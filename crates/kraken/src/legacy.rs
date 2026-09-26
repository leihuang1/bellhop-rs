use std::fs;
use std::path::{Path, PathBuf};

use crate::{Case, CaseDefinition, Diagnostic, DiagnosticReport, MAX_VECTOR_LENGTH};

const MAX_PROFILE_POINTS: usize = 100_001;

/// Load the initial Pekeris subset of a KRAKEN `.env` and FIELD `.flp` pair.
///
/// # Errors
///
/// Returns structured parse, validation, and input-file diagnostics.
#[allow(clippy::float_cmp)]
pub fn load_case(
    env_path: impl AsRef<Path>,
    flp_path: impl AsRef<Path>,
) -> Result<Case, DiagnosticReport> {
    let env_path = env_path.as_ref();
    let flp_path = flp_path.as_ref();
    let environment = parse_environment(&read_file(env_path)?, env_path)?;
    let field = parse_field(&read_file(flp_path)?, flp_path)?;

    let mut mode_sample_depths_m = environment.source_depths;
    mode_sample_depths_m.extend(environment.receiver_depths);
    mode_sample_depths_m.sort_by(f64::total_cmp);
    mode_sample_depths_m.dedup_by(|left, right| *left == *right);

    Case::from_definition(CaseDefinition {
        title: environment.title,
        frequency_hz: environment.frequency_hz,
        water_depth_m: environment.water_depth_m,
        water_sound_speed_mps: environment.water_speed,
        water_density_g_cm3: environment.water_density,
        bottom_sound_speed_mps: environment.bottom_speed,
        bottom_density_g_cm3: environment.bottom_density,
        mesh_points: environment.mesh_points,
        c_low_mps: environment.c_low,
        c_high_mps: environment.c_high,
        max_range_m: environment.max_range_m,
        mode_sample_depths_m,
        mode_limit: field.mode_limit,
        source_depths_m: field.source_depths,
        receiver_depths_m: field.receiver_depths,
        receiver_ranges_m: field.receiver_ranges_m,
        receiver_offsets_m: field.receiver_offsets_m,
    })
}

fn read_file(path: &Path) -> Result<String, DiagnosticReport> {
    fs::read_to_string(path).map_err(|error| {
        one(
            "KR0001",
            format!("unable to read input: {error}"),
            "input",
            path,
            1,
            1,
        )
    })
}

#[derive(Clone)]
struct Token {
    text: String,
    line: usize,
    column: usize,
}

struct Record {
    tokens: Vec<Token>,
    line: usize,
    slash: bool,
}

struct Reader {
    path: PathBuf,
    records: Vec<Record>,
    index: usize,
    last_line: usize,
    eof_line: usize,
}

impl Reader {
    fn new(source: &str, path: &Path) -> Result<Self, DiagnosticReport> {
        let mut records = Vec::new();
        for (index, line_text) in source.lines().enumerate() {
            let line = index + 1;
            let (tokens, slash) = tokenize(line_text, line, path)?;
            if !tokens.is_empty() || slash {
                records.push(Record {
                    tokens,
                    line,
                    slash,
                });
            }
        }
        Ok(Self {
            path: path.to_path_buf(),
            records,
            index: 0,
            last_line: 0,
            eof_line: source.lines().count() + 1,
        })
    }

    fn record(&mut self, field: &str) -> Result<Record, DiagnosticReport> {
        let Some(record) = self.records.get(self.index) else {
            return Err(one(
                "KR0101",
                "unexpected end of file",
                field,
                &self.path,
                self.eof_line,
                1,
            ));
        };
        self.index += 1;
        self.last_line = record.line;
        Ok(Record {
            tokens: record.tokens.clone(),
            line: record.line,
            slash: record.slash,
        })
    }

    fn text(&mut self, field: &str) -> Result<Token, DiagnosticReport> {
        let record = self.record(field)?;
        if let Some(token) = record.tokens.first() {
            Ok(token.clone())
        } else {
            Err(self.record_error(&record, field, "expected a value"))
        }
    }

    fn count(&mut self, field: &str) -> Result<usize, DiagnosticReport> {
        let token = self.text(field)?;
        let count = token.text.parse::<usize>().map_err(|_| {
            one(
                "KR0102",
                format!("expected an integer, got {:?}", token.text),
                field,
                &self.path,
                token.line,
                token.column,
            )
        })?;
        if count == 0 || count > MAX_VECTOR_LENGTH {
            return Err(one(
                "KR0201",
                format!("count must be in 1..={MAX_VECTOR_LENGTH}"),
                field,
                &self.path,
                token.line,
                token.column,
            ));
        }
        Ok(count)
    }

    fn scalar(&mut self, field: &str) -> Result<f64, DiagnosticReport> {
        let token = self.text(field)?;
        number(&token, &self.path, field)
    }

    fn vector(&mut self, count: usize, field: &str) -> Result<Vec<f64>, DiagnosticReport> {
        let mut result = Vec::with_capacity(count);
        while result.len() < count {
            let record = self.record(field)?;
            for token in record.tokens.iter().take(count - result.len()) {
                result.push(number(token, &self.path, field)?);
            }
            if result.len() < count && record.slash {
                return Err(self.record_error(&record, field, "slash ended an incomplete vector"));
            }
            if record.tokens.is_empty() {
                return Err(self.record_error(&record, field, "expected numeric values"));
            }
        }
        Ok(result)
    }

    fn numbers(&mut self, field: &str, count: usize) -> Result<Vec<f64>, DiagnosticReport> {
        let record = self.record(field)?;
        if record.tokens.len() != count {
            return Err(self.record_error(
                &record,
                field,
                &format!("expected exactly {count} values"),
            ));
        }
        record
            .tokens
            .iter()
            .map(|token| number(token, &self.path, field))
            .collect()
    }

    fn record_error(&self, record: &Record, field: &str, message: &str) -> DiagnosticReport {
        one(
            "KR0102",
            message,
            field,
            &self.path,
            record.line,
            record.tokens.first().map_or(1, |token| token.column),
        )
    }

    fn finish(&self) -> Result<(), DiagnosticReport> {
        if let Some(record) = self.records.get(self.index) {
            return Err(self.record_error(record, "input", "unexpected trailing input"));
        }
        Ok(())
    }
}

fn tokenize(text: &str, line: usize, path: &Path) -> Result<(Vec<Token>, bool), DiagnosticReport> {
    let chars: Vec<char> = text.chars().collect();
    let mut tokens = Vec::new();
    let mut index = 0;
    let mut slash = false;
    while index < chars.len() {
        let ch = chars[index];
        if ch.is_whitespace() || ch == ',' {
            index += 1;
        } else if ch == '!' {
            break;
        } else if ch == '/' {
            slash = true;
            break;
        } else if ch == '\'' || ch == '"' {
            let quote = ch;
            let column = index + 1;
            index += 1;
            let mut value = String::new();
            let mut closed = false;
            while index < chars.len() {
                if chars[index] == quote {
                    if chars.get(index + 1) == Some(&quote) {
                        value.push(quote);
                        index += 2;
                    } else {
                        index += 1;
                        closed = true;
                        break;
                    }
                } else {
                    value.push(chars[index]);
                    index += 1;
                }
            }
            if !closed {
                return Err(one(
                    "KR0101",
                    "unterminated quoted string",
                    "input",
                    path,
                    line,
                    column,
                ));
            }
            tokens.push(Token {
                text: value,
                line,
                column,
            });
        } else {
            let start = index;
            while index < chars.len()
                && !chars[index].is_whitespace()
                && !matches!(chars[index], ',' | '/' | '!')
            {
                index += 1;
            }
            tokens.push(Token {
                text: chars[start..index].iter().collect(),
                line,
                column: start + 1,
            });
        }
    }
    Ok((tokens, slash))
}

fn number(token: &Token, path: &Path, field: &str) -> Result<f64, DiagnosticReport> {
    token
        .text
        .replace('D', "E")
        .replace('d', "e")
        .parse::<f64>()
        .map_err(|_| {
            one(
                "KR0102",
                format!("expected a number, got {:?}", token.text),
                field,
                path,
                token.line,
                token.column,
            )
        })
}

struct Environment {
    title: String,
    frequency_hz: f64,
    water_depth_m: f64,
    water_speed: f64,
    water_density: f64,
    bottom_speed: f64,
    bottom_density: f64,
    mesh_points: usize,
    c_low: f64,
    c_high: f64,
    max_range_m: f64,
    source_depths: Vec<f64>,
    receiver_depths: Vec<f64>,
}

#[allow(clippy::float_cmp, clippy::too_many_lines)]
fn parse_environment(source: &str, path: &Path) -> Result<Environment, DiagnosticReport> {
    let mut reader = Reader::new(source, path)?;
    let title = reader.text("title")?.text;
    let frequency_hz = reader.scalar("frequency_hz")?;
    if reader.count("medium_count")? != 1 {
        return Err(reader_error(
            &reader,
            "KR0202",
            "only one water medium is supported in this slice",
            "medium_count",
        ));
    }

    let options = reader.text("top_options")?;
    let option = |index| options.text.as_bytes().get(index).copied().unwrap_or(b' ');
    if [
        option(0),
        option(1),
        option(2),
        option(3),
        option(4),
        option(5),
    ] != [b'N', b'V', b'N', b' ', b' ', b' ']
        || options
            .text
            .as_bytes()
            .iter()
            .skip(6)
            .any(|byte| !byte.is_ascii_whitespace())
    {
        return Err(one(
            "KR0202",
            "requires N interpolation, a vacuum surface, no attenuation, and one frequency",
            "top_options",
            path,
            options.line,
            options.column,
        ));
    }

    let header = reader.record("water_header")?;
    if header.tokens.len() != 3 {
        return Err(reader.record_error(
            &header,
            "water_header",
            "expected mesh count, roughness, and depth",
        ));
    }
    let mesh_points = header.tokens[0].text.parse::<usize>().map_err(|_| {
        reader.record_error(&header, "mesh_points", "mesh count must be an integer")
    })?;
    let surface_roughness = number(&header.tokens[1], path, "surface_roughness")?;
    let water_depth_m = number(&header.tokens[2], path, "water_depth_m")?;

    let mut points: Vec<[f64; 6]> = Vec::new();
    loop {
        let record = reader.record("sound_speed_profile")?;
        if !(2..=6).contains(&record.tokens.len())
            || (points.is_empty() && record.tokens.len() != 6)
        {
            return Err(reader.record_error(
                &record,
                "sound_speed_profile",
                "expected 2..=6 values; the first point must specify all 6",
            ));
        }
        let mut point = points
            .last()
            .copied()
            .unwrap_or([0.0, 0.0, 0.0, 1.0, 0.0, 0.0]);
        for (index, token) in record.tokens.iter().enumerate() {
            point[index] = number(token, path, "sound_speed_profile")?;
        }
        points.push(point);
        if points.len() > MAX_PROFILE_POINTS {
            return Err(reader.record_error(
                &record,
                "sound_speed_profile",
                "too many profile points",
            ));
        }
        if record.slash {
            break;
        }
    }
    if points.len() < 2 {
        return Err(reader_error(
            &reader,
            "KR0202",
            "a Pekeris profile needs top and interface points",
            "sound_speed_profile",
        ));
    }
    let water_speed = points[0][1];
    let water_density = points[0][3];
    if points.iter().any(|point| {
        point[1] != water_speed
            || point[2] != 0.0
            || point[3] != water_density
            || point[4] != 0.0
            || point[5] != 0.0
    }) || points[0][0] != 0.0
        || points.last().is_none_or(|point| point[0] != water_depth_m)
        || points.windows(2).any(|pair| pair[1][0] <= pair[0][0])
        || surface_roughness != 0.0
    {
        return Err(reader_error(
            &reader,
            "KR0202",
            "requires a smooth, constant, lossless fluid water column",
            "sound_speed_profile",
        ));
    }

    let bottom_option = reader.record("bottom_options")?;
    if bottom_option.tokens.len() != 2
        || bottom_option.tokens[0].text != "A"
        || number(&bottom_option.tokens[1], path, "bottom_roughness")? != 0.0
    {
        return Err(reader.record_error(
            &bottom_option,
            "bottom_options",
            "requires a smooth acoustic fluid half-space",
        ));
    }
    let bottom = reader.numbers("bottom_half_space", 6)?;
    if bottom[0] != water_depth_m || bottom[2] != 0.0 || bottom[4] != 0.0 || bottom[5] != 0.0 {
        return Err(reader_error(
            &reader,
            "KR0202",
            "bottom half-space must be fluid, lossless, and begin at the interface",
            "bottom_half_space",
        ));
    }
    let limits = reader.numbers("phase_speed_limits", 2)?;
    let max_range_m = reader.scalar("max_range_km")? * 1000.0;
    let source_depths = read_vector(&mut reader, "mode_source_depths_m")?;
    let receiver_depths = read_vector(&mut reader, "mode_receiver_depths_m")?;
    reader.finish()?;

    Ok(Environment {
        title,
        frequency_hz,
        water_depth_m,
        water_speed,
        water_density,
        bottom_speed: bottom[1],
        bottom_density: bottom[3],
        mesh_points,
        c_low: limits[0],
        c_high: limits[1],
        max_range_m,
        source_depths,
        receiver_depths,
    })
}

fn read_vector(reader: &mut Reader, field: &str) -> Result<Vec<f64>, DiagnosticReport> {
    let count = reader.count(&format!("{field}.count"))?;
    reader.vector(count, field)
}

struct Field {
    mode_limit: usize,
    source_depths: Vec<f64>,
    receiver_depths: Vec<f64>,
    receiver_ranges_m: Vec<f64>,
    receiver_offsets_m: Vec<f64>,
}

fn parse_field(source: &str, path: &Path) -> Result<Field, DiagnosticReport> {
    let mut reader = Reader::new(source, path)?;
    reader.text("field_title")?;
    let options = reader.text("field_options")?;
    let chars: Vec<char> = options.text.chars().collect();
    let option = |index| chars.get(index).copied().unwrap_or(' ');
    if option(0) != 'X'
        || option(1) != ' '
        || !matches!(option(2), ' ' | 'O')
        || !matches!(option(3), ' ' | 'C')
        || chars.iter().skip(4).any(|ch| !ch.is_whitespace())
    {
        return Err(one(
            "KR0202",
            "requires a coherent, omnidirectional line source",
            "field_options",
            path,
            options.line,
            options.column,
        ));
    }
    let mode_limit = reader.count("mode_limit")?;
    let profiles = reader.count("profile_count")?;
    let profile_ranges = reader.vector(profiles, "profile_ranges_km")?;
    if profiles != 1 || profile_ranges[0] != 0.0 {
        return Err(reader_error(
            &reader,
            "KR0202",
            "requires one range-independent profile at 0 km",
            "profile_ranges_km",
        ));
    }
    let ranges_km = read_vector(&mut reader, "receiver_ranges_km")?;
    let source_depths = read_vector(&mut reader, "field_source_depths_m")?;
    let receiver_depths = read_vector(&mut reader, "field_receiver_depths_m")?;
    let receiver_offsets = read_vector(&mut reader, "receiver_offsets_m")?;
    reader.finish()?;

    Ok(Field {
        mode_limit,
        source_depths,
        receiver_depths,
        receiver_ranges_m: ranges_km.into_iter().map(|range| range * 1000.0).collect(),
        receiver_offsets_m: receiver_offsets,
    })
}

fn reader_error(
    reader: &Reader,
    code: &'static str,
    message: &str,
    field: &str,
) -> DiagnosticReport {
    one(code, message, field, &reader.path, reader.last_line, 1)
}

fn one(
    code: &'static str,
    message: impl Into<String>,
    field: impl Into<String>,
    path: impl Into<PathBuf>,
    line: usize,
    column: usize,
) -> DiagnosticReport {
    DiagnosticReport::one(Diagnostic::new(code, message, field, path, line, column))
}

#[cfg(test)]
mod tests {
    use super::{parse_environment, parse_field};
    use std::path::Path;

    #[test]
    fn unsupported_solver_and_field_options_are_rejected() {
        let env = include_str!("../tests/fixtures/Pekeris.env").replace("'NVN'", "'NVW'");
        let error = parse_environment(&env, Path::new("Pekeris.env"))
            .err()
            .expect("attenuation option should be rejected");
        assert_eq!(error.diagnostics()[0].field, "top_options");
        assert_eq!(error.diagnostics()[0].line, 4);

        let flp = include_str!("../tests/fixtures/Pekeris.flp").replace("'X OC'", "'R OC'");
        let error = parse_field(&flp, Path::new("Pekeris.flp"))
            .err()
            .expect("point source should be rejected in this slice");
        assert_eq!(error.diagnostics()[0].field, "field_options");
        assert_eq!(error.diagnostics()[0].line, 2);
    }
}
