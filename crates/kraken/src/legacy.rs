use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};

use crate::{Case, CaseDefinition, Diagnostic, DiagnosticReport, MAX_VECTOR_LENGTH};

const MAX_PROFILE_POINTS: usize = 100_001;

/// Load the initial Pekeris subset of a KRAKEN `.env` and FIELD `.flp` pair.
///
/// # Errors
///
/// Returns structured parse, validation, and input-file diagnostics.
pub fn load_case(
    env_path: impl AsRef<Path>,
    flp_path: impl AsRef<Path>,
) -> Result<Case, DiagnosticReport> {
    let env_path = env_path.as_ref();
    let flp_path = flp_path.as_ref();
    parse_case(
        &read_file(env_path)?,
        &read_file(flp_path)?,
        env_path,
        flp_path,
    )
}

#[allow(clippy::float_cmp)]
fn parse_case(
    env_source: &str,
    flp_source: &str,
    env_path: &Path,
    flp_path: &Path,
) -> Result<Case, DiagnosticReport> {
    let environment = parse_environment(env_source, env_path)?;
    let field = parse_field(flp_source, flp_path)?;

    let mut mode_sample_depths_m = environment.source_depths;
    mode_sample_depths_m.extend(environment.receiver_depths);
    mode_sample_depths_m.sort_by(f64::total_cmp);
    mode_sample_depths_m.dedup_by(|left, right| *left == *right);

    let env_locations = environment.locations;
    let field_locations = field.locations;
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
    .map_err(|mut report| {
        for diagnostic in &mut report.diagnostics {
            let (locations, path, record) = match diagnostic.field.as_str() {
                "water_depth_m" | "mesh_points" => (&env_locations, env_path, "water_header"),
                "water_sound_speed_mps" | "water_density_g_cm3" => {
                    (&env_locations, env_path, "sound_speed_profile")
                }
                "bottom_sound_speed_mps" | "bottom_density_g_cm3" => {
                    (&env_locations, env_path, "bottom_half_space")
                }
                "max_range_m" => (&env_locations, env_path, "max_range_km"),
                "mode_sample_depths_m" => (&env_locations, env_path, "mode_receiver_depths_m"),
                "source_depths_m" => (&field_locations, flp_path, "field_source_depths_m"),
                "receiver_depths_m" => (&field_locations, flp_path, "field_receiver_depths_m"),
                "receiver_ranges_m" | "field_grid" => {
                    (&field_locations, flp_path, "receiver_ranges_km")
                }
                "receiver_offsets_m" => (&field_locations, flp_path, "receiver_offsets_m"),
                "mode_limit" => (&field_locations, flp_path, "mode_limit"),
                other => (&env_locations, env_path, other),
            };
            if let Some(&(line, column)) = locations.get(record) {
                diagnostic.path = path.to_path_buf();
                diagnostic.line = line;
                diagnostic.column = column;
            }
        }
        report
    })
}

fn read_file(path: &Path) -> Result<String, DiagnosticReport> {
    const MAX_INPUT_BYTES: u64 = 1_048_576;
    let mut source = String::new();
    File::open(path)
        .and_then(|file| file.take(MAX_INPUT_BYTES + 1).read_to_string(&mut source))
        .map_err(|error| {
            one(
                "KR0001",
                format!("unable to read input: {error}"),
                "input",
                path,
                1,
                1,
            )
        })?;
    if source.len() as u64 > MAX_INPUT_BYTES {
        return Err(one(
            "KR0201",
            "input file exceeds 1 MiB",
            "input",
            path,
            1,
            1,
        ));
    }
    Ok(source)
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
    locations: HashMap<String, (usize, usize)>,
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
            locations: HashMap::new(),
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
        self.locations.entry(field.to_owned()).or_insert((
            record.line,
            record.tokens.first().map_or(1, |token| token.column),
        ));
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
    let mut comma_pending = false;
    while index < chars.len() {
        let ch = chars[index];
        if ch.is_whitespace() {
            index += 1;
        } else if ch == ',' {
            if comma_pending || tokens.is_empty() {
                return Err(one(
                    "KR0202",
                    "Fortran null slots are not supported",
                    "input",
                    path,
                    line,
                    index + 1,
                ));
            }
            comma_pending = true;
            index += 1;
        } else if ch == '!' {
            break;
        } else if ch == '/' {
            slash = true;
            break;
        } else if ch == '\'' || ch == '"' {
            comma_pending = false;
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
            comma_pending = false;
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
    let value = token
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
        })?;
    if !value.is_finite() {
        return Err(one(
            "KR0201",
            "value must be finite",
            field,
            path,
            token.line,
            token.column,
        ));
    }
    Ok(value)
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
    locations: HashMap<String, (usize, usize)>,
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
    ] != *b"NVN   "
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
    if surface_roughness != 0.0 {
        return Err(reader.record_error(&header, "surface_roughness", "requires a smooth surface"));
    }

    let mut points: Vec<[f64; 6]> = Vec::new();
    loop {
        let record = reader.record("sound_speed_profile")?;
        if !(2..=6).contains(&record.tokens.len())
            || (points.is_empty() && record.tokens.len() != 6)
            || (record.tokens.len() < 6 && !record.slash)
        {
            return Err(reader.record_error(
                &record,
                "sound_speed_profile",
                "expected 6 values, or 2..=5 followed by / to inherit trailing values; the first point must specify all 6",
            ));
        }
        let mut point = points
            .last()
            .copied()
            .unwrap_or([0.0, 0.0, 0.0, 1.0, 0.0, 0.0]);
        for (index, token) in record.tokens.iter().enumerate() {
            point[index] = number(token, path, "sound_speed_profile")?;
        }
        if point[0] < 0.0
            || point[0] > water_depth_m
            || points
                .last()
                .is_some_and(|previous| point[0] <= previous[0])
        {
            return Err(reader.record_error(
                &record,
                "sound_speed_profile",
                "depths must increase within the water column",
            ));
        }
        if (points.is_empty() && point[0] != 0.0)
            || point[2] != 0.0
            || point[4] != 0.0
            || point[5] != 0.0
            || points
                .first()
                .is_some_and(|first| point[1] != first[1] || point[3] != first[3])
        {
            return Err(reader.record_error(
                &record,
                "sound_speed_profile",
                "requires a constant, lossless fluid water column starting at 0 m",
            ));
        }
        points.push(point);
        if points.len() > MAX_PROFILE_POINTS {
            return Err(reader.record_error(
                &record,
                "sound_speed_profile",
                "too many profile points",
            ));
        }
        if point[0] == water_depth_m {
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
        locations: reader.locations,
    })
}

fn read_vector(reader: &mut Reader, field: &str) -> Result<Vec<f64>, DiagnosticReport> {
    let count = reader.count(&format!("{field}.count"))?;
    let mut values = reader.vector(count, field)?;
    // The reference ReadVector sorts each vector, including receiver offsets.
    values.sort_by(f64::total_cmp);
    Ok(values)
}

struct Field {
    mode_limit: usize,
    source_depths: Vec<f64>,
    receiver_depths: Vec<f64>,
    receiver_ranges_m: Vec<f64>,
    receiver_offsets_m: Vec<f64>,
    locations: HashMap<String, (usize, usize)>,
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
        locations: reader.locations,
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
    use super::{parse_case, parse_environment, parse_field, read_file};
    use std::path::Path;

    #[test]
    fn oversized_files_are_rejected_before_parsing() {
        let stamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path =
            std::env::temp_dir().join(format!("kraken-input-{}-{stamp}.env", std::process::id()));
        let file = std::fs::File::create_new(&path).unwrap();
        file.set_len(1_048_577).unwrap();
        drop(file);
        let result = read_file(&path);
        std::fs::remove_file(&path).unwrap();
        assert_eq!(result.unwrap_err().diagnostics()[0].code, "KR0201");
    }

    #[test]
    fn semantic_errors_retain_file_locations() {
        let env = include_str!("../tests/fixtures/Pekeris.env").replace("50.0\n", "-50.0\n");
        let flp = include_str!("../tests/fixtures/Pekeris.flp").replace("75.0 /", "175.0 /");
        let report =
            parse_case(&env, &flp, Path::new("case.env"), Path::new("case.flp")).unwrap_err();
        for (field, path, line) in [
            ("frequency_hz", "case.env", 2),
            ("source_depths_m", "case.flp", 9),
        ] {
            let diagnostic = report
                .diagnostics()
                .iter()
                .find(|d| d.field == field)
                .unwrap();
            assert_eq!(diagnostic.path, Path::new(path));
            assert_eq!(diagnostic.line, line);
        }
    }

    #[test]
    fn legacy_vectors_follow_reference_sorting_and_reject_unsupported_syntax() {
        let flp = include_str!("../tests/fixtures/Pekeris.flp");
        let source = flp
            .replace("0.5 1.0 2.0 /", "2.0D0, 0.5d0, 1.0 /")
            .replace("0.0 0.0 0.0 /", "10.0 -10.0 0.0 /");
        let field = parse_field(&source, Path::new("case.flp")).unwrap();
        assert_eq!(field.receiver_ranges_m, [500.0, 1000.0, 2000.0]);
        assert_eq!(field.receiver_offsets_m, [-10.0, 0.0, 10.0]);
        for vector in [
            "0.5 2.0 /",
            "3*0.5 /",
            ",0.5,1.0 /",
            "0.5,,1.0 /",
            "0.5 NaN 2.0 /",
        ] {
            let source = flp.replace("0.5 1.0 2.0 /", vector);
            assert!(
                parse_field(&source, Path::new("case.flp")).is_err(),
                "accepted {vector}"
            );
        }
    }

    #[test]
    fn ssp_ends_at_interface_not_slash() {
        let env = include_str!("../tests/fixtures/Pekeris.env");
        for source in [
            env.replace(
                "0.0 1500.0 0.0 1.0 0.0 0.0\n",
                "0.0 1500.0 0.0 1.0 0.0 0.0 /\n",
            ),
            env.replace("100.0 1500.0 /", "100.0 1500.0 0.0 1.0 0.0 0.0"),
        ] {
            assert!(parse_environment(&source, Path::new("Pekeris.env")).is_ok());
        }
    }

    #[test]
    fn invalid_ssp_numbers_and_ambiguous_records_are_rejected() {
        let env = include_str!("../tests/fixtures/Pekeris.env");
        for point in [
            "NaN 1500.0 0.0 1.0 0.0 0.0",
            "inf 1500.0 0.0 1.0 0.0 0.0",
            "50.0,,1500.0,0.0,1.0,0.0,0.0",
            "50.0 1500.0", // Missing values need an explicit slash, not implicit inheritance.
        ] {
            let source = env.replace("100.0 1500.0 /", &format!("{point}\n100.0 1500.0 /"));
            assert!(
                parse_environment(&source, Path::new("Pekeris.env")).is_err(),
                "accepted {point}"
            );
        }
    }

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
