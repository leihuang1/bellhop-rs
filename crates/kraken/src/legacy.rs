use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};

use crate::{
    BottomBoundary, Case, CaseDefinition, Diagnostic, DiagnosticReport, Interpolation,
    MAX_VECTOR_LENGTH, ModeSolver, SoundSpeedPoint, SourceGeometry, SurfaceBoundary,
};

/// Maximum byte length of each legacy input, including source-based loading.
pub const MAX_INPUT_BYTES: u64 = 1_048_576;

const MAX_PROFILE_POINTS: usize = MAX_VECTOR_LENGTH;
const MAX_FREQUENCIES: usize = 1000;
const MAX_FREQUENCY_INPUT_VALUES: usize = 5_000_000;

/// Load the supported single-fluid subset of a KRAKEN `.env` and FIELD `.flp` pair.
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

/// Load a KRAKENC environment and coherent, range-independent FIELD geometry.
///
/// # Errors
///
/// Returns structured input diagnostics for unsupported configurations.
pub fn load_complex_case(
    env_path: impl AsRef<Path>,
    flp_path: impl AsRef<Path>,
) -> Result<Case, DiagnosticReport> {
    let env_path = env_path.as_ref();
    let flp_path = flp_path.as_ref();
    parse_case_with_solver(
        &read_file(env_path)?,
        &read_file(flp_path)?,
        env_path,
        flp_path,
        ModeSolver::Krakenc,
    )
}

/// Load one validated case per frequency, in legacy input order.
///
/// Each case can be passed to `solve`; results remain frequency-domain products.
/// Single-frequency inputs also work. Geometry/profile copies are bounded.
///
/// # Errors
///
/// Returns structured parse, mesh-scaling, and case validation diagnostics.
pub fn load_frequency_cases(
    env_path: impl AsRef<Path>,
    flp_path: impl AsRef<Path>,
    mode_solver: ModeSolver,
) -> Result<Vec<Case>, DiagnosticReport> {
    let env_path = env_path.as_ref();
    let flp_path = flp_path.as_ref();
    load_frequency_cases_from_sources(
        &read_file(env_path)?,
        &read_file(flp_path)?,
        env_path,
        flp_path,
        mode_solver,
    )
}

/// Load frequency cases from exact UTF-8 input snapshots, preserving source locations.
///
/// This lets adapters record hashes of the bytes actually parsed, without rereading files.
/// The same input-size and validation limits as the file loader apply.
///
/// # Errors
///
/// Returns structured input-size, parse, mesh-scaling, and case diagnostics.
pub fn load_frequency_cases_from_sources(
    env_source: &str,
    flp_source: &str,
    env_path: &Path,
    flp_path: &Path,
    mode_solver: ModeSolver,
) -> Result<Vec<Case>, DiagnosticReport> {
    for (source, path) in [(env_source, env_path), (flp_source, flp_path)] {
        check_input_size(source, path)?;
    }
    parse_frequency_cases(
        env_source,
        flp_source,
        env_path,
        flp_path,
        mode_solver,
        false,
    )
}

fn parse_case(
    env_source: &str,
    flp_source: &str,
    env_path: &Path,
    flp_path: &Path,
) -> Result<Case, DiagnosticReport> {
    parse_case_with_solver(
        env_source,
        flp_source,
        env_path,
        flp_path,
        ModeSolver::Kraken,
    )
}

fn parse_case_with_solver(
    env_source: &str,
    flp_source: &str,
    env_path: &Path,
    flp_path: &Path,
    mode_solver: ModeSolver,
) -> Result<Case, DiagnosticReport> {
    parse_frequency_cases(
        env_source,
        flp_source,
        env_path,
        flp_path,
        mode_solver,
        true,
    )
    .map(|mut cases| cases.pop().unwrap())
}

#[allow(clippy::float_cmp, clippy::too_many_lines)]
fn parse_frequency_cases(
    env_source: &str,
    flp_source: &str,
    env_path: &Path,
    flp_path: &Path,
    mode_solver: ModeSolver,
    single_frequency: bool,
) -> Result<Vec<Case>, DiagnosticReport> {
    let mut environment = parse_environment_with_solver(env_source, env_path, mode_solver)?;
    if single_frequency && environment.frequencies_hz.len() != 1 {
        let (line, column) = environment.locations["frequencies_hz"];
        return Err(one(
            "KR0202",
            "multiple frequencies require load_frequency_cases",
            "frequencies_hz",
            env_path,
            line,
            column,
        ));
    }
    let mut field = parse_field(flp_source, flp_path)?;
    // ReadSzRz stores depths in single precision. Keep an interface sample on
    // the exact validated f64 boundary even when the f32 spelling rounds upward.
    #[allow(clippy::cast_possible_truncation)]
    let boundary = f64::from(environment.water_depth_m as f32);
    for depths in [
        &mut environment.source_depths,
        &mut environment.receiver_depths,
        &mut field.source_depths,
        &mut field.receiver_depths,
    ] {
        for depth in depths {
            if *depth == boundary {
                *depth = environment.water_depth_m;
            }
        }
    }

    let mut mode_sample_depths_m = environment.source_depths;
    mode_sample_depths_m.extend(environment.receiver_depths);
    mode_sample_depths_m.sort_by(f64::total_cmp);
    mode_sample_depths_m.dedup_by(|left, right| *left == *right);

    let bottom_location = if environment.bottom_boundary == BottomBoundary::Rigid {
        "bottom_options"
    } else {
        "bottom_half_space"
    };
    let env_locations = environment.locations;
    let field_locations = field.locations;
    let definition = CaseDefinition {
        title: environment.title,
        mode_solver,
        frequency_hz: environment.frequency_hz,
        mesh_reference_frequency_hz: environment.broadband.then_some(environment.frequency_hz),
        water_depth_m: environment.water_depth_m,
        interpolation: environment.interpolation,
        surface_boundary: environment.surface_boundary,
        sound_speed_profile: environment.profile,
        water_density_g_cm3: environment.water_density,
        bottom_boundary: environment.bottom_boundary,
        bottom_sound_speed_mps: environment.bottom_speed,
        bottom_density_g_cm3: environment.bottom_density,
        bottom_attenuation_db_per_wavelength: environment.bottom_attenuation,
        source_geometry: field.source_geometry,
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
    };
    let values = definition.sound_speed_profile.len()
        + definition.mode_sample_depths_m.len()
        + definition.source_depths_m.len()
        + definition.receiver_depths_m.len()
        + definition.receiver_ranges_m.len()
        + definition.receiver_offsets_m.len();
    if values
        .checked_mul(environment.frequencies_hz.len())
        .is_none_or(|n| n > MAX_FREQUENCY_INPUT_VALUES)
    {
        let &(line, column) = env_locations
            .get("frequencies_hz")
            .unwrap_or(&env_locations["frequency_hz"]);
        return Err(one(
            "KR0201",
            "frequency cases exceed the input storage limit",
            "frequencies_hz",
            env_path,
            line,
            column,
        ));
    }
    environment
        .frequencies_hz
        .into_iter()
        .map(|frequency_hz| {
            let mut input = definition.clone();
            input.frequency_hz = frequency_hz;
            Case::from_definition(input)
                .and_then(|case| {
                    if case.mesh_reference_frequency_hz.is_some() {
                        case.mesh_points_at(1)?;
                    }
                    Ok(case)
                })
                .map_err(|mut report| {
                    for diagnostic in &mut report.diagnostics {
                        let (locations, path, record) = match diagnostic.field.as_str() {
                            "water_depth_m" | "mesh_points" => {
                                (&env_locations, env_path, "water_header")
                            }
                            "sound_speed_profile" | "water_density_g_cm3" => {
                                (&env_locations, env_path, "sound_speed_profile")
                            }
                            "bottom_sound_speed_mps"
                            | "bottom_density_g_cm3"
                            | "bottom_attenuation_db_per_wavelength" => {
                                (&env_locations, env_path, bottom_location)
                            }
                            "max_range_m" => (&env_locations, env_path, "max_range_km"),
                            "mesh_reference_frequency_hz" => {
                                (&env_locations, env_path, "frequency_hz")
                            }
                            "mode_sample_depths_m" => {
                                (&env_locations, env_path, "mode_receiver_depths_m")
                            }
                            "source_depths_m" => {
                                (&field_locations, flp_path, "field_source_depths_m")
                            }
                            "receiver_depths_m" => {
                                (&field_locations, flp_path, "field_receiver_depths_m")
                            }
                            "receiver_ranges_m" | "field_grid" => {
                                (&field_locations, flp_path, "receiver_ranges_km")
                            }
                            "receiver_offsets_m" => {
                                (&field_locations, flp_path, "receiver_offsets_m")
                            }
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
        })
        .collect()
}

fn read_file(path: &Path) -> Result<String, DiagnosticReport> {
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
    check_input_size(&source, path)?;
    Ok(source)
}

fn check_input_size(source: &str, path: &Path) -> Result<(), DiagnosticReport> {
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
    Ok(())
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
            if record.slash {
                break;
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
    frequencies_hz: Vec<f64>,
    broadband: bool,
    water_depth_m: f64,
    interpolation: Interpolation,
    surface_boundary: SurfaceBoundary,
    profile: Vec<SoundSpeedPoint>,
    water_density: f64,
    bottom_boundary: BottomBoundary,
    bottom_speed: f64,
    bottom_density: f64,
    bottom_attenuation: f64,
    mesh_points: usize,
    c_low: f64,
    c_high: f64,
    max_range_m: f64,
    source_depths: Vec<f64>,
    receiver_depths: Vec<f64>,
    locations: HashMap<String, (usize, usize)>,
}

#[cfg(test)]
fn parse_environment(source: &str, path: &Path) -> Result<Environment, DiagnosticReport> {
    parse_environment_with_solver(source, path, ModeSolver::Kraken)
}

#[allow(clippy::float_cmp, clippy::too_many_lines)]
fn parse_environment_with_solver(
    source: &str,
    path: &Path,
    mode_solver: ModeSolver,
) -> Result<Environment, DiagnosticReport> {
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
    if !matches!(option(0), b'N' | b'C' | b'P' | b'S' | b'A')
        || !matches!(option(1), b'V' | b'R')
        || !matches!(option(2), b'N' | b'W')
        || option(3) != b' '
        || (mode_solver == ModeSolver::Kraken && option(4) != b' ')
        || (mode_solver == ModeSolver::Krakenc && !matches!(option(4), b' ' | b'.'))
        || !matches!(option(5), b' ' | b'B')
        || options
            .text
            .as_bytes()
            .iter()
            .skip(6)
            .any(|byte| !byte.is_ascii_whitespace())
    {
        return Err(one(
            "KR0202",
            "requires N/C/P/S or fixed analytic A interpolation, vacuum or rigid surface, N/W attenuation without water loss, and optional B frequencies",
            "top_options",
            path,
            options.line,
            options.column,
        ));
    }

    let interpolation = match option(0) {
        b'N' => Interpolation::N2Linear,
        b'C' => Interpolation::CLinear,
        b'P' => Interpolation::Pchip,
        b'S' => Interpolation::Spline,
        _ => Interpolation::AnalyticMunk,
    };
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
    if interpolation != Interpolation::AnalyticMunk {
        loop {
            let record = reader.record("sound_speed_profile")?;
            if !(2..=6).contains(&record.tokens.len()) || (record.tokens.len() < 6 && !record.slash)
            {
                return Err(reader.record_error(
                    &record,
                    "sound_speed_profile",
                    "expected 6 values, or 2..=5 followed by / to inherit trailing values",
                ));
            }
            let mut point = points
                .last()
                .copied()
                .unwrap_or([0.0, 1500.0, 0.0, 1.0, 0.0, 0.0]);
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
                || points.first().is_some_and(|first| point[3] != first[3])
            {
                return Err(reader.record_error(
                    &record,
                    "sound_speed_profile",
                    "requires a lossless, constant-density fluid water column starting at 0 m",
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
                "a fluid profile needs top and interface points",
                "sound_speed_profile",
            ));
        }
    }
    let water_density = points.first().map_or(1.0, |point| point[3]);

    let bottom_option = reader.record("bottom_options")?;
    if bottom_option.tokens.len() != 2
        || !matches!(bottom_option.tokens[0].text.as_str(), "A" | "R")
        || number(&bottom_option.tokens[1], path, "bottom_roughness")? != 0.0
        || (bottom_option.tokens[0].text == "R" && option(2) != b'N')
    {
        return Err(reader.record_error(
            &bottom_option,
            "bottom_options",
            "requires a smooth acoustic fluid half-space or rigid bottom without loss units",
        ));
    }
    let bottom_boundary = if bottom_option.tokens[0].text == "R" {
        BottomBoundary::Rigid
    } else {
        BottomBoundary::FluidHalfSpace
    };
    let mut bottom = [0.0; 6];
    if bottom_boundary == BottomBoundary::FluidHalfSpace {
        let bottom_record = reader.record("bottom_half_space")?;
        if !(1..=6).contains(&bottom_record.tokens.len())
            || (bottom_record.tokens.len() < 6 && !bottom_record.slash)
        {
            return Err(reader.record_error(
                &bottom_record,
                "bottom_half_space",
                "expected 6 values, or trailing defaults terminated by /",
            ));
        }
        if interpolation == Interpolation::AnalyticMunk && bottom_record.tokens.len() < 4 {
            return Err(reader.record_error(
                &bottom_record,
                "bottom_half_space",
                "analytic profile requires explicit bottom sound speed and density",
            ));
        }
        bottom = points
            .last()
            .copied()
            .unwrap_or([water_depth_m, 1500.0, 0.0, 1.0, 0.0, 0.0]);
        for (index, token) in bottom_record.tokens.iter().enumerate() {
            bottom[index] = number(token, path, "bottom_half_space")?;
        }
        if bottom[0] != water_depth_m
            || bottom[2] != 0.0
            || bottom[5] != 0.0
            || (option(2) != b'W' && bottom[4] != 0.0)
        {
            return Err(reader_error(
                &reader,
                "KR0202",
                "bottom half-space must be fluid, start at the interface, and use supported loss units",
                "bottom_half_space",
            ));
        }
    }
    let limits = reader.numbers("phase_speed_limits", 2)?;
    let max_range_m = reader.scalar("max_range_km")? * 1000.0;
    let source_depths = read_vector(&mut reader, "mode_source_depths_m", true)?;
    let receiver_depths = read_vector(&mut reader, "mode_receiver_depths_m", true)?;
    let broadband = option(5) == b'B';
    let frequencies_hz = if broadband {
        let count = reader.count("frequencies_hz.count")?;
        if count > MAX_FREQUENCIES {
            return Err(reader_error(
                &reader,
                "KR0201",
                "at most 1000 frequencies are supported",
                "frequencies_hz.count",
            ));
        }
        // Unlike ReadVector, ReadfreqVec/SubTab preserves frequency order.
        let frequencies = read_vector_values(&mut reader, "frequencies_hz", count, false)?;
        if frequencies.iter().any(|&f| f <= 0.0) {
            return Err(reader_error(
                &reader,
                "KR0201",
                "frequencies must be positive",
                "frequencies_hz",
            ));
        }
        frequencies
    } else {
        vec![frequency_hz]
    };
    reader.finish()?;

    Ok(Environment {
        title,
        frequency_hz,
        frequencies_hz,
        broadband,
        water_depth_m,
        interpolation,
        surface_boundary: if option(1) == b'V' {
            SurfaceBoundary::Vacuum
        } else {
            SurfaceBoundary::Rigid
        },
        profile: points
            .into_iter()
            .map(|point| SoundSpeedPoint {
                depth_m: point[0],
                sound_speed_mps: point[1],
            })
            .collect(),
        water_density,
        bottom_boundary,
        bottom_speed: bottom[1],
        bottom_density: bottom[3],
        bottom_attenuation: bottom[4],
        mesh_points,
        c_low: limits[0],
        c_high: limits[1],
        max_range_m,
        source_depths,
        receiver_depths,
        locations: reader.locations,
    })
}

#[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
fn read_vector(
    reader: &mut Reader,
    field: &str,
    single_precision: bool,
) -> Result<Vec<f64>, DiagnosticReport> {
    let count = reader.count(&format!("{field}.count"))?;
    let mut values = read_vector_values(reader, field, count, single_precision)?;
    // The reference ReadVector sorts each vector, including receiver offsets.
    values.sort_by(f64::total_cmp);
    Ok(values)
}

#[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
fn read_vector_values(
    reader: &mut Reader,
    field: &str,
    count: usize,
    single_precision: bool,
) -> Result<Vec<f64>, DiagnosticReport> {
    let mut values = reader.vector(count, field)?;
    if single_precision {
        for value in &mut values {
            *value = f64::from(*value as f32);
        }
    }
    if values.len() < count {
        if count < 3 || !(1..=2).contains(&values.len()) {
            return Err(reader_error(
                reader,
                "KR0102",
                "incomplete vector: supply all values or 1..=2 endpoints for count >= 3",
                field,
            ));
        }
        let start = values[0];
        let end = *values.last().unwrap();
        values = if single_precision {
            let step = (end as f32 - start as f32) / (count - 1) as f32;
            (0..count)
                .map(|i| f64::from(start as f32 + i as f32 * step))
                .collect()
        } else {
            let step = (end - start) / (count - 1) as f64;
            (0..count).map(|i| start + i as f64 * step).collect()
        };
    }
    if values.iter().any(|value| !value.is_finite()) {
        return Err(reader_error(
            reader,
            "KR0201",
            "vector exceeds the numeric range",
            field,
        ));
    }
    Ok(values)
}

struct Field {
    mode_limit: usize,
    source_geometry: SourceGeometry,
    source_depths: Vec<f64>,
    receiver_depths: Vec<f64>,
    receiver_ranges_m: Vec<f64>,
    receiver_offsets_m: Vec<f64>,
    locations: HashMap<String, (usize, usize)>,
}

fn parse_field(source: &str, path: &Path) -> Result<Field, DiagnosticReport> {
    let mut reader = Reader::new(source, path)?;
    // A slash-only title means use the environment title in FIELD.
    reader.record("field_title")?;
    let options = reader.text("field_options")?;
    let chars: Vec<char> = options.text.chars().collect();
    let option = |index| chars.get(index).copied().unwrap_or(' ');
    if !matches!(option(0), 'X' | 'R')
        || !matches!(option(1), ' ' | 'A')
        || !matches!(option(2), ' ' | 'O')
        || !matches!(option(3), ' ' | 'C')
        || chars.iter().skip(4).any(|ch| !ch.is_whitespace())
    {
        return Err(one(
            "KR0202",
            "requires a coherent, omnidirectional line or point source",
            "field_options",
            path,
            options.line,
            options.column,
        ));
    }
    let mode_limit = reader.count("mode_limit")?;
    let profiles = reader.count("profile_count")?;
    let profile_ranges = reader.vector(profiles, "profile_ranges_km")?;
    if profiles != 1 || profile_ranges.len() != 1 || profile_ranges[0] != 0.0 {
        return Err(reader_error(
            &reader,
            "KR0202",
            "requires one range-independent profile at 0 km",
            "profile_ranges_km",
        ));
    }
    let ranges_km = read_vector(&mut reader, "receiver_ranges_km", false)?;
    let source_depths = read_vector(&mut reader, "field_source_depths_m", true)?;
    let receiver_depths = read_vector(&mut reader, "field_receiver_depths_m", true)?;
    let receiver_offsets = read_vector(&mut reader, "receiver_offsets_m", false)?;
    reader.finish()?;

    Ok(Field {
        mode_limit,
        source_geometry: if option(0) == 'X' {
            SourceGeometry::Line
        } else {
            SourceGeometry::Point
        },
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
    use super::{parse_case, parse_environment, parse_field, parse_frequency_cases, read_file};
    use std::path::Path;

    #[test]
    fn frequency_source_snapshots_match_files_and_enforce_size_limits() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/PekerisBroadband");
        let env = root.with_extension("env");
        let flp = root.with_extension("flp");
        let env_source = read_file(&env).unwrap();
        let flp_source = read_file(&flp).unwrap();
        let parse = |source: &str| {
            super::load_frequency_cases_from_sources(
                source,
                &flp_source,
                Path::new("snapshot.env"),
                Path::new("snapshot.flp"),
                crate::ModeSolver::Kraken,
            )
        };
        assert_eq!(
            parse(&env_source).unwrap(),
            super::load_frequency_cases(&env, &flp, crate::ModeSolver::Kraken).unwrap()
        );
        let report =
            parse(&" ".repeat(usize::try_from(super::MAX_INPUT_BYTES + 1).unwrap())).unwrap_err();
        assert_eq!(report.diagnostics()[0].path, Path::new("snapshot.env"));
        assert!(report.diagnostics()[0].message.contains("1 MiB"));
        let report = super::load_frequency_cases_from_sources(
            &env_source,
            &" ".repeat(usize::try_from(super::MAX_INPUT_BYTES + 1).unwrap()),
            &env,
            &flp,
            crate::ModeSolver::Kraken,
        )
        .unwrap_err();
        assert_eq!(report.diagnostics()[0].path, flp);
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn broadband_preserves_order_and_scales_before_mesh_rounding() {
        let env = include_str!("../tests/fixtures/PekerisBroadband.env");
        let flp = include_str!("../tests/fixtures/PekerisBroadband.flp");
        let parse = |source: &str| {
            parse_frequency_cases(
                source,
                flp,
                Path::new("case.env"),
                Path::new("case.flp"),
                crate::ModeSolver::Kraken,
                false,
            )
            .unwrap()
        };
        let cases = parse(env);
        assert_eq!(
            cases.iter().map(|c| c.frequency_hz).collect::<Vec<_>>(),
            [75.0, 50.0, 62.5]
        );
        assert_eq!(cases[0].mesh_reference_frequency_hz, Some(50.0));
        for (case, expected) in
            cases
                .iter()
                .zip([[151, 303, 606], [101, 202, 404], [126, 252, 505]])
        {
            for (multiplier, n) in [1, 2, 4].into_iter().zip(expected) {
                assert_eq!(case.mesh_points_at(multiplier).unwrap(), n);
            }
        }
        // Resolve NG=0 at freq0 (66), not at the current frequency (100).
        let automatic = parse(&env.replace("101 0.0", "0 0.0"));
        assert_eq!(automatic[0].mesh_points_at(1).unwrap(), 99);
        for frequency in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut input = cases[0].clone().into_definition();
            input.mesh_reference_frequency_hz = Some(frequency);
            let report = crate::Case::from_definition(input).unwrap_err();
            assert_eq!(report.diagnostics()[0].field, "mesh_reference_frequency_hz");
        }
        let ordered = parse(&env.replace("75.0 50.0 62.5 /", "100.0 50.0 /"));
        assert_eq!(
            ordered.iter().map(|c| c.frequency_hz).collect::<Vec<_>>(),
            [100.0, 75.0, 50.0]
        );
        assert_eq!(
            parse(&env.replace("75.0 50.0 62.5 /", "75.0 75.0 50.0 /"))[1].frequency_hz,
            75.0
        );
        let report =
            parse_case(env, flp, Path::new("case.env"), Path::new("case.flp")).unwrap_err();
        assert_eq!(report.diagnostics()[0].field, "frequencies_hz");
        assert_eq!(report.diagnostics()[0].line, 17);
    }

    #[test]
    fn broadband_rejects_bad_frequencies_and_bounded_storage() {
        let env = include_str!("../tests/fixtures/PekerisBroadband.env");
        let flp = include_str!("../tests/fixtures/PekerisBroadband.flp");
        for source in [
            env.replace("75.0 50.0 62.5 /", "0.0 50.0 62.5 /"),
            env.replace("75.0 50.0 62.5 /", "NaN 50.0 62.5 /"),
            env.replace("75.0 50.0 62.5 /", "-75.0 50.0 62.5 /"),
            env.replace("3\n75.0", "1001\n75.0"),
            env.replace("75.0 50.0 62.5 /", "75.0 50.0\n"),
            env.replace("50.0\n1\n", "0.0\n1\n"),
            env.replace("75.0 50.0 62.5 /", "0.01 50.0 62.5 /"),
        ] {
            assert!(
                parse_frequency_cases(
                    &source,
                    flp,
                    Path::new("case.env"),
                    Path::new("case.flp"),
                    crate::ModeSolver::Kraken,
                    false
                )
                .is_err(),
                "accepted {source}"
            );
        }
        let source = env.replace("3\n75.0 50.0 62.5 /", "100\n50.0 100.0 /");
        let field = flp.replace("3\n0.5 1.0 2.0 /", "100000\n0.5 2.0 /");
        let report = parse_frequency_cases(
            &source,
            &field,
            Path::new("case.env"),
            Path::new("case.flp"),
            crate::ModeSolver::Kraken,
            false,
        )
        .unwrap_err();
        assert!(
            report.diagnostics()[0]
                .message
                .contains("input storage limit")
        );
    }

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
    fn subtabulated_depths_use_reference_precision_and_keep_fractional_boundary_valid() {
        let flp = include_str!("../tests/fixtures/Pekeris.flp")
            .replace("3\n25.0 75.0 99.0 /", "10\n0.0 100.0 /")
            .replace("3\n0.0 0.0 0.0 /", "10\n0.0 /");
        let field = parse_field(&flp, Path::new("case.flp")).unwrap();
        let expected: Vec<_> = (0..10)
            .map(|i| f64::from(f32::from(u16::try_from(i).unwrap()) * (100.0_f32 / 9.0)))
            .collect();
        assert_eq!(field.receiver_depths, expected);
        assert_eq!(field.receiver_offsets_m, vec![0.0; 10]);

        let env = "'fractional boundary'\n50\n1\n'NVN'\n100 0 0.1\n0 1500 /\n0.1 1500 /\n'A' 0\n0.1 1700 /\n1400 1700\n0\n1\n0.1 /\n2\n0 0.1 /\n";
        let flp = "/\n'X OC'\n9999\n1\n0 /\n1\n1 /\n1\n0.1 /\n2\n0 0.1 /\n2\n0 0 /\n";
        let case = parse_case(env, flp, Path::new("case.env"), Path::new("case.flp")).unwrap();
        assert_eq!(case.mode_sample_depths_m, [0.0, 0.1]);
        assert_eq!(case.source_depths_m, [0.1]);
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
        for vector in ["3*0.5 /", ",0.5,1.0 /", "0.5,,1.0 /", "0.5 NaN 2.0 /"] {
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
    fn original_sduct_leaky_options_remain_unsupported() {
        let env = include_str!("../tests/fixtures/SductTrapped.env")
            .replace("'CVN'", "'CVW'")
            .replace("1450.0 1523.9", "1450.0 100000");
        assert_eq!(
            parse_case(
                &env,
                include_str!("../tests/fixtures/SductTrapped.flp"),
                Path::new("sductK.env"),
                Path::new("sductK.flp"),
            )
            .unwrap_err()
            .diagnostics()[0]
                .field,
            "phase_speed_limits"
        );
        let env = include_str!("../tests/fixtures/SductTrapped.env")
            .replace("1450.0 1523.9", "1450.0 100000");
        let report = parse_case(
            &env,
            include_str!("../tests/fixtures/SductTrapped.flp"),
            Path::new("sductK.env"),
            Path::new("sductK.flp"),
        )
        .unwrap_err();
        assert!(
            report
                .diagnostics()
                .iter()
                .any(|d| d.field == "phase_speed_limits")
        );
    }

    #[test]
    fn spline_with_no_trapped_water_is_rejected_at_legacy_boundary() {
        let env = include_str!("../tests/fixtures/Pekeris.env")
            .replace("'NVN'", "'SVN'")
            .replace(
                "100.0 1700.0 0.0 1.5 0.0 0.0 /",
                "100.0 1400.0 0.0 1.5 0.0 0.0 /",
            )
            .replace("1400.0 1700.0", "1300.0 1400.0");
        let report = parse_case(
            &env,
            include_str!("../tests/fixtures/Pekeris.flp"),
            Path::new("Pekeris.env"),
            Path::new("Pekeris.flp"),
        )
        .unwrap_err();
        assert!(report.diagnostics().iter().any(|d| {
            d.field == "bottom_sound_speed_mps" && d.path == Path::new("Pekeris.env")
        }));
    }

    #[test]
    fn rigid_bottom_has_no_half_space_record_or_loss_unit() {
        let env = include_str!("../tests/fixtures/PekerisHard.env");
        let unexpected = env.replace("'R' 0.0\n", "'R' 0.0\n100.0 1700.0 0.0 1.5 0.0 0.0 /\n");
        assert_eq!(
            parse_environment(&unexpected, Path::new("case.env"))
                .err()
                .unwrap()
                .diagnostics()[0]
                .field,
            "phase_speed_limits"
        );
        let wrong_unit = env.replace("'SVN'", "'SVW'");
        assert_eq!(
            parse_environment(&wrong_unit, Path::new("case.env"))
                .err()
                .unwrap()
                .diagnostics()[0]
                .field,
            "bottom_options"
        );
    }

    #[test]
    fn unsupported_solver_and_field_options_are_rejected() {
        let env = include_str!("../tests/fixtures/Pekeris.env").replace("'NVN'", "'NVM'");
        let error = parse_environment(&env, Path::new("Pekeris.env"))
            .err()
            .expect("unsupported attenuation unit should be rejected");
        assert_eq!(error.diagnostics()[0].field, "top_options");
        assert_eq!(error.diagnostics()[0].line, 4);
        for option in ["'NAN'", "'NFN'", "'NPN'"] {
            let env = include_str!("../tests/fixtures/Pekeris.env").replace("'NVN'", option);
            assert_eq!(
                parse_environment(&env, Path::new("Pekeris.env"))
                    .err()
                    .unwrap()
                    .diagnostics()[0]
                    .field,
                "top_options"
            );
        }
        let analytic = include_str!("../tests/fixtures/Pekeris.env").replace("'NVN'", "'AVN'");
        assert_eq!(
            parse_environment(&analytic, Path::new("Pekeris.env"))
                .err()
                .unwrap()
                .diagnostics()[0]
                .field,
            "bottom_options"
        );

        let analytic = include_str!("../tests/fixtures/MunkAnalytic.env");
        let extra_point = analytic.replace("'A' 0.0", "0 1500 /\n'A' 0.0");
        assert_eq!(
            parse_environment(&extra_point, Path::new("MunkAnalytic.env"))
                .err()
                .unwrap()
                .diagnostics()[0]
                .field,
            "bottom_options"
        );
        let missing_material = analytic.replace("5000.0  1600.00 0.0 1.8 0.8 /", "5000 1600 /");
        assert_eq!(
            parse_environment(&missing_material, Path::new("MunkAnalytic.env"))
                .err()
                .unwrap()
                .diagnostics()[0]
                .field,
            "bottom_half_space"
        );

        let lossy_water = include_str!("../tests/fixtures/MunkBottomLoss.env")
            .replace("200.0 1530.29 /", "200.0 1530.29 0.0 1.0 0.1 /");
        assert_eq!(
            parse_environment(&lossy_water, Path::new("MunkBottomLoss.env"))
                .err()
                .unwrap()
                .diagnostics()[0]
                .field,
            "sound_speed_profile"
        );
        let unsupported_bottom =
            include_str!("../tests/fixtures/Pekeris.env").replace("'A' 0.0", "'V' 0.0");
        assert_eq!(
            parse_environment(&unsupported_bottom, Path::new("Pekeris.env"))
                .err()
                .unwrap()
                .diagnostics()[0]
                .field,
            "bottom_options"
        );
        let wrong_units =
            include_str!("../tests/fixtures/MunkBottomLoss.env").replace("'NVW'", "'NVN'");
        assert_eq!(
            parse_environment(&wrong_units, Path::new("MunkBottomLoss.env"))
                .err()
                .unwrap()
                .diagnostics()[0]
                .field,
            "bottom_half_space"
        );

        let flp = include_str!("../tests/fixtures/Pekeris.flp").replace("'X OC'", "'S OC'");
        let error = parse_field(&flp, Path::new("Pekeris.flp"))
            .err()
            .expect("scaled source should be rejected in this slice");
        assert_eq!(error.diagnostics()[0].field, "field_options");
        assert_eq!(error.diagnostics()[0].line, 2);
    }
}
