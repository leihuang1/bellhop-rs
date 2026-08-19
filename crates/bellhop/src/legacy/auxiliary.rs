#![allow(clippy::cast_sign_loss, clippy::too_many_lines)]

use std::f64::consts::PI;
use std::path::Path;

use num_complex::Complex64;

use crate::diagnostic::{Diagnostic, SourceLocation};
use crate::model::{
    BoundaryInterpolation, BoundaryMaterial, BoundaryShape, BoundaryShapePoint,
    InternalReflectionCoefficientPoint, InternalReflectionCoefficientTable,
    RangeDependentSoundSpeed, ReflectionCoefficientPoint, ReflectionCoefficientTable,
    SourceBeamPattern, SourceBeamPatternPoint,
};

use super::records::{Atom, RecordReader, Slot, parse_f64, parse_i32};

const MAX_AUXILIARY_POINTS: usize = 1_000_000;
const MAX_SSP_VALUES: usize = 20_000_000;

pub(super) fn parse_range_dependent_sound_speed(
    source: &str,
    path: &Path,
    depths_m: &[f64],
) -> Result<RangeDependentSoundSpeed, Diagnostic> {
    if depths_m.len() < 2 {
        return Err(Diagnostic::error(
            "BH0201",
            "range-dependent sound speed requires at least two profile depths",
            "sound_speed.points",
            SourceLocation::file(path),
        ));
    }

    let mut reader = RecordReader::new(source, path);
    let (range_count, count_location) =
        reader.read_i32("range_dependent_sound_speed.range_count")?;
    let range_count = checked_count(
        range_count,
        2,
        "range_dependent_sound_speed.range_count",
        count_location,
    )?;
    let value_count = range_count.checked_mul(depths_m.len()).ok_or_else(|| {
        resource_error(
            "range-dependent sound-speed matrix is too large",
            "range_dependent_sound_speed.speeds",
            path,
        )
    })?;
    if value_count > MAX_SSP_VALUES {
        return Err(resource_error(
            format!(
                "range-dependent sound-speed matrix has {value_count} values; limit is {MAX_SSP_VALUES}"
            ),
            "range_dependent_sound_speed.speeds",
            path,
        ));
    }

    let range_slots = reader.read_fields("range_dependent_sound_speed.ranges", range_count)?;
    let range_atoms = required_atoms(
        range_slots,
        range_count,
        "range_dependent_sound_speed.ranges",
        path,
    )?;
    let mut ranges_m = Vec::with_capacity(range_count);
    let mut previous_range = None;
    for atom in &range_atoms {
        let profile_range_m = 1000.0 * parse_f64(atom, "range_dependent_sound_speed.ranges")?;
        require_finite(
            profile_range_m,
            "range-dependent sound-speed ranges must be finite",
            "range_dependent_sound_speed.ranges",
            &atom.location,
        )?;
        if previous_range.is_some_and(|previous| profile_range_m <= previous) {
            return Err(Diagnostic::error(
                "BH0201",
                "range-dependent sound-speed ranges must be strictly increasing",
                "range_dependent_sound_speed.ranges",
                atom.location.clone(),
            ));
        }
        previous_range = Some(profile_range_m);
        ranges_m.push(profile_range_m);
    }

    let mut speeds_mps = Vec::with_capacity(depths_m.len());
    for _ in depths_m {
        let slots = reader.read_fields("range_dependent_sound_speed.speeds", range_count)?;
        let atoms = required_atoms(
            slots,
            range_count,
            "range_dependent_sound_speed.speeds",
            path,
        )?;
        let mut row = Vec::with_capacity(range_count);
        for atom in atoms {
            let speed = parse_f64(&atom, "range_dependent_sound_speed.speeds")?;
            if !speed.is_finite() || speed <= 0.0 {
                return Err(Diagnostic::error(
                    "BH0201",
                    "sound speeds must be finite and positive",
                    "range_dependent_sound_speed.speeds",
                    atom.location,
                ));
            }
            row.push(speed);
        }
        speeds_mps.push(row);
    }

    Ok(RangeDependentSoundSpeed {
        ranges_m,
        depths_m: depths_m.to_vec(),
        speeds_mps,
    })
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum BoundarySide {
    Top,
    Bottom,
}

impl BoundarySide {
    fn name(self) -> &'static str {
        match self {
            Self::Top => "altimetry",
            Self::Bottom => "bathymetry",
        }
    }

    fn depth_error(self) -> &'static str {
        match self {
            Self::Top => "altimetry rises above the highest sound-speed profile depth",
            Self::Bottom => "bathymetry drops below the lowest sound-speed profile depth",
        }
    }
}

pub(super) fn parse_boundary_shape(
    source: &str,
    path: &Path,
    side: BoundarySide,
    profile_depth_m: f64,
) -> Result<BoundaryShape, Diagnostic> {
    let field = side.name();
    let mut reader = RecordReader::new(source, path);
    let option = reader.read_string(field)?;
    let mut characters = option.text.chars();
    let interpolation = match characters.next().unwrap_or(' ') {
        'L' => BoundaryInterpolation::PiecewiseLinear,
        'C' => BoundaryInterpolation::Curvilinear,
        other => {
            return Err(invalid_option(
                format!("unknown boundary interpolation option {other:?}; expected L or C"),
                field,
                option.location,
            ));
        }
    };
    let long_format = match characters.next().unwrap_or(' ') {
        ' ' | 'S' => false,
        'L' => true,
        other => {
            return Err(invalid_option(
                format!("unknown boundary point format {other:?}; expected blank, S, or L"),
                field,
                option.location,
            ));
        }
    };

    let (point_count, count_location) = reader.read_i32(field)?;
    let point_count = checked_count(point_count, 1, field, count_location)?;
    let values_per_point = if long_format { 7 } else { 2 };
    let mut points = Vec::with_capacity(point_count);
    let mut previous_range = None;

    for _ in 0..point_count {
        let atoms = required_atoms(
            reader.read_fields(field, values_per_point)?,
            values_per_point,
            field,
            path,
        )?;
        let range_m = 1000.0 * parse_f64(&atoms[0], field)?;
        let depth_m = parse_f64(&atoms[1], field)?;
        require_finite(
            range_m,
            "boundary ranges must be finite",
            field,
            &atoms[0].location,
        )?;
        require_finite(
            depth_m,
            "boundary depths must be finite",
            field,
            &atoms[1].location,
        )?;
        if previous_range.is_some_and(|previous| range_m <= previous) {
            return Err(Diagnostic::error(
                "BH0201",
                "boundary ranges must be strictly increasing",
                field,
                atoms[0].location.clone(),
            ));
        }
        previous_range = Some(range_m);
        let outside_profile = match side {
            BoundarySide::Top => depth_m < profile_depth_m,
            BoundarySide::Bottom => depth_m > profile_depth_m,
        };
        if outside_profile {
            return Err(Diagnostic::error(
                "BH0201",
                side.depth_error(),
                field,
                atoms[1].location.clone(),
            ));
        }

        let material = if long_format {
            let material = BoundaryMaterial {
                compressional_speed_mps: parse_f64(&atoms[2], field)?,
                shear_speed_mps: parse_f64(&atoms[3], field)?,
                density_g_cm3: parse_f64(&atoms[4], field)?,
                compressional_attenuation: parse_f64(&atoms[5], field)?,
                shear_attenuation: parse_f64(&atoms[6], field)?,
            };
            if !material.compressional_speed_mps.is_finite()
                || material.compressional_speed_mps <= 0.0
                || !material.shear_speed_mps.is_finite()
                || material.shear_speed_mps < 0.0
                || !material.density_g_cm3.is_finite()
                || material.density_g_cm3 <= 0.0
                || !material.compressional_attenuation.is_finite()
                || !material.shear_attenuation.is_finite()
            {
                return Err(Diagnostic::error(
                    "BH0201",
                    "boundary material properties contain invalid physical values",
                    field,
                    atoms[2].location.clone(),
                ));
            }
            Some(material)
        } else {
            None
        };

        points.push(BoundaryShapePoint {
            range_m,
            depth_m,
            material,
        });
    }

    Ok(BoundaryShape {
        interpolation,
        points,
    })
}

pub(super) fn parse_reflection_coefficients(
    source: &str,
    path: &Path,
    field: &'static str,
) -> Result<ReflectionCoefficientTable, Diagnostic> {
    let mut reader = RecordReader::new(source, path);
    let (point_count, count_location) = reader.read_i32(field)?;
    let point_count = checked_count(point_count, 2, field, count_location)?;
    let value_count = point_count
        .checked_mul(3)
        .ok_or_else(|| resource_error("reflection coefficient table is too large", field, path))?;
    let atoms = required_atoms(
        reader.read_fields(field, value_count)?,
        value_count,
        field,
        path,
    )?;
    let mut points = Vec::with_capacity(point_count);
    let mut previous_angle = None;
    for values in atoms.chunks_exact(3) {
        let angle_degrees = parse_f64(&values[0], field)?;
        let magnitude = parse_f64(&values[1], field)?;
        let phase_degrees = parse_f64(&values[2], field)?;
        if !angle_degrees.is_finite()
            || !magnitude.is_finite()
            || magnitude < 0.0
            || !phase_degrees.is_finite()
        {
            return Err(Diagnostic::error(
                "BH0201",
                "reflection coefficient values must be finite and magnitudes non-negative",
                field,
                values[0].location.clone(),
            ));
        }
        if previous_angle.is_some_and(|previous| angle_degrees <= previous) {
            return Err(Diagnostic::error(
                "BH0201",
                "reflection coefficient angles must be strictly increasing",
                field,
                values[0].location.clone(),
            ));
        }
        previous_angle = Some(angle_degrees);
        points.push(ReflectionCoefficientPoint {
            angle_degrees,
            magnitude,
            phase_radians: phase_degrees * PI / 180.0,
        });
    }
    Ok(ReflectionCoefficientTable { points })
}

pub(super) fn parse_internal_reflection_coefficients(
    source: &str,
    path: &Path,
) -> Result<InternalReflectionCoefficientTable, Diagnostic> {
    const FIELD: &str = "internal_reflection";
    let mut reader = RecordReader::new(source, path);
    let header = required_atoms(reader.read_fields(FIELD, 2)?, 2, FIELD, path)?;
    let title = header[0].text.clone();
    let frequency_hz = parse_f64(&header[1], FIELD)?;
    if !frequency_hz.is_finite() || frequency_hz <= 0.0 {
        return Err(Diagnostic::error(
            "BH0201",
            "internal reflection-table frequency must be finite and positive",
            FIELD,
            header[1].location.clone(),
        ));
    }

    let (point_count, count_location) = reader.read_i32(FIELD)?;
    let point_count = checked_count(point_count, 2, FIELD, count_location)?;
    let mut points = Vec::with_capacity(point_count);
    let mut previous_wavenumber = None;
    for _ in 0..point_count {
        let values = required_atoms(
            reader.read_fixed_width_fields(FIELD, &[15, 15, 15, 15, 15, 5])?,
            6,
            FIELD,
            path,
        )?;
        let horizontal_wavenumber_squared = parse_f64(&values[0], FIELD)?;
        let f = Complex64::new(parse_f64(&values[1], FIELD)?, parse_f64(&values[2], FIELD)?);
        let g = Complex64::new(parse_f64(&values[3], FIELD)?, parse_f64(&values[4], FIELD)?);
        let decimal_power = parse_i32(&values[5], FIELD)?;
        if !horizontal_wavenumber_squared.is_finite()
            || horizontal_wavenumber_squared < 0.0
            || !f.re.is_finite()
            || !f.im.is_finite()
            || !g.re.is_finite()
            || !g.im.is_finite()
        {
            return Err(Diagnostic::error(
                "BH0201",
                "internal reflection-table values must be finite and squared wavenumbers non-negative",
                FIELD,
                values[0].location.clone(),
            ));
        }
        if previous_wavenumber.is_some_and(|previous| horizontal_wavenumber_squared <= previous) {
            return Err(Diagnostic::error(
                "BH0201",
                "internal reflection-table squared wavenumbers must be strictly increasing",
                FIELD,
                values[0].location.clone(),
            ));
        }
        previous_wavenumber = Some(horizontal_wavenumber_squared);
        points.push(InternalReflectionCoefficientPoint {
            horizontal_wavenumber_squared,
            f,
            g,
            decimal_power,
        });
    }

    Ok(InternalReflectionCoefficientTable {
        title,
        frequency_hz,
        points,
    })
}

pub(super) fn parse_source_beam_pattern(
    source: &str,
    path: &Path,
) -> Result<SourceBeamPattern, Diagnostic> {
    const FIELD: &str = "source_beam_pattern";
    let mut reader = RecordReader::new(source, path);
    let (point_count, count_location) = reader.read_i32(FIELD)?;
    let point_count = checked_count(point_count, 2, FIELD, count_location)?;
    let mut points = Vec::with_capacity(point_count);
    let mut previous_angle = None;
    for _ in 0..point_count {
        let atoms = required_atoms(reader.read_fields(FIELD, 2)?, 2, FIELD, path)?;
        let angle_degrees = parse_f64(&atoms[0], FIELD)?;
        let level_db = parse_f64(&atoms[1], FIELD)?;
        if !angle_degrees.is_finite() || !level_db.is_finite() {
            return Err(Diagnostic::error(
                "BH0201",
                "source beam-pattern values must be finite",
                FIELD,
                atoms[0].location.clone(),
            ));
        }
        if previous_angle.is_some_and(|previous| angle_degrees <= previous) {
            return Err(Diagnostic::error(
                "BH0201",
                "source beam-pattern angles must be strictly increasing",
                FIELD,
                atoms[0].location.clone(),
            ));
        }
        previous_angle = Some(angle_degrees);
        points.push(SourceBeamPatternPoint {
            angle_degrees,
            level_db,
            amplitude: 10.0_f64.powf(level_db / 20.0),
        });
    }
    Ok(SourceBeamPattern { points })
}

fn checked_count(
    count: i32,
    minimum: usize,
    field: &'static str,
    location: SourceLocation,
) -> Result<usize, Diagnostic> {
    let count = usize::try_from(count).map_err(|_| {
        Diagnostic::error(
            "BH0201",
            "point count cannot be negative",
            field,
            location.clone(),
        )
    })?;
    if count < minimum {
        return Err(Diagnostic::error(
            "BH0201",
            format!("point count must be at least {minimum}"),
            field,
            location,
        ));
    }
    if count > MAX_AUXILIARY_POINTS {
        return Err(Diagnostic::error(
            "BH0203",
            format!("point count {count} exceeds limit {MAX_AUXILIARY_POINTS}"),
            field,
            location,
        ));
    }
    Ok(count)
}

fn required_atoms(
    slots: Vec<Slot>,
    expected: usize,
    field: &'static str,
    path: &Path,
) -> Result<Vec<Atom>, Diagnostic> {
    if slots.len() < expected || slots.iter().any(|slot| slot.atom.is_none()) {
        return Err(Diagnostic::error(
            "BH0104",
            format!("{expected} value(s) are required"),
            field,
            SourceLocation::file(path),
        ));
    }
    Ok(slots
        .into_iter()
        .map(|slot| slot.atom.expect("checked above"))
        .collect())
}

fn require_finite(
    value: f64,
    message: &'static str,
    field: &'static str,
    location: &SourceLocation,
) -> Result<(), Diagnostic> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(Diagnostic::error(
            "BH0201",
            message,
            field,
            location.clone(),
        ))
    }
}

fn invalid_option(
    message: impl Into<String>,
    field: &'static str,
    location: SourceLocation,
) -> Diagnostic {
    Diagnostic::error("BH0105", message, field, location)
}

fn resource_error(message: impl Into<String>, field: &'static str, path: &Path) -> Diagnostic {
    Diagnostic::error("BH0203", message, field, SourceLocation::file(path))
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;
    use std::path::Path;

    use crate::model::BoundaryInterpolation;

    use super::{
        BoundarySide, parse_boundary_shape, parse_internal_reflection_coefficients,
        parse_range_dependent_sound_speed, parse_reflection_coefficients,
        parse_source_beam_pattern,
    };

    #[test]
    fn parses_range_dependent_sound_speed() {
        let input = "3\n0 1 2\n1500 1501 1502\n1510 1511 1512\n";
        let field =
            parse_range_dependent_sound_speed(input, Path::new("case.ssp"), &[0.0, 100.0]).unwrap();
        assert_eq!(field.ranges_m, vec![0.0, 1000.0, 2000.0]);
        assert_eq!(field.speeds_mps[1], vec![1510.0, 1511.0, 1512.0]);
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn parses_short_and_long_boundary_shapes() {
        let short = parse_boundary_shape(
            "'C'\n2\n0 0\n1 10\n",
            Path::new("case.ati"),
            BoundarySide::Top,
            0.0,
        )
        .unwrap();
        assert_eq!(short.interpolation, BoundaryInterpolation::Curvilinear);
        assert_eq!(short.points[1].range_m, 1000.0);
        assert!(short.points[0].material.is_none());

        let long = parse_boundary_shape(
            "'LL'\n2\n0 100 1700 0 1.2 0.5 0\n1 90 1600 0 1.1 0.2 0\n",
            Path::new("case.bty"),
            BoundarySide::Bottom,
            100.0,
        )
        .unwrap();
        assert_eq!(long.points[0].material.as_ref().unwrap().density_g_cm3, 1.2);
    }

    #[test]
    fn parses_reflection_coefficients_and_converts_phase() {
        let table = parse_reflection_coefficients(
            "2\n0 1 0\n90 0.5 180\n",
            Path::new("case.brc"),
            "bottom_reflection",
        )
        .unwrap();
        assert!((table.points[1].phase_radians - PI).abs() < 1.0e-15);
    }

    #[test]
    fn parses_internal_reflection_impedance_functions() {
        let table = parse_internal_reflection_coefficients(
            "'generated table' 500\n2\n1 2 3 4 5 -10\n2 6 7 8 9 0\n",
            Path::new("case.irc"),
        )
        .unwrap();
        assert_eq!(table.title, "generated table");
        assert!((table.frequency_hz - 500.0).abs() < f64::EPSILON);
        assert_eq!(table.points[0].f, num_complex::Complex64::new(2.0, 3.0));
        assert_eq!(table.points[1].g, num_complex::Complex64::new(8.0, 9.0));
        assert_eq!(table.points[0].decimal_power, -10);
    }

    #[test]
    fn parses_source_pattern_and_converts_db() {
        let pattern =
            parse_source_beam_pattern("2\n-180 0\n180 -20\n", Path::new("case.sbp")).unwrap();
        assert!((pattern.points[1].amplitude - 0.1).abs() < 1.0e-15);
    }
}
