use std::f64::consts::PI;

use num_complex::Complex64;

use crate::model::{
    AttenuationUnit, BoundaryCondition, BoundaryMaterial, CurvatureCondition, EnvironmentCase,
    HalfSpace, InternalReflectionCoefficientTable, ReflectionCoefficientTable, VolumeAttenuation,
};

use super::boundary::{BoundarySegment, BoundarySide, add, dot, scale, subtract};
use super::integrator::RayState;
use super::ssp::{SegmentState, SoundSpeedModel, material_complex_speed};

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(crate) fn reflect_2d(
    incident: RayState,
    side: BoundarySide,
    boundary_segment: &BoundarySegment,
    boundary_tangent: [f64; 2],
    boundary_normal: [f64; 2],
    condition: &BoundaryCondition,
    reflection_table: Option<&ReflectionCoefficientTable>,
    internal_reflection: Option<&InternalReflectionCoefficientTable>,
    environment: &EnvironmentCase,
    sound_speed: &SoundSpeedModel,
    sound_speed_segments: &mut SegmentState,
    beam_shift: bool,
) -> Result<RayState, &'static str> {
    let tangent_component = dot(incident.tangent, boundary_tangent);
    let normal_component = dot(incident.tangent, boundary_normal);
    let mut reflected = RayState {
        tangent: subtract(
            incident.tangent,
            scale(boundary_normal, 2.0 * normal_component),
        ),
        ..incident
    };

    let sample = sound_speed.evaluate(incident.position_m, sound_speed_segments)?;
    let incident_ray_tangent = scale(incident.tangent, sample.speed_mps);
    let incident_ray_normal = [-incident_ray_tangent[1], incident_ray_tangent[0]];
    let reflected_ray_tangent = scale(reflected.tangent, sample.speed_mps);
    let reflected_ray_normal = [reflected_ray_tangent[1], -reflected_ray_tangent[0]];
    let mut curvature_change =
        2.0 * boundary_segment.curvature / sample.speed_mps.powi(2) / normal_component;
    let mut normal_gradient_jump = -dot(
        sample.gradient,
        subtract(reflected_ray_normal, incident_ray_normal),
    );
    let tangent_gradient_jump = -dot(
        sample.gradient,
        subtract(reflected_ray_tangent, incident_ray_tangent),
    );
    if side == BoundarySide::Top {
        normal_gradient_jump = -normal_gradient_jump;
        curvature_change = -curvature_change;
    }
    let incidence_tangent = tangent_component / normal_component;
    curvature_change += incidence_tangent
        * (2.0 * normal_gradient_jump - incidence_tangent * tangent_gradient_jump)
        / sample.speed_mps.powi(2);
    match environment
        .trace
        .cerveny
        .as_ref()
        .map(|options| options.curvature)
    {
        Some(CurvatureCondition::Double) => curvature_change *= 2.0,
        Some(CurvatureCondition::Zero) => curvature_change = 0.0,
        Some(CurvatureCondition::Standard) | None => {}
    }
    reflected.speed_mps = sample.speed_mps;
    reflected.p = add(incident.p, scale(incident.q, curvature_change));
    reflected.q = incident.q;

    match condition {
        BoundaryCondition::Rigid => {}
        BoundaryCondition::Vacuum => reflected.phase_radians += PI,
        BoundaryCondition::ReflectionCoefficientFile => {
            let table = reflection_table.ok_or("reflection coefficient table is not loaded")?;
            let mut incidence_degrees =
                normal_component.atan2(tangent_component).abs().to_degrees();
            if incidence_degrees > 90.0 {
                incidence_degrees = 180.0 - incidence_degrees;
            }
            let (magnitude, phase_radians) = interpolate_reflection(table, incidence_degrees);
            reflected.amplitude *= magnitude;
            reflected.phase_radians += phase_radians;
        }
        BoundaryCondition::AcoustoElastic(half_space) => {
            let properties = half_space_properties(
                environment,
                half_space,
                boundary_segment.material.as_ref(),
                incident.position_m[1],
            );
            apply_half_space_reflection(
                &mut reflected,
                incident,
                sample.density_g_cm3,
                properties,
                tangent_component,
                normal_component,
                environment.frequency_hz,
                boundary_tangent,
                beam_shift,
            );
        }
        BoundaryCondition::GrainSize { phi, .. } => {
            let properties = boundary_segment.material.as_ref().map_or_else(
                || grain_properties(environment, *phi),
                |material| {
                    segment_material_properties(environment, material, incident.position_m[1])
                },
            );
            apply_half_space_reflection(
                &mut reflected,
                incident,
                sample.density_g_cm3,
                properties,
                tangent_component,
                normal_component,
                environment.frequency_hz,
                boundary_tangent,
                beam_shift,
            );
        }
        BoundaryCondition::WriteReflectionCoefficient => {
            return Err("writing internal reflection coefficients is not supported");
        }
        BoundaryCondition::PrecalculatedReflectionCoefficient => {
            let table = internal_reflection.ok_or("internal reflection table is not loaded")?;
            let angular_frequency = 2.0 * PI * environment.frequency_hz;
            let horizontal_wavenumber_squared =
                Complex64::new((angular_frequency * tangent_component).powi(2), 0.0);
            let (f, g) = interpolate_internal_reflection(table, horizontal_wavenumber_squared)?;
            let imaginary_normal = Complex64::new(0.0, angular_frequency * normal_component);
            let coefficient = -(f - imaginary_normal * g) / (f + imaginary_normal * g);
            apply_complex_reflection(&mut reflected, incident, coefficient);
        }
    }

    match side {
        BoundarySide::Top => reflected.top_bounces += 1,
        BoundarySide::Bottom => reflected.bottom_bounces += 1,
    }
    if reflected_state_is_finite(reflected) {
        Ok(reflected)
    } else {
        Err("boundary reflection produced a non-finite ray state")
    }
}

#[derive(Clone, Copy, Debug)]
struct HalfSpaceProperties {
    compressional_speed: Complex64,
    shear_speed: Complex64,
    density_g_cm3: f64,
}

fn half_space_properties(
    environment: &EnvironmentCase,
    half_space: &HalfSpace,
    segment_material: Option<&BoundaryMaterial>,
    reflection_depth_m: f64,
) -> HalfSpaceProperties {
    if let Some(material) = segment_material {
        segment_material_properties(environment, material, reflection_depth_m)
    } else {
        HalfSpaceProperties {
            compressional_speed: material_complex_speed(
                half_space.depth_m,
                half_space.compressional_speed_mps,
                half_space.compressional_attenuation,
                environment.frequency_hz,
                environment.top_options.attenuation_unit,
                &environment.top_options.volume_attenuation,
            ),
            shear_speed: material_complex_speed(
                half_space.depth_m,
                half_space.shear_speed_mps,
                half_space.shear_attenuation,
                environment.frequency_hz,
                environment.top_options.attenuation_unit,
                &environment.top_options.volume_attenuation,
            ),
            density_g_cm3: half_space.density_g_cm3,
        }
    }
}

fn segment_material_properties(
    environment: &EnvironmentCase,
    material: &BoundaryMaterial,
    reflection_depth_m: f64,
) -> HalfSpaceProperties {
    HalfSpaceProperties {
        compressional_speed: material_complex_speed(
            reflection_depth_m,
            material.compressional_speed_mps,
            material.compressional_attenuation,
            environment.frequency_hz,
            AttenuationUnit::DbPerWavelength,
            &VolumeAttenuation::None,
        ),
        shear_speed: material_complex_speed(
            reflection_depth_m,
            material.shear_speed_mps,
            material.shear_attenuation,
            environment.frequency_hz,
            AttenuationUnit::DbPerWavelength,
            &VolumeAttenuation::None,
        ),
        density_g_cm3: material.density_g_cm3,
    }
}

fn grain_properties(environment: &EnvironmentCase, grain_size_phi: f64) -> HalfSpaceProperties {
    let speed_ratio = if (-1.0..1.0).contains(&grain_size_phi) {
        0.002_709 * grain_size_phi.powi(2) - 0.056_452 * grain_size_phi + 1.2778
    } else if (1.0..5.3).contains(&grain_size_phi) {
        -0.001_488_1 * grain_size_phi.powi(3) + 0.021_393_7 * grain_size_phi.powi(2)
            - 0.138_279_8 * grain_size_phi
            + 1.3425
    } else {
        -0.002_432_4 * grain_size_phi + 1.0019
    };
    let loss = if (-1.0..0.0).contains(&grain_size_phi) {
        0.4556
    } else if (0.0..2.6).contains(&grain_size_phi) {
        0.4556 + 0.0245 * grain_size_phi
    } else if (2.6..4.5).contains(&grain_size_phi) {
        0.1978 + 0.1245 * grain_size_phi
    } else if (4.5..6.0).contains(&grain_size_phi) {
        8.0399 - 2.5228 * grain_size_phi + 0.200_98 * grain_size_phi.powi(2)
    } else if (6.0..9.5).contains(&grain_size_phi) {
        0.9431 - 0.2041 * grain_size_phi + 0.0117 * grain_size_phi.powi(2)
    } else {
        0.0601
    };
    let real_speed = speed_ratio * 1500.0;
    let attenuation = loss * (speed_ratio / 1000.0) * 1500.0 * 10.0_f64.ln() / (40.0 * PI);
    let density_g_cm3 = if (-1.0..1.0).contains(&grain_size_phi) {
        0.007_797 * grain_size_phi.powi(2) - 0.170_57 * grain_size_phi + 2.3139
    } else if (1.0..5.3).contains(&grain_size_phi) {
        -0.016_540_6 * grain_size_phi.powi(3) + 0.229_020_1 * grain_size_phi.powi(2)
            - 1.106_903_1 * grain_size_phi
            + 3.0455
    } else {
        -0.001_297_3 * grain_size_phi + 1.1565
    };
    HalfSpaceProperties {
        compressional_speed: material_complex_speed(
            0.0,
            real_speed,
            attenuation,
            environment.frequency_hz,
            AttenuationUnit::LossParameter,
            &VolumeAttenuation::None,
        ),
        shear_speed: Complex64::new(0.0, 0.0),
        density_g_cm3,
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn apply_half_space_reflection(
    reflected: &mut RayState,
    incident: RayState,
    water_density_g_cm3: f64,
    half_space: HalfSpaceProperties,
    tangent_component: f64,
    normal_component: f64,
    frequency_hz: f64,
    boundary_tangent: [f64; 2],
    beam_shift: bool,
) {
    let angular_frequency = 2.0 * PI * frequency_hz;
    let horizontal_wavenumber = Complex64::new(angular_frequency * tangent_component, 0.0);
    let normal_wavenumber = angular_frequency * normal_component;
    let (f, g) = if half_space.shear_speed.re > 0.0 {
        let shear_vertical_squared =
            horizontal_wavenumber.powu(2) - (angular_frequency / half_space.shear_speed).powu(2);
        let compressional_vertical_squared = horizontal_wavenumber.powu(2)
            - (angular_frequency / half_space.compressional_speed).powu(2);
        let shear_vertical = shear_vertical_squared.sqrt();
        let compressional_vertical = compressional_vertical_squared.sqrt();
        let rigidity = half_space.density_g_cm3 * half_space.shear_speed.powu(2);
        let y2 = ((shear_vertical_squared + horizontal_wavenumber.powu(2)).powu(2)
            - 4.0 * shear_vertical * compressional_vertical * horizontal_wavenumber.powu(2))
            * rigidity;
        let y4 = compressional_vertical * (horizontal_wavenumber.powu(2) - shear_vertical_squared);
        (angular_frequency.powi(2) * y4, y2)
    } else {
        let mut compressional_vertical = (horizontal_wavenumber.powu(2)
            - (angular_frequency / half_space.compressional_speed).powu(2))
        .sqrt();
        if compressional_vertical.re == 0.0 && compressional_vertical.im < 0.0 {
            compressional_vertical = -compressional_vertical;
        }
        (
            compressional_vertical,
            Complex64::new(half_space.density_g_cm3, 0.0),
        )
    };
    let imaginary_normal = Complex64::new(0.0, normal_wavenumber);
    let reflection = -(water_density_g_cm3 * f - imaginary_normal * g)
        / (water_density_g_cm3 * f + imaginary_normal * g);
    let reflected_energy = reflection.norm() >= 1.0e-5;
    apply_complex_reflection(reflected, incident, reflection);
    if reflected_energy && beam_shift {
        apply_beam_shift(
            reflected,
            incident,
            half_space,
            angular_frequency,
            boundary_tangent,
        );
    }
}

fn apply_complex_reflection(reflected: &mut RayState, incident: RayState, coefficient: Complex64) {
    if coefficient.norm() < 1.0e-5 {
        reflected.amplitude = 0.0;
    } else {
        reflected.amplitude = incident.amplitude * coefficient.norm();
        reflected.phase_radians = incident.phase_radians + coefficient.arg();
    }
}

fn apply_beam_shift(
    reflected: &mut RayState,
    incident: RayState,
    half_space: HalfSpaceProperties,
    angular_frequency: f64,
    boundary_tangent: [f64; 2],
) {
    let ch = incident.speed_mps / half_space.compressional_speed.conj();
    let cosine = incident.tangent[0] * incident.speed_mps;
    let sine = incident.tangent[1] * incident.speed_mps;
    let wavenumber = angular_frequency / incident.speed_mps;
    let a = 2.0 * half_space.density_g_cm3 * (1.0 - ch.powu(2));
    let b = cosine.powi(2) - ch.powu(2);
    let d = half_space.density_g_cm3.powi(2) * sine.powi(2) + b;
    let square_root_b = b.sqrt();
    let cosine_squared = cosine.powi(2);
    let sine_squared = sine.powi(2);
    let displacement = if sine == 0.0 {
        Complex64::new(0.0, 0.0)
    } else {
        a * cosine / sine / (wavenumber * square_root_b * d)
    };
    let phase_displacement = displacement.re / (incident.speed_mps / cosine);
    let derivative = -a / (wavenumber * square_root_b * d)
        - a * cosine_squared / sine_squared / (wavenumber * square_root_b * d)
        + a * cosine_squared / (wavenumber * b * square_root_b * d)
        - a * cosine / sine / (wavenumber * square_root_b * d.powu(2))
            * (2.0 * half_space.density_g_cm3.powi(2) * sine * cosine - 2.0 * cosine * sine);
    let range_derivative = -derivative.re;
    let derivative_sign = range_derivative / range_derivative.abs();
    reflected.position_m = add(
        reflected.position_m,
        scale(boundary_tangent, displacement.re),
    );
    reflected.travel_time_s += phase_displacement;
    reflected.q = add(
        reflected.q,
        scale(
            incident.p,
            derivative_sign * range_derivative * sine * incident.speed_mps,
        ),
    );
}

fn interpolate_internal_reflection(
    table: &InternalReflectionCoefficientTable,
    x: Complex64,
) -> Result<(Complex64, Complex64), &'static str> {
    let points = &table.points;
    if x.re < points[0].horizontal_wavenumber_squared {
        return Ok((points[0].f, points[0].g));
    }
    if x.re > points[points.len() - 1].horizontal_wavenumber_squared {
        let point = &points[points.len() - 1];
        return Ok((point.f, point.g));
    }

    let right = points.partition_point(|point| point.horizontal_wavenumber_squared <= x.re);
    let left = right.saturating_sub(1).min(points.len() - 2);
    let right = (left + 2).min(points.len() - 1);
    let power = points[left].decimal_power;
    let mut abscissas = Vec::with_capacity(right - left + 1);
    let mut f_values = Vec::with_capacity(right - left + 1);
    let mut g_values = Vec::with_capacity(right - left + 1);
    for point in &points[left..=right] {
        let scale = 10.0_f64.powi(point.decimal_power - power);
        abscissas.push(Complex64::new(point.horizontal_wavenumber_squared, 0.0));
        f_values.push(point.f * scale);
        g_values.push(point.g * scale);
    }
    let f = polynomial(x, &abscissas, &f_values);
    let g = polynomial(x, &abscissas, &g_values);
    if f.re.is_finite() && f.im.is_finite() && g.re.is_finite() && g.im.is_finite() {
        Ok((f, g))
    } else {
        Err("internal reflection interpolation produced a non-finite impedance")
    }
}

fn polynomial(x: Complex64, abscissas: &[Complex64], values: &[Complex64]) -> Complex64 {
    let mut interpolated = values.to_vec();
    let offsets: Vec<Complex64> = abscissas.iter().map(|&value| value - x).collect();
    for order in 1..abscissas.len() {
        for index in 0..abscissas.len() - order {
            interpolated[index] = (offsets[index + order] * interpolated[index]
                - offsets[index] * interpolated[index + 1])
                / (offsets[index + order] - offsets[index]);
        }
    }
    interpolated[0]
}

fn interpolate_reflection(table: &ReflectionCoefficientTable, angle_degrees: f64) -> (f64, f64) {
    let points = &table.points;
    if angle_degrees < points[0].angle_degrees
        || angle_degrees > points[points.len() - 1].angle_degrees
    {
        return (0.0, 0.0);
    }
    let right = points.partition_point(|point| point.angle_degrees <= angle_degrees);
    let left = right.saturating_sub(1).min(points.len() - 2);
    let right = left + 1;
    let fraction = (angle_degrees - points[left].angle_degrees)
        / (points[right].angle_degrees - points[left].angle_degrees);
    (
        (1.0 - fraction) * points[left].magnitude + fraction * points[right].magnitude,
        (1.0 - fraction) * points[left].phase_radians + fraction * points[right].phase_radians,
    )
}

fn reflected_state_is_finite(state: RayState) -> bool {
    state.position_m.iter().all(|value| value.is_finite())
        && state.tangent.iter().all(|value| value.is_finite())
        && state.p.iter().all(|value| value.is_finite())
        && state.q.iter().all(|value| value.is_finite())
        && state.speed_mps.is_finite()
        && state.travel_time_s.re.is_finite()
        && state.travel_time_s.im.is_finite()
        && state.amplitude.is_finite()
        && state.phase_radians.is_finite()
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use crate::model::{
        InternalReflectionCoefficientPoint, InternalReflectionCoefficientTable,
        ReflectionCoefficientPoint, ReflectionCoefficientTable,
    };

    use super::{interpolate_internal_reflection, interpolate_reflection};

    #[test]
    fn internal_reflection_uses_reference_quadratic_interpolator_and_power_scaling() {
        let table = InternalReflectionCoefficientTable {
            title: "test".to_owned(),
            frequency_hz: 100.0,
            points: vec![
                InternalReflectionCoefficientPoint {
                    horizontal_wavenumber_squared: 0.0,
                    f: Complex64::new(1.0, 0.0),
                    g: Complex64::new(2.0, 0.0),
                    decimal_power: 0,
                },
                InternalReflectionCoefficientPoint {
                    horizontal_wavenumber_squared: 1.0,
                    f: Complex64::new(2.0, 0.0),
                    g: Complex64::new(3.0, 0.0),
                    decimal_power: 0,
                },
                InternalReflectionCoefficientPoint {
                    horizontal_wavenumber_squared: 2.0,
                    f: Complex64::new(5.0e-1, 0.0),
                    g: Complex64::new(6.0e-1, 0.0),
                    decimal_power: 1,
                },
            ],
        };
        let (f, g) = interpolate_internal_reflection(&table, Complex64::new(0.5, 0.0)).unwrap();
        assert!((f.re - 1.25).abs() < 1.0e-15);
        assert!((g.re - 2.25).abs() < 1.0e-15);
    }

    #[test]
    fn reflection_tables_interpolate_and_zero_outside_domain() {
        let table = ReflectionCoefficientTable {
            points: vec![
                ReflectionCoefficientPoint {
                    angle_degrees: 0.0,
                    magnitude: 1.0,
                    phase_radians: 0.0,
                },
                ReflectionCoefficientPoint {
                    angle_degrees: 90.0,
                    magnitude: 0.5,
                    phase_radians: 1.0,
                },
            ],
        };
        let midpoint = interpolate_reflection(&table, 45.0);
        assert!((midpoint.0 - 0.75).abs() < 1.0e-15);
        assert!((midpoint.1 - 0.5).abs() < 1.0e-15);
        assert_eq!(interpolate_reflection(&table, -1.0), (0.0, 0.0));
    }
}
