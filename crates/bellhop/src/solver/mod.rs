mod boundary;
mod integrator;
mod reflection;
mod ssp;

use crate::diagnostic::{Diagnostic, DiagnosticReport, SourceLocation};
use crate::model::{BoundaryCondition, Case, RunKind, SourceBeamPattern};

use boundary::{BoundaryGeometry, BoundarySide};
use integrator::{RayState, StepLimits, step_2d};
use num_complex::Complex64;
use reflection::reflect_2d;
use ssp::{SegmentState, SoundSpeedModel};

/// Resource limits for deterministic single-threaded simulation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SimulationLimits {
    pub max_rays: usize,
    pub max_steps_per_ray: usize,
    pub max_total_ray_points: usize,
}

impl Default for SimulationLimits {
    fn default() -> Self {
        Self {
            max_rays: 1_000_000,
            max_steps_per_ray: 100_000,
            max_total_ray_points: 20_000_000,
        }
    }
}

/// Results from one complete simulation.
#[derive(Clone, Debug, PartialEq)]
pub struct SimulationResult {
    pub sources: Vec<SourceRaySet>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SourceRaySet {
    pub source_depth_m: f32,
    pub rays: Vec<RayTrajectory>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RayTrajectory {
    pub launch_angle_degrees: f64,
    pub points: Vec<RayPoint>,
    pub top_bounces: u32,
    pub bottom_bounces: u32,
    pub termination: RayTermination,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RayPoint {
    pub range_m: f64,
    pub depth_m: f64,
    pub travel_time_s: f64,
    pub attenuation_time_s: f64,
    pub amplitude: f64,
    pub phase_radians: f64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RayTermination {
    ExitedTraceBox,
    LostEnergy,
    EscapedBoundary,
    SourceOutsideBoundaries,
    StepLimit,
}

/// Runs a loaded case using the deterministic compatibility path.
///
/// The current numerical milestone supports `R` (ray trace) cases. Other run
/// kinds return a structured unsupported-mode diagnostic.
///
/// # Errors
///
/// Returns diagnostics for unsupported modes, exceeded resource limits, or
/// non-finite and otherwise invalid numerical states.
#[allow(clippy::too_many_lines)]
pub fn run(case: &Case, limits: SimulationLimits) -> Result<SimulationResult, DiagnosticReport> {
    if case.environment.run.kind != RunKind::Rays {
        return Err(report(Diagnostic::error(
            "BH0301",
            "the numerical solver currently supports only R (ray trace) runs",
            "run_options.kind",
            SourceLocation::file(&case.environment.source_path),
        )));
    }
    let ray_count = selected_angles(case).len() * case.environment.positions.source_depths_m.len();
    if ray_count > limits.max_rays {
        return Err(report(Diagnostic::error(
            "BH0303",
            format!(
                "simulation requests {ray_count} rays; limit is {}",
                limits.max_rays
            ),
            "simulation_limits.max_rays",
            SourceLocation::file(&case.environment.source_path),
        )));
    }
    if limits.max_steps_per_ray < 2 || limits.max_total_ray_points < ray_count {
        return Err(report(Diagnostic::error(
            "BH0303",
            "simulation ray-point limits are too small",
            "simulation_limits",
            SourceLocation::file(&case.environment.source_path),
        )));
    }
    if matches!(
        case.environment.top_boundary.condition,
        BoundaryCondition::WriteReflectionCoefficient
    ) || matches!(
        case.environment.bottom_boundary.condition,
        BoundaryCondition::WriteReflectionCoefficient
    ) {
        return Err(report(Diagnostic::error(
            "BH0301",
            "W boundary reflection-table generation is not supported",
            "boundary.condition",
            SourceLocation::file(&case.environment.source_path),
        )));
    }

    let sound_speed = SoundSpeedModel::new(case).map_err(report)?;
    let boundaries = BoundaryGeometry::new(case);
    let step_m = if case.environment.trace.step_m == 0.0 {
        (case.environment.sound_speed.bottom_depth_m - case.environment.sound_speed.top_depth_m)
            / 10.0
    } else {
        case.environment.trace.step_m
    };
    let step_limits = StepLimits {
        base_step_m: step_m,
        max_range_m: case.environment.trace.max_range_m,
        max_depth_m: case.environment.trace.max_depth_m,
    };

    let mut total_points = 0_usize;
    let mut sources = Vec::with_capacity(case.environment.positions.source_depths_m.len());
    for &source_depth_m in &case.environment.positions.source_depths_m {
        let mut rays = Vec::with_capacity(selected_angles(case).len());
        for &launch_angle_degrees in selected_angles(case) {
            let ray = trace_ray(
                case,
                &sound_speed,
                &boundaries,
                step_limits,
                f64::from(source_depth_m),
                launch_angle_degrees,
                limits,
            )
            .map_err(|message| {
                report(Diagnostic::error(
                    "BH0302",
                    format!(
                        "ray at source depth {source_depth_m} m and launch angle {launch_angle_degrees}° failed: {message}"
                    ),
                    "ray_trace",
                    SourceLocation::file(&case.environment.source_path),
                ))
            })?;
            total_points = total_points.checked_add(ray.points.len()).ok_or_else(|| {
                report(Diagnostic::error(
                    "BH0303",
                    "total ray-point count overflowed",
                    "simulation_limits.max_total_ray_points",
                    SourceLocation::file(&case.environment.source_path),
                ))
            })?;
            if total_points > limits.max_total_ray_points {
                return Err(report(Diagnostic::error(
                    "BH0303",
                    format!(
                        "simulation produced more than {} ray points",
                        limits.max_total_ray_points
                    ),
                    "simulation_limits.max_total_ray_points",
                    SourceLocation::file(&case.environment.source_path),
                )));
            }
            rays.push(ray);
        }
        sources.push(SourceRaySet {
            source_depth_m,
            rays,
        });
    }
    Ok(SimulationResult { sources })
}

fn selected_angles(case: &Case) -> &[f64] {
    match case.environment.trace.selected_launch_angle {
        Some(index) => &case.environment.trace.launch_angles_degrees[index - 1..index],
        None => &case.environment.trace.launch_angles_degrees,
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn trace_ray(
    case: &Case,
    sound_speed: &SoundSpeedModel,
    boundaries: &BoundaryGeometry,
    step_limits: StepLimits,
    source_depth_m: f64,
    launch_angle_degrees: f64,
    limits: SimulationLimits,
) -> Result<RayTrajectory, &'static str> {
    let mut sound_segments = SegmentState::default();
    let source_position_m = [0.0, source_depth_m];
    let source_sample = sound_speed.evaluate(source_position_m, &mut sound_segments)?;
    let launch_angle_radians = launch_angle_degrees.to_radians();
    let amplitude = source_amplitude(case.source_beam_pattern.as_ref(), launch_angle_degrees);
    let geometric_q = case
        .environment
        .run
        .legacy
        .chars()
        .nth(1)
        .is_some_and(|option| option == 'G');
    let initial = RayState {
        position_m: source_position_m,
        tangent: [
            launch_angle_radians.cos() / source_sample.speed_mps,
            launch_angle_radians.sin() / source_sample.speed_mps,
        ],
        p: [1.0, 0.0],
        q: if geometric_q { [0.0, 0.0] } else { [0.0, 1.0] },
        speed_mps: source_sample.speed_mps,
        travel_time_s: Complex64::new(0.0, 0.0),
        amplitude,
        phase_radians: 0.0,
        top_bounces: 0,
        bottom_bounces: 0,
    };
    let mut states = Vec::with_capacity(limits.max_steps_per_ray.min(4096));
    states.push(initial);
    let mut top_segment = boundaries.top.segment_for_range(0.0);
    let mut bottom_segment = boundaries.bottom.segment_for_range(0.0);
    let mut beginning_top_distance = boundaries
        .top
        .signed_inside_distance(source_position_m, top_segment);
    let mut beginning_bottom_distance = boundaries
        .bottom
        .signed_inside_distance(source_position_m, bottom_segment);
    if beginning_top_distance <= 0.0 || beginning_bottom_distance <= 0.0 {
        return Ok(finish_trajectory(
            launch_angle_degrees,
            states,
            RayTermination::SourceOutsideBoundaries,
        ));
    }

    let mut consecutive_small_steps = 0_usize;
    let mut termination = RayTermination::StepLimit;
    while states.len() < limits.max_steps_per_ray {
        let incident = *states.last().expect("initial state is present");
        let stepped = step_2d(
            incident,
            sound_speed,
            &mut sound_segments,
            &boundaries.top,
            top_segment,
            &boundaries.bottom,
            bottom_segment,
            step_limits,
            &mut consecutive_small_steps,
        )?;
        states.push(stepped);

        let top_interval = boundaries.top.range_interval(top_segment);
        if stepped.position_m[0] < top_interval[0] || stepped.position_m[0] > top_interval[1] {
            top_segment = boundaries.top.segment_for_range(stepped.position_m[0]);
        }
        let bottom_interval = boundaries.bottom.range_interval(bottom_segment);
        if stepped.position_m[0] < bottom_interval[0] || stepped.position_m[0] > bottom_interval[1]
        {
            bottom_segment = boundaries.bottom.segment_for_range(stepped.position_m[0]);
        }

        let mut ending_top_distance = boundaries
            .top
            .signed_inside_distance(stepped.position_m, top_segment);
        let mut ending_bottom_distance = boundaries
            .bottom
            .signed_inside_distance(stepped.position_m, bottom_segment);

        if beginning_top_distance > 0.0 && ending_top_distance <= 0.0 {
            if states.len() >= limits.max_steps_per_ray {
                break;
            }
            let (tangent, normal) = boundaries
                .top
                .reflection_frame(stepped.position_m, top_segment);
            let reflected = reflect_2d(
                stepped,
                BoundarySide::Top,
                boundaries.top.segment(top_segment),
                tangent,
                normal,
                &case.environment.top_boundary.condition,
                case.top_reflection.as_ref(),
                &case.environment,
                sound_speed,
                &mut sound_segments,
                case.environment.run.beam_shift,
            )?;
            states.push(reflected);
            ending_top_distance = boundaries
                .top
                .signed_inside_distance(reflected.position_m, top_segment);
            ending_bottom_distance = boundaries
                .bottom
                .signed_inside_distance(reflected.position_m, bottom_segment);
        } else if beginning_bottom_distance > 0.0 && ending_bottom_distance <= 0.0 {
            if states.len() >= limits.max_steps_per_ray {
                break;
            }
            let (tangent, normal) = boundaries
                .bottom
                .reflection_frame(stepped.position_m, bottom_segment);
            let reflected = reflect_2d(
                stepped,
                BoundarySide::Bottom,
                boundaries.bottom.segment(bottom_segment),
                tangent,
                normal,
                &case.environment.bottom_boundary.condition,
                case.bottom_reflection.as_ref(),
                &case.environment,
                sound_speed,
                &mut sound_segments,
                case.environment.run.beam_shift,
            )?;
            states.push(reflected);
            ending_top_distance = boundaries
                .top
                .signed_inside_distance(reflected.position_m, top_segment);
            ending_bottom_distance = boundaries
                .bottom
                .signed_inside_distance(reflected.position_m, bottom_segment);
        }

        let current = *states.last().expect("step state is present");
        if current.position_m[0].abs() > step_limits.max_range_m
            || current.position_m[1].abs() > step_limits.max_depth_m
        {
            termination = RayTermination::ExitedTraceBox;
            break;
        }
        if current.amplitude < 0.005 {
            termination = RayTermination::LostEnergy;
            break;
        }
        if (beginning_top_distance < 0.0 && ending_top_distance < 0.0)
            || (beginning_bottom_distance < 0.0 && ending_bottom_distance < 0.0)
        {
            termination = RayTermination::EscapedBoundary;
            break;
        }
        beginning_top_distance = ending_top_distance;
        beginning_bottom_distance = ending_bottom_distance;
    }

    Ok(finish_trajectory(launch_angle_degrees, states, termination))
}

fn source_amplitude(pattern: Option<&SourceBeamPattern>, angle_degrees: f64) -> f64 {
    let Some(pattern) = pattern else {
        return 1.0;
    };
    let points = &pattern.points;
    let insertion = points.partition_point(|point| point.angle_degrees < angle_degrees);
    let left = insertion.saturating_sub(1).min(points.len() - 2);
    let fraction = (angle_degrees - points[left].angle_degrees)
        / (points[left + 1].angle_degrees - points[left].angle_degrees);
    (1.0 - fraction) * points[left].amplitude + fraction * points[left + 1].amplitude
}

fn finish_trajectory(
    launch_angle_degrees: f64,
    states: Vec<RayState>,
    termination: RayTermination,
) -> RayTrajectory {
    let final_state = *states.last().expect("ray has an initial state");
    RayTrajectory {
        launch_angle_degrees,
        points: states
            .into_iter()
            .map(|state| RayPoint {
                range_m: state.position_m[0],
                depth_m: state.position_m[1],
                travel_time_s: state.travel_time_s.re,
                attenuation_time_s: state.travel_time_s.im,
                amplitude: state.amplitude,
                phase_radians: state.phase_radians,
            })
            .collect(),
        top_bounces: final_state.top_bounces,
        bottom_bounces: final_state.bottom_bounces,
        termination,
    }
}

fn report(diagnostic: Diagnostic) -> DiagnosticReport {
    DiagnosticReport::from_diagnostic(diagnostic)
}
