use std::f64::consts::PI;

use num_complex::{Complex32, Complex64};

use crate::model::{Case, ReceiverGrid, SourceGeometry};

use super::integrator::RayState;
use super::{Arrival, RayTermination, ReceiverArrivals, ReceiverEigenrays, trajectory_from_states};

const PHASE_TOLERANCE: f64 = 0.05;

#[allow(clippy::struct_field_names)]
pub(super) struct InfluenceLimits {
    pub max_arrivals_per_receiver: usize,
    pub max_total_arrivals: usize,
    pub max_eigenrays: usize,
    pub max_eigenray_points: usize,
}

pub(super) struct InfluenceCounts {
    pub arrivals: usize,
    pub eigenrays: usize,
    pub eigenray_points: usize,
}

pub(super) enum InfluenceTarget<'a> {
    Eigenrays(&'a mut [ReceiverEigenrays]),
    Arrivals(&'a mut [ReceiverArrivals]),
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn geo_hat_cartesian(
    case: &Case,
    states: &[RayState],
    launch_angle_radians: f64,
    angular_spacing_radians: f64,
    target: &mut InfluenceTarget<'_>,
    limits: &InfluenceLimits,
    counts: &mut InfluenceCounts,
) -> Result<(), &'static str> {
    if states.len() < 2 {
        return Ok(());
    }
    let positions = &case.environment.positions;
    let ranges = &positions.receiver_ranges_m;
    let depths = &positions.receiver_depths_m;
    let rectilinear = case.environment.run.receiver_grid == ReceiverGrid::Rectilinear;
    let depth_count = if rectilinear { depths.len() } else { 1 };
    let q0 = states[0].speed_mps / angular_spacing_radians;
    let launch_angle_degrees = launch_angle_radians.to_degrees();
    let ratio = match case.environment.run.source_geometry {
        SourceGeometry::Point => launch_angle_radians.cos().abs().sqrt(),
        SourceGeometry::Line => 1.0,
    };
    let mut caustic_phase = 0.0;
    let mut q_old = states[0].q[0];
    let mut range_a = states[0].position_m[0];
    let mut receiver_index = ranges
        .iter()
        .position(|range| *range > range_a)
        .unwrap_or(0);
    if states[0].tangent[0] < 0.0 && receiver_index > 0 {
        receiver_index -= 1;
    }

    for state_index in 1..states.len() {
        let previous = states[state_index - 1];
        let current = states[state_index];
        let range_b = current.position_m[0];
        let segment = [
            current.position_m[0] - previous.position_m[0],
            current.position_m[1] - previous.position_m[1],
        ];
        let segment_length = segment[0].hypot(segment[1]);
        if segment_length < 1.0e3 * spacing(current.position_m[0]) {
            continue;
        }
        let tangent = [segment[0] / segment_length, segment[1] / segment_length];
        let normal = [-tangent[1], tangent[0]];
        let receiver_angle_degrees = tangent[1].atan2(tangent[0]).to_degrees();
        let delta_q = current.q[0] - previous.q[0];
        let delta_time = current.travel_time_s - previous.travel_time_s;
        let segment_q = previous.q[0];
        if crosses_caustic(segment_q, q_old) {
            caustic_phase += PI / 2.0;
        }
        q_old = segment_q;
        let radius_projected =
            previous.q[0].abs().max(current.q[0].abs()) / q0.abs() / tangent[0].abs();
        let (minimum_depth, maximum_depth) = if tangent[0].abs() > 0.5 {
            (
                previous.position_m[1].min(current.position_m[1]) - radius_projected,
                previous.position_m[1].max(current.position_m[1]) + radius_projected,
            )
        } else {
            (f64::NEG_INFINITY, f64::INFINITY)
        };

        loop {
            let receiver_range = ranges[receiver_index];
            if receiver_range >= range_a.min(range_b) && receiver_range < range_a.max(range_b) {
                for depth_index in 0..depth_count {
                    let receiver_depth = if rectilinear {
                        f64::from(depths[depth_index])
                    } else {
                        f64::from(depths[receiver_index])
                    };
                    if receiver_depth < minimum_depth || receiver_depth > maximum_depth {
                        continue;
                    }
                    let receiver_delta = [
                        receiver_range - previous.position_m[0],
                        receiver_depth - previous.position_m[1],
                    ];
                    let fraction = dot(receiver_delta, tangent) / segment_length;
                    let normal_distance = dot(receiver_delta, normal).abs();
                    let interpolated_q = previous.q[0] + fraction * delta_q;
                    let radius = (interpolated_q / q0).abs();
                    if normal_distance >= radius {
                        continue;
                    }
                    let delay = previous.travel_time_s + fraction * delta_time;
                    let constant = ratio
                        * (current.speed_mps / interpolated_q.abs()).sqrt()
                        * current.amplitude;
                    let hat_weight = (radius - normal_distance) / radius;
                    let amplitude = constant * hat_weight;
                    let mut phase = previous.phase_radians + caustic_phase;
                    if crosses_caustic(interpolated_q, q_old) {
                        phase += PI / 2.0;
                    }
                    let cell_index = receiver_index * depth_count + depth_index;
                    apply_contribution(
                        target,
                        cell_index,
                        state_index,
                        states,
                        ArrivalCandidate {
                            amplitude,
                            phase_radians: phase,
                            delay,
                            source_angle_degrees: launch_angle_degrees,
                            receiver_angle_degrees,
                            top_bounces: current.top_bounces,
                            bottom_bounces: current.bottom_bounces,
                        },
                        case.environment.frequency_hz,
                        limits,
                        counts,
                    )?;
                }
            }

            let next_index = if ranges[receiver_index] < range_b {
                if receiver_index + 1 >= ranges.len() {
                    break;
                }
                let next = receiver_index + 1;
                if ranges[next] >= range_b {
                    break;
                }
                next
            } else {
                if receiver_index == 0 {
                    break;
                }
                let next = receiver_index - 1;
                if ranges[next] <= range_b {
                    break;
                }
                next
            };
            receiver_index = next_index;
        }
        range_a = range_b;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn geo_hat_ray_centered(
    case: &Case,
    states: &[RayState],
    launch_angle_radians: f64,
    angular_spacing_radians: f64,
    target: &mut InfluenceTarget<'_>,
    limits: &InfluenceLimits,
    counts: &mut InfluenceCounts,
) -> Result<(), &'static str> {
    if states.len() < 2 {
        return Ok(());
    }
    let positions = &case.environment.positions;
    let ranges = &positions.receiver_ranges_m;
    let depths = &positions.receiver_depths_m;
    let rectilinear = case.environment.run.receiver_grid == ReceiverGrid::Rectilinear;
    let depth_count = if rectilinear { depths.len() } else { 1 };
    let q0 = states[0].speed_mps / angular_spacing_radians;
    let launch_angle_degrees = launch_angle_radians.to_degrees();
    let ratio = match case.environment.run.source_geometry {
        SourceGeometry::Point => launch_angle_radians.cos().abs().sqrt(),
        SourceGeometry::Line => 1.0,
    };
    let normal_depth: Vec<f64> = states
        .iter()
        .map(|state| -state.tangent[0] * state.speed_mps)
        .collect();
    let normal_range: Vec<f64> = states
        .iter()
        .map(|state| state.tangent[1] * state.speed_mps)
        .collect();

    for (depth_index, &receiver_depth) in depths.iter().take(depth_count).enumerate() {
        let receiver_depth = f64::from(receiver_depth);
        let mut caustic_phase = 0.0;
        let mut q_old = states[0].q[0];
        let (mut normal_a, mut range_a, mut receiver_a) = if normal_depth[0].abs() < 1.0e-6 {
            (1.0e10, 1.0e10, 0)
        } else {
            let normal = (receiver_depth - states[0].position_m[1]) / normal_depth[0];
            let range = states[0].position_m[0] + normal * normal_range[0];
            (normal, range, uniform_receiver_index(ranges, range))
        };

        for state_index in 1..states.len() {
            let previous = states[state_index - 1];
            let current = states[state_index];
            if normal_depth[state_index].abs() < 1.0e-10 {
                continue;
            }
            let normal_b = (receiver_depth - current.position_m[1]) / normal_depth[state_index];
            let range_b = current.position_m[0] + normal_b * normal_range[state_index];
            let receiver_b = uniform_receiver_index(ranges, range_b);
            if (current.position_m[0] - previous.position_m[0]).abs()
                < 1.0e3 * spacing(current.position_m[0])
                || receiver_a == receiver_b
            {
                normal_a = normal_b;
                range_a = range_b;
                receiver_a = receiver_b;
                continue;
            }

            let segment_q = previous.q[0];
            if crosses_caustic(segment_q, q_old) {
                caustic_phase += PI / 2.0;
            }
            q_old = segment_q;
            let receiver_angle_degrees =
                (current.tangent[1]).atan2(current.tangent[0]).to_degrees();
            let backwards = receiver_b <= receiver_a;
            let mut receiver_index = if backwards {
                receiver_a
            } else {
                receiver_a + 1
            };
            let final_index = if backwards {
                receiver_b + 1
            } else {
                receiver_b
            };
            loop {
                let fraction = (ranges[receiver_index] - range_a) / (range_b - range_a);
                let normal_distance = (normal_a + fraction * (normal_b - normal_a)).abs();
                let interpolated_q = previous.q[0] + fraction * (current.q[0] - previous.q[0]);
                let radius = interpolated_q.abs() / q0;
                if normal_distance < radius {
                    let delay = previous.travel_time_s
                        + fraction * (current.travel_time_s - previous.travel_time_s);
                    let constant = ratio * current.speed_mps.sqrt() * current.amplitude
                        / interpolated_q.abs().sqrt();
                    let hat_weight = (radius - normal_distance) / radius;
                    let amplitude = constant * hat_weight;
                    let mut phase = previous.phase_radians + caustic_phase;
                    if crosses_caustic(interpolated_q, q_old) {
                        phase += PI / 2.0;
                    }
                    let cell_index = receiver_index * depth_count + depth_index;
                    apply_contribution(
                        target,
                        cell_index,
                        state_index,
                        states,
                        ArrivalCandidate {
                            amplitude,
                            phase_radians: phase,
                            delay,
                            source_angle_degrees: launch_angle_degrees,
                            receiver_angle_degrees,
                            top_bounces: current.top_bounces,
                            bottom_bounces: current.bottom_bounces,
                        },
                        case.environment.frequency_hz,
                        limits,
                        counts,
                    )?;
                }
                if receiver_index == final_index {
                    break;
                }
                if backwards {
                    receiver_index -= 1;
                } else {
                    receiver_index += 1;
                }
            }
            normal_a = normal_b;
            range_a = range_b;
            receiver_a = receiver_b;
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn simple_gaussian_eigenrays(
    case: &Case,
    states: &[RayState],
    launch_angle_radians: f64,
    target: &mut InfluenceTarget<'_>,
    limits: &InfluenceLimits,
    counts: &mut InfluenceCounts,
) -> Result<(), &'static str> {
    if states.len() < 2 {
        return Ok(());
    }
    let ranges = &case.environment.positions.receiver_ranges_m;
    let depth_count = match case.environment.run.receiver_grid {
        ReceiverGrid::Rectilinear => case.environment.positions.receiver_depths_m.len(),
        ReceiverGrid::Irregular => 1,
    };
    let launch_angle_degrees = launch_angle_radians.to_degrees();
    let mut range_a = states[0].position_m[0];
    let mut receiver_index = 0_usize;
    for state_index in 1..states.len() {
        let current = states[state_index];
        let range_b = current.position_m[0];
        while (range_b - range_a).abs() > 1.0e3 * spacing(range_a)
            && range_b > ranges[receiver_index]
        {
            for depth_index in 0..depth_count {
                let cell_index = receiver_index * depth_count + depth_index;
                apply_contribution(
                    target,
                    cell_index,
                    state_index,
                    states,
                    ArrivalCandidate {
                        amplitude: 0.0,
                        phase_radians: 0.0,
                        delay: Complex64::new(0.0, 0.0),
                        source_angle_degrees: launch_angle_degrees,
                        receiver_angle_degrees: 0.0,
                        top_bounces: current.top_bounces,
                        bottom_bounces: current.bottom_bounces,
                    },
                    case.environment.frequency_hz,
                    limits,
                    counts,
                )?;
            }
            receiver_index += 1;
            if receiver_index >= ranges.len() {
                return Ok(());
            }
        }
        range_a = range_b;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn geo_gaussian_cartesian(
    case: &Case,
    states: &[RayState],
    launch_angle_radians: f64,
    angular_spacing_radians: f64,
    target: &mut InfluenceTarget<'_>,
    limits: &InfluenceLimits,
    counts: &mut InfluenceCounts,
) -> Result<(), &'static str> {
    if states.len() < 2 {
        return Ok(());
    }
    let positions = &case.environment.positions;
    let ranges = &positions.receiver_ranges_m;
    let depths = &positions.receiver_depths_m;
    let rectilinear = case.environment.run.receiver_grid == ReceiverGrid::Rectilinear;
    let depth_count = if rectilinear { depths.len() } else { 1 };
    let q0 = states[0].speed_mps / angular_spacing_radians;
    let launch_angle_degrees = launch_angle_radians.to_degrees();
    let ratio = match case.environment.run.source_geometry {
        SourceGeometry::Point => launch_angle_radians.cos().abs().sqrt() / (2.0 * PI).sqrt(),
        SourceGeometry::Line => 1.0 / (2.0 * PI).sqrt(),
    };
    let mut caustic_phase = 0.0;
    let mut q_old = states[0].q[0];
    let mut range_a = states[0].position_m[0];
    let mut receiver_index = ranges
        .iter()
        .position(|range| *range > range_a)
        .unwrap_or(0);
    if states[0].tangent[0] < 0.0 && receiver_index > 0 {
        receiver_index -= 1;
    }

    for state_index in 1..states.len() {
        let previous = states[state_index - 1];
        let current = states[state_index];
        let range_b = current.position_m[0];
        let segment = [
            current.position_m[0] - previous.position_m[0],
            current.position_m[1] - previous.position_m[1],
        ];
        let segment_length = segment[0].hypot(segment[1]);
        if segment_length < 1.0e3 * spacing(current.position_m[0]) {
            continue;
        }
        let tangent = [segment[0] / segment_length, segment[1] / segment_length];
        let normal = [-tangent[1], tangent[0]];
        let receiver_angle_degrees = tangent[1].atan2(tangent[0]).to_degrees();
        let delta_q = current.q[0] - previous.q[0];
        let delta_time = current.travel_time_s - previous.travel_time_s;
        let segment_q = previous.q[0];
        if crosses_caustic(segment_q, q_old) {
            caustic_phase += PI / 2.0;
        }
        q_old = segment_q;
        let wavelength = previous.speed_mps / case.environment.frequency_hz;
        let mut projected_sigma =
            previous.q[0].abs().max(current.q[0].abs()) / q0.abs() / tangent[0].abs();
        let diffraction_floor =
            (0.2 * case.environment.frequency_hz * current.travel_time_s.re).min(PI * wavelength);
        projected_sigma = projected_sigma.max(diffraction_floor);
        let radius = 4.0 * projected_sigma;
        let (minimum_depth, maximum_depth) = if tangent[0].abs() > 0.5 {
            (
                previous.position_m[1].min(current.position_m[1]) - radius,
                previous.position_m[1].max(current.position_m[1]) + radius,
            )
        } else {
            (f64::NEG_INFINITY, f64::INFINITY)
        };

        loop {
            let receiver_range = ranges[receiver_index];
            if receiver_range >= range_a.min(range_b) && receiver_range < range_a.max(range_b) {
                for depth_index in 0..depth_count {
                    let receiver_depth = if rectilinear {
                        f64::from(depths[depth_index])
                    } else {
                        f64::from(depths[receiver_index])
                    };
                    if receiver_depth < minimum_depth || receiver_depth > maximum_depth {
                        continue;
                    }
                    let receiver_delta = [
                        receiver_range - previous.position_m[0],
                        receiver_depth - previous.position_m[1],
                    ];
                    let fraction = dot(receiver_delta, tangent) / segment_length;
                    let normal_distance = dot(receiver_delta, normal).abs();
                    let interpolated_q = previous.q[0] + fraction * delta_q;
                    let mut sigma = (interpolated_q / q0).abs();
                    sigma = sigma.max(diffraction_floor);
                    if normal_distance >= 4.0 * sigma {
                        continue;
                    }
                    let spreading = (q0 / interpolated_q).abs();
                    let delay = previous.travel_time_s + fraction * delta_time;
                    let constant = ratio
                        * (current.speed_mps / interpolated_q.abs()).sqrt()
                        * current.amplitude;
                    let gaussian_weight =
                        (-0.5 * (normal_distance / sigma).powi(2)).exp() / (sigma * spreading);
                    let amplitude = constant * gaussian_weight;
                    let mut phase = previous.phase_radians + caustic_phase;
                    if crosses_caustic(interpolated_q, q_old) {
                        phase += PI / 2.0;
                    }
                    let cell_index = receiver_index * depth_count + depth_index;
                    apply_contribution(
                        target,
                        cell_index,
                        state_index,
                        states,
                        ArrivalCandidate {
                            amplitude,
                            phase_radians: phase,
                            delay,
                            source_angle_degrees: launch_angle_degrees,
                            receiver_angle_degrees,
                            top_bounces: current.top_bounces,
                            bottom_bounces: current.bottom_bounces,
                        },
                        case.environment.frequency_hz,
                        limits,
                        counts,
                    )?;
                }
            }

            let next_index = if range_b > ranges[receiver_index] {
                if receiver_index + 1 >= ranges.len() {
                    break;
                }
                let next = receiver_index + 1;
                if ranges[next] >= range_b {
                    break;
                }
                next
            } else {
                if receiver_index == 0 {
                    break;
                }
                let next = receiver_index - 1;
                if ranges[next] <= range_b {
                    break;
                }
                next
            };
            receiver_index = next_index;
        }
        range_a = range_b;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn apply_contribution(
    target: &mut InfluenceTarget<'_>,
    cell_index: usize,
    state_index: usize,
    states: &[RayState],
    candidate: ArrivalCandidate,
    frequency_hz: f64,
    limits: &InfluenceLimits,
    counts: &mut InfluenceCounts,
) -> Result<(), &'static str> {
    match target {
        InfluenceTarget::Eigenrays(receivers) => {
            let point_count = state_index + 1;
            if counts.eigenrays >= limits.max_eigenrays
                || counts.eigenray_points.saturating_add(point_count) > limits.max_eigenray_points
            {
                return Err("eigenray resource limit exceeded");
            }
            receivers[cell_index].eigenrays.push(trajectory_from_states(
                candidate.source_angle_degrees,
                &states[..=state_index],
                RayTermination::ReceiverHit,
            ));
            counts.eigenrays += 1;
            counts.eigenray_points += point_count;
            Ok(())
        }
        InfluenceTarget::Arrivals(receivers) => add_arrival(
            &mut receivers[cell_index],
            candidate,
            frequency_hz,
            limits,
            counts,
        ),
    }
}

#[allow(clippy::cast_possible_truncation)]
pub(super) fn scale_arrivals(case: &Case, receivers: &mut [ReceiverArrivals]) {
    for receiver in receivers {
        let factor = match case.environment.run.source_geometry {
            SourceGeometry::Line => 4.0 * PI.sqrt(),
            SourceGeometry::Point if receiver.range_m == 0.0 => 1.0e5,
            SourceGeometry::Point => 1.0 / receiver.range_m.sqrt(),
        } as f32;
        for arrival in &mut receiver.arrivals {
            arrival.amplitude *= factor;
        }
    }
}

struct ArrivalCandidate {
    amplitude: f64,
    phase_radians: f64,
    delay: Complex64,
    source_angle_degrees: f64,
    receiver_angle_degrees: f64,
    top_bounces: u32,
    bottom_bounces: u32,
}

fn add_arrival(
    receiver: &mut ReceiverArrivals,
    candidate: ArrivalCandidate,
    frequency_hz: f64,
    limits: &InfluenceLimits,
    counts: &mut InfluenceCounts,
) -> Result<(), &'static str> {
    let angular_frequency = 2.0 * PI * frequency_hz;
    let combine = receiver.arrivals.last().is_some_and(|last| {
        angular_frequency
            * (candidate.delay
                - Complex64::new(
                    f64::from(last.travel_time_s),
                    f64::from(last.attenuation_time_s),
                ))
            .norm()
            < PHASE_TOLERANCE
            && (f64::from(last.phase_radians) - candidate.phase_radians).abs() < PHASE_TOLERANCE
    });
    if combine {
        combine_arrival(
            receiver.arrivals.last_mut().expect("arrival was present"),
            &candidate,
        );
        return Ok(());
    }

    let arrival = candidate.into_arrival();
    if receiver.arrivals.len() >= limits.max_arrivals_per_receiver {
        let (weakest_index, weakest) = receiver
            .arrivals
            .iter()
            .enumerate()
            .min_by(|left, right| left.1.amplitude.total_cmp(&right.1.amplitude))
            .expect("arrival limit is positive");
        if arrival.amplitude > weakest.amplitude {
            receiver.arrivals[weakest_index] = arrival;
        }
        return Ok(());
    }
    if counts.arrivals >= limits.max_total_arrivals {
        return Err("total arrival resource limit exceeded");
    }
    receiver.arrivals.push(arrival);
    counts.arrivals += 1;
    Ok(())
}

#[allow(clippy::cast_possible_truncation)]
fn combine_arrival(arrival: &mut Arrival, candidate: &ArrivalCandidate) {
    let added_amplitude = candidate.amplitude as f32;
    let total_amplitude = arrival.amplitude + added_amplitude;
    let old_weight = arrival.amplitude / total_amplitude;
    let new_weight = added_amplitude / total_amplitude;
    let old_delay = Complex32::new(arrival.travel_time_s, arrival.attenuation_time_s);
    let new_delay = Complex32::new(candidate.delay.re as f32, candidate.delay.im as f32);
    let delay = old_weight * old_delay + new_weight * new_delay;
    arrival.amplitude = total_amplitude;
    arrival.travel_time_s = delay.re;
    arrival.attenuation_time_s = delay.im;
    arrival.source_angle_degrees = old_weight * arrival.source_angle_degrees
        + new_weight * candidate.source_angle_degrees as f32;
    arrival.receiver_angle_degrees = old_weight * arrival.receiver_angle_degrees
        + new_weight * candidate.receiver_angle_degrees as f32;
}

impl ArrivalCandidate {
    #[allow(clippy::cast_possible_truncation)]
    fn into_arrival(self) -> Arrival {
        Arrival {
            amplitude: self.amplitude as f32,
            phase_radians: self.phase_radians as f32,
            travel_time_s: self.delay.re as f32,
            attenuation_time_s: self.delay.im as f32,
            source_angle_degrees: self.source_angle_degrees as f32,
            receiver_angle_degrees: self.receiver_angle_degrees as f32,
            top_bounces: self.top_bounces,
            bottom_bounces: self.bottom_bounces,
        }
    }
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss
)]
fn uniform_receiver_index(ranges: &[f64], range: f64) -> usize {
    if ranges.len() == 1 {
        return 0;
    }
    let spacing = ranges[ranges.len() - 1] - ranges[ranges.len() - 2];
    let index = ((range - ranges[0]) / spacing).trunc();
    if !index.is_finite() || index <= 0.0 {
        0
    } else if index >= (ranges.len() - 1) as f64 {
        ranges.len() - 1
    } else {
        index as usize
    }
}

fn crosses_caustic(q: f64, old_q: f64) -> bool {
    (q <= 0.0 && old_q > 0.0) || (q >= 0.0 && old_q < 0.0)
}

fn dot(left: [f64; 2], right: [f64; 2]) -> f64 {
    left[0] * right[0] + left[1] * right[1]
}

fn spacing(value: f64) -> f64 {
    let magnitude = value.abs();
    if !magnitude.is_finite() {
        return f64::NAN;
    }
    if magnitude == 0.0 {
        return f64::MIN_POSITIVE * f64::EPSILON;
    }
    let bits = magnitude.to_bits();
    f64::from_bits(bits + 1) - magnitude
}

#[cfg(test)]
mod tests {
    use super::{crosses_caustic, spacing};

    #[test]
    #[allow(clippy::float_cmp)]
    fn spacing_matches_binary64_bins() {
        assert_eq!(spacing(1.0), f64::EPSILON);
        assert_eq!(spacing(-2.0), 2.0 * f64::EPSILON);
        assert_eq!(spacing(0.0), f64::from_bits(1));
    }

    #[test]
    fn caustic_crossing_uses_reference_inclusive_zero_test() {
        assert!(crosses_caustic(0.0, 1.0));
        assert!(crosses_caustic(0.0, -1.0));
        assert!(!crosses_caustic(1.0, 1.0));
    }
}
