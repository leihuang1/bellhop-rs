use std::fs;
use std::path::{Path, PathBuf};

use bellhop::legacy::load_case;
use bellhop::solver::{RayPoint, RayTrajectory, SimulationLimits, run};

#[derive(Debug)]
struct ReferenceRay {
    launch_angle_degrees: f64,
    top_bounces: u32,
    bottom_bounces: u32,
    points: Vec<[f64; 2]>,
}

#[test]
#[ignore = "requires a pinned external Fortran reference run"]
fn ray_trajectories_match_pinned_linux_reference() {
    let environment = required_path("BELLHOP_DIFFERENTIAL_ENV");
    let reference_path = required_path("BELLHOP_DIFFERENTIAL_RAY");
    let tolerance = std::env::var("BELLHOP_DIFFERENTIAL_POSITION_TOLERANCE_M")
        .map_or(Ok(1.0e-5), |value| value.parse::<f64>())
        .expect("position tolerance must be a number");

    let case = load_case(&environment).unwrap().value;
    let base_step_m = if case.environment.trace.step_m == 0.0 {
        (case.environment.sound_speed.bottom_depth_m - case.environment.sound_speed.top_depth_m)
            / 10.0
    } else {
        case.environment.trace.step_m
    };
    let minimum_step_factor = std::env::var("BELLHOP_DIFFERENTIAL_MINIMUM_STEP_FACTOR")
        .map_or(Ok(4.1), |value| value.parse::<f64>())
        .expect("minimum-step tolerance factor must be a number");
    assert!(minimum_step_factor >= 1.0);
    let minimum_step_tolerance = (minimum_step_factor * 1.0e-4 * base_step_m).max(tolerance);
    let result = run(&case, SimulationLimits::default()).unwrap();
    let actual: Vec<&RayTrajectory> = result
        .sources
        .iter()
        .flat_map(|source| &source.rays)
        .collect();
    let expected = parse_reference_rays(&reference_path);
    assert_eq!(actual.len(), expected.len(), "ray count differs");

    let mut maximum_error = 0.0_f64;
    let mut maximum_minimum_step_error = 0.0_f64;
    let mut skipped_minimum_steps = 0_usize;
    for (ray_index, (actual, expected)) in actual.iter().zip(&expected).enumerate() {
        assert!(
            (actual.launch_angle_degrees - expected.launch_angle_degrees).abs() < 1.0e-10,
            "ray {ray_index} launch angle differs"
        );
        assert_eq!(
            actual.top_bounces, expected.top_bounces,
            "ray {ray_index} top-bounce count differs"
        );
        assert_eq!(
            actual.bottom_bounces, expected.bottom_bounces,
            "ray {ray_index} bottom-bounce count differs"
        );

        if actual.points.len() == expected.points.len()
            && actual
                .points
                .iter()
                .zip(&expected.points)
                .all(|(actual, &expected)| {
                    coordinate_error(actual.range_m, actual.depth_m, expected) <= tolerance
                })
        {
            for (actual, &expected) in actual.points.iter().zip(&expected.points) {
                maximum_error =
                    maximum_error.max(coordinate_error(actual.range_m, actual.depth_m, expected));
            }
            continue;
        }

        let (branch_error, skipped) = align_compiler_sensitive_steps(
            &actual.points,
            &expected.points,
            minimum_step_tolerance,
        )
        .unwrap_or_else(|| {
            panic!(
                "ray {ray_index} cannot be aligned within the {minimum_step_tolerance:e} m minimum-step tolerance (Rust points={}, Fortran points={})",
                actual.points.len(),
                expected.points.len()
            )
        });
        maximum_minimum_step_error = maximum_minimum_step_error.max(branch_error);
        skipped_minimum_steps += skipped;
    }
    eprintln!(
        "compared {} rays; maximum aligned coordinate error {maximum_error:e} m; maximum compiler-sensitive branch error {maximum_minimum_step_error:e} m; accepted {skipped_minimum_steps} minimum-step differences",
        actual.len()
    );
}

fn align_compiler_sensitive_steps(
    actual: &[RayPoint],
    expected: &[[f64; 2]],
    tolerance: f64,
) -> Option<(f64, usize)> {
    let columns = expected.len() + 1;
    let mut parent = vec![0_u8; (actual.len() + 1).checked_mul(columns)?];
    parent[0] = 4;
    for actual_index in 0..=actual.len() {
        for expected_index in 0..=expected.len() {
            if parent[actual_index * columns + expected_index] == 0 {
                continue;
            }
            if actual_index < actual.len()
                && expected_index < expected.len()
                && coordinate_error(
                    actual[actual_index].range_m,
                    actual[actual_index].depth_m,
                    expected[expected_index],
                ) <= tolerance
            {
                let next = (actual_index + 1) * columns + expected_index + 1;
                if parent[next] == 0 {
                    parent[next] = 1;
                }
            }
            if is_short_actual_step(actual, actual_index, tolerance) {
                let next = (actual_index + 1) * columns + expected_index;
                if parent[next] == 0 {
                    parent[next] = 2;
                }
            }
            if is_short_expected_step(expected, expected_index, tolerance) {
                let next = actual_index * columns + expected_index + 1;
                if parent[next] == 0 {
                    parent[next] = 3;
                }
            }
        }
    }

    let mut actual_index = actual.len();
    let mut expected_index = expected.len();
    if parent[actual_index * columns + expected_index] == 0 {
        return None;
    }
    let mut maximum_error = 0.0_f64;
    let mut skipped = 0;
    while actual_index > 0 || expected_index > 0 {
        match parent[actual_index * columns + expected_index] {
            1 => {
                actual_index -= 1;
                expected_index -= 1;
                maximum_error = maximum_error.max(coordinate_error(
                    actual[actual_index].range_m,
                    actual[actual_index].depth_m,
                    expected[expected_index],
                ));
            }
            2 => {
                actual_index -= 1;
                skipped += 1;
            }
            3 => {
                expected_index -= 1;
                skipped += 1;
            }
            _ => return None,
        }
    }
    Some((maximum_error, skipped))
}

fn is_short_actual_step(points: &[RayPoint], index: usize, tolerance: f64) -> bool {
    index > 0
        && index + 1 < points.len()
        && coordinate_error(
            points[index].range_m,
            points[index].depth_m,
            [points[index - 1].range_m, points[index - 1].depth_m],
        ) <= tolerance
}

fn is_short_expected_step(points: &[[f64; 2]], index: usize, tolerance: f64) -> bool {
    index > 0
        && index + 1 < points.len()
        && coordinate_error(points[index][0], points[index][1], points[index - 1]) <= tolerance
}

fn coordinate_error(range_m: f64, depth_m: f64, expected: [f64; 2]) -> f64 {
    (range_m - expected[0])
        .abs()
        .max((depth_m - expected[1]).abs())
}

fn required_path(variable: &str) -> PathBuf {
    std::env::var_os(variable).map_or_else(
        || panic!("{variable} must name an input file"),
        PathBuf::from,
    )
}

fn parse_reference_rays(path: &Path) -> Vec<ReferenceRay> {
    let source = fs::read_to_string(path).expect("reference .ray file must be readable");
    let mut records = source
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty());
    for description in [
        "title",
        "frequency",
        "source counts",
        "launch-angle counts",
        "top depth",
        "bottom depth",
        "coordinate type",
    ] {
        records
            .next()
            .unwrap_or_else(|| panic!("reference .ray file is missing {description}"));
    }

    let mut rays = Vec::new();
    while let Some(angle_record) = records.next() {
        let launch_angle_degrees = parse_values::<f64>(angle_record, "launch angle")[0];
        let counts_record = records.next().expect("ray point-count record is missing");
        let counts = parse_values::<u32>(counts_record, "ray point counts");
        assert_eq!(counts.len(), 3, "ray count record must have three values");
        let point_count = usize::try_from(counts[0]).expect("point count must fit usize");
        let mut points = Vec::with_capacity(point_count);
        for _ in 0..point_count {
            let point_record = records.next().expect("ray coordinate record is missing");
            let point = parse_values::<f64>(point_record, "ray coordinates");
            assert_eq!(point.len(), 2, "2D ray point must have two coordinates");
            points.push([point[0], point[1]]);
        }
        rays.push(ReferenceRay {
            launch_angle_degrees,
            top_bounces: counts[1],
            bottom_bounces: counts[2],
            points,
        });
    }
    rays
}

fn parse_values<T: std::str::FromStr>(record: &str, description: &str) -> Vec<T> {
    record
        .split_whitespace()
        .map(|value| {
            value
                .parse::<T>()
                .unwrap_or_else(|_| panic!("invalid {description} value {value:?}"))
        })
        .collect()
}
