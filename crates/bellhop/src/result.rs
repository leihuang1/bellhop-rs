use num_complex::Complex32;

/// Results from one complete simulation.
#[derive(Clone, Debug, PartialEq)]
pub struct SimulationResult {
    pub title: String,
    pub frequency_hz: f64,
    pub legacy_run_options: String,
    pub sources: Vec<SourceRaySet>,
    pub arrival_sources: Vec<SourceArrivals>,
    pub eigenray_sources: Vec<SourceEigenrays>,
    pub field_sources: Vec<SourceField>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SourceRaySet {
    pub source_depth_m: f64,
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

#[derive(Clone, Debug, PartialEq)]
pub struct SourceArrivals {
    pub source_depth_m: f64,
    pub receivers: Vec<ReceiverArrivals>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ReceiverArrivals {
    pub range_m: f64,
    pub depth_m: f64,
    pub arrivals: Vec<Arrival>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Arrival {
    pub amplitude: f32,
    pub phase_radians: f32,
    pub travel_time_s: f32,
    pub attenuation_time_s: f32,
    pub source_angle_degrees: f32,
    pub receiver_angle_degrees: f32,
    pub top_bounces: u32,
    pub bottom_bounces: u32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SourceEigenrays {
    pub source_depth_m: f64,
    pub receivers: Vec<ReceiverEigenrays>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ReceiverEigenrays {
    pub range_m: f64,
    pub depth_m: f64,
    pub eigenrays: Vec<RayTrajectory>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SourceField {
    pub source_depth_m: f64,
    pub samples: Vec<FieldSample>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FieldSample {
    pub range_m: f64,
    pub depth_m: f64,
    pub pressure: Complex32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RayTermination {
    ExitedTraceBox,
    LostEnergy,
    EscapedBoundary,
    SourceOutsideBoundaries,
    StepLimit,
    ReceiverHit,
}
