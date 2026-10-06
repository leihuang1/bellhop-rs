use num_complex::Complex64;

#[derive(Clone, Debug, PartialEq)]
pub struct NormalMode {
    pub horizontal_wavenumber_rad_per_m: Complex64,
    pub phase_speed_mps: f64,
    pub group_speed_mps: f64,
    pub attenuation_nepers_per_m: f64,
    /// Complex pressure eigenfunction values, ordered like `ModeSet::sampled_depths_m`.
    pub eigenfunction: Vec<Complex64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ModeSet {
    pub frequency_hz: f64,
    pub sampled_depths_m: Vec<f64>,
    pub modes: Vec<NormalMode>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PressureField {
    pub source_depths_m: Vec<f64>,
    pub receiver_depths_m: Vec<f64>,
    pub receiver_ranges_m: Vec<f64>,
    pub receiver_offsets_m: Vec<f64>,
    /// Row-major values indexed as `[source_depth][receiver_depth][receiver_range]`.
    /// FIELD uses reference single-precision modal products and accumulation;
    /// final values are promoted to `Complex64` for the result model.
    pub pressure: Vec<Complex64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SimulationResult {
    pub modes: ModeSet,
    pub field: PressureField,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ProfileSimulationResult {
    /// One complete modal product per input profile, including modes above the FIELD cap.
    pub modes: Vec<ModeSet>,
    pub field: PressureField,
}
