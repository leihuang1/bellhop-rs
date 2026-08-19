use std::path::PathBuf;

use num_complex::Complex64;

/// A fully loaded and validated two-dimensional BELLHOP input case.
#[derive(Clone, Debug, PartialEq)]
pub struct Case {
    pub environment: EnvironmentCase,
    pub range_dependent_sound_speed: Option<RangeDependentSoundSpeed>,
    pub altimetry: Option<BoundaryShape>,
    pub bathymetry: Option<BoundaryShape>,
    pub bottom_reflection: Option<ReflectionCoefficientTable>,
    pub top_reflection: Option<ReflectionCoefficientTable>,
    pub internal_reflection: Option<InternalReflectionCoefficientTable>,
    pub source_beam_pattern: Option<SourceBeamPattern>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EnvironmentCase {
    pub source_path: PathBuf,
    pub title: String,
    pub frequency_hz: f64,
    pub top_options: TopOptions,
    pub top_boundary: Boundary,
    pub sound_speed: SoundSpeedInput,
    pub bottom_boundary: Boundary,
    pub positions: Positions,
    pub run: RunOptions,
    pub trace: TraceOptions,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SspInterpolation {
    N2Linear,
    CLinear,
    Pchip,
    CubicSpline,
    Quadrilateral,
    AnalyticMunk,
}

impl SspInterpolation {
    #[must_use]
    pub fn needs_range_dependent_file(self) -> bool {
        self == Self::Quadrilateral
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AttenuationUnit {
    NepersPerMeter,
    DbPerMeterKhz,
    DbPerMeter,
    DbPerWavelength,
    QualityFactor,
    LossParameter,
}

#[derive(Clone, Debug, PartialEq)]
pub enum VolumeAttenuation {
    None,
    Thorp,
    FrancoisGarrison {
        temperature_c: f64,
        salinity_psu: f64,
        ph: f64,
        mean_depth_m: f64,
    },
    Biological {
        layers: Vec<BiologicalLayer>,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub struct BiologicalLayer {
    pub top_depth_m: f64,
    pub bottom_depth_m: f64,
    pub resonance_frequency_hz: f64,
    pub quality_factor: f64,
    pub attenuation: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct TopOptions {
    pub legacy: String,
    pub interpolation: SspInterpolation,
    pub attenuation_unit: AttenuationUnit,
    pub volume_attenuation: VolumeAttenuation,
    pub has_altimetry: bool,
    pub development_options: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SoundSpeedInput {
    pub nominal_point_count: i32,
    pub surface_roughness_m: f32,
    pub top_depth_m: f64,
    pub bottom_depth_m: f64,
    pub points: Vec<SoundSpeedPoint>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SoundSpeedPoint {
    pub depth_m: f64,
    pub compressional_speed_mps: f64,
    pub shear_speed_mps: f64,
    pub density_g_cm3: f64,
    pub compressional_attenuation: f64,
    pub shear_attenuation: f64,
}

/// Sound speeds indexed as `[depth_index][range_index]`.
#[derive(Clone, Debug, PartialEq)]
pub struct RangeDependentSoundSpeed {
    pub ranges_m: Vec<f64>,
    pub depths_m: Vec<f64>,
    pub speeds_mps: Vec<Vec<f64>>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BoundaryInterpolation {
    PiecewiseLinear,
    Curvilinear,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BoundaryShape {
    pub interpolation: BoundaryInterpolation,
    pub points: Vec<BoundaryShapePoint>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BoundaryShapePoint {
    pub range_m: f64,
    pub depth_m: f64,
    pub material: Option<BoundaryMaterial>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BoundaryMaterial {
    pub compressional_speed_mps: f64,
    pub shear_speed_mps: f64,
    pub density_g_cm3: f64,
    pub compressional_attenuation: f64,
    pub shear_attenuation: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ReflectionCoefficientTable {
    pub points: Vec<ReflectionCoefficientPoint>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ReflectionCoefficientPoint {
    pub angle_degrees: f64,
    pub magnitude: f64,
    pub phase_radians: f64,
}

/// Impedance-function table stored by the Acoustics Toolbox `.irc` format.
#[derive(Clone, Debug, PartialEq)]
pub struct InternalReflectionCoefficientTable {
    pub title: String,
    pub frequency_hz: f64,
    pub points: Vec<InternalReflectionCoefficientPoint>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct InternalReflectionCoefficientPoint {
    pub horizontal_wavenumber_squared: f64,
    pub f: Complex64,
    pub g: Complex64,
    pub decimal_power: i32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SourceBeamPattern {
    pub points: Vec<SourceBeamPatternPoint>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SourceBeamPatternPoint {
    pub angle_degrees: f64,
    pub level_db: f64,
    pub amplitude: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Boundary {
    pub legacy_options: String,
    pub condition: BoundaryCondition,
    pub roughness_m: f32,
    pub has_shape_file: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BoundaryCondition {
    Vacuum,
    Rigid,
    AcoustoElastic(HalfSpace),
    GrainSize { depth_m: f64, phi: f64 },
    ReflectionCoefficientFile,
    WriteReflectionCoefficient,
    PrecalculatedReflectionCoefficient,
}

#[derive(Clone, Debug, PartialEq)]
pub struct HalfSpace {
    pub depth_m: f64,
    pub compressional_speed_mps: f64,
    pub shear_speed_mps: f64,
    pub density_g_cm3: f64,
    pub compressional_attenuation: f64,
    pub shear_attenuation: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Positions {
    /// BELLHOP stores source and receiver depths in single precision.
    pub source_depths_m: Vec<f32>,
    pub receiver_depths_m: Vec<f32>,
    pub receiver_ranges_m: Vec<f64>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RunKind {
    Rays,
    Eigenrays,
    Coherent,
    SemiCoherent,
    Incoherent,
    Arrivals,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LegacyArrivalEncoding {
    Ascii,
    Binary,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BeamFamily {
    GeometricHatCartesian,
    GeometricHatRayCentered,
    GeometricGaussianCartesian,
    GeometricGaussianRayCentered,
    SimpleGaussian,
    CervenyCartesian,
    CervenyRayCentered,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SourceGeometry {
    Point,
    Line,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReceiverGrid {
    Rectilinear,
    Irregular,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RunOptions {
    pub legacy: String,
    pub kind: RunKind,
    pub arrival_encoding: Option<LegacyArrivalEncoding>,
    pub beam_family: Option<BeamFamily>,
    pub has_source_beam_pattern: bool,
    pub source_geometry: SourceGeometry,
    pub receiver_grid: ReceiverGrid,
    pub beam_shift: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub struct TraceOptions {
    pub launch_angles_degrees: Vec<f64>,
    /// One-based to match the development option in the legacy file.
    pub selected_launch_angle: Option<usize>,
    pub step_m: f64,
    pub max_depth_m: f64,
    pub max_range_m: f64,
    pub cerveny: Option<CervenyOptions>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CervenyOptions {
    pub width: BeamWidth,
    pub curvature: CurvatureCondition,
    pub epsilon_multiplier: f64,
    pub loop_range: f64,
    pub image_count: i32,
    pub beam_window: i32,
    pub component: BeamComponent,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BeamWidth {
    SpaceFilling,
    Minimum,
    Wkb,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CurvatureCondition {
    Double,
    Standard,
    Zero,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BeamComponent {
    Pressure,
    Vertical,
    Horizontal,
    Displacement,
}
