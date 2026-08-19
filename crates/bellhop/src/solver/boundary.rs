use crate::model::{BoundaryInterpolation, BoundaryMaterial, BoundaryShape, Case};

const EXTENDED_RANGE_M: f64 = 1.340_780_792_994_259_6e149;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum BoundarySide {
    Top,
    Bottom,
}

#[derive(Clone, Debug)]
pub(crate) struct BoundaryGeometry {
    pub top: BoundaryCurve,
    pub bottom: BoundaryCurve,
}

#[derive(Clone, Debug)]
pub(crate) struct BoundaryCurve {
    interpolation: BoundaryInterpolation,
    nodes: Vec<BoundaryNode>,
    segments: Vec<BoundarySegment>,
}

#[derive(Clone, Debug)]
struct BoundaryNode {
    position_m: [f64; 2],
    tangent: [f64; 2],
    normal: [f64; 2],
    material: Option<BoundaryMaterial>,
}

#[derive(Clone, Debug)]
pub(crate) struct BoundarySegment {
    pub origin_m: [f64; 2],
    pub tangent: [f64; 2],
    pub normal: [f64; 2],
    pub length_m: f64,
    pub curvature: f64,
    pub material: Option<BoundaryMaterial>,
}

impl BoundaryGeometry {
    pub fn new(case: &Case) -> Self {
        let environment = &case.environment;
        Self {
            top: BoundaryCurve::new(
                BoundarySide::Top,
                environment.sound_speed.top_depth_m,
                case.altimetry.as_ref(),
            ),
            bottom: BoundaryCurve::new(
                BoundarySide::Bottom,
                environment.sound_speed.bottom_depth_m,
                case.bathymetry.as_ref(),
            ),
        }
    }
}

impl BoundaryCurve {
    fn new(side: BoundarySide, flat_depth_m: f64, shape: Option<&BoundaryShape>) -> Self {
        let interpolation = shape.map_or(BoundaryInterpolation::PiecewiseLinear, |shape| {
            shape.interpolation
        });
        let mut nodes = match shape {
            Some(shape) => shape
                .points
                .iter()
                .map(|point| BoundaryNode {
                    position_m: [point.range_m, point.depth_m],
                    tangent: [0.0; 2],
                    normal: [0.0; 2],
                    material: point.material.clone(),
                })
                .collect::<Vec<_>>(),
            None => vec![
                BoundaryNode {
                    position_m: [-EXTENDED_RANGE_M, flat_depth_m],
                    tangent: [0.0; 2],
                    normal: [0.0; 2],
                    material: None,
                },
                BoundaryNode {
                    position_m: [EXTENDED_RANGE_M, flat_depth_m],
                    tangent: [0.0; 2],
                    normal: [0.0; 2],
                    material: None,
                },
            ],
        };
        if shape.is_some() {
            let first = nodes[0].clone();
            let last = nodes[nodes.len() - 1].clone();
            nodes.insert(
                0,
                BoundaryNode {
                    position_m: [-EXTENDED_RANGE_M, first.position_m[1]],
                    tangent: [0.0; 2],
                    normal: [0.0; 2],
                    material: first.material,
                },
            );
            nodes.push(BoundaryNode {
                position_m: [EXTENDED_RANGE_M, last.position_m[1]],
                tangent: [0.0; 2],
                normal: [0.0; 2],
                material: last.material,
            });
        }

        let mut segments = Vec::with_capacity(nodes.len() - 1);
        let mut slopes = Vec::with_capacity(nodes.len());
        for index in 0..nodes.len() - 1 {
            let delta = subtract(nodes[index + 1].position_m, nodes[index].position_m);
            let length_m = norm(delta);
            let tangent = scale(delta, 1.0 / length_m);
            let normal = outward_normal(side, tangent);
            slopes.push(delta[1] / delta[0]);
            segments.push(BoundarySegment {
                origin_m: nodes[index].position_m,
                tangent,
                normal,
                length_m,
                curvature: 0.0,
                material: nodes[index].material.clone(),
            });
        }
        slopes.push(0.0);

        if interpolation == BoundaryInterpolation::Curvilinear {
            nodes[0].tangent = [1.0, 0.0];
            let final_node = nodes.len() - 1;
            nodes[final_node].tangent = [1.0, 0.0];
            for index in 1..final_node {
                nodes[index].tangent = scale(
                    add(segments[index - 1].tangent, segments[index].tangent),
                    0.5,
                );
            }
            for node in &mut nodes {
                node.normal = outward_normal(side, node.tangent);
            }
            for index in 0..segments.len() {
                let range_delta = nodes[index + 1].position_m[0] - nodes[index].position_m[0];
                let second_derivative = (slopes[index + 1] - slopes[index]) / range_delta;
                segments[index].curvature = second_derivative * segments[index].tangent[0].powi(3);
            }
        }

        Self {
            interpolation,
            nodes,
            segments,
        }
    }

    pub fn segment_for_range(&self, range_m: f64) -> usize {
        let insertion = self
            .nodes
            .partition_point(|node| node.position_m[0] < range_m);
        insertion.saturating_sub(1).min(self.segments.len() - 1)
    }

    pub fn segment(&self, index: usize) -> &BoundarySegment {
        &self.segments[index]
    }

    pub fn range_interval(&self, index: usize) -> [f64; 2] {
        [
            self.nodes[index].position_m[0],
            self.nodes[index + 1].position_m[0],
        ]
    }

    pub fn signed_inside_distance(&self, position_m: [f64; 2], segment: usize) -> f64 {
        -dot(
            self.segments[segment].normal,
            subtract(position_m, self.segments[segment].origin_m),
        )
    }

    pub fn reflection_frame(&self, position_m: [f64; 2], segment: usize) -> ([f64; 2], [f64; 2]) {
        if self.interpolation == BoundaryInterpolation::Curvilinear {
            let boundary_segment = &self.segments[segment];
            let proportion = dot(
                subtract(position_m, boundary_segment.origin_m),
                boundary_segment.tangent,
            ) / boundary_segment.length_m;
            (
                add(
                    scale(self.nodes[segment].tangent, 1.0 - proportion),
                    scale(self.nodes[segment + 1].tangent, proportion),
                ),
                add(
                    scale(self.nodes[segment].normal, 1.0 - proportion),
                    scale(self.nodes[segment + 1].normal, proportion),
                ),
            )
        } else {
            (
                self.segments[segment].tangent,
                self.segments[segment].normal,
            )
        }
    }
}

pub(crate) fn dot(left: [f64; 2], right: [f64; 2]) -> f64 {
    left[0] * right[0] + left[1] * right[1]
}

pub(crate) fn add(left: [f64; 2], right: [f64; 2]) -> [f64; 2] {
    [left[0] + right[0], left[1] + right[1]]
}

pub(crate) fn subtract(left: [f64; 2], right: [f64; 2]) -> [f64; 2] {
    [left[0] - right[0], left[1] - right[1]]
}

pub(crate) fn scale(vector: [f64; 2], factor: f64) -> [f64; 2] {
    [vector[0] * factor, vector[1] * factor]
}

pub(crate) fn norm(vector: [f64; 2]) -> f64 {
    dot(vector, vector).sqrt()
}

fn outward_normal(side: BoundarySide, tangent: [f64; 2]) -> [f64; 2] {
    match side {
        BoundarySide::Top => [tangent[1], -tangent[0]],
        BoundarySide::Bottom => [-tangent[1], tangent[0]],
    }
}

#[cfg(test)]
mod tests {
    use crate::model::{BoundaryInterpolation, BoundaryShape, BoundaryShapePoint};

    use super::{BoundaryCurve, BoundarySide};

    #[test]
    fn extends_piecewise_linear_boundaries_and_uses_strict_node_selection() {
        let shape = BoundaryShape {
            interpolation: BoundaryInterpolation::PiecewiseLinear,
            points: vec![
                BoundaryShapePoint {
                    range_m: 0.0,
                    depth_m: 200.0,
                    material: None,
                },
                BoundaryShapePoint {
                    range_m: 4000.0,
                    depth_m: 0.0,
                    material: None,
                },
            ],
        };
        let curve = BoundaryCurve::new(BoundarySide::Bottom, 200.0, Some(&shape));
        assert_eq!(curve.segment_for_range(0.0), 0);
        assert_eq!(curve.segment_for_range(0.001), 1);
        let sloping = curve.segment(1);
        assert!((sloping.tangent[0] - 0.998_752_338_877_844_6).abs() < 1.0e-15);
        assert!((sloping.tangent[1] + 0.049_937_616_943_892_23).abs() < 1.0e-15);
    }

    #[test]
    fn inside_distance_is_positive_between_flat_boundaries() {
        let top = BoundaryCurve::new(BoundarySide::Top, 0.0, None);
        let bottom = BoundaryCurve::new(BoundarySide::Bottom, 100.0, None);
        assert!((top.signed_inside_distance([0.0, 40.0], 0) - 40.0).abs() < 1.0e-14);
        assert!((bottom.signed_inside_distance([0.0, 40.0], 0) - 60.0).abs() < 1.0e-14);
    }
}
