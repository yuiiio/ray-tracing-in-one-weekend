use crate::vec3::{vec3_add, vec3_mul, vec3_mul_b, vec3_unit_vector_f64, Vector3};
use core::simd::num::SimdFloat;
use core::simd::Simd;

pub fn matrix_dot_3x4(
    left: &[Simd<f64, 4>; 3],
    right_transposed: &[Simd<f64, 4>; 3],
) -> [[f64; 3]; 3] {
    [
        [
            (left[0] * right_transposed[0]).reduce_sum(),
            (left[0] * right_transposed[1]).reduce_sum(),
            (left[0] * right_transposed[2]).reduce_sum(),
        ],
        [
            (left[1] * right_transposed[0]).reduce_sum(),
            (left[1] * right_transposed[1]).reduce_sum(),
            (left[1] * right_transposed[2]).reduce_sum(),
        ],
        [
            (left[2] * right_transposed[0]).reduce_sum(),
            (left[2] * right_transposed[1]).reduce_sum(),
            (left[2] * right_transposed[2]).reduce_sum(),
        ],
    ]
}

#[derive(Clone)]
pub struct Rotation {
    q: Qotation,
    oq: Qotation,
    pub left_side_matrix: [Vector3<f64>; 3],
}

impl Rotation {
    pub fn new(degrees: f64, axis: &Vector3<f64>) -> Self {
        let axis = vec3_unit_vector_f64(axis); // to unit vector
        let radians_div2 = degrees.to_radians() / 2.0;
        let cos = radians_div2.cos();
        let sin = radians_div2.sin();
        let vec3 = vec3_mul_b(&axis, sin);
        let q = Qotation { w: cos, xyz: vec3 };
        let oq = Qotation {
            w: cos,
            xyz: vec3_mul_b(&vec3, -1.0),
        };
        let q_left_matrix_3x4 = [
            Simd::from([q.xyz[0], q.w, -q.xyz[2], q.xyz[1]]),
            Simd::from([q.xyz[1], q.xyz[2], q.w, -q.xyz[0]]),
            Simd::from([q.xyz[2], -q.xyz[1], q.xyz[0], q.w]),
        ];
        let oq_right_matrix_4x3_transposed = [
            Simd::from([-oq.xyz[0], oq.w, -oq.xyz[2], oq.xyz[1]]),
            Simd::from([-oq.xyz[1], oq.xyz[2], oq.w, -oq.xyz[0]]),
            Simd::from([-oq.xyz[2], -oq.xyz[1], oq.xyz[0], oq.w]),
        ];
        // only need 3x3 matrix for rotation.
        let left_side_matrix = matrix_dot_3x4(&q_left_matrix_3x4, &oq_right_matrix_4x3_transposed);
        Rotation {
            q,
            oq,
            left_side_matrix,
        }
    }

    pub fn get_revq(&self) -> Self {
        let q = self.oq.clone();
        let oq = self.q.clone();
        let q_left_matrix_3x4 = [
            Simd::from([q.xyz[0], q.w, -q.xyz[2], q.xyz[1]]),
            Simd::from([q.xyz[1], q.xyz[2], q.w, -q.xyz[0]]),
            Simd::from([q.xyz[2], -q.xyz[1], q.xyz[0], q.w]),
        ];
        let oq_right_matrix_4x3_transposed = [
            Simd::from([-oq.xyz[0], oq.w, -oq.xyz[2], oq.xyz[1]]),
            Simd::from([-oq.xyz[1], oq.xyz[2], oq.w, -oq.xyz[0]]),
            Simd::from([-oq.xyz[2], -oq.xyz[1], oq.xyz[0], oq.w]),
        ];
        // only need 3x3 matrix for rotation.
        let left_side_matrix = matrix_dot_3x4(&q_left_matrix_3x4, &oq_right_matrix_4x3_transposed);
        Rotation {
            q,
            oq,
            left_side_matrix,
        }
    }

    pub fn rotate(&self, target: &Vector3<f64>) -> Vector3<f64> {
        // matrix dot
        // cost: 3x3= 9 mul & 2*3 = 6 add
        [
            vec3_mul(&self.left_side_matrix[0], target).iter().sum(),
            vec3_mul(&self.left_side_matrix[1], target).iter().sum(),
            vec3_mul(&self.left_side_matrix[2], target).iter().sum(),
        ]
    }
}

#[derive(Clone)]
pub struct Qotation {
    pub w: f64,
    pub xyz: [f64; 3],
}

#[allow(dead_code)]
pub fn qmul(q1: &Qotation, q2: &Qotation) -> Qotation {
    let w1 = q1.w;
    let w2 = q2.w;
    let v1 = q1.xyz;
    let v2 = q2.xyz;
    let inner = (v1[0] * v2[0]) + (v1[1] * v2[1]) + (v1[2] * v2[2]);
    let cross: Vector3<f64> = [
        (v1[1] * v2[2]) - (v1[2] * v2[1]),
        (v1[2] * v2[0]) - (v1[0] * v2[2]),
        (v1[0] * v2[1]) - (v1[1] * v2[0]),
    ];
    Qotation {
        w: (w1 * w2) - inner,
        xyz: vec3_add(
            &vec3_add(&vec3_mul_b(&v2, w1), &vec3_mul_b(&v1, w2)),
            &cross,
        ),
    }
}
