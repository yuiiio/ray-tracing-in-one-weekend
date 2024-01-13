use crate::vec3::{Vector3, vec3_add, vec3_mul_b, vec3_unit_vector_f64};

#[derive(Clone)]
pub struct Rotation {
    q: Qotation,
    oq: Qotation,
}

impl Rotation {
    pub fn new(degrees: f64, axis: &Vector3<f64>) -> Self {
        let axis = vec3_unit_vector_f64(axis); // to unit vector
        let radians: f64 = degrees * std::f64::consts::PI / 180.0 / 2.0;
        let cos = radians.cos();
        let sin = radians.sin();
        let vec3 = vec3_mul_b(&axis, sin);
        let q = Qotation{ w: cos, xyz: vec3 };
        let oq = Qotation { w:cos, xyz: vec3_mul_b(&vec3, -1.0) };
        Rotation {q, oq}
    }

    pub fn rotate(&self, target: &Vector3<f64>) -> Vector3<f64> {
        let x = qmul_q2_xyz(&self.q, target);
        qmul_ret_xyz(&x, &self.oq)
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
    let inner = (v1[0] * v2[0])
                + (v1[1] * v2[1])
                + (v1[2] * v2[2]);
    let cross: Vector3<f64> = [ (v1[1] * v2[2]) - (v1[2] * v2[1]),
                                (v1[2] * v2[0]) - (v1[0] * v2[2]),
                                (v1[0] * v2[1]) - (v1[1] * v2[0]), ];
    Qotation {
        w: (w1 * w2) - inner,
        xyz: vec3_add(&vec3_add(&vec3_mul_b(&v2, w1), &vec3_mul_b(&v1, w2)), &cross),
    }
}

pub fn qmul_q2_xyz(q1: &Qotation, v2: &[f64; 3]) -> Qotation {
    let w1 = q1.w;
    let v1 = q1.xyz;
    let inner = (v1[0] * v2[0])
                + (v1[1] * v2[1])
                + (v1[2] * v2[2]);
    let cross: Vector3<f64> = [ (v1[1] * v2[2]) - (v1[2] * v2[1]),
                                (v1[2] * v2[0]) - (v1[0] * v2[2]),
                                (v1[0] * v2[1]) - (v1[1] * v2[0]), ];
    Qotation {
        w: - inner,
        xyz: vec3_add(&vec3_mul_b(&v2, w1), &cross),
    }
}

pub fn qmul_ret_xyz(q1: &Qotation, q2: &Qotation) -> [f64; 3] {
    let w1 = q1.w;
    let w2 = q2.w;
    let v1 = q1.xyz;
    let v2 = q2.xyz;
    let cross: Vector3<f64> = [ (v1[1] * v2[2]) - (v1[2] * v2[1]),
                                (v1[2] * v2[0]) - (v1[0] * v2[2]),
                                (v1[0] * v2[1]) - (v1[1] * v2[0]), ];

    vec3_add(&vec3_add(&vec3_mul_b(&v2, w1), &vec3_mul_b(&v1, w2)), &cross)
}
