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
        let q = Qotation(cos, vec3);
        let oq = Qotation(cos, vec3_mul_b(&vec3, -1.0));
        Rotation {q, oq}
    }

    pub fn rotate(&self, target: &Vector3<f64>) -> Vector3<f64> {
        let x = qmul(&self.q, &Qotation(0.0, *target));
        let y = qmul(&x, &self.oq);
        y.1
    }
}

#[derive(Clone)]
pub struct Qotation(pub f64, pub [f64; 3]); //[w] [x, y, z]

pub fn qmul(q1: &Qotation, q2: &Qotation) -> Qotation {
    let w1 = q1.0;
    let w2 = q2.0;
    let v1 = q1.1;
    let v2 = q2.1;
    let inner = (v1[0] * v2[0])
                + (v1[1] * v2[1])
                + (v1[2] * v2[2]);
    let cross: Vector3<f64> = [ (v1[1] * v2[2]) - (v1[2] * v2[1]),
                                (v1[2] * v2[0]) - (v1[0] * v2[2]),
                                (v1[0] * v2[1]) - (v1[1] * v2[0]), ];
    Qotation((w1 * w2) - inner,
                vec3_add(&vec3_add(&vec3_mul_b(&v2, w1), &vec3_mul_b(&v1, w2)), &cross))
}
