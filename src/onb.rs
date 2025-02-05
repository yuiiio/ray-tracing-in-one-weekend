use crate::vec3::{vec3_add, vec3_mul_b, Vector3};

#[derive(Clone)]
pub struct Onb {
    pub u: Vector3<f64>,
    pub v: Vector3<f64>,
    pub w: Vector3<f64>,
}

impl Onb {
    // build_from_w should take normalized vec
    pub fn build_from_w(w: &Vector3<f64>) -> Self {
        let (v, u) = if w[0].abs() > 0.9 {
            //let b = 1.0 / (w[2].powi(2) + w[0].powi(2)).sqrt();
            // w is normalized vec so, (w[2]^2 + w[0]^2) = (1 - w[1]^2)
            let b = 1.0 / (1.0 - w[1].powi(2)).sqrt();
            let v0 = -w[2] * b;
            let v2 = w[0] * b;
            (
                [v0, 0.0, v2],
                [w[1] * v2, w[2] * v0 - w[0] * v2, -w[1] * v0],
            )
        } else {
            let b = 1.0 / (1.0 - w[0].powi(2)).sqrt();
            let v1 = w[2] * b;
            let v2 = -w[1] * b;
            (
                [0.0, v1, v2],
                [w[1] * v2 - w[2] * v1, -w[0] * v2, w[0] * v1],
            )
        };

        Onb { u, v, w: *w }
    }

    pub fn local(&self, a: &Vector3<f64>) -> Vector3<f64> {
        vec3_add(
            &vec3_mul_b(&self.u, a[0]),
            &vec3_add(&vec3_mul_b(&self.v, a[1]), &vec3_mul_b(&self.w, a[2])),
        )
    }
}
