use crate::vec3::{cross, vec3_add, vec3_mul_b, Vector3};

pub struct Onb {
    u: Vector3<f64>,
    v: Vector3<f64>,
    w: Vector3<f64>,
}

impl Onb {
    // build_from_w should take normalized vec
    pub fn build_from_w(w: &Vector3<f64>) -> Self {
        let v = if w[0].abs() > 0.9 {
            //let b = 1.0 / (w[2].powi(2) + w[0].powi(2));
            // w is normalized vec so, (w[2]^2 + w[0]^2) = (1 - w[1]^2)
            let b = 1.0 / (1.0 - w[1].powi(2));
            [
                -w[2] * b,
                0.0,
                w[0] * b,
            ]
        } else {
            let b = 1.0 / (1.0 - w[0].powi(2));
            [   
                0.0,
                w[2] * b,
                -w[1] * b,
            ]
        };
        let u = cross(&w, &v);
        Onb { u, v, w: *w }
    }

    pub fn local(&self, a: &Vector3<f64>) -> Vector3<f64> {
        vec3_add(
            &vec3_mul_b(&self.u, a[0]),
            &vec3_add(&vec3_mul_b(&self.v, a[1]), &vec3_mul_b(&self.w, a[2])),
        )
    }
}
