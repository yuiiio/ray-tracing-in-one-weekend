use crate::vec3::{cross, vec3_add, vec3_mul_b, vec3_unit_vector_f64, Vector3};

pub struct Onb {
    u: Vector3<f64>,
    v: Vector3<f64>,
    w: Vector3<f64>,
}

impl Onb {
    // build_from_w should take normalized vec
    pub fn build_from_w(w: &Vector3<f64>) -> Self {
        let wa_cross = if w[0].abs() > 0.9 {
            [
                -w[2],
                0.0,
                w[0],
            ]
        } else {
            [   
                0.0,
                w[2],
                -w[1],
            ]
        };
        let v = vec3_unit_vector_f64(&wa_cross);
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
