use crate::vec3::{vec3_add, vec3_mul_b, Vector3};

pub struct Ray {
    pub origin: Vector3<f64>,
    pub direction: Vector3<f64>,
}

impl Ray {
    pub fn point_at_parameter(&self, t: f64) -> Vector3<f64> {
        vec3_add(&self.origin, &vec3_mul_b(&self.direction, t))
    }
    pub fn get_inv_dir(&self) -> Vector3<f64> {
        [
            1.0 / self.direction[0],
            1.0 / self.direction[1],
            1.0 / self.direction[2],
        ]
    }
}
