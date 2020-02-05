use crate::vec3::{Vector3, vec3_add, vec3_mul_b};

pub struct Ray {
    pub a: Vector3<f64>,
    pub b: Vector3<f64>,
}

impl Ray {
    pub fn origin(&self) -> Vector3<f64> {
        self.a
    }
    pub fn direction(&self) -> Vector3<f64> {
        self.b
    }

    pub fn point_at_parameter(&self, t: f64) -> Vector3<f64> {
        vec3_add(self.a, vec3_mul_b(self.b, t))
    }

}
