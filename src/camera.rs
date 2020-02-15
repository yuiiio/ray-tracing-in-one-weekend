use crate::ray::{Ray};
use crate::vec3::{Vector3, vec3_add, vec3_mul_b};

pub struct Camera {
    origin: Vector3<f64>,
    lower_left_corner: Vector3<f64>,
    horizontal: Vector3<f64>,
    vertical: Vector3<f64>,
}

impl Camera {
    pub fn new() -> Camera {
        let origin = [0.0, 0.0, 0.0];
        let lower_left_corner = [-2.0, -1.0, -1.0];
        let horizontal = [4.0, 0.0, 0.0];
        let vertical = [0.0, 2.0, 0.0];
        Camera {origin, lower_left_corner, horizontal, vertical}
    }

    pub fn get_ray(&self, u: f64, v: f64) -> Ray {
        Ray::new(self.origin, vec3_add(vec3_add(self.lower_left_corner, vec3_mul_b(self.horizontal, u)), vec3_mul_b(self.vertical, v)))
    }
}
