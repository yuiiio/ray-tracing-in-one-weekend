use crate::ray::{Ray};
use crate::vec3::{Vector3, vec3_add, vec3_mul_b, vec3_sub, cross, vec3_unit_vector_f64, vec3_dot};
use std::f64;
use rand::prelude::*;

pub struct Camera {
    origin: Vector3<f64>,
    lower_left_corner: Vector3<f64>,
    horizontal: Vector3<f64>,
    vertical: Vector3<f64>,
}

impl Camera {
    pub fn new(lookfrom: Vector3<f64>, lookat: Vector3<f64>, vup: Vector3<f64>, vfov :f64, aspect :f64) -> Self {
        let theta = vfov * ((2.0 * f64::consts::PI) / 360.0);
        let half_height = (theta / 2.0).tan() * 1.0;
        let half_width = half_height * aspect;

        let w = vec3_unit_vector_f64(vec3_sub(lookfrom, lookat));
        let u = vec3_unit_vector_f64(cross(vup, w));
        let v = cross(w, u); //already length = 1.0 * 1.0 * sin(90) = 1.0

        let origin = lookfrom;
        let lower_left_corner = vec3_sub(
            vec3_add(vec3_add(origin, vec3_mul_b(u, -half_width)), vec3_mul_b(v, -half_height)) ,
            w);
        let horizontal = vec3_mul_b(u, 2.0*half_width);
        let vertical = vec3_mul_b(v, 2.0*half_height);
        Camera {origin, lower_left_corner, horizontal, vertical}
    }

    pub fn get_ray(&self, s: f64, t: f64) -> Ray {
        Ray::new(self.origin,
            vec3_sub(vec3_add(vec3_add(self.lower_left_corner, vec3_mul_b(self.horizontal, s)), vec3_mul_b(self.vertical, t)), self.origin)
        )
    }
}

fn random_in_unit_disk() -> Vector3<f64> {
    let mut rng = rand::thread_rng();
    let mut p = [0.0, 0.0, 0.0];
    loop {
        let random_x: f64 = rng.gen();
        let random_y: f64 = rng.gen();
        p = [
            (random_x * 2.0) - 1.0,
            (random_y * 2.0) - 1.0,
            0.0,
            ];
        if vec3_dot(p, p) < 1.0 {
            break;
        }
    }
    return p;
}

mod test {
    use super::*;

    #[test]
    fn random_in_disk() {
        let pos = random_in_unit_disk();
        println!("{:?}", pos);
        let result = if vec3_dot(pos, pos) < 1.0 {
            1
        } else {
            0
        };
        assert_eq!(1, result);
    }
}