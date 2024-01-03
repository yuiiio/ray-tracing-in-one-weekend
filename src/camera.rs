use crate::ray::Ray;
use crate::vec3::{cross, vec3_add, vec3_dot, vec3_mul_b, vec3_sub, vec3_unit_vector_f64, Vector3};
use rand::prelude::*;
use std::f64;

pub struct Camera {
    origin: Vector3<f64>,
    lower_left_corner: Vector3<f64>,
    horizontal: Vector3<f64>,
    vertical: Vector3<f64>,
}

impl Camera {
    pub fn new(
        lookfrom: Vector3<f64>,
        lookat: &Vector3<f64>,
        vup: &Vector3<f64>,
        vfov: f64,
        aspect: f64,
    ) -> Self {
        let theta = vfov * ((2.0 * f64::consts::PI) / 360.0);
        let half_height = (theta / 2.0).tan() * 1.0;
        let half_width = half_height * aspect;

        let w = vec3_unit_vector_f64(&vec3_sub(&lookfrom, lookat));
        let u = vec3_unit_vector_f64(&cross(vup, &w));
        let v = cross(&w, &u); //already length = 1.0 * 1.0 * sin(90) = 1.0

        let origin = lookfrom;
        let lower_left_corner = vec3_sub(
            &vec3_add(
                &vec3_add(&origin, &vec3_mul_b(&u, -half_width)),
                &vec3_mul_b(&v, -half_height),
            ),
            &w,
        );
        let horizontal = vec3_mul_b(&u, 2.0 * half_width);
        let vertical = vec3_mul_b(&v, 2.0 * half_height);
        Camera {
            origin,
            lower_left_corner,
            horizontal,
            vertical,
        }
    }

    pub fn get_ray(&self, s: f64, t: f64) -> Ray {
        Ray::new(
            self.origin,
            vec3_sub(
                &vec3_add(
                    &vec3_add(&self.lower_left_corner, &vec3_mul_b(&self.horizontal, s)),
                    &vec3_mul_b(&self.vertical, t),
                ),
                &self.origin,
            ),
        )
    }
}

fn random_in_unit_disk() -> Vector3<f64> {
    let mut rng = rand::thread_rng();
    let rand_a: f64 = rng.gen();
    let rand_b: f64 = rng.gen();

    let sqrt_a: f64 = rand_a.sqrt();
    let seeta: f64 = rand_b * 2.0 * f64::consts::PI;
    let x: f64 = seeta.cos() * sqrt_a;
    let y: f64 = seeta.sin() * sqrt_a;

    return [x, y, 0.0];
}

mod test {
    use super::*;

    #[test]
    fn random_in_disk() {
        let i = 0;
        for i in 0..1000 {
            let pos = random_in_unit_disk();
            println!("{:?}", pos);
            let result = if vec3_dot(&pos, &pos) < 1.0 { 1 } else { 0 };
            assert_eq!(1, result);
        }
    }
}
