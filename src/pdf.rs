use rand::prelude::*;
use std::f64::consts::PI;

use crate::hitable::HitRecord;
use crate::ray::Ray;
use crate::vec3::{
    cross, vec3_add, vec3_dot, vec3_mul_b, vec3_squared_length, vec3_sub, vec3_unit_vector_f64,
    Vector3,
};
use std::f64;

pub trait Pdf {
    fn value(hit_record: &HitRecord, direction: &Vector3<f64>) -> f64;
    fn generate(hit_record: &HitRecord) -> Vector3<f64>;
}

pub struct CosinePdf {}

impl Pdf for CosinePdf {
    fn value(hit_record: &HitRecord, direction: &Vector3<f64>) -> f64 {
        let n = hit_record.get_normal(); //Already normalized?
        let direction = vec3_unit_vector_f64(*direction);
        return vec3_dot(n, direction) / PI;
    }
    fn generate(hit_record: &HitRecord) -> Vector3<f64> {
        let u = hit_record.get_normal(); //Already normalized?
        let v = vec3_unit_vector_f64(cross(u, [1.0, 0.0, 0.0]));
        let w = cross(v, u);

        let rcd = random_cosine_direction();
        vec3_add(
            vec3_mul_b(v, rcd[0]),
            vec3_add(vec3_mul_b(w, rcd[1]), vec3_mul_b(u, rcd[2])),
        )
    }
}

fn random_cosine_direction() -> Vector3<f64> {
    let mut rng = rand::thread_rng();
    let r1: f64 = rng.gen();
    let r2: f64 = rng.gen();

    let a = 2.0 * PI * r1;
    let b = r2.sqrt();
    let x: f64 = a.cos() * b;
    let y: f64 = a.sin() * b;
    let z: f64 = (1.0 - r2).sqrt();
    [x, y, z]
}
