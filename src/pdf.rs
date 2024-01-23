use rand::prelude::*;
use std::f64::consts::PI;

use crate::onb::Onb;
use crate::vec3::{vec3_dot, Vector3};

pub fn cosine_pdf_value(hit_rec_normal: &Vector3<f64>, direction: &Vector3<f64>) -> f64 {
    let cosine = vec3_dot(hit_rec_normal, direction);
    if cosine.is_sign_positive() {
        cosine / PI
    } else {
        0.0
    }
}

// should return normalized vector
pub fn cosine_pdf_generate(uvw: &Onb) -> Vector3<f64> {
    let rcd = random_cosine_direction();
    uvw.local(&rcd)
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
