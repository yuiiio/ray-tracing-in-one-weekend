use rand::prelude::*;
use std::f64::consts::PI;

use crate::onb::Onb;
use crate::vec3::{vec3_dot, vec3_unit_vector_f64, Vector3};

pub fn cosine_pdf_value(hit_rec_normal: &Vector3<f64>, direction: &Vector3<f64>) -> f64 {
    let n = hit_rec_normal; //Already normalized?
    let direction = vec3_unit_vector_f64(direction); //just normalized
    let cosine = vec3_dot(&n, &direction);
    if cosine.is_sign_positive() {
        return cosine / PI;
    } else {
        return 0.0;
    };
}

// should return normalized vector
pub fn cosine_pdf_generate(hit_rec_normal: &Vector3<f64>) -> Vector3<f64> {
    let uvw = Onb::build_from_w(&hit_rec_normal);

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
