use rand::prelude::*;
use std::f64::consts::PI;
use std::sync::Arc;

use crate::hitable::{HitRecord, Hitable};
use crate::bvh_node::BvhTree;
use crate::onb::Onb;
use crate::vec3::{vec3_dot, vec3_unit_vector_f64, Vector3};

pub fn cosine_pdf_value(hit_record: &HitRecord, direction: &Vector3<f64>) -> f64 {
    let n = hit_record.normal; //Already normalized?
    let direction = vec3_unit_vector_f64(direction); //just normalized
    let cosine = vec3_dot(&n, &direction);
    if cosine.is_sign_positive() {
        return cosine / PI;
    } else {
        return 0.0;
    };
}

pub fn cosine_pdf_generate(hit_record: &HitRecord) -> Vector3<f64> {
    let uvw = Onb::build_from_w(&hit_record.normal);

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

pub fn mix_cosine_pdf_value(pdf0: &Arc<BvhTree>, hit_record: &HitRecord, direction: &Vector3<f64>) -> f64 {
    let pdf0_value = pdf0.pdf_value(&hit_record.p, direction);
    let pdf1_value = cosine_pdf_value(hit_record, direction);
    return 0.5 * pdf0_value + 0.5 * pdf1_value;
}

pub fn mix_cosine_pdf_generate(pdf0: &Arc<BvhTree>, hit_record: &HitRecord) -> Vector3<f64> {
    let mut rng = rand::thread_rng();
    let r: f64 = rng.gen();
    if r < 0.5 {
        return pdf0.random(&hit_record.p);
    } else {
        return cosine_pdf_generate(hit_record);
    }
}
