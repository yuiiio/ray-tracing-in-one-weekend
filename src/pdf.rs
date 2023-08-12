use rand::prelude::*;
use std::f64::consts::PI;
use std::sync::Arc;

use crate::hitable::{HitRecord, Hitable};
use crate::onb::Onb;
use crate::vec3::{vec3_dot, vec3_unit_vector_f64, Vector3};

pub trait Pdf {
    fn value(&self, hit_record: &HitRecord, direction: &Vector3<f64>) -> f64;
    fn generate(&self, hit_record: &HitRecord) -> Vector3<f64>;
}

pub struct CosinePdf {}

impl Pdf for CosinePdf {
    fn value(&self, hit_record: &HitRecord, direction: &Vector3<f64>) -> f64 {
        let n = hit_record.get_normal(); //Already normalized?
        let direction = vec3_unit_vector_f64(*direction); //just normalized
        let cosine = vec3_dot(n, direction);
        if cosine.is_sign_positive() {
            return cosine / PI;
        } else {
            return 0.0;
        };
    }
    fn generate(&self, hit_record: &HitRecord) -> Vector3<f64> {
        let uvw = Onb::build_from_w(&hit_record.get_normal());

        let rcd = random_cosine_direction();

        uvw.local(&rcd)
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

pub struct HitablePdf<'a, T: Hitable> {
    pub hitable: &'a Arc<T>,
}

impl<'a, T: Hitable> Pdf for HitablePdf<'a, T> {
    fn value(&self, hit_record: &HitRecord, direction: &Vector3<f64>) -> f64 {
        return self.hitable.pdf_value(&hit_record.get_p(), direction);
    }
    fn generate(&self, hit_record: &HitRecord) -> Vector3<f64> {
        self.hitable.random(&hit_record.get_p())
    }
}

pub struct MixturePdf<'a, T: Pdf> {
    pub pdf0: T,
    pub pdf1: &'a Box<dyn Pdf>,
}

impl<'a, T: Pdf> Pdf for MixturePdf<'a, T> {
    fn value(&self, hit_record: &HitRecord, direction: &Vector3<f64>) -> f64 {
        let pdf0_value = self.pdf0.value(hit_record, direction);
        let pdf1_value = self.pdf1.value(hit_record, direction);
        return 0.5 * pdf0_value + 0.5 * pdf1_value;
    }
    fn generate(&self, hit_record: &HitRecord) -> Vector3<f64> {
        let mut rng = rand::thread_rng();
        let r: f64 = rng.gen();
        if r < 0.5 {
            return self.pdf0.generate(hit_record);
        } else {
            return self.pdf1.generate(hit_record);
        }
    }
}
