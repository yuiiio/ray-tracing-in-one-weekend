use rand::prelude::*;
use std::f64::consts::PI;
use std::sync::Arc;

use crate::hitable::{HitRecord, Hitable};
use crate::ray::Ray;
use crate::vec3::{
    cross, vec3_add, vec3_dot, vec3_mul_b, vec3_squared_length, vec3_sub, vec3_unit_vector_f64,
    Vector3,
};
use std::f64;

pub trait Pdf {
    fn value(&self, hit_record: &HitRecord, direction: &Vector3<f64>) -> f64;
    fn generate(&self, hit_record: &HitRecord) -> Vector3<f64>;
}

pub struct CosinePdf {}

impl Pdf for CosinePdf {
    fn value(&self, hit_record: &HitRecord, direction: &Vector3<f64>) -> f64 {
        let n = hit_record.get_normal(); //Already normalized?
        let direction = vec3_unit_vector_f64(*direction);
        let cosine = vec3_dot(n, direction);
        if cosine > 0.0 {
            return cosine / PI;
        } else {
            return 0.0;
        };
    }
    fn generate(&self, hit_record: &HitRecord) -> Vector3<f64> {
        let u = hit_record.get_normal(); //Already normalized?
        let a: Vector3<f64>;
        if u[0].abs() > 0.9 {
            a = [0.0, 1.0, 0.0];
        } else {
            a = [1.0, 0.0, 0.0];
        }
        let v = vec3_unit_vector_f64(cross(u, a));
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
