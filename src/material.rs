use rand::prelude::*;

use crate::ray::{Ray};
use crate::vec3::{Vector3, vec3_mul_b, vec3_dot, vec3_add, vec3_unit_vector_f64, vec3_sub, vec3_squared_length};
use crate::hitable::{HitRecord};

pub struct MatRecord {
    pub scatterd: Ray,
    pub attenuation: Vector3<f64>,
}

pub trait Material {
    fn scatter(&self, r: &Ray, hit_record: &HitRecord) -> Option<MatRecord>;
}

pub struct Metal {
    pub albedo: Vector3<f64>,
}

fn reflect(v: Vector3<f64>, n: Vector3<f64>) -> Vector3<f64> {
    vec3_add(vec3_mul_b(vec3_mul_b(n, vec3_dot(vec3_mul_b(v, -1.0), n)), 2.0) , v)
}

impl Material for Metal {
    fn scatter(&self, r_in: &Ray, hit_record: &HitRecord) -> Option<MatRecord> {
        let reflected = reflect(vec3_unit_vector_f64(r_in.direction()), hit_record.normal);
        let scatterd = Ray{ a: hit_record.p, b: reflected };
        let attenuation = self.albedo;
        if vec3_dot(scatterd.direction(), hit_record.normal) > 0.0 {
            Some(MatRecord{ scatterd, attenuation }) 
        } else {
            None
        }
    }
}

pub struct Lambertion {
    pub albedo: Vector3<f64>,
}

fn random_in_unit_sphere() -> Vector3<f64> {
    let mut rng = rand::thread_rng();
    loop {
        let rand_x: f64 = rng.gen();
        let rand_y: f64 = rng.gen();
        let rand_z: f64 = rng.gen();
        let p =  vec3_sub(vec3_mul_b( [rand_x, rand_y, rand_z] , 2.0 ), [1.0, 1.0, 1.0]);
        if vec3_squared_length(p) <= 1.0 * 1.0 {
            return p
        }
    }
}

impl Material for Lambertion {
    fn scatter(&self, _r_in: &Ray, hit_record: &HitRecord) -> Option<MatRecord> {
        let temp = vec3_add(vec3_add(hit_record.p, hit_record.normal), random_in_unit_sphere());
        let scatterd = Ray{ a: hit_record.p, b: vec3_sub(temp, hit_record.p )};
        let attenuation = self.albedo;
        Some(MatRecord{ scatterd, attenuation }) 
    }
}
