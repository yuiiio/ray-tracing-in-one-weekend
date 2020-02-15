use rand::prelude::*;

use crate::ray::{Ray};
use crate::vec3::{Vector3, vec3_mul_b, vec3_dot, vec3_add, vec3_unit_vector_f64, vec3_sub, vec3_squared_length};
use crate::hitable::{HitRecord};

pub struct MatRecord {
    scatterd: Ray,
    attenuation: Vector3<f64>,
}

impl MatRecord {
    pub fn get_scatterd(&self) -> &Ray {
        &self.scatterd
    }

    pub fn get_attenuation(&self) -> Vector3<f64> {
        self.attenuation
    }
}

pub trait Material {
    fn scatter(&self, r: &Ray, hit_record: &HitRecord) -> Option<MatRecord>;
}

pub struct Metal {
    albedo: Vector3<f64>,
    fuzz: f64,
}

impl Metal {
    pub fn new(albedo: Vector3<f64>, f: f64) -> Metal {
        let fuzz = if f < 1.0 { f } else { 1.0 };
        Metal{ albedo, fuzz }
    }
}

fn reflect(v: Vector3<f64>, n: Vector3<f64>) -> Vector3<f64> {
    vec3_add(vec3_mul_b(vec3_mul_b(n, vec3_dot(vec3_mul_b(v, -1.0), n)), 2.0) , v)
}

impl Material for Metal {
    fn scatter(&self, r_in: &Ray, hit_record: &HitRecord) -> Option<MatRecord> {
        let reflected = reflect(vec3_unit_vector_f64(r_in.direction()), hit_record.get_normal());
        let scatterd = Ray::new(hit_record.get_p(), vec3_add(reflected, vec3_mul_b(random_in_unit_sphere(), self.fuzz)));
        let attenuation = self.albedo;
        if vec3_dot(scatterd.direction(), hit_record.get_normal()) > 0.0 {
            Some(MatRecord{ scatterd, attenuation }) 
        } else {
            None
        }
    }
}

pub struct Lambertian {
    albedo: Vector3<f64>,
}

impl Lambertian {
    pub fn new(albedo: Vector3<f64>) -> Lambertian {
        Lambertian{ albedo }
    }
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

impl Material for Lambertian {
    fn scatter(&self, _r_in: &Ray, hit_record: &HitRecord) -> Option<MatRecord> {
        let temp = vec3_add(vec3_add(hit_record.get_p(), hit_record.get_normal()), random_in_unit_sphere());
        let scatterd = Ray::new(hit_record.get_p(), vec3_sub(temp, hit_record.get_p()));
        let attenuation = self.albedo;
        Some(MatRecord{ scatterd, attenuation }) 
    }
}
