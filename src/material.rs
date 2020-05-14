use rand::prelude::*;

use crate::ray::{Ray};
use crate::vec3::{Vector3, vec3_mul_b, vec3_dot, vec3_add, vec3_unit_vector_f64, vec3_sub, vec3_squared_length};
use crate::hitable::{HitRecord};
use std::f64;

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

pub struct MaterialHandle(pub usize);

pub struct Materials(Vec<Box<dyn Material>>);

impl Materials {
    pub fn new() -> Materials {
        Materials(Vec::new())
    }

    pub fn add_material<M: Material + 'static>(&mut self, material: M) -> MaterialHandle {
        self.0.push(Box::new(material));
        MaterialHandle(self.0.len() -1)
    }

    pub fn get(&self, handle: MaterialHandle) -> &dyn Material {
        self.0[handle.0].as_ref()
    }
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

pub struct Dielectric {
    ref_idx: f64,
}

impl Dielectric {
    pub fn new(ref_idx: f64) -> Dielectric {
        Dielectric{ ref_idx }
    }
}

fn refract(v: Vector3<f64>, n: Vector3<f64>, ni_over_nt: f64) -> Option<Vector3<f64>> {
    let uv = vec3_unit_vector_f64(vec3_mul_b(v, -1.0));
    let dt = vec3_dot(uv, n);
    let discriminant = 1.0 - ni_over_nt * ni_over_nt * ( 1.0 - dt * dt );
    if discriminant > 0.0 {
        let refracted = vec3_sub(vec3_mul_b(vec3_mul_b(n, -1.0), discriminant.sqrt()), vec3_mul_b(vec3_sub(uv, vec3_mul_b(n, dt)), ni_over_nt));
        Some(refracted)
    } else {
        None
    }
}

fn schlick(cosine: f64, ref_idx: f64) -> f64 {
    let r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
    let r0 = r0 * r0;
    r0 + ((1.0 - r0) * (1.0 - cosine).powi(5))
}

impl Material for Dielectric {
    fn scatter(&self, _r_in: &Ray, hit_record: &HitRecord) -> Option<MatRecord> {
        let outward_normal :Vector3<f64>;
        let ni_over_nt :f64;
        let attenuation :Vector3<f64> = [1.0, 1.0, 1.0];
        let scatterd :Ray;
        let reflect_prob :f64;
        let cosine :f64;
        if vec3_dot(vec3_mul_b(_r_in.direction(), -1.0), hit_record.get_normal()) > 0.0 {
            outward_normal = hit_record.get_normal();
            ni_over_nt = 1.0 / self.ref_idx;
            cosine = 1.0 * vec3_dot(vec3_mul_b(_r_in.direction(), -1.0), hit_record.get_normal());
        } else {
            outward_normal = vec3_mul_b(hit_record.get_normal(), -1.0);
            ni_over_nt = self.ref_idx / 1.0;
            cosine = self.ref_idx * vec3_dot(_r_in.direction(), hit_record.get_normal());
        }
        let refracted = match refract(_r_in.direction(), outward_normal, ni_over_nt) {
            Some(refracted) => {
                reflect_prob = schlick(cosine, self.ref_idx);
                Some(refracted)
            },
            None => {
                reflect_prob = 1.0;
                None
            },
        };

        let mut rng = rand::thread_rng();
        let rand :f64 = rng.gen();
        if rand < reflect_prob {
                scatterd = Ray::new(hit_record.get_p(), reflect(_r_in.direction(), outward_normal));
        } else {
                scatterd = Ray::new(hit_record.get_p(), refracted.unwrap());
        }

        Some(MatRecord{ scatterd, attenuation })
    }
}
