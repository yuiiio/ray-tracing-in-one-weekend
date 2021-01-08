use rand::prelude::*;
use std::f64::consts::PI;

use crate::hitable::HitRecord;
use crate::ray::Ray;
use crate::texture::Texture;
use crate::vec3::{
    vec3_add, vec3_dot, vec3_mul_b, vec3_squared_length, vec3_sub, vec3_unit_vector_f64, Vector3,
};
use std::f64;

pub struct MatRecord {
    scatterd: Ray,
    attenuation: Vector3<f64>,
    absorabance: Vector3<f64>,
}

impl MatRecord {
    pub fn get_scatterd(&self) -> &Ray {
        &self.scatterd
    }

    pub fn get_attenuation(&self) -> Vector3<f64> {
        self.attenuation
    }

    pub fn get_absorabance(&self) -> Vector3<f64> {
        self.absorabance
    }
}

pub trait Material {
    fn scatter(&self, r: &Ray, hit_record: &HitRecord) -> Option<MatRecord>;
    fn emitted(&self, _r: &Ray, _hit_record: &HitRecord) -> Vector3<f64> {
        [0.0, 0.0, 0.0]
    }
}

#[derive(Clone, Copy)]
pub struct MaterialHandle(pub usize);

pub struct Materials(Vec<Box<dyn Material + Send + Sync>>);

impl Materials {
    pub fn new() -> Self {
        Materials(Vec::new())
    }

    pub fn add_material<M: Material + 'static + Send + Sync>(
        &mut self,
        material: M,
    ) -> MaterialHandle {
        self.0.push(Box::new(material));
        MaterialHandle(self.0.len() - 1)
    }

    pub fn get(&self, handle: MaterialHandle) -> &dyn Material {
        self.0[handle.0].as_ref()
    }
}

pub struct Metal<T: Texture> {
    fuzz: f64,
    texture: T,
}

impl<T: Texture> Metal<T> {
    pub fn new(f: f64, texture: T) -> Self {
        let fuzz = if f < 1.0 { f } else { 1.0 };
        Metal { fuzz, texture }
    }
}

fn reflect(v: Vector3<f64>, n: Vector3<f64>) -> Vector3<f64> {
    vec3_add(
        vec3_mul_b(vec3_mul_b(n, vec3_dot(vec3_mul_b(v, -1.0), n)), 2.0),
        v,
    )
}

impl<T: Texture> Material for Metal<T> {
    fn scatter(&self, r_in: &Ray, hit_record: &HitRecord) -> Option<MatRecord> {
        let reflected = reflect(
            vec3_unit_vector_f64(r_in.direction()),
            hit_record.get_normal(),
        );
        let scatterd = Ray::new(
            hit_record.get_p(),
            vec3_add(reflected, vec3_mul_b(random_in_unit_sphere(), self.fuzz)),
        );
        let attenuation =
            self.texture
                .get_value(hit_record.get_u(), hit_record.get_v(), &hit_record.get_p());
        let absorabance = [0.0, 0.0, 0.0];
        if vec3_dot(scatterd.direction(), hit_record.get_normal()) > 0.0 {
            Some(MatRecord {
                scatterd,
                attenuation,
                absorabance,
            })
        } else {
            None
        }
    }
}

pub struct Lambertian<T: Texture> {
    texture: T,
}

impl<T: Texture> Lambertian<T> {
    pub fn new(texture: T) -> Self {
        Lambertian { texture }
    }
}

fn random_in_unit_sphere() -> Vector3<f64> {
    let mut rng = rand::thread_rng();
    loop {
        let rand_x: f64 = rng.gen();
        let rand_y: f64 = rng.gen();
        let rand_z: f64 = rng.gen();
        let p = vec3_sub(vec3_mul_b([rand_x, rand_y, rand_z], 2.0), [1.0, 1.0, 1.0]);
        if vec3_squared_length(p) <= 1.0 * 1.0 {
            return p;
        }
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

impl<T: Texture> Material for Lambertian<T> {
    fn scatter(&self, _r_in: &Ray, hit_record: &HitRecord) -> Option<MatRecord> {
        let temp = vec3_add(
            vec3_add(hit_record.get_p(), hit_record.get_normal()),
            random_in_unit_sphere(),
        );
        let scatterd = Ray::new(hit_record.get_p(), vec3_sub(temp, hit_record.get_p()));
        let attenuation =
            self.texture
                .get_value(hit_record.get_u(), hit_record.get_v(), &hit_record.get_p());
        let absorabance = [0.0, 0.0, 0.0];
        Some(MatRecord {
            scatterd,
            attenuation,
            absorabance,
        })
    }
}

pub struct Dielectric {
    ref_idx: f64,
    absorabance: Vector3<f64>,
}

impl Dielectric {
    pub fn new(ref_idx: f64, absorabance: Vector3<f64>) -> Self {
        Dielectric {
            ref_idx,
            absorabance,
        }
    }
}

fn refract(v: Vector3<f64>, n: Vector3<f64>, ni_over_nt: f64) -> Option<Vector3<f64>> {
    let uv = vec3_unit_vector_f64(vec3_mul_b(v, -1.0));
    let dt = vec3_dot(uv, n);
    let discriminant = 1.0 - ni_over_nt * ni_over_nt * (1.0 - dt * dt);
    if discriminant > 0.0 {
        let refracted = vec3_sub(
            vec3_mul_b(vec3_mul_b(n, -1.0), discriminant.sqrt()),
            vec3_mul_b(vec3_sub(uv, vec3_mul_b(n, dt)), ni_over_nt),
        );
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
        let outward_normal: Vector3<f64>;
        let ni_over_nt: f64;
        let attenuation: Vector3<f64> = [1.0, 1.0, 1.0];
        let scatterd: Ray;
        let reflect_prob: f64;
        let cosine: f64;
        let mut outside_to_inside: bool = false;
        if vec3_dot(vec3_mul_b(_r_in.direction(), -1.0), hit_record.get_normal()) > 0.0 {
            outside_to_inside = true;
            outward_normal = hit_record.get_normal();
            ni_over_nt = 1.0 / self.ref_idx;
            cosine = 1.0 * vec3_dot(vec3_mul_b(_r_in.direction(), -1.0), hit_record.get_normal());
        } else {
            outward_normal = vec3_mul_b(hit_record.get_normal(), -1.0);
            ni_over_nt = self.ref_idx / 1.0;
            cosine = self.ref_idx * vec3_dot(_r_in.direction(), hit_record.get_normal());
        }
        let mut inside_to_inside: bool = false;
        let refracted = match refract(_r_in.direction(), outward_normal, ni_over_nt) {
            Some(refracted) => {
                reflect_prob = schlick(cosine, self.ref_idx);
                Some(refracted)
            }
            None => {
                reflect_prob = 1.0;
                inside_to_inside = true;
                None
            }
        };

        let mut rng = rand::thread_rng();
        let rand: f64 = rng.gen();
        let mut refracted_root: bool = false;
        if rand < reflect_prob {
            scatterd = Ray::new(
                hit_record.get_p(),
                reflect(_r_in.direction(), outward_normal),
            );
        } else {
            refracted_root = true;
            scatterd = Ray::new(hit_record.get_p(), refracted.unwrap());
        }

        let mut absorabance = [0.0, 0.0, 0.0];

        if (outside_to_inside && refracted_root) || inside_to_inside {
            absorabance = self.absorabance;
        }

        Some(MatRecord {
            scatterd,
            attenuation,
            absorabance,
        })
    }
}

pub struct DiffuseLight<T: Texture> {
    texture: T,
}

impl<T: Texture> DiffuseLight<T> {
    pub fn new(texture: T) -> Self {
        DiffuseLight { texture }
    }
}

impl<T: Texture> Material for DiffuseLight<T> {
    fn scatter(&self, _r: &Ray, _hit_record: &HitRecord) -> Option<MatRecord> {
        None
    }

    fn emitted(&self, _r: &Ray, hit_record: &HitRecord) -> Vector3<f64> {
        self.texture
            .get_value(hit_record.get_u(), hit_record.get_v(), &hit_record.get_p())
    }
}
