use rand::prelude::*;
use std::f64::consts::PI;

use crate::hitable::HitRecord;
use crate::ray::Ray;
use crate::texture::Texture;
use crate::vec3::{
    vec3_add, vec3_dot, vec3_mul_b, vec3_squared_length, vec3_sub, vec3_unit_vector_f64, Vector3,
};
use std::f64;

pub enum Scatterd {
    Ray(Ray),
    CosinePdf,
}

pub struct MatRecord {
    scatterd: Scatterd,
    attenuation: Vector3<f64>,
    absorabance: Vector3<f64>,
}

impl MatRecord {
    pub fn get_scatterd(&self) -> &Scatterd {
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
    fn scattering_pdf(&self, _r: &Ray, _hit_record: &HitRecord) -> f64 {
        return 0.0;
    }
}

#[derive(Clone, Copy)]
pub struct MaterialHandle(pub usize);

pub struct Materials {
    material_list: Vec<Box<dyn Material + Send + Sync>>,
}

impl Materials {
    pub fn new() -> Self {
        Materials{
            material_list: Vec::new(),
        }
    }

    pub fn add_material<M: Material + 'static + Send + Sync>(
        &mut self,
        material: M,
    ) -> MaterialHandle {
        self.material_list.push(Box::new(material));
        MaterialHandle(self.material_list.len() - 1)
    }

    pub fn get(&self, handle: MaterialHandle) -> &dyn Material {
        self.material_list[handle.0].as_ref()
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

fn reflect(v: &Vector3<f64>, n: &Vector3<f64>) -> Vector3<f64> {
    vec3_add(
        &vec3_mul_b(&vec3_mul_b(n, vec3_dot(&vec3_mul_b(v, -1.0), n)), 2.0),
        &v,
    )
}

impl<T: Texture> Material for Metal<T> {
    fn scatter(&self, r_in: &Ray, hit_record: &HitRecord) -> Option<MatRecord> {
        let mut reflected = reflect(
            &vec3_unit_vector_f64(&r_in.direction()),
            &hit_record.get_normal(),
        );
        if self.fuzz != 0.0 { // should return Pdf instadof Ray when fuzz != 0.0 ?
            reflected = vec3_add(&reflected, &vec3_mul_b(&random_in_unit_sphere(), self.fuzz));
        }
        let scatterd = Ray::new(hit_record.get_p(), reflected);
        let attenuation =
            self.texture
                .get_value(hit_record.get_u(), hit_record.get_v(), &hit_record.get_p());
        let absorabance = [0.0, 0.0, 0.0];
        if vec3_dot(&scatterd.direction(), &hit_record.get_normal()).is_sign_positive() {
            Some(MatRecord {
                scatterd: Scatterd::Ray(scatterd),
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
    let r1: f64 = rng.gen();
    let r2: f64 = rng.gen();
    let r3: f64 = rng.gen();

    let a = 2.0 * PI * r1;
    let b = 2.0 * (r2 * (1.0 - r2)).sqrt();

    let c = r3.cbrt();

    let x = a.cos() * b * c;
    let y = b.sin() * b * c;
    let z = (1.0 - (2.0 * r2)) * c;
    [x, y, z]
}

impl<T: Texture> Material for Lambertian<T> {
    fn scatter(&self, _r_in: &Ray, hit_record: &HitRecord) -> Option<MatRecord> {
        let attenuation =
            self.texture
                .get_value(hit_record.get_u(), hit_record.get_v(), &hit_record.get_p());
        let absorabance = [0.0, 0.0, 0.0];
        Some(MatRecord {
            scatterd: Scatterd::CosinePdf,
            attenuation,
            absorabance,
        })
    }
    fn scattering_pdf(&self, r: &Ray, _hit_record: &HitRecord) -> f64 {
        let n = _hit_record.get_normal(); //Already normalized?
        let direction = vec3_unit_vector_f64(&r.direction());
        let cosine = vec3_dot(&n, &direction);
        if cosine.is_sign_positive() {
            return cosine / PI;
        } else {
            return 0.0;
        };
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

fn refract(v: &Vector3<f64>, n: &Vector3<f64>, ni_over_nt: f64) -> Option<Vector3<f64>> {
    let uv = vec3_unit_vector_f64(&vec3_mul_b(v, -1.0));
    let dt = vec3_dot(&uv, n);
    let discriminant = 1.0 - ni_over_nt * ni_over_nt * (1.0 - dt * dt);
    if discriminant.is_sign_positive() {
        let refracted = vec3_sub(
            &vec3_mul_b(&vec3_mul_b(n, -1.0), discriminant.sqrt()),
            &vec3_mul_b(&vec3_sub(&uv, &vec3_mul_b(n, dt)), ni_over_nt),
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
    fn scatter(&self, r_in: &Ray, hit_record: &HitRecord) -> Option<MatRecord> {
        let outward_normal: Vector3<f64>;
        let ni_over_nt: f64;
        let attenuation: Vector3<f64> = [1.0, 1.0, 1.0];
        let scatterd: Ray;
        let reflect_prob: f64;
        let cosine: f64;
        let mut outside_to_inside: bool = false;
        let hit_normal = hit_record.get_normal();
        if vec3_dot(&vec3_mul_b(&r_in.direction(), -1.0), &hit_normal).is_sign_positive() {
            outside_to_inside = true;
            outward_normal = hit_record.get_normal();
            ni_over_nt = 1.0 / self.ref_idx;
            cosine = 1.0 * vec3_dot(&vec3_mul_b(&r_in.direction(), -1.0), &hit_normal);
        } else {
            outward_normal = vec3_mul_b(&hit_normal, -1.0);
            ni_over_nt = self.ref_idx / 1.0;
            cosine = self.ref_idx * vec3_dot(&r_in.direction(), &hit_normal);
        }
        let mut inside_to_inside: bool = false;
        let r_in_direction = r_in.direction();
        let refracted = match refract(&r_in_direction, &outward_normal, ni_over_nt) {
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
                reflect(&r_in_direction, &outward_normal),
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
            scatterd: Scatterd::Ray(scatterd),
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

    fn emitted(&self, r: &Ray, hit_record: &HitRecord) -> Vector3<f64> {
        if vec3_dot(&hit_record.get_normal(), &r.direction()).is_sign_negative() {
            return self.texture.get_value(
                hit_record.get_u(),
                hit_record.get_v(),
                &hit_record.get_p(),
            );
        } else {
            return [0.0, 0.0, 0.0];
        }
    }
}
