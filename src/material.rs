use rand::prelude::*;
use std::f64::consts::PI;

use crate::hitable::HitRecord;
use crate::ray::Ray;
use crate::texture::{TextureHandle, TextureList};
use crate::vec3::{
    vec3_add, vec3_dot, vec3_mul_b, vec3_sub, vec3_unit_vector_f64, Vector3,
};
use std::f64;

pub enum Scatterd {
    Ray(Ray),
    CosinePdf,
}

pub struct MatRecord {
    pub scatterd: Scatterd,
    pub attenuation: Vector3<f64>,
    pub absorabance: Vector3<f64>,
}

#[derive(Clone)]
pub enum Material {
    Metal,
    Lambertian,
    Dielectric,
    DiffuseLight,
}

#[derive(Clone)]
pub struct MaterialHandle {
    pub material_type: Material,
    position: usize, // each type
    pub needs_uv: bool,
}

pub struct MaterialList {
    metal_list: Vec<Metal>,
    lambertian_list: Vec<Lambertian>,
    dielectric_list: Vec<Dielectric>,
    diffuselight_list: Vec<DiffuseLight>,
}

impl MaterialList {
    pub fn new() -> Self {
        MaterialList{
            metal_list: Vec::new(),
            lambertian_list: Vec::new(),
            dielectric_list: Vec::new(),
            diffuselight_list: Vec::new(),
        }
    }

    pub fn add_metal_mat(&mut self, material: Metal) -> MaterialHandle {
        let pos = self.metal_list.len();
        let needs_uv = material.texture.needs_uv;
        self.metal_list.push(material);
        MaterialHandle {
            material_type: Material::Metal,
            position: pos,
            needs_uv,
        }
    }
    pub fn add_lambertian_mat(&mut self, material: Lambertian) -> MaterialHandle {
        let pos = self.lambertian_list.len();
        let needs_uv = material.texture.needs_uv;
        self.lambertian_list.push(material);
        MaterialHandle {
            material_type: Material::Lambertian,
            position: pos,
            needs_uv,
        }
    }
    pub fn add_dielectric_mat(&mut self, material: Dielectric) -> MaterialHandle {
        let pos = self.dielectric_list.len();
        self.dielectric_list.push(material);
        MaterialHandle {
            material_type: Material::Dielectric,
            position: pos,
            needs_uv: false,
        }
    }
    pub fn add_diffuselight_mat(&mut self, material: DiffuseLight) -> MaterialHandle {
        let pos = self.diffuselight_list.len();
        let needs_uv = material.texture.needs_uv;
        self.diffuselight_list.push(material);
        MaterialHandle {
            material_type: Material::DiffuseLight,
            position: pos,
            needs_uv,
        }
    }

    pub fn scatter(&self, r: &Ray, hit_record: &HitRecord, texture_list: &TextureList) -> Option<MatRecord> {
        let mat_pos = hit_record.mat_ptr.position;
        match hit_record.mat_ptr.material_type {
            Material::Metal => self.metal_list[mat_pos].scatter(r, hit_record, texture_list),
            Material::Lambertian => self.lambertian_list[mat_pos].scatter(r, hit_record, texture_list),
            Material::Dielectric => self.dielectric_list[mat_pos].scatter(r, hit_record, texture_list),
            Material::DiffuseLight => self.diffuselight_list[mat_pos].scatter(r, hit_record, texture_list),
        }
    }

    pub fn emitted(&self, r: &Ray, hit_record: &HitRecord, texture_list: &TextureList) -> Vector3<f64> {
        let mat_pos = hit_record.mat_ptr.position;
        match hit_record.mat_ptr.material_type {
            Material::Metal => [0.0, 0.0, 0.0], //self.metal_list[mat_pos].emitted(r, hit_record, texture_list),
            Material::Lambertian => [0.0, 0.0, 0.0], // self.lambertian_list[mat_pos].emitted(r, hit_record, texture_list),
            Material::Dielectric => [0.0, 0.0, 0.0], // self.dielectric_list[mat_pos].emitted(r, hit_record, texture_list),
            Material::DiffuseLight => self.diffuselight_list[mat_pos].emitted(r, hit_record, texture_list),
        }
    }

    pub fn scattering_pdf(&self, r: &Ray, hit_record: &HitRecord) -> f64 {
        let mat_pos = hit_record.mat_ptr.position;
        //should naver call on metal, Dielectric, DiffuseLight doesn't have pdf
        //scatter dpesn't return Pdf type.
        match hit_record.mat_ptr.material_type {
            Material::Metal => 0.0, //self.metal_list[mat_pos].scattering_pdf(r, hit_record),
            Material::Lambertian => self.lambertian_list[mat_pos].scattering_pdf(r, hit_record),
            Material::Dielectric => 0.0, //self.dielectric_list[mat_pos].scattering_pdf(r, hit_record),
            Material::DiffuseLight => 0.0 //self.diffuselight_list[mat_pos].scattering_pdf(r, hit_record),
        }
    }
}

pub struct Metal {
    fuzz: f64,
    texture: TextureHandle,
}

fn reflect(v: &Vector3<f64>, n: &Vector3<f64>) -> Vector3<f64> {
    vec3_add(
        &vec3_mul_b(&vec3_mul_b(n, vec3_dot(&v, n)), -2.0),
        v,
    )
}

impl Metal {
    pub fn new(f: f64, texture: TextureHandle) -> Self {
        let fuzz = if f < 1.0 { f } else { 1.0 };
        Metal { fuzz, texture }
    }

    fn scatter(&self, r_in: &Ray, hit_record: &HitRecord, texture_list: &TextureList) -> Option<MatRecord> {
        let mut reflected = reflect(
            &r_in.direction,
            &hit_record.normal,
        );
        if self.fuzz != 0.0 { // should return Pdf instadof Ray when fuzz != 0.0 ?
            reflected = vec3_unit_vector_f64(
                &vec3_add(&reflected, &vec3_mul_b(&random_in_unit_sphere(), self.fuzz))
                );
        }
        let scatterd = Ray{ origin: hit_record.p, direction: reflected };
        let attenuation = texture_list
            .get_value(hit_record.uv, &hit_record.p, &self.texture);
        let absorabance = [0.0, 0.0, 0.0];
        if vec3_dot(&scatterd.direction, &hit_record.normal).is_sign_positive() {
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

pub struct Lambertian {
    texture: TextureHandle,
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


impl Lambertian {
    pub fn new(texture: TextureHandle) -> Self {
        Lambertian { texture }
    }

    fn scatter(&self, _r_in: &Ray, hit_record: &HitRecord, texture_list: &TextureList) -> Option<MatRecord> {
        let attenuation = texture_list
            .get_value(hit_record.uv, &hit_record.p, &self.texture);
        let absorabance = [0.0, 0.0, 0.0];
        Some(MatRecord {
            scatterd: Scatterd::CosinePdf,
            attenuation,
            absorabance,
        })
    }
    fn scattering_pdf(&self, r: &Ray, _hit_record: &HitRecord) -> f64 {
        let n = _hit_record.normal;
        let cosine = vec3_dot(&n, &r.direction);
        if cosine.is_sign_positive() {
            cosine / PI
        } else {
            0.0
        }
    }
}

pub struct Dielectric {
    ref_idx: f64,
    ref_idx_seq: f64,
    nor_ref_idx: f64,
    nor_ref_idx_seq: f64,
    absorabance: Vector3<f64>,
    schlick_r0: f64,//option for realistic glass trick
    schlick_r1: f64,
}

fn refract(v: &Vector3<f64>, n: &Vector3<f64>, ni_over_nt: f64, ni_over_nt_seq: f64) -> Option<Vector3<f64>> {
    let dt = vec3_dot(&v, n);
    let discriminant = 1.0 - ni_over_nt_seq * (1.0 - dt * dt);
    if discriminant.is_sign_positive() {
        let refracted = vec3_sub(
            &vec3_mul_b(&n, -1.0 * discriminant.sqrt()),
            &vec3_mul_b(&vec3_sub(&vec3_mul_b(n, dt), &v), ni_over_nt),
        );
        Some(refracted)
    } else {
        None
    }
}

fn schlick(cosine: f64, r0: f64, r1: f64) -> f64 {
    r0 + (r1 * (1.0 - cosine).powi(5))
}

impl Dielectric {
    pub fn new(ref_idx: f64, absorabance: Vector3<f64>) -> Self {
        let nor_ref_idx = 1.0 / ref_idx;
        // for schlick
        let r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
        let schlick_r0 = r0 * r0;
        let schlick_r1 = 1.0 - schlick_r0;

        Dielectric {
            ref_idx,
            ref_idx_seq: ref_idx * ref_idx,
            nor_ref_idx,
            nor_ref_idx_seq: nor_ref_idx * nor_ref_idx,
            absorabance,
            schlick_r0,
            schlick_r1,
        }
    }

    fn scatter(&self, r_in: &Ray, hit_record: &HitRecord, _texture_list: &TextureList) -> Option<MatRecord> {
        let outward_normal: Vector3<f64>;
        let ni_over_nt: f64;
        let ni_over_nt_seq: f64;
        let attenuation: Vector3<f64> = [1.0, 1.0, 1.0];
        let reflect_prob: f64;
        let cosine: f64 = vec3_dot(&r_in.direction, &hit_record.normal);
        let mut outside_to_inside: bool = false;
        let cosine = if cosine.is_sign_negative() {
            outside_to_inside = true;
            outward_normal = hit_record.normal;
            ni_over_nt = self.nor_ref_idx;
            ni_over_nt_seq = self.nor_ref_idx_seq;
            -cosine
        } else {
            outward_normal = vec3_mul_b(&hit_record.normal, -1.0);
            ni_over_nt = self.ref_idx;// / 1.0;
            ni_over_nt_seq = self.ref_idx_seq;
            self.ref_idx * cosine
        };

        let mut refracted_root: bool = false;
        let mut inside_to_inside: bool = false;
        let scatterd: Ray = match refract(&r_in.direction, &outward_normal, ni_over_nt, ni_over_nt_seq) {
            Some(refracted) => {
                reflect_prob = schlick(cosine, self.schlick_r0, self.schlick_r1); // for real glass maty
                let mut rng = rand::thread_rng();
                let rand: f64 = rng.gen();
                if rand < reflect_prob {
                    Ray {
                        origin: hit_record.p,
                        direction: reflect(&r_in.direction, &outward_normal),
                    }
                } else {
                    refracted_root = true;
                    Ray { origin: hit_record.p, direction: refracted }
                }
            }
            None => {
                inside_to_inside = true;
                Ray {
                    origin: hit_record.p,
                    direction: reflect(&r_in.direction, &outward_normal),
                }
            }
        };

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

pub struct DiffuseLight {
    texture: TextureHandle,
}

impl DiffuseLight {
    pub fn new(texture: TextureHandle) -> Self {
        DiffuseLight { texture }
    }

    fn scatter(&self, _r: &Ray, _hit_record: &HitRecord, _texture_list: &TextureList) -> Option<MatRecord> {
        None
    }

    fn emitted(&self, r: &Ray, hit_record: &HitRecord, texture_list: &TextureList) -> Vector3<f64> {
        if vec3_dot(&hit_record.normal, &r.direction).is_sign_negative() {
            texture_list.get_value(
                hit_record.uv,
                &hit_record.p,
                &self.texture
                )
        } else {
            [0.0, 0.0, 0.0]
        }
    }
}
