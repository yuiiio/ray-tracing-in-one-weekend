use crate::ray::{Ray};
use crate::vec3::{Vector3};
use crate::material::{Material};

pub struct HitRecord<'a> {
    t: f64,
    p: Vector3<f64>,
    normal: Vector3<f64>,
    mat_ptr: &'a Box<dyn Material>
}

impl<'a> HitRecord<'a> {
    pub fn new(t: f64, p: Vector3<f64>, normal: Vector3<f64>, mat_ptr: &'a Box<dyn Material>) -> HitRecord {
        HitRecord {t, p, normal, mat_ptr}
    }

    pub fn get_mat_ptr(&self) -> &Box<dyn Material> {
        self.mat_ptr
    }

    pub fn get_t(&self) -> f64 {
        self.t
    }

    pub fn get_p(&self) -> Vector3<f64> {
        self.p
    }

    pub fn get_normal(&self) -> Vector3<f64> {
        self.normal
    }
}

pub trait Hitable {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
}
