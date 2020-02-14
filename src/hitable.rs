use crate::ray::{Ray};
use crate::vec3::{Vector3};
use crate::material::{Material};

pub struct HitRecord<'a> {
    pub t: f64,
    pub p: Vector3<f64>,
    pub normal: Vector3<f64>,
    pub mat_ptr: &'a Box<dyn Material>
}

pub trait Hitable {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
}
