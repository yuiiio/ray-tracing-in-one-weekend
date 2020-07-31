use crate::hitable::{HitRecord, Hitable};
use crate::ray::{Ray};
use crate::vec3::{Vector3, vec3_sub, vec3_add};
use crate::aabb::{Aabb};

#[derive (Clone)]
pub struct Translate {
    obj: Box<dyn Hitable + Send + Sync>,
    offset: Vector3<f64>,
}

impl Translate {
    pub fn new(obj: Box<dyn Hitable + Send + Sync>, offset: Vector3<f64>) -> Self {
        Translate {obj, offset}
    }
}

impl Hitable for Translate {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let r = Ray::new(vec3_sub(r.origin(), self.offset), r.direction());
        match self.obj.hit(&r, t_min, t_max) {
            Some(hit) => Some(HitRecord::new(hit.get_t(), hit.get_u(), hit.get_v(), vec3_add(hit.get_p(), self.offset), hit.get_normal(), hit.get_mat_ptr())),
            None => None,
        }
    }

    fn bounding_box(&self) -> Option<Aabb> {
        match self.obj.bounding_box() {
            Some(aabb) => Some(Aabb::new(vec3_add(aabb.min(), self.offset), vec3_add(aabb.max(), self.offset))),
            None => None,
        }
    }
}