use crate::hitable::{HitRecord, Hitable};
use crate::ray::{Ray};
use crate::vec3::{Vector3, vec3_sub, vec3_add};
use crate::aabb::{Aabb};
use crate::quotation::{Rotation};

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
            Some(hit) => Some(HitRecord::new(hit.get_t(), hit.get_u(), hit.get_v(),
                            vec3_add(hit.get_p(), self.offset), hit.get_normal(), hit.get_mat_ptr())),
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

#[derive(Clone)]
pub struct Rotate {
    obj: Box<dyn Hitable + Send + Sync>,
    quat: Rotation,
    revq: Rotation,
    aabb: Aabb,
}

impl Rotate {
    pub fn new(obj: Box<dyn Hitable + Send + Sync>, axis: Vector3<f64>, angle: f64) -> Self {
        let quat = Rotation::new(angle, axis);
        let revq = Rotation::new(-angle, axis);

        // let found boundingbox to enough include rotate-obj
        let bbox = obj.bounding_box().unwrap();
        let mut min = [std::f64::MAX; 3];
        let mut max = [-std::f64::MAX; 3];
        for i in [1, 0].iter() {
            for j in [1, 0].iter() {
                for k in [1, 0].iter() {
                    let x = *i as f64 * bbox.max()[0] + (1 - *i) as f64 * bbox.min()[0];
                    let y = *j as f64 * bbox.max()[1] + (1 - *j) as f64 * bbox.min()[1];
                    let z = *k as f64 * bbox.max()[2] + (1 - *k) as f64 * bbox.min()[2];

                    let rotated = quat.rotate([x, y, z]);
                    for c in 0..3 {
                        if rotated[c] > max[c] {
                            max[c] = rotated[c];
                        }
                        if rotated[c] < min[c] {
                            min[c] = rotated[c];
                        }
                    }
                }
            }
        }
        let aabb = Aabb::new(min, max);

        Rotate {obj, quat, revq, aabb}
    }
}

impl Hitable for Rotate {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let origin = self.revq.rotate(r.origin());
        let direction = self.revq.rotate(r.direction());
        let r = Ray::new(origin, direction);
        match self.obj.hit(&r, t_min, t_max) {
            Some(hit) => Some(HitRecord::new(hit.get_t(), hit.get_u(), hit.get_v(),
                            self.quat.rotate(hit.get_p()), self.quat.rotate(hit.get_normal()), hit.get_mat_ptr())),
            None => None,
        }
    }

    fn bounding_box(&self) -> Option<Aabb> {
        Some(self.aabb.clone())
    }
}