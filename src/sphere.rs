use crate::hitable::{HitRecord, Hitable};
use crate::ray::{Ray};
use crate::vec3::{Vector3, vec3_sub, vec3_dot, vec3_div_b};
use std::f64;
use crate::material::{MaterialHandle};

pub struct Sphere {
    center: Vector3<f64>,
    radius: f64,
    mat_ptr: MaterialHandle,
}

impl Sphere {
    pub fn new(center: Vector3<f64>, radius: f64, mat_ptr: MaterialHandle) -> Sphere {
        Sphere {center, radius, mat_ptr}
    }
}

impl Hitable for Sphere {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let rec: Option<HitRecord> = None;
        let oc = vec3_sub(r.origin(), self.center);
        let a = vec3_dot(r.direction(), r.direction());
        let b = 2.0 * vec3_dot(r.direction(), oc);
        let c = vec3_dot(oc, oc) - self.radius * self.radius;
        let descriminant = b * b - 4.0 * a * c;
        if descriminant >= 0.0 {
            let temp = (-b - descriminant.sqrt()) / (2.0 * a);
            if  temp < t_max && temp > t_min {
                let point = r.point_at_parameter(temp);
                let nnormal = vec3_div_b(vec3_sub(point, self.center), self.radius);
                return Some(HitRecord::new(temp, point, nnormal, MaterialHandle(self.mat_ptr.0)));
            }
            let temp = (-b + descriminant.sqrt()) / (2.0 * a);
            if  temp < t_max && temp > t_min {
                let point = r.point_at_parameter(temp);
                let nnormal = vec3_div_b(vec3_sub(point, self.center), self.radius);
                return Some(HitRecord::new(temp, point, nnormal, MaterialHandle(self.mat_ptr.0)));
            }
        }
        rec
    }
} 
