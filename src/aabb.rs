use crate::vec3::Vector3;
use crate::ray::Ray;
use crate::utils::{min, max};
use std::mem::swap;

pub struct AabbHitRecord {
}

#[derive(Clone)]
pub struct Aabb {
    b_min: Vector3<f64>,
    b_max: Vector3<f64>,
}

impl Aabb {
    pub fn new(b_min: Vector3<f64>, b_max: Vector3<f64>) -> Self {
        Aabb { b_min, b_max }
    }

    pub fn b_min(&self) -> Vector3<f64> {
        self.b_min
    }

    pub fn b_max(&self) -> Vector3<f64> {
        self.b_max
    }

    pub fn aabb_hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<AabbHitRecord> {
        let mut tmin = t_min;
        let mut tmax = t_max;
        for i in 0..3 {
            let inv_d = 1.0 / r.direction[i];
            let mut t0 = (self.b_min[i] - r.origin[i]) * inv_d;
            let mut t1 = (self.b_max[i] - r.origin[i]) * inv_d;

            if inv_d.is_sign_negative() {
                swap(&mut t0, &mut t1);
            }

            tmin = max(t0, tmin);
            tmax = min(t1, tmax);
            if tmax < tmin {
                return None
            }
        }
        return Some(AabbHitRecord{})
    }
}

mod test {
    #![allow(unused_imports)]
    use super::*;

    #[test]
    fn roop_let_test() {
        let mut s = 0;
        for i in 0..5 {
            s = s + i;
        }
        assert_eq!(s, 10);
    }

    #[test]
    fn aabb_hit_test() {
        let aabb_box = Aabb::new([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);
        let r = Ray{ origin: [0.0, 0.0, 0.0], direction: [1.5, 1.5, 1.5] };
        let result = match aabb_box.aabb_hit(&r, 0.00001, 10000.0) {
            Some(_hitrec) => true,
            None => false,
         };
        assert_eq!(true, result);
        let r = Ray{ origin: [0.0, 0.0, 0.0], direction: [1.5, 0.0, 1.5] };
        let result = match aabb_box.aabb_hit(&r, 0.00001, 10000.0) {
            Some(_hitrec) => true,
            None => false,
         };
        assert_eq!(false, result);
        let r = Ray{ origin: [3.0, 3.0, 3.0], direction: [-1.0, -1.0, -1.0] };
        let result = match aabb_box.aabb_hit(&r, 0.00001, 10000.0) {
            Some(_hitrec) => true,
            None => false,
         };
        assert_eq!(true, result);
    }
}

pub fn surrounding_box(box0: &Aabb, box1: &Aabb) -> Aabb {
    let min = [min(box0.b_min()[0], box1.b_min()[0]),
                min(box0.b_min()[1], box1.b_min()[1]),
                min(box0.b_min()[2], box1.b_min()[2])];
    let max = [max(box0.b_max()[0], box1.b_max()[0]),
                max(box0.b_max()[1], box1.b_max()[1]),
                max(box0.b_max()[2], box1.b_max()[2])];
    Aabb::new(min, max)
}
