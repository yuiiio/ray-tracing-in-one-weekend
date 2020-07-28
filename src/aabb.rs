use crate::vec3::{Vector3};
use crate::ray::{Ray};
use crate::utils::{min, max};
use crate::hitable::{Hitable, HitRecord};
use crate::material::{MaterialHandle};
use std::mem::swap;

#[derive(Clone)]
pub struct Aabb {
    min: Vector3<f64>,
    max: Vector3<f64>,
}

impl Aabb {
    pub fn new(min: Vector3<f64>, max: Vector3<f64>) -> Aabb {
        Aabb { min, max }
    }

    pub fn min(&self) -> Vector3<f64> {
        self.min
    }

    pub fn max(&self) -> Vector3<f64> {
        self.max
    }
}

impl Hitable for Aabb {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut tmin = t_min;
        let mut tmax = t_max;
        for i in 0..3 {
            let inv_d = 1.0 / r.direction()[i];
            let mut t0 = (self.min[i] - r.origin()[i]) / inv_d;
            let mut t1 = (self.max[i] - r.origin()[i]) / inv_d;

            if (inv_d < 0.0) {
                swap(&mut t0, &mut t1);
            }
        
            tmin = max(t0, tmin);
            tmax = min(t1, tmax);
            if tmax < tmin {
                return None
            }
        }
        let t = 0.0;
        let u = 0.0;
        let v = 0.0;
        let p = [0.0, 0.0, 0.0];
        let normal = [0.0, 0.0, 0.0];
        let mat_ptr = MaterialHandle(0);
        return Some(HitRecord::new(t, u, v, p, normal, mat_ptr))
    }

    fn bounding_box(&self) -> Option<Aabb> {
        Some( self.clone() )
    }
}

mod test {
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
        let r = Ray::new([0.0, 0.0, 0.0], [1.5, 1.5, 1.5]);
        let result = match aabb_box.hit(&r, 0.00001, 10000.0) {
            Some(_hitrec) => true,
            None => false,
         };
        assert_eq!(true, result);
        let r = Ray::new([0.0, 0.0, 0.0], [1.5, 0.0, 1.5]);
        let result = match aabb_box.hit(&r, 0.00001, 10000.0) {
            Some(_hitrec) => true,
            None => false,
         };
        assert_eq!(false, result);
        let r = Ray::new([3.0, 3.0, 3.0], [-1.0, -1.0, -1.0]);
        let result = match aabb_box.hit(&r, 0.00001, 10000.0) {
            Some(_hitrec) => true,
            None => false,
         };
        assert_eq!(true, result);
    }
}

pub fn surrounding_box(box0: Aabb, box1: Aabb) -> Aabb {
    let min = [min(box0.min()[0], box1.min()[0]),
                min(box0.min()[1], box1.min()[1]),
                min(box0.min()[2], box1.min()[2])];
    let max = [max(box0.max()[0], box1.max()[0]),
                max(box0.max()[1], box1.max()[1]),
                max(box0.max()[2], box1.max()[2])];
    Aabb::new(min, max)
}