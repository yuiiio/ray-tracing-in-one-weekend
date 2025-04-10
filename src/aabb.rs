use crate::ray::Ray;
use crate::utils::{max, min};
use crate::vec3::Vector3;
use std::mem::swap;

pub struct AabbHitRecord {
    pub t_max: f64,
    pub t_min: f64,
}

#[derive(Clone)]
pub struct Aabb {
    pub b_min: Vector3<f64>,
    pub b_max: Vector3<f64>,
}

impl Aabb {
    pub fn aabb_hit(
        &self,
        r: &Ray,
        r_dir_div: &Vector3<f64>,
        t_min: f64,
        t_max: f64,
    ) -> Option<AabbHitRecord> {
        let mut tmin = t_min;
        let mut tmax = t_max;
        for i in 0..3 {
            let mut t0 = (self.b_min[i] - r.origin[i]) * r_dir_div[i];
            let mut t1 = (self.b_max[i] - r.origin[i]) * r_dir_div[i];

            if r_dir_div[i].is_sign_negative() {
                swap(&mut t0, &mut t1);
            }

            tmin = max(t0, tmin);
            tmax = min(t1, tmax);
        }

        if tmax < tmin {
            return None;
        }
        Some(AabbHitRecord {
            t_max: tmax,
            t_min: tmin,
        })
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
        let aabb_box = Aabb {
            b_min: [1.0, 1.0, 1.0],
            b_max: [2.0, 2.0, 2.0],
        };
        let r = Ray {
            origin: [0.0, 0.0, 0.0],
            direction: [1.5, 1.5, 1.5],
        };
        let result = match aabb_box.aabb_hit(&r, &r.get_inv_dir(), 0.00001, 10000.0) {
            Some(_hitrec) => true,
            None => false,
        };
        assert_eq!(true, result);
        let r = Ray {
            origin: [0.0, 0.0, 0.0],
            direction: [1.5, 0.0, 1.5],
        };
        let result = match aabb_box.aabb_hit(&r, &r.get_inv_dir(), 0.00001, 10000.0) {
            Some(_hitrec) => true,
            None => false,
        };
        assert_eq!(false, result);
        let r = Ray {
            origin: [3.0, 3.0, 3.0],
            direction: [-1.0, -1.0, -1.0],
        };
        let result = match aabb_box.aabb_hit(&r, &r.get_inv_dir(), 0.00001, 10000.0) {
            Some(_hitrec) => true,
            None => false,
        };
        assert_eq!(true, result);
    }
}

pub fn surrounding_box(box0: &Aabb, box1: &Aabb) -> Aabb {
    let b_min = [
        min(box0.b_min[0], box1.b_min[0]),
        min(box0.b_min[1], box1.b_min[1]),
        min(box0.b_min[2], box1.b_min[2]),
    ];
    let b_max = [
        max(box0.b_max[0], box1.b_max[0]),
        max(box0.b_max[1], box1.b_max[1]),
        max(box0.b_max[2], box1.b_max[2]),
    ];
    Aabb { b_min, b_max }
}
