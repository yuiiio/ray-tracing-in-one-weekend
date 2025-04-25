use crate::ray::Ray;
use crate::utils::{max, min};
use crate::vec3::Vector3;
//use core::array;
use core::mem::swap;
use core::simd::Simd;
//use core::simd::cmp::SimdPartialOrd;
use core::simd::num::SimdFloat;

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
        r_origin: &Vector3<f64>,
        r_dir_div: &Vector3<f64>,
        t_min: f64,
        t_max: f64,
    ) -> Option<AabbHitRecord> {
        let mut tmin = t_min;
        let mut tmax = t_max;
        for i in 0..3 {
            let mut t0 = (self.b_min[i] - r_origin[i]) * r_dir_div[i];
            let mut t1 = (self.b_max[i] - r_origin[i]) * r_dir_div[i];

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

// BvhNode owned ?
#[derive(Clone)]
pub struct QBoxes {
    // for simd iteration
    pub bboxes: [[[f64; 4]; 3]; 2], /* quad-boxes; x,y,z; min,max*/
}

// expect compiler simd ?
pub fn aabb_hit_simd(
    /*
    aabb_one: &Aabb,
    aabb_two: &Aabb,
    aabb_three: &Aabb,
    aabb_four: &Aabb,
    */
    qbox: &QBoxes,
    r_origin: &Vector3<f64>,
    r_dir_div: &Vector3<f64>,
    t_min: f64,
    t_max: f64,
) -> (bool, bool, bool, bool) {
    /*
    let mut one_tmin = t_min;
    let mut one_tmax = t_max;
    let mut two_tmin = t_min;
    let mut two_tmax = t_max;
    let mut three_tmin = t_min;
    let mut three_tmax = t_max;
    let mut four_tmin = t_min;
    let mut four_tmax = t_max;
    */
    let mut tmin = Simd::from([t_min; 4]);
    let mut tmax = Simd::from([t_max; 4]);
    for i in 0..3 {
        let mut t0 = (Simd::from(qbox.bboxes[0][i]) - Simd::from([r_origin[i]; 4]))
            * Simd::from([r_dir_div[i]; 4]);
        let mut t1 = (Simd::from(qbox.bboxes[1][i]) - Simd::from([r_origin[i]; 4]))
            * Simd::from([r_dir_div[i]; 4]);

        /*
        let mut one_t0 = (aabb_one.b_min[i] - r_origin[i]) * r_dir_div[i];
        let mut one_t1 = (aabb_one.b_max[i] - r_origin[i]) * r_dir_div[i];
        let mut two_t0 = (aabb_two.b_min[i] - r_origin[i]) * r_dir_div[i];
        let mut two_t1 = (aabb_two.b_max[i] - r_origin[i]) * r_dir_div[i];
        let mut three_t0 = (aabb_three.b_min[i] - r_origin[i]) * r_dir_div[i];
        let mut three_t1 = (aabb_three.b_max[i] - r_origin[i]) * r_dir_div[i];
        let mut four_t0 = (aabb_four.b_min[i] - r_origin[i]) * r_dir_div[i];
        let mut four_t1 = (aabb_four.b_max[i] - r_origin[i]) * r_dir_div[i];
        */

        if r_dir_div[i].is_sign_negative() {
            swap(&mut t0, &mut t1);
            /*
            swap(&mut one_t0, &mut one_t1);
            swap(&mut two_t0, &mut two_t1);
            swap(&mut three_t0, &mut three_t1);
            swap(&mut four_t0, &mut four_t1);
            */
        }

        tmin = tmin.simd_max(t0);
        tmax = tmax.simd_min(t1);
        /*
        one_tmin = max(one_t0, one_tmin);
        one_tmax = min(one_t1, one_tmax);
        two_tmin = max(two_t0, two_tmin);
        two_tmax = min(two_t1, two_tmax);
        three_tmin = max(three_t0, three_tmin);
        three_tmax = min(three_t1, three_tmax);
        four_tmin = max(four_t0, four_tmin);
        four_tmax = min(four_t1, four_tmax);
        */
    }

    let one_result = if tmax[0] < tmin[0] { false } else { true };
    let two_result = if tmax[1] < tmin[1] { false } else { true };
    let three_result = if tmax[2] < tmin[2] { false } else { true };
    let four_result = if tmax[3] < tmin[3] { false } else { true };
    (one_result, two_result, three_result, four_result)
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
        let result = match aabb_box.aabb_hit(&r.origin, &r.get_inv_dir(), 0.00001, 10000.0) {
            Some(_hitrec) => true,
            None => false,
        };
        assert_eq!(true, result);
        let r = Ray {
            origin: [0.0, 0.0, 0.0],
            direction: [1.5, 0.0, 1.5],
        };
        let result = match aabb_box.aabb_hit(&r.origin, &r.get_inv_dir(), 0.00001, 10000.0) {
            Some(_hitrec) => true,
            None => false,
        };
        assert_eq!(false, result);
        let r = Ray {
            origin: [3.0, 3.0, 3.0],
            direction: [-1.0, -1.0, -1.0],
        };
        let result = match aabb_box.aabb_hit(&r.origin, &r.get_inv_dir(), 0.00001, 10000.0) {
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
