use crate::utils::{max, min};
use crate::vec3::{vec3_mul_b, Vector3};
//use core::array;
use core::mem::swap;
use core::simd::cmp::SimdPartialOrd;
use core::simd::num::SimdFloat;
use core::simd::{f64x4, Simd};

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
    pub fn scale(&mut self, scale_value: f64) {
        self.b_min = vec3_mul_b(&self.b_min, scale_value);
        self.b_max = vec3_mul_b(&self.b_max, scale_value);
    }
}

#[derive(Clone)]
pub struct QBoxes {
    // for simd iteration
    pub bboxes: [[[f64; 4]; 3]; 2], /* quad-boxes; x,y,z; min,max*/
}

impl QBoxes {
    pub fn encode_from_four_aabb(aabb_list: &[&Aabb; 4]) -> Self {
        QBoxes {
            bboxes: [
                [
                    [
                        aabb_list[0].b_min[0],
                        aabb_list[1].b_min[0],
                        aabb_list[2].b_min[0],
                        aabb_list[3].b_min[0],
                    ],
                    [
                        aabb_list[0].b_min[1],
                        aabb_list[1].b_min[1],
                        aabb_list[2].b_min[1],
                        aabb_list[3].b_min[1],
                    ],
                    [
                        aabb_list[0].b_min[2],
                        aabb_list[1].b_min[2],
                        aabb_list[2].b_min[2],
                        aabb_list[3].b_min[2],
                    ],
                ],
                [
                    [
                        aabb_list[0].b_max[0],
                        aabb_list[1].b_max[0],
                        aabb_list[2].b_max[0],
                        aabb_list[3].b_max[0],
                    ],
                    [
                        aabb_list[0].b_max[1],
                        aabb_list[1].b_max[1],
                        aabb_list[2].b_max[1],
                        aabb_list[3].b_max[1],
                    ],
                    [
                        aabb_list[0].b_max[2],
                        aabb_list[1].b_max[2],
                        aabb_list[2].b_max[2],
                        aabb_list[3].b_max[2],
                    ],
                ],
            ],
        }
    }
    pub fn scale(&mut self, scale_value: f64) {
        self.bboxes[0][0] = (Simd::from(self.bboxes[0][0]) * f64x4::splat(scale_value)).to_array();
        self.bboxes[0][1] = (Simd::from(self.bboxes[0][1]) * f64x4::splat(scale_value)).to_array();
        self.bboxes[0][2] = (Simd::from(self.bboxes[0][2]) * f64x4::splat(scale_value)).to_array();

        self.bboxes[1][0] = (Simd::from(self.bboxes[1][0]) * f64x4::splat(scale_value)).to_array();
        self.bboxes[1][1] = (Simd::from(self.bboxes[1][1]) * f64x4::splat(scale_value)).to_array();
        self.bboxes[1][2] = (Simd::from(self.bboxes[1][2]) * f64x4::splat(scale_value)).to_array();
    }
}

pub fn aabb_hit_simd(
    qbox: &QBoxes,
    r_origin: &Vector3<f64>,
    r_dir_div: &Vector3<f64>,
    t_min: f64,
    t_max: f64,
) -> u64 {
    let mut tmin = f64x4::splat(t_min);
    let mut tmax = f64x4::splat(t_max);
    for i in 0..3 {
        let t0 = (Simd::from(qbox.bboxes[0][i]) - f64x4::splat(r_origin[i]))
            * f64x4::splat(r_dir_div[i]);
        let t1 = (Simd::from(qbox.bboxes[1][i]) - f64x4::splat(r_origin[i]))
            * f64x4::splat(r_dir_div[i]);

        let mut t0 = &t0;
        let mut t1 = &t1;

        if r_dir_div[i].is_sign_negative() {
            swap(&mut t0, &mut t1);
        }

        tmin = tmin.simd_max(*t0);
        tmax = tmax.simd_min(*t1);
    }

    // save result as Mask<i64, 4> for each sign
    tmax.simd_ge(tmin).to_bitmask()
}

mod test {
    #![allow(unused_imports)]
    use super::*;
    use crate::ray::Ray;

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
