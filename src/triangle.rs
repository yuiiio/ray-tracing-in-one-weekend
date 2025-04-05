use rand::prelude::*;

use crate::aabb::Aabb;
use crate::hitable::{HitRecord, Hitable};
use crate::material::MaterialHandle;
use crate::onb::Onb;
use crate::quotation::Rotation;
use crate::ray::Ray;
use crate::utils::{max, min};
use crate::vec3::{
    cross, vec3_add, vec3_dot, vec3_length_f64, vec3_mul_b, vec3_sub, vec3_unit_vector_f64, Vector3,
};

#[derive(Clone)]
pub struct Triangle {
    v0: Vector3<f64>,
    v1: Vector3<f64>,
    v2: Vector3<f64>,
    mat_ptr: MaterialHandle,
    e1: Vector3<f64>,
    e2: Vector3<f64>,
    n: Vector3<f64>, // cross(e2, e1) so size is not normal
    n_norm: Vector3<f64>,
    area: f64,
    aabb_box: Aabb,
    onb_uv: (Vector3<f64>, Vector3<f64>),
}

impl Triangle {
    pub fn new(
        v0: Vector3<f64>,
        v1: Vector3<f64>,
        v2: Vector3<f64>,
        mat_ptr: MaterialHandle,
    ) -> Self {
        let e1 = vec3_sub(&v1, &v0);
        let e2 = vec3_sub(&v2, &v0);
        let n = cross(&e2, &e1);
        let area: f64 = vec3_length_f64(&n) / 2.0;

        let b_min = [
            min(min(v0[0], v1[0]), v2[0]),
            min(min(v0[1], v1[1]), v2[1]),
            min(min(v0[2], v1[2]), v2[2]),
        ];
        let b_max = [
            max(max(v0[0], v1[0]), v2[0]),
            max(max(v0[1], v1[1]), v2[1]),
            max(max(v0[2], v1[2]), v2[2]),
        ];
        let aabb_box = Aabb { b_min, b_max };

        let n_norm = vec3_unit_vector_f64(&n);
        let onb = Onb::build_from_w(&n_norm);
        Triangle {
            v0,
            v1,
            v2,
            mat_ptr,
            e1,
            e2,
            n,
            n_norm,
            area,
            aabb_box,
            onb_uv: (onb.u, onb.v),
        }
    }
}

impl Hitable for Triangle {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let d = ray.direction;
        let lal = vec3_dot(&d, &self.n);
        if lal == 0.0 {
            return None;
        }
        let nor_lal: f64 = 1.0 / lal;

        let r = vec3_sub(&ray.origin, &self.v0);
        let m = cross(&d, &r);

        let v = nor_lal * vec3_dot(&self.e1, &m);
        if v.is_sign_negative() || v > 1.0 {
            return None;
        }

        let nor_lal = -1.0 * nor_lal;

        let u = nor_lal * vec3_dot(&self.e2, &m);
        if u.is_sign_negative() || u > 1.0 {
            return None;
        }

        if u + v > 1.0 {
            return None;
        }

        let t = nor_lal * vec3_dot(&r, &self.n);
        if t < t_min || t > t_max {
            return None;
        }

        let p: Vector3<f64> = ray.point_at_parameter(t);
        Some(HitRecord {
            t,
            uv: (u, v),
            p,
            normal: self.n_norm,
            mat_ptr: &self.mat_ptr,
            onb_uv: Some(&self.onb_uv),
        })
    }

    fn bounding_box(&self) -> &Aabb {
        &self.aabb_box
    }

    fn pdf_value(&self, ray: &Ray) -> f64 {
        if let Some(rec) = self.hit(ray, 0.00001, 10000.0) {
            let distance_squared = rec.t.powi(2);
            let cosine = vec3_dot(&ray.direction, &rec.normal).abs();
            return distance_squared / (cosine * self.area);
        }
        0.0
    }

    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        let mut rng = rand::thread_rng();
        let rng_i: f64 = rng.gen();
        let rng_j: f64 = rng.gen();

        let max: f64;
        let min: f64;
        if rng_i < rng_j {
            max = rng_j;
            min = rng_i;
        } else {
            max = rng_i;
            min = rng_j;
        }

        let u = min;
        let v = 1.0 - max;
        let w = max - min;

        let random_point = vec3_add(
            &vec3_add(&vec3_mul_b(&self.v0, u), &vec3_mul_b(&self.v1, v)),
            &vec3_mul_b(&self.v2, w),
        );

        vec3_unit_vector_f64(&vec3_sub(&random_point, o)) // random should return normalized vec
    }

    fn rotate_onb(&mut self, quat: &Rotation) -> () {
        self.n_norm = quat.rotate(&self.n_norm);
        let onb = Onb::build_from_w(&self.n_norm);
        self.onb_uv = (onb.u, onb.v);
    }
}
