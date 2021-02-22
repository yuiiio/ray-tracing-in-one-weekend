use rand::prelude::*;

use crate::aabb::Aabb;
use crate::hitable::{HitRecord, Hitable};
use crate::hitablelist::HitableList;
use crate::material::MaterialHandle;
use crate::ray::Ray;
use crate::utils::{max, min};
use crate::vec3::{
    cross, vec3_add, vec3_dot, vec3_length_f64, vec3_mul_b, vec3_squared_length, vec3_sub,
    vec3_unit_vector_f64, Vector3,
};

#[derive(Clone)]
pub struct Triangle {
    v0: Vector3<f64>,
    v1: Vector3<f64>,
    v2: Vector3<f64>,
    mat_ptr: MaterialHandle,
    e1: Vector3<f64>,
    e2: Vector3<f64>,
    n: Vector3<f64>,
    area: f64,
}

impl Triangle {
    pub fn new(
        v0: Vector3<f64>,
        v1: Vector3<f64>,
        v2: Vector3<f64>,
        mat_ptr: MaterialHandle,
    ) -> Self {
        let e1 = vec3_sub(v1, v0);
        let e2 = vec3_sub(v2, v0);
        let n = cross(e2, e1);
        let area: f64 = vec3_length_f64(n) / 2.0;
        Triangle {
            v0,
            v1,
            v2,
            mat_ptr,
            e1,
            e2,
            n,
            area,
        }
    }
}

impl Hitable for Triangle {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let d = ray.direction();
        let lal = -1.0 * vec3_dot(d, self.n);
        if lal == 0.0 {
            return None;
        }
        let nor_lal: f64 = 1.0 / lal;

        let r = vec3_sub(ray.origin(), self.v0);
        let m = cross(d, r);

        let u = nor_lal * vec3_dot(self.e2, m);
        if u < 0.0 || u > 1.0 {
            return None;
        }

        let v = nor_lal * -1.0 * vec3_dot(self.e1, m);
        if v < 0.0 || v > 1.0 {
            return None;
        }

        if u + v > 1.0 {
            return None;
        }

        let t = nor_lal * vec3_dot(r, self.n);
        if t < t_min || t > t_max {
            return None;
        }

        let normal = vec3_unit_vector_f64(self.n);
        let p: Vector3<f64> = ray.point_at_parameter(t);
        let uu = u;
        let vv = v;
        return Some(HitRecord::new(t, uu, vv, p, normal, self.mat_ptr));
    }

    fn bounding_box(&self) -> Option<Aabb> {
        let min = [
            min(min(self.v0[0], self.v1[0]), self.v2[0]),
            min(min(self.v0[1], self.v1[1]), self.v2[1]),
            min(min(self.v0[2], self.v1[2]), self.v2[2]),
        ];
        let max = [
            max(max(self.v0[0], self.v1[0]), self.v2[0]),
            max(max(self.v0[1], self.v1[1]), self.v2[1]),
            max(max(self.v0[2], self.v1[2]), self.v2[2]),
        ];
        return Some(Aabb::new(min, max));
    }

    fn pdf_value(&self, o: &Vector3<f64>, v: &Vector3<f64>) -> f64 {
        match self.hit(&Ray::new(*o, *v), 0.00001, 10000.0) {
            Some(rec) => {
                let distance_squared = rec.get_t().powi(2) * vec3_squared_length(*v);
                let cosine = vec3_dot(*v, rec.get_normal()).abs() / vec3_length_f64(*v);
                return distance_squared / (cosine * self.area);
            }
            None => return 0.0,
        }
    }

    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        let mut rng = rand::thread_rng();
        let rng_i: f64 = rng.gen();
        let rng_j: f64 = rng.gen();

        let rng_i = rng_i;
        let rng_j = rng_j;

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
            vec3_add(vec3_mul_b(self.v0, u), vec3_mul_b(self.v1, v)),
            vec3_mul_b(self.v2, w),
        );

        vec3_sub(random_point, *o)
    }
}
