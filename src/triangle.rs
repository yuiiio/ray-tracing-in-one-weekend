use rand::prelude::*;

use crate::aabb::Aabb;
use crate::hitable::{HitRecord, Hitable};
use crate::hitablelist::HitableList;
use crate::material::MaterialHandle;
use crate::ray::Ray;
use crate::utils::{max, min};
use crate::vec3::{
    cross, vec3_dot, vec3_length_f64, vec3_mul_b, vec3_squared_length, vec3_sub, Vector3,
};

#[derive(Clone)]
pub struct Tri {
    v0: Vector3<f64>,
    v1: Vector3<f64>,
    v2: Vector3<f64>,
    mat_ptr: MaterialHandle,
    e1: Vector3<f64>,
    e2: Vector3<f64>,
    n: Vector3<f64>,
}

impl Tri {
    pub fn new(
        v0: Vector3<f64>,
        v1: Vector3<f64>,
        v2: Vector3<f64>,
        mat_ptr: MaterialHandle,
    ) -> Self {
        let e1 = vec3_sub(v1, v0);
        let e2 = vec3_sub(v2, v0);
        let n = cross(e2, e1);
        Tri {
            v0,
            v1,
            v2,
            mat_ptr,
            e1,
            e2,
            n,
        }
    }
}

impl Hitable for Tri {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let d = ray.direction();
        let r = vec3_sub(ray.origin(), self.v0);
        let nor_lal: f64 = 1.0 / vec3_dot(d, self.n);
        let m = cross(d, r);

        let u = nor_lal * vec3_dot(self.e2, m);
        let v = nor_lal * -1.0 * vec3_dot(self.e1, m);

        if u >= 0.0 && v >= 0.0 && u + v <= 1.0 {
            let t = nor_lal * vec3_dot(r, self.n);
            let normal = self.n;
            let p: Vector3<f64> = ray.point_at_parameter(t);
            let uu = u;
            let vv = v;
            return Some(HitRecord::new(t, uu, vv, p, normal, self.mat_ptr));
        } else {
            return None;
        }
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
}
