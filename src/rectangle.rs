
use crate::hitable::{HitRecord, Hitable};
use crate::ray::{Ray};
use crate::vec3::{Vector3, vec3_sub, vec3_dot, vec3_div_b, vec3_add_b, vec3_sub_b};
use std::f64::consts::PI;
use crate::material::{MaterialHandle};
use crate::aabb::{Aabb};

#[derive(Clone)]
pub enum AxisType {
    kXY,
    kXZ,
    kYZ,
}

#[derive(Clone)]
pub struct Rect {
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
    k: f64,
    axis: AxisType,
    mat_ptr: MaterialHandle,
}

impl Rect {
    pub fn new(x0: f64, x1: f64, y0: f64, y1: f64, k: f64, axis: AxisType, mat_ptr: MaterialHandle) -> Self {
        Rect {x0, x1, y0, y1, k, axis, mat_ptr}
    }
}

impl Hitable for Rect {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {

        let (xi, yi, zi, nnormal): (usize, usize, usize, Vector3<f64>) = match self.axis {
            AxisType::kXY => (0, 1, 2, [0.0, 0.0, 1.0]),
            AxisType::kXZ => (0, 2, 1, [0.0, 1.0, 0.0]),
            AxisType::kYZ => (1, 2, 0, [1.0, 0.0, 0.0]),
        };

        let t = (self.k - r.origin()[zi]) / r.direction()[zi];
        if t < t_min || t > t_max {
            return None
        }
        let x = r.origin()[xi] + (r.direction()[xi] * t);
        let y = r.origin()[yi] + (r.direction()[yi] * t);
        if x < self.x0 || x > self.x1 || y < self.y0 || y > self.y1 {
            return None
        }
        let u = x - self.x0 / (self.x1 - self.x0);
        let v = y - self.y0 / (self.y1 - self.y0);
        let p = r.point_at_parameter(t);
        Some(HitRecord::new(t, u, v, p, nnormal, MaterialHandle(self.mat_ptr.0)))
    }

    fn bounding_box(&self) -> Option<Aabb> {
        Some(Aabb::new([self.x0, self.y0, self.k-0.0001], [self.x1, self.y0, self.k+0.0001]))
    }
}