use rand::prelude::*;

use crate::aabb::Aabb;
use crate::hitable::{HitRecord, Hitable};
use crate::hitablelist::HitableList;
use crate::material::MaterialHandle;
use crate::ray::Ray;
use crate::vec3::{vec3_dot, vec3_length_f64, vec3_mul_b, vec3_squared_length, vec3_sub, Vector3};

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
    area: f64,
    aabb_box: Aabb,
}

impl Rect {
    pub fn new(
        x0: f64,
        x1: f64,
        y0: f64,
        y1: f64,
        k: f64,
        axis: AxisType,
        mat_ptr: MaterialHandle,
    ) -> Self {
        let area: f64 = (x1 - x0) * (y1 - y0);
        let aabb_box = match axis {
            AxisType::kXY => Aabb::new(
                [x0, y0, k - 0.0001],
                [x1, y1, k + 0.0001],
                ),
            AxisType::kXZ => Aabb::new(
                [x0, k - 0.0001, y0],
                [x1, k + 0.0001, y1],
                ),
            AxisType::kYZ => Aabb::new(
                [k - 0.0001, x0, y0],
                [k + 0.0001, x1, y1],
                ),
        };
        Rect {
            x0,
            x1,
            y0,
            y1,
            k,
            axis,
            mat_ptr,
            area,
            aabb_box,
        }
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
            return None;
        }
        let x = r.origin()[xi] + (r.direction()[xi] * t);
        let y = r.origin()[yi] + (r.direction()[yi] * t);
        if x < self.x0 || x > self.x1 || y < self.y0 || y > self.y1 {
            return None;
        }
        let u = (x - self.x0) / (self.x1 - self.x0);
        let v = (y - self.y0) / (self.y1 - self.y0);
        let p = r.point_at_parameter(t);
        Some(HitRecord::new(
            t,
            u,
            v,
            p,
            nnormal,
            self.mat_ptr,
        ))
    }

    fn bounding_box<'a>(&'a self) -> Option<&'a Aabb> {
        Some(&self.aabb_box)
    }

    fn pdf_value(&self, o: &Vector3<f64>, v: &Vector3<f64>) -> f64 {
        match self.hit(&Ray::new(*o, *v), 0.00001, 10000.0) {
            Some(rec) => {
                let distance_squared = rec.get_t().powi(2) * vec3_squared_length(v);
                let cosine = vec3_dot(v, &rec.get_normal()).abs() / vec3_length_f64(v);
                return distance_squared / (cosine * self.area);
            }
            None => return 0.0,
        }
    }

    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        let mut rng = rand::thread_rng();
        let rng_x: f64 = rng.gen();
        let rng_y: f64 = rng.gen();
        let random_point = match self.axis {
            AxisType::kXY => [
                self.x0 + rng_x * (self.x1 - self.x0),
                self.y0 + rng_y * (self.y1 - self.y0),
                self.k,
            ],
            AxisType::kXZ => [
                self.x0 + rng_x * (self.x1 - self.x0),
                self.k,
                self.y0 + rng_y * (self.y1 - self.y0),
            ],
            AxisType::kYZ => [
                self.k,
                self.x0 + rng_x * (self.x1 - self.x0),
                self.y0 + rng_y * (self.y1 - self.y0),
            ],
        };
        vec3_sub(&random_point, o)
    }
}

#[derive(Clone)]
pub struct FlipNormals {
    shape: Rect,
}

impl FlipNormals {
    pub fn new(shape: Rect) -> Self {
        FlipNormals { shape }
    }
}

impl Hitable for FlipNormals {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        match self.shape.hit(r, t_min, t_max) {
            Some(hit) => Some(HitRecord::new(
                hit.get_t(),
                hit.get_u(),
                hit.get_v(),
                hit.get_p(),
                vec3_mul_b(&hit.get_normal(), -1.0),
                hit.get_mat_ptr(),
            )),
            None => None,
        }
    }

    fn bounding_box<'a>(&'a self) -> Option<&'a Aabb> {
        Some(&self.shape.bounding_box().unwrap())
    }

    fn pdf_value(&self, o: &Vector3<f64>, v: &Vector3<f64>) -> f64 {
        self.shape.pdf_value(o, v)
    }

    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        self.shape.random(o)
    }
}

#[derive(Clone)]
pub struct Boxel {
    pmin: Vector3<f64>,
    pmax: Vector3<f64>,
    list: HitableList,
    aabb_box: Aabb,
}

impl Boxel {
    pub fn new(p0: Vector3<f64>, p1: Vector3<f64>, mat_ptr: MaterialHandle) -> Self {
        let pmin = p0;
        let pmax = p1;
        let mut list = HitableList::new();
        list.push(Rect::new(
            p0[0],
            p1[0],
            p0[1],
            p1[1],
            p1[2],
            AxisType::kXY,
            mat_ptr,
        ));
        list.push(FlipNormals::new(Rect::new(
            p0[0],
            p1[0],
            p0[1],
            p1[1],
            p0[2],
            AxisType::kXY,
            mat_ptr,
        )));

        list.push(Rect::new(
            p0[0],
            p1[0],
            p0[2],
            p1[2],
            p1[1],
            AxisType::kXZ,
            mat_ptr,
        ));
        list.push(FlipNormals::new(Rect::new(
            p0[0],
            p1[0],
            p0[2],
            p1[2],
            p0[1],
            AxisType::kXZ,
            mat_ptr,
        )));

        list.push(Rect::new(
            p0[1],
            p1[1],
            p0[2],
            p1[2],
            p1[0],
            AxisType::kYZ,
            mat_ptr,
        ));
        list.push(FlipNormals::new(Rect::new(
            p0[1],
            p1[1],
            p0[2],
            p1[2],
            p0[0],
            AxisType::kYZ,
            mat_ptr,
        )));
        let aabb_box = Aabb::new(pmin, pmax);
        Boxel { pmin, pmax, list, aabb_box}
    }
}

impl Hitable for Boxel {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        self.list.hit(r, t_min, t_max)
    }

    fn bounding_box<'a>(&'a self) -> Option<&'a Aabb> {
        Some(&self.aabb_box)
    }

    fn pdf_value(&self, o: &Vector3<f64>, v: &Vector3<f64>) -> f64 {
        self.list.pdf_value(o, v)
    }

    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        self.list.random(o)
    }
}
