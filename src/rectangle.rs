use rand::prelude::*;

use crate::aabb::Aabb;
use crate::hitable::{HitRecord, Hitable};
use crate::material::MaterialHandle;
use crate::ray::Ray;
use crate::vec3::{vec3_dot, vec3_mul_b, vec3_sub, vec3_unit_vector_f64, Vector3};

#[derive(Clone)]
pub enum AxisType {
    Kxy,
    Kxz,
    Kyz,
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
    needs_uv: bool,
    normal: Vector3<f64>,
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
        flip_normal: bool,
    ) -> Self {
        let mut normal = match axis {
            AxisType::Kxy => 
                [0.0, 0.0, 1.0],
            AxisType::Kxz => 
                [0.0, 1.0, 0.0],
            AxisType::Kyz => 
                [1.0, 0.0, 0.0],
        };
        if flip_normal == true {
            normal = vec3_mul_b(&normal, -1.0);
        }
        let needs_uv = mat_ptr.needs_uv;
        Rect {
            x0,
            x1,
            y0,
            y1,
            k,
            axis,
            mat_ptr,
            needs_uv,
            normal,
        }
    }
}

// Rect hit is valid even ray direction.
// Rect != actual Surface.

impl Hitable for Rect {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let (xi, yi, zi): (usize, usize, usize) = match self.axis {
            AxisType::Kxy => (0, 1, 2),
            AxisType::Kxz => (0, 2, 1),
            AxisType::Kyz => (1, 2, 0),
        };

        let t = (self.k - r.origin[zi]) / r.direction[zi];
        if t < t_min || t > t_max {
            return None;
        }
        let x = r.origin[xi] + (r.direction[xi] * t);
        if x < self.x0 || x > self.x1 {
            return None;
        }
        let y = r.origin[yi] + (r.direction[yi] * t);
        if y < self.y0 || y > self.y1 {
            return None;
        }

        let (u, v) = if self.needs_uv {
            let width = self.x1 - self.x0;
            let height = self.y1 - self.y0;
            (
                (x - self.x0) / width,
                (y - self.y0) / height,
            )
        } else {
            (0.0, 0.0)
        };

        let p = r.point_at_parameter(t);
        Some(HitRecord {
            t,
            uv: (u, v),
            p,
            normal: self.normal,
            mat_ptr: self.mat_ptr.clone(),
        })
    }

    fn bounding_box(&self) -> Aabb {
        match self.axis {
            AxisType::Kxy => 
                Aabb {
                    b_min: [self.x0, self.y0, self.k - 0.0001],
                    b_max: [self.x1, self.y1, self.k + 0.0001],
                },
            AxisType::Kxz => 
                Aabb {
                    b_min: [self.x0, self.k - 0.0001, self.y0],
                    b_max: [self.x1, self.k + 0.0001, self.y1],
                },
            AxisType::Kyz => 
                Aabb {
                    b_min: [self.k - 0.0001, self.x0, self.y0],
                    b_max: [self.k + 0.0001, self.x1, self.y1],
                },
        }
    }

    fn pdf_value(&self, ray: &Ray) -> f64 {
        if let Some(rec) = self.hit(ray, 0.00001, 10000.0) {
            let distance_squared = rec.t.powi(2);
            let cosine = vec3_dot(&ray.direction, &rec.normal).abs();
            let area: f64 = (self.x1 - self.x0) * (self.y1 - self.y0);
            return distance_squared / (cosine * area);
        }
        0.0
    }

    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        let mut rng = rand::rng();
        let rng_x: f64 = rng.random();
        let rng_y: f64 = rng.random();
        let width = self.x1 - self.x0;
        let height = self.y1 - self.y0;
        let random_point = match self.axis {
            AxisType::Kxy => [
                self.x0 + rng_x * width,
                self.y0 + rng_y * height,
                self.k,
            ],
            AxisType::Kxz => [
                self.x0 + rng_x * width,
                self.k,
                self.y0 + rng_y * height,
            ],
            AxisType::Kyz => [
                self.k,
                self.x0 + rng_x * width,
                self.y0 + rng_y * height,
            ],
        };
        //TODO: need all rewrite
        vec3_unit_vector_f64(&vec3_sub(&random_point, o))
    }

    fn scale(&mut self, scale_value: f64) {
        self.x0 *= scale_value;
        self.x1 *= scale_value;
        self.y0 *= scale_value;
        self.y1 *= scale_value;
        self.k *= scale_value;
    }

    fn set_material(&mut self, mat_ptr: MaterialHandle) {
        self.mat_ptr = mat_ptr;
    }
}

#[derive(Clone)]
pub struct Boxel {
    rect: [Rect; 6],
    aabb_box: Aabb,
}

impl Boxel {
    pub fn new(p0: Vector3<f64>, p1: Vector3<f64>, mat_ptr: MaterialHandle) -> Self {
        let b_min = p0;
        let b_max = p1;
        let rect: [Rect; 6] = [
            Rect::new(
                p0[0],
                p1[0],
                p0[1],
                p1[1],
                p1[2],
                AxisType::Kxy,
                mat_ptr.clone(),
                false,
            ),
            Rect::new(
                p0[0],
                p1[0],
                p0[2],
                p1[2],
                p1[1],
                AxisType::Kxz,
                mat_ptr.clone(),
                false,
            ),
            Rect::new(
                p0[1],
                p1[1],
                p0[2],
                p1[2],
                p1[0],
                AxisType::Kyz,
                mat_ptr.clone(),
                false,
            ),
            Rect::new(
                p0[0],
                p1[0],
                p0[1],
                p1[1],
                p0[2],
                AxisType::Kxy,
                mat_ptr.clone(),
                true,
            ),
            Rect::new(
                p0[0],
                p1[0],
                p0[2],
                p1[2],
                p0[1],
                AxisType::Kxz,
                mat_ptr.clone(),
                true,
            ),
            Rect::new(
                p0[1],
                p1[1],
                p0[2],
                p1[2],
                p0[0],
                AxisType::Kyz,
                mat_ptr.clone(),
                true,
            ),
        ];
        let aabb_box = Aabb { b_min, b_max };
        Boxel { rect, aabb_box }
    }
}

impl Hitable for Boxel {
    // needs to check from out-side and inside(eg: glass) rays.
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut rect_index_0 = [5, 4, 3];
        let mut rect_index_1 = [2, 1, 0];

        for i in 0..3 {
            if ray.direction[i].is_sign_negative() {
                core::mem::swap(&mut rect_index_0[i], &mut rect_index_1[i]);
            }
            if let Some(hit_rec) = self.rect[rect_index_0[i]].hit(ray, t_min, t_max) {
                return Some(hit_rec);
            }
        }
        // up is enough when only care out-side rays
        // down is needed for inside rays.
        for i in 0..3 {
            if let Some(hit_rec) = self.rect[rect_index_1[i]].hit(ray, t_min, t_max) {
                return Some(hit_rec);
            }
        }

        None
    }

    fn bounding_box(&self) -> Aabb {
        self.aabb_box.clone()
    }

    // not need checks from inside rays ?
    fn pdf_value(&self, ray: &Ray) -> f64 {
        // TODO: we needs actual pdf hit surface, now return avg all surface
        let r_dir_div = ray.get_inv_dir();
        if let Some(_aabb_hit) = self
            .aabb_box
            .aabb_hit(&ray.origin, &r_dir_div, 0.00001, 10000.0)
        {
            const DIV6: f64 = 1.0 / 6.0;
            let mut pdf_sum = 0.0;
            for i in 0..6 {
                pdf_sum += self.rect[i].pdf_value(ray);
            }
            pdf_sum * DIV6
        } else {
            0.0
        }
    }

    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        let mut rng = rand::rng();
        let rand: f64 = rng.random();
        let random_handle = (rand * 6.0) as usize;
        self.rect[random_handle].random(o)
    }

    fn scale(&mut self, scale_value: f64) {
        for i in 0..6 {
            self.rect[i].scale(scale_value);
        }
        self.aabb_box.scale(scale_value);
    }
    fn set_material(&mut self, mat_ptr: MaterialHandle) {
        for i in 0..6 {
            self.rect[i].set_material(mat_ptr.clone());
        }
    }
}
