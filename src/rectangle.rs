use rand::prelude::*;

use crate::aabb::Aabb;
use crate::hitable::{HitRecord, Hitable};
use crate::material::MaterialHandle;
use crate::onb::Onb;
use crate::quotation::Rotation;
use crate::ray::Ray;
use crate::vec3::{vec3_add, vec3_dot, vec3_mul_b, vec3_sub, vec3_unit_vector_f64, Vector3};

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
    width: f64,
    nor_width: f64,
    y0: f64,
    y1: f64,
    height: f64,
    nor_height: f64,
    k: f64,
    axis: AxisType,
    mat_ptr: MaterialHandle,
    area: f64,
    aabb_box: Aabb,
    needs_uv: bool,
    normal: Vector3<f64>,
    onb_uv: (Vector3<f64>, Vector3<f64>),
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
        let area: f64 = (x1 - x0) * (y1 - y0);
        let (aabb_box, mut normal) = match axis {
            AxisType::Kxy => (
                Aabb {
                    b_min: [x0, y0, k - 0.0001],
                    b_max: [x1, y1, k + 0.0001],
                },
                [0.0, 0.0, 1.0],
            ),
            AxisType::Kxz => (
                Aabb {
                    b_min: [x0, k - 0.0001, y0],
                    b_max: [x1, k + 0.0001, y1],
                },
                [0.0, 1.0, 0.0],
            ),
            AxisType::Kyz => (
                Aabb {
                    b_min: [k - 0.0001, x0, y0],
                    b_max: [k + 0.0001, x1, y1],
                },
                [1.0, 0.0, 0.0],
            ),
        };
        if flip_normal == true {
            normal = vec3_mul_b(&normal, -1.0);
        }
        let width = x1 - x0;
        let height = y1 - y0;
        let nor_width = 1.0 / (x1 - x0);
        let nor_height = 1.0 / (y1 - y0);
        let needs_uv = mat_ptr.needs_uv;
        let onb = Onb::build_from_w(&normal);
        Rect {
            x0,
            x1,
            width,
            nor_width,
            y0,
            y1,
            height,
            nor_height,
            k,
            axis,
            mat_ptr,
            area,
            aabb_box,
            needs_uv,
            normal,
            onb_uv: (onb.u, onb.v),
        }
    }
}

// Rect hit is valid even ray direction.
// Rect != actual Surface.

impl Rect {
    fn rect_hit(&self, r: &Ray, r_dir_div: &[f64; 3], t_min: f64, t_max: f64) -> Option<HitRecord> {
        let (xi, yi, zi): (usize, usize, usize) = match self.axis {
            AxisType::Kxy => (0, 1, 2),
            AxisType::Kxz => (0, 2, 1),
            AxisType::Kyz => (1, 2, 0),
        };

        let t = (self.k - r.origin[zi]) * r_dir_div[zi];
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
            (
                (x - self.x0) * self.nor_width,
                (y - self.y0) * self.nor_height,
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
            mat_ptr: &self.mat_ptr,
            onb_uv: Some(&self.onb_uv),
        })
    }

    fn rect_pdf_value(&self, ray: &Ray, r_dir_div: &[f64; 3]) -> f64 {
        if let Some(rec) = self.rect_hit(ray, r_dir_div, 0.00001, 10000.0) {
            let distance_squared = rec.t.powi(2);
            let cosine = vec3_dot(&ray.direction, &rec.normal).abs();
            return distance_squared / (cosine * self.area);
        }
        0.0
    }
}

impl Hitable for Rect {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        self.rect_hit(r, &r.get_inv_dir(), t_min, t_max)
    }

    fn bounding_box(&self) -> &Aabb {
        &self.aabb_box
    }

    fn pdf_value(&self, ray: &Ray) -> f64 {
        self.rect_pdf_value(ray, &ray.get_inv_dir())
    }

    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        let mut rng = rand::thread_rng();
        let rng_x: f64 = rng.gen();
        let rng_y: f64 = rng.gen();
        let random_point = match self.axis {
            AxisType::Kxy => [
                self.x0 + rng_x * self.width,
                self.y0 + rng_y * self.height,
                self.k,
            ],
            AxisType::Kxz => [
                self.x0 + rng_x * self.width,
                self.k,
                self.y0 + rng_y * self.height,
            ],
            AxisType::Kyz => [
                self.k,
                self.x0 + rng_x * self.width,
                self.y0 + rng_y * self.height,
            ],
        };
        //TODO: need all rewrite
        vec3_unit_vector_f64(&vec3_sub(&random_point, o))
    }

    fn rotate_onb(&mut self, quat: &Rotation) -> () {
        self.normal = quat.rotate(&self.normal);
        let onb = Onb::build_from_w(&self.normal);
        self.onb_uv = (onb.u, onb.v);
    }
}

#[derive(Clone)]
pub struct Boxel {
    rect: [Rect; 6],
    aabb_box: Aabb,
    center: [f64; 3],
    mat_ptr: MaterialHandle,
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
        let center = vec3_mul_b(&vec3_add(&b_min, &b_max), 0.5);
        let aabb_box = Aabb { b_min, b_max };
        Boxel {
            rect,
            aabb_box,
            center,
            mat_ptr,
        }
    }
}

impl Hitable for Boxel {
    // needs to check from out-side and inside(eg: glass) rays.
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut tmin = t_min;
        let mut tmax = t_max;

        let p0 = self.aabb_box.b_min;
        let p1 = self.aabb_box.b_max;
        let r_dir_div = ray.get_inv_dir();
        let mut norm_dir = [0.0; 3];
        for i in 0..3 {
            let mut t0 = (p0[i] - ray.origin[i]) * r_dir_div[i];
            let mut t1 = (p1[i] - ray.origin[i]) * r_dir_div[i];

            let mut flip = -1.0;
            if r_dir_div[i].is_sign_negative() {
                core::mem::swap(&mut t0, &mut t1);
                flip = 1.0;
            }
            if tmin < t0 {
                let mut norm_tmp = [0.0; 3];
                norm_tmp[i] = flip;
                norm_dir = norm_tmp;
                tmin = t0;
            }
            tmax = tmax.min(t1);
        }
        if tmax < tmin {
            return None;
        }

        let p = ray.point_at_parameter(tmin);
        Some(HitRecord {
            t: tmin,
            uv: (0.0, 0.0), // TODO
            p,
            normal: norm_dir,
            mat_ptr: &self.mat_ptr,
            onb_uv: None, //TODO //Some(&self.onb_uv),
        })

        /*
        let mut rect_index_0 = [5, 4, 3];
        let mut rect_index_1 = [2, 1, 0];

        for i in 0..3 {
            if r_dir_div[i].is_sign_negative() {
                core::mem::swap(&mut rect_index_0[i], &mut rect_index_1[i]);
            }
            if let Some(hit_rec) =
                self.rect[rect_index_0[i]].rect_hit(ray, &r_dir_div, t_min, t_max)
            {
                return Some(hit_rec);
            }
        }
        // up is enough when only care out-side rays
        // down is needed for inside rays.
        for i in 0..3 {
            if let Some(hit_rec) =
                self.rect[rect_index_1[i]].rect_hit(ray, &r_dir_div, t_min, t_max)
            {
                return Some(hit_rec);
            }
        }
        None
        */
    }

    fn bounding_box(&self) -> &Aabb {
        &self.aabb_box
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
            if r_dir_div[0] > 0.0 {
                pdf_sum += self.rect[5].rect_pdf_value(ray, &r_dir_div);
            } else {
                pdf_sum += self.rect[2].rect_pdf_value(ray, &r_dir_div);
            }
            if r_dir_div[1] > 0.0 {
                pdf_sum += self.rect[4].rect_pdf_value(ray, &r_dir_div);
            } else {
                pdf_sum += self.rect[1].rect_pdf_value(ray, &r_dir_div);
            }
            if r_dir_div[2] > 0.0 {
                pdf_sum += self.rect[3].rect_pdf_value(ray, &r_dir_div);
            } else {
                pdf_sum += self.rect[0].rect_pdf_value(ray, &r_dir_div);
            }
            pdf_sum * DIV6
        } else {
            0.0
        }
    }

    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        let mut rng = rand::thread_rng();
        let rand: f64 = rng.gen();
        let o_dir = vec3_sub(&self.center, o);
        let random_handle = (rand * 3.0) as usize;
        if o_dir[random_handle] > 0.0 {
            self.rect[5 - random_handle].random(o)
        } else {
            self.rect[2 - random_handle].random(o)
        }
    }

    fn rotate_onb(&mut self, quat: &Rotation) -> () {
        for i in 0..6 {
            self.rect[i].rotate_onb(quat);
        }
    }
}
