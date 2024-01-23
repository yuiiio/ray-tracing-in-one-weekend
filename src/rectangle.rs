use rand::prelude::*;

use crate::aabb::Aabb;
use crate::hitable::{HitRecord, Hitable};
use crate::material::MaterialHandle;
use crate::ray::Ray;
use crate::vec3::{vec3_dot, vec3_mul_b, vec3_sub, Vector3, vec3_unit_vector_f64};
use crate::onb::Onb;

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
    onb: Onb,
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
        let (aabb_box, normal) = match axis {
            AxisType::Kxy => (
                Aabb {
                b_min: [x0, y0, k - 0.0001],
                b_max: [x1, y1, k + 0.0001],
                }, 
                [0.0, 0.0, 1.0] 
            ),
            AxisType::Kxz => (
                Aabb {
                b_min: [x0, k - 0.0001, y0],
                b_max: [x1, k + 0.0001, y1],
                },
                [0.0, 1.0, 0.0]
            ),
            AxisType::Kyz => (
                Aabb {
                b_min: [k - 0.0001, x0, y0],
                b_max: [k + 0.0001, x1, y1],
                },
                [1.0, 0.0, 0.0]
            ),
        };
        let width = x1 - x0;
        let height = y1 - y0;
        let nor_width = 1.0 / (x1 - x0);
        let nor_height = 1.0 / (y1 - y0);
        let needs_uv = mat_ptr.needs_uv;
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
            onb: Onb::build_from_w(&normal),
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
            ((x - self.x0) * self.nor_width,
            (y - self.y0) * self.nor_height)
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
            onb: Some(&self.onb),
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
}

#[derive(Clone)]
pub struct FlipNormals {
    shape: Rect,
    normal: Vector3<f64>,
    onb: Onb,
}

impl FlipNormals {
    pub fn new(shape: Rect) -> Self {
        let normal = vec3_mul_b(&shape.normal, -1.0);
        FlipNormals{
            shape,
            normal,
            onb: Onb::build_from_w(&normal),
        }
    }
}

impl Hitable for FlipNormals {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        match self.shape.hit(r, t_min, t_max) {
            Some(hit) => Some(HitRecord {
                t: hit.t,
                uv: hit.uv,
                p: hit.p,
                normal: self.normal,
                mat_ptr: hit.mat_ptr,
                onb: Some(&self.onb),
            }),
            None => None,
        }
    }

    fn bounding_box(&self) -> &Aabb {
        self.shape.bounding_box()
    }

    fn pdf_value(&self, ray: &Ray) -> f64 {
        self.shape.pdf_value(ray)
    }

    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        self.shape.random(o)
    }
}

#[derive(Clone)]
pub struct Boxel {
    rect: [Rect; 3],
    flip_rect: [FlipNormals; 3],
    aabb_box: Aabb,
}

impl Boxel {
    pub fn new(p0: Vector3<f64>, p1: Vector3<f64>, mat_ptr: MaterialHandle) -> Self {
        let b_min = p0;
        let b_max = p1;
        let rect: [Rect; 3] = [
            Rect::new(
                p0[0],
                p1[0],
                p0[1],
                p1[1],
                p1[2],
                AxisType::Kxy,
                mat_ptr.clone(),
                ),
                Rect::new(
                    p0[0],
                    p1[0],
                    p0[2],
                    p1[2],
                    p1[1],
                    AxisType::Kxz,
                    mat_ptr.clone(),
                    ),
                    Rect::new(
                        p0[1],
                        p1[1],
                        p0[2],
                        p1[2],
                        p1[0],
                        AxisType::Kyz,
                        mat_ptr.clone(),
                        ),
                    ];
        let flip_rect: [FlipNormals; 3] = [
            FlipNormals::new(Rect::new(
                    p0[0],
                    p1[0],
                    p0[1],
                    p1[1],
                    p0[2],
                    AxisType::Kxy,
                    mat_ptr.clone(),
                    )),
                    FlipNormals::new(Rect::new(
                            p0[0],
                            p1[0],
                            p0[2],
                            p1[2],
                            p0[1],
                            AxisType::Kxz,
                            mat_ptr.clone(),
                            )),
                            FlipNormals::new(Rect::new(
                                    p0[1],
                                    p1[1],
                                    p0[2],
                                    p1[2],
                                    p0[0],
                                    AxisType::Kyz,
                                    mat_ptr.clone(),
                                    )),
                            ];
        let aabb_box = Aabb{b_min, b_max};
        Boxel { rect, flip_rect, aabb_box }
    }
}

impl Hitable for Boxel {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut hit_min_t = t_max;
        let mut hit_count: usize = 0;
        let mut rec: Option<HitRecord> = None;
        for i in 0..3 {
            if let Some(hit_rec) = self.rect[i].hit(r, t_min, hit_min_t) {
                hit_min_t = hit_rec.t;
                hit_count += 1;
                rec = Some(hit_rec);
            }
            if let Some(hit_rec) = self.flip_rect[i].hit(r, t_min, hit_min_t) {
                hit_min_t = hit_rec.t;
                hit_count += 1;
                rec = Some(hit_rec);
            }
            if hit_count == 2 {
                break;
            }
        }
        rec
    }

    fn bounding_box(&self) -> &Aabb {
        &self.aabb_box
    }

    fn pdf_value(&self, ray: &Ray) -> f64 {
        // TODO: we needs actual pdf hit surface, now return avg all surface
        if let Some(_aabb_hit) = self.aabb_box.aabb_hit(ray, 0.00001, 10000.0)  {
            const DIV6: f64 = 1.0 / 6.0;
            ( self.rect[0].pdf_value(ray)
              + self.rect[1].pdf_value(ray)
              + self.rect[2].pdf_value(ray)
              + self.flip_rect[0].pdf_value(ray)
              + self.flip_rect[1].pdf_value(ray)
              + self.flip_rect[2].pdf_value(ray)
            ) * DIV6
        } else {
            0.0
        }
    }

    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        let mut rng = rand::thread_rng();
        let rand: f64 = rng.gen();
        let random_handle = (rand * 3.0) as usize;
        let rand2: f64 = rng.gen();
        if rand2 < 0.5 {
            self.rect[random_handle].random(o)
        } else {
            self.flip_rect[random_handle].random(o)
        }
    }
}
