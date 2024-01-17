use rand::prelude::*;

use crate::aabb::Aabb;
use crate::hitable::{HitRecord, Hitable};
use crate::material::MaterialHandle;
use crate::ray::Ray;
use crate::vec3::{vec3_dot, vec3_mul_b, vec3_sub, Vector3, vec3_unit_vector_f64};

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
            AxisType::Kxy => Aabb {
                b_min: [x0, y0, k - 0.0001],
                b_max: [x1, y1, k + 0.0001],
            },
            AxisType::Kxz => Aabb {
                b_min: [x0, k - 0.0001, y0],
                b_max: [x1, k + 0.0001, y1],
            },
            AxisType::Kyz => Aabb {
                b_min: [k - 0.0001, x0, y0],
                b_max: [k + 0.0001, x1, y1],
            },
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
        }
    }
}

// Rect hit is valid even ray direction.
// Rect != actual Surface.
impl Hitable for Rect {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let (xi, yi, zi, nnormal): (usize, usize, usize, Vector3<f64>) = match self.axis {
            AxisType::Kxy => (0, 1, 2, [0.0, 0.0, 1.0]),
            AxisType::Kxz => (0, 2, 1, [0.0, 1.0, 0.0]),
            AxisType::Kyz => (1, 2, 0, [1.0, 0.0, 0.0]),
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

        let (u, v) = if self.needs_uv == true {
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
            normal: nnormal,
            mat_ptr: &self.mat_ptr,
        })
    }

    fn bounding_box<'a>(&'a self) -> Option<&'a Aabb> {
        Some(&self.aabb_box)
    }

    fn pdf_value(&self, o: &Vector3<f64>, v: &Vector3<f64>) -> f64 {
        match self.hit(&Ray{ origin: *o, direction: *v }, 0.00001, 10000.0) {
            Some(rec) => {
                let distance_squared = rec.t.powi(2);
                let cosine = vec3_dot(v, &rec.normal).abs();
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
}

impl FlipNormals {
    pub fn new(shape: Rect) -> Self {
        FlipNormals { shape }
    }
}

impl Hitable for FlipNormals {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        match self.shape.hit(r, t_min, t_max) {
            Some(hit) => Some(HitRecord {
                t: hit.t,
                uv: hit.uv,
                p: hit.p,
                normal: vec3_mul_b(&hit.normal, -1.0),
                mat_ptr: hit.mat_ptr,
            }),
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
        if let Some(hit_rec0) = self.rect[0].hit(r, t_min, hit_min_t) {
            hit_min_t = hit_rec0.t;
            if let Some(hit_rec1) = self.rect[1].hit(r, t_min, hit_min_t) {
                return Some(hit_rec1);
            }
            if let Some(hit_rec2) = self.rect[2].hit(r, t_min, hit_min_t) {
                return Some(hit_rec2);
            }
            if let Some(hit_rec3) = self.flip_rect[0].hit(r, t_min, hit_min_t) {
                return Some(hit_rec3);
            }
            if let Some(hit_rec4) = self.flip_rect[1].hit(r, t_min, hit_min_t) {
                return Some(hit_rec4);
            }
            if let Some(hit_rec5) = self.flip_rect[2].hit(r, t_min, hit_min_t) {
                return Some(hit_rec5);
            }
            return Some(hit_rec0); // rect[0] was min hit
        } else {
            if let Some(hit_rec1) = self.rect[1].hit(r, t_min, hit_min_t) {
                hit_min_t = hit_rec1.t;
                if let Some(hit_rec2) = self.rect[2].hit(r, t_min, hit_min_t) {
                    return Some(hit_rec2);
                }
                if let Some(hit_rec3) = self.flip_rect[0].hit(r, t_min, hit_min_t) {
                    return Some(hit_rec3);
                }
                if let Some(hit_rec4) = self.flip_rect[1].hit(r, t_min, hit_min_t) {
                    return Some(hit_rec4);
                }
                if let Some(hit_rec5) = self.flip_rect[2].hit(r, t_min, hit_min_t) {
                    return Some(hit_rec5);
                }
                return Some(hit_rec1); // rect[1] was min hit
            } else {
                if let Some(hit_rec2) = self.rect[2].hit(r, t_min, hit_min_t) {
                    hit_min_t = hit_rec2.t;
                    if let Some(hit_rec3) = self.flip_rect[0].hit(r, t_min, hit_min_t) {
                        return Some(hit_rec3);
                    }
                    if let Some(hit_rec4) = self.flip_rect[1].hit(r, t_min, hit_min_t) {
                        return Some(hit_rec4);
                    }
                    if let Some(hit_rec5) = self.flip_rect[2].hit(r, t_min, hit_min_t) {
                        return Some(hit_rec5);
                    }
                    return Some(hit_rec2); // rect[2] was min hit
                } else {
                    if let Some(hit_rec3) = self.flip_rect[0].hit(r, t_min, hit_min_t) {
                        hit_min_t = hit_rec3.t;
                        if let Some(hit_rec4) = self.flip_rect[1].hit(r, t_min, hit_min_t) {
                            return Some(hit_rec4);
                        }
                        if let Some(hit_rec5) = self.flip_rect[2].hit(r, t_min, hit_min_t) {
                            return Some(hit_rec5);
                        }
                        return Some(hit_rec3); // flip_rect[0] was min hit
                    } else {
                        if let Some(hit_rec4) = self.flip_rect[1].hit(r, t_min, hit_min_t) {
                            hit_min_t = hit_rec4.t;
                            if let Some(hit_rec5) = self.flip_rect[2].hit(r, t_min, hit_min_t) {
                                return Some(hit_rec5);
                            }
                            return Some(hit_rec4); // flip_rect[1] was min hit
                        } else {
                            if let Some(hit_rec5) = self.flip_rect[2].hit(r, t_min, hit_min_t) {
                                return Some(hit_rec5);
                            }
                            return None;
                        }
                    }
                }
            }
        }
    }

    fn bounding_box<'a>(&'a self) -> Option<&'a Aabb> {
        Some(&self.aabb_box)
    }

    fn pdf_value(&self, o: &Vector3<f64>, v: &Vector3<f64>) -> f64 {
        // TODO: we needs actual pdf hit surface, now return avg all surface
        if let Some(_aabb_hit) = self.aabb_box.aabb_hit(&Ray{ origin: *o, direction: *v }, 0.00001, 10000.0)  {

            const DIV6: f64 = 1.0 / 6.0;
            return (
                self.rect[0].pdf_value(o, v)
                + self.rect[1].pdf_value(o, v)
                + self.rect[2].pdf_value(o, v)
                + self.flip_rect[0].pdf_value(o, v)
                + self.flip_rect[1].pdf_value(o, v)
                + self.flip_rect[2].pdf_value(o, v)
                ) * DIV6
        } else {
            return 0.0;
        }
    }

    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        let mut rng = rand::thread_rng();
        let rand: f64 = rng.gen();
        let random_handle = (rand * 3.0) as usize;
        let rand2: f64 = rng.gen();
        if rand2 < 0.5 {
            return self.rect[random_handle].random(o); 
        } else {
            return self.flip_rect[random_handle].random(o); 
        }
    }
}
