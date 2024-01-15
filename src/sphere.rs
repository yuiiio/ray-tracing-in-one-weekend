use rand::prelude::*;
use std::f64::consts::PI;

use crate::aabb::Aabb;
use crate::hitable::{HitRecord, Hitable};
use crate::material::MaterialHandle;
use crate::onb::Onb;
use crate::ray::Ray;
use crate::vec3::{
    vec3_add_b, vec3_dot, vec3_mul_b, vec3_squared_length, vec3_sub, vec3_sub_b, Vector3,
};

#[derive(Clone)]
pub struct Sphere {
    center: Vector3<f64>,
    //radius: f64,
    nor_radius: f64, // for fn hit(), pre compute
    radius_sq: f64, // radius^2
    mat_ptr: MaterialHandle,
    aabb_box: Aabb,
    needs_uv: bool,
}

impl Sphere {
    pub fn new(center: Vector3<f64>, radius: f64, mat_ptr: MaterialHandle) -> Self {
        let aabb_box = Aabb {
            b_min: vec3_sub_b(&center, radius),
            b_max: vec3_add_b(&center, radius),
        };
        let nor_radius = 1.0 / radius;
        let radius_sq = radius * radius;
        let needs_uv = mat_ptr.needs_uv;
        Sphere {
            center,
            //radius,
            nor_radius,
            radius_sq,
            mat_ptr,
            aabb_box,
            needs_uv,
        }
    }
}

fn get_sphere_uv(point: Vector3<f64>) -> (f64, f64) {
    let u: f64 = 1.0 - (point[2].atan2(point[0]) + PI) / (2.0 * PI); // atan(z/x)
    let v: f64 = (point[1].asin() + (PI / 2.0)) / PI;
    (u, v)
}

impl Hitable for Sphere {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let rec: Option<HitRecord> = None;
        let oc = vec3_sub(&r.origin, &self.center);
        let a = vec3_dot(&r.direction, &r.direction);
        let b = 2.0 * vec3_dot(&r.direction, &oc);
        let c = vec3_dot(&oc, &oc) - self.radius_sq;
        let descriminant = b * b - 4.0 * a * c;
        if descriminant >= 0.0 {
            let temp = (-b - descriminant.sqrt()) / (2.0 * a);
            if temp < t_max && temp > t_min {
                let point = r.point_at_parameter(temp);
                let nnormal = vec3_mul_b(&vec3_sub(&point, &self.center), self.nor_radius);
                let uv = if self.needs_uv == true {
                    get_sphere_uv(nnormal)
                } else {
                    (0.0, 0.0)
                };
                return Some(HitRecord {
                    t: temp,
                    uv,
                    p: point,
                    normal: nnormal,
                    mat_ptr: &self.mat_ptr,
                });
            }
            let temp = (-b + descriminant.sqrt()) / (2.0 * a);
            if temp < t_max && temp > t_min {
                let point = r.point_at_parameter(temp);
                let nnormal = vec3_mul_b(&vec3_sub(&point, &self.center), self.nor_radius);
                let uv = if self.needs_uv == true {
                    get_sphere_uv(nnormal)
                } else {
                    (0.0, 0.0)
                };
                return Some(HitRecord {
                    t: temp,
                    uv,
                    p: point,
                    normal: nnormal,
                    mat_ptr: &self.mat_ptr,
                });
            }
        }
        rec
    }

    fn bounding_box<'a>(&'a self) -> Option<&'a Aabb> {
        Some(&self.aabb_box)
    }

    fn pdf_value(&self, o: &Vector3<f64>, v: &Vector3<f64>) -> f64 {
        if let Some(_aabb_hit) = self.aabb_box.aabb_hit(&Ray{ origin: *o, direction: *v }, 0.00001, 10000.0)  {

            match self.hit(&Ray{ origin: *o, direction: *v }, 0.00001, 10000.0) {
                Some(_rec) => {
                    let distabce_squared: f64 = vec3_squared_length(&vec3_sub(&self.center, o));
                    let cos_theta_max: f64 = (1.0 - (self.radius_sq / distabce_squared)).sqrt();
                    return 1.0 / (2.0 * PI * (1.0 - cos_theta_max));
                },
                None => return 0.0,
            }
        } else {
            return 0.0;
        }
    }
    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        let direction = vec3_sub(&self.center, o);
        let uvw = Onb::build_from_w(&direction);

        let distabce_squared = vec3_squared_length(&direction);
        uvw.local(&random_to_sphere(self.radius_sq, distabce_squared))
    }
}

fn random_to_sphere(radius_sq: f64, distabce_squared: f64) -> Vector3<f64> {
    let mut rng = rand::thread_rng();
    let r1: f64 = rng.gen();
    let r2: f64 = rng.gen();
    let cos_theta_max = (1.0 - (radius_sq / distabce_squared)).sqrt();
    let z = 1.0 - r2 * (1.0 - cos_theta_max);
    let a = 2.0 * PI * r1;
    let b = (1.0 - z.powi(2)).sqrt();
    let x = a.cos() * b;
    let y = a.sin() * b;
    [x, y, z]
}
