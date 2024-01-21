use rand::prelude::*;

use crate::aabb::{surrounding_box, Aabb};
use crate::hitable::{HitRecord, Hitable};
use crate::ray::Ray;
use crate::vec3::Vector3;

#[derive(Clone)]
pub struct HitableList { 
    hitable_list: Vec<Box<dyn Hitable + Send + Sync>>,
    aabb_box: Aabb,
    nor_hitable_list_len: f64,
}

impl HitableList {
    pub fn new() -> Self {
        HitableList{
            hitable_list: Vec::new(),
            aabb_box: Aabb{ b_min:[0.0, 0.0, 0.0], b_max:[0.0, 0.0, 0.0] },
            nor_hitable_list_len: 1.0,
        }
    }

    pub fn push<H: Hitable + 'static + Send + Sync>(&mut self, hitable: H) {
        self.aabb_box =
            if self.hitable_list.len() == 0 {
                hitable.bounding_box().clone()
            } else {
                surrounding_box(&self.aabb_box, hitable.bounding_box())
            };
        self.hitable_list.push(Box::new(hitable));
        self.nor_hitable_list_len = 1.0 / (self.hitable_list.len() as f64);
    }
}

impl ::std::ops::Deref for HitableList {
    type Target = Vec<Box<dyn Hitable + Send + Sync>>;

    fn deref(&self) -> &Vec<Box<dyn Hitable + Send + Sync>> {
        &self.hitable_list
    }
}

impl ::std::ops::DerefMut for HitableList {
    fn deref_mut(&mut self) -> &mut Vec<Box<dyn Hitable + Send + Sync>> {
        &mut self.hitable_list
    }
}

impl Hitable for HitableList {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut rec: Option<HitRecord> = None;
        let mut closer_so_far = t_max;
        for i in self.iter() {
            if let Some(temp_rec) = i.hit(r, t_min, closer_so_far) {
                closer_so_far = temp_rec.t;
                rec = Some(temp_rec);
            }
        }
        rec
    }

    fn bounding_box(&self) -> &Aabb {
        &self.aabb_box
    }

    fn pdf_value(&self, ray: &Ray) -> f64 {
        if let Some(_aabb_hit) = self.aabb_box.aabb_hit(ray, 0.00001, 10000.0) {
            let mut sum: f64 = 0.0;
            for i in self.iter() {
                sum += i.pdf_value(ray);
            }
            return sum * self.nor_hitable_list_len;
        }
        0.0
    }

    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        // we can clarify self.0.len() >= 1. after push some obj...
        let n = self.hitable_list.len();
        let mut rng = rand::thread_rng();
        let rand: f64 = rng.gen();

        let index: f64 = n as f64 * rand; // (1 * 0.9) as usize = 0
        self.hitable_list[index as usize].random(o)
    }
}
