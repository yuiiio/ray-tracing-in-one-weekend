use rand::prelude::*;

use crate::aabb::{surrounding_box, Aabb};
use crate::hitable::{HitRecord, Hitable};
use crate::ray::Ray;
use crate::vec3::Vector3;

#[derive(Clone)]
pub struct HitableList(Vec<Box<dyn Hitable + Send + Sync>>);

impl HitableList {
    pub fn new() -> Self {
        HitableList(Vec::new())
    }

    pub fn push<H: Hitable + 'static + Send + Sync>(&mut self, hitable: H) {
        self.0.push(Box::new(hitable))
    }

    pub fn from_vec(vec: Vec<Box<dyn Hitable + Send + Sync>>) -> Self {
        HitableList(vec)
    }
}

impl ::std::ops::Deref for HitableList {
    type Target = Vec<Box<dyn Hitable + Send + Sync>>;

    fn deref(&self) -> &Vec<Box<dyn Hitable + Send + Sync>> {
        &self.0
    }
}

impl ::std::ops::DerefMut for HitableList {
    fn deref_mut(&mut self) -> &mut Vec<Box<dyn Hitable + Send + Sync>> {
        &mut self.0
    }
}

impl Hitable for HitableList {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut rec: Option<HitRecord> = None;
        let mut closer_so_far = t_max;
        for i in self.iter() {
            if let Some(temp_rec) = i.hit(r, t_min, closer_so_far) {
                closer_so_far = temp_rec.get_t();
                rec = Some(temp_rec);
            }
        }
        rec
    }

    fn bounding_box(&self) -> Option<Aabb> {
        if self.0.len() < 1 {
            return None;
        }
        let mut temp_box: Aabb;
        match self.0[0].bounding_box() {
            Some(aabb) => temp_box = aabb,
            None => return None,
        }
        for i in self.iter().skip(1) {
            temp_box = surrounding_box(
                match i.bounding_box() {
                    Some(aabb) => aabb,
                    None => return None,
                },
                temp_box,
            );
        }
        Some(temp_box)
    }

    fn pdf_value(&self, o: &Vector3<f64>, v: &Vector3<f64>) -> f64 {
        let weight: f64 = 1.0 / self.0.len() as f64;
        let mut sum: f64 = 0.0;
        for i in self.iter() {
            sum += i.pdf_value(o, v) * weight;
        }
        sum
    }

    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        // we can clarify self.0.len() >= 1. after push some obj...
        let n = self.0.len();
        let mut rng = rand::thread_rng();
        let rand: f64 = rng.gen();

        let index: f64 = n as f64 * rand; // (1 * 0.9) as usize = 0
        self.0[index as usize].random(o)
    }
}
