use crate::hitable::{HitRecord, Hitable};
use crate::ray::{Ray};
use crate::aabb::{aabb};

pub struct HitableList(Vec<Box<dyn Hitable + Send + Sync>>);

impl HitableList {
    pub fn new() -> Self {
        HitableList(Vec::new())
    }

    pub fn push<H: Hitable + 'static + Send + Sync>(&mut self, hitable: H) {
        self.0.push(Box::new(hitable))
    }
}

impl ::std::ops::Deref for HitableList {
    type Target = Vec<Box<dyn Hitable + Send + Sync>>;

    fn deref(&self) -> &Vec<Box<dyn Hitable + Send + Sync>> {
        &self.0
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
    
    fn bounding_box(&self, t0: f64, t1: f64) -> Option<aabb> {
        if self.0.len() < 1 {
            return None
        }
        let mut temp_box: aabb;
        match self.0[0].bounding_box(t0, t1) {
            Some(aabb) => temp_box = aabb,
            None => return None,
        }
        for i in self.iter().skip(1) {
            temp_box = surrounding_box(
                match i.bounding_box(t0, t1) {
                    Some(aabb) => aabb,
                    None => return None,
                },
                temp_box);
        }
        Some(temp_box)
    }
}

fn surrounding_box(box0: aabb, box1: aabb) -> aabb {
    let min = [min(box0.min()[0], box1.min()[0]),
                min(box0.min()[1], box1.min()[1]),
                min(box0.min()[2], box1.min()[2])];
    let max = [max(box0.max()[0], box1.max()[0]),
                max(box0.max()[1], box1.max()[1]),
                max(box0.max()[2], box1.max()[2])];
    aabb::new(min, max)
}

fn max(a: f64, b: f64) -> f64 {
    if a < b {
        b
    } else {
        a
    }
}

fn min(a: f64, b: f64) -> f64 {
    if a > b {
        b
    } else {
        a
    }
}