use crate::hitable::{HitRecord, Hitable};
use crate::ray::{Ray};
use crate::aabb::{Aabb, surrounding_box};

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
            return None
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
                temp_box);
        }
        Some(temp_box)
    }
}
