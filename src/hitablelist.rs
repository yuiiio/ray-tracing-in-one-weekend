use crate::hitable::{HitRecord, Hitable};
use crate::ray::{Ray};
pub use crate::Object;

pub struct HitableList {
    pub list: Vec<Object>,
}

impl Hitable for HitableList {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut rec: Option<HitRecord> = None;
        let mut closer_so_far = t_max;
        for i in 0..self.list.len() {
            match self.list.get(i as usize).unwrap().call().unwrap().hit(r, t_min, closer_so_far) {
                Some(temp_rec) => {
                closer_so_far = temp_rec.t;
                rec = Some(temp_rec);
                },
                None => (),
            };
        };
        rec
    }
}
