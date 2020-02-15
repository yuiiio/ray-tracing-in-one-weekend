use crate::hitable::{HitRecord, Hitable};
use crate::ray::{Ray};

pub struct HitableList(Vec<Box<dyn Hitable>>);

impl HitableList {
    pub fn new() -> Self {
        HitableList(Vec::new())
    }

    pub fn push<H: Hitable + 'static>(&mut self, hitable: H) {
        self.0.push(Box::new(hitable))
    }
}

impl ::std::ops::Deref for HitableList {
    type Target = Vec<Box<dyn Hitable>>;

    fn deref(&self) -> &Vec<Box<dyn Hitable>> {
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
}
