use rand::prelude::*;

use crate::aabb::{Aabb, surrounding_box};
use crate::hitable::{Hitable, HitRecord};
use crate::ray::{Ray};
use crate::hitablelist::{HitableList};

pub struct BvhNode {
    bvh_node_box: Aabb,
    left: Box<dyn Hitable + Send + Sync>,
    right: Box<dyn Hitable + Send + Sync>,
}

impl BvhNode {
    pub fn new(hitable_list: &mut HitableList) -> BvhNode {
        let mut rng = rand::thread_rng();
        let x: f64 = rng.gen();
        let x: f64 = x * 3.0;
        let axis: usize = x as usize;
        //match axis {
        //    0 => ,
        //    1 => ,
        //    _ => ,
        //}
        let mut left_obj: Box<dyn Hitable + Send + Sync> = Box::new(Aabb::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]));
        let mut right_obj: Box<dyn Hitable + Send + Sync> = Box::new(Aabb::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]));
        let list_size = hitable_list.len();
        match list_size {
            1 => {
                left_obj = hitable_list[0];
                //right_obj = hitable_list[0];
                },
            2 => {
                //left_obj = hitable_list[0];
                //right_obj = hitable_list[1];
            },
            _ => {
                let mut a = HitableList::from_vec(hitable_list.split_off(list_size / 2));
                let mut b = hitable_list;
                left_obj = Box::new(BvhNode::new(&mut a));
                right_obj = Box::new(BvhNode::new(&mut b));
            },
        }
        let left_box = left_obj.bounding_box().expect("no bounding box in bvh_node constructor");
        let right_box = right_obj.bounding_box().expect("no bounding box in bvh_node constructor");
        BvhNode { bvh_node_box: surrounding_box(left_box, right_box),
            left: left_obj,
            right: right_obj
        }
    }
}

impl Hitable for BvhNode {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        match self.bvh_node_box.hit(r, t_min, t_max) {
            Some(_hit_rec) =>
                match self.left.hit(r, t_min, t_max) {
                    Some(left_rec) => match self.right.hit(r, t_min, t_max) {
                        Some(right_rec) => if left_rec.get_t() < right_rec.get_t() {
                            return Some(left_rec)
                        } else {
                            return Some(right_rec)
                        },
                        None => return Some(left_rec),
                    },
                    None => match self.right.hit(r, t_min, t_max) {
                        Some(right_rec) => return Some(right_rec),
                        None => return None,
                    }
                },
            None =>
                return None,
        }
    }

    fn bounding_box(&self) -> Option<Aabb> {
        Some(self.bvh_node_box)
    }
}