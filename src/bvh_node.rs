use rand::prelude::*;

use crate::aabb::{Aabb, surrounding_box};
use crate::hitable::{Hitable, HitRecord};
use crate::ray::{Ray};
use crate::hitablelist::{HitableList};
use crate::utils::{qsort};

#[derive(Clone)]
pub struct BvhNode {
    bvh_node_box: Aabb,
    left: Box<dyn Hitable + Send + Sync>,
    right: Box<dyn Hitable + Send + Sync>,
}

fn box_x_compare(a: &Box<dyn Hitable + Send + Sync>, b: &Box<dyn Hitable + Send + Sync>) -> bool {
    let box_a: Aabb = a.bounding_box().unwrap();
    let box_b: Aabb = b.bounding_box().unwrap();
    if box_a.min()[0] < box_b.min()[0] {
        return true;
    } else {
        return false;
    }
}

fn box_y_compare(a: &Box<dyn Hitable + Send + Sync>, b: &Box<dyn Hitable + Send + Sync>) -> bool {
    let box_a: Aabb = a.bounding_box().unwrap();
    let box_b: Aabb = b.bounding_box().unwrap();
    if box_a.min()[1] < box_b.min()[1] {
        return true;
    } else {
        return false;
    }
}

fn box_z_compare(a: &Box<dyn Hitable + Send + Sync>, b: &Box<dyn Hitable + Send + Sync>) -> bool {
    let box_a: Aabb = a.bounding_box().unwrap();
    let box_b: Aabb = b.bounding_box().unwrap();
    if box_a.min()[2] < box_b.min()[2] {
        return true;
    } else {
        return false;
    }
}

impl BvhNode {
    pub fn new(hitable_list: &mut HitableList) -> Self {
        let mut rng = rand::thread_rng();
        let x: f64 = rng.gen();
        let x: f64 = x * 3.0;
        let axis: usize = x as usize;
        match axis {
            0 => qsort(hitable_list, box_x_compare),
            1 => qsort(hitable_list, box_y_compare),
            _ => qsort(hitable_list, box_z_compare),
        }
        let list_size = hitable_list.len();
        let (left_obj, right_obj): (Box<dyn Hitable + Send + Sync>, Box<dyn Hitable + Send + Sync>) = match list_size {
            1 => {
                let left_obj = hitable_list[0].clone();
                let right_obj = hitable_list[0].clone();
                (left_obj, right_obj)
                },
            2 => {
                let left_obj = hitable_list[0].clone();
                let right_obj = hitable_list[1].clone();
                (left_obj, right_obj)
            },
            _ => {
                let mut a = HitableList::from_vec(hitable_list.split_off(list_size / 2));
                let mut b = hitable_list;
                let left_obj = Box::new(BvhNode::new(&mut a));
                let right_obj = Box::new(BvhNode::new(&mut b));
                (left_obj, right_obj)
            },
        };
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
            Some(_hit_rec) => {
                match self.left.hit(r, t_min, t_max) {
                    Some(left_rec) => {
                        match self.right.hit(r, t_min, t_max) {
                            Some(right_rec) => {
                                if left_rec.get_t() < right_rec.get_t() {
                                    return Some(left_rec)
                                } else {
                                    return Some(right_rec)
                                }
                            },
                            None => return Some(left_rec),
                        }
                    },
                    None => {
                        match self.right.hit(r, t_min, t_max) {
                            Some(right_rec) => return Some(right_rec),
                            None => return None,
                        }
                    },
                }
            },
            None => return None,
        }
    }

    fn bounding_box(&self) -> Option<Aabb> {
        Some(self.bvh_node_box.clone())
    }
}