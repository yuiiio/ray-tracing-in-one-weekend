use rand::prelude::*;

use crate::aabb::{surrounding_box, Aabb};
use crate::hitable::{HitRecord, Hitable};
use crate::hitablelist::HitableList;
use crate::ray::Ray;
use crate::vec3::Vector3;

#[derive(Clone)]
pub struct BvhNode {
    bvh_node_box: Aabb,
    left: Box<dyn Hitable + Send + Sync>,
    right: Box<dyn Hitable + Send + Sync>,
}

fn box_x_compare(a: &Box<dyn Hitable + Send + Sync>, b: &Box<dyn Hitable + Send + Sync>) -> bool {
    let box_a: Aabb = a.bounding_box().unwrap();
    let box_b: Aabb = b.bounding_box().unwrap();
    if box_a.b_min()[0] < box_b.b_min()[0] {
        return true;
    } else {
        return false;
    }
}

fn box_y_compare(a: &Box<dyn Hitable + Send + Sync>, b: &Box<dyn Hitable + Send + Sync>) -> bool {
    let box_a: Aabb = a.bounding_box().unwrap();
    let box_b: Aabb = b.bounding_box().unwrap();
    if box_a.b_min()[1] < box_b.b_min()[1] {
        return true;
    } else {
        return false;
    }
}

fn box_z_compare(a: &Box<dyn Hitable + Send + Sync>, b: &Box<dyn Hitable + Send + Sync>) -> bool {
    let box_a: Aabb = a.bounding_box().unwrap();
    let box_b: Aabb = b.bounding_box().unwrap();
    if box_a.b_min()[2] < box_b.b_min()[2] {
        return true;
    } else {
        return false;
    }
}

fn dmerge(
    vec: &mut Vec<usize>,
    stock_vec: &mut Vec<usize>,
    compare: fn(&Box<dyn Hitable + Send + Sync>, &Box<dyn Hitable + Send + Sync>) -> bool,
    hitable_list: &HitableList,
    left: usize,
    mid: usize,
    right: usize,
    ) {
    let mut i: usize = left;
    let mut j: usize = mid;
    let mut k: usize = 0;
    let mut l: usize = 0;

    while i < mid && j < right {
        if compare(&hitable_list[vec[i]], &hitable_list[vec[j]]) {
            stock_vec[k] = vec[i];
            i = i + 1;
        } else {
            stock_vec[k] = vec[j];
            j = j + 1;
        }
        k = k + 1;
    }

    if i == mid {
        while j < right {
            stock_vec[k] = vec[j];
            j = j + 1;
            k = k + 1;
        }
    } else {
        while i < mid {
            stock_vec[k] = vec[i];
            i = i + 1;
            k = k + 1;
        }
    }

    while l < k {
        vec[left + l] = stock_vec[l];
        l = l + 1;
    }
}

pub fn dmerge_sort(
    vec: &mut Vec<usize>,
    stock_vec: &mut Vec<usize>,
    compare: fn(&Box<dyn Hitable + Send + Sync>, &Box<dyn Hitable + Send + Sync>) -> bool,
    hitable_list: &HitableList,
    left: usize,
    right: usize,
    ) {
    if (left == right) || (left + 1 == right) {
        return;
    }
    let mid: usize = (left + right) / 2;
    dmerge_sort(vec, stock_vec, compare, hitable_list, left, mid);
    dmerge_sort(vec, stock_vec, compare, hitable_list, mid, right);
    dmerge(vec, stock_vec, compare, hitable_list, left, mid, right);
}

pub fn dmerge_sort_wrap(
    vec: &mut Vec<usize>,
    compare: fn(&Box<dyn Hitable + Send + Sync>, &Box<dyn Hitable + Send + Sync>) -> bool,
    hitable_list: &HitableList,
    ) {
    let mut stock_vec = vec.clone();
    let len = vec.len();
    dmerge_sort(vec, &mut stock_vec, compare, hitable_list, 0, len);
}

#[derive(Clone)]
enum Axis {
    X,
    Y,
    Z,
}

fn build_bvh(hitable_list: &HitableList, handle: &mut Vec<usize>, pre_sort_axis: Axis) -> BvhNode {
    let handle_size = handle.len();
    let (left_obj, right_obj): (
        Box<dyn Hitable + Send + Sync>,
        Box<dyn Hitable + Send + Sync>,
        ) = match handle_size {
        1 => {
            let left_obj = hitable_list[handle[0]].clone();
            let right_obj = hitable_list[handle[0]].clone();
            (left_obj, right_obj)
        }
        2 => {
            let left_obj = hitable_list[handle[0]].clone();
            let right_obj = hitable_list[handle[1]].clone();
            (left_obj, right_obj)
        }
        _ => {
            /*
               let mut rng = rand::thread_rng();
               let x: f64 = rng.gen();
               let x: f64 = x * 3.0;
               let axis: usize = x as usize;
               match axis {
               0 => dmerge_sort_wrap(handle, box_x_compare, hitable_list),
               1 => dmerge_sort_wrap(handle, box_y_compare, hitable_list),
               _ => dmerge_sort_wrap(handle, box_z_compare, hitable_list),
               }
               */

            let mut handle_x: Vec<usize> = handle.clone();
            let mut handle_y: Vec<usize> = handle.clone();
            let mut handle_z: Vec<usize> = handle.clone();

            match pre_sort_axis {
                Axis::X => {
                    dmerge_sort_wrap(&mut handle_y, box_y_compare, hitable_list);
                    dmerge_sort_wrap(&mut handle_z, box_z_compare, hitable_list);
                },
                Axis::Y => {
                    dmerge_sort_wrap(&mut handle_x, box_x_compare, hitable_list);
                    dmerge_sort_wrap(&mut handle_z, box_z_compare, hitable_list);
                },
                Axis::Z => {
                    dmerge_sort_wrap(&mut handle_x, box_x_compare, hitable_list);
                    dmerge_sort_wrap(&mut handle_y, box_y_compare, hitable_list);
                },
            }

            let x_max: f64 = hitable_list[handle_x[handle.len() - 1]]
                .bounding_box()
                .unwrap()
                .b_max()[0]
                - hitable_list[handle_x[0]].bounding_box().unwrap().b_min()[0];
            let y_max: f64 = hitable_list[handle_y[handle.len() - 1]]
                .bounding_box()
                .unwrap()
                .b_max()[1]
                - hitable_list[handle_y[0]].bounding_box().unwrap().b_min()[1];
            let z_max: f64 = hitable_list[handle_z[handle.len() - 1]]
                .bounding_box()
                .unwrap()
                .b_max()[2]
                - hitable_list[handle_z[0]].bounding_box().unwrap().b_min()[2];

            let sorted_axis: Axis;
            let distance_max: f64;

            let mut handle = if x_max < y_max {
                if y_max < z_max {
                    distance_max = z_max;
                    sorted_axis = Axis::Z;
                    handle_z
                } else {
                    distance_max = y_max;
                    sorted_axis = Axis::Y;
                    handle_y
                }
            } else {
                if x_max < z_max {
                    distance_max = z_max;
                    sorted_axis = Axis::Z;
                    handle_z
                } else {
                    distance_max = x_max;
                    sorted_axis = Axis::X;
                    handle_x
                }
            };

            let distance_avg: f64 = distance_max / 2.0 as f64;

            // found distance_max split point
            let mut handle_size_log: usize = handle_size.ilog2() as usize;
            let handle_size_log_pow: usize = 2usize.pow(handle_size_log as u32);
            if handle_size_log != handle_size_log_pow {
                handle_size_log = handle_size_log + 1;
            }

            let mut split_point: usize = handle_size / 2;

            for i in (0..(handle_size_log-2)).rev()
            {
                match sorted_axis {
                    Axis::X => {
                        if hitable_list[handle[split_point]].bounding_box().unwrap().b_min()[0] < distance_avg
                        {
                            split_point = split_point + 2usize.pow(i as u32);
                        } else {
                            split_point = split_point - 2usize.pow(i as u32);
                        }
                    },
                    Axis::Y => {
                        if hitable_list[handle[split_point]].bounding_box().unwrap().b_min()[1] < distance_avg
                        {
                            split_point = split_point + 2usize.pow(i as u32);
                        } else {
                            split_point = split_point - 2usize.pow(i as u32);
                        }
                    },
                    Axis::Z => {
                        if hitable_list[handle[split_point]].bounding_box().unwrap().b_min()[2] < distance_avg
                        {
                            split_point = split_point + 2usize.pow(i as u32);
                        } else {
                            split_point = split_point - 2usize.pow(i as u32);
                        }
                    },
                }
            }

            //println!("handle_size/2: {}, split_point: {}", (handle_size/2) as usize, split_point);
            if (handle_size/2) as usize != split_point {
                println!("hadnle_size/2 != split_point");
            }

            let mut a = handle.split_off(split_point);
            let mut b = handle;

            let left_obj = Box::new(build_bvh(hitable_list, &mut a, sorted_axis.clone()));
            let right_obj = Box::new(build_bvh(hitable_list, &mut b, sorted_axis));
            (left_obj, right_obj)
        }
    };
    let left_box = left_obj
        .bounding_box()
        .expect("no bounding box in bvh_node constructor");
    let right_box = right_obj
        .bounding_box()
        .expect("no bounding box in bvh_node constructor");
    BvhNode {
        bvh_node_box: surrounding_box(left_box, right_box),
        left: left_obj,
        right: right_obj,
    }
}

impl BvhNode {
    pub fn new(hitable_list: &HitableList) -> Self {
        let mut handle = Vec::with_capacity(hitable_list.len());
        for i in 0..hitable_list.len() {
            handle.push(i);
        }

        dmerge_sort_wrap(&mut handle, box_x_compare, hitable_list);
        build_bvh(hitable_list, &mut handle, Axis::X)
    }
}

impl Hitable for BvhNode {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        match self.bvh_node_box.hit(r, t_min, t_max) {
            Some(_hit_rec) => match self.left.hit(r, t_min, t_max) {
                Some(left_rec) => match self.right.hit(r, t_min, t_max) {
                    Some(right_rec) => {
                        if left_rec.get_t() < right_rec.get_t() {
                            return Some(left_rec);
                        } else {
                            return Some(right_rec);
                        }
                    }
                    None => return Some(left_rec),
                },
                None => match self.right.hit(r, t_min, t_max) {
                    Some(right_rec) => return Some(right_rec),
                    None => return None,
                },
            },
            None => return None,
        }
    }

    fn bounding_box(&self) -> Option<Aabb> {
        Some(self.bvh_node_box.clone())
    }

    fn pdf_value(&self, o: &Vector3<f64>, v: &Vector3<f64>) -> f64 {
        return (self.left.pdf_value(o, v) + self.right.pdf_value(o, v)) / 2.0; // now bvh was fill both left and right so / 2.0
    }

    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        let mut rng = rand::thread_rng();
        let rand: f64 = rng.gen();

        if rand > 0.5 {
            return self.left.random(o);
        } else {
            return self.right.random(o);
        }
    }
}
