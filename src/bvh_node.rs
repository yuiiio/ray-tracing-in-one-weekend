use rand::prelude::*;

use crate::aabb::{surrounding_box, Aabb};
use crate::hitable::{HitRecord, Hitable};
use crate::hitablelist::HitableList;
use crate::ray::Ray;
use crate::vec3::Vector3;

#[derive(Clone)]
pub struct BvhTree {
    hitable_list: HitableList,
    bvh_node_list: Vec<BvhNode>,
    aabb_box: Aabb,
}

#[derive(Clone)]
pub struct BvhNode {
    bvh_node_box: Aabb,
    left: usize,
    right: usize,
    bvh_depth: usize,// 1 is last
}

fn box_x_compare(a: &Box<dyn Hitable + Send + Sync>, b: &Box<dyn Hitable + Send + Sync>) -> bool {
    let box_a: &Aabb = a.bounding_box().unwrap();
    let box_b: &Aabb = b.bounding_box().unwrap();
    if box_a.b_min()[0] < box_b.b_min()[0] {
        return true;
    } else {
        return false;
    }
}

fn box_y_compare(a: &Box<dyn Hitable + Send + Sync>, b: &Box<dyn Hitable + Send + Sync>) -> bool {
    let box_a: &Aabb = a.bounding_box().unwrap();
    let box_b: &Aabb = b.bounding_box().unwrap();
    if box_a.b_min()[1] < box_b.b_min()[1] {
        return true;
    } else {
        return false;
    }
}

fn box_z_compare(a: &Box<dyn Hitable + Send + Sync>, b: &Box<dyn Hitable + Send + Sync>) -> bool {
    let box_a: &Aabb = a.bounding_box().unwrap();
    let box_b: &Aabb = b.bounding_box().unwrap();
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
    let length: usize = right - left;

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
        for m in k..length {
            stock_vec[m] = vec[j];
            j = j + 1;
        }
    } else {
        for m in k..length {
            stock_vec[m] = vec[i];
            i = i + 1;
        }
    }

    for m in 0..length {
        vec[left + m] = stock_vec[m];
    }
}

pub fn dmerge_sort_wrap(
    vec: &mut Vec<usize>,
    compare: fn(&Box<dyn Hitable + Send + Sync>, &Box<dyn Hitable + Send + Sync>) -> bool,
    hitable_list: &HitableList,
    ) {
    let len = vec.len();
    // stock_vec is temporary memory for merge sort.
    let mut stock_vec: Vec<usize> = Vec::with_capacity(len);
    stock_vec.resize_with(len, Default::default);

    let mut i = 0;
    while i < len {
        let next_block = i + 2;
        let right = i+1;
        if len <= right {
            break; 
        } else {
            if compare(&hitable_list[vec[right]], &hitable_list[vec[i]]) {
                vec.swap(i, right);
            }
        }
        i = next_block;
    }

    let mut k = 1 * 2;
    while k < len {// if len = 10, k => 1, 2, 4, 8
        let mut i = 0;
        while i < len { // k=1: i => 0, 2, 4, 6
                        // k=2: i => 0, 4, 8, 12
                        // k=4: i => 0, 8, 16, 24
                        // k=8: i => 0, 16
            let next_block = i + (k*2);
            // right: next_block: could over len, so need check and shrink to len
            let right = if len < next_block { len } else { next_block };
            dmerge(vec, &mut stock_vec, compare, hitable_list, i, i+k, right);

            i = next_block;
        }
        k = k*2;
    }
}

enum Axis {
    X,
    Y,
    Z,
}

#[derive(Clone)]
pub struct EmptyHitable {
}
impl EmptyHitable{
    pub fn new() -> Self {
        EmptyHitable{}
    }
}
impl Hitable for EmptyHitable {
    fn hit(&self, _r: &Ray, _t_min: f64, _t_max: f64) -> Option<HitRecord> {
        None
    }
    fn bounding_box<'a>(&'a self) -> Option<&'a Aabb> {
        None
    }
}


// push bvh_node_list and return handle
//                    15
//              7             14       <=  diff 7 (2^3 - 1)
//           3     6      10      13   <=  diff 3 (2^2 - 1)
//          1 2   4 5    8  9    11 12 <=  diff 1 (2^1 - 1)
fn build_bvh(hitable_list: &HitableList, handle: &Vec<usize>, pre_sort_axis: &Axis, bvh_node_list: &mut Vec<BvhNode>, bvh_depth: usize) -> usize {
    let handle_size = handle.len();
    match handle_size {
        1 => { // create bvh new node when item is onece,
               // but we have *not* perfect binary tree,
               // so need more checks depth in 2 => arm.
            //assert_eq!(bvh_depth, 1);
            let new_node = BvhNode {
                bvh_node_box: (hitable_list[handle[0]].bounding_box().unwrap()).clone(),
                left: handle[0],
                right: handle[0], // maybe should have empty hitable
                bvh_depth,
            };
            bvh_node_list.push(new_node);
            return bvh_node_list.len() - 1;
        },
        2 => {
            if bvh_depth == 1 {
                let new_node = BvhNode {
                    bvh_node_box: surrounding_box(hitable_list[handle[0]].bounding_box().unwrap()
                                                  , hitable_list[handle[1]].bounding_box().unwrap()),
                    left: handle[0],
                    right: handle[1],
                    bvh_depth,
                };
                bvh_node_list.push(new_node);
                return bvh_node_list.len() - 1;
            } else {
                //assert_eq!(bvh_depth, 2); // should bvh_depth 1 or 2
                
                let left_handle = build_bvh(hitable_list, &vec![handle[0]], &pre_sort_axis, bvh_node_list, bvh_depth - 1);
                let right_handle = build_bvh(hitable_list, &vec![handle[1]], &pre_sort_axis, bvh_node_list, bvh_depth - 1);
                let new_node = BvhNode {
                    bvh_node_box: surrounding_box(&bvh_node_list[left_handle].bvh_node_box
                                                  , &bvh_node_list[right_handle].bvh_node_box),
                    left: left_handle,
                    right: right_handle,
                    bvh_depth,
                };
                bvh_node_list.push(new_node);
                return bvh_node_list.len() - 1;
            };
        },
        _ => {
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
            let x_max: f64 = hitable_list[handle_x[handle_size - 1]]
                .bounding_box()
                .unwrap()
                .b_max()[0]
                - hitable_list[handle_x[0]].bounding_box().unwrap().b_min()[0];
            let y_max: f64 = hitable_list[handle_y[handle_size - 1]]
                .bounding_box()
                .unwrap()
                .b_max()[1]
                - hitable_list[handle_y[0]].bounding_box().unwrap().b_min()[1];
            let z_max: f64 = hitable_list[handle_z[handle_size - 1]]
                .bounding_box()
                .unwrap()
                .b_max()[2]
                - hitable_list[handle_z[0]].bounding_box().unwrap().b_min()[2];

            let sorted_axis: Axis;
            let mut selected_handle = if x_max < y_max {
                if y_max < z_max {
                    sorted_axis = Axis::Z;
                    handle_z
                } else {
                    sorted_axis = Axis::Y;
                    handle_y
                }
            } else {
                if x_max < z_max {
                    sorted_axis = Axis::Z;
                    handle_z
                } else {
                    sorted_axis = Axis::X;
                    handle_x
                }
            };
            let a = selected_handle.split_off(handle_size / 2);
            let b = selected_handle;

            let left_handle = build_bvh(hitable_list, &a, &sorted_axis, bvh_node_list, bvh_depth - 1);
            let right_handle = build_bvh(hitable_list, &b, &sorted_axis, bvh_node_list, bvh_depth - 1);
            let new_node = BvhNode {
                bvh_node_box: surrounding_box(&bvh_node_list[left_handle].bvh_node_box
                                              , &bvh_node_list[right_handle].bvh_node_box),
                left: left_handle,
                right: right_handle,
                bvh_depth,
            };
            bvh_node_list.push(new_node);
            return bvh_node_list.len() - 1;
        },
    };
}

impl BvhTree {
    pub fn new(hitable_list: HitableList) -> Self {
        let hitable_list_len: usize = hitable_list.len();
        let mut handle = Vec::with_capacity(hitable_list_len);
        for i in 0..hitable_list_len {
            handle.push(i);
        }

        let mut bvh_node_list: Vec<BvhNode> = Vec::new();
        bvh_node_list.push(BvhNode {
            bvh_node_box: Aabb::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
            left: 0,
            right: 0,
            bvh_depth: 0,
        }); // [0] dummy node; to actually node start at 1;
        dmerge_sort_wrap(&mut handle, box_x_compare, &hitable_list);

        let bvh_tree_depth: usize = hitable_list_len.next_power_of_two().ilog2() as usize;
        let last_node_num = build_bvh(&hitable_list, &handle, &Axis::X, &mut bvh_node_list, bvh_tree_depth);
        //println!("bvh_tree_depth: {}, last_node_num: {}", bvh_tree_depth, last_node_num);

        /*
        let mut k = 1;
        for now_depth in 0..bvh_tree_depth {
            for i in 0..k {
            }
            k = k*2;
        }
        */
        let aabb_box = bvh_node_list[last_node_num].bvh_node_box.clone();
        BvhTree {
            hitable_list,
            bvh_node_list,
            aabb_box,
        }
    }
}

impl Hitable for BvhTree {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let last_node_num: usize = self.bvh_node_list.len() - 1;
        let mut current_pos: usize = last_node_num;
        let mut min_hit_t: f64 = f64::MAX;
        let mut return_rec: Option<HitRecord> = None;
        //let top_depth = self.bvh_node_list[last_node_num].bvh_depth;
        loop {
            let mut next_pos_diff: usize = 1;
            let current_bvh_node = &self.bvh_node_list[current_pos];
            let current_depth = current_bvh_node.bvh_depth;
            if current_depth == 1 { // this node has actual item
                let right_obj = &self.hitable_list[current_bvh_node.right];
                match right_obj.bounding_box().unwrap().hit(r, t_min, t_max) { // check bounding_box
                    Some(_right_rec) => {
                        match right_obj.hit(r, t_min, t_max) { // actual hit check
                            Some(right_rec) => {
                                let left_obj = &self.hitable_list[current_bvh_node.left];
                                match left_obj.bounding_box().unwrap().hit(r, t_min, t_max) { // bounding_box
                                    Some(_left_rec) => {
                                        match left_obj.hit(r, t_min, t_max) { // acutual hit check
                                            Some(left_rec) => {
                                                let left_t = left_rec.get_t();
                                                let right_t = right_rec.get_t();
                                                if left_rec.get_t() < right_t {
                                                    if left_t < min_hit_t {
                                                        return_rec = Some(left_rec);
                                                        min_hit_t = left_t;
                                                    };
                                                } else {
                                                    if right_t < min_hit_t {
                                                        return_rec = Some(right_rec);
                                                        min_hit_t = right_t;
                                                    };
                                                };
                                            },
                                            None => {
                                                let right_t = right_rec.get_t();
                                                if right_t < min_hit_t {
                                                    return_rec = Some(right_rec);
                                                    min_hit_t = right_t;
                                                };
                                            },
                                        };
                                    },
                                    None => {
                                        let right_t = right_rec.get_t();
                                        if right_t < min_hit_t {
                                            return_rec = Some(right_rec);
                                            min_hit_t = right_t;
                                        };
                                    },
                                };
                            },
                            None => {
                                let left_obj = &self.hitable_list[current_bvh_node.left];
                                match left_obj.bounding_box().unwrap().hit(r, t_min, t_max) { // bounding_box
                                    Some(_left_rec) => {
                                        match left_obj.hit(r, t_min, t_max) { // acutual hit check
                                            Some(left_rec) => {
                                                let left_t = left_rec.get_t();
                                                if left_t < min_hit_t {
                                                    return_rec = Some(left_rec);
                                                    min_hit_t = left_t;
                                                };
                                            },
                                            None => { // nothing update
                                            },
                                        };
                                    },
                                    None => { // nothing update
                                    },
                                };
                            },
                        };
                    },
                    None => {
                        let left_obj = &self.hitable_list[current_bvh_node.left];
                        match left_obj.bounding_box().unwrap().hit(r, t_min, t_max) { // bounding_box
                            Some(_left_rec) => {
                                match left_obj.hit(r, t_min, t_max) { // acutual hit check
                                    Some(left_rec) => {
                                        let left_t = left_rec.get_t();
                                        if left_t < min_hit_t {
                                            return_rec = Some(left_rec);
                                            min_hit_t = left_t;
                                        };
                                    },
                                    None => { // nothing update
                                    },
                                };
                            },
                            None => { // nothing update
                            },
                        };
                    },
                };
                if current_pos == 0 { // *1
                    break; // last node
                }
            } else { // this node has other nodes
               match current_bvh_node.bvh_node_box.hit(r, t_min, t_max) {
                   Some(_hit_rec) => {},
                   None => {
                       next_pos_diff = 2usize.pow(current_depth as u32) - 1; // set skip number using current depth
                                                                             // (2^depth) -1
                        //println!("(skip) next_pos_diff: {}, current_pos: {}, top_depth: {}", next_pos_diff, current_pos, top_depth);
                   },
               };
            };
            //assert!(current_pos >= next_pos_diff);
            current_pos = current_pos - next_pos_diff; // 1 is always
                                                       // this_node_is_last(line)
                                                       // *1 so sould not be negative
            if current_pos == 0 {
                break; // no more hit node, ealy return;
            }
        }
        return return_rec
    }

    fn bounding_box<'a>(&'a self) -> Option<&'a Aabb> {
        Some(&self.aabb_box)
    }

    fn pdf_value(&self, o: &Vector3<f64>, v: &Vector3<f64>) -> f64 {
        let hitable_list_len = self.hitable_list.len();
        let mut pdf_sum: f64 = 0.0;
        for i in 0..hitable_list_len {
            pdf_sum = pdf_sum + self.hitable_list[i].pdf_value(o, v);
        }
        return pdf_sum / (hitable_list_len as f64);
    }

    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        let hitable_list_len = self.hitable_list.len();
        let mut rng = rand::thread_rng();
        let rand: f64 = rng.gen();
        let rand_handle = (rand * hitable_list_len as f64) as usize;
        return self.hitable_list[rand_handle].random(o);
    }
}
