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
    last_node_num: usize,
    hitable_list_num: usize,
    nor_hitable_list_num: f64,
}

#[derive(Clone)]
pub struct BvhNode {
    bvh_node_box: Aabb,
    left: usize,
    right: usize,
    next_pos_diff: usize,// 1 is last // (2^depth) -1
}

fn dmerge(
    vec: &mut Vec<usize>,
    stock_vec: &mut Vec<usize>,
    compare_axis: usize,
    center_list: &Vec<Vector3<f64>>,
    left: usize,
    mid: usize,
    right: usize,
    ) {
    let mut i: usize = left;
    let mut j: usize = mid;
    let mut k: usize = 0;
    let length: usize = right - left;

    if i < mid && j < right {
        loop {
            if center_list[vec[i]][compare_axis] < center_list[vec[j]][compare_axis] {
                stock_vec[k] = vec[i];
                i = i + 1;
                if i == mid {
                    k = k + 1;
                    break;
                }
            } else {
                stock_vec[k] = vec[j];
                j = j + 1;
                if j == right {
                    k = k + 1;
                    break;
                }
            }
            k = k + 1;
        }
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
    compare_axis: usize,
    center_list: &Vec<Vector3<f64>>,
    ) {
    let len = vec.len();
    // stock_vec is temporary memory for merge sort.
    let mut stock_vec: Vec<usize> = Vec::with_capacity(len);
    stock_vec.resize_with(len, Default::default);

    for i in 0..(len/2) { // first merge two element use swap.
        let left = i*2;
        let right = left+1;
        if center_list[vec[right]][compare_axis] < center_list[vec[left]][compare_axis] {
            vec.swap(left, right);
        }
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
            dmerge(vec, &mut stock_vec, compare_axis, center_list, i, i+k, right);

            i = next_block;
        }
        k = k*2;
    }
}

const AI_X: usize = 0;
const AI_Y: usize = 1;
const AI_Z: usize = 2;

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
fn build_bvh(hitable_list: &HitableList, handle: &Vec<usize>, pre_sort_axis: &Axis, bvh_node_list: &mut Vec<BvhNode>, bvh_depth: usize, empty_hitable_handle: usize, center_list: &Vec<Vector3<f64>>) -> usize {
    let handle_size = handle.len();
    let next_pos_diff = 2usize.pow(bvh_depth as u32) - 1;
    match handle_size {
        1 => { // create bvh new node when item is onece,
               // but we have *not* perfect binary tree,
               // so need more checks depth in 2 => arm.
            //assert_eq!(bvh_depth, 1);
            let new_node = BvhNode {
                bvh_node_box: (hitable_list[handle[0]].bounding_box().unwrap()).clone(),
                left: handle[0],
                right: empty_hitable_handle,
                next_pos_diff,
            };
            let pos = bvh_node_list.len();
            bvh_node_list.push(new_node);
            return pos;
        },
        2 => {
            if bvh_depth == 1 {
                let new_node = BvhNode {
                    bvh_node_box: surrounding_box(hitable_list[handle[0]].bounding_box().unwrap()
                                                  , hitable_list[handle[1]].bounding_box().unwrap()),
                    left: handle[0],
                    right: handle[1],
                    next_pos_diff,
                };
                let pos = bvh_node_list.len();
                bvh_node_list.push(new_node);
                return pos;
            } else {
                //assert_eq!(bvh_depth, 2); // should bvh_depth 1 or 2
                
                let left_handle = build_bvh(hitable_list, &vec![handle[0]], &pre_sort_axis, bvh_node_list, bvh_depth - 1, empty_hitable_handle, center_list);
                let right_handle = build_bvh(hitable_list, &vec![handle[1]], &pre_sort_axis, bvh_node_list, bvh_depth - 1, empty_hitable_handle, center_list);
                let new_node = BvhNode {
                    bvh_node_box: surrounding_box(&bvh_node_list[left_handle].bvh_node_box
                                                  , &bvh_node_list[right_handle].bvh_node_box),
                    left: left_handle,
                    right: right_handle,
                    next_pos_diff,
                };
                let pos = bvh_node_list.len();
                bvh_node_list.push(new_node);
                return pos;
            };
        },
        _ => {
            let mut handle_x: Vec<usize> = handle.clone();
            let mut handle_y: Vec<usize> = handle.clone();
            let mut handle_z: Vec<usize> = handle.clone();
            match pre_sort_axis {
                Axis::X => {
                    dmerge_sort_wrap(&mut handle_y, AI_Y, center_list);
                    dmerge_sort_wrap(&mut handle_z, AI_Z, center_list);
                },
                Axis::Y => {
                    dmerge_sort_wrap(&mut handle_x, AI_X, center_list);
                    dmerge_sort_wrap(&mut handle_z, AI_Z, center_list);
                },
                Axis::Z => {
                    dmerge_sort_wrap(&mut handle_x, AI_X, center_list);
                    dmerge_sort_wrap(&mut handle_y, AI_Y, center_list);
                },
            }
            /*
            let x_max: f64 = center_list[handle_x[handle_size - 1]][0]
                - center_list[handle_x[0]][0];
            let y_max: f64 = center_list[handle_x[handle_size - 1]][1]
                - center_list[handle_x[0]][1];
            let z_max: f64 = center_list[handle_x[handle_size - 1]][2]
                - center_list[handle_x[0]][2];
            */
            let x_max: f64 = hitable_list[handle_x[handle_size - 1]]
                .bounding_box()
                .unwrap()
                .b_max[0]
                - hitable_list[handle_x[0]].bounding_box().unwrap().b_min[0];
            let y_max: f64 = hitable_list[handle_y[handle_size - 1]]
                .bounding_box()
                .unwrap()
                .b_max[1]
                - hitable_list[handle_y[0]].bounding_box().unwrap().b_min[1];
            let z_max: f64 = hitable_list[handle_z[handle_size - 1]]
                .bounding_box()
                .unwrap()
                .b_max[2]
                - hitable_list[handle_z[0]].bounding_box().unwrap().b_min[2];

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

            let left_handle = build_bvh(hitable_list, &a, &sorted_axis, bvh_node_list, bvh_depth - 1, empty_hitable_handle, center_list);
            let right_handle = build_bvh(hitable_list, &b, &sorted_axis, bvh_node_list, bvh_depth - 1, empty_hitable_handle, center_list);
            let new_node = BvhNode {
                bvh_node_box: surrounding_box(&bvh_node_list[left_handle].bvh_node_box
                                              , &bvh_node_list[right_handle].bvh_node_box),
                left: left_handle,
                right: right_handle,
                next_pos_diff,
            };
            let pos = bvh_node_list.len();
            bvh_node_list.push(new_node);
            return pos;
        },
    };
}

impl BvhTree {
    pub fn new(mut hitable_list: HitableList) -> Self {
        let hitable_list_len: usize = hitable_list.len();
        let mut handle = Vec::with_capacity(hitable_list_len);
        for i in 0..hitable_list_len {
            handle.push(i);
        }

        let mut aabb_center_list = Vec::with_capacity(hitable_list_len);
        for i in 0..hitable_list_len {
            let bounding_box_max = hitable_list[i].bounding_box().unwrap().b_max;
            let bounding_box_min = hitable_list[i].bounding_box().unwrap().b_min;
            let center_point: Vector3<f64> = [
                (bounding_box_max[0] + bounding_box_min[0])*0.5,
                (bounding_box_max[1] + bounding_box_min[1])*0.5,
                (bounding_box_max[2] + bounding_box_min[2])*0.5,
            ];
            aabb_center_list.push(center_point);
        }

        let hitable_list_next_power_of_two_len = hitable_list_len.next_power_of_two();
        let mut bvh_node_list: Vec<BvhNode> = Vec::with_capacity(hitable_list_next_power_of_two_len*2); // n:(0~k), sigma(2*n)
                                                                                                        // = (2*k) - 1
                                                    // and not affect bvh_node_list.len();
                                                    // affect only bvh_node_list.capacity();
                                                    
        bvh_node_list.push(BvhNode {
            bvh_node_box: Aabb{ b_min: [0.0, 0.0, 0.0], b_max: [0.0, 0.0, 0.0] },
            left: 0,
            right: 0,
            next_pos_diff: 0,
        }); // [0] dummy node; to actually node start at 1;
        dmerge_sort_wrap(&mut handle, AI_X, &aabb_center_list);

        let empty_hitable_handle: usize =  hitable_list_len;
        let bvh_tree_depth: usize = hitable_list_next_power_of_two_len.ilog2() as usize;
        let last_node_num = build_bvh(&hitable_list, &handle, &Axis::X, &mut bvh_node_list, bvh_tree_depth, empty_hitable_handle, &aabb_center_list);
        //println!("bvh_tree_depth: {}, last_node_num: {}", bvh_tree_depth, last_node_num);

        hitable_list.push(EmptyHitable::new()); // sould add after build_bvh// because break sort.

        let hitable_list_num = hitable_list.len() - 1;
        let nor_hitable_list_num = 1.0 / (hitable_list_num as f64);
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
            last_node_num,
            hitable_list_num,
            nor_hitable_list_num,
        }
    }
}

impl Hitable for BvhTree {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut current_pos: usize = self.last_node_num;
        let mut min_hit_t: f64 = t_max; //f64::MAX;
        let mut return_rec: Option<HitRecord> = None;
        loop {
            let current_bvh_node = &self.bvh_node_list[current_pos];
            let bvh_pos_diff = current_bvh_node.next_pos_diff;
            if bvh_pos_diff == 1 { // this node has actual item
                let right_obj = &self.hitable_list[current_bvh_node.right];
                if let Some(ref right_bounding_box) = right_obj.bounding_box() {
                    if let Some(_right_rec) = right_bounding_box.aabb_hit(r, t_min, min_hit_t) { // check bounding_box
                        if let Some(right_rec) = right_obj.hit(r, t_min, min_hit_t) { // actual hit check
                            let left_obj = &self.hitable_list[current_bvh_node.left];
                            if let Some(ref left_bounding_box) = left_obj.bounding_box() {
                                if let Some(_left_rec) = left_bounding_box.aabb_hit(r, t_min, right_rec.t) { // bounding_box
                                    if let Some(left_rec) = left_obj.hit(r, t_min, right_rec.t) { // acutual hit check
                                        let left_t = left_rec.t;
                                        return_rec = Some(left_rec);
                                        min_hit_t = left_t;
                                        current_pos = current_pos - 1;
                                        if current_pos == 0 {
                                            break; // no more hit node, ealy return;
                                        } else {
                                            continue;
                                        }
                                    }; // not hit left obj or right_t < left_t
                                }; // not hit left_bounding_box
                            }; // not found left bounding_box
                            let right_t = right_rec.t;
                            return_rec = Some(right_rec);
                            min_hit_t = right_t;
                            current_pos = current_pos - 1;
                            if current_pos == 0 {
                                break; // no more hit node, ealy return;
                            } else {
                                continue;
                            }
                        }; // not hit right obj
                    }; // not hit right_bounding_box
                }; // not found right_bounding_box
                let left_obj = &self.hitable_list[current_bvh_node.left];
                if let Some(ref left_bounding_box) = left_obj.bounding_box() {
                    if let Some(_left_rec) = left_bounding_box.aabb_hit(r, t_min, min_hit_t) { // bounding_box
                        if let Some(left_rec) = left_obj.hit(r, t_min, min_hit_t) { // acutual hit check
                            let left_t = left_rec.t;
                            return_rec = Some(left_rec);
                            min_hit_t = left_t;
                            current_pos = current_pos - 1;
                            if current_pos == 0 {
                                break; // no more hit node, ealy return;
                            } else {
                                continue;
                            }
                        };
                    };
                };
                // nothing update
                current_pos = current_pos - 1;
                if current_pos == 0 {
                    break; // no more hit node, ealy return;
                } else {
                    continue;
                }
            } else { // this node has other nodes
                match current_bvh_node.bvh_node_box.aabb_hit(r, t_min, min_hit_t) {
                    Some(_hit_rec) => { // if hit, next_pos_diff => 1;
                        current_pos = current_pos - 1;
                        // perfect tree so, not need check this case.
                        /*
                        if current_pos == 0 {
                            break; // no more hit node, ealy return;
                        } else {
                            continue;
                        }
                        */
                        continue;
                    },
                    None => {
                        // next_pos_diff: set skip number using current depth
                        // (2^depth) -1
                        // println!("(skip) next_pos_diff: {}, current_pos: {}, top_depth: {}", next_pos_diff, current_pos, top_depth);
                        current_pos = current_pos - bvh_pos_diff;
                        if current_pos == 0 {
                            break; // no more hit node, ealy return;
                        } else {
                            continue;
                        }
                   },
               };
            };
        }
        return return_rec
    }

    fn bounding_box<'a>(&'a self) -> Option<&'a Aabb> {
        Some(&self.aabb_box)
    }

    fn pdf_value(&self, ray: &Ray) -> f64 {
        if let Some(_aabb_hit) = self.aabb_box.aabb_hit(ray, 0.00001, 10000.0) {
            let hitable_list_len = self.hitable_list_num; // last is empty_hitable, so needs avoid
            let mut pdf_sum: f64 = 0.0;
            for i in 0..hitable_list_len {
                pdf_sum = pdf_sum + self.hitable_list[i].pdf_value(ray);
            }
            return pdf_sum * self.nor_hitable_list_num;
        } else {
            return 0.0;
        }
    }

    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        let hitable_list_len = self.hitable_list_num; // last is empty_hitable, so needs avoid
        let mut rng = rand::thread_rng();
        let rand: f64 = rng.gen();
        let rand_handle = (rand * hitable_list_len as f64) as usize;
        return self.hitable_list[rand_handle].random(o);
    }
}
