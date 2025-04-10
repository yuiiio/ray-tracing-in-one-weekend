use rand::prelude::*;

use crate::aabb::{surrounding_box, Aabb};
use crate::hitable::{HitRecord, Hitable};
use crate::hitablelist::HitableList;
use crate::quotation::Rotation;
use crate::ray::Ray;
use crate::vec3::Vector3;

#[derive(Clone)]
pub struct BvhTree {
    hitable_list: HitableList,
    bvh_node_list: Vec<BvhNode>,
    aabb_box: Aabb,
    last_node_num: usize,
    nor_hitable_list_num: f64,
}

#[derive(Clone)]
pub struct BvhNode {
    bvh_node_box: Aabb,
    left: usize,
    right: usize,
    skip_pos: usize,
    this_node_has_hitable: bool,
    only_have_left_obj: bool,
}

fn dmerge(
    vec: &mut [usize],
    stock_vec: &mut [usize],
    compare_axis: usize,
    center_list: &[Vector3<f64>],
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
                i += 1;
                if i == mid {
                    k += 1;
                    break;
                }
            } else {
                stock_vec[k] = vec[j];
                j += 1;
                if j == right {
                    k += 1;
                    break;
                }
            }
            k += 1;
        }
    }

    if i == mid {
        for m in k..length {
            stock_vec[m] = vec[j];
            j += 1;
        }
    } else {
        for m in k..length {
            stock_vec[m] = vec[i];
            i += 1;
        }
    }

    for m in 0..length {
        vec[left + m] = stock_vec[m];
    }
}

pub fn dmerge_sort_wrap(vec: &mut [usize], compare_axis: usize, center_list: &[Vector3<f64>]) {
    let len = vec.len();
    // stock_vec is temporary memory for merge sort.
    let mut stock_vec: Vec<usize> = Vec::with_capacity(len);
    stock_vec.resize_with(len, Default::default);

    for i in 0..(len / 2) {
        // first merge two element use swap.
        let left = i * 2;
        let right = left + 1;
        if center_list[vec[right]][compare_axis] < center_list[vec[left]][compare_axis] {
            vec.swap(left, right);
        }
    }

    let mut k = 2; //1 * 2;
    while k < len {
        // if len = 10, k => 1, 2, 4, 8
        let mut i = 0;
        while i < len {
            // k=1: i => 0, 2, 4, 6
            // k=2: i => 0, 4, 8, 12
            // k=4: i => 0, 8, 16, 24
            // k=8: i => 0, 16
            let next_block = i + (k * 2);
            // right: next_block: could over len, so need check and shrink to len
            let right = if len < next_block { len } else { next_block };
            dmerge(
                vec,
                &mut stock_vec,
                compare_axis,
                center_list,
                i,
                i + k,
                right,
            );

            i = next_block;
        }
        k *= 2;
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

// push bvh_node_list and return handle
//
// now allow non perfect balance tree.
//                          19
//          3                               18
//       1     2                14                       17
//                      10              13            15    16
//                  6       9        11    12
//                4   5   7   8
// skip_node is (last right top-parent's left brother node)

fn build_bvh(
    hitable_list: &HitableList,
    handle: &[usize],
    pre_sort_axis: &Axis,
    bvh_node_list: &mut Vec<BvhNode>,
    skip_pos: usize,
    center_list: &Vec<Vector3<f64>>,
) -> usize {
    let handle_size = handle.len();
    match handle_size {
        1 => {
            // create bvh new node when item is onece,
            // but we have *not* perfect binary tree,
            // so need more checks depth in 2 => arm.
            //assert_eq!(bvh_depth, 1);
            let new_node = BvhNode {
                bvh_node_box: (hitable_list[handle[0]].bounding_box()).clone(),
                left: handle[0],
                right: handle[0], // right == left, but should skip with flag
                skip_pos,
                this_node_has_hitable: true,
                only_have_left_obj: true,
            };
            let pos = bvh_node_list.len();
            bvh_node_list.push(new_node);
            pos
        }
        2 => {
            let new_node = BvhNode {
                bvh_node_box: surrounding_box(
                    hitable_list[handle[0]].bounding_box(),
                    hitable_list[handle[1]].bounding_box(),
                ),
                left: handle[0],
                right: handle[1],
                skip_pos,
                this_node_has_hitable: true,
                only_have_left_obj: false,
            };
            let pos = bvh_node_list.len();
            bvh_node_list.push(new_node);
            pos
        }
        _ => {
            let mut handle_2: Vec<usize> = handle.to_owned();
            let mut handle_3: Vec<usize> = handle.to_owned();
            let (handle_x, handle_y, handle_z): (&[usize], &[usize], &[usize]) = match pre_sort_axis
            {
                Axis::X => {
                    dmerge_sort_wrap(&mut handle_2, AI_Y, center_list);
                    dmerge_sort_wrap(&mut handle_3, AI_Z, center_list);
                    (handle, &handle_2, &handle_3)
                }
                Axis::Y => {
                    dmerge_sort_wrap(&mut handle_2, AI_X, center_list);
                    dmerge_sort_wrap(&mut handle_3, AI_Z, center_list);
                    (&handle_2, handle, &handle_3)
                }
                Axis::Z => {
                    dmerge_sort_wrap(&mut handle_2, AI_X, center_list);
                    dmerge_sort_wrap(&mut handle_3, AI_Y, center_list);
                    (&handle_2, &handle_3, handle)
                }
            };
            /*
            let x_max: f64 = center_list[handle_x[handle_size - 1]][0]
                - center_list[handle_x[0]][0];
            let y_max: f64 = center_list[handle_x[handle_size - 1]][1]
                - center_list[handle_x[0]][1];
            let z_max: f64 = center_list[handle_x[handle_size - 1]][2]
                - center_list[handle_x[0]][2];
            */
            let x_max: f64 = hitable_list[handle_x[handle_size - 1]].bounding_box().b_max[0]
                - hitable_list[handle_x[0]].bounding_box().b_min[0];
            let y_max: f64 = hitable_list[handle_y[handle_size - 1]].bounding_box().b_max[1]
                - hitable_list[handle_y[0]].bounding_box().b_min[1];
            let z_max: f64 = hitable_list[handle_z[handle_size - 1]].bounding_box().b_max[2]
                - hitable_list[handle_z[0]].bounding_box().b_min[2];

            let sorted_axis: Axis;
            let selected_handle = if x_max < y_max {
                if y_max < z_max {
                    sorted_axis = Axis::Z;
                    handle_z
                } else {
                    sorted_axis = Axis::Y;
                    handle_y
                }
            } else if x_max < z_max {
                sorted_axis = Axis::Z;
                handle_z
            } else {
                sorted_axis = Axis::X;
                handle_x
            };

            let (a, b) = selected_handle.split_at(handle_size / 2);

            let left_handle = build_bvh(
                hitable_list,
                a,
                &sorted_axis,
                bvh_node_list,
                skip_pos,
                center_list,
            );
            let right_handle = build_bvh(
                hitable_list,
                b,
                &sorted_axis,
                bvh_node_list,
                left_handle, // for next skip_pos
                center_list,
            );
            let new_node = BvhNode {
                bvh_node_box: surrounding_box(
                    &bvh_node_list[left_handle].bvh_node_box,
                    &bvh_node_list[right_handle].bvh_node_box,
                ),
                left: left_handle,
                right: right_handle,
                skip_pos,
                this_node_has_hitable: false,
                only_have_left_obj: false,
            };
            let pos = bvh_node_list.len();
            bvh_node_list.push(new_node);
            pos
        }
    }
}

impl BvhTree {
    pub fn new(hitable_list: HitableList) -> Self {
        let hitable_list_len: usize = hitable_list.len();
        let mut handle = Vec::with_capacity(hitable_list_len);
        for i in 0..hitable_list_len {
            handle.push(i);
        }

        let mut aabb_center_list = Vec::with_capacity(hitable_list_len);
        for i in 0..hitable_list_len {
            let bounding_box_max = hitable_list[i].bounding_box().b_max;
            let bounding_box_min = hitable_list[i].bounding_box().b_min;
            let center_point: Vector3<f64> = [
                (bounding_box_max[0] + bounding_box_min[0]) * 0.5,
                (bounding_box_max[1] + bounding_box_min[1]) * 0.5,
                (bounding_box_max[2] + bounding_box_min[2]) * 0.5,
            ];
            aabb_center_list.push(center_point);
        }

        let hitable_list_next_power_of_two_len = hitable_list_len.next_power_of_two();
        let mut bvh_node_list: Vec<BvhNode> =
            Vec::with_capacity(hitable_list_next_power_of_two_len * 2); // n:(0~k), sigma(2*n)
                                                                        // = (2*k) - 1
                                                                        // and not affect bvh_node_list.len();
                                                                        // affect only bvh_node_list.capacity();

        bvh_node_list.push(BvhNode {
            bvh_node_box: Aabb {
                b_min: [0.0, 0.0, 0.0],
                b_max: [0.0, 0.0, 0.0],
            },
            left: 0,
            right: 0,
            skip_pos: 0,
            this_node_has_hitable: false,
            only_have_left_obj: false,
        }); // [0] dummy node; to actually node start at 1;
        dmerge_sort_wrap(&mut handle, AI_X, &aabb_center_list);

        //let bvh_tree_depth: usize = hitable_list_next_power_of_two_len.ilog2() as usize;
        let last_node_num = build_bvh(
            &hitable_list,
            &handle,
            &Axis::X,
            &mut bvh_node_list,
            0,
            &aabb_center_list,
        );
        //println!("bvh_tree_depth: {}, last_node_num: {}", bvh_tree_depth, last_node_num);

        let nor_hitable_list_num = 1.0 / (hitable_list_len as f64);
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
            nor_hitable_list_num,
        }
    }
}

/*
 * future idea when use wide-bvh for simd
 * skip_pos_ref_list[0..last_node_num]
 * [0 1 2 3 4 5 6 7 8 9 10 ... last_node_num]
 *  child node [o * o *] qbvh *=>miss(doesn't hit node's aabb), o=>hit
 *             [5 5 6 8] => shift hit child node's skip_pos
 *                          and fill low side as (child node's minimum skip_pos-1)
 *
 *  current_pos = skip_pos_ref_list[current_bvh_node.skip_pos];
 *
 *  but every bvh-tree hit needs allocate [usize; last_node_num]
 *  ... might slow
 *  so this stack-recursive-less traversal is not work for wide-bvh for simd
*/

impl Hitable for BvhTree {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut current_pos: usize = self.last_node_num;
        let mut min_hit_t: f64 = t_max; //f64::MAX;
        let mut return_rec: Option<HitRecord> = None;
        let r_dir_inv = &r.get_inv_dir();
        loop {
            let current_bvh_node = &self.bvh_node_list[current_pos];
            if current_bvh_node.this_node_has_hitable == true {
                // this node has actual item
                if let Some(mut aabb_hit_rec) = current_bvh_node
                    .bvh_node_box
                    .aabb_hit(r, r_dir_inv, t_min, min_hit_t)
                {
                    let left_obj = &self.hitable_list[current_bvh_node.left];
                    if let Some(left_rec) = left_obj.hit(r, aabb_hit_rec.t_min, aabb_hit_rec.t_max)
                    {
                        aabb_hit_rec.t_max = left_rec.t; // for below right_obj hit
                        min_hit_t = left_rec.t;
                        return_rec = Some(left_rec);
                    };

                    if current_bvh_node.only_have_left_obj {
                        // not need check right
                    } else {
                        let right_obj = &self.hitable_list[current_bvh_node.right];
                        if let Some(right_rec) =
                            right_obj.hit(r, aabb_hit_rec.t_min, aabb_hit_rec.t_max)
                        // aabb_hit_rec.t_max = left_rec.t(if left hit), so when hit always right_rec.t < left_rec.t
                        {
                            min_hit_t = right_rec.t;
                            return_rec = Some(right_rec);
                        };
                    }
                }
            } else {
                // this node has other nodes
                if let Some(_hit_rec) = current_bvh_node
                    .bvh_node_box
                    .aabb_hit(r, r_dir_inv, t_min, min_hit_t)
                {
                    current_pos -= 1;
                    // not hitable having node's is not current_pos = 1, so not needs index == 0 check
                    continue;
                }
            }
            // if not hit bvh_node's aabb or check's after has_hitable node.
            current_pos = current_bvh_node.skip_pos;
            if current_pos == 0 {
                break; // no more hit node, ealy return;
            } else {
                continue;
            }
        }
        return_rec
    }

    fn bounding_box(&self) -> &Aabb {
        &self.aabb_box
    }

    fn bounding_box_with_rotate(&self, quat: &Rotation) -> Aabb {
        // TODO: depth(n) trival more small aabb by child node
        // now all same as hitablelist
        self.hitable_list.bounding_box_with_rotate(quat)
    }

    fn pdf_value(&self, ray: &Ray) -> f64 {
        if let Some(_aabb_hit) = self
            .aabb_box
            .aabb_hit(ray, &ray.get_inv_dir(), 0.00001, 10000.0)
        {
            let hitable_list_len = self.hitable_list.len();
            let mut pdf_sum: f64 = 0.0;
            for i in 0..hitable_list_len {
                pdf_sum += self.hitable_list[i].pdf_value(ray);
            }
            pdf_sum * self.nor_hitable_list_num
        } else {
            0.0
        }
    }

    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        let hitable_list_len = self.hitable_list.len();
        let mut rng = rand::thread_rng();
        let rand: f64 = rng.gen();
        let rand_handle = (rand * hitable_list_len as f64) as usize;
        self.hitable_list[rand_handle].random(o)
    }

    fn rotate_onb(&mut self, quat: &Rotation) -> () {
        self.hitable_list.rotate_onb(quat);
    }
}
