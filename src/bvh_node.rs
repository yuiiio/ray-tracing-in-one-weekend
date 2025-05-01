extern crate alloc;
use rand::prelude::*;

use crate::aabb::{aabb_hit_simd, surrounding_box, Aabb, QBoxes};
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
    max_stack_size: usize,
}

#[derive(Clone)]
pub struct BvhNode {
    //bvh_node_box: Aabb,
    qbox: QBoxes,
    child: [usize; 4],
    this_node_has_hitable: bool,
    max_childs: usize,
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

fn split_heuristic<'a>(
    hitable_list: &HitableList,
    handle: &'a [usize],
    handle_2: &'a mut [usize],
    handle_3: &'a mut [usize],
    pre_sort_axis: &Axis,
    center_list: &Vec<Vector3<f64>>,
) -> (&'a [usize], &'a [usize], Axis) {
    let handle_size = handle.len();
    /*
    let mut handle_2: Vec<usize> = handle.to_owned();
    let mut handle_3: Vec<usize> = handle.to_owned();
    */
    let (handle_x, handle_y, handle_z): (&[usize], &[usize], &[usize]) = match pre_sort_axis {
        Axis::X => {
            dmerge_sort_wrap(handle_2, AI_Y, center_list);
            dmerge_sort_wrap(handle_3, AI_Z, center_list);
            (handle, handle_2, handle_3)
        }
        Axis::Y => {
            dmerge_sort_wrap(handle_2, AI_X, center_list);
            dmerge_sort_wrap(handle_3, AI_Z, center_list);
            (handle_2, handle, handle_3)
        }
        Axis::Z => {
            dmerge_sort_wrap(handle_2, AI_X, center_list);
            dmerge_sort_wrap(handle_3, AI_Y, center_list);
            (handle_2, handle_3, handle)
        }
    };

    let x_max: f64 = hitable_list[handle_x[handle_size - 1]].bounding_box().b_max[0]
        - hitable_list[handle_x[0]].bounding_box().b_min[0];
    let y_max: f64 = hitable_list[handle_y[handle_size - 1]].bounding_box().b_max[1]
        - hitable_list[handle_y[0]].bounding_box().b_min[1];
    let z_max: f64 = hitable_list[handle_z[handle_size - 1]].bounding_box().b_max[2]
        - hitable_list[handle_z[0]].bounding_box().b_min[2];

    let return_sort_axis: Axis;
    let selected_handle = if x_max < y_max {
        if y_max < z_max {
            return_sort_axis = Axis::Z;
            handle_z
        } else {
            return_sort_axis = Axis::Y;
            handle_y
        }
    } else if x_max < z_max {
        return_sort_axis = Axis::Z;
        handle_z
    } else {
        return_sort_axis = Axis::X;
        handle_x
    };

    let (handle_a, handle_b) = selected_handle.split_at(handle_size / 2);
    (handle_a, handle_b, return_sort_axis)
}

fn aabb_size(aabb: &Aabb) -> f64 {
    (aabb.b_max[0] - aabb.b_min[0])
        * (aabb.b_max[1] - aabb.b_min[1])
        * (aabb.b_max[2] - aabb.b_min[2])
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

// at traversal use return_rec(stack), push(node_pos) if aabb' hit(without most right node), and pop() if last_node or not aabb'hit.

fn build_bvh(
    hitable_list: &HitableList,
    handle: &[usize],
    pre_sort_axis: &Axis,
    bvh_node_list: &mut Vec<BvhNode>,
    center_list: &Vec<Vector3<f64>>,
) -> (usize, Aabb) {
    let handle_size = handle.len();
    match handle_size {
        1 => {
            let new_node = BvhNode {
                qbox: QBoxes::encode_from_four_aabb(&[
                    hitable_list[handle[0]].bounding_box(),
                    hitable_list[handle[0]].bounding_box(), //
                    hitable_list[handle[0]].bounding_box(), //
                    hitable_list[handle[0]].bounding_box(), // not important
                ]),
                child: [handle[0], 0, 0, 0],
                this_node_has_hitable: true,
                max_childs: 1,
            };
            let pos = bvh_node_list.len();
            bvh_node_list.push(new_node);
            (pos, (hitable_list[handle[0]].bounding_box()).clone())
        }
        2 => {
            let new_node = BvhNode {
                qbox: QBoxes::encode_from_four_aabb(&[
                    hitable_list[handle[0]].bounding_box(),
                    hitable_list[handle[1]].bounding_box(), //
                    hitable_list[handle[1]].bounding_box(), //
                    hitable_list[handle[1]].bounding_box(), // not important
                ]),
                child: [handle[0], handle[1], 0, 0],
                this_node_has_hitable: true,
                max_childs: 2,
            };
            let pos = bvh_node_list.len();
            bvh_node_list.push(new_node);
            (
                pos,
                surrounding_box(
                    hitable_list[handle[0]].bounding_box(),
                    hitable_list[handle[1]].bounding_box(),
                ),
            )
        }
        3 => {
            let new_node = BvhNode {
                qbox: QBoxes::encode_from_four_aabb(&[
                    hitable_list[handle[0]].bounding_box(),
                    hitable_list[handle[1]].bounding_box(),
                    hitable_list[handle[2]].bounding_box(),
                    hitable_list[handle[2]].bounding_box(), // not important
                ]),
                child: [handle[0], handle[1], handle[2], 0],
                this_node_has_hitable: true,
                max_childs: 3,
            };
            let pos = bvh_node_list.len();
            bvh_node_list.push(new_node);
            (
                pos,
                surrounding_box(
                    &surrounding_box(
                        hitable_list[handle[0]].bounding_box(),
                        hitable_list[handle[1]].bounding_box(),
                    ),
                    hitable_list[handle[2]].bounding_box(),
                ),
            )
        }
        4 => {
            let new_node = BvhNode {
                qbox: QBoxes::encode_from_four_aabb(&[
                    hitable_list[handle[0]].bounding_box(),
                    hitable_list[handle[1]].bounding_box(),
                    hitable_list[handle[2]].bounding_box(),
                    hitable_list[handle[3]].bounding_box(),
                ]),
                child: [handle[0], handle[1], handle[2], handle[3]],
                this_node_has_hitable: true,
                max_childs: 4,
            };
            let pos = bvh_node_list.len();
            bvh_node_list.push(new_node);
            (
                pos,
                surrounding_box(
                    &surrounding_box(
                        hitable_list[handle[0]].bounding_box(),
                        hitable_list[handle[1]].bounding_box(),
                    ),
                    &surrounding_box(
                        hitable_list[handle[2]].bounding_box(),
                        hitable_list[handle[3]].bounding_box(),
                    ),
                ),
            )
        }
        _ => {
            let mut handle_2: Vec<usize> = handle.to_owned();
            let mut handle_3: Vec<usize> = handle.to_owned();
            let (handle_a, handle_b, pre_sort_axis) = split_heuristic(
                hitable_list,
                handle,
                &mut handle_2,
                &mut handle_3,
                pre_sort_axis,
                center_list,
            );

            // handle_a split
            let mut handle_2: Vec<usize> = handle_a.to_owned();
            let mut handle_3: Vec<usize> = handle_a.to_owned();
            let (handle_c, handle_d, pre_sort_axis_a) = split_heuristic(
                hitable_list,
                handle_a,
                &mut handle_2,
                &mut handle_3,
                &pre_sort_axis,
                center_list,
            );

            // handle_b split
            let mut handle_2: Vec<usize> = handle_b.to_owned();
            let mut handle_3: Vec<usize> = handle_b.to_owned();
            let (handle_e, handle_f, pre_sort_axis_b) = split_heuristic(
                hitable_list,
                handle_b,
                &mut handle_2,
                &mut handle_3,
                &pre_sort_axis,
                center_list,
            );

            let (first_handle, first_aabb) = build_bvh(
                hitable_list,
                handle_c,
                &pre_sort_axis_a,
                bvh_node_list,
                center_list,
            );
            let (second_handle, second_aabb) = build_bvh(
                hitable_list,
                handle_d,
                &pre_sort_axis_a,
                bvh_node_list,
                center_list,
            );

            let (third_handle, third_aabb) = build_bvh(
                hitable_list,
                handle_e,
                &pre_sort_axis_b,
                bvh_node_list,
                center_list,
            );
            let (fourth_handle, fourth_aabb) = build_bvh(
                hitable_list,
                handle_f,
                &pre_sort_axis_b,
                bvh_node_list,
                center_list,
            );

            //TODO: sort child node by box-size
            //  for fast traversal.
            let aabb_sizes = [
                aabb_size(&first_aabb),
                aabb_size(&second_aabb),
                aabb_size(&third_aabb),
                aabb_size(&fourth_aabb),
            ];
            let mut sort_handle = [0, 1, 2, 3];
            if aabb_sizes[0] < aabb_sizes[1] {
                sort_handle[0] = 1;
                sort_handle[1] = 0;
            }
            if aabb_sizes[2] < aabb_sizes[3] {
                sort_handle[2] = 3;
                sort_handle[3] = 2;
            }
            if aabb_sizes[sort_handle[0]] < aabb_sizes[sort_handle[2]] {
                let tmp0 = sort_handle[0];
                sort_handle[0] = sort_handle[2];
                if aabb_sizes[tmp0] < aabb_sizes[sort_handle[3]] {
                    let tmp1 = sort_handle[1];
                    sort_handle[1] = sort_handle[3];
                    sort_handle[2] = tmp0;
                    sort_handle[3] = tmp1;
                } else {
                    let tmp1 = sort_handle[1];
                    sort_handle[1] = tmp0;
                    if aabb_sizes[tmp1] < aabb_sizes[sort_handle[3]] {
                        sort_handle[2] = sort_handle[3];
                        sort_handle[3] = tmp1;
                    } else {
                        sort_handle[2] = tmp1;
                        //sort_handle[3] = sort_handle[3];
                    }
                }
            } else {
                //sort_handle[0] = sort_handle[0];
                let tmp1 = sort_handle[1];
                if aabb_sizes[tmp1] < aabb_sizes[sort_handle[2]] {
                    sort_handle[1] = sort_handle[2];
                    if aabb_sizes[tmp1] < aabb_sizes[sort_handle[3]] {
                        sort_handle[2] = sort_handle[3];
                        sort_handle[3] = tmp1;
                    } else {
                        sort_handle[2] = tmp1;
                        //sort_handle[3] = sort_handle[3];
                    }
                } else {
                    //sort_handle[1] = tmp1;
                    //sort_handle[2] = sort_handle[2];
                    //sort_handle[3] = sort_handle[3];
                }
            }

            let handles = [first_handle, second_handle, third_handle, fourth_handle];
            let handle_aabbs = [&first_aabb, &second_aabb, &third_aabb, &fourth_aabb];

            let new_node = BvhNode {
                qbox: QBoxes::encode_from_four_aabb(&[
                    handle_aabbs[sort_handle[3]],
                    handle_aabbs[sort_handle[2]],
                    handle_aabbs[sort_handle[1]],
                    handle_aabbs[sort_handle[0]],
                ]),
                child: [
                    handles[sort_handle[3]],
                    handles[sort_handle[2]],
                    handles[sort_handle[1]],
                    handles[sort_handle[0]],
                ],
                this_node_has_hitable: false,
                max_childs: 4,
            };
            let pos = bvh_node_list.len();
            bvh_node_list.push(new_node);
            (
                pos,
                surrounding_box(
                    &surrounding_box(&first_aabb, &second_aabb),
                    &surrounding_box(&third_aabb, &fourth_aabb),
                ),
            )
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

        let mut bvh_node_list: Vec<BvhNode> = Vec::with_capacity((hitable_list_len - 1) / (4 - 1)); // qbvh min node size

        bvh_node_list.push(BvhNode {
            qbox: QBoxes {
                bboxes: [[[0.0; 4]; 3]; 2],
            },
            child: [0, 0, 0, 0],
            this_node_has_hitable: false,
            max_childs: 4,
        }); // [0] dummy node; to actually node start at 1;
        dmerge_sort_wrap(&mut handle, AI_X, &aabb_center_list);

        let bvh_tree_depth: usize = (hitable_list_len.next_power_of_two() * 2).ilog(4) as usize;
        let max_stack_size = bvh_tree_depth * (4 - 1);
        let (last_node_num, bvh_aabb) = build_bvh(
            &hitable_list,
            &handle,
            &Axis::X,
            &mut bvh_node_list,
            &aabb_center_list,
        );
        //println!("bvh_tree_depth: {}, last_node_num: {}", bvh_tree_depth, last_node_num);

        let nor_hitable_list_num = 1.0 / (hitable_list_len as f64);

        BvhTree {
            hitable_list,
            bvh_node_list,
            aabb_box: bvh_aabb,
            last_node_num,
            nor_hitable_list_num,
            max_stack_size,
        }
    }
}

impl Hitable for BvhTree {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut current_pos: usize = self.last_node_num;
        let mut min_hit_t: f64 = t_max; //f64::MAX;
        let mut return_rec: Option<HitRecord> = None;
        let r_dir_inv = &r.get_inv_dir();
        let mut return_stack: alloc::vec::Vec<usize> =
            alloc::vec::Vec::with_capacity(self.max_stack_size);
        return_stack.push(0);
        loop {
            let current_bvh_node = &self.bvh_node_list[current_pos];
            let aabb_result = aabb_hit_simd(
                &current_bvh_node.qbox,
                &r.origin,
                r_dir_inv,
                t_min,
                min_hit_t,
            );
            if current_bvh_node.this_node_has_hitable == true {
                for i in 0..current_bvh_node.max_childs {
                    if aabb_result[i] == true {
                        let child_obj = &self.hitable_list[current_bvh_node.child[i]];
                        if let Some(child_rec) = child_obj.hit(r, t_min, min_hit_t) {
                            min_hit_t = child_rec.t;
                            return_rec = Some(child_rec);
                        }
                    }
                }
            } else {
                // this node has other nodes
                if aabb_result[0] == true {
                    return_stack.push(current_bvh_node.child[0]);
                }
                if aabb_result[1] == true {
                    return_stack.push(current_bvh_node.child[1]);
                }
                if aabb_result[2] == true {
                    return_stack.push(current_bvh_node.child[2]);
                }

                if aabb_result[3] == true {
                    current_pos = current_bvh_node.child[3];
                    continue;
                }
            }
            // if not hit bvh_node's aabb or check's after has_hitable node.
            current_pos = return_stack.pop().unwrap();
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
        if let Some(_aabb_hit) =
            self.aabb_box
                .aabb_hit(&ray.origin, &ray.get_inv_dir(), 0.00001, 10000.0)
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
