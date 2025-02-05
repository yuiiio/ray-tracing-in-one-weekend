#[allow(dead_code)]
pub fn clamp(num: f64, min: f64, max: f64) -> f64 {
    if num < min {
        return min;
    }
    if max < num {
        return max;
    }
    num
}

pub fn max(a: f64, b: f64) -> f64 {
    if a < b {
        b
    } else {
        a
    }
}

pub fn min(a: f64, b: f64) -> f64 {
    if a > b {
        b
    } else {
        a
    }
}

// merge
// left, mid, right
// 2, 3, 1
//
#[allow(dead_code)]
fn merge<T: Clone>(
    vec: &mut [T],
    stock_vec: &mut [T],
    compare: fn(&T, &T) -> bool,
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
            if compare(&vec[i], &vec[j]) {
                stock_vec[k] = vec[i].clone();
                i += 1;
                if i == mid {
                    k += 1;
                    break;
                }
            } else {
                stock_vec[k] = vec[j].clone();
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
        // already left~mid is move to stock_vec
        // so mid(>j)~right move to stock_vec
        for m in k..length {
            stock_vec[m] = vec[j].clone();
            j += 1;
        }
    } else {
        for m in k..length {
            stock_vec[m] = vec[i].clone();
            i += 1;
        }
    }

    for m in 0..length {
        vec[left + m] = stock_vec[m].clone();
    }
}

#[allow(dead_code)]
pub fn merge_sort<T: Clone + std::fmt::Debug>(
    vec: &mut [T],
    stock_vec: &mut [T],
    compare: fn(&T, &T) -> bool,
    left: usize,
    right: usize,
) {
    if (left == right) || (left + 1 == right) {
        // when left+1 = right,
        // left = 1,right = 3: mid = (1+3)/2; => 2;
        // {[1, 2], [2, 3]} => return
        // merge (1, 2, 3) => mid 2 is absolutely compare
        // vs 1
        println!("ealy return: left:{}, right:{}\n", left, right);
        return;
    }
    let mid: usize = (left + right) / 2; // (right - left)/2 + left
    merge_sort(vec, stock_vec, compare, left, mid);
    println!("left-mid{:?}\n", vec);
    merge_sort(vec, stock_vec, compare, mid, right);
    println!("mid-right{:?}\n", vec);
    merge(vec, stock_vec, compare, left, mid, right);
    println!("{:?}\n", vec);
}

mod test {
    #![allow(unused_imports)]
    use super::*;

    #[test]
    fn usize_test() {
        let a: usize = (9 + 2) / 2;
        assert_eq!(a, 5);
    }

    #[test]
    fn merge_sort_test() {
        let mut vec = vec![4, 5, 2, 6, 1, 8, 3, 5, 7];
        fn compare(a: &i32, b: &i32) -> bool {
            if a < b {
                return true;
            } else {
                return false;
            }
        }

        let len = vec.len();

        let mut stock_vec: Vec<i32> = Vec::with_capacity(len);
        stock_vec.resize_with(len, Default::default);

        //merge_sort(&mut vec, &mut stock_vec, compare, 0, len);

        // if vec.len() => 9. (vec index [0, 1, 2, 3, 4, 5, 6, 7, 8])
        // 9/2 as usize => 4,
        for i in 0..(len / 2) {
            // first merge two element use swap.
            let left = i * 2;
            let right = left + 1;
            if vec[left] > vec[right] {
                println!("swap: left:{}, right:{}\n", left, right);
                vec.swap(left, right);
            }
        }

        // 2, 2, 2, 2, 1
        // 4(left:0, mid: 2, right: 4), 4(left:4, mid: 6, right: 8), 1(left: 8, mid: 10(over), right: 9(len))
        // 8(left 0, mid: 4, right: 8), 1(left: 8, mid: 12(over), right: 9(len))
        // like
        // or
        // 2, 2, 2, 2, 2
        // 4, 4, 2(left: 8, mid: 10(over), 10(len))
        // 8, 2

        let mut k = 1 * 2;
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
                merge(&mut vec, &mut stock_vec, compare, i, i + k, right);

                i = next_block;
            }
            k = k * 2;
        }

        assert_eq!(vec, [1, 2, 3, 4, 5, 5, 6, 7, 8]);
    }
}
