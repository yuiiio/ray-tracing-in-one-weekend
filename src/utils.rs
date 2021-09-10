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

fn merge<T: Clone>(
    vec: &mut Vec<T>,
    stock_vec: &mut Vec<T>,
    compare: fn(&T, &T) -> bool,
    left: usize,
    mid: usize,
    right: usize,
) {
    let mut i: usize = left;
    let mut j: usize = mid;
    let mut k: usize = 0;
    let mut l: usize = 0;

    while i < mid && j < right {
        if compare(&vec[i], &vec[j]) {
            stock_vec[k] = vec[i].clone();
            i = i + 1;
        } else {
            stock_vec[k] = vec[j].clone();
            j = j + 1;
        }
        k = k + 1;
    }

    if i == mid {
        while j < right {
            stock_vec[k] = vec[j].clone();
            j = j + 1;
            k = k + 1;
        }
    } else {
        while i < mid {
            stock_vec[k] = vec[i].clone();
            i = i + 1;
            k = k + 1;
        }
    }

    while l < k {
        vec[left + l] = stock_vec[l].clone();
        l = l + 1;
    }
}

pub fn merge_sort<T: Clone>(
    vec: &mut Vec<T>,
    stock_vec: &mut Vec<T>,
    compare: fn(&T, &T) -> bool,
    left: usize,
    right: usize,
) {
    if (left == right) || (left + 1 == right) {
        return;
    }
    let mid: usize = (left + right) / 2;
    merge_sort(vec, stock_vec, compare, left, mid);
    merge_sort(vec, stock_vec, compare, mid, right);
    merge(vec, stock_vec, compare, left, mid, right);
}

mod test {
    use super::*;

    #[test]
    fn usize_test() {
        let a: usize = (9 + 2) / 2;
        assert_eq!(a, 5);
    }

    #[test]
    fn merge_sort_test() {
        let mut vec = vec![4, 5, 2, 6, 1, 8, 3, 5];
        fn compare(a: &i32, b: &i32) -> bool {
            if a < b {
                return true;
            } else {
                return false;
            }
        };

        let mut stock_vec = vec.clone();
        let len = vec.len();
        assert_eq!(len, 8);
        merge_sort(&mut vec, &mut stock_vec, compare, 0, len);
        assert_eq!(vec, [1, 2, 3, 4, 5, 5, 6, 8]);
    }
}
