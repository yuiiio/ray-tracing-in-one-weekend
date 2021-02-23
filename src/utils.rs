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

pub fn qsort<T: Clone>(vec: &mut Vec<T>, compare: fn(&T, &T) -> bool) {
    let start = 0;
    let end = vec.len() - 1;
    qsort_partition(vec, start, end as isize, compare);
}

fn qsort_partition<T: Clone>(
    vec: &mut Vec<T>,
    start: isize,
    end: isize,
    compare: fn(&T, &T) -> bool,
) {
    if start < end && end - start >= 1 {
        let pivot = partition(vec, start as isize, end as isize, compare);
        qsort_partition(vec, start, pivot - 1, compare);
        qsort_partition(vec, pivot + 1, end, compare);
    }
}

fn partition<T: Clone>(vec: &mut Vec<T>, l: isize, h: isize, compare: fn(&T, &T) -> bool) -> isize {
    let pivot = vec[h as usize].clone();
    let mut i = l - 1;

    for j in l..h {
        if compare(&vec[j as usize], &pivot) {
            i = i + 1;
            let temp = vec[i as usize].clone();
            vec[i as usize] = vec[j as usize].clone();
            vec[j as usize] = temp;
        }
    }

    let temp = vec[(i + 1) as usize].clone();
    vec[(i + 1) as usize] = vec[h as usize].clone();
    vec[h as usize] = temp;

    i + 1
}

mod test {
    use super::*;

    #[test]
    fn qsort_test() {
        let mut vec = vec![4, 5, 2, 6, 1, 8, 3, 5];
        fn compare(a: &i32, b: &i32) -> bool {
            if a < b {
                return true;
            } else {
                return false;
            }
        };
        qsort(&mut vec, compare);
        assert_eq!(vec, [1, 2, 3, 4, 5, 5, 6, 8]);
    }
}
