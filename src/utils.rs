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

pub fn qsort<T: Copy>(vec: &mut Vec<T>, compare: fn(&T, &T) -> bool) -> Vec<T> {
    if vec.len() <= 1 {
        return vec.to_vec()
    }
    let x = vec.len() / 2;// chose random
    let pipot = vec[x];
    let mut a:Vec<T> = Vec::new();
    let mut b:Vec<T> = Vec::new();
    for i in x+1..vec.len() {
        if compare(&vec[i], &pipot) {
            vec.remove(i);
            a.push(vec[i]);
        }
    }
    for i in 0..x-1 {
        if compare(&pipot, &vec[i]) {
            vec.remove(i);
            b.push(vec[i]);
        }
    }
    let mut a_sorted = qsort(&mut a, compare);
    let mut b_sorted = qsort(&mut b, compare);
    a_sorted.append(&mut b_sorted);
    a_sorted
}

mod test {
    use super::*;

    #[test]
    fn qsort_test() {
        let mut vec = vec![4, 5, 2, 6, 1, 8, 3, 5];
        fn compare (a: &i32, b: &i32) -> bool { 
            if a <= b {
                return true;
            } else {
                return false;
            }
        };

        let sorted = qsort(&mut vec, compare);
        println!("{:?}", sorted);
        assert_eq!(vec, [1, 2, 3, 4, 5, 5, 6, 8]);
    }
}