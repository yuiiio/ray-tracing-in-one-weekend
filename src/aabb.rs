use crate::vec3::{Vector3};
use crate::ray::{Ray};
use std::mem::swap;

pub struct aabb {
    min: Vector3<f64>,
    max: Vector3<f64>,
}

impl aabb {
    pub fn new(min: Vector3<f64>, max: Vector3<f64>) -> aabb {
        aabb { min, max }
    }

    pub fn min(&self) -> Vector3<f64> {
        self.min
    }

    pub fn max(&self) -> Vector3<f64> {
        self.max
    }

    fn aabb_hit(&self, r: &Ray, t_min: f64, t_max: f64) -> bool {
        let mut tmin = t_min;
        let mut tmax = t_max;
        for i in 0..3 {
            let invD = 1.0 / r.direction()[i];
            let mut t0 = (self.min[i] - r.origin()[i]) / invD;
            let mut t1 = (self.max[i] - r.origin()[i]) / invD;

            if (invD < 0.0) {
                swap(&mut t0, &mut t1);
            }
        
            tmin = max(t0, tmin);
            tmax = min(t1, tmax);
            if tmax < tmin {
                return false
            }
        }
        return true
    }
}

fn max(a: f64, b: f64) -> f64 {
    if a < b {
        b
    } else {
        a
    }
}

fn min(a: f64, b: f64) -> f64 {
    if a > b {
        b
    } else {
        a
    }
}

mod test {
    use super::*;

    #[test]
    fn roop_let_test() {
        let mut s = 0;
        for i in 0..5 {
            s = s + i;
        }
        assert_eq!(s, 10);
    }

    #[test]
    fn aabb_hit_test() {
        let aabb_box = aabb::new([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);
        let r = Ray::new([0.0, 0.0, 0.0], [1.5, 1.5, 1.5]);
        let result = aabb_box.aabb_hit(&r, 0.00001, 10000.0);
        assert_eq!(true, result);
        let r = Ray::new([0.0, 0.0, 0.0], [1.5, 0.0, 1.5]);
        let result = aabb_box.aabb_hit(&r, 0.00001, 10000.0);
        assert_eq!(false, result);
    }
}