use crate::vec3::{Vector3, vec3_add, vec3_mul_b};

pub struct Ray {
    pub origin: Vector3<f64>,
    pub direction: Vector3<f64>,
    inv_dir: [Option<f64>; 3],
}

impl Ray {
    pub fn new(origin: Vector3<f64>, direction: Vector3<f64>) -> Self {
        Ray {
            origin, direction,
            inv_dir: [None; 3],
        }
    }

    pub fn get_inv_direction(&mut self, axis: usize) -> f64 {
        // cache inv direction
        match self.inv_dir[axis] {
            Some(inv_dir) => inv_dir,
            None => {
                let inv_dir = 1.0 / self.direction[axis];
                self.inv_dir[axis] = Some(inv_dir);
                inv_dir
            }
        }
    }

    pub fn point_at_parameter(&self, t: f64) -> Vector3<f64> {
        vec3_add(&self.origin, &vec3_mul_b(&self.direction, t))
    }
}
