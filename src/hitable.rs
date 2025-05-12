use crate::aabb::Aabb;
use crate::material::MaterialHandle;
use crate::quotation::Rotation;
use crate::ray::Ray;
use crate::vec3::{vec3_a_abs_mul_f64, vec3_add, vec3_mul_b, vec3_sub, Vector3};

pub struct HitRecord<'a> {
    pub t: f64,
    pub uv: (f64, f64),
    pub p: Vector3<f64>,
    pub normal: Vector3<f64>,
    pub mat_ptr: &'a MaterialHandle,
    pub onb_uv: Option<&'a (Vector3<f64>, Vector3<f64>)>,
}

pub trait Hitable: HitableClone {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
    fn bounding_box(&self) -> &Aabb;
    fn bounding_box_with_rotate(&self, quat: &Rotation) -> Aabb {
        let bbox = self.bounding_box();
        // https://zeux.io/2010/10/17/aabb-from-obb-with-component-wise-abs/
        //let center = vec3_mul_b(&vec3_add(&bbox.b_max, &bbox.b_min), 0.5);
        // min + ((max - min)/2) = (min + max)/2
        // extent( same distant from 8 vertex of bbox ) always positive
        // (max - min)/2
        let extent = vec3_mul_b(&vec3_sub(&bbox.b_max, &bbox.b_min), 0.5);
        let center = vec3_add(&bbox.b_min, &extent);

        // abs-rotate-matrix dot
        let abs_matrix_rotated_extent: [f64; 3] = [
            vec3_a_abs_mul_f64(&quat.left_side_matrix[0], &extent)
                .iter()
                .sum(),
            vec3_a_abs_mul_f64(&quat.left_side_matrix[1], &extent)
                .iter()
                .sum(),
            vec3_a_abs_mul_f64(&quat.left_side_matrix[2], &extent)
                .iter()
                .sum(),
        ];

        let rotated_center = quat.rotate(&center);

        let b_min = vec3_sub(&rotated_center, &abs_matrix_rotated_extent);
        let b_max = vec3_add(&rotated_center, &abs_matrix_rotated_extent);

        Aabb { b_min, b_max }
    }
    fn pdf_value(&self, _r: &Ray) -> f64 {
        0.0
    }
    fn random(&self, _o: &Vector3<f64>) -> Vector3<f64> {
        // should return normalized vector
        [1.0, 0.0, 0.0]
    }
    fn rotate_onb(&mut self, quat: &Rotation) -> (); // rotate onb and normal vec used at build time
    fn scale(&mut self, scale_value: f64);
}

pub trait HitableClone {
    fn clone_box(&self) -> Box<dyn Hitable + Send + Sync>;
}

impl<T> HitableClone for T
where
    T: 'static + Hitable + Send + Sync + Clone,
{
    fn clone_box(&self) -> Box<dyn Hitable + Send + Sync> {
        Box::new(self.clone())
    }
}

impl Clone for Box<dyn Hitable + Send + Sync> {
    fn clone(&self) -> Box<dyn Hitable + Send + Sync> {
        self.clone_box()
    }
}
