use crate::aabb::Aabb;
use crate::material::MaterialHandle;
use crate::quotation::Rotation;
use crate::ray::Ray;
use crate::vec3::Vector3;

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
        let mut b_min = [std::f64::MAX; 3];
        let mut b_max = [std::f64::MIN; 3];
        for i in [1, 0].iter() {
            for j in [1, 0].iter() {
                for k in [1, 0].iter() {
                    let x = *i as f64 * bbox.b_max[0] + (1 - *i) as f64 * bbox.b_min[0];
                    let y = *j as f64 * bbox.b_max[1] + (1 - *j) as f64 * bbox.b_min[1];
                    let z = *k as f64 * bbox.b_max[2] + (1 - *k) as f64 * bbox.b_min[2];

                    let rotated = quat.rotate(&[x, y, z]);
                    for c in 0..3 {
                        if rotated[c] > b_max[c] {
                            b_max[c] = rotated[c];
                        }
                        if rotated[c] < b_min[c] {
                            b_min[c] = rotated[c];
                        }
                    }
                }
            }
        }
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
