use crate::aabb::Aabb;
use crate::material::MaterialHandle;
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
    fn pdf_value(&self, _r: &Ray) -> f64 {
        0.0
    }
    fn random(&self, _o: &Vector3<f64>) -> Vector3<f64> { // should return normalized vector
        [1.0, 0.0, 0.0]
    }
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
