use crate::aabb::Aabb;
use crate::material::MaterialHandle;
use crate::ray::Ray;
use crate::vec3::Vector3;

pub struct HitRecord {
    t: f64,
    u: f64,
    v: f64,
    p: Vector3<f64>,
    normal: Vector3<f64>,
    mat_ptr: MaterialHandle,
}

impl HitRecord {
    pub fn new(
        t: f64,
        u: f64,
        v: f64,
        p: Vector3<f64>,
        normal: Vector3<f64>,
        mat_ptr: MaterialHandle,
    ) -> Self {
        HitRecord {
            t,
            u,
            v,
            p,
            normal,
            mat_ptr,
        }
    }

    pub fn get_mat_ptr(&self) -> MaterialHandle {
        MaterialHandle(self.mat_ptr.0)
    }

    pub fn get_t(&self) -> f64 {
        self.t
    }

    pub fn get_u(&self) -> f64 {
        self.u
    }

    pub fn get_v(&self) -> f64 {
        self.v
    }

    pub fn get_p(&self) -> Vector3<f64> {
        self.p
    }

    pub fn get_normal(&self) -> Vector3<f64> {
        self.normal
    }
}

pub trait Hitable: HitableClone {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
    fn bounding_box(&self) -> Option<Aabb>;
    fn pdf_value(&self, _o: Vector3<f64>, _v: Vector3<f64>) -> f64 {
        return 0.0;
    }
    fn random(&self, _o: Vector3<f64>) -> Vector3<f64> {
        return [1.0, 0.0, 0.0];
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
