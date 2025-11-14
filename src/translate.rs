use crate::aabb::Aabb;
use crate::hitable::{HitRecord, Hitable};
use crate::material::MaterialHandle;
use crate::quotation::Rotation;
use crate::ray::Ray;
use crate::vec3::{vec3_add, vec3_sub, Vector3};

#[derive(Clone)]
pub struct Translate {
    obj: Box<dyn Hitable + Send + Sync>,
    offset: Vector3<f64>,
}

impl Translate {
    pub fn new(obj: Box<dyn Hitable + Send + Sync>, offset: Vector3<f64>) -> Self {
        Translate {
            obj,
            offset,
        }
    }
}

impl Hitable for Translate {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let r = Ray {
            origin: vec3_sub(&r.origin, &self.offset),
            direction: r.direction,
        };
        match self.obj.hit(&r, t_min, t_max) {
            Some(hit) => Some(HitRecord {
                t: hit.t,
                uv: hit.uv,
                p: vec3_add(&hit.p, &self.offset),
                normal: hit.normal,
                mat_ptr: hit.mat_ptr,
            }),
            None => None,
        }
    }

    fn bounding_box(&self) -> Aabb {
        Aabb {
            b_min: vec3_add(&self.obj.bounding_box().b_min, &self.offset),
            b_max: vec3_add(&self.obj.bounding_box().b_max, &self.offset),
        }
    }

    fn pdf_value(&self, ray: &Ray) -> f64 {
        let on = vec3_sub(&ray.origin, &self.offset);
        self.obj.pdf_value(&Ray {
            origin: on,
            direction: ray.direction,
        }) // this use self->obj's pdf_value func
    }
    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        let on = vec3_sub(o, &self.offset);
        self.obj.random(&on)
    }

    fn scale(&mut self, scale_value: f64) {
        self.obj.scale(scale_value);

        //scale func is affected after translated(obj)
    }
    fn set_material(&mut self, mat_ptr: MaterialHandle) {
        self.obj.set_material(mat_ptr);
    }
}

#[derive(Clone)]
pub struct Rotate {
    obj: Box<dyn Hitable + Send + Sync>,
    quat: Rotation,
    revq: Rotation,
}

impl Rotate {
    pub fn new(obj: Box<dyn Hitable + Send + Sync>, quat: Rotation) -> Self {
        let revq = quat.get_revq();

        Rotate {
            obj,
            quat,
            revq,
        }
    }
}

impl Hitable for Rotate {
    fn hit(&self, input_ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let origin = self.revq.rotate(&input_ray.origin);
        let direction = self.revq.rotate(&input_ray.direction);
        let revq_r = Ray { origin, direction };
        match self.obj.hit(&revq_r, t_min, t_max) {
            Some(hit) => {
                let normal = self.quat.rotate(&hit.normal);
                Some(HitRecord {
                    t: hit.t,
                    uv: hit.uv,
                    p: input_ray.point_at_parameter(hit.t), //self.quat.rotate(&hit.p),
                    normal,
                    mat_ptr: hit.mat_ptr,
                })
            }
            None => None,
        }
    }

    fn bounding_box(&self) -> Aabb {
        // let found boundingbox to enough include rotate-obj
        self.obj.bounding_box_with_rotate(&self.quat)
    }

    fn pdf_value(&self, ray: &Ray) -> f64 {
        let ro = self.revq.rotate(&ray.origin);
        let rv = self.revq.rotate(&ray.direction);
        return self.obj.pdf_value(&Ray {
            origin: ro,
            direction: rv,
        });
    }
    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        let ro = self.revq.rotate(o);
        let rv = self.obj.random(&ro);
        self.quat.rotate(&rv)
    }

    fn scale(&mut self, scale_value: f64) {
        self.obj.scale(scale_value);

        //scale func is affected after translated(obj)
    }
    fn set_material(&mut self, mat_ptr: MaterialHandle) {
        self.obj.set_material(mat_ptr);
    }
}
