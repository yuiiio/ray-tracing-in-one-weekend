use crate::aabb::Aabb;
use crate::hitable::{HitRecord, Hitable};
use crate::quotation::Rotation;
use crate::ray::Ray;
use crate::vec3::{vec3_add, vec3_inv, vec3_sub, Vector3};

#[derive(Clone)]
pub struct Translate {
    obj: Box<dyn Hitable + Send + Sync>,
    offset: Vector3<f64>,
    aabb_box: Aabb,
}

impl Translate {
    pub fn new(obj: Box<dyn Hitable + Send + Sync>, offset: Vector3<f64>) -> Self {
        let aabb_box = Aabb {
            b_min: vec3_add(&obj.bounding_box().b_min, &offset),
            b_max: vec3_add(&obj.bounding_box().b_max, &offset),
        };
        Translate {
            obj,
            offset,
            aabb_box,
        }
    }
}

impl Hitable for Translate {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let r = Ray {
            origin: vec3_sub(&r.origin, &self.offset),
            direction: r.direction,
            inv_dir: vec3_inv(&r.direction),
        };
        match self.obj.hit(&r, t_min, t_max) {
            Some(hit) => Some(HitRecord {
                t: hit.t,
                uv: hit.uv,
                p: vec3_add(&hit.p, &self.offset),
                normal: hit.normal,
                mat_ptr: hit.mat_ptr,
                onb_uv: hit.onb_uv,
            }),
            None => None,
        }
    }

    fn bounding_box(&self) -> &Aabb {
        &self.aabb_box
    }

    fn pdf_value(&self, ray: &Ray) -> f64 {
        let on = vec3_sub(&ray.origin, &self.offset);
        self.obj.pdf_value(&Ray {
            origin: on,
            direction: ray.direction,
            inv_dir: vec3_inv(&ray.direction),
        }) // this use self->obj's pdf_value func
    }
    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        let on = vec3_sub(o, &self.offset);
        self.obj.random(&on)
    }

    fn rotate_onb(&mut self, quat: &Rotation) -> () {
        self.obj.rotate_onb(quat);
    }
}

#[derive(Clone)]
pub struct Rotate {
    obj: Box<dyn Hitable + Send + Sync>,
    quat: Rotation,
    revq: Rotation,
    aabb_box: Aabb,
}

impl Rotate {
    pub fn new(mut obj: Box<dyn Hitable + Send + Sync>, axis: &Vector3<f64>, angle: f64) -> Self {
        let quat = Rotation::new(angle, axis);
        let revq = Rotation::new(-angle, axis);

        // let found boundingbox to enough include rotate-obj

        let aabb_box = obj.bounding_box_with_rotate(&quat);

        obj.rotate_onb(&quat); // rotate obj's normal and onb

        Rotate {
            obj,
            quat,
            revq,
            aabb_box,
        }
    }
}

impl Hitable for Rotate {
    fn hit(&self, r: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let origin = self.revq.rotate(&r.origin);
        let direction = self.revq.rotate(&r.direction);
        let r = Ray {
            origin,
            direction,
            inv_dir: vec3_inv(&direction),
        };
        match self.obj.hit(&r, t_min, t_max) {
            Some(hit) => {
                let normal = match hit.onb_uv {
                    Some(_onb_uv) => hit.normal, // norm and onb is static(eg. rect, triangle)
                    // already rotated
                    None => self.quat.rotate(&hit.normal), // norm and onb is not static(eg. sphere)
                };
                Some(HitRecord {
                    t: hit.t,
                    uv: hit.uv,
                    p: self.quat.rotate(&hit.p),
                    normal,
                    mat_ptr: hit.mat_ptr,
                    onb_uv: hit.onb_uv,
                })
            }
            None => None,
        }
    }

    fn bounding_box(&self) -> &Aabb {
        &self.aabb_box
    }

    fn pdf_value(&self, ray: &Ray) -> f64 {
        if let Some(_aabb_hit) = self.aabb_box.aabb_hit(ray, 0.00001, 10000.0) {
            let ro = self.revq.rotate(&ray.origin);
            let rv = self.revq.rotate(&ray.direction);
            return self.obj.pdf_value(&Ray {
                origin: ro,
                direction: rv,
                inv_dir: vec3_inv(&rv),
            });
        }
        0.0
    }
    fn random(&self, o: &Vector3<f64>) -> Vector3<f64> {
        let ro = self.revq.rotate(o);
        let rv = self.obj.random(&ro);
        self.quat.rotate(&rv)
    }

    fn rotate_onb(&mut self, quat: &Rotation) -> () {
        self.obj.rotate_onb(quat);
    }
}
