use rand::prelude::*;

mod vec3;
mod ray;
mod hitablelist;
mod hitable;
mod sphere;
mod camera;
mod material;

use vec3::{Vector3, vec3_unit_vector_f64, vec3_mul_b, vec3_add, vec3_div_b, vec3_mul};
use ray::{Ray};
use hitable::{Hitable};
use hitablelist::{HitableList};
use sphere::{Sphere};
use camera::{Camera};
use std::f64;
use material::{Metal, Lambertion};

fn color(r: &Ray, world: &HitableList, depth: u32) -> Vector3<f64> {
    match world.hit(r, 0.0, 1000.0) {
        Some(rec) => {
            if depth < 50 {
                if let Some(mat_rec) = rec.mat_ptr.scatter(r, &rec) {
                    return vec3_mul(mat_rec.attenuation , color(&mat_rec.scatterd, &world, depth + 1))
                }
            }
            [0.0, 0.0, 0.0]
        },
        None => {
            let unit_direction = vec3_unit_vector_f64(r.direction());
            let t  = 0.5*(unit_direction[1] + 1.0);
            vec3_add(vec3_mul_b([1.0, 1.0, 1.0], 1.0 - t), vec3_mul_b([0.5, 0.7, 1.0], t))
        },
    }
}

fn main() {
    let nx: u32 = 400;
    let ny: u32 = 200;
    let ns: u32 = 100; //anti-aliasing sample-per-pixel
    let mut rng = rand::thread_rng();
    println!("P3\n {} {} \n255\n", nx, ny);
    let mut obj_list = HitableList::new();
    obj_list.push(Sphere { center: [0.0 , 0.0 , -1.0], radius: 0.5 , mat_ptr: Box::new(Lambertion{ albedo: [0.8, 0.3, 0.3] }) });
    obj_list.push(Sphere { center: [0.0, -100.5, -1.0], radius: 100.0 , mat_ptr: Box::new(Lambertion{ albedo: [0.8, 0.8, 0.0] }) });
    obj_list.push(Sphere { center: [1.0 , 0.0 , -1.0], radius: 0.5 , mat_ptr: Box::new(Metal{ albedo: [0.8, 0.6, 0.2] }) });
    obj_list.push(Sphere { center: [-1.0 , 0.0 , -1.0], radius: 0.5 , mat_ptr: Box::new(Metal{ albedo: [0.8, 0.8, 0.8] }) });

    let cam = Camera::new();

    for j in (0..ny).rev() {
        for i in 0..nx {
            let mut col = [0.0, 0.0, 0.0];
            for _s in 0..ns {
                let rand_x: f64 = rng.gen();
                let rand_y: f64 = rng.gen();
                let u: f64 = (i as f64 + rand_x )/ nx as f64; 
                let v: f64 = (j as f64 + rand_y )/ ny as f64; 
                let r = cam.get_ray(u, v);
                col = vec3_add(color(&r, &obj_list, 0), col);
            }
            col = vec3_div_b(col, ns as f64);
            col = [col[0].sqrt(), col[1].sqrt(), col[2].sqrt()];
            let ir: u32 = (255.99 * col[0]) as u32;
            let ig: u32 = (255.99 * col[1]) as u32;
            let ib: u32 = (255.99 * col[2]) as u32;
            println!("{} {} {}\n", ir, ig, ib);
            }
        }
}

