use rand::prelude::*;

mod vec3;
mod ray;
mod hitablelist;
mod hitable;
mod sphere;
mod camera;

use vec3::{Vector3, vec3_unit_vector_f64, vec3_mul_b, vec3_add, vec3_add_b, vec3_div_b};
use ray::{Ray};
use hitable::{Hitable};
use hitablelist::{HitableList};
use sphere::{Sphere};
use camera::{Camera};

fn color(r:  &Ray, world: &HitableList) -> Vector3<f64> {
    match world.hit(r, 0.0, 1000.0) {
        Some(rec) => vec3_mul_b(vec3_add_b(rec.normal, 1.0), 0.5),
        None => {
            let unit_direction = vec3_unit_vector_f64(r.direction());
            let t  = 0.5*(unit_direction[1] + 1.0);
            vec3_add(vec3_mul_b([1.0, 1.0, 1.0], 1.0 - t), vec3_mul_b([0.5, 0.7, 1.0], t))
        },
    }
}

fn main() {
    let nx: u32 = 200;
    let ny: u32 = 100;
    let ns: u32 = 100; //anti-aliasing sample-per-pixel
    let mut rng = rand::thread_rng();
    println!("P3\n {} {} \n255\n", nx, ny);
    let mut obj_list = HitableList::new();
    obj_list.push(Sphere { center: [0.0 , 0.0 , -1.0], radius: 0.5 });
    obj_list.push(Sphere { center: [0.0, -100.5, -1.0], radius: 100.0 });

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
                col = vec3_add(color(&r, &obj_list), col);
            }
            col = vec3_div_b(col, ns as f64);
            let ir: u32 = (255.99 * col[0]) as u32;
            let ig: u32 = (255.99 * col[1]) as u32;
            let ib: u32 = (255.99 * col[2]) as u32;
            println!("{} {} {}\n", ir, ig, ib);
            }
        }
}

