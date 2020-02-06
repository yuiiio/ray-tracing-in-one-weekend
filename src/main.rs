mod vec3;
mod ray;
mod hitablelist;
mod hitable;
mod sphere;

use vec3::{Vector3, vec3_unit_vector_f64, vec3_mul_b, vec3_add, vec3_add_b};
use ray::{Ray};
use hitable::{Hitable};
use hitablelist::{HitableList};
use sphere::{Sphere};

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
    println!("P3\n {} {} \n255\n", nx, ny);
    let lower_left_corner = [-2.0, -1.0, -1.0];
    let horizontal = [4.0, 0.0, 0.0];
    let vertical = [0.0, 2.0, 0.0];
    let origin = [0.0, 0.0, 0.0];
    let mut obj_list = HitableList::new();
    obj_list.push(Sphere { center: [0.0 , 0.0 , -1.0], radius: 0.5 });
    obj_list.push(Sphere { center: [0.0, -100.5, -1.0], radius: 100.0 });

    for j in (0..ny).rev() {
        for i in 0..nx {
            let u: f64 = i as f64/ nx as f64; 
            let v: f64 = j as f64/ ny as f64; 
            let r = Ray {
                a: origin,
                b: vec3_add(vec3_add(lower_left_corner, vec3_mul_b(horizontal, u)), vec3_mul_b(vertical, v))
            };
            let col = color(&r, &obj_list);
            let ir: u32 = (255.99 * col[0]) as u32;
            let ig: u32 = (255.99 * col[1]) as u32;
            let ib: u32 = (255.99 * col[2]) as u32;
            println!("{} {} {}\n", ir, ig, ib);
            }
        }
}

