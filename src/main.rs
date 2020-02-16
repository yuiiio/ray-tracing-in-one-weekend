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
use material::{Metal, Lambertian, Materials};

fn color(r: &Ray, world: &HitableList, depth: u32, material_list: &Materials) -> Vector3<f64> {
    if depth < 50 {
        match world.hit(r, 0.00001, 10000.0) {
            Some(rec) => {
                if let Some(mat_rec) = material_list.get(rec.get_mat_ptr()).scatter(r, &rec) {
                    return vec3_mul(mat_rec.get_attenuation() , color(mat_rec.get_scatterd(), &world, depth + 1, material_list))
                }
            },
            None => {
                let unit_direction = vec3_unit_vector_f64(r.direction());
                let t  = 0.5*(unit_direction[1] + 1.0);
                return vec3_add(vec3_mul_b([1.0, 1.0, 1.0], 1.0 - t), vec3_mul_b([0.5, 0.7, 1.0], t))
            },
        }
    }
    [0.0, 0.0, 0.0]
}

fn main() {
    let nx: u32 = 800;
    let ny: u32 = 400;
    let ns: u32 = 100; //anti-aliasing sample-per-pixel
    let mut rng = rand::thread_rng();
    println!("P3\n {} {} \n255\n", nx, ny);
    let mut obj_list = HitableList::new();
    let mut material_list = Materials::new();
    let mat1 = material_list.add_material(Lambertian::new([0.8, 0.3, 0.3]));
    let mat2 = material_list.add_material(Lambertian::new([0.8, 0.8, 0.0]));
    let mat3 = material_list.add_material(Metal::new([0.8, 0.6, 0.2], 0.3));
    let mat4 = material_list.add_material(Metal::new([0.8, 0.8, 0.8], 0.1));
    obj_list.push(Sphere::new([0.0 , 0.0 , -1.0], 0.5, mat1));
    obj_list.push(Sphere::new([0.0, -100.5, -1.0], 100.0, mat2));
    obj_list.push(Sphere::new([1.0 , 0.0 , -1.0], 0.5, mat3));
    obj_list.push(Sphere::new([-1.0 , 0.0 , -1.0], 0.5, mat4));

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
                col = vec3_add(color(&r, &obj_list, 0, &material_list), col);
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

