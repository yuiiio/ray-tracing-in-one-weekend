use rand::prelude::*;
use image::{RgbaImage, Rgba};

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
use material::{Metal, Lambertian, Materials, Dielectric};

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
    const NX: usize = 200;
    const NY: usize = 100;
    let mut imgbuf = vec![vec![[0, 0, 0, 255]; NY]; NX];
    const NS: usize = 100; //anti-aliasing sample-per-pixel
    let mut rng = rand::thread_rng();
    let mut obj_list = HitableList::new();
    let mut material_list = Materials::new();
    let mat1 = material_list.add_material(Lambertian::new([0.3, 0.3, 0.8]));
    let mat2 = material_list.add_material(Lambertian::new([0.8, 0.8, 0.0]));
    let mat3 = material_list.add_material(Metal::new([0.8, 0.6, 0.2], 0.3));
    let mat4 = material_list.add_material(Dielectric::new(1.5));
    let mat5 = material_list.add_material(Dielectric::new(1.5));
    obj_list.push(Sphere::new([0.0 , 0.0 , -1.0], 0.5, mat1));
    obj_list.push(Sphere::new([0.0, -100.5, -1.0], 100.0, mat2));
    obj_list.push(Sphere::new([1.0 , 0.0 , -1.0], 0.5, mat3));
    obj_list.push(Sphere::new([-1.0 , 0.0 , -1.0], 0.5, mat4));
    obj_list.push(Sphere::new([-1.0 , 0.0 , -1.0], -0.45, mat5));

    let cam = Camera::new([-2.0, 2.0, 1.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0], 90.0, (NX/NY) as f64);

    for j in 0..NY {
        for i in 0..NX {
            let mut col = [0.0 as f64; 3];
            for _s in 0..NS {
                let rand_x: f64 = rng.gen();
                let rand_y: f64 = rng.gen();
                let u: f64 = (i as f64 + rand_x )/ NX as f64; 
                let v: f64 = (j as f64 + rand_y )/ NY as f64; 
                let r = cam.get_ray(u, v);
                col = vec3_add(color(&r, &obj_list, 0, &material_list), col);
            }
            col = vec3_div_b(col, NS as f64);
            col = [col[0].sqrt(), col[1].sqrt(), col[2].sqrt()];
            let ir: u8 = (255.99 * col[0]) as u8;
            let ig: u8 = (255.99 * col[1]) as u8;
            let ib: u8 = (255.99 * col[2]) as u8;
            imgbuf[i][j][0] = ir;
            imgbuf[i][j][1] = ig;
            imgbuf[i][j][2] = ib;
        }
    }

    let mut img = RgbaImage::new(NX as u32, NY as u32);
    for x in 0..NX {
        for y in 0..NY {
            img.put_pixel(x as u32, y as u32, Rgba(imgbuf[x][NY-(y+1)]));
        }
    }
    img.save("image.png").unwrap()
}