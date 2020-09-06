use rand::prelude::*;
use image::{RgbaImage, Rgba, open};
use std::time::SystemTime;
use std::thread;
use std::sync::{Arc, Mutex};

mod vec3;
mod ray;
mod hitablelist;
mod hitable;
mod sphere;
mod camera;
mod material;
mod texture;
mod aabb;
mod utils;
mod bvh_node;
mod rectangle;
mod translate;
mod quotation;

use vec3::{Vector3, vec3_unit_vector_f64, vec3_mul_b, vec3_add, vec3_mul, vec3_div};
use ray::{Ray};
use hitable::{Hitable};
use hitablelist::{HitableList};
use sphere::{Sphere};
use camera::{Camera};
use std::f64;
use material::{Metal, Lambertian, Materials, Dielectric, DiffuseLight};
use texture::{ColorTexture, CheckerTexture, ImageTexture};
use utils::{clamp};
use bvh_node::{BvhNode};
use rectangle::{Rect, AxisType, FlipNormals, Boxel};
use translate::{Translate, Rotate};

fn color<T: Hitable>(r: &Ray, world: &Arc<T>, depth: u32, material_list: &Materials, last_absorabance: Vector3<f64>) -> Vector3<f64> {
    if depth < 50 {
        match world.hit(r, 0.00001, 10000.0) {
            Some(rec) => {
                let emitted = material_list.get(rec.get_mat_ptr()).emitted(r, &rec);
                if let Some(mat_rec) = material_list.get(rec.get_mat_ptr()).scatter(r, &rec) {
                    let absorabance = vec3_mul_b(last_absorabance, rec.get_t() * rec.get_t());
                    let absorabance = vec3_div([1.0, 1.0, 1.0], absorabance);
                    let absorabance = [ clamp(absorabance[0], 0.0, 1.0),
                                        clamp(absorabance[1], 0.0, 1.0),
                                        clamp(absorabance[2], 0.0, 1.0), ];
                    return vec3_add( emitted, vec3_mul(
                        vec3_mul(mat_rec.get_attenuation(),
                        color(mat_rec.get_scatterd(), world, depth + 1, material_list, mat_rec.get_absorabance())),
                        absorabance))
                }
                return emitted
            },
            None => {
            },
        }
    }
    [0.0, 0.0, 0.0]
}

fn main() {
    let now = SystemTime::now();
    const NX: usize = 400;
    const NY: usize = 400;
    let imgbuf = Arc::new(Mutex::new(vec![vec![[0, 0, 0, 255]; NY]; NX]));
    const NS: usize = 500; //anti-aliasing sample-per-pixel
    let mut obj_list = HitableList::new();
    let mut material_list = Materials::new();

    let red = material_list.add_material(Lambertian::new(ColorTexture::new([0.65, 0.05, 0.05])));
    let white = material_list.add_material(Lambertian::new(ColorTexture::new([0.73, 0.73, 0.73])));
    let green = material_list.add_material(Lambertian::new(ColorTexture::new([0.12, 0.45, 0.15])));
    let light = material_list.add_material(DiffuseLight::new(ColorTexture::new([15.0, 15.0, 15.0])));
    let magick = material_list.add_material(Lambertian::new(ImageTexture::new(open("./texture.png").unwrap().into_rgba())));
    let glass = material_list.add_material(Dielectric::new(2.0, [0.01, 0.01, 0.0]));
    let metal = material_list.add_material(Metal::new(0.0, ColorTexture::new([0.8, 0.85, 0.88])));

    obj_list.push(FlipNormals::new(Rect::new(0.0, 555.0 , 0.0, 555.0, 555.0, AxisType::kYZ, green)));
    obj_list.push(Rect::new(0.0, 555.0 , 0.0, 555.0, 0.0, AxisType::kYZ, red));
    obj_list.push(Rect::new(213.0, 343.0 , 227.0, 332.0, 554.0, AxisType::kXZ, light));
    obj_list.push(FlipNormals::new(Rect::new(0.0, 555.0 , 0.0, 555.0, 555.0, AxisType::kXZ, white)));
    obj_list.push(Rect::new(0.0, 555.0 , 0.0, 555.0, 0.0, AxisType::kXZ, white));
    obj_list.push(FlipNormals::new(Rect::new(0.0, 555.0 , 0.0, 555.0, 555.0, AxisType::kXY, magick)));

    obj_list.push(Translate::new(Box::new(
                    Rotate::new(Box::new(
                        Boxel::new([0.0, 0.0, 0.0], [165.0, 165.0, 165.0], white)
                    ), [0.0, 1.0, 0.0], -18.0)
                ), [130.0, 0.0, 65.0]));
    obj_list.push(Translate::new(Box::new(
                    Rotate::new(Box::new(
                        Boxel::new([0.0, 0.0, 0.0], [165.0, 330.0, 165.0], metal)
                    ), [0.0, 1.0, 0.0], 15.0)
                ), [265.0, 0.0, 295.0]));

    obj_list.push(Sphere::new([455.0, 100.0, 100.0], 100.0, glass));

    let obj_list = BvhNode::new(&mut obj_list);

    let cam = Camera::new([278.0, 278.0, -800.0], [278.0, 278.0, 0.0], [0.0, 1.0, 0.0], 40.0, (NX/NY) as f64);

    let cam = Arc::new(cam);
    let obj_list = Arc::new(obj_list);
    let material_list = Arc::new(material_list);
    let mut handles = vec![];
    for j in 0..NY {
        for i in 0..NX {
            let imgbuf_clone = Arc::clone(&imgbuf);
            let cam = Arc::clone(&cam);
            let obj_list = Arc::clone(&obj_list);
            let material_list = Arc::clone(&material_list);
            let handle = thread::spawn(move || {
                let mut col = [0.0 as f64; 3];
                let mut rng = rand::thread_rng();
                for _s in 0..NS {
                    let rand_x: f64 = rng.gen();
                    let rand_y: f64 = rng.gen();
                    let u: f64 = (i as f64 + rand_x) / NX as f64;
                    let v: f64 = (j as f64 + rand_y) / NY as f64;
                    let r = cam.get_ray(u, v);
                    col = vec3_add(color(&r, &obj_list, 0, &material_list, [0.0, 0.0, 0.0]), col);
                }
                let c = 1.0 / NS as f64;
                col = vec3_mul_b(col, c);
                col = [col[0].sqrt(), col[1].sqrt(), col[2].sqrt()];
                col = [clamp(col[0], 0.0, 1.0), clamp(col[1], 0.0, 1.0), clamp(col[2], 0.0, 1.0)];
                let ir: u8 = (255.99 * col[0]) as u8;
                let ig: u8 = (255.99 * col[1]) as u8;
                let ib: u8 = (255.99 * col[2]) as u8;
                let mut imgbuf = imgbuf_clone.lock().unwrap();
                imgbuf[i][j] = [ir, ig, ib, 255];
            });
            handles.push(handle);
        }
    }

    for handle in handles {
        handle.join().unwrap();
    }

    println!("Time elapsed: {}",now.elapsed().unwrap().as_secs_f64());

    let mut img = RgbaImage::new(NX as u32, NY as u32);
    for x in 0..NX {
        for y in 0..NY {
            let imgbuf = imgbuf.lock().unwrap();
            img.put_pixel(x as u32, y as u32, Rgba(imgbuf[x][NY-(y+1)]));
        }
    }
    img.save("image.png").unwrap()
}