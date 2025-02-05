use image::{open, Rgba, RgbaImage};
use rand::prelude::*;
use std::fs::File;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering::Relaxed;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::SystemTime;

mod aabb;
mod bvh_node;
mod camera;
mod hitable;
mod hitablelist;
mod material;
mod obj_loader;
mod onb;
mod pdf;
mod quotation;
mod ray;
mod rectangle;
mod sphere;
mod texture;
mod translate;
mod triangle;
mod utils;
mod vec3;

use bvh_node::BvhTree;
use camera::Camera;
use hitable::Hitable;
use hitablelist::HitableList;
use material::{Dielectric, DiffuseLight, Lambertian, MaterialList, Metal, Scatterd};
use obj_loader::obj_loader;
use onb::Onb;
use pdf::{cosine_pdf_generate, cosine_pdf_value};
use ray::Ray;
use rectangle::{AxisType, Boxel, FlipNormals, Rect};
use sphere::Sphere;
use std::f64;
use texture::{ColorTexture, ImageTexture, TextureList};
use translate::{Rotate, Translate};
use vec3::{vec3_add, vec3_mul, vec3_mul_b, Vector3};

const MAX_DEPTH: usize = 20;

fn color(
    ray: Ray,
    world: &BvhTree,
    light_list: &HitableList,
    texture_list: &TextureList,
    material_list: &MaterialList,
) -> Vector3<f64> {
    let mut cur_emitted: Vector3<f64> = [0.0, 0.0, 0.0];
    let mut last_throughput: Vector3<f64> = [1.0, 1.0, 1.0];
    let mut last_absorabance: Vector3<f64> = [0.0, 0.0, 0.0];
    let mut ray: Ray = ray;
    for _i in 0..MAX_DEPTH {
        //println!("ray direction length: {}", vec3::vec3_length_f64(&ray.direction));
        match world.hit(&ray, 0.00001, 10000.0) {
            Some(hit_rec) => {
                // material obj: scatter,  emitted scattering_pdf,
                // emitted is must called.
                // if scatter get Pdf then called scattering_pdf too.
                let last_emitted = material_list.emitted(&ray, &hit_rec, texture_list);
                // cur_emitted should calc before last_throughput
                cur_emitted = vec3_add(&cur_emitted, &vec3_mul(&last_throughput, &last_emitted));
                // mat_rec attenuation, absorabance scatterd(::Ray, ::Pdf)
                if let Some(mat_rec) = material_list.scatter(&ray, &hit_rec, texture_list) {
                    let distance: f64 = hit_rec.t;
                    let mut absorabance: Vector3<f64> = [1.0, 1.0, 1.0];
                    if last_absorabance[0] != 0.0 {
                        absorabance[0] = f64::exp(-(last_absorabance[0] * distance));
                    }
                    if last_absorabance[1] != 0.0 {
                        absorabance[1] = f64::exp(-(last_absorabance[1] * distance));
                    }
                    if last_absorabance[2] != 0.0 {
                        absorabance[2] = f64::exp(-(last_absorabance[2] * distance));
                    }
                    let attenuation = mat_rec.attenuation;
                    last_absorabance = mat_rec.absorabance;
                    match mat_rec.scatterd {
                        Scatterd::Ray(next_ray) => {
                            last_throughput =
                                vec3_mul(&last_throughput, &vec3_mul(&attenuation, &absorabance));
                            ray = next_ray;
                            continue;
                        }
                        Scatterd::CosinePdf => {
                            let mut rng = rand::thread_rng();
                            let rand: f64 = rng.gen();
                            let next_ray = if rand < 0.4 {
                                Ray {
                                    origin: hit_rec.p,
                                    direction: light_list.random(&hit_rec.p),
                                }
                            } else {
                                let direction = match hit_rec.onb_uv {
                                    Some(onb_uv) => cosine_pdf_generate(&Onb {
                                        u: onb_uv.0,
                                        v: onb_uv.1,
                                        w: hit_rec.normal,
                                    }),
                                    None => {
                                        cosine_pdf_generate(&Onb::build_from_w(&hit_rec.normal))
                                    }
                                };
                                Ray {
                                    origin: hit_rec.p,
                                    direction,
                                }
                            }; // next_ray direction should normalized value.

                            let light_list_pdf = light_list.pdf_value(&next_ray);
                            let cosine_pdf = cosine_pdf_value(&hit_rec.normal, &next_ray.direction);

                            let pdf_value = light_list_pdf * 0.4 + cosine_pdf * 0.6;

                            if pdf_value > 0.0 {
                                let spdf_value = material_list.scattering_pdf(&next_ray, &hit_rec);
                                let albedo = vec3_mul_b(&attenuation, spdf_value);
                                let nor_pdf_value = 1.0 / pdf_value;
                                last_throughput = vec3_mul(
                                    &last_throughput,
                                    &vec3_mul_b(&vec3_mul(&albedo, &absorabance), nor_pdf_value),
                                );
                                ray = next_ray;
                                continue;
                            } else {
                                return cur_emitted;
                            };
                        }
                    };
                };
                return cur_emitted;
            }
            None => {
                // if not hit any obj
                /*
                // sky
                let a = (ray.direction[1] + 1.0) * 0.5;
                let last_emitted = vec3_add(
                    &vec3_mul_b(&[0.5, 0.1, 0.05], 1.0 - a),
                    &vec3_mul_b(&[0.1, 0.1, 0.5], a),
                );
                */
                let last_emitted = [0.01, 0.01, 0.01];
                cur_emitted = vec3_add(&cur_emitted, &vec3_mul(&last_throughput, &last_emitted));
                return cur_emitted;
            }
        }
    }
    // when nest >= dephs
    let last_emitted = [0.5, 0.5, 0.5];
    vec3_add(&cur_emitted, &vec3_mul(&last_throughput, &last_emitted))
}

fn main() {
    let now = SystemTime::now();
    const OUTPUT_X: usize = 1920;
    const OUTPUT_Y: usize = 1080;
    const NS: usize = 8; // x^2 / per pixel sample size;
    const NX: usize = OUTPUT_X * NS;
    const NY: usize = OUTPUT_Y * NS;

    let imgbuf = Mutex::new(vec![[[0.0, 0.0, 0.0]; NY]; NX]);
    const N_TASK: usize = 1024; // task number ( expect larger than real cpu threads )
    let cpu_threads: usize = thread::available_parallelism().unwrap().get();

    let mut obj_list = HitableList::new();
    let mut light_list = HitableList::new();
    let mut material_list = MaterialList::new();

    let mut texture_list = TextureList::new();
    let red_texture = texture_list.add_color_texture(ColorTexture::new([0.65, 0.05, 0.05]));
    let white_texture = texture_list.add_color_texture(ColorTexture::new([0.73, 0.73, 0.73]));
    let green_texture = texture_list.add_color_texture(ColorTexture::new([0.12, 0.45, 0.15]));
    let light_texture = texture_list.add_color_texture(ColorTexture::new([20.0, 15.0, 10.0]));
    let magick_texture = texture_list.add_image_texture(ImageTexture::new(
        open("./texture.png").unwrap().into_rgba8(),
    ));
    let earth_texture = texture_list.add_image_texture(ImageTexture::new(
        open("./texture.jpg").unwrap().into_rgba8(),
    ));
    let metal_texture = texture_list.add_color_texture(ColorTexture::new([0.5, 0.7, 0.7]));
    //let fuzzy_metal_texture = texture_list.add_color_texture(ColorTexture::new([0.7, 0.7, 0.7]));

    let red = material_list.add_lambertian_mat(Lambertian::new(red_texture));
    let white = material_list.add_lambertian_mat(Lambertian::new(white_texture));
    let green = material_list.add_lambertian_mat(Lambertian::new(green_texture));
    let light = // light looks good on 1.0 ~ 0.0, because { emitted + (nasted result) } * accum(0.0 ~ 1.0), over flow and overflow on next path
                // but, > 1.0 can happen when powerfull light ?
        material_list.add_diffuselight_mat(DiffuseLight::new(light_texture));
    let magick = material_list.add_lambertian_mat(Lambertian::new(magick_texture));
    let earth = material_list.add_lambertian_mat(Lambertian::new(earth_texture));
    let glass = material_list.add_dielectric_mat(Dielectric::new(1.5, [0.009, 0.006, 0.0]));
    let red_glass = material_list.add_dielectric_mat(Dielectric::new(1.5, [0.005, 0.03, 0.045]));
    let metal = material_list.add_metal_mat(Metal::new(0.01, metal_texture));
    //let fuzzy_metal = material_list.add_metal_mat(Metal::new(0.1, fuzzy_metal_texture));

    /*
    obj_list.push(FlipNormals::new(Rect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        AxisType::Kyz,
        green,
    )));
    obj_list.push(Rect::new(0.0, 555.0, 0.0, 555.0, 0.0, AxisType::Kyz, red));
    */

    /*
    let light_rect = FlipNormals::new(Rect::new(
        213.0,
        343.0,
        227.0,
        332.0,
        554.0,
        AxisType::Kxz,
        light,
    ));
    obj_list.push(light_rect.clone());
    */

    /*
    obj_list.push(FlipNormals::new(Rect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        AxisType::Kxz,
        white.clone(),
    )));

    obj_list.push(FlipNormals::new(Rect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        AxisType::Kxy,
        magick,
    )));
    */

    let floor = Rect::new(
        -10000.0,
        10000.0,
        -10000.0,
        10000.0,
        0.0,
        AxisType::Kxz,
        white.clone(),
    );
    obj_list.push(floor.clone());

    let glass_box = Translate::new(
        Box::new(Rotate::new(
            Box::new(Boxel::new(
                [0.0, 0.0, 0.0],
                [165.0, 165.0, 165.0],
                red_glass,
            )),
            &[0.0, 1.0, 0.0],
            -18.0,
        )),
        [130.0, 0.0, 65.0],
    );
    obj_list.push(glass_box.clone());

    let magick_box = Translate::new(
        Box::new(Rotate::new(
            Box::new(Boxel::new([0.0, 0.0, 0.0], [130.0, 130.0, 130.0], magick)),
            &[0.5, 0.5, 0.5],
            -30.0,
        )),
        [80.0, 330.0, 90.0],
    );
    obj_list.push(magick_box.clone());

    let metal_box = Translate::new(
        Box::new(Rotate::new(
            Box::new(Boxel::new(
                [0.0, 0.0, 0.0],
                [165.0, 330.0, 165.0],
                metal.clone(),
            )),
            &[0.0, 1.0, 0.0],
            15.0,
        )),
        [265.0, 0.00, 295.0],
    );
    obj_list.push(metal_box.clone());

    let glass_sphere = Sphere::new([455.0, 100.0, 100.0], 100.0, glass.clone());
    obj_list.push(glass_sphere.clone());
    let earth_sphere = Sphere::new([500.0, 300.0, 100.0], 60.0, earth);
    obj_list.push(earth_sphere);

    let light_sphere = Sphere::new([455.0, 400.0, 100.0], 50.0, light);
    obj_list.push(light_sphere.clone());

    let bunny_list = obj_loader(&mut File::open("./lucy.obj").unwrap(), white, 0.2);

    let now1 = SystemTime::now();
    let bunny_bvh = BvhTree::new(bunny_list);
    println!(
        "BVH-1 Build Time elapsed: {}",
        now1.elapsed().unwrap().as_secs_f64()
    );

    let translated_bunny_bvh = Translate::new(
        Box::new(Rotate::new(Box::new(bunny_bvh), &[1.0, 0.5, 0.0], 270.0)),
        [130.0, 90.0, -50.0],
    );

    obj_list.push(translated_bunny_bvh.clone());

    /*
    let glass_box = obj_loader(&mut File::open("./box.obj").unwrap(), glass);

    let glass_box = Translate::new(
        Box::new(Rotate::new(
            Box::new(glass_box),
            //Box::new(Boxel::new([0.0, 0.0, 0.0], [100.0, 100.0, 100.0], glass)),
            &[1.0, 1.0, 0.0],
            45.0,
        )),
        [200.0, 300.0, 100.0],
    );

    obj_list.push(glass_box.clone());
    */

    let now2 = SystemTime::now();
    let obj_bvh = BvhTree::new(obj_list);
    println!(
        "BVH-2 Build Time elapsed: {}",
        now2.elapsed().unwrap().as_secs_f64()
    );

    //light_list.push(light_rect);
    //light_list.push(floor);
    light_list.push(light_sphere);
    light_list.push(metal_box);
    light_list.push(glass_sphere);
    light_list.push(glass_box);

    let cam = Camera::new(
        [278.0, 278.0, -800.0],
        &[278.0, 278.0, 0.0],
        &[0.0, 1.0, 0.0],
        40.0,
        NX as f64 / NY as f64,
    );

    const ALIGHN_X: usize = NX / N_TASK; //small
    const AX: usize = NX / ALIGHN_X; //large
    let mut axa: [usize; AX] = [ALIGHN_X; AX];
    const FIZZ: usize = NX % ALIGHN_X;
    if FIZZ != 0 {
        axa[AX - 1] = FIZZ;
    }

    let cam_borrowed = &cam;
    let obj_bvh_borrowed = &obj_bvh;
    let light_list_borrowed = &light_list;
    let texture_list_borrowed = &texture_list;
    let material_list_borrowed = &material_list;

    let imgbuf_arc = Arc::new(&imgbuf);

    let main_thread = &thread::current();
    let avaiable_thread = &AtomicUsize::new(cpu_threads);
    thread::scope(|s| {
        for j in 0..AX {
            let imgbuf_clone = Arc::clone(&imgbuf_arc);
            let current_avaiable_thread = avaiable_thread.load(Relaxed);
            if current_avaiable_thread == 0 {
                thread::park();
            }
            avaiable_thread.fetch_sub(1, Relaxed);
            s.spawn(move || {
                let mut img_box: Vec<[[f64; 3]; NY]> = Vec::with_capacity(axa[j]);
                for in_j in 0..axa[j] {
                    let mut img_line: [[f64; 3]; NY] = [[0.0; 3]; NY];
                    for i in 0..NY {
                        let u: f64 = (((j * ALIGHN_X) + in_j) as f64) / NX as f64;
                        let v: f64 = (i as f64) / NY as f64;
                        let r = cam_borrowed.get_ray(u, v);
                        let col = color(
                            r,
                            obj_bvh_borrowed,
                            light_list_borrowed,
                            texture_list_borrowed,
                            material_list_borrowed,
                        );
                        img_line[i] = col;
                    }
                    img_box.push(img_line);
                }
                let mut imgbuf = imgbuf_clone.lock().unwrap();
                for x in 0..axa[j] {
                    imgbuf[(ALIGHN_X * j) + x] = img_box[x];
                }
                drop(imgbuf); // explicit drop mutex lock
                avaiable_thread.fetch_add(1, Relaxed);
                main_thread.unpark();
            });
        }
    });

    println!(
        "Total Time elapsed: {}",
        now.elapsed().unwrap().as_secs_f64()
    );

    let mut output_img = RgbaImage::new(OUTPUT_X as u32, OUTPUT_Y as u32);
    const SPP_DIV: f64 = 1.0 / (NS as u32).pow(2) as f64;
    let imgbuf = imgbuf.into_inner().unwrap();
    for x in 0..OUTPUT_X {
        for y in 0..OUTPUT_Y {
            let mut accum_pixel: [f64; 3] = [0.0; 3];
            for i in 0..NS {
                for j in 0..NS {
                    accum_pixel = [
                        imgbuf[(x * NS) + i][(y * NS) + j][0] + accum_pixel[0],
                        imgbuf[(x * NS) + i][(y * NS) + j][1] + accum_pixel[1],
                        imgbuf[(x * NS) + i][(y * NS) + j][2] + accum_pixel[2],
                    ];
                }
            }
            let pixel = [
                (accum_pixel[0] * SPP_DIV).sqrt(),
                (accum_pixel[1] * SPP_DIV).sqrt(),
                (accum_pixel[2] * SPP_DIV).sqrt(),
            ];
            let ir: u8 = (255.99 * pixel[0]) as u8;
            let ig: u8 = (255.99 * pixel[1]) as u8;
            let ib: u8 = (255.99 * pixel[2]) as u8;
            output_img.put_pixel(
                x as u32,
                (OUTPUT_Y - (y + 1)) as u32,
                Rgba([ir, ig, ib, 255]),
            );
        }
    }
    output_img.save("image.png").unwrap()
}
