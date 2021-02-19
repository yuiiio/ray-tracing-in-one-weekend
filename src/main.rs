use image::{open, Rgba, RgbaImage};
use rand::prelude::*;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::SystemTime;

mod aabb;
mod bvh_node;
mod camera;
mod hitable;
mod hitablelist;
mod material;
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

use bvh_node::BvhNode;
use camera::Camera;
use hitable::Hitable;
use hitablelist::HitableList;
use material::{Dielectric, DiffuseLight, Lambertian, MaterialHandle, Materials, Metal, Scatterd};
use pdf::{CosinePdf, HitablePdf, MixturePdf, Pdf};
use ray::Ray;
use rectangle::{AxisType, Boxel, FlipNormals, Rect};
use sphere::Sphere;
use std::f64;
use texture::{CheckerTexture, ColorTexture, ImageTexture};
use translate::{Rotate, Translate};
use triangle::Triangle;
use utils::clamp;
use vec3::{
    vec3_add, vec3_div, vec3_dot, vec3_mul, vec3_mul_b, vec3_squared_length, vec3_sub,
    vec3_unit_vector_f64, Vector3,
};

fn color<T: Hitable, M: Hitable>(
    r: &Ray,
    world: &Arc<T>,
    light_list: &Arc<M>,
    depth: u32,
    material_list: &Materials,
    last_absorabance: Vector3<f64>,
) -> Vector3<f64> {
    if depth < 50 {
        match world.hit(r, 0.00001, 10000.0) {
            Some(rec) => {
                let emitted = material_list.get(rec.get_mat_ptr()).emitted(r, &rec);
                if let Some(mat_rec) = material_list.get(rec.get_mat_ptr()).scatter(r, &rec) {
                    let absorabance = vec3_mul_b(last_absorabance, rec.get_t() * rec.get_t());
                    let absorabance = vec3_div([1.0, 1.0, 1.0], absorabance);
                    let absorabance = [
                        clamp(absorabance[0], 0.0, 1.0),
                        clamp(absorabance[1], 0.0, 1.0),
                        clamp(absorabance[2], 0.0, 1.0),
                    ];
                    let scatterd = mat_rec.get_scatterd();
                    match scatterd {
                        Scatterd::Ray(next_ray) => {
                            return vec3_add(
                                emitted,
                                vec3_mul(
                                    vec3_mul(
                                        mat_rec.get_attenuation(),
                                        color(
                                            next_ray,
                                            world,
                                            light_list,
                                            depth + 1,
                                            material_list,
                                            mat_rec.get_absorabance(),
                                        ),
                                    ),
                                    absorabance,
                                ),
                            );
                        }
                        Scatterd::Pdf(pdf) => {
                            let hitable_pdf = HitablePdf {
                                hitable: light_list,
                            };

                            let mix_pdf = MixturePdf {
                                pdf0: hitable_pdf,
                                pdf1: pdf,
                            }; // mix pdf light and hitable

                            let next_ray = &Ray::new(rec.get_p(), mix_pdf.generate(&rec));
                            let pdf_value = mix_pdf.value(&rec, &next_ray.direction());

                            if pdf_value > 0.0 {
                                let spdf_value = material_list
                                    .get(rec.get_mat_ptr())
                                    .scattering_pdf(&next_ray, &rec);
                                let albedo = vec3_mul_b(mat_rec.get_attenuation(), spdf_value);
                                let nor_pdf_value = 1.0 / pdf_value;
                                return vec3_add(
                                    emitted,
                                    vec3_mul_b(
                                        vec3_mul(
                                            vec3_mul(
                                                albedo,
                                                color(
                                                    next_ray,
                                                    world,
                                                    light_list,
                                                    depth + 1,
                                                    material_list,
                                                    mat_rec.get_absorabance(),
                                                ),
                                            ),
                                            absorabance,
                                        ),
                                        nor_pdf_value,
                                    ),
                                );
                            } else {
                                return emitted;
                            }
                        }
                    };
                }
                return emitted;
            }
            None => {}
        }
    }
    [0.0, 0.0, 0.0]
}

fn main() {
    let now = SystemTime::now();
    const NX: usize = 800;
    const NY: usize = 800;
    let imgbuf = Arc::new(Mutex::new(vec![vec![[0, 0, 0, 255]; NY]; NX]));
    const NS: usize = 20; //anti-aliasing sample-per-pixel
    let mut obj_list = HitableList::new();
    let mut light_list = HitableList::new();
    let mut material_list = Materials::new();

    let red = material_list.add_material(Lambertian::new(ColorTexture::new([0.65, 0.05, 0.05])));
    let white = material_list.add_material(Lambertian::new(ColorTexture::new([0.73, 0.73, 0.73])));
    let green = material_list.add_material(Lambertian::new(ColorTexture::new([0.12, 0.45, 0.15])));
    let light =
        material_list.add_material(DiffuseLight::new(ColorTexture::new([15.0, 15.0, 15.0])));
    let magick = material_list.add_material(Lambertian::new(ImageTexture::new(
        open("./texture.png").unwrap().into_rgba(),
    )));
    let glass = material_list.add_material(Dielectric::new(2.0, [0.01, 0.01, 0.0]));
    let metal = material_list.add_material(Metal::new(0.0, ColorTexture::new([0.8, 0.85, 0.88])));

    obj_list.push(FlipNormals::new(Rect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        AxisType::kYZ,
        green,
    )));
    obj_list.push(Rect::new(0.0, 555.0, 0.0, 555.0, 0.0, AxisType::kYZ, red));

    let light_rect = FlipNormals::new(Rect::new(
        213.0,
        343.0,
        227.0,
        332.0,
        554.0,
        AxisType::kXZ,
        light,
    ));
    obj_list.push(light_rect.clone());

    obj_list.push(FlipNormals::new(Rect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        AxisType::kXZ,
        white,
    )));
    obj_list.push(Rect::new(0.0, 555.0, 0.0, 555.0, 0.0, AxisType::kXZ, white));
    obj_list.push(FlipNormals::new(Rect::new(
        0.0,
        555.0,
        0.0,
        555.0,
        555.0,
        AxisType::kXY,
        magick,
    )));

    obj_list.push(Translate::new(
        Box::new(Rotate::new(
            Box::new(Boxel::new([0.0, 0.0, 0.0], [165.0, 165.0, 165.0], white)),
            [0.0, 1.0, 0.0],
            -18.0,
        )),
        [130.0, 0.0, 65.0],
    ));
    let metal_box = Translate::new(
        Box::new(Rotate::new(
            Box::new(Boxel::new([0.0, 0.0, 0.0], [165.0, 330.0, 165.0], metal)),
            [0.0, 1.0, 0.0],
            15.0,
        )),
        [265.0, 0.00, 295.0],
    );
    obj_list.push(metal_box.clone());

    let glass_sphere = Sphere::new([455.0, 100.0, 100.0], 100.0, glass);
    obj_list.push(glass_sphere.clone());

    let triangle = Triangle::new(
        [200.0, 200.0, 200.0],
        [200.0, 400.0, 200.0],
        [400.0, 200.0, 200.0],
        white,
    );
    obj_list.push(triangle.clone());

    let obj_list = BvhNode::new(&mut obj_list);

    light_list.push(light_rect);
    light_list.push(metal_box);
    light_list.push(glass_sphere);

    let cam = Camera::new(
        [278.0, 278.0, -800.0],
        [278.0, 278.0, 0.0],
        [0.0, 1.0, 0.0],
        40.0,
        (NX / NY) as f64,
    );

    let cam = Arc::new(cam);
    let obj_list = Arc::new(obj_list);
    let light_list = Arc::new(light_list);
    let material_list = Arc::new(material_list);
    let mut handles = vec![];
    for j in 0..NY {
        let imgbuf_clone = Arc::clone(&imgbuf);
        let cam = Arc::clone(&cam);
        let obj_list = Arc::clone(&obj_list);
        let light_list = Arc::clone(&light_list);
        let material_list = Arc::clone(&material_list);
        let handle = thread::spawn(move || {
            for i in 0..NX {
                let mut col = [0.0 as f64; 3];
                let mut rng = rand::thread_rng();
                for _s in 0..NS {
                    let rand_x: f64 = rng.gen();
                    let rand_y: f64 = rng.gen();
                    let u: f64 = (i as f64 + rand_x) / NX as f64;
                    let v: f64 = (j as f64 + rand_y) / NY as f64;
                    let r = cam.get_ray(u, v);
                    col = vec3_add(
                        color(
                            &r,
                            &obj_list,
                            &light_list,
                            0,
                            &material_list,
                            [0.0, 0.0, 0.0],
                        ),
                        col,
                    );
                }
                let c = 1.0 / NS as f64;
                col = vec3_mul_b(col, c);
                col = [col[0].sqrt(), col[1].sqrt(), col[2].sqrt()];
                col = [
                    clamp(col[0], 0.0, 1.0),
                    clamp(col[1], 0.0, 1.0),
                    clamp(col[2], 0.0, 1.0),
                ];
                let ir: u8 = (255.99 * col[0]) as u8;
                let ig: u8 = (255.99 * col[1]) as u8;
                let ib: u8 = (255.99 * col[2]) as u8;
                let mut imgbuf = imgbuf_clone.lock().unwrap();
                imgbuf[i][j] = [ir, ig, ib, 255];
            }
        });
        handles.push(handle);
    }

    for handle in handles {
        handle.join().unwrap();
    }

    println!("Time elapsed: {}", now.elapsed().unwrap().as_secs_f64());

    let mut img = RgbaImage::new(NX as u32, NY as u32);
    for x in 0..NX {
        for y in 0..NY {
            let imgbuf = imgbuf.lock().unwrap();
            img.put_pixel(x as u32, y as u32, Rgba(imgbuf[x][NY - (y + 1)]));
        }
    }
    img.save("image.png").unwrap()
}
