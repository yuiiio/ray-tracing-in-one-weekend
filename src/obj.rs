use std::fs::File;
use std::io::prelude::*;

use crate::aabb::Aabb;
use crate::hitable::{HitRecord, Hitable};
use crate::hitablelist::HitableList;
use crate::material::MaterialHandle;
use crate::ray::Ray;
use crate::triangle::Triangle;
use crate::vec3::{vec3_dot, vec3_length_f64, vec3_mul_b, vec3_squared_length, vec3_sub, Vector3};

pub fn obj_loader(file: &mut File) -> HitableList {
    let mut contents = String::new();
    file.read_to_string(&mut contents);

    let vertex_line: Vec<&str> = contents
        .lines()
        .filter(|s| s.get(0..1).unwrap() == "v")
        .collect();

    let mut vertex: Vec<Vector3<f64>> = vec![];

    for s in vertex_line {
        let attr: Vec<&str> = s.split(" ").skip(1).collect(); //remove "v" at front

        let xattr: Vec<&str> = attr.get(0).unwrap().split("e-").collect();
        let x: f64 = xattr.get(0).unwrap().parse().unwrap();
        let xe: f64 = xattr.get(1).unwrap().parse().unwrap();
        let minus_xe = -1.0 * xe;
        let xi = x * minus_xe.exp();

        let yattr: Vec<&str> = attr.get(1).unwrap().split("e-").collect();
        let y: f64 = yattr.get(0).unwrap().parse().unwrap();
        let ye: f64 = yattr.get(1).unwrap().parse().unwrap();
        let minus_ye = -1.0 * ye;
        let yi = y * minus_ye.exp();

        let zattr: Vec<&str> = attr.get(2).unwrap().split("e-").collect();
        let z: f64 = zattr.get(0).unwrap().parse().unwrap();
        let ze: f64 = zattr.get(1).unwrap().parse().unwrap();
        let minus_ze = -1.0 * ze;
        let zi = z * minus_ze.exp();

        //scale
        let xi = xi * 200.0;
        let yi = yi * 200.0;
        let zi = zi * 200.0;

        let v: Vector3<f64> = [xi, yi, zi];
        vertex.push(v);
    }

    let face_line: Vec<&str> = contents
        .lines()
        .filter(|s| s.get(0..1).unwrap() == "f")
        .collect();

    let mut hitablelist = HitableList::new();

    for s in face_line {
        let attr: Vec<&str> = s.split(" ").skip(1).collect(); //remove "f" at front

        let i0: usize = attr.get(0).unwrap().parse().unwrap();
        let i1: usize = attr.get(1).unwrap().parse().unwrap();
        let i2: usize = attr.get(2).unwrap().parse().unwrap();

        let v0 = vertex[i0 - 1];
        let v1 = vertex[i1 - 1];
        let v2 = vertex[i2 - 1];

        let triangle = Triangle::new(v0, v1, v2, MaterialHandle(0)); //0:red, 1:white
        hitablelist.push(triangle);
    }

    return hitablelist;
}
