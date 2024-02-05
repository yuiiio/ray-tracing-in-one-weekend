use std::fs::File;
use std::io::prelude::*;

use crate::hitablelist::HitableList;
use crate::material::MaterialHandle;
use crate::triangle::Triangle;
use crate::vec3::Vector3;

pub fn obj_loader(file: &mut File, material_handle: MaterialHandle, scale: f64) -> HitableList {
    let mut contents = String::new();
    file.read_to_string(&mut contents)
        .expect("Failed to read to string obj-file");

    let vertex_line: Vec<&str> = contents
        .lines()
        .filter(|s| s.get(0..1).unwrap() == "v")
        .collect();

    let mut vertex: Vec<Vector3<f64>> = Vec::with_capacity(vertex_line.len());

    for s in vertex_line {
        let attr: Vec<&str> = s.split(&" ").skip(1).collect(); //remove "v" at front

        let x: f64 = attr[0].parse().unwrap();
        let y: f64 = attr[1].parse().unwrap();
        let z: f64 = attr[2].parse().unwrap();

        //scale
        let x = x * scale;
        let y = y * scale;
        let z = z * scale;

        let v: Vector3<f64> = [x, y, z];
        vertex.push(v);
    }

    let face_line: Vec<&str> = contents
        .lines()
        .filter(|s| s.get(0..1).unwrap() == "f")
        .collect();

    let mut hitablelist = HitableList::with_capacity(face_line.len());

    for s in face_line {
        let attr: Vec<&str> = s.split(&" ").skip(1).collect(); //remove "f" at front

        let first: Vec<&str> = attr[0].split(&"/").collect();
        let i0: usize = first[0].parse().unwrap();

        let second: Vec<&str> = attr[1].split(&"/").collect();
        let i1: usize = second[0].parse().unwrap();

        let third: Vec<&str> = attr[2].split(&"/").collect();
        let i2: usize = third[0].parse().unwrap();

        let v0 = vertex[i0 - 1];
        let v1 = vertex[i1 - 1];
        let v2 = vertex[i2 - 1];

        let triangle = Triangle::new(v0, v2, v1, material_handle.clone()); 
        hitablelist.push(triangle);
    }

    hitablelist
}
