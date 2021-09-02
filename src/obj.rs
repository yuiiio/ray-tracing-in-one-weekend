use std::fs::File;
use std::io::prelude::*;

use crate::hitablelist::HitableList;
use crate::material::MaterialHandle;
use crate::triangle::Triangle;
use crate::vec3::Vector3;

pub fn obj_loader(file: &mut File) -> HitableList {
    let mut contents = String::new();
    file.read_to_string(&mut contents)
        .expect("Failed to read to string obj-file");

    let vertex_line: Vec<&str> = contents
        .lines()
        .filter(|s| s.get(0..1).unwrap() == "v")
        .collect();

    let mut vertex: Vec<Vector3<f64>> = vec![];

    for s in vertex_line {
        let attr: Vec<&str> = s.split(" ").skip(1).collect(); //remove "v" at front

        let x: f64 = attr.get(0).unwrap().parse().unwrap();
        let y: f64 = attr.get(1).unwrap().parse().unwrap();
        let z: f64 = attr.get(2).unwrap().parse().unwrap();

        //scale
        let x = x * 200.0;
        let y = y * 200.0;
        let z = z * 200.0;

        let v: Vector3<f64> = [x, y, z];
        vertex.push(v);
    }

    let face_line: Vec<&str> = contents
        .lines()
        .filter(|s| s.get(0..1).unwrap() == "f")
        .collect();

    let mut hitablelist = HitableList::new();

    for s in face_line {
        let attr: Vec<&str> = s.split(" ").skip(1).collect(); //remove "f" at front

        let first: Vec<&str> = attr.get(0).unwrap().split("/").collect();
        let i0: usize = first.get(0).unwrap().parse().unwrap();

        let second: Vec<&str> = attr.get(1).unwrap().split("/").collect();
        let i1: usize = second.get(0).unwrap().parse().unwrap();

        let third: Vec<&str> = attr.get(2).unwrap().split("/").collect();
        let i2: usize = third.get(0).unwrap().parse().unwrap();

        let v0 = vertex[i0 - 1];
        let v1 = vertex[i1 - 1];
        let v2 = vertex[i2 - 1];

        let triangle = Triangle::new(v0, v2, v1, MaterialHandle(5)); //0:red, 1:white
        hitablelist.push(triangle);
    }

    return hitablelist;
}
