use core::f64;
use core::ops::{Add, Div, Mul, Sub};

pub type Vector3<T> = [T; 3];

#[inline(always)]
pub fn vec3_sub<T>(a: &Vector3<T>, b: &Vector3<T>) -> Vector3<T>
where
    T: Copy + Sub<T, Output = T>,
{
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline(always)]
pub fn vec3_add<T>(a: &Vector3<T>, b: &Vector3<T>) -> Vector3<T>
where
    T: Copy + Add<T, Output = T>,
{
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline(always)]
pub fn vec3_mul<T>(a: &Vector3<T>, b: &Vector3<T>) -> Vector3<T>
where
    T: Copy + Mul<T, Output = T>,
{
    [a[0] * b[0], a[1] * b[1], a[2] * b[2]]
}

#[allow(dead_code)]
#[inline(always)]
pub fn vec3_div_f64(a: &Vector3<f64>, b: &Vector3<f64>) -> Vector3<f64> {
    // assert!(b[0] != 0.0 && b[1] != 0.0 && b[2] != 0.0);
    [a[0] / b[0], a[1] / b[1], a[2] / b[2]]
}

#[inline(always)]
pub fn vec3_sub_b<T>(a: &Vector3<T>, b: T) -> Vector3<T>
where
    T: Copy + Sub<T, Output = T>,
{
    [a[0] - b, a[1] - b, a[2] - b]
}

#[inline(always)]
pub fn vec3_add_b<T>(a: &Vector3<T>, b: T) -> Vector3<T>
where
    T: Copy + Add<T, Output = T>,
{
    [a[0] + b, a[1] + b, a[2] + b]
}

#[inline(always)]
pub fn vec3_mul_b<T>(a: &Vector3<T>, b: T) -> Vector3<T>
where
    T: Copy + Mul<T, Output = T>,
{
    [a[0] * b, a[1] * b, a[2] * b]
}

#[inline(always)]
pub fn vec3_dot<T>(a: &Vector3<T>, b: &Vector3<T>) -> T
where
    T: Copy + Add<T, Output = T> + Mul<T, Output = T>,
{
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline(always)]
pub fn vec3_length_f64(a: &Vector3<f64>) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

#[inline(always)]
pub fn vec3_squared_length<T>(a: &Vector3<T>) -> T
where
    T: Copy + Add<T, Output = T> + Mul<T, Output = T>,
{
    a[0] * a[0] + a[1] * a[1] + a[2] * a[2]
}

#[inline(always)]
pub fn vec3_unit_vector_f64(a: &Vector3<f64>) -> Vector3<f64> {
    let b: f64 = vec3_length_f64(a);
    //assert!(b != 0.0);
    let c = 1.0 / b;
    vec3_mul_b(a, c)
}

#[inline(always)]
pub fn cross(a: &Vector3<f64>, b: &Vector3<f64>) -> Vector3<f64> {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
