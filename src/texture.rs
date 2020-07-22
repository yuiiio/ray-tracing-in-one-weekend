use crate::vec3::Vector3;

pub trait Texture {
    fn get_value(&self, u: f64, v: f64, p: &Vector3<f64>) -> Vector3<f64>;
}

pub struct ColorTexture {
    m_color: Vector3<f64>
}

impl ColorTexture {
    pub fn new(m_color: Vector3<f64>) -> ColorTexture {
        ColorTexture { m_color }
    }
}

impl Texture for ColorTexture {
    fn get_value(&self, u: f64, v: f64, p: &Vector3<f64>) -> Vector3<f64> {
        return self.m_color
    }
}