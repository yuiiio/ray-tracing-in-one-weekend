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

pub struct CheckerTexture<T: Texture> {
    m_odd: T,
    m_even: T,
    m_freq: f64,
}

impl<T: Texture> CheckerTexture<T> {
    pub fn new(m_odd: T, m_even: T, m_freq: f64) -> CheckerTexture<T> {
        CheckerTexture { m_odd, m_even, m_freq }
    }
}

impl<T: Texture> Texture for CheckerTexture<T> {
    fn get_value(&self, u: f64, v: f64, p: &Vector3<f64>) -> Vector3<f64> {
        let sines: f64 = (self.m_freq * p[0]).sin() * (self.m_freq * p[1]).sin() * (self.m_freq * p[2]).sin();
        if sines < 0.0 {
            self.m_even.get_value(u, v, p)
        } else {
            self.m_odd.get_value(u, v, p)
        }
    }
}