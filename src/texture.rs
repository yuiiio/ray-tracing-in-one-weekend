use crate::utils::min;
use crate::vec3::Vector3;
use image::RgbaImage;

pub trait Texture {
    fn get_value(&self, u: f64, v: f64, p: &Vector3<f64>) -> Vector3<f64>;
}

pub struct ColorTexture {
    m_color: Vector3<f64>,
}

impl ColorTexture {
    pub fn new(m_color: Vector3<f64>) -> ColorTexture {
        ColorTexture { m_color }
    }
}

impl Texture for ColorTexture {
    fn get_value(&self, _u: f64, _v: f64, _p: &Vector3<f64>) -> Vector3<f64> {
        return self.m_color;
    }
}

pub struct CheckerTexture<T: Texture> {
    m_odd: T,
    m_even: T,
    m_freq: f64,
}

impl<T: Texture> CheckerTexture<T> {
    pub fn new(m_odd: T, m_even: T, m_freq: f64) -> Self {
        CheckerTexture {
            m_odd,
            m_even,
            m_freq,
        }
    }
}

impl<T: Texture> Texture for CheckerTexture<T> {
    fn get_value(&self, u: f64, v: f64, p: &Vector3<f64>) -> Vector3<f64> {
        let sines: f64 =
            (self.m_freq * p[0]).sin() * (self.m_freq * p[1]).sin() * (self.m_freq * p[2]).sin();
        if sines.is_sign_negative() {
            self.m_even.get_value(u, v, p)
        } else {
            self.m_odd.get_value(u, v, p)
        }
    }
}

pub struct ImageTexture {
    teximage: RgbaImage,
}

impl ImageTexture {
    pub fn new(teximage: RgbaImage) -> Self {
        ImageTexture { teximage }
    }
}

impl Texture for ImageTexture {
    fn get_value(&self, u: f64, v: f64, _p: &Vector3<f64>) -> Vector3<f64> {
        let width = self.teximage.width();
        let height = self.teximage.height();
        let x = min(width as f64 * u, width as f64 - 1.0);
        let y = min(height as f64 * (1.0 - v), height as f64 - 1.0);
        let pixel = self.teximage.get_pixel(x as u32, y as u32);
        [
            pixel[0] as f64 / 255.99,
            pixel[1] as f64 / 255.99,
            pixel[2] as f64 / 255.99,
        ]
    }
}
