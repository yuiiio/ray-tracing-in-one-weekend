use crate::utils::min;
use crate::vec3::Vector3;
use image::RgbaImage;

pub enum Texture {
    ColorTexture,
    //CheckerTexture,
    ImageTexture,
}

pub struct TextureHandle {
    texture_type: Texture,
    position: usize, // each type
    pub needs_uv: bool,
}

pub struct TextureList {
    color_texture_list: Vec<ColorTexture>,
    //checker_texture_list: Vec<CheckerTexture>,
    image_texture_list: Vec<ImageTexture>,
}

impl TextureList {
    pub fn new() -> TextureList {
        TextureList {
            color_texture_list: Vec::new(),
            //checker_texture_list: Vec::new(),
            image_texture_list: Vec::new(),
        }
    }

    pub fn add_color_texture(&mut self, texture: ColorTexture) -> TextureHandle {
        let pos = self.color_texture_list.len();
        self.color_texture_list.push(texture);
        TextureHandle {
            texture_type: Texture::ColorTexture,
            position: pos,
            needs_uv: false,
        }
    }

    pub fn add_image_texture(&mut self, texture: ImageTexture) -> TextureHandle {
        let pos = self.image_texture_list.len();
        self.image_texture_list.push(texture);
        TextureHandle {
            texture_type: Texture::ImageTexture,
            position: pos,
            needs_uv: true,
        }
    }

    pub fn get_value(&self, uv: (f64, f64), p: &Vector3<f64>, texture_handle: &TextureHandle) -> Vector3<f64> {
        let texture_pos = texture_handle.position;
        match texture_handle.texture_type {
            Texture::ColorTexture => self.color_texture_list[texture_pos].get_value(),
            Texture::ImageTexture => self.image_texture_list[texture_pos].get_value(uv, p),
        }
    }
}

pub struct ColorTexture {
    m_color: Vector3<f64>,
}

impl ColorTexture {
    pub fn new(m_color: Vector3<f64>) -> ColorTexture {
        ColorTexture { m_color }
    }

    fn get_value(&self) -> Vector3<f64> {
        self.m_color
    }
}

/* disable recursion type Texture now.
pub struct CheckerTexture {
    m_odd: TextureHandle,
    m_even: TextureHandle,
    m_freq: f64,
}

impl CheckerTexture {
    pub fn new(m_odd: TextureHandle, m_even: TextureHandle, m_freq: f64) -> Self {
        CheckerTexture {
            m_odd,
            m_even,
            m_freq,
        }
    }

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
*/

pub struct ImageTexture {
    teximage: RgbaImage,
}

impl ImageTexture {
    pub fn new(teximage: RgbaImage) -> Self {
        ImageTexture { teximage }
    }

    fn get_value(&self, uv: (f64, f64), _p: &Vector3<f64>) -> Vector3<f64> {
        let (u, v) = (uv.0, uv.1);
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
