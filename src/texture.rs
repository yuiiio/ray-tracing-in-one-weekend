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
    teximage: Vec<Vec<Vector3<f64>>>,
    width: f64, // usize as f64 for uv calc
    height: f64,
}

impl ImageTexture {
    pub fn new(teximage: RgbaImage) -> Self {
        let width = teximage.width();
        let height = teximage.height();
        let mut pixels: Vec<Vec<Vector3<f64>>> = Vec::with_capacity(height as usize);
        for i in 0..height {
            let mut line: Vec<Vector3<f64>> = Vec::with_capacity(width as usize);
            for j in 0..width {
                let pixel = teximage.get_pixel(j, i);
                line.push([
                          pixel[0] as f64 / 255.99,
                          pixel[1] as f64 / 255.99,
                          pixel[2] as f64 / 255.99,
                ])
            }
            pixels.push(line);
        }
        ImageTexture { teximage: pixels, width: width as f64, height: height as f64 }
    }

    fn get_value(&self, uv: (f64, f64), _p: &Vector3<f64>) -> Vector3<f64> {
        let (u, v) = (uv.0, uv.1);
        let x = min(self.width * u, self.width - 1.0);
        let y = min(self.height * (1.0 - v), self.height - 1.0);
        self.teximage[y as usize][x as usize]
    }
}
