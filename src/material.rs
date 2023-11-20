use crate::math::{Float2, Float3};
use crate::surface::Surface;

pub trait Material
{
    fn get_color(&self, uv: &Float2) -> Float3;
}

pub struct LinearColorMaterial
{
    color: Float3
}

impl Material for LinearColorMaterial
{
    fn get_color(&self, uv: &Float2) -> Float3
    {
        return self.color;
    }
}

impl LinearColorMaterial
{
    pub fn new(color: Float3) -> Self
    {
        LinearColorMaterial
        {
            color
        }
    }
}

pub struct CheckerBoardMaterial
{

}

impl Material for CheckerBoardMaterial
{
    fn get_color(&self, uv: &Float2) -> Float3 {
        let mut ix = (uv.x * 2.0 + 96.01) as i32;
        let mut iz = (uv.y * 2.0 + 96.01) as i32;

        if ix == 98 && iz == 98
        {
            ix = (uv.x * 32.01) as i32;
            iz = (uv.y * 32.01) as i32;
        }
        if ix == 94 && iz == 98
        {
            ix = (uv.x * 64.01) as i32;
            iz = (uv.y * 64.01) as i32;
        }

        if ((ix + iz) & 1) != 0
        {
            return Float3::from_a(1.0);
        }
        return Float3::from_a(0.3);
    }
}

impl CheckerBoardMaterial
{
    pub fn new() -> Self { CheckerBoardMaterial{} }
}

pub struct LogoMaterial
{
    texture: Surface
}

impl Material for LogoMaterial
{
    fn get_color(&self, uv: &Float2) -> Float3
    {
        let ix = ((uv.x + 4.0) * (128.0 / 8.0)) as i32;
        let iy = ((2.0 - uv.y) * (64.0 / 3.0)) as i32;
        let pixel = self.texture.pixels[((ix & 127) + (iy & 63) * 128) as usize];
        let r = ((pixel >> 16) & 255) as f32;
        let g = ((pixel >> 8) & 255) as f32;
        let b = (pixel & 255) as f32;
        Float3::from_xyz(r,g,b) * Float3::from_a(1.0 / 255.0)
    }
}

impl LogoMaterial
{
    pub fn new() -> Self
    {
        LogoMaterial
        {
            texture: Surface::load_from_file(std::path::Path::new("./assets/logo.png"))
        }
    }
}

fn get_color_from_pixel(pixel: u32) -> Float3
{
    let r = ((pixel >> 16) & 255) as f32;
    let g = ((pixel >> 8) & 255) as f32;
    let b = (pixel & 255) as f32;
    Float3::from_xyz(r, g, b) * Float3::from_a(1.0 / 255.0)
}

fn get_red_and_blue_color(uv: &Float2, texture: &Surface) -> Float3
{
    let ix = ((uv.x - 4.0) * (512.0 / 7.0)) as i32;
    let iy = ((2.0 - uv.y) * (512.0 / 3.0)) as i32;

    let pixel = texture.pixels[((ix & 511) + (iy & 511) * 512) as usize];
    get_color_from_pixel(pixel)
}


pub struct RedMaterial
{
    texture: Surface
}

impl Material for RedMaterial
{
    fn get_color(&self, uv: &Float2) -> Float3
    {
        get_red_and_blue_color(uv, &self.texture)
    }
}

impl RedMaterial
{
    pub fn new() -> Self
    {
        RedMaterial
        {
            texture: Surface::load_from_file(std::path::Path::new("./assets/red.png"))
        }
    }
}

pub struct BlueMaterial
{
    texture: Surface
}

impl Material for BlueMaterial
{
    fn get_color(&self, uv: &Float2) -> Float3
    {
        get_red_and_blue_color(uv, &self.texture)
    }
}

impl BlueMaterial
{
    pub fn new() -> Self
    {
        BlueMaterial
        {
            texture: Surface::load_from_file(std::path::Path::new("./assets/blue.png"))
        }
    }
}