#![allow(dead_code)]

use crate::material::UVMaterial::{CheckerboardMaterial, LogoMaterial, JaccoMaterial};
use crate::math::{Float2, Float3};
use crate::surface::Surface;

pub enum UVMaterial
{
    CheckerboardMaterial,
    LogoMaterial (Surface),
    JaccoMaterial (Surface)
}

pub enum SimpleMaterial
{
    BlackMaterial,
    LinearColorMaterial (Float3),
    UV (UVMaterial)
}

pub enum Material
{
    Simple (SimpleMaterial),
    FullyReflectiveMaterial (SimpleMaterial),
    ReflectiveMaterial (SimpleMaterial, f32  /* reflectivity */),
    RefractiveMaterial (SimpleMaterial /* material */, f32 /* n1/n2 */, f32 /* absorption rate */),
    RefractiveMaterialWithExit (SimpleMaterial /* material */, f32 /* n1/n2 */, f32 /* absorption rate */),
}

fn get_color_from_pixel(pixel: u32) -> Float3
{
    let r = ((pixel >> 16) & 255) as f32;
    let g = ((pixel >> 8) & 255) as f32;
    let b = (pixel & 255) as f32;
    Float3::from_xyz(r, g, b) * Float3::from_a(1.0 / 255.0)
}

fn get_jacco_color(uv: &Float2, texture: &Surface) -> Float3
{
    let ix = ((uv.x - 4.0) * (512.0 / 7.0)) as i32;
    let iy = ((2.0 - uv.y) * (512.0 / 3.0)) as i32;

    let pixel = texture.pixels[((ix & 511) + (iy & 511) * 512) as usize];
    get_color_from_pixel(pixel)
}

pub fn get_color_from_uv_material(material: &UVMaterial, uv: &Float2) -> Float3
{
    match material
    {
        UVMaterial::CheckerboardMaterial =>
        {
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
        },
        UVMaterial::LogoMaterial(texture) =>
        {
            let ix = ((uv.x + 4.0) * (128.0 / 8.0)) as i32;
            let iy = ((2.0 - uv.y) * (64.0 / 3.0)) as i32;
            let pixel = texture.pixels[((ix & 127) + (iy & 63) * 128) as usize];
            let r = ((pixel >> 16) & 255) as f32;
            let g = ((pixel >> 8) & 255) as f32;
            let b = (pixel & 255) as f32;
            Float3::from_xyz(r,g,b) * Float3::from_a(1.0 / 255.0)
        },
        UVMaterial::JaccoMaterial(texture) =>
        {
            get_jacco_color(&uv, &texture)
        }
    }
}

pub fn get_simple_material(material: &Material) -> &SimpleMaterial
{
    match material
    {
        Material::Simple(simple_material) => simple_material,
        Material::ReflectiveMaterial(simple_material, _) => simple_material,
        Material::FullyReflectiveMaterial(simple_material) => simple_material,
        Material::RefractiveMaterial(simple_material, _, _) => simple_material,
        Material::RefractiveMaterialWithExit(simple_material, _, _) => simple_material,
    }
}

pub fn linear_color_simple(color: Float3) -> SimpleMaterial
{
    SimpleMaterial::LinearColorMaterial(color)
}

pub fn linear_color(color: Float3) -> Material
{
    Material::Simple(linear_color_simple(color))
}

pub fn checkerboard_material_simple() -> SimpleMaterial
{
    SimpleMaterial::UV(CheckerboardMaterial)
}

pub fn checkerboard_material() -> Material
{
    Material::Simple(checkerboard_material_simple())
}

pub fn blue_material_simple() -> SimpleMaterial
{
    SimpleMaterial::UV(JaccoMaterial(Surface::load_from_file(std::path::Path::new("./assets/blue.png"))))
}

pub fn blue_material() -> Material
{
    Material::Simple(blue_material_simple())
}

pub fn red_material_simple() -> SimpleMaterial
{
    SimpleMaterial::UV(JaccoMaterial(Surface::load_from_file(std::path::Path::new("./assets/red.png"))))
}

pub fn red_material() -> Material
{
    Material::Simple(red_material_simple())
}

pub fn logo_material_simple() -> SimpleMaterial
{
    SimpleMaterial::UV(LogoMaterial(Surface::load_from_file(std::path::Path::new("./assets/logo.png"))))
}

pub fn logo_material() -> Material
{
    Material::Simple(logo_material_simple())
}

pub fn reflective_material(material: SimpleMaterial, reflectivity: f32) -> Material
{
    Material::ReflectiveMaterial(material, reflectivity)
}

pub fn refractive_material(material: SimpleMaterial, n1n2: f32, absorption: f32) -> Material
{
    Material::RefractiveMaterial(material, n1n2, absorption)
}

pub fn refractive_material_with_exit(material: SimpleMaterial, n2: f32, absorption: f32) -> Material
{
    Material::RefractiveMaterialWithExit(material, n2, absorption)
}

pub fn fully_reflective_material(material: SimpleMaterial) -> Material
{
    Material::FullyReflectiveMaterial(material)
}