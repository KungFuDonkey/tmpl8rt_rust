use crate::camera::Camera;
use crate::math::*;
use crate::surface::*;
use crate::material::*;
use crate::scene::{Ray, Scene};
use rayon::prelude::*;

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum RenderMode
{
    Standard,
    Normals,
    Distance
}

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum LightingMode
{
    None,
    SoftShadows,
    HardShadows
}

pub struct RaytracingSettings
{
    pub max_bounces: i32,
    pub lighting_mode: LightingMode,
    pub area_sample_size: i32,
}


pub struct Renderer
{
    //accumulator: Vec<Float4>,
    pub random_seeds: Vec<u32>,
    pub render_target: Surface,
    pub render_mode: RenderMode,
    pub ray_tracing_settings: RaytracingSettings,
    seed: u32,
}


impl Renderer
{
    pub fn new() -> Renderer
    {
        let seed_base: u32 = 0x123456;
        let mut seed = init_seed(seed_base);

        Renderer{
            //accumulator: vec![Float4::zero(); SCRWIDTH * SCRHEIGHT],
            random_seeds: vec![0; SCRWIDTH * SCRHEIGHT],
            render_target: Surface::new(),
            render_mode: RenderMode::Standard,
            ray_tracing_settings: RaytracingSettings {
                max_bounces: 1,
                lighting_mode: LightingMode::None,
                area_sample_size: 1,
            },
            seed
        }
    }

    fn render_normals(ray: &mut Ray, scene: &Scene) -> Float3
    {
        scene.intersect_scene(ray);
        if ray.obj_idx == -1
        {
            return Float3::zero();
        }
        let intersection = ray.intersection_point();
        let normal = scene.get_normal(ray, &intersection, &ray.direction);

        return (normal + 1.0) * 0.5;
    }

    fn render_distances(ray: &mut Ray, scene: &Scene) -> Float3
    {
        scene.intersect_scene(ray);
        if ray.obj_idx == -1
        {
            return Float3::zero();
        }
        return Float3::from_xyz(ray.t, ray.t, ray.t) * 0.1;
    }

    fn get_lighting_color(ray: &Ray, scene: &Scene, intersection: &Float3, material: &Material, render_settings: &RaytracingSettings, seed: &mut u32) -> Float3
    {
        let lighting: Float3 = match render_settings.lighting_mode {
            LightingMode::None => Float3::from_a(1.0),
            LightingMode::HardShadows => scene.direct_lighting_hard(&intersection, &scene.get_normal(ray, &intersection, &ray.direction)),
            LightingMode::SoftShadows => scene.direct_lighting_soft(&intersection, &scene.get_normal(ray, &intersection, &ray.direction), render_settings.area_sample_size as usize, seed)
        };
        return lighting * Self::get_color_from_material(ray, scene, &intersection, material);
    }

    fn trace(ray: &mut Ray, scene: &Scene, render_settings: &RaytracingSettings, seed: &mut u32, bounces: i32) -> Float3
    {
        scene.intersect_scene(ray);
        if ray.obj_idx == -1
        {
            return Float3::zero();
        }
        let intersection = ray.intersection_point();
        let material = scene.get_material(ray);

        if bounces > render_settings.max_bounces
        {
            return Renderer::get_lighting_color(ray, scene, &intersection, material, render_settings, seed);
        }

        match material
        {
            Material::ReflectiveMaterial(reflection_material, reflectivity) =>
            {
                let normal = scene.get_normal(ray, &intersection, &ray.direction);
                let reflected_ray_direction = reflect(&ray.direction, &normal);
                let mut new_ray = Ray::directed(intersection + reflected_ray_direction * EPSILON, reflected_ray_direction);
                let material_color = Self::get_color_from_simple_material(ray, scene, &intersection, reflection_material);
                let object_color = Renderer::get_lighting_color(ray, scene, &intersection, material, render_settings, seed);
                return (material_color) * (*reflectivity * Renderer::trace(&mut new_ray, scene, render_settings, seed, bounces + 1) + (1.0 - reflectivity) * object_color);
            }
            Material::FullyReflectiveMaterial(reflection_material) =>
            {
                let normal = scene.get_normal(ray, &intersection, &ray.direction);
                let reflected_ray_direction = reflect(&ray.direction, &normal);
                let mut new_ray = Ray::directed(intersection + reflected_ray_direction * EPSILON, reflected_ray_direction);
                let material_color = Self::get_color_from_simple_material(ray, scene, &intersection, reflection_material);
                return (material_color) * Renderer::trace(&mut new_ray, scene, render_settings, seed, bounces + 1);
            }
            _ =>
            {
                return Renderer::get_lighting_color(ray, scene, &intersection, material, render_settings, seed);
            }
        }
    }

    pub fn render(&mut self, scene: &Scene, camera: &Camera)
    {
        for i in 0..(SCRWIDTH * SCRHEIGHT)
        {
            self.random_seeds[i] = random_uint_s(&mut self.seed);
        }

        match self.render_mode {
            RenderMode::Standard => {
                self.render_target.pixels.par_iter_mut().enumerate().for_each(|(index, value)|
                {
                    let mut seed = self.random_seeds[index];
                    let mut ray = camera.get_primary_ray_indexed(index);
                    let ray_color = Renderer::trace(&mut ray, &scene, &self.ray_tracing_settings, &mut seed, 1);
                    *value = rgbf32_to_rgb8_f3(&ray_color);
                });
            },
            RenderMode::Normals => {
                self.render_target.pixels.par_iter_mut().enumerate().for_each(|(index, value)|
                {
                    let mut ray = camera.get_primary_ray_indexed(index);
                    let ray_color = Renderer::render_normals(&mut ray, &scene);
                    *value = rgbf32_to_rgb8_f3(&ray_color);
                });
            },
            RenderMode::Distance => {
                self.render_target.pixels.par_iter_mut().enumerate().for_each(|(index, value)|
                {
                    let mut ray = camera.get_primary_ray_indexed(index);
                    let ray_color = Renderer::render_distances(&mut ray, &scene);
                    *value = rgbf32_to_rgb8_f3(&ray_color);
                });
            }
        }
    }

    fn get_color_from_material(ray: &Ray, scene: &Scene, intersection: &Float3, material: &Material) -> Float3
    {
        Renderer::get_color_from_simple_material(ray, scene, intersection, get_simple_material(material))
    }

    fn get_color_from_simple_material(ray: &Ray, scene: &Scene, intersection: &Float3, simple_material: &SimpleMaterial) -> Float3
    {
        match simple_material
        {
            SimpleMaterial::LinearColorMaterial(color) => { return *color; }
            SimpleMaterial::UV(uv_material) =>
            {
                let uv = scene.get_uv(ray, &intersection);
                return get_color_from_uv_material(uv_material, &uv);
            }
            SimpleMaterial::BlackMaterial => { return Float3::zero(); }
        }
    }
}