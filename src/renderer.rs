use crate::camera::Camera;
use crate::math::*;
use crate::surface::*;
use crate::material::*;
use crate::scene::{MeshIntersectionSetting, Scene};
use crate::ray::*;
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

pub struct RenderSettings
{
    pub max_bounces: i32,
    pub lighting_mode: LightingMode,
    pub area_sample_size: i32,
    pub mesh_intersection_setting: MeshIntersectionSetting
}


pub struct Renderer
{
    accumulator: Vec<Float3>,
    pub random_seeds: Vec<u32>,
    pub render_target: Surface,
    pub render_mode: RenderMode,
    pub render_settings: RenderSettings,
    seed: u32,
    accumulated_pixels: u32,
    prev_render_mode: RenderMode,
    prev_lighting_mode: LightingMode
}


impl Renderer
{
    pub fn new() -> Renderer
    {
        let seed_base: u32 = 0x123456;
        let mut seed = init_seed(seed_base);

        Renderer{
            accumulator: vec![Float3::zero(); SCRWIDTH * SCRHEIGHT],
            random_seeds: vec![0; SCRWIDTH * SCRHEIGHT],
            render_target: Surface::new(),
            render_mode: RenderMode::Normals,
            render_settings: RenderSettings {
                max_bounces: 1,
                lighting_mode: LightingMode::None,
                area_sample_size: 1,
                mesh_intersection_setting: MeshIntersectionSetting::Grid
            },
            seed,
            accumulated_pixels: 0,
            prev_render_mode: RenderMode::Normals,
            prev_lighting_mode: LightingMode::None
        }
    }

    fn render_normals(ray: &mut Ray, scene: &Scene, render_settings: &RenderSettings) -> Float3
    {
        scene.intersect_scene(ray, &render_settings.mesh_intersection_setting);
        if ray.obj_idx == usize::MAX
        {
            return Float3::zero();
        }
        let intersection = ray.intersection_point();
        let normal = scene.get_normal(ray, &intersection, &ray.direction);

        return (normal + 1.0) * 0.5;
    }

    fn render_distances(ray: &mut Ray, scene: &Scene, render_settings: &RenderSettings) -> Float3
    {
        scene.intersect_scene(ray, &render_settings.mesh_intersection_setting);
        if ray.obj_idx == usize::MAX
        {
            return Float3::zero();
        }
        return Float3::from_xyz(ray.t, ray.t, ray.t) * 0.1;
    }

    fn get_lighting_color(ray: &Ray, scene: &Scene, intersection: &Float3, material: &Material, render_settings: &RenderSettings, seed: &mut u32) -> Float3
    {
        let lighting: Float3 = match render_settings.lighting_mode {
            LightingMode::None => Float3::from_a(1.0),
            LightingMode::HardShadows => scene.direct_lighting_hard(&intersection, &scene.get_normal(ray, &intersection, &ray.direction), &render_settings.mesh_intersection_setting),
            LightingMode::SoftShadows => scene.direct_lighting_soft(&intersection, &scene.get_normal(ray, &intersection, &ray.direction), render_settings.area_sample_size as usize, seed, &render_settings.mesh_intersection_setting)
        };
        return lighting * Self::get_color_from_material(ray, scene, &intersection, material);
    }

    fn trace(ray: &mut Ray, scene: &Scene, render_settings: &RenderSettings, seed: &mut u32, bounces: i32) -> Float3
    {
        scene.intersect_scene(ray, &render_settings.mesh_intersection_setting);
        if ray.obj_idx == usize::MAX
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
            Material::RefractiveMaterialWithExit(refractive_material, n2, absorption) =>
            {
                // assumes no collisions with any other objects, so a full glass ball with no objects inside

                let normal = scene.get_normal(ray, &intersection, &ray.direction);
                let reversed_dir = -ray.direction;
                let theta1 = dot( &normal, &reversed_dir);
                let theta1_2 = theta1 * theta1;
                let n1n2 = 1.0 / n2; // div by 1
                let n1n2_2 = n1n2 * n1n2;

                let k = 1.0 - n1n2_2 * (1.0 - theta1_2);

                if k < 0.0
                {
                    // do reflection
                }

                let refract_direction = (n1n2 * ray.direction) + (normal * (n1n2 * theta1_2 - k.sqrt()));
                let refract_direction = normalize(&refract_direction);
                let mut refract_ray = Ray::directed(intersection + EPSILON * refract_direction, refract_direction);

                // only do self intersection
                scene.intersect_object(&mut refract_ray, ray.obj_idx as usize, ray.obj_type);
                let refract_normal = scene.get_normal(&refract_ray, &intersection, &refract_direction);
                let reversed_refract_dir = -refract_direction;
                let theta1 = dot (&refract_normal, &reversed_refract_dir);
                let theta1_2 = theta1 * theta1;
                let n2n1 = *n2;
                let n2n1_2 = n2n1 * n2n1;

                let k = 1.0 - n2n1_2 * (1.0 - theta1_2);

                let new_direction = (n2n1 * refract_direction) + (refract_normal * (n2n1 * theta1_2 - k.sqrt()));
                let new_direction = normalize(&new_direction);
                let refract_intersection = refract_ray.intersection_point();

                let mut new_ray = Ray::directed(refract_intersection + EPSILON * new_direction, new_direction);
                let material_color = Float3::from_a(1.0) - Self::get_color_from_simple_material(ray, scene, &intersection, refractive_material);

                let absorption_vector = Float3::from_xyz((-material_color.x * refract_ray.t * absorption).exp(), (-material_color.y * refract_ray.t * absorption).exp(), (-material_color.z * refract_ray.t * absorption).exp());

                return (absorption_vector) * Renderer::trace(&mut new_ray, scene, render_settings, seed, bounces + 1);
            }
            _ =>
            {
                return Renderer::get_lighting_color(ray, scene, &intersection, material, render_settings, seed);
            }
        }
    }

    pub fn render(&mut self, scene: &Scene, camera: &Camera)
    {
        if self.prev_render_mode != self.render_mode || self.render_settings.lighting_mode != self.prev_lighting_mode
        {
            self.reset_accumulator();
        }
        self.prev_render_mode = self.render_mode;
        self.prev_lighting_mode = self.render_settings.lighting_mode;

        for i in 0..(SCRWIDTH * SCRHEIGHT)
        {
            self.random_seeds[i] = random_uint_s(&mut self.seed);
        }

        match self.render_mode {
            RenderMode::Standard => {
                self.accumulated_pixels += 1;

                self.accumulator.par_iter_mut().enumerate().for_each(|(index, value)|
                {
                    let mut seed = self.random_seeds[index];
                    let mut ray = camera.get_primary_ray_indexed(index);
                    *value += Renderer::trace(&mut ray, &scene, &self.render_settings, &mut seed, 1);
                });

                self.render_target.pixels.par_iter_mut().enumerate().for_each(|(index, value)|
                {
                    let acc_color = self.accumulator[index] / (self.accumulated_pixels as f32);
                    *value = rgbf32_to_rgb8_f3(&acc_color);
                });
            },
            RenderMode::Normals => {
                self.render_target.pixels.par_iter_mut().enumerate().for_each(|(index, value)|
                {
                    let mut ray = camera.get_primary_ray_indexed(index);
                    let ray_color = Renderer::render_normals(&mut ray, &scene, &self.render_settings);
                    *value = rgbf32_to_rgb8_f3(&ray_color);
                });
            },
            RenderMode::Distance => {
                self.render_target.pixels.par_iter_mut().enumerate().for_each(|(index, value)|
                {
                    let mut ray = camera.get_primary_ray_indexed(index);
                    let ray_color = Renderer::render_distances(&mut ray, &scene, &self.render_settings);
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

    pub fn reset_accumulator(&mut self)
    {
        self.accumulator.par_iter_mut().for_each(|value|{
            *value = Float3::zero();
        });
        self.accumulated_pixels = 0;
    }
}