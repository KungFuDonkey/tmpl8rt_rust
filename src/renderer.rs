use crate::camera::Camera;
use crate::math::*;
use crate::surface::*;
use crate::material::*;
use crate::scene::{Scene};
use crate::ray::*;
use rayon::prelude::*;
use crate::objects::mesh::MeshIntersectionSetting;

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum RenderMode
{
    Standard,
    Normals,
    Distance,
    Complexity,
    RelativeComplexity
}

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum LightingMode
{
    None,
    SoftShadows,
    HardShadows
}
#[derive(PartialEq, Copy, Clone, Debug)]
pub enum ComplexityMode
{
    Primary,
    Shadow
}

#[derive(PartialEq, Copy, Clone, Debug)]
pub struct RenderSettings
{
    pub max_bounces: i32,
    pub render_mode: RenderMode,
    pub lighting_mode: LightingMode,
    pub area_sample_size: i32,
    pub mesh_intersection_setting: MeshIntersectionSetting,
    pub max_expected_intersection_tests: i32,
    pub complexity_mode: ComplexityMode
}

pub struct RayInfo
{
    pub primary_ray_steps: u32,
    pub primary_ray_triangle_tests: u32,
    pub shadow_ray_steps: u32,
    pub shadow_triangle_tests: u32,
}

pub struct Renderer
{
    accumulator: Vec<Float3>,
    pub random_seeds: Vec<u32>,
    pub render_target: Surface,
    pub render_settings: RenderSettings,
    pub prev_render_settings: RenderSettings,
    seed: u32,
    accumulated_pixels: u32,
    pub complexity_max: u32,
    pub complexity_min: u32,
    pub avg_complexity: f32,
}


impl Renderer
{
    pub fn new() -> Renderer
    {
        let seed_base: u32 = 0x123456;
        let seed = init_seed(seed_base);

        let render_settings = RenderSettings {
            max_bounces: 1,
            render_mode: RenderMode::Complexity,
            lighting_mode: LightingMode::None,
            area_sample_size: 1,
            mesh_intersection_setting: MeshIntersectionSetting::BvhSpatial128,
            max_expected_intersection_tests: 1000,
            complexity_mode: ComplexityMode::Primary
        };

        Renderer{
            accumulator: vec![Float3::zero(); SCRWIDTH * SCRHEIGHT],
            random_seeds: vec![0; SCRWIDTH * SCRHEIGHT],
            render_target: Surface::new(),
            render_settings: render_settings,
            prev_render_settings: render_settings,
            seed,
            accumulated_pixels: 0,
            complexity_max: 0,
            complexity_min: 0,
            avg_complexity: 0.0
        }
    }

    fn render_normals(ray: &mut Ray, scene: &Scene, render_settings: &RenderSettings) -> Float3
    {
        scene.intersect_scene(ray);
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
        scene.intersect_scene(ray);
        if ray.obj_idx == usize::MAX
        {
            return Float3::zero();
        }
        return Float3::from_xyz(ray.t, ray.t, ray.t) * 0.1;
    }

    fn render_complexity(ray: &mut Ray, scene: &Scene, render_settings: &RenderSettings) -> u32
    {
        scene.intersect_scene(ray);

        return ray.intersection_tests;
    }

    fn get_lighting_color(ray: &Ray, scene: &Scene, intersection: &Float3, material: &Material, render_settings: &RenderSettings, seed: &mut u32) -> (Float3, u32, u32)
    {
        let (lighting, intersection_tests, triangle_intersections): (Float3, u32, u32) = match render_settings.lighting_mode {
            LightingMode::None => (Float3::from_a(1.0), 0, 0),
            LightingMode::HardShadows => scene.direct_lighting_hard(&intersection, &scene.get_normal(ray, &intersection, &ray.direction)),
            LightingMode::SoftShadows => scene.direct_lighting_soft(&intersection, &scene.get_normal(ray, &intersection, &ray.direction), render_settings.area_sample_size as usize, seed)
        };
        return (lighting * Self::get_color_from_material(ray, scene, &intersection, material), intersection_tests, triangle_intersections);
    }

    fn trace(ray: &mut Ray, scene: &Scene, render_settings: &RenderSettings, seed: &mut u32, bounces: i32) -> (Float3, RayInfo)
    {
        scene.intersect_scene(ray);
        let mut ray_info = RayInfo
        {
            primary_ray_steps: ray.intersection_tests,
            primary_ray_triangle_tests: ray.triangle_intersection_tests,
            shadow_ray_steps: 0,
            shadow_triangle_tests: 0
        };
        if ray.obj_idx == usize::MAX
        {
            return (Float3::zero(), ray_info);
        }
        let intersection = ray.intersection_point();
        let material = scene.get_material(ray);

        if bounces > render_settings.max_bounces
        {
            let (color, traversal_steps, triangle_intersections) = Renderer::get_lighting_color(ray, scene, &intersection, material, render_settings, seed);
            ray_info.shadow_ray_steps = traversal_steps;
            ray_info.shadow_triangle_tests = triangle_intersections;
            return (color, ray_info);
        }

        match material
        {
            Material::ReflectiveMaterial(reflection_material, reflectivity) =>
            {
                let normal = scene.get_normal(ray, &intersection, &ray.direction);
                let reflected_ray_direction = reflect(&ray.direction, &normal);
                let mut new_ray = Ray::directed(intersection + reflected_ray_direction * EPSILON, reflected_ray_direction);
                let material_color = Self::get_color_from_simple_material(ray, scene, &intersection, reflection_material);
                let (object_color, traversal_steps, triangle_intersections) = Renderer::get_lighting_color(ray, scene, &intersection, material, render_settings, seed);
                ray_info.shadow_ray_steps = traversal_steps;
                ray_info.shadow_triangle_tests = triangle_intersections;
                let (reflection_color, _) = Renderer::trace(&mut new_ray, scene, render_settings, seed, bounces + 1);
                return ((material_color) * (*reflectivity * reflection_color + (1.0 - reflectivity) * object_color), ray_info);
            }
            Material::FullyReflectiveMaterial(reflection_material) =>
            {
                let normal = scene.get_normal(ray, &intersection, &ray.direction);
                let reflected_ray_direction = reflect(&ray.direction, &normal);
                let mut new_ray = Ray::directed(intersection + reflected_ray_direction * EPSILON, reflected_ray_direction);
                let material_color = Self::get_color_from_simple_material(ray, scene, &intersection, reflection_material);

                let (reflection_color, _) = Renderer::trace(&mut new_ray, scene, render_settings, seed, bounces + 1);

                return ((material_color) * reflection_color, ray_info);
            }
            Material::RefractiveMaterialWithExit(refractive_material, n2, absorption) =>
            {
                // currently broken when inside object

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
                if refract_ray.obj_idx == usize::MAX
                {
                    println!("no int");
                    return (Float3::zero(), ray_info);
                }
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

                let (refraction_color, _) = Renderer::trace(&mut new_ray, scene, render_settings, seed, bounces + 1);

                return ((absorption_vector) * refraction_color, ray_info);
            }
            _ =>
            {
                let (color, traversal_steps, triangle_intersections) = Renderer::get_lighting_color(ray, scene, &intersection, material, render_settings, seed);
                ray_info.shadow_ray_steps = traversal_steps;
                ray_info.shadow_triangle_tests = triangle_intersections;
                return (color, ray_info);
            }
        }
    }

    pub fn render(&mut self, scene: &Scene, camera: &Camera)
    {
        if self.render_settings.render_mode != self.prev_render_settings.render_mode || self.render_settings.lighting_mode != self.prev_render_settings.lighting_mode
        {
            self.reset_accumulator();
        }
        self.prev_render_settings = self.render_settings;

        for i in 0..(SCRWIDTH * SCRHEIGHT)
        {
            self.random_seeds[i] = random_uint_s(&mut self.seed);
        }

        match self.render_settings.render_mode {
            RenderMode::Standard => {
                self.accumulated_pixels += 1;

                self.accumulator.par_iter_mut().enumerate().for_each(|(index, value)|
                {
                    let mut seed = self.random_seeds[index];
                    let mut ray = camera.get_primary_ray_indexed(index);
                    let (color, _) = Renderer::trace(&mut ray, &scene, &self.render_settings, &mut seed, 1);
                    *value += color;
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
            },
            RenderMode::Complexity =>
            {
                self.render_target.pixels.par_iter_mut().enumerate().for_each(|(index, value)|
                {
                    let mut seed = self.random_seeds[index];
                    let mut ray = camera.get_primary_ray_indexed(index);
                    let (_, info) = Renderer::trace(&mut ray, &scene, &self.render_settings, &mut seed, 1);
                    if self.render_settings.complexity_mode == ComplexityMode::Primary
                    {
                        *value = info.primary_ray_triangle_tests + info.primary_ray_steps;
                        return;
                    }
                    if self.render_settings.complexity_mode == ComplexityMode::Shadow
                    {
                        *value = info.shadow_triangle_tests + info.shadow_ray_steps;
                        return;
                    }
                });

                self.complexity_max = *self.render_target.pixels.par_iter().max().unwrap();
                self.complexity_min = *self.render_target.pixels.par_iter().min().unwrap();
                self.avg_complexity = (self.render_target.pixels.par_iter().sum::<u32>() as f32) / ((SCRWIDTH * SCRHEIGHT) as f32);

                let max = self.render_settings.max_expected_intersection_tests;

                self.render_target.pixels.par_iter_mut().for_each(|value|
                {
                    let val = *value;
                    let mut scaled_value = ((val as f32) / (max as f32)).min(1.0);
                    if scaled_value < 0.5
                    {
                        scaled_value *= 2.0;
                        *value = rgbf32_to_rgb8_f3(&Float3::from_xyz(0.0, scaled_value, 1.0 - scaled_value));
                        return;
                    }
                    scaled_value = (scaled_value - 0.5) * 2.0;
                    *value = rgbf32_to_rgb8_f3(&Float3::from_xyz(scaled_value, 1.0 - scaled_value, 0.0));

                });
            }
            RenderMode::RelativeComplexity =>
            {
                self.render_target.pixels.par_iter_mut().enumerate().for_each(|(index, value)|
                {
                    let mut seed = self.random_seeds[index];
                    let mut ray = camera.get_primary_ray_indexed(index);
                    let (_, info) = Renderer::trace(&mut ray, &scene, &self.render_settings, &mut seed, 1);

                    *value = info.primary_ray_steps;
                    if self.render_settings.complexity_mode == ComplexityMode::Primary
                    {
                        *value = info.primary_ray_triangle_tests + info.primary_ray_steps;
                        return;
                    }
                    if self.render_settings.complexity_mode == ComplexityMode::Shadow
                    {
                        *value = info.shadow_triangle_tests + info.shadow_ray_steps;
                        return;
                    }
                });

                let max = *self.render_target.pixels.par_iter().max().unwrap();
                self.complexity_max = max;
                self.complexity_min = *self.render_target.pixels.par_iter().min().unwrap();
                self.avg_complexity = (self.render_target.pixels.par_iter().sum::<u32>() as f32) / ((SCRWIDTH * SCRHEIGHT) as f32);

                self.render_target.pixels.par_iter_mut().for_each(|value|
                {
                    let val = *value;
                    let mut scaled_value = ((val as f32) / (max as f32)).min(1.0);
                    if scaled_value < 0.5
                    {
                        scaled_value *= 2.0;
                        *value = rgbf32_to_rgb8_f3(&Float3::from_xyz(0.0, scaled_value, 1.0 - scaled_value));
                        return;
                    }
                    scaled_value = (scaled_value - 0.5) * 2.0;
                    *value = rgbf32_to_rgb8_f3(&Float3::from_xyz(scaled_value, 1.0 - scaled_value, 0.0));

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

    pub fn export_frame(&self, delta_time: f32)
    {
        let intersection_name = match self.render_settings.mesh_intersection_setting
        {
            MeshIntersectionSetting::Raw => "Raw",
            MeshIntersectionSetting::Bvh4 => "BVH4",
            MeshIntersectionSetting::Bvh128 => "BVH128",
            MeshIntersectionSetting::BvhSpatial4 => "BVHSpatial4",
            MeshIntersectionSetting::BvhSpatial128 => "BVHSpatial128",
            MeshIntersectionSetting::Grid16 => "Grid16",
            MeshIntersectionSetting::Grid32 => "Grid32",
            MeshIntersectionSetting::Grid64 => "Grid64",
            MeshIntersectionSetting::KDTree24 => "KDTree24",
            MeshIntersectionSetting::KDTree8 => "KDTree8",
            MeshIntersectionSetting::KDTree16 => "KDTree16",
        };

        let mut complexity_mode = "";

        if self.render_settings.render_mode == RenderMode::Complexity || self.render_settings.render_mode == RenderMode::RelativeComplexity
        {
            complexity_mode = match self.render_settings.complexity_mode
            {
                ComplexityMode::Primary => "_Primary",
                ComplexityMode::Shadow => "_Shadow"
            };
        }

        let mut complexity_value: String = String::new();
        if self.render_settings.render_mode == RenderMode::Complexity
        {
            complexity_value = self.render_settings.max_expected_intersection_tests.to_string();
        }

        let render_mode = match self.render_settings.render_mode
        {
            RenderMode::Standard => "Standard",
            RenderMode::Complexity => "Complexity",
            RenderMode::RelativeComplexity => "Relative_Complexity",
            RenderMode::Distance => "Distance",
            RenderMode::Normals => "Normals"
        };

        let lighting_mode = match self.render_settings.lighting_mode {
            LightingMode::None => "None",
            LightingMode::HardShadows => "Hard",
            LightingMode::SoftShadows => "Soft"
        };

        let mut name = String::from("./output/");
        name += intersection_name;
        name += "_";
        name += render_mode;
        name += complexity_value.as_str();
        name += complexity_mode;
        name += "_";
        name += lighting_mode;
        name += "_";
        name += delta_time.to_string().as_str();
        name += ".png";

        let mut img = image::RgbaImage::new(SCRWIDTH as u32, SCRHEIGHT as u32);

        img.pixels_mut().enumerate().for_each(|(index, pixel)|{
            let value = self.render_target.pixels[index];
            let r = value >> 16;
            let g = (value & 0x0000ff00) >> 8;
            let b = value & 0x000000ff;
            pixel[0] = r as u8;
            pixel[1] = g as u8;
            pixel[2] = b as u8;
            pixel[3] = 255;
        });

        img.save(name).expect("Failed to save file to disk");
    }
}