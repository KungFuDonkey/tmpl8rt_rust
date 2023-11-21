use crate::camera::Camera;
use crate::math::*;
use crate::surface::*;
use crate::scene::{Ray, Scene};
use rayon::prelude::*;

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum RenderMode
{
    Standard,
    HardShadows,
    SoftShadows,
    Normals,
    Distance
}

pub struct Renderer
{
    //accumulator: Vec<Float4>,
    pub random_seeds: Vec<u32>,
    pub render_target: Surface,
    pub render_mode: RenderMode,
    seed: u32
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
        let normal = scene.get_normal(ray.obj_idx, &ray.obj_type, &intersection, &ray.direction);

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

    fn trace(ray: &mut Ray, scene: &Scene) -> Float3
    {
        scene.intersect_scene(ray);
        if ray.obj_idx == -1
        {
            return Float3::zero();
        }
        let intersection = ray.intersection_point();

        return scene.get_albedo(ray.obj_idx, &ray.obj_type, &intersection);
    }

    fn trace_soft_shadow(ray: &mut Ray, scene: &Scene, seed: &mut u32) -> Float3
    {
        scene.intersect_scene(ray);
        if ray.obj_idx == -1
        {
            return Float3::zero();
        }
        let intersection = ray.intersection_point();
        let normal = scene.get_normal(ray.obj_idx, &ray.obj_type, &intersection, &ray.direction);

        return scene.direct_lighting_soft(&intersection, &normal, seed) * scene.get_albedo(ray.obj_idx, &ray.obj_type, &intersection);
    }

    fn trace_hard_shadow(ray: &mut Ray, scene: &Scene) -> Float3
    {
        scene.intersect_scene(ray);
        if ray.obj_idx == -1
        {
            return Float3::zero();
        }
        let intersection = ray.intersection_point();
        let normal = scene.get_normal(ray.obj_idx, &ray.obj_type, &intersection, &ray.direction);

        return scene.direct_lighting_hard(&intersection, &normal) * scene.get_albedo(ray.obj_idx, &ray.obj_type, &intersection);
    }

    pub fn render(&mut self, scene: &Scene, camera: &Camera)
    {
        if self.render_mode == RenderMode::SoftShadows
        {
            for i in 0..(SCRWIDTH * SCRHEIGHT)
            {
                self.random_seeds[i] = random_uint_s(&mut self.seed);
            }
        }

        match self.render_mode {
            RenderMode::Standard => {
                self.render_target.pixels.par_iter_mut().enumerate().for_each(|(index, value)|
                {
                    let mut ray = camera.get_primary_ray_indexed(index);
                    let ray_color = Renderer::trace(&mut ray, &scene);
                    *value = rgbf32_to_rgb8_f3(&ray_color);
                });
            },
            RenderMode::SoftShadows => {
                self.render_target.pixels.par_iter_mut().enumerate().for_each(|(index, value)|
                    {
                        let mut seed = self.random_seeds[index];
                        let mut ray = camera.get_primary_ray_indexed(index);
                        let ray_color = Renderer::trace_soft_shadow(&mut ray, &scene, &mut seed);
                        *value = rgbf32_to_rgb8_f3(&ray_color);
                    });
            },
            RenderMode::HardShadows => {
                self.render_target.pixels.par_iter_mut().enumerate().for_each(|(index, value)|
                    {
                        let mut ray = camera.get_primary_ray_indexed(index);
                        let ray_color = Renderer::trace_hard_shadow(&mut ray, &scene);
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
}