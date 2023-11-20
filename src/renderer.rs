use crate::camera::Camera;
use crate::math::*;
use crate::surface::*;
use crate::scene::{Ray, Scene};
use rayon::prelude::*;

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum RenderMode
{
    Standard,
    Normals,
    Distance
}

pub struct Renderer
{
    //accumulator: Vec<Float4>,
    pub render_target: Surface,
    pub render_mode: RenderMode,
}

impl Renderer
{
    pub fn new() -> Renderer
    {
        Renderer{
            //accumulator: vec![Float4::zero(); SCRWIDTH * SCRHEIGHT],
            render_target: Surface::new(),
            render_mode: RenderMode::Standard
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
        let normal = scene.get_normal(ray.obj_idx, &ray.obj_type, &intersection, &ray.direction);

        return scene.get_albedo(ray.obj_idx, &ray.obj_type, &intersection);
    }

    pub fn render(&mut self, scene: &Scene, camera: &Camera)
    {
        match self.render_mode {
            RenderMode::Standard => {
                self.render_target.pixels.par_iter_mut().enumerate().for_each(|(index, value)|
                {
                    let mut ray = camera.get_primary_ray_indexed(index);
                    let ray_color = Renderer::trace(&mut ray, &scene);
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