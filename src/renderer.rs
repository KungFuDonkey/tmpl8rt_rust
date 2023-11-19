use std::ffi::{CStr, CString};
use std::ops::{Deref, Index};
use rayon::prelude::*;
use crate::camera::Camera;
use crate::math::*;
use crate::scene::{Ray, Scene};
use crate::screen::*;
use crate::timer::{FrameTimer, Timer};
use imgui_glfw_rs::imgui::Ui;
use imgui_glfw_rs::imgui::ImStr;
use imgui_glfw_rs::imgui::ImString;

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum RenderMode
{
    Standard,
    Normals,
    Distance
}

pub struct Renderer
{
    accumulator: Vec<Float4>,
    pub screen: Screen,
    camera: Camera,
    scene: Scene,
    mouse_pos: Int2,
    animating: bool,
    animation_time: f32,
    render_mode: RenderMode,
    internal_timer: FrameTimer
}

impl Renderer
{
    pub fn new() -> Self
    {
        Renderer {
            accumulator: vec![Float4::zero(); SCRWIDTH * SCRHEIGHT],
            screen: Screen::new(),
            camera: Camera::new(),
            scene: Scene::new(),
            mouse_pos: Int2::zero(),
            animating: true,
            animation_time: 0.0,
            render_mode: RenderMode::Normals,
            internal_timer: FrameTimer::new()
        }
    }

    pub fn init(&mut self)
    {
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

    pub fn tick(&mut self, delta_time: f32)
    {
        if self.animating
        {
            self.animation_time += delta_time;
            self.scene.set_time(self.animation_time);
        }

        self.internal_timer.reset();

        match self.render_mode {
            RenderMode::Standard => {
                self.screen.pixels.par_iter_mut().enumerate().for_each(|(index, value)|
                {
                    let mut ray = self.camera.get_primary_ray_indexed(index);
                    let ray_color = Renderer::trace(&mut ray, &self.scene);
                    *value = rgbf32_to_rgb8_f3(&ray_color);
                });
            },
            RenderMode::Normals => {
                self.screen.pixels.par_iter_mut().enumerate().for_each(|(index, value)|
                {
                    let mut ray = self.camera.get_primary_ray_indexed(index);
                    let ray_color = Renderer::render_normals(&mut ray, &self.scene);
                    *value = rgbf32_to_rgb8_f3(&ray_color);
                });
            },
            RenderMode::Distance => {
                self.screen.pixels.par_iter_mut().enumerate().for_each(|(index, value)|
                {
                    let mut ray = self.camera.get_primary_ray_indexed(index);
                    let ray_color = Renderer::render_distances(&mut ray, &self.scene);
                    *value = rgbf32_to_rgb8_f3(&ray_color);
                });
            }
        }

        self.internal_timer.print_frame_time();
    }

    pub fn ui(&mut self, ui: &mut Ui)
    {
        ui.checkbox(ImString::new("Animate scene").deref(), &mut self.animating);
        ui.text(ImString::new("Render mode:").deref());
        ui.radio_button(ImString::new("Ray tracing").deref(), &mut self.render_mode, RenderMode::Standard);
        ui.radio_button(ImString::new("Normals").deref(), &mut self.render_mode, RenderMode::Normals);
        ui.radio_button(ImString::new("Distance").deref(), &mut self.render_mode, RenderMode::Distance);
    }

    pub fn shutdown(&mut self)
    {

    }

    pub fn mouse_button_down(&mut self, button: i32)
    {
        todo!()
    }

    pub fn mouse_button_up(&mut self, button: i32)
    {
        todo!()
    }

    pub fn mouse_move(&mut self, x: i32, y: i32)
    {
        todo!()
    }

    pub fn mouse_wheel(&mut self, y: f32)
    {
        todo!()
    }

    pub fn key_up(&mut self, key: i32)
    {
        todo!()
    }

    pub fn key_down(&mut self, key: i32)
    {
        todo!()
    }
}