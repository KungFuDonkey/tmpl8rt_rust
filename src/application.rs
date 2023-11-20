
use std::ops::{Deref};
use rayon::prelude::*;
use crate::camera::Camera;
use crate::math::*;
use crate::scene::{Ray, Scene};
use crate::surface::*;
use crate::timer::{FrameTimer, Timer};
use imgui_glfw_rs::imgui::Ui;
use imgui_glfw_rs::imgui::ImString;
use crate::input::Input;
use crate::renderer::{Renderer, RenderMode};


pub struct Application
{
    pub renderer: Renderer,
    camera: Camera,
    scene: Scene,
    mouse_pos: Int2,
    animating: bool,
    animation_time: f32,
    internal_timer: FrameTimer
}

impl Application
{
    pub fn new() -> Self
    {
        Application {
            renderer: Renderer::new(),
            camera: Camera::new(),
            scene: Scene::new(),
            mouse_pos: Int2::zero(),
            animating: true,
            animation_time: 0.0,
            internal_timer: FrameTimer::new()
        }
    }

    pub fn init(&mut self)
    {
    }

    pub fn tick(&mut self, delta_time: f32, input: &Input)
    {
        if self.animating
        {
            self.animation_time += delta_time;
            self.scene.set_time(self.animation_time);
        }

        self.internal_timer.reset();

        self.renderer.render(&self.scene, &self.camera);

        self.internal_timer.print_frame_time();

        self.camera.handle_input(&input, delta_time);
    }

    pub fn ui(&mut self, ui: &mut Ui)
    {
        ui.checkbox(ImString::new("Animate scene").deref(), &mut self.animating);
        ui.text(ImString::new("Render mode:").deref());
        ui.radio_button(ImString::new("Ray tracing").deref(), &mut self.renderer.render_mode, RenderMode::Standard);
        ui.radio_button(ImString::new("Normals").deref(), &mut self.renderer.render_mode, RenderMode::Normals);
        ui.radio_button(ImString::new("Distance").deref(), &mut self.renderer.render_mode, RenderMode::Distance);
    }

    pub fn shutdown(&mut self)
    {

    }
}