use std::ops::{Deref};
use crate::camera::Camera;
use crate::scene::{MeshIntersectionSetting, Scene};
use crate::timer::{FrameTimer};
use imgui_glfw_rs::imgui::Ui;
use imgui_glfw_rs::imgui::ImString;
use crate::input::Input;
use crate::renderer::{LightingMode, Renderer, RenderMode};


pub struct Application
{
    pub renderer: Renderer,
    camera: Camera,
    scene: Scene,
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
            animating: false,
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

        if self.camera.handle_input(&input, delta_time)
        {
            self.renderer.reset_accumulator();
        }
    }

    pub fn ui(&mut self, ui: &mut Ui)
    {
        ui.text(ImString::new("Render mode:").deref());
        ui.radio_button(ImString::new("Ray tracing").deref(), &mut self.renderer.render_mode, RenderMode::Standard);
        ui.radio_button(ImString::new("Normals").deref(), &mut self.renderer.render_mode, RenderMode::Normals);
        ui.radio_button(ImString::new("Distance").deref(), &mut self.renderer.render_mode, RenderMode::Distance);
        ui.radio_button(ImString::new("Complexity").deref(), &mut self.renderer.render_mode, RenderMode::Complexity);
        ui.radio_button(ImString::new("Complexity - Relative").deref(), &mut self.renderer.render_mode, RenderMode::RelativeComplexity);

        ui.text(ImString::new("").deref());
        ui.text(ImString::new("Acceleration structure:").deref());
        ui.radio_button(ImString::new("None (not recommended)").deref(), &mut self.renderer.render_settings.mesh_intersection_setting, MeshIntersectionSetting::Raw);
        ui.radio_button(ImString::new("Grid").deref(), &mut self.renderer.render_settings.mesh_intersection_setting, MeshIntersectionSetting::Grid);
        ui.radio_button(ImString::new("KDTree").deref(), &mut self.renderer.render_settings.mesh_intersection_setting, MeshIntersectionSetting::KDTree);
        ui.radio_button(ImString::new("BVH4").deref(), &mut self.renderer.render_settings.mesh_intersection_setting, MeshIntersectionSetting::Bvh4);
        ui.radio_button(ImString::new("BVH128").deref(), &mut self.renderer.render_settings.mesh_intersection_setting, MeshIntersectionSetting::Bvh128);
        ui.radio_button(ImString::new("BVHSpatial4").deref(), &mut self.renderer.render_settings.mesh_intersection_setting, MeshIntersectionSetting::BvhSpatial4);
        ui.radio_button(ImString::new("BVHSpatial128").deref(), &mut self.renderer.render_settings.mesh_intersection_setting, MeshIntersectionSetting::BvhSpatial128);

        ui.text(ImString::new("").deref());
        ui.text(ImString::new("Scene settings:").deref());
        ui.checkbox(ImString::new("Animate scene").deref(), &mut self.animating);

        ui.text(ImString::new("").deref());
        ui.text(ImString::new("Raytracing Settings:").deref());
        ui.slider_int(ImString::new("Bounces").deref(), &mut self.renderer.render_settings.max_bounces, 0, 8 ).build();

        ui.text(ImString::new("").deref());
        ui.text(ImString::new("Shadows:").deref());
        ui.radio_button(ImString::new("No shadows").deref(), &mut self.renderer.render_settings.lighting_mode, LightingMode::None);
        ui.radio_button(ImString::new("Ray tracing + Hard Shadows").deref(), &mut self.renderer.render_settings.lighting_mode, LightingMode::HardShadows);
        ui.radio_button(ImString::new("Ray tracing + Soft Shadows").deref(), &mut self.renderer.render_settings.lighting_mode, LightingMode::SoftShadows);

        ui.text(ImString::new("").deref());
        ui.slider_int(ImString::new("Light sample size").deref(), &mut self.renderer.render_settings.area_sample_size, 1, 8 ).build();

    }

    pub fn shutdown(&mut self)
    {

    }
}