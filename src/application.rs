use std::ops::{Deref};
use crate::camera::Camera;
use crate::scene::{MeshIntersectionSetting, Scene};
use crate::timer::{FrameTimer,Timer};
use imgui_glfw_rs::imgui::Ui;
use imgui_glfw_rs::imgui::ImString;
use crate::input::Input;
use crate::renderer::{LightingMode, Renderer, RenderMode};


pub struct Application
{
    pub renderer: Renderer,
    camera: Camera,
    scene: Scene,
    is_animating: bool,
    is_rendering: bool,
    animation_time: f32,
    internal_timer: FrameTimer,
    manual_frames_to_render: i32,
    manual_redraw_screen: bool,
    static_frame_timer: Timer
}

impl Application
{
    pub fn new() -> Self
    {
        Application {
            renderer: Renderer::new(),
            camera: Camera::new(),
            scene: Scene::new(),
            is_animating: false,
            is_rendering: true,
            animation_time: 0.0,
            internal_timer: FrameTimer::new(),
            static_frame_timer: Timer::new(),
            manual_frames_to_render: 0,
            manual_redraw_screen: false
        }
    }

    pub fn init(&mut self)
    {
    }

    pub fn tick(&mut self, delta_time: f32, input: &Input)
    {
        self.handle_input(input, delta_time);

        if !self.is_rendering
        {
            return;
        }

        self.internal_timer.reset();
        self.renderer.render(&self.scene, &self.camera);
        self.internal_timer.print_frame_time();

        if self.manual_frames_to_render > 0
        {
            self.is_animating = false;
            self.manual_frames_to_render -= 1;
            if self.manual_frames_to_render == 0
            {
                self.is_rendering = false;
                println!();
                println!("FINISHED RENDER OF STATIC FRAME");
                println!("total time: {} seconds", self.static_frame_timer.elapsed_seconds())
            }
        }

        if self.is_animating
        {
            self.animation_time += delta_time;
            self.scene.set_time(self.animation_time);
        }
    }

    fn handle_input(&mut self, input: &Input, delta_time: f32)
    {
        if input.is_key_pressed(glfw::Key::Space)
        {
            self.is_rendering = !self.is_rendering;
            if self.is_rendering && self.manual_frames_to_render > 0
            {
                if self.manual_redraw_screen
                {
                    self.renderer.reset_accumulator();
                }
                self.static_frame_timer.reset();
                self.renderer.render_mode = RenderMode::Standard;

                println!("STARTED RENDER OF STATIC FRAME");
            }
        }

        if !self.is_rendering
        {
            return;
        }

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
        ui.checkbox(ImString::new("Animate scene").deref(), &mut self.is_animating);

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

        ui.text(ImString::new("").deref());
        ui.slider_int(ImString::new("Static frames to render").deref(), &mut self.manual_frames_to_render, 0, 10000 ).build();

        ui.checkbox(ImString::new("Redraw screen?").deref(), &mut self.manual_redraw_screen);
    }

    pub fn shutdown(&mut self)
    {

    }
}