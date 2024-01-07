use std::ops::{Deref};
use crate::camera::Camera;
use crate::scene::{GPUScene, Scene};
use crate::timer::{FrameTimer,Timer};
use imgui_glfw_rs::imgui::Ui;
use imgui_glfw_rs::imgui::ImString;
use crate::gpu_renderer::{GPURenderer, GPURenderMode};
use crate::input::Input;
use crate::renderer::{ComplexityMode, LightingMode, Renderer, RenderMode};
use crate::objects::mesh::MeshIntersectionSetting;
use crate::opencl::OpenCL;


pub struct Application
{
    pub renderer: Renderer,
    pub gpu_renderer: GPURenderer,
    cl: OpenCL,
    camera: Camera,
    scene: Scene,
    gpu_scene: GPUScene,
    is_animating: bool,
    is_rendering: bool,
    animation_time: f32,
    internal_timer: FrameTimer,
    manual_frames_to_render: i32,
    manual_redraw_screen: bool,
    static_frame_timer: Timer,
    ms: f32,
    fps: f32,
    rps: f32,
    pub use_gpu: bool
}

impl Application
{
    pub fn new() -> Self
    {
        let renderer = Renderer::new();
        let mut scene = Scene::new();
        scene.change_intersection_setting(&renderer.render_settings.mesh_intersection_setting);

        let camera = Camera::new();

        let cl = OpenCL::init();
        let gpu_scene = GPUScene::from_scene(&cl, &scene);
        let mut gpu_renderer =  GPURenderer::new(&cl);
        gpu_renderer.set_scene(&gpu_scene);
        gpu_renderer.set_camera(&camera);

        Application {
            renderer,
            cl,
            gpu_renderer,
            camera,
            scene,
            gpu_scene,
            is_animating: false,
            is_rendering: true,
            animation_time: 0.0,
            internal_timer: FrameTimer::new(),
            static_frame_timer: Timer::new(),
            manual_frames_to_render: 0,
            manual_redraw_screen: false,
            ms: 0.0,
            fps: 0.0,
            rps: 0.0,
            use_gpu: true
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
        if self.use_gpu
        {
            self.gpu_scene.fluid_system.simulate(&self.cl, delta_time);
            self.gpu_renderer.render(&self.cl, &self.gpu_scene, &self.camera);
        }
        else
        {
            self.renderer.render(&self.scene, &self.camera);
        }

        (self.ms, self.fps, self.rps) = self.internal_timer.get_frame_time();

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
        if self.renderer.render_settings.mesh_intersection_setting != self.renderer.prev_render_settings.mesh_intersection_setting
        {
            self.scene.change_intersection_setting(&self.renderer.render_settings.mesh_intersection_setting);
        }

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
                self.renderer.render_settings.render_mode = RenderMode::Standard;

                println!("STARTED RENDER OF STATIC FRAME");
            }
        }

        if input.is_key_pressed(glfw::Key::P)
        {
            self.gpu_scene.fluid_system.settings.paused = !self.gpu_scene.fluid_system.settings.paused;
        }

        if !self.is_rendering
        {
            return;
        }

        if input.is_key_pressed(glfw::Key::C)
        {
            if self.renderer.render_settings.render_mode == RenderMode::Complexity || self.renderer.render_settings.render_mode == RenderMode::RelativeComplexity
            {
                println!("{} {} {}", self.renderer.complexity_max, self.renderer.complexity_min, self.renderer.avg_complexity);
            }
            self.renderer.export_frame(delta_time);
        }

        if self.camera.handle_input(&input, delta_time)
        {
            if self.use_gpu
            {
                self.gpu_renderer.set_camera(&self.camera);
            }
            else
            {
                self.renderer.reset_accumulator();
            }
        }
    }

    pub fn ui(&mut self, ui: &mut Ui)
    {
        if self.is_rendering
        {
            let render_string = format!("ms : {}\nfps: {}\nrps: {}\n", self.ms, self.fps, self.rps);

            ui.text(ImString::new(render_string).deref());
        }
        else
        {
            ui.text(ImString::new("paused").deref());
        }


        if self.use_gpu
        {
            ui.text(ImString::new("Render mode:").deref());
            ui.radio_button(ImString::new("Debug").deref(), &mut self.gpu_renderer.render_mode, GPURenderMode::Debug);
            ui.radio_button(ImString::new("Path Tracing").deref(), &mut self.gpu_renderer.render_mode, GPURenderMode::PathTracer);

            ui.text(ImString::new("Fluid").deref());
            ui.slider_float(ImString::new("Time Scale").deref(), &mut self.gpu_scene.fluid_system.settings.time_scale, 0.0, 1.0 ).build();
            ui.slider_float(ImString::new("Gravity").deref(), &mut self.gpu_scene.fluid_system.settings.gravity, -10.0, 10.0 ).build();
            ui.slider_float(ImString::new("Collision Damping").deref(), &mut self.gpu_scene.fluid_system.settings.collision_damping, 0.001, 1.0).build();
            ui.slider_float(ImString::new("Smoothing Radius").deref(), &mut self.gpu_scene.fluid_system.settings.smoothing_radius, 0.001, 1.0).build();
            ui.slider_float(ImString::new("Target Density").deref(), &mut self.gpu_scene.fluid_system.settings.target_density, 0.001, 5000.0).build();
            ui.slider_float(ImString::new("Pressure Multiplier").deref(), &mut self.gpu_scene.fluid_system.settings.pressure_multiplier, 0.001, 4.0).build();
            ui.slider_float(ImString::new("Near Pressure Multiplier").deref(), &mut self.gpu_scene.fluid_system.settings.near_pressure_multiplier, 0.001, 4.0).build();
            ui.slider_float(ImString::new("Viscosity strength").deref(), &mut self.gpu_scene.fluid_system.settings.viscosity_strength, 0.001, 0.075).build();

            return;
        }
        ui.text(ImString::new("Render mode:").deref());
        ui.radio_button(ImString::new("Ray tracing").deref(), &mut self.renderer.render_settings.render_mode, RenderMode::Standard);
        ui.radio_button(ImString::new("Normals").deref(), &mut self.renderer.render_settings.render_mode, RenderMode::Normals);
        ui.radio_button(ImString::new("Distance").deref(), &mut self.renderer.render_settings.render_mode, RenderMode::Distance);
        ui.radio_button(ImString::new("Complexity").deref(), &mut self.renderer.render_settings.render_mode, RenderMode::Complexity);
        ui.radio_button(ImString::new("Complexity - Relative").deref(), &mut self.renderer.render_settings.render_mode, RenderMode::RelativeComplexity);

        ui.slider_int(ImString::new("Max intersection Tests").deref(), &mut self.renderer.render_settings.max_expected_intersection_tests, 100, 1000 ).build();

        ui.text(ImString::new("Complexity mode:").deref());
        ui.radio_button(ImString::new("Primary rays").deref(), &mut self.renderer.render_settings.complexity_mode, ComplexityMode::Primary);
        ui.radio_button(ImString::new("Shadow rays").deref(), &mut self.renderer.render_settings.complexity_mode, ComplexityMode::Shadow);

        ui.text(ImString::new("").deref());
        ui.text(ImString::new("Acceleration structure:").deref());
        ui.radio_button(ImString::new("None (not recommended)").deref(), &mut self.renderer.render_settings.mesh_intersection_setting, MeshIntersectionSetting::Raw);
        ui.radio_button(ImString::new("Grid16").deref(), &mut self.renderer.render_settings.mesh_intersection_setting, MeshIntersectionSetting::Grid16);
        ui.radio_button(ImString::new("Grid32").deref(), &mut self.renderer.render_settings.mesh_intersection_setting, MeshIntersectionSetting::Grid32);
        ui.radio_button(ImString::new("Grid64").deref(), &mut self.renderer.render_settings.mesh_intersection_setting, MeshIntersectionSetting::Grid64);
        ui.radio_button(ImString::new("KDTree8").deref(), &mut self.renderer.render_settings.mesh_intersection_setting, MeshIntersectionSetting::KDTree8);
        ui.radio_button(ImString::new("KDTree16").deref(), &mut self.renderer.render_settings.mesh_intersection_setting, MeshIntersectionSetting::KDTree16);
        ui.radio_button(ImString::new("KDTree24").deref(), &mut self.renderer.render_settings.mesh_intersection_setting, MeshIntersectionSetting::KDTree24);
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