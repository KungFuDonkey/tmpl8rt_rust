use std::path;
use crate::camera::Camera;
use crate::opencl::*;
use crate::math::*;
use crate::scene::GPUScene;
use crate::surface::{SCRHEIGHT, SCRWIDTH};
use crate::blue_noise::load_blue_noise_from_file;

pub struct GPURenderer
{
    generate_rays_kernel: OpenCLKernel,
    extend_kernel: OpenCLKernel,
    shade_kernel: OpenCLKernel,
    connect_kernel: OpenCLKernel,
    finalize_kernel: OpenCLKernel,
    num_rays: OpenCLBuffer<u32>,
    write_back_ids: OpenCLBuffer<u32>,
    ts: OpenCLBuffer<f32>,
    origins: OpenCLBuffer<Float3>,
    directions: OpenCLBuffer<Float3>,
    normals: OpenCLBuffer<Float3>,
    materials: OpenCLBuffer<u64>,
    ray_colors: OpenCLBuffer<Float3>,
    ray_lights: OpenCLBuffer<Float3>,
    accumulator: OpenCLBuffer<Float3>,
    shadow_write_back_ids: OpenCLBuffer<u32>,
    shadow_ray_ts: OpenCLBuffer<f32>,
    shadow_ray_origins: OpenCLBuffer<Float3>,
    shadow_ray_directions: OpenCLBuffer<Float3>,
    shadow_light_colors: OpenCLBuffer<Float3>,
    blue_noise_texture: OpenCLBuffer<u8>,
    seed: u32,
    rendered_frames: u32,
    num_primary_rays: usize,
    num_bounces: u32,
    pub output_buffer: OpenCLBuffer<u32>,
}


impl GPURenderer
{
    pub fn new(cl: &OpenCL) -> Self
    {
        let generate_rays_program = OpenCLProgram::from_file(cl, path::Path::new("./src/kernels/generate_rays.cl"));
        let generate_rays_kernel = OpenCLKernel::from_program(cl, &generate_rays_program, "generate_rays");

        let extend_program = OpenCLProgram::from_file(cl, path::Path::new("./src/kernels/extend.cl"));
        let extend_kernel = OpenCLKernel::from_program(cl, &extend_program, "extend");

        let shade_program = OpenCLProgram::from_file(cl, path::Path::new("./src/kernels/shade.cl"));
        let shade_kernel = OpenCLKernel::from_program(cl, &shade_program, "shade");

        let connect_program = OpenCLProgram::from_file(cl, path::Path::new("./src/kernels/connect.cl"));
        let connect_kernel = OpenCLKernel::from_program(cl, &connect_program, "connect");

        let finalize_program = OpenCLProgram::from_file(cl, path::Path::new("./src/kernels/finalize.cl"));
        let finalize_kernel = OpenCLKernel::from_program(cl, &finalize_program, "finalize");

        let num_primary_rays = SCRWIDTH * SCRHEIGHT;
        let num_bounces = 10;

        let mut num_rays: Vec<u32> = Vec::with_capacity(num_bounces * 2);
        for _ in 0..(num_bounces + 1)
        {
            num_rays.push(0);
            num_rays.push(0);
        }
        num_rays[0] = num_primary_rays as u32; // set initial rays

        let mut write_back_ids: Vec<u32> = Vec::with_capacity(num_primary_rays);
        let mut ts: Vec<f32> = Vec::with_capacity(num_primary_rays);
        let mut origins: Vec<Float3> = Vec::with_capacity(num_primary_rays);
        let mut directions: Vec<Float3> = Vec::with_capacity(num_primary_rays);
        let mut normals: Vec<Float3> = Vec::with_capacity(num_primary_rays);
        let mut output_buffer: Vec<u32> = Vec::with_capacity(num_primary_rays);
        let mut materials: Vec<u64> = Vec::with_capacity(num_primary_rays);
        let mut accumulator: Vec<Float3> = Vec::with_capacity(num_primary_rays);
        let mut ray_colors: Vec<Float3> = Vec::with_capacity(num_primary_rays);
        let mut ray_lights: Vec<Float3> = Vec::with_capacity(num_primary_rays);
        let mut shadow_write_back_ids: Vec<u32> = Vec::with_capacity(num_primary_rays);
        let mut shadow_ray_ts: Vec<f32> = Vec::with_capacity(num_primary_rays);
        let mut shadow_ray_origins: Vec<Float3> = Vec::with_capacity(num_primary_rays);
        let mut shadow_ray_directions: Vec<Float3> = Vec::with_capacity(num_primary_rays);
        let mut shadow_light_colors: Vec<Float3> = Vec::with_capacity(num_primary_rays);

        for _ in 0..num_primary_rays
        {
            write_back_ids.push(0);
            ts.push(0.0);
            origins.push(Float3::zero());
            directions.push(Float3::zero());
            normals.push(Float3::zero());
            output_buffer.push(0);
            materials.push(0);
            accumulator.push(Float3::zero());
            ray_colors.push(Float3::zero());
            ray_lights.push(Float3::zero());
            shadow_write_back_ids.push(0);
            shadow_ray_ts.push(0.0);
            shadow_ray_origins.push(Float3::zero());
            shadow_ray_directions.push(Float3::zero());
            shadow_light_colors.push(Float3::zero());
        }

        let num_rays = OpenCLBuffer::read_write(cl, num_rays);
        let write_back_ids = OpenCLBuffer::read_write(cl, write_back_ids);
        let ts = OpenCLBuffer::read_write(cl, ts);
        let origins = OpenCLBuffer::read_write(cl, origins);
        let directions = OpenCLBuffer::read_write(cl, directions);
        let normals = OpenCLBuffer::read_write(cl, normals);
        let materials = OpenCLBuffer::read_write(cl, materials);
        let output_buffer = OpenCLBuffer::write_only(cl, output_buffer);
        let accumulator = OpenCLBuffer::read_write(cl, accumulator);
        let ray_colors = OpenCLBuffer::read_write(cl, ray_colors);
        let ray_lights = OpenCLBuffer::read_write(cl, ray_lights);
        let shadow_write_back_ids = OpenCLBuffer::read_write(cl, shadow_write_back_ids);
        let shadow_ray_ts = OpenCLBuffer::read_write(cl, shadow_ray_ts);
        let shadow_ray_origins = OpenCLBuffer::read_write(cl, shadow_ray_origins);
        let shadow_ray_directions = OpenCLBuffer::read_write(cl, shadow_ray_directions);
        let shadow_light_colors = OpenCLBuffer::read_write(cl, shadow_light_colors);
        let blue_noise_texture = load_blue_noise_from_file(cl, std::path::PathBuf::from("./assets/blue_noise.png"));

        generate_rays_kernel.set_argument(1, SCRWIDTH as u32);
        generate_rays_kernel.set_argument(2, SCRHEIGHT as u32);
        generate_rays_kernel.set_argument(7, &write_back_ids);
        generate_rays_kernel.set_argument(8, &ts);
        generate_rays_kernel.set_argument(9, &origins);
        generate_rays_kernel.set_argument(10, &directions);
        generate_rays_kernel.set_argument(11, &normals);
        generate_rays_kernel.set_argument(12, &materials);
        generate_rays_kernel.set_argument(13, &ray_colors);
        generate_rays_kernel.set_argument(14, &ray_lights);

        extend_kernel.set_argument(0, num_bounces as u32);
        extend_kernel.set_argument(1, &num_rays);
        extend_kernel.set_argument(2, &ts);
        extend_kernel.set_argument(3, &origins);
        extend_kernel.set_argument(4, &directions);
        extend_kernel.set_argument(5, &normals);
        extend_kernel.set_argument(6, &materials);

        shade_kernel.set_argument(3, SCRWIDTH as u32);
        shade_kernel.set_argument(4, &blue_noise_texture);
        shade_kernel.set_argument(5, &num_rays);
        shade_kernel.set_argument(6, &write_back_ids);
        shade_kernel.set_argument(7, &ts);
        shade_kernel.set_argument(8, &origins);
        shade_kernel.set_argument(9, &directions);
        shade_kernel.set_argument(10, &normals);
        shade_kernel.set_argument(11, &materials);
        shade_kernel.set_argument(12, &ray_colors);
        shade_kernel.set_argument(13, &ray_lights);
        shade_kernel.set_argument(14, &shadow_write_back_ids);
        shade_kernel.set_argument(15, &shadow_ray_ts);
        shade_kernel.set_argument(16, &shadow_ray_origins);
        shade_kernel.set_argument(17, &shadow_ray_directions);
        shade_kernel.set_argument(18, &shadow_light_colors);

        connect_kernel.set_argument(0, num_bounces as u32);
        connect_kernel.set_argument(1, &num_rays);
        connect_kernel.set_argument(2, &shadow_write_back_ids);
        connect_kernel.set_argument(3, &ray_lights);
        connect_kernel.set_argument(4, &shadow_ray_ts);
        connect_kernel.set_argument(5, &shadow_ray_origins);
        connect_kernel.set_argument(6, &shadow_ray_directions);
        connect_kernel.set_argument(7, &shadow_light_colors);

        finalize_kernel.set_argument(1, &ray_lights);
        finalize_kernel.set_argument(2, &accumulator);
        finalize_kernel.set_argument(3, &output_buffer);

        let seed = init_seed(random_uint());
        println!("generated kernels");

        GPURenderer
        {
            generate_rays_kernel,
            extend_kernel,
            shade_kernel,
            connect_kernel,
            finalize_kernel,
            write_back_ids,
            ts,
            origins,
            directions,
            normals,
            materials,
            accumulator,
            output_buffer,
            num_rays,
            ray_colors,
            ray_lights,
            shadow_write_back_ids,
            shadow_ray_ts,
            shadow_ray_origins,
            shadow_ray_directions,
            shadow_light_colors,
            seed,
            num_primary_rays,
            blue_noise_texture,
            num_bounces: num_bounces as u32,
            rendered_frames: 1
        }
    }

    pub fn set_scene(&mut self, scene: &GPUScene)
    {
        self.extend_kernel.set_argument(7, scene.sphere_positions.host_buffer.len() as u32);
        self.extend_kernel.set_argument(8, &scene.sphere_positions);
        self.extend_kernel.set_argument(9, &scene.sphere_radi);
        self.extend_kernel.set_argument(10, &scene.sphere_materials);
        self.extend_kernel.set_argument(11, scene.plane_normals.host_buffer.len() as u32);
        self.extend_kernel.set_argument(12, &scene.plane_normals);
        self.extend_kernel.set_argument(13, &scene.plane_distances);
        self.extend_kernel.set_argument(14, &scene.plane_materials);
        self.extend_kernel.set_argument(15, scene.quad_sizes.host_buffer.len() as u32);
        self.extend_kernel.set_argument(16, &scene.quad_sizes);
        self.extend_kernel.set_argument(17, &scene.quad_inv_transforms);
        self.extend_kernel.set_argument(18, &scene.quad_materials);
        self.extend_kernel.set_argument(19, scene.mesh_offsets.host_buffer.len() as u32);
        self.extend_kernel.set_argument(20, &scene.mesh_offsets);
        self.extend_kernel.set_argument(21, &scene.mesh_triangle_offsets);
        self.extend_kernel.set_argument(22, &scene.mesh_inv_transforms);
        self.extend_kernel.set_argument(23, &scene.mesh_min_bounds);
        self.extend_kernel.set_argument(24, &scene.mesh_max_bounds);
        self.extend_kernel.set_argument(25, &scene.mesh_tri_counts);
        self.extend_kernel.set_argument(26, &scene.mesh_left_firsts);
        self.extend_kernel.set_argument(27, &scene.mesh_triangles);
        self.extend_kernel.set_argument(28, &scene.mesh_triangle_normals);
        self.extend_kernel.set_argument(29, &scene.mesh_materials);

        self.shade_kernel.set_argument(19, scene.quad_sizes.host_buffer.len() as u32);
        self.shade_kernel.set_argument(20, &scene.quad_sizes);
        self.shade_kernel.set_argument(21, &scene.quad_inv_transforms);
        self.shade_kernel.set_argument(22, &scene.quad_materials);

        self.connect_kernel.set_argument(8, scene.sphere_positions.host_buffer.len() as u32);
        self.connect_kernel.set_argument(9, &scene.sphere_positions);
        self.connect_kernel.set_argument(10, &scene.sphere_radi);
        self.connect_kernel.set_argument(11, scene.quad_sizes.host_buffer.len() as u32);
        self.connect_kernel.set_argument(12, &scene.quad_sizes);
        self.connect_kernel.set_argument(13, &scene.quad_inv_transforms);
        self.connect_kernel.set_argument(14, scene.mesh_offsets.host_buffer.len() as u32);
        self.connect_kernel.set_argument(15, &scene.mesh_offsets);
        self.connect_kernel.set_argument(16, &scene.mesh_triangle_offsets);
        self.connect_kernel.set_argument(17, &scene.mesh_inv_transforms);
        self.connect_kernel.set_argument(18, &scene.mesh_min_bounds);
        self.connect_kernel.set_argument(19, &scene.mesh_max_bounds);
        self.connect_kernel.set_argument(20, &scene.mesh_tri_counts);
        self.connect_kernel.set_argument(21, &scene.mesh_left_firsts);
        self.connect_kernel.set_argument(22, &scene.mesh_triangles);
    }

    pub fn set_camera(&mut self, camera: &Camera)
    {
        self.generate_rays_kernel.set_argument(3, &camera.position);
        self.generate_rays_kernel.set_argument(4, &camera.top_left);
        self.generate_rays_kernel.set_argument(5, &camera.bottom_left);
        self.generate_rays_kernel.set_argument(6, &camera.top_right);
        self.rendered_frames = 1;
    }

    pub fn render(&mut self, cl: &OpenCL, scene: &GPUScene, camera: &Camera)
    {
        self.num_rays.copy_to_device(cl);

        self.generate_rays_kernel.set_argument(0, self.seed);
        random_uint_s(&mut self.seed);

        self.generate_rays_kernel.run2d(cl, SCRWIDTH, SCRHEIGHT);
        self.shade_kernel.set_argument(2, self.rendered_frames);
        self.finalize_kernel.set_argument(0, self.rendered_frames);

        for bounces in 0u32..self.num_bounces
        {
            self.extend_kernel.set_argument(0, bounces);
            self.extend_kernel.run(cl, self.num_primary_rays);
            self.shade_kernel.set_argument(0, self.seed);
            self.shade_kernel.set_argument(1, bounces);
            random_uint_s(&mut self.seed);
            self.shade_kernel.run(cl, self.num_primary_rays);

            self.connect_kernel.set_argument(0, bounces);
            self.connect_kernel.run(cl, self.num_primary_rays);
        }

        self.finalize_kernel.run(cl, self.num_primary_rays);

        cl.flush_queue();

        self.output_buffer.copy_from_device(cl);

        self.rendered_frames += 1;
    }
}