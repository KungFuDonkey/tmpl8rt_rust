use std::path;
use crate::camera::Camera;
use crate::opencl::*;
use crate::math::*;
use crate::scene::GPUScene;
use crate::surface::{SCRHEIGHT, SCRWIDTH};

pub struct GPURenderer
{
    generate_rays_kernel: OpenCLKernel,
    extend_kernel: OpenCLKernel,
    shade_kernel: OpenCLKernel,
    finalize_kernel: OpenCLKernel,
    ts: OpenCLBuffer<f32>,
    origins: OpenCLBuffer<Float3>,
    directions: OpenCLBuffer<Float3>,
    normals: OpenCLBuffer<Float3>,
    intersections: OpenCLBuffer<Float3>,
    obj_ids: OpenCLBuffer<u32>,
    accumulator: OpenCLBuffer<Float3>,
    num_rays: usize,
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

        let finalize_program = OpenCLProgram::from_file(cl, path::Path::new("./src/kernels/finalize.cl"));
        let finalize_kernel = OpenCLKernel::from_program(cl, &finalize_program, "finalize");

        let num_rays = SCRWIDTH * SCRHEIGHT;

        let mut ts: Vec<f32> = Vec::with_capacity(num_rays);
        let mut origins: Vec<Float3> = Vec::with_capacity(num_rays);
        let mut directions: Vec<Float3> = Vec::with_capacity(num_rays);
        let mut normals: Vec<Float3> = Vec::with_capacity(num_rays);
        let mut intersections: Vec<Float3> = Vec::with_capacity(num_rays);
        let mut output_buffer: Vec<u32> = Vec::with_capacity(num_rays);
        let mut obj_ids: Vec<u32> = Vec::with_capacity(num_rays);
        let mut accumulator: Vec<Float3> = Vec::with_capacity(num_rays);

        for _ in 0..num_rays
        {
            ts.push(0.0);
            origins.push(Float3::zero());
            directions.push(Float3::zero());
            normals.push(Float3::zero());
            intersections.push(Float3::zero());
            output_buffer.push(0);
            obj_ids.push(0);
            accumulator.push(Float3::zero());
        }

        let ts = OpenCLBuffer::read_write(cl, ts);
        let origins = OpenCLBuffer::read_write(cl, origins);
        let directions = OpenCLBuffer::read_write(cl, directions);
        let normals = OpenCLBuffer::read_write(cl, normals);
        let intersections = OpenCLBuffer::read_write(cl, intersections);
        let obj_ids = OpenCLBuffer::read_write(cl, obj_ids);
        let output_buffer = OpenCLBuffer::write_only(cl, output_buffer);
        let accumulator = OpenCLBuffer::read_write(cl, accumulator);

        generate_rays_kernel.set_argument(0, SCRWIDTH as u32);
        generate_rays_kernel.set_argument(1, SCRHEIGHT as u32);
        generate_rays_kernel.set_argument(6, &ts);
        generate_rays_kernel.set_argument(7, &origins);
        generate_rays_kernel.set_argument(8, &directions);
        generate_rays_kernel.set_argument(9, &intersections);
        generate_rays_kernel.set_argument(10, &normals);
        generate_rays_kernel.set_argument(11, &obj_ids);

        extend_kernel.set_argument(0, num_rays as u32);
        extend_kernel.set_argument(1, &ts);
        extend_kernel.set_argument(2, &origins);
        extend_kernel.set_argument(3, &directions);
        extend_kernel.set_argument(4, &intersections);
        extend_kernel.set_argument(5, &normals);
        extend_kernel.set_argument(6, &obj_ids);

        shade_kernel.set_argument(0, num_rays as u32);
        shade_kernel.set_argument(1, &ts);
        shade_kernel.set_argument(2, &origins);
        shade_kernel.set_argument(3, &directions);
        shade_kernel.set_argument(4, &intersections);
        shade_kernel.set_argument(5, &normals);
        shade_kernel.set_argument(6, &obj_ids);
        shade_kernel.set_argument(9, &accumulator);

        finalize_kernel.set_argument(0, &accumulator);
        finalize_kernel.set_argument(1, &output_buffer);

        ts.copy_to_device(cl);
        origins.copy_to_device(cl);
        directions.copy_to_device(cl);
        obj_ids.copy_to_device(cl);
        output_buffer.copy_to_device(cl);
        accumulator.copy_to_device(cl);

        println!("generated kernels");

        GPURenderer
        {
            generate_rays_kernel,
            extend_kernel,
            shade_kernel,
            finalize_kernel,
            ts,
            origins,
            directions,
            normals,
            intersections,
            obj_ids,
            accumulator,
            output_buffer,
            num_rays
        }
    }

    pub fn render(&mut self, cl: &OpenCL, scene: &GPUScene, camera: &Camera)
    {
        self.generate_rays_kernel.set_argument(2, &camera.position);
        self.generate_rays_kernel.set_argument(3, &camera.top_left);
        self.generate_rays_kernel.set_argument(4, &camera.bottom_left);
        self.generate_rays_kernel.set_argument(5, &camera.top_right);

        self.generate_rays_kernel.run2d(cl, SCRWIDTH, SCRHEIGHT);

        self.extend_kernel.set_argument(7, scene.sphere_positions.host_buffer.len() as u32);
        self.extend_kernel.set_argument(8, &scene.sphere_positions);
        self.extend_kernel.set_argument(9, &scene.sphere_radi);
        self.extend_kernel.run(cl, self.num_rays);

        self.shade_kernel.set_argument(7, &scene.sphere_positions);
        self.shade_kernel.set_argument(8, &scene.sphere_radi);
        self.shade_kernel.run(cl, self.num_rays);

        self.finalize_kernel.run(cl, self.num_rays);

        cl.flush_queue();

        self.output_buffer.copy_from_device(cl);
    }
}