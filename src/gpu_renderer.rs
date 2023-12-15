use std::path;
use crate::opencl::*;
use crate::math::*;

pub struct GPURenderer
{
    generate_rays_kernel: OpenCLKernel,
    extend_kernel: OpenCLKernel,
    finalize_kernel: OpenCLKernel
    /*ts: OpenCLBuffer<f32>,
    origins: OpenCLBuffer<Float3>,
    directions: OpenCLBuffer<Float3>,
    output_buffer: OpenCLBuffer<u32>*/
}


impl GPURenderer
{
    pub fn new(cl: &OpenCL) -> Self
    {
        let generate_rays_program = OpenCLProgram::from_file(cl, path::Path::new("./src/kernels/generate_rays.cl"));
        let generate_rays_kernel = OpenCLKernel::from_program(cl, &generate_rays_program, "generate_rays");

        let extend_program = OpenCLProgram::from_file(cl, path::Path::new("./src/kernels/extend.cl"));
        let extend_kernel = OpenCLKernel::from_program(cl, &extend_program, "extend");

        let finalize_program = OpenCLProgram::from_file(cl, path::Path::new("./src/kernels/finalize.cl"));
        let finalize_kernel = OpenCLKernel::from_program(cl, &finalize_program, "finalize");

        println!("generated kernels");

        GPURenderer
        {
            generate_rays_kernel,
            extend_kernel,
            finalize_kernel
        }
    }
}