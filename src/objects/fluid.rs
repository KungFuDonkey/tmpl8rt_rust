
use crate::math::*;
use crate::objects::aabb::{AABB, intersect_aabb};
use crate::opencl::*;

const MAX_PARTICLES : usize = 4096;

pub struct FluidSystem
{
    pub min_bounds: OpenCLBuffer<Float3>,
    pub max_bounds: OpenCLBuffer<Float3>,
    pub particle_counts: OpenCLBuffer<u32>,
    pub particle_offsets: OpenCLBuffer<u32>,
    pub particle_ids: OpenCLBuffer<u32>,
    pub particle_mass: OpenCLBuffer<f32>,
    pub particle_pressures: OpenCLBuffer<f32>,
    pub particle_stiffness: OpenCLBuffer<f32>,
    pub particle_rest_densities: OpenCLBuffer<f32>,
    pub particle_densities: OpenCLBuffer<f32>,
    pub particle_viscosities: OpenCLBuffer<f32>,
    pub particle_positions: OpenCLBuffer<Float3>,
    pub particle_old_positions: OpenCLBuffer<Float3>,
    pub particle_velocities: OpenCLBuffer<Float3>,
    pub particle_accelerations: OpenCLBuffer<Float3>,
    pub particle_rest_coefficients: OpenCLBuffer<f32>,
    pub particle_support_radi: OpenCLBuffer<f32>,
    pub particle_thresholds: OpenCLBuffer<f32>,
    pub particle_surface_tensions: OpenCLBuffer<f32>,
    pub particle_colors: OpenCLBuffer<Float3>,
    pub particle_opacities: OpenCLBuffer<f32>
}

impl FluidSystem
{

    pub fn new(cl: &OpenCL) -> FluidSystem
    {
        let initial_particles: u32 = 512;

        let mut particle_offsets: Vec<u32> = Vec::with_capacity(1);
        particle_offsets.push(0);
        let mut min_bounds: Vec<Float3> = Vec::with_capacity(1);
        min_bounds.push(Float3::from_a(-1.0));
        let mut max_bounds: Vec<Float3> = Vec::with_capacity(1);
        max_bounds.push(Float3::from_a(1.0));
        let mut particle_counts: Vec<u32> = Vec::with_capacity(1);
        particle_counts.push(initial_particles);

        let mut particle_ids : Vec<u32> = Vec::with_capacity(MAX_PARTICLES);
        let mut particle_mass : Vec<f32> = Vec::with_capacity(MAX_PARTICLES);
        let mut particle_pressure : Vec<f32> = Vec::with_capacity(MAX_PARTICLES);
        let mut particle_stiffness : Vec<f32> = Vec::with_capacity(MAX_PARTICLES);
        let mut particle_rest_density : Vec<f32> = Vec::with_capacity(MAX_PARTICLES);
        let mut particle_density : Vec<f32> = Vec::with_capacity(MAX_PARTICLES);
        let mut particle_viscosity : Vec<f32> = Vec::with_capacity(MAX_PARTICLES);
        let mut particle_position : Vec<Float3> = Vec::with_capacity(MAX_PARTICLES);
        let mut particle_old_position : Vec<Float3> = Vec::with_capacity(MAX_PARTICLES);
        let mut particle_velocity : Vec<Float3> = Vec::with_capacity(MAX_PARTICLES);
        let mut particle_acceleration : Vec<Float3> = Vec::with_capacity(MAX_PARTICLES);
        let mut particle_rest_coeff : Vec<f32> = Vec::with_capacity(MAX_PARTICLES);
        let mut particle_support_radius : Vec<f32> = Vec::with_capacity(MAX_PARTICLES);
        let mut particle_threshold : Vec<f32> = Vec::with_capacity(MAX_PARTICLES);
        let mut particle_surface_tension : Vec<f32> = Vec::with_capacity(MAX_PARTICLES);
        let mut particle_color : Vec<Float3> = Vec::with_capacity(MAX_PARTICLES);
        let mut particle_opacity : Vec<f32> = Vec::with_capacity(MAX_PARTICLES);

        let mut seed: u32 = 934845939;
        for _ in 0..initial_particles
        {
            particle_ids.push(particle_mass.len() as u32);
            particle_mass.push(0.02);
            particle_pressure.push(0.0);
            particle_stiffness.push(3.0);
            particle_rest_density.push(998.29);
            particle_density.push(0.0);
            particle_viscosity.push(3.5);
            let x: f32 = 0.2 * (random_float_s(&mut seed) * 2.0 - 1.0);
            let y: f32 = 0.2 * (random_float_s(&mut seed) * 2.0 - 1.0) - 0.5;
            let z: f32 = 0.2 * (random_float_s(&mut seed) * 2.0 - 1.0);
            particle_position.push(Float3::from_xyz(x, y, z));
            particle_old_position.push(Float3::from_xyz(x, y, z));
            particle_support_radius.push(0.457); // should be the same as grid box size?
            particle_velocity.push(Float3::zero());
            particle_opacity.push(0.8);
            particle_threshold.push(7.065);
            particle_surface_tension.push(0.0728);
            particle_color.push(Float3::from_xyz(0.0, 0.0, 0.8));
            particle_rest_coeff.push(0.9);
            particle_acceleration.push(Float3::zero());
        }

        let min_bounds = OpenCLBuffer::read_write(cl, min_bounds);
        let max_bounds = OpenCLBuffer::read_write(cl, max_bounds);
        let particle_counts = OpenCLBuffer::read_write(cl, particle_counts);
        let particle_offsets = OpenCLBuffer::read_write(cl, particle_offsets);
        let particle_ids = OpenCLBuffer::read_write(cl, particle_ids);
        let particle_mass = OpenCLBuffer::read_write(cl, particle_mass);
        let particle_pressures = OpenCLBuffer::read_write(cl, particle_pressure);
        let particle_stiffness = OpenCLBuffer::read_write(cl, particle_stiffness);
        let particle_rest_densities = OpenCLBuffer::read_write(cl, particle_rest_density);
        let particle_densities = OpenCLBuffer::read_write(cl, particle_density);
        let particle_viscosities = OpenCLBuffer::read_write(cl, particle_viscosity);
        let particle_positions = OpenCLBuffer::read_write(cl, particle_position);
        let particle_old_positions = OpenCLBuffer::read_write(cl, particle_old_position);
        let particle_velocities = OpenCLBuffer::read_write(cl, particle_velocity);
        let particle_accelerations = OpenCLBuffer::read_write(cl, particle_acceleration);
        let particle_rest_coefficients = OpenCLBuffer::read_write(cl, particle_rest_coeff);
        let particle_support_radi = OpenCLBuffer::read_write(cl, particle_support_radius);
        let particle_thresholds = OpenCLBuffer::read_write(cl, particle_threshold);
        let particle_surface_tensions = OpenCLBuffer::read_write(cl, particle_surface_tension);
        let particle_colors = OpenCLBuffer::read_write(cl, particle_color);
        let particle_opacities = OpenCLBuffer::read_write(cl, particle_opacity);

        min_bounds.copy_to_device(cl);
        max_bounds.copy_to_device(cl);
        particle_counts.copy_to_device(cl);
        particle_offsets.copy_to_device(cl);
        particle_ids.copy_to_device(cl);
        particle_mass.copy_to_device(cl);
        particle_pressures.copy_to_device(cl);
        particle_stiffness.copy_to_device(cl);
        particle_rest_densities.copy_to_device(cl);
        particle_densities.copy_to_device(cl);
        particle_viscosities.copy_to_device(cl);
        particle_positions.copy_to_device(cl);
        particle_old_positions.copy_to_device(cl);
        particle_velocities.copy_to_device(cl);
        particle_accelerations.copy_to_device(cl);
        particle_rest_coefficients.copy_to_device(cl);
        particle_support_radi.copy_to_device(cl);
        particle_thresholds.copy_to_device(cl);
        particle_surface_tensions.copy_to_device(cl);
        particle_colors.copy_to_device(cl);
        particle_opacities.copy_to_device(cl);

        let mut fluid_system = FluidSystem
            {
                min_bounds,
                max_bounds,
                particle_counts,
                particle_offsets,
                particle_ids,
                particle_mass,
                particle_pressures,
                particle_stiffness,
                particle_rest_densities,
                particle_densities,
                particle_viscosities,
                particle_positions,
                particle_old_positions,
                particle_velocities,
                particle_accelerations,
                particle_rest_coefficients,
                particle_support_radi,
                particle_thresholds,
                particle_surface_tensions,
                particle_colors,
                particle_opacities,
            };




        return fluid_system;
    }


    // this kinda prevents multiple fluids unless I start pre-allocating memory
    pub fn add_water_particle(&mut self, seed: &mut u32)
    {
        self.particle_ids.host_buffer.push(self.particle_mass.host_buffer.len() as u32);
        self.particle_mass.host_buffer.push(0.02);
        self.particle_pressures.host_buffer.push(0.0);
        self.particle_stiffness.host_buffer.push(3.0);
        self.particle_rest_densities.host_buffer.push(998.29);
        self.particle_densities.host_buffer.push(0.0);
        self.particle_viscosities.host_buffer.push(3.5);
        let x: f32 = 0.2 * (random_float_s(seed) * 2.0 - 1.0);
        let y: f32 = 0.2 * (random_float_s(seed) * 2.0 - 1.0) - 0.5;
        let z: f32 = 0.2 * (random_float_s(seed) * 2.0 - 1.0);
        self.particle_positions.host_buffer.push(Float3::from_xyz(x, y, z));
        self.particle_old_positions.host_buffer.push(Float3::from_xyz(x, y, z));
        self.particle_support_radi.host_buffer.push(0.457); // should be the same as grid box size?
        self.particle_velocities.host_buffer.push(Float3::zero());
        self.particle_opacities.host_buffer.push(0.8);
        self.particle_thresholds.host_buffer.push(7.065);
        self.particle_surface_tensions.host_buffer.push(0.0728);
        self.particle_colors.host_buffer.push(Float3::from_xyz(0.0, 0.0, 0.8));
        self.particle_rest_coefficients.host_buffer.push(0.9);
        self.particle_accelerations.host_buffer.push(Float3::zero());
    }
}