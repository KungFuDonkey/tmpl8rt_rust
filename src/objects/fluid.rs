use std::path;
use crate::math::*;
use crate::objects::aabb::{AABB, intersect_aabb};
use crate::opencl::*;

const MAX_PARTICLES : usize = 16384;


#[derive(Clone, Copy)]
pub struct FluidSystemSettings
{
    pub paused: bool,
    pub local_to_world: Mat4,
    pub world_to_local: Mat4,
    pub time_scale: f32,
    pub num_particles: u32,
    num_sort_stages: u32,
    pub gravity: f32,
    pub collision_damping: f32,
    pub smoothing_radius: f32,
    pub target_density: f32,
    pub pressure_multiplier: f32,
    pub near_pressure_multiplier: f32,
    pub viscosity_strength: f32,
}

pub struct FluidSystem
{
    sort_spatial_indices_kernel: OpenCLKernel,
    calculate_offsets_kernel: OpenCLKernel,
    external_forces_kernel: OpenCLKernel,
    update_spatial_hash_kernel: OpenCLKernel,
    calculate_densities_kernel: OpenCLKernel,
    calculate_pressure_force_kernel: OpenCLKernel,
    calculate_viscosity_kernel: OpenCLKernel,
    update_positions_kernel: OpenCLKernel,

    pub settings: FluidSystemSettings,
    prev_settings: FluidSystemSettings,

    pub particle_positions: OpenCLBuffer<Float3>,
    pub predicted_particle_positions: OpenCLBuffer<Float3>,
    pub particle_velocities: OpenCLBuffer<Float3>,
    pub particle_densities: OpenCLBuffer<Float2>,
    pub spatial_indices: OpenCLBuffer<Uint3>,
    pub spatial_offsets: OpenCLBuffer<u32>
}

impl FluidSystem
{

    pub fn new(cl: &OpenCL) -> FluidSystem
    {
        let fluid_sim_program = OpenCLProgram::from_file(cl, path::Path::new("./src/kernels/fluid_simulation/simulate.cl"));
        let sort_spatial_indices_kernel = OpenCLKernel::from_program(cl, &fluid_sim_program, "sort_spatial_indices");
        let calculate_offsets_kernel = OpenCLKernel::from_program(cl, &fluid_sim_program, "calculate_offsets");
        let external_forces_kernel = OpenCLKernel::from_program(cl, &fluid_sim_program, "external_forces");
        let update_spatial_hash_kernel = OpenCLKernel::from_program(cl, &fluid_sim_program, "update_spatial_hash");
        let calculate_densities_kernel = OpenCLKernel::from_program(cl, &fluid_sim_program, "calculate_densities");
        let calculate_pressure_force_kernel = OpenCLKernel::from_program(cl, &fluid_sim_program, "calculate_pressure_force");
        let calculate_viscosity_kernel = OpenCLKernel::from_program(cl, &fluid_sim_program, "calculate_viscosity");
        let update_positions_kernel = OpenCLKernel::from_program(cl, &fluid_sim_program, "update_positions");

        let initial_particles: u32 = 4096;

        let mut particle_positions : Vec<Float3> = Vec::with_capacity(MAX_PARTICLES);
        let mut predicted_particle_positions : Vec<Float3> = Vec::with_capacity(MAX_PARTICLES);
        let mut particle_velocities : Vec<Float3> = Vec::with_capacity(MAX_PARTICLES);
        let mut particle_densities : Vec<Float2> = Vec::with_capacity(MAX_PARTICLES);
        let mut spatial_indices : Vec<Uint3> = Vec::with_capacity(MAX_PARTICLES);
        let mut spatial_offsets : Vec<u32> = Vec::with_capacity(MAX_PARTICLES);

        let mut seed: u32 = 934845939;
        for _ in 0..initial_particles
        {
            let x: f32 = 0.2 * (random_float_s(&mut seed) * 2.0 - 1.0);
            let y: f32 = 0.2 * (random_float_s(&mut seed) * 2.0 - 1.0);
            let z: f32 = 0.2 * (random_float_s(&mut seed) * 2.0 - 1.0);
            particle_positions.push(Float3::from_xyz(x, y, z));
            predicted_particle_positions.push(Float3::from_xyz(x, y, z));
            particle_velocities.push(Float3::zero());
            particle_densities.push(Float2::zero());
            spatial_indices.push(Uint3::zero());
            spatial_offsets.push(0);
        }

        let particle_positions = OpenCLBuffer::read_write(cl, particle_positions);
        let predicted_particle_positions = OpenCLBuffer::read_write(cl, predicted_particle_positions);
        let particle_velocities = OpenCLBuffer::read_write(cl, particle_velocities);
        let particle_densities = OpenCLBuffer::read_write(cl, particle_densities);
        let spatial_indices = OpenCLBuffer::read_write(cl, spatial_indices);
        let spatial_offsets = OpenCLBuffer::read_write(cl, spatial_offsets);

        particle_positions.copy_to_device(cl);
        predicted_particle_positions.copy_to_device(cl);
        particle_velocities.copy_to_device(cl);
        particle_densities.copy_to_device(cl);
        spatial_indices.copy_to_device(cl);
        spatial_offsets.copy_to_device(cl);

        let num_particles = initial_particles;
        let gravity: f32 = -9.81;
        let collision_damping: f32 = 0.8;
        let smoothing_radius: f32 = 0.5;
        let target_density: f32 = 900.0;
        let pressure_multiplier: f32 = 0.4;
        let near_pressure_multiplier: f32 = 0.4;
        let viscosity_strength: f32 = 0.001;

        let local_to_world: Mat4 = Mat4::identity_matrix();
        let world_to_local: Mat4 = Mat4::identity_matrix();

        sort_spatial_indices_kernel.set_argument(0, num_particles);
        sort_spatial_indices_kernel.set_argument(4, &spatial_indices);

        calculate_offsets_kernel.set_argument(0, num_particles);
        calculate_offsets_kernel.set_argument(1, &spatial_indices);
        calculate_offsets_kernel.set_argument(2, &spatial_offsets);

        external_forces_kernel.set_argument(1, num_particles);
        external_forces_kernel.set_argument(2, gravity);
        external_forces_kernel.set_argument(3, &particle_positions);
        external_forces_kernel.set_argument(4, &predicted_particle_positions);
        external_forces_kernel.set_argument(5, &particle_velocities);

        update_spatial_hash_kernel.set_argument(0, num_particles);
        update_spatial_hash_kernel.set_argument(1, smoothing_radius);
        update_spatial_hash_kernel.set_argument(2, &predicted_particle_positions);
        update_spatial_hash_kernel.set_argument(3, &spatial_indices);
        update_spatial_hash_kernel.set_argument(4, &spatial_offsets);

        calculate_densities_kernel.set_argument(0, num_particles);
        calculate_densities_kernel.set_argument(1, smoothing_radius);
        calculate_densities_kernel.set_argument(2, &predicted_particle_positions);
        calculate_densities_kernel.set_argument(3, &particle_densities);
        calculate_densities_kernel.set_argument(4, &spatial_indices);
        calculate_densities_kernel.set_argument(5, &spatial_offsets);

        calculate_pressure_force_kernel.set_argument(1, num_particles);
        calculate_pressure_force_kernel.set_argument(2, smoothing_radius);
        calculate_pressure_force_kernel.set_argument(3, target_density);
        calculate_pressure_force_kernel.set_argument(4, pressure_multiplier);
        calculate_pressure_force_kernel.set_argument(5, near_pressure_multiplier);
        calculate_pressure_force_kernel.set_argument(6, &predicted_particle_positions);
        calculate_pressure_force_kernel.set_argument(7, &particle_velocities);
        calculate_pressure_force_kernel.set_argument(8, &particle_densities);
        calculate_pressure_force_kernel.set_argument(9, &spatial_indices);
        calculate_pressure_force_kernel.set_argument(10, &spatial_offsets);

        calculate_viscosity_kernel.set_argument(1, num_particles);
        calculate_viscosity_kernel.set_argument(2, smoothing_radius);
        calculate_viscosity_kernel.set_argument(3, viscosity_strength);
        calculate_viscosity_kernel.set_argument(4, &predicted_particle_positions);
        calculate_viscosity_kernel.set_argument(5, &particle_velocities);
        calculate_viscosity_kernel.set_argument(6, &spatial_indices);
        calculate_viscosity_kernel.set_argument(7, &spatial_offsets);

        update_positions_kernel.set_argument(1, num_particles);
        update_positions_kernel.set_argument(2, collision_damping);
        update_positions_kernel.set_argument(3, &local_to_world);
        update_positions_kernel.set_argument(4, &world_to_local);
        update_positions_kernel.set_argument(5, &particle_positions);
        update_positions_kernel.set_argument(6, &particle_velocities);

        let settings = FluidSystemSettings
        {
            paused: false,
            local_to_world,
            world_to_local,
            time_scale: 1.0,
            num_particles,
            num_sort_stages: FluidSystem::calc_num_sort_stages(num_particles),
            gravity,
            collision_damping,
            smoothing_radius,
            target_density,
            pressure_multiplier,
            near_pressure_multiplier,
            viscosity_strength
        };

        let mut fluid_system = FluidSystem
        {
            sort_spatial_indices_kernel,
            calculate_offsets_kernel,
            external_forces_kernel,
            update_spatial_hash_kernel,
            calculate_densities_kernel,
            calculate_pressure_force_kernel,
            calculate_viscosity_kernel,
            update_positions_kernel,
            settings,
            prev_settings: settings,
            particle_positions,
            predicted_particle_positions,
            particle_velocities,
            particle_densities,
            spatial_indices,
            spatial_offsets
        };

        return fluid_system;
    }

    fn calc_num_sort_stages(num_particles: u32) -> u32
    {
        for i in 0u32..32u32
        {
            if (num_particles >> i) == 0
            {
                return i;
            }
        }
        return 0;
    }

    fn update_settings(&mut self)
    {
        let mut changed = false;
        if self.settings.local_to_world != self.prev_settings.local_to_world
        {
            changed = true;
            self.update_positions_kernel.set_argument(3, &self.settings.local_to_world);
            self.update_positions_kernel.set_argument(4, &self.settings.world_to_local);
        }
        if self.settings.num_particles != self.prev_settings.num_particles
        {
            changed = true;
            self.sort_spatial_indices_kernel.set_argument(0, self.settings.num_particles);
            self.calculate_offsets_kernel.set_argument(0, self.settings.num_particles);
            self.external_forces_kernel.set_argument(1, self.settings.num_particles);
            self.update_spatial_hash_kernel.set_argument(0, self.settings.num_particles);
            self.calculate_densities_kernel.set_argument(0, self.settings.num_particles);
            self.calculate_pressure_force_kernel.set_argument(1, self.settings.num_particles);
            self.calculate_viscosity_kernel.set_argument(1, self.settings.num_particles);
            self.update_positions_kernel.set_argument(1, self.settings.num_particles);
            self.settings.num_sort_stages = FluidSystem::calc_num_sort_stages(self.settings.num_particles);
        }
        if self.settings.gravity != self.prev_settings.gravity
        {
            changed = true;
            self.external_forces_kernel.set_argument(2, self.settings.gravity);
        }
        if self.settings.collision_damping != self.prev_settings.collision_damping
        {
            changed = true;
            self.update_positions_kernel.set_argument(2, self.settings.collision_damping);
        }
        if self.settings.smoothing_radius != self.prev_settings.smoothing_radius
        {
            changed = true;
            self.update_spatial_hash_kernel.set_argument(1, self.settings.smoothing_radius);
            self.calculate_densities_kernel.set_argument(1, self.settings.smoothing_radius);
            self.calculate_pressure_force_kernel.set_argument(2, self.settings.smoothing_radius);
            self.calculate_viscosity_kernel.set_argument(2, self.settings.smoothing_radius);
        }
        if self.settings.target_density != self.prev_settings.target_density
        {
            changed = true;
            self.calculate_pressure_force_kernel.set_argument(3, self.settings.target_density);
        }
        if self.settings.pressure_multiplier != self.prev_settings.pressure_multiplier
        {
            changed = true;
            self.calculate_pressure_force_kernel.set_argument(4, self.settings.pressure_multiplier);
        }
        if self.settings.near_pressure_multiplier != self.prev_settings.near_pressure_multiplier
        {
            changed = true;
            self.calculate_pressure_force_kernel.set_argument(5, self.settings.near_pressure_multiplier);
        }
        if self.settings.viscosity_strength != self.prev_settings.viscosity_strength
        {
            changed = true;
            self.calculate_viscosity_kernel.set_argument(3, self.settings.viscosity_strength);
        }

        if changed
        {
            self.prev_settings = self.settings;
        }
    }

    pub fn simulate(&mut self, cl: &OpenCL, delta_time: f32)
    {
        if self.settings.paused
        {
            return;
        }

        self.update_settings();
        let delta_time = delta_time / self.settings.time_scale;
        self.external_forces_kernel.set_argument(0, delta_time);
        self.calculate_pressure_force_kernel.set_argument(0, delta_time);
        self.calculate_viscosity_kernel.set_argument(0, delta_time);
        self.update_positions_kernel.set_argument(0, delta_time);

        // pre calc some forces
        self.external_forces_kernel.run(cl, MAX_PARTICLES);

        // update the spatial hashes
        self.update_spatial_hash_kernel.run(cl, MAX_PARTICLES);

        // sort
        let num_stages = self.settings.num_sort_stages;
        for stage_index in 0..num_stages
        {
            for step_index in 0..(stage_index + 1)
            {
                let group_width = 1 << (stage_index - step_index);
                let group_height = 2 * group_width - 1;
                self.sort_spatial_indices_kernel.set_argument(1, group_width);
                self.sort_spatial_indices_kernel.set_argument(2, group_height);
                self.sort_spatial_indices_kernel.set_argument(3, step_index);
                self.sort_spatial_indices_kernel.run(cl, MAX_PARTICLES);
            }
        }

        // calc offsets
        self.calculate_offsets_kernel.run(cl, MAX_PARTICLES);

        // calc densities
        self.calculate_densities_kernel.run(cl, MAX_PARTICLES);

        self.particle_densities.copy_from_device(cl);
        // calc pressures
        self.calculate_pressure_force_kernel.run(cl, MAX_PARTICLES);

        self.particle_velocities.copy_from_device(cl);

        // calc viscosity
        self.calculate_viscosity_kernel.run(cl, MAX_PARTICLES);

        self.particle_velocities.copy_from_device(cl);

        // update positions
        self.update_positions_kernel.run(cl, MAX_PARTICLES);
    }


    // this kinda prevents multiple fluids unless I start pre-allocating memory
    pub fn add_water_particle(&mut self, seed: &mut u32)
    {

    }
}