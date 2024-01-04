#include "src/kernels/tools/constants.cl"


#define GRAVITY ((float3)(0.0, -9.81, 0.0))

float compute_density(
    uint idx,
    uint particle_count,
    float* fluid_particle_mass,
    float* fluid_particle_densities,
    float3* fluid_particle_positions,
    float* fluid_particle_support_radi
)
{
    float density = 0.0f;
    float3 particle_position = fluid_particle_positions[idx];
    float particle_support_radius = fluid_particle_support_radi[idx];
    float particle_support_radius_2 = particle_support_radius * particle_support_radius;

    // todo grid based operations
    for (uint p = 0; p < particle_count; p++)
    {
        float3 neighbour_position = fluid_particle_positions[p];
        float3 diff = particle_position - neighbour_position;
        float magnitude_2 = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
        if (magnitude_2 <= particle_support_radius_2)
        {
            float pow_mag = particle_support_radius_2 - magnitude_2;
            pow_mag = pow_mag * pow_mag * pow_mag;
            float h_mag = particle_support_radius_2 * particle_support_radius_2;
            h_mag = h_mag * h_mag * particle_support_radius;
            density += ((315.0f * pow_mag) / (64.0f * PI * h_mag)) * fluid_particle_mass[p];
        }
    }

    if (density == 0.0f)
    {
        density = fluid_particle_mass[idx] / .00000001f;
    }

    return density;
}

__kernel void compute_densities(
    uint num_fluids,
    __global uint* fluid_particle_offsets,
    __global uint* fluid_particle_counts,
    __global float* fluid_particle_mass,
    __global float* fluid_particle_densities,
    __global float3* fluid_particle_positions,
    __global float* fluid_particle_support_radi
)
{
    uint idx = get_global_id(0);

    for (uint i = 0; i < num_fluids; i++)
    {
        uint particle_count = fluid_particle_counts[i];
        if (idx >= particle_count)
        {
            continue;
        }

        uint particle_offset = fluid_particle_offsets[i];

        float density = compute_density(idx, particle_count, fluid_particle_mass + particle_offset, fluid_particle_densities + particle_offset, fluid_particle_positions + particle_offset, fluid_particle_support_radi + particle_offset);

        fluid_particle_densities[idx + particle_offset] = density;
    }
}

float compute_pressure(
    uint idx,
    uint particle_count,
    float* fluid_particle_pressures,
    float* fluid_particle_stiffness,
    float* fluid_particle_rest_densities,
    float* fluid_particle_densities
)
{
    float stiffness = fluid_particle_stiffness[idx];
    float density = fluid_particle_densities[idx];
    float rest_density = fluid_particle_rest_densities[idx];
    float mag = (density / rest_density);
    float mag_2 = mag * mag;
    float mag_7 = mag_2 * mag_2 * mag_2 * mag;
    float pressure = stiffness * (mag_7 - 1.0);
    return pressure;
}

__kernel void compute_pressures(
    uint num_fluids,
    __global uint* fluid_particle_offsets,
    __global uint* fluid_particle_counts,
    __global float* fluid_particle_pressures,
    __global float* fluid_particle_stiffness,
    __global float* fluid_particle_rest_densities,
    __global float* fluid_particle_densities
)
{
    uint idx = get_global_id(0);

    for (uint i = 0; i < num_fluids; i++)
    {
        uint particle_count = fluid_particle_counts[i];

        if (idx >= particle_count)
        {
            continue;
        }

        uint particle_offset = fluid_particle_offsets[i];

        float pressure = compute_pressure(idx, particle_count, fluid_particle_pressures + particle_offset, fluid_particle_stiffness + particle_offset, fluid_particle_rest_densities + particle_offset, fluid_particle_densities + particle_offset);

        fluid_particle_pressures[idx + particle_offset] = pressure;
    }
}

float3 compute_acceleration(
    uint idx,
    uint particle_count,
    float* fluid_particle_densities,
    float3* fluid_particle_positions,
    float* fluid_particle_support_radi,
    float* fluid_particle_mass,
    float* fluid_particle_pressures,
    float* fluid_particle_viscosities,
    float3* fluid_particle_velocities
)
{
    float particle_density = fluid_particle_densities[idx];
    float3 particle_position = fluid_particle_positions[idx];
    float particle_support_radius = fluid_particle_support_radi[idx];
    float particle_support_radius_2 = particle_support_radius * particle_support_radius;
    float particle_pressure = fluid_particle_pressures[idx];
    float3 particle_velocity = fluid_particle_velocities[idx];
    float particle_viscosity = fluid_particle_viscosities[idx];

    float3 pressure;
    float3 viscosity;
    float3 diff;
    float p_coeff;
    float3 v_coeff;

    // todo grid based operations
    for (uint i = 0; i < particle_count; i++)
    {
        if (i == idx)
        {
            continue;
        }

        float3 other_position = fluid_particle_positions[i];
        float3 diff = particle_position - other_position;
        float magnitude_2 = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
        if (magnitude_2 <= particle_support_radius_2)
        {
            float other_particle_pressure = fluid_particle_pressures[i];
            float other_particle_mass = fluid_particle_mass[i];
            float other_particle_density = fluid_particle_densities[i];
            float3 other_particle_velocity = fluid_particle_velocities[i];
            float other_particle_volume = other_particle_mass / other_particle_density;
            p_coeff = (particle_pressure + other_particle_pressure) / 2.0f * other_particle_volume;
            float diff_mag = length(diff);
            float pow_mag = (particle_support_radius - diff_mag);
            pow_mag = pow_mag * pow_mag;
            float h_mag = particle_support_radius_2 * particle_support_radius_2 * particle_support_radius_2;
            pressure += diff * ((-45.0f * pow_mag) / (PI * h_mag)) * p_coeff;

            v_coeff = (other_particle_velocity - particle_velocity) * other_particle_volume;
            viscosity += v_coeff * ((45.0f * (particle_support_radius - diff_mag)) / (PI * h_mag));
        }
    }
    viscosity = viscosity * particle_viscosity;

    float3 force = GRAVITY * particle_density + (viscosity - pressure);
    float3 acceleration = force / particle_density;
    return acceleration;
}

__kernel void compute_forces(
    uint num_fluids,
    __global uint* fluid_particle_offsets,
    __global uint* fluid_particle_counts,
    __global float* fluid_particle_densities,
    __global float3* fluid_particle_positions,
    __global float* fluid_particle_support_radi,
    __global float* fluid_particle_mass,
    __global float* fluid_particle_pressures,
    __global float* fluid_particle_viscosities,
    __global float3* fluid_particle_velocities,
    __global float3* fluid_particle_accelerations
)
{
    uint idx = get_global_id(0);

    for (uint i = 0; i < num_fluids; i++)
    {
        uint particle_count = fluid_particle_counts[i];

        if (idx >= particle_count)
        {
            continue;
        }

        uint particle_offset = fluid_particle_offsets[i];

        float3 acceleration = compute_acceleration(idx, particle_count, fluid_particle_densities + particle_offset, fluid_particle_positions + particle_offset, fluid_particle_support_radi + particle_offset, fluid_particle_mass + particle_offset, fluid_particle_pressures + particle_offset, fluid_particle_viscosities + particle_offset, fluid_particle_velocities + particle_offset);

        fluid_particle_accelerations[idx + particle_offset] = acceleration;
    }
}

void compute_position
(
    float dt,
    uint idx,
    uint particle_count,
    float3 min_bounds,
    float3 max_bounds,
    float3* fluid_particle_positions,
    float3* fluid_particle_old_positions,
    float3* fluid_particle_velocities,
    float3* fluid_particle_velocities_half,
    float3* fluid_particle_accelerations
)
{
    float3 velocity_half = fluid_particle_velocities_half[idx];
    float3 acceleration = fluid_particle_accelerations[idx];
    float3 position = fluid_particle_positions[idx];
    float rest_coefficient = 0.9f;

    float3 temp_velocity_half = velocity_half + acceleration * dt;
    float3 temp_position = position + (temp_velocity_half * dt);
    float3 temp_velocity = velocity_half;

    if (temp_position.x > max_bounds.x)
    {
        temp_position.x = max_bounds.x;
        temp_velocity.x = -rest_coefficient * temp_velocity_half.x;
        temp_velocity_half.x = -rest_coefficient * temp_velocity_half.x;
    }
    else if (temp_position.x < min_bounds.x)
    {
        temp_position.x = min_bounds.x;
        temp_velocity.x = -rest_coefficient * temp_velocity_half.x;
        temp_velocity_half.x = -rest_coefficient * temp_velocity_half.x;
    }
    if (temp_position.y > max_bounds.y)
    {
        temp_position.y = max_bounds.y;
        temp_velocity.y = -rest_coefficient * temp_velocity_half.y;
        temp_velocity_half.y = -rest_coefficient * temp_velocity_half.y;
    }
    else if (temp_position.y < min_bounds.x)
    {
        temp_position.y = min_bounds.y;
        temp_velocity.y = -rest_coefficient * temp_velocity_half.y;
        temp_velocity_half.y = -rest_coefficient * temp_velocity_half.y;
    }
    if (temp_position.z > max_bounds.z)
    {
        temp_position.z = max_bounds.z;
        temp_velocity.z = -rest_coefficient * temp_velocity_half.z;
        temp_velocity_half.z = -rest_coefficient * temp_velocity_half.z;
    }
    else if (temp_position.z < min_bounds.z)
    {
        temp_position.z = min_bounds.z;
        temp_velocity.z = -rest_coefficient * temp_velocity_half.z;
        temp_velocity_half.z = -rest_coefficient * temp_velocity_half.z;
    }

    fluid_particle_old_positions[idx] = position;
    fluid_particle_velocities_half[idx] = temp_velocity_half;
    fluid_particle_positions[idx] = temp_position;
    fluid_particle_velocities[idx] = temp_velocity;
}


__kernel void compute_positions
(
    float dt,
    uint num_fluids,
    __global uint* fluid_particle_offsets,
    __global uint* fluid_particle_counts,
    __global float3* fluid_min_bounds,
    __global float3* fluid_max_bounds,
    __global float3* fluid_particle_positions,
    __global float3* fluid_particle_old_positions,
    __global float3* fluid_particle_velocities,
    __global float3* fluid_particle_velocities_half,
    __global float3* fluid_particle_accelerations
)
{
    uint idx = get_global_id(0);
    for (uint i = 0; i < num_fluids; i++)
    {
        uint particle_count = fluid_particle_counts[i];

        if (idx >= particle_count)
        {
            continue;
        }

        uint particle_offset = fluid_particle_offsets[i];
        compute_position(dt, idx, particle_count, fluid_min_bounds[i], fluid_max_bounds[i], fluid_particle_positions + particle_offset, fluid_particle_old_positions + particle_offset, fluid_particle_velocities + particle_offset, fluid_particle_velocities_half + particle_offset, fluid_particle_accelerations + particle_offset);
    }
}


__kernel void simulate(
    uint num_fluids,
    __global uint* fluid_particle_offsets,
    __global float3* fluid_min_bounds,
    __global float3* fluid_max_bounds,
    __global uint* fluid_particle_counts,
    __global uint* fluid_particle_ids,
    __global float* fluid_particle_mass,
    __global float* fluid_particle_pressures,
    __global float* fluid_particle_stiffness,
    __global float* fluid_particle_rest_densities,
    __global float* fluid_particle_densities,
    __global float* fluid_particle_viscosities,
    __global float3* fluid_particle_positions,
    __global float3* fluid_particle_old_positions,
    __global float3* fluid_particle_velocities,
    __global float3* fluid_particle_accelerations,
    __global float* fluid_particle_support_radi,
    __global float* fluid_particle_thresholds,
    __global float* fluid_particle_surface_tensions
)
{

}