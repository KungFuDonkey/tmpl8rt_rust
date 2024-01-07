#include "src/kernels/tools/constants.cl"
#include "src/kernels/types/mat4.cl"
#include "src/kernels/fluid_simulation/spatial_hashing.cl"

// ----------------------------------------- SORTING -------------------------------------------------------

__kernel void sort_spatial_indices
(
    uint num_particles,
    uint group_width,
    uint group_height,
    uint step_index,
    __global uint3* spatial_indices
)
{
    uint idx = get_global_id(0);
    uint h_index = idx & (group_width - 1);
    uint index_left = h_index + (group_height + 1) * (idx / group_width);
    uint right_step_size = step_index == 0 ? group_height - 2 * h_index : (group_height + 1) / 2;
    uint index_right = index_left + right_step_size;

    if (index_right >= num_particles)
    {
        return;
    }

    uint3 left_value = spatial_indices[index_left];
    uint3 right_value = spatial_indices[index_right];

    if (left_value.z <= right_value.z)
    {
        return;
    }

    spatial_indices[index_left] = right_value;
    spatial_indices[index_right] = left_value;
}

__kernel void calculate_offsets(
    uint num_particles,
    __global uint3* spatial_indices,
    __global uint* spatial_offsets
)
{
    uint idx = get_global_id(0);
    if (idx >= num_particles)
    {
        return;
    }

    uint3 value = spatial_indices[idx];
    uint key_prev = idx == 0 ? num_particles : spatial_indices[idx - 1].z;

    if (value.z == key_prev) return;

    spatial_offsets[value.z] = idx;
}


// -------------------------------------------- SIMULATION -------------------------------------------------

__kernel void external_forces(
    float delta_time,
    uint num_particles,
    float gravity,
    __global float3* particle_positions,
    __global float3* predicted_particle_positions,
    __global float3* particle_velocities
)
{
    uint idx = get_global_id(0);
    if (idx >= num_particles)
    {
        return;
    }

    float3 velocity = particle_velocities[idx];
    float3 position = particle_positions[idx];
    velocity += (float3)(0, gravity, 0) * delta_time;

    // mul by constant framerate for consistent simulations
    float3 predicted_position = position + velocity * (1 / 120.0f);

    predicted_particle_positions[idx] = predicted_position;
    particle_velocities[idx] = velocity;
}


__kernel void update_spatial_hash
(
    uint num_particles,
    float smoothing_radius,
    __global float3* predicted_particle_positions,
    __global uint3* spatial_indices,
    __global uint* spatial_offsets
)
{
    uint idx = get_global_id(0);
    if (idx >= num_particles)
    {
        return;
    }

    float3 predicted_position = predicted_particle_positions[idx];

    int3 cell = get_cell(predicted_position, smoothing_radius);
    uint hash = hash_cell(cell);
    uint key = key_from_hash(hash, num_particles);

    spatial_indices[idx] = (uint3)(idx, hash, key);

    // reset offsets
    spatial_offsets[idx] = num_particles;
}

__kernel void calculate_densities(
    uint num_particles,
    float smoothing_radius,
    __global float3* predicted_particle_positions,
    __global float2* particle_densities,
    __global uint3* spatial_indices,
    __global uint* spatial_offsets
)
{
    uint idx = get_global_id(0);
    if (idx >= num_particles)
    {
        return;
    }

    float3 position = predicted_particle_positions[idx];
    int3 origin_cell = get_cell(position, smoothing_radius);
    float radius_2 = smoothing_radius * smoothing_radius;
    float density = 0;
    float near_density = 0;

    for (int i = 0; i < 27; i++)
    {
        uint hash = hash_cell(origin_cell + sp_offsets[i]);
        uint key = key_from_hash(hash, num_particles);
        uint current_index = spatial_offsets[key];

        while (current_index < num_particles)
        {
            uint3 index_data = spatial_indices[current_index];
            current_index++;

            // Exit if no longer looking at correct bin
            if (index_data.z != key) break;
            // Skip if hash does not match
            if (index_data.y != hash) continue;

            uint neighbour_index = index_data.x;

        }
    }

    UNROLLED_NEIGHBOURHOOD(
        float3 neighbour_position = predicted_particle_positions[neighbour_index];
        float3 diff = neighbour_position - position;
        float distance_2 = dot(diff, diff);

        // Skip if not within radius
        if (distance_2 > radius_2) continue;

        // Calculate density and near density
        float distance = sqrt(distance_2);
        float scale = 15 / (2 * PI * radius_2 * radius_2 * smoothing_radius); // radius ^ 5
        float near_scale = 15 / (PI * radius_2 * radius_2 * radius_2); // radius ^ 6
        float v = smoothing_radius - distance;
        float v_2 = v * v;
        density += v_2 * scale;
        near_density += v_2 * v * near_scale;);

    particle_densities[idx] = (float2)(density, near_density);
}

__kernel void calculate_pressure_force(
    float delta_time,
    uint num_particles,
    float smoothing_radius,
    float target_density,
    float pressure_multiplier,
    float near_pressure_multiplier,
    __global float3* predicted_particle_positions,
    __global float3* particle_velocities,
    __global float2* particle_densities,
    __global uint3* spatial_indices,
    __global uint* spatial_offsets
)
{
    uint idx = get_global_id(0);
    if (idx >= num_particles)
    {
        return;
    }

    float2 densities = particle_densities[idx];
    float3 position = predicted_particle_positions[idx];
    float3 velocity = particle_velocities[idx];

    float density = densities.x;
    float near_density = densities.y;

    float pressure = (density - target_density) * pressure_multiplier;
    float near_pressure = near_density * near_pressure_multiplier;

    float3 pressure_force = (float3)0.0f;

    int3 origin_cell = get_cell(position, smoothing_radius);
    float radius_2 = smoothing_radius * smoothing_radius;


    UNROLLED_NEIGHBOURHOOD(
        // skip self
        if (neighbour_index == idx) continue;

        float3 neighbour_position = predicted_particle_positions[neighbour_index];
        float3 diff = neighbour_position - position;
        float distance_2 = dot(diff, diff);

        // Skip if not within radius
        if (distance_2 > radius_2) continue;

        // calc pressure_force
        float2 neighbour_densities = particle_densities[neighbour_index];
        float neighbour_density = neighbour_densities.x;
        float neighbour_near_density = neighbour_densities.y;
        float neighbour_pressure = (neighbour_density - target_density) * pressure_multiplier;
        float neighbour_near_pressure = neighbour_near_density * near_pressure_multiplier;

        float shared_pressure = (pressure + neighbour_pressure) / 2;
        float shared_near_pressure = (near_pressure + neighbour_near_pressure) / 2;

        float distance = sqrt(distance_2);

        // go up if points on same position
        float3 direction = distance > 0 ? diff / distance : (float3)(0.0f, 1.0f, 0.0f);

        float radius_4 = radius_2 * radius_2;
        float pressure_scale = 15.0f / (radius_4 * smoothing_radius * PI);
        float near_pressure_scale = 45.0f / (radius_4 * radius_2 * PI);
        float v = smoothing_radius - distance;

        pressure_force += direction * (-v * pressure_scale) * shared_pressure / neighbour_density;
        pressure_force += direction * (-v * v * near_pressure_scale) * shared_near_pressure / neighbour_near_density;);

    float3 acceleration = pressure_force / density;
    particle_velocities[idx] = velocity + acceleration * delta_time;
}

__kernel void calculate_viscosity(
    float delta_time,
    uint num_particles,
    float smoothing_radius,
    float viscosity_strength,
    __global float3* predicted_particle_positions,
    __global float3* particle_velocities,
    __global uint3* spatial_indices,
    __global uint* spatial_offsets
)
{
    uint idx = get_global_id(0);
    if (idx >= num_particles)
    {
        return;
    }

    float3 position = predicted_particle_positions[idx];
    float3 velocity = particle_velocities[idx];

    int3 origin_cell = get_cell(position, smoothing_radius);
    float radius_2 = smoothing_radius * smoothing_radius;

    float3 viscosity_force = 0.0f;

    UNROLLED_NEIGHBOURHOOD(
        // skip self
        if (neighbour_index == idx) continue;

        float3 neighbour_position = predicted_particle_positions[neighbour_index];
        float3 diff = neighbour_position - position;
        float distance_2 = dot(diff, diff);

        // Skip if not within radius
        if (distance_2 > radius_2) continue;

        float3 neighbour_velocity = particle_velocities[neighbour_index];

        float radius_4 = radius_2 * radius_2;

        float scale = 315.0f / (64 * PI * radius_4 * radius_4 * smoothing_radius);
        float v = radius_2 - distance_2;
        viscosity_force += (neighbour_velocity - velocity) * v * v * v * scale;);

    particle_velocities[idx] = velocity + viscosity_force * viscosity_strength * delta_time;
}

__kernel void update_positions
(
    float delta_time,
    uint num_particles,
    float collision_damping,
    struct mat4 local_to_world,
    struct mat4 world_to_local,
    __global float3* particle_positions,
    __global float3* particle_velocities
)
{
    uint idx = get_global_id(0);
    if (idx >= num_particles)
    {
        return;
    }

    float3 position = particle_positions[idx];
    float3 velocity = particle_velocities[idx];

    position += velocity * delta_time;

    float3 position_local = transform_position(&position, &world_to_local);
    float3 velocity_local = transform_vector(&velocity, &world_to_local);

    float3 half_size = (float3)0.5f;
    float3 edge_distance = half_size - fabs(position_local);

    if (edge_distance.x <= 0.0f)
    {
        position_local.x = half_size.x * sign(position_local.x);
        velocity_local.x *= -1.0f * collision_damping;
    }
    if (edge_distance.y <= 0.0f)
    {
        position_local.y = half_size.y * sign(position_local.y);
        velocity_local.y *= -1.0f * collision_damping;
    }
    if (edge_distance.z <= 0.0f)
    {
        position_local.z = half_size.z * sign(position_local.z);
        velocity_local.z *= -1.0f * collision_damping;
    }


    position = transform_position(&position_local, &local_to_world);
    velocity = transform_vector(&velocity_local, &local_to_world);

    particle_positions[idx] = position;
    particle_velocities[idx] = velocity;
}
