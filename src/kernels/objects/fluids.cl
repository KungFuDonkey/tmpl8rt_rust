#include "src/kernels/objects/aabb.cl"
#include "src/kernels/objects/spheres.cl"

#define PARTICLE_RADIUS 0.006f
#define PARTICLE_RADIUS 0.006f
#define PARTICLE_RADIUS_2 (PARTICLE_RADIUS * PARTICLE_RADIUS)
#define PARTICLE_INV_RADIUS (1.0f / PARTICLE_RADIUS)

void intersect_fluid(
    float* ray_t,
    float3* ray_origin,
    float3* ray_direction,
    float3* ray_normal,
    ulong* ray_material,
    uint fluid_num_particles,
    struct mat4* local_to_world,
    struct mat4* world_to_local,
    float3* particle_positions,
    float3* particle_velocities,
    float3* particle_densities
)
{
    //{
    //    float t = intersect_aabb(ray_t, ray_origin, ray_direction, fluid_min_bound, fluid_max_bound);
    //    if (t == 1e30)
    //    {
    //        //return false;
    //    }
    //}

    float particle_radius_2 = PARTICLE_RADIUS_2;

    uint particle_idx = 0;
    // todo do grid based intersection
    bool intersected = false;
    for (uint i = 0; i < fluid_num_particles; i++)
    {
        float3 particle_position = particle_positions[i];
        if (intersect_sphere(ray_t, ray_origin, ray_direction, &particle_position, &particle_radius_2))
        {
            intersected = true;
            particle_idx = i;
        }
    }


    if (intersected)
    {
        float3 particle_position = particle_positions[particle_idx];

        float3 velocity = particle_velocities[particle_idx];
        float speed = length(velocity);

        if (speed <= 0)
        {
            *ray_material = (128 << 16) + (128 << 8) + 255;
        }
        else
        {
            *ray_material = (255 << 16);
        }

        float3 intersection = *ray_t * *ray_direction + *ray_origin;
        *ray_normal = (intersection - particle_position) * PARTICLE_INV_RADIUS;
    }
}


/*void intersect_fluids(
    float* ray_t,
    float3* ray_origin,
    float3* ray_direction,
    float3* ray_normal,
    ulong* ray_material,
    uint num_fluids,
    uint* fluid_particle_offsets,
    float3* fluid_min_bounds,
    float3* fluid_max_bounds,
    uint* fluid_particle_counts,
    uint* fluid_particle_ids,
    float3* fluid_particle_positions,
    float3* fluid_particle_colors
    )
{
    uint ray_particle_idx;
    for (int i = 0; i < num_fluids; i++)
    {
        uint particle_offset = fluid_particle_offsets[i];
        float3 fluid_min_bound = fluid_min_bounds[i];
        float3 fluid_max_bound = fluid_max_bounds[i];
        uint fluid_particle_count = fluid_particle_counts[i];

        if (!intersect_fluid(ray_t, ray_origin, ray_direction, &ray_particle_idx, &fluid_min_bound, &fluid_max_bound, fluid_particle_count, fluid_particle_ids + particle_offset, fluid_particle_positions + particle_offset))
        {
            continue;
        }

        float3 particle_position = fluid_particle_positions[particle_offset + ray_particle_idx];
        *ray_material = (200 << 24) + (128 << 16) + (128 << 8) + 255;
        float3 intersection = *ray_t * *ray_direction + *ray_origin;
        *ray_normal = (intersection - particle_position) * PARTICLE_INV_RADIUS;
    }

}*/