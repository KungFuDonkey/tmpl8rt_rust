#include "src/kernels/tools/constants.cl"
#include "src/kernels/tools/material.cl"
#include "src/kernels/tools/random.cl"
#include "src/kernels/objects/spheres.cl"
#include "src/kernels/objects/quads.cl"

__kernel void shade(
    uint glob_seed,
    uint num_bounces,
    volatile __global uint* num_rays,
    __global uint* write_back_ids,
    __global float* ts,
    __global float3* origins,
    __global float3* directions,
    __global float3* normals,
    __global ulong* materials,
    __global float3* ray_colors,
    __global float3* ray_lights,
    __global uint* shadow_write_back_ids,
    __global float* shadow_ray_ts,
    __global float3* shadow_ray_origins,
    __global float3* shadow_ray_directions,
    __global float3* shadow_light_colors,
    uint num_quads,
    __global float* quad_sizes,
    __global struct mat4* quad_inv_transforms,
    __global ulong* quad_materials
    )
{
    uint idx = get_global_id(0);
    uint max_idx = num_rays[num_bounces * 2];

    if (idx >= max_idx)
    {
        return;
    }

    uint seed = init_seed(glob_seed + idx);
    float t = ts[idx];

    // no hit
    if (t == 1e30 || t == -1.0f)
    {
        return;
    }

    ulong material = materials[idx];
    float3 color = (float3)0;
    float emission_strength = 0;
    get_material_properties(material, &color, &emission_strength);

    float3 ray_color = ray_colors[idx];
    uint write_back_idx = write_back_ids[idx];

    // hit a light
    if (emission_strength > 0.0f)
    {
        ray_lights[write_back_idx] += ray_color * color * emission_strength;
        return;
    }

    float3 BRDF = color * INV_PI;

    float3 normal = normals[idx];
    float3 ray_origin = origins[idx];
    float3 ray_direction = directions[idx];
    float3 intersection = ray_origin + ray_direction * (t - EPSILON);

    // todo create shadow rays for nee
    {
        uint random_quad = random_uint(&seed) % num_quads;
        float quad_size = quad_sizes[random_quad];
        struct mat4 quad_inv_transform = quad_inv_transforms[random_quad];
        struct mat4 quad_transform = invert_mat4(&quad_inv_transform);
        float3 quad_point = random_point_on_quad(&quad_inv_transform, &quad_size, &seed);
        float3 shadow_ray_dir = quad_point - intersection;
        float3 quad_normal = (float3)(-quad_transform.cell[1], -quad_transform.cell[5], -quad_transform.cell[9]);
        if (dot(normal, shadow_ray_dir) > 0.0f && dot(quad_normal, -shadow_ray_dir) > 0.0f)
        {
            float shadow_ray_t = length(shadow_ray_dir);
            shadow_ray_dir = normalize(shadow_ray_dir);
            ulong quad_material = quad_materials[random_quad];
            float3 quad_color = (float3)0;
            float quad_emission_strength = 0;
            get_material_properties(quad_material, &quad_color, &quad_emission_strength);

            uint shadow_ray_idx = atomic_inc(num_rays + num_bounces * 2 + 1);

            shadow_ray_origins[shadow_ray_idx] = intersection;
            shadow_ray_directions[shadow_ray_idx] = shadow_ray_dir;

            float lightPDF = (shadow_ray_t * shadow_ray_t) / (dot(quad_normal, -shadow_ray_dir)); // todo change
            shadow_light_colors[shadow_ray_idx] = ray_color * quad_color * quad_emission_strength * BRDF * (dot(normal, shadow_ray_dir) / lightPDF);

            shadow_ray_ts[shadow_ray_idx] = shadow_ray_t;
            shadow_write_back_ids[shadow_ray_idx] = write_back_idx;
        }
    }

    float survival_chance = 1.0f;

    // russian roulette, start after 2 bounces
    if (num_bounces > 1)
    {
        survival_chance = clamp(max(max(ray_color.x, ray_color.y), ray_color.z), 0.1f, 0.9f);
        if (survival_chance < random_float(&seed))
        {
            return;
        }
    }

    float3 new_direction = normalize(normal + random_direction(&seed));

    //float3 new_direction = normalize(random_hemisphere_direction(&normal, &seed));

    uint new_ray_idx = atomic_inc(num_rays + (num_bounces * 2) + 2);

    write_back_ids[new_ray_idx] = write_back_idx;
    ts[new_ray_idx] = 1e30;
    origins[new_ray_idx] = intersection;
    directions[new_ray_idx] = new_direction;
    ray_colors[new_ray_idx] = ray_color * (1.0f / survival_chance) * PI * BRDF;
}