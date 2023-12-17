#include "src/kernels/tools/constants.cl"
#include "src/kernels/tools/material.cl"
#include "src/kernels/tools/random.cl"
#include "src/kernels/objects/spheres.cl"

__kernel void shade(
    uint glob_seed,
    uint num_bounces,
    uint num_rays,
    __global float* ts,
    __global float3* origins,
    __global float3* directions,
    __global float3* normals,
    __global ulong* materials,
    __global float3* ray_colors,
    __global float3* ray_lights)
{
    int idx = get_global_id(0);
    if (idx >= num_rays)
    {
        return;
    }

    uint seed = init_seed(glob_seed + idx);
    float t = ts[idx];

    // no hit
    if (t == 1e30 || t == -1.0f)
    {
        // todo remove when wavefront
        ts[idx] = -1.0f;
        return;
    }

    ulong material = materials[idx];
    float3 color = (float3)0;
    float emission_strength = 0;
    get_material_properties(material, &color, &emission_strength);

    float3 emitted_light = color * emission_strength;
    float3 ray_color = ray_colors[idx];
    float3 ray_light = ray_lights[idx] + emitted_light * ray_color;
    ray_lights[idx] = ray_light;

    // hit a light
    if (emission_strength > 0.0f)
    {
        // todo remove when wavefront
        ts[idx] = -1.0f;
        return;
    }

    // todo create shadow rays for nee

    ray_color *= color;
    // russian roulette, start after 2 bounces
    if (num_bounces > 1)
    {
        float survival_chance = clamp(max(max(ray_color.x, ray_color.y), ray_color.z), 0.1f, 0.9f);
        if (survival_chance < random_float(&seed))
        {
            // todo remove when wavefront
            ts[idx] = -1.0f;
            return;
        }
    }

    float3 normal = normals[idx];
    float3 ray_origin = origins[idx];
    float3 ray_direction = directions[idx];
    float3 intersection = ray_origin + ray_direction * (t - EPSILON);
    float3 new_direction = normalize(normal + random_direction(&seed));

    ts[idx] = 1e30;
    origins[idx] = intersection;
    directions[idx] = new_direction;
    ray_colors[idx] = ray_color;
}