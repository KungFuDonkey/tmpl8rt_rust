#include "src/kernels/objects/spheres.cl"
#include "src/kernels/tools/random.cl"

__kernel void generate_rays(
    uint glob_seed,
    uint screen_width,
    uint screen_height,
    float3 cam_position,
    float3 cam_top_left,
    float3 cam_bottom_left,
    float3 cam_top_right,
    __global float* ts,
    __global float3* origins,
    __global float3* directions,
    __global float3* normals,
    __global ulong* materials,
    __global float3* ray_colors,
    __global float3* ray_lights)
{
    uint x = get_global_id(0);
    uint y = get_global_id(1);

    uint seed = init_seed(glob_seed + x + y * screen_width);

    float pixel_width = 1.0f / screen_width;
    float pixel_height = 1.0f / pixel_height;


    //float u = ((float)x + pixel_width * random_float(&seed)) / (float)screen_width;
    //float v = ((float)y + pixel_height * random_float(&seed)) / (float)screen_height;

    float u = ((float)x + random_float(&seed)) / (float)screen_width;
    float v = ((float)y + random_float(&seed)) / (float)screen_height;

    float3 p = cam_top_left + (cam_top_right - cam_top_left) * u + (cam_bottom_left - cam_top_left) * v;
    float3 d = normalize(p - cam_position);

    uint idx = x + y * screen_width;

    ts[idx] = 1e30f;
    origins[idx] = cam_position;
    directions[idx] = d;
    normals[idx] = d; // set to something
    materials[idx] = MAX_ULONG;
    ray_colors[idx] = (float3)1;
    ray_lights[idx] = (float3)0;
}