#include "src/kernels/objects/spheres.cl"

__kernel void generate_rays(
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

    float u = (float)x / (float)screen_width;
    float v = (float)y / (float)screen_height;

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