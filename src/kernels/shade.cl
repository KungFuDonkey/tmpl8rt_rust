#include "src/kernels/tools/constants.cl"
#include "src/kernels/tools/color_conversion.cl"
#include "src/kernels/objects/spheres.cl"

__kernel void shade(
    uint num_rays,
    __global float* ts,
    __global float3* origins,
    __global float3* directions,
    __global float3* intersections,
    __global float3* normals,
    __global uint* materials,
    __global float3* sphere_positions,
    __global float* sphere_radi,
    __global float3* accumulator)
{
    int idx = get_global_id(0);
    if (idx >= num_rays)
    {
        return;
    }

    float t = ts[idx];
    if (t == 1e30)
    {
        accumulator[idx] = (float3)0;
        return;
    }

    uint material = materials[idx];
    float3 color = from_uint_to_float3(material);
    accumulator[idx] = color;
}