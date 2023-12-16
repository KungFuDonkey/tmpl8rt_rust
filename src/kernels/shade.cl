#include "src/kernels/tools/constants.cl"
#include "src/kernels/objects/spheres.cl"

__kernel void shade(
    uint num_rays,
    __global float* ts,
    __global float3* origins,
    __global float3* directions,
    __global float3* intersections,
    __global float3* normals,
    __global uint* obj_ids,
    __global float3* sphere_positions,
    __global float* sphere_radi,
    __global float3* accumulator)
{
    int idx = get_global_id(0);
    if (idx >= num_rays)
    {
        return;
    }

    uint ray_obj_id = obj_ids[idx];
    if (ray_obj_id == MAX_UINT)
    {
        accumulator[idx] = (float3)0;
        return;
    }

    if (ray_obj_id == 0)
    {
        accumulator[idx] = (float3)1;
    }
    if (ray_obj_id == 1)
    {
        accumulator[idx] = (float3)(0,1,0);
    }
}