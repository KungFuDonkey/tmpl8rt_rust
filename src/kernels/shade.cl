#include "src/kernels/tools/constants.cl"
#include "src/kernels/objects/spheres.cl"

float3 get_normal(
    float3* intersection,
    uint obj_id,
    __global float3* sphere_positions,
    __global float* sphere_inv_radi)
{
    uint obj_type = obj_id & OBJ_TYPE_MASK;
    if (obj_type == SPHERE_OBJ_TYPE)
    {
        return get_sphere_normal(intersection, obj_id, sphere_positions, sphere_inv_radi);
    }
    uint real_obj_id = obj_id & OBJ_ID_MASK;

    return (float3)0;
}

__kernel void shade(
    uint num_rays,
    __global float* ts,
    __global float3* origins,
    __global float3* directions,
    __global uint* obj_ids,
    __global float3* sphere_positions,
    __global float* sphere_radi2,
    __global float* sphere_inv_radi,
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

    return;

    float ray_t = ts[idx];
    float3 ray_origin = origins[idx];
    float3 ray_direction = directions[idx];

    float3 intersection = ray_origin + ray_t * ray_direction;
}