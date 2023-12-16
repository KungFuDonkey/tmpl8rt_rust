#include "src/kernels/objects/spheres.cl"

__kernel void extend(
    uint num_rays,
    __global float* ts,
    __global float3* origins,
    __global float3* directions,
    __global float3* intersections,
    __global float3* normals,
    __global uint* obj_ids,
    uint num_spheres,
    __global float3* sphere_positions,
    __global float* sphere_radi)
{
    int idx = get_global_id(0);
    if (idx >= num_rays)
    {
        return;
    }

    float ray_t = ts[idx];
    float3 ray_origin = origins[idx];
    float3 ray_direction = directions[idx];
    float3 ray_intersection = intersections[idx];
    float3 ray_normal = normals[idx];
    uint ray_obj_id = obj_ids[idx];

    intersect_spheres(&ray_t, &ray_origin, &ray_direction, &ray_intersection, &ray_normal, &ray_obj_id, num_spheres, sphere_positions, sphere_radi);

    ts[idx] = ray_t;
    origins[idx] = ray_origin;
    directions[idx] = ray_direction;
    intersections[idx] = ray_intersection;
    normals[idx] = ray_normal;
    obj_ids[idx] = ray_obj_id;
}