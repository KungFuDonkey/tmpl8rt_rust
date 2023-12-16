#include "src/kernels/objects/spheres.cl"
#include "src/kernels/objects/planes.cl"
#include "src/kernels/objects/quads.cl"

__kernel void extend(
    uint num_rays,
    __global float* ts,
    __global float3* origins,
    __global float3* directions,
    __global float3* intersections,
    __global float3* normals,
    __global uint* materials,
    uint num_spheres,
    __global float3* sphere_positions,
    __global float* sphere_radi,
    __global uint* sphere_materials,
    uint num_planes,
    __global float3* plane_normals,
    __global float* plane_distances,
    __global uint* plane_materials,
    uint num_quads,
    __global float* quad_sizes,
    __global struct mat4* quad_inv_transforms,
    __global struct uint* quad_materials)
{
    int idx = get_global_id(0);
    if (idx >= num_rays)
    {
        return;
    }

    float ray_t = ts[idx];
    float original_ray_t = ray_t;
    float3 ray_origin = origins[idx];
    float3 ray_direction = directions[idx];
    float3 ray_intersection = (float3)0;
    float3 ray_normal = (float3)0;
    uint ray_material = 0;

    intersect_spheres(&ray_t, &ray_origin, &ray_direction, &ray_intersection, &ray_normal, &ray_material, num_spheres, sphere_positions, sphere_radi, sphere_materials);
    intersect_planes(&ray_t, &ray_origin, &ray_direction, &ray_intersection, &ray_normal, &ray_material, num_planes, plane_normals, plane_distances, plane_materials);
    intersect_quads(&ray_t, &ray_origin, &ray_direction, &ray_intersection, &ray_normal, &ray_material, num_quads, quad_sizes, quad_inv_transforms, quad_materials);

    if (ray_t >= original_ray_t)
    {
        return;
    }

    ts[idx] = ray_t;
    origins[idx] = ray_origin;
    directions[idx] = ray_direction;
    intersections[idx] = ray_intersection;
    normals[idx] = ray_normal;
    materials[idx] = ray_material;
}