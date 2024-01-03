#include "src/kernels/objects/spheres.cl"
#include "src/kernels/objects/planes.cl"
#include "src/kernels/objects/quads.cl"
#include "src/kernels/objects/meshes.cl"
#include "src/kernels/objects/fluids.cl"

__kernel void extend(
    uint num_bounces,
    volatile __global uint* num_rays,
    __global float* ts,
    __global float3* origins,
    __global float3* directions,
    __global float3* normals,
    __global ulong* materials,
    uint num_spheres,
    __global float3* sphere_positions,
    __global float* sphere_radi,
    __global ulong* sphere_materials,
    uint num_planes,
    __global float3* plane_normals,
    __global float* plane_distances,
    __global ulong* plane_materials,
    uint num_quads,
    __global float* quad_sizes,
    __global struct mat4* quad_inv_transforms,
    __global ulong* quad_materials,
    uint num_meshes,
    __global uint* mesh_offsets,
    __global uint* mesh_triangle_offsets,
    __global struct mat4* mesh_inv_transforms,
    __global float3* mesh_min_bounds,
    __global float3* mesh_max_bounds,
    __global uint* mesh_tri_counts,
    __global uint* mesh_left_firsts,
    __global struct triangle* mesh_triangles,
    __global float3* mesh_triangle_normals,
    __global ulong* mesh_materials,
    uint num_fluids,
    __global uint* fluid_particle_offsets,
    __global float3* fluid_min_bounds,
    __global float3* fluid_max_bounds,
    __global uint* fluid_particle_counts,
    __global uint* fluid_particle_ids,
    __global float3* fluid_particle_positions,
    __global float3* fluid_particle_colors)
{
    uint idx = get_global_id(0);
    uint max_idx = num_rays[num_bounces * 2];
    if (idx >= max_idx)
    {
        return;
    }

    float ray_t = ts[idx];

    float original_ray_t = ray_t;
    float3 ray_origin = origins[idx];
    float3 ray_direction = directions[idx];
    float3 ray_normal = (float3)0;
    ulong ray_material = 0;

    intersect_spheres(&ray_t, &ray_origin, &ray_direction, &ray_normal, &ray_material, num_spheres, sphere_positions, sphere_radi, sphere_materials);
    intersect_planes(&ray_t, &ray_origin, &ray_direction, &ray_normal, &ray_material, num_planes, plane_normals, plane_distances, plane_materials);
    intersect_quads(&ray_t, &ray_origin, &ray_direction, &ray_normal, &ray_material, num_quads, quad_sizes, quad_inv_transforms, quad_materials);
    intersect_meshes(&ray_t, &ray_origin, &ray_direction, &ray_normal, &ray_material, num_meshes, mesh_offsets, mesh_triangle_offsets, mesh_inv_transforms, mesh_min_bounds, mesh_max_bounds, mesh_tri_counts, mesh_left_firsts, mesh_triangles, mesh_triangle_normals, mesh_materials);
    intersect_fluids(&ray_t, &ray_origin, &ray_direction, &ray_normal, &ray_material, num_fluids, fluid_particle_offsets, fluid_min_bounds, fluid_max_bounds, fluid_particle_counts, fluid_particle_ids, fluid_particle_positions, fluid_particle_colors);

    if (ray_t >= original_ray_t)
    {
        return;
    }

    ts[idx] = ray_t;
    origins[idx] = ray_origin;
    directions[idx] = ray_direction;
    normals[idx] = ray_normal;
    materials[idx] = ray_material;
}