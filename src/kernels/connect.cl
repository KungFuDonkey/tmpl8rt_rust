#include "src/kernels/tools/constants.cl"
#include "src/kernels/objects/spheres.cl"
#include "src/kernels/objects/meshes.cl"


void __kernel connect(
    uint num_shadow_rays,
    __global float3* ray_colors,
    __global float3* ray_lights,
    __global float3* ray_normals,
    __global float3* ray_BRDFs,
    __global float* shadow_ray_ts,
    __global float3* shadow_ray_origins,
    __global float3* shadow_ray_directions,
    __global float3* shadow_light_normals,
    __global float* shadow_light_areas,
    __global float3* shadow_light_colors,
    uint num_spheres,
    __global float3* sphere_positions,
    __global float* sphere_radi,
    uint num_quads,
    __global float* quad_sizes,
    __global struct mat4* quad_inv_transforms,
    uint num_meshes,
    __global uint* mesh_offsets,
    __global uint* mesh_triangle_offsets,
    __global struct mat4* mesh_inv_transforms,
    __global float3* mesh_min_bounds,
    __global float3* mesh_max_bounds,
    __global uint* mesh_tri_counts,
    __global uint* mesh_left_firsts,
    __global struct triangle* mesh_triangles
)
{
    int idx = get_global_id(0);
    if (idx >= num_shadow_rays)
    {
        return;
    }

    float ray_t = shadow_ray_ts[idx];

    if (ray_t < 0)
    {
        // todo remove when wavefront
        return;
    }
    float3 ray_origin = shadow_ray_origins[idx];
    float3 ray_direction = shadow_ray_directions[idx];


    if (occlude_spheres(&ray_t, &ray_origin, &ray_direction, num_spheres, sphere_positions, sphere_radi))
    {
        return;
    }
    if (occlude_meshes(&ray_t, &ray_origin, &ray_direction, num_meshes, mesh_offsets, mesh_triangle_offsets, mesh_inv_transforms, mesh_min_bounds, mesh_max_bounds, mesh_tri_counts, mesh_left_firsts, mesh_triangles))
    {
        return;
    }

    float distance2 = ray_t * ray_t;
    float light_area = shadow_light_areas[idx];
    float3 light_normal = shadow_light_normals[idx];
    float3 light_color = shadow_light_colors[idx];
    float3 ray_normal = ray_normals[idx];
    float3 ray_light = ray_lights[idx];
    float3 ray_color = ray_colors[idx];
    float3 ray_BRDF = ray_BRDFs[idx];

    float solidAngle = (dot(light_normal, -ray_direction) * light_area) / distance2;
    float lightPDF = 1 / solidAngle;

    ray_light += ray_color * (dot(ray_normal, ray_direction) / lightPDF) * ray_BRDF * light_color;

    ray_lights[idx] = ray_light;
}