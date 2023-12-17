#include "src/kernels/tools/constants.cl"

bool intersect_sphere(
    float* ray_t,
    float3* ray_origin,
    float3* ray_direction,
    float3* sphere_position,
    float* sphere_radius2)
{
    float3 oc = *ray_origin - *sphere_position;
    float b = dot(oc, *ray_direction);
    float c = dot(oc, oc) - *sphere_radius2;
    float d = b * b - c;
    if (d <= 0.0f)
    {
        return false;
    }

    d = sqrt(d);
    float t = -b - d;
    if (t < *ray_t && t > 0.0f)
    {
        *ray_t = t;
        return true;
    }
    if (c > 0.0f)
    {
        return false;
    }

    t = d - b;
    if (t < *ray_t && t > 0.0f)
    {
        *ray_t = t;
        return true;
    }
    return false;
}

void intersect_spheres(
    float* ray_t,
    float3* ray_origin,
    float3* ray_direction,
    float3* ray_normal,
    ulong* ray_material,
    uint num_spheres,
    float3* sphere_positions,
    float* sphere_radi,
    ulong* sphere_materials)
{
    for (uint i = 0; i < num_spheres; i++)
    {
        float3 sphere_position = sphere_positions[i];
        float sphere_radius = sphere_radi[i];
        float sphere_radius2 = sphere_radius * sphere_radius;

        if (!intersect_sphere(ray_t, ray_origin, ray_direction, &sphere_position, &sphere_radius2))
        {
            continue;
        }

        float sphere_inv_radius = 1.0f / sphere_radius;

        *ray_material = sphere_materials[i];
        float3 intersection = *ray_t * *ray_direction + *ray_origin;
        float3 normal = (intersection - sphere_position) * sphere_inv_radius;
        *ray_normal = dot(normal, *ray_direction) > 0.0f ? -normal : normal;
    }
}

bool occlude_sphere
(
    float* ray_t,
    float3* ray_origin,
    float3* ray_direction,
    float3* sphere_position,
    float* sphere_radius2
)
{
    float3 oc = *ray_origin - *sphere_position;
    float b = dot(oc, *ray_direction);
    float c = dot(oc, oc) - *sphere_radius2;

    float d = b * b - c;
    if (d <= 0.0)
    {
        return false;
    }

    float t = -b - sqrt(d);
    return t < *ray_t && t > 0.0f;
}

bool occlude_spheres(
    float* ray_t,
    float3* ray_origin,
    float3* ray_direction,
    uint num_spheres,
    float3* sphere_positions,
    float* sphere_radi
)
{
    for (uint i = 0; i < num_spheres; i++)
    {
        float3 sphere_position = sphere_positions[i];
        float sphere_radius = sphere_radi[i];
        float sphere_radius2 = sphere_radius * sphere_radius;

        if (occlude_sphere(ray_t, ray_origin, ray_direction, &sphere_position, &sphere_radius2))
        {
            return true;
        }
    }
    return false;
}