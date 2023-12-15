#include "src/kernels/constants.cl"

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
    uint* ray_obj_id,
    uint num_spheres,
    float3* sphere_positions,
    float* sphere_radi2)
{
    for (uint i = 0; i < num_spheres; i++)
    {
        float3 sphere_position = sphere_positions[i];
        float sphere_radius2 = sphere_radi2[i];

        if (!intersect_sphere(ray_t, ray_origin, ray_direction, &sphere_position, &sphere_radius2))
        {
            continue;
        }

        *ray_obj_id = i;
    }
}

float3 get_sphere_normal(
    float3* intersection,
    uint obj_id,
    float3* sphere_positions,
    float* sphere_inv_radi)
{
    float3 sphere_position = sphere_positions[obj_id];
    float sphere_inv_radius = sphere_inv_radi[obj_id];
    return (*intersection - sphere_position) * sphere_inv_radius;
}
