

bool intersect_sphere(float* ray_t, float3* ray_origin, float3* ray_direction, float3* sphere_position, float* sphere_radius2)
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

void intersect_spheres(float* ray_t, float3* ray_origin, float3* ray_direction, uint* ray_obj_id, uint* ray_obj_type, uint num_spheres, float3* sphere_positions, float* sphere_radi2, int* sphere_obj_ids)
{
    for (uint i = 0; i < num_spheres; i++)
    {
        float3 sphere_position = sphere_positions[i];
        float sphere_radius2 = sphere_radi2[i];

        if (!intersect_sphere(ray_t, ray_origin, ray_direction, &sphere_position, &sphere_radius2))
        {
            continue;
        }

        *ray_obj_id = sphere_obj_ids[i];
        *ray_obj_type = 0; // sphere
    }
}

__kernel void intersect_spheres_kernel(
    uint num_rays,
    __global float* ts,
    __global float3* origins,
    __global float3* directions,
    __global uint* obj_ids,
    __global uint* obj_types,
    uint num_spheres,
    __global float3* sphere_positions,
    __global float* sphere_radi2,
    __global uint* sphere_obj_ids)
{
    int idx = get_global_id(0);
    if (idx >= num_rays)
    {
        return;
    }

    float ray_t = ts[idx];
    float3 ray_origin = origins[idx];
    float3 ray_direction = directions[idx];
    uint ray_obj_id = obj_ids[idx];
    uint ray_obj_type = obj_types[idx];

    intersect_spheres(&ray_t, &ray_origin, &ray_direction, &ray_obj_id, &ray_obj_type, num_spheres, sphere_positions, sphere_radi2, sphere_obj_ids);

    ts[idx] = ray_t;
    origins[idx] = ray_origin;
    directions[idx] = ray_direction;
    obj_ids[idx] = ray_obj_id;
    obj_types[idx] = ray_obj_type;
}