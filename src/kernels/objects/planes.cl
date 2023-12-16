
bool intersect_plane(
    float* ray_t,
    float3* ray_origin,
    float3* ray_direction,
    float3* plane_normal,
    float* plane_distance)
{
    float t = -(dot(*ray_origin, *plane_normal) + *plane_distance) / (dot(*ray_direction, *plane_normal));

    if (t < *ray_t && t > 0.0f)
    {
        *ray_t = t;
        return true;
    }
    return false;
}

void intersect_planes(
    float* ray_t,
    float3* ray_origin,
    float3* ray_direction,
    float3* ray_intersection,
    float3* ray_normal,
    uint* ray_material,
    uint num_planes,
    float3* plane_normals,
    float* plane_distances,
    uint* plane_materials)
{
    for (uint i = 0; i < num_planes; i++)
    {
        float3 plane_normal = plane_normals[i];
        float plane_distance = plane_distances[i];

        if (!intersect_plane(ray_t, ray_origin, ray_direction, &plane_normal, &plane_distance))
        {
            continue;
        }

        *ray_material = plane_materials[i];
        *ray_intersection = *ray_t * *ray_direction + *ray_origin;
        *ray_normal = plane_normal;
    }
}
