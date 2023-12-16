#include "src/kernels/types/mat4.cl"


bool intersect_quad(
    float* ray_t,
    float3* ray_origin,
    float3* ray_direction,
    struct mat4* quad_inv_transform,
    float* quad_size)
{
    float o_y = quad_inv_transform->cell[4] * ray_origin->x + quad_inv_transform->cell[5] * ray_origin->y + quad_inv_transform->cell[6] * ray_origin->z + quad_inv_transform->cell[7];
    float d_y = quad_inv_transform->cell[4] * ray_direction->x + quad_inv_transform->cell[5] * ray_direction->y + quad_inv_transform->cell[6] * ray_direction->z;
    float t = o_y / -d_y;
    if (t < *ray_t && t > 0.0)
    {
        float o_x = quad_inv_transform->cell[0] * ray_origin->x + quad_inv_transform->cell[1] * ray_origin->y + quad_inv_transform->cell[2] * ray_origin->z + quad_inv_transform->cell[3];
        float o_z = quad_inv_transform->cell[8] * ray_origin->x + quad_inv_transform->cell[9] * ray_origin->y + quad_inv_transform->cell[10] * ray_origin->z + quad_inv_transform->cell[11];
        float d_x = quad_inv_transform->cell[0] * ray_direction->x + quad_inv_transform->cell[1] * ray_direction->y + quad_inv_transform->cell[2] * ray_direction->z;
        float d_z = quad_inv_transform->cell[8] * ray_direction->x + quad_inv_transform->cell[9] * ray_direction->y + quad_inv_transform->cell[10] * ray_direction->z;
        float i_x = o_x + t * d_x;
        float i_z = o_z + t * d_z;
        if (i_x > -*quad_size && i_x < *quad_size && i_z > -*quad_size && i_z < *quad_size)
        {
            *ray_t = t;
            return true;
        }
    }
    return false;
}


void intersect_quads(
    float* ray_t,
    float3* ray_origin,
    float3* ray_direction,
    float3* ray_intersection,
    float3* ray_normal,
    uint* ray_material,
    uint num_quads,
    float* quad_sizes,
    struct mat4* quad_inv_transforms,
    uint* quad_materials)
{
    for (uint i = 0; i < num_quads; i++)
    {
        float quad_size = quad_sizes[i];
        struct mat4 quad_inv_transform = quad_inv_transforms[i];

        if (!intersect_quad(ray_t, ray_origin, ray_direction, &quad_inv_transform, &quad_size))
        {
            continue;
        }

        struct mat4 quad_transform = invert_mat4(&quad_inv_transform);
        *ray_material = quad_materials[i];
        *ray_intersection = *ray_t * *ray_direction + *ray_origin;
        *ray_normal = (float3)(quad_transform.cell[1], quad_transform.cell[5], quad_transform.cell[9]);
    }
}

