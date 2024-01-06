#include "src/kernels/types/mat4.cl"
#include "src/kernels/tools/random.cl"

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
    float3* ray_normal,
    ulong* ray_material,
    uint num_quads,
    float* quad_sizes,
    struct mat4* quad_inv_transforms,
    ulong* quad_materials)
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
        float3 normal = (float3)(-quad_transform.cell[1], -quad_transform.cell[5], -quad_transform.cell[9]);
        *ray_normal = dot(normal, *ray_direction) > 0.0f ? -normal : normal;
    }
}

float3 quad_centre(
    struct mat4* quad_inv_transform,
    float* quad_size
)
{
    struct mat4 quad_transform = invert_mat4(quad_inv_transform);
    float3 sized_corner1 = (float3)(-*quad_size, 0.0, -*quad_size);
    float3 sized_corner2 = (float3)(*quad_size, 0.0, -*quad_size);
    float3 sized_corner3 = (float3)(*quad_size, 0.0, *quad_size);
    float3 corner1 = transform_position(&sized_corner1, &quad_transform);
    float3 corner2 = transform_position(&sized_corner2, &quad_transform);
    float3 corner3 = transform_position(&sized_corner3, &quad_transform);
    float r1 = 0.5f;
    float r2 = 0.5f;
    return corner1 + r2 * (corner2 - corner1) + r1 * (corner3 - corner1);
}

float3 random_point_on_quad(
    struct mat4* quad_inv_transform,
    float* quad_size,
    uint* seed)
{
    struct mat4 quad_transform = invert_mat4(quad_inv_transform);
    float3 sized_corner1 = (float3)(-*quad_size, 0.0, -*quad_size);
    float3 sized_corner2 = (float3)(*quad_size, 0.0, -*quad_size);
    float3 sized_corner3 = (float3)(*quad_size, 0.0, *quad_size);
    float3 corner1 = transform_position(&sized_corner1, &quad_transform);
    float3 corner2 = transform_position(&sized_corner2, &quad_transform);
    float3 corner3 = transform_position(&sized_corner3, &quad_transform);
    float r1 = random_float(seed);
    float r2 = random_float(seed);
    return corner1 + r2 * (corner2 - corner1) + r1 * (corner3 - corner1);
}

float3 random_point_on_quad_blue_noise(
    struct mat4* quad_inv_transform,
    float* quad_size,
    __constant uchar* blue_noise_texture,
    uint frame_idx,
    uint pixel_x,
    uint pixel_y
)
{
    struct mat4 quad_transform = invert_mat4(quad_inv_transform);
    float3 sized_corner1 = (float3)(-*quad_size, 0.0, -*quad_size);
    float3 sized_corner2 = (float3)(*quad_size, 0.0, -*quad_size);
    float3 sized_corner3 = (float3)(*quad_size, 0.0, *quad_size);
    float3 corner1 = transform_position(&sized_corner1, &quad_transform);
    float3 corner2 = transform_position(&sized_corner2, &quad_transform);
    float3 corner3 = transform_position(&sized_corner3, &quad_transform);

    float2 random_point = random_blue_noise_point(blue_noise_texture, frame_idx, pixel_x, pixel_y);
    return corner1 + random_point.x * (corner2 - corner1) + random_point.y * (corner3 - corner1);
}