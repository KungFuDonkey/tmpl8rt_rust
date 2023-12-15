

void intersect_sphere(float& ray_t, float3& ray_origin, float3& ray_direction, float3& sphere_position, float& sphere_radius2)
{

}



__kernel void intersect_spheres(__global int num_rays, __global float* ts, __global float3* origins, __global float3* directions, __global uint* obj_ids, __global uint* obj_types, __global int num_spheres, __global float3* sphere_positions, __global float* sphere_radi2, __global int* sphere_obj_ids)
{
    int idx = get_global_id(0);
    if idx >= num_rays
    {
        return;
    }

    float ray_t = ts[idx];
    float3 ray_origin = origins[idx];
    float3 ray_direction = directions[idx];
    uint obj_id = obj_ids[idx];
    uint obj_type = obj_types[idx];
}