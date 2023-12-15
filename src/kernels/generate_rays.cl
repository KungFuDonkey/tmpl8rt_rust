
__kernel void generate_rays(uint screen_width, uint screen_height, float3 cam_position, float3 cam_top_left, float3 cam_bottom_left, float3 cam_top_right, __global float* ts, __global float3* origins, __global float3* directions)
{
    uint idx = get_global_id(0);

    uint y = idx / screen_width;
    uint x = idx - (screen_width * y);

    float u = (float)x / (float)screen_width;
    float v = (float)y / (float)screen_height;

    float3 p = cam_top_left + (cam_top_right - cam_top_left) * u + (cam_bottom_left - cam_top_left) * v;
    float3 d = normalize(p - cam_position);

    ts[idx] = 1e30f;
    origins[idx] = cam_position;
    directions[idx] = d;
}