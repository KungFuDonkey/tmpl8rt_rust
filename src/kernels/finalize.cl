

__kernel void finalize(__global float3* colors, __global uint* output_buffer)
{
    uint idx = get_global_id(0);

    float3 one = {1.0f, 1.0f, 1.0f};
    float3 ranged_color = min(colors[idx], one) * 255.0f;

    uint r = (uint)ranged_color.x;
    uint g = (uint)ranged_color.y;
    uint b = (uint)ranged_color.z;

    output_buffer[idx] = (r << 16) + (g << 8) + b;
}