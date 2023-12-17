

__kernel void finalize(
    uint rendered_frames,
    __global float3* colors,
    __global float3* accumulator,
    __global uint* output_buffer)
{
    uint idx = get_global_id(0);

    float3 accumulated_color = accumulator[idx];
    float3 rendered_color = colors[idx];

    float weight = 1.0f / (float)(rendered_frames);
    float3 new_color = accumulated_color * (1 - weight) + rendered_color * weight;

    float3 one = (float3)1;
    float3 ranged_color = min(new_color, one) * 255.0f;

    uint r = (uint)ranged_color.x;
    uint g = (uint)ranged_color.y;
    uint b = (uint)ranged_color.z;

    output_buffer[idx] = (r << 16) + (g << 8) + b;
    accumulator[idx] = new_color;
}