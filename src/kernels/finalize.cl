

__kernel void finalize(
    uint rendered_frames,
    __global float3* colors,
    __global float3* accumulator,
    __global uint* output_buffer)
{
    uint idx = get_global_id(0);

    float3 accumulated_color = accumulator[idx];
    float3 rendered_color = colors[idx];

    float3 one = (float3)1;
    float3 ranged_color = min(rendered_color, one) * 255.0f;

    float weight = 1.0f / (float)(rendered_frames);
    float3 new_color = accumulated_color * (1 - weight) + ranged_color * weight;

    uint r = (uint)new_color.x;
    uint g = (uint)new_color.y;
    uint b = (uint)new_color.z;

    // todo remove for accumulated rendering
    //r = (uint)ranged_color.x;
    //g = (uint)ranged_color.y;
    //b = (uint)ranged_color.z;

    output_buffer[idx] = (r << 16) + (g << 8) + b;
    accumulator[idx] = new_color;
}