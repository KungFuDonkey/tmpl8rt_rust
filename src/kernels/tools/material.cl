#include "src/kernels/tools/color_conversion.cl"
#include "src/kernels/tools/constants.cl"


void get_material_properties(ulong material, float3* color, float* emission_strength)
{
    ulong mat_type = material >> 61;
    uint mat_color_description = (uint)(material & MAX_UINT);
    if (mat_type == MAT_TYPE_LINEAR_COLOR || mat_type == MAT_TYPE_LIGHT)
    {
        *color = from_uint_to_float3(mat_color_description);
    }
    if (mat_type == MAT_TYPE_LIGHT)
    {
        *emission_strength = 4.0f;
    }
}