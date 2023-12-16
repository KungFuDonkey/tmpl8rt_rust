
#define INV_BYTE_VAL 0.00392156862f

float3 from_uint_to_float3(uint value)
{
    return (float3)(
        (float)(value >> 16) * INV_BYTE_VAL,
        (float)((value >> 8) & 255) * INV_BYTE_VAL,
        (float)(value & 255) * INV_BYTE_VAL);
}