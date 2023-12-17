#pragma once
#include "src/kernels/tools/constants.cl"

uint wang_hash(uint seed)
{
    uint v = (seed ^ 61) ^ (seed >> 16);
    v *= 9;
    v = v ^ (v >> 4);
    v *= 0x27d4eb2d;
    v = v ^ (v >> 15);
    return v;
}

uint init_seed(uint seed_base)
{
    return wang_hash((seed_base + 1) * 17);
}

uint random_uint(uint* seed)
{
    *seed ^= *seed << 13;
    *seed ^= *seed >> 17;
    *seed ^= *seed << 5;
    return *seed;
}

float random_float(uint* seed)
{
    return (float)(random_uint(seed)) * 2.3283064365387e-10;
}

float random_float_normal_distribution(uint* seed)
{
    float theta = 2.0f * PI * random_float(seed);
    float rho = sqrt(-2.0f * log(random_float(seed)));
    return rho * cos(theta);
}

float3 random_direction(uint* seed)
{
    float x = random_float_normal_distribution(seed);
    float y = random_float_normal_distribution(seed);
    float z = random_float_normal_distribution(seed);
    return normalize((float3)(x,y,z));
}

float3 random_hemisphere_direction(float3* normal, uint* seed)
{
    float3 dir = random_direction(seed);
    return dir * sign(dot(*normal, dir));
}