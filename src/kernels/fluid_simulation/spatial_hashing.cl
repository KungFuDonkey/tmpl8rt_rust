#pragma once

// ----------------------------------------- SPATIAL HASHING -----------------------------------------------
const int3 sp_offsets[27] =
{
    (int3)(-1, -1, -1),
	(int3)(-1, -1, 0),
	(int3)(-1, -1, 1),
	(int3)(-1, 0, -1),
	(int3)(-1, 0, 0),
	(int3)(-1, 0, 1),
	(int3)(-1, 1, -1),
	(int3)(-1, 1, 0),
	(int3)(-1, 1, 1),
	(int3)(0, -1, -1),
	(int3)(0, -1, 0),
	(int3)(0, -1, 1),
	(int3)(0, 0, -1),
	(int3)(0, 0, 0),
	(int3)(0, 0, 1),
	(int3)(0, 1, -1),
	(int3)(0, 1, 0),
	(int3)(0, 1, 1),
	(int3)(1, -1, -1),
	(int3)(1, -1, 0),
	(int3)(1, -1, 1),
	(int3)(1, 0, -1),
	(int3)(1, 0, 0),
	(int3)(1, 0, 1),
	(int3)(1, 1, -1),
	(int3)(1, 1, 0),
	(int3)(1, 1, 1)
};

#define INTERNAL_NEIGHBOURHOOD(offset, ops)                    \
    {                                                          \
        uint hash = hash_cell(origin_cell + offset);           \
        uint key = key_from_hash(hash, num_particles);         \
        uint current_index = spatial_offsets[key];             \
        while (current_index < num_particles)                  \
        {                                                      \
            uint3 index_data = spatial_indices[current_index]; \
            current_index++;                                   \
            if (index_data.z != key) break;                    \
            if (index_data.y != hash) continue;                \
            uint neighbour_index = index_data.x;               \
            ops                                                \
        }                                                      \
    }

#define UNROLLED_NEIGHBOURHOOD(ops) \
    {\
        INTERNAL_NEIGHBOURHOOD(((int3)(-1, -1, -1)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(-1, -1, 0)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(-1, -1, 1)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(-1, 0, -1)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(-1, 0, 0)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(-1, 0, 1)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(-1, 1, -1)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(-1, 1, 0)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(-1, 1, 1)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(0, -1, -1)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(0, -1, 0)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(0, -1, 1)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(0, 0, -1)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(0, 0, 0)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(0, 0, 1)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(0, 1, -1)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(0, 1, 0)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(0, 1, 1)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(1, -1, -1)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(1, -1, 0)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(1, -1, 1)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(1, 0, -1)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(1, 0, 0)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(1, 0, 1)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(1, 1, -1)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(1, 1, 0)), ops) \
        INTERNAL_NEIGHBOURHOOD(((int3)(1, 1, 1)), ops) \
    }\

const uint hash1 = 15823;
const uint hash2 = 9737333;
const uint hash3 = 440817757;

// Convert floating point position into an integer cell coordinate
int3 get_cell(float3 position, float radius)
{
    float3 floored = floor(position / radius);
	return convert_int3(floored);
}

// Hash cell coordinate to a single unsigned integer
uint hash_cell(int3 cell)
{
    uint3 cell2 = convert_uint3(cell);
	return (cell2.x * hash1) + (cell2.y * hash2) + (cell2.z * hash3);
}

uint key_from_hash(uint hash, uint tableSize)
{
	return hash % tableSize;
}