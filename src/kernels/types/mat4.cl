
struct mat4
{
    float cell[16];
};

#define IDENTITY_MATRIX {{1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f}}

struct mat4 invert_mat4(struct mat4* other)
{
    float inv[16] = {
        other->cell[5] * other->cell[10] * other->cell[15] - other->cell[5] * other->cell[11] * other->cell[14] - other->cell[9] * other->cell[6] * other->cell[15] +
            other->cell[9] * other->cell[7] * other->cell[14] + other->cell[13] * other->cell[6] * other->cell[11] - other->cell[13] * other->cell[7] * other->cell[10],
        -other->cell[1] * other->cell[10] * other->cell[15] + other->cell[1] * other->cell[11] * other->cell[14] + other->cell[9] * other->cell[2] * other->cell[15] -
            other->cell[9] * other->cell[3] * other->cell[14] - other->cell[13] * other->cell[2] * other->cell[11] + other->cell[13] * other->cell[3] * other->cell[10],
        other->cell[1] * other->cell[6] * other->cell[15] - other->cell[1] * other->cell[7] * other->cell[14] - other->cell[5] * other->cell[2] * other->cell[15] +
            other->cell[5] * other->cell[3] * other->cell[14] + other->cell[13] * other->cell[2] * other->cell[7] - other->cell[13] * other->cell[3] * other->cell[6],
        -other->cell[1] * other->cell[6] * other->cell[11] + other->cell[1] * other->cell[7] * other->cell[10] + other->cell[5] * other->cell[2] * other->cell[11] -
            other->cell[5] * other->cell[3] * other->cell[10] - other->cell[9] * other->cell[2] * other->cell[7] + other->cell[9] * other->cell[3] * other->cell[6],
        -other->cell[4] * other->cell[10] * other->cell[15] + other->cell[4] * other->cell[11] * other->cell[14] + other->cell[8] * other->cell[6] * other->cell[15] -
            other->cell[8] * other->cell[7] * other->cell[14] - other->cell[12] * other->cell[6] * other->cell[11] + other->cell[12] * other->cell[7] * other->cell[10],
        other->cell[0] * other->cell[10] * other->cell[15] - other->cell[0] * other->cell[11] * other->cell[14] - other->cell[8] * other->cell[2] * other->cell[15] +
            other->cell[8] * other->cell[3] * other->cell[14] + other->cell[12] * other->cell[2] * other->cell[11] - other->cell[12] * other->cell[3] * other->cell[10],
        -other->cell[0] * other->cell[6] * other->cell[15] + other->cell[0] * other->cell[7] * other->cell[14] + other->cell[4] * other->cell[2] * other->cell[15] -
            other->cell[4] * other->cell[3] * other->cell[14] - other->cell[12] * other->cell[2] * other->cell[7] + other->cell[12] * other->cell[3] * other->cell[6],
        other->cell[0] * other->cell[6] * other->cell[11] - other->cell[0] * other->cell[7] * other->cell[10] - other->cell[4] * other->cell[2] * other->cell[11] +
            other->cell[4] * other->cell[3] * other->cell[10] + other->cell[8] * other->cell[2] * other->cell[7] - other->cell[8] * other->cell[3] * other->cell[6],
        other->cell[4] * other->cell[9] * other->cell[15] - other->cell[4] * other->cell[11] * other->cell[13] - other->cell[8] * other->cell[5] * other->cell[15] +
            other->cell[8] * other->cell[7] * other->cell[13] + other->cell[12] * other->cell[5] * other->cell[11] - other->cell[12] * other->cell[7] * other->cell[9],
        -other->cell[0] * other->cell[9] * other->cell[15] + other->cell[0] * other->cell[11] * other->cell[13] + other->cell[8] * other->cell[1] * other->cell[15] -
            other->cell[8] * other->cell[3] * other->cell[13] - other->cell[12] * other->cell[1] * other->cell[11] + other->cell[12] * other->cell[3] * other->cell[9],
        other->cell[0] * other->cell[5] * other->cell[15] - other->cell[0] * other->cell[7] * other->cell[13] - other->cell[4] * other->cell[1] * other->cell[15] +
            other->cell[4] * other->cell[3] * other->cell[13] + other->cell[12] * other->cell[1] * other->cell[7] - other->cell[12] * other->cell[3] * other->cell[5],
        -other->cell[0] * other->cell[5] * other->cell[11] + other->cell[0] * other->cell[7] * other->cell[9] + other->cell[4] * other->cell[1] * other->cell[11] -
            other->cell[4] * other->cell[3] * other->cell[9] - other->cell[8] * other->cell[1] * other->cell[7] + other->cell[8] * other->cell[3] * other->cell[5],
        -other->cell[4] * other->cell[9] * other->cell[14] + other->cell[4] * other->cell[10] * other->cell[13] + other->cell[8] * other->cell[5] * other->cell[14] -
            other->cell[8] * other->cell[6] * other->cell[13] - other->cell[12] * other->cell[5] * other->cell[10] + other->cell[12] * other->cell[6] * other->cell[9],
        other->cell[0] * other->cell[9] * other->cell[14] - other->cell[0] * other->cell[10] * other->cell[13] - other->cell[8] * other->cell[1] * other->cell[14] +
            other->cell[8] * other->cell[2] * other->cell[13] + other->cell[12] * other->cell[1] * other->cell[10] - other->cell[12] * other->cell[2] * other->cell[9],
        -other->cell[0] * other->cell[5] * other->cell[14] + other->cell[0] * other->cell[6] * other->cell[13] + other->cell[4] * other->cell[1] * other->cell[14] -
            other->cell[4] * other->cell[2] * other->cell[13] - other->cell[12] * other->cell[1] * other->cell[6] + other->cell[12] * other->cell[2] * other->cell[5],
        other->cell[0] * other->cell[5] * other->cell[10] - other->cell[0] * other->cell[6] * other->cell[9] - other->cell[4] * other->cell[1] * other->cell[10] +
            other->cell[4] * other->cell[2] * other->cell[9] + other->cell[8] * other->cell[1] * other->cell[6] - other->cell[8] * other->cell[2] * other->cell[5]
    };
    float det = other->cell[0] * inv[0] + other->cell[1] * inv[4] + other->cell[2] * inv[8] + other->cell[3] * inv[12];
    struct mat4 ret_val = IDENTITY_MATRIX;
    if (det != 0.0f)
    {
        float inv_det = 1.0 / det;
        for (int i = 0; i < 16; i++)
        {
            ret_val.cell[i] = inv[i] * inv_det;
        }
    }
    return ret_val;
}