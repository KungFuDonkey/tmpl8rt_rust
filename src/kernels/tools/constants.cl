
#define MAX_UINT 4294967295

#define SPHERE_OBJ_TYPE (0)
#define PLANE_OBJ_TYPE (1 << 30)
#define QUAD_OBJ_TYPE (2 << 30)
#define MESH_OBJ_TYPE (3 << 30)

#define OBJ_TYPE_MASK (3 << 30)
#define OBJ_ID_MASK (MAX_UINT - OBJ_TYPE_MASK)