use crate::math::*;
use crate::ray::*;

pub struct Quad
{
    pub size: f32,
    pub t: Mat4,
    pub inv_t: Mat4,
    pub obj_idx: usize,
    pub mat_idx: usize,
}

impl Quad
{
    pub fn new(obj_idx: usize, mat_idx: usize, size: f32, transform: &Mat4) -> Quad
    {
        Quad{
            obj_idx,
            mat_idx,
            size: size * 0.5,
            t: transform.clone(),
            inv_t: transform.inverted_no_scale()
        }
    }

    pub fn random_points(&self, num: usize, seed: &mut u32) -> Vec<Float3>
    {
        let mut vector: Vec<Float3> = Vec::with_capacity(num);
        let size = self.size;
        let corner1 = transform_position(&Float3::from_xyz(-size, 0.0, -size), &self.t);
        let corner2 = transform_position(&Float3::from_xyz(size, 0.0, -size), &self.t);
        let corner3 = transform_position(&Float3::from_xyz(-size, 0.0, size), &self.t);

        for _ in 0..num
        {
            let r1 = random_float_s(seed);
            let r2 = random_float_s(seed);
            vector.push(corner1 + r2 * (corner2 - corner1) + r1 * (corner3 - corner1));
        }
        return vector;
    }

    pub fn center_point(&self) -> Float3
    {
        let r1 = 0.5;
        let r2 = 0.5;
        let size = self.size;
        let corner1 = transform_position(&Float3::from_xyz(-size, 0.0, -size), &self.t);
        let corner2 = transform_position(&Float3::from_xyz(size, 0.0, -size), &self.t);
        let corner3 = transform_position(&Float3::from_xyz(-size, 0.0, size), &self.t);
        return corner1 + r2 * (corner2 - corner1) + r1 * (corner3 - corner1);
    }
}

impl RayHittableObject for Quad
{
    fn intersect(&self, ray: &mut Ray) {
        let o_y: f32 = self.inv_t.cell[4] * ray.origin.x + self.inv_t.cell[5] * ray.origin.y + self.inv_t.cell[6] * ray.origin.z + self.inv_t.cell[7];
        let d_y: f32 = self.inv_t.cell[4] * ray.direction.x + self.inv_t.cell[5] * ray.direction.y + self.inv_t.cell[6] * ray.direction.z;
        let t = o_y / -d_y;
        if t < ray.t && t > 0.0
        {
            let o_x: f32 = self.inv_t.cell[0] * ray.origin.x + self.inv_t.cell[1] * ray.origin.y + self.inv_t.cell[2] * ray.origin.z + self.inv_t.cell[3];
            let o_z: f32 = self.inv_t.cell[8] * ray.origin.x + self.inv_t.cell[9] * ray.origin.y + self.inv_t.cell[10] * ray.origin.z + self.inv_t.cell[11];
            let d_x: f32 = self.inv_t.cell[0] * ray.direction.x + self.inv_t.cell[1] * ray.direction.y + self.inv_t.cell[2] * ray.direction.z;
            let d_z: f32 = self.inv_t.cell[8] * ray.direction.x + self.inv_t.cell[9] * ray.direction.y + self.inv_t.cell[10] * ray.direction.z;
            let i_x = o_x + t * d_x;
            let i_z = o_z + t * d_z;
            if i_x > -self.size && i_x < self.size && i_z > -self.size && i_z < self.size
            {
                ray.t = t;
                ray.obj_idx = self.obj_idx;
                ray.obj_type = RayHittableObjectType::Quad
            }
        }
    }

    fn get_normal(&self, _: &Ray, i: &Float3) -> Float3 {
        Float3::from_xyz(-self.t.cell[1], -self.t.cell[5], -self.t.cell[9])
    }

    fn get_uv(&self, i: &Float3) -> Float2 {
        Float2::zero()
    }

    fn is_occluded(&self, ray: &Ray) -> bool {
        let o_y: f32 = self.inv_t.cell[4] * ray.origin.x + self.inv_t.cell[5] * ray.origin.y + self.inv_t.cell[6] * ray.origin.z + self.inv_t.cell[7];
        let d_y: f32 = self.inv_t.cell[4] * ray.direction.x + self.inv_t.cell[5] * ray.direction.y + self.inv_t.cell[6] * ray.direction.z;
        let t = o_y / -d_y;
        if t < ray.t && t > 0.0
        {
            let o_x: f32 = self.inv_t.cell[0] * ray.origin.x + self.inv_t.cell[1] * ray.origin.y + self.inv_t.cell[2] * ray.origin.z + self.inv_t.cell[3];
            let o_z: f32 = self.inv_t.cell[8] * ray.origin.x + self.inv_t.cell[9] * ray.origin.y + self.inv_t.cell[10] * ray.origin.z + self.inv_t.cell[11];
            let d_x: f32 = self.inv_t.cell[0] * ray.direction.x + self.inv_t.cell[1] * ray.direction.y + self.inv_t.cell[2] * ray.direction.z;
            let d_z: f32 = self.inv_t.cell[8] * ray.direction.x + self.inv_t.cell[9] * ray.direction.y + self.inv_t.cell[10] * ray.direction.z;
            let i_x = o_x + t * d_x;
            let i_z = o_z + t * d_z;
            return i_x > -self.size && i_x < self.size && i_z > -self.size && i_z < self.size;
        }

        return false;
    }
}