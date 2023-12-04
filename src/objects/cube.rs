use crate::math::*;
use crate::ray::*;

pub struct Cube
{
    pub obj_idx: i32,
    pub mat_idx: i32,
    pub b: [Float3; 2],
    pub m: Mat4,
    pub inv_m: Mat4
}

impl RayHittableObject for Cube
{
    fn intersect(&self, ray: &mut Ray) {
        let origin = transform_position(&ray.origin, &self.inv_m);
        let direction = transform_vector(&ray.direction, &self.inv_m);
        let rdx = 1.0 / direction.x;
        let rdy = 1.0 / direction.y;
        let rdz = 1.0 / direction.z;
        let signx = (direction.x < 0.0) as usize;
        let signy = (direction.y < 0.0) as usize;
        let signz = (direction.z < 0.0) as usize;

        let mut tmin = (self.b[signx].x - origin.x) * rdx;
        let mut tmax = (self.b[1 - signx].x - origin.x) * rdx;
        let tymin = (self.b[signy].y - origin.y) * rdy;
        let tymax = (self.b[1 - signy].y - origin.y) * rdy;

        if tmin > tymax || tymin > tmax
        {
            return;
        }

        tmin = tmin.max(tymin);
        tmax = tmax.min(tymax);
        let tzmin = (self.b[signz].z - origin.z) * rdz;
        let tzmax = (self.b[1 - signz].z - origin.z) * rdz;

        if tmin > tzmax || tzmin > tmax
        {
            return;
        }

        tmin = tmin.max(tzmin);
        tmax = tmax.min(tzmax);
        if tmin > 0.0
        {
            if tmin < ray.t
            {
                ray.t = tmin;
                ray.obj_idx = self.obj_idx;
                ray.obj_type = RayHittableObjectType::Cube;
            }
        }
        else if tmax > 0.0
        {
            if tmax < ray.t
            {
                ray.t = tmax;
                ray.obj_idx = self.obj_idx;
                ray.obj_type = RayHittableObjectType::Cube;
            }
        }
    }

    fn get_normal(&self, i: &Float3) -> Float3 {
        let obj_i = transform_position(&i, &self.inv_m);
        let mut n = Float3::from_xyz(-1.0, 0.0, 0.0);
        let d0 = (obj_i.x - self.b[0].x).abs();
        let d1 = (obj_i.x - self.b[1].x).abs();
        let d2 = (obj_i.y - self.b[0].x).abs();
        let d3 = (obj_i.y - self.b[1].x).abs();
        let d4 = (obj_i.z - self.b[0].x).abs();
        let d5 = (obj_i.z - self.b[1].x).abs();
        let mut min_dist = d0;
        if d1 < min_dist
        {
            min_dist = d1;
            n.x = 1.0;
        }
        if d2 < min_dist
        {
            min_dist = d2;
            n = Float3::from_xyz(0.0, -1.0, 0.0);
        }
        if d3 < min_dist
        {
            min_dist = d3;
            n = Float3::from_xyz(0.0, 1.0, 0.0);
        }
        if d4 < min_dist
        {
            min_dist = d4;
            n = Float3::from_xyz(0.0, 0.0, -1.0);
        }
        if d5 < min_dist
        {
            n = Float3::from_xyz(0.0, 0.0, 1.0);
        }

        return transform_vector(&n, &self.m);
    }

    fn get_uv(&self, i: &Float3) -> Float2
    {
        Float2::zero()
    }

    fn is_occluded(&self, ray: &Ray) -> bool {
        let origin = transform_position(&ray.origin, &self.inv_m);
        let direction = transform_vector(&ray.direction, &self.inv_m);
        let rdx = 1.0 / direction.x;
        let rdy = 1.0 / direction.y;
        let rdz = 1.0 / direction.z;
        let t1 = (self.b[0].x - origin.x) * rdx;
        let t2 = (self.b[1].x - origin.x) * rdx;
        let t3 = (self.b[0].y - origin.y) * rdy;
        let t4 = (self.b[1].y - origin.y) * rdy;
        let t5 = (self.b[0].z - origin.z) * rdz;
        let t6 = (self.b[1].z - origin.z) * rdz;
        let tmin = t1.min( t2 ).max( t3.min( t4 ) ).max( t5.min( t6 ) );
        let tmax = t1.max( t2 ).min( t3.max( t4 ) ).min( t5.max( t6 ) );
        return tmax > 0.0 && tmin < tmax && tmin < ray.t;
    }
}

impl Cube
{
    pub fn new(obj_idx: i32, mat_idx: i32, pos: Float3, size: Float3) -> Self
    {
        Cube{
            obj_idx,
            mat_idx,
            b:[
                pos - size * 0.5,
                pos + size * 0.5
            ],
            m: Mat4::identity_matrix(),
            inv_m: Mat4::identity_matrix().inverted()
        }
    }
}