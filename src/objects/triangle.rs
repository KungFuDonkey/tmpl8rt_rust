use crate::math::*;
use crate::ray::*;


#[derive(Clone, Copy)]
pub struct Triangle
{
    pub tri_idx: u32,
    pub vertex0: Float3,
    pub vertex1: Float3,
    pub vertex2: Float3
}


pub fn intersect_triangle(triangle: &Triangle, ray: &mut Ray) -> bool
{
    ray.triangle_intersection_tests += 1;

    let edge1 = triangle.vertex1 - triangle.vertex0;
    let edge2 = triangle.vertex2 - triangle.vertex0;
    let h = cross(&ray.direction, &edge2);
    let a = dot(&edge1, &h);
    if a  > -0.0001 && a < 0.0001
    {
        return false;
    }
    let s = ray.origin - triangle.vertex0;
    let f = 1.0 / a;
    let u = f * dot(&s, &h);
    if u < 0.0 || u > 1.0
    {
        return false;
    }
    let q = cross(&s, &edge1);
    let v = f * dot(&ray.direction, &q);
    if v < 0.0 || u + v > 1.0
    {
        return false;
    }
    let t = f * dot(&edge2, &q);
    if t > 0.0001 && t < ray.t
    {
        ray.t = t;
        ray.sub_obj_idx = triangle.tri_idx as usize;
        return true;
    }
    return false;
}