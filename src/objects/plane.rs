use crate::math::*;
use crate::ray::*;

pub struct PlaneUVFunction
{
    f: fn(i: &Float3) -> Float2
}

impl PlaneUVFunction
{
    pub fn empty() -> Self
    {
        PlaneUVFunction
        {
            f: |_: &Float3| Float2::zero()
        }
    }

    pub fn xz() -> Self
    {
        PlaneUVFunction
        {
            f: |i: &Float3| {
                Float2::from_xy(i.x, i.z)
            }
        }
    }

    pub fn xy() -> Self
    {
        PlaneUVFunction
        {
            f: |i: &Float3| {
                Float2::from_xy(i.x, i.y)
            }
        }
    }

    pub fn zy() -> Self
    {
        PlaneUVFunction
        {
            f: |i: &Float3| {
                Float2::from_xy(i.z, i.y)
            }
        }
    }
}

pub struct Plane
{
    pub obj_idx: i32,
    pub mat_idx: i32,
    pub normal: Float3,
    pub distance: f32,
    pub uv_function: PlaneUVFunction
}

impl RayHittableObject for Plane
{
    fn intersect(&self, ray: &mut Ray) {
        let t = -(dot(&ray.origin, &self.normal) + self.distance) / (dot(&ray.direction, &self.normal));
        if t < ray.t && t > 0.0
        {
            ray.t = t;
            ray.obj_idx = self.obj_idx;
            ray.obj_type = RayHittableObjectType::Plane;
        }
    }

    fn get_normal(&self, _: &Float3) -> Float3 {
        return self.normal;
    }

    fn get_uv(&self, i: &Float3) -> Float2
    {
        return (self.uv_function.f)(i)
    }

    fn is_occluded(&self, ray: &Ray) -> bool {
        return false;
    }
}

impl Plane
{
    pub fn new(obj_idx: i32, mat_idx: i32, normal: Float3, distance: f32, uv_function: PlaneUVFunction) -> Self
    {
        Plane{
            obj_idx,
            mat_idx,
            normal,
            distance,
            uv_function
        }
    }
}