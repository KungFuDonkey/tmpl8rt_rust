use crate::math::*;

#[derive(Debug, Clone, Copy)]
pub enum RayHittableObjectType
{
    None,
    Sphere,
    Plane,
    Cube,
    Torus,
    Quad,
    Mesh,
    Bvh
}

#[repr(C, align(64))]
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    pub origin: Float3,
    pub direction: Float3,
    pub r_direction: Float3,
    pub t: f32,
    pub obj_idx: usize,
    pub sub_obj_idx: usize,
    pub obj_type: RayHittableObjectType,
    pub inside: bool
}

impl Ray {

    #![allow(dead_code)]
    pub fn new() -> Self
    {
        Ray {
            origin: Float3::zero(),
            direction: Float3::zero(),
            r_direction: Float3::zero(),
            t: 1.0e34,
            obj_idx: usize::MAX,
            sub_obj_idx: 0,
            obj_type: RayHittableObjectType::None,
            inside: false
        }
    }

    pub fn directed(origin: Float3, direction: Float3) -> Self
    {
        Ray {
            origin: Float3::from_xyz(origin.x, origin.y, origin.z),
            direction: Float3::from_xyz(direction.x, direction.y, direction.z),
            r_direction: Float3::from_xyz(1.0 / direction.x, 1.0 / direction.y, 1.0 / direction.z),
            t: 1.0e34,
            obj_idx: usize::MAX,
            sub_obj_idx: 0,
            obj_type: RayHittableObjectType::None,
            inside: false
        }
    }

    #[allow(dead_code)]
    pub fn directed_distance(origin: Float3, direction: Float3, distance: f32) -> Self
    {
        Ray {
            origin: Float3::from_xyz(origin.x, origin.y, origin.z),
            direction: Float3::from_xyz(direction.x, direction.y, direction.z),
            r_direction: Float3::from_xyz(1.0 / direction.x, 1.0 / direction.y, 1.0 / direction.z),
            t: distance,
            obj_idx: usize::MAX,
            sub_obj_idx: 0,
            obj_type: RayHittableObjectType::None,
            inside: false
        }
    }

    pub fn intersection_point(&self) -> Float3
    {
        return self.origin + self.direction * self.t;
    }
}

pub trait RayHittableObject
{
    fn intersect(&self, ray: &mut Ray);

    fn get_normal(&self, ray: &Ray, i: &Float3) -> Float3;

    fn get_uv(&self, i: &Float3) -> Float2;

    fn is_occluded(&self, ray: &Ray) -> bool;
}

pub trait RayOccluder
{
    fn is_occluded(&self, ray: &Ray);
}