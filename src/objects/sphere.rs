use crate::math::*;
use crate::ray::*;

pub struct Sphere
{
    pub position: Float3,
    pub radius: f32,
    pub radius2: f32,
    pub inv_radius: f32,
    pub obj_idx: i32,
    pub mat_idx: i32,
}

impl RayHittableObject for Sphere
{
    fn intersect(&self, ray: &mut Ray) {
        let oc = ray.origin - self.position;
        let b = dot(&oc, &ray.direction);
        let c = dot(&oc, &oc) - self.radius2;
        let mut d = b * b - c;
        if d <= 0.0
        {
            return;
        }

        d = d.sqrt();
        let mut t = -b - d;
        let mut hit = t < ray.t && t > 0.0;
        if hit
        {
            ray.t = t;
            ray.obj_idx = self.obj_idx;
            ray.obj_type = RayHittableObjectType::Sphere;
            return;
        }
        if c > 0.0
        {
            return;
        }

        t = d - b;
        hit = t < ray.t && t > 0.0;
        if hit
        {
            ray.t = t;
            ray.obj_idx = self.obj_idx;
            ray.obj_type = RayHittableObjectType::Sphere;
            return;
        }
    }

    fn get_normal(&self, i: &Float3) -> Float3 {
        return (*i - self.position) * self.inv_radius;
    }

    fn get_uv(&self, i: &Float3) -> Float2
    {
        Float2::zero()
    }

    fn is_occluded(&self, ray: &Ray) -> bool {
        let oc = ray.origin - self.position;
        let b = dot(&oc, &ray.direction);
        let c = dot(&oc, &oc) - self.radius2;

        let d = b * b - c;
        if d <= 0.0
        {
            return false
        }

        let t = -b - d.sqrt();
        return t < ray.t && t > 0.0;
    }
}

impl Sphere
{
    pub fn new(obj_idx: i32, mat_idx: i32, position: Float3, radius: f32) -> Self
    {
        Sphere
        {
            obj_idx,
            mat_idx,
            position,
            radius,
            radius2: radius * radius,
            inv_radius: 1.0 / radius
        }
    }
}