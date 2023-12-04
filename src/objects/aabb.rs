use crate::math::*;
use crate::ray::*;


#[derive(Clone, Copy)]
pub struct AABB
{
    pub min_bound: Float3,
    pub max_bound: Float3
}

impl AABB
{
    pub fn from_empty() -> Self
    {
        AABB{
            min_bound: Float3::from_a(1e34),
            max_bound: Float3::from_a(-1e34)
        }
    }

    pub fn from_bounds(min_bound: &Float3, max_bound: &Float3) -> Self
    {
        AABB{
            min_bound: *min_bound,
            max_bound: *max_bound
        }
    }

    pub fn reset(&mut self)
    {
        self.min_bound = Float3::from_a(1e34);
        self.max_bound = Float3::from_a(-1e34);
    }

    pub fn contains(&self, vector: &Float3) -> bool
    {
        let p = *vector;
        let va = p - self.min_bound;
        let vb = self.max_bound - p;
        return va.x >= 0.0 && va.y >= 0.0 && va.z >= 0.0 && vb.x >= 0.0 && vb.y >= 0.0 && vb.z >= 0.0;
    }

    pub fn grow_aabb(&mut self, bb: &AABB)
    {
        self.min_bound = self.min_bound.min(&bb.min_bound);
        self.max_bound = self.max_bound.max(&bb.max_bound);
    }

    pub fn grow(&mut self, vector: &Float3)
    {
        self.min_bound = self.min_bound.min(vector);
        self.max_bound = self.max_bound.max(vector);
    }

    pub fn union(&self, other: &AABB) -> Self
    {
        let min_bound = self.min_bound.min(&other.min_bound);
        let max_bound = self.max_bound.max(&other.max_bound);

        AABB{
            min_bound,
            max_bound
        }
    }

    pub fn intersection(&self, other: &AABB) -> Self
    {
        let min_bound = self.min_bound.max(&other.min_bound);
        let max_bound = self.max_bound.min(&other.max_bound);

        AABB{
            min_bound,
            max_bound
        }
    }

    pub fn extend_x(&self) -> f32
    {
        return self.max_bound.x - self.min_bound.x;
    }

    pub fn extend_y(&self) -> f32
    {
        return self.max_bound.y - self.min_bound.y;
    }

    pub fn extend_z(&self) -> f32
    {
        return self.max_bound.z - self.min_bound.z;
    }

    pub fn area(&self) -> f32
    {
        let e = self.max_bound - self.min_bound;
        let zero: f32 = 0.0;
        return zero.max(e.x * e.y + e.x * e.z + e.y * e.z);
    }

    pub fn set_bounds(&mut self, max_bound: &Float3, min_bound: &Float3)
    {
        self.min_bound = *min_bound;
        self.max_bound = *max_bound;
    }

    pub fn center(&self, axis: usize) -> f32
    {
        match axis
        {
            0 => (self.min_bound.x + self.max_bound.x) * 0.5,
            1 => (self.min_bound.y + self.max_bound.y) * 0.5,
            2 => (self.min_bound.z + self.max_bound.z) * 0.5,
            _ => panic!("Outside of axis range")
        }
    }

    pub fn minimal(&self, axis: usize) -> f32
    {
        match axis
        {
            0 => self.min_bound.x,
            1 => self.min_bound.y,
            2 => self.min_bound.z,
            _ => panic!("Outside of axis range")
        }
    }

    pub fn maximal(&self, axis: usize) -> f32
    {
        match axis
        {
            0 => self.max_bound.x,
            1 => self.max_bound.y,
            2 => self.max_bound.z,
            _ => panic!("Outside of axis range")
        }
    }
}

pub fn intersect_aabb(ray: &Ray, aabb: &AABB) -> f32
{
    let tx1 = (aabb.min_bound.x - ray.origin.x) * ray.r_direction.x;
    let tx2 = (aabb.max_bound.x - ray.origin.x) * ray.r_direction.x;
    let mut tmin = tx1.min(tx2);
    let mut tmax = tx1.max(tx2);
    let ty1 = (aabb.min_bound.y - ray.origin.y) * ray.r_direction.y;
    let ty2 = (aabb.max_bound.y - ray.origin.y) * ray.r_direction.y;
    tmin = tmin.max(ty1.min(ty2));
    tmax = tmax.min(ty1.max(ty2));
    let tz1 = (aabb.min_bound.z - ray.origin.z) * ray.r_direction.z;
    let tz2 = (aabb.max_bound.z - ray.origin.z) * ray.r_direction.z;
    tmin = tmin.max(tz1.min(tz2));
    tmax = tmax.min(tz1.max(tz2));
    if tmax >= tmin && tmin < ray.t && tmax > 0.0
    {
        return tmin;
    }
    return 1e30;
}