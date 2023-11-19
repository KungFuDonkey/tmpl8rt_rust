use crate::math::*;
use crate::scene::Ray;
use crate::surface::{SCRHEIGHT, SCRWIDTH};

pub struct Camera
{
    pub position: Float3,
    pub target: Float3,
    pub top_left: Float3,
    pub top_right: Float3,
    pub bottom_left: Float3,
    pub aspect_ratio: f32
}

impl Camera
{
    pub fn new() -> Self
    {
        let aspect = (SCRWIDTH as f32) / (SCRHEIGHT as f32);
        Camera {
            position: Float3::from_xyz(0.0,0.0,-2.0),
            target: Float3::from_xyz(0.0, 0.0, -1.0),
            top_left: Float3::from_xyz(-aspect, 1.0, 0.0),
            top_right: Float3::from_xyz(aspect, 1.0, 0.0),
            bottom_left: Float3::from_xyz(-aspect, -1.0, 0.0),
            aspect_ratio: aspect
        }
    }

    pub fn set_aspect_ratio(&mut self, aspect: f32)
    {
        self.aspect_ratio = aspect;
        self.top_left = Float3::from_xyz(-aspect, 1.0, 0.0);
        self.top_right = Float3::from_xyz(aspect, 1.0, 0.0);
        self.bottom_left = Float3::from_xyz(-aspect, -1.0, 0.0);
    }

    pub fn get_primary_ray(&self, x: f32, y: f32) -> Ray
    {
        let u = x / (SCRWIDTH as f32);
        let v = y / (SCRHEIGHT as f32);
        let p = self.top_left + (self.top_right - self.top_left) * u + (self.bottom_left - self.top_left) * v;
        let d = p - self.position;
        return Ray::directed(self.position, normalize_f3(&d))
    }

    pub fn get_primary_ray_indexed(&self, index: usize) -> Ray
    {
        let y = index / SCRWIDTH;
        let x = index - y * SCRWIDTH;
        return self.get_primary_ray(x as f32, y as f32);
    }
}