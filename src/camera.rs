use crate::input::Input;
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

    pub fn handle_input(&mut self, input: &Input, delta_time: f32)
    {
        if !input.window_has_focus()
        {
            return;
        }

        let speed = 2.5 * delta_time;
        let dir = self.target - self.position;
        let mut ahead = normalize_f3( &dir );
        let tmp_up = Float3::from_xyz(0.0, 1.0, 0.0);
        let mut right = normalize_f3(&cross(&tmp_up, &ahead));
        let mut up = normalize_f3(&cross( &ahead, &right ));
        let mut changed = false;

        if input.is_key_down(glfw::Key::A)
        {
            self.position -= right * speed * 2.0;
            changed = true;
        }
        if input.is_key_down(glfw::Key::D)
        {
            self.position += right * speed * 2.0;
            changed = true;
        }
        if input.is_key_down(glfw::Key::W)
        {
            self.position += ahead * speed * 2.0;
            changed = true;
        }
        if input.is_key_down(glfw::Key::S)
        {
            self.position -= ahead * speed * 2.0;
            changed = true;
        }
        if input.is_key_down(glfw::Key::R)
        {
            self.position += up * speed * 2.0;
            changed = true;
        }
        if input.is_key_down(glfw::Key::F)
        {
            self.position -= up * speed * 2.0;
            changed = true;
        }

        self.target = self.position + ahead;
        if input.is_key_down(glfw::Key::Up)
        {
            self.target += up * speed;
            changed = true;
        }
        if input.is_key_down(glfw::Key::Down)
        {
            self.target -= up * speed;
            changed = true;
        }
        if input.is_key_down(glfw::Key::Left)
        {
            self.target -= right * speed;
            changed = true;
        }
        if input.is_key_down(glfw::Key::Right)
        {
            self.target += right * speed;
            changed = true;
        }

        if !changed
        {
            return;
        }

        ahead = normalize_f3( &(self.target - self.position) );
        up = normalize_f3( &cross( &ahead, &right ) );
        right = normalize_f3( &cross( &up, &ahead ) );
        self.top_left = self.position + ahead * 2.0 - right * self.aspect_ratio  + up;
        self.top_right = self.position + ahead * 2.0  + right * self.aspect_ratio + up;
        self.bottom_left = self.position + ahead * 2.0 - right * self.aspect_ratio - up;
    }
}