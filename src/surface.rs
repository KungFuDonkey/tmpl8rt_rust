

pub const SCRWIDTH: usize = 1024;
pub const SCRHEIGHT: usize = 640;

pub struct Surface
{
    pub pixels: Vec<u32>
}

impl Surface
{
    pub fn new() -> Self
    {
        Surface {
            pixels: vec![0; SCRWIDTH * SCRHEIGHT]
        }
    }

    pub fn set_pixel(&mut self, x: usize, y: usize, value: u32)
    {
        self.pixels[y * SCRWIDTH + x] = value;
    }
}