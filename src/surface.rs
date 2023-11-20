use image::GenericImageView;

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

    pub fn load_from_file(file: &std::path::Path) -> Self
    {
        let img = image::open(file).expect("Image not found");

        let mut pixels: Vec<u32> = vec![0; (img.width() * img.height()) as usize];
        for (x,y, value) in img.pixels()
        {
            let r = value[0] as u32;
            let g = value[1] as u32;
            let b = value[2] as u32;
            pixels[(y * img.height() + x) as usize] = (r << 16) + (g << 8) + b;
        }

        return Surface { pixels }
    }

    pub fn set_pixel(&mut self, x: usize, y: usize, value: u32)
    {
        self.pixels[y * SCRWIDTH + x] = value;
    }
}