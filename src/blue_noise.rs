use image::GenericImageView;
use crate::opencl::*;


pub fn load_blue_noise_from_file(cl: &OpenCL, file: std::path::PathBuf) -> OpenCLBuffer<u8>
{
    let img = image::open(file).expect("Blue noise not found");

    let mut pixels: Vec<u8> = vec![0; (img.width() * img.height()) as usize];
    for (x,y, value) in img.pixels()
    {
        let r = value[0] as u8;
        pixels[(y * img.height() + x) as usize] = r;
    }

    return OpenCLBuffer::read_only(cl, pixels);
}