use std::f32::consts::PI;
use crate::math::*;
use crate::material::*;
use crate::ray::*;
use crate::objects::sphere::*;
use crate::objects::plane::*;
use crate::objects::cube::*;
use crate::objects::torus::*;
use crate::objects::quad::*;
use crate::objects::mesh::*;
use crate::objects::bvh::*;

pub struct Scene
{
    spheres: Vec<Sphere>,
    planes: Vec<Plane>,
    cubes: Vec<Cube>,
    tori: Vec<Torus>,
    quads: Vec<Quad>,
    meshes: Vec<Mesh>,
    bvhs: Vec<BVH>,
    materials: Vec<Material>,
    animation_time: f32,
}

impl Scene
{
    pub fn new() -> Self
    {
        let mut torus = Torus::new(0, 0, 0.8, 0.25);
        let translation = Float3::from_xyz(-0.25, 0.0, 2.0);
        torus.t = Mat4::translate(&translation) * Mat4::rotate_x(PI / 4.0);
        torus.inv_t = torus.t.inverted();

        let transform = Mat4::translate( &Float3::from_xyz(0.0, 0.0, 2.5)) * Mat4::scale(0.5);
        let triangle_mesh = Mesh::triangle(0, 6, transform);
        let triangle_bvh = BVH::from_mesh(&triangle_mesh);

        let transform = Mat4::translate( &Float3::from_xyz(0.0, 0.0, -1.0)) * Mat4::scale(0.5);
        let unity_mesh = Mesh::from_tri_file(1, 8, transform, std::path::Path::new("./assets/unity.tri"));
        let unity_bvh = BVH::from_mesh(&unity_mesh);

        Scene{
            spheres: vec![
                Sphere::new(0, 11, Float3::from_a(0.0), 0.3),
                Sphere::new(1, 0, Float3::from_xyz( 0.0, 2.5, -3.07 ), 8.0)
            ],
            planes: vec![
                Plane::new(0, 4, Float3::from_xyz(1.0, 0.0, 0.0), 3.0, PlaneUVFunction::zy()),
                Plane::new(1, 5, Float3::from_xyz(-1.0, 0.0, 0.0), 2.99, PlaneUVFunction::zy()),
                Plane::new(2, 2, Float3::from_xyz(0.0, 1.0, 0.0), 1.0, PlaneUVFunction::xz()),
                Plane::new(3, 1, Float3::from_xyz(0.0, -1.0, 0.0), 2.0, PlaneUVFunction::empty()),
                Plane::new(4, 1, Float3::from_xyz(0.0, 0.0, 1.0), 3.0, PlaneUVFunction::empty()),
                Plane::new(5, 3, Float3::from_xyz(0.0, 0.0, -1.0), 3.99, PlaneUVFunction::xy()),
            ],
            cubes: vec![
                Cube::new(0, 10, Float3::zero(), Float3::from_a(1.15))
            ],
            tori: vec![
                torus
            ],
            quads: vec![
                Quad::new(0, 6, 0.5, &Mat4::translate(&Float3::from_xyz(-1.0, 1.5, -1.0))),
                Quad::new(1, 7, 0.5, &Mat4::translate(&Float3::from_xyz(1.0, 1.5, -1.0))),
                Quad::new(2, 8, 0.5, &Mat4::translate(&Float3::from_xyz(1.0, 1.5, 1.0))),
                Quad::new(3, 0, 0.5, &Mat4::translate(&Float3::from_xyz(-1.0, 1.5, 1.0))),
            ],
            meshes: vec![
                //triangle_mesh
            ],
            bvhs: vec![
                triangle_bvh,
                unity_bvh
            ],
            materials: vec![
                linear_color(Float3::from_a(1.0)),
                linear_color(Float3::from_a(0.93)),
                checkerboard_material(),
                logo_material(),
                red_material(),
                blue_material(),
                linear_color(Float3::from_xyz(1.0, 0.0, 0.0)),
                linear_color(Float3::from_xyz(0.0, 1.0, 0.0)),
                linear_color(Float3::from_xyz(0.0, 0.0, 1.0)),
                reflective_material(linear_color_simple(Float3::from_a(1.0)), 0.02),
                fully_reflective_material(linear_color_simple(Float3::from_xyz(1.0, 0.0, 0.0))),
                refractive_material_with_exit(linear_color_simple(Float3::from_xyz(1.0, 0.0, 0.0)), 1.5, 10.0),
            ],
            animation_time: 0.0
        }
    }

    pub fn intersect_scene(&self, ray: &mut Ray)
    {
        for sphere in &self.spheres
        {
            sphere.intersect(ray);
        }

        for plane in &self.planes
        {
            plane.intersect(ray);
        }

        for cube in &self.cubes
        {
            cube.intersect(ray);
        }

        for torus in &self.tori
        {
            torus.intersect(ray);
        }

        for quad in &self.quads
        {
            quad.intersect(ray);
        }

        for mesh in &self.meshes
        {
            mesh.intersect(ray);
        }

        for bvhs in &self.bvhs
        {
            bvhs.intersect(ray);
        }
    }

    pub fn is_occluded(&self, ray: &Ray) -> bool
    {
        for sphere in &self.spheres
        {
            if sphere.is_occluded(ray)
            {
                return true;
            }
        }

        // skip planes

        for cube in &self.cubes
        {
            if cube.is_occluded(ray)
            {
                return true;
            }
        }

        for quad in &self.quads
        {
            if quad.is_occluded(ray)
            {
                return true;
            }
        }

        for torus in &self.tori
        {
            if torus.is_occluded(ray)
            {
                return true;
            }
        }

        for mesh in &self.meshes
        {
            if mesh.is_occluded(ray)
            {
                return true;
            }
        }

        return false;
    }


    fn get_quad_light_color(&self, quad: &Quad, light_point: &Float3) -> Float3
    {
        let simple_material = get_simple_material(&self.materials[quad.mat_idx as usize]);
        match simple_material
        {
            SimpleMaterial::LinearColorMaterial(color) => *color,
            SimpleMaterial::UV(uv_material) => get_color_from_uv_material(uv_material, &quad.get_uv(&light_point)),
            SimpleMaterial::BlackMaterial => Float3::zero()
        }
    }

    pub fn direct_lighting_soft(&self, point: &Float3, normal: &Float3, area_sample_size: usize, seed: &mut u32) -> Float3
    {
        let total_sample_points = self.quads.len() * area_sample_size;
        let sample_strength = 1.0 / (total_sample_points as f32);
        let mut lighting = Float3::zero();
        for quad in &self.quads
        {
            let light_points = quad.random_points(area_sample_size, seed);
            for light_point in light_points
            {
                let ray_dir = light_point - *point;
                let ray_dir_n = normalize(&ray_dir);
                let origin = (*point) + (EPSILON * ray_dir_n);
                let ray = Ray::directed_distance(origin, ray_dir_n, length(&ray_dir) - (2.0 * EPSILON));

                if self.is_occluded(&ray)
                {
                    continue;
                }

                let light_color = self.get_quad_light_color(&quad, &light_point);
                // take distance to the light?
                lighting += sample_strength * dot(&ray_dir_n, &normal) * light_color;
            }
        }

        return lighting;
    }

    pub fn direct_lighting_hard(&self, point: &Float3, normal: &Float3) -> Float3
    {
        let total_sample_points = (self.quads.len() as u32);
        let light_strength = 1.0 / (total_sample_points as f32);
        let mut lighting = Float3::zero();
        for quad in &self.quads
        {
            let light_point = quad.center_point();
            let ray_dir = light_point - *point;
            let ray_dir_n = normalize(&ray_dir);
            let origin = (*point) + (EPSILON * ray_dir_n);
            let ray = Ray::directed_distance(origin, ray_dir_n, length(&ray_dir) - (2.0 * EPSILON));

            if self.is_occluded(&ray)
            {
                continue;
            }

            let light_color = self.get_quad_light_color(&quad, &light_point);

            // take distance to the light?
            lighting += light_strength * dot(&ray_dir_n, &normal) * light_color;
        }

        return lighting;
    }


    pub fn get_normal(&self, ray: &Ray, i: &Float3, wo: &Float3) -> Float3
    {
        let obj_idx = ray.obj_idx;
        if obj_idx == -1
        {
            println!("ERROR: obj_idx not set or no object was hit");
            return Float3::zero();
        }

        let obj_idx = obj_idx as usize;
        let obj_type = ray.obj_type;

        let normal = match obj_type
        {
            RayHittableObjectType::None =>
            {
                println!("ERROR: tried to get normal of non object");
                Float3::zero()
            },
            RayHittableObjectType::Sphere =>
            {
                self.spheres[obj_idx].get_normal(i)
            },
            RayHittableObjectType::Plane =>
            {
                self.planes[obj_idx].get_normal(i)
            },
            RayHittableObjectType::Cube =>
            {
                self.cubes[obj_idx].get_normal(i)
            },
            RayHittableObjectType::Torus =>
            {
                self.tori[obj_idx].get_normal(i)
            },
            RayHittableObjectType::Mesh =>
            {
                self.meshes[obj_idx].triangle_normals[ray.sub_obj_idx]
            },
            RayHittableObjectType::Quad =>
            {
                self.quads[obj_idx].get_normal(i)
            }
            RayHittableObjectType::Bvh =>
            {
                self.bvhs[obj_idx].get_normal(i)
            }
        };

        if dot(&normal, &wo) > 0.0
        {
            return -normal;
        }
        return normal;
    }

    pub fn get_material(&self, ray: &Ray) -> &Material
    {
        let obj_idx = ray.obj_idx as usize;
        let obj_type = ray.obj_type;
        let mat_idx = match obj_type
        {
            RayHittableObjectType::None =>
            {
                -1
            },
            RayHittableObjectType::Sphere =>
            {
                self.spheres[obj_idx as usize].mat_idx
            },
            RayHittableObjectType::Plane =>
            {
                self.planes[obj_idx as usize].mat_idx
            },
            RayHittableObjectType::Cube =>
            {
                self.cubes[obj_idx as usize].mat_idx
            },
            RayHittableObjectType::Quad =>
            {
                self.quads[obj_idx as usize].mat_idx
            },
            RayHittableObjectType::Mesh =>
            {
                self.meshes[obj_idx as usize].mat_idx
            },
            RayHittableObjectType::Torus =>
            {
                self.tori[obj_idx as usize].mat_idx
            },
            RayHittableObjectType::Bvh =>
            {
                self.bvhs[obj_idx as usize].mat_idx
            }
        };
        &self.materials[mat_idx as usize]
    }

    pub fn get_uv(&self, ray: &Ray, i: &Float3) -> Float2
    {
        let obj_idx = ray.obj_idx as usize;
        let obj_type = ray.obj_type;
        match obj_type
        {
            RayHittableObjectType::None =>
            {
                Float2::zero()
            },
            RayHittableObjectType::Sphere =>
            {
                self.spheres[obj_idx].get_uv(i)
            },
            RayHittableObjectType::Plane =>
            {
                self.planes[obj_idx].get_uv(i)
            },
            RayHittableObjectType::Cube =>
            {
                self.cubes[obj_idx].get_uv(i)
            },
            RayHittableObjectType::Quad =>
            {
                self.quads[obj_idx].get_uv(i)
            },
            RayHittableObjectType::Mesh =>
            {
                // should get uv from triangle
                self.meshes[obj_idx].get_uv(i)
            },
            RayHittableObjectType::Torus =>
            {
                self.tori[obj_idx].get_uv(i)
            },
            RayHittableObjectType::Bvh =>
            {
                self.bvhs[obj_idx].get_uv(i)
            }
        }
    }

    pub fn intersect_object(&self, ray: &mut Ray, obj_idx: usize, obj_type: RayHittableObjectType)
    {
        let obj_idx = obj_idx as usize;
        let obj_type = obj_type;
        match obj_type
        {
            RayHittableObjectType::None =>
            {

            },
            RayHittableObjectType::Sphere =>
            {
                self.spheres[obj_idx].intersect(ray);
            },
            RayHittableObjectType::Plane =>
            {
                self.planes[obj_idx].intersect(ray)
            },
            RayHittableObjectType::Cube =>
            {
                self.cubes[obj_idx].intersect(ray)
            },
            RayHittableObjectType::Quad =>
            {
                self.quads[obj_idx].intersect(ray)
            },
            RayHittableObjectType::Mesh =>
            {
                // should get uv from triangle
                self.meshes[obj_idx].intersect(ray)
            },
            RayHittableObjectType::Torus =>
            {
                self.tori[obj_idx].intersect(ray)
            },
            RayHittableObjectType::Bvh =>
            {
                self.bvhs[obj_idx].intersect(ray)
            }
        }
    }

    pub fn set_time(&mut self, time: f32)
    {
        self.animation_time = time;

        let m2_base = Mat4::rotate_x(PI / 4.0) * Mat4::rotate_z(PI / 4.0);
        let translation = Float3::from_xyz(1.8, 0.0, 2.5);
        let m2 = Mat4::translate(&translation) * Mat4::rotate_y(self.animation_time * 0.5) * m2_base;
        self.cubes[0].m = m2;
        self.cubes[0].inv_m = m2.inverted_no_scale();

        // bouncing sphere
        let modulo = (self.animation_time % 2.0) - 1.0;
        let modulo2 = modulo * modulo;
        let tm = 1.0 - modulo2;
        self.spheres[0].position = Float3::from_xyz(-1.8, -0.4 * tm, 2.0);
    }
}