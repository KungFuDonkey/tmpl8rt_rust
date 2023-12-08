use std::f32::consts::PI;
use crate::math::*;
use crate::material::*;
use crate::obj_loader::load_obj;
use crate::ray::*;
use crate::objects::sphere::*;
use crate::objects::plane::*;
use crate::objects::cube::*;
use crate::objects::torus::*;
use crate::objects::quad::*;
use crate::objects::mesh::*;
use crate::objects::bvh::*;
use crate::objects::grid::*;
use crate::objects::kd_tree::*;

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum MeshIntersectionSetting
{
    Raw,
    Bvh4,
    Bvh128,
    BvhSpatial4,
    BvhSpatial128,
    Grid,
    KDTree
}


pub struct Scene
{
    spheres: Vec<Sphere>,
    planes: Vec<Plane>,
    cubes: Vec<Cube>,
    tori: Vec<Torus>,
    quads: Vec<Quad>,
    meshes: Vec<Mesh>,
    bvh_4: Vec<BVH>,
    bvh_128: Vec<BVH>,
    bvh_spatial_4: Vec<BVH>,
    bvh_spatial_128: Vec<BVH>,
    grids: Vec<Grid>,
    kd_trees: Vec<KDTree>,
    materials: Vec<Material>,
    animation_time: f32,
}

impl Scene
{
    pub fn new() -> Self
    {
        let mut materials = vec![
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
        ];

        let mut torus = Torus::new(0, 0, 0.8, 0.25);
        let translation = Float3::from_xyz(-0.25, 0.0, 2.0);
        torus.t = Mat4::translate(&translation) * Mat4::rotate_x(PI / 4.0);
        torus.inv_t = torus.t.inverted();

        let mut meshes: Vec<Mesh> = Vec::new();
        //let transform = Mat4::translate( &Float3::from_xyz(0.0, 0.0, 2.5)) * Mat4::scale(0.5);
        //meshes.push(Mesh::triangle(0, 6, transform));

        //let transform = Mat4::translate( &Float3::from_xyz(0.0, 0.0, 1.0)) * Mat4::scale(0.5);
        //meshes.push(Mesh::from_tri_file(1, 8, transform, std::path::Path::new("./assets/unity.tri")));

        let transform = Mat4::translate( &Float3::from_xyz(0.0, 0.0, 0.5)) * Mat4::scale(0.5);
        let (msh, mts) = load_obj(&std::path::Path::new("./assets/suzanne.obj"), meshes.len(), materials.len(), &transform);

        for mesh in msh
        {
            meshes.push(mesh);
        }

        for material in mts
        {
            materials.push(material);
        }

        let mut scene = Scene{
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
                //torus
            ],
            quads: vec![
                Quad::new(0, 6, 0.5, &Mat4::translate(&Float3::from_xyz(-1.0, 1.5, -1.0))),
                Quad::new(1, 7, 0.5, &Mat4::translate(&Float3::from_xyz(1.0, 1.5, -1.0))),
                Quad::new(2, 8, 0.5, &Mat4::translate(&Float3::from_xyz(1.0, 1.5, 1.0))),
                Quad::new(3, 0, 0.5, &Mat4::translate(&Float3::from_xyz(-1.0, 1.5, 1.0))),
            ],
            meshes,
            bvh_4: Vec::new(),
            bvh_spatial_4: Vec::new(),
            bvh_128: Vec::new(),
            bvh_spatial_128: Vec::new(),
            grids: Vec::new(),
            kd_trees: Vec::new(),
            materials,
            animation_time: 0.0
        };

        scene.construct_acceleration_structures();

        return scene;
    }

    fn construct_acceleration_structures(&mut self)
    {
        for mesh in &self.meshes
        {
            self.grids.push(Grid::from_mesh(mesh, 64, 64, 64));
            self.kd_trees.push(KDTree::from_mesh(mesh, 16, 2, 128));
            self.bvh_4.push(BVH::from_mesh(mesh, 4));
            self.bvh_128.push(BVH::from_mesh(mesh, 128));
            self.bvh_spatial_4.push(BVH::from_mesh_spatial(mesh, 4));
            self.bvh_spatial_128.push(BVH::from_mesh_spatial(mesh, 128));
        }
    }

    pub fn intersect_scene(&self, ray: &mut Ray, mesh_setting: &MeshIntersectionSetting)
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

        match mesh_setting {
            MeshIntersectionSetting::Raw =>
            {
                for mesh in &self.meshes
                {
                    mesh.intersect(ray);
                }
            },
            MeshIntersectionSetting::Bvh4 =>
            {
                for bvh in &self.bvh_4
                {
                    bvh.intersect(ray);
                }
            },
            MeshIntersectionSetting::Bvh128 =>
            {
                for bvh in &self.bvh_128
                {
                    bvh.intersect(ray);
                }
            },
            MeshIntersectionSetting::BvhSpatial4 =>
            {
                for bvh in &self.bvh_spatial_4
                {
                    bvh.intersect(ray);
                }
            },
            MeshIntersectionSetting::BvhSpatial128 =>
            {
                for bvh in &self.bvh_spatial_128
                {
                    bvh.intersect(ray);
                }
            },
            MeshIntersectionSetting::Grid =>
            {
                for grid in &self.grids
                {
                    grid.intersect(ray);
                }
            },
            MeshIntersectionSetting::KDTree =>
            {
                for kd_tree in &self.kd_trees
                {
                    kd_tree.intersect(ray);
                }
            }
        }
    }

    pub fn is_occluded(&self, ray: &Ray, mesh_setting: &MeshIntersectionSetting) -> bool
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

        match mesh_setting {
            MeshIntersectionSetting::Raw => {
                for mesh in &self.meshes
                {
                    if mesh.is_occluded(ray)
                    {
                        return true;
                    }
                }
            }
            MeshIntersectionSetting::Bvh4 => {
                for bvh in &self.bvh_4
                {
                    if bvh.is_occluded(ray)
                    {
                        return true;
                    }
                }
            }
            MeshIntersectionSetting::BvhSpatial4 => {
                for bvh in &self.bvh_spatial_4
                {
                    if bvh.is_occluded(ray)
                    {
                        return true;
                    }
                }
            }
            MeshIntersectionSetting::Bvh128 => {
                for bvh in &self.bvh_128
                {
                    if bvh.is_occluded(ray)
                    {
                        return true;
                    }
                }
            }
            MeshIntersectionSetting::BvhSpatial128 => {
                for bvh in &self.bvh_spatial_128
                {
                    if bvh.is_occluded(ray)
                    {
                        return true;
                    }
                }
            }
            MeshIntersectionSetting::Grid => {
                for grid in &self.grids
                {
                    if grid.is_occluded(ray)
                    {
                        return true;
                    }
                }
            }
            MeshIntersectionSetting::KDTree =>
            {
                for kd_tree in &self.kd_trees
                {
                    if kd_tree.is_occluded(ray)
                    {
                        return true;
                    }
                }
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

    pub fn direct_lighting_soft(&self, point: &Float3, normal: &Float3, sample_size: usize, seed: &mut u32, mesh_setting: &MeshIntersectionSetting) -> Float3
    {
        let total_sample_points = sample_size;
        let sample_strength = 1.0 / (total_sample_points as f32);
        let mut lighting = Float3::zero();

        for _ in 0..sample_size
        {
            let quad_index = random_uint_s(seed) % (self.quads.len() as u32);
            let quad = &self.quads[quad_index as usize];
            let light_point = quad.random_point(seed);
            let ray_dir = light_point - *point;
            let ray_dir_n = normalize(&ray_dir);
            let origin = (*point) + (EPSILON * ray_dir_n);
            let ray = Ray::directed_distance(origin, ray_dir_n, length(&ray_dir) - (2.0 * EPSILON));

            if self.is_occluded(&ray, mesh_setting)
            {
                continue;
            }

            let light_color = self.get_quad_light_color(&quad, &light_point);
            // take distance to the light?
            lighting += sample_strength * dot(&ray_dir_n, &normal) * light_color;
        }

        return lighting;
    }

    pub fn direct_lighting_hard(&self, point: &Float3, normal: &Float3, mesh_setting: &MeshIntersectionSetting) -> Float3
    {
        let total_sample_points = self.quads.len() as u32;
        let light_strength = 1.0 / (total_sample_points as f32);
        let mut lighting = Float3::zero();
        for quad in &self.quads
        {
            let light_point = quad.center_point();
            let ray_dir = light_point - *point;
            let ray_dir_n = normalize(&ray_dir);
            let origin = (*point) + (EPSILON * ray_dir_n);
            let ray = Ray::directed_distance(origin, ray_dir_n, length(&ray_dir) - (2.0 * EPSILON));

            if self.is_occluded(&ray, mesh_setting)
            {
                continue;
            }

            let light_color = self.get_quad_light_color(&quad, &light_point);

            // take distance to the light?
            lighting += light_strength * dot(&ray_dir_n, &normal) * light_color;
        }

        return lighting;
    }


    pub fn get_normal(&self, ray: &Ray, i: &Float3, wo: &Float3, mesh_intersection_setting: &MeshIntersectionSetting) -> Float3
    {
        let obj_idx = ray.obj_idx;
        if obj_idx == usize::MAX
        {
            println!("ERROR: obj_idx not set or no object was hit in get_normal()");
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
                self.spheres[obj_idx].get_normal(ray, i)
            },
            RayHittableObjectType::Plane =>
            {
                self.planes[obj_idx].get_normal(ray, i)
            },
            RayHittableObjectType::Cube =>
            {
                self.cubes[obj_idx].get_normal(ray, i)
            },
            RayHittableObjectType::Torus =>
            {
                self.tori[obj_idx].get_normal(ray, i)
            },
            RayHittableObjectType::Mesh =>
            {
                self.meshes[obj_idx].get_normal(ray, i)
            },
            RayHittableObjectType::Quad =>
            {
                self.quads[obj_idx].get_normal(ray, i)
            }
            RayHittableObjectType::Bvh =>
            {
                match mesh_intersection_setting
                {
                    MeshIntersectionSetting::Bvh4 => self.bvh_4[obj_idx].get_normal(ray, i),
                    MeshIntersectionSetting::Bvh128 => self.bvh_128[obj_idx].get_normal(ray, i),
                    MeshIntersectionSetting::BvhSpatial4 => self.bvh_spatial_4[obj_idx].get_normal(ray, i),
                    MeshIntersectionSetting::BvhSpatial128 => self.bvh_spatial_128[obj_idx].get_normal(ray, i),
                    _ => Float3::zero()
                }
            }
            RayHittableObjectType::Grid =>
            {
                self.grids[obj_idx].get_normal(ray, i)
            }
            RayHittableObjectType::KDTree =>
            {
                self.kd_trees[obj_idx].get_normal(ray, i)
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
                usize::MAX
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
                self.meshes[obj_idx as usize].mat_idx
            },
            RayHittableObjectType::Grid =>
            {
                self.grids[obj_idx as usize].mat_idx
            },
            RayHittableObjectType::KDTree =>
            {
                self.kd_trees[obj_idx].mat_idx
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
                self.meshes[obj_idx].get_uv(i)
            },
            RayHittableObjectType::Grid =>
            {
                self.grids[obj_idx].get_uv(i)
            },
            RayHittableObjectType::KDTree =>
            {
                self.kd_trees[obj_idx].get_uv(i)
            }
        }
    }

    pub fn intersect_object(&self, ray: &mut Ray, obj_idx: usize, obj_type: RayHittableObjectType, mesh_intersection_setting: &MeshIntersectionSetting)
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
                match mesh_intersection_setting
                {
                    MeshIntersectionSetting::Bvh4 => self.bvh_4[obj_idx].intersect(ray),
                    MeshIntersectionSetting::Bvh128 => self.bvh_128[obj_idx].intersect(ray),
                    MeshIntersectionSetting::BvhSpatial4 => self.bvh_spatial_4[obj_idx].intersect(ray),
                    MeshIntersectionSetting::BvhSpatial128 => self.bvh_spatial_128[obj_idx].intersect(ray),
                    _ => panic!("cannot reach this")
                }
            },
            RayHittableObjectType::Grid =>
            {
                self.grids[obj_idx].intersect(ray)
            },
            RayHittableObjectType::KDTree =>
            {
                self.kd_trees[obj_idx].intersect(ray)
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