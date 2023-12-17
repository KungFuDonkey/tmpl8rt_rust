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
use crate::objects::triangle::*;
use crate::obj_loader::*;
use crate::opencl::{OpenCL, OpenCLBuffer};

pub struct Scene
{
    spheres: Vec<Sphere>,
    planes: Vec<Plane>,
    cubes: Vec<Cube>,
    tori: Vec<Torus>,
    quads: Vec<Quad>,
    meshes: Vec<Mesh>,
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
            reflective_material(linear_color_simple(Float3::from_xyz(1.0, 0.0, 0.0)), 0.8),
        ];

        let mut torus = Torus::new(0, 0, 0.8, 0.25);
        let translation = Float3::from_xyz(-0.25, 0.0, 2.0);
        torus.t = Mat4::translate(&translation) * Mat4::rotate_x(PI / 4.0);
        torus.inv_t = torus.t.inverted();

        let mut meshes: Vec<Mesh> = Vec::new();

        let transform = Mat4::translate( &Float3::from_xyz(0.0, 0.0, 1.0)) * Mat4::scale(0.5);
        meshes.push(Mesh::from_tri_file(0, 8, transform, std::path::Path::new("./assets/unity.tri")));

        let transform = Mat4::translate( &Float3::from_xyz(1.0, 0.0, 0.5)) * Mat4::scale(0.5);
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
                Plane::new(4, 1, Float3::from_xyz(0.0, 0.0, 1.0), 8.0, PlaneUVFunction::empty()),
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
            materials,
            animation_time: 0.0
        };

        scene.set_time(0.0);

        return scene;
    }

    pub fn change_intersection_setting(&mut self, mesh_setting: &MeshIntersectionSetting)
    {
        /*for mesh in &mut self.meshes
        {
            mesh.mesh_intersection_setting = *mesh_setting;
        }*/
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
    }

    pub fn is_occluded(&self, ray: &mut Ray) -> bool
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

    pub fn direct_lighting_soft(&self, point: &Float3, normal: &Float3, sample_size: usize, seed: &mut u32) -> (Float3, u32, u32)
    {
        let total_sample_points = sample_size;
        let sample_strength = 1.0 / (total_sample_points as f32);
        let mut lighting = Float3::zero();
        let mut traversal_steps: u32 = 0;
        let mut triangle_tests: u32 = 0;
        for _ in 0..sample_size
        {
            let quad_index = random_uint_s(seed) % (self.quads.len() as u32);
            let quad = &self.quads[quad_index as usize];
            let light_point = quad.random_point(seed);
            let ray_dir = light_point - *point;
            let ray_dir_n = normalize(&ray_dir);
            let origin = (*point) + (EPSILON * ray_dir_n);
            let mut ray = Ray::directed_distance(origin, ray_dir_n, length(&ray_dir) - (2.0 * EPSILON));

            let is_occluded = self.is_occluded(&mut ray);
            traversal_steps += ray.intersection_tests;
            triangle_tests += ray.triangle_intersection_tests;
            if is_occluded
            {
                continue;
            }

            let light_color = self.get_quad_light_color(&quad, &light_point);
            // take distance to the light?
            lighting += sample_strength * dot(&ray_dir_n, &normal) * light_color;
        }

        return (lighting, traversal_steps, triangle_tests);
    }

    pub fn direct_lighting_hard(&self, point: &Float3, normal: &Float3) -> (Float3, u32, u32)
    {
        let total_sample_points = self.quads.len() as u32;
        let light_strength = 1.0 / (total_sample_points as f32);
        let mut lighting = Float3::zero();
        let mut traversal_steps: u32 = 0;
        let mut triangle_tests: u32 = 0;
        for quad in &self.quads
        {
            let light_point = quad.center_point();
            let ray_dir = light_point - *point;
            let ray_dir_n = normalize(&ray_dir);
            let origin = (*point) + (EPSILON * ray_dir_n);
            let mut ray = Ray::directed_distance(origin, ray_dir_n, length(&ray_dir) - (2.0 * EPSILON));

            let is_occluded = self.is_occluded(&mut ray);
            traversal_steps += ray.intersection_tests;
            triangle_tests += ray.triangle_intersection_tests;
            if is_occluded
            {
                continue;
            }

            let light_color = self.get_quad_light_color(&quad, &light_point);

            // take distance to the light?
            lighting += light_strength * dot(&ray_dir_n, &normal) * light_color;
        }

        return (lighting, traversal_steps, triangle_tests);
    }


    pub fn get_normal(&self, ray: &Ray, i: &Float3, wo: &Float3) -> Float3
    {
        let obj_idx = ray.obj_idx;
        if obj_idx == usize::MAX
        {
            panic!("ERROR: obj_idx not set or no object was hit in get_normal()");
        }

        let obj_idx = obj_idx as usize;
        let obj_type = ray.obj_type;

        let normal = match obj_type
        {
            RayHittableObjectType::None =>
            {
                panic!("ERROR: tried to get normal of non object");
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
                self.meshes[obj_idx].get_uv(i)
            },
            RayHittableObjectType::Torus =>
            {
                self.tori[obj_idx].get_uv(i)
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
                {},
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

pub struct GPUScene
{
    pub sphere_positions: OpenCLBuffer<Float3>,
    pub sphere_radi: OpenCLBuffer<f32>,
    pub sphere_materials: OpenCLBuffer<u64>,
    pub plane_normals: OpenCLBuffer<Float3>,
    pub plane_distances: OpenCLBuffer<f32>,
    pub plane_materials: OpenCLBuffer<u64>,
    pub quad_sizes: OpenCLBuffer<f32>,
    pub quad_inv_transforms: OpenCLBuffer<Mat4>,
    pub quad_materials: OpenCLBuffer<u64>,
    pub mesh_offsets: OpenCLBuffer<u32>,
    pub mesh_triangle_offsets: OpenCLBuffer<u32>,
    pub mesh_inv_transforms: OpenCLBuffer<Mat4>,
    pub mesh_min_bounds: OpenCLBuffer<Float3>,
    pub mesh_max_bounds: OpenCLBuffer<Float3>,
    pub mesh_tri_counts: OpenCLBuffer<u32>,
    pub mesh_left_firsts: OpenCLBuffer<u32>,
    pub mesh_triangles: OpenCLBuffer<Triangle>,
    pub mesh_triangle_normals: OpenCLBuffer<Float3>,
    pub mesh_materials: OpenCLBuffer<u64>
}

impl GPUScene
{
    pub fn new(cl: &OpenCL) -> Self
    {
        let scene = Scene::new();

        return GPUScene::from_scene(cl, &scene);
    }

    pub fn from_scene(cl: &OpenCL, scene: &Scene) -> Self
    {
        let mut sphere_positions: Vec<Float3> = Vec::with_capacity(scene.spheres.len());
        let mut sphere_radi: Vec<f32> = Vec::with_capacity(scene.spheres.len());
        let mut sphere_materials: Vec<u64> = Vec::with_capacity(scene.spheres.len());
        for sphere in &scene.spheres
        {
            sphere_positions.push(sphere.position);
            sphere_radi.push(sphere.radius);
            sphere_materials.push(sphere.mat_idx as u64);
        }
        // todo convert materials
        sphere_materials[0] = (255 << 16);
        sphere_materials[1] = 255;

        let mut plane_normals: Vec<Float3> = Vec::with_capacity(scene.planes.len());
        let mut plane_distances: Vec<f32> = Vec::with_capacity(scene.planes.len());
        let mut plane_materials: Vec<u64> = Vec::with_capacity(scene.planes.len());
        for plane in &scene.planes
        {
            plane_normals.push(plane.normal);
            plane_distances.push(plane.distance);
            plane_materials.push(plane.mat_idx as u64);
        }

        plane_materials[0] = (255 << 16) + (128 << 8) + 255;
        plane_materials[1] = (255 << 16) + (0   << 8) + 255;
        plane_materials[2] = (0   << 16) + (255 << 8) + 255;
        plane_materials[3] = (255 << 16) + (255 << 8) + 0  ;
        plane_materials[4] = (0   << 16) + (255 << 8) + 0  ;
        plane_materials[5] = (128 << 16) + (128 << 8) + 128;

        let mut quad_sizes: Vec<f32> = Vec::with_capacity(scene.quads.len());
        let mut quad_inv_transforms: Vec<Mat4> = Vec::with_capacity(scene.quads.len());
        let mut quad_materials: Vec<u64> = Vec::with_capacity(scene.quads.len());
        for quad in &scene.quads
        {
            quad_sizes.push(quad.size);
            quad_inv_transforms.push(quad.inv_t);
            quad_materials.push(quad.mat_idx as u64);
        }

        quad_materials[0] = (255 << 16) + (255 << 8) + 255 + (1 << 61);
        quad_materials[1] = (255 << 16) + (255 << 8) + 255 + (1 << 61);
        quad_materials[2] = (255 << 16) + (255 << 8) + 255 + (1 << 61);
        quad_materials[3] = (255 << 16) + (255 << 8) + 255 + (1 << 61);

        let mut mesh_offsets: Vec<u32> = Vec::with_capacity(scene.meshes.len());
        let mut mesh_triangle_offsets: Vec<u32> = Vec::with_capacity(scene.meshes.len());
        let mut mesh_inv_transforms: Vec<Mat4> = Vec::with_capacity(scene.meshes.len());
        let mut mesh_min_bounds: Vec<Float3> = Vec::with_capacity(scene.meshes.len());
        let mut mesh_max_bounds: Vec<Float3> = Vec::with_capacity(scene.meshes.len());
        let mut mesh_tri_counts: Vec<u32> = Vec::with_capacity(scene.meshes.len());
        let mut mesh_left_firsts: Vec<u32> = Vec::with_capacity(scene.meshes.len());
        let mut mesh_triangles: Vec<Triangle> = Vec::with_capacity(scene.meshes.len());
        let mut mesh_triangle_normals: Vec<Float3> = Vec::with_capacity(scene.meshes.len());
        let mut mesh_materials: Vec<u64> = Vec::with_capacity(scene.meshes.len());

        let mut mesh_offset = 0;
        let mut triangle_offset = 0;
        for mesh in &scene.meshes
        {
            mesh_offsets.push(mesh_offset);
            mesh_offset += mesh.bvh_4.bvh_nodes.len() as u32;
            mesh_inv_transforms.push(mesh.inv_t);
            for bvh_node in &mesh.bvh_4.bvh_nodes
            {
                mesh_min_bounds.push(bvh_node.bounds.min_bound);
                mesh_max_bounds.push(bvh_node.bounds.max_bound);
                mesh_tri_counts.push(bvh_node.tri_count as u32);
                mesh_left_firsts.push(bvh_node.left_first as u32);
            }
            for id in &mesh.bvh_4.triangle_idx
            {
                mesh_triangles.push(mesh.bvh_4.triangles[*id].internal_triangle);
            }
            for normal in &mesh.triangle_normals
            {
                mesh_triangle_normals.push(*normal);
            }
            mesh_triangle_offsets.push(triangle_offset);
            triangle_offset += mesh.bvh_4.triangle_idx.len() as u32;
            mesh_materials.push(mesh.mat_idx as u64);
        }

        mesh_materials[0] = (255 << 16);
        mesh_materials[1] = (255 << 8);

        let sphere_positions = OpenCLBuffer::read_write(cl, sphere_positions);
        let sphere_radi = OpenCLBuffer::read_write(cl, sphere_radi);
        let sphere_materials = OpenCLBuffer::read_write(cl, sphere_materials);
        let plane_normals = OpenCLBuffer::read_write(cl, plane_normals);
        let plane_distances = OpenCLBuffer::read_write(cl, plane_distances);
        let plane_materials = OpenCLBuffer::read_write(cl, plane_materials);
        let quad_sizes = OpenCLBuffer::read_write(cl, quad_sizes);
        let quad_inv_transforms = OpenCLBuffer::read_write(cl, quad_inv_transforms);
        let quad_materials = OpenCLBuffer::read_write(cl, quad_materials);
        let mesh_offsets = OpenCLBuffer::read_write(cl, mesh_offsets);
        let mesh_triangle_offsets = OpenCLBuffer::read_write(cl, mesh_triangle_offsets);
        let mesh_inv_transforms = OpenCLBuffer::read_write(cl, mesh_inv_transforms);
        let mesh_min_bounds =  OpenCLBuffer::read_write(cl, mesh_min_bounds);
        let mesh_max_bounds = OpenCLBuffer::read_write(cl, mesh_max_bounds);
        let mesh_tri_counts = OpenCLBuffer::read_write(cl, mesh_tri_counts);
        let mesh_left_firsts = OpenCLBuffer::read_write(cl, mesh_left_firsts);
        let mesh_triangles = OpenCLBuffer::read_write(cl, mesh_triangles);
        let mesh_triangle_normals = OpenCLBuffer::read_write(cl, mesh_triangle_normals);
        let mesh_materials = OpenCLBuffer::read_write(cl, mesh_materials);

        sphere_positions.copy_to_device(cl);
        sphere_radi.copy_to_device(cl);
        sphere_materials.copy_to_device(cl);
        plane_normals.copy_to_device(cl);
        plane_distances.copy_to_device(cl);
        plane_materials.copy_to_device(cl);
        quad_sizes.copy_to_device(cl);
        quad_inv_transforms.copy_to_device(cl);
        quad_materials.copy_to_device(cl);
        mesh_offsets.copy_to_device(cl);
        mesh_triangle_offsets.copy_to_device(cl);
        mesh_inv_transforms.copy_to_device(cl);
        mesh_min_bounds.copy_to_device(cl);
        mesh_max_bounds.copy_to_device(cl);
        mesh_tri_counts.copy_to_device(cl);
        mesh_left_firsts.copy_to_device(cl);
        mesh_triangles.copy_to_device(cl);
        mesh_triangle_normals.copy_to_device(cl);
        mesh_materials.copy_to_device(cl);

        GPUScene
        {
            sphere_positions,
            sphere_radi,
            sphere_materials,
            plane_normals,
            plane_distances,
            plane_materials,
            quad_sizes,
            quad_inv_transforms,
            quad_materials,
            mesh_inv_transforms,
            mesh_offsets,
            mesh_triangle_offsets,
            mesh_min_bounds,
            mesh_max_bounds,
            mesh_tri_counts,
            mesh_left_firsts,
            mesh_triangles,
            mesh_triangle_normals,
            mesh_materials
        }
    }
}