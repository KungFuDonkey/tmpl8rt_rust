use crate::math::*;
use crate::ray::*;
use crate::objects::aabb::{AABB, intersect_aabb};
use crate::objects::bvh::*;
use crate::objects::kd_tree::*;
use crate::objects::grid::*;
use crate::timer::Timer;

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum MeshIntersectionSetting
{
    Raw,
    Bvh4,
    Bvh128,
    BvhSpatial4,
    BvhSpatial128,
    Grid16,
    Grid32,
    Grid64,
    KDTree8,
    KDTree16,
    KDTree24,
}

pub struct Mesh
{
    pub vertices: Vec<Float3>,
    pub triangles: Vec<[usize; 3]>,
    pub triangle_normals: Vec<Float3>,
    pub obj_idx: usize,
    pub mat_idx: usize,
    pub t: Mat4,
    pub inv_t: Mat4,
    pub bounds: AABB,
    pub bvh_4: BVH,
/*    pub bvh_128: BVH,
    pub bvh_spatial_4: BVH,
    pub bvh_spatial_128: BVH,
    pub grid_16: Grid,
    pub grid_32: Grid,
    pub grid_64: Grid,
    pub kd_tree_8: KDTree,
    pub kd_tree_16: KDTree,
    pub kd_tree_24: KDTree,
    pub mesh_intersection_setting: MeshIntersectionSetting*/
}

impl Mesh
{
    #[allow(dead_code)]
    pub fn from_tri_file(obj_idx: usize, mat_idx: usize, transform: Mat4, path: &std::path::Path) -> Self
    {
        let mut vertices = Vec::new();
        let mut triangles = Vec::new();
        let contents = std::fs::read_to_string(path).expect("Could not find file");
        let lines = contents.split("\n");

        for line in std::fs::read_to_string(path).expect("Could not find file").lines() {
            let floats = line.split(" ");

            let mut raw_floats: [f32; 9] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
            let mut index = 0;
            for float in floats
            {
                raw_floats[index] = float.parse::<f32>().unwrap();
                index += 1;
            }

            let vertices_index = vertices.len();
            vertices.push(Float3::from_xyz(raw_floats[0], raw_floats[1], raw_floats[2]));
            vertices.push(Float3::from_xyz(raw_floats[3], raw_floats[4], raw_floats[5]));
            vertices.push(Float3::from_xyz(raw_floats[6], raw_floats[7], raw_floats[8]));

            triangles.push([vertices_index, vertices_index + 1, vertices_index + 2])
        }

        let triangle_normals = Mesh::compute_normals(&triangles, &vertices, &transform);

        let bounds = Mesh::compute_bounds(&vertices);

        let mut timer = Timer::new();
        timer.reset();

        let bvh_4 = BVH::from_mesh(&triangles, &vertices, &bounds, 4);
        println!("bvh_4: {}", timer.elapsed());
        timer.reset();

        /*let bvh_128 = BVH::from_mesh(&triangles, &vertices, &bounds, 128);
        println!("bvh_128: {}", timer.elapsed());
        timer.reset();

        let bvh_spatial_4 = BVH::from_mesh_spatial(&triangles, &vertices, &bounds, 4);
        println!("bvh_spatial_4: {}", timer.elapsed());
        timer.reset();

        let bvh_spatial_128 = BVH::from_mesh_spatial(&triangles, &vertices, &bounds, 128);
        println!("bvh_spatial_128: {}", timer.elapsed());
        timer.reset();

        let kd_tree_8 = KDTree::from_mesh(&triangles, &vertices, &bounds, 8, 16, 128);
        println!("kd_tree_8 {}", timer.elapsed());
        timer.reset();

        let kd_tree_16 = KDTree::from_mesh(&triangles, &vertices, &bounds, 16, 16, 128);
        println!("kd_tree_16 {}", timer.elapsed());
        timer.reset();

        let kd_tree_24 = KDTree::from_mesh(&triangles, &vertices, &bounds, 24, 16, 128);
        println!("kd_tree_24 {}", timer.elapsed());
        timer.reset();

        let grid_16 = Grid::from_mesh(&triangles, &vertices, &bounds, 16, 16, 16);
        println!("grid_16 {}", timer.elapsed());
        timer.reset();

        let grid_32 = Grid::from_mesh(&triangles, &vertices, &bounds, 32, 32, 32);
        println!("grid_32 {}", timer.elapsed());
        timer.reset();

        let grid_64 = Grid::from_mesh(&triangles, &vertices, &bounds, 64, 64, 64);
        println!("grid_64 {}", timer.elapsed());
        timer.reset();*/

        Mesh
        {
            vertices,
            triangles,
            triangle_normals,
            obj_idx,
            mat_idx,
            t: transform,
            inv_t: transform.inverted(),
            bounds,
            //mesh_intersection_setting: MeshIntersectionSetting::BvhSpatial128,
            bvh_4,
            /*bvh_128,
            bvh_spatial_4,
            bvh_spatial_128,
            kd_tree_24,
            kd_tree_8,
            kd_tree_16,
            grid_16,
            grid_32,
            grid_64*/
        }
    }

    pub fn triangle(obj_idx: usize, mat_idx: usize, transform: Mat4) -> Self
    {
        let vertices = vec![Float3::from_xyz(0.0, 0.0, 0.0), Float3::from_xyz(1.0, 0.0, 0.0), Float3::from_xyz(0.0, 1.0, 0.3)];
        let triangles = vec![[0,2,1]];

        let triangle_normals = Mesh::compute_normals(&triangles, &vertices, &transform);
        let bounds = Mesh::compute_bounds(&vertices);

        let mut timer = Timer::new();
        timer.reset();

        let bvh_4 = BVH::from_mesh(&triangles, &vertices, &bounds, 4);
        println!("bvh_4: {}", timer.elapsed());
        timer.reset();
/*
        let bvh_128 = BVH::from_mesh(&triangles, &vertices, &bounds, 128);
        println!("bvh_128: {}", timer.elapsed());
        timer.reset();

        let bvh_spatial_4 = BVH::from_mesh_spatial(&triangles, &vertices, &bounds, 4);
        println!("bvh_spatial_4: {}", timer.elapsed());
        timer.reset();

        let bvh_spatial_128 = BVH::from_mesh_spatial(&triangles, &vertices, &bounds, 128);
        println!("bvh_spatial_128: {}", timer.elapsed());
        timer.reset();

        let kd_tree_8 = KDTree::from_mesh(&triangles, &vertices, &bounds, 8, 16, 128);
        println!("kd_tree_8 {}", timer.elapsed());
        timer.reset();

        let kd_tree_16 = KDTree::from_mesh(&triangles, &vertices, &bounds, 16, 16, 128);
        println!("kd_tree_16 {}", timer.elapsed());
        timer.reset();

        let kd_tree_24 = KDTree::from_mesh(&triangles, &vertices, &bounds, 24, 16, 128);
        println!("kd_tree_24 {}", timer.elapsed());
        timer.reset();

        let grid_16 = Grid::from_mesh(&triangles, &vertices, &bounds, 16, 16, 16);
        println!("grid_16 {}", timer.elapsed());
        timer.reset();

        let grid_32 = Grid::from_mesh(&triangles, &vertices, &bounds, 32, 32, 32);
        println!("grid_32 {}", timer.elapsed());
        timer.reset();

        let grid_64 = Grid::from_mesh(&triangles, &vertices, &bounds, 64, 64, 64);
        println!("grid_64 {}", timer.elapsed());
        timer.reset();*/

        Mesh
        {
            vertices,
            triangles,
            triangle_normals,
            obj_idx,
            mat_idx,
            t: transform,
            inv_t: transform.inverted(),
            bounds,
            //mesh_intersection_setting: MeshIntersectionSetting::BvhSpatial128,
            bvh_4,
            /*bvh_128,
            bvh_spatial_4,
            bvh_spatial_128,
            kd_tree_24,
            kd_tree_8,
            kd_tree_16,
            grid_16,
            grid_32,
            grid_64*/
        }
    }

    pub fn from_data(obj_idx: usize, mat_idx: usize, transform: Mat4, vertices: Vec<Float3>, triangles: Vec<[usize; 3]>) -> Self
    {
        let triangle_normals = Mesh::compute_normals(&triangles, &vertices, &transform);

        let bounds = Mesh::compute_bounds(&vertices);

        let mut timer = Timer::new();
        timer.reset();

        let bvh_4 = BVH::from_mesh(&triangles, &vertices, &bounds, 4);
        println!("bvh_4: {}", timer.elapsed());
        timer.reset();

        /*let bvh_128 = BVH::from_mesh(&triangles, &vertices, &bounds, 128);
        println!("bvh_128: {}", timer.elapsed());
        timer.reset();

        let bvh_spatial_4 = BVH::from_mesh_spatial(&triangles, &vertices, &bounds, 4);
        println!("bvh_spatial_4: {}", timer.elapsed());
        timer.reset();

        let bvh_spatial_128 = BVH::from_mesh_spatial(&triangles, &vertices, &bounds, 128);
        println!("bvh_spatial_128: {}", timer.elapsed());
        timer.reset();

        let kd_tree_8 = KDTree::from_mesh(&triangles, &vertices, &bounds, 8, 16, 128);
        println!("kd_tree_8 {}", timer.elapsed());
        timer.reset();

        let kd_tree_16 = KDTree::from_mesh(&triangles, &vertices, &bounds, 16, 16, 128);
        println!("kd_tree_16 {}", timer.elapsed());
        timer.reset();

        let kd_tree_24 = KDTree::from_mesh(&triangles, &vertices, &bounds, 24, 16, 128);
        println!("kd_tree_24 {}", timer.elapsed());
        timer.reset();

        let grid_16 = Grid::from_mesh(&triangles, &vertices, &bounds, 16, 16, 16);
        println!("grid_16 {}", timer.elapsed());
        timer.reset();

        let grid_32 = Grid::from_mesh(&triangles, &vertices, &bounds, 32, 32, 32);
        println!("grid_32 {}", timer.elapsed());
        timer.reset();

        let grid_64 = Grid::from_mesh(&triangles, &vertices, &bounds, 64, 64, 64);
        println!("grid_64 {}", timer.elapsed());
        timer.reset();*/

        Mesh
        {
            vertices,
            triangles,
            triangle_normals,
            obj_idx,
            mat_idx,
            t: transform,
            inv_t: transform.inverted(),
            bounds,
            //smesh_intersection_setting: MeshIntersectionSetting::BvhSpatial128,
            bvh_4,
            /*bvh_128,
            bvh_spatial_4,
            bvh_spatial_128,
            kd_tree_24,
            kd_tree_8,
            kd_tree_16,
            grid_16,
            grid_32,
            grid_64*/
        }
    }

    fn compute_normals(triangles: &Vec<[usize; 3]>, vertices: &Vec<Float3>, transform: &Mat4) -> Vec<Float3>
    {
        let mut triangle_normals: Vec<Float3> = Vec::with_capacity(triangles.len());
        for triangle in triangles
        {
            let v0 = vertices[triangle[0]];
            let v1 = vertices[triangle[1]];
            let v2 = vertices[triangle[2]];

            let v0v1 = v1 - v0;
            let v0v2 = v2 - v0;
            triangle_normals.push(normalize(&transform_vector(&cross(&v0v1, &v0v2), transform)));
        }
        return triangle_normals;
    }

    fn compute_bounds(vertices: &Vec<Float3>) -> AABB
    {
        // calculate bounding box
        let mut bounds = AABB::from_empty();
        for v in vertices
        {
            bounds.grow(v);
        }
        return bounds;
    }

    fn intersect_triangle(&self, ray: &mut Ray, vertex0: usize, vertex1: usize, vertex2: usize, triangle_index: usize) -> bool
    {
        ray.triangle_intersection_tests += 1;

        let v0 = self.vertices[vertex0];
        let v1 = self.vertices[vertex1];
        let v2 = self.vertices[vertex2];
        let v0v1 = v1 - v0;
        let v0v2 = v2 - v0;
        let pvec = cross(&ray.direction, &v0v2);
        let det = dot(&v0v1, &pvec);

        // cull backsides
        if det < EPSILON
        {
            return false;
        }

        let inv_det = 1.0 / det;

        let tvec = ray.origin - v0;
        let u = dot(&tvec, &pvec) * inv_det;
        if u < 0.0 || u > 1.0
        {
            return false;
        }

        let qvec = cross( &tvec, &v0v1);
        let v = dot(&ray.direction, &qvec) * inv_det;
        if v < 0.0 || u + v > 1.0
        {
            return false;
        }

        let t = dot(&v0v2, &qvec) * inv_det;
        if t > 0.0 && t < ray.t
        {
            ray.t = dot(&v0v2, &qvec) * inv_det;
            ray.sub_obj_idx = triangle_index;
            return true;
        }
        return false;
    }

    fn is_occluded_triangle(&self, ray: &Ray, vertex0: usize, vertex1: usize, vertex2: usize) -> bool
    {
        let v0 = self.vertices[vertex0];
        let v1 = self.vertices[vertex1];
        let v2 = self.vertices[vertex2];
        let v0v1 = v1 - v0;
        let v0v2 = v2 - v0;
        let pvec = cross(&ray.direction, &v0v2);
        let det = dot(&v0v1, &pvec);

        // do not cull backsides
        if det.abs() < EPSILON
        {
            return false;
        }

        let inv_det = 1.0 / det;

        let tvec = ray.origin - v0;
        let u = dot(&tvec, &pvec) * inv_det;
        if u < 0.0 || u > 1.0
        {
            return false;
        }

        let qvec = cross( &tvec, &v0v1);
        let v = dot(&ray.direction, &qvec) * inv_det;
        if v < 0.0 || u + v > 1.0
        {
            return false;
        }

        let t = dot(&v0v2, &qvec) * inv_det;
        return t > 0.0 && t < ray.t;
    }
}

impl RayHittableObject for Mesh
{
    fn intersect(&self, ray: &mut Ray) {
        let origin = transform_position(&ray.origin, &self.inv_t );
        let direction = transform_vector(&ray.direction, &self.inv_t );

        let mut ray_t = Ray::directed_distance(origin, direction, ray.t);
        ray_t.obj_idx = ray.obj_idx;

        // check for intersection with boundary of mesh
        let t = intersect_aabb(&mut ray_t, &self.bounds);
        if t == 1e30
        {
            return;
        }

        let mut intersected = self.bvh_4.intersect(&mut ray_t);
        /*match self.mesh_intersection_setting
        {
            MeshIntersectionSetting::Raw =>
            {
                for (triangle_index, triangle) in self.triangles.iter().enumerate()
                {
                    intersected = self.intersect_triangle(&mut ray_t, triangle[0], triangle[1], triangle[2], triangle_index) || intersected;
                }
            },
            MeshIntersectionSetting::Bvh4 =>
            {
                intersected = self.bvh_4.intersect(&mut ray_t);
            },
            MeshIntersectionSetting::Bvh128 =>
            {
                intersected = self.bvh_128.intersect(&mut ray_t);
            },
            MeshIntersectionSetting::BvhSpatial4 =>
            {
                intersected = self.bvh_spatial_4.intersect(&mut ray_t);
            },
            MeshIntersectionSetting::BvhSpatial128 =>
            {
                intersected = self.bvh_spatial_128.intersect(&mut ray_t);
            }
            MeshIntersectionSetting::KDTree24 =>
            {
                intersected = self.kd_tree_24.intersect(&mut ray_t, t, &self.bounds);
            }
            MeshIntersectionSetting::KDTree8 =>
            {
                intersected = self.kd_tree_8.intersect(&mut ray_t, t, &self.bounds);
            }
            MeshIntersectionSetting::KDTree16 =>
            {
                intersected = self.kd_tree_16.intersect(&mut ray_t, t, &self.bounds);
            }
            MeshIntersectionSetting::Grid16 =>
            {
                intersected = self.grid_16.intersect(&mut ray_t, t);
            }
            MeshIntersectionSetting::Grid32 =>
            {
                intersected = self.grid_32.intersect(&mut ray_t, t);
            }
            MeshIntersectionSetting::Grid64 =>
            {
                intersected = self.grid_64.intersect(&mut ray_t, t);
            }
        }*/


        ray.intersection_tests += ray_t.intersection_tests;
        ray.triangle_intersection_tests += ray_t.triangle_intersection_tests;
        if intersected
        {
            ray.t = ray_t.t;
            ray.obj_idx = self.obj_idx;
            ray.sub_obj_idx = ray_t.sub_obj_idx;
            ray.obj_type = RayHittableObjectType::Mesh;
        }
    }

    fn get_normal(&self, ray: &Ray, i: &Float3) -> Float3 {
        // needs to interpolate between normals
        self.triangle_normals[ray.sub_obj_idx]
    }

    fn get_uv(&self, i: &Float3) -> Float2 {
        Float2::zero()
    }

    fn is_occluded(&self, ray: &mut Ray) -> bool {
        let mut shadow = ray.clone();
        shadow.t = 1e34;
        shadow.obj_idx = usize::MAX;
        self.intersect(&mut shadow);
        ray.intersection_tests += shadow.intersection_tests;
        ray.triangle_intersection_tests += shadow.triangle_intersection_tests;
        return shadow.obj_idx != usize::MAX;
    }
}
