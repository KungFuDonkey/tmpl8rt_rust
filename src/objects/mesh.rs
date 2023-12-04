use crate::math::*;
use crate::ray::*;

pub struct Mesh
{
    pub vertices: Vec<Float3>,
    pub triangles: Vec<[usize; 3]>,
    pub triangle_normals: Vec<Float3>,
    pub obj_idx: i32,
    pub mat_idx: i32,
    pub t: Mat4,
    pub inv_t: Mat4
}

impl Mesh
{
    pub fn from_tri_file(obj_idx: i32, mat_idx: i32, transform: Mat4, path: &std::path::Path) -> Self
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

        Mesh
        {
            vertices,
            triangles,
            triangle_normals,
            obj_idx,
            mat_idx,
            t: transform,
            inv_t: transform.inverted()
        }
    }

    pub fn triangle(obj_idx: i32, mat_idx: i32, transform: Mat4) -> Self
    {
        let vertices = vec![Float3::from_xyz(0.0, 0.0, 0.0), Float3::from_xyz(1.0, 0.0, 0.0), Float3::from_xyz(0.0, 1.0, 0.3)];
        let triangles = vec![[0,2,1]];

        let triangle_normals = Mesh::compute_normals(&triangles, &vertices, &transform);

        Mesh
        {
            vertices,
            triangles,
            triangle_normals,
            obj_idx,
            mat_idx,
            t: transform,
            inv_t: transform.inverted()
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

    fn intersect_triangle(&self, ray: &mut Ray, vertex0: usize, vertex1: usize, vertex2: usize, triangle_index: usize) -> bool
    {
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

        for (triangle_index, triangle) in self.triangles.iter().enumerate()
        {
            if self.intersect_triangle(&mut ray_t, triangle[0], triangle[1], triangle[2], triangle_index)
            {
                ray.obj_idx = self.obj_idx;
                ray.sub_obj_idx = ray_t.sub_obj_idx;
                ray.t = ray_t.t;
                ray.obj_type = RayHittableObjectType::Mesh;
            }
        }
    }

    fn get_normal(&self, ray: &Ray, i: &Float3) -> Float3 {
        // needs to interpolate between normals
        self.triangle_normals[ray.sub_obj_idx]
    }

    fn get_uv(&self, i: &Float3) -> Float2 {
        Float2::zero()
    }

    fn is_occluded(&self, ray: &Ray) -> bool {
        for triangle in &self.triangles
        {
            if self.is_occluded_triangle(ray, triangle[0], triangle[1], triangle[2])
            {
                return true;
            }
        }
        return false;
    }
}
