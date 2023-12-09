use crate::math::*;
use crate::objects::aabb::*;
use crate::objects::triangle::*;
use crate::objects::mesh::*;
use crate::ray::{Ray, RayHittableObject, RayHittableObjectType};

struct KDTriangle
{
    pub internal_triangle: Triangle
}

struct KDNode
{
    pub plane_value: f32,
    pub left_first: usize,
    pub triangle_count: usize,
    pub is_leaf: bool
}

impl KDNode
{
    pub fn is_leaf(&self) -> bool
    {
        self.is_leaf
    }
}

pub struct KDTree
{
    triangles: Vec<KDTriangle>,
    kd_nodes: Vec<KDNode>,
    pub triangle_ids: Vec<usize>,
    pub t: Mat4,
    pub inv_t: Mat4,
    pub normals: Vec<Float3>,
    pub obj_idx: usize,
    pub mat_idx: usize,
    pub max_depth: u32,
    pub max_triangles_in_leaf: usize,
    pub scan_samples: u32,
    pub bounds: AABB
}

impl KDTree
{
    pub fn from_mesh(mesh: &Mesh, max_depth: u32, max_triangles_in_leaf: usize, scan_samples: u32) -> KDTree
    {
        let triangle_count = mesh.triangles.len();
        let mut triangles: Vec<KDTriangle> = Vec::with_capacity(triangle_count);
        let triangle_ids: Vec<usize> = Vec::with_capacity(triangle_count * 2);

        let mut start_triangle_ids: Vec<usize> = Vec::with_capacity(triangle_count);

        let mut tri_idx = 0;
        for triangle in &mesh.triangles
        {
            let vertex0 = mesh.vertices[triangle[0]];
            let vertex1 = mesh.vertices[triangle[1]];
            let vertex2 = mesh.vertices[triangle[2]];

            triangles.push(KDTriangle{
                internal_triangle: Triangle {
                    vertex0,
                    vertex1,
                    vertex2,
                    tri_idx
                }
            });

            start_triangle_ids.push(tri_idx as usize);
            tri_idx += 1;
        }

        let mut kd_nodes: Vec<KDNode> = Vec::with_capacity(triangle_count * 4);

        let root = KDNode
        {
            left_first: 0,
            triangle_count: 0,
            plane_value: 0.0,
            is_leaf: false
        };

        kd_nodes.push(root);

        let mut kd_tree = KDTree
        {
            triangles,
            triangle_ids,
            kd_nodes,
            normals: mesh.triangle_normals.to_vec(),
            t: mesh.t,
            inv_t: mesh.inv_t,
            obj_idx: mesh.obj_idx,
            mat_idx: mesh.mat_idx,
            max_depth,
            max_triangles_in_leaf,
            scan_samples,
            bounds: mesh.bounds
        };

        kd_tree.build_tree(&start_triangle_ids);

        return kd_tree;
    }

    fn build_tree(&mut self, start_ids: &Vec<usize>)
    {
        let bounds = self.bounds;
        self.build(&bounds, 0, start_ids, 0);
    }

    fn find_split_plane(&self, bounds: &AABB, triangle_ids: &Vec<usize>, axis: usize) -> f32
    {
        let min_value = bounds.min_bound.get_axis(axis);
        let max_value = bounds.max_bound.get_axis(axis);
        let range = max_value - min_value;
        let step_size = range / (self.scan_samples as f32);
        let mut best_value: f32 = 1e30;
        let mut best_score: usize = 0;

        let mut current_value = min_value;
        for i in 0..self.scan_samples
        {
            current_value += step_size;

            let (left_indices, right_indices) = self.partition_ids(current_value, axis, triangle_ids);
            let same_indices = left_indices.len() + right_indices.len() - triangle_ids.len();
            let split_score = (left_indices.len() - same_indices) * (right_indices.len() - same_indices);
            if split_score > best_score
            {
                best_score = split_score;
                best_value = current_value;
            }
        }

        // could not find any split with a score > 0, can happen when partitioning on y, but the object is flat
        if best_value == 1e30
        {
            // take center
            best_value = (max_value + min_value) / 2.0;
        }

        return best_value;
    }

    fn build(&mut self, bounds: &AABB, current_node_id: usize, current_triangle_ids: &Vec<usize>, depth: u32)
    {
        if depth == self.max_depth || current_triangle_ids.len() < self.max_triangles_in_leaf
        {
            self.build_leaf_node(current_node_id, current_triangle_ids);
            return;
        }

        let axis = (depth % 3) as usize;
        let best_plane = self.find_split_plane(bounds, current_triangle_ids, axis);

        let (left_triangles, right_triangles) = self.partition_ids(best_plane, axis, current_triangle_ids);

        self.kd_nodes[current_node_id].left_first = self.kd_nodes.len();
        self.kd_nodes[current_node_id].plane_value = best_plane;
        self.kd_nodes.push(KDNode
        {
            plane_value: 0.0,
            left_first: 0,
            triangle_count: 0,
            is_leaf: false
        });

        self.kd_nodes.push(KDNode
        {
            plane_value: 0.0,
            left_first: 0,
            triangle_count: 0,
            is_leaf: false
        });

        let mut left_bounds = *bounds;
        left_bounds.max_bound.set_axis(axis, best_plane);
        let mut right_bounds = *bounds;
        right_bounds.min_bound.set_axis(axis, best_plane);

        self.build(&left_bounds, self.kd_nodes[current_node_id].left_first, &left_triangles, depth + 1);
        self.build(&right_bounds, self.kd_nodes[current_node_id].left_first + 1, &right_triangles, depth + 1);
    }

    fn build_leaf_node(&mut self, current_node_id: usize, current_triangle_ids: &Vec<usize>)
    {
        self.kd_nodes[current_node_id].triangle_count = current_triangle_ids.len();
        self.kd_nodes[current_node_id].left_first = self.triangle_ids.len();
        self.kd_nodes[current_node_id].is_leaf = true;
        self.triangle_ids.extend(current_triangle_ids);
    }

    fn partition_ids(&self, plane_value: f32, axis: usize, triangle_ids: &Vec<usize>) -> (Vec<usize>, Vec<usize>)
    {
        let mut left_triangles: Vec<usize> = Vec::new();
        let mut right_triangles: Vec<usize> = Vec::new();

        for triangle_id in triangle_ids
        {
            let id = *triangle_id;
            let triangle = &self.triangles[id].internal_triangle;

            let v0 = triangle.vertex0.get_axis(axis);
            let v1 = triangle.vertex1.get_axis(axis);
            let v2 = triangle.vertex2.get_axis(axis);
            let min_value = v0.min(v1.min(v2));
            let max_value = v0.max(v1.max(v2));

            if max_value < plane_value - EPSILON
            {
                left_triangles.push(id);
            }
            else if min_value > plane_value + EPSILON
            {
                right_triangles.push(id);
            }
            else
            {
                left_triangles.push(id);
                right_triangles.push(id);
            }
        }

        return (left_triangles, right_triangles)
    }

    fn intersect_kd_node(&self, ray: &mut Ray, current_t: f32, node_id: usize, depth: u32, bounds: &AABB) -> bool
    {
        if current_t >= ray.t
        {
            return false;
        }

        ray.intersection_tests += 1;

        let current_node = &self.kd_nodes[node_id];

        // if node is a leaf check if it intersects with the geometry
        if current_node.is_leaf()
        {
            if current_node.triangle_count == 0
            {
                return false;
            }
            let mut intersected = false;
            for i in &self.triangle_ids[current_node.left_first..(current_node.left_first + current_node.triangle_count)]
            {
                let triangle = &self.triangles[*i].internal_triangle;
                intersected = intersect_triangle(triangle, ray) || intersected;
            }

            return intersected;
        }

        let axis = (depth % 3) as usize;

        // check if there is a plane intersection WITHIN the bounds of the current node (and positive)
        let current_origin_point = ray.origin + ray.direction * current_t;
        let ray_direction_value = ray.direction.get_axis(axis);
        let ray_origin_value = current_origin_point.get_axis(axis);
        let distance = current_node.plane_value - ray_origin_value;
        let t = distance / ray_direction_value;
        let mut plane_intersection = false;
        if t > 0.0
        {
            let intersection_point = current_origin_point + ray.direction * t;
            let a1 = (axis + 1) % 3;
            let a2 = (axis + 2) % 3;

            let i_a1 = intersection_point.get_axis(a1);
            let i_a2 = intersection_point.get_axis(a2);
            let min_a1 = bounds.min_bound.get_axis(a1);
            let min_a2 = bounds.min_bound.get_axis(a2);
            let max_a1 = bounds.max_bound.get_axis(a1);
            let max_a2 = bounds.max_bound.get_axis(a2);
            plane_intersection = min_a1 <= i_a1 && i_a1 <= max_a1 && min_a2 <= i_a2 && i_a2 <= max_a2;
        }

        let left_ordering = ray_origin_value <= current_node.plane_value;

        // if there is none -> check on which side of the plane the ray started and take that side
        if !plane_intersection
        {
            if left_ordering
            {
                let mut left_bounds = *bounds;
                left_bounds.max_bound.set_axis(axis, current_node.plane_value + EPSILON);
                return self.intersect_kd_node(ray, current_t, current_node.left_first, depth + 1, &left_bounds);
            }
            let mut right_bounds = *bounds;
            right_bounds.min_bound.set_axis(axis, current_node.plane_value - EPSILON);
            return self.intersect_kd_node(ray, current_t, current_node.left_first + 1, depth + 1, &right_bounds);
        }

        let t = current_t + t;

        let mut left_bounds = *bounds;
        left_bounds.max_bound.set_axis(axis, current_node.plane_value + EPSILON);
        let mut right_bounds = *bounds;
        right_bounds.min_bound.set_axis(axis, current_node.plane_value - EPSILON);

        if left_ordering
        {
            let intersection = self.intersect_kd_node(ray, current_t, current_node.left_first, depth + 1, &left_bounds);

            if intersection && ray.t < t
            {
                return true;
            }

            return self.intersect_kd_node(ray, t, current_node.left_first + 1, depth + 1, &right_bounds) || intersection;
        }

        let intersection = self.intersect_kd_node(ray, current_t, current_node.left_first + 1, depth + 1, &right_bounds);

        if intersection && ray.t < t
        {
            return true;
        }
        return self.intersect_kd_node(ray, t, current_node.left_first, depth + 1, &left_bounds) || intersection;
    }
}

impl RayHittableObject for KDTree
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

        let mut bounds = self.bounds;
        bounds.max_bound += EPSILON_VECTOR;
        bounds.min_bound -= EPSILON_VECTOR;
        if self.intersect_kd_node(&mut ray_t, t,0, 0, &bounds)
        {
            ray.t = ray_t.t;
            ray.obj_type = RayHittableObjectType::KDTree;
            ray.obj_idx = self.obj_idx;
            ray.sub_obj_idx = ray_t.sub_obj_idx;
        }

        ray.intersection_tests += ray_t.intersection_tests;
    }

    fn get_normal(&self, ray: &Ray, i: &Float3) -> Float3 {
        return self.normals[ray.sub_obj_idx];
    }

    fn get_uv(&self, i: &Float3) -> Float2 {
        return Float2::zero();
    }

    fn is_occluded(&self, ray: &Ray) -> bool {
        let mut shadow = ray.clone();
        shadow.t = 1e34;
        shadow.obj_idx = usize::MAX;
        self.intersect(&mut shadow);
        return shadow.obj_idx != usize::MAX;
    }
}