use crate::math::*;
use crate::objects::aabb::{AABB, intersect_aabb, intersect_aabb_triangle, triangle_is_in_aabb};
use crate::objects::triangle::{intersect_triangle, Triangle};
use crate::objects::mesh::Mesh;
use crate::ray::*;

#[derive(Clone, Copy)]
struct GridTriangle
{
    pub bounds: AABB,               // 24 bytes
    pub internal_triangle: Triangle // 40 bytes
}                                   // 64 bytes


pub struct Grid
{
    pub bounds: AABB,
    pub triangles: Vec<GridTriangle>,
    pub grid: Vec<Vec<usize>>,
    pub t: Mat4,
    pub inv_t: Mat4,
    pub normals: Vec<Float3>,
    pub res_x: usize,
    pub res_y: usize,
    pub res_z: usize,
    pub cell_size: Float3,
    pub obj_idx: usize,
    pub mat_idx: usize
}


struct Clipped
{
    pub vertices: u32,
    pub v: [Float3; 9],
    pub bounds: AABB
}

impl Clipped
{
    pub fn new(triangle: &GridTriangle, bounding_box: &AABB) -> Self
    {
        let mut v = [Float3::zero(); 9];
        v[0] = triangle.internal_triangle.vertex0;
        v[1] = triangle.internal_triangle.vertex1;
        v[2] = triangle.internal_triangle.vertex2;
        let mut vertices: usize = 3;

        for a in 0..3
        {
            let mut ntmp = 0;
            let mut tmp = [Float3::zero(); 9];
            let mut C = Float3::zero();
            let mut v0 = v[vertices - 1];
            let mut v1 = Float3::zero();
            let mut plane = bounding_box.minimal(a);
            let mut x = (v0.get_axis(a) - plane) >= 0.0;
            for i in 0..vertices
            {
                v1 = v[i];
                let d0 = v0.get_axis(a) - plane;
                let d1 = v1.get_axis(a) - plane;
                if x && d1 >= 0.0
                {
                    /* both in */
                    tmp[ntmp] = v1;
                    ntmp += 1;
                }
                else if !x && d1 > 0.0
                {
                    // coming in: emit C and d1
                    C = v0 + (d0 / (d0 - d1)) * (v1 - v0);
                    C.set_axis(a, plane);
                    tmp[ntmp] = C;
                    ntmp += 1;
                    tmp[ntmp] = v1;
                    ntmp += 1;
                    x = true;
                }
                else if x && d1 < 0.0 // going out: emit C
                {
                    C = v0 + (d0 / (d0 - d1)) * (v1 - v0);
                    C.set_axis(a, plane);
                    tmp[ntmp] = C;
                    ntmp += 1;
                    x = false;
                }

                v0 = v1;
            }
            vertices = 0;
            if ntmp < 3
            {
                break;
            }
            v0 = tmp[ntmp - 1];
            plane = bounding_box.maximal(a);
            x = (plane - v0.get_axis(a)) >= 0.0;
            for i in 0..ntmp
            {
                v1 = tmp[i];
                let d0 = plane - v0.get_axis(a);
                let d1 = plane - v1.get_axis(a);
                if x && d1 >= 0.0
                {
                    /* both in */
                    v[vertices] = v1;
                    vertices += 1;
                }
                else if !x && d1 > 0.0
                {
                    // coming in: emit C and d1
                    C = v0 + (d0 / (d0 - d1)) * (v1 - v0);
                    C.set_axis(a,plane);
                    v[vertices] = C;
                    vertices += 1;
                    v[vertices] = v1;
                    vertices += 1;
                    x = true;
                }
                else if x && d1 < 0.0
                {
                    // going out: emit C
                    C = v0 + (d0 / (d0 - d1)) * (v1 - v0);
                    C.set_axis(a,plane);
                    v[vertices] = C;
                    vertices += 1;
                    x = false;
                }

                v0 = v1;
            }
            if vertices < 3
            {
                // nothing or degenerate
                break;
            }
        }
        // calculate bounding box
        let mut bounds = AABB::from_empty();
        for i in 0..vertices
        {
            bounds.grow(&v[i]);
        }

        Clipped
        {
            vertices: vertices as u32,
            v,
            bounds
        }
    }
}

impl Grid
{
    pub fn from_mesh(mesh: &Mesh, res_x: usize, res_y: usize, res_z: usize) -> Grid
    {
        let mut grid = Vec::with_capacity((res_x * res_y * res_z) as usize);
        for i in 0..(res_x * res_y * res_z)
        {
            grid.push(Vec::with_capacity(0));
        }

        let mut triangles: Vec<GridTriangle> = Vec::with_capacity(mesh.triangles.len());

        let mut tri_idx = 0;
        for triangle in &mesh.triangles
        {
            let vertex0 = mesh.vertices[triangle[0]];
            let vertex1 = mesh.vertices[triangle[1]];
            let vertex2 = mesh.vertices[triangle[2]];

            triangles.push(GridTriangle
            {
                bounds: AABB::from_empty(),
                internal_triangle: Triangle {
                    vertex0,
                    vertex1,
                    vertex2,
                    tri_idx,
                },
            });
            tri_idx += 1;
        }

        triangles.iter_mut().for_each(|triangle|
        {
            triangle.bounds = AABB::from_bounds(&triangle.internal_triangle.vertex0, &triangle.internal_triangle.vertex0);
            triangle.bounds.grow(&triangle.internal_triangle.vertex1);
            triangle.bounds.grow(&triangle.internal_triangle.vertex2);
        });

        let mut bounds: AABB = AABB::from_empty();

        for triangle in &triangles
        {
            bounds.min_bound = bounds.min_bound.min(&triangle.bounds.min_bound);
            bounds.max_bound = bounds.max_bound.max(&triangle.bounds.max_bound);
        }

        let cell_size = Float3::from_xyz(
            (bounds.max_bound.x - bounds.min_bound.x) / (res_x as f32),
            (bounds.max_bound.y - bounds.min_bound.y) / (res_y as f32),
            (bounds.max_bound.z - bounds.min_bound.z) / (res_z as f32));


        let mut grid = Grid
        {
            bounds,
            triangles,
            grid,
            t: mesh.t,
            inv_t: mesh.inv_t,
            normals: mesh.triangle_normals.to_vec(),
            obj_idx: mesh.obj_idx,
            mat_idx: mesh.mat_idx,
            res_x,
            res_y,
            res_z,
            cell_size
        };
        grid.construct_grid();
        return grid;
    }

    fn construct_grid(&mut self)
    {
        for triangle in &self.triangles
        {
            let mut cells_0x: usize = 0;
            let mut cells_0y: usize = 0;
            let mut cells_0z: usize = 0;
            let mut cells_1x: usize = 0;
            let mut cells_1y: usize = 0;
            let mut cells_1z: usize = 0;
            let mut cells_2x: usize = 0;
            let mut cells_2y: usize = 0;
            let mut cells_2z: usize = 0;

            self.get_cell_f3(&triangle.internal_triangle.vertex0, &mut cells_0x, &mut cells_0y, &mut cells_0z);
            self.get_cell_f3(&triangle.internal_triangle.vertex1, &mut cells_1x, &mut cells_1y, &mut cells_1z);
            self.get_cell_f3(&triangle.internal_triangle.vertex2, &mut cells_2x, &mut cells_2y, &mut cells_2z);

            let max_x: i32 = (cells_0x.max(cells_1x.max(cells_2x)) as i32) + 1;
            let max_y: i32 = (cells_0y.max(cells_1y.max(cells_2y)) as i32) + 1;
            let max_z: i32 = (cells_0z.max(cells_1z.max(cells_2z)) as i32) + 1;
            let min_x: i32 = (cells_0x.min(cells_1x.min(cells_2x)) as i32) - 1;
            let min_y: i32 = (cells_0y.min(cells_1y.min(cells_2y)) as i32) - 1;
            let min_z: i32 = (cells_0z.min(cells_1z.min(cells_2z)) as i32) - 1;

            let mut assigned = false;

            for z in min_z..(max_z + 1)
            {
                for y in min_y..(max_y + 1)
                {
                    for x in min_x..(max_x + 1)
                    {
                        if x < 0 || x >= self.res_x as i32 ||
                            y < 0 || y >= self.res_y as i32 ||
                            z < 0 || z >= self.res_z as i32
                        {
                            continue;
                        }

                        let x : usize = x as usize;
                        let y : usize = y as usize;
                        let z : usize = z as usize;
                        /*let cell = self.get_cell_aabb(x, y, z);
                        if !(triangle_is_in_aabb(&cell, &triangle.internal_triangle, &self.normals[triangle.internal_triangle.tri_idx as usize]))
                        {
                            continue;
                        }*/
                        self.grid[x + y * self.res_x + z * self.res_x * self.res_y].push(triangle.internal_triangle.tri_idx as usize);
                        assigned = true;
                    }
                }
            }

/*            if !assigned
            {
                println!("triangle was not assigned: {}, {}, {}, {}, {}, {}", min_x, max_x, min_y, max_y, min_z, max_z);
                self.get_cell_f3(&triangle.internal_triangle.vertex0, &mut cells_0x, &mut cells_0y, &mut cells_0z);
                self.get_cell_f3(&triangle.internal_triangle.vertex1, &mut cells_1x, &mut cells_1y, &mut cells_1z);
                self.get_cell_f3(&triangle.internal_triangle.vertex2, &mut cells_2x, &mut cells_2y, &mut cells_2z);
                let cell = self.get_cell_aabb(min_x,min_y,min_z);
                if !(triangle_is_in_aabb(&cell, &triangle.internal_triangle, &self.normals[triangle.internal_triangle.tri_idx as usize]))
                {
                    return;
                }
                self.grid[min_x + min_y * self.res_x + min_z * self.res_x * self.res_y].push(triangle.internal_triangle.tri_idx as usize);
                assigned = true;
            }*/
            if !assigned
            {
                println!("triangle was not assigned: {}, {}, {}, {}, {}, {}", min_x, max_x, min_y, max_y, min_z, max_z);
            }
        }
    }

    #[inline(always)]
    fn intersect_triangles_in_cell(&self, ray: &mut Ray, x: usize, y: usize, z: usize) -> bool
    {
        let mut intersected = false;
        for triangle in &self.grid[x + y * self.res_x + z * self.res_x * self.res_y]
        {
            intersected = intersect_triangle(&self.triangles[*triangle].internal_triangle, ray) || intersected;
        }
        return intersected;
    }

    fn get_cell_aabb(&self, x: usize, y: usize, z: usize) -> AABB
    {
        let cell_index = Float3::from_xyz(x as f32, y as f32, z as f32);
        let min_bound = cell_index * self.cell_size + self.bounds.min_bound;
        let max_bound = min_bound + self.cell_size;
        return AABB{
            min_bound,
            max_bound
        }
    }

    fn get_cell_f3(&self, position: &Float3, x: &mut usize, y: &mut usize, z: &mut usize)
    {
        let rel_position = (*position - self.bounds.min_bound).max(&Float3::zero());
        *x = (rel_position.x / self.cell_size.x).floor() as usize;
        *y = (rel_position.y / self.cell_size.y).floor() as usize;
        *z = (rel_position.z / self.cell_size.z).floor() as usize;
        
        if *x >= self.res_x
        {
            *x = self.res_x - 1;
        }
        if *y >= self.res_y
        {
            *y = self.res_y - 1;
        }
        if *z >= self.res_z
        {
            *z = self.res_z - 1;
        }
    }
}

impl RayHittableObject for Grid
{
    fn intersect(&self, ray: &mut Ray) {
        let origin = transform_position(&ray.origin, &self.inv_t );
        let direction = transform_vector(&ray.direction, &self.inv_t );

        let mut ray_t = Ray::directed_distance(origin, direction, ray.t);
        ray_t.obj_idx = ray.obj_idx;
        ray_t.obj_type = ray.obj_type;

        // check for intersection with boundary of grid
        let t = intersect_aabb(&ray_t, &self.bounds);
        if t == 1e30
        {
            return;
        }

        // draws bounding box
        /*ray.t = ray_t.t;
        ray.obj_idx = self.obj_idx;
        ray.obj_type = RayHittableObjectType::Grid;
        ray.sub_obj_idx = 0;
        return;*/

        //

        // draws bounding cell of first vertex of triangle
        /*let mut x: usize = 0;
        let mut y: usize = 0;
        let mut z: usize = 0;
        self.get_cell_f3(&self.triangles[0].internal_triangle.vertex0, &mut x, &mut y, &mut z);
        let aabb = self.get_cell_aabb(x,y,z);
        let t = intersect_aabb(&ray_t, &aabb);
        if t != 1e30
        {
            ray.t = t;
            ray.obj_idx = self.obj_idx;
            ray.obj_type = RayHittableObjectType::Grid;
            ray.sub_obj_idx = 0;
            return;
        }

        return;*/

        // draws all bounding boxes that satisfy the condition
        /*for z in 0..self.res_z
        {
            for y in 0..self.res_y
            {
                for x in 0..self.res_x
                {
                    let aabb = self.get_cell_aabb(x,y,z);
                    let t = intersect_aabb(&ray_t, &aabb);
                    if t != 1e30 && self.grid[x + y * self.res_x + z * self.res_x * self.res_z].contains(&0)
                    {
                        ray.t = t;
                        ray.obj_idx = self.obj_idx;
                        ray.obj_type = RayHittableObjectType::Grid;
                        ray.sub_obj_idx = 0;
                        return;
                    }
                }
            }
        }
        return;*/

        let start_position = ray_t.origin + ray_t.direction * (t + 0.0001);

        let ray_orig_grid = ray_t.origin - self.bounds.min_bound;

        let mut delta_t = Float3::zero();

        let mut x: usize = 0;
        let mut y: usize = 0;
        let mut z: usize = 0;
        self.get_cell_f3(&start_position, &mut x, &mut y, &mut z);

        let mut x: i32 = x as i32;
        let mut y: i32 = y as i32;
        let mut z: i32 = z as i32;

        let mut x_inc: i32 = 1;
        let mut y_inc: i32 = 1;
        let mut z_inc: i32 = 1;

        if ray_t.direction.x < 0.0
        {
            delta_t.x = -self.cell_size.x / ray_t.direction.x;
            x_inc = -1;
        }
        else {
            delta_t.x = self.cell_size.x / ray_t.direction.x;
        }
        if ray_t.direction.y < 0.0
        {
            delta_t.y = -self.cell_size.y / ray_t.direction.y;
            y_inc = -1;
        }
        else {
            delta_t.y = self.cell_size.y / ray_t.direction.y;
        }
        if ray_t.direction.z < 0.0
        {
            delta_t.z = -self.cell_size.z / ray_t.direction.z;
            z_inc = -1;
        }
        else {
            delta_t.z = self.cell_size.z / ray_t.direction.z;
        }

        let mut t_x = t + delta_t.x;
        let mut t_y = t + delta_t.y;
        let mut t_z = t + delta_t.z;

        let mut t: f32 = 0.0;

        let mut intersected= false;

        loop
        {
            if self.intersect_triangles_in_cell(&mut ray_t, x as usize, y as usize, z as usize)
            {
                intersected = true;
            }

            if t_x < t_y && t_x < t_z
            {
                t = t_x;
                t_x += delta_t.x;
                x += x_inc;
                if x < 0 || x >= self.res_x as i32
                {
                    break;
                }
            }
            else if t_y < t_z
            {
                t = t_y;
                t_y += delta_t.y;
                y += y_inc;
                if y < 0 || y >= self.res_y as i32
                {
                    break;
                }
            }
            else
            {
                t = t_z;
                t_z += delta_t.z;
                z += z_inc;
                if z < 0 || z >= self.res_z as i32
                {
                    break;
                }
            }

            if intersected
            {
                break;
            }
        }


        if intersected
        {
            ray.t = ray_t.t;
            ray.obj_idx = self.obj_idx;
            ray.obj_type = RayHittableObjectType::Grid;
            ray.sub_obj_idx = ray_t.sub_obj_idx;
        }

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