use crate::bitvec::BitVector;
use crate::math::*;
use crate::objects::aabb::{AABB, intersect_aabb};
use crate::objects::triangle::{intersect_triangle, Triangle};
use crate::objects::mesh::Mesh;
use crate::ray::*;

#[derive(Clone, Copy)]
struct GridTriangle
{
    pub internal_triangle: Triangle // 40 bytes
}


pub struct Grid
{
    pub bounds: AABB,
    triangles: Vec<GridTriangle>,
    pub grid: Vec<Vec<usize>>,
    pub res_x: usize,
    pub res_y: usize,
    pub res_z: usize,
    pub cell_size: Float3,
}

impl Grid
{
    pub fn from_mesh(triangles: &Vec<[usize; 3]>, vertices: &Vec<Float3>, bounds: &AABB, res_x: usize, res_y: usize, res_z: usize) -> Grid
    {
        let mut grid = Vec::with_capacity((res_x * res_y * res_z) as usize);
        for i in 0..(res_x * res_y * res_z)
        {
            grid.push(Vec::with_capacity(0));
        }

        let mut grid_triangles: Vec<GridTriangle> = Vec::with_capacity(triangles.len());

        let mut tri_idx = 0;
        for triangle in triangles
        {
            let vertex0 = vertices[triangle[0]];
            let vertex1 = vertices[triangle[1]];
            let vertex2 = vertices[triangle[2]];

            grid_triangles.push(GridTriangle
            {
                internal_triangle: Triangle {
                    vertex0,
                    vertex1,
                    vertex2,
                    tri_idx,
                },
            });
            tri_idx += 1;
        }

        let bounds = *bounds;

        let cell_size = Float3::from_xyz(
            (bounds.max_bound.x - bounds.min_bound.x) / (res_x as f32),
            (bounds.max_bound.y - bounds.min_bound.y) / (res_y as f32),
            (bounds.max_bound.z - bounds.min_bound.z) / (res_z as f32));


        let mut grid = Grid
        {
            bounds,
            triangles: grid_triangles,
            grid,
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
    fn intersect_triangles_in_cell(&self, ray: &mut Ray, x: usize, y: usize, z: usize, mailbox: &mut BitVector) -> bool
    {
        let mut intersected = false;
        for triangle in &self.grid[x + y * self.res_x + z * self.res_x * self.res_y]
        {
            let id = *triangle;
            if mailbox.get_unchecked(id)
            {
                continue;
            }
            mailbox.set_true(id);
            intersected = intersect_triangle(&self.triangles[id].internal_triangle, ray) || intersected;
        }
        return intersected;
    }

    #[allow(dead_code)]
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

    pub fn intersect(&self, ray_t: &mut Ray, t: f32) -> bool {

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

        let mut mailbox = BitVector::new(self.triangles.len());

        let start_position = ray_t.origin + ray_t.direction * (t + 0.0001);

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

        let mut t: f32;

        let mut intersected= false;

        loop
        {
            if self.intersect_triangles_in_cell(ray_t, x as usize, y as usize, z as usize, &mut mailbox)
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

            if intersected && ray_t.t < t
            {
                break;
            }
        }

        return intersected;
    }
}