use std::f32::consts::PI;
use crate::material::*;
use crate::math::*;

#[derive(Debug, Clone, Copy)]
pub enum RayHittableObjectType
{
    None,
    Sphere,
    Plane,
    Cube,
    Torus,
    Quad,
    Mesh,
    Bvh
}

#[repr(C, align(64))]
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    pub origin: Float3,
    pub direction: Float3,
    pub  r_direction: Float3,
    pub t: f32,
    pub obj_idx: i32,
    pub sub_obj_idx: usize,
    pub obj_type: RayHittableObjectType,
    pub inside: bool
}

impl Ray {

    #![allow(dead_code)]
    pub fn new() -> Self
    {
        Ray {
            origin: Float3::zero(),
            direction: Float3::zero(),
            r_direction: Float3::zero(),
            t: 1.0e34,
            obj_idx: -1,
            sub_obj_idx: 0,
            obj_type: RayHittableObjectType::None,
            inside: false
        }
    }

    pub fn directed(origin: Float3, direction: Float3) -> Self
    {
        Ray {
            origin: Float3::from_xyz(origin.x, origin.y, origin.z),
            direction: Float3::from_xyz(direction.x, direction.y, direction.z),
            r_direction: Float3::from_xyz(1.0 / direction.x, 1.0 / direction.y, 1.0 / direction.z),
            t: 1.0e34,
            obj_idx: -1,
            sub_obj_idx: 0,
            obj_type: RayHittableObjectType::None,
            inside: false
        }
    }

    #[allow(dead_code)]
    pub fn directed_distance(origin: Float3, direction: Float3, distance: f32) -> Self
    {
        Ray {
            origin: Float3::from_xyz(origin.x, origin.y, origin.z),
            direction: Float3::from_xyz(direction.x, direction.y, direction.z),
            r_direction: Float3::from_xyz(1.0 / direction.x, 1.0 / direction.y, 1.0 / direction.z),
            t: distance,
            obj_idx: -1,
            sub_obj_idx: 0,
            obj_type: RayHittableObjectType::None,
            inside: false
        }
    }

    pub fn intersection_point(&self) -> Float3
    {
        return self.origin + self.direction * self.t;
    }
}

trait RayHittableObject
{
    fn intersect(&self, ray: &mut Ray);

    fn get_normal(&self, i: &Float3) -> Float3;

    fn get_uv(&self, i: &Float3) -> Float2;

    fn is_occluded(&self, ray: &Ray) -> bool;
}

trait RayOccluder
{
    fn is_occluded(&self, ray: &Ray);
}

pub struct Sphere
{
    pub position: Float3,
    pub radius: f32,
    pub radius2: f32,
    pub inv_radius: f32,
    pub obj_idx: i32,
    pub mat_idx: i32,
}

impl RayHittableObject for Sphere
{
    fn intersect(&self, ray: &mut Ray) {
        let oc = ray.origin - self.position;
        let b = dot(&oc, &ray.direction);
        let c = dot(&oc, &oc) - self.radius2;
        let mut d = b * b - c;
        if d <= 0.0
        {
            return;
        }

        d = d.sqrt();
        let mut t = -b - d;
        let mut hit = t < ray.t && t > 0.0;
        if hit
        {
            ray.t = t;
            ray.obj_idx = self.obj_idx;
            ray.obj_type = RayHittableObjectType::Sphere;
            return;
        }
        if c > 0.0
        {
            return;
        }

        t = d - b;
        hit = t < ray.t && t > 0.0;
        if hit
        {
            ray.t = t;
            ray.obj_idx = self.obj_idx;
            ray.obj_type = RayHittableObjectType::Sphere;
            return;
        }
    }

    fn get_normal(&self, i: &Float3) -> Float3 {
        return (*i - self.position) * self.inv_radius;
    }

    fn get_uv(&self, i: &Float3) -> Float2
    {
        Float2::zero()
    }

    fn is_occluded(&self, ray: &Ray) -> bool {
        let oc = ray.origin - self.position;
        let b = dot(&oc, &ray.direction);
        let c = dot(&oc, &oc) - self.radius2;

        let d = b * b - c;
        if d <= 0.0
        {
            return false
        }

        let t = -b - d.sqrt();
        return t < ray.t && t > 0.0;
    }
}

impl Sphere
{
    fn new(obj_idx: i32, mat_idx: i32, position: Float3, radius: f32) -> Self
    {
        Sphere
        {
            obj_idx,
            mat_idx,
            position,
            radius,
            radius2: radius * radius,
            inv_radius: 1.0 / radius
        }
    }
}

pub struct PlaneUVFunction
{
    f: fn(i: &Float3) -> Float2
}

impl PlaneUVFunction
{
    pub fn empty() -> Self
    {
        PlaneUVFunction
        {
            f: |_: &Float3| Float2::zero()
        }
    }

    pub fn xz() -> Self
    {
        PlaneUVFunction
        {
            f: |i: &Float3| {
                Float2::from_xy(i.x, i.z)
            }
        }
    }

    pub fn xy() -> Self
    {
        PlaneUVFunction
        {
            f: |i: &Float3| {
                Float2::from_xy(i.x, i.y)
            }
        }
    }

    pub fn zy() -> Self
    {
        PlaneUVFunction
        {
            f: |i: &Float3| {
                Float2::from_xy(i.z, i.y)
            }
        }
    }
}

pub struct Plane
{
    pub obj_idx: i32,
    pub mat_idx: i32,
    pub normal: Float3,
    pub distance: f32,
    pub uv_function: PlaneUVFunction
}

impl RayHittableObject for Plane
{
    fn intersect(&self, ray: &mut Ray) {
        let t = -(dot(&ray.origin, &self.normal) + self.distance) / (dot(&ray.direction, &self.normal));
        if t < ray.t && t > 0.0
        {
            ray.t = t;
            ray.obj_idx = self.obj_idx;
            ray.obj_type = RayHittableObjectType::Plane;
        }
    }

    fn get_normal(&self, _: &Float3) -> Float3 {
        return self.normal;
    }

    fn get_uv(&self, i: &Float3) -> Float2
    {
        return (self.uv_function.f)(i)
    }

    fn is_occluded(&self, ray: &Ray) -> bool {
        return false;
    }
}

impl Plane
{
    pub fn new(obj_idx: i32, mat_idx: i32, normal: Float3, distance: f32, uv_function: PlaneUVFunction) -> Self
    {
        Plane{
            obj_idx,
            mat_idx,
            normal,
            distance,
            uv_function
        }
    }
}

pub struct Cube
{
    pub obj_idx: i32,
    pub mat_idx: i32,
    pub b: [Float3; 2],
    pub m: Mat4,
    pub inv_m: Mat4
}

impl RayHittableObject for Cube
{
    fn intersect(&self, ray: &mut Ray) {
        let origin = transform_position(&ray.origin, &self.inv_m);
        let direction = transform_vector(&ray.direction, &self.inv_m);
        let rdx = 1.0 / direction.x;
        let rdy = 1.0 / direction.y;
        let rdz = 1.0 / direction.z;
        let signx = (direction.x < 0.0) as usize;
        let signy = (direction.y < 0.0) as usize;
        let signz = (direction.z < 0.0) as usize;

        let mut tmin = (self.b[signx].x - origin.x) * rdx;
        let mut tmax = (self.b[1 - signx].x - origin.x) * rdx;
        let tymin = (self.b[signy].y - origin.y) * rdy;
        let tymax = (self.b[1 - signy].y - origin.y) * rdy;

        if tmin > tymax || tymin > tmax
        {
            return;
        }

        tmin = tmin.max(tymin);
        tmax = tmax.min(tymax);
        let tzmin = (self.b[signz].z - origin.z) * rdz;
        let tzmax = (self.b[1 - signz].z - origin.z) * rdz;

        if tmin > tzmax || tzmin > tmax
        {
            return;
        }

        tmin = tmin.max(tzmin);
        tmax = tmax.min(tzmax);
        if tmin > 0.0
        {
            if tmin < ray.t
            {
                ray.t = tmin;
                ray.obj_idx = self.obj_idx;
                ray.obj_type = RayHittableObjectType::Cube;
            }
        }
        else if tmax > 0.0
        {
            if tmax < ray.t
            {
                ray.t = tmax;
                ray.obj_idx = self.obj_idx;
                ray.obj_type = RayHittableObjectType::Cube;
            }
        }
    }

    fn get_normal(&self, i: &Float3) -> Float3 {
        let obj_i = transform_position(&i, &self.inv_m);
        let mut n = Float3::from_xyz(-1.0, 0.0, 0.0);
        let d0 = (obj_i.x - self.b[0].x).abs();
        let d1 = (obj_i.x - self.b[1].x).abs();
        let d2 = (obj_i.y - self.b[0].x).abs();
        let d3 = (obj_i.y - self.b[1].x).abs();
        let d4 = (obj_i.z - self.b[0].x).abs();
        let d5 = (obj_i.z - self.b[1].x).abs();
        let mut min_dist = d0;
        if d1 < min_dist
        {
            min_dist = d1;
            n.x = 1.0;
        }
        if d2 < min_dist
        {
            min_dist = d2;
            n = Float3::from_xyz(0.0, -1.0, 0.0);
        }
        if d3 < min_dist
        {
            min_dist = d3;
            n = Float3::from_xyz(0.0, 1.0, 0.0);
        }
        if d4 < min_dist
        {
            min_dist = d4;
            n = Float3::from_xyz(0.0, 0.0, -1.0);
        }
        if d5 < min_dist
        {
            n = Float3::from_xyz(0.0, 0.0, 1.0);
        }

        return transform_vector(&n, &self.m);
    }

    fn get_uv(&self, i: &Float3) -> Float2
    {
        Float2::zero()
    }

    fn is_occluded(&self, ray: &Ray) -> bool {
        let origin = transform_position(&ray.origin, &self.inv_m);
        let direction = transform_vector(&ray.direction, &self.inv_m);
        let rdx = 1.0 / direction.x;
        let rdy = 1.0 / direction.y;
        let rdz = 1.0 / direction.z;
        let t1 = (self.b[0].x - origin.x) * rdx;
        let t2 = (self.b[1].x - origin.x) * rdx;
        let t3 = (self.b[0].y - origin.y) * rdy;
        let t4 = (self.b[1].y - origin.y) * rdy;
        let t5 = (self.b[0].z - origin.z) * rdz;
        let t6 = (self.b[1].z - origin.z) * rdz;
        let tmin = t1.min( t2 ).max( t3.min( t4 ) ).max( t5.min( t6 ) );
        let tmax = t1.max( t2 ).min( t3.max( t4 ) ).min( t5.max( t6 ) );
        return tmax > 0.0 && tmin < tmax && tmin < ray.t;
    }
}

impl Cube
{
    pub fn new(obj_idx: i32, mat_idx: i32, pos: Float3, size: Float3) -> Self
    {
        Cube{
            obj_idx,
            mat_idx,
            b:[
                pos - size * 0.5,
                pos + size * 0.5
            ],
            m: Mat4::identity_matrix(),
            inv_m: Mat4::identity_matrix().inverted()
        }
    }
}

pub struct Quad
{
    pub size: f32,
    pub t: Mat4,
    pub inv_t: Mat4,
    pub obj_idx: i32,
    pub mat_idx: i32,
}

impl Quad
{
    pub fn new(obj_idx: i32, mat_idx: i32, size: f32, transform: &Mat4) -> Quad
    {
        Quad{
            obj_idx,
            mat_idx,
            size: size * 0.5,
            t: transform.clone(),
            inv_t: transform.inverted_no_scale()
        }
    }

    pub fn random_points(&self, num: usize, seed: &mut u32) -> Vec<Float3>
    {
        let mut vector: Vec<Float3> = Vec::with_capacity(num);
        let size = self.size;
        let corner1 = transform_position(&Float3::from_xyz(-size, 0.0, -size), &self.t);
        let corner2 = transform_position(&Float3::from_xyz(size, 0.0, -size), &self.t);
        let corner3 = transform_position(&Float3::from_xyz(-size, 0.0, size), &self.t);

        for _ in 0..num
        {
            let r1 = random_float_s(seed);
            let r2 = random_float_s(seed);
            vector.push(corner1 + r2 * (corner2 - corner1) + r1 * (corner3 - corner1));
        }
        return vector;
    }

    pub fn center_point(&self) -> Float3
    {
        let r1 = 0.5;
        let r2 = 0.5;
        let size = self.size;
        let corner1 = transform_position(&Float3::from_xyz(-size, 0.0, -size), &self.t);
        let corner2 = transform_position(&Float3::from_xyz(size, 0.0, -size), &self.t);
        let corner3 = transform_position(&Float3::from_xyz(-size, 0.0, size), &self.t);
        return corner1 + r2 * (corner2 - corner1) + r1 * (corner3 - corner1);
    }
}

impl RayHittableObject for Quad
{
    fn intersect(&self, ray: &mut Ray) {
        let o_y: f32 = self.inv_t.cell[4] * ray.origin.x + self.inv_t.cell[5] * ray.origin.y + self.inv_t.cell[6] * ray.origin.z + self.inv_t.cell[7];
        let d_y: f32 = self.inv_t.cell[4] * ray.direction.x + self.inv_t.cell[5] * ray.direction.y + self.inv_t.cell[6] * ray.direction.z;
        let t = o_y / -d_y;
        if t < ray.t && t > 0.0
        {
            let o_x: f32 = self.inv_t.cell[0] * ray.origin.x + self.inv_t.cell[1] * ray.origin.y + self.inv_t.cell[2] * ray.origin.z + self.inv_t.cell[3];
            let o_z: f32 = self.inv_t.cell[8] * ray.origin.x + self.inv_t.cell[9] * ray.origin.y + self.inv_t.cell[10] * ray.origin.z + self.inv_t.cell[11];
            let d_x: f32 = self.inv_t.cell[0] * ray.direction.x + self.inv_t.cell[1] * ray.direction.y + self.inv_t.cell[2] * ray.direction.z;
            let d_z: f32 = self.inv_t.cell[8] * ray.direction.x + self.inv_t.cell[9] * ray.direction.y + self.inv_t.cell[10] * ray.direction.z;
            let i_x = o_x + t * d_x;
            let i_z = o_z + t * d_z;
            if i_x > -self.size && i_x < self.size && i_z > -self.size && i_z < self.size
            {
                ray.t = t;
                ray.obj_idx = self.obj_idx;
                ray.obj_type = RayHittableObjectType::Quad
            }
        }
    }

    fn get_normal(&self, i: &Float3) -> Float3 {
        Float3::from_xyz(-self.t.cell[1], -self.t.cell[5], -self.t.cell[9])
    }

    fn get_uv(&self, i: &Float3) -> Float2 {
        Float2::zero()
    }

    fn is_occluded(&self, ray: &Ray) -> bool {
        let o_y: f32 = self.inv_t.cell[4] * ray.origin.x + self.inv_t.cell[5] * ray.origin.y + self.inv_t.cell[6] * ray.origin.z + self.inv_t.cell[7];
        let d_y: f32 = self.inv_t.cell[4] * ray.direction.x + self.inv_t.cell[5] * ray.direction.y + self.inv_t.cell[6] * ray.direction.z;
        let t = o_y / -d_y;
        if t < ray.t && t > 0.0
        {
            let o_x: f32 = self.inv_t.cell[0] * ray.origin.x + self.inv_t.cell[1] * ray.origin.y + self.inv_t.cell[2] * ray.origin.z + self.inv_t.cell[3];
            let o_z: f32 = self.inv_t.cell[8] * ray.origin.x + self.inv_t.cell[9] * ray.origin.y + self.inv_t.cell[10] * ray.origin.z + self.inv_t.cell[11];
            let d_x: f32 = self.inv_t.cell[0] * ray.direction.x + self.inv_t.cell[1] * ray.direction.y + self.inv_t.cell[2] * ray.direction.z;
            let d_z: f32 = self.inv_t.cell[8] * ray.direction.x + self.inv_t.cell[9] * ray.direction.y + self.inv_t.cell[10] * ray.direction.z;
            let i_x = o_x + t * d_x;
            let i_z = o_z + t * d_z;
            return i_x > -self.size && i_x < self.size && i_z > -self.size && i_z < self.size;
        }

        return false;
    }
}

pub struct Torus
{
    pub rt2: f32,
    pub rc2: f32,
    pub r2: f32,
    pub obj_idx: i32,
    pub mat_idx: i32,
    pub t: Mat4,
    pub inv_t: Mat4
}

impl Torus
{
    pub fn new(obj_idx: i32, mat_idx: i32, a: f32, b: f32) -> Self
    {
        Torus
        {
            obj_idx,
            mat_idx,
            rc2: a * a,
            rt2: b * b,
            t: Mat4::identity_matrix(),
            inv_t: Mat4::identity_matrix(),
            r2: (a + b).sqrt()
        }
    }

    fn cbrt_fast(n: f64) -> f64
    {
        let mut x1 = n / 10.0;
        let mut x2 = 1.0;
        let mut turn = 0;
        while (x1 - x2).abs() > 0.00000001 && turn < 100
        {
            turn += 1;
            x1 = x2;
            x2 = (2.0 / 3.0 * x1) + (n / (3.0 * x1 * x1));
        }
        return x2;
    }

    fn cbrt_fast_f32(n: f32) -> f32
    {
        let mut x1 = n / 10.0;
        let mut x2 = 1.0;
        let mut turn = 0;
        while (x1 - x2).abs() > 0.00000001 && turn < 100
        {
            turn += 1;
            x1 = x2;
            x2 = (2.0 / 3.0 * x1) + (n / (3.0 * x1 * x1));
        }
        return x2;
    }
}

impl RayHittableObject for Torus
{
    fn intersect(&self, ray: &mut Ray) {
        // via: https://www.shadertoy.com/view/4sBGDy
        let origin: Float3 = transform_position(&ray.origin, &self.inv_t );
        let distance: Float3 = transform_vector(&ray.direction, &self.inv_t);

        let r2 = self.r2 as f64;
        let rt2 = self.rt2 as f64;
        let rc2 = self.rc2 as f64;

        // extension rays need double precision for the quadratic solver!
        let mut po = 1.0;
        let m = dot(&origin, &origin) as f64;
        let mut k3 = dot(&origin, &distance) as f64;
        let mut k32 = k3 * k3;

        // bounding sphere test
        let v = k32 - m + r2;
        if v < 0.0
        {
            return;
        }

        // setup torus intersection
        let k = (m - rt2 - rc2) * 0.5;
        let mut k2 = k32 + rc2 * ((distance.z * distance.z) as f64) + k;
        let mut k1 = k * k3 + rc2 * ((origin.z * distance.z) as f64);
        let mut k0 = k * k + rc2 * ((origin.z * origin.z) as f64) - rc2 * rt2;
        // solve quadratic equation
        if ( k3 * (k32 - k2) + k1 ).abs() < 0.0001
        {
            let temp = k1;
            k1 = k3;
            k3 = temp;
            po = -1.0;
            k0 = 1.0 / k0;
            k1 = k1 * k0;
            k2 = k2 * k0;
            k3 = k3 * k0;
            k32 = k3 * k3;
        }
        let mut c2 = 2.0 * k2 - 3.0 * k32;
        let mut c1 = k3 * (k32 - k2) + k1;
        let mut c0 = k3 * (k3 * (-3.0 * k32 + 4.0 * k2) - 8.0 * k1) + 4.0 * k0;
        c2 *= 0.33333333333;
        c1 *= 2.0;
        c0 *= 0.33333333333;
        let q = c2 * c2 + c0;
        let r = 3.0 * c0 * c2 - c2 * c2 * c2 - c1 * c1;
        let mut h = r * r - q * q * q;
        let mut z: f64;
        if h < 0.0
        {
            let s_q = q.sqrt();
            z = 2.0 * s_q * ((r / (s_q * q)).acos() * 0.33333333333).cos();
        }
        else
        {
            let s_q = Torus::cbrt_fast(h.sqrt() + r.abs());
            z = ( s_q + q / s_q).abs().copysign(r);
        }
        z = c2 - z;
        let mut d1 = z - 3.0 * c2;
        let mut d2 = z * z - 3.0 * c0;
        if d1.abs() < 1.0e-8
        {
            if d2 < 0.0
            {
                return;
            }
            d2 = d2.sqrt();
        }
        else
        {
            if d1 < 0.0
            {
                return;
            }
            d1 = ( d1 * 0.5 ).sqrt();
            d2 = c1 / d1;
        }
        let mut t = 1e20;
        h = d1 * d1 - z + d2;
        if h > 0.0
        {
            h = h.sqrt();
            let mut t1 = -d1 - h - k3;
            let mut t2 = -d1 + h - k3;

            if po < 0.0 {
                t1 = 2.0 / t1;
                t2 = 2.0 / t2;
            }
            if t1 > 0.0
            {
                t = t1;
            }
            if t2 > 0.0
            {
                t = t.min(t2);
            }
        }
        h = d1 * d1 - z - d2;
        if h > 0.0
        {
            h = h.sqrt();
            let mut t1 = d1 - h - k3;
            let mut t2 = d1 + h - k3;
            if po < 0.0 {
                t1 = 2.0 / t1;
                t2 = 2.0 / t2;
            }
            if t1 > 0.0
            {
                t = t.min(t1);
            }
            if t2 > 0.0
            {
                t = t.min(t2);
            }
        }
        let ft = t as f32;
        if ft > 0.0 && ft < ray.t
        {
            ray.t = ft;
            ray.obj_idx = self.obj_idx;
            ray.obj_type = RayHittableObjectType::Torus;
        }
    }

    fn get_normal(&self, i: &Float3) -> Float3 {
        let l = transform_position(&i, &self.inv_t );
        let x = l * (- self.rc2 * Float3::from_xyz(1.0, 1.0, -1.0 ) + dot(&l, &l) - self.rt2);
        let n = normalize( &x );
        return transform_vector(&n, &self.t );
    }

    fn get_uv(&self, i: &Float3) -> Float2 {
        Float2::zero()
    }

    fn is_occluded(&self, ray: &Ray) -> bool {
        // via: https://www.shadertoy.com/view/4sBGDy
        let origin = transform_position(&ray.origin, &self.inv_t );
        let direction = transform_vector(&ray.direction, &self.inv_t );
        let mut po = 1.0;
        let mut m = dot(&origin, &origin);
        let mut k3 = dot(&origin, &direction);
        let mut k32 = k3 * k3;

        let r2 = self.r2;
        let rt2 = self.rt2;
        let rc2 = self.rc2;

        // bounding sphere test
        let v = k32 - m + r2;
        if v < 0.0
        {
            return false;
        }

        // setup torus intersection
        let mut k = (m - rt2 - rc2) * 0.5;
        let mut k2 = k32 + rc2 * direction.z * direction.z + k;
        let mut k1 = k * k3 + rc2 * origin.z * direction.z;
        let mut k0 = k * k + rc2 * origin.z * origin.z - rc2 * rt2;
        // solve quadratic equation
        if (k3 * (k32 - k2) + k1 ).abs() < 0.01
        {
            let temp = k1;
            k1 = k3;
            k3 = temp;
            po = -1.0;
            k0 = 1.0 / k0;
            k1 = k1 * k0;
            k2 = k2 * k0;
            k3 = k3 * k0;
            k32 = k3 * k3;
        }
        let mut c2 = 2.0 * k2 - 3.0 * k32;
        let mut c1 = k3 * (k32 - k2) + k1;
        let mut c0 = k3 * (k3 * (-3.0 * k32 + 4.0 * k2) - 8.0 * k1) + 4.0 * k0;
        c2 *= 0.33333333333;
        c1 *= 2.0;
        c0 *= 0.33333333333;
        let Q = c2 * c2 + c0;
        let R = 3.0 * c0 * c2 - c2 * c2 * c2 - c1 * c1;
        let mut h = R * R - Q * Q * Q;
        let mut z: f32;
        if h < 0.0
        {
            let s_q = Q.sqrt();
            z = 2.0 * s_q * ((R / (s_q * Q)).acos() * 0.33333333333).cos();
        }
        else
        {
            let s_q = Torus::cbrt_fast_f32(h.sqrt() + R.abs());
            z = ( s_q + Q / s_q).abs().copysign(R);
        }
        z = c2 - z;
        let mut d1 = z - 3.0 * c2;
        let mut d2 = z * z - 3.0 * c0;
        if d1.abs() < 1.0e-4
        {
            if d2 < 0.0
            {
                return false;
            }
            d2 = d2.sqrt();
        }
        else
        {
            if d1 < 0.0
            {
                return false;
            }
            d1 = ( d1 * 0.5 ).sqrt();
            d2 = c1 / d1;
        }
        let t = 1e20;
        h = d1 * d1 - z + d2;
        if h > 0.0
        {
            let mut t1 = -d1 - h.sqrt() - k3;
            if po < 0.0
            {
                t1 = 2.0 / t1;
            }
            if t1 > 0.0 && t1 < ray.t
            {
                return true;
            }
        }
        h = d1 * d1 - z - d2;
        if h > 0.0
        {
            let mut t1 = d1 - h.sqrt() - k3;
            if po < 0.0
            {
                t1 = 2.0 / t1;
            }
            if t1 > 0.0 && t1 < ray.t
            {
                return true;
            }
        }
        return false;
    }
}

#[derive(Clone, Copy)]
pub struct Triangle
{
    pub bounds: AABB,     // 24 bytes
    pub tri_idx: i32,     //  4 bytes
    pub vertex0: Float3,  // 12 bytes
    pub vertex1: Float3,  // 12 bytes
    pub vertex2: Float3   // 12 bytes
}                         // 64 bytes

#[derive(Clone, Copy)]
pub struct BVHNode
{
    pub bounds: AABB,
    pub tri_count: usize,
    pub left_first: usize,
}

impl BVHNode
{
    pub fn is_leaf(&self) -> bool
    {
        return self.tri_count > 0;
    }
}

pub struct Clipped
{
    pub vertices: u32,
    pub v: [Float3; 9],
    pub bounds: AABB
}

impl Clipped
{
    pub fn new(triangle: &Triangle, bounding_box: &AABB) -> Self
    {
        let mut v = [Float3::zero(); 9];
        v[0] = triangle.vertex0;
        v[1] = triangle.vertex1;
        v[2] = triangle.vertex2;
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

pub struct BVH
{
    pub triangles: Vec<Triangle>,
    pub bvh_nodes: Vec<BVHNode>,
    pub triangle_idx: Vec<usize>,
    pub triangle_tmp: Vec<usize>,
    pub root_node_idx: usize,
    pub nodes_used: usize,
    pub triangle_ptr: usize,
    pub spatial_splits: usize,
    pub t: Mat4,
    pub inv_t: Mat4,
    pub normals: Vec<Float3>,
    pub obj_idx: i32,
    pub mat_idx: i32
}

const BVH_BINS: usize = 128;

impl BVH
{
    pub fn from_mesh(mesh: &Mesh) -> Self
    {
        let prim_count = mesh.triangles.len();
        let mut triangles: Vec<Triangle> = Vec::with_capacity(prim_count * 2);
        let mut triangle_idx: Vec<usize> = Vec::with_capacity(prim_count * 2);
        let mut triangle_tmp: Vec<usize> = Vec::with_capacity(prim_count * 2);

        let mut tri_idx = 0;
        for triangle in &mesh.triangles
        {
            let vertex0 = mesh.vertices[triangle[0]];
            let vertex1 = mesh.vertices[triangle[1]];
            let vertex2 = mesh.vertices[triangle[2]];

            triangles.push(Triangle
            {
                bounds: AABB::from_empty(),
                vertex0,
                vertex1,
                vertex2,
                tri_idx,
            });
            triangle_idx.push(tri_idx as usize);
            triangle_tmp.push(0);
            tri_idx += 1;
        }

        // build bvh
        let mut bvh_nodes: Vec<BVHNode> = Vec::with_capacity(prim_count * 4);

        triangles.iter_mut().for_each(|triangle|
        {
            triangle.bounds = AABB::from_bounds(&triangle.vertex0, &triangle.vertex0);
            triangle.bounds.grow(&triangle.vertex1);
            triangle.bounds.grow(&triangle.vertex2);
        });

        let mut root = BVHNode {
            left_first: 0,
            tri_count: prim_count,
            bounds: AABB::from_empty()
        };

        for triangle in &triangles
        {
            root.bounds.min_bound = root.bounds.min_bound.min(&triangle.bounds.min_bound);
            root.bounds.max_bound = root.bounds.max_bound.max(&triangle.bounds.max_bound);
        }

        bvh_nodes.push(root);
        for _ in 1..prim_count * 4
        {
            bvh_nodes.push(
                BVHNode {
                    left_first: 0,
                    tri_count: 0,
                    bounds: AABB::from_empty()
                }
            );
        }

        let root_node_idx: usize = 0;

        let spatial_splits: usize = 0;

        for _ in &mesh.triangles
        {
            triangles.push(Triangle
            {
                bounds: AABB::from_empty(),
                vertex0: Float3::zero(),
                vertex1: Float3::zero(),
                vertex2: Float3::zero(),
                tri_idx: 0,
            });
            triangle_idx.push(0);
            triangle_tmp.push(0);
        }

        let mut bvh = BVH
        {
            triangles,
            bvh_nodes,
            triangle_idx,
            triangle_tmp,
            root_node_idx,
            nodes_used: 2,
            triangle_ptr: prim_count,
            spatial_splits,
            t: mesh.t,
            inv_t: mesh.inv_t,
            normals: mesh.triangle_normals.to_vec(),
            obj_idx: mesh.obj_idx,
            mat_idx: mesh.mat_idx
        };

        bvh.subdivide(root_node_idx, prim_count);
        bvh.finalize_sbvh(0);

        println!("SAH:            {}", bvh.get_total_sah(0));
        println!("nodes:          {}", bvh.get_node_count(0));
        println!("leafs:          {}", bvh.get_leaf_count(0));
        println!("spatial splits: {}", bvh.spatial_splits);
        println!("SAH cost:       {}", bvh.sah_cost(0));

        return bvh;
    }

    fn get_total_sah(&self, node_idx: usize) -> f32
    {
        let node = &self.bvh_nodes[node_idx];
        if node.is_leaf()
        {
            return node.bounds.area();
        }
        return self.get_total_sah(node.left_first) + self.get_total_sah(node.left_first + 1);
    }

    fn get_node_count(&self, node_idx: usize) -> u32
    {
        let node = &self.bvh_nodes[node_idx];
        if node.is_leaf()
        {
            return 1;
        }
        return self.get_node_count(node.left_first) + self.get_node_count(node.left_first + 1) + 1;
    }

    fn get_leaf_count(&self, node_idx: usize) -> u32
    {
        let node = &self.bvh_nodes[node_idx];
        if node.is_leaf()
        {
            return 1;
        }
        return self.get_leaf_count(node.left_first) + self.get_leaf_count(node.left_first + 1);
    }

    fn sah_cost(&self, node_idx: usize) -> f32
    {
        let node = &self.bvh_nodes[node_idx];
        let area = node.bounds.area();
        if node.is_leaf()
        {
            return area * (node.tri_count as f32);
        }
        let mut cost: f32 = 0.0;
        cost += self.sah_cost(node.left_first);
        cost += self.sah_cost(node.left_first + 1);
        cost += area;
        if (node_idx == 0)
        {
            cost *= 1.0 / area;
        }
        return cost;
    }


    fn intersect_triangle(&self, ray: &mut Ray, tri_idx: usize) -> bool
    {
        let triangle = &self.triangles[tri_idx];
        let edge1 = triangle.vertex1 - triangle.vertex0;
        let edge2 = triangle.vertex2 - triangle.vertex0;
        let h = cross(&ray.direction, &edge2);
        let a = dot(&edge1, &h);
        if a  > -0.0001 && a < 0.0001
        {
            return false;
        }
        let s = ray.origin - triangle.vertex0;
        let f = 1.0 / a;
        let u = f * dot(&s, &h);
        if u < 0.0 || u > 1.0
        {
            return false;
        }
        let q = cross(&s, &edge1);
        let v = f * dot(&ray.direction, &q);
        if v < 0.0 || u + v > 1.0
        {
            return false;
        }
        let t = f * dot(&edge2, &q);
        if t > 0.0001 && t < ray.t
        {
            ray.t = t;
            ray.sub_obj_idx = triangle.tri_idx as usize;
            return true;
        }
        return false;
    }

    fn intersect_aabb(&self, ray: &Ray, aabb: &AABB) -> f32
    {
        let tx1 = (aabb.min_bound.x - ray.origin.x) * ray.r_direction.x;
        let tx2 = (aabb.max_bound.x - ray.origin.x) * ray.r_direction.x;
        let mut tmin = tx1.min(tx2);
        let mut tmax = tx1.max(tx2);
        let ty1 = (aabb.min_bound.y - ray.origin.y) * ray.r_direction.y;
        let ty2 = (aabb.max_bound.y - ray.origin.y) * ray.r_direction.y;
        tmin = tmin.max(ty1.min(ty2));
        tmax = tmax.min(ty1.max(ty2));
        let tz1 = (aabb.min_bound.z - ray.origin.z) * ray.r_direction.z;
        let tz2 = (aabb.max_bound.z - ray.origin.z) * ray.r_direction.z;
        tmin = tmin.max(tz1.min(tz2));
        tmax = tmax.min(tz1.max(tz2));
        if tmax >= tmin && tmin < ray.t && tmax > 0.0
        {
            return tmin;
        }
        return 1e30;
    }

    fn finalize_sbvh(&mut self, node_idx: usize)
    {
        let node = self.bvh_nodes[node_idx];
        if node.is_leaf()
        {
            for i in 0..node.tri_count
            {
                let idx = self.triangle_idx[node.left_first + i];
                self.triangle_idx[node.left_first + i] = self.triangles[idx].tri_idx as usize;
            }
            return;
        }
        self.finalize_sbvh(node.left_first);
        self.finalize_sbvh(node.left_first + 1);
    }

    fn find_best_object_split_plane(&mut self, node: &BVHNode, axis: &mut usize, split_pos: &mut f32, l_box: &mut AABB, r_box: &mut AABB) -> f32
    {
        let mut best_cost: f32 = 1e30;
        for a in 0..3
        {
            let mut bounds_min: f32 = 1e30;
            let mut bounds_max: f32 = -1e30;
            for i in 0..node.tri_count
            {
                let triangle = &self.triangles[self.triangle_idx[node.left_first + i]];

                let triangle_center = triangle.bounds.center( a );
                bounds_min = bounds_min.min( triangle_center );
                bounds_max = bounds_max.max( triangle_center );
            }
            if bounds_min == bounds_max
            {
                continue;
            }
            // loop over split plane candidates
            let bin_width = (bounds_max - bounds_min) / (BVH_BINS as f32);
            for b in 1..BVH_BINS
            {
                let plane = bounds_min + bin_width * (b as f32);

                let mut left_box = AABB::from_empty();
                let mut right_box = AABB::from_empty();

                let mut left_count: u32 = 0;
                let mut right_count: u32 = 0;
                for i in 0..node.tri_count
                {
                    let triangle = &self.triangles[self.triangle_idx[node.left_first + i]];
                    if triangle.bounds.center( a ) < plane
                    {
                        left_box.grow_aabb( &triangle.bounds );
                        left_count += 1;
                    }
                    else
                    {
                        right_box.grow_aabb( &triangle.bounds );
                        right_count += 1;
                    }
                }
                let plane_cost = (left_count as f32) * left_box.area() + (right_count as f32) * right_box.area();
                if plane_cost < best_cost
                {
                    *axis = a;
                    *split_pos = plane;
                    best_cost = plane_cost;
                    *l_box = left_box;
                    *r_box = right_box;
                }

            }
        }

        return best_cost;
    }

    fn find_best_spatial_split_plane(&mut self, node: &BVHNode, axis: &mut usize, split_pos: &mut f32, n_left: &mut usize, n_right: &mut usize, bounds_left: &mut AABB, bounds_right: &mut AABB, splitted: &mut i32) -> f32
    {
        let mut best_cost = 1e30;
        for a in 0..3
        {
            let bounds_min = node.bounds.min_bound.get_axis(a);
            let bounds_max = node.bounds.max_bound.get_axis(a);
            if bounds_min == bounds_max
            {
                continue;
            }

            // calculate cost of spatial splits
            let bin_extend = (bounds_max - bounds_min) / (BVH_BINS as f32);
            for b in 1..BVH_BINS
            {
                // calculate spatial split plane position
                let pos = bounds_min + (b as f32) * bin_extend;
                // construct left and right bounding box
                let mut left_bounds = AABB::from_empty();
                let mut right_bounds = AABB::from_empty();
                let mut left_box = AABB::from_bounds(&node.bounds.min_bound, &node.bounds.max_bound);
                let mut right_box = AABB::from_bounds(&node.bounds.min_bound, &node.bounds.max_bound);

                left_box.max_bound.set_axis(a, pos);
                right_box.min_bound.set_axis(a, pos);

                let mut left_count: usize = 0;
                let mut right_count: usize = 0;
                // loop over triangles
                for i in 0..node.tri_count
                {
                    // fetch a triangle
                    let triangle = &self.triangles[self.triangle_idx[node.left_first + i]];

                    // extend left and right side
                    let left_triangle = Clipped::new(&triangle, &left_box);

                    if left_triangle.vertices > 2
                    {
                        left_bounds.grow_aabb( &left_triangle.bounds );
                        left_count += 1;
                    }

                    let right_triangle = Clipped::new(&triangle, &right_box);

                    if right_triangle.vertices > 2
                    {
                        right_bounds.grow_aabb( &right_triangle.bounds );
                        right_count += 1;
                    }
                }
                // calculate cost for this split plane
                if left_count > 0 && right_count > 0
                {
                    let zero: f32 = 0.0;
                    let cost = zero.max(left_bounds.area()) * (left_count as f32) + zero.max(right_bounds.area()) * (right_count as f32);

                    if cost < best_cost
                    {
                        *axis = a;
                        best_cost = cost;
                        *split_pos = pos;

                        *splitted = ((right_count + left_count) as i32) - (node.tri_count as i32);
                        *n_left = left_count;
                        *n_right = right_count;
                        *bounds_left = left_bounds;
                        *bounds_right = right_bounds;
                    }

                }
            }
        }
        return best_cost;
    }

    fn calculate_node_cost(node: &BVHNode) -> f32
    {
        let e = node.bounds.max_bound - node.bounds.min_bound;
        let area = e.x * e.y + e.y * e.z + e.z * e.x;
        return (node.tri_count as f32) * area;
    }

    fn subdivide(&mut self, node_idx: usize, slack: usize)
    {
        let node = self.bvh_nodes[node_idx].clone();
        let mut obj_split_axis: usize = 0;
        let mut spatial_split_axis: usize = 0;
        let mut splitted: i32 = 0;
        let mut left_box = AABB::from_empty();
        let mut right_box = AABB::from_empty();

        let mut obj_split_pos: f32 = 0.0;
        let mut spatial_split_pos: f32 = 0.0;
        let no_split_cost: f32 = BVH::calculate_node_cost(&node);

        let mut obj_split_cost = self.find_best_object_split_plane(&node, &mut obj_split_axis, &mut obj_split_pos, &mut left_box, &mut right_box);

        let root_area = self.bvh_nodes[0].bounds.area();
        let lambda = left_box.intersection( &right_box ).area() / root_area;

        if lambda > 1e-5
        {
            let mut n_left: usize = 0;
            let mut n_right: usize = 0;
            let mut bounds_left = AABB::from_empty();
            let mut bounds_right = AABB::from_empty();

            let mut spatial_split_cost = self.find_best_spatial_split_plane(&node, &mut spatial_split_axis, &mut spatial_split_pos, &mut n_left, &mut n_right, &mut bounds_left, &mut bounds_right, &mut splitted);

            if spatial_split_cost < obj_split_cost && splitted < (slack as i32)
            {
                if spatial_split_cost >= no_split_cost
                {
                    return; // don't split, not worth it
                }
                let mut left_of_split = node.bounds.clone();
                let mut right_of_split = node.bounds.clone();
                left_of_split.max_bound.set_axis(spatial_split_axis, spatial_split_pos);
                right_of_split.min_bound.set_axis(spatial_split_axis, spatial_split_pos);

                let mut left_pos = node.left_first;
                let mut left_count: usize = 0;
                let mut right_pos = node.left_first + node.tri_count + slack;
                let mut right_count: usize = 0;

                for i in 0..node.tri_count
                {
                    let mut idx = self.triangle_idx[node.left_first + i];
                    let mut left_part = Clipped::new(&self.triangles[idx], &left_of_split);
                    let mut right_part = Clipped::new(&self.triangles[idx], &right_of_split);
                    let mut in_left = left_part.vertices >= 3;
                    let mut in_right = right_part.vertices >= 3;

                    if in_left && in_right
                    {
                        let c1 = bounds_left.union(&self.triangles[idx].bounds).area() * (n_left as f32) + bounds_right.area() * ((n_right - 1) as f32);
                        let c2 = bounds_left.area() * ((n_left - 1) as f32) + bounds_right.union( &self.triangles[idx].bounds ).area() * (n_right as f32);
                        if c1 < spatial_split_cost || c2 < spatial_split_cost
                        {
                            if c1 < c2
                            {
                                spatial_split_cost = c1;
                                n_right -= 1;
                                bounds_left.grow_aabb( &self.triangles[idx].bounds );
                                left_part.bounds = self.triangles[idx].bounds;
                                in_right = false; // undo clip
                            }
                            else {
                                spatial_split_cost = c2;
                                n_left -= 1;
                                bounds_right.grow_aabb( &self.triangles[idx].bounds );
                                right_part.bounds = self.triangles[idx].bounds;
                                in_left = false; // undo clip
                            }
                        }
                    }

                    if in_left
                    {
                        self.triangle_idx[left_pos] = idx;
                        left_pos += 1;
                        self.triangles[idx].bounds = left_part.bounds;
                        left_count += 1;
                    }
                    if in_right
                    {
                        if in_left
                        {
                            self.triangles[self.triangle_ptr] = self.triangles[idx].clone();
                            idx = self.triangle_ptr;
                            self.triangle_ptr += 1;
                        }
                        right_pos -= 1;
                        self.triangle_tmp[right_pos] = idx;
                        self.triangles[idx].bounds = right_part.bounds;
                        right_count += 1;
                    }
                }

                let slack = ((slack as i32) - splitted) as usize;
                let half_slack = slack / 2;

                let triangle_idx_adr = node.left_first + left_count + half_slack;
                let triangle_tmp_adr = right_pos;

                for i in 0..right_count
                {
                    self.triangle_idx[triangle_idx_adr + i] = self.triangle_tmp[triangle_tmp_adr + i];
                }

                let left_child_idx = self.nodes_used;
                self.nodes_used += 1;
                let right_child_idx = self.nodes_used;
                self.nodes_used += 1;

                self.bvh_nodes[left_child_idx].left_first = node.left_first;
                self.bvh_nodes[left_child_idx].tri_count = left_count;
                self.bvh_nodes[right_child_idx].left_first = node.left_first + left_count + half_slack;
                self.bvh_nodes[right_child_idx].tri_count = right_count;

                self.bvh_nodes[node_idx].left_first = left_child_idx;
                self.bvh_nodes[node_idx].tri_count = 0;

                self.bvh_nodes[left_child_idx].bounds = bounds_left;
                self.bvh_nodes[right_child_idx].bounds = bounds_right;
                self.spatial_splits += 1;
                self.subdivide(left_child_idx, half_slack);
                self.subdivide(right_child_idx, half_slack);
                return;
            }
        }
        if obj_split_cost >= no_split_cost
        {
            return;
        }

        let mut i = node.left_first;
        let mut j: i32 = ((i + node.tri_count) as i32) - 1;

        while (i as i32) <= j
        {
            let jx = j as usize;
            let triangle = &self.triangles[self.triangle_idx[i]];
            if triangle.bounds.center(obj_split_axis) < obj_split_pos
            {
                i += 1;
            }
            else {
                if i != jx
                {
                    self.triangle_idx.swap(i, jx);
                }
                j -= 1;
            }
        }

        let half_slack = slack / 2;
        let left_count = i - node.left_first;
        if left_count == 0 || left_count == node.tri_count
        {
            return;
        }

        if half_slack > 0
        {
            let lb = i;
            let ub = i + (node.tri_count - left_count);
            self.triangle_idx.copy_within((lb..ub), i + half_slack);
        }

        let left_child_idx = self.nodes_used;
        self.nodes_used += 1;
        let right_child_idx = self.nodes_used;
        self.nodes_used += 1;
        self.bvh_nodes[left_child_idx].left_first = node.left_first;
        self.bvh_nodes[left_child_idx].tri_count = left_count;
        self.bvh_nodes[right_child_idx].left_first = i + half_slack;
        self.bvh_nodes[right_child_idx].tri_count = node.tri_count - left_count;

        self.bvh_nodes[node_idx].left_first = left_child_idx;
        self.bvh_nodes[node_idx].tri_count = 0;

        self.bvh_nodes[left_child_idx].bounds = left_box;
        self.bvh_nodes[right_child_idx].bounds = right_box;
        self.subdivide(left_child_idx, half_slack);
        self.subdivide(right_child_idx, half_slack);
    }
}

impl RayHittableObject for BVH
{
    fn intersect(&self, ray: &mut Ray) {
        let O = transform_position( &ray.origin, &self.inv_t );
        let D = transform_vector( &ray.direction, &self.inv_t );

        let mut ray_t = Ray::directed_distance(O, D, ray.t);
        ray_t.obj_idx = ray.obj_idx;

        let mut node = &self.bvh_nodes[self.root_node_idx];
        let mut stack = [BVHNode{
            bounds: AABB::from_empty(),
            tri_count: 0,
            left_first: 0
        }; 64];

        let mut stack_ptr: usize = 0;
        loop
        {
            if node.is_leaf()
            {
                for i in 0..node.tri_count
                {
                    if self.intersect_triangle(&mut ray_t, self.triangle_idx[node.left_first + i])
                    {
                        ray.obj_idx = self.obj_idx;
                        ray.obj_type = RayHittableObjectType::Bvh;
                    }
                }

                if stack_ptr == 0
                {
                    break;
                }
                else
                {
                    stack_ptr -= 1;
                    node = &stack[stack_ptr];
                }
                continue;
            }
            let mut child1 = &self.bvh_nodes[node.left_first];
            let mut child2 = &self.bvh_nodes[node.left_first + 1];
            let mut dist1 = self.intersect_aabb( &ray_t, &child1.bounds );
            let mut dist2 = self.intersect_aabb( &ray_t, &child2.bounds );
            if dist1 > dist2
            {
                let tmp = dist1;
                dist1 = dist2;
                dist2 = tmp;
                let tmp = child1;
                child1 = child2;
                child2 = tmp;
            }
            if dist1 == 1e30
            {
                if stack_ptr == 0
                {
                    break;
                }
                else
                {
                    stack_ptr -= 1;
                    node = &stack[stack_ptr];
                }
                continue;
            }
            node = child1;
            if dist2 != 1e30
            {
                stack[stack_ptr] = *child2;
                stack_ptr += 1;
            }
        }
        ray.t = ray_t.t;
        ray.sub_obj_idx = ray_t.sub_obj_idx;
    }

    fn get_normal(&self, i: &Float3) -> Float3 {
        return self.normals[0];
    }

    fn get_uv(&self, i: &Float3) -> Float2 {
        return Float2::zero();
    }

    fn is_occluded(&self, ray: &Ray) -> bool {
        let mut shadow = ray.clone();
        shadow.t = 1e34;
        shadow.obj_idx = i32::MAX;
        self.intersect(&mut shadow);
        return shadow.obj_idx == i32::MAX;
    }
}


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

    fn get_normal(&self, i: &Float3) -> Float3 {
        // needs to interpolate between normals
        self.triangle_normals[0]
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