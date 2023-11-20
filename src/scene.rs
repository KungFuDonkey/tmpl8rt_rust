use std::f32::consts::PI;
use crate::material::{BlueMaterial, CheckerBoardMaterial, LinearColorMaterial, LogoMaterial, Material, RedMaterial};
use crate::math::*;

#[derive(Debug, Clone, Copy)]
pub enum RayHittableObjectType
{
    None,
    Sphere,
    Plane,
    Cube,
    Torus,
    Quad
}

#[repr(C, align(64))]
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    pub origin: Float3,
    pub direction: Float3,
    pub  r_direction: Float3,
    pub t: f32,
    pub obj_idx: i32,
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
}

impl RayHittableObject for Quad
{
    fn intersect(&self, ray: &mut Ray) {
        let Oy: f32 = self.inv_t.cell[4] * ray.origin.x + self.inv_t.cell[5] * ray.origin.y + self.inv_t.cell[6] * ray.origin.z + self.inv_t.cell[7];
        let Dy: f32 = self.inv_t.cell[4] * ray.direction.x + self.inv_t.cell[5] * ray.direction.y + self.inv_t.cell[6] * ray.direction.z;
        let t = Oy / -Dy;
        if t < ray.t && t > 0.0
        {
            let Ox: f32 = self.inv_t.cell[0] * ray.origin.x + self.inv_t.cell[1] * ray.origin.y + self.inv_t.cell[2] * ray.origin.z + self.inv_t.cell[3];
            let Oz: f32 = self.inv_t.cell[8] * ray.origin.x + self.inv_t.cell[9] * ray.origin.y + self.inv_t.cell[10] * ray.origin.z + self.inv_t.cell[11];
            let Dx: f32 = self.inv_t.cell[0] * ray.direction.x + self.inv_t.cell[1] * ray.direction.y + self.inv_t.cell[2] * ray.direction.z;
            let Dz: f32 = self.inv_t.cell[8] * ray.direction.x + self.inv_t.cell[9] * ray.direction.y + self.inv_t.cell[10] * ray.direction.z;
            let Ix = Ox + t * Dx;
            let Iz = Oz + t * Dz;
            if Ix > -self.size && Ix < self.size && Iz > -self.size && Iz < self.size
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
}


pub struct Scene
{
    spheres: Vec<Sphere>,
    planes: Vec<Plane>,
    cubes: Vec<Cube>,
    tori: Vec<Torus>,
    quads: Vec<Quad>,
    materials: Vec<Box<dyn Material + Sync>>,
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

        Scene{
            spheres: vec![
                Sphere::new(0, 0, Float3::from_a(0.0), 0.6),
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
                Cube::new(0, 0, Float3::zero(), Float3::from_a(1.15))
            ],
            tori: vec![
                torus
            ],
            quads: vec![
                Quad::new(0, 6, 0.5, &Mat4::translate(&Float3::from_xyz(-1.0, 1.5, -1.0))),
                Quad::new(1, 6, 0.5, &Mat4::translate(&Float3::from_xyz(1.0, 1.5, -1.0))),
                Quad::new(2, 6, 0.5, &Mat4::translate(&Float3::from_xyz(1.0, 1.5, 1.0))),
                Quad::new(3, 6, 0.5, &Mat4::translate(&Float3::from_xyz(-1.0, 1.5, 1.0))),
            ],
            materials: vec![
                Box::new(LinearColorMaterial::new(Float3::from_a(1.0))),
                Box::new(LinearColorMaterial::new(Float3::from_a(0.93))),
                Box::new(CheckerBoardMaterial::new()),
                Box::new(LogoMaterial::new()),
                Box::new(RedMaterial::new()),
                Box::new(BlueMaterial::new()),
                Box::new(LinearColorMaterial::new(Float3::from_xyz(1.0, 0.0, 0.0))),
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
    }


    pub fn get_normal(&self, obj_idx: i32, obj_type: &RayHittableObjectType, i: &Float3, wo: &Float3) -> Float3
    {
        if obj_idx == -1
        {
            println!("ERROR: obj_idx not set or no object was hit");
            return Float3::zero();
        }

        let normal = match obj_type
        {
            RayHittableObjectType::None =>
            {
                println!("ERROR: tried to get normal of non object");
                Float3::zero()
            },
            RayHittableObjectType::Sphere =>
            {
                self.spheres[obj_idx as usize].get_normal(i)
            },
            RayHittableObjectType::Plane =>
            {
                self.planes[obj_idx as usize].get_normal(i)
            },
            RayHittableObjectType::Cube =>
            {
                self.cubes[obj_idx as usize].get_normal(i)
            },
            RayHittableObjectType::Torus =>
            {
                self.tori[obj_idx as usize].get_normal(i)
            },
            RayHittableObjectType::Quad =>
            {
                self.quads[obj_idx as usize].get_normal(i)
            }
        };

        if dot(&normal, &wo) > 0.0
        {
            return -normal;
        }
        return normal;
    }

    pub fn get_albedo(&self, obj_idx: i32, obj_type: &RayHittableObjectType, i: &Float3) -> Float3
    {
        let (uv, mat_idx) = match obj_type
        {
            RayHittableObjectType::None =>
            {
                (Float2::zero(), -1)
            },
            RayHittableObjectType::Sphere =>
            {
                let s = &self.spheres[obj_idx as usize];
                (s.get_uv(i), s.mat_idx)
            },
            RayHittableObjectType::Plane =>
            {
                let p = &self.planes[obj_idx as usize];
                (p.get_uv(i), p.mat_idx)
            },
            RayHittableObjectType::Cube =>
            {
                let c = &self.cubes[obj_idx as usize];
                (c.get_uv(i), c.mat_idx)
            },
            RayHittableObjectType::Quad =>
            {
                let q = &self.quads[obj_idx as usize];
                (q.get_uv(i), q.mat_idx)
            },
            RayHittableObjectType::Torus =>
            {
                let t = &self.tori[obj_idx as usize];
                (t.get_uv(i), t.mat_idx)
            }
        };

        if mat_idx == -1
        {
            return Float3::zero();
        }

        return self.materials[mat_idx as usize].get_color(&uv);
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