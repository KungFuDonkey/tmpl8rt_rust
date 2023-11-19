use std::f32::consts::PI;
use crate::math::*;

#[derive(Debug, Clone, Copy)]
pub enum RayHittableObjectType
{
    None,
    Sphere,
    Plane,
    Torus
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

    pub fn new() -> Self
    {
        let r: Ray;
        unsafe
        {
            r = Ray {
                origin: Float3::zero(),
                direction: Float3::zero(),
                r_direction: Float3::zero(),
                t: 1.0e34,
                obj_idx: -1,
                obj_type: RayHittableObjectType::None,
                inside: false
            };
        }
        return r;
    }

    pub fn directed(origin: Float3, direction: Float3) -> Self
    {
        let r: Ray;
        unsafe
            {
                r = Ray {
                    origin: Float3::from_xyz(origin.x, origin.y, origin.z),
                    direction: Float3::from_xyz(direction.x, direction.y, direction.z),
                    r_direction: Float3::from_xyz(1.0 / direction.x, 1.0 / direction.y, 1.0 / direction.z),
                    t: 1.0e34,
                    obj_idx: -1,
                    obj_type: RayHittableObjectType::None,
                    inside: false
                };
            }
        return r;
    }

    pub fn directed_distance(origin: Float3, direction: Float3, distance: f32) -> Self
    {
        let r: Ray;
        unsafe
            {
                r =  Ray {
                    origin: Float3::from_xyz(origin.x, origin.y, origin.z),
                    direction: Float3::from_xyz(direction.x, direction.y, direction.z),
                    r_direction: Float3::from_xyz(1.0 / direction.x, 1.0 / direction.y, 1.0 / direction.z),
                    t: distance,
                    obj_idx: -1,
                    obj_type: RayHittableObjectType::None,
                    inside: false
                };
            }
        return r;
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

    fn get_albedo(&self, i: &Float3) -> Float3;
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
    pub obj_idx: i32
}

impl RayHittableObject for Sphere
{
    fn intersect(&self, ray: &mut Ray) {
        let oc = ray.origin - self.position;
        let b = dot_f3(&oc, &ray.direction);
        let c = dot_f3(&oc, &oc) - self.radius2;
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

    fn get_albedo(&self, i: &Float3) -> Float3
    {
        return Float3::from_a(0.93)
    }
}

impl Sphere
{
    fn new(obj_idx: i32, position: Float3, radius: f32) -> Self
    {
        Sphere
        {
            obj_idx,
            position,
            radius,
            radius2: radius * radius,
            inv_radius: 1.0 / radius
        }
    }
}

pub struct Plane
{
    pub obj_idx: i32,
    pub normal: Float3,
    pub distance: f32
}

impl RayHittableObject for Plane
{
    fn intersect(&self, ray: &mut Ray) {
        let t = -(dot_f3(&ray.origin, &self.normal) + self.distance) / (dot_f3(&ray.direction, &self.normal));
        if t < ray.t && t > 0.0
        {
            ray.t = t;
            ray.obj_idx = self.obj_idx;
            ray.obj_type = RayHittableObjectType::Plane;
        }
    }

    fn get_normal(&self, i: &Float3) -> Float3 {
        return self.normal;
    }

    fn get_albedo(&self, i: &Float3) -> Float3 {
        return Float3::from_xyz(1.0, 0.0, 0.0);
    }
}

impl Plane
{
    pub fn new(obj_idx: i32, normal: Float3, distance: f32) -> Self
    {
        Plane{
            obj_idx,
            normal,
            distance
        }
    }
}

pub struct Cube
{
    pub obj_idx: i32,
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

        tmin = fmaxf(tmin, tymin);
        tmax = fminf(tmax, tymax);
        let tzmin = (self.b[signz].z - origin.z) * rdz;
        let tzmax = (self.b[1 - signz].z - origin.z) * rdz;

        if tmin > tzmax || tzmin > tmax
        {
            return;
        }

        tmin = fmaxf(tmin, tzmin);
        tmax = fminf(tmax, tzmax);
        if tmin > 0.0
        {
            if tmin < ray.t
            {
                ray.t = tmin;
                ray.obj_idx = self.obj_idx;
                ray.obj_type = RayHittableObjectType::Plane;
            }
        }
        else if tmax > 0.0
        {
            if tmax < ray.t
            {
                ray.t = tmax;
                ray.obj_idx = self.obj_idx;
                ray.obj_type = RayHittableObjectType::Plane;
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
            min_dist = d5;
            n = Float3::from_xyz(0.0, 0.0, 1.0);
        }

        return transform_vector(&n, &self.m);
    }

    fn get_albedo(&self, i: &Float3) -> Float3 {
        return Float3::from_a(1.0);
    }
}

impl Cube
{
    pub fn new(obj_idx: i32, pos: Float3, size: Float3) -> Self
    {
        Cube{
            obj_idx,
            b:[
                pos - size * 0.5,
                pos + size * 0.5
            ],
            m: Mat4::identity_matrix(),
            inv_m: Mat4::identity_matrix().inverted()
        }
    }
}

pub struct Torus
{
    pub rt2: f32,
    pub rc2: f32,
    pub r2: f32,
    pub obj_idx: i32,
    pub t: Mat4,
    pub inv_t: Mat4
}

impl Torus
{
    pub fn new(obj_idx: i32, a: f32, b: f32) -> Self
    {
        Torus
        {
            obj_idx,
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
        let O: Float3 = transform_position( &ray.origin, &self.inv_t );
        let D: Float3 = transform_vector(&ray.direction, &self.inv_t);

        let r2 = self.r2 as f64;
        let rt2 = self.rt2 as f64;
        let rc2 = self.rc2 as f64;

        // extension rays need double precision for the quadratic solver!
        let mut po = 1.0;
        let mut m = dot_f3(&O, &O) as f64;
        let mut k3 = dot_f3(&O, &D) as f64;
        let mut k32 = k3 * k3;

        // bounding sphere test
        let v = k32 - m + r2;
        if v < 0.0
        {
            return;
        }

        // setup torus intersection
        let mut k = (m - rt2 - rc2) * 0.5;
        let mut k2 = k32 + rc2 * ((D.z * D.z) as f64) + k;
        let mut k1 = k * k3 + rc2 * ((O.z * D.z) as f64);
        let mut k0 = k * k + rc2 * ((O.z * O.z) as f64) - rc2 * rt2;
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
        let mut Q = c2 * c2 + c0;
        let mut R = 3.0 * c0 * c2 - c2 * c2 * c2 - c1 * c1;
        let mut h = R * R - Q * Q * Q;
        let mut z: f64;
        if h < 0.0
        {
            let mut sQ = Q.sqrt();
            z = 2.0 * sQ * ((R / (sQ * Q)).acos() * 0.33333333333).cos();
        }
        else
        {
            let mut sQ = Torus::cbrt_fast(h.sqrt() + R.abs());
            z = ( sQ + Q / sQ ).abs().copysign(R);
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
        let mut ft = t as f32;
        if ft > 0.0 && ft < ray.t
        {
            ray.t = ft;
            ray.obj_idx = self.obj_idx;
            ray.obj_type = RayHittableObjectType::Torus;
        }
    }

    fn get_normal(&self, i: &Float3) -> Float3 {
        let l = transform_position(&i, &self.inv_t);
        let x = -Float3::from_xyz(1.0, 1.0, -1.0) * self.rc2 + l * (dot_f3(&l, &l) - self.rt2);
        let n = normalize_f3(&x);
        return transform_vector(&n, &self.t);
    }

    fn get_albedo(&self, i: &Float3) -> Float3 {
        return Float3::from_a(1.0);
    }
}


pub struct Scene
{
    spheres: Vec<Sphere>,
    planes: Vec<Plane>,
    cubes: Vec<Cube>,
    tori: Vec<Torus>,
    animation_time: f32,
}

impl Scene
{
    pub fn new() -> Self
    {
        let mut torus = Torus::new(0, 0.8, 0.25);
        let translation = Float3::from_xyz(-0.25, 0.0, 2.0);
        torus.t = Mat4::translate(&translation) * Mat4::rotate_x(PI / 4.0);
        torus.inv_t = torus.t.inverted();

        Scene{
            spheres: vec![
                Sphere::new(0, Float3::from_a(0.0), 0.6),
                Sphere::new(1, Float3::from_xyz( 0.0, 2.5, -3.07 ), 8.0)
            ],
            planes: vec![
                Plane::new(0, Float3::from_xyz(1.0, 0.0, 0.0), 3.0),
                Plane::new(1, Float3::from_xyz(-1.0, 0.0, 0.0), 2.99),
                Plane::new(2, Float3::from_xyz(0.0, 1.0, 0.0), 1.0),
                Plane::new(3, Float3::from_xyz(0.0, -1.0, 0.0), 2.0),
                Plane::new(4, Float3::from_xyz(0.0, 0.0, 1.0), 3.0),
                Plane::new(5, Float3::from_xyz(0.0, 0.0, -1.0), 3.99),
            ],
            cubes: vec![
                Cube::new(0, Float3::zero(), Float3::from_a(1.15))
            ],
            tori: vec![
                torus
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
    }


    pub fn get_normal(&self, obj_idx: i32, obj_type: &RayHittableObjectType, i: &Float3, wo: &Float3) -> Float3
    {
        if obj_idx == -1
        {
            println!("ERROR: obj_idx not set or no object was hit");
            return Float3::zero();
        }

        match obj_type
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
            RayHittableObjectType::Torus =>
            {
                self.tori[obj_idx as usize].get_normal(i)
            }
        }
    }

    pub fn get_albedo(&self, obj_idx: i32, obj_type: &RayHittableObjectType, i: &Float3) -> Float3
    {
        match obj_type
        {
            RayHittableObjectType::None =>
            {
                Float3::zero()
            },
            RayHittableObjectType::Sphere =>
            {
                self.spheres[obj_idx as usize].get_albedo(i)
            }
            RayHittableObjectType::Plane =>
            {
                self.planes[obj_idx as usize].get_albedo(i)
            }
            RayHittableObjectType::Torus =>
            {
                self.tori[obj_idx as usize].get_albedo(i)
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