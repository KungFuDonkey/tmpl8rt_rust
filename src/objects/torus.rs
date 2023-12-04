use crate::math::*;
use crate::ray::*;

pub struct Torus
{
    pub rt2: f32,
    pub rc2: f32,
    pub r2: f32,
    pub obj_idx: usize,
    pub mat_idx: usize,
    pub t: Mat4,
    pub inv_t: Mat4
}

impl Torus
{
    pub fn new(obj_idx: usize, mat_idx: usize, a: f32, b: f32) -> Self
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

    fn get_normal(&self, _: &Ray, i: &Float3) -> Float3 {
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