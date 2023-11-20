#![allow(dead_code)]

use std::cmp::min;
use std::cmp::max;
use std::ops::*;

#[repr(C, align(8))]
#[derive(Debug, Clone, Copy)]
pub struct Int2 {
    pub x: i32,
    pub y: i32
}

impl Int2 {
    pub fn zero() -> Self {
        Int2 { x: 0, y: 0 }
    }

    pub fn from_xy(x: i32, y: i32) -> Self
    {
        Int2 { x, y }
    }

    pub fn from_a(a: i32) -> Self {
        Int2 { x:a, y:a }
    }
}

#[repr(C, align(8))]
#[derive(Debug, Clone, Copy)]
pub struct Uint2 {
    pub x: u32,
    pub y: u32
}

impl Uint2 {
    pub fn zero() -> Self {
        Uint2 { x: 0, y: 0 }
    }

    pub fn from_xy(x: u32, y: u32) -> Self
    {
        Uint2 { x, y }
    }

    pub fn from_a(a: u32) -> Self {
        Uint2 { x:a, y:a }
    }
}

#[repr(C, align(8))]
#[derive(Debug, Clone, Copy)]
pub struct Float2 {
    pub x: f32,
    pub y: f32
}

impl Float2 {
    pub fn zero() -> Self {
        Float2 { x: 0.0, y: 0.0 }
    }

    pub fn from_xy(x: f32, y: f32) -> Self
    {
        Float2 { x, y }
    }

    pub fn from_a(a: f32) -> Self {
        Float2 { x:a, y:a }
    }
}

#[repr(C, align(16))]
#[derive(Debug, Clone, Copy)]
pub struct Int3 {
    pub x: i32,
    pub y: i32,
    pub z: i32
}

impl Int3 {
    pub fn zero() -> Self {
        Int3 { x: 0, y: 0, z: 0 }
    }

    pub fn from_xyz(x: i32, y: i32, z: i32) -> Self
    {
        Int3 { x, y, z }
    }

    pub fn from_a(a: i32) -> Self {
        Int3 { x: a, y: a, z: a }
    }
}

#[repr(C, align(16))]
#[derive(Debug, Clone, Copy)]
pub struct Uint3 {
    pub x: u32,
    pub y: u32,
    pub z: u32
}

impl Uint3 {
    pub fn zero() -> Self {
        Uint3 { x: 0, y: 0, z: 0 }
    }

    pub fn from_xyz(x: u32, y: u32, z: u32) -> Self
    {
        Uint3 { x, y, z }
    }

    pub fn from_a(a: u32) -> Self {
        Uint3 { x: a, y: a, z: a }
    }

}

#[repr(C, align(16))]
#[derive(Debug, Clone, Copy)]
pub struct Float3 {
    pub x: f32,
    pub y: f32,
    pub z: f32
}

impl Float3 {
    pub fn zero() -> Self {
        Float3 { x: 0.0, y: 0.0, z: 0.0 }
    }

    pub fn from_xyz(x: f32, y: f32, z: f32) -> Self
    {
        Float3 { x, y, z }
    }

    pub fn from_a(a: f32) -> Self {
        Float3 { x: a, y: a, z: a }
    }
}

#[repr(C, align(16))]
#[derive(Debug, Clone, Copy)]
pub struct Int4 {
    pub x: i32,
    pub y: i32,
    pub z: i32,
    pub w: i32
}

impl Int4 {
    pub fn zero() -> Self {
        Int4 { x: 0, y: 0, z: 0, w: 0 }
    }

    pub fn from_xyzw(x: i32, y: i32, z: i32, w: i32) -> Self
    {
        Int4 { x, y, z, w }
    }

    pub fn from_a(a: i32) -> Self {
        Int4 { x: a, y: a, z: a, w: a }
    }
}

#[repr(C, align(8))]
#[derive(Debug, Clone, Copy)]
pub struct Uint4 {
    pub x: u32,
    pub y: u32,
    pub z: u32,
    pub w: u32
}

impl Uint4 {
    pub fn zero() -> Self {
        Uint4 { x: 0, y: 0, z: 0, w: 0 }
    }

    pub fn from_xyzw(x: u32, y: u32, z: u32, w: u32) -> Self
    {
        Uint4 { x, y, z, w }
    }

    pub fn from_a(a: u32) -> Self {
        Uint4 { x: a, y: a, z: a, w: a }
    }
}

#[repr(C, align(8))]
#[derive(Debug, Clone, Copy)]
pub struct Float4 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub w: f32
}

impl Float4 {
    pub fn zero() -> Self {
        Float4 { x: 0.0, y: 0.0, z: 0.0, w: 0.0 }
    }

    pub fn from_xyzw(x: f32, y: f32, z: f32, w: f32) -> Self
    {
        Float4 { x, y, z, w }
    }

    pub fn from_a(a: f32) -> Self {
        Float4 { x: a, y: a, z: a, w: a }
    }
}

// swap
pub fn swap<T>(a: &mut T, b: &mut T) {
    core::mem::swap(a,b);
}

// random numbers
static mut SEED: u32 = 0x12345678;

fn wang_hash(s: u32) -> u32
{
    let mut v = (s ^ 61) ^ (s >> 16);
    v *= 9;
    v = v ^ (v >> 4);
    v *= 0x27d4eb2d;
    v = v ^ (v >> 15);
    return v;
}

pub fn init_seed(seed_base: u32) -> u32
{
    return wang_hash((seed_base + 1) * 17)
}

pub fn random_uint() -> u32
{
    let v: u32 = 0;
    unsafe {
        SEED ^= SEED << 13;
        SEED ^= SEED >> 17;
        SEED ^= SEED << 5;
    }
    return v;
}

pub fn random_uint_s(seed: &mut u32) -> u32
{
    *seed ^= *seed << 13;
    *seed ^= *seed >> 17;
    *seed ^= *seed << 5;
    return *seed;
}

pub fn random_float() -> f32
{
    (random_uint() as f32) * 2.3283064365387e-10
}

pub fn random_float_s(seed: &mut u32) -> f32
{
    (random_uint_s(seed) as f32) * 2.3283064365387e-10
}

pub fn rand(range: f32) -> f32
{
    random_float() * range
}

static NUMX: i32 = 512;
static NUMY: i32 = 512;
static NUMOCTAVES: i32 = 7;
static PRIMEINDEX: i32 = 0;
static PERSISTENCE: f32 = 0.5;

static PRIMES: [[i32; 3]; 10] = [
    [ 995615039, 600173719, 701464987 ],
    [ 831731269, 162318869, 136250887 ],
    [ 174329291, 946737083, 245679977 ],
    [ 362489573, 795918041, 350777237 ],
    [ 457025711, 880830799, 909678923 ],
    [ 787070341, 177340217, 593320781 ],
    [ 405493717, 291031019, 391950901 ],
    [ 458904767, 676625681, 424452397 ],
    [ 531736441, 939683957, 810651871 ],
    [ 997169939, 842027887, 423882827 ]
];

pub fn noise(i: usize, x: i32, y: i32) -> f32
{
    let mut n = x + y * 57;
    n = (n << 13) ^ n;
    let a = PRIMES[i][0];
    let b = PRIMES[i][1];
    let c = PRIMES[i][2];
    let t: i32 = (n * (n * n * a + b) + c) & 0x7fffffff;
    return 1.0 - (t as f32) / 1073741824.0;
}

pub fn expf( a: &Float3) -> Float3
{
    Float3::from_xyz(fast_math::exp(a.x), fast_math::exp(a.y), fast_math::exp(a.z))
}

macro_rules! impl_default_operator_2
{
    ($type:ty, $trait:ty, $trait_fn:ident) => {
        impl $trait for $type
        {
            type Output = $type;
            fn $trait_fn(self, rhs: Self) -> Self::Output
            {
                <$type>::from_xy(self.x.$trait_fn(rhs.x), self.y.$trait_fn(rhs.y))
            }
        }
    }
}

macro_rules! impl_default_assign_operator_2
{
    ($type:ty, $trait:ty, $trait_fn:ident) => {
        impl $trait for $type
        {
            fn $trait_fn(&mut self, rhs: Self)
            {
                self.x.$trait_fn(rhs.x);
                self.y.$trait_fn(rhs.y);
            }
        }
    }
}

macro_rules! impl_single_operator_2 {
    ($type:ty, $subtype:ty, $trait:ty, $trait_t:ty, $trait_fn:ident) => {
        impl $trait for $type
        {
            type Output = $type;
            fn $trait_fn(self, rhs: $trait_t) -> Self::Output
            {
                <$type>::from_xy(self.x.$trait_fn(rhs as $subtype), self.y.$trait_fn(rhs as $subtype))
            }
        }
    };
}

macro_rules! impl_single_assign_operator_2 {
    ($type:ty, $subtype:ty, $trait:ty, $trait_t:ty, $trait_fn:ident) => {
        impl $trait for $type
        {
            fn $trait_fn(&mut self, rhs: $trait_t)
            {
                self.x.$trait_fn(rhs as $subtype);
                self.y.$trait_fn(rhs as $subtype);
            }
        }
    };
}

macro_rules! impl_operators_2 {
    ($type:ty, $subtype:ty) => {

        impl_default_operator_2!($type, Sub, sub);
        impl_default_assign_operator_2!($type, SubAssign, sub_assign);
        impl_single_operator_2!($type, $subtype, Sub<i32>, i32, sub);
        impl_single_assign_operator_2!($type, $subtype, SubAssign<i32>, i32, sub_assign);
        impl_single_operator_2!($type, $subtype, Sub<u32>, u32, sub);
        impl_single_assign_operator_2!($type, $subtype, SubAssign<u32>, u32, sub_assign);
        impl_single_operator_2!($type, $subtype, Sub<f32>, f32, sub);
        impl_single_assign_operator_2!($type, $subtype, SubAssign<f32>, f32, sub_assign);

        impl_default_operator_2!($type, Add, add);
        impl_default_assign_operator_2!($type, AddAssign, add_assign);
        impl_single_operator_2!($type, $subtype, Add<i32>, i32, add);
        impl_single_assign_operator_2!($type, $subtype, AddAssign<i32>, i32, add_assign);
        impl_single_operator_2!($type, $subtype, Add<u32>, u32, add);
        impl_single_assign_operator_2!($type, $subtype, AddAssign<u32>, u32, add_assign);
        impl_single_operator_2!($type, $subtype, Add<f32>, f32, add);
        impl_single_assign_operator_2!($type, $subtype, AddAssign<f32>, f32, add_assign);

        impl_default_operator_2!($type, Mul, mul);
        impl_default_assign_operator_2!($type, MulAssign, mul_assign);
        impl_single_operator_2!($type, $subtype, Mul<i32>, i32, mul);
        impl_single_assign_operator_2!($type, $subtype, MulAssign<i32>, i32, mul_assign);
        impl_single_operator_2!($type, $subtype, Mul<u32>, u32, mul);
        impl_single_assign_operator_2!($type, $subtype, MulAssign<u32>, u32, mul_assign);
        impl_single_operator_2!($type, $subtype, Mul<f32>, f32, mul);
        impl_single_assign_operator_2!($type, $subtype, MulAssign<f32>, f32, mul_assign);
    };
}

macro_rules! impl_additional_operators_2 {
    ($type:ty) => {
        impl Neg for $type
        {
            type Output = $type;

            fn neg(self) -> Self::Output {
                <$type>::from_xy(-self.x, -self.y)
            }
        }
    };
}

impl_operators_2!(Int2, i32);
impl_additional_operators_2!(Int2);
impl_operators_2!(Float2, f32);
impl_additional_operators_2!(Float2);
impl_operators_2!(Uint2, u32);
impl_default_operator_2!(Float2, Div, div);
impl_single_operator_2!(Float2, f32, Div<f32>, f32, div);
impl_single_operator_2!(Int2, i32, Shl<i32>, i32, shl);
impl_single_operator_2!(Int2, i32, Shr<i32>, i32, shr);


macro_rules! impl_default_operator_3
{
    ($type:ty, $trait:ty, $trait_fn:ident) => {
        impl $trait for $type
        {
            type Output = $type;
            fn $trait_fn(self, rhs: Self) -> Self::Output
            {
                <$type>::from_xyz(self.x.$trait_fn(rhs.x), self.y.$trait_fn(rhs.y), self.z.$trait_fn(rhs.z))
            }
        }
    }
}

macro_rules! impl_default_assign_operator_3
{
    ($type:ty, $trait:ty, $trait_fn:ident) => {
        impl $trait for $type
        {
            fn $trait_fn(&mut self, rhs: Self)
            {
                self.x.$trait_fn(rhs.x);
                self.y.$trait_fn(rhs.y);
                self.z.$trait_fn(rhs.z);
            }
        }
    }
}

macro_rules! impl_single_operator_3 {
    ($type:ty, $subtype:ty, $trait:ty, $trait_t:ty, $trait_fn:ident) => {
        impl $trait for $type
        {
            type Output = $type;
            fn $trait_fn(self, rhs: $trait_t) -> Self::Output
            {
                <$type>::from_xyz(self.x.$trait_fn(rhs as $subtype), self.y.$trait_fn(rhs as $subtype), self.z.$trait_fn(rhs as $subtype))
            }
        }
    };
}

macro_rules! impl_single_assign_operator_3 {
    ($type:ty, $subtype:ty, $trait:ty, $trait_t:ty, $trait_fn:ident) => {
        impl $trait for $type
        {
            fn $trait_fn(&mut self, rhs: $trait_t)
            {
                self.x.$trait_fn(rhs as $subtype);
                self.y.$trait_fn(rhs as $subtype);
                self.z.$trait_fn(rhs as $subtype);
            }
        }
    };
}

macro_rules! impl_operators_3 {
    ($type:ty, $subtype:ty) => {

        impl_default_operator_3!($type, Sub, sub);
        impl_default_assign_operator_3!($type, SubAssign, sub_assign);
        impl_single_operator_3!($type, $subtype, Sub<i32>, i32, sub);
        impl_single_assign_operator_3!($type, $subtype, SubAssign<i32>, i32, sub_assign);
        impl_single_operator_3!($type, $subtype, Sub<u32>, u32, sub);
        impl_single_assign_operator_3!($type, $subtype, SubAssign<u32>, u32, sub_assign);
        impl_single_operator_3!($type, $subtype, Sub<f32>, f32, sub);
        impl_single_assign_operator_3!($type, $subtype, SubAssign<f32>, f32, sub_assign);

        impl_default_operator_3!($type, Add, add);
        impl_default_assign_operator_3!($type, AddAssign, add_assign);
        impl_single_operator_3!($type, $subtype, Add<i32>, i32, add);
        impl_single_assign_operator_3!($type, $subtype, AddAssign<i32>, i32, add_assign);
        impl_single_operator_3!($type, $subtype, Add<u32>, u32, add);
        impl_single_assign_operator_3!($type, $subtype, AddAssign<u32>, u32, add_assign);
        impl_single_operator_3!($type, $subtype, Add<f32>, f32, add);
        impl_single_assign_operator_3!($type, $subtype, AddAssign<f32>, f32, add_assign);

        impl_default_operator_3!($type, Mul, mul);
        impl_default_assign_operator_3!($type, MulAssign, mul_assign);
        impl_single_operator_3!($type, $subtype, Mul<i32>, i32, mul);
        impl_single_assign_operator_3!($type, $subtype, MulAssign<i32>, i32, mul_assign);
        impl_single_operator_3!($type, $subtype, Mul<u32>, u32, mul);
        impl_single_assign_operator_3!($type, $subtype, MulAssign<u32>, u32, mul_assign);
        impl_single_operator_3!($type, $subtype, Mul<f32>, f32, mul);
        impl_single_assign_operator_3!($type, $subtype, MulAssign<f32>, f32, mul_assign);
    };
}

macro_rules! impl_additional_operators_3 {
    ($type:ty) => {
        impl Neg for $type
        {
            type Output = $type;

            fn neg(self) -> Self::Output {
                <$type>::from_xyz(-self.x, -self.y, -self.z)
            }
        }
    };
}

impl_operators_3!(Int3, i32);
impl_additional_operators_3!(Int3);
impl_operators_3!(Float3, f32);
impl_additional_operators_3!(Float3);
impl_operators_3!(Uint3, u32);
impl_default_operator_3!(Float3, Div, div);
impl_single_operator_3!(Float3, f32, Div<f32>, f32, div);
impl_single_operator_3!(Int3, i32, Shl<i32>, i32, shl);
impl_single_operator_3!(Int3, i32, Shr<i32>, i32, shr);

macro_rules! impl_default_operator_4
{
    ($type:ty, $trait:ty, $trait_fn:ident) => {
        impl $trait for $type
        {
            type Output = $type;
            fn $trait_fn(self, rhs: Self) -> Self::Output
            {
                <$type>::from_xyzw(self.x.$trait_fn(rhs.x), self.y.$trait_fn(rhs.y), self.z.$trait_fn(rhs.z), self.w.$trait_fn(rhs.w))
            }
        }
    }
}

macro_rules! impl_default_assign_operator_4
{
    ($type:ty, $trait:ty, $trait_fn:ident) => {
        impl $trait for $type
        {
            fn $trait_fn(&mut self, rhs: Self)
            {
                self.x.$trait_fn(rhs.x);
                self.y.$trait_fn(rhs.y);
                self.z.$trait_fn(rhs.z);
                self.w.$trait_fn(rhs.w);
            }
        }
    }
}

macro_rules! impl_single_operator_4 {
    ($type:ty, $subtype:ty, $trait:ty, $trait_t:ty, $trait_fn:ident) => {
        impl $trait for $type
        {
            type Output = $type;
            fn $trait_fn(self, rhs: $trait_t) -> Self::Output
            {
                <$type>::from_xyzw(self.x.$trait_fn(rhs as $subtype), self.y.$trait_fn(rhs as $subtype), self.z.$trait_fn(rhs as $subtype), self.w.$trait_fn(rhs as $subtype))
            }
        }
    };
}

macro_rules! impl_single_assign_operator_4 {
    ($type:ty, $subtype:ty, $trait:ty, $trait_t:ty, $trait_fn:ident) => {
        impl $trait for $type
        {
            fn $trait_fn(&mut self, rhs: $trait_t)
            {
                self.x.$trait_fn(rhs as $subtype);
                self.y.$trait_fn(rhs as $subtype);
                self.z.$trait_fn(rhs as $subtype);
                self.w.$trait_fn(rhs as $subtype);
            }
        }
    };
}

macro_rules! impl_operators_4 {
    ($type:ty, $subtype:ty) => {

        impl_default_operator_4!($type, Sub, sub);
        impl_default_assign_operator_4!($type, SubAssign, sub_assign);
        impl_single_operator_4!($type, $subtype, Sub<i32>, i32, sub);
        impl_single_assign_operator_4!($type, $subtype, SubAssign<i32>, i32, sub_assign);
        impl_single_operator_4!($type, $subtype, Sub<u32>, u32, sub);
        impl_single_assign_operator_4!($type, $subtype, SubAssign<u32>, u32, sub_assign);
        impl_single_operator_4!($type, $subtype, Sub<f32>, f32, sub);
        impl_single_assign_operator_4!($type, $subtype, SubAssign<f32>, f32, sub_assign);

        impl_default_operator_4!($type, Add, add);
        impl_default_assign_operator_4!($type, AddAssign, add_assign);
        impl_single_operator_4!($type, $subtype, Add<i32>, i32, add);
        impl_single_assign_operator_4!($type, $subtype, AddAssign<i32>, i32, add_assign);
        impl_single_operator_4!($type, $subtype, Add<u32>, u32, add);
        impl_single_assign_operator_4!($type, $subtype, AddAssign<u32>, u32, add_assign);
        impl_single_operator_4!($type, $subtype, Add<f32>, f32, add);
        impl_single_assign_operator_4!($type, $subtype, AddAssign<f32>, f32, add_assign);

        impl_default_operator_4!($type, Mul, mul);
        impl_default_assign_operator_4!($type, MulAssign, mul_assign);
        impl_single_operator_4!($type, $subtype, Mul<i32>, i32, mul);
        impl_single_assign_operator_4!($type, $subtype, MulAssign<i32>, i32, mul_assign);
        impl_single_operator_4!($type, $subtype, Mul<u32>, u32, mul);
        impl_single_assign_operator_4!($type, $subtype, MulAssign<u32>, u32, mul_assign);
        impl_single_operator_4!($type, $subtype, Mul<f32>, f32, mul);
        impl_single_assign_operator_4!($type, $subtype, MulAssign<f32>, f32, mul_assign);
    };
}

macro_rules! impl_additional_operators_4 {
    ($type:ty) => {
        impl Neg for $type
        {
            type Output = $type;

            fn neg(self) -> Self::Output {
                <$type>::from_xyzw(-self.x, -self.y, -self.z, -self.w)
            }
        }
    };
}

impl_operators_4!(Int4, i32);
impl_additional_operators_4!(Int4);
impl_operators_4!(Float4, f32);
impl_additional_operators_4!(Float4);
impl_operators_4!(Uint4, u32);
impl_default_operator_4!(Float4, Div, div);
impl_single_operator_4!(Float4, f32, Div<f32>, f32, div);
impl_single_operator_4!(Int4, i32, Shl<i32>, i32, shl);
impl_single_operator_4!(Int4, i32, Shr<i32>, i32, shr);

pub fn min_f2(a: &Float2, b: &Float2 ) -> Float2
{
    return Float2::from_xy( a.x.min(b.x) , a.y.min(b.y) );
}

pub fn min_f3(a: &Float3, b: &Float3 ) -> Float3
{
    return Float3::from_xyz( a.x.min(b.x), a.y.min(b.y), a.z.min( b.z ) );
}

pub fn min_f4(a: &Float4, b: &Float4 ) -> Float4
{
    return Float4::from_xyzw( a.x.min(b.x), a.y.min(b.y), a.z.min( b.z ) , a.w.min( b.w ) )
}

pub fn min_i2(a: &Int2, b: &Int2 ) -> Int2
{
    return Int2::from_xy( min( a.x, b.x ), min( a.y, b.y ) );
}

pub fn min_i3(a: &Int3, b: &Int3 ) -> Int3
{
    return Int3::from_xyz( min( a.x, b.x ), min( a.y, b.y ), min( a.z, b.z ) );
}

pub fn min_i4(a: &Int4, b: &Int4 ) -> Int4
{
    return Int4::from_xyzw( min( a.x, b.x ), min( a.y, b.y ), min( a.z, b.z ), min( a.w, b.w ) )
}

pub fn min_u2(a: &Uint2, b: &Uint2 ) -> Uint2
{
    return Uint2::from_xy( min( a.x, b.x ), min( a.y, b.y ) );
}

pub fn min_u3(a: &Uint3, b: &Uint3 ) -> Uint3
{
    return Uint3::from_xyz( min( a.x, b.x ), min( a.y, b.y ), min( a.z, b.z ) );
}

pub fn min_u4(a: &Uint4, b: &Uint4 ) -> Uint4
{
    return Uint4::from_xyzw( min( a.x, b.x ), min( a.y, b.y ), min( a.z, b.z ), min( a.w, b.w ) )
}

pub fn max_f2(a: &Float2, b: &Float2 ) -> Float2
{
    return Float2::from_xy( a.x.max( b.x ), a.y.max( b.y ) );
}

pub fn max_f3(a: &Float3, b: &Float3 ) -> Float3
{
    return Float3::from_xyz( a.x.max( b.x ), a.y.max( b.y ) , a.z.max( b.z ) );
}

pub fn max_f4(a: &Float4, b: &Float4 ) -> Float4
{
    return Float4::from_xyzw( a.x.max( b.x ), a.y.max( b.y ) , a.z.max( b.z ) , a.w.max(b.w ) )
}

pub fn max_i2(a: &Int2, b: &Int2 ) -> Int2
{
    return Int2::from_xy( max( a.x, b.x ), max( a.y, b.y ) );
}

pub fn max_i3(a: &Int3, b: &Int3 ) -> Int3
{
    return Int3::from_xyz( max( a.x, b.x ), max( a.y, b.y ), max( a.z, b.z ) );
}

pub fn max_i4(a: &Int4, b: &Int4 ) -> Int4
{
    return Int4::from_xyzw( max( a.x, b.x ), max( a.y, b.y ), max( a.z, b.z ), max( a.w, b.w ) )
}

pub fn max_u2(a: &Uint2, b: &Uint2 ) -> Uint2
{
    return Uint2::from_xy( max( a.x, b.x ), max( a.y, b.y ) );
}

pub fn max_u3(a: &Uint3, b: &Uint3 ) -> Uint3
{
    return Uint3::from_xyz( max( a.x, b.x ), max( a.y, b.y ), max( a.z, b.z ) );
}

pub fn max_u4(a: &Uint4, b: &Uint4 ) -> Uint4
{
    return Uint4::from_xyzw( max( a.x, b.x ), max( a.y, b.y ), max( a.z, b.z ), max( a.w, b.w ) )
}

pub fn lerp(a: f32, b: f32, t: f32) -> f32
{
    return a + t * (b - a);
}

pub fn lerp_f2(a: &Float2, b: &Float2, t: f32) -> Float2
{
    return *a + (*b - *a) * t;
}

pub fn lerp_f3(a: &Float3, b: &Float3, t: f32) -> Float3
{
    return *a + (*b - *a) * t;
}

pub fn lerp_f4(a: &Float4, b: &Float4, t: f32) -> Float4
{
    return *a + (*b - *a) * t;
}

pub fn clamp_f(f: f32, a: f32, b: f32) -> f32
{
    return a.max(f.min( b ) );
}

pub fn clamp_i(f: i32, a: i32, b: i32) -> i32
{
    return max( a, min( f, b ) );
}

pub fn clamp_u(f: u32, a: u32, b: u32) -> u32
{
    return max( a, min( f, b ) );
}

pub fn clamp_f2(f: &Float2, a: f32, b: f32) -> Float2
{
    return Float2::from_xy(clamp_f(f.x, a, b), clamp_f(f.y, a, b));
}

pub fn clamp_i2(f: &Int2, a: i32, b: i32) -> Int2
{
    return Int2::from_xy(clamp_i(f.x, a, b), clamp_i(f.y, a, b));
}

pub fn clamp_u2(f: &Uint2, a: u32, b: u32) -> Uint2
{
    return Uint2::from_xy(clamp_u(f.x, a, b), clamp_u(f.y, a, b));
}

pub fn clamp_f3(f: &Float3, a: f32, b: f32) -> Float3
{
    return Float3::from_xyz(clamp_f(f.x, a, b), clamp_f(f.y, a, b), clamp_f(f.z, a, b));
}

pub fn clamp_i3(f: &Int3, a: i32, b: i32) -> Int3
{
    return Int3::from_xyz(clamp_i(f.x, a, b), clamp_i(f.y, a, b), clamp_i(f.z, a, b));
}

pub fn clamp_u3(f: &Uint3, a: u32, b: u32) -> Uint3
{
    return Uint3::from_xyz(clamp_u(f.x, a, b), clamp_u(f.y, a, b), clamp_u(f.z, a, b));
}

pub fn clamp_f4(f: &Float4, a: f32, b: f32) -> Float4
{
    return Float4::from_xyzw(clamp_f(f.x, a, b), clamp_f(f.y, a, b), clamp_f(f.z, a, b), clamp_f(f.w, a, b));
}

pub fn clamp_i4(f: &Int4, a: i32, b: i32) -> Int4
{
    return Int4::from_xyzw(clamp_i(f.x, a, b), clamp_i(f.y, a, b), clamp_i(f.z, a, b), clamp_i(f.w, a, b));
}

pub fn clamp_u4(f: &Uint4, a: u32, b: u32) -> Uint4
{
    return Uint4::from_xyzw(clamp_u(f.x, a, b), clamp_u(f.y, a, b), clamp_u(f.z, a, b), clamp_u(f.w, a, b));
}

pub fn dot_f2 (a: &Float2, b: &Float2) -> f32
{
    return a.x * b.x + a.y * b.y;
}

pub fn dot_f3 (a: &Float3, b: &Float3) -> f32
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

pub fn dot_f4 (a: &Float4, b: &Float4) -> f32
{
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

pub fn dot_i2 (a: &Int2, b: &Int2) -> i32
{
    return a.x * b.x + a.y * b.y;
}

pub fn dot_i3 (a: &Int3, b: &Int3) -> i32
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

pub fn dot_i4 (a: &Int4, b: &Int4) -> i32
{
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

pub fn dot_u2 (a: &Uint2, b: &Uint2) -> u32
{
    return a.x * b.x + a.y * b.y;
}

pub fn dot_u3 (a: &Uint3, b: &Uint3) -> u32
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

pub fn dot_u4 (a: &Uint4, b: &Uint4) -> u32
{
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

pub fn sqr_length_f2( v: &Float2) -> f32
{
    return dot_f2(v, v);
}

pub fn sqr_length_f3( v: &Float3) -> f32
{
    return dot_f3(v, v);
}

pub fn sqr_length_f4( v: &Float4) -> f32
{
    return dot_f4(v, v);
}

pub fn length_f2(v: &Float2) -> f32
{
    return dot_f2(v,v).sqrt();
}

pub fn length_f3(v: &Float3) -> f32
{
    return dot_f3(v,v).sqrt();
}

pub fn length_f4(v: &Float4) -> f32
{
    return dot_f4(v,v).sqrt();
}

pub fn length_i2(v: &Int2) -> f32
{
    return (dot_i2(v,v) as f32).sqrt();
}

pub fn length_i3(v: &Int3) -> f32
{
    return (dot_i3(v,v) as f32).sqrt();
}

pub fn length_i4(v: &Int4) -> f32
{
    return (dot_i4(v,v) as f32).sqrt();
}

pub fn normalize_f2(v: &Float2) -> Float2
{
    return *v * (1.0 / length_f2(v));
}

pub fn normalize_f3(v: &Float3) -> Float3
{
    return *v * (1.0 / length_f3(v));
}

pub fn normalize_f4(v: &Float4) -> Float4
{
    return *v * (1.0 / length_f4(v));
}

pub fn reflect(i: &Float3, n: &Float3) -> Float3
{
    return *i - *n * dot_f3( n, i ) * 2.0;
}

pub fn cross(a: &Float3, b: &Float3) -> Float3
{
    return Float3::from_xyz( a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x );
}

#[repr(C, align(64))]
#[derive(Debug, Clone, Copy)]
pub struct Mat4
{
    pub cell: [f32; 16]
}

impl Mat4
{
    pub fn identity_matrix() -> Self
    {
        Mat4 { cell: [ 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0 ] }
    }

    pub fn zero_matrix() -> Self
    {
        Mat4 { cell: [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ] }
    }

    pub fn get_translation(&self) -> Float3
    {
        Float3::from_xyz(self.cell[3], self.cell[7], self.cell[11])
    }

    pub fn rotate_x(a: f32) -> Self
    {
        let mut r = Mat4::identity_matrix();
        let ca = a.cos();
        let sa = a.sin();
        r.cell[5] = ca;
        r.cell[6] = -sa;
        r.cell[9] = sa;
        r.cell[10] = ca;
        return r;
    }

    pub fn rotate_y(a: f32) -> Self
    {
        let mut r = Mat4::identity_matrix();
        let ca = a.cos();
        let sa = a.sin();
        r.cell[0] = ca;
        r.cell[2] = sa;
        r.cell[8] = -sa;
        r.cell[10] = ca;
        return r;
    }

    pub fn rotate_z(a: f32) -> Self
    {
        let mut r = Mat4::identity_matrix();
        let ca = a.cos();
        let sa = a.sin();
        r.cell[0] = ca;
        r.cell[1] = -sa;
        r.cell[4] = sa;
        r.cell[5] = ca;
        return r;
    }

    pub fn translate(p: &Float3) -> Self
    {
        let mut ret_val = Mat4::identity_matrix();
        ret_val.cell[3] = p.x;
        ret_val.cell[7] = p.y;
        ret_val.cell[11] = p.z;
        return ret_val;
    }

    pub fn inverted(&self) -> Mat4
    {
        let inv = [
            self.cell[5] * self.cell[10] * self.cell[15] - self.cell[5] * self.cell[11] * self.cell[14] - self.cell[9] * self.cell[6] * self.cell[15] +
                self.cell[9] * self.cell[7] * self.cell[14] + self.cell[13] * self.cell[6] * self.cell[11] - self.cell[13] * self.cell[7] * self.cell[10],
            -self.cell[1] * self.cell[10] * self.cell[15] + self.cell[1] * self.cell[11] * self.cell[14] + self.cell[9] * self.cell[2] * self.cell[15] -
                self.cell[9] * self.cell[3] * self.cell[14] - self.cell[13] * self.cell[2] * self.cell[11] + self.cell[13] * self.cell[3] * self.cell[10],
            self.cell[1] * self.cell[6] * self.cell[15] - self.cell[1] * self.cell[7] * self.cell[14] - self.cell[5] * self.cell[2] * self.cell[15] +
                self.cell[5] * self.cell[3] * self.cell[14] + self.cell[13] * self.cell[2] * self.cell[7] - self.cell[13] * self.cell[3] * self.cell[6],
            -self.cell[1] * self.cell[6] * self.cell[11] + self.cell[1] * self.cell[7] * self.cell[10] + self.cell[5] * self.cell[2] * self.cell[11] -
                self.cell[5] * self.cell[3] * self.cell[10] - self.cell[9] * self.cell[2] * self.cell[7] + self.cell[9] * self.cell[3] * self.cell[6],
            -self.cell[4] * self.cell[10] * self.cell[15] + self.cell[4] * self.cell[11] * self.cell[14] + self.cell[8] * self.cell[6] * self.cell[15] -
                self.cell[8] * self.cell[7] * self.cell[14] - self.cell[12] * self.cell[6] * self.cell[11] + self.cell[12] * self.cell[7] * self.cell[10],
            self.cell[0] * self.cell[10] * self.cell[15] - self.cell[0] * self.cell[11] * self.cell[14] - self.cell[8] * self.cell[2] * self.cell[15] +
                self.cell[8] * self.cell[3] * self.cell[14] + self.cell[12] * self.cell[2] * self.cell[11] - self.cell[12] * self.cell[3] * self.cell[10],
            -self.cell[0] * self.cell[6] * self.cell[15] + self.cell[0] * self.cell[7] * self.cell[14] + self.cell[4] * self.cell[2] * self.cell[15] -
                self.cell[4] * self.cell[3] * self.cell[14] - self.cell[12] * self.cell[2] * self.cell[7] + self.cell[12] * self.cell[3] * self.cell[6],
            self.cell[0] * self.cell[6] * self.cell[11] - self.cell[0] * self.cell[7] * self.cell[10] - self.cell[4] * self.cell[2] * self.cell[11] +
                self.cell[4] * self.cell[3] * self.cell[10] + self.cell[8] * self.cell[2] * self.cell[7] - self.cell[8] * self.cell[3] * self.cell[6],
            self.cell[4] * self.cell[9] * self.cell[15] - self.cell[4] * self.cell[11] * self.cell[13] - self.cell[8] * self.cell[5] * self.cell[15] +
                self.cell[8] * self.cell[7] * self.cell[13] + self.cell[12] * self.cell[5] * self.cell[11] - self.cell[12] * self.cell[7] * self.cell[9],
            -self.cell[0] * self.cell[9] * self.cell[15] + self.cell[0] * self.cell[11] * self.cell[13] + self.cell[8] * self.cell[1] * self.cell[15] -
                self.cell[8] * self.cell[3] * self.cell[13] - self.cell[12] * self.cell[1] * self.cell[11] + self.cell[12] * self.cell[3] * self.cell[9],
            self.cell[0] * self.cell[5] * self.cell[15] - self.cell[0] * self.cell[7] * self.cell[13] - self.cell[4] * self.cell[1] * self.cell[15] +
                self.cell[4] * self.cell[3] * self.cell[13] + self.cell[12] * self.cell[1] * self.cell[7] - self.cell[12] * self.cell[3] * self.cell[5],
            -self.cell[0] * self.cell[5] * self.cell[11] + self.cell[0] * self.cell[7] * self.cell[9] + self.cell[4] * self.cell[1] * self.cell[11] -
                self.cell[4] * self.cell[3] * self.cell[9] - self.cell[8] * self.cell[1] * self.cell[7] + self.cell[8] * self.cell[3] * self.cell[5],
            -self.cell[4] * self.cell[9] * self.cell[14] + self.cell[4] * self.cell[10] * self.cell[13] + self.cell[8] * self.cell[5] * self.cell[14] -
                self.cell[8] * self.cell[6] * self.cell[13] - self.cell[12] * self.cell[5] * self.cell[10] + self.cell[12] * self.cell[6] * self.cell[9],
            self.cell[0] * self.cell[9] * self.cell[14] - self.cell[0] * self.cell[10] * self.cell[13] - self.cell[8] * self.cell[1] * self.cell[14] +
                self.cell[8] * self.cell[2] * self.cell[13] + self.cell[12] * self.cell[1] * self.cell[10] - self.cell[12] * self.cell[2] * self.cell[9],
            -self.cell[0] * self.cell[5] * self.cell[14] + self.cell[0] * self.cell[6] * self.cell[13] + self.cell[4] * self.cell[1] * self.cell[14] -
                self.cell[4] * self.cell[2] * self.cell[13] - self.cell[12] * self.cell[1] * self.cell[6] + self.cell[12] * self.cell[2] * self.cell[5],
            self.cell[0] * self.cell[5] * self.cell[10] - self.cell[0] * self.cell[6] * self.cell[9] - self.cell[4] * self.cell[1] * self.cell[10] +
                self.cell[4] * self.cell[2] * self.cell[9] + self.cell[8] * self.cell[1] * self.cell[6] - self.cell[8] * self.cell[2] * self.cell[5]
        ];
        let det = self.cell[0] * inv[0] + self.cell[1] * inv[4] + self.cell[2] * inv[8] + self.cell[3] * inv[12];
        let mut ret_val = Mat4::identity_matrix();
        if det != 0.0
        {
            let inv_det = 1.0 / det;
            for i in 0..16
            {
                ret_val.cell[i] = inv[i] * inv_det;
            }
        }
        return ret_val;
    }

    pub fn inverted_no_scale(&self) -> Mat4
    {
        let mut r = Mat4::identity_matrix();
        r.cell[0] = self.cell[0];
        r.cell[1] = self.cell[4];
        r.cell[2] = self.cell[8];
        r.cell[4] = self.cell[1];
        r.cell[5] = self.cell[5];
        r.cell[6] = self.cell[9];
        r.cell[8] = self.cell[2];
        r.cell[9] = self.cell[6];
        r.cell[10] = self.cell[10];
        r.cell[3] = -(self.cell[3] * r.cell[0] + self.cell[7] * r.cell[1] + self.cell[11] * r.cell[2]);
        r.cell[7] = -(self.cell[3] * r.cell[4] + self.cell[7] * r.cell[5] + self.cell[11] * r.cell[6]);
        r.cell[11] = -(self.cell[3] * r.cell[8] + self.cell[7] * r.cell[9] + self.cell[11] * r.cell[10]);
        return r;
    }
}

impl Mul<&Mat4> for Float4
{
    type Output = Float4;

    fn mul(self, rhs: &Mat4) -> Self::Output {
        Float4::from_xyzw(
            rhs.cell[0] * self.x + rhs.cell[1] * self.y + rhs.cell[2] * self.z + rhs.cell[3] * self.w,
            rhs.cell[4] * self.x + rhs.cell[5] * self.y + rhs.cell[6] * self.z + rhs.cell[7] * self.w,
            rhs.cell[8] * self.x + rhs.cell[9] * self.y + rhs.cell[10] * self.z + rhs.cell[11] * self.w,
            rhs.cell[12] * self.x + rhs.cell[13] * self.y + rhs.cell[14] * self.z + rhs.cell[15] * self.w
        )
    }
}

impl Mul<&Mat4> for Mat4
{
    type Output = Mat4;

    fn mul(self, rhs: &Mat4) -> Self::Output {
        let mut r = Mat4::identity_matrix();
        for i in (0..16).step_by(4) {
            for j in 0..4 {
                r.cell[i + j] =
                    (self.cell[i + 0] * rhs.cell[j + 0]) +
                        (self.cell[i + 1] * rhs.cell[j + 4]) +
                        (self.cell[i + 2] * rhs.cell[j + 8]) +
                        (self.cell[i + 3] * rhs.cell[j + 12]);
            }
        }
        return r;
    }
}

impl Mul<Mat4> for Mat4
{
    type Output = Mat4;

    fn mul(self, rhs: Mat4) -> Self::Output {
        return self * &rhs;
    }
}

pub fn transform_position(a: &Float3, m: &Mat4) -> Float3
{
    let f4 = Float4::from_xyzw(a.x, a.y, a.z, 1.0) * m;
    Float3::from_xyz(f4.x, f4.y, f4.z)
}

pub fn transform_vector(a: &Float3, m: &Mat4) -> Float3
{
    let f4 = Float4::from_xyzw(a.x, a.y, a.z, 0.0) * m;
    Float3::from_xyz(f4.x, f4.y, f4.z)
}

pub fn rgbf32_to_rgb8_f4(value: &Float4) -> u32
{
    let r: u32 = (255.0 * value.x.min(1.0)) as u32;
    let g: u32 = (255.0 * value.y.min(1.0)) as u32;
    let b: u32 = (255.0 * value.z.min(1.0)) as u32;
    return (r << 16) + (g << 8) + b;
}

pub fn rgbf32_to_rgb8_f3(value: &Float3) -> u32
{
    let r: u32 = (255.0 * value.x.min(1.0)) as u32;
    let g: u32 = (255.0 * value.y.min(1.0)) as u32;
    let b: u32 = (255.0 * value.z.min(1.0)) as u32;
    return (r << 16) + (g << 8) + b;
}