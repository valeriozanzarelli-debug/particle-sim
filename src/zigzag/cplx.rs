//! Minimal f64 complex numbers for the zig-zag walk.

use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct C {
    pub re: f64,
    pub im: f64,
}

impl C {
    pub const ZERO: C = C { re: 0.0, im: 0.0 };
    pub const ONE: C = C { re: 1.0, im: 0.0 };
    pub const I: C = C { re: 0.0, im: 1.0 };
    #[inline] pub const fn new(re: f64, im: f64) -> Self { Self { re, im } }
    #[inline] pub const fn real(re: f64) -> Self { Self { re, im: 0.0 } }
    /// e^{iθ}
    #[inline] pub fn expi(t: f64) -> Self { let (s, c) = t.sin_cos(); Self { re: c, im: s } }
    #[inline] pub fn conj(self) -> Self { Self { re: self.re, im: -self.im } }
    #[inline] pub fn norm2(self) -> f64 { self.re * self.re + self.im * self.im }
    #[inline] pub fn abs(self) -> f64 { self.norm2().sqrt() }
    #[inline] pub fn arg(self) -> f64 { self.im.atan2(self.re) }
    #[inline] pub fn scale(self, k: f64) -> Self { Self { re: self.re * k, im: self.im * k } }
}

impl Add for C { type Output = C; #[inline] fn add(self, o: C) -> C { C::new(self.re + o.re, self.im + o.im) } }
impl Sub for C { type Output = C; #[inline] fn sub(self, o: C) -> C { C::new(self.re - o.re, self.im - o.im) } }
impl Neg for C { type Output = C; #[inline] fn neg(self) -> C { C::new(-self.re, -self.im) } }
impl AddAssign for C { #[inline] fn add_assign(&mut self, o: C) { self.re += o.re; self.im += o.im; } }
impl SubAssign for C { #[inline] fn sub_assign(&mut self, o: C) { self.re -= o.re; self.im -= o.im; } }
impl Mul for C {
    type Output = C;
    #[inline] fn mul(self, o: C) -> C { C::new(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re) }
}
impl MulAssign for C { #[inline] fn mul_assign(&mut self, o: C) { *self = *self * o; } }
impl Mul<f64> for C { type Output = C; #[inline] fn mul(self, k: f64) -> C { self.scale(k) } }
impl Div for C {
    type Output = C;
    #[inline] fn div(self, o: C) -> C { let d = o.norm2(); (self * o.conj()).scale(1.0 / d) }
}
