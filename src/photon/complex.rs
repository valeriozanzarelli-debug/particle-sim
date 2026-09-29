//! Minimal complex arithmetic (phases, Jones vectors, two-photon states).

use std::ops::{Add, AddAssign, Mul, Neg, Sub};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct C {
    pub re: f32,
    pub im: f32,
}

impl C {
    pub const ZERO: C = C { re: 0.0, im: 0.0 };
    pub const ONE:  C = C { re: 1.0, im: 0.0 };
    pub const I:    C = C { re: 0.0, im: 1.0 };

    #[inline] pub const fn new(re: f32, im: f32) -> Self { Self { re, im } }
    #[inline] pub fn real(re: f32) -> Self { Self { re, im: 0.0 } }
    /// e^{iθ}
    #[inline] pub fn expi(theta: f32) -> Self { let (s, c) = theta.sin_cos(); Self { re: c, im: s } }
    #[inline] pub fn conj(self) -> Self { Self { re: self.re, im: -self.im } }
    #[inline] pub fn norm2(self) -> f32 { self.re * self.re + self.im * self.im }
    #[inline] pub fn scale(self, k: f32) -> Self { Self { re: self.re * k, im: self.im * k } }
}

impl Add for C { type Output = C; #[inline] fn add(self, o: C) -> C { C::new(self.re + o.re, self.im + o.im) } }
impl Sub for C { type Output = C; #[inline] fn sub(self, o: C) -> C { C::new(self.re - o.re, self.im - o.im) } }
impl Neg for C { type Output = C; #[inline] fn neg(self) -> C { C::new(-self.re, -self.im) } }
impl AddAssign for C { #[inline] fn add_assign(&mut self, o: C) { self.re += o.re; self.im += o.im; } }
impl Mul for C {
    type Output = C;
    #[inline] fn mul(self, o: C) -> C {
        C::new(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    }
}
impl Mul<f32> for C { type Output = C; #[inline] fn mul(self, k: f32) -> C { self.scale(k) } }

/// Polarization state in a beam's transverse basis (e1, e2).
pub type Jones = [C; 2];

/// ⟨a|b⟩
#[inline]
pub fn inner(a: Jones, b: Jones) -> C {
    a[0].conj() * b[0] + a[1].conj() * b[1]
}

#[inline]
pub fn jones_norm2(j: Jones) -> f32 {
    j[0].norm2() + j[1].norm2()
}

/// Linear polarization at angle `a` from e1.
#[inline]
pub fn linear(a: f32) -> Jones {
    [C::real(a.cos()), C::real(a.sin())]
}
