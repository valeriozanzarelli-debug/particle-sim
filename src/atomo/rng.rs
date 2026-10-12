//! xoshiro256++ (same generator as src/molecole/dmc.rs): random positions for drawing balls.

#[derive(Clone, Debug)]
pub struct Rng([u64; 4]);

impl Rng {
    pub fn new(seed: u64) -> Self {
        let mut z = seed;
        let mut s = [0u64; 4];
        for v in s.iter_mut() {
            z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            *v = x ^ (x >> 31);
        }
        Rng(s)
    }
    pub fn next_u64(&mut self) -> u64 {
        let s = &mut self.0;
        let r = s[0].wrapping_add(s[3]).rotate_left(23).wrapping_add(s[0]);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        r
    }
    pub fn uni(&mut self) -> f64 { (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64) }
    pub fn gauss(&mut self) -> f64 {
        let u1 = 1.0 - self.uni();
        let u2 = self.uni();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }
}

/// A fixed pseudo-random number in [0, 1) for an integer (the same cell always gets the same jitter).
#[inline]
pub fn hash01(i: u64) -> f32 {
    let mut x = i.wrapping_mul(0x9E37_79B9_7F4A_7C15).wrapping_add(0x632B_E59B_D9B4_E019);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^= x >> 31;
    (x >> 40) as f32 / (1u64 << 24) as f32
}
