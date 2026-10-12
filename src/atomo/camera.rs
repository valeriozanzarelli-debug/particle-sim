//! Camera with a LOGARITHMIC zoom (one wheel notch = the same factor at every scale, from 100 Å down to
//! 0.1 fm): left drag orbits, wheel zooms, smooth exponential approach to the target. The focus is kept in
//! f64 (physics coordinates, bohr) and the scene is drawn relative to it (render.rs), so the femtometre
//! view stays exact even when the nucleus has moved by ångströms.

use bevy::{input::mouse::{MouseMotion, MouseScrollUnit, MouseWheel}, prelude::*, render::camera::Projection};

/// What the camera looks at.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Bersaglio { Atomo, Nucleo, Elettrone }

#[derive(Resource)]
pub struct Cam {
    pub yaw: f32,
    pub pitch: f32,
    /// ln of the camera distance (bohr) and its target
    pub ln_r: f64,
    pub ln_r_target: f64,
    /// focus (bohr, physics coordinates) and where it goes
    pub focus: [f64; 3],
    pub bersaglio: Bersaglio,
    pub fov: f32,
    pub altezza_px: f32,
}

pub const R_MIN: f64 = 2.0e-6; // ≈ 0.1 fm
pub const R_MAX: f64 = 400.0; // ≈ 200 Å

impl Cam {
    pub fn new(r: f64, b: Bersaglio) -> Self {
        Cam { yaw: 0.6, pitch: 0.35, ln_r: r.ln(), ln_r_target: r.ln(), focus: [0.0; 3], bersaglio: b, fov: 0.8, altezza_px: 900.0 }
    }
    pub fn r(&self) -> f64 { self.ln_r.exp() }
    pub fn vai(&mut self, b: Bersaglio, r: f64) {
        self.bersaglio = b;
        self.ln_r_target = r.clamp(R_MIN, R_MAX).ln();
    }
    /// World length (bohr) of one pixel at the focus.
    pub fn bohr_per_px(&self) -> f64 { 2.0 * self.r() * (self.fov as f64 * 0.5).tan() / self.altezza_px as f64 }
}

#[derive(Resource, Default)]
pub struct PointerOverHud(pub bool);

#[allow(clippy::too_many_arguments)]
pub fn camera_sys(
    time: Res<Time>,
    mut cam: ResMut<Cam>,
    over: Res<PointerOverHud>,
    mut motion: EventReader<MouseMotion>,
    mut wheel: EventReader<MouseWheel>,
    buttons: Res<ButtonInput<MouseButton>>,
    lab: Res<crate::Lab>,
    auto: Res<crate::Auto>,
    windows: Query<&Window>,
    mut q: Query<(&mut Transform, &mut Projection), With<Camera3d>>,
) {
    let dt = if auto.attivo { 1.0 / 60.0 } else { time.delta_secs_f64().min(0.1) };
    if auto.attivo { motion.clear(); wheel.clear(); }
    if !over.0 && buttons.pressed(MouseButton::Left) {
        for ev in motion.read() {
            cam.yaw -= ev.delta.x * 0.005;
            cam.pitch = (cam.pitch + ev.delta.y * 0.005).clamp(-1.5, 1.5);
        }
    } else {
        motion.clear();
    }
    for ev in wheel.read() {
        if over.0 { continue; }
        let notches = match ev.unit { MouseScrollUnit::Line => ev.y as f64, MouseScrollUnit::Pixel => ev.y as f64 / 60.0 };
        cam.ln_r_target = (cam.ln_r_target - 0.18 * notches).clamp(R_MIN.ln(), R_MAX.ln());
    }
    // smooth exponential approach (zoom in log space; the focus in a straight line)
    let a = 1.0 - (-6.0 * dt).exp();
    cam.ln_r += (cam.ln_r_target - cam.ln_r) * a;
    let target = lab.punto_bersaglio(cam.bersaglio);
    // the focus moves at a speed proportional to the view size, so long flights stay smooth
    for c in 0..3 { cam.focus[c] += (target[c] - cam.focus[c]) * a; }
    if let Ok(w) = windows.get_single() { cam.altezza_px = w.height().max(100.0); }
    let r = cam.r() as f32;
    let pos = Vec3::new(cam.yaw.sin() * cam.pitch.cos(), cam.pitch.sin(), cam.yaw.cos() * cam.pitch.cos()) * r;
    for (mut t, mut p) in q.iter_mut() {
        *t = Transform::from_translation(pos).looking_at(Vec3::ZERO, Vec3::Y);
        if let Projection::Perspective(pp) = &mut *p {
            pp.near = r * 1e-3;
            pp.far = r * 1e4;
            pp.fov = cam.fov;
        }
    }
}
