//! Orbit camera: LMB rotates, RMB pans, scroll zooms.

use bevy::{input::mouse::{MouseMotion, MouseWheel}, prelude::*};
use crate::config::{CAM_RADIUS_MIN, CAM_RADIUS_MAX, CAM_RADIUS_DEFAULT};

#[derive(Resource)]
pub struct OrbitState {
    pub yaw:    f32,
    pub pitch:  f32,
    pub radius: f32,
    pub focus:  Vec3,
}

impl Default for OrbitState {
    fn default() -> Self {
        Self { yaw: 0.3, pitch: 0.35, radius: CAM_RADIUS_DEFAULT, focus: Vec3::ZERO }
    }
}

pub fn orbit_camera(
    mut state: ResMut<OrbitState>,
    mut motion: EventReader<MouseMotion>,
    mut wheel: EventReader<MouseWheel>,
    buttons: Res<ButtonInput<MouseButton>>,
    mut camera: Query<&mut Transform, With<Camera3d>>,
) {
    if buttons.pressed(MouseButton::Left) {
        for ev in motion.read() {
            state.yaw   -= ev.delta.x * 0.005;
            state.pitch  = (state.pitch + ev.delta.y * 0.005).clamp(-1.4, 1.4);
        }
    } else if buttons.pressed(MouseButton::Right) {
        // Pan
        let cos_p = state.pitch.cos();
        let cam_pos = orbit_pos(&state);
        let forward = (state.focus - cam_pos).normalize_or_zero();
        if forward.length_squared() > 1e-8 {
            let right = forward.cross(Vec3::Y).normalize_or_zero();
            let up    = right.cross(forward).normalize();
            let pan   = state.radius * 0.0018;
            for ev in motion.read() {
                state.focus += right * (-ev.delta.x * pan) + up * (ev.delta.y * pan);
            }
        }
        let _ = cos_p;
    } else {
        motion.clear();
    }

    for ev in wheel.read() {
        let factor = if ev.y > 0.0 { 0.88 } else { 1.0 / 0.88 };
        state.radius = (state.radius * factor).clamp(CAM_RADIUS_MIN, CAM_RADIUS_MAX);
    }

    let pos = orbit_pos(&state);
    for mut cam in camera.iter_mut() {
        cam.translation = pos;
        cam.look_at(state.focus, Vec3::Y);
    }
}

fn orbit_pos(s: &OrbitState) -> Vec3 {
    let cos_p = s.pitch.cos();
    s.focus + Vec3::new(
        s.radius * cos_p * s.yaw.sin(),
        s.radius * s.pitch.sin(),
        s.radius * cos_p * s.yaw.cos(),
    )
}
