//! Drawing the photon: grains as points (colour = phase), E-field arrows,
//! optical elements and detector clicks.

use bevy::{
    prelude::*,
    render::{
        mesh::{PrimitiveTopology, VertexAttributeValues},
        render_asset::RenderAssetUsages,
    },
};

use crate::config::*;
use crate::physics::{Experiment, TAU};
use crate::LabState;

#[derive(Resource)]
pub struct CloudMeshes {
    grains: Handle<Mesh>,
    clicks: Handle<Mesh>,
}

/// Where dead grains are parked (outside any view).
const PARKED: [f32; 3] = [1.0e7, 1.0e7, 1.0e7];

fn point_mesh(n: usize) -> Mesh {
    let mut mesh = Mesh::new(
        PrimitiveTopology::PointList,
        RenderAssetUsages::MAIN_WORLD | RenderAssetUsages::RENDER_WORLD,
    );
    resize(&mut mesh, n);
    mesh
}

/// Every attribute must have the same length, or Bevy draws the shortest.
fn resize(mesh: &mut Mesh, n: usize) {
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, vec![PARKED; n]);
    mesh.insert_attribute(Mesh::ATTRIBUTE_COLOR,    vec![[0.0f32; 4]; n]);
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL,   vec![[0.0f32, 1.0, 0.0]; n]);
    mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0,     vec![[0.0f32; 2]; n]);
}

pub fn setup_render(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let mat = materials.add(StandardMaterial { base_color: Color::WHITE, unlit: true, ..default() });
    let grains = meshes.add(point_mesh(1));
    let clicks = meshes.add(point_mesh(1));
    commands.spawn((Mesh3d(grains.clone()), MeshMaterial3d(mat.clone()), Transform::default()));
    commands.spawn((Mesh3d(clicks.clone()), MeshMaterial3d(mat), Transform::default()));
    commands.insert_resource(CloudMeshes { grains, clicks });
}

/// Phase colour wheel: the colour of a grain is the angle of its arrow.
fn phase_color(angle: f32) -> [f32; 4] {
    let t = angle.rem_euclid(TAU) / TAU;
    // Cyclic palette: cyan → blue → violet → magenta → gold → cyan
    const STOPS: [[f32; 3]; 5] = [
        [0.20, 0.95, 1.00],
        [0.25, 0.45, 1.00],
        [0.70, 0.30, 1.00],
        [1.00, 0.35, 0.70],
        [1.00, 0.85, 0.30],
    ];
    let x = t * 5.0;
    let i = (x.floor() as usize) % 5;
    let s = x - x.floor();
    let (a, b) = (STOPS[i], STOPS[(i + 1) % 5]);
    [a[0] + (b[0] - a[0]) * s, a[1] + (b[1] - a[1]) * s, a[2] + (b[2] - a[2]) * s, 1.0]
}

pub fn update_meshes(
    state: Res<LabState>,
    handles: Option<Res<CloudMeshes>>,
    mut meshes: ResMut<Assets<Mesh>>,
) {
    let Some(h) = handles else { return };
    let sim = &state.sim;
    let k = sim.lab.k();

    // ── grains ──
    let n: usize = sim.photons.iter().map(|p| p.grains.len()).sum::<usize>().max(1);
    if let Some(mesh) = meshes.get_mut(&h.grains) {
        if mesh.count_vertices() != n { resize(mesh, n); }
        let mut pos: Vec<[f32; 3]> = Vec::with_capacity(n);
        let mut col: Vec<[f32; 4]> = Vec::with_capacity(n);
        for ph in &sim.photons {
            for g in &ph.grains {
                if !g.alive {
                    pos.push(PARKED);
                    col.push([0.0; 4]);
                    continue;
                }
                pos.push(g.pos.to_array());
                let mut c = phase_color(g.phase(k) - k * sim.clock);
                if !ph.pol_defined {
                    // Entangled, unmeasured: phase yes, polarization none → pale
                    for v in c.iter_mut().take(3) { *v = 0.55 + 0.45 * *v * 0.5; }
                }
                col.push(c);
            }
        }
        pos.resize(n, PARKED);
        col.resize(n, [0.0; 4]);
        set_points(mesh, pos, col);
    }

    // ── clicks on the screen ──
    let clicks = &sim.stats.clicks;
    let m = clicks.len().min(MAX_CLICKS_DRAWN).max(1);
    if let Some(mesh) = meshes.get_mut(&h.clicks) {
        if mesh.count_vertices() != m { resize(mesh, m); }
        let start = clicks.len().saturating_sub(MAX_CLICKS_DRAWN);
        let mut pos: Vec<[f32; 3]> = clicks[start..].iter().map(|c| c.to_array()).collect();
        let mut col: Vec<[f32; 4]> = vec![[0.75, 1.0, 0.6, 1.0]; pos.len()];
        pos.resize(m, PARKED);
        col.resize(m, [0.0; 4]);
        set_points(mesh, pos, col);
    }
}

fn set_points(mesh: &mut Mesh, pos: Vec<[f32; 3]>, col: Vec<[f32; 4]>) {
    if let Some(VertexAttributeValues::Float32x3(p)) = mesh.attribute_mut(Mesh::ATTRIBUTE_POSITION) { *p = pos; }
    if let Some(VertexAttributeValues::Float32x4(c)) = mesh.attribute_mut(Mesh::ATTRIBUTE_COLOR) { *c = col; }
}

/// Optical bench and E-field arrows.
pub fn draw_bench(state: Res<LabState>, mut gizmos: Gizmos) {
    let sim = &state.sim;
    let lab = &sim.lab;
    let dim   = Color::srgba(0.55, 0.6, 0.7, 0.6);
    let gold  = Color::srgb(1.0, 0.8, 0.3);
    let green = Color::srgb(0.6, 1.0, 0.5);

    let rect = |g: &mut Gizmos, z: f32, x0: f32, x1: f32, y0: f32, y1: f32, c: Color| {
        let p = [Vec3::new(x0, y0, z), Vec3::new(x1, y0, z), Vec3::new(x1, y1, z), Vec3::new(x0, y1, z)];
        for i in 0..4 { g.line(p[i], p[(i + 1) % 4], c); }
    };
    let polarizer = |g: &mut Gizmos, z: f32, angle: f32, c: Color| {
        let r = 10.0;
        let iso = Isometry3d::new(Vec3::new(0.0, 0.0, z), Quat::IDENTITY);
        g.circle(iso, r, c);
        let a = Vec3::new(angle.cos(), angle.sin(), 0.0) * r;
        g.line(Vec3::new(0.0, 0.0, z) - a, Vec3::new(0.0, 0.0, z) + a, c);
    };

    // Source
    if lab.experiment != Experiment::Bell {
        gizmos.circle(Isometry3d::new(Vec3::new(0.0, 0.0, lab.source_z - 12.0), Quat::IDENTITY), 2.0, gold);
    } else {
        gizmos.sphere(Isometry3d::from_translation(Vec3::ZERO), 1.5, gold);
    }

    match lab.experiment {
        Experiment::FreeBeam => {}
        Experiment::DoubleSlit => {
            let (c, h) = lab.slits();
            let (x0, x1, yh) = (-40.0, 40.0, 25.0);
            let z = lab.wall_z;
            // Wall pieces between/around the slits
            rect(&mut gizmos, z, x0, c[0] - h, -yh, yh, dim);
            rect(&mut gizmos, z, c[0] + h, c[1] - h, -yh, yh, dim);
            rect(&mut gizmos, z, c[1] + h, x1, -yh, yh, dim);
            // Which-path polarizers on the slits
            for (s, f) in lab.slit_filters().iter().enumerate() {
                if let Some(f) = f {
                    let a = f[1].re.atan2(f[0].re);
                    let d = Vec3::new(a.cos(), a.sin(), 0.0) * 3.0;
                    let p = Vec3::new(c[s], 0.0, z - 0.5);
                    gizmos.line(p - d, p + d, gold);
                }
            }
        }
        Experiment::Polarizers => {
            for p in lab.polarizers.iter().filter(|p| p.enabled) {
                polarizer(&mut gizmos, p.z, p.angle, gold);
            }
        }
        Experiment::Bell => {
            let last = sim.stats.last_bell;
            for (side, d) in [(0usize, lab.bell.dist_a), (1, lab.bell.dist_b)] {
                let z = if side == 0 { d } else { -d };
                let angle = match (last, lab.bell.random_chsh) {
                    (Some((a, _, b, _)), true) => if side == 0 { a } else { b },
                    _ => if side == 0 { lab.bell.angle_a } else { lab.bell.angle_b },
                };
                let c = if side == 0 { Color::srgb(0.4, 0.9, 1.0) } else { Color::srgb(1.0, 0.6, 0.3) };
                polarizer(&mut gizmos, z, angle, c);
                let zd = if side == 0 { d + 15.0 } else { -d - 15.0 };
                rect(&mut gizmos, zd, -6.0, 6.0, -6.0, 6.0, dim);
            }
        }
    }

    // Detection screen
    if lab.experiment != Experiment::Bell {
        rect(&mut gizmos, lab.screen_z, -lab.screen_half, lab.screen_half, -25.0, 25.0, green.with_alpha(0.5));
    }

    // E-field arrows on a subset of grains: rotate + translate made visible
    if state.show_arrows {
        let k = lab.k();
        for ph in sim.photons.iter().filter(|p| p.pol_defined) {
            let alive: Vec<_> = ph.grains.iter().filter(|g| g.alive).collect();
            let stride = (alive.len() / MAX_ARROWS).max(1);
            for g in alive.iter().step_by(stride) {
                let e = g.e_field(k, sim.clock) * ARROW_LEN;
                gizmos.line(g.pos, g.pos + e, Color::srgba(1.0, 1.0, 1.0, 0.55));
            }
        }
    }
}
