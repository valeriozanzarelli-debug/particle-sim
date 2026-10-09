//! Drawing: the balls (one point per electron per ball), the nuclei, the surfer, the node surfaces.

use bevy::{
    prelude::*,
    render::{mesh::PrimitiveTopology, render_asset::RenderAssetUsages},
};

use crate::dmc::Vincolo;
use crate::scena::{punti, Scelta};
use crate::Lab;

/// The entity drawing the ball cloud.
#[derive(Component)]
pub struct Nuvola;

/// A point cloud with the given positions and colours.
fn cloud_mesh(pos: Vec<[f32; 3]>, col: Vec<[f32; 4]>) -> Mesh {
    let n = pos.len();
    let mut mesh = Mesh::new(PrimitiveTopology::PointList, RenderAssetUsages::RENDER_WORLD);
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, pos);
    mesh.insert_attribute(Mesh::ATTRIBUTE_COLOR, col);
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, vec![[0.0f32, 1.0, 0.0]; n]);
    mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, vec![[0.0f32; 2]; n]);
    mesh
}

pub fn setup_render(mut commands: Commands, mut meshes: ResMut<Assets<Mesh>>, mut materials: ResMut<Assets<StandardMaterial>>) {
    let mat = materials.add(StandardMaterial { base_color: Color::WHITE, unlit: true, ..default() });
    let m = cloud_mesh(vec![[1e7; 3]], vec![[0.0; 4]]);
    commands.spawn((Nuvola, Mesh3d(meshes.add(m)), MeshMaterial3d(mat), Transform::default()));
}

/// Colours: electron 1 (spin up) cyan, electron 2 (spin down) orange; brighter near the nuclei.
fn colore(e: usize, r: f32) -> [f32; 4] {
    let k = (1.0 / (1.0 + 0.35 * r)).clamp(0.25, 1.0);
    let base = if e == 0 { [0.35, 0.85, 1.0] } else { [1.0, 0.62, 0.25] };
    [base[0] * (0.45 + 0.55 * k), base[1] * (0.45 + 0.55 * k), base[2] * (0.45 + 0.55 * k), 1.0]
}

/// Replace the cloud mesh with the current balls.
///
/// Bevy 0.15 never frees the GPU allocation of a mesh that is *modified* in place, so the asset is
/// swapped (new handle, old one removed) — the same rule as `swap_mesh` in src/photon/render.rs.
pub fn aggiorna_nuvola(lab: Res<Lab>, mut meshes: ResMut<Assets<Mesh>>, mut q: Query<&mut Mesh3d, With<Nuvola>>) {
    let mut pos = vec![];
    let mut col = vec![];
    for (p, e) in punti(&lab.dmc) {
        if !lab.mostra_e[e.min(1)] { continue; }
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        pos.push([p[0], p[2], p[1]]); // physics z (molecular axis) -> screen y (up)
        col.push(colore(e, r));
    }
    if pos.is_empty() {
        pos.push([1e7; 3]);
        col.push([0.0; 4]);
    }
    for mut m in q.iter_mut() {
        let old = m.0.id();
        m.0 = meshes.add(cloud_mesh(pos.clone(), col.clone()));
        meshes.remove(old);
    }
}

fn v(p: [f64; 3]) -> Vec3 { Vec3::new(p[0] as f32, p[2] as f32, p[1] as f32) }

/// Nuclei, the surfer and the node surfaces.
pub fn disegna(lab: Res<Lab>, mut g: Gizmos) {
    let gold = Color::srgb(1.0, 0.82, 0.3);
    let white = Color::srgb(1.0, 1.0, 1.0);
    let node = Color::srgba(0.7, 0.5, 1.0, 0.35);
    for n in &lab.dmc.sys.nuclei {
        let c = v(n.pos);
        let r = 0.06 * n.z.sqrt() as f32;
        g.sphere(Isometry3d::from_translation(c), r, gold);
        g.sphere(Isometry3d::from_translation(c), r * 0.5, gold);
    }
    // the surfer: one ball (a whole configuration of the electrons), at rest
    if lab.mostra_surfer {
        for p in &lab.surfer {
            let c = v(*p);
            g.sphere(Isometry3d::from_translation(c), 0.09, white);
            for d in [Vec3::X, Vec3::Y, Vec3::Z] { g.line(c - d * 0.3, c + d * 0.3, white); }
        }
        if lab.surfer.len() == 2 { g.line(v(lab.surfer[0]), v(lab.surfer[1]), Color::srgba(1.0, 1.0, 1.0, 0.35)); }
    }
    // node surfaces of the excited states
    let s = lab.scelta;
    let ext = 12.0;
    let piano = |g: &mut Gizmos, n: Vec3, a: Vec3, b: Vec3| {
        for i in -4..=4 {
            let t = i as f32 / 4.0 * ext;
            g.line(a * t - b * ext, a * t + b * ext, node);
            g.line(b * t - a * ext, b * t + a * ext, node);
        }
        let _ = n;
    };
    let mut vincoli: Vec<Vincolo> = lab.dmc.sys.vincoli.clone();
    match s {
        Scelta::H2s => vincoli.push(Vincolo::SferaFuori(lab.res.f("h2s_node").unwrap_or(2.0))),
        Scelta::H3p => {
            vincoli.push(Vincolo::PianoZ);
            vincoli.push(Vincolo::SferaFuori(lab.res.f("h3p_node").unwrap_or(6.0)));
        }
        _ => {}
    }
    for c in vincoli {
        match c {
            // physics z -> screen Y; physics x -> X; physics y -> Z
            Vincolo::PianoZ => piano(&mut g, Vec3::Y, Vec3::X, Vec3::Z),
            Vincolo::PianoX => piano(&mut g, Vec3::X, Vec3::Y, Vec3::Z),
            Vincolo::PianoY => piano(&mut g, Vec3::Z, Vec3::X, Vec3::Y),
            Vincolo::SferaDentro(r) | Vincolo::SferaFuori(r) => {
                g.sphere(Isometry3d::IDENTITY, r as f32, node);
            }
        }
    }
}
