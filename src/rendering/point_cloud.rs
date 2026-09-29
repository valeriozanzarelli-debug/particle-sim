//! Point-cloud renderer: one PointList Mesh updated every frame.
//! Each particle maps to a single vertex (= one pixel when rendered).
//!
//! Colour encodes local density (neighbour count):
//!   blue  → cyan → green → yellow → red → white
//! This reveals emerging structures: cold lonely particles vs hot dense clusters.

use bevy::{
    prelude::*,
    render::{
        mesh::{PrimitiveTopology, VertexAttributeValues},
        render_asset::RenderAssetUsages,
    },
};
use crate::simulation::world::ParticleWorld;
use crate::config::COLOR_SCALE_N;

// ─── Handle resource ──────────────────────────────────────────────────────────

#[derive(Resource)]
pub struct ParticleMesh(pub Handle<Mesh>);

// ─── Setup: create mesh + material once ───────────────────────────────────────

pub fn setup_point_cloud(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    world: Res<ParticleWorld>,
) {
    let n = world.n;

    let positions: Vec<[f32; 3]> = world.positions.iter().map(|p| p.to_array()).collect();
    let colors:    Vec<[f32; 4]> = world.neighbor_counts.iter().map(|&c| density_color(c)).collect();
    // Dummy normals and UVs required by Bevy's default shader
    let normals: Vec<[f32; 3]> = vec![[0.0, 1.0, 0.0]; n];
    let uvs:     Vec<[f32; 2]> = vec![[0.0, 0.0]; n];

    let mut mesh = Mesh::new(
        PrimitiveTopology::PointList,
        RenderAssetUsages::MAIN_WORLD | RenderAssetUsages::RENDER_WORLD,
    );
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, positions);
    mesh.insert_attribute(Mesh::ATTRIBUTE_COLOR,    colors);
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL,   normals);
    mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0,     uvs);

    let handle = meshes.add(mesh);
    commands.insert_resource(ParticleMesh(handle.clone()));

    let mat = materials.add(StandardMaterial {
        base_color: Color::WHITE,
        unlit: true,
        ..default()
    });

    commands.spawn((Mesh3d(handle), MeshMaterial3d(mat), Transform::default()));
}

// ─── Per-frame update: positions + colours in-place ──────────────────────────

pub fn update_point_cloud(
    world: Res<ParticleWorld>,
    handle: Option<Res<ParticleMesh>>,
    mut meshes: ResMut<Assets<Mesh>>,
) {
    let Some(handle) = handle else { return };
    let Some(mesh)   = meshes.get_mut(&handle.0) else { return };

    // Particle count changed (after reset): every attribute must be resized,
    // or Bevy draws only as many vertices as the shortest attribute holds.
    if mesh.count_vertices() != world.n {
        mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, vec![[0.0f32, 1.0, 0.0]; world.n]);
        mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0,   vec![[0.0f32, 0.0]; world.n]);
    }

    // Update positions
    if let Some(VertexAttributeValues::Float32x3(ref mut pos)) =
        mesh.attribute_mut(Mesh::ATTRIBUTE_POSITION)
    {
        // Grow or shrink if particle count changed (after reset)
        if pos.len() != world.n {
            *pos = world.positions.iter().map(|p| p.to_array()).collect();
        } else {
            for (slot, p) in pos.iter_mut().zip(world.positions.iter()) {
                *slot = p.to_array();
            }
        }
    }

    // Update colours
    if let Some(VertexAttributeValues::Float32x4(ref mut col)) =
        mesh.attribute_mut(Mesh::ATTRIBUTE_COLOR)
    {
        if col.len() != world.n {
            *col = world.neighbor_counts.iter().map(|&c| density_color(c)).collect();
        } else {
            for (slot, &c) in col.iter_mut().zip(world.neighbor_counts.iter()) {
                *slot = density_color(c);
            }
        }
    }
}

// ─── Density colour map ───────────────────────────────────────────────────────

/// Maps neighbour count to an RGBA colour.
/// 0 neighbours = deep blue (free particle).
/// High density  = white-hot (bound cluster).
pub fn density_color(n: u8) -> [f32; 4] {
    let t = (n as f32 / COLOR_SCALE_N).clamp(0.0, 1.0);
    // 5-stop gradient: blue → cyan → green → yellow → red → white
    let (r, g, b) = if t < 0.2 {
        let s = t / 0.2;
        (0.0,       s * 0.5,  1.0)
    } else if t < 0.4 {
        let s = (t - 0.2) / 0.2;
        (0.0,       0.5 + s * 0.5, 1.0 - s)
    } else if t < 0.6 {
        let s = (t - 0.4) / 0.2;
        (s,         1.0,       0.0)
    } else if t < 0.8 {
        let s = (t - 0.6) / 0.2;
        (1.0,       1.0 - s,  0.0)
    } else {
        let s = (t - 0.8) / 0.2;
        (1.0,       s,         s)    // red → white
    };
    [r, g, b, 1.0]
}
