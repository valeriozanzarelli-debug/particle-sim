//! Drawing of the zig-zag scenes (free photon, live atom) and the animated proton knot. Same soft-ball
//! shader as render.rs; meshes replaced, never modified in place.

use bevy::{
    pbr::{NotShadowCaster, NotShadowReceiver},
    prelude::*,
    render::view::NoFrustumCulling,
};

use std::sync::Arc;

use crate::camera::Cam;
use crate::istanze::{impacca, quadrato, StratoIstanze};
use crate::oit_palline::{vortice_disegnato, OitPalline};
use crate::render::{quads, quads_s, schermo, BallMaterial, BallParams};
use crate::{Lab, Scena};

#[derive(Component, Clone, Copy, PartialEq, Eq, Debug)]
pub enum ZzStrato { Palline, Nodo, Onde }

/// The layers drawn as GPU instances.
#[derive(Component, Clone, Copy, PartialEq, Eq, Debug)]
pub enum ZzStratoI { Palline, Onde }

/// Drawing path: GPU instances (default) or the old rebuilt meshes (`--mesh`, kept for the comparison R2).
#[derive(Resource, Clone, Copy)]
pub struct ModoDisegno {
    pub istanze: bool,
    /// the photon's energy balls with weighted blended OIT (oit_palline.rs); `--senza-oit` = the old alpha blending
    pub oit: bool,
}

#[derive(Resource)]
pub struct ZzVisuals { mats: [Handle<BallMaterial>; 3], versione: u64, scena: Option<Scena>, nodo_fatto: bool }

pub fn setup(mut commands: Commands, mut meshes: ResMut<Assets<Mesh>>, mut mats: ResMut<Assets<BallMaterial>>) {
    let p = BallParams { size: 0.1, gain: 1.0, soft: 3.0, core: 0.3, blend: 0.0 };
    let m0 = mats.add(BallMaterial { params: p, blend: false });
    let m1 = mats.add(BallMaterial { params: p, blend: false });
    let m2 = mats.add(BallMaterial { params: BallParams { blend: 1.0, ..p }, blend: true });
    // the blended wave balls first (drawn behind), the glowing layers over them
    for (s, m) in [(ZzStrato::Onde, m2.clone()), (ZzStrato::Palline, m0.clone()), (ZzStrato::Nodo, m1.clone())] {
        commands.spawn((s, Mesh3d(meshes.add(quads(&[]))), MeshMaterial3d(m), Transform::default(),
                        NoFrustumCulling, NotShadowCaster, NotShadowReceiver));
    }
    // the same two layers as GPU instances (exam R2 of docs/MOTORE_GRAFICO.md); used unless `--mesh`
    for (s, alfa) in [(ZzStratoI::Onde, true), (ZzStratoI::Palline, false)] {
        commands.spawn((s, Mesh3d(meshes.add(quadrato())), StratoIstanze { alfa, ..default() }, Transform::default(),
                        NoFrustumCulling, NotShadowCaster, NotShadowReceiver));
    }
    commands.insert_resource(ZzVisuals { mats: [m0, m1, m2], versione: u64::MAX, scena: None, nodo_fatto: false });
}

/// Hue from a direction (the pion isovector of the hedgehog points along r̂): a colour wheel in azimuth,
/// brighter towards the poles.
fn colore_direzione(p: [f64; 3], d: f32) -> [f32; 4] {
    let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt().max(1e-30);
    let az = p[1].atan2(p[0]) as f32;
    let cz = (p[2] / r) as f32;
    let h = |o: f32| 0.5 + 0.5 * (az + o).cos();
    let base = [h(0.0), h(-2.094), h(2.094)];
    let w = 0.35 + 0.65 * d;
    let mix = |c: f32| (c * (1.0 - 0.4 * cz.abs()) + 0.4 * cz.abs()) * w + (1.0 - w) * 0.9;
    [mix(base[0]).min(1.0), (mix(base[1]) * 0.85).min(1.0), mix(base[2]).min(1.0), 0.07 + 0.22 * d]
}

#[allow(clippy::too_many_arguments)]
pub fn aggiorna(
    time: Res<Time>,
    mut lab: ResMut<Lab>,
    cam: Res<Cam>,
    mut vis: ResMut<ZzVisuals>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut mats: ResMut<Assets<BallMaterial>>,
    mut q: Query<(&ZzStrato, &mut Mesh3d, &mut Transform)>,
    mut qi: Query<(&ZzStratoI, &mut StratoIstanze)>,
    modo: Res<ModoDisegno>,
    mut oit: ResMut<OitPalline>,
    q_cam: Query<&GlobalTransform, With<Camera3d>>,
) {
    let lab = &mut *lab;
    // GPU path of the photon (gpu_fotone.rs): the energy balls are already instances on the GPU
    let onde_gpu = if lab.scena == Scena::Fotone { lab.fot.as_ref().and_then(|c| c.lock().unwrap().out.onde_gpu.clone()) } else { None };
    let r_cam = cam.r();
    let focus = cam.focus;
    let scena = lab.scena;
    // what the main layer shows
    type Strati = (Vec<crate::render::Ball>, Vec<crate::render::Ball>, Vec<[f32; 3]>);
    let (items, versione): (Option<Strati>, u64) = match scena {
        Scena::Fotone => {
            let out = lab.fot.as_ref().map(|c| c.lock().unwrap().out.clone());
            match out {
                Some(o) if vis.versione != o.versione || vis.scena != Some(scena) => (Some(((*o.palline).clone(), (*o.onde).clone(), (*o.allunga).clone())), o.versione),
                Some(o) => (None, o.versione),
                None => (None, 0),
            }
        }
        Scena::Atomo => {
            let out = lab.zz.as_ref().map(|c| c.lock().unwrap().out.clone());
            match out {
                Some(o) if vis.versione != o.versione || vis.scena != Some(scena) => (Some(((*o.palline).clone(), (*o.onde).clone(), vec![])), o.versione),
                Some(o) => (None, o.versione),
                None => (None, 0),
            }
        }
        Scena::Orbitali => {
            let out = lab.orb.as_ref().map(|c| c.lock().unwrap().out.clone());
            match out {
                Some(o) if vis.versione != o.versione || vis.scena != Some(scena) => (Some(((*o.palline).clone(), (*o.onde).clone(), vec![])), o.versione),
                Some(o) => (None, o.versione),
                None => (None, 0),
            }
        }
    };
    let unita_fm = lab.metri_per_unita() / 1e-15;
    for (s, mut m3, mut tr) in q.iter_mut() {
        match s {
            ZzStrato::Palline => {
                tr.translation = schermo([-focus[0], -focus[1], -focus[2]]);
                if let Some((it, onde, _)) = &items {
                    let old = m3.0.id();
                    m3.0 = meshes.add(quads(if modo.istanze { &[] } else { it }));
                    meshes.remove(old);
                    lab.palline_disegnate = it.len() + onde.len() + onde_gpu.as_ref().map_or(0, |g| g.n as usize);
                }
            }
            ZzStrato::Onde => {
                tr.translation = schermo([-focus[0], -focus[1], -focus[2]]);
                if let Some((_, it, st)) = &items {
                    let old = m3.0.id();
                    let st0: Vec<[f32; 3]> = it.iter().zip(st.iter()).map(|(b, v)| if b.1 < 0.95 { [0.0; 3] } else { *v }).collect();
                    m3.0 = if modo.istanze { meshes.add(quads(&[])) } else { meshes.add(quads_s(it, &st0)) };
                    meshes.remove(old);
                }
            }
            ZzStrato::Nodo => {
                if !vis.nodo_fatto {
                    // balls drawn once in fm (random radii from the baryon density, random directions: no lattice)
                    let items: Vec<crate::render::Ball> = lab.nodo_palline.iter()
                        .map(|(p, d)| (schermo(*p).to_array(), 1.0, colore_direzione(*p, d.clamp(0.0, 1.0)))).collect();
                    let old = m3.0.id();
                    m3.0 = meshes.add(quads(&items));
                    meshes.remove(old);
                    vis.nodo_fatto = true;
                }
                // the knot lives at the nucleus; its own clock is slowed so that its motions can be seen:
                // it SPINS (for the hedgehog an isorotation is a rotation: the colours turn) and BREATHES
                let t = lab.t_nucleo;
                let ang = (lab.omega_rot * t) as f32;
                let resp = 1.0 + 0.07 * (lab.omega_resp * t).sin() as f32;
                let k = (1.0 / unita_fm) as f32; // fm → scene units
                *tr = Transform::from_translation(schermo([-focus[0], -focus[1], -focus[2]]))
                    .with_rotation(Quat::from_rotation_y(ang))
                    .with_scale(Vec3::splat(k * resp));
            }
        }
    }
    if modo.istanze {
        let off = schermo([-focus[0], -focus[1], -focus[2]]).to_array();
        let con_oit = modo.oit && scena == Scena::Fotone;
        for (s, mut si) in qi.iter_mut() {
            si.offset = off;
            si.gpu = if *s == ZzStratoI::Onde && !con_oit { onde_gpu.clone() } else { None };
            if let Some((it, onde, st)) = &items {
                // without OIT the tails are drawn round, as before (their stretch is read only by the OIT shader)
                let st_vecchio = || -> Vec<[f32; 3]> { onde.iter().zip(st.iter()).map(|(b, v)| if b.1 < 0.95 { [0.0; 3] } else { *v }).collect() };
                let dati = Arc::new(match s {
                    ZzStratoI::Palline => impacca(it, &[]),
                    ZzStratoI::Onde => if con_oit { impacca(onde, st) } else { impacca(onde, &st_vecchio()) },
                });
                // with OIT the energy balls go to the OIT layer, the instanced alpha layer stays empty
                if *s == ZzStratoI::Onde && con_oit { oit.dati = dati; si.dati = Arc::new(vec![]); } else { si.dati = dati; }
            }
        }
        oit.attivo = con_oit;
        if con_oit {
            // the aesthetic knobs can be tried without rebuilding (OIT_PROFONDITA, OIT_ESPOSIZIONE, OIT_CODA_LATO, OIT_GUADAGNO)
            let env = |k: &str, d: f32| std::env::var(k).ok().and_then(|v| v.parse::<f32>().ok()).unwrap_or(d);
            oit.profondita = env("OIT_PROFONDITA", crate::oit_palline::DEF_PROFONDITA);
            oit.esposizione = env("OIT_ESPOSIZIONE", crate::oit_palline::DEF_ESPOSIZIONE);
            oit.coda_lato = env("OIT_CODA_LATO", crate::oit_palline::DEF_CODA_LATO);
            oit.guadagno = env("OIT_GUADAGNO", crate::oit_palline::DEF_GUADAGNO);
            oit.gpu = onde_gpu.clone();
            oit.offset = off;
            if let Some(o) = lab.fot.as_ref().map(|c| c.lock().unwrap().out.clone()) {
                let par = o.par;
                let lam = par.cpl() as f32;
                let d = crate::fotone3d::norma(crate::fotscena::DIREZIONI[par.dir % crate::fotscena::DIREZIONI.len()].1);
                oit.asse = schermo(d).normalize_or_zero().to_array();
                oit.centro = (schermo(o.centro) + Vec3::from(off)).to_array();
                oit.lambda = lam;
                oit.w0 = (par.w0_lambda as f32) * lam;
                // energy ∝ |F|²: the length of the pulse in energy is σ/√2
                oit.sigma = (crate::fotscena::SIGMA_LAMBDA as f32) * lam / std::f32::consts::SQRT_2;
                let teste = crate::fotscena::PALLINE_QUAL[par.qual.min(crate::fotscena::PALLINE_QUAL.len() - 1)].max(crate::fotdisegno::N_SURF) as f32;
                let n_ist = oit.gpu.as_ref().map(|g| g.n as usize).unwrap_or(oit.dati.len()) as f32;
                oit.teste = teste;
                oit.ns = (n_ist / teste - 1.0).round().max(0.0);
                // exam V2 (docs/MOTORE_GRAFICO_OIT.md): with OIT_VORTICE=file.json, the swirl as drawn at the frozen tick
                let fermo = crate::fotscena::TICK_FERMO.load(std::sync::atomic::Ordering::Relaxed);
                if let (Ok(path), Some(gt)) = (std::env::var("OIT_VORTICE"), q_cam.iter().next()) {
                    if fermo > 0 && o.ticks == fermo && !oit.dati.is_empty() && oit.gpu.is_none() {
                        let m = gt.compute_matrix();
                        let (right, up, fwd, eye) = (m.x_axis.truncate().normalize(), m.y_axis.truncate().normalize(), -m.z_axis.truncate().normalize(), m.w_axis.truncate());
                        let ns = oit.ns as usize;
                        let (c, l3, n) = vortice_disegnato(&oit.dati, ns, &oit, eye, fwd, right, up);
                        let ll = Vec3::from(oit.asse).dot(fwd);
                        let j = format!("{{\"tick\": {}, \"polar\": \"{:?}\", \"qual\": {}, \"ns\": {}, \"teste\": {}, \"C_rad\": {:.6}, \"L_asse_3d\": {:.6e}, \"asse_dot_sguardo\": {:.4}, \"regola\": \"{}\"}}\n",
                            o.ticks, par.polar, par.qual, ns, n, c, l3, ll, o.regola);
                        let _ = std::fs::write(&path, j);
                    }
                }
            }
        }
    }
    vis.versione = versione;
    vis.scena = Some(scena);
    // sizes: screen-constant, capped
    let fm = 1e-15 / lab.metri_per_unita();
    let g_nodo = if lab.mostra_protone && scena != Scena::Fotone {
        let x = (r_cam / (300.0 * fm)).ln();
        (1.0 - x / 2.0).clamp(0.0, 1.0) as f32
    } else { 0.0 };
    let (s_main, g_main) = match scena {
        Scena::Fotone => ((0.0032 * r_cam).min(0.9), lab.luce),
        Scena::Atomo => ((0.0035 * r_cam).min(0.12), lab.luce * (1.0 - 0.8 * g_nodo)),
        // balls about half a grid step at the full view (the grid step is 2R/64), smaller when zooming in
        Scena::Orbitali => ((0.0045 * r_cam).min(0.018 * lab.orb_r), lab.luce * (1.0 - 0.8 * g_nodo)),
    };
    if let Some(m) = mats.get_mut(&vis.mats[0]) { m.params = BallParams { size: (s_main * lab.grandezza as f64) as f32, gain: g_main, soft: 3.0, core: 0.3, blend: 0.0 }; }
    if let Some(m) = mats.get_mut(&vis.mats[1]) { m.params = BallParams { size: ((0.004 * r_cam).min(0.06 * fm)) as f32, gain: 1.3 * g_nodo, soft: 3.0, core: 0.4, blend: 0.0 }; }
    // wave balls: a little smaller than the glowing layers, as in the photon laboratory
    if let Some(m) = mats.get_mut(&vis.mats[2]) { m.params = BallParams { size: (0.8 * s_main * lab.grandezza as f64) as f32, gain: lab.luce.min(1.5), soft: 2.0, core: 0.0, blend: 1.0 }; }
    if let Some(m) = mats.get(&vis.mats[2]) { oit.size = m.params.size; oit.gain = m.params.gain; oit.soft = m.params.soft; }
    if modo.istanze {
        for (s, mut si) in qi.iter_mut() {
            let k = match s { ZzStratoI::Palline => 0, ZzStratoI::Onde => 2 };
            if let Some(m) = mats.get(&vis.mats[k]) {
                let p = m.params;
                si.size = p.size; si.gain = p.gain; si.soft = p.soft; si.core = p.core;
            }
        }
    }
    let _ = time;
}
