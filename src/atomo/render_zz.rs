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
use crate::render::{quads, quads_s, schermo, BallMaterial, BallParams};
use crate::{Lab, Scena};

#[derive(Component, Clone, Copy, PartialEq, Eq, Debug)]
pub enum ZzStrato { Palline, Nodo, Onde }

/// The layers drawn as GPU instances.
#[derive(Component, Clone, Copy, PartialEq, Eq, Debug)]
pub enum ZzStratoI { Palline, Onde, Scie }

/// Drawing path: GPU instances (default) or the old rebuilt meshes (`--mesh`, kept for the comparison R2).
#[derive(Resource, Clone, Copy)]
pub struct ModoDisegno { pub istanze: bool }

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
    for (s, alfa) in [(ZzStratoI::Onde, true), (ZzStratoI::Scie, false), (ZzStratoI::Palline, false)] {
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
    graf: Res<crate::grafica::Grafica>,
    budget: Res<crate::grafica::Budget>,
) {
    let lab = &mut *lab;
    // GPU path of the photon (gpu_fotone.rs): the energy balls are already instances on the GPU
    let onde_gpu = if lab.scena == Scena::Fotone { lab.fot.as_ref().and_then(|c| c.lock().unwrap().out.onde_gpu.clone()) } else { None };
    let r_cam = cam.r();
    let focus = cam.focus;
    let bpp = cam.bohr_per_px();
    // the eye seen from the focus, physics coordinates (camera.rs: r·(sin yaw cos pitch, sin pitch, cos yaw cos pitch) on
    // screen; screen (x, y, z) = physics (x, z, −y))
    // the light behind the cut fades over at most 0.6 × the camera distance (close up only a thin layer stays: no haze)
    let taglio_prof = graf.taglio_prof.min(0.6 * r_cam as f32).max(4.0);
    let occhio_dir = { let (yw, pt) = (cam.yaw as f64, cam.pitch as f64); [yw.sin() * pt.cos(), -(yw.cos() * pt.cos()), pt.sin()] };
    // the atom scene's level of detail (grafica.rs): its surfers are one group of glows
    let (lod_atomo, scie_atomo, n_testa) = if lab.scena == Scena::Atomo {
        let o = lab.zz.as_ref().map(|c| { let g = c.lock().unwrap(); (g.out.r_nuvola, g.out.passo_medio, g.out.n_testa, g.out.scie.clone(), g.out.scie_allunga.clone(), g.out.versione) });
        match o {
            Some((r, pm, nt, sc, sa, v)) => {
                let k = ((graf.scia_s as f64) * 25.0).round().min(15.0);
                // the scene graph of today: one node, the atom, whose glows are its surfers
                let nodo = crate::grafica::Nodo { nome: "atomo", centro: [0.0; 3], figli: vec![],
                    // the disc that holds most of the cloud: 1.5 × its rms radius
                    gruppo: crate::grafica::Gruppo { n: nt.max(1), raggio: 1.5 * r.max(0.1), scia: pm * k } };
                let mut dec = Vec::with_capacity(1);
                crate::grafica::percorri(&nodo, bpp, 3.7 * lab.grandezza as f64, graf.lod as f64 * budget.fattore, 1.0, &mut dec);
                let l = dec[0].lod;
                if std::env::var("ATOMO_DEBUG_LOD").is_ok() && lab.frame_count % 60 == 0 {
                    eprintln!("lod: r_cam {r_cam:.1} bpp {bpp:.4} n {nt} r_nuvola {r:.2} passo {pm:.3} scie {} -> {l:?}", sc.len());
                }
                (Some(l), Some((sc, sa, v)), nt)
            }
            None => (None, None, 0),
        }
    } else { (None, None, 0) };
    if lab.scena == Scena::Atomo {
        if let Some(c) = &lab.zz {
            c.lock().unwrap().vista = crate::grafica::VistaWorker {
                occhio: occhio_dir, r_cam, bpp, lod_f: lod_atomo.map(|l| l.f).unwrap_or(1.0), scia_s: graf.scia_s as f64,
                taglio: graf.taglio, taglio_prof: taglio_prof as f64, fetta: graf.fetta, zoom_atomo: graf.zoom_atomo, budget: graf.budget_luce(),
            };
        }
    }
    // the photon's worker sorts its glows by distance from the eye (physics coordinates, cells): camera.rs puts the
    // eye at r·(sin yaw cos pitch, sin pitch, cos yaw cos pitch) on screen, relative to the focus
    if lab.scena == Scena::Fotone {
        let (yw, pt) = (cam.yaw as f64, cam.pitch as f64);
        let e = [yw.sin() * pt.cos() * r_cam, pt.sin() * r_cam, yw.cos() * pt.cos() * r_cam];
        let occhio = [e[0] + focus[0], -e[2] + focus[1], e[1] + focus[2]];
        if let Some(c) = &lab.fot { c.lock().unwrap().occhio = Some(occhio); }
    }
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
        Scena::Tavola => {
            let out = lab.tav.out();
            match out {
                Some(o) if vis.versione != o.versione || vis.scena != Some(scena) => (Some(((*o.palline).clone(), (*o.onde).clone(), vec![])), o.versione),
                Some(o) => (None, o.versione),
                None => (None, 0),
            }
        }
        Scena::Elettrone => {
            let out = lab.el.as_ref().map(|c| c.lock().unwrap().out.clone());
            match out {
                Some(o) if vis.versione != o.versione || vis.scena != Some(scena) => (Some(((*o.palline).clone(), vec![], vec![])), o.versione),
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
                    m3.0 = if modo.istanze { meshes.add(quads(&[])) } else { meshes.add(quads_s(it, st)) };
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
        for (s, mut si) in qi.iter_mut() {
            si.offset = off;
            si.gpu = if *s == ZzStratoI::Onde { onde_gpu.clone() } else { None };
            if let Some((it, onde, st)) = &items {
                match s {
                    ZzStratoI::Palline => si.dati = Arc::new(impacca(it, &[])),
                    ZzStratoI::Onde => si.dati = Arc::new(impacca(onde, st)),
                    ZzStratoI::Scie => si.dati = Arc::new(match &scie_atomo { Some((sc, sa, _)) => impacca(sc, sa), None => vec![] }),
                }
            }
        }
    }
    vis.versione = versione;
    lab.versione_disegnata = versione;
    vis.scena = Some(scena);
    // sizes: screen-constant, capped
    let fm = 1e-15 / lab.metri_per_unita();
    let g_nodo = if lab.mostra_protone && scena != Scena::Fotone && scena != Scena::Elettrone && scena != Scena::Tavola {
        let x = (r_cam / (300.0 * fm)).ln();
        (1.0 - x / 2.0).clamp(0.0, 1.0) as f32
    } else { 0.0 };
    let (s_main, g_main) = match scena {
        Scena::Fotone => ((0.0032 * r_cam).min(0.9), lab.luce),
        // the atom: the surfer glows keep ~3.7 px at every distance (grafica::lod thins and grows them when far)
        Scena::Atomo => (3.7 * bpp, lab.luce * (1.0 - 0.8 * g_nodo)),
        Scena::Elettrone => ((0.0035 * r_cam).min(0.12), lab.luce),
        // the periodic-table scene is in bohr: glows ~1/30 of the view, at most 0.08 a0 (drawing choice)
        Scena::Tavola => ((0.004 * r_cam).min(0.08 * lab.tav.raggio_cam.max(1.0)), lab.luce),
        // balls about half a grid step at the full view (the grid step is 2R/64), smaller when zooming in
        Scena::Orbitali => ((0.0045 * r_cam).min(0.018 * lab.orb_r), lab.luce * (1.0 - 0.8 * g_nodo)),
    };
    // the atom's surfers: crisper (tighter glow, hotter core)
    let (soft, core) = if scena == Scena::Atomo { (4.0, 0.9) } else { (3.0, 0.3) };
    if let Some(m) = mats.get_mut(&vis.mats[0]) { m.params = BallParams { size: (s_main * lab.grandezza as f64) as f32, gain: g_main, soft, core, blend: 0.0 }; }
    if let Some(m) = mats.get_mut(&vis.mats[1]) { m.params = BallParams { size: ((0.004 * r_cam).min(0.06 * fm)) as f32, gain: 1.3 * g_nodo, soft: 3.0, core: 0.4, blend: 0.0 }; }
    // wave balls: a little smaller than the glowing layers, as in the photon laboratory
    // (atom scene: fixed in world units, S_ONDE, so the light glows belong to the light's cells at every zoom)
    let s_onde = if scena == Scena::Atomo { crate::zzscena::S_ONDE } else { 0.8 * s_main };
    if let Some(m) = mats.get_mut(&vis.mats[2]) { m.params = BallParams { size: (s_onde * lab.grandezza as f64) as f32, gain: lab.luce.min(1.5), soft: 2.0, core: 0.0, blend: 1.0 }; }
    if modo.istanze {
        for (s, mut si) in qi.iter_mut() {
            let k = match s { ZzStratoI::Palline | ZzStratoI::Scie => 0, ZzStratoI::Onde => 2 };
            if let Some(m) = mats.get(&vis.mats[k]) {
                let p = m.params;
                si.size = p.size; si.gain = p.gain; si.soft = p.soft; si.core = p.core;
            }
            // level of detail, depth and the cut (istanze.wgsl); off in the other scenes
            let px_mondo = (2.0 * (cam.fov as f64 * 0.5).tan() / cam.altezza_px as f64) as f32;
            si.n_lod = 0; si.lod = [1.0, 1.0, 1.0, 0.0]; si.prof = [0.0; 4]; si.taglio = [0.0; 4];
            if let Some(l) = lod_atomo {
                let prof = [r_cam as f32, 2.0 * zz_r_nuvola(lab) as f32, graf.profondita, px_mondo];
                match s {
                    ZzStratoI::Palline => {
                        si.n_lod = n_testa as u32;
                        si.lod = [l.f, l.alfa, l.cresci, 1.6];
                        si.prof = prof;
                        si.gain *= graf.bagliori;
                    }
                    ZzStratoI::Scie => {
                        // already thinned by the worker with the same keys; same growth and brightness as the heads
                        si.n_lod = u32::MAX;
                        si.lod = [1.0, l.alfa, l.cresci, 1.2];
                        si.prof = prof;
                        si.gain *= graf.bagliori;
                        // a tail is a smooth light: no hot core in each piece (it would make beads), a broad profile
                        si.core = 0.0;
                        si.soft = 2.0;
                    }
                    ZzStratoI::Onde => {
                        si.lod = [1.0, 1.0, 1.0, 1.5];
                        si.n_lod = u32::MAX;
                        si.prof = [r_cam as f32, (r_cam as f32).max(60.0), 0.5 * graf.profondita, px_mondo];
                        si.taglio = [0.0, 1.0, 1.0, 14.0];
                        if graf.taglio && !graf.fetta { si.taglio = [1.0, crate::zzscena::W_TAGLIO as f32, taglio_prof, 14.0]; }
                    }
                }
            }
        }
    }
    let _ = time;
}

fn zz_r_nuvola(lab: &Lab) -> f64 { lab.zz.as_ref().map(|c| c.lock().unwrap().out.r_nuvola).unwrap_or(3.0).max(0.5) }
