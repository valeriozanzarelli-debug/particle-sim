//! LABORATORIO ATOMO — one carbon atom from our rules, from ångströms down to the proton.
//!
//! • the electron cloud: 6 electrons as balls of the electron-ball rule (fixed-node diffusion Monte
//!   Carlo per ball, carbonio.rs), guide from our own hydrogen-like orbitals; inputs = physical constants,
//!   Z = 6 and the nuclear mass;
//! • zoom into the nucleus: the proton as the B = 1 knot of the pion field (skyrme.rs, from its ODE);
//! • zoom onto an electron: the zig-zag Dirac walk (src/zigzag/walk.rs), with its Zitterbewegung;
//! • photons: absorbed only at energies our own DMC computed (2p → 3s, ionization), else they pass;
//! • reversible time (tempo.rs) and a temperature that moves the nucleus (Maxwell–Boltzmann).
//!
//! Viewer:   cargo run --bin atomo [-- --zoom atomo|nucleo|elettrone]
//! Measures: cargo run --bin atomo -- --misura [veloce|tutto]      (numbers -> assets/atomo/risultati.txt)
//! Tests:    cargo test --bin atomo        (long one: cargo test --bin atomo -- --ignored)

use bevy::{
    core_pipeline::{bloom::Bloom, tonemapping::Tonemapping},
    prelude::*,
    window::{PresentMode, WindowResolution},
};
use bevy_egui::{EguiContexts, EguiPlugin};

mod camera;
mod carbonio;
mod fenscena;
mod luce;
mod render_zz;
mod zz;
mod zzscena;
#[path = "../molecole/costanti.rs"]
mod costanti;
#[path = "../zigzag/cplx.rs"]
mod cplx;
mod misura;
mod mondo;
mod render;
mod risultati;
mod skyrme;
mod tempo;
mod ui;
#[path = "../zigzag/walk.rs"]
#[allow(dead_code)]
mod walk;
#[cfg(test)]
mod tests;

use camera::{Bersaglio, Cam, PointerOverHud};
use cplx::C;
use mondo::{Input, Livelli, Mondo};
use risultati::Risultati;
use tempo::Storia;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Modo { Avanti, Pausa, Indietro }

/// The three scenes. The first two are the zig-zag rule live (light; electron + light); the third is the
/// older electron-ball rule (DMC) for carbon, kept as an optional view.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scena { Fenditura, Atomo, Carbonio }

/// The nucleus is shown on its own slowed clock: spin angular velocity on screen (rad per screen second).
pub const OMEGA_ROT_SCHERMO: f64 = 1.1;

/// Zig-zag electron: cells per side, spacing (λ_C), packet width (λ_C), ticks before a restart.
pub const N_DIRAC: usize = 48;
pub const A_DIRAC: f64 = 0.25;
pub const SIGMA_DIRAC: f64 = 2.0;
pub const TICK_DIRAC_MAX: u64 = 100;

#[derive(Resource)]
pub struct Lab {
    pub scena: Scena,
    pub fen: Option<fenscena::ScenaFenditura>,
    pub zz: Option<std::sync::Arc<std::sync::Mutex<zzscena::Condiviso>>>,
    /// photon energy for the zig-zag atom (keV) and the worker pace (ticks per second)
    pub e_kev: f64,
    pub ritmo: f64,
    pub passi_frame: f64,
    acc_passi: f64,
    /// the nucleus clock (screen seconds) and its angular frequencies on screen; real ones in `nodo_info2`
    pub t_nucleo: f64,
    pub omega_rot: f64,
    pub omega_resp: f64,
    /// [Λ, M_Δ−M_N (MeV), ω_rot (1/s), ω_breath (1/s), ħω_breath (MeV)]
    pub nodo_info2: [f64; 5],
    pub r: Risultati,
    pub lv: Livelli,
    pub mondo: Mondo,
    pub storia: Storia,
    /// frame shown (== storia.ultimo() when live)
    pub vista: u64,
    pub modo: Modo,
    /// physical seconds per screen second
    pub velocita: f64,
    pub t_kelvin: f64,
    pub lambda_nm: f64,
    /// direction index of the next photon (0..6)
    pub dir_sparo: usize,
    pub spara: bool,
    pub panel_open: bool,
    pub mostra_protone: bool,
    pub mostra_elettrone: bool,
    pub mostra_nuvola: bool,
    pub bagliore: bool,
    pub grandezza: f32,
    pub luce: f32,
    /// shell boundaries found from the balls (bohr) and the accumulated radial histogram
    pub confini: Vec<f64>,
    pub istogramma: Vec<f64>,
    pub versione_colori: u64,
    /// the proton: balls (fm, density) and its numbers
    pub nodo_palline: Vec<([f64; 3], f32)>,
    pub nodo_info: [f64; 5],
    /// the zig-zag electron
    pub dirac: Option<walk::Dirac3D>,
    pub dirac_media: Vec<[f64; 3]>,
    pub versione_dirac: u64,
    pub zbw_ampiezza: f64,
    pub punto_elettrone: [f64; 3],
    pub frame_count: u64,
    pub palline_disegnate: usize,
    pub e_media: f64,
    pub fps: f64,
}

pub const DR_ISTO: f64 = 0.025;
pub const R_ISTO: f64 = 14.0;

impl Lab {
    pub fn new() -> Self {
        let r = Risultati::carica();
        let lv = Livelli::da(&r);
        let mondo = Mondo::new(&r, 2026);
        let storia = Storia::new(&mondo);
        let prof = skyrme::Profilo::calcola();
        let info = [prof.a, prof.numero_barionico(), prof.energia_12pi2(), prof.raggio_rms_fm(), prof.massa_mev()];
        let nodo_palline = prof.palline(40_000, 12);
        let (w_resp, e_resp) = prof.respiro();
        let w_rot = prof.omega_rot_s();
        let info2 = [prof.inerzia(), prof.delta_n_mev(), w_rot, w_resp, e_resp];
        let lambda_nm = lv.lambda_ecc_nm();
        let mut lab = Lab {
            scena: Scena::Fenditura, fen: Some(fenscena::ScenaFenditura::new()), zz: None, e_kev: 0.0, ritmo: 25.0,
            passi_frame: 1.0, acc_passi: 0.0, t_nucleo: 0.0, omega_rot: OMEGA_ROT_SCHERMO,
            omega_resp: OMEGA_ROT_SCHERMO * w_resp / w_rot, nodo_info2: info2,
            r, lv, mondo, storia, vista: 0, modo: Modo::Avanti, velocita: 3.0e-18, t_kelvin: 293.0, lambda_nm, dir_sparo: 0,
            spara: false, panel_open: false, mostra_protone: true, mostra_elettrone: false, mostra_nuvola: true, bagliore: true,
            grandezza: 1.0, luce: 1.0, confini: vec![], istogramma: vec![0.0; (R_ISTO / DR_ISTO) as usize],
            versione_colori: 0, nodo_palline, nodo_info: info, dirac: None, dirac_media: vec![], versione_dirac: 0,
            zbw_ampiezza: 0.0, punto_elettrone: [1.2, 0.0, 0.3], frame_count: 0, palline_disegnate: 0, e_media: 0.0, fps: 0.0,
        };
        lab.e_media = lab.mondo.dmc.e_mix;
        lab.aggiorna_confini(true);
        lab
    }

    pub fn dal_vivo(&self) -> bool { self.vista == self.storia.ultimo() }

    /// Metres per scene unit.
    pub fn metri_per_unita(&self) -> f64 {
        match self.scena {
            Scena::Carbonio => costanti::bohr_m(),
            Scena::Atomo => zz::LAMBDA_C_M,
            // the photon laboratory's 500 nm light: 12 cells per wavelength
            Scena::Fenditura => 500e-9 / 12.0,
        }
    }
    /// Seconds per tick of the scene's rule.
    pub fn secondi_per_tick(&self) -> f64 {
        match self.scena {
            Scena::Atomo => zz::A * zz::T_UNIT_S,
            Scena::Fenditura => 500e-9 / 12.0 / 299_792_458.0,
            Scena::Carbonio => 0.0,
        }
    }
    /// Seconds of the nucleus per screen second (its slowed clock).
    pub fn rallentamento_nucleo(&self) -> f64 { OMEGA_ROT_SCHERMO / self.nodo_info2[2] }

    pub fn cambia_scena(&mut self, s: Scena) {
        if s == self.scena && (s != Scena::Atomo || self.zz.is_some()) { return; }
        self.scena = s;
        match s {
            Scena::Fenditura => { if self.fen.is_none() { self.fen = Some(fenscena::ScenaFenditura::new()); } }
            Scena::Atomo => { if self.zz.is_none() { self.zz = Some(zzscena::avvia()); } }
            Scena::Carbonio => {}
        }
        if let Some(z) = &self.zz {
            let mut g = z.lock().unwrap();
            g.comando = if s == Scena::Atomo && self.modo == Modo::Avanti { zzscena::Comando::Avanti } else { zzscena::Comando::Pausa };
        }
    }
    pub fn punto_scena(&self, b: Bersaglio) -> [f64; 3] {
        match self.scena {
            Scena::Carbonio => self.punto_bersaglio(b),
            Scena::Fenditura => [0.0; 3],
            Scena::Atomo => match b {
                Bersaglio::Elettrone => self.zz.as_ref().map(|z| z.lock().unwrap().out.punto).unwrap_or([0.0; 3]),
                _ => [0.0; 3],
            },
        }
    }
    pub fn distanza_tipica(&self, b: Bersaglio) -> f64 {
        let u = self.metri_per_unita();
        match (self.scena, b) {
            (_, Bersaglio::Nucleo) => 3.0e-15 / u,
            (Scena::Carbonio, Bersaglio::Elettrone) => 30.0 * costanti::alpha_fs(),
            (Scena::Atomo, Bersaglio::Elettrone) => 3.0,
            (Scena::Carbonio, _) => 9.0,
            (Scena::Atomo, _) => 11.0,
            (Scena::Fenditura, _) => 230.0,
        }
    }

    /// Where the camera looks (bohr).
    pub fn punto_bersaglio(&self, b: Bersaglio) -> [f64; 3] {
        let n = self.storia.frame(self.vista).nucleo;
        match b {
            Bersaglio::Atomo | Bersaglio::Nucleo => n,
            Bersaglio::Elettrone => [n[0] + self.punto_elettrone[0], n[1] + self.punto_elettrone[1], n[2] + self.punto_elettrone[2]],
        }
    }

    /// Shell colours: boundaries from the radial histogram of the balls (accumulated over ticks).
    fn aggiorna_confini(&mut self, reset: bool) {
        let nb = self.istogramma.len();
        if reset { self.istogramma.iter_mut().for_each(|v| *v = 0.0); }
        for v in self.istogramma.iter_mut() { *v *= 0.97; }
        for c in self.mondo.dmc.x.chunks(3) {
            let r = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
            let i = (r / DR_ISTO) as usize;
            if i < nb { self.istogramma[i] += 1.0; }
        }
        let c = carbonio::confini(&self.istogramma, DR_ISTO);
        let changed = c.len() != self.confini.len() || c.iter().zip(&self.confini).any(|(a, b)| (a - b).abs() > 0.04);
        if changed {
            self.confini = c;
            self.versione_colori += 1;
        }
    }

    /// Pick the point where we look at "one real electron": one ball of the cloud, an L-shell electron.
    pub fn scegli_elettrone(&mut self) {
        let ne = self.mondo.dmc.g.ne();
        let n = self.mondo.dmc.n();
        let mut rng = carbonio::Rng::new(self.frame_count + 77);
        for _ in 0..1000 {
            let i = (rng.next_u64() % n as u64) as usize;
            let e = 2 + (rng.next_u64() % (ne as u64 - 2)) as usize;
            let p = self.mondo.dmc.elettrone(e).nth(i).unwrap();
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            if (0.8..2.0).contains(&r) { self.punto_elettrone = p; return; }
        }
    }

    pub fn reset_dirac(&mut self) {
        let n = N_DIRAC;
        let mut d = walk::Dirac3D::new(n, A_DIRAC, A_DIRAC);
        // zig-zag start: a purely right-handed spin-up packet (R₀ = 1): the mass flip turns it left-handed and
        // back at 2mc²/ħ, and ⟨r⟩ trembles along the spin axis (the Zitterbewegung)
        let mut norm = 0.0;
        for i in 0..n * n * n {
            let x = [d.x(i / (n * n)), d.x((i / n) % n), d.x(i % n)];
            let r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
            let g = (-r2 / (4.0 * SIGMA_DIRAC * SIGMA_DIRAC)).exp();
            d.psi[i] = [C::ZERO, C::ZERO, C::real(g), C::ZERO];
            norm += g * g;
        }
        let kk = 1.0 / norm.sqrt();
        for p in d.psi.iter_mut() { for c in p.iter_mut() { *c = c.scale(kk); } }
        self.dirac = Some(d);
        self.dirac_media.clear();
        self.versione_dirac += 1;
    }

    fn passo_dirac(&mut self) {
        if self.dirac.as_ref().map(|d| d.ticks >= TICK_DIRAC_MAX).unwrap_or(true) { self.reset_dirac(); }
        let d = self.dirac.as_mut().unwrap();
        d.tick();
        let n = d.n;
        let mut s = [0.0f64; 4];
        for (i, p) in d.psi.iter().enumerate() {
            let w: f64 = p.iter().map(|c| c.norm2()).sum();
            s[0] += w * d.x(i / (n * n));
            s[1] += w * d.x((i / n) % n);
            s[2] += w * d.x(i % n);
            s[3] += w;
        }
        self.dirac_media.push([s[0] / s[3], s[1] / s[3], s[2] / s[3]]);
        if self.dirac_media.len() > 30 {
            let z: Vec<f64> = self.dirac_media.iter().rev().take(30).map(|v| v[2]).collect();
            let (lo, hi) = z.iter().fold((f64::MAX, f64::MIN), |(a, b), &v| (a.min(v), b.max(v)));
            self.zbw_ampiezza = 0.5 * (hi - lo);
        }
        self.versione_dirac += 1;
    }

    /// Balls of the zig-zag electron (positions in bohr relative to its centre).
    pub fn palline_dirac(&self) -> Vec<render::Ball> {
        let Some(d) = self.dirac.as_ref() else { return vec![] };
        let n = d.n;
        let s = costanti::alpha_fs(); // bohr per λ_C
        let rl: Vec<(f32, f32)> = d.psi.iter().map(|p| ((p[0].norm2() + p[1].norm2()) as f32, (p[2].norm2() + p[3].norm2()) as f32)).collect();
        let rmax = rl.iter().map(|v| v.0 + v.1).fold(0.0, f32::max).max(1e-30);
        let mut out = Vec::new();
        for (i, &(l, r)) in rl.iter().enumerate() {
            let (a, b, c) = (i / (n * n), (i / n) % n, i % n);
            if (a + b + c) % 2 == 1 || l + r < 0.04 * rmax { continue; }
            // a random point inside the cell (fixed per cell): no lattice pattern on screen
            let mut hr = carbonio::Rng::new(i as u64 * 0x9E37 + 1);
            let p = [(d.x(a) + d.a * (hr.uni() - 0.5) * 2.0) * s, (d.x(b) + d.a * (hr.uni() - 0.5) * 2.0) * s, (d.x(c) + d.a * (hr.uni() - 0.5) * 2.0) * s];
            let h = r / (l + r);
            let col = [0.2 + 0.8 * h, 0.55, 1.0 - 0.88 * h];
            out.push((render::schermo(p).to_array(), 1.0, [col[0], col[1], col[2], 0.09 * ((l + r) / rmax).powf(0.5)]));
        }
        // ⟨r⟩: a trail of white balls (the trembling centre)
        for (j, m) in self.dirac_media.iter().enumerate() {
            let a = 0.3 + 0.7 * j as f32 / self.dirac_media.len().max(1) as f32;
            out.push((render::schermo([m[0] * s, m[1] * s, m[2] * s]).to_array(), 0.6, [1.0, 1.0, 1.0, a]));
        }
        out
    }

    /// One frame of the laboratory.
    fn avanza(&mut self, dt: f64) {
        self.frame_count += 1;
        let dt = dt.min(0.1);
        if dt > 0.0 { self.fps += (1.0 / dt - self.fps) * 0.05; }
        match self.modo {
            Modo::Avanti => self.t_nucleo += dt,
            Modo::Indietro => self.t_nucleo -= dt,
            Modo::Pausa => {}
        }
        match self.scena {
            Scena::Carbonio => self.avanza_carbonio(dt),
            Scena::Fenditura => {
                let Some(f) = self.fen.as_mut() else { return };
                self.acc_passi += self.passi_frame;
                while self.acc_passi >= 1.0 {
                    self.acc_passi -= 1.0;
                    match self.modo {
                        Modo::Avanti => f.avanti(),
                        Modo::Indietro => { if f.f.ticks == 0 { self.modo = Modo::Pausa; } else { f.indietro(); } }
                        Modo::Pausa => {}
                    }
                }
                self.spara = false;
            }
            Scena::Atomo => {
                if let Some(z) = &self.zz {
                    let mut g = z.lock().unwrap();
                    g.comando = match self.modo { Modo::Avanti => zzscena::Comando::Avanti, Modo::Indietro => zzscena::Comando::Indietro, Modo::Pausa => zzscena::Comando::Pausa };
                    g.ritmo = self.ritmo;
                    if self.spara {
                        self.spara = false;
                        if self.e_kev <= 0.0 { self.e_kev = g.out.omega_res * zz::MC2_KEV; }
                        g.spara = Some(self.e_kev / zz::MC2_KEV);
                    }
                }
            }
        }
    }

    fn avanza_carbonio(&mut self, dt: f64) {
        match self.modo {
            Modo::Pausa => {
                if self.spara { self.riprendi(); self.passo_vivo(0.0); }
            }
            Modo::Indietro => {
                self.spara = false;
                if self.vista > self.storia.inizio { self.vista -= 1; } else { self.modo = Modo::Pausa; }
            }
            Modo::Avanti => {
                if !self.dal_vivo() { self.riprendi(); }
                self.passo_vivo(dt);
                if self.mostra_elettrone { self.passo_dirac(); }
            }
        }
    }

    /// From a past frame: rebuild it exactly and drop the future.
    fn riprendi(&mut self) {
        if !self.dal_vivo() {
            self.mondo = self.storia.riprendi_da(self.vista, &self.r, &self.lv);
            self.vista = self.storia.ultimo();
        }
    }

    fn passo_vivo(&mut self, dt: f64) {
        let dirs: [[f64; 3]; 6] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.6, 0.6, 0.53], [-0.6, 0.48, -0.64]];
        let spara = if self.spara {
            self.spara = false;
            let d = dirs[self.dir_sparo % 6];
            let n = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            Some((self.lambda_nm, [d[0] / n, d[1] / n, d[2] / n]))
        } else { None };
        let inp = Input { dt_s: dt * self.velocita, dt_schermo: dt, ticks: if dt > 0.0 { 1 } else { 0 }, t_kelvin: self.t_kelvin, spara };
        let stato0 = self.mondo.p.stato;
        self.mondo.passo(&inp, &self.r, &self.lv);
        self.storia.registra(&self.mondo, inp);
        self.vista = self.storia.ultimo();
        let changed = self.mondo.p.stato != stato0;
        if changed { self.e_media = self.mondo.dmc.e_mix; }
        self.e_media += 0.01 * (self.mondo.dmc.e_mix - self.e_media);
        if changed || self.mondo.p.tick % 15 == 0 { self.aggiorna_confini(changed); }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if let Some(i) = args.iter().position(|a| a == "--misura") {
        misura::main(args.get(i + 1).map(|s| s.as_str()).unwrap_or("tutto"));
        return;
    }
    let zoom = args.iter().position(|a| a == "--zoom").and_then(|i| args.get(i + 1)).map(|s| s.as_str()).unwrap_or("atomo");
    let scena_arg = args.iter().position(|a| a == "--scena").and_then(|i| args.get(i + 1)).map(|s| s.as_str()).unwrap_or("");
    let mut lab = Lab::new();
    let scena = match (scena_arg, zoom) {
        ("fenditura", _) => Scena::Fenditura,
        ("carbonio", _) => Scena::Carbonio,
        ("atomo", _) | ("", "nucleo") | ("", "protone") | ("", "elettrone") => Scena::Atomo,
        _ => Scena::Fenditura,
    };
    lab.cambia_scena(scena);
    let arg = |k: &str| args.iter().position(|a| a == k).and_then(|i| args.get(i + 1)).cloned();
    let mut auto = Auto::default();
    if let Some(f) = arg("--foto") {
        auto.attivo = true;
        auto.foto = f;
        auto.frames = arg("--frames").and_then(|v| v.parse().ok()).unwrap_or(240);
        auto.spara_al = arg("--spara-al").and_then(|v| v.parse().ok());
        auto.lambda = arg("--lambda").and_then(|v| v.parse().ok());
        lab.panel_open = args.iter().any(|a| a == "--pannello");
    }
    if let Some(t) = arg("--temperatura").and_then(|v| v.parse().ok()) { lab.t_kelvin = t; }
    if let Some(v) = arg("--velocita").and_then(|v| v.parse().ok()) { lab.velocita = v; }
    let b = match zoom { "nucleo" | "protone" => Bersaglio::Nucleo, "elettrone" => Bersaglio::Elettrone, _ => Bersaglio::Atomo };
    if scena == Scena::Carbonio && b == Bersaglio::Elettrone {
        lab.mostra_elettrone = true;
        lab.scegli_elettrone();
        lab.reset_dirac();
    }
    let cam = Cam::new(lab.distanza_tipica(b), b);
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "Laboratorio Atomo — carbonio dalle palline".into(),
                resolution: WindowResolution::new(1600.0, 900.0),
                present_mode: PresentMode::AutoVsync,
                ..default()
            }),
            ..default()
        }))
        .add_plugins(EguiPlugin)
        .add_plugins(render::BallsPlugin)
        .insert_resource(ClearColor(Color::srgb(0.008, 0.01, 0.02)))
        .insert_resource(lab)
        .insert_resource(auto)
        .insert_resource(cam)
        .init_resource::<PointerOverHud>()
        .add_systems(Startup, (setup_scene, render_zz::setup))
        .add_systems(Update, (
            ui::draw_hud,
            hud_pointer.after(ui::draw_hud),
            run_lab.after(hud_pointer),
            camera::camera_sys.after(run_lab),
            render::aggiorna.after(camera::camera_sys),
            render_zz::aggiorna.after(camera::camera_sys),
            bagliore.after(run_lab),
            auto_sys.after(run_lab),
        ))
        .run();
}

fn setup_scene(mut commands: Commands) {
    commands.spawn((
        Camera3d::default(),
        Camera { hdr: true, ..default() },
        Tonemapping::TonyMcMapface,
        Bloom { intensity: 0.18, ..Bloom::NATURAL },
        Transform::from_xyz(0.0, 3.0, 9.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));
}

fn bagliore(lab: Res<Lab>, mut q: Query<&mut Bloom>) {
    for mut b in q.iter_mut() {
        let want = if lab.bagliore { 0.18 } else { 0.0 };
        if b.intensity != want { b.intensity = want; }
    }
}

fn hud_pointer(mut ctx: EguiContexts, mut over: ResMut<PointerOverHud>) {
    let c = ctx.ctx_mut();
    over.0 = c.wants_pointer_input() || c.is_pointer_over_area();
}

fn run_lab(time: Res<Time>, mut lab: ResMut<Lab>, auto: Res<Auto>) {
    // automatic runs use a fixed frame time, so the pictures are reproducible
    let dt = if auto.attivo { 1.0 / 60.0 } else { time.delta_secs_f64() };
    lab.avanza(dt);
}

/// Unattended run for checks: `--foto file.png [--frames N] [--spara-al F] [--lambda nm] [--pannello]`.
/// Mouse and keyboard are ignored; a screenshot is saved at frame N and the program exits.
#[derive(Resource, Default)]
pub struct Auto {
    pub attivo: bool,
    pub foto: String,
    pub frames: u64,
    pub spara_al: Option<u64>,
    pub lambda: Option<f64>,
    pub fatto: bool,
}

fn auto_sys(mut commands: Commands, mut auto: ResMut<Auto>, mut lab: ResMut<Lab>, mut exit: EventWriter<AppExit>) {
    if !auto.attivo { return; }
    let f = lab.frame_count;
    if Some(f) == auto.spara_al {
        if let Some(l) = auto.lambda { lab.lambda_nm = l; }
        lab.spara = true;
    }
    if f == auto.frames && !auto.fatto {
        auto.fatto = true;
        let path = auto.foto.clone();
        commands.spawn(bevy::render::view::screenshot::Screenshot::primary_window())
            .observe(bevy::render::view::screenshot::save_to_disk(path));
        let zz = lab.zz.as_ref().map(|z| { let o = z.lock().unwrap().out.clone(); format!("tick {} p1s {:.3} dip {:+.3} {:.0} ms/tick", o.ticks, o.p1s, o.dipolo, o.ms_tick) }).unwrap_or_default();
        eprintln!("foto al frame {f}: scena {:?}, {} palline, {:.0} fps {zz}", lab.scena, lab.palline_disegnate, lab.fps);
    }
    if f >= auto.frames + 20 { exit.send(AppExit::Success); }
}
