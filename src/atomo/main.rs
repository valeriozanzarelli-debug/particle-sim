//! LABORATORIO ATOMO — light and matter from ONE rule: the zig-zag walk (src/zigzag/walk.rs).
//!
//! • scene "fotone": one free photon in a 3-D box, made of spin-1 balls (F = E + icB) — fotone3d.rs,
//!   fotscena.rs: Gaussian pulse, chosen polarization and wavelength, helicity, diffraction, wavefronts;
//! • scene "atomo": the one-electron ion Nd⁵⁹⁺ live with the same rule plus mass and the Coulomb phase (zz.rs,
//!   zzscena.rs); zoom into the nucleus: the proton as the B = 1 knot of the pion field (skyrme.rs).
//!
//! Viewer: cargo run --bin atomo [-- --scena fotone|atomo] [--zoom nucleo]
//! Tests:  cargo test --bin atomo        (long ones: cargo test --bin atomo -- --ignored)

use bevy::{
    core_pipeline::{bloom::Bloom, tonemapping::Tonemapping},
    prelude::*,
    window::{PresentMode, WindowResolution},
};
use bevy_egui::{EguiContexts, EguiPlugin};

mod accoppiato;
mod camera;
#[path = "../zigzag/cplx.rs"]
#[allow(dead_code)]
mod cplx;
mod fotdisegno;
mod fotone3d;
mod istanze;
mod oit_palline;
mod fotscena;
#[allow(dead_code)]
mod gpu_regola;
#[allow(dead_code)]
mod gpu_fotone;
mod misura;
mod orbitali;
mod render;
mod render_zz;
mod rng;
mod skyrme;
mod ui;
#[path = "../zigzag/walk.rs"]
#[allow(dead_code)]
mod walk;
mod zz;
mod zzscena;
#[cfg(test)]
mod tests;
#[cfg(test)]
mod gpu_esami;
#[cfg(test)]
mod test_battimento;

use camera::{Bersaglio, Cam, PointerOverHud};
use fotone3d::Polar;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Modo { Avanti, Pausa, Indietro }

/// The two scenes, both the zig-zag rule live: the free photon (light alone) and the atom (electron + light).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scena { Fotone, Atomo, Orbitali }

/// The nucleus is shown on its own slowed clock: spin angular velocity on screen (rad per screen second).
pub const OMEGA_ROT_SCHERMO: f64 = 1.1;

#[derive(Resource)]
pub struct Lab {
    pub scena: Scena,
    pub fot: Option<std::sync::Arc<std::sync::Mutex<fotscena::Condiviso>>>,
    /// the photon the panel will launch
    pub fot_par: fotscena::Parametri,
    pub ritmo_fot: f64,
    /// the camera follows the photon
    pub segui: bool,
    pub vista_fot: fotscena::Vista,
    /// the real wavelength (m): only the size of a cell depends on it (the rule has no length of its own)
    pub lambda_m: f64,
    pub zz: Option<std::sync::Arc<std::sync::Mutex<zzscena::Condiviso>>>,
    /// the orbitals (n = 1..4, s p d f, superpositions) from the radial form of the rule
    pub orb: Option<std::sync::Arc<std::sync::Mutex<orbitali::Condiviso>>>,
    /// what the orbital panel will show
    pub orb_scelta: orbitali::Scelta,
    pub ritmo_orb: f64,
    /// measurement and surfer options (all scenes), a measurement asked, clear the clicks
    pub opz: misura::Opzioni,
    /// atom scene: light coupled live (Maxwell–Dirac) instead of the prescribed pulse; photon energy set for it
    pub accoppiato: bool,
    pub e_kev_acc: bool,
    pub misura: bool,
    pub azzera_clic: bool,
    /// radius of the orbital view (λ_C) the camera was last set for
    pub orb_r: f64,
    /// photon energy for the zig-zag atom (keV) and the worker pace (ticks per second)
    pub e_kev: f64,
    pub ritmo: f64,
    /// the nucleus clock (screen seconds) and its angular frequencies on screen; real ones in `nodo_info2`
    pub t_nucleo: f64,
    pub omega_rot: f64,
    pub omega_resp: f64,
    /// [Λ, M_Δ−M_N (MeV), ω_rot (1/s), ω_breath (1/s), ħω_breath (MeV)]
    pub nodo_info2: [f64; 5],
    pub modo: Modo,
    pub spara: bool,
    pub panel_open: bool,
    pub mostra_protone: bool,
    pub bagliore: bool,
    pub grandezza: f32,
    pub luce: f32,
    /// the proton: balls (fm, density) and its numbers
    pub nodo_palline: Vec<([f64; 3], f32)>,
    pub nodo_info: [f64; 5],
    pub frame_count: u64,
    pub palline_disegnate: usize,
    pub fps: f64,
}

impl Lab {
    pub fn new(par: fotscena::Parametri) -> Self {
        let prof = skyrme::Profilo::calcola();
        let info = [prof.a, prof.numero_barionico(), prof.energia_12pi2(), prof.raggio_rms_fm(), prof.massa_mev()];
        let nodo_palline = prof.palline(40_000, 12);
        let (w_resp, e_resp) = prof.respiro();
        let w_rot = prof.omega_rot_s();
        let info2 = [prof.inerzia(), prof.delta_n_mev(), w_rot, w_resp, e_resp];
        Lab {
            scena: Scena::Fotone, fot: None, fot_par: par, ritmo_fot: 20.0, segui: false, vista_fot: fotscena::Vista::default(), lambda_m: 500e-9, zz: None, orb: None, orb_scelta: orbitali::Scelta::default(), ritmo_orb: 30.0, opz: misura::Opzioni::default(), accoppiato: true, e_kev_acc: false, misura: false, azzera_clic: false, orb_r: 25.0, e_kev: 0.0, ritmo: 25.0,
            t_nucleo: 0.0, omega_rot: OMEGA_ROT_SCHERMO, omega_resp: OMEGA_ROT_SCHERMO * w_resp / w_rot, nodo_info2: info2,
            modo: Modo::Avanti, spara: false, panel_open: false, mostra_protone: true, bagliore: true, grandezza: 1.0, luce: 1.0,
            nodo_palline, nodo_info: info, frame_count: 0, palline_disegnate: 0, fps: 0.0,
        }
    }

    /// Metres per scene unit.
    pub fn metri_per_unita(&self) -> f64 {
        match self.scena {
            Scena::Atomo | Scena::Orbitali => zz::LAMBDA_C_M,
            Scena::Fotone => self.cella_m(),
        }
    }
    /// Seconds per tick of the scene's rule.
    pub fn secondi_per_tick(&self) -> f64 {
        match self.scena {
            Scena::Atomo => zz::A * zz::T_UNIT_S,
            Scena::Orbitali => zz::T_UNIT_S,
            Scena::Fotone => self.cella_m() / 299_792_458.0,
        }
    }
    /// Size of one cell of the photon box (m): λ / (cells per λ) of the photon now flying.
    pub fn cella_m(&self) -> f64 {
        let cpl = self.fot.as_ref().map(|f| f.lock().unwrap().out.par.cpl()).unwrap_or(self.fot_par.cpl());
        self.lambda_m / cpl
    }
    /// Seconds of the nucleus per screen second (its slowed clock).
    pub fn rallentamento_nucleo(&self) -> f64 { OMEGA_ROT_SCHERMO / self.nodo_info2[2] }

    pub fn cambia_scena(&mut self, s: Scena) {
        self.scena = s;
        match s {
            Scena::Fotone => { if self.fot.is_none() { self.fot = Some(fotscena::avvia(self.fot_par)); } }
            Scena::Atomo => { if self.zz.is_none() { self.zz = Some(zzscena::avvia()); } }
            Scena::Orbitali => { if self.orb.is_none() { self.orb = Some(orbitali::avvia()); } }
        }
        if let Some(o) = &self.orb {
            let mut g = o.lock().unwrap();
            g.comando = if s == Scena::Orbitali && self.modo == Modo::Avanti { orbitali::Comando::Avanti } else { orbitali::Comando::Pausa };
        }
        // the hidden scene's worker rests
        if let Some(z) = &self.zz {
            let mut g = z.lock().unwrap();
            g.comando = if s == Scena::Atomo && self.modo == Modo::Avanti { zzscena::Comando::Avanti } else { zzscena::Comando::Pausa };
        }
        if let Some(f) = &self.fot {
            let mut g = f.lock().unwrap();
            g.comando = if s == Scena::Fotone && self.modo == Modo::Avanti { fotscena::Comando::Avanti } else { fotscena::Comando::Pausa };
        }
    }
    pub fn punto_scena(&self, b: Bersaglio) -> [f64; 3] {
        match self.scena {
            Scena::Fotone => {
                if self.segui { self.fot.as_ref().map(|f| f.lock().unwrap().out.centro).unwrap_or([0.0; 3]) } else { [0.0; 3] }
            }
            Scena::Atomo => match b {
                Bersaglio::Elettrone => self.zz.as_ref().map(|z| z.lock().unwrap().out.punto).unwrap_or([0.0; 3]),
                _ => [0.0; 3],
            },
            Scena::Orbitali => [0.0; 3],
        }
    }
    pub fn distanza_tipica(&self, b: Bersaglio) -> f64 {
        let u = self.metri_per_unita();
        match (self.scena, b) {
            (_, Bersaglio::Nucleo) => 3.0e-15 / u,
            (Scena::Atomo, Bersaglio::Elettrone) => 3.0,
            (Scena::Atomo, _) => 11.0,
            (Scena::Orbitali, Bersaglio::Elettrone) => 6.0,
            (Scena::Orbitali, _) => 2.4 * self.orb.as_ref().map(|o| o.lock().unwrap().out.r).unwrap_or(25.0),
            (Scena::Fotone, _) => 1.0 * self.fot.as_ref().map(|f| f.lock().unwrap().out.par.n()).unwrap_or(self.fot_par.n()) as f64,
        }
    }

    /// Launch a new photon with the panel's parameters.
    pub fn lancia_fotone(&mut self) {
        if let Some(f) = &self.fot { f.lock().unwrap().lancia = Some(self.fot_par); }
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
        let (misura, azzera) = (self.misura, self.azzera_clic);
        match self.scena {
            Scena::Fotone => {
                if self.spara { self.spara = false; self.lancia_fotone(); }
                if let Some(f) = &self.fot {
                    let mut g = f.lock().unwrap();
                    g.comando = match self.modo { Modo::Avanti => fotscena::Comando::Avanti, Modo::Indietro => fotscena::Comando::Indietro, Modo::Pausa => fotscena::Comando::Pausa };
                    g.ritmo = self.ritmo_fot;
                    g.vista = self.vista_fot;
                    g.opz = self.opz;
                    if self.misura { g.misura = true; }
                    if self.azzera_clic { g.azzera = true; }
                }
            }
            Scena::Atomo => {
                if let Some(z) = &self.zz {
                    let mut g = z.lock().unwrap();
                    g.comando = match self.modo { Modo::Avanti => zzscena::Comando::Avanti, Modo::Indietro => zzscena::Comando::Indietro, Modo::Pausa => zzscena::Comando::Pausa };
                    g.ritmo = self.ritmo;
                    g.opz = self.opz;
                    g.accoppiato = self.accoppiato;
                    if self.misura { g.misura = true; }
                    if self.azzera_clic { g.azzera = true; }
                    if self.spara {
                        self.spara = false;
                        if self.e_kev <= 0.0 { self.e_kev = g.out.omega_res * zz::MC2_KEV; }
                        g.spara = Some(self.e_kev / zz::MC2_KEV);
                    }
                }
            }
            Scena::Orbitali => {
                self.spara = false;
                if let Some(o) = &self.orb {
                    let mut g = o.lock().unwrap();
                    g.comando = match self.modo { Modo::Avanti => orbitali::Comando::Avanti, Modo::Indietro => orbitali::Comando::Indietro, Modo::Pausa => orbitali::Comando::Pausa };
                    g.ritmo = self.ritmo_orb;
                    g.opz = self.opz;
                    if self.misura { g.misura = true; }
                    if self.azzera_clic { g.azzera = true; }
                }
            }
        }
        let _ = (misura, azzera);
        self.misura = false;
        self.azzera_clic = false;
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let arg = |k: &str| args.iter().position(|a| a == k).and_then(|i| args.get(i + 1)).cloned();
    let zoom = arg("--zoom").unwrap_or_else(|| "intero".into());
    let scena_arg = arg("--scena").unwrap_or_default();
    // the photon's parameters from the command line (for unattended pictures)
    let mut par = fotscena::Parametri::default();
    if let Some(v) = arg("--qual").and_then(|v| v.parse().ok()) { par.qual = v; }
    if let Some(q) = arg("--qualita").and_then(|v| fotscena::preset(&v)) { par.qual = q; }
    // where the light rule runs: gpu (default for alta/ultra) or cpu (default for bassa/media)
    if let Some(r) = arg("--regola") { gpu_fotone::REGOLA.store(if r == "gpu" { 1 } else if r == "cpu" { 2 } else { 0 }, std::sync::atomic::Ordering::Relaxed); }
    if let Some(v) = arg("--w0").and_then(|v| v.parse().ok()) { par.w0_lambda = v; }
    if let Some(v) = arg("--dir").and_then(|v| v.parse().ok()) { par.dir = v; }
    if let Some(p) = arg("--polar") {
        par.polar = match p.as_str() { "h" => Polar::Orizzontale, "v" => Polar::Verticale, "c-" => Polar::CircolareMeno, _ => Polar::CircolarePiu };
    }
    let mut lab = Lab::new(par);
    if let Some(v) = arg("--stati") { if let Some(s) = orbitali::scelta_da(&v) { lab.orb_scelta = s; } }
    if let Some(v) = arg("--lambda-m").and_then(|v| v.parse().ok()) { lab.lambda_m = v; }
    if let Some(v) = arg("--strati") { lab.vista_fot = strati(&v); }
    if args.iter().any(|a| a == "--surfisti") { lab.opz.surfisti = true; }
    if args.iter().any(|a| a == "--collassa") { lab.opz.collassa = true; }
    let scena = match (scena_arg.as_str(), zoom.as_str()) {
        ("atomo", _) | ("", "nucleo") | ("", "protone") | ("", "elettrone") => Scena::Atomo,
        ("orbitali", _) => Scena::Orbitali,
        _ => Scena::Fotone,
    };
    lab.cambia_scena(scena);
    if scena == Scena::Orbitali { if let Some(o) = &lab.orb { o.lock().unwrap().scegli = Some(lab.orb_scelta.clone()); } }
    let mut auto = Auto::default();
    if let Some(f) = arg("--foto") {
        auto.attivo = true;
        auto.foto = f;
        auto.frames = arg("--frames").and_then(|v| v.parse().ok()).unwrap_or(240);
        auto.spara_al = arg("--spara-al").and_then(|v| v.parse().ok());
        auto.misure = arg("--misure").and_then(|v| v.parse().ok()).unwrap_or(0);
        lab.panel_open = args.iter().any(|a| a == "--pannello");
        if let Some(r) = arg("--ritmo").and_then(|v| v.parse().ok()) { lab.ritmo_fot = r; }
        if let Some(t) = arg("--tick").and_then(|v| v.parse::<u64>().ok()) { fotscena::TICK_FERMO.store(t, std::sync::atomic::Ordering::Relaxed); }
    }
    if let Some(y) = arg("--yaw").and_then(|v| v.parse::<f32>().ok()) { auto.yaw = Some(y); }
    if let Some(p) = arg("--pitch").and_then(|v| v.parse::<f32>().ok()) { auto.pitch = Some(p); }
    let b = match zoom.as_str() { "nucleo" | "protone" => Bersaglio::Nucleo, "elettrone" => Bersaglio::Elettrone, _ => Bersaglio::Atomo };
    let mut cam = Cam::new(lab.distanza_tipica(b), b);
    if scena == Scena::Fotone { cam.yaw = 0.42; cam.pitch = 0.28; }
    if let Some(y) = auto.yaw { cam.yaw = y; }
    if let Some(p) = auto.pitch { cam.pitch = p; }
    if let Some(d) = arg("--distanza").and_then(|v| v.parse::<f64>().ok()) { cam.ln_r = d.ln(); cam.ln_r_target = d.ln(); }
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "Laboratorio Atomo — la regola zig-zag".into(),
                resolution: WindowResolution::new(1600.0, 900.0),
                present_mode: PresentMode::AutoVsync,
                ..default()
            }),
            ..default()
        }).set(bevy::render::RenderPlugin {
            // the discrete GPU when there is one (GTX 1650 Ti on this laptop, not the Intel UHD)
            render_creation: bevy::render::settings::WgpuSettings {
                power_preference: bevy::render::settings::PowerPreference::HighPerformance,
                ..default()
            }.into(),
            ..default()
        }))
        .add_plugins(EguiPlugin)
        .add_plugins(render::BallsPlugin)
        .add_plugins(istanze::IstanzePlugin)
        .add_plugins(oit_palline::OitPlugin)
        .insert_resource(render_zz::ModoDisegno { istanze: !args.iter().any(|a| a == "--mesh"), oit: !args.iter().any(|a| a == "--senza-oit") })
        .insert_resource(ClearColor(Color::srgb(0.008, 0.01, 0.02)))
        .insert_resource(lab)
        .insert_resource(auto)
        .insert_resource(cam)
        .init_resource::<PointerOverHud>()
        .add_systems(Startup, (setup_scene, render_zz::setup, gpu_fotone::prendi_dispositivo))
        .add_systems(Update, (
            ui::draw_hud,
            hud_pointer.after(ui::draw_hud),
            run_lab.after(hud_pointer),
            camera::camera_sys.after(run_lab),
            render_zz::aggiorna.after(camera::camera_sys),
            bagliore.after(run_lab),
            auto_sys.after(run_lab),
        ))
        .run();
}

/// Layers from a string of letters: t tracers, f wavefronts, l balls stretched along E, E E-streaks, B B-streaks, s Poynting,
/// o spin, x cross-section.
fn strati(s: &str) -> fotscena::Vista {
    fotscena::Vista {
        traccianti: s.contains('t'), fronti: s.contains('f'), allunga: s.contains('l'), e_vettori: s.contains('E'), b_vettori: s.contains('B'),
        poynting: s.contains('s'), spin: s.contains('o'), sezione: s.contains('x'),
    }
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

/// Unattended run for checks: `--foto file.png [--frames N] [--spara-al F] [--pannello] [--yaw a --pitch b]`.
/// Mouse and keyboard are ignored; a screenshot is saved at frame N and the program exits.
#[derive(Resource, Default)]
pub struct Auto {
    pub attivo: bool,
    pub foto: String,
    pub frames: u64,
    pub spara_al: Option<u64>,
    /// unattended: one measurement every 3 frames, this many in all
    pub misure: u64,
    pub yaw: Option<f32>,
    pub pitch: Option<f32>,
    pub fatto: bool,
}

fn auto_sys(mut commands: Commands, mut auto: ResMut<Auto>, mut lab: ResMut<Lab>, mut exit: EventWriter<AppExit>) {
    if !auto.attivo { return; }
    let f = lab.frame_count;
    if Some(f) == auto.spara_al { lab.spara = true; }
    if auto.misure > 0 && f >= 30 && f % 3 == 0 && (f - 30) / 3 < auto.misure { lab.misura = true; }
    if f == auto.frames && !auto.fatto {
        auto.fatto = true;
        let path = auto.foto.clone();
        commands.spawn(bevy::render::view::screenshot::Screenshot::primary_window())
            .observe(bevy::render::view::screenshot::save_to_disk(path));
        let info = match lab.scena {
            Scena::Fotone => lab.fot.as_ref().map(|z| {
                let o = z.lock().unwrap().out.clone();
                format!("tick {} v {:?} elicita {:?} w {:.2} theta {:?} {:.0} ms/tick", o.ticks, o.v, o.elicita, o.w, o.theta, o.ms_tick)
            }).unwrap_or_default(),
            Scena::Atomo => lab.zz.as_ref().map(|z| { let o = z.lock().unwrap().out.clone(); format!("tick {} p1s {:.3} dip {:+.3} {:.0} ms/tick", o.ticks, o.p1s, o.dipolo, o.ms_tick) }).unwrap_or_default(),
            Scena::Orbitali => lab.orb.as_ref().map(|z| { let o = z.lock().unwrap().out.clone(); format!("t {:.1} R {:.1} livelli {:?} {:.0} ms", o.t, o.r, o.livelli.iter().map(|l| (l.0.nome(), l.1, l.2)).collect::<Vec<_>>(), o.ms) }).unwrap_or_default(),
        };
        eprintln!("foto al frame {f}: scena {:?}, {} palline, {:.0} fps {info}", lab.scena, lab.palline_disegnate, lab.fps);
    }
    if f >= auto.frames + 20 { exit.send(AppExit::Success); }
}
