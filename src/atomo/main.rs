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

mod camera;
#[path = "../zigzag/cplx.rs"]
#[allow(dead_code)]
mod cplx;
mod fotone3d;
mod fotscena;
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

use camera::{Bersaglio, Cam, PointerOverHud};
use fotone3d::Polar;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Modo { Avanti, Pausa, Indietro }

/// The two scenes, both the zig-zag rule live: the free photon (light alone) and the atom (electron + light).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scena { Fotone, Atomo }

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
            scena: Scena::Fotone, fot: None, fot_par: par, ritmo_fot: 20.0, segui: false, vista_fot: fotscena::Vista::default(), lambda_m: 500e-9, zz: None, e_kev: 0.0, ritmo: 25.0,
            t_nucleo: 0.0, omega_rot: OMEGA_ROT_SCHERMO, omega_resp: OMEGA_ROT_SCHERMO * w_resp / w_rot, nodo_info2: info2,
            modo: Modo::Avanti, spara: false, panel_open: false, mostra_protone: true, bagliore: true, grandezza: 1.0, luce: 1.0,
            nodo_palline, nodo_info: info, frame_count: 0, palline_disegnate: 0, fps: 0.0,
        }
    }

    /// Metres per scene unit.
    pub fn metri_per_unita(&self) -> f64 {
        match self.scena {
            Scena::Atomo => zz::LAMBDA_C_M,
            Scena::Fotone => self.cella_m(),
        }
    }
    /// Seconds per tick of the scene's rule.
    pub fn secondi_per_tick(&self) -> f64 {
        match self.scena {
            Scena::Atomo => zz::A * zz::T_UNIT_S,
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
        }
    }
    pub fn distanza_tipica(&self, b: Bersaglio) -> f64 {
        let u = self.metri_per_unita();
        match (self.scena, b) {
            (_, Bersaglio::Nucleo) => 3.0e-15 / u,
            (Scena::Atomo, Bersaglio::Elettrone) => 3.0,
            (Scena::Atomo, _) => 11.0,
            (Scena::Fotone, _) => 1.45 * self.fot.as_ref().map(|f| f.lock().unwrap().out.par.n()).unwrap_or(self.fot_par.n()) as f64,
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
        match self.scena {
            Scena::Fotone => {
                if self.spara { self.spara = false; self.lancia_fotone(); }
                if let Some(f) = &self.fot {
                    let mut g = f.lock().unwrap();
                    g.comando = match self.modo { Modo::Avanti => fotscena::Comando::Avanti, Modo::Indietro => fotscena::Comando::Indietro, Modo::Pausa => fotscena::Comando::Pausa };
                    g.ritmo = self.ritmo_fot;
                    g.vista = self.vista_fot;
                }
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
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let arg = |k: &str| args.iter().position(|a| a == k).and_then(|i| args.get(i + 1)).cloned();
    let zoom = arg("--zoom").unwrap_or_else(|| "intero".into());
    let scena_arg = arg("--scena").unwrap_or_default();
    // the photon's parameters from the command line (for unattended pictures)
    let mut par = fotscena::Parametri::default();
    if let Some(v) = arg("--qual").and_then(|v| v.parse().ok()) { par.qual = v; }
    if let Some(v) = arg("--w0").and_then(|v| v.parse().ok()) { par.w0_lambda = v; }
    if let Some(v) = arg("--dir").and_then(|v| v.parse().ok()) { par.dir = v; }
    if let Some(p) = arg("--polar") {
        par.polar = match p.as_str() { "h" => Polar::Orizzontale, "v" => Polar::Verticale, "c-" => Polar::CircolareMeno, _ => Polar::CircolarePiu };
    }
    let mut lab = Lab::new(par);
    if let Some(v) = arg("--lambda-m").and_then(|v| v.parse().ok()) { lab.lambda_m = v; }
    if let Some(v) = arg("--strati") { lab.vista_fot = strati(&v); }
    let scena = match (scena_arg.as_str(), zoom.as_str()) {
        ("atomo", _) | ("", "nucleo") | ("", "protone") | ("", "elettrone") => Scena::Atomo,
        _ => Scena::Fotone,
    };
    lab.cambia_scena(scena);
    let mut auto = Auto::default();
    if let Some(f) = arg("--foto") {
        auto.attivo = true;
        auto.foto = f;
        auto.frames = arg("--frames").and_then(|v| v.parse().ok()).unwrap_or(240);
        auto.spara_al = arg("--spara-al").and_then(|v| v.parse().ok());
        lab.panel_open = args.iter().any(|a| a == "--pannello");
        if let Some(r) = arg("--ritmo").and_then(|v| v.parse().ok()) { lab.ritmo_fot = r; }
    }
    if let Some(y) = arg("--yaw").and_then(|v| v.parse::<f32>().ok()) { auto.yaw = Some(y); }
    if let Some(p) = arg("--pitch").and_then(|v| v.parse::<f32>().ok()) { auto.pitch = Some(p); }
    let b = match zoom.as_str() { "nucleo" | "protone" => Bersaglio::Nucleo, "elettrone" => Bersaglio::Elettrone, _ => Bersaglio::Atomo };
    let mut cam = Cam::new(lab.distanza_tipica(b), b);
    if scena == Scena::Fotone { cam.yaw = 0.0; cam.pitch = 0.2; }
    if let Some(y) = auto.yaw { cam.yaw = y; }
    if let Some(p) = auto.pitch { cam.pitch = p; }
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "Laboratorio Atomo — la regola zig-zag".into(),
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
            render_zz::aggiorna.after(camera::camera_sys),
            bagliore.after(run_lab),
            auto_sys.after(run_lab),
        ))
        .run();
}

/// Layers from a string of letters: e E-cloud, b B-cloud, E E-vectors, B B-vectors, s Poynting, h spin, x section, p polarizer.
fn strati(s: &str) -> fotscena::Vista {
    fotscena::Vista {
        e_nuvola: s.contains('e'), b_nuvola: s.contains('b'), e_vettori: s.contains('E'), b_vettori: s.contains('B'),
        poynting: s.contains('s'), spin: s.contains('h'), sezione: s.contains('x'), polarizzatore: s.contains('p'),
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
    pub yaw: Option<f32>,
    pub pitch: Option<f32>,
    pub fatto: bool,
}

fn auto_sys(mut commands: Commands, mut auto: ResMut<Auto>, mut lab: ResMut<Lab>, mut exit: EventWriter<AppExit>) {
    if !auto.attivo { return; }
    let f = lab.frame_count;
    if Some(f) == auto.spara_al { lab.spara = true; }
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
        };
        eprintln!("foto al frame {f}: scena {:?}, {} palline, {:.0} fps {info}", lab.scena, lab.palline_disegnate, lab.fps);
    }
    if f >= auto.frames + 20 { exit.send(AppExit::Success); }
}
