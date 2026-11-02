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
mod bagliori;
mod camera;
mod cascata;
mod elscena;
mod esportazioni;
#[path = "../zigzag/cplx.rs"]
#[allow(dead_code)]
mod cplx;
mod fotdisegno;
mod fotone3d;
mod istanze;
mod fotscena;
mod grafica;
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
mod tavola;
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
#[cfg(test)]
mod test_emissione;
#[cfg(test)]
mod test_sottoreticolo;
#[cfg(test)]
mod test_diffusa;
#[cfg(test)]
mod test_bagliori;

use camera::{Bersaglio, Cam, PointerOverHud};
use fotone3d::Polar;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Modo { Avanti, Pausa, Indietro }

/// The two scenes, both the zig-zag rule live: the free photon (light alone) and the atom (electron + light).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scena { Fotone, Atomo, Orbitali, Elettrone, Tavola }

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
    /// atom scene: the scattered view (docs/ATOMO.md §3.8)
    pub diffusa: bool,
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
    /// the free electron (or muon) drawn with glows (elscena.rs) and the particle the panel shows
    pub el: Option<std::sync::Arc<std::sync::Mutex<elscena::Condiviso>>>,
    pub el_part: elscena::Particella,
    pub ritmo_el: f64,
    /// version of the scene output last handed to the drawing (render_zz), for frame-by-frame sequences
    pub versione_disegnata: u64,
    /// the periodic-table mode (tavola.rs): the table, its scene worker, what was pressed
    pub tav: tavola::Pannello,
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
            scena: Scena::Fotone, fot: None, fot_par: par, ritmo_fot: 20.0, segui: false, vista_fot: fotscena::Vista::default(), lambda_m: 500e-9, zz: None, orb: None, orb_scelta: orbitali::Scelta::default(), ritmo_orb: 30.0, opz: misura::Opzioni::default(), accoppiato: true, diffusa: true, e_kev_acc: false, misura: false, azzera_clic: false, orb_r: 25.0, e_kev: 0.0, ritmo: 25.0,
            t_nucleo: 0.0, omega_rot: OMEGA_ROT_SCHERMO, omega_resp: OMEGA_ROT_SCHERMO * w_resp / w_rot, nodo_info2: info2,
            modo: Modo::Avanti, spara: false, panel_open: false, mostra_protone: true, bagliore: true, grandezza: 1.0, luce: 1.0,
            nodo_palline, nodo_info: info, frame_count: 0, palline_disegnate: 0, fps: 0.0,
            el: None, el_part: elscena::Particella::Elettrone, ritmo_el: 12.0, versione_disegnata: 0, tav: tavola::Pannello::default(),
        }
    }

    /// Metres per scene unit.
    pub fn metri_per_unita(&self) -> f64 {
        match self.scena {
            Scena::Atomo | Scena::Orbitali => zz::LAMBDA_C_M,
            Scena::Fotone => self.cella_m(),
            Scena::Elettrone => zz::LAMBDA_C_M / self.el_part.massa(),
            Scena::Tavola => tavola::A0_M,
        }
    }
    /// Seconds per tick of the scene's rule.
    pub fn secondi_per_tick(&self) -> f64 {
        match self.scena {
            Scena::Atomo => zz::A * zz::T_UNIT_S,
            Scena::Orbitali => zz::T_UNIT_S,
            Scena::Fotone => self.cella_m() / 299_792_458.0,
            Scena::Elettrone => elscena::THETA * zz::T_UNIT_S / self.el_part.massa(),
            Scena::Tavola => tavola::T_AU_S,
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
            Scena::Elettrone => { if self.el.is_none() { self.el = Some(elscena::avvia(self.el_part)); } }
            Scena::Tavola => self.tav.avvia(),
        }
        self.tav.comanda(s == Scena::Tavola && self.modo == Modo::Avanti, 30.0);
        if let Some(e) = &self.el {
            let mut g = e.lock().unwrap();
            g.comando = if s == Scena::Elettrone && self.modo == Modo::Avanti { elscena::Comando::Avanti } else { elscena::Comando::Pausa };
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
            Scena::Orbitali | Scena::Elettrone | Scena::Tavola => [0.0; 3],
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
            (Scena::Elettrone, _) => 16.0,
            (Scena::Tavola, _) => 2.6 * self.tav.out().map(|o| o.raggio).unwrap_or(3.0),
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
                    g.diffusa = self.diffusa;
                    if self.misura { g.misura = true; }
                    if self.azzera_clic { g.azzera = true; }
                    // the coupled atom's own 1s → 2p (the uncoupled radial value is not its resonance)
                    if g.out.accoppiato && !self.e_kev_acc { self.e_kev = g.out.omega_res * zz::MC2_KEV; self.e_kev_acc = true; }
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
            Scena::Elettrone => {
                if let Some(e) = &self.el {
                    let mut g = e.lock().unwrap();
                    g.comando = if self.modo == Modo::Avanti { elscena::Comando::Avanti } else { elscena::Comando::Pausa };
                    g.ritmo = self.ritmo_el;
                    if self.spara { self.spara = false; g.ricomincia = true; }
                    if g.out.particella != self.el_part && g.particella.is_none() { g.particella = Some(self.el_part); }
                }
            }
            Scena::Tavola => {
                self.tav.comanda(self.modo == Modo::Avanti, 30.0);
                // F (or "lancia"): hydrogen gets a photon on the chosen line (cascata.rs); elsewhere nothing is simulated
                if self.spara && self.tav.idrogeno() { let r = self.tav.riga_lancio(); self.tav.richiedi(|q| q.lancia = Some(r)); }
                self.spara = false;
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
    // the atom on one colour of the cube (docs/SOTTORETICOLO.md: 4× less work, but one colour alone is anisotropic):
    // `--reticolo colore|pieno`, or `--stati colore` (the name the bridge already accepts)
    match (arg("--reticolo").as_deref(), arg("--stati").as_deref()) {
        (Some("colore"), _) | (None, Some("colore")) => zz::RETICOLO.store(1, std::sync::atomic::Ordering::Relaxed),
        (Some("pieno"), _) => zz::RETICOLO.store(2, std::sync::atomic::Ordering::Relaxed),
        _ => {}
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
        ("elettrone", _) | ("muone", _) => Scena::Elettrone,
        ("tavola", _) => Scena::Tavola,
        _ => Scena::Fotone,
    };
    if scena_arg == "muone" || arg("--particella").as_deref() == Some("muone") { lab.el_part = elscena::Particella::Muone; }
    lab.cambia_scena(scena);
    if scena == Scena::Orbitali { if let Some(o) = &lab.orb { o.lock().unwrap().scegli = Some(lab.orb_scelta.clone()); } }
    // the periodic-table mode: `--elemento He` (or a molecule export name, e.g. H2), `--buca` for the Coulomb hole
    if scena == Scena::Tavola {
        if let Some(e) = arg("--elemento") { if let Some(s) = lab.tav.scegli_nome(&e) { lab.cambia_scena(s); } }
        lab.tav.buca = args.iter().any(|a| a == "--buca");
        // hydrogen: the line the photon is tuned to, by its upper level (`--riga 2p` Lyman-α, `--riga 3p` Lyman-β)
        if let Some(r) = arg("--riga") { lab.tav.riga = cascata::LIVELLI.iter().position(|x| x.0 == r).unwrap_or(2); }
    }
    let mut auto = Auto::default();
    if let Some(f) = arg("--foto") {
        auto.attivo = true;
        auto.foto = f;
        auto.frames = arg("--frames").and_then(|v| v.parse().ok()).unwrap_or(240);
        auto.spara_al = arg("--spara-al").and_then(|v| v.parse().ok());
        auto.misure = arg("--misure").and_then(|v| v.parse().ok()).unwrap_or(0);
        // atom scene: fire the pulse as soon as the coupled atom is prepared, and take the photo `--foto-dopo T`
        // (λ_C/c) after it (then `--frames` is only an upper limit)
        auto.spara_pronto = args.iter().any(|a| a == "--spara-pronto");
        auto.foto_dopo = arg("--foto-dopo").and_then(|v| v.parse().ok());
        if args.iter().any(|a| a == "--senza-diffusa") { lab.diffusa = false; }
        lab.panel_open = args.iter().any(|a| a == "--pannello");
        if let Some(r) = arg("--ritmo").and_then(|v| v.parse().ok()) { lab.ritmo_fot = r; }
        if let Some(r) = arg("--ritmo").and_then(|v| v.parse().ok()) { lab.ritmo_el = r; }
        if let Some(t) = arg("--tick").and_then(|v| v.parse::<u64>().ok()) {
            fotscena::TICK_FERMO.store(t, std::sync::atomic::Ordering::Relaxed);
            elscena::TICK_FERMO.store(t, std::sync::atomic::Ordering::Relaxed);
        }
        // `--sequenza K [--passo P]`: K pictures, one every P ticks of the scene's rule starting from `--tick`
        // (default 1), each frame saved when the drawing shows exactly that tick: a slow-motion clip, frame by frame
        if let Some(k) = arg("--sequenza").and_then(|v| v.parse::<u64>().ok()) {
            let t0 = arg("--tick").and_then(|v| v.parse::<u64>().ok()).unwrap_or(1).max(1);
            auto.seq = Some(Sequenza { k, passo: arg("--passo").and_then(|v| v.parse().ok()).unwrap_or(1).max(1), t0, fatte: 0, attesa: 0 });
            fotscena::TICK_FERMO.store(t0, std::sync::atomic::Ordering::Relaxed);
            elscena::TICK_FERMO.store(t0, std::sync::atomic::Ordering::Relaxed);
        }
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
        .insert_resource(render_zz::ModoDisegno { istanze: !args.iter().any(|a| a == "--mesh") })
        .insert_resource(ClearColor(Color::srgb(0.008, 0.01, 0.02)))
        .insert_resource(lab)
        .insert_resource(auto)
        .insert_resource(cam)
        .insert_resource(grafica::Grafica::da_argomenti(&args))
        .init_resource::<grafica::Budget>()
        .init_resource::<PointerOverHud>()
        .add_systems(Startup, (setup_scene, render_zz::setup, gpu_fotone::prendi_dispositivo))
        .add_systems(Update, (
            ui::draw_hud,
            grafica::pannello.after(ui::draw_hud).before(hud_pointer),
            grafica::budget_sys,
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
    /// atom scene: fire when the coupled atom is ready; photo this long (λ_C/c) after the pulse
    pub spara_pronto: bool,
    pub sparato: bool,
    pub foto_dopo: Option<f64>,
    pub frame_foto: u64,
    /// atom scene: the frame at which the pulse really left (it waits for the coupled atom to be prepared)
    pub frame_impulso: Option<u64>,
    /// frame-by-frame sequence (`--sequenza`)
    pub seq: Option<Sequenza>,
}

/// A slow-motion sequence: picture `fatte` shows tick t0 + fatte·passo.
#[derive(Clone, Copy, Debug)]
pub struct Sequenza { pub k: u64, pub passo: u64, pub t0: u64, pub fatte: u64, pub attesa: u64 }

/// (ticks, version) of the output of the scene shown.
fn tick_versione(lab: &Lab) -> Option<(u64, u64)> {
    match lab.scena {
        Scena::Fotone => lab.fot.as_ref().map(|z| { let o = &z.lock().unwrap().out; (o.ticks, o.versione) }),
        Scena::Atomo => lab.zz.as_ref().map(|z| { let o = &z.lock().unwrap().out; (o.ticks, o.versione) }),
        // the orbitals have no tick to stop at: no sequences there
        Scena::Orbitali => None,
        Scena::Elettrone => lab.el.as_ref().map(|z| { let o = &z.lock().unwrap().out; (o.ticks, o.versione) }),
        Scena::Tavola => None,
    }
}

fn auto_sys(mut commands: Commands, mut auto: ResMut<Auto>, mut lab: ResMut<Lab>, mut exit: EventWriter<AppExit>) {
    if !auto.attivo { return; }
    if let Some(mut q) = auto.seq {
        // wait until the drawing shows the wanted tick, let two frames pass, save, then ask for the next tick
        let voluto = q.t0 + q.fatte * q.passo;
        if q.fatte >= q.k {
            q.attesa += 1;
            if q.attesa > 20 { exit.send(AppExit::Success); }
        } else if let Some((t, v)) = tick_versione(&lab) {
            if t == voluto && lab.versione_disegnata == v {
                q.attesa += 1;
                if q.attesa == 3 {
                    let base = auto.foto.trim_end_matches(".png").to_string();
                    let path = format!("{base}_{:04}.png", q.fatte);
                    commands.spawn(bevy::render::view::screenshot::Screenshot::primary_window())
                        .observe(bevy::render::view::screenshot::save_to_disk(path));
                    if let Some(info) = info_scena(&lab) { eprintln!("sequenza {}/{}: {info}", q.fatte + 1, q.k); }
                }
                if q.attesa >= 5 {
                    q.fatte += 1;
                    q.attesa = 0;
                    let t = q.t0 + q.fatte * q.passo;
                    fotscena::TICK_FERMO.store(t, std::sync::atomic::Ordering::Relaxed);
                    elscena::TICK_FERMO.store(t, std::sync::atomic::Ordering::Relaxed);
                }
            } else { q.attesa = 0; }
        }
        auto.seq = Some(q);
        return;
    }
    let f = lab.frame_count;
    if Some(f) == auto.spara_al { lab.spara = true; }
    let o_zz = if lab.scena == Scena::Atomo { lab.zz.as_ref().map(|z| z.lock().unwrap().out.clone()) } else { None };
    if auto.spara_pronto && !auto.sparato {
        if let Some(o) = &o_zz { if o.pronto && o.accoppiato { lab.spara = true; auto.sparato = true; } }
    }
    if auto.frame_impulso.is_none() { if let Some(o) = &o_zz { if o.t_dal_impulso.is_some() { auto.frame_impulso = Some(f); } } }
    let pronta = match (auto.foto_dopo, &o_zz) {
        (Some(t), Some(o)) => o.t_dal_impulso.map(|x| x >= t).unwrap_or(false),
        // atom scene with --spara-al: the pulse waits for the prepared atom, the photo keeps the same distance in frames
        // from the pulse (frames − spara-al) as asked; never later than frame 20000
        (None, Some(_)) => match (auto.spara_al, auto.frame_impulso) {
            (Some(sa), Some(fi)) => f >= fi + auto.frames.saturating_sub(sa),
            _ => false,
        },
        _ => false,
    };
    // the table's hydrogen: `--spara-pronto` fires when its spectrum is computed, `--foto-dopo T` = T screen seconds after
    let o_tav = if lab.scena == Scena::Tavola { lab.tav.out() } else { None };
    if auto.spara_pronto && !auto.sparato { if let Some(o) = &o_tav { if o.spettro.is_some() { lab.spara = true; auto.sparato = true; } } }
    let pronta = pronta || matches!((auto.foto_dopo, o_tav.as_ref().and_then(|o| o.dal_lancio)), (Some(t), Some(x)) if x >= t);
    let atomo_attende = (o_zz.is_some() && auto.foto_dopo.is_none() && auto.spara_al.is_some()) || (o_tav.is_some() && auto.foto_dopo.is_some());
    if auto.misure > 0 && f >= 30 && f % 3 == 0 && (f - 30) / 3 < auto.misure { lab.misura = true; }
    if ((f == auto.frames && !atomo_attende) || pronta || f == 20_000) && !auto.fatto {
        auto.fatto = true;
        auto.frame_foto = f;
        let path = auto.foto.clone();
        commands.spawn(bevy::render::view::screenshot::Screenshot::primary_window())
            .observe(bevy::render::view::screenshot::save_to_disk(path));
        let info = info_scena(&lab).unwrap_or_default();
        eprintln!("foto al frame {f}: scena {:?}, {} palline, {:.0} fps {info}", lab.scena, lab.palline_disegnate, lab.fps);
    }
    if auto.fatto && f >= auto.frame_foto + 20 { exit.send(AppExit::Success); }
}

fn info_scena(lab: &Lab) -> Option<String> {
        let info = match lab.scena {
            Scena::Fotone => lab.fot.as_ref().map(|z| {
                let o = z.lock().unwrap().out.clone();
                format!("tick {} v {:?} elicita {:?} w {:.2} theta {:?} {:.0} ms/tick", o.ticks, o.v, o.elicita, o.w, o.theta, o.ms_tick)
            }).unwrap_or_default(),
            Scena::Atomo => lab.zz.as_ref().map(|z| { let o = z.lock().unwrap().out.clone(); format!("tick {} p1s {:.3} dip {:+.3} {:.0} ms/tick p_exc {:.3} t_dal_impulso {:?} luce_diffusa_x {:?} e_diffusa {:.3e}", o.ticks, o.p1s, o.dipolo, o.ms_tick, o.p_exc, o.t_dal_impulso, o.amp_diffusa, o.e_diffusa) }).unwrap_or_default(),
            Scena::Orbitali => lab.orb.as_ref().map(|z| { let o = z.lock().unwrap().out.clone(); format!("t {:.1} R {:.1} livelli {:?} {:.0} ms", o.t, o.r, o.livelli.iter().map(|l| (l.0.nome(), l.1, l.2)).collect::<Vec<_>>(), o.ms) }).unwrap_or_default(),
            Scena::Elettrone => lab.el.as_ref().map(|z| { let o = z.lock().unwrap().out.clone();
                format!("{} tick {} t {:.3} centro bagliori ({:+.4}, {:+.4}, {:+.4}) centro onda ({:+.4}, {:+.4}, {:+.4}) {:.0} ms/tick", o.particella.nome(), o.ticks, o.t,
                        o.centro_palline[0], o.centro_palline[1], o.centro_palline[2], o.centro_onda[0], o.centro_onda[1], o.centro_onda[2], o.ms_tick) }).unwrap_or_default(),
            Scena::Tavola => lab.tav.out().map(|o| format!("{} · {} · {} palline {:.0} ms", o.righe.0, o.righe.1, o.palline.len() + o.onde.len(), o.ms)).unwrap_or_default(),
        };
        Some(info)
}
