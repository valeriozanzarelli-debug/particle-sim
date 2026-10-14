//! PERIODIC TABLE mode of the Laboratorio Atomo (docs/ATOMO.md §10): the whole table, every element with its
//! status read from data (`assets/tavola_stato.json` + the export folders of `esportazioni.rs`). An element can be
//! pressed only when real simulation data stand behind it; everything else stays grey, "in arrivo", never faked.
//!
//! The scene of this mode shows the element (or molecule) that was pressed:
//!   • an EXPORT (He today, H₂ and the others as the studies hand them over): glows drawn from the exported
//!     photographs of all the electrons together (|Ψ|² samples), fresh ones every few frames, each one fading in
//!     and out where it was drawn — the motion is the sampler, NOT a trajectory (declared). "Buca di Coulomb":
//!     one electron is turned onto a fixed direction (allowed by the symmetry of the state) and the others are
//!     drawn around it: they avoid it.
//!   • hydrogen: the absorption and cascade scene of `cascata.rs`.
//!   • the ion Nd⁵⁹⁺: the live zig-zag atom scene (Scena::Atomo).

use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use bevy_egui::egui;
use egui::{Color32, RichText};
use serde_json::Value;

use crate::camera::{Bersaglio, Cam};
use crate::esportazioni::{self, esami_da, Esame, Esportazione, Simmetria};
use crate::render::{schermo, Ball};
use crate::rng::Rng;
use crate::{Lab, Scena};

/// Bohr radius (m): the unit of length of this scene.
pub const A0_M: f64 = 5.291_772_109_03e-11;
/// Atomic unit of time (s).
pub const T_AU_S: f64 = 2.418_884_326_585_7e-17;

pub const SIMBOLI: [&str; 118] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

pub fn simbolo(z: u32) -> &'static str { SIMBOLI.get(z.saturating_sub(1) as usize).copied().unwrap_or("?") }
pub fn z_di(s: &str) -> Option<u32> { SIMBOLI.iter().position(|x| *x == s).map(|i| i as u32 + 1) }

/// Row and column (0-based) of element Z in the usual layout: 7 periods, then the lanthanides (row 8) and the
/// actinides (row 9) apart. Only a layout, no physics.
pub fn posizione(z: u32) -> (usize, usize) {
    let z = z as usize;
    match z {
        1 => (0, 0),
        2 => (0, 17),
        3..=10 | 11..=18 => {
            let (p, i) = if z <= 10 { (1, z - 3) } else { (2, z - 11) };
            (p, if i < 2 { i } else { i + 10 })
        }
        19..=36 => (3, z - 19),
        37..=54 => (4, z - 37),
        55..=56 => (5, z - 55),
        57..=71 => (8, z - 57 + 2),
        72..=86 => (5, z - 72 + 3),
        87..=88 => (6, z - 87),
        89..=103 => (9, z - 89 + 2),
        _ => (6, z.saturating_sub(104) + 3),
    }
}

/// The status of one table entry from `assets/tavola_stato.json`.
#[derive(Clone, Debug, Default)]
pub struct Voce {
    pub disponibile: bool,
    pub livello: Option<String>,
    pub scena: Option<String>,
    pub validazione: String,
    pub fonte: String,
    pub esami: Vec<Esame>,
    pub nota: String,
}

fn voce_da(v: &Value) -> Voce {
    let s = |k: &str| v[k].as_str().unwrap_or("").to_string();
    Voce {
        disponibile: v["disponibile"].as_bool().unwrap_or(false),
        livello: v["livello"].as_str().map(String::from),
        scena: v["scena"].as_str().map(String::from),
        validazione: s("validazione"), fonte: s("fonte"), esami: esami_da(&v["esami"]), nota: s("nota"),
    }
}

#[derive(Clone, Debug)]
pub struct Ione { pub nome: String, pub z: u32, pub voce: Voce }

/// How an element can be shown.
#[derive(Clone, Debug, PartialEq)]
pub enum Accesso { Scena(String), Esportazione(usize), No }

/// The whole table: statuses, ions, molecule notes, exports (loaded once).
#[derive(Clone, Debug, Default)]
pub struct Tavola {
    pub elementi: HashMap<String, Voce>,
    pub ioni: Vec<Ione>,
    pub molecole_note: Vec<(String, String)>,
    pub esportazioni: Vec<Arc<Esportazione>>,
    pub errori: Vec<String>,
}

impl Tavola {
    pub fn da_testo(stato: &str, esp: Vec<Esportazione>, mut errori: Vec<String>) -> Self {
        let v: Value = match serde_json::from_str(stato) { Ok(v) => v, Err(e) => { errori.push(format!("tavola_stato.json: {e}")); Value::Null } };
        let elementi = v["elementi"].as_object().map(|o| o.iter().map(|(k, x)| (k.clone(), voce_da(x))).collect()).unwrap_or_default();
        let ioni = v["ioni"].as_array().map(|a| a.iter().map(|x| Ione {
            nome: x["nome"].as_str().unwrap_or("?").into(), z: x["Z"].as_u64().unwrap_or(0) as u32, voce: { let mut w = voce_da(x); w.disponibile = true; w },
        }).collect()).unwrap_or_default();
        let molecole_note = v["molecole"].as_object().map(|o| o.iter().map(|(k, x)| (k.clone(), x["nota"].as_str().unwrap_or("in arrivo").to_string())).collect()).unwrap_or_default();
        Tavola { elementi, ioni, molecole_note, esportazioni: esp.into_iter().map(Arc::new).collect(), errori }
    }

    /// Read the status file and every export folder.
    pub fn carica() -> Self {
        let dir = esportazioni::cartella_assets();
        let (esp, mut err) = esportazioni::tutte();
        let stato = std::fs::read_to_string(dir.join("tavola_stato.json")).unwrap_or_else(|e| {
            err.push(format!("tavola_stato.json non letto ({e}): uso la copia compilata"));
            include_str!("../../assets/tavola_stato.json").to_string()
        });
        Tavola::da_testo(&stato, esp, err)
    }

    /// The neutral-atom export of element Z, if any.
    pub fn esportazione_atomo(&self, z: u32) -> Option<usize> {
        self.esportazioni.iter().position(|e| !e.molecola() && e.z == z && e.carica == 0)
    }

    /// How element Z can be shown: its own scene (status file) or an export; else not at all.
    pub fn accesso(&self, z: u32) -> Accesso {
        if let Some(v) = self.elementi.get(simbolo(z)) {
            if v.disponibile { if let Some(s) = &v.scena { return Accesso::Scena(s.clone()); } }
        }
        match self.esportazione_atomo(z) { Some(i) => Accesso::Esportazione(i), None => Accesso::No }
    }

    /// Tier, exams, source and validation of element Z, from its export or the status file.
    pub fn scheda(&self, z: u32) -> Voce {
        match self.accesso(z) {
            Accesso::Esportazione(i) => {
                let e = &self.esportazioni[i];
                Voce { disponibile: true, livello: Some(e.livello.clone()), scena: None, validazione: e.metodo.clone(), fonte: e.fonte.clone(),
                       esami: e.esami.clone(), nota: e.note.join(" · ") }
            }
            _ => self.elementi.get(simbolo(z)).cloned().unwrap_or_default(),
        }
    }

    pub fn molecole(&self) -> impl Iterator<Item = (usize, &Arc<Esportazione>)> { self.esportazioni.iter().enumerate().filter(|(_, e)| e.molecola()) }

    /// "Costruisci molecole": the molecule export made of exactly these nuclei, if one exists. Never a fake bond.
    pub fn costruisci(&self, banco: &[u32]) -> Option<usize> {
        let mut b = banco.to_vec();
        b.sort_unstable();
        if b.len() < 2 { return None; }
        self.molecole().find(|(_, e)| e.composizione() == b).map(|(i, _)| i)
    }
}

/// Formula from a list of Z (Hill-like order: C, H first, then alphabetical), e.g. [1, 8, 1] → "H2O".
pub fn formula(banco: &[u32]) -> String {
    let mut c: Vec<(String, usize)> = vec![];
    for z in banco { let s = simbolo(*z).to_string(); if let Some(x) = c.iter_mut().find(|x| x.0 == s) { x.1 += 1 } else { c.push((s, 1)) } }
    let ha_c = c.iter().any(|x| x.0 == "C");
    c.sort_by_key(|x| (if ha_c && x.0 == "C" { 0 } else if ha_c && x.0 == "H" { 1 } else { 2 }, x.0.clone()));
    c.iter().map(|(s, n)| if *n > 1 { format!("{s}{n}") } else { s.clone() }).collect()
}

// ───────────────────────────── the scene worker ─────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Comando { Avanti, Pausa }

#[derive(Clone)]
pub enum Oggetto { Esportazione(Arc<Esportazione>), Idrogeno }

#[derive(Clone)]
pub struct Uscita {
    pub palline: Arc<Vec<Ball>>,
    pub onde: Arc<Vec<Ball>>,
    pub versione: u64,
    pub pronto: bool,
    /// what is shown and how large it is (bohr)
    pub nome: String,
    pub raggio: f64,
    pub righe: (String, String),
    pub ms: f64,
}

pub struct Condiviso {
    pub comando: Comando,
    pub ritmo: f64,
    pub metti: Option<Oggetto>,
    pub buca: bool,
    pub fine: bool,
    pub out: Uscita,
}

/// A random rotation (uniform on SO(3)) from a random unit quaternion.
fn rotazione(rng: &mut Rng) -> [[f64; 3]; 3] {
    let (mut a, mut b, mut c, mut d) = (rng.gauss(), rng.gauss(), rng.gauss(), rng.gauss());
    let n = (a * a + b * b + c * c + d * d).sqrt().max(1e-12);
    a /= n; b /= n; c /= n; d /= n;
    [[a * a + b * b - c * c - d * d, 2.0 * (b * c - a * d), 2.0 * (b * d + a * c)],
     [2.0 * (b * c + a * d), a * a - b * b + c * c - d * d, 2.0 * (c * d - a * b)],
     [2.0 * (b * d - a * c), 2.0 * (c * d + a * b), a * a - b * b - c * c + d * d]]
}
fn rot_z(phi: f64) -> [[f64; 3]; 3] { let (s, c) = phi.sin_cos(); [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]] }
fn applica(m: &[[f64; 3]; 3], p: [f64; 3]) -> [f64; 3] {
    [m[0][0] * p[0] + m[0][1] * p[1] + m[0][2] * p[2], m[1][0] * p[0] + m[1][1] * p[1] + m[1][2] * p[2], m[2][0] * p[0] + m[2][1] * p[1] + m[2][2] * p[2]]
}
/// The rotation that takes the unit vector u onto +z (Rodrigues).
fn verso_z(u: [f64; 3]) -> [[f64; 3]; 3] {
    let (x, y, z) = (u[0], u[1], u[2]);
    if z < -0.999_999 { return [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]; }
    // axis u × z = (y, −x, 0), cos = z
    let k = 1.0 / (1.0 + z);
    [[1.0 - x * x * k, -x * y * k, -x], [-x * y * k, 1.0 - y * y * k, -y], [x, y, z]]
}

/// One photograph slot: a configuration, how it is turned, its age and its life (steps).
#[derive(Clone, Copy)]
struct Slot { i: usize, rot: [[f64; 3]; 3], eta: u32, vita: u32, rif: Option<usize> }

/// Electron glows per frame (all configurations together).
pub const N_BAGLIORI: usize = 5000;
const VITA: u32 = 8;
/// Coulomb hole: the reference electron is kept at this distance from the nucleus (bohr), on +z.
pub const BUCA_R: (f64, f64) = (0.35, 0.75);

pub struct Nuvola { e: Arc<Esportazione>, slot: Vec<Slot>, rng: Rng }

impl Nuvola {
    pub fn new(e: Arc<Esportazione>) -> Self {
        let k = (N_BAGLIORI / e.n_el.max(1)).clamp(100, 4000).min(e.n.max(1));
        let rng = Rng::new(31);
        let slot = (0..k).map(|j| Slot { i: 0, rot: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], eta: (j as u32) % VITA, vita: VITA, rif: None }).collect();
        let mut n = Nuvola { e, slot, rng };
        for j in 0..n.slot.len() { let eta = n.slot[j].eta; n.slot[j] = n.pesca(false); n.slot[j].eta = eta; }
        n
    }

    /// A fresh photograph: a random configuration turned by a rotation the symmetry allows. With the hole on, one
    /// electron is turned onto +z (spherical) or onto the x–z half plane x > 0 (axial), and kept only if it sits
    /// where the reference is held.
    fn pesca(&mut self, buca: bool) -> Slot {
        let e = &self.e;
        for _ in 0..200 {
            let i = (self.rng.next_u64() % e.n as u64) as usize;
            if !buca || e.n_el < 2 || e.simmetria == Simmetria::Nessuna {
                let rot = match e.simmetria { Simmetria::Sferica => rotazione(&mut self.rng), Simmetria::AssialeZ => rot_z(std::f64::consts::TAU * self.rng.uni()), Simmetria::Nessuna => rot_z(0.0) };
                return Slot { i, rot, eta: 0, vita: VITA, rif: None };
            }
            let r = (self.rng.next_u64() % e.n_el as u64) as usize;
            let p = e.el(i, r);
            match e.simmetria {
                Simmetria::Sferica => {
                    let d = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                    if d < BUCA_R.0 || d > BUCA_R.1 { continue; }
                    let mut rot = verso_z([p[0] / d, p[1] / d, p[2] / d]);
                    // any turn about z afterwards is still allowed
                    rot = mul(&rot_z(std::f64::consts::TAU * self.rng.uni()), &rot);
                    return Slot { i, rot, eta: 0, vita: VITA, rif: Some(r) };
                }
                _ => {
                    // axial: reference on the right half (z > z_max/4 of the nuclei, at least 0.35) and near the axis
                    let zc = e.nuclei.iter().map(|n| n.1[2].abs()).fold(0.0, f64::max).max(0.7) * 0.5;
                    let rho = (p[0] * p[0] + p[1] * p[1]).sqrt();
                    if p[2] < zc || rho > 0.8 { continue; }
                    return Slot { i, rot: rot_z(-p[1].atan2(p[0])), eta: 0, vita: VITA, rif: Some(r) };
                }
            }
        }
        Slot { i: 0, rot: rot_z(0.0), eta: 0, vita: VITA, rif: None }
    }

    pub fn passo(&mut self, buca: bool) {
        for j in 0..self.slot.len() {
            let s = &mut self.slot[j];
            s.eta += 1;
            if s.eta >= s.vita || (s.rif.is_some() != buca) { self.slot[j] = self.pesca(buca); }
        }
    }

    /// The glows: electrons (cyan-white; with the hole, the reference electron gold and the others blue), nuclei.
    pub fn disegna(&self, buca: bool) -> Vec<Ball> {
        let e = &self.e;
        let mut out = Vec::with_capacity(self.slot.len() * e.n_el + e.nuclei.len());
        for s in &self.slot {
            let w = (std::f32::consts::PI * (s.eta as f32 + 0.5) / s.vita as f32).sin();
            for k in 0..e.n_el {
                let p = applica(&s.rot, e.el(s.i, k));
                let (col, sz) = match (buca, s.rif) {
                    (true, Some(r)) if r == k => ([1.0, 0.82, 0.3, 0.55 * w], 0.6),
                    (true, Some(_)) => ([0.35, 0.6, 1.0, 0.45 * w], 1.0),
                    _ => ([0.55, 0.85, 1.0, 0.42 * w], 1.0),
                };
                out.push((schermo(p).to_array(), sz, col));
            }
        }
        for (z, p) in &e.nuclei {
            // the nucleus: a warm glow, larger for a larger charge (drawing choice, its real size is 10⁻⁵ of this)
            out.push((schermo(*p).to_array(), 1.5 + 0.3 * (*z as f32).sqrt(), [1.0, 0.6, 0.25, 0.9]));
        }
        out
    }

    /// Radius (bohr) that holds 99 % of the electrons of the first 2000 photographs.
    pub fn raggio(&self) -> f64 {
        let e = &self.e;
        let mut r: Vec<f64> = (0..e.n.min(2000)).flat_map(|i| (0..e.n_el).map(move |k| { let p = e.el(i, k); (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt() })).collect();
        r.sort_by(|a, b| a.total_cmp(b));
        r.get((r.len() as f64 * 0.99) as usize).copied().unwrap_or(3.0)
    }
}

fn mul(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut o = [[0.0; 3]; 3];
    for i in 0..3 { for j in 0..3 { for k in 0..3 { o[i][j] += a[i][k] * b[k][j]; } } }
    o
}

fn righe_esportazione(e: &Esportazione, buca: bool) -> (String, String) {
    let (p, f) = e.n_pass();
    (format!("{} ({}) · livello {} · {} elettroni · {} fotografie di |Ψ|² · esami {} PASS / {} FAIL", e.nome, e.simbolo, e.livello.to_uppercase(), e.n_el, e.n, p, f),
     format!("{}{}", e.energia.map(|x| format!("E = {:.6} Eh{} · ", x, e.riferimento.map(|r| format!(" (rif. {:.6})", r)).unwrap_or_default())).unwrap_or_default(),
             if buca { "BUCA DI COULOMB: un elettrone (oro) tenuto fermo, gli altri (blu) lo evitano" } else { "bagliori = fotografie fresche dello stato (il moto è il campionatore, non una traiettoria)" }))
}

enum Mondo { Vuoto, Nuvola(Nuvola) }

pub fn avvia() -> Arc<Mutex<Condiviso>> {
    let c = Arc::new(Mutex::new(Condiviso {
        comando: Comando::Avanti, ritmo: 30.0, metti: None, buca: false, fine: false,
        out: Uscita { palline: Arc::new(vec![]), onde: Arc::new(vec![]), versione: 0, pronto: false, nome: String::new(), raggio: 3.0,
                      righe: ("Tavola periodica: premi un elemento acceso".into(), String::new()), ms: 0.0 },
    }));
    let cc = c.clone();
    std::thread::spawn(move || {
        let mut mondo = Mondo::Vuoto;
        let mut versione = 0u64;
        let mut ms = 0.0;
        let mut ultimo = Instant::now();
        loop {
            let t0 = Instant::now();
            let (cmd, ritmo, metti, buca, fine) = {
                let mut g = cc.lock().unwrap();
                (g.comando, g.ritmo, g.metti.take(), g.buca, g.fine)
            };
            if fine { break; }
            let _dt_schermo = t0.duration_since(ultimo).as_secs_f64().min(0.1);
            ultimo = t0;
            let mut nome = None;
            if let Some(o) = metti {
                mondo = match o {
                    Oggetto::Esportazione(e) => { nome = Some(e.nome.clone()); Mondo::Nuvola(Nuvola::new(e)) }
                    Oggetto::Idrogeno => { nome = Some("idrogeno".into()); Mondo::Vuoto }
                };
            }
            let avanti = cmd == Comando::Avanti;
            let disegno = match &mut mondo {
                Mondo::Vuoto => None,
                Mondo::Nuvola(n) => {
                    if avanti { n.passo(buca); }
                    Some((n.disegna(buca), vec![], n.raggio(), righe_esportazione(&n.e, buca)))
                }
            };
            if let Some((p, o, r, righe)) = disegno {
                versione += 1;
                let el = t0.elapsed().as_secs_f64() * 1e3;
                ms += 0.15 * (el.min(500.0) - ms);
                let mut g = cc.lock().unwrap();
                let nome = nome.unwrap_or_else(|| g.out.nome.clone());
                g.out = Uscita { palline: Arc::new(p), onde: Arc::new(o), versione, pronto: true, nome, raggio: r, righe, ms };
            }
            let want = Duration::from_secs_f64(1.0 / ritmo.max(0.5));
            let el = t0.elapsed();
            if el < want { std::thread::sleep(want - el); }
            if !avanti { std::thread::sleep(Duration::from_millis(15)); }
        }
    });
    c
}

// ───────────────────────────── the panel ─────────────────────────────

/// Everything the periodic-table mode keeps in the Lab.
#[derive(Default)]
pub struct Pannello {
    pub tavola: Option<Arc<Tavola>>,
    pub mondo: Option<Arc<Mutex<Condiviso>>>,
    /// the element (Z) or molecule export shown, and what was last pressed
    pub scelto: Option<u32>,
    pub molecola: Option<usize>,
    pub buca: bool,
    /// "costruisci molecole": the nuclei put on the bench
    pub banco: Vec<u32>,
    pub messaggio: String,
    /// radius (bohr) the camera was last set for
    pub raggio_cam: f64,
}

impl Pannello {
    pub fn avvia(&mut self) {
        if self.tavola.is_none() { self.tavola = Some(Arc::new(Tavola::carica())); }
        if self.mondo.is_none() { self.mondo = Some(avvia()); }
    }
    pub fn out(&self) -> Option<Uscita> { self.mondo.as_ref().map(|m| m.lock().unwrap().out.clone()) }
    pub fn comanda(&self, avanti: bool, ritmo: f64) {
        if let Some(m) = &self.mondo { let mut g = m.lock().unwrap(); g.comando = if avanti { Comando::Avanti } else { Comando::Pausa }; g.ritmo = ritmo; g.buca = self.buca; }
    }
    fn metti(&self, o: Oggetto) { if let Some(m) = &self.mondo { m.lock().unwrap().metti = Some(o); } }
    /// What the scene shows is hydrogen (its own photon and cascade).
    pub fn idrogeno(&self) -> bool { self.scelto == Some(1) && self.molecola.is_none() }

    /// Press element Z: returns the scene to switch to, if it lives in another scene (the ion → Scena::Atomo).
    pub fn premi(&mut self, z: u32) -> Option<Scena> {
        let Some(t) = self.tavola.clone() else { return None };
        match t.accesso(z) {
            Accesso::Scena(s) if s == "idrogeno" => { self.scelto = Some(z); self.molecola = None; self.metti(Oggetto::Idrogeno); self.messaggio.clear(); None }
            Accesso::Scena(s) if s == "atomo" => Some(Scena::Atomo),
            Accesso::Scena(s) => { self.messaggio = format!("scena '{s}' sconosciuta a questa versione del programma"); None }
            Accesso::Esportazione(i) => { self.scelto = Some(z); self.molecola = None; self.metti(Oggetto::Esportazione(t.esportazioni[i].clone())); self.messaggio.clear(); None }
            Accesso::No => { self.messaggio = format!("{}: in arrivo (nessun dato calcolato ancora)", simbolo(z)); None }
        }
    }
    pub fn mostra_molecola(&mut self, i: usize) {
        let Some(t) = self.tavola.clone() else { return };
        if let Some(e) = t.esportazioni.get(i) { self.molecola = Some(i); self.scelto = None; self.metti(Oggetto::Esportazione(e.clone())); self.messaggio.clear(); }
    }
    /// From the command line: an element symbol ("He", "H") or a molecule export name ("H2").
    pub fn scegli_nome(&mut self, s: &str) -> Option<Scena> {
        if let Some(z) = z_di(s) { return self.premi(z); }
        let t = self.tavola.clone()?;
        if let Some((i, _)) = t.molecole().find(|(_, e)| e.simbolo == s || e.nome == s) { self.mostra_molecola(i); }
        None
    }
}

const VERDE: Color32 = Color32::from_rgb(60, 140, 90);
const AMBRA: Color32 = Color32::from_rgb(170, 125, 40);
const GRIGIO: Color32 = Color32::from_rgb(38, 40, 48);
const ORO: Color32 = Color32::from_rgb(255, 205, 80);
const TENUE: Color32 = Color32::from_rgb(150, 155, 170);

fn colore_livello(l: Option<&str>) -> Color32 { match l { Some("esatto") => VERDE, Some(_) => AMBRA, None => VERDE } }
fn badge(l: Option<&str>) -> &'static str { match l { Some("esatto") => "E", Some("approssimato") => "A", _ => "" } }

fn mostra_esami(ui: &mut egui::Ui, esami: &[Esame]) {
    for e in esami {
        let (s, c) = if e.pass { ("PASS", Color32::from_rgb(140, 230, 140)) } else { ("FAIL", Color32::from_rgb(255, 120, 100)) };
        ui.horizontal_wrapped(|ui| { ui.label(RichText::new(s).monospace().strong().color(c)); ui.label(RichText::new(format!("{} — {}", e.id, e.cosa)).small()); });
    }
}

fn scheda(ui: &mut egui::Ui, titolo: &str, v: &Voce) {
    ui.horizontal(|ui| {
        ui.label(RichText::new(titolo).strong().color(ORO));
        if let Some(l) = &v.livello {
            ui.label(RichText::new(format!(" {} ", l.to_uppercase())).strong().color(Color32::WHITE).background_color(colore_livello(Some(l))));
        }
    });
    if !v.validazione.is_empty() { ui.label(RichText::new(&v.validazione).small()); }
    if !v.fonte.is_empty() { ui.label(RichText::new(format!("dati: {}", v.fonte)).small().color(TENUE)); }
    mostra_esami(ui, &v.esami);
    if !v.nota.is_empty() { ui.label(RichText::new(&v.nota).small().color(TENUE)); }
}

/// The periodic-table window (scene "tavola").
pub fn pannello(ctx: &egui::Context, lab: &mut Lab, cam: &mut Cam) {
    lab.tav.avvia();
    let Some(t) = lab.tav.tavola.clone() else { return };
    // the camera follows the size of what is shown
    if let Some(o) = lab.tav.out() {
        if o.pronto && (o.raggio - lab.tav.raggio_cam).abs() > 1e-6 && cam.bersaglio == Bersaglio::Atomo {
            lab.tav.raggio_cam = o.raggio;
            cam.vai(Bersaglio::Atomo, 2.6 * o.raggio);
        }
    }
    if !lab.panel_open { return; }
    let mut open = true;
    let mut vai_a: Option<Scena> = None;
    egui::Window::new("Tavola periodica — solo ciò che è calcolato si accende").open(&mut open).default_pos([12.0, 80.0]).default_width(760.0).vscroll(true).show(ctx, |ui| {
        ui.label(RichText::new("Premi un elemento acceso: appare nella scena. Verde = livello ESATTO, ambra = APPROSSIMATO (HF + correlazione \
            imparata). Grigio = in arrivo: nessun dato, quindi niente disegno finto.").small().color(TENUE));
        // the grid
        let w = 38.0;
        egui::Grid::new("tavola").spacing([3.0, 3.0]).show(ui, |ui| {
            for riga in 0..10 {
                for col in 0..18 {
                    let z = (1..=118u32).find(|z| posizione(*z) == (riga, col));
                    match z {
                        None => { ui.allocate_exact_size(egui::vec2(w, if riga == 7 { 6.0 } else { 34.0 }), egui::Sense::hover()); }
                        Some(z) => {
                            let acc = t.accesso(z);
                            let v = t.scheda(z);
                            let on = acc != Accesso::No;
                            let fill = if on { colore_livello(v.livello.as_deref()) } else { GRIGIO };
                            let txt = RichText::new(format!("{}\n{}{}", simbolo(z), z, if on { format!(" {}", badge(v.livello.as_deref())) } else { String::new() }))
                                .size(11.0).color(if on { Color32::WHITE } else { Color32::from_gray(95) });
                            let b = egui::Button::new(txt).fill(fill).min_size(egui::vec2(w, 34.0)).selected(lab.tav.scelto == Some(z));
                            let r = ui.add(b).on_hover_ui(|ui| {
                                if on { scheda(ui, &format!("{} (Z = {z})", simbolo(z)), &v); }
                                else { ui.label(format!("{} (Z = {z}): in arrivo{}", simbolo(z), if v.nota.is_empty() { String::new() } else { format!(" — {}", v.nota) })); }
                            });
                            if r.clicked() {
                                if let Some(s) = lab.tav.premi(z) { vai_a = Some(s); }
                            }
                            if r.secondary_clicked() { lab.tav.banco.push(z); }
                        }
                    }
                }
                ui.end_row();
            }
        });
        if !lab.tav.messaggio.is_empty() { ui.label(RichText::new(&lab.tav.messaggio).color(Color32::from_rgb(255, 170, 120))); }
        ui.separator();
        // ions
        ui.horizontal_wrapped(|ui| {
            ui.label(RichText::new("IONI A UN ELETTRONE").strong().color(ORO));
            for io in &t.ioni {
                if ui.add(egui::Button::new(RichText::new(format!("{} ({})", io.nome, io.voce.livello.as_deref().unwrap_or("?"))).color(Color32::WHITE)).fill(VERDE))
                    .on_hover_ui(|ui| scheda(ui, &io.nome, &io.voce)).clicked() {
                    if io.voce.scena.as_deref() == Some("atomo") { vai_a = Some(Scena::Atomo); }
                }
            }
        });
        ui.separator();
        // what is shown
        if let Some(z) = lab.tav.scelto { scheda(ui, &format!("IN SCENA: {} (Z = {z})", simbolo(z)), &t.scheda(z)); }
        if let Some(i) = lab.tav.molecola { let e = &t.esportazioni[i]; scheda(ui, &format!("IN SCENA: molecola {}", e.nome), &Voce {
            disponibile: true, livello: Some(e.livello.clone()), scena: None, validazione: e.metodo.clone(), fonte: e.fonte.clone(), esami: e.esami.clone(), nota: e.note.join(" · ") }); }
        // projectiles
        ui.separator();
        ui.label(RichText::new("LANCIA").strong().color(ORO));
        let idro = lab.tav.idrogeno();
        ui.horizontal_wrapped(|ui| {
            let r = ui.add_enabled(idro, egui::Button::new(RichText::new(" fotone (F) ").strong()).fill(Color32::from_rgb(40, 60, 110)));
            if r.clicked() { lab.spara = true; }
            if !idro { r.on_disabled_hover_text("per questo oggetto la luce non è ancora simulata (idrogeno: sì; Nd⁵⁹⁺: nella scena 'atomo zig-zag')"); }
            for n in ["elettrone", "protone", "neutrone"] {
                ui.add_enabled(false, egui::Button::new(n)).on_disabled_hover_text("non ancora simulato: nessuna regola calcolata per l'urto");
            }
            ui.label(RichText::new("elettrone/protone: non ancora simulato").small().color(TENUE));
        });
        // the Coulomb hole
        if let Some(o) = lab.tav.molecola.or_else(|| lab.tav.scelto.and_then(|z| t.esportazione_atomo(z))) {
            let e = &t.esportazioni[o];
            ui.add_enabled(e.n_el >= 2 && e.simmetria != Simmetria::Nessuna, egui::Checkbox::new(&mut lab.tav.buca,
                "buca di Coulomb: tieni fermo un elettrone e guarda dove stanno gli altri"));
        }
        // molecules
        ui.separator();
        ui.label(RichText::new("MOLECOLE").strong().color(ORO));
        ui.horizontal_wrapped(|ui| {
            for (i, e) in t.molecole() {
                let (p, f) = e.n_pass();
                if ui.add(egui::Button::new(RichText::new(format!("{} [{}]", e.nome, badge(Some(&e.livello)))).color(Color32::WHITE)).fill(colore_livello(Some(&e.livello))).selected(lab.tav.molecola == Some(i)))
                    .on_hover_text(format!("{} · {} PASS / {} FAIL · {}", e.metodo, p, f, e.fonte)).clicked() { lab.tav.mostra_molecola(i); }
            }
            for (n, nota) in &t.molecole_note {
                if t.molecole().any(|(_, e)| &e.simbolo == n || &e.nome == n) { continue; }
                ui.add_enabled(false, egui::Button::new(n.as_str())).on_disabled_hover_text(nota.as_str());
            }
        });
        ui.label(RichText::new("COSTRUISCI MOLECOLE (clic destro su un elemento = mettilo sul banco)").small().strong().color(TENUE));
        ui.horizontal_wrapped(|ui| {
            if let Some(z) = lab.tav.scelto { if ui.button(format!("+ {}", simbolo(z))).clicked() { lab.tav.banco.push(z); } }
            ui.label(RichText::new(if lab.tav.banco.is_empty() { "banco vuoto".to_string() } else { format!("banco: {}", formula(&lab.tav.banco)) }).monospace());
            match t.costruisci(&lab.tav.banco) {
                Some(i) => { if ui.add(egui::Button::new(RichText::new(format!(" unisci → {} ", t.esportazioni[i].nome)).strong()).fill(VERDE)).clicked() { lab.tav.mostra_molecola(i); } }
                None if lab.tav.banco.len() >= 2 => { ui.add_enabled(false, egui::Button::new("unisci")); ui.label(RichText::new("non ancora calcolata: nessun legame disegnato").color(Color32::from_rgb(255, 170, 120))); }
                None => {}
            }
            if ui.button("svuota").clicked() { lab.tav.banco.clear(); }
        });
        if !t.errori.is_empty() {
            ui.separator();
            for e in &t.errori { ui.label(RichText::new(format!("dato non letto: {e}")).small().color(Color32::from_rgb(255, 120, 100))); }
        }
        ui.horizontal(|ui| {
            ui.add(egui::Slider::new(&mut lab.grandezza, 0.3..=3.0).text("palline"));
            ui.add(egui::Slider::new(&mut lab.luce, 0.2..=4.0).text("luce"));
        });
    });
    if !open { lab.panel_open = false; }
    if let Some(s) = vai_a {
        lab.cambia_scena(s);
        let r = lab.distanza_tipica(Bersaglio::Atomo);
        cam.vai(Bersaglio::Atomo, r);
        cam.ln_r = cam.ln_r_target;
    }
}

/// The two lines at the top of the screen.
pub fn stato(lab: &Lab) -> (String, String) {
    match lab.tav.out() {
        Some(o) if o.pronto => o.righe,
        _ => ("TAVOLA PERIODICA · premi un elemento acceso (Tab apre la tavola)".into(), String::new()),
    }
}

/// The time bar of the table scene.
pub fn barra(ui: &mut egui::Ui, lab: &mut Lab) {
    let o = lab.tav.out();
    ui.label(RichText::new(o.map(|o| format!("{} · {:.0} ms", if o.nome.is_empty() { "nessun elemento" } else { o.nome.as_str() }, o.ms)).unwrap_or_default()).monospace().color(Color32::from_rgb(90, 220, 255)));
    if ui.add_enabled(lab.tav.idrogeno(), egui::Button::new(RichText::new(" Lancia (F) ").strong()).fill(Color32::from_rgb(40, 60, 110))).clicked() { lab.spara = true; }
}

#[cfg(test)]
mod test {
    use super::*;

    /// Layout: 118 distinct places, no element on the spacer row.
    #[test]
    fn tavola_posizioni() {
        let mut v: Vec<(usize, usize)> = (1..=118).map(posizione).collect();
        assert!(v.iter().all(|p| p.0 != 7 && p.0 < 10 && p.1 < 18));
        v.sort();
        v.dedup();
        assert_eq!(v.len(), 118);
        assert_eq!(posizione(26), (3, 7));
        assert_eq!(posizione(60), (8, 5));
        assert_eq!(formula(&[8, 1, 1]), "H2O");
        assert_eq!(formula(&[6, 1, 1, 1, 1]), "CH4");
    }

    /// The status file and the He export load; He is lit by its export, Li stays grey, a bench without an
    /// export never makes a molecule.
    #[test]
    fn tavola_carica() {
        let t = Tavola::carica();
        assert!(t.errori.is_empty(), "{:?}", t.errori);
        assert_eq!(t.accesso(1), Accesso::Scena("idrogeno".into()));
        let he = match t.accesso(2) { Accesso::Esportazione(i) => t.esportazioni[i].clone(), a => panic!("He: {a:?}") };
        assert_eq!((he.n_el, he.n, he.livello.as_str()), (2, 16000, "esatto"));
        let r: f64 = (0..he.n).flat_map(|i| (0..2).map(move |k| (i, k))).map(|(i, k)| { let p = he.el(i, k); (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt() }).sum::<f64>() / (2 * he.n) as f64;
        assert!((r / 0.929472 - 1.0).abs() < 0.01, "<r> = {r}");
        assert_eq!(t.accesso(3), Accesso::No);
        assert!(t.ioni.iter().any(|i| i.z == 60));
        assert_eq!(t.costruisci(&[8, 1, 1]), None);
        let h2 = t.costruisci(&[1, 1]).expect("H2 export");
        assert_eq!((t.esportazioni[h2].nome.as_str(), t.esportazioni[h2].simmetria), ("H2", Simmetria::AssialeZ));
        // the Coulomb hole: with one electron held on +z at 0.35-0.75 a0 the other is farther from it than without
        let mut n = Nuvola::new(he);
        let mut d = [0.0f64; 2];
        for (j, buca) in [false, true].into_iter().enumerate() {
            let mut s = 0.0;
            let mut c = 0;
            for _ in 0..400 {
                let sl = n.pesca(buca);
                let k = sl.rif.unwrap_or(0);
                let a = applica(&sl.rot, n.e.el(sl.i, k));
                let b = applica(&sl.rot, n.e.el(sl.i, 1 - k));
                if buca { assert!(a[2] > 0.34 && a[0].abs() < 1e-9 && a[1].abs() < 1e-9); }
                s += ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt();
                c += 1;
            }
            d[j] = s / c as f64;
        }
        assert!(d[0] > 1.0 && d[1] > 0.8, "{d:?}");
    }
}
