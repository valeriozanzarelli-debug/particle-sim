//! EXPORTS of atoms and molecules computed by the studies (thread "Debiti aperti" and others), in ONE generic
//! format so a new element or molecule lights up in the periodic table just by adding its folder.
//!
//! A folder `assets/esportazioni/<nome>/` holds:
//!   manifest.json   — what it is (see `docs/ATOMO.md` §10.2 for every field)
//!   campioni.f32    — n configurations × n_elettroni × 3 float32 little endian (bohr): each row is ONE
//!                     "photograph" of all the electrons together, drawn from |Ψ(r₁…r_N)|² by the study.
//! Optional extra search folder: env `ATOMO_ESPORTAZIONI` (e.g. /mnt/project-files/debiti/export).
//!
//! Nothing here computes physics: it reads what a study exported and keeps its tier ("esatto" /
//! "approssimato") and its declared PASS/FAIL so the screen can show them.

use std::path::{Path, PathBuf};

use serde_json::Value;

#[derive(Clone, Debug, PartialEq)]
pub struct Esame { pub id: String, pub cosa: String, pub pass: bool }

/// Which random rotations leave the exported state unchanged (used to draw "fresh" photographs from a few
/// thousand samples): a closed-shell S atom is spherical, a linear molecule along z is axial.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Simmetria { Sferica, AssialeZ, Nessuna }

#[derive(Clone, Debug)]
pub struct Esportazione {
    pub cartella: PathBuf,
    pub nome: String,
    pub simbolo: String,
    /// "atomo" or "molecola"
    pub tipo: String,
    /// atoms: Z of the nucleus (molecules: 0)
    pub z: u32,
    pub carica: i32,
    /// "esatto" or "approssimato"
    pub livello: String,
    pub metodo: String,
    pub fonte: String,
    /// (Z, position in bohr)
    pub nuclei: Vec<(u32, [f64; 3])>,
    pub n_el: usize,
    /// n × n_el × 3, bohr
    pub campioni: Vec<f32>,
    pub n: usize,
    pub simmetria: Simmetria,
    pub energia: Option<f64>,
    pub riferimento: Option<f64>,
    pub esami: Vec<Esame>,
    pub note: Vec<String>,
}

impl Esportazione {
    pub fn molecola(&self) -> bool { self.tipo == "molecola" }
    /// Electron e of configuration i (bohr).
    pub fn el(&self, i: usize, e: usize) -> [f64; 3] {
        let k = (i * self.n_el + e) * 3;
        [self.campioni[k] as f64, self.campioni[k + 1] as f64, self.campioni[k + 2] as f64]
    }
    /// Sorted list of the nuclear charges (what "costruisci molecole" compares).
    pub fn composizione(&self) -> Vec<u32> {
        let mut v: Vec<u32> = self.nuclei.iter().map(|n| n.0).collect();
        v.sort_unstable();
        v
    }
    pub fn n_pass(&self) -> (usize, usize) {
        let p = self.esami.iter().filter(|e| e.pass).count();
        (p, self.esami.len() - p)
    }
}

pub fn esami_da(v: &Value) -> Vec<Esame> {
    v.as_array().map(|a| a.iter().map(|e| Esame {
        id: e["id"].as_str().unwrap_or("?").to_string(),
        cosa: e["cosa"].as_str().unwrap_or("").to_string(),
        pass: e["PASS"].as_bool().unwrap_or(false),
    }).collect()).unwrap_or_default()
}

fn testo(v: &Value) -> String { v.as_str().unwrap_or("").to_string() }

/// Read one export folder. Errors are returned as text (shown in the panel, never hidden).
pub fn leggi(cartella: &Path) -> Result<Esportazione, String> {
    let man = std::fs::read_to_string(cartella.join("manifest.json")).map_err(|e| format!("{}: {e}", cartella.display()))?;
    let m: Value = serde_json::from_str(&man).map_err(|e| format!("{}: manifest.json: {e}", cartella.display()))?;
    if !m["formato"].as_str().unwrap_or("").starts_with("particle-sim/esportazione-v1") {
        return Err(format!("{}: formato sconosciuto {:?}", cartella.display(), m["formato"]));
    }
    let n_el = m["n_elettroni"].as_u64().ok_or("n_elettroni mancante")? as usize;
    let file = m["campioni"]["file"].as_str().unwrap_or("campioni.f32");
    let byte = std::fs::read(cartella.join(file)).map_err(|e| format!("{}: {file}: {e}", cartella.display()))?;
    if n_el == 0 || byte.len() % (12 * n_el) != 0 { return Err(format!("{}: {file}: {} byte non multiplo di 12·{n_el}", cartella.display(), byte.len())); }
    let campioni: Vec<f32> = byte.chunks_exact(4).map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]])).collect();
    let n = campioni.len() / (3 * n_el);
    let nuclei = m["nuclei"].as_array().map(|a| a.iter().map(|x| {
        let p = &x["pos"];
        (x["Z"].as_u64().unwrap_or(0) as u32, [p[0].as_f64().unwrap_or(0.0), p[1].as_f64().unwrap_or(0.0), p[2].as_f64().unwrap_or(0.0)])
    }).collect()).unwrap_or_default();
    let simmetria = match m["simmetria"].as_str().unwrap_or("") { "sferica" => Simmetria::Sferica, "assiale_z" => Simmetria::AssialeZ, _ => Simmetria::Nessuna };
    Ok(Esportazione {
        cartella: cartella.to_path_buf(),
        nome: testo(&m["nome"]), simbolo: testo(&m["simbolo"]), tipo: testo(&m["tipo"]),
        z: m["Z"].as_u64().unwrap_or(0) as u32, carica: m["carica"].as_i64().unwrap_or(0) as i32,
        livello: testo(&m["livello"]), metodo: testo(&m["metodo"]), fonte: testo(&m["fonte"]),
        nuclei, n_el, campioni, n, simmetria,
        energia: m["energia_Eh"].as_f64(), riferimento: m["riferimento_Eh"].as_f64(),
        esami: esami_da(&m["esami"]),
        note: m["note"].as_array().map(|a| a.iter().filter_map(|s| s.as_str().map(String::from)).collect()).unwrap_or_default(),
    })
}

/// The assets folder: `assets/` next to where the program runs, else the one of the source tree.
pub fn cartella_assets() -> PathBuf {
    let p = PathBuf::from("assets");
    if p.join("tavola_stato.json").exists() { p } else { PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("assets") }
}

/// Every export found (sorted by name), and the errors of the folders that could not be read.
pub fn tutte() -> (Vec<Esportazione>, Vec<String>) {
    let mut cartelle = vec![cartella_assets().join("esportazioni")];
    if let Ok(x) = std::env::var("ATOMO_ESPORTAZIONI") { cartelle.push(PathBuf::from(x)); }
    let (mut ok, mut err) = (vec![], vec![]);
    for c in cartelle {
        let Ok(rd) = std::fs::read_dir(&c) else { continue };
        let mut dirs: Vec<PathBuf> = rd.filter_map(|e| e.ok().map(|e| e.path())).filter(|p| p.join("manifest.json").exists()).collect();
        dirs.sort();
        for d in dirs {
            match leggi(&d) {
                // the same name twice: the first folder (the repository) wins
                Ok(e) => if !ok.iter().any(|x: &Esportazione| x.nome == e.nome) { ok.push(e) },
                Err(s) => err.push(s),
            }
        }
    }
    (ok, err)
}
