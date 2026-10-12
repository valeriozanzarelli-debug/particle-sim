# Contesto per la sessione in cloud: indice di riferimento

Scritto il 12 ottobre 2026 (sera) per una sessione Claude Code nuova, che non ricorda niente della
conversazione sul portatile. Non è un racconto, è un indice: dove sta cosa, cosa è stato misurato,
cosa è fallito, cosa resta da pagare.

- I numeri sono copiati dai documenti degli studi (`docs/*.md`). Se un numero non è qui, va letto
  nel documento: non va ricostruito a memoria.
- Quello che viene dallo stato di git o dei file, e non da un documento, è marcato **(osservato)**.
- Lo stato esatto di ogni cantiere al momento del passaggio è in `docs/PASSAGGIO_CLOUD.md`.

---

## 1. Cosa leggere per primo

| # | File | Perché |
|---|------|--------|
| 1 | `CLAUDE.md` | regole vincolanti, modo di lavorare, firma dei commit (viene letto in automatico) |
| 2 | `docs/PIANO.md` | obiettivo, regole, tabella di stato, 4 cantieri aperti (§4), prossimi passi (§5), idee (§6), setup (§7) |
| 3 | `docs/CONTESTO_CLOUD.md` (questo file) | mappa del repo, indice di tutti gli studi con i numeri, debiti, strumenti, lezioni |
| 4 | `docs/PONTE.md` | come usare il portatile di Valerio come laboratorio GPU |

Poi, secondo il lavoro:
- `docs/PASSAGGIO_CLOUD.md`: stato di ogni cantiere al passaggio.
- `docs/MOLTI_ELETTRONI_PIANO.md`: piano già giudicato per gli atomi a più elettroni (da He a Ne).
  **Non è stato eseguito niente.** I numeri "development observations" vanno riverificati nella
  pre-registrazione. I risultati andranno in `docs/MOLTI_ELETTRONI.md`.
- `docs/PROMPT_CLOUD.md`: il prompt di avvio di questa sessione.
- Il documento dello studio su cui lavori (`docs/<STUDIO>.md`) e la sua pre-registrazione
  (`lab/out/<studio>_prereg.json`).

### Stato in poche righe

- **Una sola strada:** la regola zig-zag (ogni pallina corre a c, poi ruota; l'angolo θ = m c a/ħ
  è la massa, la fase di Coulomb è la carica). Si sale di scala solo con l'**agglomeratore
  automatico** (blocca, impara, controlla, ripeti).
- **Idrogeno dalla regola** (forma radiale): 1s 13,60575 eV contro Dirac 13,60587.
  **Agglomerato a scala di Bohr:** 1s −13,59967 eV contro −13,60569 (−0,044 %).
- **Quattro cantieri aperti** (PIANO §4): agglomeratore sulla zig-zag, colla nucleare, reticolo
  dell'universo, orbitali del Laboratorio Atomo.
- **Luce + atomo accoppiati** (`src/atomo/accoppiato.rs`, in `docs/ATOMO.md`): v2 E2v2 PASS
  (−0,6 %), E1v2 FAIL (+7,1 %). **(osservato)** La tabella di stato di PIANO.md non elenca ancora
  questo studio (commit ed14b66, 635839f, c7271d7, d782a44).
- **(osservato)** Ultimo commit alla stesura: d782a44 (luce accoppiata dal vivo nella scena
  dell'atomo). Non committati: `lab/colla_figure.py`, `lab/colla_summary.py`,
  `lab/out/colla_units.json`. `docs/COLLA.md` era non tracciato e poi è stato aggiunto all'indice
  di git. Se nel clone manca `docs/COLLA.md`, il cantiere Colla va ripreso dalla sua
  pre-registrazione (commit 78d76a7), come dice PIANO §4.
- **(osservato)** Il ponte ha fatto un giro completo:
  - file `ponte/risultati/esempio_ping/esito.json`: `ok: true`, 2026-10-12T22:47Z, revisione 88ff896;
  - GPU elencate: NVIDIA GeForce GTX 1650 Ti e Intel UHD Graphics; 8 core; Python 3.14.2.

---

## 2. Mappa del repo

Setup (PIANO §7):

    git clone https://github.com/valeriozanzarelli-debug/particle-sim.git
    pip install numpy scipy torch numba ripser gudhi matplotlib
    cargo build --bin atomo        # profilo dev con opt-level 3
    cargo test --bin atomo

### 2.1 Binari Rust (`Cargo.toml`)

Dipendenze principali: bevy 0.15, bevy_egui 0.31, wgpu 23.

| bin | sorgente | cosa | stato |
|-----|----------|------|-------|
| `atomo` | `src/atomo/` | **Laboratorio Atomo**: luce e materia dalla regola zig-zag. Contiene: scena del fotone 3-D; atomo vivo Nd59+ (Z=60) con protone fatto da un nodo di Skyrme; orbitali n=1..4; surfisti e misura; motore accoppiato luce-atomo | **corrente, laboratorio principale** |
| `zigzag` | `src/zigzag/` | viewer della regola (modi luce, elettrone, idrogeno). `walk.rs` è la regola in Rust, forma `Radial` compresa; anche `atomo` la usa via `#[path]` | corrente (riferimento della regola, 8 test) |
| `fotone` | `src/photon/` | Laboratorio Fotone: regola a due ruoli (`pilot.rs`), `base3d.rs`, `gpu.rs`, `timeline.rs`, `haze.rs`, `entangle.rs` | storico; resta riferimento di **stile grafico** (PIANO regola 9) |
| `molecole` | `src/molecole/` | DMC per palline: H, H2+, H2, He | storico |
| `scontro` | `src/scontro/` | viewer della collisione fotone-elio | storico |
| `compton` | `src/compton/` | scattering Compton | storico |
| `particle_sim` | `src/main.rs`, `src/simulation/`, `src/rendering/`, `src/config.rs` | primo prototipo: Primordial Particle System 3-D, calcolo con wgpu | legacy |

**(osservato)**
- `src/adaptive/` non è usato da nessun binario.
- `src/rendering/camera.rs` è condiviso da `zigzag` e `fotone`.

Note pratiche:
- **Release:** la build con LTO dura circa 35 min.
- **Cartella di lavoro:** i programmi grafici vanno lanciati dalla cartella del repo (gli shader
  sono in `assets/`).
- **Server senza schermo:** i programmi grafici non si aprono. Su Linux servono `libasound2-dev`,
  `libudev-dev`, `pkg-config`.
- **Eseguibili:** `Laboratorio Atomo.exe` e `Laboratorio Fotone.exe` stanno fuori dal repo, nella
  cartella madre. Si rifanno con `cargo build` e si copiano.
- **Test grafici e foto:** si fanno solo sul portatile, tramite il ponte (vedi §5).

### 2.2 `src/atomo/` (laboratorio principale)

| file | ruolo |
|------|-------|
| `fotone3d.rs`, `fotscena.rs`, `fotdisegno.rs` | scena del fotone (luce 3-D della regola) |
| `zz.rs`, `zzscena.rs`, `render_zz.rs` | atomo vivo sulla regola zig-zag |
| `orbitali.rs` | orbitali n=1..4 dalla forma radiale della regola |
| `skyrme.rs` | nodo di Skyrme (il protone) |
| `misura.rs` | misura e surfisti |
| `accoppiato.rs` | accoppiamento luce-atomo (Maxwell-Dirac semiclassico) |
| `tests.rs` | test (veloci + `--ignored`) |
| `main.rs`, `ui.rs`, `camera.rs`, `render.rs`, `rng.rs`, `balls.wgsl` | supporto |
| `../zigzag/walk.rs`, `../zigzag/cplx.rs` | regola e numeri complessi, inclusi con `#[path]` |

Opzioni da riga di comando (ATOMO.md):
- scena: `--scena fotone|atomo|orbitali [--zoom nucleo] [--stati 3d | 2p+3d]`;
- fotone: `--qual 0..3 --w0 --dir 0..4 --polar h|v|c+|c- --lambda-m --strati tflEBsox`;
- foto senza finestra: `--foto file.png --frames N [--spara-al F] [--yaw a --pitch b] [--ritmo r]`.

Test:
- `cargo test --bin atomo`: 7 test veloci (F1a, F2, F4, regressione dell'apertura, livelli degli
  orbitali, atomo, nodo).
- `cargo test --bin atomo -- --ignored`: test lunghi o con FAIL registrato (`zigzag_risonanza`,
  direzioni oblique, `bench_fotone`, `accoppiato_e1_e2`, `accoppiato_v2`).

### 2.3 `lab/`: studi in Python, per prefisso

| prefisso | documento | sulla strada |
|----------|-----------|--------------|
| `zz_` | ZIGZAG.md (regola congelata in `zz_core.py`) | sì |
| `agg_zz_` | AGGLOMERATO_ZIGZAG.md | sì |
| `reticolo_` | RETICOLO_UNIVERSO.md | sì |
| `zoom_` | ZOOM.md (algoritmo generico dell'agglomeratore) | sì |
| `colla_` | COLLA.md (cantiere aperto) | sì |
| `knot_`, `knot2_` | NODO_3D.md | sì |
| `nuclei_` | NUCLEI_DA_NODI.md | sì |
| `zoomknot_` | ZOOM_NODO.md | sì |
| `redo_` | IDEE_RIPRESE.md | sì |
| `box_` | COPPIE_IN_SCATOLA.md | sì |
| `multi_e_` | MOLTI_ELETTRONI_PIANO.md (file da creare, non esistono ancora) | sì |
| `press_` | PRESSIONE.md | no |
| `glue_` | GLUEBALL.md | no |
| `kk_` (+ `kk_shape_`, `kk_speed_`) | CARICA_NASCOSTA.md, FORMA_NASCOSTA.md | no |
| `geo_` | FORMA_DALLE_PALLINE.md | no |
| `relax_` | EQUILIBRIO_QUANTISTICO.md | no |
| `pre_` | PRE_ESPERIMENTI.md | no |
| `big_` | GEOMETRIA_AUTOCOSTRUITA.md | no |
| `orb_` | ORBIFOLD.md | no |
| `eye_`, `flow_` | OCCHIO_E_FLUSSO.md | no |
| `kz_` | SINCRONIA.md | no |
| `fluid_` | FLUIDO_IN_NEGATIVO.md | no |
| `trap_` | LUCE_INTRAPPOLATA.md | no |
| `search_` | RICERCA_REGOLE.md | no |
| `dirac_` | SURFISTA_DIRAC.md | no |
| `h2vib_` | H2_VIBRAZIONE.md | no |
| `identical_` | REGOLA_PALLINE_IDENTICHE.md | no |
| `base_` | REGOLA_BASE.md | no |
| `multi_` (`multi_core`, `multi_gra`, …) | MULTIFOTONE.md | no |
| `img_` | IMMAGINI.md | no |
| `compton_check.py` | COMPTON.md | no |
| `molecole_misura.sh` | MOLECOLE.md | no |
| `optics.py`, `exact_guide*.py`, `photon_rule.py`, `neighbour_study.py`, `validate.py`, `isolated_ball.py`, `rules.py`, `pilot*.py` | REGOLA_FOTONE.md | no (`optics.py` resta riferimento esatto RS) |

Attenzione al prefisso: il piano a più elettroni usa `multi_e_`, mentre `multi_` è già usato da
MULTIFOTONE. Non confonderli.

### 2.4 `assets/`

| file | cosa | prodotto da |
|------|------|-------------|
| `assets/atomo/agglomerato_zz.json` | livelli agglomerati della regola zig-zag (catena A, livello B = 256 per il porto in Rust) | `lab/agg_zz_export.py` |
| `assets/atomo/reticolo_spec.json` | formato del reticolo sparso multilivello per Rust (con tutti i 40 casi R2) | `lab/reticolo_spec.py` |
| `assets/atomo/risultati.txt` | numeri archiviati del carbonio DMC (codice vecchio al commit a84692d) | storico |
| `assets/molecole/risultati.txt` | i 33 valori DMC usati anche da H2_VIBRAZIONE | `cargo run --bin molecole -- --misura …` |
| `assets/shaders/photon_grains.wgsl` | shader del port GPU (GPU_PORT) | |
| `assets/shaders/pps_compute.wgsl` | shader di `particle_sim` | |

Gli shader delle palline (`balls.wgsl`) stanno in `src/<bin>/`.

### 2.5 `lab/out/` e file grandi

- **Ignorato da git.** `.gitignore`: `/target`, `*.npy`, `__pycache__/`, `lab/out/`. Figure e JSON
  di risultato vanno aggiunti con `git add -f`.
- **(osservato)** Sono già tracciati: 24 pre-registrazioni `*_prereg.json`, molti log, alcuni
  `.npz`/`.pkl` piccoli (es. `agg_zz_3d_A.pkl`, `agg_zz_3d_B.pkl`) e i tre script
  `lab/out/agg_zz_run*.sh`.
- **Non sono nel repo e vanno rigenerati:**
  - `.pt` (checkpoint di Colla), `.npz` dei campi grandi, `.npy`;
  - `lab/out/glue_raw/` (~420 MB);
  - `lab/out/multi_cache/` (~95 MB);
  - `lab/out/zoom_ising_mc.npz`;
  - i campi `.npz` e i checkpoint di Colla.
- **Log scritti di continuo:** non vanno tracciati, perché rompono `git pull --rebase --autostash`.

### 2.6 Altro

- `ponte/`: `ponte.py` (gira sul portatile), `richieste/`, `risultati/`. Vedi PONTE.md.
- `docs/fig/`: figure dei documenti. `docs/RELAZIONE.html`, `docs/tre_punti.html`,
  `docs/universo.html` sono pagine HTML presenti nel repo.

---

## 3. Indice degli studi

Stato: *validato*, *parziale*, *in corso*, *fallito*, *storico*, *strumento*, *visione*, *piano*.
"Strada" vuol dire: serve alla strada attuale (zig-zag + agglomeratore + cantieri di PIANO).

| # | File | Cosa | Stato | Numeri chiave (nostro vs riferimento) | Strada |
|---|------|------|-------|---------------------------------------|--------|
| 1 | ZIGZAG.md | regola zig-zag (Dirac quantum walk) su luce, elettrone libero, idrogeno radiale; regola congelata il 2026-10-09 in `lab/zz_core.py` | parziale (Python 22 PASS / 18 FAIL su 40; Rust 8/8) | 1s 13,60575 eV vs Dirac 13,60587; struttura fine 10,94367 GHz vs misurata 10,96905 (−0,23 %, manca QED); A(2p→1s) 6,249e8 vs 6,2649e8 s⁻¹ (−0,25 %); doppia fenditura λ=8a corr 0,595 (FAIL), λ=24a 0,9988 | sì |
| 2 | AGGLOMERATO_ZIGZAG.md | agglomeratore automatico sulla zig-zag: Schrödinger + (v/c)², 1s a scala di Bohr su 192³, zoom ibrido vicino al protone, piano a più elettroni | parziale | 1s −13,59967 vs −13,60569 eV (−0,044 %, PASS); 2p −0,29 % PASS; 2s +1,003 % FAIL (soglia 1 %); celle×tick 4,3e9 volte meno della regola fine | sì |
| 3 | RETICOLO_UNIVERSO.md | reticolo sparso multilivello h_k = 2^k h_0, stessa regola di luce su ogni livello, cuciture tra livelli | parziale | velocità −0,162 % sui livelli 0, 1, 10, 20 (PASS); passaggio fine→grosso T 0,99993 a 16 celle/λ; mondo a 21 livelli 717 568 celle, 34,7 MB; radio al livello 20 non fatta | sì |
| 4 | ATOMO.md | Laboratorio Atomo (Rust/Bevy) + motore accoppiato luce-atomo | parziale | orbitali 1s1/2 +39,9 %, 4f7/2 +1,3 % vs Dirac; apertura del fascio −20,9 % (8 celle/λ), −6,0 % (16); E2v2 −0,6 % PASS; E1v2 +7,1 % FAIL | sì |
| 5 | ZOOM.md | algoritmo generico dell'agglomeratore (blocca→impara→passo→controlla→ripeti) su catena di molle, Ising 2-D critico, fasci di luce | strumento | Ising y_h 1,8738 vs 1,875; η 0,252 vs 0,25; catena livello 3 K_eff +1,3 %; onde corte falliscono oltre ~0,4 del Nyquist | sì |
| 6 | MOLTI_ELETTRONI_PIANO.md | piano giudicato (3 progetti, 3 giudici) per atomi da He a Ne dalla regola agglomerata | piano, nulla eseguito | nessun risultato; numeri di sviluppo da riverificare | sì |
| 7 | NODO_3D.md | nodo di Skyrme B=1 sul reticolo 3-D: una sola massa, Lorentz, annichilazione; parte 2 sugli stencil (c4) | parziale | E/12π² +0,77 % vs 1,231446; c4 v=0,6 h=0,15 velocità −2,2 % PASS; v=0,8 FAIL; annichilazione: energia uscita 98,4 % / 99,1 % (S4) | sì |
| 8 | NUCLEI_DA_NODI.md | nuclei B=1..4 (toro, tetraedro, cubo) dai nodi | parziale | E/B entro 0,26–0,41 % della letteratura; le forme nascono da sole (N2); legame deutone 72 MeV vs 2,22 reale (×32), ⁴He 310 vs 28,3 (×11) | sì |
| 9 | ZOOM_NODO.md | livelli annidati di agglomerati sul nodo, regola di blocco imparata | validato | p(0) −0,08 %, E −0,16 %, Λ +0,9 % con 195 112 palline (24 389 fini); griglia uniforme: p(0) +67,7 % | sì |
| 10 | COLLA.md | la "colla" mancante del legame nucleare: massa del pione, quantizzazione, punto zero, GHS, termine sestico | in corso (documento non finito) | massa del pione: deutone 68 MeV; quantizzazione rigida: 160 MeV vs Braaten–Carson 158; + punto zero: 78 MeV (35× il reale); sestico c₆=200 **alza** il legame (C3b-i FAIL) | sì |
| 11 | IDEE_RIPRESE.md | R1 Schwinger/Sauter sulla zig-zag 1+1D; R2 flusso quantizzato contro le isole; R3 surfista deterministico | parziale | R1a reticolo/Schwinger 0,98–1,04 (PASS); Sauter 0,992–0,996 (PASS); R2 FAIL; R3 PASS ma fragile (controllo fine FAIL) | sì |
| 12 | COPPIE_IN_SCATOLA.md | γ → e⁺e⁻ tra due lastre: zig-zag 2+1D in una striscia | parziale | m_eff segue p d = π − arctan(p/m) entro l'1 % per d ≥ 32; crossover d* = 2,25 lunghezze Compton (PASS); scatole strette contate male (A1a FAIL) | sì |
| 13 | PIANO.md | piano maestro e stato | in corso | idrogeno 13,60575 eV; Schwinger 0,3 %, Sauter 0,5 %; Lorentz da fluido 1,5 %; pressione nel protone 0,605 vs 0,586 fm; Kibble–Zurek 0,507 vs 0,573 (fallito) | sì |
| 14 | PONTE.md | ponte cloud↔portatile (portatile = laboratorio GPU) | strumento | controllo ogni 90 s; timeout 1 ora per azione; id fino a 64 caratteri | sì |
| 15 | VISIONE.md | nota iniziale del "Blender atomico" (30 set): livelli di zoom, osservatore che accetta una regola solo sotto una soglia d'errore | visione | nessun numero misurato | sì (antenato dell'agglomeratore) |
| 16 | VISIONE_SISTEMA.md | visione di prodotto + aggiunta del 12 ott: motore tipo Unreal, qualità grafica = livelli di agglomerato | visione | cita "0,08% con 7 livelli" (ZOOM_NODO) | sì |
| 17 | PRESSIONE.md | pressione, taglio, D-term dentro il nodo B=1 vs JLab | parziale (13 PASS, 3 FAIL) | r₀ 0,605 fm vs 0,586 JLab; p(0) 6,5·10³⁴ Pa vs tripolo JLab 1,9·10³⁵ (0,34, FAIL previsto); nucleo sul reticolo h=0,2 sbagliato del 57 % | no |
| 18 | GLUEBALL.md | SU(2) puro 2+1D: confinamento e glueball | parziale | a√σ entro 0,6–1,6 %; m(0++)/√σ nel continuo 4,93(10) vs 4,737; β=12 0++ +5,8 % (FAIL); 0++* FAIL | no |
| 19 | LUCE_INTRAPPOLATA.md | la materia è luce intrappolata? | storico | la luce nel vuoto lineare scappa nel tempo dei raggi; skyrmion M/Bogomolny 1,2315 vs 1,232; g_A 0,640 vs 1,275 | no |
| 20 | RICERCA_REGOLE.md | ricerca automatica (CMA-ES) di regole che fanno particelle dalla luce | storico | CMA-ES vs caso: 28 vs 2 ammissibili; finalista f_EM 0,38 %; zero coincidenze entro il 5 % (a caso se ne aspettano il 61 %) | no |
| 21 | CARICA_NASCOSTA.md | carica = avvolgimento su un cerchio nascosto (Kaluza-Klein sulla zig-zag) | parziale, chiuso | carica intera conservata (deriva 0); massa sbagliata di ~10²¹ (10¹⁷ GeV) | no |
| 22 | FORMA_NASCOSTA.md | dimensione nascosta come forma (ovale, triangolo, stadio): masse = modi | fallito | m₂/m₁ = 1,40–1,81 vs 206,768; nessuna famiglia a 0–1 parametri passa | no |
| 23 | FORMA_DALLE_PALLINE.md | ricostruire lo spazio dall'informazione mutua del vuoto; inversa dalle masse | parziale, chiuso | raggio nascosto / lunghezza 0,0792 vs 0,0819; la regola raddoppia lo spazio (due fogli); Koide: 1 su 5 | no |
| 24 | EQUILIBRIO_QUANTISTICO.md | regola di Born come equilibrio raggiunto intorno ai vortici (nel continuo) | parziale | M=25: τ 2,4; senza vortici nessun rilassamento; 3b 0,864 → 0,980 vs ideale 0,9976 (FAIL) | no |
| 25 | PRE_ESPERIMENTI.md | P1–P5: zig-zag con pesi sui link e surfisti stocastici | parziale | P1c mescolamento finto (M=1 veloce come M=16); P2 falsi positivi 60–100 %; P4 U(k+(π,π)) = U(k) a 5e-16 | no |
| 26 | GEOMETRIA_AUTOCOSTRUITA.md | lettore riparato (R1), mescolamento finto (R2), link rinforzati (B) | fallito | R1 PASS (0/25 falsi buchi); R2 FAIL (rapporti 0,65 e 1,46 vs ≥ 3); B: isole, 0 buchi | no |
| 27 | ORBIFOLD.md | palline su orbifold imposti a mano (K3, T2/Z3, T6/Z3) | parziale | O1–O3 PASS (conteggi uguali ai libri, 0 parametri continui); "3 famiglie" solo scegliendo M e settore (17 scelte su 72) | no |
| 28 | OCCHIO_E_FLUSSO.md | occhio (cavità, punte coniche) e flusso di Ollivier-Ricci | fallito | punte T4/Z2 16/16; cavità di T2/Z3 0,230 < soglia 0,237 (FAIL); il flusso non appiattisce (3,0–3,8× invece di 10×) | no |
| 29 | SINCRONIA.md | metronomi = Kibble–Zurek nell'XY 3-D | fallito | α = 0,507 [0,493; 0,523] vs 0,573; ν = 0,677 ± 0,018 vs 0,6717 (PASS) | no |
| 30 | FLUIDO_IN_NEGATIVO.md | particelle come difetti in un fluido di palline | storico | vortice GP A = 2,078 vs 2,070; kink: Lorentz entro l'1,45 % fino a 0,8c; la massa del buco non è un numero solo (FAIL) | no |
| 31 | REGOLA_FOTONE.md | regola del fotone a due ruoli ("onda pilota"), porto in Rust, poi base3d | storico | fascio libero 0,2 %; doppia fenditura corr 0,995; Rust base3d Airy +1,7 %; 29 test di `fotone` PASS | no |
| 32 | REGOLA_BASE.md | regola base congelata (η=2, ε=10⁻³, δ=0,02 λ, Θ=90°); v2 con zona di bordo λ/2 | storico | v2: 56 righe tutte PASS; Young corr 0,9962; Arago 0,0294 vs 0,0293 | no (la sua doppia fenditura è il banco di prova della zig-zag) |
| 33 | ESPERIMENTI_FOTONE.md | catalogo di 96 esperimenti col fotone | storico | 96 schede; 8 validate al commit 83194fa | no |
| 34 | MULTIFOTONE.md | n fotoni + livelli di zoom della luce | storico | 96/96 PASS; CHSH S = 2,8322 ± 0,0141 vs 2,8284; raggi validi all'1 % per N_F > 4400 | no |
| 35 | FOTONE_UI.md | interfaccia del Laboratorio Fotone, tempo reversibile | strumento (storico) | keyframe ~190 MB su un budget di 192 MB | no |
| 36 | GPU_PORT.md | port wgpu della fisica a grani classica | strumento (storico) | CPU vs GPU corr 0.9999995; 10⁶ grani: GPU 11–13 ms/tick vs CPU 16.5 ms/passo | no |
| 37 | IMMAGINI.md | come i fotoni fanno un'immagine (sensore, lente, risoluzione) | parziale | primo anello scuro 3,169 ± 0,028 λ vs 3,142 λ; R5b FAIL (esponente −2,72) | no |
| 38 | REGOLA_PALLINE_IDENTICHE.md | fotone fatto di N palline identiche (Many-Interacting-Worlds) | fallito | fascio libero −0,05 % (N=1000); doppia fenditura al massimo ~0,60 vs ~0,97 richiesto | no |
| 39 | SURFISTA_DIRAC.md | elettrone con la vecchia regola, surfista Bohm–Dirac | storico (34 PASS, 2 FAIL) | g = 1,9962 (σ=8); struttura fine 10,94369 GHz; surfista 1s con \|v\| = α c sin θ | no (la formula del surfista è il passo 3 della zig-zag) |
| 40 | COMPTON.md | Compton su elettrone libero | storico | 9/9 test Rust, 19/19 controlli Python; Klein–Nishina è un INGRESSO | no |
| 41 | MOLECOLE.md | DMC per palline: H, H2+, H2, He | storico | H −13,60567 vs −13,605693 eV; H2 D0 4,4721 vs 4,47813 eV; riga He 58,437 vs 58,4334 nm | no |
| 42 | H2_VIBRAZIONE.md | diagnosi del −1,4 % sulla vibrazione di H2 | storico | causa: il fit a 5 esponenziali; corretto ω_e 4393 ± 36 vs 4401,21 (−0,18 %) | no |

### 3.1 Strada attuale: zig-zag + agglomeratore + reticolo + Laboratorio Atomo

**ZIGZAG.md**
- **FAIL:**
  - L1a: dispersione lungo la diagonale −3,53 % a λ=8a. La tolleranza del 2 % era troppo stretta
    per un errore che va come (a/λ)².
  - L5: doppia fenditura a λ=8a, corr 0,595. La diffrazione è più debole di un fattore ka·cot(ka);
    serve λ ≥ 16a.
  - H1: livelli s 7,0e-5 contro una tolleranza di 1e-5 (errore di contatto, θ² ln θ).
  - A1: regola cubica 3-D a Z=60, −11,6 % dopo l'estrapolazione. Per stare entro lo 0,5 % serve
    a ≤ 0,05 λ_C, cioè circa 480³ celle.
  - A2: il surfista 3-D non è stazionario.
  - S3: 2s dall'ascolto −3,0 %.
  - F5: fattore g a θ=0,5 vale 1,820.
- **Debiti:**
  - Ingressi messi a mano: c, ħ, spin, m_e (= θ), α = 1/137,035999, Z=1 puntiforme, massa del
    protone (solo nella massa ridotta).
  - L'idrogeno vero è fatto solo in forma radiale (in 3-D servirebbero ~3e11 celle).
  - I livelli vengono dal risolvente; confrontati con l'ascolto vero solo per l'1s.
  - QED assente: Lamb, g−2, i 25 MHz della struttura fine.
  - La regola non è derivata: è il Dirac quantum walk noto.
- **Rifare:**
  - `python lab/zz_core.py` (autotest, 2 s);
  - `timeout 590 python lab/zz_hydrogen.py` (249 s), `zz_hydrogen_s.py` (46 s), `zz_light.py` (72 s),
    `zz_slit_conv.py` (318 s), `zz_free.py` (284 s), `zz_atom3d.py` (257 s);
  - `cargo test --bin zigzag`;
  - `cargo run --bin zigzag -- --modo luce|elettrone|idrogeno`.

**AGGLOMERATO_ZIGZAG.md**
- **PASS:**
  - A1a: la regola grossa riproduce 1000 tick fini.
  - A1b: m_kin/m 1,00018 … 0,99995.
  - A1c: a B = 32, 64, 192.
  - A2a: 1s −0,044 %.
  - A2b: 2p −0,29 %.
  - A3a: ibrido −0,071 %.
  - A3c: costo 7,2e4 volte più basso.
- **FAIL:**
  - A1c a B=384: coefficiente di k⁴ 1,50, cioè sbagliato del 50 %.
  - A2b: 2s +1,003 %. La cella a0/8 è troppo grossa per la cuspide.
  - A2c: struttura fine non disponibile, perché il kernel imparato non ha spin (1e-10).
  - A3b: Zitterbewegung non vista. Il picco a 2,2008 è la 3ª armonica del passo grosso, non 2E.
- **Catena 3-D A:** B=16 (a0/256), B=256 (8,56 λ_C = a0/16), B=512 (a0/8).
- **Debiti:**
  - Spin-orbita assente al livello grosso.
  - Il 2s all'1 % richiede a0/16 su ±12 a0, cioè 384³.
  - Il potenziale grosso (media di −α/r sulla cella) è scritto a mano.
  - Manca il rinculo del protone.
  - I molti elettroni sono solo un piano.
  - Deviazioni dichiarate: onde piane invece di pacchetti; molte corse di sviluppo prima della
    pre-registrazione.
- **Passo dichiarato (§10): porto in Rust**
  - livello B=256: cella 8,5647 λ_C, τ = 8,5647 ħ/mc²;
  - G(k) per FFT dello stencil;
  - 128³ (±4 a0) dal vivo.
- **Rifare:**
  - `python lab/agg_zz_prereg.py` (pre-registrazione committata in b714d9e);
  - `python lab/agg_zz_1d.py` (~1 min);
  - `sh lab/out/agg_zz_run3d.sh` (~15 min), `sh lab/out/agg_zz_runatom.sh` (~20 min),
    `sh lab/out/agg_zz_runzoom.sh` (~15 min);
  - `python lab/agg_zz_export.py && python lab/agg_zz_figure.py`.

**RETICOLO_UNIVERSO.md**
- **PASS:**
  - R1: kernel = `zz_core` a 3e-14.
  - R2a: sugli assi; dopo la cucitura la velocità è 1,000000 c.
  - R2b: a N=16.
  - R2c: sugli assi.
  - R3: mondo a 21 livelli.
  - R3b e R4.
- **FAIL:**
  - R2a a 30°/45°: −0,81 % / −1,02 % (previsto; servono ≥ 32 celle/λ).
  - R2c a 30°: fedeltà 0,9979 a 16 celle/λ.
  - R3: radio al livello 20 non fatta (~15 ore in Python). Sostituita post hoc da R3b al livello 10.
  - Fantasma del gamma sui livelli 2–7.
- **Debiti:**
  - Motore cambiato dopo la pre-registrazione: cucitura v2. Le corse v1 sono in
    `lab/out/reticolo_v1/`.
  - Energia non esatta all'interfaccia: −0,12 % / +0,04 % a 30°. Serve una cucitura conservativa.
  - La restrizione fine→grosso deve diventare conservativa.
  - Solo luce: palline massive e agglomerati non ancora collegati.
  - 3-D provato solo all'interfaccia.
  - Python lento: 52 ms per tick del gamma. Il Rust deve essere 10–100 volte più veloce.
- **Stime per la Via Lattea:**
  - impulso radio λ = 1 m ovunque: 4,4e7 celle, 2 GB;
  - radio + gamma da 1 MeV con 40 livelli di camino: 2,5e8 celle, 12 GB.
- **Rifare:**
  - `python lab/reticolo_prereg.py`;
  - `python lab/reticolo_r2.py 2 3 4 6 8 12 16` (anche `--ang 30`, `--d3`);
  - `python lab/reticolo_r2a.py 0` e `30 45`;
  - `python lab/reticolo_r3.py ladder|both|frames`;
  - `python lab/reticolo_r4.py`;
  - `python lab/reticolo_spec.py`.

**ATOMO.md**
- **PASS:**
  - F1a: 0,99647 c.
  - F1b a 16 celle/λ: −0,085 %.
  - F2: elicità.
  - F4: 1,9e-13.
  - E2v2.
- **FAIL:**
  - F1b a 8 celle/λ: −0,58 %.
  - F3 apertura del fascio: −20,9 % / −6,0 %. La velocità di gruppo trasversale è q⊥/tan(ka) invece
    di q⊥/(ka).
  - Luce obliqua: 0,900 c lungo (1,1,0) e 0,876 c lungo (1,1,1), a 8 celle/λ.
  - Orbitali s troppo legati per colpa della cella: 1s +39,9 %.
  - Accoppiato v1: E1 −2,01 %; E2 con potenza negativa.
  - Accoppiato E1v2 +7,1 %. La riga vera è 1s→2p3/2 = 0,0778 più un 1s sporco del reticolo;
    il confronto non è stato cambiato dopo la corsa.
  - Doppia fenditura −12,7 %: scena tolta.
- **Debiti:**
  - Il fotone singolo è classico: non c'è il "clic".
  - ≥ 24 celle/λ a ≥ 30 tick/s non entra nel portatile (servono ~7 milioni di celle per tick;
    alla qualità massima si va a ~5 tick/s).
  - L'atomo è Nd59+ invece di idrogeno.
  - Accoppiamento semiclassico: uno stato eccitato puro non decade. Serve la luce quantizzata.
  - L'accoppiamento è K = 10 volte quello reale.
  - L'1s 3-D è mescolato con livelli del reticolo (0,8405 e 0,882).
  - Gli orbitali vengono dalla forma radiale (per n ≥ 3 servirebbero ~350³ celle).
  - I due numeri di Skyrme sono presi da Adkins–Nappi–Witten.
- **Costo per tick** (i5 a 4 core, build dev): ~9–14 ms a 96³, ~31 ms a 128³, ~41 ms a 144³,
  ~105 ms a 192³.
- **Rifare:** `cargo test --bin atomo` (veloci) e `-- --ignored` (lunghi). Le foto si fanno col ponte.

**ZOOM.md** (l'algoritmo)
- **FAIL:**
  - Onde corte: ogni regola grossa senza memoria sbaglia oltre ~0,4 del suo Nyquist (manca la
    memoria di Mori–Zwanzig).
  - La regola grossa è a lungo raggio (in 3-D circa 1150 vicini per super-pallina).
  - Luce: raggruppare butta via la fase.
  - Ising: nessun guadagno di velocità.
- **Debiti:**
  - Iperparametri scelti a mano: raggio 6,5 h, nodi 0,25 h, cluster 1,5 h, blocco k=2.
  - Descrittori fissi, senza termine di memoria.
- **Rifare:** `timeout 590 python lab/zoom_chain.py`, `zoom_ising.py` (scrive un `.npz` grande),
  `zoom_ising_rg.py`, `zoom_light.py`.

**MOLTI_ELETTRONI_PIANO.md** (dal file, osservato)
- **Milestone:**
  - M0: pre-registrazione `lab/multi_e_prereg.py`, controlli sulla regola, ECG, campo medio.
  - M1: He.
  - M2: Li.
  - M3: Be.
  - M4 (opzionale): da B a F.
  - M5: Ne.
  - M6: verifica 3-D R7 e documento.
- **Nuovo codice:** binario Rust `elettroni` in `src/elettroni/`; export
  `assets/atomo/molti_elettroni.json`.
- **Costo dichiarato nel riassunto del piano:** circa 30–40 ore di PC in 5–6 sessioni.
- **Decisione per Valerio:** il campo medio con scambio (in pratica Hartree–Fock ricavato dalla
  regola) serve solo per i nodi, ma MOLECOLE.md diceva "niente Hartree–Fock".
- **(osservato) Attenzione:** gli script di sviluppo che il piano dice di copiare in
  `lab/multi_e_dev/` stanno nella cartella temporanea del portatile. **Non sono nel repo** e in
  cloud non esistono.

### 3.2 Nuclei da nodi (cantiere Colla)

**NODO_3D.md**
- **FAIL:**
  - A v=0,6 e v=0,8 il nodo frena e si gonfia (attrito di Peierls–Nabarro). A v=0,8 fallisce con
    ogni stencil e passo h.
  - Nodo + antinodo fermi non si annichilano.
  - Λ sul reticolo −23 % (la coda 1/R è tagliata dalla scatola); ancora −4,6 % con la correzione
    della coda.
  - μp/μn: +6 % sul reticolo, −4,3 % in S3.
  - g_A −52 % e F_π −30 % (attesi).
  - Disegno in negativo: 1,19 % di triangoli contro un limite dell'1 %.
  - Fit spettrale c = 0,37 / 0,46: lo stimatore è sbagliato.
- **Debiti:**
  - Niente elettroni, carica o fotoni nel reticolo.
  - c delle onde di pione = c della luce è un'assunzione.
  - Spin ½ citato dalla letteratura, non mostrato.
  - Rotore rigido semiclassico.
  - m_π = 0.
  - Servono ≥ 4 palline per raggio.
  - L-BFGS invece del flusso di Newton.
- **Rifare:** `bash lab/knot_run.sh knot_kN.py ARGS` e `bash lab/knot2_run.sh SCRIPT ARGS` (codice
  d'uscita 3 = checkpoint, si rilancia).

**NUCLEI_DA_NODI.md**
- **PASS:**
  - N1: forme ed energie.
  - N2: due nodi diventano il toro (uguale a N1 a 10⁻⁸); quattro nodi diventano il cubo.
- **FAIL:**
  - Legame ×11–32 (noto da decenni).
  - L'asse previsto del toro era sbagliato (è venuto lungo z, non y).
- **Debiti:**
  - Solo forme classiche, niente spin/isospin.
  - Il deutone vero è largo (~2 fm), il toro è compatto.
  - Pioni senza massa.
  - N2 è una discesa in energia, non una collisione nel tempo.
  - La simmetria di N1 viene dal punto di partenza.
  - `tet4` fermato a 126 iterazioni.
- **Rifare:**
  - `bash lab/nuclei_drive.sh JOB...` (oppure `bash lab/nuclei_drive2.sh NMAX JOB...`);
  - `python nuclei_run.py JOB --finalize`;
  - `nuclei_summary.py`, `nuclei_figure.py`.
  - Pre-registrazione: commit cd859e2.

**ZOOM_NODO.md** (validato)
- **FAIL:**
  - Z1b: l'errore nel nucleo cresce allargando la zona fine (domina la parete a 4,8).
  - Z1c: E −1,01 %.
  - Motore v1: "gradino" alla cucitura (corretto dopo la pre-registrazione; corse v1 in
    `lab/out/zoomknot_v1bug/`).
  - Al livello h=0,2 imparare non aiuta.
  - Z4 (nodo in moto) saltato.
- **Regola imparata:** punto fisso da h=0,8 in su, k2 ≈ 0,941 e k4 ≈ 0,805; 14 numeri imparati.
- **Debiti:**
  - Termine cinetico degli agglomerati non rinormalizzato.
  - La dinamica su griglie annidate richiede un passo temporale locale.
- **Rifare:**
  - `bash lab/zoomknot_drive.sh "CONFIG MODE" ...`;
  - `zoomknot_learn.py`, `zoomknot_summary.py`, `zoomknot_figure.py`.
  - Pre-registrazione: commit bab46be.

**COLLA.md** (cantiere aperto, documento non finito)
- **Fatto:**
  - C1a: rifit AN e = 4,8426, F_π = 108,26 MeV (PASS).
  - C1c: forme mantenute (PASS).
  - C1d: la sola massa del pione non basta (deutone, ³He, ⁴He: 68 / 170 / 300 MeV).
  - C2a: la quantizzazione rigida peggiora il legame (160 / 344 / 610 MeV). Deutone 160 MeV contro
    Braaten–Carson 158 (PASS).
  - C2b: punto zero di B=2 81,7 MeV; deutone 78 MeV (35× il reale).
  - C3b: il termine sestico alza il legame (previsione sbagliata).
- **FAIL:**
  - C1b: −1,2 % (soglia 1 %).
  - C2a-i: Λ −5,6 % (soglia 3 %).
  - C3b-i.
- **Da fare:**
  - Scansione C3a (GHS, α = 0, 0,5, 0,8, 0,95). **(osservato)** JSON ancora mancanti alla stesura:
    `ghs95_B3`, `ghs80_B3`, `ghs95_B2pair`, `ghs95_B4tet` (`ghs50_B4` è arrivato).
  - Esami C3a-ii (monotonia del legame in α), C3a-iii (forme a α=0,95), C3a-iv (fit di α* su
    ⁴He = 28,3 MeV, poi prevedere d e ³He entro un fattore 2).
  - Tabella C4.
  - Riempire le sezioni segnaposto (`<!-- TABELLONE -->`, `<!-- C3a -->`, `<!-- C4 -->`,
    `<!-- DEBITI -->`).
  - Rilanciare `colla_summary.py` e `colla_figure.py`.
  - Commit con `git add -f` per `lab/out`.
- **Debiti:**
  - C2b è una mezza stima: manca il punto zero del continuo (Casimir).
  - Il deutone a 10 dimensioni di Leese–Manton–Schroers non è stato rifatto (in letteratura
    6,2 MeV, 2,18 fm).
  - Il segno di w per ³He è stato corretto dopo aver visto i numeri (dichiarato).
  - C3a e C3b aggiungono 1 parametro libero ciascuno.
  - `.pt` e `.npz` non sono in git.
- **Rifare:**
  - `bash lab/colla_drive.sh JOB...`;
  - `python colla_modes.py an_B1 12 --wall 540`;
  - `colla_quant.py`;
  - `colla_summary.py` → `lab/out/colla_results.json`.
  - **(osservato)** `colla_drive_a.sh`, `colla_drive_b.sh` e `colla_drive_bc.sh` aspettano i
    processi con `tasklist` di Windows: su Linux vanno adattati.
- **Da chiarire:** PIANO §4 mette come obiettivo "deutone quantizzato ~6 MeV secondo
  Braaten–Carson 1988". COLLA C2a-iv invece confronta con Braaten–Carson 158 MeV (rigido), e
  i 6,2 MeV sono di Leese–Manton–Schroers 1995.

**PRESSIONE.md** (non sulla strada)
- **FAIL:**
  - P2b2 e P2c: previsti.
  - P4: nucleo del reticolo sbagliato del 57 % a h=0,2. Il problema è stato poi risolto da ZOOM_NODO
    (p(0) −0,08 %).
- **Debiti:** il M² del tripolo JLab è citato a memoria; massa classica; m_π = 0.
- **Rifare:** `python lab/press_ode.py`.

### 3.3 Vuoto e coppie (base del passo 6 di PIANO: γ → e⁺e⁻ vicino a un nucleo pesante)

**IDEE_RIPRESE.md**
- **R1 (PASS):**
  - I doppioni danno P_π/P_0 = 1,000000.
  - Contando su un periodo di Bloch intero si ottiene 2× Schwinger (previsto).
- **R2 FAIL:** il flusso di Ollivier-Ricci ha il punto piatto instabile per qualunque passo; il freno
  di flusso non salva la forma.
- **R3 fragile:**
  - I risultati dipendono dal numero di thread (modi degeneri di ARPACK).
  - Il controllo sul reticolo fine fallisce.
- **Debiti:**
  - Unruh non fatto.
  - Nessuna retroazione del campo: le coppie non schermano E.
  - La formula di Sauter era scritta capovolta all'inizio; l'ha trovato il controllo ODE.
- **Rifare:** `lab/redo_r1.py` (~20 s), `lab/redo_r2.py 0|1|2|4|diag`, `lab/redo_r3.py main|fine`.

**COPPIE_IN_SCATOLA.md**
- **Emendamento dichiarato:** pareti a specchio sostituite da pareti di massa, perché gli specchi
  portano uno stato di bordo (errore 10–35 %).
- **FAIL:**
  - A1a: le scatole strette (d m < 1,5) sono contate male (fondo ~2e-4).
  - A1b: la formula senza spin è sbagliata per Dirac fino a un fattore 2.
  - A2: a d=128.
  - A4(a): errore nella pre-registrazione.
- **Rifare:** `lab/box_run.py`, `lab/box_threshold.py D`, `lab/box_formation.py`.

### 3.4 Piano, visione, ponte

- **PIANO.md:**
  - Fonte delle regole e dei 4 cantieri.
  - Obiettivi pre-registrati per i cantieri:
    - agglomeratore: 1s entro lo 0,1 %, massa dal flip entro l'1 %, −k⁴/8m³ entro il 10 %,
      costo ridotto ≥ b²;
    - reticolo: c su ogni livello (0,5 %), passaggio ≥ 99 % con ≥ 16 celle/λ, gamma al livello 0
      e radio al livello 20 nello stesso mondo.
  - **(osservato)** La tabella di stato non ha l'accoppiato.
- **VISIONE_SISTEMA.md:** la proposta del 3 ott (integrare QuTiP, Geant4, PySCF, LAMMPS…) è in
  tensione con la regola 1 di PIANO (niente codici esterni). Vale PIANO.
- **PONTE.md:**
  - Azioni ammesse: `ping`, `build`, `test` (`filtro`, `ignored`), `foto` (argomenti controllati).
  - Richieste in `ponte/richieste/<id>.json`, risultati in `ponte/risultati/<id>/` (`esito.json`,
    `log.txt`, `foto.png`).
  - Un id già eseguito non si ripete.
  - **(osservato)** Il doc parla di tempi misurati sulla Intel UHD, ma il ping elenca anche una
    GTX 1650 Ti. Il commit d782a44 chiede la GPU ad alte prestazioni.

### 3.5 Storico: fotone a due ruoli (non si riapre senza Valerio)

- **REGOLA_FOTONE.md:**
  - La moneta di flusso ai bordi è una correzione, non viene da primi principi.
  - I due "primi passi dietro un oggetto" non sono stati unificati.
  - Rifare: `cargo test --bin fotone -- --test-threads=1` (329 s).
- **REGOLA_BASE.md:**
  - v1 Arago FAIL, quindi v2 con zona di bordo (correzione del modello di bordo).
  - Il contrasto di Arago 3,3 vs 4,1 non è spiegato.
  - Rifare: `bash lab/base_run_all.sh`.
- **ESPERIMENTI_FOTONE.md:** stato fermo al commit 83194fa; non è aggiornato con REGOLA_BASE e
  MULTIFOTONE.
- **MULTIFOTONE.md:**
  - La non-località richiede un "adesso" comune.
  - Rifare: `timeout 570 python lab/multi_<x>.py`. La cache di ~95 MB non è versionata.
- **FOTONE_UI.md:** il tempo reversibile (keyframe + piccole operazioni) è uno schema riusabile.
  La foschia termica è solo visiva.
- **GPU_PORT.md:**
  - sin/cos nativi di WGSL falliscono il test; si usa un sincos polinomiale.
  - Il port non è collegato a Sim né al render.
- **IMMAGINI.md:**
  - R5b FAIL.
  - L'atomo è un segnaposto (idrogeno di Schrödinger esatto).
  - Rifare: `PYTHONIOENCODING=utf-8 timeout 590 python lab/img_<x>.py`.
- **REGOLA_PALLINE_IDENTICHE.md:** l'ingrediente minimo mancante (pesi μ e un operatore locale che
  conservi la probabilità su punti irregolari) non è stato trovato.

### 3.6 Storico: Dirac a due ruoli, DMC, molecole

- **SURFISTA_DIRAC.md:**
  - FAIL: Z1d (−5,10 % contro tolleranza 5 %) e G2 (−13 %).
  - La regola del surfista (v = c ψ†αψ/ψ†ψ) è il passo 3 della zig-zag.
- **COMPTON.md:** le sezioni d'urto di Klein–Nishina sono un INGRESSO. La mappa rango→angolo è una
  scelta, non la traiettoria di Bohm a 6-D.
- **MOLECOLE.md:**
  - Acqua non fatta: problema del segno, nodi fissi.
  - §7.2 e §9 dicono ancora "causa non trovata" (superato da H2_VIBRAZIONE).
  - Rifare: `cargo run --bin molecole -- --misura <target>` e `sh lab/molecole_misura.sh`.
- **H2_VIBRAZIONE.md:**
  - Le barre d'errore Rust erano sottostimate di ~2×.
  - `Curva::fit` in `src/molecole/misura.rs` usa ancora 5 esponenziali (va portato a 6).

### 3.7 Storico: materia dalla luce

- **LUCE_INTRAPPOLATA.md e RICERCA_REGOLE.md:**
  - La luce da sola non si lega. Teorema: nessun campo statico localizzato per un'elettrodinamica
    non lineare convessa.
  - Quando una regola locale fa materia stabile, la materia è fatta del campo di "struttura", non
    di luce.
  - Nessuna previsione riuscita.
- **GLUEBALL.md:**
  - I dati grezzi `lab/out/glue_raw/` (~420 MB) non sono nel repo.
  - Rifare: `lab/glue_run.py BETA [budget]` (exit 3 = altri pezzi da fare). Semi diversi per ogni
    pezzo da 9 minuti.

### 3.8 Chiusi: dimensioni nascoste, forme, Born sul reticolo

- **CARICA_NASCOSTA.md:**
  - Cambiare la velocità c_w sul cerchio non aiuta: c_w è il raggio travestito.
  - Con una rigidità K propria, α e le masse diventano numeri messi a mano.
- **FORMA_NASCOSTA.md:** 18 previsioni a zero parametri tutte sbagliate. Una "copia fantasma"
  (doppione del walk) accoppia ogni livello.
- **FORMA_DALLE_PALLINE.md:**
  - La regola zig-zag raddoppia lo spazio: due fogli intrecciati.
  - Il bordo aperto salda i fogli.
- **EQUILIBRIO_QUANTISTICO.md:**
  - Nel continuo il rilassamento a Born funziona con i vortici.
  - Le corse sono durate 4π invece degli 8π pre-registrati.
- **PRE_ESPERIMENTI.md, GEOMETRIA_AUTOCOSTRUITA.md, IDEE_RIPRESE R3:** sul reticolo il salto casuale
  dei surfisti crea mescolamento finto. Non è stato rimosso, né a θ = 0,125 né a 0,0625.
- **ORBIFOLD.md:**
  - Il quoziente è imposto a mano.
  - La lisciatura delle punte è presa dai libri.
- **OCCHIO_E_FLUSSO.md:**
  - La regola v2 per le punte è stata scritta dopo la calibrazione (emendamento dichiarato).
  - Il flusso di Ollivier-Ricci spezza il grafo in isole.

### 3.9 Altre idee chiuse

- **SINCRONIA.md:**
  - L'incertezza del protocollo (~0,06) è quasi grande quanto l'effetto cercato (0,073).
  - I vortici dei metronomi restano inchiodati al reticolo.
  - Idea per rifarlo (PIANO §6): un campo con ampiezza (Ginzburg–Landau).
- **FLUIDO_IN_NEGATIVO.md:**
  - In un fluido viscoso 2-D la massa del buco dipende dalla scala di tempo (paradosso di Stokes).
  - Il consiglio del documento (difetti topologici in un mezzo 3-D) è stato seguito con i nodi di
    Skyrme.

---

## 4. Debiti aperti (lista unica)

**Regola zig-zag**
1. Ingressi messi a mano: c, ħ, spin, m_e (θ), α, Z, massa del protone (massa ridotta). La regola è
   il Dirac quantum walk noto, non derivata.
2. QED assente: Lamb, g−2, i 25 MHz della struttura fine.
3. Raddoppio dello spazio (due fogli): U(k+(π,π)) = U(k). La correzione di Wilson esiste solo in
   spazio k (non è un automa locale vero).
4. Dispersione trasversale ka/tan(ka) e rallentamento obliquo. Servono ≥ 24–32 celle/λ, oppure una
   correzione dentro la zig-zag che non è stata trovata.
5. Errore di contatto dei livelli s (θ² ln θ). L'1S–2S non vede ancora il Lamb.
6. Idrogeno vero in 3-D non fatto (~3e11 celle). La regola 3-D a Z=60 dà −11,6 %.
7. Surfisti:
   - il lettore va migliorato (equivarianza 1,26e-2);
   - non esiste una guida deterministica che conservi |ψ|² esattamente sul reticolo;
   - il salto casuale crea mescolamento finto.

**Agglomeratore**
8. Spin-orbita al livello grosso: imparare l'accoppiamento della cella con un campo elettrico.
9. 2s all'1 % (a0/16 su 384³, o con lo zoom). Potenziale grosso scritto a mano. Rinculo del protone.
10. Molti elettroni: c'è solo il piano (MOLTI_ELETTRONI_PIANO.md); gli script di sviluppo non sono
    nel repo.
11. Livelli più alti (a0/4, a0/2, …) per molecole e gas. Porto in Rust del livello B = 256.
12. Problemi generali emersi in ZOOM e ZOOM_NODO:
    - memoria per le onde corte;
    - coda a lungo raggio;
    - iperparametri a mano;
    - termine cinetico non rinormalizzato;
    - dinamica su griglie annidate (passo locale).

**Reticolo dell'universo**
13. Porto in Rust, 10–100 volte più veloce.
14. Cucitura conservativa del flusso; restrizione conservativa (contro il fantasma del gamma).
15. Mondo 3-D completo con camino; direzione obliqua fuori dal piano.
16. Palline massive (θ raddoppia a ogni livello) e aggancio agli agglomerati. Radio al livello 20.

**Laboratorio Atomo e accoppiato**
17. Emissione spontanea: serve la luce quantizzata (oggi è semiclassico).
18. Accoppiamento K = 10.
19. Nd59+ invece di H.
20. 1s 3-D sporco.
21. Rifare il confronto E1 con la riga 1s→2p3/2 e un 1s misurato meglio (annotato, non rigiudicato).
22. Cantiere orbitali (PIANO §4):
    - selettore n = 1..4, l = s..f, m e sovrapposizioni;
    - gusci morbidi da lontano, dinamica pallina per pallina da vicino;
    - niente orbitali da tabella.
23. Idrogeno vero a scala di Bohr nel laboratorio (`agglomerato_zz.json` + `reticolo_spec.json`),
    con zoom fino al tremolio di Dirac.
24. Doppia fenditura fatta di atomi agglomerati: confronto col riferimento continuo della stessa
    geometria, con ≥ 24 celle/λ.
25. Prestazioni: ≥ 24 celle/λ a ≥ 30 tick/s non entra nel portatile. Preset di qualità pensati per
    una GPU da 16 GB.

**Nuclei**
26. Legame ×11–32. Finire Colla (C3a, C4, tabellone, debiti, figura, commit).
27. Quantizzazione solo rigida; punto zero del continuo; deutone LMS non rifatto.
28. Effetto di m_π sulla pressione e sul D-term. M² di JLab da verificare sulla fonte primaria.
29. Nodo a v = 0,8 (h ≈ 0,12 o stencil di ordine 6). Misurare direttamente il freno di
    Peierls–Nabarro. Coda di Λ.
30. Da chiarire la cifra "~6 MeV secondo Braaten–Carson" in PIANO §4 (vedi §3.2 di questo file, voce COLLA.md).

**Vuoto e coppie**
31. Unruh; retroazione del campo; conteggio nelle scatole strette (d m < 1,5).
32. Passo 6 di PIANO: γ → e⁺e⁻ vicino a un nucleo pesante (vuoto di coppie in 2-D HD).

**Documenti e infrastruttura**
33. PIANO.md: aggiungere l'accoppiato alla tabella di stato.
34. COLLA.md: completare e committare.
35. MOLECOLE.md: aggiornare §7.2 e §9 con la causa trovata.
36. `Curva::fit` (5 esponenziali) e bootstrap in `src/molecole/misura.rs`.
37. ESPERIMENTI_FOTONE.md: stato fermo al commit 83194fa.
38. PONTE.md: parla di Intel UHD, ma c'è anche la GTX 1650 Ti.
39. `colla_drive_a/b/bc.sh` dipendono da `tasklist` di Windows.
40. File grandi da rigenerare (vedi §2.5).
41. Idee di PIANO §6, non sulla strada:
    - geone (serve la gravità);
    - sincronia con Ginzburg–Landau;
    - glueball eccitato e β = 12;
    - problema delle masse (muone 207, tau 3477);
    - energia del vuoto (10⁴⁴–10¹²⁰);
    - anomalie (vita del neutrone 878 vs 888 s, Hubble);
    - la mela (ossidazione del ferro).

---

## 5. Strumenti riutilizzabili

| strumento | dove | a cosa serve |
|-----------|------|--------------|
| Regola zig-zag di riferimento | `lab/zz_core.py` (+ `zz_fast.py` numba); Rust `src/zigzag/walk.rs` (anche `Radial`) | regola congelata, con autotest; base di quasi tutti gli studi Python sulla strada |
| Agglomeratore sulla zig-zag | `lab/agg_zz_core.py` (`auto_level`, `kcheck` con `energy_mode="state"`) | sceglie da solo blocco b, stencil w e numeri per cella; impara la regola grossa per minimi quadrati; controlla. Livelli già imparati in `assets/atomo/agglomerato_zz.json` e `lab/out/agg_zz_3d_{A,B}.pkl` |
| Algoritmo generico dello zoom | `lab/zoom_core.py` | blocca→impara→passo→controlla→ripeti (usato anche da `zoomknot_`) |
| Zoom annidato sul nodo | `lab/zoomknot_core.py` | 8 livelli 2×2×2, banda di cucitura con interpolazione cubica, k2 e k4 imparati per livello |
| Motore dei nodi | `lab/knot_core.py`, `lab/knot2_core.py` (stencil `old`/`c2`/`c4`) | Skyrme 3-D sul reticolo. c4: errore di fase 0,17 % e di gruppo 0,9 % per kh ≤ π/4 |
| Nuclei | `lab/nuclei_core.py` | mappe razionali (B=1..4), misure di forma |
| Colla | `lab/colla_core.py`, `colla_modes.py` | energia generalizzata, metrica cinetica, tensori d'inerzia, Hessiana; modi con LOBPCG + precondizionatore FFT |
| Reticolo sparso multilivello | `lab/reticolo_core.py` + formato `assets/atomo/reticolo_spec.json` | livelli h_k = 2^k h_0, blocchi 16² (8³ in 3-D), camino, passo Berger–Oliger, cucitura v2 |
| Occhio (lettore di forma) | `lab/eye_core.py` (anche `big_core.py`, `geo_core.py`) | informazione mutua → forma. Soglia buchi: persistenza/diametro > 0,4. Cavità H2: soglia 0,237. Punte: regola v2. Testare sempre sui nulli prima |
| Riferimenti ottici esatti | `lab/optics.py` | Rayleigh–Sommerfeld esatto, punteggi (usato anche da ZIGZAG §3.4) |
| Walk con pesi sui link + surfisti | `lab/pre_core.py`, `lab/relax_core.py` | surfisti tipo Bell (manopola λ), rilassamento, H-bar |
| Cerchio nascosto | `lab/kk_core.py` | zig-zag con una dimensione compatta |
| SU(2) su reticolo | `lab/glue_core.py` | heat-bath + over-relaxation, correlatori, GEVP |
| Camera e PSF | `lab/img_core.py` | modello di sensore, PSF, Planck |
| Accoppiamento luce-atomo | `src/atomo/accoppiato.rs` | Maxwell–Dirac semiclassico sulla zig-zag |
| Tempo reversibile | `src/photon/timeline.rs` | keyframe + piccole operazioni per fotogramma (non ancora usato in `src/atomo`) |
| sincos polinomiale GPU | `src/photon/gpu.rs`, `assets/shaders/photon_grains.wgsl` | evita l'errore dei sin/cos nativi di WGSL; accumulatore a punto fisso i64 |
| **Pre-registrazione** | `lab/<studio>_prereg.py` → `lab/out/<studio>_prereg.json` | esami, soglie, configurazioni committati prima delle corse; emendamenti dichiarati in un file a parte (es. `box_amend.json`, `eye_amend.json`); i cambi del motore dopo la pre-registrazione si dichiarano e si tengono i log vecchi (es. `reticolo_v1/`, `zoomknot_v1bug/`) |
| **Corse a pezzi** | `lab/knot_run.sh`, `knot2_run.sh`, `nuclei_drive.sh`, `zoomknot_drive.sh`, `colla_drive.sh`, `glue_run.py` | ogni script ≤ 10 min (`timeout 590`/`570`), checkpoint, exit code 3 = "rilanciami" |
| **Ponte** | `ponte/ponte.py`, `docs/PONTE.md` | il cloud scrive `ponte/richieste/<id>.json` e fa push; il portatile esegue `ping`/`build`/`test`/`foto`; risultati in `ponte/risultati/<id>/` dopo il pull |

---

## 6. Lezioni apprese

1. **La risoluzione del reticolo è fisica.**
   - La luce trasversale va col fattore ka/tan(ka): 0,785 (8 celle/λ), 0,907 (12), 0,948 (16),
     0,977 (24), 0,987 (32).
   - Per errori di pochi % servono ≥ 24 celle per λ; per la luce obliqua tra livelli ≥ 32.
   - Le tolleranze si scrivono partendo da come scala l'errore. La L1a di ZIGZAG è fallita per una
     soglia del 2 % su un errore (a/λ)².
2. **Almeno 4 palline per raggio di un difetto.** Sotto ~3 si rompe perfino la topologia. Il nodo in
   moto regge con 3,4–4,0 palline per raggio contratto.
3. **"Il forte diventa più forte" fa isole, non strade.** È successo:
   - coi link rinforzati (GEOMETRIA_AUTOCOSTRUITA);
   - col flusso di Ollivier-Ricci (OCCHIO_E_FLUSSO);
   - col freno di flusso (IDEE_RIPRESE R2).
   Lo stato piatto era instabile. Prima di chiedere emergenza serve una dinamica con lo stato piatto
   stabile.
4. **I surfisti sono un'interpretazione.**
   - Stessa statistica, storie diverse: la mappa rango→angolo di COMPTON, l'"adesso" comune di
     MULTIFOTONE.
   - Sul reticolo il salto casuale crea mescolamento finto (PRE P1c, GEOMETRIA R2).
   - Il surfista di Bohm–Dirac gira con raggio uguale alla larghezza del pacchetto, non ħ/2mc.
5. **Un fit non è una previsione.** Contare le manopole:
   - Koide passa 1 volta su 5;
   - 2 parametri = fit (FORMA_NASCOSTA);
   - "3 famiglie" in 17 scelte su 72 (ORBIFOLD);
   - probabilità di trovarlo per caso altrove 21 % (RICERCA_REGOLE).
6. **Spesso il colpevole è il fit o lo stimatore, non le palline.** Esempi:
   - vibrazione di H2: il fit a 5 esponenziali;
   - stimatore a parabola di IMMAGINI;
   - fit spettrale di NODO_3D S4;
   - rimbalzi sul bordo della scatola in `zz_listen` (risolti con la spugna).
   Le barre d'errore bootstrap erano ottimistiche di ~2×.
7. **Un lettore di forma vede forme nel rumore** finché non viene provato sui nulli (P2: 60–100 % di
   falsi positivi). Poi il lettore riparato non vede più le cose piccole.
8. **Continuo ≠ reticolo.** Esempi:
   - Born per mescolamento funziona nel continuo, sul reticolo è fragile;
   - la zig-zag raddoppia lo spazio, e i conteggi su un periodo di Bloch intero valgono 2×;
   - le pareti a specchio portano stati di bordo.
9. **La cucitura fine/grosso è il punto delicato** degli agglomerati: gradino v1 in ZOOM_NODO,
   cucitura v2 nel RETICOLO. Nelle zone a legge di potenza la regola imparata va a un punto fisso.
10. **Una regola grossa senza memoria sbaglia sulle onde corte** (oltre ~0,4 del Nyquist). La regola
    grossa è a lungo raggio.
11. **Confrontare col riferimento giusto:**
    - la riga dell'atomo 3-D è 1s→2p3/2 della regola 3-D, non ω21 della forma radiale;
    - la doppia fenditura va confrontata col continuo della stessa geometria, non con la formula a
      due punti.
12. **Discesa ≠ dinamica, e la simmetria può venire dal punto di partenza** (NUCLEI N1). La prova vera
    è partire da pezzi separati (N2).
13. **Riproducibilità.**
    - I modi degeneri di ARPACK cambiano col numero di thread: fissare il vettore iniziale e
      scegliere i modi in modo canonico.
    - Semi diversi per pezzo cambiano i bit se la corsa viene spezzata diversamente.
14. **Il tempo di calcolo è un muro reale.** Esempi: radio al livello 20 ~15 ore in Python; v=0,8 del
    nodo non fatto. Script ≤ 10 min con checkpoint; il lavoro pesante va in Rust o via ponte.
15. **La luce da sola non fa materia.** Con il vuoto di Euler–Heisenberg mancano 10⁴–10⁵, e dove
    sarebbe abbastanza forte crea coppie. La materia stabile richiede topologia più un termine
    anti-collasso (Skyrme).
16. **Semiclassico ≠ quantistico.** Uno stato eccitato puro non decade senza luce quantizzata.
17. **Git.**
    - `git add` con percorsi espliciti; `git add -f` per figure e JSON di `lab/out`.
    - Mai `.npz`/`.pt`/`.npy`/`.pkl` grandi.
    - I log tracciati che vengono scritti di continuo rompono l'autostash di `git pull --rebase`.
    - Mai force-push.
