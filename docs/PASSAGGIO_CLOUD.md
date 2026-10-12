# Passaggio al cloud — stato esatto dei cantieri (13 ottobre 2026, notte)

Questo file dice alla sessione cloud **da dove riprendere** ogni lavoro aperto al momento del
passaggio dal portatile. Il quadro generale è in `docs/PIANO.md`, l'indice degli studi in
`docs/CONTESTO_CLOUD.md`.

## Hardware del portatile (dal ponte, ping del 12/10 22:47 UTC)
- CPU Intel i5 (famiglia 6 modello 165), 8 thread logici, Windows 11, Python 3.14, Rust con cargo.
- GPU: **Intel UHD Graphics** (1 GB condivisa) **e NVIDIA GeForce GTX 1650 Ti (4 GB)**.
- Problema aperto: wgpu vede solo la Intel (Vulkan e DX12), anche chiedendo l'adattatore ad alte
  prestazioni. Valerio deve impostare `Laboratorio Atomo.exe` (e `target/debug/atomo.exe` in
  `../ponte_lab`) su "Prestazioni elevate" in Impostazioni di Windows → Schermo → Grafica, oppure nel
  pannello NVIDIA. Poi verificare con un `foto` via ponte: il nome della scheda compare sulla barra
  di scala.

## 1. Laboratorio Atomo (`src/atomo/**`, `docs/ATOMO.md`)
Commit pulito `d782a44`: build verde, 7 test passati (12 ignorati: lunghi o FAIL registrati).

**Fatto**
- Scena fotone libero 3D con palline che cavalcano il campo, elica di E, surfisti spenti di default.
- Scena orbitali (n = 1..4, s/p/d/f, sovrapposizioni che battono) — livelli: 1s +39,9%, 2p3/2 +3,5%,
  3d5/2 +1,8%, 4f7/2 +1,3% rispetto a Dirac esatto (celle grosse vicino al nucleo).
- Tasto MISURA (M) in tutte le scene, con "collassa" (§3.4 di ATOMO.md).
- Accoppiamento luce–atomo nei due versi (`accoppiato.rs`), semi-classico, accoppiamento ×10
  (debito), emissione solo da sovrapposizione (debito). Pre-registrazione v2: ampiezza della luce
  rimandata vs dipolo −0,6% **PASS**; frequenza +7,1% vs 1s→2p½ **FAIL** (la riga emessa è la
  1s→2p3/2, +3,5%, e l'1s 3-D sul reticolo non è un livello pulito). Decadimento non visibile dal vivo
  (~8900 ħ/mc² anche con ×10).

**Da fare, in ordine**
1. Controllo dal vivo dell'impulso accoppiato (via ponte, azione `foto`):
   `{"azione":"foto","args":{"scena":"atomo","frames":4500,"spara-al":3600,"distanza":300,"pitch":0.5}}`
   (gli argomenti ammessi dal ponte sono in `FOTO_ARGS` di `ponte/ponte.py`: se il laboratorio
   ne guadagna di nuovi, aggiungerli lì con il loro controllo).
2. Facoltativo: E1v3 contro 1s→2p3/2 con ascolto più lungo dell'1s, **nuova** pre-registrazione.
3. Task B: luminosità ∝ ρ, moto lungo j/ρ, scie corte in tutte le scene; sovrapposizioni fino a 4
   stati in atomo e orbitali (`fotdisegno.rs`, `zzscena.rs`, `orbitali.rs`).
4. Task C: strato "lunga esposizione" (|ψ|² o |F|² mediato su un tempo scelto) al posto dei surfisti.
5. GPU: preset Bassa/Media/Alta/Ultra con budget di memoria calcolati; port WGSL della regola (luce,
   elettrone, sorgente accoppiata) con test GPU vs CPU ~1e-6 e coefficienti di rotazione calcolati
   sulla CPU; `--qualita ultra --foto N` → `lab/out/atomo_ultra_*.png`. Media/Alta dimensionati sulla
   GTX 1650 Ti, Ultra per una scheda da 16 GB.

Test di sviluppo ignorati da conoscere: `accoppiato_sviluppo`, `accoppiato_post_hoc`,
`accoppiato_righe`, `accoppiato_livelli_3d`.

## 2. Colla dei nuclei (`lab/colla_*`, `docs/COLLA.md`)
Pre-registrazione in `78d76a7`. Agente fermato al passaggio (13/10, ~00:10). Il documento
`docs/COLLA.md` è **quasi completo**: C1 (massa del pione), C2a (trottola quantistica, deuterone
160 MeV come Braaten–Carson 158), C2b (punto zero, deuterone 78 MeV, ancora ×35), C3a (GHS: legame
giù di ~20 volte da α = 0 a 0,95; α* = 0,894 dal ⁴He è un **fit**), C3b (sestico: il legame **sale**,
FAIL), C4 tabella e raccomandazione per l'agglomeratore, debiti.

**Mancano**
1. La corsa `ghs95_B3` (interrotta a metà, ripartire con `lab/colla_run.py ghs95_B3 --wall 480`
   dal suo checkpoint se esiste, altrimenti da capo) e, se in coda, `ghs0_B3`, `ghs50_B3`.
2. La figura `lab/out/colla_figure.png` (`lab/colla_figure.py`) e il riepilogo
   (`lab/colla_summary.py`), poi aggiornare la tabella C4 di `docs/COLLA.md` con il B = 3.
3. Il debito più importante dal lato "0 parametri": la quantizzazione del deuterone con le vibrazioni
   accoppiate (non solo corpo rigido), vedi la sezione "Debiti" del documento.

## 3. Atomi a più elettroni (nuovo, non iniziato)
Piano giudicato in `docs/MOLTI_ELETTRONI_PIANO.md` (espansione a cluster come spina dorsale). Primo
passo: audit del nucleo di diffusione e permesso del tempo immaginario (`lab/multi_e_kernel.py`),
poi solutore esatto a pochi corpi (`lab/multi_e_ecg.py`), poi elio. Va pre-registrato prima di
qualunque corsa di produzione.

## 4. Già chiusi oggi (solo per riferimento)
Agglomeratore zig-zag (`AGGLOMERATO_ZIGZAG.md`), reticolo universo (`RETICOLO_UNIVERSO.md`), zoom del
nodo, nuclei da nodi, pressione, nodo 3D bis, glueball, sincronia, coppie in scatola, H₂.

## 5. Pulizia git da fare una volta sola
- In `git stash list` ci sono 3–4 voci `autostash`: sono vecchie copie dei log della colla e di
  `lab/reticolo_*.py`, più vecchie dei file committati. **Non applicarle.** Si possono eliminare
  con `git stash clear` dopo aver verificato che `docs/COLLA.md` e i file colla sono committati.
  (Sono solo sul portatile, non su GitHub: il cloud non le vede.)
- Non tracciare log scritti di continuo in `lab/out/` (causa di tutto questo).

## 6. Ponte
`python ponte/ponte.py` sul portatile, in esecuzione continua (Ctrl+C per fermarlo). Lavora nella
copia `../ponte_lab`. Prova riuscita: `ponte/risultati/esempio_ping/`.
