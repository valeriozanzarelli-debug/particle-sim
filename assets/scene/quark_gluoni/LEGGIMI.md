# Scene di quark e gluoni per il programma (dati calcolati, studio `qg_`)

Questa cartella contiene i dati **calcolati** dallo studio `docs/QUARK_GLUONI.md`, pronti per essere caricati da una
scena del programma unico e disegnati con i bagliori (`src/atomo/bagliori.rs`: un bagliore per pallina, coda corta con
`scia`). Niente è disegnato a mano: le palline sono campioni del campo misurato (quark) o pezzetti di energia del campo
portati dalla sua stessa corrente (gluoni). Scritti da `lab/qg_esporta.py`, ricontrollati da `lab/qg_export_check.py`
(esito in `lab/out/qg_export_check.json`). Formato e disposizione dei byte sono **gli stessi** degli esportatori di
riferimento del thread delle scene (`lab/sc_esporta_quark.py`, `lab/sc_esporta_gluoni.py`, ramo
`claude/project-thread-sxhn28`: `particle-sim/scena-quark-v1` e `particle-sim/scena-gluoni-v1`). Due differenze
dichiarate: nei quark `quota` e `livello` sono i valori esatti del campionamento (non ricalcolati alla posizione
spostata a caso); nei gluoni si tiene **un fotogramma salvato su 3** (più l'ultimo) per stare sotto ~25 MB.

Tutti i numeri sono **little endian**. Ogni cartella ha un `manifest.json` che dice tutto (forma, colonne, tempi,
voti, commit sorgente). Niente `.npz/.npy`: solo byte grezzi.

| cartella | file | cosa | byte |
|---|---|---|---|
| `quark_Y/` | `palline.f32` | 60 000 palline del tubo a Y | 1 200 000 |
| `quark_Y/` | `campo.f32` | il campo misurato sul reticolo 12³ | 27 648 |
| `gluoni_video_grande/` | `palline.bin` | film, ampiezza 1 (si urtano), 48 fotogrammi | 7 680 000 |
| `gluoni_video_piccolo/` | `palline.bin` | film, ampiezza 0,05 (si attraversano), 48 fotogrammi | 7 680 000 |
| `gluoni_video_tre/` | `palline.bin` | film, tre gluoni da x, y, z (aggiunta dichiarata), 48 fotogrammi | 7 680 000 |

## 1. Tre quark e il tubo a Y (`quark_Y/`)

**Cos'è.** Un'istantanea in **tempo immaginario**: SU(3) su reticolo 12⁴, β = 5,7 (a ≈ 0,17 fm), 600 misure Monte
Carlo, tre quark **fermi e infinitamente pesanti** ai vertici di un triangolo (d = 4, lato 5,66 passi ≈ 1 fm), loop di
Wilson barionico con T = 4. **Non è un film**: nella scena la nuvola sta ferma (si può far girare la camera o far
"respirare" la luminosità, ma questo è disegno, non fisica). Voti: A placchetta PASS, B tensione di corda PASS,
C forma Y PASS (R_Y = 1,23 ± 0,05).

**`palline.f32`**: `n × 5` float32, riga dopo riga (20 byte a pallina):

| colonna | significato | unità |
|---|---|---|
| x, y, z | posizione | passi di reticolo dal **centro del triangolo** (1 passo = 0,17 fm) |
| quota | quota elettrica dell'azione espulsa Δp_t/(Δp_t+Δp_s), 0…1 | colore (tavolozza zig-zag: ciano = magnetico, viola = metà, oro = elettrico) |
| livello | f/f_max, 0…1 | luminosità 0,1 + 0,9·livello^0,35 (come `zzscena`) |

Le palline sono in **ordine casuale**: accendere le prime k dà un sottoinsieme fedele. Le tre posizioni dei quark sono
in `manifest.json` → `quark` (tre bagliori bianchi grandi). Nasce con densità ∝ f (spline cubica dal reticolo), solo
dove f supera 2 errori jackknife e tolto il fondo uniforme (mediana di f a più di 5 passi dal centro): **scelta di
disegno dichiarata**; gli esami usano f com'è.

**`campo.f32`**: `4 × 12 × 12 × 12` float32, ordine `[canale][x][y][z]` (z il più veloce). Canali: `f`, `f_errore`
(jackknife), `dpt`, `dps` (parte temporale e spaziale; f = (dpt+dps)/2). Indice = sito del reticolo periodico; il punto
dove si toccano le linee del loop è il sito (0,0,0), i quark sono nei siti (4,0,0), (0,4,0), (0,0,4); il centro del
triangolo è (4/3, 4/3, 4/3). Nessuna sottrazione del fondo. Unità: reticolo (densità d'azione adimensionale).

## 2. Gluoni (`gluoni_<corsa>/`)

**Cos'è.** Campo dei gluoni **classico** (Yang–Mills SU(3) in tempo reale, salto della rana, dt = 0,1), scatola
periodica (`manifest.forma`, in passi di reticolo). Ogni pallina è un pezzetto fisso dell'energia: nasce con densità ∝
energia e viene **portata dal vettore di Poynting** del campo (v = S/ε) dentro la simulazione. **Non è quantistico,
non c'è confinamento**, i colori iniziali dipendono dalla gauge (il fatto che si urtino no).

**`palline.bin`**: fotogramma dopo fotogramma, ognuno di `8·n` byte (n = `manifest.n` = 20 000):

```
q   n × 3 uint16   posizione / forma * 65535, intrecciati: x y z x y z ...   (6n byte)
rb  n uint8        quota elettrica E²/(E²+B²) * 255                         (n byte)  -> colore
lv  n uint8        255 * min(ε/ε_max, 1)^0,35                               (n byte)  -> luminosità 0,1 + 0,9 lv/255
```

Sono **le stesse palline nello stesso ordine in ogni fotogramma** (le prime 20 000 delle 160 000/200 000 della corsa,
nate in ordine casuale), quindi le code (`scia`) sono vere traiettorie. ε_max = 99,5° percentile di ε alle palline nel
fotogramma 0 (`manifest.eps_max`). Tempi in `manifest.t` (unità di reticolo, c = 1); fra due fotogrammi tenuti passano
1,2 unità (3 salvataggi da 0,4): la scena interpoli linearmente, tenendo conto che la scatola è **periodica** (una
pallina che esce da un lato rientra dall'altro: non interpolare attraverso il bordo). Fedeltà G7 (correlazione palline
del film intero / energia del campo per cella 8³) per fotogramma in `manifest.G7_fedelta`; S8 = quanto le 20 000 tenute
rappresentano tutte le palline.

## 3. Come leggere

Python:
```python
import json, numpy as np
d = "assets/scene/quark_gluoni/"
q = np.fromfile(d + "quark_Y/palline.f32", "<f4").reshape(-1, 5)        # x y z quota livello
m = json.load(open(d + "gluoni_video_tre/manifest.json")); n = m["n"]
b = np.fromfile(d + "gluoni_video_tre/palline.bin", "u1").reshape(m["n_fotogrammi"], 8 * n)
k = 10                                                                    # fotogramma
pos = b[k, :6 * n].copy().view("<u2").reshape(n, 3) / 65535 * np.array(m["forma"])
quota, lv = b[k, 6 * n:7 * n] / 255, b[k, 7 * n:] / 255
```

Rust:
```rust
let raw = std::fs::read("assets/scene/quark_gluoni/quark_Y/palline.f32")?;
let f: Vec<f32> = raw.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect();
for p in f.chunks_exact(5) { /* p[0..3] posizione, p[3] quota, p[4] livello */ }
let g = std::fs::read("assets/scene/quark_gluoni/gluoni_video_tre/palline.bin")?;
let (n, forma) = (20000usize, [64.0f32, 64.0, 64.0]);          // dal manifest
let fr = &g[k * 8 * n..(k + 1) * 8 * n];
let u = |i: usize| u16::from_le_bytes([fr[2 * i], fr[2 * i + 1]]) as f32 / 65535.0;
for i in 0..n { let pos = [u(3*i) * forma[0], u(3*i+1) * forma[1], u(3*i+2) * forma[2]];
                let (quota, lv) = (fr[6*n + i] as f32 / 255.0, fr[7*n + i] as f32 / 255.0); }
```

## 4. Misurato o scelto

* **Misurato** (esce dalla simulazione con β o con le condizioni iniziali): f, dpt, dps e i loro errori; posizioni,
  energia e quota elettrica delle palline dei gluoni; i tempi; G7.
* **Scelto** (disegno, dichiarato): quante palline (60 000 / 20 000), la soglia 2σ e la sottrazione del fondo dei
  quark, la legge della luminosità (^0,35), ε_max, un fotogramma su 3, il seme del campionamento (7).
* **Messo a mano** (debito della fisica): le posizioni dei quark (statici, infinitamente pesanti); le condizioni
  iniziali dei pacchetti di gluoni.

## 5. Stato onesto

* Quark: **tempo immaginario, quark statici, niente quark dinamici (quenched), reticolo grossolano (2 fm di scatola)**.
  La Y è un risultato noto della QCD su reticolo: qui è validazione, non scoperta.
* Gluoni: **classici**. L'urto (20 % dell'energia spostata ad ampiezza 1) è fisica nota delle equazioni di Yang–Mills.
  G7 (palline fedeli al campo, r ≥ 0,9 in ogni fotogramma, votato sulle due corse pre-registrate): **FAIL** di un
  soffio. `video_grande` min r = 0,8998 (un solo fotogramma sotto 0,9, a t = 50,0; mediana 0,974); `video_piccolo` min
  r = 0,9093 (mediana 0,974); `video_tre` (aggiunta, non votata) min r = 0,934. Le palline seguono una corrente
  d'energia sul reticolo non esattamente conservata. Le 20 000 tenute rappresentano tutte le palline con r ≥ 0,993 (S8).
* Le corse sono state rifatte in parte con gli stessi semi il 2/11 e riproducono **bit per bit**
  (`lab/out/qg_repro.json`).
* Controllo di lettura: `python3 lab/qg_export_check.py` → `lab/out/qg_export_check.json` (tutto ok: posizioni u16
  identiche, quota entro 0,5/255, luminosità entro 0,5/255, palline e campo dei quark entro l'arrotondamento float32).
* Il thread delle scene usava la cartella `assets/scene_calcolate/`: questi file sono gli stessi formati, messi qui.
