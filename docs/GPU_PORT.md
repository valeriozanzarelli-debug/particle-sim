# Port GPU dei grani del fotone

Stato: modulo autonomo `src/photon/gpu.rs` + shader `assets/shaders/photon_grains.wgsl`,
compilato nel binario `fotone` ma **non ancora collegato** a `Sim` né al rendering.
Obiettivo: fotoni fatti di milioni di grani, con la fisica identica a `physics.rs`.

## 1. Cosa fa la GPU

Riproduce `advance()` di `physics.rs` grano per grano:

- traslazione a velocità c lungo `dir`; fase = `phase0 + k·|pos − origin|`;
- gli elementi ottici del fotone (max 4, ordinati lungo il verso di volo, stessa
  lista di `Sim::emit`) vengono attraversati **in ordine e con incroci di piano
  esatti** anche con passi grandi (`t = (z − pos.z)/dir.z`, `pos.z = z` esatto);
- **Wall** → riemissione di Huygens dalla fenditura: l'angolo usa lo **stesso
  hash per grano** di `physics.rs` (`hash3(seed, id, gi)`, splitmix a 64 bit
  emulato in WGSL con coppie di `u32`), quindi CPU e GPU scelgono lo stesso angolo;
  filtri di polarizzazione sulle fenditure, peso `1/U_inc`, obliquità, fase −π/4;
- **Screen** → il grano muore e deposita l'ampiezza complessa nel suo bin;
- **Polarizer / Wall a livello di fotone** (regola di Born): calcolate sulla GPU
  (vedi pass `probe` e `decide`) con un numero casuale estratto dall'host al
  lancio; nessun viaggio GPU→CPU durante il volo.
- **BellSplitter**: per ora attraversato senza effetto; il collasso entangled
  (`entangle::measure`) resta sulla CPU (servirà una lettura di `status` e un
  `set_jones` per slot — non implementato).

Precisione: le `sin/cos` native WGSL garantiscono solo 2⁻¹¹ e sul driver Intel UHD
l'errore reale è ~1.4·10⁻⁵ sull'angolo di riemissione (→ 1.2·10⁻³ sulle posizioni
70 unità dopo: il test falliva). Lo shader usa quindi una `sincos` polinomiale
(riduzione a quadranti + Taylor fino a r¹⁰, errore ~10⁻⁷) e riduce la fase a
[−π, π] prima della trigonometria.

## 2. Layout dei dati

Tutte le struct sono `#[repr(C)] Pod` in Rust con `assert!` statico sulla
dimensione, e documentate nell'header dello shader.

| Buffer | Tipo | Dimensione | Contenuto |
|---|---|---|---|
| `params` | uniform | 96 B | `dist, k, λ, w0, θmax, screen_half, n_bins, n_grains, slit_c0, slit_c1, slit_half, filter_mask, filter0: vec4, filter1: vec4, amp_scale, self_scale, prob_scale, n_slots` |
| `grains` | storage rw (+VERTEX) | **64 B/grano** | `pos: vec3, amp, dir: vec3, phase0, origin: vec3, state: u32, jones: vec4 (re e1, im e1, re e2, im e2)` |
| `photons` | storage rw | 112 B/slot | `elem_z, elem_kind, elem_arg, rand, op_jones: vec4; n_elems, key_lo, key_hi, first_grain, n_grains, front, op, fate` |
| `status` | storage atomic | 32 B/slot | `reach_mask, alive, den (i64), num (i64), fate, front` |
| `screen` | storage atomic | 48 B/bin/slot | `Re a₁, Im a₁, Re a₂, Im a₂, Σ|a|²` come **i64 a virgola fissa** (lo+hi), `count`, pad |

`state` del grano: bit 0‑7 = prossimo elemento, bit 8 = vivo, bit 9 = diffratto,
bit 16‑31 = slot del fotone. Più fotoni in volo condividono un solo buffer di
grani; ogni fotone occupa un intervallo che **inizia a un multiplo di 64**
(dimensione del workgroup), così un workgroup non contiene mai due fotoni e le
riduzioni in memoria condivisa sono corrette. I grani di riempimento sono morti.

### Deposito sullo schermo: virgola fissa, nessun limite

Ogni deposito è convertito in `i32` (`round(x·scale)`) e sommato con
`atomicAdd` sulla parola bassa, con propagazione del riporto sulla parola alta:
un accumulatore **i64 esatto** costruito con atomici a 32 bit.

- L'ordine dei thread non conta (somma intera associativa): risultato
  deterministico bit a bit.
- Scale scelte dall'host (`set_lab`) come potenze di due tali che il deposito
  massimo possibile, `w_max = L / (cos²θmax·√(λL)) / U_inc(bordo fenditura)`, stia
  in 2³⁰. Per la doppia fenditura di default: `amp_scale = 2²⁵`, `self_scale = 2²⁰`.
  Con 64 bit restano ~2.7·10¹¹ depositi massimi prima dell'overflow: ben oltre i
  milioni di grani.
- Nessun contenitore a capienza fissa: al contrario del vecchio `max_per_cell`,
  qui non si può perdere nulla.

Le "hits" (la y di ogni atterraggio, usata da `finalize` per posizionare il click
in verticale) non sono ancora portate. Proposta: per bin, tenere la y del grano
con hash minimo (`atomicMin` su una chiave `hash<<32|indice`), cioè un campione
uniforme deterministico: 2 parole per bin invece di un vettore che cresce.

## 3. Pass per tick

Un tick = un `tick(dist, n)`; `n` tick consecutivi vanno in **un solo command
buffer**, senza nessuna lettura:

1. `clear_buffer(status)`
2. **`probe`** (1 thread/grano): per i grani vivi calcola quali elementi il
   fronte del fotone raggiunge entro `dist` (maschera di bit) e le somme della
   regola di Born per l'elemento `front`: per il Wall Σ|J|² e Σ attraverso le
   fenditure (con filtri), per il Polarizer Σ|J|² e Σ|⟨a|J⟩|². Riduzione ad
   albero in memoria di workgroup → **1 atomico per workgroup**, non per grano.
3. **`decide`** (1 thread/fotone): come `front_events` di `physics.rs`: per ogni
   elemento raggiunto applica la regola di Born con il numero casuale estratto al
   lancio, scrive `op = KILL | SET_JONES`, `fate`, `front`.
4. **`advance`** (1 thread/grano): applica `op` (uccide tutti i grani o imposta
   la loro Jones), poi traslazione, incroci di piano, Huygens, deposito;
   conta i sopravvissuti (riduzione di workgroup → `status.alive`).

Dispatch 2D quando i workgroup superano 65 535 (> 4.19 M grani).

Lancio di un fotone (`launch`): upload dei grani (64 B ciascuno), del record del
fotone e azzeramento del suo schermo. Letture (piccole, solo su richiesta):
`read_status()` 32 B/fotone, `read_screen(slot)` 48 B × bin (11.5 KB per 240 bin).
I grani **non** vengono mai letti (solo nel test, `read_grains`, `#[cfg(test)]`).

## 4. Rendering senza readback per frame

Il rendering deve leggere direttamente `grain_buffer()` (creato con uso
`STORAGE | VERTEX`: stride 64, `pos` a offset 0 `Float32x3`, `phase0` a 28,
`origin` a 32, `state` a 44, `jones` a 48). Un buffer però vive in **un solo
device**; il vecchio `gpu_pps.rs` usava un device dedicato, e per questo doveva
copiare tutto sulla CPU e ricaricarlo nelle mesh di Bevy a ogni frame.

Piano:

1. **Stesso device di Bevy.** Aggiungere `GrainEngine::from_device(device, queue)`
   (oggi `new()` crea un device dedicato, comodo per i test). Bevy 0.15 usa
   wgpu 23, la stessa versione del crate: `RenderDevice::wgpu_device()` e
   `RenderQueue` (Arc di `wgpu::Queue`) sono direttamente utilizzabili.
2. **Motore nel render world.** Una risorsa del render world possiede il
   `GrainEngine`; il main world estrae ogni frame solo i comandi (dist del frame,
   lanci con le Jones/numeri casuali, lab cambiato) via `ExtractResource`.
   Un sistema in `RenderSet::Prepare` (o un nodo del render graph prima del main
   pass) codifica `probe/decide/advance` nello stesso encoder del frame.
3. **Disegno.** Una pipeline custom (come l'esempio `custom_shader_instancing`
   di Bevy) disegna un quad/punto per istanza con `grain_buffer()` come instance
   buffer o storage letto nel vertex shader (vertex pulling). Il colore (freccia
   di fase, campo E = Re[J·e^{i(φ−ωt)}], grani morti nascosti tramite `state`) si
   calcola nel vertex shader con gli stessi dati: `render.rs` oggi fa lo stesso
   sulla CPU.
4. **Letture asincrone.** `status` e `screen` si leggono con `map_async` e si
   raccolgono al frame successivo (`device.poll(Maintain::Poll)`): 1 frame di
   latenza, zero stalli. Si legge lo schermo solo quando `alive == 0`.

Alternativa più semplice ma peggiore: device dedicato + `copy_buffer_to_buffer`
verso un buffer di Bevy non è possibile (device diversi); servirebbe comunque
passare dalla CPU. Quindi la strada è il device condiviso.

Passo successivo utile: un kernel `emit` che genera i grani direttamente in GPU
(Box–Muller sull'hash per grano, deterministico) per evitare l'upload di 64 MB
per fotone da un milione di grani.

## 5. Memoria per milione di grani

| Voce | Per 10⁶ grani |
|---|---|
| grani (64 B) | 64 MB (61 MiB) |
| fotoni / status / schermo | trascurabili (112 + 32 B + 11.5 KB per fotone) |
| **Totale oggi** | **≈ 64 MB** |
| R4: copia ordinata ping-pong dei grani | +64 MB |
| R4: chiave di cella + indice ordinato (2 × u32) | +8 MB |
| R4: ψ locale stimata per grano (2 complessi, vec4) | +16 MB |
| R4: griglia (count, start, cursor) per M celle | 12·M B (es. M = 2²⁰ → 12 MB) |
| **Totale con R4** | **≈ 165 MB** |

Limiti wgpu: di default `max_storage_buffer_binding_size` = 128 MiB → 2 M grani
per binding. `GrainEngine::new` legge `adapter.limits()` e chiede di più solo se
serve e solo quanto l'adapter consente (altrimenti errore chiaro). Oltre si
spezza il buffer in più binding.

## 6. Regola locale R4 e griglia a counting sort (senza limite per cella)

R4: ogni grano stima il campo locale ψ(x_i) dai vicini entro una finestra h,
ψ̂(x_i) = Σ_j a_j K_h(x_i − x_j): una ricerca di vicini per ogni grano.

Il vecchio `pps_compute.wgsl` usava una griglia con `max_per_cell` slot fissi:
nelle celle dense i grani in più venivano **scartati in silenzio**, esattamente
dove la densità (cioè |ψ|²) è maggiore. Per R4 si usa invece un
**counting sort** completo, con memoria proporzionale a N e non a celle × capienza:

1. **key**: `cell = hash(floor(pos/h))` (tabella di M = potenza di 2 ≥ 2N, o una
   griglia densa se il volume è limitato).
2. **count**: `atomicAdd(count[cell], 1)` — nessun limite.
3. **scan**: prefisso esclusivo su `count` → `start[cell]` (Blelloch a più
   livelli o decoupled look‑back; M fino a milioni).
4. **scatter**: `dst = start[cell] + atomicAdd(cursor[cell], 1)`; si scrive
   l'indice (o direttamente il grano nel buffer ping-pong: migliore coerenza di
   memoria nella fase di query).
5. **query**: per le 27 celle (3D) intorno al grano si scorre
   `start[c] .. start[c+1]`: lunghezza variabile, costo proporzionale alla densità
   reale, mai troncata.
6. **update**: legge lo stato vecchio (buffer A), scrive il nuovo (buffer B),
   poi scambio (ping-pong): R4 legge i vicini, quindi non può aggiornare sul posto.

Riproducibilità: l'ordine dentro una cella dipende dall'ordine degli atomici;
per somme bit‑identiche si ordina la cella per indice di grano (ranking locale nel
workgroup, oppure radix sort su `(cell, indice)`) oppure si accumula ψ̂ in
virgola fissa come per lo schermo. Celle molto dense: un workgroup per cella
che carica i vicini a blocchi in memoria condivisa (tiling) invece di un thread
per grano che scorre liste lunghe.

## 7. Risultati del test CPU vs GPU

`cargo test --bin fotone gpu -- --nocapture` (Intel UHD Graphics, Vulkan, build dev).

`gpu_matches_cpu_double_slit`: doppia fenditura di default, 60 000 grani, seed 1
(primo seed in cui il fotone passa la parete), 33 passi da 7.3 (passo non
commensurabile: incroci a metà passo). Stessi grani iniziali presi da `Sim`,
stesso hash per grano.

| Grandezza | CPU (`Sim`) | GPU |
|---|---|---|
| p(passa la parete) | 0.239500 | 0.239500 |
| grani vivi dopo 15 passi | 14 370 | 14 370 (0 differenze vivo/diffratto) |
| max \|Δpos\| dopo 15 passi | — | **3.05·10⁻⁵** (soglia 10⁻³) |
| grani arrivati sullo schermo | — | 11 450 |
| correlazione profili d'intensità | — | **0.9999995** (soglia 0.999) |
| max \|ΔI\| (profili normalizzati) | picco 2.15·10⁻² | 7.19·10⁻⁵ |

Con le `sin/cos` native lo stesso test dava max |Δpos| = 1.19·10⁻³ e correlazione
0.99991: da qui la `sincos` polinomiale.

`gpu_photon_level_decisions`: con numero casuale 0.999999 il fotone viene
assorbito intero alla parete (`fate = WALL`, 0 grani vivi, 0 conteggi sullo
schermo); luce H contro un polarizzatore a 90° → `fate = POLARIZER`, 0 vivi.

`gpu_million_grains_throughput`: 10⁶ grani (61 MiB), 100 tick da 2.0 in un solo
command buffer, poi 5 passi di `Sim::step` sugli stessi grani (rayon):

| Condizione | GPU (probe+decide+advance) | CPU `Sim::step` |
|---|---|---|
| macchina scarica | **11–13 ms/tick** | 16.5 ms/passo |
| macchina carica (compilazioni in parallelo) | 17–39 ms/tick | 61–172 ms/passo |

Lettura onesta: su questa GPU integrata il volo libero è limitato dalla banda di
memoria (64 B letti + 64 B scritti per grano in `advance`, più 64 B letti in
`probe`: ~190 MB per tick con 10⁶ grani), quindi il guadagno sulla sola
traslazione è piccolo. Il vantaggio vero arriva (a) dall'eliminazione della
copia per frame verso il renderer, che oggi sulla CPU costa quanto il passo, e
(b) da R4, dove la ricerca dei vicini è calcolo e non banda. Ottimizzazioni
ovvie per il volo libero: saltare `probe` quando nessun fotone ha un elemento
entro `dist` (l'host conosce z del fronte), e dividere il grano in parte calda
(`pos`, `state`: 16 B) e parte fredda (letta solo agli incroci) → ~4× meno
traffico.

Nel test di confronto (60 000 grani, 33 passi) i tempi non sono significativi:
la GPU include due letture complete dei grani fatte apposta per il test.
