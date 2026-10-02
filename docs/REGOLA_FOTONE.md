# La regola del fotone: una sola regola locale per ogni pallina

Laboratorio: cartella `lab/` (Python). Tutti i numeri qui sotto vengono da esecuzioni reali degli
script elencati in fondo (seed fissi, 8 core). Nessun "PASS" è dichiarato senza una corsa che lo mostri.

---

## 1. In una frase

Un fotone è fatto di **tante palline-onda** e di **una pallina-fotone**. Le palline-onda volano dritte
e portano con sé da dove sono partite e con che ampiezza. La pallina-fotone, a ogni passo, guarda
**solo le palline-onda che le passano accanto**, somma le loro onde e si muove lungo la direzione in cui
"scorre" l'onda in quel punto. Dove sbatte sullo schermo, lì c'è il clic.

Con questa regola, e senza nessuna scelta globale, si ottengono: il fascio gaussiano libero (errore 0,2%),
la doppia fenditura (correlazione 0,995 con la teoria esatta, errore sulle altezze delle frange 4%),
il cancellino "which-path" con polarizzatori, la legge di Malus e la costruzione del disegno un clic
alla volta.

---

## 2. Il quadro (per chi non è fisico)

* L'**onda** di un fotone viene rappresentata da palline. Ogni pallina è un pezzetto di un'"ondina"
  circolare (onda di Huygens) partita da un punto preciso: la sua **origine**. La pallina sa da dove
  è partita, in che direzione vola, con che ampiezza e che polarizzazione (un vettore di Jones).
* Il **fotone vero e proprio** — la cosa che alla fine fa "clic" in un solo punto — è una pallina
  speciale. Non ha energia propria: si fa "trasportare" dall'onda, come un tappo di sughero dalla
  corrente di un fiume (è la meccanica di de Broglie–Bohm, detta anche "onda pilota").
* L'unica casualità del modello è **dove nasce** la pallina-fotone alla sorgente (più probabile dove
  l'onda è più intensa) e un **lancio di moneta locale** quando la pallina-fotone tocca un oggetto
  che assorbe (bordo di una fenditura, polarizzatore assorbente). Tutto il resto è deterministico.
* Niente decisioni globali: nessuno guarda "tutto il disegno" per scegliere dove far comparire il clic.

---

## 3. La regola, passo per passo

### 3.1 Cosa fa una pallina-onda (sempre la stessa regola)

1. **Vola dritta** alla velocità della luce. La sua fase cresce di 2π ogni lunghezza d'onda λ di strada
   (fase = `phase0 + k·|pos − origin|`).
2. Porta con sé, fissati al momento della nascita: `origin`, `dir`, `amp` (ampiezza complessa),
   `jones` (polarizzazione).
3. **Quando tocca un oggetto** (sorgente, parete con fenditure, polarizzatore, beam-splitter):
   * se tocca una parte opaca **muore**;
   * altrimenti l'oggetto la **ri-emette**: diventa l'origine di un ventaglio di nuove palline
     con direzioni regolari θ_j = −Θ + (j+½)·dθ che coprono tutto il semipiano in avanti (Θ ≈ 90°).
     L'ampiezza del ventaglio è `a = ψ_in(punto) · T(punto) · ds`, dove ψ_in è l'onda che arriva lì
     (letta localmente) e T è la trasmissione dell'oggetto in quel punto (matrice di Jones × fase).
   * Le origini sono su una **griglia regolare** di passo `ds` lungo l'oggetto (non a caso!).

### 3.2 Cosa fa la pallina-fotone a ogni passo

1. **Ascolta i vicini.** Per ogni origine *i* dell'ultimo oggetto attraversato, cerca le **due palline
   di quell'origine i cui raggi le passano a cavallo** (una a sinistra, una a destra). Il peso è
   "a cappello": lineare tra le due, `1−f` e `f` (f = posizione frazionaria tra i due raggi).
   Se nessuna pallina di quell'origine le passa accanto (morte contro un ostacolo, o fuori dal
   ventaglio) quell'origine non conta.
2. **Somma le ondine.** Ogni pallina ascoltata contribuisce con l'ondina della sua origine valutata
   nella posizione del fotone: `peso · a_i · G(P − origin_i)` e il suo gradiente. G è l'ondina di
   Huygens esatta in 2D (Rayleigh–Sommerfeld): `G = (ik/2)·Δz·H1(kR)/R`. Le due componenti di
   polarizzazione restano separate: ψ = (ψ_x, ψ_y).
3. **Calcola la corrente** `j = Σ_pol Im(ψ* ∇ψ)` (somma sulle polarizzazioni: così due onde con
   polarizzazioni ortogonali non interferiscono, come deve essere).
4. **Si muove alla velocità della luce lungo j**: in coordinate `dx/dz = j_x / j_z`
   (integrazione Runge–Kutta 4, passo `dz = min(2λ, max(0,01λ, 0,1·distanza dall'oggetto))`).
5. Se non sente **nessun** vicino (o j_z ≤ 0) **continua nella sua direzione**: un fotone fatto di
   una sola pallina è un raggio dritto.
6. **Oggetti:**
   * parete opaca → la pallina-fotone muore (fotone assorbito, tutto intero);
   * **regola del flusso** (oggetti che assorbono: bordi di fenditura, polarizzatori assorbenti):
     subito dietro l'oggetto (a distanza δ) la pallina-fotone sopravvive con probabilità
     `p = min(1, Φ_out / Φ_in)`, dove Φ_in è il flusso in avanti dell'onda in arrivo e Φ_out quello
     dell'onda uscente, **entrambi letti nella sua posizione**. Per un polarizzatore p = cos²θ;
     per una fenditura p ≈ 1 al centro e scende verso i bordi (vedi §5);
   * oggetti che non assorbono (beam-splitter polarizzatore, prisma): nessuna moneta; il primo
     passo dietro l'oggetto si fa già lungo la corrente dell'onda uscente;
   * schermo → **clic** nel punto d'impatto.

### 3.3 Parametri usati

| parametro | valore | ruolo |
|---|---|---|
| λ | 1 | unità di lunghezza |
| `ds` (passo delle origini) | 0,02 λ su fenditure con bordi netti; 0,1 λ su sorgenti/oggetti lisci | risoluzione dell'integrale di Huygens |
| `δ` (distanza della regola del flusso) | 5·ds (0,1 λ fenditure; 0,5 λ oggetti lisci) | vicino al bordo, ma ≫ ds |
| `dθ` (passo delle direzioni) | 0,01 rad (ma 0,1 dà lo stesso risultato, variante H) | quante palline per origine |
| Θ (apertura del ventaglio) | π/2 − 10⁻³ | tutto il semipiano avanti: indispensabile nel campo vicino |
| finestra | "cappello" = le 2 palline a cavallo per origine (larghezza = spaziatura locale dei raggi R·dθ) | adattiva, numero fisso di vicini |
| palline per fotone (doppia fenditura) | 200 origini × 314 direzioni = 62.800 (con dθ=0,1: 6.200) | |
| passo del fotone | RK4, dz = min(2λ, max(0,01λ, 0,1·Δz)) | |
| geometria doppia fenditura | fascio gaussiano w0=7, fenditure larghe 2 a ±4, schermo a 160 | come `pilot.py` |

---

## 4. Compito 1 — perché la guida con il campo ESATTO dava solo 0,95

Script: `lab/exact_guide2.py` (6000 clic; campionatore ideale con 6000 clic: corr 0,9972±0,0010,
errore frange rms 3,2%, max 5,4%). Riferimento: flusso j_z esatto sullo schermo (Rayleigh–Sommerfeld).
"corr" = correlazione grossolana (bin da 8), come prima.

| passo | cosa cambia | corr | err. frange rms / max | tempo |
|---|---|---|---|---|
| 0 | diagnosi vecchia (`exact_guide.py`): campo approssimato, legge `v=Im(ψ*∂xψ)/(k|ψ|²)`, partenza \|U_inc\|² a z=0,05, RK2 dz=0,5 | 0,955 (0,951 vs parassiale) | 15,5% / 27,6% | 47 s |
| 1 | + legge del flusso dx/dz = j_x/j_z | 0,952 | 17,5% / 31,4% | 117 s |
| 2 | + campo esatto (Hankel) e RK4 adattivo; partenza ancora \|U_inc\|² | 0,969 | 14,3% / 27,6% | |
| 3 | + partenza **equivariante** j_z(x, z0=0,05) | **0,9959** | 2,4% / 3,5% | |
| 4 | partenza j_z a z0=0,5 | 0,9968 | 2,3% / 3,4% | |
| 5 | partenza j_z a z0=2 | 0,9968 | 2,7% / 4,4% | 52 s |
| 6 | come 5 con passi dimezzati (convergenza) | 0,9981 | 2,3% / 4,4% | 101 s |
| 7 | **partenza \|U_inc\|² + regola del flusso locale** a z0=0,05 | **0,9978** | 4,2% / 6,7% | 156 s |
| 8 | come 5 ma legge vecchia `j_x/(k|ψ|²)` | 0,9940 | 2,2% / 3,7% | 52 s |

**Diagnosi.**
1. **La causa principale è la distribuzione di partenza** (passo 2 → 3: da 0,969 a 0,996).
   Subito dietro una fenditura con bordi netti il flusso di probabilità j_z **non** è |U_inc|²:
   vicino ai bordi è impoverito (fino a 0,26 volte) e al centro leggermente arricchito (fino a 1,05),
   e nel modello di Kirchhoff/Rayleigh–Sommerfeld solo il 92% del flusso incidente passa davvero.
   Anche uno scarto del 2% sulla distribuzione iniziale diventa un 20–30% su una singola frangia,
   perché le traiettorie di Bohm non si incrociano: la frangia centrale è alimentata proprio dalle
   traiettorie che partono vicino ai bordi interni delle fenditure.
2. Il passo grosso (dz=0,5 partendo da 0,05) e il campo approssimato vicino alle fenditure contano
   anche loro (passi 0–1).
3. La legge vecchia `j_x/(k|ψ|²)` è sbagliata solo al secondo ordine nell'angolo (0,994 invece di 0,997):
   la legge giusta per contare i clic su uno schermo piatto è `dx/dz = j_x/j_z`.
4. Il vecchio riferimento **parassiale** non è la verità: con quello esatto ha correlazione 0,9917,
   quindi da solo impedisce di superare ~0,99. Tutti i confronti ora usano il flusso esatto.

La punta di partenza non si può "aggiustare" globalmente (servirebbe conoscere la distribuzione
di tutto il fascio). Il rimedio **locale** è la regola del flusso (§3.2.6): la pallina-fotone entra
con |U_inc|², e subito dietro la fenditura sopravvive con probabilità min(1, j_z/(k|U_inc|²)) letta
nella sua posizione. Così la popolazione sopravvissuta è proprio j_z (passo 7: 0,9978).
Fisicamente: la parte d'onda che vicino al bordo non scorre in avanti (onde evanescenti) non può
trasportare il fotone; la moneta rappresenta l'assorbimento/riflessione sul bordo, che un modello
di Kirchhoff non descrive.

**Obiettivo ≥ 0,98 con ≥ 3000 clic: raggiunto** (0,9978 con 6000 clic in 156 s).

---

## 5. Compito 2 — stimare l'onda dai soli vicini

Script: `lab/neighbour_study.py`, stessa doppia fenditura, riferimento esatto.
Campionatore ideale: 3000 clic → 0,9945±0,0018 (rms 4,4%, max 7,5%); 6000 clic → 0,9972±0,0010 (3,2%, 5,4%).

| variante | palline/fotone | clic | corr | χ²/gdl | err. frange rms / max | tempo |
|---|---|---|---|---|---|---|
| A palline a caso (come `pilot.py`/Rust oggi), finestra gaussiana h=1, onda piana dalla pallina (E1), senza regola del flusso | 4.900 | 2000 | 0,931 | 5,7 | 15,7% / 20,6% | 283 s |
| B palline a caso, finestra h=1, ondina esatta dall'origine (E2), senza regola del flusso | 4.900 | 1400 | 0,965 | 2,4 | 13,0% / 24,7% | 265 s |
| D reticolo, cappello, E1 (solo fase+direzione della pallina) | 125.600 | 4386 | 0,9575 | 6,0 | 16,0% / 29,6% | 185 s |
| E reticolo, cappello, **E2**, ds=0,01, δ=0,05 | 125.600 | 5905 | **0,9954** | 1,26 | **3,4% / 4,2%** | 224 s |
| F come E ma **senza** regola del flusso | 125.600 | 6000 | 0,9766 | 5,1 | 12,0% / 21,4% | 263 s |
| G come E, ds=0,04, δ=0,2 | 31.400 | 5741 | 0,9938 | 1,31 | 4,8% / 7,0% | 38 s |
| H come E, dθ=0,1 (31 direzioni) | 12.400 | 5926 | 0,9951 | 1,33 | 3,7% / 4,8% | 262 s |
| **I (regola finale)**: E2, ds=0,02, δ=0,1 | 62.800 | 5858 | **0,9949** | 1,34 | **3,9% / 4,7%** | 98–135 s |
| J come I, primo passo lungo la corrente uscente | 62.800 | 5804 | 0,9932 | 1,47 | 5,2% / 8,8% | 131 s |

Nota A/B: con palline a caso e finestra fissa h=1 la regola del flusso ucciderebbe quasi tutti i
fotoni (il ventaglio vicino alla fenditura è più stretto della finestra, l'onda stimata è ~0), quindi
A e B girano senza. Sono lente (molti vicini per fotone), perciò usano meno fotoni. Corse di prova con
3000 fotoni: A 0,938 (RK4), B 0,965 (RK2) — stesso quadro. Il vecchio `pilot_batch.py` dava 0,88–0,89
(legge `j_x/(k|ψ|²)`, Eulero, 1900 clic).

**Cosa si impara.**
* **Il rumore delle palline a caso** è il primo nemico (A, B): con poche migliaia di palline
  messe a caso la somma di fasi è rumorosa, e una finestra fissa h=1 **sottostima di molto** l'onda
  vicino alle fenditure (il ventaglio è più stretto della finestra). Rimedio: **origini e direzioni
  su un reticolo regolare** + finestra a cappello larga quanto la spaziatura locale dei raggi
  (= sempre 2 vicini per origine). Con questa scelta i pesi di ogni origine sommano esattamente a 1
  e la stima coincide con l'integrale di Huygens discretizzato con passo ds (errore < 10⁻⁴,
  verificato punto per punto).
* **L'onda piana della pallina non basta vicino alle fenditure** (D): lì R è di pochi λ e un'onda
  circolare non è un'onda piana; anche il flusso j_z per la regola del bordo viene sbagliato
  (84% assorbiti invece di 78%). Nel campo lontano E1 ≈ E2 (errore 0,3% sull'ampiezza a z=50).
  Usare l'**ondina esatta dall'origine** (E2) — informazione che la pallina già porta (`origin`) —
  risolve il problema.
* Il numero di direzioni non conta per E2 (H: 12.400 palline = stesso risultato di 125.600).
  Conta il passo delle origini `ds` (G un po' peggiore).
* La regola del flusso è indispensabile (F: senza, 0,977 e frangia centrale +21%).
* **Obiettivo ≥ 0,97 e errore frange < 10%: raggiunto** (I: 0,9949, rms 3,9%, max 4,7%, 5858 clic,
  statisticamente indistinguibile dal campionatore ideale; χ²/gdl = 1,34 su ~28 gradi di libertà).

---

## 6. Compito 3 e 4 — la STESSA regola sugli altri esperimenti

Script: `lab/validate.py`, `lab/isolated_ball.py`.

| esperimento | misurato | teoria | esito |
|---|---|---|---|
| Fascio libero w0=3, larghezza a z = zR | 4,247 | 4,243 (RS esatto 4,252) | **PASS** (+0,10%) |
| … a z = 2 zR | 6,718 | 6,708 (6,725) | **PASS** (+0,15%) |
| … a z = 5 zR | 15,325 | 15,297 (15,340) | **PASS** (+0,18%) |
| Doppia fenditura (regola finale, 5858 clic) | corr 0,9949; frange rms 3,9%, max 4,7% | ideale 0,9972±0,0010 | **PASS** |
| Frazione assorbita alla parete | 78,9% | 76,2% geometrica; 78,1% modello RS | ok (vedi §7) |
| Which-path, polarizzatori paralleli H,H | V = 0,925; corr 0,9953 | V = 0,907 | **PASS** (> 0,9) |
| Which-path, polarizzatori ortogonali H,V | V = −0,023; corr 0,9851 (ideale 0,990±0,003) | V = 0,029 | **PASS** (< 0,15) |
| Malus M1, polarizzatore assorbente, θ = 0 / 22,5 / 45 / 67,5 / 90° | 0,9962 / 0,8565 / 0,4885 / 0,1452 / 0,0000 | 1 / 0,8536 / 0,5 / 0,1464 / 0 | **PASS** (entro 3σ o 1%) |
| Malus M2, beam-splitter + assorbitore (deterministico) | 0,9993 / 0,8470 / 0,5040 / 0,1430 / 0,0000 | idem | **PASS** |
| Polarizzazione dell'onda dove si trova il fotone dopo il polarizzatore | \|⟨a\|ψ⟩\|²/\|ψ\|² = 1,0000 | 1 | **PASS** |
| Costruzione un clic alla volta: corr dopo 50 / 200 / 1000 / 3000 / 5909 fotoni | 0,57 / 0,88 / 0,975 / 0,992 / 0,998 | ideale 0,78±0,08 / 0,93±0,02 / 0,984±0,006 / 0,994±0,002 / 0,997±0,001 | **PASS** (solo rumore statistico) |
| Pallina isolata (N=1), 50 palline a caso, E2 ed E1 | scarto dal raggio dritto ≤ 4·10⁻¹³ λ | 0 | **PASS** |
| Pallina senza nessun vicino | scarto 0 | 0 | **PASS** |

Dettagli:
* **Fascio libero.** La sorgente è la vita del fascio: emette il ventaglio a reticolo con ampiezza
  gaussiana, e la pallina-fotone nasce con probabilità |U|². Le traiettorie si allargano come
  w(z)/w0: l'errore sulla larghezza è 0,1–0,2%, molto sotto il 3% richiesto.
* **Which-path.** La luce arriva polarizzata a 45°; su ogni fenditura un polarizzatore. Le palline che
  passano dalla fenditura *i* portano il vettore di Jones proiettato. La pallina-fotone usa la corrente
  **sommata sulle polarizzazioni**: con polarizzatori ortogonali le due ondine non interferiscono e
  le frange spariscono (V ≈ 0), con quelli paralleli restano (V = 0,925). La metà dei fotoni è
  assorbita dai polarizzatori con la stessa regola del flusso (p = 1/2). Visibilità misurata come
  (I_centro − I_buio)/(I_centro + I_buio) su finestre di ±2,5 λ.
* **Malus, due versioni:**
  * **M1 (polarizzatore assorbente):** la regola del flusso dà p = Φ_out/Φ_in = |⟨a|J⟩|² = cos²θ
    (misurato r = 0,849–0,851 a 22,5°). È un lancio di moneta **della sola pallina-fotone**,
    con la sua polarizzazione locale. Il fotone passa tutto intero o non passa.
  * **M2 (beam-splitter polarizzatore + assorbitore sull'uscita scartata):** nessuna moneta.
    L'oggetto divide l'onda in due fasci con polarizzazioni a e b che si separano (±0,25 rad). La
    pallina-fotone segue la corrente e finisce in uno dei due fasci **a seconda di dove si trovava**
    nel fascio: è deterministico, e poiché la posizione iniziale è distribuita come |ψ|², la frazione
    che finisce nel fascio "a" è esattamente cos²θ. È la risposta "alla Bohm" completa.
    Per farlo funzionare il primo passo dietro l'oggetto va fatto lungo la corrente uscente:
    volando dritti per δ il fascio inclinato si sposta di δ·α e la divisione si sbilancia
    (misurato: 0,835 invece di 0,854 prima della correzione).

**Perché un polarizzatore assorbente non si può fare con una regola locale deterministica "da
ferma".** Se la decisione passa/non passa fosse una funzione fissa dello stato della sola pallina
(posizione, Jones), per ottenere cos²θ per ogni θ bisognerebbe sapere in che quantile del fascio si
trova la pallina: informazione globale. La parte assorbita dell'onda sparisce, quindi non c'è una
"strada" per una traiettoria. Le due vie locali sono: (M1) una moneta locale della pallina, oppure
(M2) modellare l'assorbitore come un canale a parte (dove la parte scartata dell'onda va davvero), e
allora è la dinamica stessa a decidere, in modo deterministico e locale. Entrambe validate.

---

## 7. Cosa non funziona, o funziona solo in parte (onestà)

1. **La regola del flusso al bordo è una correzione, non fisica "prima".** Nasce dal fatto che il modello
   di Kirchhoff/Rayleigh–Sommerfeld mette un campo che salta di colpo al bordo della fenditura: solo il
   92% del flusso incidente passa, e al centro della fenditura il rapporto supera 1 (fino a 1,05–1,09).
   Una moneta locale può solo *togliere* fotoni, non aggiungerli: dove il rapporto supera 1 viene
   tagliato a 1. Effetto: frazione assorbita 78,9% invece di 78,1% e un residuo di qualche % sulla
   frangia centrale (−3…−5%, entro le fluttuazioni ma sistematicamente dello stesso segno).
   Un modello di bordo più realistico (bordi morbidi, o l'onda riflessa dalla parete) renderebbe la
   moneta superflua o più piccola.
2. **Due varianti del "primo passo" dietro un oggetto.** Per le fenditure (incidenza normale, bordi
   netti) funziona meglio volare dritti per δ e poi applicare la moneta (I: 0,9949) che fare il primo
   passo lungo la corrente uscente (J: 0,9932). Per gli oggetti che inclinano il fascio (beam-splitter)
   serve invece il primo passo lungo la corrente uscente. Sono due approssimazioni coerenti dello
   stesso requisito (la popolazione a distanza δ deve essere j_z); non sono riuscito a unificarle in
   una sola formula che vada bene per entrambi con δ finito.
3. **E1 (solo fase e direzione della pallina) fallisce vicino agli oggetti** (0,958). La regola finale
   usa l'origine della pallina per calcolare l'ondina esatta. Localmente è lecito (è stato della pallina
   stessa, già presente in `Grain`), ma vuol dire che le palline fanno da **messaggeri**: in spazio
   libero ogni origine ha sempre due palline accanto al fotone, quindi il fotone di fatto somma le
   ondine di **tutte** le origini dell'ultimo oggetto. Il costo è O(numero di origini) per passo.
   Il beneficio della località si vede quando qualcosa uccide delle palline (un ostacolo tra
   fenditura e fotone spegne le origini che non "vedono" il fotone).
4. **Palline a caso = rumore.** L'emissione casuale di oggi in Rust (`huygens()`: θ uniforme in ±0,45,
   posizioni da |U|²) non basta: va sostituita da un reticolo regolare di origini e direzioni, e il
   ventaglio deve coprire tutto il semipiano (con ±0,45 rad il campo vicino è sbagliato).
5. **Solo 2D.** Le fenditure sono lunghe in y. In 3D le origini diventano una griglia 2D sull'apertura
   (costo ∝ numero di origini²) e l'ondina diventa e^{ikR}/R.
6. **Malus con polarizzatore assorbente** usa una moneta locale (M1). La versione deterministica (M2)
   richiede di modellare il canale dell'assorbitore.
7. **Entanglement** (coppie di fotoni, Bell) resta non locale per forza (teorema di Bell) ed è gestito
   in `src/photon/entangle.rs`; non fa parte di questo lavoro.
8. MIW (R2/R3) non è stato ripreso: con la regola a onda pilota non serve.

---

## 8. Specifica per il porting in Rust (`src/photon/physics.rs`)

### 8.1 Dati

```rust
pub struct Grain {            // pallina-onda (già esistente, con 2 campi in più)
    pos: Vec3, dir: Vec3,
    origin: Vec3,             // punto dell'ultima (ri)emissione
    phase0: f32,              // non serve più per la guida (la fase viene da |P-origin|), resta per il disegno
    amp: C,                   // NUOVO significato: ampiezza della propria origine (a_i), uguale per tutte le palline dell'origine
    jones: Jones,             // polarizzazione portata (già proiettata da eventuali filtri)
    alive: bool, next: u8,
    fan: u16,                 // NUOVO: indice dell'oggetto che l'ha emessa (0 = sorgente, 1 = parete, ...)
    src: u32,                 // NUOVO: indice dell'origine nel reticolo di quell'oggetto
    dir_idx: u16,             // NUOVO: indice j della direzione nel ventaglio
}
pub struct Pilot {            // la pallina-fotone
    pos: Vec3, slope: f32,    // slope = dx/dz corrente
    fan: u16,                 // ventaglio che la guida ora
    alive: bool,
}
pub struct Fan { ze: f32, origins: Vec<f32> /* x_i */, amp: Vec<[C;2]> /* a_i ⊗ jones */, th: f32, dth: f32 }
```

### 8.2 Emissione (a ogni oggetto, incluso la sorgente)

* Origini `x_i = x_lo + (i+½)·ds` sulle parti trasparenti; `ds = 0,02λ` se l'oggetto ha bordi netti,
  `0,1λ` se è liscio. Direzioni `θ_j = −Θ + (j+½)·dθ`, `Θ = π/2 − 10⁻³`, `dθ = 0,01` (0,1 va bene per la
  guida; più fitto serve solo a disegnare).
* `a_i = ψ_in(x_i) · T(x_i) · ds` (vettore di Jones complesso). Alla sorgente ψ_in è il profilo
  gaussiano; a un oggetto a valle ψ_in(x_i) si ottiene con la stessa somma di §8.3 fatta nel punto
  x_i (le palline che arrivano lì).
  Filtro di polarizzazione f sulla fenditura: `a_i = f ⟨f|J_in⟩ · inc(x_i) · ds`.
  Beam-splitter a angolo θ con deviazione ±α: `a_i = (1/√cosα)(a⟨a|ψ⟩e^{+ik sinα x} + b⟨b|ψ⟩e^{−ik sinα x}) ds`.
* La pallina-fotone nasce alla sorgente in `x ~ |U_src(x)|²` (unico numero casuale globale del modello,
  ma è un evento *della sorgente*, locale).

### 8.3 Guida (per ogni passo della pallina-fotone, sotto-passi RK4 in z)

```
fn field(P, fan) -> ([C;2] psi, [C;2] dpsi_dx, [C;2] dpsi_dz):
  dz = P.z - fan.ze
  for i in origins:
     t = atan2(P.x - x_i, dz)
     q = (t + Θ)/dθ - 0.5 ; j = floor(q) ; f = q - j
     w = (1-f)·alive(i,j) + f·alive(i,j+1)          // le 2 palline a cavallo (0 se morte o assenti)
     if w == 0 { continue }
     R = hypot(P.x - x_i, dz) ; u = k R
     H0, H1 = hankel(u)          // u > 12: serie asintotica a 4 termini (errore < 2e-5);
                                 // u ≤ 12: J0,J1,Y0,Y1 (Abramowitz–Stegun 9.4) — serve solo entro 2λ dall'oggetto
     F  = H1/R ; Fp = (k H0 - 2 H1/R)/R
     G  = (i k/2) dz F
     Gx = (i k/2) dz Fp (P.x - x_i)/R
     Gz = (i k/2) (F + dz Fp dz/R)
     psi += w G a_i ; dpsi_dx += w Gx a_i ; dpsi_dz += w Gz a_i
slope(P) = Σ_pol Im(conj(psi)·dpsi_dx) / Σ_pol Im(conj(psi)·dpsi_dz)   // se il denominatore ≤ 0 o psi≈0: slope invariato
passo:  dz = min(2λ, max(0,01λ, 0,1·(P.z − fan.ze)))  ; RK4 su x(z)
```

* **Nella versione a impulso** (tutte le palline di un fotone emesse insieme, come ora in Rust) le
  palline di un'origine stanno su un arco di raggio c·t: per E2 conta solo che il **raggio** della
  pallina passi accanto al fotone, non dove si trova lungo il raggio. La ricerca dei vicini deve quindi
  accettare palline entro ±(differenza di cammino massima + lunghezza dell'impulso) lungo il raggio
  (qui ~3–5λ). Con gli indici `(fan, src, dir_idx)` la ricerca è una divisione, non una ricerca spaziale.
* Costo: O(N_origini) per sotto-passo: 200 origini × ~600 valutazioni per fotone nella doppia
  fenditura ≈ 10⁵ ondine per fotone, banale per rayon.

### 8.4 Eventi agli oggetti

```
parete (parte opaca)           -> pilot.alive = false  (fotone assorbito intero)
oggetto che assorbe            -> vola dritto per δ = 5·ds; r = Φ_out/Φ_in con
                                  Φ_in  = Σ_pol Im(conj ψ_in · ∂zψ_in) nel punto d'arrivo (ventaglio vecchio)
                                        (per un fascio che arriva normale: k·|ψ_in|²)
                                  Φ_out = stesso con il ventaglio nuovo a (x, ze+δ)
                                  sopravvive se rng_locale < min(1, r)
oggetto che non assorbe (PBS)  -> primo passo implicito: x1 = x + δ·slope_nuovo(x1), 2 iterazioni; nessuna moneta
schermo                        -> clic in pilot.pos (sostituisce il campionamento globale di finalize())
```

Da togliere in Rust: `wall_decision()` (probabilità di passare calcolata su tutto il fotone),
`measure_single()` (moneta globale del polarizzatore) e il campionamento del clic da tutto il
disegno in `finalize()`. Da tenere: `entangle.rs` per le coppie.

---

## 9. Esecuzioni finali e come rieseguire

Tutti gli script stanno in `lab/`, usano seed fissi, stampano con `flush=True`, usano 8 processi
(`multiprocessing`). I risultati (clic) sono salvati in `lab/out/*.npy`.

| comando | cosa fa | tempo misurato |
|---|---|---|
| `python lab/exact_guide2.py` | compito 1, passi 5 e 7 (`... 0 3` per altri passi) | ~3,5 min |
| `python lab/photon_rule.py` | la regola finale sulla doppia fenditura | 135 s |
| `python lab/neighbour_study.py` | variante I (`A B D E F G H J` una per volta) | 98–135 s (A 283 s, B 265 s, E/F/H 220–265 s) |
| `python lab/validate.py free` | fascio libero | 83 s |
| `python lab/validate.py which` | which-path | 186–277 s |
| `python lab/validate.py malus1` / `malus2` | Malus M1 / M2 | 187 s / 306 s |
| `python lab/validate.py build` | costruzione un clic alla volta | 138 s |
| `python lab/isolated_ball.py` | pallina isolata | < 5 s |

File: `lab/optics.py` (riferimenti esatti e punteggi: campo di Rayleigh–Sommerfeld 2D, flusso,
correlazione, errori sulle frange, χ², visibilità, campionatore ideale), `lab/photon_rule.py`
(la regola), `lab/exact_guide2.py`, `lab/neighbour_study.py`, `lab/validate.py`,
`lab/isolated_ball.py`. I vecchi `rules.py`, `pilot.py`, `pilot_batch.py`, `exact_guide.py` non sono
stati modificati.

---

## 10. Regola pilota nel programma (`fotone`, Rust)

La regola dei §3 e §8 è ora una modalità del programma: `src/photon/pilot.rs` (la regola),
`src/photon/pilot_validation.rs` (i test). Si sceglie nel pannello nascosto (Tab / ☰: "rule:
pilot (local, 2-D) / classic") o da riga di comando `--regola pilota|classica`; è la scelta
predefinita per fascio libero, doppia fenditura e polarizzatori (l'entanglement resta classico).
La modalità classica è intatta: `Lab::preset` la usa ancora e i suoi 14 test passano invariati.

**Cosa si vede.** Le palline-onda (colore = fase, come prima) partono dalla sorgente come un ventaglio
a reticolo; quando il fronte arriva a un oggetto nasce il ventaglio nuovo (200 origini × 314 direzioni
= 62.800 palline dietro le fenditure). La pallina-fotone è la stella bianca con la sua scia: viaggia
dentro la sua onda, muore intera sulla parete o al bordo (moneta locale), e il clic è dove tocca lo
schermo. Il pannello scrive "2-D": la regola vive nel piano x–z; lo spessore in y delle palline-onda
(±1,5) e dei clic (±12) è **solo disegno**. Le palline che escono di lato oltre |x| = 100 non si
disegnano più (non cambia nulla per la pallina-fotone, vedi sotto). Il tempo reversibile funziona
(stato clonabile, `Arc` per i ventagli condivisi, stima della dimensione dei keyframe aggiornata).
Nessuna passata globale in questa modalità: `wall_decision`, `measure_single` e il campionamento del
clic da tutto il disegno in `finalize` non vengono mai chiamati.

**Test Rust** (`cargo test --bin fotone pilot -- --nocapture --test-threads=1`, profilo dev
opt-level 3, seed fissi, i5-10300H 4 core):

| test | Rust | Python (§5–§6) | soglia | tempo |
|---|---|---|---|---|
| Hankel H0, H1 contro scipy (32 punti, u = 0,01…1000) | errore relativo max 1,85·10⁻⁷ / 2,18·10⁻⁷ (a u = 12) | — | < 10⁻⁶ | < 0,1 s |
| pallina isolata: nessun vicino / una sola origine | scarto dal raggio 1,2·10⁻¹⁴ λ / 2,1·10⁻¹⁴ λ | ≤ 4·10⁻¹³ λ | < 10⁻⁶ λ | < 0,1 s |
| fascio libero w0 = 3, 20.000 fotoni, larghezza a z = 0,5 / zR / 2zR / 5zR | +0,31% / +0,38% / +0,42% / +0,44% (stat. ±0,5%) | +0,10% / +0,15% / +0,18% (4000 fotoni) | 3% | 15 s |
| doppia fenditura, 16.384 fotoni, 3210 clic | corr 0,9930; frange rms 2,8%, max 5,1% | 0,9949; 3,9% / 4,7% (5858 clic) | ≥ 0,98; ≤ 8% | 18 s |
| … frazione assorbita alla parete | 78,9% (persi alla sorgente 0,04%) | 78,9% | ±2 punti | |
| which-path paralleli / ortogonali (1689 / 1682 clic) | V = 0,924 / −0,043 (teoria 0,907 / 0,029) | 0,925 / −0,023 | > 0,85 / < 0,15 | 34 s |
| Malus, polarizzatore assorbente 0° / 45° / 90° (≈ 2000 fotoni ciascuno) | 1,0000 / 0,5145 / 0,0000 (σ a 45° = 0,011: 1,3σ) | 0,9962 / 0,4885 / 0 | max(3σ, 1%) | 7 s |
| corsa "live" (palline-onda memorizzate) = corsa senza palline, 20 fotoni | identiche bit per bit (1133 passi) | — | | 0,2 s |
| tempo all'indietro e ripetizione, 300 fotogrammi | esatti | — | | 0,5 s |
| costo di un passo a velocità ×1 (62.800 + 62.744 palline) | 0,30–0,37 ms medio, 4–5 ms peggiore | — | < 50 ms | |

Il programma vero (`fotone.exe --esperimento fenditure --fps`, build dev, velocità ×1) gira a
92–121 fotogrammi/s con 62.744–125.544 palline-onda in volo (tutto incluso: fisica, disegno,
interfaccia), nessun panic; memoria privata ogni 5 s: 614 → 676 → 677 → 681 MB (e in un'altra
corsa 612 → 676 → 681 → 678 → 675 MB): piatta dopo il riempimento della storia del tempo.

**Differenze rispetto al Python, dette chiaramente.**
1. **Φ_in letto nello stesso punto di Φ_out.** Il §8.4 dice "Φ_in nel punto d'arrivo"; così il primo
   Malus in Rust dava 0,9890 a 0° (fallito): il polarizzatore qui è a z = 20 con la vita del fascio a
   z = 0, e in δ = 0,5 il fascio si allarga dello 0,8% (d ln w/dz · δ), quindi Φ_out/Φ_in < 1 per
   tutti (nel Python il polarizzatore era a z = 5: 0,3%, da cui 0,9962). Ora Φ_in è il flusso delle
   palline in arrivo (quelle che toccano l'oggetto) valutato **nella posizione della pallina-fotone
   dietro l'oggetto**, come Φ_out — è la lettera del §3.2 ("entrambi letti nella sua posizione").
   Alle fenditure (alla vita del fascio) non cambia nulla: stessi numeri prima e dopo.
2. **Hankel a 6 termini** per kR ≥ 12 (il Python ne usa 4: errore ~7·10⁻⁶ a u = 12, non passerebbe
   10⁻⁶); sotto 12 Abramowitz–Stegun 9.4.1–9.4.6 (9.4.3/9.4.6 tra 3 e 12).
3. **La sorgente non è alla vita del fascio** (doppia fenditura: sorgente a z = −40, vita sulla parete;
   polarizzatori: −20 e 0). Il ventaglio della sorgente porta il fascio gaussiano 2-D complesso su
   quel piano; l'onda sulla parete è quindi propagata dalle palline, non scritta a mano: la sua
   figura sullo schermo coincide con quella esatta del Python (correlazione 1,00000). La
   pallina-fotone nasce con probabilità |a_i|² (più uno spostamento uniforme entro ds) e alla
   sorgente fa la moneta del flusso con Φ_in = k|U_src|², come `validate.py` (perdite 0,04–0,2%).
4. **Impulso e "chi sente chi".** Un ventaglio nasce quando il fronte sull'asse raggiunge l'oggetto.
   Come chiede il §8.3 conta solo che il *raggio* della pallina passi accanto al fotone; il raggio
   è "chiuso" solo se la sua pallina è morta contro un oggetto **prima del piano z del fotone**
   (gli oggetti sono piani ⊥ z; tolleranza 0,01 λ perché le palline muoiono proprio sull'oggetto
   dove il fotone legge Φ_in). Un primo tentativo con il criterio "raggio di morte ≥ R" sbagliava
   con ventagli radi (dθ = 0,1) e il test live = senza palline lo ha scoperto. In questi esperimenti
   nessun raggio viene mai chiuso prima del fotone, quindi la corsa senza palline (test, velocità
   ≥ ×50) è identica bit per bit a quella disegnata.
5. **Palline per origine** = palline per fotone del pannello / numero di origini (minimo 31 = dθ 0,1,
   variante H): la guida (E2) non dipende da dθ.
6. **Somme deterministiche**: la somma sulle origini è fatta a blocchi fissi di 32 (in parallelo con
   rayon nella corsa live) e sommata in ordine, quindi il risultato non dipende dal numero di thread
   (serve al tempo reversibile).
7. Restano valide le note del §7 (moneta al bordo come correzione di Kirchhoff, solo 2-D, Malus con
   moneta locale; l'entanglement resta in `entangle.rs`).
