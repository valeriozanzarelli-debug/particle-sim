# Coppie in scatola: se la scatola è troppo stretta, la coppia nasce?

La domanda di Valerio: *se il fotone gamma passa fra due lastre vicinissime, la coppia elettrone–positrone non si
forma, perché l'onda dell'elettrone non ha spazio?*

Risposta breve: **sì, in pratica sì**. Stringere la scatola rende l'elettrone **più pesante** (m_eff), e la
creazione di coppie dipende in modo esponenziale dal quadrato della massa: con la scatola larga meno di una
lunghezza Compton le coppie in un campo forte diminuiscono di **quasi mille volte** (e oltre: il calcolo esatto dice di
miliardi di volte; il reticolo arriva a misurare solo fino a 10⁻⁴). Per un fotone serve più energia: la soglia sale
da 2m a 2m_eff. La scala a cui succede è proprio la lunghezza Compton (d ≈ 2,2/m). Però la formula del compito
(m_eff con π/d) **non** è quella giusta per un elettrone: la parete gli lascia un po' più di spazio (vedi A1b).

Codice: `lab/box_core.py` (motore), `lab/box_run.py` (A1, A3), `lab/box_threshold.py` (A2),
`lab/box_formation.py` (A4), `lab/box_posthoc.py` (controllo a posteriori, non valutato), figura
`lab/box_figure.py` → `lab/out/box_figure.png`. Pre-registrazione: `lab/out/box_prereg.json` (commit 64d0afd) e un
emendamento dichiarato prima delle corse valutate, `lab/out/box_amend.json` (commit dc120f1). Risultati:
`lab/out/box_run.json`, `box_threshold_d*.json`, `box_formation.json`, `box_posthoc.json`.

---

## 1. La regola

È la stessa regola zig-zag di `lab/redo_r1.py` (quella che rifà Schwinger entro il 2–4 % e Sauter entro l'1 %,
`docs/IDEE_RIPRESE.md` R1), con **una direzione in più**, y, larga d palline:

* ogni tick la pallina fa il passo in x (R a destra, L a sinistra, con la fase del campo elettrico lungo x);
* poi fa un passo in y, su o giù secondo il suo "spin lungo y" (in mezzo alla striscia è l'onda di Dirac in 2+1 D,
  stessa velocità della luce nelle due direzioni);
* poi la moneta della R1 (cos θ, −i sin θ): θ = 0,1 è la massa.

Le pareti. All'inizio (pre-registrato) le pareti erano "specchi": la pallina che esce torna indietro. **Non
funziona**: quello specchio porta con sé uno *stato di bordo* (un'onda che corre lungo la parete senza massa,
come negli isolanti topologici), che chiude il buco fra mare e particelle; il conteggio delle coppie veniva
sbagliato del 10–35 % (`A1a_old_reflecting_wall` nel json). L'ho dichiarato e cambiato **prima** delle corse
valutate: la parete diventa una "parete di massa", 12 palline fuori dalla striscia con moneta θ_w = 1,5 (massa
grandissima, con lo stesso segno: è la parete "a massa infinita" dei libri, senza stati di bordo). Controllato: il
buco minimo resta a kx = 0, nessuna banda lo attraversa.

Unità: ħ = c = passo = tick = carica = 1; lunghezza Compton 1/m = 10 palline.

---

## 2. Risultati, voce per voce

| voce | che cosa | esito |
|---|---|---|
| A1a | coppie nel campo uniforme = Σⱼ exp(−π mⱼ²/E), mⱼ masse misurate dei modi della striscia, entro 10 % per ogni d, x | **FAIL** (9 casi su 27 fuori; tutti i d ≥ 24 dentro) |
| A1b | m_eff = √(m² + (π/d)²) entro 10 % | **FAIL** (atteso e dichiarato) |
| A1c | a campo fisso le coppie calano stringendo, fattore ≥ 10 a d ≤ 8 | **FAIL** sul "calano sempre" (fondo di misura), PASS sul fattore 10 |
| A2 | soglia del fotone = 2 m_eff entro 10 %, e sale stringendo | **FAIL** (d = 128: il massimo cade sul 2° modo); d = 8, 16, 32 entro 3–6 % |
| A3 | la scala di passaggio è la lunghezza Compton: 0,5 ≤ d* m ≤ 5 | **PASS** (d* m = 2,25) |
| A4 | lunghezza di formazione: (a) → 1 entro 15 % per ℓ ≥ 8 ℓ_f; (b) ≤ 0,5 per ℓ ≤ ℓ_f | **(a) FAIL** (errore di pre-registrazione: il gradino esce dalla banda), **(b) PASS** |

### A1 — Schwinger nella striscia

Campo uniforme E lungo x, tutto il mare della striscia evoluto esattamente (kx da −1 a +1). Il rapporto
reticolo / previsione:

| d (d·m) | x = E/m_eff² = 1 | 2 | 4 |
|---|---|---|---|
| 4 (0,4) | 0,84 | 0,55 | 1,02 |
| 6 (0,6) | 0,70 | 0,86 | 0,88 |
| 8 (0,8) | 0,72 | 0,94 | 0,87 |
| 12 (1,2) | 0,88 | 0,95 | 1,00 |
| 16 (1,6) | 0,89 | 0,95 | 0,95 |
| 24 (2,4) | 0,93 | 0,97 | 1,00 |
| 32 (3,2) | 0,97 | 0,98 | 0,97 |
| 64 (6,4) | 0,99 | 1,00 | 0,99 |
| 128 (12,8) | 0,99 | 0,99 | 1,01 |

Quando la scatola è larga la legge exp(−π m_eff²/E) è giusta all'1 %. Quando è stretta il numero ballava: il
controllo a posteriori (`box_posthoc.py`, **non valutato**) rifà lo stesso conto con finestre kx = ±0,6, ±0,8, ±1,3
e trova a d ≤ 16 valori fra 0,3 e 1,5 per lo stesso caso: il conteggio "tutto il mare" in una finestra finita è
sporcato dagli stati che entrano ed escono dalla finestra (le bande della parete stanno a |w| ≥ 1,5, i modi trasversi
alti ci arrivano vicino). Non è la fisica che sbaglia, è il mio modo di contare; ma il criterio era quello, quindi
FAIL. Con la scatola stretta m_eff arriva a 0,3–0,4 (non più "piccola" per il reticolo), e anche questo pesa.

### A1b — quanta massa dà la scatola

| d | d·m | m_eff misurata / m | con π/d (compito) | con π/2d | Dirac, pareti a massa infinita |
|---|---|---|---|---|---|
| 4 | 0,4 | 3,78 | 7,92 | 4,05 | 4,59 |
| 8 | 0,8 | 2,42 | 4,05 | 2,20 | 2,65 |
| 16 | 1,6 | 1,64 | 2,20 | 1,40 | 1,70 |
| 32 | 3,2 | 1,25 | 1,40 | 1,11 | 1,27 |
| 64 | 6,4 | 1,085 | 1,114 | 1,030 | 1,088 |
| 128 | 12,8 | 1,025 | 1,030 | 1,008 | 1,026 |

La formula "π/d" è quella di una particella **senza spin** in una scatola con le pareti dure (l'onda deve stare a
zero sulle due pareti: mezza lunghezza d'onda nella scatola). L'elettrone è un'onda di Dirac: la parete non lo
obbliga a zero, gli chiede solo che una certa combinazione delle sue due componenti si annulli. Così l'onda "sporge"
un po' e l'impulso trasverso è più piccolo: per la parete a massa infinita la regola giusta è
p·d = π − arctan(p/m), che va da π/(2d) (scatola strettissima) a π/d (scatola larga). Le nostre palline seguono
questa curva entro l'8 % (d = 8), il 3 % (d = 16), l'1 % (d ≥ 32): la parete vera (12 palline, θ_w = 1,5) lascia
penetrare l'onda un poco, quindi la massa viene un filo più bassa. La formula del compito sbaglia fino a un fattore 2
quando la scatola è stretta: FAIL, come dichiarato prima.

### A1c — a campo fisso

E = 2 E_crit (senza scatola: 0,208 coppie per ogni modo kx che passa). Con la scatola:

| d·m | 12,8 | 6,4 | 3,2 | 2,4 | 1,6 | 1,2 | 0,8 | 0,6 | 0,4 |
|---|---|---|---|---|---|---|---|---|---|
| coppie (reticolo) | 0,54 | 0,24 | 0,088 | 0,049 | 0,013 | 0,0025 | 2,7·10⁻⁴ | 3,0·10⁻⁴ | 1,9·10⁻⁴ |
| previsione Σ exp(−π mⱼ²/E) | 0,54 | 0,24 | 0,089 | 0,051 | 0,014 | 0,0033 | 1,0·10⁻⁴ | 1,8·10⁻⁶ | 1,8·10⁻¹⁰ |

(Con la scatola larga i numeri superano 0,208 perché ci sono **più modi trasversi**: la striscia larga è tanti
canali 1+1 D in parallelo.) Sotto d·m ≈ 1 il reticolo misura un **fondo** di ~2·10⁻⁴ (la coda non adiabatica agli
estremi della finestra kx), e lì la curva non scende più: 3,0·10⁻⁴ a d = 6 è più di 2,7·10⁻⁴ a d = 8. Quindi
"calano sempre" FAIL per colpa del fondo; il fattore di soppressione rispetto a 1+1 D è ≥ 690 a d ≤ 8 (PASS).
La previsione dice che il vero calo è molto più grande (10⁻¹⁰ a d·m = 0,4).

### A2 — la soglia del fotone

Un "fotone" = un campo debole che oscilla, A(t) = 0,01 sin(Ωt), acceso per 600 tick; misuro quanto mare viene
eccitato per ogni frequenza Ω (spettro di assorbimento, figura in basso a sinistra).

| d | 2 m_eff | Ω* (massimo) | Ω*/2m_eff | inizio (10 % del massimo) |
|---|---|---|---|---|
| 8 | 0,484 | 0,498 | 1,03 | 0,450 |
| 16 | 0,329 | 0,343 | 1,04 | 0,293 |
| 32 | 0,251 | 0,265 | 1,06 | 0,215 |
| 128 | 0,205 | 0,258 | **1,26** | 0,183 |

Senza scatola la soglia è 2m = 0,2. Con la scatola **sale**: 0,25 (d = 32), 0,33 (d = 16), 0,48 (d = 8), sempre
a 2 m_eff. Un fotone con energia fra 2m e 2m_eff passa fra le lastre **senza poter fare la coppia**: è proprio
l'idea di Valerio. A d = 128 (scatola larga 13 lunghezze Compton) i modi trasversi sono fitti e il massimo dello
spettro cade sul secondo-terzo modo, non sul primo: con la definizione pre-registrata ("il massimo") è FAIL.
L'inizio dell'assorbimento sta un po' sotto 2 m_eff per l'allargamento dovuto ai 600 tick (larghezza 2π/600 = 0,01).
Il passo della griglia in kx per d = 32 e 128 è 0,008 / 0,016 invece di 0,004 (per stare nei 6 minuti; dichiarato
nel file, effetto trascurabile rispetto alla larghezza 0,01).

### A3 — la scala

m_eff = √2·m (l'impulso trasverso uguale alla massa) a **d* = 22,5 palline = 2,25 lunghezze Compton**: PASS.
Sopra qualche lunghezza Compton la scatola non conta quasi nulla (d·m = 6,4: massa +8,5 %); sotto una lunghezza
Compton conta tutto.

### A4 — la lunghezza di formazione

Un campo che agisce solo su un tratto lungo ℓ (gradino di Sauter, campo massimo E = 2 E_crit fisso). La coppia
deve prendere dal campo l'energia 2m: le serve almeno un tratto ℓ_f = 2m/E = 10 palline. Rapporto fra le coppie vere
e la stima "locale" (come se ogni pezzo di campo facesse coppie per conto suo):

| ℓ/ℓ_f | 0,2–0,5 | 0,7 | 1 | 1,4 | 2 | 2,8 | 4 | 5,6 | 8 | 11 | 16 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| senza scatola | 0 | 0,31 | 0,64 | 0,81 | 0,90 | 0,95 | 0,97 | 0,985 | 0,991 | 0,60* | 0* |
| striscia d = 16 (m → m_eff) | 0 | 0 | 0,11 | 0,43 | 0,68 | 0,82 | 0,90 | 0,94 | 0,96 | 0,71* | 0* |

\* il gradino è più alto di π: nel reticolo l'energia è un angolo, un gradino più alto di mezzo giro "si
riavvolge" e il conto non ha senso. È un errore della mia pre-registrazione (ℓ fino a 160 con E fisso), non della
fisica: (a) FAIL. Fino a ℓ = 8 ℓ_f (gradino 3,2, ancora dentro la banda) la stima locale torna, 0,99.
(b) PASS: un campo più corto della lunghezza di formazione **non fa coppie** (zero esatto sotto ℓ_f/2, dove il
gradino è più basso di 2m), e a ℓ = ℓ_f ne fa il 64 %. È l'analogo dell'effetto LPM (Landau–Pomeranchuk–Migdal):
un processo ha bisogno di un tratto minimo per "formarsi", e se il campo è più corto viene soppresso. Con la
scatola la lunghezza di formazione cresce (2 m_eff/E), la curva arancione si sposta a destra: un'altra strada per
cui le lastre vicine ostacolano le coppie.

---

## 3. Che cosa resta vero e che cosa no

* **Vero (forte):** stringere la scatola alza la massa efficace, e le coppie seguono exp(−π m_eff²/E) all'1 %
  quando il conteggio è pulito (d ≥ 24). La soglia per un fotone sale a 2 m_eff. La scala è la lunghezza Compton.
* **Non vero:** la formula m_eff = √(m² + (π/d)²) per l'elettrone. Quella è per una particella senza spin. Per
  l'elettrone (Dirac) l'impulso trasverso minimo va da π/(2d) a π/d; con pareti a massa infinita
  p·d = π − arctan(p/m).
* **Debole:** con la scatola strettissima (d·m < 1,5) il mio conteggio ha errori del 10–50 % (finestra, fondo
  2·10⁻⁴); il PASS/FAIL A1a è FAIL lì.
* **Onesto sulle pareti:** la parete pre-registrata ("specchio") aveva uno stato di bordo e l'ho cambiata prima delle
  corse valutate (emendamento nel json, con i numeri di prova già visti). Le pareti vere di un condensatore non sono
  né l'una né l'altra: questo è un modello.
* Sulla domanda fisica vera: per un fotone gamma da 1 MeV la lunghezza Compton dell'elettrone è 0,4 pm. Lastre a
  quella distanza non esistono (sono più piccole di un atomo): nella realtà l'effetto è irrilevante, ma la logica di
  Valerio è giusta.
