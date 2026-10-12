# Nodo 3D: una particella come nodo del fluido di palline (modello di Skyrme su reticolo)

**L'idea (di Valerio).** Le palline sono i "pixel" della materia: tantissime, legate forte fra vicine, e **sono** il
vuoto. Una particella non è una pallina in più: è un **nodo** nel modo in cui le palline sono orientate, un nodo che
non si può sciogliere senza strappare il fluido. Per questo dovrebbe avere **una sola massa**, sempre la stessa.

**Come l'abbiamo provato.** Ogni pallina porta una piccola "freccia" a 4 componenti di lunghezza 1
(φ = (σ, π₁, π₂, π₃), cioè un'orientazione SU(2)). Le vicine vogliono essere allineate (accoppiamento forte), e c'è un
secondo termine, quartico, che impedisce al nodo di stringersi fino a sparire (teorema di Derrick). È esattamente il
**modello di Skyrme** (1961) messo su un reticolo cubico di palline. Il nodo più semplice è il "riccio" (hedgehog):
lontano tutte le frecce puntano in su (il vuoto), al centro puntano in giù, e in mezzo girano a seconda della
direzione. Il numero di avvolgimenti B (numero barionico) è un intero e non può cambiare in modo continuo.

> **Attenzione, da dire chiaro.** Gli skyrmioni descrivono i **barioni** (protone, neutrone, nuclei), non gli
> elettroni. Il campo di Skyrme è il campo dei **pioni**; la sua "velocità della luce" c = 1 è la velocità dei pioni
> senza massa. L'elettrone (leggero, puntiforme fino a 10⁻¹⁸ m, carica senza struttura forte) **non** è un
> skyrmione, e qui non c'è nessun elettrone.

Pre-registrazione: `lab/knot_prereg.py` → `lab/out/knot_prereg.json`, **committata (d1bebbf) prima delle corse di
produzione**. Le prove di sviluppo (dev) viste prima sono dichiarate nel suo docstring.
Codice: `lab/knot_core.py` (motore), `lab/knot_k1.py` … `lab/knot_k5.py`, `lab/knot_summary.py`, `lab/knot_figure.py`.
Numeri: `lab/out/knot_*.json`, log in `lab/out/knot_logs/`. Figura: `lab/out/knot_figure.png`.

## Il modello in una riga, e quanti parametri

Unità di Skyrme (energia F_π/4e, lunghezza 2/(e F_π), c = 1):

* energia statica per pallina: |∂φ|² + ½[(tr D)² − tr D²], con D_ij = ∂_iφ·∂_jφ;
* energia cinetica: |φ̇|²(1 + |∂φ|²) − Σ_i (φ̇·∂_iφ)²  → il modello è invariante di Lorentz con c = 1;
* limite inferiore (Faddeev–Bogomolny): E ≥ 12π²|B|.

**Parametri.** La lagrangiana ne ha **2** (F_π ed e); la massa del pione l'abbiamo messa a zero (sarebbe il 3°).
In unità di Skyrme **non resta nessun parametro libero**: forma del nodo, E/12π², rapporti, contrazione di Lorentz
sono numeri puri (0 parametri). Per passare a MeV e fm servono 2 misure (noi, come ANW: M_N e M_Δ). Le scelte
numeriche (passo h, scatola, dt, spugna) non sono fisica: si controllano facendole variare.

Reticolo: differenze in avanti/indietro simmetrizzate (il termine di Skyrme è la media delle 4 combinazioni
avanti/indietro). Con le differenze centrali del 4° ordine il riccio si **srotolava** del tutto (modo a scacchiera a
energia zero, lo stesso guaio già visto in 2D in `LUCE_INTRAPPOLATA.md`): dichiarato nel dev. Bordi: φ = vuoto (Dirichlet).
Dinamica: forma hamiltoniana (π = 2Mφ̇, M inversa calcolata localmente), RK4 + proiezione sul vincolo |φ| = 1 a ogni passo.

## Tabellone

| prova | cosa | esito pre-registrato |
|---|---|---|
| K1a | E/12π² (Richardson in h², metà scatola 4.7) entro 1 % da 1.2315 | **PASS** (1.2410, +0.77 %) |
| K1b | stessa massa al cambiare della scatola (L = 48 → 64) entro 0.5 % | **PASS** (0.26 %) |
| K1c | la massa converge come potenza 1/L³, non come ln L | **PASS** (rapporto 3.47; ln L darebbe 1.4) |
| K1d | B = 1 entro 1 %; partenza "sbagliata" → stessa E entro 0.2 % | **PASS** (|B−1| ≤ 0.005; 6·10⁻¹²) |
| K2a | E = γE₀ entro 3 % (v = 0.3, 0.6, 0.8) | **PASS** (+0.2 %, +1.0 %, +2.8 %) |
| K2b | contrazione di Lorentz σx/σ⊥ = 1/γ entro 3 % | v = 0.3 **PASS** (−0.9 %), 0.6 PASS (+1.8 %, ma vedi sotto), 0.8 **FAIL** (+47 %) |
| K2c | il nodo contratto mantiene la forma (3 %) | v = 0.3 **PASS** (−0.7 %), 0.6 **FAIL** (+14 %), 0.8 **FAIL** (×3.6) |
| K2d | velocità del centro = v entro 3 % | v = 0.3 **PASS** (0.299), 0.6 **FAIL** (0.48), 0.8 **FAIL** (0.51) |
| K2e | quante palline per raggio servono | **3.8** passa (h = 0.2); 3.2 si rompe; a h = 0.15, 4.4 palline: E e contrazione ok, velocità −5 % |
| K2 (h = 0.15, v = 0.8) | ripetizione pre-registrata sul reticolo fine | **FAIL** ma meglio (velocità 0.62, contrazione +24 %) |
| K3a | nodo + antinodo: ∫\|b\| scende sotto 0.2, a ogni velocità | **FAIL** (v = 0.3, 0.6: sì, 0.06–0.07; v = 0: mai) |
| K3b | ≥ 70 % dell'energia esce come onde di pioni dalla sfera r < 4.5 entro t = 12 | **FAIL** (57 %, 62 %; a v = 0: 1.6 %) |
| K3c | nodo + nodo: B = 2 resta, toro B = 2 o urto a 90° | **PASS** (B = 1.993; esce a 90°) |
| controllo | nodo da solo: B costante entro 1 % | **PASS** (0.9959 → 0.9948) |
| K4a | momento d'inerzia Λ del reticolo entro 2 % dall'ODE | **FAIL** (−23 %: coda 1/R della scatola) |
| K4b | \|μp/μn\| entro 3 % da 1.460 (e fissato da M_Δ/M_N) | reticolo **FAIL** (1.548, +6 %); ODE esatta **PASS** (1.417, −2.9 %) |
| K4c | fallimenti attesi: g_A e F_π entro 30 % | **FAIL come previsto** (g_A 0.61 contro 1.27; F_π 129 contro 184 MeV) |
| K5a | disegno in negativo: triangoli / palline < 1 % | **FAIL** (1.19 %: 3120 triangoli per 262 144 palline) |
| K5b | il disegno ritrova il nodo (centro < h, anisotropia < 5 %) | **PASS** (centro 2·10⁻⁶ h, anisotropia 2.1 %) |

---

## K1 — Una sola massa, che non dipende dalla scatola

**Riferimento esatto.** Il riccio a simmetria sferica, φ = (cos f(r), sin f(r) r̂), risolto come ODE radiale.
Lo abbiamo rifatto con collocazione (`radial_bvp`, condizione esatta della coda f = C/r²): **E/12π² = 1.231446**,
E₂ = E₄ = 72.9233 (Derrick, uguali a 10⁻⁹), coda **C = 2.160**. Il valore di letteratura per pioni senza massa
è 1.232 (ANW: M = 36.5 F_π/e = 1.2315·12π²): **verificato**.

**Reticolo.** Riccio rilassato (minimo dell'energia, B = 1) su cubi L³:

| h (passo) | L | metà scatola | E/12π² | B |
|---|---|---|---|---|
| 0.2  | 32 | 3.1 | 1.21645 | 0.9950 |
| 0.2  | 48 | 4.7 | 1.20578 | 0.9956 |
| 0.2  | 64 | 6.3 | 1.20271 | 0.9958 |
| 0.15 | 32 | 2.3 | 1.25434 | 0.9984 |
| 0.15 | 48 | 3.5 | 1.22898 | 0.9987 |
| 0.15 | 64 | 4.7 | 1.22119 | 0.9988 |
| 0.25 (info) | 38 | 4.6 | 1.18365 | 0.9857 |

* Il reticolo **sottostima** l'energia di qualche % (O(h²): −2 % a h = 0.2, −0.8 % a h = 0.15).
  Estrapolando in h² (stessa scatola, 4.7): **1.2410**. È più del valore infinito (1.2315) di +0.77 % (PASS), e
  quasi tutto quello che resta è l'effetto della scatola: lo stesso riccio chiuso in una **palla** di raggio 4.8
  con f = 0 sul bordo vale 1.2442 (ODE). Cioè: reticolo e continuo concordano entro ~0.3 % a parità di scatola.
* Ingrandendo la scatola la massa si **ferma**: da L = 48 a 64 cambia di 0.26 %. Le differenze successive stanno
  nel rapporto 3.47 (h = 0.2) e 3.25 (h = 0.15): è una convergenza a potenza (coda ∝ 1/R³, perché f ∝ 1/r²), **non**
  la crescita logaritmica del vortice 2D di `FLUIDO_IN_NEGATIVO.md` (F2: E ∝ ln R). In 3D il nodo ha **una massa vera**.
* Partendo da un riccio "sbagliato" (f = π e^(−r/0.6), troppo piccolo) si arriva alla **stessa** energia fino a
  6·10⁻¹²: il nodo trova da solo la sua forma.
* Quante palline servono: con h = 0.2 il raggio a metà (f = π/2) è ~4 palline. Nel dev, a h ≥ 0.3 (< 3 palline
  per raggio) il nodo si **srotolava** passando fra le palline (B → 0): con troppo poche palline la topologia non
  protegge niente. A h = 0.25 (3.5 palline) in produzione regge (B = 0.986).

**Deviazione dichiarata.** La pre-registrazione diceva "arrested Newton flow, dt = 0.15 h". Al primo tentativo
(h = 0.2, L = 48) è arrivato a 1.205782 e poi un modo rigido del reticolo è esploso. Abbiamo usato L-BFGS sulla
stessa energia (φ = ψ/|ψ|); il minimo raggiunto è lo stesso (1.2057818).

## K2 — Lorentz in 3D: il nodo che corre si schiaccia?

Prendiamo il nodo rilassato del reticolo (h = 0.2), lo schiacciamo lungo x di 1/γ, gli diamo la velocità giusta
e lo lasciamo correre con la dinamica completa del 2° ordine (scatola 112×48×48 palline, energia conservata a
10⁻⁵–2·10⁻³). Misuriamo l'energia, la forma della densità barionica e la velocità del centro.

| v | γ | raggio contratto / h | E/(γE₀) | σx/σ⊥ ÷ (1/γ) | forma fine/inizio | velocità misurata |
|---|---|---|---|---|---|---|
| 0 | 1 | 4.0 | 1.000 | 1.000 (riposo) | 0.994 | 0.005 |
| 0.3 | 1.048 | 3.8 | **1.002** | **0.991** | **0.993** | **0.299** |
| 0.6 | 1.25 | 3.2 | 1.010 | 1.018 | 1.14 | 0.48 (rallenta: 0.58 → 0.46) |
| 0.8 | 1.667 | 2.4 | 1.028 | 1.47 | 3.6 | 0.51 (0.68 → 0.49) |
| 0.8, h = 0.15 | 1.667 | 3.3 | 1.016 | 1.24 | 2.9 | 0.62 (0.73 → 0.61) |
| 0.6, h = 0.15 (extra) | 1.25 | 4.4 | **1.004** | **0.994** | 0.96 | 0.568 (quasi costante, 0.587 → 0.565) |

* A **v = 0.3 tutto torna**: E = γE₀ entro 0.2 %, contrazione entro 1 %, il nodo contratto è una soluzione che
  viaggia senza deformarsi, a velocità 0.299. Il fluido di palline "sa" la relatività, come in 1D (FK, F3).
* A v = 0.6 e 0.8 il nodo **frena** e si **gonfia**: perde impulso nel reticolo (attrito di Peierls–Nabarro e
  radiazione), mentre l'energia totale resta conservata. A v = 0.6 il rapporto σx/σ⊥ passa il test solo per caso
  (il nodo rallentato oscilla attorno alla forma giusta per la velocità iniziale); forma e velocità falliscono.
* **Quante palline servono (K2e):** a h = 0.2 passa tutto solo con **3.8 palline per raggio** (contratto); a 3.2 si
  rompe. Un nodo veloce si schiaccia, quindi a parità di reticolo **esiste una velocità massima "di reticolo"** ben
  sotto c. Il reticolo più fine (h = 0.15) lo conferma: a v = 0.8 la velocità sale da 0.51 a 0.62 e il nodo si gonfia
  meno; a v = 0.6 (4.4 palline per raggio, corsa in più non pre-registrata) energia e contrazione tornano entro 1 %,
  B resta 0.998, e il nodo corre quasi senza frenare, ma a **0.568 invece di 0.6** (−5 %) e si allarga del 4 %.
  Quindi anche con 4.4 palline per raggio la velocità non è ancora giusta entro 3 %: in 3D serve più di quanto
  bastava in 1D (FK, F3: ~4 palline davano 1.5 %).
* **Possibile causa in più (debito).** Nel nostro reticolo il termine cinetico usa le differenze centrali e quello
  statico la media avanti/indietro: le due "velocità della luce" del reticolo coincidono solo per onde lunghe.
  Un nodo schiacciato contiene onde corte, che vanno più piano. Una discretizzazione con cinetico e statico costruiti
  allo stesso modo potrebbe migliorare; non l'abbiamo provata.

**Deviazione dichiarata.** dt = 0.25 h invece di 0.1 h (costo). Test di sviluppo a v = 0.8: deriva dell'energia
2·10⁻⁵ per unità di tempo; in produzione la deriva massima è 0.19 % < 0.5 % (corse valide).

## K3 — Nodo + antinodo → onde

**Esperimento.** Su un cubo di 64³ palline (h = 0.2, metà lato 6.3, una "spugna" che assorbe le onde negli ultimi
1.4) mettiamo un nodo (B = +1) e un antinodo (B = −1, il riccio con le frecce rovesciate) a distanza 4 (prodotto
U₁U₂), e li lanciamo uno contro l'altro. Misuriamo ∫|b| (quanto "nodo" c'è, positivo o negativo) e l'energia rimasta
dentro la sfera r < 4.5. Il numero totale B = 0 resta 0 a 10⁻³ per tutto il tempo.

| corsa | v (ciascuno) | ∫\|b\| min | quando sotto 0.2 | energia uscita da r < 4.5 a t = 12 | esito |
|---|---|---|---|---|---|
| nodo + antinodo | 0   | 1.983 | mai | 1.6 % | K3a **FAIL**, K3b **FAIL** |
| nodo + antinodo | 0.3 | 0.058 | t = 7.75 | **57 %** | K3a PASS, K3b **FAIL** |
| nodo + antinodo | 0.6 | 0.074 | t = 7.25 | **62 %** | K3a PASS, K3b **FAIL** |
| nodo da solo (controllo) | 0 | 0.992 | — | 2.3 % | B costante: **PASS** |
| nodo + nodo (canale attrattivo) | 0.5 | B = 1.993 sempre | — | 15 % | K3c **PASS**: urto a 90° |

* **Il nodo si scioglie davvero.** A v = 0.3 e 0.6 i due nodi si toccano (t ≈ 3), restano un po' attaccati, e fra
  t = 5 e t = 8 ∫|b| crolla da 2 a 0.06: il nodo e l'antinodo si **annullano**, non rimbalzano. Questa è la
  differenza vera con il kink 1D di `FLUIDO_IN_NEGATIVO.md` (F3), dove kink + antikink non si annichilavano.
* **Diventano onde, ma non abbastanza in fretta per il nostro esame.** Dopo l'annichilazione un guscio di onde di pioni
  esce a velocità c: l'energia dentro la sfera scende di 40–90 unità per unità di tempo proprio alla fine
  (t = 10 → 12: 257 → 129 a v = 0.3; 233 → 139 a v = 0.6). Al tempo pre-registrato t = 12 ne è uscito il 57 % e il
  62 %, sotto il 70 %: **FAIL**. Andando avanti di ~1 unità di tempo si passerebbe il 70 %, ma questo è dopo i
  fatti e non conta. Lo sbaglio di progetto: abbiamo stimato male quanto dura lo "scioglimento" (3–4 unità di tempo).
* **A riposo non succede niente.** Nell'orientazione scelta (antinodo = riccio con le frecce rovesciate, nessuna
  rotazione relativa) nodo e antinodo fermi **non si attirano**, anzi si allontanano un poco (centri da 2.00 a 2.13
  in 12 unità). Topologicamente potrebbero annullarsi (B = 0), ma c'è una barriera di energia. L'annichilazione
  richiede velocità (o un'altra orientazione relativa, che non abbiamo provato).
* **Nodo + nodo (B = 2).** Nel canale attrattivo (secondo nodo girato di 180° nell'isospin attorno a z) a v = 0.5:
  i due nodi arrivano lungo x, a t ≈ 4 formano per un attimo una configurazione a **anello** (densità al centro
  0.2 del massimo, il buco del toro B = 2), poi escono **lungo y, a 90°**, a velocità ~0.3 (il resto dell'energia
  va in onde). B resta 1.993 (0.4 %). È proprio lo "scattering a 90°" della letteratura (Battye–Sutcliffe; Manton–
  Sutcliffe, *Topological Solitons*, cap. 9): non rimbalzano e non si attraversano, girano la direzione di 90°.
  Il toro **legato** non si forma con questa velocità (troppa energia cinetica); non l'abbiamo cercato rilassando.

**Deviazione dichiarata.** dt = 0.25 h invece di 0.1 h; il controllo col nodo da solo tiene B entro 0.1 %.

## K4 — Spin ½ e numeri senza parametri (Adkins–Nappi–Witten 1983)

**Spin ½: citato, non derivato.** Il nodo B = 1 è un oggetto rigido che può ruotare. Che i suoi stati quantistici
abbiano spin ½ (fermione) viene da Finkelstein–Rubinstein (1968): lo spazio delle configurazioni ha un anello
(rotazione di 2π) che non si può contrarre (π₁ = Z₂), e si può scegliere la quantizzazione in cui quella rotazione
cambia segno alla funzione d'onda; Witten (1983) ha mostrato che con N_c = 3 colori la scelta è proprio quella.
**Noi non l'abbiamo dimostrato numericamente**: lo citiamo.

**Momento d'inerzia.** Ruotando rigidamente il nodo nello spazio dell'isospin con velocità angolare ω,
T = ½Λω². ODE esatta: **Λ = 106.77** in unità di Skyrme (= (2π/3)·50.98 nelle unità di ANW, che riportano 50.9:
verificato). *Correzione dichiarata:* la pre-registrazione scriveva 103.90, calcolato con un raccordo della coda
troppo presto (C = 1.99 invece di 2.16): era sbagliato del 2.7 %.

**K4a sul reticolo: FAIL.** Λ(h=0.2) = 63.3 / 72.3 / 77.1 per metà scatola 3.1 / 4.7 / 6.3; Λ(h=0.15) = 58.0 / 69.5 / 76.5.
Richardson (scatola 4.7): **81.9**, cioè −23 %. Il motivo: Λ pesa r² sin²f ≈ C²/r², una coda che converge solo come
**1/R** (non 1/R³ come la massa). Post-hoc (non pre-registrato): estrapolando Λ(R) = Λ∞ − c/R a ogni h e poi in h²
si ottiene 104.0 (−2.6 %). Quindi **la massa è locale, l'inerzia di rotazione no**: la "nuvola di pioni" lontana
conta molto per come il nodo gira.

**Numeri di ANW (rifatti qui).** Due ingressi, M_N = 938.9 MeV e M_Δ = 1232 MeV (rotore quantistico:
M = M₀ + J(J+1)/2Λ). Con l'ODE esatta:

| grandezza | tipo | noi (ODE) | ANW 1983 | misura |
|---|---|---|---|---|
| e | **fit** (da M_Δ/M_N) | 5.44 | 5.45 | — |
| F_π | **fit** (da M_N, dato e) | 129 MeV | 129 MeV | 184 MeV (−30 %) |
| raggio isoscalare ⟨r²⟩^½ | previsione (1 ingresso) | 0.59 fm | 0.59 fm | 0.72 fm |
| μ_p | previsione | 1.88 | 1.87 | 2.79 |
| μ_n | previsione | −1.33 | −1.31 | −1.91 |
| **\|μp/μn\|** | previsione, **dipende solo da e** | **1.417** | 1.43 | 1.460 (**−2.9 %**) |
| g_A | previsione | 0.61 | 0.61 | 1.27 (−52 %) |
| μ_NΔ | previsione | 2.27 | 2.3 | 3.3 |

Numeri **a zero ingressi** (puri, in unità di Skyrme): E/12π² = 1.2314; Λ = 106.8; r_E/r_B = 1.377; rapporti del
rotore (M_{5/2}−M_{3/2})/(M_{3/2}−M_{1/2}) = 5/3; μ_NΔ/(μp−μn) = 1/√2 (identità della quantizzazione collettiva).
I rapporti come μp/μn non dipendono da F_π ma solo da e, quindi sono "un ingresso" (M_Δ/M_N), non zero.

* K4b col Λ del **reticolo** (pre-registrato): \|μp/μn\| = 1.548 (+6 %) **FAIL**, perché il Λ del reticolo è troppo
  piccolo (scatola). Con l'ODE esatta: 1.417 (−2.9 %) **PASS**, appena.
* K4c: come la letteratura aveva detto, g_A (0.61) e F_π (129 MeV) sono sbagliati di −52 % e −30 % (**FAIL atteso**).
  Il modello di Skyrme è un modello "al 30 %" dei barioni.

## K5 — Disegnare in negativo: solo il nodo

Disegniamo solo la superficie dove la densità barionica vale il 10 % del massimo (marching tetrahedra), niente palline.
Nodo rilassato h = 0.2, cubo 64³ = 262 144 palline: **3120 triangoli** (1.19 %, **FAIL** sul limite 1 %, per poco).
Il disegno è fedele: centro a 2·10⁻⁶ h dal vero, sfera con anisotropia 2.1 % (il taglio dei cubi in tetraedri), e
la superficie contiene B = 0.72. Il numero di triangoli non dipende dalla scatola ma dal passo (area/h²): con L = 96
lo stesso nodo starebbe sotto 0.4 %. Nella figura (pannello e) c'è il nodo disegnato così.

## Debiti e cose che mancano

* **Elettroni e leptoni: niente.** Lo skyrmione è un barione (raggio ~0.6 fm, massa ~1 GeV, struttura). Un elettrone
  non è un nodo di pioni. Serve un'altra idea (o un altro campo) per i leptoni.
* **Carica elettrica: non c'è.** Q = B/2 + I₃ appare solo dopo la quantizzazione collettiva (la usiamo per μp, μn),
  ma nel reticolo non c'è nessun campo elettromagnetico, nessun fotone, nessuna forza di Coulomb.
* **c dei fononi = c della luce: non dimostrato.** Qui c = 1 è la velocità delle onde del campo di palline (i pioni
  senza massa), messa nella lagrangiana. Che sia la stessa c della luce è un'ipotesi, non un risultato.
* **Spin ½**: citato (Finkelstein–Rubinstein, Witten), non dimostrato col reticolo.
* **Quantizzazione semiclassica**: rotore rigido; ANW stessi la considerano una stima al 30 %.
* **Massa del pione = 0**: con m_π ≠ 0 (3° parametro) la coda diventa esponenziale e K4a sarebbe molto più facile.
* **Reticolo**: servono ≥ 4 palline per raggio; a velocità alte il nodo si schiaccia e il reticolo lo frena.
  "Le palline sono i pixel" ha un prezzo: un pixel troppo grosso rompe la relatività e, sotto ~3 palline, anche la topologia.
* **Parametri:** 2 nella lagrangiana (F_π, e), 0 per forme e rapporti in unità di Skyrme, 1 ingresso (M_Δ/M_N) per
  μp/μn, 2 ingressi per i numeri in MeV e fm.
