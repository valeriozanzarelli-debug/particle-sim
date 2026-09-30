# Un fotone fatto di palline tutte uguali: cosa funziona e cosa no

Laboratorio: `lab/identical_*.py` (Python, numpy/scipy, seed fissi). Ogni numero qui sotto viene da una
corsa fatta davvero; il comando è indicato accanto. Riferimenti e punteggi (flusso esatto di
Rayleigh–Sommerfeld sullo schermo, correlazione "grossolana", errore sulle altezze delle frange)
sono quelli di `lab/optics.py`, come in `docs/REGOLA_FOTONE.md`.

---

## 1. La risposta in breve

* **Il requisito:** un fotone = N palline **identiche**. Ogni pallina segue **la stessa regola**, guarda
  solo le palline vicine, e la **densità delle palline** deve essere |ψ|².
* **Fascio libero: funziona.** Due regole di questo tipo danno la larghezza giusta ben entro il 3%:
  * **regola A** (la forza "Many Interacting Worlds" di Hall–Deckert–Wiseman): errore −0,49% con 100
    palline, −0,05% con 1000;
  * **regola C** (ogni pallina porta anche una fase): la dinamica è esatta (errore della regola ≤ 0,005%);
    l'unico errore è quello di come sono state messe le palline all'inizio.
  * **Quante palline:** con palline messe in modo regolare ne bastano **50** (errore 1,3%). Con palline
    messe a caso l'errore è il rumore di partenza, circa 1/√(2N): servono **qualche centinaio**.
* **Doppia fenditura: non funziona.** Nessuna regola a palline identiche arriva vicino a 0,97.
  * Regola A: correlazione **0,42–0,60** per N = 50…1000. Nessuna frangia: le palline riempiono le
    zone buie (18% delle palline contro il 2% della teoria). Non si incrociano mai.
  * Regola C: si **rompe** nel momento in cui le due metà del fotone si incontrano (a z ≈ 4–8).
    Punteggio finale tra −0,01 e 0,45. Più palline **non** aiutano: l'errore al momento
    dell'incontro **cresce** con N (0,6% con 100 palline, 35% con 3000).
  * Un fotone "perfetto" di N palline (palline messe nei quantili della figura giusta) farebbe 0,96
    con 50 palline e 0,99 con 500. Quindi il problema non è il numero: è la regola.
* **Perché (in parole semplici):** le frange scure nascono dove l'onda di sinistra e quella di destra
  si cancellano. Per cancellarsi devono esserci entrambe. Ma le palline stanno dove c'è luce: tra le due
  fenditure, all'inizio, non c'è quasi nessuna pallina. L'onda che viene da destra dovrebbe arrivare alle
  palline di sinistra **attraversando una zona vuota**, e nessuno la porta. Quando i due gruppi si
  toccano, ognuno ha solo la propria metà d'onda: manca proprio il pezzo che fa le frange.
* **Ingrediente minimo che manca:** (1) ogni pallina deve portare una **fase** (senza fase è impossibile,
  §6.1); (2) non basta: servono palline che portino l'onda **anche dove non c'è luce**. Allora la densità
  delle palline non può essere |ψ|² dappertutto: ogni pallina deve portare anche un **peso** (quanto
  "fotone" porta). Le palline possono restare tutte uguali (stessa regola, numeri diversi). La regola già
  validata in `REGOLA_FOTONE.md` (0,995) è il caso estremo: peso 0 per le palline-onda, 1 per una sola.
* **Polarizzazione:** ogni pallina porta un vettore di Jones. Con un **beam-splitter polarizzatore** la
  legge di Malus esce **senza nessuna pallina speciale e senza monete**: la frazione di palline che finisce
  nel fascio "a" è 0,8540 / 0,5000 / 0,1460 contro cos²θ = 0,8536 / 0,5 / 0,1464. Il singolo clic però
  richiede di dire **quale** pallina è "la nostra" (§7).

---

## 2. Le regole provate: cosa fa UNA pallina a ogni passo

Modello: 2D parassiale. x = posizione trasversale, z = distanza lungo il fascio (fa da "tempo").
Il fotone in x si muove come una particella di massa k = 2π/λ con ħ = 1: i ∂ψ/∂z = −(1/2k) ∂²ψ/∂x².
Ogni pallina avanza in z alla velocità della luce; la regola decide quanto si sposta di lato.
λ = 1 in tutto il documento.

### Regola A — "Many Interacting Worlds" (HDW 2014) — la migliore tra quelle *senza fase*

A ogni passo, **una pallina**:
1. misura la distanza dalle 2 palline alla sua sinistra e dalle 2 alla sua destra (g = distanze tra
   vicine consecutive; oltre l'ultima pallina la distanza conta come infinita);
2. calcola la spinta quantistica F = −∂U/∂x della "energia" U = (1/8k) Σ_n (1/g_n − 1/g_{n−1})²
   (dipende solo da quelle 4 distanze: formula a 5 punti di Hall–Deckert–Wiseman);
3. aggiorna la velocità, v ← v + (F/k)·dz, e la posizione, x ← x + v·dz (Verlet).

Il passo dz deve essere piccolo dove le palline sono fitte: dz ≤ 0,2·k·g_min² (esplicito). Oppure si usa
il punto medio implicito (vedi §8): stesso risultato, 20–60 volte più veloce sul fascio libero.

### Regola C — ogni pallina porta un'ampiezza complessa ψ (fase + ampiezza) — la migliore in assoluto

A ogni passo, **una pallina** i:
1. guarda le sue **2m palline più vicine** (per distanza; m = 3, oppure m = 8 se le palline sono state
   messe a caso), e per ognuna legge posizione x_j e ψ_j;
2. adatta ai rapporti log(ψ_j/ψ_i) una parabola a₁d + a₂d² (d = x_j − x_i): è un piccolo "fascetto
   gaussiano" locale E(d) = e^{a₁d + a₂d²};
3. adatta ai resti ψ_j/E(d_j) − ψ_i una parabola c₁d + c₂d²;
4. ne ricava la pendenza e la curvatura della sua onda: ψ' = a₁ψ + c₁, ψ'' = (2a₂ + a₁²)ψ + 2a₁c₁ + 2c₂;
5. **si muove lungo la corrente**: v = Im(ψ*·ψ')/(k|ψ|²), x ← x + v·dz;
6. **aggiorna la sua onda** come la vede chi viaggia con lei: ψ ← ψ + dz·[(i/2k)ψ'' + v·ψ'].

La pallina non legge mai una densità. Le palline restano distribuite come |ψ|² perché ognuna segue la
corrente (è la proprietà di "equivarianza" delle traiettorie di Bohm). Il passo 6 si integra con la regola
del trapezio (Crank–Nicolson), usando la stessa velocità nel moto e nell'onda. La versione esplicita pura
(ogni pallina fa solo il suo conto, Runge–Kutta 4, passi piccoli) dà gli stessi numeri:
N = 100: −0,64% / −0,64% / −0,64% (`python lab/identical_free.py Ce 100`);
N = 300: −0,21% in tutti e tre i punti, come C (`... Ce 300` e `... C 300`).

### Varianti provate e scartate (stessa regola per tutte le palline)

| sigla | cosa fa la pallina | esito sul fascio libero (N = 100, partenza regolare) | comando |
|---|---|---|---|
| B | densità da una curva adattata alle posizioni dei vicini in funzione del loro "numero d'ordine" (Savitzky–Golay, grado 4, m = 4); Q = −(1/2k)(√ρ)''/√ρ e la sua derivata dalla stessa curva | **si blocca** a z = 0,88: nelle code la forza ha il segno sbagliato, le palline si ammassano e il passo diventa minuscolo | `python lab/identical_free.py B 100 q tmax=60` |
| BH | come A ma "lisciata": U = (1/8k) Σ J'²/J⁴ con J, J' dalla curva sui 2m+1 vicini | m = 1: −4,04% / −6,14% / −7,28%; m = 2: **si blocca** a z = 10,9 | `... BH 100 q m=1`, `... BH 100 q m=2 tmax=60` |
| P | porta ampiezza R e fase S separate; v = S'/k, dS/dz = kv²/2 − Q con Q = −(1/2k)R''/R | −4,20% / −13,55% / **si rompe** a z = 70,9 | `... P 100 q tmax=100` |
| H | porta ψ; le palline ferme si scambiano onda solo con le 2 vicine con "molle" simmetriche (conservano la probabilità esattamente), poi ognuna si sposta e rilegge l'onda dai vicini | +2,69% / −4,82% / −51,9% (N = 1000); con partenza casuale si rompe a z = 8,6 | `... H 1000`, `... H 1000 1` |

Note sul percorso (senza numeri, perché ottenuti con versioni intermedie del codice): la prima versione
della regola C (senza il fascetto gaussiano del punto 2, parabola di grado 4, vicini presi per numero
d'ordine) sbagliava le code e poi esplodeva; il fascetto gaussiano l'ha resa esatta sul fascio libero.
Un primo integratore "predictor-corrector" esplodeva per N ≥ 1000 e per passi *più piccoli*: la
causa era l'accoppiamento tra moto e onda; la versione con la stessa velocità nel moto e nei due
operatori (`cn2`) è stabile.

---

## 3. Compito 1 — fascio gaussiano libero (w0 = 3, λ = 1, zR = 28,27)

Larghezza = 2·deviazione standard delle N palline di UN fotone. Teoria w(z) = w0·√(1+(z/zR)²).
"Errore della partenza": quanto la larghezza delle palline di partenza già differisce da w0 (con N
palline messe nei quantili è −1,26% per N = 50 e scende come 1/N). "Errore della regola": larghezza
misurata / (larghezza di partenza × w(z)/w0) − 1, cioè quello che aggiunge la dinamica.

| regola | N | partenza | z = zR | z = 2zR | z = 5zR | errore della regola (5zR) | comando |
|---|---|---|---|---|---|---|---|
| A (esplicita) | 100 | regolare | −0,49% | −0,40% | −0,35% | +0,29% | `python lab/identical_free.py A 100` |
| A (esplicita) | 1000 | regolare | −0,05% | −0,04% | −0,03% | +0,04% | `... A 1000` (117 s) |
| A (implicita) | 1000 | regolare | −0,047% | −0,036% | −0,030% | +0,035% | `python lab/identical_nstudy.py free A 1000 q` (2 s) |
| A | 1000 | casuale (seed 1) | **esplode** | | | | `... A 1000 1 tmax=60` |
| C | 100 | regolare | −0,64% | −0,63% | −0,63% | +0,003% | `... C 100` |
| C | 1000 | regolare | −0,065% | −0,064% | −0,064% | +0,001% | `python lab/identical_nstudy.py free C 1000 q` |
| C (m = 8) | 1000 | casuale (seed 1) | −1,37% | −1,37% | −1,82% | −0,46% | `... C 1000 1 m=8` |

**Tutte le corse riuscite sono entro il 3% richiesto.** La regola C segue le traiettorie esatte di Bohm
pallina per pallina: con N = 1000 lo scarto massimo di una pallina dalla sua traiettoria esatta è
4,5·10⁻⁶, 1,8·10⁻⁴ e 5,0·10⁻⁴ volte la larghezza del fascio a zR, 2zR, 5zR (controllo in linea con
`run_psi(..., scheme='cn2')` di `identical_core.py`, m = 3). Il fascio libero è però un caso
facile per C: il suo "fascetto gaussiano" è esatto per un fascio gaussiano.

---

## 4. Compito 2 — doppia fenditura (L = 160, d = 8, a = 2, w0 = 7)

Le palline partono a z0 = 1 subito dietro le fenditure, distribuite come il flusso esatto j_z(x, z0)
(la partenza "a flusso" dello studio precedente); ogni pallina riceve il valore esatto dell'onda lì
(è la fenditura che la emette). Da z0 in poi tutte seguono la stessa regola. Punteggio: posizioni delle
palline sullo schermo contro il flusso esatto RS.

**Tetto parassiale:** anche un'evoluzione parassiale perfetta (FFT) dalla stessa partenza fa 0,9918
contro il flusso esatto. Quindi 0,97 è raggiungibile in linea di principio.

| regola | N | corr. | err. frange rms / max | palline nelle zone buie (teoria 2,0%) | incroci tra palline | comando |
|---|---|---|---|---|---|---|
| A | 50 | 0,4194 | 52,9% / 75,3% | 20,0% | 0 | `python lab/identical_nstudy.py slit A,C 50,100,200,500,1000 q` |
| A | 100 | 0,5898 | 42,5% / 51,7% | 22,0% | 0 | idem |
| A | 200 | 0,5454 | 37,9% / 58,1% | 18,5% | 0 | idem |
| A | 500 | 0,5987 | 46,0% / 64,5% | 17,0% | 0 | idem |
| A | 1000 | 0,5705 | 46,9% / 65,1% | 18,0% | 0 | idem |
| A | 2000 | non finita: limite di 3·10⁶ passi raggiunto a z = 47,8 (227 s) | | | | `python lab/identical_nstudy.py slit A 2000 q` |
| A | 300 | partenza casuale (seed 1): **esplode** (tutte le palline escono dallo schermo) | | | | `python lab/identical_slit.py A 300 1` |
| C | 50 | 0,3652 | 70,6% / 100% | 14,0% | 2.278 | `identical_nstudy.py slit A,C ...` |
| C | 100 | 0,4032 | 49,6% / 73,5% | 21,0% | 5.860 | idem |
| C | 200 | 0,2607 | 45,5% / 69,5% | 13,0% | 17.621 | idem |
| C | 500 | 0,4509 | 36,9% / 55,1% | 5,4% | 65.445 | idem |
| C | 1000 | −0,0063 | 48,3% / 67,9% | 5,0% | 197.895 | idem |
| C | 3000 | 0,2804 | 45,1% / 59,9% | 1,8% | 810.825 | `python lab/identical_slit.py C 3000 q` |

Per confronto, un **fotone perfetto** di N palline (palline nei quantili della figura parassiale giusta):
N = 50: 0,9574; 100: 0,9827; 200: 0,9881; 500: 0,9910; 1000: 0,9914 (stesse corse). Un campionatore
casuale ideale con 300 clic fa 0,9509 ± 0,0155.

**Le palline attraversano i nodi?**
* Regola A: **no, mai** (0 incroci in tutte le corse: la spinta HDW respinge le vicine). Ma non si formano
  proprio i nodi: le palline restano una nuvola liscia, senza frange, e riempiono le zone buie
  (17–22% contro 2%). La figura è quella di un "fluido classico" che si allarga.
* Regola C: sì, e in massa (da 2·10³ a 8·10⁵ scambi d'ordine): la regola perde il controllo proprio quando
  le due metà del fotone si toccano, e da lì in poi le palline si mescolano.

---

## 5. Compito 3 — quante palline fanno un fotone?

### 5.1 Fascio libero

Partenza **regolare** (palline nei quantili di |ψ|²). Errore sulla larghezza (a zR / 2zR / 5zR) e,
tra parentesi, l'errore della sola regola a 5zR.
Comandi: `python lab/identical_nstudy.py free A,C 50,100,200,500,1000,2000,5000 q`
(per A a 5000: `... free A 500,1000,2000,5000 q`).

| N | errore della partenza | A | C |
|---|---|---|---|
| 50 | −1,262% | −0,99 / −0,82 / −0,73% (+0,54%) | −1,26 / −1,26 / −1,26% (+0,005%) |
| 100 | −0,637% | −0,49 / −0,40 / −0,35% (+0,29%) | −0,64 / −0,63 / −0,63% (+0,003%) |
| 200 | −0,321% | −0,24 / −0,19 / −0,17% (+0,15%) | −0,32 / −0,32 / −0,32% (+0,002%) |
| 500 | −0,129% | −0,10 / −0,07 / −0,06% (+0,066%) | −0,13 / −0,13 / −0,13% (+0,001%) |
| 1000 | −0,065% | −0,05 / −0,04 / −0,03% (+0,035%) | −0,065 / −0,064 / −0,064% (+0,001%) |
| 2000 | −0,033% | −0,02 / −0,02 / −0,01% (+0,018%) | −0,033 / −0,032 / −0,032% (+0,000%) |
| 5000 | −0,013% | −0,009 / −0,007 / −0,005% (+0,008%) | −0,013 / −0,012 / −0,013% (+0,000%) |

* **Come scala:** l'errore della partenza va come **1/N** (N × 100 → errore / 100). L'errore della
  regola A va anch'esso come **1/N** (+0,54% → +0,008%). Quello della regola C è ~0: la regola non
  aggiunge nulla.
* **Soglia:** non c'è un N sotto cui il fascio "si sfalda": anche con 50 palline la larghezza è giusta
  entro 1,3%. Il risultato smette di cambiare entro 0,15% da **N ≈ 500** ed entro 0,07% da N ≈ 1000.

Partenza **casuale** (palline estratte da |ψ|², seed 1, 2, 3), regola C con m = 8.
Comando: `python lab/identical_nstudy.py free C 50,100,200,500,1000,2000,5000 1,2,3`.

| N | errore della partenza (seed 1 / 2 / 3) | errore della regola a 2zR | a 5zR |
|---|---|---|---|
| 50 | −11,9 / −2,0 / +9,2% | +0,004 / +0,003 / −0,87% | −3,4 / +0,005 / −7,0% |
| 100 | −14,8 / −4,9 / +6,6% | +0,003 / +0,003 / −0,29% | +0,004 / +0,004 / −2,9% |
| 200 | −7,6 / −3,8 / +3,0% | ≤ 0,003% | +0,002 / +0,004 / −0,93% |
| 500 | −8,7 / +1,2 / +0,07% | ≤ 0,002% | +0,06% / **esplode** / **esplode** |
| 1000 | −1,4 / +1,3 / +0,7% | ≤ 0,001% | −0,46 / +0,001 / 0,000% |
| 2000 | +0,7 / −0,2 / −0,7% | ≤ 0,001% | −0,30 / +0,001 / −0,10% |
| 5000 | +0,06 / +0,14 / +0,55% | ≤ 0,001% | ≤ 0,001% |

* Con palline casuali l'errore dominante è il **rumore della partenza**, che va come **1/√N**
  (circa 1/√(2N): ±7% con 100 palline, ±2% con 1000). Per stare entro il 3% servono **alcune centinaia
  di palline**; ~1000 per stare tranquilli.
* La regola C su palline disordinate è esatta fino a 2zR ma poi, in 2 casi su 21, si rompe (a N = 500)
  e in altri 7 devia di 0,1–7%. Il motivo è numerico: su punti irregolari il "fit" locale dà un operatore
  che non conserva la probabilità (ha modi che crescono). Con partenza regolare non succede.
* La regola A su palline casuali **esplode subito** (una coppia di palline vicinissime riceve una spinta
  enorme): la regola HDW richiede palline già ben ordinate.

### 5.2 Doppia fenditura

Tabella del §4: la correlazione della regola A **non migliora con N** (0,42 → 0,59 → 0,55 → 0,60 → 0,57
per N = 50 … 1000). Non c'è una soglia oltre la quale "nasce" l'interferenza, almeno fino a 1000 (2000
e 3000 non sono finite nel limite di passi). La regola C si rompe per tutti gli N.

### 5.3 Dove si rompe: due metà "morbide" che si incontrano

Per togliere ogni dubbio sui bordi netti delle fenditure: due fenditure **gaussiane** (larghezza 1,
a x = ±4) con campo analitico esatto. Regola C, partenza regolare, errore della pallina peggiore
rispetto all'onda esatta nel punto in cui si trova (`python lab/identical_meet.py meeting 100,300,1000,3000`,
poi `... meeting 10000` e `... meeting 30000`):

| N | pallina più interna all'inizio | errore su ψ a z = 4 | a z = 6 | a z = 8 |
|---|---|---|---|---|
| 100 | \|x\| = 2,84 | 0,0004% | 0,59% | 18% |
| 300 | 2,64 | 0,0017% | 2,4% | 49% |
| 1000 | 2,46 | 0,01% | 10% | 44% |
| 3000 | 2,30 | 0,05% | 35% | 35% |
| 10000 | 2,14 | 0,24% | 34% | 27% |
| 30000 | 2,01 | 0,89% | 29% | esplosa |

Finché le due metà sono separate la regola è esatta. **Quando si toccano l'errore salta al 10–50%, e
con più palline arriva prima e più grande.** L'errore è esattamente il pezzo d'onda che viene dall'altra
fenditura e che doveva arrivare attraverso la zona vuota: a z = 6 e x = −0,67 quel pezzo vale circa il 10%
dell'onda locale, ed è proprio l'errore misurato con N = 1000. Le palline in più finiscono dove la luce
c'è già; la pallina più interna si avvicina al centro solo come √(ln N) (da 2,84 a 2,01 passando da 100
a 30.000 palline), quindi i due gruppi si toccano prima, in una zona ancora più povera di palline.
Le corse complete fino allo schermo si rompono tutte: N = 300 finisce con ψ fuori scala e 70.533 incroci
(corr 0,48), N = 1000 si ferma a z = 47,1, N = 3000 a z = 18,8
(`python lab/identical_meet.py soft N q m=3 P=2 dz=0.02 vmax=0`).

**Risposta a "quante palline fanno un fotone":**
* per un fascio che si propaga da solo: **~50 palline** messe in ordine (o qualche centinaio messe a caso)
  bastano; l'errore scende come 1/N (in ordine) o 1/√N (a caso);
* per l'interferenza: **nessun numero basta** con palline identiche distribuite come |ψ|². Non è una
  questione di risoluzione: più palline peggiorano l'incontro delle due metà.

---

## 6. Compito 4 — perché la doppia fenditura non riesce, e cosa manca

### 6.1 Senza fase è impossibile (regole A, B, BH)

Una regola che usa solo posizioni e velocità dei vicini non può sapere **dove** mettere le frange.
Prova (`python lab/identical_phase_argument.py`): due doppie fenditure che differiscono solo per una fase
π sulla fenditura di destra.
* Le figure sullo schermo sono opposte: frangia centrale chiara contro scura (probabilità nei 4 bin
  centrali 0,0572 contro 0,0004; correlazione tra le due figure −0,364).
* Ma le palline partono quasi identiche: per il 90% della luce la densità di partenza differisce di meno
  dell'1,8% e la velocità di meno di 0,015. L'informazione sulla fase relativa sta solo in queste
  increspature di pochi percento su scala λ, dovute al poco di luce che dalle due fenditure si
  sovrappone già a z0 = 1.

Una regola senza fase deve ricostruire le frange da increspature di ~1% nella densità: le palline non ci
riescono. È il risultato della regola A: una nuvola liscia, senza frange (§4).

### 6.2 Con la fase non basta (regola C)

La regola C porta la fase ed è esatta finché l'onda sta dove stanno le palline. Ma l'equazione di
Schrödinger sposta pezzi d'onda **attraverso** zone quasi buie: l'onda di una fenditura arriva nella zona
dell'altra prima delle sue palline. Le palline, distribuite come |ψ|², in quelle zone non ci sono, e
quindi nessuno trasporta quel pezzo d'onda. Lo stesso difetto si vede nel beam-splitter (§7): dopo la
separazione un piccolo resto d'onda "sbagliato" non riesce a uscire dal fascio (non ci sono palline tra
i due fasci che lo portino via), si accumula, e a z ≈ 30–35 la regola si rompe.

In parole semplici: **l'onda deve vivere anche dove la luce non c'è, ma le palline stanno solo dove la
luce c'è.**

### 6.3 L'ingrediente minimo

1. **Ogni pallina porta una fase** (un'ampiezza complessa, cioè ψ; con la polarizzazione un vettore di
   Jones). Necessario per §6.1.
2. **Le palline devono coprire anche le zone buie.** Allora la loro densità non può essere |ψ|² ovunque:
   ogni pallina deve portare anche un **peso** μ (la sua quota di fotone), fissato alla nascita. Se si
   muove lungo la corrente, il suo μ resta costante. La figura sullo schermo è la densità **pesata**.
   Tutte le palline seguono ancora la stessa regola: sono uguali, portano solo numeri diversi.
3. Il clic singolo richiede di scegliere **una** pallina con probabilità μ: con pesi uguali è la
   risposta "Many Interacting Worlds" (siamo in uno dei mondi, a caso).

Stato delle prove: la versione estrema di questo schema è la regola di `REGOLA_FOTONE.md`: palline-onda
con peso 0 disposte uniformemente, più una pallina con peso 1. Lì la doppia fenditura dà 0,995. **Una
versione "morbida" con palline che si muovono tutte lungo la corrente e un 30% di palline in più nella
zona buia l'ho provata e non ha funzionato**: si rompe a z = 35,5
(`python lab/identical_meet.py soft 1000 q m=3 P=2 dz=0.02 vmax=0 floor=0.3`). Le palline buie hanno
onda quasi zero, velocità mal determinate e punti sempre più irregolari. Farla funzionare richiede un
operatore locale che conservi esattamente la probabilità su punti irregolari; qui non ci sono riuscito.
Quindi il punto 2 è **motivato ma non dimostrato** con palline che si muovono.

---

## 7. Compito 5 — polarizzazione

**Come la porta una pallina:** ψ_i diventa un vettore di Jones complesso (ψ_x, ψ_y). La regola C lavora
componente per componente. La velocità usa la corrente **sommata** sulle due componenti, così due onde
con polarizzazioni ortogonali non interferiscono.

**Beam-splitter polarizzatore (Malus senza monete, senza pallina speciale).** Fascio w0 = 3 polarizzato
H; a z = 0 ogni pallina trasforma la sua ψ in a⟨a|H⟩e^{+ikαx} + b⟨b|H⟩e^{−ikαx} (α = 0,25). Le palline
si dividono da sole nei due fasci; conto quelle con x > 0 a z = 25 (N = 1000, `python lab/identical_polar.py 1000`):

| θ | frazione nel fascio a | cos²θ | polarizzazione media nel fascio a, \|⟨a\|ψ⟩\|²/\|ψ\|² |
|---|---|---|---|
| 0° | 0,9990 | 1 | 1,0000 |
| 22,5° | 0,8540 | 0,8536 | 0,9970 |
| 45° | 0,5000 | 0,5 | 0,9974 |
| 67,5° | 0,1460 | 0,1464 | 0,9870 |
| 90° | 0,0010 | 0 | — |

(a 0° e 90° una pallina di coda sta dalla parte sbagliata del taglio x = 0; la granularità è 1/N = 0,001.)
In quale fascio va una pallina lo decide **dove si trova** nel fascio: tutto è deterministico e locale.
Limite: oltre z ≈ 30 la regola si rompe (§6.2); il conteggio a z = 25 è già a fasci separati.

**Polarizzatore che assorbe.** Ogni pallina proietta la sua ψ su a. Perché la densità resti |ψ|², una
frazione sin²θ delle palline deve sparire: ogni pallina lancia una moneta **sua** con probabilità
|⟨a|ψ_i⟩|²/|ψ_i|². Sopravvivono 0,8450 / 0,5180 / 0,1570 delle palline contro 0,8536 / 0,5 / 0,1464
(stessa corsa, N = 1000; lo scarto è il rumore binomiale, ±0,016).

**Chi decide il singolo clic?** Con palline identiche il fotone è "la nuvola che resta". Ma un rivelatore
fa **un** clic in **un** punto, o non lo fa. Le palline da sole danno solo la *frazione* (la probabilità),
non il singolo evento. Serve un'informazione in più: quale pallina è "la nostra". È la risposta
Many-Interacting-Worlds: tutti i mondi sono uguali e noi siamo in uno a caso. Questa scelta si può fare
una volta sola alla sorgente, in modo uniforme, perché le palline hanno tutte lo stesso peso. Non è una
pallina con una regola diversa: è un'etichetta. Poi
* beam-splitter: il clic è nel fascio dove va la nostra pallina (probabilità cos²θ, per quanto sopra);
* polarizzatore che assorbe: il fotone passa se la nostra pallina sopravvive alla sua moneta
  (probabilità cos²θ).

---

## 8. Limiti e note oneste

1. **Solo parassiale e 2D.** La regola evolve l'equazione parassiale. Il tetto di questo modello sulla
   doppia fenditura è 0,9918 (FFT). In 3D i "vicini più prossimi" restano definiti, ma i fit diventano
   superfici.
2. **Il passo in z è comune a tutte le palline.** La versione implicita (Crank–Nicolson per C, punto medio
   implicito per A) risolve a ogni passo un sistema lineare "a banda" su tutte le palline. È un artificio
   numerico per fare passi lunghi: la versione esplicita pura, dove ogni pallina fa solo il suo conto
   locale (`Ce`, `A`), dà gli stessi numeri (§2, §3), ma con passi 10–10⁴ volte più corti.
   Il punto medio implicito sulla doppia fenditura non converge: si ferma a z = 9,0 con N = 1000,
   a z = 34,8 con N = 2000 e a z = 6,6 con N = 5000 (`python lab/identical_slit.py Ai N q dzmax=0.05`).
   Lì uso la versione esplicita.
3. **La partenza della doppia fenditura** (flusso esatto e onda esatta a z0 = 1) è il modo in cui la
   fenditura emette le palline, come nello studio precedente. Non contiene ancora l'interferenza lontana.
4. Il limite dei 3·10⁶ passi ha fermato la regola A esplicita a N = 2000 (z = 47,8) e a N = 3000 (z = 4,45,
   `python lab/identical_slit.py A 3000 q`).
5. `identical_meet.py` con partenza "soft" usa un riferimento FFT per la correlazione e il campo
   analitico per il test "meeting".

---

## 9. File e comandi

| file | cosa fa |
|---|---|
| `lab/identical_core.py` | le regole A, B, BH, P, C (e H), gli integratori, le partenze, i riferimenti |
| `lab/identical_free.py RULE N [seed\|q] [chiave=valore]` | fascio libero, una corsa |
| `lab/identical_slit.py RULE N [seed\|q] [chiave=valore]` | doppia fenditura, una corsa (punteggi contro RS e FFT) |
| `lab/identical_nstudy.py free\|slit RULES Ns [seeds]` | studio in N, 4 processi |
| `lab/identical_meet.py CASE N …` / `meeting Ns` | due metà gaussiane che si incontrano; errore all'incontro |
| `lab/identical_phase_argument.py` | fenditure con fase 0 e π: partenze quasi uguali, figure opposte |
| `lab/identical_polar.py [N]` | beam-splitter polarizzatore (Malus deterministico) e polarizzatore che assorbe |

Tempi misurati (i5-10300H, macchina condivisa): fascio libero C con N = 5000, 37–51 s; A implicita con
N = 5000, 83 s; doppia fenditura A esplicita con N = 1000, 84 s; studio in N del fascio libero, < 3 min.
Nessun file esistente in `src/` o in `lab/` è stato modificato.
