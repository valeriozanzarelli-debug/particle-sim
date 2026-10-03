# Molti fotoni: la stessa regola "a due ruoli" per 2, 3, 4 fotoni e per i fasci

Laboratorio: cartella `lab/`, file `multi_*.py` (Python 3.14 + numpy/scipy). Tutti i numeri qui sotto
vengono da esecuzioni reali con seed fissi su una macchina a 4 core occupata da altri lavori; i tempi
sono "a freddo" (tabelle cancellate prima di ogni script) e sono riportati in §9.
La regola del singolo fotone è quella di `docs/REGOLA_FOTONE.md`; qui viene estesa, **una volta sola**,
a n fotoni (`lab/multi_core.py`), e poi usata **senza ritocchi** in tutti gli esperimenti.

---

## 1. In una frase

Ogni fotone ha la sua **pallina-fotone**; le **palline-onda** sono le stesse di sempre. Quando i fotoni
non hanno legami (fotoni indipendenti, o tanti fotoni nello stesso fascio) ogni pallina-fotone ascolta
solo le palline-onda che le passano accanto, esattamente come nella regola del singolo fotone. Quando
i fotoni sono **entangled** o **identici e mescolati** (beam-splitter), l'onda è una **piccola somma di
prodotti** e la velocità di una pallina-fotone dipende anche da **dove si trovano le altre** palline-fotone:
è l'unica eccezione alla località, ed è obbligatoria (teorema di Bell). Il clic è dove la pallina atterra.

Con questa regola fissa: Grangier–Roger–Aspect (α = 0), Mach–Zehnder a un fotone (V = 1,000),
Hong–Ou–Mandel (visibilità 0,9993), Hanbury Brown–Twiss (g2 = 0 / 1,03 / 1,98), Bell–CHSH con il
fotone che sceglie l'uscita del PBS (S = 2,832 ± 0,014 su 40.000 coppie), GHZ (quattro correlazioni ±1
esatte), teletrasporto (successo 0,253, fedeltà 1,000) e scambio di entanglement (fedeltà 1,000),
cancellino quantistico a scelta ritardata (frange/anti-frange solo dopo la selezione). **96 test su 96
PASS** negli 8 script di esperimento (più l'autotest della regola, PASS), compresi i test dei "livelli di zoom".

---

## 2. La regola, scritta come "cosa fa ogni pallina"

### 2.1 Palline-onda (identiche alla regola del singolo fotone)

* Volano dritte; quando toccano un oggetto vengono ri-emesse come **ventaglio a reticolo** (origini
  ogni `DS = 0,1 λ` sugli oggetti lisci, tutte le direzioni in avanti). Ampiezza di un'origine =
  onda in arrivo letta lì (dalle palline che arrivano) × matrice locale dell'oggetto.
* Portano un **vettore interno** = (etichetta di percorso) ⊗ (polarizzazione di Jones). Un oggetto è
  una **matrice di Jones** su questo vettore più una **rampa di fase** locale che fissa l'inclinazione
  d'uscita:
  * BS 50:50: `a → t·c + r·d`, `b → r·c + t·d`, t = 1/√2, r = i/√2; c esce inclinata di +α, d di −α;
  * PBS nella base (e1, e2): `|e1⟩⟨e1|` → uscita "+" (+α), `|e2⟩⟨e2|` → uscita "−" (−α)
    (base lineare a un angolo, oppure circolare R/L: è la lamina λ/4 + PBS);
  * specchi/prismi: cambiano solo l'inclinazione (rampa) e, se serve, una fase e^{iφ}.
  Fattore di flusso √(cos θ_in / cos θ_out), così l'oggetto non crea né distrugge fotoni.
* L'insieme delle palline-onda partite da **una** sorgente e passate per **un** percorso di oggetti
  è un **canale** (etichetta, inclinazione, ventaglio). Canali con la stessa etichetta interferiscono;
  canali con etichette diverse no (come due polarizzazioni ortogonali nella regola singola).
  Perché le due uscite di un BS hanno etichette diverse: in un BS vero escono in **direzioni diverse
  dello spazio 3-D**; nel nostro piano 2-D questa separazione è rappresentata dall'etichetta.

### 2.2 L'onda di n fotoni senza esplosione esponenziale

    Ψ(x1, …, xn) = Σ_k c_k · ψ_k1(x1) ⊗ ψ_k2(x2) ⊗ … ⊗ ψ_kn(xn)

* `ψ_ki` è un'onda **di un solo fotone** (i canali della sua sorgente, con il suo vettore di
  polarizzazione e un'**etichetta spettrale** τ_ki — il "colore/istante" del fotone, che gli oggetti non
  toccano e che serve per ritardi e distinguibilità).
* Fotoni indipendenti: **1 termine**. Coppia entangled |Φ+⟩: **2 termini** (decomposizione di Schmidt).
  Fotoni identici che si incontrano su un BS: termini **simmetrizzati** (bosoni).
* Le palline-onda **non** si moltiplicano con i termini: i termini differiscono solo per vettori
  interni e coefficienti, quindi tutti i fotoni di un braccio leggono **gli stessi canali**.
  Costo = (termini) × (fotoni) letture di onde a un fotone, mai una griglia in 2n dimensioni.

| esperimento | fotoni n | termini K | palline-onda in tutto |
|---|---|---|---|
| GRA, HBT (un ingresso) | 1 per volta (o n nello stesso modo) | 1 | 78.120 (braccio con BS) |
| Mach–Zehnder | 1 | 1 | 2 BS + specchi, stesso ordine |
| HOM | 2 nello stesso braccio | 2 (simmetrizzati) | 78.120 |
| Bell/CHSH | 2 in due bracci | 2 (Schmidt) | 2 × 37.603 |
| GHZ | 3 in tre bracci | 2 | 3 × 37.603 |
| teletrasporto | 3 (1,2 sul BS; 3 al PBS) | 4 | 78.120 + 37.603 |
| scambio di entanglement | 4 (2,3 sul BS) | 8 | 78.120 + 2 × 37.603 |
| cancellino a scelta ritardata | 2 | 2 | 136.152 |

### 2.3 Pallina-fotone (la sola cosa nuova)

A ogni passo la pallina-fotone i legge **solo le palline-onda accanto a sé** (per ogni termine k e
ogni canale) e ottiene ψ_ki(x_i) e il suo gradiente. Da qui due numeri per ogni coppia di termini (k, l):

* `Q^i_kl` = "sovrapposizione di flusso" del fotone i tra i termini k e l, letta **dove si trova**
  (per k = l è il flusso j_z della sua onda, diviso k);
* `X^i_kl` = la stessa cosa per la corrente trasversale j_x.

La pallina-fotone i si muove alla velocità della luce con pendenza

    dx_i/dz = Re Σ_kl c_k* c_l · X^i_kl · Π_{j≠i} Q^j_kl  /  Re Σ_kl c_k* c_l · Π_j Q^j_kl

* **Un solo fotone, o fotoni non legati (un termine):** i fattori degli altri si semplificano e resta
  esattamente `j_x / j_z` della sua onda: la regola del singolo fotone. Verificato: differenza
  5,6·10⁻¹⁷ (autotest, §4).
* **Fotoni legati:** compare `Π_{j≠i} Q^j_kl`, cioè **dove sono adesso gli altri fotoni**, anche
  lontanissimi. **Questa è l'unica eccezione alla località** e non è una scelta: per il teorema di Bell
  nessuna regola locale può riprodurre le correlazioni misurate (S > 2). Le palline-onda restano locali;
  la non località è solo in questo prodotto.
* **Perché funziona sempre:** le "correnti di transizione" tra due onde che soddisfano Helmholtz hanno
  divergenza nulla, quindi la densità `ρ = Re Σ c_k* c_l Π_j Q^j_kl` si conserva lungo il moto
  (equivarianza esatta): fotoni nati con la densità di Born restano distribuiti secondo Born.
  Non c'è nessun collasso, nessuna moneta ai PBS e ai BS.
* **Nascita:** la sorgente piazza **tutte** le palline-fotone insieme con densità ρ (unico numero
  casuale: è l'evento della sorgente, per esempio la coppia che nasce nel cristallo).
* **Oggetti che non assorbono** (BS, PBS, specchi): nessuna moneta; il primo passo dietro il piano
  (lunghezza `DELTA = 0,5 λ`) si fa lungo la corrente uscente (Eulero implicito, come nella regola
  singola).
* **Rivelatore:** la pallina si ferma sul piano del rivelatore; il clic è l'uscita (o la posizione) dove
  è atterrata. Una pallina già rivelata resta ferma ma continua a contare, con la sua posizione, nella
  guida delle altre (è così che "la misura di Alice decide per Bob" senza collasso).
* **Oggetti che assorbono** (non usati in questi esperimenti): moneta locale p = Φ_out/Φ_in con il
  flusso della propria onda condizionata (stessa formula con Q prima e dopo l'oggetto).
* **Orologio comune:** tutte le palline-fotone avanzano insieme in z (passi RK4,
  `dz = min(2, max(0,05, 0,2·distanza dall'ultimo piano))`). Serve un tempo comune perché la regola non
  locale ha bisogno di un "adesso"; le statistiche osservabili non dipendono da questa scelta
  (verificato in §5.4), i singoli esiti sì.

### 2.4 Parametri (gli stessi per tutti gli esperimenti)

| parametro | valore |
|---|---|
| λ | 1 |
| passo origini `DS` | 0,1 λ (oggetti lisci, REGOLA §3.3) |
| primo passo dietro un oggetto `DELTA` | 0,5 λ (= 5 DS) |
| passi della pallina-fotone | RK4, dz = min(2, max(0,05, 0,2·Δz)) |
| fasci | gaussiani, w0 = 5 λ (solo le due sorgenti strette del segnale nel cancellino: 1,5 λ) |
| deflessione di BS e PBS | α = ±0,25 rad |
| geometria standard | oggetto a z = 60, rivelatori a z = 160 (le uscite sono separate di ≈ 3 larghezze) |
| uscita di un clic | etichetta del fascio uscente il cui centro è più vicino al punto d'atterraggio |
| tabelle | griglia che segue il fascio: Δu = 0,25 λ (Hermite), Δz = 0,5 λ (lineare) |

Nessun parametro della regola è cambiato da un esperimento all'altro; cambiano solo la geometria
(dove stanno sorgenti, oggetti, rivelatori) e lo stato iniziale (i termini).

### 2.5 Accelerazione (non cambia la regola)

La somma dei vicini di un canale è la stessa funzione per tutti i fotoni che la leggono, quindi viene
**tabulata una volta** su una griglia che segue il fascio (cache su disco in `lab/out/multi_cache/`,
non versionata). Verifiche (autotest): errore massimo tabella / somma diretta 9,7·10⁻⁶ su campo e
gradienti; lettura "a cappello" della REGOLA = somma diretta (differenza 0); 40 traiettorie con tabella
e con somma diretta: |Δx| ≤ 2,4·10⁻⁴ λ al rivelatore (fasci larghi ~8 λ). La tabella è valida finché
nessuna pallina-onda viene uccisa tra l'oggetto e la pallina-fotone (vero in tutti questi esperimenti:
gli oggetti sono piani interi). Con la tabella una pallina-fotone costa ~0,4 ms per tutto il suo volo.
Anche le **coppie di canali** che possono interferire sono precalcolate: ogni passo è un solo prodotto
di matrici per fotone (questo ha portato lo scambio di entanglement da >9 min a 2 min).

---

## 3. Tabella dei risultati

σ = errore statistico; la tolleranza è 3σ (o 3/N quando la teoria è 0 o 1 e σ → 0), salvo dove
indicato. Tutti i conteggi sono eventi completi (coppie, terne, quaterne).

| # | esperimento | misurato | teoria | tolleranza | esito | conteggi | tempo |
|---|---|---|---|---|---|---|---|
| a1 | GRA: anticorrelazione α, fotone singolo | **0,000** (0 coincidenze; P_c = 0,5024) | 0 (GRA 1986: 0,18 ± 0,06; classico ≥ 1) | α < 0,05 | PASS | 20.000 inneschi | |
| a2 | GRA: impulsi laser μ = 0,5 | α = 0,975 ± 0,022 | 1 | 3σ | PASS | 40.000 inneschi, 19.903 fotoni | |
| a3 | Mach–Zehnder a un fotone, P(e) a φ = 0…315° (8 fasi) | 1,0000 / 0,8564 / 0,5016 / 0,1508 / 0,0000 / 0,1392 / 0,4888 / 0,8516 | cos²(φ/2) | 3σ | 8 × PASS | 8 × 2500 | |
| a4 | … visibilità delle frange a un fotone | 1,0000 | 1 | > 0,98 | PASS | | |
| b1 | HOM: coincidenze a ritardo 0 | 0,0003 (1 su 3000) | 0 (con la forma reale dei fasci: 6·10⁻⁵) | 3/N | PASS | 3000 coppie | |
| b2 | HOM: coincidenze a τσω = 0,5 / 1 / 1,5 / 2 / 3 / 5 | 0,1103 / 0,3110 / 0,4557 / 0,4913 / 0,4953 / 0,4980 | 0,1106 / 0,3161 / 0,4473 / 0,4908 / 0,4999 / 0,5000 | 3σ | 6 × PASS | 6 × 3000 | |
| b3 | HOM: visibilità del buco | **0,9993** | 1 | ≥ 0,95 | PASS | | |
| b4 | HOM: polarizzazioni a θ = 0 / 30 / 45 / 60 / 90° | 0,0000 / 0,1310 / 0,2660 / 0,3587 / 0,4917 | ½ sin²θ: 0 / 0,125 / 0,25 / 0,375 / 0,5 | 3σ | 5 × PASS (45° e 60° a 2σ) | 5 × 3000 | |
| c1 | HBT fotone singolo g2(0) | **0,000** | 0 | 3σ | PASS | 15.000 finestre | |
| c2 | HBT laser g2(0) | 1,026 ± 0,015 | 1 | 3σ | PASS | 60.000 finestre | |
| c3 | HBT termica g2(0) (M = 50 emettitori) | 1,980 ± 0,034 | 1,945 per questo campo (ideale 2 − 1/M = 1,98) | 3σ | PASS | 60.000 finestre | |
| c4 | HBT termica g2(τ), τ = 0,5 / 1 / 1,5 / 2 / 3 | 1,717 / 1,294 / 1,060 / 0,988 / 1,028 | 1,720 / 1,297 / 1,048 / 1,008 / 1,017 | 3σ | 5 × PASS | | |
| d1 | Bell/CHSH, pallina-fotone decide al PBS, PBS simultanei | **S = 2,8322 ± 0,0141** | 2√2 = 2,8284 | 3σ | PASS | 40.000 coppie | 55 s |
| d2 | E(a,b) per (0°,22,5°) (0°,67,5°) (45°,22,5°) (45°,67,5°) | +0,7100 / −0,7040 / +0,7058 / +0,7124 (± 0,007) | ±0,7071 | 3σ | 4 × PASS | 4 × 10.000 | |
| d3 | nessun segnale: P(A+) al variare di b, P(B+) al variare di a | differenze +0,005 / +0,011 / −0,010 / +0,004 | 0 (e 0,5) | 3σ = 0,021 | 4 × PASS | | |
| d4 | E(0°, b), b = 0…90° a passi di 15° | +1,000 / +0,863 / +0,492 / +0,001 / −0,522 / −0,860 / −1,000 | cos 2b | 3σ | 7 × PASS | 7 × 2000 | |
| d5 | CHSH con il PBS di A (z = 60) prima di quello di B (z = 100) | S = 2,8570 ± 0,0313 | 2,8284 | 3σ | PASS | 8000 coppie | 41 s |
| d6 | modello a collasso di `entangle.rs` (stessi 40.000) | S = 2,840 ± 0,014 | 2,8284 (Rust: 2,840) | — | riferimento | 40.000 | |
| e | GHZ: E(XYY), E(YXY), E(YYX), E(XXX) | −1,0000 / −1,0000 / −1,0000 / **+1,0000** (0 terne di segno sbagliato) | −1 / −1 / −1 / +1 (realismo locale: E(XXX) = −1) | 3/N | 5 × PASS | 4 × 4000 terne | 37 s |
| f1 | teletrasporto, ingressi H, V, +45, −45, R, L: successo della misura di Bell | 0,2607 / 0,2663 / 0,2487 / 0,2433 / 0,2437 / 0,2573 (media 0,2533) | 1/4 | 3σ | 6 × PASS | 6 × 3000 | |
| f2 | … fedeltà dato il successo | **1,0000** per tutti e 6 (730–799 eventi ciascuno) | 1 (limite classico 2/3) | 3σ | 6 × PASS | | |
| f3 | … fedeltà se la misura di Bell fallisce / senza selezione | 0,327–0,346 / 0,498–0,508 | 1/3 / 1/2 | — | coerente | | |
| f4 | scambio di entanglement: E_ZZ, E_XX, E_YY dei fotoni 1-4 dato il successo | −1,0000 / −1,0000 / −1,0000 (senza selezione −0,024 / −0,003 / −0,009) | −1 (senza selezione 0) | 3/N | 3 × PASS | 3 × 6000 quaterne | |
| f5 | … successo / fedeltà con \|Ψ−⟩ | 0,2570 / 0,2525 / 0,2522; F = **1,0000** | 1/4; 1 (separabile ≤ 1/2) | 3σ; > 0,97 | 4 × PASS | ~1500 successi per base | |
| g1 | cancellino a scelta ritardata: tutti i clic del segnale | V = 0,020 ± 0,016 (nessuna frangia) | 0,001 | 3σ | PASS | 8000 coppie | |
| g2 | … segnale dato D1 / D2 | V = 0,801 / 0,860 (± 0,032), fasi −1,61 / +1,61 | 0,857 / 0,857, fasi −1,57 / +1,57 | 3σ | 2 × PASS | ~2000 ciascuno | |
| g3 | … segnale dato D3 / D4 (which-path) | V = 0,077 / 0,045 | 0,002 / 0,002 | 3σ | 2 × PASS | ~2000 ciascuno | |
| g4 | … frange D1 contro anti-frange D2 | differenza di fase 3,07 rad | π | ± 0,3 | PASS | | |
| g5 | … frazioni D1/D2/D3/D4 | 0,2497 / 0,2520 / 0,2485 / 0,2497 | 0,2501 / 0,2501 / 0,2499 / 0,2499 | 3σ | 4 × PASS | | |

(La colonna "tempo" completa per script è in §9.)

---

## 4. Autotest della regola (`python lab/multi_core.py`)

| controllo | risultato |
|---|---|
| tabella contro somma diretta delle palline (campo e gradienti, 3 canali, 1200 punti) | errore relativo max 9,7·10⁻⁶ |
| lettura "a cappello" della REGOLA contro somma diretta | identiche |
| due fotoni in uno stato prodotto contro due regole a un fotone | |Δv| max 5,6·10⁻¹⁷ |
| 40 traiettorie con tabella contro somma diretta | |Δx| max 2,4·10⁻⁴ λ |

---

## 5. Gli esperimenti (che cosa si vede)

### 5.1 Grangier–Roger–Aspect e Mach–Zehnder (`multi_gra.py`)
Un fotone = una pallina-fotone: dopo il BS l'onda va in entrambe le uscite, la pallina in una sola,
scelta dalla sua posizione nel fascio (metà sinistra/destra, in pratica). Quindi **zero coincidenze**:
α = 0 (GRA misurarono 0,18 per i fotoni spuri della sorgente). Con impulsi laser (tante palline nella
stessa onda) α = 0,975 ± 0,022: le palline di un impulso sono indipendenti. Se le due strade si
ricombinano su un secondo BS (Mach–Zehnder) la pallina, guidata dall'onda che è passata da **entrambe**
le strade, esce da "e" con probabilità cos²(φ/2): interferenza di un fotone solo, visibilità 1,000.

### 5.2 Hong–Ou–Mandel (`multi_hom.py`)
Due palline-fotone nello stesso braccio, onda simmetrizzata. Dopo il BS l'onda a due fotoni è
`c(x1)c(x2) − d(x1)d(x2)`: non ha supporto dove un fotone è in c e l'altro in d, e poiché le palline
seguono il flusso e non attraversano i nodi, **finiscono sempre insieme**. Un ritardo (fase e^{iωτ}
sull'etichetta spettrale) o una polarizzazione ruotata riducono la sovrapposizione e le coincidenze
risalgono come ½(1 − |⟨τ1|τ2⟩|²|⟨p1|p2⟩|²). Visibilità 0,9993 (1 coincidenza su 3000 a τ = 0:
condizioni iniziali nelle code dei fasci; la teoria con la forma reale dei fasci dà 6·10⁻⁵).

### 5.3 Hanbury Brown–Twiss (`multi_hbt.py`)
Stessa onda, stesso BS; cambia solo **quante** palline-fotone la sorgente mette nell'onda in ogni
finestra: 1 (fotone singolo, g2 = 0), Poisson(μ) (laser, g2 = 1,03 ± 0,02), Poisson(μ·|A(t)|²) con A
somma di 50 emettitori a fase e frequenza casuali (termica, g2 = 1,98 ± 0,03, poi g2(τ) scende a 1 dopo
il tempo di coerenza). Il raggruppamento termico **nasce nella sorgente** (l'intensità dell'onda
fluttua), non nel BS: lì ogni pallina sceglie da sola. È esattamente la spiegazione di HBT.

### 5.4 Bell/CHSH (`multi_bell.py`)
Lo stato |Φ+⟩ ha due termini; ogni fotone va a un PBS nel suo braccio; **l'esito è l'uscita da cui la
pallina-fotone esce dal PBS**. Nessun collasso, nessuna moneta: le uscite sono decise dalle posizioni
iniziali (xA, xB) e dalla dinamica non locale. Risultato S = 2,8322 ± 0,0141 (40.000 coppie), contro
2,840 del modello a collasso (`entangle.rs` / `validation.rs`, 40.000 coppie) e 2√2 = 2,8284.
Marginali sempre 0,50 ± 0,005: **nessun segnale**.
Cosa fa la non località, misurato con le **stesse posizioni iniziali** cambiando solo l'angolo di Bob:
* PBS simultanei: l'esito di **Alice** cambia nel 35% delle coppie (la sua pallina "sente" l'apparato di
  Bob mentre attraversa il suo PBS) — ma la statistica di Alice non cambia;
* PBS di Alice prima: l'esito di Alice cambia nello **0,00%** delle coppie, quello di Bob nel 72%: chi
  arriva prima decide da solo (con la sua posizione), chi arriva dopo è guidato da dove l'altro è
  già atterrato. Anche così S = 2,857 ± 0,031.

### 5.5 GHZ (`multi_ghz.py`)
Tre bracci, (HHH + VVV)/√2, PBS in base X (±45°) o Y (R/L). Le palline danno E(XYY) = E(YXY) = E(YYX) =
−1 e E(XXX) = +1 su **ogni** terna (zero eccezioni su 16.000), mentre il realismo locale, dai primi
tre, prevede E(XXX) = −1. Ogni singolo fotone resta 50/50 (marginali 0,49–0,52).

### 5.6 Teletrasporto e scambio di entanglement (`multi_teleport.py`)
La misura di Bell lineare è un BS: successo = le palline 1 e 2 escono da porte diverse. Con le palline
questo avviene nel 25,3% dei casi (teoria 1/4) e **solo** per la parte antisimmetrica |Ψ−⟩; in quei
casi la pallina 3, al suo PBS, sceglie sempre l'uscita "stato teletrasportato": fedeltà 1,0000 per H, V,
±45°, R, L (limite classico 2/3). Quando le palline 1 e 2 escono insieme la fedeltà è 1/3, senza
selezione 1/2: nessuna informazione viaggia senza il canale classico. Nello scambio, i fotoni 1 e 4
(mai incontrati) mostrano E_ZZ = E_XX = E_YY = −1 **solo** sugli eventi in cui 2 e 3 escono separati
(F = 1,000 con |Ψ−⟩), e correlazioni ~0 sugli altri.

### 5.7 Cancellino quantistico a scelta ritardata (`multi_eraser.py`)
Il fotone segnale (da A o da B) arriva allo schermo a z = 50 e fa clic; l'idler, nel suo braccio, passa
dai "rubinetti" which-path a z = 30 e dal BS che cancella a z = 60, **dopo** il clic del segnale.
Tutti i clic del segnale: nessuna frangia (V = 0,02). Selezionati con D1: frange (V = 0,80; teoria
0,857); con D2: anti-frange (fase opposta, 3,07 rad); con D3/D4: niente frange. Nel quadro delle
palline nulla va indietro nel tempo: il segnale atterra per primo, guidato da un'onda condizionata che
non mostra frange (l'idler è già in uno dei due fasci a o b); poi è l'idler a essere guidato verso D1 o
D2 in modo correlato a **dove il segnale è già atterrato**.

---

## 6. MOLTI fotoni: i "livelli di zoom" della luce (`multi_beams.py`)

| livello | che cosa si simula | quando basta (criterio) | costo | errore misurato |
|---|---|---|---|---|
| **1 — quantistico** | palline-onda + una pallina-fotone per fotone, onda a n fotoni come somma di K prodotti | sempre; **necessario** per luce non classica: fotoni singoli, coppie, entanglement, HOM, sub-Poisson | K × n letture per passo; K cresce come n! se n fotoni identici in modi diversi si mescolano | — (riferimento) |
| **2 — onda condivisa** | **una** sola onda per modo, n palline-fotone indipendenti | fotoni nello stesso modo (laser, luce termica in un modo): stato prodotto, la regola a n fotoni = n regole a un fotone (|Δv| = 6·10⁻¹⁷) | palline-onda + n (non n × palline-onda) | clic Poisson (Fano 0,977 ± 0,033), indipendenti (g2 = 1,007 ± 0,036), figura = quella a un fotone (corr 0,9931) |
| **3 — solo onda (Maxwell)** | niente palline-fotone; i clic, se servono, sono estrazioni di Poisson del flusso j_z | luce classica e domanda sulla **media**: l'unica differenza col livello 2 è il rumore granulare √((1−p)/(n p)) per cella; con n̄ ≫ 1/ε² clic per cella il livello 2 non aggiunge nulla | indipendente da n | intensità dell'onda contro Rayleigh–Sommerfeld esatto: 0,02%; clic estratti dall'onda contro clic delle palline: χ²/gdl = 1,69 (stessa distribuzione entro 3σ) |
| **4 — raggi** | niente palline-onda; la pallina-fotone va dritta (si piega solo sui bordi) | aperture ≫ λ: numero di Fresnel N_F = a²/(λL) ≫ 1 | O(1) per fotone per oggetto | luce fuori dal fascio geometrico ≈ 0,11/√N_F; distanza totale (0,231 + 0,051 ln N_F)/√N_F: **1% per N_F > 4400, 10% per N_F > 13** |

### 6.1 Livello 1 → 2: tanti fotoni nello stesso modo (laser)
4000 impulsi laser con μ = 4 (15.966 palline-fotone emesse, 3832 entrano nelle fenditure, 3374 clic)
sulla doppia fenditura della REGOLA (stessa regola a un fotone, variante G veloce: ds = 0,04, δ = 0,2),
**tutte con lo stesso ventaglio di 31.400 palline-onda**:
* figura di tutti i clic contro il flusso esatto: corr 0,9931, χ²/gdl 1,84 (REGOLA variante G: 0,9938);
* clic per impulso: Fano (varianza/media) 0,977 ± 0,033 sullo schermo, 1,009 ± 0,057 nella frangia
  centrale, 0,974 ± 0,041 a destra (Poisson: 1); g2 sinistra/destra 1,007 ± 0,036: i fotoni di un
  impulso non si "parlano".
* **Costo misurato.** Far volare le 31.400 palline-onda fino allo schermo: 3 ms; una pallina-fotone che
  le legge (somma esatta su 100 origini, ~150 passi RK4): 36 ms. Memoria: 1,5 MB per l'onda condivisa,
  contro 1,5 GB per 1000 fotoni con un'onda ciascuno (1,5 TB per un milione). Il tempo di calcolo è
  dominato dalla lettura del fotone (∝ origini); con l'onda **tabulata** (`multi_core`) una pallina costa
  0,36–0,40 ms per tutto il volo attraverso un BS e il tempo cresce come n^0,96 (1000, 4000, 16.000
  palline: 0,4 / 1,6 / 5,8 s): **costo ∝ palline-onda (una volta) + n**.

### 6.2 Livello 2 → 3: fascio brillante, solo onda
* L'onda da sola (flusso j_z delle palline-onda, nessuna pallina-fotone) riproduce il flusso esatto di
  Rayleigh–Sommerfeld entro lo 0,02% sui bin luminosi.
* L'errore delle palline-fotone rispetto all'onda, per bin luminosi di 8 λ: 52% / 31% / 20% / 15% con
  100 / 300 / 1000 / 3374 clic. Fino a ~300 clic è puro rumore granulare (54% / 31%); oltre compare un
  **pavimento sistematico dell'11,8% per bin** (4,5% rms con la metrica delle frange della REGOLA):
  è l'errore della moneta al bordo delle fenditure della variante veloce (REGOLA §7.1), che il livello 3
  **non** ha. Il rumore da solo arriva all'1% in ogni bin luminoso con ~9·10⁵ clic.
* Clic estratti a caso dal flusso dell'onda (livello 3) contro clic delle palline (livello 2):
  χ²/gdl a due campioni 1,69 (soglia 3σ: 1,81) → stessa distribuzione.
* **Criterio.** Le palline-fotone servono solo per luce **non classica**. Per luce classica (laser,
  termica) la scelta tra livello 2 e 3 dipende solo dal rumore che si vuole vedere: si passa al livello 3
  quando il numero medio di fotoni per modo e per cella di rivelazione n̄ ≫ 1 (in pratica ≫ 1/ε² per
  un errore ε). Esempi: laser He-Ne 1 mW con riga 1 MHz: n̄ = 3·10⁹ fotoni per modo → solo onda;
  luce solare a 500 nm: n̄ = 0,007 per modo — pochissimi fotoni per modo, ma è luce classica: i clic
  sono estrazioni di Poisson da un'onda che fluttua (termica), mai serve il livello 1.

### 6.3 Livello 3 → 4: raggi
Fenditura di semilarghezza a illuminata da un'onda piana, schermo a distanza L, soluzione esatta di
Fresnel:

| N_F | 0,25 | 1 | 4 | 16 | 64 | 256 | 1024 | 4096 | 16.384 |
|---|---|---|---|---|---|---|---|---|---|
| distanza totale onda–raggi | 26,1% | 22,4% | 15,1% | 9,3% | 5,6% | 3,2% | 1,83% | 1,03% | 0,57% |
| luce fuori dal fascio geometrico | 25,9% | 12,6% | 5,8% | 2,8% | 1,39% | 0,69% | 0,35% | 0,17% | 0,086% |

* La luce fuori dal fascio scende come N_F^−0,505 (atteso −1/2: la zona di bordo è larga √(λL)).
* La distanza totale scende più lentamente: N_F^−0,40 se si forza una legge di potenza; la forma giusta
  è (0,231 + 0,051 ln N_F)/√N_F (scarto ≤ 0,2%), perché le increspature di Fresnel dentro il fascio
  calano solo come 1/(distanza dal bordo). **Raggi buoni all'1% per N_F > 4400, al 10% per N_F > 13.**
* Con le palline (fenditura a = 4 λ, regola del singolo fotone, 1500 fotoni): le palline che si
  scostano dal raggio dritto di più di 0,1 a sono il 12,0% / 3,0% / 0,0% a N_F = 1 / 4 / 16 (scarto
  mediano 0,081 a / 0,010 a / 0,002 a), e l'istogramma dei clic dista dai raggi quanto l'onda esatta
  (24,6% vs 22,4%, 15,6% vs 15,1%, 10,2% vs 9,3%). **La pallina va dritta, tranne vicino ai bordi.**

### 6.4 Raggi gamma
λ = 10⁻¹² m, apertura 1 mm (a = 0,5 mm), schermo a 1 m: **N_F = 2,5·10⁵**. Errore dei raggi (Fresnel
esatto): distanza totale 0,17%, luce fuori dal fascio geometrico 0,022%; zona di bordo √(λL) = 1 µm su
500 µm. Con la regola a palline-onda: origini ogni 0,1 λ → in 2-D 10¹⁰ origini (3·10¹¹ palline), in 3-D
**8·10¹⁹ origini e 2,4·10²¹ palline** (≈ 10²³ byte), e ogni passo di ogni fotone leggerebbe 8·10¹⁹
ondine; anche solo la zona di bordo in 3-D ha ~10¹⁸ origini. Con i raggi: un raggio per fotone, costo
O(1). **Per i gamma il livello 4 è obbligatorio** (e lì contano altre fisiche — Compton, produzione di
coppie — che sono interazioni da particella, non diffrazione).

### 6.5 Regola di scelta (riassunto)
1. La luce è **non classica** (fotoni singoli, coppie, entanglement, HOM)? → livello 1, sui soli
   fotoni coinvolti; K termini piccoli.
2. Altrimenti: n fotoni nello stesso modo → livello 2 (una sola onda + n palline), oppure livello 3 se
   interessa l'intensità media o n̄ per cella ≫ 1/ε².
3. In ogni caso, per ogni oggetto: se N_F = a²/(λL) > ~4400 (1%) o > 13 (10%) → raggi (livello 4)
   per quel tratto; le palline-onda solo dove N_F è piccolo (bordi, aperture strette).

---

## 7. Che cosa non funziona, o funziona solo con condizioni (onestà)

1. **Due tentativi falliti nei test dei livelli di zoom.** (a) Avevo previsto che l'errore delle palline
   contro l'onda scendesse come n^−1/2: misurato n^−0,35, perché sotto il rumore c'è il pavimento
   sistematico dell'11,8% della variante veloce (vedi §6.2); il test è stato cambiato in "fino a 300
   clic = rumore puro" e il pavimento è riportato. (b) Avevo previsto distanza totale ∝ N_F^−1/2:
   misurato −0,40, per il logaritmo delle increspature di Fresnel; la legge −1/2 vale per la luce fuori
   dal fascio (−0,505). Nessuna tolleranza è stata allargata: sono state corrette due previsioni
   sbagliate, e i fallimenti sono scritti qui.
2. **La non località richiede un "adesso" comune** (qui: la coordinata z, uguale per tutte le palline).
   Le statistiche non ne dipendono (S uguale con PBS simultanei o sfalsati), ma **i singoli esiti sì**:
   con le stesse posizioni iniziali l'esito di Alice cambia nel 35% delle coppie se si cambia l'angolo
   di Bob quando i PBS sono simultanei. È il prezzo noto di ogni teoria a variabili nascoste non locale;
   non permette segnali (marginali 0,50 ± 0,005).
3. **Le uscite di BS/PBS sono etichette ortogonali.** Nel piano 2-D due fasci che si incrociano
   interferirebbero davvero; un BS vero li manda in direzioni diverse nello spazio. L'etichetta
   rappresenta questa terza dimensione. Il porting 3-D (agente `base_*`) potrà usare direzioni vere.
4. **Fotoni identici in molti modi diversi**: la simmetrizzazione richiede fino a n! termini (permanenti:
   è il problema del "boson sampling", difficile anche per la natura). Il livello 1 va usato solo sui
   pochi fotoni davvero indistinguibili e mescolati; i fasci classici passano ai livelli 2–4.
5. **Il raggruppamento termico (HBT) è messo dalla sorgente** (intensità fluttuante, numero di palline
   di Poisson): la dinamica delle palline al BS non lo crea né lo distrugge. È fisicamente corretto
   (g2 è proprietà del campo), ma non è un test della guida.
6. **HOM: una coincidenza su 3000 a ritardo zero** (attese 0,2 con la forma vera dei fasci): errore
   numerico di integrazione per partenze nelle code dei due fasci. Con polarizzazioni a 45° e 60° gli
   scarti sono a 2σ (dentro la tolleranza ma segnalati).
7. **Tabelle**: valide solo se nessuna pallina-onda viene uccisa fra l'oggetto e la pallina-fotone
   (oggetti a piano intero). Per fenditure/ostacoli si torna alla somma diretta (REGOLA), più lenta.
8. **Moneta degli oggetti assorbenti** a più fotoni: definita (flusso condizionato) ma non esercitata in
   questi test; tutti gli oggetti qui sono non assorbenti e i rivelatori stanno sulle uscite.
9. **Rivelatori ideali** (efficienza 1, niente buio, niente conteggi accidentali): i numeri sono quelli
   "ideali" del catalogo `ESPERIMENTI_FOTONE.md`; GRA misurò α = 0,18 per imperfezioni della sorgente.
10. **Solo 2-D** (come la REGOLA). Il pavimento sistematico della §6.2 viene dalla variante veloce
    della regola a un fotone (ds = 0,04); con la regola finale (ds = 0,02) è più basso ma 4 volte più lento.

---

## 8. Interfaccia (compatibile con la regola 3-D `base_*`)

* **Pallina-onda**: come nella REGOLA (`origin`, `dir`, `amp`, `jones`, `fan`, `src`, `dir_idx`) con il
  vettore di Jones esteso a (etichetta di percorso) ⊗ (polarizzazione). Non serve altro.
* **Oggetto**: per ogni etichetta d'ingresso una lista di (etichetta d'uscita, inclinazione d'uscita,
  matrice di Jones 2×2 o scalare). La rampa di fase e il fattore di flusso seguono dall'inclinazione.
* **Pallina-fotone**: posizione, pendenza, braccio, fermata (bool). Nessuno stato quantistico proprio.
* **Stato a n fotoni**: lista di termini `(c_k, [(sorgente, jones, etichetta_spettrale) per fotone])`.
  L'unica operazione non locale è `slopes()`: raccoglie per ogni fotone le matrici K×K Q e X (lette
  localmente) e fa il prodotto sugli altri fotoni.
* Per il Rust: `entangle.rs` (collasso) può restare come modalità classica; la modalità pilota usa
  `Pair = {termini, palline-fotone}` e nessuna `measure()`: l'esito è la porta d'uscita.

---

## 9. Come rieseguire (tempi a freddo, tabelle cancellate prima di ogni script)

Tutti gli script: `timeout 570 python lab/<script>.py`, seed fissi, stampa con `flush=True`, riga finale
`RESULT: n/n PASS`. Log in `lab/out/<script>.log`.

| script | cosa fa | test | tempo a freddo (parete) |
|---|---|---|---|
| `multi_core.py` | la regola + autotest | SELFTEST PASS | 43 s |
| `multi_gra.py` | GRA (α), impulsi laser, Mach–Zehnder | 11/11 | 204 s |
| `multi_hom.py` | HOM: ritardo e polarizzazione | 13/13 | 119 s |
| `multi_hbt.py` | HBT: fotone singolo, laser, termica, g2(τ) | 9/9 | 91 s |
| `multi_bell.py` | CHSH (40.000 coppie), E(0,b), nessun segnale, PBS sfalsati, dipendenza dall'angolo lontano, modello a collasso | 17/17 | 111 s |
| `multi_ghz.py` | GHZ a tre fotoni | 5/5 | 37 s |
| `multi_teleport.py` | teletrasporto (6 stati) e scambio di entanglement | 19/19 | 229 s |
| `multi_eraser.py` | cancellino a scelta ritardata | 11/11 | 127 s |
| `multi_beams.py` | livelli di zoom, raggi, gamma (4 processi) | 11/11 | 115 s |

Il tempo a freddo include la costruzione delle tabelle (da 20 s a ~110 s per il cancellino); con le
tabelle già in `lab/out/multi_cache/` (~95 MB, non versionate) gli script sono più rapidi. Tutti sotto
i 5 minuti con la macchina occupata. Le esecuzioni a freddo hanno dato esattamente gli stessi numeri di
quelle precedenti (seed fissi, calcolo deterministico).
