# Molecole dalle palline: il primo livello (atomi → molecole)

Codice: `src/molecole/` (Rust). La regola è in `dmc.rs`, le misure in `misura.rs`, il visore nel binario
`molecole`. Tutti i numeri di questo documento vengono da esecuzioni reali con seed fissi
(`cargo run --bin molecole -- --misura <bersaglio>`, tutto insieme: `sh lab/molecole_misura.sh`) e sono
salvati in `assets/molecole/risultati.txt` (il visore e i test li leggono da lì). Macchina: i5-10300H,
4 core, condivisa con altri lavori: i tempi sono indicativi.

**Niente PySCF, niente Hartree–Fock, niente DFT, nessun programma di chimica esterno.** Le molecole
vengono dalle stesse palline del resto del progetto, con una regola semplicissima per pallina.

---

## 1. In una frase

Ogni elettrone è fatto di **palline-onda**; ogni pallina, a ogni tick, fa **un passo a caso** e poi
**si duplica o sparisce** a seconda dell'energia del punto in cui si trova. Dopo un po' la nuvola delle
palline **è** la funzione d'onda dello stato fondamentale, e il numero che tiene ferma la popolazione
**è** l'energia. Facendo questo per tante distanze tra i nuclei si ottiene la curva E(R): il suo minimo è
la lunghezza del legame, la sua curvatura (con le masse dei nuclei) è la vibrazione, la differenza tra
il fondo della curva e gli atomi separati è l'energia di legame. Le differenze di energia tra livelli
sono i **fotoni** che l'atomo emette.

Con questa sola regola: idrogeno −13,6057 eV (errore 0,02 meV), Lyman-α 121,565 nm (NIST 121,567),
Balmer-α 656,54 nm (NIST 656,46), polarizzabilità 4,60 ± 0,12 a₀³ (esatta 4,5), H₂⁺ e H₂ con lunghezza
di legame entro lo 0,1 % (H₂) e lo 0,4 % (H₂⁺), D₀ entro lo 0,25 %, vibrazioni entro l'1,5 %, elio
−78,994 ± 0,013 eV (misura −79,005).
L'acqua **non** è stata fatta: spiego perché al §8.

---

## 2. Gli unici ingressi

| ingresso | valore (CODATA 2018) | a cosa serve |
|---|---|---|
| ħ | 1,054571817·10⁻³⁴ J s | |
| e | 1,602176634·10⁻¹⁹ C | |
| mₑ | 9,1093837015·10⁻³¹ kg | |
| ε₀ | 8,8541878128·10⁻¹² F/m | |
| c | 299 792 458 m/s | solo per trasformare energie in lunghezze d'onda e per la velocità di emissione |
| carica dei nuclei | Z = 1 (protone), 2 (elio) | |
| massa dei nuclei | protone 1,67262192369·10⁻²⁷ kg, particella α 6,6446573357·10⁻²⁷ kg | vibrazioni, massa ridotta |

Da questi il programma **ricava** (`costanti.rs`) le unità atomiche: raggio di Bohr
a₀ = 4πε₀ħ²/(mₑe²) = 0,5291772103 Å, hartree Eₕ = ħ²/(mₑa₀²) = 27,21138628 eV, il tempo atomico
ħ/Eₕ = 2,42·10⁻¹⁷ s, la costante di struttura fine α = e²/(4πε₀ħc) = 1/137,036. Non c'è nessuna tabella
di proprietà di atomi o molecole. Le misure (NIST, Huber & Herzberg) servono **solo** per il confronto
finale.

---

## 3. La regola, scritta come "cosa fa ogni pallina"

### 3.1 Che cos'è una pallina

Una pallina è **un punto nello spazio di tutti gli elettroni** del sistema. Per l'idrogeno (1 elettrone)
è un punto in 3-D. Per H₂ o He (2 elettroni) è **una coppia** di punti che si muovono insieme: 6 numeri.
Per l'acqua sarebbe un punto in 30 dimensioni (10 elettroni).

### 3.2 Cosa fa a ogni tick (durata τ)

1. **Tremolio.** Fa un passo a caso, gaussiano, con varianza 2Dτ per coordinata, dove
   **D = ħ/2mₑ** (in unità atomiche D = ½). È la stessa diffusione dell'equazione di Schrödinger
   "in tempo immaginario".
2. **Vita.** Calcola l'energia potenziale V nel punto dove si trova (attrazione −Ze²/4πε₀r di ogni
   nucleo, repulsione e²/4πε₀r tra gli elettroni, repulsione tra i nuclei) e diventa
   ⌊w + u⌋ palline (0, 1, 2, …) con **w = exp(−τ (V − E_ref))** e u un numero a caso tra 0 e 1.
   Dove l'energia è bassa (vicino ai nuclei) le palline si moltiplicano, dove è alta muoiono.
3. **E_ref** è un solo numero globale che si aggiusta per tenere costante il numero di palline.
   Quando la nuvola si è stabilizzata, **E_ref è l'energia dello stato fondamentale.**

Questo è il metodo che i fisici chiamano *diffusion Monte Carlo*: qui lo scriviamo come regola per pallina.

### 3.3 La regola dei nodi (stati eccitati e fermioni)

Uno stato eccitato ha dei **nodi**: superfici dove l'onda vale zero e cambia segno. Regola: **una pallina
che attraversa un nodo muore**; se gli passa vicino muore con la probabilità che il suo cammino a caso
l'abbia toccato, exp(−2dd′/τ) (d, d′ = distanze dal nodo prima e dopo il passo). Così le palline restano
in una "tasca" e trovano lo stato più basso **con quel nodo**.

* **Nodi dalla simmetria** (esatti): lo stato 2p ha il piano z = 0 (è dispari per lo specchio z → −z);
  lo stato 3d_xy ha i due piani x = 0 e y = 0.
* **Nodi radiali** (2s, 3p): sono sfere, e il loro raggio **non** lo diamo noi: lo trovano le palline.
  Regola: per un nodo vero le due tasche (dentro e fuori dalla sfera) hanno **la stessa energia**
  (ogni tasca di uno stato vero ha la sua stessa energia). Facciamo vivere le palline separatamente
  dentro e fuori per tre raggi, interpoliamo le due energie e prendiamo il raggio dove si incrociano.
* **Fermioni.** Con 2 elettroni di spin opposto (H₂, He) lo stato fondamentale non ha nodi:
  **nessun problema di segno, il risultato è esatto** (a parte rumore e tick finito). Con più elettroni
  dello stesso spin (l'acqua ne ha 5 per spin) servirebbe la regola del nodo fisso con nodi costruiti
  dagli orbitali del livello sotto: è un'**approssimazione** (vedi §8).

### 3.4 La guida (seconda generazione di palline)

Le palline "nude" del §3.2 sono rumorose. Il livello sotto ci dà un aiuto: dalla nuvola di una prima
generazione di palline nude misuriamo la forma dell'orbitale (per l'idrogeno: la densità delle palline
decade come e^(−a r) con **a = 0,9987 /a₀**, misurato, contro 1 esatto). Una seconda generazione usa
questa forma ψ_g come **guida**:

* oltre al tremolio, ogni pallina **scivola** di τ·∇ln ψ_g verso dove la guida è più grande;
* per la vita usa l'energia locale E_loc = (Hψ_g)/ψ_g invece di V (più liscia, molto meno rumore).

La guida **non cambia la risposta** (salvo attraverso i suoi nodi): cambia solo il rumore. Con la guida la
densità delle palline diventa ψ_g·ψ ≈ |ψ|², cioè proprio la **densità elettronica**.
Le guide usate: l'orbitale 1s misurato dalle palline (per un nucleo di carica Z: lo stesso scalato di Z,
che è una proprietà matematica esatta del potenziale di Coulomb); per gli stati eccitati il fattore
angolare dato dalla simmetria (z, xy) per e^(−a r) con a misurato sulla nuvola nuda; per le molecole la
somma degli 1s sui due protoni, con il decadimento fissato dalla **condizione di cuspide**
(∂lnψ/∂r = −Z sul nucleo, conseguenza della legge di Coulomb); con due elettroni, il fattore
exp(r₁₂/2/(1+r₁₂)) che la stessa legge impone quando due elettroni si toccano.

### 3.5 Una regola in più per le molecole: accetta o rifiuta il passo

Per H₂⁺, H₂ e He le palline nude o solo guidate hanno un errore di tick troppo grande (per H₂⁺ a
τ = 0,01: −0,60120 invece di −0,60263 Eₕ). Aggiungo una regola **per pallina** (Reynolds et al. 1982):
dopo il passo la pallina confronta la guida nel punto nuovo e in quello vecchio e **tiene il passo** con
probabilità min(1, ψ_g²(nuovo)·G(vecchio←nuovo) / ψ_g²(vecchio)·G(nuovo←vecchio)), altrimenti resta
ferma. È locale, usa solo la guida, e toglie quasi tutto l'errore di tick (H₂⁺ a τ = 0,02: −0,60245,
a τ = 0,01: −0,60250, esatto −0,60263). Per l'idrogeno e i suoi stati eccitati **non** è usata.

### 3.6 Il surfer

Per uno stato stazionario **reale** (come tutti questi) la fase dell'onda è costante, quindi la velocità
di Bohm ∇S/m è **zero**: il surfer (la pallina di Bohm) **sta fermo** nel punto in cui è nato, scelto a
caso dalla nuvola. Nel visore è la stella bianca. Con due elettroni il surfer è **una configurazione
intera** (due stelle unite da un filo): in questo metodo non c'è un surfer per elettrone, perché la
pallina stessa vive nello spazio di entrambi.

---

## 4. Livello 1: l'atomo di idrogeno

`--misura idrogeno` (≈ 80–180 s). 20 000 palline.

| misura | calcolato | riferimento | errore |
|---|---|---|---|
| E, palline nude, τ = 0,02 / 0,01 / 0,005 | −0,49994 / −0,50054 / −0,49947 Eₕ (± 0,0007) | | |
| E, palline nude, τ → 0 | **−0,49966 ± 0,00077 Eₕ = −13,596 ± 0,021 eV** | −13,6057 eV (Schrödinger esatto) | −0,07 % (0,5 σ) |
| decadimento della nuvola nuda (= ψ) | e^(−0,9987 r/a₀) | e^(−r/a₀) | −0,13 % |
| E, palline guidate, τ → 0 | **−0,4999992 ± 0,0000009 Eₕ = −13,60567 eV** | −13,605693 eV | 0,02 meV |
| decadimento della nuvola guidata (= \|ψ\|²) | e^(−1,992 r/a₀) | e^(−2r/a₀) | −0,4 % |
| profilo radiale guidato vs e^(−2r)/π, 0,2 < r < 4 a₀ | scarto massimo 2,7 % | | |

Vicino al protone (r < 0,1 a₀) la nuvola è un po' più bassa del giusto (−8 % nel primo intervallo a
τ = 0,005, −20 % a τ = 0,01): è l'errore del tick finito e sparisce con τ → 0.
Con la massa finita del protone (massa ridotta, ingresso permesso) E = −13,5983 eV; l'energia di
ionizzazione misurata è 13,59843 eV (NIST): il resto (≈ 0,03 meV) è relatività e QED, che la nostra regola
non contiene.

---

## 5. Livelli eccitati e la luce emessa

`--misura eccitati` (2p, 3d, ≈ 150 s), `--misura nodo2s` (≈ 500 s), `--misura nodo3p` (≈ 460 s).

### 5.1 I livelli

| stato | nodi | E palline nude | E palline guidate (τ → 0) | esatto −1/(2n²) |
|---|---|---|---|---|
| 1s | nessuno | −0,49966 ± 0,00077 | −0,4999992 ± 0,0000009 | −0,5 |
| 2p | piano z = 0 (simmetria) | −0,12487 ± 0,00028 | −0,1249895 ± 0,0000011 | −0,125 |
| 2s | sfera r = **1,9970 ± 0,0023 a₀** (trovata dalle palline; esatta 2) | | −0,12512 ± 0,00008 | −0,125 |
| 3d | piani x = 0, y = 0 (simmetria) | −0,05542 ± 0,00022 | −0,0555520 ± 0,0000006 | −0,0555556 |
| 3p | piano z = 0 + sfera r = **5,997 ± 0,014 a₀** (esatta 6) | | −0,05558 ± 0,00008 | −0,0555556 |

**La degenerazione emerge da sola**: E(2s) − E(2p) = −0,00013 ± 0,00008 Eₕ, E(3p) − E(3d) =
−0,00003 ± 0,00008 Eₕ (zero entro gli errori). Nessuno ha detto alle palline che i livelli dipendono solo
da n.

Errori onesti: le barre d'errore statistiche sottostimano un piccolo errore sistematico del tick
(2p: 1,05·10⁻⁵ Eₕ sopra l'esatto, 9 σ; 3d: 3,6·10⁻⁶, 6 σ). Per 2s e 3p l'errore è dominato da dove sta
il nodo.

### 5.2 I fotoni emessi (λ = hc/ΔE, nel vuoto, con la massa ridotta)

| riga | transizione | calcolata | NIST (vuoto) | errore | colore |
|---|---|---|---|---|---|
| Lyman-α | 2p → 1s | **121,5653 ± 0,0004 nm** | 121,567 nm | −0,0017 nm (−1,4·10⁻⁵) | ultravioletto, invisibile |
| Lyman-β | 3p → 1s | **102,580 ± 0,017 nm** | 102,572 nm | +0,008 nm | ultravioletto |
| Balmer-α | 3d → 2p | **656,536 ± 0,012 nm** | 656,461 nm | +0,075 nm (+1,1·10⁻⁴) | rosso |

Il rosso dell'idrogeno (la riga Hα delle nebulose) esce dalle palline al decimo di nanometro. Lo scarto di
Balmer-α (6 σ) è il piccolo errore sistematico di tick del 2p (§5.1). La struttura fine (≈ 0,016 nm su
Hα) è relatività: la regola non la contiene.

---

## 6. Regole per la luce dell'atomo di idrogeno

Quanti fotoni manda indietro un atomo di idrogeno illuminato. Questi numeri sono pensati per il livello
"formazione delle immagini" (`docs/IMMAGINI.md` §5, che finora usava la soluzione esatta come segnaposto).
`--misura campo` (≈ 250 s) e `--misura riga` (≈ 450 s).

### 6.1 Polarizzabilità statica (come si deforma la nuvola in un campo elettrico)

Alle palline si aggiunge solo l'energia dell'elettrone nel campo uniforme F (eFz) e si misura
E(F) = E(0) − ½αF² − γF⁴/24 a F = 0; 0,01; 0,02; 0,03; 0,04 a.u. (1 a.u. = 5,14·10¹¹ V/m), poi un fit
con E(0), α, γ liberi (cioè F → 0).

| F (a.u.) | 0 | 0,01 | 0,02 | 0,03 | 0,04 |
|---|---|---|---|---|---|
| E (Eₕ) | −0,4999828(8) | −0,500185(16) | −0,500916(21) | −0,502178(53) | −0,503823(60) |

**α = 4,60 ± 0,12 a₀³** (esatta 4,5; +2,3 %, 0,8 σ). Il dipolo indotto letto direttamente dallo spostamento
della nuvola (2⟨z⟩/F, stima "estrapolata") dà 4,1–4,9: è più rumoroso e dipende da F; uso quello
dall'energia.

### 6.2 Forza della riga Lyman-α

Elemento di matrice ⟨1s|z|2p₀⟩ dalle due nuvole guidate (1s e 2p). Come: la densità delle palline guidate
è ψ_g·ψ, e ψ/ψ_g dipende solo da r (la guida porta già il fattore angolare esatto, imposto dalla
simmetria); quindi ogni nuvola **misura la correzione radiale** R(r) = ψ/ψ_g contando le palline in gusci
sferici: R ∝ (palline nel guscio)/∫ψ_g² sul guscio. Poi ψ = ψ_g R, normalizzazione e integrale con gli
angoli fatti esattamente. Due coppie di tick (τ₁ₛ, τ₂ₚ) = (0,01; 0,02) e (0,005; 0,01), poi τ → 0:
0,7315 ± 0,0031 e 0,7427 ± 0,0009 → **⟨1s|z|2p₀⟩ = 0,754 ± 0,004 a₀** (esatto 0,7449; +1,2 %, 2,4 σ).
È la parte più debole: la densità vicino al nucleo risente del tick e l'estrapolazione lineare in τ
probabilmente esagera un po' (il valore al tick più piccolo, 0,7427, è a −0,3 %).

| grandezza | palline | esatto | errore |
|---|---|---|---|
| ⟨1s\|z\|2p₀⟩ | 0,754 ± 0,004 a₀ | 0,7449 a₀ | +1,2 % |
| forza d'oscillatore f(1s→2p) = 2ΔE\|z\|² | 0,427 ± 0,004 | 0,4162 | +2,5 % |
| A(2p→1s) = 4ω³\|d\|²/3c³ | (6,42 ± 0,06)·10⁸ s⁻¹ | 6,2649·10⁸ s⁻¹ | +2,4 % |
| vita media del 2p | 1,558 ns | 1,596 ns | −2,4 % |

### 6.3 Le tre regole per la luce (da usare nel livello immagini)

| regola | calcolata (palline) | esatta | errore |
|---|---|---|---|
| sezione d'urto di Rayleigh a 500 nm, σ = (8π/3) k⁴ α², α come volume (α·a₀³) | **9,72·10⁻³² m²** (con α statica) | 9,29·10⁻³² m² (α = 4,5 statica) | +4,6 % |
| sezione d'urto risonante Lyman-α, 3λ²/2π | **7,056·10⁻¹⁵ m²** | 7,056·10⁻¹⁵ m² | < 10⁻⁴ |
| massimo ritmo di diffusione risonante, Γ/2 = A/2 | **3,21·10⁸ fotoni/s** | 3,13·10⁸ fotoni/s | +2,4 % |

Nota per chi usa la sezione di Rayleigh: a 500 nm la polarizzabilità **dinamica** è più grande di quella
statica. Con la sola riga Lyman-α (la nostra f e il nostro ΔE) α(ω) = α(0) + f[1/(ΔE² − ω²) − 1/ΔE²]
cresce di almeno +4 %; il valore esatto usato in `IMMAGINI.md` è 4,741 a₀³ (+5,4 %), che porta σ a
1,03·10⁻³¹ m². Le nostre palline per ora danno solo α statica.

---

## 7. Livello 2: le prime molecole, H₂⁺ e H₂

`--misura h2piu`, `--misura h2` (riprendono dai punti già salvati; ≈ 10 + 30 min con la macchina condivisa).
Per ogni distanza R tra i protoni (14 valori per H₂⁺, 18 per H₂) le palline danno E(R) (10 000 palline,
τ = 0,02, 300 / 200 unità di tempo per punto). La curva è interpolata con E(R) = E∞ + Σₖ₌₁..₅ cₖ e^(−kβR) (E∞ = gli atomi separati,
calcolati dalle palline; β scelto dal χ²). Poi i nuclei: equazione di Schrödinger radiale (J = 0) sulla
curva, con la massa dei nuclei (+ elettroni) → livelli vibrazionali, punto zero, ν(0→1), D₀ = E∞ − E(v=0).
Le barre d'errore vengono da 200 ri-fit con i punti spostati del loro errore.

### 7.1 H₂⁺ (due protoni, un elettrone)

| grandezza | palline | misura | errore | fonte |
|---|---|---|---|---|
| distanza di legame Rₑ | 1,0557 ± 0,0024 Å | 1,052 Å | +0,35 % | Huber & Herzberg 1979 |
| ωₑ (armonica) | 2302 ± 30 cm⁻¹ | 2321,7 cm⁻¹ | −0,9 % | idem |
| ν(0→1) | 2170 ± 23 cm⁻¹ | 2191,1 cm⁻¹ | −0,9 % | idem (ΔG½) |
| D₀ | 2,6443 ± 0,0045 eV | 2,6507 eV | −0,24 % | idem |
| Dₑ | 2,7856 ± 0,0060 eV | | | |

(La teoria esatta a nuclei fermi dà Rₑ = 1,0569 Å: il valore "sperimentale" 1,052 Å di Huber & Herzberg
è più vecchio e meno preciso.) Controllo del tick a R = 2 a₀: −0,602446 ± 0,000164 (τ = 0,02) e
−0,602500 ± 0,000135 (τ = 0,01; esatto −0,602634). χ²/punto del fit = 3,0: le barre d'errore dei singoli
punti sono ottimiste, e nei ri-fit sono moltiplicate per 2,3 (errori in tabella già gonfiati).

### 7.2 H₂ (due protoni, due elettroni: esatto, nessun problema di segno)

| grandezza | palline | misura | errore | fonte |
|---|---|---|---|---|
| distanza di legame Rₑ | **0,7420 ± 0,0008 Å** | 0,74144 Å | +0,08 % | Huber & Herzberg 1979 |
| D₀ (con il punto zero) | **4,4721 ± 0,0020 eV** | 4,47813 eV | −0,14 % (3 σ) | Liu et al. 2009 |
| ωₑ (armonica, dalla curvatura) | **4339 ± 18 cm⁻¹** (fit locale: 4349) | 4401,21 cm⁻¹ | −1,4 % (3,5 σ) | Huber & Herzberg 1979 |
| ν(0→1) (fondamentale) | **4129 ± 16 cm⁻¹** | 4161,17 cm⁻¹ | −0,8 % | Dickenson et al. 2013 |
| punto zero (ZPE) | 2155 ± 9 cm⁻¹ | | | |
| Dₑ (fondo della curva) | 4,7392 ± 0,0026 eV | (teoria esatta a nuclei fermi: 4,7477 eV) | −0,18 % | |

La curva: 18 distanze da 0,9 a 5,5 a₀; vicino al minimo (1,2–1,6 a₀) tre corse indipendenti per punto.
Energia a R = 1,4 a₀: −1,174332 ± 0,000096 Eₕ (teoria esatta −1,174476; 1,5 σ). Controllo del tick a
1,4 a₀: τ = 0,01 dà −1,173991 ± 0,000133, compatibile. Il fit ha χ²/punto = 3,5: le barre d'errore dei
punti sono troppo piccole di un fattore ~2, quindi nei ri-fit le ho moltiplicate per 2,3 (gli errori in
tabella sono già "gonfiati"). Controllo indipendente della curvatura: una quartica locale sugli 11 punti
vicini al minimo dà ωₑ = 4349 cm⁻¹, d'accordo con il fit globale.

**Cosa non torna.** La curvatura (ωₑ, ν(0→1)) è un po' bassa, −1,4 % e −0,8 %, oltre le barre d'errore.
Non ho trovato la causa: il tick non la spiega (controllo a τ/2) e il fit locale dà lo stesso numero. È la
grandezza più difficile, perché è una differenza seconda di energie rumorose: ±0,2 a₀ attorno al minimo
l'energia cambia solo di ~0,003 Eₕ, e ogni punto ha ~1·10⁻⁴ Eₕ di rumore.

### 7.3 La luce delle molecole

Una molecola biatomica omonucleare non ha dipolo, quindi la riga di vibrazione si vede solo debolmente
(quadrupolo, Raman). Dalle palline:

| molecola | riga di vibrazione 0→1 dalle palline | misura | errore |
|---|---|---|---|
| H₂ | ν = 4129 ± 16 cm⁻¹ → **λ = 2,422 µm** (infrarosso) | 4161,17 cm⁻¹ → 2,4032 µm | +0,8 % in λ |
| H₂⁺ | ν = 2170 ± 23 cm⁻¹ → **λ = 4,61 µm** | 2191,1 cm⁻¹ → 4,564 µm | +1,0 % in λ |

Sono i fotoni che queste molecole assorbono o emettono quando cambiano stato di vibrazione: è il ponte
verso il livello del fotone per le molecole. (L'intensità di queste righe, che per H₂ è piccolissima,
non l'ho calcolata.)

---

## 8. Elio (2 elettroni) e perché l'acqua non c'è ancora

### 8.1 Elio

`--misura elio` (≈ 3 min). Una pallina = una coppia di punti (spin ↑ e ↓), guida e^(−2r₁)e^(−2r₂) (il
nostro 1s scalato per Z = 2, cuspide esatta) × fattore elettrone–elettrone; regola "accetta o rifiuta";
10 000 palline, 150 unità di tempo per tick.

| tick τ | 0,02 | 0,01 | 0,005 | τ → 0 |
|---|---|---|---|---|
| E (Eₕ) | −2,90453 ± 0,00046 | −2,90412 ± 0,00048 | −2,90362 ± 0,00035 | **−2,90337 ± 0,00047** |

| grandezza | palline | riferimento | errore |
|---|---|---|---|
| E, nucleo fermo | −79,005 ± 0,013 eV | −79,0145 eV (teoria esatta, −2,903724 Eₕ) | +0,01 eV (0,8 σ) |
| E con la massa del nucleo (×μ/mₑ) | **−78,994 ± 0,013 eV** | **−79,0052 eV** (NIST: 24,5874 + 54,4178 eV) | +0,011 eV (0,9 σ) |

Due elettroni di spin opposto: nessun nodo, nessun problema di segno: il risultato è esatto entro il
rumore. (La correzione di massa trascura la "polarizzazione di massa", ~0,002 eV, e la relatività.)

### 8.2 L'acqua: non fatta, e perché

Il compito diceva: H₂O solo se 1–4 passano. Passano, ma l'acqua **non** l'ho calcolata, per due motivi
onesti:

1. **Il problema del segno.** L'acqua ha 10 elettroni, 5 per spin. La funzione d'onda di 5 elettroni
   uguali cambia segno quando due di loro si scambiano: le palline, che sono solo positive, da sole
   cadrebbero nello stato "bosonico" (sbagliato, molto più basso). Serve la regola del nodo fisso con una
   superficie nodale in 30 dimensioni, costruita (come vuole il progetto) solo dagli orbitali del livello
   sotto: 1s, 2s, 2p dell'ossigeno idrogenoidi (cioè i nostri orbitali dell'idrogeno scalati con Z) e 1s
   degli idrogeni, antisimmetrizzati (determinante di Slater). Il risultato è **un'approssimazione**:
   le palline trovano lo stato più basso con *quel* nodo, che è sopra il vero. Con nodi
   da orbitali buoni l'errore sull'energia totale dell'acqua è dell'ordine di 0,5 eV; con nodi da orbitali
   idrogenoidi "nudi" sarà più grande, e non lo so stimare senza farlo.
2. **Il costo.** Gli elettroni 1s dell'ossigeno (Z = 8) hanno energie di ~ −30 Eₕ e si muovono 8 volte
   più in fretta: servono tick ~ 64 volte più piccoli di quelli dell'idrogeno. La geometria richiede
   differenze di energia di ~10⁻⁴ Eₕ su un totale di 76 Eₕ. Stima: ~10⁸–10⁹ passi-pallina per punto
   della superficie, ore per punto su questa macchina in Rust su 4 core, e servono ~20 punti
   (due distanze, un angolo). Non sta nel limite di 10 minuti per esecuzione senza un lavoro in più
   (campionamento correlato tra geometrie, determinanti aggiornati per un elettrone alla volta, ecc.).

---

## 9. Cosa non va / limiti

* **Barre d'errore.** Sono calcolate a blocchi (il più grande di 32/16/8 blocchi). Per gli stati estesi
  (n = 3) e per i punti delle curve delle molecole (χ²/punto ≈ 3) sono ottimiste di un fattore ~2; nei
  fit delle molecole le ho gonfiate di conseguenza, altrove no.
* **Vibrazione di H₂** −1,4 % (ωₑ), −0,8 % (ν(0→1)): oltre l'errore dichiarato, causa non trovata (§7.2).
* **Tick finito.** Senza la regola "accetta o rifiuta" l'errore del tick è il limite principale
  (estrapolato τ → 0 con una retta, che lascia ~10⁻⁵ Eₕ sistematici sul 2p).
* **Elemento di matrice 1s–2p** (+1,2 %) e quindi A e Γ/2 (+2,4 %): è la misura più debole, perché usa
  la forma della nuvola e non solo l'energia.
* **Approssimazioni fisiche** che la regola contiene per costruzione: nuclei fermi per gli elettroni
  (Born–Oppenheimer; poi la massa dei nuclei entra nelle vibrazioni), niente relatività, niente QED,
  niente struttura fine/iperfine.
* **Polarizzabilità**: solo statica.
* **Fermioni con più di 2 elettroni**: non fatti (vedi §8). Quando si faranno, il nodo fisso sarà
  un'approssimazione e l'errore andrà dichiarato contro le misure.

---

## 10. Vedere le palline: il binario `molecole`

`cargo run --release --bin molecole` (opzione `--sistema h2|h2piu|he|2p`). Schermo pulito: la nuvola,
i nuclei (sfere dorate), il surfer (stella bianca). **Tab** (o ☰) apre il pannello:

* scelta del sistema: H 1s, 2p, 2s, 3d, 3p, H₂⁺, H₂, He;
* le palline girano **dal vivo** con la stessa regola delle misure (8000 palline, guida + accetta/rifiuta);
  ogni punto è un elettrone di una pallina; con due elettroni azzurro = spin ↑, arancio = spin ↓
  (si possono nascondere separatamente);
* superfici viola = i nodi (piani dalla simmetria, sfere trovate dalle palline); per 2p e 3d il visore
  mostra una sola tasca (il lobo z > 0, o il quadrante x, y > 0): le altre sono la sua immagine allo
  specchio, con la stessa energia;
* tabella calcolato / misurato / differenza; per l'idrogeno le righe di luce con il loro colore
  (Balmer-α rosso, Lyman in UV: riquadro vuoto = invisibile) e le "regole per la luce";
* per H₂ e H₂⁺: vibrazione dei nuclei alla frequenza calcolata (rallentata ~10¹⁴ volte, ampiezza del
  punto zero), le palline seguono i nuclei;
* il pannello dice chiaramente che la nuvola è la **densità** elettronica (quello che misura la
  diffrazione X), non tanti elettroni, e che con più elettroni non c'è un surfer per elettrone.

Memoria: la nuvola viene ridisegnata sostituendo la mesh (mai modificandola), come `swap_mesh` del
fotone, perché Bevy 0.15 non libera una mesh modificata sul posto.

---

## 11. Il prossimo livello: due molecole (il legame a idrogeno tra due acque)

Che cosa servirà, in ordine:

1. **L'atomo di ossigeno** con le palline: 8 elettroni, 5 + 3 per spin → primo vero caso di nodo fisso.
   Nodi dagli orbitali idrogenoidi (i nostri, scalati con Z) e una prima stima di quanto sbagliano, contro
   le energie di ionizzazione misurate (NIST).
2. **L'acqua**: determinanti con orbitali del livello sotto, fattore di cuspide, mosse di un elettrone
   alla volta (per non ricalcolare tutto), campionamento correlato tra geometrie vicine (le stesse palline
   valutate a due geometrie: l'errore sulla *differenza* è molto più piccolo). Geometria contro 0,9578 Å e
   104,48°, frequenze contro 3657 / 1595 / 3756 cm⁻¹, dipolo contro 1,855 D.
3. **Due acque**: il legame a idrogeno vale ~5 kcal/mol ≈ 0,008 Eₕ = 0,2 eV su 152 Eₕ di energia totale.
   Gli errori di nodo si cancellano in buona parte nella differenza (stessi orbitali nei due casi), ma
   bisognerà dimostrarlo, contro la misura dell'energia del dimero (D₀ ≈ 13,2 kJ/mol, Rocher-Casterline
   et al. 2011).
4. **Potenza di calcolo**: le palline sono indipendenti tra un passo e l'altro → la GPU (come il fotone).
   Senza questo il dimero dell'acqua non sta in tempi ragionevoli su questa macchina.

---

## 12. Riproduzione

```
sh lab/molecole_misura.sh                        # tutte le misure (ognuna < 10 min), log in lab/out/
cargo test --bin molecole                        # 9 test (regola, algebra della guida, nodi, elettroni
                                                 # per pallina, livelli nucleari, misure vs esperimento)
cargo run --release --bin molecole               # il visore
```
