# La regola base del fotone: una regola fissa, due ruoli, in 2D e in 3D

Laboratorio: `lab/base_*.py` (Python) e `lab/base_kernel/` (le somme delle ondine scritte in Rust e chiamate
da Python: stessa aritmetica della versione numpy di riferimento, solo più veloce).
Tutti i numeri di questo documento vengono da esecuzioni reali (seed fissi, 4 thread, i5-10300H condiviso con
altri processi). Registri: `lab/out/base_run_all_v1.log` (regola v1), `lab/out/base_run_all_v2_batch1.log`
(regola v2, prima corsa completa); risultati per esperimento: `lab/out/base_<esperimento>_<versione>.json`
(quelli della prima corsa v2 in `lab/out/base_v2_batch1/`); tabelle generate con `python lab/base_table.py v2`.

Il punto di partenza è la regola "pilota" di `docs/REGOLA_FOTONE.md` (2D). Qui diventa **una sola regola con
parametri fissi**, con il **numero di palline-onda dato dalla fisica** (i gradi di libertà del campo, cioè il
numero di Shannon) ed è estesa al **3D vettoriale**. La regola è stata congelata (v1) prima di far girare gli
esperimenti; un esperimento (la macchia di Arago) è fallito, la regola è stata cambiata **per tutti** (v2) e
**tutti** gli esperimenti sono stati rifatti. Nessuna tolleranza è stata cambiata.

---

## 1. La regola (cosa fa ogni pallina a ogni passo)

Un fotone = **tante palline-onda** + **una pallina-fotone**.

### 1.1 Cosa fa una pallina-onda

1. **Nasce su un oggetto** (sorgente, apertura, polarizzatore, beam-splitter, specchio, lamina di fase).
   Ogni oggetto fa sempre la stessa cosa: prende l'onda che gli arriva e la **ri-emette come un ventaglio
   regolare** di palline:
   * **origini** su una griglia regolare sulla parte illuminata dell'oggetto (il baricentro della parte
     trasparente di ogni cella); ampiezza di un'origine `a = B(origine) × area della cella`, dove
     `B = T·ψ(lato da cui arriva la luce) + R·ψ(lato verso cui va)` è l'onda uscente sul bordo e T, R sono le
     matrici di Jones di trasmissione e di riflessione dell'oggetto (§3);
   * **direzioni** su una griglia regolare che copre **tutto** il semispazio in avanti (Θ = 90°).
2. **Vola dritta** alla velocità della luce; la sua fase cresce di 2π per ogni λ di strada.
3. **Muore** quando incontra un altro oggetto (lì viene ri-emessa da quell'oggetto) o una parte opaca.
4. **Quante sono** non si sceglie: è il numero di gradi di libertà del campo (§2).

### 1.2 Cosa fa la pallina-fotone a ogni passo

1. **Nasce** sulla sorgente in un punto scelto con probabilità |B|² (dove c'è la luce).
2. **Ascolta** le palline-onda che le passano accanto: quelle di ogni ventaglio il cui raggio
   (origine → fotone) arriva fino a lei senza essere fermato da un oggetto. In 3D gli oggetti sono piani
   perpendicolari al fascio: sente il ventaglio dell'ultimo oggetto attraversato. In 2D (dove ci sono specchi e
   beam-splitter a 45°) il raggio viene controllato contro ogni oggetto; un oggetto che viene inserito o tolto
   (esperimento di Wheeler) ferma i raggi e ri-emette palline solo mentre c'è.
3. **Somma le ondine esatte** delle origini che sente, nella propria posizione:
   * 2D: `G = (ik/2)·(R·n)·H1(kR)/R` (Rayleigh–Sommerfeld esatto in 2D) e il suo gradiente;
   * 3D, Rayleigh–Sommerfeld **vettoriale**: `E_x = −2 Σ B_x ∂_n g`, `E_y = −2 Σ B_y ∂_n g`,
     `E_z = 2 Σ (B_x ∂_x g + B_y ∂_y g)`, con `g = e^{ikR}/(4πR)` (quindi `−2∂_z g` è proprio la ∂/∂z
     dell'ondina sferica), e il gradiente di ogni componente. La componente E_z è quella che rende il campo
     trasverso (∇·E = 0, verificato a 5·10⁻¹⁶): è il modo coerente di fare la "proiezione trasversa" in 3D;
   * ogni origine pesa `w_ℓ(R)` secondo il suo livello (§2.2); i pesi sommano esattamente a 1.
4. **Calcola la corrente** `j = Σ_pol Im(ψ*∇ψ)` (2D: le due componenti di Jones; 3D: E_x, E_y, E_z).
5. **Si muove alla velocità della luce lungo j** (Dormand–Prince 5(4), errore ≤ 10⁻⁴ λ per passo,
   passo ≤ 2 λ). Se non sente nessuna onda va dritta.
6. **Quando tocca un oggetto:**
   * parte opaca → muore (fotone assorbito, tutto intero);
   * parte trasparente → passa e riparte a δ = 0,02 λ oltre l'oggetto;
   * oggetto che **assorbe** (pareti con fenditure, aperture, polarizzatori, la sorgente stessa) →
     moneta locale `p = min(1, Φ_out/Φ_in)`, Φ = flusso normale n·j letto nella propria posizione
     (Φ_in subito prima dell'oggetto, Φ_out a δ dopo);
   * **(v2)** se lo attraversa a meno di **λ/2 da un bordo**, non vola dritta: sul segmento lungo λ che
     attraversa il bordo nel suo punto legge il flusso che arriva (Φ_in) e quello che esce a δ (j_z, anche
     dalla parte in ombra), e **conserva il suo posto in fila nel flusso**: se prima del bordo aveva una
     frazione F del flusso entrante, dopo si mette dove il flusso uscente ha la stessa frazione F; la moneta
     è il rapporto dei due flussi totali su quel segmento. È ciò che fa la dinamica di Bohm in una dimensione
     (le traiettorie non si incrociano) e serve perché un'onda con bordi netti manda subito una parte del
     flusso nell'ombra geometrica (§5);
   * oggetto che **non assorbe** (beam-splitter, lamina di fase) → nessuna moneta: dall'altra parte segue la
     corrente e finisce in un'uscita sola;
   * specchio → dietro non c'è onda: la corrente la gira prima che lo tocchi (se lo tocca è un errore
     numerico; lo contiamo: 0–7 fotoni su 1500);
   * schermo → **clic** dove arriva.

L'unica casualità: dove nasce e una moneta per ogni oggetto che assorbe. Nessuna decisione globale.

---

## 2. Quante palline: i gradi di libertà del campo

### 2.1 Livello 0: le palline "propaganti" (numero di Shannon)

Un'onda che esce da un oggetto largo W e può andare in tutte le direzioni in avanti (Θ → 90°: serve nel campo
vicino, come già visto in REGOLA_FOTONE) ha frequenze spaziali trasverse |f| ≤ sinΘ/λ = 1/λ. I suoi gradi di
libertà indipendenti sono, per ogni dimensione trasversa,

    N = 2·W·sinΘ/λ          (numero di Shannon / prodotto spazio-banda)

Campionamento di Nyquist con **un solo fattore di sovracampionamento η**:

* origini: `ds = λ/(2η·sinΘ) = λ/(2η)`  → `η·N` origini per dimensione;
* direzioni: `d(sinθ) = λ/(η·W)`       → `η·N` direzioni per dimensione (3D: griglia in (sinθ_x, sinθ_y)
  dentro il cerchio unitario, ≈ π/4·(2ηW_x/λ)(2ηW_y/λ) direzioni per origine).

**η = 2** fissato prima delle prove: la pallina-fotone non legge il campo ma il prodotto "campo × ondina",
la cui banda è la somma delle due (1/λ + 1/λ); Nyquist per il prodotto dà `ds = λ/4`, cioè η = 2.
Con η < 1 le copie (alias) della griglia diventano onde vere (lobi di reticolo).

**W** è l'estensione della parte dell'oggetto dove |B| ≥ ε·max|B|, **ε = 10⁻³** (sotto questa ampiezza una
componente si considera assente). W è un risultato del calcolo, non un input.

### 2.2 Livelli ℓ ≥ 1: le palline "evanescenti" (gradi di libertà vicino all'oggetto)

Vicino a un bordo netto il campo ha anche componenti evanescenti (|f| > 1/λ) che si spengono come
e^{−2πfz}. A distanza R da un'origine sono ancora sopra ε solo le frequenze fino a `f(R) = ln(1/ε)/(2πR)`.
Lo stesso Nyquist (stesso η) dà quindi, a chi legge a distanza R, `ds(R) = 1/(2η·f(R))`. Si realizza con
livelli di palline `ds_ℓ = ds₀/2^ℓ` che si sentono solo fino a `R < R₀·2^{2−ℓ}`, `R₀ = ln(1/ε)λ/(2π) = 1,10 λ`
(pesi "smoothstep" in log R tra livelli vicini, somma esattamente 1). Le palline fini sono **a vita breve**:
è il modo in cui la fisica dice che le onde evanescenti hanno corto raggio. Il livello più fine serve alla
distanza più piccola a cui una pallina-fotone legge un ventaglio, δ: `ℓ_max = ⌈log₂(R₀/δ)⌉ + 1 = 7`.
Il valore B delle palline fini si ottiene da quello che portano le palline in arrivo (campionato a λ/16 e
interpolato, Catmull–Rom) × la maschera dell'oggetto, calcolata esattamente sulla cella.

### 2.3 Parametri congelati

| parametro | valore | da dove viene |
|---|---|---|
| Θ (apertura dei ventagli) | 90° | tutto il semispazio (campo vicino) |
| η (sovracampionamento di Nyquist) | 2 | banda del prodotto campo × ondina (§2.1) |
| ε (soglia di ampiezza) | 10⁻³ | definisce W e la portata delle palline evanescenti |
| ds₀, d(sinθ) | λ/4, λ/(2W) | Nyquist con η = 2 |
| R₀, ℓ_max | 1,10 λ, 7 | ln(1/ε)λ/2π; ⌈log₂(R₀/δ)⌉+1 |
| δ (passaggio di un oggetto, lettura di Φ_out) | 0,02 λ | il più piccolo pratico (v1, invariato) |
| zona di bordo (v2) | λ/2 per lato | larghezza di transizione di un gradino limitato alla banda propagante (1/λ) |
| integratore della pallina-fotone | Dormand–Prince 5(4), ≤ 10⁻⁴ λ per passo, passo ≤ 2 λ | |

Il **numero di palline-onda per fotone** è un **risultato** (colonna della tabella §4):
`Σ_oggetti (origini × direzioni)`. Diamo quelle del livello 0 (i gradi di libertà propaganti) e il totale con
i livelli evanescenti (`N₀·Σ_ℓ 2^{ℓ·d}` origini, d = 1 in 2D, 2 in 3D). Il totale è enorme (fino a 10¹²) ma
quasi tutte le palline fini non vengono mai sentite: la pallina-fotone, vicino a un bordo, ne legge qualche
migliaio per passo (in Rust si creano solo quando servono e si riusano per tutti i fotoni).
Le direzioni, con l'ondina esatta dall'origine (stimatore "E2" di REGOLA_FOTONE) e il cappello normalizzato,
**non entrano nella guida**: dicono solo chi sente chi. Il loro numero è quello di Nyquist, riportato per
completezza; per costruzione la guida non cambia se le si infittisce o le si dirada (finché coprono il
semispazio), come già trovato in REGOLA_FOTONE (variante H).

### 2.4 Prova di convergenza (il numero di palline non è arbitrario)

Script `lab/base_conv.py` (regola v2, cambiando **un solo** numero alla volta).

**Onda** (deterministico): errore delle palline rispetto all'integrale di Rayleigh–Sommerfeld esatto
(quadratura fine indipendente), doppia fenditura 2D (max|U−U_RS|/max|U_RS|, a distanza z dalle fenditure) e
apertura circolare 3D D = 6 (max errore relativo di |E_x|):

| η | origini nelle fenditure | z = 0,02 | z = 0,1 | z = 0,5 | z = 2 | z = 10 | z = 160 | 3D z = 0,05 | 3D z = 2 | 3D z = 50 |
|---|---|---|---|---|---|---|---|---|---|---|
| 0,5 (sotto Nyquist) | 4 | 2,9e-2 | 6,2e-2 | 3,0e-1 | **6,0e-1** | 2,1e-1 | 6,7e-2 | 1,3e-1 | **4,4e-1** | 1,3e-2 |
| 1 (Nyquist del campo) | 8 | 3,1e-3 | 7,0e-3 | 2,6e-2 | 4,3e-2 | 3,3e-2 | 1,6e-2 | 6,2e-3 | 1,0e-1 | 3,3e-3 |
| **2 (regola)** | 16 | 3,1e-4 | 8,5e-4 | 3,7e-3 | 9,4e-3 | 7,5e-3 | 3,9e-3 | 2,0e-3 | 2,6e-2 | 9,4e-4 |
| 4 | 32 | 7,0e-5 | 1,6e-4 | 8,3e-4 | 2,3e-3 | 1,8e-3 | 9,6e-4 | 8,4e-4 | 6,6e-3 | 2,9e-4 |

L'errore scende come 1/η² sopra Nyquist; sotto (η = 0,5) l'onda è sbagliata del 30–60%.

**Fotoni** (doppia fenditura 2D, 12 000 fotoni emessi, ~2 570 clic; campionatore ideale con lo stesso numero
di clic: correlazione 0,9937 ± 0,0024, errore frange rms 4,8%; frazione assorbita esatta 0,7806):

| caso | palline/fotone (liv. 0 · tutti) | correlazione | χ²/gdl | err. frange rms / max | assorbiti |
|---|---|---|---|---|---|
| η = 0,5 | 1 409 · 3,6·10⁵ | **0,571** | **44,8** | **68% / 106%** | 0,856 |
| η = 1 | 5 636 · 1,4·10⁶ | 0,9912 | 1,39 | 10,7% / 23,1% | 0,785 |
| **η = 2 (regola)** | 22 249 · 5,7·10⁶ | 0,9949 | 0,96 | 6,6% / 13,1% | 0,786 |
| η = 4 | 89 585 · 2,3·10⁷ | 0,9919 | 1,31 | 6,6% / 11,7% | 0,786 |
| ε = 10⁻² (η = 2) | 15 040 · 3,8·10⁶ | 0,9940 | 1,14 | 5,7% / 12,1% | 0,786 |
| ε = 10⁻⁴ (η = 2, ℓ_max = 8) | 29 540 · 1,5·10⁷ | 0,9944 | 1,01 | 5,2% / 9,4% | 0,785 |

Lettura onesta: con i fotoni, sopra Nyquist (η ≥ 2) i risultati non cambiano più (le differenze sono rumore
statistico); a η = 1 si vede un peggioramento (χ² 1,39, una frangia +23%), sotto Nyquist (η = 0,5) il disegno è
distrutto. Con 2 500 clic la statistica non basta a separare η = 1 da η = 2 in modo netto: lo fa l'onda
(tabella sopra, 4–10× meno errore). ε cambia poco tra 10⁻² e 10⁻⁴.

---

## 3. Gli oggetti (definiti una volta: matrici di Jones + geometria)

Ogni oggetto è un segmento (2D) o un piano ⊥ z (3D) con una **maschera** (dove agisce) e due matrici di Jones
2×2, **T** (trasmissione) e **R** (riflessione). Se riceve onda dai due lati emette un ventaglio per lato:
`B₊ = T·ψ₋ + R·ψ₊`, `B₋ = T·ψ₊ + R·ψ₋`. Base di Jones: 2D (s = E_y, p = H_y, ognuna obbedisce all'equazione di
Helmholtz scalare), 3D (E_x, E_y tangenti al piano; E_z la ricostruisce l'ondina vettoriale).

| oggetto | maschera | T | R | assorbe (moneta)? |
|---|---|---|---|---|
| sorgente | tutto o un'apertura | fascio gaussiano (w₀, posizione della vita; w₀ = ∞ onda piana) × Jones | – | sì |
| parete con fenditure, apertura circolare, disco opaco, semipiano, reticolo | intervalli / cerchio dentro / cerchio fuori / rettangoli | 𝟙 | – | sì (+ zona di bordo v2) |
| polarizzatore lineare ad angolo θ | tutto | `a a†`, a = (cosθ, sinθ) | – | sì: p = \|⟨a\|ψ⟩\|²/\|ψ\|² |
| polarizzatore su una fenditura (which-path) | la fenditura | `f f†` | – | sì |
| beam-splitter 50:50 a 45° | segmento | 𝟙/√2 | i𝟙/√2 (stesse su entrambe le facce) | no |
| specchio metallico a 45° | segmento | – | diag(−1, +1) | no |
| lamina di fase | segmento | e^{iφ}·𝟙 | – | no |
| ostacolo / "bomba" | segmento | – | – | opaco: chi lo tocca muore |
| schermo, rivelatore | segmento o piano | – | – | clic |
| piano passivo (solo misura) | – | – | – | registra il passaggio, non tocca nulla |

"Sorgente con apertura" (Airy, Arago, fenditura lunga 3D; fenditura singola, reticolo, Talbot, bordo in 2D):
la pallina-fotone nasce direttamente sull'apertura con probabilità |B|²; è come una sorgente che illumina
l'apertura, senza simulare i fotoni che muoiono sulla parte opaca (la loro frazione è calcolata e riportata).
Doppia fenditura e which-path invece hanno sorgente e parete separate, come in REGOLA_FOTONE.

---

## 4. Risultati (regola v2, quella finale)

Tolleranze scritte prima delle corse (`TOL` in ogni script), mai cambiate. "σ" = errore statistico
(bootstrap o binomiale). "Campionatore ideale" = lo stesso stimatore applicato a clic estratti dal disegno
esatto con lo stesso numero di clic: dice quanto rumore ci si deve aspettare e se lo stimatore stesso ha un bias.
"Palline-onda per fotone" = somma su tutti i ventagli, livello 0 · tutti i livelli. "Tempo" = tempo reale
dello script su 4 thread, macchina condivisa.

| # | esperimento / grandezza | misurato | teoria | tolleranza | esito | fotoni / clic | palline-onda per fotone (livello 0 · tutti i livelli) | tempo |
|---|---|---|---|---|---|---|---|---|
| 1 | 2-D free beam width at z=1 zR | 4.2438 | 4.2426 | max(3σ=0.015, 0.02) | **PASS** | 20000 / 19975 | 3.969 · 1e+06 | 19 s |
| 2 | 2-D free beam width at z=2 zR | 6.7123 | 6.7082 | max(3σ=0.015, 0.02) | **PASS** | 20000 / 19975 (idem) | idem | 19 s |
| 3 | 2-D free beam width at z=5 zR | 15.3098 | 15.2971 | max(3σ=0.015, 0.02) | **PASS** | 20000 / 19975 (idem) | idem | 19 s |
| 4 | 3-D free beam width at z=1 zR | 4.3102 | 4.2426 | max(3σ=0.018, 0.02) | **PASS** | 7000 / 6982 | 9.743.750 · 2.1e+11 | 182 s |
| 5 | 3-D free beam width at z=2 zR | 6.8208 | 6.7082 | max(3σ=0.018, 0.02) | **PASS** | 7000 / 6982 (idem) | idem | 182 s |
| 6 | 3-D free beam width at z=5 zR | 15.5621 | 15.2971 | max(3σ=0.018, 0.02) | **PASS** | 7000 / 6982 (idem) | idem | 182 s |
| 7 | 3-D free beam, same photons: expansion w(1 zR)/w(1λ) | 1.4154 | 1.4133 | ±0.02 rel. | **PASS** | 7000 / 6982 (idem) | idem | 182 s |
| 8 | 3-D free beam, same photons: expansion w(2 zR)/w(1λ) | 2.2398 | 2.2347 | ±0.02 rel. | **PASS** | 7000 / 6982 (idem) | idem | 182 s |
| 9 | 3-D free beam, same photons: expansion w(5 zR)/w(1λ) | 5.1103 | 5.0958 | ±0.02 rel. | **PASS** | 7000 / 6982 (idem) | idem | 182 s |
| 10 | Gouy phase on axis, ζ(+zR) − ζ(−zR) (read by a photon ball on the axis) | 1.5722 rad | π/2 = 1.5708 (paraxial); exact RS 1.5707 | ±0.03 | **PASS** | 0 / 0 | 193.057.128 · 4.2e+12 | 3 s |
| 11 | Gouy phase on axis, ζ(+2zR) − ζ(−2zR) | 2.2080 rad | 2·atan(2) = 2.2143 | ±0.03 | **PASS** | 0 / 0 (idem) | idem | 3 s |
| 12 | Young double slit: correlation with exact flux | 0.9962 | 1 (ideal sampler 0.9971±0.0011) | >= 0.98 | **PASS** | 27000 / 5844 | 22.249 · 5.7e+06 | 35 s |
| 13 | Young double slit: bright-fringe population error rms / max | 4.3% / 7.6% | 0 (ideal sampler rms 3.2%) | rms <= 8% | **PASS** | 27000 / 5844 (idem) | idem | 35 s |
| 14 | Young double slit: absorbed fraction | 0.7836 | 0.7806 | ±0.02 | **PASS** | 27000 / 5844 (idem) | idem | 35 s |
| 15 | Taylor build-up: correlation vs number of clicks | N=10: 0.166 (ideal 0.487±0.139); N=50: 0.818 (ideal 0.783±0.072); N=200: 0.942 (ideal 0.929±0.024); N=1000: 0.986 (ideal 0.984±0.006); N=3000: 0.996 (ideal 0.995±0.002); N=5844: 0.996 (ideal 0.997±0.001) | ideal sampler | within 3σ of the ideal sampler at every N | **PASS** | 27000 / 5844 (idem) | idem | 35 s |
| 16 | 3-D long double slit (x-marginal) vs exact 2-D flux: correlation | 0.9952 | 1 (ideal sampler 0.9971±0.0010); 2-D rule: see Young | >= 0.98 | **PASS** | 5895 / 5635 | 3.557.772 · 7.8e+10 | 137 s |
| 17 | 3-D long double slit: bright-fringe error rms / max | 3.9% / 6.6% | 0 (ideal rms 3.2%) | rms <= 8% | **PASS** | 5895 / 5635 (idem) | idem | 137 s |
| 18 | 3-D long double slit: absorbed fraction | 0.7723 | 0.7806 (exact 2-D RS) | ±0.02 | **PASS** | 5895 / 5635 (idem) | idem | 137 s |
| 19 | single slit a=4.0: minimum m=1 (mean of ±) | 104.17 ± 0.52 | 103.28 = L·tan(asin(mλ/a)) [paraxial mλL/a = 100.0; exact RS 103.28] | max(3σ, 2%) | **PASS** | 58172 / 55197 | 256 · 6.5e+04 | 88 s |
| 20 | single slit a=4.0: minimum m=2 (mean of ±) | 231.49 ± 1.62 | 230.94 = L·tan(asin(mλ/a)) [paraxial mλL/a = 200.0; exact RS 230.95] | max(3σ, 2%) | **PASS** | 58172 / 55197 (idem) | idem | 88 s |
| 21 | grating N=5, d=8.0: principal maximum m=1 (centroid, mean of ±) | 25.204 ± 0.011 | 25.198 = L·tan(asin(mλ/d)) [paraxial 25.00] | max(3σ, 1%) | **PASS** | 37681 / 34196 | 5.440 · 1.4e+06 | 67 s |
| 22 | grating N=5, d=8.0: principal maximum m=2 (centroid, mean of ±) | 51.612 ± 0.019 | 51.640 = L·tan(asin(mλ/d)) [paraxial 50.00] | max(3σ, 1%) | **PASS** | 37681 / 34196 (idem) | idem | 67 s |
| 23 | grating N=5: secondary maxima between order 0 and ±1 | 3 (+), 3 (−) | 3 | exact | **PASS** | 37681 / 34196 (idem) | idem | 67 s |
| 24 | Airy disk D=6.0: first zero radius | 41.63 ± 0.24 | 41.52 = L·tan(asin(1.22λ/D)) [paraxial 1.22λL/D = 40.67; exact RS 41.52] | max(3σ, 2%) | **PASS** | 9799 / 9356 | 219.252 · 4.8e+09 | 120 s |
| 25 | Arago spot (disk a=3.0, beam w=5.0, L=12.0): fraction of clicks within r < 0.77 | 0.0294 ± 0.0022 | 0.0293 (exact RS) | max(3σ, 10%) | **PASS** | 6829 / 6627 | 78.418.960 · 1.7e+12 | 181 s |
| 26 | Arago spot: density in the spot / density in the shadow ring | 3.3 | 4.1 (exact RS) | > 3.0 | **PASS** | 6829 / 6627 (idem) | idem | 181 s |
| 27 | Talbot self-image: z of best image (clicks) | 125.4 | 128.0 = 2d²/λ  [exact RS best plane 125.4] | ±2.6 | **PASS** | 18791 / 17111 | 181.440 · 4.6e+07 | 123 s |
| 28 | Talbot half-distance: the image is shifted by d/2 (contrast minimum near zT/2) | 61.4 | 64.0 | ±2.6 | **PASS** | 18791 / 17111 (idem) | idem | 123 s |
| 29 | Fresnel edge: I(0)/I0 at the geometric shadow edge | 0.2610 ± 0.0110 | 0.25 | max(3σ, 0.02) | **PASS** | 59724 / 59255 | 44.100 · 1.1e+07 | 126 s |
| 30 | Fresnel edge: first maximum position | 8.438 ± 0.219 | 8.748 (exact RS, this beam) [plane-wave Fresnel u=1.2172: 8.607] | max(3σ, 2%) | **PASS** | 59724 / 59255 (idem) | idem | 126 s |
| 31 | Fresnel edge: first maximum height I/I0 (window ±0.1·sqrt(λL/2)) | 1.4516 ± 0.0166 | 1.4465 (exact RS, this beam) [plane-wave Fresnel 1.3615] | max(3σ, 0.02) | **PASS** | 59724 / 59255 (idem) | idem | 126 s |
| 32 | Fresnel edge: first minimum position | 13.423 ± 0.097 | 13.347 (exact RS, this beam) [plane-wave Fresnel u=1.8725: 13.241] | max(3σ, 2%) | **PASS** | 59724 / 59255 (idem) | idem | 126 s |
| 33 | Fresnel edge: first minimum height I/I0 (window ±0.1·sqrt(λL/2)) | 0.6637 ± 0.0143 | 0.6893 (exact RS, this beam) [plane-wave Fresnel 0.7920] | max(3σ, 0.02) | **PASS** | 59724 / 59255 (idem) | idem | 126 s |
| 34 | Fresnel edge: second maximum position | 16.679 ± 0.137 | 16.720 (exact RS, this beam) [plane-wave Fresnel u=2.3445: 16.578] | max(3σ, 2%) | **PASS** | 59724 / 59255 (idem) | idem | 126 s |
| 35 | Fresnel edge: second maximum height I/I0 (window ±0.1·sqrt(λL/2)) | 1.3985 ± 0.0273 | 1.3755 (exact RS, this beam) [plane-wave Fresnel 1.1822] | max(3σ, 0.02) | **PASS** | 59724 / 59255 (idem) | idem | 126 s |
| 36 | which-path, parallel polarizers: fringe visibility | 0.920 ± 0.015 | 0.907 | > 0.85 | **PASS** | 60000 / 6502 | 21.801 · 5.6e+06 | 72 s |
| 37 | which-path, orthogonal polarizers: fringe visibility | -0.009 ± 0.042 | 0.029 | < 0.15 | **PASS** | 60000 / 6502 (idem) | idem | 72 s |
| 38 | Malus (absorbing polarizer, 3-D) θ=0°: P(pass) | 1.0000 | 1.0000 = cos²θ | max(3σ=0.0020, 0.01) | **PASS** | 7500 / 3727 | 4.757.193 · 1e+11 | 53 s |
| 39 | Malus (absorbing polarizer, 3-D) θ=22.5°: P(pass) | 0.8554 | 0.8536 = cos²θ | max(3σ=0.0274, 0.01) | **PASS** | 7500 / 3727 (idem) | idem | 53 s |
| 40 | Malus (absorbing polarizer, 3-D) θ=45°: P(pass) | 0.4900 | 0.5000 = cos²θ | max(3σ=0.0388, 0.01) | **PASS** | 7500 / 3727 (idem) | idem | 53 s |
| 41 | Malus (absorbing polarizer, 3-D) θ=67.5°: P(pass) | 0.1475 | 0.1464 = cos²θ | max(3σ=0.0275, 0.01) | **PASS** | 7500 / 3727 (idem) | idem | 53 s |
| 42 | Malus (absorbing polarizer, 3-D) θ=90°: P(pass) | 0.0000 | 0.0000 = cos²θ | max(3σ=0.0020, 0.01) | **PASS** | 7500 / 3727 (idem) | idem | 53 s |
| 43 | three polarizers 0/45/90°, light polarized at 0°: P(pass) | 0.2490 | 0.2500 | max(3σ=0.0206, 0.01) | **PASS** | 8000 / 1482 | 22.298.321 · 4.9e+11 | 130 s |
| 44 | two crossed polarizers 0/90° (the middle one removed): P(pass) | 0.0000 | 0 | 0 exactly | **PASS** | 8000 / 1482 (idem) | idem | 130 s |
| 45 | three polarizers 0/45/90°, unpolarized light: P(pass) | 0.1239 | 0.1250 | max(3σ=0.0157, 0.01) | **PASS** | 8000 / 1482 (idem) | idem | 130 s |
| 46 | Mach–Zehnder φ=0°: P(D) | 0.9967 | 1.0000 | max(3σ=0.0020, 0.02) | **PASS** | 3000 / 2983 | 157.564 · 4e+07 | 239 s |
| 47 | Mach–Zehnder φ=90°: P(D) | 0.4993 | 0.5000 | max(3σ=0.0388, 0.02) | **PASS** | 3000 / 2983 (idem) | idem | 239 s |
| 48 | Mach–Zehnder φ=180°: P(D) | 0.0000 | 0.0000 | max(3σ=0.0020, 0.02) | **PASS** | 3000 / 2980 | 157.564 · 4e+07 | 217 s |
| 49 | Mach–Zehnder φ=270°: P(D) | 0.5043 | 0.5000 | max(3σ=0.0387, 0.02) | **PASS** | 3000 / 2980 (idem) | idem | 217 s |
| 50 | bomb tester: P(explosion) | 0.5088 | 0.5000 | max(3σ=0.0274, 0.02) | **PASS** | 3000 / 1471 | 151.909 · 3.9e+07 | 110 s |
| 51 | bomb tester: P(C) = bomb detected without interaction (dark port) | 0.2341 | 0.2500 | max(3σ=0.0237, 0.02) | **PASS** | 3000 / 1471 (idem) | idem | 110 s |
| 52 | bomb tester: P(D) | 0.2571 | 0.2500 | max(3σ=0.0237, 0.02) | **PASS** | 3000 / 1471 (idem) | idem | 110 s |
| 53 | same interferometer without bomb: P(C) (dark port) | 0.0000 | 0.0000 | max(3σ=0.0030, 0.02) | **PASS** | 3000 / 1471 (idem) | idem | 110 s |
| 54 | Wheeler, BS2 inserted after the photon passed BS1: P(D) | 0.9946 | 1.0000 | max(3σ=0.0040, 0.02) | **PASS** | 1500 / 1490 | 157.564 · 4e+07 | 147 s |
| 55 | Wheeler, BS2 removed after the photon passed BS1: P(D) | 0.5053 | 0.5000 | max(3σ=0.0548, 0.02) | **PASS** | 1500 / 1490 (idem) | idem | 147 s |
| 56 | Wheeler, BS2 removed after the photon passed BS1: fraction of photons whose detector is the straight continuation of their arm (C⇔A, D⇔B) | 0.0000 | 0.0000 | max(3σ=0.0040, 0.02) | **PASS** | 1500 / 1490 (idem) | idem | 147 s |

Note:
* **Fascio libero 3D**: lo scarto di +1,6…1,7% sulla larghezza assoluta c'è già a z = 1 λ (3,045 invece di
  3,009 attesi dal flusso j_z, che per il campo vettoriale è più largo di |E|² dello 0,3%): è il campione di
  nascita (+1,2%, 2σ). Con il seed 12 i primi 4 000 fotoni sono gli stessi in tutte le corse (v1: +2,0% con
  4 000 fotoni; prima corsa v2: +1,0% con 10 000), per questo lo scarto si ripete. Con 100 000 nascite la
  larghezza a z = δ è 3,0079 contro 3,0084 del flusso esatto. L'allargamento delle stesse traiettorie
  (w(z)/w(1 λ)) segue la teoria entro 0,1–0,3%.
* **Gouy**: è una proprietà dell'onda letta da una pallina-fotone sull'asse (nessuna statistica). Le palline
  coincidono con il Rayleigh–Sommerfeld esatto a 1,7·10⁻³ rad; lo scarto fino a 0,048 rad dalla curva
  parassiale atan(z/zR) lontano dal fuoco è fisico (fascio non parassiale, w₀ = 2λ).
* **Talbot**: il miglior autoritratto è a 125,4 λ invece di 2d²/λ = 128: è la vera distanza per un reticolo
  con d = 8λ (il Rayleigh–Sommerfeld esatto dà lo stesso 125,4); la formula 2d²/λ è parassiale. Passa al
  limite della tolleranza (2%).
* **Bordo di Fresnel**: le altezze dei massimi si confrontano con il Rayleigh–Sommerfeld esatto **per lo stesso
  fascio gaussiano** (w = 20 λ = 2,8 scale di Fresnel); i numeri dell'onda piana (tra parentesi) non valgono
  per questo fascio (§5, A3). I(0)/I₀ = 1/4 vale in entrambi i casi.
* **Wheeler**: con BS2 tolto, le traiettorie non si incrociano dove i due fasci si incrociano (posto vuoto di
  BS2): un fotone del braccio A finisce nel rivelatore "di B" e viceversa (traiettorie "surrealistiche",
  Englert–Scully–Süssmann–Walther 1992). La meccanica quantistica non dice nulla su un cammino non osservato;
  le statistiche dei clic sono quelle quantistiche.
* **Bomba**: il fotone del braccio con la bomba muore se la tocca (esplosione); quello dell'altro braccio
  sente a BS2 solo la propria onda e si divide 50/50: 1/4 dei fotoni clicca nel rivelatore buio senza aver
  toccato la bomba.
* **Fenditura lunga 3D**: è la doppia fenditura 2D fatta con il motore 3D (fenditure lunghe 120 λ, fascio largo
  7 in x e 4 in y): la distribuzione in x dei clic coincide con quella 2D esatta (0,9952; rms 3,9%),
  e il nucleo 3D sommato lungo una riga di origini in y riproduce l'ondina 2D di Hankel a 10⁻⁴
  (`lab/base_test.py`).

---

## 5. Registro delle versioni della regola

### v1 — congelata il 2026-10-03 (η = 2, ε = 10⁻³, δ = 0,02 λ, DP45 10⁻⁴ λ)

Prima della corsa ufficiale v1 sono stati trovati e corretti (senza toccare regola né tolleranze):
* **A1** (analisi) fenditura singola: il minimo stimato con una parabola simmetrica era spostato (107 invece
  di 103,3). Lo stesso stimatore su clic ideali estratti dal disegno esatto dava 107,3: era lo stimatore, non la
  regola. Sostituito con una cubica pesata, verificata sul campionatore ideale (bias < 0,4σ).
* **A2** (analisi) reticolo: il contatore dei massimi secondari lisciava su 3,75 λ e ne vedeva 1 invece di 3,
  mentre l'istogramma dei clic li mostrava tutti e 3 al posto giusto. Sostituito con un contatore a
  prominenza (3σ) che sul campionatore ideale trova sempre 3.
* **A3** (progetto dell'esperimento) bordo di Fresnel: le altezze dei massimi erano confrontate con la formula
  dell'onda piana (1,370 / 0,778 / 1,199), ma il fascio usato è gaussiano con w = 20 λ, solo 2,8 scale di
  Fresnel: per questo fascio il Rayleigh–Sommerfeld esatto dà 1,447 / 0,689 / 1,376 e i clic danno
  1,454 / 0,670 / 1,400. La prima versione dava FAIL sulle altezze (file `base_edge_v1_firstdesign.json`);
  la teoria corretta per l'esperimento fatto è quella esatta per lo stesso fascio.
* **B1, B2** (programma) gli specchi non emettevano il ventaglio (maschera letta dalla trasmissione invece che
  dall'estensione dell'oggetto); gli istanti di inserimento/rimozione di BS2 erano indicizzati male. Trovati
  con prove piccole prima della corsa.

Risultati v1 (corsa completa, `lab/out/base_run_all_v1.log`):

| # | esperimento / grandezza | misurato | teoria | tolleranza | esito | fotoni / clic | palline-onda per fotone (livello 0 · tutti i livelli) | tempo |
|---|---|---|---|---|---|---|---|---|
| 1 | 2-D free beam width at z=1 zR | 4.2438 | 4.2426 | max(3σ=0.015, 0.02) | **PASS** | 20000 / 19975 | 3.969 · 1e+06 | 21 s |
| 2 | 2-D free beam width at z=2 zR | 6.7123 | 6.7082 | max(3σ=0.015, 0.02) | **PASS** | 20000 / 19975 (idem) | idem | 21 s |
| 3 | 2-D free beam width at z=5 zR | 15.3098 | 15.2971 | max(3σ=0.015, 0.02) | **PASS** | 20000 / 19975 (idem) | idem | 21 s |
| 4 | 3-D free beam width at z=1 zR | 4.3402 | 4.2426 | max(3σ=0.024, 0.02) | **PASS** | 4000 / 3987 | 9.743.750 · 2.1e+11 | 103 s |
| 5 | 3-D free beam width at z=2 zR | 6.8684 | 6.7082 | max(3σ=0.024, 0.02) | **FAIL** | 4000 / 3987 (idem) | idem | 103 s |
| 6 | 3-D free beam width at z=5 zR | 15.6709 | 15.2971 | max(3σ=0.024, 0.02) | **FAIL** | 4000 / 3987 (idem) | idem | 103 s |
| 7 | Gouy phase on axis, ζ(+zR) − ζ(−zR) (read by a photon ball on the axis) | 1.5722 rad | π/2 = 1.5708 (paraxial); exact RS 1.5707 | ±0.03 | **PASS** | 0 / 0 | 193.057.128 · 4.2e+12 | 3 s |
| 8 | Gouy phase on axis, ζ(+2zR) − ζ(−2zR) | 2.2080 rad | 2·atan(2) = 2.2143 | ±0.03 | **PASS** | 0 / 0 (idem) | idem | 3 s |
| 9 | Young double slit: correlation with exact flux | 0.9965 | 1 (ideal sampler 0.9973±0.0010) | >= 0.98 | **PASS** | 27000 / 5809 | 22.249 · 5.7e+06 | 34 s |
| 10 | Young double slit: bright-fringe population error rms / max | 4.0% / 7.1% | 0 (ideal sampler rms 3.2%) | rms <= 8% | **PASS** | 27000 / 5809 (idem) | idem | 34 s |
| 11 | Young double slit: absorbed fraction | 0.7849 | 0.7806 | ±0.02 | **PASS** | 27000 / 5809 (idem) | idem | 34 s |
| 12 | Taylor build-up: correlation vs number of clicks | N=10: 0.247 (ideal 0.487±0.139); N=50: 0.843 (ideal 0.783±0.072); N=200: 0.959 (ideal 0.929±0.024); N=1000: 0.989 (ideal 0.984±0.006); N=3000: 0.993 (ideal 0.995±0.002); N=5809: 0.997 (ideal 0.997±0.001) | ideal sampler | within 3σ of the ideal sampler at every N | **PASS** | 27000 / 5809 (idem) | idem | 34 s |
| 13 | 3-D long double slit (x-marginal) vs exact 2-D flux: correlation | 0.9907 | 1 (ideal sampler 0.9970±0.0010); 2-D rule: see Young | >= 0.98 | **PASS** | 6000 / 5551 | 3.557.772 · 7.8e+10 | 94 s |
| 14 | 3-D long double slit: bright-fringe error rms / max | 6.5% / 9.9% | 0 (ideal rms 3.2%) | rms <= 8% | **PASS** | 6000 / 5551 (idem) | idem | 94 s |
| 15 | 3-D long double slit: absorbed fraction | 0.7397 | 0.7806 (exact 2-D RS) | ±0.02 | **FAIL** | 6000 / 5551 (idem) | idem | 94 s |
| 16 | single slit a=4.0: minimum m=1 (mean of ±) | 103.91 ± 0.42 | 103.28 = L·tan(asin(mλ/a)) [paraxial mλL/a = 100.0; exact RS 103.28] | max(3σ, 2%) | **PASS** | 60000 / 56852 | 256 · 6.5e+04 | 83 s |
| 17 | single slit a=4.0: minimum m=2 (mean of ±) | 230.34 ± 1.47 | 230.94 = L·tan(asin(mλ/a)) [paraxial mλL/a = 200.0; exact RS 230.95] | max(3σ, 2%) | **PASS** | 60000 / 56852 (idem) | idem | 83 s |
| 18 | grating N=5, d=8.0: principal maximum m=1 (centroid, mean of ±) | 25.192 ± 0.010 | 25.198 = L·tan(asin(mλ/d)) [paraxial 25.00] | max(3σ, 1%) | **PASS** | 40000 / 36324 | 5.440 · 1.4e+06 | 53 s |
| 19 | grating N=5, d=8.0: principal maximum m=2 (centroid, mean of ±) | 51.616 ± 0.018 | 51.640 = L·tan(asin(mλ/d)) [paraxial 50.00] | max(3σ, 1%) | **PASS** | 40000 / 36324 (idem) | idem | 53 s |
| 20 | grating N=5: secondary maxima between order 0 and ±1 | 3 (+), 3 (−) | 3 | exact | **PASS** | 40000 / 36324 (idem) | idem | 53 s |
| 21 | Airy disk D=6.0: first zero radius | 41.76 ± 0.23 | 41.52 = L·tan(asin(1.22λ/D)) [paraxial 1.22λL/D = 40.67; exact RS 41.52] | max(3σ, 2%) | **PASS** | 9952 / 9459 | 219.252 · 4.8e+09 | 80 s |
| 22 | Arago spot (disk a=3.0, beam w=5.0, L=12.0): fraction of clicks within r < 0.77 | 0.0143 ± 0.0013 | 0.0293 (exact RS) | max(3σ, 10%) | **FAIL** | 9957 / 9489 | 78.418.960 · 1.7e+12 | 177 s |
| 23 | Arago spot: density in the spot / density in the shadow ring | 1.6 | 4.1 (exact RS) | > 3.0 | **FAIL** | 9957 / 9489 (idem) | idem | 177 s |
| 24 | Talbot self-image: z of best image (clicks) | 125.4 | 128.0 = 2d²/λ  [exact RS best plane 125.4] | ±2.6 | **PASS** | 20000 / 18117 | 181.440 · 4.6e+07 | 114 s |
| 25 | Talbot half-distance: the image is shifted by d/2 (contrast minimum near zT/2) | 61.4 | 64.0 | ±2.6 | **PASS** | 20000 / 18117 (idem) | idem | 114 s |
| 26 | Fresnel edge: I(0)/I0 at the geometric shadow edge | 0.2628 ± 0.0110 | 0.25 | max(3σ, 0.02) | **PASS** | 60000 / 59504 | 44.100 · 1.1e+07 | 125 s |
| 27 | Fresnel edge: first maximum position | 8.550 ± 0.232 | 8.748 (exact RS, this beam) [plane-wave Fresnel u=1.2172: 8.607] | max(3σ, 2%) | **PASS** | 60000 / 59504 (idem) | idem | 125 s |
| 28 | Fresnel edge: first maximum height I/I0 (window ±0.1·sqrt(λL/2)) | 1.4544 ± 0.0163 | 1.4465 (exact RS, this beam) [plane-wave Fresnel 1.3615] | max(3σ, 0.02) | **PASS** | 60000 / 59504 (idem) | idem | 125 s |
| 29 | Fresnel edge: first minimum position | 13.403 ± 0.098 | 13.347 (exact RS, this beam) [plane-wave Fresnel u=1.8725: 13.241] | max(3σ, 2%) | **PASS** | 60000 / 59504 (idem) | idem | 125 s |
| 30 | Fresnel edge: first minimum height I/I0 (window ±0.1·sqrt(λL/2)) | 0.6702 ± 0.0139 | 0.6893 (exact RS, this beam) [plane-wave Fresnel 0.7920] | max(3σ, 0.02) | **PASS** | 60000 / 59504 (idem) | idem | 125 s |
| 31 | Fresnel edge: second maximum position | 16.794 ± 0.128 | 16.720 (exact RS, this beam) [plane-wave Fresnel u=2.3445: 16.578] | max(3σ, 2%) | **PASS** | 60000 / 59504 (idem) | idem | 125 s |
| 32 | Fresnel edge: second maximum height I/I0 (window ±0.1·sqrt(λL/2)) | 1.4002 ± 0.0279 | 1.3755 (exact RS, this beam) [plane-wave Fresnel 1.1822] | max(3σ, 0.02) | **PASS** | 60000 / 59504 (idem) | idem | 125 s |
| 33 | which-path, parallel polarizers: fringe visibility | 0.924 ± 0.014 | 0.907 | > 0.85 | **PASS** | 60000 / 6513 | 21.801 · 5.6e+06 | 80 s |
| 34 | which-path, orthogonal polarizers: fringe visibility | 0.000 ± 0.045 | 0.029 | < 0.15 | **PASS** | 60000 / 6513 (idem) | idem | 80 s |
| 35 | Malus (absorbing polarizer, 3-D) θ=0°: P(pass) | 1.0000 | 1.0000 = cos²θ | max(3σ=0.0020, 0.01) | **PASS** | 7500 / 3727 | 4.757.193 · 1e+11 | 51 s |
| 36 | Malus (absorbing polarizer, 3-D) θ=22.5°: P(pass) | 0.8554 | 0.8536 = cos²θ | max(3σ=0.0274, 0.01) | **PASS** | 7500 / 3727 (idem) | idem | 51 s |
| 37 | Malus (absorbing polarizer, 3-D) θ=45°: P(pass) | 0.4900 | 0.5000 = cos²θ | max(3σ=0.0388, 0.01) | **PASS** | 7500 / 3727 (idem) | idem | 51 s |
| 38 | Malus (absorbing polarizer, 3-D) θ=67.5°: P(pass) | 0.1475 | 0.1464 = cos²θ | max(3σ=0.0275, 0.01) | **PASS** | 7500 / 3727 (idem) | idem | 51 s |
| 39 | Malus (absorbing polarizer, 3-D) θ=90°: P(pass) | 0.0000 | 0.0000 = cos²θ | max(3σ=0.0020, 0.01) | **PASS** | 7500 / 3727 (idem) | idem | 51 s |
| 40 | three polarizers 0/45/90°, light polarized at 0°: P(pass) | 0.2490 | 0.2500 | max(3σ=0.0206, 0.01) | **PASS** | 8000 / 1482 | 22.298.321 · 4.9e+11 | 126 s |
| 41 | two crossed polarizers 0/90° (the middle one removed): P(pass) | 0.0000 | 0 | 0 exactly | **PASS** | 8000 / 1482 (idem) | idem | 126 s |
| 42 | three polarizers 0/45/90°, unpolarized light: P(pass) | 0.1239 | 0.1250 | max(3σ=0.0157, 0.01) | **PASS** | 8000 / 1482 (idem) | idem | 126 s |
| 43 | Mach–Zehnder φ=0°: P(D) | 0.9967 | 1.0000 | max(3σ=0.0020, 0.02) | **PASS** | 3000 / 2983 | 157.564 · 4e+07 | 262 s |
| 44 | Mach–Zehnder φ=90°: P(D) | 0.4993 | 0.5000 | max(3σ=0.0388, 0.02) | **PASS** | 3000 / 2983 (idem) | idem | 262 s |
| 45 | Mach–Zehnder φ=180°: P(D) | 0.0000 | 0.0000 | max(3σ=0.0020, 0.02) | **PASS** | 3000 / 2980 | 157.564 · 4e+07 | 241 s |
| 46 | Mach–Zehnder φ=270°: P(D) | 0.5043 | 0.5000 | max(3σ=0.0387, 0.02) | **PASS** | 3000 / 2980 (idem) | idem | 241 s |
| 47 | bomb tester: P(explosion) | 0.5088 | 0.5000 | max(3σ=0.0274, 0.02) | **PASS** | 3000 / 1471 | 151.909 · 3.9e+07 | 139 s |
| 48 | bomb tester: P(C) = bomb detected without interaction (dark port) | 0.2341 | 0.2500 | max(3σ=0.0237, 0.02) | **PASS** | 3000 / 1471 (idem) | idem | 139 s |
| 49 | bomb tester: P(D) | 0.2571 | 0.2500 | max(3σ=0.0237, 0.02) | **PASS** | 3000 / 1471 (idem) | idem | 139 s |
| 50 | same interferometer without bomb: P(C) (dark port) | 0.0000 | 0.0000 | max(3σ=0.0030, 0.02) | **PASS** | 3000 / 1471 (idem) | idem | 139 s |
| 51 | Wheeler, BS2 inserted after the photon passed BS1: P(D) | 0.9946 | 1.0000 | max(3σ=0.0040, 0.02) | **PASS** | 1500 / 1490 | 157.564 · 4e+07 | 218 s |
| 52 | Wheeler, BS2 removed after the photon passed BS1: P(D) | 0.5053 | 0.5000 | max(3σ=0.0548, 0.02) | **PASS** | 1500 / 1490 (idem) | idem | 218 s |
| 53 | Wheeler, BS2 removed after the photon passed BS1: fraction of photons whose detector is the straight continuation of their arm (C⇔A, D⇔B) | 0.0000 | 0.0000 | max(3σ=0.0040, 0.02) | **PASS** | 1500 / 1490 (idem) | idem | 218 s |

Fallimenti v1 e diagnosi:
* **Arago** (righe "Arago"): 1,4% dei clic nella macchia invece di 2,9%, contrasto 1,6 invece di 4,1.
  Diagnosi in tre passi: (1) l'onda delle palline coincide con quella esatta sullo schermo (entro 1%);
  (2) le traiettorie integrate coincidono con un RK4 fine indipendente (stessi punti d'arrivo a 10⁻³ λ);
  (3) le traiettorie che finiscono nella macchia partono, a z = δ, **dentro l'ombra geometrica** (r < 3 λ,
  fino a 0,4 λ dentro il disco). L'onda di Rayleigh–Sommerfeld di un bordo netto manda nell'ombra circa il
  2% del flusso già a z = δ, **per ogni δ** (0,0211 / 0,0202 / 0,0196 per δ = 0,02 / 0,005 / 0,002): la
  moneta può solo togliere fotoni, non metterli dove non sono passati. Ridurre δ non aiuta; serve una regola
  che sposti il fotone lungo il flusso vicino al bordo → **v2**.
  (Nella stessa indagine è emerso anche **B3**: le nascite vicino ai bordi erano campionate a celle intere di
  λ/16, troppo poche nelle celle tagliate dal bordo; corretto in v2. Con B3 corretto ma la regola v1 la
  macchia dà ancora 1,8% invece di 2,9%: FAIL. Con la regola v2: 3,1%. Quindi serve proprio la v2.)
* **Fascio libero 3D** (2 righe su 3 FAIL, +2,4%): fluttuazione del campione di nascita (larghezza iniziale
  +2,0% con 4 000 fotoni, 2,4σ); le traiettorie stesse scalano esattamente. In v2: 10 000 fotoni e la riga
  aggiuntiva "allargamento delle stesse traiettorie".
* **Fenditura lunga 3D, frazione assorbita** (0,740 invece di 0,781): errore di contabilità, la frazione di
  fascio sulla parete era calcolata solo dentro il riquadro delle fenditure (B4); corretto in v2.
* Le righe Mach–Zehnder fase 180/270, bomba e Wheeler della v1 sono state scritte dal programma quando
  l'etichetta era già passata a "v2"; il percorso di calcolo (nessun bordo nel fascio) è identico nelle due
  versioni, e la v2 le ha comunque rifatte.

### v2 — 2026-10-03: zona di bordo

Unico cambiamento della regola: **zona di bordo di λ/2** (§1.2 punto 6). Più le correzioni B3, B4 del
programma. **Tutti** gli esperimenti rifatti (tabella §4): 56 righe, tutte PASS. La doppia fenditura cambia di poco
(correlazione 0,9965 → 0,9962, rms 4,0% → 4,3%), la macchia di Arago passa da 1,4% a 2,9% (teoria 2,93%).

---

## 6. Cosa non funziona, o funziona solo in parte (onestà)

1. **La zona di bordo è una correzione del modello di bordo, non fisica "prima".** Serve perché il bordo di
   Kirchhoff/Rayleigh–Sommerfeld è singolare: il 2% del flusso passa nell'ombra in uno strato sottile come δ.
   La regola v2 lo tratta con un riordino monotono in una dimensione (attraverso il bordo). Limiti noti:
   * non sposta fotoni **lungo** il bordo. In 3D con luce polarizzata il flusso vicino al bordo si concentra
     dove E è perpendicolare al bordo (fino a 1,6 volte il flusso incidente) e si svuota dove è parallelo:
     dove il segmento ha guadagno (> 1) la moneta è tagliata a 1 e la parte esterna della zona resta un po'
     scarica. Effetto residuo visto: frazione nella macchia di Arago giusta (2,94% e 2,86% nelle due corse,
     esatto 2,93%), ma contrasto macchia/anello d'ombra 3,3 contro 4,1 esatto (passa la soglia 3 fissata
     prima, però è più basso; non ho isolato la causa);
   * la larghezza λ/2 è giustificata (larghezza di un gradino ridotto alla banda propagante) ma è un numero
     in più nella regola;
   * legge l'onda su un segmento lungo λ: è locale alla scala λ, non puntuale.
2. **La moneta toglie soltanto.** Dove Φ_out/Φ_in > 1 (lontano dai bordi, raro) viene tagliata a 1.
3. **Le direzioni non contano per la guida.** Con l'ondina esatta dall'origine e il cappello normalizzato il
   numero di direzioni è solo contabilità: lo riportiamo (Nyquist) ma la prova di convergenza su di esso
   sarebbe piatta per costruzione. Il numero che conta è quello delle origini (η, ε), ed è lì che la
   convergenza si vede.
4. **Con i fotoni, η = 1 e η = 2 si distinguono poco** con 2 500 clic; la differenza netta è nell'onda.
5. **Palline evanescenti**: il loro numero totale è enorme (10⁶–10¹² per fotone); è il conto dei gradi di
   libertà fino a δ, ma in pratica se ne leggono poche migliaia per passo. Il costo vero è il campo vicino in 3D.
6. **Geometria**: in 3D solo oggetti piani perpendicolari al fascio; il Mach–Zehnder con specchi e
   beam-splitter a 45° è fatto in 2D (motore 2D generale). Il beam-splitter è un foglio sottile ideale.
7. **Specchi**: 0–7 fotoni su 1 500 (≤ 0,5%) toccano uno specchio per errore numerico vicino alla superficie
   (onda stazionaria) e vengono contati come persi; 1 fotone per corsa si ferma al limite di passi.
8. **Tempi**: nella prima corsa completa v2 tre script hanno superato i 5 minuti (fascio libero 3D 363 s,
   Arago 421 s, bordo 360 s, macchina condivisa). Per rispettare il limite: fascio libero 3D e Arago con 7 000
   fotoni invece di 10 000 e analisi del bordo più rapida (riferimento esatto su una griglia meno fitta). La
   seconda corsa completa v2 (stesso codice, solo queste tre modifiche) è quella della tabella §4: tutti gli
   altri numeri sono identici bit per bit alla prima (seed fissi). La seconda corsa è partita per errore con
   l'uscita su schermo scartata (numeri e tempi sono nei file json); il registro della prima corsa v2 è
   `base_run_all_v2_batch1.log`, i suoi json in `lab/out/base_v2_batch1/` (Arago 0,0286 ± 0,0017 con 10 000
   fotoni, fascio libero 3D +1,1…1,3%). Anche alcune corse della prova di convergenza superano i 5 minuti se
   fatte insieme (η = 0,5: 212 s da sola, perché l'onda sbagliata rende le traiettorie lente).
9. **Entanglement** e più fotoni: non fanno parte di questo lavoro (vedi `docs/MULTIFOTONE.md`, altro lavoro).
10. **Il motore Rust in `src/photon/` non è stato toccato**: la regola base vive in `lab/`.

---

## 7. Come rieseguire

`bash lab/base_run_all.sh` esegue tutto in fila (prove della macchina, convergenza, 2D, 3D, Mach–Zehnder),
ogni script con `timeout 570`. Singolarmente:

| comando | cosa fa |
|---|---|
| `python lab/base_test.py` | Rust = numpy (2·10⁻⁷ 2D, 2·10⁻¹⁵ 3D), ∇·E = 0, palline vs Rayleigh–Sommerfeld, 3D → 2D, velocità |
| `python lab/base_conv.py wave \| photons 0.5 \| photons 1 2 4 \| eps` | convergenza (§2.4) |
| `python lab/base_exp2d.py free \| double \| which \| single \| grating \| talbot \| edge` | esperimenti 2D (double include la costruzione un clic alla volta) |
| `python lab/base_exp3d.py free \| gouy \| slit3d \| airy \| arago \| malus \| three` | esperimenti 3D |
| `python lab/base_mz.py phase 0 90 \| phase 180 270 \| bomb \| wheeler` | Mach–Zehnder, bomba, Wheeler |
| `python lab/base_table.py v2` | tabella dei risultati dai file json |

File: `lab/base_rule.py` (la regola, gli oggetti, la pallina-fotone, il riferimento numpy delle somme),
`lab/base_kernel/` (somme in Rust, compilate da `base_rule.py` con `cargo build --release` se manca la dll),
`lab/base_exp2d.py`, `lab/base_exp3d.py`, `lab/base_mz.py`, `lab/base_conv.py`, `lab/base_test.py`,
`lab/base_table.py`, `lab/base_run_all.sh`. Uscite in `lab/out/base_*`.
