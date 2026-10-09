# La carica come giro attorno a una dimensione nascosta

Laboratorio: `lab/kk_*.py` (Python 3.14 + numpy), costruito sulla regola zig-zag di `lab/zz_core.py`
(`docs/ZIGZAG.md`). Registri `lab/out/kk_*.log`, numeri `lab/out/kk_*.json`, figura `lab/out/kk_force.png`.
Prototipo veloce (circa un'ora): tutti i numeri qui sotto vengono da esecuzioni reali, deterministiche
(seed 20261009). Unità del reticolo: passo a = tick = 1, c = ħ = 1.

## 1. L'idea (di Valerio)

> "La carica non è una cosa dentro la particella: è una **circolazione attorno a una dimensione
> nascosta** di cui vediamo solo metà. Tre circolazioni a 120° si bilanciano e danno il neutro; se ne
> manca una, la circolazione che resta la vediamo come carica, e il suo verso distingue elettrone e
> protone."

Chi ci ha già pensato: **Kaluza (1921) e Klein (1926)**: con una quinta dimensione arrotolata su un
cerchio piccolissimo, la carica è la quantità di moto attorno al cerchio e il campo elettromagnetico è la
"torsione" del cerchio (un pezzo della geometria). **Wheeler (1955-57)**: "carica senza carica", linee di
campo intrappolate in un buco della geometria. **La tripla del colore (trialità Z₃)**: i quark hanno tre
"colori", e solo le combinazioni con colore totale nullo (tre quark, o quark + antiquark) esistono libere.

## 2. Che cosa abbiamo fatto

La regola zig-zag di sempre (ogni pallina va a c e gira il suo stato; la massa è l'angolo di flip θ),
con **una sola aggiunta**: oltre a x (e y) c'è un terzo asse **w**, un cerchio chiuso di N siti, raggio
R = N/2π. La pallina salta anche lungo w **con la stessa regola** (traslazione + rotazione). Niente forza
messa a mano: l'unica cosa in più, per il test della forza, è che il cerchio nascosto può essere
**storto**: quando la pallina fa un passo lungo y, scivola anche di un tratto D attorno al cerchio
(un termine fuori diagonale della metrica, g_yw: è il potenziale vettore A = D/R scritto come geometria).

## 3. Che cosa è emerso (numeri)

**Test 1 - la carica è intera da sola, e si conserva** (`kk_charge_mass.py`)
* Uno stato **a caso** sul cerchio (N = 32, θ = 0,1), lasciato andare per 4096 tick: lo spettro in
  frequenza ha **17 righe**, una per ogni |n| = 0…16, tutte sulle energie dei modi interi
  (distanza massima 7·10⁻⁴, risoluzione 1,5·10⁻³); la riga più vicina a un "n semintero" è a 4·10⁻²,
  cioè 26 volte la risoluzione. Il cerchio chiuso permette solo giri interi: k_w = n/R.
* Un pacchetto con n = +1, −1, +2 in un potenziale casuale e con il cerchio storto a caso (2+1 D,
  300 tick): la carica ⟨n⟩ resta **+1, −1, +2 esatti** (deriva 0 in doppia precisione).
* Controprova: se il cerchio è **bitorzoluto** (potenziale che dipende da w) la carica **non** si
  conserva più (deriva 1,35 / 0,85 / 2,25). La conservazione della carica = il cerchio è liscio
  (simmetria per rotazione, teorema di Noether).

**Test 2 - quanto costa girare (il problema di Kaluza-Klein)**
Energia a riposo misurata sulla regola, θ = m = 0,1 (λ_C = 10 celle):

| R / λ_C | E(n=1)/m | E(n=2)/m | errore rispetto a √(m² + (n/R)²) |
|---|---|---|---|
| 0,25 | 4,046 | 7,904 | −1,6·10⁻³ |
| 0,51 | 2,201 | 4,046 | −1,3·10⁻³ |
| 1,00 | 1,411 | 2,228 | −8·10⁻⁴ |
| 2,01 | 1,117 | 1,411 | −3·10⁻⁴ |
| 4,01 | 1,031 | 1,117 | −1·10⁻⁴ |

E = √(m²c⁴ + (nħc/R)²): girare costa. Perché il giro dell'elettrone non costi più della sua massa
serve **R ≳ ħ/(mₑc) = 3,86·10⁻¹³ m**. Ma Kaluza-Klein lega R alla forza della carica:
α = 4Għ/(c³R²) ⇒ per α = 1/137 serve **R = 2 l_P/√α = 3,8·10⁻³⁴ m**, e allora il primo modo carico
pesa **5·10¹⁷ GeV ≈ 10²¹ mₑ**. Al contrario, con R = 3,86·10⁻¹³ m la carica sarebbe α ≈ **7·10⁻⁴⁵**
(forte come la gravità). **La tensione è netta: in questa forma semplice l'elettrone non può essere
un "giro" attorno al cerchio di Kaluza-Klein.** (È il problema storico; in più Witten, 1981: KK non dà
facilmente particelle sinistrorse/destrorse diverse.)

**Test 3 - la forza viene dalla geometria** (`kk_cyclotron.py`, `kk_ab.py`, figura)
* Cerchio storto in modo uniforme (2+1 D, 128×128 × N = 24, p = 0,3, b = 0,015): **n = +1 gira in
  senso antiorario, n = −1 in senso orario, n = 0 va dritto** (deviazione laterale 0,65 celle su 96).
  Raggio del picco di densità 18,7 e 18,0 celle contro p/(|n|b) = 20; periodo 197 e 184 tick contro
  204. Accordo al 10 %, non di più: il pacchetto si allarga e si sfasa lungo l'orbita; non l'ho fatto
  convergere (prototipo).
* **Aharonov-Bohm** su un anello (400 siti × N = 32): i due mezzi pacchetti fanno i due lati dell'anello
  e si incontrano dall'altra parte. Sull'anello il cerchio nascosto **non ha curvatura** (nessun
  "campo"): la torsione totale Φ è tutta nel buco. Eppure la fase dell'interferenza si sposta di
  **n·Φ** (errore massimo 1,5·10⁻¹¹ su 24 casi, n = −1, 0, +1, +2); con Φ = 2π (una scivolata di un
  giro intero) non cambia niente. Torsione su un solo legame o spalmata su tutti: densità uguali a
  3·10⁻¹⁴ (conta solo il totale nel buco). Le frange misurate direttamente con un fit: stesso spostamento
  a 0,6 % (n = 1).

**Test 4 - il "tre" di Valerio** (`kk_z3.py`): cerchio di **N = 3 siti** (120°), giri n ∈ {0, +1, −1}.
Fase AB misurata con una torsione di **un sito** (una rotazione di 120°): n = ±1 → ∓2,0944 rad
(= 2π/3, esatto), n = 0 → 0. Neutri (fase 0 mod 2π):
* 1 pallina: solo (0);
* 2 palline: (0,0) e **(+1,−1)**;
* 3 palline: (0,0,0), **(0,+1,−1)**, **(+1,+1,+1)**, **(−1,−1,−1)** (4 su 10 combinazioni).
* Togliendone una da (0,+1,−1): resta (+1,−1) neutro, oppure (0,−1) con fase +2π/3 o (0,+1) con
  −2π/3: **la carica che resta ha verso opposto a seconda di quale manca**, come diceva Valerio.
* Con una torsione **non** di un sito intero (0,37 siti) (+1,+1,+1) **non** è più neutro (fase 3·0,775):
  la neutralità dei tre è una proprietà Z₃ (rotazioni a scatti di 120°), non della carica elettrica U(1).

## 4. Che cosa NON funziona (onestamente)

* **Massa** (test 2): con la carica osservata il giro costa ~10¹⁷ GeV. Non è un dettaglio.
* **Il "tre" è il colore, non la carica elettrica.** I neutri Z₃ sono quelli della trialità dei quark
  (barioni qqq, mesoni q q̄). Ma: (a) qui un cerchio di 3 siti conosce solo n mod 3 (n = +1 e n = −2 sono
  la stessa cosa), quindi **nessuna carica frazionaria** 2/3, −1/3; nel Modello Standard i terzi vengono
  dall'ipercarica (Q = T₃ + Y), e il legame con il colore è il quoziente Z₆ del gruppo, che qui **non**
  si ricava; (b) **l'atomo di idrogeno è neutro con due costituenti** (e⁻, p⁺): per la carica elettrica
  basta la coppia (+1, −1), il "tre" non serve.
* Il cerchio di 3 siti è lontanissimo dal continuo: energia del giro 1,059 invece di 2,104, velocità di
  gruppo negativa (artefatto del reticolo). Robusta è solo la fase a 120°.
* Artefatto trovato: con k_w·a = π/4 (N = 8, n = 1) la regola a passi divisi ha banda **piatta** nel
  piano (il pacchetto non si muove). Servono cerchi di molti siti (qui 24-32) per avere fisica pulita.
* I campi sono **di sfondo**: il cerchio storto lo metto io. Due palline cariche non si attraggono né si
  respingono ancora.

## 5. Il prossimo passo (lungo)

Rendere **dinamica** la geometria nascosta: la torsione D (cioè A) deve avere una sua regola (la luce
di `zz_core.py`, spin 1, F = E + icB) e deve essere **sorgente** delle palline cariche (la corrente
di n). Allora due pacchetti con n = +1 e n = −1 dovrebbero attirarsi da soli (e due +1 respingersi),
con la legge 1/r² in 3+1 D, e la costante α sarebbe fissata da R. È Kaluza-Klein completo sul reticolo
(con la gravità di sfondo congelata): difficile, e porta dritto al problema del test 2.

## 6. File

`lab/kk_core.py` (regola con il cerchio nascosto), `lab/kk_charge_mass.py` (test 1-2, 21 s),
`lab/kk_cyclotron.py` (test 3a, 123 s), `lab/kk_ab.py` (test 3b, 72 s), `lab/kk_z3.py` (test 4, 15 s),
`lab/kk_figure.py` → `lab/out/kk_force.png`.
