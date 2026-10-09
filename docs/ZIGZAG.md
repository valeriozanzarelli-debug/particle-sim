# Zig-zag: una sola regola per la luce, l'elettrone e l'atomo di idrogeno

Laboratorio: `lab/zz_*.py` (Python 3.14 + numpy/scipy/numba). Visore: `cargo run --bin zigzag`
(sorgenti in `src/zigzag/`). Unità naturali ħ = c = mₑ = 1: lunghezza λ_C = ħ/(mₑc) = 3,8616·10⁻¹³ m,
tempo ħ/(mₑc²) = 1,2881·10⁻²¹ s, energia mₑc² = 510 998,95 eV. Tutti i numeri qui sotto vengono da
esecuzioni reali (seed fissi, macchina a 4 core condivisa con altri due agenti). Registri in
`lab/out/zz_*.log`, risultati in `lab/out/zz_*.json`, figure `lab/out/zz_*.png`, immagini del visore
`lab/out/zigzag_*.png`.

---

## 0. In breve

**L'idea.** Una sola regola costruisce tutte le particelle con le stesse palline: ogni pallina-onda
**si sposta alla velocità della luce e gira il suo stato interno**. Cambia solo il "tipo" di pallina
(lo spin), un angolo di "flip" (la massa) e una fase (la carica). Il surfista legge le palline e si
muove con v = c ψ†αψ/ψ†ψ (la regola di `SURFISTA_DIRAC.md`).

**Che cosa abbiamo trovato (tutto in questo documento, con i numeri e i fallimenti):**

* **Luce** (palline di spin 1 che portano F = E + icB): lungo gli assi del reticolo la regola è **esatta**
  (ω = c|k| per ogni lunghezza d'onda); in diagonale l'errore va come (a/λ)² (−3,5 % a λ = 8a, −0,21 % a
  32a). La luce resta trasversale (energia nel modo "fermo" 4·10⁻³ a 8a, 3·10⁻⁴ a 16a) e ha spin 1
  (fase di Berry 2π·(1 − cos 60°) = π su un cono, come misurato da Tomita e Chiao nel 1986).
  **Ma** la diffrazione è più debole del vero di un fattore ka·cot(ka): con λ = 8a la doppia fenditura
  non riproduce le frange (correlazione 0,60, **FALLITO**); con λ = 16a sì (0,980 a bin fini, 0,996 a bin
  grossi), con 24a 0,996 / 0,999. Serve λ ≳ 16a.
* **Elettrone libero** (palline di Weyl sinistrorse e destrorse mescolate dal flip θ = mₑca/ħ): energia
  a riposo **esatta**, E² = p²c² + m²c⁴ con errore 0,19 % fino a p = mc per θ = 0,1 (scende come θ²);
  Zitterbewegung con ω = 2,0095 mc²/ħ e ampiezza 0,491 ħ/mc, **uguale** al calcolo del continuo
  (`dirac_free.py`); in 3-D ⟨r⟩ gira su un cerchio di raggio ~ħ/2mc; lo spinore ha periodo 4π (una
  rotazione del reticolo di 90° vuole lo spinore girato di 45°); g → 2 dalla corrente di Dirac (con un
  errore del reticolo ∝ θ² che abbiamo misurato e spiegato).
* **Idrogeno** (lo stesso elettrone nella fase di Coulomb di un protone puntiforme): i livelli sono le
  frequenze a cui "suona" la regola. Estrapolati a reticolo nullo (gli s con i reticoli fini del §5.4):
  1s₁/₂ = **13,60575 eV** (Dirac
  13,60587; Bohr 13,60569), tutti i livelli p e d a 10⁻⁸ dalla formula di Dirac, i livelli s a ~10⁻⁵ (un
  errore del reticolo vicino al nucleo, §5.4). Struttura fine 2p₃/₂ − 2p₁/₂ = **10,9437 GHz** (Dirac
  10,9437; misurato 10,9690: i 25 MHz mancanti sono QED). Lyman-α **121,568 nm** (misurato 121,567),
  Balmer-α **656,467 nm** (656,46). Orbitali 1s/2s/2p uguali a quelli veri a 10⁻⁶. Il surfista dell'1s
  gira attorno allo spin a **α c sinθ = 2,188·10⁶ m/s** (come in `SURFISTA_DIRAC.md`). Vita del 2p dalle
  funzioni d'onda della regola: 1,600 ns (A = 6,249·10⁸ s⁻¹ contro 6,265·10⁸). L'1s lo trovano anche le
  palline da sole, "ascoltando" (8·10⁻⁷ dal risolvente). La regola cubica 3-D completa lega l'elettrone
  (provata a Z = 60) ma sui reticoli possibili qui sbaglia ancora del 10 % (**FALLITO**, §5.6): l'idrogeno
  vero è calcolato nella forma radiale della stessa regola.
* **Ingressi** (messi a mano): c, ħ (unità), massa dell'elettrone (= angolo di flip), carica (α), Z = 1,
  massa del protone (solo nella massa ridotta), spin delle palline. **Non** spiega: perché quelle masse e
  quella carica, lo spostamento di Lamb e g − 2 (QED), il protone.

Esito dei test (tolleranze fissate nel codice prima di ogni prima esecuzione): vedi la tabella del §8.
Molti fallimenti sono **test scritti male da me** (tolleranze troppo strette, un pacchetto confrontato
con il valore sbagliato, due errori di codice): sono tutti elencati, con la causa, nel §9.

---

## 1. Un po' di storia: chi ha già pensato allo "zig-zag"

| anno | chi | che cosa |
|---|---|---|
| 1928 | P. A. M. Dirac, *Proc. R. Soc. A* 117, 610 | l'equazione dell'elettrone: la "velocità" c α ha autovalori ±c. Ogni misura istantanea della velocità darebbe ±c |
| 1929 | H. Weyl, *Z. Phys.* 56, 330 | senza massa l'equazione si spezza in due pezzi indipendenti: sinistrorso (L) e destrorso (R), ciascuno va sempre a c |
| 1930 | E. Schrödinger, *Sitzungsber. Preuss. Akad. Wiss.* 24, 418 | Zitterbewegung: la posizione media trema a 2mc²/ħ con ampiezza ħ/2mc |
| 1946–1965 | R. P. Feynman (Feynman e Hibbs, *Quantum Mechanics and Path Integrals*, 1965, probl. 2-6) | la "scacchiera": in 1+1 D l'elettrone va sempre a c, a destra o a sinistra; a ogni svolta un fattore i·m·ε. Sommando su tutti i cammini a zig-zag si ottiene Dirac |
| 1957 | C. S. Wu et al., *Phys. Rev.* 105, 1413 | la natura distingue sinistra e destra (violazione della parità): L e R sono davvero due oggetti diversi |
| 1964 / 2012 | P. Higgs, *PRL* 13, 508; F. Englert, R. Brout, *PRL* 13, 321; ATLAS e CMS, *Phys. Lett. B* 716, 1 e 30 | la massa dell'elettrone nasce dall'accoppiamento che trasforma L in R (e viceversa) con il campo di Higgs |
| 2004 | R. Penrose, *The Road to Reality*, cap. 25 | l'elettrone come "zig" (L) e "zag" (R) che si scambiano di continuo; la massa è il ritmo dello scambio |
| 1993 | Y. Aharonov, L. Davidovich, N. Zagury, *PRA* 48, 1687 | le "passeggiate quantistiche" (quantum walks) |
| 1994 | I. Bialynicki-Birula, *PRD* 49, 6920 | Weyl, Dirac **e Maxwell** sul reticolo come automi cellulari unitari; per la luce usa F = E + icB, che obbedisce a un'equazione "alla Weyl" con spin 1 |
| 1996 | D. A. Meyer, *J. Stat. Phys.* 85, 551 | automi cellulari quantistici: niente di non banale senza uno stato interno ("da automi cellulari a gas reticolari quantistici") |
| 2006 | F. W. Strauch, *PRA* 73, 054302 | limite relativistico delle passeggiate quantistiche |
| 2012–2014 | G. Di Molfetta, F. Debbasch (*J. Math. Phys.* 53, 123302); P. Arrighi, V. Nesme, M. Forets (*J. Phys. A* 47, 465302) | passeggiate quantistiche a "passi divisi" che convergono a Dirac in 3-D; spostamento lungo x, y, z uno alla volta |
| 2010 | R. Gerritsma et al., *Nature* 463, 68 | Zitterbewegung vista in un simulatore (uno ione intrappolato); sull'elettrone libero non è mai stata osservata |

Quindi la regola non è inventata qui: è la passeggiata quantistica di Dirac a passi divisi, con la luce
scritta nella stessa forma (Bialynicki-Birula). Il nostro lavoro: scriverla come "che cosa fa ogni
pallina a ogni tick", controllarla contro le misure e vedere fin dove arriva senza aggiungere altro.

---

## 2. La regola (congelata il 2026-10-09 in `lab/zz_core.py`, prima di ogni test)

Reticolo cubico di passo a; un tick dura a/c.

**1. SPOSTA.** Lungo x, poi y, poi z (nei tick dispari l'ordine è z, y, x: due tick formano un
palindromo e la regola diventa precisa al secondo ordine): lo stato interno della pallina si divide sugli
autovettori della sua matrice di spin S_j. La parte con autovalore +1 si sposta di **una cella avanti**,
quella con −1 di **una cella indietro**, quella con 0 (solo la luce) **resta**. Una cella per tick = c:
nessuna pallina va più veloce o più piano.

| tipo di pallina | stato interno | S_j | chi è |
|---|---|---|---|
| Weyl destrorsa R | 2 numeri complessi | σ_j | metà dell'elettrone |
| Weyl sinistrorsa L | 2 numeri complessi | −σ_j | l'altra metà |
| luce (spin 1) | F = E + icB, 3 numeri complessi | (S_k)_ij = −i ε_kij | il fotone di Bialynicki-Birula |

**2. RUOTA.** In ogni cella, dopo lo spostamento:

    L' = e^{−iVΔt} (cos θ · L − i sin θ · R)
    R' = e^{−iVΔt} (cos θ · R − i sin θ · L)

* θ = **m c² Δt/ħ = m c a/ħ**: è la **massa** (il "flip" di Higgs). θ = 0 → senza massa (la luce, i neutrini
  di Weyl). Per l'elettrone θ = a/λ_C: con a = 0,1 λ_C l'angolo è 0,1 rad per tick. **Ingresso: mₑ.**
* V = −Zαħc/r: la **carica** del nucleo (puntiforme, Z = 1) come fase per tick. **Ingressi: α, Z.**

**3. SURFISTA.** v = c (R†σR − L†σL)/(|L|² + |R|²): è la regola di Bohm–Dirac, scritta con L e R. Vale
|v| ≤ c per **qualunque** spinore (disuguaglianza di Cauchy–Schwarz). In una dimensione: v = c (|R|² − |L|²)/(|R|²
+ |L|²): il surfista va a destra quanto la pallina è destrorsa.

**Limite continuo** (a → 0): i ∂ψ/∂t = (c α·p + β mc² + V)ψ (Dirac, rappresentazione chirale) e
i ∂F/∂t = c ∇×F (Maxwell nel vuoto). Il passo a **non è fisica**: ogni numero è dato per più valori di a ed
estrapolato a a → 0.

**Forma radiale (idrogeno).** Per uno stato di momento angolare dato (numero quantico di Dirac κ) la stessa
regola si scrive lungo il raggio: un camminatore χ₊ che va **fuori** a c e uno χ₋ che va **dentro** a c;
la rotazione è exp(−iΔt (m σ_z + κ/r σ_x − Zα/r)) su (g, f) = componenti radiali grande e piccola; in r = 0
chi va dentro torna fuori (g(0) = 0). È esattamente la regola 3-D per quello stato (nel continuo).

**Perché la forma radiale.** La regola 3-D deve risolvere la lunghezza ħ/mc (θ ≪ 1) **e** contenere qualche
raggio di Bohr (a₀ = 137 ħ/mc): servirebbero ~(5 a₀ / 0,1 λ_C)³ ≈ 3·10¹¹ celle. La regola 3-D sul problema
di Coulomb è comunque provata, con Zα più grande (§5.6).

**Come si trovano i livelli: "ascoltare" le palline.** Uno stato stazionario è una frequenza a cui la
regola suona. Se si fa girare la regola e si sommano le ampiezze delle palline con una fase e^{iωt} (un
lock-in, come in un ricevitore radio), sopravvive solo lo stato con frequenza ω; più a lungo si ascolta,
più è pulito. Ascoltare per un tempo infinito è il risolvente (1 − e^{iωΔt}W)⁻¹ = Σₙ (e^{iωΔt}W)ⁿ, che si
calcola in un colpo con l'eliminazione di Gauss ("iterazione inversa"). `zz_listen.py` controlla che
l'ascolto vero converga a quel valore. **Perché non il tempo immaginario (o la DMC di `MOLECOLE.md`)?**
Lo spettro di Dirac non ha un fondo (il "mare" a energia negativa): un cammino in tempo immaginario
affonderebbe nel limite del reticolo, −πħ/Δt, invece di fermarsi all'1s (il collasso variazionale di
Brown e Ravenhall, *Proc. R. Soc. A* 208, 552 (1951)).

---

## 3. T1 — La luce (`zz_light.py`, `zz_slit.py`, `zz_slit_conv.py`)

### 3.1 Dispersione del reticolo (analisi di Bloch esatta dei due tick)

Errore relativo ω/(c|k|) − 1:

| direzione | λ = 4a | 6a | 8a | 12a | 16a | 32a |
|---|---|---|---|---|---|---|
| asse (1,0,0) | 0 | 0 | 0 | 0 | 0 | 0 |
| faccia (1,1,0) | −12,6 % | −5,02 % | −2,71 % | −1,17 % | −0,65 % | −0,16 % |
| corpo (1,1,1) | −15,0 % | −6,40 % | −3,53 % | −1,54 % | −0,86 % | −0,21 % |

Lungo gli assi lo spostamento è esatto: la luce va a c per ogni λ. In diagonale la regola è un po' più
lenta, con errore ∝ (a/λ)² (rapporto 8a/16a = 4,09). Nel piano x–z vale una formula chiusa:
**cos(ωΔt) = cos(k_x a) cos(k_z a)** (la stessa forma di cos(ωΔt) = cos θ cos(ka) dell'elettrone).

### 3.2 Velocità di un pacchetto

Pacchetto 3-D (larghezza 5a, elicità +1, 30 tick). La velocità del centro **misurata** coincide con quella
che prevede la dispersione della regola entro lo 0,4 %; rispetto a Maxwell nel continuo (stesso pacchetto)
è uguale lungo l'asse e più lenta in diagonale:

| λ, direzione | misurata | dalla dispersione della regola | Maxwell, stesso pacchetto |
|---|---|---|---|
| 8a, asse | 0,9794 c | 0,9802 c | 0,9841 c |
| 8a, faccia | 0,8941 c | 0,8975 c | 0,9841 c |
| 8a, corpo | 0,8705 c | 0,8743 c | 0,9841 c |
| 16a, asse | 0,9360 c | 0,9365 c | 0,9400 c |
| 16a, faccia | 0,9168 c | 0,9183 c | 0,9400 c |
| 16a, corpo | 0,9104 c | 0,9124 c | 0,9400 c |

(Un pacchetto largo 5a contiene molte direzioni, quindi il suo centro va più piano di c anche nel
continuo. La velocità di gruppo di un'onda piana, dalla dispersione: 1 sull'asse; 0,916/0,980 in faccia e
0,892/0,974 in diagonale a 8a/16a.)

### 3.3 Trasversalità e spin 1

* Una luce fisica (div F = 0) che entra nel reticolo mette nel modo "fermo" (frequenza 0, la parte 0 di S_j
  che non si muove) una frazione **4,26·10⁻³** dell'energia a λ = 8a, **3,46·10⁻⁴** a 16a (∝ (a/λ)⁴):
  quasi tutta la luce resta trasversale e si muove.
* **Fase di Berry.** Si porta la direzione di propagazione lungo un cono e si guarda di quanto ruota la fase
  dello stato che si muove: per spin s la teoria dà s·(angolo solido). Cerchio massimo (angolo solido 2π):
  2π → **+1** (esatto). Cono di 60° (angolo solido π): atteso π, trovato π − 0,335 / 0,082 / 0,020 / 0,0051
  rad a λ = 8a / 16a / 32a / 64a (∝ (a/λ)²). È l'effetto misurato sulla luce in una fibra avvolta a elica
  (Tomita, Chiao, *PRL* 57, 937 (1986)).

### 3.4 Doppia fenditura (la stessa di `REGOLA_BASE.md`)

λ = 1, fenditure larghe 2 con centri a ±4, fascio gaussiano largo 7, schermo a L = 160; muro di una cella
che assorbe le palline; una sorgente piana a z = −3 aggiunge ogni tick un'onda che va avanti (cos ωt);
bordi "a spugna". Si misura il flusso di Poynting medio sullo schermo e si confronta con il flusso esatto di
Rayleigh–Sommerfeld (`lab/optics.py`).

| λ/a | schermo | correlazione (bin grossi) | correlazione (bin fini) | errore frange rms |
|---|---|---|---|---|
| 6 | 160 | 0,377 | 0,332 | 62 % |
| 8 | 160 | **0,595** | 0,576 | **45 %** |
| 8 | 40 | 0,944 | 0,709 | — |
| 16 | 40 | **0,9955** | **0,980** | — |
| 24 | 40 | **0,9988** | **0,9961** | — |

**FALLITO** a λ = 8a (soglia 0,98, la regola base del fotone fa 0,9962). La causa, trovata dopo: in
diffrazione conta la **curvatura** della dispersione. Da cos ω = cos k_x a · cos k_z a, un fascio lungo z si
allarga con d²ω/dk_x² = cot(k_z a)/a invece di 1/k_z: **la diffrazione è troppo debole del fattore
ka·cot(ka)** (0,785 a λ = 8a, 0,948 a 16a). Un fascio largo 1,5 λ a λ = 8a si allarga a 11,6 invece di
13,9 dopo 120 tick (controllo spettrale esatto). È lo stesso errore che nell'elettrone dà la "massa
cinetica" (1 + θ²/3) m. Per la luce bisogna dunque lavorare con λ ≥ 16a (meglio 24a).

---

## 4. T2 — L'elettrone libero (`zz_free.py`, `zz_gfactor.py`)

### 4.1 E² = p²c² + m²c⁴

Errore relativo dell'energia della regola (analisi di Bloch esatta dei due tick), direzione (1,1,1):

| θ = a/λ_C | p = 0 | 0,1 mc | 0,5 mc | 1 mc | 2 mc |
|---|---|---|---|---|---|
| 0,4 | 1·10⁻¹⁶ | 2,7·10⁻⁴ | 7,1·10⁻³ | 3,1·10⁻² | 1,4·10⁻¹ |
| 0,2 | 2·10⁻¹⁶ | 6,7·10⁻⁵ | 1,8·10⁻³ | 7,8·10⁻³ | 3,4·10⁻² |
| 0,1 | 4·10⁻¹⁶ | 1,7·10⁻⁵ | 4,4·10⁻⁴ | 1,9·10⁻³ | 8,5·10⁻³ |

* L'energia a riposo è **esatta** (mc², a 10⁻¹⁶): il flip θ **è** la massa.
* Lungo un asse vale cos(ωΔt) = cos θ cos(ka) (Feynman–Bialynicki-Birula). Per p piccolo l'energia
  cinetica è p²/(2 m_kin) con m_kin = **(1 + θ²/3) m**: misurato 1,0570 / 1,01355 / 1,003347 per
  θ = 0,4 / 0,2 / 0,1 (formula 1,0533 / 1,01333 / 1,003333). È l'errore principale del reticolo e va a zero
  come θ².
* La velocità di gruppo più grande su tutta la zona di Brillouin (20 000 punti, θ = 0,2) è 0,980 c: nessuna
  onda della regola va più veloce della luce.

### 4.2 Zitterbewegung

**1-D** (stesso pacchetto di `dirac_free.py`: σ = 5 λ_C, spinore (1,1)/√2 nella rappresentazione di Dirac,
che qui è un pacchetto **tutto destrorso**):

| θ | ω [mc²/ħ] | ampiezza [ħ/mc] |
|---|---|---|
| 0,2 | 2,00944 | 0,49405 |
| 0,1 | 2,00951 | 0,49139 |
| 0,05 | 2,00951 | 0,49072 |
| continuo (`dirac_free.py`) | 2,00952 | 0,49050 |

ω = 2,0095 mc²/ħ = 1,560·10²¹ rad/s (2mc²/ħ = 1,553·10²¹), ampiezza 0,491 ħ/mc = 1,90·10⁻¹³ m (ħ/2mc = 1,93·10⁻¹³ m).
Lo stesso pacchetto proiettato sugli stati a energia positiva **della regola** non trema: residuo 1·10⁻¹⁵.
Nella regola si vede *perché* trema: le palline vanno sempre a c, il flip le gira da R a L (cambiano verso)
a un ritmo fissato dalla massa.

**3-D** (σ = 3, θ = 0,25, spinore (L, R) = ((0, i), (1, 0))/√2): ⟨r⟩ gira su un **cerchio** nel piano x–z:
regola raggio 0,433/0,438 ħ/mc, ω = 2,071; Dirac nel continuo (stesso pacchetto) 0,428/0,431, ω = 2,073.
(Per un pacchetto largo il raggio tende a ħ/2mc = 0,5; a σ = 3 è più piccolo, come in `SURFISTA_DIRAC.md`.)

### 4.3 Spin ½ (periodo 4π)

* **Fase di Berry** dello stato R che si muove, portando la direzione lungo un cerchio massimo:
  π → **e^{iπ} = −1**: un giro completo dà un segno meno (spin ½). Su un cono di 60°: atteso π/2, errore
  −0,168 / −0,041 / −0,0102 rad a λ = 16a / 32a / 64a (∝ (a/λ)²).
* **Rotazione del reticolo.** Ruoto il reticolo di 90° attorno a z e lo spinore di un angolo χ
  (exp(−iχΣ_z/2)), faccio 40 tick e confronto con "prima 40 tick, poi ruoto". Il disaccordo è minimo
  per **χ = 90°** (0,068, contro 0,67 a χ = 180° e 0,74 a χ = 0): la regola vuole lo spinore girato con la
  legge del "mezzo angolo". Quattro quarti di giro danno exp(−iπΣ_z) = **−1**, otto danno +1: periodo 4π.

### 4.4 Momento magnetico: g dalla corrente

Pacchetto a energia positiva (proiettato sugli stati della regola), spin su, fermo; μ = ½∫ r × j con
j = c ψ†αψ; g = 2μ/S.

| | σ = 4 | σ = 8 |
|---|---|---|
| regola, θ = 0,5 | 1,797 | 1,820 |
| Dirac nel continuo | 1,966 | 1,989 |

**FALLITO** (soglie 2·10⁻³ e 1 %). Ipotesi scritta prima della prova successiva: è la massa cinetica del
reticolo (la corrente di spin va come 1/m_kin). Prova (`zz_gfactor.py`, σ = 3, continuo g = 1,9442):

| θ | g regola | g/g_continuo − 1 | −(θ²/3)/(1 + θ²/3) |
|---|---|---|---|
| 0,5 | 1,7739 | −0,0876 | −0,0769 |
| 0,25 | 1,9023 | −0,0216 | −0,0204 |
| 0,125 | 1,9338 | −0,0054 | −0,0052 |

Confermato: g della regola → g di Dirac come θ², e g di Dirac → 2 al crescere del pacchetto (1,966 → 1,989
da σ = 4 a 8). Il g − 2 della QED (0,00116) non c'è, come in Dirac.

### 4.5 Surfisti

* |v| ≤ c a ogni stadio RK4 di 2000 surfisti (max 1,000000000000000).
* Equivarianza in 1-D (ogni surfista dovrebbe tenere il suo "posto in fila" nella densità): errore massimo
  1,26·10⁻² su 40 ħ/mc² a θ = 0,1 — **FALLITO** di poco (soglia 10⁻²). I surfisti leggono le palline con
  un'interpolazione lineare e un passo di due tick: è un errore del lettore, non della regola (con lo
  spettrale di `dirac_free.py` era 1,5·10⁻⁵).

---

## 5. T3 — L'idrogeno (`zz_hydrogen.py`, `zz_hydrogen_s.py`, `zz_listen.py`, `zz_atom3d.py`)

### 5.1 Ascoltare l'atomo (la dinamica delle palline)

Regola radiale per gli stati s (κ = −1) su un reticolo grosso, θ = a/λ_C = **1,0** (così la fase di legame per
tick è grande abbastanza), scatola di 20 raggi di Bohr con una "spugna" sul 20 % esterno (le palline che
si allontanano svaniscono). L'elettrone parte come una macchia gaussiana g = r e^{−r²/(2(2a₀)²)}, f = 0:
un miscuglio di 1s, 2s, 3s, … e di stati liberi. Niente formula, niente stato di prova "giusto".

1. **Spettroscopia.** 10⁶ tick; si registra quanto della macchia iniziale c'è ancora, C(n) = ⟨χ₀|χₙ⟩.
   Lo spettro di C(n) ha due righe (picchi a 21,117 e 5,578 eV di legame): la regola "suona" a quelle
   frequenze. Il risolvente sullo stesso reticolo dà 1s = 21,135156 eV e 2s = 5,290770 eV (su questo
   reticolo grossolano le energie sono più grandi del vero: massa cinetica (1 + θ²/3) e il resto, §5.2): picchi a
   0,018 e 0,29 eV, dentro la risoluzione 2πħ/(N Δt) = 3,2 eV (test S1 **PASS**).
2. **Lock-in.** Altri 10⁶ tick, sommando le ampiezze delle palline con la fase e^{iωt} dei due picchi sotto
   finestre gaussiane:

| finestra (deviazione standard) | 1s dall'ascolto | differenza dal risolvente |
|---|---|---|
| 41 667 tick | 18,837740 eV | 1,1·10⁻¹ |
| 83 333 tick | 21,118026 eV | 8,1·10⁻⁴ |
| 166 667 tick | **21,135139 eV** | **8,1·10⁻⁷** |

   Più a lungo si ascolta, più ci si avvicina al risolvente: lo stato fondamentale lo trovano le palline
   da sole (test S2 **PASS**). Il 2s dall'ascolto è 5,132 eV (−3,0 %, test S3 **FALLITO**): con 10⁶ tick il 2s
   e il 3s (a ~2,9 eV su questo reticolo) non sono separati, servirebbe ascoltare più a lungo.
   (Le prime due esecuzioni non sono finite in 590 s; la seconda ha anche mostrato righe false equidistanti
   di 5,06 eV: la parte libera della macchia rimbalzava sul bordo della scatola. Da qui la spugna.
   Registri `zz_listen_run1_killed.log`, `zz_listen_run2_killed.log`.)

### 5.2 I livelli (risolvente = ascolto infinito), quattro reticoli

Energia di legame (eV, nucleo infinitamente pesante):

| livello | θ = 0,4 | 0,2 | 0,1 | 0,05 | a → 0 | Dirac (formula) | errore relativo |
|---|---|---|---|---|---|---|---|
| 1s₁/₂ | 14,365109 | 13,782649 | 13,647771 | 13,615535 | 13,604922 | 13,605874 | 7,0·10⁻⁵ |
| 2s₁/₂ | 3,593302 | 3,446623 | 3,412409 | 3,404106 | 3,401361 | 3,401480 | 3,5·10⁻⁵ |
| 3s₁/₂ | 1,597315 | 1,531968 | 1,516689 | 1,512962 | 1,511729 | 1,511764 | 2,3·10⁻⁵ |
| 2p₁/₂ | 3,464117 | 3,416729 | 3,405268 | 3,402425 | 3,4014799 | 3,4014799 | 6·10⁻⁹ |
| 2p₃/₂ | 3,464072 | 3,416684 | 3,405222 | 3,402380 | 3,4014346 | 3,4014346 | 6·10⁻⁹ |
| 3p₁/₂ | 1,558679 | 1,523090 | 1,514571 | 1,512464 | 1,5117638 | 1,5117638 | 2·10⁻⁸ |
| 3p₃/₂ | 1,558665 | 1,523076 | 1,514558 | 1,512451 | 1,5117504 | 1,5117504 | 2·10⁻⁸ |
| 3d₃/₂ | 1,528376 | 1,515812 | 1,512760 | 1,512002 | 1,5117504 | 1,5117504 | 3·10⁻⁹ |
| 3d₅/₂ | 1,528372 | 1,515808 | 1,512756 | 1,511998 | 1,5117459 | 1,5117459 | 3·10⁻⁹ |

La formula di Dirac (Sommerfeld) è usata **solo per controllare**. Il test H1 (tutti entro 10⁻⁵) è
**FALLITO** per i livelli s; i p e i d sono esatti a 10⁻⁸.

### 5.3 Le grandezze misurabili

| grandezza | regola | misura | differenza | da dove viene |
|---|---|---|---|---|
| energia di ionizzazione 1s (massa ridotta) | 13,59834 eV (s raffinati, §5.4) | 13,598435 eV (NIST) | −7·10⁻⁶ | errore del reticolo sugli s (1·10⁻⁵) + Lamb 1s (+2·10⁻⁶, QED) |
| struttura fine 2p₃/₂ − 2p₁/₂ | **10,94367 GHz** (massa inf. 10,94963) | 10,96905 GHz | **−25,4 MHz (−0,23 %)** | QED (soprattutto il momento magnetico anomalo); Dirac dà 10,94966 GHz |
| Lyman-α (2p → 1s, baricentro 2:1) | 121,5780 nm (corsa 1); **121,5681 nm** con gli s raffinati | 121,567 nm | +90 ppm; **+9 ppm** | errore del reticolo sull'1s |
| Balmer-α (n = 3 → 2, baricentro) | **656,4674 nm** | 656,46 nm | +11 ppm | la larghezza fine della riga (le 7 componenti vanno da 656,445 a 656,498 nm) e il peso delle componenti |
| senza il rinculo del protone | Lyman-α 121,512 nm | | −455 ppm | la massa ridotta serve davvero |
| 1S–2S | 2,466042·10¹⁵ Hz | 2,466061·10¹⁵ Hz | −19 GHz | errore del reticolo sugli s (~30 GHz) più grande del Lamb (+7 GHz): qui la regola non basta |
| A(2p → 1s) dalle funzioni della regola | 6,249·10⁸ s⁻¹ (vita 1,600 ns) | 6,2649·10⁸ s⁻¹ (NIST) | −0,25 % | |

Il baricentro di Balmer-α pesa le componenti 3p→2s, 3s→2p, 3d→2p con popolazioni statistiche e con i
coefficienti A calcolati dalle funzioni d'onda della regola (A(3p→2s) = 2,24·10⁷, A(3s→2p) = 6,34·10⁶,
A(3d→2p) = 6,46·10⁷ s⁻¹, vedi `zz_hydrogen.json`) e le forze relative dentro i multipletti (algebra dei
momenti angolari, non un adattamento).

### 5.4 Perché gli s convergono peggio (`zz_hydrogen_s.py`)

Solo gli stati s hanno densità diversa da zero sul nucleo. Lì l'angolo di rotazione per tick κ/r·Δt è
grande (2 rad sulla prima cella) e l'errore del reticolo contiene un termine "di contatto" ∝ |ψ(0)|²:
infatti l'errore estrapolato scala come 1/n³ (1s : 2s : 3s = 9,5·10⁻⁴ : 1,2·10⁻⁴ : 3,5·10⁻⁵ eV). Con reticoli più
fini:

| θ | 0,2 | 0,1 | 0,05 | 0,025 | 0,0125 |
|---|---|---|---|---|---|
| 1s, errore relativo | +1,30·10⁻² | +3,08·10⁻³ | +7,10·10⁻⁴ | +1,54·10⁻⁴ | +3,07·10⁻⁵ |
| 2s | +1,33·10⁻² | +3,21·10⁻³ | +7,72·10⁻⁴ | +1,81·10⁻⁴ | +4,14·10⁻⁵ |

Estrapolazione dai tre più fini: 1s −9,2·10⁻⁶ (polinomio in θ²), −3,0·10⁻⁶ aggiungendo un termine
θ² ln θ; 2s −4,6·10⁻⁶ / −1,5·10⁻⁶. L'altra condizione al centro (f(0) = 0 invece di g(0) = 0) cambia l'1s a
θ = 0,05 di +1,3·10⁻⁴ (13,61726 invece di 13,61554 eV): anche questo sparisce con a → 0.

### 5.5 Orbitali e surfista

* Orbitali: 1 − |⟨regola|vero⟩|² = 4·10⁻⁷ (1s), 1,4·10⁻⁶ (2s), 1,6·10⁻⁷ (2p₁/₂), 1,5·10⁻⁷ (2p₃/₂), contro le
  funzioni R_nl dell'idrogeno; i massimi di r²|ψ|² a 0,999 a₀ (1s), 4,0 a₀ (2p), 5,23 a₀ (2s). Peso della
  componente piccola: 1,33·10⁻⁵ (1s), 3,3·10⁻⁶ (n = 2) ≈ α²/4 per l'1s.
* Surfista dell'1s: |v| all'equatore = **1,000649 α c** a θ = 0,05 (θ = 0,4 / 0,2 / 0,1: 1,0553 / 1,0127 /
  1,0029), estrapolato 0,99990 α c (H5b **FALLITO**: soglia 10⁻⁴; anche questo è l'errore degli s vicino al
  nucleo). f/g = −0,0036526 (Dirac −0,0036487). Il surfista gira attorno all'asse dello spin a
  2,188·10⁶ m/s · sinθ: è lo stesso risultato di `SURFISTA_DIRAC.md`, ora ottenuto dalle palline che fanno
  zig-zag.

### 5.6 La regola 3-D vera sul problema di Coulomb (`zz_atom3d.py`)

La regola cubica completa (palline a 4 componenti, spostamenti lungo x, y, z, flip, fase di Coulomb) non
può contenere l'idrogeno vero (§2), ma il problema di Dirac–Coulomb dipende solo da Zα. Allora: ione
idrogenoide di **neodimio, Z = 60** (Zα = 0,4378, raggio di Bohr λ_C/Zα = 2,28 λ_C), scatola periodica di
24 λ_C, nucleo fra 8 celle, reticoli a = 0,4 / 0,3 / 0,24 λ_C (60³, 80³, 100³ celle). Partenza: una macchia
esponenziale grande quanto il raggio di Bohr (dagli ingressi), proiettata sugli stati a energia positiva
della regola libera; poi lock-in alla sua frequenza media con una finestra di 27 ħ/mc². Ingressi: mₑ, α,
Z = 60. Dirac 1s: legame 0,100948 mc² = 51,584 keV.

| a [λ_C] | legame dalla regola 3-D | regola radiale, stesso a | 3-D / radiale | v_φ/(Zα c sinθ) | v_r rms /(Zα c) |
|---|---|---|---|---|---|
| 0,4 | 0,130043 mc² (+28,8 %) | 0,125310 (+24,1 %) | +3,8 % | 0,603 | 0,58 |
| 0,3 | 0,120024 (+18,9 %) | 0,114134 (+13,1 %) | +5,2 % | 0,845 | 0,50 |
| 0,24 | 0,111273 (+10,2 %) | 0,109349 (+8,3 %) | +1,8 % | 0,942 | 0,37 |
| 0,1 | — | 0,102584 (+1,62 %) | | | |
| 0,05 | — | 0,101430 (+0,48 %) | | | |
| 0,025 | — | 0,101093 (+0,14 %) | | | |

* La regola 3-D **lega** l'elettrone al nucleo con la sola fase di Coulomb, e su ogni reticolo dà quasi lo
  stesso numero della regola radiale (entro 2–5 %): la forma radiale usata per l'idrogeno è davvero la
  stessa regola.
* Ma i reticoli 3-D possibili qui sono troppo grossi: l'estrapolazione in a² sbaglia di −11,6 % (−5,2 % con
  i due più fini), test A1 **FALLITO** (soglia 1 %). Il motivo (`zz_atom3d_radial.py`, scritto dopo): con
  Zα = 0,44 l'onda di Dirac va come r^{γ−1} vicino al nucleo (γ = 0,899) e l'errore del reticolo va come
  a^{2γ} = a^{1,8}, non a²; e diventa piccolo (≤ 0,5 %) solo per a ≤ 0,05 λ_C, cioè 480³ celle.
* Il surfista dello stato 3-D gira attorno allo spin con v_φ che sale verso Zα c sinθ (0,60 → 0,84 → 0,94),
  ma lo stato filtrato ha ancora una velocità radiale grande (0,37 Zα c): non è ancora stazionario
  (ascolto troppo breve, miscuglio con il 2s e con gli stati liberi). Test A2 **FALLITO**.

---

## 6. T4 — Tabella di onestà

| | che cosa | come entra |
|---|---|---|
| **ingressi** | c, ħ | le unità: una cella per tick = c |
| | spin delle palline (½ per l'elettrone, 1 per la luce) | il tipo di pallina (le matrici S_j) |
| | massa dell'elettrone mₑ | l'angolo di flip θ = mₑca/ħ |
| | carica: α = 1/137,035999 | la fase per tick −Zα ħc/r |
| | Z = 1, protone puntiforme | la stessa fase |
| | massa del protone | solo nella massa ridotta μ (righe e struttura fine) |
| **escono** | velocità della luce per ogni λ lungo gli assi; trasversalità; spin 1 della luce | |
| | E² = p²c² + m²c⁴, Zitterbewegung (2mc²/ħ, ħ/2mc), spin ½ (4π), g → 2 | |
| | livelli dell'idrogeno con la struttura fine di Dirac, orbitali, righe, vite, il surfista che gira a αc sinθ | |
| **non spiega** | perché mₑ e α hanno quei valori | sono ingressi |
| | spostamento di Lamb (1057,8 MHz), g − 2 (0,00116), i 25 MHz della struttura fine | QED: il vuoto quantistico |
| | il protone (quark, gluoni, dimensione 0,84 fm, struttura iperfine) | è solo una carica puntiforme |
| | il rinculo oltre la massa ridotta; più elettroni | non fatti |
| | la creazione di coppie (l'energia negativa è solo "mare") | una regola a una particella non può |
| | perché la regola è proprio questa | è la passeggiata di Dirac nota; non l'abbiamo derivata da qualcosa di più semplice |

---

## 7. Il visore `zigzag`

    cargo run --bin zigzag                         (oppure target\debug\zigzag.exe dalla cartella del progetto)
    cargo run --bin zigzag -- --modo luce|elettrone|idrogeno [--stato 0..3] [--pannello] [--diagonale] [--cerchio]

* **1 luce**: un pacchetto di palline di spin 1 (λ = 8 celle) che corre a c lungo l'asse o la diagonale;
  si vedono le creste di (Re F)²; il pannello dà la velocità del centro misurata in quel momento.
* **2 elettrone**: palline arancioni (R) e blu (L), reticolo θ = 0,25; il pacchetto parte tutto destrorso e
  cambia colore a 2mc²/ħ; i surfisti (verdi) e ⟨r⟩ (bianco) fanno lo zig-zag lungo lo spin; con
  "spinore misto" girano in tondo. Nel pannello la **scacchiera di Feynman**: una sola pallina R al centro,
  il tempo in su. Metà della nube (davanti) è tagliata per vedere dentro.
* **3 idrogeno**: la nube dello stato scelto (1s, 2s, 2p₁/₂, 2p₃/₂), estratta da |ψ|² con le funzioni
  radiali **calcolate dal visore stesso** all'avvio (regola radiale, θ = 0,2 · 0,1 · 0,05, estrapolazione;
  ~12 s in un thread a parte); colore = fase dell'onda, che gira con la frequenza del livello (rallentata:
  1 s di schermo = 2·10⁻¹⁷ s); nucleo rosso; surfisti con le scie attorno all'asse dello spin. Il pannello
  mostra livelli, righe, struttura fine e surfista contro le misure, e gli ingressi.

Comandi: Tab pannello · Spazio pausa · 1/2/3 modo · R ricomincia · mouse: sinistro ruota, destro sposta,
rotella zoom. Le palline sono rotonde e morbide (shader proprio, `src/zigzag/balls.wgsl`); la mesh viene
**sostituita** a ogni fotogramma (Bevy 0.15 perde memoria se la si modifica sul posto).

Prove sul visore (Windows, `target\debug\zigzag.exe`, ~20 s per modo, memoria privata ogni 5 s):
| modo | memoria privata ogni 5 s (MB) | panic | immagine |
|---|---|---|---|
| luce (con pannello) | 578, 521, 521, 520 | no | `lab/out/zigzag_luce.png` |
| elettrone, zig-zag (con pannello) | 637, 638, 583, 579 | no | `lab/out/zigzag_elettrone.png` |
| elettrone, cerchio (60 s) | 578, 583, 588, 583, 580, 583, 584, 585, 584, 585, 583, 583 | no | `lab/out/zigzag_elettrone_cerchio.png` |
| idrogeno 1s (con pannello) | 643, 523, 523, 523 | no | `lab/out/zigzag_idrogeno.png` |
| idrogeno 2p₃/₂ | 642, 642, 522, 522 | no | `lab/out/zigzag_idrogeno_2p.png` |

La memoria resta piatta (il primo valore più alto è il calcolo dei livelli all'avvio). Il visore e il
laboratorio danno gli stessi livelli: 2p₁/₂ a θ = 0,2 = 3,4167292091 eV (Rust) contro 3,4167292090 (Python).
(Per arrivarci il prodotto scalare del quoziente di Rayleigh in Rust usa la somma compensata di Kahan:
con la somma semplice l'errore di arrotondamento, ~10⁻¹⁴ sulla fase, spostava i livelli di 2·10⁻⁸ eV e la
struttura fine di 6 MHz.)

`cargo test --bin zigzag` (8 test, 0 avvisi): unitarietà (regola 3-D di Dirac con potenziale, luce, radiale:
norma conservata a 10⁻¹²), matrice di Bloch = ciò che fa la regola a un'onda piana (10⁻¹²), relazione di
dispersione del reticolo (cos ω = cos θ cos ka esatta, energia a riposo esatta, errore ∝ θ²), velocità
≤ c (10⁵ spinori a caso; una pallina destrorsa con spin lungo x va esattamente una cella per tick),
idrogeno: 1s a θ = 0,2 uguale al laboratorio Python a 10⁻⁷ (13,7826491 eV) e struttura fine estrapolata =
Dirac a 10⁻³, armoniche sferiche spinoriali normalizzate.

---

## 8. Tabella dei test

Tolleranze nel dizionario `TOL` in testa a ogni script, fissate prima della prima esecuzione e mai
cambiate. "corsa" = quale esecuzione (le precedenti sono conservate nei `*_run*.log`).

| # | test | calcolato | atteso | tolleranza | esito | script (corsa, tempo) |
|---|---|---|---|---|---|---|
| L1a | dispersione della luce a λ = 8a | −3,53 % (diagonale) | 0 | 2 % | **FALLITO** | zz_light (3ª, 72 s) |
| L1b | errore ∝ (a/λ)² | rapporto 8a/16a = 4,09 | 4 | [3; 5] | PASS | zz_light |
| L2 ×6 | velocità del centro di un pacchetto largo 5a | 0,87–0,98 c | 1 e la velocità di gruppo dell'onda piana | 2 % / 0,5 % e 10⁻² | **FALLITO ×6** (test mal posto, §9) | zz_light |
| L3a | luce trasversale a 8a | 4,26·10⁻³ nel modo fermo | 0 | < 10⁻² | PASS | zz_light |
| L3b | migliora con λ/a | ×12,3 da 8a a 16a | — | ≥ 8 | PASS | zz_light |
| L4a | fase di Berry, cerchio massimo (spin 1) | e^{iγ} = +1,0000 | +1 | 10⁻² rad | PASS | zz_light |
| L4b | fase di Berry, cono di 60° (λ = 16a) | π + 0,082 | π | 10⁻² rad | **FALLITO** (0,005 a 64a) | zz_light |
| L5a | doppia fenditura, correlazione (λ = 8a, schermo 160) | 0,595 | ≥ 0,98 | | **FALLITO** | zz_slit (1ª, 295 s); a 6a in zz_light 1ª: 0,377 |
| L5b | errore delle frange chiare | 45 % rms | ≤ 8 % | | **FALLITO** | zz_slit |
| F1a | energia a riposo | 4,4·10⁻¹⁶ | 0 | 10⁻¹² | PASS | zz_free (3ª, 284 s) |
| F1b | E(p) fino a p = mc, θ = 0,1 | 0,194 % | 0 | 1 % | PASS | zz_free |
| F1c | errore ∝ θ² | rapporto 4,00 | 4 | [3; 5] | PASS | zz_free |
| F2a | frequenza ZBW 1-D (θ = 0,1) | 2,00951 | 2,00952 (continuo) | 1 % | PASS | zz_free |
| F2b | ampiezza ZBW 1-D | 0,49139 | 0,49050 | 2 % | PASS | zz_free |
| F2c | niente ZBW a energia positiva | 1,1·10⁻¹⁵ | 0 | 10⁻⁶ | PASS | zz_free |
| F3 | cerchio ZBW 3-D, regola contro continuo | 0,433/0,438 e 2,071 | 0,428/0,431 e 2,073 | 2 % | PASS | zz_free |
| F4a | fase di Berry spin ½, cerchio massimo | e^{iγ} = −1 (10⁻¹⁵) | −1 | 10⁻² rad | PASS | zz_free |
| F4a' | fase di Berry spin ½, cono 60° (64a) | −1,02·10⁻² rad | 0 | 10⁻² rad | **FALLITO** (di poco) | zz_free |
| F4b | legge del mezzo angolo | minimo a 90°, d(90)/d(180) = 0,101 | minimo a 90°, < 0,1 | | **FALLITO** (di poco) | zz_free |
| F5a | g regola contro continuo (σ = 8, θ = 0,5) | 1,820 contro 1,989 | uguali | 2·10⁻³ | **FALLITO** (θ², §4.4) | zz_free |
| F5b | g = 2 (σ = 8) | 1,820 | 2 | 1 % | **FALLITO** | zz_free |
| F6a | surfisti, velocità ≤ c | 1,000000000000000 | ≤ 1 | 10⁻¹² | PASS | zz_free |
| F6b | equivarianza 1-D | 1,26·10⁻² | 0 | 10⁻² | **FALLITO** | zz_free |
| H1 | 9 livelli contro Dirac (a → 0) | p, d ≤ 2·10⁻⁸; s 7,0·10⁻⁵ | 0 | 10⁻⁵ | **FALLITO** (gli s, §5.4) | zz_hydrogen (1ª, 249 s) |
| H2a | struttura fine contro Dirac | 10,94963 GHz | 10,94966 | 10⁻³ | PASS | zz_hydrogen |
| H2b | struttura fine contro misura | 10,94367 GHz | 10,96905 | 0,5 % | PASS (−0,23 %, QED) | zz_hydrogen |
| H3a | Lyman-α | 121,5780 nm | 121,567 | 10⁻⁴ | PASS (+90 ppm) | zz_hydrogen |
| H3b | Balmer-α | 656,4674 nm | 656,46 | 10⁻⁴ | PASS (+11 ppm) | zz_hydrogen |
| H4 | orbitali 1s, 2s, 2p | 1 − sovrapposizione ≤ 1,4·10⁻⁶ | 0 | 10⁻⁴ | PASS | zz_hydrogen |
| H5a | surfista 1s, θ = 0,05 | 1,000649 αc | αc | 10⁻³ | PASS | zz_hydrogen |
| H5b | surfista 1s, a → 0 | 0,9998957 αc | αc | 10⁻⁴ | **FALLITO** | zz_hydrogen |
| H6 | A(Lyman-α) | 6,2494·10⁸ s⁻¹ | 6,2649·10⁸ | 1 % | PASS | zz_hydrogen |
| S1 | lo spettro dell'ascolto ha le righe 1s e 2s | 0,018 / 0,29 eV dal risolvente | 0 | risoluzione 3,2 eV | PASS | zz_listen (3ª, 370 s) |
| S2 | 1s dall'ascolto = risolvente | 8,1·10⁻⁷ | 0 | 10⁻⁴ e decrescente | PASS | zz_listen |
| S3 | 2s dall'ascolto | −3,0·10⁻² | 0 | 10⁻² | **FALLITO** | zz_listen |
| A1 | regola 3-D, Z = 60, legame (a → 0) | −11,6 % | 0 | 1 % | **FALLITO** (§5.6) | zz_atom3d (1ª, 257 s) |
| A2 | regola 3-D, surfista | v_φ 0,94, v_r 0,37 | 1 e 0 | 5 % | **FALLITO** | zz_atom3d |
| R1–8 | `cargo test --bin zigzag` | 8 test | | | PASS ×8 | src/zigzag/tests.rs |

Totale nel laboratorio Python: **22 PASS, 18 FAIL** su 40 (contando L2 come 6). Più gli 8 test Rust, tutti PASS.

---

## 9. Che cosa non funziona, o funziona solo in parte

1. **La diffrazione della luce è troppo debole** (L5). È la proprietà più importante trovata qui: per
   questa regola a passi divisi cos ωΔt = cos k_x a · cos k_z a, quindi la curvatura trasversale della
   dispersione, che decide quanto un fascio si allarga, è sbagliata del fattore ka·cot(ka) (−21 % a λ = 8a,
   −5 % a 16a, −2,3 % a 24a). Lungo gli assi la luce va a c esattamente, ma una fenditura diffrange troppo
   poco finché λ < ~16a. La regola base del fotone (`REGOLA_BASE.md`) non ha questo problema perché usa
   le ondine esatte; qui ogni dimezzamento di a costa 8× in 2-D e 16× in 3-D. Lo stesso errore, per
   l'elettrone, è la massa cinetica (1 + θ²/3) m.
2. **Test scritti male da me.** L2 confrontava la velocità del centro di un pacchetto largo 5a con c, ma un
   pacchetto così stretto contiene molte direzioni e va più piano di c anche nel continuo; il confronto
   giusto (aggiunto dopo, senza cambiare il test) mostra che la regola fa esattamente ciò che dice la sua
   dispersione (0,4 %) e che rispetto a Maxwell è uguale lungo l'asse e più lenta del 9 % (8a) / 2,4 % (16a)
   in diagonale. L1a, L4b, F4a', F4b: tolleranze troppo strette per un errore del reticolo che va come
   (a/λ)²; F4b manca di 0,001 con il minimo nel posto giusto. F5: ho provato g con θ = 0,5 senza pensare
   alla massa cinetica; `zz_gfactor.py` mostra che g(regola) → g(Dirac) come θ².
3. **Gli stati s vicino al nucleo** (H1, H5b): l'errore del reticolo ha un termine di contatto θ² ln θ;
   con reticoli fino a θ = 0,0125 l'1s arriva a 9,2·10⁻⁶ (3,0·10⁻⁶ con il termine logaritmico). Per questo
   il 1S–2S della regola (−19 GHz) non può ancora vedere lo spostamento di Lamb (+7 GHz).
4. **La regola 3-D sull'atomo** (A1, A2) è provata solo a Z = 60 e su reticoli grossi: lega l'elettrone e
   concorda con la forma radiale entro 2–5 %, ma per arrivare all'1 % servirebbero ~480³ celle.
   **L'idrogeno vero (Z = 1) è calcolato solo nella forma radiale** (esatta per simmetria nel continuo).
5. **Ascolto**: l'1s esce dall'ascolto delle palline a 8·10⁻⁷ dal risolvente; il 2s no (−3 %) con 10⁶ tick.
6. **Il lettore dei surfisti** (F6b): interpolazione lineare + RK4 su due tick: errore 1,3·10⁻² sul "posto in
   fila". Il surfista di `SURFISTA_DIRAC.md` legge un campo spettrale e fa 1,5·10⁻⁵.
7. **Errori di codice, corretti (registri conservati).** (a) `zz_light` 1ª corsa: pacchetto costruito con
   autovettori di fase casuale (un pacchetto "sparpagliato"): corretto con il proiettore; (b) `zz_free` 1ª
   corsa: `dirac_velocity` restituiva la corrente invece della velocità (g ≈ 0); (c) `zz_free` 2ª corsa: due
   autovettori degeneri non ortonormali usati come proiettore (S = 0,473 invece di 0,498); (d) `zz_listen`
   1ª e 2ª corsa oltre i 10 minuti, e righe false dovute ai rimbalzi sul bordo: reticolo più grosso e spugna;
   (e) `zz_slit_conv` 1ª corsa: il punto λ = 24a è stato interrotto dal limite di tempo prima di scrivere il
   json. La doppia fenditura era nata dentro `zz_light.py` (1ª corsa: λ = 6a, correlazione 0,377), poi
   spostata in `zz_slit.py` (numba) perché a 8a in numpy servivano 753 s.
8. **Non fatto.** L'accoppiamento **dinamico** della carica ai camminatori della luce (fasi di Peierls sui
   collegamenti): c'è solo il potenziale statico del nucleo. Quindi niente emissione "vera" di luce
   dall'atomo nella regola: le righe sono differenze di livelli e le vite vengono dalla formula di dipolo
   applicata alle funzioni d'onda della regola. Niente campo magnetico esterno, niente più elettroni
   (servirebbe una fogliazione, come in `MULTIFOTONE.md`), niente rinculo oltre la massa ridotta.
9. **Che cosa la regola non può dare comunque**: QED (Lamb, g − 2, i 25 MHz della struttura fine), la
   creazione di coppie, il perché delle masse e di α, il protone. La Zitterbewegung dell'elettrone libero
   non è mai stata osservata direttamente.

---

## 10. Come rieseguire

Dalla cartella del progetto (seed fissi; tempi sulla macchina condivisa):

    python lab/zz_core.py            # autotest della macchina (2 s)
    python lab/zz_fast.py            # autotest dei nuclei numba contro numpy
    timeout 590 python lab/zz_light.py          # T1: dispersione, velocità, trasversalità, Berry (72 s)
    timeout 590 python lab/zz_slit.py           # T1: doppia fenditura a lambda = 8a (295 s)
    timeout 590 python lab/zz_slit_conv.py      # T1: convergenza della doppia fenditura, schermo a 40 (318 s)
    timeout 590 python lab/zz_free.py           # T2 (284 s)
    timeout 590 python lab/zz_gfactor.py        # T2: g contro theta (40 s)
    timeout 590 python lab/zz_hydrogen.py       # T3: livelli, righe, orbitali, surfista (249 s)
    timeout 590 python lab/zz_hydrogen_s.py     # T3: gli s su reticoli più fini (46 s)
    timeout 590 python lab/zz_listen.py         # T3: ascoltare l'atomo (370 s)
    timeout 590 python lab/zz_atom3d.py         # T3: la regola 3-D con Z = 60 (257 s)
    python lab/zz_atom3d_radial.py               # T3: la stessa regola radiale a Z = 60 (5 s)
    cargo test --bin zigzag
    cargo run --bin zigzag

Nota: sotto Windows `timeout` ferma il lanciatore `python.exe` di WindowsApps ma non sempre il processo
figlio; i tempi riportati sono quelli scritti dagli script.
