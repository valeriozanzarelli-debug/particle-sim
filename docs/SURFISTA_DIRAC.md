# Il surfista di Dirac: l'elettrone con la regola "a due ruoli"

Laboratorio: `lab/dirac_*.py` (Python 3.14 + numpy/scipy). Unità naturali ħ = c = m_e = 1, cioè
lunghezza = λ̄_C = ħ/(m c) = 3,8616·10⁻¹³ m, tempo = ħ/(m c²) = 1,2881·10⁻²¹ s, energia = m c² = 510 998,95 eV.
Tutti i numeri di questo documento vengono da esecuzioni reali (seed fissi, macchina a 4 core condivisa con
altri processi). Registri in `lab/out/dirac_*.log`, risultati in `lab/out/dirac_*.json`, figure
`lab/out/dirac_zbw_1d.png`, `lab/out/dirac_zbw_loop.png`, `lab/out/dirac_h1s_circulation.png` (in tutto 0,34 MB).

---

## 0. In breve

**L'idea di Valerio.** L'elettrone = la regola a due ruoli (palline-onda + una pallina "surfista" guidata dalla
corrente), più un fattore geometrico che fa cadere il surfista in uno stato semi-stabile in cui gira in tondo
velocissimo; l'onda fa un'orbita strettissima e così, in media, il surfista è "quasi fermo".

**Che cosa ne dice la fisica (e i nostri numeri).**

* Il pezzo "gira in tondo velocissimo" **esiste già** ed è dentro l'equazione di Dirac: è la
  *Zitterbewegung* (Schrödinger 1930). Con la regola di Bohm–Dirac `v = ψ†αψ / ψ†ψ` l'abbiamo vista
  direttamente: surfisti che girano su cerchi di raggio **0,486 ħ/mc = 1,88·10⁻¹³ m** a velocità
  **0,99993 c**, con frequenza **2,06 mc²/ħ ≈ 1,6·10²¹ rad/s** (atteso: ħ/2mc = 1,93·10⁻¹³ m, c, 2mc²/ħ = 1,55·10²¹ rad/s).
  La deriva media è la velocità di gruppo, la velocità non supera mai c, e una popolazione di surfisti
  distribuita come |ψ|² resta distribuita come |ψ|².
* **Ma**: questi giri esistono solo se l'onda contiene componenti a **energia negativa** (un elettrone
  localizzato in meno di ~ħ/mc). Un elettrone "normale" (energia positiva) non ha Zitterbewegung: misurato
  residuo < 10⁻¹⁵.
* **Nell'atomo di idrogeno** (stato 1s, risolto esattamente con Dirac) il surfista **non** è quasi fermo e
  **non** fa un'orbita strettissima: gira lentamente attorno all'asse dello spin, su cerchi di raggio
  r sinθ (≈ raggio di Bohr), a velocità **esattamente α c sinθ = 2,19·10⁶ m/s · sinθ**, uguale a ogni
  distanza dal nucleo. Nella vecchia regola non relativistica (Schrödinger) lo stesso surfista è fermo.
* **Spin e g = 2** escono dalla corrente di Dirac (calcolato: g = 1,996 per un pacchetto di 8 ħ/mc, → 2 al
  crescere del pacchetto). **Non** escono dal "cerchietto classico" da solo: una carica che gira a c sul
  raggio ħ/2mc ha momento angolare ħ/2 ma momento magnetico eħ/4m, cioè g = 1 (§4). Il g − 2 della QED
  (0,00115965…) non c'è in nessuno dei due quadri.
* **Idrogeno**: energie di Dirac riprodotte dal nostro risolutore numerico con errore relativo ≤ 5·10⁻¹¹;
  struttura fine 2p₃/₂ − 2p₁/₂ = **10,9437 GHz** contro **10,9690 GHz** misurati (−0,23 %: i 25 MHz che
  mancano sono QED); Lamb shift 2s₁/₂ − 2p₁/₂: Dirac dà **0** (misurato 1057,8 MHz).

Esito dei test: **34 PASS, 2 FAIL** su 36 (tutti e due "fallimenti di previsione" spiegati in §8; nessuna
tolleranza cambiata dopo aver visto i risultati).

---

## 1. Da dove viene l'idea, e chi l'ha già avuta

| idea | chi | che cosa dice |
|---|---|---|
| Zitterbewegung ("moto tremolante") | E. Schrödinger, *Sitzungsber. Preuss. Akad. Wiss.* 24, 418 (1930) | nell'equazione di Dirac la "velocità" dell'elettrone è c·α, i cui autovalori sono ±c; la posizione media oscilla a 2mc²/ħ con ampiezza ~ħ/2mc |
| il surfista di Dirac | D. Bohm, *Prog. Theor. Phys.* 9, 273 (1953); P. R. Holland, *The Quantum Theory of Motion* (Cambridge 1993), cap. 12 | la particella si muove con v = c ψ†αψ/ψ†ψ; |v| ≤ c sempre |
| il caso a molte particelle | Dürr, Goldstein, Münch-Berndl, Zanghì, *Phys. Rev. A* 60, 2729 (1999) | serve una "fogliazione" dello spazio-tempo (un "adesso" comune), come nella nostra regola a molti fotoni |
| l'elettrone come carica che gira a c | D. Hestenes, *Found. Phys.* 20, 1213 (1990); *Found. Phys.* 40, 1 (2010) (arXiv:0802.2728) | "zitter": elica percorsa a velocità c, raggio λ_e = ħ/2mc ≈ 1,9308·10⁻³ Å, frequenza 2mc²/ħ |
| la corrente si divide in "moto" + "spin" | W. Gordon, *Z. Phys.* 50, 630 (1928) | j = (parte convettiva) + (1/2m)∇×(densità di spin): la seconda è una corrente che circola attorno allo spin |

Quindi la regola da costruire non va inventata: è la dinamica di Bohm–Dirac. Il nostro lavoro è (1) scriverla
come "cosa fa ogni pallina a ogni tick", (2) controllare che faccia davvero quello che si dice, (3) vedere che
cosa dà nell'idrogeno, (4) confrontare con le misure.

---

## 2. La regola (congelata il 2026-10-09 in `lab/dirac_core.py`, prima di ogni test)

### 2.1 Palline-onda: portano uno spinore

Nella regola del fotone ogni pallina-onda porta un vettore di Jones (2 numeri complessi). Per l'elettrone porta
uno **spinore di Dirac**: 4 numeri complessi in 3+1 D; 2 numeri complessi in 1+1 D e in 2+1 D (le dimensioni
che simuliamo qui). L'onda evolve con l'equazione di Dirac libera

    i ∂ψ/∂t = (α·p + β m) ψ          1+1 D: α = σx, β = σz;   2+1 D: α = (σx, σy), β = σz

Numericamente: **split-operator spettrale** (FFT). Per ogni onda piana k il propagatore è esatto,
`U(k,t) = cos(E t) − i sin(E t) H(k)/E`, `E = √(k² + m²)`. Niente "fermion doubling": ogni k della griglia ha
il suo E(k) vero (il raddoppio nasce dalle derivate alle differenze finite, che qui non usiamo).
Autotest: unitarietà 2·10⁻¹⁶, U(1,2)U(2,5) − U(3,7) = 1·10⁻¹⁶, proiettori P₊ + P₋ = 1 a 7·10⁻¹⁷.

### 2.2 Il surfista: che cosa fa a ogni tick

1. **Nasce** in un punto scelto con probabilità ψ†ψ (l'unico numero casuale della regola).
2. **Legge lo spinore** ψ = (a, b) dove si trova (interpolazione cubica delle componenti dalle palline-onda
   vicine).
3. **Calcola la sua velocità**

       v = c · ψ†αψ / ψ†ψ        in 2+1 D:  vx = 2 Re(a* b)/(|a|²+|b|²),  vy = 2 Im(a* b)/(|a|²+|b|²)

   È sempre |v| ≤ c, per **qualunque** spinore (disuguaglianza di Cauchy–Schwarz): non serve nessun freno.
   Controllato su 10⁵ spinori casuali: max |v| = 0,99999999996.
4. **Si muove** con un passo RK4 di durata DT = 0,02 ħ/mc² (157 passi per periodo di Zitterbewegung),
   leggendo l'onda negli istanti esatti degli stadi intermedi.
5. Fine. Non c'è altro: niente monete, niente "fattore geometrico".

Parametri congelati: DT = 0,02, interpolazione cubica, seed 20261009. Due **accelerazioni** aggiunte dopo
(non cambiano i numeri: verificate nell'autotest): coefficienti spline calcolati una volta per campo
(differenza 0,0) e propagazione a mezzi passi U(dt/2) invece che da t = 0 (differenza 3,8·10⁻¹⁶).

### 2.3 Che cosa cambia rispetto alla regola non relativistica, e perché è la stessa a bassa velocità

| | regola non relativistica (Schrödinger) | regola di Dirac |
|---|---|---|
| pallina-onda porta | 1 numero complesso | uno spinore (2 in 2+1 D, 4 in 3+1 D) |
| velocità del surfista | v = (ħ/m) Im(ψ*∇ψ)/|ψ|² | v = c ψ†αψ/ψ†ψ |
| limite di velocità | nessuno | |v| ≤ c automatico |
| stato fermo (onda reale) | surfista **fermo** | surfista che **gira** attorno allo spin |

A bassa velocità (v ≪ c) la componente "piccola" dello spinore vale χ ≈ (σ·p)φ/(2mc), e sostituendola in
ψ†αψ si trova (decomposizione di Gordon, fatta a mano in 2+1 D):

    v = (ħ/m) [ Im(φ*∇φ) + ½ ∇×(|φ|² ẑ) ] / |φ|²          (regola di Pauli)

cioè **la vecchia regola + un termine di spin**, un vortice attorno all'asse dello spin che non porta via
nulla (ha divergenza nulla, quindi non cambia dove arrivano i surfisti in media) ma li fa ruotare.
Test (`dirac_pauli.py`): stessi punti di partenza, tre regole, ε = v/c = 1/σ:

| ε = k₀ = 1/σ | tempo T = σ²/2 | |x_Dirac − x_Pauli| / σ (RMS) | |x_Dirac − x_Schrödinger senza spin| / σ |
|---|---|---|---|
| 0,200 | 12,5 | 3,21·10⁻² | 0,334 |
| 0,141 | 25 | 1,73·10⁻² | 0,356 |
| 0,100 | 50 | 8,87·10⁻³ | 0,345 |

La differenza Dirac–Pauli si dimezza a ogni passo √2 di ε (rapporti 1,86 e 1,95; atteso 2): è di
**secondo ordine in v/c**, come deve. La differenza con la regola di Schrödinger senza spin **non** scende:
resta un terzo della larghezza del pacchetto. Conclusione: la regola di Dirac sostituisce quella vecchia,
e a bassa velocità si riduce alla vecchia **più il termine di spin**, che è necessario.

---

## 3. Elettrone libero: Zitterbewegung e surfisti (`dirac_free.py`, `dirac_surf2d.py`)

Pacchetto gaussiano di larghezza σ (σ = 5 ħ/mc in 1+1 D, σ = 4 in 2+1 D), spinore (1, 1)/√2. Questo spinore a
impulso zero è **metà energia positiva e metà negativa** (misurato: peso positivo 0,5000). Lo stesso pacchetto
proiettato con P₊ = (1 + H/E)/2 è un elettrone "normale".

### 3.1 Zitterbewegung di ⟨x⟩ (figura `dirac_zbw_1d.png`, sinistra)

* misto, 1+1 D: oscillazione con **ω = 2,0095 mc²/ħ = 1,560·10²¹ rad/s** (2⟨E⟩ = 2,0099 per questo pacchetto;
  2mc²/ħ = 2) e **ampiezza 0,4905 ħ/mc = 1,894·10⁻¹³ m** (ħ/2mc = 1,931·10⁻¹³ m);
* solo energia positiva: residuo **2,2·10⁻¹⁵** (zero numerico);
* 2+1 D, misto fermo: ⟨r⟩(t) percorre un **cerchio** (figura `dirac_zbw_loop.png`, sinistra): ampiezze
  0,4752 (x) e 0,4745 (y), ω = 2,029. È proprio il "cerchietto di Hestenes" di raggio ħ/2mc. L'ampiezza è un
  po' sotto 0,5 perché il pacchetto è stretto (σ = 4); al crescere di σ va a 0,5 (`dirac_zbw_scan.py`):

  | σ [ħ/mc] | 2 | 4 | 8 | 16 |
  |---|---|---|---|---|
  | ampiezza | 0,3695 | 0,4752 | 0,4954 | 0,4990 |
  | ω [mc²/ħ] | 2,075 | 2,029 | 2,008 | 2,002 |

### 3.2 I surfisti

* **Velocità ≤ c**: su tutti gli stadi RK4 di tutti i surfisti (2000 × 2000 passi in 1D, 3000 × 1200 in 2D e
  le altre corse) max |v| = 1 − (≤ 4·10⁻¹⁶): il limite viene toccato ma mai superato oltre l'arrotondamento.
* **Anelli**: i 20 surfisti nati più vicini al centro del pacchetto misto 2D, nei primi due periodi, girano su
  cerchi di raggio **0,486 ħ/mc = 1,88·10⁻¹³ m**, con ω = 2,06 mc²/ħ e velocità media **0,99993 c**
  (figura `dirac_zbw_loop.png`, centro). In 1D il surfista centrale va avanti e indietro con ampiezza 0,494 e
  ω = 2,019 (zig-zag nella figura `dirac_zbw_1d.png`, destra). Con un impulso k₀ = 0,5 mc i cerchi si
  trascinano: cicloidi nel piano, eliche nello spazio-tempo (figura `dirac_zbw_loop.png`, destra).
  **Questo è esattamente il "surfista che gira velocissimo in un'orbita strettissima" dell'idea.**
* **Deriva = velocità di gruppo**: pacchetti a energia positiva con k₀ = 0,5 / 1 / 2 (1D): deriva dei
  surfisti 0,44656 / 0,70518 / 0,89394 c contro ⟨k/E⟩ = 0,44656 / 0,70521 / 0,89393 c (v_g(k₀) = 0,447 /
  0,707 / 0,894). In 2D (k₀ = 0,5): 0,4448 contro 0,4436. Pacchetti misti: la deriva è la media
  pesata dei due rami (il ramo a energia negativa va all'indietro), 0,192 / 0,476 / 0,803 c contro
  0,201 / 0,498 / 0,799 c in 1D e 0,198 contro 0,200 in 2D (la differenza è il resto della
  Zitterbewegung della media, ≤ 0,5/T, dichiarato nella tolleranza).
* **Equivarianza** (i surfisti nati con |ψ|² restano distribuiti come |ψ|²):
  * 1D: in Bohm 1D le traiettorie non si incrociano, quindi ogni surfista deve conservare il suo "posto in
    fila" F = ∫_{−∞}^{x} |ψ|². Misurato: max |F(x_i(t), t) − F(x_i(0), 0)| = **1,5·10⁻⁵** su 2000 surfisti
    per 40 ħ/mc² (13 periodi di Zitterbewegung). Media dei surfisti contro ⟨x⟩(t) del campo: scarto massimo
    0,051 (errore statistico 0,147).
  * 2D: χ² dell'istogramma dei 3000 surfisti contro |ψ(t)|²: t = 0: 102,5/87 (p = 0,12); t = 12:
    73,2/89 (p = 0,89); t = 24: 101,7/110 (p = 0,70). La media dei surfisti segue il cerchio di ⟨r⟩(t) con
    scarto massimo 0,087 (errore statistico 0,090).

### 3.3 Che cosa vuol dire

Il surfista di Dirac fa davvero i giri veloci dell'idea, con il raggio (ħ/2mc) e la frequenza (2mc²/ħ) giusti,
e in media va alla velocità di gruppo, cioè sta "quasi fermo" se il pacchetto è fermo. Ma la causa dei giri
non è un fattore geometrico esterno: è la **interferenza tra energia positiva e negativa** dentro l'onda.
Togli l'energia negativa (P₊) e i giri spariscono (residuo 10⁻¹⁵). Per avere energia negativa in quantità
serve un elettrone compresso sotto ~ħ/mc = 3,9·10⁻¹³ m; in un atomo non succede (§5).

---

## 4. Spin e momento magnetico dal "cerchietto"

### 4.1 La derivazione, in parole semplici

Prendi una carica e che gira a velocità c su un cerchio di raggio r = ħ/2mc.

1. **Frequenza.** Un giro è lungo 2πr, quindi ω = c/r = **2mc²/ħ**: è la frequenza della Zitterbewegung.
   Coerente.
2. **Momento angolare.** Se ciò che gira porta l'intera energia mc² (come un fotone, che a velocità c ha
   impulso p = E/c = mc), allora L = r·p = (ħ/2mc)·mc = **ħ/2**. È lo spin dell'elettrone. Coerente.
3. **Momento magnetico.** Una carica che gira è una spira di corrente I = eω/2π, di area πr²:
   μ = I·A = eωr²/2 = e·c·r/2 = **eħ/4m**, cioè **metà** del magnetone di Bohr eħ/2m. Il rapporto
   g = (2m/e)(μ/L) = **1**, non 2.

Quindi la frase "una carica che gira a c sul raggio ħ/2mc ha spin ħ/2 e momento eħ/2m" è giusta per lo spin
ma **sbagliata di un fattore 2 per il momento magnetico**, se si usa la fisica classica della spira. Per avere
g = 2 con lo stesso cerchio bisogna assumere qualcos'altro: per esempio che solo metà dell'energia mc² stia
nel moto circolare (allora p = mc/2 e servirebbe r = ħ/mc per avere ancora L = ħ/2, ma allora la velocità
sarebbe 2c), oppure mettere il momento magnetico a mano. Hestenes lo dice esplicitamente nel lavoro del 2010:
la sua lagrangiana del modello "zitter" ha costanti di accoppiamento separate per il potenziale A e per il
campo F, e "appears to presume the electron magnetic moment" (Hestenes 2010, arXiv:0802.2728); nella teoria di Dirac
invece g = 2 è una conseguenza.

**Dove sta davvero il fattore 2.** Nella corrente di Dirac. La decomposizione di Gordon dice

    j = j_moto + (1/m) ∇×s          s = densità di spin (ħ/2 per elettrone)

e per una distribuzione localizzata ∫ r × (∇×s) d³r = **2** ∫ s d³r (integrazione per parti; vale anche in
2D). Quindi il momento magnetico della corrente di spin è μ = (e/2)∫ r × j = (e/m)·S = 2·(e/2m)·S: **g = 2**.
Il fattore 2 viene dall'integrazione per parti di un vortice, non da una spira: il vortice dello spin
"pesa" due volte nel momento magnetico rispetto a quanto peserebbe come massa che gira. (Corollario curioso:
se si interpretasse la corrente di spin come massa in moto, porterebbe momento angolare ħ, non ħ/2: la densità
di impulso di Dirac **non** è m·j.)

### 4.2 Il conto numerico (2+1 D, `dirac_free.py`, test G1)

Pacchetto a energia positiva, spinore (1, 0), fermo; S = ⟨σz⟩/2, L = ⟨−i(x∂y − y∂x)⟩, μ = ½∫(x j_y − y j_x)
(con e = m = 1), g = 2μ/S:

| σ [ħ/mc] | ⟨S_z⟩ | ⟨L_z⟩ | J = S + L | μ [eħ/m] | g |
|---|---|---|---|---|---|
| 2 | 0,47386 | 0,02614 | 0,500000 | 0,46200 | 1,9499 |
| 4 | 0,49258 | 0,00742 | 0,500000 | 0,48898 | 1,9854 |
| 8 | 0,49807 | 0,00193 | 0,500000 | 0,49712 | **1,9962** |

J = ½ esatto (conservato); g → 2 come 1/σ² (per pacchetti più stretti una parte del momento angolare passa da
spin a orbitale). I surfisti di questo pacchetto **ruotano** attorno al centro (test G2): velocità
angolare 0,0272 mc²/ħ, velocità 0,109 c a r = σ (3,3·10⁷ m/s); vedi §8 per il confronto con la previsione.

### 4.3 Che cosa questo quadro NON dà: g − 2

L'elettrone vero ha g/2 = 1,001 159 652 180 59 (13) (Fan et al., *Phys. Rev. Lett.* 130, 071801 (2023)),
cioè l'anomalia a_e = (g − 2)/2 = **0,00115965**. Né il cerchietto né Dirac la producono: Dirac dà g = 2
esatto. L'anomalia viene dalla **elettrodinamica quantistica**: l'elettrone emette e riassorbe fotoni virtuali
(correzione di vertice); il primo termine è α/2π = 0,0011614 (Schwinger, *Phys. Rev.* 73, 416 (1948)), il resto
da termini in (α/π)², (α/π)³, … e da piccoli contributi adronici ed elettrodeboli. Nel linguaggio del progetto:
serve il **vuoto quantistico** (punto 4 dell'ordine di lavoro in `VISIONE.md`), non una geometria migliore.

---

## 5. Idrogeno con la regola di Dirac — il test decisivo (`dirac_hydrogen.py`)

### 5.1 Energie

Equazioni radiali di Dirac–Coulomb (Z = 1), risolte con tiro (shooting) dall'origine e dall'infinito e
raccordo del wronskiano, DOP853 in ln r con tolleranza 10⁻¹³. La formula di Sommerfeld–Dirac è usata **solo**
per controllare, mai come ingresso.

| livello | E − mc² numerico [eV] | formula di Dirac [eV] | differenza relativa |
|---|---|---|---|
| 1s₁/₂ | −13,605874258 | −13,605874258 | 5,2·10⁻¹² |
| 2s₁/₂ | −3,401479886 | −3,401479886 | 7,1·10⁻¹² |
| 2p₁/₂ | −3,401479886 | −3,401479886 | 7,1·10⁻¹² |
| 2p₃/₂ | −3,401434601 | −3,401434601 | 8,6·10⁻¹² |
| 3d₅/₂ | −1,511745917 | −1,511745917 | 4,5·10⁻¹¹ |

Confronto con le misure (massa ridotta: E × μ/m, μ/m = 1/(1 + m_e/m_p)):

| quantità | Dirac (questo lavoro) | misura | differenza | da dove viene la differenza |
|---|---|---|---|---|
| energia di legame 1s, massa infinita | 13,605874 eV = 13,605693 (Bohr) + **1,811·10⁻⁴** (relatività) | — | — | — |
| energia di ionizzazione 1s, massa ridotta | 13,598468 eV | 13,598 434 600 eV (NIST) | +3,37·10⁻⁵ eV = **8,15 GHz** | Lamb shift dell'1s (8,17 GHz, QED) |
| 2p₃/₂ − 2p₁/₂ | 10,94965 GHz (massa inf.), **10,94369 GHz** (ridotta) | **10,96905 GHz** (9911,200 MHz [Hagley–Pipkin 1994] + 1057,845 MHz [Lundeen–Pipkin 1981]) | −25,4 MHz (−0,23 %) | QED: soprattutto il momento magnetico anomalo a_e, che sposta 2p₁/₂ e 2p₃/₂ in direzioni opposte |
| 2s₁/₂ − 2p₁/₂ (Lamb shift) | **0** (3,7·10⁻⁶ MHz, zero numerico) | 1057,845 MHz (Lundeen–Pipkin 1981; 1057,829 Bezginov 2019) | tutta | **non è in Dirac**: autoenergia e polarizzazione del vuoto (Lamb–Retherford 1947) |
| 1S–2S | 2,466 068 519·10¹⁵ Hz | 2,466 061 413 187·10¹⁵ Hz (Parthey 2011) | +7,11 GHz | Lamb shift 1s − Lamb shift 2s ≈ 8,17 − 1,05 = 7,1 GHz |

Nota onesta: la frase "Dirac dà 10,969 GHz" non è esatta. Dirac dà 10,944–10,950 GHz; i 10,969 GHz misurati
contengono già 25 MHz di QED. La regola di Dirac spiega il 99,77 % della struttura fine; il resto, e tutto il
Lamb shift, richiedono il vuoto quantistico.

### 5.2 Il surfista nello stato 1s (figura `dirac_h1s_circulation.png`)

Per lo stato 1s con spin su (κ = −1, m_j = +½), con le funzioni radiali G (grande) e F (piccola), la regola
dà una velocità **solo azimutale**:

    v = −c · 2GF/(G² + F²) · sinθ · φ̂

Per l'1s il rapporto F/G è costante: F/G = −√((1−γ)/(1+γ)), γ = √(1 − α²) (numerico: −0,0036487249 a
tutte le r, teoria −0,0036487249). Allora

    |v| = c √(1 − γ²) sinθ = **α c sinθ**          (esatto, indipendente dalla distanza dal nucleo)

Numerico: |v|/(αc) all'equatore tra 0,99999999999994 e 1,0000000000000013 su r da 10⁻⁴ a 4,8 raggi di Bohr.

| | valore |
|---|---|
| traiettoria | cerchi attorno all'asse dello spin (z e distanza dall'asse conservati: 0 e 5,6·10⁻¹¹ a₀ su un giro, 2000 surfisti RK4) |
| velocità | α c sinθ = **2,188·10⁶ m/s × sinθ** (la velocità di Bohr!), media sulla nube α c π/4 = 0,0057 c |
| raggio dell'orbita | r sinθ: tipicamente un raggio di Bohr, 5,3·10⁻¹¹ m (non ħ/2mc = 1,9·10⁻¹³ m) |
| velocità angolare | α c / r: a r = a₀ **4,13·10¹⁶ rad/s** (un giro in 1,52·10⁻¹⁶ s); sale come 1/r verso il nucleo |
| senso | antiorario attorno a +z per spin su |
| regola di Pauli (Schrödinger 1s + termine di spin) | |v| = α c sinθ: **identica** |
| regola di Schrödinger senza spin | v = 0: il surfista è **fermo** (funzione d'onda reale) |
| stato 2s | stessa forma, ma vicino al nodo di G (r ≈ 2 a₀) dove |G| = |F| la velocità sale fino a c sinθ (sulla nostra griglia: 0,935 c) |

**Risposta diretta all'idea.** Con la regola di Dirac il surfista dell'elettrone nell'atomo **non è fermo**: gira
attorno all'asse dello spin a velocità α c sinθ. Il moto è dovuto allo **spin**, non al momento angolare orbitale
(che nell'1s è zero), ed è lo stesso "vortice di spin" della §4 visto dentro un atomo. Ma **non** è l'orbita
ultra-stretta della Zitterbewegung: niente raggio ħ/2mc, niente 10²¹ rad/s, perché uno stato legato stazionario
non contiene energia negativa. Il "quasi fermo in media" è vero in un senso preciso: la velocità media della
nube è zero (sono cerchi), e il surfista non va mai da nessuna parte.

---

## 6. Il "fattore geometrico": che cosa dovrebbe fare senza l'equazione di Dirac

Se si volesse ottenere il cerchietto stabile solo da Maxwell + geometria, il fattore geometrico dovrebbe:
(1) tenere una carica su un'orbita chiusa di raggio ħ/2mc a velocità c senza che irraggi (una carica classica
che gira a 2mc²/ħ irraggerebbe tutta la sua energia in ~10⁻²³ s); (2) fissare il momento angolare a ħ/2 per
**ogni** elettrone, qualunque cosa succeda; (3) dare g = 2 e non 1 (§4.1); (4) produrre la statistica di Fermi
e l'antiparticella. Due tentativi noti della letteratura vanno in questa direzione. La **soluzione di
Kerr–Newman** della relatività generale (buco nero carico e rotante) ha g = 2 esattamente (B. Carter,
*Phys. Rev.* 174, 1559 (1968)); con i parametri dell'elettrone il momento angolare per unità di massa è
a = ħ/2mc (proprio il raggio della zitter), molto più grande della massa geometrica, quindi niente orizzonte ma
una singolarità ad anello nuda di raggio ħ/2mc: è il programma di A. Burinskii (per es. *Grav. Cosmol.* 14, 109
(2008)), che però non deriva la quantizzazione né l'equazione di Dirac. La **teoria di Einstein–Cartan** (con
torsione) accoppia la torsione dello spazio-tempo alla densità di spin dei fermioni (F. W. Hehl, P. von der
Heyde, G. D. Kerlick, J. M. Nester, *Rev. Mod. Phys.* 48, 393 (1976)), ma presuppone già i campi di Dirac e i
suoi effetti diventano importanti solo a densità enormi (~10⁵⁴ g/cm³ per gli elettroni). In nessuno dei due
casi il fattore geometrico sostituisce l'equazione di Dirac: nel primo dà il g giusto ma non la quantistica,
nel secondo parte dalla quantistica.

---

## 7. Tabella di validazione

Tolleranze fissate nel codice (dizionario `TOL` o in testa allo script) prima della prima corsa. Comando:
`timeout 590 python lab/<script>` dalla cartella del progetto. Tempi su macchina condivisa.

| # | test | calcolato | atteso / misurato | tolleranza | esito | script (tempo) |
|---|---|---|---|---|---|---|
| Z1a | frequenza ZBW 1D | 2,0095 mc²/ħ = 1,560·10²¹ rad/s | 2mc²/ħ = 1,553·10²¹ rad/s (2⟨E⟩ = 2,0099) | 2 % | PASS | dirac_free.py (203 s) |
| Z1b | ampiezza ZBW 1D | 0,4905 ħ/mc = 1,894·10⁻¹³ m | ħ/2mc = 1,931·10⁻¹³ m | 5 % | PASS | dirac_free.py |
| Z1c | niente ZBW con sola energia positiva, 1D | 2,2·10⁻¹⁵ | 0 | 10⁻⁶ | PASS | dirac_free.py |
| Z1d | ZBW 2D, pacchetto fermo (cerchio) | A = 0,4752 / 0,4745, ω = 2,029 | 0,5, 2 | 5 %, 2 % | **FAIL** (ampiezza x −4,96 %, y −5,10 %; ω +1,4 %) | dirac_free.py |
| Z1d' | ZBW 2D con k₀ = 0,5 | ω = 2,233 | ≈ 2⟨E⟩ (solo controllo di presenza) | 20 % | PASS | dirac_free.py |
| Z1e | niente ZBW con sola energia positiva, 2D | 1,0·10⁻¹⁵ | 0 | 10⁻⁶ | PASS | dirac_free.py |
| S1a,b | |v| ≤ c in ogni stadio RK4 | 1 + 4·10⁻¹⁶ al massimo | ≤ 1 | 10⁻¹² | PASS (×8 corse) | free, surf2d |
| S2a | posto in fila conservato, 1D | 1,5·10⁻⁵ | 0 | 10⁻³ | PASS | dirac_free.py |
| S2b | media surfisti = ⟨x⟩ del campo, 1D | 0,051 | errore stat. 0,147 | 4σ | PASS | dirac_free.py |
| S2c | χ² surfisti vs |ψ|², 2D, t = 0 / 12 / 24 | p = 0,12 / 0,89 / 0,70 | — | p > 10⁻³ | PASS ×3 | dirac_surf2d.py (246 s) |
| S2d | media surfisti = ⟨r⟩ del campo, 2D | 0,087 | errore stat. 0,090 | 4σ | PASS | dirac_surf2d.py |
| S3 | deriva = velocità di gruppo, 1D, k₀ = 0,5/1/2, positivo | 0,44656 / 0,70518 / 0,89394 | 0,44656 / 0,70521 / 0,89393 | 1 % + 4σ | PASS ×3 | dirac_free.py |
| S3 | deriva, 1D, misto | 0,1919 / 0,4758 / 0,8032 | 0,2014 / 0,4975 / 0,7991 | 1 % + 4σ + 0,5/T | PASS ×3 | dirac_free.py |
| S3 | deriva, 2D, k₀ = 0,5, positivo / misto | 0,4448 / 0,1985 | 0,4436 / 0,1999 | idem | PASS ×2 | dirac_surf2d.py |
| G1 | g dalla corrente di Dirac (σ = 8) | 1,9962 (J = 0,500000) | 2 | 1 % | PASS | dirac_free.py |
| G2 | rotazione dei surfisti, pacchetto con spin | Ω = 0,02721 | 1/(2σ²) = 0,03125 | 5 % | **FAIL** (−13 %; vedi §8) | dirac_surf2d.py |
| P1 | Dirac → Pauli al secondo ordine in v/c | rapporti 1,86 e 1,95 | 2 | [1,5; 2,5] | PASS | dirac_pauli.py (148 s) |
| P2 | Dirac = Pauli a ε = 0,1 | 8,9·10⁻³ σ | 0 | 2·10⁻² σ | PASS | dirac_pauli.py |
| P3 | serve il termine di spin | 0,33–0,36 σ | O(1) | > 0,1 σ | PASS | dirac_pauli.py |
| H1 | 5 livelli di Dirac numerici vs formula | ≤ 4,5·10⁻¹¹ relativo | 0 | 10⁻⁹ | PASS ×5 | dirac_hydrogen.py (253 s) |
| H2 | struttura fine 2p₃/₂ − 2p₁/₂ | 10,94369 GHz | 10,96905 GHz (misura) | 0,5 % (QED attesa ~α/π) | PASS (−0,23 %) | dirac_hydrogen.py |
| — | Lamb shift 2s₁/₂ − 2p₁/₂ | 0 | 1057,845 MHz | — | **non riprodotto** (non è in Dirac) | dirac_hydrogen.py |
| — | ionizzazione 1s | 13,598468 eV | 13,598435 eV | — | +8,15 GHz = Lamb 1s | dirac_hydrogen.py |
| H3 | 1s: F/G costante | −0,0036487249 | −0,0036487249 | 10⁻⁶ | PASS | dirac_hydrogen.py |
| H4 | 1s: |v| = α c sinθ | 1 ± 1,3·10⁻¹⁵ (in unità αc) | 1 | 10⁻⁶ | PASS | dirac_hydrogen.py |
| H5 | 1s: 2000 surfisti RK4, |v| ≤ c, cerchi | max |v| = α, Δz = 0, Δr⊥ = 5,6·10⁻¹¹ a₀ | — | — | PASS | dirac_hydrogen.py |

Totale: 34 PASS, 2 FAIL (Z1d, G2) su 36 test, più il Lamb shift dichiarato in anticipo come non riproducibile.
Supplemento senza test: `dirac_zbw_scan.py` (135 s), ampiezza ZBW 2D in funzione di σ (§3.1).

---

## 8. Che cosa non funziona, o funziona solo in parte (onestà)

1. **Z1d (FAIL).** Ampiezza della Zitterbewegung 2D 0,4752 (x, −4,96 %) e 0,4745 (y, −5,10 %) contro 0,5 con
   tolleranza 5 %: la componente y è fuori di 0,1 punti percentuali.
   Non è un errore della regola: per un pacchetto di larghezza finita l'ampiezza è ≈ ½(1 − 3/(4σ²)) = 0,477
   per σ = 4, e tende a 0,5 per σ grande (0,4954 a σ = 8, 0,4990 a σ = 16, §3.1). La tolleranza era stata
   scelta senza tenere conto di questa correzione, che in 2D è doppia che in 1D. Non è stata cambiata.
2. **G2 (FAIL).** Previsione congelata Ω = 1/(2σ²) = 0,03125; misurato 0,02721. La previsione dimenticava che il
   pacchetto si allarga: σ(t)² = σ²(1 + (t/2σ²)²) e Ω(t) = 1/(2σ(t)²). Calcolata **dopo** (non è un test),
   la media su 0 < t < 20 è 0,02793, e il misurato è il 97,4 % di questa; il 2,6 % che resta è compatibile con le
   correzioni relativistiche a σ = 4 (⟨S_z⟩ = 0,4926 invece di 0,5).
3. **Errori nel codice dei test, corretti (registri della prima corsa conservati).**
   (a) Z1e 2D confrontava ⟨r⟩(t) con una costante, ma il pacchetto proiettato **deriva** (0,0147 c: P₊ dello
   spinore (1,1) preferisce k_x > 0); corretto togliendo il trend, esattamente come già nel test 1D Z1c.
   (b) Il χ² 2D assegnava la massa di una cella tagliata da un bordo di bin tutta a un bin: falliva già a t = 0,
   prima di ogni movimento (χ² = 132/81, `lab/out/dirac_surf2d_run1.log`); bordi allineati alle celle.
   (c) Il test di Pauli v1 confrontava la posizione iniziale di Dirac (la funzione di registrazione salvava solo
   ogni N passi): corretto; inoltre la v1 (ε = 0,2 / 0,1 / 0,05) richiedeva ~40 minuti, quindi la v2 usa
   ε = 0,2 / 0,141 / 0,1 con rapporto atteso 2 (tolleranze v2 fissate prima della sua corsa; il punto ε = 0,2
   della v1 aveva dato 3,27·10⁻², quello della v2 3,21·10⁻²).
   (d) La prima corsa di `dirac_free.py` superava i 10 minuti sulla macchina condivisa: i test dei surfisti 2D
   sono stati spostati in `dirac_surf2d.py` con griglia 256² e T = 24 invece di 384² e T = 30.
4. **Niente Lamb shift, niente g − 2.** Sono effetti del vuoto quantistico (QED). La regola di Dirac da sola non
   può darli; è il limite più importante per l'idrogeno "vero" (8 GHz sull'1s, 1 GHz sul 2s, 25 MHz sulla
   struttura fine).
5. **Dimensioni.** Le simulazioni con le palline-onda sono in 1+1 e 2+1 D (spinori a 2 componenti). L'idrogeno è
   in 3+1 D vero ma risolto con le equazioni radiali (metodo standard ed esatto), non con il motore a palline su
   una griglia 3D; il moto dei surfisti nell'1s usa il campo di velocità ricavato da quella soluzione. Il motore
   3+1 D (spinori a 4 componenti) non è stato scritto.
6. **Energia negativa.** La Zitterbewegung esiste solo se l'onda ha componenti a energia negativa, che in teoria
   dei campi sono positroni (il "mare di Dirac"). Una regola a un solo elettrone non sa trattare la creazione di
   coppie; per pacchetti più stretti di ħ/mc la regola è fuori dal suo campo di validità fisico, anche se
   numericamente fa quello che dice l'equazione. La Zitterbewegung di un elettrone libero **non è mai stata
   osservata direttamente**; è stata vista in simulatori quantistici (ione intrappolato: Gerritsma et al.,
   *Nature* 463, 68 (2010)).
7. **La corrente non è unica.** Alla corrente di Dirac si può aggiungere un rotore senza cambiare nessuna
   statistica osservabile (Holland 1993, §12); le traiettorie singole (per esempio i cerchi dell'1s) dipendono da
   questa scelta, le distribuzioni dei clic no. Abbiamo usato la corrente di Dirac standard ψ†αψ, la scelta
   più semplice e invariante di Lorentz.
8. **Più elettroni.** Per due o più elettroni la regola di Bohm–Dirac richiede un tempo comune (una
   fogliazione, Dürr et al. 1999), come già nella regola a molti fotoni (`MULTIFOTONE.md` §2.3). Non è stato fatto.

---

## 9. La regola per il progetto (da copiare nel motore)

**Palline-onda** (come quelle del fotone, con una differenza): portano uno spinore (2 numeri complessi in 2D,
4 in 3D) invece di un vettore di Jones; si propagano con il propagatore di Dirac libero (spettrale, esatto); un
potenziale V (per esempio il nucleo) entra come fase e^{−iVdt} nel passo di split-operator.

**Surfista**, a ogni tick:
1. legge lo spinore ψ dove si trova (interpolazione cubica dalle palline-onda vicine);
2. calcola v = c ψ†αψ / ψ†ψ (in 2D: vx = 2Re(a*b)/|ψ|², vy = 2Im(a*b)/|ψ|²);
3. fa un passo RK4 con dt = 0,02 ħ/mc² (o il passo più grande che risolve la frequenza più alta presente:
   2mc²/ħ se c'è energia negativa, altrimenti la frequenza delle energie in gioco);
4. nasce con probabilità ψ†ψ; non c'è nessuna moneta.

**Quando usarla al posto della regola non relativistica**: quando servono spin, struttura fine o velocità vicine a c.
A bassa velocità la regola di Dirac coincide (al secondo ordine in v/c) con la regola di Pauli = regola di
Schrödinger + termine di spin ½∇×(|φ|²ŝ)/(m|φ|²); il termine di spin non cambia le statistiche di arrivo ma
cambia le traiettorie (§2.3), quindi se si vogliono traiettorie "giuste" va tenuto anche nella regola lenta.

---

## 10. Come rieseguire

Dalla cartella del progetto (ogni script < 10 minuti, seed fissi):

    timeout 590 python lab/dirac_core.py        # autotest della macchina (0,2 s)
    timeout 590 python lab/dirac_free.py        # Z1, S 1D, S3 1D, G1           (203 s)
    timeout 590 python lab/dirac_surf2d.py      # surfisti 2D, S3 2D, G2, figura (246 s)
    timeout 590 python lab/dirac_zbw_scan.py    # ampiezza ZBW 2D contro sigma  (135 s)
    timeout 590 python lab/dirac_pauli.py       # limite di Pauli                (148 s)
    timeout 590 python lab/dirac_hydrogen.py    # idrogeno di Dirac, surfista 1s (253 s)
