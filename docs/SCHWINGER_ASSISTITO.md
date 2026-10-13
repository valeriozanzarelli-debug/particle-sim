# Schwinger assistito: coppie dal vuoto sotto soglia con due campi insieme

Studio dei "problemi aperti". Prefisso dei file: `assist_`. Codice in `lab/assist_*.py`, risultati in
`lab/out/assist_*`. Pre-registrazione: `lab/out/assist_prereg.json` (generata da `lab/assist_prereg.py`).

**In breve.** La regola zig-zag 1+1 ritrova Schwinger (0,1%), il passaggio tunnel → multifotone di Keldysh
(0,02–0,16%, ordine a² pulito) e le regole multifotone (Δₙ ∝ ξⁿ, pendenze 1,005 / 2,980 / 4,958). Con due campi
ritrova la **soglia netta** di Schützhold–Gies–Dunne: γ_crit = **1,606** contro 1,592 previsto (+0,9%), e l'esponente
dell'istantone entro lo 0,8–1,8% per γ ≤ 6. Fuori dalle formule (vicino alla soglia, ω ~ m, ε finito) dà numeri
esatti del modello: lì il WKB a ε finito sbaglia di fattori da 2 a 10⁴, e lo sviluppo "ordine 0 + 1 in ε" sbaglia di
un fattore 2–3 vicino alla soglia (a ε = 0,1). **È validazione + esame di formule note, non una scoperta.** Due esami
pre-registrati falliscono per motivi numerici (spiegati sotto). Per i laser veri la strada resta lontanissima:
esponenti ≥ 30 anche con fotoni da 100 keV.

---

## 1. Il problema, in parole semplici

Il vuoto della regola zig-zag è il "mare" pieno di elettroni a energia negativa (la banda bassa). Un campo
elettrico molto forte può strappare un elettrone dal mare: nasce una coppia elettrone + positrone. È
l'**effetto Schwinger**. Il guaio: serve un campo enorme, E_S = m²c³/(eħ) ≈ 1,3·10¹⁸ V/m. La probabilità va
come exp(−π E_S/E): con E = E_S/10 è già e^(−31) ≈ 3·10⁻¹⁴, con i laser di oggi (E ≈ 10⁻³ E_S) è e^(−3000),
cioè zero. **Non è mai stato visto con un campo di laser.**

L'idea (Schützhold, Gies, Dunne 2008): sommare al campo **forte e lento** un campo **debole e veloce**. Il campo
veloce da solo non basta (i suoi "fotoni" hanno energia ω ≪ 2m), quello forte da solo nemmeno; insieme
l'esponente crolla. Immagine: il campo forte inclina il mare e crea una "barriera" da attraversare in
tunnel; il campo veloce dà un calcio di energia a metà strada e accorcia il tunnel.

## 2. Cosa è noto (letteratura)

Fonti lette online in questa sessione (arXiv/ar5iv), salvo dove scritto "a memoria".

* **Schützhold, Gies, Dunne, PRL 101, 130404 (2008)** (arXiv:0807.0754). Campo
  E(t) = E/cosh²(Ωt) + ε/cosh²(ωt), con ε ≪ E e Ω ≪ ω ≪ m. Parametro di Keldysh **combinato**
  γ = mω/(eE) (frequenza del campo veloce, intensità di quello forte). Con il metodo delle traiettorie
  istantone (worldline instanton) trovano una **soglia netta**: per γ < γ_crit = π/2 l'esponente resta
  quello di Schwinger; per γ > π/2 diventa
  𝒜 = (m²/eE)·[2 arcsin(π/2γ) + (π/2γ²)·√(4γ² − π²)]  (eq. 15),
  **indipendente da ε** (ε entra solo nel prefattore). Per γ grande 𝒜 → 2πm/ω.
* **Dunne, Gies, Schützhold, PRD 80, 111301 (2009)** "Catalysis of Schwinger vacuum pair production" (a memoria):
  stessa idea, istantoni per ε finito.
* **Orthaber, Hebenstreit, Alkofer, PLB 698, 80 (2011)** (arXiv:1102.2182): numerica con l'equazione cinetica
  quantistica (spinori, campo omogeneo) per due impulsi di Sauter, ℰ₁ = 0,25 E_S, τ₁ = 10⁻⁴ eV⁻¹ (≈ 51/m),
  ℰ₂ = 0,025 E_S, τ₂ = 7·10⁻⁶ eV⁻¹ (≈ 3,6/m). Trovano un aumento del numero di coppie e spettri in impulso
  strutturati. Nel testo **non ci sono numeri** (solo figure): non posso confrontare cifra per cifra.
  Nota: con la loro definizione γ = m/(eℰ₁τ₂) ≈ 1,1 < π/2, cioè sono **sotto** la soglia degli istantoni.
* **Linder, Schneider, Sicking, Szpak, Schützhold, PRD 92, 085009 (2015)** (arXiv:1505.05685): **1+1 dimensioni**,
  E(t) = E₁ + E₂·g(ω₂t) con g = 1/cosh² (Sauter), gaussiana, cos. Con WKB e istantoni:
  - Sauter: γ_crit = π/2 indipendente da ε;
  - gaussiana: γ_crit ~ √|ln ε|; coseno: γ_crit ~ |ln ε| (la soglia dipende dalla forma dell'impulso!);
  - coseno al primo ordine in ε (eq. 37): esponente = (πE_S/E₁)·[1 − 2ε cos(pγ) I₁(γ)/γ].
  Nessuna numerica esatta in quell'articolo.
* **Schneider, Schützhold, JHEP (2016)** (arXiv:1603.00864) "Prefactor in the dynamically assisted Sauter-Schwinger
  effect": prefattore con gli istantoni (QED scalare 1+1), confronto con l'integrazione numerica di Riccati
  "eccellente" a E = E_S/5, si stacca a E_S/3.
* **Torgrimsson, Schneider, Oertel, Schützhold, JHEP 06, 043 (2017)** (arXiv:1703.09203): sviluppo in potenze del
  campo debole (ε): per l'impulso di Sauter **ordine zero + primo ordine** bastano e concordano con la numerica;
  per gaussiana/coseno gli ordini alti possono dominare. Numeri della loro numerica: non estratti (sez. 6 non
  leggibile dal testo).
* **Varianti e verifiche** (a memoria, non verificate online): Akal, Villalba-Chávez, Müller PRD 90, 113004 (2014)
  (due frequenze); Hebenstreit, Fillion-Gourdeau PLB 739 (2014) (ottimizzazione della forma); Otto, Seipt, Blaschke,
  Kämpfer, Smolyansky (2015) (assistenza dinamica con soglia multifotone).
* **Identità utile** (a memoria, standard: Schmidt et al. 1998, Hebenstreit tesi 2011): per un campo omogeneo
  lineare l'equazione cinetica della QED a 3+1 dimensioni con spinori è **la stessa** dell'equazione di Dirac
  1+1 con m → m_⊥ = √(m² + p_⊥²), più un fattore 2 di spin. Quindi la 1+1 non è un giocattolo: è la fetta
  p_⊥ = 0 della risposta vera. Cambiano i prefattori (integrale su p_⊥), non gli esponenti.
* **Riferimenti per la validazione:** Schwinger/Nikishov in 1+1: numero medio di coppie per lunghezza e tempo
  (eE/2π)·exp(−πm²/eE). Impulso singolo di Sauter E sech²(t/τ): spettro esatto (Narozhny–Nikishov 1970),
  f(p) = [cosh(2πeEτ²) − cosh(πτ(ω₊ − ω₋))] / [2 sinh(πτω₊) sinh(πτω₋)], ω± = √(m² + (p ± eEτ)²); copre tutto il
  passaggio tunnel (γ = m/(eEτ) ≪ 1) → multifotone (γ ≫ 1). Campo oscillante: risonanze a n fotoni
  nω = 2q̄ (massa "vestita" dal campo), a p = 0 solo n dispari per la simmetria di coniugazione di carica
  (Mocken, Ruf, Müller, Keitel PRA 81, 022122 (2010), a memoria).

**Cosa è aperto.** (1) L'effetto non è mai stato osservato. (2) Le formule semiclassiche valgono per ω ≪ m e
ε piccolo ma "non troppo": vicino alla soglia e per ω ~ m servono numeri esatti, che esistono solo per alcuni
parametri e in forma di figure. (3) Quanto la forma dell'impulso conti (soglia netta solo per Sauter).

**Cosa possiamo fare noi, onestamente.** In un campo **omogeneo** la regola zig-zag 1+1 dà il numero di coppie
**esatto del modello**, modo per modo, come già in R1 (Schwinger 0,98–1,04, Sauter 0,992–0,996). Qui non c'è
fisica nuova da scoprire nella QED continua: la QED 1+1 in campo esterno è un problema risolto in linea di
principio (equazione cinetica). Quello che possiamo dare:
* **validazione**: la regola ritrova Schwinger, il passaggio di Keldysh e le regole multifotone;
* **esame delle formule semiclassiche** (soglia γ_crit, esponente) con numeri esatti, convergenti in a → 0;
* **numeri dove la semiclassica non vale** (vicino alla soglia, ω ~ m): previsioni esatte del modello 1+1,
  che però **coincidono** con quelle della QED 1+1 continua (la regola converge a Dirac). Non sono "scoperte":
  sono numeri controllati, confrontabili con la numerica di letteratura quando i parametri coincidono.

---

## 3. Previsione prima dei dati (fase 2)

Pre-registrazione: `lab/out/assist_prereg.json` (script `lab/assist_prereg.py`), committata **prima** delle corse di
produzione. Parametri liberi: **0**. Le prove di sviluppo fatte prima (parametri diversi da quelli di produzione:
E1 = 0,125 con ε = 0,05; impulsi singoli (0,15; 22), (0,05; 2), (0,3; 3); Floquet a ξ = 0,05–0,3; fit sul solo
Dirac continuo a ε = 0,01) sono dichiarate nello script. Mi hanno insegnato due cose che hanno cambiato gli esami:

1. **Il WKB "esatto" a ε finito non è affidabile** per l'impulso di Sauter: il punto di svolta sta attaccato al polo
   dell'impulso veloce e lì il passaggio è brusco, non adiabatico (parametro locale c = ε m²/(γ² E1) < 1 ai campi
   raggiungibili). Quindi la soglia si esamina sulla **risposta a ε → 0**,
   Q = [f(+ε) + f(−ε) − 2 f(0)]/(2ε²) (la parte "un fotone del campo veloce + tunnel nel campo forte"), il cui
   esponente è proprio l'azione dell'istantone di Schützhold–Gies–Dunne.
2. Le separazioni di Floquet sotto ~10⁻⁸ sono al pavimento numerico → ξ = 0,1 e 0,2.

**Previsioni (dalla teoria, calcolate prima):**
* soglia: con il campo lento sech²(Ωt), γ_s = mΩ/E1 = 0,2, la soglia degli istantoni sale da π/2 = 1,571 a
  **γ_crit = 1,5915**; esponente × E1: 3,111 (sotto soglia) → 2,775 (γ = 2) → 1,994 (γ = 3) → 1,035 (γ = 6) →
  0,626 (γ = 10). Quindi a E1 = 0,1 E_S l'aumento previsto a γ = 6 è e^(20,8) ≈ 10⁹ in esponente (prefattore a parte).
* Schwinger (impulso lento γ_s = 0,05): N per lunghezza, regola locale (LCFA) 2,234·10⁻¹⁴ contro esatto 2,283·10⁻¹⁴
  a E = 0,1 (differenza 2%: entro la soglia del 10%).
* Coseno debole: esponente WKB e formula di Linder (primo ordine) coincidono entro l'1% dove la correzione
  2εI₁(γ)/γ ≤ 0,05; si separano oltre (γ ≳ |ln ε|).

| esame | cosa | soglia |
|---|---|---|
| V1a | spettro Schwinger (impulso lento) regola vs esatto | 3% (θ = 0,0125), Richardson 1% |
| V1b | N regola vs tasso locale di Schwinger | 10% |
| V2 | passaggio di Keldysh γ = 0,1…10, f(0) e N vs formula esatta | 3%, Richardson 1% |
| V3a | regola di selezione: a p = 0 niente risonanze pari | Δ_pari/Δ_dispari ≤ 10⁻² |
| V3b | legge multifotone Δ_n ∝ ξⁿ | pendenza entro 5% di n |
| V3c | posizione risonanze regola vs Dirac | 10⁻³ |
| A1 | assistito Sauter: regola vs Dirac continuo | 3% (θ = 0,025), Richardson 1%; Q entro 5% |
| A2 | esponente di Q vs istantone, γ ≤ 6 | 3% |
| A3 | soglia γ_crit dal fit | 5% di 1,5915 |
| O1 | assistito coseno: regola vs Dirac | 3%, Richardson 1% |
| O2 | esponente vs WKB (dove 2εI₁/γ ≤ 0,3) | 3% |
| O3 | esponente vs Linder eq. 37 (dove 2εI₁/γ ≤ 0,05) | 3% |
| A4, O4 | numeri vicino a soglia e a ω ~ m: **previsione** del modello | solo riportati |

---

## 4. Cosa abbiamo fatto (fase 3)

**Il metodo è quello di R1.** In un campo elettrico omogeneo E(t) (dipende dal tempo, non dal posto) ogni onda piana
del mare resta un'onda piana: il mare si evolve **modo per modo**, una matrice 2×2 per ogni impulso p. Il campo entra
come fase di Peierls sul salto (corri a c), poi la rotazione di angolo θ = m a (la massa). All'inizio il modo è nella
banda bassa (vuoto), alla fine si misura quanto è finito nella banda alta: è il numero di coppie f(p) in quel modo.
Il numero di coppie per lunghezza è N = ∫ dp/2π f(p). I doppioni del reticolo (buco a k = π, vedi R1b) non si
toccano mai: ogni modo usato resta lontano da k = π (controllo nel codice).

**Tre modi di fare lo stesso conto** (`lab/assist_core.py`):
1. la **regola** zig-zag a quattro passi reticolari a = 0,2 / 0,1 / 0,05 / 0,025 (in unità 1/m), più
   l'estrapolazione di Richardson a → 0 (l'errore va come a², verificato: il rapporto degli errori fra passi
   successivi è 3,97–4,00);
2. il **Dirac continuo** 1+1 (quello che risolve l'equazione cinetica quantistica), integrato con Magnus al 4º ordine;
   indipendente dalla regola; controllato contro la formula esatta dell'impulso di Sauter a 10⁻⁸;
3. la **semiclassica**: WKB/Dykhne col punto di svolta complesso (= istantone per campi che dipendono dal tempo), e le
   formule chiuse di letteratura (Schützhold–Gies–Dunne eq. 15; Linder et al. eq. 37).

**I campi** (dati come potenziale A(t), così E si spegne liscio):
* forte e lento: E1 sech²(Ωt), con γ_s = mΩ/E1 = 0,2 (praticamente costante sul tempo del tunnel, ~m/E1);
* debole e veloce, tipo Sauter: ε E1 sech²(ωt) (come nell'articolo originale);
* debole e veloce, oscillante: ε E1 cos(ωt) con la stessa busta lenta.
Parametro di Keldysh combinato: **γ = mω/(eE1)**.

## 5. Risultati

### 5.1 Validazione (ritrovare cose note)

**V1 – Schwinger** (impulso lento γ_s = 0,05, `lab/out/assist_val.json`):

| E/E_S | f(p) regola/esatto, peggior modo (a = 0,0125) | Richardson | N regola / N esatto | N regola / Schwinger locale |
|---|---|---|---|---|
| 0,25 | 0,05% | 4·10⁻⁷ | 0,9997 | 1,010 |
| 0,16 | 0,07% | 8·10⁻⁷ | 0,9995 | 1,014 |
| 0,10 | 0,10% | 5·10⁻⁶ | 0,9992 | 1,021 |

PASS (soglie 3%, 1%, 10%). La differenza dell'1–2% con la formula di Schwinger "locale" è fisica (l'impulso ha
γ_s = 0,05), non un errore: con la formula esatta dell'impulso l'accordo è 0,1%.

**V2 – Keldysh, un impulso solo** E sech²(t/τ), da tunnel a multifotone:

| γ = m/(eEτ) | E/E_S | f(0) esatto | regola/esatto a = 0,025 | N regola/esatto | Richardson |
|---|---|---|---|---|---|
| 0,1 | 0,2 | 1,57·10⁻⁷ | 0,9984 | 0,9983 | 1 − 4·10⁻⁶ |
| 0,3 | 0,2 | 2,11·10⁻⁷ | 0,9985 | 0,9984 | 1 − 4·10⁻⁶ |
| 1 | 0,2 | 2,23·10⁻⁶ | 0,9991 | 0,9991 | 1 − 10⁻⁶ |
| 3 | 0,1 | 2,77·10⁻⁷ | 0,9996 | 0,9997 | 1 − 2·10⁻⁷ |
| 10 | 0,05 | 5,89·10⁻⁶ | 0,99998 | 1,0001 | 1 − 4·10⁻⁹ |

PASS. Nota: a γ = 10 Schwinger darebbe 5·10⁻²⁸, il valore vero è 6·10⁻⁶: è il regime multifotone, e la regola lo
segue senza saperne niente.

**V3 – multifotone, campo oscillante** A = ξ m sin ωt, p = 0, risonanze dai livelli di Floquet:
* V3b PASS: separazione alla risonanza a n fotoni Δₙ ∝ ξⁿ, pendenze 1,005 / 2,980 / 4,958 (n = 1, 3, 5; Dirac uguale
  alla quarta cifra).
* V3c PASS: posizioni delle risonanze regola/Dirac entro 2·10⁻⁶. La formula adiabatica della massa vestita
  nω = 2q̄ le dà entro lo 0,33% (non giudicato).
* **V3a FAIL (tecnico).** La regola di selezione (a p = 0 niente risonanze pari) era da verificare come
  Δ_pari/Δ_dispari ≤ 10⁻²: a ξ = 0,1, n = 4 viene 0,030. Il motivo è il minimizzatore (si ferma a Δ ~ 10⁻⁸, il
  pavimento numerico), non la fisica. **Diagnosi a posteriori dichiarata** (`lab/assist_v3diag.py`): stringendo la
  griglia in ω il minimo pari scende da 10⁻⁸ a 10⁻¹² (incrocio vero), quello dispari resta fermo a 3,8·10⁻⁷; il
  rapporto vero è < 3·10⁻⁶. Il verdetto pre-registrato resta FAIL.

### 5.2 Il punto aperto: due campi insieme (impulso di Sauter veloce)

**A2 – esponente e soglia.** Per ogni γ ho misurato la risposta Q (la parte proporzionale a ε²: un "calcio" del
campo veloce + tunnel nel campo forte) a cinque campi forti E1 = 0,16 … 0,0625 E_S, e ho estratto l'esponente col fit
ln Q = −S/E1 + ν ln(1/E1) + c. Previsione: azione dell'istantone (ε → 0).

| γ | 0,5 | 1,0 | 1,5 | 1,6 | 1,7 | 1,8 | 2,0 | 2,5 | 3 | 4 | 6 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S regola (× E1) | 3,128 | 3,123 | 3,125 | 3,108 | 3,071 | 3,008 | 2,823 | 2,340 | 1,988 | 1,522 | 1,027 | 0,603 |
| S istantone | 3,111 | 3,111 | 3,111 | 3,109 | 3,051 | 2,965 | 2,775 | 2,335 | 1,994 | 1,529 | 1,035 | 0,626 |
| rapporto | 1,005 | 1,004 | 1,004 | 1,000 | 1,007 | 1,014 | 1,018 | 1,002 | 0,997 | 0,995 | 0,992 | **0,964** |

**A2 PASS** (tutti entro 3% per γ ≤ 6). A γ = 10 la frequenza veloce è ω = 0,6–1,6 m: la semiclassica (che vuole
ω ≪ m) comincia a sbagliare, del 3,6%. Il Dirac continuo dà gli stessi S alla terza cifra.

**A3 PASS – la soglia.** Fit a un solo parametro della curva degli istantoni sui punti della regola:
**γ_crit = 1,606**, previsto 1,5915 (= π/2 corretto per la busta lenta γ_s = 0,2), differenza **+0,9%**.
Sotto soglia l'esponente non si muove (rapporto 1,000–1,007); sopra crolla, come detto da Schützhold, Gies e Dunne.

**A1 FAIL (numerico, dichiarato).** Regola contro Dirac continuo, f(0) in 390 punti: il caso peggiore è 3,8% a
a = 0,025 e 21% dopo Richardson (soglie 3% e 1%). **Tutti** i punti cattivi sono a ε < 0 vicino alla soglia, dove
l'ampiezza di Schwinger e quella del calcio veloce si **cancellano** (f crolla fino a 10⁻²⁵, vicino al pavimento
dell'arrotondamento): lì un piccolo errore in ampiezza diventa un grande errore relativo. Sui 260 punti con ε ≥ 0 il
peggiore è 0,8% (a = 0,025) e 0,3% (Richardson); Q, che è quello che entra in A2–A3, va entro lo 0,6%. Senza i
44 punti con f < 10⁻²⁰ (controllo a posteriori, dichiarato) il Richardson va entro lo 0,15%.

**A4 – numeri esatti del modello dove le formule non bastano** (E1 = 0,1 E_S; `lab/out/assist_analysis.json`):

| γ | 1,0 | 1,5 | 1,8 | 2,0 | 2,5 | 3 | 4 | 6 | 10 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| aumento di f(0), ε = 0,01 | 1,5 | 2,0 | 3,6 | 6,7 | 62 | 577 | 1,4·10⁴ | 3,5·10⁵ | 2,6·10⁶ | – |
| aumento di f(0), ε = 0,1 | 30 | 97 | 313 | 792 | 9,0·10³ | 7,6·10⁴ | 1,6·10⁶ | 3,6·10⁷ | 2,6·10⁸ | – |
| N insieme / (N forte + N debole), ε = 0,01 | 1,33 | 1,31 | 1,41 | 1,75 | **5,8** | 3,6 | 1,7 | 1,19 | 1,05 | 1,01 |
| N insieme / (N forte + N debole), ε = 0,1 | 17 | 30 | 73 | **114** | 18 | 5,1 | 2,0 | 1,23 | 1,05 | 1,01 |
| f esatto / f (ordine 0 + 1 in ε), ε = 0,1 | 2,3 | 3,0 | 3,1 | 2,9 | 1,9 | 1,5 | 1,15 | 1,03 | 1,004 | – |
| N WKB (ε finito) / N regola, ε = 0,1 | 1,00 | 1,06 | 1,35 | 1,70 | 2,8 | 4,5 | 10 | 31 | 31 | 316 |

Cosa si legge:
* **L'aumento per il modo centrale è enorme** (fino a 10⁸), ma l'aumento del **numero totale** di coppie è al
  massimo ~100 (ε = 0,1, γ ≈ 2) e ~6 (ε = 0,01, γ ≈ 2,5). Perché: oltre γ ≈ 3 il campo veloce **da solo** (a
  ω = γ E1 ≳ 0,3 m) produce già più coppie del campo forte, e la "collaborazione" diventa piccola. La finestra utile è
  **subito sopra la soglia**, 1,6 < γ < 3. Questo è coerente con la letteratura, ma i numeri esatti di questa
  tabella (1+1 D, questi impulsi) non li ho trovati pubblicati.
* **Il WKB a ε finito sbaglia** di fattori 2–300 sopra soglia a ε = 0,1 (fino a 3·10⁴ a ε = 0,01): il punto di svolta sta attaccato al polo dell'impulso
  veloce e lì il passaggio non è adiabatico. Era la lezione della fase 2: l'esponente è giusto, il prefattore no.
* **Ordine 0 + 1 in ε** (Torgrimsson et al. 2017 dicono che per Sauter basta): vero lontano dalla soglia (γ ≥ 6:
  entro il 3%), ma a ε = 0,1 vicino alla soglia gli ordini alti valgono un **fattore 2–3**. A ε = 0,01 la stessa
  verifica è fatta in A2 tramite Q. Non contraddice l'articolo (che parla di ε piccolo), lo quantifica.
* A γ < γ_crit l'aumento a ε = 0,1 (×17–30) non è "assistenza dinamica": è solo che il campo totale al centro è
  E1(1 + ε) — l'esponente π/(1,1·E1) invece di π/E1.

### 5.3 Campo debole oscillante (coseno)

* **O1 PASS**: regola/Dirac entro 0,8% (a = 0,025), 0,11% (Richardson), tutti i 135 punti.
* **O2 PASS**: dove la correzione 2εI₁(γ)/γ ≤ 0,3, esponente regola / WKB fra 0,9993 e 1,0092.
* **O3 PASS**: dove la correzione ≤ 0,05, esponente relativo / formula di Linder (eq. 37) fra 1,0000 e 1,0098.
  Oltre, la formula al primo ordine si stacca come previsto (1,036 a γ = 5, 1,13 a γ = 6 con ε = 0,01).
* **O4 – oltre la soglia del coseno** (γ ≳ |ln ε|, previsione del modello): ε = 0,01, γ = 8 / 10 / 12 → aumento del
  numero totale ×236 / ×45 / ×428; ε = 0,001, γ = 10 / 12 → ×42 / ×789. Qui la dipendenza da E1 non è più un solo
  esponente (residui del fit 0,05–0,6) e il WKB col punto di svolta dominante sbaglia l'esponente del 4–26%: ci sono
  molti punti di svolta che interferiscono (Dumlu–Dunne). Non c'è una formula chiusa da confrontare: i numeri esatti
  sono nel JSON (`lab/out/assist_analysis.json`, voci specO). Il numero non sale in modo monotono (γ = 10 < γ = 8 a
  ε = 0,01) perché a ω = m il campo veloce da solo produce già molte coppie.

## 6. Validazione, scoperta, limite

* **Validazione:** Schwinger, Keldysh, regole multifotone, soglia γ_crit e esponente dell'istantone, formula di
  Linder: tutto ritrovato dalla regola, con convergenza a² pulita e 0 parametri.
* **Scoperta:** **nessuna** in senso stretto. La QED 1+1 in campo omogeneo è un problema risolto in linea di principio
  e la regola converge a Dirac; i numeri della sez. 5.2–5.3 sono esatti e controllati ma qualunque codice cinetico
  li darebbe. Al massimo è una **quantificazione** non trovata pubblicata: (i) a ε = 0,1 vicino alla soglia
  l'ordine 0+1 sbaglia ×2–3; (ii) la finestra utile per il numero totale è stretta, 1,6 < γ < 3.
* **Limiti:**
  - *1+1 contro 3+1.* Per campi omogenei la QED 3+1 con spin è la 1+1 con m → m_⊥ e un fattore 2 (a memoria,
    Schmidt et al. 1998). Fatto il conto con la regola sommando su m_⊥ (`lab/assist_limits.py`, E1 = 0,1, ε = 0,1,
    γ = 3): aumento 1+1 = 5,13, aumento 3+1 = **4,77** (regola e Dirac uguali entro 0,03%). Gli esponenti non
    cambiano, i prefattori di qualche %–decine di %. Ma questo usa la regola 1-D a massa m_⊥: la regola zig-zag 3-D
    vera avrebbe la sua dispersione trasversa (debito).
  - *Campi veri.* Un laser non è un campo elettrico omogeneo: ha campo magnetico, è focalizzato, ha onde stazionarie.
    Tutto questo qui manca (debito), come in quasi tutta la letteratura citata.
  - *Quanto lontano sono i laser.* (Esponente dell'istantone, forma di Sauter, ε → 0: la stima più ottimista.)

| campo forte | E/E_S | Schwinger da solo | + XFEL 10 keV (ω = 0,020 m) | + 25 keV | + 100 keV | + 300 keV |
|---|---|---|---|---|---|---|
| 10²³ W/cm² (classe 10 PW, oggi) | 6,6·10⁻⁴ | e^(−4789) | e^(−321) | e^(−128) | e^(−32) | e^(−11) |
| 10²⁴ W/cm² (prossima generazione) | 2,1·10⁻³ | e^(−1515) | e^(−320) | e^(−128) | e^(−32) | e^(−11) |
| 10²⁵ W/cm² (oltre i piani) | 6,6·10⁻³ | e^(−479) | e^(−306) | e^(−127) | e^(−32) | e^(−11) |

  Il guaio si vede subito: con i laser veri γ è enorme (30–900), e lì l'esponente diventa **2πm/ω**, cioè non
  dipende quasi più dal campo forte. È lo stesso esponente dell'impulso veloce **da solo** (la coda di Fourier di
  sech²). Detto semplice: con questi numeri il campo forte non aiuta più, e servono fotoni "veloci" da ~100 keV o più
  (raggi gamma, non XFEL). Per un XFEL vero (quasi monocromatico, forma più simile al coseno) non c'è soglia netta: la
  soglia γ ~ |ln ε| dipende dall'intensità del campo veloce, e con ε piccolo il guadagno è molto minore (sez. 5.3).
  Mi sembra in linea con quello che si legge in letteratura (a memoria, non verificato: servono campi forti di
  qualche 10⁻³–10⁻² E_S *e* fotoni veloci ben sopra le decine di keV); qui è rifatto con numeri nel modello 1+1.
  Anche e^(−11) non vuol dire "visto": il numero assoluto dipende da volume, durata e prefattore.

## 7. Debiti (cose messe a mano)

1. **Campo esterno dato**: E(t) omogeneo, forma sech² o coseno, scelto da me (è l'esperimento simulato, non un
   parametro della regola). Nessuna reazione delle coppie sul campo (niente schermo).
2. **1+1 dimensioni**; la stima 3+1 usa l'identità m → m_⊥ (citata a memoria) con la regola 1-D, non la regola 3-D.
3. **Busta lenta** γ_s = 0,2 per accendere/spegnere il campo forte: sposta la soglia da π/2 a 1,592 (calcolato e
   previsto, non un aggiustamento).
4. Il fit dell'esponente usa una legge con potenza (ν ln(1/E1)) per il prefattore: è una scelta di forma, dichiarata
   nella pre-registrazione.
5. Numeri di letteratura delle numeriche (Orthaber 2011, Torgrimsson 2017) **non estratti** (solo figure): il
   confronto con la "numerica di letteratura" è fatto con il mio Dirac continuo indipendente, non con le loro cifre.
6. Diagnosi a posteriori dichiarate: `assist_v3diag.py` (V3a) e `A1_posthoc` in `assist_analysis.py` (A1).
   Una correzione di codice dopo la pre-registrazione: finestra temporale del solo campo oscillante debole
   (`Field.T`), necessaria per farlo girare; nessun esame cambiato.

**Parametri liberi: 0.** Ingressi: m (tramite θ = m a), e i campi (l'esperimento).

## 8. Come rifare le corse (≈ 15 minuti in tutto, 1 core)

```
export OMP_NUM_THREADS=1
python lab/assist_prereg.py          # pre-registrazione (teoria)
python lab/assist_val.py V1; python lab/assist_val.py V2; python lab/assist_val.py V3
python lab/assist_v3diag.py          # diagnosi dichiarata di V3a
python lab/assist_run.py A; python lab/assist_run.py O
python lab/assist_run.py specA; python lab/assist_run.py specO   # spettri (JSON grandi, non committati)
python lab/assist_run.py wkb         # ~5 min
python lab/assist_analysis.py; python lab/assist_limits.py; python lab/assist_figure.py
```

Figura: `lab/out/assist_figura.png` (in alto a sinistra la regola sopra la formula esatta da tunnel a multifotone;
in alto a destra la soglia; in basso gli spettri e l'aumento del numero di coppie).

## 9. Bilancio PASS/FAIL contro la pre-registrazione

| esame | esito |
|---|---|
| V1a, V1b Schwinger | PASS, PASS |
| V2 Keldysh | PASS |
| V3a regola di selezione | **FAIL** (pavimento del minimizzatore; diagnosi a posteriori: rapporto vero < 3·10⁻⁶) |
| V3b, V3c multifotone | PASS, PASS |
| A1 regola vs Dirac (Sauter) | **FAIL** (punti di cancellazione a ε < 0; per ε ≥ 0 entro 0,8%) |
| A2 esponente istantone | PASS |
| A3 soglia γ_crit | PASS (1,606 contro 1,592) |
| O1, O2, O3 (coseno) | PASS, PASS, PASS |
