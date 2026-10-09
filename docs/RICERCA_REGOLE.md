# Ricerca automatica di regole: può la luce, da sola o con una "struttura", fare materia?

Laboratorio: `lab/search_*.py` (Python 3.14 + numpy/scipy; PyTorch 2.14 CPU installato via pip per i
gradienti automatici). Tutti i numeri vengono da esecuzioni reali, con semi fissati, su una macchina condivisa
(i5-10300H, al massimo 2 processi nostri). Registri e risultati: `lab/out/search_*.log|json`, figure
`lab/out/search_*.png`. Come rieseguire: §8.

Prima di leggere: questo documento continua `docs/LUCE_INTRAPPOLATA.md` (che cosa era già stato provato: vuoto
lineare, Euler–Heisenberg, Kerr, Derrick, skyrmioni). Qui la domanda è più ambiziosa: **cercare
automaticamente** una regola locale che faccia nascere dalla luce oggetti stabili, fermi, contabili come
particelle — e poi **controllare se quella regola predice qualcosa che non le abbiamo detto**.

---

## 0. In breve (il verdetto prima dei dettagli)

* **La luce da sola non forma oggetti fermi** con nessuna regola locale che lasci la luce ordinaria com'è
  (teorema per i campi statici, §3.1; per quelli oscillanti l'irraggiamento, `LUCE_INTRAPPOLATA.md`).
* **Il meccanismo di Valerio — un vuoto esterno a energia più bassa che spinge e chiude la luce — è il modello
  "a sacco"** (MIT 1974, Friedberg–Lee 1977). Lo abbiamo simulato in 2D con un campo di struttura σ e tre modi
  di accoppiarlo alla luce, partendo da un urto γγ nel vuoto normale e anche da bolle già pronte: si forma al
  più una bolla effimera (4 tempi di attraversamento dopo l'urto; 15 nel caso migliore preparato a mano);
  **nessuna** tiene la luce per 100. Funziona per i gluoni, non per la luce: §2.
* **La ricerca automatica** (CMA-ES contro ricerca casuale, fitness e vincoli pre-registrati) trova regole che
  fanno oggetti stabili con carica intera e una piccola torre di stati — ma **ritrova il modello baby-Skyrme e
  spegne la luce** (0,4% dell'energia è elettromagnetica). Con la luce accoppiata davvero (8%) gli oggetti
  esistono ma non superano tutti i vincoli nel nostro budget: §3.
* **Predizioni: nessuna riuscita.** Rapporti di massa 1,5–2,9 contro 206,8 (μ/e), 1836 (p/e), 3477 (τ/e); in
  3+1 D il monopolo con g = e pesa 137–245 volte un bosone carico di pochi MeV, escluso: §4–5.
* In una frase: **quando una regola locale fa materia stabile, la materia è fatta della "struttura", non della
  luce**; e la luce libera che osserviamo vieta proprio gli ingredienti (massa, vuoto che confina,
  auto-interazione, carica conservata) che servirebbero per chiuderla.

---

## 1. Come non barare

Una regola che viene **aggiustata** finché dà un numero (per esempio il rapporto muone/elettrone 206,77) non
spiega niente: con abbastanza manopole qualunque numero esce. Le regole del gioco che ci siamo dati, **prima**
di far girare la ricerca:

1. **Spazio delle regole finito e scritto in anticipo** (§3.1): un candidato = 8 numeri.
2. **Vincoli duri e punteggio fissati in anticipo** in un file di pre-registrazione,
   `lab/out/search_rules_prereg.json`, la cui impronta SHA-256
   (`5aafa18d93fc24fac06fca76f47530ef269b20791744239fb5e657b3b0c98ee7`) è copiata in ogni checkpoint della
   ricerca. (Onestà: il codice di valutazione è stato *messo a punto* su tre teorie note — baby-Skyrme, sigma
   O(3) con gauge, Higgs abeliano — prima di scrivere il file; nessuna valutazione della ricerca vera è stata
   fatta prima.)
3. **Base di confronto**: una ricerca casuale con lo stesso numero di valutazioni.
4. **Test di predizione dopo**, sui migliori, senza ritoccare: una sola massa misurata fissa la scala, tutto il
   resto (rapporti di massa, dimensioni) è predizione. Se qualcosa coincide, si stima la probabilità che sia
   un caso (effetto "guarda altrove", *look-elsewhere*).

---

## 2. Parte A — il meccanismo di Valerio: la luce chiusa da un vuoto che spinge

### 2.1 L'idea, detta in fisica

"Una struttura procedurale, condizionata dalle onde del fotone, che dopo un urto ad altissima energia
(γγ) forma un sistema chiuso; il vuoto esterno, avendo energia **negativa** rispetto all'interno, crea una
pressione che impedisce al fotone di riaprirsi e disperdersi."

Questa idea ha un nome preciso nella fisica delle particelle: è il **modello a sacco** (*bag model*).

* **Sacco del MIT** (Chodos, Jaffe, Johnson, Thorn, Weisskopf, *Phys. Rev. D* 9, 3471 (1974)): dentro una
  bolla di raggio R i campi sono liberi; fuori c'è un vuoto con densità di energia **più bassa** di una
  quantità B (la "costante del sacco"). Avere una bolla costa B·volume; i quanti dentro hanno energia
  x/R (più piccola è la bolla, più sono "schiacciati"). L'energia totale
  E(R) = Σ xᵢ ħc/R + (4π/3) B R³ (− Z₀ ħc/R, termine di punto zero)
  ha un **minimo**: R* = (S/4πB)^{1/4}, E* = (4/3) S/R*, con S = Σxᵢ − Z₀. Esattamente la "pressione esterna
  che tiene chiuso il fotone".
* **Solitone non topologico di Friedberg–Lee** (*Phys. Rev. D* 15, 1694 (1977); rassegna Lee & Pang,
  *Phys. Rep.* 221, 251 (1992)): il "vuoto diverso" fuori non è messo a mano, è un **campo di struttura** σ con
  due minimi (σ = 0 dentro, σ = σ_v fuori). Il sacco del MIT è il suo limite di parete sottile.
* **Dielettrico di colore / superconduttore duale** ('t Hooft, *Nucl. Phys. B* 138, 1 (1978); Mandelstam,
  *Phys. Rep.* 23, 245 (1976); Friedberg & Lee 1977): il vuoto esterno ha costante dielettrica ε → 0 per il
  campo di colore, e il campo viene espulso dentro i sacchi.

**Per i gluoni (e i quark) funziona.** Per la luce c'è un problema che diciamo subito e poi misuriamo: nel
**nostro** vuoto la luce si propaga liberamente (è il vincolo C1). Quindi la fase "che confina" non può essere
il vuoto di tutto l'universo; potrebbe esistere solo **localmente**, come un guscio creato dall'urto. È
esattamente quello che abbiamo provato.

### 2.2 Il sacco del MIT fatto di luce (`search_bag.py mit`)

Radici di Bessel calcolate (condizioni del sacco n·E = 0, n×B = 0 per i campi vettoriali):

| modo | x (radice) |
|---|---|
| fotone/gluone TE (1⁺), il più basso | 2,7437 |
| TE con l = 2 | 3,8702 |
| TM (1⁻) | 4,4934 |
| quark senza massa (MIT) | 2,0428 |

**Controllo sugli adroni** (dove il sacco funziona): fissando B con la massa del nucleone (3 quark nel modo più
basso), senza termine di punto zero: B^{1/4} = 96 MeV, R_N = 1,72 fm; con Z₀ = 1,84 (DeGrand, Jaffe, Johnson,
Kiskis, *Phys. Rev. D* 12, 2060 (1975)): B^{1/4} = 125,5 MeV, R_N = 1,20 fm. La **glueball** di due gluoni TE
viene allora 864 MeV (Z₀ = 0) o 832 MeV (Z₀ = 1,84): predizione ~0,85 GeV contro 1,73 GeV della QCD su reticolo
per lo 0⁺⁺ (Morningstar & Peardon 1999): **fuori del 50%**. Il rapporto Δ/N viene **1** senza il termine
"colore-magnetico" (lo scambio di un gluone tra i quark): nel modello del MIT Δ/N = 1,31 si ottiene solo
aggiungendo quel termine con una costante α_c **fittata**. Quindi nemmeno lì Δ/N è una predizione pura; ma il
quadro complessivo (masse degli adroni leggeri entro ~10–20% con 4 parametri) è una descrizione riconosciuta.

**Sacco di fotoni.** Un fotone nel modo TE, scala fissata da una massa:

| massa che fissa la scala | B^{1/4} | raggio R* | nota |
|---|---|---|---|
| elettrone (0,511 MeV) | 0,0955 MeV | 1 413 fm = 1,4·10⁻¹² m | l'elettrone è puntiforme fino a ~10⁻²⁰ m: **escluso di 10⁸** |
| elettrone, con Casimir di Boyer | 0,0943 MeV | 1 436 fm | il Casimir sposta dell'1,7% |
| protone (938,3 MeV) | 175 MeV | 0,77 fm | dimensione giusta (r_p = 0,84 fm), ma vedi sotto |

Tre ostacoli non numerici, ma decisivi:

1. **Spin**: un numero qualsiasi di fotoni (bosoni di spin 1) in qualsiasi modo ha momento angolare **intero**.
   Un sacco di fotoni non può essere un elettrone, un muone o un protone (spin ½). Questo vale per *qualunque*
   regola che metta insieme solo fotoni senza topologia.
2. **Carica**: i fotoni sono neutri; un sacco di fotoni è neutro.
3. **Il vuoto esterno**: il sacco del MIT richiede che **fuori** la luce non possa esistere (condizione al bordo
   o fase che confina). Nel nostro vuoto la luce esiste e viaggia: §2.3–2.4 misurano che cosa succede se la fase
   che confina è solo un guscio locale.

**Il termine di punto zero (Casimir).** Nei fit del MIT si usa −Z₀/R con Z₀ > 0 (attrattivo, aiuta a tenere
piccolo il sacco). Ma il calcolo vero per una sfera conduttrice perfetta (Boyer, *Phys. Rev.* 174, 1764
(1968); Milton, DeRaad, Schwinger, *Ann. Phys.* 115, 388 (1978)) dà E_C = **+0,0462 ħc/R**: **repulsivo**
(la sfera vuole allargarsi), cioè Z₀ = −0,046, l'opposto del segno che servirebbe; per le condizioni del
sacco (duali) il valore è lo stesso. Effetto numerico: sposta R* dell'1,7% e **non** aiuta il confinamento. Il
Z₀ = 1,84 del MIT è un parametro fenomenologico, non un Casimir calcolato.

**Predizioni del sacco di fotoni e "guarda altrove".** Rapporti di massa predetti (E ∝ S^{3/4}):
2 fotoni/1 = 1,682; 3/1 = 2,280; TM/TE = 1,448; TE(l=2)/TE(l=1) = 1,294; (TE+TM)/TE = 2,070. Confronto: il più
vicino a Δ/N = 1,312 è TE(l=2)/TE = 1,294 (**−1,4%**!). È una coincidenza? Probabilità che 5 rapporti messi a
caso (log-uniformi tra 1 e 4000) ne abbiano almeno uno entro il 5% di uno dei 4 bersagli: **21%**. Non è
significativo. E se il numero N di fotoni è libero, E_N/E_1 = N^{3/4} "indovina" tutto: N = 1223 dà
206,81 (μ/e, scarto 0,02%), N = 22 484 dà 1836,14 (p/e, 0,0008%), N = 52 680 dà 3477,24 (τ/e) — perché a quegli
N i rapporti possibili distano tra loro meno dello 0,06%. **Esempio da manuale di numerologia**: una famiglia
abbastanza fitta contiene ogni numero.

### 2.3 Simulazione 2D: urto γγ nel vuoto normale con una struttura σ "condizionata"
(`search_bag.py run col_*`)

**Il modello** (2+1 D, modo TE: il campo è ψ = A_z, con E_z = −∂ₜψ, B = ∇ψ×ẑ):

L = ½ a(σ) ψₜ² − ½ b(σ) |∇ψ|² − ½ m²(σ) ψ² + ½ σₜ² − ½ |∇σ|² − U(σ),
U(σ) = ½ m_σ² σ²(1−σ)² + ΔU·s(σ),  s(σ) = 3σ² − 2σ³ (passa da 0 a 1 tra i due minimi).

σ = 0 è il vuoto normale (a = b = 1, m = 0: **Maxwell esatto**, C1 rispettato). σ = 1 è la fase "sacco", con
energia ΔU = +0,05 **sopra** il vuoto normale (il vuoto esterno ha energia più bassa: la richiesta di Valerio).
m_σ = 2 (spessore della parete ~ λ/2, tensione τ = m_σ/6 = 0,333). Tre modi in cui la struttura "sente" la luce:

| sigla | accoppiamento | che cosa fa nella fase σ = 1 | Lorentz-invariante? |
|---|---|---|---|
| D | dielettrico: a = b = ε(σ), ε(1) = ε₁ (20 o 400) | luce con la stessa velocità c ma impedenza diversa: le pareti **riflettono** | sì (accoppia σ a F = (E²−B²)/2: un'onda piana singola ha F = 0 e **non** tocca σ, solo luce che si incrocia sì) |
| U | a = 1 + g s(σ), b = 1 − g s(σ), g = 0,6 (= termine + g·u·s(σ) con u = (E²+B²)/2) | luce più lenta (indice 2): riflessione totale interna, "galleria dei sussurri" | **no** (u dipende dal sistema di riferimento) |
| M | massa del fotone m(σ)² = 64 s(σ) | per ω < 8 la luce non entra: specchio perfetto (Meissner/Higgs) | sì nel modo TE (A_z con campo carico), ma **dà massa al fotone** dove σ = 1 |

**Nota importante sul dielettrico.** Per la sola radiazione, un ε **uniforme** è un cambio di unità
(A → A/√ε): la fisica della luce libera non cambia. Conta solo il **rapporto** ε_dentro/ε_fuori. "ε piccolo
fuori" (il dielettrico di colore) e "ε grande dentro" sono quindi la stessa cosa; abbiamo messo ε = 1 fuori
(vuoto normale) e ε₁ ≫ 1 dentro.

**Urto.** Due pacchetti di luce (λ = 1, lunghi ~1λ, larghi ~2λ) partono da x = ±6 uno contro l'altro, con
σ = 0 ovunque (più un rumore 10⁻³ in σ: le "fluttuazioni"); griglia 160² con passo λ/6, RK4, spugna
assorbente ai bordi. Misura: energia di **luce** nel disco r < 3λ dopo che i pacchetti, se nulla li ferma, sono
già usciti (t = 30 = t_uscita + 2 tempi di attraversamento).

Soglia attesa per l'accoppiamento U: σ = 0 diventa instabile dove 6g·u > m_σ², cioè u > 1,1 (e₀ ≳ 1).
Per D la soglia è bassissima ((ε₁ − 1)·3(E² − B²) > m_σ²), ma solo dove i fasci si incrociano (F ≠ 0).

| accoppiamento | campo e₀ dei pacchetti | σ massimo raggiunto | bolla (σ medio nel disco > 0,5)? | luce rimasta nel disco a t = 30 | luce nel disco alla fine |
|---|---|---|---|---|---|
| D, ε₁ = 20 | 0,1 / 0,5 / 1,5 | 0,004 / 0,004 / 0,38 | no | 1·10⁻⁷ / 2·10⁻⁷ / 6·10⁻⁴ | < 10⁻⁴ |
| D, ε₁ = 400 | 0,1 / 0,5 / 1,5 | 0,004 / 0,09 / 0,15 | no | 3·10⁻⁶ / 7·10⁻³ / 4·10⁻² | 4·10⁻⁷ / 1,4·10⁻³ / 1,9·10⁻³ (t = 71) |
| U, g = 0,6 | 0,5 / 1,5 / 4 | 0,004 / 0,004 / **1,43** | solo per e₀ = 4: **sì, da t = 7,7 a t = 32,6** | 2·10⁻⁷ / 2·10⁻⁷ / 1,4·10⁻² | 2·10⁻⁴ (t = 75) |
| M, m₀ = 8 | 0,5 / 1,5 / 4 | 0,004 (= rumore) | no: la luce spinge σ verso 0 | 2·10⁻⁷ | < 10⁻¹¹ |

(Le corse si fermano da sole quando il disco resta sotto l'1% dell'energia iniziale per 3 tempi di
attraversamento: per questo finiscono a t ≈ 54–75 invece di t = 618.) In sintesi: **con la regola
relativistica (D) e con quella "Higgs" (M) l'urto non forma niente**; con la regola non relativistica U e un
urto abbastanza forte si forma davvero una bolla della fase σ = 1 — l'idea di Valerio "funziona" per un
istante — ma la bolla vive circa 25 unità di tempo (≈ 4 tempi di attraversamento) e trattiene ~1% della luce,
poi collassa e lascia uscire tutto.

### 2.4 Il caso più favorevole: la bolla già pronta (`search_bag.py run seed_*`)

Per separare "si forma?" da "regge?", abbiamo messo a mano una bolla già formata (σ = 1 in un disco di raggio
3, pareti morbide) con dentro luce in un modo che gira (ω = π), e misurato quanta energia **di luce** resta nel
disco r < 5 (tempo di attraversamento t_c = 6).

Frazione dell'energia **di luce** nel disco rispetto all'inizio (per M include il termine di massa ½m²ψ², per
questo può superare 1 mentre la bolla si stringe):

| configurazione | 1 t_c | 5 t_c | 10 t_c | 25 t_c | 100 t_c | che cosa succede |
|---|---|---|---|---|---|---|
| D, ε₁ = 20, bolla con luce | 0,25 | 1·10⁻³ | 9·10⁻⁵ | 2·10⁻⁶ | < 10⁻¹² | la luce esce a ogni rimbalzo, la bolla collassa |
| D, ε₁ = 400 (registro: corsa interrotta a t = 168 per tempo) | 0,13 (t = 10) | 0,10 | 0,067 | 7·10⁻⁴ (t = 148) | — | la bolla collassa in ~1,6 t_c (σ max 1 → 0,32 a t = 10), poi la luce si disperde |
| U, g = 0,6, modo "galleria dei sussurri" (ℓ = 15) | 1,06 | 1,03 | 3·10⁻⁴ | 2·10⁻⁸ | 6·10⁻¹⁵ | regge ~5 t_c, poi la bolla (vuoto falso) collassa e la luce esce |
| idem, luce × 2,25 (amp 0,45) (registro: interrotta a t = 538) | — | 0,77 (t = 29) | 0,56 (t = 58) | 4·10⁻⁴ (t = 96) | 0 | **la più longeva**: ~15 t_c, poi collassa |
| idem, luce × 4 (amp 0,6) | 0,95 | 0,19 | 0,24 | < 10⁻⁵ | 0 | troppa luce: deforma e rompe la bolla |
| M, guscio di fase massiva (vuoto normale dentro e fuori) | 1,50 | 1·10⁻³ | 4·10⁻⁴ | 5·10⁻⁶ | 10⁻¹¹ | il guscio (vuoto falso) si chiude in ~1 t_c |
| M, Friedberg–Lee vero: fase massiva **ovunque fuori**, ΔU < 0 | 1,86 | 1,73 | 0,089 | 5·10⁻³ | 1·10⁻³ | regge ~5 t_c, poi il vuoto esterno invade la bolla |
| idem, luce × 2,25 | 0,75 | 0,44 | 0,57 | 4·10⁻³ | 4·10⁻⁴ | ~10 t_c |
| idem, luce × 4 | 0,55 | 0,63 | 9·10⁻³ | 7·10⁻⁴ | 1·10⁻⁴ | ~5–10 t_c |

**Nessuna** configurazione tiene il 90% della luce per 100 tempi di attraversamento: C2 fallisce ovunque,
anche partendo dal caso più favorevole.

**Coleman (bolle senza luce).** Teoria a parete sottile in 2D: E(R) = 2πRτ ± πR²|ΔU|. Una bolla di vuoto
"falso" (più energia dentro) ha entrambi i termini che spingono a chiudersi: **collassa sempre**. Una bolla di
vuoto "vero" (meno energia dentro) si allarga senza fine se R > R_c = τ/|ΔU| = 6,67 (Coleman, *Phys. Rev. D* 15,
2929 (1977)). Verificato: bolla falsa R = 3 → collassa (σ medio nel disco 0,37 → 0,06); bolla vera R = 4 < R_c
→ collassa (0,45 → 0,01); bolla vera R = 8 > R_c → si espande e riempie tutta la scatola (energia della regione
interna da +16,5 a −21,7, σ = 1 ovunque). Quindi un guscio "locale" di una fase confinante è sempre di vuoto
falso (se no l'universo intero sarebbe già in quella fase, e la luce non si propagherebbe — contro C1) e deve
essere **tenuto aperto dalla pressione della luce**, che a sua volta deve essere **tenuta dentro dal guscio**.

Un dettaglio che si vede in tutte le corse con σ: dopo che la bolla è scomparsa resta nel disco un grumo di
energia del solo campo σ (2,2 su 19 per U, 2,3 su 11 per la bolla falsa di Coleman) che oscilla e si spegne
lentamente: è un **oscillone** (Gleiser, *Phys. Rev. D* 49, 2978 (1994)), un oggetto quasi stabile del campo di
struttura, **senza luce** e senza carica conservata.

### 2.5 Perché il sacco di luce non regge (i conti)

1. **Fuga attraverso il guscio.** Con il dielettrico, a ogni rimbalzo sulla parete la luce passa con
   probabilità 1 − R = 4ε₁/(1 + ε₁)² ≈ 4/ε₁ (incidenza normale): 18% per ε₁ = 20, 1% per ε₁ = 400. Per tenere il
   90% dell'energia per 100 tempi di attraversamento serve 1 − R ≲ 10⁻³, cioè ε₁ ≳ 4000. Con la massa (M) la
   luce passa il guscio per effetto tunnel, e il guscio deve essere spesso molte volte 1/√(m² − ω²).
2. **Il guscio non è tenuto su.** Con l'accoppiamento Lorentz-invariante ε(σ)F, la forza della luce su σ è
   ½ε'(σ)(E² − B²), e per un modo stazionario in una cavità ⟨E²⟩ = ⟨B²⟩ in media: la luce **non** spinge in
   media sulla bolla (lo abbiamo visto: la bolla dielettrica collassa in ~1,6 t_c). La bolla falsa (ΔU > 0) più
   la tensione della parete la chiudono; chiudendosi, la frequenza della luce sale e la luce scappa.
   L'accoppiamento U (alla densità di energia u) spinge davvero (la bolla regge da 5 a ~15 t_c con dentro la luce), ma
   non abbastanza e **rompe la relatività**: u non è lo stesso in tutti i sistemi di riferimento, e ogni fotone
   singolo molto energetico (per esempio un fotone della radiazione cosmica di fondo visto da un protone dei
   raggi cosmici) attiverebbe la regola. Vincoli sulla violazione di Lorentz nel settore dei fotoni sono
   strettissimi (Coleman & Glashow, *Phys. Rev. D* 59, 116008 (1999); Kostelecký & Russell, *Rev. Mod. Phys.* 83,
   11 (2011)).
3. **Niente numero intero conservato.** Il solitone di Friedberg–Lee esiste perché il campo dentro è
   **complesso** con una carica di Noether conservata Q = ∫ ω|φ|² (il numero di particelle): non può
   disperdersi senza portare via la carica. Il campo della luce è **reale**: il "numero di fotoni" non è
   conservato da un'interazione non lineare (si possono convertire due fotoni in uno di frequenza doppia, o
   cedere energia a σ). Al massimo è un invariante adiabatico. Quindi C3 (un intero conservato) **fallisce per
   costruzione** per un sacco di luce; un sacco di **gluoni** ha invece il colore (e i quark il numero
   barionico).
4. **Anche il caso "vero Friedberg–Lee"** (fotone massivo **ovunque** fuori, vuoto esterno più basso) nelle
   nostre corse regge solo 5–10 t_c: con un campo reale la luce scambia energia con σ e le pareti, la bolla si
   restringe, la frequenza sale finché la luce supera la soglia m₀ e esce (punto 3: non c'è una carica che lo
   impedisca). E comunque è una teoria in cui la luce **non** si propaga nel vuoto: violazione diretta di C1
   (limite sperimentale sulla massa del fotone m_γ < 10⁻¹⁸ eV, Particle Data Group).

**Verdetto della parte A.** Il meccanismo "vuoto esterno che spinge" è fisica vera e funziona per i **gluoni**
(sacco del MIT, Friedberg–Lee, dielettrico di colore: masse degli adroni leggeri a ~10–20%, raggio del protone
giusto), perché lì il vuoto di QCD *è* una fase che confina, creata dall'auto-interazione dei gluoni. Per la
**luce** non funziona, per una ragione che i numeri rendono esplicita: il nostro vuoto lascia passare la luce,
quindi una fase che la chiude può essere solo un guscio locale di vuoto falso, che (i) perde luce a ogni
rimbalzo, (ii) non è sostenuto in media dalla luce se la regola è relativistica, (iii) non conserva nessun
intero. Nelle nostre simulazioni l'urto γγ forma **al più** una bolla effimera (accoppiamento U, campo e₀ = 4:
la bolla compare a t ≈ 8 e sparisce a t ≈ 33, circa 4 tempi di attraversamento), mai un oggetto che regge 100.

---

## 3. Parte B — la ricerca automatica (`search_rules.py`)

### 3.1 Lo spazio delle regole

Si lavora in 2+1 D (un piano più il tempo: costa poco), c = 1. Campi:

* **luce**: il potenziale A_μ, con la libertà di gauge U(1) (la regola deve dare la stessa fisica per A e
  A + ∇χ: è ciò che rende il fotone senza massa); campo magnetico B (in 2D è un numero), campo elettrico E;
* **struttura** (la "procedura generativa che aiuta le palline a prendere forme"): un vettore unitario
  n = (n₁, n₂, n₃), |n| = 1, in ogni punto. La luce "vede" la struttura in due modi: ruota la coppia
  (n₁, n₂) (cioè ψ = n₁ + i n₂ ha carica g) e può cambiare la costante dielettrica del vuoto.

Energia statica (per una configurazione ferma, A₀ = 0):

E = ∫ d²x [ ½|D n|² + V(n) + (κ²/2)|D_x n × D_y n|² + ½(1 + ζ(1 − n₃))B² + a₄B⁴ + a₆B⁶ ],
V(n) = μ₁²(1 − n₃) + μ₂²(1 − n₃)² + ν n₃²,   D = ∇ − i g A sulla parte carica ψ,

e la dinamica completa è quella della lagrangiana relativistica corrispondente (termine cinetico ½|Dₜn|²,
parte temporale del termine di Skyrme, elettrodinamica L_EM = −ε(n)𝓕 − 4a₄𝓕² − 8a₆𝓕³ con 𝓕 = (B² − E²)/2).

**Un candidato = 8 numeri** θ = (g, ζ, μ₁, μ₂, ν, κ, a₄, a₆) nella scatola
g ∈ [0, 3], ζ ∈ [−0,45, 3], μ₁, μ₂ ∈ [0, 1,5], ν ∈ [−1, 1], κ ∈ [0, 1,5], a₄ ∈ [−0,3, 0,3], a₆ ∈ [0, 0,3].

| pezzo | che cosa rappresenta | teoria nota se è l'unico acceso |
|---|---|---|
| a₄, a₆ (polinomio negli invarianti della luce) | auto-interazione della luce, tipo Euler–Heisenberg / Born–Infeld troncati | elettrodinamica non lineare |
| g | accoppiamento luce ↔ struttura (gauge) | modello sigma O(3) con gauge (Schroers, *Phys. Lett. B* 356, 291 (1995)) |
| ζ | la struttura cambia il dielettrico per la luce | "dielettrico di colore" (Friedberg–Lee) |
| μ₁, μ₂, ν | potenziale della struttura: vuoto al polo (n₃ = 1) o all'equatore | baby-Skyrme (Piette, Schroers, Zakrzewski, *Z. Phys. C* 65, 165 (1995)); con ν > 0 vuoto all'equatore = Higgs abeliano nel limite di London (vortici di Nielsen–Olesen, *Nucl. Phys. B* 61, 45 (1973)) |
| κ | termine di Skyrme (quartico) | baby-Skyrme; con g: Gladikowski, Piette, Schroers, *Phys. Rev. D* 53, 844 (1996) |

**Escluso dallo spazio, con il motivo (analitico):**

* **Solo luce, con qualunque non linearità** (qualsiasi L(𝓕, 𝒢) con energia convessa): *teorema* — non
  esiste nessun campo statico, localizzato, senza sorgenti, diverso da zero. Dimostrazione in due righe: statico
  e senza correnti ⇒ ∇×H = 0 ⇒ H = ∇ψ; allora ∫H·B = ∫∇ψ·B = −∫ψ ∇·B = 0; ma con energia convessa H·B > 0 dove
  B ≠ 0, quindi B = 0 ovunque (lo stesso per E, con ∇×E = 0 e ∇·D = 0). Le uniche "palle di luce" possibili
  oscillano, e oscillando irradiano perché la luce non ha soglia di massa (`LUCE_INTRAPPOLATA.md` §5.3). Per
  questo nello spazio la non linearità della luce è sempre **accompagnata** dalla struttura n.
* **Termine di Chern–Simons** (esiste solo in 2+1 D): dà massa al fotone ⇒ viola C1.
* **Simmetria di gauge non abeliana senza struttura** (Yang–Mills puro, i "gluoni" classici): teorema di
  Deser (*Phys. Lett. B* 64, 463 (1976)) e di Coleman ("There are no classical glueballs", *Commun. Math.
  Phys.* 55, 113 (1977)): nessuna soluzione statica di energia finita. Le glueball sono un effetto
  **quantistico**. Con una struttura nella rappresentazione aggiunta si ottiene il monopolo di 't Hooft–Polyakov:
  lo trattiamo in 3+1 D nella parte C.
* **Invariante 𝒢 = E·B**: in 2+1 D non esiste come scalare; torna in 3+1 D (§4.2).

### 3.2 Vincoli duri e fitness (fissati prima)

Testo completo in `lab/out/search_rules_prereg.json` (SHA-256 `5aafa18d…c98ee7`). Riassunto:

| vincolo | come si misura | perché |
|---|---|---|
| **C0** buona posizione | energia della luce convessa in B (se a₄ < 0: 144a₄² < 120a₆) | se no le equazioni non hanno soluzioni sensate |
| **C1** luce libera invariata | nel vuoto la parte quadratica in A è **esattamente** Maxwell: massa del fotone g·ρ_vuoto < 10⁻⁶ (ρ = modulo di ψ nel vuoto); niente vuoti degeneri con ρ > 0 | propagazione, dispersione, interferenza della luce ordinaria identiche; la scala E* sotto cui i termini extra sono < 1% è riportata per ogni candidato e confrontata con PVLAS/ATLAS in §4.2 |
| **C2a** oggetto fermo | rilassamento (flusso di Newton arrestato, 1500 passi) stazionario (variazione relativa dell'energia negli ultimi 300 passi < 5·10⁻³, forza massima < 0,1) e ≥ 90% dell'energia nel disco di raggio L/4 = 3,6 | "fermo e localizzato" |
| **C2b** robusto | perturbazione casuale del 5% di tutti i campi, nuovo rilassamento: stessa energia entro l'1% e stessa carica | un minimo vero (teorema di Lagrange–Dirichlet: un minimo isolato dell'energia è stabile nella dinamica); la dinamica vera per 100 tempi di attraversamento è fatta sui finalisti (§3.6) |
| **C3** intero conservato | grado topologico, in valore assoluto, = 1 e conservato (vuoto al polo: grado di n, formula esatta di Berg–Lüscher sul reticolo; vuoto all'equatore: avvolgimento al bordo) e oggetto risolto (raggio della densità topologica > 2 passi) | "contabile come una particella" |

**Punteggi morbidi** (ognuno tra 0 e 1), pesi fissati: F = [C0…C3 tutti veri] × (0,2·S_a + 0,2·S_b + 0,3·S_c + 0,3·S_d).

* S_a (Derrick): tanh(E''(1)/E)·exp(−|E'(1)|/(0,1E)), dove E(s) è l'energia dell'oggetto ingrandito di un
  fattore s; ogni termine scala con una sua potenza (gradienti s⁰, potenziale s², Skyrme e Maxwell s⁻², B⁴ s⁻⁶,
  B⁶ s⁻¹⁰). Premia un minimo netto rispetto alla scala.
* S_b (spin): 1 se la carica è del tipo π₂(S²) (vuoto al polo). In 2+1 D questo ammette un termine di Hopf
  (Wilczek & Zee, *Phys. Rev. Lett.* 51, 2250 (1983)) che rende lo spin del solitone **qualsiasi**, anche ½
  (anioni); 0 per i vortici di π₁(S¹). È un'affermazione sulla topologia, non un calcolo numerico.
* S_c (spettro): metà del numero di settori Q = 2, 3 che danno a loro volta un oggetto fermo, localizzato, con
  carica conservata e **legato** (E_Q < Q·E₁·1,005, entro la tolleranza numerica).
* S_d (luce): min(1, 4 f_EM), con f_EM = frazione dell'energia dell'oggetto Q = 1 che è campo elettromagnetico:
  misura "quanto l'oggetto è fatto di luce".

Per ordinare i candidati nella CMA-ES (anche quando F = 0) si usa rango = (numero di vincoli superati) + F.
Budget: CMA-ES 10 generazioni × 8 = 80 valutazioni (seme 2026), ricerca casuale uniforme 80 valutazioni
(seme 2027).

### 3.3 Numerica e controlli

* Griglia 48 × 48, passo h = 0,3 (scatola L = 14,4), n fissato al vuoto sul bordo, A libero. Energia sul
  reticolo con derivate covarianti sui collegamenti (gauge esatta sul reticolo: ψ del vicino viene trasportato
  con la fase e^{−igAh}); gradienti con la differenziazione automatica di PyTorch; tutte le configurazioni di
  una generazione rilassate insieme in un unico lotto.
* **Controllo (baby-skyrmione noto, κ = 1, V = 0,25(1 − n₃))**: sul reticolo E = 22,892, carica −1,0000 (il
  segno è l'orientazione), raggio 1,68; il profilo radiale 1D molto fine dà 23,182: errore di griglia −1,25%
  (`search_report.py babycheck`). Gli errori sui rapporti E₂/E₁ sono dello stesso ordine.
* **Due trappole numeriche trovate e corrette (istruttive).**
  (i) Con un campo di struttura a modulo libero (|φ| ≠ 1, potenziale λ(|φ|² − 1)², λ = 10) il solitone si
  srotola facendo passare |φ| per zero nel centro: la carica topologica **non** è protetta se il bersaglio non
  è davvero una sfera. Per questo la struttura è un vettore unitario.
  (ii) Il termine di Skyrme scritto con le sole differenze "avanti-avanti" ha un modo spurio che abbassa
  l'energia rendendo quasi antiparalleli i vicini: il solitone resta 400 passi a E = 22,35 e poi crolla a zero.
  Con la media sulle 4 combinazioni di lati della plaquette il crollo sparisce (E = 22,892 stabile per 2500
  passi, forza residua 2·10⁻⁸).
* La dinamica vera (§3.6) usa l'hamiltoniana completa nella gauge temporale A₀ = 0 (energia cinetica della
  struttura con la parte temporale del termine di Skyrme, campo D coniugato ad A), RK4; la legge di Gauss è
  soddisfatta all'inizio (tutto fermo) e conservata dalla simmetria di gauge; spugna assorbente vicino al bordo.

### 3.4 Storia della ricerca: CMA-ES contro ricerca casuale

Figura: `lab/out/search_history_bag.png` (pannello di sinistra).

| | CMA-ES (seme 2026) | ricerca casuale (seme 2027) |
|---|---|---|
| valutazioni (budget uguale) | 72 (9 generazioni × 8) | 72 (le prime 72 delle 80 fatte) |
| candidati che passano tutti i vincoli duri | **28** | **2** |
| primo candidato fattibile | valutazione 18 | valutazione 1 (fortuna) |
| fitness migliore | **0,654** | 0,467 |
| fitness mediana dei fattibili | 0,621 | 0,345 |
| bocciati da C1 (fotone con massa) | 8 | 28 |
| bocciati da C3 (carica persa) | 23 | 63 |
| g massimo tra i fattibili | 0,21 | 0,09 |
| f_EM massimo tra i fattibili | **0,0074** | 0,0007 |

Nota onesta sul budget: la pre-registrazione prevedeva 10 generazioni. La decima (8 candidati su 8 fattibili,
32 rilassamenti) richiede più dei 590 s concessi a ogni esecuzione ed è stata interrotta due volte; ci siamo
fermati a 9 generazioni e confrontiamo a budget uguale (72 contro 72). La fitness migliore della CMA-ES era già
0,654 alla valutazione 24 e non è più cambiata.

**Che cosa ha imparato la ricerca, in parole.** La CMA-ES ha portato in fretta la struttura al "vuoto al polo"
(C1) e il termine di Skyrme al massimo permesso (κ → 1,5), e ha **spento la luce**: g → 0,07, f_EM → 0,4%. Le
regole con luce davvero accoppiata (g ≳ 0,3) cadono in due modi, visibili nei registri:
(i) con g ≳ 1 e κ piccolo il solitone si restringe e la carica si perde (Q → 0: 20 casi);
(ii) con g ≈ 0,4–1,4 e κ grande l'oggetto c'è, conserva Q = 1, contiene il 4–9% di energia elettromagnetica, ma
**non** si ferma nei 1500 passi di rilassamento (energia che scende ancora dello 0,5–8%, forza 0,06–0,5): il
campo magnetico si sistema lentamente. Il vincolo C2a pre-registrato (con un budget finito) li ha esclusi.
È un limite del nostro test, e lo diciamo: per questo §3.5 ne riprende uno con un rilassamento più lungo.

### 3.5 Le regole migliori e le teorie note ritrovate

Rivalutati con 4000 passi di rilassamento (`search_report.py finalists`, figura `lab/out/search_objects.png`):

**Finalista 1 — la migliore regola trovata (CMA-ES, F = 0,654 in ricerca).**
θ = (g = 0,071, ζ = 0,081, μ₁ = 0,224, μ₂ = 0,346, ν = −0,108, κ = 1,5, a₄ = −0,062, a₆ = 0,115).
Vuoto al polo, fotone senza massa (C1 ✓); oggetto Q = 1: E₁ = 28,19, raggio 1,98, carica −1,0000, 95% dell'energia
nel disco; perturbato del 5% torna a 28,186 (−0,02%). E₂/E₁ = 1,878 (legato), E₃/E₁ = 2,834 (legato, ma con
4000 passi l'anello Q = 3 non è ancora fermo: drift 1,1%, quindi nella rivalutazione S_c scende a 0,5 e F a
0,50). Frazione di luce **f_EM = 0,38%**: è un **baby-skyrmione** (Piette, Schroers, Zakrzewski 1995) con un
filo di campo magnetico attaccato. Tutti i 28 fattibili della CMA-ES e i 2 della ricerca casuale sono di questo
tipo (g ≤ 0,21).

**Finalista 2 — la migliore regola con luce vera (scelta POST-HOC, dichiarata).** Nessun candidato fattibile ha
f_EM ≥ 5%; prendiamo quello con f_EM ≥ 5% che conserva la carica e supera più vincoli (CMA-ES, valutazione della
generazione 3): θ = (g = 0,759, ζ = −0,146, μ₁ = 0,043, μ₂ = 0,147, ν = −0,130, κ = 1,49, a₄ = 0,006,
a₆ = 0,174). Con 4000 passi: E₁ = 18,65, Q = −1,0000, **f_EM = 8,3%**, flusso magnetico 6,35 (non quantizzato:
vuoto al polo, nessun superconduttore), E₂/E₁ = 1,53, E₃/E₁ = 1,99 (molto più legati del caso senza luce). Ma
la rivalutazione dà ancora drift 0,23% e, dopo la perturbazione del 5%, un'energia più alta del 2% non ancora
ferma: **C2b non superato** anche con un budget 2,7 volte maggiore. È il **baby-skyrmione con gauge** di
Gladikowski, Piette e Schroers (1996). La sua stabilità vera è giudicata dalla dinamica (§3.6).

**Migliore della ricerca casuale**: F = 0,467, g = 0,024, κ = 1,06: ancora un baby-skyrmione senza luce.

### 3.6 Dinamica vera dei finalisti (C2 completo)

Hamiltoniana completa, RK4 (passo 0,033), oggetto Q = 1 rilassato e lasciato **fermo**, 100 tempi di
attraversamento (t_c = 2 r_E: 3,96 e 4,20), disco fisso di raggio 3,6 attorno al baricentro iniziale
(`search_rules.py dyn`, risultati `lab/out/search_rules_dyn_*.json`).

| finalista | perturbazione | energia nel disco (min / fine, rispetto all'inizio) | energia totale H (inizio → fine) | carica Q | C2 dinamico |
|---|---|---|---|---|---|
| 1 (baby-Skyrme, luce 0,4%) | nessuna | 0,991 / 1,014 | 28,192 → 28,150 (−0,15%, spugna) | −1 nel 92,5% dei campionamenti, salti isolati di ±2 | **sì** (≥ 90% per 100 t_c) |
| 2 (con luce 8%, post-hoc) | nessuna | 0,997 / 1,002 | 18,710 → 18,706 | −1 solo nel 42% dei campionamenti | energia sì, **carica no** |
| 1 | 5% | 0,12 / 0,56 | 64,7 → 29,9 | −1 nel **100%** | **no** secondo il criterio |
| 2 | 5% | — | 2,2·10⁵ all'inizio → esplode (NaN a t = 2) | — | **no** (inconcludente) |

Come leggerlo, onestamente:

* Senza perturbazione il finalista 1 è un oggetto fermo che resta fermo: energia nel disco costante entro
  l'1,4% per 100 tempi di attraversamento, energia totale conservata entro lo 0,15%. I salti di Q di ±2 in
  campionamenti isolati (Q torna subito a −1) sono un difetto della misura: nella dinamica il modulo |φ| non è
  esattamente 1 e la formula di Berg–Lüscher sbaglia dove due vicini sono quasi opposti.
* Il finalista 2 (con luce) conserva l'energia ma la sua carica misurata salta continuamente, e la figura
  `search_objects.png` mostra struttura a scacchiera nel nucleo: il nucleo **non è risolto** dalla griglia.
  Non possiamo dire che sia un oggetto stabile (né che non lo sia) senza una griglia più fine.
* La "perturbazione del 5%" pre-registrata (rumore moltiplicativo su φ, rumore additivo su A pari al 5% del
  massimo di A) nella dinamica vera inietta **più energia dell'oggetto stesso** (finalista 1: da 28 a 65, per il
  rumore su |φ| e sul campo magnetico a scala di griglia; finalista 2: il rumore su B moltiplicato dal termine
  a₆B⁶ dà 2·10⁵ e la simulazione esplode). È un difetto della nostra ricetta di perturbazione, non della regola:
  il test perturbato dinamico è quindi **inconcludente**. Nota a favore del finalista 1: anche così la carica
  resta −1 in tutti i 201 campionamenti e l'energia totale torna verso quella dell'oggetto (29,9 contro 28,2)
  mentre la spugna assorbe il rumore; l'energia nel disco fisso scende perché l'oggetto, spinto dal rumore, si
  sposta (il disco non lo segue).
* La stabilità del baby-Skyrme è comunque nota in letteratura (Piette, Schroers, Zakrzewski 1995); quella del
  baby-Skyrme con gauge anche (Gladikowski, Piette, Schroers 1996). Il nostro risultato non aggiunge né toglie
  a questo: conferma il primo, non riesce a verificare il secondo con questa griglia.

---

## 4. Parte C — promozione a 3+1 D e il vincolo C1 con numeri veri

### 4.1 Che cosa diventa in 3+1 D la regola "luce + struttura con vuoto al polo"

In 2+1 D la carica topologica della struttura è il grado della mappa piano → sfera, π₂(S²). In 3+1 D gli
oggetti puntiformi con lo stesso meccanismo (la luce resta senza massa nel vuoto, gli altri bosoni di gauge
diventano massivi, una carica topologica intera all'infinito) sono i **monopoli di 't Hooft–Polyakov**
('t Hooft, *Nucl. Phys. B* 79, 276 (1974); Polyakov, *JETP Lett.* 20, 194 (1974)) del modello di
Georgi–Glashow: gauge SU(2) (la luce più due bosoni carichi W± che **interagiscono tra loro** come i gluoni) e
una struttura nella rappresentazione aggiunta. Il versante "nodo" (carica di Hopf π₃(S²), termine di Skyrme in
3D) è il modello di Faddeev–Niemi (*Nature* 387, 58 (1997)), con energie che crescono come Q^{3/4}
(Vakulenko–Kapitanski; Battye & Sutcliffe, *Phys. Rev. Lett.* 81, 4798 (1998)): non l'abbiamo ricalcolato.

Abbiamo risolto (`search_monopole.py`) l'equazione radiale del monopolo (ansatz a riccio, 900 punti,
L-BFGS con gradiente esatto): M = (4πv/g)·f, m_W = g v, m_H = √(2λ) v.

| m_H/m_W | f calcolato (+ coda coulombiana oltre il taglio) | raggio che contiene metà energia |
|---|---|---|
| 0 (limite BPS) | 1,0008 (esatto: 1) | 2,39/m_W |
| 1 | 1,239 | 1,45/m_W |
| 2 | 1,352 | 1,22/m_W |
| 4,47 | 1,496 | 1,03/m_W |
| 14,1 | 1,659 | 0,89/m_W |

(limite m_H/m_W → ∞: f = 1,787, Kirkman & Zachos, *Phys. Rev. D* 24, 999 (1981).)

**Che cosa predice, senza ritocchi.** Se la luce della regola è la luce vera, l'accoppiamento è quello
elettrico: g = e, e²/4π = α = 1/137,036. Allora M/m_W = 4πf/e² = f/α ∈ [137, 245]: **il monopolo pesa 137–245
volte il bosone carico**, e ha **carica magnetica** 4π/e = 68,5 e (non elettrica).

* Scala fissata dal protone: m_W = 4,1–6,8 MeV, raggio ~40–70 fm (50–80 volte il protone). Un bosone vettore
  di carica e e massa ~5 MeV contribuirebbe al g − 2 dell'elettrone per ~(α/π)²(m_e/m_W)²/45 ≈ 7·10⁻¹⁰ (stima
  d'ordine di grandezza), 5000 volte l'incertezza sperimentale (1,3·10⁻¹³), e sarebbe prodotto in coppia in
  qualunque collisione e⁺e⁻ sopra ~10 MeV. **Escluso.**
* Scala fissata dall'elettrone: m_W = 2–4 keV, raggio ~10⁻¹⁰ m. **Escluso** ancora più duramente.
* Protoni ed elettroni hanno carica **elettrica**, non magnetica; le ricerche di monopoli (MoEDAL, IceCube, …)
  non ne hanno mai trovato uno.
* **Una trappola di numerologia**: per m_H/m_W ≈ 4,6 il rapporto M/m_W vale 206,8 = m_μ/m_e. Ma m_H/m_W è libero
  e M/m_W percorre con continuità [137, 245]: "ritrovare" 206,77 sarebbe un **fit** a un parametro, non una
  predizione (e il muone ha carica elettrica e spin ½).
* Spin: il monopolo da solo è bosonico; legato a un campo carico può avere spin semintero ("spin
  dall'isospin": Jackiw & Rebbi, *Phys. Rev. Lett.* 36, 1116 (1976); Hasenfratz & 't Hooft, *ibid.* 36, 1119
  (1976)). Spettro: nel limite BPS i multi-monopoli pesano esattamente Q·M (nessun legame); per m_H > 0 monopoli
  dello stesso segno si respingono: **nessuna torre** di stati legati di carica più alta.

### 4.2 Il vincolo C1 con i numeri veri della luce

* **Birifrangenza del vuoto (PVLAS)**: Δn/B² = (+12 ± 17)·10⁻²³ T⁻² misurato (Ejlli et al., *Phys. Rep.* 871,
  1 (2020)) contro 4,0·10⁻²⁴ T⁻² predetto dalla QED. Un termine extra c₁(E² − B²)² dà Δn = −4c₁B²: a 2σ
  |c₁| < 87 volte il coefficiente di Euler–Heisenberg della QED. Tradotto: **qualunque** non linearità della
  luce permessa da PVLAS è sotto l'1% di Maxwell per campi E < E* ≈ 2·10¹⁸ V/m (`search_monopole.py`), oltre
  1000 volte i laser più intensi (~10¹⁵ V/m a 10²³ W/cm²).
* **Diffusione luce-luce (ATLAS)**: osservata nelle collisioni ultraperiferiche Pb–Pb (*Nature Phys.* 13, 852
  (2017); *Phys. Rev. Lett.* 123, 052001 (2019), 8,2σ), compatibile con il Modello Standard; CMS, *Phys. Lett.
  B* 797, 134826 (2019). Non lascia spazio a un'auto-interazione della luce di ordine 1 a energie del GeV.
* **Massa del fotone**: m_γ < 10⁻¹⁸ eV (Particle Data Group). Esclude ogni regola con vuoto "di Higgs" per la
  luce (vuoto all'equatore con g > 0): nella ricerca sono scartate da C1.
* Le regole che passano C1 hanno vuoto al polo: al livello classico **nessuna** interazione luce-luce (A compare
  solo moltiplicato da coppie di quanti carichi della struttura); a un anello i quanti carichi di massa m_q danno
  un termine di Euler–Heisenberg "scalare" ∝ α²/m_q⁴ (coefficiente di (E²−B²)² pari a 7/16 di quello
  dell'elettrone a parità di massa; Dunne, hep-th/0406216), trascurabile se m_q ≫ m_e — ma allora devono
  esistere **nuove particelle cariche**, che i collisori escludono fino a ~100 GeV.

---

## 5. Test di predizione

Fatto **dopo** la ricerca, sulle regole finali, senza toccare nessun parametro (`search_report.py predict`,
`lab/out/search_predictions.json`). Numeri adimensionali che la regola produce da sola:

| regola | E₂/E₁ | E₃/E₁ | raggio × massa del quanto carico r_E·m_q | E₁/4π (sopra il limite topologico) |
|---|---|---|---|---|
| finalista 1 (migliore, senza luce) | 1,878 | 2,834 | 1,024 | 2,243 |
| finalista 2 (con luce, post-hoc) | 1,532 | 1,991 | 1,080 | 1,484 |
| baby-Skyrme di controllo | 1,877 | 2,868 | 0,862 | 1,822 |
| baby-Skyrme con gauge (GPS) di controllo | 1,618 | 2,170 | 0,769 | 1,398 |

Bersagli: m_μ/m_e = 206,768, m_τ/m_e = 3477,2, m_p/m_e = 1836,15, m_Δ/m_N = 1,312, e r_p·m_π± = 0,595
(raggio del protone 0,8409 fm × massa del pione carico, l'analogo misurato di r_E·m_q).

* **Rapporti di massa**: tutti tra 1,5 e 2,9, perché E_Q ≈ Q·E₁ meno un'energia di legame. Una torre di oggetti
  topologici con energia ~Q non può mai dare 206 o 1836 con Q = 1, 2, 3: **nessun accordo** con μ/e, τ/e, p/e
  (scarti > 99%). Il più vicino a Δ/N = 1,312 è E₂/E₁ = 1,53 (+17%) — ma Δ e N hanno lo **stesso** numero
  barionico, mentre qui Q cambia: il confronto non ha nemmeno senso fisico.
* **Dimensioni**: r_E·m_q = 0,77–1,08 contro 0,595 per il protone (+29% … +81%); e un oggetto 2+1 D non è un
  protone 3+1 D.
* **Look-elsewhere**: 16 numeri predetti × 5 bersagli. Se fossero messi a caso (log-uniformi tra 1 e 4000), la
  probabilità di almeno una coincidenza entro l'1% sarebbe 17%, entro il 5% 61%, entro il 10% 84%. Abbiamo
  **zero** coincidenze entro il 5% e una entro il 10% (E₁/4π = 1,398 per il GPS contro Δ/N, +6,5%, un numero che
  non è nemmeno un rapporto di masse). **Risultato: nessuna predizione riuscita**, e le coincidenze sono meno
  di quelle attese dal caso.
* **Spin, carica, g**: gli oggetti trovati sono neutri (nessuna carica elettrica: la luce li attraversa con un
  flusso magnetico non quantizzato); lo spin ½ è *permesso* dalla topologia (termine di Hopf) ma è una scelta,
  non una predizione; un momento magnetico g non è definito senza carica. In 3+1 D (§4) il monopolo ha carica
  **magnetica** 68,5 e e massa 137–245 volte quella del bosone carico: escluso.

Confronto utile: lo skyrmione 3D (ricalcolato in `LUCE_INTRAPPOLATA.md` §5.4) con **due** numeri fissati
predice μ_p/μ_n al 2% e le altre grandezze al 30–50%. È il livello di un modello serio; nessuna delle regole
"di luce" qui lo raggiunge.

---

## 6. Teorie note ritrovate (riassunto)

| che cosa è uscito | teoria nota | riferimento |
|---|---|---|
| il vincitore della ricerca (luce spenta, κ massimo, vuoto al polo) | modello **baby-Skyrme** | Piette, Schroers, Zakrzewski, *Z. Phys. C* 65, 165 (1995); Leese, Peyrard, Zakrzewski, *Nonlinearity* 3, 773 (1990) |
| oggetti con luce accoppiata e carica conservata | **baby-Skyrme con gauge** (flusso non quantizzato) | Gladikowski, Piette, Schroers, *Phys. Rev. D* 53, 844 (1996) |
| sigma O(3) con gauge senza termine di Skyrme (controllo): collassa in questo punto dei parametri | modello di **Schroers** (solitoni di Bogomolny solo con potenziale e accoppiamento "giusti") | Schroers, *Phys. Lett. B* 356, 291 (1995) |
| vuoto all'equatore con g > 0 (bocciato da C1) | **Higgs abeliano**, vortici di Nielsen–Olesen (il fotone prende massa: superconduttore) | Nielsen & Olesen, *Nucl. Phys. B* 61, 45 (1973) |
| promozione 3+1 D della regola vincente con luce | **monopolo di 't Hooft–Polyakov** (Georgi–Glashow) | 't Hooft 1974; Polyakov 1974 |
| bolla di vuoto falso tenuta su dalla luce | **sacco del MIT / solitone di Friedberg–Lee / dielettrico di colore** | Chodos et al. 1974; Friedberg & Lee 1977; Lee & Pang 1992 |
| grumo del solo campo σ dopo il collasso delle bolle | **oscillone** | Gleiser 1994 |
| solo luce con qualunque non linearità | nessun oggetto statico (teorema §3.1); in mezzi di Kerr solo fasci che corrono a c | `LUCE_INTRAPPOLATA.md` §5.3 |

Ritrovare teorie note è un esito onesto e anche rassicurante: la ricerca automatica, partendo da un menu di 8
manopole, converge sulle stesse strutture che i fisici hanno trovato a mano negli anni '70–'90 — e sulle loro
stesse limitazioni.

---

## 7. Verdetto

**Che cosa abbiamo trovato.**

1. **La luce da sola non fa materia ferma**, con nessuna regola locale ragionevole: è un teorema per i campi
   statici (§3.1), e per quelli che oscillano vale l'argomento della soglia di massa (`LUCE_INTRAPPOLATA.md`).
   La ricerca lo conferma per esclusione: nessuna regola senza struttura entra nello spazio dei candidati
   possibili.
2. **Il "sacco" di Valerio** (vuoto esterno più basso che spinge) è fisica vera — è il modello del MIT e di
   Friedberg–Lee — e funziona per i **gluoni** (fase del vuoto che confina, creata dalla loro auto-interazione).
   Per la **luce** l'abbiamo provato in tutte le varianti (dielettrico relativistico, accoppiamento alla
   densità di energia, fotone massivo nel guscio, Friedberg–Lee con fotone massivo fuori): l'urto γγ forma al
   più una bolla effimera (~4 tempi di attraversamento), la bolla migliore preparata a mano regge ~15, mai 100.
   Le ragioni sono quantitative: perdita a ogni rimbalzo ≈ 4ε_fuori/ε_dentro, nessuna pressione media della
   luce con una regola relativistica, Coleman (un guscio locale di vuoto confinante è sempre vuoto falso),
   nessun intero conservato per un campo reale. E il caso che funzionerebbe (fase confinante ovunque) spegne la
   luce libera: viola C1.
3. **La ricerca automatica** (CMA-ES, 72 valutazioni, contro 72 casuali; vincitrice 28 a 2 nei candidati
   fattibili, fitness 0,654 contro 0,467) ha trovato regole che fanno oggetti fermi, stabili (il migliore anche nella
   dinamica vera senza perturbazione: §3.6; il test perturbato è inconcludente), con una carica topologica intera e una piccola torre Q = 1, 2, 3 legata. **Ma ha
   ritrovato il baby-Skyrme e ha spento la luce** (f_EM = 0,4%): l'oggetto è fatto del campo di struttura, non
   di luce. Le varianti con luce accoppiata (baby-Skyrme con gauge, f_EM ≈ 8%) esistono ma, nel nostro budget
   numerico, si assestano lentamente e non superano tutti i vincoli. In 3+1 D l'analogo con luce è il monopolo
   di 't Hooft–Polyakov.
4. **Predizioni**: nessuna. I rapporti di massa delle torri topologiche sono ~Q (1,5–2,9), lontanissimi da
   μ/e = 206,8, τ/e = 3477, p/e = 1836; le dimensioni non tornano; le poche "coincidenze" sono meno di quelle
   attese dal caso. Il monopolo, con la sola ipotesi g = e, pesa 137–245 volte il bosone carico e porta carica
   magnetica: escluso dagli esperimenti. Due esempi espliciti di numerologia (N^{3/4} dei sacchi di fotoni,
   M/m_W del monopolo) mostrano come si "indovinerebbe" 206,77 barando.

**Verdetto semplice.** Le "regole che cerchiamo" per far diventare la luce materia, nella forma "luce chiusa da
una struttura locale", **non le abbiamo trovate**, e abbiamo trovato perché: la luce libera che vediamo ogni
giorno (C1) toglie esattamente gli ingredienti che servono (massa, vuoto che confina, auto-interazione, carica
conservata). Quando una regola produce materia stabile, la materia è fatta della **struttura** (pioni nello
skyrmione, campo n qui, gluoni nel sacco), non della luce. La strada che funziona in natura è un'altra: la
luce si **trasforma** in coppie di particelle che esistono come campi propri (γγ → e⁺e⁻, Breit–Wheeler), e i
campi che si auto-confinano sono quelli non abeliani (gluoni → protone). Il passo serio successivo, per il
simulatore a palline, è quindi il campo di **gluoni** su reticolo (confinamento, legge dell'area) e lo
skyrmione 3D completo, non un'altra regola per la luce.

---

## 8. Come rieseguire

```
python -m pip install torch --index-url https://download.pytorch.org/whl/cpu     # una volta
# Parte A (sacchi)
python lab/search_bag.py mit                                   # istantaneo
python lab/search_bag.py list                                  # elenco delle configurazioni
timeout 590 python lab/search_bag.py run col_U0.6_e4           # ~2-8 min ciascuna (griglia 160^2, RK4)
python lab/search_report.py bag                                # tabella riassuntiva
# Parte B (ricerca)
python lab/search_rules.py prereg                              # pre-registro (non lo riscrive se esiste)
python lab/search_rules.py known                               # teorie note dentro lo spazio (~3 min)
timeout 590 python lab/search_rules.py cma                     # ripetere fino a 9-10 generazioni (riprende)
timeout 590 python lab/search_rules.py rand                    # idem per la ricerca casuale
timeout 590 python lab/search_report.py finalists              # ~9 min: finalisti con 4000 passi + figura
timeout 590 python lab/search_rules.py dyn lab/out/search_finalist_assoluto.json assoluto 0     # ripetere: riprende
python lab/search_report.py predict                            # test di predizione + look-elsewhere
python lab/search_report.py babycheck                          # controllo del reticolo (profilo radiale)
python lab/search_report.py figs                               # figure
# Parte C
timeout 590 python lab/search_monopole.py                      # ~2 min: monopolo + numeri PVLAS
```

File: `lab/search_bag.py`, `lab/search_rules.py`, `lab/search_monopole.py`, `lab/search_report.py`;
risultati in `lab/out/search_*.json`, registri `lab/out/search_*.log`, figure `lab/out/search_objects.png`,
`lab/out/search_history_bag.png`.
