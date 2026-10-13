# Vuoto che polarizza: il mare di Dirac della zig-zag attorno a un nucleo

Studio con prefisso `vuoto_`. Codice: `lab/vuoto_core.py` (motore), `lab/vuoto_run1d.py`, `lab/vuoto_run3d.py`,
`lab/vuoto_analyze3d.py`, `lab/vuoto_levels.py`, `lab/vuoto_runB.py`, `lab/vuoto_posthocB.py`, `lab/vuoto_figure.py`.
Pre-registrazione: `lab/out/vuoto_prereg.json` (commit 9773fdb, **prima** di ogni corsa di produzione), emendamento
`lab/out/vuoto_amend.json` (commit 77eac84, prima dei numeri V3/V4), esplorazione B `lab/out/vuoto_prereg_B.json`
(commit 9d4c1db). Risultati: `lab/out/vuoto_1d.json`, `vuoto_3d_exams.json` (tabelle grezze `vuoto_3d.json` +
`vuoto_3d_b.json`), `vuoto_levels.json`, `vuoto_B.json`, `vuoto_posthocB.json`, `vuoto_posthoc_Hfine.json`.
Figura: `lab/out/vuoto_figure.png`.

---

## 0. In breve (per Valerio)

**La domanda.** Il vuoto della regola zig-zag è il "mare di Dirac": tutti gli stati a energia negativa pieni. Se
metti un nucleo, il mare si deve deformare (le coppie virtuali si orientano: gli elettroni del mare un po' verso il
nucleo, i "buchi" un po' via). È la **polarizzazione del vuoto**, il primo pezzo di QED che non ha bisogno di fotoni
quantizzati. Nella natura dà il potenziale di **Uehling** e sposta il 2S dell'idrogeno di −27 MHz.

**Cosa è uscito, in quattro righe.**
1. **Il mare intero della regola, così com'è, non polarizza bene.** In 1+1 D non polarizza **per niente** (esatto,
   a ogni ordine: il "doppione" della regola è un vuoto capovolto che cancella tutto). In 3-D la regola ha 32 copie
   pesanti normali e 32 copie **senza massa capovolte** (non le conoscevamo): il mare intero dà una polarizzazione
   finita ma **sbagliata** (anti-schermo che cresce a grande distanza). Previsto prima e misurato (−0,3 %).
2. **Contando solo la specie fisica** (correzione dichiarata, come R1 "conta solo il buco k = 0") la regola dà la
   QED vera: il logaritmo della rinormalizzazione della carica **2α/3π a −0,07 %**, la forma di Uehling entro
   **0,1 %** (fino a q ≈ 20 m), il potenziale di Uehling nello spazio entro **0,2 %**. Zero parametri adattati.
3. **Livelli:** Nd⁵⁹⁺ 1s −11,61 eV contro −11,67 (**−0,5 %, PASS**). Idrogeno 2S/1S: **FAIL** pre-registrato
   (+14 %): non è il vuoto, è il reticolo dell'**atomo** (un picco della funzione d'onda sul primo sito). Con reticoli
   più fini (a posteriori, non valutato) 2S = **−26,85 MHz** contro −26,89 e 1S −214,9 contro −215,1 (−0,1 %).
4. **Idea geometrica (provata):** se la massa (l'angolo di flip) "respira" fra tick pari e dispari (θ−μ, θ+μ), i
   doppioni senza massa prendono massa μ e diventano **regolatori di Pauli–Villars**: allora è il **mare intero, senza
   proiettore**, a dare Uehling (0,03–0,3 %) e una carica rinormalizzata **finita**, Π(0) = (2α/3π) ln(M/m) (0,5 %).
   Le soglie pre-registrate di B però erano per μ grandi e sono **fallite**; il successo è a posteriori.

**Quello che NON è uscito e non poteva uscire:** l'autoenergia (≈ +1085 MHz del Lamb 2S–2P, quasi tutti i 7,1 GHz
dell'1S–2S), g−2, i 25 MHz della struttura fine: servono la luce quantizzata e il loop elettrone–fotone (piano al §8).

**Validazione, non scoperta.** Uehling è del 1935: qui la regola lo **ritrova**. Le cose nuove sono sulla regola: il
teorema "il mare intero 1-D non polarizza", i 32 coni senza massa del 3-D, l'anisotropia del vuoto del reticolo, e il
trucco della moneta che respira.

---

## 1. Che cosa si calcola (spiegato semplice)

* **Il vuoto.** La regola fa un passo (tick) alla volta. Gli stati stazionari di una regola a passi non hanno
  un'energia ma una **quasi-energia** Ω (un angolo per periodo). Il **mare** = tutti gli stati con Ω fra −π e 0.
  In 1-D il periodo è 1 tick; in 3-D sono **2 tick** (la regola alterna l'ordine xyz e zyx).
* **Il nucleo.** È la stessa fase della regola, e^{−iV} dopo ogni moneta (`zz_core.Walk`). Per il vuoto bastano
  potenziali deboli: il risultato è al **primo ordine in Zα**, che è esattamente l'ordine di Uehling.
* **La risposta.** Per un potenziale V(x) = V_q e^{iqx} calcolo **esattamente** come cambia il proiettore sul mare
  (formula di Riesz per un operatore unitario, solo coppie pieno–vuoto, sicura anche con livelli degeneri) e quindi la
  densità indotta δn(q) = −χ(q) V_q. Tutto in spazio k con le matrici di Bloch di `zz_core`.
* **Controlli del motore** (`lab/vuoto_dev_check.py`, sviluppo): la formula contro differenze finite (5·10⁻¹⁰); il
  codice in spazio k contro il **mare vero** ottenuto diagonalizzando la regola (`zz_core.Walk`) su un anello 1-D di
  64 siti e su una scatola 3-D 8³ (accordo 10⁻⁸).
* **Unità**: ħ = c = mₑ = 1; sul reticolo a = 1, quindi m = θ. Per il 3-D: Π(q) = 4πα·χ(q)/q² per specie; la carica
  rinormalizzata usa lo **schema on-shell** (la carica vera è quella vista da lontano, q → 0), lo stesso ovunque:
  Π_R(q) = Π(0) − Π(q). Uehling: Π_R = (2α/π)∫x(1−x)ln(1+x(1−x)q²/m²)dx.

---

## 2. I doppioni: chi c'è nel vuoto della regola

**1-D (periodo 1 tick).** Due specie: k = 0 (massa θ, buco a quasi-energia 0) e k = π (massa θ, buco a
quasi-energia π). Il mare "Ω ∈ (−π, 0)" riempie la specie k = 0 normalmente, ma la specie k = π **al contrario**
(pieno sopra il suo buco, vuoto sotto). Esiste una simmetria esatta: con S = (−1)^x, **S U S = −U** per **qualunque**
potenziale locale. Conseguenza (teorema): la densità del mare intero è **uniforme per ogni potenziale**, a tutti gli
ordini. Il doppione capovolto cancella la specie fisica, sito per sito.

**3-D (periodo 2 tick).** W(k + π e_j) = W(k). In ogni cella di lato π ci sono:
* **4 specie con massa θ** in (0,0,0), (π/2,π/2,0), (π/2,0,π/2), (0,π/2,π/2) (stessa massa e velocità, misurate);
* **4 coni senza massa** a quasi-energia π in (π/2,0,0), (0,π/2,0), (0,0,π/2), (π/2,π/2,π/2), velocità cos θ
  (misurata 0,995 c a θ = 0,1). Il mare li riempie **capovolti**.

In k-spazio è un **reticolo di salgemma** (NaCl): un reticolo FCC di specie pesanti e un FCC di coni senza massa
negli interstizi ottaedrici. Su tutta la zona: 32 + 32. **Questi coni senza massa sono una cosa nuova sulla regola**
(il debito #3 parlava solo del raddoppio). Per la fisica a una particella con potenziali lisci non danno fastidio (non
si mescolano), ma il **vuoto** li vede tutti.

**La correzione dichiarata** (come R1): un proiettore liscio sulla specie fisica, Q = Σ_k √w(|k|) |k⟩⟨k|, con
w = 1 per |k|a ≤ K1 = 0,6 e 0 per |k|a ≥ K2 = 1,05 (gradino C^∞), **lo stesso** per ogni θ, q, 1-D e 3-D. È anche il
regolatore (il "cutoff"). Cambiandolo in (0,5; 0,9) i risultati rinormalizzati cambiano di 0,2–0,35 % (V2g).

---

## 3. Risultati contro la pre-registrazione

| esame | che cosa | misurato | soglia | esito |
|---|---|---|---|---|
| V1a | mare intero 1-D, potenziale forte V(0) = −0,8 per tick | spread della densità 4,6·10⁻¹⁵; χ in spazio k 7·10⁻¹⁶ | < 10⁻¹² | **PASS** |
| V1b | specie fisica 1-D contro Dirac continuo 1+1 D | θ = 0,025: +0,021 %, +0,020 %, −0,008 %, −0,46 % (q/m = 0,3; 1; 3; 10); errore ∝ θ² (×3,9–4,3) | ≤ 1 % e ≥ ×3 | **PASS** |
| V1c | specie fisica con periodo 1 o 2 tick | 2,4·10⁻¹⁵ | 10⁻⁹ | **PASS** |
| V2a | mare intero 3-D: log a grande distanza | pendenza 0,016916 contro 1/(6π² cos θ) = 0,016972 (−0,33 %); Π < 0 ovunque | 2 % | **PASS** |
| V2b | log della rinormalizzazione della carica | fit L = 0,016874 contro 1/(6π²) = 0,016887 (**−0,07 %**); ultimo passo −0,14 % | 2 % | **PASS** |
| V2c | forma di Uehling (θ = 0,0125) | −0,055 %, −0,057 %, −0,071 %, −0,096 % (q/m = 0,3; 1; 3; 10) | 0,5 % / 1,5 % | **PASS** |
| V2d | coefficiente q² (quello che sposta i livelli dell'idrogeno) | 0,99944 | 1 % | **PASS** |
| V2e | isotropia di Π_R (direzione (1,1,1) contro media x,y,z) | +0,073 % | 1 % | **PASS** |
| V2f | anisotropia di Π(0): costante del reticolo, non errore θ² | (z−x)/L = 0,37 → 0,41 → 0,42 → 0,43; cambia +4,4 % da 0,05 a 0,0125 | < 30 % | **PASS** (cioè: l'anisotropia **resta**) |
| V2g | indipendenza dal regolatore (K1,K2) | Π_R +0,22 % / +0,35 %; L −0,08 % | 0,5 % / 2 % | **PASS** |
| V3 | potenziale di Uehling nello spazio (nucleo smussato uguale per i due) | max scarto 0,21 % per r = 0,05…3 λ_C | 2 % / 3 % | **PASS** |
| V4a | H 2S₁/₂ | −30,62 MHz contro −26,89 (+13,9 %) | 1 % | **FAIL** |
| V4b | H 1S₁/₂ | −245,0 MHz contro −215,1 (+13,9 %) | 1 % | **FAIL** |
| V4c | H 2P₁/₂ | +0,005 MHz (rif. −0,00008) | < 0,01 MHz | **PASS** (ma è il rumore del metodo, vedi §5) |
| V4d | Nd⁵⁹⁺ 1s | **−11,605 eV** contro −11,667 (−0,53 %); a 3 punti −0,43 % | 2 % | **PASS** |
| V4e | Nd⁵⁹⁺ 1s, potenziale messo davvero nella regola contro primo ordine | rapporto 1,00058 | 1 % | **PASS** |

Totale: **15 PASS, 2 FAIL** (gli esami dell'idrogeno). Convergenza della quadratura 3-D: 5·10⁻⁷.

Esplorazione B (moneta che respira, pre-registrata a parte):

| esame | misurato | esito |
|---|---|---|
| B1 Uehling dal mare intero (θ = 0,025) | μ = 0,3: −4,4 / −3,8 / −1,7 %; μ = 0,6: −16 / −14 / −4 % | **FAIL** (soglia 2 %) |
| B2 log della carica | μ = 0,3: +2,3 %; μ = 0,6: +12,8 % | **FAIL** (soglia 3 %; μ = 0,3 da solo passerebbe) |
| B3 log della massa del regolatore | +35 %, +45 % | **FAIL** |
| B4 carica rinormalizzata finita e schermante | sì in tutti i casi | PASS |

---

## 4. Il vuoto della regola: spiegazione dei numeri

**Il logaritmo (V2b).** Più fine è il reticolo, più coppie virtuali "vede" il vuoto, e la carica nuda va schermata di
più: Π(0) = (2α/3π) ln(Λ/m) + costante, con Λ ~ 1/a. Il reticolo dà la pendenza esatta (fig. b): 0,016874 contro
0,016887. Prima di avvicinarsi al continuo la pendenza sbaglia di −1,4 % (θ 0,1→0,05), −0,45 %, −0,14 %.

**La forma (V2c, V2d, V3).** Dopo aver tolto Π(0) resta la parte fisica, uguale a Uehling entro 0,1 % da q = 0,02 m a
q ≈ 22 m; sale a +0,4 % a 32 m e +1,0 % all'ultimo punto (48 m, dove il regolatore comincia a tagliare) (fig. c), e nello spazio entro 0,2 % da 0,05 a 3 λ_C (fig. d). La coda e^{−2r} c'è.

**Il mare intero (V2a).** Senza proiettore: specie pesanti meno coni senza massa. Per q ≫ m le due cose si
cancellano (niente polarizzazione a corta distanza!), per q ≪ m resta il log dei coni con il segno sbagliato:
l'interazione di Coulomb diventerebbe **più forte a grande distanza**, come ln(r). È fisica sbagliata, e la regola lo
fa esattamente come previsto (0,3 %).

**L'anisotropia (V2f) — una cosa da sapere.** La regola non ha la simmetria del cubo: alterna xyz e zyx, e la densità
si guarda alla fine dei 2 tick. Il vuoto risulta un po' "birifrangente": la costante della carica nuda lungo z supera
quella lungo x di 0,43·(1/6π²), e non sparisce affinando il reticolo. In pratica la carica vista da un nucleo, in una
sola fetta di tempo, dipende dalla direzione di 4πα·0,0073 ≈ **7·10⁻⁴**; mediando sulle due fette di tempo (che si
scambiano x e z) resta **4·10⁻⁵** fra y e (x, z). Sugli stati s non conta (conta solo la media sulle direzioni), ma è
un difetto del vuoto della regola. **Nuovo debito.**

---

## 5. I livelli

**Metodo.** Spostamento al primo ordine: ΔE = ⟨V_U⟩ = −(2Zα/π)∫Π_R(q) ρ̃(q) dq, con Π_R del vuoto della regola
(θ_vac = 0,0125, fino a q = 48 m) e ρ̃(q) dalla densità (g² + f²) dello stato della **regola radiale**
(`zz_core.RadialWalk`). Riferimenti (solo esame): H con la densità di Schrödinger e Uehling analitico (−215,10 e
−26,887 MHz; la formula di testa dà −217,2 e −27,15; la differenza −1 % è il termine (Zα)⁵ di portata finita);
Nd⁵⁹⁺ con la densità esatta di Dirac 1s (−11,667 eV; formula di testa −11,632). Nucleo puntiforme, massa infinita.

**Nd⁵⁹⁺ (lo Z del Laboratorio Atomo): PASS.** θ_atom = 0,05 / 0,025 / 0,0125: −12,977 / −11,975 / −11,711 eV;
estrapolato con la legge a^{2γ} (γ = 0,899): **−11,605 eV** (−0,53 %). Messo dentro la regola come potenziale vero,
l'autovalore si sposta di −10,725 eV contro −10,719 al primo ordine (stesso nucleo smussato): 0,06 %.

**Idrogeno: FAIL pre-registrato, e perché.** Con θ_atom = 0,05 e 0,025 viene −356,8 e −272,9 MHz (1S), cioè la
convergenza non va come θ² e l'estrapolazione pre-registrata sbaglia di +14 %. Il colpevole l'ho trovato dopo: la
funzione d'onda della regola radiale ha sul **primo sito** (r = a/2) una densità **37 volte** (θ = 0,05) e **108
volte** (θ = 0,025) quella vera. È l'errore di contatto degli stati s (debito #5) visto da vicino: un pezzetto di
"doppione" dell'elettrone incastrato vicino al nucleo. Pesa pochissimo nella norma, ma Uehling vive proprio lì (r ≲
λ_C/2) e lo amplifica. A posteriori (**non valutato**):
* reticoli più fini (`vuoto_posthoc_Hfine.json`): 1S −233,4 (θ 0,0125), −219,5 (0,00625) → **−214,9 MHz** (−0,10 %);
  2S −29,17, −27,43 → **−26,85 MHz** (−0,13 %);
* lo stesso proiettore della specie fisica applicato anche all'elettrone (il picco sparisce): 1S −213,6 (−0,7 %), 2S
  −26,70 (−0,7 %) già con θ = 0,05 e 0,025.

**2P₁/₂.** Il metodo in spazio q ha un pavimento numerico di ~0,005 MHz (con densità e Uehling **analitici** dà
+0,0052 invece di −0,00008): il PASS di V4c dice solo "piccolo", il segno non significa niente.

---

## 6. Bilancio del Lamb (onesto)

Numeri di riferimento della QED (solo esame, arrotondati, Eides–Grotch–Shelyuto, Phys. Rep. 342, 63 (2001)):

| pezzo | 2S–2P (Lamb, 1057,8 MHz) | 1S–2S (la regola/Dirac sta +7,11 GHz sopra la misura, `NUCLEO_1S.md`) | chi lo può dare |
|---|---|---|---|
| **polarizzazione del vuoto (Uehling)** | **−27,1 MHz** | 1S −215,1, 2S −26,9 → **+188 MHz** (va nel verso opposto) | **la regola, col mare di Dirac: fatto qui** (2S −26,85, 1S −214,9 → +188,0 MHz, a posteriori sui reticoli fini; Nd⁵⁹⁺ PASS) |
| **autoenergia** (elettrone che emette e riassorbe un fotone, compreso il momento magnetico anomalo ≈ +68 MHz) | ≈ **+1085 MHz** | ≈ **−7,3 GHz** (quasi tutto) | serve la **luce quantizzata** + il loop elettrone–fotone: **non c'è** (piano al §8) |
| rinculo oltre la massa ridotta, due loop, ordini alti | qualche MHz | qualche decina di MHz | rinculo: serve il protone dinamico; due loop: dopo l'autoenergia |
| dimensione del protone (0,84 fm) | +0,14 MHz (sul 2S) | 1S +1,1, 2S +0,14 → −1,0 MHz | dal nodo del protone (studi COLLA/NUCLEO): non qui |

In parole: dei ~7,1 GHz mancanti la polarizzazione del vuoto ne spiega **+0,19 GHz con il segno contrario** (quindi
il buco da riempire con l'autoenergia è in realtà ≈ 7,3 GHz). Del Lamb 2S–2P ne spiega **−27 MHz su 1058** (2,6 %).
Il resto è autoenergia: senza fotoni virtuali la regola non può averla.

---

## 7. Idee geometriche e visive (richiesta di Valerio del 13/10)

Per ognuna: cosa si vedrebbe, perché funziona, costo (Intel UHD, ~42 ms per passo con ~33 800 palline), regole.

**7.1 La moneta che respira (provata, la più promettente).** *Cosa si vede:* il vuoto come una scacchiera nel tempo:
nei tick pari le palline girano un po' meno (θ−μ), nei dispari un po' di più (θ+μ); la massa media resta θ.
*Perché funziona:* i coni senza massa stanno a quasi-energia π, dove contano le **differenze** fra tick pari e
dispari; la differenza 2μ apre lì un buco di massa μ (misurato: 0,300 e 0,600 per μ = 0,3, 0,6) e lascia la massa
fisica intatta (0,1000). I doppioni capovolti diventano esattamente i **regolatori di Pauli–Villars** (particelle
fantasma con segno opposto e massa grande). A posteriori, con μ = 0,1 e θ = 0,0125 (M/m = 8): il **mare intero, senza
proiettore**, dà Uehling regolato alla Pauli–Villars entro −0,32 % / −0,04 % (q = 1, 3 m) e Π(0) = (1/6π²) ln(M/m)
entro **0,4 %**. Con μ grandi (0,3–0,6, quelli pre-registrati) i regolatori sono deformati dal reticolo (−3 / −15 %):
B FAIL. *Costo:* zero (un numero in più per tick). *Regole:* **cambia la regola** (un parametro μ in più, che fa da
cutoff fisico): serve la decisione di Valerio. Vantaggio: niente proiettore, niente correzione a mano, la carica
rinormalizzata finita esce da sola. Svantaggio: μ è una manopola (anche se i risultati a q ≪ μ/a non dipendono da
lei).

**7.2 La mappa delle specie: un cristallo di salgemma nello spazio degli impulsi.** *Cosa si vede:* il cubo di
Brillouin con palline luminose: blu (specie pesanti normali, FCC) e rosse (coni senza massa capovolti, negli
interstizi): il "vuoto" è un cristallo NaCl di particelle e anti-vuoti. Con la moneta che respira le rosse si
"gonfiano" (prendono massa). *Perché:* §2, misurato. *Costo:* trascurabile (64 palline). *Regole:* è solo una vista.

**7.3 La scacchiera di Feynman capovolta (1-D).** *Cosa si vede:* il mare 1-D su una scacchiera; colorando i siti
pari e dispari, il doppione è la stessa figura con bianco e nero scambiati (S = (−1)^x). Mettendo un nucleo fortissimo
la densità del mare resta **piatta** (10⁻¹⁵). È il modo più chiaro di far vedere "perché il vuoto intero non si
polarizza". *Costo:* ~10³ palline, nulla. *Regole:* sì.

**7.4 L'alone di Uehling attorno al nucleo.** *Cosa si vede:* attorno al nucleo di Nd⁵⁹⁺ un guscio morbido di palline
di carica indotta, dentro la nube 1s (raggio di Bohr 2,3 λ_C, alone a ~0,3–1 λ_C): carica che si accumula vicino al
nucleo e un velo opposto un po' più fuori (somma zero dopo la rinormalizzazione). Le palline dell'alone potrebbero
"respirare" fra le due fette di tempo (l'anisotropia x↔z del §4 è visibile: un ellissoide che si gira a ogni tick).
*Perché:* il profilo ρ_ind(r) esce dalla tabella Π_R già calcolata (una trasformata radiale). *Costo:* un profilo
radiale calcolato una volta; qualche migliaio di palline: niente per la Intel UHD. *Regole:* sì, è densità della
regola. È **rendering**: va come proposta al thread "Modello standard a palline".

**7.5 Rompere la birifrangenza con l'ordine degli assi.** *Cosa si vedrebbe:* un vuoto con carica uguale in tutte le
direzioni. *Perché dovrebbe funzionare:* l'anisotropia nasce dall'ordine xyz/zyx; ciclando sulle 6 permutazioni (o
sulle 3 rotazioni cicliche xyz, yzx, zxy con i loro palindromi) il tensore di rango 2 diventa per forza isotropo
(simmetria cubica ripristinata sul periodo). *Costo:* periodo di 6 tick invece di 2; da controllare che non nascano
nuovi coni. *Regole:* cambia la regola (l'ordine degli assi): decisione di Valerio. Non provato.

**7.6 La cipolla dell'atomo (lega con NUCLEO_1S).** Il FAIL dell'idrogeno viene dal primo sito vicino al nucleo.
Celle più fini vicino al nucleo (gusci) o lo stesso proiettore di specie anche sull'elettrone (provato a posteriori:
−0,7 %) tolgono il picco. *Cosa si vedrebbe:* la nube dell'1s senza la "spina" luminosa al centro.

---

## 8. Piano per l'autoenergia dentro la regola (con fattibilità)

Il Lamb vero è quasi tutto autoenergia: l'elettrone zig-zag emette e riassorbe una pallina di luce zig-zag. Due
strade, entrambe con la stessa regola per elettrone e luce (accoppiamento minimo e = √(4πα), come in `emiss_core.py`
dello studio EMISSIONE, modo "qe" a un'eccitazione):

1. **Il loop in spazio k (la strada di questo studio, consigliata).** Lo stesso motore di `vuoto_core` (matrici di
   Bloch, formula di Riesz) applicato al **settore a un fotone**: propagatore dell'elettrone-walk × propagatore della
   luce-walk, integrato su k in 3-D. Dà (a) lo spostamento della massa δm ∝ (3α/4π) m ln(Λ²/m²) → il θ nudo va
   corretto (un **debito dichiarato**: θ_nudo ≠ θ_fisico); (b) il fattore di vertice → **g−2 = α/2π** come esame
   pulito; (c) la parte ad alta energia dell'autoenergia legata. *Costo:* integrali 3-D come quelli fatti qui (~10 s
   per punto), qualche ora di CPU.
2. **Il logaritmo di Bethe dalla regola radiale (parte a bassa energia).** ln k₀ = media logaritmica delle energie di
   eccitazione pesata con |⟨n|p|m⟩|²: si ottiene dal risolvente della regola radiale nei canali κ = +1 e −2 integrando
   sull'energia del fotone. *Fattibilità misurata:* la fattorizzazione LU della regola radiale per H fino a 200 a₀ a
   θ = 0,1 (274 000 siti) costa **1,4 s** e ogni soluzione 0,04 s: 200 energie × 2 canali × 4 reticoli ≈ 40 minuti.
   Esame: ln k₀(1S) = 2,984, ln k₀(2S) = 2,812.

Perché non nello spazio reale come EMISSIONE: per l'autoenergia servono insieme fotoni da ~10 eV (λ ~ 10⁻⁷ m) e da
~mₑc² (λ ~ λ_C): un fattore 5·10⁴ in lunghezza, impossibile in una scatola 3-D (anche con il reticolo sparso del
RETICOLO_UNIVERSO servirebbero ~16 livelli). Il loop in k e il risolvente radiale non hanno questo muro.

---

## 9. Debiti

**Pagati (in parte):**
* #2 QED assente → la **polarizzazione del vuoto** ora emerge dal mare della regola (Uehling, log della carica,
  livelli di Nd⁵⁹⁺). Restano autoenergia, g−2, i 25 MHz.
* #5 "l'1S–2S non vede il Lamb": la parte di vuoto (+188 MHz) si calcola; l'1S–2S della regola resta dominata
  dall'autoenergia mancante.

**Nuovi o rimasti:**
1. **Proiettore di specie** (correzione dichiarata, 2 numeri K1 = 0,6, K2 = 1,05; effetto sul rinormalizzato
   ≤ 0,35 %). Alternativa senza proiettore: la moneta che respira (§7.1), che cambia la regola.
2. **La regola 3-D ha coni senza massa a quasi-energia π** (32 su tutta la zona). Il mare intero della regola dà una
   polarizzazione sbagliata (log a grande distanza). Da aggiungere al debito #3.
3. **Anisotropia del vuoto** (birifrangenza del reticolo): 7·10⁻⁴ per fetta di tempo, 4·10⁻⁵ mediando; non sparisce
   col reticolo fine (§7.5 è una cura possibile).
4. Solo **primo ordine in Zα**: niente Wichmann–Kroll (per Z = 60 è qualche % di Uehling, di segno opposto).
5. Nucleo puntiforme, massa infinita del nucleo, nucleo smussato su ~1/q_c per il potenziale nello spazio
   (emendamento A1: q_c = 0,15/θ invece di 0,3/θ, per un problema numerico dell'integrale, deciso prima dei numeri).
6. **Idrogeno**: il picco di contatto della regola radiale sul primo sito; servono θ_atom ≤ 0,00625 o il proiettore
   anche sull'elettrone. Pre-registrato FAIL.
7. Il metodo in spazio q ha un pavimento ~0,005 MHz (stati p).
8. **Parametri liberi: 0 adattati.** Ingressi: θ (= mₑ), α, Z. K1, K2: regolatore dichiarato (variato).

---

## 10. Patch proposte per altri file (non toccati)

* `docs/CONTESTO_CLOUD.md` §4: debito #2 → "QED: polarizzazione del vuoto fatta (VUOTO_POLARIZZA.md: Uehling 0,1 %,
  log 2α/3π −0,07 %, Nd⁵⁹⁺ 1s −0,5 %); restano autoenergia, g−2, 25 MHz". Debito #3 → aggiungere "la regola 3-D ha
  anche 32 coni senza massa a quasi-energia π; il mare intero 1-D non polarizza (S U S = −U)". Nuovo debito:
  anisotropia del vuoto.
* `docs/PIANO.md` §3: riga "Vuoto che polarizza (Uehling dal mare di Dirac): fatto, idrogeno solo a posteriori".
* `src/atomo/` (Laboratorio Atomo, altro thread): proposta di aggiungere a Nd⁵⁹⁺ lo spostamento di Uehling della
  regola (−11,6 eV sull'1s) e, come vista, l'alone del §7.4.
* Ponte: nessuna richiesta necessaria. Se si volesse una foto dell'alone: richiesta `foto` della scena Nd⁵⁹⁺ dopo che
  il thread di rendering avesse aggiunto il profilo.

---

## 11. Dove l'intuizione di Valerio può aiutare

1. **La moneta che respira.** Ti piace l'idea che la massa "vibri" fra tick pari e dispari? Fisicamente vuol dire che
   l'elettrone ha due masse che si alternano e una media. È l'unico modo che ho trovato per far nascere **da solo** un
   vuoto sano. Lo vuoi come regola nuova (allora decidiamo μ), o preferisci tenere il proiettore dichiarato?
2. **I coni senza massa.** La regola 3-D, oltre all'elettrone, contiene particelle di luce "nascoste" nel vuoto
   (coni senza massa a quasi-energia π). Ti fanno pensare a qualcosa (neutrini? fotoni del reticolo?), o sono da
   eliminare?
3. **Il vuoto birifrangente.** La regola vede direzioni privilegiate per colpa dell'ordine dei passi. Ti piacerebbe
   un ordine "a giro" (xyz, yzx, zxy) che le renda tutte uguali?
4. **L'autoenergia.** Per il Lamb vero l'elettrone deve emettere e riassorbire palline di luce. Preferisci provarlo
   "a conti" (loop in k, §8.1, veloce) o "a palline" nello spazio (più bello da vedere, ma ~5·10⁴ di scala da coprire)?

---

## 12. Come rieseguire

    python lab/vuoto_core.py                 # autotest (Riesz contro differenze finite)
    python lab/vuoto_dev_check.py            # sviluppo: k-spazio contro mare vero in spazio reale (~1 min)
    python lab/vuoto_run1d.py                # V1 (~10 s)
    sh lab/vuoto_drive.sh  (e, in parallelo, sh lab/vuoto_drive2.sh)   # tabelle 3-D, a pezzi da 10 min (~2,5 h in tutto)
    python lab/vuoto_analyze3d.py            # V2
    timeout 590 python lab/vuoto_levels.py   # V3, V4 (~5 min)
    sh lab/vuoto_driveB.sh; sh lab/vuoto_drivePH.sh                     # esplorazione B e sguardo a posteriori
    python lab/vuoto_figure.py
