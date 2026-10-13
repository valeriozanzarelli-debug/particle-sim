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

**Quello che in fase 1 NON era uscito:** l'autoenergia (≈ +1084 MHz del Lamb 2S–2P, quasi tutti i 7,3 GHz
dell'1S–2S), g−2, i 25 MHz della struttura fine: servono la luce quantizzata e il loop elettrone–fotone (piano al §8;
fatto in parte nella fase 2, §13: δm e g−2 sì, il Lamb completo non ancora).

**Fase 2 (§13), in cinque righe.**
1. **Autoenergia dalla regola** (elettrone-walk × luce-walk, loop in k): la massa vestita cresce come
   (3α/2π) m ln(1/ma), il coefficiente di Weisskopf, entro lo **0,06 % (PASS)**.
2. **g−2 della regola = α/2π entro il 2,4 % (PASS)**, ma dipende dall'orientazione (0,94–0,99, non converge: G3
   FAIL) perché il "c" dell'elettrone vestito è **anisotropo** (nuovo: viene dall'ordine dei passi xyz/zyx).
3. **Logaritmo di Bethe** dalla regola radiale: 1S uguale alla QM fino a 300 Ry entro 0,006 (ibrido ln k₀ = 2,978
   contro 2,984), ma **FAIL pre-registrato** (2S a θ = 0,05 rotto, due emendamenti dichiarati non lo curano).
4. **Moneta che respira** con μ = √θ (principio): il **mare intero, senza correzioni**, dà carica finita e Uehling
   (B5–B8 **PASS**), ma l'elettrone rallenta a c·cos μ (rottura di Lorentz θ/2, trovata anche da REGOLA_LOCALE):
   esplorazione, la regola **non** è cambiata.
5. **Idrogeno col proiettore anche sull'elettrone: PASS**, 1S–2S di vuoto +188,57 MHz contro +188,21 (+0,19 %).
Le 4 specie massive di questo studio sono **le 4 valli** di BATTIMENTO_1S; i 4 coni sono in più.

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

Numeri di riferimento della QED (solo esame, arrotondati). **Verificati sulla fonte in fase 2** (tabella dei
contributi CODATA-2018 nel materiale supplementare di arXiv:2211.05411, in kHz): autoenergia a un loop 1S
**8 383 339**, 2S **1 071 281**; polarizzazione del vuoto 1S **−215 168**, 2S **−26 897**; dimensione del protone
(termine di testa) 1S 1 107, 2S 138; Lamb totale 1S 8 172 847, 2S 1 045 009. **Non verificati** (scritti a
memoria, restano approssimati): Lamb 2S–2P ≈ 1057,8 MHz, autoenergia del 2P₁/₂ ≈ −12,8 MHz, parte del momento
anomalo ≈ +68 MHz. Le righe sotto usano i valori verificati dove ci sono:

| pezzo | 2S–2P (Lamb, 1057,8 MHz) | 1S–2S (la regola/Dirac sta +7,11 GHz sopra la misura, `NUCLEO_1S.md`) | chi lo può dare |
|---|---|---|---|
| **polarizzazione del vuoto (Uehling)** | **−27,1 MHz** | 1S −215,1, 2S −26,9 → **+188 MHz** (va nel verso opposto) | **la regola, col mare di Dirac: fatto qui** (2S −26,85, 1S −214,9 → +188,0 MHz, a posteriori sui reticoli fini; Nd⁵⁹⁺ PASS) |
| **autoenergia** (elettrone che emette e riassorbe un fotone, compreso il momento magnetico anomalo ≈ +68 MHz, non verificato) | ≈ **+1084 MHz** (2S: +1071,28 verificato; 2P₁/₂ ≈ −12,8 non verificato) | 1S +8383,34, 2S +1071,28 → **−7312,06 MHz** (verificato) | serve la **luce quantizzata** + il loop elettrone–fotone: **non c'è** (piano al §8) |
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

**7.7 (fase 2) Il cono di luce vestito, schiacciato.** *Cosa si vedrebbe:* un elettrone fermo che emette palline di
luce; il suo "fronte" di velocità massima non è una sfera ma un ellissoide leggermente schiacciato lungo x e allungato
lungo z, che si gira di 90° (x ↔ z) a ogni tick. *Perché:* G4, i δv_j misurati (−0,005; 0,044; 0,072) α. È piccolo
(~5·10⁻⁴), quindi andrebbe esagerato ×1000 e dichiarato come tale. *Costo:* nullo (una forma da densità). *Regole:* sì,
è una proprietà della regola. Solo proposta per il thread del motore grafico.

**7.8 (fase 2) Un solo rimedio per due birifrangenze.** Il vuoto (§4) e l'elettrone vestito (G4) sono anisotropi per
la stessa ragione: l'ordine dei passi. L'ordine "a giro" del §7.5 (6 permutazioni) li curerebbe tutti e due; la prova
più economica è il g−2 con le 6 permutazioni (lo stesso codice di `vuoto_f2_core.g2_all`, periodo 6 invece di 2).
Non provato: cambia la regola.

**7.9 (fase 2) Lo specchio dello spettro.** La regola radiale ha ogni livello ripetuto a fase +π (lo stesso fenomeno
delle valli): si potrebbe disegnare l'atomo con due "piani" sovrapposti, quello fisico e il suo gemello capovolto,
che lampeggiano in controfase. Utile per far vedere da dove vengono i problemi del logaritmo di Bethe (§13.2).

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
9. **Fase 2** (dettagli al §13.5): θ_nudo ≠ θ_fisico (δm/m ≈ 1,55 α); controtermine di velocità δv_j anisotropo;
   massa del fotone λ → 0; nucleo di Coulomb come ingresso nella parte istantanea; filtro di banda e due emendamenti
   numerici nel logaritmo di Bethe (2S a θ = 0,05 FAIL aperto); g−2 su una sola fetta di tempo. Il debito #6
   (picco di contatto dell'idrogeno) è **pagato** dal proiettore sull'elettrone (§13.4, nessuna manopola nuova).

---

## 10. Patch proposte per altri file (non toccati)

* `docs/CONTESTO_CLOUD.md` §4: debito #2 → "QED: polarizzazione del vuoto fatta (VUOTO_POLARIZZA.md: Uehling 0,1 %,
  log 2α/3π −0,07 %, Nd⁵⁹⁺ 1s −0,5 %); restano autoenergia, g−2, 25 MHz". Debito #3 → aggiungere "la regola 3-D ha
  anche 32 coni senza massa a quasi-energia π; il mare intero 1-D non polarizza (S U S = −U)". Nuovo debito:
  anisotropia del vuoto.
* `docs/PIANO.md` §3: riga "Vuoto che polarizza (Uehling dal mare di Dirac): fatto, idrogeno solo a posteriori".
* `src/atomo/` (Laboratorio Atomo, altro thread): proposta di aggiungere a Nd⁵⁹⁺ lo spostamento di Uehling della
  regola (−11,6 eV sull'1s) e, come vista, l'alone del §7.4.
* (fase 2) `docs/CONTESTO_CLOUD.md`: debito #2 → aggiungere "autoenergia nella regola: δm (3α/2π) ln, 0,06 %; g−2 =
  α/2π −2,4 % con dipendenza dall'orientazione; Bethe 1S ok fino a 300 Ry, 2S aperto". Nuovo debito: c vestito
  anisotropo (G4). `docs/PIANO.md`: riga "Vuoto fase 2: δm, g−2 PASS; Bethe FAIL; moneta che respira esplorata
  (non adottata: rompe Lorentz di θ/2)".
* (fase 2) `src/atomo/`: l'idrogeno col proiettore sull'elettrone dà il pezzo di vuoto del 1S–2S (+188,6 MHz):
  si può mostrare come numero accanto alla riga 1S–2S.
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

5. **(fase 2) La moneta che respira, il verdetto.** Con μ = √θ il vuoto intero viene giusto **da solo** (B5–B7
   PASS), ma l'elettrone corre a c·cos μ mentre la luce va a c. Accetteresti una regola dove luce ed elettrone hanno
   "c" diversi di θ/2 (che spariscono solo nel continuo)? O aspettiamo la variante in cui respira anche la luce
   (REGOLA_LOCALE fase 5)? Io voterei per aspettare.
6. **(fase 2) Il controtermine di velocità.** Per avere g−2 giusto devo dire all'elettrone vestito "corri a c". È
   un'aggiunta a mano, fissata da una condizione fisica. Ti sembra accettabile, o pensi che la regola dovrebbe
   garantirlo da sola (cioè: luce ed elettrone sullo stesso reticolo dovrebbero vestirsi senza cambiare c)?
7. **(fase 2) L'ordine dei passi.** Due anisotropie (vuoto e vestito) vengono dall'ordine xyz/zyx. Vale la pena
   provare l'ordine a giro (6 permutazioni)? Costa poco nei conti, ma è un cambio di regola.

---

## 12. Come rieseguire

    python lab/vuoto_core.py                 # autotest (Riesz contro differenze finite)
    python lab/vuoto_dev_check.py            # sviluppo: k-spazio contro mare vero in spazio reale (~1 min)
    python lab/vuoto_run1d.py                # V1 (~10 s)
    sh lab/vuoto_drive.sh  (e, in parallelo, sh lab/vuoto_drive2.sh)   # tabelle 3-D, a pezzi da 10 min (~2,5 h in tutto)
    python lab/vuoto_analyze3d.py            # V2
    timeout 590 python lab/vuoto_levels.py   # V3, V4 (~5 min)
    sh lab/vuoto_driveB.sh; sh lab/vuoto_drivePH.sh                     # esplorazione B e sguardo a posteriori
    # fase 2 (pre-registrazioni lab/out/vuoto_f2_prereg.json e vuoto_f2_prereg_L.json)
    python lab/vuoto_f2_nrref.py             # riferimento d'esame (QM non relativistica, ~5 min)
    sh lab/vuoto_f2_drive_se.sh              # δm e g−2 (~1 h 30)
    sh lab/vuoto_f2_drive_bethe.sh           # logaritmo di Bethe pre-registrato (~15 min)
    VUOTO_AMEND_L=1 sh lab/vuoto_f2_drive_bethe.sh   # emendamento A-L1;  VUOTO_AMEND_L=2 → A-L2
    sh lab/vuoto_f2_drive_breath.sh          # moneta che respira (~2 h)
    sh lab/vuoto_f2_drive_hproj.sh           # idrogeno col proiettore
    python lab/vuoto_f2_analyze.py           # → lab/out/vuoto_f2_exams.json
    python lab/vuoto_figure.py

---

## 13. Fase 2: il loop dell'autoenergia, la moneta che respira, l'idrogeno col proiettore

Pre-registrazioni: `lab/out/vuoto_f2_prereg.json` (commit 4be9feb, **prima** delle corse S, G, B, H) e
`lab/out/vuoto_f2_prereg_L.json` (commit d892181, **prima** del logaritmo di Bethe sull'idrogeno). Codice:
`lab/vuoto_f2_core.py` (motori), `lab/vuoto_f2_se.py`, `vuoto_f2_bethe.py`, `vuoto_f2_breath.py`, `vuoto_f2_hproj.py`
(corse a pezzi da ≤ 10 min, con i driver `lab/vuoto_f2_drive_*.sh`), `lab/vuoto_f2_nrref.py` (riferimento solo
d'esame), `lab/vuoto_f2_analyze.py` → `lab/out/vuoto_f2_exams.json`. Figura: pannelli (e), (f) di
`lab/out/vuoto_figure.png`.

**Le 4 valli di BATTIMENTO_1S sono le mie 4 specie massive.** Le 4 specie pesanti che il §2 trova in ogni cella π
del 3-D (k = 0 e k = (π/2, π/2, 0) con le permutazioni) sono **esattamente** le 4 valli di `docs/BATTIMENTO_1S.md`:
il nucleo puntiforme le mescola, quindi l'1s del Laboratorio si apre in un multipletto e ogni livello esatto sta
25/25/25/25 sulle 4 valli. I 4 coni senza massa a quasi-energia π (§2) sono **in più**: non entrano nel multipletto
dell'1s, ma sono loro che rovinano il vuoto intero (§3). La moneta che respira (sotto) lascia le 4 valli a massa θ
**esatta** (B8) e dà massa μ solo ai coni: quindi **non** cura il battimento, cura solo il vuoto.

### 13.1 Il loop dell'autoenergia nella regola (spazio k)

**Che cosa si calcola (semplice).** Un elettrone fermo emette una pallina di luce e la riprende. Nel frattempo
l'elettrone vola con impulso −k e la luce con +k. Si sommano tutti i k. Gli ingredienti sono **tutti della
regola**: le bande e gli spinori dell'elettrone sono quelli di Floquet della zig-zag 3-D (due tick); i fotoni, la
loro velocità e le loro due polarizzazioni vengono dal **cammino della luce** (`bloch_batch(spin=1)`), anche lui a
due tick; il vertice è "l'elettrone sente il potenziale vettore sullo spigolo" (Peierls: dH_eff/dk a metà strada,
H_eff = i log W / 2). Il taglio a grande k è lo stesso proiettore di specie della fase 1 (K1 = 0,6, K2 = 1,05).
Il nucleo di Coulomb 4πα/k² è un ingresso (come in fase 1). Formule nel gauge di Coulomb, forma di Feynman
(teoria dei buchi): parte di Coulomb istantanea + parte trasversa con i fotoni veri. **Prima** le stesse formule
sono state provate sul Dirac del continuo (devono dare il risultato del libro, 3α/2π · ln).

**δm (S1, S2): PASS.** δm/m = A + B ln(1/θ) + Cθ su θ = 0,05 … 0,003125:

| | B_Coulomb/(α/π) | B_trasv/(α/2π) | B_tot/(3α/2π) | soglia | esito |
|---|---|---|---|---|---|
| S1 continuo (controllo delle formule) | 1,0002 | 1,0003 | 1,0003 | 0,5 % | PASS |
| **S2 regola (elettrone-walk × luce-walk)** | **1,0004** | **1,0010** | **1,0006** | 2 / 3 / 2 % | **PASS** |

Cioè: la massa dell'elettrone della regola, vestita dalla sua luce, cresce come **(3α/2π) m ln(1/(m a))**, il
risultato di Weisskopf (1939), con 1/π dalla parte di Coulomb e 1/2π dai fotoni trasversi. **È una validazione**:
il coefficiente del logaritmo nasce tra m e 1/a, dove la regola è già il Dirac del continuo. La costante A
dipende dal taglio (proiettore) e non ha significato fisico: è il **debito θ_nudo ≠ θ_fisico** (δm/m ≈ 1,55 α a
θ = 0,05, cioè l'1,1 %: per avere la massa vera bisogna partire da un θ nudo un po' più piccolo).


**g−2 (G0–G4).** Stesso loop, ma l'elettrone viene **colpito** da un fotone esterno di impulso Q molto piccolo
(sistema di Breit, p = ∓Q/2). La parte "magnetica" del vertice dà F₂(0), il momento anomalo. Il libro dice
F₂ = α/2π (Schwinger 1948). Massa del fotone λ come regolatore infrarosso (λ/m = 10⁻³ e 3·10⁻⁴, poi λ → 0
lineare), poi θ → 0 con a + bθ + cθ².

*Due lezioni prima di arrivarci (sul continuo, prima della pre-registrazione).* Nel gauge di Coulomb con un taglio
sull'impulso a 3 dimensioni il vertice da solo **non** basta: dà 0,664. Servono altri due pezzi che di solito non si
vedono: (1) le **gambe esterne** (l'autoenergia sul guscio mescola l'elettrone con gli stati a energia negativa allo
stesso impulso, dopo aver tolto δm): +1,00; (2) un **controtermine di velocità** δv: il taglio rompe Lorentz, e
l'elettrone vestito correrebbe con un "c" un po' diverso. Lo fisso chiedendo che il vestito abbia E² − p² = m²: −0,33.
Con i tre pezzi il continuo dà **1,0000**.

| | θ = 0,05 | 0,025 | 0,0125 | θ → 0 | soglia | esito |
|---|---|---|---|---|---|---|
| G0 continuo (controllo) | 1,0034 | 1,0009 | 1,0002 | **0,99999** | 0,5 % | PASS |
| **G1 regola** (media delle 6 coppie asse-Q / asse-corrente) | 0,943 | 0,965 | 0,972 | **0,976** | 3 % | **PASS** |
| G2 regola a θ = 0,0125 | | | 0,972 | | 6 % | PASS |
| G3 dispersione fra le 6 orientazioni | 0,063 | 0,056 | 0,054 | | deve dimezzare | **FAIL** |

(valori in unità di α/2π, dopo λ → 0). Controllo di griglia (θ = 0,025, griglia più fine): 0,961004 contro
0,961002, cioè le somme su k sono convergenti a 2·10⁻⁶.

**Cosa vuol dire.** In media il **momento magnetico anomalo dell'elettrone della regola è α/2π entro il 2,4 %**:
l'elettrone-walk vestito dalla luce-walk ha il g−2 di Schwinger. È una **validazione** (numero del 1948), non una
scoperta. **Ma** la media nasconde una cosa: le orientazioni non sono uguali. A θ = 0,0125 le 4 coppie con la
corrente lungo x o z danno 0,983–0,994 (quasi giuste), le 2 con la corrente lungo y danno 0,94. E questa differenza
**non cala** quando θ → 0 (G3 FAIL): non è un errore del reticolo che sparisce, è una proprietà del loop della regola
vicino al taglio.

**G4, la scoperta collegata: il "c" vestito della regola è anisotropo.** Il controtermine di velocità che serve sulla
regola è diverso lungo i tre assi: δv/α = (−0,005; 0,044; 0,072) sulla fetta di tempo 0 (contro 0,053 uguale in
tutte le direzioni sul continuo), quasi indipendente da θ. Sulla fetta 1 x e z si scambiano: (0,074; 0,053; 0,007).
Cioè: **la luce della regola veste l'elettrone in modo diverso lungo gli assi, perché i passi vanno in ordine xyz e
poi zyx.** È lo stesso tipo di birifrangenza già trovato nel vuoto (fase 1, §4), ma qui l'anisotropia è grande quanto
l'effetto stesso (δv cambia di ~150 % fra gli assi, mentre nel vuoto era 7·10⁻⁴ in relativo); in assoluto resta
piccola: ~5·10⁻⁴ sulla velocità dell'elettrone. La media sulle due fette, (0,035; 0,049; 0,040), è molto più isotropa; il g−2 è
stato calcolato solo sulla fetta 0 (debito, sotto). La cura naturale è la stessa del §7.5: l'ordine degli assi a giro
(6 permutazioni), che renderebbe isotropi sia il vuoto sia il vestito.

### 13.2 Il logaritmo di Bethe dal risolvente della regola radiale (L1–L3: FAIL pre-registrato)

**Che cosa si calcola (semplice).** La parte grossa del Lamb è "a bassa energia": l'elettrone legato emette un
fotone morbido e passa per un momento su un altro stato. Quanto costa dipende da una media logaritmica delle energie
di salto, pesata con quanto è facile saltare: il **logaritmo di Bethe** ln k₀ (1947). Qui lo prendo dal
**risolvente della regola radiale**: si parte dallo stato s della regola, si dà un calcio di dipolo (canali p,
κ' = +1 e −2), si tengono solo gli stati a energia positiva della banda fisica, e si somma su tutte le energie del
fotone fino a un taglio K.

**Come si giudica.** Il valore dei libri (2,984 per l'1S, 2,812 per il 2S, K → ∞) **non** è un esame pulito per un
Dirac sul reticolo (sopra ~m la teoria non relativistica non vale più, e il reticolo ha la sua coda). Quindi
l'esame pre-registrato confronta la regola con la **meccanica quantistica non relativistica allo stesso taglio K**
(`lab/vuoto_f2_nrref.py`, Schrödinger a differenze finite, due griglie + Richardson; solo esame). Soglie: L1
|Δ| ≤ 0,02 a 100 e 300 Ry (1S e 2S, θ = 0,05); L2 |L(0,05) − L(0,1)| ≤ 0,01; L3 la "forza di dipolo" D_K entro 1 %.

| ln k₀(K) regola − QM | 30 Ry | 100 Ry | 300 Ry | 1000 Ry | 3000 Ry | D_K regola/QM a 300 Ry |
|---|---|---|---|---|---|---|
| **1S, θ = 0,05** | +0,0003 | −0,0015 | −0,0063 | −0,019 | −0,041 | 0,968 |
| 1S, θ = 0,1 | +0,008 | +0,012 | +0,014 | +0,011 | −0,002 | 0,986 |
| 2S, θ = 0,1 | +0,011 | +0,017 | +0,026 | +0,051 | +0,13 | 0,990 |
| 2S, θ = 0,05 | +0,11 | +0,37 | **+1,00** | +2,3 | +3,8 | 1,47 (rotto) |

**Esito: L1, L2, L3 FAIL.** Il pre-registrato resta questo.
* L'**1S** a θ = 0,05 segue la QM fino a 300 Ry entro 0,006 (sotto la soglia), ma D_K è basso del 3 % (L3 FAIL) e la
  differenza fra θ = 0,05 e 0,1 a 300 Ry è 0,020 (L2 FAIL). Non so perché D_K peggiori col reticolo più fine: lo
  dichiaro come non capito.
* Il **2S a θ = 0,05 è rotto**: a grandi K la somma esplode (D_K cresce come ~K^0,6). Le regole dicono che il
  risultato pre-registrato non si tocca; ho provato due emendamenti **dichiarati, scritti prima di ciascuna corsa e
  non valutati**:
  - **A-L1** (`lab/out/vuoto_f2_amend_L.json`): quadratura più stretta del proiettore sulle energie positive
    (ipotesi: una perdita verso il mare). Il 2S a θ = 0,1 migliora (+0,13 → +0,045 a 3000 Ry), il 2S a θ = 0,05
    **peggiora** (+1,91 a 300 Ry). Ipotesi sbagliata o incompleta.
  - **A-L2** (`lab/out/vuoto_f2_amend_L2.json`): misurando lo spettro ho trovato che la regola radiale ha uno
    **spettro a specchio**: lo stato degenere col 2S (il 2P₁/₂, distanza 1,6·10⁻¹⁰) ha un gemello a fase +π, che
    cade proprio sul polo della trasformata di Cayley usata nel risolvente (amplificazione ~10¹⁰). Tolti
    esplicitamente i 6 stati più vicini al polo: **non basta** (+1,54 a 300 Ry). Il gemello è vero (ed è lo stesso
    fenomeno delle valli/doppioni), ma non è l'unico colpevole.
  Mi fermo qui: un terzo emendamento sarebbe "pescare" il risultato. Il 2S a θ = 0,05 resta **un FAIL aperto**.

**Cosa resta di buono (a posteriori, non valutato).** Unendo la regola sotto K = 300 Ry con la teoria dei libri
sopra K ("ibrido"): **ln k₀(1S) = 2,978** contro 2,984 (−0,006, cioè −6,5 MHz sul Lamb dell'1S, che vale ~1085 MHz
per unità di ln k₀) e ln k₀(2S) = 2,836 (θ = 0,1) contro 2,812 (+3,4 MHz). Quindi la **fisica a bassa energia del
Lamb dell'1S c'è nella regola**, al livello dello 0,2 %. Onestà: la parte ad alta energia (sopra K) qui è presa dai
libri, non dalla regola; per chiudere il conto serve unire questo pezzo al loop in k del §13.1 (stesso taglio), lavoro
non fatto.

### 13.3 La moneta che respira con μ scelto da un principio (esplorazione, la regola NON cambia)

**Il principio (pre-registrato): μ = √θ.** La massa dei regolatori M = μ/a è la **media geometrica** fra la massa
dell'elettrone m = θ/a e il taglio del reticolo 1/a: il regolatore sta "a metà strada" (in scala logaritmica) fra
la fisica e i quadretti. Così si bilanciano i due errori della moneta che respira, tutti e due di ordine θ: il
termine di Uehling del regolatore stesso, (m/M)² = θ, e la deformazione della specie fisica, di ordine μ² = θ. Non
è stato scelto guardando le soglie.

**Che cosa fa alle bande (B8: PASS, errore 4·10⁻¹⁷).** Le 4 valli pesanti restano a massa θ **esatta**; i 4 coni
prendono massa μ esatta (0,1581 a θ = 0,025; 0,1118 a θ = 0,0125). Ma c'è un prezzo che in fase 1 non avevo visto:
**l'elettrone rallenta.** La velocità limite della specie fisica diventa **cos μ** (misurata: v² = 0,98918 con
μ = 0,1, cos²μ = 0,99003 a meno dell'errore θ² del reticolo), uguale in tutte le direzioni, mentre la luce resta a c
(il cammino della luce non ha moneta). Con μ = √θ: 1 − v ≈ θ/2.
**Questo è una rottura di Lorentz** (l'elettrone e la luce non hanno più lo stesso c, differenza μ²/2): lo studio
REGOLA_LOCALE (fase 4, `docs/REGOLA_LOCALE.md` §13) l'ha trovata in modo indipendente (μ = 0,1 → elettrone a
0,995 c). Con μ **fisso** non sparisce raffinando il reticolo. Con il principio di questo studio, μ = √θ, invece
1 − v ≈ θ/2 va a zero come la massa in unità di reticolo, cioè **sparisce nel limite del continuo**, ma lentamente
(O(θ), non O(θ²) come gli altri errori della regola). La variante in cui anche la luce "respira" dello stesso cos μ
(così i due c tornano uguali) la prova REGOLA_LOCALE in fase 5; se funziona, il prezzo sugli atomi qui sotto andrebbe
rifatto.

**Risultati pre-registrati (B5, B6, B7: tutti PASS).** Mare intero, **senza proiettore**, μ = √θ, media su x, y, z:

| θ | μ = √θ | B5: Π(0) / (1/6π²) ln(μ/θ) | B6: Π_R / Uehling-PV, q = 1 m | q = 3 m | B7: specie fisica con/senza respiro, q = 1 m | q = 3 m |
|---|---|---|---|---|---|---|
| 0,025 | 0,158 | 1,0104 | 0,9925 | 0,9995 | 0,9946 | 1,0020 |
| 0,0125 | 0,112 | 1,0050 | 0,9960 | 0,9995 | 0,9966 | 1,0003 |
| soglia | | 3 % | 2 % | 2 % | 1 % | 1 % |

Cioè: con la moneta che respira e μ scelto **prima** dal principio, il vuoto intero della regola dà da solo una
carica rinormalizzata **finita** (errore 1,0 % → 0,5 %: cala come θ, come detto dal principio) e la forma di Uehling
regolata alla Pauli–Villars entro lo 0,4–0,8 % a q = m, 0,05 % a 3 m. Il prezzo sulla specie fisica (B7) è lo
0,3–0,5 % a q = m, come previsto dalla fase di sviluppo. In fase 1 le soglie B erano fallite perché i μ erano troppo
grandi (0,3–0,6); con il principio μ = √θ passano. **Questo è il primo vuoto sano della regola senza correzioni a
mano**, ma a una condizione: cambiare la regola (sotto, il prezzo).

**Attenzione: la deformazione non è solo "velocità".** Ho provato a spiegare B7 come Uehling con c → cos μ
(predizione Π_R(q cos μ)/cos μ): **sbagliato**, non torna con i numeri misurati (sviluppo, a θ = 0,05). La specie fisica con il respiro è un
Dirac deformato in modo più complicato (anche gli spinori cambiano). Lo dichiaro: il numero B7 è misurato, non capito
fino in fondo.

**Cosa cambierebbe per gli altri studi se diventasse regola** (non l'ho cambiata; un cambio di regola si coordina
con tutti i thread):
* **Vuoto (questo studio):** niente proiettore di specie; la carica rinormalizzata esce finita da sola,
  Π(0) = (2α/3π) ln(M/m) con M/m = 1/√θ (a θ = 0,0125: ln 8,9 = 2,19). Prezzo: errore O(θ) invece di O(θ²).
* **Atomi (ATOMO, NUCLEO_1S, BATTIMENTO, Laboratorio):** con velocità v = cos μ la massa cinetica diventa m/v²
  e le energie di legame crescono di 1/v² − 1 ≈ μ² = **θ in relativo** (a θ = 0,05: **+5 %** sui legami!). Andrebbe compensato (per esempio ritoccando la regola) o si perde la precisione attuale degli
  atomi. **Le 4 valli restano** (massa θ esatta): il multipletto dell'1s di BATTIMENTO non cambia.
* **Regola radiale e 1-D:** in 1-D il doppione a k = π **non** prende massa (a k = π le due monete commutano:
  C(θ+μ)C(θ−μ) = C(2θ)); cambia solo la velocità (cos μ). Nessun guadagno, solo il prezzo.
* **Simulazioni in spazio reale (EMISSIONE, LOCALE, motore):** i modi spuri senza massa a quasi-energia π (i coni)
  diventano massivi: meno "sfarfallio" stroboscopico a lunga distanza. Costo zero (un numero per tick).
* **Luce:** nessun cambiamento.

**Verdetto onesto.** La moneta che respira è un'idea bella e funziona sul vuoto, ma rallenta l'elettrone di
1 − cos μ: per i conti di precisione sugli atomi costa più di quanto rende. Il proiettore dichiarato resta la scelta
migliore finché non si trova una "respirazione" che lasci intatta la velocità.

### 13.4 Idrogeno col proiettore anche sull'elettrone (H1–H4: tutti PASS)

**Che cosa cambia rispetto alla fase 1.** In fase 1 il FAIL dell'idrogeno (+14 %) veniva dal "picco di contatto":
la funzione d'onda della regola radiale ha una spina sul primo sito vicino al nucleo, e lo spostamento di Uehling
(che vive proprio lì, entro un λ_C) la sente tutta. Ora lo **stesso** proiettore di specie del vuoto (K1, K2: nessun
numero nuovo) si applica anche allo stato dell'elettrone prima di calcolare ⟨V_Uehling⟩. Pre-registrato (H1–H4)
prima della corsa.

| | θ_atom = 0,05 | 0,025 | 0,0125 | estrapolato (E₀ + cθ + dθ²) | riferimento | esito |
|---|---|---|---|---|---|---|
| 1S (MHz) | −211,65 | −213,13 | −214,20 (−0,42 %) | **−215,50 (+0,19 %)** | −215,10 | H1, H3 PASS |
| 2S (MHz) | −26,456 | −26,639 | −26,772 (−0,43 %) | **−26,934 (+0,17 %)** | −26,887 | H2, H3 PASS |
| 1S–2S (MHz) | | | | **+188,57 (+0,19 %)** | +188,21 | H4 PASS |

Soglie: 1 % sul reticolo più fine (H1, H2), 0,5 % sull'estrapolato (H3), 1 % sull'1S–2S (H4). I riferimenti sono la
parte di polarizzazione del vuoto di CODATA-2018 (arXiv:2211.05411, verificati in tabella: −215,168 e −26,897 MHz;
il valore usato è l'Uehling puro a nucleo puntiforme, la differenza è < 0,04 %).

**Onestà sulla convergenza.** Le differenze fra un reticolo e il successivo calano di un fattore ~1,37 per ogni
dimezzamento di θ_atom, cioè come **θ^0,46** (non θ, non θ²). Il fit E₀ + cθ + dθ² su tre punti quindi è un po'
ottimista: se le differenze continuassero a calare di 1,37 (serie geometrica) l'1S arriverebbe a −217,1 MHz, cioè
+0,9 % oltre il riferimento. La verità sta fra i due estrapolati. Il risultato robusto è il punto più fine (−0,42 %),
già entro la soglia; l'estrapolato va preso con un'incertezza di ±0,7 %. Lo dichiaro come avvertenza, non cambia il PASS
pre-registrato. Il debito #6 del §9 (il picco di contatto) è **pagato** con un correttivo già dichiarato (nessuna
manopola nuova).

### 13.5 Riassunto fase 2, debiti, manopole

| esame | che cosa | esito |
|---|---|---|
| S1, S2 | δm: pendenza del logaritmo 3α/2π (Weisskopf), continuo e regola | PASS (0,03 % e 0,06 %) |
| G0 | g−2 sul continuo (controllo delle formule) | PASS (0,001 %) |
| G1, G2 | g−2 della regola, θ → 0 e a θ = 0,0125 | PASS (−2,4 % e −2,8 %) |
| G3 | le orientazioni convergono | **FAIL** (dispersione 0,063 → 0,054) |
| G4 | c vestito anisotropo | trovato (nuovo debito) |
| L1, L2, L3 | logaritmo di Bethe contro QM allo stesso taglio | **FAIL** (2S a θ = 0,05 rotto; 1S buono fino a 300 Ry) |
| B5–B8 | moneta che respira, μ = √θ | PASS (ma rallenta l'elettrone: rottura di Lorentz θ/2) |
| H1–H4 | idrogeno col proiettore anche sull'elettrone | PASS (1S–2S +188,57 contro +188,21 MHz) |

**Debiti nuovi della fase 2** (anche nel §9):
* **θ_nudo ≠ θ_fisico**: la massa vestita è più grande di δm/m ≈ 1,55 α (a θ = 0,05); controtermine di massa.
* **Controtermine di velocità δv_j** (uno per asse): serve perché il taglio rompe Lorentz; sulla regola è anisotropo
  (G4). Il g−2 è calcolato sulla sola fetta di tempo 0.
* **Massa del fotone λ** come regolatore infrarosso (estrapolata a zero; dichiarata).
* Il **nucleo di Coulomb 4πα/k²** della parte istantanea è un ingresso (come in fase 1), non viene dalla luce-walk.
* Logaritmo di Bethe: **filtro di banda** (Ω < Ω_n + 0,6) e proiettore di specie sulla sorgente; due emendamenti
  numerici dichiarati; il pezzo sopra K dei libri nell'ibrido.
* Moneta che respira: **un numero μ in più** se diventasse regola (scelto da un principio, non adattato).

**Manopole adattate: 0.** Ingressi: θ (= mₑ in unità del reticolo), α, Z. Regolatori dichiarati: K1, K2 (ereditati),
δm e δv_j (fissati da condizioni fisiche: massa e dispersione del vestito), λ → 0, μ = √θ (solo esplorazione).

**Cosa resta da fare.**
1. Il 2S del logaritmo di Bethe a θ = 0,05 (un FAIL aperto: lo specchio del 2P non basta a spiegarlo).
2. Unire il pezzo a bassa energia (§13.2) e quello ad alta energia (§13.1) con lo stesso taglio → il Lamb 2S–2P
   (≈ +1058 MHz, valore non riverificato) tutto dalla regola. È il passo che manca per dire "Lamb emerso".
3. Il g−2 su tutte e due le fette di tempo (e con l'ordine degli assi a giro, §7.5) per vedere se G3 si chiude.
4. Coordinare con REGOLA_LOCALE la variante "respira anche la luce" (fase 5 di quello studio).

**Valori di autoenergia verificati sulla fonte** (CODATA-2018, arXiv:2211.05411, tabelle in kHz): autoenergia 1S
8 383 339,47 e 2S 1 071 281,00; polarizzazione del vuoto 1S −215 167,77 e 2S −26 897,00; raggio del nucleo 1S
1107,15 e 2S 138,39; Lamb totale 1S 8 172 847 e 2S 1 045 009,1. **Non riverificati** (scritti a memoria e quindi
dichiarati tali): Lamb 2S–2P 1057,8 MHz, autoenergia del 2P₁/₂ −12,8 MHz, parte del momento anomalo +68 MHz, e i
valori dei libri di ln k₀ (2,984; 2,812).
