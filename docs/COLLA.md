# La colla dei nuclei: che cosa manca perché i nodi si leghino quanto i nuclei veri

**Da dove partiamo.** In `NUCLEI_DA_NODI.md` i nodi del fluido di palline (skyrmioni) si mettono insieme da soli
in forme bellissime e giuste (ciambella, tetraedro, cubo) con E/B entro lo 0.4 % dalla letteratura. Ma si **legano
troppo**: il deuterone di 72 MeV invece di 2.2, il ⁴He di 310 invece di 28.3. **La domanda di Valerio:** "allora
mancano dei pezzi: completa il quadro onestamente, così l'agglomeratore sa cosa decidere a livello dei nuclei".

Qui abbiamo provato, **uno per uno e con gli esami scritti prima**, i pezzi che la letteratura indica:

| strada | che cosa aggiunge | parametri nuovi |
|---|---|---|
| C1 | la **massa del pione** (138 MeV, una costante della natura) | 0 liberi (le 2 unità si rifittano su M_N, M_Δ) |
| C2a | la **trottola quantistica**: spin e isospin del nucleo come rotazioni quantizzate | 0 |
| C2b | il **tremolio di punto zero**: le vibrazioni del nodo non possono stare ferme | 0 |
| C3a | un potenziale in più che rende i nodi "leggermente legati" (Gillard–Harland–Speight 2015) | **1** (α) |
| C3b | il termine "sestico" (quello del mesone ω) | **1** (c₆) |

Pre-registrazione: `lab/colla_prereg.py` → `lab/out/colla_prereg.json`, **committata (78d76a7) prima delle corse**
(cosa avevamo visto prima è scritto nel suo docstring). Codice: `lab/colla_core.py` (motore: energia generalizzata,
metrica cinetica, tensori d'inerzia, hessiano), `lab/colla_units.py` (unità), `lab/colla_run.py` (rilassamenti a pezzi
da 8 minuti con checkpoint), `lab/colla_quant.py` (C2a), `lab/colla_modes.py` (C2b), `lab/colla_summary.py` (esami),
`lab/colla_figure.py`. Numeri: `lab/out/colla_*.json`, log in `lab/out/colla_logs/`. Figura: `lab/out/colla_figure.png`.
Motore dei nodi (`knot_core.py`, `knot2_core.py`, `nuclei_core.py`) usato in sola lettura, schema c4, h = 0.2.

<!-- TABELLONE -->

## Le cose da sapere per capire (in breve)

* **Energia di legame** = (somma delle masse dei pezzi separati) − (massa del nucleo). Se è positiva il nucleo sta
  insieme. In tutto lo studio i due termini sono calcolati **sullo stesso reticolo e con lo stesso modello**, così gli
  errori del reticolo (≈ 1 %) si cancellano quasi del tutto nella differenza.
* **Unità.** Il modello in "unità di Skyrme" non ha numeri liberi; per avere MeV servono 2 misure. Come Adkins, Nappi
  e Witten usiamo la massa del nucleone (939 MeV) e della Δ (1232 MeV). Senza massa del pione vengono F_π = 129 MeV,
  e = 5.44 (ANW 1983); con la massa del pione vengono **e = 4.843, F_π = 108.3 MeV** (Adkins–Nappi 1984: 4.84 e 108).
* **Perché la Δ?** Il nucleone e la Δ sono **lo stesso nodo** che gira in due modi diversi (spin ½ e spin 3/2):
  M = E_classica + ħ² J(J+1)/(2Λ), con Λ il momento d'inerzia del nodo. La differenza Δ − N fissa Λ e quindi
  quanto "pesa" la trottola. **Questo è il punto chiave di C2a**: nei 939 MeV del nucleone ci sono già circa 75 MeV di
  trottola di spin.

## C1 — La massa del pione (costante della natura, 0 parametri liberi)

Il pione vero ha una massa (138 MeV): nel modello è un termine 2m²(1 − σ) in più, con m = 2 m_π/(e F_π) = 0.526 in
unità di Skyrme. I nodi diventano più piccoli e la loro coda muore in fretta (e^(−m r) invece di 1/r²).

Unità rifittate su M_N e M_Δ, in modo autoconsistente (m dipende da e e F_π): **e = 4.8426, F_π = 108.26 MeV**,
m = 0.5265, unità di energia F_π/4e = 5.589 MeV, di lunghezza 0.753 fm, ħ = 2e² = 46.9 (Adkins–Nappi 1984: e = 4.84,
F_π = 108 MeV → **C1a PASS**, scarti +0.05 % e +0.24 %).

Poi abbiamo rilassato i nodi B = 1…4 con il pione massivo, partendo dai campi senza massa di `NUCLEI_DA_NODI` più un
piccolo rumore liscio (3 %, per lasciare cadere eventuali selle simmetriche):

| B | E/12π²B (reticolo) | B misurato | raggio barionico | forma |
|---|---|---|---|---|
| 1 | 1.2920 | 0.9946 | 0.66 fm | sfera |
| 2 | 1.2404 | 1.9918 | 0.91 fm | **ciambella** (simmetrie 30°/45°/90° attorno all'asse: asym ≤ 0.003; centro vuoto 1 %) |
| 3 | 1.2063 | 2.9889 | 1.06 fm | **tetraedro** (rotazioni di 180° e 120°: asym < 10⁻⁴; inversione: 0.85, quindi non un cubo) |
| 4 | 1.1789 | 3.9878 | 1.18 fm | **cubo** (90°, 120°, inversione: asym < 10⁻⁴) |

* **C1c PASS**: con la massa vera del pione le forme restano le stesse (come dice Battye–Sutcliffe 2006 fino a B ≈ 8).
* **C1b FAIL (dichiarato)**: il B = 1 sul reticolo viene −1.2 % sotto il valore esatto 1D (soglia 1 %). Abbiamo
  verificato che è l'errore del passo h = 0.2: con la massa il nodo è più piccolo e i gradienti più ripidi (il
  profilo esatto messo sul reticolo dà −1.05 % a h = 0.2 e −0.48 % a h = 0.15). Nel legame l'errore si cancella
  (stesso reticolo per tutti i B).
* **Legame classico con il pione massivo** (unità AN): **deuterone 68 MeV, ³He 170 MeV, ⁴He 300 MeV** (contro 72,
  178, 310 senza massa). Cambia pochissimo: **C1d PASS** (la previsione era "la massa del pione da sola non basta").
  In frazione della massa: 4.0 %, 6.6 %, 8.8 % (veri: 0.12 %, 0.28 %, 0.76 %).

## C2a — La trottola quantistica (0 parametri)

Un nodo non è una pallina ferma: può **ruotare nello spazio** (spin J) e **ruotare nello spazio dell'isospin**
(isospin I: protone ↔ neutrone). In meccanica quantistica queste rotazioni sono quantizzate, e quanto costano dipende
dai **tensori d'inerzia** U (isospin), V (spin), W (misto), che abbiamo misurato direttamente sui campi rilassati del
reticolo (`colla_core.inertia`). Le regole (Finkelstein–Rubinstein, già in letteratura) dicono quali stati sono
permessi:

* nucleone: J = I = ½ → ħ²·3/(8Λ);
* deuterone: lo stato più basso permesso della ciambella è **J = 1, I = 0** (Braaten–Carson 1988) → ħ²/V_⊥
  (V_⊥ = inerzia di spin attorno a un asse che giace nel piano della ciambella);
* ³He/³H: J = I = ½ con "spin + isospin del corpo = 0" (come il riccio) → ¾(u + v + 2w)/(2(uv − w²)) ħ²;
* ⁴He: **J = I = 0 → nessuna energia di rotazione**.

**Il controllo dei tensori.** Per la ciambella la teoria dice U₁₁ = U₂₂, V₁₁ = V₂₂, V₃₃ = 4U₃₃, |W₃₃| = 2U₃₃,
W₁₁ = 0 (girare la ciambella attorno al suo asse equivale a girarla due volte nell'isospin). Sul reticolo:
U₁₁/U₂₂ − 1 = 3·10⁻⁴, V₃₃/U₃₃ = 3.9925, |W₃₃|/U₃₃ = 1.9977, W₁₁/U₁₁ = 2·10⁻⁶: **C2a-ii PASS**. Per il B = 1 invece
il momento d'inerzia sul reticolo (Λ = 59.4) è il **5.6 % sotto** quello esatto 1D (62.9): **C2a-i FAIL** (soglia 3 %):
l'inerzia "sente" la coda del nodo, che nella scatola (mezza larghezza 4.7, cioè 3.5 fm) è un po' tagliata, e il passo
del reticolo. Effetto sui numeri sotto: qualche MeV, non cambia nessuna conclusione.

| | inerzie (unità di Skyrme) | energia di rotazione |
|---|---|---|
| nucleone (B=1) | Λ = 59.4 (U = V = −W, il riccio) | **77.6 MeV** (M_N sul reticolo 932.8, M_Δ 1243) |
| deuterone (B=2, J=1, I=0) | V_⊥ = 193.3, V_asse = 342.0, U_⊥ = 124.0, U_asse = 85.7 | **63.6 MeV** |
| ³He (B=3, J=I=½) | u = 157.8, v = 507.1, w = +105.9 | **58.8 MeV** |
| ⁴He (B=4, J=I=0) | — | **0** |

**Il risultato, e perché.** Legame = B·M_N − M_nucleo, con M = energia classica + rotazione:

| | classico (C1) | con la trottola quantistica (C1 + C2a) | vero |
|---|---|---|---|
| deuterone | 68 MeV | **160 MeV** | 2.2 |
| ³He | 170 MeV | **344 MeV** | 7.7 |
| ⁴He | 300 MeV | **610 MeV** | 28.3 |

La quantizzazione "a corpo rigido" **peggiora** il legame (**C2a-iii PASS**: lo avevamo previsto). Il motivo è
semplice: ogni nucleone libero porta con sé ~78 MeV di trottola (è dentro i 939 MeV), mentre il nucleo che ruota come
un corpo unico ne porta poca (deuterone: 64 MeV per tutti e due) o niente (⁴He, spin 0). Quindi separare i pezzi costa
ancora di più. Il deuterone a 160 MeV coincide con quello di **Braaten–Carson 1988 (158 MeV)** — **C2a-iv PASS** —
che fu proprio il primo a trovare questo problema.

**Nota onesta sul ³He (correzione dopo aver visto i numeri).** Il segno di w conta. Il primo script assumeva w < 0
"come il riccio" e dava 30.4 MeV di rotazione. Ma per il lemma di Schur W = c·R (R una rotazione propria), e il segno
giusto è c = segno(det W)·|W|: per il riccio viene c = −Λ (giusto), per il tetraedro c = **+105.9**. Con il segno
giusto la rotazione del ³He è 58.8 MeV e il legame 344 MeV (con il segno sbagliato sarebbe 372). Lo stato con
"spin + isospin del corpo = 1" non è permesso dalla simmetria del tetraedro (darebbe 363 MeV): non cambia il verdetto.

**Che cosa ci dice la letteratura (verificata).** Il problema non è la quantizzazione in sé ma la "rigidità".
Leese–Manton–Schroers 1995 (hep-ph/9502405) hanno lasciato muovere il deuterone su uno spazio di 10 dimensioni (i due
nucleoni possono allontanarsi, girarsi, ecc.) e risolto l'equazione di Schrödinger lì sopra: **legame 6.2 MeV, raggio
2.18 fm** (il vero: 2.22 MeV e 2.14 fm). Il deuterone vero **non è una ciambella rigida**: è una nuvola quantistica di
due nucleoni quasi liberi che passano solo una parte del tempo vicini. Noi **non** l'abbiamo rifatto (vedi debiti).

## C2b — Il tremolio di punto zero (0 parametri)

Ogni modo di vibrazione con frequenza ω ha un'energia minima ħω/2 (non può stare fermo). Abbiamo cercato i modi
normali del B = 1 e del B = 2 risolvendo il problema agli autovalori K ξ = ω² (2M) ξ sul reticolo (K = hessiano
dell'energia sulla sfera, M = metrica cinetica del modello), con LOBPCG e un precondizionatore FFT (`colla_modes.py`).
Contano come "modi legati" quelli sotto la soglia del pione (ω < m): sopra c'è il continuo (pioni che scappano).

Risultati (pione massivo, unità AN, 48³, h = 0.2):

| | modi zero (≈ 0) | modi legati (ω < m_π) | primo modo sopra soglia |
|---|---|---|---|
| B = 1 | **6** (3 traslazioni, ω² = 0.0008 m²; 3 rotazioni, 0) | **nessuno** | ω² = 1.86 m² |
| B = 2 | **8** (3 traslazioni, 3 isorotazioni, 2 rotazioni dell'asse: queste a 0.013 m², il reticolo cubico le frena un po') | **1 doppietto**: ω² = 0.350 m², cioè ħω = 81.6 MeV | ω² = 1.68 m² |

* **C2b-i PASS** (6 e 8 modi zero, come deve essere: per la ciambella girarla attorno al suo asse è la stessa cosa che
  isoruotarla, quindi 3 + 3 + 3 − 1 = 8).
* Il nucleone con il pione vero **non ha vibrazioni legate** (il suo "respiro" sta nel continuo). La ciambella ha un
  doppietto legato: è uno dei modi in cui i due nucleoni dentro la ciambella "si muovono uno rispetto all'altro".
* Punto zero: B = 2 → ½·2·81.6 = **81.7 MeV**; due nucleoni → 0. Il deuterone si **slega di 82 MeV**: dopo C2a + C2b
  il legame è **78 MeV**, ancora 35 volte il vero → **C2b-ii PASS** (avevamo previsto "riduce, ma non basta").
* **Attenzione, è una stima a metà.** L'energia di punto zero vera è la somma su *tutti* i modi, anche quelli del
  continuo, che sono infiniti e vanno "rinormalizzati" (servono nuovi parametri: i controtermini). Per il B = 1 la
  letteratura (Moussallam) stima −200 MeV circa; per la ciambella 3D **nessuno l'ha mai calcolata**
  (Walliser–Holzwarth 1999). Quindi il segno e la grandezza della correzione completa al deuterone **non sono noti**.

## C3a — Nodi "leggermente legati" (1 parametro: α)

**L'idea (Gillard–Harland–Speight 2015, NPB 895, 272).** Si sposta un po' di peso dal termine "di gradiente" a un
potenziale α(1 − σ)⁴, che si annulla nel vuoto e cresce molto al centro del nodo:
E = ∫ (1−α)(|∂φ|² + 2m²(1−σ)) + Skyrme + α(1−σ)⁴, con m = 1 come loro. Per α → 1 c'è un limite inferiore
all'energia che solo il B = 1 raggiunge: i nodi grandi **quasi** non si legano più. α è **un parametro in più**.
Calibrazione come GHS: E₁ (classica) := 938.27 MeV (loro, per far tornare anche il raggio, hanno bisogno di
F_π = 29 MeV e m_π = 303 MeV, lontanissimi dai valori veri 186 MeV, nella loro convenzione, e 138 MeV: un debito che
loro stessi dichiarano).

**Controllo con la letteratura.** α = 0.95, B = 1 in 1D: E/8π² = 1.1093 (GHS 1.1082, +0.1 %); sul reticolo il
legame per nucleone del B = 2 è 1.58 MeV (GHS 1.84) e del B = 4 4.15 MeV (GHS 4.56): **C3a-i PASS** (soglia 30 %).
(Il punto α = 0 è il modello standard con m = 1: il cubo viene E/12π²B = 1.279 contro 1.307 di Battye–Manton–Sutcliffe:
**C1e FAIL**, −2.1 % per l'errore del reticolo a h = 0.2, che con m = 1 è più grande; il rapporto (E₄/4)/E₁ invece
è 0.9192 contro 0.9230, −0.4 %.)

**La scansione** (rilassamenti dal campo senza massa + rumore liscio; 52³ palline, h = 0.2):

| α | deuterone | ⁴He | forma B=2 | forma B=4 |
|---|---|---|---|---|
| 0 | 67.5 MeV (3.6 %) | 303 MeV (8.1 %) | ciambella | cubo |
| 0.5 | 35.2 MeV (1.9 %) | 163 MeV (4.3 %) | **2 grumi** (q₂ = 0.42) | **4 grumi a tetraedro** (q₃ = 0.23, niente inversione) |
| 0.8 | 14.3 MeV (0.76 %) | 69.3 MeV (1.8 %) | 2 grumi | 4 grumi a tetraedro |
| 0.95 | 3.2 MeV (0.17 %) | 16.6 MeV (0.44 %) | 2 grumi a 2.7 fm | 4 grumi a tetraedro (q₃ = 0.49) |
| vero | 2.22 MeV | 28.3 MeV | nuvola larga ~2 fm | — |

* **C3a-ii PASS**: il legame scende con α in modo regolare (di un fattore ~20 da 0 a 0.95).
* **C3a-iii FAIL** (metà giusta, metà sbagliata): avevamo previsto che a α = 0.5 ciambella e cubo restassero e che a
  0.95 si rompessero in grumi. A 0.95 sì, ma **già a α = 0.5 la ciambella è diventata due grumi e il cubo quattro
  grumi su un tetraedro**. Le forme "belle" di `NUCLEI_DA_NODI` spariscono appena il legame scende sotto il ~2 %.
  Cioè: **forme platoniche e legame giusto, in questo modello, non stanno insieme**.
* **C3a-iv (un fit, poi una previsione).** Scegliamo α perché il ⁴He venga 28.3 MeV (interpolazione log-lineare fra
  0.8 e 0.95): **α\* = 0.894 (FIT, dichiarato)**. Con quello stesso α il modello *prevede*: deuterone **5.5 MeV**
  (vero 2.22: 2.5 volte troppo) e ³He **⟨3He⟩**. Esame: entrambi entro un fattore 2 → **⟨VERDETTO⟩**.
  (Corse dirette a α\* = 0.894 per controllare l'interpolazione: messe in coda ma **non fatte**, vedi debiti.)

## C3b — Il termine sestico (1 parametro: c₆)

Il termine c₆·B₀² (B₀ = densità barionica) è quello che in letteratura viene dal mesone ω: "costa" avere tanti
avvolgimenti nello stesso posto. Lo abbiamo **aggiunto** al modello con il pione massivo (unità AN tenute fisse,
c₆ = 50 e 200: nel B = 1 il termine pesa il 5 % e il 13 % dell'energia).

| c₆ | E₁ (MeV) | deuterone | ⁴He | forme |
|---|---|---|---|---|
| 0 (C1) | 855 | 68 MeV (4.0 %) | 300 MeV (8.8 %) | ciambella, cubo |
| 50 | 900 | **84 MeV (4.7 %)** | **352 MeV (9.8 %)** | ciambella, cubo |
| 200 | 987 | **106 MeV (5.4 %)** | **429 MeV (10.9 %)** | ciambella, cubo |

**C3b-i FAIL**: avevamo previsto che il legame scendesse (poco). Invece **sale**. Il perché, a posteriori, è semplice:
B₀² punisce chi è *concentrato*, e il nodo più concentrato di tutti è il nucleone singolo (tutto il suo avvolgimento
in una pallina di 0.66 fm); i nuclei più grandi sono gusci vuoti (ciambella, gabbie) con la densità sparsa, quindi il
termine li penalizza meno. Il sestico aiuta solo **vicino al limite "BPS"** (quando il sestico e un potenziale
*sostituiscono* i termini di Skyrme: modelli di Adam–Sánchez-Guillén–Wereszczyński), e anche lì Gillard–Harland–Speight
(modello sestico senza termine di Skyrme) trovano ancora 3.6–6 % di legame. Come pezzo mancante da solo: **no**.

## C4 — Il quadro completo e che cosa deve decidere l'agglomeratore

### C4a — Tabella riassuntiva

| strada | costanti/fit usati | parametri liberi in più | deuterone | ³He | ⁴He | forme |
|---|---|---|---|---|---|---|
| vero | — | — | **2.22** | **7.72** | **28.3** | d: nuvola larga |
| C0 classico senza massa | M_N, M_Δ (2 fit) | 0 | 72 | 178 | 310 | ciambella/tetraedro/cubo |
| C1 + massa del pione | m_π = 138 (costante) + M_N, M_Δ | 0 | 68 | 170 | 300 | **tenute** |
| C1 + C2a trottola rigida | idem | 0 | 160 | 344 | 610 | tenute |
| C1 + C2a + C2b punto zero (modi legati) | idem | 0 | 78 | — | — | tenute |
| C3a GHS α = 0.95 | E₁ := 938 MeV | 1 (α, valore di GHS) | 3.2 | ⟨3He95⟩ | 16.6 | **perse** (grumi) |
| C3a GHS α\* = 0.894 | E₁ := 938 MeV + ⁴He | 1 (α, **fit** sul ⁴He) | 5.5 (previsto) | ⟨3HeS⟩ (previsto) | 28.3 (fit) | perse |
| C3b sestico c₆ = 200 | come C1 | 1 (c₆) | 106 | — | 429 | tenute |

(MeV. "—" = non calcolato. Le righe C1–C3b sono di questo studio, sullo stesso reticolo; C0 da `NUCLEI_DA_NODI`.)

### C4b — Che cosa abbiamo capito

1. **Nessun pezzo "a costo zero" ripara il legame.** La massa del pione (costante della natura) è giusta da mettere
   ma cambia il legame di pochi %. La quantizzazione rigida lo **peggiora** (×2). Il punto zero dei modi legati lo
   riduce un po' (−82 MeV per il deuterone), ma resta 35 volte troppo e il pezzo più grosso (il continuo) è ignoto.
2. **Il legame giusto arriva solo pagando un parametro** (C3a), e allora i nuclei diventano **grumi di nucleoni quasi
   intatti** che si toccano appena — che è proprio come sono i nuclei veri (il deuterone vero è grande e lasco). Le
   forme platoniche (ciambella, cubo) erano il segno che il modello classico lega troppo: fondono i nucleoni in un
   oggetto nuovo, mentre in natura i nucleoni dentro il nucleo restano quasi se stessi.
3. **La trottola va attaccata al pezzo giusto.** Nel modello standard il nucleo è un corpo rigido unico e la trottola
   dei nucleoni (78 MeV l'uno) si perde → legame enorme. Se i nucleoni restano grumi quasi separati (C3a), ognuno tiene
   la sua trottola, e nel legame si cancella. (Non l'abbiamo calcolato qui: è il motivo per cui GHS confrontano con
   E₁ classica; è anche il cuore del risultato di Leese–Manton–Schroers, 6 MeV per il deuterone.)
4. Il sestico, aggiunto al modello standard, **non aiuta** (anzi lega di più).

### C4c — Raccomandazione per l'agglomeratore (livello "nuclei")

**Che cosa deve contenere la regola del livello nuclei:**

* **Le unità del livello sono i nucleoni (nodi B = 1), non il campo di palline.** Il nodo B = 1 si prende dal livello
  sotto (forma, massa, raggio, momento d'inerzia, già verificati in `NODO_3D` / `ZOOM_NODO`), con la massa del pione
  vera dentro (0 parametri in più).
* **Le etichette quantistiche vengono dalla topologia (0 parametri, robuste):** numero barionico intero, spin e
  isospin permessi (regole di Finkelstein–Rubinstein: deuterone J = 1, I = 0; ⁴He J = I = 0; ³He J = I = ½). Questo
  il modello lo fa giusto anche dove sbaglia le energie.
* **L'energia fra nucleoni va presa "leggermente legata":** o il modello GHS con **1 parametro fittato (α\* ≈ 0.89,
  sul ⁴He)**, dichiarato come debito, o — più onestamente — un'interazione efficace nucleone–nucleone da fittare sui
  nuclei leggeri. In tutti e due i casi il legame del ⁴He **non è una previsione**; deuterone e ³He sì, e lì il modello
  sbaglia ancora di un fattore ~2.5 (⟨3He commento⟩).
* **La trottola di spin/isospin si quantizza per nucleone**, non per nucleo rigido, quando il nucleo è "a grumi"; la
  quantizzazione a corpo rigido va usata solo per stati in cui i nucleoni sono davvero fusi (e lì sbaglia).
* **Non** usare le energie dei nuclei classici di Skyrme (ciambella/cubo) come energie di legame: sbagliano di 10–30
  volte. Si possono usare le loro **forme e simmetrie** come "stampo" per gli stati eccitati e per le regole di spin.

**Debiti che restano aperti** (sotto): la dinamica quantistica vera del deuterone (nuvola, non corpo rigido), le
correzioni del continuo (loop), la calibrazione GHS con F_π e m_π non fisici, e le corse dirette a α\*.

## Debiti e cose non fatte (dette chiare)

* **Il deuterone come nuvola quantistica non è fatto.** La strada giusta in letteratura (Leese–Manton–Schroers 1995:
  Schrödinger su uno spazio di 10 dimensioni di configurazioni di due nucleoni, legame 6.2 MeV, raggio 2.18 fm) qui
  non è rifatta. È il pezzo più importante che manca dal lato "0 parametri".
* **Il continuo delle vibrazioni (correzioni a un loop, "Casimir") non è calcolato**: serve rinormalizzare, cioè
  nuovi parametri. Per il B = 1 vale circa −200 MeV (Moussallam); per la ciambella 3D non è noto a nessuno. C2b conta
  solo i modi legati: è una stima parziale, il segno della correzione completa non è noto.
* **Reticolo h = 0.2.** Con il pione massivo il nodo è più piccolo e l'errore cresce: E₁ −1.2 % (m = 0.53), −1.7 %
  (m = 1) → **C1b FAIL** e **C1e FAIL**; momento d'inerzia −5.6 % → **C2a-i FAIL**. Nei legami (differenze sullo
  stesso reticolo) l'errore si cancella quasi tutto; nelle energie di rotazione resta (qualche MeV).
* **Calibrazione GHS non fisica** (F_π = 29 MeV, m_π = 303 MeV in GHS); qui il GHS è confrontato con E₁ classica
  := 938 MeV, senza trottola quantistica (vedi C4b punto 3: per nuclei a grumi è la scelta sensata, ma non è calcolata).
* **Corse non fatte (coda tagliata per il trasloco nel cloud, 23:54):** corse dirette a α\* = 0.894 (`ghsS_B1..B4`,
  modello già in `colla_units.json`; `ghsS_B3` era partito ed è stato fermato), partenze a coppia/tetraedro per
  α = 0.95 (`ghs95_B2pair`, `ghs95_B4tet`, pre-registrate come controllo della forma: con il rumore il campo si era già
  rotto in grumi da solo), B = 3 a α = 0 e 0.5 (`ghs0_B3`, `ghs50_B3`, solo per la tabella), sestico B = 3.
  **Come riprendere:** da `lab/`, `COLLA_THREADS=2 ./colla_drive.sh ghsS_B1 ghsS_B2 ghsS_B4 ghsS_B3 ghs95_B2pair ghs95_B4tet ghs0_B3 ghs50_B3`
  (riparte dai checkpoint `lab/out/colla_ckpt_*.pt` se ci sono), poi `python colla_summary.py` e `python colla_figure.py`.
* **Le forme a α = 0.5 vengono da partenze rumorose**: non abbiamo confrontato con un cubo "tenuto simmetrico" per
  vedere se il cubo è ancora un minimo locale (più alto). Il verdetto "il cubo si rompe" vuol dire "dal cubo leggermente
  disturbato il campo scende nei grumi".
* **Cambi dopo la pre-registrazione (tutti prima di vedere i numeri che toccano, salvo l'ultimo):** rumore di partenza
  2 % casuale → 3 % liscio e localizzato (21:26, per velocità: il rumore a scacchiera rallentava L-BFGS); tolleranza
  1e-8 → 1e-7 dalle 22:25 (an_B1..B4 e ghs95_B1 a 1e-8); tetto LOBPCG 150 → 60 iterazioni per il B = 2 (autovalori sotto
  soglia già fermi a 10⁻⁴); **segno di w per il ³He corretto dopo aver visto i numeri** (vedi C2a, con la ragione).
* **Log:** alcune righe dei log in `lab/out/colla_logs/` sono andate perse quando i log erano ancora tracciati da git
  (stash/ripristino automatici di altri agenti); i JSON dei risultati sono la fonte.
* **Il sestico** è stato provato solo sul modello con il pione massivo, unità AN tenute fisse (niente rifit), senza il
  suo contributo all'inerzia.
