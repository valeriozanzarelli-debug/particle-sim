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

<!-- C3a -->

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

<!-- C4 -->

## Debiti e cose non fatte (dette chiare)

<!-- DEBITI -->
