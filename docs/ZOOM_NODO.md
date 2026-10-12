# Zoom nel nodo: palline fini dove serve, agglomerati altrove

**L'idea (di Valerio).** "Dove zoomi vedi la verità (palline fini); altrove metti agglomerati precisi che introducono un
errore piccolo e controllato." Qui la applichiamo al **nodo di Skyrme** (B = 1) di `NODO_3D.md` e `PRESSIONE.md`.

**Il problema da risolvere.** In `PRESSIONE.md` (prova P4) il nodo fatto di palline tutte uguali (passo h = 0.2) sbagliava
la **pressione al centro del 57 %**. Per scendere al 5 % con palline tutte uguali servirebbe h ≈ 0.06, cioè circa
160³ = 4 milioni di palline: troppo per questo computer (4 core condivisi con altri due agenti).
**Obiettivo:** p(0) entro il 5 % e l'energia E/12π² entro lo 0.3 % dei valori esatti (E/12π² = 1.231446,
p(0) = 12.24 in unità di Skyrme), spendendo **non più di un reticolo uniforme 64³ = 262 144 palline**.

Pre-registrazione: `lab/zoomknot_prereg.py` → `lab/out/zoomknot_prereg.json`, **committata (bab46be) prima delle corse**,
insieme al motore `lab/zoomknot_core.py`. Corse: `lab/zoomknot_run.py` (+ `zoomknot_drive.sh`), apprendimento degli
agglomerati: `lab/zoomknot_learn.py`, verdetti: `lab/zoomknot_summary.py` → `lab/out/zoomknot_results.json`,
figura: `lab/zoomknot_figure.py` → `lab/out/zoomknot_figure.png`. Diagnosi post-hoc: `lab/zoomknot_ballcheck.py`.

## Come funziona, in parole semplici

Immagina il nodo come una cipolla. Al centro le frecce (le "palline-frecce", ognuna è un vettore a 4 componenti di
lunghezza 1) girano di 180° in pochissimo spazio: lì serve guardare da vicino. Fuori, il campo cambia lentamente
(va a zero come 1/r²): lì basta guardare da lontano.

Allora costruiamo **scatole una dentro l'altra**, come le matrioske:

| livello | passo h | la scatola arriva a (mezzo lato) | ruolo |
|---|---|---|---|
| 0 | 0.05 | 0.6 | palline fini: "la verità" nel cuore del nodo |
| 1 | 0.1 | 1.2 | palline medie |
| 2 | 0.2 | 2.4 | il reticolo "vecchio" di P4 |
| 3 → 7 | 0.4, 0.8, 1.6, 3.2, 6.4 | 4.8, 9.6, 19.2, 38.4, 89.6 | **agglomerati**: ogni pallina vale 2×2×2 di quelle del livello sotto |

Ogni livello ha **lo stesso numero di palline** (29³ = 24 389): è lo stesso disegno ripetuto a scale diverse, un
frattale. In totale 195 112 palline, meno di un 64³, ma la scatola arriva a **89.6** invece di 4.8.

Tre regole tengono insieme i livelli:
1. **Una sola verità per ogni punto.** Dove c'è un livello fine, il livello grosso copia i suoi valori (restrizione);
   l'energia di ogni pezzetto di spazio è contata **una volta sola** (pesi "a trapezio" che si sommano a 1).
2. **Il bordo del livello fine è "incollato" al livello grosso.** Le palline fini dell'ultima fascia (larga una
   pallina grossa) non sono libere: sono l'interpolazione cubica del livello grosso. Così il campo non può fare uno
   "scalino" fra fine e grosso (vedi sotto: senza questa regola tutto andava storto).
3. **Si rilassano tutti i livelli insieme** (L-BFGS, come nel P4), partendo dal profilo esatto.

### Gli agglomerati imparano la loro regola (`zoomknot_learn.py`)

Un agglomerato grosso, se usa la stessa formula dell'energia delle palline fini, **sbaglia**: in un campo che va come
1/r² le differenze finite grosse **sovrastimano** l'energia del 6–7 %. Allora facciamo come in `ZOOM.md`
(raggruppa, impara, controlla):

* **raggruppa:** prendiamo il nodo esatto ingrandito o rimpicciolito (fattori 0.7, 0.85, 1.2, 1.4: il fattore 1 è
  tenuto da parte per il controllo), lo mettiamo su palline 4 volte più fini, e lo **raggruppiamo** in agglomerati
  (due passi di blocco 2×2×2, pesi 1/4–1/2–1/4 per asse);
* **impara:** per ogni livello cerchiamo 2 numeri, k2 (per il termine σ) e k4 (per il termine di Skyrme), tali che
  l'energia degli agglomerati sia uguale a quella delle palline fini (minimi quadrati);
* **controlla:** sul nodo vero (fattore 1, mai visto in allenamento) confrontiamo l'energia di ogni guscio con
  l'integrale esatto.

| livello | h | guscio | k2 | k4 | errore della regola nuda | errore della regola imparata |
|---|---|---|---|---|---|---|
| 1 | 0.1 | 0.6–1.2 | 1.004 | 1.005 | −0.58 % | −0.14 % |
| 2 | 0.2 | 1.2–2.4 | 1.002 | 0.997 | +0.96 % | **+1.09 %** (peggio) |
| 3 | 0.4 | 2.4–4.8 | 0.965 | 0.888 | +5.9 % | +2.1 % |
| 4 | 0.8 | 4.8–9.6 | 0.942 | 0.809 | +6.6 % | +0.44 % |
| 5 | 1.6 | 9.6–19.2 | 0.941 | 0.805 | +6.6 % | +0.35 % |
| 6 | 3.2 | 19.2–38.4 | 0.941 | 0.805 | +6.6 % | +0.35 % |
| 7 | 6.4 | 38.4–89.6 | 0.943 | 0.805 | +6.4 % | +0.33 % |

La cosa bella: **da 0.8 in su la regola imparata è sempre la stessa** (k2 ≈ 0.941, k4 ≈ 0.805). È un **punto fisso**,
come il magnete di Ising a T_c in `ZOOM.md`: lontano dal nodo il campo è una legge di potenza (1/r²), che vista da più
lontano ha lo stesso aspetto, quindi la regola degli agglomerati non cambia più. 14 numeri imparati in tutto.
Il residuo +0.35 % è l'errore della "verità" usata per imparare (palline 4 volte più fini, non infinite).

## Tabellone

| prova | cosa | previsione | esito |
|---|---|---|---|
| Z0 | il codice riproduce il P4: reticolo uniforme h = 0.2, E/12π² entro 0.2 % da 1.20578 e nucleo sbagliato > 30 % | PASS | **PASS** (1.205475, −0.03 %; nucleo sbagliato del **68 %**, p(0) = 20.5) |
| Z1a | griglie annidate (fine 0.05 nel cuore, 0.2 fuori): p(0) e nucleo entro 5 % con ≤ 64³ palline | FAIL | **PASS** (R_f = 0.6: +1.0 %, nucleo 2.5 %; R_f = 1.0: +3.4 %, nucleo 3.5 %) |
| Z1b | l'errore del nucleo scende allargando la scatola fine (0.6 → 1.0 → 1.4) | PASS | **FAIL** (2.5 % → 3.5 % → 8.6 %: sale invece di scendere) |
| Z1c | migliore Z1 (R_f = 0.6): E entro 0.3 % | FAIL | **FAIL** (−1.01 %; per curiosità R_f = 1.0 dà +0.14 %, ma per compensazione) |
| Z2a | agglomerati con regola imparata: p(0) e nucleo entro 5 % | PASS (incerto) | **PASS** (p(0) = 12.23, **−0.08 %**; nucleo 0.3 %) |
| Z2b | E entro 2 % (obiettivo) / entro 0.3 % (traguardo) | PASS / 50-50 | **PASS / PASS** (E/12π² = 1.229512, **−0.16 %**) |
| Z2c | momento d'inerzia Λ entro 3 % | PASS | **PASS** (107.8 contro 106.8, +0.9 %) |
| Z2d | imparare aiuta: (i) ≥ 6 livelli su 7 migliori nel controllo; (ii) E migliore della regola nuda | PASS | **PASS** (6/7; E −0.16 % contro −0.23 %) |
| Z2e | errore statico totale degli agglomerati < 0.3 % di E | PASS | **PASS** (+0.21 %; regola nuda +0.13 %, vedi sotto) |
| Z3 | bilancio degli errori, "principio dello zoom" | — | vedi tabella sotto: **verificato** |
| Z4 | nodo in moto (v = 0.3) sulle griglie annidate | SKIPPED | **SKIPPED** (budget) |

## Il principio dello zoom, in una tabella (Z3)

Errori rispetto all'ODE esatta. "fine" = palline del livello più fine; "agglomerati" = tutte le altre (comprese le
palline fantasma e quelle coperte: è un conto per eccesso del lavoro).

| configurazione | palline totali | fini | grosse/agglomerati | errore E | errore p(0) | nucleo (max) | errore Λ |
|---|---|---|---|---|---|---|---|
| U0: tutte uguali h = 0.2, scatola 4.8 | 117 649 | — | 117 649 | −2.11 % | **+67.7 %** | 67.7 % | −32 % |
| Z1: fine 0.05 fino a 0.6 + 0.2 fuori | 142 038 | 24 389 | 117 649 | −1.01 % | +0.96 % | 2.5 % | −24 % |
| Z1: fine 0.05 fino a 1.0 + 0.2 fuori | 208 774 | 91 125 | 117 649 | +0.14 % | +3.40 % | 3.5 % | −22 % |
| Z1: fine 0.05 fino a 1.4 + 0.2 fuori | 344 630 (oltre budget) | 226 981 | 117 649 | +0.42 % | +8.64 % | 8.6 % | −23 % |
| Z1: fine 1/15 fino a 1.0 + 0.2 fuori | 160 524 | 42 875 | 117 649 | +0.12 % | +6.48 % | 6.4 % | −22 % |
| Z2: 8 livelli, regola nuda | 195 112 | 24 389 | 170 723 | −0.23 % | +1.01 % | 1.2 % | −2.1 % |
| **Z2: 8 livelli, regola imparata** | **195 112** | 24 389 | 170 723 | **−0.16 %** | **−0.08 %** | **0.3 %** | **+0.9 %** |

Come si legge:
* **Tutto grosso → 68 %** di errore sul centro. **Palline fini solo nel cuore → 1–3 %** (con h = 0.05; con h = 1/15 il 6.5 %: l'errore scala circa come h², rapporto 1.9 contro 1.8). È lo zoom: la verità sta dove
  si guarda da vicino. Costo in più: appena 24 389 palline fini (un quinto del reticolo grosso).
* Ma con la sola scatola fine (Z1) **E e Λ restano sbagliati** per un motivo che non ha niente a che fare col cuore:
  la scatola finisce a 4.8 e il muro (campo forzato a zero) schiaccia il nodo. Lo mostra il conto esatto in una sfera
  di raggio 4.8 (`zoomknot_ballcheck.py`, post-hoc): il solo muro dà E +1.0 %, p(0) +15 %, Λ −28 %. Quindi i buoni
  numeri di Z1 sono in parte una **compensazione** fra il muro (che spinge p(0) in su) e il reticolo grosso fuori
  (che lo tira giù). Per questo Z1b fallisce: allargando la scatola fine si toglie errore del reticolo grosso, ma
  resta quello del muro, e l'equilibrio si sposta.
* **Gli agglomerati (Z2) spostano il muro a 89.6** con sole 5 × 24 389 palline in più, e tutto torna: E −0.16 %,
  p(0) −0.08 %, Λ +0.9 %. È la seconda metà dell'idea di Valerio: lontano non servono palline fini, servono
  **agglomerati con la regola giusta**.
* **Imparare la regola conta soprattutto per Λ** (−2.1 % → +0.9 %) e un po' per E (−0.23 % → −0.16 %) e p(0)
  (+1.0 % → −0.08 %). Λ dipende molto dalla coda lontana (∫ r² sin²f), che è proprio dove la regola nuda sbaglia del 6 %.

Altri numeri di Z2 (regola imparata): B = 0.9966, zero della pressione r₀ = 1.078 (ODE 1.076), virial
(E₄−E₂)/E = −0.5 % (esatto 0), von Laue ∫p/∫|p| = −0.012 (P4 sul reticolo uniforme: 0.21). Rilassamento convergito
(12 passi esterni, 8.5 minuti su 1 core).

## Il fallimento che ci ha insegnato di più: lo scalino al confine

La **prima versione** del motore (quella committata con la pre-registrazione) lasciava libere anche le palline fini
sulla faccia della scatola fine; le palline "fantasma" appena fuori erano l'interpolazione del livello grosso.
Risultato (`lab/out/zoomknot_v1bug/`): Z1 con R_f = 1.0 dava E −1.2 %, p(0) −6.3 %, e **B = 0.967** invece di 1;
Z2 dava E −1.6 % e p(0) −7.4 %. Guardando dentro: le ultime due file di palline fini si scostavano dal profilo esatto
di 0.04–0.07, mentre le fantasma subito fuori di solo 0.003. Cioè il campo faceva **uno scalino** sul confine, e
quello scalino costava quasi zero: il "ponte" fra la faccia e le fantasma pesa solo 1/4 nell'energia, e il livello
grosso non lo vede (vede solo una pallina fine ogni 2 o 4). Una crepa nel vetro: il nodo ci si "infilava" per
abbassare l'energia in modo finto.

**Correzione** (prima delle corse buone, non pre-registrata — la dichiariamo): l'ultima fascia del livello fine (larga
una pallina grossa) non è più libera ma è **l'interpolazione cubica del livello grosso** (regola 2 sopra). Le corse
nella tabella sono tutte fatte con la versione corretta; U0 (un solo livello) non cambia. Lezione per lo zoom in
generale: **il punto delicato non è l'agglomerato, è la cucitura** fra agglomerato e palline fini.

## Debiti e limiti, detti chiari

* **Correzione del motore dopo la pre-registrazione** (vedi sopra). I criteri e le configurazioni sono quelli
  pre-registrati; le corse con il motore difettoso sono conservate in `lab/out/zoomknot_v1bug/`.
* **Una corsa pilota** (Z1, R_f = 1.0, motore v1) era stata fatta prima della pre-registrazione e dichiarata lì.
* **Z2e passa ma la regola nuda fa meglio nel totale statico** (+0.13 % contro +0.21 %): la regola nuda sbaglia il
  livello 1 in giù e gli altri in su, e gli errori si compensano per caso. Livello per livello la regola imparata è
  migliore 6 volte su 7. Il livello 2 (h = 0.2, guscio 1.2–2.4) è quello dove imparare **non** aiuta: lì il nodo non
  è ancora una legge di potenza e un solo numero per termine non basta a descrivere configurazioni di taglie diverse.
* **Z1b FAIL**: l'errore nel cuore non scende allargando la scatola fine, perché domina il muro a 4.8 (diagnosi
  post-hoc, `zoomknot_ballcheck.py`).
* **Λ**: il termine cinetico non è rinormalizzato (k = 1 per l'energia cinetica); il buon accordo di Z2 (+0.9 %)
  contiene una compensazione fra la coda tagliata a 89.6 (−1.5 % nell'ODE) e il resto.
* **Z4 (nodo in moto) non fatto**: le dinamiche sulle griglie annidate richiedono passi di tempo fissati dalla
  pallina più fine; fuori budget.
* Costo contato per eccesso (tutte le palline, comprese fantasma e coperte). Tempo reale: 8–12 minuti per corsa su
  1 core, contro i ~2.5 minuti (2 core) del reticolo uniforme h = 0.2 che sbaglia del 68 %.

## Quanti parametri

* Fisica: **0 parametri liberi** (unità di Skyrme, come in P1–P4).
* Metodo: la geometria delle scatole (8 livelli × 3 numeri: passo, mezzo lato, palline per lato) scelta prima delle
  corse; **14 numeri imparati** (k2, k4 per i 7 livelli non fini) dalle palline fini, mai riaggiustati sui risultati.

## In una frase

Con le palline fini **solo nel cuore** e **agglomerati che imparano la loro regola** fuori, il nodo di Skyrme ha la
pressione al centro giusta entro lo **0.1 %** (prima: 68 %) e l'energia entro lo **0.16 %**, con meno palline di un
reticolo 64³: lo zoom funziona, a patto di **cucire bene** i bordi fra fine e grosso.
