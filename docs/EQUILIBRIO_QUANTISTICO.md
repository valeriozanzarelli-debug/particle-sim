# Equilibrio quantistico: la regola di Born si raggiunge mescolando?

**Idea di Valerio.** La regola di Born (i surfisti si distribuiscono come |ψ|²) non sarebbe una legge,
ma un **equilibrio**: come lo zucchero nel caffè, che si spande da solo se qualcuno gira il cucchiaino.
Il cucchiaino sono i **vortici**: i punti dove l'onda si annulla (ψ = 0) e la fase gira attorno.
Lì vicino i surfisti girano velocissimi e si mescolano (Valentini 1991, "teorema H sub-quantistico";
Efthymiopoulos–Contopoulos: il caos nasce dai vortici ai nodi; Towler–Russell–Valentini 2012:
il rilassamento è quasi esponenziale e più veloce con più modi).

**Il ricordo che ci ha portato qui.** La doppia fenditura dava 0,89 invece di 0,995 quando i surfisti
partivano dalla distribuzione sbagliata. Lo abbiamo risolto *facendoli partire giusti*. Domanda di oggi:
se partono sbagliati, si aggiustano da soli?

## Come si legge

* **H̄(t)** = "quanto la folla dei surfisti è diversa da |ψ|²", misurata a celle (griglia 16×16 salvo
  dove detto). H̄ = 0 vuol dire equilibrio. Con 20 000 surfisti c'è un **pavimento** di rumore
  F = (celle − 1)/(2N) ≈ 0,0064: sotto non si scende.
* **τ** = tempo in cui (H̄ − F) scende di un fattore e. Velocità di mescolamento = 1/τ.
* **Vortici** = punti dove la fase di ψ fa un giro intero (contati su griglia 128², media nel tempo).
* Unità: ħ = m = 1, scatola [0, π]², ψ = somma di M modi sin(nx)sin(my) con fasi a caso.

**Pre-registrazione.** I criteri sono stati scritti e committati *prima* delle corse:
`lab/relax_prereg.py` → `lab/out/relax_prereg.json` (commit 782880f, 2026-10-09 18:12).
Deviazioni dichiarate: (1) durata tagliata da 8π a 4π per stare nei tempi (una prova di 1 unità di
tempo, fatta per misurare la velocità del codice, mi ha già mostrato un pezzetto di H̄ per M = 25);
(2) nel test 3a ho aggiunto *dopo* due partenze più sbagliate, perché quella pre-registrata era troppo
poco sbagliata (vedi sotto). Quelle aggiunte sono marcate "post-hoc".

## Risultati

Figura: `lab/out/relax_figure.png`. Numeri: `lab/out/relax_results.json`.

### Test 1 — Il rilassamento noto si riproduce? **Sì, quasi tutto.**

| M (modi) | vortici medi | partenza | H̄(0) | H̄(4π) | τ (16×16) |
|---|---|---|---|---|---|
| 1 | 0 | \|ψ\|⁴ | 0,136 | 0,136 | ∞ (fermo) |
| 1 | 0 | uniforme | 1,154 | 1,154 | ∞ |
| 4 | 1,0 | \|ψ\|⁴ | 0,148 | 0,111 | > 4π |
| 4 | 1,0 | uniforme | 0,861 | 0,790 | > 4π |
| 16 | 5,5 | \|ψ\|⁴ | 0,233 | 0,022 | 4,3 |
| 16 | 5,5 | uniforme | 0,789 | 0,113 | 7,2 |
| 25 | 9,9 | \|ψ\|⁴ | 0,198 | 0,008 | 2,4 |
| 25 | 9,9 | uniforme | 0,609 | 0,024 | 3,2 |
| 16 controllo | 5,5 | \|ψ\|² | 0,0070 | 0,0060 | — (resta sul pavimento) |

| criterio pre-registrato | esito |
|---|---|
| C1 M = 1 non rilassa (H̄ finale/iniziale > 0,9) | **PASS** (1,00 e 1,00: senza vortici i surfisti sono fermi) |
| C2 partendo da \|ψ\|² si resta in equilibrio (max H̄ ≤ 3F) | **PASS** (0,0078 ≤ 0,019): il motore non "inventa" disordine |
| C3 M = 16 e 25 scendono sotto il 10% in 4π | **FAIL** per un caso: M=16 uniforme resta al 13,7% (gli altri: 6,7%, 0,7%, 2,9%). Con la pendenza misurata arriverebbe sotto il 10% poco dopo, ma non l'ho corso |
| C4 più modi → più vortici → più veloce (stesso ordine) | **PASS**: vortici 1 / 5,5 / 9,9; velocità ≈0 / 0,23 / 0,42 (\|ψ\|⁴) e ≈0 / 0,14 / 0,31 (uniforme) |
| C5 decadimento circa esponenziale (R² ≥ 0,9) | **PASS** (R² 0,975–0,994) |

In parole: **con molti vortici la folla si mescola e va a |ψ|² da sola, in modo quasi esponenziale;
senza vortici non succede nulla**. Con un solo vortice (M = 4) il mescolamento è lentissimo e parziale.
Raddoppiare i vortici (5,5 → 9,9) quasi raddoppia la velocità. Solo 3 punti, però: è una tendenza, non una legge.

### Test 2 — Le scale grandi restano indietro? **Metà sì, metà no.**

* **2a (stessa corsa, celle grandi contro celle piccole): FAIL, anzi il contrario.** A celle 4×4 il
  disordine sparisce *prima* (τ = 1,1) che a celle 16×16 (τ = 4,3) per M = 16; per M = 25 1,95 contro
  2,39. Le celle grosse non vedono i filamenti fini, quindi vedono l'equilibrio prima.
* **2b (perturbazione di lunghezza d'onda lunga contro corta, celle 32×32): PASS.** Una "gobba" larga
  (cos 2x) si cancella con τ = 6,9; una fitta (cos 12x) con τ = 2,8: la lunga è **2,5 volte più lenta**.
* **Onestà.** Nell'argomento di Valentini sul fondo cosmico il ritardo nasce da modi che escono
  dall'orizzonte (e quindi smettono del tutto di evolvere), non dalla lunghezza d'onda dentro una
  scatola. Qui ho solo un'analogia debole: "le disuguaglianze larghe si spianano più lentamente di quelle
  fini". Non dice niente sul cosmo.

### Test 3 — La nostra doppia fenditura si aggiusta da sola? **No, se prima delle fenditure non ci sono vortici.**

**3a, una dimensione trasversa (come il nostro esperimento classico).** Prima delle fenditure c'è un
fascio gaussiano: **zero vortici**. In 1D poi le traiettorie di Bohm **non si possono incrociare**:
l'ordine dei surfisti è congelato, quindi l'errore di partenza arriva intatto sullo schermo.

| schermo z | 0,5 | 2 | 5 | 10 | 20 |
|---|---|---|---|---|---|
| campionatore ideale (media ± σ) | 0,9985±0,0006 | 0,9988±0,0004 | 0,9992±0,0004 | 0,9993±0,0004 | 0,9991±0,0005 |
| partenza \|ψ\|² | 0,9939 | 0,9990 | 0,9993 | 0,9998 | 0,9996 |
| partenza uniforme (pre-registrata) | 0,9921 | 0,9958 | 0,9984 | 0,9982 | 0,9981 |
| post-hoc: 2/3 a sinistra | 0,927 | 0,940 | 0,943 | 0,943 | 0,939 |
| distanza KS, 2/3 a sinistra | 0,163 | 0,157 | 0,159 | 0,161 | 0,162 |

* Criterio pre-registrato (uniforme sotto la banda a ogni z): **FAIL**, ma per un motivo noioso: con
  bordi morbidi e fascio largo "uniforme nelle fenditure" è quasi uguale a |ψ|², quindi l'errore è troppo
  piccolo per vederlo con la correlazione a bin larghi (a z ≥ 5 cade dentro la banda).
* Post-hoc, errore grande (2/3 dei surfisti nella fenditura di sinistra): la distanza KS resta
  **0,16 da z = 0,5 a z = 20**, la correlazione resta 0,94. **Nessun rilassamento**: le frange restano
  storte per sempre. Questo è il punto del vecchio 0,89: in quella geometria l'errore non si sarebbe
  mai corretto da solo.

**3b, due dimensioni trasverse con una "guida multimodo" prima delle fenditure.** I surfisti partono
uniformi (sbagliati) in una guida quadrata con M modi, ci viaggiano per una lunghezza Z0, poi passano
la doppia fenditura e arrivano allo schermo (z = 1,5 dopo le fenditure).

| guida | Z0 | partenza sbagliata | partenza \|ψ\|² | campionatore ideale |
|---|---|---|---|---|
| M = 16 (5,5 vortici) | 0 | **0,864** | 0,998 | 0,9974±0,0010 |
| M = 16 | π | **0,931** | 0,998 | 0,9975±0,0009 |
| M = 16 | 4π | **0,980** | 0,995 | 0,9976±0,0007 |
| M = 1 (0 vortici) | 0, π, 4π | **0,970, 0,970, 0,970** | 0,999 | 0,9989±0,0005 |

* La correlazione **sale da sola** 0,86 → 0,93 → 0,98 quando c'è mescolamento a monte; con M = 1
  resta **ferma** a 0,970 qualunque sia la lunghezza.
* Criterio pre-registrato (M = 16, Z0 = 4π dentro 3σ dal campionatore ideale): **FAIL**. 0,980 non è
  0,995: in 4π il mescolamento ha tolto ~86% del disordine (coerente con H̄ M=16 uniforme, 13,7% residuo),
  ma non tutto. Servirebbe una guida più lunga o con più modi (M = 25 è ~2× più veloce).

### Test 4 (CHSH con coppie entangled) — **saltato** per tempo.

## Cosa vuol dire, in breve

1. **Sì**: nel motore surfista la regola di Born *si comporta come un equilibrio*. Partendo sbagliati,
   con vortici, si arriva a |ψ|² in modo circa esponenziale; più vortici, più veloce. Partendo giusti si
   resta giusti.
2. **Ma serve il cucchiaino.** Senza vortici (un solo modo, un fascio gaussiano, qualsiasi cosa in 1D)
   l'errore resta congelato per sempre. La nostra doppia fenditura storica è proprio così: quindi far
   partire i surfisti giusti non era un trucco, era l'unica strada in quella geometria.
3. Questo non dimostra che la natura funzioni così: dimostra che **la nostra regola è compatibile** con
   l'idea, e riproduce i numeri noti in letteratura.

## Debiti e limiti

* Durata 4π invece di 8π (tempo). C3 e 3b falliscono di poco forse solo per questo: non verificato.
* Solo 3 valori di M con vortici: la relazione "velocità ∝ numero di vortici" è un indizio.
* Integrazione: RK4 con sotto-passi vicino ai nodi (passo massimo 0,02 in spazio); il controllo C2 dice
  che non crea disordine finto, ma vicino ai vortici gli errori numerici aiutano un po' il mescolamento.
* 3b: lo schermo è a z = 1,5 (in campo vicino-intermedio) e la velocità a valle è interpolata su griglia.
  Lo script 3b dura 319 s, più dei 3 minuti previsti.
* 2a/2b sono analogie, non il meccanismo cosmologico di Valentini.
* La geometria è una scatola continua con modi analitici, **non** la regola su reticolo (zig-zag):
  il mescolamento sul reticolo non è ancora provato.

## File

* `lab/relax_prereg.py` → `lab/out/relax_prereg.json` (pre-registrazione)
* `lab/relax_core.py` (scatola, campi, RK4 con sotto-passi, H̄, conta vortici)
* `lab/relax_box.py` (test 1 e 2; ~15–180 s per corsa) → `lab/out/relax_box_*.npz`
* `lab/relax_slit.py 1d|2d` (test 3; 27 s / 319 s) → `lab/out/relax_slit_{1d,2d}.json`
* `lab/relax_figure.py` → `lab/out/relax_figure.png`, `lab/out/relax_results.json`
