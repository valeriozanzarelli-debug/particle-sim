# Pressione nel nodo: quanto "spinge" il vuoto dentro una particella

**L'idea (di Valerio).** La massa sarebbe "quanto il vuoto preme la bolla" dentro una particella.
C'è una versione **misurata** di questa idea: la **distribuzione di pressione dentro il protone**
(Burkert, Elouadrhiri, Girod, *Nature* 557, 396 (2018), dati JLab). Lì dentro la pressione è **positiva al centro**
(spinge verso fuori, ~10³⁵ Pa, più che nel cuore di una stella di neutroni) e **negativa sul bordo** (tira verso dentro,
"confina"). Le due parti si cancellano esattamente: ∫ r² p(r) dr = 0 (condizione di **von Laue**: un oggetto fermo e
stabile non può avere pressione netta, altrimenti esploderebbe o collasserebbe). Il numero che riassume tutto si chiama
**D-term** (D < 0 per un oggetto stabile).

Qui prendiamo il **nodo di Skyrme** già costruito in `NODO_3D.md` (il "riccio" di palline-frecce, B = 1) e ne
calcoliamo pressione p(r), forza di taglio s(r), densità di energia T₀₀(r) e D. Poi li mettiamo accanto al protone.

Pre-registrazione: `lab/press_prereg.py` → `lab/out/press_prereg.json`, **committata (4331d44) prima delle corse**,
con le previsioni scritte prima (anche quelle che sapevamo sarebbero FAIL).
Codice: `lab/press_ode.py` (P1, P2, P3, P5), `lab/press_lattice.py` (P4), `lab/press_lattice_posthoc.py` (diagnosi
**post-hoc**, non pre-registrata), `lab/press_figure.py`. Numeri: `lab/out/press_ode.json`, `press_lattice.json`,
`press_posthoc_h*_L*.json`. Figura: `lab/out/press_figure.png`.

## La fisica in tre righe

Il campo ha due termini di energia: e₂ = |∂φ|² (il "termine σ", quello che vuole le frecce allineate) ed
e₄ = termine di Skyrme (quartico). Dal tensore energia-impulso viene, punto per punto:

* **p = (e₄ − e₂)/3.** Il termine di Skyrme dà pressione **positiva** p₄ = +e₄/3, esattamente come un gas di luce
  (radiazione: p = ρ/3); il termine σ dà pressione **negativa** p₂ = −e₂/3 (tira, come un elastico).
* Al centro vince e₄ (spinge), fuori vince e₂ (tira): **una sola** inversione di segno, a r₀.
* Il nodo è in equilibrio quando le due spinte si bilanciano in totale: ∫p d³r = 0 ⇔ E₂ = E₄ (teorema di **Derrick**).

Formule (unità di Skyrme): s = 2(f′² − sin²f/r²)(1 + sin²f/r²), T₀₀ = e₂ + e₄,
D = M ∫d³r r² p = −(4/15) M ∫d³r r² s (Polyakov–Schweitzer 2018). Sono le stesse di Cebulla, Goeke, Ossmann,
Schweitzer, NPA 794 (2007) 87: lo verifichiamo riproducendo la loro Tabella 2 (P1e).

## Tabellone

| prova | cosa | esito pre-registrato |
|---|---|---|
| P1a | ∫T₀₀ d³r = M (massa del risolutore) entro 10⁻⁶ | **PASS** (3·10⁻¹¹) |
| P1b | von Laue: \|∫r²p dr\| / ∫r²\|p\| dr < 10⁻³ | **PASS** (6·10⁻¹⁰) |
| P1c | p(0) > 0, un solo cambio di segno, coda negativa fino a r = 50 | **PASS** (r₀ = 1.076 unità Skyrme) |
| P1d | D da p e D da s uguali entro 2 % | **PASS** (2·10⁻⁹) |
| P1e | ai parametri di Cebulla (e = 4.628, F_π = 131.3 MeV) riprodurre T₀₀(0), p(0), r₀, d₁ entro 2 % | **PASS** (tutti entro 0.2 %: 1.577 vs 1.58 GeV/fm³; 0.3168 vs 0.317; 0.699 vs 0.70 fm; −6.603 vs −6.6) |
| P2a | r₀ entro 25 % da 0.6 fm | **PASS** (0.605 fm; tripolo JLab 0.586 fm) |
| P2b1 | p(0) entro un fattore 2 da 10³⁵ Pa | **PASS** (6.5·10³⁴ Pa) |
| P2b2 | p(0) entro un fattore 2 dal p(0) del tripolo JLab | **FAIL** (6.5·10³⁴ vs 1.9·10³⁵ Pa: rapporto 0.34) — *previsto FAIL* |
| P2c | D nella banda della letteratura Skyrme [−5.28, −3.58] | **FAIL** (D = −2.77) — *previsto FAIL* (D ∝ 1/e⁴, e di ANW più grande) |
| P3a | Derrick E₂ = E₄ entro 10⁻⁴ | **PASS** (2·10⁻¹⁰) |
| P3b | ∫T^μ_μ = 2E₂ = M; traccia del termine di Skyrme = 0 punto per punto | **PASS** (2·10⁻¹⁰; 7·10⁻¹⁵) |
| P3c | massa di Tolman ∫(T₀₀ + 3p) = M | **PASS** (2·10⁻¹⁰) |
| P4a | reticolo (h = 0.2, L = 48): nel nucleo r < r₀/2, \|p_ret − p_ODE\| / p(0) < 5 % | **FAIL** (57 %) |
| P5a | G M²/R ÷ Mc² (M = 939 MeV, R = 0.6 fm) in [10⁻⁴⁰, 10⁻³⁸] | **PASS** (2.1·10⁻³⁹) |
| P5b | luce: massa attiva di Tolman ρ + 3p/c² = 2ρ | **PASS** (fattore 2 esatto) |

13 PASS, 3 FAIL (due previsti prima di correre, uno no: P4).

## P1 — Il nodo esatto: la pressione c'è, ed è fatta come quella del protone

Dal profilo esatto f(r) (collocazione `radial_bvp`, coda analitica f = C/r²) viene la figura in alto a destra:
la pressione è positiva fino a r₀, poi negativa e piccola fino all'infinito (coda −2C²/r⁶); la forza di taglio s è
sempre positiva. Tutti i controlli interni passano a livello di 10⁻⁹–10⁻¹¹, e soprattutto **i numeri di Cebulla et al.
(2007) escono uguali entro 0.2 %**: il nostro calcolo e il loro sono lo stesso calcolo. In unità di Skyrme (zero
parametri): p(0) = 12.24, T₀₀(0) = 60.9, r₀ = 1.076, D·e⁴ = −2423.

## P2 — In unità vere: accanto al protone

Per passare a fm e Pa servono 2 numeri: e = 5.44, F_π = 129 MeV (**FIT** di Adkins–Nappi–Witten su massa del nucleone
e della Δ, non una previsione). Allora 1 unità di lunghezza = 0.562 fm, 1 unità di energia = 5.93 MeV, e la massa
classica del nodo è 865 MeV.

| | nodo (ANW) | protone (JLab, tripolo) |
|---|---|---|
| zero della pressione r₀ | **0.605 fm** | 0.586 fm ("vicino a 0.6 fm") |
| picco di r²p (spinta) | 0.31 fm | 0.25 fm |
| minimo di r²p (presa) | 0.87 fm | 0.92 fm |
| p(0) | **6.5·10³⁴ Pa** (0.41 GeV/fm³) | 1.9·10³⁵ Pa (1.18 GeV/fm³) |
| D | **−2.77** | −1.63 (solo quark, μ² = 1.4 GeV²) |

Come si legge (figura in alto a sinistra): **la forma è quasi la stessa** — spinta al centro, presa sul bordo, zero
nello stesso punto entro il 3 %. Il nodo però è **più "morbido" al centro** (p(0) tre volte più piccolo del tripolo) e
spinge un po' più in fuori. Il "10³⁵ Pa" che si legge nei giornali è un numero arrotondato: contro di lui il nodo passa
(fattore 0.65); contro la formula vera del fit JLab no (fattore 0.34). Lo avevamo previsto prima di correre.

Sul D: la letteratura Skyrme (Cebulla, e = 4.628) dà D da −5.3 a −3.6; con e = 5.44 il nodo dà −2.77, **fuori dalla
banda** (FAIL, previsto: D scala come 1/e⁴). Curiosamente è più vicino ai dati: JLab −1.63 (quark soli; il dato è
discusso — Kumerički, *Nature* 570 (2019) E1, dice che è compatibile anche con zero) e QCD su reticolo
−3.87 ± 0.97 (quark + gluoni, Hackett–Pefkou–Shanahan, PRL 132 (2024) 251904): siamo a 1.1 σ da quest'ultimo.
Non è un successo da vendere: abbiamo solo spostato e, e D si muove come e⁻⁴.

**Debito dichiarato.** Il parametro M² = 1.02 GeV² del tripolo JLab è citato a memoria (il testo primario non era
raggiungibile in questa sessione). Controllo di coerenza fatto prima delle corse: con quel valore lo zero cade a
3/M = 0.586 fm e il picco di r²p a (3−√3)/M = 0.248 fm, come scritto nel testo verificato (zero "vicino a 0.6 fm",
picco vicino a 0.25 fm).

## P3 — "Massa = pressione"? Cosa dice davvero il conto

Qui sta la risposta onesta all'idea di Valerio, in tre fatti esatti (nessun fit):

1. **La pressione totale è zero.** ∫p d³r = 0 (von Laue). Quindi la massa **non** è "la pressione" sommata: sommata fa
   zero. La massa è il prezzo di un **equilibrio** fra chi spinge e chi tira.
2. **Metà e metà.** Nel nodo E₂ = E₄ = M/2 esattamente (Derrick). Il termine che spinge (Skyrme, p = +e₄/3, si comporta
   come luce) e quello che tira (σ, p = −e₂/3) portano ciascuno **il 50 %** della massa: 432 MeV + 432 MeV.
3. **La massa è la traccia.** T^μ_μ = T₀₀ − 3p = 2e₂, e il suo integrale è **esattamente M**. Il termine di Skyrme ha
   traccia zero (come la luce, come i gluoni classici): da solo non saprebbe "fissare una scala". La massa viene tutta
   dal termine che rompe l'invarianza di scala. Nel protone succede la stessa cosa in forma quantistica: M = ⟨traccia⟩,
   e la traccia è fatta di **anomalia di traccia** dei gluoni + massa dei quark.

Confronto **solo qualitativo** con la decomposizione della massa del protone (χQCD, Yang et al., PRL 121 (2018) 212001,
reticolo, MS-bar 2 GeV; la decomposizione è di Ji, PRL 74 (1995) 1071): massa dei quark ~9 %, energia dei quark ~32 %,
energia del campo dei gluoni ~36 %, anomalia di traccia ~23 %. La lezione comune: **la massa è quasi tutta energia di
campo in equilibrio, non "roba"**. Ma le caselle non si corrispondono una a una (il nodo non ha quark né gluoni):
nessun numero del nodo va confrontato con un numero del protone qui.

P3c: la "massa di Tolman" ∫(T₀₀ + 3p) — quella che sente la gravità — è uguale a M **proprio perché** ∫p = 0. La
pressione dentro c'è, enorme, ma per la gravità si cancella.

## P4 — Il nodo di palline (reticolo): FAIL nel nucleo

Abbiamo rilassato il nodo sul reticolo (h = 0.2, L = 48, metà scatola 4.7, L-BFGS come K1: E/12π² = 1.205782, la
stessa di `NODO_3D.md`) e misurato p(r) con derivate centrali del 4° ordine, mediate su gusci attorno al centro.

* **Fuori dal nucleo va bene**: lo zero della pressione è a r₀ = 1.049 (ODE 1.076, −2.5 %), la parte negativa segue l'ODE.
* **Nel nucleo no**: nel primo guscio (r ≈ 0.17, cioè meno di una pallina dal centro) p = 18.0 contro 11.0 dell'ODE.
  Scarto massimo 57 % di p(0) contro il 5 % richiesto: **FAIL**.
* Anche il bilancio globale soffre: virial (E₄−E₂)/E = +9 %, von Laue sul reticolo 0.21 (invece di 0) — il reticolo e
  la scatola con bordo fisso tagliano la coda 1/r⁶ e alterano l'equilibrio. Dichiarato come "non giudicato" nella
  pre-registrazione.

**Diagnosi post-hoc (non pre-registrata)**, `press_lattice_posthoc.py`, stessa scatola 4.7:

| h | L | pendenza al centro a (ODE 2.0084) | errore nel nucleo, derivate centrali | errore, differenze dell'energia | virial |
|---|---|---|---|---|---|
| 0.2  | 48 | 2.098 (+4.5 %) | 57 % | 48 % | +7.4 % |
| 0.15 | 64 | 2.014 (+0.3 %) | 33 % | 29 % | +4.6 % |

L'errore scende come h² (rapporto 1.72; h² darebbe 1.78): è **discretizzazione**, non fisica sbagliata. La pressione
al centro è la quantità più cattiva da misurare: p(0) = a²(a²−1) va come la **quarta potenza** della pendenza, e vicino
al centro del riccio le frecce girano di 180° in una o due palline. Per arrivare al 5 % servirebbe h ≈ 0.06
(≈ 160³ palline), fuori dal budget di questa sessione (4 core condivisi). Nota: a h = 0.15 il rilassamento si è
fermato per tempo (E variava di 10⁻⁷ relativo), non per il criterio di convergenza.

## P5 — Riquadro di onestà: la gravità non c'entra

* Energia gravitazionale del protone su se stesso: G M²/R con M = 939 MeV, R = 0.6 fm, divisa per Mc²:
  **2·10⁻³⁹**. La gravità è 39 ordini di grandezza troppo debole per fare la massa di una particella. La "pressione"
  di cui parliamo qui è quella del **campo forte** (pioni nel nodo; quark e gluoni nel protone), non dello spazio-tempo.
* Curiosità vera (Tolman): un gas di luce, che ha p = ρc²/3, gravita come ρ + 3p/c² = **2ρ**, il doppio della sua
  energia. Il termine di Skyrme dentro il nodo è "luce" in questo senso (peso 2), il termine σ ha peso 0; sommati
  sull'intero nodo, grazie a von Laue, danno esattamente M. La pressione **conta** per la gravità punto per punto,
  ma in un oggetto stabile si cancella.

## Quanti parametri

* Unità di Skyrme: **0 parametri liberi** (tutti i risultati di P1, P3, P4 sono numeri puri).
* Fm, MeV, Pa: **2 parametri di FIT** (e, F_π da ANW), non riaggiustati qui.
* Protone: tripolo JLab con 3 numeri presi dall'articolo (d₁ = −2.04, M² = 1.02 GeV², α = 3), non rifittati;
  D = (4/5) d₁. Massa del pione = 0 (limite chirale): un 3° parametro che non abbiamo usato.

## Debiti e limiti, detti chiari

* M² = 1.02 GeV² del fit JLab citato a memoria (vedi P2); il resto dei numeri JLab verificato sul testo.
* Massa classica del nodo (large-N_c, come Cebulla), senza correzioni quantistiche di rotazione: lo spin/isospin del
  protone non è nel conto.
* Limite chirale (m_π = 0): con m_π = 138 MeV la coda della pressione diventa esponenziale e D cambia (Cebulla: d₁ da
  −6.6 a −4.5 a e fisso).
* Il D di JLab è solo dei quark e a una scala; il D del nodo è "tutto". Il confronto di P2 è di forma, non di precisione.
* P4: il nodo di palline non riproduce la pressione del nucleo a h = 0.2; serve un reticolo ~3 volte più fine.

## In una frase

Il nodo di Skyrme ha **la stessa architettura di pressione del protone** — spinge al centro, tira sul bordo, si annulla
in totale, zero a 0.6 fm — e mostra in modo esatto che la massa non è "pressione sommata" (che fa zero) ma **l'energia
di un equilibrio** fra un termine che spinge come la luce e uno che tira come un elastico, metà e metà; la gravità in
tutto questo pesa 10⁻³⁹.
