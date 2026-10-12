# Glueball: materia fatta solo di "luce" che si attacca a se stessa

L'idea di Valerio: *esiste una materia fatta solo di luce che si attacca a se stessa?*

Con la luce vera (i fotoni) no: due fotoni si attraversano senza toccarsi. Esiste però una "luce" che **si tocca
da sola**: il campo dei gluoni, cioè la colla che tiene insieme i quark. Le palline di questa luce si trascinano a
vicenda. Il risultato, anche senza nessun quark, è che il campo forma **grumi con massa**: i *glueball* ("palle di colla").
Nella regola non c'è nessuna massa: c'è solo un accoppiamento. La massa la produce il campo da solo.

L'abbiamo verificato sul portatile, nella versione più semplice che si può fare per davvero: **SU(2) in 2+1
dimensioni** (due di spazio e una di tempo), su reticolo. Il risultato in breve:

- **Confinamento: sì.** La "corda" fra due cariche ha una tensione σ costante, e a√σ coincide con la letteratura
  entro lo 0,6–1,6% a tutti e tre i β (G1 PASS).
- **Grumi con massa: sì.** Lo stato più leggero (0++) pesa m ≈ 4,7–4,9 √σ, contro 4,74 della letteratura. Ai due
  reticoli più grossolani coincide entro l'1%. A quello più fine (β = 12) esce **più alto del 5,8%**, e questo è un
  **FAIL** pre-registrato (G2a), con la causa non trovata.
- **Rapporto 2++/0++** (un numero puro, senza parametri): 1,53 ± 0,07 al reticolo più fine, contro **1,639**:
  **PASS** secondo la regola fissata prima. Ma è basso del 7%, e l'estrapolazione al continuo (1,45 ± 0,09) è
  bassa dell'11%: passa solo perché l'errore è grande.
- **Primo eccitato 0++\***: **FAIL** al reticolo più fine (1,71 ± 0,07 contro 1,449). La nostra base di operatori
  non lo separa bene dagli stati più su.
- **La massa è energia (G4): non concluso.** La prova con il glueball in movimento è troppo rumorosa.

Codice: `lab/glue_core.py` (motore), `lab/glue_run.py` (Monte Carlo a pezzi da ≤ 9 minuti, riprendibile),
`lab/glue_analysis.py` (correlatori, metodo variazionale, regola del plateau, jackknife), `lab/glue_verdict.py`
(voti), `lab/glue_posthoc.py` (diagnostica a posteriori, non valutata), `lab/glue_figure.py` →
`lab/out/glue_figure.png`. Pre-registrazione: `lab/glue_prereg.py` → `lab/out/glue_prereg.json`, nel commit
**74d3747**, fatto prima di tutte le corse di produzione. Risultati: `lab/out/glue_results.json`,
`lab/out/glue_verdict.json`, `lab/out/glue_posthoc.json`.

---

## 1. Che cos'è un glueball (come a un ventenne)

Pensa alla luce come a tante freccette (il campo) attaccate ai lati di una griglia. Nell'elettromagnetismo una
freccetta è solo un angolo, cioè una fase. Le fasi si sommano e basta: due onde di luce si attraversano senza
accorgersi l'una dell'altra.

Nella teoria dei gluoni la freccetta è una piccola **rotazione** in uno spazio interno, una matrice SU(2), cioè un
punto su una sfera a tre dimensioni. Le rotazioni **non commutano**: girare prima così e poi cosà dà un risultato
diverso dal contrario. Per questo l'energia del campo contiene termini in cui il campo agisce su se stesso. La luce
di gluoni **porta la stessa carica a cui è sensibile**, e quindi si attira da sola.

Che cosa succede allora:
1. le linee di forza non si allargano come nell'elettromagnetismo, ma si stringono in un **tubo** (una corda);
2. la corda ha un'energia per metro costante (la **tensione di corda σ**), e quindi due cariche non si separano
   mai: è il **confinamento**;
3. anche senza cariche, il campo può chiudersi su se stesso in un anello o in un grumo che resta in piedi da solo.
   Quel grumo ha un'energia minima, cioè **una massa**: è il **glueball**.

I nomi: 0++ vuol dire spin 0, parità +, coniugazione di carica + (il più leggero). 0++\* è il suo primo stato
eccitato, la stessa cosa ma "vibrante". 2++ ha spin 2: un grumo che gira, con una forma a "ellisse".

## 2. La regola (e perché non c'è nessuna massa)

- Reticolo T × L × L di punti, periodico. Su ogni lato c'è una matrice SU(2) U, salvata come quaternione
  unitario (4 numeri reali).
- Energia (azione di Wilson): per ogni quadratino (placchetta) si prende il prodotto dei 4 lati attorno,
  U_P, e si paga β·(1 − ½ tr U_P). Se il campo è "liscio" U_P ≈ 1 e non si paga niente.
- **L'unico numero della regola è β = 4/(a g²)**: la forza dell'accoppiamento in unità del passo di reticolo a.
  **Non c'è nessun termine di massa.** Classicamente queste onde sono senza massa, come la luce.
- Dinamica Monte Carlo: heat-bath di Kennedy–Pendleton più 3 sweep di over-relaxation. Il primo pesca ogni
  lato dalla sua distribuzione esatta, il secondo è un "rimbalzo" che non cambia l'energia e decorrela più in fretta.

**G3: da dove viene allora la scala?** In 2+1 dimensioni l'accoppiamento g² ha le dimensioni di una massa: è
**l'unica scala della teoria**, e ogni massa deve essere (numero puro) × g². Il reticolo serve solo a calcolare quel
numero puro. Al variare di β il passo a cambia, ma il prodotto β·a√σ/4 = √σ/g² deve andare verso una costante.
Noi troviamo β·a√σ = 1,516, 1,453, 1,445 a β = 6, 9, 12, che si estrapola a **√σ/g² = 0,342 ± 0,007**, contro
**0,3353** di Teper (entro il 2%): **PASS**. Onestà: in 3+1 dimensioni (il mondo vero) g è adimensionale e la
scala nasce dalla "trasmutazione dimensionale" (la corsa dell'accoppiamento). In 2+1 D la scala è g² stesso.
In tutti e due i casi però **nessuno mette una massa a mano**: i grumi pesanti nascono da un campo senza massa.
Questo è il punto di Valerio, e regge.

## 3. Che cosa abbiamo misurato

Tre reticoli, scelti prima e scritti nella pre-registrazione:

| β | reticolo | misure | a√σ (Athenodorou–Teper 2017) |
|---|---|---|---|
| 6 | 20² × 24 | 8000 | 0,251 |
| 9 | 28² × 32 | 6000 | 0,161 |
| 12 | 40² × 32 | 7500 (5000 + 2500, vedi §6) | 0,118 |

Il lato spaziale è circa 5/√σ dappertutto, abbastanza grande perché due "corde avvolte" costino più del 2++.

**Operatori.** Per "vedere" un glueball si crea un disturbo del campo a un certo tempo e si guarda quanto ne resta
più tardi. Il disturbo è un anello di Wilson nel piano spaziale (quadrati l×l e rettangoli 2l×l), costruito con
lati "levigati" (smearing APE a 3 livelli) e sommato su tutto il piano, così il grumo è fermo (impulso zero).
- 0++: quadrati, più rettangolo orizzontale **+** verticale;
- 2++: rettangolo orizzontale **−** verticale (cambia segno ruotando di 90°: è la "forma a ellisse" di spin 2);
- corda: un anello che fa tutto il giro della scatola (il *torelone*, una corda avvolta). Dalla sua energia
  m_P = σL·√(1 − π/(3σL²)) (formula di Nambu–Goto) si ricava σ.

**Correlatore ed energia.** C(t) = ⟨O(t) O(0)⟩ − ⟨O⟩². Per tempi lunghi C(t) ∝ e^(−m t): l'esponente **è
l'energia** dello stato più leggero con quei numeri quantici. A impulso zero quell'energia è proprio la massa.
Si combinano tutti gli operatori di un canale (metodo variazionale, GEVP a t = 0/1) per trovare la miscela che
assomiglia di più al grumo. Poi si guarda la "massa efficace" m_eff(t) = ln C(t)/C(t+1) (forma cosh).

**Regola del plateau, fissata prima:** t_min è il primo t ≥ 1 in cui m_eff(t) e m_eff(t+1) hanno ciascuna un
errore relativo ≤ 20% e coincidono entro 1σ. Il fit va da t_min fino all'ultimo t con segnale/rumore ≥ 3. Se non
esiste un t così, si usa t = 1 (ripiego segnalato). Gli errori vengono dal jackknife su 40 blocchi.

## 4. I numeri contro la letteratura

Letteratura controllata prima della pre-registrazione:
- **Athenodorou & Teper 2017** (JHEP 02 (2017) 015, arXiv:1609.03873), continuo SU(2): m(0++)/√σ = **4,7367(55)**,
  m(0++\*)/√σ = **6,861(14)**, m(2++)/√σ = **7,762(10)**. Quindi **2++/0++ = 1,639** e **0++\*/0++ = 1,449**.
- **Teper 1998** (PRD 59 014512): 4,718(43), 6,83(10), 7,82(14) (2++/0++ = 1,657). √σ/g² = 0,3353(18).

### G1 – confinamento
| β | a√σ noi (torelone) | a√σ AT2017 | Creutz χ(6,6) − π/(24·5,5²), diviso σa² |
|---|---|---|---|
| 6 | 0,2526(26) | 0,2510 | 1,02 |
| 9 | 0,1615(15) | 0,1608 | 1,08 |
| 12 | 0,1204(10) | 0,1185 | 0,99 |

G1a (a√σ entro il 3%): **PASS** a tutti e tre. G1b (rapporti di Creutz ≈ σ più la correzione di Lüscher entro il
15%): **PASS** a tutti e tre. I rapporti di Creutz scendono con R e si appoggiano sulla tensione della corda: è il
confinamento visto direttamente (pannello (c) della figura).

### G2 – spettro
| β | a m(0++) | AT2017 | a m(0++\*) | a m(2++) | 2++/0++ | 0++\*/0++ | m(0++)/√σ | m(2++)/√σ |
|---|---|---|---|---|---|---|---|---|
| 6 | 1,193(30) | 1,191 | 1,88(10)† | 2,12(14)† | 1,78(13) | 1,58(9) | 4,72(13) | 8,40(57) |
| 9 | 0,770(13) | 0,762 | 1,205(25) | 1,213(30)† | 1,575(47) | 1,564(41) | 4,77(9) | 7,51(20) |
| 12 | 0,592(8) | 0,559 | 1,010(39) | 0,907(40) | **1,532(70)** | **1,706(70)** | 4,92(8) | 7,53(33) |
| continuo (fit in (a√σ)²) | | | | | **1,450(85)** | 1,650(73) | 4,93(10) | 7,11(39) |
| letteratura | | | | | 1,639 | 1,449 | 4,737 | 7,762 |

† plateau non trovato: preso a t = 1 (il ripiego della regola). Al β più grossolano gli stati pesanti muoiono in
2 passi di reticolo.

Voti (regola fissata prima: |noi − lett.| ≤ 2σ + tolleranza):
- G2a massa 0++ entro il 5% di AT2017: β = 6 **PASS**, β = 9 **PASS**, β = 12 **FAIL** (+5,8%).
- G2b **2++/0++ al β più fine: PASS** (1,53 ± 0,07 contro 1,64; tolleranza 2σ + 5%).
- G2c **0++\*/0++ al β più fine: FAIL** (1,71 ± 0,07 contro 1,45).
- G2d al β più fine: m(0++)/√σ **PASS** (4,92 contro 4,74), m(2++)/√σ **PASS** (7,53 contro 7,76),
  m(0++\*)/√σ **FAIL** (8,39 contro 6,86).
- G2e continuo: 2++/0++ **PASS**, 0++\*/0++ **PASS**, m(0++)/√σ **PASS**, m(2++)/√σ **PASS**. Questi PASS vanno
  letti con cautela. Tre punti con errori grandi danno un continuo con errori grandi: 2++/0++ = 1,45 è basso
  dell'11%, e 0++\*/0++ = 1,65 "passa" anche se a β = 9 e 12 il valore è alto (+8%, +18%).

**Il messaggio fisico regge comunque.** Il 2++ è circa **1,5–1,6 volte** più pesante dello 0++, ed entrambi stanno
in rapporti fissi con la tensione di corda (m ≈ 4,8 √σ e ≈ 7,5 √σ), senza nessuna massa nella regola. Con la
precisione di un portatile riproduciamo i numeri dei professionisti entro il 3–7% per 0++, 2++ e σ. Non li
riproduciamo per lo 0++\*.

### G4 – la massa è energia
Una cosa vera ma quasi banale: la massa del glueball **è** un'energia. È l'esponente con cui decade C(t), cioè
la differenza fra l'energia del grumo e quella del vuoto (la matrice di trasferimento del reticolo è e^(−aH)).
Nessuno ha messo "materia": c'è solo energia del campo, ferma in un grumo, e quella energia, divisa c², è la massa.

La prova non banale che avevamo pre-registrato: un grumo in **movimento** (impulso p = 2π/L) deve avere
E(p)² = m² + c²p², con c = 1 (la velocità della luce del reticolo). Risultato: c² = 1,10 ± 0,45 a β = 9 e
−0,77 ± 0,65 a β = 12. Errori troppo grandi per decidere: **INCONCLUSIVE** a tutti e due. A β = 12 l'energia del
grumo in moto, E(p) = 0,576(11), esce **più bassa** della massa a riposo misurata, 0,592(8): impossibile per un grumo
vero. È lo stesso sintomo del FAIL di G2a (vedi §5). Notiamo, a posteriori, che E(p) è proprio quello che ci si
aspetta con la massa della letteratura (√(0,559² + 0,157²) = 0,581).

## 5. Fallimenti, onestamente

1. **G2a a β = 12: la massa 0++ esce 0,592(8) invece di 0,559 (+5,8%, circa 4σ).** Diagnostica a posteriori
   (`glue_posthoc.py`, non valutata):
   - le due catene indipendenti danno m_eff(1) = 0,594 e 0,572, una differenza di 2–3σ;
   - il tempo di autocorrelazione è piccolo (τ_int ≈ 1 misura), quindi gli errori del jackknife non sembrano
     sottostimati per colpa delle autocorrelazioni;
   - il canale con impulso dà un valore compatibile con la letteratura;
   - a β = 6 e 9 la stessa catena di codice dà la letteratura entro l'1%.

   Non abbiamo trovato la causa. Tra i candidati: una fluttuazione statistica più grande degli errori stimati, oppure
   il fatto che a β = 12 il grumo è grande (raggio di circa 8 passi) e i nostri operatori più grandi (8×8) lo
   coprono appena. Un test pulito sarebbe rifare β = 12 con operatori più grandi e più statistica. Non lo facciamo
   qui, perché sarebbe una scelta fatta guardando il risultato.
2. **0++\* (G2c, G2d):** il secondo stato del canale 0++ esce sempre più pesante del valore vero (+8% a β = 9,
   +18% a β = 12). I nostri operatori sono quasi tutti la stessa cosa a scale diverse (correlazioni fino a
   0,9999), quindi la base non ha abbastanza "forme diverse" per isolare lo stato eccitato. Probabilmente
   vediamo un misto di 0++\* e 0++\*\* (8,38 √σ).
3. **G4 non concluso**: serviva una precisione sull'energia dello 0,3%, e non l'avevamo.
4. **β = 6**: gli stati pesanti si leggono solo a t = 1 (ripiego della regola). Lì il reticolo è grossolano.

## 6. Cose dichiarate (deviazioni e note di processo)

- **Pilota** prima della pre-registrazione: β = 6, 400 misure, solo per provare il codice. Ha mostrato che una
  regola del plateau basata solo sull'accordo entro 1σ accetta rumore, e per questo è stato aggiunto il requisito
  "errore relativo ≤ 20%". È scritto nel JSON della pre-registrazione.
- **Accelerazioni del codice durante le corse.** Dopo il primo pezzo di β = 12 abbiamo riscritto in numba il
  calcolo degli anelli, e la placchetta ora viene dalla densità d'azione già calcolata. I risultati sono identici
  bit per bit a quelli della versione numpy (differenza 0,0, verificato), quindi la misura non cambia. numba 0.68
  era già installato. Ogni processo usa un solo thread.
- **β = 12: seconda catena indipendente** (seme diverso, 2500 misure) lanciata quando un core si è liberato, per
  avere più statistica. La catena principale ha completato le 5000 misure pre-registrate. L'analisi usa 5000 + 2500
  = 7500 misure, e la regola di arresto ("finire entrambe le catene") è stata fissata prima di guardare qualsiasi
  numero di β = 12 oltre a due controlli intermedi. Onestà: i controlli intermedi ci sono stati (a 1331 e 2401
  misure, am0 ≈ 0,60) e non hanno cambiato niente del piano: le catene sono state portate alla lunghezza decisa.
- I dati grezzi (`lab/out/glue_raw/`, circa 420 MB) non sono nel repository: si rigenerano con `glue_run.py` e
  gli stessi semi. I semi della catena cambiano a ogni pezzo da 9 minuti: è sempre una catena di Markov valida,
  ma non è riproducibile bit per bit se la si spezza in punti diversi.

## 7. Debiti (che cosa NON abbiamo fatto)

- **2+1 D, non 3+1 D.** Il mondo vero ha tre dimensioni di spazio. In 3+1 D i glueball di SU(3) hanno m(0++)
  ≈ 1,6–1,7 GeV e m(2++)/m(0++) ≈ 1,4. Lì serve un reticolo 4D, troppo per questo portatile in poche ore.
- **SU(2), non SU(3).** I gluoni veri hanno 3 colori (SU(3)). SU(2) è il fratello più piccolo, con le stesse idee
  (non abeliano, confinamento, gap di massa). I rapporti di massa in 2+1 D dipendono poco da N (Teper 1998).
- **Niente quark.** È la teoria "pura colla". Nel mondo vero i glueball si mescolano con i mesoni fatti di quark,
  e nessuno li ha ancora identificati con certezza negli esperimenti.
- **Artefatti di reticolo.** Tre soli β e una sola dimensione di scatola per β. L'estrapolazione al continuo con
  tre punti ha errori grandi. Gli effetti di volume finito non li abbiamo misurati: ci siamo affidati alla
  scelta L√σ ≈ 5.
- **Parametri messi a mano: zero oltre al reticolo.** β, L, T e i semi sono scelte di calcolo, non fisica. Nella
  regola non c'è nessuna massa. Gli unici fit sono gli esponenziali dei correlatori e la retta per il continuo.

## 8. Figura

`lab/out/glue_figure.png`:
(a) masse efficaci a β = 12 per 0++, 0++\* e 2++, con le bande del fit e i valori della letteratura;
(b) m/√σ ai tre β con l'estrapolazione al continuo, contro Athenodorou–Teper;
(c) rapporti di Creutz divisi per σ, che scendono verso 1 (confinamento);
(d) √σ/g² contro 1/β: l'unica scala è g² (G3), più i numeri del test di dispersione (G4).
