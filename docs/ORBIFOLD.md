# Le forme nascoste "alla Calabi–Yau" fatte nel modo più semplice: gli orbifold

Laboratorio: `lab/orb_*.py` (nucleo piccolo `lab/orb_core.py`). Numeri: `lab/out/orb_o1.json`, `orb_o2.json`,
`orb_o3.json`. Figura: `lab/out/orb_figure.png`.
Pre-registrazione: `lab/orb_prereg.py` → `lab/out/orb_prereg.json`, **messa in git (commit 61486a9) prima di
qualsiasi corsa**. Tempo totale delle corse: circa 2,5 minuti.

## In due righe

Un orbifold è il modo più semplice di fare una forma nascosta "alla Calabi–Yau": prendi un toro piatto (il foglio
con i bordi incollati) e dichiari uguali i punti che una rotazione porta l'uno nell'altro. Ho messo le palline
zig-zag su questi orbifold e ho contato tutto **solo con i loro operatori**. **I numeri dei libri escono tutti, esatti**
(16 e 27 punti fissi, χ = 24 e 72, 36 famiglie, spinore che torna uguale solo nei casi giusti). Con un flusso nel toro
nascosto più una "torsione" (linea di Wilson) si arriva anche a **3 modi leggeri**, ma **solo scegliendo a mano** il
flusso e il settore. Debito grande: **il quoziente lo impongo io**, nessuna dinamica lo sceglie.

---

## 0. Controllo dei fatti (cosa era giusto, cosa no)

| affermazione | verdetto |
|---|---|
| T⁴/Z2 (x → −x) = Kummer, 16 punti fissi, χ = 24 | **giusto, con una precisazione**: il quoziente "grezzo" (con le 16 punte) ha χ = **8** e numeri di Betti (1,0,6,0,1). Diventa K3 (χ = 24, b2 = 22) solo dopo aver **lisciato** ogni punta con una sferetta (ogni punta: χ da 1 a 2, quindi +16). Il "χ di orbifold" dei fisici, (1/|G|)·Σ χ dei punti fissi delle coppie (g,h), dà già 24 senza lisciare. |
| T²/Z3 ha 3 punti fissi | **giusto**. Ma **T²/Z3 NON è un Calabi–Yau**: è una **sfera** con tre punte da 120° (χ = 2). Un Calabi–Yau di dimensione 1 è solo il toro. T²/Z3 è un mattone, non la casa. |
| T⁶/Z3: 27 punti fissi, h11 = 36, h21 = 0, χ = 72, 36 generazioni (embedding standard) | **giusto**. In dettaglio: 9 dal settore "non torto" + 27 dai punti fissi = 36; χ = 2(h11 − h21) = 72. Il quoziente grezzo ha χ = 18; lisciando ogni punta con un piano proiettivo (+2) si arriva a 18 + 54 = 72. |
| Le linee di Wilson riducono a 3 famiglie (Ibáñez–Kim–Nilles–Quevedo 1987) | **giusto** (Phys. Lett. B 191 (1987) 282; prima Ibáñez–Nilles–Quevedo, PLB 187 (1987) 25). Precisazione: su T²/Z3 la linea di Wilson **non è libera**: deve essere uguale sui due lati del rombo e di ordine 3 (W³ = 1). Il "3" delle famiglie, nel conto dei punti fissi, viene dai **3 punti fissi del toro non toccato** (vedi O3c). Quali classi portano davvero famiglie leggere lo decidono i dati del gruppo di gauge E8×E8, che qui non ho. |
| CY ⇔ spinore covariantemente costante (olonomia ⊂ SU(n)) | **giusto in dimensione 4 e 6** (le nostre). In dimensione 7 e 8 ci sono anche G2 e Spin(7), che hanno spinori costanti ma non sono CY. Attenzione: anche il **toro piatto** ha spinori costanti (è un CY banale), quindi lo spinore costante non garantisce una forma "complicata". |
| "regola doppia" ↔ Z2, "tre punti a 120°" ↔ Z3 | **analogie**, non identità: il Z2 di Kummer agisce su **tutte e 4** le direzioni nascoste insieme; il Z3 di T⁶/Z3 agisce **sui tre tori insieme**. Se agisce su uno solo, lo spinore non torna (vedi i controlli di O1). |

---

## O1 — L'olonomia vista dalle palline

**Come.** (a) Pallina zig-zag in 4 dimensioni nascoste (rete n⁴ con n = 4 e 6): trasla lungo ogni asse secondo le 4
matrici di spin, poi gira (massa θ = 0,3). (b) Per i casi con 120° la regola a passi su rete quadrata non ha la
simmetria di 120°, quindi uso la **versione continua nel tempo** della stessa pallina su rete esagonale (Dirac con
termine di Wilson, exp(−iH dt)). Dichiarato prima.

Non ho scritto io come la rotazione agisce sullo spinore: l'ho **ricavato dalla pallina**: cerco tutte le matrici S
tali che "muovere e poi ruotare" = "ruotare e poi muovere". Poi chiedo che S sia **la rotazione vera della pallina**
(fatta con i suoi generatori), e conto gli spinori con S·χ = χ: sono quelli che **tornano uguali dopo un giro intorno a
un punto fisso**.

| caso | la pallina accetta la rotazione? | spinori che tornano uguali |
|---|---|---|
| **K3 = T⁴/Z2** (n = 4 e 6) | sì, scarto 2,7·10⁻¹⁶ | **2** (stessa chiralità) ✔ |
| **T⁶/Z3** (tre tori, 120° insieme) | sì, scarto 1,6·10⁻¹⁴ | **2** ✔ |
| T²/Z3 (previsto 0: non è CY) | sì | **0** ✔ |
| controllo: Z2 su una sola coordinata complessa di T⁴ | sì | **0** ✔ |
| controllo: Z3 su un solo toro di T⁶ | sì | **0** ✔ |
| controllo: riflessione x1 → −x1 con massa | **no**: nessuna S esiste (spazio delle soluzioni = 0) | – ✔ |

Un pacchetto messo nel "settore orbifold" di K3 e fatto camminare 60 tick **non esce** dal settore (perdita 2·10⁻¹⁵).
**O1: PASSA** (6/6 come pre-registrato).

**Due cose oneste.**
* Su un toro, le rotazioni intorno ai 16 (o 3) punti fissi differiscono solo per una traslazione di un giro intero,
  che sul toro è l'identità: quindi "il giro intorno a ogni punto fisso" è **lo stesso operatore**. Ho controllato
  l'olonomia una volta per rete (K3) e 9 combinazioni di centri (T⁶/Z3), tutte uguali.
* **La scappatoia della fase.** La pallina fissa S solo a meno di una fase. Se permetto una **fase U(1) in più** (cioè
  una torsione del cerchio nascosto, una carica KK), allora anche i controlli "non CY" hanno spinori che tornano uguali
  (Z2 su una coordinata: 2; Z3 su un solo toro: 4; T²/Z3: 1; barre a righe in figura). Quello non è più geometria: è uno spinore **carico** sotto un altro
  cerchio. Il test vale solo se la fase è vietata.

## O2 — Contare con gli operatori delle palline

**Come.** Punti fissi: contati uno per uno sulla rete. Numeri di Betti (quanti pezzi, anelli, cavità...): dalla rete
fatta a triangoli (T², L = 6 e 9) o a cubetti (T⁴, n = 4), tenendo solo le "catene" che la rotazione lascia uguali (sono
i modi zero del Laplaciano delle celle sul quoziente). T⁶/Z3: dal modo in cui la rotazione agisce sugli anelli di
**un** toro esagonale, moltiplicato per tre tori (Künneth). Modi zero di Dirac: la pallina di Wilson-Dirac ristretta
agli stati che la rotazione lascia uguali.

| | punti fissi | Betti del quoziente grezzo | χ grezzo | χ dopo lisciatura | altro |
|---|---|---|---|---|---|
| T²/Z3 | **3** | (1, 0, 1) = sfera | **2** | 2 (in 1-D non serve) | Dirac: **0** modi zero |
| T⁴/Z2 | **16** | (1, 0, 6, 0, 1) | **8** | **24** = 8 + 16, b2 = 22 | χ orbifold 24; Dirac: **2** modi zero, stessa chiralità (= Â(K3) = 2) |
| T⁶/Z3 | **27** | (1, 0, 9, 2, 9, 0, 1) | **18** | **72** = 18 + 2·27 | h11 = 9 + 27 = **36**, h21 = **0**, **36** famiglie |

**O2: PASSA**, con una **deviazione dalla pre-registrazione**: per T⁴/Z2 i Betti completi li ho calcolati solo con
n = 4 (con n = 6 i ranghi densi erano troppo lenti per il tempo); con n = 6 ho controllato punti fissi (16), χ grezzo
(8, per media di Lefschetz) e modi zero di Dirac (2).

**Debito dichiarato.** La lisciatura (+1 per ogni punta di K3, +2 per ogni punta di T⁶/Z3) e i 27 contributi "torti"
a h11 **li aggiungo io** dai libri: le palline vedono il quoziente grezzo, non la versione liscia. Le palline puntiformi
non hanno un "settore torto" (quello è una cosa delle stringhe).

## O3 — La torsione del cerchio nascosto (linea di Wilson) cambia i modi leggeri

**(a) Quali torsioni sono permesse.** Sul T²/Z3 senza flusso ho provato 16 coppie di fasi sui due lati. La pallina
accetta la rotazione di 120° **solo** per (0,0), (1/3, 2/3), (2/3, 1/3) (in giri): esattamente le Z3 con W³ = 1
(nella mia base dei lati, W₂ = W₁⁻¹). ✔

**(b) Flusso + torsione.** Metto un flusso AB uniforme di M quanti attraverso il toro nascosto (come i nostri flussi
nei buchi). Il livello più basso (= modi zero del Dirac chirale) ha **esattamente M stati** per ogni M (con un salto
netto sopra). La rotazione di 120° li divide in tre settori (n0, n1, n2):

| M | W = 0 | W = (1/3, 2/3) | M dispari: W = (1/2, 1/2) |
|---|---|---|---|
| 2 | (1,1,0) | (1,0,1) | – |
| 4 | (1,1,2) | (2,1,1) | – |
| 6 | (**3**,1,2) | (2,2,2) | – |
| 7 | – | – | (2,2,**3**) |
| 8 | (**3**,**3**,2) | (**3**,2,**3**) | – |
| 9 | – | – | (4,2,**3**) |
| 10 | (**3**,**3**,4) | (4,**3**,**3**) | – |
| 11 | – | – | (4,4,**3**) |
| 12 | (5,**3**,4) | (4,4,4) | – |

Per M dispari serve la torsione (1/2, 1/2): senza, la rotazione non è compatibile col flusso. (Questo è lo stesso tipo
di regola che trova la letteratura sui "magnetized orbifolds" di Abe, Kobayashi e collaboratori; non ho confrontato
la loro tabella numero per numero.)

**Dove stanno i modi.** Con M = 6 e W = 0 i 3 modi del settore 0 pesano su **tutti e tre** i punti fissi, gli altri
settori ci si annullano. Accendendo la torsione, ogni settore si sposta su **un solo** punto fisso diverso (pannello in
basso a destra della figura). Quindi sì: **la torsione cambia quali punti fissi portano i modi leggeri.**

**Il "3" arriva?** Sì: 17 scelte su 72 (M ≤ 12, torsione compatibile, settore) danno esattamente 3. **Ma è una scelta,
non una previsione**: ogni settore ha circa M/3 modi, e 3 esce perché scegliamo M tra 6 e 12 e il settore giusto.

**(c) T⁶/Z3: i 27 punti fissi.** Con 0, 1, 2, 3 torsioni indipendenti i 27 punti si dividono in 1×27, 3×9, **9×3**,
27×1. Con due torsioni ogni classe ha **3** punti: i tre punti fissi del toro non toccato. È questo il motivo
"geometrico" per cui le famiglie vengono a gruppi di 3 nei modelli di Ibáñez e collaboratori. **Quale** classe porta
famiglie senza massa dipende dai vettori di E8×E8 e dalla condizione di livello: **non fatto**, non ho quei dati.

## O4 — Emergenza

**Saltato**, come pre-registrato (opzionale). Motivo in più, scoperto qui: T²/Z3 è una **sfera** (b1 = 0, b2 = 1), e il
lettore di forma riparato (`big_core`) vede solo anelli (H1), **non le cavità**: non potrebbe nemmeno riconoscere un
T²/Z3 se emergesse. E nel passo grande di GEOMETRIA_AUTOCOSTRUITA le palline hanno fatto solo isole.

---

## Conteggio dei parametri

* Scelte discrete fatte a mano: il **gruppo** (Z2 o Z3) e **come agisce** (su tutte le coordinate insieme); per O3
  anche il flusso **M**, la **torsione** (3 scelte, 1 per M dispari) e il **settore** (3 scelte).
* Parametri continui adattati: **zero**. La massa θ e il passo dt non cambiano nessun conteggio.
* Aggiunti dai libri, non dalle palline: le lisciature (+1, +2 per punta) e i 27 contributi torti a h11.

## Debiti

1. **Il quoziente è imposto a mano.** Nessuna regola delle palline sceglie "120°" o "x → −x".
2. Le palline vedono solo il quoziente grezzo; la forma liscia (K3 vero, CY vero) richiede la lisciatura, messa a mano.
3. La simmetria di 120° vale esattamente solo per la versione continua nel tempo della pallina, non per il tick a
   passi su rete quadrata.
4. La scappatoia della fase U(1): senza un motivo per vietarla, "spinore che torna uguale" non distingue CY da non-CY.
5. Il "3" di O3 è scelto (M e settore), e il vero meccanismo delle famiglie ha bisogno dei dati del gruppo di gauge.

## Verdetto per uno studio grande "orbifold emergente"

**Non partire ora.** La macchina funziona (tutti i numeri dei libri escono dalle palline), ma la domanda interessante,
se la forma con le punte **nasce da sola**, non si può fare finché (i) il lettore di forma non vede le cavità e le punte
(b2, o curvatura concentrata in punti), e (ii) non c'è una dinamica che abbia almeno una possibilità di rompere la
simmetria verso Z3 (il rinforzo ha dato solo isole). Passo piccolo sensato prima: aggiungere al lettore la
cavità/le punte e provarlo su T²/Z3 e T⁴/Z2 disegnati.
