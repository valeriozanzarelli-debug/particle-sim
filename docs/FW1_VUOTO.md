# Il vuoto di FW1: polarizzazione senza proiettore

Studio `fw1vuoto_`. Variante in prova **FW1** (`docs/SENZA_DOPPIONI.md` §5.3) accanto alla regola congelata, che resta
il riferimento. Codice: `lab/fw1vuoto_core.py` (motore), pre-registrazione `lab/fw1vuoto_prereg.py` →
`lab/out/fw1vuoto_prereg.json` (commit prima di ogni corsa di produzione).

*(Documento in costruzione: sezioni 1–2 scritte prima delle corse.)*

---

## 1. Il problema in parole semplici

Il vuoto della regola zig-zag è il "mare di Dirac": tutti gli stati a energia negativa pieni. Un nucleo lo deforma
(polarizzazione del vuoto): la carica vista da lontano è più piccola di quella nuda di un fattore che cresce col
logaritmo del taglio (coefficiente 2α/3π), e a distanze ~λ_C compare il potenziale di **Uehling**.

Il thread Debiti (`docs/VUOTO_POLARIZZA.md` su origin/master, fasi 1–4) ha trovato che con la regola congelata:

- il **mare intero** non va bene: in ogni cella della zona ci sono 4 specie pesanti normali e **4 coni senza massa**
  a quasi-energia π, riempiti "al contrario"; il risultato è un anti-schermo che cresce a grande distanza;
- contando **solo la specie fisica** (un proiettore dichiarato, K1 = 0,6, K2 = 1,05) la regola ridà la QED:
  log 2α/3π a −0,07 %, Uehling entro 0,1–0,2 %;
- in fase 4 la luce non si poteva vestire col mare intero (i coni a π chiudono la zona di Floquet), e il taglio di
  specie rompe la gauge (massa del fotone da togliere a mano, c della luce sporcato di 0,007–0,05 α).

FW1 (studio SENZA_DOPPIONI) toglie i doppioni leggeri: tutti i gemelli stanno a ≥ 0,37/a da energia 0 e ≥ 0,40/a dal
bordo π; le coppie dal vuoto contate su tutta la zona danno 0,99 di una specie. Se è vero che i gemelli pesanti sono
"innocui", il **mare intero** di FW1 deve polarizzare come la QED **senza nessun proiettore**: i gemelli possono
cambiare solo la costante della carica nuda (un controtermine), non il logaritmo né la forma di Uehling.

## 2. Cosa si sa (letteratura)

Non ho potuto verificare online in questa sessione: cito a memoria, con autori e anno.

- **Uehling (1935), Serber (1935)**: il potenziale di polarizzazione del vuoto al primo ordine in α. Il coefficiente
  del log della rinormalizzazione della carica, 2α/3π per specie, è di Dirac–Heisenberg–Weisskopf (anni '30).
- **Fermioni di Wilson sul reticolo (Wilson 1974; Karsten–Smit 1981; Kawai–Nakayama–Seo 1981)**: nel reticolo
  euclideo la polarizzazione del vuoto con fermioni di Wilson è invariante di gauge e i doppioni (masse ~ r/a)
  contribuiscono solo a costanti e a termini O(a²): il log viene da una sola specie. È il risultato standard che qui
  si vuole **ritrovare** per una regola a tempo discreto: validazione, non scoperta.
- **Massa ritorta (twisted mass, Frezzotti–Grassi–Sint–Weisz 2001, Frezzotti–Rossi 2004)**: un termine di Wilson
  "ruotato" rispetto alla massa (come il Wilson lungo Y di FW1, ortogonale alla moneta lungo X) rompe parità e
  inversione del tempo ma conserva la loro combinazione; a torsione massima gli errori O(a) si cancellano.
- **Tempo discreto contro reticolo 4-D**: nel reticolo euclideo la simmetria ipercubica vieta termini di dimensione 4
  che rompono Lorentz nel settore del fotone. Una regola a tempo discreto con solo simmetria cubica **nello spazio**
  (anzi, la congelata ha solo D4h con y speciale) non ha questa protezione: la velocità della luce vestita può
  differire da c di O(α) senza calare con a (Collins, Perez, Sudarsky, Urrutia, Vucetich, PRL 93, 191301 (2004)).
- **Termini dispari**: Carroll–Field–Jackiw (PRD 41, 1231 (1990)): un termine di Chern–Simons b_μ ε^{μνρσ}A_ν F_ρσ
  (CPT-dispari) rende il vuoto birifrangente; Jackiw–Kostelecký (PRL 82, 3572 (1999)) e altri: un tale termine
  indotto dai fermioni è ambiguo/dipende dal regolatore. Settore del fotone del Modello Standard Esteso
  (Colladay–Kostelecký 1997–98, Kostelecký–Mewes 2002): tutti i coefficienti di dimensione 3 (CPT-dispari) e 4
  (CPT-pari, tra cui il termine magneto-elettrico E·κ·B, P-dispari e T-dispari). Effetto magneto-elettrico/assione:
  Wilczek (1987).
- **Passeggiate quantistiche**: non conosco lavori che calcolino la polarizzazione del vuoto del mare di Floquet di un
  automa di Dirac 3-D (D'Ariano–Perinotti, Arrighi, Bisio studiano la dinamica a una particella). Il thread Debiti è,
  per quanto so, il primo; questo studio ne è la continuazione con FW1.

**Cosa è aperto qui** (per la nostra regola, non per la fisica): se FW1 dia un vuoto sano senza correzioni a mano, e
quanto la sua inversione del tempo rotta (coppie di Kramers spaccate 5·10⁻²) sporchi il vuoto.
