# Tempo di tunnel e tempo di arrivo sulla regola zig-zag

Studio `tunnel_` (file: `lab/tunnel_*.py`, `lab/out/tunnel_*`). Regola: quella congelata di `lab/zz_core.py`
(1+1 D, `dirac1`): ogni pallina corre di una cella a c, poi ruota (flip θ = m c a/ħ) e prende la fase del
potenziale. Unità ħ = c = mₑ = 1 (lunghezza λ_C = 3,86·10⁻¹³ m, tempo ħ/mₑc² = 1,29·10⁻²¹ s).

Commit: ricerca `039613b`, pre-registrazione `2fedf0a`, corse `eab1b64`, documento e figura: vedi `git log`.
Figura: `lab/out/tunnel_figura.png` (copia in `/mnt/project-files/problemi-aperti/tunnel_figura.png`).

## 0. In breve (per chi ha fretta)

* **Il cono di luce della regola è rigido e lo abbiamo visto esatto**: in ogni corsa (6 corse, fino a 22 000 tick,
  7 copie) oltre il "raggio di luce" del pacchetto iniziale la densità è **zero esatto** (non 10⁻³⁰: zero), e
  l'ultima cella accesa avanza di **una cella per tick esatta** (velocità del fronte = c, errore 0). Il fronte
  passa la barriera **senza essere attenuato** (solo i cammini che non ruotano mai arrivano fin lì).
* **Effetto Hartman: sì, anche nella regola.** Il tempo di fase satura con lo spessore (tra L = 8 e 16 cambia
  meno di 8·10⁻⁴) e coincide con Dirac a 4·10⁻⁶. Velocità apparenti L/τ_φ fino a **4c** già a L = 4 ħ/mc.
  Con un pacchetto ultra-relativistico (E = V0 = 5 mc², L = 10) il **picco** trasmesso esce **4,05 ħ/mc
  davanti** al raggio di luce partito dal picco in arrivo (previsto 4,05). Non è una violazione: il fronte vero è
  ~135 ħ/mc più avanti; il picco è fatto dalle palline della **coda davanti** del pacchetto (lo mostrano i
  surfisti). **VALIDAZIONE** (Krekora–Su–Grobe 2001, Winful 2006).
* **Orologio di Larmor dentro la regola = dwell** (barriera simmetrica) a 2·10⁻⁷; nel limite non relativistico
  la regola ritrova Büttiker 1983 a **3·10⁻⁵**; con una barriera gaussiana ritrova la forma "più lenti → più
  veloci" di Ramos et al. 2020 (entro 0,5 %, correzioni relativistiche V0/mc² = 1 %). Il valore 0,62 ms di
  Ramos esce se la barriera è alta ~7 nK, ma l'altezza vera **non sono riuscito a leggerla**: confronto solo di
  forma. **VALIDAZIONE.**
* **L'orologio che conta le rotazioni (domanda 3) è nuovo per quanto abbiamo trovato.** Contare quante volte i
  cammini a zig-zag ruotano dentro la barriera (pesati con l'ampiezza) è *esattamente* un orologio di Larmor
  attaccato alla **massa** invece che al potenziale, e legge il **tempo proprio**. Non coincide con nessuno dei
  tempi noti. Nel tunnel spesso:
  - parte "nel piano": τ_my = (ε/m)·τ_dwell esatto (ε = E − V0): **va all'indietro** quando V0 > E;
  - parte "inclinata": cresce come **m L/(ħκ) ≥ L/c**: mai più veloce della luce, e **esattamente L/c** al centro
    del gap (V0 = E), dove invece Büttiker–Landauer dà 0 (tunnel "istantaneo", come l'attoclock);
  - nel limite non relativistico torna Larmor (Büttiker).
  Tutti PASS (6 esami pre-registrati, E11a–f). È una **SCOPERTA piccola**: matematicamente semplice, segue da Dirac; il
  nuovo è la lettura (rotazioni = tempo proprio) e la previsione falsificabile per un simulatore di Dirac.
* **Surfisti (Bohm):** tutti i trasmessi stanno sotto c (minimo 75 ħ/mc² per L = 10), ma il loro tempo nella
  barriera (~115) **non è una proprietà della barriera**: cresce come ≈ 5,2 σ + 6 con la larghezza del pacchetto
  (esplorativo, 4 larghezze).
* **Tempo di arrivo: FALLITO quello che avevo previsto.** Il rivelatore "ideale" della regola (toglie tutto ciò
  che corre oltre X_D, la condizione "solo onda uscente") **riflette** una frazione (1 − v)/(1 + v) dell'onda
  (misurato 0,07188 contro 0,07180 per v = 0,866) e con la barriera forma una cavità: conta ±10 % di
  probabilità rispetto a quella che passa davvero (esame E13, soglia 10⁻⁶: FAIL). Il tempo di arrivo dipende dal
  rivelatore: limite onesto, coerente con Alonso–Sala Mayato–Leavens (prima-arrivata ben definita solo per v → c).
* **Esami: 18 pre-registrati, 17 PASS, 1 FAIL (E13).** Parametri liberi: **0**.


---

## 1. Il problema e cosa si sa (ricerca, 13/10/2026)

**La domanda.** Una particella passa una barriera "per effetto tunnel" (energia sotto il bordo della
barriera). Quanto tempo ci ha passato dentro? In meccanica quantistica il tempo non è un operatore (Pauli,
1933: un operatore tempo autoaggiunto coniugato a un'energia limitata dal basso non esiste), quindi ci sono
**tante definizioni diverse**, che danno numeri diversi. Le principali:

| tempo | definizione | cosa dice | chi |
|---|---|---|---|
| **fase (Wigner, gruppo)** | ritardo del picco: τ_φ = ħ ∂(arg t + kL)/∂E | per barriere spesse **satura** (non cresce con L): "effetto Hartman"; la velocità apparente L/τ_φ supera c | Wigner 1955; Hartman, *J. Appl. Phys.* 33, 3427 (1962) |
| **permanenza (dwell)** | τ_D = ∫_barriera \|ψ\|² dx / flusso incidente | tempo medio passato dentro, senza distinguere chi passa da chi torna indietro; anche lui satura | Smith 1960; Büttiker 1983 |
| **Larmor (Baz', Rybachenko, Büttiker)** | un orologio: un piccolo campo magnetico solo dentro la barriera fa ruotare lo spin; l'angolo misura il tempo. Due componenti: τ_y (rotazione nel piano) = −ħ ∂arg t/∂V, τ_z (inclinazione) = −ħ ∂ln\|t\|/∂V | per una barriera simmetrica τ_y = τ_D; τ_z per barriere spesse → m L/(ħκ) (tempo "semiclassico" di Büttiker–Landauer, cresce con L) | Baz' 1967; Büttiker, *PRB* 27, 6178 (1983); Büttiker e Landauer, *PRL* 49, 1739 (1982) |
| **valore debole** | τ_y − i τ_z = valore debole (complesso) del "tempo passato dentro" sui soli trasmessi | parte reale = tempo vero, parte immaginaria = disturbo della misura | Steinberg, *PRL* 74, 2405 (1995) |
| **Bohm** | traiettorie che seguono la corrente (v = j/ρ); in 1-D non si incrociano | i trasmessi sono quelli partiti davanti; tempo finito, mai sopra c nella versione di Dirac | Leavens (anni '90); Hauge e Støvneng, *RMP* 61, 917 (1989, rassegna) |
| **attoclock** | ionizzazione di atomi da laser ellittico: angolo dell'elettrone uscito | in idrogeno tempo di tunnel compatibile con **zero** (< 1,8 as) | Torlina et al., *Nat. Phys.* 11, 503 (2015); Sainadh et al., *Nature* 568, 75 (2019) |

**Le misure.** Ramos, Spierings, Racicot, Steinberg, *Nature* 583, 529 (2020) (arXiv 1907.13523, letto
l'abstract online): atomi di ⁸⁷Rb condensati (energia ~1 nK), barriera ottica spessa 1,3 µm, orologio di
Larmor fatto con una transizione Raman solo dentro la barriera: **tempo di attraversamento 0,62(7) ms**
(un riassunto di corso, Illinois PHYS596 2022, riporta τ_y = 0,61(7) ms al picco della distribuzione di
velocità e la previsione "controintuitiva" che **energie più basse danno tempi più corti**). Altezza della
barriera e larghezza gaussiana esatta: **non sono riuscito a leggerle** (il testo completo non era
raggiungibile: limite di accesso del proxy); le uso solo come forma, vedi §6. Con la luce: Steinberg, Kwiat,
Chiao, *PRL* 71, 708 (1993) (fotoni, ritardo ~1,5 fs "superluminale"); Spielmann et al., *PRL* 73, 2308
(1994) (impulsi in specchi a strati).

**Relatività (il nostro terreno).** Con l'equazione di Dirac la storia cambia poco ma in modo importante:
* Krekora, Su, Grobe, *PRA* 63, 032107 (2001) (abstract letto online): pacchetti di Dirac attraverso una
  barriera; picco e fronte possono sembrare superluminali, ma **Dirac non viola la causalità** (Schrödinger sì,
  perché non ha un cono di luce). Anche *PRA* 64, 022105 (2001) e *PRL* 92, 040406 (2004) (Klein).
* Petrillo e Janner, *PRA* 67, 012110 (2003); De Leo e Rotelli, *EPJ C* 51, 241 (2007); De Leo, *J. Phys. A*
  46, 155306 (2013) (letto online): nel limite opaco il tempo di fase dà velocità apparenti oltre c; De Leo
  conclude che il tempo di fase "non è il tempo di tunnel giusto" per l'effetto filtro (passano solo le
  componenti più energetiche) e propone come lavoro futuro le **correzioni relativistiche a Larmor**.
* Li, Liang, Kobe (arXiv quant-ph/0111117, letto l'abstract): Larmor relativistico per una particella neutra di
  spin ½: il tempo di precessione è uguale al dwell (come in Schrödinger).
* Flores e Galapon (arXiv 2207.09040, 2022, letto l'abstract): operatore "tempo di arrivo" relativistico per
  spin 0: le componenti sotto barriera passano in tempo zero.
* Winful, *Phys. Rep.* 436, 1 (2006) (non verificato online): il ritardo di gruppo è un tempo di "accumulo e
  rilascio" di energia, non di transito; la saturazione viene dalla saturazione dell'energia immagazzinata
  (niente superluminalità vera).

**Tempo di arrivo.** Nessun operatore accettato (Pauli). Candidati: flusso j(X,t) nel punto del rivelatore
(può essere negativo: *backflow*); regola del bordo assorbente di Tumulka (2016–2022, non verificato online:
per Dirac la condizione al bordo "solo onda uscente" n·α ψ = ψ); prima-arrivata sui cammini a scacchiera di
Feynman: Alonso, Sala Mayato, Leavens (arXiv quant-ph/0206071, letto l'abstract): con i cammini di Feynman la
distribuzione di prima-arrivata **non** si definisce in modo coerente per l'interferenza tra prime e
successive arrivate (tranne vicino a v → c).

**Il conteggio delle svolte.** Jacobson e Schulman, *J. Phys. A* 17, 375 (1984) (riferimento verificato, il
contenuto no): nella scacchiera di Feynman il numero di svolte è legato al passaggio relativistico →
non relativistico (in tempo immaginario le svolte sono un processo di Poisson con ritmo ∝ m c²/ħ).

**Cosa è aperto.** (i) Quale tempo è "quello vero" (probabilmente nessuno: domande diverse, risposte
diverse); (ii) se l'effetto Hartman dica qualcosa sul transito o solo sul rimodellamento del pacchetto
(Winful: il secondo); (iii) un tempo di arrivo con una distribuzione di probabilità accettata; (iv) cosa
misura davvero l'attoclock.

**Cosa può dire la nostra regola** (e cosa no). La zig-zag ha due cose che Schrödinger non ha: un **cono di
luce rigido** (una pallina non può fare più di una cella per tick) e un **orologio interno naturale**: le
rotazioni (flip) dei cammini a scacchiera. Non ha (per costruzione) niente di nuovo da dire sulle definizioni
standard: quelle sono proprietà dell'equazione di Dirac che la regola ritrova nel limite a → 0 (sarà
VALIDAZIONE). Il punto dove potremmo dire qualcosa di non banale è la domanda 3: **un orologio che conta le
rotazioni** dentro la barriera. Lo definiamo così: ogni cammino della scacchiera che attraversa la barriera fa
un certo numero f di rotazioni dentro la barriera; la "media" sui cammini trasmessi pesata con le ampiezze
(valore debole) è f_w = Σ f A_cammino / Σ A_cammino = ∂ln t/∂ln z, con z peso extra per ogni rotazione dentro.
Si dimostra (identità esatta sul reticolo, verificata a 10⁻⁹) che
  f_w = sinθ cosθ ∂ln t/∂θ_B + sin²θ n_w,
cioè contare le rotazioni equivale a un **orologio di Larmor accoppiato alla massa** (una piccola variazione di
m solo dentro la barriera) invece che al potenziale. Nel continuo questo orologio legge ψ̄ψ invece di ψ†ψ: è il
**tempo proprio** (per una particella libera: T_m = (L/v)(m/E), la dilatazione del tempo). Non abbiamo trovato
in letteratura questo orologio applicato al tunnel (ricerca limitata: 4 ricerche web; "relativistic
corrections to Larmor" sono indicate come lavoro futuro da De Leo 2013).

---

## 2. Cosa abbiamo fatto (in parole semplici)

**La regola** (`lab/zz_core.py`, 1+1 D, nessuna modifica): due tipi di palline, R (corre a destra a c) e L (corre
a sinistra a c). A ogni tick tutte fanno un passo, poi in ogni cella R e L si mescolano un po' (rotazione di
angolo θ = m c a/ħ: è la massa) e prendono la fase e^{−iVa} del potenziale. L'unica variante (`tunnel_core.step`)
è che θ e V possono cambiare da cella a cella: una **barriera** = N celle con V = V0.

**Quando è tunnel.** Fuori dalla barriera l'elettrone ha energia E > mc². Dentro ha ε = E − V0. Se
|ε| < mc², dentro la barriera l'onda non corre: si spegne come e^{−κx}, κ = √(m² − ε²) (in unità ħ = c = 1). È
la "zona di tunnel" di Dirac. Abbiamo scelto E = 2 mc² e V0 = 1,2 … 2,8 (ε da +0,8 a −0,8), e un caso
ultra-relativistico E = V0 = 5 mc².

**Due modi di calcolare, sempre con la stessa regola.**
1. *Stazionario esatto*: la regola ripetuta per sempre su un'onda di quasi-energia fissata si risolve cella per
   cella con una matrice 2×2 (`lattice_scatter`). Dà t (ampiezza trasmessa), r, la densità dentro. Ogni tempo è una
   derivata esatta di ln t: tempo di fase (∂/∂E), Larmor (∂/∂V dentro), rotazioni (∂/∂m dentro). Abbiamo dovuto
   scrivere la soluzione in forma stabile (det = 1 esatto, densità propagata all'indietro), altrimenti per
   barriere spesse (T ~ 10⁻¹⁴) i numeri erano spazzatura (§9).
2. *Pacchetti in tempo reale*: un pacchetto gaussiano (σ = 20 ħ/mc, solo energia positiva, supporto finito
   tagliato a 7σ) corre contro la barriera tick per tick, con **7 copie** insieme: base, Larmor V0 ± δ/2, massa
   θ ± δa/2 (δ = 10⁻⁴), senza barriera, e con un rivelatore assorbente. Gli orologi qui sono **veri**: un grado di
   libertà in più (come lo spin degli atomi di Steinberg) che ruota solo dentro la barriera; si legge l'angolo
   della parte trasmessa.

**I tempi misurati** (tutti dentro la stessa simulazione):

| nome | come si legge nella regola |
|---|---|
| fase τ_φ (Wigner/Hartman) | ∂arg t/∂E + L/v: ritardo del picco |
| dwell τ_D | Σ_tick P(barriera)·a: quanto è stata dentro in media l'onda (trasmessa + riflessa) |
| Larmor τ_y, τ_z | angolo dell'orologio sul potenziale, nel piano e verso l'asse (τ_y − iτ_z = i ∂ln t/∂V) |
| rotazioni τ_my, τ_mz | lo stesso con l'orologio sulla massa (T_m = i ∂ln t/∂m); = contatore delle rotazioni |
| Bohm (surfisti) | in 1-D i surfisti non si incrociano: il surfista q è il punto con probabilità q alla sua destra; tempo passato nella barriera dai trasmessi |
| arrivo | flusso al rivelatore J(t) oppure probabilità assorbita da un rivelatore ideale F(t) |

**Il contatore di rotazioni.** Nella scacchiera di Feynman ogni cammino è una successione di "dritto" (peso cos θ)
e "giro" (peso −i sin θ). Diamo a ogni giro fatto *dentro la barriera* un peso extra z: t(z) = Σ_f A_f z^f. Il
numero medio di giri, pesato con le ampiezze dei cammini trasmessi, è f_w = ∂ln t/∂ln z (un valore debole: può
essere complesso). Identità esatta (derivando i pesi cos θ, sin θ):

  f_w = sin θ cos θ · ∂ln t/∂θ_barriera + sin²θ · n_w ,   n_w = numero (debole) di tick dentro,

verificata sul reticolo a **2,5·10⁻⁹** contando direttamente con z (`flips_direct`). Per a → 0 dà
i f_w → T_m = τ_my − iτ_mz (con errore O(a), visto: −7,900, −7,953, −7,978, −7,989 → −8,0000 al centro del gap,
L = 8). Quindi **contare le rotazioni = orologio di massa**. Per una particella libera T_m = (L/v)(m/E) = il
**tempo proprio** (la dilatazione del tempo): esame E11a, regola 2,30938 contro 2,30940 (8·10⁻⁶).

## 3. Risultati contro la pre-registrazione

Pre-registrazione: `lab/out/tunnel_prereg.json` (commit `2fedf0a`, prima delle corse). Passi reticolari
a = 0,1 / 0,05 / 0,025 / 0,0125 ħ/mc (stazionario) e 0,05 / 0,025 / 0,0125 (pacchetti); estrapolazione
c0 + c1 a + c2 a².

| esame | cosa | numero della regola | riferimento | soglia | esito |
|---|---|---|---|---|---|
| E1 | cono di luce, 6 corse × 7 copie | densità oltre il cono = 0 esatto; fronte 1 cella/tick | c | 10⁻¹⁰ | **PASS** |
| E2 | τ_φ regola (a→0) vs Dirac, 5 V0 × 6 L | errore max 3,9·10⁻⁶ | Dirac continuo | 10⁻³ | **PASS** |
| E3 | saturazione τ_φ(16) vs τ_φ(8) | cambio max 7,6·10⁻⁴ (V0 = 1,2: κ piccola) | satura | 10⁻³ | **PASS** (di poco) |
| E4 | velocità apparente L/τ_φ a L = 4 | 1,54c … 3,98c | > c | > 1 | **PASS** |
| E5 | spostamento del picco trasmesso | PK1 9,0018 / PK2 7,0242 | L − vτ_φ = 9,0000 / 7,0000 | 0,02 L | **PASS** (2·10⁻⁴ L, 3·10⁻³ L) |
| E5 | picco "prima della luce" | PK1 +4,052, PK2 −29,82 | +4,051 / −29,82 | segno | **PASS** |
| E6 | Larmor = dwell, stazionario | 2,3·10⁻⁷ (trasmessi e riflessi) | Büttiker | 10⁻³ | **PASS** |
| E6 | P_T τ_y^T + P_R τ_y^R = dwell, pacchetto | 1,5·10⁻⁵ | identità | 10⁻⁴ | **PASS** |
| E7 | orologio vero (pacchetto) vs stazionario | τ_y 1,8·10⁻⁴, τ_z 9·10⁻⁶ | media spettrale | 10⁻³ / 10⁻² | **PASS** |
| E8 | Larmor regola vs Dirac | 1,7·10⁻⁵ | Dirac | 10⁻³ | **PASS** |
| E9 | limite NR (V0/m → 0), rettangolo κ0L = 3 | 3,4·10⁻⁵ | Büttiker (Schrödinger) | 1 % | **PASS** |
| E10 | barriera gaussiana, V0/m = 0,01 | 0,5 % max; forma come prevista | Schrödinger | 1 % | **PASS** |
| E11a | rotazioni in zona libera = tempo proprio; identità del contatore | 8·10⁻⁶; 2,5·10⁻⁹ | L m/k | 10⁻³; 10⁻⁶ | **PASS** |
| E11b | τ_my/τ_y (L = 16) | = ε/m a 6·10⁻⁷ (5 V0) | ε/m | 1 % | **PASS** |
| E11c | pendenze dτ_mz/dL, dτ_z/dL | m/κ e ε/κ a 10⁻⁴ | formule | 1 % | **PASS** |
| E11d | centro del gap: τ_mz/L | 0,99967 (L=4), 1,000002 (L ≥ 8); τ_my/τ_y < 0,003 | 1 | 1 % | **PASS** |
| E11e | segno di τ_my = segno di ε | sì (L = 8, 12, 16) | — | — | **PASS** |
| E11f | barriere sottili: τ_mz < L/c | 0,8509 L (L = 1), 0,981 L (L = 2) | Dirac 0,8509 | — | **PASS** |
| E12 | surfisti trasmessi ≥ L/c − 2a; decomposizione del dwell | min 74,8 (L = 10), 69,4 (L = 8); 3·10⁻¹¹ | — | 10⁻³ | **PASS** (vedi nota) |
| E13 | rivelatore ideale: probabilità assorbita = probabilità passata | rapporto 0,902 (PK1), 1,078 (PK2) | 1 | 10⁻⁶ | **FAIL** |

Nota su E12: τ_B dei riflessi l'ho ricavato dalla densità (parte della barriera dietro la separatrice), quindi
la decomposizione è quasi un'identità; il controllo vero è traiettorie contro densità per i trasmessi: 4·10⁻⁵.

## 4. I risultati spiegati

### 4.1 Il cono di luce e il fronte
Ogni pallina fa una cella per tick: dopo n tick niente può stare oltre n celle dal bordo del pacchetto di
partenza. Nella regola questo non è "quasi": è **zero esatto** (le celle non sono mai state toccate). E il primo
pezzo che arriva è fatto dal solo cammino che non ruota mai: ampiezza cosⁿθ × (bordo iniziale) × fasi.
**La barriera non può attenuarlo** (cambia solo la fase): è il "precursore" di Sommerfeld–Brillouin, qui visibile
cella per cella (in figura, pannello a, la coda rossa a destra a ~10⁻¹⁷ che arriva fino a x ≈ 235).

### 4.2 Hartman e il picco "superluminale"
Per barriere spesse t ≈ F(ε) e^{−κL} e la fase di F non dipende da L: τ_φ smette di crescere. Valori saturati
della regola (a → 0, E = 2): τ_φ = 2,694 / 1,512 / 1,155 / 1,008 / 1,155 ħ/mc² per ε = +0,8 … −0,8 (Dirac:
uguali a 4·10⁻⁶). A L = 16 il tempo "per attraversare" è ~1 ħ/mc² mentre la luce ci mette 16: velocità
apparente 6–16 c.
Pacchetti: PK1 (E = V0 = 5 mc², v = 0,98 c, L = 10, T = 8·10⁻⁹): il picco trasmesso è **9,002 ħ/mc avanti** al
picco del pacchetto libero, ed esce **4,05 ħ/mc davanti** al raggio di luce partito dal picco in arrivo. PK2 (E =
V0 = 2, v = 0,866): avanti di 7,02 ma 29,8 *dietro* la luce (la particella è troppo lenta perché il vantaggio
basti). **Come si spiega** (e la regola lo fa vedere): i surfisti trasmessi (pannello e) partono tutti dalla coda
davanti del pacchetto, entrano nella barriera presto, ci restano a lungo (sono lentissimi lì dentro) ed escono;
nessuno va mai oltre c. Il picco trasmesso non è "il picco in arrivo che ha attraversato": è una gobba nuova fatta
dalla coda davanti (Winful 2006: accumulo e rilascio; Krekora–Su–Grobe 2001 con Dirac). **Validazione.**

### 4.3 Larmor, dwell, Büttiker, Ramos
* τ_y (Larmor nel piano) = dwell per trasmessi e riflessi, come Büttiker 1983 e, nel caso relativistico, Li–Liang–Kobe:
  2·10⁻⁷ sulla regola.
* Limite non relativistico (V0/mc² = 0,04, 0,02, 0,01 → 0, κ0L = 3): la regola dà τ_y V0/ħ = 0,51365 / 0,81648 /
  1,13478 / 1,48673 per E/V0 = 0,2 / 0,4 / 0,6 / 0,8; Büttiker 0,51367 / 0,81650 / 1,13479 / 1,48673. Anche τ_z
  (Büttiker–Landauer) e τ_φ a ≤ 1,5·10⁻⁵. La correzione relativistica al primo ordine in V0/mc² è misurata
  (pendenze nel JSON `lab/out/tunnel_r3.json`, es. τ_y: +0,14 … +0,45 per unità di V0/mc²).
* **Ramos 2020, solo forma.** Con barriera gaussiana V0 e^{−2x²/w²} e l'orologio proporzionale al profilo (la
  Raman di Ramos la fa la luce della barriera): per κ0w = 2, τ_y V0/ħ = 0,930 → 1,355 → 1,716 → 1,902 per E/V0 =
  0,3 → 0,9: **energia più bassa, tempo più corto**, come la previsione citata per Ramos. Per κ0w = 1 la curva è
  quasi piatta con un massimo (0,818 a E/V0 ≈ 0,5–0,7). In millisecondi per ⁸⁷Rb (ipotesi: w = 1,3 µm come raggio
  1/e², E = 1 nK): τ_y = 4,2 ms (V0 = 1,25 nK) … 0,80 ms (6 nK) … 0,52 ms (8 nK): **0,62 ms esce per V0 ≈ 7 nK**
  (trasmissione ~1 %). Non posso dire se è l'altezza vera: **non è un esame**, solo ordine di grandezza giusto.

### 4.4 L'orologio delle rotazioni (la parte nuova)
**Cosa misura.** L'orologio di Larmor si attacca al potenziale: "conta" ψ†ψ, la densità. L'orologio delle
rotazioni si attacca alla massa: conta ψ̄ψ (densità "scalare"). Per una particella che corre ψ̄ψ/ψ†ψ = m/E = √(1 −
v²/c²): è il ritmo del **tempo proprio**. Quindi le rotazioni della zig-zag sono l'orologio da polso della
particella (Jacobson–Schulman 1984 avevano già legato il numero di svolte della scacchiera al limite non
relativistico; qui lo usiamo come orologio dentro la barriera).

**Cosa succede nel tunnel** (regola a → 0, barriera spessa, tutti pre-registrati e PASS):
1. Nel piano: τ_my = (ε/m)·τ_y **esatto** (6·10⁻⁷). Dentro un'onda che si spegne ψ̄ψ/ψ†ψ = ε/m (non m/E!). Quando la
   barriera è più alta dell'energia totale (ε < 0) l'orologio di rotazione **gira all'indietro**: −0,315 ħ/mc² per
   ε = −0,4 mentre Larmor dice +0,787.
2. Inclinazione: τ_mz cresce con pendenza **m/κ** (Larmor: ε/κ). È il tempo proprio del cammino "a velocità
   immaginaria" v = iκc/ε di Büttiker–Landauer: γ⁻¹ = √(1 − v²) = m/ε, e (Lε/κ)·(m/ε) = Lm/κ. Siccome κ ≤ m,
   **dτ_mz/dL ≥ 1/c**: per barriere spesse l'orologio delle rotazioni non è mai superluminale.
3. Al centro del gap (V0 = E): Büttiker–Landauer dà **0** (pendenza ε/κ = 0: tunnel istantaneo, come
   l'attoclock), Larmor nel piano dà τ_y = v·ħ/mc² (0,866), il tempo di fase satura a 1,155; l'orologio delle
   rotazioni dà **τ_mz = L/c esatto** (0,999997·L per L ≥ 8) e τ_my = 0. *Il tunnel, contato in rotazioni, dura
   esattamente quanto la luce.*
4. Ma non è un limite rigido: per barriere sottili τ_mz = 0,851 L (L = 1) e 0,981 L (L = 2): sotto L/c (previsto da
   Dirac e ritrovato a 4·10⁻⁶).
5. Limite non relativistico: ε → m, ε/m → 1, m/κ → ε/κ: torna Larmor (Büttiker). Le differenze sono effetti
   relativistici di ordine V0/mc², non misurabili con atomi freddi (V0/mc² ~ 10⁻³⁶), misurabili in un
   **simulatore di Dirac**.

**Previsione falsificabile.** In un simulatore quantistico dell'equazione di Dirac 1+1 D dove la "massa" è un
accoppiamento controllabile (ioni intrappolati: Gerritsma et al. hanno simulato Zitterbewegung, *Nature* 463, 68
(2010), e il paradosso di Klein, *PRL* 106, 060503 (2011) — quest'ultimo riferimento da memoria, non verificato;
guide d'onda fotoniche binarie, Dreisow et al. ~2012, non verificato), rendere la massa dipendente da un qubit
"orologio" solo dentro la barriera. Prevediamo: angolo nel piano/angolo dell'orologio di potenziale = (E − V0)/mc²,
cambio di segno a V0 = E; inclinazione per unità di spessore = m/κ (≥ 1/c), uguale a 1/c a V0 = E.

**Quanto è nuovo.** Le formule seguono da Dirac in due righe una volta scelto di derivare rispetto alla massa
dentro la barriera: nessuna fisica nuova. Nuovo (per quel che abbiamo trovato, 4 ricerche web) è: (i) l'identità
"contare le rotazioni della zig-zag = orologio di massa" come lettura di un tempo di tunnel; (ii) le tre proprietà
sopra (ε/m, m/κ ≥ 1/c, L/c al centro). Se esiste in letteratura un "Larmor clock coupled to the mass term" non
l'abbiamo trovato: va cercato meglio prima di dirlo in pubblico.

### 4.5 I surfisti (Bohm)
I surfisti trasmessi sono i più avanti del pacchetto (quantile < T ≈ 10⁻⁸): entrano nella barriera presto e
aspettano lì che passi il grosso del pacchetto. Tempo medio nella barriera: 114,8 (PK1) e 112,2 (PK2) ħ/mc²,
minimo 74,8 e 69,4: **sempre ≥ L/c**, nessuno supera c (con la regola è automatico: un quantile non può spostarsi
più di una cella per tick). **Esplorativo (non pre-registrato)**: per PK2 con σ = 10, 20, 40, 80 il tempo dei
surfisti trasmessi è 58, 112, 217, 420: **circa 5,2 σ + 6** (lineare). Nel tunnel spesso il tempo di Bohm misura il pacchetto,
non la barriera (nel limite stazionario va a τ_D/T → ∞, Leavens). Gli altri orologi invece non cambiano con σ
(τ_y = 0,8684 / 0,8666 / 0,8662 / 0,8661).

### 4.6 Tempo di arrivo
Rivelatore a X_D = L + 20 ħ/mc. Il flusso J(t) attraverso il legame (|R|² − |L|²) è la lettura "senza
disturbo": il pacchetto passato dalla barriera arriva **prima** di quello libero (tempo medio 174,5 contro 183,6
in PK1; 197,4 contro 205,5 in PK2: anticipo = (L − vτ_φ)/v). Backflow (J < 0) presente ma minuscolo (6·10⁻⁷ del
massimo). **Il rivelatore ideale della regola** (togliere tutto quello che entra in x ≥ X_D: è la condizione "solo
onda uscente" L = 0 al bordo, la stessa forma della regola del bordo assorbente di Tumulka per Dirac) **non è
innocuo**: una particella con massa che corre ha sempre un po' di L dentro, e tagliarlo crea un'onda riflessa con
probabilità (1 − v)/(1 + v) (post-hoc: pacchetto 0,07188, stazionario della regola 0,07185, formula 0,07180 per
v = 0,866; 0,01026 per v = 0,98). Con la barriera dietro (che riflette quasi tutto) si forma una **cavità** e il
conteggio cambia di −10 % (PK1) e +8 % (PK2); il tempo medio di arrivo assorbito differisce di −1,8 e +4,4 ħ/mc²
da quello del flusso. Avevo previsto (E13) che il rivelatore ideale contasse esattamente ciò che passa: **FAIL**.
È un limite onesto e istruttivo: anche nella regola il "tempo di arrivo" dipende dal rivelatore, e il rivelatore
perfetto esiste solo per v → c (massa zero), come concludevano Alonso, Sala Mayato e Leavens per la prima-arrivata
sulla scacchiera.

## 5. Validazione, scoperta, limite

* **VALIDAZIONE** (la regola ritrova cose note): cono di luce e precursore; Hartman e saturazione; picco
  "superluminale" solo come rimodellamento; Larmor = dwell; Büttiker nel limite non relativistico; forma di
  Ramos (energie basse → tempi corti per barriere opache); Bohm causale.
* **SCOPERTA (piccola)**: contare le rotazioni della zig-zag = orologio di massa = tempo proprio; nel tunnel spesso
  τ_my/τ_y = ε/m (con segno), dτ_mz/dL = m/κ ≥ 1/c, τ_mz = L/c esatto a V0 = E. Controllata con convergenza in a
  (4 passi), identità del contatore a 10⁻⁹, orologio vero nei pacchetti a 10⁻⁵. Previsione per simulatori di Dirac.
  Non verificata la novità in modo esaustivo.
* **LIMITI**: (1) 1+1 D e una particella sola: niente coppie (V0 ≥ 2mc² in teoria dei campi creerebbe coppie,
  vedi Schwinger/Sauter in `IDEE_RIPRESE.md`); (2) l'attoclock (ionizzazione in 3-D con un laser) è fuori portata
  qui; (3) il tempo di arrivo con rivelatore dipende dal rivelatore (E13 FAIL); (4) il valore di Ramos in ms non è
  verificabile senza l'altezza della barriera; (5) la regola non sceglie quale tempo sia "vero": dà tutti i tempi
  noti più uno nuovo, ognuno risponde a una domanda diversa.

## 6. Debiti (cose messe a mano)

| debito | perché |
|---|---|
| barriera di potenziale rettangolare/gaussiana V(x) data a mano | è l'"esperimento", ma non nasce da cariche vere |
| orologi (V ± δ/2, θ ± δa/2) messi solo dentro la barriera | è l'apparato di misura; nella realtà sarebbe un campo/una Raman |
| rivelatore assorbente (cancellazione delle ampiezze oltre X_D) | apparato di misura idealizzato |
| pacchetto iniziale preparato con gli autovettori di Bloch e tagliato a 7σ | stato iniziale scelto |
| parametri di Ramos (w = 1,3 µm come raggio 1/e², E = 1 nK, V0 scansionata) | non letti dall'articolo |
| massa del ⁸⁷Rb | ammessa (massa di un nucleo/atomo), usata solo per convertire in ms |

**Parametri liberi: 0.** E, V0, L, σ sono condizioni dell'esperimento; δ è infinitesimo (derivata).

## 7. Come rifare tutto (≈ 6 minuti su un core)

```bash
cd lab
export OMP_NUM_THREADS=1
python3 tunnel_core.py                 # controlli: t = 1 senza barriera, T + R = 1, convergenza in a
python3 tunnel_prereg.py               # pre-registrazione (già committata) -> out/tunnel_prereg.json
python3 tunnel_r1_stat.py              # R1 stazionario esatto (1 s)  -> out/tunnel_r1.json
for c in PK1 PK2; do for a in 0.05 0.025 0.0125; do python3 tunnel_r2_pack.py $c $a; done; done   # R2 (~2,5 min)
python3 tunnel_r2_verdict.py           # giudizio R2 -> out/tunnel_r2_verdict.json
for s in 10 40 80; do python3 tunnel_r2_pack.py PK2 0.05 $s; done   # esplorativo (larghezze)
python3 tunnel_r3_nr.py                # R3 limite NR, Ramos, assorbitore (15 s) -> out/tunnel_r3.json
python3 tunnel_figure.py               # figura
```
I JSON dei pacchetti a a = 0,05 e 0,025 e quelli esplorativi completi non sono nel repo (300 kB–1 MB l'uno): si
rigenerano; `tunnel_r2x_summary.json` riassume gli esplorativi.

## 8. Errori, cambi e cose post-hoc (dichiarati)

1. **Prima versione del calcolo stazionario instabile** per barriere spesse (cancellazione tra esponenziali e^{±κL}
   e det calcolato numericamente): corretta *prima* della pre-registrazione (det M = 1 esatto, densità propagata
   all'indietro). Le previsioni del continuo nel JSON sono quelle corrette.
2. **Prima esecuzione di E10 con un bug**: l'orologio sulla barriera gaussiana era un campo uniforme su ±3w invece
   che proporzionale al profilo (come il riferimento di Schrödinger e come in Ramos): valori ×10 sbagliati.
   Corretto e rieseguito; la definizione dell'esame non è cambiata.
3. **E10 ritoccato prima delle corse** (nello stesso commit della pre-registrazione): la prima stesura diceva
   "τ_y cresce con E/V0" per tutte le larghezze; il calcolo di Schrödinger fatto per scrivere le previsioni
   mostrava che per κ0w = 1 non è monotono: la previsione registrata è quella giusta.
4. **E12**: il tempo dei riflessi è ricavato dalla densità, quindi la decomposizione del dwell è quasi
   un'identità; il controllo vero (traiettorie contro densità) è riportato.
5. **Post-hoc / esplorativo, non pre-registrato**: dipendenza dalla larghezza σ (§4.5), riflessione del
   rivelatore (1 − v)/(1 + v) (§4.6), valori di Ramos in ms.
6. Non ho potuto leggere il testo completo di Ramos et al. 2020 (limite di accesso del proxy).

## 9. Unità fisiche (per l'elettrone)
ħ/mc = 3,86·10⁻¹³ m, ħ/mc² = 1,29·10⁻²¹ s. Una barriera di L = 10 ħ/mc è 3,9 pm; il tempo di fase saturato
(~1 ħ/mc²) è ~1,3 zeptosecondi; la luce ci mette 13 zs. Barriere così (V0 ~ MeV su pm) esistono solo vicino ai
nuclei: per l'elettrone il caso pratico è il limite non relativistico (§4.3), dove tutti gli orologi della
regola tornano quelli di Büttiker.
