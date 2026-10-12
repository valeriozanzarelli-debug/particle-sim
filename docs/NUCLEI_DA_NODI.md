# Nuclei dai nodi: quando i nodi del fluido di palline si mettono insieme

**L'idea (di Valerio).** Le palline si aggregano in atomi, e gli atomi hanno proprietà che le palline da sole non
hanno. Allo stesso modo i **nodi** del fluido di palline (le particelle di `NODO_3D.md`) dovrebbero aggregarsi in
particelle composte, con **forme nuove** e proprietà che non sono "la somma dei pezzi".

Nel modello di Skyrme questo è letterale. Il nodo con B = 1 è il nucleone (protone/neutrone). Mettendo insieme 2, 3, 4
nodi, il campo di palline si risistema da solo in forme precise e sorprendenti:

| B | nucleo "corrispondente" | forma (letteratura) |
|---|---|---|
| 1 | nucleone | sfera (riccio) |
| 2 | deuterone | **ciambella** (toro) |
| 3 | ³He / ³H | **tetraedro** (gabbia a 4 facce) |
| 4 | ⁴He | **cubo** (gabbia a 6 facce) |

(Braaten–Townsend–Carson 1990; Battye–Sutcliffe 1997 e 2002; Houghton–Manton–Sutcliffe 1998.)
Qui lo abbiamo **rifatto sul nostro reticolo di palline**, controllato con esami scritti prima, e abbiamo guardato
anche se le forme **nascono da sole** partendo da nodi B = 1 separati.

Pre-registrazione: `lab/nuclei_prereg.py` → `lab/out/nuclei_prereg.json`, **committata (cd859e2) prima delle corse**.
Le prove di sviluppo viste prima (mappe razionali e campi iniziali non rilassati) sono dichiarate nel suo docstring
(`lab/nuclei_dev.py` → `lab/out/nuclei_dev.json`).
Codice: `lab/nuclei_core.py` (mappe razionali, misure di forma), `lab/nuclei_run.py` (rilassamento a pezzi da 8 min
con checkpoint), `lab/nuclei_summary.py` (esami), `lab/nuclei_figure.py`. Motore: `knot_core.py` e `knot2_core.py`
(schema c4) usati in sola lettura. Numeri: `lab/out/nuclei_*.json`, log in `lab/out/nuclei_logs/`.
Figura: `lab/out/nuclei_figure.png`.

## Come funziona, in breve

* **Stesso modello di `NODO_3D.md`**: ogni pallina porta una freccia a 4 componenti di lunghezza 1; energia di Skyrme;
  pioni senza massa. Unità di Skyrme. Reticolo 48³ palline, passo h = 0.2 (mezza scatola 4.7), bordi = vuoto.
  Derivate dello schema "c4" (4° ordine sfalsato, `knot2_core.py`): a h = 0.2 il nodo B = 1 viene E/12π² = 1.2272,
  solo −0.4 % dal valore esatto 1.2315.
* **Rilassare** = far scendere l'energia finché non si ferma (L-BFGS). Non è una dinamica nel tempo: è "lasciare che il
  campo si sistemi nella forma di energia minima vicina".
* **Punto di partenza "mappe razionali"** (HMS 1998): un trucco per scrivere un campo con B avvolgimenti e la simmetria
  giusta. Abbiamo verificato le mappe: danno B = 1, 2, 3, 4 esatti e gli integrali angolari I = 1, 5.808, 13.577,
  20.650 (letteratura 1, 5.81, 13.58, 20.65).
* **Misure di forma**: densità barionica b(x) (dove "sta" il numero di nodi); punteggio di simmetria
  asym(g) = ‖b − g·b‖/‖b‖ (0 = perfettamente simmetrico sotto la rotazione g); multipoli q_l (q₂ ≠ 0 = schiacciato o
  allungato, q₃ ≠ 0 = forma "a tetraedro", q₄ = forma "a cubo"); tensore dei momenti secondi; superficie della densità
  al 25 % del massimo e il suo numero di buchi.

**Parametri.** 0 parametri liberi per forme ed E/B (in unità di Skyrme, pioni senza massa). Solo per i MeV servono le
2 misure di ANW (M_N, M_Δ → F_π = 129 MeV, e = 5.44): unità di **fit**, dichiarate.

## Tabellone

| prova | cosa | esito |
|---|---|---|
| N1a B=2 | resta una ciambella perfetta | **PASS** |
| N1a B=3 | resta un tetraedro (non un cubo) | **PASS** |
| N1a B=4 | resta un cubo | **PASS** |
| N1b | E/12π²B entro 2 % dalla letteratura (B=2,3,4) | **PASS** (−0.26 %, −0.29 %, −0.30 %) |
| N1c | E/B scende con B | **PASS** |
| N1d | B intero entro 1 % | **PASS** (1.9949, 2.9932, 3.9926) |
| N2 B=2 | da due nodi B=1 nasce la ciambella da sola | **PASS** (distanza 2.5 e 3.0) |
| N2 previsione asse | asse della ciambella lungo y | **FAIL** (viene lungo z) |
| N2 B=4 tetraedro | da quattro nodi B=1 nasce il cubo | **PASS** (fermato per tempo, quasi a convergenza) |
| N2 B=4 quadrato | da quattro nodi B=1 nasce il cubo | **PASS** |
| N3a | B=2 schiacciato, somma ingenua allungata | **PASS** |
| N3b | B=2 ha un anello, somma ingenua ha 2 picchi | **PASS** |
| N3c | legame del deuterone entro un fattore 2 | **FAIL (atteso)**: 72 MeV contro 2.2 MeV, **32 volte** |
| N4 | superfici: sfera (χ = 2) e ciambella (χ = 0) | **PASS** (B=3: 3 buchi, B=4: 5 buchi) |

## N1 — Le forme dei nuclei di Skyrme (partendo dalla forma giusta)

| B | E/12π²B reticolo | letteratura (BS 2002) | scarto | B misurato | raggio barionico (lett. BS 1997) |
|---|---|---|---|---|---|
| 1 | 1.2272 | 1.2322 | −0.41 % | 0.9965 | 1.00 (1.03) |
| 2 | 1.1760 | 1.1791 | −0.26 % | 1.9949 | 1.37 (1.42) |
| 3 | 1.1429 | 1.1462 | −0.29 % | 2.9932 | 1.61 (1.64) |
| 4 | 1.1167 | 1.1201 | −0.30 % | 3.9926 | 1.80 (1.86) |

* L'energia per nucleone **scende**: 1.227 → 1.176 → 1.143 → 1.117. Mettere insieme nodi conviene: è il "legame".
* I rapporti sullo stesso reticolo (E_B/B)/E₁ = 0.9583, 0.9313, 0.9100 contro la letteratura 0.9569, 0.9302, 0.9090:
  d'accordo al 0.15 %. Gli errori del reticolo (≈ −0.3 %, e la scatola un po' stretta) si cancellano quasi tutti.
* Il punto di partenza (mappa razionale) aveva E/12π²B = 1.219, 1.190, 1.139: rilassando scende del 2–4 %, cioè la
  mappa razionale è una buona approssimazione ma **non** è la forma esatta; il campo la corregge da solo.
* Forme: la ciambella è simmetrica sotto rotazioni di 30°, 45°, 90° attorno al suo asse (asym ≤ 0.003) e ha un buco
  (densità al centro / massimo = 0.007); il tetraedro ha le 3 rotazioni di 180° e quella di 120° (asym ~ 5·10⁻⁶) ma
  **non** l'inversione (asym 0.88): è un tetraedro, non un cubo; il cubo ha le rotazioni di 90°, di 120° e l'inversione.
  Multipoli: q₂ = 0.25 (B=2, schiacciato), q₃ = 0.20 (B=3), q₄ = 0.147 (B=4) e tutti gli altri ~0.
* **Onestà:** in N1 la simmetria **è messa dentro** dal punto di partenza, e il reticolo cubico la conserva. Quindi N1
  verifica che il rilassamento **non rovina** la forma e che l'energia è giusta; **non** che la forma nasca da sola.
  Per quello c'è N2.

## N2 — La forma nasce da sola? Sì

**Due nodi → ciambella.** Due ricci B = 1 a distanza 2.5, il secondo girato di 180° nello spazio dell'isospin (attorno
all'asse z: è il "canale attrattivo", in cui i due nodi si attirano). All'inizio la densità è chiaramente **due palle**
(2 picchi, forma allungata: momenti 0.37, 0.51, 1.85; E/12π²B = 1.316). Rilassando:

* E/12π²B scende a **1.175997**, uguale alla ciambella di N1 fino al **10⁻⁸** (relativo);
* i due picchi spariscono: la densità diventa un **anello** (centro/anello = 0.007, variazione lungo l'anello 0.2 %),
  simmetrico sotto rotazioni di 30°, 45°, 90° attorno al suo asse (asym ≤ 0.003);
* la superficie al 25 % ha χ = 0: **un buco**, una ciambella.

**Questo è il punto dell'idea di Valerio**: due "palle" messe vicine non restano due palle attaccate. Il fluido si
risistema in un oggetto nuovo, con una forma (un anello) che nessuno dei due pezzi aveva.

**Previsione sbagliata (FAIL dichiarato).** Avevamo scritto prima che l'asse della ciambella sarebbe venuto lungo y
(perpendicolare sia alla separazione x sia all'asse z della rotazione). È venuto lungo **z** (|asse_z| = 0.99999): la
ciambella giace nel piano che contiene la linea che univa i due nodi, e il suo asse è l'asse della rotazione relativa.
Il fatto (la ciambella si forma) è giusto; la nostra intuizione sulla sua orientazione era sbagliata.

**Quattro nodi → cubo, da due partenze diverse.**

| partenza | E/12π²B inizio → fine | q₂, q₃, q₄ (cubo di N1: 0, 0, 0.147) | superficie | converge? |
|---|---|---|---|---|
| tetraedro (lato 2.5, rotazioni 1, π_x, π_y, π_z: tutte le coppie attrattive) | 1.355 → **1.11615** | 0.001, 0.0001, 0.148 | χ = −8 (5 buchi) | no: fermato per tempo a 126 iterazioni, l'energia cambiava ancora di 10⁻⁶ |
| quadrato (lato 2.5, rotazioni 1, π_z, 1, π_z) | 1.353 → **1.11642** | 0.0006, 0.00001, 0.148 | χ = −8 (5 buchi) | sì (46 iterazioni) |

In tutti e due i casi **nasce il cubo**: stessa energia del cubo di N1 entro lo 0.05 %, stessi multipoli, 8 picchi
(i vertici). I cubi nati da soli sono **ruotati** rispetto al reticolo (non hanno le rotazioni di 90° degli assi
x, y, z ma hanno l'inversione), e hanno un'energia sul reticolo un po' **più bassa** (−0.05 % e −0.03 %) di quello
allineato agli assi di N1: un piccolo effetto del reticolo (l'energia dipende un poco da come la forma è orientata
rispetto alle palline). Il quadrato passa per una fase intermedia (E/12π²B ≈ 1.14 dopo 10 iterazioni) prima di
"chiudersi" nel cubo.

**Controlli in più (pre-registrati "se c'è tempo", fatti).**
* Due nodi a distanza **3.0** invece di 2.5: stessa storia, ciambella con E/12π²B = 1.175997 (uguale a 10⁻⁸), asse
  di nuovo lungo z. La ciambella non dipende dalla distanza di partenza (in questo intervallo).
* Cubo in una scatola più grande (64³, mezza scatola 6.3): E/12π²B = **1.11434** (−0.21 % rispetto a 48³; −0.51 % dalla
  letteratura, sempre dentro il 2 %), B = 3.9929, raggio 1.82 (lett. 1.86). La scatola 48³ stringe un po' la coda del
  campo: la parte del raggio che manca (−3 %) è soprattutto questo.

## N3 — Il composto non è la somma dei pezzi

**Forma.** Prendiamo due nodi B = 1 veri (rilassati) e li mettiamo uno accanto all'altro alla distanza (1.88) che dà
**lo stesso raggio** del B = 2 vero (1.374). Stessa taglia, ma:

| | B = 2 vero | due B = 1 affiancati |
|---|---|---|
| momenti secondi (dal più piccolo) | 0.317, 0.785, 0.785 → **schiacciato** (disco) | 0.334, 0.334, 1.219 → **allungato** (sigaro) |
| picchi della densità | anello (8 massimi tutti uguali, è il reticolo) | **2** picchi |
| q₂ | 0.25 | 0.47 |

**Taglia.** I raggi barionici 1.00, 1.37, 1.61, 1.80 crescono più in fretta di B^(1/3) (che darebbe 1.26, 1.44, 1.59
per una goccia compatta): i nuclei di Skyrme piccoli sono **gusci vuoti** (anello, gabbie), non palline piene.
Momenti d'inerzia (della densità barionica): B = 2 ha due valori uguali e uno diverso (1.10, 1.10, 1.57: una trottola
simmetrica), B = 3 e B = 4 hanno tre valori uguali (1.73 ×3 e 2.16 ×3: sfere "per l'inerzia" anche se non sono sfere, perché
tetraedro e cubo hanno troppe simmetrie per avere un asse privilegiato).

**Energia di legame (il FALLIMENTO noto).** Legame = B·E₁ − E_B. In unità di ANW (F_π = 129 MeV, e = 5.44, **fit**;
1 unità di 12π² = 702 MeV, E₁ = 862 MeV):

| nucleo | legame Skyrme (reticolo) | legame Skyrme (letteratura) | vero | troppo forte di |
|---|---|---|---|---|
| deuterone (B=2) | 72 MeV | 75 MeV | 2.22 MeV | **×32** |
| ³He/³H (B=3) | 178 MeV | 181 MeV | 7.7 / 8.5 MeV | ×23 |
| ⁴He (B=4) | 310 MeV | 315 MeV | 28.3 MeV | ×11 |

Il modello di Skyrme classico lega i nodi **molto troppo**: il 4–9 % della massa invece dello 0.1–0.75 % vero.
Per ⁴He è il famoso "fattore ~10"; per il deuterone è peggio (×32). È un fallimento noto da decenni; i rimedi studiati in
letteratura (massa del pione, modelli "quasi-BPS", energia di vibrazione e rotazione quantistica, mesoni ρ/ω) non li
abbiamo provati qui.

## N4 — Disegnare in negativo

Solo la superficie dove la densità barionica vale il 25 % del suo massimo, niente palline (pannelli in alto nella
figura; in negativo vuol dire: si disegna il nodo, non le palline). Il numero di buchi si legge dalla caratteristica di Eulero χ = V − E + F della superficie (una superficie chiusa
con g buchi ha χ = 2 − 2g):

| | B=1 | B=2 | B=3 | B=4 | B=2 da 2 nodi | B=4 da tetraedro / quadrato |
|---|---|---|---|---|---|---|
| χ | 2 | 0 | −4 | −8 | 0 (anche d = 3) | −8 / −8 |
| buchi | 0 (sfera) | 1 (ciambella) | 3 | 5 | 1 | 5 / 5 |
| triangoli | 1 896 | 5 952 | 9 804 | 15 168 | 5 952 | 14 836 |

Il tetraedro è una gabbia con 4 finestre (genere 3), il cubo una gabbia con 6 finestre (genere 5): la topologia della
superficie conta le "facce bucate" del poliedro.

## Debiti e cose che mancano

* **Il legame è sbagliato di 10–30 volte** (N3c). Le forme sono belle, l'energia di legame no.
* **Solo forme classiche.** Spin, isospin, il fatto che il deuterone abbia spin 1 e ⁴He spin 0: servono la
  quantizzazione dei modi collettivi (letteratura: Braaten–Carson per il deuterone, Walet, Manton e altri). Non fatto.
* **Il deuterone vero è grande e debole** (≈ 2 fm, due nucleoni quasi separati); la ciambella è compatta. Anche questo
  è una differenza nota dal mondo reale.
* **Pioni senza massa.** Con la massa del pione (3° parametro) per B più grandi le forme cambiano (catene di cubi, ecc.);
  fino a B = 4 restano queste.
* **N2 è una discesa di energia, non una collisione nel tempo.** Mostra che la forma è il minimo vicino alle partenze,
  non come ci si arriva davvero urtando. Le collisioni vere (due nodi che si scontrano e diffondono a 90°) sono un
  esperimento futuro (servono la dinamica di `knot2_core` e scatole più grandi).
* **Reticolo e scatola.** h = 0.2, mezza scatola 4.7: energie −0.3 % rispetto al continuo, B −0.2 %, raggi −3 % (la
  coda che esce dalla scatola). I rapporti fra B diversi sono giusti al 0.15 %.
* **tet4 fermato per budget** a 126 iterazioni (6 pezzi da 8 minuti), energia ferma a 10⁻⁶: il verdetto "cubo" non
  dipende da questo, ma è dichiarato. Incidente tecnico (dichiarato): per una ventina di minuti due "driver" hanno
  fatto girare in parallelo lo stesso pezzo di tet4 (stesso punto di partenza, stesso checkpoint); un pezzo è stato
  buttato. Nessun effetto sui numeri, solo tempo perso.
* **Previsione sull'asse della ciambella sbagliata** (vedi N2).
* **Parametri:** 0 per forme ed E/B; 2 ingressi (M_N, M_Δ) per i MeV.
