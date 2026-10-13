# Motore ottimizzato: quale livello simulare, dove, e quando (lato fisico)

**La domanda di Valerio.** "Un ottimo sistema di ottimizzazione grafica come Unreal Engine", con idee geniali, visive e
geometriche. Nel progetto le impostazioni grafiche **sono** i livelli di agglomerato (`VISIONE_SISTEMA.md`). E poi il
principio, detto con la sua immagine: **una mela in 4K e una in 480p devono marcire tutte e due.** Sulla Intel UHD del
portatile si accetta meno risoluzione e più errore, ma la **dinamica** deve restare la stessa.

**Perimetro (cambiato dal coordinatore a metà lavoro).** Qui c'è solo il **lato fisico** dello schema multilivello:
1. quale livello di agglomerato usa ogni regione, entro un limite di errore fisico (anche **proiettato sullo schermo**,
   come criterio di scelta, alla Nanite);
2. come si passa da un livello all'altro **senza scatti**, conservando carica ed energia;
3. il **passo di tempo locale** (fino a "orologio fermo" per le regioni che nessuno guarda);
4. che cosa vuol dire, in modo misurabile, **"dinamica conservata"** ai livelli grossi, con esami Ultra contro Intel.

Il disegno (splat gaussiani, preset grafici, regola sulla GPU) è del thread "Modello standard a palline": le idee che
avevo le trovi in fondo, nell'appendice da passare a loro. La GTX 1650 Ti è rotta: ogni budget qui è sulla **Intel UHD**
(e su una scheda da 16 GB come obiettivo).

File: motore `lab/motore_core.py`; pre-registrazione `lab/motore_prereg.py` → `lab/out/motore_prereg.json`
(**committata in e540d7f prima delle corse**); corse `lab/motore_m1.py` (A, B1, B2), `lab/motore_m2.py` (C),
`lab/motore_m3.py` (D), `lab/motore_m4.py` (EL, EE, giudizio); emendamento dichiarato `lab/motore_amend.py`;
figura `lab/motore_figure.py` → `lab/out/motore_figure.png`.

![figura](../lab/out/motore_figure.png)

---

## 0. In breve

| esame | cosa | risultato | esito |
|---|---|---|---|
| A1 | il modello d'errore (dispersione di Bloch della regola, pesata sullo spettro del pacchetto) prevede deriva e fase della luce contro Maxwell vero | 16 casi su 16 (deriva entro 5 % + 0,01 λ, fase entro 0,05 rad) | **PASS** |
| A2 | lo stesso con un solo k₀ (il modello economico che userebbe il motore) | 16/16 | **PASS** |
| B1 | scena "errore sullo schermo": 8 pacchetti da 0,1 a 20 λ per pixel, 100 λ di volo | errore misurato ≤ **0,43 px** su tutti i pacchetti non segnalati; costo **÷28** (Intel) e ÷8 (Ultra) rispetto a tutto a 64 celle/λ; ma i 2 pacchetti segnalati misurati sotto 0,5 px per un difetto della **mia misura** (scatola troppo piccola, vedi §3) | **FAIL** (sotto-criterio; i punti principali passano) |
| B2 | l'errore della fisica **ha memoria**: camera che si avvicina | ingenua (alla Nanite) **14,7 px** di errore finale; col "passaporto d'errore" **0,26 px** | **PASS** (la ingenua doveva fallire, e fallisce) |
| C1 | congelare tutto fuori dal **cono di luce all'indietro** della regione guardata | regione guardata identica **al bit** (luce 2-D e Dirac 2-D con nucleo Z = 60); con un cono più stretto di 1 cella errore 0,58–1,45 | **PASS** |
| C2 | "orologio pigro": recuperare 1024 tick con un salto esatto in Fourier | differenza 1,9·10⁻¹³; 38 volte più veloce dei tick uno per uno | **PASS** |
| C3 | la regola agglomerata imparata ha un cono? | quasi: coda fuori dal cono fisico 10⁻⁴–10⁻² per una sola cella accesa, fuori da 2× lo stencil ≤ 6·10⁻⁶ dopo un passo e 10⁻¹⁵ dopo 4 | descrittivo |
| D1 | cambio di livello della luce (N → N/2 → N) | **spettrale**: energia conservata a 10⁻¹¹, nessun salto; cubico (quello di oggi nel reticolo) 3·10⁻⁵…9·10⁻³; media a blocchi 2·10⁻³…4·10⁻² | **PASS** |
| D2 | cambio di livello dell'elettrone (fine → livello imparato → fine) | blocco dell'agglomeratore + prolungamento spettrale (che "disfa" il blocco): carica e energia cinetica tornano a 10⁻¹⁵; prolungamento a gradini: energia cinetica **+105 %…+3800 %** | **PASS** |
| EL | dinamica della luce, Ultra (32 celle/λ) / Intel (12) / sotto il pavimento (8), contro Maxwell | Ultra **PASS**, Intel **PASS** (diffrazione −9,3 %: al limite del 10 %), sotto il pavimento **FAIL** come previsto (diffrazione −21 %, cucitura T = 0,983, anisotropia 8,5 %) | come pre-registrato |
| EE | dinamica dell'elettrone a scala atomica sui livelli imparati (celle a₀/256, a₀/21, a₀/11) | velocità, allargamento, energie cinetiche entro 0,4 %, norma 10⁻¹³ a **tutti e tre** i livelli | **PASS** (3/3) |

**Il colpo d'occhio.** Le tre idee che reggono il motore sono uscite tutte **esatte o prevedibili**:
il **cono di luce** permette di non calcolare niente fuori da quello che potrà influenzare ciò che guardi (esatto al
bit, perché la pallina zig-zag non va mai più veloce di una cella per tick); il **cambio di livello spettrale**
conserva energia e carica a 10⁻¹¹; e l'**errore di ogni livello si calcola in anticipo** dalla regola stessa (16 casi su
16). La lezione più importante è un **no**: nella fisica, a differenza di Nanite, l'errore **si accumula**. Scegliere il
livello solo in base a quanto è lontano adesso porta a 14,7 pixel di errore quando ti avvicini. Serve un "passaporto"
dell'errore per ogni regione, oppure (per la luce libera) il salto esatto nel tempo.

---

## 1. L'idea in parole semplici

Nanite (Unreal) disegna ogni pezzo di montagna col numero di triangoli che basta perché l'errore **sullo schermo** stia
sotto un pixel: lontano pochi triangoli, vicino tanti. Noi vogliamo la stessa cosa ma per la **fisica**: lontano celle
grosse (livelli di agglomerato alti), vicino celle fini, e l'errore detto in pixel.

Ma ci sono due differenze di fondo, e lo studio gira attorno a queste:

1. **Una montagna non si muove; la fisica sì.** Se simuli un fotone con celle grosse, sbaglia un po' la velocità; dopo
   100 lunghezze d'onda ha accumulato un ritardo. Quando la camera si avvicina, il ritardo **è ancora lì**: non sparisce
   passando a celle fini. L'errore ha memoria (esame B2).
2. **Sotto un certo livello la fisica cambia, non si sfoca soltanto.** Un'onda con meno di ~12 celle per lunghezza
   d'onda non diventa "più sfocata": cambia **diffrazione** (−21 % a 8 celle/λ), **velocità di sbieco** (−8,5 %) e
   attraversa male le cuciture. È il **pavimento dinamico**: l'hardware debole può scendere fino lì, non sotto
   (esami EL). Per non scendere sotto, sulla Intel si **rimpicciolisce la scena**, non la risoluzione.

## 2. Il modello d'errore (A1, A2): sapere in anticipo quanto sbaglia un livello

La regola zig-zag della luce su un livello con N celle per λ ha una dispersione nota esattamente (la matrice di Bloch
di due tick, `zz_core.bloch_batch`). Da lì escono, **senza simulare**:

* errore di frequenza δω/ω (la fase che si accumula: le frange scivolano di δφ λ/2π);
* errore della velocità di gruppo (l'inviluppo arriva in ritardo, e di sbieco devia: "walk-off");
* fattore di diffrazione (quanto si apre il fascio rispetto alla luce vera).

| celle per λ | velocità a 45° | frequenza a 45° | diffrazione a 0° | diffrazione a 45° |
|---|---|---|---|---|
| 8 | −8,43 % | −2,71 % | 0,785 | 1,135 |
| 12 | −3,57 % | −1,17 % | 0,907 | 1,058 |
| 16 | −1,97 % | −0,65 % | 0,948 | 1,033 |
| 32 | −0,48 % | −0,16 % | 0,987 | 1,008 |
| 64 | −0,12 % | −0,04 % | 0,997 | 1,002 |

(Lungo gli assi la luce è un puro spostamento: errore zero di velocità e frequenza, ma la diffrazione è il fattore
ka/tan ka già noto. L'errore scala come N⁻².)

**A1/A2.** Ho fatto volare 16 pacchetti (N = 8, 12, 16, 32 × angoli 0°, 22,5°, 30°, 45°) per 50 λ con la regola vera
(evoluzione esatta in Fourier) e con **Maxwell continuo** (la luce vera, l'esame), dagli stessi dati iniziali.
La deriva misurata cade sulla retta della previsione per 4 ordini di grandezza (figura, pannello a). 16/16 sia col
modello pesato sullo spettro, sia con quello a un solo k₀.

**Una lezione di sviluppo (dichiarata).** Nella prima prova la deriva di traverso a 30° era sbagliata del 45 %: non era
il modello, era la **scatola periodica troppo piccola** (le code del pacchetto rientravano dall'altro lato). Con ±4,5 σ
torna tutto. Questa stessa trappola mi ha poi morso in B1 (§3).

## 3. La scelta del livello con l'errore sullo schermo (B1)

**La regola del motore.** Per ogni regione: impronta del pixel p (quante λ copre un pixel a quella distanza), orizzonte
T (quanto a lungo vivrà la cosa), pavimento del livello (32 per Ultra, 12 per Intel). L'errore visibile in pixel è

* **inviluppo**: deriva / p — sempre;
* **frange**: (δφ/2π)·λ / p — **solo se le frange sono risolte** sullo schermo (almeno 2 pixel per λ); se un pixel
  copre più di mezza λ, il pixel fa la media delle frange e la fase non si vede.

Si sceglie il livello più grosso con errore previsto ≤ 0,5 px (margine 2). Se nessun livello della scala basta, la
regione viene **segnalata**.

| impronta [λ/px] | 0,1 | 0,2 | 0,5 | 1 | 2 | 5 | 10 | 20 |
|---|---|---|---|---|---|---|---|---|
| angolo | 30° | 30° | 0° | 22,5° | 45° | 60° | 15° | 0° |
| livello Intel | segn. | segn. | 12 | 32 | 24 | 16 | 12 | 12 |
| errore misurato Intel [px] | — | — | 0,013 | 0,395 | 0,431 | 0,365 | 0,214 | 0,000 |
| livello Ultra | segn. | segn. | 32 | 32 | 32 | 32 | 32 | 32 |
| errore misurato Ultra [px] | — | — | 0,002 | 0,395 | 0,241 | 0,091 | 0,029 | 0,000 |

* **Errore ≤ 1 px dappertutto**: massimo 0,43 px; previsione e misura coincidono entro il 3 %.
* **Costo** (celle × tick, in 2-D ∝ N³) sui pacchetti non segnalati: **÷28** su Intel rispetto a tenerli tutti a 64
  celle/λ (soglia ÷4: PASS), ÷8 su Ultra (il pavimento di Ultra a 32 costa).
* Notare il pacchetto a 1 λ/px a 22,5°: anche su Intel vuole 32 celle/λ. Le frange non si vedono, ma in 100 λ
  l'inviluppo di sbieco si sposta più di mezzo pixel. **L'angolo conta quanto la distanza.**

**Il FAIL (dichiarato).** Il pre-registro chiedeva anche che i due pacchetti segnalati, misurati al livello più fine
(64), sbagliassero davvero più di 0,5 px. Misurati: 0,605 e **0,302** px. Il secondo "smentiva" la segnalazione. Ma per
non riempire la memoria li avevo fatti con σ = 2 λ, e in 100 λ un pacchetto così stretto si apre per diffrazione
fino al bordo della scatola (±9 λ): le code rientrano e falsano il centro, lo stesso difetto di §2. **Emendamento A1**
(`lab/motore_amend.py`, dopo la corsa, dichiarato): stessi due casi in una scatola di ±7 σ → **1,05 px e 0,52 px**:
la segnalazione era giusta. Il verdetto di B1 resta FAIL perché così era scritto; il motore non ha sbagliato, la mia
misura sì.

## 4. L'errore ha memoria (B2): il passaporto

La camera parte lontana (5 λ/px) e arriva vicina (0,1 λ/px) in 100 λ di volo, pacchetto a 30°.

* **Politica ingenua** (Nanite puro: ogni 10 λ sceglie il livello che basta **adesso**): usa 12 celle/λ quasi tutto il
  tempo (da lontano basta), passa a 16 alla fine. Quando la camera è vicina il pacchetto ha già 1,5 λ di ritardo e
  5,2 rad di fase accumulati: **14,7 px** di errore. Cambiare livello alla fine non serve: l'errore è già nella fisica.
* **Passaporto d'errore**: ogni regione porta con sé l'errore **accumulato** (previsto dal modello di §2) e il motore
  sceglie il livello pensando alla **impronta più piccola che la camera potrà raggiungere** (la camera ha una velocità
  massima: anche lei ha un "cono"). Risultato: 128 poi 96 celle/λ, errore finale **0,26 px** (previsto: fase 0,080,
  misurata 0,081 rad).
* **Il prezzo**: il passaporto costa 3,3·10¹⁰ celle×tick contro 6,4·10⁷ dell'ingenua, **500 volte di più**. È il prezzo
  vero di vedere da vicino una cosa che ha viaggiato a lungo.

**Come non pagarlo (per la luce libera).** Nel vuoto la regola è lineare e uguale dappertutto: lo stato dopo m tick si
calcola **con un salto solo** (W(k)^m in Fourier, esame C2: identico ai tick uno per uno a 2·10⁻¹³). Quindi un fotone che
nessuno guarda **non si simula affatto**: si tiene il suo fotogramma chiave e, quando la camera arriva, si salta
direttamente al presente al livello fine. Costo: una trasformata, invece di 10¹⁰ celle×tick. Per la **materia che
interagisce** il salto non esiste: lì il passaporto dice "rifai da un fotogramma chiave al livello fine" (il costo si
rimanda, non si toglie).

## 5. Il tempo locale (C1, C2, C3)

**C1: il cono di luce come "frustum" causale.** La regola zig-zag sposta ogni pallina di **al massimo una cella per
asse per tick**. Quindi quello che vedrai fra T tick nella finestra guardata dipende solo dalle celle entro distanza T
(nella metrica "a scacchiera", L∞). Tutto il resto si può **congelare**. Prova: mondo 192², finestra 48², 64 tick, campo
dappertutto (tre pacchetti + fondo casuale), sia luce sia elettrone di Dirac con nucleo Z = 60:

| | errore nella finestra | con cono 1 cella più stretto |
|---|---|---|
| luce 2-D | **0 (identico al bit)** | 0,58 |
| Dirac 2-D + Coulomb | **0 (identico al bit)** | 1,45 |

Il cono è **esatto e minimo**: una cella in meno e il risultato è sbagliato. Celle aggiornate: il 37 % (risparmio ×2,7
su 64 tick); ×12 su 8 tick, ×5,6 su 32 (in 3-D i risparmi vanno al cubo invece che al quadrato).

**C2: l'orologio pigro.** Una regione congelata accumula un "debito di tick". Per la luce nel vuoto lo si paga con un
salto in Fourier: 1024 tick in 0,09 s invece di 3,4 s (n = 128), differenza 1,9·10⁻¹³.

**C3: e le regole imparate dall'agglomeratore?** La regola grossa imparata (1-D, livelli B = 16, 64, 192) non ha un cono
perfetto: è il fattore unitario di uno stencil, quindi ha una coda. Partendo da **una sola cella accesa** (il caso
peggiore: contenuto fino al limite del livello) dopo un passo l'1 % circa della norma sta fuori dal cono fisico
(1 cella grossa per passo), ma fuori da **2 volte lo stencil** l'ampiezza è ≤ 6·10⁻⁶ dopo un passo e 10⁻¹⁵ dopo 4.
Regola pratica per il motore: sui livelli imparati, congelare fuori da un cono di velocità **2w celle per passo** (w =
raggio dello stencil), con un errore ≤ 10⁻⁵ per passo, controllabile.

**Berger–Oliger (calcolo, non esame).** Nel mondo R3 del reticolo (21 livelli, 717 568 celle) ogni livello fa metà dei
tick di quello sotto: costo per tick fine 127 357 contro 717 568, **×5,6**. Il cono e l'orologio pigro sono il caso
estremo dello stesso principio: il passo di tempo di una regione che nessuno può vedere è **infinito**.

## 6. Cambiare livello senza scatti (D1, D2)

**Luce (D1).** Pacchetto a 30°, proiettato sulla banda giusta del reticolo fine, portato al livello grosso e
ritorno:

| celle/λ sul grosso | spettrale: ΔE/E | cubico (oggi nel reticolo): ΔE/E andata e ritorno | media a blocchi | banda sbagliata sul grosso |
|---|---|---|---|---|
| 32 | −2·10⁻¹² | −3,5·10⁻⁵ | −2,4·10⁻³ | 3,5·10⁻⁶ |
| 16 | −4·10⁻¹² | −5,5·10⁻⁴ | −9,6·10⁻³ | 5,7·10⁻⁵ |
| 12 | −6·10⁻¹² | −1,7·10⁻³ | −1,7·10⁻² | 1,9·10⁻⁴ |
| 8 | −9·10⁻¹² | −8,6·10⁻³ | −3,8·10⁻² | 1,0·10⁻³ |

* Lo **spettrale** (tieni i coefficienti di Fourier sotto il limite del livello grosso, Parseval) conserva l'energia a
  10⁻¹¹, il centro non salta (10⁻¹⁵ λ), il giro di andata e ritorno ha fedeltà 1. **È il cambio "senza pop".**
* La parte che finisce sulla **banda sbagliata** del reticolo grosso (elicità o modo fermo) non dipende dal metodo: è
  il fatto che i "fotoni del reticolo" a 12 celle/λ sono un po' diversi da quelli a 24. 1,9·10⁻⁴ a 12 celle/λ,
  sotto la soglia 10⁻³; a 8 celle/λ arriva a 10⁻³ (un altro segnale del pavimento).
* La media a blocchi (il modo più ovvio) **perde l'1–4 % dell'energia** a ogni cambio: un motore che cambia livello
  spesso così farebbe "svanire" la luce.

**Elettrone (D2).** Il descrittore dell'agglomeratore è la media sul blocco. Quella media **non** conserva
esattamente la norma (−3,7·10⁻⁴ a B = 192, k = 1/a₀: è un filtro). La soluzione esatta:
* la **carica vestita** del livello grosso, Q = Σ_K |ψ̂(K)|²/|S(K)|² (S = il filtro del blocco), è uguale alla carica
  fine a 10⁻¹⁶ ed è conservata dalla regola grossa (che è unitaria e uguale dappertutto);
* il **prolungamento spettrale** (dividere per S in Fourier: "disfare" il blocco) riporta stato, norma ed energia
  cinetica al fine a 10⁻¹⁵.

Il prolungamento a gradini (ripetere il valore della cella grossa) conserva la norma ma **inventa energia cinetica**:
+105 % a B = 16, +1850 % a B = 192, +3800 % a B = 384. È il "pop" fisico: gli spigoli dei gradini sono onde corte
finte.

## 7. La mela in 4K e in 480p: che cosa deve restare (EL, EE)

**La definizione (pre-registrata).**

| deve restare uguale | come si misura | tolleranza Ultra / Intel |
|---|---|---|
| leggi di conservazione | norma (carica), energia sum\|F\|² nel livello | **10⁻¹⁰** (tutti) |
| tempi caratteristici e frequenze | velocità, frequenza, apertura per diffrazione, energie cinetiche | 2 % / 10 % |
| sequenza degli eventi | due eventi separati più di tol·t restano nello stesso ordine (≡ errore di velocità < tol/2) | 1 % / 5 % sulle velocità |
| simmetrie | purezza di elicità, isotropia (max − min della velocità sugli angoli) | 2 % / 10 % |
| attraversare le cuciture fra livelli | T ≥ 1 − tol/20, energia persa e elicità sbagliata ≤ tol/20, velocità dopo entro tol | 0,999 / 0,995 |

**Può peggiorare:** la nitidezza (la cella), il contenuto più corto di ~12 celle, la Zitterbewegung.

**Luce** (riferimento: Maxwell continuo; la cucitura col mondo a blocchi di `reticolo_core`):

| | Ultra 32 celle/λ | Intel 12 | sotto il pavimento 8 |
|---|---|---|---|
| velocità peggiore (45°) | −0,48 % | −3,58 % | −8,47 % |
| frequenza peggiore (45°) | −0,16 % | −1,17 % | −2,70 % |
| diffrazione 0° / 45° (misura; previsione) | 0,987 / 1,008 (0,987 / 1,008) | **0,907** / 1,057 (0,907 / 1,058) | **0,786** / 1,131 |
| cucitura T peggiore (30°) / energia persa | 0,99993 / 7·10⁻⁵ | 0,9962 / 3,8·10⁻³ | **0,983** / 1,7·10⁻² |
| anisotropia | 0,48 % | 3,6 % | **8,5 %** |
| energia nel livello | 1,5·10⁻¹² | 5·10⁻¹³ | 3·10⁻¹³ |
| verdetto | **PASS** | **PASS** | **FAIL** (diffrazione, cucitura, sequenza) |

**Elettrone** (livelli 1-D imparati da `agg_zz_1d`, riferimento: la regola fine esatta; pacchetto atomico k₀ = 1/a₀,
σ = a₀/2, fino al tempo di allargamento 2mσ² ≈ 9400 ħ/mc²):

| | Ultra B = 16 (a₀/256) | Intel B = 192 (a₀/21) | sotto B = 384 (a₀/11) |
|---|---|---|---|
| velocità | +1,4·10⁻⁵ | +7,9·10⁻⁴ | +7,2·10⁻⁴ |
| allargamento | +7·10⁻⁶ | −2,1·10⁻⁴ | −2,2·10⁻⁴ |
| energia cinetica a k = 0,5 / 1 / 2 (1/a₀) | 3·10⁻⁵ / 2·10⁻⁵ / 1·10⁻⁵ | 3,6·10⁻³ / 1,7·10⁻³ / 7,8·10⁻⁴ | 3,4·10⁻³ / 1,5·10⁻³ / 6,8·10⁻⁴ |
| norma | 2·10⁻¹³ | 2·10⁻¹³ | 1·10⁻¹⁴ |
| verdetto | **PASS** | **PASS** | **PASS** |

* **La luce ha un pavimento a ~12 celle/λ**, fissato dalla diffrazione (ka/tan ka), non dalla velocità. Intel a 12 passa,
  ma di poco (−9,3 % contro 10 %). Sotto (8) la luce **cambia fisica**: si apre il 21 % in meno, perde l'1,7 % alle
  cuciture, e due fotoni che arrivano a 5 % di distanza possono scambiarsi l'ordine.
* **L'elettrone a scala atomica non tocca il pavimento** nemmeno a a₀/11: le regole imparate sono già Schrödinger con la
  massa giusta e l'errore è di qualche 10⁻³. (Il pavimento per la materia viene dal nucleo: la cuspide degli stati s,
  `AGGLOMERATO_ZIGZAG.md`: 2s +1,0 % ad a₀/8.)
* **Zitterbewegung (EE5, descrittivo).** Ad Ultra (passo 0,54 ħ/mc²) il tremolio si vede alla frequenza giusta (2,03
  contro 2). A Intel (passo 6,4) non è "sparito": la regola imparata con i due fogli (R, L) lo **ricampiona** e mostra un
  battito **lento finto** a 0,12 mc²/ħ, identico a quello del fine guardato agli stessi istanti. Per il motore è una
  trappola visiva: un'oscillazione che sembra fisica e non lo è. Ai livelli grossi va usato il descrittore di **Pauli**
  (solo energie positive), come già fa la catena 3-D: lì il tremolio non c'è e non c'è nemmeno il fantasma.

## 8. Il budget sull'hardware (calcolo, da misurare col ponte)

**Dati.** Misurato nel Laboratorio (dal coordinatore): ~42 ms per passo fisico con ~33 800 palline, sulla Intel UHD
(oggi la fisica gira sulla CPU; il porto WGSL della regola non c'è ancora). Dati di catalogo, **non misurati**: Intel
UHD della famiglia del portatile ≈ 0,4 TFLOPS, memoria condivisa DDR4, banda pratica ~25 GB/s; scheda da 16 GB ≈
288 GB/s (fascia media) – 717 GB/s (fascia alta).

**Conto.** La regola della luce è limitata dalla banda di memoria: 24 byte per cella (3 complessi a 32 bit), da 48 byte
per cella-tick (passata fusa) a 144 (tre passate per asse). Celle×tick al secondo:

| | Intel UHD | 16 GB fascia media | 16 GB fascia alta |
|---|---|---|---|
| celle×tick/s | 1,7·10⁸ – 5,2·10⁸ | 2,0·10⁹ – 6,0·10⁹ | 5,0·10⁹ – 1,5·10¹⁰ |
| livello | pavimento: 12 celle/λ | Ultra: 32 | Ultra: 32 |
| scatola 3-D con il fotone che avanza di 1 λ al secondo (costo (NL)³·N) | **L ≈ 20–29 λ** | L ≈ 12–18 λ | L ≈ 17–24 λ |
| memoria (due copie del campo) | (12·25)³ celle → 1,3 GB | (32·20)³ → 12,6 GB | idem |

**La cosa bella:** perché il costo va come N⁴ in 3-D, il rapporto Ultra/Intel (32/12)⁴ ≈ 50 è dello stesso ordine del
rapporto di banda (12–29). Quindi **la stessa scena, grande ~20 λ, gira sulla Intel a 12 celle/λ e sulla scheda da
16 GB a 32 celle/λ**: è esattamente "la mela in 480p e in 4K". Se la scena è più grande, sulla Intel si **taglia la
scatola** (o si congela fuori dal cono, §5), mai la risoluzione sotto il pavimento.

**Atomo** (regola imparata di Pauli, due FFT per passo, ~100 byte per cella): Intel a₀/8 (96³, ±6 a₀) ~3,5 ms per passo;
Ultra a₀/16 (192³) ~2,4 ms sulla scheda media. Ma un battito 1s–2p dura 18 000–37 000 passi: **circa un minuto** per
battito su entrambe. Per vedere battere un atomo isolato dal vivo serve il salto esatto (idea 4 sotto), non più GPU.

**Richieste di ponte proposte** (formato `docs/PONTE.md`, **non create**):
1. `{"azione": "ping"}` — conferma che wgpu vede la Intel UHD e quanta RAM condivisa c'è.
2. `{"azione": "test", "args": {"filtro": "motore_banco_regola", "ignored": true}}` — dopo la patch §10.3: tempo per tick
   della regola della luce sulla CPU a 96³ e 128³, per tarare i byte per cella-tick.
3. Dopo il porto WGSL della regola (thread "Modello standard"): lo stesso banco sulla GPU, con il nome dell'adattatore,
   per sostituire le stime della tabella con numeri veri.

## 9. Idee geometriche e visive

| # | idea | cosa si vedrebbe | perché funziona | costo | regole |
|---|---|---|---|---|---|
| 1 | **Il cono di luce come frustum** (provata, C1) | niente di diverso a schermo: è il punto. In una vista di controllo, le palline fuori dal cono restano ferme e si rimettono in moto quando entrano | la pallina zig-zag non supera mai una cella per tick: dipendenza esatta | zero; risparmio ×2,7–12 in 2-D, al cubo in 3-D | fisica pura, nessun fake |
| 2 | **Clipmap di livelli attorno all'occhio** | gusci concentrici, ognuno grande il doppio e con celle il doppio, centrati sulla camera; mentre ti muovi i gusci scorrono | è il "camino" del reticolo (R4: 25 610 celle per livello), ma attorno all'occhio invece che attorno a un oggetto: costo ∝ numero di livelli = log(mondo) | 20 livelli × 64³ ≈ 5·10⁶ celle per un mondo 2²⁰ volte la cella fine | si scala solo con l'agglomeratore; il contenuto vince sull'occhio (il pavimento §7) |
| 3 | **Passaporto d'errore come colore** | in vista di controllo, ogni regione ha un alone del colore del suo errore accumulato in pixel (verde < 0,5, rosso > 1) | il passaporto è già calcolato (§4) | nullo | è solo un indicatore; niente linee |
| 4 | **Atomo-cartolina spettrale** | un atomo lontano che "respira": le palline seguono \|Σ cₙ ψₙ e^{−iEₙt}\|², battono alle frequenze vere | un atomo isolato è lineare: nella base dei suoi stati fermi (calcolati dalla regola stessa, come in A2 di AGGLOMERATO_ZIGZAG) il tempo è un giro di fasi esatto; il salto di §4 per la materia | pochi numeri per atomo; errore = solo il peso degli stati tagliati, misurabile | gli stati vengono dalla regola, nessuna tabella; si torna alla griglia quando arriva la luce |
| 5 | **Reticolo più rotondo (BCC/FCC)** | la luce di sbieco non rallenta: fronti d'onda davvero sferici | il difetto obliquo (−3,6 % a 12 celle/λ) è il cubo; esistono automi di Weyl sul reticolo cubico a corpo centrato (D'Ariano–Perinotti) più isotropi | da calcolare: dispersione e anisotropia a 12 celle/λ | **cambia la regola**: non si fa senza Valerio (una sola strada) |
| 6 | **Il pavimento come "risoluzione onesta"** | sulla Intel la stessa scena, più piccola o più sgranata, ma con la luce che si apre e arriva come sulla scheda grossa | §7–§8: sotto ~12 celle/λ cambia la fisica, sopra cambia solo la nitidezza | — | è la regola del motore |

La più promettente nel mio perimetro (1) l'ho provata: esatta al bit.

## 10. Patch proposte per altri file (non toccati)

1. **`src/atomo/` (thread Laboratorio Atomo)** — modulo nuovo `livelli.rs`:
   * `TabellaErrori`: per N ∈ {8, 12, 16, 24, 32, 48, 64} e angoli a passi di 5°, δω/ω, errore di velocità, walk-off,
     fattore di diffrazione, calcolati all'avvio dalla matrice di Bloch (porto di `motore_core.level_errors`; con il
     `zigzag/walk.rs` esistente).
   * `Passaporto { livello, deriva_acc: [f32; 3], fase_acc: f32, t_chiave, stato_chiave }` per regione; scelta del
     livello con l'impronta minima raggiungibile dalla camera (velocità massima della camera).
   * `cono.rs`: maschera L∞ delle celle da aggiornare (distanza dalla finestra ≤ tick rimanenti); sui livelli
     imparati velocità 2w.
   * `trasferisci.rs`: cambio di livello spettrale (rustfft, già nell'ecosistema) per luce ed elettrone, con
     deconvoluzione del blocco.
2. **`lab/reticolo_core.py` / `reticolo_spec.json`** — sostituire, per il cambio di livello di un'intera regione, la
   restrizione cubica (−1, 9, 9, −1)/16 con quella spettrale (debito 14: chiude anche il "fantasma del gamma" per il
   contenuto risolto). La cucitura locale fra blocchi resta cubica; quella conservativa del flusso è ancora da fare.
3. **`src/atomo/tests.rs`** — test ignorato `motore_banco_regola`: 200 tick della regola della luce a 96³ e 128³,
   stampa ms/tick e byte/s stimati (per la richiesta di ponte 2).
4. **`docs/PIANO.md`, `docs/CONTESTO_CLOUD.md`** — aggiungere questo studio; in CONTESTO §4 debito 25 ("≥ 24 celle/λ a
   ≥ 30 tick/s non entra nel portatile"): sostituire con il pavimento 12 celle/λ e la tabella di §8.
5. **Agglomeratore 1-D** — ai livelli grossi (passo > π/2 ħ/mc²) usare il descrittore a energia positiva, per evitare
   il battito finto della Zitterbewegung ricampionata (EE5).

## 11. Fallimenti, debiti, cambi rispetto alla pre-registrazione

**FAIL:** B1 (sotto-criterio dei pacchetti segnalati: difetto della mia misura, emendamento A1 lo mostra); EL "sotto il
pavimento" (atteso, è la prova che il pavimento esiste).

**Cambi e correzioni dichiarati:**
* Prima della pre-registrazione: esplorazioni di sviluppo elencate nel JSON.
* Nel pre-registro, prima delle misure: il costo di B1 contato sui soli pacchetti non segnalati (la prima stesura
  contava i segnalati come N = 64).
* B2: la scala del passaporto estesa a 96 e 128 celle/λ (nessun livello fino a 64 tiene 0,5 px a 0,1 λ/px per
  100 λ: lo dicono già le segnalazioni di B1). B2 riorganizzato per far girare di fila gli intervalli con lo stesso
  livello (stessa fisica, solo più veloce: la prima versione non stava nei 10 minuti).
* D2: un errore di arrotondamento nel mio script (indici di Fourier `astype(int)` invece di `rint`) dava −11 % a
  B = 16; corretto e rifatto, l'esame non è cambiato.
* Emendamento A1 (dopo B1): i due pacchetti segnalati rifatti in una scatola più grande; il verdetto non cambia.

**Debiti nuovi o rimasti:**
1. Le regioni di B1 sono pacchetti in scatole separate (la luce non interagisce con la luce); la convivenza di più
   livelli in un mondo è quella del reticolo (R3b, 4,8·10⁻⁷), non rifatta qui.
2. Il cambio di livello spettrale è globale sulla regione (FFT): per i blocchi del reticolo serve una versione locale.
3. Il passaporto è provato solo per la luce libera; per la materia (atomo in un campo) l'errore per livello va preso dai
   controlli dell'agglomeratore (8,8·10⁻⁵ … 2,8·10⁻⁴ di energia per livello in 3-D), non verificato in una scena.
4. Solo 2-D (luce) e 1-D (elettrone). Il 3-D è un conto (§8): l'errore della luce è lo stesso (la dispersione), il
   costo va al cubo/quarta.
5. Tutti i tempi GPU sono stime di catalogo: da misurare col ponte.
6. L'atomo-cartolina (idea 4) e la clipmap (idea 2) non sono provate.

**Quanti numeri:** **0 parametri di fisica** (c e la regola sono quelli di `ZIGZAG.md`; i livelli dell'elettrone sono quelli
già imparati in `agg_zz_1d`). **Manopole di progetto** (scelte, non fisica): tolleranze 2 % / 10 %, margine 0,5 px,
scala dei livelli, intervallo di controllo 10 λ, pavimenti 32 / 12, regola del cono 2w sui livelli imparati.

## 12. Dove l'intuizione di Valerio può aiutare

1. **Quanto deve "marcire uguale" la mela?** Ho messo 10 % per la Intel. Per te conta di più che i **tempi** siano
   giusti (il fotone arriva quando deve) o che le **forme** lo siano (quanto si apre il fascio)? Sulla Intel a 12 celle/λ
   i tempi sbagliano del 3,6 %, le forme del 9 %.
2. **L'errore che ha memoria.** Accettiamo che una cosa vista da lontano per tanto tempo, quando ti avvicini, "si
   ricalcoli" (pausa di un istante, rifatta dal fotogramma chiave)? Oppure preferisci che il motore la simuli sempre
   abbastanza bene "per sicurezza"? La prima costa 500 volte meno nel caso di B2.
3. **Il cono di luce come regola del mondo.** Fisicamente dice: quello che non può raggiungerti non conta *ancora*. Ti
   piace che il motore **fermi l'orologio** delle zone che nessuno può vedere? Si vede solo se lo si mostra apposta.
4. **Reticolo più rotondo.** La luce di sbieco rallenta perché le palline vivono su un cubo. Un reticolo a "corpo
   centrato" (come il ferro) sarebbe più rotondo, ma è un'altra regola. Vuoi che lo valuti (solo calcolo, senza cambiare
   strada)?
5. **L'atomo lontano come cartolina che respira.** Un atomo da solo è "una musica": pochi stati che girano. Ti sembra
   onesto disegnarlo così da lontano (stati calcolati dalla regola), e tornare alle palline vere quando arriva la luce?

## 13. Cosa resta da fare

1. Ponte: `ping` e, dopo la patch, il banco della regola (richieste di §8).
2. Passaporto per la materia: un atomo agglomerato colpito da un impulso, Ultra (a₀/16) contro Intel (a₀/8), con le
   stesse definizioni di "dinamica conservata" (righe, tempi, conservazione).
3. Atomo-cartolina spettrale: provare che il salto nella base degli stati fermi riproduce la griglia, e quanto costa
   tornare alla griglia quando arriva la luce.
4. Cambio di livello spettrale **locale** per i blocchi del reticolo.

## Appendice — idee di disegno da passare al thread "Modello standard a palline"

* I descrittori della super-pallina di `ZOOM.md` (n, X, s²: quante, dove, quanto larga) **sono già uno splat gaussiano**:
  ogni agglomerato si disegna come una pallina morbida con centro, larghezza e peso suoi.
* Lontano, una pallina per impronta di pixel che porta la **somma dell'energia** delle celle che copre (un mipmap
  dell'energia): la luminosità totale si conserva esattamente quando si cambia livello, quindi niente "pop" di luce.
* Le frange più strette di 2 pixel non si disegnano (il pixel ne fa la media, come in §3): si disegna l'inviluppo.
* "Lunga esposizione" = quello che misura un rivelatore vero (energia integrata su un tempo): è fisica, e le palline
  che la campionano possono continuare a muoversi lungo il flusso d'energia.
* Vista di controllo: palline colorate per livello e per errore del passaporto (idee 1 e 3), mai linee.
* Budget di §8: la Intel regge ~20 λ di scatola a 12 celle/λ; una 16 GB la stessa scatola a 32.

## Come rifarlo

    python lab/motore_prereg.py
    python lab/motore_m1.py A ; python lab/motore_m1.py B1 ; python lab/motore_m1.py B2   # ~1, 2, 5 min
    python lab/motore_m2.py                                                               # ~1 min
    python lab/motore_m3.py                                                               # ~3 min
    python lab/motore_m4.py EE ; python lab/motore_m4.py EL4 ; python lab/motore_m4.py EL3 ; python lab/motore_m4.py judge
    python lab/motore_amend.py ; python lab/motore_figure.py
