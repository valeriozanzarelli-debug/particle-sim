# Bagliori: il fotone che si legge anche da davanti

Studio del 13/10/2026. Prefisso `gl_`. Pre-registrazione `lab/out/gl_prereg.json` (commit 051f1c7, prima di
ogni misura), emendamenti `lab/out/gl_amend.json`, numeri `lab/out/gl_risultati.json`.
Clip: `lab/out/gl_fotone_davanti.mp4/.gif` (di fronte), `lab/out/gl_fotone_lato.mp4/.gif` (di lato).

## 0. In breve

Regola di Valerio: **solo palline e bagliori**, niente linee, scie, palline allungate o surfisti speciali; la
dinamica la fanno i bagliori che si muovono. Il fotone visto **di fronte** era un disco bianco uniforme. Adesso:

* ogni pallina di energia è un **bagliore tondo** la cui **grandezza e luce crescono con la densità di energia**
  u = E² + c²B² nel punto (dove la luce è forte i bagliori sono più grandi e più accesi);
* il **colore resta la direzione di E**, ma la tavolozza (scritta in sRGB) ora è convertita in colore lineare
  prima del disegno HDR: prima ogni colore usciva pastello, e migliaia di pastelli di tutti i colori fanno bianco;
* le palline si disegnano **dalla più lontana alla più vicina** (CPU: ordinamento; GPU: 4096 "cassetti" di
  profondità, un ordinamento per conteggio che regge anche il milione di palline dell'ultra);
* `--sequenza K [--passo P]` (nuova opzione del Laboratorio): salva un PNG per tick, esattamente al tick voluto,
  per fare clip al rallentatore anche senza schermo (`lab/gl_clip.sh`).

Tutti gli esami del fotone **passano** (con un emendamento dichiarato su I3, sotto).

## 1. Cosa si vede

**Di lato** (`gl_fotone_lato`): bande di colore che scorrono a c: è E che gira attorno all'asse (luce circolare:
tutti i colori una volta per λ). **Di fronte** (`gl_fotone_davanti`): un nucleo più pieno e luminoso e un bordo
rado; le macchioline colorate del bordo **girano attorno all'asse** mentre il fotone arriva: è il vortice
dell'energia, cioè lo spin del fotone. Il giro è lento (~80° in tutto il passaggio): nella clip si vede guardando
il bordo.

## 2. I numeri

| esame | cosa | misurato | soglia | esito |
|---|---|---|---|---|
| S1 | verso del vortice al bordo (0,7–1,3 w0), dalle palline disegnate | +0,02414 (circ+) / −0,02414 (circ−) rad/tick | segni opposti | PASS |
| S2 | luce lineare | +0,00001 rad/tick (0,04 % della circolare) | < 15 % | PASS |
| S3 | contro il fascio gaussiano: Ω = 2σc/(k w0²) = 0,02653 | rapporto 0,910 | 0,6–1,4 | PASS |
| I1 | pixel saturi nella foto di fronte | 0 % | ≤ 5 % | PASS |
| I2 | luce del nucleo / luce della corona (1,2–1,6 R50) | 2,91 | ≥ 2 | PASS |
| I3 | saturazione del colore (non nebbia bianca) | 0,359 (prima dell'emendamento 0,213) | ≥ 0,35 | PASS dopo em. 2 |
| I5 | stessa foto sul cammino GPU | I1 0 %, I3 0,359, I2 uguale alla CPU | | PASS |

**Validazione, non scoperta.** S3 ritrova la formula nota del flusso di spin di un fascio (v_φ = −(σc/2k) ∂_r ln u);
il 9 % in meno viene da un fascio focalizzato e corto (σ = λ), non gaussiano infinito. Il vortice c'era già
(lo spin ±0,999 ħ di FOTONE_PALLINE F2); qui si vede e si misura sulle palline che vanno sullo schermo.

## 3. Fallimenti ed emendamenti

* **I3 fallito alla prima foto (0,213).** Causa: la tavolozza in sRGB passata come colore lineare (e non la
  profondità, come pensavo: l'ordinamento da solo non cambiava niente). Emendamento 2: conversione ~x^2,2. Soglia
  non toccata. Le manopole del bagliore sono rimaste quelle scritte prima delle foto.
* Con 20 000 palline (preset bassa) di fronte il colore NON può essere uniforme: per ogni pixel la pallina davanti
  sta a una profondità a caso entro ~λ, e il colore cambia ogni λ/5. Si vedono coriandoli colorati, non un colore
  solo. Con molte più palline (alta/ultra, GPU seria) lo strato davanti diventa coerente.

## 4. Manopole (scelte di disegno, non fisica): 7

`fotdisegno::BAGLIORE`: grandezza 0,45 + 0,85 √u, luce min(0,03 + 0,42 u, 0,5), fattore colore 1; più
l'esponente della conversione sRGB (2,2) e i 4096 cassetti. `ATOMO_BAGLIORE="s0,s1,a0,a1,p,amax,c"` le cambia per
provare senza ricompilare.

## 5. Debiti

* Lo stretch lungo E e la scia di 4 campioni sono ancora nel codice (spenti), per gli esami vecchi della GPU.
* L'ordinamento è rifatto quando la camera si sposta più dell'1 % della distanza: girando molto in fretta a fotone
  fermo l'ordine può restare indietro di un fotogramma.

## 6. Prossimo passo: l'elettrone libero (e il muone) a bagliori

Fermato per priorità (gluoni e quark). Già fatto, NON finito: `src/atomo/elscena.rs` (scena `--scena elettrone`,
`--scena muone`), pacchetto di F3 di ZIGZAG.md, 40 000 bagliori sulla corrente. Prima passata degli esami
(`src/atomo/test_bagliori.rs`): il tremito (Zitterbewegung) si ritrova dai bagliori disegnati, **E1 PASS**
(ω 2,081 / 2,079 contro 2,071) ed **E2 PASS** (raggio 0,424 / 0,427 contro 0,433 / 0,438); **E3 FAIL**: il
centro dei bagliori anticipa quello dell'onda (scarto max 0,285 ħ/mc, soglia 0,06). Con la corrente esatta per
link (emendamento 1, `elscena::flussi`, controllata a 1e-9 da `bagliori_flussi_esatti`) lo sfasamento sparisce ma
resta 0,12 (deriva dell'interpolazione lineare dentro la cella: serve il tracciamento di Pollock di
SURFISTI_ESATTI). Il lettore rifà le proiezioni della regola: debito, la patch `docs/patch/walk_tick_flusso.patch`
(non ancora compilata) lo toglierebbe. Mancano: Pollock, clip, foto, il muone (massa × 206,768 = ingresso: senza
nucleo è lo stesso disegno 207 volte più piccolo e veloce, nessun numero nuovo).
