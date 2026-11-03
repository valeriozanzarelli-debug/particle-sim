# Tavola periodica, parte 2: dal boro in su sulla griglia agglomerata — PARZIALE (ossigeno: HF sì, correlazione no; fase 3 di coppia: FAIL)

Studio `tavg_`. Stato al 2/11/2026: **PARZIALE**. Pre-registrazione della fase 1: `lab/out/tavg_prereg.json` (commit
69414b4). Emendamento fase 2 (correlazione, bersaglio ossigeno): `lab/out/tavg_prereg2.json` (commit ec4a4ba, prima
di ogni fit). Verdetti di O (con N, C, H): `lab/out/tavg_verdetti_ossigeno.json` (`python lab/tavg_verdetti_o.py`).

## Ossigeno, in breve (2/11)

- **Hartree–Fock dell'ossigeno funziona sui tre livelli.** Fine: O −74,8093985 (chiave −74,809398, 5·10⁻⁷), O⁺
  −74,3726057 (3·10⁻⁷). Griglia 3-D a a₀/4 col blocco imparato (r_c 1,2/1,175 a₀): O −74,812490, O⁺ −74,372471.
- **La rottura di simmetria di −3,1 mEh è fisica di HF, non la griglia** (prova `lab/tavg_simm.py`): con i tre p
  costretti ad avere la stessa forma radiale (quella del livello fine) la griglia ridà il radiale entro +0,11 mEh; se li
  lascio liberi scende di −3,09 mEh a a₀/4 e −3,05 a a₀/6, uguale. Il p doppiamente occupato si allarga (⟨r⟩ 1,276
  contro 1,204 a₀). È HF ristretto con orbitali non equivalenti (noto in letteratura), non UHF (che dà −10,4 mEh).
- Costo: il −3,1 mEh sta solo nel neutro, così l'IP HF della griglia (11,974 eV) supera quello del livello fine
  (11,886) di 0,088 eV: **G4 FAIL** (soglia 0,05). È la stessa fisica, ma la soglia è quella e resta.
- **La correlazione imparata NON regge sull'ossigeno.** Imparata solo dagli esatti della regola (H⁻, He, He 2³S, Li⁺,
  Li, Li 2P, Be⁺, ioni tipo He Z = 6–10; 12 punti), forma LYP a 4 manopole (scelta da sola per validazione incrociata:
  1,6 mEh contro 2,4 della forma a 2), prevista alla cieca: E_c(O) = −172 mEh contro −258 (**−33 %, C1 FAIL**);
  IP(O) = 13,204 eV contro 13,618 (**−3,04 %, C2 FAIL** per un soffio; con l'HF del livello fine 13,116, −3,7 %).
  Va invece bene la riga: IP(N) 14,139 contro 14,534 (−2,7 %, C5 PASS), caduta N→O 0,935 eV contro 0,916 (C3, C4 PASS).
- **Perché**: nell'allenamento ci sono solo coppie del guscio interno (1s², 1s2s, 1s2p) e lì E_c è quasi costante
  (42–47 mEh). Nell'ossigeno ~¾ della correlazione sta nelle coppie del guscio 2 (2s², 2s2p, 2p²), che nessun dato
  ha mai mostrato. Prova (post-hoc, non contata): le costanti LYP pubblicate (prese da He soltanto) danno −265 mEh
  su O ma sbagliano i NOSTRI esatti (Be⁺ −61 invece di −47): la forma a densità non distingue i tipi di coppia, il
  risultato giusto su O della letteratura è in parte fortuna. Togliere gli ioni He-like o cambiare i pesi non cambia
  niente (−171/−174 mEh).
- **Export di O: NON fatto** (G4 e C1, C2 falliscono). Il blocco O (`lab/out/tavg_blocchi/Z08.json`, HF) è usabile
  con `hf_griglia` con errore dichiarato (−3,1 mEh sul totale, +0,09 eV sull'IP HF).
Tier di questi risultati: **"approssimato"** (Hartree–Fock + blocchi imparati; non esatto come H, He, H₂, Li, Be).

## Cosa si fa, detto semplice

L'ECG esatto non va oltre 4–5 elettroni. Per B…Kr la strada che scala è a **tre livelli dell'agglomeratore**:

1. **Livello fine = zoom radiale** (`lab/tavg_radial.py`): l'atomo intero, tutti gli elettroni, con l'operatore della
   regola (p²/2m_kin − Z/r + 1/r_ij), Hartree–Fock ristretto del termine di Hund. Nessun orbitale da tabella: si parte
   dagli orbitali del nucleo nudo e il campo medio si costruisce da solo.
2. **Blocco-nucleo** (`lab/tavg_core.py`, `lab/tavg_blocks.py`): nucleo + gusci interni diventano UNA pallina. La
   sua regola (come spinge gli elettroni esterni, un potenziale per ogni l) si impara invertendo l'equazione del livello
   fine; la grandezza del blocco r_c la sceglie la macchina (il più piccolo blocco che la griglia grossa riesce a
   portare), e la cella h la sceglie la macchina (la più grossa fra a₀/4, a₀/6, a₀/8 che trasferisce il catione entro
   0,8 mEh). In linguaggio standard è uno pseudopotenziale a conservazione di norma generato da noi: lo dichiaro.
3. **Griglia 3-D** (`lab/tavg_grid.py`): gli elettroni esterni su celle cubiche, con i blocchi in posizioni qualsiasi.

## Numeri che ci sono (produzione, contro la pre-registrazione)

| esame | cosa | risultato | esito |
|---|---|---|---|
| F1 (B, C) | ROHF fine contro limite HF | B −24,5290607 (chiave −24,529061, +3·10⁻⁷); C −37,6886190 (chiave −37,688619, ~0) | **PASS** (2 su 6) |
| F2 (B⁺, C⁺) | catione | B⁺ −24,2375752 (−24,237575); C⁺ −37,2922238 (−37,292224) | **PASS** (2 su 6) |
| F4 (B, C) | aufbau: la configurazione più bassa fra i candidati è quella vera | B 2s²2p, B⁺ 2s², C 2s²2p², C⁺ 2s²2p | **PASS** (4 casi su 40) |
| IP HF (fine) | ΔSCF | B 7,932 eV, C 10,786 eV (sperimentali 8,298 e 11,260: manca la correlazione, fase 2) | solo riportato |
| REL | stima relativistica su IP | B −0,0023 eV | trascurabile |
| B2 (B) | blocco → catione | +0,37 mEh a a₀/4 | riportato |
| B3 (B) | promozione 2s→2p | −1,51 mEh (soglia 1,6) | **PASS** |
| F1, F2 (N, O) | fine contro chiave | N −2·10⁻⁷, N⁺ −1·10⁻⁸, O −5·10⁻⁷, O⁺ +3·10⁻⁷ | **PASS** |
| F3, F4 (N, O) | gusci, Kato, aufbau | 2 massimi; Kato 6,9993 / 7,9992; aufbau giusto per N, N⁺, O, O⁺ | **PASS** |
| REL (N, O) | stima relativistica su IP | N −0,007 eV, O −0,010 eV | sotto la bandiera |
| B3 (N, O) | promozione 2s→2p | N −3,05, O −3,40 mEh | **FAIL** (atteso) |
| G1 (N, O) | griglia contro radiale, stesso blocco | N −0,02, N⁺ −0,08; O −3,09 (finestra non sferica), O⁺ +0,07 mEh | **PASS** |
| G2 (O) | a₀/4 contro a₀/6, L = 12 | 0,06 mEh | **PASS** |
| G3 (N, O) | totale griglia contro limite HF | N −0,02, N⁺ +0,24; O −3,09, O⁺ +0,14 mEh | **PASS** |
| G4 (N, O) | IP HF griglia contro fine | N +0,007 eV; O +0,088 eV | N **PASS**, O **FAIL** |
| G5 (O) | stencil imparato a₀/4 | −5,5 mEh | **FAIL** (atteso) |
| G6 (C, O) | uova: atomo spostato di (0,37, 0,21, 0,11) celle | C 0,098 mEh; O 0,165 mEh (a a₀/6: 6·10⁻⁷) | C **PASS**, O **FAIL** |
| G7 | molecole: H₂, N₂ | H₂ −0,014 mEh; N₂ +3,62 mEh (soglia 3) | H₂ **PASS**, N₂ **FAIL** |
| C1–C5 | correlazione imparata, cieca | vedi sopra | C1, C2 **FAIL**; C3, C4, C5 **PASS** |

Sviluppo (dichiarato, prima della pre-registrazione): livello fine Ne +4·10⁻¹⁰ Eh, Kr +9·10⁻⁹, Fe +2·10⁻⁶ dai limiti HF;
griglia 3-D contro il radiale con lo stesso blocco: B +0,001 mEh, C −0,11, Ne +0,08, O −3,0 mEh (a a₀/4 e a₀/6: è
rottura di simmetria vera, il p doppiamente occupato prende un raggio diverso); una corsa di griglia dura 5–20 s.
Ne ha bisogno della cella a₀/8 (catione: 24 mEh a a₀/4, 4,1 a a₀/6, 0,74 a a₀/8). C: la promozione 2s→2p sbaglia di
~4,5 mEh qualunque r_c (integrali valenza–valenza dentro il blocco): FAIL atteso per B3 su C, N, O.
Lo stencil imparato a₀/4 di aggmol vale solo fino a k ≈ 6/a₀ (−2,4 % a 7,9, −8 % a 9,4): per C…Ne la griglia usa la
banda della regola k²/2m spettrale (debito: ri-imparare il livello con l'osservabile dei gusci 2p).

## Funzione per gli altri studi (molecole)

`lab/tavg_grid.py`: `hf_griglia(nuclei=[(Z, (x, y, z)), ...], carica=0, spin=None, h=None, L=16.0, restricted=True)`
→ dict con `E` (totale, Eh), `E_val`, `orbitali`, `occ`, `eps`, `densita`, `griglia`. h = None prende la cella più
fine scelta dai blocchi presenti. Richiede i blocchi in `lab/out/tavg_blocchi/Z??.json` (oggi solo **Z05**, boro;
oggi Z01 H, Z05 B, Z06 C, Z07 N, Z08 O; gli altri con `./tavg_drive.sh tavg_blocks_prod.py <Z...>`).
**Bug trovato e corretto il 2/11**: in `_build_ext` il nome `Y` (griglia) veniva riusato per le armoniche sferiche:
con due nuclei e proiettori KB (N₂, e quindi H₂O) la funzione si rompeva. Un atomo solo e H₂ non erano toccati.
Per l'acqua: la O–H è a 1,81 a₀ e il blocco O ha coda fino a 1,2 + 2,75 a₀; N₂ (2,07 a₀) sbaglia di +3,6 mEh: aspettarsi
un errore di qualche mEh da sovrapposizione dei blocchi finché non si misura su H₂O stessa.

## Debiti

Forma TM del blocco e coda a carica puntiforme oltre r_c + 2,75 a₀ (manopole: bilancio di coda 2·10⁻⁴ Eh, 0,8 π/h,
0,8 mEh per scegliere h, canale locale p); nucleo congelato; Born–Oppenheimer; determinanti a orbitali reali (per d², d³,
d⁷, d⁸ non sono termini puri); cinetico della regola spettrale invece dello stencil imparato oltre k ≈ 6; relatività
solo stimata. Parametri liberi di fisica in HF: zero (ingressi Z e costanti). **Correlazione**: forma LYP presa dalla
letteratura (debito), 4 manopole fittate sui nostri 12 esatti (a = 0,0248, ab = 6,9·10⁻⁴, c = 0,30, d = 0,05 — d al
bordo della griglia pre-registrata: i dati non lo fissano). Esclusi dall'allenamento: Be (studio non chiuso), H₂ e LiH
(manca la loro densità HF sullo stesso livello).

## Come stringere la regola (proposta, invece di abbassare le soglie)

Il problema non è la forma: è che l'agglomeratore non ha mai visto una coppia del guscio 2. Regola proposta,
"correlazione per coppia, imparata nel campo vero":
1. **Coppie di valenza esatte sul blocco**: due elettroni nel campo del blocco O già imparato (1s² congelato) —
   cioè O⁶⁺ con la valenza 2s², 2s2p (³P, ¹P), 2p² (³P, ¹D, ¹S) — con l'ECG della regola (2 elettroni: costo piccolo).
   Pezzo nuovo da scrivere: il potenziale separabile del blocco dentro l'ECG. Uguale per C⁴⁺ e N⁵⁺ per vedere come la
   coppia cambia con la densità.
2. **Be quando chiude** (2s² con la quasi-degenerazione 2s/2p, ~45 mEh): primo esatto di una coppia del guscio 2 a
   densità "atomo neutro".
3. **Legge per tipo di coppia**: E_c = Σ_coppie ε_tipo(densità della coppia), con le manopole contate per tipo
   (s², sp, p² singoletto/tripletto, interguscio), poi un controllo a 3 elettroni (O⁵⁺ 2s²2p) per vedere se le coppie
   si sommano. Solo dopo, di nuovo O alla cieca con le stesse soglie (C1 5 %, C2 3 %).
Diagnosi gratis per la regola di h: anche le "uova" (G6) devono stare nel bilancio; O a a₀/6 passa (6·10⁻⁷).

## Come riprendere

```
cd lab
./tavg_drive.sh tavg_fine.py 5 36                       # livello fine (riprende dal checkpoint, B e C fatti)
./tavg_drive.sh tavg_blocks_prod.py 6 7 8 9 10 1 3 4 11 ... 36
./tavg_drive.sh tavg_run_grid.py atom 5 6 7 8 9 10       # poi conv, learned, shift, mol H2 N2
./tavg_drive.sh tavg_simm.py                             # prova della rottura di simmetria di O
python tavg_corr.py train && python tavg_corr.py predict # fase 2 (dopo tavg_prereg2.json)
python tavg_verdetti_o.py && python tavg_figure.py       # verdetti e figura lab/out/tavg_figure.png
```
Poi: verdetti G1–G7, export in `/mnt/project-files/debiti/dati/<nome>/` (tier "approssimato"), figura, fase 2
(correlazione imparata dagli esatti He-like/Li-like/Be-like della regola, `lab/tavg_ecg.py` copia di multi_e_ecg).

## Dove l'intuizione di Valerio può aiutare / Idee geometriche e visive

- **Il blocco-nucleo come pallina visibile**: a schermo il nucleo con i gusci interni è una sfera morbida di raggio r_c
  (0,8–1,5 a₀) che "sa" respingere gli elettroni esterni; zoomando dentro, la sfera si apre nei gusci veri del livello
  fine. Costo zero (i dati ci sono), rispetta le regole (è un livello dell'agglomeratore). Domanda: ti va che lo zoom
  sull'atomo mostri proprio questo passaggio blocco → gusci?
- **Le coppie come palline doppie**: la proposta sopra è vicina all'idea dei livelli: due elettroni che si evitano
  sono un "mini-agglomerato" con la sua regola imparata. Domanda per te: ti torna che ogni coppia di elettroni sia un
  oggetto con la sua correzione, come una molla che si impara una volta e si riusa nell'acqua?
- **La rottura di simmetria dell'ossigeno**: il p doppiamente occupato è più largo (⟨r⟩ 1,275 contro 1,204 a₀): a
  schermo una nube leggermente schiacciata lungo un asse. È fisica vera di HF, non un difetto.

Figura: `lab/out/tavg_figure.png` (a: correlazione imparata, allenamento e cieco; b: IP di N e O; c: la prova della
rottura di simmetria).

## Fase 3 — correlazione per coppia nel campo vero (CHIUSA il 3/11: **FAIL**, niente export di O)

Pre-registrazione: `lab/out/tavg_prereg_f3.json` (commit 94ecce4, prima di ogni corsa e di ogni fit), emendamento 1
`lab/out/tavg_prereg_f3_em1.json` (commit 2c90f3c: corse a 3 e- da K 250 a 600 perché non erano convergiute; fatto
prima di calcolare qualunque previsione; soglie intatte: E_c entro 5 %, IP entro 3 %). Verdetti:
`lab/out/tavg_f3_verdetti.json` (`python lab/tavg_f3_judge.py`), figura `lab/out/tavg_f3_figura.png`.

**Il collasso dell'SVM è risolto** (`lab/tavg_pecg2.py`). La causa: il proiettore del blocco scritto come 40 gaussiane
con coefficienti fino a 1,3·10⁵ che si cancellano (fattore 3,5·10⁵), più la sottrazione A_jj − A_ij²/A_ii. Rimedio:
formula riscritta senza quella sottrazione e calcolata in long double (80 bit): errore contro un riferimento a 45 cifre
2·10⁻¹¹ Eh invece di 8·10⁻⁸; in più l'SVM accetta un candidato solo se è nuovo per almeno 10⁻⁷ (norma residua) e usa
l'ortogonalizzazione canonica. Prima collassava a K ≈ 50 (energie di −10⁹ Eh); ora arriva a K = 240 senza mai
scendere sotto l'esatto. (Provata e scartata la quadratura diretta: errori 10⁻⁷, la griglia non risolve le funzioni
molto correlate.)

Corse (4 processi in parallelo, spezzoni da 590 s, `lab/tavg_f3_drive.py` + `lab/tavg_f3_loop.sh`, ~3 ore):
- coppie a 2 elettroni nei blocchi B, C, N, O con una CI nuova a onde parziali (`lab/tavg_ci2.py`, l ≤ 8, 24 funzioni
  di Laguerre per l, coda in l). Provata su He nudo: HF −2,8616800 (esatto), 1s² a l ≤ 4 entro 0,05 mEh dal limite.
  E_c delle coppie (mEh, B⁺→O⁴⁺): 2s² −65 → −115; 2s2p ³P −7,7 → −8,7; ¹P −51 → −63; 2p² ³P ≈ −9,6; ¹D ≈ −53;
  barre ≤ 0,09 mEh;
- ECG esatti nel blocco O per controllare la CI; ioni tipo Li nudi Z = 7, 8 (core-valenza); 2s²2p a 3 elettroni nei
  blocchi C e O (verifica di additività).

| esame | cosa | risultato | esito |
|---|---|---|---|
| M1 | macchina: oltre K = 200, senza collasso, monotona | He nudo dE 1,4·10⁻⁷ ok; 2 e- senza repulsione nel blocco O dE +4,8·10⁻⁶, mai sotto l'esatto, ma risale di 1,7·10⁻⁶ ai riavvii (taglio canonico a 10⁻¹²) | **FAIL** (solo "monotona") |
| V1 | ECG contro CI nel blocco O, entro 1 mEh | 2s² −0,08, 2s2p ³P −0,05, ¹P −0,10 mEh | PASS |
| A1 | additività a 3 e- (2s²2p), entro 10 % | C: coppie −116,3 contro esatto −88,7 mEh (+31 %); O: −154,1 contro −114,9 (+34 %) | **FAIL** |
| C1 | E_c(O) entro 5 % di −257,9 mEh | −492,2 mEh (+91 %) | **FAIL** |
| C2 | IP(O) entro 3 % di 13,618 eV | 16,17 eV (+18,7 %) | **FAIL** |
| C3, C4 | IP(N) > IP(O), caduta 0,916 ± 0,3 eV | IP(N) 15,59 < IP(O): caduta −0,58 eV | **FAIL** |
| C5 | IP(N) entro 3 % di 14,534 | 15,59 eV (+7,3 %) | **FAIL** |
| R1 | E_c(N) (riportato) | −309,3 contro −188,8 mEh | — |

Leggi scelte da sole (lascia-uno-fuori): lineare in x = 1/√(⟨r⟩ₐ⟨r⟩_b) per 2s², ¹P, ¹D, ³P(sp); 1/x per 2p² ³P;
10 manopole, tutte fit. O neutro sta dentro l'intervallo di allenamento (x 0,8–0,9 contro 0,52–1,13): non è
un'estrapolazione. Pezzi di O (mEh): 1s² −45,4; core-valenza −42,3; valenza −404,4 (2s² −96,6, 2s2p −163,6,
2p² −144,3).

**Perché fallisce (detto semplice)**: in uno ione con 2 soli elettroni di valenza, la coppia ha tutti gli orbitali del
guscio 2 liberi per "schivarsi" (2s² ↔ 2p², 2s2p ¹P ↔ 2p3d, …): da lì vengono i −50/−110 mEh delle coppie singoletto.
Nell'atomo vero quegli orbitali sono già occupati dagli altri elettroni, e ogni coppia non può più usarli. Sommare
coppie prese da ioni a 2 elettroni conta la stessa "via di fuga" più volte. Lo dice già A1, con 3 elettroni esatti:
+31/34 %. Con 6 elettroni di valenza l'errore cresce fino a quasi raddoppiare E_c. Anche il core-valenza preso dallo
ione tipo Li (5–8 mEh per elettrone) è troppo, per la stessa ragione (2s compatto). È la fisica nota delle teorie di
coppia "non chiuse" (Sinanoğlu): la correlazione non è additiva per coppie fra gusci aperti. Il metodo nuovo funziona
(V1: ECG e CI d'accordo entro 0,1 mEh), è la regola di somma che non regge.

Debiti: blocco HF congelato; trasferimento ione→atomo via x; pesi di configurazione media; ε(2p² ¹S) := ε(¹D)
(il ¹S di 2p² non è la radice più bassa: E_c da +37 a −20 mEh, solo riportato); core-valenza senza scala; additività.

**Export di O: NON fatto** (C1, C2 falliscono; anche M1 per la monotonia). Soglie non toccate.

Cosa resta aperto (proposta, non provata): la correlazione deve venire da sistemi con il guscio 2 già pieno quanto
nell'atomo, non da coppie isolate: ECG/CI a 3–4 elettroni nel blocco (A1 ora gira a K = 600 in ~2 ore), o una
correzione di esclusione per coppia imparata da A1 (dichiarata come nuova manopola, da pre-registrare).

Per rifare:
```
cd lab; export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MKL_NUM_THREADS=1
for g in g1 g2 g3 g4 g5 g6; do nohup ./tavg_f3_loop.sh $g >/dev/null 2>&1 & done   # basi in out/tavg_f3_base/ (non tracciate)
timeout 590 python tavg_f3_judge.py
```
