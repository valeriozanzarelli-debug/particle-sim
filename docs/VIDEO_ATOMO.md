# Video atomo: l'atomo che assorbe la luce e si svuota, al rallentatore (BOZZA)

> **Questo taglio è una BOZZA** (scritto anche a schermo). I dati sono l'esportazione con impronta
> `12322e1817` (sha256 di manifest.json + fotogrammi.csv, manifest del 13/10/2026 01:46 UTC). Il gruppo Debiti ha
> trovato uno smorzamento finto ×2,2 in `accoppiato.rs` (docs/EMISSIONE.md §7): la corsa verrà rifatta e lo
> svuotamento che si vede (velocità del calo di p_exc) **può cambiare**. Tutti i numeri di p_exc qui sotto sono di
> questa esportazione; il video li legge sempre dal file, non sono scritti a mano nel codice. Per rifare tutto con
> la nuova esportazione: `lab/vl_tutto.sh` (riconosce l'impronta nuova, cancella gli intermedi e rifà campi,
> palline, esami, fotogrammi, video, figura; `VL_ETICHETTA=""` toglie la scritta BOZZA).

Studio `vl_` del 13/10/2026. Il primo video lento della corsa accoppiata luce + atomo (1s esatto della regola
3-D, dati del commit addffd7, `/mnt/project-files/ripresa/dati/atomo_accoppiato/`).
Codice: `lab/vl_campi.py` (campi che muovono e colorano le palline), `lab/vl_palline.py` (moto delle palline, a
pezzi con checkpoint), `lab/vl_esami.py` (esami), `lab/vl_render.py` (disegno), `lab/vl_figura.py` (figura e
fermo-immagine), `lab/vl_video.sh` (montaggio). Pre-registrazione: `lab/out/vl_prereg.json` (**commit cf95ca6,
prima di muovere qualunque pallina**) + emendamento 1 `lab/out/vl_amend1.json` (commit edd687b, prima di rifare le
palline). Risultati: `lab/out/vl_campi.json`, `lab/out/vl_esami.json`,
`lab/out/vl_serie.json` (prima versione, col giro finto: `*_v1.json`).

Video e fermo-immagine: `/mnt/project-files/modello-standard/video_atomo_1080p.mp4`, `video_atomo_720p.mp4`,
`video_atomo_assorbe.png`, `video_atomo_oscilla.png`, `video_atomo_riemette.png`, `video_atomo_si_svuota.png`.

![figura](fig/video_atomo.png)

---

## 0. In breve

* **Cosa si vede.** A sinistra, grande, *dentro l'atomo*: 120 000 palline che sono l'elettrone. All'inizio sono
  messe una volta sola dove |ψ|² è grande, poi si muovono **solo** col flusso che serve a portare |ψ|² da un
  fotogramma all'altro (continuità): ferme in un livello puro, in moto solo quando la luce crea la sovrapposizione. A destra in alto la luce intorno all'atomo (scatola di 416 λ_C, l'atomo è il puntino al centro, in scala
  vera): palline di energia della luce. In basso a destra il tempo, la curva p_exc e il dipolo ⟨z⟩.
* **La storia in tre atti** (49 s di video, 24 fotogrammi al secondo):
  1. arriva l'impulso (bande rosa/verdi = E su e giù);
  2. la nuvola dell'elettrone viene scossa: **oro in alto e blu in basso, poi il contrario**, al ritmo
     ω₂₁ = 0,0752 (periodo 83,6 λ_C/c): è il dipolo che oscilla, cioè la sovrapposizione 1s + stato eccitato;
  3. l'impulso esce a destra; la sovrapposizione 1s + 2p creata dalla luce **continua a oscillare** e irraggia
     **come un'antenna**: dal centro escono **gusci di luce** (amplificata ×1000 in energia, dichiarato). **Non è
     decadimento spontaneo**: un 2p puro qui non emetterebbe (manca la luce quantizzata, Emissione fase 3, 57de9e9): p_exc scende da 0,120 (t = 220) a 0,064
     (t = 1600) mentre l'oscillazione si smorza.
* **Esami** (pre-registrati): le palline restano fedeli a |ψ|² per tutta la corsa senza mai essere rinate
  (correlazione minima **0,9996**), il dipolo letto dalle palline coincide con quello del calcolo (**0,999997**,
  scarto 0,015 λ_C), la luce ricostruita è coerente (0,96) e **l'energia della luce che esce dalla sfera di 40 λ_C
  è il 91 % dell'energia persa dall'atomo**. Un esame è **fallito come previsto** (p_exc letto dalle palline con la
  sola fase esportata: +0,05 di troppo). Dettagli in §4.
* **Validazione, non scoperta.** Tutto quello che si vede (dipolo che oscilla a ω₂₁, emissione a forma di
  ciambella attorno a z, svuotamento lento) è fisica nota; il valore è che **esce dalle palline della regola**,
  senza nulla disegnato.

## 1. Cosa si vede, pannello per pannello (come a un ventenne)

**Dentro l'atomo (sinistra).** Il nucleo è al centro (non disegnato). Le palline sono la "nebbia" dell'elettrone:
fitte dove la probabilità di trovarlo è alta. Il colore dice **come cambia la nuvola rispetto all'1s di partenza**:

* bianco-lavanda = come l'1s;
* **oro** = lì adesso c'è *più* elettrone di prima;
* **blu** = lì ce n'è *meno*.

Quando la luce scuote l'atomo, una parte dell'elettrone (al massimo il 12 %) passa allo stato eccitato. Le due
onde (1s ed eccitato) hanno energie diverse e quindi **battono**: la nuvola si sposta su e giù lungo z (la
direzione del campo elettrico dell'impulso). Si vede come **oro sopra / blu sotto**, poi **blu sopra / oro sotto**,
ogni 42 λ_C/c. Il centro diventa un po' blu in media: l'elettrone eccitato sta più largo, il centro si svuota un po'.
Prima che arrivi la luce le palline **stanno ferme**: in un livello puro questo reticolo non porta corrente
(scritto a schermo; §5 punto 3). Si muovono solo su e giù (e un po' dentro-fuori nelle code) quando c'è la
sovrapposizione. **Nessuna rotazione**: quella di Dirac attorno all'asse dello spin qui non c'è (debito aperto).

**La luce (destra in alto).** Ogni pallina è un pezzetto dell'energia della luce, (E² + B²)/2. Si muovono con la
velocità dell'energia (mai più di c). Colore = verso di E: rosa = E in su, verde-acqua = E in giù. La vista è
"tagliata": si disegna una fetta 0 < y < 72 λ_C, così si vede dentro. Prima l'impulso passa da sinistra a destra
come una fila di bande. Dopo che è uscito (t ≈ 240 vicino all'atomo) si vedono nascere dal centro **gusci** che
si allontanano a velocità c: è la luce del dipolo oscillante (1s + 2p guidati dalla luce), un'antenna, non un
decadimento spontaneo. Ogni pallina della luce dell'atomo vale 1/1000 di una
pallina dell'impulso: senza questo fattore non si vedrebbe (l'energia rimandata è 1,8·10⁻³ di quella
dell'impulso, sparsa su un guscio grande). Col tempo i gusci si fanno più radi: l'atomo si svuota.

**Il pannello in basso.** t in λ_C/c e in zeptosecondi (1 λ_C/c = 1,288·10⁻²¹ s; tutto il video sono 2 attosecondi),
il fattore di rallentamento (5·10¹⁹ volte nella parte lenta, 1,6·10¹⁹ nella parte veloce), la curva **p_exc** (quota
fuori dall'1s, **dal calcolo**, fotogrammi.csv) e il dipolo **⟨z⟩** (punti = palline, riga = calcolo).

## 2. Come è fatto (numeri)

### 2.1 Elettrone

* **Smussamento ~2 λ_C (dichiarato a schermo).** L'1s esatto del reticolo cubico è un miscuglio delle 4 valli del
  cammino: la sua densità ha una trama a scacchiera di periodo 4 celle = 2 λ_C (BATTIMENTO_1S §2.4). Densità e
  corrente sono smussate con una gaussiana σ = 2 celle = 1 λ_C (larghezza a metà altezza 2,35 λ_C). La componente
  di periodo 4 celle scende da 1,65 a 8,5·10⁻⁵ del totale (attenuazione 5·10⁻⁵).
* **Semina, una volta sola**: 120 000 palline, cella scelta con probabilità ∝ ρ_s·a³, posizione uniforme dentro la
  cella (seme 20261013). **Nessuna rinascita** in tutta la corsa.
* **Il flusso.** La corrente esportata (R†σR − L†σL per cella, mediata sui 16 tick) **non rispetta la continuità**:
  |∇·j| ≈ 0,30 contro |∂ρ/∂t| ≈ 0,03–0,06, correlazione tra −0,66 e +0,31. Non è la corrente sulle facce che usano i
  surfisti esatti (SURFISTI_ESATTI), è una lettura per cella che contiene una parte "di compressione" finta.
  **Prima versione** (pre-registrata): si teneva la j esportata e si correggeva col cambio più piccolo
  (J = j_s − ρ_s∇χ, ∇·(ρ_s∇χ) = ∇·j_s + ∂ρ_s/∂t): la correzione era l'84 % della j, e quello che restava era per
  l'81 % un giro attorno a z (le palline facevano ~4 giri, 0,12 di Bohm–Dirac). **Poi Debiti ha mostrato** che
  nessuno stato stazionario del reticolo 3-D porta corrente netta sui collegamenti (≤ 5·10⁻¹³ per coppia di tick):
  quel giro era una lettura finta della j per cella. **Emendamento 1** (`lab/out/vl_amend1.json`, prima di rifare):
  la j esportata non si usa più; il flusso è quello **di energia cinetica minima imposto dalla continuità**:
  J = −ρ_s∇χ con ∇·(ρ_s∇χ) = ∂ρ_s/∂t, sulle facce delle celle, per ogni intervallo tra due fotogrammi (gradiente
  coniugato, residuo < 10⁻⁷, ≤ 192 iterazioni). Risultato: rotazione netta attorno a z **7·10⁻⁸** (era 0,80);
  velocità media 0,017 c (era 0,05). Dentro la cella le palline usano i flussi sulle facce (RT0, come i surfisti)
  e la densità della cella lineare nel tempo: una nuvola uniforme in una cella resta uniforme. Residuo di
  continuità 4·10⁻⁴. Prima della luce le palline si spostano di 0,04 λ_C per fotogramma (rumore del float16
  esportato nelle code); spostamento radiale rms su tutta la corsa 0,43 λ_C (era 2,7 col giro finto).
* **Integratore**: punto medio, passo ≤ 0,5 λ_C/c; velocità tagliata a 0,95 c in 311 valutazioni su ~7·10⁸
  (code esterne); 0 urti contro il bordo della scatola (24 λ_C).

### 2.2 Luce

* Il dato esportato è solo E (64³, celle da 8 λ_C, media 2×2×2). **B si ricostruisce** con Faraday,
  ∂B/∂t = −∇×E (rotore spettrale, trapezi nel tempo), partendo da B₀ = x̂ × E₀ (l'impulso va lungo +x).
  Primo tentativo: dopo il passaggio dell'impulso resta un **B fermo e finto** grande 5–15 volte il campo
  dell'atomo (errore di integrazione + la parte non trasversa dell'impulso focalizzato). Siccome è fermo, si
  toglie la sua media sui fotogrammi 60–200. Dopo: B coincide con r̂ × E (radiazione) con correlazione media
  **0,963** (minimo 0,832); sull'impulso 0,999.
* **Palline dell'impulso**: 150 000, messe a t = 0 ∝ u = (E² + B²)/2 dentro |x|,|y|,|z| < 208 λ_C (fuori dalla
  spugna), mosse con v = 2E×B/(E²+B²) (Catmull-Rom nel tempo). Muoiono nella spugna.
* **Palline della luce dell'atomo**: nascono sulla sfera R₀ = 40 λ_C (circa mezza λ), dal momento in cui
  l'impulso ha lasciato la sfera (t = 240, dal flusso misurato), con probabilità ∝ flusso uscente (S·n̂)⁺; poi si
  muovono come le altre. Ognuna vale **1/1000** dell'energia di una pallina dell'impulso: **amplificazione ×1000
  in energia (×31,6 in campo), dichiarata a schermo**. Ne nascono 251 271 in tutto; vive insieme da ~80 000
  (t ≈ 400) a ~28 000 (t = 1600).
* La luce dell'atomo **durante** l'impulso non si separa: lì è mescolata all'impulso (è proprio l'interferenza
  che toglie energia all'impulso). Non si è sottratto un "impulso libero" perché l'esportazione non ha una corsa
  senza atomo: scelta dichiarata, si mostra solo dopo t = 240.

### 2.3 Video

* 1121 fotogrammi di disegno per 201 fotogrammi di dati: **12 per fotogramma di dati fino a t = 280** (assorbimento
  e primi battiti), poi raccordo, **4 dopo t = 360**; più 1,5 s fermi all'inizio e 3 s alla fine. 1920×1080 e
  1280×720, H.264.
* Disegno su CPU (Python, 4 core condivisi con altri lavori): ~10 s per fotogramma 1080p.

## 3. Che cosa è amplificato o smussato (tutto dichiarato a schermo)

| cosa | quanto | perché |
|---|---|---|
| densità e corrente dell'elettrone | smussate σ = 1 λ_C (~2 λ_C di larghezza) | togliere la scacchiera delle 4 valli del reticolo |
| luce dell'atomo | ×1000 in energia (×31,6 in campo) | è 1,8·10⁻³ dell'impulso, sparsa su gusci grandi |
| accoppiamento luce-elettrone (dai dati, non dal video) | K = 10 (e → 10e), tassi ×100 | vedere lo svuotamento in una corsa corta |
| colore dell'elettrone | pieno quando il cambio locale di densità è ~25 % | perché si legga il dipolo anche dove ci sono poche palline |
| corrente dell'elettrone | non usata: flusso minimo dalla continuità (emendamento 1) | la j esportata per cella non è il flusso del reticolo (§2.1) |

## 4. Esami (pre-registrati in `lab/out/vl_prereg.json`)

| esame | soglia | risultato | esito |
|---|---|---|---|
| E1 palline fedeli a |ψ|² smussata, 201 fotogrammi, senza rinascite | Pearson ≥ 0,97 sempre; TV ≤ TV(0)+0,05 | min **0,9996** (fot. 21); TV 0,022 → max 0,029 | **PASS** |
| E2 dipolo ⟨z⟩ dalle palline contro il calcolo | Pearson ≥ 0,98; rms ≤ 0,10 λ_C | **0,999997**; rms **0,015** λ_C; ampiezza max 1,169 contro 1,166 | **PASS** |
| E3a p_exc dalle palline = stesso stimatore sul reticolo | max diff ≤ 0,02 | 0,013 | **PASS** |
| E3b p_exc dalle palline contro p_exc vero | max diff ≤ 0,02 | 0,068 (media +0,050) | **FAIL (previsto)** |
| E4 B ricostruito coerente con la radiazione | media ≥ 0,90, min ≥ 0,80 | 0,963 / 0,832 | PASS (controllo già visto) |
| E5 bilancio di energia: luce uscita da R₀ / energia persa dall'atomo | in [0,67; 1,5] | **0,906** (3,69·10⁻³ contro 4,08·10⁻³ mc²) | PASS (controllo già visto) |
| E6 palline dell'impulso fedeli a u (fot. 10, 20, 30) | Pearson ≥ 0,90 | 0,985 / 0,984 / 0,980 | **PASS** |
| E7 dichiarazioni a schermo | leggibili nei fermo-immagine | smussamento, ×1000, K = 10, rallentatore, "dal calcolo", BOZZA, impronta dei dati, "in un livello puro... le palline stanno ferme" | PASS |

Note oneste sugli esami:

* E1 ed E2 passano **anche perché il flusso è costruito** per rispettare la continuità: misurano l'integratore e
  l'interpolazione, non la corrente esportata (che da sola non li passerebbe, §2.1). Prima versione (giro finto):
  E1 0,9995, E2 0,99999, E3a 0,017, E3b 0,071 — stessi esiti.
* E4 ed E5 erano già stati misurati prima della pre-registrazione (scritto nella pre-registrazione): sono controlli,
  non previsioni. E5 però è un buon numero di fisica: l'energia che la luce porta fuori è il 91 % di quella che
  l'atomo perde (p_exc da 0,1198 a 0,0640 × ω₂₁). Il 9 % mancante: luce ancora vicino all'atomo, la quota di
  p_exc che finisce nei livelli del multipletto del reticolo (non a ω₂₁), la sfera campionata su un reticolo grosso.
* E3b: con l'esportazione si ha solo la fase della componente dominante dello spinore; lo stimatore che la usa
  sbaglia di +0,04…0,06 già sul reticolo. Per questo **la curva p_exc a schermo è quella del calcolo**, ed è scritto.

## 5. Fallimenti e sorprese

1. **La corrente esportata non è un flusso** (l'84 % è compressione finta e il resto è un giro che il reticolo non
   ha). Con una prima correzione non pesata dalla densità le palline andavano sopra c in 75 000 valutazioni su 32
   fotogrammi e urtavano il bordo 26 000 volte; con la correzione pesata giravano per finta (punto 3). Serve esportare **la corrente sulle facce**
   della regola (quella dei surfisti esatti): debito 1.
2. **B non c'è nell'esportazione**: ricostruito con Faraday + un residuo fermo tolto (§2.2). Funziona (0,96) ma è
   una ricostruzione; esportare B è facile: debito 2.
3. **Il giro dell'1s: prima versione sbagliata, buttata (debito aperto, scritto a schermo).** Nella prima versione
   le palline dell'1s giravano attorno a z (~4 giri in tutta la corsa, 0,12 di Bohm–Dirac Zα c sinθ): veniva dalla
   j per cella, che non è il flusso del reticolo. Debiti (SURFISTI_ESATTI B5 e controlli di spin) ha mostrato che
   nessun livello puro dell'atomo 3-D porta corrente netta (≤ 5·10⁻¹³ per collegamento e coppia di tick): nella
   valle k = 0 l'elettrone gira a 0,87–0,97 di Zα c sinθ, ma le altre 3 copie del reticolo lo cancellano. Quindi nel
   video rifatto **le palline dell'1s stanno ferme** e si muovono solo nella sovrapposizione creata dalla luce: è la
   regola com'è oggi. Scritto a schermo: "in un livello puro questo reticolo non porta corrente: le palline stanno
   ferme (debito aperto: le copie del reticolo, da provare)". Nessuna rotazione disegnata o aggiunta. La rotazione
   di spin che Valerio vuole vedere arriverà solo quando il reticolo la porta davvero.
4. Il flusso "minimo" della continuità è **una scelta** (come Bohm sceglie la sua): tra tutti i flussi che portano
   |ψ|² giusto, quello senza giri. Se il reticolo avesse giri veri a divergenza zero, questo flusso non li vedrebbe.
   Per saperlo servono i flussi sui collegamenti nell'esportazione (debito 1).
5. Le palline della luce dell'atomo nascono solo fuori da R₀ = 40 λ_C: dentro (zona vicina) non ci sono palline di
   luce, anche se il campo c'è.
6. Righe deboli in più nella luce (se si vedono): secondo Emissione fase 3 vengono da un pezzo di elettrone strappato e
   tenuto dalla scatola periodica: **artefatto della scatola**. Nel video non le ho riconosciute a occhio.
7. Il colore dell'elettrone mostra soprattutto il **dipolo**; il fatto che l'eccitato sia più largo si vede solo
   come un centro un po' blu.

## 6. Debiti

0. **Prossimo taglio** (quando arriva l'esportazione nuova: 1s partito nella valle k = 0, nucleo imparato, senza lo
   smorzamento ×2,2): `lab/vl_tutto.sh` riconosce l'impronta, sceglie da solo il flusso (emendamento 2: la j esportata
   se la correzione per la continuità è ≤ 25 %) e scrive a schermo il giro dell'elettrone col rapporto con Dirac
   del passo usato (Debiti, surfisti fase 4: 0,867 / 0,910 / 0,943 / 0,958 a a = 0,5 / 0,4 / 0,3 / 0,25 λ_C),
   dichiarato come **errore del passo del reticolo**; nessuna scala per chiudere lo scarto.
0b. **Taglio BOZZA**: la corsa sotto ha uno smorzamento finto ×2,2 (EMISSIONE §7); da rifare con `lab/vl_tutto.sh`
   sulla nuova esportazione.
1. Flusso dell'elettrone ricostruito dalla sola continuità: servono i **flussi sui collegamenti** nell'esportazione.
1b. Rotazione di spin dell'1s assente sul reticolo (le 4 copie la cancellano): debito della regola, da provare.
2. B ricostruito con Faraday: esportare B (o F complesso) dal Laboratorio.
3. Smussamento σ = 1 λ_C per la scacchiera delle 4 valli: è un debito della regola 3-D (raddoppio dei fermioni),
   non del video.
4. Luce dell'atomo ×1000 e separata dall'impulso solo nel tempo (da t = 240): manca una corsa "senza atomo" da
   sottrarre.
5. K = 10, campo classico, ione Z = 60 (debiti della corsa, dal manifest).
6. Le palline sono una **lettura** (traccianti della densità e dell'energia), non le palline della regola: è
   un'interpretazione, come i surfisti.
7. Disegno su CPU: il porto su GPU è il lavoro del motore grafico (stesse palline, più palline e più risoluzione
   su una GPU seria; sulla Intel UHD meno palline e fotogrammi più piccoli, ma la stessa dinamica).

**Parametri fisici liberi: 0.** **Manopole estetiche (non toccano gli esami): 12** — grandezza delle palline
(elettrone, luce), bagliore (2), scia (2), soglia del colore, allungamento lungo E, guadagni di esposizione (2),
camere (2 percorsi), più la palette e lo spessore della fetta tagliata.

## 7. Come rifarlo

Un comando solo: `lab/vl_tutto.sh` (passi sotto, uno per uno).

```
python3 lab/vl_campi.py                  # ~6 min: campi in SCR (non in git) + lab/out/vl_campi.json
lab/vl_drive.sh elettrone; lab/vl_drive.sh luce     # palline, a pezzi da <= 8 min
python3 lab/vl_esami.py                  # esami
python3 lab/vl_render.py refs 0 50       # esposizione fissa per tutto il video
lab/vl_render_drive.sh 0 1121 1.0        # fotogrammi 1080p (0.6667 per 720p)
lab/vl_video.sh                          # mp4
python3 lab/vl_figura.py                 # figura e fermo-immagine
```
