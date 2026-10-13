# Video idrogeno: dal grande al piccolo (BOZZA)

Studio `vl_idro` del 13/10/2026. Zoom sull'idrogeno vero (Z = 1, α vera) dello studio IDROGENO_BOHR: da lontano la
nuvola 1s + 2p_z che oscilla su un periodo di Lyman-α, poi lo 1s fermo e lo zoom fino al protone, dove lavora la
zig-zag vera. Dati: `/mnt/project-files/ripresa/dati/idrogeno_bohr/` (manifest commit 918a71e, impronta
`cd7be29798`). Codice: `lab/vl_idro.py` (campi, palline, esami, disegno, video), `lab/vl_idro_drive.sh`.
Pre-registrazione: `lab/out/vl_idro_prereg.json` (commit 958dc5a, **prima di ogni pallina**). Risultati:
`lab/out/vl_idro_esami.json`, `lab/out/vl_idro_profilo.json` (controllo dopo, dichiarato).

Video e fermo-immagine: `/mnt/project-files/modello-standard/video_idrogeno_bohr_1080p.mp4`, `..._720p.mp4`,
`video_idrogeno_bohr_lontano.png`, `video_idrogeno_bohr_zoom.png`, `video_idrogeno_bohr_protone.png`.

---

## 0. In breve

* **Tre atti (30 s + pause, 24 fps).**
  1. **Lontano** (celle a0/8 = 6,6 pm, campo visivo 12 a0): 100 000 palline nate dove |ψ|² è grande, mosse dalla
     corrente della regola imparata. La nuvola va su e giù lungo z: è il dipolo della riga di Lyman-α (122 nm).
     Tre periodi (0,41 fs l'uno), rallentati ~10¹⁶ volte. Oro = più densità della media sul periodo, blu = meno.
  2. **Zoom** (8 s): si passa allo 1s fermo; le palline della griglia a0/16 (3,3 pm) svaniscono quando diventano più
     grandi del campo visivo, e compaiono quelle del blocco fine a0/512 (0,10 pm) attorno al protone. Campo visivo
     da 12 a0 a 0,22 a0 (da 635 pm a 11,6 pm): fattore 55. Esposizione automatica (dichiarata).
  3. **Sul protone** (32 tick della zig-zag, 6 fotogrammi per tick): la sfera r < 0,094 a0 del blocco fine; il
     colore è |R|² − |L|² sulla densità (quanto della pallina corre in un verso o nell'altro), amplificato ~500
     volte (pieno a 2·10⁻³).
* **La cosa onesta più importante.** Il manifest dice che la chiralità "cambia segno da un tick all'altro (il
  tremolio della zig-zag)". **Nei dati non è così**: è un dipolo quasi fermo lungo z (correlazione con z −0,95 a ogni
  tick, ±10⁻³ della densità), con una fluttuazione del 14 % e solo l'1,5–6,7 % delle celle che cambia segno fra due
  tick. Quindi a schermo si vede un "mezzo oro, mezzo blu" quasi fermo, e c'è scritto (esame H6 FAIL, previsto).
  Il tremolio a c della pallina della zig-zag *dentro* un tick non è nei dati esportati (si esporta un fotogramma
  per tick): debito.
* **Le palline del livello puro stanno ferme** (1s grosso e fine): coerente con Debiti (nessuna corrente netta nei
  livelli puri). Nessuna rotazione disegnata.

## 1. Come è fatto

* **Lontano.** ρ e j hanno solo le armoniche 0 e 1 del periodo (le altre pesano 2·10⁻⁹): si interpolano in modo
  esatto nel tempo per il colore. La j esportata ("velocità della regola dG/dk") è **quasi un flusso vero**: la
  correzione minima per rispettare la continuità è il **2,2 %** (mediana; massimo 11 %), contro l'84 % dell'atomo
  accoppiato. Rotazione netta attorno a z 2·10⁻⁸ (giusto: 2p_z ha m = 0). Palline sui flussi delle facce (RT0),
  passo 0,044 u.a., 351 velocità tagliate a 20 αc (code). La scatola ±6 a0 contiene il 92 % della probabilità
  (le code del 2p escono): le palline non escono dal bordo.
* **Zoom.** 100 000 palline ∝ ρ della griglia a0/16 (±2 a0, fuori dal blocco fine, probabilità 0,854) e 100 000 ∝ ρ
  del blocco fine (probabilità 0,00175). La luce di ogni pallina è proporzionale alla probabilità che porta: la
  stessa densità fa la stessa luce ai due livelli. La grandezza di una pallina è metà della sua cella.
* **Sul protone.** Densità quasi uniforme nel blocco (0,23–0,318 a0⁻³: il blocco è lungo 0,19 a0), con la cuspide
  al centro 0,317 ≈ 1/π.

## 2. Esami (pre-registrati)

| esame | soglia | risultato | esito |
|---|---|---|---|
| H1 palline lontane fedeli a ρ(t), 73 istanti su 3 periodi | Pearson ≥ 0,97; TV ≤ TV(0)+0,05 | min **0,996**; TV 0,097 → 0,103 | **PASS** |
| H2 dipolo dalle palline | Pearson ≥ 0,99; ampiezza entro 3 % | **0,9999999**; 0,7456 contro 0,7465 (griglia); file 0,7489 | **PASS** |
| H3 la j della regola è già quasi un flusso | correzione ≤ 25 % | 2,2 % (max 11 %) | **PASS** |
| H4 ritorno dopo 1, 2, 3 periodi | riportato | rms 0,11 / 0,19 / 0,28 a0 | — |
| H5 palline fini fedeli a ρ_fine; nucleo | Pearson ≥ 0,97; nucleo entro 10 % | Pearson **0,25**; nucleo 0,324 contro 0,316 (+2,7 %) | **FAIL** |
| H6 tremolio tick per tick (come dice il manifest) | ≥ 50 % celle cambiano segno | 1,5–6,7 % | **FAIL (previsto)** |
| H7 dichiarazioni a schermo | cambi di scala, amplificazione, palline ferme, taglio, esposizione | presenti | PASS |

* **H5 è un esame scritto male da me**: nel blocco la densità è quasi piatta (varia del 28 %), e 100 000 palline su
  96³ celle sono 0,1 per cella: la correlazione cella per cella misura il rumore di conteggio, non la fedeltà.
  Controllo fatto **dopo** (dichiarato, non cambia l'esito): profilo radiale in 9 gusci, palline contro dati, scarto
  massimo **2,7 %** (`lab/out/vl_idro_profilo.json`).
* H4: le palline non tornano esattamente dove erano dopo un periodo (0,11 a0 rms, cresce piano): il flusso di
  Bohm di una sovrapposizione non è periodico per le singole traiettorie anche se la densità lo è. Non è un errore.

## 3. Amplificato, smussato, scelto (tutto a schermo)

| cosa | quanto |
|---|---|
| colore della chiralità sul protone | pieno a |R|²−|L|² = 2·10⁻³ della densità (~×500) |
| colore lontano | pieno a ~25 % di cambio locale rispetto alla media sul periodo |
| esposizione | automatica fotogramma per fotogramma (la densità per pixel cambia di 10⁴ nello zoom) |
| vista | sfera r < 0,094 a0 nel blocco fine (il blocco è un cubo) |
| svanire delle palline grosse | fra campo visivo 1,2 a0 e 0,28 a0 (sono più grandi di quel che si guarda) |
| smussamento | nessuno |

## 4. Debiti

1. Il tremolio della zig-zag dentro un tick non è esportato (un fotogramma per tick; la direzione del passo non c'è).
2. La frase del manifest sul cambio di segno a ogni tick non è confermata dai dati (H6): da chiarire con chi l'ha
   esportato.
3. Dal manifest: nucleo fisso (massa ridotta), spin della regola grossa scartato, corrente spuria dalla parte
   dispari della regola imparata (j_grosso ≤ 0,0065, non usata: le palline dell'1s stanno ferme), vista lontana a
   a0/8 (1s −0,35 %).
4. Le palline sono una lettura (traccianti di |ψ|²), come i surfisti.
5. Esame H5 mal posto (vedi sopra).

**Parametri fisici liberi: 0.** Manopole estetiche: grandezza delle palline (metà cella), bagliore, scia (3 passi),
camera (percorso dello zoom), soglie di colore (2), svanire delle palline grosse.

## 5. Come rifarlo

```
python3 lab/vl_idro.py campi      # a pezzi da <= 8 min: rilanciare finché scrive campi.json
python3 lab/vl_idro.py palline
python3 lab/vl_idro.py esami
lab/vl_idro_drive.sh              # fotogrammi
python3 lab/vl_idro.py video      # mp4 + fermo-immagine
```
