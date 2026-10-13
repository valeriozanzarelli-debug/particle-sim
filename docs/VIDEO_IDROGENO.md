# Video idrogeno: dal grande al piccolo

Studio `vl_idro` del 13/10/2026. Zoom sull'idrogeno vero (Z = 1, α vera) dello studio IDROGENO_BOHR: da lontano la
nuvola 1s + 2p_z che oscilla su un periodo di Lyman-α, poi lo 1s fermo e lo zoom fino al protone, dove lavora la
zig-zag vera. **Dati della fase 2** (`/mnt/project-files/ripresa/dati/idrogeno_bohr/fase2/`, manifest → `fase2`,
commit 056007a, docs/IDROGENO_BOHR.md §11): regola T con parità esatta (mescolamento s–p 2,2·10⁻⁹), colla con reflux,
Lyman-α a 122,21 nm, blocco fine 7³ celle grosse a a0/256. Impronta dei dati `85de02bc47` (sha256 di manifest.json +
fase2/dipolo_lontano.txt). Codice: `lab/vl_idro.py`, `lab/vl_idro_drive.sh` (tutto con un comando).
Pre-registrazione: `lab/out/vl_idro_prereg.json` (commit 958dc5a) + **emendamento 1 `lab/out/vl_idro_amend1.json`
(commit 580668c, prima di ogni pallina sui dati della fase 2)**. Risultati: `lab/out/vl_idro_esami.json`.

Video e fermo-immagine: **in rendering** (verranno copiati in `/mnt/project-files/modello-standard/` appena finiti;
questo paragrafo sarà aggiornato con i nomi verificati).

---

## 0. In breve

* **Aspetto.** Ogni pallina è un **bagliore morbido** che si somma agli altri (come la scena zig-zag del Laboratorio);
  **niente scie, niente righe**: il moto si vede perché i bagliori si spostano. La fase dell'onda qui non è
  esportata, quindi il colore lontano è il cambio di densità (oro = più della media sul periodo, blu = meno) e sul
  protone è la chiralità |R|² − |L|².
* **Tre atti (30 s + pause, 24 fps).**
  1. **Lontano** (celle a0/8 = 6,6 pm, campo visivo 12 a0): 100 000 palline nate dove |ψ|² è grande, mosse dalla
     corrente della regola (corretta al minimo per la continuità: **0,6 %**). La nuvola va su e giù lungo z: è il
     dipolo della riga di Lyman-α (122,21 nm). Tre periodi (0,41 fs l'uno), rallentati ~10¹⁶ volte.
  2. **Zoom** (8 s): si passa allo 1s fermo; le palline della griglia a0/16 (3,3 pm) svaniscono quando diventano più
     grandi del campo visivo, e compaiono quelle del blocco fine a0/256 (0,21 pm) attorno al protone. Campo visivo da
     12 a0 a 0,44 a0 (da 635 pm a 23 pm). Esposizione automatica (dichiarata).
  3. **Sul protone** (16 tick della zig-zag, 12 fotogrammi per tick, 1 tick = 6,9·10⁻²² s): la sfera r < 0,21 a0 del
     blocco fine; colore = |R|² − |L|² sulla densità, pieno a 4·10⁻³ (amplificato ~250 volte, dichiarato).
* **Scritte di onestà a schermo** (richieste e vere):
  - lo 1s ibrido **non converge** con la grandezza del blocco: **−0,043 ± 0,004 %** (§11.2 di IDROGENO_BOHR);
  - la corrente residua **1,7·10⁻³** dei dati viene dal transitorio (si parte dallo stato della griglia grossa, non
    da quello dell'ibrido): **non è fisica** e non si usa. Le palline dello 1s **stanno ferme** (in un livello puro
    non c'è corrente): quindi nella vista dello 1s nulla si muove per colpa di quella corrente;
  - la chiralità nei dati è un dipolo quasi fermo lungo z, **non** un tremolio che cambia segno a ogni tick (H6).

## 1. Come è fatto

* **Lontano.** ρ e j hanno solo le armoniche 0 e 1 del periodo: interpolazione esatta nel tempo per il colore. La j
  esportata è **quasi un flusso vero**: la correzione minima per la continuità è lo **0,6 %** (fase 1: 2,2 %; con la
  regola T a parità esatta è migliorata). Rotazione netta attorno a z 3·10⁻¹⁰ (giusto: 2p_z ha m = 0). Palline sui
  flussi delle facce (RT0), 278 velocità tagliate a 20 αc (code).
* **Zoom.** 100 000 palline ∝ ρ della griglia a0/16 fuori dalle 7³ celle coperte dal blocco fine (probabilità 0,838)
  e 100 000 ∝ ρ del blocco fine (probabilità 0,0176). La luce di ogni pallina è proporzionale alla probabilità che
  porta: la stessa densità fa la stessa luce ai due livelli.
* **Sul protone.** Densità nel blocco 0,150–0,318 a0⁻³, al centro 0,314 (1/π = 0,318).

## 2. Esami (pre-registrati; H5 cambiato dall'emendamento 1 prima delle corse)

| esame | soglia | risultato | esito |
|---|---|---|---|
| H1 palline lontane fedeli a ρ(t), 73 istanti su 3 periodi | Pearson ≥ 0,97; TV ≤ TV(0)+0,05 | min **0,996**; TV 0,097 → 0,103 | **PASS** |
| H2 dipolo dalle palline | Pearson ≥ 0,99; ampiezza entro 3 % | **0,9999997**; 0,7429 contro 0,7471 (griglia); file 0,7495 | **PASS** |
| H3 la j della regola è già quasi un flusso | correzione ≤ 25 % | **0,6 %** | **PASS** |
| H4 ritorno dopo 1, 2, 3 periodi | riportato | rms 0,014 / 0,016 / 0,016 a0 | — |
| H5 palline fini contro ρ_fine: profilo radiale, 9 gusci | scarto max ≤ 5 % | **4,8 %** (vecchio Pearson cella per cella 0,30, riportato) | **PASS** (di poco) |
| H6 tremolio tick per tick (come dice il manifest) | ≥ 50 % celle cambiano segno | 0–1,2 % | **FAIL (previsto)** |
| H7 dichiarazioni a schermo | scale, amplificazione, palline ferme, corrente spuria, non convergenza, taglio, esposizione | presenti | PASS |

* H4: con la regola T le palline tornano quasi dove erano dopo un periodo (0,014 a0, era 0,11 nella fase 1): il flusso
  ora è quasi periodico anche per le singole traiettorie.
* H5 passa di poco (4,8 % contro 5 %): è rumore di conteggio nei gusci interni (pochi bagliori vicino al centro).

## 3. Amplificato, smussato, scelto (tutto a schermo)

| cosa | quanto |
|---|---|
| colore della chiralità sul protone | pieno a \|R\|²−\|L\|² = 4·10⁻³ della densità (~×250) |
| colore lontano | pieno a ~25 % di cambio locale rispetto alla media sul periodo |
| esposizione | automatica fotogramma per fotogramma (la densità per pixel cambia di 10⁴ nello zoom) |
| vista | sfera r < 0,21 a0 nel blocco fine (il blocco è un cubo) |
| svanire delle palline grosse | fra mezzo campo visivo 1,2 a0 e 0,35 a0 |
| smussamento | nessuno |
| scie | nessuna (Valerio, 13/10) |

## 4. Debiti

1. Il tremolio della zig-zag dentro un tick non è esportato; la frase del manifest sul cambio di segno a ogni tick
   non è confermata dai dati (H6).
2. Lo 1s ibrido non converge con la grandezza del blocco (−0,043 ± 0,004 %), picco al blocco 7 non capito (§11.2).
3. Corrente spuria 1,7·10⁻³ dal transitorio: partire dall'autostato dell'ibrido (prossimo passo di IDROGENO_BOHR).
4. La fase dell'onda non è esportata: il colore non può essere la fase come nella scena del Laboratorio.
5. Dal manifest: nucleo fisso (massa ridotta), spin della regola grossa scartato, vista lontana a a0/8 (1s −0,385 %),
   massa del cammino fine (+0,034 % nel 2p).
6. Le palline sono una lettura (traccianti di |ψ|²), come i surfisti.

**Parametri fisici liberi: 0.** Manopole estetiche: grandezza dei bagliori (metà cella), bagliore largo, camera
(percorso dello zoom), soglie di colore (2), svanire delle palline grosse.

## 5. Come rifarlo

```
lab/vl_idro_drive.sh      # campi (a pezzi), palline, esami, fotogrammi, mp4 e fermo-immagine; riprende da dove era
```

## 6. Storia

Un primo taglio sui dati della fase 1 (blocco a0/512, 32 tick, impronta `cd7be29798`) non è mai stato montato: superato
dai dati della fase 2. Esami di allora: H1 0,996, H2 0,9999999, H3 2,2 %, H5 Pearson 0,25 (mal posto; profilo 2,7 %),
H6 1,5–6,7 % (FAIL previsto).
