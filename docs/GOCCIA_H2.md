# Goccia di H₂: test frattale del ciclo dei livelli — PARZIALE

Studio `aggmol_frat_`. Idea di Valerio: una regola è una regola solo se è uguale a ogni zoom; si sceglie solo quanto tiene
ogni passo. La salita è molecola → pallina-molecola → ammasso → pallina-blob → goccia, con il ciclo di
`CICLO_LIVELLI.md`. **L'H₂ liquido sta al posto dell'acqua** finché l'ossigeno non passa la regola a molti elettroni.
I nuclei sono classici, quindi legano troppo (manca il moto di punto zero): va dichiarato in ogni esame.

**Stato: PARZIALE.** Ci sono solo corse di sviluppo. La pre-registrazione **non** è ancora committata come JSON
(lo script è una bozza) e non è partita nessuna corsa di produzione. Fermato per la pausa di Valerio.

## Fatto finora (sviluppo, dichiarato)

**Livello fine anisotropo completo.**
- Calcolata la quinta orientazione S (parallele, inclinate di 45°) con la macchina della fase 3
  (`lab/aggmol_frat_aniso.py` → `lab/out/aggmol_frat_aniso.json`, committato).
- Con le 5 orientazioni la forza anisotropa di rango 2 è determinata del tutto (`lab/aggmol_frat_pot.py`).
- C9 (tre corpi, Axilrod–Teller) = 40,14 u.a., dagli stessi pseudo-stati.

**Passo 1, molecola → pallina (sviluppo, R = 6,5 a₀).**
- Pallina V₀: −33,41 K. Rotori quantistici para (J pari): −33,50 K (−0,26 %).
- Rotori classici a 20 K: −35,60 K. Minimo classico (forma T): −52,13 K.
- Trimero equilatero: rotori quantistici −100,50 K contro pallina −100,24 K.
- Cioè: per l'H₂ para la pallina sferica è quasi esatta. Un livello fine *classico* sbaglierebbe del 7 % e oltre.

**Motore MD** (`lab/aggmol_frat_md.py`).
- Minimi di Mackay 13 / 55 / 147: −112,7 / −164,2 / −190,6 K per molecola.
- Cubo di 1000 molecole:
  - a 20 K si arrotonda in una goccia (asfericità 0,004–0,011) con un po' di vapore, Epot/N ≈ −156 K;
  - a 15 K: −176 K, niente evapora.
- Costo: 8 ms per passo.

**Livello blob (il punto dove si rompe).**
- Blob = gruppi di ~13 molecole scelti da k-means. Forza fra blob imparata per force matching.
- Il modello a coppie spiega solo **l'1–4 %** della forza sui blob dell'ammasso di 147.
- Con una finestra temporale (10–25 fotogrammi) la convalida incrociata dà errore ≥ 1, cioè nessun potere predittivo.
- Le gocce grosse collassano (Rg dei blob 2–8 Å contro ~15 della fine).
- Lezione provvisoria: in un liquido i centri dei blob **non sono variabili lente**, perché le molecole cambiano blob
  in ~5 ps. Il ciclo vuole "conservate + lente": qui sono densità e quantità di moto (idrodinamica), non i gruppi.

## Per riprendere
1. Decidere la variabile del livello blob. Due candidati:
   - densità e quantità di moto per cella (idrodinamica fluttuante);
   - inversione di Boltzmann della g(r) dei blob.
   Poi correggere `lab/aggmol_frat_prereg.py` (esami L1/L2) e fare il commit di `lab/out/aggmol_frat_prereg.json`
   **prima** delle corse.
2. Corse (un thread, ≤ 10 min ciascuna, con checkpoint):
```
cd lab
python aggmol_frat_prereg.py                      # poi git add -f lab/out/aggmol_frat_prereg.json, commit
python aggmol_frat_run.py step1                   # passo 1: dimeri, trimeri, B2
python aggmol_frat_run.py clusters                # passo 2: 13/55/147 a 5-20 K
python aggmol_frat_run.py learn                   # livello blob (da rifare con la variabile scelta)
python aggmol_frat_run.py drop 20                 # ripetere finché non scrive più "checkpoint"; poi drop 15
python aggmol_frat_run.py coarse 20; python aggmol_frat_run.py coarse 15
python aggmol_frat_run.py analyze                 # esporta anche lab/out/aggmol_frat_drop*_frames.json per i bagliori
```
   Il monomero `lab/out/aggmol_h2h2_mono.npz` (non nel repo) si rifà con `python aggmol_h2h2.py mono`.
