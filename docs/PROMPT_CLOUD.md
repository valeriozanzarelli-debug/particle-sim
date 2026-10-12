# Prompt per avviare la sessione in cloud

Copia il testo qui sotto (dal riquadro) come **primo messaggio** della sessione cloud aperta sul repo
`valeriozanzarelli-debug/particle-sim`, branch `master`. Il file `CLAUDE.md` alla radice del repo viene
letto automaticamente; questo prompt dice cosa fare per primo.

---

```text
Ciao, sono Valerio. Riprendi il progetto particle_sim da dove l'abbiamo lasciato sul portatile
(sessione del 9–13 ottobre 2026). Da ora il lavoro vive qui in cloud; il mio portatile fa solo da
laboratorio GPU tramite il ponte (docs/PONTE.md).

1. Leggi, in ordine: CLAUDE.md, docs/PIANO.md, docs/CONTESTO_CLOUD.md, docs/PONTE.md,
   docs/PASSAGGIO_CLOUD.md (stato esatto di ogni cantiere al momento del passaggio), e se esiste
   docs/MOLTI_ELETTRONI_PIANO.md (il piano già giudicato per gli atomi a più elettroni).
   Rispetta tutte le regole di CLAUDE.md, in particolare: niente fake, tutto emerge,
   pre-registrazione, una sola strada (regola zig-zag + agglomeratore automatico), commit con
   percorsi espliciti e firma Opus 5.5.

2. Prepara l'ambiente: installa i pacchetti Python che servono (numpy, scipy, torch CPU, numba,
   ripser, gudhi, matplotlib), verifica che la toolchain Rust compili `cargo build --bin atomo`
   se possibile (in cloud non c'è GPU né schermo: i test grafici si fanno col ponte). Controlla
   quanti core hai e adatta i thread degli studi. Manda una richiesta `ping` sul ponte
   (ponte/richieste/<data>_ping.json) per verificare che il portatile risponda.

3. Chiudi i cantieri lasciati aperti, seguendo docs/PASSAGGIO_CLOUD.md: per ognuno controlla
   cosa è committato, cosa manca, e riprendi dalla sua pre-registrazione senza rifare ciò che è
   già fatto e senza cambiare le soglie.

4. Poi avanza il piano, nell'ordine di docs/PIANO.md §5:
   a) atomi a più elettroni dalla regola agglomerata (elio prima: −2,903724 Eₕ, ionizzazione
      24,587 eV; poi litio, berillio, … fino al neon), seguendo docs/MOLTI_ELETTRONI_PIANO.md;
   b) idrogeno vero a scala di Bohr nel Laboratorio Atomo usando assets/atomo/agglomerato_zz.json
      e il reticolo sparso (assets/atomo/reticolo_spec.json), con zoom fino al nucleo;
   c) eccitazione → salto di orbitale → riemissione accoppiata elettrone–luce, e i preset di qualità
      fino a Ultra per una GPU da 16 GB, con il calcolo della regola su GPU (verifica sul portatile
      via ponte, azioni build/test/foto);
   d) molecole solo dopo la regola a molti elettroni.

5. Lavora come abbiamo sempre fatto: subagenti Opus in background, uno per studio; corse a pezzi
   con checkpoint; commit dopo ogni fase; documento in italiano semplice per ogni studio con
   numeri, fallimenti, debiti, parametri e una figura; aggiorna docs/PIANO.md quando uno studio si
   chiude. Usa il ponte per ogni cosa che richiede GPU o schermo, e mandami le immagini che
   tornano.

6. Quando mi scrivi: italiano informale, numeri veri, distingui validazione da scoperta, dimmi
   quando esce qualcosa di figo e anche quando non esce niente.

Inizia dicendomi in poche righe cosa hai trovato (stato dei cantieri, ambiente, ping del ponte) e
qual è il primo studio che lanci.
```
