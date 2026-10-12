# particle_sim — istruzioni per ogni sessione di Claude

Progetto di Valerio: un simulatore in cui **tutto emerge da palline con una regola semplice**
(luce → elettrone → atomo → molecola → materiali → oggetti → universo). Prodotto finale: un motore
"come Unreal" con fisica vera sotto, dove le impostazioni grafiche sono i livelli di agglomerato.

## Leggere all'avvio, in quest'ordine
1. `docs/PIANO.md` — regole, stato, cantieri aperti, prossimi passi.
2. `docs/CONTESTO_CLOUD.md` — indice di tutti gli studi, numeri, debiti, strumenti, lezioni.
3. `docs/PONTE.md` — come usare il portatile di Valerio come laboratorio GPU.
4. Il documento dello studio su cui si lavora (`docs/<STUDIO>.md`) e la sua pre-registrazione
   (`lab/out/<studio>_prereg.json`).

## Regole vincolanti
- **Niente fake**: ingressi ammessi = costanti fisiche, cariche e masse dei nuclei (e costanti di
  natura come la massa del pione). Niente tabelle di materiali, niente codici di chimica esterni.
  Gli esperimenti servono solo come esame.
- **Tutto emerge**: ogni cosa messa a mano è un **debito**, scritto nel documento dello studio.
- **Pre-registrazione**: prima di ogni corsa di produzione si committa
  `lab/out/<studio>_prereg.json` con esami, soglie, configurazioni. Cambi dopo: dichiarati.
- **Manopole contate**; un fit non è una previsione. **Fallimenti dichiarati** come i successi.
- **Una sola strada**: la regola zig-zag (corri a c, poi ruota). Si sale di scala solo con
  l'**agglomeratore automatico** (blocca, impara, controlla, ripeti). Molecole solo dopo la regola
  a molti elettroni. Le strade storiche (fotone a due ruoli, DMC separato, forme disegnate) non si
  riaprono senza che Valerio lo chieda.
- **HD**: la risoluzione è fisica (es. ≥ 24–32 celle per λ). Grafica pensata per una GPU seria da
  16 GB; il portatile serve per controllare (anche lento o a fotogrammi singoli).
- **Grafica**: palline morbide che si muovono (stile del vecchio Laboratorio Fotone), forme che
  nascono da densità e corrente; niente linee/fili, niente nuvole ferme, niente geometria aggiunta.

## Come si lavora
- Lavoro lungo → subagenti **Opus** in background, uno per studio, con file dedicati
  (`lab/<prefisso>_*.py`, `docs/<STUDIO>.md`, `lab/out/<prefisso>_*`).
- Ogni studio: pre-registrazione → corse a pezzi con checkpoint (≤ 10 min per script) →
  documento in **italiano semplice** ("come a un ventenne") con numeri, fallimenti, debiti,
  numero di parametri, una figura.
- Commit dopo ogni fase finita. `git add` con **percorsi espliciti**; PNG/JSON di `lab/out/` con
  `git add -f`; mai i file grandi (`.npz`, `.pt`, `.npy`, `.pkl`). `git pull --rebase --autostash`
  poi push; **mai force-push**. Non tracciare file di log scritti di continuo (rompono l'autostash).
- Firma dei commit: `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Aggiornare `docs/PIANO.md` (stato e cantieri) quando uno studio si chiude o cambia.

## Laboratorio GPU (ponte)
Il cloud non ha GPU né schermo. Per compilare/provare il Laboratorio Atomo sul portatile: scrivere
`ponte/richieste/<id>.json` (azioni `ping`, `build`, `test`, `foto`), push, poi leggere
`ponte/risultati/<id>/` dopo il pull. Dettagli in `docs/PONTE.md`.

## Con Valerio
Italiano informale ("bro/fra"), spiegazioni semplici ma oneste. Avvisarlo quando esce qualcosa di
bello **e** quando non esce niente, sempre con i numeri. Non gonfiare: distinguere validazione
(ritrovare numeri noti) da scoperta. Le sue intuizioni fisiche vanno tradotte in regole e test, dicendo
cosa corrisponde a fisica nota, cosa è nuovo, cosa è sbagliato. Chiavi di servizi esterni (es. Meshy)
non vanno mai copiate in questo repo.
