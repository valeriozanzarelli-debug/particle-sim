# Fotone — interfaccia

## Interfaccia

La schermata di default è pulita: solo la scena 3D, una barra del tempo in
basso e il pulsante ☰ in alto a sinistra. Tutto il resto sta in un unico
pannello, nascosto all'avvio.

### Barra del tempo (in basso)

| Controllo | Funzione |
|---|---|
| ⏮ | riavvolge: il tempo scorre all'indietro attraverso la storia (di nuovo ⏮ = pausa) |
| ▶ / ⏸ | avanti verso il presente / pausa |
| ⏭ | avanti veloce (×10 sul fattore di velocità; di nuovo ⏭ = velocità normale) |
| ×fattore | velocità da ×0.1 a ×1000. A ×1 un fotone attraversa il banco in ~4 s. Da ×50 in su i fotoni vengono eseguiti interi, senza animazione (come il vecchio "turbo") |
| cursore storia | trascinandolo si scorre la storia conservata (mette in pausa) |
| ct = … λ · … fs | orologio della luce: distanza percorsa dalla luce dall'emissione, in lunghezze d'onda e in femtosecondi (1 unità = 500 nm, c = 299 792 458 m/s) |
| T … K | temperatura dell'ambiente (0–3000 K, default 293 K), con il picco di Wien |

Il mouse sopra la barra o il pannello non muove la camera.

### Pannello (Tab oppure ☰)

Scelta dell'esperimento (default: doppia fenditura), emissione di un fotone,
modalità continua, azzeramento dei conteggi, frecce del campo E, lunghezza
d'onda, numero di grani, waist, polarizzazione, parametri delle fenditure,
polarizzatori, impostazioni Bell, grafici e verifiche di validazione, e lo
stato della memoria del tempo.

### Temperatura (solo visiva, per ora)

I fotoni termici dell'ambiente sono disegnati come una foschia tenue: il
numero di punti segue la densità di fotoni di corpo nero (∝ T³, quindi a
temperatura ambiente sono pochissimi), il colore segue lo spettro di Planck
(picco di Wien λ = 2898 µm·K / T: infrarosso scuro sotto ~800 K, rosso →
arancio → bianco sopra). I fotoni della foschia viaggiano a c e seguono la
barra del tempo (fermi in pausa, all'indietro in riavvolgimento). Non
interagiscono ancora con il fotone studiato: l'accoppiamento con la fisica
verrà dopo (`src/photon/haze.rs`).

### Tempo reversibile: snapshot e memoria

La simulazione è deterministica (RNG `StdRng` con seme dentro `Sim`, hash
per grano, passi esatti di qualsiasi lunghezza). Quindi lo stato al frame n è
determinato da uno stato precedente più la lista delle operazioni applicate.

Salvare ogni frame costerebbe 60 000 grani × 64 B × 600 frame ≈ 2,3 GB.
La storia (`src/photon/timeline.rs`) conserva invece:

- un'operazione piccola per frame (`Op`: emetti / avanza di d / esegui n
  fotoni interi / completa il volo), pochi byte;
- un keyframe completo (clone di `Sim`: grani, coppia entangled, statistiche,
  stato dell'RNG) ogni k frame, con k = ⌈600 · dimensione_keyframe / 192 MB⌉
  (60 000 grani ≈ 4,3 MB → k = 14 → ~43 keyframe ≈ 190 MB);
- un keyframe dopo ogni operazione costosa da rigiocare (fotoni interi ad
  alta velocità), che è piccolo perché non c'è nessun fotone in volo.

Si tengono almeno 600 frame; i keyframe più vecchi vengono scartati quando la
finestra è coperta o si supera il budget di 192 MB, quindi la memoria resta
limitata in ogni caso.

- **Riavvolgere / trascinare il cursore**: si clona l'ultimo keyframe ≤ n e si
  rigiocano i passi uniti in un solo passo (un passo di qualsiasi lunghezza è
  fisica esatta e incontra gli stessi eventi con le stesse estrazioni casuali
  nello stesso ordine): stessa storia a meno dell'arrotondamento, al costo di
  circa un passo per frame.
- **Ripartire in avanti da un punto passato**: lo stato viene ricostruito
  esattamente (tutte le operazioni registrate, bit per bit) e le operazioni
  registrate vengono rieseguite: il futuro è identico a quello già visto.
- **Cambiare qualcosa nel passato** (un parametro, emettere un fotone,
  azzerare i conteggi): quel frame diventa il nuovo presente, il futuro
  registrato viene scartato e la simulazione riprende dal vivo da lì.

Il test `timeline::tests::rewind_and_replay_are_exact` verifica che
riavvolgere e riavanzare riproduca esattamente le posizioni dei grani e che le
anteprime coincidano entro 10⁻³.
