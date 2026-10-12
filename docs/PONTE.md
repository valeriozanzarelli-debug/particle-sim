# Il ponte cloud ↔ portatile

Il lavoro vive in una sessione **cloud** (codice, studi Python, documenti, commit). Il cloud però
non ha scheda video né schermo. Il **portatile** fa da laboratorio: compila, esegue i test grafici,
scatta le "fotografie" del Laboratorio Atomo e misura i tempi sulla sua GPU. Il ponte è **GitHub**.

## Come funziona

1. La sessione cloud scrive una richiesta in `ponte/richieste/<id>.json` e fa push.
2. Sul portatile gira `python ponte/ponte.py`. Ogni 90 secondi scarica il repo in una **copia
   separata** (`../ponte_lab`, accanto al repo), esegue le richieste nuove e fa push dei risultati
   in `ponte/risultati/<id>/` (`esito.json`, `log.txt`, `foto.png`).
3. La sessione cloud fa `git pull` e legge i risultati.

La copia separata serve a non toccare mai la cartella in cui lavori a mano.

## Azioni ammesse (e solo queste)

| azione | cosa fa | argomenti (`"args"`) |
|---|---|---|
| `ping` | revisione git, scheda video, CPU, RAM | nessuno |
| `build` | `cargo build --bin atomo` | nessuno |
| `test` | `cargo test --bin atomo` | `filtro` (nome del test), `ignored` (true/false) |
| `foto` | compila e lancia il Laboratorio Atomo in modalità `--foto`, restituisce l'immagine | `scena`, `zoom`, `stati`, `qual`, `qualita`, `polar`, `dir`, `w0`, `lambda-m`, `strati`, `frames`, `spara-al`, `misure`, `ritmo`, `yaw`, `pitch`, e i flag `surfisti`, `collassa`, `pannello` |

Ogni valore viene controllato (scelte fisse o espressioni regolari) prima di essere usato. Il ponte
**non** esegue comandi liberi scritti nel repo. Resta vero che compilare ed eseguire il codice del
repo significa eseguire quel codice: va puntato solo al nostro repository.

## Esempi di richiesta

```json
{"azione": "foto", "args": {"scena": "orbitali", "stati": "2p+3d", "frames": 300}}
```

```json
{"azione": "test", "args": {"filtro": "orbitali_livelli"}}
```

Il nome del file è l'id: lettere, cifre, `_` e `-`, massimo 64 caratteri (es.
`2026-10-13_foto_orbitali.json`). Una richiesta già eseguita (cartella con `esito.json`) non
viene ripetuta: per rifarla si usa un id nuovo.

## Avvio sul portatile

```bash
python ponte/ponte.py
```

Si ferma con Ctrl+C. Una sola passata: `python ponte/ponte.py --una`. Al primo avvio clona il repo
in `../ponte_lab` e la prima compilazione richiede parecchi minuti.

## Limiti

- Il portatile deve essere acceso, con il ponte in esecuzione.
- I tempi misurati sono quelli della scheda Intel UHD del portatile, non di una scheda seria.
- Una richiesta lunga (es. `foto` a qualità ultra) blocca le successive finché non finisce
  (timeout 1 ora per azione).
