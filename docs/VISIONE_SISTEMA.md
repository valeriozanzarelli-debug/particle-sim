# Visione: il sistema operativo per creare e simulare

Briefing di Valerio, 3 ottobre 2026. Questo documento fissa l'idea e la
prima struttura; i dettagli arriveranno nei briefing successivi.

## L'idea in una frase

Un sistema operativo unico (desktop, mobile, tablet, realtà aumentata) con al
centro un simulatore che va dal fotone all'oggetto industriale, guidato da un
agente AI serio che accompagna la persona dall'idea al progetto certificato.
Lo stesso sistema gestisce una comunità di 150–200 persone senza denaro.

## Le tre parti

### 1. Simulazione e creazione (il cuore)

Livelli di zoom, dal basso verso l'alto:

| Livello | Cosa si simula | Esempi di strumenti |
|---|---|---|
| Fotone | luce: interferenza, rifrazione, polarizzazione, entanglement | il nostro Photon Lab |
| Quantistico | particelle, atomi, chip quantistici, acceleratori | laboratorio quantistico |
| Chimico | molecole, materiali, reazioni | simulazione di materiali |
| Fisico | fluidi, vento, calore, strutture, circuiti, elettronica | CFD, elementi finiti, circuiti |
| Oggetto | macchine, CPU, RAM, computer, architettura | progettazione 3D, stampa 3D |
| Creativo | immagini, video, modelli | "Photoshop", video editing, generazione 3D (Meshy) |

Principi:
- **Un solo motore grafico** con dettaglio adattivo (tipo Nanite), condiviso
  da tutti i livelli: zoomando si cambia livello fisico, non programma.
- **Niente tabella infinita di strumenti**: l'agente propone solo quello che
  serve al momento, in finestre che compaiono dove servono.
- **Realtà aumentata** come modo principale di lavorare, lo schermo come
  alternativa.
- **Dall'idea al progetto certificato**: schizzo o descrizione, l'AI propone,
  si itera insieme, si arriva a un progetto ingegneristico verificato.

### 2. L'agente

Un agente serio alla base del sistema, che conosce gli strumenti e li usa per
la persona. Il costo di calcolo dell'agente è un problema aperto da risolvere
(da discutere).

### 3. Gestione della comunità

Sezioni del sistema per una comunità di 150–200 persone che:
- non usa denaro;
- gestisce cibo e risorse con un sistema condiviso;
- lavora circa tre ore al giorno nella manutenzione di una città
  tecnologicamente avanzata, dove l'automazione toglie il più lavoro possibile.

Dettagli da definire nei prossimi briefing.

## Scelte di architettura (proposta)

1. **Non reinventare i solutori che esistono già e sono validati.** A ogni
   livello di zoom il sistema orchestra strumenti affermati e liberi, e il
   nostro lavoro è collegarli e passare i dati tra i livelli:

   | Livello | Solutori esistenti da integrare (candidati) |
   |---|---|
   | Fotone, ottica | il nostro motore (validato su `docs/ESPERIMENTI_FOTONE.md`) |
   | Quantistico | QuTiP, Qiskit (chip quantistici), Geant4 (acceleratori, rivelatori) |
   | Chimico | PySCF, LAMMPS, modelli di reazione |
   | Fisico | OpenFOAM (fluidi, vento), CalculiX/FEniCS (strutture, calore), ngspice (circuiti) |
   | Oggetto | CAD parametrico (es. OpenCASCADE), slicing per stampa 3D |

2. **Il valore nostro sono i livelli di zoom**: regole ricavate dal livello
   sotto e verificate contro misure reali, con l'errore dichiarato a schermo.
3. **Un formato di scena unico** che tutti i livelli leggono e scrivono, così
   l'agente può passare un oggetto da un livello all'altro.
4. **Motore grafico unico** per tutte le piattaforme: Rust + wgpu (già in uso)
   gira su Windows, macOS, Linux, Android, iOS e nel browser.

## Primi tasselli (in ordine)

1. Fotone validato sugli esperimenti del catalogo (in corso).
2. Formato di scena unico e "guscio" dell'app: una finestra 3D, il tempo,
   la temperatura, e un pannello per l'agente.
3. Primo collegamento tra due livelli: il fotone che attraversa un materiale
   descritto dal livello chimico (indice di rifrazione, assorbimento).
4. Primo solutore esterno integrato (circuiti con ngspice, il più semplice).
5. Modulo comunità: primo briefing dedicato.

## Domande aperte

- Costo dell'agente AI: modello locale, cloud, o misto.
- Hardware per la realtà aumentata.
- Certificazione ingegneristica: quali norme e chi firma.
- Chiavi di servizi esterni (Meshy e altri): restano negli altri repository,
  mai copiate qui.
