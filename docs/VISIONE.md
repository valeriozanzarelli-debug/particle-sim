# Visione: un "Blender atomico"

Appunti di progetto, 30 settembre 2026. Scopo: fissare dove vogliamo arrivare
e in che ordine, così ogni pezzo costruito adesso serve anche dopo.

## Cosa vogliamo costruire

Un simulatore di materiali in cui si può costruire e provare qualsiasi cosa —
un materiale, una mela, un chip, un computer intero — con fedeltà fisica, dal
livello quantistico in su. Lo strumento di misura è la luce: si simula come la
luce interagisce con la materia, fotone per fotone dove serve.

## L'idea centrale: livelli di zoom

"A questo livello di zoom cento palline si muovono come venti."

In fisica questa idea esiste e ha un nome: **coarse-graining**, e la teoria che
dice come cambiano le regole quando si cambia scala è il **gruppo di
rinormalizzazione**. È lo stesso principio di Nanite per le mesh (dettaglio
solo dove lo guardi), applicato alle regole fisiche invece che ai triangoli.

Ogni livello ha le sue palline e la sua regola locale. La regola di un livello
non si inventa: si **ricava** dal livello sotto e si **verifica** contro misure
reali.

| Livello | "Pallina" | Regola locale | Si verifica contro |
|---|---|---|---|
| 0 | grano del fotone / dell'elettrone | ruota e trasla, fase dai vicini | ottica, doppia fenditura, Bell |
| 1 | atomo | livelli di energia, emissione e assorbimento, vuoto | spettri atomici |
| 2 | legame, molecola | forze tra atomi, reazioni | energie di legame, spettri molecolari |
| 3 | ammasso (pochi pixel) | assorbimento e diffusione della luce, diffusione chimica, calore | spettri di materiali, curve di reazione |
| 4 | oggetto | meccanica, calore, ottica dell'insieme | misure macroscopiche |
| 5 | mondo, pianeta | gravità, atmosfera, clima | osservazioni |

## Il sistema che osserva gli ammassi

Un osservatore (per ora Claude, poi un modello dedicato) che:
1. fa girare la simulazione fine su un campione piccolo;
2. ne ricava la regola dell'ammasso (quante palline servono, come si muovono);
3. misura l'errore rispetto alla simulazione fine e ai dati reali;
4. accetta la regola solo se l'errore sta sotto una soglia dichiarata.

Senza il punto 3 il risultato è grafica, non simulazione.

## Limiti da tenere presenti

- **Entanglement:** non può venire da regole solo locali (teorema di Bell).
  Serve un legame non locale tra le palline gemelle.
- **Computer quantistici:** simulare esattamente molti qubit costa memoria che
  raddoppia a ogni qubit. Oltre qualche decina di qubit servono
  approssimazioni, e bisogna dire quali.
- **Fedeltà:** a ogni salto di livello si perde informazione. Il programma deve
  mostrare a quale livello sta lavorando e con quale errore stimato.
- **Hardware:** la GPU integrata Intel basta per sviluppare, non per scene grandi.

## Ordine di lavoro

1. Fotone: regola locale per pallina, validata su tutti gli esperimenti.
2. Elettrone: stessa idea con massa (pacchetto che si allarga, effetto tunnel).
3. Atomo di idrogeno: elettrone legato, livelli di energia, emissione di un fotone.
4. Vuoto quantistico dove conta: emissione spontanea, Lamb shift, Casimir.
5. Legami e molecole semplici (H2, acqua).
6. Primo salto di livello: ammassi, con l'osservatore che ricava le regole.
7. Primo oggetto: la fetta di mela che imbrunisce (ossigeno, temperatura, tempo).
8. Scelta della piattaforma e costruzione dell'editor.

La scelta della piattaforma viene dopo il punto 6: prima di allora non si sa
quali strutture dati servono davvero.
