# Il reticolo dell'universo: palline solo dove c'è qualcosa, a tutte le scale

**L'idea (di Valerio).** "Il nostro reticolo deve poter contenere tutta la Via Lattea o l'universo conosciuto; per
vedere intera un'onda radio devi zoomare molto indietro, ma noi la vogliamo vedere."

**La risposta in una riga.** Si può, se il reticolo è **rado** (le celle esistono solo dove il campo non è zero) e
**a livelli** (celle di 1, 2, 4, 8, ... volte più grandi), con **la stessa regola della luce su ogni livello**.
Il costo per *tenere* un'onda non dipende da quanto è grande l'universo né da quanto è grande l'onda: dipende solo da
**quante lunghezze d'onda** è lunga. Il vero muro non è la memoria, è **il tempo**: un raggio gamma fa un miliardo
di miliardi di tick mentre un'onda radio ne fa uno.

Pre-registrazione: `lab/reticolo_prereg.py` → `lab/out/reticolo_prereg.json`, **committata (de7b1af) prima delle
corse**, insieme al motore `lab/reticolo_core.py`. Corse: `lab/reticolo_r2.py` (interfaccia), `lab/reticolo_r2a.py`
(velocità), `lab/reticolo_r3.py` (il mondo con gamma e radio), `lab/reticolo_r4.py` (quanto costa), figure
`lab/reticolo_figure.py`, formato per il Rust `lab/reticolo_spec.py` → `assets/atomo/reticolo_spec.json`.

## 1. Come è fatto, spiegato semplice

**La regola della luce non ha un metro suo.** La pallina di luce (spin 1, porta F = E + icB, `docs/ZIGZAG.md`) a ogni
tick si sposta di una cella e gira il suo stato: non sa quanto è grande la cella. Un fotone radio di 1 m e un fotone
di 1 µm sono **lo stesso disegno di palline**, cambia solo il righello. Per questo possiamo usare la **stessa regola**
su celle grandi h₀, 2h₀, 4h₀, ..., 2ᵏh₀ (il "livello k"), con il tick lungo quanto la cella (dt = h/c).
Lo abbiamo verificato al bit: lo stesso pacchetto messo al livello 0, 1, 10 o 20 dà **gli stessi identici numeri**
(R2a).

**Rado: solo i blocchi che servono.** Ogni livello è un insieme di **blocchi** di 16×16 celle (8×8×8 in 3-D) con un
indirizzo intero. Dove non c'è un blocco il campo è **zero esatto** e non costa niente. Un blocco del livello k+1
contiene esattamente 2×2 (2×2×2) blocchi del livello k: è un albero (quadtree/octree).

**Quale livello tiene cosa.** Un'onda con lunghezza λ "vive bene" sul livello dove ha **16–32 celle per λ**
(il suo livello di casa). Un livello più fine serve solo dove il campo ha dettagli che il livello più grosso non sa
tenere: si misura con la differenza fra il campo fine e quello interpolato dal grosso (se supera una soglia, il
blocco resta, più un blocco di margine). Il livello più grosso tiene i blocchi dove |F| supera una soglia.
Ogni livello fine deve stare **dentro** il livello sopra (con 2 blocchi di margine): attorno a un oggetto piccolo,
tutti i livelli fino a quello di un oggetto grande formano un **camino** di blocchi (si vede nella figura sotto).

**Il tempo.** Il livello k+1 fa un tick mentre il k ne fa due (un tick grosso dura il doppio). Si fa come
Berger e Oliger (1984): prima avanza il grosso, poi il fine lo raggiunge con due tick.

**La cucitura fra due livelli** (come in `ZOOM_NODO.md`):
1. dove il fine esiste (nella sua parte interna) il grosso **copia** il fine (interpolazione di quart'ordine a metà
   strada, pesi −1, 9, 9, −1 su 16);
2. la **fascia** più esterna del fine (4 celle fini = 2 grosse) viene **riscritta** con l'interpolazione cubica del
   grosso; fra una cucitura e l'altra la fascia corre con la regola fine (gli zeri da fuori sporcano al massimo 4
   celle della fascia in 4 tick, mai l'interno);
3. ogni punto dello spazio è contato **una volta sola** nell'energia.

**Una correzione dopo la pre-registrazione (dichiarata).** La prima versione (v1, quella pre-registrata) cuciva i
livelli **dopo ogni tick grosso**. Ma un tick solo (x poi y) non è simmetrico: la regola è precisa al second'ordine
solo dopo una **coppia** di tick (x,y poi y,x: un palindromo). Cucendo a metà coppia, i due livelli si scambiavano un
errore di ordine (k h)² che a 30° **creava o distruggeva circa l'1 % di energia** (a 16 celle/λ: +0,92 % passando da
grosso a fine, −0,95 % da fine a grosso; l'errore scendeva come N⁻², la firma di questo difetto). La v2 cuce **solo
dopo ogni coppia palindroma**: a 30° e 16 celle/λ l'errore scende a −0,12 % / +0,04 % e la riflessione da
2,8·10⁻⁵ a 6,7·10⁻⁶. A incidenza normale non cambia niente. Tutte le corse sotto sono della v2; quelle v1 sono in
`lab/out/reticolo_v1/`.

## 2. Tabellone

| prova | cosa | previsione | esito |
|---|---|---|---|
| R1 | reticolo rado a livelli, 2+1 D e 3+1 D, regola della luce su ogni livello | — | **fatto** (kernel = regola di `zz_core` entro 3·10⁻¹⁴ in 2-D e 3-D) |
| R2a asse | velocità entro 0,5 % su ogni livello (λ = 24 celle) | PASS | **PASS**: −0,162 % su livelli 0, 1, 10, 20, **identici al bit** (Maxwell, stesso pacchetto: −0,158 %) |
| R2a 30°/45° | idem in obliquo | FAIL (−0,65/−0,87 % dalla regola) | **FAIL** come previsto: −0,81 % / −1,02 % (servono ≥ 32 celle/λ) |
| R2a cucitura | velocità dopo l'interfaccia (asse) | PASS | **PASS**: 1,000000 c (entro 3·10⁻¹⁰) |
| R2b | interfaccia a ≥ 16 celle/λ (grosso): T ≥ 99 %, R ≤ 1 %, nei due versi | PASS | **PASS** (sotto) |
| R2c asse | fase ≤ 0,01 rad, fedeltà ≥ 0,999, elicità sbagliata + modo fermo ≤ 10⁻³ | PASS | **PASS**: fase ≤ 3·10⁻⁶ rad, fedeltà 1 − 10⁻¹⁰, elicità sbagliata 2·10⁻¹⁰ |
| R2c 30° | fedeltà ≥ 0,999, elicità | PASS | **FAIL a 16 celle/λ** (fedeltà 0,9979), PASS da 24 (0,99961): è la diversa velocità dei due livelli, non la cucitura (§4) |
| R2 3-D | 3-D, 16 celle/λ, asse | PASS | **PASS** (stessi numeri del 2-D: a incidenza normale il problema è 1-D) — obliquo 3-D in §3 |
| R3 | un mondo, gamma al livello 0 e radio al livello 20, celle ≤ 10⁶, memoria ≤ 100 MB | PASS | **PASS**: 717 568 celle, 34,7 MB (massimo 945 408 celle, 45,6 MB in corsa) |
| R3 gamma | gamma dopo i suoi tick: fedeltà ≥ 0,999 | PASS | **PASS**: 0,99973 dopo 2048 tick, fase 8·10⁻⁷ |
| R3 radio | radio dopo ≥ 1 tick radio (2²⁰ tick gamma) | PASS | **NON FATTO al livello 20** (servirebbero ~15 ore in Python, e il gamma intanto si allarga per diffrazione); fatto al livello 10: §5 |
| R4 | costo indipendente dalla scala; camino lineare nei livelli | PASS | **PASS**: 106 496 celle identiche ai livelli 0…40; +25 610 celle per livello di camino |
| R5 | formato per il Rust + questo documento | — | **fatto** |

## 3. R2: cosa succede quando un'onda attraversa il confine fra due livelli

Un pacchetto di luce (elicità +1, largo 3 λ) passa da un livello fine (h₀) a uno grosso (2h₀) e viceversa. N = celle
per λ **sul livello grosso** (sul fine sono 2N). T = energia che passa, R = energia che torna indietro.

**Incidenza normale** (lungo l'asse x):

| N | fine→grosso T | grosso→fine T | R (tutti) |
|---|---|---|---|
| 2 | 0,784 | 0,324 | ≤ 2·10⁻³ |
| 3 | 0,949 | 0,775 | ~0 |
| 4 | 0,983 | 0,849 | ~0 |
| 6 | 0,9965 | 0,963 | ~0 |
| 8 | 0,99889 | 0,98778 | ~0 |
| 12 | 0,99978 | 0,99750 | ~0 |
| **16** | **0,99993** | **0,99920** | **4·10⁻¹⁴** |
| 24 | 0,999986 | 0,99984 | 10⁻¹³ |
| 32 | 0,999996 | 0,99995 | 2·10⁻¹³ |

**A 30°:**

| N | fine→grosso T | R | grosso→fine T | R |
|---|---|---|---|---|
| 2 | 0,00004 | 0,03 | 0 | **1,00** (tutto riflesso) |
| 3 | 0,37 | 0,11 | 0,05 | 0,23 |
| 4 | 0,866 | 9·10⁻³ | 0,988 | 1,2·10⁻² |
| 8 | 0,983 | 1,5·10⁻⁴ | 1,006 | 1,6·10⁻⁴ |
| 12 | 0,9962 | 2,3·10⁻⁵ | 1,0013 | 2,4·10⁻⁵ |
| **16** | **0,9988** | **6,7·10⁻⁶** | **1,0004** | **6,8·10⁻⁶** |
| 24 | 0,99976 | 1,2·10⁻⁶ | 1,00008 | 1,2·10⁻⁶ |
| 32 | 0,99993 | 3,8·10⁻⁷ | 1,00003 | 3,8·10⁻⁷ |

Come si legge:
* **Da 16 celle per λ in su la cucitura è quasi invisibile**: passa il 99,9 % (normale) o il 99,88 % (30°), torna
  indietro meno di 10⁻⁵. La parte che manca va come N⁻⁴ (normale, è l'errore dell'interpolazione cubica) o N⁻²
  (30°, è la diversa dispersione dei due livelli, §4).
* **Le onde corte rispetto alla cella grossa NON passano, ed è giusto così.** Un livello grosso non può portare
  un'onda di 2–4 celle: a incidenza normale l'onda si **perde** (T = 0,32 a N = 2), di sbieco **torna indietro**
  (R = 1 a N = 2 da grosso a fine: l'onda resta intrappolata nel grosso). È la fisica di ogni reticolo a più
  livelli: per questo il livello giusto si sceglie con il criterio "≥ 16 celle per λ".
* **A incidenza normale la riflessione è zero** (10⁻¹³–10⁻³⁰): lungo un asse la regola è un puro spostamento e la
  parte che va avanti non si mescola con quella che va indietro. La riflessione nasce solo di sbieco.
* A 30° da grosso a fine T supera 1 di 4·10⁻⁴ (a 16 celle/λ): la cucitura non conserva l'energia **esattamente**;
  dopo la correzione v2 l'errore è piccolo (prima era ±0,9 %), ma non è zero. È un debito (§7).
* **Fase ed elicità** (R2c): a incidenza normale il pacchetto che esce è **lo stesso** di quello che non ha
  attraversato niente: fase spostata di ≤ 3·10⁻⁶ rad, elicità sbagliata 2·10⁻¹⁰, modo fermo 10⁻¹². A 30° l'elicità
  resta pulita (sbagliata 2·10⁻⁷, modo fermo ≤ 9,5·10⁻⁴ a 16 celle/λ).
* **3-D**: a incidenza normale (N = 8 e 16) i numeri sono identici al 2-D, come deve essere (il problema è 1-D).
  OBLIQUO_3D

## 4. Perché di sbieco la fedeltà è 0,998 a 16 celle/λ (R2c 30°: FAIL dichiarato)

Ogni livello porta la luce di sbieco **un po' più lenta di c** (è un difetto noto della regola, `ZIGZAG.md` §3.1:
a 30° la velocità di gruppo è −1,48 % a 16 celle/λ e −0,36 % a 32). Il livello grosso è più "lento" del fine:
il pacchetto che ha fatto metà strada sul grosso arriva un po' indietro (circa 5 celle fini su 500) e con la fase
girata di 0,37 rad rispetto a uno che ha viaggiato tutto sul fine. Il conto della dispersione dà 0,36 rad e una
perdita di fedeltà di ~9·10⁻⁴: **la cucitura non c'entra**, è la regola stessa ad avere velocità diverse a scale
diverse quando λ non è molto più grande della cella. Da 24 celle/λ la fedeltà passa (0,99961), a 32 è 0,99988.
**Lezione pratica: per la luce obliqua ogni livello va usato con ≥ 24–32 celle per λ**, non 16.

## 5. R3: un solo mondo con un "gamma" e una "radio"

**Il mondo.** 21 livelli. Al livello 0 un impulso "gamma" con λ = 24 h₀ (largo 1 λ lungo il volo e 3 λ di traverso,
va lungo +x). Al livello 20 un impulso "radio" con λ = 24 celle grosse = 24·2²⁰ h₀ ≈ 2,5·10⁷ h₀ (largo 1,5 λ, va a
60°). Il gamma sta dentro la radio, a 1,4 larghezze dal centro. In mezzo i livelli 1–19 esistono solo come **camino**
attorno al gamma. (Se il livello 20 fosse un'onda radio di 1 m, h₀ sarebbe 40 nm e il "gamma" avrebbe λ = 0,95 µm:
un infrarosso. Un vero gamma da 1 MeV contro una radio di 1 m sono 40 livelli, non 20: costa 20 piani di camino
in più, vedi R4.)

**Celle e memoria (all'inizio):**

| livello | 0 | 1 | 2 | 3 | 4 … 18 | 19 | 20 | totale |
|---|---|---|---|---|---|---|---|---|
| blocchi 16×16 | 337 | 196 | 140 | 120 | 100 ciascuno | 119 | 391 | 2803 |
| celle | 86 272 | 50 176 | 35 840 | 30 720 | 25 600 ciascuno | 30 464 | 100 096 | **717 568** |

Memoria: **34,7 MB** (numeri complessi a 64 bit; a 32 bit sarebbero 17 MB). Un reticolo **pieno** che tenesse la
radio intera alla risoluzione del gamma avrebbe (352·2²⁰)² ≈ 1,4·10¹⁷ celle: **2·10¹¹ volte di più**.

**La figura** (`lab/out/reticolo_zoom_ladder.png`, e un'immagine per livello `reticolo_zoom_L00…L20.png`; la radio
intera con i suoi blocchi in `reticolo_zoom_radio_whole.png`): ogni pannello è 2 volte più lontano del precedente.
Al livello 0 si vede il gamma (fronti verticali); salendo il gamma si rimpicciolisce fino a un puntino; dal livello 8
in su si vede un pezzetto della radio quasi uniforme, poi i suoi fronti inclinati, e al livello 20 la radio **intera**.
Si vede anche un difetto: ai livelli 2–7 resta un **fantasma** del gamma (ampiezza 0,4 al livello 4, 0,03 al livello
8), perché la copia fine→grosso è un'interpolazione che, per onde troppo corte, **ripiega** l'onda invece di farne la
media (§7).

**Il gamma corre** (2048 tick del livello 0, tutti i livelli insieme con lo schema di Berger–Oliger): confrontato con
la soluzione esatta della regola su un reticolo uniforme (calcolata in Fourier), **fedeltà 0,99973**, ampiezza
+0,026 %, fase 8·10⁻⁷ rad. Il camino segue il gamma (514 ricostruzioni dei blocchi in 2048 tick). Tempo: **52 ms per
tick del gamma** (Python + numba, 1 thread, CPU condivisa).

**La radio al livello 20 non l'abbiamo fatta muovere, e il motivo è il risultato più importante.** Un tick della radio
sono 2²⁰ ≈ un milione di tick del gamma: a 52 ms l'uno sono **15 ore** in Python (in Rust da qualche minuto a qualche ora,
secondo quanto si allarga il gamma). E in quel milione di tick il gamma, stretto 3 λ, si allarga per diffrazione fino a ~50 000 celle: un ventaglio
di 10⁷–10⁸ celle. La memoria regge, **il tempo no**. Per vedere davvero tutti e due muoversi abbiamo messo la radio al
livello 10 (1024 tick del gamma per tick radio):

R3B_TABLE

## 6. R4: quanto costa tenere le cose

* **Un impulso da solo** (λ = 24 celle del suo livello, largo 1,5 λ): **106 496 celle (416 blocchi) al livello 0, 1,
  5, 10, 20, 30, 40**, sempre le stesse. Il lato della regione va da 352 h₀ a 3,9·10¹⁴ h₀. **La scala assoluta non
  costa niente.**
* **Due oggetti a m livelli di distanza** (gamma + radio): il camino costa **25 610 celle per livello** (100 blocchi,
  il punto fisso della regola di annidamento con 2 blocchi di margine) + 206 106: m = 2, 5, 10 dal motore, m = 20 dal
  mondo di R3 (717 568 contro 718 306 della retta: 0,1 %).
* **Quindi**: costo ≈ (celle per tenere ogni oggetto al suo livello, ∝ (L/λ)^d) + (livelli fra le scale) ×
  (blocchi del camino). Il lato della scatola non entra.

**Stime oneste per la Via Lattea** (diametro ~10⁵ anni luce ≈ 10²¹ m; in 3-D, blocchi 8³, 48 byte per cella):

| cosa | celle | memoria | ci sta? |
|---|---|---|---|
| un impulso radio (λ = 1 m, h = 4 cm) ovunque nella galassia | 4,4·10⁷ | 2 GB | **sì** |
| la stessa radio + un gamma da 1 MeV (40 livelli di camino da 4·10⁶ celle) | 2,5·10⁸ | 12 GB (6 a 32 bit) | **sì**, su una macchina grossa |
| indirizzi: il gamma a risoluzione piena in qualunque punto della galassia | 114 bit per asse | — | sì con interi a 128 bit (universo osservabile: 134 bit, **non** basta i128 al livello 0 del gamma; basta salire di 7 livelli) |
| un'onda radio di 1 m che **riempie** la galassia | 1,2·10⁶⁷ | 6·10⁶⁸ byte | **no** (10⁶⁷: più degli atomi della Terra, ~10⁵⁰) |
| un'onda gamma che riempie la galassia | 6·10¹⁰² | — | **no**, mai |
| un fotone radio che **attraversa** la galassia | 2,3·10²² tick, e il fascio largo almeno √(Dλ/π) ≈ 1,7·10¹⁰ m = 4·10¹¹ celle | — | **no**: il tempo |
| un gamma che attraversa la galassia | 1,8·10³⁴ tick, fascio ≥ 19 km = 4·10¹⁷ celle | — | **no**, di 10²⁰ volte e più |
| tick del gamma per ogni tick della radio (1 MeV contro 1 m) | 2⁴⁰ ≈ 1,1·10¹² | — | è il muro del tempo |

**In una frase**: il reticolo rado a livelli può **contenere** la Via Lattea e **guardare** un'onda radio intera
accanto a un gamma (fotografia: milioni di celle). Non può **far vivere** a piena risoluzione cose piccole per tempi
lunghi su distanze grandi: lì serve un'altra idea (gli agglomerati che imparano la regola del livello, o un salto
analitico per la luce libera, che nel vuoto è nota esattamente).

## 7. Debiti e limiti, detti chiari

* **Correzione del motore dopo la pre-registrazione** (sincronizzazione a coppie di tick, v2), dichiarata in §1;
  le corse v1 sono conservate.
* **R3 al livello 20 non fatta muovere** (tempo e diffrazione, §5); sostituita, dopo, da R3b al livello 10. Anche
  questa è una deviazione dalla pre-registrazione.
* **Il fantasma del gamma** ai livelli 2–7: la restrizione fine→grosso è un'interpolazione (ottima per onde con
  ≥ 16 celle/λ), ma per contenuto troppo fine per il livello grosso ripiega l'onda invece di mediarla. Gonfia i
  blocchi dei livelli 2–4 (da 140→355 blocchi in 2048 tick). Va sostituita, per il contenuto non risolto, con una
  media conservativa — o meglio con un agglomerato che sa cosa rappresenta.
* **Energia non conservata esattamente all'interfaccia**: −0,12 % / +0,04 % a 30° e 16 celle/λ (v2). Una cucitura
  a flusso conservativo sarebbe la versione seria.
* **Velocità diverse dei livelli in obliquo** (§4): per la luce obliqua ≥ 24–32 celle per λ.
* **Sincronizzazione nel tempo**: i livelli grossi "aspettano avanti" (Berger–Oliger). L'energia totale del mondo ha
  senso solo quando tutti i livelli sono allo stesso istante (ogni 2^(K+1) tick fini).
* **Solo luce**. Le palline con massa hanno θ = m c h/ħ che raddoppia a ogni livello: oltre pochi livelli servono gli
  agglomerati con la regola imparata (`lab/agg_zz_*`), qui non collegati. Il formato ha il "gancio" per farlo.
* **3-D**: provato solo sull'interfaccia (R2 3-D), non il mondo intero con il camino.
* **Python lento**: 52 ms per tick del gamma con ~10⁶ celle in 21 livelli; il Rust deve fare 10–100 volte meglio.

**Numero di parametri**: **0 di fisica** (c e la regola sono quelli di `ZIGZAG.md`), **0 imparati**; di ingegneria:
blocchi 16 (8 in 3-D), fascia 4 celle, cucitura ogni 2 tick, interpolazioni cubiche, soglie 10⁻⁵ (dettaglio) e
10⁻⁴ (livello in cima), margine 1 blocco, camino con 2 blocchi, ricostruzione ogni 8 tick.

## 8. Il formato per il Rust

`assets/atomo/reticolo_spec.json`: livelli (h = 2ᵏ h₀, dt = h/c), blocchi (16² o 8³ celle, ordine x più veloce, 3
complessi per cella, blocco mancante = zero), la regola di un asse, lo schema dei tempi in pseudocodice, i pesi della
cucitura (prolungamento cubico ai quarti, restrizione −1, 9, 9, −1 su 16), la fascia, la regola dell'energia, la
ricostruzione dei blocchi, i limiti noti e tutti i 40 casi di R2 per controllare l'implementazione.

## Come rifarlo

    python lab/reticolo_prereg.py                 # pre-registrazione
    python lab/reticolo_r2.py 2 3 4 6 8 12 16     # interfaccia, incidenza normale (poi 24 32; --ang 30; --d3)
    python lab/reticolo_r2a.py 0 ; python lab/reticolo_r2a.py 30 45
    python lab/reticolo_r3.py ladder              # mondo a 21 livelli (4 min)
    python lab/reticolo_r3.py both                # radio al livello 10 + gamma (ripetere: riparte dal checkpoint)
    python lab/reticolo_r3.py frames              # immagini del mondo a t = 0
    python lab/reticolo_r4.py ; python lab/reticolo_figure.py ; python lab/reticolo_spec.py
