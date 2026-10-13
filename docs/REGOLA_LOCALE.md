# Regola locale pulita: doppioni e luce di traverso

Studio sui **debiti #3 e #4** di `CONTESTO_CLOUD.md` §4:

- **#3** la zig-zag "raddoppia lo spazio" (U(k+(π,π)) = U(k)), e l'unica cura (Wilson) esisteva solo in spazio k,
  cioè non era un automa locale vero;
- **#4** la luce di traverso si allarga col fattore ka/tan(ka) e in obliquo va piano: servono ≥ 24–32 celle per λ.

La domanda: esiste una versione **locale** di "corri a c, poi ruota" (stessa pallina, stessi ingressi c, ħ, m) che
tolga i doppioni e/o richieda meno celle per λ? Oppure si può dimostrare che non esiste?

**La regola congelata (`lab/zz_core.py`) non è stata toccata.** Qui ci sono solo proposte, con i test da rifare.

File: motore `lab/locale_core.py`, candidati `lab/locale_cands.py`, pre-registrazione `lab/locale_prereg.py` →
`lab/out/locale_prereg.json` (commit **9d4f0bc**, prima delle corse), corse `lab/locale_r1.py`, `locale_r2.py`,
`locale_r3.py`, emendamento dichiarato `lab/locale_amend.py`, esplorazione dopo la pre-registrazione
`lab/locale_hexs4.py`, figura `lab/locale_figure.py` → `lab/out/locale_figure.png`.

---

## 0. In breve (per Valerio)

1. **I "due fogli" erano in realtà 8 (in 3-D) o 4 (in 2-D) universi identici e intrecciati.** Ogni pallina in un tick
   si sposta di ±1 su *ogni* asse, quindi dopo due tick tutte le coordinate sono cambiate di un numero pari. Il reticolo
   cubico si divide in 8 sotto-reticoli (2Z³ spostato) che non si parlano mai. La regola congelata calcola **8 volte
   la stessa fisica**. Basta tenerne uno: **stessa identica fisica con 1/8 delle celle** (verificato: differenza
   3·10⁻⁷ celle sul centro e sulla larghezza di un pacchetto, A2). È locale ed esatto, niente spazio k.
2. **Dentro ognuno di questi universi restano però dei doppioni veri** (fermioni "fantasma" a bassa energia nei punti
   N della zona di Brillouin, con chiralità opposta). **Teorema (alla Nielsen–Ninomiya, per i passeggi):** per *ogni*
   regola fatta a passi lungo una direzione alla volta, il "numero di avvolgimento" W₃ è zero, quindi la somma delle
   chiralità a energia 0 è zero: un elettrone senza massa a k = 0 si porta dietro per forza dei gemelli. Verificato su
   14 regole (tutte 0,0000) e con un controllo positivo che il calcolo sa vedere 1 e −2.
3. **Per l'elettrone con massa il teorema non vieta niente**, e ho trovato la **cura locale**: una piccola rotazione
   "sui collegamenti a due passi" (la parte sinistrorsa in x ruota con la destrorsa in x ± 2), attorno a un asse
   **perpendicolare** a quello della massa ("Wilson ritorto"). Risultato: i doppioni salgono da E = m a E ≥ 0,39/a
   (anche quando m → 0), la massa a riposo resta **esatta** (4·10⁻¹⁶) e la massa cinetica è **identica** a quella della
   regola congelata (differenza 4·10⁻¹⁰). Costo: raggio 2 celle, ~10 volte più operazioni per tick, ma su 1/8 dei siti.
4. **Luce: con soli passi in avanti non si può fare meglio** (teorema di Sheng–Suzuki: niente ordine > 2; ricerca su
   14 parole: la migliore è proprio la zig-zag). **Con un passo all'indietro sì**: la sequenza di blocchi
   **avanti 4, indietro 2, avanti 4** (S4: 1,1,1,1,−2,1,1,1,1; somma dei cubi = 0) è di **quarto ordine**. Bastano
   **14,6 celle/λ** per l'1% (contro 36,3), il costo in 3-D scende a **×0,10**. Sull'elettrone la massa cinetica passa
   da 3,3·10⁻³ a **1,8·10⁻⁵** (θ = 0,1). **Ma**: pagando con onde non risolte più veloci della luce (fino a 1,66 c) e un
   leggerissimo sorpasso di c (1,001 c) anche per onde risolte. Per questo l'esame pre-registrato M3 è **FAIL**.
5. Il reticolo **esagonale** in 2-D (luce lungo 3 direzioni a 60°) dà 21,7 celle/λ (costo ×0,41) **senza** nessun
   sorpasso di c: è la cura "pulita" per il 2-D, solo positiva. In 3-D nessun reticolo (fcc, bcc, misto) aiuta.

Tabellone: r1 **66 PASS / 1 FAIL**, r2 **33 PASS / 16 FAIL**, r3 **20 PASS / 12 FAIL**, emendamento 6/0.
Molti FAIL sono **previsti nella pre-registrazione** (la regola congelata che fallisce gli esami dei doppioni); quelli
inattesi sono spiegati in §6.

![figura](../lab/out/locale_figure.png)

*(a) energia più bassa lungo Γ–N–H–Γ–P (θ = 0,1): la regola congelata ha una valle finta a E = m nel punto N (un
doppione); con il Wilson locale la valle sale a 0,8. (b) errore sulla diffrazione contro celle per λ. (c) errore della
massa cinetica: S4 scende come θ⁴. (d) celle per λ per l'1% e costo relativo.*

---

## 1. Come ho lavorato

Ogni candidato è un **prodotto di pezzi locali**: spostamenti (la parte di autovalore ±1 di s·n si sposta di ±v),
rotazioni sul posto, rotazioni su un collegamento. Per una regola uguale dappertutto, l'effetto su un'onda piana è una
matrice W(k) (la "matrice di Bloch"): da lì si leggono energie, velocità, doppioni, avvolgimenti. Ho controllato che la
matrice coincida con la regola fatta girare sulle palline (G4, ≤ 2·10⁻¹⁴ per tutti) e, per la congelata, con
`zz_core.bloch_batch` (4·10⁻¹⁶).

Esplorazioni di sviluppo (dichiarate nella pre-registrazione, voce `dev_seen_before_prereg`): ho visto i numeri della
luce congelata, dell'esagonale, di S4 e del Wilson prima di scrivere le soglie; le soglie non sono state aggiustate su
quei numeri.

---

## 2. I fogli: 8 universi invece di 1 (debito #3, parte facile)

**Perché.** In un tick ogni componente si sposta di (±1, ±1, ±1). Dopo due tick lo spostamento ha tutte le coordinate
pari. La rotazione (massa) e il potenziale agiscono sul posto. Quindi un sito con coordinate tutte pari ai tick pari
parla solo con siti "tutti dispari" ai tick dispari, e viceversa: **8 classi** in 3-D (4 in 2-D) che non si toccano.
In formula: W(k + π e_j) = W(k) (T1a: 7·10⁻¹⁶). Lo stesso ragionamento vale per qualsiasi passeggio di spin ½ con
rotazioni sul posto: il numero di copie è almeno 2^d.

**Prova.** Pacchetto 3-D (larghezza 3 celle, k = 0,3) fatto girare 32 tick sul reticolo intero e su una sola classe:

| esame | risultato | soglia | esito |
|---|---|---|---|
| T1b sulla classe resta zero fuori | 0 esatto | 0 | PASS |
| T1b stessa fisica, scatola 32³ | centro diverso di 1,18 celle | 1% della larghezza | **FAIL** (pacchetto più largo della scatola: il "centro" in una scatola periodica non ha senso) |
| A2 (emendamento) stessa fisica, scatola 96³ | 3,4·10⁻⁷ / 1,3·10⁻⁷ celle (larghezza 27,1) | 0,27 | PASS |

Il FAIL di T1b resta nel tabellone: è un errore mio nella scelta della scatola, dichiarato e rifatto in A2.

**Conseguenza pratica (grossa).** Tutti i codici che usano la zig-zag per l'elettrone in 3-D (Laboratorio Atomo,
`zz_atom3d`, accoppiato) possono girare su **1/8 dei siti** con la stessa precisione. L'idrogeno 3-D "da ~3·10¹¹ celle"
del debito #6 diventa ~4·10¹⁰.

---

## 3. I doppioni veri e il teorema (debito #3, parte difficile)

Dentro una sola classe la regola vive su un reticolo cubico di passo 2 con zona di Brillouin [−π/2, π/2)³. Lì, per
l'elettrone senza massa (Weyl), i punti dove la regola ha un cono sono (r1, T2, `find_nodes`):

| energia | punti (k in unità di π) | chiralità | somma |
|---|---|---|---|
| 0 | Γ (0,0,0); spigoli N: (½,½,0), (½,0,½), (0,½,½) | −1; +1, −1, +1 | **0** |
| π (per periodo) | facce (½,0,0), (0,½,0), (0,0,½); vertice (½,½,½) | +1, −1, +1; −1 | **0** |

**Una regola geometrica bella:** in coordinate q = 2k sono gli 8 vertici del cubo {0, π}³. Quelli con un numero
**pari** di π stanno a energia 0, quelli con un numero **dispari** stanno a energia π. È il raddoppio di
Nielsen–Ninomiya (2³ = 8 doppioni), con metà spostati a energia π dal tempo discreto.

**Teorema (versione per i passeggi).** Per una matrice 2×2 di determinante 1, la somma delle chiralità dei coni a energia
0 è uguale al grado W₃ della mappa k → W(k) sulla sfera S³ (e così quella a energia π). W₃ di un prodotto è la somma
dei W₃ dei fattori (la classe di grado 3 di un gruppo di Lie è "primitiva"). Ogni fattore della zig-zag (spostamento
lungo v, rotazione su un collegamento v, rotazione sul posto) dipende da **una sola** combinazione k·v: passa per un
cerchio, e una mappa che passa per un cerchio ha W₃ = 0. Quindi **ogni regola a passi lungo una direzione alla volta
ha W₃ = 0**: un cono a Γ ha per forza dei gemelli di chiralità opposta a energia 0. Inoltre W(k) = e^{−iH(k)T} con H
continua ha anch'esso W₃ = 0 (si deforma fino all'identità).

Verifiche: W₃ = 0,0000 per U(xyz), U(zyx), W congelato, S4 e 10 regole casuali (con direzioni, assi di spin e rotazioni
su collegamenti casuali): 14/14 PASS. Controllo positivo (emendamento A1): su una mappa non locale nota il calcolo dà
+1,000000, −2,000000, +1,000000, 0 come atteso.

**Onestà sui limiti del teorema.** È dimostrato per i prodotti di pezzi a una direzione (tutto quello che si può
costruire "corri lungo un asse, ruota"). Per un automa locale generico in 3-D che non si scompone così non l'ho
dimostrato: congettura, collegata al teorema di Dubail–Read (proiettori a raggio finito hanno numero di Chern zero).
Per un neutrino (Weyl senza massa) quindi **i doppioni non si tolgono** con nessuna regola di questo tipo.

### 3.1 La cura locale per l'elettrone con massa: il Wilson "ritorto"

Con 4 componenti (L e R) il teorema non vale: basta dare ai doppioni una **massa grande**, che sia zero a Γ.

**La regola proposta (F2).** Dopo gli spostamenti di ogni tick:

1. metà del termine di Wilson: per ogni v in {±2e_x, ±2e_y, ±2e_z} la parte L in x e la parte R in x+v ruotano
   l'una nell'altra di un angolo piccolo (h/2 diviso sui 6 collegamenti, ordine palindromo), attorno all'asse **Y** di
   chiralità (perpendicolare all'asse X della massa); poi si toglie la parte "sul posto", così a k = 0 non succede
   niente;
2. la rotazione della massa θ come oggi;
3. l'altra metà del Wilson;
4. nei tick dispari l'ordine dei collegamenti è specchiato (−v prima di +v).

In spazio k al primo ordine il termine è −h Σ_v (1 − cos k·v) Y: zero a Γ, grande nei doppioni. Raggio 2 celle
(resta sulla stessa classe del §2).

Perché "ritorto": se l'asse fosse X (lo stesso della massa) la massa efficiente diventerebbe θ + w(k) e la massa
cinetica sbaglierebbe di O(h) (F2b: +9,1% a h = 0,1). Con Y si ha E² = p² + θ² + w(k)²: nessun termine misto, e w² ~
(h p²)² è piccolissimo. Il trucco dello specchio nei tick dispari è servito perché con l'ordine fisso la particella a
riposo **scivolava** a 2,4·10⁻³ c (∝ h³, dall'ordine dei collegamenti); con lo specchio la deriva è 3,5·10⁻¹² c.
Nella fisica dei reticoli è il "twisted mass" di Frezzotti–Rossi: là lo usano per lo stesso motivo.

| esame (θ = 0,1, 3-D) | congelata F0 | **Wilson Y, F2** | Wilson X, F2b | soglia |
|---|---|---|---|---|
| D1 E(0)/m − 1 | 4·10⁻¹⁶ | 4·10⁻¹⁶ | 4·10⁻¹⁶ | 10⁻¹² |
| D2 dispersione p ≤ m | 1,94·10⁻³ | 1,55·10⁻³ | **2,2·10⁻² FAIL** | 1% |
| D3 massa cinetica | +3,35·10⁻³ | +3,35·10⁻³ | **+9,1·10⁻² FAIL** | 1% |
| D4 deriva a riposo | 4·10⁻¹² c | 3,5·10⁻¹² c | 2·10⁻¹² c | 10⁻⁹ |
| D5 doppioni, θ = 0,1 | **0,100 FAIL** (punto N) | 0,404 PASS | 0,389 PASS | ≥ 0,3 |
| D5 doppioni, θ = 0,02 | **0,020 FAIL** | 0,392 PASS | 0,389 PASS | ≥ 0,3 |
| D7 velocità massima nella zona | 0,995 c PASS | **1,030 c FAIL** | 1,028 c FAIL | ≤ c |

In 2-D: F0 D5 FAIL (0,100 e 0,020 nel punto (½,½)π), F2 PASS (0,405 e 0,394), D7 F2 2-D 1,031 c FAIL.

Il valore 0,40 di D5 è il cono di Γ stesso al bordo |k| = 0,4: **sotto 0,39/a non c'è nient'altro**, e questo non dipende
da m (i doppioni non tornano giù quando la massa va a zero).

D6 (indipendenza da h): la massa cinetica non cambia con h (differenze 1,6·10⁻¹⁰, 4·10⁻¹⁰, 8·10⁻⁸ a h = 0,05; 0,1;
0,2) ma la dispersione a p = m cambia di 1·10⁻⁴ … 1,2·10⁻³ per il termine (h p²)². La soglia pre-registrata (10⁻⁶ su
entrambe) era troppo stretta: **FAIL** ×3. Con h = 0,1 l'effetto (4·10⁻⁴) è più piccolo dell'errore che la regola
congelata ha già (1,9·10⁻³), e va a zero come a⁴.

D7: il Wilson allarga il "cono di luce" del reticolo (collegamenti di raggio 2): onde corte con |k| ≥ 0,62 (3-D) / 0,43
(2-D) vanno fino a 1,03 c. Sono onde non risolte (meno di ~10 celle per lunghezza d'onda), ma è un difetto: dichiarato.

**Costo.** Per tick: 3 spostamenti + 24 rotazioni di collegamento + 3 rotazioni sul posto (contro 3 + 1). Su 1/8 dei
siti: in totale circa **lo stesso costo** della regola congelata di oggi, ma senza doppioni.

**Parametri.** Uno solo nuovo, h (regolatore senza fisica dentro: la massa cinetica non lo vede). Gli ingressi fisici
restano c, ħ, m.

---

## 4. La luce: meno celle per λ (debito #4)

### 4.1 Perché la zig-zag sbaglia di traverso, e perché non si può aggiustare con soli passi in avanti

Lungo un asse lo spostamento di una cella **è** la propagazione esatta per un tempo a/c. Quindi l'unico errore è quello
di "fare un asse alla volta" (errore di Trotter), con un passo di tempo grande (a/c). La zig-zag alterna l'ordine
(xyz, zyx): è la formula di Strang, errore ∝ (ka)². Teorema di Sheng (1989) e Suzuki (1991): **qualsiasi** formula
di questo tipo di ordine > 2 ha bisogno di passi negativi. Con passi interi in avanti, la ricerca su tutte le parole
palindrome fino a 12 lettere (2-D) dà come migliore proprio "xyyx", cioè la zig-zag (T3: 36,33 contro 36,33, PASS).

Celle per λ per l'1% (r3, M3; 2-D: 46 direzioni, 3-D: 160):

| candidato | fase | gruppo | diffrazione | peggiore | siti/λᵈ | costo ×congelata | onde risolte > c? |
|---|---|---|---|---|---|---|---|
| L0 congelata 3-D | 14,8 | 25,7 | **36,3** | 36,3 | 4,8·10⁴ | 1 | no |
| L3 ordini ciclati 3-D | 15,8 | 26,8 | 37,1 | 37,1 | 5,1·10⁴ | 1,06 | no |
| L4a fcc (diagonali di faccia) | 29,5 | 50,8 | 75,6 | 75,6 | 2,2·10⁵ | 3,19 | no |
| L4b bcc (diagonali di cubo) | 33,8 | 58,4 | 81,4 | 81,4 | 1,4·10⁵ | 1,63 | no |
| L5 assi + bcc (anisotropia cancellata al primo ordine) | 19,9 | 32,9 | 37,6 | 37,6 | 5,3·10⁴ | 0,98 | no |
| **L2 S4 3-D (avanti 4, indietro 2, avanti 4)** | 12,5 | 14,5 | **14,6** | **14,6** | 3,1·10³ | **0,10** | **1,001 c** |
| L0 congelata 2-D | 13,0 | 22,3 | 36,3 | 36,3 | 1320 | 1 | no |
| **L1 esagonale 2-D (6 ordini ciclati)** | 13,0 | 20,9 | 21,7 | 21,7 | 545 | **0,41** | no |
| **L2 S4 2-D** | 12,4 | 14,3 | 14,4 | 14,4 | 207 | **0,23** | **1,001 c** |

Esame M3 pre-registrato ("debito ridotto" se peggiore ≤ 20 **e** costo ≤ 1 **e** nessuna onda risolta più veloce di c):
**nessun candidato passa.** S4 manca solo per il terzo criterio (1,001 c), l'esagonale solo per il primo (21,7 > 20).

### 4.2 S4: il passo del gambero

La composizione palindroma di blocchi a due tick con passi interi m = (1,1,1,1,−2,1,1,1,1) ha Σm = 6 e Σm³ = 0: le
condizioni di Yoshida/Suzuki per il quarto ordine, con numeri interi. (Ricerca esaustiva: non esiste per lunghezze ≤ 7
con |m| ≤ 7; a lunghezza 9 questa è quella coi passi più piccoli.) Il blocco con m = −2 è la regola fatta **all'indietro
nel tempo**: le palline si spostano di 2 celle contro il loro spin e il flip torna indietro di 2θ. Resta "corri a c,
poi ruota", con la stessa pallina e gli stessi ingressi.

- Luce: c esatta lungo gli assi come prima; |ω/ck − 1| a 1000 celle/λ = 8·10⁻¹¹ (congelata 2·10⁻⁶).
- Elettrone (F3): massa a riposo esatta (D1 PASS), dispersione 2,5·10⁻⁵ (D2), massa cinetica −1,8·10⁻⁵ (D3), ordine
  in θ: rapporto 14,6 tra θ = 0,1 e 0,05 (D8 PASS, atteso 16). Per gli atomi è un guadagno enorme: l'errore θ²/3 del
  reticolo diventa ~θ⁴.
- **Prezzi**: (i) sotto ~10 celle/λ si rompe di colpo (a 8 celle: errore 200%, contro il 21% della congelata);
  (ii) il cono di luce del reticolo si allarga a 20 celle in 12 tempi = **1,66 c**: il 4% delle onde della zona va più
  veloce della luce (D7 FAIL 1,64 c per l'elettrone); (iii) anche onde risolte superano c di 10⁻³ (M3 FAIL);
  (iv) l'ordine misurato tra 16 e 32 celle è 5,7–6,0 invece di 12–20 (T3 FAIL: a 16 celle il quarto ordine non è
  ancora puro; tra 32 e 48 il rapporto è ~6, cioè (3/2)⁴ = 5,1 — coerente col quarto ordine, ma l'esame era su 16→32);
  (v) il periodo è 12 tick: le energie si "ripiegano" ogni π/12a, quindi serve E·a < 0,26 (θ ≤ ~0,1);
  (vi) il raggio di un pezzo resta 2 celle, ma in un periodo una componente può andare avanti/indietro fino a 20 celle.

### 4.3 Esagonale (2-D)

Con tre direzioni a 60° la simmetria è di ordine 6, e un errore "di grado 4" con simmetria 6 è per forza uguale in tutte
le direzioni: la luce sbaglia allo stesso modo in ogni direzione, e l'errore di diffrazione scende da 36 a 21,7 celle/λ.
Il costo per area è ×0,41 e nessuna onda supera c. Quello che non va: il modo "fermo" (longitudinale) è un po' più
mescolato (difetto 1,0·10⁻² a 16 celle contro 4·10⁻⁴), e serve un reticolo diverso da quello degli altri codici.

### 4.4 Cose provate che non funzionano (3-D)

Ciclare gli ordini (L3): nessun guadagno, l'errore dipende dal passo, non dall'ordine. fcc/bcc: peggio (passi più
lunghi). Mescolare assi e bcc per cancellare l'anisotropia al primo ordine (L5): l'anisotropia va via ma resta un errore
isotropo grande: 37,6. Esagonale + S4 (`locale_hexs4.py`, esplorazione dopo la pre-registrazione): 18,9 celle/λ,
peggio di S4 sul quadrato (14,4).

---

## 5. Tabellone

| corsa | PASS | FAIL | FAIL previsti nella pre-registrazione | FAIL non previsti |
|---|---|---|---|---|
| r1 (`locale_r1.json`, 65 s) | 66 | 1 | — | T1b (scatola troppo piccola; rifatto in A2: PASS) |
| r2 (`locale_r2.json`, 189 s) | 33 | 16 | D5 F0 ×4, D2/D3 F2b, D6 F2b ×3, D7 F2/F3 ×3 | D6 F2 ×3 (soglia troppo stretta sulla dispersione); D7 F2b |
| r3 (`locale_r3.json`, 158 s) | 20 | 12 | M3 L0 ×2 (la congelata) | M3 di tutti gli altri; M1 L4b (1,1·10⁻⁵ contro 10⁻⁵); T3 S4 ×2 (ordine 5,7–6,0) |
| emendamento A1–A2 (`locale_amend.json`) | 6 | 0 | | |

---

## 6. Fallimenti e debiti

**Fallimenti dichiarati.**

1. Nessuna regola locale toglie i doppioni del **neutrino** (Weyl senza massa): è il teorema del §3. Per l'elettrone
   sì, con il Wilson ritorto.
2. Nessun candidato passa l'esame M3 come scritto: S4 per 10⁻³ di velocità sopra c; l'esagonale per 1,7 celle.
3. Il Wilson e S4 allargano il cono di luce del reticolo (1,03 c e 1,66 c per onde non risolte).
4. Il mio T1b era scritto male (scatola), il mio D6 aveva una soglia irrealistica per la dispersione, il mio T3 per S4
   guardava un intervallo dove il quarto ordine non è ancora puro.
5. I numeri dei doppioni a energia π (il "buco π", rilevante per il conto delle coppie di `IDEE_RIPRESE.md` R1b) non
   migliorano: distanza minima da π 0,075 (F0) e 0,078 rad (F2). Il Wilson non ripara il buco π.

**Debiti chiusi o ridotti.**

- **#3 per l'elettrone**: chiuso in proposta (F2 locale, raggio 2), con i test da rifare (§7). Per i fogli: chiuso
  (1 classe su 8, esatto).
- **#3 per il neutrino / Weyl senza massa**: non chiudibile in questa classe di regole (teorema). Resta come proprietà.
- **#4**: ridotto in proposta (S4: 14,6 celle/λ, costo ×0,10 in 3-D; esagonale in 2-D 21,7, costo ×0,41), con il
  prezzo dichiarato del cono di luce. Non chiuso secondo l'esame pre-registrato.

**Debiti nuovi.**

- D-L1: h del Wilson è una manopola (regolatore). La massa cinetica non la vede, la dispersione sì a O((hp²)²).
- D-L2: S4 ha onde corte superluminali (fino a 1,66 c) e un sorpasso di 10⁻³ per onde risolte; con un accoppiamento non
  lineare o una sorgente stretta potrebbero essere eccitate. Va misurato in un esperimento vero.
- D-L3: il teorema dei doppioni è dimostrato per i prodotti di pezzi a una direzione, non per ogni automa locale.
- D-L4: la potenziale sulla singola classe di siti cambia leggermente la media rispetto a oggi (prima era la media di 8
  campionamenti spostati di una cella).

**Parametri liberi.** Fogli: 0. Wilson: 1 regolatore (h = 0,1). S4: 0 (una sequenza intera, non un numero reale).
Esagonale: 0.

---

## 7. Proposta (non applicata) e test da rifare

Proposta in tre pezzi indipendenti, da adottare solo se Valerio lo decide:

1. **Una sola classe di siti** per l'elettrone (subito, nessun rischio): rifare `zz_free` (F1–F5), `zz_atom3d` (A1–A2) e
   l'accoppiato con il potenziale campionato sulla classe; atteso: stessi numeri a meno di ~10⁻⁶, tempo ÷8.
2. **Wilson ritorto (F2)** dove contano i modi a bassa energia della zona intera: vuoto e coppie (`IDEE_RIPRESE` R1,
   `COPPIE_IN_SCATOLA`), lettore della forma del vuoto (`FORMA_DALLE_PALLINE`, i "due fogli"), agglomeratore (che oggi
   può imparare i doppioni). Da rifare: Schwinger/Sauter (attese: stesse probabilità per modo a k = 0), conteggio dei
   doppioni nel vuoto, 1s di `zz_atom3d`.
3. **S4** per la luce dove servono fronti d'onda belli con poche celle (doppia fenditura, fotone nel Laboratorio Atomo)
   e per gli atomi dove serve precisione (massa cinetica θ⁴). Da rifare: L1–L5 di `zz_light`, `zz_slit` a 8 e 16
   celle/λ (attesa: a 16 celle la correlazione deve superare 0,995; a 8 celle **peggiora**), il fronte di un impulso
   stretto (cercare precursori superluminali), `zz_free` F1–F3, struttura fine con la forma radiale (S4 radiale da
   scrivere).

### Patch proposte per altri file

- `lab/zz_core.py`: nessuna modifica (congelato). Al massimo una classe nuova `WalkS4` / opzione `coset=True` in un file
  nuovo che importa `zz_core`.
- `src/zigzag/walk.rs` e `src/atomo/**` (altro thread): indicizzazione su una sola classe 2Z³ (all'indice (i,j,k)
  della classe corrisponde il sito (2i+p, 2j+p, 2k+p), p = parità del tick); opzione S4 come sequenza di tick con
  passo m e flip mθ.
- `docs/CONTESTO_CLOUD.md` §4: #3 → "fogli: 8 copie esatte, basta una (REGOLA_LOCALE §2); doppioni: teorema W₃ = 0 per
  Weyl; cura locale per Dirac: Wilson ritorto (§3.1)". #4 → "S4 (passo indietro) 14,6 celle/λ con cono 1,66 c; esagonale
  2-D 21,7".
- Richiesta di ponte proposta (NON creata): `ponte/richieste/locale_s4_foto.json` con `{"azione": "foto", "scena":
  "fotone libero", "regola": "S4", "celle_per_lambda": 12}` per vedere se a schermo la differenza si nota (serve prima il
  porto in Rust).

---

## 8. Idee geometriche e visive

Ognuna: cosa si vedrebbe, perché funzionerebbe, costo, regole.

1. **Gli 8 universi a colori.** A schermo: colorare le palline secondo la loro classe (8 colori, come una scacchiera 3-D);
   un pacchetto liscio appare come 8 copie sovrapposte che si muovono insieme; spegnendone 7, la forma non cambia.
   Perché: §2, esatto. Costo: nessuno (anzi ÷8). Regole: rispettate, nessun fake. È anche un bel modo per spiegare cosa
   vuol dire "doppione".
2. **Il cubo dei doppioni.** A schermo: la zona di Brillouin come un cubo trasparente; al centro e su 3 spigoli
   brillano 4 coni a energia 0 (i gemelli, colori per chiralità ±1), sulle 3 facce e su un vertice i 4 coni a energia π.
   La regola "numero pari di π → energia 0, dispari → π" si vede. Con il Wilson acceso i coni sugli spigoli si
   spengono (la valle si alza, figura (a)). Perché: r1 T2, tabella §3. Costo: solo grafica di un diagramma. Regole: è
   un'immagine della regola, non della natura.
3. **Il "fiore" dell'anisotropia.** A schermo: un anello di luce che parte da un punto. Con la zig-zag il fronte è un
   po' quadrato (più lento in diagonale: 4 petali); sull'esagonale è rotondo (simmetria 6 → errore uguale in tutte le
   direzioni); con S4 è rotondo e più preciso. Perché: §4.3, tabella M3. Costo: esagonale ×0,41 in 2-D. Regole: ok;
   l'esagonale è un reticolo diverso, non un nuovo ingrediente.
4. **Il passo del gambero (S4).** A schermo: ogni 18 tick le palline fanno 8 tick in avanti, 2 all'indietro (le scie
   dei surfisti tornano indietro un attimo), 8 in avanti; il fronte d'onda alla fine è più pulito. È visivamente strano
   ma onesto: la regola "ripassa" all'indietro per cancellare l'errore. Perché: Σm³ = 0. Costo ÷10 in 3-D a parità di
   errore. Regole: stessa pallina, stessa regola fatta anche all'indietro; il cono 1,66 c va mostrato, non nascosto.
5. **Collegamenti che respirano (Wilson ritorto).** A schermo: tra palline a due passi di distanza un piccolo
   scambio di colore L↔R, forte solo dove il campo cambia segno da una cella all'altra (onde corte), spento dove il campo
   è liscio. È il Wilson che "spegne i gemelli". Perché: §3.1. Costo: ~10× operazioni per tick, compensato dal ÷8.
   Regole: è un pezzo locale nuovo con una manopola (h), dichiarata.
6. **Esplorazione fatta (la più promettente tra le geometriche nel mio studio): esagonale + S4.** Risultato: 18,9
   celle/λ, **peggio** di S4 sul quadrato (14,4) e costo 1240 contro 620 (2-D). Quindi: esagonale per la luce senza passi
   indietro, quadrato per S4. Non si sommano.

---

## 9. Dove l'intuizione di Valerio può aiutare

1. **Il passo all'indietro è fisica o trucco?** Una regola che ogni tanto va indietro nel tempo per correggersi:
   ti sembra accettabile per "una sola strada"? Se no, il debito #4 resta a ~36 celle/λ in 3-D (teorema di
   Sheng–Suzuki), e l'unico guadagno pulito è l'esagonale in 2-D.
2. **Il cono di luce del reticolo.** La congelata non supera mai c. S4 sì, solo per onde corte non risolte (fino a
   1,66 c). Per te conta più "mai più veloce della luce, nemmeno per errore" o "giusto con meno celle"?
3. **I gemelli del neutrino.** Il teorema dice che una particella senza massa e di una sola mano porta sempre gemelli.
   In natura il neutrino ha una massa piccolissima: la regola "suggerisce" che ogni particella abbia almeno una massa
   minima o un partner? È un'osservazione sul reticolo, non una previsione: si può trasformare in un test?
4. **Gli 8 universi intrecciati.** Fisicamente sono 8 copie che non si parlano. Ti piacerebbe usarle (ad esempio 8
   "sapori" o 8 elettroni indipendenti nello stesso reticolo, come nella fogliazione di `MULTIFOTONE`), oppure buttarle?

---

## 10. Come rifare tutto

```
OMP_NUM_THREADS=1 python lab/locale_prereg.py      # pre-registrazione (già committata)
OMP_NUM_THREADS=1 timeout 590 python lab/locale_r1.py     # G1-G4, T1, T2 (65 s)
OMP_NUM_THREADS=1 timeout 590 python lab/locale_r2.py     # fermioni D1-D8 (189 s)
OMP_NUM_THREADS=1 timeout 590 python lab/locale_amend.py  # emendamento A1-A2 (55 s)
OMP_NUM_THREADS=1 timeout 590 python lab/locale_r3.py     # luce M1-M4, T3 (158 s; exit 3 = rilanciare)
OMP_NUM_THREADS=1 python lab/locale_hexs4.py       # esplorazione esagonale + S4
OMP_NUM_THREADS=1 python lab/locale_figure.py
```

Letteratura citata (da verificare sulle fonti prima di usarla altrove): Nielsen–Ninomiya (1981); D'Ariano–Perinotti,
PRA 90, 062106 (2014) (la regola a un tick è il loro automa di Weyl su bcc: cos ω = c_x c_y c_z + s_x s_y s_z);
Sheng (1989), Suzuki (1991) per l'ordine dei metodi a passi; Yoshida (1990); Frezzotti–Rossi (twisted mass, 2004);
Dubail–Read (2015) per i proiettori a raggio finito.
