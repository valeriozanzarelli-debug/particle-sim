# La vibrazione di H₂: perché veniva −1,4 %

L'unica anomalia che non capivamo (`docs/MOLECOLE.md` §7.2): dalle palline (DMC) la molecola H₂ aveva la lunghezza
di legame giusta (0,7420 Å) e l'energia di legame giusta (D₀ 4,472 eV), ma la vibrazione bassa:
ωₑ = 4339 ± 18 cm⁻¹ contro 4401,21 (−1,4 %, 3,5 σ) e ν(0→1) = 4129 ± 16 contro 4161,17 (−0,8 %).

**Causa trovata: il fit della curva, non le palline.** La formula usata per interpolare i 18 punti E(R)
(5 esponenziali, β su una griglia a passi di 0,04, E∞ fisso) non riesce a seguire bene il fondo della buca. Se le
diamo i valori **esatti** dell'energia (Pachucki 2010) alle stesse 18 distanze, senza nessun rumore, restituisce
ωₑ = 4333,6 invece di 4404,5: **−1,6 %**, e ν(0→1) **−0,86 %**. È esattamente la nostra anomalia. I punti delle
palline sono buoni; era il righello sbagliato.

Corretto (fit a 6 esponenziali, scelto **solo** sui dati esatti, poi applicato alle palline):
**ωₑ = 4393 ± 36 cm⁻¹ (−0,18 %, 0,2 σ)** e **ν(0→1) = 4152 ± 20 cm⁻¹ (−0,21 %, 0,4 σ)**. Anomalia risolta;
resto −0,2 %, dentro l'errore. L'errore vero però è **il doppio** di quello dichiarato prima (vedi B4).

Codice: `lab/h2vib_run.py` (lista B0–B7), `lab/h2vib_fix.py` (numero corretto), `lab/h2vib_pachucki.py` (tabella
esatta, usata solo come riferimento), figura `lab/h2vib_figure.py` → `lab/out/h2vib_figure.png`.
Pre-registrazione: `lab/out/h2vib_prereg.json` (commit 64d0afd, prima delle corse). Risultati:
`lab/out/h2vib_run.json`, `lab/out/h2vib_fix.json`. Nessuna nuova corsa di palline: uso i 33 valori già salvati in
`assets/molecole/risultati.txt`.

---

## 1. Il riferimento

K. Pachucki, *Born–Oppenheimer potential for H₂*, Phys. Rev. A 82, 032509 (2010), tabella IV: E(R) esatta a
10⁻¹⁵ per R da 0,1 a 20 bohr. Interpolata con una spline cubica (i punti della tabella sono fitti: 0,01–0,02 bohr
vicino al minimo). Con questa curva e lo stesso risolutore dei nuclei usato dalle palline:

| | curva esatta, massa ridotta (mₚ+mₑ)/2 | con mₚ/2 | misura |
|---|---|---|---|
| ωₑ armonica (curvatura) | 4404,5 | 4405,7 | 4401,21 (Huber & Herzberg) |
| ν(0→1) | 4162,2 | 4163,3 | 4161,17 (Dickenson et al. 2013) |
| Rₑ | 0,74143 Å | | 0,74144 Å |

Quindi con la curva **giusta** il nostro metodo dà già la misura entro +0,05 % (ν) e +0,08 % (ωₑ): il resto
(−1 / −2 cm⁻¹) sono le correzioni oltre Born–Oppenheimer (adiabatica, non adiabatica, relatività, QED), che la
letteratura dà proprio di questa grandezza (Komasa, Piszczatowski, Lach, Przybytek, Jeziorski, Pachucki,
J. Chem. Theory Comput. 7, 3105 (2011)).

---

## 2. La lista pre-registrata

| voce | ipotesi | risultato | causa? |
|---|---|---|---|
| B0 | il codice Python rifà i numeri del Rust | ωₑ 4338,84, ν 4128,81, Rₑ 0,74201 Å: identici | (prerequisito: PASS) |
| B1 | confronto armonica contro ωₑ "spettroscopico" | curva esatta: armonica 4405,7 vs 4401,21 = +0,10 % | **no** |
| B2 | massa ridotta (protone vs atomo) | cambia ωₑ di +0,03 % | **no** |
| B3 | correzioni oltre Born–Oppenheimer | ν(BO esatta) − misura = +2,1 cm⁻¹ = +0,05 % | **no** |
| B4 | rumore statistico | vedi sotto: il rumore non serve a spiegarla | **no** (vedi nota) |
| B5 | forma del fit / intervallo | fit sui punti **esatti**: ωₑ −1,61 %, ν −0,86 % | **SÌ** |
| B6 | passo della griglia dei nuclei | ν cambia di 0,016 % fra 1500 e 12000 punti | **no** |
| B7 | errore di tick che dipende da R | curvatura dei residui → ωₑ +0,03 ± 0,29 % | **no** |

### B1 — armonica e anarmonica

Il valore sperimentale ωₑ = 4401,21 viene da un fit dei livelli G(v) = ωₑ(v+½) − ωₑxₑ(v+½)² + ωₑyₑ(v+½)³.
Controllo: 4401,21 − 2·121,33 + (13/4)·0,8129 = 4161,19, e la misura del fondamentale è 4161,17: la tabella è
coerente. Sulla curva esatta l'armonica (dalla curvatura) è 4405,7, e un fit G(v) con v = 0..3 dà 4403,6 / 121,5
(v = 0..5: 4400,0 / 119,4): la differenza fra "armonica" e "ωₑ spettroscopico" è di pochi cm⁻¹, non di 62.
Con il vecchio fit, le palline davano anche ωₑxₑ = 101 (invece di ~120): il vecchio fit schiacciava tutta la
forma della buca. Con il fit corretto: ωₑ = 4386, ωₑxₑ = 117 (curva esatta, stesso conto: 4399 / 119).

### B4 — il rumore

Ho preso l'energia esatta alle 18 distanze, aggiunto rumore gaussiano grande come le nostre barre d'errore × 2,3
(il fattore che il χ² chiede) e rifatto lo stesso fit 2000 volte. Risultato: ωₑ = 4338 ± 35 (mediana ± metà
dell'intervallo 16–84 %), ν = 4130 ± 21. Le palline danno 4339 e 4129: **esattamente il centro**. Cioè: le palline
si comportano come la curva esatta più rumore, e lo spostamento di −62 cm⁻¹ è tutto distorsione del fit.
Il criterio pre-registrato "colpo di sfortuna statistico" risulta formalmente soddisfatto (0,02 σ dalla
distribuzione sintetica), ma la distribuzione sintetica è quella **già spostata** dal fit: non è sfortuna, è il fit.
Due note oneste:

* l'errore dichiarato prima, ±18 cm⁻¹, era **metà** di quello vero (±35–39 con quel fit): il ri-fit "bootstrap"
  del Rust lo sottostimava;
* i residui palline − esatto hanno χ²/punto = 4,9 con le barre d'errore nude (2,6 togliendo lo spostamento medio):
  le barre dei singoli punti sono ottimiste di un fattore ~1,6–2,2, come già scritto in `MOLECOLE.md`. Il pre-registrato
  "≤ 4" è soddisfatto solo dopo aver tolto lo spostamento medio.

### B5 — il colpevole

La figura (in alto a sinistra) mostra fit − esatto quando il fit riceve i punti esatti: il vecchio fit fa un'onda
di ±4·10⁻⁴ hartree, con un **dosso a R = 1,4 e due buche a 1,05 e 1,9**: cioè il fondo della buca viene più
largo e più piatto, la curvatura scende del 3,2 %, ωₑ (che va come la radice) dell'1,6 %. Con 18 punti e 6 parametri
quasi uguali fra loro (esponenziali con β, 2β, ... 5β) la formula non ha abbastanza libertà vicino al minimo, e
il β scelto dal χ² (0,68) è un compromesso con i punti lontani (0,9 e 5,5 bohr).

Il "controllo indipendente" che c'era in `MOLECOLE.md` (quartica locale: 4349) non era indipendente abbastanza:
anche lui era sotto, per il rumore (sui punti esatti la quartica locale dà 4405,8, sui nostri 4371 nella mia
versione).

| modello (giudicato sui dati esatti) | distorsione ωₑ | distorsione ν | spread ωₑ | palline: ωₑ | palline: ν |
|---|---|---|---|---|---|
| vecchio (5 exp, β griglia 0,04) | −1,61 % | −0,86 % | 39 | 4338,8 | 4128,8 |
| 5 exp, β continuo | −1,22 % | −0,63 % | 24 | 4341,1 | 4130,2 |
| **6 exp, β griglia 0,04** | **−0,03 %** | **−0,08 %** | **36** | **4396,6** | **4153,6** |
| 6 exp, β continuo | +0,02 % | −0,05 % | 37 | 4408,8 | 4160,3 |
| 5 exp, E∞ libero | −0,08 % | −0,07 % | 43 | 4395,4 | 4154,8 |
| 5 exp, E∞ libero, β continuo | −0,13 % | −0,10 % | 41 | 4401,7 | 4158,1 |

La regola di scelta (scritta in `h2vib_fix.py` prima di eseguirlo; onestà: in B5 avevo già visto due varianti sulle palline, 6 exp → 4397 ed E∞ libero → 4395): distorsione ≤ 0,3 % su ωₑ e ν,
con e senza rumore; fra questi, lo spread più piccolo → "6 exp, β griglia". Tutti i modelli non distorti danno
4395–4409 e 4154–4160: la scelta non cambia la conclusione.

---

## 3. Il numero corretto

Dalle palline (curva Born–Oppenheimer, fit a 6 esponenziali): ωₑ = 4396,6 ± 36,5, ν(0→1) = 4153,6 ± 20,1 cm⁻¹,
Rₑ = 0,74132 Å (era 0,74201; misura 0,74144), D₀ = 4,4716 eV (era 4,4721: non cambia).
Correzione Born–Oppenheimer → misura (presa dalla curva esatta con la stessa massa ridotta, §1): −3,3 e −1,1 cm⁻¹.

| | prima | adesso | misura | scarto |
|---|---|---|---|---|
| ωₑ | 4339 ± 18 | **4393 ± 36** | 4401,21 | −0,18 % (0,2 σ) |
| ν(0→1) | 4129 ± 16 | **4152 ± 20** | 4161,17 | −0,21 % (0,4 σ) |

**Risolto** (criterio pre-registrato: entrambi entro 2 σ). Resto −0,2 %, compatibile con zero. Per andare oltre
servono punti più precisi vicino al minimo (oggi ~1,5·10⁻⁴ hartree per punto: per un ±0,2 % su ωₑ servirebbe
~4 volte meno rumore, ~16 volte più tempo di palline) — non lo faccio qui.

Da correggere nel codice (non è un mio file, lo segnalo e basta): in `src/molecole/misura.rs` `Curva::fit` usa
5 esponenziali; con 6 (o E∞ libero) la distorsione sparisce. Anche l'errore del bootstrap del Rust (±18) andrebbe
controllato: lo stesso fit su punti esatti + rumore della stessa taglia ha uno spread di ±39. `docs/MOLECOLE.md` §7.2 e
§9 dicono ancora "causa non trovata": la causa è qui.
