# Atomo supercritico (Z oltre 137): l'1s che si tuffa nel mare di Dirac

Studio `supercrit_` — file: `lab/supercrit_*.py`, `lab/out/supercrit_*`. Stato: **chiuso** (13/10/2026). Figura: `lab/out/supercrit_figura.png`.

## 1. Il problema (cosa è noto, cosa è aperto)

- Con un nucleo **puntiforme** l'equazione di Dirac dà E_1s = mc² √(1 − (Zα)²): a Zα = 1 (Z ≈ 137) l'energia va
  a zero e oltre diventa immaginaria: la formula si rompe ("caduta sul centro").
- Con un nucleo **esteso** (raggio di qualche fm) l'1s continua a scendere e tocca −mc² (il bordo del "mare di
  Dirac") a una carica critica Z_cr ≈ 170–173: Pomeranchuk–Smorodinsky (J. Phys. USSR 9, 97, 1945, stima troppo
  alta); Gershtein–Zeldovich (Lett. Nuovo Cimento 1, 835, 1969); Pieper–Greiner (Z. Phys. 218, 327, 1969);
  Müller–Rafelski–Greiner (Z. Phys. 257, 62 e 183, 1972); Zeldovich–Popov (Sov. Phys. Usp. 14, 673, 1972).
  Valori moderni: Z_cr ≈ 173 (Maltsev et al., PRL 123, 113401, 2019, con modello di nucleo realistico);
  ≈ 170 con R = 1,2 (2,5 Z)^{1/3} fm (Popov et al., arXiv:2209.00125).
- Oltre Z_cr lo stato non è più legato: diventa una **risonanza** dentro il continuo negativo. Se il guscio K è
  vuoto, un elettrone del mare lo riempie e il buco che resta scappa: **emissione spontanea di un positrone**;
  il vuoto diventa "carico". Le larghezze Γ della risonanza sono piccole (eV–keV): esempio, tabella esatta di
  Godunov, Machet, Vysotsky (EPJC 77, 782, 2017) per un nucleo "a guscio" di raggio R = 0,031 ħ/mc (12 fm):
  Z_cr ≈ 177,5, Γ(184) = 5,2·10⁻⁴ mc² = 265 eV, Γ(190) = 6,1·10⁻³ mc², Γ(232) = 0,19 mc².
  Vicino alla soglia Γ è soppressa dalla barriera di Coulomb vista dal positrone: Γ ∝ exp(−2πη),
  η = Zα E/p (Zeldovich–Popov).
- **Mai osservato.** Nessun nucleo ha Z > 118; si usano collisioni lente U+U (Z = 184 per ~10⁻²¹ s). I picchi di
  positroni visti a GSI negli anni '80 (EPOS, ORANGE) non sono stati confermati (APEX, anni '90). Oggi si
  propongono collisioni U–U a energie sotto la barriera con traiettorie scelte (Maltsev et al. 2019; Popov et
  al.), esperimenti previsti a GSI/FAIR e HIAF. È quindi una previsione standard della QED in campo forte
  mai verificata: il problema aperto è **sperimentale**, non teorico.
- **Analogo in materiali**: il "collasso atomico" in grafene (Shytov–Katsnelson–Levitov 2007; Pereira–Nilsson–
  Castro Neto 2007) è stato osservato (Wang et al., Science 340, 734, 2013). Lì il reticolo di carbonio fa da
  taglio vicino all'impurità, proprio come farà il nostro reticolo.
- **Quantum walk / automi di Dirac**: ci sono walk con campo elettrico e con Coulomb (Di Molfetta–Debbasch–
  Brachet 2014; walk di Dirac con Coulomb per la ricerca spaziale, arXiv:2210.13920), automi di Dirac 1-D con
  scattering da potenziale (Bisio–D'Ariano–Tosini 2013), simulazioni su ioni intrappolati (Huerta Alderete et al.
  2020). **Non abbiamo trovato** un lavoro che faccia con un walk/automa di Dirac il tuffo dell'1s oltre Z_cr e il
  decadimento della risonanza nel tempo. Le simulazioni tempo-dipendenti della risonanza esistono, ma con
  l'equazione di Dirac ordinaria (Müller, Reinhardt, Greiner anni '70–'80; Maltsev, Shabaev, Tupitsyn 2015–2019).
  Quindi: se funziona, è **nuovo come metodo, non come fisica**.

(Letteratura controllata online il 13/10/2026: arXiv:2008.05005, arXiv:2209.00125, arXiv:1707.07497, tesi di
R.V. Popov 2024 — U+U: R_cr = 32,7 fm, larghezza naturale < 2 keV a R = 16 fm. Gli altri riferimenti storici sono
citati a memoria.)

## 2. Previsione prima dei dati

Vedi `lab/out/supercrit_prereg.json` (generato da `lab/supercrit_prereg.py`). In breve:

1. Nucleo puntiforme, Z ≤ 137: la regola ritrova Dirac, ma con errore che si chiude come a^{2γ} (γ = √(1−(Zα)²)):
   a Z = 137, 2γ = 0,05, cioè **in pratica non converge** (limite atteso).
2. Nucleo puntiforme, Zα > 1: il reticolo non ha punti dentro r = a/2, quindi fa da **nucleo finto grande ~a**.
   Previsione: l'1s tocca −mc² a un Z_cr(a) finito, che scende verso 137 come 1 + π²/(2 ln²(1/(c·a))) con c fisso.
   È un artefatto (dichiarato): non è fisica.
3. Nucleo esteso (DEBITO: il raggio): Z_cr = 177,46 per il guscio R = 0,031 (letteratura), 170–173 per la sfera.
4. Oltre Z_cr: la risonanza si vede sia dalla fase (salto di π) sia dal tempo (lo stato si svuota con legge
   esponenziale e un'onda di carica positiva esce). Γ entro il 10 % della tabella di Godunov et al.


## 3. In due righe (per Valerio)

- **Con un nucleo vero (esteso) la regola fa tutto giusto**: l'1s si tuffa nel mare di Dirac a **Z_cr = 177,463**
  (letteratura 177,46 per lo stesso nucleo a guscio di 12 fm), a **169,706** per una sfera uniforme di
  1,2 (2,5Z)^{1/3} fm (letteratura ≈ 170), e oltre Z_cr diventa una risonanza con la larghezza giusta, da 1 eV
  (Z = 180) a 120 keV (Z = 232), **entro lo 0,2–1 %** della tabella esatta di Godunov et al. 2017.
- **La parte nuova**: lasciando andare la regola nel tempo si vede davvero l'emissione: lo stato si svuota con
  legge esponenziale, un guscio di carica esce a ~0,8 c (il **positrone**, +e) e al centro resta il "vuoto
  carico" (−e). Il ritmo dello svuotamento e lo spettro dell'onda uscita danno la stessa Γ della fase (entro lo
  0,3 %, ma con un'analisi corretta dopo aver visto i dati, vedi §5).
- **Il limite trovato**: con nucleo **puntiforme** la regola non sa cosa fare vicino a Zα = 1. Il reticolo fa
  tuffare l'1s già a **Z ≈ 131,4** (non a 137!), e questo numero **non si muove** facendo il reticolo più fine.
  È un artefatto della regola, non fisica. Per Z ≤ 120 invece ritrova Dirac a 10⁻⁵ o meglio.
- È **validazione**, non scoperta di fisica: la fisica del tuffo è nota dal 1969–1972. Nuovo è il **metodo**
  (un walk di Dirac che fa tutto il film: tuffo, risonanza, positrone che esce) e il **limite** della regola a
  Zα → 1 con nucleo puntiforme.

## 4. Cosa abbiamo fatto

La regola è quella congelata (`zz_core.RadialWalk`, onda parziale s₁/₂, κ = −1): corri a c (χ₊ verso
l'esterno, χ₋ verso l'interno, rimbalzo in r = 0), poi ruota con exp(−ia(mσ_z + κ/r σ_x + V)). L'unica cosa
cambiata è il potenziale V(r) (fase per tick): puntiforme −Zα/r, guscio (−Zα/R dentro R), sfera uniforme.
`lab/supercrit_core.py` controlla che con il nucleo puntiforme sia **identica** a RadialWalk (differenza 0,0).

Strumenti (tutti nuovi, in `supercrit_core.py`):
1. **Livelli**: autovettori dell'operatore di un tick W (shift-invert), si segue l'1s in Z finché tocca −mc².
2. **Scansione di fase**: l'equazione stazionaria W χ = e^{−iωa} χ si risolve come ricorrenza dall'origine
   verso fuori (soluzione regolare); a r = 25 si scompone sulle due onde di Bloch locali; la fase relativa salta
   di 2π attraverso la risonanza → fit di Breit–Wigner (posizione E_r, larghezza Γ).
3. **Film**: lo stato quasi-legato (la regola in una scatola chiusa a metà della barriera di Coulomb) si mette
   nella scatola grande (r < 70, spugna assorbente da r = 40) e si lascia andare: si misura il peso dentro,
   l'onda che passa a r = 30 e il suo spettro.
4. **Riferimento** (solo per controllo): l'equazione di Dirac radiale nel continuo, integrata con scipy, con lo
   stesso V(r). Ritrova la tabella di Godunov (Z_cr 177,463; E(170) = −0,66392; Γ(190) = 6,07·10⁻³).

Tutto in unità ħ = c = m_e = 1 (lunghezza 386 fm, energia 511 keV). Nessun parametro libero nella regola.

## 5. Numeri ed esami (contro `lab/out/supercrit_prereg.json`; tutto in `lab/out/supercrit_exams.json`)

### 5.1 V1 — nucleo puntiforme, Z ≤ 137, contro la formula di Dirac (validazione)

Passi a = 0,03125 … 0,00049 (7 reticoli), estrapolazione E₀ + K a^{2γ} + B a² sui 3 più fini.

| Z | errore relativo (sul legame) al reticolo più fine | dopo a → 0 | esame |
|---|---|---|---|
| 20 | −6·10⁻⁸ | −3·10⁻¹⁰ | V1a PASS |
| 60 | −1,5·10⁻⁶ | −1,2·10⁻⁹ | V1a PASS |
| 92 | −3,2·10⁻⁵ | −2,9·10⁻⁸ | V1a PASS |
| 110 | −3,9·10⁻⁴ | −2,7·10⁻⁷ | V1a PASS |
| 120 | −2,8·10⁻³ | +9,7·10⁻⁶ | V1b PASS |
| 130 | −0,128 | +2,1·10⁻³ | V1c PASS (soglia 10⁻²) |
| 137 | +0,21 (sul livello giusto, vedi sotto) | non estrapolabile | V1c **FAIL** (previsto: limite) |

A Z = 137 (2γ = 0,05) la regola converge così lentamente che non si può estrapolare; in più (scoperto qui) il
livello più basso del reticolo non è l'1s di Dirac: per a ≤ 0,016 si è già tuffato sotto −mc² (§5.2) e il livello
che tende all'1s vero (0,0229 mc²) sta ancora a 0,231 con a = 0,00049.

### 5.2 C1 — nucleo puntiforme oltre Zα = 1: il reticolo come nucleo finto (LIMITE)

Previsione pre-registrata: il reticolo fa da "taglio" di dimensione ~a, Z_cr(a) > 137 e scende verso 137 come
nel continuo (1 + π²/(2 ln²R), verificato col nostro integratore: a R = 10⁻⁶ 140,63). **Sbagliata.**

| a [ħ/mc] | 0,1 | 0,03 | 0,01 | 0,003 | 0,001 | 3·10⁻⁴ | 10⁻⁴ | 3·10⁻⁵ |
|---|---|---|---|---|---|---|---|---|
| Z_cr(a) | 146,75 | 137,80 | 134,43 | 132,78 | 132,08 | 131,71 | 131,55 | 131,45 |

Z_cr(a) **scende sotto 137** e si ferma: fit Z* + c a^q → **Z* = 131,37** (ζ* = Z*α = 0,9587), q = 0,62.
C1a FAIL, C1b FAIL (per a ≤ 0,01 non esiste un nucleo esteso equivalente: un nucleo esteso sposta Z_cr *sopra*
137, mai sotto).

**Perché** (spiegazione, controllata solo in parte): per Zα > √3/2 (Z > 118,7) entrambe le soluzioni di Dirac
vicino all'origine, r^{+γ} e r^{−γ}, sono a quadrato integrabile: la fisica non è fissata dall'equazione da
sola ma da cosa succede nel nucleo. Il reticolo, vicino all'origine, è un problema senza scala che sceglie un
miscuglio fisso β(ζ) delle due soluzioni (lo stesso numero β misurato dallo studio `nucleo1s_`, che cresce con
Z). Se β ha un polo a ζ*, il livello tocca −mc² a ζ − ζ* ∝ a^{2γ*}: esponente previsto 2γ* = 0,569, misurato 0,62
(9 % di scarto). Tra Z* e 137 il reticolo ha un livello **in più** che si tuffa, e i livelli veri (1s, 2s, …)
arrivano da sopra, lentissimi. Morale: **con nucleo puntiforme e Z > 131 la regola dà fisica sbagliata a
qualunque passo**; non importa per gli elementi veri (Z ≤ 118), importa per chiunque voglia usare la regola
"nuda" in campo supercritico. Con nucleo esteso risolto (R ≫ a) il problema sparisce (sotto).
Che un reticolo sposti la carica critica è noto in grafene (dipende dal taglio); il valore Z* = 131,4 di questa
regola è nostro.

### 5.3 C2, C3 — nucleo esteso: la carica critica (validazione)

| nucleo | regola a → 0 (a = 0,004 … 0,00025) | nostro Dirac continuo | letteratura | esame |
|---|---|---|---|---|
| guscio R = 0,031 ħ/mc = 12 fm, Z_cr | **177,463** | 177,463 | 177,46 (Godunov et al. 2017, interp. tra 177 e 178) | C2a PASS |
| guscio, E_1s a Z = 170 | −0,66392 mc² | −0,66392 | −0,664 | C2b PASS |
| sfera uniforme R = 1,2 (2,5Z)^{1/3} fm (9,0 fm), Z_cr | **169,706** | 169,706 | ≈ 170 (Popov et al.); ≈ 173 con modelli di nucleo realistici | C3 PASS |

Convergenza: l'errore va come a (il potenziale ha uno spigolo in R), es. guscio 177,734 / 177,605 / 177,549 /
177,507 / 177,485.

### 5.4 G1, G3 — la risonanza oltre Z_cr dalla scansione di fase (validazione)

Γ in unità mc² e in eV (regola, a → 0 da a = 0,004 … 0,0005), contro la tabella "esatta" di Godunov et al.:

| Z | E_r regola | E_r lett. | Γ regola | Γ regola [eV] | Γ lett. | scarto Γ |
|---|---|---|---|---|---|---|
| 180 | −1,1270 | −1,127 | 1,850·10⁻⁶ | 0,95 | 1,846·10⁻⁶ | +0,2 % |
| 181 | −1,1785 | −1,178 | 1,915·10⁻⁵ | 9,8 | 1,895·10⁻⁵ | +1,1 % |
| 182 | −1,2311 | −1,231 | 8,446·10⁻⁵ | 43 | 8,336·10⁻⁵ | +1,3 % |
| 184 | −1,3389 | −1,338 | 5,195·10⁻⁴ | 265 | 5,182·10⁻⁴ | +0,3 % |
| 186 | −1,4505 | −1,450 | 1,588·10⁻³ | 812 | 1,581·10⁻³ | +0,4 % |
| 188 | −1,5660 | −1,565 | 3,427·10⁻³ | 1 751 | 3,414·10⁻³ | +0,4 % |
| 190 | −1,6854 | −1,685 | 6,079·10⁻³ | 3 107 | 6,060·10⁻³ | +0,3 % |
| 195 | −2,0014 | −2,001 | 1,618·10⁻² | 8 267 | 1,606·10⁻² | +0,7 % (+0,19 %*) |
| 200 | −2,3420 | −2,341 | 3,088·10⁻² | 15 780 | 3,034·10⁻² | +1,8 % (+0,17 %*) |
| 210 | −3,0961 | −3,094 | 7,698·10⁻² | 39 340 | 6,928·10⁻² | **+11 %** (+0,18 %*) |
| 232 | −5,0662 | −5,057 | 0,2355 | 120 300 | 0,1928 | **+22 %** (+0,21 %*) |

(*) finestra di fit ±3Γ invece di ±12Γ: **post-hoc**, dichiarato (`lab/supercrit_posthoc_wide.py`).

- G1a (Γ entro 10 % per Z ≥ 182): **FAIL** con il fit pre-registrato a Z = 210 e 232. Non è la regola: lo stesso
  fit (±12Γ, fondo quadratico) applicato alla fase del **Dirac continuo** dà gli stessi scarti (0,0740 a Z = 210);
  per risonanze larghe la "Breit–Wigner su finestra larga" non misura il polo. Con ±3Γ regola e continuo
  coincidono con Godunov entro lo 0,2 %.
- G1b (Z = 180, Γ = 1 eV, entro 30 %): PASS (+0,2 %). G1c (E_r entro 1 %): PASS (≤ 0,18 %).
- G3 (legge di Gamow vicino alla soglia, ln Γ = c − s·2πη, η = Zα|E|/p, Z = 180…188): s = **1,12**, PASS (soglia
  ±0,2). Il positrone deve passare sotto la barriera di Coulomb: è questo che rende Γ così piccola vicino a Z_cr.

### 5.5 G2, D1 — il film: lo stato si svuota e il positrone esce (la parte nuova)

Z = 190, 200, 210, 232, a = 0,002. Pesi e spettri:

| Z | Γ dal tempo (pre-reg.) / Γ dalla fase | Γ dal tempo con "pavimento" (post-hoc) / Γ fase ±3Γ | larghezza dello spettro dell'onda uscita / Γ fase | picco dello spettro vs E_r |
|---|---|---|---|---|
| 190 | 0,964 | **0,998** | 1,000 | −1,67424 vs −1,67419 |
| 200 | 0,936 | **0,999** | 1,000 | −2,32742 vs −2,32697 |
| 210 | 0,854 | **0,997** | 0,999 | −3,0775 vs −3,0769 |
| 232 | 0,779 | **0,999** | 0,999 | −5,0357 vs −5,0376 |

- G2 (pre-registrato: pendenza di ln P tra 0,7 e 0,05, entro 5 %): **FAIL**. Due cause, trovate guardando i
  dati: (i) lo stato di partenza ha lo 0,5 % di peso sui livelli **legati veri** (2s, 3s, …), che non decadono e
  fanno un "pavimento" che piega la coda; (ii) per Z = 210, 232 il riferimento di fase era quello del fit largo
  (§5.4). Con P(t) = A e^{−Γt} + C (post-hoc, `lab/supercrit_posthoc_decay.py`) tempo e fase coincidono entro lo 0,3 %.
- D1 (Z = 200, spettro dell'onda a r = 30: picco entro 1 %, larghezza entro 15 %): **PASS** (picco 0,02 %,
  larghezza −1,5 % rispetto al fit largo, 0,0 % rispetto a quello stretto).

Che cosa si vede (figura, pannelli d): un guscio di carica parte dal nucleo e corre verso fuori a ~0,8 c, con
l'energia giusta del positrone emesso (|E_r| − mc² = 1,34 mc² ≈ 685 keV di energia cinetica a Z = 200); il peso
interno cala come e^{−Γt}. Nella lingua della teoria delle buche: lo stato tuffato vuoto è un buco nel mare;
l'onda che esce **è** il positrone (+e) e il nucleo si ritrova circondato da un elettrone in più (il vuoto carico,
−e). È la stessa equazione a una particella usata da Müller–Rafelski–Greiner; non c'è il calcolo a molti corpi
(lo dichiariamo).

### 5.6 Riassunto PASS/FAIL

| esame | esito |
|---|---|
| V1a, V1b, V1c(130) | PASS |
| V1c(137) | FAIL (previsto, limite) |
| C1a, C1b | **FAIL** (previsione sbagliata: Z_cr(a) → 131,37, non 137) |
| C2a, C2b, C3 | PASS |
| G1a | FAIL a Z = 210, 232 col fit pre-registrato (difetto del fit, non della regola: post-hoc 0,2 %) |
| G1b, G1c, G3, D1 | PASS |
| G2 | FAIL pre-registrato (−4…−22 %); post-hoc 0,3 % |

## 6. Validazione, scoperta, limite

- **Validazione**: Dirac fino a Z = 130; Z_cr 177,463 e 169,706; E_r e Γ della risonanza da Z = 180 (1 eV) a
  Z = 232 (120 keV) entro lo 0,2–1 %; legge di Gamow. Tutti numeri noti.
- **Nuovo come metodo** (non come fisica): il tuffo e l'emissione spontanea fatti per intero da un walk/automa
  di Dirac a tempo discreto (cercato in letteratura: walk con Coulomb e campo elettrico esistono, il tuffo
  supercritico no, per quanto abbiamo trovato). Tre misure indipendenti di Γ (fase, tempo, spettro) concordano.
- **Limite** (il risultato più interessante per noi): la regola con nucleo puntiforme ha una sua carica critica
  finta **Z* = 131,37 ± ~0,05** (indipendente dal passo), sotto 137. Per Z tra 131 e 137 i livelli s della regola
  non convergono a Dirac in tempi umani. Per usare la regola in campo forte serve il nucleo esteso (che però è un
  debito).
- **Non fatto**: QED (polarizzazione del vuoto, autoenergia: spostano Z_cr di circa uno; c'è ora
  `docs/VUOTO_POLARIZZA.md` di un altro thread, non usato qui), schermo degli altri elettroni, nucleo che si
  muove (collisioni U+U: la cosa che serve davvero agli esperimenti), calcolo a molti corpi del vuoto carico.

## 7. Debiti (cose messe a mano)

1. **Raggio del nucleo** (non è un ingresso ammesso): guscio R = 0,031 ħ/mc copiato dal lavoro di Godunov et al.
   per poterli confrontare; sfera R = 1,2 A^{1/3} fm con **A = 2,5 Z** (nessun nucleo con Z > 118 esiste: anche la
   massa è inventata). La regola non sa fare il nucleo.
2. Forma del nucleo (guscio o sfera uniforme): a mano.
3. Solo l'onda s₁/₂ (κ = −1) nella forma radiale: la simmetria sferica è imposta. (L'avviso sui 4 doppioni della
   regola 3-D a cubo, `docs/BATTIMENTO_1S.md`, non riguarda questo studio: qui si usa solo la forma radiale.)
4. Stato iniziale del film = stato della regola in una scatola chiusa a metà barriera (scelta di comodo).
5. Spugna assorbente (non unitaria) oltre r = 40.
6. Nessuna QED, nessuno schermo elettronico, nucleo fermo.

**Parametri liberi**: 0 nella regola. Modelli d'errore per a → 0 fissati prima (a^{2γ}, a², a). Il fit
Z* + c a^q del §5.2 ha 3 parametri ed è post-hoc (descrittivo).

## 8. Come rifare le corse (1 core, ognuna < 10 min)

    python lab/supercrit_core.py              # controlli (identità con RadialWalk)
    python lab/supercrit_prereg.py            # pre-registrazione (commit 972f466)
    python lab/supercrit_run_v1.py            # V1 (1 min)
    python lab/supercrit_run_c1.py            # C1 (a = 3e-5 da solo ~8 min; rilanciare per continuare)
    python lab/supercrit_run_c2.py shell; python lab/supercrit_run_c2.py sphere; python lab/supercrit_run_c2.py cont
    python lab/supercrit_run_g1.py            # G1 (2 min)
    python lab/supercrit_run_g2.py 190; ... 200; ... 210; ... 232   # film (190: 5 min)
    python lab/supercrit_posthoc_wide.py; python lab/supercrit_posthoc_decay.py   # post-hoc dichiarati
    cd lab && python supercrit_exams.py && python supercrit_figure.py

Emendamenti rispetto alla pre-registrazione (dichiarati): A1 passo massimo 0,5 nel tracciamento in Z di C1 (il
primo passaggio con passo 4 saltava sopra il tuffo a a = 10⁻⁴, 3·10⁻⁵: risultati 150,45 e 147,43 scartati, nel
JSON sotto `first_pass_discarded`); G1 calcolato a tutti e 4 i passi anche per Z = 181…189 (pre-registrato: 2);
seme dell'energia nella scansione = stato in scatola; il valore continuo di C3/C2 è stato copiato dal log nel JSON
(scrittura persa per una corsa concorrente). I due controlli post-hoc sono in file separati.
