# Perché la regola fa tuffare l'1s a Z* = 131,37 e non a 137

Studio `zstar_` — file: `lab/zstar_*.py`, `lab/out/zstar_*`. Stato: **chiuso** (13/10/2026, giro 2).
Figura: `lab/out/zstar_figura.png` (copia in `/mnt/project-files/problemi-aperti/zstar_figura.png`).
Punto di partenza: `docs/SUPERCRITICO.md` §5.2 (giro 1).

## 0. In due righe (per Valerio)

- **Trovato il perché, con una formula senza parametri che dà Z* = 131,3614** (la regola, estrapolata, dava
  131,37 ± 0,05). Vicino al nucleo la regola non ha scala: nella cella i gli angoli sono κ/(i+½) e Zα/(i+½), e il
  passo a sparisce. Quel pezzetto di reticolo sceglie da solo **quale soluzione di Dirac** usare vicino al centro
  (la buona r^{+γ} o la cattiva r^{−γ}). A Z* = 131,36 sceglie **solo la cattiva**: lì il reticolo ha uno "stato di
  contatto" a energia zero, e da lì in su l'1s si tuffa. La formula di raccordo (reticolo dentro, Dirac fuori)
  prevede anche **tutta la curva** Z_cr(a) (regola base: scarto 0,004 a a = 10⁻⁴, 0,001 a 3·10⁻⁵; a passo grosso
  sbaglia dove contano i termini O(a)).
- **Non è universale: è un artefatto della condizione al centro.** Cambiando un dettaglio che la fisica non fissa
  (il segno/la fase del rimbalzo in r = 0, o dove cadono le celle rispetto al nucleo) Z* va da 95 a 131, oppure
  **sparisce**. Le previsioni fatte prima delle corse (cieche) tornano al limite a → 0 entro 10⁻⁶–3·10⁻³ in Z (un esame, F2, fallisce a passo grosso: §5).
- **Rimedio dentro le regole**: mettere le celle a r = (i+1)a (il nucleo "su una cella", come un nucleo che sta
  su un nodo del reticolo 3-D) o a r = (i+√3/2)a (prima cella alla distanza dei primi vicini del cubo 3-D, che ha il
  nucleo al centro di un cubetto). Allora non c'è nessun Z* finto: l'1s si tuffa **sopra** 137 e scende verso 137
  come un nucleo di raggio ~a (145,46 a a = 10⁻⁴; 138,7 a a = 10⁻⁹). È il limite giusto, ma lentissimo
  (logaritmico): **137 si raggiunge solo a a → 0**, come per qualunque nucleo che si rimpicciolisce. La scelta della
  posizione è un debito (una scelta geometrica).
- **Grafene**: stessa matematica (Dirac + 1/r, la carica critica è dove γ = √(κ² − ζ²) diventa immaginario, il
  reticolo fissa la condizione al centro), fisica diversa nel segno: lì il taglio del reticolo **alza** la carica
  critica (con la gap) o non la sposta (senza gap), qui il reticolo la **abbassa** perché fa da attrazione di
  contatto. Dettagli e numeri in §6.
- È **validazione + spiegazione di un limite** (teoria nota: estensioni autoaggiunte, "caduta sul centro",
  anomalia di scala), applicata alla nostra regola con una formula nuova e controllata. Nessuna fisica nuova.

## 1. Il problema

Con il nucleo puntiforme l'equazione di Dirac dà E_1s = mc² √(1 − ζ²), ζ = Zα: a ζ = 1 (Z = 137,036) l'energia va
a zero e oltre la formula si rompe. Nel giro 1 la regola radiale congelata (`zz_core.RadialWalk`, κ = −1) ha fatto
una cosa diversa: l'1s tocca −mc² a una carica Z_cr(a) che **scende sotto 137** e si ferma a **Z* = 131,37**
quando il passo a → 0 (146,75 a a = 0,1 … 131,45 a a = 3·10⁻⁵). L'ipotesi del giro 1: il reticolo, vicino
all'origine, sceglie un miscuglio fisso delle soluzioni r^{+γ} e r^{−γ}, e quel miscuglio ha un polo. Qui la
controllo, la trasformo in una formula, la metto alla prova con varianti della regola mai girate prima.

## 2. Cosa è noto (letteratura)

- **Potenziali singolari e condizione al contorno.** Case (Phys. Rev. 80, 797, 1950): con potenziali che vanno come
  1/r² (Schrödinger) o 1/r (Dirac) "troppo forti", l'equazione da sola non decide la fisica vicino al centro: serve
  una condizione in più (un'estensione autoaggiunta dell'Hamiltoniana). Per Dirac–Coulomb: per ζ ≤ √3/2 (Z ≤ 118,7)
  c'è una sola scelta; per √3/2 < ζ < 1 entrambe le soluzioni r^{±γ} sono a quadrato integrabile e c'è una
  famiglia di scelte a un parametro, con una "distinta" (la r^{+γ}, quella che si ottiene da un nucleo esteso che
  si rimpicciolisce); per ζ > 1 nessuna scelta privilegiata ("caduta sul centro"). Testi: Gitman, Tyutin, Voronov,
  *Self-adjoint Extensions in Quantum Mechanics* (Birkhäuser 2012), cap. su Dirac–Coulomb; Esteban, Lewin, Séré
  (Rev. Mat. Iberoam. 2019 e lavori seguenti) sull'estensione distinta per ζ < 1; Hogreve (J. Phys. A 46, 025301,
  2013) sull'operatore oltre ζ = 1.
- **Anomalia di scala / "conformality lost".** Per 1/r² (Schrödinger) e 1/r (Dirac senza massa) il problema non ha
  scala; la condizione al centro è un accoppiamento di contatto g che "corre" con la scala. Per ζ < 1 ci sono due
  punti fissi (r^{+γ} stabile, r^{−γ} instabile), che si fondono a ζ = 1 e spariscono (Kaplan, Lee, Son,
  Stephanov, PRD 80, 125005, 2009; prima Beane et al. PRA 64, 042103, 2001; Camblong, Ordóñez et al. PRL 87,
  220402, 2001). Un taglio a corta distanza mette il valore "nudo" di g: se cade **oltre** il punto fisso instabile
  nasce uno stato legato anche sotto la soglia — è quello che succede alla nostra regola (§4).
- **Grafene.** Gli elettroni del grafene sono Dirac senza massa in 2-D; per un'impurezza carica la soglia è
  β = Zα_eff = ½ (canale j = ½), "collasso atomico" oltre (Shytov, Katsnelson, Levitov PRL 99, 236801, 2007;
  Pereira, Nilsson, Castro Neto PRL 99, 166802, 2007). Con una gap (massa) il taglio del reticolo **alza** la
  carica critica: β̃_c > ½ (Pereira, Kotov, Castro Neto PRB 78, 085101, 2008: potenziale tagliato a r₀ = distanza
  C–C, "necessariamente g̃_c > g_c"). Misure: Wang et al. (Science 340, 734, 2013), dimeri di Ca: 1–2 dimeri sotto
  soglia, 3 in transizione, 4–5 supercritici (Z/Z_c = 1,4 ± 0,2 a 3 dimeri), e cambiare il taglio r₀ tra 0,5 e
  1,5 nm non cambia i fit; Mao et al. (Nature Phys. 12, 545, 2016), vacanza carica: risonanze di collasso per
  β ≈ 0,95–1,29, sotto soglia β ≈ 0,45; la vacanza ha anche un suo modo a energia zero, che dà picchi in più ma
  secondo gli autori non sposta la soglia ½. (Letti online il 13/10/2026: arXiv 0803.4195, 0706.2872, 1510.02890,
  1508.07667; il resto citato a memoria.)
- **Nel nostro repo**: `docs/NUCLEO_1S.md` (thread principale, su origin/master) ha già misurato il miscuglio β(Z)
  della regola per Z ≤ 92 (β = 0,10301 a Z = 60) e mostrato che l'errore dell'1s va come β a^{2γ}. Questo studio
  lo estende fino al polo.
- **Aperto?** Il fatto generale (un reticolo fissa la condizione al centro e può spostare la carica critica) è
  noto. Il numero Z* = 131,36 della nostra regola e la formula che lo dà sono nostri; non ho trovato un calcolo
  equivalente per un walk/automa di Dirac.

## 3. La teoria (prima dei dati nuovi)

**Vicino al nucleo la regola non ha scala.** Nella cella i (r_i = (i + s)a, s = ½ nella regola congelata) il giro
per tick ha angolo centrifugo κa/r_i = κ/(i+s) e fase di Coulomb ζa/r_i = ζ/(i+s): **a sparisce**. La massa
(ma) e l'energia (ωa) per r ≪ 1 sono trascurabili. Quindi le prime celle sono uguali per ogni reticolo.

**Le matrici di trasferimento.** A energia zero e massa zero lo stato stazionario della regola si propaga cella per
cella con una matrice 2×2 (calcolata a mano dalla moneta, `lab/zstar_core.py`):

    v_{i+1} = T_i v_i,   T_i = [[e^{iφ_i}/cos θ_i, −tan θ_i], [−tan θ_i, e^{−iφ_i}/cos θ_i]],
    θ_i = κ/(i+s),  φ_i = ζ/(i+s),   v_0 = (refl0, 1)   (refl0 = −1: il rimbalzo χ− → −χ+ in r = 0)

det T_i = 1 e T_i sta in SU(1,1): per questo il risultato ha una struttura reale (conservazione della corrente) e
il miscuglio è un numero reale. Lontano dal nucleo T_i → 1 + K/ρ con K = [[iζ, −κ], [−κ, −iζ]], autovalori ±γ:
la soluzione diventa A ρ^{+γ} + B ρ^{−γ}, cioè **Dirac con la condizione al contorno τ = β a^{2γ}**, β = B/A.

**Perché la regola converge di solito**: τ = β a^{2γ} → 0 per a → 0 (errore dell'1s ∝ a^{2γ}, come in NUCLEO_1S).
**Perché a Z* no**: se a un certo ζ* il reticolo sceglie **solo** r^{−γ} (A(ζ*) = 0, β = ∞), allora vicino a ζ*,
β ≈ b/(ζ* − ζ) e τ resta finito anche per a → 0 quando ζ − ζ* ~ a^{2γ*}. Il livello si tuffa a ζ_cr(a) = ζ* +
O(a^{2γ*}) → **ζ* per a → 0**.

**La formula per ζ*** (equazione trascendente, senza parametri):

    F(ζ*) = lim_{n→∞} n^{−γ} ℓ₊ · T_{n−1} ⋯ T_1 T_0 (refl0, 1)ᵀ = 0     (ℓ₊ = autovettore sinistro di K per +γ)

Si calcola come "wronskiano" tra la soluzione dal centro e quella che decade (ricorrenza all'indietro, stabile) a
2·10⁵ celle: stabile a 10⁻¹³ cambiando la lunghezza (5·10⁴ … 10⁶) e la cella di valutazione. **Non ho trovato una
forma chiusa**: il prodotto infinito di matrici che non commutano non si somma (al primo ordine, trascurando la non
commutatività delle prime celle, si otterrebbe refl0 = (γ − iζ)/κ, cioè ζ* = 0 per la regola: sbagliato; sono
proprio le prime 2–3 celle a fare tutto).

**La curva intera Z_cr(a) (raccordo)**: dentro, la soluzione del reticolo (A₀, B₀ dal fit su ρ ∈ [40, 10⁵]);
fuori, Dirac del continuo con m = 1 a E = −1 (solo come strumento di raccordo); Z_cr(a) è dove la soluzione
raccordata ha lo stato di soglia. Trascura i termini O(a) (massa ed energia nelle prime celle).

**Corollari esatti** (dimostrati, poi controllati): (i) T(−κ) = σ_z T(κ) σ_z, quindi **Z*(κ, refl0) = Z*(−κ, −refl0)**;
(ii) l'ordine dei due mezzi passi (W = M S oppure S M) non cambia nessun autovalore (matrici simili);
(iii) spezzare la moneta (Trotter) non cambia nulla a m = 0: Z* identico, solo O(a) diverso; (iv) il rimbalzo
esterno refl1 (a r = 40) non conta.

Controllo di partenza: β(Z = 20, 60, 92) = 0,06628, 0,10301, 0,22037 — gli stessi di NUCLEO_1S (0,06629,
0,10301, 0,22037), calcolati da un altro thread in un altro modo.

## 4. La previsione (pre-registrata, commit 448d9c8) e il perché di 131,37

`lab/out/zstar_prereg.json` (da `lab/zstar_prereg.py`):

| variante | cosa cambia | ζ* (formula) | **Z*** | 2γ* |
|---|---|---|---|---|
| base (regola congelata) | — | 0,9585907 | **131,361** | 0,569 |
| R+1 | rimbalzo χ− → +χ+ | 0,918721 | 125,898 | 0,79 |
| R+i / R−i | rimbalzo con fase ±i | 0,69414 / 0,94246 | 95,122 / 129,151 | 1,44 / 0,67 |
| S0.25 | celle a (i + ¼)a | 0,756697 | 103,695 | 1,31 |
| S0.75, S0.866, S1 | celle a (i + ¾)a, (i + √3/2)a, (i + 1)a | **nessuna** | — | — |
| K+1 (p½) | κ = +1 | 0,918721 | 125,898 (= R+1, dualità) | 0,79 |
| K+1R+1 | κ = +1, rimbalzo +1 | 0,9585907 | 131,361 (= base) | 0,57 |
| K−2 (p3/2) | κ = −2 | 1,503406 | 206,021 | 2,64 |
| K+2 (d3/2) | κ = +2 | 0,306877 | **42,053** | 3,90 |

Scansioni (figura, pannello c): con s = ½ **ogni** fase del rimbalzo dà un Z* (da 95 a 131); con refl0 = −1 il Z*
sparisce per s > s_c ≈ 0,618 (la radice esce da ζ < 1). Per κ = ±½ (l'analogo radiale del grafene, §6) nessuna
radice sotto ζ = ½.

**Il perché di 131,37, detto semplice**: la regola fa rimbalzare l'onda in r = 0 e le prime due-tre celle,
dove l'angolo centrifugo è grande (κ/(½) = 2 rad, κ/(3/2) = 0,67 rad), fanno da "nucleo finto" senza dimensione. A
ζ* = 0,95859 quel nucleo finto è esattamente risonante a energia zero: ha uno stato che decade come r^{−γ} senza
pezzo di r^{+γ} (pannello a). Sopra ζ* quello stato si stacca e scende fino a −mc², in una finestra di Z che si
stringe come a^{2γ*} = a^{0,57}: per questo Z_cr(a) → 131,36 lentamente e non si muove più.

## 5. Numeri ed esami (tutto in `lab/out/zstar_exams.json`, `lab/out/zstar_runs.json`)

### 5.1 Regola base: la formula ritrova il giro 1 (F1, retrospettivo, già visto)

| a | 3·10⁻³ | 10⁻³ | 3·10⁻⁴ | 10⁻⁴ | 3·10⁻⁵ | 10⁻⁵ | a → 0 |
|---|---|---|---|---|---|---|---|
| regola (Z_cr) | 132,783 | 132,084 | 131,712 | 131,545 | 131,453 | 131,410 | fit: **131,3612** (3 più fini del preregistrato), **131,3615** (10⁻⁴…10⁻⁵) |
| formula di raccordo | 132,908 | 132,125 | 131,724 | 131,549 | 131,454 | — | **131,3614** (ζ* esatto) |

La differenza regola − raccordo va come a (0,125 → 0,041 → 0,012 → 0,004 → 0,001): sono i termini O(a) (massa ed
energia nelle prime celle) che il raccordo trascura. F1 PASS (−0,0011 a a = 3·10⁻⁵). Il fit "Z* + c a^{2γ*} + c' a"
sulla sola regola dà 131,3612 / 131,3615 contro 131,3614 della formula: **Z* = 131,3614 ± 0,0003**.

### 5.2 Varianti cieche (previsioni nel commit 448d9c8, corse dopo)

Z_cr(a) della regola (metodo del segno, emendamento A3; dove il tracker pre-registrato funzionava i due metodi
coincidono a 10⁻⁷) contro il raccordo, e Z_inf = estrapolazione della **sola regola** (a = 10⁻³, 3·10⁻⁴, 10⁻⁴):

| variante | Z* formula | regola a = 10⁻³ / 10⁻⁴ / 3·10⁻⁵ | raccordo 10⁻³ / 10⁻⁴ | Z_inf regola | esami |
|---|---|---|---|---|---|
| R+1 | 125,8979 | 126,268 / 125,963 / 125,924 | 126,364 / 125,971 | **125,8987** | F2 ✓ F3 ✓ |
| R−i | 129,1507 | 129,693 / 129,267 / 129,203 | 129,754 / 129,272 | — | F2 ✓ |
| R+i | 95,1219 | 95,176 / 95,127 / 95,124 | 96,795 / 95,195 | — | **F2 ✗** |
| S0.25 | 103,6947 | 103,727 / 103,698 / 103,696 | 103,900 / 103,705 | **103,69466** | F2 ✓ F3 ✓ |
| K+1 (p½) | 125,8979 | 126,703 / 126,041 / 125,955 | 126,991 / 126,069 | **125,9005** | **F2 ✗** (a=10⁻³) F3 ✓ |
| K+1R+1 | 131,3614 | 132,707 / 131,697 / 131,527 | 132,832 / 131,709 | **131,3600** | F2 ✓ |

- **F2 FAIL** (pre-registrato, dichiarato): R+i a 10⁻³ (−1,62) e 10⁻⁴ (−0,068 > 0,05); K+1 a 10⁻³ (−0,29 > 0,2).
  Il motivo è lo stesso: dove 2γ* > 1 (R+i: 1,44) o l'approccio è veloce, i termini O(a) che il raccordo trascura
  dominano su quelli a^{2γ*} che il raccordo tiene. Il valore limite è comunque giusto: R+i la regola va a
  95,124 (3·10⁻⁵) contro 95,122.
- **F3 PASS**: estrapolando la sola regola, Z_inf coincide con la formula entro 0,0008 (R+1), 10⁻⁶ (S0.25),
  0,0026 (K+1).
- **F4 PASS (dualità)**: Z_inf(K+1, refl0 = +1) = 131,3600 = base 131,3612; Z_inf(K+1) = 125,9005 = R+1 125,8987.
  A passo finito sono diversi (la massa e il continuo distinguono κ = ±1), al limite coincidono come previsto.
- **U1 (universalità)**: lo spread di Z* fra le varianti è **36,2** (95,1 … 131,4): **Z* non è universale**, è un
  artefatto della condizione al centro. Previsione confermata.

### 5.3 Altri canali: anche p3/2 e d3/2 hanno un tuffo finto (K2, K3)

- **p3/2 (κ = −2)**: la regola fa tuffare il livello più basso a Z_cr = 206,192 / 206,078 / 206,038 / 206,026
  (a = 3·10⁻³ … 10⁻⁴) → **206,02** previsto. Nel continuo il p3/2 col nucleo puntiforme non si tuffa mai prima di
  ζ = 2 (Z = 274). K2 PASS (con il metodo del segno; il tracker pre-registrato aveva saltato la finestra e trovato
  il tuffo "vero" del p3/2 sopra 274: 285,5 … 278,5, emendamento A2).
- **d3/2 (κ = +2)**: uno stato di contatto (grande ~a, energia ~ (Z − Z*)/a) attraversa tutta la gap a
  Z_x = 42,075 / 42,060 / **42,055** (a = 10⁻³, 3·10⁻⁴, 10⁻⁴; energia zero), previsto **42,053**. K3 PASS.
  Attenzione: Z = 42 è il molibdeno. Per la regola radiale in onda d3/2, vicino a Z = 42,05 (in una finestra larga
  ~0,1·(a/10⁻³) in Z) c'è uno stato in più che non esiste in natura. Fuori dalla finestra non disturba
  (2γ* = 3,9: il suo effetto sugli altri livelli si spegne come a^{3,9}).
- Con le celle a (i+1)a: K−2 a 261,4 e K+2 a 249,5 (formula), cioè oltre ogni elemento.

### 5.4 Dettagli che NON contano (O1–O3, tutti PASS)

- Ordine dei passi (W = M S oppure S M): autovalori identici a 4·10⁻¹⁶ (sono matrici simili).
- Moneta spezzata (Trotter: fase, poi centrifugo, poi massa): Z_cr cambia di −0,067, −0,021, −0,006, −0,002,
  −0,0005, −0,0002 (a = 3·10⁻³ … 10⁻⁵), cioè O(a): **stesso Z***.
- Rimbalzo esterno refl1 = +1 invece di −1: ΔZ_cr = −1,6·10⁻⁸.

### 5.5 Riassunto PASS/FAIL

| esame | esito |
|---|---|
| F1 (retrospettivo), F3, F4, U1, R1, K2, K3, O1, O2, O3 | PASS |
| F2 | **FAIL** (R+i e K+1 a passo grosso: termini O(a) non nel raccordo; il limite a → 0 è giusto) |
| tracker pre-registrato su S0.25, R+i (a ≤ 3·10⁻⁴), K−2, K+2 | **FAIL dello strumento** (finestre strette): emendamenti A1–A3 |


## 6. Confronto col grafene (punto 3)

**Stessa matematica.** L'equazione radiale di Dirac 2-D senza massa del grafene ha la stessa forma della nostra con
κ → j = ±½: γ = √(j² − β²), soglia a β = ½. In entrambi i casi la soglia è dove γ diventa immaginario, e il reticolo
(o il taglio r₀) fissa la condizione al centro.

**Fisica diversa, con i numeri.**
1. *Grafene senza gap* (Wang 2013, Mao 2016): non c'è una scala di energia, la soglia resta β = ½ qualunque sia
   il taglio; r₀ cambia solo l'energia delle risonanze (Wang: r₀ da 0,5 a 1,5 nm non cambia i fit; transizione a 3
   dimeri, Z/Z_c = 1,4 ± 0,2 a ±15 %, troppo impreciso per vedere uno spostamento del 4 %).
2. *Grafene con gap / atomo con nucleo esteso*: un taglio "regolare" (potenziale piatto dentro r₀, o nucleo di
   raggio R) **alza** la soglia: Pereira–Kotov–Castro Neto g̃_c > ½; per l'atomo Z_cr(R) = 137 (1 + π²/(2 ln²R))
   circa, 177,5 a 12 fm (SUPERCRITICO.md).
3. *La nostra regola* (s = ½): la soglia scende a ζ* = 0,9586, il 4,1 % **sotto**. Non è un taglio regolare: le
   prime celle fanno un'**attrazione di contatto** che a ζ* è risonante (nel linguaggio di Kaplan et al. 2009 il
   valore nudo dell'accoppiamento di contatto attraversa il punto fisso instabile). In grafene l'analogo non è il
   taglio ma un difetto risonante a corta distanza: la vacanza ha davvero un modo a energia zero (Mao 2016), ma lì
   viene dalla rottura dei sottoreticoli, non dalla carica, e non sposta la soglia ½.
4. *L'analogo radiale del grafene nella nostra regola*: con |κ| = ½ (stessa regola, s = ½, refl0 = −1) la formula
   **non trova nessun ζ* sotto ½**: con quella forma il reticolo non abbasserebbe la soglia. (Non è il reticolo a
   nido d'ape vero: non simulato.)

Risposta: **stessa matematica (condizione al centro fissata dal reticolo), meccanismo diverso** (contatto
risonante, non taglio), segno dello spostamento opposto. Un confronto numerico diretto con il grafene non è
possibile: lì l'impurezza e il reticolo sono diversi e le misure sono al 15 %.


## 7. Rimedio (punto 4)

**Basta spostare le celle.** Con le celle a r_i = (i + s)a la formula dà ζ* < 1 solo per s < s_c ≈ 0,618 (con
refl0 = −1; con s = ½ invece **nessuna** fase del rimbalzo salva: Z* va da 95 a 131 ma c'è sempre). Per s ≥ s_c lo
stato di contatto non c'è e l'1s si tuffa **sopra** 137, verso 137 come un nucleo di raggio ~a:

| celle | regola a = 3·10⁻³ / 10⁻³ / 3·10⁻⁴ / 10⁻⁴ / 3·10⁻⁵ | raccordo 10⁻⁴ | raccordo 10⁻⁶ / 10⁻⁹ |
|---|---|---|---|
| (i + 1)a ("nucleo su una cella") | 156,785 / 151,530 / 147,779 / 145,456 / — | 145,4565 | 140,80 / 138,69 |
| (i + √3/2)a (primi vicini del cubo 3-D) | 153,230 / 149,111 / 146,133 / 144,263 / 142,795 | 144,2629 | 140,41 / 138,57 |
| (i + ¾)a | 148,918 / 145,974 / 143,858 / 142,528 / 141,478 | 142,528 | 139,74 / 138,34 |

R1 PASS: sempre sopra 137,036, in discesa, e la formula di raccordo coincide con la regola a 2·10⁻⁴.

Cosa vuol dire, onestamente:
- **Funziona senza fake**: non entra il raggio del nucleo, entra solo dove si mettono le celle. Ma è una **scelta
  a mano** (debito D1): la regola 3-D (`zz_core.Walk`, griglia (i − n/2 + ½)a) mette il nucleo al centro di un
  cubetto, con i primi vicini a a√3/2: è questo l'argomento per s = √3/2, ma la riduzione radiale della regola 3-D
  **non** è RadialWalk, quindi non è una deduzione.
- **137 non si raggiunge a passo finito**: si arriva da sopra come 1 + π²/(2 ln²(·)), cioè lentissimo (138,6 a
  a = 10⁻⁹). Questo però è il comportamento **giusto** (lo stesso di un nucleo vero che si rimpicciolisce); il
  valore fisico vero (Z_cr ≈ 170–177) richiede il nucleo esteso (SUPERCRITICO, debito del raggio).
- Il prezzo: per s = 1 restano tuffi finti in p3/2 a Z = 261 e in d3/2 a 249 (oltre ogni elemento reale).
- Una scelta della fase del rimbalzo "dedotta dalla regola 3-D" non esiste: con s = ½ tutte le fasi danno un Z*.


## 8. Validazione, scoperta, limite

- **Validazione**: β(Z) della regola uguale a NUCLEO_1S (altro thread, altro metodo); il quadro "estensione
  autoaggiunta scelta dal reticolo" (Case 1950, Kaplan et al. 2009) spiega tutto il comportamento.
- **Spiegazione controllata (nuova per noi, piccola)**: una formula senza parametri (prodotto di matrici SU(1,1)
  delle prime celle) che dà Z* = 131,3614 e, raccordata con Dirac, tutta la curva Z_cr(a); previsioni cieche per 11
  varianti tornate entro 10⁻⁶–3·10⁻³ al limite a → 0. Non ho trovato in letteratura lo stesso calcolo per un walk
  di Dirac; è nuovo come analisi, non come fisica.
- **Limiti**: (1) con le celle a (i+½)a e nucleo puntiforme la regola ha tuffi finti: s₁/₂ e p₁/₂ a 131,36 e
  125,90, p3/2 a 206,02, **d3/2 a 42,05** (dentro la tavola periodica, in una finestra stretta); (2) il rimedio
  porta a 137 solo nel limite a → 0, logaritmicamente; (3) nessuna forma chiusa per ζ*; (4) il raccordo a passo
  grosso sbaglia dove contano i termini O(a) (F2 FAIL).


## 9. Debiti e parametri

1. **D1 (solo se si adotta il rimedio)**: la posizione delle celle s (= 1 o √3/2) è una scelta geometrica a mano.
   Nello studio è una variante, non un ingresso.
2. **D2**: refl0 (rimbalzo) e s = ½ della regola congelata sono anche loro scelte a mano già presenti: questo
   studio mostra che fissano Z*.
3. Il raccordo usa Dirac del continuo fuori dalle prime celle come **strumento di analisi** (non entra nella
   regola). Il fit ρ ∈ [40, 10⁵] con 3 correzioni 1/ρ è una scelta numerica (stabile a 10⁻⁵ su β).
4. Solo forma radiale (simmetria sferica imposta); la regola 3-D a cubo non è studiata qui.
5. Nucleo puntiforme, nessuna QED, nessuno schermo.

**Parametri liberi: 0** (nella regola e nella formula). I fit di estrapolazione Z_inf + c₁a^{2γ*} + c₂a hanno 2
coefficienti descrittivi con l'esponente fissato dalla teoria. Emendamenti (dichiarati, `lab/out/zstar_amend.json`):
A1 tracker a passo fine (fallito per tempo, ucciso a 10 min), A2 griglia fine per K±2, A3 metodo del segno (validato
contro il tracker a 10⁻⁷ prima di usarlo). Non fatto: a = 3·10⁻⁵ per S1 e a = 10⁻⁵ per le varianti (limite dei 10
minuti); s_c calcolato solo a 1,5·10⁻⁵ con il criterio ζ* < 0,999 (0,6181: vicino al rapporto aureo 0,6180, ma non
ho controllato se è un caso — **non** lo dico).


## 10. Compito secondario: il film del tuffo esportato a fotogrammi

Fatto: `lab/zstar_film.py` → `lab/out/zstar_film.json` (**2,9 MB**, committato) + 6 anteprime
`/mnt/project-files/problemi-aperti/zstar_film_1…6.png` (t = 0, 10, 30, 60, 120, 240 λ_C/c).

Su origin/master l'unico export di fotogrammi è quello dell'atomo accoppiato (`lab/atomo_export_leggi.py`: cartella
con manifest.json + file float16 + fotogrammi.csv, 3-D). Il film qui è 1-D e chiesto in JSON, quindi ho tenuto le
stesse parole (manifest, campi con "forma", scalari per fotogramma) in un solo file:

    manifest: Z, nucleo, regola, unità, carica, fotogrammi (301), dt_fotogramma (0,8 λ_C/c = 400 tick), griglia,
              campi {rho, j, Q: forma [301, 384]}, risonanza {E_r, Γ, Γ dal tempo}, controllo_continuita_max
    r_bordi [385], r_centri [384]   (passo 0,025 fino a r = 2, poi 0,125 fino a r = 40; λ_C)
    t [301]
    rho[k][b]  densità radiale di probabilità dell'elettrone (per λ_C, già integrata sull'angolo)
    j[k][e]    corrente radiale attraverso la sfera r_bordi[e+1] (probabilità per λ_C/c), media esatta sul
               fotogramma (flusso del reticolo |χ+_i|² − |χ−_{i+1}|² per tick); j[0] = 0
    Q[k][e]    probabilità dentro r_bordi[e+1]
    fotogrammi_scalari: t, P_in_rw (peso dentro r_w = 0,699), Q_entro_40, passato_r30 (carica passata a r = 30)

Valori a 4 cifre significative. Continuità controllata: Q(t_k) − Q(t_{k−1}) + j·Δt = 0 entro 1,4·10⁻¹³ (prima di
arrotondare). Lettura in carica (teoria delle buche): numeri dell'elettrone × (−e); quello che esce è il positrone
(+e), quello che resta il vuoto carico (−e). Stessa corsa del giro 1 (Z = 200, guscio 12 fm, a = 0,002): Γ dal
tempo 0,0284 (finestra pre-registrata del giro 1) contro 0,0302 dalla fase, come in SUPERCRITICO §5.5.


## 11. Come rifare le corse (1 core, ognuna < 10 min)

    python lab/zstar_core.py          # controlli: beta(Z) = NUCLEO_1S, D(zeta) cambia segno a 0,9586
    python lab/zstar_prereg.py        # formula + previsioni (7 min) -> lab/out/zstar_prereg.json (commit 448d9c8)
    lab/zstar_drive.sh                # corse pre-registrate (~25 min in tutto, ogni chiamata < 10 min)
    python lab/zstar_amend.py; lab/zstar_drive3.sh    # emendamento A3 (metodo del segno, ~60 min a pezzi)
    python lab/zstar_film.py          # film (35 s)
    python lab/zstar_exams.py; python lab/zstar_figure.py
