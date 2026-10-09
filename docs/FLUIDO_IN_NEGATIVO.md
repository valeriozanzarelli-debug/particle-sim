# Fluido in negativo

**L'idea (di Valerio).** Le palline sono i "pixel" della materia: sono tante, si spingono forte e *sono* lo spazio.
Il vuoto è un fluido pieno di palline. Una particella è un posto dove le palline **mancano** o girano in modo
diverso: un buco, un vortice, un difetto. Per disegnare il mondo basta disegnare i difetti, non le palline
("in negativo"), e la GPU lavora poco. La massa dovrebbe venire dalle palline: quante ne sposti (o metti in
rotazione) per la velocità del suono al quadrato, cioè E = m c² con c = velocità del suono del fluido.

Abbiamo fatto quattro prove (F1–F4). I criteri sono stati scritti **prima** delle corse in
`lab/out/fluid_prereg.json` (commit `5bea426`). Le prove di sviluppo (dev) viste prima sono dichiarate lì.

Figura: `lab/out/fluid_figure.png`. Numeri: `lab/out/fluid_*.json`. Codice: `lab/fluid_*.py`.

## Riassunto in una tabella

| prova | cosa | esito pre-registrato |
|---|---|---|
| F1-E1 | la massa del buco sta nella "forbice" [s = 2c/R, s = c/(2R)] attorno a m·N_mancanti | **PASS** (3/3 taglie) |
| F1-E2 | massa a s* = c/R entro ±25 % di m·N_mancanti | **FAIL** (0.90, 0.80, **0.74**) |
| F1-E3 | la massa scala come N_mancanti (pendenza 0.85–1.15) | **PASS** (0.91) |
| F1-E4 | non dipende dalla massa "nuda" del buco (< 10 %) | **PASS** (2 %) |
| F2-E2a | energia del vortice ∝ ln R con pendenza π | **PASS** (0.9993 π) |
| F2-E2b | costante del nucleo come in letteratura (2 %) | **PASS** (2.078 vs 2.070) |
| F2-E2c | energia del vortice che gira = (palline mancanti)·c² + costante | **PASS** (deriva 0.008 su π) |
| F3-E3a | contrazione di Lorentz: larghezza·γ costante entro 5 % fino a 0.8 c | **PASS** (peggiore −1.5 %) |
| F3-E3b | E = γ E₀ entro 5 % fino a 0.8 c | **PASS** (peggiore +1.2 %) |
| F3-E3c | velocità limite = c **e** formula di McLaughlin–Scott entro 2 % | **FAIL** (v < c sempre sì, ma −4.4 % a F = 0.2) |
| F3-E3d | previsione: kink + antikink **non** si annichilano in fononi | previsione **confermata** (cioè: niente "e⁺e⁻ → luce") |
| F4-E4a | disegno in negativo: primitive / palline < 1 % | **FAIL** (2.9 %: 168 segmenti per 5844 palline) |
| F4-E4b | il disegno ritrova il buco (centro < 1σ, conteggio entro 20 %) | **PASS** (0.07σ, −0.9 %) |

## F1 — Quanto pesa un buco?

**Esperimento.** 2D, palline morbide repulsive (WCA, massa 1, diametro 1), densità 0.75, temperatura 1:
un liquido denso, non un gas. Il "buco" è un disco libero di raggio R_w che respinge le palline con una parete
liscia (solo forza radiale). Gli diamo un calcio (+dV e −dV dallo stesso stato, poi sottraiamo: così il rumore
termico si cancella in gran parte) e guardiamo quanto rallenta. Velocità del suono misurata (via termodinamica):
**c_s = 6.21** (isoterma 4.82).

**Palline mancanti.** Dal profilo di densità attorno al buco: R_w = 3, 6, 10 → **N_mancanti = 30.4, 102.3, 268.2**.

**Cosa ci si aspetta.** In un fluido ideale incomprimibile un cilindro in 2D si porta dietro una massa aggiunta
ρπR², cioè **esattamente la massa delle palline che mancano** (Lamb, *Hydrodynamics*; per la sfera in 3D è ½ρV).
Quindi un buco con massa nuda zero dovrebbe pesare m·N_mancanti. È l'idea di Valerio, in forma precisa.

**Cosa è successo.** La velocità del disco cala sempre, senza un plateau. In 2D il fluido viscoso continua a
trascinare sempre più palline (è il "paradosso di Stokes" 2D): **non esiste una massa unica**, esiste una massa che
dipende dalla scala di tempo. Lo abbiamo visto già nel dev e l'abbiamo dichiarato. Misuriamo allora la massa
M(s) a ogni "frequenza" s (trasformata di Laplace della risposta):

| R_w | N_mancanti | M a s=c/(2R) | M a s*=c/R | M a s=2c/R | rapporto a s* |
|---|---|---|---|---|---|
| 3  | 30.4  | 46.6  | 27.2 ± 0.8 | 14.7  | 0.90 |
| 6  | 102.3 | 128.8 | 82.3 ± 1.2 | 48.3  | 0.80 |
| 10 | 268.2 | 297.4 | 197.9 ± 3.5 | 119.4 | 0.74 |

- La massa delle palline mancanti sta **sempre dentro la forbice** (E1 pass).
- Scala con il numero di palline mancanti: pendenza **0.91** (E3 pass).
- Con massa nuda 4 volte più piccola cambia del **2 %** (E4 pass): l'inerzia è del fluido, non del disco.
- Ma al punto pre-registrato s* = c/R il rapporto scende con la taglia: 0.90 → 0.80 → **0.74** (E2 **fail**, per poco).

**Perché scende (post-hoc, dichiarato, zero parametri).** A s = c/R il fluido non è più incomprimibile: il buco
manda via onde sonore. Per un cilindro in un fluido ideale *comprimibile* la massa vale
ρπR²·K₁(x)/[x(K₀(x)+K₁(x)/x)] con x = sR/c; a x = 1 fa **0.588**. Più il buco è grande, meno conta la viscosità
(che alza la massa), e il rapporto scende verso 0.59. Nella figura (pannello a) la curva teorica a zero parametri
sta sotto i dati e i dati ci si avvicinano al crescere di R; ad alte frequenze le curve coincidono. Quindi: il
risultato è coerente con l'idrodinamica, ma **la nostra scelta s* = c/R era sbagliata come punto di test**; il
plateau "massa = palline mancanti" esiste solo per R molto grandi (finestra viscosa/comprimibile larga
~ (ν/(cR))^(1/3) ≈ 20–25 % alle nostre taglie).

**In parole semplici:** un buco nel fluido ha davvero un'inerzia, e quell'inerzia è circa il numero di palline che
mancano. Ma "circa" va preso sul serio: dipende da quanto in fretta lo spingi. Una particella vera ha **una** massa,
qualunque sia il tempo. Un buco in un liquido classico no.

## F2 — Quanto pesa un vortice (palline che girano)?

Gross–Pitaevskii 2D (un superfluido = palline nel limite quantistico), unità ħ = m = g n = 1, quindi ξ = 1,
c_s = 1. Vortice al centro, energia dentro un raggio R = 10…160.

- E = π ln(A R/ξ) con pendenza **0.9993 π** e **A = 2.078**. Letteratura (Fetter, Rev. Mod. Phys. 81, 647, 2009):
  1.464 con ξ = ħ/√(2mgn), che nelle nostre unità è 1.464·√2 = **2.070**. Scarto 0.4 %. (E2a, E2b pass)
- **Il risultato più vicino all'idea di Valerio:** anche le palline mancanti crescono come π ln R (attorno al
  vortice la densità cala come 1 − ξ²/(2r²)). E vale **energia del giro = m·c²·N_mancanti − 3.98**, con la
  differenza che *non* cresce con R (pendenza 0.008, cioè 0.2 % di π). È la legge di Bernoulli: dove il fluido gira,
  la pressione e quindi la densità cala, e (palline mancanti)×c² = energia cinetica del giro. È una **E = m c²
  vera, ma di origine idrodinamica**, con c = velocità del suono. (E2c pass)
- Attenzione: questa "massa" diverge con ln R (dipende dalla grandezza del sistema). Un elettrone non ha massa che
  dipende dalla grandezza della stanza.

**Massa di un vortice: la letteratura NON è d'accordo** (citazioni da memoria, da ricontrollare sui testi):
- Popov (1973), Duan (Phys. Rev. B 49, 12381, 1994): massa "relativistica" M ≈ E/c², che **diverge** come ln R.
- Baym & Chandler (J. Low Temp. Phys. 50, 57, 1983): massa legata alla comprimibilità / al nucleo, piccola.
- Thouless & Anglin (Phys. Rev. Lett. 99, 105301, 2007): a basse frequenze la massa efficace è **finita**,
  dell'ordine del nucleo.
- Kopnin (superfluidi fermionici / superconduttori): massa grande dagli stati legati nel nucleo.
- Simula (Phys. Rev. A 97, 023609, 2018) e lavori successivi su vortici con nucleo riempito (Richaud et al.):
  massa dipendente dal contesto.
Il punto comune: come in F1, **la massa dipende dalla scala di tempo / frequenza** con cui misuri. Noi **non**
abbiamo misurato la massa dinamica del vortice (debito).

## F3 — Relatività dal fluido (catena di Frenkel–Kontorova)

Catena di palline unite da molle su un substrato periodico (in unità sine-Gordon: φ̈ₙ = (φₙ₊₁−2φₙ+φₙ₋₁)/h² − sin φₙ).
c = 1 è la velocità del suono della catena. h = distanza tra palline / larghezza del kink (h = 0.25: ~4 palline per
larghezza). Un kink è un posto dove la catena "scivola" di un passo: una palline in più (o in meno) rispetto al
substrato. Un difetto.

**Contrazione di Lorentz ed E = γ m c²** (h = 0.25):

| v₀ | v misurata | larghezza·γ (rel.) | E/(γE₀) (rel.) |
|---|---|---|---|
| 0.2 | 0.1995 | −0.02 % | +0.01 % |
| 0.4 | 0.3988 | −0.09 % | +0.04 % |
| 0.6 | 0.5975 | −0.31 % | +0.18 % |
| 0.8 | 0.7936 | −1.45 % | +1.22 % |
| 0.9 | 0.8693 | −4.3 % | +12 % (qui si rompe) |

Fino a 0.8c **pass** (E3a, E3b). Con palline più rade si rompe prima: h = 0.5 regge fino a ~0.6c; con h = 1 (una
pallina per larghezza) il kink irradia e **si ferma a v ≈ 0.40** anche se lo lanci a 0.9. **La relatività emerge
solo se il difetto è molto più grande dei pixel.** Fatto noto per le dislocazioni nei cristalli (Frank 1949,
Eshelby 1949: contrazione "relativistica" delle dislocazioni a vite con c = velocità del suono trasversale).

**Velocità limite.** Spinto da una forza F con attrito α = 0.05, il kink non supera mai c: a F = 0.8 va a 0.993
(0.9989 nella corsa lunga). Ma il criterio pre-registrato chiedeva anche l'accordo al 2 % con la formula del
continuo (McLaughlin–Scott 1978): **FAIL** a F = 0.2 (−4.4 %) e F = 0.4 (−4.1 %). Ipotesi "non era arrivato a
regime": testata post-hoc con tempo 4× → stessi numeri, **ipotesi respinta**. È un effetto del reticolo
(probabile aggancio a risonanze con i fononi del reticolo, attrito per irraggiamento Cherenkov). Non verificato.

*Bug dichiarato:* nella prima corsa il kink spinto partiva vicino al bordo sinistro e ci andava a sbattere
(velocità senza senso, 1.44). Corretto la posizione di partenza (al centro) e il tracciamento; la prima corsa è in
`lab/out/fluid_f3_run1_bug.json`. Il criterio non è stato toccato.

**Kink + antikink → luce? No.** (previsione pre-registrata, confermata)
- h = 0.25: a tutte le velocità (0.05…0.8) si **attraversano** e ripartono: il sine-Gordon continuo è integrabile,
  le "particelle" non muoiono.
- h = 1 (catena molto discreta), v ≤ 0.4: si **legano** in un "breather" (una coppia che oscilla), che perde energia
  in fononi lentissimamente: tasso ~5·10⁻⁵ per unità di tempo (vita ~20 000). Il 31–37 % dell'energia è andato in
  fononi, il resto resta nella coppia. Come un "positronio" che non decade quasi mai.
- h = 1, v = 0.8: si attraversano, 24 % in fononi.
Quindi **nessuna annichilazione completa in "luce"**. E c'è un'altra differenza: i "fotoni" qui sono fononi con
gap (massa 1 per via del substrato), non luce senza massa.

## F4 — Disegnare in negativo

Dallo stesso istante del buco R_w = 10: invece di 5844 palline disegniamo solo la linea dove la densità
(media su 1.5 diametri) scende a metà. Risultato: **un solo anello, 168 segmenti** (×35 di compressione,
2.9 % → criterio < 1 % **FAIL**). Ma il disegno ritrova il buco benissimo: centro sbagliato di 0.07 diametri,
palline mancanti contate dal campo 265.7 contro 268.2 (−0.9 %) → **PASS**.
Vortice: 40 000 celle di griglia → 36 segmenti (×1111). Kink: 1600 palline → 1 marcatore (×1600).
Nota onesta: la compressione dipende da come si conta; un anello si potrebbe descrivere con 5–10 numeri
(post-hoc, non testato). E il costo vero non sparisce: per sapere *dove* sono i difetti devi comunque simulare
le palline (o avere una teoria efficace dei difetti, che è proprio la fisica delle particelle).

## Conteggio dei parametri e debiti

- F1: fisici 2 (densità 0.75, temperatura 1) + la taglia (variabile). Numerici: dt, skin, dV, attese, finestra +6
  per contare le mancanti. Massa nuda: formula fissata prima.
- F2: **0** parametri liberi.
- F3: 1 fisico (h, la discretezza) + α, F solo per E3c.
- F4: 3 scelte di disegno (larghezza 1.5, griglia 0.5, livello ½).

Debiti: (1) massa dinamica del vortice non misurata; (2) la finestra s di F1 era mal scelta (lo sapevamo a
metà dal dev, pre-registrato comunque); (3) la curva comprimibile di F1 è post-hoc; (4) la deviazione da
McLaughlin–Scott in F3 non è spiegata; (5) una sola temperatura/densità per il liquido; (6) tutto in 1D/2D.

## Cosa servirebbe per particelle vere

1. **Una massa sola, non una massa per ogni scala di tempo.** F1 e F2 danno inerzie che dipendono dalla frequenza
   e (in 2D) dalla grandezza del sistema. Nella relatività del kink (F3) invece la massa è unica: lì il difetto è
   *topologico* e il fluido è "rigido" (onde di tipo Klein-Gordon). Quella è la strada: difetti topologici in un
   mezzo con simmetria di Lorentz delle piccole onde, non buchi in un liquido viscoso.
2. **Spin ½.** Un difetto classico di palline è, per default, un bosone. Ci sono vie note: Finkelstein–Rubinstein
   (1968) — un solitone può comportarsi da fermione se lo spazio delle configurazioni ha un "giro" che
   cambia segno (π₁ = Z₂); è il meccanismo degli Skyrmioni (nucleoni come solitoni). Nel libro di Volovik i fermioni
   ci sono già nel fluido (³He è fatto di fermioni): gli spinori emergono vicino ai "punti di Fermi", non dalle
   palline classiche. Nessuno dei nostri difetti ha spin.
3. **Carica.** I vortici hanno una carica topologica (il numero di giri) e in 2D si attraggono/respingono come
   cariche 2D (energia ∝ ln r). In 3D però le linee di vortice non sono punti, e non viene una legge 1/r² pulita
   con una carica che si conserva come quella elettrica. Servirebbe un campo di gauge emergente (in ³He-A c'è,
   secondo Volovik, ma non è l'elettromagnetismo con le sue costanti).
4. **Perché la c del suono dovrebbe essere la c della luce.** Qui c'è un solo tipo di onda e la sua velocità fa da
   "c" per tutti i difetti di *quella* catena. Ma in un fluido vero ci sono suono longitudinale, trasversale,
   secondo suono, ognuno con la sua velocità. Nella natura **tutte** le particelle vedono la stessa c con
   precisione enorme (test di Lorentz ~10⁻²⁰ e oltre). In un fluido questo richiede una "sintonia fine" o un punto
   fisso che lo imponga (Barceló–Liberati–Visser, *Analogue Gravity*, Living Rev. Rel.; Volovik discute che
   l'universalità emerge solo a basse energie vicino a un punto di Fermi).
5. **I problemi aperti di Volovik** (dal libro del 2003, riassunti onestamente): la gravità che emerge nel fluido
   **non è Einstein** (la metrica efficace non obbedisce alle equazioni giuste); la simmetria di Lorentz si rompe ad
   alta energia (alla "scala Planck" del fluido, cioè alla scala delle palline: lo abbiamo visto in F3 con h = 1);
   il problema della costante cosmologica (Volovik propone che un vuoto in equilibrio non gravita, una proposta
   interessante ma non dimostrata per il nostro universo); campi di gauge non abeliani e chiralità non escono
   gratis.
6. **Kelvin ci aveva già provato.** Gli "atomi-vortice" (Kelvin 1867, "On vortex atoms") volevano gli atomi come
   nodi di vortice nell'etere. Non funzionò (stabilità, spettri), ma diede origine alla teoria dei nodi.
   **Dirac** (1930): un *buco* nel mare di elettroni si comporta come una particella (il positrone) — la parte
   "in negativo" dell'idea di Valerio ha un precedente illustre e corretto.

## Verdetto e consiglio

- **Funziona:** un difetto in un mezzo di palline ha inerzia ≈ palline mancanti (F1, a meno di ~20–25 %), e nel
  superfluido energia del giro = palline mancanti × c² esattamente (F2, Bernoulli). In una catena, il difetto
  obbedisce a Lorentz ed E = γmc² con c = velocità del suono fino a 0.8 c (F3), purché sia largo almeno ~4 palline.
- **Non funziona (ancora):** la massa del buco in un liquido non è un numero unico; i difetti non si annichilano
  in "luce"; niente spin, niente carica vera; la compressione del disegno è ×35, non ×100.
- **Consiglio:** se si continua, abbandonare il "buco in un liquido classico" e passare a **difetti topologici in un
  mezzo elastico/superfluido 3D** (kink → vortici-anello → Skyrmioni), con un solo tipo di onda veloce, e misurare
  per prima cosa se la loro massa è *una sola* a tutte le frequenze. Il rendering in negativo è un buon trucco
  grafico, ma non risparmia il calcolo delle palline: risparmia solo il disegno.
