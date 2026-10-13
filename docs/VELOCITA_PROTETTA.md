# Velocità protetta: si può costruire una regola in cui i loop non separano la c dell'elettrone da quella della luce?

Studio `protetta_` — 13/10/2026 (giro 2 dei problemi aperti). Regola `lab/zz_core.py` congelata, non toccata.
**Parametri liberi: 0.** Pre-registrazione: `lab/protetta_prereg.py` → `lab/out/protetta_prereg.json` (commit cdb7d59,
**prima** delle corse). Codice: `lab/protetta_libero.py` (livello libero), `lab/protetta_vestito.py` +
`lab/protetta_drive.sh` (livello vestito), `lab/protetta_analisi.py` (esami e figura). Il motore del loop è quello di
Debiti (`vuoto_f3_core.py`, origin/master 35afcc5), copiato **in sola lettura** in `lab/protetta_vf3.py`,
`lab/protetta_vf2.py`, `lab/protetta_vcore.py` (cambiano solo gli import). Figura: `lab/out/protetta_figura.png`.

## In una riga

Nessuno dei meccanismi noti (supersimmetria, rotazione tempo↔spazio del reticolo euclideo, flusso verso Lorentz
nell'infrarosso, settori fortemente accoppiati, fotone indotto) funziona per la nostra regola così com'è: o non c'è
nella regola, o con α guadagna un fattore ~0,7 dove ne serve 10⁻¹⁶. L'idea "luce e materia sono la stessa pallina"
invece **regge a metà, e in modo misurabile**: (1) la luce della regola **è già esattamente** una coppia di palline
di Weyl con metà impulso ciascuna (identità fra matrici, 10⁻¹⁵); (2) se i costituenti sono palline **con carica**,
la c vestita che si portano dietro è la stessa dell'elettrone entro **≤ 2,5·10⁻⁴ α** (in pratica zero entro la griglia;
contro i 0,045 α del gap con la luce nuda); quindi al livello delle autoenergie il fotone composto **eredita** la c
dell'elettrone e il gap sparisce. **Ma non è risolto:** manca il termine di scambio fra i due costituenti (stesso
ordine α, nessuna simmetria lo azzera), servono costituenti carichi e senza massa (che non vediamo), e il requisito
sperimentale è 3·10⁻¹⁸ α, cioè 14 ordini sotto la precisione di questo conto. Se i costituenti sono neutri, il gap
resta tutto (+0,045 α).

---

## 1. Il problema, semplice

Lo studio `lvloop_` (docs/LORENTZ_LOOP.md §10) e la fase 3 di Debiti (VUOTO_POLARIZZA.md §14 su master) hanno trovato
che l'elettrone **vestito** della regola corre con una velocità limite diversa da quella della luce **nuda**:
+0,045 α (taglio K12/3), −0,067 α (taglio K12), +0,044 α con l'ordine dei passi a giro. Questa differenza **non cala**
quando il reticolo diventa fine (p = 0). È l'effetto di Collins, Perez, Sudarsky, Urrutia, Vucetich (PRL 93, 191301
(2004), verificato online): un taglio che sceglie un sistema di riferimento porta, tramite i loop, una violazione di
Lorentz di ordine α giù fino alle basse energie. Gli esperimenti vogliono |c_e − c_γ|/c sotto ~10⁻²⁰ (elettroni del
Granchio, fotoni che non decadono negli sciami) e 5·10⁻¹⁵ (LEP, il limite più robusto). Diviso per α: serve un gap
sotto **~3·10⁻¹⁸ α**. La regola dà ~5·10⁻² α. Mancano 16 ordini.

La domanda di questo studio: **esiste un meccanismo che leghi la c dell'elettrone a quella della luce, così che i loop
non possano separarle?** In altre parole: una "simmetria custode" della velocità limite.

Nota sul coordinamento: la c della luce **vestita** della regola la calcola Debiti nella fase 4 (non duplicata qui).
Tutti i gap "contro la luce" di questo documento sono quindi **contro la luce nuda** (riserva già scritta in
LORENTZ_LOOP §10). Il fotone **composto** (coppia di walker) è un oggetto diverso, ed è il soggetto di questo studio.

## 2. Cosa è noto (ricerca)

Fonti verificate online in questa sessione dove indicato; le altre sono **a memoria** e lo dico.

| meccanismo | idea | funziona in natura/letteratura? | funziona per la nostra regola? |
|---|---|---|---|
| **Supersimmetria** (Groot Nibbelink–Pospelov PRL 94, 081601 (2005); Polchinski CQG 29, 088001 (2012), a memoria) | la SUSY vieta gli operatori di dimensione 4 che violano Lorentz; dopo la rottura resta ~(m_soft/M)² | sì in principio (unica protezione nota secondo Polchinski) | **no**: la regola non ha partner bosonici dei walker; la SUSY esatta su reticolo esiste solo per una sottoalgebra nilpotente (Catterall–Kaplan–Ünsal, Phys. Rep. 484, 71 (2009), a memoria), e per un automa a tempo discreto non c'è nemmeno H da scrivere come {Q, Q†}. Nessuna costruzione nota. |
| **Rotazione tempo↔spazio del reticolo euclideo** (simmetria ipercubica H(4) della QCD su reticolo) | la rotazione di 90° nel piano (t, x) manda c₀₀ in c_xx: l'unico tensore simmetrico invariante è δ_μν, quindi tutte le particelle hanno la stessa c nel limite del continuo, senza sintonia | sì: per questo la QCD su reticolo **isotropo** non deve sintonizzare nessuna "velocità della luce" | **no**: la regola è un automa **in tempo reale** (passi unitari): tick e cella non sono legati da nessuna simmetria. È il caso estremo del **reticolo anisotropo**, dove in QCD la velocità dei quark e quella dei gluoni vanno sintonizzate **separatamente** (Klassen, Nucl. Phys. B 533, 557 (1998); Edwards–Joó–Lin, PRD 78, 054501 (2008), a memoria). Collins, nella pratica della QCD su reticolo, è proprio questo. |
| **Discretezza senza sistema privilegiato** (insiemi causali: Bombelli–Henson–Sorkin, Mod. Phys. Lett. A 24, 2579 (2009), verificato che esiste) | una distribuzione casuale (Poisson) di punti nello spazio-tempo non sceglie nessun riferimento: niente Collins | sì per costruzione, ma la dinamica è non locale | **fuori strada**: niente reticolo cubico, niente "corri a c e ruota". Si riapre solo se lo chiede Valerio. |
| **Flusso nell'infrarosso verso Lorentz** (Chadha–Nielsen, Nucl. Phys. B 217, 125 (1983), a memoria; Anber–Donoghue, arXiv 1102.0789 (PRD 83, 105027 (2011)), verificato) | le velocità di specie che interagiscono convergono a una sola andando verso le basse energie | sì, ma **logaritmicamente**: Anber e Donoghue scrivono "le differenze normalmente svaniscono solo logaritmicamente" e stimano che per avere differenze del 10 % la scala di partenza dovrebbe essere "oltre 10^(10¹³) GeV" | **no, quantificato in §2.1**: da Planck all'elettrone si guadagna un fattore ~0,7, ne serve 10⁻¹⁶ |
| **Settore fortemente accoppiato** (Bednik–Pujolàs–Sibiryakov, arXiv 1305.0011 (JHEP 2013), verificato) | se il flusso va a un punto fisso fortemente accoppiato, la violazione cala come una **potenza** del rapporto fra le scale | sì nei modelli olografici (N grande) | **no**: la nostra regola è debolmente accoppiata (α). Servirebbe un esponente ≥ 0,72 su 22 decadi (§2.1) |
| **Fotone indotto** (Zel'dovich/Sakharov; fotone come stato legato, Bjorken, Ann. Phys. 24, 174 (1963); Volovik, *The Universe in a Helium Droplet* (2003); a memoria) | se il termine cinetico della luce nasce tutto dal loop dei fermioni, il cono della luce è quello dei fermioni al **logaritmo dominante** | eredità al log dominante | **no guadagno** (§2.2): la violazione sta nei termini di soglia (senza log), che restano O(α) |
| **Fotone composto da walker** (Bisio–D'Ariano–Perinotti, "Quantum cellular automaton theory of light", arXiv 1407.6928 (Ann. Phys. 368, 177 (2016), rivista a memoria); verificato: F^μ(k) = ∫ f_k(q) φᵀ(k/2−q) σ^μ ψ(k/2+q), ω(k) = 2|n_{k/2}|, cioè 2 volte la frequenza del Weyl a k/2; f concentrata a q ≪ k) | la luce è una coppia di fermioni senza massa: la sua c è **quella dei costituenti** | vero cinematicamente | **testato qui** (§3–4) |
| **Teoria del neutrino della luce** (de Broglie; Jordan, anni '30; problema di Pryce 1938, a memoria) | come sopra, con costituenti neutri | storicamente abbandonata (statistica e trasversalità) | caso particolare del precedente (§4.3) |
| **Covarianza di Lorentz discreta dei quantum walk** (Arrighi–Facchini–Forets, NJP 16, 093007 (2014); Bisio–D'Ariano–Perinotti, relatività deformata dell'automa (2016); a memoria) | una trasformazione del reticolo che fa da "boost" | solo 1+1 o solo libero | **non mostrata con le interazioni**: candidata generatrice (§5) |
| **Sintonia** (una velocità nuda dell'elettrone ≠ c, es. la "moneta che respira") | si sceglie il parametro nudo perché il gap vestito sia zero | funziona sempre | **sì, ma** a 1 parte su 10¹⁶–10¹⁷: è proprio il fine tuning di Collins. Non è una protezione. |

**Cosa è aperto in letteratura:** nessuno conosce un meccanismo che protegga l'universalità della velocità limite in
una teoria con taglio che sceglie un riferimento, a parte la supersimmetria (Polchinski 2012) e i settori fortemente
accoppiati con flussi a potenza. È un problema aperto vero (lo dicono gli stessi Collins et al.: "serve un meccanismo
che preservi Lorentz automaticamente", verificato sull'abstract).

### 2.1 Quanto basterebbe il flusso (Chadha–Nielsen, Anber–Donoghue): conto

Il gap δ fra due specie che interagiscono con accoppiamento α scorre come dδ/d ln μ = b α δ, con b di ordine 1/π
(Anber–Donoghue: i coefficienti sono g²/(6π²)–g²/(8π²) per i fermioni e i bosoni). Integrando dalla scala del reticolo Λ
a m_e: δ(m_e) = δ(Λ) · (m_e/Λ)^{bα}.
- Con Λ = massa di Planck: ln(Λ/m_e) = ln(2,4·10²²) = 51,5. Con b = 1 (generoso): fattore e^{−0,38} = **0,69**.
- Per portare 3·10⁻⁴ (gap della regola) a 2,5·10⁻²⁰ serve un fattore 8·10⁻¹⁷ = e^{−37}: b α ≈ 0,72, cioè **b ≈ 100**.
  Nessuna teoria debolmente accoppiata ha b così.
- Variante "a potenza" (Bednik et al.): δ ∝ (m_e/Λ)^γ. Serve γ ≥ ln(8·10⁻¹⁷)/ln(4,2·10⁻²³) = **0,72**, cioè una
  dimensione anomala di ordine 1 per l'operatore che viola Lorentz su 22 decadi. Serve un settore fortemente
  accoppiato fra Planck e l'elettrone. La regola non ce l'ha.

**Verdetto:** il flusso nell'infrarosso è reale ma, con α, dà un guadagno del 30 %. Ne serve 10¹⁶. Non basta.

### 2.2 Il fotone "indotto" eredita solo il logaritmo

Se il termine cinetico della luce è tutto generato dal loop di fermioni (Zel'dovich–Sakharov), il suo coefficiente è
(1/12π²)·ln(Λ²/m²) con il cono dei fermioni, più termini di soglia C_LV (senza log) che vengono dalla zona p ~ Λ, dove
il reticolo rompe Lorentz al 100 %. Lo spostamento relativo della c della luce è C_LV/ln(Λ²/m²). Ma nel fotone indotto
proprio 1/e² = ln(Λ²/m²)/(12π²), cioè ln(Λ²/m²) = 3π/α (per una specie). Quindi lo spostamento è C_LV·α/(3π): **lo
stesso ordine α** del caso con il fotone elementare. L'eredità vale solo al log dominante; la violazione vive nei
termini di soglia, che è esattamente dove vive Collins. Nessun guadagno.

## 3. L'idea candidata: luce e materia sono la stessa pallina

Valerio: "la luce e la materia sono la stessa pallina". Bisio, D'Ariano e Perinotti costruiscono il fotone come
coppia di fermioni di Weyl dell'automa. Se è così, la c della luce non è un parametro suo: è **ereditata** da quella
dei costituenti. La previsione (pre-registrata, T1–T6):

1. **Livello libero:** la luce della regola è già esattamente una coppia di walker di Weyl a impulso relativo zero.
2. **Livello vestito:** se i costituenti sono le stesse palline (con carica) dell'elettrone senza flip, il fotone
   composto prende la c vestita del costituente; resta da vedere se la c vestita dell'elettrone ad alto impulso (che è
   quella del costituente senza massa) è uguale a quella a basso impulso.
3. **Cosa non copre:** il termine di scambio fra i due costituenti (vertice), richiesto dall'invarianza di gauge.

### 3.1 Livello libero: risultato esatto (L1–L4, tutti PASS)

**L1. La luce della regola è la coppia.** Per ogni asse, exp(−i k S_j) (spin 1) = exp(−i k/2 σ_j) ⊗ exp(−i k/2 σ_j)
ristretto al tripletto. È un'identità tra matrici: vale per ogni k, per l'ordine xyz|zyx e per l'ordine a giro.
Numeri (2000 k a caso nella zona di Brillouin): spettri uguali entro **1,3·10⁻¹⁵** (regola) e 1,8·10⁻¹⁵ (a giro);
operatore uguale (con un intertwiner fisso) entro **2,3·10⁻¹⁵** e 5,9·10⁻¹⁵. In parole: la pallina della luce di
Bialynicki-Birula **è** due palline di Weyl che camminano insieme, ognuna con metà dell'impulso. Il centro di massa
fa (±1 ± 1)/2 = ±1 oppure 0 celle: esattamente le tre mosse della luce.

**L2. Quindi c_luce = c_Weyl a ogni k.** ω_luce(k) = 2 ω_Weyl(k/2) e la velocità di gruppo della luce a k è quella
del Weyl a k/2, lungo x, y, z, (1,1,0), (1,1,1), da |k|a = 10⁻³ a 0,5: differenza ≤ **3·10⁻¹⁰** (derivata numerica).
Nota: questo dice che la regola **già** soddisfa "luce = coppia di palline". Non lo avevamo scritto da nessuna parte.

**L3. Il prezzo del composto.** Se i due costituenti non sono esattamente paralleli (impulso relativo q_⊥), la coppia
ha energia k + 2q_⊥²/k: è una particella di massa 2q_⊥, con velocità 1 − 2(q_⊥/k)². Misurato a ka = 0,05: il rapporto
(1 − v)/(2(q_⊥/k)²) vale 1,0005 / 0,9996 / 0,993 per q_⊥/k = 0,01 / 0,02 / 0,05. Per non dare al fotone una massa sopra
il limite (m_γ < 10⁻¹⁸ eV, PDG, a memoria) i costituenti devono avere q_⊥ < 5·10⁻¹⁹ eV, cioè essere correlati su
distanze di ~4·10¹¹ m (qualche unità astronomica). È un vincolo forte sulla funzione f_k(q) di Bisio et al., non un
problema di Collins.

**L4. Le palline con il flip non possono fare la luce.** Se i costituenti girano (θ > 0, cioè **sono** elettroni),
la coppia ferma ha quasi-energia 2θ per tick (misurato: 2E/(2θ) = 1 entro 2·10⁻¹⁵): un "fotone" di massa 2 m_e.
Quindi la luce non può essere una coppia di **elettroni**. Può essere una coppia di palline **senza flip**.

## 4. Livello vestito: il fotone composto eredita la c vestita?

**La domanda numerica (D1–D3).** Il fotone composto viaggia alla velocità dei suoi costituenti (L1–L2). Se i
costituenti sono palline **senza flip ma con carica**, si vestono con la luce come l'elettrone. La c vestita di una
pallina senza massa è la c vestita dell'elettrone **ad alto impulso** (p ≫ m, dove la massa non conta). Quindi il gap
del fotone composto, al livello delle autoenergie, è

  gap = c_γ(composto) − c_e = δv(p ≫ m) − δv(p ≪ m) ≡ D(p/m).

Se D → 0 il composto eredita la c dell'elettrone. L'ho misurato con il motore esatto a tempo discreto di Debiti
(copiato, non modificato; riproduce il suo numero 0,048221 α a θ = 0,003125 alla sesta cifra), facendo variare
l'impulso dell'elettrone in δv = [E δE(p) − m δm]/p².

### 4.1 Risultati (D1, regola congelata, due tagli, assi x e y, 5 passi del reticolo)

| taglio, asse | D a p = 1 m, a → 0 | D a p = 3 m, a → 0 | D a p = 10 m (griglia di produzione / fine, EM1) |
|---|---|---|---|
| K12, x | −1,1·10⁻⁶ α | −1,4·10⁻⁶ α | −1,2·10⁻⁴ / **+8·10⁻⁶** α |
| K12, y | −1,3·10⁻⁶ α | −3,0·10⁻⁶ α | −1,2·10⁻⁴ / **+1,8·10⁻⁵** α |
| K12/3, x | −1,1·10⁻⁶ α | −1,4·10⁻⁶ α | −2,5·10⁻⁴ / **−3,5·10⁻⁵** α |
| K12/3, y | −1,3·10⁻⁶ α | −3,0·10⁻⁶ α | −2,5·10⁻⁴ / **−2,8·10⁻⁵** α |

(a → 0: fit quadratico in θ a p/m fisso; 5 punti per p = 1, 3 m, 3 punti per p = 10 m. La griglia fine è
nr = 128, nct = 40, nph = 64 al θ più piccolo, emendamento EM1 dichiarato dopo la corsa: l'emissione quasi collineare
ad alto p/m ha bisogno di più angoli.)

A p fisso in unità di m, D cala come θ² (= (p a)²: figura b, pendenza 2) fino a un pavimento di ~2·10⁻⁶ α, che è
la precisione della griglia. **Il c vestito dell'elettrone non dipende dall'impulso**: dove la regola è fine, una
pallina veloce e una lenta portano la stessa velocità vestita. È quello che dice la teoria efficace (l'operatore di
dimensione 4 c_μν non sa niente della massa), ora visto sulla regola. Le parti separate invece si muovono molto
(trasversa −1,05 → −0,69 α, Coulomb +1,10 → +0,74 α fra p = 0,05 m e p = 10 m): la somma no.

**Limite senza massa (D3):** θ = 3,9·10⁻⁴, p = 26, 51, 77 m: |D| ≤ 2,8·10⁻⁴ α (produzione), ≤ 5·10⁻⁴ α (griglia
fine, a p = 77 m): il segno cambia con la griglia, quindi è rumore numerico dell'angolo collineare, non un effetto.
**L'ordine a giro (D2):** D = −1,5…−2,8·10⁻⁴ α a p = 10 m (griglia di produzione), stessa cosa.

Anche l'**anisotropia** si eredita: D_x ≈ D_y, quindi il composto ha la stessa c lungo y dell'elettrone. Con il
fotone composto, la differenza **fra** c_e lungo y e c_γ lungo y va a zero come il gap isotropo. (Un'anisotropia
comune a tutto è un cambio di coordinate: nessuno la vede.)

### 4.2 Cosa vuol dire, e cosa NO

- **Sì:** al livello delle autoenergie, il fotone fatto di palline cariche senza flip **eredita** la c vestita
  dell'elettrone: gap ≤ 2,5·10⁻⁴ α (limite di griglia; il valore a griglia fine è ≤ 3,5·10⁻⁵ α a p = 10 m), contro
  0,045–0,067 α del gap contro la luce nuda. Un guadagno di almeno **200 volte**, probabilmente molto di più: il
  numero è compatibile con zero.
- **No:** (1) **il termine di scambio** fra i due costituenti (la coppia carica si scambia luce mentre viaggia) è dello
  stesso ordine α e, per l'invarianza di gauge, va **insieme** alle autoenergie (è la stessa struttura del loop di
  polarizzazione del vuoto a due loop). Non l'ho calcolato. Nel continuo, per due cariche parallele alla velocità della
  luce, la forza elettrica e quella magnetica si cancellano (fattore 1 − v₁·v₂ = 0): è un motivo per sperare che sia
  piccolo, ma sul reticolo vicino al taglio la cancellazione non è garantita da nessuna simmetria. (2) **Precisione**:
  serve 3·10⁻¹⁸ α, il conto arriva a ~10⁻⁵–10⁻⁴ α. "Compatibile con zero a 10⁻⁴" non è "zero a 10⁻¹⁸". Un'eredità
  **esatta** richiederebbe un teorema (un'identità per ogni ordine), non un numero. (3) **Fisica**: costituenti
  carichi senza massa sarebbero fermioni carichi senza massa, che dovrebbero comparire liberi (coppie prodotte a ogni
  energia; corsa di α diversa). Non li vediamo. Con il flip, la coppia ha massa 2 m_e (L4).

### 4.3 Costituenti neutri (teoria del neutrino della luce)

Se le palline del fotone sono neutre, non si vestono con la luce a ordine α: c_γ = c nuda = 1 (L1–L2) e il gap è
quello di prima, **+0,045 α** (K12/3; −0,067 α con K12). Nessun guadagno: è la barra rossa della figura.

### 4.4 Riassunto delle varianti (figura a)

| variante | gap c_γ − c_e (a → 0) | nota |
|---|---|---|
| regola, luce nuda, K12 | +0,067 α | Debiti 35afcc5 |
| regola, luce nuda, K12/3 | −0,045 α | Debiti |
| ordine a giro, luce nuda, K12/3 | −0,044 α | Debiti |
| continuo con lo stesso taglio | −0,053 α | Debiti |
| fotone composto da palline **neutre** | −0,045 α | §4.3 |
| fotone composto da palline **cariche**, solo autoenergie | **|gap| ≤ 2,5·10⁻⁴ α** (compatibile con 0) | questo studio, D1–D3 |
| fotone composto da **elettroni** (con flip) | non è luce: massa 2 m_e | L4 |
| richiesto dagli esperimenti | < 3·10⁻¹⁸ α (Granchio, indicativo); < 7·10⁻¹³ α (LEP) | §1 |

## 5. Alternativa: una scelta della regola che annulla il gap isotropo per simmetria?

**Le simmetrie della regola (elenco completo noto):**
1. traslazioni di una cella (spazio) e di un periodo (2 tick, oppure 6 con l'ordine a giro);
2. gruppo puntuale: **D₄h** (16 elementi, asse y distinto) per xyz|zyx; **T_h** (24) per l'ordine a giro
   (misurato in `lvloop_`, confermato da Debiti);
3. P (L ↔ R, k → −k); T antiunitaria con T² = −1 (Surfisti); coniugazione di carica C (usata da Debiti nella
   combinazione simmetrica dE₊ − dE₋);
4. U(1) globale della carica; U(1) di gauge quando la luce entra con la fase di Peierls;
5. **accidentali del caso libero:** la dispersione è simmetrica per tutte le permutazioni degli assi; e (nuovo,
   §3.1) la luce è esattamente il prodotto tensoriale di due Weyl a k/2. Né l'una né l'altra sono simmetrie della
   teoria con interazioni.

**Perché nessuna di queste può vietare il gap isotropo.** Il gap isotropo è c₀₀ − (1/3)Σ c_jj (per l'elettrone) meno
la stessa cosa per la luce: è uno **scalare** sotto qualunque rotazione dello spazio, pari sotto P, T e C. Tutte le
simmetrie 1–4 agiscono solo sullo spazio, o su spazio e tempo separatamente (T rovescia il tempo ma non lo mescola con
lo spazio). Nessuna mescola tick e cella. Quindi **nessuna** vieta il gap. Un conto di teoria dei gruppi, non numerico.

**Quale generatore bisognerebbe aggiungere.** Uno che mandi la componente tempo-tempo di un tensore in una componente
spazio-spazio:
- (a) **una rotazione tempo↔spazio** come nel reticolo euclideo. Per un automa unitario vorrebbe dire che il circuito
  letto "di lato" (lo spazio come tempo) è ancora unitario: sono i circuiti **duali-unitari** (Bertini–Kos–Prosen,
  PRL 123, 210601 (2019), a memoria). In 1+1 il walker senza massa (uno SWAP) lo è; la moneta θ lo rompe. **Ma** Debiti
  misura un gap che non va a zero per θ → 0 (+0,045 α, piatto in θ): quindi la regola con la luce **già a θ = 0** non
  ha questa simmetria. Se l'avesse, il gap sarebbe ∝ θ^p con p ≥ 1 (rottura morbida) e il passo di Planck basterebbe.
  In 3+1, con il reticolo cubico e la luce, non conosco nessun automa duale-unitario con un elettrone e un fotone.
- (b) **un boost discreto** (Arrighi–Facchini–Forets 2014 in 1+1; la "relatività deformata" dell'automa di Bisio et
  al.). Se una teoria **con interazioni** fosse covariante sotto un boost deformato che agisce allo stesso modo su
  tutte le specie, l'universalità sarebbe protetta. Nessuno l'ha mostrato con le interazioni (nemmeno per l'automa di
  Thirring, Bisio et al. 2018, a memoria).
- (c) **la supersimmetria** (vedi §2): un generatore fermionico che lega la palla di spin ½ a una palla di spin 0 o 1.
  Darebbe il legame che cerchiamo fra luce ed elettrone solo se fossero nello **stesso supermultipletto**, e non lo
  sono nemmeno in natura (il fotone sta con il fotino). Nella regola non c'è.
- (d) **"stessa pallina" come principio dinamico**, cioè che l'interazione vesta **la pallina** e non le specie:
  allora ogni particella fatta di palline erediterebbe la stessa c. Il §4 misura fin dove questo regge.

**Scelte provate o già misurate (ordine dei passi, moneta, simmetria fra assi):**
- ordine a giro: toglie l'anisotropia (T_h), il gap isotropo resta (+0,044 α);
- moneta che respira: aggiunge una velocità nuda, cioè **una manopola**: si può sintonizzare il gap a zero, ma a
  1 su 10¹⁶ (non è una protezione);
- simmetria fra assi completa (O_h): non tocca uno scalare.

Nessuna scelta dentro la strada annulla il gap isotropo **per simmetria**. Serve un generatore nuovo che mescoli tick e
cella (a, b), oppure la supersimmetria (c). Tutti e tre sono fuori dalla regola attuale.

## 6. PASS/FAIL contro la pre-registrazione

| esame | soglia | trovato | esito |
|---|---|---|---|
| L1 luce = coppia di Weyl (tripletto) | < 10⁻¹² | 1,3·10⁻¹⁵ (R), 1,8·10⁻¹⁵ (C); operatore 2,3·10⁻¹⁵ / 5,9·10⁻¹⁵ | **PASS** |
| L2 v_luce(k) = v_Weyl(k/2) | < 10⁻⁹ | 2,9·10⁻¹⁰ (R), 9·10⁻¹¹ (C) | **PASS** |
| L3 coppia non collineare: 1 − v = 2(q_⊥/k)² | entro 5 % | 1,0005 / 0,9996 / 0,993 | **PASS** |
| L4 coppia con flip: energia a riposo 2θ | entro 10⁻³ | 1 entro 2·10⁻¹⁵ | **PASS** (il composto di elettroni è massivo, come previsto) |
| D1 \|d₀\| < 10⁻³ α (ogni taglio, asse, p) | 10⁻³ α | ≤ 1,1–3·10⁻⁶ α (p = 1, 3 m); ≤ 2,5·10⁻⁴ α (p = 10 m, 3 punti) | **PASS** |
| D2 ordine a giro, \|D\| < 10⁻³ α | 10⁻³ α | ≤ 2,8·10⁻⁴ α | **PASS** |
| D3 limite senza massa, \|D\| < 2·10⁻³ α | 2·10⁻³ α | ≤ 2,8·10⁻⁴ α (produzione), ≤ 5·10⁻⁴ α (griglia fine) | **PASS** |

Emendamenti: **EM1** (dopo la corsa, non cambia nessuna soglia): controllo di griglia sui punti ad alto p/m
(`lab/protetta_griglia.py` → `lab/out/protetta_griglia.json`). Nella pre-registrazione avevo scritto il gap del
composto come −D; con la convenzione gap = c_γ − c_e è +D (solo il segno, il valore assoluto non cambia).
Fallimenti: nessun esame pre-registrato è fallito. **Il fallimento vero è di portata**: la domanda "esiste una regola
in cui i loop **non possono** separare le due c" resta senza risposta positiva (§5), e la parte che regge (§4) è
un'eredità al livello delle autoenergie, a precisione 10⁻⁴, non un teorema.

## 7. Validazione, scoperta, limite

- **Validazione:** (i) il problema e l'elenco dei rimedi sono della letteratura (Collins et al. 2004, Polchinski 2012,
  Groot Nibbelink–Pospelov 2005, Anber–Donoghue 2011, Bednik et al. 2013). (ii) ω_luce(k) = 2 ω_Weyl(k/2) è il risultato
  di Bisio, D'Ariano, Perinotti (2014/2016) per il loro automa; qui ritrovato per il nostro. (iii) Il c vestito che non
  dipende dall'impulso è ciò che la teoria efficace con operatori di dimensione 4 prevede.
- **Nuovo per la nostra regola (controllato, non una scoperta di fisica):** (a) la luce di Bialynicki-Birula della
  regola **è esattamente** la coppia simmetrica di due walker di Weyl R della regola a k/2, per ogni k e per tutti e
  due gli ordini dei passi (10⁻¹⁵): "luce e materia sono la stessa pallina" è già vero, alla lettera, al livello libero;
  (b) misurato sulla regola che la velocità vestita dell'elettrone è la stessa da p = 0,05 m a p = 77 m entro
  ≤ 5·10⁻⁴ α (≤ 3·10⁻⁶ α fino a 3 m), quindi un fotone composto di palline cariche ne erediterebbe la c (al livello
  delle autoenergie); (c) l'elenco delle simmetrie della regola e la dimostrazione (di teoria dei gruppi) che nessuna
  può vietare il gap isotropo, con i generatori che servirebbero.
- **Limite:** la regola **non ha** una protezione della velocità limite. L'unica strada interna (fotone composto)
  sposta il problema sul termine di scambio e chiede costituenti carichi senza massa, che non esistono in natura.
  Le strade che funzionerebbero (rotazione tempo↔spazio, boost discreto con interazioni, supersimmetria, settore
  fortemente accoppiato, discretezza casuale) sono tutte fuori dalla regola attuale. Senza una di queste, la regola
  come reticolo fondamentale richiede una sintonia a 1 parte su ~10¹⁶–10¹⁷ (Collins).

## 8. Debiti (cose messe a mano o non fatte)

1. **Motore del loop di Debiti** (vertice di Peierls a ogni salto, fotone = cammino della luce, Coulomb 4πα/q² dopo
   ogni moneta, peso di specie K12 = (0,6; 1,05) e K12/3, regola delle buche): tutto ereditato, con i loro debiti.
2. **Il fotone composto non è dinamico qui:** ho usato la cinematica (L1–L4) e la velocità vestita dei costituenti.
   Non ho costruito un campo bosonico composto con la funzione f_k(q) di Bisio et al. né controllato la statistica.
3. **Termine di scambio fra i costituenti: non calcolato.** È il debito principale (richiesta §9).
4. **Costituenti carichi e senza massa:** ipotesi di lavoro, in contrasto con l'osservazione (nessun fermione carico
   senza massa). Dichiarata come il prezzo dell'idea.
5. **Gap "contro la luce nuda"** (riserva di LORENTZ_LOOP): la luce vestita elementare è di Debiti fase 4.
6. **Precisione:** i gap del composto sono limitati dalla griglia angolare (EM1), ~10⁻⁵–5·10⁻⁴ α. Il requisito
   sperimentale è 3·10⁻¹⁸ α.
7. **Riferimenti a memoria** (non verificati online in questa sessione): Polchinski 2012, Groot Nibbelink–Pospelov
   (rivista), Chadha–Nielsen (contenuto), Klassen 1998, Edwards–Joó–Lin 2008, Catterall–Kaplan–Ünsal 2009,
   Bertini–Kos–Prosen 2019, Arrighi–Facchini–Forets 2014, de Broglie/Jordan/Pryce, Volovik 2003, limite sulla massa
   del fotone (PDG), rivista di Bisio et al. (Ann. Phys. 368, 177 (2016)). Verificati: Collins et al. (abstract),
   Anber–Donoghue (arXiv 1102.0789), Bednik–Pujolàs–Sibiryakov (arXiv 1305.0011), Bombelli–Henson–Sorkin
   (gr-qc/0605006), Bisio–D'Ariano–Perinotti (arXiv 1407.6928: formula di F^μ, ω = 2|n_{k/2}|, condizioni su f).
8. **b = 1 e γ** nel conto del §2.1: stime d'ordine (i coefficienti veri sono ~1/π, quindi il guadagno è ancora
   più piccolo).

**Parametri liberi: 0.** Ingressi: θ = m a, α. Manopole adattate: nessuna.

## 9. Richieste ad altri thread

1. **Debiti (fase 4, già in corso):** la c della luce **vestita** della regola (isotropa e lungo y). Con quel numero
   la riga "luce nuda" della figura diventa la riga vera, e si sa se il vestito della luce va nella direzione giusta
   (stesso segno dell'elettrone) o peggiora.
2. **Debiti (o uno studio nuovo):** il termine di **scambio dentro il fotone composto** (i due costituenti carichi
   che si scambiano la luce mentre viaggiano paralleli). È un conto a due loop, che per l'invarianza di gauge va
   insieme alle autoenergie. È il numero che decide se il fotone composto eredita la c anche oltre il livello
   delle autoenergie.
3. **Valerio:** se la strada "stessa pallina" va presa sul serio, bisogna decidere **di quali palline è fatta la
   luce**: palline senza flip e **con carica** (ma allora esisterebbero fermioni carichi senza massa, che non vediamo),
   oppure **neutre** (e allora il §4.3 dice che non si guadagna niente).

## 10. Come rifare le corse

Figura: `lab/out/protetta_figura.png` (copia in /mnt/project-files/problemi-aperti/). (a) gap per variante, scala
logaritmica, con i limiti sperimentali; (b) |D| contro il passo del reticolo: cala come θ² fino al pavimento della griglia.


```
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
python3 lab/protetta_prereg.py                # pre-registrazione (già committata)
cd lab && python3 protetta_libero.py          # L1-L4, pochi secondi
./protetta_drive.sh                           # D1-D3, ~15 min a pezzi da <= 10 min (checkpoint protetta_vestito.json)
python3 protetta_griglia.py                   # EM1, controllo di griglia (~6 min)
python3 protetta_analisi.py                   # esami, gap, figura
```
