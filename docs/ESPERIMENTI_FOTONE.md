# Catalogo degli esperimenti sul fotone

Scopo: elencare gli esperimenti **reali e storici** con cui un modello di fotone deve misurarsi, con
un numero che possa fare da criterio **PASS/FAIL** in una simulazione. È il seguito naturale di
`docs/REGOLA_FOTONE.md` (la regola a onda pilota, livello 0 di `docs/VISIONE.md`).

Stesura: 2 ottobre 2026. I riferimenti sono al lavoro originale (autori, rivista, anno). I numeri
sono quelli riportati nei lavori citati o le previsioni teoriche standard. Dove non sono sicuro del
valore misurato scrivo **"valore indicativo"**: in quel caso il criterio PASS va preso dalla
**previsione teorica** (che è esatta e calcolabile), non dal numero storico.

---

## 0. Come leggere il catalogo

Ogni scheda ha questa forma:

* **Rif.** — autori, anno, rivista.
* **Mostra** — cosa ha insegnato sul fotone, in una frase.
* **PASS** — l'osservabile e il valore (con tolleranza) da usare come criterio di simulazione.
* **Tag** — riga fissa `Gruppo | Geometria | Elementi`, usata per i conteggi del §6.
* **Note** — avvertenze; in particolare i flag:
  * **[DET]** il risultato dipende dal modello di **rivelatore** (efficienza, conteggi di buio, tempo
    morto, finestra di coincidenza) che oggi non abbiamo: oggi il clic è ideale (efficienza 1).
  * **[MAT]** serve un **modello di materiale** (indice di rifrazione, dispersione, anisotropia,
    assorbimento) che oggi non abbiamo. Se il parametro del materiale viene *dato* (n, tensore,
    potere rotatorio) l'esperimento torna di ottica lineare passiva.
  * **[SRC]** la sorgente (cristallo non lineare, atomo, punto quantico) non si simula: si modella
    come **uno stato quantistico dato** (es. coppia entangled, fotone singolo annunciato,
    luce termica). La richiesta del compito esclude le sorgenti non lineari dal gruppo L se sono
    modellabili così.
  * **[VALIDATO]** già coperto dai test esistenti (vedi §0.2).

### 0.1 Gruppi

| Gruppo | Significato |
|---|---|
| **O** | ottica di **un** fotone nel vuoto o con ottica lineare passiva (fenditure, lenti, specchi, BS, polarizzatori, lamine, interfacce dielettriche con n dato) |
| **M** | **più fotoni**: statistica quantistica, interferenza a due fotoni, entanglement, con ottica lineare (sorgente = stato dato) |
| **L** | luce–**materia**: servono elettroni/atomi (assorbimento, emissione, diffusione, vuoto quantistico) |
| **G** | **gravità / relatività** (moto della sorgente o dell'apparato, campo gravitazionale) |

### 0.2 Cosa è già validato nel programma (stato al commit `83194fa`)

| Test | File | Regola | Esperimento del catalogo |
|---|---|---|---|
| velocità c e fase 2π per λ, nessuna dispersione nel vuoto | `validation.rs` | classica | G1 (in parte) |
| elicità della polarizzazione circolare, campo trasverso | `validation.rs` | classica | O29 (in parte) |
| fascio gaussiano w(z) e Δx·Δθ = λ/4π (errore < 3%) | `validation.rs`, `pilot_validation.rs` | entrambe | O15 |
| doppia fenditura: passo λL/d (< 3%), visibilità > 0,9, corr ≥ 0,98 col flusso esatto | entrambi | entrambe | O1 |
| which-path con polarizzatori sulle fenditure: V > 0,85–0,9 / < 0,15 | entrambi | entrambe | O17 |
| un fotone alla volta costruisce le frange | `validation.rs`, `REGOLA_FOTONE.md` §6 | entrambe | O2 |
| Malus cos²θ (entro max(3σ, 1%)), tre polarizzatori 0/45/90 → 1/4 | entrambi | entrambe (3 pol.: classica) | O18 |
| CHSH: S = 2√2 ± 0,06, modello a variabili nascoste ≤ 2, nessun segnale (marginali 0,5 ± 0,02) | `validation.rs` | classica, `entangle.rs` (non locale per costruzione) | M10, M12 (versione ideale) |
| collasso istantaneo a distanza | `validation.rs` | classica | M11 (in parte) |

Tutto il resto del catalogo è **da fare**. La regola pilota oggi è solo 2-D, con oggetti piani ⊥ z
(sorgente, parete con fenditure, polarizzatore, beam-splitter polarizzatore, schermo).

### 0.3 Sigle degli elementi ottici

| Sigla | Elemento | Oggi nel simulatore? |
|---|---|---|
| `SCH` | schermo / rivelatore di posizione (clic) | sì |
| `F1` | fenditura singola | sì (parete con aperture) |
| `F2` | doppia fenditura / più fenditure | sì |
| `RET` | reticolo (molte fenditure periodiche) | quasi (parete con N aperture) |
| `OST` | ostacolo opaco (striscia, disco, bordo) | quasi (parete con parti opache) |
| `AC` | apertura circolare / disco (3-D) | no |
| `POL` | polarizzatore assorbente | sì |
| `PBS` | beam-splitter polarizzatore | sì (versione M2) |
| `BS` | beam-splitter non polarizzatore (50:50 o variabile) | no |
| `S` | specchio (riflessione totale) | no |
| `WP` | lamina d'onda (λ/2, λ/4) / rotatore di polarizzazione | no (si fa con matrici di Jones) |
| `RIT` | linea di ritardo / sfasatore (differenza di cammino regolabile) | no |
| `L` | lente (o specchio curvo) | no |
| `DI` | interfaccia dielettrica con indice n (rifrazione, riflessione parziale, TIR) | no |
| `MB` | mezzo birifrangente / anisotropo / chirale (parametro dato) | no |
| `CAV` | cavità (due specchi parzialmente riflettenti) | no |
| `EOM` | commutatore veloce (cella di Pockels) per scelte ritardate | no |
| `COINC` | conteggio di coincidenze tra più rivelatori (contatore a due/tre clic) | parziale (`entangle.rs`) |
| `SRC2` | sorgente di stato a più fotoni dato (coppia entangled, fotoni annunciati, luce termica) | parziale (`entangle.rs`: coppie in polarizzazione) |
| `MOV` | apparato in moto / rotazione / campo gravitazionale | no |
| `MAT` | modello di materia (atomi, elettroni) | no |

---

## 1. Gruppo O — un fotone, ottica lineare passiva

### O1 — Doppia fenditura (Young)
* **Rif.** T. Young, "Experiments and calculations relative to physical optics", *Phil. Trans. R. Soc.* 94, 1 (1804).
* **Mostra** La luce interferisce: il fotone passa "da entrambe le fenditure" come onda.
* **PASS** Passo delle frange Δx = λL/d entro 1% (oggi 3%); visibilità V > 0,9 per fenditure uguali e luce coerente; profilo correlato ≥ 0,98 con l'integrale di Rayleigh–Sommerfeld esatto.
* **Tag** `O | 2-D | F2, SCH`
* **Note** [VALIDATO] (classica e pilota: corr 0,993, frange rms 2,8%).

### O2 — Interferenza con luce debolissima (Taylor)
* **Rif.** G. I. Taylor, "Interference fringes with feeble light", *Proc. Camb. Phil. Soc.* 15, 114 (1909).
* **Mostra** Le frange restano identiche anche quando la luce è così debole che in media c'è meno di un fotone alla volta nell'apparato (esposizione fino a ~2000 ore).
* **PASS** Con un fotone alla volta, l'istogramma dei clic converge al profilo esatto: correlazione dentro il rumore statistico del campionatore ideale (es. 0,984 ± 0,006 a 1000 clic, 0,997 ± 0,001 a ~6000); nessuna dipendenza della forma dal flusso.
* **Tag** `O | 2-D | F2, SCH`
* **Note** [VALIDATO] (`REGOLA_FOTONE.md` §6). Taylor usò un ago (diffrazione dal bordo), non due fenditure: l'esito è lo stesso.

### O3 — Diffrazione da una fenditura (Grimaldi, Fraunhofer)
* **Rif.** F. M. Grimaldi, *Physico-mathesis de lumine* (1665); J. Fraunhofer, *Denkschr. Akad. Wiss. München* 8, 1 (1821–22).
* **Mostra** Un fascio che passa da un'apertura si allarga: Δθ ≈ λ/a (incertezza posizione–impulso).
* **PASS** Campo lontano I(θ) ∝ sinc²(πa sinθ/λ): minimi a sinθ = mλ/a entro 1%; primo massimo laterale 0,0472·I₀ entro 5%; 90,3% del flusso nel lobo centrale (fenditura lunga, 2-D).
* **Tag** `O | 2-D | F1, SCH`
* **Note** Gratis con gli elementi attuali: è il primo test da aggiungere.

### O4 — Biprisma e doppio specchio di Fresnel
* **Rif.** A. Fresnel, "Mémoire sur la diffraction de la lumière", *Ann. Chim. Phys.* (1819; memoria premiata 1818). Doppio specchio 1816.
* **Mostra** Le frange compaiono anche senza fenditure, con due immagini virtuali della stessa sorgente: l'interferenza non è un effetto dei bordi.
* **PASS** Passo Δx = λ D/d (d = distanza tra le sorgenti virtuali, D = distanza sorgenti–schermo) entro 1%; V > 0,9.
* **Tag** `O | 2-D | DI, S, SCH`
* **Note** Biprisma = due interfacce dielettriche inclinate [MAT con n dato]; il doppio specchio richiede solo `S`.

### O5 — Specchio di Lloyd
* **Rif.** H. Lloyd, "On a new case of interference of the rays of light", *Trans. R. Irish Acad.* 17, 171 (1837; letto 1834).
* **Mostra** La riflessione su un mezzo più denso cambia la fase di π.
* **PASS** Frangia **scura** sul bordo dello specchio (differenza di cammino zero); passo λD/(2h) entro 1% (h = altezza della sorgente sopra lo specchio).
* **Tag** `O | 2-D | S, SCH`
* **Note** Test diretto del segno della fase in riflessione.

### O6 — Macchia di Poisson–Arago
* **Rif.** Predizione di S. D. Poisson nella commissione per la memoria di Fresnel (1818); osservazione di F. Arago (1818–19), riportata da Fresnel, *Ann. Chim. Phys.* 11 (1819).
* **Mostra** Dietro un disco opaco, al centro dell'ombra, c'è un punto luminoso: prova decisiva della teoria ondulatoria contro quella corpuscolare di allora.
* **PASS** Onda piana su disco di raggio R a distanza z (numero di Fresnel R²/λz ≲ 10): intensità sull'asse = intensità senza disco entro 5%; profilo ≈ J₀²(kRr/z), primo zero a r = 0,383 λz/R entro 2%.
* **Tag** `O | 3-D | AC, SCH`
* **Note** In 2-D (striscia opaca lunga) c'è ancora una frangia centrale chiara, ma non vale la legge "intensità = incidente": il vero test è 3-D. Sensibile alla rotondità del bordo.

### O7 — Principio di Babinet
* **Rif.** J. Babinet, *C. R. Acad. Sci.* 4, 638 (1837).
* **Mostra** Uno schermo e il suo complementare (fenditura ↔ striscia opaca) danno la stessa figura di diffrazione fuori dal fascio geometrico.
* **PASS** Fuori dal fascio diretto, |ψ_fenditura|² = |ψ_striscia|² punto per punto entro 2%.
* **Tag** `O | 2-D | F1, OST, SCH`
* **Note** Test di linearità (somma di aperture complementari = onda libera).

### O8 — Disco di Airy e criterio di Rayleigh
* **Rif.** G. B. Airy, "On the diffraction of an object-glass with circular aperture", *Trans. Camb. Phil. Soc.* 5, 283 (1835); Lord Rayleigh, *Phil. Mag.* 8, 261 (1879).
* **Mostra** Un'apertura circolare non può focalizzare in un punto: limite di risoluzione di ogni strumento ottico.
* **PASS** Primo anello scuro a sinθ = 1,220 λ/D entro 1%; 83,8% dell'energia nel disco centrale entro 1%; picco del primo anello 1,75% del centrale. Due sorgenti incoerenti separate di 1,22 λ/D: avvallamento al centro ≈ 73,5% del massimo.
* **Tag** `O | 3-D | AC, L, SCH`
* **Note** La lente serve solo a portare il campo lontano a distanza finita; senza lente basta uno schermo molto lontano.

### O9 — Reticolo di diffrazione
* **Rif.** J. Fraunhofer, *Denkschr. Akad. Wiss. München* 8, 1 (1821–22).
* **Mostra** Molte fenditure periodiche separano le lunghezze d'onda: d sinθ = mλ.
* **PASS** Posizione degli ordini entro 0,5%; larghezza dei picchi ∝ 1/N; potere risolutivo λ/Δλ = mN entro 10%; ordini mancanti dove coincidono con i minimi della singola fenditura (es. d = 2a: ordini pari spenti, intensità < 1% del primo).
* **Tag** `O | 2-D | RET, SCH`

### O10 — Effetto Talbot
* **Rif.** H. F. Talbot, "Facts relating to optical science. No. IV", *Phil. Mag.* 9, 401 (1836); spiegazione di Lord Rayleigh, *Phil. Mag.* 11, 196 (1881).
* **Mostra** Dietro un reticolo l'onda si "auto-ricostruisce" a distanze regolari: campo vicino puramente ondulatorio.
* **PASS** Immagine del reticolo ripetuta a z_T = 2d²/λ entro 1% (paraxiale; esatta z_T = λ/(1 − √(1 − λ²/d²))); a z_T/2 immagine spostata di d/2; correlazione immagine/reticolo > 0,9 a z_T (reticolo con molti periodi).
* **Tag** `O | 2-D | RET, SCH`
* **Note** Ottimo test della regola pilota nel campo vicino (onde evanescenti trascurabili se d ≫ λ).

### O11 — Interferometro di Fabry–Pérot
* **Rif.** C. Fabry, A. Pérot, "Théorie et applications d'une nouvelle méthode de spectroscopie interférentielle", *Ann. Chim. Phys.* 16, 115 (1899).
* **Mostra** Riflessioni multiple: la luce che rimbalza molte volte interferisce con se stessa (risonanza).
* **PASS** Trasmissione T = 1/(1 + F sin²(δ/2)), F = 4R/(1−R)²; intervallo spettrale libero Δν = c/(2nL) entro 0,5%; finezza π√R/(1−R) entro 5%; T = 1 in risonanza (specchi senza perdite).
* **Tag** `O | 2-D | CAV, SCH`
* **Note** Per un fotone singolo: probabilità di uscita = T. Richiede riflessione parziale e molti passaggi: il ventaglio di Huygens di oggi non riflette all'indietro.

### O12 — Anelli di Newton
* **Rif.** I. Newton, *Opticks*, libro II (1704).
* **Mostra** Interferenza in un sottile strato d'aria tra lente e lastra: lo spessore si legge in unità di λ.
* **PASS** Raggi degli anelli scuri in riflessione r_m = √(mλR) entro 1%; centro scuro (sfasamento π).
* **Tag** `O | 3-D | DI, L, SCH`
* **Note** [MAT con n dato]. In 2-D con superficie cilindrica si ottengono frange rettilinee con la stessa legge.

### O13 — Onde stazionarie di Wiener
* **Rif.** O. Wiener, "Stehende Lichtwellen und die Schwingungsrichtung polarisirten Lichtes", *Ann. Phys.* 276, 203 (1890).
* **Mostra** Davanti a uno specchio la luce forma onde stazionarie; è il campo **elettrico** (non il magnetico) che annerisce l'emulsione.
* **PASS** Nodi di E sullo specchio e ogni λ/2 (entro 1%); ventri a λ/4 dallo specchio; con polarizzazione p a 45° di incidenza il contrasto si annulla (E dei due fasci ortogonali).
* **Tag** `O | 2-D | S, SCH`
* **Note** [DET] La lettura "E e non B" dipende dal rivelatore (accoppiamento al campo elettrico).

### O14 — Fase di Gouy
* **Rif.** L. G. Gouy, *C. R. Acad. Sci.* 110, 1251 (1890).
* **Mostra** Un'onda che passa per un fuoco guadagna una fase extra rispetto all'onda piana.
* **PASS** Fase di Gouy −arctan(z/z_R) in 3-D (totale π attraverso il fuoco), **−½ arctan(z/z_R) in 2-D** (totale π/2), entro 2%. Si misura come spostamento delle frange tra il fascio focalizzato e un'onda di riferimento.
* **Tag** `O | 2-D | L, SCH`
* **Note** Il caso 2-D è un controllo prezioso: la regola 2-D deve dare π/2, non π.

### O15 — Fascio gaussiano e principio di indeterminazione
* **Rif.** H. Kogelnik, T. Li, "Laser beams and resonators", *Appl. Opt.* 5, 1550 (1966) (formulazione standard).
* **Mostra** Un fascio di larghezza finita diverge; alla vita Δx·Δθ = λ/4π è il minimo di Heisenberg.
* **PASS** w(z) = w₀√(1 + (z/z_R)²) entro 3% (oggi 0,2–0,4%); Δx·Δθ = λ/4π entro 3%.
* **Tag** `O | 2-D | SCH`
* **Note** [VALIDATO].

### O16 — Mach–Zehnder a singolo fotone
* **Rif.** L. Zehnder, *Z. Instrumentenkd.* 11, 275 (1891); L. Mach, *Z. Instrumentenkd.* 12, 89 (1892). Con fotoni singoli: P. Grangier, G. Roger, A. Aspect, *Europhys. Lett.* 1, 173 (1986).
* **Mostra** Un solo fotone, diviso in due bracci, sceglie l'uscita secondo la differenza di fase: P = cos²(φ/2).
* **PASS** P(uscita 1) = cos²(φ/2) entro max(3σ, 1%) per φ in 8 punti; visibilità > 0,98 (Grangier 1986 riporta visibilità ≈ 98%, valore indicativo).
* **Tag** `O | 2-D | BS, S, RIT, SCH`
* **Note** È **l'esperimento base** di quasi tutto il gruppo M: va fatto subito dopo il `BS`.

### O17 — Leggi di Fresnel–Arago (which-path con polarizzazione)
* **Rif.** F. Arago, A. Fresnel, "Sur l'action que les rayons de lumière polarisée exercent les uns sur les autres", *Ann. Chim. Phys.* 10, 288 (1819).
* **Mostra** Due fasci polarizzati ortogonalmente non interferiscono; tornano a interferire se si proietta su una polarizzazione comune (cancellino).
* **PASS** Polarizzatori paralleli sulle fenditure V > 0,85; ortogonali V < 0,15; ortogonali + analizzatore a 45° davanti allo schermo: V torna > 0,85 (frange metà intensità).
* **Tag** `O | 2-D | F2, POL, SCH`
* **Note** [VALIDATO] i primi due casi; il terzo (il cancellino) **manca** ed è facile.

### O18 — Legge di Malus e tre polarizzatori
* **Rif.** É.-L. Malus, *Mém. Soc. Arcueil* 2, 143 (1809).
* **Mostra** La probabilità che un fotone attraversi un polarizzatore è cos²θ; un polarizzatore in mezzo a due incrociati fa passare 1/4.
* **PASS** cos²θ entro max(3σ, 1%); 0°/45°/90° → 0,25 ± 0,01; senza il polarizzatore centrale → 0.
* **Tag** `O | 2-D | POL, PBS, SCH`
* **Note** [VALIDATO].

### O19 — Effetto Zeno quantistico con la polarizzazione
* **Rif.** P. Kwiat, H. Weinfurter, T. Herzog, A. Zeilinger, M. A. Kasevich, "Interaction-free measurement", *Phys. Rev. Lett.* 74, 4763 (1995) (parte Zeno); idea: B. Misra, E. C. G. Sudarshan, *J. Math. Phys.* 18, 756 (1977).
* **Mostra** Misurare spesso "congela" l'evoluzione: N rotazioni di π/2N, ciascuna seguita da un polarizzatore H, lasciano il fotone quasi sempre H.
* **PASS** P(sopravvive) = cos^{2N}(π/2N): N = 1 → 0; N = 2 → 0,25; N = 6 → 0,660; N = 20 → 0,884; entro max(3σ, 1%).
* **Tag** `O | 2-D | WP, POL, SCH`
* **Note** Con i soli `POL` attuali si può già fare (polarizzatori ruotati in sequenza = stessa matematica).

### O20 — Misura senza interazione (bomba di Elitzur–Vaidman)
* **Rif.** A. C. Elitzur, L. Vaidman, *Found. Phys.* 23, 987 (1993); P. Kwiat et al., *Phys. Rev. Lett.* 74, 4763 (1995); P. G. Kwiat et al., "High-efficiency quantum interrogation measurements via the quantum Zeno effect", *Phys. Rev. Lett.* 83, 4725 (1999).
* **Mostra** Si può sapere che c'è un oggetto assorbente in un braccio senza che il fotone lo tocchi.
* **PASS** Mach–Zehnder 50:50 bilanciato con assorbitore in un braccio: P(assorbito) = 1/2, P(uscita "buia") = 1/4, P(uscita "chiara") = 1/4 entro max(3σ, 1%); efficienza η = P(buia)/(P(buia)+P(assorbito)) = 1/3. Versione Zeno (Kwiat 1999): η fino a 73% misurato (teoria → 1 per N → ∞).
* **Tag** `O | 2-D | BS, S, OST, SCH`
* **Note** **Decisivo per la regola pilota**: la pallina-fotone che passa nel braccio libero deve "sapere" che l'altro braccio è chiuso. Nella regola attuale lo sa perché le palline-onda del braccio chiuso muoiono (§7.3 di `REGOLA_FOTONE.md`): è il caso in cui la località delle palline-onda conta davvero.

### O21 — Coefficienti di Fresnel (riflessione parziale)
* **Rif.** A. Fresnel, "Mémoire sur la loi des modifications que la réflexion imprime à la lumière polarisée" (letta 1823), *Mém. Acad. Sci.* 11, 393 (1832).
* **Mostra** All'interfaccia tra due mezzi il fotone viene riflesso o trasmesso con probabilità date da n e dall'angolo; per un fotone singolo è un lancio di moneta con queste probabilità.
* **PASS** Incidenza normale aria/vetro n = 1,5: R = ((n−1)/(n+1))² = 0,040 entro 1% relativo; R_s(θ), R_p(θ) di Fresnel su 10 angoli entro 1% assoluto; R + T = 1 entro 0,1%.
* **Tag** `O | 2-D | DI, SCH`
* **Note** [MAT con n dato]. Elemento chiave: serve a O4, O12, O22–O28.

### O22 — Angolo di Brewster (polarizzazione per riflessione)
* **Rif.** É.-L. Malus, *Nouv. Bull. Sci. Soc. Philomat.* 1, 266 (1809) (scoperta, 1808); D. Brewster, *Phil. Trans. R. Soc.* 105, 125 (1815).
* **Mostra** A un angolo preciso la luce riflessa è tutta polarizzata s.
* **PASS** tanθ_B = n₂/n₁ (n = 1,5 → 56,31°) entro 0,2°; R_p(θ_B) < 10⁻³; R_s(θ_B) = 0,148 entro 2%.
* **Tag** `O | 2-D | DI, SCH`
* **Note** [MAT con n dato].

### O23 — Rifrazione (legge di Snell)
* **Rif.** W. Snell (1621, inedito); R. Descartes, *La Dioptrique* (1637).
* **Mostra** Il fotone cambia direzione entrando in un mezzo: n₁ sinθ₁ = n₂ sinθ₂ (in un mezzo la fase avanza n volte più in fretta).
* **PASS** Direzione del fascio trasmesso entro 0,1° su 10 angoli; nel mezzo λ' = λ/n entro 0,5% (passo delle frange di un interferometro con un braccio nel mezzo).
* **Tag** `O | 2-D | DI, SCH`
* **Note** [MAT con n dato].

### O24 — Riflessione totale interna
* **Rif.** J. Kepler, *Dioptrice* (1611) (osservazione); trattazione quantitativa: Fresnel (1823, vedi O21).
* **Mostra** Oltre l'angolo critico il fotone è sempre riflesso, ma l'onda "sporge" oltre la superficie (onda evanescente).
* **PASS** θ_c = arcsin(n₂/n₁) (1,5 → aria: 41,81°) entro 0,2°; R = 1 entro 10⁻³ oltre θ_c; profondità di penetrazione 1/κ, κ = k√(n₁² sin²θ − n₂²), entro 2%.
* **Tag** `O | 2-D | DI, SCH`
* **Note** [MAT con n dato]. La regola pilota usa oggi solo il semipiano in avanti: la riflessione richiede ventagli all'indietro.

### O25 — Riflessione totale frustrata (tunnel ottico)
* **Rif.** I. Newton, *Opticks* (1704) (osservazione qualitativa); E. E. Hall, "The penetration of totally reflected light into the rarer medium", *Phys. Rev.* (Series I) 15, 73 (1902); J. C. Bose (1897, con microonde).
* **Mostra** Avvicinando un secondo prisma a meno di ~λ, il fotone "salta" il vuoto: l'analogo ottico dell'effetto tunnel.
* **PASS** Trasmissione T(d) dalla formula esatta a due interfacce; per d ≫ 1/κ, T ∝ exp(−2κd) con κ come in O24 entro 5% sulla pendenza di ln T.
* **Tag** `O | 2-D | DI, SCH`
* **Note** [MAT con n dato]. Ottimo test dei campi vicini/evanescenti, che la regola pilota oggi taglia (Θ = π/2 − 10⁻³, niente onde evanescenti).

### O26 — Spostamento di Goos–Hänchen
* **Rif.** F. Goos, H. Hänchen, "Ein neuer und fundamentaler Versuch zur Totalreflexion", *Ann. Phys.* 436, 333 (1947).
* **Mostra** Un fascio totalmente riflesso esce spostato lungo la superficie, come se penetrasse un po' oltre.
* **PASS** Spostamento laterale D = −∂φ_r/∂(k sinθ) (formula di Artmann, *Ann. Phys.* 1948) entro 10%; ordine di grandezza: qualche λ vicino a θ_c (valore indicativo), diverso per s e p.
* **Tag** `O | 2-D | DI, SCH`
* **Note** [MAT con n dato]. Goos e Hänchen moltiplicarono l'effetto con decine di riflessioni.

### O27 — Effetto Hall di spin della luce (Imbert–Fedorov)
* **Rif.** F. I. Fedorov, *Dokl. Akad. Nauk SSSR* 105, 465 (1955); C. Imbert, *Phys. Rev. D* 5, 787 (1972); O. Hosten, P. Kwiat, "Observation of the spin Hall effect of light via weak measurements", *Science* 319, 787 (2008).
* **Mostra** All'interfaccia un fotone circolare destro e uno sinistro si spostano in direzioni opposte, **fuori** dal piano di incidenza: accoppiamento spin–orbita del fotone.
* **PASS** Segno opposto per elicità ±1 e valore della teoria di Bliokh–Bliokh entro 10%; ordine di grandezza sotto la lunghezza d'onda (Hosten–Kwiat: risoluzione ~1 Å con misura debole; valore indicativo).
* **Tag** `O | 3-D | DI, WP, SCH`
* **Note** [MAT con n dato]. Serve la 3-D (lo spostamento è perpendicolare al piano).

### O28 — Rombo di Fresnel
* **Rif.** A. Fresnel (1817), in *Œuvres complètes*, vol. 1.
* **Mostra** Due riflessioni totali danno uno sfasamento tra s e p: si trasforma luce lineare in circolare senza materiali birifrangenti.
* **PASS** n = 1,51, θ = 54,6°: sfasamento 45° ± 0,5° per riflessione, 90° per due; uscita circolare con grado di polarizzazione circolare > 0,99.
* **Tag** `O | 2-D | DI, SCH`
* **Note** [MAT con n dato]. Bella verifica della fase di Fresnel nella TIR.

### O29 — Birifrangenza (doppia rifrazione)
* **Rif.** E. Bartholin, *Experimenta crystalli Islandici disdiaclastici* (1669); C. Huygens, *Traité de la lumière* (1690).
* **Mostra** In un cristallo anisotropo un fotone si divide in due fasci polarizzati ortogonalmente (ordinario e straordinario); un fotone singolo esce da uno solo.
* **PASS** Calcite a 589 nm: n_o = 1,658, n_e = 1,486 (dati); angolo di separazione (walk-off) dal calcolo con il tensore entro 1%; probabilità per fascio = cos²/sin² della polarizzazione di ingresso entro max(3σ, 1%).
* **Tag** `O | 3-D | MB, SCH`
* **Note** [MAT con tensore dato]. In 2-D si può fare con l'asse ottico nel piano. Il `PBS` attuale ne è la versione semplificata.

### O30 — Lamine d'onda e tomografia di polarizzazione (Stokes)
* **Rif.** G. G. Stokes, *Trans. Camb. Phil. Soc.* 9, 399 (1852); D. F. V. James, P. G. Kwiat, W. J. Munro, A. G. White, "Measurement of qubits", *Phys. Rev. A* 64, 052312 (2001).
* **Mostra** Lo stato di polarizzazione di un fotone è un punto sulla sfera di Poincaré; con λ/4, λ/2 e polarizzatore si ricostruisce.
* **PASS** Per 6 stati di prova (H, V, D, A, R, L) e stati arbitrari: parametri di Stokes ricostruiti entro 3σ; fedeltà > 0,99.
* **Tag** `O | 2-D | WP, POL, SCH`
* **Note** Strumento di misura indispensabile per il gruppo M.

### O31 — Attività ottica (potere rotatorio)
* **Rif.** F. Arago (1811); J.-B. Biot, *Mém. Acad. Sci.* (1812–1815).
* **Mostra** In un mezzo chirale (quarzo, zucchero) il piano di polarizzazione ruota: le due elicità viaggiano a velocità diverse.
* **PASS** Quarzo a 589 nm: 21,7°/mm (dato); rotazione ∝ spessore entro 1%; rotazione = differenza di fase tra elicità/2.
* **Tag** `L | 2-D | MB, POL, SCH`
* **Note** [MAT]. Con potere rotatorio dato diventa ottica lineare (rotatore di Jones).

### O32 — Fase geometrica (Pancharatnam–Berry)
* **Rif.** S. Pancharatnam, *Proc. Indian Acad. Sci. A* 44, 247 (1956); A. Tomita, R. Y. Chiao, "Observation of Berry's topological phase by use of an optical fiber", *Phys. Rev. Lett.* 57, 937 (1986).
* **Mostra** Un fotone il cui stato di polarizzazione (o direzione) percorre un circuito chiuso acquista una fase uguale a metà dell'angolo solido racchiuso (o all'angolo solido nello spazio delle direzioni).
* **PASS** Sequenza di lamine λ/4–λ/2(α)–λ/4 in un Mach–Zehnder: fase geometrica = 2α entro 1°; fibra elicoidale: rotazione della polarizzazione = angolo solido entro 5% (valore indicativo per l'esperimento originale).
* **Tag** `O | 3-D | WP, BS, S, SCH`
* **Note** La versione con lamine è 2-D; quella con la fibra è 3-D.

### O33 — Momento angolare di spin della luce (Beth)
* **Rif.** R. A. Beth, "Mechanical detection and measurement of the angular momentum of light", *Phys. Rev.* 50, 115 (1936).
* **Mostra** Ogni fotone circolare porta momento angolare ±ħ: una lamina λ/2 che ne inverte l'elicità riceve 2ħ per fotone.
* **PASS** Coppia sulla lamina = 2 P/ω (P = potenza) con segno corretto; accordo del lavoro originale entro ~10% (valore indicativo).
* **Tag** `O | 3-D | WP, MOV`
* **Note** Richiede il bilancio di momento angolare tra fotone ed elemento ottico (oggi gli elementi sono fissi e non ricevono nulla).

### O34 — Momento angolare orbitale (fasci a vortice)
* **Rif.** L. Allen, M. W. Beijersbergen, R. J. C. Spreeuw, J. P. Woerdman, *Phys. Rev. A* 45, 8185 (1992); J. Leach, M. J. Padgett, S. M. Barnett, S. Franke-Arnold, J. Courtial, "Measuring the orbital angular momentum of a single photon", *Phys. Rev. Lett.* 88, 257901 (2002).
* **Mostra** Un fotone con fase e^{iℓφ} porta momento angolare ℓħ indipendente dallo spin; un fotone singolo si smista in base a ℓ.
* **PASS** Fascio di Laguerre–Gauss: zero di intensità sull'asse, fase che gira di 2πℓ; interferenza con onda piana = spirale a ℓ bracci; smistatore di Leach (Mach–Zehnder con prismi di Dove a 90°) separa ℓ pari/dispari con contrasto > 0,9 (valore indicativo).
* **Tag** `O | 3-D | L, BS, S, SCH`
* **Note** Impossibile in 2-D.

### O35 — Traiettorie medie dei fotoni (misura debole)
* **Rif.** S. Kocsis, B. Braverman, S. Ravets, M. J. Stevens, R. P. Mirin, L. K. Shalm, A. M. Steinberg, "Observing the average trajectories of single photons in a two-slit interferometer", *Science* 332, 1170 (2011).
* **Mostra** Le "traiettorie medie" dei fotoni dietro una doppia fenditura, ricostruite con misure deboli della direzione, coincidono con le traiettorie di Bohm e **non attraversano** l'asse di simmetria.
* **PASS** Campo di velocità trasversale misurato = j_x/j_z della teoria (cioè la legge di guida della regola pilota) entro 10% sui piani misurati; nessuna traiettoria attraversa l'asse centrale.
* **Tag** `O | 2-D | F2, L, MB, SCH`
* **Note** **Il test più diretto della regola pilota**: le traiettorie della pallina-fotone si possono confrontare con quelle pubblicate. La misura debole (lamina di calcite sottile) è [MAT con tensore dato]; il confronto con la teoria si può fare anche senza simularla. [SRC] fotoni da punto quantico.

### O36 — Scelta ritardata di Wheeler con un fotone singolo
* **Rif.** J. A. Wheeler (1978); V. Jacques, E. Wu, F. Grosshans, F. Treussart, P. Grangier, A. Aspect, J.-F. Roch, "Experimental realization of Wheeler's delayed-choice gedanken experiment", *Science* 315, 966 (2007).
* **Mostra** Decidere *dopo* che il fotone è entrato nell'interferometro se chiuderlo (frange) o lasciarlo aperto (quale via) dà sempre la statistica della configurazione finale: il fotone non "decide" in anticipo se essere onda o particella.
* **PASS** Configurazione chiusa: P = cos²(φ/2) con V > 0,9 (originale ≈ 0,94, valore indicativo); aperta: P = 1/2 indipendente da φ (V < 0,05) e via riconoscibile > 0,99; il risultato non dipende dall'istante della scelta (anche quando la scelta avviene a fotone già nei bracci).
* **Tag** `O | 2-D | BS, PBS, EOM, RIT, SCH`
* **Note** **Decisivo**. In un modello a onda pilota il risultato è corretto purché la seconda divisione agisca sull'onda già presente in entrambi i bracci. [SRC] centro NV.

### O37 — Complementarità quantitativa V² + D² ≤ 1
* **Rif.** B.-G. Englert, *Phys. Rev. Lett.* 77, 2154 (1996); V. Jacques et al., "Delayed-choice test of quantum complementarity with interfering single photons", *Phys. Rev. Lett.* 100, 220402 (2008).
* **Mostra** Più si sa del cammino, meno frange: la somma dei quadrati di visibilità e distinguibilità non supera 1, e vale 1 per stati puri.
* **PASS** Variando la riflettività del secondo BS: V² + D² = 1 entro 0,03 per tutti i punti.
* **Tag** `O | 2-D | BS, S, EOM, RIT, SCH`
* **Note** Generalizza O16, O17, O36.

### O38 — Contestualità con un solo fotone (KCBS)
* **Rif.** A. A. Klyachko, M. A. Can, S. Binicioğlu, A. S. Shumovsky, *Phys. Rev. Lett.* 101, 020403 (2008); R. Łapkiewicz, P. Li, C. Schaeff, N. K. Langford, S. Ramelow, M. Wieśniak, A. Zeilinger, "Experimental non-classicality of an indivisible quantum system", *Nature* 474, 490 (2011).
* **Mostra** Anche un **singolo** fotone (in tre cammini: un "qutrit") non può avere risultati prefissati indipendenti dal contesto della misura.
* **PASS** Somma KCBS ≥ −3 per ogni modello non contestuale; QM = 5 − 4√5 = −3,944; misurato −3,893 ± 0,006. PASS: valore simulato < −3,8 e entro 0,02 da −3,944 (ottica ideale).
* **Tag** `O | 2-D | BS, PBS, WP, RIT, SCH`
* **Note** Teorema di Kochen–Specker: il fotone singolo è già non classico, senza entanglement.

### O39 — Kochen–Specker indipendente dallo stato (Peres–Mermin)
* **Rif.** E. Amselem, M. Rådmark, M. Bourennane, A. Cabello, "State-independent quantum contextuality with single photons", *Phys. Rev. Lett.* 103, 160405 (2009).
* **Mostra** La contraddizione con la non contestualità vale per **qualunque** stato del fotone (polarizzazione × cammino).
* **PASS** ⟨χ_PM⟩ ≤ 4 classico, QM = 6 per ogni stato; misurato ≈ 5,4–5,5 (valore indicativo). PASS: simulato 6 ± 0,05 per 4 stati diversi.
* **Tag** `O | 2-D | BS, PBS, WP, SCH`

---

## 2. Gruppo M — più fotoni, statistica, entanglement

### M1 — Hanbury Brown–Twiss (raggruppamento dei fotoni)
* **Rif.** R. Hanbury Brown, R. Q. Twiss, "Correlation between photons in two coherent beams of light", *Nature* 177, 27 (1956); "A test of a new type of stellar interferometer on Sirius", *Nature* 178, 1046 (1956).
* **Mostra** I fotoni di una sorgente termica arrivano "a grappoli": le coincidenze a due rivelatori dietro un BS sono il doppio di quelle casuali a ritardo zero.
* **PASS** g⁽²⁾(0) = 2 per luce termica (polarizzata, monomodo, risoluzione temporale ≪ tempo di coerenza) entro 5%; g⁽²⁾(τ) → 1 per τ ≫ τ_c; g⁽²⁾ = 1 per luce laser.
* **Tag** `M | 2-D | BS, COINC, SRC2`
* **Note** [SRC] luce caotica come miscela di stati coerenti con fase casuale. [DET] Nell'originale il contrasto misurato è ridotto dalla banda dei rivelatori. Il diametro di Sirio è un valore indicativo (~6–7 millisecondi d'arco).

### M2 — Statistica dei conteggi: laser contro luce termica
* **Rif.** F. T. Arecchi, "Measurement of the statistical distribution of Gaussian and laser sources", *Phys. Rev. Lett.* 15, 912 (1965).
* **Mostra** Il numero di fotoni contati in una finestra è poissoniano per un laser e di Bose–Einstein per la luce termica.
* **PASS** Varianza ⟨Δn²⟩ = ⟨n⟩ (laser) e ⟨n⟩ + ⟨n⟩² (termica, un modo) entro 3%; test χ² della distribuzione con p > 0,01.
* **Tag** `M | 2-D | COINC, SRC2`
* **Note** [SRC] [DET] (efficienza < 1 rende più poissoniana qualunque statistica).

### M3 — Interferenza tra fasci indipendenti a livello di fotone singolo
* **Rif.** R. L. Pfleegor, L. Mandel, "Interference of independent photon beams", *Phys. Rev.* 159, 1084 (1967).
* **Mostra** Due laser indipendenti producono frange (rilevate come correlazioni tra conteggi) anche quando nell'apparato c'è circa un fotone alla volta: "ogni fotone interferisce solo con se stesso" (Dirac) non è la frase giusta.
* **PASS** Anticorrelazione tra i conteggi in due zone dello schermo distanti mezza frangia; coefficiente di correlazione con il segno e il modulo della teoria entro 20% (valore indicativo).
* **Tag** `M | 2-D | SCH, COINC, SRC2`
* **Note** [SRC] due stati coerenti con fase relativa casuale ma stabile nella finestra. [DET].

### M4 — Anticorrelazione del fotone singolo a un beam splitter
* **Rif.** P. Grangier, G. Roger, A. Aspect, "Experimental evidence for a photon anticorrelation effect on a beam splitter: a new light on single-photon interferences", *Europhys. Lett.* 1, 173 (1986). Precedente: J. F. Clauser, *Phys. Rev. D* 9, 853 (1974).
* **Mostra** Un fotone annunciato, arrivato su un BS, fa clic in **una sola** uscita: non si divide. Ma lo stesso fotone in un Mach–Zehnder interferisce (O16).
* **PASS** Parametro di anticorrelazione α = P_c/(P_t P_r) (≈ g⁽²⁾(0) condizionato): classico α ≥ 1; misurato α = 0,18 ± 0,06; ideale 0. PASS: α < 0,05 con sorgente ideale (nessuna coincidenza vera in 10⁴ fotoni annunciati).
* **Tag** `M | 2-D | BS, COINC, SRC2`
* **Note** **Decisivo.** [SRC] fotoni annunciati (cascata atomica). [DET] il valore storico dipende dai fotoni doppi della sorgente. Nella regola pilota è automatico: c'è una sola pallina-fotone.

### M5 — Interferenza di Hong–Ou–Mandel
* **Rif.** C. K. Hong, Z. Y. Ou, L. Mandel, "Measurement of subpicosecond time intervals between two photons by interference", *Phys. Rev. Lett.* 59, 2044 (1987).
* **Mostra** Due fotoni identici che arrivano insieme sulle due facce di un BS 50:50 escono **sempre dalla stessa parte**: le coincidenze scompaiono. È un effetto puramente bosonico, senza analogo per due "palline" classiche.
* **PASS** Coincidenze C(τ) con avvallamento a τ = 0: visibilità (C_max − C_min)/C_max → 1 con fotoni identici (moderni > 0,99; l'originale era inferiore, valore indicativo); larghezza dell'avvallamento = tempo di coerenza (originale ≈ 16 µm di cammino FWHM, ~50 fs); onde classiche: visibilità ≤ 0,5.
* **Tag** `M | 2-D | BS, RIT, COINC, SRC2`
* **Note** **Decisivo.** Richiede che due fotoni condividano un'**onda a due particelle simmetrica** (bosoni). Con due palline-fotone guidate ciascuna dalla propria onda indipendente il risultato è 1/2 di coincidenze: FAIL. È il prossimo vero ostacolo concettuale del modello.

### M6 — Hong–Ou–Mandel con distinguibilità e cancellino
* **Rif.** P. G. Kwiat, A. M. Steinberg, R. Y. Chiao, "Observation of a 'quantum eraser': a revival of coherence in a two-photon interference experiment", *Phys. Rev. A* 45, 7729 (1992).
* **Mostra** Se i due fotoni sono resi distinguibili (polarizzazione ruotata di θ) l'avvallamento si riduce come il loro sovrapporsi; proiettando dopo il BS su polarizzazioni a 45° l'avvallamento torna.
* **PASS** Visibilità V(θ) = cos²θ entro 0,03 (θ = 0, 30°, 45°, 60°, 90°); con polarizzatori a 45° dopo il BS la visibilità torna al valore di θ = 0.
* **Tag** `M | 2-D | BS, WP, POL, RIT, COINC, SRC2`
* **Note** Decisivo insieme a M5.

### M7 — Interferometria di Franson (tempo–energia)
* **Rif.** J. D. Franson, "Bell inequality for position and time", *Phys. Rev. Lett.* 62, 2205 (1989); P. G. Kwiat, A. M. Steinberg, R. Y. Chiao, "High-visibility interference in a Bell-inequality experiment for energy and time", *Phys. Rev. A* 47, R2472 (1993).
* **Mostra** Due fotoni, ognuno in un interferometro sbilanciato (cammino corto/lungo), mostrano frange nelle coincidenze che dipendono dalla **somma** delle fasi, anche se nessuno dei due da solo interferisce.
* **PASS** Coincidenze ∝ 1 + V cos(φ₁ + φ₂); singoli piatti (V_singoli < 0,02); V > 1/√2 = 0,707 necessaria per violare Bell; misurato 0,804 ± 0,006. PASS ideale: V = 1 ± 0,02 (dopo la selezione delle coincidenze).
* **Tag** `M | 2-D | BS, S, RIT, COINC, SRC2`
* **Note** [SRC] coppia con tempo di nascita incerto ed energie correlate.

### M8 — Coerenza indotta senza emissione indotta
* **Rif.** X. Y. Zou, L. J. Wang, L. Mandel, "Induced coherence and indistinguishability in optical interference", *Phys. Rev. Lett.* 67, 318 (1991). Applicazione: G. B. Lemos et al., "Quantum imaging with undetected photons", *Nature* 512, 409 (2014).
* **Mostra** Il segnale di due sorgenti interferisce solo se i loro fotoni "gemelli" (idler) sono indistinguibili; bloccare l'idler, che non viene mai rivelato, spegne le frange del segnale.
* **PASS** Visibilità del segnale = |T| (trasmissione **in ampiezza** dell'oggetto sul cammino dell'idler), cioè lineare in |T|, entro 0,03; idler bloccato → V < 0,02. La frequenza di conteggio del segnale non cambia.
* **Tag** `M | 2-D | BS, OST, RIT, SCH, SRC2`
* **Note** [SRC] modellabile come stato dato (|s₁i⟩ + e^{iφ}|s₂i⟩, con l'idler del primo cristallo allineato al secondo). Test molto severo di "chi porta l'onda".

### M9 — Freedman–Clauser (prima violazione di Bell)
* **Rif.** S. J. Freedman, J. F. Clauser, "Experimental test of local hidden-variable theories", *Phys. Rev. Lett.* 28, 938 (1972).
* **Mostra** Le correlazioni di polarizzazione di coppie di fotoni violano la disuguaglianza di Bell (forma CHSH/Freedman).
* **PASS** δ = |R(22,5°)/R₀ − R(67,5°)/R₀| − 1/4 ≤ 0 per i modelli locali; misurato δ = 0,050 ± 0,008, in accordo con la previsione QM per le loro efficienze (≈ 0,051; valore indicativo). Ideale con polarizzatori perfetti: δ = (√2 − 1)/4 = 0,104.
* **Tag** `M | 2-D | POL, COINC, SRC2`
* **Note** [SRC] cascata del calcio. [DET] (ipotesi di campionamento equo).

### M10 — Aspect–Grangier–Roger (CHSH con polarizzatori a due canali)
* **Rif.** A. Aspect, P. Grangier, G. Roger, "Experimental realization of Einstein–Podolsky–Rosen–Bohm Gedankenexperiment: a new violation of Bell's inequalities", *Phys. Rev. Lett.* 49, 91 (1982).
* **Mostra** Con analizzatori a due uscite (PBS) la violazione di CHSH è netta.
* **PASS** S = 2,697 ± 0,015 misurato; QM 2,70 ± 0,05 (con le aperture reali); ideale 2√2 = 2,828. PASS: S = 2,828 ± 0,03 ideale; S ≤ 2 per qualunque modello locale del simulatore.
* **Tag** `M | 2-D | PBS, WP, COINC, SRC2`
* **Note** [VALIDATO] nella versione ideale (`validation.rs`, modalità classica, `entangle.rs` non locale per costruzione). [SRC].

### M11 — Aspect–Dalibard–Roger (analizzatori che cambiano durante il volo)
* **Rif.** A. Aspect, J. Dalibard, G. Roger, "Experimental test of Bell's inequalities using time-varying analyzers", *Phys. Rev. Lett.* 49, 1804 (1982).
* **Mostra** La violazione resta anche quando l'orientazione dei polarizzatori cambia mentre i fotoni sono in volo (separazione di tipo spazio, 12 m).
* **PASS** S' = 0,101 ± 0,020 contro il limite locale 0 (QM 0,112); ideale come M9. PASS: violazione invariata quando la scelta avviene dopo l'emissione; marginali di ciascun lato indipendenti dalla scelta dell'altro entro 0,02.
* **Tag** `M | 2-D | PBS, EOM, COINC, SRC2`
* **Note** [VALIDATO] in parte (collasso istantaneo, nessun segnale).

### M12 — Weihs et al. (località stretta)
* **Rif.** G. Weihs, T. Jennewein, C. Simon, H. Weinfurter, A. Zeilinger, "Violation of Bell's inequality under strict Einstein locality conditions", *Phys. Rev. Lett.* 81, 5039 (1998).
* **Mostra** Scelte casuali veloci e rivelatori a 400 m: la violazione non si spiega con segnali tra i due lati.
* **PASS** S = 2,73 ± 0,02 (30 σ). PASS ideale come M10, con scelte generate da sorgenti casuali indipendenti per i due lati e decise dopo l'emissione.
* **Tag** `M | 2-D | PBS, EOM, COINC, SRC2`
* **Note** [SRC] coppie da conversione parametrica. [DET].

### M13 — Bell senza scappatoie con fotoni (Giustina, Shalm)
* **Rif.** M. Giustina et al., "Significant-loophole-free test of Bell's theorem with entangled photons", *Phys. Rev. Lett.* 115, 250401 (2015); L. K. Shalm et al., "Strong loophole-free test of local realism", *Phys. Rev. Lett.* 115, 250402 (2015).
* **Mostra** Chiudendo insieme la scappatoia della località e quella del rivelatore (efficienza > 2/3), il realismo locale è escluso.
* **PASS** Disuguaglianza di Eberhard/CH: J ≤ 0 per i modelli locali; p-value 3,74 × 10⁻³¹ (Giustina), 2,3 × 10⁻⁷ (Shalm). PASS nel simulatore: con un modello di rivelatore di efficienza η (≈ 0,75–0,8, valore indicativo) e stato non massimamente entangled ottimizzato, J < 0 con significatività > 5σ; con η < 2/3 nessuna violazione possibile.
* **Tag** `M | 2-D | PBS, EOM, COINC, SRC2`
* **Note** **Decisivo**, ma [DET]: richiede un modello esplicito di efficienza e conteggi di buio dei rivelatori.

### M14 — Gomma quantistica a scelta ritardata (Kim et al.)
* **Rif.** Y.-H. Kim, R. Yu, S. P. Kulik, Y. Shih, M. O. Scully, "Delayed 'choice' quantum eraser", *Phys. Rev. Lett.* 84, 1 (2000).
* **Mostra** Il fotone segnale, rivelato **prima** del gemello, mostra frange solo se si selezionano i casi in cui il gemello ha perso l'informazione sul cammino; la somma di tutti i casi non ha mai frange.
* **PASS** D₀ totale: nessuna frangia (V < 0,05). Coincidenze D₀∧D₁ e D₀∧D₂: frange con sfasamento π tra loro (frange e anti-frange), V > 0,9 ideale; D₀∧D₃, D₀∧D₄ (cammino noto): nessuna frangia. Ordine temporale dei clic irrilevante.
* **Tag** `M | 2-D | F2, BS, S, SCH, COINC, SRC2`
* **Note** **Decisivo.** [SRC] coppia nata in una delle due regioni del cristallo (stato dato). Il caso ideale ha V = 1; le visibilità storiche sono più basse (valore indicativo).

### M15 — Gomma quantistica con doppia fenditura e lamine λ/4
* **Rif.** S. P. Walborn, M. O. Terra Cunha, S. Pádua, C. H. Monken, "Double-slit quantum eraser", *Phys. Rev. A* 65, 033818 (2002).
* **Mostra** Lamine λ/4 sulle fenditure marcano il cammino del fotone segnale (frange spente); un polarizzatore sul gemello entangled le riaccende, anche se messo dopo.
* **PASS** Senza cancellino V < 0,1; con polarizzatore sul gemello a ±45°: frange e anti-frange, V > 0,9 ideale.
* **Tag** `M | 2-D | F2, WP, POL, SCH, COINC, SRC2`

### M16 — Scelta ritardata quantistica (sovrapposizione di onda e particella)
* **Rif.** A. Peruzzo, P. Shadbolt, N. Brunner, S. Popescu, J. L. O'Brien, "A quantum delayed-choice experiment", *Science* 338, 634 (2012); F. Kaiser, T. Coudreau, P. Milman, D. B. Ostrowsky, S. Tanzilli, "Entanglement-enabled delayed-choice experiment", *Science* 338, 637 (2012).
* **Mostra** Se la "scelta" (chiudere o no l'interferometro) è essa stessa un qubit in sovrapposizione, il fotone mostra un comportamento intermedio continuo tra onda e particella.
* **PASS** Con l'ancilla nello stato cosα|chiuso⟩ + sinα|aperto⟩: I(φ) = cos²α·cos²(φ/2) + sin²α·½, visibilità cos²α entro 0,03.
* **Tag** `M | 2-D | BS, PBS, WP, RIT, COINC, SRC2`

### M17 — Scambio di entanglement a scelta ritardata
* **Rif.** X.-S. Ma, S. Zotter, J. Kofler, R. Ursin, T. Jennewein, Č. Brukner, A. Zeilinger, "Experimental delayed-choice entanglement swapping", *Nature Physics* 8, 479 (2012).
* **Mostra** Se due fotoni già rivelati fossero entangled o no si decide **dopo**, con una misura su altri due fotoni.
* **PASS** Sottoinsieme "misura di Bell": correlazioni dello stato entangled (testimone di entanglement < 0); sottoinsieme "misura separabile": correlazioni nulle nelle basi coniugate; indipendente dall'ordine temporale.
* **Tag** `M | 2-D | BS, PBS, EOM, COINC, SRC2`
* **Note** Quattro fotoni.

### M18 — Teletrasporto quantistico
* **Rif.** D. Bouwmeester, J.-W. Pan, K. Mattle, M. Eibl, H. Weinfurter, A. Zeilinger, "Experimental quantum teleportation", *Nature* 390, 575 (1997).
* **Mostra** Lo stato di polarizzazione di un fotone viene trasferito a un altro fotone lontano usando entanglement e due bit classici, senza che il primo viaggi.
* **PASS** Fedeltà dello stato teletrasportato > 2/3 (limite classico); ideale 1 nei casi selezionati (con BS lineare si identifica solo |ψ⁻⟩: 25% dei casi; con l'analizzatore completo lineare al massimo 50%). Storico ≈ 0,8 (valore indicativo).
* **Tag** `M | 2-D | BS, PBS, WP, COINC, SRC2`
* **Note** Il limite del 50% per la misura di Bell con sola ottica lineare (Lütkenhaus, Calsamiglia, Suominen, *Phys. Rev. A* 59, 3295, 1999) è esso stesso un test.

### M19 — Scambio di entanglement
* **Rif.** J.-W. Pan, D. Bouwmeester, H. Weinfurter, A. Zeilinger, "Experimental entanglement swapping: entangling photons that never interacted", *Phys. Rev. Lett.* 80, 3891 (1998).
* **Mostra** Due fotoni che non si sono mai incontrati diventano entangled.
* **PASS** Condizionando su |ψ⁻⟩ al BS centrale: correlazioni dei fotoni esterni = |ψ⁻⟩, visibilità ideale 1 in tutte le basi; visibilità storica ≈ 0,65 (valore indicativo).
* **Tag** `M | 2-D | BS, PBS, COINC, SRC2`

### M20 — GHZ a tre fotoni
* **Rif.** D. M. Greenberger, M. A. Horne, A. Zeilinger (1989); J.-W. Pan, D. Bouwmeester, M. Daniell, H. Weinfurter, A. Zeilinger, "Experimental test of quantum nonlocality in three-photon Greenberger–Horne–Zeilinger entanglement", *Nature* 403, 515 (2000).
* **Mostra** Con tre fotoni il realismo locale fallisce non statisticamente ma "a colpo singolo": predice il segno opposto a quello osservato.
* **PASS** Misure XYY, YXY, YYX danno prodotto −1 ⇒ il realismo locale predice XXX = −1, QM +1 (convenzioni del lavoro: ruoli scambiati). Frazione di eventi conformi alla QM: misurata 0,85 ± 0,04; ideale 1. Mermin: M ≤ 2 locale, QM 4; misurato ≈ 2,84 (da visibilità 71%).
* **Tag** `M | 2-D | BS, PBS, WP, COINC, SRC2`

### M21 — Codifica densa
* **Rif.** K. Mattle, H. Weinfurter, P. G. Kwiat, A. Zeilinger, "Dense coding in experimental quantum communication", *Phys. Rev. Lett.* 76, 4656 (1996).
* **Mostra** Agendo su un solo fotone di una coppia entangled si trasmettono più di 1 bit.
* **PASS** Capacità del canale log₂3 = 1,58 bit con analizzatore di Bell lineare (3 classi distinguibili su 4) entro 0,05 (ideale).
* **Tag** `M | 2-D | BS, PBS, WP, COINC, SRC2`

### M22 — Paradosso di Hardy
* **Rif.** L. Hardy, *Phys. Rev. Lett.* 68, 2981 (1992) e 71, 1665 (1993); J. R. Torgerson, D. Branning, C. H. Monken, L. Mandel, *Phys. Lett. A* 204, 323 (1995); W. T. M. Irvine, J. F. Hodelin, C. Simon, D. Bouwmeester, "Realization of Hardy's thought experiment with photons", *Phys. Rev. Lett.* 95, 030401 (2005).
* **Mostra** Un evento che il realismo locale dichiara impossibile accade con probabilità finita: "non località senza disuguaglianze".
* **PASS** Tre probabilità che devono essere 0 (entro 0,005) e una che vale (5√5 − 11)/2 = 0,0902 al massimo, per lo stato ottimale non massimamente entangled; entro 0,005.
* **Tag** `M | 2-D | PBS, WP, COINC, SRC2`

### M23 — Disuguaglianza di Leggett (non località realistica)
* **Rif.** A. J. Leggett, *Found. Phys.* 33, 1469 (2003); S. Gröblacher, T. Paterek, R. Kaltenbaek, Č. Brukner, M. Żukowski, M. Aspelmeyer, A. Zeilinger, "An experimental test of non-local realism", *Nature* 446, 871 (2007); T. Paterek et al., *Phys. Rev. Lett.* 99, 210406 (2007).
* **Mostra** Esclude anche una classe di modelli **non locali** in cui ogni fotone ha una polarizzazione definita.
* **PASS** Violazione della disuguaglianza di Leggett dal valore QM esatto dell'articolo; storico: pochi σ (Gröblacher), ~80 σ senza ipotesi di simmetria (Paterek) (valori indicativi).
* **Tag** `M | 3-D | PBS, WP, COINC, SRC2`
* **Note** Richiede polarizzazioni ellittiche e direzioni fuori da un grande cerchio della sfera di Poincaré (la geometria 3-D è sulla sfera, non nello spazio: il fascio resta 2-D). Interessante per il nostro `entangle.rs`, che è un modello non locale.

### M24 — Tempo di tunnel di un fotone singolo
* **Rif.** A. M. Steinberg, P. G. Kwiat, R. Y. Chiao, "Measurement of the single-photon tunneling time", *Phys. Rev. Lett.* 71, 708 (1993).
* **Mostra** Il picco del pacchetto di un fotone che attraversa uno specchio multistrato (barriera) arriva **prima** di uno che viaggia nel vuoto (ritardo 1,47 ± 0,21 fs, "velocità" apparente ≈ 1,7 c), senza violare la causalità.
* **PASS** Anticipo misurato con un HOM tra fotone attraversante e fotone di riferimento = tempo di fase di Wigner della barriera entro 15%.
* **Tag** `M | 2-D | DI, BS, RIT, COINC, SRC2`
* **Note** [MAT con n dato] (multistrato). Richiede M5.

### M25 — Diffrazione a due fotoni (litografia quantistica)
* **Rif.** M. D'Angelo, M. V. Chekhova, Y. Shih, "Two-photon diffraction and quantum lithography", *Phys. Rev. Lett.* 87, 013602 (2001).
* **Mostra** Una coppia di fotoni entangled che attraversa insieme una doppia fenditura dà frange con passo **dimezzato** (lunghezza d'onda di de Broglie di gruppo λ/2).
* **PASS** Coincidenze tra due rivelatori nello stesso punto: passo λL/(2d) entro 2%; conteggi singoli: nessuna frangia di quel passo.
* **Tag** `M | 2-D | F2, SCH, COINC, SRC2`

### M26 — Stati NOON e lunghezza d'onda di de Broglie a più fotoni
* **Rif.** P. Walther, J.-W. Pan, M. Aspelmeyer, R. Ursin, S. Gasparoni, A. Zeilinger, "De Broglie wavelength of a non-local four-photon state", *Nature* 429, 158 (2004); M. W. Mitchell, J. S. Lundeen, A. M. Steinberg, "Super-resolving phase measurements with a multiphoton entangled state", *Nature* 429, 161 (2004).
* **Mostra** N fotoni in uno stato NOON oscillano N volte più veloce con la fase.
* **PASS** Periodo delle coincidenze a N fotoni in funzione di φ = 2π/N entro 1%; visibilità ideale 1.
* **Tag** `M | 2-D | BS, PBS, WP, RIT, COINC, SRC2`

### M27 — Immagine "fantasma" a due fotoni
* **Rif.** T. B. Pittman, Y. H. Shih, D. V. Strekalov, A. V. Sergienko, "Optical imaging by means of two-photon quantum entanglement", *Phys. Rev. A* 52, R3429 (1995).
* **Mostra** L'immagine di un oggetto attraversato da un fotone compare nelle coincidenze con un rivelatore mobile sul fotone gemello, che non ha mai toccato l'oggetto.
* **PASS** Immagine nelle coincidenze con ingrandimento dato dall'equazione di lente "a due fotoni" entro 2%; conteggi singoli uniformi.
* **Tag** `M | 2-D | F2, L, SCH, COINC, SRC2`

### M28 — Campionamento di bosoni
* **Rif.** M. A. Broome et al., "Photonic boson sampling in a tunable circuit", *Science* 339, 794 (2013); J. B. Spring et al., "Boson sampling on a photonic chip", *Science* 339, 798 (2013).
* **Mostra** N fotoni identici in una rete di BS escono con probabilità date dai **permanenti** della matrice: l'interferenza a molti bosoni.
* **PASS** Per 3 fotoni in una rete 6×6 nota: frequenze di uscita = |Perm(U_S)|² entro 3σ (distanza di variazione totale < 0,05).
* **Tag** `M | 2-D | BS, RIT, COINC, SRC2`
* **Note** Estensione di M5; costo che esplode con N (vedi limiti in `VISIONE.md`).

### M29 — Traiettorie di Bohm "surreali" e non locali
* **Rif.** D. H. Mahler, L. Rozema, K. Fisher, L. Vermeyden, K. J. Resch, H. M. Wiseman, A. Steinberg, "Experimental nonlocal and surreal Bohmian trajectories", *Science Advances* 2, e1501466 (2016).
* **Mostra** Con coppie entangled, la traiettoria media di un fotone (misura debole) dipende dalla misura fatta sull'altro lontano; le traiettorie possono sembrare "surreali" rispetto a un rivelatore di cammino.
* **PASS** Velocità trasversale media condizionata alla polarizzazione del gemello = previsione di Bohm per lo stato a due particelle entro 10% (valore indicativo).
* **Tag** `M | 2-D | F2, MB, POL, SCH, COINC, SRC2`
* **Note** Il test che direbbe se la regola pilota si estende in modo coerente alle coppie.

### M30 — Tomografia di uno stato a due fotoni
* **Rif.** D. F. V. James, P. G. Kwiat, W. J. Munro, A. G. White, *Phys. Rev. A* 64, 052312 (2001).
* **Mostra** La matrice densità 4×4 di una coppia si ricostruisce da 16 misure di coincidenza.
* **PASS** Fedeltà con lo stato preparato > 0,99; concurrence di |ψ⁻⟩ = 1 ± 0,02.
* **Tag** `M | 2-D | WP, PBS, COINC, SRC2`
* **Note** Strumento di misura per M9–M23.

---

## 3. Gruppo L — luce e materia

Tutti questi richiedono un modello di materia che oggi non esiste (livello 1 e oltre di `VISIONE.md`).
Si elencano perché definiscono la "particella" fotone (energia hν, impulso h/λ, statistica) e perché
saranno i test dei livelli successivi.

### L1 — Dispersione nel prisma (Newton)
* **Rif.** I. Newton, "A letter ... containing his new theory about light and colors", *Phil. Trans. R. Soc.* 6, 3075 (1672).
* **Mostra** La luce bianca è una miscela di colori; ogni colore ha un suo indice di rifrazione.
* **PASS** Deviazione minima del prisma per 3 lunghezze d'onda da n(λ) dato (formula di Sellmeier) entro 0,05°.
* **Tag** `L | 2-D | DI, SCH, MAT`
* **Note** [MAT]. Con n(λ) dato diventa O; dedurre n(λ) dagli atomi è un test del livello 3.

### L2 — Effetto fotoelettrico
* **Rif.** P. Lenard, *Ann. Phys.* 313, 149 (1902); A. Einstein, *Ann. Phys.* 322, 132 (1905); R. A. Millikan, "A direct photoelectric determination of Planck's h", *Phys. Rev.* 7, 355 (1916).
* **Mostra** L'energia scambiata con un elettrone viene in pacchetti hν, indipendenti dall'intensità.
* **PASS** Potenziale d'arresto eV = hν − W: pendenza h entro 0,5% (Millikan: h = 6,57 × 10⁻³⁴ J·s); nessuna soglia di intensità, nessun ritardo misurabile all'accensione.
* **Tag** `L | 3-D | MAT, SCH`

### L3 — Effetto Compton
* **Rif.** A. H. Compton, "A quantum theory of the scattering of X-rays by light elements", *Phys. Rev.* 21, 483 (1923).
* **Mostra** Il fotone urta un elettrone come una particella con impulso h/λ.
* **PASS** Δλ = (h/m_e c)(1 − cosθ), λ_C = 2,426 pm, entro 1% su 3 angoli; sezione d'urto di Klein–Nishina entro 5%.
* **Tag** `L | 3-D | MAT, SCH`

### L4 — Effetto Raman
* **Rif.** C. V. Raman, K. S. Krishnan, "A new type of secondary radiation", *Nature* 121, 501 (1928).
* **Mostra** Un fotone diffuso da una molecola cede o acquista un quanto di vibrazione.
* **PASS** Righe Stokes/anti-Stokes a ±ν_vib (es. CCl₄: 459 cm⁻¹, la più intensa) entro 1 cm⁻¹; rapporto anti-Stokes/Stokes ≈ exp(−hν_vib/kT) × fattore ν⁴ entro 10%.
* **Tag** `L | 3-D | MAT, SCH`

### L5 — Radiazione del corpo nero
* **Rif.** O. Lummer, E. Pringsheim, *Verh. Dtsch. Phys. Ges.* 1, 23 (1899); M. Planck, *Verh. Dtsch. Phys. Ges.* 2, 237 (1900).
* **Mostra** Il fotone come quanto di energia hν dà lo spettro termico (e la statistica di Bose–Einstein).
* **PASS** Picco di Wien λ_max T = 2,898 × 10⁻³ m·K entro 0,5%; legge di Stefan–Boltzmann entro 1%.
* **Tag** `L | 3-D | CAV, MAT`

### L6 — Antiraggruppamento in fluorescenza di risonanza
* **Rif.** H. J. Kimble, M. Dagenais, L. Mandel, "Photon antibunching in resonance fluorescence", *Phys. Rev. Lett.* 39, 691 (1977).
* **Mostra** Un singolo atomo emette un fotone alla volta: g⁽²⁾(0) < 1, impossibile per un'onda classica.
* **PASS** g⁽²⁾(0) < 1 (moderno con atomi/punti quantici < 0,01; valore indicativo); g⁽²⁾(τ) = 1 − e^{−3Γτ/2}(cos Ωτ + ...) per l'atomo a due livelli (formula di Mollow/Kimble) entro 5%.
* **Tag** `L | 2-D | BS, COINC, MAT`
* **Note** [DET]. Con sorgente data come fotone singolo diventa M4.

### L7 — Luce compressa (squeezing)
* **Rif.** R. E. Slusher, L. W. Hollberg, B. Yurke, J. C. Mertz, J. F. Valley, *Phys. Rev. Lett.* 55, 2409 (1985); L.-A. Wu, H. J. Kimble, J. L. Hall, H. Wu, *Phys. Rev. Lett.* 57, 2520 (1986); H. Vahlbruch, M. Mehmet, K. Danzmann, R. Schnabel, *Phys. Rev. Lett.* 117, 110801 (2016).
* **Mostra** Il rumore di una quadratura del campo scende sotto quello del vuoto (con l'altra che cresce).
* **PASS** Rumore di rivelazione omodina sotto il livello dello sparo; record 15 dB (Vahlbruch 2016). PASS con stato dato: varianza della quadratura = e^{−2r}/4 entro 2% dopo un BS di perdita nota.
* **Tag** `L | 2-D | BS, RIT, MAT`
* **Note** [SRC] il cristallo non lineare non si simula; con stato compresso dato è un test di gruppo M che richiede una rappresentazione a numero di fotoni variabile (la regola attuale ha un fotone per evento).

### L8 — Emissione spontanea (vita media)
* **Rif.** V. Weisskopf, E. Wigner, *Z. Phys.* 63, 54 (1930).
* **Mostra** Un atomo eccitato emette un fotone anche nel vuoto, con legge esponenziale.
* **PASS** Vita dello stato 2p dell'idrogeno 1,596 ns entro 1%; decadimento esponenziale (scarti < 1% fino a 5 vite).
* **Tag** `L | 3-D | MAT`

### L9 — Effetto Purcell (emissione modificata in cavità)
* **Rif.** E. M. Purcell, *Phys. Rev.* 69, 681 (1946) (abstract); P. Goy, J. M. Raimond, M. Gross, S. Haroche, *Phys. Rev. Lett.* 50, 1903 (1983); R. G. Hulet, E. S. Hilfer, D. Kleppner, *Phys. Rev. Lett.* 55, 2137 (1985).
* **Mostra** La rapidità con cui un atomo emette dipende dai modi del campo intorno: una cavità la accelera o la blocca.
* **PASS** Fattore di Purcell F_P = (3/4π²)(λ/n)³ Q/V entro 10%; emissione inibita sotto la frequenza di taglio (Hulet: vita allungata di oltre un ordine di grandezza, valore indicativo).
* **Tag** `L | 3-D | CAV, MAT`

### L10 — Spostamento di Lamb
* **Rif.** W. E. Lamb, R. C. Retherford, "Fine structure of the hydrogen atom by a microwave method", *Phys. Rev.* 72, 241 (1947).
* **Mostra** Il vuoto del campo elettromagnetico sposta i livelli atomici.
* **PASS** 2S₁/₂ − 2P₁/₂ dell'idrogeno = 1057,8 MHz entro 0,1% (1947: ≈ 1000 MHz).
* **Tag** `L | 3-D | MAT`

### L11 — Forza di Casimir
* **Rif.** H. B. G. Casimir, *Proc. K. Ned. Akad. Wet.* 51, 793 (1948); S. K. Lamoreaux, *Phys. Rev. Lett.* 78, 5 (1997); U. Mohideen, A. Roy, *Phys. Rev. Lett.* 81, 4549 (1998).
* **Mostra** Due specchi nel vuoto si attraggono perché limitano i modi del campo tra loro.
* **PASS** Pressione π²ħc/(240 d⁴) tra piani perfetti; accordo sperimentale 5% (Lamoreaux), 1% (Mohideen–Roy) con correzioni di conducibilità finita e rugosità.
* **Tag** `L | 3-D | CAV, MOV, MAT`

### L12 — Assorbimento risonante senza rinculo (Mössbauer)
* **Rif.** R. L. Mössbauer, *Z. Phys.* 151, 124 (1958).
* **Mostra** In un cristallo il rinculo va all'intero reticolo: righe gamma larghe ~10⁻¹² dell'energia (valore indicativo per ⁵⁷Fe, 14,4 keV).
* **PASS** Frazione senza rinculo dal fattore di Debye–Waller entro 5%; larghezza naturale della riga.
* **Tag** `L | 3-D | MAT`
* **Note** Strumento di G6.

### L13 — Pressione di radiazione
* **Rif.** P. N. Lebedev, *Ann. Phys.* 311, 433 (1901); E. F. Nichols, G. F. Hull, *Phys. Rev.* (Series I) 13, 307 (1901) e 17, 26 (1903).
* **Mostra** La luce porta impulso E/c: spinge.
* **PASS** Pressione = (1 + R)·I/c entro 1% (accordo storico di qualche %, valore indicativo).
* **Tag** `L | 3-D | S, MOV, MAT`

### L14 — Effetto Faraday
* **Rif.** M. Faraday, "On the magnetization of light and the illumination of magnetic lines of force", *Phil. Trans. R. Soc.* 136, 1 (1846).
* **Mostra** Un campo magnetico nella materia fa ruotare la polarizzazione (non reciproco: raddoppia tornando indietro).
* **PASS** Rotazione = V·B·L con costante di Verdet data entro 1%; raddoppio in andata e ritorno (diverso dall'attività ottica O31, che si annulla).
* **Tag** `L | 2-D | MB, S, POL, MAT`

### L15 — Bell senza scappatoie con spin elettronici (Hensen)
* **Rif.** B. Hensen et al., "Loophole-free Bell inequality violation using electron spins separated by 1.3 kilometres", *Nature* 526, 682 (2015).
* **Mostra** Due spin di centri NV a 1,3 km, resi entangled tramite fotoni (scambio di entanglement), violano CHSH senza scappatoie.
* **PASS** S = 2,42 ± 0,20 su 245 prove (p = 0,039); ideale 2√2.
* **Tag** `L | 2-D | BS, COINC, MAT`
* **Note** I fotoni servono solo come messaggeri (misura di Bell al BS, come M19).

### L16 — Fotoni indistinguibili da un punto quantico
* **Rif.** C. Santori, D. Fattal, J. Vučković, G. S. Solomon, Y. Yamamoto, "Indistinguishable photons from a single-photon device", *Nature* 419, 594 (2002).
* **Mostra** Fotoni emessi in tempi diversi dallo stesso emettitore sono identici (HOM tra fotoni successivi).
* **PASS** Sovrapposizione media ≈ 0,8 (valore indicativo); con emettitore ideale 1 (riduce a M5).
* **Tag** `L | 2-D | BS, RIT, COINC, MAT`

### L17 — Righe spettrali (Fraunhofer, Kirchhoff–Bunsen)
* **Rif.** J. Fraunhofer (1814–15); G. Kirchhoff, R. Bunsen, *Ann. Phys.* 186, 161 (1860).
* **Mostra** Ogni atomo emette e assorbe fotoni a frequenze precise: il fotone è un'impronta digitale della materia.
* **PASS** Righe D del sodio 589,0 / 589,6 nm entro 0,01 nm; serie di Balmer dell'idrogeno (656,3; 486,1; 434,0 nm) entro 0,01%.
* **Tag** `L | 2-D | RET, SCH, MAT`
* **Note** Primo test del livello 1 (atomo di idrogeno), punto 3 dell'ordine di lavoro in `VISIONE.md`.

---

## 4. Gruppo G — gravità e relatività

### G1 — Velocità della luce (Fizeau, Foucault)
* **Rif.** H. Fizeau, *C. R. Acad. Sci.* 29, 90 (1849); L. Foucault, *C. R. Acad. Sci.* 55, 501 (1862).
* **Mostra** La luce ha velocità finita (Fizeau ≈ 3,13 × 10⁸ m/s, +5%; Foucault 2,98 × 10⁸ m/s, −0,6%).
* **PASS** Tempo di volo / distanza = c entro 10⁻⁶ (nel simulatore è costruita dentro la regola).
* **Tag** `G | 2-D | S, MOV, SCH`
* **Note** [VALIDATO] in parte (velocità c di ogni pallina).

### G2 — Aberrazione stellare
* **Rif.** J. Bradley, *Phil. Trans. R. Soc.* 35, 637 (1729).
* **Mostra** La direzione apparente della luce dipende dal moto dell'osservatore: costante di aberrazione 20,5″ ≈ v/c.
* **PASS** Angolo di aberrazione = arctan(γ v/c) (formula relativistica) entro 0,1%.
* **Tag** `G | 2-D | MOV, SCH`

### G3 — Trascinamento di Fresnel nell'acqua in moto (Fizeau)
* **Rif.** H. Fizeau, *C. R. Acad. Sci.* 33, 349 (1851); A. A. Michelson, E. W. Morley, *Am. J. Sci.* 31, 377 (1886).
* **Mostra** La luce in un mezzo in moto è trascinata solo in parte: coefficiente 1 − 1/n² (poi spiegato dalla composizione relativistica delle velocità).
* **PASS** Acqua (n = 1,333): coefficiente 0,437; Michelson–Morley 1886: 0,434 ± 0,02. PASS entro 2%.
* **Tag** `G | 2-D | DI, BS, S, MOV, SCH`
* **Note** [MAT con n dato] e mezzo in moto.

### G4 — Michelson–Morley
* **Rif.** A. A. Michelson, E. W. Morley, "On the relative motion of the Earth and the luminiferous ether", *Am. J. Sci.* 34, 333 (1887).
* **Mostra** La velocità della luce non dipende dalla direzione del moto della Terra: niente etere.
* **PASS** Spostamento atteso con l'etere 0,4 frange; osservato < 0,02 (meno di 1/20, "probabilmente meno di 1/40"). PASS: spostamento delle frange ruotando l'interferometro < 10⁻³ frange.
* **Tag** `G | 2-D | BS, S, MOV, SCH`

### G5 — Doppler trasverso (Ives–Stilwell)
* **Rif.** H. E. Ives, G. R. Stilwell, "An experimental study of the rate of a moving atomic clock", *J. Opt. Soc. Am.* 28, 215 (1938).
* **Mostra** La frequenza della luce emessa da un atomo in moto è ridotta del fattore γ (dilatazione del tempo).
* **PASS** Spostamento medio delle righe avanti/indietro = (γ − 1)ν entro 10% (valore indicativo per l'originale; moderni 10⁻⁹).
* **Tag** `G | 2-D | MOV, MAT, SCH`

### G6 — Spostamento verso il rosso gravitazionale (Pound–Rebka)
* **Rif.** R. V. Pound, G. A. Rebka, "Apparent weight of photons", *Phys. Rev. Lett.* 4, 337 (1960); R. V. Pound, J. L. Snider, *Phys. Rev.* 140, B788 (1965).
* **Mostra** Un fotone che sale in un campo gravitazionale perde energia: Δν/ν = gh/c².
* **PASS** 22,5 m: gh/c² = 2,46 × 10⁻¹⁵; rapporto misura/teoria 1,05 ± 0,10 (1960), 0,9990 ± 0,0076 (1965). PASS entro 1%.
* **Tag** `G | 2-D | MOV, MAT, SCH`
* **Note** Richiede L12 per misurarlo come allora; nel simulatore basta una fase dipendente dal potenziale.

### G7 — Deflessione della luce del Sole (Eddington)
* **Rif.** F. W. Dyson, A. S. Eddington, C. Davidson, *Phil. Trans. R. Soc. A* 220, 291 (1920) (eclisse del 1919).
* **Mostra** La gravità devia la luce del doppio della previsione newtoniana.
* **PASS** Al bordo del Sole 1,75″ (RG; 0,87″ newtoniano). Misurato 1,98″ ± 0,12″ (Sobral), 1,61″ ± 0,30″ (Príncipe). Moderni (VLBI): parametro γ = 1 entro ~10⁻⁴ (valore indicativo). PASS entro 0,1%.
* **Tag** `G | 3-D | MOV`

### G8 — Ritardo di Shapiro
* **Rif.** I. I. Shapiro, "Fourth test of general relativity", *Phys. Rev. Lett.* 13, 789 (1964); B. Bertotti, L. Iess, P. Tortora, "A test of general relativity using radio links with the Cassini spacecraft", *Nature* 425, 374 (2003).
* **Mostra** La luce che passa vicino a una massa impiega più tempo (fino a ~200 µs per radar su Venere/Marte).
* **PASS** Cassini: γ − 1 = (2,1 ± 2,3) × 10⁻⁵. PASS: ritardo (1 + γ)(GM/c³) ln(...) con γ = 1 entro 10⁻³.
* **Tag** `G | 3-D | MOV`

### G9 — Effetto Sagnac e rotazione terrestre
* **Rif.** G. Sagnac, *C. R. Acad. Sci.* 157, 708 e 1410 (1913); A. A. Michelson, H. G. Gale, *Astrophys. J.* 61, 140 (1925).
* **Mostra** In un interferometro che ruota i due versi di percorrenza non sono equivalenti: Δφ = 8πAΩ/(λc).
* **PASS** Michelson–Gale: 0,230 ± 0,005 frange contro 0,236 previste. PASS: Δφ entro 1% per tre velocità angolari.
* **Tag** `G | 2-D | BS, S, MOV, SCH`

### G10 — Onde gravitazionali con un interferometro (LIGO)
* **Rif.** B. P. Abbott et al. (LIGO Scientific Collaboration and Virgo Collaboration), "Observation of gravitational waves from a binary black hole merger", *Phys. Rev. Lett.* 116, 061102 (2016).
* **Mostra** Un Michelson con bracci di 4 km misura deformazioni dello spazio di 10⁻²¹; il rumore quantistico dei fotoni (sparo, pressione di radiazione) è il limite.
* **PASS** Picco di deformazione h ≈ 1,0 × 10⁻²¹ (GW150914). Nel simulatore: sfasamento = 2kLh per un braccio che cambia lunghezza di hL, entro 1%; rumore di sparo ∝ 1/√N_fotoni.
* **Tag** `G | 2-D | BS, S, CAV, MOV, SCH`

---

## 5. Totale

39 (O) + 30 (M) + 17 (L) + 10 (G) = **96 schede**. Nota: O31 (attività ottica) è nella sezione O
perché è un elemento di polarizzazione, ma il suo tag è **L**: senza un potere rotatorio dato serve
la materia. Contando per tag: **O 38, M 30, L 18, G 10**.

---

## 6. Tabelle riassuntive

### 6.1 Per gruppo e geometria (contati dai tag)

| Gruppo | Esperimenti | di cui 2-D sufficiente | di cui 3-D necessario | già validati (anche in parte) |
|---|---|---|---|---|
| O | 38 | 30 | 8 | 5 (O1, O2, O15, O17, O18) |
| M | 30 | 29 | 1 | 2 (M10 ideale, M11 in parte) |
| L | 18 | 8 | 10 | 0 |
| G | 10 | 8 | 2 | 1 in parte (G1) |
| **Totale** | **96** | **75** | **21** | **8** |

Gli esperimenti 3-D del gruppo O sono: O6 (macchia di Poisson), O8 (Airy), O12 (anelli di Newton),
O27 (Imbert–Fedorov), O29 (birifrangenza generale), O32 (fibra elicoidale), O33 (Beth), O34 (OAM);
più in parte O14 (fase di Gouy = π in 3-D, π/2 in 2-D: si testano entrambe). Quasi tutto
l'entanglement è 2-D: la "terza dimensione" di Bell è la polarizzazione, non lo spazio.

### 6.2 Per elemento ottico (quanti esperimenti lo richiedono)

| Elemento | N | Esperimenti |
|---|---|---|
| `SCH` schermo/rivelatore di posizione | 58 | quasi tutto O, più M3, M8, M14, M15, M25, M27, M29, L1–L4, L17, G1–G6, G9, G10 |
| `BS` beam-splitter non polarizzatore | 32 | O16, O20, O32, O34, O36–O39, M1, M4–M8, M14, M16–M21, M24, M26, M28, L6, L7, L15, L16, G3, G4, G9, G10 |
| `COINC` coincidenze | 32 | tutto M, L6, L15, L16 |
| `SRC2` stato a più fotoni dato | 30 | tutto M |
| `MAT` materia | 19 | tutto L tranne O31 (tag L ma nella sezione O), G5, G6 |
| `PBS` beam-splitter polarizzatore | 18 | O18, O36, O38, O39, M10–M13, M16–M23, M26, M30 |
| `WP` lamina d'onda / rotatore | 18 | O19, O27, O30, O32, O33, O38, O39, M6, M10, M15, M16, M18, M20–M23, M26, M30 |
| `S` specchio | 17 | O4, O5, O13, O16, O20, O32, O34, O37, M7, M14, L13, L14, G1, G3, G4, G9, G10 |
| `RIT` ritardo / sfasatore | 14 | O16, O36–O38, M5–M8, M16, M24, M26, M28, L7, L16 |
| `DI` interfaccia dielettrica (n dato) | 13 | O4, O12, O21–O28, M24, L1, G3 |
| `MOV` moto / rotazione / gravità | 13 | O33, L11, L13, tutto G |
| `POL` polarizzatore assorbente | 10 | O17–O19, O30, O31, M6, M9, M15, M29, L14 |
| `F2` doppia fenditura | 9 | O1, O2, O17, O35, M14, M15, M25, M27, M29 |
| `L` lente | 6 | O8, O12, O14, O34, O35, M27 |
| `EOM` commutatore veloce | 6 | O36, O37, M11–M13, M17 |
| `MB` mezzo anisotropo / chirale | 5 | O29, O31, O35, M29, L14 |
| `CAV` cavità | 5 | O11, L5, L9, L11, G10 |
| `RET` reticolo | 3 | O9, O10, L17 |
| `OST` ostacolo opaco | 3 | O7, O20, M8 |
| `F1` fenditura singola | 2 | O3, O7 |
| `AC` apertura circolare | 2 | O6, O8 |

(I conteggi sono fatti sui tag; un esperimento conta una volta per elemento. Sono indicativi di
priorità, non una specifica: alcuni elementi si possono sostituire con altri, es. `PBS` + `POL`.)

---

## 7. Ordine consigliato di implementazione

Il criterio: prima ciò che **sblocca più esperimenti** con **meno fisica nuova**, e prima i test
**decisivi** (§8) che possono far cadere il modello.

1. **Gratis con gli elementi di oggi (2-D, parete con aperture, POL, PBS).**
   O3 (singola fenditura), O7 (Babinet), O9 (reticolo), O10 (Talbot), O17 terzo caso (cancellino con
   analizzatore a 45°), O19 (Zeno con polarizzatori in sequenza). Nessun elemento nuovo; sei
   esperimenti in più, due dei quali (Talbot, Zeno) mettono alla prova il campo vicino e la moneta
   locale ripetuta.
2. **Specchio `S` e beam-splitter `BS` (+ `RIT`).** È il salto più redditizio: sblocca O5 (Lloyd),
   O13 (Wiener), O16 (Mach–Zehnder), O20 (bomba di Elitzur–Vaidman), O37, e in combinazione col
   punto 4 quasi tutto il gruppo M. Richiede che la regola pilota gestisca **ventagli riflessi**
   (oggi il ventaglio copre solo il semipiano in avanti) e oggetti non perpendicolari a z.
   O20 è il test in cui la località delle palline-onda pesa davvero.
3. **Lamine d'onda `WP` e commutatore `EOM`.** Solo matrici di Jones (già presenti) e un oggetto che
   cambia nel tempo. Sblocca O30 (tomografia), O36 (Wheeler, decisivo), O37, O38–O39
   (contestualità di un fotone solo: decisivi e ancora a un fotone).
4. **Due fotoni con onda bosonica comune (`SRC2` generale + `COINC`).** Il punto più difficile
   concettualmente: oggi `entangle.rs` gestisce solo la polarizzazione di coppie e in modo non
   locale per costruzione; la regola pilota ha una pallina-fotone per onda. Per M5 (HOM) serve
   un'onda **a due particelle simmetrica** nello spazio delle configurazioni (o una regola
   equivalente): senza, HOM fallisce con 1/2 di coincidenze. Sblocca M4–M8, M14–M30, cioè ~25
   esperimenti, tutti in 2-D. Ordine interno: M4 → M5 → M6 → M7 → M14 → M8 → M25/M26 → M18–M21 →
   M28.
5. **Modello di rivelatore** (efficienza η, conteggi di buio, tempo morto, finestra di
   coincidenza). Necessario per i numeri storici di M1, M2, M9, M13 e per dire onestamente se
   una violazione di Bell è "senza scappatoie". Piccolo da scrivere, ma va dichiarato.
6. **Interfaccia dielettrica `DI` con n dato** (rifrazione, riflessione parziale, TIR, onde
   evanescenti). Sblocca O4, O12, O21–O28, M24, G3. Richiede di rivedere il taglio Θ = π/2 − 10⁻³
   (le onde evanescenti di O24–O26 vivono proprio lì) e la moneta del flusso alle superfici
   parzialmente riflettenti.
7. **Lenti `L` e cavità `CAV`.** O8, O11, O14, M27, G10. La cavità richiede molti passaggi della
   stessa onda (ventagli che si rigenerano).
8. **Passaggio al 3-D** (origini su una griglia 2-D dell'apertura, ondina e^{ikR}/R): O6, O8, O27,
   O34 e la fase di Gouy completa. Costo ∝ (origini)² (vedi §7.5 di `REGOLA_FOTONE.md`).
9. **Mezzi anisotropi `MB`** con tensore dato (O29, O35 completo, M29).
10. **Moto e gravità `MOV`** (gruppo G). G4 e G9 sono interferometri (vengono quasi gratis dopo il
    punto 2, se l'apparato può muoversi o ruotare); G6–G8 richiedono una metrica.
11. **Gruppo L**, che coincide con i livelli 1–4 di `VISIONE.md`: L17 e L8 per primi (atomo di
    idrogeno), poi L10–L11 (vuoto), L2–L3 (elettrone libero).

Con i punti 1–4 si coprono circa **45 dei 96** esperimenti (circa 19 del gruppo O, compresi i 5 già validati, e 27 del gruppo M: tutti tranne M24, M27, M29 che chiedono `DI`, `L`, `MB`), tutti in 2-D. Con i punti 5–7 si arriva a circa 65.

---

## 8. Esperimenti decisivi

Se il modello ne fallisce anche uno solo, **non sta simulando un fotone** (o almeno non uno
quantistico). In ordine di implementazione:

| # | Esperimento | Che cosa non può fallire | Stato |
|---|---|---|---|
| 1 | O1/O2 Young + Taylor | frange da un fotone alla volta, con clic singoli | VALIDATO |
| 2 | O17 Fresnel–Arago + cancellino | la polarizzazione ortogonale spegne le frange, il proiettore le riaccende | 2/3 validato |
| 3 | O18 Malus | cos²θ per un fotone **intero** (passa tutto o niente) | VALIDATO |
| 4 | O16 Mach–Zehnder a fotone singolo | P = cos²(φ/2) con un fotone solo | da fare |
| 5 | **M4 Grangier–Roger–Aspect** | un fotone a un BS fa **un solo** clic (α ≈ 0), *e* lo stesso fotone interferisce nel MZ | da fare |
| 6 | O20 Elitzur–Vaidman | un fotone sa che l'altro braccio è chiuso senza passarci | da fare |
| 7 | O36/O37 Wheeler a scelta ritardata, V² + D² ≤ 1 | il risultato dipende solo dall'apparato al momento della misura | da fare |
| 8 | O38 KCBS con un fotone | contestualità anche senza entanglement | da fare |
| 9 | **M5/M6 Hong–Ou–Mandel** | due fotoni identici escono insieme: coincidenze → 0, e la distinguibilità riduce la buca come cos²θ | da fare: **ostacolo principale** |
| 10 | M10–M12 CHSH | S = 2√2 e nessun segnale tra i lati | VALIDATO (ideale, non locale per costruzione) |
| 11 | M13 Bell senza scappatoie | violazione anche con efficienza realistica > 2/3 | da fare, [DET] |
| 12 | M14/M15 gomma quantistica a scelta ritardata | il D₀ totale non ha mai frange; i sottoinsiemi sì, con segni opposti | da fare |
| 13 | M8 Zou–Wang–Mandel | V del segnale = \|T\| dell'idler mai rivelato | da fare |
| 14 | M20 GHZ | contraddizione "a colpo singolo" | da fare |

Test **specifici per la regola pilota** (non decisivi per "è un fotone", ma decisivi per "è questa
regola"): **O35 (Kocsis)** — le traiettorie della pallina-fotone devono coincidere con quelle medie
misurate e non attraversare l'asse; **M29 (Mahler)** — la regola estesa a coppie deve dare traiettorie
che dipendono dalla misura lontana; **O10 (Talbot) e O25 (TIR frustrata)** — campo vicino ed
evanescente.

---

## 9. Onestà: cosa dipende da cose che non abbiamo

* **Rivelatori.** Oggi il clic è ideale. I numeri storici di M1, M2, M3, M4, M9, M12, M13, L6, O13
  dipendono da efficienza, rumore e risoluzione temporale. Per questi il PASS va dato sulla
  **previsione ideale**; per confrontarsi con il numero storico serve il modello di rivelatore (§7,
  punto 5). M13 *senza* modello di rivelatore non ha senso: la scappatoia del rivelatore è il suo
  contenuto.
* **Sorgenti.** Conversione parametrica, cascate atomiche, punti quantici, centri NV non si
  simulano: sono **stati dati** [SRC]. Tutto il gruppo M è quindi un test della propagazione e della
  misura, non della generazione.
* **Materiali.** Tutto ciò che contiene `DI`, `MB` usa un n o un tensore **dato**; ricavarli dagli
  atomi è il lavoro dei livelli 1–3. Il criterio PASS di quelle schede vale per n dato.
* **Entanglement e località.** `entangle.rs` è non locale per costruzione (necessario per il teorema
  di Bell, vedi `VISIONE.md`). Passare M10–M13 con esso dimostra che la statistica è giusta, non
  che una regola locale la produce (non può).
* **Numeri "indicativi".** Sono marcati nelle schede: O16 (visibilità 98%), M3, M5 (visibilità
  storica), M14, M18, M19, M23, M29, O26, O27, O32–O34, O36, O39, L6, L9, L12, L13, G5, G7 (VLBI),
  M1 (diametro di Sirio). In tutti questi casi il PASS è la previsione teorica esatta.
* **Numero di fotoni variabile.** La regola attuale ha un fotone per evento (o una coppia). Luce
  compressa (L7), statistica laser/termica (M2), HBT (M1) richiedono stati con numero di fotoni
  non definito: vanno rappresentati come miscele/sovrapposizioni di eventi a N fotoni, cosa che oggi
  non esiste.
