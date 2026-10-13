# Protocollo sperimentale per l'orologio di massa nel tunnel

Studio `tunproto_` (file: `lab/tunproto_*.py`, `lab/out/tunproto_*`). Segue `docs/TEMPO_TUNNEL.md` (giro 1).
Commit: ricerca `f2932a2`, pre-registrazione `d3936f6`, corse e giudizio `ae1e80a`, documento e figura: vedi `git log`.
Figura: `lab/out/tunproto_figura.png` (copia in `/mnt/project-files/problemi-aperti/tunproto_figura.png`).

## 0. In breve (per chi ha fretta)

* **La misura esiste ed è semplice**: un qubit "ancella" che cambia la **massa** di ±δ/2 solo dentro la barriera,
  letto sui soli trasmessi, dà esattamente l'orologio di massa del giro 1: inclinazione s_z = −tanh(δτ_mz), angolo
  nel piano = −δτ_my. Dimostrato tick per tick sulla regola (uguale alla previsione a **3·10⁻⁶**), sul Dirac
  continuo (split-step, **8·10⁻⁴**), e nel limite debole/stazionario (A1–A5 tutti PASS). La parte inclinata non
  ha nemmeno bisogno dell'ancella: è −½ ∂ln T/∂m_B.
* **Piattaforma 1, cammino quantistico in fibra = la regola stessa** (θ = 0,5, 40 passi, moneta programmabile
  per posto e passo): previsione nativa della regola, esatta e nuova di questo studio: al centro del gap
  l'orologio di massa cresce di **1/cos θ = 1/v_max passi per sito** (il tunnel dura quanto la luce *del reticolo*);
  il rapporto nel piano è **tan(εa)/tan θ**; tutto verificato a 10⁻⁶. Con un pacchetto di 4 siti la pendenza
  "misurata" è 1,153 (orologio di potenziale: 0,142). **Ripetizioni per 5σ**: ~1 200 trasmessi per punto per
  orologio (pendenza massa ≠ pendenza potenziale), ~200 con un solo spessore; il segno negativo a V0 = 1,2 E
  vuole ~30 000–40 000 trasmessi (il Monte Carlo a 29 000 dà 4,3σ: **FAIL** del mio conto analitico).
* **Il vero nemico è il rumore, non la statistica** (scoperta post-hoc, dichiarata): un rumore di fase che cambia
  da passo a passo non conserva l'energia e fa passare un po' di pacchetto **sopra** la barriera; con T ~ 10⁻³ basta
  pochissimo. Per stare entro il 10 % servono **≤ 0,003 rad rms** di fase e **≤ 0,005 rad** di moneta per sito e
  per passo. Ai livelli che avevo pre-registrato (0,05 e 0,01 rad) l'esperimento **fallisce** (pendenza perfino
  negativa). La perdita dipendente dalla polarizzazione (3 %/passo) invece è innocua.
* **Piattaforma 2, guide d'onda binarie** (polarizzazione = ancella, birifrangenza alternata Δn ≈ 1,5·10⁻⁶ dentro
  la barriera, 12 cm): l'orologio di massa segue L/c con pendenza 1,11; disordine di fabbricazione (0,4 %) e
  accoppiamento dipendente dalla polarizzazione (−7 %, sottraibile a 2·10⁻⁴) sono gestibili. Ma il reticolo delle
  guide **non è Dirac continuo** a L = 4 Compton (scarto fino al 13 %: FAIL della mia soglia del 10 %): la
  previsione per l'esperimento va fatta col modello delle guide, non con Dirac.
* **Onestà sulla novità**: l'idea "orologio attaccato alla massa" per il tunnel esiste già, non relativistica
  (Schach e Giese, *Sci. Adv.* 2024). Nostro: il caso di Dirac, dove l'orologio di massa e quello di potenziale
  danno numeri diversi (L/c contro ~0 al centro del gap), e la versione esatta sul reticolo. Un simulatore di Dirac
  ben fatto **non può smentire Dirac**: l'esperimento controllerebbe la lettura (rotazioni = tempo proprio) e la
  qualità del simulatore, non una legge nuova.
* **Esami: 21 pre-registrati, 14 PASS, 7 FAIL** (B1, C1, C4, tre C5, D1); i FAIL sono spiegati in §4. Parametri
  liberi: **0**.


## 1. Ricerca (13/10/2026)

### 1.1 Da dove partiamo
Nel giro 1 (`TEMPO_TUNNEL.md`, §4.4) abbiamo visto che contare le rotazioni della zig-zag dentro la barriera,
pesate con le ampiezze, è esattamente un **orologio di Larmor attaccato alla massa** (T_m = i ∂ln t/∂m_B, con m_B
la massa solo dentro la barriera), e che nel continuo di Dirac questo orologio legge ψ̄ψ, cioè il ritmo del
**tempo proprio**. Previsioni (barriera spessa, |ε| < m, ε = E − V0, κ = √(m² − ε²), ħ = c = 1):
* parte "nel piano": τ_my = (ε/m) τ_y (va all'indietro per V0 > E);
* parte "inclinata": dτ_mz/dL = m/κ ≥ 1/c, quindi τ_mz → L/c esatto a V0 = E, dove Büttiker–Landauer dà 0;
* barriere sottili: τ_mz < L/c (0,851 L a L = 1).

Qui la domanda è: **si può misurare davvero?** Con quale apparato, quali numeri, quante ripetizioni.

### 1.2 Simulatori di Dirac 1+1 D (cosa esiste davvero)

| piattaforma | Dirac come | m ↔ | c ↔ | potenziale V(x) | stato | fonte |
|---|---|---|---|---|---|---|
| **ione intrappolato** ⁴⁰Ca⁺ | spinore = 2 livelli interni, x = moto (modo assiale) | ħΩ (portante), Ω ≤ 2π·13 kHz | c = 2ηΩ̃Δ = 0,052 Δ/µs (η = 0,06, Ω̃ = 2π·68 kHz, ω_ax = 2π·1,36 MHz) | **solo lineare** (e quadratico) in x, con un secondo ione ancella | Zitterbewegung misurata fino a 150 µs; 10⁴–3·10⁴ misure per punto | Gerritsma et al., *Nature* 463, 68 (2010) (arXiv 0909.0674, letto); Casanova et al., arXiv 1004.5400 (letto); Gerritsma et al., *PRL* 106, 060503 (2011) (abstract letto) |
| **guide d'onda binarie** (silice, scritte col laser fs) | due sottoreticoli = due componenti, z (lunghezza di propagazione) = tempo | σ = metà del salto alternato delle costanti di propagazione | κ (accoppiamento tra guide): c = 2κ guide per unità di z | indice di rifrazione cambiato su un gruppo di guide (gradino/barriera) | Zitterbewegung (Dreisow et al., *PRL* 105, 143902 (2010)); **tunnel di Klein osservato** con un gradino, "diverse masse" (Dreisow, Keil, Tünnermann, Nolte, Longhi, Szameit, *EPL* 97, 10008 (2012), abstract letto) | Longhi, *Opt. Lett.* 35, 235 (2010) (arXiv 0912.5071, letto: passo 10 µm, λ = 633 nm, lunghezza 10 cm); Longhi, *PRB* 81, 075102 (2010); rassegna arXiv 1111.3461 (letta) |
| **cammino quantistico fotonico in fibra** (multiplazione nel tempo) | **è la nostra regola**: posto = intervallo di tempo, R/L = polarizzazione H/V (o anello corto/lungo), moneta = rotazione | angolo della moneta θ (lamina o EOM) | 1 posto per passo (fisso) | fase di un modulatore elettro-ottico (EOM) per posto e passo | Schreiber et al., *PRL* 104, 050502 (2010): 5 passi, moneta da lamina; Schreiber et al., *PRL* 106, 180403 (2011): 28 passi, fasi per posto/passo con EOM (rapporto φ_V/φ_H ≈ 3,5 fisso), > 10⁴ eventi per passo in ≤ 1 h; Regensburger et al., *Nature* 488, 167 (2012): due anelli, accoppiatore fisso 50:50 (θ = π/4), fasi con modulatore al niobato di litio, 70–80 passi con amplificatori; Boutari et al. (arXiv 1607.00891): 62 passi; Zhang et al. (Nanjing, arXiv 2202.09630): moneta programmabile per intervallo di tempo (Sagnac + EOM), 11 passi | letti online (abstract/testo) |
| **atomi freddi con accoppiamento spin-orbita** (⁸⁷Rb, Raman) | spinore = 2 stati iperfini | Ω_Raman/2 | ħk_R/M (~ mm/s) | barriere ottiche (le stesse di Ramos et al. 2020) | Zitterbewegung: LeBlanc et al., *NJP* 15, 073011 (2013); Qu et al., *PRA* 88, 021604 (2013) (**da memoria, non verificati online**) | — |
| **circuiti superconduttori** | catena di qubit nel sottospazio a un'eccitazione (dimerizzata: SSH/Dirac) | frequenze alternate | accoppiamento J | frequenze dei qubit (flusso) | cammini quantistici su 12–62 qubit (Yan et al., *Science* 364, 753 (2019); Gong et al., *Science* 372, 948 (2021); **da memoria**) | — |

**Lettura onesta della tabella.**
* Gli **ioni** sono il simulatore di Dirac più pulito, ma il potenziale è una funzione *lineare* di x (un campo
  elettrico): **una barriera rettangolare non è stata proposta né fatta** (Casanova et al. trattano lineare e
  quadratico). La "barriera" degli ioni sarebbe la zona di gap di un potenziale lineare (regime di Sauter): bordi
  lisci per forza, larghezza 2m/g. Si può fare, ma l'orologio va acceso solo nella zona di gap e oggi non c'è un
  modo dimostrato. Li tengo come seconda scelta.
* Le **guide d'onda binarie** hanno già fatto proprio il tunnel di Klein contro un gradino, con masse diverse: la
  barriera c'è. Luce classica (niente rumore di conteggio), ma ogni campione è un pezzo di vetro (difetti di
  fabbricazione) e il "tempo" totale è la lunghezza del campione (10–15 cm).
* Il **cammino in fibra** è letteralmente la regola zig-zag: corri (ritardo H/V di un posto) e ruota (moneta).
  Qui la previsione si fa sulla regola **a θ finito**, senza bisogno del limite del continuo. Il limite è il
  numero di passi: 28 con moneta+EOM programmabili (Schreiber 2011), 62–80 senza moneta variabile.
* Atomi con spin-orbita: il termine p²/2M rovina il Dirac puro; circuiti: plausibili ma niente di dimostrato
  su questo. Non li scelgo.

### 1.3 Come si legge un orologio (operativamente)
**Steinberg / Ramos.** Ramos, Spierings, Racicot, Steinberg, *Nature* 583, 529 (2020): l'orologio è un grado di
libertà interno (due stati iperfini del ⁸⁷Rb) che ruota solo dentro la barriera (Raman accoppiata dalla luce
della barriera stessa); si misura lo stato interno dei soli atomi trasmessi: angolo nel piano (τ_y) e
inclinazione (τ_z). Valore debole complesso: Steinberg, *PRL* 74, 2405 (1995).

**La nostra versione.** Un'**ancella** (qubit) con stati |0⟩, |1⟩ cambia la **massa** dentro la barriera di ± δ/2:
H = σ_z p + [m + (δ/2) s_z χ_B(x)] σ_x + V(x) (s_z = Pauli dell'ancella, χ_B = 1 dentro la barriera). Si prepara
l'ancella in (|0⟩ + |1⟩)/√2 e si misura sui soli trasmessi. Nel regime stazionario (dimostrazione in §2):
  s_z = −tanh(δ τ_mz),   s_x + i s_y = e^{−iδ τ_my} / cosh(δ τ_mz),
cioè l'angolo nel piano dà τ_my e l'inclinazione dà τ_mz, **esattamente** come l'orologio di Büttiker ma con la
massa al posto del potenziale. La parte inclinata è anche leggibile **senza ancella**: τ_mz = −½ ∂ln T/∂m_B
(quanto la probabilità trasmessa cambia se si cambia la massa dentro la barriera).

### 1.4 È davvero nuovo? (cosa ho trovato)
* **Schach e Giese, *Sci. Adv.* 10, eadl6078 (2024)** (arXiv 2404.14382, letto su ar5iv): orologio di **Ramsey**
  per il tunnel con **difetto di massa** (stati interni di massa diversa, m_{e/g} = m̄(1 ± Δω/2ω̄)): è l'idea
  "orologio attaccato alla massa" più vicina alla nostra. Differenze: (i) loro sono non relativistici (Schrödinger
  con m c² + p²/2m), noi Dirac pieno; (ii) il loro difetto di massa è ovunque (si sottrae la dilatazione del
  tempo con un riferimento), il nostro è solo dentro la barriera; (iii) loro notano che nel limite non
  relativistico il difetto di massa equivale a una barriera dipendente dallo stato (entra il prodotto m̄V0): è
  proprio il nostro limite NR, dove l'orologio di massa torna Larmor. Concludono "niente tunnel istantaneo né
  superluminale", nello spirito del nostro τ_mz = L/c al centro del gap. **Quindi l'idea di un orologio di massa
  per il tunnel non è nostra**; nostro è il caso di Dirac (dove l'orologio di massa è davvero diverso da quello di
  potenziale) con le tre proprietà (ε/m, m/κ, L/c).
* Larmor relativistico: Li, Liang, Kobe (quant-ph/0111117) e Li, Liang (quant-ph/0109001): particella neutra di
  spin ½ con momento magnetico: orologio sul **potenziale** (dà il dwell), non sulla massa.
* Tempi di transito per Dirac con pacchetti: Krekora, Su, Grobe, *PRA* 63, 032107 (2001); De Leo (arXiv
  1304.1561, solo abstract: velocità di transito ∝ spessore nel limite opaco); De Leo e Rotelli; nessuno usa la
  massa come orologio.
* Sokolovski e Akhmatskaya, *Sci. Rep.* 11, 10040 (2021) (letto): ogni tempo di Larmor è la risposta
  dell'ampiezza a una piccola perturbazione (un valore debole), **non un intervallo di tempo fisico**. Vale anche
  per il nostro: τ_m è la risposta alla massa. Lo dico chiaro nel documento.
* Feynman a scacchiera: Jacobson e Schulman, *J. Phys. A* 17, 375 (1984): il numero di svolte è legato al tempo
  proprio/limite NR (da `TEMPO_TUNNEL.md`). Collegamento noto in generale, non applicato al tunnel (per quanto
  trovato).
* Orologi quantistici e tempo proprio: Zych, Costa, Pikovski, Brukner, *Nat. Commun.* 2, 505 (2011); Smith e
  Ahmadi, *Nat. Commun.* 11, 5360 (2020) (**da memoria**): un orologio interno accoppiato via E = mc² legge il
  tempo proprio. La nostra lettura "orologio di massa = tempo proprio" è la stessa idea, nel tunnel.

**Conclusione della ricerca.** La misura proposta è un orologio di Larmor/Ramsey con l'accoppiamento spostato dal
potenziale alla massa, solo dentro la barriera, in un simulatore di Dirac. L'idea generale (orologio di massa)
esiste già (Schach–Giese 2024, NR); il caso di Dirac con le sue tre previsioni non l'ho trovato (8 ricerche web,
non esaustive). Attenzione: una piattaforma che simula bene Dirac **non può** smentire Dirac; l'esperimento
controlla (a) che il simulatore è Dirac anche nell'evanescente, (b) la lettura "rotazioni = tempo proprio", e
(c) che l'orologio di massa e quello di potenziale danno numeri diversi (L/c contro 0).

---

## 2. L'osservabile operativo (e perché ridà τ_m)

**L'apparato, in parole semplici.** Prendi la particella di Dirac e un qubit in più (l'ancella). Dentro la
barriera, e solo lì, la massa vale m + δ/2 se l'ancella è in |0⟩ e m − δ/2 se è in |1⟩. Prepari l'ancella "a metà"
((|0⟩+|1⟩)/√2), mandi il pacchetto, e alla fine guardi l'ancella **solo quando la particella è passata**.

**Il conto (due righe).** Ogni ramo dell'ancella vede una barriera leggermente diversa, quindi passa con ampiezza
t(m ± δ/2) = t·e^{∓iδT_m/2}, con T_m = τ_my − iτ_mz = i ∂ln t/∂m_B (la definizione del giro 1). La coerenza tra i
rami e le popolazioni danno
  s_x + i s_y = e^{−iδτ_my}/cosh(δτ_mz),   s_z = −tanh(δτ_mz),
quindi τ_my = −(angolo nel piano)/δ e τ_mz = −atanh(s_z)/δ. Nel piano l'ancella "ticchetta" alla frequenza δ per un
tempo τ_my (orologio che gira); in z si inclina perché il ramo più pesante passa meno (filtro). È lo schema di
Büttiker/Steinberg/Ramos con la massa al posto del potenziale. Nel linguaggio di Steinberg, T_m è il **valore
debole** di ∫β χ_B dt sui trasmessi (β = la matrice della massa, σ_x qui): "tempo proprio passato nella barriera".

**Le prove numeriche (pre-registrate, sezione A, E0 = 2 mc², L = 4 ħ/mc, σ = 10, δ = 0,05):**

| esame | cosa | risultato | soglia | esito |
|---|---|---|---|---|
| A1 | regola tick per tick con l'ancella (2 rami = 2 copie) = previsione dallo stazionario della regola | 3,4·10⁻⁶ | 10⁻⁵ | **PASS** |
| A2 | regola a a → 0 (0,05 / 0,025 / 0,0125) = Dirac continuo | 6,5·10⁻⁴ (estrapolazione lineare, come registrato); 3,4·10⁻⁶ con quella quadratica (post-hoc: l'errore della regola qui è O(a²)) | 10⁻³ | **PASS** |
| A3 | Dirac continuo indipendente (split-step di Fourier, dx 0,05, dt 0,0125) = Dirac spettrale | τ: 8,4·10⁻⁴; P_T: 8,6·10⁻⁴ rel. | 5·10⁻³; 10⁻³ | **PASS** |
| A4 | δ → 0 e pacchetto largo (σ = 80) → τ stazionario a E0 | scarto ≤ 10⁻² e decrescente con σ | — | **PASS** |
| A5 | la lettura conferma le previsioni del giro 1 | τ_mz/L = **1,0005** al centro (potenziale: τ_z/L = 0,127); a V0 = 2,4 τ_my < 0 e τ_y > 0; τ_my/τ_y = 0,411 e −0,396 (ε/m = ±0,4) | vedi prereg | **PASS** |

Un dettaglio utile: τ_z (Büttiker–Landauer) al centro del gap **non è zero** a L finito: è ≈ 0,5 ħ/mc² e non cresce
con L. Quindi la discriminazione pulita è sulla **pendenza**: dτ_mz/dL = 1/c contro dτ_z/dL = 0.

## 3. Il protocollo pre-registrato (`lab/out/tunproto_prereg.json`)

### 3.1 Piattaforma 1: cammino quantistico fotonico in fibra (= la regola)

| cosa | nella regola | nell'apparato |
|---|---|---|
| posto x | una cella | un intervallo di tempo (time bin) |
| R / L | due palline | polarizzazione H / V (o anello corto/lungo) |
| corri a c | passo di una cella | ritardo differenziale di un bin per giro |
| ruota (massa) | θ = 0,5 rad | moneta: lamine + EOM (tipo Sagnac di Zhang et al. 2022), programmabile per bin |
| potenziale V0 | fase e^{−iV0 a} per passo | fase dell'EOM nei bin della barriera: V0 a = 1,0 rad (al centro del gap) |
| barriera | N celle | N = 2, 4, 6, 8 bin |
| orologio di massa | θ ± δa/2 nella barriera | moneta ± **0,0375 rad** solo nei bin della barriera |
| orologio di potenziale | V ± δ/2 | fase EOM ± 0,0375 rad |
| ancella | 2 copie | un qubit a "macro-bin" (copia anticipata/ritardata del treno) con EOM diverso nei due rami; ricombinazione interferometrica alla fine |
| pacchetto | gaussiano, σ = 4 siti, solo energia positiva | treno di impulsi con ampiezza e polarizzazione programmate, tagliato a ±4σ (16 bin) |
| durata | 40 passi | 40 giri (oggi: 28 con moneta programmabile, 62–80 a moneta fissa) |

**Previsioni native della regola (barriera spessa, esatte, verificate a 10⁻⁶ su θ = 0,3 / 0,5 / 0,7):**
con cosh(κa) = cos(εa)/cos θ,
* pendenza orologio di massa dτ_mz/dL = cos(εa) sin θ / (cos²θ sinh κa) → al centro del gap **1/cos θ**, cioè
  **1/v_max** (v_max = c cos θ è la velocità massima dei pacchetti della regola): *il tunnel, contato in rotazioni,
  dura quanto la luce del reticolo*. Nel continuo (θ → 0) torna L/c;
* pendenza orologio di potenziale sin(εa)/(cos θ sinh κa) → 0 al centro;
* rapporto nel piano τ_my/τ_y = **tan(εa)/tan θ** → ε/m nel continuo.

**Previsioni per l'esperimento (pacchetto vero, 40 passi; tempi in passi, barre = solo statistica a N₅σ):**

| V0/E | N | P_T | τ_mz (massa) | τ_z (potenziale) | τ_my (massa) | τ_y (potenziale) |
|---|---|---|---|---|---|---|
| 1,0 | 2 | 0,31 | 2,01 | 0,82 | 0,53 | 1,57 |
| 1,0 | 4 | 0,047 | 4,56 | 1,29 | 0,30 | 1,83 |
| 1,0 | 6 | 6,8·10⁻³ | 6,96 | 1,57 | 0,26 | 1,91 |
| 1,0 | 8 | 1,0·10⁻³ | 9,17 | 1,86 | 0,38 | 2,12 |
| 1,2 | 8 | 1,1·10⁻³ | 9,29 | −2,61 | **−0,60 ± 0,14** | 1,75 |
| 0,8 | 8 | 5,9·10⁻³ | 10,85 | 6,76 | 3,80 | 5,21 |

(Tutti i valori in `lab/out/tunproto_rC_pulito.json`; τ in passi = τ[ħ/mc²]/a.) Pendenza al centro (N = 4..8):
**1,153** contro 1/cos θ = 1,139 (+1,2 % per il pacchetto finito); potenziale **0,142**.

**Quante ripetizioni (contando solo il rumore di proiezione dell'ancella):**

| discriminazione | N₅σ (trasmessi rivelati) | Monte Carlo 2 000 esperimenti | esito |
|---|---|---|---|
| D1: pendenza massa ≠ pendenza potenziale (V0 = E, N = 4, 6, 8), δ = 0,15 | 1 174 per punto per orologio | z mediano 4,85 (16° perc. 3,85) | C3 **PASS** |
| D1, δ = 0,6 | 768 | 4,81 | C3b **PASS** |
| D1b: un solo spessore N = 8, τ_mz contro τ_z | 206 per orologio | (analitico) | — |
| D2: τ_my < 0 a V0 = 1,2 E, N = 8, δ = 0,15 | 28 959 | **4,30** | C4 **FAIL** |
| D2, δ = 0,6 | 76 014 | 4,99 | C4b **PASS** |

Perché C4 fallisce: il mio N₅σ usava la previsione "spettrale" (pacchetto gaussiano intero, τ_my = −0,353 in
unità ħ/mc²), ma il pacchetto tagliato a 4σ dà −0,301: segnale più piccolo del 15 %, servono ~ 39 000 trasmessi.
Con P_T ≈ 10⁻³ vuol dire ~4·10⁷ particelle mandate per punto. Con luce laser attenuata (equivalente al singolo fotone
in un apparato lineare) e ~10⁴ eventi/ora come Schreiber 2011, D1 richiede minuti, D2 qualche ora per punto; con
luce classica e fotodiodi (Regensburger) la statistica non è un problema.
Notare: δ grande non aiuta nel piano (a δ = 0,6 l'inclinazione δτ_mz ≈ 2,8 rad spegne il contrasto, 1/cosh).

**Rumore (pre-registrato C5, poi scansione post-hoc):**

| sistematico | livello registrato | effetto sulla pendenza di τ_mz | esito |
|---|---|---|---|
| fase bianca per bin e passo, uguale nei due rami | 0,05 rad | −148 % (pendenza negativa!) | **FAIL** |
| fase bianca, diversa nei due rami | 0,05 rad | −157 %, contrasto 0,05 | **FAIL** |
| angolo della moneta bianco | 0,01 rad | −31 % | **FAIL** |
| perdita dipendente dalla polarizzazione | 3 % per passo | −0,6 %; segno e contrasto ok | **PASS** |

Scansione post-hoc (`tunproto_rC_scan.json`, pannelli e–f): pendenza entro il 10 % per fase ≤ **0,002–0,003 rad**
e moneta ≤ **0,005 rad** rms per sito e per passo; l'effetto va come (rumore)². Meccanismo: il rumore che cambia
nel tempo dà energia al pacchetto; la parte "riscaldata" passa sopra la barriera (P_T a N = 8 cresce da 1,0·10⁻³ a
1,5·10⁻² con 0,05 rad) e non porta il segnale dell'orologio. Contromisure: misurare P_T (deve tornare la previsione
del tunnel), filtrare in tempo di volo (la parte riscaldata arriva prima), barriere più sottili (T più grande).
**Questo è il requisito vero del protocollo**, più della statistica.

### 3.2 Piattaforma 2: guide d'onda binarie (silice scritta al laser)

| cosa | Dirac | guide |
|---|---|---|
| tempo | t | lunghezza di propagazione z (12 cm) |
| c | 1 | 2κ = 9 guide/cm (κ = 4,5 cm⁻¹, ipotesi) |
| massa | m | σ = 1,5 cm⁻¹ (metà del salto alternato della costante di propagazione); Compton = 6 guide |
| E0 | 2 mc² | fascio a β = 3 cm⁻¹ sopra il centro di banda, inclinato vicino all'angolo di Bragg, solo banda alta (fascio sagomato) |
| barriera | V0, L | indice aumentato su 12/18/24 guide, Φ = 2,4 / 3,0 / 3,6 cm⁻¹ (Δn ≈ 2·10⁻⁵) |
| ancella | qubit | **polarizzazione** H/V del fascio |
| orologio di massa | m ± δ/2 | birifrangenza **alternata** solo nelle guide della barriera: δ = 0,3 cm⁻¹ (Δn ≈ 1,5·10⁻⁶) |
| lettura | s_x, s_y, s_z | parametri di Stokes della luce uscita dalle guide oltre la barriera |

Risultati (modi accoppiati = modello dell'apparato, contro Dirac continuo con le stesse unità):

| V0/E | L (Compton) | τ_mz guide | τ_mz Dirac | τ_z guide | τ_z Dirac |
|---|---|---|---|---|---|
| 1 | 2 | 1,93 | 1,99 | 0,54 | 0,59 |
| 1 | 3 | 3,04 | 3,14 | 0,77 | 0,90 |
| 1 | 4 | 4,16 | 4,36 | 1,19 | 1,62 |
| 1,2 | 4 | 4,38 | 5,01 | −1,88 | −2,98 |

(tempi in 1/σ = 0,67 cm di propagazione). Pendenza al centro **1,11** (soglia 1 ± 15 %: ok) ma scarto massimo
**12,7 %** (D1 **FAIL**, soglia 10 %): per barriere spesse il reticolo delle guide filtra lo spettro in modo diverso
dal continuo. Disordine di fabbricazione 0,05 cm⁻¹: media −0,26 %, dispersione tra campioni 0,4 % (D3 **PASS**).
Accoppiamento dipendente dalla polarizzazione dell'1 %: sposta τ_mz del −7 %, e la corsa a orologio spento lo
toglie fino a 1,6·10⁻⁴ (D2 **PASS**). Statistica: luce classica, il limite è la polarimetria (~10⁻³ sui Stokes →
~0,01 su τ), non il conteggio. Rumore temporale: il vetro è fermo, quindi niente "riscaldamento": vantaggio vero.

### 3.3 Le altre piattaforme (perché no, per ora)
* **Ioni** (Gerritsma 2010: c = 0,052 Δ/µs, mc²/h ≤ 13 kHz → Compton ≈ 0,64 Δ, ~12 tempi di massa in 150 µs):
  Dirac ottimo ma potenziali solo lineari/quadratici; il tunnel sarebbe quello di Sauter (bordi lisci, larghezza
  2m/g) e l'orologio andrebbe acceso solo nella zona di gap: non c'è un modo dimostrato. I bordi lisci di per sé
  non sono un problema (B1r: pendenza entro 0,1 % per w ≤ 0,5 ħ/mc, 1,5 % per w = 1, se l'orologio copre il gap).
* **Atomi con spin-orbita**: il termine p²/2M rovina Dirac; **circuiti superconduttori**: plausibili (frequenze
  alternate, ancella con accoppiamento ZZ di segno alternato) ma nessun dato su cui appoggiarsi.

## 4. Tutti gli esami contro la pre-registrazione

| esame | esito | nota |
|---|---|---|
| A1–A5 | 5 **PASS** | la misura con l'ancella = τ_m (regola e Dirac) |
| B1 | **FAIL** (2,4·10⁻⁴ > 10⁻⁴, bordi lisci w = 1 col profilo) | colpa della *previsione*: il calcolo stazionario tagliava la coda della tanh a ±4w; allargando a ±8w (post-hoc) l'accordo è 2·10⁻⁷ |
| B1r | **PASS** | bordi lisci innocui per w ≤ 0,5 |
| B2 | **PASS** | δ finito: scarto 3·10⁻⁵ a δτ = 0,8; 0,1 % anche a δτ = 3,2 (il conto a differenza centrale è molto robusto) |
| C1 | **FAIL** (7 letture su 48, tutte a N = 8) | il pacchetto tagliato a ±4σ ha code spettrali che passano sopra la barriera quando T ~ 10⁻³; l'inclinazione di massa (l'osservabile principale) resta sempre entro tolleranza. Lezione (vista già in fase di progetto: con taglio a 3σ la lettura a N = 8 era sbagliata del 45 %): **lo stato preparato deve avere potenza sopra la barriera ≪ T** |
| C2, C3, C3b, C4b | **PASS** | pendenza 1,153; 5σ con ~1 200 e ~800 trasmessi |
| C4 | **FAIL** | 4,3σ invece di 5 (N₅σ sottostimato del ~35 %) |
| C5 | 1 **PASS** (pdl), 3 **FAIL** | il rumore temporale è il requisito critico (§3.1) |
| D1 | **FAIL** | guide ≠ Dirac continuo fino al 13 % |
| D2, D3 | **PASS** | contaminazione e disordine sotto controllo |

Totale **14 PASS su 21**. Nessuna soglia cambiata dopo le corse.

## 5. È davvero diverso da cose note? (onestà)

* **Dal Larmor di Büttiker (τ_y, τ_z).** Matematicamente è lo stesso schema (derivata di ln t rispetto a un
  parametro dentro la barriera), con un parametro diverso: la massa invece del potenziale. In Schrödinger le due
  cose coincidono (nel limite NR il difetto di massa equivale a una barriera, lo notano Schach e Giese); in Dirac
  no: pendenze m/κ contro ε/κ, rapporto nel piano ε/m. La parte inclinata è "solo" −½ ∂ln T/∂m_B: chi ha misurato
  il tunnel di Klein con masse diverse (Dreisow et al. 2012) ha in mano qualcosa di vicino, ma con la massa
  cambiata ovunque, non solo nella barriera.
* **Dai valori deboli di Steinberg.** È un valore debole (di ∫βχ_B dt invece che di ∫χ_B dt). Vale la critica di
  Sokolovski–Akhmatskaya (2021): è la risposta dell'ampiezza a una perturbazione, **non un intervallo di tempo in
  senso classico**. Chiamarlo "tempo proprio" è una lettura, giustificata dal caso libero (dà esattamente
  la dilatazione del tempo) e dal conteggio delle rotazioni della scacchiera.
* **Dal tempo proprio dei cammini di Feynman a scacchiera.** Che il numero di svolte sia legato al tempo proprio è
  noto (Jacobson–Schulman 1984); l'applicazione al tunnel e il valore 1/cos θ = 1/v_max sul reticolo non li ho
  trovati.
* **Da Krekora–Su–Grobe.** Loro fanno pacchetti di Dirac e tempi di picco/fronte; nessun orologio di massa.
* **Qualcuno l'ha già proposto?** L'orologio di massa per il tunnel sì, in versione non relativistica e con il
  difetto di massa ovunque (Schach e Giese, *Sci. Adv.* 10, eadl6078 (2024)), con la stessa conclusione "niente
  tunnel istantaneo". La versione di Dirac con l'ancella solo dentro la barriera, le tre previsioni e il protocollo
  in fibra/guide non li ho trovati (ricerca web limitata, ~15 interrogazioni). Prima di dirlo in pubblico va
  cercato meglio (Ord, Kauffman–Noyes sulla scacchiera; letteratura sul "proper time" in Dirac).
* **Validazione, scoperta, limite.** VALIDAZIONE: la misura ridà τ_m (A1–A5), Büttiker nel limite noto.
  SCOPERTA (piccola, controllata a 10⁻⁶): le formule esatte della regola 1/cos θ e tan(εa)/tan θ, cioè la previsione
  per un cammino quantistico reale senza passare dal continuo; il requisito quantitativo sul rumore (≤ 0,003 rad).
  LIMITI: un simulatore di Dirac non può smentire Dirac; le guide non sono Dirac a L = 4 Compton; il protocollo in
  fibra chiede 40 passi con moneta programmabile (dimostrati 28) e rumore 10–20 volte sotto quello che avevo
  ipotizzato.

## 6. Debiti (cose messe a mano)

| debito | perché |
|---|---|
| barriera rettangolare (o tanh), orologi solo nella barriera | è l'apparato |
| ancella a macro-bin in fibra, birifrangenza alternata nelle guide | proposte mie, **non dimostrate** in laboratorio |
| θ = 0,5, 40 passi, σ = 4 siti; κ = 4,5 cm⁻¹, σ = 1,5 cm⁻¹, 12 cm | scelti da me: plausibili, non verificati come fattibili |
| modello a modi accoppiati delle guide (primi vicini, Longhi 2010) | modello dell'apparato, non la regola |
| rumore gaussiano bianco per sito e per passo | il rumore vero è colorato (deriva lenta della fibra è quasi innocua: fase globale) |
| precisione polarimetrica 10⁻³ | ipotesi |

**Parametri liberi: 0** (θ, V0, L, σ, δ sono impostazioni dell'esperimento).

## 7. Come rifare tutto (≈ 6 minuti su un core)

```bash
cd lab; export OMP_NUM_THREADS=1
python3 tunproto_prereg.py              # pre-registrazione (già committata) -> out/tunproto_prereg.json
python3 tunproto_rA.py regola           # A1, A2, A5 (~45 s)
python3 tunproto_rA.py dirac            # A3 split-step (~2,5 min)
python3 tunproto_rA.py bordi            # B1, B2 (~20 s)
python3 tunproto_rC.py pulito; python3 tunproto_rC.py mc; python3 tunproto_rC.py rumore
python3 tunproto_rC.py scan             # post-hoc: soglia di rumore (~10 s)
python3 tunproto_rD.py                  # guide d'onda (~2 s)
python3 tunproto_verdict.py             # giudizio -> out/tunproto_verdict.json
python3 tunproto_figure.py              # figura
```

## 8. Cose post-hoc e cambi (dichiarati)
1. Estrapolazione quadratica in A2 (l'errore della regola qui è O(a²)): solo come informazione; l'esame usa la lineare.
2. B1: la previsione aveva una coda di tanh tagliata; ricalcolata post-hoc con zona ±8w (accordo 2·10⁻⁷). L'esame resta FAIL.
3. Scansione del rumore (`tunproto_rC_scan.json`) non pre-registrata.
4. Nel progetto (prima della pre-registrazione) ho visto che un taglio del pacchetto a 3σ rovina le letture a N = 8;
   ho registrato 4σ, che basta per l'inclinazione ma non del tutto per il piano (C1).
5. La figura, pannello (a): banda = statistica ⊕ scarto tra pacchetto tagliato e previsione spettrale; pannello (b):
   lo zero di τ_my cade a V0/E ≈ 1,05 e non a 1 per il pacchetto finito (σ = 4 siti), effetto previsto dalla stessa
   regola.
