# Protocollo sperimentale per l'orologio di massa nel tunnel

Studio `tunproto_` (file: `lab/tunproto_*.py`, `lab/out/tunproto_*`). Segue `docs/TEMPO_TUNNEL.md` (giro 1).
Stato: **fase 1 (ricerca) scritta**; pre-registrazione, corse e risultati nelle sezioni successive.

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
