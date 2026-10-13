# Idrogeno vero a scala di Bohr: regola imparata lontano, zig-zag vicino al protone

**Studio:** `idrogeno_bohr` (PIANO §5 passo 2, richiesta 4b di Valerio).
**Pre-registrazione:** `lab/out/idrogeno_bohr_prereg.json`. È stata committata prima delle corse di produzione (commit 6327b86, con un cambio dichiarato in c84f198).
**Codice:**
- `src/atomo/idrogeno.rs`: il motore.
- `src/atomo/idrogeno_main.rs`: il binario headless `idrogeno`.
- `lab/idro_bohr_figura.py`: la figura.

**Numeri:** `lab/out/idrogeno_bohr_risultati.json`.
**Figura:** `lab/out/idrogeno_bohr_figure.png`.
**Dati per il video:** `/mnt/project-files/ripresa/dati/idrogeno_bohr/` (manifest.json).

---

## 0. In breve

- **Che cosa abbiamo costruito.** Un atomo di idrogeno vero (Z = 1, protone, elettrone, α vera) in 3-D, che evolve nel tempo e gira su CPU.
- **Come è fatto.**
  - Lontano dal nucleo lavora la **regola imparata dall'agglomeratore**, con celle da a0/16, cioè la regola grossa.
  - Vicino al protone c'è un **blocco annidato** in cui lavora la **regola zig-zag vera**: celle da a0/512, e la pallina corre a c e poi ruota.
  - In mezzo c'è una **colla**: un prolungamento conservativo e una media di blocco.
- **Risultato principale.** Tutti gli 8 esami pre-registrati passano.

| Esame | Che cosa misura | Valore | Soglia | Esito |
|---|---|---|---|---|
| E1 | energia 1s, ibrido (blocco 3³), scatola ±8 a0 | **−0,066 %** (−13,58955 eV) | ≤ 0,10 % | PASS |
| E2 | 2s e 2p (regola grossa a0/8, scatola ±18 a0) | 2s **−0,138 %**, 2p **+0,039 %**, tripletto aperto 3·10⁻⁵ | ≤ 0,25 %, ≤ 10⁻⁴ | PASS |
| E3 | 1s fermo nel tempo (12 → 24 passi) | densità grossa 4,5·10⁻⁴, blocco fine 9·10⁻⁴ | ≤ 10⁻³, ≤ 10⁻² | PASS |
| E4 | risonanza 1s→2p dal tempo (150 passi veri) | **10,1498 eV = 122,15 nm** (−0,48 %); coerenza con gli autovalori 2,5·10⁻⁸ | ≤ 0,6 %, ≤ 10⁻⁴ | PASS |
| E5 | costo | **3,8 s per passo grosso** (32 tick fini dentro), 0,7 s senza blocco; 2 thread | ≤ 10 s | PASS |
| E6 | quota della colla (pozzo liscio, blocco 3³) | 6,7·10⁻⁵ | ≤ 10⁻⁴ | PASS |
| E7 | nessuna riflessione (pacchetto che attraversa il blocco) | passa 0,99999998, riflesso 5·10⁻¹⁰ | ≥ 0,99, ≤ 10⁻² | PASS |
| E8 | carica ed energia attraverso la colla | carica 7,8·10⁻⁸ per passo, ⟨H⟩ 8·10⁻⁸, pacchetto 6·10⁻⁸ | ≤ 3·10⁻⁷, ≤ 10⁻⁵, ≤ 10⁻⁴ | PASS |

- **Validazione o scoperta?** È **validazione**: i numeri sono quelli noti di Dirac, e li ritroviamo.
- **Le cose nuove sono tre**, e stanno nella diagnosi:
  1. La scatola periodica è un **cristallo di atomi** che sposta i livelli. Fa +13 % sul 2p se la scatola è troppo stretta.
  2. Il cammino zig-zag su una cella "grande" (a_f = 0,27 lunghezze Compton) ha una **massa efficace sbagliata del 2,4 %**.
     - L'angolo giusto si **deriva**, non si fitta.
     - È stata la causa principale della debolezza della colla vista in `NUCLEO_1S.md`.
  3. L'errore residuo della regola si spiega in **due pezzi**:
     - la massa del cammino di base: **+0,037 %**, che si vede da sola nel 2p (+0,039 %);
     - un termine di contatto proporzionale a |ψ(0)|². Tocca solo gli stati s ed è la metà nel 2s rispetto all'1s, come previsto.
- **Che cosa non è venuto.** L'ibrido **non è ancora convergente** nella grandezza del blocco fine.
  - Con blocchi 1³, 3³ e 5³ l'1s dà −0,097 %, −0,066 % e −0,043 %.
  - La dispersione è 5,4·10⁻⁴, contro una colla perfetta che darebbe 0. Vedi §5.

---

## 1. L'idea in parole semplici

L'elettrone dell'idrogeno vive in una nuvola larga qualche a0. Il protone però è un punto, e lì la nuvola ha una punta (la "cuspide").

Se simulassimo tutto con la zig-zag vera alla risoluzione che serve al protone, ci vorrebbero (a0 ⋅ 16 / (a0/512))³ ≈ 10¹⁵ celle. Impossibile.

Quindi si fa come una macchina fotografica con lo zoom:

1. **Lontano: la regola grossa.**
   - Il cubo è diviso in celle da a0/16 (256³ celle, scatola ±8 a0).
   - Su ogni cella agisce la regola che l'agglomeratore ha imparato guardando la zig-zag (`agglomerato_zz.json`, catena A, livello B = 256).
   - Un passo della regola = un tempo τ = a0/16 / c. Dentro c'è già tutto: cinetica, massa, la fase di riposo.
2. **Vicino al protone: la zig-zag vera.**
   - Un cubetto di 3 × 3 × 3 celle grosse (lato 0,19 a0) è diviso ancora in 32 × 32 × 32 per cella. Sono 96³ celle fini da a0/512.
   - Lì l'elettrone è la pallina della zig-zag: in un tick corre di una cella a c, poi ruota. Il Coulomb è una fase a ogni tick.
   - In un passo grosso ci sono 32 tick fini.
3. **La colla, ad ogni passo grosso.**
   - Il bordo del cubetto fine (lo "alone", spesso una cella grossa) viene ricostruito dalla griglia grossa. La ricostruzione è conservativa: la media di ogni blocco fine ricostruito è esattamente il valore grosso (R∘P = I).
   - Dentro il cubetto la zig-zag evolve per 32 tick.
   - Poi le celle grosse coperte prendono la media dei blocchi fini.
   - Il resto della griglia fa un passo della regola grossa.
   - Il "calcio" del Coulomb nelle celle coperte viene dal fine. Così la regola grossa vede il nucleo come lo vede il fine.

Che cosa vede l'occhio:
- **da lontano**, la nuvola che oscilla, cioè la sovrapposizione 1s + 2p che "respira" alla frequenza di Lyman-alfa;
- **da vicino**, la cuspide sul protone e il tremolio delle palline che corrono a c avanti e indietro.

## 2. Che cosa è validazione e che cosa è nuovo

| Validazione (numeri noti ritrovati) | Nuovo (non c'era prima) |
|---|---|
| E(1s), E(2s), E(2p) entro 0,07–0,14 % di Dirac × massa ridotta | regola grossa + zig-zag annidata nello **stesso** atomo 3-D, dipendente dal tempo, in Rust su CPU |
| densità al nucleo dal blocco fine 0,3176 a0⁻³ (esatto 1/π = 0,3183) | la massa del cammino su celle grandi (θ cot θ) come causa principale dell'errore della colla, e l'angolo derivato che la toglie |
| elemento di matrice ⟨1s|z|2p⟩ = 0,7489 a0 (esatto 0,7449) | la scatola periodica come cristallo, misurata (1s: ±4 a0 → +1,7 %; 2p: ±8 a0 → +13 %) |
| Lyman-alfa 122,15 nm (Dirac × μ: 121,57 nm) | la decomposizione dell'errore: massa di base +0,037 % + contatto ∝ \|ψ(0)\|² |

## 3. Unità e confronto

- **Unità.** ħ = c = μ = 1, con μ la massa ridotta. Il problema di Dirac-Coulomb con massa ridotta è esattamente questo. Il rinculo entra solo nella conversione in eV.
- **Debito.** Il nucleo è fermo: è la massa ridotta, non un protone che si muove.
- **Riferimenti** (Dirac × μ):
  - 1s = −13,598468 eV;
  - 2s½ = −3,399628 eV;
  - baricentro del 2p (2·2p½ + 4·2p3/2)/6 = −3,399598 eV;
  - 1s→2p = 10,198870 eV.
- **Lo spin.** La regola grossa non ha spin: la parte di spin imparata è ~10⁻¹⁰ ed è scartata. Per questo la confronto con il baricentro.

## 4. Gli esami uno per uno

### E1 — energia dell'1s (PASS)

- **Ibrido**, blocco 3³, scatola ±8 a0, 24 passi. L'energia si legge dalla pendenza della fase di ⟨1s|ψ(t)⟩ dal passo 8.
  - **−13,589554 eV, −0,066 %.**
- **Solo grosso:** −13,587576 eV, −0,080 %.
- Il blocco fine sposta di +1,5·10⁻⁴ verso Dirac.

### E2 — stati n = 2 (PASS)

**Perché a0/8.** Gli stati n = 2 sono 4 volte più larghi dell'1s, e la scatola deve essere ±16–20 a0:
- a ±8 a0 il 2p esce **+13 %** (vedi §6.1);
- a a0/16 servirebbero 640³ celle, troppe per questa CPU.

Quindi si usa il livello imparato successivo, a0/8, con 288³ celle e scatola ±18 a0.

| stato | E regola (eV) | contro Dirac × μ | parità |
|---|---|---|---|
| 1s | −13,550798 | −0,351 % | +1,0000 |
| 2p_z / 2p_x / 2p_y | −3,401028 / −3,400878 / −3,400838 | +0,042 / +0,038 / +0,037 % | −0,9995 / −0,9943 / −0,9998 |
| 2p (media) | −3,400915 | **+0,039 %** | |
| 2s | −3,394925 | **−0,138 %** | +0,9936 |

- **Il tripletto 2p** è aperto di 3,3·10⁻⁵ relativo. È l'anisotropia del reticolo cubico imparato.
- **La parità** non è ±1 esatta: s e p si mescolano un po' (vedi i debiti). Il motivo è che la regola imparata ha una piccola parte dispari (|Im g| fino a 1,8·10⁻³).
- **2s sta sopra il 2p di 6,0 meV.**
  - In Dirac il 2s½ sta sotto il baricentro di 0,03 meV (il Lamb shift qui non c'è).
  - Quindi quei 6 meV sono errore della regola: il termine di contatto dello stato s, vedi §6.3.

### E3 — l'1s sta fermo (PASS)

Tra il passo 12 e il passo 24:
- la densità sulla griglia grossa cambia al massimo di 4,5·10⁻⁴ del suo massimo;
- nel blocco fine di 9,1·10⁻⁴.

Il blocco fine parte dalla ricostruzione liscia (senza cuspide) e la cuspide si forma da sola nei primi passi: è il transitorio.

### E4 — la risonanza 1s → 2p (PASS)

**Il metodo.**
- Lo stato (1s + 2p_z)/√2 viene evoluto con **150 passi veri** della regola grossa (a0/8, ±18 a0).
- ω si legge dalla pendenza della differenza di fase delle due proiezioni.

**I risultati.**
- ω = **10,149770 eV**, cioè λ = **122,15 nm**.
- Contro Dirac × μ (10,198870 eV): **−0,48 %**.
- Contro la differenza degli autovalori della regola: 2,5·10⁻⁸. La dinamica fa esattamente quello che dice lo spettro.
- La "cartolina" (i due autostati con le loro fasi) e lo stato evoluto davvero coincidono con fedeltà 1 − 10⁻¹⁰.
- La norma si conserva a 10⁻⁹.

**Da dove viene l'errore.** Il −0,48 % è tutto dell'1s a a0/8 (−0,35 %); vedi §6.3. Per avere la risonanza al livello dell'1s ibrido servirebbe il 2p a a0/16 in una scatola ±18 a0, cioè 576³ celle: non sta in 15 GB con il solutore attuale.

**Il periodo** della risonanza è 16,84 u.a. = 0,41 fs = 18 457 passi della regola a0/8.

### E5 — costo (PASS, su 2 thread)

Misurato su 2 thread, perché gli altri 2 core erano di altri agenti.

| configurazione | s per passo grosso |
|---|---|
| solo grosso a0/16, 256³ | 0,72 |
| ibrido, blocco 1³ (patch 96³) | 1,19 |
| **ibrido, blocco 3³ (patch 160³)** | **3,81** |
| ibrido, blocco 5³ (patch 224³) | 8,06 |
| solo grosso a0/8, 288³ | 1,5 |

- **Un tick fine** costa (3,81 − 0,72)/32 ≈ 0,1 s per 160³ × 4 componenti.
- **Su 4 thread** (stima, non misura) mi aspetto ~2,2 s per passo grosso: FFT e tick sono paralleli quasi perfetti.
- **Per il tempo atomico** bisogna fare i conti:
  - 1 u.a. = 2193 passi a0/16;
  - il periodo di Lyman-alfa è 36 900 passi a0/16.
  - Quindi il "film" lungo non si fa passo per passo: si fa con la cartolina (vedi §7).

### E6 — quanto pesa la colla (PASS sulla parte pre-registrata)

**(a) Pozzo liscio.** Al posto del Coulomb si mette un pozzo gaussiano (8 Ha, largo a0/2): senza cuspide, l'ibrido dovrebbe dare **esattamente** il grosso.
- Blocco 3³: E_ibrido/E_grosso − 1 = **−6,7·10⁻⁵** (soglia 10⁻⁴), PASS.
- Blocco 1³: −5·10⁻⁶.
- Blocco 5³ (sviluppo): −2,3·10⁻⁴.

**(b) Dispersione con la grandezza del blocco (Coulomb, ±8 a0).**

| blocco | 1s | ibrido/grosso − 1 |
|---|---|---|
| nessuno | −0,080 % | 0 |
| 1³ | −0,097 % | −1,7·10⁻⁴ |
| 3³ | −0,066 % | +1,5·10⁻⁴ |
| 5³ | −0,043 % | +3,7·10⁻⁴ |

Interpretazione onesta: §5.

### E7 — nessuna riflessione (PASS)

Un pacchetto libero (σ = 3 celle) attraversa il blocco fine 3³:

| impulso | passa (fedeltà con il solo grosso) | riflesso | carica |
|---|---|---|---|
| 0,5 / cella | 0,99999998 | 4,9·10⁻¹⁰ | −2,3·10⁻⁸ |
| 0,25 / cella | 1,00000005 | 2,5·10⁻⁹ | +6,2·10⁻⁸ |

Il criterio del reticolo era "passa ≥ 99 %": qui si sta a 10⁻⁸ dal 100 %.

### E8 — carica ed energia attraverso la colla (PASS)

Era richiesta del coordinatore, confrontata con `docs/MOTORE_OTTIMIZZATO.md` §6 D2. Lì:
- la media di blocco perde fino a 9·10⁻³;
- il prolungamento a gradini inventa +105…+3800 % di energia cinetica;
- quello spettrale conserva a 10⁻¹⁵.

Qui:
- **carica:** +7,8·10⁻⁸ per passo grosso con l'1s e blocco 3³ (−2,5·10⁻⁸ con blocco 1³);
- **⟨H⟩ della griglia grossa:** −7,7·10⁻⁸ nella seconda metà della corsa;
- **pacchetto che attraversa:** ≤ 6·10⁻⁸.

La nostra colla non inventa energia cinetica, perché il prolungamento è polinomiale di ordine 9 conservativo, non a gradini. Però non è esatta come quella spettrale.

Su un periodo di Lyman-alfa (36 900 passi) la deriva di carica arriverebbe a ~3·10⁻³. **Per i film lunghi serve la colla a flusso conservativo** (debito).

## 5. La colla: che cosa abbiamo capito e che cosa no

### Che cosa si è aggiustato

`NUCLEO_1S.md` §5.2 aveva una "debolezza della colla ~0,1 %". In sviluppo, con il cammino "nudo" (angolo di rotazione θ = a_f = 0,268):
- blocco 1³: −2,5·10⁻⁴;
- blocco 3³: +1,1·10⁻³.

**Causa trovata.** Un cammino che ruota di θ a ogni cella a ha, per piccoli impulsi, massa efficace m/(a cot θ). Con θ = a = 0,268 viene θ cot θ = 0,976: l'elettrone nel blocco è **2,4 % più pesante** che nella regola grossa. La regola grossa viene dalla catena con θ = 0,0335, quindi è più leggera di così.

**Rimedio (derivato, non fittato).**
- Si sceglie l'angolo in modo che a cot(angolo) = κ, dove κ è il coefficiente cinetico della regola grossa, letto da G(k) per k → 0.
- Si rimette a posto la fase di riposo con un fattore uniforme e^{−i(a − angolo)} per tick.
- Il cammino resta esatto: confrontato con la matrice di Bloch a 1·10⁻¹⁵.

**Effetto.** La dispersione blocco 1/blocco 3 scende da 1,4·10⁻³ a 3,2·10⁻⁴.

### Che cosa resta

1. **Il residuo della colla con un pozzo liscio.** Vale −5·10⁻⁶, −6,7·10⁻⁵ e −2,3·10⁻⁴ per blocchi 1, 3 e 5. Cresce con il blocco, e va nel verso sbagliato rispetto alla soluzione esatta del pozzo (−42,6416 eV): il grosso dà −42,6133, l'ibrido 3³ −42,6104.
   - Anche un pacchetto **fermo e libero** nel blocco perde carica: 3·10⁻⁷, 2,5·10⁻⁶ e 6·10⁻⁶ per passo per blocchi 1, 3 e 5.
   - Alzare l'ordine della ricostruzione da 5 a 9 toglie il 40 %; 13 non toglie di più.
   - Quindi **non** è la ricostruzione. È che la colla non fa combaciare i **flussi**:
     - la regola grossa (larga 9 celle) muove carica dentro e fuori dalle celle coperte;
     - la zig-zag nel blocco ne muove un'altra;
     - nessuno le pareggia.
   - È lo stesso problema che negli adattivi si risolve con il "reflux".
2. **La regola grossa e la zig-zag non hanno la stessa dispersione a k grandi.** *(Corretto in fase 2, §11.)* In fase 1 avevo scritto "4 volte la correzione relativistica": era sbagliato. Gran parte di quell'effetto era la **parte dispari** della regola, perché G(k) lungo +x non è uguale a G(−k). Tolta la parte dispari, resta un errore vero solo ad alto k:
   - regola vecchia simmetrizzata: −5·10⁻³ a k = 0,8 per cella lungo l'asse;
   - regola T: −3,8·10⁻³.

   È l'errore ad alto k ancora aperto in `AGGLOMERATO_COMPLETO.md`.
3. **Il blocco fine non è convergente in grandezza** (Coulomb: −0,097, −0,066, −0,043 % per blocchi 1, 3, 5).
   - Il blocco fine corregge davvero la fisica vicino al nucleo: la forma vicino al protone torna entro 1,5·10⁻⁴ dalla media esatta, contro il 2·10⁻³ del solo grosso.
   - Però la parte di errore che la regola grossa fa vicino al nucleo è sparsa su ≥ 5 celle (la regola è larga ±4 celle).
   - **La quota della colla** sull'errore 3-D dell'ibrido 3³ è ≈ 7·10⁻⁵ (dal pozzo liscio) su un errore totale di 6,6·10⁻⁴: circa il 10 %.
   - **Il resto** è la regola grossa fuori dal blocco: il contatto e la massa di base, §6.3.

## 6. Le tre cose trovate in sviluppo

### 6.1 La scatola è un cristallo

La griglia è periodica. L'atomo vede le sue copie a distanza 2L, e le code delle nuvole si toccano: è un piccolo cristallo di idrogeno.

| | ±4 a0 | ±6 a0 | ±8 a0 | ±12 a0 | ±18 a0 |
|---|---|---|---|---|---|
| 1s (a0/16) | +1,70 % | −0,035 % | −0,080 % | | |
| 2p (a0/8) | | | +13,3 % | | +0,042 % |
| 2s (a0/8) | | | −5 % (non conv.) | +1,0 % | −0,138 % |

- Anche il codice Python (`agg_zz_atom`) dà −13,8366 eV a ±4 a0: non è un errore del Rust.
- Il "−0,044 %" di uno studio precedente a ±6 a0 era in parte aiutato da questo effetto.
- **Regola pratica:** scatola ≥ ±8 a0 per n = 1 e ≥ ±18 a0 per n = 2.

### 6.2 La massa del cammino dipende dalla cella

Vale m_eff = m/(θ cot θ), con θ = cella/lunghezza Compton:

| cella | θ cot θ | massa in più |
|---|---|---|
| catena A (base, a0/4096) | 0,999627 | +0,037 % |
| blocco fine (a0/512) | 0,976 | +2,4 % |

È fisica nota: la massa del reticolo di Dirac rinormalizzata. Qui però si vede in numeri dentro il nostro motore.

### 6.3 L'errore della regola in due pezzi

- **Pezzo 1, massa di base.** La regola imparata eredita la massa del cammino di base: +0,037 % di legame in più su **ogni** stato.
  - Il **2p**, che al nucleo ha densità zero, sbaglia di **+0,039 %**: la massa di base, e quasi nient'altro.
- **Pezzo 2, contatto.** Gli stati s hanno in più un termine di contatto ∝ |ψ(0)|². Relativo all'energia scala come 1/n (|ψ(0)|² ∝ 1/n³, E ∝ 1/n²). Quindi:
  - 1s a a0/8: −0,351 % = +0,037 (massa) − 0,388 (contatto);
  - 2s previsto: +0,037 − 0,388/2 = −0,157 %; misurato **−0,138 %**. Entro il 12 %, con dentro il mescolamento s–p.
- **Scala con h².** Il contatto a a0/16 dovrebbe essere 1/4: −0,097 %, quindi 1s ≈ −0,06 %. Il solo grosso a ±8 a0 dà −0,080 %.

## 7. Dati per il video (`/mnt/project-files/ripresa/dati/idrogeno_bohr/`, 273 MB)

Formato float16 LE come `atomo_accoppiato`, descritto in `manifest.json`.

- **Vista lontana** (`rho_lontano.f16` [24, 96, 96, 96], `j_lontano.f16` [24, 3, 96, 96, 96], `dipolo_lontano.txt`):
  - (1s + 2p_z)/√2 su un periodo di Lyman-alfa, 24 fotogrammi, griglia a0/8, ±6 a0.
  - La nuvola sale e scende lungo z; il dipolo oscilla con ampiezza 0,749 a0.
  - È fatta con la cartolina: i due autostati della regola con le loro energie. Il 2p lungo z è ricostruito dal tripletto, perché il solutore ha trovato tre 2p ruotati.
  - È esatta per la regola, ed è controllata contro 150 passi veri.
- **Zoom sul protone** (`rho_fine.f16` e `chir_fine.f16` [32, 96, 96, 96]):
  - Blocco fine, cella a0/512, lato 0,19 a0, ogni tick per un passo grosso (1 tick = 3,4·10⁻²² s; 32 tick = 1,1·10⁻²⁰ s).
  - Al centro la densità vale 0,3176 a0⁻³, contro l'esatto 1/π.
  - `chir` = |R|² − |L|²: quanta pallina corre "su" e quanta "giù" a c. Cambia a ogni tick: è il tremolio della zig-zag.
- **Griglia grossa dello stesso 1s** (`rho_grosso.f16`, `j_grosso.f16`, 64³ a a0/16).
  - Attenzione: la corrente dovrebbe essere zero. Il residuo (0,0065 u.a.) viene dalla parte dispari della regola imparata. È un debito.

## 8. Fallimenti, debiti, manopole

**Fallimenti e limiti dichiarati**
1. L'ibrido non converge nella grandezza del blocco: dispersione 5,4·10⁻⁴ fra blocchi 1, 3 e 5 (§5).
2. Colla non a flusso conservativo: la carica deriva di ~10⁻⁷ per passo. Su un periodo di Lyman-alfa fa ~3·10⁻³.
3. Gli stati n = 2 e la risonanza sono a a0/8, non con il blocco fine. La risonanza sbaglia di −0,48 % per colpa dell'1s grossolano a quel livello.
4. Il "passo a blocchi" (m passi liberi esatti + Coulomb a metà, per saltare tempo) costa troppo in precisione. A m = 64 l'1s sbaglia di 4,7·10⁻³. Quindi il film lungo si fa con la cartolina, non con la dinamica vera.
5. La soglia di E5 (costo) è stata misurata su 2 thread, non 4: vincolo della sessione.

**Debiti** (cose messe a mano o non emerse)
- Nucleo fermo: massa ridotta, non protone dinamico.
- Spin della regola grossa scartato: lontano non c'è struttura fine.
- Parte dispari della regola imparata: mescola s e p (parità 0,994) e dà una corrente spuria. Viene dall'apprendimento, non dalla fisica.
- Dispersione della regola imparata a k grandi 4 volte la relativistica (§5): da segnalare all'agglomeratore.
- Scelte di configurazione (non fit):
  - a_f = h/32;
  - blocco 3³ e alone di 1 cella;
  - ordine di ricostruzione 9;
  - scatole ±8 a0 e ±18 a0.
- L'angolo del cammino fine è **derivato** dalla regola grossa, con la fase di riposo uniforme. Non è un parametro, ma è un pezzo di colla scritto a mano.
- La vista lontana usa la cartolina (autostati), non 18 457 passi veri.

**Manopole fittate: 0.** Gli ingressi sono α, m_p/m_e, Z = 1 e la regola imparata.

**Cambi dopo la pre-registrazione**
1. κ estrapolato a k → 0 dai due k più piccoli (dichiarato nel file, prima delle corse).
2. La corsa P4 (spettro) è partita prima del commit della pre-registrazione: al commit avevo letto solo la riga dell'1s.
3. L'export lontano è stato rifatto con il 2p ricostruito lungo z (il primo export aveva il "2p_z" del solutore, che era ruotato: dipolo 0,15 invece di 0,75). Le energie e l'esame non cambiano.

## 9. Dove l'intuizione di Valerio può aiutare

1. **"Lo zoom deve far vedere la stessa cosa."**
   - A Valerio piace l'idea che da vicino e da lontano si veda la stessa pallina.
   - Qui abbiamo trovato che non è vero gratis: la pallina vicino al protone era il 2,4 % più pesante.
   - Un test visivo semplice: la stessa nuvola vista con il blocco fine acceso e spento non deve "saltare". Diventa un esame per ogni zoom futuro.
2. **Il cristallo.** L'atomo che vede le sue copie è proprio l'inizio di un solido.
   - Valerio potrebbe chiedere: "e se le copie le metto apposta più vicine?".
   - È la strada verso H₂ e i cristalli di idrogeno, con la stessa regola, senza aggiungere niente.
3. **Il flusso che non combacia.** Un'immagine aiuterebbe: palline che passano dalla zona grossa a quella fine come persone tra due stanze, dove ognuno conta chi entra.
   - Oggi le due stanze contano in modo diverso, e qualcuno si perde (10⁻⁷ per passo).
   - La regola giusta è che chi esce da una stanza entra nell'altra: lo stesso conteggio da entrambe le parti, ed è il "reflux".
4. **Il film lungo.** Per vedere l'atomo "respirare" servono 18 000–37 000 passi per periodo.
   - Valerio ha l'intuizione del rallenty: la cartolina è proprio questo. Si calcolano due stati fermi e si fa scorrere solo la loro fase.
   - Da chiedere a lui: ci basta questo per i video, o vuole vedere la dinamica vera con una luce che arriva? Nel secondo caso serve la colla conservativa.

## 10. Prossimi passi

1. **Colla a flusso conservativo** (reflux) tra la regola grossa larga 9 celle e il blocco fine. Esame: pacchetto fermo libero con carica ≤ 10⁻¹⁰ per passo e pozzo liscio ≤ 10⁻⁶.
2. **Blocco fine più largo e più grosso** (a_f = a0/256, blocco 7³–9³) per arrivare alla convergenza dell'1s ibrido. Ci si aspetta di andare verso +0,037 %, il limite della massa di base.
3. **Agglomeratore:** segnalare la dispersione a k grandi (4 volte la relativistica) e la parte dispari |Im g| ~ 3·10⁻³. Con una regola simmetrica per parità spariscono il mescolamento s–p e la corrente spuria.
4. **2p a a0/16** con un solutore a blocchi su griglia composita (a0/8 fuori, a0/16 dentro), per portare la risonanza dove è l'1s.
5. **Collegare al Laboratorio Atomo** (grafica): usare `idrogeno.rs` come motore e i dati esportati per il primo video.

## Come rifarlo

```
cargo build --release -j 2 --bin idrogeno
B=target/release/idrogeno; export RAYON_NUM_THREADS=2
$B spettro 512 288 0 /tmp/stati                       # E2 (≈ 22 min)
IDRO_ORD=9 $B esame_ibrido 256 256 3 24 /tmp/zoom      # E1 E3 E5 E8 + export zoom (≈ 3 min)
IDRO_ORD=9 $B esame_ibrido 256 256 1 24                # E6b
IDRO_ORD=9 $B esame_ibrido 256 256 5 16                # E6b
IDRO_ORD=9 IDRO_POT=gauss $B esame_ibrido 256 128 3 16 # E6a (esatto del pozzo: −42,641613 eV, solutore radiale)
IDRO_ORD=9 $B attraversa 64 3 0.5 3 -20 80             # E7
IDRO_ORD=9 $B attraversa 64 3 0.25 3 -20 160           # E7
$B esame_risonanza 512 288 /tmp/stati 150 0 24 /tmp/lontano   # E4 + export lontano
python3 lab/idro_bohr_figura.py
```

Controlli interni: `$B testfine` (R∘P = I esatto, cammino contro matrice di Bloch 1·10⁻¹⁵), `$B fft`, `$B massa` (massa efficace del cammino contro θ cot θ).


---

## 11. Fase 2: togliere le imprecisioni (13/10)

**Pre-registrazione:** `lab/out/idrogeno_bohr_fase2_prereg.json`.
- Primo commit 714b17f.
- Emendamenti 1–3 in 2111160, fe4487d, f37da41, tutti committati prima delle loro corse.

**Numeri:** `lab/out/idrogeno_bohr_risultati.json` → `fase2`. **Figura:** `lab/out/idrogeno_bohr_fase2_figure.png`.

### 11.1 Colla che pareggia i flussi ("reflux"): PASS

**L'idea.** C'è una sola porta fra le due stanze, e la si conta da una parte sola.
- La carica che entra nella zona fine o ne esce è quella mossa dalla regola grossa.
- La zig-zag dentro decide solo **dove** sta: la cuspide, il tremolio.
- In pratica, a ogni passo le celle coperte ricevono in totale Q(prima) − Q(celle scoperte, dopo), e il blocco fine è riscalato dello stesso fattore reale.
- Nessun parametro.

| esame | soglia | misura |
|---|---|---|
| R1a pacchetto fermo, deriva per passo | ≤ 10⁻¹¹ | **4·10⁻¹⁶** (fase 1: 2,5·10⁻⁶) |
| R1b 1s ibrido, deriva per passo | ≤ 10⁻¹¹ | **1,9·10⁻¹⁶**; E1 −0,0666 % |
| R1c pacchetto che attraversa | fedeltà ≥ 0,99, riflesso ≤ 10⁻² | 1,00000000, 2,9·10⁻¹⁰, ΔQ 5·10⁻¹⁵ |
| R1d pozzo liscio | ≤ 10⁻⁴ | 6,5·10⁻⁵ |

**Conseguenza.** Su un periodo di Lyman-alfa (36 900 passi) la deriva di carica è ~10⁻¹¹, cioè ≪ 10⁻⁴. Il debito "colla non conservativa" è **chiuso**.

**Limite.** Il reflux non tocca la fase. Un pacchetto fermo accumula comunque −2·10⁻⁶ rad in 8 passi, perché regola grossa e cammino non hanno la stessa dispersione ad alto k.

### 11.2 Convergenza nella grandezza del blocco: FAIL (come registrato)

Celle fini a0/256, reflux, scatola ±8 a0. L'1s contro Dirac × μ (%):

| blocco | 3 | 5 | 7 | 9 | 11 | 13 |
|---|---|---|---|---|---|---|
| regola vecchia (72 passi) | −0,0614 | −0,0280 | −0,0158 | −0,0381 | | |
| regola T (48 passi) | −0,0675 | −0,0489 | −0,0212 | −0,0380 | −0,0469 | −0,0409 |

- **Il criterio pre-registrato è mancato** (|E5−E3| e |E7−E5| < 10⁻⁴):
  - regola vecchia: 3,3·10⁻⁴ e 1,2·10⁻⁴;
  - regola T: 1,9·10⁻⁴ e 2,8·10⁻⁴.
- **Il blocco 7 è un picco**, con tutte e due le regole. Dal blocco 9 in poi i valori stanno fra −0,047 e −0,038 %. Fra 11 e 13 la differenza è 6·10⁻⁵, e il criterio "a blocco grande" dell'emendamento 3 passa, ma di poco.
- **La stima onesta dell'1s ibrido** è **−0,043 ± 0,004 %**. È il valore stabile per blocchi da 9 a 13, non un fit.
- **La lettura dell'energia non è il problema.** Su 72 passi le finestre [8,72], [24,72], [36,72] e [48,72] concordano entro 10⁻⁵.
- **Che cosa non capisco ancora.**
  - Il picco al blocco 7 non dipende dalla regola, quindi viene dalla geometria della colla: per esempio l'alone che cade sul bordo dello stencil della regola, largo ±4 celle. Non l'ho dimostrato.
  - Il valore a blocco grande (−0,04 %) è ancora lontano da +0,037 %, il limite atteso con la sola massa del cammino di base.
  - Il solo grosso cambia di 4,5·10⁻⁴ fra scatola ±6 e ±8 a0, quindi la scatola o la regola ad alto k (§11.3) pesano ancora a questo livello.

### 11.3 Regola con parità esatta

- **Il mio tentativo** (simmetrizzazione O_h, `assets/atomo/agglomerato_zz_simm_v1.json`, `lab/idro_simm_regola.py`) è venuto (R3a):
  - parità 1 − |P| ≤ 3·10⁻¹⁴;
  - tripletto 2p aperto 4·10⁻¹³.

  Su richiesta del coordinatore è stato **superato** dalla regola dell'agglomeratore completo.
- **Regola T** (emendamento 2):
  - È il metodo di `AGGLOMERATO_COMPLETO.md` §11 (learner "T": massa e k⁴ fissati al cammino fine, parità esatta) applicato ai livelli dell'idrogeno, a0/16 e a0/8, imparati direttamente dal cammino fine.
  - Script `lab/idro_regola_T.py`, nuovo file `assets/atomo/agglomerato_zz_T_v1.json`.
  - `agglomerato_zz.json` non è stato toccato.

| R4a (a0/8, ±18 a0) | soglia | regola vecchia | regola T |
|---|---|---|---|
| 1 − \|P\| (mescolamento s–p), peggiore | ≤ 10⁻⁸ | 6,4·10⁻³ (2s) | **2,2·10⁻⁹** PASS |
| tripletto 2p aperto | ≤ 10⁻⁷ | 3,3·10⁻⁵ | **6,3·10⁻⁸** PASS |
| 2p | ≤ 0,06 % | +0,039 % | **+0,034 %** PASS |
| 2s | ≤ 0,25 % | −0,138 % | **−0,163 %** PASS |
| 1s a a0/8 | (riportato) | −0,351 % | −0,385 % |
| 1s, solo grosso a0/16 | (riportato) | −0,080 % | −0,075 % |

**Una correzione al messaggio del coordinatore.** La regola T ha la massa uguale a quella del **cammino fine** (1,000373), non a quella di Dirac. Il learner fissa la massa al "maestro".
- Il +0,037 % della fase 1 è proprio la massa del cammino fine a θ = a0/4096 (θ²/3), e **resta**.
- Si vede nel 2p, che sta a +0,034 %.
- Per toglierlo serve un cammino di base più fine, oppure la regola della fase 4 dell'agglomeratore, con massa esatta a 10⁻¹³. Quella non l'ho adottata, per non rifare corse già finite: la propongo come prossimo emendamento.

**R4c** (regola T):
- E4: ω = 10,14535 eV = **122,21 nm**. È −0,52 % da Dirac × μ (soglia 0,6 %, PASS); la coerenza con gli autovalori è 2,4·10⁻⁸.
- Dipolo ⟨1s|z|2p⟩ = 0,7495 a0 (esatto 0,7449).
- E1 con blocco 7: −0,0216 % (PASS ≤ 0,10 %). Il blocco 7 era stato scelto come "migliore" prima di vedere 11 e 13: è il picco, e il valore onesto è quello di §11.2.
- **Corrente spuria** ≤ 10⁻⁶: **FAIL**, misurata 1,7·10⁻³ (fase 1: 6,5·10⁻³).
  - La regola ora è pari e reale, quindi la corrente non viene più da lì.
  - Viene dal transitorio dell'ibrido: lo stato iniziale è l'autostato grosso, non quello dell'ibrido, e nei primi passi "respira".
  - L'avevo attribuita tutta alla regola: era sbagliato. Dichiarato.

**Errore ad alto k della regola imparata** (−3,8·10⁻³ a k = 0,8 per cella lungo l'asse): **FAIL dichiarato**, finché `AGGLOMERATO_COMPLETO` non lo chiude.

### 11.4 Dati per il video (fase 2)

`/mnt/project-files/ripresa/dati/idrogeno_bohr/fase2/` (250 MB). Le voci sono in `manifest.json` → `fase2`.
- Vista lontana con la regola T: 24 fotogrammi su un periodo.
- Zoom: blocco 7³ a a0/256, 112³ celle, 16 tick.

### 11.5 Conto

**Esami**
- R1 4/4 PASS.
- R2, R2bis, R4b FAIL.
- Emendamento 3 PASS di poco.
- R3a PASS (superato).
- R4a 4/4 PASS.
- R4c: E1 e E4 PASS, corrente spuria FAIL.
- Alto k FAIL dichiarato.

**Manopole fittate: 0.** I file nuovi:
- `agglomerato_zz_simm_v1.json` e `agglomerato_zz_T_v1.json` (asset versionati);
- `lab/idro_simm_regola.py`, `lab/idro_regola_T.py`.

**Prossimo**
1. Capire il picco al blocco 7: un alone di 2 celle al posto di 1, e lo stesso con la regola larga w = 2.
2. Partire dall'autostato dell'ibrido invece che da quello grosso, per togliere il transitorio e la corrente.
3. Adottare la regola della fase 4 dell'agglomeratore, con massa esatta: deve far sparire il +0,034 % del 2p.

## 12. Fase 3: il picco del blocco 7, l'autostato dell'ibrido, la regola della fase 4 (13/10)

Pre-registrazione: `lab/out/idrogeno_bohr_fase3_prereg.json`, committata **prima** delle corse (commit 79fa535).
Emendamento 1 (commit ce8660f) committato prima della corsa di controllo: aggiunge la diagnostica della corrente e il
controllo con la sola regola grossa. Il criterio di convergenza è quello della fase 2, **non allentato**.

Codice nuovo:
- `IDRO_NH`: lo spessore dell'alone, in celle grosse.
- `IDRO_FILTRO=K`: un filtro spettrale. Fa K passi dell'ibrido e somma lo stato intero (campo grosso + 4 componenti fini),
  pesato con una finestra di Hann e "girato all'indietro" con la fase dell'1s. Si riparte da questa somma. Quello che
  non gira alla frequenza dell'1s si cancella. Manopola: K = 48.
- `lab/idro_regola_P4.py` → `assets/atomo/agglomerato_zz_P4_v1.json`.

### 12.1 Il picco del blocco 7 (P1, P1b)

Configurazione: errore dell'1s rispetto a Dirac×μ, regola T, celle fini a0/256, reflux, 48 passi, finestra [24,48].

| blocco | fase 2 (alone 1) | P1: alone 2 | P1b: alone 1, partenza filtrata |
|---|---|---|---|
| 3 | −0,0675 % | −0,0648 % | −0,0674 % |
| 5 | −0,0489 % | −0,0403 % | −0,0490 % |
| 7 | −0,0212 % | **−0,0071 %** | −0,0214 % |
| 9 | −0,0380 % | −0,0135 % | −0,0371 % |

- **P1b FAIL** su convergenza e picco. Partire dall'autostato dell'ibrido **non cambia nulla**: le energie sono uguali
  alla fase 2 entro 1·10⁻⁵. Quindi il picco non viene dallo stato di partenza.
- **P1 FAIL** su convergenza (|E5−E3| = 2,5·10⁻⁴) e picco (E7 − (E5+E9)/2 = 2,0·10⁻⁴; soglie 1·10⁻⁴).
  Però **l'alone conta**: i blocchi 7 e 9 si spostano di 1,4·10⁻⁴ e 2,5·10⁻⁴.
- Detto semplice: quanto vale l'energia dipende **da come cuciamo** il blocco fine alla regola grossa (alone e
  ricostruzione), più che da quanto è grande il blocco. Finché la cucitura conta per 10⁻⁴, il blocco non può
  "convergere" a 10⁻⁴.
- Il valore vero dell'1s resta incerto fra −0,01 % e −0,05 %. Tutti i numeri della serie stanno **dentro 0,07 %**.
- La prova con la regola stretta w = 2, scritta come idea nella §11.5, **non è stata fatta**: era fuori dalla
  pre-registrazione e il budget era stretto.

### 12.2 L'autostato dell'ibrido e la corrente spuria (P2 + emendamento 1)

Configurazione: blocco 7, alone 1, partenza filtrata.

| esame | risultato | soglia | esito |
|---|---|---|---|
| E1 | −0,0214 % | 0,10 % | PASS |
| stabilità della lettura (finestre [8,48]…[32,48]) | 5·10⁻⁷ | 1·10⁻⁵ | PASS |
| stazionarietà grosso / fine | 3,0·10⁻⁴ / 5,3·10⁻⁴ (fase 2 ~1·10⁻³) | 1·10⁻⁴ | FAIL |
| corrente spuria max\|j\| | 1,64·10⁻³ α c/a0³ | 1·10⁻⁶ | FAIL |

**Cosa si capisce.**
- Il filtro sposta lo stato di appena 1,7·10⁻⁹ in probabilità. L'ibrido partiva già praticamente dal suo autostato.
- La corrente **non è un transitorio**, come invece avevo ipotizzato in fase 2: correzione dichiarata.
- La diagnostica dell'emendamento 1 mostra dove sta la corrente, per zone:

  | zona | P2 | fase 2 |
  |---|---|---|
  | celle grosse coperte dal blocco | 1,66·10⁻³ | 1,72·10⁻³ |
  | appena fuori dal blocco | 2,7·10⁻⁴ | 3,3·10⁻⁴ |
  | oltre 0,5 a0 | 1,3·10⁻⁴ | 2,3·10⁻⁴ |

- Nelle celle coperte il campo grosso è solo la media del campo fine. La sua "corrente", letta con la regola grossa,
  non è un flusso vero: è il segno della **cucitura**, la stessa che sposta l'energia nella §12.1.
- Era radiale, verso il nucleo, e in uno stato stazionario non può esserci un flusso netto: con la carica conservata
  a 10⁻¹⁵ per passo, non scorre niente davvero.

**Controllo con la sola regola grossa** (emendamento 1, nessun blocco, 48 passi):

| misura | sola regola grossa | ibrido |
|---|---|---|
| corrente | 2,9·10⁻⁵ | 1,64·10⁻³ |
| stazionarietà | 2,2·10⁻⁵ | 3,0·10⁻⁴ |

- La corrente della sola regola grossa è il **pavimento** della misura, che viene dal residuo del solutore.
- Quindi la soglia 1·10⁻⁶ della fase 2 era sotto il pavimento di circa 30 volte: va dichiarato.
- L'ibrido però sta **57 volte sopra** il pavimento, quindi il FAIL è vero e resta.

### 12.3 La regola della fase 4 dell'agglomeratore (P3)

**Costruzione.** `assets/atomo/agglomerato_zz_P4_v1.json` (nuovo, versionato) rifà il metodo di
AGGLOMERATO_COMPLETO §12:
- c2, c4 e c6 fissati al livello sotto, parità esatta;
- B = 512 imparato dal B = 256;
- passa il controllo kcheck, con w = 4 a entrambi i livelli.

Massa assoluta: 1,000373 = 1 + θ²/3 del cammino fine.

**Previsione scritta prima delle corse.** La "massa esatta" della fase 4 è esatta rispetto al **maestro** (il cammino
fine), non rispetto a Dirac. Quindi il +0,034 % del 2p **non sparisce**.
- Esce **+0,0338 %**: previsione giusta.
- Correggo la §11.5 punto 3, che diceva il contrario.

**Regola.**

| controllo | B = 256 | B = 512 |
|---|---|---|
| c6 (vincolo) / Dirac | 1,00018 | 0,706 |
| banda / cammino a k = 0,785 per cella sull'asse | +1,3·10⁻³ | +1,5·10⁻³ |

L'errore di banda della regola T al primo livello era −3,8·10⁻³.

**Spettro a a0/8** (`spettro 512 288`):

| esame | regola P4 | regola T | soglia | esito |
|---|---|---|---|---|
| parità 1−\|P\| | ≤ 5·10⁻⁷ (stampa a 6 cifre) | | ≤ 1·10⁻⁶ | PASS |
| allargamento del tripletto 2p | 8,5·10⁻⁸ | | ≤ 1·10⁻⁶ | PASS |
| 2p | +0,0338 % | | +0,029…+0,039 % | PASS |
| 2s | −0,1615 % | −0,163 % | | |
| 1s(a0/8) | −0,383 % | −0,385 % | | |
| Lyman-α dagli autovalori | 10,1457 eV (−0,52 %) | | | |

**Ibrido** al blocco 7 con la regola P4 e la partenza filtrata:
- E1 **−0,0276 %**: PASS entro 0,10 %, contro −0,0214 % con la regola T;
- corrente 1,62·10⁻³ e stazionarietà 2,7·10⁻⁴: FAIL, come P2.

La regola della fase 4 è un po' migliore ad alto k, ma **non cambia** gli stati legati oltre 0,01 %.

### 12.4 Conto della fase 3

**Esami.**
- **FAIL:** P1 (convergenza e picco), P1b (convergenza e picco), P2 su stazionarietà e corrente.
- **PASS:** P2 su energia e lettura, P3 su costruzione, parità, tripletto, 2p (previsione giusta) e ibrido.

**Manopole.** Una nuova: K = 48 del filtro, scelta prima. 0 manopole fittate.

**Debiti aperti.**
- **La cucitura.** Pesa ~1–3·10⁻⁴ sull'energia e dà la corrente nelle celle coperte. È il prossimo cantiere vero:
  serve una restrizione/prolungamento che rispetti anche la **fase** (la corrente), non solo la carica come fa il reflux.
- **La massa del maestro.** Il +0,037 % viene dal cammino fine a a0/4096, θ²/3. Si toglie solo con un cammino più fine
  o con la massa corretta nel maestro, non con la regola grossa.

**Come rifarlo.**
```
IDRO_ASSET=assets/atomo/agglomerato_zz_T_v1.json IDRO_REFLUX=1 IDRO_ORD=9 IDRO_B=16 IDRO_NH=2 target/release/idrogeno esame_ibrido 256 256 7 48
IDRO_ASSET=... IDRO_REFLUX=1 IDRO_ORD=9 IDRO_B=16 IDRO_FILTRO=48 target/release/idrogeno esame_ibrido 256 256 7 48 <dir>
python3 lab/idro_regola_P4.py; IDRO_ASSET=assets/atomo/agglomerato_zz_P4_v1.json target/release/idrogeno spettro 512 288 0 <dir>
```

**Video.** Nessun nuovo dato: i fotogrammi della fase 2 restano validi, perché l'ibrido filtrato è uguale entro 10⁻⁵.

## 13. Fase 4: la cucitura, la massa del maestro, il grosso di default (13/10)

Pre-registrazione: `lab/out/idrogeno_bohr_fase4_prereg.json`, committata **prima** delle corse di produzione
(commit 9be2ca3), con le corse di sviluppo già fatte scritte dentro. Emendamento 1 (commit 9bf3421): solo il salvataggio e
la ripresa delle corse lunghe, per stare nei 10 minuti per script; soglie invariate.

Tre lavori:
1. **la cucitura** fra il blocco fine (zig-zag) e la regola grossa;
2. **la massa del maestro** (il +0,034 % del 2p), da sistemare con l'agglomeratore e non a mano;
3. **il grosso di default**: regola T + bordi ritorti (0, π) + coda di Born sulla cella del nucleo, come emendamento
   dichiarato (i numeri del filone debiti, rifatti con la regola nuova).

Figura: `lab/out/idrogeno_bohr_fase4_figura.png`.

### 13.1 La cucitura

**Il bug trovato.** Il blocco fine prende il suo bordo (l'alone) dalla regola grossa. Per la componente piccola del
campo (χ = Σ G_j(−i∂_j)φ) serve una derivata. All'ultima cella grossa del patch il codice faceva la differenza centrata
con l'indice "bloccato" al bordo: in pratica **metà derivata**. Quindi al bordo dell'alone la corrente era sbagliata.
- Effetto: il blocco **perdeva carica** di ~5·10⁻⁶ per passo, sempre nello stesso verso. È proprio la perdita che in
  fase 2 avevamo tappato col reflux.
- Correzione: al bordo una differenza unilatera del secondo ordine, (−3, 4, −1)/2. Zero parametri.
  (Il vecchio comportamento resta con `IDRO_BORDO_VECCHIO=1`.)
- Con la correzione la perdita scende a **3·10⁻⁸ per passo** (mediana). Il reflux resta acceso ma non ha quasi più niente
  da fare.

**La partenza.** Partire dall'1s grosso vuol dire partire con un po' di stati eccitati mescolati, perché l'ibrido è
più legato dell'1s grosso. Questi battono su tempi di ~1 unità atomica, cioè migliaia di passi: un filtro corto non può
toglierli (il filtro della fase 3 toglieva solo roba a centinaia di Hartree).
Nuova partenza (`IDRO_PARTENZA=veff:M:R`), **senza regola fittata**:
- si fanno M passi dell'ibrido;
- su ogni cella coperta dal blocco si misura quanta fase l'ibrido aggiunge rispetto al passo grosso semplice: è un
  potenziale locale V_c (sul nucleo +50 eV, cioè una cella sola che "spinge" diversamente);
- si prende lo stato più basso della regola grossa + V + V_c, e si ripete R volte (converge in 2 giri: 1−|⟨prima|dopo⟩|² = 3·10⁻¹⁴).
Cambia **solo lo stato di partenza**: dopo, la dinamica è l'ibrido senza aggiunte. M = 32, R = 2–3, dichiarati.

**Esami della cucitura** (N = 128, ±4 a0, blocco 7, 240 passi; pavimento = la stessa misura con la sola regola grossa):

| esame | fase 3 | fase 4 | pavimento | soglia | esito |
|---|---|---|---|---|---|
| C1 perdita di carica per passo (mediana) | ~5·10⁻⁶ | **3,1·10⁻⁸** | 0 | ≤ 2·10⁻⁷ | PASS |
| C2 stazionarietà delle celle grosse fuori dal blocco | 3·10⁻⁴ (tutte) | **1,7·10⁻⁵** | 1,9·10⁻⁶ | ≤ 10 × pavimento | PASS (8,9 ×) |
| C2b stazionarietà grosso (tutte le celle) | 3,0·10⁻⁴ | **1,1·10⁻⁵** | | ≤ 1·10⁻⁴ | PASS |
| C3 corrente fuori dal blocco (α c/a0³) | 2,7·10⁻⁴ | **1,24·10⁻⁴** | 7,1·10⁻⁷ | ≤ 10 × pavimento | **FAIL** (175 ×, previsto) |
| corrente sulle celle coperte | 1,66·10⁻³ | 1,48·10⁻³ | | descrittivo | |
| stazionarietà del blocco fine | 5,3·10⁻⁴ | 3,0·10⁻⁴ | | descrittivo | |

- **Stazionarietà: 30 volte meglio**, ora a 9 volte il pavimento.
- **Corrente: non risolta.** Fuori dal blocco resta 1,2·10⁻⁴ ed è **costante** nel tempo, quindi non è un transitorio.
  È un flusso fermo della cucitura: il grosso "vede" le celle coperte (che sono solo la media del campo fine) con una
  fase un po' diversa dalla sua, e la legge come una corrente.
  - In velocità vale j/ρ ≈ 4·10⁻⁴ α c, cioè 4 decimillesimi della velocità dell'elettrone nell'1s.
  - Previsione scritta prima: FAIL. Giusta.

**Convergenza nel blocco** (N = 128, 144 passi, finestra [72,144], partenza veff; ibrido/grosso − 1 in 10⁻⁴):

| blocco | 3 | 5 | 7 | 9 | 11 |
|---|---|---|---|---|---|
| fase 4 (bordo giusto) | 0,18 | 2,05 | **5,58** | **5,30** | **5,56** |
| sviluppo, bordo vecchio (48 passi) | 0,72 | 2,37 | 4,88 | 3,44 | 2,60 |

- **V2 plateau: PASS.** 7, 9 e 11 stanno dentro ±1,5·10⁻⁵ (|E11−E9| = 2,6·10⁻⁵, soglia 3·10⁻⁵; |E9−E7| = 2,8·10⁻⁵, soglia
  1·10⁻⁴). In fase 3 i blocchi 7 e 9 differivano di 1,7·10⁻⁴ e la serie col bordo vecchio **scende** dopo il 7 invece
  di fermarsi.
- **V1 monotona: FAIL.** Dal 7 al 9 l'energia risale di 2,8·10⁻⁵ (soglia 1·10⁻⁵). Il residuo della cucitura è sceso
  da ~2·10⁻⁴ a ~3·10⁻⁵, ma non è zero.
- **V3 "picco": FAIL** (1,9·10⁻⁴). Questo esame però misura la curvatura: una serie che sale e poi si ferma al 7 lo
  sbaglia per costruzione. La previsione "V3 PASS" era sbagliata e lo dico.
- Detto semplice: con il bordo giusto **il blocco converge** (plateau a 5,5·10⁻⁴ sopra il grosso). La bolla al 7 è
  rimasta ma è 6 volte più piccola.

**Energia a N = 256 (C4/R)** (±8 a0, blocco 9, 120 passi, partenza veff, regola nuova): **−13,591198 eV, −0,0535 % da
Dirac×μ**, PASS (soglia 0,10 %). Lettura stabile fra le finestre entro 2·10⁻⁶.
- Il grosso da solo (motore Rust, potenziale come media di cella) sta a −0,1116 %; il blocco fine porta +0,058 %.
- In fase 3 l'ibrido stava a −0,021 %, ma lì aiutava la massa sbagliata del maestro (+0,037 %). Tolta quella, il numero
  è **più onesto e un po' più lontano**: il −0,05 % che resta è il potenziale "media di cella" del grosso del Rust fuori
  dal blocco (vedi §13.3: nel Rust non c'è ancora il campo di Galerkin).

### 13.2 La massa del maestro

Il maestro è il cammino zig-zag fine da cui l'agglomeratore impara la regola grossa. Con cella θ = a0/4096 la sua massa
cinetica è tan θ/θ = 1 + θ²/3 = **1,000373**. La regola la eredita, e si vede da sola nel 2p (+0,034 %).

Non ho corretto la massa a mano. Ho dato all'agglomeratore **un controllo in più**: "il maestro è abbastanza fine?".
- Rifà la catena con maestri sempre più fini, θ_j = θ_A/2^j. Stessa regola zig-zag, stessa cella grossa, stesso learner
  della fase 3 (c2, c4, c6 fissati, parità esatta, kcheck a tempo fisico fisso).
- Si ferma al primo j in cui la massa della regola a0/16 cambia meno di **1·10⁻⁵** (budget pre-registrato).

| maestro | a0/4096 | a0/8192 | a0/16384 | a0/32768 | **a0/65536** |
|---|---|---|---|---|---|
| m − 1 della regola | 3,73·10⁻⁴ | 9,33·10⁻⁵ | 2,33·10⁻⁵ | 5,83·10⁻⁶ | **1,46·10⁻⁶** |
| \|Δm\| | | 2,8·10⁻⁴ | 7,0·10⁻⁵ | 1,7·10⁻⁵ | **4,4·10⁻⁶** |

| esame | valore | soglia | esito |
|---|---|---|---|
| M1 il controllo si ferma | a j = 4 (a0/65536), \|Δm\| = 4,4·10⁻⁶ | < 1·10⁻⁵ entro j ≤ 4 | PASS |
| M2 massa della regola | m − 1 = 1,46·10⁻⁶ (a0/16), 1,47·10⁻⁶ (a0/8) | ≤ 2·10⁻⁵ | PASS |
| M3 2p nel motore Rust (a0/8, ±16 a0, media bordi 0 e π) | **+0,0004 %** (prima +0,0338 %) | ±0,005 % | PASS |
| M4 kcheck su ogni livello | packet ≤ 4·10⁻²¹, energia ≤ 3,4·10⁻⁵ | come fase 3 | PASS |

- Il passo da 1,000373 a 1,0000015 è esattamente θ²/3 che scala come 4⁻ʲ, come previsto. È **validazione**.
- Costa poco: ogni maestro è un conto in k, ~6 s. Il cammino più fine non costa niente a run-time, perché la regola
  grossa ha la stessa forma.
- Nuovo asset: `lab/out/idrogeno_bohr_fase4_regola.json` (stesso formato degli `agglomerato_zz_*.json`; livelli con
  B_fine_cells = 4096 per a0/16 e 8192 per a0/8).

### 13.3 Il grosso di default (emendamento)

Da ora il livello grosso di default è:
- la regola T (quella nuova del maestro a0/65536);
- i **bordi ritorti**: media degli stati con θ = 0 e θ = π su tutti e tre gli assi (toglie i "vicini" della scatola);
- la **coda di Born** sulla cella del nucleo: ΔE_B = −2m(4πZα)²ρ_c J e ρ_B = ρ_c/(1 − 8πmZα I)², con I e J costanti
  del cubo di Brillouin. **0 parametri.**

Motore: `lab/agg2s_core.py` in sola lettura (campo di Galerkin, a0/8, N = 256, L = 32 a0), più il dV del primo ordine
dell'accoppiamento come in AGGLOMERATO_2S. Script: `lab/idrogeno_bohr_fase4_grosso.py`.

| livello | prima (maestro a0/4096) | fase 4 | senza Born | soglia | esito |
|---|---|---|---|---|---|
| 1s1/2 | +0,026 % (+0,028 % T) | **+0,0029 %** | −0,0089 % | ±0,005 % | PASS |
| 2s1/2 | +0,0347 % | **+0,0034 %** | −0,0025 % | ±0,005 % | PASS |
| 2p1/2 | +0,038 % | **+0,0006 %** | | ±0,005 % | PASS |
| 2p3/2 | +0,038 % | **+0,0006 %** | | ±0,005 % | PASS |
| densità sul nucleo, 1s | 0,916 → 0,998 | grezza 0,918 → Born **1,0007** | | entro 0,5 % | PASS (G2) |

- **G1 e G2 PASS.** I quattro livelli n = 1, 2 stanno entro **0,0034 %** da Dirac. Prima stavano a +0,026…+0,038 %,
  e quasi tutto era la massa del maestro.
- Quello che resta sugli s (+0,003 %) è la coda di Born che corregge un po' troppo: senza Born gli s stanno a
  −0,009 % (1s) e −0,0025 % (2s).
- **Nell'ibrido** la coda di Born non si usa: lì il nucleo lo risolve il blocco fine (dichiarato).
- **Debito:** il motore Rust ha i bordi ritorti (`IDRO_TORSIONE`), ma non il campo di Galerkin. Quindi Born nel Rust
  non c'è ancora, e il grosso del Rust (media di cella del potenziale) resta più "basso" sugli s di ~0,1 %. È per questo
  che il blocco fine nell'ibrido guadagna 5,5·10⁻⁴ (in gran parte è questo).

### 13.4 Conto della fase 4

**Esami.**
- **PASS (11):** C1, C2, C2b, C4 (energia N = 256), V2, M1, M2, M3, M4, G1, G2.
- **FAIL (3):** C3 corrente fuori dal blocco (175 × il pavimento, previsto), V1 monotonia (+2,8·10⁻⁵ dal 7 al 9),
  V3 picco (esame di curvatura, previsione sbagliata).
- G3 è descrittivo.

**Manopole.** Nessuna fittata. Impostazioni dichiarate: M = 32 e R della partenza veff (solo lo stato iniziale), budget
del maestro 1·10⁻⁵ e JMAX = 4.

**Che cosa non ha funzionato.**
- **Il Δ imparato sulle celle coperte** (comando `impara`: una correzione della regola grossa con le classi di simmetria
  del cubo, dal passo del blocco fine). Il fit non regge: residuo 0,65–1,26, e l'energia salta di 1,1·10⁻³. Scartato.
- Le varianti di cucitura "fase" (reflux complesso) e "cella" (ampiezza grossa e fase fine) **non** tolgono
  la corrente delle celle coperte (1,5–1,7·10⁻³). La causa vera della perdita era il bug del bordo.
- La corrente ferma al bordo del blocco resta.

**Debiti aperti.**
- **Corrente della cucitura**: ~1·10⁻⁴ fuori dal blocco, ferma. Serve una restrizione che consegni al grosso anche la
  fase "giusta per lui", non solo la media. Candidato: restringere il campo fine **dopo** averlo portato nella base
  della regola grossa (oggi la media di blocco è quella della componente grande).
- **Bolla di 3·10⁻⁵** al blocco 7.
- **Born e campo di Galerkin nel Rust** (vedi §13.3).
- Il dV del primo ordine dell'accoppiamento è ancora quello imparato in AGGLOMERATO_COMPLETO col maestro vecchio. Pesa
  ~10⁻⁵ eV, letto e non rifatto.

**Come rifarlo.**
```
python3 lab/idrogeno_bohr_fase4_regola.py J        # J = 0..4, poi:
python3 lab/idrogeno_bohr_fase4_regola.py scegli   # -> lab/out/idrogeno_bohr_fase4_regola.json
A=lab/out/idrogeno_bohr_fase4_regola.json
IDRO_ASSET=$A IDRO_TORSIONE=000 target/release/idrogeno livello 8192 256 2p     # e 111
for j in 1s_000 1s_111 2s_000 2s_111 2p_000 2p_111; do python3 lab/idrogeno_bohr_fase4_grosso.py $j; done
python3 lab/idrogeno_bohr_fase4_grosso.py esame
export IDRO_ASSET=$A IDRO_ORD=9 IDRO_B=16 IDRO_REFLUX=1
IDRO_DIAG=48 target/release/idrogeno esame_ibrido 4096 128 0 240                          # pavimento
IDRO_DIAG=48 IDRO_PARTENZA=veff:32:3 target/release/idrogeno esame_ibrido 4096 128 7 240  # cucitura
IDRO_PARTENZA=veff:32:2 IDRO_CKPT=<file> IDRO_TEMPO=540 target/release/idrogeno esame_ibrido 4096 128 <blocco> 144
python3 lab/idrogeno_bohr_fase4_esame.py <cartella dei log>; python3 lab/idrogeno_bohr_fase4_figura.py
```
