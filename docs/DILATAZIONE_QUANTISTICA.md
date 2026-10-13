# Dilatazione quantistica del tempo con l'orologio interno della regola

Studio `dilat_` (file: `lab/dilat_*.py`, `lab/out/dilat_*`). Regola: quella congelata di `lab/zz_core.py`
(`Walk('dirac1')`, 1+1 D, libera). Unità ħ = c = m = 1.

Commit: ricerca + pre-registrazione `909edb6`, corse + post-hoc + figura `2a2c4e9`, documento: vedi `git log`.
Figura: `lab/out/dilat_figura.png` (copia in `/mnt/project-files/problemi-aperti/dilat_figura.png`).

## 0. In breve (per chi ha fretta)

* **La regola ha un orologio suo, ed è giusto.** Contare le rotazioni (o, che è lo stesso, confrontare due copie
  con massa m ± δ/2) dà il tempo proprio: per un pacchetto a impulso p, τ/t → ⟨m/E⟩ = 1/γ con errore che va
  **esattamente** come a² (rapporti 4,00–4,17) e con il coefficiente che avevamo previsto a carta e penna
  (−0,2502 contro −0,2501 a p = mc). VALIDAZIONE.
* **Formula esatta della regola (nuova, piccola):** il ritmo dell'orologio è **dτ/dt = √(1 − v²) · cos(ka)**,
  v = velocità di gruppo della regola. Il fattore cos(ka) è il "difetto" del reticolo: a impulso fissato sparisce
  come 1 − (pa)²/2, ma vicino al **doppione** (ka = π) vale −1 e **l'orologio gira all'indietro**.
* **La dilatazione quantistica: la regola dà la stessa di Smith–Ahmadi**, non una diversa. Nel regime lento
  (p ≤ 0,1 mc) la regola coincide con la loro eq. 25 entro 1,8 % (la differenza è tutta l'approssimazione non
  relativistica della formula; con il continuo esatto la regola sta a 4·10⁻⁴). È VALIDAZIONE, e c'è un motivo
  semplice: per una particella libera il ritmo dipende solo da |ψ(p)|², e la regola conserva k come il continuo.
* **Due cose che la formula di Smith–Ahmadi non dice e la regola (come il continuo esatto) sì:**
  (1) per p ~ mc la correzione è **più piccola e può cambiare segno**: con pesi uguali Δ_Q ∝ −f″(p̄),
  f = m/E, quindi **si annulla a p̄ = mc/√2** (v = c/√3). Misurato: zero a p̄ = 0,7103 (regola, a = 1/32),
  0,7100 (continuo), flesso 0,7071; S–A dice un valore costante. (2) Con i pacchetti spostati di d nello spazio
  la correzione oscilla come una frangia cos(p̄d + φ) e muore come e^{−Δ²d²/4}; ma **una volta preparata non
  cambia più nel tempo**, anche quando i pacchetti si allontanano (ritmo costante a 3·10⁻¹⁰ nella regola vera).
  Sono conseguenze del continuo relativistico, non della regola: **non le chiamiamo scoperte**.
* **Dove la regola dice qualcosa in più, perché l'orologio è la massa di un elettrone di Dirac (non un orologio
  scalare):** (a) se si prepara il qubit-orologio "di colpo" (stesso spinore per le due masse), nasce una
  **Zitterbewegung nella fase dell'orologio** grande ≈ 4–5 (δp/4E²)², dello stesso ordine dell'effetto
  "non classico" di Paige–Plato–Kim finché δ·t ≲ 2; a tempi lunghi torna esattamente PPK (|D − PPK| = 3·10⁻⁶
  contro una separazione PPK–massa media di 8,5·10⁻⁵). (b) Uno stato di sola chiralità fa tremare l'orologio e il
  ramo a energia negativa conta **all'indietro**: ritmo medio ⟨v·m/E⟩ = 0,4292 (regola 0,4290).
* **Esami pre-registrati: 10 PASS, 6 FAIL.** I FAIL sono tutti soglie di precisione troppo strette a a = 1/32 per
  impulsi alti (errore O(a²) perfetto ma coefficiente grande), un pavimento numerico (3·10⁻¹² contro 10⁻¹²), e
  l'esame E3 dove avevamo previsto con il continuo **scalare** e la regola segue invece il continuo di **Dirac**
  (vedi §4.4). Nessun FAIL è un errore della regola rispetto alla fisica; tutti spiegati in §5 con controlli
  post-hoc dichiarati.
* **Esperimento:** in uno ione che simula Dirac (Gerritsma 2010) l'orologio di massa è ⟨σ_z⟩, che si legge con
  la fluorescenza. Con lo stato fondamentale (Δp = 0,31 mc) e p̄ = 0,5 e 0,94 mc la correzione quantistica è
  **0,060** (6 %) sul ritmo: ~2·10⁴ ripetizioni per 5σ. Visibilità dell'orologio a 244 µs: 0,58 (sovrapposizione)
  contro 0,28 (miscuglio).

## 1. Il problema e cosa è noto (ricerca, 13/10/2026)

**L'idea.** Un orologio che si muove va più piano: τ = t/γ. Se però l'orologio è una particella quantistica,
può muoversi in una **sovrapposizione** di due velocità. Che tempo segna? La risposta "classica" sarebbe la
media dei due tempi (come se avessimo un miscuglio: metà delle volte va piano, metà veloce). La letteratura
dice che c'è un pezzo in più, che dipende dalla **coerenza** tra i due pacchetti: la "dilatazione quantistica".

Fonti lette (online, 13/10/2026):

| lavoro | cosa dice | numeri |
|---|---|---|
| Zych, Costa, Pikovski, Brukner, *Nat. Commun.* 2, 505 (2011) | orologio (livello interno) in un interferometro: la differenza di tempo proprio tra i due bracci si vede come **perdita di visibilità** | proposta, mai fatta col tempo proprio vero |
| Smith & Ahmadi, *Nat. Commun.* 11, 5360 (2020) | orologio = gradi interni di una particella, H = √(p²c² + M²c⁴), M = m + H_clock/c². Pacchetto stretto → dilatazione classica. Sovrapposizione cos θ\|p̄⟩ + e^{iφ} sin θ\|p̄′⟩ → ⟨T⟩ = (γ_C⁻¹ + γ_Q⁻¹) τ | eq. 24: γ_C⁻¹ = 1 − (p̄² cos²θ + p̄′² sin²θ)/2m²c² − Δ²/4m²c² (il miscuglio); eq. 25: γ_Q⁻¹ = cos φ sin 2θ [(p̄′−p̄)² − 2(p̄′²−p̄²) cos 2θ] / {8m²c² [cos φ sin 2θ + e^{(p̄′−p̄)²/4Δ²}]}. Per ⁸⁷Rb a 5 e 15 m/s: γ_Q⁻¹ ~ 10⁻¹⁵ |
| Paige, Plato, Kim, *PRL* 124, 160602 (2020) | un orologio in sovrapposizione di livelli interni ha una sovrapposizione di **masse**: la fase relativa tra i livelli n, m si dilata di 1 − p²/(2MₙMₘc²) (spinta in impulso), non con un γ classico; con spinta in velocità torna il fattore classico | correzione ~10⁻¹⁰ dello spostamento già visto negli ioni |
| Khandelwal, Lock, Woods, *Quantum* 4, 309 (2020) | orologi delocalizzati: correzioni quantistiche universali alla dilatazione | (non letto in dettaglio) |
| Grochowski, Smith, Dragan, Dębski, *PRR* 3, 023053 (2021) | dilatazione quantistica negli spettri atomici (tasso di emissione di un atomo in sovrapposizione di impulsi) | |
| Paczos, Dębski, Grochowski, Smith, Dragan, *Quantum* 8, 1338 (2024) | versione gravitazionale (sovrapposizione di quote) | spostamento frazionario alla portata degli orologi atomici |
| Dębski et al., *CQG* 41, 135014 (2024) | "universalità" della dilatazione quantistica | (solo titolo) |
| arXiv:2509.09573 (PRL 2026) "Quantum signatures of proper time in optical ion clocks" | **proposta** con ²⁷Al⁺: lo spostamento in frequenza da solo NON distingue quantistico da semiclassico; il testimone vero è la **perdita di visibilità** da entanglement moto-orologio (stati compressi) | V ≈ 0,93 a t ≈ 1 s con compressione r = 2,26 |
| Gerritsma et al., *Nature* 463, 68 (2010) | simulazione di Dirac 1+1 D con uno ione ⁴⁰Ca⁺: H = c p σ_x + mc² σ_z, c = 2ηΩ̃Δ, mc² = ħΩ | η = 0,06, ν = 2π·1,36 MHz, Ω̃ = 2π·68 kHz, 0 < Ω ≤ 2π·13 kHz |
| Lan et al., *Science* 339, 554 (2013) (+ critiche di Wolf et al.) | "orologio di Compton": la massa come orologio a frequenza mc²/ħ | controverso: la fase assoluta non si misura |

**Cosa è aperto.** La dilatazione quantistica non è mai stata misurata (solo proposte, fino al 2026). In
più, quasi tutta la teoria usa una particella scalare con H = √(p² + M²), al primo ordine in H_clock/mc² e
spesso al secondo ordine in p/mc (la formula di Smith–Ahmadi è non relativistica).

**Un punto chiave (lo diciamo subito, non è nostro).** Per una particella libera l'operatore "ritmo
dell'orologio" mc²/H(p) è diagonale in impulso. Quindi il tempo letto dipende **solo da |ψ(p)|²**. La
"correzione quantistica" di Smith–Ahmadi è l'interferenza dei due pacchetti dentro |ψ(p)|²: esiste solo se i
due pacchetti si **sovrappongono in impulso** (fattore e^{−(Δp)²/4Δ²}) e, come aggiungiamo noi in §4, anche
nello spazio (fattore e^{−Δ²d²/4} e una frangia cos(p̄d + φ)). Se i pacchetti non si toccano in impulso, il
miscuglio e la sovrapposizione danno lo stesso tempo.

## 2. Perché la regola c'entra, e la previsione (prima dei dati)

Nella regola zig-zag la massa è l'angolo di rotazione θ = m a. Il giro 1 (`docs/TEMPO_TUNNEL.md`) ha mostrato
che contare le rotazioni (pesate con le ampiezze) è un orologio attaccato alla massa e legge il tempo proprio.
Qui usiamo due versioni dello stesso orologio:

* **orologio a qubit** (sonda dichiarata, come l'orologio di Larmor del giro 1): lo stesso stato evolve con
  massa m − δ/2 e m + δ/2 (un livello interno che pesa un po' di più). Si legge O(t) = ⟨ψ₋(t)|ψ₊(t)⟩:
  tempo proprio τ = −arg O/δ, visibilità V = |O|. È esattamente l'orologio di Zych / Smith–Ahmadi / Paige–Plato–Kim.
* **contatore di rotazioni** (dentro la regola, nessun grado di libertà in più): ogni flip ha un peso z;
  f_w = ∂ ln⟨ψ(n)|ψ_z(n)⟩/∂ ln z. Identità esatta: Re f_w = n sin²θ, Im f_w = −(sin 2θ/2a) τ. Contare i flip =
  orologio di massa con δ → 0.

**Il ritmo dell'orologio della regola** (nostro calcolo): dalla dispersione cos ωa = cos θ cos ka,

  r(k) = ∂ω/∂m = sin θ cos ka / sin ωa  →  m/E (a → 0),   r(k) ≈ (m/E)[1 − (k²a²/3)(1 − m²/2E²)].

Quindi **previsione**: la regola dà la stessa dilatazione (classica e quantistica) del continuo relativistico
esatto, cioè ⟨m/E⟩ su |ψ(k)|², con un errore O(a²) dal coefficiente sopra. Nel regime non relativistico deve
coincidere con Smith–Ahmadi eq. 25; per p ~ mc no (la formula di S–A è al secondo ordine in p). Con δ finito
deve dare Paige–Plato–Kim (1 − p²/2M₊M₋), non la massa media. Cose specifiche della regola che prevediamo:
(i) r² + v² = 1 − sin²ka sin²θ / sin²ωa < 1: **sul reticolo il tempo proprio non è √(1 − v²)** (errore O(a²));
(ii) attorno al **doppione** k = π/a l'orologio **gira all'indietro** (r → −m/E_q); (iii) uno stato con energie
di due segni (es. chiralità pura) fa la **Zitterbewegung anche nell'orologio**, e la parte a energia negativa
conta il tempo all'indietro.

Pre-registrazione: `lab/out/dilat_prereg.json` (generata da `lab/dilat_prereg.py`), esami E0–E6 con soglie.

## 3. Cosa abbiamo fatto (in parole semplici)

**La regola** (`zz_core.Walk('dirac1')`, nessuna modifica): palline R che corrono a destra a c, palline L a
sinistra, e a ogni tick in ogni cella si mescolano un po' (angolo θ = m a). Niente potenziali.

**Due motori, stessa regola.** (1) *Nel tempo*: la regola vera tick per tick, con tre copie (massa m − δ/2,
m + δ/2, e il contatore di flip con peso z). (2) *Spettrale*: la regola libera non mescola impulsi diversi, quindi
per ogni k un tick è una matrice 2×2 W(k) ∈ SU(2) e W(k)ⁿ si scrive in forma chiusa. I due motori coincidono a
1,7·10⁻¹³ (esame E0); lo spettrale permette di fare migliaia di configurazioni in secondi.

**Gli stati.** Pacchetti gaussiani di sola energia positiva (proiettore di banda della regola), alla Smith–Ahmadi:
g(p) = e^{−(p−p̄)²/2Δ²}. La sovrapposizione è cos θ g₁ + e^{iφ} sin θ g₂ (θ, φ come S–A). Il "miscuglio" è la
media cos²θ·(corsa 1) + sin²θ·(corsa 2).

**Il passo reticolare.** a = 1/4, 1/8, 1/16, 1/32 ħ/mc (fino a 200 celle per lunghezza d'onda a p = mc). Scatola
1024 ħ/mc (32 768 celle a a = 1/32).

## 4. Risultati

### 4.1 Validazione: dilatazione classica (E1, R1)

| p̄ [mc] | ⟨m/E⟩ (continuo) | regola a = 1/32 | errore rel. | rapporti errore a→a/2 | coeff. a² misurato / previsto |
|---|---|---|---|---|---|
| 0 | 0,999377 | 0,999376 | −2,6·10⁻⁷ | 4,03 4,03 3,11 | (≈ 0) |
| 0,5 | 0,894248 | 0,894204 | −4,9·10⁻⁵ | 4,04 4,01 4,00 | −0,0503 / −0,0503 |
| 1 | 0,707217 | 0,707044 | −2,4·10⁻⁴ | 4,06 4,01 4,00 | −0,2502 / −0,2501 |
| 2 | 0,447292 | 0,446768 | −1,2·10⁻³ | 4,10 4,03 4,01 | −1,2005 / −1,1998 |
| 3 | 0,316261 | 0,315380 | −2,8·10⁻³ | 4,17 4,04 4,01 | −2,852 / −2,850 |

Estrapolando a a → 0 (Richardson su 1/16 e 1/32, post-hoc) l'errore scende a ≤ 1·10⁻⁵ per tutti.

**La formula esatta.** Dalla dispersione della regola cos ωa = cos θ cos ka seguono tre righe di algebra:
√(1 − v²) = sin θ / sin ωa, e quindi **r = ∂ω/∂m = √(1 − v²) · cos(ka)** (verificata a 10⁻¹⁴, `dilat_r_R1.json`).
Cioè: l'orologio della regola segna "il tempo proprio della velocità con cui si muove davvero il pacchetto",
moltiplicato per cos(ka). A p = mc e a = 1/8 il fattore è 0,9922: **l'orologio della regola va più piano di
√(1 − v²) dello 0,8 %**; a a = 1/32 dello 0,05 %. Sparisce come (pa)²/2: è un difetto di Lorentz del reticolo,
non fisica.

### 4.2 La dilatazione quantistica (E2)

**Regime lento (E2a, contro Smith–Ahmadi eq. 25).** Δ = 0,02 mc, p̄₁ = 0,02 mc, p̄₂ = p̄₁ + sΔ. Esempi (θ = π/8):

| s | φ | regola (a = 1/32) | continuo esatto | S–A eq. 25 |
|---|---|---|---|---|
| 0,5 | 0 | −3,0250·10⁻⁵ | −3,0242·10⁻⁵ | −3,0290·10⁻⁵ |
| 0,5 | π | 1,4995·10⁻⁴ | 1,4991·10⁻⁴ | 1,5015·10⁻⁴ |
| 1 | π | 1,9846·10⁻⁴ | 1,9840·10⁻⁴ | 1,9872·10⁻⁴ |
| 2 | π | 1,2842·10⁻⁴ | 1,2839·10⁻⁴ | 1,2857·10⁻⁴ |
| 4 | π | 1,1772·10⁻⁵ | 1,1768·10⁻⁵ | 1,1770·10⁻⁵ |

Su 45 configurazioni: regola contro S–A al massimo 1,8 % (soglia 3 %: PASS), regola contro continuo esatto
4,3·10⁻⁴ (soglia 10⁻³: PASS). Con φ = π/2 la correzione è zero (S–A: cos φ); la regola dà 3·10⁻¹² con la fase
(FAIL contro la soglia 10⁻¹², pavimento numerico) e 2·10⁻¹⁶ con la media spettrale esatta (post-hoc).

Quanto è grande? Per questi parametri (velocità ~0,02–0,1 c) la correzione è 10⁻⁵–2·10⁻⁴, cioè dello stesso
ordine della dilatazione classica stessa (la dilatazione classica qui è 1 − 0,9996 = 4·10⁻⁴): **la parte
quantistica non è una piccola correzione della classica, è dello stesso ordine** quando i pacchetti si
sovrappongono in impulso (s ≲ 2). Muore come e^{−s²/4}.

**Regime relativistico (E2b).** Δ = 0,3 mc, p̄₁ = 0,5 mc. Qui la formula di S–A (al secondo ordine in p/mc)
sbaglia di un fattore 1,6–100 e anche il segno; la regola segue il continuo esatto:

| s | θ | φ | regola (1/32) | continuo esatto | S–A (NR) |
|---|---|---|---|---|---|
| 1 | π/8 | π | 4,3902·10⁻² | 4,3877·10⁻² | 7,07·10⁻² |
| 1 | π/4 | π | −5,2695·10⁻³ | −5,2576·10⁻³ | −3,96·10⁻² |
| 2 | π/4 | 0 | −2,868·10⁻⁵ | −3,166·10⁻⁵ | +1,21·10⁻² |
| 3 | π/4 | π | 9,845·10⁻⁴ | 9,867·10⁻⁴ | −1,19·10⁻² |

La soglia di precisione (10⁻³ relativo a a = 1/32) fallisce dove la correzione quasi si annulla (s = 2, θ = π/4:
errore 9 %); l'ordine a² è perfetto (16 rapporti tutti 4,00–4,01) e con Richardson l'errore massimo è 2,7·10⁻⁴.

**Perché cambia segno (esplorativo, `lab/dilat_segno.py`, non pre-registrato).** Con pesi uguali la correzione è
Δ_Q ≈ −[O/(1+O)] f″(p̄) (Δp)²/8 con f(p) = m/E e O = e^{−(Δp)²/4Δ²} (per f = 1 − p²/2m² è esattamente S–A).
La dilatazione quantistica **misura la curvatura della curva 1/γ(p)**. Questa curva ha un flesso a p = mc/√2:
lì la correzione si annulla e sopra cambia segno. Scansione (Δ = Δp = 0,05 mc): zero a p̄ = 0,7100 (continuo),
0,7103 (regola, a = 1/32), flesso 0,7071; la formula reggeva entro 1 % su tutta la scansione 0–2 mc, mentre S–A
dà un valore costante (1,37·10⁻⁴). Non abbiamo trovato questa affermazione nei lavori letti (S–A, PPK, Paczos
2024), ma è una conseguenza di due righe del continuo: **la classifichiamo come osservazione del continuo,
non come scoperta della regola**.

**Pacchetti separati nello spazio (E2c, PASS).** Stessa algebra di S–A con g₂ spostato di d: il termine incrociato
prende e^{−Δ²d²/4}·e^{−ip̄d} e un centro complesso p̄ − iΔ²d/2. Regola contro questa formula: entro 1,3 % (il
resto è la NR). La correzione cambia segno con d (frangia cos(p̄d + φ)):
d·Δ = 0, 0,5, 1, 2, 3 → Δ_Q = 5,34, 8,20, 8,44, −12,82, 1,07 (·10⁻⁵, φ = 0).

**Persistenza (E2d, PASS, regola vera tick per tick).** p̄ = 0,5 e 0,8, Δ = 0,3, a = 1/16, 3200 tick: i centri
si allontanano di 35 ħ/mc (15 volte la larghezza iniziale), le frange della densità cambiano (figura, pannello a),
ma il ritmo dell'orologio della sovrapposizione è costante a **2,9·10⁻¹⁰** e Δ_Q = 6,60·10⁻⁴ in ogni finestra.
Quello che conta non è se i pacchetti si toccano *adesso*, ma l'interferenza in impulso fissata alla preparazione.

### 4.3 Visibilità dell'orologio (E4, PASS)

Orologio con δ = 0,05 mc², p̄ = 0,5 e 1,0, Δ = 0,3, θ = π/4. La visibilità |⟨ψ₋|ψ₊⟩| (Zych) scende perché le
componenti di impulso segnano tempi diversi. Regola (a = 1/32) contro continuo: 5,8·10⁻⁴ al massimo (errori
3,2·10⁻², 7,9·10⁻³, 2,0·10⁻³, 4,9·10⁻⁴: O(a²)).

| t [ħ/mc²] | 0 | 50 | 100 | 200 | 400 | 800 |
|---|---|---|---|---|---|---|
| V sovrapposizione | 1 | 0,967 | 0,874 | 0,563 | 0,010 | 0,027 |
| V miscuglio | 1 | 0,960 | 0,846 | 0,476 | 0,144 | 0,041 |

La differenza tra sovrapposizione e miscuglio è grande (0,56 contro 0,48 a t = 200; 0,01 contro 0,14 a t = 400):
è lo stesso effetto (|ψ(p)|² diverso), visto nella visibilità invece che nel ritmo. La proposta per ioni del 2026
(arXiv:2509.09573) dice proprio che la visibilità è il testimone giusto.

### 4.4 Orologio a qubit con δ finito: Paige–Plato–Kim e la Zitterbewegung dell'orologio (E3)

Con δ finito la fase relativa delle due masse cresce, al primo ordine in p², con 1 − p²/(2M₊M₋) (PPK), non con la
massa media 1 − p²/2m². La pre-registrazione confrontava la regola a t = 2 con il continuo **scalare**: FAIL
(errore 9 % a δ = 0,5). La causa, trovata dopo (post-hoc, dichiarato): l'elettrone della regola è di **Dirac**.
Lo spinore di energia positiva dipende dalla massa; uno stato preparato per la massa m non è di energia positiva
per m ± δ/2, quindi l'orologio "acceso di colpo" fa una **Zitterbewegung nella sua fase**:

* la regola segue il continuo di Dirac con la stessa preparazione: 1,238750·10⁻³ contro 1,238414·10⁻³ (a = 1/32,
  errore O(a²));
* la fase "tremante" ha ampiezza ≈ 4–5 ε², ε = δp/4E² (misurato su 6 casi: rapporto 4,1–5,0);
* il ritmo a tempi lunghi (fit su t ∈ [50, 250], Richardson) torna lo scalare a 3·10⁻⁵ relativo, e il test PPK
  passa: |D − PPK| = 3,3·10⁻⁶ contro |PPK − massa media| = 8,5·10⁻⁵ (p̄ = 0,05, δ = 0,5).

Messaggio: il termine "non classico" di PPK accumula una fase ≈ δ³p²t/8, il tremolio di Dirac è limitato a
≈ δ²p²/4: il primo vince solo per **δ·t ≳ 2**, cioè dopo che l'orologio ha fatto qualche giro. Per orologi atomici
veri (δ/mc² ~ 10⁻¹⁰) è irrilevante; per un simulatore (δ ~ mc²) no.

### 4.5 Cose solo della regola: doppione e chiralità (E5, E6)

* **Doppione (E5, PASS).** Un pacchetto attorno a k = π/a − q ha ritmo −⟨m/E_q⟩: −0,894204 contro −0,894248
  (q = 0,5), −0,707044 contro −0,707217 (q = 1). L'orologio gira **all'indietro** (il doppione è un "elettrone"
  con massa di segno opposto per l'orologio). Una sovrapposizione elettrone + doppione ha ritmo 2·10⁻⁹: l'orologio
  si ferma. In 1-D è un artefatto noto della regola (stare lontani da k = π); qui ha una firma netta.
* **Chiralità pura (E6).** Stato di sole palline R (σ_x = 1, p̄ = 1): contiene energie + e −. L'orologio trema alla
  frequenza 2E (Zitterbewegung, figura pannello e) e il ritmo medio è ⟨(1+v)/2 − (1−v)/2⟩·m/E = ⟨p m/E²⟩ = 0,4292:
  il pezzo a energia negativa **conta il tempo all'indietro**. Regola: 0,4290 (PASS, 6·10⁻⁴). La traccia τ(t)
  fallisce la soglia 2·10⁻³ a a = 1/32 (1,0·10⁻²), errore O(a²) netto (0,66, 0,16, 0,041, 0,010); con Richardson
  9·10⁻⁶.

### 4.6 Numeri per un esperimento con uno ione (domanda 4, `lab/dilat_ione.py`)

Gerritsma et al. 2010: ⁴⁰Ca⁺, H = c p σ_x + mc² σ_z, c = 2ηΩ̃x₀, mc² = ħΩ, η = 0,06, Ω̃ = 2π·68 kHz,
Ω ≤ 2π·13 kHz. Con Ω = 2π·13 kHz: mc = 1,59 ħ/x₀, unità di tempo ħ/mc² = 12,2 µs, e lo stato fondamentale del
moto ha σ_p = 0,31 mc (Δ = 0,44 mc alla S–A): **un pacchetto di larghezza paragonabile a mc**, il regime dove la
correzione quantistica è grande. Lì l'orologio di massa è σ_z (operatore di massa): per uno stato di energia
positiva il ritmo è ⟨σ_z⟩ (Hellmann–Feynman), che si misura con la fluorescenza.

| p̄₁, p̄₂ [mc] | θ, φ | ⟨σ_z⟩ sovrapp. | ⟨σ_z⟩ miscuglio | Δ_Q (continuo / regola a = 1/16) | ripetizioni per 5σ |
|---|---|---|---|---|---|
| 0,5; 0,94 | π/8, π | 0,9171 | 0,8571 | **+0,0599** / 0,0601 | ~2·10⁴ |
| 0,5; 0,94 | π/8, 0 | 0,8398 | 0,8571 | −0,0174 / −0,0174 | ~2,5·10⁵ |
| 0; 0,44 | π/4, π | 0,8733 | 0,9264 | −0,0532 / −0,0533 | ~2,7·10⁴ |
| 0,5; 1,39 | π/8, π | 0,8734 | 0,8368 | +0,0366 / 0,0367 | ~5,6·10⁴ |

(rumore di proiezione ≤ 1/√N per ⟨σ_z⟩, tre misure). Serve preparare la sovrapposizione di due impulsi **di sola
energia positiva**: Gerritsma ha già preparato uno spinore di una sola energia al 99 %; la sovrapposizione di due
spinte è il passo nuovo (spinta dipendente dallo stato + disaccoppiamento). Con un vero orologio a qubit (un livello
in più che "pesa" Ω ± δ/2, δ = 0,5 mc²) la visibilità a t = 20 ħ/mc² = 244 µs sarebbe V² = 0,58 (sovrapposizione)
contro 0,28 (miscuglio) per la prima riga. **Attenzione**: questo è un simulatore. Misura la matematica della
dilatazione quantistica (|ψ(p)|² e coerenza), non il tempo proprio di un oggetto che si muove davvero a 0,7 c.
Non abbiamo messo riscaldamento, decoerenza e errori di preparazione (DEBITO).

## 5. Esami contro la pre-registrazione

| esame | cosa | misurato | soglia | verdetto |
|---|---|---|---|---|
| E0 | regola nel tempo = spettrale; contatore di flip | 1,7·10⁻¹³; Re 3·10⁻¹²; Im 4,9·10⁻⁷ | 10⁻¹⁰; 10⁻⁹; 10⁻⁶ | **PASS** |
| E1a | τ/t = ⟨m/E⟩ a a = 1/32 | 2,6·10⁻⁷ … 2,8·10⁻³ | 10⁻³ | **FAIL** (p̄ = 2, 3) |
| E1b | ordine a² | rapporti 4,00–4,17 | 3,5–4,5 | **PASS** |
| E1c | coefficiente a² previsto | entro 0,1 % | 5 % | **PASS** |
| E2a | contro S–A eq. 25 / continuo | 1,8 % / 4,3·10⁻⁴ | 3 % / 10⁻³ | **PASS** |
| E2a′ | zero per φ = π/2 | 3,1·10⁻¹² | 10⁻¹² | **FAIL** (pavimento numerico) |
| E2b | relativistico: precisione | 9,4·10⁻² (dove Δ_Q ≈ 0) | 10⁻³ | **FAIL** |
| E2b′ | relativistico: ordine | 4,00–4,01 | 3–5 | **PASS** |
| E2c | pacchetti spostati | 1,3 % | 3 % | **PASS** |
| E2d | ritmo costante nel tempo | 2,9·10⁻¹⁰ | 10⁻⁶ | **PASS** |
| E3a | δ finito, t = 2, contro continuo scalare | 9 % | 1 % | **FAIL** (la regola segue Dirac, §4.4) |
| E3b | PPK contro massa media (a t = 2) | 1,2·10⁻⁴ contro 8,5·10⁻⁵ | < 0,2× | **FAIL** (idem; post-hoc a tempi lunghi: 3,3·10⁻⁶, passa) |
| E4 | visibilità | 5,8·10⁻⁴ | 10⁻³ | **PASS** |
| E5 | doppione all'indietro | 2,4·10⁻⁴ | 2·10⁻³ | **PASS** |
| E6a | traccia τ(t) chiralità pura | 1,0·10⁻² | 2·10⁻³ | **FAIL** (O(a²), Richardson 9·10⁻⁶) |
| E6b | ritmo medio ⟨pm/E²⟩ | 6·10⁻⁴ | 2 % | **PASS** |

Totale: 10 PASS, 6 FAIL (contando E2a′, E2b′, E6a/b separati). Post-hoc dichiarati in `lab/dilat_posthoc.py` →
`lab/out/dilat_emendamenti.json` (A1 Richardson, A2 media spettrale, A3 Dirac contro scalare); non cambiano i
verdetti. Lezione: le soglie di precisione a un solo a vanno scalate con il coefficiente a² previsto (che avevamo
calcolato! ma non usato per le soglie), e per orologi di Dirac con δ finito si prevede con Dirac, non con lo scalare.

## 6. Validazione, scoperta, limite (onesto)

**Il confronto è legittimo?** Sì, al primo ordine: nel modello di Smith–Ahmadi l'orologio entra come
M = m + H_clock/c², e al primo ordine il suo ritmo è ∂H/∂m = m/E. Il nostro orologio a qubit è esattamente questo
(due livelli con masse m ± δ/2), e il contatore di rotazioni ne è il limite δ → 0 (identità esatta, E0). Quindi
"l'elettrone della regola è un orologio" nel senso preciso di Hellmann–Feynman: il suo ritmo è ⟨β⟩ = ⟨ψ̄ψ⟩.
**Dove no:** (i) la regola non ha un orologio *separato*: per leggerlo serve la sonda a due masse (debito) o un
valore debole (il contatore di flip, non una misura proiettiva); la fase assoluta dei flip (frequenza di Compton,
Lan et al.) non è osservabile, solo le differenze (la critica di Wolf et al. vale anche qui); (ii) l'orologio è
legato allo spin: con δ finito o stati a due energie compaiono tremolio e conteggio all'indietro, assenti nel
modello scalare; (iii) il doppione ha l'orologio rovesciato; (iv) siamo in 1+1 D e liberi: niente gravità, niente
potenziale.

* **VALIDAZIONE:** dilatazione classica (convergenza a², coefficiente previsto); dilatazione quantistica = S–A nel
  regime lento; PPK a tempi lunghi; visibilità = continuo; Zitterbewegung dell'orologio = Dirac.
* **OSSERVAZIONI DEL CONTINUO (non della regola, non trovate nei lavori letti, piccole):** Δ_Q ∝ −f″(p̄) con zero a
  p̄ = mc/√2 (la formula NR sbaglia segno sopra); frangia spaziale cos(p̄d + φ)e^{−Δ²d²/4}; persistenza dopo la
  separazione; per orologi di Dirac il tremolio ~(δp/2E²)² compete con PPK finché δt ≲ 2.
* **DELLA REGOLA (piccola, esatta):** dτ/dt = √(1 − v²)·cos(ka); orologio all'indietro al doppione.
* **LIMITE:** nessuna previsione diversa dalla letteratura per la particella libera, e non poteva esserci: per una
  particella libera la dilatazione quantistica è una proprietà di |ψ(p)|², e la regola la conserva come il
  continuo. L'unica differenza è il reticolo (cos ka), che va a zero. **Non abbiamo una scoperta sulla dilatazione
  quantistica.** Il posto dove la regola potrebbe dire qualcosa di nuovo è un orologio in un campo (nucleo,
  potenziale, gravità simulata), dove ∂E/∂m non è più diagonale in impulso — non fatto qui.

## 7. Debiti e parametri

* Parametri liberi: **0** (m, c, ħ; tutto il resto è configurazione dell'esame).
* Debiti: (1) la sonda a due masse (qubit orologio con splitting δ) è messa a mano, come l'orologio di Larmor del
  giro 1; (2) stati iniziali gaussiani di sola energia positiva preparati con il proiettore di banda; (3) il
  miscuglio è una media di corse separate, non una dinamica; (4) per lo ione: niente decoerenza/riscaldamento,
  preparazione ideale.

## 8. Come rifare le corse

```
cd lab
export OMP_NUM_THREADS=1
python3 dilat_prereg.py          # previsioni + soglie -> out/dilat_prereg.json (gia' committata)
python3 dilat_run.py all         # E0..E6, R1 -> out/dilat_r_*.json  (~20 s in tutto)
python3 dilat_posthoc.py         # A1-A3 -> out/dilat_emendamenti.json (~90 s)
python3 dilat_segno.py           # esplorativo: zero a p = mc/sqrt2 -> out/dilat_r_segno.json
python3 dilat_ione.py            # numeri per lo ione -> out/dilat_ione.json
python3 dilat_figure.py          # out/dilat_figura.png
```

Fonti online: [Smith & Ahmadi 2020](https://www.nature.com/articles/s41467-020-18264-4),
[Paige, Plato, Kim 2020](https://spiral.imperial.ac.uk/handle/10044/1/79739),
[Paczos et al. 2024](https://quantum-journal.org/papers/q-2024-05-07-1338/),
[arXiv:2509.09573](https://arxiv.org/html/2509.09573v1), [Gerritsma et al. 2010](https://arxiv.org/abs/0909.0674).
