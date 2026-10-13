# Dilatazione quantistica del tempo con l'orologio interno della regola

Studio `dilat_` (file: `lab/dilat_*.py`, `lab/out/dilat_*`). Regola: quella congelata di `lab/zz_core.py`
(`Walk('dirac1')`, 1+1 D, libera). Unità ħ = c = m = 1.

*(documento in costruzione: fase 1 ricerca, fase 2 previsione; i risultati arrivano dopo le corse)*

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
