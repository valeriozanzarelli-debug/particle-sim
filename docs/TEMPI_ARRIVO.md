# Tempi di arrivo delle particelle di spin 1/2 sulla regola zig-zag

Studio `arrivo_` (file: `lab/arrivo_*.py`, `lab/out/arrivo_*`). Regola: quella congelata di `lab/zz_core.py`
(1+1 D `dirac1` e 2+1 D `dirac` con assi y, z): ogni pallina corre di una cella a c, poi ruota (θ = m c a/ħ).
Unità della regola ħ = c = mₑ = 1 (lunghezza λ_C = 3,86·10⁻¹³ m); unità di Das–Dürr ("DD") ħ = m = d = 1
(d = lunghezza della scatola iniziale, tempo m d²/ħ).

> Stato: **fase 1 (ricerca)**. Le sezioni sui risultati vengono aggiunte dopo la pre-registrazione e le corse.

---

## 1. Il problema, in parole semplici

Lanci un elettrone in un tubo e metti un rivelatore a distanza L. **Quando** fa clic? La meccanica quantistica
standard non ha una risposta unica: il tempo non è un'osservabile come la posizione (teorema di Pauli: non c'è
un operatore "tempo di arrivo" autoaggiunto e covariante). Ci sono più ricette, e danno numeri diversi:

* **flusso** (corrente di probabilità nel punto del rivelatore): semplice, ma può diventare negativa
  (*backflow*) e allora non è una probabilità;
* **Bohm**: le particelle hanno traiettorie; il tempo di arrivo è quando la traiettoria tocca il rivelatore la
  prima volta;
* **bordo assorbente** (Tumulka): il rivelatore è un bordo con una condizione al contorno che lascia solo uscire;
  la probabilità del clic è il flusso uscente, sempre ≥ 0;
* **potenziali immaginari**, misure ripetute, modelli di rivelatore espliciti…

## 2. Cosa è noto (letteratura, verificata online il 13/10/2026 dove indicato)

**Das & Dürr, *Sci. Rep.* 9, 2242 (2019), arXiv:1802.07141** (letto il PDF tramite WebFetch). Elettrone in una
guida d'onda semi-infinita (asse z) con confinamento trasverso armonico ω, chiuso da un muro a z = 0, inizialmente
nello stato fondamentale di una scatola 0 < z < d: ψ₀ = √(2ω/π) sin(πz/d) e^{−ωρ²/2} χ, con spinore
χ = (cos α/2, e^{iβ} sin α/2). A t = 0 si toglie il muro a z = d. Rivelatore: il piano z = L. Moto bohmiano con la
**corrente di Pauli** (termine convettivo + termine di spin (ħ/2m) ∇×(ψ†σψ)/ψ†ψ). Lo stato si fattorizza,
ψ = φ₀(x,y) e^{−iωt} W(z,t) χ, e le equazioni del moto (loro Eq. 13, scritte da noi in cartesiane per β = 0) sono

    Ẋ = −ω Y cos α,   Ẏ = Re(W′/W)(Z,t) sin α + ω X cos α,   Ż = Im(W′/W)(Z,t) + ω Y sin α .

* spin lungo la guida (α = 0 o π): il termine di spin non tocca Z, la distribuzione dei tempi = **flusso**
  Π_qf(τ) = Im(W* W′)(L,τ), coda ∝ 4(L/π)³/τ⁴, indipendente da ω;
* spin perpendicolare (α = π/2): distribuzione "esotica" con un **tempo massimo netto** τ_max: nessuna traiettoria
  arriva dopo. Per L = 100 d e ω = 10³ (unità DD): **τ_max ≈ 42,9** (≈ 1 ms per d = 50 µm); τ_max cresce
  linearmente con L, scende con ω e tende a una costante per ω ≫ 1; per ω → 0 l'effetto sparisce.
* esperimento proposto (mai fatto): elettroni in trappole di Penning, d = 50 µm, L ≈ 5 mm, ω ≈ 4,6·10⁷ rad/s
  (= 10³ DD), tempo unità m d²/ħ ≈ 21,7 µs. Il semiclassico dà invece una coda ∝ τ⁻².

**Das, *Detlef Dürr, arrival-time distributions, and spin in Bohmian mechanics* (arXiv:2309.15815, 2023)**
(letto): riassunto; la guida a piatti paralleli sopprime gli arrivi su un intervallo di α più largo: è la
geometria "più promettente". Nota 35: alcune condizioni al bordo assorbente darebbero dipendenza dallo spin, ma
perché lo spin è messo esplicitamente nella condizione.

**Goldstein, Tumulka, Zanghì, *Sci. Rep.* (febbraio 2024, PMC10866986), arXiv:2309.11835** (letto l'arXiv; il
numero dell'articolo nella rivista non l'ho verificato). Argomento semplice e forte: le statistiche di **qualunque** esperimento
quantistico vengono da una POVM, quindi P_n(t) = ⟨n|O(t)|n⟩ con O(t) operatore 2×2 sullo spin. Allora
P_n + P_{−n} = Tr O non dipende da n; con la simmetria assiale P↑ + P↓ = 2 P→. Le distribuzioni bohmiane di DD
hanno P↑ = P↓ ≠ P→ (fortemente): **non vengono da nessuna POVM**, quindi **nessun rivelatore** le può misurare;
l'errore di qualunque misura approssimata è ≥ ‖P→ − P↑‖/2. Con la simmetria per riflessione (chirale) un
rivelatore simmetrico non vede **nessuna** dipendenza dallo spin.

**Das & Aristarhov, arXiv:2312.01802 (2023)**: commento a GTZ (letto solo l'abstract: rispondono, il contenuto
non l'ho potuto leggere).

**Tumulka, regola del bordo assorbente (ABR)**: arXiv:1601.03715 ("Distribution of the time at which an ideal
detector clicks", poi *J. Stat. Phys.*), Schrödinger con ∂ₙψ = iκψ sul rivelatore (κ > 0 numero d'onda
caratteristico del rivelatore); densità del clic = (ħκ/m)|ψ|² sul bordo. Un'onda piana e^{ikz} è riflessa con
|(k−κ)/(k+κ)|²: κ = k la assorbe tutta. **Per Dirac**: arXiv:1601.04571 (letto): condizione **ideale**
n·α ψ = ψ ("ABCD"), la particella esce a velocità c; la versione semi-ideale (n·α + θβ)ψ = √(1+θ²) ψ ha limite non
relativistico la condizione di Pauli σ·∇φ = iκ (n·σ) φ con **κ = (2mc/ħ)(√(1+θ²) − θ)** (sua Eq. 66). Per θ = 0
(la condizione ideale) κ = 2mc/ħ, che diverge per c → ∞: Tumulka scrive che l'ABCD ideale **"non ha limite non
relativistico"**.

**Jozani & Tumulka, arXiv:2603.22044 (2026)** (letto l'HTML): proprio il set-up di DD con tre modelli di rivelatore
— potenziale immaginario, ABC scalare ∂_zψ = iκψ, ABC spinoriale (quella di Pauli sopra) — con κ = π. Trovano:
segni di **riflessione parziale** al rivelatore (oscillazioni), frazione rivelata 90 % (L = 100) e 96 % (L = 10)
per l'ABC scalare, distribuzioni **senza dipendenza dallo spin** (per l'ABC spinoriale "così debole da non essere
notabile"), ma l'ABC spinoriale dipende dalla larghezza della guida (𝔼T = A + B√ω). Confermano τ_max ≈ 43 di DD
(per ω = 500). **Quindi il punto (ii) del nostro compito (bordo assorbente sullo stato di DD) è già stato fatto:
lo ripetiamo solo come validazione.**

**Cosa resta aperto.** (1) Quale distribuzione misurerebbe davvero l'esperimento proposto (dipende dal
rivelatore, nessuno l'ha fatto). (2) Se il tempo bohmiano sia "fisico" (GTZ: non misurabile). (3) Quale κ abbia un
rivelatore vero: l'ABR lascia κ libero.

## 3. Perché la nostra regola c'entra, e cosa può dire di nuovo

* **I surfisti** della regola hanno v = c ψ†αψ/ψ†ψ: è la velocità di Bohm–Dirac. Nel limite non relativistico la
  corrente di Dirac contiene da sola (decomposizione di Gordon) il termine di spin ∇×(ψ†σψ)/2m: quello che DD
  aggiungono a mano a Schrödinger. In 2-D (y, z) con lo spin lungo x la regola deve far nascere l'effetto di DD
  **senza** che nessuno scriva il termine di spin.
* **Il rivelatore naturale** della regola: una cella che toglie tutto quello che ci entra (studio `tunnel_`,
  `docs/TEMPO_TUNNEL.md` §4.6). Niente rientra da oltre il piano z = L: le componenti che corrono verso −z sono
  zero lì, cioè α_z ψ = ψ: **è esattamente l'ABCD ideale di Tumulka**, realizzata come misura ripetuta a ogni tick.
  Riflette (1−v)/(1+v) (misurato in `tunnel_`: 0,07188 contro 0,07180 per v = 0,866). Nel limite non relativistico
  equivale all'ABR di Schrödinger con **κ = 2mc/ħ** (Tumulka Eq. 66 con θ = 0): per un elettrone lento è quasi uno
  specchio (assorbe ≈ 2v/c per urto).
* **GTZ in pratica**: ogni rivelatore fatto con la regola è lineare (una POVM) e la regola è (quasi) simmetrica per
  riflessione: la previsione è che **i clic non dipendano dallo spin**, mentre i surfisti senza rivelatore sì.

Cosa ci aspettiamo che sia **validazione**: le curve di DD (flusso, τ_max), il limite NR della regola, l'ABR di
Tumulka/Jozani. Cosa può essere **nuovo**: il rivelatore della regola come ABCD ideale con κ = 2mc/ħ fissato (zero
parametri) e la sua previsione numerica per l'esperimento proposto; i surfisti della regola che danno da soli il
termine di spin di Pauli. Cosa sarà probabilmente un **limite**: la regola non sa costruire un rivelatore
"morbido" (servirebbero atomi), e il suo rivelatore ideale è irrealistico per particelle lente.

## 4. Prove esplorative prima della pre-registrazione (dichiarate)

Per scegliere griglie e soglie ho fatto prove piccole (N ≤ 2000 traiettorie) prima di scrivere la
pre-registrazione. Le dichiaro perché hanno influenzato le scelte:
* riferimento NR (`arrivo_nr.py`): τ_max(L = 100, ω = 10³) ≈ 41,4 con 2000 traiettorie (il massimo del campione
  cresce col campione); a L = 10: ω = 40 → τ_max ≈ 6,0, ω = 10³ → 4,15;
* la velocità "di un sito" della regola (|R|² − |L|²)/ρ (`zz_core.dirac_velocity`) ha un errore O(a)∇ρ: vicino al
  muro dava v sbagliata di un fattore 10. Uso la **corrente conservata del reticolo centrata sul sito** (media dei
  due flussi di legame); con questa l'errore scende a O(a²);
* il filtro spettrale 2-D centrato sul quoziente di Rayleigh sbagliava stato (l'1 % di energia negativa sposta il
  quoziente di 0,009): ora filtro largo → Rayleigh → filtro stretto (due volte);
* rivelatore della regola (1-D, D = 20, L = 10): F(20) = 0,465 / 0,503 / 0,511 per a = 0,5 / 0,25 / 0,125, contro
  l'ABR con κ = 2mc/ħ: 0,517; ABR con κ = π: 0,965 (Jozani–Tumulka: 96 %).
