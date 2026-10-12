# Sincronia: l'universo come un campo di metronomi

Laboratorio: `lab/kz_*.py` (motore `lab/kz_core.py`). Numeri: `lab/out/kz_*.json` / `.log`. Figura: `lab/out/kz_figure.png`.
Pre-registrazione: `lab/kz_prereg.py` → `lab/out/kz_prereg.json`, **committata (d8a3c2e) prima delle corse**.
Prima di scriverla ho fatto solo prove di tempi e due "corse di sviluppo" con un seme, dichiarate nel file
(compreso l'unico numero simile ad alpha che avevo visto: una pendenza a due punti di circa 0,6).

## L'idea, in parole semplici
Immagina un'enorme stanza piena di **metronomi**, uno per ogni punto dello spazio. All'inizio battono a caso,
ognuno per conto suo (fa "caldo": tutto è agitato). Poi la stanza si raffredda e i metronomi vicini cominciano
a mettersi d'accordo, a sincronizzarsi. Ma la notizia "andiamo a tempo così" viaggia a velocità finita: zone
lontane si mettono d'accordo ciascuna a modo suo. Dove tre o più zone si incontrano con fasi incompatibili, la
sincronia **non si può completare**: resta un "nodo" nel campo, un punto (in 2D) o una linea (in 3D) attorno a cui
la fase gira di un giro intero. L'idea di Valerio: **le particelle sono questi nodi congelati**.

Questo meccanismo esiste davvero ed è famoso: si chiama **Kibble–Zurek (KZ)**. Kibble lo pensò per l'universo
primordiale, Zurek per l'elio superfluido. Fa una previsione precisa: più raffreddi piano, meno difetti restano.
Se raffreddi in un tempo τ_Q, la densità di difetti scende come

  n ∝ τ_Q^(−α),   α = (D − d)·ν / (1 + z·ν)

dove D = 3 (dimensioni dello spazio), d = 1 (i difetti sono linee), ν dice quanto cresce la "zona d'accordo"
vicino alla temperatura critica, z quanto rallenta il sistema lì vicino. Ho verificato la formula: è quella
standard (la zona congelata ha taglia ξ̂ ∝ τ_Q^(ν/(1+zν)), e le linee hanno densità 1/ξ̂², quindi α = 2ν/(1+zν)).

* **Campo medio** (la teoria "da lavagna"): ν = 1/2, z = 2 → **α = 0,5**.
* **Vera classe 3D XY** (quella dei metronomi in 3D, e dell'elio): ν = 0,6717. Con la nostra dinamica
  (tipo "modello A", z ≈ 2) → **α = 0,573**. (Per l'elio vero, modello F, z = 3/2 → α ≈ 0,67: non lo testiamo.)

La **sfida**: misurare α così bene da dire "0,573, non 0,5". Sarebbe una previsione che richiede calcolo vero,
non solo la teoria da lavagna. Le prove "usa e getta" di prima avevano dato −0,22 (anello 1D), −0,62 (vortici 2D)
e −0,95 in 3D, ma quest'ultima contava troppo tardi (i difetti si erano già mangiati a vicenda).

## Il simulatore
* Reticolo cubico L×L×L di metronomi (angoli θ), energia −Σ cos(θ_i − θ_j) tra vicini, bordi periodici.
* Dinamica: Metropolis a scacchiera (ogni metronomo prova a ruotare di ±1 radiante e accetta con la regola
  termica). È "modello A": niente conservazione, z ≈ 2. Un "sweep" = un'unità di tempo.
* Raffreddamento lineare: T(t) = Tc(1 − t/τ_Q), con Tc = 2,20184 (valore noto), partendo da 1,2·Tc.
* Conto dei difetti: per ogni quadratino (plaquette) sommo i salti di fase attorno; se fanno un giro intero,
  una linea di vortice lo attraversa.

**Il problema che ho trovato prima di registrare** (dichiarato): subito sotto Tc il campo è pieno di piccoli
anelli di vortice **termici** (rumore all'equilibrio, grandi come un sito). Contati "crudi", sono 0,5–0,6 per sito e
quasi non dipendono da τ_Q. Soluzione registrata: prima di contare, liscio una **copia** del campo a T = 0
("gradient flow") per un tempo b·t̂, e conto al tempo a·t̂ dopo Tc, con t̂ = τ_Q^0,573 (il tempo di congelamento
della teoria). Siccome con z = 2 il tempo t̂ va come ξ̂², lisciare per un tempo ∝ t̂ cancella gli anellini senza
rovinare la scala giusta: se KZ vale, il conto è proprio ∝ ξ̂^−2.

## K1 – L'esponente α (il test principale)
Corsa principale registrata: L = 48, 8 semi per τ_Q, τ_Q = 64…4096 (1,8 decadi, 7 valori), a = 2, b = 0,5.
Fit di ln n contro ln τ_Q, intervallo di confidenza al 95% col bootstrap sui semi.

| τ_Q | 64 | 128 | 256 | 512 | 1024 | 2048 | 4096 |
|---|---|---|---|---|---|---|---|
| linea di vortice per sito (L=48) | 0,0481 | 0,0324 | 0,0238 | 0,0167 | 0,0124 | 0,0084 | 0,0056 |

**α = 0,507, IC 95% [0,493 ; 0,523].**

Criterio pre-registrato: l'IC deve **escludere 0,5 e includere 0,573**. Qui fa esattamente il contrario: include
0,5 ed esclude 0,573. **K1: FAIL.** Quindi ROBUST non si pone (FAIL).

Controlli (tutti registrati):

| variante | α | IC 95% |
|---|---|---|
| 9 combinazioni a ∈ {1,2,4}, b ∈ {0,25; 0,5; 1} | da 0,503 a 0,527 | (scarto max dal principale 0,020) |
| solo τ_Q ≥ 256, L = 48 | 0,519 | [0,493 ; 0,546] |
| L = 32 | 0,485 | [0,451 ; 0,521] |
| L = 64, τ_Q ≥ 256 (4 semi) | 0,499 | [0,464 ; 0,544] |
| conteggio "da campo medio", t̂ = τ_Q^0,5 | 0,449 | [0,438 ; 0,461] |

Cosa ci dicono:
* **Taglia finita**: L = 32, 48, 64 danno lo stesso numero entro gli errori (a τ_Q = 4096 la distanza media tra
  linee è 13 siti; con L = 64 la densità è 0,0062 ± 0,0005 contro 0,0056 a L = 48). Non è un artefatto di scatola.
* **Quando/quanto lisciare** (a, b): conta poco, ±0,02.
* **La regola per scegliere l'istante del conteggio conta molto**: se conto a t ∝ τ_Q^0,5 invece di τ_Q^0,573,
  α passa da 0,507 a 0,449. Cioè l'incertezza "di protocollo" è circa **0,06**, quasi quanto la differenza che
  volevamo vedere (0,573 − 0,5 = 0,073). Il motivo: dopo il congelamento i difetti continuano lentamente ad
  annichilarsi, e spostare l'istante del conteggio cambia la pendenza.
* Le pendenze locali (tra τ_Q vicini) sono 0,57 / 0,45 / 0,51 / 0,43 / 0,57 / 0,58: rumorose. Le ultime due
  sembrano salire verso 0,57, ma è un'osservazione **a posteriori** e L = 64 non la conferma. Da non usare.

**Verdetto onesto**: il simulatore misura la legge di potenza KZ (pulita, 1,8 decadi, densità che scende di 9 volte)
con α ≈ 0,51 ± 0,02 (statistica) ± ~0,06 (protocollo). **Non è una previsione oltre il campo medio**: il numero
coincide col campo medio, e comunque l'errore sistematico è troppo grande per separare 0,5 da 0,573.

## K2 – Il simulatore conosce la sua classe? (controllo all'equilibrio)
Algoritmo a cluster (Swendsen–Wang con il trucco di riflessione di Wolff) + over-relaxation, L = 8, 12, 16, 24, 32
esattamente a Tc. Le derivate d ln⟨|m|⟩/dβ e d ln⟨m²⟩/dβ crescono come L^(1/ν).

* **ν = 0,677 ± 0,018** (entrambi gli stimatori 0,677; togliendo L = 8 viene 0,700 ± 0,034). Atteso 0,6717.
  Lontanissimo dal campo medio 0,5.
* Incrocio delle curve di Binder U4(T): T = 2,1989 (8–16), 2,2044 (12–24), **2,2041 (16–32)**, cioè 0,10% dal valore
  noto 2,20184. Il valore di U4 all'incrocio (1,252) è vicino a quello di letteratura (~1,243).

**K2: PASS.** All'equilibrio il simulatore "sa" di essere 3D XY. Quindi il fallimento di K1 non viene da un
simulatore sbagliato: viene dal modo di contare i difetti fuori equilibrio (e/o da correzioni di scala alle τ_Q
raggiungibili).

## K3 – I difetti nati si muovono?
Qui cambio dinamica: dai metronomi "nella melassa" (primo ordine) passo a metronomi con **inerzia**
(θ'' = −γθ' + Σ sin(θ_j − θ_i)), dove le onde viaggiano a velocità c = 1 (come la luce nel nostro mondo-giocattolo).

* **K3a, coppia vortice–antivortice a distanza 24 (pre-registrata): FAIL.** Non si annichila: in 300 unità di tempo
  la distanza passa solo da 24 a 22. L'energia si conserva benissimo (1,4·10⁻⁵), quindi non è un errore numerico:
  i vortici sono **incastrati nel reticolo** (barriera di Peierls–Nabarro). Il nucleo di un vortice nei puri
  metronomi è grande come un solo sito, quindi "sente" i gradini della griglia; l'attrazione 2π/d a d = 24 è troppo
  debole per farli saltare.
* **A posteriori (non registrato)**: provo distanze più piccole. d = 6, 8, 12 si annichilano (t ≈ 10, 15, 85), d = 16 resta
  bloccata. Rifaccio tutta l'analisi a **d = 8**: si annichila a t = 13; velocità massima dei vortici ~0,4 (< c);
  energia conservata a 7·10⁻⁴; dopo l'annichilazione il **96%** dell'energia della coppia esce come un **anello d'onda**
  il cui fronte corre a **0,96 c** (il raggio del 90% dell'energia), e l'energia si divide metà cinetica e metà
  elastica (onde vere). Tutti i criteri di K3a passano, ma solo per d = 8, scelto dopo.
* **K3b, "calcio" alla coppia** (velocità iniziale 0,5 data con la trasformazione di Lorentz, solo da riportare):
  la coppia **non parte**. Il centro si muove a 0,026 nei primi 20 tempi e poi si ferma; l'energia del calcio se ne va
  in onde. Quindi **i vortici dei metronomi su reticolo non sono "particelle alla Lorentz"**: il reticolo li frena.
  (Il kink di sine-Gordon in 1D passava perché è largo molti siti.)
* **K3c, i difetti nati da un raffreddamento 2D** (L = 128, τ_Q = 256, fermato a T = 0,2, poi dinamica inerziale con
  γ = 0,02): da 68 vortici a 14 in t = 200 (−79%), e intanto compaiono onde (energia cinetica fino a 800). **PASS.**
  Qui le coppie vicine si trovano e si annichilano.
* **K3d, 3D** (configurazione finale di K1, L = 48, τ_Q = 1024): la lunghezza delle linee scende da 0,0099 a 0,0018 per
  sito in t = 100. Le linee si muovono e si accorciano.

## Bilancio
| test | esito |
|---|---|
| K1 α oltre il campo medio | **FAIL**: α = 0,507 [0,493 ; 0,523]; include 0,5, esclude 0,573 |
| K2 ν e Tc all'equilibrio | **PASS**: ν = 0,677 ± 0,018; Tc entro 0,1% |
| K3a coppia d = 24 | **FAIL** (bloccata dal reticolo); d = 8 a posteriori: passa tutto |
| K3b calcio di Lorentz | non si muove (solo riportato) → dinamica non "alla Lorentz" |
| K3c difetti 2D nati dal quench | **PASS** (−79%) |

**L'idea dei metronomi regge nella parte qualitativa**: raffreddando si formano davvero linee di vortice congelate,
tante quanto dice la legge di KZ (legge di potenza netta su quasi due decadi), e poi si muovono e si annichilano
liberando onde. **Ma non abbiamo una previsione oltre la teoria da lavagna**: il nostro α è indistinguibile da 0,5.

## Debiti (cosa resta da fare)
1. **Protocollo di conteggio**: la scelta dell'istante del conteggio vale ~0,06 in α. Serve un metodo che non
   dipenda da un esponente scelto prima, per esempio misurare la lunghezza di correlazione congelata ξ̂ direttamente,
   o un "collasso di dati" su molti istanti (n·ξ̂² funzione di t/t̂) fittando insieme l'esponente.
2. **Più decadi**: τ_Q fino a 10⁵ e L ≥ 96 (ore di calcolo, non minuti) per vedere se le pendenze locali salgono
   davvero verso 0,57.
3. **Vortici "lisci"**: nei puri metronomi il nucleo è largo un sito e il reticolo li incastra. Un campo con ampiezza
   (Ginzburg–Landau / "Higgs", nucleo largo molti siti) è necessario per avere difetti che si muovono come particelle.
4. **Elio vero** (modello F, z = 3/2, α ≈ 0,67): richiede una dinamica che conservi energia e quantità di moto.
5. K3a/K3b a d = 24 sono falliti come registrati; il successo a d = 8 è a posteriori.

## File
* `lab/kz_core.py` motore; `lab/kz_prereg.py` pre-registrazione; `lab/kz_k1.py` (corse K1, coda `lab/kz_run.sh`),
  `lab/kz_k2.py`, `lab/kz_k3.py` (coda `lab/kz_run2.sh`); `lab/kz_summary.py` analisi K1; `lab/kz_figure.py`.
* `lab/out/kz_k1.json`, `kz_k2.json`, `kz_k3.json`, log `kz_k*.log`, figura `kz_figure.png`.
