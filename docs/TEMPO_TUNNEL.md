# Tempo di tunnel e tempo di arrivo sulla regola zig-zag

Studio `tunnel_` (file: `lab/tunnel_*.py`, `lab/out/tunnel_*`). Regola: quella congelata di `lab/zz_core.py`
(1+1 D, `dirac1`): ogni pallina corre di una cella a c, poi ruota (flip θ = m c a/ħ) e prende la fase del
potenziale. Unità ħ = c = mₑ = 1 (lunghezza λ_C = 3,86·10⁻¹³ m, tempo ħ/mₑc² = 1,29·10⁻²¹ s).

*(documento in costruzione: §1 scritta prima delle corse; il resto dopo la pre-registrazione)*

---

## 1. Il problema e cosa si sa (ricerca, 13/10/2026)

**La domanda.** Una particella passa una barriera "per effetto tunnel" (energia sotto il bordo della
barriera). Quanto tempo ci ha passato dentro? In meccanica quantistica il tempo non è un operatore (Pauli,
1933: un operatore tempo autoaggiunto coniugato a un'energia limitata dal basso non esiste), quindi ci sono
**tante definizioni diverse**, che danno numeri diversi. Le principali:

| tempo | definizione | cosa dice | chi |
|---|---|---|---|
| **fase (Wigner, gruppo)** | ritardo del picco: τ_φ = ħ ∂(arg t + kL)/∂E | per barriere spesse **satura** (non cresce con L): "effetto Hartman"; la velocità apparente L/τ_φ supera c | Wigner 1955; Hartman, *J. Appl. Phys.* 33, 3427 (1962) |
| **permanenza (dwell)** | τ_D = ∫_barriera \|ψ\|² dx / flusso incidente | tempo medio passato dentro, senza distinguere chi passa da chi torna indietro; anche lui satura | Smith 1960; Büttiker 1983 |
| **Larmor (Baz', Rybachenko, Büttiker)** | un orologio: un piccolo campo magnetico solo dentro la barriera fa ruotare lo spin; l'angolo misura il tempo. Due componenti: τ_y (rotazione nel piano) = −ħ ∂arg t/∂V, τ_z (inclinazione) = −ħ ∂ln\|t\|/∂V | per una barriera simmetrica τ_y = τ_D; τ_z per barriere spesse → m L/(ħκ) (tempo "semiclassico" di Büttiker–Landauer, cresce con L) | Baz' 1967; Büttiker, *PRB* 27, 6178 (1983); Büttiker e Landauer, *PRL* 49, 1739 (1982) |
| **valore debole** | τ_y − i τ_z = valore debole (complesso) del "tempo passato dentro" sui soli trasmessi | parte reale = tempo vero, parte immaginaria = disturbo della misura | Steinberg, *PRL* 74, 2405 (1995) |
| **Bohm** | traiettorie che seguono la corrente (v = j/ρ); in 1-D non si incrociano | i trasmessi sono quelli partiti davanti; tempo finito, mai sopra c nella versione di Dirac | Leavens (anni '90); Hauge e Støvneng, *RMP* 61, 917 (1989, rassegna) |
| **attoclock** | ionizzazione di atomi da laser ellittico: angolo dell'elettrone uscito | in idrogeno tempo di tunnel compatibile con **zero** (< 1,8 as) | Torlina et al., *Nat. Phys.* 11, 503 (2015); Sainadh et al., *Nature* 568, 75 (2019) |

**Le misure.** Ramos, Spierings, Racicot, Steinberg, *Nature* 583, 529 (2020) (arXiv 1907.13523, letto
l'abstract online): atomi di ⁸⁷Rb condensati (energia ~1 nK), barriera ottica spessa 1,3 µm, orologio di
Larmor fatto con una transizione Raman solo dentro la barriera: **tempo di attraversamento 0,62(7) ms**
(un riassunto di corso, Illinois PHYS596 2022, riporta τ_y = 0,61(7) ms al picco della distribuzione di
velocità e la previsione "controintuitiva" che **energie più basse danno tempi più corti**). Altezza della
barriera e larghezza gaussiana esatta: **non sono riuscito a leggerle** (il testo completo non era
raggiungibile: limite di accesso del proxy); le uso solo come forma, vedi §6. Con la luce: Steinberg, Kwiat,
Chiao, *PRL* 71, 708 (1993) (fotoni, ritardo ~1,5 fs "superluminale"); Spielmann et al., *PRL* 73, 2308
(1994) (impulsi in specchi a strati).

**Relatività (il nostro terreno).** Con l'equazione di Dirac la storia cambia poco ma in modo importante:
* Krekora, Su, Grobe, *PRA* 63, 032107 (2001) (abstract letto online): pacchetti di Dirac attraverso una
  barriera; picco e fronte possono sembrare superluminali, ma **Dirac non viola la causalità** (Schrödinger sì,
  perché non ha un cono di luce). Anche *PRA* 64, 022105 (2001) e *PRL* 92, 040406 (2004) (Klein).
* Petrillo e Janner, *PRA* 67, 012110 (2003); De Leo e Rotelli, *EPJ C* 51, 241 (2007); De Leo, *J. Phys. A*
  46, 155306 (2013) (letto online): nel limite opaco il tempo di fase dà velocità apparenti oltre c; De Leo
  conclude che il tempo di fase "non è il tempo di tunnel giusto" per l'effetto filtro (passano solo le
  componenti più energetiche) e propone come lavoro futuro le **correzioni relativistiche a Larmor**.
* Li, Liang, Kobe (arXiv quant-ph/0111117, letto l'abstract): Larmor relativistico per una particella neutra di
  spin ½: il tempo di precessione è uguale al dwell (come in Schrödinger).
* Flores e Galapon (arXiv 2207.09040, 2022, letto l'abstract): operatore "tempo di arrivo" relativistico per
  spin 0: le componenti sotto barriera passano in tempo zero.
* Winful, *Phys. Rep.* 436, 1 (2006) (non verificato online): il ritardo di gruppo è un tempo di "accumulo e
  rilascio" di energia, non di transito; la saturazione viene dalla saturazione dell'energia immagazzinata
  (niente superluminalità vera).

**Tempo di arrivo.** Nessun operatore accettato (Pauli). Candidati: flusso j(X,t) nel punto del rivelatore
(può essere negativo: *backflow*); regola del bordo assorbente di Tumulka (2016–2022, non verificato online:
per Dirac la condizione al bordo "solo onda uscente" n·α ψ = ψ); prima-arrivata sui cammini a scacchiera di
Feynman: Alonso, Sala Mayato, Leavens (arXiv quant-ph/0206071, letto l'abstract): con i cammini di Feynman la
distribuzione di prima-arrivata **non** si definisce in modo coerente per l'interferenza tra prime e
successive arrivate (tranne vicino a v → c).

**Il conteggio delle svolte.** Jacobson e Schulman, *J. Phys. A* 17, 375 (1984) (riferimento verificato, il
contenuto no): nella scacchiera di Feynman il numero di svolte è legato al passaggio relativistico →
non relativistico (in tempo immaginario le svolte sono un processo di Poisson con ritmo ∝ m c²/ħ).

**Cosa è aperto.** (i) Quale tempo è "quello vero" (probabilmente nessuno: domande diverse, risposte
diverse); (ii) se l'effetto Hartman dica qualcosa sul transito o solo sul rimodellamento del pacchetto
(Winful: il secondo); (iii) un tempo di arrivo con una distribuzione di probabilità accettata; (iv) cosa
misura davvero l'attoclock.

**Cosa può dire la nostra regola** (e cosa no). La zig-zag ha due cose che Schrödinger non ha: un **cono di
luce rigido** (una pallina non può fare più di una cella per tick) e un **orologio interno naturale**: le
rotazioni (flip) dei cammini a scacchiera. Non ha (per costruzione) niente di nuovo da dire sulle definizioni
standard: quelle sono proprietà dell'equazione di Dirac che la regola ritrova nel limite a → 0 (sarà
VALIDAZIONE). Il punto dove potremmo dire qualcosa di non banale è la domanda 3: **un orologio che conta le
rotazioni** dentro la barriera. Lo definiamo così: ogni cammino della scacchiera che attraversa la barriera fa
un certo numero f di rotazioni dentro la barriera; la "media" sui cammini trasmessi pesata con le ampiezze
(valore debole) è f_w = Σ f A_cammino / Σ A_cammino = ∂ln t/∂ln z, con z peso extra per ogni rotazione dentro.
Si dimostra (identità esatta sul reticolo, verificata a 10⁻⁹) che
  f_w = sinθ cosθ ∂ln t/∂θ_B + sin²θ n_w,
cioè contare le rotazioni equivale a un **orologio di Larmor accoppiato alla massa** (una piccola variazione di
m solo dentro la barriera) invece che al potenziale. Nel continuo questo orologio legge ψ̄ψ invece di ψ†ψ: è il
**tempo proprio** (per una particella libera: T_m = (L/v)(m/E), la dilatazione del tempo). Non abbiamo trovato
in letteratura questo orologio applicato al tunnel (ricerca limitata: 4 ricerche web; "relativistic
corrections to Larmor" sono indicate come lavoro futuro da De Leo 2013).
