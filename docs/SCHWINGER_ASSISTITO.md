# Schwinger assistito: coppie dal vuoto sotto soglia con due campi insieme

Studio dei "problemi aperti". Prefisso dei file: `assist_`. Codice in `lab/assist_*.py`, risultati in
`lab/out/assist_*`. Pre-registrazione: `lab/out/assist_prereg.json` (generata da `lab/assist_prereg.py`).

*(Documento in costruzione: questa è la fase 1, la ricerca. Le fasi 2–4 si aggiungono sotto.)*

---

## 1. Il problema, in parole semplici

Il vuoto della regola zig-zag è il "mare" pieno di elettroni a energia negativa (la banda bassa). Un campo
elettrico molto forte può strappare un elettrone dal mare: nasce una coppia elettrone + positrone. È
l'**effetto Schwinger**. Il guaio: serve un campo enorme, E_S = m²c³/(eħ) ≈ 1,3·10¹⁸ V/m. La probabilità va
come exp(−π E_S/E): con E = E_S/10 è già e^(−31) ≈ 3·10⁻¹⁴, con i laser di oggi (E ≈ 10⁻³ E_S) è e^(−3000),
cioè zero. **Non è mai stato visto con un campo di laser.**

L'idea (Schützhold, Gies, Dunne 2008): sommare al campo **forte e lento** un campo **debole e veloce**. Il campo
veloce da solo non basta (i suoi "fotoni" hanno energia ω ≪ 2m), quello forte da solo nemmeno; insieme
l'esponente crolla. Immagine: il campo forte inclina il mare e crea una "barriera" da attraversare in
tunnel; il campo veloce dà un calcio di energia a metà strada e accorcia il tunnel.

## 2. Cosa è noto (letteratura)

Fonti lette online in questa sessione (arXiv/ar5iv), salvo dove scritto "a memoria".

* **Schützhold, Gies, Dunne, PRL 101, 130404 (2008)** (arXiv:0807.0754). Campo
  E(t) = E/cosh²(Ωt) + ε/cosh²(ωt), con ε ≪ E e Ω ≪ ω ≪ m. Parametro di Keldysh **combinato**
  γ = mω/(eE) (frequenza del campo veloce, intensità di quello forte). Con il metodo delle traiettorie
  istantone (worldline instanton) trovano una **soglia netta**: per γ < γ_crit = π/2 l'esponente resta
  quello di Schwinger; per γ > π/2 diventa
  𝒜 = (m²/eE)·[2 arcsin(π/2γ) + (π/2γ²)·√(4γ² − π²)]  (eq. 15),
  **indipendente da ε** (ε entra solo nel prefattore). Per γ grande 𝒜 → 2πm/ω.
* **Dunne, Gies, Schützhold, PRD 80, 111301 (2009)** "Catalysis of Schwinger vacuum pair production" (a memoria):
  stessa idea, istantoni per ε finito.
* **Orthaber, Hebenstreit, Alkofer, PLB 698, 80 (2011)** (arXiv:1102.2182): numerica con l'equazione cinetica
  quantistica (spinori, campo omogeneo) per due impulsi di Sauter, ℰ₁ = 0,25 E_S, τ₁ = 10⁻⁴ eV⁻¹ (≈ 51/m),
  ℰ₂ = 0,025 E_S, τ₂ = 7·10⁻⁶ eV⁻¹ (≈ 3,6/m). Trovano un aumento del numero di coppie e spettri in impulso
  strutturati. Nel testo **non ci sono numeri** (solo figure): non posso confrontare cifra per cifra.
  Nota: con la loro definizione γ = m/(eℰ₁τ₂) ≈ 1,1 < π/2, cioè sono **sotto** la soglia degli istantoni.
* **Linder, Schneider, Sicking, Szpak, Schützhold, PRD 92, 085009 (2015)** (arXiv:1505.05685): **1+1 dimensioni**,
  E(t) = E₁ + E₂·g(ω₂t) con g = 1/cosh² (Sauter), gaussiana, cos. Con WKB e istantoni:
  - Sauter: γ_crit = π/2 indipendente da ε;
  - gaussiana: γ_crit ~ √|ln ε|; coseno: γ_crit ~ |ln ε| (la soglia dipende dalla forma dell'impulso!);
  - coseno al primo ordine in ε (eq. 37): esponente = (πE_S/E₁)·[1 − 2ε cos(pγ) I₁(γ)/γ].
  Nessuna numerica esatta in quell'articolo.
* **Schneider, Schützhold, JHEP (2016)** (arXiv:1603.00864) "Prefactor in the dynamically assisted Sauter-Schwinger
  effect": prefattore con gli istantoni (QED scalare 1+1), confronto con l'integrazione numerica di Riccati
  "eccellente" a E = E_S/5, si stacca a E_S/3.
* **Torgrimsson, Schneider, Oertel, Schützhold, JHEP 06, 043 (2017)** (arXiv:1703.09203): sviluppo in potenze del
  campo debole (ε): per l'impulso di Sauter **ordine zero + primo ordine** bastano e concordano con la numerica;
  per gaussiana/coseno gli ordini alti possono dominare. Numeri della loro numerica: non estratti (sez. 6 non
  leggibile dal testo).
* **Varianti e verifiche** (a memoria, non verificate online): Akal, Villalba-Chávez, Müller PRD 90, 113004 (2014)
  (due frequenze); Hebenstreit, Fillion-Gourdeau PLB 739 (2014) (ottimizzazione della forma); Otto, Seipt, Blaschke,
  Kämpfer, Smolyansky (2015) (assistenza dinamica con soglia multifotone).
* **Identità utile** (a memoria, standard: Schmidt et al. 1998, Hebenstreit tesi 2011): per un campo omogeneo
  lineare l'equazione cinetica della QED a 3+1 dimensioni con spinori è **la stessa** dell'equazione di Dirac
  1+1 con m → m_⊥ = √(m² + p_⊥²), più un fattore 2 di spin. Quindi la 1+1 non è un giocattolo: è la fetta
  p_⊥ = 0 della risposta vera. Cambiano i prefattori (integrale su p_⊥), non gli esponenti.
* **Riferimenti per la validazione:** Schwinger/Nikishov in 1+1: numero medio di coppie per lunghezza e tempo
  (eE/2π)·exp(−πm²/eE). Impulso singolo di Sauter E sech²(t/τ): spettro esatto (Narozhny–Nikishov 1970),
  f(p) = [cosh(2πeEτ²) − cosh(πτ(ω₊ − ω₋))] / [2 sinh(πτω₊) sinh(πτω₋)], ω± = √(m² + (p ± eEτ)²); copre tutto il
  passaggio tunnel (γ = m/(eEτ) ≪ 1) → multifotone (γ ≫ 1). Campo oscillante: risonanze a n fotoni
  nω = 2q̄ (massa "vestita" dal campo), a p = 0 solo n dispari per la simmetria di coniugazione di carica
  (Mocken, Ruf, Müller, Keitel PRA 81, 022122 (2010), a memoria).

**Cosa è aperto.** (1) L'effetto non è mai stato osservato. (2) Le formule semiclassiche valgono per ω ≪ m e
ε piccolo ma "non troppo": vicino alla soglia e per ω ~ m servono numeri esatti, che esistono solo per alcuni
parametri e in forma di figure. (3) Quanto la forma dell'impulso conti (soglia netta solo per Sauter).

**Cosa possiamo fare noi, onestamente.** In un campo **omogeneo** la regola zig-zag 1+1 dà il numero di coppie
**esatto del modello**, modo per modo, come già in R1 (Schwinger 0,98–1,04, Sauter 0,992–0,996). Qui non c'è
fisica nuova da scoprire nella QED continua: la QED 1+1 in campo esterno è un problema risolto in linea di
principio (equazione cinetica). Quello che possiamo dare:
* **validazione**: la regola ritrova Schwinger, il passaggio di Keldysh e le regole multifotone;
* **esame delle formule semiclassiche** (soglia γ_crit, esponente) con numeri esatti, convergenti in a → 0;
* **numeri dove la semiclassica non vale** (vicino alla soglia, ω ~ m): previsioni esatte del modello 1+1,
  che però **coincidono** con quelle della QED 1+1 continua (la regola converge a Dirac). Non sono "scoperte":
  sono numeri controllati, confrontabili con la numerica di letteratura quando i parametri coincidono.
