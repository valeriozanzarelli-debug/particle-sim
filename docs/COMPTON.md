# Compton: un fotone contro un elettrone libero

Programma: `cargo run --bin compton` (codice in `src/compton/`). Controllo indipendente in Python:
`lab/compton_check.py`. Tutti i numeri qui sotto vengono da esecuzioni reali con seed fissi.

---

## 1. In una frase

Un pacchetto-onda di fotone (palline-onda, la nostra regola) arriva su un pacchetto-onda di elettrone fermo.
L'onda che esce è una **somma di coppie**: per ogni direzione θ del fotone c'è un pezzo d'onda del fotone che
va lì **e** un pezzo d'onda dell'elettrone che rincula nella direzione "coniugata". Solo un pezzettino minuscolo
dell'onda viene deviato; il resto passa dritto. I due surfisti (fotone ed elettrone) scelgono **una** coppia,
oppure niente: il fotone passa e l'elettrone resta fermo. Dopo l'urto i due sono **entangled**.

## 2. Che cosa succede, passo per passo

**Prima.** L'onda del fotone: palline-onda che volano a c lungo il fascio (asse z), più fitte sulle creste
(distanti λ). L'onda dell'elettrone: un pacchetto fermo (palline azzurre). Dentro ognuna c'è il suo surfista.
Il surfista dell'elettrone è **fermo**: è il surfista di Schrödinger/Dirac di un pacchetto a energia positiva
in quiete (in `docs/SURFISTA_DIRAC.md` si vede che senza componenti a energia negativa non c'è Zitterbewegung).

**L'urto: l'onda a due particelle.** Dopo l'incontro l'onda è

    Ψ(x_γ, x_e) ≈ Ψ_passa(x_γ, x_e)  +  Σ_θ,φ  f(θ, φ) · [onda del fotone verso (θ, φ)] ⊗ [onda dell'elettrone verso (θₑ, φ + π)]

È la stessa forma "piccola somma di prodotti" di `docs/MULTIFOTONE.md`, con tanti termini (uno per direzione).
Per ogni direzione del fotone, **energia e impulso conservati** fissano tutto il resto (questo è derivato, esatto):

* lunghezza d'onda diffusa: **λ' − λ = (h/mₑc)(1 − cos θ)**, h/mₑc = 2,42631 pm;
* energia del fotone: E' = E / (1 + (E/mₑc²)(1 − cos θ));
* elettrone: prende esattamente l'impulso mancante, p_e = p − p'; il suo angolo soddisfa
  **tan θₑ = cot(θ/2) / (1 + E/mₑc²)**, la sua energia cinetica è T = E − E'.

**Quanta** onda va in ogni direzione, invece, **non** lo ricaviamo dalle palline: lo prendiamo dalla QED, la
formula di **Klein–Nishina** (1929), con la dipendenza dalla polarizzazione:

    dσ/dΩ = (r_e²/2) (E'/E)² [E'/E + E/E' − 2 sin²θ cos²φ]      (polarizzazione lineare lungo x)
    dσ/dΩ = (r_e²/2) (E'/E)² [E'/E + E/E' − sin²θ]               (non polarizzata)

È un **INGRESSO**, come lo sono h, c, mₑ e r_e. Integrata dà σ_T = 6,652·10⁻²⁹ m² a bassa energia (Thomson) e meno
ad alta energia.

**Quanto è piccolo il pezzo deviato.** La probabilità che un fotone del pacchetto venga deviato è
σ × (densità di probabilità della distanza trasversa relativa b in b = 0) = σ / (2π s²), con s la larghezza
trasversa relativa dei due pacchetti. Con un fotone largo 10 λ e un elettrone largo 1 nm (INGRESSI, scelti da noi):
p ≈ 6,6·10⁻¹² a 17,48 keV, 4·10⁻¹⁹ per la luce visibile. Il resto dell'onda passa dritto, e lo si vede andare avanti.

**I surfisti scelgono una coppia (la regola).** Quando il surfista-fotone raggiunge il piano del surfista-elettrone:

1. si guarda la sua distanza trasversa **b** dal surfista-elettrone e il suo "posto in fila" nel flusso entrante,
   u = F(b) = 1 − exp(−b²/2s²) (la frazione di probabilità relativa entro il raggio b);
2. se **u < p** (cioè b cade nel disco di area σ attorno all'elettrone, raggio √(σ/π) ≈ 4,5 fm) il fotone entra
   nella parte deviata: il rango v = u/p sceglie θ (v = 0 → θ = 180°, cioè b più piccolo → più all'indietro),
   l'angolo di b attorno all'asse sceglie l'azimut φ (con la sua distribuzione condizionata di Klein–Nishina);
   l'elettrone prende la direzione e la velocità coniugate;
3. altrimenti il fotone passa dritto e l'elettrone resta fermo.

È la stessa idea della regola base al bordo (`docs/REGOLA_BASE.md` v2): il surfista **conserva il suo posto in fila
nel flusso**. Siccome i surfisti nascono con densità |ψ|², questa mappa conserva la probabilità (equivarianza): gli
angoli escono con la statistica di |ampiezza|² = Klein–Nishina. **Onestà:** la vera traiettoria di Bohm nell'onda
completa a 6 dimensioni sarebbe un'altra mappa che conserva la probabilità: stesse frequenze, singole traiettorie
diverse. Non abbiamo risolto quell'onda; la mappa "a ranghi" è una scelta, non una derivazione.

**Dopo.** I due surfisti volano via, il fotone a c, l'elettrone a βc lungo p_e. Sono entangled: la direzione
dell'uno fissa quella dell'altro (linea rosa). Ogni coppia somma **esattamente** l'impulso iniziale, evento per
evento. Le linee rosa tenui uniscono altre coppie dell'onda: l'onda le contiene tutte, i surfisti ne seguono una.

## 3. Che cosa è derivato e che cosa è ingresso

| cosa | da dove |
|---|---|
| λ' − λ, E', angolo ed energia dell'elettrone | **derivato**: conservazione di energia e impulso (relatività), esatto |
| quanta onda in ogni direzione, polarizzazione | **ingresso**: Klein–Nishina (QED) |
| σ totale | ingresso (integrale di Klein–Nishina; tabella numerica = formula chiusa a 10⁻¹²) |
| larghezza del pacchetto-fotone (10 λ) e dell'elettrone (1 nm) | **ingresso**, scelto da noi; cambia solo p |
| la mappa rango → angolo (b piccolo → indietro) | **scelta** della regola; cambia le singole storie, non le frequenze |
| costanti | CODATA 2018: mₑc² = 510,99895 keV, hc = 1239,841984 keV·pm, r_e = 2,8179403262 fm |

## 4. Che cosa si vede (`cargo run --bin compton`)

* **Sezione nel piano dell'urto** (vista iniziale): una fetta sottile dell'onda nel piano (fascio, fotone, elettrone).
  Il disegno è ruotato attorno al fascio in modo che il piano dell'urto guardi la telecamera (solo un cambio di
  punto di vista). Le creste del guscio diffuso sono anelli distanti λ'(θ): **all'indietro le bande sono più larghe**.
  Colore = allungamento della lunghezza d'onda: giallo chiaro (niente) → arancio → rosso → magenta (2 h/mₑc).
* **Cono di rinculo**: la nuvola azzurra sono i pezzi d'onda dell'elettrone che rincula (uno per ogni direzione del
  fotone); le curve azzurre sono il loro fronte. Nella vista 3-D (a metà o intera) si vede la forma tridimensionale.
* **L'onda che passa**: il pacchetto entrante continua dritto (quasi tutta l'onda).
* **Surfisti**: pallina bianca (fotone) e azzurra (elettrone) con la scia; dopo l'urto la linea rosa della coppia.
* **Pannello** (☰ o Tab): energia (luce visibile 2,48 eV, Mo Kα 17,48 keV, 100 keV, Cs-137 662 keV, Co-60 1,17 e
  1,33 MeV, o cursore logaritmico), polarizzazione, vista, velocità; numeri del colpo (θ, λ', E', T, angolo
  dell'elettrone contro la formula, residuo di impulso); "lancia 20 000 colpi" = istogrammi contro Klein–Nishina.
  Spazio = pausa, N = nuovo colpo.
* "Mostra solo i colpi in cui il fotone devia" (acceso di base): i surfisti nascono con la stessa |ψ|² ma ristretta
  al disco di diffusione, altrimenti bisognerebbe aspettare ~10¹¹ colpi per vederne uno. Spento: i colpi sono veri
  (il fotone passa quasi sempre).

**Non in scala, detto chiaramente.** 1 unità della scena = λ del fotone entrante. Il pacchetto-elettrone è disegnato
largo 1,6 λ (quello vero, 1 nm, è 14 λ ai raggi X e 5·10⁻⁴ λ per la luce visibile). Il guscio diffuso e la nuvola di
rinculo sono disegnati circa 10¹⁰–10¹⁸ volte più luminosi del loro peso vero (p per fotone, nel pannello).
**Rallentamento:** c = 10 λ al secondo sullo schermo, quindi il fattore è c/(10 λ):

| energia | λ | λ' a 90° | λ' a 180° | Δλ/λ max | σ/σ_T | p per fotone | bordo Compton | rallentamento |
|---|---|---|---|---|---|---|---|---|
| luce 2,48 eV | 499,94 nm | +2,43 pm | +4,85 pm | 9,7·10⁻⁶ | 1,0000 | 4,2·10⁻¹⁹ | 2,4·10⁻⁵ eV | 6,0·10¹³ |
| Mo Kα 17,48 keV | 70,929 pm | 73,356 pm | 75,782 pm | 0,068 | 0,9372 | 6,6·10⁻¹² | 1,119 keV | 4,2·10¹⁷ |
| 100 keV | 12,398 pm | 14,825 pm | 17,251 pm | 0,39 | 0,7407 | 7,7·10⁻¹² | 28,13 keV | 2,4·10¹⁸ |
| Cs-137 662 keV | 1,8729 pm | 4,2992 pm | 6,7255 pm | 2,59 | 0,3850 | 4,1·10⁻¹² | 477,65 keV | 1,6·10¹⁹ |
| Co-60 1,1732 MeV | 1,0568 pm | 3,4831 pm | 5,9094 pm | 4,59 | 0,2932 | 3,1·10⁻¹² | 963,39 keV | 2,8·10¹⁹ |
| Co-60 1,3325 MeV | 0,9305 pm | 3,3568 pm | 5,7831 pm | 5,22 | 0,2746 | 2,9·10⁻¹² | 1118,11 keV | 3,2·10¹⁹ |

Con la luce visibile l'allungamento è ~10⁻⁶ λ (limite di Thomson): sullo schermo non si vede nessuna differenza di
bande o colore se non attraverso il colore per angolo (che indica (1 − cos θ), non la grandezza assoluta: il valore vero
è nel pannello). Il pacchetto-elettrone si allarga in τ = 2mₑs²/ħ = 1,7·10⁻¹⁴ s; la scena intera dura 3·10⁻¹⁷ s ai
raggi X e meno ai gamma, quindi l'allargamento è trascurato (per la luce visibile la scena dura 2·10⁻¹³ s: lì non lo
è, ma lì l'elettrone non si muove comunque).

## 5. Validazione

Comandi (ognuno < 10 minuti, seed fissi):

    cargo test --bin compton                      # 9 test, 0,2 s
    target/debug/compton.exe --misura             # tabella dei numeri, 6 s
    timeout 590 python lab/compton_check.py       # controllo indipendente (KN riscritta con scipy), 10 s

| # | test | risultato | soglia |
|---|---|---|---|
| 1 | λ' dalle quantità di moto contro λ + λ_C(1 − cos θ), 8 energie × 721 angoli × 4 azimut | max err. rel. **2,1·10⁻¹²** | < 10⁻⁹ |
| 1 | stesso, evento per evento sui 60 000 eventi dei file CSV (Python) | max **1,8·10⁻¹²** | < 10⁻⁹ |
| 1 | angolo dell'elettrone contro tan θₑ = cot(θ/2)/(1+k) | max **6·10⁻¹⁴ rad** | < 10⁻⁹ |
| – | σ totale: tabella numerica contro formula chiusa di Klein–Nishina (6 energie) | ≤ **2,2·10⁻¹²** | < 10⁻⁹ |
| 2 | distribuzione angolare, 20 000 eventi, regola dei surfisti, 40 classi in cos θ (Rust) | χ² = **50,1 / 46,8 / 28,0** su 39 (17,48 / 662 / 1332,5 keV) | < 39 + 4√78 = 74,3 |
| 2 | stesso, controllo Python con KN indipendente (altri seed) | χ² = **33,0 / 37,2 / 32,6** su 39, p = 0,74 / 0,55 / 0,76 | p > 10⁻³ |
| 2 | χ² su 8 seed da 10 000 eventi (deve oscillare attorno a 39) | media **41,8 ± 1,6** (17,48 keV), **40,0 ± 2,1** (662 keV), KS p = 0,13 / 0,37 | KS p > 10⁻³ |
| 2 | azimut con polarizzazione lineare, 80° < θ < 100°, 12 classi | χ² = **11,0 / 11**; rapporto conteggi φ≈0 / φ≈90° = **0,107** | p > 10⁻³; < 0,25 |
| 2 | frazione deviata senza condizionamento (s = 20 fm, solo per avere p contabile) | **0,02439** contro p = **0,02481** (200 000 colpi, 1,2 σ) | < 4 σ |
| 3 | bordo Compton a 662 keV (formula) | **477,650 keV** | 477,65 |
| 3 | elettrone più energetico su 20 000 eventi a 662 keV | **477,629–477,647 keV**, nessuno sopra il bordo | ≤ bordo, > 0,99 bordo |
| 3 | spettro di T dell'elettrone contro KN, 50 classi | χ² = **53,2/47** (Rust), **44,8/47** (Python) a 662 keV; 69,0/47 a 17,48 keV (Rust) | < 47 + 4√94 = 85,8 |
| 4 | Compton 1923, Mo Kα a 90° | misurato **2,2 pm**, nostro **2,4263 pm** (−9 %) | vedi §6 |
| 5 | impulso conservato evento per evento: max \|p'_γ + p_e − p\|/p | **2,0·10⁻¹⁶** (tutti gli eventi, 3 energie) | < 10⁻⁹ |
| 5 | energia conservata evento per evento | **3,1·10⁻¹⁶** | < 10⁻⁹ |

Esito: **9/9** test Rust passati, **19/19** controlli Python passati. Le soglie dei χ² sono state fissate prima di
guardare i risultati. Il test 1 e il test 5 sono "esatti per costruzione" (l'elettrone prende l'impulso mancante):
controllano che il codice della cinematica non abbia errori, non dicono niente di nuovo sulla fisica.
Il test 2 è quello vero della regola: prova che la mappa "rango nel flusso → angolo" dei surfisti riproduce |ampiezza|²
(la statistica di Klein–Nishina). Un errore trovato così: la prima versione dell'inversione dell'azimut usciva dal ciclo
un passo troppo tardi quando trovava lo zero esatto, e la polarizzazione dava χ² = 121,8/11; corretto, 11,0/11.

Figura: `lab/out/compton_validation.png` (istogrammi contro Klein–Nishina, spettro con il bordo, Δλ(θ) con il punto
di Compton). Immagini del programma: `lab/out/compton_prima.png`, `compton_durante.png`, `compton_dopo.png`,
`compton_pannello.png`, `compton_3d.png`.

## 6. Il confronto storico (Compton 1923)

A. H. Compton, *A quantum theory of the scattering of X-rays by light elements*, Phys. Rev. 21, 483 (1923): raggi X
Mo Kα (λ = 0,708 Å) diffusi dalla grafite a 90°; la riga spostata è a 0,730 Å, cioè Δλ = **0,022 Å = 2,2 pm**, contro
**2,426 pm** dalla formula (allora 0,0242 Å). Differenza −9 %, circa 1 volta l'incertezza dello spettrometro
dell'epoca (~0,002 Å). **Attenzione:** questi valori li riportiamo come vengono citati dall'articolo; in questa sessione
non siamo riusciti a ricontrollarli sul PDF originale (accesso negato). Le misure a 45° e 135° del secondo articolo
(Phys. Rev. 22, 409, 1923) non sono in tabella per lo stesso motivo; i nostri valori sono λ' = 71,641 / 73,356 /
75,072 pm a 45° / 90° / 135° per λ = 70,93 pm.

Il bordo Compton a 477,65 keV per i gamma di 662 keV del Cs-137 è quello che si vede in ogni spettro misurato con uno
scintillatore (il valore preciso con la riga vera di 661,657 keV è 477,33 keV).

## 7. Che cosa non funziona, o funziona solo con condizioni (onestà)

* **Klein–Nishina è un ingresso.** Le nostre palline non sanno calcolare quanta onda va in ogni direzione: non abbiamo
  un modello dell'accoppiamento elettrone–campo fatto di palline. Senza questo ingresso l'esperimento non si fa.
* **La mappa dei surfisti è una scelta.** Riproduce le frequenze (verificato), ma non è la traiettoria di Bohm
  dell'onda completa a 6 dimensioni (non l'abbiamo calcolata). In particolare "b piccolo → all'indietro" è
  un'ipotesi (ragionevole, come nell'urto classico), non un risultato.
* **Il surfista-fotone prima dell'urto vola dritto** (limite del pacchetto libero largo); la regola 3-D completa di
  `base3d.rs` non è usata qui: con un pacchetto largo 10 λ e senza oggetti dà la stessa linea retta.
* **Elettrone libero e fermo.** Nella grafite di Compton gli elettroni sono legati e si muovono: questo allarga la riga
  (profilo Compton) e dà anche la riga "non spostata" (diffusione sull'atomo intero). Qui non c'è nessuno dei due.
* **Niente spin dell'elettrone né polarizzazione del fotone uscente**: Klein–Nishina è sommata su questi; il programma
  non li segue.
* **Niente correzioni radiative** (ordine α della QED oltre Klein–Nishina: ~1 % ad alta energia) né diffusione doppia.
* **Il disegno non è in scala** (§4) e il guscio diffuso è enormemente amplificato in luminosità.
* Il χ² dello spettro dell'elettrone a 17,48 keV (69,0/47, p ≈ 0,02) è un po' alto ma sotto la soglia fissata; con altri
  seed (Python) i χ² oscillano attorno al valore atteso.

## 8. Come rieseguire

    cargo build --bin compton
    cargo test --bin compton
    target/debug/compton.exe --misura
    timeout 590 python lab/compton_check.py
    target/debug/compton.exe --eventi 20000 --kev 662 --seed 1 --pol u --out lab/out/eventi.csv
    cargo run --bin compton -- --kev 17.48 --pol x --vista meta   # opzioni: --seed, --ferma-a T, --pannello,
                                                                   #   --statistica, --velocita V, --vista sezione|meta|3d
