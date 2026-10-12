# Laboratorio Atomo: luce e materia con una sola regola di palline

Binario: `atomo` (`src/atomo/`). Eseguibile: `Laboratorio Atomo.exe` sul Desktop (nella cartella
"Programma di simulazione"). Numeri: `assets/atomo/risultati.txt` (copiati anche dentro il programma),
log della corsa lunga `lab/out/atomo_misura.log`, script `lab/atomo_misura.sh`.

## 0. In una frase

Si apre su un fotone fatto di palline che passa una **doppia fenditura**; poi si può passare a un **atomo con un
elettrone** dove le palline si muovono tutte, a ogni tick, e un fotone le fa saltare in uno stato misto che
oscilla e rimanda luce; zoomando nel nucleo c'è il **protone come nodo** che gira e respira. Il carbonio con la
regola DMC resta come terza scena (la "nuvola media").

## 1. I tre motori, detti chiari

| motore | che cosa fa | dove si vede qui | dinamica vera? |
|---|---|---|---|
| **regola zig-zag** (`src/zigzag`, `zz.rs`, `luce.rs`) | ogni pallina avanza a c e gira il suo stato; senza massa è la luce (spin 1, F = E + icB), con la massa (giro L ↔ R di θ = a) è l'elettrone di Dirac; la fase di Coulomb lega l'elettrone al nucleo | scena 1 (doppia fenditura) e scena 2 (atomo) | **sì**: tempo reale, tick dopo tick, e all'indietro esatto |
| regola delle palline-elettrone (DMC, `carbonio.rs`) | palline che tremano, si moltiplicano o muoiono: trovano lo stato fondamentale (tempo immaginario) | scena 3 (carbonio) | **no**: dà forme ed energie di stati fermi; il tremolio è il metodo, non un moto |
| regola del fotone del laboratorio (`src/photon`) | grani che ruotano e traslano | non usata qui (la luce qui è la zig-zag) | — |

La regola zig-zag è **l'unico motore in cui luce e materia sono la stessa regola di palline** (la luce è la
stessa regola senza massa) e in cui la luce che diventa materia funziona (creazione di coppie, `docs/ZIGZAG.md`
e `docs/COPPIE_IN_SCATOLA.md`). Per questo è la strada scelta da Valerio: il resto è "eredità".

## 2. Le scene

### 2.1 Doppia fenditura (la prima prova del laboratorio nuovo)

* Palline a spin 1 su una scatola 176 × 144 × 2 celle (periodica e sottile in z: le fenditure sono lunghe, il
  problema è piano, come quello dei libri). λ = 12 celle (la luce verde del laboratorio del fotone, 500 nm:
  1 cella = 41,7 nm). Un muro di 2 celle con due fenditure (distanza 36 celle, larghe 6), una spugna sui bordi
  che assorbe, uno schermo a 106 celle dal muro che somma |E|².
* Tutto dal vivo: le creste (palline oro) avanzano, passano le fenditure, si aprono e si incrociano; sullo
  schermo i "clic" (palline azzurre) si accumulano dove arriva più energia. Un fotone dopo l'altro (stesso
  pacchetto), lo schermo continua a contare.
* **Verifica**: prima frangia luminosa contro la geometria esatta (le due strade dalle fenditure differiscono
  di una λ). Misurato (test `fenditura_frange`): **33,5 celle contro 38,3 (−12,7 %)**. Con 8 celle per λ
  −25 %, con 12 −13 %, con 16 −11 %: si avvicina piano. Il motivo: il reticolo cubico non è isotropo, in
  diagonale le palline della luce vanno più piano e le frange si stringono. **Le frange ci sono, la loro
  posizione non è ancora giusta entro 5 %**: è un debito del reticolo, non della fisica.
* Indietro nel tempo: muro e spugna assorbono, quindi la regola non si può "disfare"; si torna indietro con
  fotografie ogni 8 tick + rigioco (esatto).

### 2.2 L'atomo zig-zag (dal vivo)

* **Perché non l'idrogeno.** Il raggio dell'idrogeno è 137 ħ/mc e le celle devono essere ≲ 0,5 ħ/mc: servirebbero
  ~10¹⁰ celle. Il problema di Dirac–Coulomb dipende solo da Zα: prendiamo lo **ione idrogenoide di neodimio,
  Nd⁵⁹⁺ (Z = 60, un solo elettrone)**, raggio di Bohr 2,28 ħ/mc ≈ 880 fm. È l'atomo più grande (in Zα) che questo
  PC tiene dal vivo: 48³ celle da 0,5 ħ/mc, ~35 ms per tick, su un thread a parte (la finestra resta a 60 fps).
* **Partenza**: l'1s della forma radiale della stessa regola (stesso a, stesso Zα) messo sul reticolo cubico.
  Sulle celle grosse è solo ~60 % l'1s della regola 3-D: il resto si vede come tremolio in più (anche
  Zitterbewegung). Legame 1s della regola radiale con a = 0,5: 0,141 mc² (Dirac esatto 0,101: +40 %, celle
  grosse; con celle più fini torna, `docs/ZIGZAG.md` §5.6).
* **Ogni pallina si muove a ogni tick**: 20 000 surfisti trasportati dalla corrente dell'onda (v = j/ρ ≤ c); la
  1s di Dirac ha una corrente che gira attorno allo spin, quindi girano. Il **colore è la fase** dell'onda nel
  punto (tolta la rotazione della massa a riposo e^{−imc²t/ħ}): uno stato fermo ha densità ferma ma il suo
  colore gira alla frequenza di Bohr. La luminosità di ogni pallina segue la densità nel suo punto. Sei
  surfisti hanno una scia di palline luminose (niente linee, niente frecce).
* **Il fotone**: un impulso con l'energia scelta (keV). All'atomo il suo campo entra nella fase come e E·r
  (forma di dipolo: giusta, perché λ ≈ 55 ħ/mc è 25 volte l'atomo). Il pacchetto che arriva è disegnato (6 volte
  più corto). **La risonanza esce dalla regola**: con l'energia 1s → 2p della regola (58,5 keV, ħω = 0,1146 mc²)
  la popolazione 1s scende e la nuvola resta in una **sovrapposizione 1s + 2p** che va su e giù lungo z (il
  dipolo ⟨z⟩ misurato sull'onda). Fuori risonanza (×1,3) lo stesso impulso muove il dipolo 4 volte meno
  (0,17 contro 0,69, prova di sviluppo; test `zigzag_risonanza`, lungo, `--ignored`).
* **La luce che esce**: gusci di palline (oro = dipolo in su, azzurro = in giù) con la direzione sin²θ del dipolo
  e la luminosità del dipolo **dell'onda stessa** al tempo ritardato t − r/c.
* **Indietro nel tempo esatto**: ogni operazione della regola è unitaria, quindi il tick si disfa: 60 tick avanti
  con l'impulso e 60 indietro ridanno lo stato di partenza a |Δψ|² = 1,5·10⁻²⁷ (test `zigzag_indietro_esatto`).
  I surfisti tornano indietro con −v (approssimato, non esatto).

### 2.3 Il nucleo: il protone come nodo che gira e respira

* Il nodo B = 1 del campo dei pioni (modello di Skyrme, `docs/NODO_3D.md`), **calcolato qui** dalla sua equazione
  (tiro con bisezione): B = 1,00000, E/12π² = 1,23150 (letteratura 1,2315), raggio rms 0,59 fm, massa 863 MeV.
* Disegnato con 40 000 palline estratte a caso dalla densità barionica (raggi e direzioni casuali: mai i punti di
  un reticolo); colore = direzione del campo dei pioni (nel riccio punta lungo r̂), luminosità = densità.
* **Gira**: per il riccio un'isorotazione è una rotazione, quindi si vede girare il colore. Λ dalla stessa
  soluzione (con la coda 1/r² fuori dalla griglia): M_Δ − M_N = 3ħ²/2Λ = **319 MeV** (misura 293, +9 %).
  **Respira**: modo di scala (Derrick), ħω = 216 MeV (stima nostra; la risonanza di Roper è ~500 MeV sopra il
  nucleone: la stima è grezza). Orologio del nucleo rallentato (1 s di schermo ≈ 4·10⁻²⁴ s), ampiezza del
  respiro scelta per vederla (7 %).

### 2.4 Carbonio (regola DMC, nuvola media)

Sei elettroni come palline nello spazio a 18 dimensioni, nodi fissi da una guida fatta dai nostri orbitali
idrogenoidi (un determinante per spin, ³P scritto come M_S = 1: 4 su + 2 giù; la richiesta diceva 3 + 3, ma
M_S = 0 avrebbe bisogno di due determinanti). Gli esponenti vengono dal minimo della nostra energia variazionale.
Fotoni: assorbiti solo alle energie che le palline hanno calcolato (2p → 3s, ionizzazione), altrimenti passano.
Temperatura: il nucleo si muove come un atomo di gas a 1 atm (Maxwell–Boltzmann, urti al tempo libero medio);
si vede solo a velocità del tempo alte (ps al secondo). È **la vista statica** che Valerio ha scartato come scena
principale: resta per i numeri.

## 3bis. Risultati della corsa di produzione (pre-registrata sotto, §3)

| prova | risultato | esito |
|---|---|---|
| V1 energia di C | **−37,82587 ± 0,00106 Eₕ** contro −37,8450 (+0,051 %) | **PASS** |
| V2 ionizzazione | **11,334 ± 0,032 eV** contro 11,2603 (+0,66 %) | **PASS** |
| V3 2p → 3s | **8,036 ± 0,041 eV** contro 7,482 (+7,4 %) | **FAIL** |
| V4 protone | B = 1,00000, E/12π² = 1,23150, rms 0,592 fm, M = 863,1 MeV | come previsto |

V3: lo stato eccitato con un solo determinante e un 3s idrogenoide con due parametri sta 0,55 eV troppo alto
(nella corsa di sviluppo, con esponenti diversi, era 7,55 eV): il nodo dello stato eccitato è il punto debole.
La scena del carbonio usa comunque il numero calcolato (λ = 154,3 nm), non quello misurato. Il guscio K/L le
palline lo trovano da sole al minimo della distribuzione radiale: 0,64 bohr (0,34 Å).

## 4. Che cosa è vero, che cosa è disegnato, che cosa è un debito

**Dinamica vera (la regola, tick dopo tick)**: la luce della doppia fenditura; l'onda dell'elettrone dell'atomo
zig-zag (e i surfisti che porta); la risposta alla luce (risonanza, sovrapposizione, dipolo); il ritorno indietro.
**Precalcolato**: il profilo del nodo (poi animato con le sue frequenze); le energie DMC del carbonio; lo stato
di partenza 1s (dalla forma radiale della stessa regola).
**Disegnato**: il pacchetto del fotone che arriva all'atomo (la sua forza nel calcolo è il campo E(t)); i gusci
della luce uscente (la loro luminosità è il dipolo calcolato); la foschia termica.

Debiti:
1. Il campo del fotone sull'atomo è un **impulso forte** (come un laser), non un fotone solo.
2. **Emissione spontanea** (vuoto quantistico) assente: la sovrapposizione non decade da sola.
3. Celle grosse: legame 1s +40 %, frange −13 %, partenza solo 60 % 1s della regola 3-D.
4. Nd⁵⁹⁺ invece dell'idrogeno (stesso Zα·fisica, Z diverso); il nucleo di Nd ha 144 barioni, si vede un protone.
5. I due numeri del modello di Skyrme (F_π = 129 MeV, e = 5,45) presi da Adkins–Nappi–Witten (1983).
6. Temperatura: nella scena zig-zag non ha effetti visibili (in 10⁻¹⁹ s un nucleo a 300 K si sposta di 10⁻⁵ fm);
   la foschia è solo visiva. Nel carbonio muove il nucleo (gas a 1 atm, urti come sfere dure).
7. Carbonio: nodi fissi; V3 FAIL; forza della riga e vita media non calcolate (attesa compressa).

## 5. Come si usa

* Doppio clic su `Laboratorio Atomo.exe`. Schermo pulito: barra del tempo in basso (⏮ indietro, ▶/⏸, ⏭,
  ritmo), scala in basso a sinistra, scene e zoom in alto a sinistra, pannello con Tab o ☰.
* Rotella = zoom logaritmico (dalla scala della luce, micrometri, fino a 0,1 fm), trascinare = girare attorno. Tasti: Spazio
  pausa, F spara il fotone, 1 vista intera, 2 nucleo, 3 elettrone.
* Riga di comando: `"Laboratorio Atomo.exe" --scena fenditura|atomo|carbonio [--zoom nucleo|elettrone]`.
  **Atomo aperto dentro il nucleo**: `"Laboratorio Atomo.exe" --scena atomo --zoom nucleo`.
* Misure: `atomo --misura tutto` (≈ 11 min) o `sh lab/atomo_misura.sh`; test: `cargo test --bin atomo`
  (7 test, ~15 s) e `cargo test --bin atomo -- --ignored` (lunghi).
* Prove senza mani: `--foto file.png --frames N [--spara-al F]` salva un'immagine ed esce (mouse e tastiera
  ignorati).

## 3. Pre-registrazione (prima della corsa di produzione)

Corse di sviluppo già viste (dichiarate): VMC + DMC brevi (1000 palline, 2000–3000 tick, τ = 0,005) hanno dato
E(C) = −37,830 ± 0,007 e −37,819 ± 0,004 Eₕ; E(C⁺) = −37,412 ± 0,004; E(C*, 2p3s ³P°) = −37,542 ± 0,004;
quindi ionizzazione ≈ 11,08 eV ed eccitazione 2p→3s ≈ 7,55 eV. Protone (nodo B = 1): E/12π² = 1,23150,
B = 1,0000, raggio rms 0,59 fm.

Corsa di produzione (fissata ora): `atomo --misura tutto` = VMC a campionamento correlato (30 000
configurazioni, 3 giri) per gli esponenti di C, C⁺, C*; poi DMC a nodi fissi con 2000 palline, τ = 0,005 Eₕ⁻¹·ħ,
4000 tick di equilibrio e 40 000 tick misurati per ciascuno stato; errore dai blocchi (il più grande tra 32, 16, 8
blocchi).

| prova | cosa | soglia |
|---|---|---|
| V1 | energia totale di C (³P) contro l'esatta non relativistica −37,8450 Eₕ (Chakravorty et al., Phys. Rev. A 47, 3649, 1993) | entro **1 %** |
| V2 | energia di ionizzazione E(C⁺) − E(C) contro 11,2603 eV (NIST ASD) | entro **3 %** |
| V3 | eccitazione 2s²2p² ³P → 2s²2p3s ³P° contro 7,482 eV (NIST ASD, media sui J) | entro **3 %** |
| V4 | il protone: B = 1 entro 10⁻³, E/12π² entro 0,2 % da 1,2315 | informativa (già vista nel dev) |

Le soglie V2 e V3 non vengono cambiate dopo la corsa. Se V3 fallisce, la scena di risonanza del carbonio usa
comunque il numero calcolato (non quello misurato) e il fallimento viene scritto qui.
