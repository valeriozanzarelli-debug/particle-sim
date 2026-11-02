# Il Modello Standard fatto di palline: cosa si può, cosa no, in che ordine

Scritto il 12/10/2026 dopo la domanda di Valerio: "riusciamo a rappresentare tutte le particelle del modello
standard, solo con le nostre palline, in HD, mostrando i flussi e dove l'energia è più densa?".

Questo è un **parere ragionato**, non un risultato: dice per ogni particella che cosa la regola zig-zag fa
già, che cosa potrebbe fare con lavoro vero, e che cosa oggi sarebbe solo un disegno (quindi vietato dalle
regole del progetto). Il primo pezzo concreto è il fotone senza spirale: `FOTONE_PALLINE.md`.

---

## 1. La regola in una riga (per ricordarsi da dove partiamo)

Ogni pallina-onda **corre a c e poi gira il suo stato interno** (`ZIGZAG.md` §2). Cambia solo:

* il **tipo** di stato interno (spin 1/2 sinistrorso o destrorso, spin 1, …);
* l'**angolo di flip** θ = m c a / ħ, cioè la **massa**;
* la **fase** che prende in un campo, cioè la **carica**.

Tutto quello che non si scrive in queste tre manopole va aggiunto alla regola, e ogni aggiunta va contata.

## 2. Particella per particella

| particella | cosa serve nella regola | stato | onesto? |
|---|---|---|---|
| **fotone** γ | palline spin 1, F = E + icB, θ = 0 | **fatto** (ZIGZAG §3, ATOMO §2, 4.1) | sì: ω = c\|k\| esatto lungo gli assi, trasversale, spin 1 (fase di Berry) |
| **elettrone, positrone** | Weyl L + R, flip θ = mₑca/ħ | **fatto** (idrogeno 13,60575 eV, struttura fine, Zitterbewegung) | sì, con la massa come ingresso |
| **neutrini** (3) | Weyl L con θ minuscolo (o zero) | fattibile subito | la massa e il mescolamento (3 angoli + 1 fase) sono ingressi; le oscillazioni escono dalla regola se si mettono 3 flip diversi |
| **muone, tau** | lo stesso elettrone con θ × 207 e × 3477 | fattibile subito | la massa **non emerge**: è il debito "problema delle masse" (PIANO §6). Il decadimento del muone no (serve la forza debole quantizzata) |
| **quark** (6) | zig-zag come l'elettrone + uno stato di **colore** a 3 componenti, carica 2/3 o −1/3 | fattibile, lavoro grosso | il moto sì; il **confinamento no** in tempo reale (vedi §3) |
| **gluoni** (8) | palline spin 1 come la luce, ma con colore: si urtano tra loro (non lineari) | fattibile come campo classico su reticolo (Yang–Mills classico) | sì come campo classico; masse delle glueball solo in tempo immaginario (GLUEBALL.md, entro 4 % in 2+1 D SU(2)) |
| **Higgs** | pallina spin 0 (due camminatori che vanno avanti e indietro, come Klein–Gordon), con il suo "cappello messicano" | fattibile | i numeri del potenziale sono costanti di natura (ammessi), ma vanno contati |
| **W±, Z** | "luce" spin 1 che si appoggia al Higgs e diventa pesante (meccanismo di Higgs) | fattibile come simulazione classica (elettrodebole su reticolo, sfaleroni) | sì; il **decadimento** in 3·10⁻²⁵ s no (serve la luce quantizzata) |
| **protone, neutrone** | oggi: **nodo di Skyrme** (NODO_3D, PRESSIONE, NUCLEI) | fatto come campo dei pioni | è un modello efficace, non i quark; il legame dei nuclei è ancora × 11–32 (COLLA) |

## 3. Le tre cose che oggi **non** si possono far emergere (e perché)

1. **Il confinamento dei quark in tempo reale.** La QCD su reticolo funziona in "tempo immaginario"
   (Monte Carlo, come le nostre glueball). In tempo reale, con molti quark, il calcolo esatto cresce in modo
   esponenziale (il "problema del segno"). Un protone fatto di tre quark che restano chiusi da soli,
   guardato mentre si muove, oggi non lo sa fare nessuno. Disegnarlo sarebbe un fake.
2. **I decadimenti** (muone → elettrone + neutrini, W → coppia, Z → coppia). Servono campi quantizzati
   (creazione e distruzione di particelle una alla volta). La nostra regola fa onde; la luce quantizzata è già
   un debito aperto (CONTESTO_CLOUD §4, voce 17). Le **coppie dal vuoto** invece le facciamo già (Schwinger allo 0,3 %).
3. **Le masse** (perché il muone pesa 207 elettroni). Nessuna regola provata qui le fa uscire; restano ingressi.

## 4. "Urti e deformazioni": dove sono veri

Nella regola zig-zag lineare (luce libera, elettrone libero) le palline **non si urtano**: si sommano
(interferenza). Le forme che vediamo lì sono forme del flusso, non di scontri. Gli urti e le deformazioni veri
arrivano solo dove c'è **non linearità**:

* i gluoni tra loro (Yang–Mills);
* il nodo del protone (Skyrme: due nodi che si scontrano si deformano, si fondono in toro, cubo, …: NUCLEI_DA_NODI);
* l'elettrone nel campo che lui stesso cambia (Maxwell–Dirac accoppiato: `accoppiato.rs`);
* le coppie che nascono dal vuoto.

Questo va detto quando si fa vedere qualcosa: un fotone che "rimbalza" su un altro fotone nel vuoto non esiste
(salvo la diffusione luce-luce, 10⁻⁶⁵ cm² a energie ottiche).

## 5. Come lo farei vedere (GPU, HD)

Lo stesso schema per tutte le particelle, dal piccolo al grande:

* **Il campo vero** gira sulla GPU come texture 3-D (la regola è locale: un passo è uno shader per asse).
  Risoluzione fisica: 24–32 celle per λ (o per ħ/mc).
* **Le palline che si vedono** sono pezzetti della quantità che si conserva: energia per la luce, probabilità
  (o carica) per l'elettrone. Nascono una volta sola con densità proporzionale a quella quantità e poi vengono
  **portate dalla sua corrente** (Poynting per la luce, corrente di Dirac per l'elettrone). Siccome la quantità
  si conserva, la nuvola resta fedele al campo senza ricrearla: si può controllare (esame F1 del fotone).
* **Forma della pallina** = un'informazione locale del campo (per la luce: allungata lungo E; per l'elettrone:
  lungo lo spin). **Colore** = una fase o una direzione. **Luminosità e bagliore** = densità di energia.
* **Trasparenza** "pesata per profondità" (weighted blended OIT): milioni di palline senza ordinarle.
* Bilancio su una GPU da 16 GB: campo 256³–384³ (×3 componenti complesse) + 2–4 milioni di palline a
  60 fotogrammi; il portatile (solo Intel UHD: la GTX 1650 Ti è rotta) controlla a qualità bassa o a fotogrammi singoli.
* Dal piccolo al grande: a ogni livello dell'agglomeratore le palline sono le **stesse quantità conservate**
  del livello sopra, sommate per blocco. Così l'immagine non cambia natura quando si fa zoom indietro.

## 6. Ordine proposto (nessuna molecola prima che atomi e quantistica siano a posto)

1. **Fotone a palline** (in corso, `FOTONE_PALLINE.md`): via la spirale disegnata, spin letto dalle palline.
2. **Elettrone a palline**: stessa idea con la corrente di Dirac (Zitterbewegung visibile come tremolio
   della nuvola); esame: ritrovare 2mc²/ħ e ħ/2mc già misurati in ZIGZAG.
3. **Muone, tau, neutrini**: la stessa scena con θ diversi (validazione, nessuna scoperta: la massa è un ingresso).
4. **Higgs e W/Z classici** su reticolo piccolo (meccanismo di Higgs visto con le palline).
5. **Gluoni classici** (SU(3) di Yang–Mills in tempo reale): qui si vedono gli urti veri tra palline di colore.
6. Quark e confinamento: solo dopo aver chiuso i debiti degli atomi, e solo con un'idea che non sia un disegno.
