# Piano del progetto e stato dei lavori

Aggiornato: 12 ottobre 2026, sera. Questo file serve a **riprendere il lavoro da qualunque
macchina** (anche dal cloud) senza la conversazione originale. Leggerlo per primo.

---

## 1. Lo scopo

Un simulatore in cui **tutto emerge da palline con una regola semplice**: luce → elettrone →
atomo → molecola → materiali → oggetti (la mela che marcisce) → pianeti → universo. Destinazione
finale: un motore tipo Unreal per simulazioni vere (vedi `VISIONE.md`, `VISIONE_SISTEMA.md`),
dove le impostazioni grafiche sono i livelli di agglomerato e il motore scala con l'hardware.

## 2. Le regole del lavoro (vincolanti)

1. **Niente fake.** Ingressi ammessi: costanti fisiche, numero di protoni e masse dei nuclei
   (e costanti di natura come la massa del pione). Niente tabelle di proprietà dei materiali,
   niente codici di chimica esterni. Gli esperimenti servono solo come **esame**.
2. **Tutto deve emergere.** Ogni cosa messa a mano è un **debito**: va scritta nel documento
   dello studio e sostituita col tempo.
3. **Pre-registrazione.** Prima di ogni corsa di produzione si committa un JSON con soglie e
   configurazioni (`lab/out/<studio>_prereg.json`). I risultati si giudicano contro quello.
   Cambi dopo la pre-registrazione: dichiarati.
4. **Manopole contate.** Ogni risultato riporta quanti parametri liberi usa. Un fit non è
   una previsione.
5. **Fallimenti dichiarati** come i successi. Mai gonfiare.
6. **Una sola strada: la regola zig-zag** (decisione del 12/10/2026). "Corri a c, poi ruota":
   luce se non rimbalza, elettrone se rimbalza (massa = angolo di flip θ = m c a/ħ). È l'unica
   regola in cui luce e materia sono la stessa pallina e in cui la luce diventa materia
   (coppie dal vuoto, Schwinger allo 0,3%). Le altre regole (fotone a due ruoli, DMC) sono
   storia: si tengono nel repo ma non si sviluppano più.
7. **Si sale di scala solo con l'agglomeratore**: blocca, impara la regola dei blocchi dai
   dati, controlla, ripeti; scelta automatica del raggruppamento che perde meno errore. Le
   palline non cambiano mai. **Le molecole arrivano solo dopo** che l'agglomeratore ha
   imparato la regola a molti elettroni.
8. **Tutto in HD**, anche per i calcoli: la risoluzione è fisica (la luce di traverso sul
   reticolo va con il fattore ka/tan ka; servono ≥ 24 celle per λ per errori di qualche %).
9. **Grafica**: niente linee/fili. Palline morbide che si muovono, nello stile del
   Laboratorio Fotone (palline d'onda che scorrono + surfisti luminosi con scia). Ogni cosa
   sullo schermo deve muoversi.

## 3. Dove siamo (gradino per gradino)

| Gradino | Stato | Documento |
|---|---|---|
| Fotone (regola a due ruoli, 96 esperimenti) | fatto, storico | `ESPERIMENTI_FOTONE.md`, `REGOLA_BASE.md`, `MULTIFOTONE.md` |
| Regola zig-zag: luce, elettrone, idrogeno 13,60575 eV, struttura fine | fatto | `ZIGZAG.md` |
| Fotone libero 3D zig-zag nel Laboratorio Atomo (E, B, elica, Poynting) | fatto, grafica rifatta | `ATOMO.md` |
| Fotone a palline senza spirale (13/10): palline di energia portate da Poynting, spin letto dalle palline 1,0012 come atteso; nel Laboratorio via l'elica | fatto (validazione); debito: vista da dietro nebbia bianca | `FOTONE_PALLINE.md`, `MODELLO_STANDARD_PALLINE.md` |
| Motore grafico a palline (13/10): istanze da 32 byte, regola della luce sulla GPU (= CPU bit per bit), preset bassa/media/alta/ultra; Intel UHD a bassa 60 fps e 32 ms/tick | fatto; debiti: regola GPU non collegata, Ultra non provato, livelli di disegno solo fino a λ/6 | `MOTORE_GRAFICO.md` |
| Motore grafico GPU (13/10): alta e ultra con regola, palline, misure e istanze sulla GPU (niente campo riletto); G1-G4 PASS (GPU = CPU bit per bit su llvmpipe, foto PSNR 104 dB), G5 stessi fenomeni a bassa CPU e alta GPU: flusso, segno dello spin, elicità PASS; spin dalle palline 0,91 ad alta (FAIL 5 %) e diffrazione a bassa −12,6 % (FAIL 10 %), uguali sui due cammini | fatto; debiti: spin 0,91 attorno al fuoco, rinascita diversa, Ultra non provato | `MOTORE_GRAFICO_GPU.md` |
| Atomo vivo zig-zag (ione a 1 elettrone Nd⁵⁹⁺, risonanza 1s→2p dalla regola); orbitali n ≤ 4 sulla mesh a cipolla dello zoom (13/10): tutti entro 0,36 % di Dirac (1s +0,05 %, prima +39,9 %) | fatto; l'atomo 3-D dal vivo resta con celle cubiche da 0,5 (1s +40 %, debito 12 di `NUCLEO_1S.md`) | `ATOMO.md`, `NUCLEO_1S.md` |
| Atomo vivo su un colore del cubo (13/10): il cubo = 4 reticoli indipendenti; un colore è esatto come algoritmo (4·10⁻¹⁵), regola ×4,09, passo con luce ×1,99 | fatto come opzione, **spento di default**: un colore da solo è anisotropo (E2 FAIL, dipolo 0,44 %); proposta due colori + rotazione di 180° (×2) | `SOTTORETICOLO.md` |
| Atomo vivo: accoppiamento luce–atomo centrato + trasverso (14/10): attrito finto tolto (ritmo 1,22 → 1,004 della formula, K = 1 1,006), nucleo imparato + 1s in valle k = 0: palline che girano (v_φ 0,88 di Dirac), scacchiera sparita (C4 0,90 → 10⁻⁵); export rigenerato | fatto, default; FAIL dichiarati: L1 (−0,45 %), e con il nucleo imparato p_exc/ritmo/righe non puliti (1s deriva 0,035 in 4000 tick) | `ATOMO.md` §3.7, `EMISSIONE.md` |
| Idrogeno vero a scala di Bohr (13/10): regola imparata a0/16 lontano + zig-zag a0/512 sul protone, 3-D, nel tempo, Rust su CPU (3,8 s/passo); 1s −0,066 %, 2p +0,039 %, 2s −0,138 %, Lyman-α 122,15 nm (−0,48 %), 8/8 esami PASS | fatto (validazione); debiti: ibrido non convergente nel blocco (5·10⁻⁴), colla non a flusso conservativo (10⁻⁷/passo), n = 2 solo a a0/8 | `IDROGENO_BOHR.md` |
| Vuoto di coppie (Schwinger 0,3%, Sauter 0,5%, scatola stretta) | fatto | `IDEE_RIPRESE.md`, `COPPIE_IN_SCATOLA.md` |
| Molecole dalla regola DMC: H, He, H₂ (vibrazione ora −0,2%) | fatto, storico | `MOLECOLE.md`, `H2_VIBRAZIONE.md` |
| Relatività da un fluido di palline (Lorentz 1,5%) | fatto | `FLUIDO_IN_NEGATIVO.md` |
| Protone come nodo (Skyrme): una massa, Lorentz fino a 0,6c, annichilazione 98–99% | fatto | `NODO_3D.md` |
| Pressione dentro il protone (cambio di segno 0,605 vs 0,586 fm) | fatto | `PRESSIONE.md` |
| Nuclei da nodi: toro, tetraedro, cubo nascono da soli; legame ×11–32 troppo forte | fatto | `NUCLEI_DA_NODI.md` |
| Colla dei nuclei: massa del pione, trottola quantistica, punto zero (0 parametri) non bastano (deuterone 78–160 MeV); GHS con 1 parametro: d 3.2, ³He 8.75, ⁴He 16.6 MeV a α = 0.95, ma nessun α mette a posto tutti e tre e le forme platoniche si perdono | fatto (13/10), debito: deuterone quantistico con vibrazioni accoppiate | `COLLA.md` |
| Zoom ad agglomerati sul nodo: p(0) −0,08%, E −0,16% con meno celle | fatto | `ZOOM_NODO.md` |
| Glueball (materia da luce autointeragente, 2+1 D SU(2)): massa base entro 4% | fatto, dettagli falliti | `GLUEBALL.md` |
| Carica come giro su cerchio nascosto: carica sì, massa no (×10²¹) | fatto, fallito sulla massa | `CARICA_NASCOSTA.md` |
| Forme nascoste, forma dalle palline, orbifold, occhio e flusso | fatto, emersione fallita (isole) | `FORMA_*.md`, `ORBIFOLD.md`, `OCCHIO_E_FLUSSO.md`, `GEOMETRIA_AUTOCOSTRUITA.md` |
| Equilibrio quantistico (Born come equilibrio da vortici) | continuo sì, reticolo fragile | `EQUILIBRIO_QUANTISTICO.md`, `IDEE_RIPRESE.md` |
| Sincronia / Kibble–Zurek 3D (esponente elio) | fallito: 0,507 vs 0,573 | `SINCRONIA.md` |
| Immagine complessiva dell'universo | scritta | `universo.html` (pagina), `RELAZIONE.html` |

## 4. Cantieri aperti al 12/10/2026 sera

Quattro studi erano in corso quando è stato scritto questo file. Se uno risulta interrotto
(file non committati o documento mancante), si riprende rileggendo la sua pre-registrazione.

| Cantiere | File | Obiettivo | Esami principali |
|---|---|---|---|
| **Agglomeratore sulla zig-zag** | `lab/agg_zz_*`, `docs/AGGLOMERATO_ZIGZAG.md`, `assets/atomo/agglomerato_zz.json` | imparare dai dati la regola dei blocchi della zig-zag (deve uscire Schrödinger + correzioni (v/c)²), scelta automatica dei livelli con budget di errore; idrogeno vero a scala di Bohr con poche celle; zoom fine solo vicino al protone | 1s entro 0,1% di −13,6057 eV; massa dal flip entro 1%; −k⁴/8m³ entro 10%; costo ridotto ≥ b²; legge dei costi per N elettroni e piano per imparare la regola a molti elettroni da 2–3 elettroni esatti |
| **Molti elettroni (elio, poi Li/Be)** | `lab/multi_e_*`, `docs/MOLTI_ELETTRONI.md`, `lab/out/multi_e_prereg.json` | DMC e ECG sull'operatore della regola agglomerata (m_kin e D dal livello 1); capire l'elio fino in fondo prima di Li/Be | 13/10: G0 R1 R2 R3 R4b G3 D7 G1' M1 M2 M3 E1 E2 D1 Q1 PASS; R4a G1 G5 D6 R5 R6 FAIL con cause trovate; ECG elio −2,903761 ± 0,000080, IE 24,5885 eV (E1 PASS); v3 (livello 1 T/M di AGGLOMERATO_COMPLETO): m0 = 1 a 1e-9, M1 +0,003 mEh, banda = Dirac cinetico, alto k FAIL dichiarato; elio DMC −2,903757 ± 0,000071 Eh (esatto −2,903724), ionizzazione 24,5885 ± 0,0022 eV (esp. 24,5874), 2³S −2,175328 ± 0,000072 (esatto −2,175229), solo riportati; prossimo: decisioni su G1 e D6, poi ECG elio, poi Li |
| ~~Colla dei nuclei~~ **chiuso 13/10** | `lab/colla_*`, `docs/COLLA.md` | trovare i pezzi mancanti del legame troppo forte | chiuso: vedi §3 e `COLLA.md` (C3a-iv FAIL: a α\* = 0.894 deuterone ×2.5–3.4; resta il debito del deuterone come nuvola quantistica, idea "valle del deuterone" scritta nei debiti, non corsa; sestico B = 3 non fatto) |
| **Reticolo universo** | `lab/reticolo_*`, `docs/RETICOLO_UNIVERSO.md`, `assets/atomo/reticolo_spec.json` | reticolo sparso a livelli (h_k = 2^k h_0), stessa regola della luce su ogni livello, cucitura tra livelli | velocità c su ogni livello (0,5%); passaggio fine↔grosso ≥ 99% con ≥ 16 celle/λ (curva di riflessione onesta); gamma al livello 0 e radio al livello 20 nello stesso mondo; stima onesta per una scatola grande come la Via Lattea |
| **Laboratorio Atomo: grafica + orbitali** | `src/atomo/**` (in particolare `fotscena.rs`, `render_zz.rs`, `orbitali.rs`), `docs/ATOMO.md` | (fatto) fotone con palline che cavalcano il campo e surfisti; (in corso) atomo a un elettrone con selettore n = 1..4, l = s..f, m e sovrapposizioni; da lontano forme degli orbitali come gusci morbidi, da vicino dinamica pallina per pallina | grafica almeno bella quanto il Laboratorio Fotone vecchio; niente orbitali da tabella |

## 5. Prossimi passi, in ordine

1. Chiudere i quattro cantieri sopra.
2. ~~Portare nel Laboratorio Atomo l'idrogeno **vero** a scala di Bohr~~ **motore fatto 13/10**
   (`IDROGENO_BOHR.md`, `src/atomo/idrogeno.rs`, dati per il video in `ripresa/dati/idrogeno_bohr/`).
   Fase 2 (13/10): reflux PASS (carica 10⁻¹⁶/passo), parità esatta con la regola T (AGGLOMERATO_COMPLETO §11),
   convergenza nel blocco FAIL (picco al blocco 7; blocchi 9–13: 1s −0,043 ± 0,004 %). Resta: picco del blocco 7,
   partenza dall'autostato dell'ibrido, regola della fase 4 (massa esatta), 2p a a0/16, grafica del Laboratorio.
3. **Elio e litio esatti** sulla zig-zag (o sulla regola agglomerata), poi imparare da loro la
   regola effettiva a molti elettroni. Esame: energie di ionizzazione e righe.
4. Atomi uno a uno, numero atomico per numero atomico. Dopo: legami, H₂, acqua.
5. Parete di atomi agglomerati → doppia fenditura vera → fotone che la attraversa (servono
   ≥ 24 celle per λ e un riferimento di continuo, vedi `ATOMO.md` §5).
6. Scena gamma → elettrone + positrone vicino a un nucleo pesante (vuoto di coppie in 2D HD).

## 6. Idee in lista (non dimenticare)

- **Geone** (luce tenuta dalla propria gravità): serve la gravità nel motore.
- **Sincronia / metronomi** (le particelle come difetti di una sincronia non finita): rifare
  con campo dotato di ampiezza (Ginzburg–Landau) per avere vortici larghi e mobili, conteggio
  senza esponente scelto a priori.
- **Glueball**: sistemare lo stato eccitato (base di operatori) e l'errore sistematico a β=12.
- **Il problema delle masse** (muone 207, tau 3477): nessuna forma semplice lo dà; servono
  posti di "incastro" che nascano da soli.
- **Vuoto che non preme l'universo** (10⁴⁴–10¹²⁰): collegato alla pressione del protone.
- **Anomalie reali** su cui puntare per una scoperta: vita del neutrone (878 vs 888 s), Hubble.
- **La mela**: ossidazione del ferro come primo oggetto che "invecchia" per motivi veri.

## 7. Come continuare dal cloud (o da un altro PC)

```bash
git clone https://github.com/valeriozanzarelli-debug/particle-sim.git
cd particle-sim
pip install numpy scipy torch numba ripser gudhi matplotlib
cargo build --bin atomo            # profilo dev con opt-level 3; la release LTO richiede ~35 min
cargo test  --bin atomo
```

- I programmi grafici (Bevy) vanno lanciati con la cartella del repo come working dir (gli
  shader sono in `assets/`). Su un server senza schermo non si aprono: lì si lavora sugli
  studi Python in `lab/` e sui test Rust.
- Su Linux Bevy richiede i pacchetti di sistema per audio e input (es. `libasound2-dev`,
  `libudev-dev`, `pkg-config`).
- `lab/out/` è ignorato da git: le figure e i JSON dei risultati sono aggiunti con
  `git add -f`; i file grandi (`.npz`, `.pt`, dati grezzi delle glueball) **non** sono nel
  repo e vanno rigenerati.
- Gli eseguibili (`Laboratorio Atomo.exe`, `Laboratorio Fotone.exe`) stanno fuori dal repo,
  nella cartella superiore: si ricreano con `cargo build` e si copiano.
- Commit: `git add` con percorsi espliciti, `git pull --rebase --autostash`, mai force-push.

## 8. Come lavora Valerio (per chi riprende)

- Parla italiano informale; vuole spiegazioni semplici ma oneste ("come a un ventenne").
- Lavoro lungo delegato a subagenti (modello Opus) in background, per non perdere il filo.
- Quando un programma funziona, aprirglielo subito perché veda "che forme prendono".
- Avvisarlo quando esce qualcosa di bello **e** quando non esce niente, sempre con i numeri.
- Propone spesso intuizioni fisiche: tradurle in regole e test, dire con chiarezza cosa
  corrisponde a fisica nota, cosa è nuovo, cosa è sbagliato.
- Chiavi di servizi esterni (es. Meshy) stanno in altri repository: non copiarle qui.
