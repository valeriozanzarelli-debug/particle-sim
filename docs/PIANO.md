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
| Video lenti a bagliori (13/10): atomo Z = 60 che gira (~0,83 di Dirac, errore del passo), assorbe e rimanda luce come antenna (esportazione c7c193e6bf); idrogeno vero dal grande al piccolo (dati fase 2). Palline = bagliori additivi, niente scie | fatto; FAIL dichiarati: atomo E1 (fedeltà 0,958), E5 (1,61, deriva dell'1s), E3b; idrogeno H6 (previsto) | `VIDEO_ATOMO.md`, `VIDEO_IDROGENO.md`, `modello-standard/video_atomo_bagliori_*`, `video_idrogeno_bohr_*` |
| Video supercritico a bagliori (13/10): l'1s di Z = 200 si tuffa, il guscio (positrone, lettura delle buche) esce a 0,90 c; dati radiali del film ZSTAR (fc8a08c), palline a coordinata di massa | fatto; V1/V3/V4/V5 PASS, V2 FAIL (export a 4 cifre: continuità 1e-4 per fotogramma) | `VIDEO_SUPERCRITICO.md`, `modello-standard/video_supercritico_*` |
| Motore grafico GPU (13/10): alta e ultra con regola, palline, misure e istanze sulla GPU (niente campo riletto); G1-G4 PASS (GPU = CPU bit per bit su llvmpipe, foto PSNR 104 dB), G5 stessi fenomeni a bassa CPU e alta GPU: flusso, segno dello spin, elicità PASS; spin dalle palline 0,91 ad alta (FAIL 5 %) e diffrazione a bassa −12,6 % (FAIL 10 %), uguali sui due cammini | fatto; debiti: spin 0,91 attorno al fuoco, rinascita diversa, Ultra non provato | `MOTORE_GRAFICO_GPU.md` |
| Atomo vivo zig-zag (ione a 1 elettrone Nd⁵⁹⁺, risonanza 1s→2p dalla regola); orbitali n ≤ 4 sulla mesh a cipolla dello zoom (13/10): tutti entro 0,36 % di Dirac (1s +0,05 %, prima +39,9 %) | fatto; l'atomo 3-D dal vivo resta con celle cubiche da 0,5 (1s +40 %, debito 12 di `NUCLEO_1S.md`) | `ATOMO.md`, `NUCLEO_1S.md` |
| Atomo vivo su un colore del cubo (13/10): il cubo = 4 reticoli indipendenti; un colore è esatto come algoritmo (4·10⁻¹⁵), regola ×4,09, passo con luce ×1,99 | fatto come opzione, **spento di default**: un colore da solo è anisotropo (E2 FAIL, dipolo 0,44 %); proposta due colori + rotazione di 180° (×2) | `SOTTORETICOLO.md` |
| Atomo vivo: accoppiamento luce–atomo centrato + trasverso (14/10): attrito finto tolto (ritmo 1,22 → 1,004 della formula, K = 1 1,006), nucleo imparato + 1s in valle k = 0: palline che girano (v_φ 0,88 di Dirac), scacchiera sparita (C4 0,90 → 10⁻⁵); export rigenerato | fatto, default; FAIL dichiarati: L1 (−0,45 %), e con il nucleo imparato p_exc/ritmo/righe non puliti (1s deriva 0,035 in 4000 tick) | `ATOMO.md` §3.7, `EMISSIONE.md` |
| Motore grafico di base a bagliori (13/10): atomo in 3-D (luce aperta verso la camera, profondità, rotazione lenta), dettaglio continuo (da lontano meno bagliori e più grandi, niente palla bianca; la luce passa da bagliori di cella a scintille fini zoomando), surfisti più chiari con scia di 0,3 s, finestra Impostazioni (G), grafo dei livelli + budget per fotogramma; interfaccia unica `bagliori.rs` per tutte le scene | fatto nella scena atomo; da portare: orbitali, tavola, molecole, quark/gluoni | `ATOMO.md` §3.9 |
| Idrogeno vero a scala di Bohr (13/10): regola imparata a0/16 lontano + zig-zag a0/512 sul protone, 3-D, nel tempo, Rust su CPU (3,8 s/passo); 1s −0,066 %, 2p +0,039 %, 2s −0,138 %, Lyman-α 122,15 nm (−0,48 %), 8/8 esami PASS | fatto (validazione). **Fase 4 (13/10, §13):** bug del bordo dell'alone trovato (perdita 5·10⁻⁶ → 3·10⁻⁸/passo), blocchi 7/9/11 a plateau entro 3·10⁻⁵; massa del maestro sistemata dall'agglomeratore (maestro a0/65536, m−1 = 1,5·10⁻⁶; 2p +0,0004 %); grosso di default T + bordi ritorti + Born: 1s/2s/2p entro 0,0034 % da Dirac. Debiti: corrente ferma al bordo del blocco (~1·10⁻⁴, FAIL), bolla 3·10⁻⁵ al blocco 7, Born e campo di Galerkin non ancora nel Rust | `IDROGENO_BOHR.md` |
| Vuoto di coppie (Schwinger 0,3%, Sauter 0,5%, scatola stretta) | fatto | `IDEE_RIPRESE.md`, `COPPIE_IN_SCATOLA.md` |
| Molecole dalla regola DMC: H, He, H₂ (vibrazione ora −0,2%) | fatto, storico | `MOLECOLE.md`, `H2_VIBRAZIONE.md` |
| Relatività da un fluido di palline (Lorentz 1,5%) | fatto | `FLUIDO_IN_NEGATIVO.md` |
| Protone come nodo (Skyrme): una massa, Lorentz fino a 0,6c, annichilazione 98–99% | fatto | `NODO_3D.md` |
| Pressione dentro il protone (cambio di segno 0,605 vs 0,586 fm) | fatto | `PRESSIONE.md` |
| Nuclei da nodi: toro, tetraedro, cubo nascono da soli; legame ×11–32 troppo forte | fatto | `NUCLEI_DA_NODI.md` |
| Colla dei nuclei: massa del pione, trottola quantistica, punto zero (0 parametri) non bastano (deuterone 78–160 MeV); GHS con 1 parametro: d 3.2, ³He 8.75, ⁴He 16.6 MeV a α = 0.95, ma nessun α mette a posto tutti e tre e le forme platoniche si perdono | fatto (13/10), debito: deuterone quantistico con vibrazioni accoppiate | `COLLA.md` |
| Zoom ad agglomerati sul nodo: p(0) −0,08%, E −0,16% con meno celle | fatto | `ZOOM_NODO.md` |
| Glueball (materia da luce autointeragente, 2+1 D SU(2)): massa base entro 4% | fatto, dettagli falliti | `GLUEBALL.md` |
| Quark e gluoni (13–14/10): SU(3) 3+1 D in tempo immaginario, tre quark statici, tubo a Y misurato (R_Y = 1,23 ± 0,05; placchetta e a√σ = 0,400 PASS) + Yang–Mills SU(3) classico in tempo reale (si attraversano a piccola ampiezza, si urtano a grande, D ∝ amp^2,04) | fatto; G7 FAIL (fedeltà palline 0,8998); film in montaggio | `QUARK_GLUONI.md` |
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
| **Molti elettroni (elio, poi Li/Be)** | `lab/multi_e_*`, `docs/MOLTI_ELETTRONI.md`, `lab/out/multi_e_prereg.json` | DMC e ECG sull'operatore della regola agglomerata (m_kin e D dal livello 1); capire l'elio fino in fondo prima di Li/Be | 13/10: G0 R1 R2 R3 R4b G3 D7 G1' M1 M2 M3 E1 E2 D1 Q1 PASS; R4a G1 G5 D6 R5 R6 FAIL con cause trovate; ECG elio −2,903761 ± 0,000080, IE 24,5885 eV (E1 PASS); v3/v5 (livello 1 di AGGLOMERATO_COMPLETO, fase 5: k⁶ fissato): m0 = 1, M1 +0,003 mEh, Pauli Z=10 6/6 semi (L1c_v5 PASS), banda = Dirac cinetico a 0,03 %, DEBT-4/9 sbagliato da +0,41 a +0,0002 mEh; D6' (cuspide del misto) FAIL di poco: 0,447 ± 0,029 contro 0,50 ± 0,05; elio DMC −2,903757 ± 0,000071 Eh (esatto −2,903724), ionizzazione 24,5885 ± 0,0022 eV (esp. 24,5874), 2³S −2,175328 ± 0,000072 (esatto −2,175229), solo riportati; elio chiuso e documentato; prossimo (non in questo giro): litio con l'ECG, D2 fermo per G5 |
| ~~Colla dei nuclei~~ **chiuso 13/10** | `lab/colla_*`, `docs/COLLA.md` | trovare i pezzi mancanti del legame troppo forte | chiuso: vedi §3 e `COLLA.md` (C3a-iv FAIL: a α\* = 0.894 deuterone ×2.5–3.4; resta il debito del deuterone come nuvola quantistica, idea "valle del deuterone" scritta nei debiti, non corsa; sestico B = 3 non fatto) |
| **Reticolo universo** | `lab/reticolo_*`, `docs/RETICOLO_UNIVERSO.md`, `assets/atomo/reticolo_spec.json` | reticolo sparso a livelli (h_k = 2^k h_0), stessa regola della luce su ogni livello, cucitura tra livelli | velocità c su ogni livello (0,5%); passaggio fine↔grosso ≥ 99% con ≥ 16 celle/λ (curva di riflessione onesta); gamma al livello 0 e radio al livello 20 nello stesso mondo; stima onesta per una scatola grande come la Via Lattea |
| **Laboratorio Atomo: grafica + orbitali** | `src/atomo/**` (in particolare `fotscena.rs`, `render_zz.rs`, `orbitali.rs`), `docs/ATOMO.md` | (fatto) fotone con palline che cavalcano il campo e surfisti; (in corso) atomo a un elettrone con selettore n = 1..4, l = s..f, m e sovrapposizioni; da lontano forme degli orbitali come gusci morbidi, da vicino dinamica pallina per pallina | grafica almeno bella quanto il Laboratorio Fotone vecchio; niente orbitali da tabella |
| **Bagliori** (13/10) | `docs/BAGLIORI.md`, `lab/out/gl_*`, `src/atomo/elscena.rs` | fotone leggibile di fronte (bagliori che seguono u, tavolozza lineare, ordine di profondità), `--sequenza` per clip al rallentatore; poi elettrone/muone libero a bagliori | fotone: S1 S2 S3 I1 I2 I3 I5 PASS (I3 dopo emendamento); elettrone (fermato, prossimo passo): E1 E2 PASS, E3 FAIL |
| **Modello standard a palline** (in pausa 14/10, 01:00) | `docs/VIDEO_H2.md`, `docs/VIDEO_ELIO.md`, `docs/QUARK_GLUONI.md`, `lab/h2v_*`, `lab/he_v_*`, `lab/qg_*`, `lab/h2v_zzlook.py` (aspetto zig-zag copiato da `zzscena.rs`/`render_zz.rs`/`balls.wgsl`) | H₂ ed elio a bagliori (8/9 PASS + V7 FAIL dichiarato; 5/5 PASS), tubo a Y dei tre quark finale (600/600, A B C PASS, `docs/fig/qg_tre_quark_Y.png`); clip H₂ ed elio **parziali** (612/912 e 260/432 fotogrammi: disco pieno e cartella condivisa in errore I/O), film dei gluoni: urto ad ampiezza grande **montato** (`docs/fig/qg_gluoni_urto.gif`), ampiezza piccola 54 % e tre gluoni 78 % fermi (checkpoint solo nello scratch del container: se persi si riparte da zero, comandi in `QUARK_GLUONI.md` §7); G7 FAIL 0,8998 | **alla ripresa, primo passo**: portare H₂, elio, tubo a Y e gluoni come scene dentro il programma unico sul renderer zig-zag (`bagliori.rs`, interfaccia unica del thread principale), usando le funzioni `carica()`/`disegna()` già separate; poi finire le clip (`python3 lab/h2v_render.py batch 0 912 --scala 0.6667`, `python3 lab/he_v_render.py batch 0 432 --scala 0.6667`, poi ffmpeg) e il film dei gluoni dai checkpoint; poi nuclei di nucleoni (B = 1–4 da `NUCLEI_DA_NODI.md` a bagliori, poi ⁶Li/⁷Li, ¹²C, ¹⁶O con pre-registrazione) e la modalità Modello standard (particella per particella, con stato onesto: dal vivo e⁻ e⁺ fotone muone; gluoni solo classici; quark solo statici; forza debole spenta) |

## 4bis. Pausa del 14/10/2026 (01:00 UTC): dove si è fermato il programma unico

Valerio ha messo in pausa tutto. Nessun agente attivo. Stato del **programma unico** `atomo` (Laboratorio Atomo):
- **Su master e compilato sul portatile** (ponte `2026-10-14_build_tavola`, git 278b739):
  - motore a bagliori 3-D (`src/atomo/bagliori.rs`, `grafica.rs`, ATOMO.md §3.9, commit 3937246): camera girabile, dettaglio
    continuo con lo zoom (niente palla bianca), scie dei surfisti 0,3 s, finestra Impostazioni grafica (G);
  - tavola periodica (`tavola.rs`, `esportazioni.rs`, §10, e7052b5): accesi H, He, H₂, Nd⁵⁹⁺; formato delle esportazioni
    `assets/esportazioni/<nome>/manifest.json + campioni.f32` (§10.2), pronto per Li e le molecole dei Debiti;
  - cascata dell'idrogeno (`cascata.rs`, §10.4, ed6254d + d739593): 1s…3d entro 0,03 % di Dirac, A(2p) +0,07 %, ramo 3p→2s
    11,83 %; C4 (2s a due fotoni) non calcolato.
- **Comandi sul portatile** (cartella `ponte_lab`): `target\debug\atomo.exe --scena atomo --surfisti --distanza 60`;
  `--scena tavola`; `--scena tavola --elemento He --buca`.
- **Il ponte.py che gira sul portatile è vecchio**: le foto con `scena tavola` o `spara-pronto` falliscono finché Valerio non lo
  riavvia (`2026-10-14_foto_*` FALLITE per questo, non per il programma).

**Prossimi passi, in ordine** (per chi riprende a freddo):
1. **Vista da lontano (LOD)**: a distanza 300 l'atomo è un puntino sopra i gusci di luce; le palline bianche di flusso sono
   invadenti a distanza 60. Sistemare in `render_zz.rs`/`zzscena.rs` restando sulla stessa rappresentazione (§3.9, debiti).
2. **Tavola e cascata sul motore nuovo**: oggi disegnano col vecchio percorso; portarle su `bagliori.rs` (regola di Valerio:
   la grafica dell'atomo zig-zag è la base di tutto, nessuno stile parallelo). Poi orbitali, elettrone e fotone liberi.
3. **Budget di fotogramma sull'Intel UHD**: provato solo in cloud (lavapipe). Foto/misure dal ponte dopo il riavvio di ponte.py.
4. **Clip a bagliori dentro il programma come scene**: i video (atomo che assorbe/riemette, idrogeno di Bohr, elio, H₂, quark
   del Modello standard) devono diventare scene di `atomo` sul motore comune, non file a parte.
5. Accendere nella tavola gli elementi/molecole che Debiti esporta (Li chiuso, LiH in corso).
6. Idrogeno di Bohr: campo di Galerkin nel motore Rust (il −0,05 % residuo dell'ibrido, IDROGENO_BOHR §13).

## 4ter. Pausa del 14/10/2026 (01:00 UTC): dove si è fermata la fisica dei Debiti

Nessun agente attivo. Tutto su master (commit con firma, file grandi .npy/.npz non tracciati; basi ECG salvate in JSON).
- **Chiusi**: litio esatto 11/0 (`docs/TAVOLA_LITIO_BERILLIO.md`, ionizzazione 5,39133 eV, riga 670,99 nm, vita 2p 27,107 ns);
  LiH esatto 5/1 (`docs/MOLECOLE_PICCOLE.md`: R_e −0,015 %, D_e −22 meV, ω_e +0,24 %, dipolo +0,73 %; FAIL L5 base del Li
  +0,56 mEh); H₂ 6/1 e forza H₂–H₂ 4/0 (`docs/AGGLOMERATO_DUE_ATOMI.md` §12–§15); atomo ibrido entro 1 % (AGGLOMERATO_COMPLETO
  §14–§17); 2s +0,034 % (AGGLOMERATO_2S); surfisti limite 0,5 % (SURFISTI_ESATTI §15); ciclo tra livelli fasi 1–2
  (`docs/CICLO_LIVELLI.md`: punto fisso κ², variabili grosse automatiche, φ perde contro b=3).
- **Export per la tavola**: `lab/out/tav_dati_litio/manifest.json` + `litio_livelli.json`, `lab/out/mol_dati_lih/manifest.json`;
  campioni da rigenerare dalle basi JSON (`lab/out/tav_basi/LEGGIMI.txt`, `lab/out/mol_basi/`) con
  `python lab/tav_export.py li`. Da convertire nel formato `assets/esportazioni/<nome>/` (§4bis) e accendere Li e LiH.

**Prossimi passi della fisica, in ordine** (per chi riprende a freddo):
1. **Berillio**: corse ferme a metà (1S K=700, +0,43 mEh non giudicato; ¹P K=228/450). Ricostruire i .npz dalle basi JSON,
   `lab/tav_drive.sh be1S` e `be1P` (~2–3 h), poi `python lab/tav_analyze.py be` e `python lab/tav_export.py be`.
2. **LiH**: ~3 punti in più a R = 2,1–2,4 a₀ (~30 min) per togliere il muro a 2,5 a₀ (salto 0→1 1481 contro 1360 cm⁻¹).
3. **Boro → neon** (`docs/TAVOLA_GRIGLIA.md`, `lab/tavg_*`, PARZIALE): HF fine B, C PASS; fare N–Ne, blocchi-nucleo per
   tutti (serve anche H e N per le molecole), esami di griglia G1–G7, poi fase 2 correlazione imparata dai casi esatti
   (H₂, He, Li, Be) con ionizzazioni alla cieca. Atteso FAIL B3 (promozione 2s→2p) per C, N, O.
4. **Acqua → goccia** (stella polare, vedi memoria): H₂O sulla griglia, dimero e trimero (`lab/mol_acqua_geom.py`,
   scan in `lab/out/mol_acqua_scan.json`), forza acqua–acqua imparata (niente TIP4P), goccia di ~1000 molecole
   con nuclei classici; esami pre-registrati: densità, picco O–O 2,8 Å, ~3,5 legami H, tensione superficiale.
5. **Ciclo tra livelli fase 3** (PARZIALE, `lab/ciclo_p3_*`): scrivere il codice misto con variabili lisciate (s=3),
   fissare le soglie B1–B5 e committarle prima delle corse, poi produzione e pannello f della figura. Fase 4 = test a due
   vie sulla goccia (solo piano, esami C1–C4).
6. **H₂ liquido**: scambio-dispersione, C6 entro 3 %, palline-molecola e B(T), poi path integral.

## 5. Prossimi passi, in ordine

1. Chiudere i quattro cantieri sopra.
2. ~~Portare nel Laboratorio Atomo l'idrogeno **vero** a scala di Bohr~~ **motore fatto 13/10**
   (`IDROGENO_BOHR.md`, `src/atomo/idrogeno.rs`, dati per il video in `ripresa/dati/idrogeno_bohr/`).
   Fase 2 (13/10): reflux PASS (carica 10⁻¹⁶/passo), parità esatta con la regola T (AGGLOMERATO_COMPLETO §11),
   convergenza nel blocco FAIL (picco al blocco 7; blocchi 9–13: 1s −0,043 ± 0,004 %). Resta: picco del blocco 7,
   partenza dall'autostato dell'ibrido, regola della fase 4 (massa esatta), 2p a a0/16, grafica del Laboratorio.
   Fase 3 (13/10): il picco del blocco 7 viene dalla cucitura (alone 2 lo sposta di 1–3·10⁻⁴), non dallo stato iniziale;
   convergenza ancora FAIL, corrente spuria non transitoria (artefatto della cucitura nelle celle coperte, 57× il pavimento);
   regola fase 4 adottata (asset P4_v1): 2p +0,0338 % come previsto (massa = quella del maestro). Prossimo: cucitura che
   rispetti la fase, maestro con massa corretta.
3. **Elio e litio esatti** sulla zig-zag (o sulla regola agglomerata), poi imparare da loro la
   regola effettiva a molti elettroni. Esame: energie di ionizzazione e righe.
4. Atomi uno a uno, numero atomico per numero atomico. Dopo: legami, H₂, acqua.
5. Parete di atomi agglomerati → doppia fenditura vera → fotone che la attraversa (servono
   ≥ 24 celle per λ e un riferimento di continuo, vedi `ATOMO.md` §5).
6. Scena gamma → elettrone + positrone vicino a un nucleo pesante (vuoto di coppie in 2D HD).
7. **Programma unico** (richiesta di Valerio 13–14/10): tutto sul renderer dell'atomo zig-zag (bagliori con scia corta), modalità tavola periodica + Impostazioni grafiche + modalità Modello standard; le scene di molecole, elio, quark e gluoni vanno portate lì (vedi riga Modello standard a palline in §4).

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
