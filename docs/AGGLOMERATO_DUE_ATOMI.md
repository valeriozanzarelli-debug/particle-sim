# Il livello delle due atomi: l'agglomeratore pronto per la scala di H₂

Studio `aggmol_`. Continua `AGGLOMERATO_ZIGZAG.md` e `AGGLOMERATO_COMPLETO.md`. **Qui non si calcola il legame di H₂**
(le molecole vengono dopo la regola a molti elettroni e dopo l'elio, decisione di Valerio): si **prepara e si collauda la
macchina** che servirà, alla scala di due atomi (protoni a 1,4 a₀, nube su 5–10 a₀).

File: motore `lab/aggmol_core.py` (campo di più nuclei, Hamiltoniano grosso con settori di parità), soluzione esatta del
problema a due centri (solo esame) `lab/aggmol_exact.py`, apprendimento `lab/aggmol_learn.py`, un atomo
`lab/aggmol_atom.py`, due centri `lab/aggmol_h2p.py`, costi `lab/aggmol_costs.py`, esiti `lab/aggmol_analyze.py`,
figura `lab/aggmol_figure.py`; prove dopo il FAIL `lab/aggmol_posthoc16.py`; sviluppo `lab/aggmol_dev.py`.
Copie congelate del codice di `aggfull_` (che lavora in parallelo sui suoi file): `lab/aggmol_fz_*.py` (commit a70a372).
Pre-registrazione `lab/aggmol_prereg.py` → `lab/out/aggmol_prereg.json`, **commit e79e57a, prima delle corse di
produzione**. Esiti: `lab/out/aggmol_verdict.json`, `aggmol_learn.json`, `aggmol_atom.json`, `aggmol_h2p.json`,
`aggmol_costs.json`, `aggmol_posthoc16.json`.

![figura](../lab/out/aggmol_figure.png)

---

## 0. In breve

**Risultato: 14 PASS / 1 FAIL.**

1. **La macchina ha scelto da sola il livello a₀/4** (cella di un quarto di raggio di Bohr) come "livello delle due
   atomi", e ha **rifiutato a₀/2**: lì la cella è così grossa che la regola non riesce più a imitare il livello sotto per
   gli impulsi dell'1s (errore 1,0–1,3·10⁻³ contro il budget 3·10⁻⁴, da tutte e due le strade). Il pavimento per gli
   atomi è a₀/4.
2. **Un atomo con a₀/4** (stessi vincoli esatti della fase 3: Taylor fissato, parità; positività impossibile come
   previsto dal limite dei momenti): 1s −0,088 % dalla regola fine, 2s +0,011 %, 2p +0,008 %; struttura fine del 2p
   44,94 μeV contro 45,28 di Dirac (−0,77 %), rapporto di Landé 0,5000. Darwin e spin-orbita riletti dalla cella imparata
   entro l'1 %.
3. **Collaudo a due centri** (H₂⁺: **un solo elettrone** fra due protoni fermi; collaudo della macchina, non chimica)
   contro la soluzione esatta: a 64³ celle l'energia è giusta entro **0,38 mEh** a R = 2 a₀ e **0,61 mEh** a R = 1,4 a₀
   (sotto l'accuratezza chimica, 1,59 mEh); legame D_e 102,52 contro 102,63 mEh; distanza di equilibrio +0,15 %; tre stati
   eccitati entro 0,04 %; la molecola si può spostare e ruotare sulla griglia cambiando l'energia di ≤ 4·10⁻⁵.
4. **La cosa più bella (idea geometrica, provata):** il nodo dello stato 2σg esce un **ellissoide con i due protoni nei
   fuochi**, come vuole la soluzione esatta (separabile in coordinate sferoidali): sul nodo ξ = (r_a + r_b)/R varia dello
   **0,08 %**, mentre il raggio varia del 5,7 %. Nessuno ha messo ellissi nella regola cubica: la geometria nascosta
   del problema a due centri emerge da sola.
5. **FAIL (C8):** la scala a₀/4 → a₀/8 → a₀/16 non scende in modo monotono: il punto a₀/16 l'ho messo in una scatola ±6 a₀
   per stare nella memoria, e la scatola periodica (il "cristallo di atomi") lega di +4,4·10⁻⁴ in più. Errore di disegno
   dell'esame, verificato dopo (§5).
6. **La legge dei costi** (§1) dice la cosa che conta per H₂ vero: con un elettrone la scatola di due atomi a a₀/4 è
   **64³ = 2,6·10⁵ celle** (22 MB, ~20 ms per passo stimati sulla Intel UHD); con **due elettroni esatti** lo stato vive
   in 6 dimensioni e a a₀/4 sono 275 GB. **Non entra in nessuna scheda**: serve la regola a molti elettroni (§6).

---

## 1. Compito 1: la legge dei costi per una scatola di due atomi

Scatola ±8 a₀ (lo stato fondamentale con i due protoni a ±0,7 a₀; la nube di 5–10 a₀ ci sta con margine) e ±12 a₀ (stati
eccitati). Tempi **misurati** qui (un core, numpy) per un'applicazione dell'Hamiltoniano grosso (cinetico imparato + campo
a 8 lanterne); stime su hardware con il modello di banda di `MOTORE_OTTIMIZZATO.md` §8 (Intel UHD ~25 GB/s, schede da 16 GB
288–717 GB/s: numeri di catalogo, **non misurati**; un passo con il campo ≈ 18 FFT + 9 passate, ~1,9 kB per cella).

| cella | scatola | N³ celle | memoria (dinamica / autostati) | H, CPU 1 core | Intel UHD (stima) | 16 GB alta (stima) | passi per unità atomica di tempo | H₂ esatto a 2 elettroni (6-D, singoletto) |
|---|---|---|---|---|---|---|---|---|
| a₀/2 | ±8 a₀ | 32³ = 3,3·10⁴ | 3 MB / 6 MB | 6 ms | 2,5 ms | 0,09 ms | 274 | 4,3 GB |
| **a₀/4** | **±8 a₀** | **64³ = 2,6·10⁵** | 22 MB / 44 MB | 60 ms | **20 ms** | 0,7 ms | 548 | **275 GB** |
| a₀/4 | ±12 a₀ | 96³ = 8,8·10⁵ | 74 MB / 0,15 GB | 0,24 s | 66 ms | 2,3 ms | 548 | 3,1 TB |
| a₀/8 | ±8 a₀ | 128³ = 2,1·10⁶ | 0,18 / 0,35 GB | 0,61 s | 157 ms | 5,5 ms | 1096 | 18 TB |
| a₀/8 | ±12 a₀ | 192³ = 7,1·10⁶ | 0,6 / 1,2 GB | 3,1 s | 530 ms | 18 ms | 1096 | 200 TB |
| a₀/16 | ±8 a₀ | 256³ = 1,7·10⁷ | 1,4 / 2,8 GB | (~7 s) | 1,3 s | 44 ms | 2193 | 1,1 PB |
| a₀/16 | ±12 a₀ | 384³ = 5,7·10⁷ | 4,8 / 9,5 GB | — | 4,2 s | 148 ms | 2193 | 13 PB |

Lettura semplice:
* **Con un elettrone il costo per tempo fisico va come h⁻⁴** (celle h⁻³ per passi h⁻¹). Il battito σg↔σu di H₂⁺ a
  R = 2 a₀ dura 14,4 unità atomiche (0,35 fs): a a₀/4 sono ~7 900 passi nativi → **~2,5 minuti sulla Intel UHD**,
  ~5 s su una scheda alta da 16 GB. A a₀/8: 16 volte di più. È la "mela a 480p e in 4K" per gli atomi: Intel a₀/4,
  scheda grossa a₀/8 (o a₀/16 solo vicino ai nuclei, con lo zoom).
* **Sotto a₀/4 non si scende** (la macchina lo rifiuta, §2): è il pavimento degli atomi, come le 12 celle/λ lo sono per
  la luce.
* **Due elettroni esatti** sono un altro mondo: lo stato è una funzione di 6 coordinate, M² numeri. Anche alla cella più
  grossa ammessa (a₀/4) sono 275 GB. L'unica casella sotto i 16 GB è a₀/2 ±8 a₀ (4,3 GB), che però la macchina ha
  rifiutato. **Quindi H₂ esatto sulla griglia non si fa a forza bruta**: serve la regola a molti elettroni imparata (coppie,
  Jastrow, camminatori), come pianificato in `MOLTI_ELETTRONI_PIANO.md`.

## 2. Compito 2: imparare il livello successivo con la macchina automatica

Si parte dal livello a₀/8 di `aggfull_` (cinetico T: c₂, c₄, c₆ fissati esattamente al livello sotto, parità) e dal suo
accoppiamento cella-campo imparato. Candidati, dal più grosso: b = 4 (a₀/2), poi b = 2 (a₀/4); poi **si ripete** dal
livello scelto. Per ogni candidato: impara (stesso codice della fase 3–4), controlla il moto libero, impara come la cella
sente il campo (Duhamel dal livello **imparato** sotto), controlla l'1s contro il livello sotto (accuratezza chimica).

| esame | cosa | risultato | esito |
|---|---|---|---|
| L1 | la macchina sceglie a₀/4 e rifiuta a₀/2 (previsto dallo sviluppo) | a₀/4: errore energia 2,2·10⁻⁴ ≤ 3·10⁻⁴, 1s +0,44 mEh dal livello sotto: ammesso. a₀/2: 0,97 / 1,30·10⁻³ (da a₀/8, w = 4/5) e 1,02 / 1,31·10⁻³ (da a₀/4) | **PASS** |
| L2 | positività (regola di Markov) impossibile, dal limite dei momenti | deviazione dei vincoli: 283 (a₀/4), 1133 (a₀/2); servirebbe \|c₄\|/c₂ = 2·10⁻⁴ (a₀/4), i salti positivi danno ≥ 1/12 | **PASS** |
| L3 | Darwin e spin-orbita riletti dalla cella a₀/4 entro 2 % | stencil a stella scelto da solo; g_k −0,12618 (+0,94 %), g_x +0,25100 (+0,40 %), g_so +0,25256 (+1,02 %) | **PASS** |

**Perché a₀/2 non va (fisica, non sfortuna).** La cella grossa vede solo impulsi fino a π/h. A a₀/2 il bordo della zona
è a k = 2π/a₀, e la "coda" dell'1s (densità degli impulsi ∝ 1/(1+k²a₀²)⁴) ha ancora peso lì. Il livello sotto, visto
attraverso b passi, ha una banda che **non è periodica** con il periodo della nuova cella: vicino al bordo la regola
grossa non può seguirla, e l'errore sull'energia di moto pesata sull'1s sale a 10⁻³. Con il cinetico esatto (sviluppo,
dichiarato): 1s −0,68 %, H₂⁺ −0,42 % a a₀/2, contro −0,094 % e −0,057 % a a₀/4.

**Vincoli esatti (come la fase 3).** Massa: m_kin = tanθ/θ = 1,000373 portata esatta (per costruzione, non conta come
esame). c₄ della regola fine riportato giusto (1,0001 di Dirac); il c₆ letto dal livello a₀/8 è −2,4 volte Dirac: a
queste celle il k⁶ vale ~10⁻¹⁰ dell'energia, è rumore del livello sotto (dichiarato, senza effetto misurabile).

## 3. Esami su un atomo con il livello a₀/4

Scatola ±16 a₀ (128³), regola grossa + primo ordine del resto dell'accoppiamento imparato (Darwin, spin-orbita). Confronto
con la regola fine radiale allo stesso θ (`aggfull_fineref.json`) e con Dirac.

| stato | a₀/4 (eV) | contro la regola fine | contro Dirac |
|---|---|---|---|
| 1s½ | −13,597926 | −0,088 % (+0,44 mEh) | −0,058 % (+0,29 mEh) |
| 2s½ | −3,402979 | +0,011 % (−0,013 mEh) | +0,044 % |
| 2p½ | −3,402188 | +0,008 % | +0,021 % |
| 2p3/2 | −3,402143 | +0,008 % | +0,021 % |

| esame | soglia | esito |
|---|---|---|
| A1 | i quattro stati entro 0,25 % e 1 kcal/mol della regola fine | **PASS** |
| A2 | idem contro Dirac | **PASS** |
| A3 | struttura fine 2p entro 2 % di Dirac, rapporto 0,5 entro 1 % | 44,94 μeV contro 45,28 (−0,77 %), 0,5000: **PASS** |

L'1s è il più sbagliato (−0,088 %): è la cuspide del nucleo, che una cella di a₀/4 liscia (stessa lezione di
`AGGLOMERATO_COMPLETO` §4.3). Gli stati n = 2, che al nucleo pesano 8 volte meno, sono a 10⁻⁴.

## 4. Compito 3: collaudo a due centri (H₂⁺, un elettrone, protoni fermi)

**Dichiarazione.** H₂⁺ ha un solo elettrone: è ammesso **come collaudo della macchina a due centri**. Non è il legame di
H₂ (che ha due elettroni e la loro repulsione e correlazione: §6). I protoni sono fermi (Born–Oppenheimer): il loro moto è
un debito. L'unico ingresso nuovo rispetto all'atomo è **la posizione del secondo protone**: il campo grosso è la somma
dei due campi di Galerkin (−α/r di ciascun protone), le posizioni possono essere ovunque, non agganciate alla griglia.

**Chiave d'esame.** La soluzione esatta non relativistica del problema a due centri, che si separa in coordinate
sferoidali (ξ, η, φ) ed è calcolata qui con due equazioni a una variabile (Chebyshev). Controllata: R = 2, 1σg =
−1,1026342144949 Eh (il numero classico), 1σu −0,6675343922, 2σg −0,3608648753, 1πu −0,4287718199.

| esame | cosa | risultato | esito |
|---|---|---|---|
| C1 | 1σg a R = 1,4 e 2,0 a₀, a₀/4 (64³), entro 1 kcal/mol e 0,1 % | +0,61 mEh (−0,048 %) e +0,38 mEh (−0,034 %) | **PASS** |
| C2 | D_e a R = 2 (stesso livello e scatola per l'atomo) entro 1 kcal/mol | 102,517 contro 102,634 mEh (−0,12 mEh) | **PASS** |
| C3 | curva R = 1,6…2,4: R_e entro 0,5 %, curvatura entro 2 % (stesso fit sui valori esatti) | R_e +0,15 %, k_e −0,95 % | **PASS** |
| C4 | 1σu, 1πu, 2σg a R = 2 (±12 a₀) entro 0,1 % | −0,041 %, +0,031 %, +0,011 % | **PASS** |
| C5 | spostare la coppia di frazioni di cella: ≤ 10⁻⁴ | 3,9 / 2,7 / 1,9·10⁻⁵ | **PASS** |
| C6 | coppia lungo la diagonale (1,1,1) invece che lungo z: ≤ 10⁻⁴ | 1,0·10⁻⁵ | **PASS** |
| C7 | nodo del 2σg = ellissoide confocale | spread di ξ 0,081 % (del raggio 5,7 %); ξ₀ = 1,90745 contro 1,90781 (−0,04 % di ξ₀−1) | **PASS** |
| C8 | scala a₀/4 → a₀/8 → a₀/16 contro l'esatto con la massa della regola: errore decrescente, ≤ 0,02 % a a₀/16; D_e entro 1 kcal/mol ovunque | −5,5·10⁻⁴, −2,1·10⁻⁵, **+4,1·10⁻⁵** (non monotono); D_e −0,12 / +0,02 / −0,17 mEh | **FAIL** |
| C9 | salto dell'elettrone fra i protoni: σu − σg a R = 4 e 6 a₀ entro 1 % | +0,055 % e +0,13 %: periodo 62,46 contro 62,50 e 294,2 contro 294,6 unità atomiche (1,5 e 7,1 fs) | **PASS** |

L'errore della curva (inserto del pannello b) è quasi costante, +0,3…+0,6 mEh: è soprattutto l'errore di cuspide dei due
nuclei, che è lo stesso dell'atomo e in gran parte si cancella nel legame (D_e sbaglia di 0,12 mEh soltanto).

## 5. Fallimenti (dichiarati)

**C8.** Il punto a₀/16 l'ho messo, nella pre-registrazione, in una scatola ±6 a₀ (192³) per restare sotto i 2,5 GB.
Prova dopo il FAIL (`aggmol_posthoc16.py`, non pre-registrata): con il cinetico **esatto** k²/2m_kin, stessa scatola e
stesso campo, l'1s a a₀/16 sta a +4,4·10⁻⁴ dalla chiave; e a a₀/8 lo stesso cinetico esatto dà +3,3·10⁻⁴ a ±6 a₀ e
−1,1·10⁻⁴ a ±8 a₀. Cioè ~4,5·10⁻⁴ sono **la scatola periodica** (il cristallo di atomi di `IDROGENO_BOHR.md`), non il
livello: il cinetico imparato a a₀/16 contribuisce solo −3,4·10⁻⁵ Eh. Rifare a a₀/16 in ±8 a₀ (256³) vuole ~2,8 GB: sopra
il limite di questo container (provato, `MemoryError`). Il FAIL resta. Lezione (già nota, ripetuta): **la scatola è un
livello del passaporto d'errore**.

Errore mio di sviluppo, trovato e corretto prima della pre-registrazione: nel cercatore di nodi l'indice della
riempitura di zeri veniva arrotondato male (spostava il nodo del 5 % in un piano).

## 6. Compito 4: cosa manca per H₂ vero

H₂ ha **due elettroni**. Rispetto al collaudo di questo studio mancano:

1. **La regola a molti elettroni** (thread principale, `MOLTI_ELETTRONI_PIANO.md`): lo stato di due elettroni vive in 6-D;
   a a₀/4 sono 275 GB (§1). Serve la regola di coppia imparata da corse esatte di 2 elettroni (elio, H⁻) e controllata a
   campione. **La macchina a due centri di qui è pronta a riceverla**: il campo di più nuclei, i settori di simmetria, la
   covarianza per traslazione e rotazione (C5, C6) sono collaudati.
2. **La repulsione elettrone-elettrone sulla griglia grossa**: va imparata come la cella sente il campo di un'altra carica
   che si muove (stessa tecnica dell'accoppiamento cella-campo, ma con un campo che non è fisso). Oggi c'è solo il campo
   dei nuclei.
3. **I nodi** non servono per H₂ nello stato fondamentale (singoletto senza nodi: la parte spaziale è simmetrica e
   positiva), quindi H₂ è il primo caso dove i camminatori sono esatti. Ma i camminatori vogliono un nucleo **positivo**, e
   L2 dice che a a₀/4 una regola positiva non può portare il pezzo relativistico: per H₂ (non relativistico a 10⁻⁵) si può
   usare un nucleo positivo che porta solo la massa, dichiarando l'errore k⁴ del reticolo; va misurato.
4. **I protoni che si muovono**: qui sono fermi. `AGGLOMERATO_COMPLETO` ha già il protone come pallina pesante imparata
   (rinculo entro 0,6 %): per la vibrazione di H₂ serve la coppia protone-protone con la stessa macchina.
5. **La cuspide**: a a₀/4 l'errore principale è ai nuclei (+0,3…0,6 mEh). Per la chimica fine serve lo zoom ai nuclei
   (`nucleo1s_`, `IDROGENO_BOHR.md`), con la colla già provata lì.
6. **La scatola**: a ±6 a₀ costa 4·10⁻⁴; per due atomi ±8 a₀ è il minimo.

## 7. Debiti

**Chiusi o avanzati:** CONTESTO §4 voce 11 ("livelli più alti a₀/4, a₀/2 per molecole e gas") → **a₀/4 imparato e
collaudato; a₀/2 rifiutato dalla macchina per gli atomi** (resta possibile per oggetti più larghi dell'1s: gas, stati di
Rydberg, molecole viste da lontano).

**Aperti / nuovi:**
- protoni fermi (Born–Oppenheimer) nel collaudo a due centri;
- iper-parametri ereditati (peso uniforme 3·10⁻⁴, budget 3·10⁻⁴, osservabile "1s di taglia a₀", soglia dell'accuratezza
  chimica per il controllo dell'atomo: quest'ultima è una scala esterna nota, non un fit);
- il campo dei due protoni è la **somma** dei due campi di Galerkin: corretto per linearità, ma l'accoppiamento cella-campo
  è imparato per un campo debole e applicato al primo ordine (Darwin e spin-orbita) solo sull'atomo, non su H₂⁺
  (effetti ~10⁻⁵, dichiarati);
- c₆ del livello a₀/8 non di Dirac (−2,4×; irrilevante qui);
- C8: rifare a₀/16 in ±8 a₀ su una macchina con più memoria.

**Parametri liberi:** nessuno nuovo di fisica. Ingressi: m_e (θ), α, Z = 1 per protone, **posizioni dei protoni**. Numeri
imparati: livello a₀/4, 729 posti × 2×2 complessi del cinetico (fissati 40 vincoli esatti dal livello sotto) + stella
2×2×7×7 dell'accoppiamento. Nessun fit contro l'esatto o contro Dirac: le chiavi esatte sono solo esami.

## 8. Idee geometriche e visive

**V1 — I nodi confocali (provata: C7 PASS).** *A schermo:* gli stati eccitati della coppia di protoni mostrano superfici
dove la nube si annulla; sono **ellissoidi e iperboloidi con i due protoni nei fuochi**, come le orbite dei pianeti
attorno a due soli. Disegnate come palline che "scivolano" lungo la superficie, si vedrebbero ellissi perfette che nessuno
ha disegnato. *Perché:* il problema a due centri ha una simmetria nascosta (una costante del moto in più, che lo rende
separabile in coordinate sferoidali); la regola cubica non la conosce, ma la ritrova: ξ costante entro 0,08 %. *Costo:*
zero (si leggono gli stati che già ci sono). *Regole:* sì, è fisica emersa, nessuna forma disegnata.

**V2 — L'elettrone che salta fra i protoni (calcolata: C9 PASS).** *A schermo:* metti l'elettrone su un protone: la nube
passa all'altro e torna, avanti e indietro, con periodo 62 unità atomiche a R = 4 a₀ e 294 a R = 6 a₀ (1,5 e 7,1 fs):
più lontani, più lento (esponenzialmente). Rallentato al ritmo che vuole Valerio, è il "tunnel" visto con gli occhi. *Perché:*
è il battimento fra σg e σu, la cui differenza esce entro 0,13 %. *Costo:* a a₀/4 con 80³ celle, sulla Intel ~1 minuto di
calcolo per un salto a R = 4 (o la "cartolina" a due stati: istantanea). *Regole:* sì.

**V3 — Molecola che scivola sulla griglia senza scatti (provata: C5, C6).** *A schermo:* la molecola che si sposta o ruota
non "scatta" da una cella all'altra. *Perché:* le 8 lanterne per cella rendono il campo covariante: spostare la coppia di
una frazione di cella cambia l'energia di ≤ 4·10⁻⁵ (40 μEh), ruotarla di 10⁻⁵. È la condizione per far muovere i nuclei
(vibrazioni, urti) senza forze finte di griglia: con 4·10⁻⁵ Eh di "scatola delle uova" la forza finta è ~10⁻⁴ Eh/a₀,
mille volte sotto la forza di legame. *Costo:* già pagato. *Regole:* sì.

**V4 — Livelli come lenti a fuoco variabile (proposta).** *A schermo:* lontano dai nuclei celle a₀/4, attorno a ogni
nucleo un "occhio" a₀/16 o lo zoom zig-zag; fra due protoni il ponte resta a₀/4. *Perché:* l'errore di C1 è tutto nelle
cuspidi (stesso +0,3 mEh per atomo e per molecola); lo zoom locale lo toglie senza moltiplicare le celle (a₀/16 dappertutto
costa 64 volte). *Costo:* la colla di `IDROGENO_BOHR.md`. *Regole:* sì.

**V5 — La bolla di legame (proposta per il thread grafico).** *A schermo:* la densità di 1σg meno quella dei due atomi
separati: una ciambella luminosa fra i protoni dove l'elettrone si accumula quando lega (pannello c mostra la nube; la
differenza ne mostra il "perché"). *Costo:* una sottrazione. *Regole:* sì, è un'immagine di un dato calcolato.

## 9. Dove l'intuizione di Valerio può aiutare

1. **Il pavimento a₀/4.** Per la "mela a 480p" degli atomi la macchina non scende sotto a₀/4. Ti va bene che sulla Intel
   la differenza sia la **grandezza della scena** (meno atomi) e non la risoluzione? O vuoi che proviamo a₀/2 solo per
   oggetti grandi (gas, atomi lontani) dove l'1s non si guarda?
2. **Due elettroni.** La regola di coppia: come immagini la repulsione fra due palline-elettrone nella regola zig-zag? Come
   due correnti che si respingono (campo dell'una che passa per le celle dell'altra), o come un "costo" quando stanno nella
   stessa cella? Le due cose corrispondono a pezzi diversi di fisica nota (Coulomb lungo raggio vs cuspide di coppia).
3. **Il salto fra protoni (V2).** È la versione più semplice della chimica che "si vede". Vuoi che sia la prima scena del
   laboratorio a due atomi, prima di H₂?
4. **I protoni che si muovono.** Con V3 la griglia non ostacola il moto dei nuclei. Vuoi H₂⁺ che vibra (protoni come
   palline pesanti della stessa regola) come collaudo successivo, sempre con un elettrone?

## 10. Patch proposte per altri file (non toccati)

- `docs/CONTESTO_CLOUD.md` §4 voce 11 → "a₀/4 imparato e collaudato (AGGLOMERATO_DUE_ATOMI: 1s −0,088 %, H₂⁺ entro
  0,6 mEh); a₀/2 rifiutato dalla macchina per gli atomi". Indice: aggiungere AGGLOMERATO_DUE_ATOMI.
- `docs/PIANO.md` §3: riga "Livello delle due atomi (macchina a₀/4, collaudo a due centri a un elettrone)".
- Per il Laboratorio Atomo (Rust, non toccato): il livello a₀/4 è in `lab/out/aggmol_levels.pkl` (non nel repo; si rifà con
  `python lab/aggmol_learn.py`, ~3 minuti); un passo = due FFT del cinetico + 8 lanterne per il campo, 64³ celle.
- **Ponte (proposta, non creata):** `ponte/richieste/aggmol_salto.json`, azione `foto`, scena H₂⁺ a R = 4 a₀ con
  l'elettrone che parte su un protone (V2), una foto ogni 8 unità atomiche.

## 11. Come rieseguire

```
cd lab
python aggmol_prereg.py                       # (già committato, e79e57a)
python aggmol_learn.py                        # livelli a0/4 (e il rifiuto di a0/2), ~3 min
python aggmol_atom.py s p pert_s pert_p       # un atomo a a0/4
python aggmol_h2p.py atom_4 atom_8 curve_4_R2.0 ...   # vedi l'elenco in aggmol_h2p.py / lab/out/aggmol_prod.log
python aggmol_h2p.py gap_4_R4.0 gap_4_R6.0
python aggmol_costs.py; python aggmol_posthoc16.py
python aggmol_analyze.py; python aggmol_figure.py
```
