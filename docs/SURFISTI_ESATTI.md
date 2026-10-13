# Surfisti esatti: palline che cavalcano l'onda senza barare

Studio del 13/10/2026. Prefisso `surf_`. Pre-registrazione: `lab/surf_prereg.py` → `lab/out/surf_prereg.json`
(commit 2067d27, prima di ogni corsa di produzione). Emendamenti dichiarati: `lab/out/surf_amend.json`.
Figura: `lab/out/surf_figure.png`. Tabellone: `lab/out/surf_results.json`.

## 0. In breve

**Il debito** (CONTESTO_CLOUD §4, #7): i surfisti (le palline che mostrano "dove sta" l'elettrone) (a) avevano un
lettore impreciso (equivarianza 1,26·10⁻²), (b) non c'era una guida deterministica che conservasse |ψ|²
*esattamente* sul reticolo, (c) il salto a caso creava mescolamento finto.

**La soluzione trovata, in una frase.** La regola zig-zag sposta la probabilità solo attraverso i link, e sa
*esattamente* quanta ne passa da ogni link a ogni tick (la parte che "va avanti" meno quella che "torna indietro").
Usiamo quella corrente esatta come un fiume che scorre tra le celle, e muoviamo ogni surfista *dentro* la cella con
la soluzione analitica di quel fiume (il metodo di Pollock dei geologi che seguono le gocce d'acqua nelle falde,
in "tempo di flusso"). Risultato: **|ψ|² è conservato esattamente, con zero numeri casuali**, e la regola si può
**riavvolgere** esattamente all'indietro.

**I numeri principali** (tutti contro la pre-registrazione):

| cosa | vecchio | nuovo (flusso esatto) |
|---|---|---|
| F6 (1-D, il test che falliva): errore sul "posto in fila" | 1,26·10⁻² (riprodotto: 1,260·10⁻²) | **2,9·10⁻¹⁵** |
| idrogeno 1s (Z=1) fermo, 10⁴ tick: quanto si muove un surfista | lettore vecchio: z = 104 (deriva), Δq 6,8·10⁻² | **Δq 7·10⁻¹⁴, zero salti** |
| atomo 3-D vivo (Z=60, 48³), 1s+2p che batte, 10 periodi: \|z\| massimo | **60 327** (Rust attuale), 28 949 con le rinascite | **2,2** (soglia 4) |
| onda libera 3-D: \|z\| massimo | 24 218 | **2,5** |
| riavvolgimento 20 tick avanti + 20 indietro | non reversibile | **≤ 1,3·10⁻¹² celle** |
| mescolamento finto (partenza sbagliata in uno stato fermo, 1s 3-D) | Bell: H̄ cala all'**15 %** in 1000 tick | **98 %** (flusso a due tick): non mescola |

Esiti: **PASS 47, FAIL 1** (con l'emendamento A3 sul riferimento del TV; senza: 44 PASS, 4 FAIL — §6.2) (dettagli in §4). Il FAIL vero è uno: la velocità di rotazione del surfista
dell'1s nel reticolo 3-D grosso (0,21 invece di ~1 volte Zα c sinθ). Non è colpa della regola dei surfisti: è lo
stato 1s del reticolo a θ = 0,5 che ha lo spin mezzo "sciolto" (⟨Σz⟩ = 0,51 invece di 0,99) — §6.

**Validazione o scoperta?** Validazione e ingegneria. La matematica usata è nota (correnti di Bell 1984,
processo di Vink 1993, tracciamento di Pollock 1988, elementi di Raviart–Thomas 1977). È nuovo, per quanto so,
metterli insieme sul Dirac quantum walk: la corrente esatta per link della zig-zag + flusso RT0 in tempo di flusso
dà surfisti deterministici, esatti e reversibili. Nessuna previsione fisica nuova.

---

## 1. Il problema, come a un ventenne

Un surfista è una pallina che sta "sopra" l'onda dell'elettrone. Se ne lanci tante, si devono distribuire come
|ψ|² (la regola di Born), *e restarci* mentre l'onda si muove. Questa proprietà si chiama **equivarianza**.

Nel continuo (de Broglie–Bohm) basta dare ai surfisti la velocità v = j/ρ (corrente diviso densità). Sul reticolo
non è così semplice:

* il vecchio lettore (Python F6 e Rust `zz.rs::muovi_surfisti`) calcola v = ψ†αψ/ψ†ψ dallo spinore
  **interpolato** nel punto e fa un passo. Ma la zig-zag non è un'onda liscia: in ogni tick metà dell'onda salta
  avanti di una cella e metà indietro. La corrente vera vive **sui link fra due celle**, sfasata di mezza cella e
  di mezzo tick rispetto a dove il lettore la va a leggere. L'errore è di ordine a, e si accumula;
* i surfisti a salti casuali di `pre_core.py` erano esatti in media ma **rimescolavano** la folla da soli (una
  catena di Markov si mescola anche se l'onda è ferma): "mescolamento finto";
* la mappa monotona per asse (`redo_r3.py`, R3) era esatta ma, facendo x poi y poi z a ogni mezzo passo,
  faceva "ballare" i surfisti al ritmo dello zig-zag e lasciava un mescolamento lento.

## 2. L'idea: il fiume tra le celle

### 2.1 La corrente esatta della zig-zag

In un mezzo passo lungo l'asse j, ogni cella x divide il suo stato: la parte A(x) che va avanti attraversa il link
(x, x+1), la parte B(x+1) della cella vicina lo attraversa all'indietro. Il moneta (flip) agisce *dentro* la cella e
non sposta |ψ|². Quindi la **corrente netta per link**

    F_j(x) = A_j(x) − B_j(x+e_j)

è esatta, e la densità soddisfa ρ_dopo = ρ_prima − div F **a precisione macchina** (verificato: 9·10⁻¹⁶). È la
corrente di Bell (1984) per la nostra regola. Sommando i tre mezzi passi (o i sei di una coppia di tick, il periodo
vero della regola: ordine x,y,z poi z,y,x) si ha il flusso netto del tick.

### 2.2 Dentro la cella: un fiume che si può risolvere a mano

Il flusso è noto sulle sei facce della cella. Dentro, lo riempiamo nel modo più semplice compatibile: ogni
componente varia linearmente lungo il suo asse (è l'elemento di **Raviart–Thomas** di ordine più basso). La densità
della cella passa linearmente da ρ_t a ρ_{t+1}. Due fatti geometrici:

1. **dentro ogni cella** la continuità è esatta (la divergenza del campo RT0 è costante e uguale a −dρ/dt);
2. **sulle facce** il flusso normale è continuo (la faccia è condivisa).

Quindi la legge di conservazione vale ovunque, e il moto dei surfisti trasporta la misura "uniforme a pezzi" ρ_t
**esattamente** su ρ_{t+1}. Niente approssimazioni.

**Il trucco del tempo di flusso.** Con τ = ∫ ds/ρ(s) ogni coordinata si muove indipendentemente:
q_j(τ) = q_j(0)·e^{d_j τ}, dove q_j è il flusso nel punto e d_j la differenza tra le due facce. Si risolve con un
logaritmo il tempo in cui il surfista tocca una faccia, si passa alla cella vicina, e si continua fino a fine tick.
È esattamente l'algoritmo di Pollock (1988) con il "magazzino" che cambia nel tempo. Costo: qualche operazione per
surfista per tick (0,06 s per 2·10⁵ surfisti in un solo thread).

**Proprietà che vengono gratis:** deterministico (nessun numero casuale dopo la nascita da |ψ|²); in 1-D è la mappa
monotona (unica), quindi coincide con R3; **reversibile**: con −F da ρ_{t+1} a ρ_t ritorna esattamente indietro;
**nessuna manopola** (0 parametri continui).

### 2.3 Gli schemi messi a confronto

| sigla | nome | casuale? | esatto? |
|---|---|---|---|
| S0 | lettore vecchio (Rust `zz.rs`) | no | **no** |
| S0+ | lettore vecchio + 0,2 % di rinascite per tick (come `zzscena.rs`) | sì | no |
| S1 | Bell balistico (salti con i flussi lordi: tutti saltano come l'onda) | molto | sì |
| S2 | Bell minimo per mezzo passo (`pre_core` λ = 0) | sì | sì |
| S2b | Bell in tempo di flusso (salti continui con il flusso netto del tick) | sì | sì |
| S3 | mappa monotona per asse (R3) | no | sì |
| **S4** | **flusso esatto, 1 tick** | **no** | **sì** |
| **S4-2** | **flusso esatto sul periodo di 2 tick** (consigliato) | **no** | **sì** |
| S5 | orbitali analitici: corrente di Dirac sulle facce + correzione di Poisson minima, poi S4 | no | sì |

## 3. Esami pre-registrati e corse

Tutte le corse: un thread, pezzi ≤ 10 min con checkpoint. L'atomo 3-D è **proprio quello del Laboratorio Atomo**
(`src/atomo/zz.rs`: Z = 60, 48³ celle, a = 0,5 λ_C, θ = 0,5). Gli stati 1s e 2p½ sono stati ottenuti *ascoltando* la
regola (filtro di Hann, due passate da 4096 tick, solo tick pari perché il periodo della regola è di due tick):
1s a ω = 0,89945 (residuo su due tick 1,5·10⁻⁴), 2p½ a ω = 0,97553 (6,6·10⁻⁴). Battimento 1s+2p: 82,6 λ_C/c.

Statistica: istogramma dei surfisti contro |ψ|² su blocchi 2³ (3-D), z = (χ² − gradi)/√(2·gradi), 41 controlli per
corsa; un campionatore ideale dà |z| di ordine 1.

## 4. Risultati contro la pre-registrazione

| esame | cosa | soglia | risultato | esito |
|---|---|---|---|---|
| E0a | cammino con corrente = `zz_core.Walk`; ρ₁ = ρ₀ − div F | 10⁻¹⁴ / 10⁻¹³ | 0 / 9·10⁻¹⁶ | PASS |
| E0b | Liouville su un tick (ρ·J costante), stati 3-D 1s / 1s+2p / libero | mediana ≤ 10⁻⁵, 99 % ≤ 10⁻³ | 8·10⁻⁹ / 4,7·10⁻⁷ / 3,7·10⁻⁹ (99 %: 5·10⁻⁵ / 8·10⁻⁵ / 6·10⁻⁸) | PASS ×3 |
| E0c | riavvolgimento 20+20 tick (celle) | mediana ≤ 10⁻¹², 99 % ≤ 10⁻⁹ | ≤ 3,6·10⁻¹⁵ / 99 % ≤ 1,3·10⁻¹² | PASS ×3 |
| E1 | F6: 1-D libero θ = 0,1, 40 ħ/mc², 2000 surfisti, errore sul posto in fila | ≤ 10⁻¹² | **2,9·10⁻¹⁵** (vecchio: 1,260·10⁻²) | PASS |
| E2 | limite continuo 1-D (a = 0,2/0,1/0,05): distanza dai surfisti di Bohm–Dirac continui | decrescente, ordine ≥ 0,8 | 0,251 / 0,111 / 0,052 λ_C, ordine **1,10**; identità dei quantili di Bohm continuo 4,6·10⁻⁷ | PASS ×3 |
| E3 | idrogeno vero (Z = 1) 1s radiale, 10⁴ tick: surfisti fermi | Δq ≤ 10⁻⁹; salti ≤ 10⁻⁶ per tick | Δq **7·10⁻¹⁴**; Bell minimo e in tempo di flusso: **0 salti** (balistico: 1 salto per tick; vecchio: Δq 0,068, z = 104) | PASS ×3 |
| E4 | radiale Z = 60, 1s+2s, 20 battimenti (33 160 tick) | quantili ≤ 10⁻¹⁰; \|z\| ≤ 4 | 2,1·10⁻¹³; \|z\| 0,94 / 2,63 / 1,85 / 3,28 (flusso / Bell flusso / Bell / balistico); vecchio **37 282** | PASS ×5 |
| E5 | atomo 3-D 1s ascoltato, 1000 tick: equivarianza 6 schemi | \|z\| ≤ 4, TV ≤ 1,2 | \|z\| 2,2–3,0, TV 0,97–1,04; vecchio 15 933 | PASS ×12 |
| E5 | "quasi fermi": crescita della deriva radiale 100→1000 tick (flusso 2 tick) | ≤ 1,5 | **1,04** (flusso 1 tick 2,45; Bell 1,7; CDF 2,1) | PASS |
| E5 | partenza sbagliata (\|ψ\|⁴): H̄ finale/iniziale (flusso 2 tick) | ≥ 0,9 | **0,979** | PASS |
| E5 | v_r RMS / Zα (flusso 2 tick) | ≤ 0,25 | 0,057 | PASS |
| E5 | **v_φ /(Zα sinθ) (flusso 2 tick)** | 0,75–1,25 | **0,21** | **FAIL** |
| E6 | 3-D 1s+2p, 1670 tick (10 battimenti) | \|z\| ≤ 4, TV ≤ 1,2 | \|z\| 2,2 / 2,4 / 2,2; TV 1,01 / 1,01 / 1,00; vecchio **60 327**, con rinascite 28 949 | PASS ×6 |
| E7 | 3-D onda libera, 300 tick | \|z\| ≤ 4, TV ≤ 1,2 | \|z\| 2,5 / 2,2 / 2,6; TV 1,00 / 1,01 / 1,00 (emend. A3; riferimento a t = 0: 1,70) ; vecchio 24 218 | PASS ×6 |
| E8 | radiale: partenza sbagliata, H̄ finale/iniziale (flusso) | ≥ 0,999 | **1,000** | PASS |
| E9 | orbitali analitici (come `orbitali.rs`), 10 battimenti, S5 | \|z\| ≤ 4 | 3,86 (vecchio 7,3; con rinascite 2,8) | PASS |

Nota sul lettore vecchio: nella vista dall'alto (figura, pannello in alto al centro) i suoi surfisti dell'1s
**cadono verso il nucleo** e ci restano aggrovigliati: la velocità letta dallo spinore nel punto, a θ = 0,5, ha una
componente radiale finta verso il centro. Per questo z cresce a 10⁴.

## 5. Mescolamento finto (punto c del compito)

Prova: si fanno partire i surfisti **sbagliati** (da |ψ|⁴, più concentrati) in stati dove de Broglie–Bohm nel continuo
**non mescola**, e si guarda se la folla "si aggiusta da sola" (H̄ relativo; 1 = nessun mescolamento, 0 = equilibrio).

| schema | radiale 1s+2s, 20 battimenti (1-D: l'ordine è congelato) | 1s 3-D, 1000 tick (r si conserva) | salti casuali per surfista per tick (1s 3-D) | bit casuali per tick |
|---|---|---|---|---|
| flusso (1 tick) | **1,000** | 0,694 | 0 | 0 |
| **flusso (2 tick)** | — (in 1-D coincide col flusso a 1 tick) | **0,979** | **0** | **0** |
| CDF per asse (R3) | — | 0,811 | 0 | 0 |
| Bell in tempo di flusso | 0,879 | 0,148 (τ ≈ 124 tick) | 0,81 | — (tempi continui) |
| Bell minimo per asse | 0,866 | 0,148 (τ ≈ 150 tick) | 0,82 | 1,87 |
| Bell balistico | 0,079 | 0,140 (τ ≈ 24 tick) | 3,00 | 2,52 |
| de Broglie–Bohm continuo | 1 | 1 (r conservato) | 0 | 0 |

In parole:
* I **Bell** (tutte e tre le versioni) sono esatti in media ma **mescolano**: la folla sbagliata va all'equilibrio in
  poche centinaia di tick anche dove Bohm la lascerebbe sbagliata per sempre. Il mescolamento finto non dipende
  dal Bell "minimo" o "balistico": dipende dal fatto che c'è un dado.
* Il **flusso a 1 tick** in 3-D mescola un po' (H̄ al 69 %): perché il tick pari (x,y,z) e il dispari (z,y,x) fanno
  "respirare" la densità di una cella avanti e indietro; i surfisti seguono quel respiro e le due mappe non si
  annullano. È il residuo già visto in R3.
* Il **flusso a 2 tick** (il periodo vero della regola) **non mescola** (98 %): sull'1s il flusso di una coppia di
  tick è a divergenza nulla (§8, idea 2: la parte "gradiente" pesa 2,6·10⁻⁸), quindi i surfisti girano su linee
  chiuse come in Bohm. Il raggio quadratico medio dei loro spostamenti non cresce (1,04 fra 100 e 1000 tick: niente
  diffusione).
* **Tempo di correlazione** (H̄ scende a 1/e): Bell balistico ~24 tick, Bell minimo ~150, Bell in tempo di flusso
  ~124; flusso (1 e 2 tick) e CDF: mai entro i 1000 tick.

## 6. Fallimenti e sorprese (onestamente)

1. **E5 velocità di rotazione dell'1s 3-D: FAIL.** Bohm–Dirac esatto: v_φ = Zα c sinθ (identità: 2x/(1+x²) = Zα con
   x = (1−γ)/Zα). Sul reticolo del laboratorio: 0,21 (flusso), 0,09 (lettore vecchio). Diagnosi (dopo la corsa,
   dichiarata): l'1s *del reticolo* a θ = 0,5 ha ⟨Σz⟩ = 0,51 invece di ~0,99, e la sua stessa corrente dà 0,23;
   la formula di Bohm–Dirac sullo stesso spinore dà 0,36. Il reticolo grosso "scioglie" metà dello spin
   (zitterbewegung a scala di cella, i due fogli). Lo stato radiale messo sul cubo (non autostato) dava 0,82. Il
   surfista è fedele alla regola; è la regola a θ = 0,5 che non è ancora il continuo. Da rifare con a più piccolo.
2. **TV relativo (E5–E7): errore mio di implementazione.** Il riferimento era calcolato a t = 0, mentre la
   pre-registrazione dice "campionatore ideale della stessa densità". Per l'onda libera, che si allarga, dava 1,70
   per *tutti* gli schemi esatti allo stesso modo (segno che era il riferimento). Ricalcolato alla densità finale
   (`surf_tvnull.py`, emendamento A3): 1,00–1,01 per tutti gli schemi esatti, 4,4–5,5 per il lettore vecchio. Riporto entrambi.
3. **La corsa radiale 1s+2s** è stata uccisa dal limite di 10 minuti al primo tentativo (non aveva checkpoint);
   rifatta da zero con checkpoint, stesso seme (emendamento A2).
4. **Il lettore vecchio è molto peggio di quanto si credeva.** Con l'onda vera della zig-zag, a θ = 0,5, z arriva a
   10⁴–6·10⁴: la folla del Laboratorio Atomo non è |ψ|² (le rinascite dello 0,2 % per tick la tengono "vicino" solo
   in apparenza). Sugli orbitali analitici lisci (`orbitali.rs`) invece va quasi bene (z 7,3; 2,8 con le rinascite):
   lì l'onda è liscia e non c'è lo zig-zag.
5. **Limite continuo di ordine 1, non 2.** In 1-D la distanza fra i surfisti del reticolo e quelli di Bohm–Dirac del
   continuo scende come a¹·¹ (0,251 → 0,111 → 0,052 λ_C per a = 0,2 → 0,1 → 0,05). È la regola (spostamento poi
   moneta, non palindroma in 1-D) a essere del primo ordine nel punto dove la guardiamo, non i surfisti, che
   seguono la densità del reticolo esattamente. I surfisti di Bell si allargano attorno a quelli deterministici
   come √a (2,56 → 1,91 → 1,32 λ_C): il processo di Bell tende a Bohm quando a → 0 (Vink 1993), ma lentamente.
6. **Orbitali (S5)**: |z| massimo 3,86, sotto la soglia 4 ma vicino. Su 41 controlli correlati non è un allarme, ma
   lo segnalo.

## 7. Debiti

**Chiusi** (debito #7 di CONTESTO_CLOUD):
* "il lettore va migliorato (equivarianza 1,26·10⁻²)": **chiuso**, 2,9·10⁻¹⁵ sullo stesso test.
* "non esiste una guida deterministica che conservi |ψ|² esattamente sul reticolo": **chiuso** (S4/S4-2; in 3-D
  Liouville per tick a 10⁻⁸–5·10⁻⁷ con differenze finite, riavvolgimento a 10⁻¹²).
* "il salto casuale crea mescolamento finto": **chiuso per costruzione** (zero salti); misurato: nessun mescolamento
  nei due stati di controllo con S4-2 (1-D: H̄ esattamente 1,000 ; 3-D: 0,98).

**Aperti / nuovi:**
1. La rotazione dell'1s sul reticolo grosso (0,21 Zα c sinθ): dipende dal reticolo (θ = 0,5), non dai surfisti.
   Serve a più piccolo o la regola agglomerata.
2. Scelte discrete (non manopole continue, ma scelte): flusso netto invece di lordo; periodo di 2 tick; campo RT0
   dentro la cella; ρ lineare nel tempo dentro il tick; per S5 la correzione di Poisson L2-minima. Altre scelte
   (per esempio un campo dentro la cella di ordine più alto) darebbero traiettorie un po' diverse con la stessa
   esattezza: **le traiettorie restano un'interpretazione**, come già detto in CONTESTO §6.4.
3. Vortici sugli spigoli: dove la corrente gira attorno a uno spigolo del reticolo un surfista vicino allo spigolo
   può fare molti giri in un tick (visto su stati casuali: fino a centinaia di attraversamenti). Esatto ma costoso
   in casi patologici; limite di sicurezza 10⁵ attraversamenti (mai raggiunto nelle corse vere).
4. Le patch in Rust non sono applicate (§9): servono un altro thread e il ponte.

**Parametri liberi: 0.**

## 8. Idee geometriche e visive

Ogni idea: cosa si vedrebbe, perché funziona, costo (Intel UHD, l'unica GPU rimasta: ~42 ms per passo fisico con
~33 800 palline), regole.

**1. Il fiume tra le celle è il teorema di Stokes del reticolo** (provata: è il cuore dello studio). La densità vive
nelle celle (volumi), la corrente sulle facce (superfici). "Quanto entra meno quanto esce = quanto cambia dentro" è
il teorema di Stokes discreto, ed è esatto perché il bordo di un bordo è vuoto. *A schermo:* surfisti che non si
"perdono" mai, nemmeno in 10⁴ tick. *Costo:* nullo in più (la corrente è già calcolata dalla regola). *Regole:* sì,
è la corrente della zig-zag stessa.

**2. Anelli di spin e respiro: colorare i surfisti con la decomposizione di Hodge** (provata, `surf_hodge.py`). Il
flusso di due tick si spezza in una parte "gradiente" (quella che cambia la densità: il respiro) e una parte
"rotazionale" (a divergenza zero: le correnti che girano, nell'1s la corrente di spin di Gordon). Numeri:
* 1s 3-D: la parte gradiente pesa **2,6·10⁻⁸** del totale: il moto è *tutto* rotazione (anelli chiusi);
* 1s+2p che batte: la parte gradiente pesa da 6·10⁻⁴ a **7,8·10⁻³**, e sale e scende col battimento;
* la parte rotazionale è per il 69 % azimutale attorno all'asse dello spin (85 % sui blocchi 2³: il resto sono
  piccoli giri alla scala della cella, i "due fogli").
*A schermo:* surfisti blu-ghiaccio quando girano (spin), dorati quando vengono spinti dal respiro; nell'1s si vede un
toro di anelli tutto blu; quando l'atomo assorbe luce e batte, onde dorate pulsano lungo z al ritmo del dipolo.
*Costo:* un Poisson per fotogramma su 48³ (gradiente coniugato con partenza dal fotogramma prima, ~10–20 iterazioni:
stimato 5–10 ms sulla Intel UHD) — oppure solo ogni 4 fotogrammi. *Regole:* sì: nessuna geometria aggiunta, è la
corrente della regola divisa in due parti; colore = numero calcolato, non disegnato.

**3. Surfisti agglomerati: l'esattezza vale a ogni livello di blocco** (provata, `surf_macro.py`). Sommando la
corrente delle facce fini sulle facce di un blocco 2³ si ottiene la corrente *esatta* del blocco (Stokes:
errore 3·10⁻¹⁷). Quindi i surfisti possono vivere al livello grosso con la stessa esattezza (|z| 2,6) e:
* fanno anelli **più tondi**: deriva radiale RMS 2,3 % del raggio contro 9,1 % sulle celle fini (400 tick, 1s);
* mescolano ancora meno: H̄ 0,996 contro 0,978;
* costano 8 volte meno lavoro di griglia.
*A schermo:* è il "livello di agglomerato" come impostazione grafica che Valerio voleva: da lontano surfisti su blocchi
grossi (anelli lisci), da vicino sulle celle fini (si vede il tremolio dello zig-zag). *Costo:* minore del lettore
attuale. *Regole:* sì, ed è proprio la filosofia dell'agglomeratore (blocca, controlla l'errore: qui l'errore di
equivarianza è zero a ogni livello).

**4. Il tasto "riavvolgi" diventa vero.** Flusso esatto + regola reversibile = le scie dei surfisti tornano indietro
esattamente (10⁻¹² celle in 20+20 tick). *A schermo:* premendo Indietro le palline ripercorrono le loro scie al
contrario, senza "rinascite" che le teletrasportano. *Costo:* zero. *Regole:* sì.

**5. L'orologio del tempo di flusso: vedere i vortici senza disegnare fili.** Il tempo di flusso τ = ∫ds/ρ corre veloce
dove la densità è piccola: vicino ai nodi (i vortici, i punti dove l'onda si annulla) i surfisti girano velocissimi.
*A schermo:* colorare la scia con la velocità di attraversamento delle facce: i filamenti dei vortici si accendono da
soli come anelli di palline che sfrecciano, senza nessuna linea disegnata. *Perché:* è la stessa grandezza che il
motore calcola per muovere i surfisti. *Costo:* zero. *Regole:* sì (i vortici sono quelli dell'onda).

**6. Reticolo diverso (bcc, a diagonali) per anelli più rotondi** — solo proposta. Il cubo semplice rende la
rotazione dell'1s quadrata a scala di cella (azimutale 69 %). Un reticolo cubico a corpo centrato con 4 direzioni
diagonali è più isotropo. *Ma* è una regola nuova (cambia lo spostamento): **non rispetta "una sola strada"** senza
che Valerio lo chieda; l'idea 3 (blocchi) ottiene parte dello stesso effetto restando sulla strada.

## 9. Patch proposte per altri file (solo proposta, non applicate)

### 9.1 `src/atomo/zz.rs` (atomo vivo): corrente esatta + surfisti a flusso

La funzione Rust qui sotto è stata **compilata e confrontata** con la versione numba su 20 000 surfisti in un tick
dell'atomo 48³: celle identiche e posizioni identiche bit per bit (differenza 0,0).

```rust
/// A surfer: its cell and where it is inside the cell (0..1 per axis; 0.5 = centre, as x_of()).
#[derive(Clone, Copy, Debug)]
pub struct Surfista { pub c: [usize; 3], pub xi: [f64; 3] }
impl Surfista {
    pub fn pos(&self) -> [f64; 3] { [0, 1, 2].map(|k| x_of(self.c[k]) + (self.xi[k] - 0.5) * A) }
}

/// log1p(x)/x (→ 1 at 0, +∞ at x ≤ −1)
#[inline] fn phi(x: f64) -> f64 { if x <= -1.0 { f64::INFINITY } else if x.abs() < 1e-8 { 1.0 - 0.5 * x + x * x / 3.0 } else { x.ln_1p() / x } }
/// expm1(x)/x (→ 1 at 0)
#[inline] fn ee(x: f64) -> f64 { if x.abs() < 1e-8 { 1.0 + 0.5 * x + x * x / 6.0 } else { x.exp_m1() / x } }

/// Move one surfer through one step, EXACTLY keeping |ψ|² (Pollock 1988 in flow time τ = ∫ds/ρ).
pub fn pollock(p: &mut Surfista, n: usize, rho0: &[f64], rho1: &[f64], fl: &[Vec<f64>; 3]) -> usize {
    let idx = |c: [usize; 3]| (c[0] * n + c[1]) * n + c[2];
    let mut s = 0.0f64;
    let mut cross = 0;
    for _ in 0..100_000 {
        let c = idx(p.c);
        let (r0, r1) = (rho0[c], rho1[c]);
        let d = r1 - r0;
        let rs = (r0 + s * d).max(1e-300);
        let tau_end = (1.0 - s) / rs * phi(d * (1.0 - s) / rs);
        let (mut tau, mut jb, mut side) = (tau_end, usize::MAX, 0i32);
        let mut ab = [[0.0f64; 2]; 3];
        for j in 0..3 {
            let mut m = p.c; m[j] = (m[j] + n - 1) % n;
            let (a, b) = (fl[j][idx(m)], fl[j][c]);
            ab[j] = [a, b];
            let (dj, x0) = (b - a, p.xi[j]);
            let q0 = a + dj * x0;
            if q0 > 0.0 && b > 0.0 {
                let t = (1.0 - x0) / q0 * phi(dj * (1.0 - x0) / q0);
                if t < tau { tau = t; jb = j; side = 1; }
            } else if q0 < 0.0 && a < 0.0 {
                let t = (-x0 / q0) * phi(-dj * x0 / q0);
                if t < tau { tau = t; jb = j; side = -1; }
            }
        }
        if !tau.is_finite() { tau = 0.0; }
        for j in 0..3 {
            let [a, b] = ab[j];
            let dj = b - a;
            let q0 = a + dj * p.xi[j];
            p.xi[j] = (p.xi[j] + q0 * tau * ee(dj * tau)).clamp(0.0, 1.0);
        }
        if jb == usize::MAX { break; }
        s = (s + rs * tau * ee(d * tau)).min(1.0);
        cross += 1;
        if side > 0 { p.c[jb] = (p.c[jb] + 1) % n; p.xi[jb] = 0.0; } else { p.c[jb] = (p.c[jb] + n - 1) % n; p.xi[jb] = 1.0; }
    }
    cross
}
```

Modifiche a `Atomo`:

```rust
// nuovi campi
pub flusso: [Vec<f64>; 3],   // flusso netto della coppia di tick, link (s, s+e_j)
pub rho_coppia: Vec<f64>,    // densità all'inizio della coppia (tick pari)

// in shift(j, verso), dopo aver riempito tmp (t[0]=⟨u+|L⟩, t[1]=⟨u−|L⟩, t[2]=⟨u+|R⟩, t[3]=⟨u−|R⟩):
let tmp = &self.tmp;
self.flusso[j].par_iter_mut().enumerate().for_each(|(s, f)| {
    let (t, tn) = (tmp[s], tmp[nb(s, j, 1, n)]);
    // "va avanti" in questo verso: avanti L u− e R u+; indietro L u+ e R u−
    let (piu, meno) = if verso > 0 { (t[1].norm2() + t[2].norm2(), tn[0].norm2() + tn[3].norm2()) }
                      else { (t[0].norm2() + t[3].norm2(), tn[1].norm2() + tn[2].norm2()) };
    *f += piu - meno;
});

// in avanti(): se ticks è pari, prima di tutto: rho_coppia = densità; flusso = 0.
//              dopo il secondo tick della coppia: muovi i surfisti
pub fn muovi_surfisti_esatti(&self, s: &mut [Surfista]) {
    let rho1: Vec<f64> = self.psi.par_iter().map(|p| p.iter().map(|c| c.norm2()).sum()).collect();
    s.par_iter_mut().for_each(|p| { pollock(p, N, &self.rho_coppia, &rho1, &self.flusso); });
}
// in indietro(): la stessa cosa al contrario: il flusso accumulato disfacendo i due tick è −F e la densità va da
// ρ(2m+2) a ρ(2m): i surfisti ripercorrono esattamente la loro strada.
```

In `zzscena.rs::passo`: con i surfisti accesi **togliere le 40 rinascite per tick** (non servono più e rompono la
reversibilità); tenerle solo dopo `collassa` / `imposta` (l'onda cambia di colpo e i surfisti vanno ricampionati
con `nascite`, che va adattata a restituire `Surfista {c, xi}`). Per il disegno a ogni tick (la coppia ne muove due
alla volta) si può interpolare linearmente la posizione fra due aggiornamenti, oppure usare la versione a 1 tick
(esatta anch'essa, ma con un po' di "respiro" pari/dispari, §5). Costo: 20 000 surfisti × ~0,2 attraversamenti di
faccia ogni 2 tick: trascurabile rispetto al tick (stimato < 1 ms per tick in Rust con rayon).

### 9.2 `src/atomo/orbitali.rs` (orbitali analitici): corrente proiettata (S5)

Lì non c'è la zig-zag ma ψ(t) analitico su 64³: per ogni passo, P0 = |ψ(t)|², P1 = |ψ(t+dt)|² sulle celle,
corrente di Dirac ψ†αψ a metà passo mediata sulle facce, poi la correzione minima F ← F − ∇φ con Lφ = div F −
(P0 − P1) (gradiente coniugato senza dipendenze nuove, partenza dal φ del passo prima), poi `pollock`. Provato in
Python (E9): |z| ≤ 3,86 su 10 battimenti. La correzione è piccola (RMS 2·10⁻⁶ per passo) e toglie l'errore del
lettore. Togliere le 25 rinascite per passo.

### 9.3 Richiesta di ponte proposta (NON creata)

```json
{"id": "surfisti-esatti-1", "azione": "test", "note": "dopo l'applicazione della patch 9.1: cargo test --bin atomo; aggiungere un test che fa 200 tick avanti e 200 indietro con 2000 surfisti e controlla posizione finale = iniziale entro 1e-9, e il chi2 dei surfisti contro |psi|^2 dopo 400 tick (|z| < 4)"}
{"id": "surfisti-esatti-foto", "azione": "foto", "note": "scena zig-zag 1s con surfisti accesi, 600 tick, prima e dopo la patch: anelli attorno all'asse dello spin"}
```

## 10. Dove l'intuizione di Valerio può aiutare

1. **Lo spin "sciolto" dell'1s grosso.** Sul reticolo a θ = 0,5 l'1s ha metà dello spin (⟨Σz⟩ = 0,51) e i surfisti
   girano a 1/5 della velocità di Bohm–Dirac. Domanda: preferisci vedere nel laboratorio la *verità del reticolo*
   (giro lento) o andare subito a reticoli più fini / alla regola agglomerata per vedere il giro giusto?
2. **Quale storia raccontano i surfisti?** Le traiettorie esatte non sono uniche (si può riempire la cella in altri
   modi con la stessa esattezza). Bohm sceglie "la corrente locale"; noi scegliamo "la corrente più semplice dentro
   la cella". C'è un principio fisico (minima azione? minima dissipazione? simmetria del cubo?) che secondo te
   dovrebbe sceglierne una?
3. **Il dado o il fiume?** Bell (salti) è esatto ma mescola; il flusso non mescola. Nella tua idea di "equilibrio
   quantistico come zucchero nel caffè", vuoi che il cucchiaino sia solo l'onda (flusso, nessun dado) o anche un
   piccolo rumore vero? Con il flusso esatto possiamo finalmente rifare il test di Valentini sulla regola a reticolo
   senza il mescolamento finto.
4. **I due fogli.** La rotazione rotazionale fine (31 % non azimutale) sono piccoli giri alla scala della cella che
   spariscono sui blocchi 2³. Sono "fisica" (un pezzo dello spin a scala λ_C) o solo il raddoppio dei fermioni? Una
   tua idea su come distinguerli si potrebbe tradurre in un test.

## 11. File e come rifare

* `lab/surf_core.py` — cammini con corrente (1-D, radiale, 3-D), mover esatto (numba), Bell, CDF per asse, lettore
  vecchio, statistica. `python lab/surf_core.py` (autotest, 2 s).
* `lab/surf_dev.py` — controlli di sviluppo (dichiarati).
* `lab/surf_states.py 1s|2p` — ascolto degli stati 3-D (driver `lab/surf_drive_states.sh`, ~20 min).
* `lab/surf_prereg.py` → `lab/out/surf_prereg.json`; `lab/surf_amend.py` → `lab/out/surf_amend.json`.
* `lab/surf_run.py f6|cont|h1s|beat|a3d free|a3d 1s|a3d 1s2p|orb` (driver `lab/surf_drive_all.sh`,
  `lab/surf_drive_beat.sh`) → `lab/out/surf_<corsa>.json`.
* `lab/surf_tvnull.py` (emendamento A3), `lab/surf_hodge.py`, `lab/surf_macro.py` (idee §8).
* `lab/surf_figure.py` → `lab/out/surf_figure.png`, `lab/out/surf_results.json`.
* Gli stati `.npy` (lab/out/surf_ckpt/) non sono nel repo: si rigenerano con `surf_states.py`.

![surfisti esatti](../lab/out/surf_figure.png)
