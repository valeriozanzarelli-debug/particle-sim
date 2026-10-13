# Riflusso quantistico sulla regola zig-zag

(in corso: questa e' la parte di ricerca, scritta prima delle corse)

## 1. Il problema in due righe

Prendi una particella libera fatta **solo** di onde che vanno a destra (impulsi positivi). Ti aspetti
che la probabilita' di trovarla a sinistra di un punto possa solo calare. In meccanica quantistica no:
per un po' di tempo puo' **crescere**. E' il *riflusso quantistico* (quantum backflow, Allcock 1969;
Bracken e Melloy 1994). Il massimo di probabilita' che puo' tornare indietro in un tempo qualsiasi T e'
un numero puro, la costante di Bracken–Melloy **c_bm = 0,0384517** (Penz, Grübl, Kreidl, Wagner,
J. Phys. A 39, 423, 2006), noto solo numericamente: non c'e' una formula chiusa.

## 2. Cosa e' noto (ricerca, ottobre 2026)

| chi | cosa |
|---|---|
| Allcock, Ann. Phys. 53, 311 (1969) | prima osservazione: corrente negativa con impulsi positivi |
| Bracken & Melloy, J. Phys. A 27, 2197 (1994) | il massimo e' un autovalore di un operatore integrale, ~0,04, indipendente da T, m, ħ |
| Penz et al., J. Phys. A 39, 423 (2006) | c_bm = 0,0384517 |
| Eveson, Fewster, Verch, Ann. H. Poincaré 6, 1 (2005) | limiti al flusso istantaneo (illimitato dal basso nel continuo) |
| Bostelmann, Cadamuro, Fewster, PRA 96, 012112 (2017) | riflusso e scattering |
| Melloy & Bracken, Found. Phys. 28, 505 (1998) | elettrone di Dirac (energia e impulso positivi): il massimo dipende da ε = √(4ħ/(mc²T)) e cala da c_bm al crescere di ε |
| Su & Chen, Mod. Phys. Lett. A 33, 1850186 (2018) | Dirac spin 1/2 |
| Ashfaque, Lynch, Strange, Phys. Scr. 94, 125107 (2019), arXiv:1903.07618 | tabella c_rbf(ε): 0,03686 (ε=0,1), 0,03088 (0,5), 0,02722 (0,8), 0,02498 (1), 0,01947 (1,6), 0,01660 (2), 0,01372 (2,5); fit c_bm·exp[−(4ε/9)(1−4αε)] |
| Bialynicki-Birula et al., J. Phys. A 55, 255702 (2022) e Commento di Barbier, Fewster, Goussev, Morozov, Srivastava | riflusso nelle equazioni relativistiche; il commento distingue riflusso ondulatorio e riflusso quantistico |
| Paccoia, Panella, Roy, arXiv:2201.12916 (2022) | Dirac su un anello |
| Goussev, Quinque, Joo, Burbanks, PRA 110, 022216 (2024) | anello: limiti ottimali, curva frattale di dimensione 7/4 |
| Torres Arvizu, Ortega, Larralde, arXiv:2511.16867 (2025–26) | **reticolo tight-binding** (tempo continuo): massimo 0,0764734 ≈ 2 c_bm (a ν = τT/ħ ≈ 50); su anello 0,1313 |
| Eliezer, Zacharias, Bahabad, Optica 7, 72 (2020), e altri lavori ottici successivi | solo analoghi ottici; il riflusso della materia non e' mai stato misurato |

Aperto: (1) una forma chiusa di c_bm; (2) per Dirac la curva c_rbf(ε) e' solo numerica (tabella fino a
ε = 2,5); la legge per ε grande (regime ultrarelativistico) non l'ho trovata in letteratura (ho letto
solo abstract e la tabella di Ashfaque et al.: non e' una ricerca esaustiva); (3) nessuno studio trovato
sul riflusso in una **camminata quantistica a tempo discreto** (quantum walk di Dirac), che e' esattamente
la nostra regola.

## 3. Come si scrive il riflusso nella regola

Nella regola 1-D ogni tick sposta R di un sito a destra e L di uno a sinistra, poi ruota (R, L) di θ.
Quindi la probabilita' che attraversa il taglio tra i siti 0 e 1 in un tick e' **esatta**:
|R(0)|² va a destra, |L(1)|² va a sinistra. Il riflusso in N tick e' Δ_N = Σ (|L_n(1)|² − |R_n(0)|²).
Gli stati ammessi sono quelli della banda a energia positiva (cos ω = cos θ cos k) con k in una
finestra: k ∈ (0, π/2) e' la "finestra fisica" (dove la regola somiglia a Dirac), k ∈ (0, π) la banda
intera (contiene il "doppione" vicino a k = π, un artefatto del reticolo). Il massimo e' l'autovalore piu'
grande di una matrice, calcolato in `lab/rifl_core.py`; il nucleo e' stato controllato contro la regola
vera (`zz_core.Walk`) a 10⁻¹⁶.
