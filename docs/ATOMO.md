# Laboratorio Atomo: un atomo di carbonio dalle palline

(Documento in costruzione: questa prima versione contiene solo la **pre-registrazione**, scritta e
committata prima della corsa di produzione `atomo --misura tutto`.)

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
