"""motore_f2_prereg.py - PRE-REGISTRATION of phase 2 of the "motore ottimizzato" study (written and committed before
the production runs motore_f2_run.py). Development done before this file (declared):
 (i)  octant solver at the Intel level (a0/8, octant 96^3, box +-12 a0): 1s -13.54753, 2p -3.39145, 2s -3.43550 eV,
      against the full-box values of AGGLOMERATO_ZIGZAG (same level, same box, unsymmetrized stencil): -13.54776,
      -3.39159, -3.43552 -> the cube symmetrization changes the levels by <= 4e-5 (relative);
 (ii) the local cosine transform: sign convention (+1 gives fast decay), bell = sine iterated twice, blocks of 32 fine
      cells (16 coarse = the block of reticolo_core), half-overlap 16 fine cells: high-frequency leakage of a smooth
      packet at 24 fine cells per lambda 1.3e-8 (1-D).
Run: python lab/motore_f2_prereg.py -> lab/out/motore_f2_prereg.json"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import motore_f2_core as f2

Q = 0.5
p2p = f2.exact_kick_2p(Q)
p1s = (16 / (4 + Q ** 2) ** 2) ** 2
RY = 0.5 * f2.zz.ALPHA ** 2 * f2.MC2
prereg = dict(
    study="motore_ottimizzato fase 2", date="2026-10-13",
    tiers=dict(ultra=dict(cell="a0/16 (chain A, B = 256)", octant=192, tol=0.02),
               intel=dict(cell="a0/8 (chain A, B = 512)", octant=96, tol=0.10)),
    scene=("atomo di idrogeno (regola agglomerata imparata, nucleo su uno spigolo, scatola +-12 a0) nello stato 1s, colpito "
           "da un impulso brevissimo lungo z (calcio: psi -> e^{i q z} psi, q = 0.5/a0, il limite di un impulso molto piu' "
           "corto del periodo atomico). Poi l'atomo 'batte': il dipolo <z>(t) oscilla alla riga 1s->2p."),
    exams=dict(
        G1a="riga 1s->2p (E2p - E1s letta dalla regola) entro la tolleranza del livello contro Schroedinger 10.2043 eV",
        G1b="probabilita' dopo il calcio: P(2p0) contro l'esatto %.6f, P(1s) contro l'esatto %.6f, entro la tolleranza" % (p2p, p1s),
        G1c="conservazione: norma dello stato colpito sotto la regola per M passi (100 Intel, 40 Ultra) <= 1e-10",
        G1d="sequenza: E1s < E2p e E1s < E2s; P(1s) > P(2p0) > P(2s); i primi 4 zeri del dipolo a due livelli entro tol*t dagli esatti",
        G1e="Intel contro Ultra: riga e P(2p0) entro il 10 % (la stessa dinamica a due risoluzioni)",
        G2=("cartolina spettrale: gli stati 1s, 2s, 2p0 sono fermi sotto la regola: residuo r = |W psi - lambda psi|; lo stato "
            "colpito proiettato e fatto girare con le fasi coincide con la griglia dopo M passi, nel sottospazio, entro "
            "M*r*2 + 1e-12; peso tagliato (fuori da 1s, 2s, 2p0) riportato; costo per passo griglia/cartolina riportato"),
        G3a=("cambio di livello LOCALE (trasformata a coseni locali, blocchi 32 fini / 16 grossi, sovrapposizione 16/8, campana "
             "seno iterato 2 volte) per la luce 2-D, pacchetti a 30 gradi come D1: |dE/E| alla restrizione e nel giro <= 1e-6 "
             "per N grosso >= 12; 8 riportato; confronto con cubico (D1)"),
        G3b="localita': una perturbazione dentro un blocco fine cambia il campo grosso solo entro il blocco +- la sovrapposizione; oltre: <= 1e-14",
        G3c="carica dell'elettrone (1-D, regola fine theta = a0/4096): 4 dimezzamenti locali (fine -> 16x), |dN/N| <= 1e-6 e giro <= 1e-6, energia cinetica <= 1e-6",
        G4=("descrittore a energia positiva al livello Intel 1-D (B = 192): (a) peso sulla banda negativa grossa di un pacchetto "
            "atomico positivo dopo il blocco (riportato); (b) ampiezza del battito finto (componente alla frequenza aliasata del "
            "tremolio) senza proiezione (riportata) e con proiezione <= 1e-10; (c) perdita fra bande per passo <= 1e-12; "
            "(d) dinamica lenta (velocita', larghezza) cambiata dalla proiezione <= 1e-6; (e) uno stato con peso negativo > 1e-3 "
            "(R puro fermo) viene SEGNALATO e resta al livello fine"),
        G5="gusci di livelli attorno all'occhio: calcolo del costo (celle e celle x tick con passo locale) contro il mondo pieno; nessuna soglia"),
    predictions=dict(
        G1a="Intel: 2p - 1s = -3.39145 + 13.54753 = 10.156 eV (-0.47 %) PASS; Ultra entro 0.3 % PASS",
        G1b="Intel ed Ultra entro qualche %: PASS", G2="residui 1e-7..1e-8: PASS", G3a="1e-8..1e-7: PASS per N >= 12",
        G4="PASS", exact_P2p=p2p, exact_P1s=p1s, line_exact_eV=RY * 0.75),
    configs=dict(q_a0=Q, box_a0=12, M_intel=100, M_ultra=40, lot=dict(L=32, e=16, bell_iter=2), seed=20261013))
json.dump(prereg, open(os.path.join(f2.OUT, "motore_f2_prereg.json"), "w"), indent=1)
print("P2p exact", p2p, "P1s exact", p1s, "line", RY * 0.75)
