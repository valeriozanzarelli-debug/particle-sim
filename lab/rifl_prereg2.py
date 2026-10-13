"""rifl_prereg2.py - seconda pre-registrazione (Q1b), scritta DOPO aver visto la prima corsa Q1 e PRIMA delle corse Q1b.
Genera lab/out/rifl_prereg2.json.

Cosa e' successo: la prima corsa Q1 (eps fino a 100) aveva il taglio in impulso troppo basso (P <= 1000 << 1/T per
eps >= 20, griglia con un pavimento di fase sbagliato per T piccoli): quei numeri sono SBAGLIATI (c*eps^2 -> ~0.30).
Con una griglia geometrica corretta, due prove di sviluppo: eps = 20 -> 6.772e-4, eps = 100 -> 3.1173e-5 (taglio 800/T).

Previsione analitica (ricavata dopo la prima corsa, prima di Q1b): per m c^2 T << hbar il nucleo di Dirac
K = (1/pi) sin((E-E')T/2)/(E-E') (L L' - R R') per p << 1/T diventa (T/2pi)(|L><L| - |R><R|), di rango 2.
Il suo autovalore positivo tende a (T/2pi) int_0^inf L(p)^2 dp = (T/2pi) (m/2), perche'
int_0^inf (1 - p/sqrt(p^2+m^2))/2 dp = m/2 (e int R^2 -> inf spegne il termine misto).
Quindi:   c_rbf  ->  m c^2 T / (4 pi hbar) = 1/(pi eps^2),   cioe'  c_rbf * eps^2 -> 1/pi = 0.318310.
"""
import json, os
pre = {
    "studio": "RIFLUSSO, domanda Q1b (regime ultrarelativistico di Dirac)",
    "data": "2026-10-13",
    "previsione": "lim_{eps->inf} c_rbf(eps) * eps^2 = 1/pi = 0.3183099 (c_rbf ~ m c^2 T/(4 pi hbar))",
    "configurazioni": {"eps": [2.5, 5, 10, 20, 50, 100, 200, 500, 1000, 3000, 10000],
                       "taglio": "P = X/T, X in {50,100,200,400,800,1600}, estrapolato in 1/P",
                       "griglia": "dirac_grid(geo=True), cosc = 4, controllo a cosc = 2"},
    "esami": [
        {"id": "Q1b-1", "cosa": "delta(eps) = 1 - pi c eps^2 per eps in {200,500,1000,3000,10000}",
         "soglia": "delta decresce in modo monotono e |delta(10000)| <= 1%"},
        {"id": "Q1b-2", "cosa": "la regola zig-zag a theta piccolo nello stesso regime (theta N = m T fisso piccolo)",
         "soglia": "regola(a->0) entro 1% del Dirac continuo a eps = 10 e eps = 5 (completa E3)"},
    ],
    "sviluppo_gia_visto": {"eps20": 6.772e-4, "eps100": 3.1173e-5},
    "parametri_liberi": 0,
}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "rifl_prereg2.json")
json.dump(pre, open(out, "w"), indent=1, ensure_ascii=False)
print("scritto", out)
