"""nucleo1s_f5_amend2.py -> lab/out/nucleo1s_f5_amend2.json : SECOND DECLARED AMENDMENT of phase 5, written after
the learned s0 at a = 0.5 for Z = 20, 40, 80 were seen (0.641, 0.618, 0.579; Z = 60: 0.555) and BEFORE computing the
new predictions below. Input from docs/ZSTAR.md (commit fc8a08c, other thread): near the nucleus the radial rule has no
scale (angles kappa/(i+s), Z alpha/(i+s)); the first cells fix the mixing beta of r^{+gamma} and r^{-gamma}; the 1s
error goes as beta a^{2 gamma}; sites at (i + sqrt3/2) a (the radius of the 3-D nearest neighbours) remove the fake dive.
New predictions (no learning, no fit):
 P_C1  scale-free limit of P_B1: the factor s on the Coulomb angle of the FIRST radial site (r = a/2) that makes the
       radial 1s of uniform cells equal to the resolved one, at a = 0.01 (where O(a) terms are gone and only the
       contact mixing beta is left: s cancels beta). Root rule: the largest root <= 1 (the bare lattice over-binds:
       the contact acts as an extra attraction, so the nucleus must pull less). a = 0.02 reported as a stability check.
 P_C2  same with radial sites at (i + sqrt3/2) a (first site at the 3-D nearest-neighbour radius).
Use with s1 = 1 (the 6 next sites unchanged).
Exams (new data, never computed before this commit):
 S5: P_C1 (primary) s0 at Z = 20, 40, 60, 80 within 0.05 of the learned s0 at a = 0.5 (ALREADY SEEN: declared as a
     comparison, not blind).
 S6 (blind): with s = (P_C1(Z), 1) at Z = 60, a = 0.3 and 0.2, and Z = 40, 80 at a = 0.3: valleys closed (spread of
     the 1s lines, weight >= 0.05, <= 0.002) at every point; at Z = 60 also |1s - Dirac| <= 1 %.
 The 3-D 9x9 junction (point 2) is lower priority (coordinator, 13/10): it is finished only if budget remains."""
import json, os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
json.dump(dict(amends="nucleo1s_f5_prereg.json (89e99f2), nucleo1s_f5_amend1.json (1d8467f)",
               seen_before=dict(learned_s0_a05={"20": 0.641461, "40": 0.617755, "60": 0.554593, "80": 0.578939}),
               predictions=dict(PC1=dict(a=0.01, offset=0.5, root="largest <= 1"),
                                PC2=dict(a=0.01, offset=0.8660254037844386, root="largest <= 1")),
               exams=dict(S5=dict(tol=0.05), S6=dict(points=[[60, 0.3], [60, 0.2], [40, 0.3], [80, 0.3]], spread=0.002,
                                                     tol_Z60=0.01))),
          open(os.path.join(OUT, "nucleo1s_f5_amend2.json"), "w"), indent=1)
print("written")
