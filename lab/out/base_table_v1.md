| # | esperimento / grandezza | misurato | teoria | tolleranza | esito | fotoni / clic | palline-onda per fotone (livello 0 · tutti i livelli) | tempo |
|---|---|---|---|---|---|---|---|---|
| 1 | 2-D free beam width at z=1 zR | 4.2438 | 4.2426 | max(3σ=0.015, 0.02) | **PASS** | 20000 / 19975 | 3.969 · 1e+06 | 21 s |
| 2 | 2-D free beam width at z=2 zR | 6.7123 | 6.7082 | max(3σ=0.015, 0.02) | **PASS** | 20000 / 19975 (idem) | idem | 21 s |
| 3 | 2-D free beam width at z=5 zR | 15.3098 | 15.2971 | max(3σ=0.015, 0.02) | **PASS** | 20000 / 19975 (idem) | idem | 21 s |
| 4 | 3-D free beam width at z=1 zR | 4.3402 | 4.2426 | max(3σ=0.024, 0.02) | **PASS** | 4000 / 3987 | 9.743.750 · 2.1e+11 | 103 s |
| 5 | 3-D free beam width at z=2 zR | 6.8684 | 6.7082 | max(3σ=0.024, 0.02) | **FAIL** | 4000 / 3987 (idem) | idem | 103 s |
| 6 | 3-D free beam width at z=5 zR | 15.6709 | 15.2971 | max(3σ=0.024, 0.02) | **FAIL** | 4000 / 3987 (idem) | idem | 103 s |
| 7 | Gouy phase on axis, ζ(+zR) − ζ(−zR) (read by a photon ball on the axis) | 1.5722 rad | π/2 = 1.5708 (paraxial); exact RS 1.5707 | ±0.03 | **PASS** | 0 / 0 | 193.057.128 · 4.2e+12 | 3 s |
| 8 | Gouy phase on axis, ζ(+2zR) − ζ(−2zR) | 2.2080 rad | 2·atan(2) = 2.2143 | ±0.03 | **PASS** | 0 / 0 (idem) | idem | 3 s |
| 9 | Young double slit: correlation with exact flux | 0.9965 | 1 (ideal sampler 0.9973±0.0010) | >= 0.98 | **PASS** | 27000 / 5809 | 22.249 · 5.7e+06 | 34 s |
| 10 | Young double slit: bright-fringe population error rms / max | 4.0% / 7.1% | 0 (ideal sampler rms 3.2%) | rms <= 8% | **PASS** | 27000 / 5809 (idem) | idem | 34 s |
| 11 | Young double slit: absorbed fraction | 0.7849 | 0.7806 | ±0.02 | **PASS** | 27000 / 5809 (idem) | idem | 34 s |
| 12 | Taylor build-up: correlation vs number of clicks | N=10: 0.247 (ideal 0.487±0.139); N=50: 0.843 (ideal 0.783±0.072); N=200: 0.959 (ideal 0.929±0.024); N=1000: 0.989 (ideal 0.984±0.006); N=3000: 0.993 (ideal 0.995±0.002); N=5809: 0.997 (ideal 0.997±0.001) | ideal sampler | within 3σ of the ideal sampler at every N | **PASS** | 27000 / 5809 (idem) | idem | 34 s |
| 13 | 3-D long double slit (x-marginal) vs exact 2-D flux: correlation | 0.9907 | 1 (ideal sampler 0.9970±0.0010); 2-D rule: see Young | >= 0.98 | **PASS** | 6000 / 5551 | 3.557.772 · 7.8e+10 | 94 s |
| 14 | 3-D long double slit: bright-fringe error rms / max | 6.5% / 9.9% | 0 (ideal rms 3.2%) | rms <= 8% | **PASS** | 6000 / 5551 (idem) | idem | 94 s |
| 15 | 3-D long double slit: absorbed fraction | 0.7397 | 0.7806 (exact 2-D RS) | ±0.02 | **FAIL** | 6000 / 5551 (idem) | idem | 94 s |
| 16 | single slit a=4.0: minimum m=1 (mean of ±) | 103.91 ± 0.42 | 103.28 = L·tan(asin(mλ/a)) [paraxial mλL/a = 100.0; exact RS 103.28] | max(3σ, 2%) | **PASS** | 60000 / 56852 | 256 · 6.5e+04 | 83 s |
| 17 | single slit a=4.0: minimum m=2 (mean of ±) | 230.34 ± 1.47 | 230.94 = L·tan(asin(mλ/a)) [paraxial mλL/a = 200.0; exact RS 230.95] | max(3σ, 2%) | **PASS** | 60000 / 56852 (idem) | idem | 83 s |
| 18 | grating N=5, d=8.0: principal maximum m=1 (centroid, mean of ±) | 25.192 ± 0.010 | 25.198 = L·tan(asin(mλ/d)) [paraxial 25.00] | max(3σ, 1%) | **PASS** | 40000 / 36324 | 5.440 · 1.4e+06 | 53 s |
| 19 | grating N=5, d=8.0: principal maximum m=2 (centroid, mean of ±) | 51.616 ± 0.018 | 51.640 = L·tan(asin(mλ/d)) [paraxial 50.00] | max(3σ, 1%) | **PASS** | 40000 / 36324 (idem) | idem | 53 s |
| 20 | grating N=5: secondary maxima between order 0 and ±1 | 3 (+), 3 (−) | 3 | exact | **PASS** | 40000 / 36324 (idem) | idem | 53 s |
| 21 | Airy disk D=6.0: first zero radius | 41.76 ± 0.23 | 41.52 = L·tan(asin(1.22λ/D)) [paraxial 1.22λL/D = 40.67; exact RS 41.52] | max(3σ, 2%) | **PASS** | 9952 / 9459 | 219.252 · 4.8e+09 | 80 s |
| 22 | Arago spot (disk a=3.0, beam w=5.0, L=12.0): fraction of clicks within r < 0.77 | 0.0143 ± 0.0013 | 0.0293 (exact RS) | max(3σ, 10%) | **FAIL** | 9957 / 9489 | 78.418.960 · 1.7e+12 | 177 s |
| 23 | Arago spot: density in the spot / density in the shadow ring | 1.6 | 4.1 (exact RS) | > 3.0 | **FAIL** | 9957 / 9489 (idem) | idem | 177 s |
| 24 | Talbot self-image: z of best image (clicks) | 125.4 | 128.0 = 2d²/λ  [exact RS best plane 125.4] | ±2.6 | **PASS** | 20000 / 18117 | 181.440 · 4.6e+07 | 114 s |
| 25 | Talbot half-distance: the image is shifted by d/2 (contrast minimum near zT/2) | 61.4 | 64.0 | ±2.6 | **PASS** | 20000 / 18117 (idem) | idem | 114 s |
| 26 | Fresnel edge: I(0)/I0 at the geometric shadow edge | 0.2628 ± 0.0110 | 0.25 | max(3σ, 0.02) | **PASS** | 60000 / 59504 | 44.100 · 1.1e+07 | 125 s |
| 27 | Fresnel edge: first maximum position | 8.550 ± 0.232 | 8.748 (exact RS, this beam) [plane-wave Fresnel u=1.2172: 8.607] | max(3σ, 2%) | **PASS** | 60000 / 59504 (idem) | idem | 125 s |
| 28 | Fresnel edge: first maximum height I/I0 (window ±0.1·sqrt(λL/2)) | 1.4544 ± 0.0163 | 1.4465 (exact RS, this beam) [plane-wave Fresnel 1.3615] | max(3σ, 0.02) | **PASS** | 60000 / 59504 (idem) | idem | 125 s |
| 29 | Fresnel edge: first minimum position | 13.403 ± 0.098 | 13.347 (exact RS, this beam) [plane-wave Fresnel u=1.8725: 13.241] | max(3σ, 2%) | **PASS** | 60000 / 59504 (idem) | idem | 125 s |
| 30 | Fresnel edge: first minimum height I/I0 (window ±0.1·sqrt(λL/2)) | 0.6702 ± 0.0139 | 0.6893 (exact RS, this beam) [plane-wave Fresnel 0.7920] | max(3σ, 0.02) | **PASS** | 60000 / 59504 (idem) | idem | 125 s |
| 31 | Fresnel edge: second maximum position | 16.794 ± 0.128 | 16.720 (exact RS, this beam) [plane-wave Fresnel u=2.3445: 16.578] | max(3σ, 2%) | **PASS** | 60000 / 59504 (idem) | idem | 125 s |
| 32 | Fresnel edge: second maximum height I/I0 (window ±0.1·sqrt(λL/2)) | 1.4002 ± 0.0279 | 1.3755 (exact RS, this beam) [plane-wave Fresnel 1.1822] | max(3σ, 0.02) | **PASS** | 60000 / 59504 (idem) | idem | 125 s |
| 33 | which-path, parallel polarizers: fringe visibility | 0.924 ± 0.014 | 0.907 | > 0.85 | **PASS** | 60000 / 6513 | 21.801 · 5.6e+06 | 80 s |
| 34 | which-path, orthogonal polarizers: fringe visibility | 0.000 ± 0.045 | 0.029 | < 0.15 | **PASS** | 60000 / 6513 (idem) | idem | 80 s |
| 35 | Malus (absorbing polarizer, 3-D) θ=0°: P(pass) | 1.0000 | 1.0000 = cos²θ | max(3σ=0.0020, 0.01) | **PASS** | 7500 / 3727 | 4.757.193 · 1e+11 | 51 s |
| 36 | Malus (absorbing polarizer, 3-D) θ=22.5°: P(pass) | 0.8554 | 0.8536 = cos²θ | max(3σ=0.0274, 0.01) | **PASS** | 7500 / 3727 (idem) | idem | 51 s |
| 37 | Malus (absorbing polarizer, 3-D) θ=45°: P(pass) | 0.4900 | 0.5000 = cos²θ | max(3σ=0.0388, 0.01) | **PASS** | 7500 / 3727 (idem) | idem | 51 s |
| 38 | Malus (absorbing polarizer, 3-D) θ=67.5°: P(pass) | 0.1475 | 0.1464 = cos²θ | max(3σ=0.0275, 0.01) | **PASS** | 7500 / 3727 (idem) | idem | 51 s |
| 39 | Malus (absorbing polarizer, 3-D) θ=90°: P(pass) | 0.0000 | 0.0000 = cos²θ | max(3σ=0.0020, 0.01) | **PASS** | 7500 / 3727 (idem) | idem | 51 s |
| 40 | three polarizers 0/45/90°, light polarized at 0°: P(pass) | 0.2490 | 0.2500 | max(3σ=0.0206, 0.01) | **PASS** | 8000 / 1482 | 22.298.321 · 4.9e+11 | 126 s |
| 41 | two crossed polarizers 0/90° (the middle one removed): P(pass) | 0.0000 | 0 | 0 exactly | **PASS** | 8000 / 1482 (idem) | idem | 126 s |
| 42 | three polarizers 0/45/90°, unpolarized light: P(pass) | 0.1239 | 0.1250 | max(3σ=0.0157, 0.01) | **PASS** | 8000 / 1482 (idem) | idem | 126 s |
| 43 | Mach–Zehnder φ=0°: P(D) | 0.9967 | 1.0000 | max(3σ=0.0020, 0.02) | **PASS** | 3000 / 2983 | 157.564 · 4e+07 | 262 s |
| 44 | Mach–Zehnder φ=90°: P(D) | 0.4993 | 0.5000 | max(3σ=0.0388, 0.02) | **PASS** | 3000 / 2983 (idem) | idem | 262 s |
| 45 | Mach–Zehnder φ=180°: P(D) | 0.0000 | 0.0000 | max(3σ=0.0020, 0.02) | **PASS** | 3000 / 2980 | 157.564 · 4e+07 | 241 s |
| 46 | Mach–Zehnder φ=270°: P(D) | 0.5043 | 0.5000 | max(3σ=0.0387, 0.02) | **PASS** | 3000 / 2980 (idem) | idem | 241 s |
| 47 | bomb tester: P(explosion) | 0.5088 | 0.5000 | max(3σ=0.0274, 0.02) | **PASS** | 3000 / 1471 | 151.909 · 3.9e+07 | 139 s |
| 48 | bomb tester: P(C) = bomb detected without interaction (dark port) | 0.2341 | 0.2500 | max(3σ=0.0237, 0.02) | **PASS** | 3000 / 1471 (idem) | idem | 139 s |
| 49 | bomb tester: P(D) | 0.2571 | 0.2500 | max(3σ=0.0237, 0.02) | **PASS** | 3000 / 1471 (idem) | idem | 139 s |
| 50 | same interferometer without bomb: P(C) (dark port) | 0.0000 | 0.0000 | max(3σ=0.0030, 0.02) | **PASS** | 3000 / 1471 (idem) | idem | 139 s |
| 51 | Wheeler, BS2 inserted after the photon passed BS1: P(D) | 0.9946 | 1.0000 | max(3σ=0.0040, 0.02) | **PASS** | 1500 / 1490 | 157.564 · 4e+07 | 218 s |
| 52 | Wheeler, BS2 removed after the photon passed BS1: P(D) | 0.5053 | 0.5000 | max(3σ=0.0548, 0.02) | **PASS** | 1500 / 1490 (idem) | idem | 218 s |
| 53 | Wheeler, BS2 removed after the photon passed BS1: fraction of photons whose detector is the straight continuation of their arm (C⇔A, D⇔B) | 0.0000 | 0.0000 | max(3σ=0.0040, 0.02) | **PASS** | 1500 / 1490 (idem) | idem | 218 s |

* 2-D free beam width at z=1 zR: err +0.03%  n=19975
* 2-D free beam width at z=2 zR: err +0.06%  n=19975
* 2-D free beam width at z=5 zR: err +0.08%  n=19975
* 3-D free beam width at z=1 zR: err +2.30%  n=3987
* 3-D free beam width at z=2 zR: err +2.39%  n=3987
* 3-D free beam width at z=5 zR: err +2.44%  n=3987
* Gouy phase on axis, ζ(+zR) − ζ(−zR) (read by a photon ball on the axis): max |ζ − atan(z/zR)| over ±2zR: 0.0476 rad; vs exact RS 1.74e-03
* Young double slit: correlation with exact flux: REGOLA_FOTONE: 0.9949
* Young double slit: bright-fringe population error rms / max: REGOLA_FOTONE: 3.9% / 4.7%
* Young double slit: absorbed fraction: REGOLA_FOTONE: 0.789
* 3-D long double slit: bright-fringe error rms / max: per fringe [3, 4, -9, -1, 10]
* 3-D long double slit: absorbed fraction: wall 0.7187
* single slit a=4.0: minimum m=1 (mean of ±): err +0.61%; same estimator on ideal samples of the exact pattern: 103.47 ± 0.42
* single slit a=4.0: minimum m=2 (mean of ±): err -0.26%; same estimator on ideal samples of the exact pattern: 231.03 ± 1.59
* grating N=5, d=8.0: principal maximum m=1 (centroid, mean of ±): err -0.02%
* grating N=5, d=8.0: principal maximum m=2 (centroid, mean of ±): err -0.05%
* grating N=5: secondary maxima between order 0 and ±1: zeros every λL/(Nd) = 5.0; same counter on 20 ideal samples of the exact pattern: 3.00 (min 3)
* Airy disk D=6.0: first zero radius: err +0.56%; same estimator on ideal samples of the exact pattern 41.60 ± 0.22; fraction inside the first zero 0.879 (exact 0.916 within 2.2 r1)
* Arago spot (disk a=3.0, beam w=5.0, L=12.0): fraction of clicks within r < 0.77: 136 clicks in the spot; exact I(0)/I_free(0) = 0.444
* Talbot self-image: z of best image (clicks): image fraction in the slit windows at best plane 0.822 (open fraction a/d = 0.25; exact 0.825)
* Talbot half-distance: the image is shifted by d/2 (contrast minimum near zT/2): contrast there 0.000
* Fresnel edge: I(0)/I0 at the geometric shadow edge: (window |x|<0.25; ideal sampler of the exact pattern, normalised to its own clicks: 0.2566 ± 0.0113)
* Fresnel edge: first maximum position: ideal sampler of the exact pattern 8.725 ± 0.197
* Fresnel edge: first minimum position: ideal sampler of the exact pattern 13.347 ± 0.115
* Fresnel edge: second maximum position: ideal sampler of the exact pattern 16.784 ± 0.119
* which-path, parallel polarizers: fringe visibility: corr with exact 0.9935; 3273 clicks
* which-path, orthogonal polarizers: fringe visibility: corr with exact 0.9865; 3240 clicks
* Malus (absorbing polarizer, 3-D) θ=0°: P(pass): 1497/1497 photons
* Malus (absorbing polarizer, 3-D) θ=22.5°: P(pass): 1278/1494 photons
* Malus (absorbing polarizer, 3-D) θ=45°: P(pass): 732/1494 photons
* Malus (absorbing polarizer, 3-D) θ=67.5°: P(pass): 220/1492 photons
* Malus (absorbing polarizer, 3-D) θ=90°: P(pass): 0/1488 photons
* three polarizers 0/45/90°, light polarized at 0°: P(pass): 989/3972
* two crossed polarizers 0/90° (the middle one removed): P(pass): 0/1987
* three polarizers 0/45/90°, unpolarized light: P(pass): 493/3978 (H: 493/2008, V: 0/1970)
* Mach–Zehnder φ=0°: P(D): 1493/1498 (C: 0)
* Mach–Zehnder φ=90°: P(D): 748/1498 (C: 742)
* Mach–Zehnder φ=180°: P(D): 0/1492 (C: 1488)
* Mach–Zehnder φ=270°: P(D): 756/1499 (C: 736)
* bomb tester: P(explosion): 1524/2995 
* bomb tester: P(C) = bomb detected without interaction (dark port): 701/2995 
* bomb tester: P(D): 770/2995 
* same interferometer without bomb: P(C) (dark port): 0/999 
* Wheeler, BS2 inserted after the photon passed BS1: P(D): 743/747 (C: 0) — QM with BS2 present: cos²(0) = 1
* Wheeler, BS2 removed after the photon passed BS1: P(D): 379/750 (C: 368) — QM with BS2 absent: 1/2
* Wheeler, BS2 removed after the photon passed BS1: fraction of photons whose detector is the straight continuation of their arm (C⇔A, D⇔B): 0/747 — dBB prediction 0 (no crossing of trajectories); QM: no prediction
