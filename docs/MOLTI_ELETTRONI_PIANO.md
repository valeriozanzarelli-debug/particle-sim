# Piano: atomi a più elettroni dalla regola zig-zag agglomerata

Prodotto il 12–13 ottobre 2026 da un panel: tre progetti indipendenti (palline in tempo immaginario, espansione a cluster, campo medio + correlazione), tre giudici avversari (fisico severo, esperto di calcolo, avvocato delle regole di PIANO.md), una sintesi. Punteggi totali: {"imaginary-time-balls": 66.5, "cluster-expansion": 71.5, "mean-field-then-correlation": 67}.

Il testo sotto è il brief per l'agente che implementa (in inglese). Il documento dei risultati sarà `docs/MOLTI_ELETTRONI.md`, in italiano. **Non è ancora stato eseguito nulla**: è un piano, e i numeri citati come 'development observations' vanno riverificati nella pre-registrazione.

---

# Implementation brief: many-electron atoms (He to Ne) from the agglomerated zig-zag

Language of the code comments and of `docs/MOLTI_ELETTRONI.md`: Italian, like the rest of the repo. This brief is in English for the implementing agent.

---

## 0. How the three designs were combined

- **Backbone: "cluster-expansion" (score 71.5).** It keeps:
  - exact few-electron solutions (ECG) of the rule's own operator as the agglomerator's CHECK level;
  - a cluster rule (pair and triple terms) learned on Z ≤ 6 and tested on held-out ions;
  - the node error measured exactly at 3 and 4 electrons;
  - the full learned dispersion fed in through the F(s) identity, which gives exam R6 (the alpha² content of the rule).
- **Taken from "imaginary-time-balls":**
  - an automatic, audited choice of the walker level (kernel positivity per level);
  - the three-way check "same rule, two time directions" (R4);
  - He 2³S with its exact node r1 = r2 (D2), and bosonic Li (D3);
  - exact Coulomb scaling E ∝ m_kin, so theta → 0 costs nothing;
  - a mechanical configuration generator (the Layzer complex, built from the rule's own n-shell degeneracy);
  - listing the development observations in the pre-registration.
- **Taken from "mean-field-then-correlation":**
  - honest chain-of-custody accounting;
  - mean field "listened" from random starts, to prove no orbitals are hidden;
  - the mean field treated as a candidate descriptor that the agglomerator formally rejects and keeps only as the node carrier;
  - H⁻ unbound in the mean field and bound only with correlation;
  - the 3-D learned-lattice check of the radial ("angular agglomeration") reduction (R7);
  - pre-registered expected ΔSCF failures;
  - the warning about the A2a error cancellation;
  - pair-rule learning from walker histograms, kept for the later molecule step.
- **Dropped:**
  - "agglomeration in time" as a claim (it is only the central-limit theorem);
  - the walker-density self-consistency loop (mixed-estimator bias, selection bias);
  - the same-spin cusp exam (ill-posed);
  - exams that pass by construction (E1c, E2b, E8 report-only), which are now recorded but not scored.

---

## 1. Goal and honest scope

**Goal.** Total energies and observables of He, Li, Be (then optionally B to F) and Ne, plus their cations, from the agglomerated zig-zag rule:
- no tabulated orbitals or exponents;
- no external chemistry code;
- every imposed piece labelled as a DEBT;
- pre-registered exams.

Molecules come after this, starting with H₂, which is nodeless and therefore exact.

**What the zig-zag actually gives to atoms.** State this verbatim in the document. At the non-relativistic level the agglomerated rule contributes:
1. a 2-component (Pauli) descriptor whose generator is bounded below. This is the permission for imaginary time; the fine Dirac walk does not have it;
2. one number, m_kin(theta), which goes to 1 as theta → 0 (tanθ/θ law);
3. spin as a conserved label (spin part 4e-11);
4. kernel positivity at β above a measured β_min;
5. the alpha² mass-velocity content of its band, with no Darwin term.

Everything else is quantum mechanics of identical fermions with the rule's operator, plus declared debts. Quantitatively the scalable part is fixed-node DMC with Slater-Jastrow guides. What is new is the traceability of each ingredient and the exact CHECK at N ≤ 4, not the method. Say so plainly to pre-empt "old DMC relabelled".

**Rule 7 (PIANO.md) status.**
- The learned many-electron rule is the cluster factor J. It sets VMC energies and DMC cost. It cannot move nodes.
- The nodes are imposed (DEBT-7). Their error is measured exactly only up to N = 4.
- So rule 7 is met in the weak sense: the cheap rule is learned and checked, and the accurate rule rests on imposed nodes.
- Before molecules beyond H₂, Valerio must decide whether that is enough.

---

## 2. Derivation chain

Tags: GIVEN, LEARNED, DERIVED, EXACT, COMPUTED, CHOSEN (automatic, by a pre-registered criterion), DEBT-n (imposed by hand; to be retired later).

**L0 GIVEN.** Fine zig-zag (`lab/zz_core.py`, `lab/agg_zz_core.py`).
- 4-component Dirac walk, flip θ = m c a/ħ, Coulomb phase −Zα/r per tick.
- Inputs: m_e (through θ), α, Z.
- Imaginary time is illegal here. e^{−βH_eff} with H_eff = i log W / Δt sinks toward the lattice bottom of the negative band (about −πħ/Δt, ZIGZAG.md §2) instead of stopping at the 1s. Note: "no bottom" in the design text is wrong; the bottom is the lattice limit.

**L1 LEARNED (exists).** Chain A, level 1 (B = 16, cell a0/256, coarse tick 0.5353 ħ/mc² = 2.85e-5 a.u.).
- Descriptor `pauli2`, positive band only.
- Hermitian generator with H(k) ≥ H(0) (neg_gen = 0).
- Spin part 4.3e-11.
- Levels a0/16 and a0/8 exist but are only valid as real-time rules (see L2).

**L1a CHECK + DEBT-8 (parity).** The exported stencil has a parity-odd piece:
- (H(k) − H(−k))/H = 3e-4 at kh = 0.01, about 3e-2 at kh = 1.
- The fine Bloch matrices are checked at ±k (gate G0).
- If they are even, the odd piece is a learning artifact and the stencil is symmetrised: H_j → Re(H_j + H_−j)/2.
- If not, the odd part is kept as a complex weight. Either way it is declared.
- All numbers below are from the symmetrised stencil.

**L1b LEARNED.** m_kin(θ) from the k → 0 curvature of the symmetrised dispersion.
- Chain A: 1.000384 (symmetrised ratio 0.999614 at k = 1/a0, isotropic).
- Chain B (θ = 0.2676): 1.024518.
- Fine-walk law tanθ/θ gives 1.000373 and 1.024574, so the agglomerator adds about 1e-5.
- A new chain at θ_A/2 gives a third, independent point.
- m_kin(0) comes from a fit in θ² with the θ⁴ term allowed.

**L1c LEARNED (new run, same algorithm and hyper-parameters, only the observable changes).** Re-learn level 1 with the Z = 10 1s as observable.
- Call: `auto_level(..., energy_mode="state", k_obs = 10/a0, eps_energy = 3e-5)`, at θ_A and θ_A/2.
- Reason: the existing level-1 band, checked only on an a0-sized state, carries just 24–48% of the Dirac k⁴ term at k = 10/a0 (axis/face/body after removing the m_kin part). This is a development observation.
- Outputs:
  - spherically averaged band ε̄_L(k) and an anisotropy table;
  - F(s) = E_{p∼N(0, s·I₃)}[ε̄_L(|p|)];
  - F4(s) = the same average of |p|⁴.
- Beyond the level-1 zone (|k| > π/h ≈ 804/a0) the fine walk's own positive band is spliced in (gate G4).

**L2 DERIVED (Wick rotation).** Rotating the level-1 Strang step to imaginary time (τ → −iβ) gives
φ ← e^{−βV/2} K_β ∗ (e^{−βV/2} φ), with K_β = IFFT[e^{−β(H(k) − H(0))}].
- H bounded below makes e^{−βH} a projector onto the ground state of each symmetry sector.
- The rest phase becomes a constant absorbed into E_ref (no mod-2π ambiguity).

**L2a LEARNED, with an honest limit (positivity).**
- The generator itself is not stoquastic: 42 off-diagonal taps have the wrong sign (largest +1580 Eh against a diagonal of 87800 Eh).
- K_β has 0.6–3% negative weight for β ≤ 2e-5 a.u.
- It becomes a probability (negative weight about 1e-14) only from β_min ≈ 5e-5 a.u. (about 1.75 coarse ticks), with ΣK = 1 − 1e-11 (balls conserved).
- So what is learned is a bounded-below band. Positivity is the Gaussianisation of a nearly quadratic band above β_min. Do not claim more.
- Walker level, chosen automatically: the coarsest level whose kernel has negative weight < 1e-10 over the working β range. That is level 1.
- a0/16 and a0/8 are rejected (2–30% negative weight at β = 1e-3…1e-2, because their learned band is wrong for kh > 0.5). Declared, already seen.

**L2b DERIVED + DEBT-9 (continuum kernel).**
- Production walkers use the continuum limit of K_β: Gaussian jitter with D = ħ/2m_kin = 0.499808 a.u. at θ_A, read from the export. No literal 1/2 anywhere in the code.
- That is the only learned number in the walkers, and it goes to 1/2 as θ → 0. This is the finding, not a defect to hide.
- The lattice-to-continuum step is checked for one electron (R4) and optionally for He (D0b).
- Branching e^{−β[(V(x) + V(y))/2 − E_ref]} is the Strang split; with importance sampling it becomes the local energy E_L.

**L3 EXACT.** The free step of N balls is the tensor product of one-ball steps, so the generator is Σ_i ε(p_i). Spin labels are conserved (L1).

**DEBT-1.** The electron–electron interaction is the instantaneous Coulomb phase +α/r_ij, the same law and constant as the nucleus phase. It is not derived from exchanged light balls. Missing: Breit and retardation, roughly 1e-5–1e-4 Eh in He and order 1e-2 Eh in Ne.

**DEBT-2.** Antisymmetry under exchange of same-spin labels (spin-statistics) is a postulate. Exam D3 shows it is not free: bosonic Li sits below −8 Eh.

**DEBT-3.** Continuum production uses point-nucleus Coulomb. The lattice rule uses the cell average, itself AGGLOMERATO debt 3. Nuclear mass is infinite; reduced-mass scaling (from allowed nuclear masses) is applied only for experiment exams; mass polarisation is neglected.

**DEBT-4.** Non-relativistic projection ε → k²/2m_kin for all energy exams except R5/R6. Darwin and spin-orbit are absent from the rule (AGGLOMERATO debt 1). Mixing mass-velocity without Darwin would shift Ne by about −0.6 Eh, so the full-rule numbers are reported separately.

**L4 EXACT few-body level (N ≤ 4, L = 0).** ECG with stochastic variational basis growth (SVM), using the rule's operator.
- The kinetic term enters through the identity ⟨φ_A|ε(p_i)|φ_B⟩ = S_AB · F(s_i), with s_i = 2[A(A+B)⁻¹B]_ii (per-component momentum variance).
- This is exact for centred ECGs, because the marginal momentum distribution is isotropic. Cubic anisotropy enters only at second order, through L = 4 admixture; estimate it and report.
- Coulomb integrals are analytic.
- Antisymmetry: sum over same-spin permutations only (N↑!·N↓! signed terms) with a primitive spin function of fixed M_S. That is 1 term for He ¹S, 2 for He ³S and Li, 4 for Be — not 24.
- Scope: L = 0 atoms and ions only. P states and molecules go to walkers.
- Certified as the rule's exact level by the chain R4 → (D0b) → D1.

**L5 COMPUTED + DEBT-5 (mean-field carrier).**
- Each orbital is the stationary state of the rule's non-relativistic operator in the field of the nucleus plus the other electrons: Hartree from Gauss shells, exchange from antisymmetry via Y^k with exact s/p angular factors. In standard language this is Hartree–Fock (spherically averaged ROHF for open shells), and the document must call it that.
- The solver starts from random band-limited fields, with no orbital forms anywhere.
- Formal CHECK against exact 2-electron He: error about 1.45%, against an energy budget of 3e-4, so REJECTED as a final level. It is kept only as the node carrier.
- Conflict: MOLECOLE.md says "niente Hartree–Fock". Declare it, and flag it as a decision for Valerio.
- The radial reduction is checked against the 3-D learned lattice (R7).
- A Hartree-only variant is reported at Li and Be, where node errors can be measured exactly.

**L6 CHOSEN + DEBT-10 (configurations).**
- Candidate generator (mechanical, the same for every atom): all CSFs with the reference's L, S and parity inside the n-shell that is degenerate in the rule's own spectrum (2s/2p degenerate to 0.044 eV in A2b). This gives:
  - 1 CSF: He, Li, N, O, F, Ne;
  - 2 CSFs: Be (2s² + 2p²), B, C.
- Correlating orbital: the improved virtual orbital, listened in the (N−1)-electron field.
- CI coefficient from the rule's own 2×2 CI: 0 free parameters.
- Keep rule: a candidate CSF stays only if a pilot fixed-node DMC lowers the energy by more than max(1 mEh, 3σ_pilot).
- Production runs use fresh seeds, independent of the pilot (this fixes the min-of-noisy selection bias).

**L7 LEARNED (cluster rule) + DEBT-6 (functional family).**
- J = Σ χ(r_iI) + Σ u_σσ′(r_ij) + Σ f(r_iI, r_jI, r_ij) [+ Σ g(r_ij, r_ik, r_jk)].
- Cusps are DERIVED from the operator (Kato): antiparallel m_kin/2, parallel m_kin/4. The electron–nucleus cusp is exact in the radial orbitals, so χ has zero slope at 0.
- Learning: force matching of ∇J to ∇ln|Ψ_ECG| − ∇ln|Φ_mf| on samples of |Ψ_ECG|².
- Training set: He-like Z = 1…6 (¹S), Z = 2…6 (³S); Li-like Z = 3…6; Be-like Z = 5, 6, added only after exam T1.
- Held out: He-like and Li-like Z = 7…10, and Be (Z = 4).
- Cluster order and functional family are chosen by leave-one-Z-out cross-validation inside the training set.

**DEBT-7.** Fixed node for two or more same-spin electrons.
- Measured exactly at N = 3 (Li) and N = 4 (Be).
- Not measurable for N > 4 on this laptop: 5-electron L = 0 ground states do not exist for B-like ions, and released node is hopeless for Ne (boson–fermion gap about 130 Eh).

**L8 READ-OUT.**
- τ → 0 by the pre-registered rule: a linear fit through 3 τ values; quadratic only if χ² gives p < 0.05, and declared.
- θ → 0 by the exact scaling E(m_kin) = m_kin · E(1) (non-relativistic). Report both the θ_A and θ → 0 values.
- Reduced-mass scaling μ/m_e for experiment exams only.

---

## 3. Algorithm

### 3.1 Kernel audit and permission (`lab/multi_e_kernel.py`, Python, 2 BLAS threads)

1. **Load and convert.** Load `assets/atomo/agglomerato_zz.json`.
   - For each level of chains A and B: H_j = ((G₀₀ + G₁₁)/2)/τ/α² in Eh.
   - Patterns are in the scratchpad dev scripts `kcheck*.py` and `judge2_kbeta.py`. Copy them into `lab/multi_e_dev/`.
2. **G0 parity.** Evaluate the fine 3-D two-tick Bloch matrices at ±k on 200 random k with |k| ≤ 0.5/cell. Do it by applying `fine3_tick_fast` to plane waves; `fine1_bloch` exists for 1-D. Then symmetrise or keep the odd part (L1a).
3. **K_β scan.** On 64³ and 128³ grids, for β ∈ {1e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 1e-3, 1e-2, 3e-2}, record:
   - negative weight, ΣK, second-moment D, excess kurtosis;
   - β_min;
   - the automatic walker-level choice.
4. **m_kin.** Curvature fit at |k| ∈ [1e-3, 4e-3] cells⁻¹ for chains A, B and the new A/2. Fit m_kin(θ) = 1 + aθ² + bθ⁴. Use two learning seeds per chain for the noise estimate.
5. **R1, 1-D.**
   - Fine 1-D walk in a soft Coulomb well: diagonalise W, form H_eff with eigenphases in (−π, π]/Δt, project e^{−βH_eff}, and show the energy sinks below −mc².
   - The learned 1-D Pauli level (`chain_1d`) converges to its LOBPCG bound state.
6. **Export** `lab/out/multi_e_kernel.json`: D(θ_A), m_kin(θ_A, θ_B, θ_A/2), m_kin(0) ± noise, β_min, chosen level, stencil hash.

### 3.2 Re-learned level and dispersion tables (`lab/multi_e_disp.py`)

1. Run `auto_level` from the fine rule with energy_mode "state", the Z = 10 1s observable, energy budget 3e-5, and otherwise the unchanged AGGLOMERATO §9 hyper-parameters. Do it at θ_A and θ_A/2. Plane-wave learning takes minutes. Report the level the machine chooses.
2. Symmetrise (if G0 passes) and compute:
   - ε̄_L(k) on a radial table up to the zone edge, fine-band splice beyond it;
   - the anisotropy table;
   - F(s) and F4(s) by 1-D radial quadrature over isotropic Gaussians.
3. **G4 gate.** For the ECG bases actually used, the part of F from beyond the level-1 zone must be < 1e-6 relative, or the splice must be continuous to 1e-8.
4. Export `lab/out/multi_e_disp.json`.

### 3.3 Two time directions on the lattice (`lab/multi_e_lattice.py`; numba or a small Rust helper)

1. Z = 16 ion on the symmetrised level-1 lattice, 192³ cells (box ±6 a_Z, cell a0/256). Reuse `cell_avg_inv_r` and the LOBPCG path from `lab/agg_zz_atom.py` with the level-1 stencil.
2. **Deterministic three-way check:**
   - (a) LOBPCG eigenvalue of G/τ + V;
   - (b) FFT power iteration of T_β = e^{−βV/2} K_β e^{−βV/2} at β ∈ {2e-4, 1e-4, 5e-5}, then −ln λ_max/β extrapolated β → 0 (quadratic in β);
   - (c) real-time −arg⟨ψ|W|ψ⟩/τ of the same lattice rule.
3. **Lattice walkers.**
   - Integer positions; moves sampled from an alias table of K_β (support ±6σ, about 1e5 entries); branching from the Strang weights.
   - Mixed estimator Σw (H_lat ψ_T)/Σw ψ_T with a NON-exact smooth guide ψ_T = e^{−Z m_kin r}, so the test is not tautological.
   - β ∈ {2e-4, 1e-4}.
4. Optional D0b: He 1¹S with lattice walkers (β ≥ β_min) against continuum Gaussian walkers at the same τ. Diagnostic only.

### 3.4 Exact few-body solver (`lab/multi_e_ecg.py`, numpy, 2 processes × 2 threads)

1. **Basis.** φ_k(r) = exp(−rᵀ(A_k ⊗ I₃)r), with r the electron positions relative to the nucleus.
   - A_k = Σ_{i≤j} α_ij w_ij w_ijᵀ (N(N+1)/2 nonlinear parameters).
   - Positive-definiteness is guarded.
2. **Antisymmetry.** Primitive spin function α…αβ…β with given M_S. Matrix elements are Σ_P sign(P) over P ∈ S_N↑ × S_N↓ (B → PᵀBP). The variational energy is an upper bound to the lowest state of that M_S sector with L = 0 and even parity, which is the ground state for He ¹S, He ³S (M_S = 1), Li ²S, Be ¹S and their ions. Optional: explicit spatial symmetrisation for singlets, for speed.
3. **Matrix elements** (C = A + B′):
   - S = (π^N/det C)^{3/2};
   - Coulomb ⟨1/|wᵀr|⟩ = S · (2/√π) · (wᵀC⁻¹w)^{−1/2};
   - kinetic, non-relativistic: S · (3/m_kin) · tr[A C⁻¹ B′];
   - kinetic, full rule: S · Σ_i F(s_i), with s_i = 2[A C⁻¹ B′]_ii.
4. **SVM.**
   - Grow one function at a time: 20–50 random trials, scaled by Z; keep the best by rank-1 secular update.
   - 3 refinement sweeps.
   - Conditioning guard: reject a candidate whose normalised overlap with the current span exceeds 1 − 1e-10.
   - Store E against K and extrapolate.
5. **Runs.**
   - Non-relativistic at m_kin(θ_A), then exact rescaling.
   - Full ε_L only where R5/R6 need it: He-like Z = 1, 2, 6…10.
6. **K budgets.** 150–200 for 2e (H⁻ needs 200), 400 for 3e, 600 for 4e.
7. **G1 dev validation (uses no exam system):**
   - 1-electron H, −0.5 to 1e-8;
   - Hooke's atom 2e at ω = 1/2, exact E = 2 Eh;
   - non-interacting Li-like and Be-like sums (−10.125 and −20 Eh at Z = 3 and 4), which test the antisymmetriser;
   - a Monte Carlo check of random 1/r_ij elements to 1e-4.

### 3.5 Mean-field carrier (`lab/multi_e_meanfield.py`)

1. Radial log grid (pattern of scratchpad `he_hf2.py`):
   - r_min = 1e-9/Z, r_max = 40–60 a0, 3000–12000 points, with doubling;
   - shift-invert eigensolver;
   - kinetic 1/(2m_kin);
   - exact point nucleus;
   - Hartree by Gauss shells, exchange by Y^k (k = 0, 1, 2) with exact s/p angular factors;
   - ROHF (spherical average) for open shells.
2. **Starts and occupations.**
   - Random band-limited starts; the fixed point must match the start from the rule's own hydrogen orbital to 1e-9 Eh (E1b in the old numbering).
   - Occupation enumerator over the lowest listened levels; the minimum is kept (aufbau emerges).
3. **Hartree-only variant** for Li and Be (report).
4. **Configurations.** IVO 2p in the (N−1) field; 2×2 CI from radial Slater integrals.
5. Export spline orbitals and CI to `lab/out/multi_e_mf.json`.

### 3.6 Cluster-rule learning (`lab/multi_e_learn.py`)

1. **Samples.** Metropolis sampling of |Ψ_ECG|²: 2e4 decorrelated samples per training state.
2. **Target.** t_i = ∇_i ln|Ψ_ECG| − ∇_i ln|Φ_mf|, with Φ_mf the same-symmetry determinant(s) for that ion. Exclude:
   - r_ij < 0.05/Z (ECG has no cusp there);
   - r_iI < 0.02/Z;
   - |Φ_mf| < 1e-3 × median (near mean-field nodes).
3. **Fit.** Linear least squares on the J coefficients, Huber loss, ridge 1e-8. Basis: polynomials in r̄ = x/(1 + x).
4. **Candidate families (DEBT-6):**
   - (A) u(r_ij) with coefficients c0 + c1/Z, no density scaling;
   - (B) u(r_ij / r_s(ρ_mf at the pair midpoint)).
5. **Candidate orders:** {u, χ}, {u, χ, f}, {u, χ, f, g}.
6. **Choice.** Smallest set within 1 percentage point of the best leave-one-Z-out VMC correlation fraction, using training Z only.
7. **Fallback.** If force matching gives a worse cross-validation score than matching smooth exact densities (intracule and electron–nucleus density, by iterative log-ratio update u ← u + ln(I_exact/I_VMC)), use the latter. Decided by cross-validation, pre-registered.
8. Exam T1 runs before Be-like Z = 5, 6 enter training.
9. Export `lab/out/multi_e_rule.json`.

### 3.7 Walker engine (Rust, new bin `elettroni`, `src/elettroni/`)

1. **Reuse.**
   - `Rng`, `drift_capped`, `media_errore` from `src/molecole/dmc.rs`;
   - determinant inverse (5×5) and fixed-node logic from `git show a84692d:src/atomo/carbonio.rs` (removed from HEAD in ce6adb7).
   - Do not edit `src/atomo` or `src/molecole`. Another agent's work tree is live; add only a `[[bin]]` entry to `Cargo.toml`.
2. **Guide.** Ψ_T = Σ_c C_c D↑_c D↓_c · e^J.
   - At most 4 determinants, at most 5×5 per spin.
   - Orbitals: radial splines × real Y_lm (l ≤ 1).
   - J from JSON.
   - D = 1/(2m_kin) from `multi_e_kernel.json`.
3. **Moves.**
   - Electron-by-electron, Sherman–Morrison updates;
   - UNR drift cap;
   - Reynolds accept/reject;
   - a move across a node of Ψ_T is rejected.
4. **Branching and estimators.**
   - Branching with the average of E_L before and after the move; population control T_pop = 1 a.u.
   - Mixed energy; blocking errors (largest of 32/16/8 blocks).
   - Spin-resolved intracule histograms, mixed and VMC, giving the extrapolated estimator 2I_mix − I_VMC.
   - Optional n(k) estimator.
5. **Modes.** VMC mode, DMC mode, unguided bosonic mode (D3), checkpoint/restart, `--misura <target>`.
6. **Tests (G2).** Local energies identical to the Python VMC at fixed configurations to 1e-10 for He and Li; hydrogen DMC −0.5.

### 3.8 3-D check of the radial reduction (`lab/multi_e_lattice3d.py`, R7)

1. He mean field on the learned a0/16 lattice and on a learned a0/32 level. Learn a0/32 from level 1 with the unchanged algorithm, about 10 min.
2. Isolated Poisson by zero-padded rfft: 512³ padding, about 1 GB per array, 4–6 GB peak.
3. Richardson extrapolation with the measured order. Unextrapolated values reported.
4. Optional: H₂ at R = 1.4 on a0/16 (the molecule path).

---

## 4. Pre-registration (`lab/multi_e_prereg.py` → `lab/out/multi_e_prereg.json`)

Commit it (`git add -f`) **before any production run**, including the R-class exams.

**It must contain:**
- every gate and exam below, with its expected outcome;
- K, τ triplets, walker counts, run lengths, grids and budgets;
- the configuration and keep rules;
- the J families and the training/held-out split;
- the τ-extrapolation rule;
- the pilot/production independence rule.

**Development observations already seen (list them as such):**
- level-1 negative weight about 1e-14 for β ≥ 5e-5, and 0.6–3% for β ≤ 2e-5;
- 42 wrong-sign taps;
- a0/16 and a0/8 negative weight 2–30%;
- parity-odd part 3e-4 at kh = 0.01, 3e-2 at kh = 1;
- symmetrised dispersion 0.999614 at k = 1/a0, and 0.99930 / 0.99903 / 0.99898 at k = 10/a0 (24–48% of the Dirac k⁴ term);
- D 0.499808 / 0.488034; m_kin(0) = 1.000001 from chains A and B;
- He HF −2.8616803 Eh in 8 s;
- the scratchpad scripts `kcheck.py`…`kcheck4.py`, `judge1_*.py`, `judge2_*.py`, `he_hf*.py`, `disp*.py`, copied into `lab/multi_e_dev/`.

**Pre-registered mean-field expectations (record, do not score):**
- E1c rejection, 1.45% on He;
- ΔSCF ionisation energies of about 23.45 / 5.34 / 8.04 / 19.85 eV (He / Li / Be / Ne), all expected to fail;
- H⁻ unbound in the mean field (−0.48793 > −0.5).

All energies are in Eh, non-relativistic and with infinite nuclear mass unless an exam says otherwise. σ is the statistical error after τ-extrapolation.

### Gates (must pass before production; a failure stops the work and is declared)

- **G0 Fine-rule parity.**
  - Target: Bloch matrices at k and −k equal to 1e-12.
  - Why: this justifies symmetrisation (DEBT-8); otherwise keep the odd part.
- **G1 ECG machinery.**
  - Target: H −0.5 to 1e-8; Hooke's atom 2.0 to 1e-8; non-interacting Li-like and Be-like to 1e-8; Monte Carlo 1/r_ij check to 1e-4.
  - Why: unit tests on non-exam systems, so no peeking.
- **G2 Rust vs Python.** Local energies equal to 1e-10.
- **G3 Mean-field implementation.**
  - Target with m_kin = 1: HF limits He −2.8616800, Li ROHF −7.4327269, Be −14.5730232, Ne −128.5470981, H⁻ −0.4879297, to 2e-5 relative. Random-start independence to 1e-9 Eh.
  - Why: literature HF values are used only as implementation unit tests, declared; they are never physics exams.
- **G4 Dispersion splice.** Beyond-zone fraction < 1e-6.
- **G5 Monotonic orbital ratio.** 2s/1s must be monotonic on the grid used for He 2³S, so the determinant's node is exactly r1 = r2.

### R-class: exams that touch the zig-zag rule itself

- **R1 Permission (1-D).**
  - Target: fine walk energy < −1 mc² by β = 50 ħ/mc² (it sinks toward the band bottom). Pauli level |E − E_LOBPCG| < 1e-6 relative.
  - Why: a qualitative contrast; 1e-6 is the solver precision.
- **R2 Walker-kernel audit.** Targets:
  - level 1: β_min ≤ 1e-4 a.u., negative weight < 1e-10 for all β ∈ [β_min, 1e-2], |ΣK − 1| < 1e-9;
  - D within 2e-5 of the curvature value;
  - |excess kurtosis| < 0.02 at β = 1e-2;
  - a0/16 and a0/8 recorded as non-walker levels (negative weight > 1e-3 at β ≤ 0.03; expected, already seen).

  Why: 1e-10 keeps sign effects far below any energy tolerance. The working τ (≥ 1e-3) is at least 20 β_min. The failure below β_min is disclosed, not hidden.
- **R3 m_kin(θ → 0).**
  - Target: 1, with tolerance 5e-5.
  - Why: learning noise is about 1e-5 per chain and the θ² extrapolation amplifies it about 1.7×; 5e-5 is about 3σ. The value is already roughly known (dev observation), so this is a confirmation with a fresh chain at θ_A/2.
- **R4 Same rule, two time directions (Z = 16, level 1).** Targets:
  - (a) power iteration (β → 0) = LOBPCG to 1e-6 relative; real-time readout = LOBPCG to 3e-5 relative;
  - (b) lattice walkers = power iteration at the same β within 3σ, σ ≤ 1e-4 relative.

  Why: (a) is deterministic, with the Trotter error at τ = 2.85e-5 a.u. below 3e-5. (b) σ is the about 20-minute budget; a non-exact guide avoids a tautology.
- **R5 ECG represents the rule's operator (one electron).** Targets:
  - non-relativistic Z = 1, 2, 10: E = −m_kin Z²/2 to 1e-6 relative;
  - full ε_L at Z = 10: (E_full − E_NR)/(−5α²Z⁴/8) = 1 ± 0.05.

  Why: higher-order (Zα)² effects are about 0.5%, and lattice θ² effects are 1e-3.
- **R6 α² content of the many-electron rule.**
  - Target: He-like Z = 6…10, (E_full − E_NR)/(−α²⟨Σp_i⁴⟩/8) = 1 ± 0.10, on the re-learned level.
  - Also reported on the existing level: expected FAIL (24–48%).
  - Why: ±10% is the A1c precedent. The Darwin term is absent by construction and is not tested.
- **R7 Angular agglomeration.**
  - Target: He mean field on learned 3-D lattices (a0/16 and a0/32, Richardson) = radial mean field to 1e-3 relative.
  - Why: a0/16 anisotropy at He momenta is about 1.5e-3 unextrapolated, so the Richardson value must reach 1e-3.

### M-class: many-body machinery on the rule's non-relativistic operator, θ → 0 by exact scaling

Without the θ → 0 conversion, M1 would miss by 1.1 mEh, so these exams also test the m_kin chain. They are not physics discoveries; label them as machinery.

| Exam | Target | Tolerance | Why |
|---|---|---|---|
| M1 He 1¹S (ECG) | −2.903724377 | ±0.29 mEh (1e-4 relative) | project standard; ECG converges to 1e-7 |
| M2 He 2³S (ECG) | −2.175229378 | ±0.22 mEh | same |
| M3 H⁻ (ECG) | −0.527751017; must also be below −0.5 | ±0.05 mEh | diffuse but cheap at K = 200 |
| M4 Li (ECG) | −7.478060323 | ±0.75 mEh | K = 400 expected to give 2e-5 to 5e-5 |
| M5 Be (ECG) | −14.6673565 | ±1.5 mEh | K = 600 expected to give 1e-4 to 5e-4 |
| Li⁺, Be⁺ (ECG, used for ionisation energies) | −7.279913413, −14.324763177 | 1e-4 relative | same reasoning |

If M5 fails to converge, the 4-electron CHECK is declared downgraded (D5 and C2 then use the ECG upper bound with its convergence uncertainty).

### D-class: projector and fermions

- **D1** He fixed-node DMC (nodeless, Wick-rotated rule, D from the export) = ECG He on the same operator. |Δ| ≤ max(3σ, 0.1 mEh), σ ≤ 1e-4. Why: it certifies that walkers and ECG represent the same rule.
- **D2** He 2³S DMC with the exact node r1 = r2 = ECG, within max(3σ, 0.1 mEh). Why: zero nodal error calibrates the fermion machinery.
- **D3** Bosonic Li (unguided, no node): E < −8.0 Eh. Why: the expected value is about −8.46 to −8.5 against −7.478, so the sign of the statistics is unmistakable.
- **D4** Li node error (chosen carrier, single determinant): −3σ ≤ E_DMC − E_ECG ≤ 0.3 mEh + 3σ, with σ ≤ 1e-4. Why: the expected error is about 0.03–0.1 mEh; the lower bound is the variational check.
- **D5** Be node error (automatically chosen node): −3σ ≤ ΔE ≤ 1.5 mEh + 3σ, with σ ≤ 2e-4. Single-determinant node error reported; expected about 8–12 mEh, which would fail. Why: the 2s/2p near-degeneracy is the textbook case where node choice matters.
- **D6** Kato cusp emerges in the walkers. He, guide with u = 0 (no electron–electron Jastrow); extrapolated intracule binned at 0.01 a0; fit h(r) = h0(1 + c·r + d·r²) on r < 0.25 a0. Target c = m_kin ± 0.10, opposite spin only. Why: about 7e6 pair samples fall below 0.25 a0, and the extrapolated-estimator bias is second order (it affects d, not c). The same-spin case (h ∼ r²) is not resolvable and is dropped.
- **D7** Variational sanity on every fixed-node energy: E_DMC ≥ E_exact − 3σ (for Ne, E_est − 1 mEh − 3σ). Why: a value below the exact energy means a bug.
- **D0b** (optional diagnostic) Lattice-walker He = continuum He within 3σ, σ ≤ 1e-3. The expected lattice correction is about 1e-5.

### S-class: selection emergence

- **S1.**
  - Li ground state 1s²2s, below 1s²2p by more than 3σ in DMC.
  - Mean-field enumeration gives 1s²2s²2p⁶ for Ne.
  - The configuration keep rule adds 2p² for Be and nothing for He, Li or Ne.
  - Why: all choices are made by a pre-registered internal criterion and are checked here.
- **S2** (optional, if B–F are run): Hund's rule emerges.
  - C ³P < ¹D, N ⁴S < ²D, O ³P < ¹D. The partner states are single-determinant M_L = 2 states.
  - Splittings 1.264 / 2.384 / 1.967 eV, each within ±0.15 eV.
  - Why: fixed-node errors of the two states differ by up to about 0.1 eV.

### C-class: the learned cluster rule (VMC)

Correlation fraction = (E_mf − E_VMC)/(E_mf − E_exact), with the same node in both.

- **C1** Held-out He-like Z = 7…10: ≥ 95%. Held-out Li-like Z = 7…10: ≥ 90%. Why: interpolation in 1/Z of a smooth family should be good; failing it means the cluster rule does not transfer even within its class.
- **T1** L–L transfer test, run before Be-like enters training. Be-like Z = 5, 6 with the J trained on 2e/3e data only: ≥ 70%. Why: it tests the hypothesis that the pair hole transfers by density, and it is genuinely falsifiable.
- **C2** Be (Z = 4, held out), final J, chosen node: ≥ 80%.
- **C3** Ne cheap rule (VMC, transferred J, no refit): ≥ 60% of E_mf − E_est (about 0.39 Eh). A real chance of failing is pre-registered as such.

### E-class: experiment and multi-electron exact values

Non-relativistic values with reduced-mass scaling; mass polarisation neglected.

| Exam | Target | Tolerance | Why |
|---|---|---|---|
| E1 He first ionisation energy | 24.587389 eV (NIST) | ±0.010 eV | prediction 24.58834; residual about +1 meV from relativity, QED and mass polarisation |
| E2 H⁻ electron affinity | 0.754195 eV (Lykke 1991) | ±0.005 eV | prediction about 0.75473 |
| E3 Li ionisation energy | 5.391715 eV | ECG ±0.010 eV; DMC route ±0.02 eV | ECG prediction 5.39138 |
| E4 Li 2s → 2p | 1.8478 eV (D-line centroid) | ±0.02 eV | both states fixed-node; σ_diff ≈ 4 meV |
| E5 Be ionisation energy | 9.322699 eV | ECG ±0.015 eV; DMC route ±0.05 eV | ECG prediction 9.32187 |

- **N1 Ne total** (DMC, chosen carrier, single determinant, θ → 0), σ ≤ 4 mEh.
  - (a) E ≥ −128.9376 − 1 mEh − 3σ. Sharp, variational.
  - (b) At least 90% of E_mf − E_est recovered, i.e. within about 39 mEh.
  - Why: the Chakravorty 1993 estimate is uncertain by ±0.5–1 mEh. The expected node error of 15–25 mEh comes from general fixed-node DMC experience and cannot be derived inside the project, so (b) is declared WEAK.
- **N2 Ne ionisation energy:** 21.564540 eV, ±0.25 eV, with σ_IE ≤ 0.08 eV.
  - Why: 3σ is about 0.24 eV, and node errors of Ne and Ne⁺ only partly cancel.
  - Declared as a ~1% test. Do not loosen it after the run; extend the checkpointed runs instead.
- **B1** (optional, only if added before those runs): ionisation energies of B, C, N, O and F (8.298019, 11.260288, 14.53413, 13.618055, 17.42282 eV), each ±0.15 eV, with the Be → B and N → O drops reproduced in sign.

**Report-only:**
- full-rule α² content for Ne, first order from the VMC n(k) folded with ε_L − ε_NR;
- Hartree-only carrier node errors for Li and Be;
- the A2a cancellation note (+0.034% from the mass against about −0.08% from the cell-averaged Coulomb at Z = 1, which will not cancel at Z = 10).

---

## 5. Cost on this laptop

Assumptions: i5, 4 cores, 16 GB, shared with the colla agent. Python jobs use 2 BLAS threads. Rust uses rayon on 4 threads when the machine is free.

Calibration from the repo: the old carbon run (2000 walkers × 44k ticks at τ = 0.005) took about 200 s per state, about 2.3 µs per walker-step wall, and gave σ ≈ 1.1 mEh. The old He run (10k walkers × 150 a.u.) took about 3 min for 3 τ values at σ ≈ 4.7e-4.

| Atom / step | Work | Wall time | RAM |
|---|---|---|---|
| M0 kernel, dispersion, R1–R4 | FFTs 64³–192³, 3 re-learnings, LOBPCG + power iteration on 192³, lattice walkers | 2–3 h | ≤ 3 GB |
| He | ECG: 19 non-relativistic states + 7 full-ε_L states, K = 150–200 | about 1 h | < 1 GB |
| He | DMC: He, He 2³S, H⁻ at σ ≤ 1e-4 (3 τ, 4000–10000 walkers) + D6 cusp run | 2–3 h | < 0.5 GB |
| Li | ECG Li-like Z = 3…10, K = 400, 2 terms per element | 3–5 h background, 2 processes | < 2 GB |
| Li | mean field (seconds), J learning on He + Li sets (about 1 h), VMC C1 (20 min) | about 1.5 h | |
| Li | DMC Li, Li 2²P, bosonic Li, Li⁺ check; σ ≤ 1e-4 for Li | 3–4 h | |
| Be | ECG Be-like Z = 4, 5, 6, K = 600, 4 terms per element | 4–8 h background | < 3 GB |
| Be | T1, J retrain, configuration pilots, DMC Be (single determinant and chosen node), Be⁺; σ ≤ 2e-4 | 5–6 h | |
| B–F (optional) | DMC atom + cation + Hund partner, 2000–4000 walkers | 2–3 h per atom, about 12–15 h | < 1 GB |
| Ne | VMC C3 (20 min); DMC Ne and Ne⁺ at τ 0.004 / 0.002 / 0.001, 2000–4000 walkers, plus a walker-count check | 6–10 h (1–2 nights, checkpointed) | < 1 GB |
| R7 | 3-D lattice He at a0/16 and a0/32 (+ H₂ optional) | 1.5 h | 4–6 GB |
| Write-up | figure, export, document | 1–2 h | |

How the Ne figure was obtained:
- Electron-by-electron moves cost about 7–12 µs per walker-step wall.
- About 4.4e5 walker-a.u. per τ gives about 2.5 mEh per τ.
- The extrapolated σ is about 4 mEh.

**Totals.** About 30–40 machine-hours for the required path, plus 12–15 h optional, over 5–6 sessions with two overnight runs. ECG runs go in the background alongside DMC.

For contrast, the grid route: (192³)^N = 5e13 (He) and about 1e68 (Ne).

---

## 6. Milestones, in order

1. **M0, pre-registration and the rule.**
   - Write `multi_e_prereg.py` and commit the JSON.
   - Then G0, R2, R3, R1, R4 (`multi_e_kernel.py`, `multi_e_lattice.py`) and the re-learned level plus F(s) (`multi_e_disp.py`, G4).
   - Implement the ECG and pass G1.
   - Implement the mean field and pass G3 and G5.
2. **M1, He.**
   - ECG He-like series (training and held-out), full-ε_L runs: M1, M2, M3, R5, R6.
   - Build the Rust bin `elettroni` and pass G2.
   - DMC: D1, D2, D6.
   - Exams E1, E2.
   - Record the mean-field CHECK rejection and the H⁻ result.
3. **M2, Li.**
   - ECG Li-like Z = 3…10: M4, Li⁺.
   - Mean field for Li (S1).
   - Learn J on He-like + Li-like Z ≤ 6 (cross-validation picks family and order). Exam C1.
   - DMC Li, Li 2²P, bosonic Li: D3, D4, E3, E4.
   - Hartree-only carrier report.
   - Optional: released-node diagnostic on Li (gap about 1 Eh, t ≤ 3 a.u.).
4. **M3, Be.**
   - ECG Be-like Z = 5, 6, then T1, then retrain J including Be-like.
   - ECG Be (M5) and Be⁺.
   - Configuration pilots (keep rule), then fresh production DMC: D5, C2, E5, S1 (Be part).
   - Optional M3b: learned backflow (about 6 counted parameters, VMC-optimised on Li-like Z ≤ 6 and Be-like Z = 5, 6, then frozen). Test on held-out Be: node error reduced by ≥ 30%. This is the only learned piece that moves nodes.
5. **M4 (optional), B, C, N, O, F.**
   - Mean field (spherical average declared for open p shells), configuration rule (B and C get 2 CSFs), DMC atom + cation + Hund partner: S2, B1.
   - Only if those exams were in the committed pre-registration before the runs; otherwise report-only.
6. **M5, Ne.**
   - Mean field; VMC with transferred J (C3).
   - Overnight DMC of Ne ¹S and Ne⁺ ²P (M_L = 0 hole), with a walker-count check at τ = 0.002: N1, N2.
   - n(k) α² report. If σ is above target, extend the runs rather than loosening tolerances.
7. **M6, R7 and write-up.**
   - R7 3-D check (optional H₂ on the lattice).
   - `docs/MOLTI_ELETTRONI.md` in Italian: PASS/FAIL table against the pre-registration, every FAIL, debts DEBT-1 to DEBT-10, parameter count, the "what the zig-zag gives" box, the HF decision, rerun instructions.
   - Figure `lab/out/multi_e_figure.png`.
   - Export `assets/atomo/molti_elettroni.json` (ECG A-matrices and coefficients for He/Li/Be for live sampling of |Ψ|² in the Laboratorio Atomo, orbitals, J, energies).
   - Update PIANO.md §3 and §4.
   - Commit with explicit paths.
8. **After this plan (not part of it): molecules.**
   - H₂ (R = 1.4, −1.174475714) as nodeless 2-electron walkers: an exact level below that ECG cannot reach (shifted Gaussians break the F(s) identity).
   - Learn pair rules for molecules from walker histograms (u ← u + ln(I_exact/I_T), extrapolated estimator).
   - Then LiH and H₂O, only after Valerio's decision on rule 7.

---

## 7. Files to create

**Python (`lab/`):**
- `multi_e_prereg.py` → `lab/out/multi_e_prereg.json`
- `multi_e_kernel.py` → `lab/out/multi_e_kernel.json` (G0, R1, R2, R3)
- `multi_e_disp.py` → `lab/out/multi_e_disp.json` (re-learned level, F(s), F4(s), G4)
- `multi_e_lattice.py` (R4, optional D0b)
- `multi_e_ecg.py` → `lab/out/multi_e_ecg.json` (G1, M-class, R5, R6; E against K)
- `multi_e_meanfield.py` → `lab/out/multi_e_mf.json` (G3, G5, S1 enumeration, CI, IVO, Hartree-only variant)
- `multi_e_learn.py` → `lab/out/multi_e_rule.json` (J, cross-validation, T1, C1)
- `multi_e_vmc.py` (Python reference VMC for G2 and small checks)
- `multi_e_lattice3d.py` (R7)
- `multi_e_summary.py` → `lab/out/multi_e_results.json` (judging against the pre-registration)
- `multi_e_figure.py` → `lab/out/multi_e_figure.png`
- `multi_e_dev/` (copied development scripts from the scratchpad, listed in the pre-registration)
- `lab/out/multi_e_run_*.sh` (checkpointed drivers)

**Rust:** `src/elettroni/main.rs`, `motore.rs` (VMC/DMC), `guida.rs` (multi-determinant + Jastrow), `orbitali.rs` (splines × Y_lm), `misura.rs`, `tests.rs`, plus one `[[bin]] name = "elettroni"` entry in `Cargo.toml`.

**Docs and assets:** `docs/MOLTI_ELETTRONI.md`, `assets/atomo/molti_elettroni.json`.

`lab/out` is git-ignored: add result files with `git add -f`. Never commit `.npz` or other large files.

---

## 8. The judges' fatal flaws and how they are resolved

**Backbone (cluster-expansion):**
1. **J cannot move nodes, so the Ne energies are set by imposed nodes.** This cannot be fixed inside this plan: it is the fermion sign problem.
   - Node error measured exactly at N = 3 and 4 (D4, D5).
   - Ne node error declared (DEBT-7); N1(b) labelled weak.
   - Optional learned backflow (M3b) is the only learned node-moving piece.
   - Released node is impossible for Ne (signal ∝ e^{−130 t}).
2. **No L–L pairs in training; the r_s scaling is imposed.** Partly fixed.
   - Exact 4-electron Be-like Z = 5, 6 added to training (2s² pairs and the 2s/2p near-degeneracy), after a transfer exam (T1); Be held out.
   - Density scaling becomes DEBT-6, competing against a family without it, chosen by cross-validation.
   - Remaining gap, stated plainly: the 2p–2p pairs of Ne's 2p⁶ are seen only through the 2p² CSF of Be-like ions.
3. **The exact exams test ECG Schrödinger, not the zig-zag.** Fixed by relabelling them as M-class machinery exams, separating out R-class rule exams (R1–R7), and certifying the ECG as the rule's exact level through R4 → D1. The full ε_L is used in R5 and R6.
4. **4-electron ECG risk.**
   - The cost was overestimated: the primitive spin function needs N↑!N↓! = 4 terms, not 24.
   - Staged gate G1 on non-exam systems.
   - Pre-registered downgrade if M5 does not converge.
5. **Force matching inherits oscillating ECG derivatives.** Fixed with exclusion zones, a robust loss, validation against intracules and VMC energy, and a density-matching fallback chosen by cross-validation.
6. **F(s) only covers centred L = 0 states.** Scope stated: P states and molecules go to walkers.
7. **R6 hinges on re-learning; the current k⁴ content is 24–48%.** The existing-level R6 is pre-registered as an expected FAIL; the real exam is on the re-learned level.
8. **R3 noise amplification; quoted numbers include the parity-odd artifact; θ/tanθ typo.** Fixed: G0 then symmetrise, curvature fit, two seeds, tolerance 5e-5, tanθ/θ.
9. **ROHF conflicts with "niente Hartree–Fock".** Cannot be removed without a learned node rule. Mitigated:
   - the mean field is listened from the rule's own operator, from random starts;
   - it never supplies a reported energy;
   - it is formally rejected by the agglomerator's CHECK;
   - DEBT-5, with an explicit decision for Valerio.

**Grafts from imaginary-time-balls:**
1. **Positivity is not learned at the generator level.** Reframed (L2a): the learned item is a bounded-below band; positivity is measured above β_min ≈ 5e-5, and the failure below is disclosed in R2.
2. **"No bottom" was wrong.** Corrected to "sinks to the lattice limit −π/Δt".
3. **Agglomeration in time is decorative.** Dropped.
4. **X3 σ unreachable or tautological.** Now a deterministic three-way check plus a walker leg with a non-exact guide (R4).
5. **X4 doubtful.** Now optional (D0b, σ ≤ 1e-3).
6. **Layzer complex is hindsight.** The generator is mechanical from the rule's measured degeneracy, applied to every atom (it also gives B and C two CSFs); keep/drop is by an internal pilot energy gain; Be is held out of training.
7. **Min-E_FN selection bias and the mixed-density SCF loop.** Pilot selection with independent production; no walker-density SCF loop.
8. **Ne σ optimistic.** Re-budgeted to 1–2 nights, with explicit σ targets and a rule to extend rather than loosen.
9. **Ne tolerance sized from the literature.** Split into a sharp variational part (N1a) and a declared-weak part (N1b).

**Grafts from mean-field-then-correlation:**
1. **Ill-posed cusp exam.** Now opposite spin only, extrapolated estimator, fine bins (D6); same-spin dropped.
2. **E1c vacuous; E1 against HF limits.** E1c is recorded, not scored; E1 is now gate G3, declared as a unit test.
3. **E2b qualitative; E8 report-only.** Replaced by S1, S2 and the C-class thresholds.
4. **Lax ionisation-energy tolerances.** ECG-based ±0.010–0.015 eV; Ne ±0.25 eV with its reason stated.
5. **"Capacity 2 per orbital is not imposed" overstated.** Corrected: 2 internal states come from pauli2; exclusion comes from the imposed antisymmetry.
6. **a0/32 memory.** 4–6 GB fits in 16 GB; H₂ optional.

---

## 9. Risks

1. **Nodes beyond 4 electrons (main risk).** The Ne result rests on mean-field nodes with no exact check. If they are worse than expected, N1(b) and N2 can fail with no remedy inside this plan. Declare the result; do not tune.
2. **4-electron ECG convergence in numpy.** Conditioning and K = 600 may not reach 1.5 mEh in one background session. Fallback: the declared downgrade.
3. **Transfer of the cluster rule** (T1, C2, C3). The density-scaling hypothesis may fail. The VMC numbers are then honest but the many-electron rule is "learned" only in the weak sense; say so.
4. **Re-learned level.** It may not pass its own CHECK at the 3e-5 budget, or R6 may fail. The non-relativistic exams are unaffected; the alpha² claim is withdrawn.
5. **Statistics and time step at Z = 10.** Core time-step bias and the population bias are checked with three τ values and two walker counts. Budget slips of 2× are plausible.
6. **Credibility.** It is quantitatively HF + fixed-node DMC + ECG. The value lies in chain of custody, R-class exams, measured node errors and declared debts. Write that first, not last.
7. **Peeking.** All targets are known values. Mitigations: tolerances and expected outcomes pre-registered, development observations listed, unit tests on non-exam systems, changes after the pre-registration declared.
8. **Shared machine and work tree.** The colla agent is active (uncommitted `lab/colla_*`, `docs/COLLA.md`). Use explicit-path commits, `git pull --rebase --autostash`, never force-push, a new bin only, and run heavy Rust jobs when the cores are free.
9. **Decision dependency.** If Valerio rejects the HF-like carrier (DEBT-5), N ≥ 3 nodes fall back to the Hartree-only carrier, whose cost in accuracy is measured at Li and Be.

---

## 10. Parameter count (to report in the document)

- **Physical inputs:** m_e (through θ), α, Z = 1…10. Nuclear masses of ¹H, ⁴He, ⁷Li, ⁹Be, ²⁰Ne (and B–F) only for μ scaling.
- **Learned from the fine rule:**
  - level-1 stencil (125 taps, symmetrised);
  - D (= 1/2m_kin);
  - m_kin(0);
  - the re-learned Z = 10 level (stencil + F tables), used only for R5/R6 and reports.
- **Learned from exact few-body data:** J coefficients (about 70–80; about 6 more if the three-electron term is selected); optional backflow (about 6).
- **Computed, 0 free:** orbitals, IVO, CI coefficient, Kato cusps, configurations (by the keep rule).
- **Variational / counted:** none in the main path. A VMC-optimised CI coefficient, if reported, counts as 1 parameter, used only as a side check.
- **Method knobs (pre-registered):** K, SVM trials and sweeps, grids, τ triplets, walkers, run lengths, T_pop, thresholds, ridge, exclusion cuts, cross-validation margin.
- **Fitted to any exam quantity:** 0.

---

## 11. Riassunto per Valerio

1. Elio, litio e berillio li risolviamo "esatti" con gaussiane correlate, usando la regola agglomerata (massa imparata, banda imparata) più l'antisimmetria messa a mano e dichiarata: sono il livello di controllo.
2. Da questi casi esatti l'agglomeratore impara la regola a grappoli tra elettroni (coppie, terne) e la prova su ioni che non ha mai visto e sul berillio.
3. Per il neon l'energia la danno le palline in tempo immaginario della regola agglomerata (lecito solo al livello Pauli, e lo verifichiamo); i nodi non sono imparati: vengono dal campo medio, e il loro errore lo misuriamo esatto solo fino a 4 elettroni.
4. Onestà: senza relatività la zig-zag porta agli atomi la massa (che tende a 1), lo spinore a 2 componenti, il permesso del tempo immaginario e il termine (v/c)²; il resto è meccanica quantistica dei fermioni, e lo scriviamo.
5. Costo: circa 30–40 ore di PC in 5–6 sessioni (il neon di notte). Una decisione spetta a te: il campo medio con scambio (in pratica Hartree–Fock ricavato dalla regola) lo usiamo solo per portare i nodi, ma MOLECOLE.md diceva "niente Hartree–Fock".

---

Relevant paths (repo root `C:\Users\valer\Desktop\Programma di simulazione\particle_sim`):
- `docs\PIANO.md`, `docs\AGGLOMERATO_ZIGZAG.md`, `docs\MOLECOLE.md`, `docs\ZIGZAG.md`
- `assets\atomo\agglomerato_zz.json`, `lab\agg_zz_core.py` (`auto_level`, `kcheck` with `energy_mode="state"`), `lab\agg_zz_atom.py`
- `src\molecole\dmc.rs`, `git show a84692d:src/atomo/carbonio.rs`
- Development scripts to copy into `lab\multi_e_dev\`: `C:\Users\valer\AppData\Local\Temp\claude\C--Users-valer-Desktop-Programma-di-simulazione-particle-sim\2978b0e8-1bc3-443b-89d8-2e1ba7531c9b\scratchpad\` (`kcheck.py`…`kcheck4.py`, `judge1_check.py`, `judge1_kpos.py`, `judge2_kbeta.py`, `judge2_disp.py`, `judge2_raw.py`, `he_hf.py`, `he_hf2.py`, `disp.py`, `disp2.py`)
