"""Independent check of the Compton binary (docs/COMPTON.md).

Runs `target/debug/compton.exe --eventi N` (surfer rule, events born inside the scattering disk),
reads the CSV and checks, with Klein-Nishina written again here (scipy.integrate.quad):
  1. Compton's formula for every event (lambda' - lambda vs lambda_C (1 - cos theta)), < 1e-9
  2. angular distribution vs Klein-Nishina (chi2, 40 bins in cos theta) at 17.48, 662, 1332.5 keV
     and the same at several seeds (mean chi2 / dof), polarized azimuth at 90 deg
  3. electron kinetic-energy spectrum and Compton edge at 662 keV (477.65 keV)
  4. Compton 1923: Mo K-alpha at 90 deg
  5. momentum conservation per event (< 1e-9), electron angle tan(phi) = cot(theta/2)/(1+k)
Figure: lab/out/compton_validation.png.   Run:  timeout 590 python lab/compton_check.py
"""
import os
import subprocess
import sys
import time

import numpy as np
from scipy import integrate, stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXE = os.path.join(ROOT, "target", "debug", "compton.exe")
OUT = os.path.join(ROOT, "lab", "out")
os.makedirs(OUT, exist_ok=True)

ME = 510.99895000          # keV
HC = 1239.841984332        # keV pm
LC = 2.42631023867         # pm
N = 20000

results = []


def check(name, ok, detail):
    results.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")


def kn_dmu(k, mu):
    p = 1.0 / (1.0 + k * (1.0 - mu))
    return np.pi * p * p * (p + 1.0 / p - (1.0 - mu * mu))


def kn_dmu_dphi_pol(k, mu, phi):
    p = 1.0 / (1.0 + k * (1.0 - mu))
    return 0.5 * p * p * (p + 1.0 / p - 2.0 * (1.0 - mu * mu) * np.cos(phi) ** 2)


def run(e, seed, pol="u", n=N):
    f = os.path.join(OUT, f"compton_ev_{e:g}_{pol}_{seed}.csv")
    subprocess.run([EXE, "--eventi", str(n), "--kev", repr(e), "--seed", str(seed), "--pol", pol, "--out", f],
                   check=True, cwd=ROOT, capture_output=True)
    d = np.loadtxt(f, delimiter=",", comments="#", skiprows=2)
    cols = "theta,phi,e_out,lam_out,pgx,pgy,pgz,pex,pey,pez,t_e,theta_e,beta_e,bx,by".split(",")
    return {c: d[:, i] for i, c in enumerate(cols)}


def chi2_mu(ev, e, nb=40):
    k = e / ME
    edges = np.linspace(-1, 1, nb + 1)
    tot = integrate.quad(lambda m: kn_dmu(k, m), -1, 1, epsabs=0, epsrel=1e-12)[0]
    exp = np.array([integrate.quad(lambda m: kn_dmu(k, m), a, b, epsabs=0, epsrel=1e-12)[0]
                    for a, b in zip(edges[:-1], edges[1:])]) / tot * len(ev["theta"])
    h, _ = np.histogram(np.cos(ev["theta"]), edges)
    c = float(((h - exp) ** 2 / exp).sum())
    return c, nb - 1, h, exp, tot


t0 = time.time()
if not os.path.exists(EXE):
    sys.exit("build first: cargo build --bin compton")

# ---------- 1, 5: kinematics per event ----------
data = {}
for e, seed in [(17.48, 1), (662.0, 2), (1332.5, 3)]:
    ev = run(e, seed)
    data[e] = ev
    lam = HC / e
    k = e / ME
    # 1. Compton's formula, with lambda' from the photon momentum
    pg = np.sqrt(ev["pgx"] ** 2 + ev["pgy"] ** 2 + ev["pgz"] ** 2)
    lam2 = HC / pg
    err1 = np.max(np.abs((lam2 - lam) - LC * (1 - np.cos(ev["theta"]))) / lam2)
    check(f"formula di Compton per evento, {e} keV", err1 < 1e-9, f"max err rel {err1:.2e}")
    # 5. momentum and energy per event
    r = np.sqrt((ev["pgx"] + ev["pex"]) ** 2 + (ev["pgy"] + ev["pey"]) ** 2 + (ev["pgz"] + ev["pez"] - e) ** 2) / e
    pe2 = ev["pex"] ** 2 + ev["pey"] ** 2 + ev["pez"] ** 2
    ee = np.sqrt(pe2 + ME ** 2)
    re = np.abs(e + ME - pg - ee) / (e + ME)
    check(f"impulso ed energia conservati evento per evento, {e} keV", r.max() < 1e-9 and re.max() < 1e-9,
          f"max |dp|/p {r.max():.1e}, max |dE|/E {re.max():.1e}")
    th = ev["theta"]
    m = (th > 1e-6) & (th < np.pi - 1e-6)
    phe = np.arctan2(np.cos(th[m] / 2), np.sin(th[m] / 2) * (1 + k))
    err5 = np.max(np.abs(phe - ev["theta_e"][m]))
    check(f"angolo dell'elettrone tan phi = cot(theta/2)/(1+k), {e} keV", err5 < 1e-9, f"max {err5:.1e} rad")

# ---------- 2: angular distribution ----------
for e in data:
    c, dof, *_ = chi2_mu(data[e], e)
    p = stats.chi2.sf(c, dof)
    check(f"distribuzione angolare vs Klein-Nishina, {e} keV, 20 000 eventi", p > 1e-3, f"chi2 = {c:.1f}/{dof}, p = {p:.3f}")
# several seeds: the chi2 must scatter around dof
for e in [17.48, 662.0]:
    cs = [chi2_mu(run(e, 1000 + s, n=10000), e)[0] for s in range(8)]
    m, se = np.mean(cs), np.std(cs, ddof=1) / np.sqrt(len(cs))
    ks = stats.kstest(stats.chi2.cdf(cs, 39), "uniform").pvalue
    check(f"chi2 su 8 seed (10 000 eventi), {e} keV", ks > 1e-3, f"media {m:.1f} +- {se:.1f} (atteso 39), KS p = {ks:.2f}")
# polarized azimuth at 90 deg
ev = run(17.48, 7, pol="x", n=40000)
k = 17.48 / ME
band = (ev["theta"] > np.radians(80)) & (ev["theta"] < np.radians(100))
edges = np.linspace(0, 2 * np.pi, 13)
w = np.array([integrate.dblquad(lambda ph, mu: kn_dmu_dphi_pol(k, mu, ph), np.cos(np.radians(100)), np.cos(np.radians(80)),
                                a, b)[0] for a, b in zip(edges[:-1], edges[1:])])
exp = w / w.sum() * band.sum()
h, _ = np.histogram(np.mod(ev["phi"][band], 2 * np.pi), edges)
c = float(((h - exp) ** 2 / exp).sum())
p = stats.chi2.sf(c, 11)
check("azimut con polarizzazione lineare, 80-100 gradi", p > 1e-3, f"chi2 = {c:.1f}/11, p = {p:.3f}")

# ---------- 3: electron spectrum, Compton edge ----------
ev = data[662.0]
k = 662.0 / ME
edge = 662.0 * 2 * k / (1 + 2 * k)
check("bordo Compton a 662 keV (formula)", abs(edge - 477.65) < 0.005, f"{edge:.3f} keV")
check("nessun elettrone sopra il bordo, il piu' energetico vicino al bordo",
      ev["t_e"].max() <= edge * (1 + 1e-12) and ev["t_e"].max() > 0.995 * edge, f"T max = {ev['t_e'].max():.3f} keV")
# spectrum vs KN: dsigma/dT = dsigma/dmu * dmu/dT
tb = np.linspace(0, 1.05 * edge, 51)
def mu_of_t(t):
    t = np.minimum(t, edge)
    e2 = 662.0 - t
    return np.clip(1 - (662.0 / e2 - 1) / k, -1, 1)
tot = integrate.quad(lambda m: kn_dmu(k, m), -1, 1, epsrel=1e-12)[0]
expt = np.array([integrate.quad(lambda m: kn_dmu(k, m), mu_of_t(b), mu_of_t(a), epsrel=1e-12)[0]
                 for a, b in zip(tb[:-1], tb[1:])]) / tot * len(ev["t_e"])
ht, _ = np.histogram(ev["t_e"], tb)
sel = expt >= 5
c = float(((ht[sel] - expt[sel]) ** 2 / expt[sel]).sum())
p = stats.chi2.sf(c, sel.sum() - 1)
check("spettro di energia dell'elettrone a 662 keV vs Klein-Nishina", p > 1e-3, f"chi2 = {c:.1f}/{sel.sum() - 1}, p = {p:.3f}")

# ---------- 4: Compton 1923 ----------
# Phys. Rev. 21, 483 (1923): Mo K-alpha (0.708 A) scattered by graphite at 90 deg -> 0.730 A (shift 0.022 A)
meas = 2.2
ours = LC * (1 - np.cos(np.pi / 2))
check("Compton 1923, Mo Ka a 90 gradi (precisione storica ~0.2 pm)", abs(meas - ours) < 0.3,
      f"misurato {meas:.1f} pm, nostro {ours:.4f} pm ({100 * (meas - ours) / ours:+.1f} %)")

# ---------- figure ----------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for e, col in [(17.48, "C0"), (662.0, "C1"), (1332.5, "C3")]:
        c, dof, h, exp, _ = chi2_mu(data[e], e)
        x = np.linspace(-1, 1, 41)
        xm = 0.5 * (x[1:] + x[:-1])
        ax[0].step(xm, h, where="mid", color=col, label=f"{e:g} keV (surfisti)")
        ax[0].plot(xm, exp, "k--", lw=0.8)
    ax[0].set_xlabel("cos θ del fotone"); ax[0].set_ylabel("eventi / classe"); ax[0].legend(fontsize=8)
    ax[0].set_title("angoli: surfisti vs Klein–Nishina (tratteggio)")
    tm = 0.5 * (tb[1:] + tb[:-1])
    ax[1].step(tm, ht, where="mid", label="surfisti, 662 keV")
    ax[1].plot(tm, expt, "k--", lw=0.8, label="Klein–Nishina")
    ax[1].axvline(477.65, color="m", lw=1, label="bordo Compton 477,65 keV")
    ax[1].set_xlabel("T elettrone (keV)"); ax[1].legend(fontsize=8); ax[1].set_title("spettro dell'elettrone")
    th = np.linspace(0, 180, 200)
    ax[2].plot(th, LC * (1 - np.cos(np.radians(th))), label="λ' − λ = λc(1 − cos θ)")
    evx = data[17.48]
    ax[2].plot(np.degrees(evx["theta"][::40]), (evx["lam_out"] - HC / 17.48)[::40], ".", ms=2, label="eventi (17,48 keV)")
    ax[2].errorbar([90], [2.2], yerr=[0.2], fmt="s", color="r", label="Compton 1923, Mo Kα")
    ax[2].set_xlabel("θ (gradi)"); ax[2].set_ylabel("Δλ (pm)"); ax[2].legend(fontsize=8)
    ax[2].set_title("spostamento di lunghezza d'onda")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "compton_validation.png"), dpi=110)
    print("figure -> lab/out/compton_validation.png")
except ImportError:
    print("matplotlib missing: no figure")

npass = sum(ok for _, ok, _ in results)
print(f"\n{npass} PASS / {len(results)} test   ({time.time() - t0:.1f} s)")
