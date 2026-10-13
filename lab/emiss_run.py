"""emiss_run.py - the coupled light + atom runs of the study "Emissione vera" (docs/EMISSIONE.md), in pieces of
<= ~9 minutes with checkpoints (exit code 3 = "run me again").

Modes
  sc : semiclassical (as src/atomo/accoppiato.rs): one electron wave, one real light field; the electron feels
       e E.r (E = light at the atom + an optional prescribed drive), the light is sourced by -e d<r>/dt.
  qe : one-excitation quantum (no RWA on the walks, RWA on the coupling): vacuum sector psi (starts as a pure
       excited state), photon sector = channel states phi_c (lattice 1s, evolved by the same walk) times a
       complex photon amplitude field per channel (two light walks: real and imaginary part).
       psi += -i a e (Phi_c(atom) . r) phi_c ;  source of Phi_c: J_c = -e d<phi_c| r |psi>/dt.
  light: no atom: a prescribed dipole current I0 cos(w t) z (smooth start), to check the light rule's radiation.

Usage: python lab/emiss_run.py CONFIG_NAME   (configs in emiss_configs.py)
"""
import sys, os, time, json
import numpy as np
from emiss_core import *


class Runner:
    def __init__(self, cfg):
        self.cfg = cfg
        self.name = cfg["name"]
        self.ck = os.path.join(CK, f"run_{self.name}.npz")
        c = cfg
        self.e = E_CHARGE * c.get("K", 1.0)
        self.mode = c["mode"]
        self.b = c["b"]
        self.nL = c["nL"]
        self.sub = c.get("sub", 1)
        if self.mode != "light":
            self.atom = c["atom_factory"]()
            self.a = self.atom.a
            assert abs(self.sub * self.a - self.b) < 1e-9 * self.b, "b must be sub * a"
        else:
            self.a = self.b / self.sub
        nf = 1 if self.mode in ("sc", "light") else 2 * len(c["channels"](self.atom))
        self.L = [Light3D(self.nL, self.b, sponge=c.get("sponge", 12), sigma=c.get("sigma", 1.0)) for _ in range(nf)]
        self.t_ticks = 0
        self.ser = {k: [] for k in ("t", "norm", "pexc", "eph", "abs", "probe", "dip", "Dc", "pgnd")}
        self.dser = []             # (t, <r>) every atom tick (sc)
        self.probe_idx = (self.nL // 2 - c.get("probe", 40), self.nL // 2, self.nL // 2)
        if c.get("transverse"):
            for L in self.L:
                L.make_transverse()
        self.Jhist = [[] for _ in range(max(1, nf // 2))]

    # ---------------------------------------------------------------- state
    def init_state(self):
        c = self.cfg
        if self.mode == "light":
            return
        self.psi = c["psi0_factory"](self.atom)
        self.excited = [p.copy() for p in c.get("excited", lambda at: [])(self.atom)]
        self.ground = c.get("ground", lambda at: None)(self.atom)
        if self.mode == "qe":
            self.phis = [f.copy() for f in c["channels"](self.atom)]
            self.record(np.array([self.atom.dip_matrix(p, self.psi) for p in self.phis]))
        else:
            self.record(self.atom.position(self.psi))

    def save(self):
        d = dict(t_ticks=self.t_ticks, F=np.array([L.F for L in self.L]), absorbed=np.array([L.absorbed for L in self.L]),
                 ltick=np.array([L.tick_count for L in self.L]), dser=np.array(self.dser))
        for k, v in self.ser.items():
            d["ser_" + k] = np.array(v)
        if self.mode != "light":
            d["psi"] = self.psi
            d["excited"] = np.array(self.excited) if self.excited else np.zeros(0)
            if self.ground is not None:
                d["ground"] = self.ground
            if self.mode == "qe":
                d["phis"] = np.array(self.phis)
        tmp = self.ck + ".tmp.npz"
        np.savez(tmp, **d)
        os.replace(tmp, self.ck)

    def load(self):
        d = np.load(self.ck, allow_pickle=False)
        self.t_ticks = int(d["t_ticks"])
        for i, L in enumerate(self.L):
            L.F = d["F"][i].copy()
            L.absorbed = float(d["absorbed"][i])
            L.tick_count = int(d["ltick"][i])
        self.dser = [tuple(x) for x in d["dser"]] if d["dser"].size else []
        for k in self.ser:
            self.ser[k] = list(d["ser_" + k])
        if self.mode != "light":
            self.psi = d["psi"].copy()
            self.excited = list(d["excited"]) if d["excited"].size else []
            self.ground = d["ground"].copy() if "ground" in d else None
            if self.mode == "qe":
                self.phis = list(d["phis"])

    # ---------------------------------------------------------------- one light tick
    def drive(self, t):
        dr = self.cfg.get("drive")
        if dr is None:
            return np.zeros(3)
        x = (t - dr["t0"]) / dr["tau"]
        if abs(x) > 4:
            return np.zeros(3)
        env = dr["eE0"] * np.exp(-x * x)
        ph = dr["omega"] * (t - dr["t0"])
        if dr["pol"] == "z":
            return np.array([0.0, 0.0, env * np.cos(ph)])
        if dr["pol"] == "sigma-":          # E = Re[(x - i y) e^{-i w t}] : drives Delta m = -1
            return np.array([env * np.cos(ph), -env * np.sin(ph), 0.0])
        raise ValueError(dr["pol"])

    def light_tick(self):
        a, e, b = self.a, self.e, self.b
        if self.mode == "light":
            c = self.cfg
            t = self.t_ticks * a
            ramp = 1 - np.exp(-(t / c["ramp"])**2) if c.get("ramp") else 1.0
            I = np.array([0, 0, c["I0"] * ramp * np.cos(c["omega"] * t)])
            self.t_ticks += self.sub
            self.L[0].tick()
            self.L[0].deposit(I)
            self.L[0].sponge()
            self.record(I)
            return
        if self.mode == "sc":
            r0 = self.atom.position(self.psi)
            centered = self.cfg.get("centered", False)
            L0 = self.L[0]
            if centered:
                L0.tick()
                h = self.Jhist[0]
                if len(h) >= 2 and np.dot(h[-2], h[-2]) > 0:
                    # real current: extrapolate linearly (J ~ cos): J_n ~ 2 cos(wb) J_{n-1} - J_{n-2}, cos from history
                    jp = 2 * h[-1] - h[-2] if len(h) < 3 else self.sc_predict(h)
                else:
                    jp = h[-1] if h else np.zeros(3)
                Elight = L0.field_at_atom() + 0.5 * self.selff(L0, jp)
            else:
                Elight = L0.field_at_atom()
            for s in range(self.sub):
                t = self.t_ticks * a
                eE = e * Elight + self.drive(t)
                self.atom.tick(self.psi, self.t_ticks, eE)
                self.t_ticks += 1
                if self.cfg.get("dipole_every_tick"):
                    self.dser.append((self.t_ticks * a,) + tuple(self.atom.position(self.psi)))
            r1 = self.atom.position(self.psi)
            I = -e * (r1 - r0) / b
            if centered:
                dE = 0.5 * self.selff(L0, I - jp)
                # correction of the self field: a phase e dE.r applied once for the whole light tick
                self.atom.tick_phase(self.psi, e * dE * b)
            else:
                L0.tick()
            self.Jhist[0].append(I)
            self.dep(L0, I)
            L0.sponge()
            self.record(r1)
            return
        # qe
        nc = len(self.phis)
        D0 = [self.atom.dip_matrix(p, self.psi) for p in self.phis]
        centered = self.cfg.get("centered", False)
        if centered:
            for L in self.L:              # free propagation first: F' = U F_n
                L.tick()
            fields, Jp = [], []
            for c in range(nc):
                h = self.Jhist[c]
                if len(h) >= 2 and np.vdot(h[-2], h[-2]).real > 0:
                    r = np.vdot(h[-2], h[-1]) / np.vdot(h[-2], h[-2])
                    r = r / max(1.0, abs(r))
                    jp = h[-1] * r
                elif len(h) == 1:
                    jp = h[-1]
                else:
                    jp = np.zeros(3, complex)
                Jp.append(jp)
                LA, LB = self.L[2 * c], self.L[2 * c + 1]
                fields.append(LA.field_at_atom() + 0.5 * self.selff(LA, jp.real)
                              + 1j * (LB.field_at_atom() + 0.5 * self.selff(LB, jp.imag)))
        else:
            fields = [self.L[2 * c].field_at_atom() + 1j * self.L[2 * c + 1].field_at_atom() for c in range(nc)]
        if hasattr(self.atom, "qe_block"):
            self.atom.qe_block(self.psi, self.phis, fields, self.sub, a, e)
            self.t_ticks += self.sub
        else:
            for s in range(self.sub):
                self.atom.tick(self.psi, self.t_ticks)
                for p in self.phis:
                    self.atom.tick(p, self.t_ticks)
                for c in range(nc):
                    self.atom.add_rphi(self.psi, -1j * a * e * fields[c], self.phis[c])
                self.t_ticks += 1
        D1 = [self.atom.dip_matrix(p, self.psi) for p in self.phis]
        for c in range(nc):
            J = -e * (D1[c] - D0[c]) / b
            if centered:
                LA, LB = self.L[2 * c], self.L[2 * c + 1]
                dE = 0.5 * (self.selff(LA, (J - Jp[c]).real) + 1j * self.selff(LB, (J - Jp[c]).imag))
                self.atom.add_rphi(self.psi, -1j * b * e * dE, self.phis[c])
            self.Jhist[c].append(J)
            for q, part in ((0, J.real), (1, J.imag)):
                L = self.L[2 * c + q]
                if not centered:
                    L.tick()
                self.dep(L, part)
                L.sponge()
        self.record(np.array(D1))

    def sc_predict(self, h):
        """linear prediction J_n = 2c J_{n-1} - J_{n-2} with c fitted on the last three (cos(w b))"""
        x, y, z = h[-3], h[-2], h[-1]
        den = 2 * np.dot(y, y)
        c = np.dot(z + x, y) / den if den > 0 else 1.0
        c = min(max(c, -1.0), 1.0)
        return 2 * c * z - y

    def selff(self, L, I):
        if self.cfg.get("transverse"):
            return L.self_field(I)
        return -np.sum(L.w**2) * np.asarray(I) / L.b**2

    def dep(self, L, I):
        if self.cfg.get("transverse"):
            L.deposit_T(I)
        else:
            L.deposit(I)

    def record(self, extra):
        a = self.a
        t = self.t_ticks * a
        s = self.ser
        s["t"].append(t)
        eph = sum(L.energy() for L in self.L)
        ab = sum(L.absorbed for L in self.L)
        s["eph"].append(eph)
        s["abs"].append(ab)
        pr = []
        for L in self.L:
            pr += list(L.F[:, self.probe_idx[0], self.probe_idx[1], self.probe_idx[2]].real)
        s["probe"].append(pr)
        if self.mode == "light":
            s["dip"].append(list(extra))
            return
        s["norm"].append(norm2(self.psi))
        ev = self.t_ticks % 2 == 0
        s["pexc"].append([abs(ovl(x, self.psi))**2 if ev else np.nan for x in self.excited])
        s["pgnd"].append(abs(ovl(self.ground, self.psi))**2 if (self.ground is not None and ev) else np.nan)
        if self.mode == "sc":
            s["dip"].append(list(extra))
        else:
            s["Dc"].append([complex(v) for v in np.ravel(extra)])

    # ---------------------------------------------------------------- run in pieces
    def run(self, budget=540):
        t0 = time.time()
        if os.path.exists(self.ck):
            self.load()
            print(f"[{self.name}] resumed at t = {self.t_ticks * self.a:.1f}", flush=True)
        else:
            self.init_state()
        T = self.cfg["T"]
        last = time.time()
        while self.t_ticks * self.a < T - 1e-9:
            self.light_tick()
            if time.time() - last > 60:
                last = time.time()
                s = self.ser
                print(f"[{self.name}] t = {s['t'][-1]:.0f}/{T}  norm {s['norm'][-1] if s['norm'] else 0:.8f}  "
                      f"eph+abs {s['eph'][-1] + s['abs'][-1]:.4e}  ({time.time() - t0:.0f} s)", flush=True)
            if time.time() - t0 > budget:
                self.save()
                print(f"[{self.name}] checkpoint at t = {self.t_ticks * self.a:.1f}; run again", flush=True)
                return False
        self.save()
        print(f"[{self.name}] done ({time.time() - t0:.0f} s)", flush=True)
        return True


if __name__ == "__main__":
    from emiss_configs import CONFIGS
    name = sys.argv[1]
    r = Runner(CONFIGS[name])
    ok = r.run(float(sys.argv[2]) if len(sys.argv) > 2 else 540)
    sys.exit(0 if ok else 3)
