"""Ponte cloud <-> portatile: the laptop as the GPU laboratory of the cloud session.

The cloud session cannot see the laptop's screen or GPU. It writes a request file in
`ponte/richieste/<id>.json` and pushes it to GitHub. This script, running on the laptop, pulls
the repo every minute or two, executes each new request from a FIXED list of allowed actions,
and pushes the result to `ponte/risultati/<id>/` (esito.json, log.txt, images).

Safety: the runner never executes free-form commands from the repo. Only these actions exist:
    ping   machine info (git revision, GPU, CPU, RAM)
    build  cargo build --bin atomo
    test   cargo test --bin atomo [-- filter] [--ignored]
    foto   run the atomo binary in --foto mode with validated arguments, return the PNG
Every argument is checked against a whitelist / regex below. Building and running the repo's own
code is still running the repo's code: only point this at your own repository.

It works in a SEPARATE clone (default: ../ponte_lab next to the repo), so it never touches the
working copy you edit by hand.

Run:   python ponte/ponte.py            (loop forever, poll every 90 s)
       python ponte/ponte.py --una      (one pass and exit)
Stop:  Ctrl+C
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO_URL = "https://github.com/valeriozanzarelli-debug/particle-sim.git"
BRANCH = "master"
ATTRIBUTION = "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"

ID_RE = re.compile(r"^[A-Za-z0-9_\-]{1,64}$")
FILTER_RE = re.compile(r"^[A-Za-z0-9_:]{0,80}$")

# foto arguments: name -> validator
def _choice(*vals):
    return lambda v: str(v) in vals

def _regex(rx):
    r = re.compile(rx)
    return lambda v: bool(r.fullmatch(str(v)))

def _num(lo, hi, integer=False):
    def ok(v):
        try:
            x = int(v) if integer else float(v)
        except (TypeError, ValueError):
            return False
        return lo <= x <= hi
    return ok

FOTO_ARGS = {
    "scena": _choice("fotone", "atomo", "orbitali", "elettrone", "tavola"),
    "elemento": _regex(r"[A-Za-z0-9]{1,8}"),
    "riga": _choice("2p", "3p"),
    "zoom": _choice("intero", "nucleo", "protone", "elettrone"),
    "stati": _regex(r"[0-9a-z+]{1,24}"),
    "qual": _num(0, 3, integer=True),
    "qualita": _choice("bassa", "media", "medio-alta", "alta", "ultra"),
    "polar": _choice("h", "v", "c+", "c-"),
    "dir": _num(0, 4, integer=True),
    "w0": _num(0.1, 100.0),
    "lambda-m": _num(1e-15, 1e5),
    "strati": _regex(r"[A-Za-z]{1,16}"),
    "frames": _num(1, 20000, integer=True),
    "spara-al": _num(0, 20000, integer=True),
    "misure": _num(0, 100000, integer=True),
    "ritmo": _num(0.0, 1000.0),
    "yaw": _num(-10.0, 10.0),
    "pitch": _num(-10.0, 10.0),
    "distanza": _num(0.01, 1e6),
    "foto-dopo": _num(0.0, 5000.0),
    "scia": _num(0.0, 2.0),
    "lod": _num(0.0, 5.0),
}
FOTO_FLAGS = {"surfisti", "collassa", "pannello", "spara-pronto", "senza-diffusa", "buca", "fetta", "zoom-atomo", "senza-taglio"}

TIMEOUT = {"ping": 60, "build": 3600, "test": 3600, "foto": 3600}


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def run(cmd: list[str], cwd: Path, timeout: int) -> tuple[int, str]:
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout,
                           encoding="utf-8", errors="replace")
        return p.returncode, (p.stdout or "") + (p.stderr or "")
    except subprocess.TimeoutExpired as e:
        out = (e.stdout or "") if isinstance(e.stdout, str) else ""
        return 124, out + f"\n[ponte] TIMEOUT after {timeout} s"
    except FileNotFoundError as e:
        return 127, f"[ponte] command not found: {e}"


def git(args: list[str], cwd: Path, timeout: int = 300) -> tuple[int, str]:
    return run(["git", *args], cwd, timeout)


def ensure_clone(lab: Path) -> None:
    if (lab / ".git").exists():
        return
    lab.parent.mkdir(parents=True, exist_ok=True)
    code, out = run(["git", "clone", "--branch", BRANCH, REPO_URL, str(lab)], lab.parent, 1800)
    if code != 0:
        sys.exit(f"[ponte] clone failed:\n{out}")


def sync(lab: Path) -> bool:
    code, out = git(["pull", "--rebase", "--autostash", "origin", BRANCH], lab)
    if code != 0:
        print(f"[ponte] pull failed, retry next cycle:\n{out[-800:]}")
        git(["rebase", "--abort"], lab)
        return False
    return True


def gpu_info() -> str:
    if platform.system() == "Windows":
        code, out = run(["powershell", "-NoProfile", "-Command",
                         "Get-CimInstance Win32_VideoController | ForEach-Object { $_.Name + ' | ' + $_.AdapterRAM }"],
                        Path.cwd(), 60)
        return out.strip()
    code, out = run(["sh", "-c", "lspci | grep -i -E 'vga|3d' || true"], Path.cwd(), 60)
    return out.strip()


def cargo() -> str:
    return shutil.which("cargo") or str(Path.home() / ".cargo" / "bin" / "cargo")


def exe_path(lab: Path) -> Path:
    name = "atomo.exe" if platform.system() == "Windows" else "atomo"
    return lab / "target" / "debug" / name


def act_ping(req: dict, lab: Path, out_dir: Path) -> tuple[bool, str, dict]:
    _, rev = git(["rev-parse", "HEAD"], lab)
    info = {
        "git": rev.strip(),
        "gpu": gpu_info(),
        "cpu": platform.processor(),
        "cores": os.cpu_count(),
        "os": platform.platform(),
        "python": sys.version.split()[0],
        "cargo": cargo(),
    }
    return True, json.dumps(info, indent=2), info


def act_build(req: dict, lab: Path, out_dir: Path) -> tuple[bool, str, dict]:
    code, out = run([cargo(), "build", "--bin", "atomo", "-j", "2"], lab, TIMEOUT["build"])
    return code == 0, out, {"exit": code}


def act_test(req: dict, lab: Path, out_dir: Path) -> tuple[bool, str, dict]:
    filt = str(req.get("args", {}).get("filtro", ""))
    if not FILTER_RE.fullmatch(filt):
        return False, f"[ponte] filtro non valido: {filt!r}", {}
    cmd = [cargo(), "test", "--bin", "atomo", "-j", "2", "--", "--test-threads=1"]
    if filt:
        cmd.append(filt)
    if req.get("args", {}).get("ignored"):
        cmd.append("--ignored")
    code, out = run(cmd, lab, TIMEOUT["test"])
    return code == 0, out, {"exit": code}


def act_foto(req: dict, lab: Path, out_dir: Path) -> tuple[bool, str, dict]:
    ok_b, out_b, _ = act_build(req, lab, out_dir)
    if not ok_b:
        return False, "[ponte] build failed before foto\n" + out_b[-4000:], {}
    args = req.get("args", {})
    cmd = [str(exe_path(lab))]
    for k, v in args.items():
        if k in FOTO_FLAGS:
            if v:
                cmd.append(f"--{k}")
            continue
        if k not in FOTO_ARGS:
            return False, f"[ponte] argomento non ammesso: {k!r}", {}
        if not FOTO_ARGS[k](v):
            return False, f"[ponte] valore non valido per {k}: {v!r}", {}
        cmd += [f"--{k}", str(v)]
    png = out_dir / "foto.png"
    cmd += ["--foto", str(png)]
    t0 = time.time()
    code, out = run(cmd, lab, TIMEOUT["foto"])
    extra = {"exit": code, "secondi": round(time.time() - t0, 1), "comando": " ".join(cmd[1:])}
    ok = code == 0 and png.exists()
    if not png.exists():
        out += "\n[ponte] nessuna immagine prodotta"
    return ok, out, extra


ACTIONS = {"ping": act_ping, "build": act_build, "test": act_test, "foto": act_foto}


def handle(req_file: Path, lab: Path) -> str | None:
    rid = req_file.stem
    if not ID_RE.fullmatch(rid):
        return None
    out_dir = lab / "ponte" / "risultati" / rid
    if (out_dir / "esito.json").exists():
        return None
    try:
        req = json.loads(req_file.read_text(encoding="utf-8"))
    except Exception as e:  # malformed request: record it, do not crash
        req = {"azione": "?", "errore": str(e)}
    out_dir.mkdir(parents=True, exist_ok=True)
    azione = str(req.get("azione", "?"))
    start = now()
    t0 = time.time()
    if azione in ACTIONS:
        try:
            ok, log, extra = ACTIONS[azione](req, lab, out_dir)
        except Exception as e:
            ok, log, extra = False, f"[ponte] eccezione: {e!r}", {}
    else:
        ok, log, extra = False, f"[ponte] azione sconosciuta: {azione!r} (ammesse: {sorted(ACTIONS)})", {}
    _, rev = git(["rev-parse", "HEAD"], lab)
    esito = {
        "id": rid, "azione": azione, "ok": ok, "inizio": start, "fine": now(),
        "durata_s": round(time.time() - t0, 1), "git": rev.strip(), "macchina": platform.node(),
        **extra,
    }
    (out_dir / "esito.json").write_text(json.dumps(esito, indent=2), encoding="utf-8")
    lines = log.splitlines()
    (out_dir / "log.txt").write_text("\n".join(lines[-600:]) + "\n", encoding="utf-8")
    return f"{rid} ({azione}, {'ok' if ok else 'FALLITO'})"


def publish(lab: Path, done: list[str]) -> None:
    git(["add", "-f", "ponte/risultati"], lab)
    msg = "ponte: risultati " + ", ".join(done) + "\n\n" + ATTRIBUTION
    code, out = git(["commit", "-m", msg], lab)
    if code != 0:
        print(f"[ponte] nothing to commit or commit failed:\n{out[-400:]}")
        return
    for _ in range(3):
        code, out = git(["push", "origin", BRANCH], lab)
        if code == 0:
            print(f"[ponte] pushed: {', '.join(done)}")
            return
        sync(lab)
    print(f"[ponte] push failed:\n{out[-800:]}")


def one_pass(lab: Path) -> None:
    if not sync(lab):
        return
    req_dir = lab / "ponte" / "richieste"
    if not req_dir.exists():
        return
    done = []
    for f in sorted(req_dir.glob("*.json")):
        r = handle(f, lab)
        if r:
            print(f"[ponte] {now()} done {r}")
            done.append(r)
    if done:
        publish(lab, done)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    here = Path(__file__).resolve().parents[1]
    ap.add_argument("--lab", default=str(here.parent / "ponte_lab"), help="separate clone used by the runner")
    ap.add_argument("--ogni", type=int, default=90, help="seconds between polls")
    ap.add_argument("--una", action="store_true", help="one pass and exit")
    a = ap.parse_args()
    lab = Path(a.lab)
    ensure_clone(lab)
    print(f"[ponte] lab clone: {lab}  poll every {a.ogni} s  actions: {sorted(ACTIONS)}")
    while True:
        try:
            one_pass(lab)
        except KeyboardInterrupt:
            raise
        except Exception as e:
            print(f"[ponte] errore nel ciclo: {e!r}")
        if a.una:
            break
        time.sleep(a.ogni)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[ponte] fermato.")
