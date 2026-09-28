"""
200c_m3_tepe.py — KALEM 200-C: Kapı M3c (tepe/hump bulucu çözünürlüğü)
=================================================================================

Yakın çift olaylarında (δ̃<0.2) M_n = max|Z| tepesinin 129-nokta ızgara+parabol
ile 513-nokta ızgara+parabol arasındaki farkı — YALNIZ FARK raporlanır (mutlak
M dağılımı yok). Pencere başına 200 rastgele olay (tohum 3103). Gate: |ΔM|/M <= 1e-6.

Çıktı: M3c_200C.json.
"""
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def _load(modname):
    if modname in sys.modules:
        return sys.modules[modname]
    spec = importlib.util.spec_from_file_location(modname, HERE / f"{modname}.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[modname] = m
    spec.loader.exec_module(m)
    return m


ortak = _load("200c_ortak")
motor = _load("200c_motor")

SEED_EV = 3103
N_CHECK = 200
EPS_FOR_M3 = 0.2
GRID_A = 129
GRID_B = 513
GATE_REL = 1e-6


def main():
    t_start = time.time()
    out = {"gate_rel_dM": GATE_REL, "grid_a": GRID_A, "grid_b": GRID_B,
           "seed": SEED_EV, "n_check_per_window": N_CHECK, "windows": {}}

    overall_pass = True
    for w in ortak.WINDOWS:
        win = ortak.load_window(w)
        tab = ortak.get_anchor_table(win, n_workers=8)
        Z_func = ortak.Z_func_for(win, tab)

        ev = ortak.events_from_window(win, eps_max=EPS_FOR_M3)
        n_ev = len(ev["n"])
        rng = np.random.default_rng(SEED_EV)
        keep = rng.choice(n_ev, size=min(N_CHECK, n_ev), replace=False)
        ga, gb = ev["ga"][keep], ev["gb"][keep]

        t0 = time.time()
        Ma = ortak.ortak_b.peak_grid_parabola(Z_func, ga, gb, n_grid=GRID_A)
        Mb = ortak.ortak_b.peak_grid_parabola(Z_func, ga, gb, n_grid=GRID_B)
        dt_calc = time.time() - t0

        rel = np.abs(Ma - Mb) / np.maximum(Mb, 1e-300)
        gate_pass = bool(np.max(rel) <= GATE_REL)
        overall_pass = overall_pass and gate_pass

        out["windows"][w] = {
            "n_events_eps02": int(n_ev), "n_checked": int(len(ga)),
            "rel_dM_max": float(np.max(rel)), "rel_dM_p99": float(np.quantile(rel, 0.99)),
            "rel_dM_mean": float(np.mean(rel)), "gate_pass": gate_pass,
            "calc_time_s": dt_calc,
        }
        print(f"{w}: n_events(eps<0.2)={n_ev} checked={len(ga)} "
              f"rel_dM max={np.max(rel):.3e} p99={np.quantile(rel,0.99):.3e} -> "
              f"{'GECTI' if gate_pass else 'KALDI'} ({dt_calc:.1f}s)")

    out["overall_pass"] = overall_pass
    out["runtime_s"] = time.time() - t_start
    with open(HERE / "M3c_200C.json", "w") as f:
        json.dump(out, f, indent=1)
    print(f"\nM3c GENEL: {'GECTI' if overall_pass else 'KALDI'}  (toplam {out['runtime_s']:.1f} s)")
    print(f"Kaydedildi: {HERE / 'M3c_200C.json'}")


if __name__ == "__main__":
    main()
