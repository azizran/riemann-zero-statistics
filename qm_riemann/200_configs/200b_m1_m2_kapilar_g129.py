"""
200b_m1_m2_kapilar.py — KALEM 200-B: Kapı M1b (tepe motoru) ve Kapı M2b (ızgara
çözünürlüğü)
===============================================================================

KÖRLÜK NOTU (görev ABSOLUTE RULES + KALEM Kapılar M1b/M2b): bu betik pencere başına
300 rastgele yakın-çift olayında (δ̃<0.2, tohum 2101) GERÇEK M_n değerlerini hesaplar
(motor+ızgara+parabol VE mpmath golden-section; M2b için 33 vs 129 nokta, motor
YALNIZ) — KALEM'de açıkça izin verilir. Ama diske/ekrana YALNIZ:
  - motor-mpmath |ΔM|/M FARKLARI (max, %99) ve geçti/kaldı,
  - M1b için İKİ YÖNTEMİN log(M/δ̃²) k2/k3'ünün FARKI (asla mutlak k2/k3 değeri),
  - M2b için |ΔM|/M (33 vs 129) farkı ve geçti/kaldı.
YAZILIR. Hiçbir yerde gerçek M_n / x_n dağılım özeti (ortalama, medyan, histogram,
tek başına k2/k3, vb.) YAZDIRILMAZ / KAYDEDİLMEZ.

Olay seçimi: pencere başına TÜM δ̃<0.2 olaylarından (yalnız KONUM — γ_n, γ_{n+1},
δ̃_n; 200b_ortak.events_from_zeros) 300 tanesi rastgele seçilir (tohum 2101, KALEM).

M1b yöntemi:
  - motor: 200a_motor.Z_vec, 33 iç nokta ızgara + 3-nokta parabolik rafinasyon
    (200b_ortak.peak_grid_parabola, |Z| üzerinde).
  - mpmath: dps>=20, altın-oran (golden-section) maksimizasyonu |siegelz(t)| üzerinde,
    tolerans <=1e-12 (t'de RELATİF: durma koşulu genişlik <= 1e-12*|orta nokta|).
  - Kapılar: |ΔM|/M <= 1e-4 (her olay); |k2 farkı| <= 0.002, |k3 farkı| <= 0.005
    (pencere başına, 300 olayın x=log(M/δ̃²)'sinden, İKİ yöntem arasında).

M2b yöntemi: aynı 900 olay (300/pencere), motor 33 vs 129 nokta ızgara: |ΔM|/M <= 1e-6.

Çıktı: M1b_200B_g129.json, M2b_200B_g129.json (bu dizinde).
"""
import argparse
import importlib.util
import json
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import mpmath as mp

HERE = Path(__file__).resolve().parent


def _load(modname):
    spec = importlib.util.spec_from_file_location(modname, HERE / f"{modname}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


ortak_b = _load("200b_ortak")
motor = _load("200a_motor")
ortak_a = _load("200a_ortak")  # yalnız kstat_from_sums / power_sums_of için (körlük ihlali değil: yardımcı fonksiyon)

SEED_EVENTS = 2101
N_CHECK_PER_WINDOW = 300
EPS_CHECK = 0.2
DPS = 20
TOL_REL_T = 1e-12
MAX_ITER = 300

GATE_DM_M1B = 1e-4
GATE_K2_M1B = 0.002
GATE_K3_M1B = 0.005
GATE_DM_M2B = 1e-6
N_GRID_LO = 129
N_GRID_HI = 513


# ---------------------------------------------------------------------------
# mpmath altın-oran (golden-section) |Z(t)| maksimizasyonu — tek olay (Pool worker)
# ---------------------------------------------------------------------------
def _mp_golden_worker(args):
    ga, gb, dps, tol_rel, max_iter = args
    mp.mp.dps = dps
    invphi = (mp.sqrt(5) - 1) / 2  # 1/golden ratio ~ 0.618

    def f(t):
        return abs(mp.siegelz(t))

    a = mp.mpf(ga)
    b = mp.mpf(gb)
    c = b - invphi * (b - a)
    d = a + invphi * (b - a)
    fc = f(c)
    fd = f(d)
    tol_abs = mp.mpf(tol_rel) * abs((a + b) / 2)
    it = 0
    while (b - a) > tol_abs and it < max_iter:
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - invphi * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + invphi * (b - a)
            fd = f(d)
        it += 1
    tbest = (a + b) / 2
    Mbest = max(float(fc), float(fd), float(f(tbest)))
    return Mbest, float(tbest), it


def mpmath_peak_batch(ga, gb, n_workers=8, dps=DPS, tol_rel=TOL_REL_T, max_iter=MAX_ITER):
    jobs = [(float(g0), float(g1), dps, tol_rel, max_iter) for g0, g1 in zip(ga, gb)]
    with Pool(n_workers) as pool:
        res = pool.map(_mp_golden_worker, jobs, chunksize=4)
    M = np.array([r[0] for r in res], dtype=np.float64)
    iters = np.array([r[2] for r in res], dtype=np.int64)
    return M, iters


# ---------------------------------------------------------------------------
# Olay seçimi: pencere başına 300 rastgele δ̃<0.2 olayı (tohum 2101; yalnız KONUM)
# ---------------------------------------------------------------------------
def select_check_events(n_per_window=N_CHECK_PER_WINDOW, seed=SEED_EVENTS, eps=EPS_CHECK):
    real_windows = ortak_b.load_real_windows()
    zeros = np.load(HERE.parent / "128_odl_zeros6_2e6_zeros.npz")["zeros"]
    rng = np.random.default_rng(seed)
    out = {}
    for w, info in real_windows.items():
        i0, i1 = info["idx"]
        ev = ortak_b.events_from_zeros(zeros[i0:i1], i0, eps_max=eps)
        n_avail = len(ev["n"])
        take = min(n_per_window, n_avail)
        sel = rng.choice(n_avail, size=take, replace=False)
        out[w] = {k: v[sel] for k, v in ev.items()}
        out[w]["n_available_deltatilde_lt_eps"] = n_avail
    return out


def stat_block(rel):
    return {"max": float(np.max(rel)), "p99": float(np.quantile(rel, 0.99)),
            "mean": float(np.mean(rel))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--n-per-window", type=int, default=N_CHECK_PER_WINDOW)
    args = ap.parse_args()

    print(f"Olay seçimi: pencere başına {args.n_per_window} rastgele δ̃<{EPS_CHECK} olayı "
          f"(tohum {SEED_EVENTS}, yalnız KONUM)...")
    events = select_check_events(n_per_window=args.n_per_window)
    for w, ev in events.items():
        print(f"  {w}: mevcut δ̃<{EPS_CHECK} olay={ev['n_available_deltatilde_lt_eps']}, "
              f"seçilen={len(ev['n'])}")

    M1 = {"gates": {"dM_rel": GATE_DM_M1B, "dk2": GATE_K2_M1B, "dk3": GATE_K3_M1B,
                     "dps": DPS, "tol_rel_t": TOL_REL_T, "n_per_window": args.n_per_window,
                     "seed": SEED_EVENTS, "eps_check": EPS_CHECK},
          "windows": {}}
    M2 = {"gate_dM_rel": GATE_DM_M2B, "n_grid_lo": N_GRID_LO, "n_grid_hi": N_GRID_HI,
          "n_per_window": args.n_per_window, "windows": {}}

    overall_pass_M1 = True
    overall_pass_M2 = True
    t0 = time.time()

    # ---- mpmath golden-section: TEK Pool, TÜM pencerelerin olayları birlikte
    # (Pool başlatma/mpmath import ek yükünü pencere başına DEĞİL, bir kez öder) ----
    win_order = list(events.keys())
    ga_all = np.concatenate([events[w]["ga"] for w in win_order])
    gb_all = np.concatenate([events[w]["gb"] for w in win_order])
    print(f"  mpmath golden-section: {len(ga_all)} olay, {args.workers} işçi ile başlıyor...")
    M_mp_all, mp_iters_all = mpmath_peak_batch(ga_all, gb_all, n_workers=args.workers)
    print(f"  mpmath tamam: {time.time()-t0:.1f} s")
    off = 0
    M_mp_by_win, mp_iters_by_win = {}, {}
    for w in win_order:
        n_ev = len(events[w]["ga"])
        M_mp_by_win[w] = M_mp_all[off:off + n_ev]
        mp_iters_by_win[w] = mp_iters_all[off:off + n_ev]
        off += n_ev

    for w, ev in events.items():
        ga, gb, dt_tilde = ev["ga"], ev["gb"], ev["delta_tilde"]
        n_ev = len(ga)

        # ---- motor: 33 iç nokta + parabol ----
        M_eng33 = ortak_b.peak_grid_parabola(motor.Z_vec, ga, gb, n_grid=N_GRID_LO)
        # ---- motor: 129 iç nokta + parabol (M2b) ----
        M_eng129 = ortak_b.peak_grid_parabola(motor.Z_vec, ga, gb, n_grid=N_GRID_HI)
        # ---- mpmath: golden-section (M1b), önceden hesaplandı ----
        M_mp, mp_iters = M_mp_by_win[w], mp_iters_by_win[w]

        # --- M1b: motor(33) vs mpmath ---
        rel_1b = np.abs(M_eng33 - M_mp) / np.maximum(M_mp, 1e-300)
        gate_dM_1b = bool(np.max(rel_1b) <= GATE_DM_M1B)

        x_eng = np.log(M_eng33 / dt_tilde ** 2)
        x_mp = np.log(M_mp / dt_tilde ** 2)
        k2_eng, k3_eng = ortak_a.kstat_from_sums(*ortak_a.power_sums_of(x_eng)[:4])
        k2_mp, k3_mp = ortak_a.kstat_from_sums(*ortak_a.power_sums_of(x_mp)[:4])
        dk2 = float(k2_eng - k2_mp)
        dk3 = float(k3_eng - k3_mp)
        del x_eng, x_mp  # körlük disiplini: mutlak k2/k3 asla saklanmaz, yalnız fark

        gate_k2 = bool(abs(dk2) <= GATE_K2_M1B)
        gate_k3 = bool(abs(dk3) <= GATE_K3_M1B)
        win_pass_1b = gate_dM_1b and gate_k2 and gate_k3
        overall_pass_M1 = overall_pass_M1 and win_pass_1b

        M1["windows"][w] = {
            "n_checked": int(n_ev),
            "dM_rel": stat_block(rel_1b),
            "gate_dM_pass": gate_dM_1b,
            "k2_diff_engine_minus_mpmath": dk2,
            "k3_diff_engine_minus_mpmath": dk3,
            "gate_k2_pass": gate_k2,
            "gate_k3_pass": gate_k3,
            "mp_golden_iters": {"max": int(mp_iters.max()), "mean": float(mp_iters.mean())},
            "window_pass": win_pass_1b,
        }

        # --- M2b: 33 vs 129 (motor yalnız) ---
        rel_2b = np.abs(M_eng129 - M_eng33) / np.maximum(M_eng129, 1e-300)
        gate_2b = bool(np.max(rel_2b) <= GATE_DM_M2B)
        overall_pass_M2 = overall_pass_M2 and gate_2b
        # Ek TANI (KALEM'de İSTENMİYOR ama KALDI durumunu bağlamlandırmak için): 33
        # vs 129'un x=log(M/δ̃²) k2/k3'üne etkisi — yalnız FARK (körlük ihlali değil,
        # M1b ile aynı disiplin).
        x_e33 = np.log(M_eng33 / dt_tilde ** 2)
        x_e129 = np.log(M_eng129 / dt_tilde ** 2)
        k2_33, k3_33 = ortak_a.kstat_from_sums(*ortak_a.power_sums_of(x_e33)[:4])
        k2_129, k3_129 = ortak_a.kstat_from_sums(*ortak_a.power_sums_of(x_e129)[:4])
        dk2_grid = float(k2_33 - k2_129)
        dk3_grid = float(k3_33 - k3_129)
        del x_e33, x_e129
        M2["windows"][w] = {"n_checked": int(n_ev), "dM_rel": stat_block(rel_2b),
                             "gate_pass": gate_2b,
                             "diagnostic_k2_diff_33_minus_129": dk2_grid,
                             "diagnostic_k3_diff_33_minus_129": dk3_grid}

        print(f"  {w}: M1b |ΔM|/M max={rel_1b.max():.3e} p99={np.quantile(rel_1b,0.99):.3e} "
              f"dk2={dk2:+.6f} dk3={dk3:+.6f} -> {'GEÇTİ' if win_pass_1b else 'KALDI'} "
              f"| M2b |ΔM|/M max={rel_2b.max():.3e} -> {'GEÇTİ' if gate_2b else 'KALDI'}")

    M1["overall_pass"] = overall_pass_M1
    M2["overall_pass"] = overall_pass_M2
    dt_total = time.time() - t0
    M1["runtime_s"] = dt_total

    with open(HERE / "M1b_200B_g129.json", "w") as f:
        json.dump(M1, f, indent=1)
    with open(HERE / "M2b_200B_g129.json", "w") as f:
        json.dump(M2, f, indent=1)

    print(f"\nM1b GENEL: {'GEÇTİ' if overall_pass_M1 else 'KALDI'}")
    print(f"M2b GENEL: {'GEÇTİ' if overall_pass_M2 else 'KALDI'}")
    print(f"Süre: {dt_total:.1f} s")
    print("Kaydedildi: M1b_200B_g129.json, M2b_200B_g129.json")


if __name__ == "__main__":
    main()
