"""
200c_m1_m2_kapilar.py — KALEM 200-C: Kapı M1c (motor doğruluğu) + Kapı M2c (sıfır sağlaması)
================================================================================================

KÖRLÜK NOTU: KALEM'in açıkça izin verdiği gibi ("Engine checks may compute real
Z/Z' values"), bu betik gerçek çapalı motorumuzla VE mpmath.siegelz ile GERÇEK
Z(t), Z'(t) değerlerini hesaplar — ama yalnız İKİ MOTORUN FARKLARI (|ΔZ|, |ΔZ'|,
k2/k3 FARKI) ve geçti/kaldı sayıları raporlanır/kaydedilir. Hiçbir mutlak k2/k3
ya da Z/Z' dağılım özeti YAZDIRILMAZ.

M1c: pencere başına 200 rastgele t (tohum 3101; adayların YARISI çapaya EN UZAK
noktalardan seçilir — geniş bir aday havuzundan |dt| en büyük 100 nokta + 100
sıradan rastgele nokta) ve 200 rastgele sıfır (tohum 3102), motor vs mpmath
(siegelz, dps>=25). Gates: |ΔZ|<=2e-5, |ΔZ'|<=2e-4; iki metodun k2/k3 FARKI
(x=log|Z| t-noktalarında, x=log(|Z'|*2pi/L) sıfırlarda) <= 0.002/0.005.

M2c: aynı 200 sıfırda |Z(gamma_n)| <= 1e-4*|Z'(gamma_n)| (yalnız motor).

Çıktı: M1c_200C.json, M2c_200C.json.
"""
import importlib.util
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import mpmath as mp

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
ortak_a = _load("200a_ortak")

TWO_PI = 2 * np.pi
DPS = 25
SEED_T = 3101
SEED_Z = 3102
N_CHECK = 200
N_POOL_FARTHEST = 2000  # aday havuzu (en-uzak yarı secimi icin)

GATE_DZ = 2e-5
GATE_DZP = 2e-4
GATE_K2 = 0.002
GATE_K3 = 0.005
GATE_M2_RATIO = 1e-4


def _mp_worker(args):
    T0, t_off, dps = args
    mp.mp.dps = dps
    tm = mp.mpf(T0) + mp.mpf(float(t_off))
    z = float(mp.siegelz(tm))
    zp = float(mp.siegelz(tm, derivative=1))
    return z, zp


def pick_t_points(win, tab, seed, n_check, n_pool):
    """n_check nokta: yarisi cok GENIS bir havuzdan capaya EN UZAK |dt|'ye sahip
    olanlar, yarisi siradan rastgele (KALEM M1c kurali)."""
    rng = np.random.default_rng(seed)
    pool = rng.uniform(win.t_a_off, win.t_b_off, n_pool)
    aidx = np.clip(np.round((pool - tab.offsets[0]) / tab.spacing).astype(np.int64),
                    0, len(tab.offsets) - 1)
    dt = np.abs(pool - tab.offsets[aidx])
    order = np.argsort(-dt)  # en uzaktan en yakina
    n_half = n_check // 2
    farthest = pool[order[:n_half]]
    regular = rng.uniform(win.t_a_off, win.t_b_off, n_check - n_half)
    return np.concatenate([farthest, regular])


def main():
    t_start = time.time()
    M1 = {"gates": {"dZ": GATE_DZ, "dZp": GATE_DZP, "dk2": GATE_K2, "dk3": GATE_K3,
                     "dps": DPS, "n_check_per_window": N_CHECK,
                     "seed_t": SEED_T, "seed_zeros": SEED_Z},
          "windows": {}}
    M2 = {"gate_ratio": GATE_M2_RATIO, "n_check_per_window": N_CHECK, "windows": {}}

    per_window_data = {}
    mp_jobs = []
    job_index = {}  # (w, kind) -> (start, n)

    for w in ortak.WINDOWS:
        win = ortak.load_window(w)
        tab = ortak.get_anchor_table(win, n_workers=8)

        ts = pick_t_points(win, tab, SEED_T, N_CHECK, N_POOL_FARTHEST)
        rng_z = np.random.default_rng(SEED_Z)
        zi = rng_z.choice(win.n_read - 1, size=N_CHECK, replace=False)  # icerideki indeks (0-based)
        zt = win.off[zi]
        Ls = ortak.L_of_offset(win, zt)

        Z_eng_t, Zp_eng_t = motor.Z_and_Zprime(ts, tab)
        Z_eng_z, Zp_eng_z = motor.Z_and_Zprime(zt, tab)

        per_window_data[w] = dict(win=win, tab=tab, ts=ts, zt=zt, Ls=Ls,
                                   Z_eng_t=Z_eng_t, Zp_eng_t=Zp_eng_t,
                                   Z_eng_z=Z_eng_z, Zp_eng_z=Zp_eng_z)

        start = len(mp_jobs)
        for tv in ts:
            mp_jobs.append((win.T0, tv, DPS))
        job_index[(w, "t")] = (start, N_CHECK)
        start = len(mp_jobs)
        for tv in zt:
            mp_jobs.append((win.T0, tv, DPS))
        job_index[(w, "z")] = (start, N_CHECK)

    print(f"toplam mpmath noktasi: {len(mp_jobs)} (dps={DPS}) -- havuzla hesaplaniyor...")
    t0 = time.time()
    with Pool(8) as pool:
        mp_res = pool.map(_mp_worker, mp_jobs, chunksize=4)
    print(f"  mpmath tamam: {time.time()-t0:.0f} s")
    mp_res = np.array(mp_res)

    overall_pass_M1 = True
    for w in ortak.WINDOWS:
        d = per_window_data[w]
        s_t, n_t = job_index[(w, "t")]
        s_z, n_z = job_index[(w, "z")]
        Z_mp_t, Zp_mp_t = mp_res[s_t:s_t + n_t, 0], mp_res[s_t:s_t + n_t, 1]
        Z_mp_z, Zp_mp_z = mp_res[s_z:s_z + n_z, 0], mp_res[s_z:s_z + n_z, 1]

        dZ_t = np.abs(d["Z_eng_t"] - Z_mp_t)
        dZp_t = np.abs(d["Zp_eng_t"] - Zp_mp_t)
        dZ_z = np.abs(d["Z_eng_z"] - Z_mp_z)
        dZp_z = np.abs(d["Zp_eng_z"] - Zp_mp_z)

        relZ_t = dZ_t / np.maximum(np.abs(Z_mp_t), 1e-3)
        relZp_t = dZp_t / np.maximum(np.abs(Zp_mp_t), 1e-3)
        relZ_z = dZ_z / np.maximum(np.abs(Z_mp_z), 1e-3)
        relZp_z = dZp_z / np.maximum(np.abs(Zp_mp_z), 1e-3)

        def stat_block(dd, rel):
            return {"max": float(np.max(dd)), "p99": float(np.quantile(dd, 0.99)),
                    "rel_max": float(np.max(rel)), "rel_p99": float(np.quantile(rel, 0.99))}

        gate_t_Z = bool(np.max(dZ_t) <= GATE_DZ)
        gate_t_Zp = bool(np.max(dZp_t) <= GATE_DZP)
        gate_z_Z = bool(np.max(dZ_z) <= GATE_DZ)
        gate_z_Zp = bool(np.max(dZp_z) <= GATE_DZP)

        x_t_eng = np.log(np.abs(d["Z_eng_t"]))
        x_t_mp = np.log(np.abs(Z_mp_t))
        x_z_eng = np.log(np.abs(d["Zp_eng_z"]) * TWO_PI / d["Ls"])
        x_z_mp = np.log(np.abs(Zp_mp_z) * TWO_PI / d["Ls"])

        k2_t_eng, k3_t_eng = ortak_a.kstat_from_sums(*ortak_a.power_sums_of(x_t_eng)[:4])
        k2_t_mp, k3_t_mp = ortak_a.kstat_from_sums(*ortak_a.power_sums_of(x_t_mp)[:4])
        k2_z_eng, k3_z_eng = ortak_a.kstat_from_sums(*ortak_a.power_sums_of(x_z_eng)[:4])
        k2_z_mp, k3_z_mp = ortak_a.kstat_from_sums(*ortak_a.power_sums_of(x_z_mp)[:4])

        dk2_t = float(k2_t_eng - k2_t_mp)
        dk3_t = float(k3_t_eng - k3_t_mp)
        dk2_z = float(k2_z_eng - k2_z_mp)
        dk3_z = float(k3_z_eng - k3_z_mp)
        del x_t_eng, x_t_mp, x_z_eng, x_z_mp

        gate_k2_t = bool(abs(dk2_t) <= GATE_K2)
        gate_k3_t = bool(abs(dk3_t) <= GATE_K3)
        gate_k2_z = bool(abs(dk2_z) <= GATE_K2)
        gate_k3_z = bool(abs(dk3_z) <= GATE_K3)

        win_pass = all([gate_t_Z, gate_t_Zp, gate_z_Z, gate_z_Zp,
                         gate_k2_t, gate_k3_t, gate_k2_z, gate_k3_z])
        overall_pass_M1 = overall_pass_M1 and win_pass

        M1["windows"][w] = {
            "anchor_spacing": d["tab"].spacing, "n_anchors": int(len(d["tab"].offsets)),
            "random_t": {"dZ": stat_block(dZ_t, relZ_t), "dZp": stat_block(dZp_t, relZp_t),
                         "gate_dZ_pass": gate_t_Z, "gate_dZp_pass": gate_t_Zp,
                         "k2_diff_engineC_minus_mpmath": dk2_t,
                         "k3_diff_engineC_minus_mpmath": dk3_t,
                         "gate_k2_pass": gate_k2_t, "gate_k3_pass": gate_k3_t},
            "random_zeros": {"dZ": stat_block(dZ_z, relZ_z), "dZp": stat_block(dZp_z, relZp_z),
                              "gate_dZ_pass": gate_z_Z, "gate_dZp_pass": gate_z_Zp,
                              "k2_diff_engineC_minus_mpmath": dk2_z,
                              "k3_diff_engineC_minus_mpmath": dk3_z,
                              "gate_k2_pass": gate_k2_z, "gate_k3_pass": gate_k3_z},
            "window_pass": win_pass,
        }

        ratio = np.abs(d["Z_eng_z"]) / np.abs(d["Zp_eng_z"])
        pass_mask = ratio <= GATE_M2_RATIO
        M2["windows"][w] = {"pass_fraction": float(np.mean(pass_mask)),
                             "worst_ratio": float(np.max(ratio)), "n": int(n_z)}

    M1["overall_pass"] = overall_pass_M1
    M2["overall_pass"] = bool(all(M2["windows"][w]["pass_fraction"] == 1.0 for w in ortak.WINDOWS))
    M2["overall_pass_worst_ratio"] = bool(all(M2["windows"][w]["worst_ratio"] <= GATE_M2_RATIO
                                               for w in ortak.WINDOWS))
    M1["runtime_s"] = time.time() - t_start
    M2["runtime_s"] = M1["runtime_s"]

    with open(HERE / "M1c_200C.json", "w") as f:
        json.dump(M1, f, indent=1)
    with open(HERE / "M2c_200C.json", "w") as f:
        json.dump(M2, f, indent=1)

    print("\n=== M1c (motor dogrulugu) ===")
    for w in ortak.WINDOWS:
        rt = M1["windows"][w]["random_t"]
        rz = M1["windows"][w]["random_zeros"]
        print(f"{w}: spacing={M1['windows'][w]['anchor_spacing']} n_anchors={M1['windows'][w]['n_anchors']} "
              f"| t |dZ| max={rt['dZ']['max']:.2e} |dZp| max={rt['dZp']['max']:.2e} "
              f"dk2={rt['k2_diff_engineC_minus_mpmath']:+.6f} dk3={rt['k3_diff_engineC_minus_mpmath']:+.6f} "
              f"-> pencere gecti={M1['windows'][w]['window_pass']}")
        print(f"    zero |dZ| max={rz['dZ']['max']:.2e} |dZp| max={rz['dZp']['max']:.2e} "
              f"dk2={rz['k2_diff_engineC_minus_mpmath']:+.6f} dk3={rz['k3_diff_engineC_minus_mpmath']:+.6f}")
    print(f"M1c GENEL: {'GECTI' if overall_pass_M1 else 'KALDI'}")

    print("\n=== M2c (sifir saglamasi) ===")
    for w in ortak.WINDOWS:
        m2w = M2["windows"][w]
        print(f"{w}: gecme orani={m2w['pass_fraction']*100:.2f}%  en kotu oran={m2w['worst_ratio']:.2e}")
    print(f"M2c GENEL: {'GECTI' if M2['overall_pass'] else 'KALDI'}")
    print(f"\nToplam sure: {M1['runtime_s']:.1f} s")


if __name__ == "__main__":
    main()
