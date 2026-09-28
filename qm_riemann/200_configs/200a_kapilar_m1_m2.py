"""
200a_kapilar_m1_m2.py — KALEM 200-A: Kapı M1 (motor doğruluğu) ve Kapı M2 (sıfır sağlaması)
===============================================================================================

KÖRLÜK NOTU: bu betik gerçek sıfır konumlarında GERÇEK Z ve Z' değerlerini hesaplar
(M1/M2'ye KALEM tarafından açıkça izin verilir: "Engine checks (M1, M2) may compute
real Z and Z' values"). Ama şu ANCAK raporlanır:
  - motor-mpmath FARKLARI (|ΔZ|, |ΔZ'|; max, %99, relatif) ve geçti/kaldı sayıları,
  - M1 için İKİ MOTORUN k2/k3'ünün FARKI (asla mutlak k2/k3 değeri değil),
  - M2 için geçme oranı ve en kötü oran.
Hiçbir yerde gerçek Z / Z' değerlerinin dağılım özeti (ortalama, medyan, histogram,
tek başına k2/k3 vb.) YAZDIRILMAZ / KAYDEDİLMEZ.

Çıktı: M1_200A.json, M2_200A.json (bu dizinde).
"""
import importlib.util
import json
import time
from pathlib import Path
from multiprocessing import Pool

import numpy as np
import mpmath as mp

HERE = Path(__file__).resolve().parent


def _load(modname):
    spec = importlib.util.spec_from_file_location(modname, HERE / f"{modname}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


motor = _load("200a_motor")
ortak = _load("200a_ortak")

TWO_PI = 2 * np.pi
DPS = 20
SEED_T = 2001
SEED_Z = 2002
N_CHECK = 2000

GATE_DZ = 2e-5
GATE_DZP = 2e-4
GATE_K2 = 0.002
GATE_K3 = 0.005


def _mp_worker(t):
    mp.mp.dps = DPS
    tm = mp.mpf(t)
    z = float(mp.siegelz(tm))
    zp = float(mp.siegelz(tm, derivative=1))
    return z, zp


def main():
    tahmin = ortak.load_tahmin()
    zeros = np.load(HERE.parent / "128_odl_zeros6_2e6_zeros.npz")["zeros"]

    rng_t = np.random.default_rng(SEED_T)
    rng_z = np.random.default_rng(SEED_Z)

    per_window = {}
    all_t_points = []   # (window, kind, t_value) flat list for the mpmath pool
    for w in ortak.WINDOWS:
        idx0, idx1 = tahmin[w]["idx"]
        t_a, t_b = tahmin[w]["t"]
        ts = rng_t.uniform(t_a, t_b, N_CHECK)
        zi = rng_z.choice(np.arange(idx0, idx1), size=N_CHECK, replace=False)
        gammas = zeros[zi]
        Ls = np.log(gammas / TWO_PI)
        per_window[w] = dict(ts=ts, gammas=gammas, Ls=Ls)
        all_t_points.extend(ts.tolist())
        all_t_points.extend(gammas.tolist())

    print(f"toplam mpmath noktası: {len(all_t_points)} (dps={DPS}) — havuzla hesaplanıyor...")
    t0 = time.time()
    with Pool(8) as pool:
        mp_res = pool.map(_mp_worker, all_t_points, chunksize=8)
    print(f"  mpmath tamam: {time.time()-t0:.0f} s")
    mp_res = np.array(mp_res)  # (2*3*N_CHECK, 2) columns: Z_mp, Zp_mp

    M1 = {"gates": {"dZ": GATE_DZ, "dZp": GATE_DZP, "dk2": GATE_K2, "dk3": GATE_K3,
                     "dps": DPS, "n_check_per_window": N_CHECK,
                     "seed_t": SEED_T, "seed_zeros": SEED_Z},
          "windows": {}}
    M2 = {"gate_ratio": 1e-4, "n_check_per_window": N_CHECK, "windows": {}}

    off = 0
    overall_pass_M1 = True
    for w in ortak.WINDOWS:
        ts = per_window[w]["ts"]
        gammas = per_window[w]["gammas"]
        Ls = per_window[w]["Ls"]
        n = len(ts)

        Z_mp_t, Zp_mp_t = mp_res[off:off + n, 0], mp_res[off:off + n, 1]
        off += n
        Z_mp_z, Zp_mp_z = mp_res[off:off + n, 0], mp_res[off:off + n, 1]
        off += n

        Z_eng_t, Zp_eng_t = motor.Z_and_Zprime(ts)
        Z_eng_z, Zp_eng_z = motor.Z_and_Zprime(gammas)

        dZ_t = np.abs(Z_eng_t - Z_mp_t)
        dZp_t = np.abs(Zp_eng_t - Zp_mp_t)
        dZ_z = np.abs(Z_eng_z - Z_mp_z)
        dZp_z = np.abs(Zp_eng_z - Zp_mp_z)

        # göreli farklar (payda sıfırdan uzak tutulur: |mpmath| ile normalize)
        relZ_t = dZ_t / np.maximum(np.abs(Z_mp_t), 1e-3)
        relZp_t = dZp_t / np.maximum(np.abs(Zp_mp_t), 1e-3)
        relZ_z = dZ_z / np.maximum(np.abs(Z_mp_z), 1e-3)
        relZp_z = dZp_z / np.maximum(np.abs(Zp_mp_z), 1e-3)

        def stat_block(d, rel):
            return {"max": float(np.max(d)), "p99": float(np.quantile(d, 0.99)),
                    "rel_max": float(np.max(rel)), "rel_p99": float(np.quantile(rel, 0.99))}

        gate_t_Z = bool(np.max(dZ_t) <= GATE_DZ)
        gate_t_Zp = bool(np.max(dZp_t) <= GATE_DZP)
        gate_z_Z = bool(np.max(dZ_z) <= GATE_DZ)
        gate_z_Zp = bool(np.max(dZp_z) <= GATE_DZP)

        # --- k2/k3 FARKI (yalnız fark saklanır) ---
        x_t_eng = np.log(np.abs(Z_eng_t))
        x_t_mp = np.log(np.abs(Z_mp_t))
        x_z_eng = np.log(np.abs(Zp_eng_z) * TWO_PI / Ls)
        x_z_mp = np.log(np.abs(Zp_mp_z) * TWO_PI / Ls)

        k2_t_eng, k3_t_eng = ortak.kstat_from_sums(*ortak.power_sums_of(x_t_eng)[:4])
        k2_t_mp, k3_t_mp = ortak.kstat_from_sums(*ortak.power_sums_of(x_t_mp)[:4])
        k2_z_eng, k3_z_eng = ortak.kstat_from_sums(*ortak.power_sums_of(x_z_eng)[:4])
        k2_z_mp, k3_z_mp = ortak.kstat_from_sums(*ortak.power_sums_of(x_z_mp)[:4])

        dk2_t = float(k2_t_eng - k2_t_mp)
        dk3_t = float(k3_t_eng - k3_t_mp)
        dk2_z = float(k2_z_eng - k2_z_mp)
        dk3_z = float(k3_z_eng - k3_z_mp)
        # değişkenler artık gerekmiyor: bellekten at (M9 ruhu, blindness disiplini)
        del x_t_eng, x_t_mp, x_z_eng, x_z_mp

        gate_k2_t = bool(abs(dk2_t) <= GATE_K2)
        gate_k3_t = bool(abs(dk3_t) <= GATE_K3)
        gate_k2_z = bool(abs(dk2_z) <= GATE_K2)
        gate_k3_z = bool(abs(dk3_z) <= GATE_K3)

        win_pass = all([gate_t_Z, gate_t_Zp, gate_z_Z, gate_z_Zp,
                         gate_k2_t, gate_k3_t, gate_k2_z, gate_k3_z])
        overall_pass_M1 = overall_pass_M1 and win_pass

        M1["windows"][w] = {
            "random_t": {"dZ": stat_block(dZ_t, relZ_t), "dZp": stat_block(dZp_t, relZp_t),
                         "gate_dZ_pass": gate_t_Z, "gate_dZp_pass": gate_t_Zp,
                         "k2_diff_engineB_minus_mpmath": dk2_t,
                         "k3_diff_engineB_minus_mpmath": dk3_t,
                         "gate_k2_pass": gate_k2_t, "gate_k3_pass": gate_k3_t},
            "random_zeros": {"dZ": stat_block(dZ_z, relZ_z), "dZp": stat_block(dZp_z, relZp_z),
                              "gate_dZ_pass": gate_z_Z, "gate_dZp_pass": gate_z_Zp,
                              "k2_diff_engineB_minus_mpmath": dk2_z,
                              "k3_diff_engineB_minus_mpmath": dk3_z,
                              "gate_k2_pass": gate_k2_z, "gate_k3_pass": gate_k3_z},
            "window_pass": win_pass,
        }

        # --- M2: sıfır sağlaması (yalnız motor, mpmath yok) ---
        ratio = np.abs(Z_eng_z) / np.abs(Zp_eng_z)
        pass_mask = ratio <= 1e-4
        M2["windows"][w] = {
            "pass_fraction": float(np.mean(pass_mask)),
            "worst_ratio": float(np.max(ratio)),
            "n": int(n),
        }

    M1["overall_pass"] = overall_pass_M1
    M2["overall_pass"] = bool(all(M2["windows"][w]["pass_fraction"] == 1.0 for w in ortak.WINDOWS))
    # M2'nin KALEM eşiği: 2000 sıfırın HEPSİ oranı sağlamalı (aksi hâlde motor tablo ile
    # uyumsuz sayılır); daha yumuşak bir "en kötü oran <= eşik" testi de ekleyelim:
    M2["overall_pass_worst_ratio"] = bool(all(M2["windows"][w]["worst_ratio"] <= 1e-4 for w in ortak.WINDOWS))

    with open(HERE / "M1_200A.json", "w") as f:
        json.dump(M1, f, indent=1)
    with open(HERE / "M2_200A.json", "w") as f:
        json.dump(M2, f, indent=1)

    print("\n=== M1 (motor doğruluğu) ===")
    for w in ortak.WINDOWS:
        rt = M1["windows"][w]["random_t"]
        rz = M1["windows"][w]["random_zeros"]
        print(f"{w}: t-örnek |dZ| max={rt['dZ']['max']:.2e} p99={rt['dZ']['p99']:.2e} "
              f"| |dZp| max={rt['dZp']['max']:.2e} | dk2={rt['k2_diff_engineB_minus_mpmath']:+.5f} "
              f"dk3={rt['k3_diff_engineB_minus_mpmath']:+.5f} -> pencere geçti={M1['windows'][w]['window_pass']}")
        print(f"    sıfır-örnek |dZ| max={rz['dZ']['max']:.2e} | |dZp| max={rz['dZp']['max']:.2e} "
              f"| dk2={rz['k2_diff_engineB_minus_mpmath']:+.5f} dk3={rz['k3_diff_engineB_minus_mpmath']:+.5f}")
    print(f"M1 GENEL: {'GEÇTİ' if overall_pass_M1 else 'KALDI'}")

    print("\n=== M2 (sıfır sağlaması) ===")
    for w in ortak.WINDOWS:
        m2w = M2["windows"][w]
        print(f"{w}: geçme oranı={m2w['pass_fraction']*100:.2f}%  en kötü oran={m2w['worst_ratio']:.2e}")
    print(f"M2 GENEL (hepsi geçti): {'GEÇTİ' if M2['overall_pass'] else 'KALDI'}  "
          f"(en kötü oran <= 1e-4: {'GEÇTİ' if M2['overall_pass_worst_ratio'] else 'KALDI'})")


if __name__ == "__main__":
    main()
