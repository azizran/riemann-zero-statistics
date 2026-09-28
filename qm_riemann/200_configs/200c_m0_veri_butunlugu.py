"""
200c_m0_veri_butunlugu.py — KALEM 200-C: Kapı M0c (veri bütünlüğü, YALNIZ KONUM)
====================================================================================

GERÇEKTEN ÇALIŞTIRILIR (görev kuralı: M0c yalnız sıfır KONUMLARINI kontrol eder,
hiçbir Z/Z'/M değeri hesaplamaz — körlük ihlali değil). Her dosya için:
  - sha256
  - okunan blok sayısı (1.5e6 sıfıra ulaşana kadar)
  - ilk/son sıfır (mutlak t, taban+ofset)
  - kesin artanlık (strict monotonicity)
  - açılmamış (raw) ortalama aralık VE açılmış (unfolded, *L/2pi) ortalama aralık
    (KALEM: "ortalama aralık 2π/L ile %1 içinde" -> unfolded ortalama ~1.0 olmalı)
  - Riemann-von Mangoldt sayım kontrolü: N_dosya = ilk_blok_N0 + n_max (okunan
    sıfırların GLOBAL kümülatif indeksi) vs N_RvM(T_son) = theta(T_son)/pi + 1;
    |fark| <= 5 (S(T) mertebesi; KALEM'in pilot ile bulduğu +0.32/+0.84/+0.45/+0.24
    ile aynı mertebede beklenir).

Çıktı: M0c_200C.json.
"""
import json
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
import sys
sys.path.insert(0, str(HERE))
import importlib.util


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

GATE_COUNT = 5.0


def main():
    out = {"gate_count_diff": GATE_COUNT, "windows": {}}
    t_start = time.time()
    for w in ortak.WINDOWS:
        t0 = time.time()
        path = ortak.data_path(w)
        sha = ortak.sha256_of(path)
        win = ortak.load_window(w)
        off = win.off
        n_read = win.n_read
        n_blocks_read = len(win.bloklar)

        strict_mono = bool(np.all(np.diff(off) > 0))
        raw_gap_mean = float(np.mean(np.diff(off)))

        L_all = ortak.L_of_offset(win, off)
        Lbar = float(np.mean(L_all))
        unfolded_gap_mean = float(np.mean(np.diff(off) * 0.5 * (L_all[:-1] + L_all[1:]) / ortak.TWO_PI))

        t_first = win.T0 + off[0]
        t_last = win.T0 + off[-1]

        first_block_N0 = int(win.bloklar[0][2])  # (t0,t1,N0,N1)
        global_index_end = first_block_N0 + n_read  # kumulatif GLOBAL sifir sayisi, T_last'a kadar

        N_rvm = float(motor.theta_f64(np.array([t_last]))[0] / np.pi + 1.0)
        count_diff = global_index_end - N_rvm

        rec = {
            "file": ortak.FILES[w],
            "sha256": sha,
            "T0": win.T0,
            "n_blocks_read": n_blocks_read,
            "n_zeros_read": n_read,
            "t_first": float(t_first),
            "t_last": float(t_last),
            "strict_monotonic": strict_mono,
            "raw_gap_mean": raw_gap_mean,
            "unfolded_gap_mean": unfolded_gap_mean,
            "Lbar_pilot": Lbar,
            "L_nominal_KALEM": ortak.L_NOMINAL[w],
            "first_block_N0": first_block_N0,
            "global_zero_index_end": global_index_end,
            "N_RvM_main_term": N_rvm,
            "count_diff_S_T_order": count_diff,
            "gate_count_pass": bool(abs(count_diff) <= GATE_COUNT),
            "gate_unfolded_gap_pass": bool(abs(unfolded_gap_mean - 1.0) <= 0.01),
            "read_time_s": time.time() - t0,
        }
        out["windows"][w] = rec
        print(f"{w}: sha256={sha[:16]}... n_blocks={n_blocks_read} n_zeros={n_read} "
              f"t=[{t_first:.3f},{t_last:.3f}] mono={strict_mono} "
              f"unfolded_gap={unfolded_gap_mean:.6f} Lbar={Lbar:.4f} (KALEM {ortak.L_NOMINAL[w]}) "
              f"count_diff={count_diff:+.3f} -> {'GECTI' if rec['gate_count_pass'] else 'KALDI'}")

    out["overall_pass"] = bool(all(
        out["windows"][w]["gate_count_pass"] and out["windows"][w]["strict_monotonic"]
        and out["windows"][w]["gate_unfolded_gap_pass"] for w in ortak.WINDOWS))
    out["runtime_s"] = time.time() - t_start

    with open(HERE / "M0c_200C.json", "w") as f:
        json.dump(out, f, indent=1)
    print(f"\nM0c GENEL: {'GECTI' if out['overall_pass'] else 'KALDI'}  (toplam {out['runtime_s']:.1f} s)")
    print(f"Kaydedildi: {HERE / 'M0c_200C.json'}")


if __name__ == "__main__":
    main()
