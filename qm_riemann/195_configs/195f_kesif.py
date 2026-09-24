# -*- coding: utf-8 -*-
"""
195f — ÖN-KAYITSIZ KEŞİF (hükümler görüldükten SONRA; HÜKÜM DIŞI)
================================================================
Soru: B'de yasak uydularda görülen küçük POZİTİF Re K̃ (~+0.004), uydu-dışı
"boş" konumlarda da var mı (evrensel taban)? Boş konum: Δω ∈ [−1.0, 2.5],
a,b ≤ 10 sade rasyonellerin log'larına (Δω=0 dahil) uzaklığı ≥ 0.045.
Aynı makine (195o.cekirdek_blok), aynı bloklar; adalar (üst bölge) + zeta düşük
ızgarası (17 blok) + zeta son (193 blokları, ALT=3).
Çıktı: scratchpad/195/K_kesif_195.json
"""
import importlib.util
import json
import time
from math import gcd
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("o195", HERE / "195o_ortak.py")
o = importlib.util.module_from_spec(spec)
spec.loader.exec_module(o)


def log(*a):
    print(*a, flush=True)


_W = {}


def _is(arg):
    ad, j, idx_f, idx_a, Lb, adlar = arg
    D = _W[ad]
    r = o.cekirdek_blok(D["g"], D["mid"], idx_f, idx_a, Lb, D["Q"], D["A"], D["CHI"], adlar)
    return ad, j, r["num"], r["den"], r["N_tam"]


if __name__ == "__main__":
    T0 = time.time()
    rats = sorted({np.log(a / b) for a in range(1, 11) for b in range(1, 11) if gcd(a, b) == 1})
    aday = np.round(np.arange(-1.0, 2.5001, 0.005), 3)
    bos = [x for x in aday if np.min(np.abs(np.array(rats) - x)) >= 0.045]
    # birbirinden ≥ 0.1 uzak olanları seç (seyrek, dağınık)
    sec = []
    for x in bos:
        if not sec or x - sec[-1] >= 0.1:
            sec.append(float(x))
    adlar = []
    for x in sec:
        et = f"bos{x:+.3f}"
        o.UYDULAR[et] = (float(np.exp(x)), 1.0)
        adlar.append(et)
    log(f"boş konumlar ({len(sec)}): {sec}")
    tab, _ = o.karakter_tablolari()
    isler = []
    Q, _, A = o.asal_kuvvetler(np.exp(12.1 + 2.5 + 0.04))
    for ad in o.ADALAR:
        K = o.ada_kin(ad)
        _W[ad] = dict(g=K["g"], mid=K["mid"], Q=Q, A=A,
                      CHI=o.chi_dizi(tab[ad], K["k"], Q.astype(np.int64)))
        for (j, lo, hi, idx) in o.dL_bloklari(K["Lm"], o.L_UST):
            isler.append((ad, j, idx, idx, float(K["Lm"][idx].mean()), adlar))
    Kd = o.zeta_kin("dusuk")
    _W["zeta_dusuk"] = dict(g=Kd["g"], mid=Kd["mid"], Q=Q, A=A, CHI=np.ones(len(Q)))
    for (j, lo, hi, idx) in o.dL_bloklari(Kd["Lm"], o.L_UST):
        isler.append(("zeta_dusuk", j, idx, idx, float(Kd["Lm"][idx].mean()), adlar))
    Ks = o.zeta_kin("son")
    Qs, _, As = o.asal_kuvvetler(np.exp(12.1 + 2.5 + 0.04))
    _W["zeta_son"] = dict(g=Ks["g"], mid=Ks["mid"], Q=Qs, A=As, CHI=np.ones(len(Qs)))
    alt = np.arange(0, len(Ks["mid"]), 3)
    for b, idx in enumerate(o.esit_sayim_bloklari(len(Ks["mid"]))):
        ia = alt[(alt >= idx[0]) & (alt <= idx[-1])]
        isler.append(("zeta_son", b, idx, ia, float(np.log(Ks["mid"][idx] / o.TWO_PI).mean()), adlar))
    import multiprocessing as mp
    R = {}
    with mp.get_context("fork").Pool(8) as pool:
        for ad, j, num, den, N in pool.imap_unordered(_is, isler):
            R.setdefault(ad, []).append((j, num, den, N))
    out = {"konumlar": sec, "not": "ÖN-KAYITSIZ KEŞİF — hüküm dışı", "sonuc": {}}
    for ad in o.ADALAR + ["zeta_dusuk", "zeta_son"]:
        rr = sorted(R[ad], key=lambda x: x[0])
        num = np.array([x[1] for x in rr])
        den = np.array([x[2] for x in rr])
        Nw = np.array([x[3] for x in rr])
        K, reps = o.havuz(num, den, Nw)
        s = o.ozet(K, reps)
        out["sonuc"][ad] = {"re": s["re"].tolist(), "se_re": s["se_re"].tolist(),
                            "z_re": s["z_re"].tolist(),
                            "ort_re": float(np.mean(s["re"])),
                            "ort_se": float(np.sqrt(np.mean(s["se_re"] ** 2)) / np.sqrt(len(sec)))}
        log(f"  {ad:>10}: boş konum Re K̃ ort {np.mean(s['re']):+.5f}  aralık "
            f"[{np.min(s['re']):+.5f}, {np.max(s['re']):+.5f}]  z: "
            + " ".join(f"{z:+.1f}" for z in s["z_re"]))
    out["zaman"] = time.ctime()
    json.dump(out, open(o.S195 / "K_kesif_195.json", "w"), indent=1, ensure_ascii=False)
    log(f"-> K_kesif_195.json ({time.time()-T0:.0f}s)")
