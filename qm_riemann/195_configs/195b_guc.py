# -*- coding: utf-8 -*-
"""
195b — A: YAPI ÇARPANI GÜCÜ |Ĝ_b(Δω)|² (adalar + zeta kontrolleri)
=================================================================
ONKAYIT_195 dondurdu (sha denetimi). Her blok b: Ĝ_b(Δω) = mean e^{i(L_b+Δω)m_n}
(MUTLAK m_n), Δω ∈ [−2.0, 2.6] adım 0.001. Uydu n: pencere |Δω − log n| < 0.03,
halka 0.06–0.25 (a,b ≤ 6 rasyonel komşuları ±0.035 dışlanmış).
D_b = mean_pencere|Ĝ_b|² − mean_halka|Ĝ_b|²; D = Σ N_b D_b/Σ N_b; z = D/se_jk(D)
(blok-loo). Adalar: 8.5-çapalı ΔL=0.05 blokları (ön-kayıt listesiyle assert).
Zeta kontrolleri: son (193 blokları), düşük (ΔL ızgarası).
Çıktı: scratchpad/195/A_guc.json, A_profiller.npz
"""
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("o195", HERE / "195o_ortak.py")
o = importlib.util.module_from_spec(spec)
spec.loader.exec_module(o)


def log(*a):
    print(*a, flush=True)


def onkayit():
    ONK = json.load(open(o.S195 / "ONKAYIT_195.json"))
    s = hashlib.sha256((HERE / "195a_onkayit.py").read_bytes()).hexdigest()
    if s != ONK["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {s} != {ONK['sha256']}")
    return ONK


_W = {}


def _is(arg):
    ad, et, idx, Lb = arg
    mid = _W[ad]
    return ad, et, o.G_guc_blok(mid[idx], Lb), len(idx), Lb


def blok_listesi(ad, ONK):
    if ad in o.ADALAR:
        K = o.ada_kin(ad)
        bl = o.dL_bloklari(K["Lm"], o.L_UST)
        B = ONK["ada_bloklari"][ad]
        assert [x[0] for x in bl] == B["j"] and [len(x[3]) for x in bl] == B["N_b"], ad
        return K["mid"], [(f"j{j}", idx, float(K["Lm"][idx].mean())) for (j, lo, hi, idx) in bl]
    if ad == "zeta_son":
        K = o.zeta_kin("son")
        Lb193 = json.load(open(o.S193 / "ONKAYIT_193.json"))["bloklar"]["L_b"]
        out = []
        for b, idx in enumerate(o.esit_sayim_bloklari(len(K["mid"]))):
            Lb = float(np.log(K["mid"][idx] / o.TWO_PI).mean())
            assert abs(Lb - Lb193[b]) < 1e-12
            out.append((f"b{b}", idx, Lb))
        return K["mid"], out
    K = o.zeta_kin("dusuk")
    return K["mid"], [(f"j{j}", idx, float(K["Lm"][idx].mean()))
                      for (j, lo, hi, idx) in o.dL_bloklari(K["Lm"], o.L_UST)]


if __name__ == "__main__":
    T0 = time.time()
    ONK = onkayit()
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    log("=" * 78)
    log(f"195b / A YAPI ÇARPANI GÜCÜ  [on-kayit {ONK['zaman']} sha {ONK['sha256'][:12]}]")
    log("=" * 78)
    kumeler = o.ADALAR + ["zeta_son", "zeta_dusuk"]
    isler, meta = [], {}
    for ad in kumeler:
        mid, bl = blok_listesi(ad, ONK)
        _W[ad] = mid
        meta[ad] = [x[0] for x in bl]
        for (et, idx, Lb) in bl:
            isler.append((ad, et, idx, Lb))
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    guc = {ad: {} for ad in kumeler}
    with ctx.Pool(nw) as pool:
        for ad, et, P, n, Lb in pool.imap_unordered(_is, sorted(isler, key=lambda x: -len(x[2]))):
            guc[ad][et] = (P, n, Lb)
    log(f"  |Ĝ|² hesapları bitti ({time.time()-T0:.0f}s; {len(isler)} blok)")

    maskeler = {u: o.pencere_halka_maskeleri(u) for u in o.A_LISTE + o.A_KAYIT}
    out = {"sha_onkayit": ONK["sha256"], "sonuc": {}}
    prof = {}
    for ad in kumeler:
        ets = meta[ad]
        P = np.array([guc[ad][e][0] for e in ets])
        Nb = np.array([guc[ad][e][1] for e in ets], float)
        Lb = np.array([guc[ad][e][2] for e in ets])
        prof[f"{ad}_P"] = (Nb[:, None] * P).sum(0) / Nb.sum()        # ağırlıklı |Ĝ|²
        prof[f"{ad}_Pn"] = (Nb[:, None] * (Nb[:, None] * P)).sum(0) / Nb.sum()  # gürültü tabanı birimi
        prof[f"{ad}_Lb"] = Lb
        prof[f"{ad}_Nb"] = Nb
        k = {"zeta_son": 1, "zeta_dusuk": 1}.get(ad) or ONK["K0b"]["adalar"][ad]["k"]
        res = {}
        for u in o.A_LISTE + o.A_KAYIT:
            W, R = maskeler[u]
            Db = P[:, W].mean(1) - P[:, R].mean(1)
            D = float((Nb * Db).sum() / Nb.sum())
            reps = np.array([((Nb * Db).sum() - Nb[j] * Db[j]) / (Nb.sum() - Nb[j])
                             for j in range(len(Nb))])
            se = float(o.jk_se(reps))
            # gürültü-tabanı biriminde (bilgi): N_b·D_b
            Dn = float((Nb * (Nb * Db)).sum() / Nb.sum())
            res[u] = {"D": D, "se": se, "z": D / se, "D_taban_biriminde": Dn,
                      "yasak": bool(o.yasak_mi(u, k)) if k > 1 else False,
                      "pencere_ort": float((Nb * P[:, W].mean(1)).sum() / Nb.sum()),
                      "halka_ort": float((Nb * P[:, R].mean(1)).sum() / Nb.sum())}
        out["sonuc"][ad] = {"k": k, "n_blok": len(ets), "N": int(Nb.sum()), "uydular": res}
        log(f"\n  [{ad}] k={k} blok={len(ets)} N={int(Nb.sum())}")
        for u in o.A_LISTE + o.A_KAYIT:
            r = res[u]
            log(f"     {u:>10} {'YASAK ' if r['yasak'] else 'izinli'}: D={r['D']:+.3e}±{r['se']:.2e}"
                f"  z={r['z']:+7.2f}  (N_b·D={r['D_taban_biriminde']:+.2f} taban)")
    np.savez_compressed(o.S195 / "A_profiller.npz", dw=o.DW_IZGARA, **prof)
    out["zaman"] = time.ctime()
    json.dump(out, open(o.S195 / "A_guc.json", "w"), indent=1, ensure_ascii=False)
    log(f"\n-> A_guc.json, A_profiller.npz ({time.time()-T0:.0f}s)")
