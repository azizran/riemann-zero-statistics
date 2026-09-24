# -*- coding: utf-8 -*-
"""
195c — B: ÇEKİRDEK K̃ (adalar; genelleştirilmiş makine 195o, K0c'de zeta'da kapılandı)
=====================================================================================
ONKAYIT_195 dondurdu (sha denetimi). Her ada, her üst-bölge bloğu b (8.5-çapalı
ΔL=0.05; ön-kayıt listesiyle assert):
  pencere çizgileri τ = log q/L_b ∈ [0.45, 0.86), χ(q) ≠ 0, a^χ = χ(q)a_q;
  c^öz_{q,b} (tam blok); her uydu: |log q' − L_b − log n| < 0.03, χ(q') ≠ 0 →
  seri → izdüşüm c^{(uydu)}_{q,b}; num_b = Σ_q c^{(uydu)} conj(c^öz), den_b = Σ|c^öz|².
K̃ = Σ N_b num_b / Σ N_b den_b (üst bölge HAVUZ; bant [10.0, U_χ) HAVUZ), blok-loo jk.
Pencere-içi uydu çizgisi sayımı (−log5 dışlaması denetimi).
Çıktı: scratchpad/195/B_cekirdek.json, B_bloklar.npz
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
    ad, j, idx, Lb = arg
    D = _W[ad]
    t0 = time.time()
    r = o.cekirdek_blok(D["g"], D["mid"], idx, idx, Lb, D["Q"], D["A"], D["CHI"],
                        o.UYDU_SIRA)
    lq = np.log(D["Q"])
    ic = []
    for u in o.UYDU_SIRA:
        sel = (np.abs(lq - Lb - o.log_n(u)) < o.DELTA) & (D["CHI"] != 0)
        ic.append(int(np.sum(sel & (lq / Lb < o.TAU_HI))))
    return dict(ad=ad, j=j, num=r["num"], den=r["den"], N=r["N_tam"],
                ncz_pencere=r["ncz_pencere"], ncz_uydu=r["ncz_uydu"],
                ncz_ic=np.array(ic), Lb=Lb, sure=time.time() - t0)


def ozetle(num, den, Nw):
    K, reps = o.havuz(num, den, Nw)
    s = o.ozet(K, reps)
    return {u: {"re": float(s["re"][i]), "im": float(s["im"][i]),
                "se_re": float(s["se_re"][i]), "se_im": float(s["se_im"][i]),
                "z_re": float(s["z_re"][i]), "z_im": float(s["z_im"][i]),
                "aci": float(s["aci"][i])} for i, u in enumerate(o.UYDU_SIRA)}


if __name__ == "__main__":
    T0 = time.time()
    ONK = onkayit()
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    log("=" * 78)
    log(f"195c / B ÇEKİRDEK K̃ (adalar)  [on-kayit {ONK['zaman']} sha {ONK['sha256'][:12]}]")
    log("=" * 78)
    tab, _ = o.karakter_tablolari()
    isler = []
    Lmax_all = 0.0
    kin = {}
    for ad in o.ADALAR:
        K = o.ada_kin(ad)
        kin[ad] = K
        Lmax_all = max(Lmax_all, float(K["Lm"].max()))
    qmax = np.exp(Lmax_all + np.log(10) + o.DELTA + 0.01)
    t0 = time.time()
    Q, LAM, A = o.asal_kuvvetler(qmax)
    log(f"  asal-kuvvetler q ≤ {qmax:.0f}: {len(Q)} ({time.time()-t0:.0f}s)")
    for ad in o.ADALAR:
        K = kin[ad]
        k = K["k"]
        CHI = o.chi_dizi(tab[ad], k, Q.astype(np.int64))
        _W[ad] = dict(g=K["g"], mid=K["mid"], Q=Q, A=A, CHI=CHI)
        bl = o.dL_bloklari(K["Lm"], o.L_UST)
        B = ONK["ada_bloklari"][ad]
        assert [x[0] for x in bl] == B["j"] and [len(x[3]) for x in bl] == B["N_b"], ad
        for (j, lo, hi, idx) in bl:
            Lb = float(K["Lm"][idx].mean())
            assert abs(Lb - B["L_b"][B["j"].index(j)]) < 1e-12
            isler.append((ad, j, idx, Lb))
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    R = {ad: {} for ad in o.ADALAR}
    with ctx.Pool(nw) as pool:
        for r in pool.imap_unordered(_is, sorted(isler, key=lambda x: -x[3])):
            R[r["ad"]][r["j"]] = r
    log(f"  {len(isler)} blok bitti ({time.time()-T0:.0f}s)")

    out = {"sha_onkayit": ONK["sha256"], "adalar": {}}
    npz = {}
    for ad in o.ADALAR:
        B = ONK["ada_bloklari"][ad]
        js = B["j"]
        num = np.array([R[ad][j]["num"] for j in js])
        den = np.array([R[ad][j]["den"] for j in js])
        Nw = np.array([R[ad][j]["N"] for j in js])
        ncz = np.array([R[ad][j]["ncz_uydu"] for j in js])
        ic = np.array([R[ad][j]["ncz_ic"] for j in js])
        ust = ozetle(num, den, Nw)
        bm = np.array([j in B["bant_bloklari_j"] for j in js])
        bant = ozetle(num[bm], den[bm], Nw[bm])
        for i, u in enumerate(o.UYDU_SIRA):
            ust[u]["n_cizgi"] = int(ncz[:, i].sum())
            ust[u]["n_pencere_ici"] = int(ic[:, i].sum())
            bant[u]["n_cizgi"] = int(ncz[bm, i].sum())
        out["adalar"][ad] = {"k": ONK["K0b"]["adalar"][ad]["k"],
                             "parite": ONK["K0b"]["adalar"][ad]["parite"],
                             "n_blok": len(js), "N": int(Nw.sum()),
                             "n_blok_bant": int(bm.sum()), "N_bant": int(Nw[bm].sum()),
                             "ust": ust, "bant": bant,
                             "pencere_cizgi_blok_ort": float(np.mean([R[ad][j]["ncz_pencere"] for j in js]))}
        npz[f"{ad}_num"] = num
        npz[f"{ad}_den"] = den
        npz[f"{ad}_N"] = Nw
        npz[f"{ad}_Lb"] = np.array([R[ad][j]["Lb"] for j in js])
        npz[f"{ad}_ncz"] = ncz
        log(f"\n  [{ad}] k={out['adalar'][ad]['k']} parite={out['adalar'][ad]['parite']} "
            f"blok={len(js)} N={int(Nw.sum())}  (bant: {int(bm.sum())} blok, N={int(Nw[bm].sum())})")
        for u in o.UYDU_SIRA:
            v, w = ust[u], bant[u]
            ys = "YASAK " if o.yasak_mi(u, out['adalar'][ad]['k']) else "izinli"
            log(f"     {u:>10} {ys}: ÜST Re={v['re']:+.5f}±{v['se_re']:.5f} (z={v['z_re']:+6.2f}) "
                f"Im={v['im']:+.5f} (z={v['z_im']:+5.2f}) ∠{v['aci']:+6.1f}°  çizgi={v['n_cizgi']} "
                f"p-içi={v['n_pencere_ici']} | BANT Re={w['re']:+.5f}±{w['se_re']:.5f} (z={w['z_re']:+6.2f})")
    np.savez_compressed(o.S195 / "B_bloklar.npz", uydular=np.array(o.UYDU_SIRA), **npz)
    out["zaman"] = time.ctime()
    json.dump(out, open(o.S195 / "B_cekirdek.json", "w"), indent=1, ensure_ascii=False)
    log(f"\n-> B_cekirdek.json, B_bloklar.npz ({time.time()-T0:.0f}s)")
