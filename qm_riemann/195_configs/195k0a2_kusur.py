# -*- coding: utf-8 -*-
"""
195k0a2 — K0a EKİ: dosya bütünlüğü (tamsayı indeks sıçramaları) + düzeltilmiş kalibrasyon
======================================================================================
195k0ab'de χ₃ ve χ₈ₒ için ortalama-yuvarlama j0 (116/147) ile bağımsız sayım ve
tek-nokta (114/145) 2 farklı çıktı. Bu betik:
  (1) her adada d_i = N̄_χ(m_i) − i − j0_bağımsız dizisinin TAMSAYI sıçramalarını
      bulur (yerel medyan, yuvarlanmış) — eksik sıfır = +1 sıçrama;
  (2) sıçramayı taşıyan aralıkta mpmath Hardy-Z_χ işaret değişimlerini sayar
      (adım 0.002) — eksik sıfırların GERÇEK konumları;
  (3) kusur düzeltilmiş indeksle (sıçramadan sonra n += eksik sayısı)
      ortalama[N̄_χ(m_n) − n] (TÜM ve ÜST L≥8.5) — C_χ kapısı (|·|<0.05);
  (4) LİTERAL sayılar (tek-nokta j0, düzeltmesiz) de KAYIT olarak yazılır.
Not: A ve B ölçümleri n'yi KULLANMAZ (yalnız mutlak m_n, g_n); kusur(lar)
üst bölgenin (L ≥ 8.5) dışındaysa ölçümler etkilenmez.
Çıktı: scratchpad/195/K0a_kusur.json
"""
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location(
    "o195", Path(__file__).resolve().parent / "195o_ortak.py")
o = importlib.util.module_from_spec(spec)
spec.loader.exec_module(o)


def log(*a):
    print(*a, flush=True)


if __name__ == "__main__":
    import mpmath as mp
    mp.mp.dps = 20
    K0 = json.load(open(o.S195 / "K0ab.json"))
    tab, _ = o.karakter_tablolari()
    out = {}
    gecti = True
    for ad in o.ADALAR:
        K = o.ada_kin(ad)
        z, mid, k, par = K["z"], K["mid"], K["k"], K["a"]
        C = o.C_chi(par)
        j0b = K0["K0a"][ad]["j0_bagimsiz_mpmath"]
        i = np.arange(1, len(mid) + 1)
        d = o.Nbar(mid, k, C) - i - j0b
        # yerel medyan (41 nokta) → yuvarlanmış tamsayı seviye
        w = 41
        pad = np.pad(d, w // 2, mode="edge")
        from numpy.lib.stride_tricks import sliding_window_view
        med = np.median(sliding_window_view(pad, w), axis=1)
        seviye = np.round(med).astype(int)
        sic = np.where(np.diff(seviye) != 0)[0]
        kusurlar = []
        chi = [tab[ad][r] for r in range(k)]

        def Z(t):
            s = mp.mpc(0.5, t)
            th = t / 2 * mp.log(mp.mpf(k) / mp.pi) + mp.im(mp.loggamma((s + par) / 2))
            return float(mp.re(mp.exp(1j * th) * mp.dirichlet(s, chi)))

        gorulen = set()
        for s_ in sic:
            # sıçramayı taşıyan en büyük boşluk (±w/2 içinde); aynı boşluk bir kez
            lo_, hi_ = max(0, s_ - w // 2), min(len(mid) - 1, s_ + w // 2)
            gix = lo_ + int(np.argmax(np.diff(z)[lo_:hi_ + 1]))
            if gix in gorulen:
                continue
            gorulen.add(gix)
            a_, b_ = z[gix], z[gix + 1]
            ts = np.arange(a_ - 0.3, b_ + 0.3, 0.002)
            v = np.array([Z(float(t)) for t in ts])
            ch = ts[np.where(np.sign(v[:-1]) * np.sign(v[1:]) < 0)[0]] + 0.001
            # dosya sıfırlarıyla eşleşmeyen (±0.03) mpmath işaret değişimleri = eksik
            zyakin = z[(z > ts[0]) & (z < ts[-1])]
            ic = np.array([c for c in ch if np.min(np.abs(zyakin - c)) > 0.03])
            assert len(ch) - len(ic) == len(zyakin), (ad, "eşleşme sayısı")
            kusurlar.append({"orta_nokta_idx": int(gix), "z_sol": float(a_), "z_sag": float(b_),
                             "bosluk": float(b_ - a_),
                             "L": float(np.log(k * a_ / o.TWO_PI)),
                             "seviye_sicrama": int(seviye[s_ + 1] - seviye[s_]),
                             "eksik_sifirlar_mpmath": [float(x) for x in ic]})
        # düzeltilmiş indeks
        n = j0b + i.copy()
        for kk in kusurlar:
            n[kk["orta_nokta_idx"] + 1:] += len(kk["eksik_sifirlar_mpmath"])
        gecerli = np.ones(len(mid), bool)
        for kk in kusurlar:        # kusurlu boşluğun "orta noktası" gerçek orta nokta değil
            gecerli[kk["orta_nokta_idx"]] = False
        r_all = float(np.mean((o.Nbar(mid, k, C) - n)[gecerli]))
        ust = (K["Lm"] >= o.L_UST) & gecerli
        r_ust = float(np.mean((o.Nbar(mid, k, C) - n)[ust]))
        lit_j0 = K0["K0a"][ad]["j0_tek_nokta_Nbar_z1"]
        r_lit = float(np.mean(o.Nbar(mid, k, C) - (lit_j0 + i)))
        ok = abs(r_all) < 0.05 and abs(r_ust) < 0.05
        gecti &= ok
        kusur_ust_bolgede = any(kk["L"] >= o.L_UST for kk in kusurlar)
        out[ad] = {"j0_bagimsiz": j0b, "kusurlar": kusurlar,
                   "duzeltilmis_ort_TUM": r_all, "duzeltilmis_ort_UST": r_ust,
                   "literal_tek_nokta_j0": lit_j0, "literal_ort_TUM_duzeltmesiz": r_lit,
                   "kusur_ust_bolgede_mi": kusur_ust_bolgede, "kapi_0.05": ok}
        log(f"  {ad:>6}: j0(bağımsız)={j0b}; kusur sayısı={len(kusurlar)} "
            f"{[(round(kk['z_sol'],3), round(kk['z_sag'],3), kk['eksik_sifirlar_mpmath']) for kk in kusurlar]}")
        log(f"          düzeltilmiş mean[N̄−n]: TÜM={r_all:+.5f} ÜST={r_ust:+.5f} "
            f"(literal düzeltmesiz: {r_lit:+.4f}) → {'GEÇTİ' if ok else 'KALDI'}; "
            f"kusur üst bölgede mi: {kusur_ust_bolgede}")
    out["K0a_duzeltilmis_GECTI"] = bool(gecti)
    out["zaman"] = time.ctime()
    json.dump(out, open(o.S195 / "K0a_kusur.json", "w"), indent=1, ensure_ascii=False)
    log(f"-> K0a_kusur.json; düzeltilmiş K0a: {'GEÇTİ' if gecti else 'KALDI'}")
