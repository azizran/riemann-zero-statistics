# -*- coding: utf-8 -*-
"""
187g — İKİNCİL: DERİN SONDA (1.20,1.30] — 100k alt-örneklem
============================================================
ONKAYIT_187: "gerçek, (1.20,1.30] katmanı 100k alt-örneklem (her 3. nokta
n%3==0 — pencere boyu korunur), ayrı se, İKİNCİL damgası; ana hükümlere
GİRMEZ." Karşılaştırılabilirlik için m^ya(1.20) ve m^ya(1.30) AYNI
alt-örneklem üzerinde izdüşürülür (kayıtlı katman serileri yeniden
kullanılır; yalnız yeni katman serisi alt-örneklemde inşa edilir).

Çıktı: 187/IKINCIL_derin_sonda.json
"""
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S155, S184, S185, S186, S187 = (SCR / d for d in
                                ["155", "184", "185", "186", "187"])
NJACK = 8
PBLOK = 2048
LBLOK = 4096

spec = importlib.util.spec_from_file_location(
    "b185", QM / "185_configs" / "185b_oz_muhasebe.py")
b185 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b185)
spec7 = importlib.util.spec_from_file_location(
    "b187", QM / "187_configs" / "187b_katman_defteri.py")
b187 = importlib.util.module_from_spec(spec7)
spec7.loader.exec_module(b187)

if __name__ == "__main__":
    ONK = b187.onkayit()
    print("=" * 78)
    print(f"187g / İKİNCİL DERİN SONDA (1.20,1.30]  [on-kayit sha "
          f"{ONK['sha256'][:12]}]  — ana hükümlere GİRMEZ")
    print("=" * 78, flush=True)
    L = ONK["L"]
    izg = ONK["katman_izgara_gercek"]
    NKAT = len(izg) - 1

    G = np.load(S184 / "K1_gercek.npz")
    OZg = np.load(S185 / "OZ_gercek.npz")
    w_win = np.asarray(G["w"], float)
    tau = np.asarray(G["tau"], float)
    aq = np.asarray(G["aq"], float)
    ae = np.asarray(G["aq_eff"], float)
    g, mid, Lg = b185.kinematik("gercek")
    N = len(mid)
    alt = np.arange(0, N, 3)          # n%3==0 → 100000 nokta
    kenar_jk = np.linspace(0, N, NJACK + 1).astype(int)

    # yeni katman çizgileri (1.20, 1.30]
    qmax = np.exp(1.30 * L)
    qa, lam, taua, aqa = b187.asal_kuvvetler(qmax, L)
    idx = np.where((taua > 1.20) & (taua <= 1.30))[0]
    print(f"  (1.20,1.30]: {len(idx)} çizgi (q≤{qmax:.0f})", flush=True)

    # yeni katman serisi (yalnız alt-örneklem noktalarında)
    yol = S187 / "katman_8sub_gercek.npy"
    if not yol.exists():
        dds8 = b187.katman_serisi(g[alt], mid[alt], np.log(qa[idx]),
                                  aqa[idx], "derin-sonda (alt-örneklem)")
        np.save(yol, dds8)
    dds8 = np.load(yol)

    # alt-örneklemde tüm serilerin izdüşümü (taban + 7 katman + yeni)
    seriler = [np.load(S187 / f"katman_{i}_gercek.npy")[alt]
               for i in range(NKAT + 1)] + [dds8]
    NS = len(seriler)
    re = np.zeros((NS, NJACK, len(w_win)))
    im = np.zeros((NS, NJACK, len(w_win)))
    nb_alt = np.zeros(NJACK, int)
    t0 = time.time()
    for b in range(NJACK):
        m_blk = alt[(alt >= kenar_jk[b]) & (alt < kenar_jk[b + 1])]
        nb_alt[b] = len(m_blk)
        for s0 in range(0, len(m_blk), 1500):
            sl = m_blk[s0:s0 + 1500]
            P = np.outer(mid[sl], w_win)
            cP = np.cos(P)
            sP = np.sin(P)
            del P
            j0 = np.searchsorted(alt, sl[0])
            j1 = j0 + len(sl)
            for i in range(NS):
                v = seriler[i][j0:j1]
                re[i, b] += v @ cP
                im[i, b] -= v @ sP
            del cP, sP
        print(f"  [izdüşüm-alt] jk-blok {b+1}/8  ({time.time()-t0:.0f}s)",
              flush=True)
    Nalt = int(nb_alt.sum())

    kenar_b = json.load(open(S186 / "ONKAYIT_186.json"))["kenar"]
    maskeler = [(tau >= kenar_b[b]) & (tau < kenar_b[b + 1]) for b in range(8)]
    maskeler.append((tau >= kenar_b[0]) & (tau < kenar_b[-1]))
    et = [f"{kenar_b[b]:.2f}-{kenar_b[b+1]:.2f}" for b in range(8)] + ["HAVUZ"]

    def cplx(i, disari):
        if disari < 0:
            return 2.0 * (re[i].sum(0) + 1j * im[i].sum(0)) / Nalt
        return 2.0 * ((re[i].sum(0) - re[i][disari]) +
                      1j * (im[i].sum(0) - im[i][disari])) / \
            (Nalt - nb_alt[disari])

    def defter(disari):
        a_oz = np.abs(b185.c_oz(OZg, aq, disari))
        cya20 = sum(cplx(i, disari) for i in range(NKAT + 1))
        cya30 = cya20 + cplx(NS - 1, disari)
        out = np.zeros((9, 2))
        for k, m in enumerate(maskeler):
            sae = ae[m].sum()
            soz = a_oz[m].sum()
            out[k, 0] = (np.abs(cya20[m]).sum() - soz) / sae
            out[k, 1] = (np.abs(cya30[m]).sum() - soz) / sae
        return out

    tam = defter(-1)
    reps = np.array([defter(i) for i in range(NJACK)])
    se = b185.jk_se(reps)
    dm = tam[:, 1] - tam[:, 0]
    se_dm = b185.jk_se(reps[:, :, 1] - reps[:, :, 0])

    K1 = json.load(open(S187 / "K1_katman_defteri.json"))
    m_olc = np.array(K1["m_olc"])
    print("\nİKİNCİL DEFTER (alt-örneklem; ±jk se) — ana hükümlere GİRMEZ:")
    print(f"{'bant':>10} {'m^ya_alt(1.20)':>15} {'m^ya_alt(1.30)':>15} "
          f"{'Δm':>18} {'m_ölç':>8}")
    for k in range(9):
        print(f"{et[k]:>10} {tam[k,0]:8.4f}±{se[k,0]:.4f} "
              f"{tam[k,1]:8.4f}±{se[k,1]:.4f} "
              f"{dm[k]:+.4f}±{se_dm[k]:.4f} {m_olc[k]:8.4f}")
    acik20 = tam[8, 0] - m_olc[8]
    acik30 = tam[8, 1] - m_olc[8]
    print(f"\n  HAVUZ açık: 1.20 → {acik20:+.4f}, 1.30 → {acik30:+.4f} "
          f"(Δm = {dm[8]:+.4f} ± {se_dm[8]:.4f})")

    sonuc = {"damga": "İKİNCİL — ana hükümlere girmez",
             "sha_onkayit": ONK["sha256"], "etiket": et,
             "cizgi_sayisi": int(len(idx)), "n_alt": Nalt,
             "m_ya_alt_120": tam[:, 0].tolist(),
             "m_ya_alt_130": tam[:, 1].tolist(),
             "se": se.tolist(), "dm": dm.tolist(),
             "se_dm": se_dm.tolist()}
    json.dump(sonuc, open(S187 / "IKINCIL_derin_sonda.json", "w"), indent=1,
              ensure_ascii=False)
    print(f"\n-> {S187/'IKINCIL_derin_sonda.json'}  BİTTİ", flush=True)
