# -*- coding: utf-8 -*-
"""
192b — K1: B SEÇİM KURALI KATALOĞU (mevcut Δω profillerinden; yeni harita YOK)
=============================================================================
ONKAYIT_192 dondurdu. İki profil:
 BİRİNCİL — DÜŞÜK pencere ω-dilim (0.025) havuz profili: 190b D çıktıları
   (scratchpad/190/omega/b<b>_c<k>.npz; blok b, dilim j re/im blok toplamları).
   T_{b,j}^{(v)} = Σ_{q∈HAVUZ} 2(re+i·im)_{b,j}(q)·conj(mix_q^{(v)}) / Σ|mix^{(v)}|²,
   K_j^{(v)} = Σ_{b∈B_j^{(v)}} T_{b,j}^{(v)} / Σ_{b∈B_j^{(v)}} n_b; κ = −Re K.
   v = tam / loo-i (B_j^{(i)} = B_j \\ {i}; mix loo — 188b K_hesap AYNEN).
   Mühür: T_{b,j}^{tam}/n_b ≡ 190 harita_omega_dusuk K[b, HAVUZ, j].
 İKİNCİL (KAYIT) — SON penceresi: 188/harita_proj_gercek.npz blok τ'-dilim
   κ_b(s) (190c tau_blok_profili AYNEN; mix loo eklenir), yoğunluk κ_b/(0.005·L_son)
   Δω_b ekseninden 0.025 ızgarasına np.interp (190e son_yog AYNEN), kapsayan
   bloklarda n_b-ağırlıklı ortalama. Mühür: κ_b^{tam} ≡ profiller_190 ts_kap.
Her katalog konumu: κ_tepe (±0.03 arama), jk se, yerel taban (0.05<|Δ|≤0.15
halkası; medyan, 1.4826·MAD), hüküm yanar / söner / belirsiz (+ KAYIT'lar).
Çıktı: 192/K1_katalog.json, 192/profiller_192.npz + ekran.
"""
import hashlib
import importlib.util
import json
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
           "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S = {d: SCR / d for d in ["155", "184", "185", "186", "188", "190", "192"]}
ZD = S["190"] / "zincir_dusuk"
OM = S["190"] / "omega"
TWO_PI = 2 * np.pi
NJACK = 8
DW = 0.025
EPS = 1e-9


def yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec)
    sys.modules[ad] = m
    spec.loader.exec_module(m)
    return m


b188 = yukle("b188", QM / "188_configs" / "188b_harita.py")


def onkayit():
    o = json.load(open(S["192"] / "ONKAYIT_192.json"))
    s = hashlib.sha256((QM / "192_configs" / "192a_onkayit.py").read_bytes()).hexdigest()
    if s != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {s} != {o['sha256']}")
    return o


def jk(r):
    r = np.asarray(r, float)
    return float(np.sqrt((NJACK - 1) / NJACK * np.sum((r - r.mean()) ** 2)))


def mixler(K1, OZ, P):
    tau = np.asarray(K1["tau"], float)
    win = (tau >= 0.45) & (tau < 0.86)
    aq = np.asarray(K1["aq"], float)
    tam = b188.karisim(P, OZ, aq)[win]
    reps = [b188.karisim(P, OZ, aq, j)[win] for j in range(NJACK)]
    return tam, reps, tau[win]


# ---------------------------------------------------------------------------
def dusuk_profil(ONK, ONK188):
    L = float(ONK["kinematik"]["L_dusuk"])
    Lb = np.array(ONK["kinematik"]["L_b_dusuk"])
    K1 = np.load(ZD / "K1_gercek_dusuk.npz")
    mix_t, mix_r, tau_win = mixler(K1, np.load(ZD / "OZ_gercek_dusuk.npz"),
                                   np.load(ZD / "G1_proj_gercek_dusuk.npz"))
    mask, et = b188.bant_maskeleri(tau_win, ONK188)
    hv = mask[et.index("HAVUZ")]
    H = np.load(S["190"] / "harita_omega_dusuk.npz")
    J = H["J"]
    jmin = int(J.min())
    nJ = len(J)
    T = np.zeros((NJACK + 1, NJACK, nJ), complex)     # v=0 tam, v=1..8 loo(i=v−1)
    nb = np.zeros(NJACK, int)
    mixes = [mix_t] + mix_r
    for b in range(NJACK):
        parca = sorted(OM.glob(f"b{b}_c*.npz"), key=lambda p: int(p.stem.split("_c")[1]))
        for p in parca:
            d = np.load(p)
            C = 2.0 * (d["re"] + 1j * d["im"])          # (nj, nwin), /n_b YOK
            ix = d["j"] - jmin
            for v, mx in enumerate(mixes):
                T[v, b, ix] = (C[:, hv] @ np.conj(mx[hv])) / np.sum(np.abs(mx[hv]) ** 2)
            nb[b] = int(d["nb"])
    assert np.array_equal(nb, H["nb"])
    iH = list(H["etiket"]).index("HAVUZ")
    muhur = float(np.max(np.abs(T[0] / nb[:, None] - H["K"][:, iH, :])))
    merkez = J * DW
    kap = np.zeros((NJACK, nJ), bool)
    for b in range(NJACK):
        kap[b] = ((merkez - 0.5 * DW >= 0.86 * L - Lb[b] - EPS) &
                  (merkez + 0.5 * DW <= 1.30 * L - Lb[b] + EPS))
    Kt = np.zeros(nJ, complex)
    Kr = np.zeros((NJACK, nJ), complex)
    Kt_sum = (T[0] * kap).sum(0)
    nsum = (nb[:, None] * kap).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        Kt = Kt_sum / nsum
        for i in range(NJACK):
            kk = kap.copy()
            kk[i] = False
            Kr[i] = (T[1 + i] * kk).sum(0) / (nb[:, None] * kk).sum(0)
    return {"merkez": merkez, "J": J, "kap": -Kt.real, "kap_reps": -Kr.real,
            "K": Kt, "K_reps": Kr, "nkap": kap.sum(0), "kapsama": kap, "nb": nb,
            "muhur_K_om": muhur, "ncz": H["ncz"]}


def son_profil(ONK, ONK188, merkez):
    Ls = float(ONK["kinematik"]["L_son"])
    Lbs = np.array(ONK["kinematik"]["L_b_son"])
    K1 = np.load(S["184"] / "K1_gercek.npz")
    mix_t, mix_r, tau_win = mixler(K1, np.load(S["185"] / "OZ_gercek.npz"),
                                   np.load(S["186"] / "G1_proj_gercek.npz"))
    PR = np.load(S["188"] / "harita_proj_gercek.npz")
    k = PR["kenar"]
    o = 0.5 * (k[:-1] + k[1:])
    nb = PR["nb"]
    mixes = [mix_t] + mix_r
    kapb = np.zeros((NJACK + 1, NJACK, len(o)))
    dw = np.zeros((NJACK, len(o)))
    for b in range(NJACK):
        C = 2 * (PR["re"][b] + 1j * PR["im"][b]) / nb[b]
        dw[b] = o * Ls - Lbs[b]
        for v, mx in enumerate(mixes):
            kapb[v, b] = -(C @ np.conj(mx)).real / np.sum(np.abs(mx) ** 2)
    P190 = np.load(S["190"] / "profiller_190.npz")
    muhur = float(np.max(np.abs(kapb[0] - P190["ts_kap"])))
    muhur_dw = float(np.max(np.abs(dw - P190["ts_dw"])))
    kap = np.zeros((NJACK, len(merkez)), bool)
    for b in range(NJACK):
        kap[b] = ((merkez - 0.5 * DW >= o[0] * Ls - Lbs[b] - EPS) &
                  (merkez + 0.5 * DW <= o[-1] * Ls - Lbs[b] + EPS))
    Y = np.zeros((NJACK + 1, NJACK, len(merkez)))
    for v in range(NJACK + 1):
        for b in range(NJACK):
            Y[v, b] = np.interp(merkez, dw[b], kapb[v, b] / (0.005 * Ls),
                                left=0.0, right=0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        yt = (nb[:, None] * kap * Y[0]).sum(0) / (nb[:, None] * kap).sum(0)
        yr = np.zeros((NJACK, len(merkez)))
        for i in range(NJACK):
            kk = kap.copy()
            kk[i] = False
            yr[i] = (nb[:, None] * kk * Y[1 + i]).sum(0) / (nb[:, None] * kk).sum(0)
    return {"kap": yt, "kap_reps": yr, "nkap": kap.sum(0), "nb": nb,
            "muhur_ts_kap": muhur, "muhur_ts_dw": muhur_dw, "kap_blok": kapb[0], "dw_blok": dw}


# ---------------------------------------------------------------------------
def degerlendir(ad, kat, prof, J, katalog, profil_adi):
    kapA = prof["kap"]
    rep = prof["kap_reps"]
    nk = prof["nkap"]
    jmin = int(J[0])
    ara = [j for j in kat["arama_j"]]
    if any(nk[j - jmin] < 4 for j in ara):
        return {"hukum": "erişilemedi", "neden": "arama diliminde kapsama < 4"}
    vals = np.array([kapA[j - jmin] for j in ara])
    js = ara[int(np.argmax(vals))]
    kt = float(kapA[js - jmin])
    se = jk(rep[:, js - jmin])
    hal = [j for j in kat["halka_j"] if nk[j - jmin] >= 4]
    hv = np.array([kapA[j - jmin] for j in hal])
    med = float(np.median(hv))
    yay = float(1.4826 * np.median(np.abs(hv - med)))
    ust = med + 2.0 * yay
    yanar = (kt > 3 * se) and (kt > ust)
    soner = (abs(kt) < 2 * se) or (abs(kt - med) <= 2.0 * yay)
    cukur = (kt < med - 2.0 * yay) and (abs(kt) >= 2 * se)
    hk = "yanar" if yanar else ("söner" if soner else "belirsiz")
    # KAYIT: çakışmasız arama
    ev = {katalog[x]["ev_dilimi_j"] for x in katalog if x != ad}
    ara2 = [j for j in ara if j not in ev]
    if ara2:
        v2 = np.array([kapA[j - jmin] for j in ara2])
        js2 = ara2[int(np.argmax(v2))]
        kt2 = float(kapA[js2 - jmin])
        se2 = jk(rep[:, js2 - jmin])
        y2 = (kt2 > 3 * se2) and (kt2 > ust)
        s2 = (abs(kt2) < 2 * se2) or (abs(kt2 - med) <= 2.0 * yay)
        cak = {"j": js2, "merkez": js2 * DW, "kappa": kt2, "se": se2,
               "hukum": "yanar" if y2 else ("söner" if s2 else "belirsiz"),
               "yanar_4se": bool((kt2 > 4 * se2) and (kt2 > ust))}
    else:
        cak = None
    return {"profil": profil_adi, "tepe_j": js, "tepe_merkez": js * DW,
            "kappa_tepe": kt, "se": se, "kappa_bolu_se": kt / se if se > 0 else None,
            "arama_min_kappa": float(vals.min()),
            "taban_medyan": med, "taban_yayilim": yay, "taban_ust": ust,
            "halka_n": len(hal), "yanar": bool(yanar), "soner": bool(soner),
            "cukur": bool(cukur), "hukum": hk,
            "yanar_4se": bool((kt > 4 * se) and (kt > ust)),
            "yalin_4se": bool(kt > 4 * se),
            "cakismasiz_KAYIT": cak if set(ara2) != set(ara) else "aynı (çakışma yok)"}


if __name__ == "__main__":
    ONK = onkayit()
    ONK188 = json.load(open(S["188"] / "ONKAYIT_188.json"))
    KAT = ONK["B"]["katalog"]
    print("=" * 78)
    print(f"192b / K1 B KATALOĞU  [on-kayit {ONK['zaman']} sha {ONK['sha256'][:12]}]")
    print("=" * 78)
    D = dusuk_profil(ONK, ONK188)
    print(f"  düşük: {len(D['J'])} dilim; mühür max|T/n_b − K_om(HAVUZ)| = "
          f"{D['muhur_K_om']:.1e}; nb = {D['nb'].tolist()}")
    # kapsama ön-kayıtla aynı mı
    jmin = int(D["J"][0])
    for ad, k in KAT.items():
        for js, bl in k["kapsama_dusuk"].items():
            assert list(np.where(D["kapsama"][:, int(js) - jmin])[0]) == bl, (ad, js)
    Sn = son_profil(ONK, ONK188, D["merkez"])
    print(f"  son: mühür max|κ_b − ts_kap(190)| = {Sn['muhur_ts_kap']:.1e}, "
          f"max|Δω_b − ts_dw| = {Sn['muhur_ts_dw']:.1e}")
    for ad, k in KAT.items():
        for js, bl in k["kapsama_son"].items():
            assert int(Sn["nkap"][int(js) - jmin]) == len(bl), (ad, js)

    out = {"sha_onkayit": ONK["sha256"], "muhur": {
        "dusuk_T_bolu_nb_vs_K_om": D["muhur_K_om"],
        "son_kappa_b_vs_ts_kap": Sn["muhur_ts_kap"],
        "son_dw_vs_ts_dw": Sn["muhur_ts_dw"]}, "dusuk": {}, "son": {}}
    for pad, prof in (("dusuk", D), ("son", Sn)):
        print(f"\n  [{pad.upper()}] {'hedef':>10} {'liste':>7} {'Δω':>7} | tepe  "
              f"{'κ_tepe':>10} {'jk-se':>9} {'κ/se':>6} | taban med±yay       | hüküm")
        for ad, k in KAT.items():
            r = degerlendir(ad, k, prof, D["J"], KAT, pad)
            out[pad][ad] = r
            if r["hukum"] == "erişilemedi":
                print(f"   {ad:>10} {k['liste']:>7} {k['delta_omega']:+.3f} | ERİŞİLEMEDİ")
                continue
            ek = (" ÇUKUR" if r["cukur"] else "") + (" [4se-YANAR]" if r["yanar_4se"] else "")
            cak = r["cakismasiz_KAYIT"]
            if isinstance(cak, dict):
                ek += f"  (çakışmasız: {cak['merkez']:.3f} κ={cak['kappa']:.2e}±{cak['se']:.1e} → {cak['hukum']})"
            print(f"   {ad:>10} {k['liste']:>7} {k['delta_omega']:+.3f} | {r['tepe_merkez']:+.3f} "
                  f"{r['kappa_tepe']:+.3e} {r['se']:.2e} {r['kappa_bolu_se']:+6.1f} | "
                  f"{r['taban_medyan']:+.2e}±{r['taban_yayilim']:.1e} | {r['hukum']}{ek}")
    np.savez_compressed(
        S["192"] / "profiller_192.npz", merkez=D["merkez"], J=D["J"],
        d_kap=D["kap"], d_kap_reps=D["kap_reps"], d_nkap=D["nkap"], d_nb=D["nb"],
        d_K=D["K"], d_K_reps=D["K_reps"],
        s_kap=Sn["kap"], s_kap_reps=Sn["kap_reps"], s_nkap=Sn["nkap"], s_nb=Sn["nb"])
    json.dump(out, open(S["192"] / "K1_katalog.json", "w"), indent=1, ensure_ascii=False)
    print(f"\n-> K1_katalog.json, profiller_192.npz  BİTTİ")
