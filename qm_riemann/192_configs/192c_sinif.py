# -*- coding: utf-8 -*-
"""
192c — K2: A1 KALINTI-SINIFI AYRIŞTIRMASI (SON penceresi; 188 makinesi AYNEN)
============================================================================
ONKAYIT_192 dondurdu. 188b_harita.py importlib ile yüklenir (DOSYA DÜZENLENMEZ);
seri_ve_G, izdusum, c_proj, karisim, K_matris, bant_maskeleri, jk_kenar AYNEN.

 M   MAKİNE MÜHRÜ: 188/dilim_44.npz (τ'∈(1.080,1.085]) aynı liste + seri_ve_G ile
     dds/Gre/Gim/q/a/τ BİT-BİT; tek-dilim K (50 bant) vs 188/harita_K_gercek.npz
     K[:,44] (maks|Δ|); doğrusallık: dilim-44 çizgilerinin q' mod 3 sınıf
     serileri toplamı = dilim serisi.
 A   beş uydu penceresi × sınıf (q' mod a; 'bölünen' ayrı; +log2 için ek mod-4
     bölmesi) → sınıf serileri (alt-örneklem her 3. nokta; kontrol noktalı).
 B   tek geçişte izdüşüm (alt-örneklem) → K_r(HAVUZ) tam + 8 loo.
 C   tam-örneklem kontrolü: +log3 r=1 ve −log3 r=1 (N = 299 999).
Çıktı: 192/seri_<uydu>_<sinif>.npz, 192/muhur_192.json, 192/sinif_proj.npz,
       192/tam_<..>.npz, 192/K2_sinif.json, 192/K2_sinif.npz
Koşu: nice -n 19, TEK süreç, OMP/BLAS = 1.
"""
import hashlib
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
           "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S184, S185, S186, S188, S192 = (SCR / d for d in ["184", "185", "186", "188", "192"])
NJACK = 8


def yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec)
    sys.modules[ad] = m
    spec.loader.exec_module(m)
    return m


b188 = yukle("b188", QM / "188_configs" / "188b_harita.py")
b185 = b188.b185
b187 = b188.b187


def onkayit():
    o = json.load(open(S192 / "ONKAYIT_192.json"))
    s = hashlib.sha256((QM / "192_configs" / "192a_onkayit.py").read_bytes()).hexdigest()
    if s != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {s} != {o['sha256']}")
    return o


def wrap(x):
    return (np.asarray(x, float) + 180.0) % 360.0 - 180.0


def jk(r):
    r = np.asarray(r, float)
    return float(np.sqrt((NJACK - 1) / NJACK * np.sum((r - r.mean()) ** 2)))


def aci_se(tam_deg, reps_deg):
    """188b se_aci deseni: sarılı sapmalar."""
    d = wrap(np.asarray(reps_deg) - tam_deg)
    return float(np.sqrt((NJACK - 1) / NJACK * np.sum((d - d.mean()) ** 2)))


def log(*a):
    print(*a, flush=True)


def seri_kontrol(yol, g_, mid_, bid_, w_l, a_l, etiket):
    if yol.exists():
        d = np.load(yol)
        log(f"    [{etiket}] kontrol noktası var ({len(a_l)} çizgi)")
        return d["dds"]
    t0 = time.time()
    dds, _, _ = b188.seri_ve_G(g_, mid_, bid_, w_l, a_l)
    tmp = str(yol) + ".tmp.npz"
    np.savez(tmp, dds=dds, w=w_l, a=a_l)
    os.replace(tmp, yol)
    log(f"    [{etiket}] seri bitti: {len(a_l)} çizgi, {time.time()-t0:.0f}s")
    return dds


if __name__ == "__main__":
    T0 = time.time()
    ONK = onkayit()
    ONK188 = json.load(open(S188 / "ONKAYIT_188.json"))
    L = float(ONK188["L"])
    A1 = ONK["A1"]["uydular"]
    log("=" * 78)
    log(f"192c / K2 A1 SINIF AYRIŞTIRMASI (son)  [on-kayit {ONK['zaman']} sha "
        f"{ONK['sha256'][:12]}]")
    log("=" * 78)

    # ---- kinematik + pencere (188b.kos AYNEN) ----
    G = np.load(S184 / "K1_gercek.npz")
    OZ = np.load(S185 / "OZ_gercek.npz")
    P = np.load(S186 / "G1_proj_gercek.npz")
    w_all = np.asarray(G["w"], float)
    tau_all = np.asarray(G["tau"], float)
    aq_all = np.asarray(G["aq"], float)
    win = (tau_all >= 0.45) & (tau_all < 0.86)
    w_win, tau_win = w_all[win], tau_all[win]
    g, mid, Lg = b185.kinematik("gercek")
    N = len(mid)
    kenar_jk = b188.jk_kenar(N)
    assert np.array_equal(np.diff(kenar_jk), G["nb"])
    bid_tam = np.searchsorted(kenar_jk, np.arange(N), side="right") - 1
    alt = np.arange(0, N, b188.ALT)
    g_a, mid_a, bid_a = g[alt], mid[alt], bid_tam[alt]
    maskeler, et = b188.bant_maskeleri(tau_win, ONK188)
    iH = len(et) - 1
    assert et[iH] == "HAVUZ"
    log(f"  N={N} alt={len(alt)} L={L:.9f} pencere çizgi={int(win.sum())}")

    t0 = time.time()
    qa, lam, taua, aqa = b187.asal_kuvvetler(np.exp(1.30 * L), L)
    log(f"  asal-kuvvetler: {len(qa)} ({time.time()-t0:.0f}s)")
    qi = qa.astype(np.int64)
    mix_reps = [b188.karisim(P, OZ, aq_all, j)[win] for j in range(NJACK)]
    mix_tam = b188.karisim(P, OZ, aq_all)[win]

    def K_hav(re, im, nb, Nn, disari):
        C = b188.c_proj(re, im, nb, Nn, disari)
        mix = mix_tam if disari < 0 else mix_reps[disari]
        return b188.K_matris(C, mix, maskeler)[0]

    muhur = {}
    # ================= M: makine mührü (dilim 44) =================
    kenar = np.array(ONK188["dilim_izgara_gercek"])
    i44 = 44
    idx = np.where((taua > kenar[i44]) & (taua <= kenar[i44 + 1]))[0]
    t0 = time.time()
    dd, Gre, Gim = b188.seri_ve_G(g_a, mid_a, bid_a, np.log(qa[idx]), aqa[idx])
    R = np.load(S188 / f"dilim_{i44}.npz")
    bit = {"dds": bool(np.array_equal(dd, R["dds"])),
           "Gre": bool(np.array_equal(Gre, R["Gre"])),
           "Gim": bool(np.array_equal(Gim, R["Gim"])),
           "q": bool(np.array_equal(qa[idx], R["q"])),
           "a": bool(np.array_equal(aqa[idx], R["a"])),
           "tau": bool(np.array_equal(taua[idx], R["tau"]))}
    bit["HEPSI"] = all(bit.values())
    log(f"  MÜHÜR dilim_44 ({len(idx)} çizgi, ({kenar[i44]:.3f},{kenar[i44+1]:.3f}]) "
        f"bit-bit: {bit}  ({time.time()-t0:.0f}s)")
    muhur["dilim44_bitbit"] = bit
    # tek-dilim K vs 188 harita_K_gercek
    re1, im1, nb_a, N_a = b188.izdusum(dd[:, None], mid_a, bid_a, w_win)
    K1t = K_hav(re1, im1, nb_a, N_a, -1)[:, 0]
    K1r = np.array([K_hav(re1, im1, nb_a, N_a, j)[:, 0] for j in range(NJACK)])
    HK = np.load(S188 / "harita_K_gercek.npz")
    dK = float(np.max(np.abs(K1t - HK["K"][:, i44])))
    dKr = float(np.max(np.abs(K1r - HK["K_reps"][:, :, i44])))
    muhur["dilim44_K_maks_fark_tam"] = dK
    muhur["dilim44_K_maks_fark_reps"] = dKr
    muhur["dilim44_K_HAVUZ"] = [float(K1t[iH].real), float(K1t[iH].imag)]
    log(f"  MÜHÜR tek-dilim K (50 bant) vs 188 harita_K: maks|Δ| tam {dK:.1e}, "
        f"loo {dKr:.1e}; K_HAVUZ = {K1t[iH]:.6f}")
    # doğrusallık: q' mod 3 sınıf serileri toplamı
    q44 = qi[idx]
    top = np.zeros_like(dd)
    for sel in ((q44 % 3 == 1), (q44 % 3 == 2), (np.gcd(q44, 3) > 1)):
        if sel.any():
            d_, _, _ = b188.seri_ve_G(g_a, mid_a, bid_a, np.log(qa[idx][sel]),
                                      aqa[idx][sel])
            top += d_
    dlin = float(np.max(np.abs(top - dd)))
    muhur["dilim44_dogrusallik_maks"] = dlin
    muhur["dilim44_dds_maks"] = float(np.max(np.abs(dd)))
    log(f"  DOĞRUSALLIK Σ_r seri_r − seri (dilim 44, mod 3) maks|Δ| = {dlin:.1e} "
        f"(maks|seri| = {np.max(np.abs(dd)):.3f})")
    json.dump(muhur, open(S192 / "muhur_192.json", "w"), indent=1, ensure_ascii=False)
    del dd, Gre, Gim, R, re1, im1, top, HK

    # ================= A: sınıf serileri =================
    siniflar = []   # (uydu, sınıf etiketi, seçim)
    for ad, w in A1.items():
        m = np.where((taua >= w["tau"][0]) & (taua <= w["tau"][1]))[0]
        q = qi[m]
        a = w["a"]
        if ad == "+log2":
            parca = [("r1mod4", q % 4 == 1), ("r3mod4", q % 4 == 3)]
        else:
            parca = [(f"r{r}", q % a == r) for r in w["siniflar"]]
        bol = np.gcd(q, a) > 1
        parca.append(("bolunen", bol))
        kap = np.zeros(len(q), bool)
        for _, s in parca:
            assert not (kap & s).any()
            kap |= s
        assert kap.all(), ad
        assert len(m) == w["cizgi_sayisi"], (ad, len(m), w["cizgi_sayisi"])
        for et_, s in parca:
            if s.any():
                siniflar.append((ad, et_, m[s]))
    log(f"  sınıf serileri: {len(siniflar)} (toplam çizgi "
        f"{sum(len(x[2]) for x in siniflar)})")
    D = np.zeros((len(alt), len(siniflar)))
    for k, (ad, et_, ix) in enumerate(siniflar):
        yol = S192 / f"seri_{ad.replace('+', 'p').replace('-', 'm')}_{et_}.npz"
        D[:, k] = seri_kontrol(yol, g_a, mid_a, bid_a, np.log(qa[ix]), aqa[ix],
                               f"{ad} {et_}")
        log(f"      ({k+1}/{len(siniflar)}; geçen {time.time()-T0:.0f}s)")

    # ================= B: izdüşüm =================
    yolB = S192 / "sinif_proj.npz"
    if not yolB.exists():
        t0 = time.time()
        re, im, nb_a, N_a = b188.izdusum(D, mid_a, bid_a, w_win)
        np.savez_compressed(yolB, re=re, im=im, nb=nb_a, N=N_a,
                            ad=np.array([x[0] for x in siniflar]),
                            sinif=np.array([x[1] for x in siniflar]))
        log(f"  [izdüşüm] -> {yolB.name} ({time.time()-t0:.0f}s)")
    PR = np.load(yolB)
    re, im, nb_a, N_a = PR["re"], PR["im"], PR["nb"], int(PR["N"])
    Kt = K_hav(re, im, nb_a, N_a, -1)[iH]                       # (nsınıf,)
    Kr = np.array([K_hav(re, im, nb_a, N_a, j)[iH] for j in range(NJACK)])

    # ================= C: tam-örneklem kontrolü =================
    tam = {}
    for ad, et_ in (("+log3", "r1"), ("-log3", "r1")):
        k = [i for i, x in enumerate(siniflar) if x[0] == ad and x[1] == et_][0]
        ix = siniflar[k][2]
        tag = f"{ad.replace('+', 'p').replace('-', 'm')}_{et_}"
        yolT = S192 / f"tam_{tag}.npz"
        if not yolT.exists():
            t0 = time.time()
            dT = seri_kontrol(S192 / f"tamseri_{tag}.npz", g, mid, bid_tam,
                              np.log(qa[ix]), aqa[ix], f"TAM {ad} {et_}")
            reT, imT, nbT, NT = b188.izdusum(dT[:, None], mid, bid_tam, w_win)
            np.savez_compressed(yolT, re=reT, im=imT, nb=nbT, N=NT)
            log(f"  [tam kontrol {ad} {et_}] ({time.time()-t0:.0f}s)")
        T = np.load(yolT)
        KT = K_hav(T["re"], T["im"], T["nb"], int(T["N"]), -1)[iH, 0]
        KTr = np.array([K_hav(T["re"], T["im"], T["nb"], int(T["N"]), j)[iH, 0]
                        for j in range(NJACK)])
        se_alt = float(np.hypot(jk(Kr[:, k].real), jk(Kr[:, k].imag)))
        se_tam = float(np.hypot(jk(KTr.real), jk(KTr.imag)))
        fark = abs(Kt[k] - KT)
        tam[f"{ad} {et_}"] = {
            "K_alt": [float(Kt[k].real), float(Kt[k].imag)],
            "K_tam": [float(KT.real), float(KT.imag)],
            "se_alt": se_alt, "se_tam": se_tam, "fark": float(fark),
            "fark_bolu_se_alt": float(fark / se_alt),
            "aci_alt": float(np.degrees(np.angle(Kt[k]))),
            "aci_tam": float(np.degrees(np.angle(KT)))}
        log(f"  TAM-ÖRNEKLEM {ad} {et_}: K_alt={Kt[k]:.5f} K_tam={KT:.5f} "
            f"|Δ|={fark:.2e} (= {fark/se_alt:.2f} se_alt; se_tam {se_tam:.1e})")

    # ================= sonuçlar =================
    out = {"sha_onkayit": ONK["sha256"], "muhur": muhur, "tam_kontrol": tam,
           "uydular": {}}
    for ad, w in A1.items():
        ks = [i for i, x in enumerate(siniflar) if x[0] == ad]
        top_t = Kt[ks].sum()
        top_r = Kr[:, ks].sum(1)
        th_top = float(np.degrees(np.angle(top_t)))
        th_top_r = np.degrees(np.angle(top_r))
        U = {"K_top": [float(top_t.real), float(top_t.imag)],
             "mod_top": float(abs(top_t)), "se_mod_top": jk(np.abs(top_r)),
             "kappa_top": float(-top_t.real), "se_kappa_top": jk(-top_r.real),
             "aci_top": th_top, "se_aci_top": aci_se(th_top, th_top_r),
             "siniflar": {}}
        for i in ks:
            et_ = siniflar[i][1]
            kt, kr = Kt[i], Kr[:, i]
            th = float(np.degrees(np.angle(kt)))
            thr = np.degrees(np.angle(kr))
            rel = float(wrap(th - th_top))
            relr = wrap(thr - th_top_r)
            oran = float(abs(kt) / abs(top_t))
            oranr = np.abs(kr) / np.abs(top_r)
            U["siniflar"][et_] = {
                "n_cizgi": int(len(siniflar[i][2])),
                "sum_a": float(aqa[siniflar[i][2]].sum()),
                "K": [float(kt.real), float(kt.imag)], "mod": float(abs(kt)),
                "se_mod": jk(np.abs(kr)), "aci": th, "se_aci": aci_se(th, thr),
                "aci_goreli_top": rel, "se_aci_goreli": aci_se(rel, relr),
                "oran_top": oran, "se_oran": jk(oranr)}
        # sınıf farkları
        def fark(r1, r2):
            i1 = [i for i in ks if siniflar[i][1] == r1][0]
            i2 = [i for i in ks if siniflar[i][1] == r2][0]
            d = float(wrap(np.degrees(np.angle(Kt[i2])) - np.degrees(np.angle(Kt[i1]))))
            dr = wrap(np.degrees(np.angle(Kr[:, i2])) - np.degrees(np.angle(Kr[:, i1])))
            return {"fark": d, "se": aci_se(d, dr)}
        ff = {}
        if ad in ("+log3", "-log3"):
            ff["1->2"] = fark("r1", "r2")
        elif ad == "+log5":
            for r in (1, 2, 3):
                ff[f"{r}->{r+1}"] = fark(f"r{r}", f"r{r+1}")
            ff["1->4"] = fark("r1", "r4")
        elif ad == "+log6":
            ff["1->5"] = fark("r1", "r5")
        elif ad == "+log2":
            ff["r1mod4->r3mod4"] = fark("r1mod4", "r3mod4")
        U["farklar"] = ff
        # +log2: mod-2 tek sınıfı = mod-4 sınıflarının toplamı (bölünen yok)
        if ad == "+log2":
            ko = [i for i in ks if siniflar[i][1] in ("r1mod4", "r3mod4")]
            t1 = Kt[ko].sum()
            t1r = Kr[:, ko].sum(1)
            th1 = float(np.degrees(np.angle(t1)))
            rel = float(wrap(th1 - th_top))
            U["tek_sinif_mod2"] = {
                "K": [float(t1.real), float(t1.imag)], "oran_top": float(abs(t1) / abs(top_t)),
                "se_oran": jk(np.abs(t1r) / np.abs(top_r)), "aci_goreli_top": rel,
                "se_aci_goreli": aci_se(rel, wrap(np.degrees(np.angle(t1r)) - th_top_r))}
        out["uydular"][ad] = U
        log(f"\n  [{ad}] τ'{w['tau']}: K_top = {top_t:.5f} |K| {abs(top_t):.5f}±"
            f"{U['se_mod_top']:.5f} ∠{th_top:+.1f}°±{U['se_aci_top']:.1f}  "
            f"κ = {-top_t.real:.5f}±{U['se_kappa_top']:.5f}")
        for et_, c in U["siniflar"].items():
            log(f"     sınıf {et_:>8} ({c['n_cizgi']:5d} çizgi): |K| {c['mod']:.5f}±"
                f"{c['se_mod']:.5f} ∠{c['aci']:+7.1f}°±{c['se_aci']:.1f} | göreli "
                f"{c['aci_goreli_top']:+7.1f}°±{c['se_aci_goreli']:.1f} | oran "
                f"{c['oran_top']:.3f}±{c['se_oran']:.3f}")
        for k_, v in ff.items():
            log(f"     Δ_{k_} = {v['fark']:+.1f}°±{v['se']:.1f}")
        if ad == "+log2":
            t = U["tek_sinif_mod2"]
            log(f"     mod-2 tek sınıfı: oran {t['oran_top']:.4f} göreli "
                f"{t['aci_goreli_top']:+.2f}°")
    np.savez_compressed(S192 / "K2_sinif.npz", Kt=Kt, Kr=Kr,
                        ad=np.array([x[0] for x in siniflar]),
                        sinif=np.array([x[1] for x in siniflar]))
    json.dump(out, open(S192 / "K2_sinif.json", "w"), indent=1, ensure_ascii=False)
    log(f"\n-> K2_sinif.json, K2_sinif.npz  BİTTİ ({time.time()-T0:.0f}s)")
