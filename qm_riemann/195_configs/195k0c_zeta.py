# -*- coding: utf-8 -*-
"""
195k0c — K0c ZETA KAPISI (genelleştirilmiş kod, χ ≡ 1, k = 1)
=============================================================
Ölçümden (adalar) ÖNCE. Biri tutmazsa DUR.

 M1  seri fonksiyonu (195o.seri) 188'in dilim_44.npz'ini (son, alt-örneklem,
     τ'∈(1.080,1.085]) BİT-BİT üretir (dds, q, a, τ).
 M2  195o.izdusum_blok ≡ 190b.izdusum_blok (bit-bit) ve 195o.c_oz_blok ≡
     185b.c_oz (düşük zinciri OZ_gercek_dusuk, tüm pencere; göreli maks fark).
 R   ESKİ YOL referansı (mix-normalize κ = −Re K_HAVUZ; 193b.blok_yerel_D +
     188b.c_proj/K_matris + karışım AYNEN, a=1 tek sınıf): +log2, +log3, +log6,
     +log10; son (193 blokları, ALT=3) ve düşük (190 blokları, ALT=1).
     +log3/+log10 (son) 193 K1_sinif κ_top(+bölünen) ile hane hane.
 Y   YENİ YOL (195o.cekirdek_blok; öz-terime normalize K̃, blok-yerel pencere
     çizgileri): 10 uydu; son (193 blokları, ALT=3), düşük (190 blokları,
     ALT=1), düşük ΔL=0.05 ızgarası (8.5 çapalı; 17 blok + ada-kısmi bloklar).
 G   ORAN KAPISI: ρ(n) = K(+log n)/K(+log2), n ∈ {3,6,10}: |ρ_yeni/ρ_eski − 1|
     ≤ 0.10 (son VE düşük; 6 oran).
 İ   İŞARET REFERANSI: sign Re K̃_ζ (son); düşük ile uyum KAYIT.
Çıktı: scratchpad/195/K0c.json, K0c_bloklar.npz
Kullanım: 195k0c_zeta.py [işçi]   (nohup)
"""
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

HERE = Path(__file__).resolve().parent


def yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec)
    sys.modules[ad] = m
    spec.loader.exec_module(m)
    return m


o = yukle("o195", HERE / "195o_ortak.py")
b193 = yukle("b193", o.QM / "193_configs" / "193b_sinif.py")
b188 = b193.b188
b190 = b193.b190
b185 = b188.b185
b187 = b188.b187
S184, S185, S186 = (o.SCR / d for d in ["184", "185", "186"])
NJ = 8
REF_UYDU = ["+log2", "+log3", "+log6", "+log10"]


def log(*a):
    print(*a, flush=True)


# ------------------------------------------------------------ ortak işçi durumu
_W = {}


def _init(d):
    _W.update(d)


def _yeni_is(arg):
    et, idx_full, idx_an, Lb = arg
    t0 = time.time()
    r = o.cekirdek_blok(_W["g"], _W["mid"], idx_full, idx_an, Lb, _W["Q"], _W["A"],
                        _W["CHI"], o.UYDU_SIRA)
    r["et"] = et
    r["sure"] = time.time() - t0
    # τ' < 0.86 (pencere-içi) uydu çizgisi payı
    lq = np.log(_W["Q"])
    ic = []
    for ad in o.UYDU_SIRA:
        sel = np.abs(lq - Lb - o.log_n(ad)) < o.DELTA
        ic.append(int(np.sum(sel & (lq / Lb < o.TAU_HI))))
    r["ncz_pencere_ici"] = np.array(ic)
    del r["c_uydu"]
    return r


def _eski_is(arg):
    pen, ad = arg
    E = _W["eski"][pen]
    b193.w_win_G[0] = E["w_win"]
    h = {"a": 1, "siniflar": [0], "log_n": o.log_n(ad)}
    t0 = time.time()
    RE, IM, ncz, sin_ad, tau_ar = b193.blok_yerel_D(
        ad, h, E["L_B"], E["qi"], E["qa"], E["taua"], E["aqa"], E["kenar_jk_a"],
        E["g_a"], E["mid_a"], E["bid_a"])
    return pen, ad, RE, IM, ncz, tau_ar, time.time() - t0


def eski_hazirla(pen):
    ONK188 = json.load(open(o.S188 / "ONKAYIT_188.json"))
    if pen == "son":
        G = np.load(S184 / "K1_gercek.npz")
        OZ = np.load(S185 / "OZ_gercek.npz")
        P = np.load(S186 / "G1_proj_gercek.npz")
        g, mid, _ = b185.kinematik("gercek")
        ALT = 3
        L_B = np.array(json.load(open(o.S193 / "ONKAYIT_193.json"))["bloklar"]["L_b"])
        Lref = float(ONK188["L"])
    else:
        G = np.load(o.ZD / "K1_gercek_dusuk.npz")
        OZ = np.load(o.ZD / "OZ_gercek_dusuk.npz")
        P = np.load(o.ZD / "G1_proj_gercek_dusuk.npz")
        g, mid, _ = b190.kinematik_eta(o.ZD / "eta_dusuk_t0.4_c4000.npz")
        ALT = 1
        O190 = json.load(open(o.S190 / "ONKAYIT_190.json"))
        L_B = np.array(O190["bloklar"]["L_b"])
        Lref = float(O190["pencere"]["L"])
    w_all = np.asarray(G["w"], float)
    tau_all = np.asarray(G["tau"], float)
    aq_all = np.asarray(G["aq"], float)
    win = (tau_all >= 0.45) & (tau_all < 0.86)
    w_win, tau_win = w_all[win], tau_all[win]
    N = len(mid)
    kenar_jk = b188.jk_kenar(N)
    assert np.array_equal(np.diff(kenar_jk), G["nb"])
    L_chk = np.array([float(np.log(mid[kenar_jk[b]:kenar_jk[b + 1]] / o.TWO_PI).mean())
                      for b in range(NJ)])
    assert np.allclose(L_chk, L_B, atol=1e-12), (pen, "L_B")
    bid = np.searchsorted(kenar_jk, np.arange(N), side="right") - 1
    alt = np.arange(0, N, ALT)
    g_a, mid_a, bid_a = g[alt], mid[alt], bid[alt]
    kenar_jk_a = np.searchsorted(bid_a, np.arange(NJ + 1))
    nb_a = np.diff(kenar_jk_a)
    maskeler, et = b188.bant_maskeleri(tau_win, ONK188)
    iH = len(et) - 1
    assert et[iH] == "HAVUZ"
    qa, lam, taua, aqa = b187.asal_kuvvetler(np.exp(1.30 * Lref), Lref)
    mix_tam = b188.karisim(P, OZ, aq_all)[win]
    mix_reps = [b188.karisim(P, OZ, aq_all, j)[win] for j in range(NJ)]
    return dict(w_win=w_win, L_B=L_B, qi=qa.astype(np.int64), qa=qa, taua=taua,
                aqa=aqa, kenar_jk_a=kenar_jk_a, g_a=g_a, mid_a=mid_a, bid_a=bid_a,
                nb_a=nb_a, N_a=len(alt), maskeler=maskeler, iH=iH, mix_tam=mix_tam,
                mix_reps=mix_reps, g=g, mid=mid, kenar_jk=kenar_jk, ALT=ALT, Lref=Lref)


def eski_kappa(E, RE, IM):
    def K_hav(disari):
        C = b188.c_proj(RE, IM, E["nb_a"], E["N_a"], disari)
        mix = E["mix_tam"] if disari < 0 else E["mix_reps"][disari]
        return b188.K_matris(C, mix, E["maskeler"])[0][E["iH"]]
    Kt = K_hav(-1).sum()          # r0 + bölünen (boş)
    Kr = np.array([K_hav(j).sum() for j in range(NJ)])
    return complex(Kt), Kr


if __name__ == "__main__":
    T0 = time.time()
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    log("=" * 78)
    log("195k0c — K0c ZETA KAPISI (genelleştirilmiş kod, χ≡1, k=1)")
    log("=" * 78)
    out = {"zaman_baslangic": time.ctime()}

    # ================= M1: seri bit-bit (188 dilim_44) =================
    ONK188 = json.load(open(o.S188 / "ONKAYIT_188.json"))
    Lson = float(ONK188["L"])
    Ks = o.zeta_kin("son")
    g_s, mid_s = Ks["g"], Ks["mid"]
    g185, mid185, _ = b185.kinematik("gercek")
    assert np.array_equal(g_s, g185) and np.array_equal(mid_s, mid185), "son kinematik"
    alt = np.arange(0, len(mid_s), 3)
    t0 = time.time()
    Qs, LAMs, As = o.asal_kuvvetler(np.exp(1.30 * Lson))
    taus = np.log(Qs) / Lson
    log(f"  asal-kuvvetler (son evreni q≤e^(1.30L)={np.exp(1.30*Lson):.0f}): {len(Qs)} "
        f"({time.time()-t0:.0f}s)")
    kenar = np.array(ONK188["dilim_izgara_gercek"])
    idx = np.where((taus > kenar[44]) & (taus <= kenar[45]))[0]
    dd = o.seri(g_s[alt], mid_s[alt], np.log(Qs[idx]), 1.0 * As[idx])
    R = np.load(o.S188 / "dilim_44.npz")
    M1 = {"dds": bool(np.array_equal(dd, R["dds"])),
          "q": bool(np.array_equal(Qs[idx], R["q"])),
          "a": bool(np.array_equal(As[idx], R["a"])),
          "tau": bool(np.array_equal(taus[idx], R["tau"])),
          "n_cizgi": int(len(idx)), "maks_fark_dds": float(np.max(np.abs(dd - R["dds"])))}
    M1["HEPSI"] = all(M1[k] for k in ["dds", "q", "a", "tau"])
    log(f"  M1 seri vs 188 dilim_44 ({len(idx)} çizgi, alt-örneklem) bit-bit: {M1}")
    out["M1_seri_bitbit"] = M1
    # M2a: izdüşüm bit-bit (190b.izdusum_blok)
    Dt = np.stack([dd[:5000], dd[5000:10000]], axis=1)
    wt = np.log(Qs[(taus >= 0.45) & (taus < 0.86)])[:400]
    r1, i1 = o.izdusum_blok(Dt, mid_s[alt][:5000], wt)
    r2, i2 = b190.izdusum_blok(Dt, mid_s[alt][:5000], wt)
    M2a = bool(np.array_equal(r1, r2) and np.array_equal(i1, i2))
    log(f"  M2a izdusum_blok ≡ 190b.izdusum_blok bit-bit: {M2a}")
    out["M2a_izdusum_bitbit"] = M2a
    del dd, R, Dt
    # M2b: c^öz ≡ 185b.c_oz (düşük zinciri, tüm pencere, tüm çizgi evreni)
    t0 = time.time()
    Gd = np.load(o.ZD / "K1_gercek_dusuk.npz")
    OZd = np.load(o.ZD / "OZ_gercek_dusuk.npz")
    Kd0 = o.zeta_kin("dusuk")
    coz_yeni = o.c_oz_blok(Kd0["g"], Kd0["mid"], np.asarray(Gd["w"], float),
                           np.asarray(Gd["aq"], float))
    coz_185 = b185.c_oz(OZd, np.asarray(Gd["aq"], float))
    M2b = float(np.max(np.abs(coz_yeni - coz_185)) / np.max(np.abs(coz_185)))
    log(f"  M2b c_oz_blok vs 185b.c_oz (düşük, {len(coz_185)} çizgi, tüm pencere): "
        f"maks|Δ|/maks|c| = {M2b:.1e}  ({time.time()-t0:.0f}s)")
    out["M2b_coz_goreli_maks_fark"] = M2b
    M2b_ok = M2b < 1e-10
    del coz_yeni, coz_185, Kd0

    # ================= R: eski yol hazırlığı =================
    t0 = time.time()
    ESKI = {pen: eski_hazirla(pen) for pen in ["son", "dusuk"]}
    log(f"  eski yol hazırlığı ({time.time()-t0:.0f}s)")

    # ================= Y: yeni yol blok listeleri =================
    Kd = o.zeta_kin("dusuk")
    ONK193 = json.load(open(o.S193 / "ONKAYIT_193.json"))
    ONK190 = json.load(open(o.S190 / "ONKAYIT_190.json"))
    isler_son, isler_dus = [], []
    bl = o.esit_sayim_bloklari(len(mid_s))
    for b, idx_f in enumerate(bl):
        Lb = float(np.log(mid_s[idx_f] / o.TWO_PI).mean())
        assert abs(Lb - ONK193["bloklar"]["L_b"][b]) < 1e-12
        idx_an = alt[(alt >= idx_f[0]) & (alt <= idx_f[-1])]
        isler_son.append((f"son_b{b}", idx_f, idx_an, Lb))
    bl = o.esit_sayim_bloklari(len(Kd["mid"]))
    for b, idx_f in enumerate(bl):
        Lb = float(np.log(Kd["mid"][idx_f] / o.TWO_PI).mean())
        assert abs(Lb - ONK190["bloklar"]["L_b"][b]) < 1e-12
        isler_dus.append((f"d190_b{b}", idx_f, idx_f, Lb))
    # ΔL ızgarası (8.5 çapalı) + ada-kısmi üst bloklar
    Lm_d = Kd["Lm"]
    grid = o.dL_bloklari(Lm_d, o.L_UST)
    isler_grid = []
    for (j, lo, hi, idx_f) in grid:
        isler_grid.append((f"dg_j{j}_{lo:.2f}_{hi:.4f}", idx_f, idx_f,
                           float(Lm_d[idx_f].mean())))
    Umax = {}
    for ad in o.ADALAR:
        K = o.ada_kin(ad)
        Umax[ad] = float(K["Lm"].max())
    kismi_et = {}
    for ad in o.ADALAR:
        U = min(o.BANT[1], Umax[ad])
        for (j, lo, hi, idx_f) in o.dL_bloklari(Lm_d, o.L_UST, ust_sinir=U):
            if lo < o.BANT[0] - 1e-9:
                continue
            et = f"dg_j{j}_{lo:.2f}_{hi:.4f}"
            if et not in [x[0] for x in isler_grid]:
                if et not in kismi_et:
                    isler_grid.append((et, idx_f, idx_f, float(Lm_d[idx_f].mean())))
                kismi_et.setdefault(et, []).append(ad)
    log(f"  yeni yol blokları: son {len(isler_son)}, düşük-190 {len(isler_dus)}, "
        f"düşük-ızgara {len(isler_grid)} (kısmi ek: {kismi_et})")

    # ================= işçiler =================
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    Qd, _, Ad = o.asal_kuvvetler(np.exp(1.30 * Kd["L"]))
    sonuc_bloklar = {}
    for pen, isler, Q, A, K in [("son", isler_son, Qs, As, Ks),
                                ("dusuk", isler_dus + isler_grid, Qd, Ad, Kd)]:
        t0 = time.time()
        _W.clear()
        _W.update(g=K["g"], mid=K["mid"], Q=Q, A=A, CHI=np.ones(len(Q)))
        with ctx.Pool(nw) as pool:
            for r in pool.imap_unordered(_yeni_is, sorted(isler, key=lambda x: -len(x[2]))):
                sonuc_bloklar[r["et"]] = r
                log(f"    [yeni {pen}] {r['et']} N={r['N_tam']} (an {r['N_an']}) "
                    f"pencere={r['ncz_pencere']} uydu çizgi={r['ncz_uydu'].tolist()} "
                    f"({r['sure']:.0f}s)")
        log(f"  [yeni {pen}] bitti ({time.time()-t0:.0f}s)")

    # eski yol (8 iş: 2 pencere × 4 uydu)
    t0 = time.time()
    _W.clear()
    _W["eski"] = ESKI
    eski = {}
    with ctx.Pool(min(nw, 8)) as pool:
        for pen, ad, RE, IM, ncz, tau_ar, dt in pool.imap_unordered(
                _eski_is, [(p, a) for p in ["son", "dusuk"] for a in REF_UYDU]):
            Kt, Kr = eski_kappa(ESKI[pen], RE, IM)
            kap = -Kt.real
            se = float(o.jk_se(-Kr.real))
            eski[(pen, ad)] = dict(K=[Kt.real, Kt.imag], kappa=kap, se_kappa=se,
                                   kappa_reps=(-Kr.real).tolist(),
                                   n_cizgi=int(ncz.sum()), tau_aralik=list(tau_ar))
            log(f"    [eski {pen}] {ad}: κ={kap:+.6f}±{se:.6f} ∠K={np.degrees(np.angle(Kt)):+.1f}° "
                f"çizgi={int(ncz.sum())} ({dt:.0f}s)")
    log(f"  [eski] bitti ({time.time()-t0:.0f}s)")

    # 193 hane-hane kontrolü
    K193 = json.load(open(o.S193 / "K1_sinif.json"))["hedefler"]
    h193 = {}
    for ad in ["+log3", "+log10"]:
        ref = K193[ad]["kappa_top"] + K193[ad]["bolunen"]["kappa"]
        v = eski[("son", ad)]["kappa"]
        h193[ad] = {"eski_yol": v, "193_kappa_top_arti_bolunen": ref,
                    "goreli_fark": abs(v - ref) / abs(ref)}
        log(f"  193 yeniden üretimi {ad}: eski yol κ={v:.8f}  193={ref:.8f}  "
            f"göreli fark {abs(v-ref)/abs(ref):.1e}")
    out["R_193_yeniden_uretim"] = h193

    # ================= havuzlar =================
    def havuzla(etiketler):
        num = np.array([sonuc_bloklar[e]["num"] for e in etiketler])
        den = np.array([sonuc_bloklar[e]["den"] for e in etiketler])
        Nw = np.array([sonuc_bloklar[e]["N_tam"] for e in etiketler])
        K, reps = o.havuz(num, den, Nw)
        s = o.ozet(K, reps)
        return {ad: {"re": float(s["re"][i]), "im": float(s["im"][i]),
                     "se_re": float(s["se_re"][i]), "se_im": float(s["se_im"][i]),
                     "z_re": float(s["z_re"][i]), "z_im": float(s["z_im"][i]),
                     "aci": float(s["aci"][i]),
                     "n_cizgi": int(sum(sonuc_bloklar[e]["ncz_uydu"][i] for e in etiketler)),
                     "n_pencere_ici_cizgi": int(sum(sonuc_bloklar[e]["ncz_pencere_ici"][i]
                                                    for e in etiketler))}
                for i, ad in enumerate(o.UYDU_SIRA)}, K, reps

    Y = {}
    Y["son"], Kson, Rson = havuzla([x[0] for x in isler_son])
    Y["dusuk190"], Kd190, Rd190 = havuzla([x[0] for x in isler_dus])
    std_grid = [f"dg_j{j}_{lo:.2f}_{hi:.4f}" for (j, lo, hi, _) in grid]
    Y["dusuk_izgara_tum"], _, _ = havuzla(std_grid)
    # bant referansları (ada-eşli aralık)
    bant_ref = {}
    for ad in o.ADALAR:
        U = min(o.BANT[1], Umax[ad])
        ets = []
        for (j, lo, hi, idx_f) in o.dL_bloklari(Lm_d, o.L_UST, ust_sinir=U):
            if lo < o.BANT[0] - 1e-9:
                continue
            ets.append(f"dg_j{j}_{lo:.2f}_{hi:.4f}")
        bant_ref[ad] = {"U": U, "bloklar": ets}
        bant_ref[ad]["K"], _, _ = havuzla(ets)
    Y["bant_ref"] = bant_ref
    out["Y_yeni_Ktilde"] = Y

    log("\n  YENİ YOL K̃_ζ (öz-normalize; Re±se, z_Re; Im, z_Im; açı):")
    for pen in ["son", "dusuk190", "dusuk_izgara_tum"]:
        log(f"   [{pen}]")
        for ad in o.UYDU_SIRA:
            v = Y[pen][ad]
            log(f"     {ad:>10}: Re={v['re']:+.5f}±{v['se_re']:.5f} (z={v['z_re']:+.1f})  "
                f"Im={v['im']:+.5f} (z={v['z_im']:+.1f})  ∠{v['aci']:+.1f}°  "
                f"çizgi={v['n_cizgi']} pencere-içi={v['n_pencere_ici_cizgi']}")

    # ================= G: oran kapısı =================
    G = {}
    gecti = True
    for pen, yk in [("son", "son"), ("dusuk", "dusuk190")]:
        G[pen] = {}
        for ad in ["+log3", "+log6", "+log10"]:
            r_eski = eski[(pen, ad)]["kappa"] / eski[(pen, "+log2")]["kappa"]
            r_yeni = Y[yk][ad]["re"] / Y[yk]["+log2"]["re"]
            sap = r_yeni / r_eski - 1
            ok = abs(sap) <= 0.10
            gecti &= ok
            G[pen][ad] = {"rho_eski": r_eski, "rho_yeni": r_yeni, "sapma": sap, "ok": ok}
            log(f"  ORAN KAPISI {pen:>5} {ad:>7}/+log2: eski {r_eski:.4f}  yeni {r_yeni:.4f}"
                f"  sapma {sap:+.3f} → {'TUTTU' if ok else 'TUTMADI'}")
    out["G_oran_kapisi"] = G
    out["G_GECTI"] = bool(gecti)
    out["R_eski_kappa"] = {f"{p}|{a}": v for (p, a), v in eski.items()}

    # ================= İ: işaret referansı =================
    ISR = {}
    for ad in o.UYDU_SIRA:
        s_son = Y["son"][ad]
        s_dus = Y["dusuk190"][ad]
        ISR[ad] = {"isaret_son": int(np.sign(s_son["re"])), "z_son": s_son["z_re"],
                   "isaret_dusuk190": int(np.sign(s_dus["re"])), "z_dusuk190": s_dus["z_re"],
                   "isaret_dusuk_izgara": int(np.sign(Y["dusuk_izgara_tum"][ad]["re"])),
                   "z_dusuk_izgara": Y["dusuk_izgara_tum"][ad]["z_re"],
                   "uyum": bool(np.sign(s_son["re"]) == np.sign(s_dus["re"]))}
    out["ISARET_REFERANSI"] = ISR
    out["M_hepsi"] = bool(M1["HEPSI"] and M2a and M2b_ok)
    out["K0c_GECTI"] = bool(M1["HEPSI"] and M2a and M2b_ok and gecti)
    out["Umax_ada"] = Umax
    out["zaman_bitis"] = time.ctime()
    np.savez_compressed(
        o.S195 / "K0c_bloklar.npz",
        et=np.array(list(sonuc_bloklar.keys())),
        num=np.array([sonuc_bloklar[e]["num"] for e in sonuc_bloklar]),
        den=np.array([sonuc_bloklar[e]["den"] for e in sonuc_bloklar]),
        N=np.array([sonuc_bloklar[e]["N_tam"] for e in sonuc_bloklar]),
        N_an=np.array([sonuc_bloklar[e]["N_an"] for e in sonuc_bloklar]),
        Lb=np.array([sonuc_bloklar[e]["Lb"] for e in sonuc_bloklar]),
        ncz=np.array([sonuc_bloklar[e]["ncz_uydu"] for e in sonuc_bloklar]),
        uydular=np.array(o.UYDU_SIRA))
    json.dump(out, open(o.S195 / "K0c.json", "w"), indent=1, ensure_ascii=False,
              default=lambda x: x.tolist() if hasattr(x, "tolist") else str(x))
    log(f"\nK0c: M1 {M1['HEPSI']}  M2a {M2a}  ORAN KAPISI {'GEÇTİ' if gecti else 'KALDI — DUR'}")
    log(f"-> K0c.json, K0c_bloklar.npz  ({time.time()-T0:.0f}s)")
