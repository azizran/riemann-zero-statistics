# -*- coding: utf-8 -*-
"""
190c — K2: HÜKÜM (H-190a/b/c; ONKAYIT_190 eşikleri AYNEN, kurtarma yok) + KAYIT
==============================================================================
 (a) H-190a: düşük pencere, blok başına HAVUZ κ_b(Δω) (ω-dilim 0.025, 190b D),
     tepe_b = hedef ±0.15 içindeki dilim merkezlerinde κ maksimumu; blok-medyan
     (+ blok-loo jk se); hedef/rakip karşılaştırması. KAYIT: −log2/−log3,
     tepe_b'nin L_b'ye eğimi (0 evrensel / −1 τ'-sabit), τ'-dilim profili
     çapraz kontrolü, SON penceresi profili (188 dosyalarından; yeniden koşu yok).
 (b) H-190b: 188d.analiz AYNEN (harita_K_dusuk) → corr(Re v, S_Re) ± jk.
 (c) H-190c: |corr(Re v, W)| ≥ |corr(Re v, S_Re)| ?
 (d) KAYIT: ζ_g düşük (190k0) vs son (187 K2), Σ_{τ'≤1.20/1.30} K, rang-1 payı.
Çıktı: 190/HUKUM_190.json, 190/profiller_190.npz (figür için) + ekran.
"""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S = {d: SCR / d for d in ["155", "184", "185", "186", "187", "188", "190"]}
ZD = S["190"] / "zincir_dusuk"
TWO_PI = 2 * np.pi
NJACK = 8


def yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec)
    sys.modules[ad] = m
    spec.loader.exec_module(m)
    return m


b185 = yukle("b185", QM / "185_configs" / "185b_oz_muhasebe.py")
b188d = yukle("b188d", QM / "188_configs" / "188d_analiz.py")


def onkayit():
    o = json.load(open(S["190"] / "ONKAYIT_190.json"))
    sha = hashlib.sha256(
        (QM / "190_configs" / "190a_onkayit.py").read_bytes()).hexdigest()
    if sha != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {sha} != {o['sha256']}")
    return o


def jk(r):
    r = np.asarray(r, float)
    return float(np.sqrt((NJACK - 1) / NJACK * np.sum((r - r.mean()) ** 2)))


def mix_havuz(K1, OZ, P):
    """188f(c) AYNEN: mix = c^kesik − c^öz (tam örneklem), HAVUZ pencere."""
    ck = 2 * (P["re_bir"].sum(0) + 1j * P["im_bir"].sum(0)) / int(P["N"])
    coz = b185.c_oz(OZ, K1["aq"])
    tau = K1["tau"]
    win = (tau >= 0.45) & (tau < 0.86)
    return (ck - coz)[win]


def L_bloklar(mid):
    kj = np.linspace(0, len(mid), 9).astype(int)
    return np.array([float(np.log(mid[kj[b]:kj[b + 1]] / TWO_PI).mean())
                     for b in range(8)]), kj


def tau_blok_profili(PR, mix, L, Lb):
    """188f(c) AYNEN: yalnız-blok K_HAVUZ (−Re) τ'-dilimlerinde; Δω = τ'·L − L_b."""
    k = PR["kenar"]
    o = 0.5 * (k[:-1] + k[1:])
    kap = np.zeros((8, len(o)))
    dw = np.zeros((8, len(o)))
    for bb in range(8):
        C = 2 * (PR["re"][bb] + 1j * PR["im"][bb]) / PR["nb"][bb]
        kap[bb] = -(C @ np.conj(mix)).real / np.sum(np.abs(mix) ** 2)
        dw[bb] = o * L - Lb[bb]
    return kap, dw


def tepe(x, y, hedef, pen=0.15, gecerli=None):
    m = np.abs(x - hedef) <= pen + 1e-9
    if gecerli is not None:
        m &= gecerli
    if not m.any():
        return float("nan"), float("nan")
    i = np.where(m)[0][np.argmax(y[m])]
    return float(x[i]), float(y[i])


def medyan_jk(t):
    t = np.asarray(t, float)
    loo = np.array([np.median(np.delete(t, i)) for i in range(len(t))])
    return float(np.median(t)), jk(loo)


def egim(Lb, t):
    """tepe_b = a + s·L_b (en küçük kareler) + blok-loo jk se."""
    def s_(x, y):
        return float(np.polyfit(x, y, 1)[0])
    s = s_(Lb, t)
    loo = [s_(np.delete(Lb, i), np.delete(t, i)) for i in range(len(t))]
    return s, jk(loo)


if __name__ == "__main__":
    ONK = onkayit()
    ONK188 = json.load(open(S["188"] / "ONKAYIT_188.json"))
    print("=" * 78)
    print(f"190c / K2 HÜKÜM  [on-kayit {ONK['zaman']} sha {ONK['sha256'][:12]}]")
    print("=" * 78)
    L = float(ONK["pencere"]["L"])
    L_SON = float(ONK["L_son"])
    HED = ONK["hedef_hukum"]
    HEDK = ONK["hedef_kayit"]
    RAK = ONK["rakip_tau_sabit"]
    PEN = 0.15
    TOL = 0.05
    out = {"sha_onkayit": ONK["sha256"]}

    # ================= (a) düşük ω-dilim profili =================
    H = np.load(S["190"] / "harita_omega_dusuk.npz")
    et = list(H["etiket"])
    iH = et.index("HAVUZ")
    merkez = H["merkez"]
    ncz = H["ncz"]
    Lb = H["L_b"]
    assert np.array_equal(Lb, np.array(ONK["bloklar"]["L_b"]))
    kap_om = -H["K"][:, iH, :].real            # (8, nJ)
    nb = H["nb"]
    N = int(nb.sum())
    print("\n(a) H-190a — blok başına tepe konumu (ω-dilim 0.025; HAVUZ κ; "
          "hedef ±0.15):")
    tablo = {}
    tum_hedef = dict(HED)
    tum_hedef.update(HEDK)
    for ad, h in tum_hedef.items():
        tp, kp = [], []
        for b in range(8):
            x, y = tepe(merkez, kap_om[b], h, PEN, ncz[b] > 0)
            tp.append(x)
            kp.append(y)
        tp = np.array(tp)
        med, se = medyan_jk(tp)
        s, se_s = egim(Lb, tp)
        r = RAK.get(ad)
        d_h = abs(med - h)
        d_r = abs(med - r) if r is not None else None
        kosul = (d_h <= TOL + 1e-12) and (r is None or d_h < d_r)
        tablo[ad] = {"hedef": h, "rakip": r, "tepe_blok": tp.tolist(),
                     "kappa_tepe": kp, "medyan": med, "se_medyan_jk": se,
                     "fark_hedef": med - h,
                     "fark_rakip": (med - r) if r is not None else None,
                     "kosul": bool(kosul), "egim_Lb": s, "se_egim": se_s,
                     "hukum_hedefi": ad in HED}
        print(f"  {ad:>6} hedef {h:+.4f} rakip "
              f"{('%+.4f' % r) if r is not None else '   —   '} | tepe_b = "
              f"{' '.join('%+.3f' % v for v in tp)} | medyan {med:+.4f}±{se:.4f} "
              f"(Δh {med-h:+.4f}{'' if r is None else ', Δr %+.4f' % (med-r)}) "
              f"{'GEÇTİ' if kosul else 'GEÇMEDİ'}  eğim(L_b) {s:+.3f}±{se_s:.3f}"
              f"{'' if ad in HED else '  [KAYIT]'}")
    out["H190a_tablo"] = tablo

    # τ'-dilim profili çapraz kontrolü (düşük) ve SON profili (188 dosyaları)
    K1d = np.load(ZD / "K1_gercek_dusuk.npz")
    mix_d = mix_havuz(K1d, np.load(ZD / "OZ_gercek_dusuk.npz"),
                      np.load(ZD / "G1_proj_gercek_dusuk.npz"))
    PRd = np.load(S["190"] / "harita_proj_dusuk.npz")
    kap_td, dw_td = tau_blok_profili(PRd, mix_d, L, Lb)
    mid_s = np.load(S["155"] / "eta_son_t0.4_c4000.npz")["mid"]
    Lb_s, _ = L_bloklar(mid_s)
    K1s = np.load(S["184"] / "K1_gercek.npz")
    mix_s = mix_havuz(K1s, np.load(S["185"] / "OZ_gercek.npz"),
                      np.load(S["186"] / "G1_proj_gercek.npz"))
    PRs = np.load(S["188"] / "harita_proj_gercek.npz")
    kap_ts, dw_ts = tau_blok_profili(PRs, mix_s, L_SON, Lb_s)
    print("\n  KAYIT — τ'-dilim (0.005) blok profilleri, aynı tepe kuralı "
          "(düşük: 0.052 ω; son: 0.060 ω; son 188 dosyalarından):")
    capraz = {}
    for ad, h in tum_hedef.items():
        td = np.array([tepe(dw_td[b], kap_td[b], h)[0] for b in range(8)])
        ts = np.array([tepe(dw_ts[b], kap_ts[b], h)[0] for b in range(8)])
        md, sd = medyan_jk(td)
        ms, ss = medyan_jk(ts)
        eg_s, eg_se = egim(Lb_s, ts)
        capraz[ad] = {"dusuk_tau_tepe_blok": td.tolist(), "dusuk_tau_medyan": md,
                      "dusuk_tau_se": sd, "son_tepe_blok": ts.tolist(),
                      "son_medyan": ms, "son_se": ss, "son_egim_Lb": eg_s,
                      "son_se_egim": eg_se}
        print(f"  {ad:>6} hedef {h:+.4f} | düşük(τ'-dilim) medyan {md:+.4f}±{sd:.4f}"
              f" | SON medyan {ms:+.4f}±{ss:.4f} (eğim {eg_s:+.3f}±{eg_se:.3f})")
    out["KAYIT_tau_dilim_ve_son"] = capraz
    out["L_b_son"] = Lb_s.tolist()

    # H-190a hükmü
    g2 = tablo["+log2"]["kosul"]
    hepsi = all(tablo[a]["kosul"] for a in HED)
    if hepsi:
        hA = "MÜHÜR"
    elif not g2:
        hA = "ÖLDÜ"
    else:
        hA = "KAYIT (kısmi: " + ", ".join(a for a in HED if not tablo[a]["kosul"]) \
            + " geçmedi)"

    # ================= (b)/(c) 188d.analiz AYNEN =================
    def analiz_kos(yol, Lx):
        Hk = np.load(yol)
        onk = json.loads(json.dumps(ONK188))
        onk["uydular"]["tau_p"] = {p: float(np.log(int(p)) / Lx)
                                  for p in ("2", "3", "5", "7")}
        tam = b188d.analiz(Hk["K"], Hk["pay"], Hk["S"], Hk["SRe"], Hk["W"],
                           Hk["orta"], Hk["kenar"], Hk["ince_kenar"], onk)
        reps = [b188d.analiz(Hk["K_reps"][j], Hk["pay_reps"][j], Hk["S_reps"][j],
                             Hk["SRe_reps"][j], Hk["W"], Hk["orta"], Hk["kenar"],
                             Hk["ince_kenar"], onk) for j in range(NJACK)]
        W = Hk["W"]
        c = {"pay1": tam["pay1"], "se_pay1": jk([r["pay1"] for r in reps]),
             "corr_RevSRe": tam["corr_RevSRe"],
             "se_RevSRe": jk([r["corr_RevSRe"] for r in reps]),
             "corr_RevW": b188d.pearson(tam["v"].real, W),
             "se_RevW": jk([b188d.pearson(r["v"].real, W) for r in reps]),
             "corr_vS": tam["corr_vS"], "se_vS": jk([r["corr_vS"] for r in reps]),
             "corr_vW": tam["corr_vW"], "se_vW": jk([r["corr_vW"] for r in reps])}
        fark = [abs(b188d.pearson(r["v"].real, W)) - abs(r["corr_RevSRe"])
                for r in reps]
        c["absfark_W_eksi_SRe"] = abs(c["corr_RevW"]) - abs(c["corr_RevSRe"])
        c["se_absfark"] = jk(fark)
        return c, tam, Hk

    cd, tam_d, Hd = analiz_kos(S["190"] / "harita_K_dusuk.npz", L)
    cs, tam_s, Hs = analiz_kos(S["188"] / "harita_K_gercek.npz", L_SON)
    print("\n(b)/(c) 188d.analiz AYNEN (88 τ'-dilim; ± 8-blok loo jk):")
    for ad, c in (("düşük", cd), ("son (188 yeniden)", cs)):
        print(f"  [{ad}] rang-1 payı {c['pay1']:.4f}±{c['se_pay1']:.4f} | "
              f"corr(Re v,S_Re) {c['corr_RevSRe']:+.4f}±{c['se_RevSRe']:.4f} | "
              f"corr(Re v,W) {c['corr_RevW']:+.4f}±{c['se_RevW']:.4f} | "
              f"corr(|v|,S) {c['corr_vS']:+.4f}±{c['se_vS']:.4f} | "
              f"corr(|v|,W) {c['corr_vW']:+.4f}±{c['se_vW']:.4f}")
    rb = cd["corr_RevSRe"]
    Hc = abs(cd["corr_RevW"]) >= abs(rb)
    if Hc:
        hB = "ÖLDÜ (H-190c tuttu: tarak resmi ölür)"
        hA = "ÖLDÜ (H-190c tuttu: tarak resmi ölür)"
        hC = "MÜHÜR (yapısız sıfır)"
    else:
        hC = "ÖLDÜ"
        if rb <= -0.80:
            hB = "MÜHÜR"
        elif abs(rb) < 0.60:
            hB = "ÖLDÜ"
        else:
            hB = "KAYIT" + (" (işaret ters)" if rb > 0 else "")
    out["H190b_c"] = {"dusuk": cd, "son_188_yeniden": cs}

    # ================= (d) KAYIT: ζ_g =================
    zd = json.load(open(ZD / "zincir_gercek_dusuk.json"))
    zs = json.load(open(S["187"] / "K2_zeta.json"))["gercek"]
    K1h = json.load(open(S["190"] / "K1_harita_dusuk.json"))
    K1h_s = json.load(open(S["188"] / "K1_harita_gercek.json"))
    print("\n(d) KAYIT — ζ_g (τ_c = 0.86; 187c AYNEN; ± jk):")
    print(f"{'bant':>10} {'düşük |ζ|':>16} {'açı°':>15} | {'son |ζ|':>16} {'açı°':>15}")
    for k, e in enumerate(zd["etiket187"]):
        print(f"{e:>10} {zd['zeta_mod'][k]:.4f}±{zd['se_mod'][k]:.4f} "
              f"{zd['zeta_aci'][k]:+8.2f}±{zd['se_aci'][k]:.2f} | "
              f"{zs['zeta_mod'][k]:.4f}±{zs['se_mod'][k]:.4f} "
              f"{zs['zeta_aci'][k]:+8.2f}±{zs['se_aci'][k]:.2f}")
    zh = zd["zeta_mod"][8]
    for ad in ("le120", "le130"):
        a = K1h[f"HAVUZ_toplam_{ad}"]
        b = K1h_s[f"HAVUZ_toplam_{ad}"]
        print(f"  HAVUZ Σ K ({ad}): düşük {a['mod']:.4f}±{a['se_mod']:.4f} "
              f"∠{a['aci']:+.2f}° (ζ_g'nin %{100*a['mod']/zh:.1f}'i) | son "
              f"{b['mod']:.4f}±{b['se_mod']:.4f} ∠{b['aci']:+.2f}° "
              f"(%{100*b['mod']/zs['zeta_mod'][8]:.1f})")
    out["KAYIT_zeta"] = {"dusuk": {k: zd[k] for k in ("etiket187", "zeta_mod", "se_mod",
                                                      "zeta_aci", "se_aci",
                                                      "genlik_okuma", "se_genlik")},
                         "son": {k: zs[k] for k in ("zeta_mod", "se_mod", "zeta_aci",
                                                    "se_aci")},
                         "sigmaK_dusuk": {a: K1h[f"HAVUZ_toplam_{a}"]
                                          for a in ("le120", "le130")},
                         "sigmaK_son": {a: K1h_s[f"HAVUZ_toplam_{a}"]
                                        for a in ("le120", "le130")}}
    out["makine"] = {"M1": K1h["M1"], "M2": K1h["M2"], "M3_maks": max(K1h["M3"])}

    # ================= HÜKÜM =================
    out["hukum"] = {"H-190a": hA, "H-190b": hB, "H-190c": hC,
                    "KAYIT": "ζ_g düşük + son profili + eğimler"}
    print("\nHÜKÜMLER (ön-kayıt eşikleri; kurtarma yok):")
    for a in HED:
        t = tablo[a]
        print(f"  H-190a [{a}] medyan {t['medyan']:+.4f} hedef {t['hedef']:+.4f} "
              f"(|Δ|={abs(t['fark_hedef']):.4f} ≤0.05?"
              f"{'' if t['rakip'] is None else ' ; rakip %.4f |Δr|=%.4f' % (t['rakip'], abs(t['fark_rakip']))})"
              f" → {'GEÇTİ' if t['kosul'] else 'GEÇMEDİ'}")
    print(f"  H-190a → {hA}")
    print(f"  H-190b: corr(Re v, S_Re) = {rb:+.4f}±{cd['se_RevSRe']:.4f} "
          f"(≤ −0.80 MÜHÜR; |corr|<0.60 ÖLÜM) → {hB}")
    print(f"  H-190c: |corr(Re v,W)| = {abs(cd['corr_RevW']):.4f} ≥ |corr(Re v,S_Re)| "
          f"= {abs(rb):.4f}? {'EVET' if Hc else 'HAYIR'} → {hC}")

    # ================= figür profilleri =================
    np.savez_compressed(
        S["190"] / "profiller_190.npz",
        om_merkez=merkez, om_kap=kap_om, om_ncz=ncz, om_nb=nb, Lb=Lb,
        td_kap=kap_td, td_dw=dw_td, ts_kap=kap_ts, ts_dw=dw_ts, Lb_s=Lb_s,
        nb_s=PRs["nb"], nb_d=PRd["nb"],
        v_d=tam_d["v"], SRe_d=Hd["SRe"], W_d=Hd["W"], orta=Hd["orta"],
        v_s=tam_s["v"], SRe_s=Hs["SRe"], W_s=Hs["W"])
    json.dump(out, open(S["190"] / "HUKUM_190.json", "w"), indent=1,
              ensure_ascii=False)
    print(f"\n-> {S['190']/'HUKUM_190.json'}, profiller_190.npz  BİTTİ")
