# -*- coding: utf-8 -*-
"""
188d — K3: ANALİZ + HÜKÜMLER (H-188a..d, K3 yeniden kurulum)
============================================================
ONKAYIT_188 AYNEN:
 (i)   41 ince bant × 88 dilim karmaşık K, ağırlıksız SVD → rang-1 payı,
       u(τ_b), v(τ') (faz: Σv reel-pozitif)
 (ii)  corr(|v|,S), corr(|v|,W)  [İKİNCİL: gürültü-düzeltmeli pay; corr(Re v,S_Re)]
 (iii) (1.05,1.10] κ-katkısının 1+τ_2 / 1+τ_3 uydu pencerelerindeki payı
 (iv)  rang-1 artığında köşegen şerit kontrastı (τ_b ≥ 0.70) + yükseliş payı
 (v)   H-188d: κ-profili gerçek vs ikiz (34 dilim × 0.01)
 (vi)  8-bant ζ^{≤1.30} yeniden kurulumu (u vs artık), bant-bant
Tüm se: 8-blok loo-jackknife (her replika baştan: SVD, korelasyon, pay).

Çıktı: 188/K3_analiz.json + ekran defteri.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S187, S188 = SCR / "187", SCR / "188"
NJACK = 8
NI = 41   # ince bant sayısı


def onkayit():
    o = json.load(open(S188 / "ONKAYIT_188.json"))
    sha = hashlib.sha256(
        (QM / "188_configs" / "188a_onkayit.py").read_bytes()).hexdigest()
    if sha != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {sha} != {o['sha256']}")
    return o


def jk_se(reps):
    reps = np.asarray(reps, float)
    return np.sqrt((NJACK - 1) / NJACK * np.sum((reps - reps.mean(0)) ** 2, 0))


def pearson(x, y):
    x = np.asarray(x, float) - np.mean(x)
    y = np.asarray(y, float) - np.mean(y)
    return float(np.sum(x * y) / np.sqrt(np.sum(x * x) * np.sum(y * y)))


def svd1(M):
    U, s, Vh = np.linalg.svd(M, full_matrices=False)
    v = Vh[0].copy()
    u = s[0] * U[:, 0]
    ph = np.angle(v.sum())
    v *= np.exp(-1j * ph)
    u *= np.exp(1j * ph)
    K1 = np.outer(u, v)
    return dict(pay=float(s[0] ** 2 / np.sum(s ** 2)), s=s, u=u, v=v, K1=K1)


def analiz(K, pay, S, SRe, W, orta, kenar, ince_kenar, ONK, tauH=None,
           KH=None):
    """Tek veri kümesi (tam ya da loo replikası) için tüm nicelikler."""
    out = {}
    M = K[:NI]
    r = svd1(M)
    out["pay1"] = r["pay"]
    v, u, K1 = r["v"], r["u"], r["K1"]
    out["absv"] = np.abs(v)
    out["v"] = v
    out["u"] = u
    out["corr_vS"] = pearson(np.abs(v), S)
    out["corr_vW"] = pearson(np.abs(v), W)
    out["corr_RevSRe"] = pearson(v.real, SRe)
    # (iii) (1.05,1.10] uydu payı
    kap = -K[-1].real                       # κ(s) = −Re K_HAVUZ
    out["kappa"] = kap
    tp = ONK["uydular"]["tau_p"]
    pen = ONK["uydular"]["pencere"]
    k510 = (kenar[:-1] >= 1.05 - 1e-9) & (kenar[1:] <= 1.10 + 1e-9)
    uy = np.zeros(len(orta), bool)
    for p in ("2", "3"):
        t = 1.0 + tp[p]
        uy |= np.abs(orta - t) <= pen
    out["uydu_pay_510"] = float(kap[k510 & uy].sum() / kap[k510].sum())
    out["uydu_sifir_510"] = float((k510 & uy).sum() / k510.sum())
    # (iv) köşegen şerit kontrastı (artık)
    R = M - K1
    tb = 0.5 * (ince_kenar[:-1] + ince_kenar[1:])
    bsel = ince_kenar[:-1] >= 0.70 - 1e-9
    ser = (np.abs(orta[None, :] - (2.0 - tb[:, None])) <= 0.02 + 1e-12) | \
          (np.abs(orta[None, :] - tb[:, None]) <= 0.02 + 1e-12)
    A = np.abs(R[bsel])
    sm = ser[bsel]
    out["kontrast"] = float(A[sm].mean() / A[~sm].mean())
    out["serit_hucre"] = int(sm.sum())
    # ayrı ayrı iki şerit (bilgi)
    sk = (np.abs(orta[None, :] - (2.0 - tb[:, None])) <= 0.02 + 1e-12)[bsel]
    sy = (np.abs(orta[None, :] - tb[:, None]) <= 0.02 + 1e-12)[bsel]
    out["kontrast_katlanma"] = float(A[sk].mean() / A[~sm].mean()) \
        if sk.any() else float("nan")
    out["kontrast_yakin"] = float(A[sy].mean() / A[~sm].mean()) \
        if sy.any() else float("nan")
    # (vi) 8-bant yeniden kurulum (Σ|mix|²-ağırlıklı ince-bant toplamı)
    b8 = ONK["bant8_kenar"]
    z8 = np.zeros(8, complex)
    z8_r1 = np.zeros(8, complex)
    z8_dog = np.zeros(8, complex)
    for B in range(8):
        sel = (ince_kenar[:-1] >= b8[B] - 1e-9) & \
              (ince_kenar[1:] <= b8[B + 1] + 1e-9)
        w = pay[:NI][sel] / pay[:NI][sel].sum()
        z8[B] = np.sum(w[:, None] * M[sel])
        z8_r1[B] = np.sum(w[:, None] * K1[sel])
        z8_dog[B] = K[NI + B].sum()
    out["z8"] = z8
    out["z8_r1"] = z8_r1
    out["z8_art"] = z8 - z8_r1
    out["z8_dogrudan_fark"] = float(np.max(np.abs(z8 - z8_dog)))
    z = -z8.real
    za = -(z8 - z8_r1).real
    Bmin = int(np.argmin(z))
    out["Bmin"] = Bmin
    out["yukselis"] = float(z[7] - z[Bmin])
    out["art_pay_yukselis"] = float((za[7] - za[Bmin]) / (z[7] - z[Bmin])) \
        if z[7] - z[Bmin] > 0 else float("nan")
    # (v) H-188d
    if KH is not None:
        m120 = kenar[1:] <= 1.20 + 1e-9
        kr = kap[m120].reshape(-1, 2).sum(1)
        kh = -KH[-1].real
        out["kappa01_gercek"] = kr
        out["kappa01_ikiz"] = kh
        out["corr_d"] = pearson(kr, kh)
        Kr01 = K[-1][m120].reshape(-1, 2).sum(1)
        out["corr_d_abs"] = pearson(np.abs(Kr01), np.abs(KH[-1]))
    return out


if __name__ == "__main__":
    ONK = onkayit()
    print("=" * 78)
    print(f"188d / K3 ANALİZ  [on-kayit {ONK['zaman']} sha {ONK['sha256'][:12]}]")
    print("=" * 78)
    H = np.load(S188 / "harita_K_gercek.npz")
    Hk = np.load(S188 / "harita_K_Hkeskin.npz")
    K, Kr = H["K"], H["K_reps"]
    pay, payr = H["pay"], H["pay_reps"]
    S, Sr, SRe, SRer, W = H["S"], H["S_reps"], H["SRe"], H["SRe_reps"], H["W"]
    orta, kenar = H["orta"], H["kenar"]
    ince_kenar = H["ince_kenar"]
    KH, KHr = Hk["K"], Hk["K_reps"]

    tam = analiz(K, pay, S, SRe, W, orta, kenar, ince_kenar, ONK, KH=KH)
    reps = [analiz(Kr[j], payr[j], Sr[j], SRer[j], W, orta, kenar,
                   ince_kenar, ONK, KH=KHr[j]) for j in range(NJACK)]

    def se(ad):
        return float(jk_se([r[ad] for r in reps]))

    # gürültü-düzeltmeli rang-1 payı (İKİNCİL)
    M = K[:NI]
    trN = float(np.sum(jk_se(Kr[:, :NI].real) ** 2 + jk_se(Kr[:, :NI].imag) ** 2))
    s = np.linalg.svd(M, compute_uv=False)
    pay_nc = float(s[0] ** 2 / (np.sum(s ** 2) - trN))
    pay_gurultu = float(trN / np.sum(s ** 2))

    print("\n(i) RANG-1 (41×88 karmaşık K, ağırlıksız SVD):")
    print(f"   rang-1 varyans payı = {tam['pay1']:.4f} ± {se('pay1'):.4f}"
          f"   [ikinci σ²-payı {s[1]**2/np.sum(s**2):.4f}]")
    print(f"   İKİNCİL: gürültü (tr N) payı = {pay_gurultu:.4f}; "
          f"gürültü-düzeltmeli rang-1 payı = {pay_nc:.4f}")
    print("\n(ii) KORELASYONLAR (88 dilim):")
    print(f"   corr(|v|,S) = {tam['corr_vS']:.4f} ± {se('corr_vS'):.4f}")
    print(f"   corr(|v|,W) = {tam['corr_vW']:.4f} ± {se('corr_vW'):.4f}")
    dfr = [r["corr_vS"] - r["corr_vW"] for r in reps]
    print(f"   fark        = {tam['corr_vS']-tam['corr_vW']:.4f} ± "
          f"{float(jk_se(dfr)):.4f}")
    print(f"   İKİNCİL corr(Re v, S_Re) = {tam['corr_RevSRe']:.4f} ± "
          f"{se('corr_RevSRe'):.4f}")
    print("\n(iii) (1.05,1.10] κ-katkısının 1+τ_2/1+τ_3 pencerelerindeki payı:")
    print(f"   pay = {tam['uydu_pay_510']:.4f} ± {se('uydu_pay_510'):.4f}  "
          f"(sıfır beklentisi = {tam['uydu_sifir_510']:.2f})")
    print("\n(iv) KÖŞEGEN ŞERİT KONTRASTI (artık R, τ_b≥0.70, ±0.02):")
    print(f"   kontrast = {tam['kontrast']:.4f} ± {se('kontrast'):.4f} "
          f"({tam['serit_hucre']} şerit hücresi)  [katlanma "
          f"{tam['kontrast_katlanma']:.3f}, yakın {tam['kontrast_yakin']:.3f}]")
    print(f"   yükseliş Δz = {tam['yukselis']:.4f} ± {se('yukselis'):.4f} "
          f"(B_min={tam['Bmin']}); artık payı = {tam['art_pay_yukselis']:.4f}"
          f" ± {se('art_pay_yukselis'):.4f}")
    print("\n(v) H-188d (κ-profili gerçek vs ikiz, 34 dilim):")
    print(f"   corr = {tam['corr_d']:.4f} ± {se('corr_d'):.4f}   "
          f"İKİNCİL |K|-profil corr = {tam['corr_d_abs']:.4f} ± "
          f"{se('corr_d_abs'):.4f}")

    # (vi) K3 tablosu
    K2 = json.load(open(S187 / "K2_zeta.json"))
    zf = np.array(K2["gercek"]["zeta_mod"][:8])
    zf_se = np.array(K2["gercek"]["se_mod"][:8])
    b8 = ONK["bant8_kenar"]
    z = -tam["z8"].real
    zr1 = -tam["z8_r1"].real
    za = -tam["z8_art"].real
    z_se = jk_se([-r["z8"].real for r in reps])
    zr1_se = jk_se([-r["z8_r1"].real for r in reps])
    za_se = jk_se([-r["z8_art"].real for r in reps])
    print("\n(vi) K3 — 8-BANT ζ^{≤1.30} YENİDEN KURULUMU (z = −Re ζ):")
    print(f"{'bant':>10} {'|ζ| tam(187c)':>15} {'|ζ≤1.30|':>9} {'açı°':>8} "
          f"{'z':>14} {'z_rang1':>14} {'z_artık':>14} {'z/|ζ|tam':>9}")
    k3 = []
    for B in range(8):
        a = np.degrees(np.angle(tam["z8"][B]))
        print(f"{b8[B]:.2f}-{b8[B+1]:.2f} {zf[B]:8.4f}±{zf_se[B]:.4f} "
              f"{abs(tam['z8'][B]):9.4f} {a:+8.2f} "
              f"{z[B]:7.4f}±{z_se[B]:.4f} {zr1[B]:7.4f}±{zr1_se[B]:.4f} "
              f"{za[B]:7.4f}±{za_se[B]:.4f} {z[B]/zf[B]:9.3f}")
        k3.append({"bant": f"{b8[B]:.2f}-{b8[B+1]:.2f}",
                   "zeta_tam_187c": float(zf[B]), "se_zeta_tam": float(zf_se[B]),
                   "zeta130_mod": float(abs(tam["z8"][B])), "zeta130_aci": float(a),
                   "z": float(z[B]), "se_z": float(z_se[B]),
                   "z_r1": float(zr1[B]), "se_z_r1": float(zr1_se[B]),
                   "z_art": float(za[B]), "se_z_art": float(za_se[B])})
    print(f"   ince→8 bant toplama kontrolü maks|Δ| = {tam['z8_dogrudan_fark']:.1e}")

    # ---- u(τ_b) profili, |v| tepeleri ----
    tb = 0.5 * (ince_kenar[:-1] + ince_kenar[1:])
    u = tam["u"]
    tp = ONK["uydular"]["tau_p"]
    etk = {"Bragg": 1.0}
    for p in ("2", "3", "5", "7"):
        etk[f"1+τ{p}"] = 1 + tp[p]
        etk[f"1−τ{p}"] = 1 - tp[p]

    def en_yakin(t):
        ad = min(etk, key=lambda k: abs(etk[k] - t))
        return ad, etk[ad] - t

    absv = tam["absv"]
    sira = np.argsort(-absv)[:12]
    print("\n|v(τ')| EN BÜYÜK 12 DİLİM (en yakın uydu):")
    tepe = []
    for i in sira:
        ad, d = en_yakin(orta[i])
        print(f"   τ'={orta[i]:.4f} |v|={absv[i]:.4f} ∠{np.degrees(np.angle(tam['v'][i])):+7.1f}°"
              f"  κ={tam['kappa'][i]:+.5f}  S={S[i]:.3e} W={W[i]:.3e}  "
              f"[{ad} {d:+.4f}]")
        tepe.append({"tau": float(orta[i]), "absv": float(absv[i]),
                     "kappa": float(tam["kappa"][i]), "uydu": ad,
                     "uzaklik": float(d)})
    print("\nu(τ_b) (rang-1 bant çarpanı; |u|, açı):")
    for B in range(0, NI, 5):
        print(f"   τ_b={tb[B]:.3f} |u|={abs(u[B]):.4f} ∠{np.degrees(np.angle(u[B])):+.1f}°")

    # ---- HÜKÜMLER ----
    Ha = ONK["H_188a"]
    k_i = tam["pay1"] >= 0.80
    k_ii = (tam["corr_vS"] >= 0.80) and \
        (tam["corr_vS"] - tam["corr_vW"] >= 0.15)
    k_iii = tam["uydu_pay_510"] >= 0.60
    olum_a = (tam["pay1"] < 0.60) or (tam["corr_vS"] < 0.50)
    Hc = tam["corr_vW"] >= tam["corr_vS"]
    if Hc:
        hA = "ÖLDÜ (H-188c tuttu: tarak köprüsü ölür)"
    elif k_i and k_ii and k_iii:
        hA = "MÜHÜR"
    elif olum_a:
        hA = "ÖLDÜ"
    else:
        hA = "KAYIT"
    hB_k = tam["kontrast"] >= 3.0
    hB_p = (not np.isnan(tam["art_pay_yukselis"])) and \
        tam["art_pay_yukselis"] >= 0.60
    hB = "MÜHÜR" if (hB_k and hB_p) else ("ÖLDÜ" if not hB_k else "KAYIT")
    hC = "MÜHÜR (yayvan taban)" if Hc else "ÖLDÜ"
    hD = "MÜHÜR" if tam["corr_d"] >= 0.80 else "ÖLDÜ"
    print("\nHÜKÜMLER (ön-kayıt eşikleri; kurtarma yok):")
    print(f"  H-188a: (i) {tam['pay1']:.3f}≥0.80 {'EVET' if k_i else 'HAYIR'};"
          f" (ii) {'EVET' if k_ii else 'HAYIR'}; (iii) {tam['uydu_pay_510']:.3f}"
          f"≥0.60 {'EVET' if k_iii else 'HAYIR'}; ölüm koşulu "
          f"{'EVET' if olum_a else 'HAYIR'} → {hA}")
    print(f"  H-188b: kontrast {tam['kontrast']:.3f}≥3 {'EVET' if hB_k else 'HAYIR'};"
          f" artık payı {tam['art_pay_yukselis']:.3f}≥0.60 "
          f"{'EVET' if hB_p else 'HAYIR'} → {hB}")
    print(f"  H-188c: corr(|v|,W)={tam['corr_vW']:.3f} ≥ corr(|v|,S)="
          f"{tam['corr_vS']:.3f}? {'EVET' if Hc else 'HAYIR'} → {hC}")
    print(f"  H-188d: corr={tam['corr_d']:.3f}≥0.80? → {hD}")

    def J(x):
        return [float(t) for t in np.ravel(x)]

    sonuc = {
        "sha_onkayit": ONK["sha256"],
        "i": {"pay1": tam["pay1"], "se": se("pay1"),
              "sigma2_pay": J(s ** 2 / np.sum(s ** 2)),
              "ikincil_gurultu_pay": pay_gurultu,
              "ikincil_pay_nc": pay_nc},
        "ii": {"corr_vS": tam["corr_vS"], "se_vS": se("corr_vS"),
               "corr_vW": tam["corr_vW"], "se_vW": se("corr_vW"),
               "fark": tam["corr_vS"] - tam["corr_vW"],
               "se_fark": float(jk_se(dfr)),
               "ikincil_corr_RevSRe": tam["corr_RevSRe"],
               "se_RevSRe": se("corr_RevSRe")},
        "iii": {"pay": tam["uydu_pay_510"], "se": se("uydu_pay_510"),
                "sifir": tam["uydu_sifir_510"]},
        "iv": {"kontrast": tam["kontrast"], "se": se("kontrast"),
               "katlanma": tam["kontrast_katlanma"],
               "yakin": tam["kontrast_yakin"],
               "yukselis": tam["yukselis"], "se_yukselis": se("yukselis"),
               "Bmin": tam["Bmin"], "art_pay": tam["art_pay_yukselis"],
               "se_art_pay": se("art_pay_yukselis")},
        "v": {"corr_d": tam["corr_d"], "se": se("corr_d"),
              "ikincil_corr_abs": tam["corr_d_abs"],
              "se_abs": se("corr_d_abs"),
              "kappa01_gercek": J(tam["kappa01_gercek"]),
              "kappa01_ikiz": J(tam["kappa01_ikiz"])},
        "vi": k3,
        "tepeler_v": tepe,
        "u_mod": J(np.abs(u)), "u_aci": J(np.degrees(np.angle(u))),
        "v_mod": J(absv), "v_aci": J(np.degrees(np.angle(tam["v"]))),
        "kappa": J(tam["kappa"]), "S": J(S), "W": J(W), "SRe": J(SRe),
        "orta": J(orta), "tb": J(tb),
        "hukum": {"H-188a": hA, "H-188b": hB, "H-188c": hC, "H-188d": hD,
                  "K3": "KAYIT"},
    }
    json.dump(sonuc, open(S188 / "K3_analiz.json", "w"), indent=1,
              ensure_ascii=False)
    np.savez_compressed(S188 / "K3_analiz.npz", K1=svd1(K[:NI])["K1"],
                        u=u, v=tam["v"])
    print(f"\n-> {S188/'K3_analiz.json'}  BİTTİ")
