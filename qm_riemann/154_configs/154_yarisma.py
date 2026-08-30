"""
154 — KAPALI-FORM YARIŞI
========================
Girdi: 154_gercek.py'nin iki penceresi × iki bant ızgarası (JSON).
Her adayda EN ÇOK 2 serbest parametre (8 veri noktası; aşırı-fit yasak).

HATA MODELİ (dürüst, muhafazakâr):
  σ_i = max( jackknife_i , |R_son − R_orta|/2 , 0.02 )
  Jackknife bant-içi çizgi-seçim saçılması; iki pencere farkı bağımsız
  üst sınır; 0.02 taban ise R-ızgarasının kendi çözünürlüğü.
SKOR: χ² = Σ((R_i − f_i)/σ_i)² , AIC = χ² + 2k.

İKİ TAHMİNCİ: R_emp (ampirik ağırlıklı faz eşlemesi — 150'nin tanımı) ve
R_gauss = φ/(A²σΔ²) (Gauss limiti). Yüksek τ'da ayrışırlar; ayrışma
sistematik hatanın kendisidir ve rapora girer.
"""
import json
import numpy as np
from pathlib import Path
from scipy.special import erfcx
from scipy.optimize import minimize

OUT = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/154")
TWO_PI = 2 * np.pi
FIT_TAU = [0.5375, 0.585, 0.66, 0.74]     # 0.815 fite GİRMEZ (gürültülü)


def R_erfc(tau, tc=0.68, dl=0.125, p=2.0):
    """−(1/2π)·dln[½erfc((τ−tc)/Δ)]^p/dτ = (p/(π^{3/2}Δ))/erfcx(u)."""
    u = (np.asarray(tau, float) - tc) / dl
    return p / (np.pi**1.5 * dl) / erfcx(u)


def yukle(pen, bs):
    j = json.loads((OUT / f"R_gercek_{pen}_{bs}.json").read_text())
    b = [x for x in j["bantlar"] if x.get("olculdu")]
    return j, {round(x["tau"], 4): x for x in b}


def Rg(rec, sA2):
    """Gauss-limit tahmincisi φ/(A²σΔ²)."""
    A = TWO_PI * rec["tau"]
    return rec["phi"] / (A * A * sA2)


# ---------- ρ'dan türev: ince ızgaraya log-kübik, analitik d/dτ ----------
class RhoTurev:
    def __init__(self, tau, rho, derece=3):
        self.t = np.asarray(tau, float)
        self.c = np.polyfit(self.t, np.log(np.asarray(rho, float)), derece)
        self.d = np.polyder(self.c)
        self.derece = derece

    def dln(self, tau):
        return np.polyval(self.d, np.asarray(tau, float))

    def lnrho(self, tau):
        return np.polyval(self.c, np.asarray(tau, float))


def yarist(tau, R, sig, basl, adaylar):
    print(f"\n=== {basl}  (n={len(R)}) ===")
    print(f"{'aday':54s} {'k':>2s} {'χ²':>9s} {'AIC':>9s} {'RMS':>6s}"
          f"   parametreler")
    sat = []
    for ad, k, f, x0s, pad in adaylar:
        if k == 0:
            th = np.array([]); r = (R - f(th, tau)) / sig
            chi = float(r @ r)
        else:
            best = None
            for x0 in x0s:
                res = minimize(lambda th: float(np.sum(
                    ((R - f(th, tau)) / sig)**2)), np.array(x0, float),
                    method="Nelder-Mead",
                    options=dict(xatol=1e-11, fatol=1e-13, maxiter=40000,
                                 maxfev=40000))
                if best is None or res.fun < best.fun:
                    best = res
            th = best.x; chi = float(best.fun)
        art = R - f(th, tau)
        ps = "  ".join(f"{n}={v:.4g}" for n, v in zip(pad, th)) or "—"
        print(f"{ad:54s} {k:2d} {chi:9.2f} {chi+2*k:9.2f} "
              f"{float(np.sqrt(np.mean(art**2))):6.3f}   {ps}")
        sat.append(dict(ad=ad, k=k, chi2=chi, aic=chi + 2 * k,
                        rms=float(np.sqrt(np.mean(art**2))),
                        theta=[float(v) for v in th],
                        artik=[float(v) for v in art]))
    print(f"\n  artık (R_ölç − model), τ = "
          f"{['%.4f' % t for t in tau[:len(tau)//2 or len(tau)]]} …")
    for s in sat:
        print(f"   {s['ad'][:46]:46s} " +
              " ".join(f"{a:+6.3f}" for a in s["artik"]))
    return sat


if __name__ == "__main__":
    js, S = yukle("son", "std"); jo, O = yukle("orta", "std")
    jsi, Si = yukle("son", "ince"); joi, Oi = yukle("orta", "ince")
    sA2s, sA2o = js["sA2"], jo["sA2"]

    print("=" * 78)
    print("1) İKİ PENCERE — L-DEĞİŞMEZLİK (std bantlar)")
    print(f"   L: son {js['L']:.4f}   orta {jo['L']:.4f}   "
          f"(ΔL = {js['L']-jo['L']:.3f})")
    print(f"{'τ̄':>7s} {'R_son':>7s} {'±jk':>6s} {'R_orta':>7s} {'±jk':>6s} "
          f"{'|fark|':>7s} | {'Rg_son':>7s} {'Rg_orta':>8s} | "
          f"{'ρ_son':>7s} {'ρ_orta':>7s} | {'n_eff':>7s}")
    tum = FIT_TAU + [0.815]
    for t in tum:
        a, b = S[t], O[t]
        print(f"{t:7.4f} {a['R']:7.3f} {a['sR_jk']:6.3f} {b['R']:7.3f} "
              f"{b['sR_jk']:6.3f} {abs(a['R']-b['R']):7.3f} | "
              f"{Rg(a,sA2s):7.3f} {Rg(b,sA2o):8.3f} | {a['rho']:7.4f} "
              f"{b['rho']:7.4f} | {a['neff']:7.0f}")

    print("\n2) İNCE IZGARA — R(τ) ve ρ(τ) aynı koşudan")
    print(f"{'τ̄':>7s} {'R_son':>7s} {'±jk':>6s} {'R_orta':>7s} {'±jk':>6s} "
          f"{'ρ_son':>8s} {'±':>7s} {'ρ_orta':>8s} {'±':>7s}")
    ft = sorted(Si)
    for t in ft:
        a, b = Si[t], Oi.get(t)
        print(f"{t:7.4f} {a['R']:7.3f} {a['sR_jk']:6.3f} "
              f"{b['R']:7.3f} {b['sR_jk']:6.3f} {a['rho']:8.4f} "
              f"{a['sRho_jk']:7.4f} {b['rho']:8.4f} {b['sRho_jk']:7.4f}")

    RT_s = RhoTurev(ft, [Si[t]["rho"] for t in ft])
    RT_o = RhoTurev(ft, [Oi[t]["rho"] for t in ft])
    RT_b = RhoTurev(ft, [0.5 * (Si[t]["rho"] + Oi[t]["rho"]) for t in ft])
    print(f"\n   lnρ(τ) kübik: son  {np.poly1d(RT_s.c)}")
    print(f"   lnρ(τ) kübik: orta {np.poly1d(RT_o.c)}")

    print("\n3) ADAY (f): R = −(1/π)·dlnρ/dτ = 2·(−dlnρ/dA)  [0 PARAMETRE]")
    print("   (2 çarpanı: bond İKİ adımlıdır — dsΔ=(ds_n+ds_{n+1})/2 ⇒")
    print("    adım başına e^{−A·R₁·ds} iki kez ⇒ ölçülen R = 2R₁)")
    print(f"{'τ':>7s} {'R_ölç(s)':>9s} {'R_ölç(o)':>9s} {'−dlnρ/dA(s)':>12s} "
          f"{'×2':>7s} {'oran R/(−dlnρ/dA)':>19s}")
    for t in ft:
        d1s = -RT_s.dln(t) / TWO_PI
        print(f"{t:7.4f} {Si[t]['R']:9.3f} {Oi[t]['R']:9.3f} {d1s:12.3f} "
              f"{2*d1s:7.3f} {Si[t]['R']/d1s:19.3f}")

    # ---------------- adaylar ----------------
    ADAYLAR = [
        ("a  R = c  (sabit)", 1,
         lambda th, t: np.full_like(np.asarray(t, float), th[0]),
         [[1.0], [2.0]], ["c"]),
        ("b  R = c·(τ−τ₀)  (doğrusal)", 2,
         lambda th, t: th[0] * (np.asarray(t, float) - th[1]),
         [[8.0, 0.49], [5.0, 0.45], [12.0, 0.52]], ["c", "τ₀"]),
        ("c  R = R∞·(1−e^{−(τ−0.52)/w})", 2,
         lambda th, t: th[0] * (1 - np.exp(-(np.asarray(t, float) - 0.52)
                                           / max(th[1], 1e-4))),
         [[3.0, 0.12], [2.5, 0.08], [5.0, 0.25], [2.2, 0.05]], ["R∞", "w"]),
        ("d  −dlnρ_erfc/dA (0.68/0.125, p=2)   [0 PAR]", 0,
         lambda th, t: R_erfc(t), None, []),
        ("d2 2·(−dlnρ_erfc/dA)  (iki-adım)     [0 PAR]", 0,
         lambda th, t: 2 * R_erfc(t), None, []),
        ("d′ s·(−dlnρ_erfc/dA)  (ölçek serbest)", 1,
         lambda th, t: th[0] * R_erfc(t), [[1.0], [2.0], [0.5]], ["s"]),
        ("d″ erfc ama τ_c,Δ serbest (p=2)", 2,
         lambda th, t: R_erfc(t, th[0], max(abs(th[1]), 1e-3)),
         [[0.68, 0.125], [0.85, 0.25], [1.0, 0.4], [0.75, 0.2],
          [1.3, 0.6]], ["τ_c", "Δ"]),
        ("e  R = k·A²σΔ² = k(2πτ)²σΔ²", 1,
         lambda th, t: th[0] * (TWO_PI * np.asarray(t, float))**2 * sA2s,
         [[0.6], [1.0]], ["k"]),
        ("f  −(1/π)·dlnρ_ÖLÇÜLEN/dτ           [0 PAR]", 0,
         lambda th, t: -RT_b.dln(t) / np.pi, None, []),
        ("f1 −(1/2π)·dlnρ_ÖLÇÜLEN/dτ (tek adım)[0 PAR]", 0,
         lambda th, t: -RT_b.dln(t) / TWO_PI, None, []),
        ("f′ s·(−(1/2π)dlnρ_ÖLÇÜLEN/dτ)  (s serbest)", 1,
         lambda th, t: -th[0] * RT_b.dln(t) / TWO_PI,
         [[2.0], [1.0]], ["s"]),
        ("g  R = c(τ−τ₀)/τ²   [φ doğrusal + kinematik A⁻²]", 2,
         lambda th, t: th[0] * (np.asarray(t, float) - th[1])
         / np.asarray(t, float)**2, [[4.0, 0.50], [3.0, 0.45],
                                     [5.0, 0.52]], ["c", "τ₀"]),
        ("g½ R = c(τ−1/2)/τ²  (τ₀ ≡ 1/2 SABİT)", 1,
         lambda th, t: th[0] * (np.asarray(t, float) - 0.5)
         / np.asarray(t, float)**2, [[4.0], [3.0]], ["c"]),
        ("c′ R∞(1−e^{−(τ−τ₀)/w}), τ₀=0.5139 ÖLÇÜLEN", 2,
         lambda th, t: th[0] * (1 - np.exp(-(np.asarray(t, float) - 0.5139)
                                           / max(th[1], 1e-4))),
         [[2.5, 0.10], [2.2, 0.06], [4.0, 0.2]], ["R∞", "w"]),
        ("gτ R = c(τ−τ₀)/τ²  τ₀=0.5139 ÖLÇÜLEN", 1,
         lambda th, t: th[0] * (np.asarray(t, float) - 0.5139)
         / np.asarray(t, float)**2, [[5.0], [4.0]], ["c"]),
    ]

    tau = np.array(FIT_TAU)
    Rs = np.array([S[t]["R"] for t in FIT_TAU])
    Ro = np.array([O[t]["R"] for t in FIT_TAU])
    Gs = np.array([Rg(S[t], sA2s) for t in FIT_TAU])
    Go = np.array([Rg(O[t], sA2o) for t in FIT_TAU])
    jsk = np.array([S[t]["sR_jk"] for t in FIT_TAU])
    jok = np.array([O[t]["sR_jk"] for t in FIT_TAU])
    fark = np.abs(Rs - Ro)
    ss = np.maximum(np.maximum(jsk, fark / 2), 0.02)
    so = np.maximum(np.maximum(jok, fark / 2), 0.02)

    print("\n" + "=" * 78)
    s8 = yarist(np.concatenate([tau, tau]), np.concatenate([Rs, Ro]),
                np.concatenate([ss, so]), "8 NOKTA — R_emp (150 tanımı)",
                ADAYLAR)
    fg = np.abs(Gs - Go)
    sg = np.maximum(np.maximum(jsk, fg / 2), 0.02)
    s8g = yarist(np.concatenate([tau, tau]), np.concatenate([Gs, Go]),
                 np.concatenate([sg, sg]),
                 "8 NOKTA — R_gauss = φ/(A²σΔ²) (ikinci tahminci)", ADAYLAR)
    # ince ızgara: 9 τ × 2 pencere = 18 nokta, bağımsız çapraz sınav
    tf = np.array(ft)
    Rfs = np.array([Si[t]["R"] for t in ft]); Rfo = np.array([Oi[t]["R"]
                                                              for t in ft])
    jfs = np.array([Si[t]["sR_jk"] for t in ft])
    jfo = np.array([Oi[t]["sR_jk"] for t in ft])
    ff = np.abs(Rfs - Rfo)
    sfs = np.maximum(np.maximum(jfs, ff / 2), 0.02)
    sfo = np.maximum(np.maximum(jfo, ff / 2), 0.02)
    s18 = yarist(np.concatenate([tf, tf]), np.concatenate([Rfs, Rfo]),
                 np.concatenate([sfs, sfo]),
                 "18 NOKTA — ince ızgara (çapraz sınav)", ADAYLAR)

    # ---- 6) SİSTEMATİKLİ hata modeli: taban taramasının kendi yayılımı ----
    print("\n" + "=" * 78)
    print("6) SİSTEMATİKLİ YARIŞ — σ² = σ_jk² + σ_pencere² + σ_TABAN²")
    print("   σ_TABAN(τ) = |R@0.58 − R@0.46|/2  (ölçülmüş konvansiyon")
    print("   duyarlılığı; τ<0.60'ta 0.46↔0.52 farkının iki katı).")
    print("   Bant τ>0.76 ELENDİ (çözünmeyen çizgi bölgesi — 152 uyarısı).")
    try:
        TBj = {k: {round(x["tau"], 4): x for x in json.loads(
            (OUT / f"R_taban_{k}.json").read_text())["bantlar"]
            if x.get("olculdu")} for k in ("0.46", "0.52", "0.58")}
        tk = [t for t in ft if t <= 0.76]
        syst = []
        for t in tk:
            if t in TBj["0.58"] and t in TBj["0.46"]:
                syst.append(abs(TBj["0.58"][t]["R"] - TBj["0.46"][t]["R"]) / 2)
            else:
                syst.append(2 * abs(TBj["0.52"][t]["R"] - TBj["0.46"][t]["R"]))
        syst = np.array(syst)
        tk = np.array(tk)
        Rk_s = np.array([Si[t]["R"] for t in tk])
        Rk_o = np.array([Oi[t]["R"] for t in tk])
        jk_s = np.array([Si[t]["sR_jk"] for t in tk])
        jk_o = np.array([Oi[t]["sR_jk"] for t in tk])
        pen = np.abs(Rk_s - Rk_o) / 2
        sg_s = np.sqrt(jk_s**2 + pen**2 + syst**2)
        sg_o = np.sqrt(jk_o**2 + pen**2 + syst**2)
        print(f"{'τ':>7s} {'R':>7s} {'σ_jk':>6s} {'σ_pen':>6s} "
              f"{'σ_taban':>8s} {'σ_top':>7s} {'%':>5s}")
        for i, t in enumerate(tk):
            print(f"{t:7.3f} {Rk_s[i]:7.3f} {jk_s[i]:6.3f} {pen[i]:6.3f} "
                  f"{syst[i]:8.3f} {sg_s[i]:7.3f} "
                  f"{100*sg_s[i]/Rk_s[i]:5.1f}")
        sS = yarist(np.concatenate([tk, tk]),
                    np.concatenate([Rk_s, Rk_o]),
                    np.concatenate([sg_s, sg_o]),
                    "16 NOKTA — sistematikli (τ≤0.76)", ADAYLAR)
        print(f"\n   serbestlik derecesi: n=16;  χ²/dof kabul çizgisi ≈ 1")
        for s in sorted(sS, key=lambda x: x["aic"])[:6]:
            print(f"   {s['ad'][:50]:50s} χ²/dof = "
                  f"{s['chi2']/(16-s['k']):6.2f}")
    except FileNotFoundError as e:
        sS = None
        print(f"   (taban JSON yok: {e})")

    # ---- φ düzeyinde: ilkel gözlenebilir R değil φ'dir ----
    print("\n" + "=" * 78)
    print("4) İLKEL GÖZLENEBİLİR φ_Γ(τ) — sıfırı nerede?")
    print("   (a) tüm ince ızgaraya doğrusal fit — φ EĞRİ olduğu için")
    print("       τ₀ aralığa bağlı; YANILTICI, yalnız kayıt:")
    for et, D in (("son ", Si), ("orta", Oi)):
        t = np.array(sorted(D)); p = np.array([D[x]["phi"] for x in sorted(D)])
        c1_, c0_ = np.polyfit(t, p, 1)
        art = p - (c1_ * t + c0_)
        print(f"       {et}: φ = {c1_:.3f}·(τ − {-c0_/c1_:.4f})   "
              f"rms artık {np.sqrt(np.mean(art**2)):.4f}")
    print("   (b) YEREL sıfır — en alttaki üç banda parabol (dürüst tahmin):")
    for et, D in (("son(0.52) ", Si), ("orta(0.52)", Oi)):
        t = np.array(sorted(D))[:3]; p = np.array(
            [D[x]["phi"] for x in sorted(D)])[:3]
        k = np.polyfit(t, p, 2); r = np.roots(k)
        r = [x.real for x in r if abs(x.imag) < 1e-9 and 0.40 < x.real < 0.56]
        print(f"       {et}: τ₀ = {r[0]:.4f}  (bantlar "
              f"{t[0]:.3f}/{t[1]:.3f}/{t[2]:.3f}, EKSTRAPOLASYON)")
    try:
        D = {round(x["tau"], 4): x for x in json.loads(
            (OUT / "R_taban_0.46.json").read_text())["bantlar"]
            if x.get("olculdu")}
        t = np.array(sorted(D))[:3]; p = np.array(
            [D[x]["phi"] for x in sorted(D)])[:3]
        k = np.polyfit(t, p, 2); r = np.roots(k)
        r = [x.real for x in r if abs(x.imag) < 1e-9 and 0.40 < x.real < 0.56]
        print(f"       taban0.46 : τ₀ = {r[0]:.4f}  (bantlar "
              f"{t[0]:.3f}/{t[1]:.3f}/{t[2]:.3f} — φ burada NEGATİF "
              f"[{p[0]:+.3f},{p[1]:+.3f}], yani DOĞRUDAN ÖLÇÜM)")
        print(f"       ⇒ φ'nin sıfırı TABANLA KAYMIYOR: 0.46 ve 0.52 "
              f"tabanları aynı τ₀'ı veriyor.")
        # R'nin kendi sıfırı: R=0'da M=M_emp, arg M_emp ≠ 0 ⇒ R'nin
        # sıfırı φ'nin sıfırından biraz YUKARIDA. Doğrudan ölçüm:
        tt3 = np.array(sorted(D))[:4]
        rr3 = np.array([D[x]["R"] for x in sorted(D)])[:4]
        if np.any(rr3 < 0) and np.any(rr3 > 0):
            k = np.polyfit(tt3, rr3, 2); rt = np.roots(k)
            rt = [x.real for x in rt if abs(x.imag) < 1e-9
                  and 0.44 < x.real < 0.58]
            print(f"   (c) R'NİN KENDİ SIFIRI (doğrudan, taban 0.46): "
                  f"τ_R0 = {rt[0]:.4f}")
            print(f"       ölçülen R: " + "  ".join(
                f"{a:.3f}@{b:.3f}" for b, a in zip(tt3, rr3)))
            print(f"       [karşılaştır: aday (g)'nin τ≥0.5375 bantlarından "
                  f"FİT ettiği τ₀]")
    except FileNotFoundError:
        pass

    print("\n5) TABAN TARAMASI — R(τ) (τ−taban)'ın mı fonksiyonu?")
    try:
        TB = {}
        for tb in ("0.46", "0.52", "0.58"):
            jj = json.loads((OUT / f"R_taban_{tb}.json").read_text())
            TB[tb] = {round(x["tau"], 4): x for x in jj["bantlar"]
                      if x.get("olculdu")}
        ortak = sorted(set(TB["0.46"]) & set(TB["0.52"]) & set(TB["0.58"]))
        print(f"{'τ':>7s} " + "".join(f"{'R@'+k:>9s}{'ρ@'+k:>9s}"
                                      for k in TB))
        for t in ortak:
            print(f"{t:7.4f} " + "".join(
                f"{TB[k][t]['R']:9.3f}{TB[k][t]['rho']:9.4f}" for k in TB))
        print("   (τ−taban) motifi doğruysa taban 0.58 sütunu τ≈0.60'ta "
              "R≈0.4 vermeliydi.")
        for tb in TB:
            t = np.array(sorted(TB[tb])); p = np.array(
                [TB[tb][x]["phi"] for x in sorted(TB[tb])])
            m = t >= 0.53
            if m.sum() >= 3:
                c1_, c0_ = np.polyfit(t[m], p[m], 1)
                print(f"   taban {tb}: φ sıfırı τ₀ = {-c0_/c1_:.4f} "
                      f"(eğim {c1_:.3f}, τ≥0.53 bantlarından)")
    except FileNotFoundError as e:
        print(f"   (taban taraması JSON'u yok: {e})")

    json.dump(dict(
        tau=FIT_TAU, R_son=Rs.tolist(), R_orta=Ro.tolist(),
        Rg_son=Gs.tolist(), Rg_orta=Go.tolist(),
        sig_son=ss.tolist(), sig_orta=so.tolist(),
        ince_tau=ft, ince_R_son=Rfs.tolist(), ince_R_orta=Rfo.tolist(),
        ince_sig_son=sfs.tolist(), ince_sig_orta=sfo.tolist(),
        ince_rho_son=[Si[t]["rho"] for t in ft],
        ince_rho_orta=[Oi[t]["rho"] for t in ft],
        lnrho_kubik_son=RT_s.c.tolist(), lnrho_kubik_orta=RT_o.c.tolist(),
        skor8=s8, skor8_gauss=s8g, skor18=s18, skor_syst=sS),
        open(OUT / "yarisma.json", "w"), indent=1)
    print("\n-> yarisma.json")
