"""
173f — EK TANILAR (ön-mühürsüz KOKLAMALAR; iddia değil, not)
=============================================================
Yeni ölçüm YOK; `173/EGRI.json` + `167/insa_*.json` üstünde çalışır.
Bu betiğin çıktıları **ÖN-KAYITLI DEĞİLDİR** ve hüküm sayılmazlar;
173'ün ölçtükleri arasında göze çarpan düzenlilikleri 174'e NOT olarak
bırakır. (172 §6.4'ün "ön-mühürsüz koklamalar açıkça etiketlenir"
kuralı.)

T1  M(λ)'nın SAĞLIKLI SEKİZ NOKTALI log-λ kuadratiği (L145 dışarıda —
    169'un R_bant filtresini düşürdüğü için, gevşetme değil dışlama).
T2  σ marjinallerinin iç oranları ve yerel kuvvet üsteli p(λ).
T3  SADAKAT SINIRI: R_bant(λ) 1'i nerede kesiyor?
T4  θ(λ) = 1 geçişi (θ ≡ KALİB/g_cal; θ = 1 ⟺ artık-alan eşleşmesi
    SIFIR, 172 §G2.0).
T5  c(λ) ve 4/π².

Kullanım: 173f_ek_tanilar.py
Çıktı:    scratchpad/173/EK.json
"""
import glob
import json
import os

import numpy as np

SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
       "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
SCR173 = SCR + "/173"
C_HIP = 4.0 / np.pi ** 2


def main():
    E = json.load(open(SCR173 + "/EGRI.json"))
    D9 = E["defter"]
    sira = [g for g in ["L040", "L050", "L060", "L070", "L085", "Hkeskin",
                        "L115", "L130", "L140", "L145"] if g in D9]
    lam = np.array([D9[g]["lam"] for g in sira])
    sag = np.array([D9[g]["saglikli"] for g in sira])
    M = np.array([D9[g]["M"] for g in sira])
    sM = np.array([D9[g]["sM"] for g in sira])
    c = np.array([D9[g]["c"] for g in sira])
    sc = np.array([D9[g]["s_tot"] for g in sira])
    th = np.array([D9[g]["theta"] for g in sira])
    Rlo = np.array([min(D9[g]["R_bant"]) for g in sira])
    Rhi = np.array([max(D9[g]["R_bant"]) for g in sira])
    sds = np.array([D9[g]["sigds"] for g in sira])
    sX = np.array([D9[g]["sigX"] for g in sira])
    sC = np.array([D9[g]["sigC"] for g in sira])
    out = {}

    print("=" * 96)
    print("173f — EK TANILAR  [ÖN-MÜHÜRSÜZ KOKLAMALAR — iddia değil]")
    print("=" * 96)

    # ---- T1: sağlıklı gazların M kuadratiği ---------------------------
    print("\nT1 — M(λ)'nın SAĞLIKLI noktalarla log-λ kuadratiği")
    out["T1"] = {}
    for et, m in (("λ ≤ 1.30 (sağlıklı)", sag & (lam <= 1.31)),
                  ("BÜTÜN sağlıklılar", sag)):
        if m.sum() < 4:
            continue
        p = np.polyfit(np.log(lam[m]), np.log(M[m]), 2)
        art = 100 * (M[m] / np.exp(np.polyval(p, np.log(lam[m]))) - 1)
        p1 = np.polyfit(np.log(lam[m]), np.log(M[m]), 1)
        art1 = 100 * (M[m] / np.exp(np.polyval(p1, np.log(lam[m]))) - 1)
        tepe = float(np.exp(-p[1] / (2 * p[0])))
        print("  [%s]  uyum noktaları: %s"
              % (et, ", ".join("%.2f" % x for x in lam[m])))
        print("    log M = %+.5f (logλ)² %+.5f logλ %+.5f   "
              "(durgun nokta λ = %.4f, M = %.4f)"
              % (p[0], p[1], p[2], tepe,
                 float(np.exp(np.polyval(p, np.log(tepe))))))
        print("    artıklar (%%): %s   rms %.2f%%   (ölçüm hatası ~%.2f%%)"
              % (" ".join("%+5.2f" % x for x in art),
                 float(np.sqrt((art ** 2).mean())),
                 float(100 * np.mean(sM[m] / M[m]))))
        print("    [karş.] saf kuvvet yasası M = %.5f·λ^(%+.5f): rms %.2f%%"
              % (np.exp(p1[1]), p1[0], float(np.sqrt((art1 ** 2).mean()))))
        d = dict(kat=[float(x) for x in p], artik=[float(x) for x in art],
                 rms=float(np.sqrt((art ** 2).mean())), tepe=tepe,
                 kuvvet=[float(np.exp(p1[1])), float(p1[0])],
                 kuvvet_rms=float(np.sqrt((art1 ** 2).mean())),
                 lam=[float(x) for x in lam[m]])
        out["T1"][et] = d
        if et.startswith("λ ≤"):        # figür bu uyumu kullanır
            out["T1"].update(d)

    # ---- T2: marjinallerin iç oranları --------------------------------
    print("\nT2 — MARJİNALLERİN İÇ ORANLARI ve YEREL KUVVET ÜSTELİ")
    print("  λ      σ_ds      σ_X̃/σ_ds   σ_Ĉ/σ_ds   p = dlogσ_ds/dlogλ")
    pp = np.diff(np.log(sds)) / np.diff(np.log(lam))
    for i, l in enumerate(lam):
        s = ("%+.4f (→%.2f)" % (pp[i], lam[i + 1])) if i < len(pp) else ""
        print("  %.2f   %.5f   %.5f    %.5f   %s"
              % (l, sds[i], sX[i] / sds[i], sC[i] / sds[i], s))
    out["T2"] = dict(oran_X=[float(x) for x in sX / sds],
                     oran_C=[float(x) for x in sC / sds],
                     p=[float(x) for x in pp])

    # ---- T3: sadakat sınırı -------------------------------------------
    print("\nT3 — SADAKAT SINIRI: R_bant(λ)")
    print("  λ      R_bant (min – maks, 5 hüküm bandı)   169 filtresi")
    for i, l in enumerate(lam):
        print("  %.2f   %.4f – %.4f                       %s"
              % (l, Rlo[i], Rhi[i], "GEÇTİ" if sag[i] else "**DÜŞTÜ**"))
    kes = []
    for i in range(len(lam) - 1):
        for nm, R in (("alt bant", Rlo), ("üst bant", Rhi)):
            if (R[i] - 1) * (R[i + 1] - 1) < 0:
                kes.append((nm, float(lam[i] + (lam[i + 1] - lam[i])
                                      * (R[i] - 1) / (R[i] - R[i + 1]))))
    for nm, k in kes:
        print("  R_bant = 1 geçişi (%s): λ ≈ %.4f" % (nm, k))
    kes98 = []
    for i in range(len(lam) - 1):
        if (Rhi[i] - 0.98) * (Rhi[i + 1] - 0.98) < 0:
            kes98.append(float(lam[i] + (lam[i + 1] - lam[i])
                               * (Rhi[i] - 0.98) / (Rhi[i] - Rhi[i + 1])))
    print("  169 filtresinin (0.98) üst-bant geçişi: λ ≈ %s"
          % ["%.4f" % k for k in kes98])
    out["T3"] = dict(Rlo=[float(x) for x in Rlo], Rhi=[float(x) for x in Rhi],
                     kesis1=[(n, k) for n, k in kes], kesis098=kes98)

    # ---- T4: θ = 1 geçişi ----------------------------------------------
    print("\nT4 — θ(λ) = 1 GEÇİŞİ  (θ ≡ KALİB_orta/g_cal; 172 §G2.0'a göre")
    print("      θ = 1 ⟺ beş ARTIK-ALAN teriminin toplamı SIFIR)")
    print("  λ:    %s" % "  ".join("%.2f" % x for x in lam))
    print("  θ:    %s" % "  ".join("%.4f" % x for x in th))
    kt = [float(lam[i] + (lam[i + 1] - lam[i]) * (th[i] - 1)
                / (th[i] - th[i + 1]))
          for i in range(len(lam) - 1) if (th[i] - 1) * (th[i + 1] - 1) < 0]
    print("  ** θ = 1 geçişi: λ ≈ %s   (θ(0.40) = %.4f > 1 İLK KEZ)"
          % (["%.4f" % k for k in kt], th[0]))
    out["T4"] = dict(theta=[float(x) for x in th], kesis=kt)

    # ---- T5: c ve 4/π² -------------------------------------------------
    print("\nT5 — c(λ) ve 4/π² = %.7f" % C_HIP)
    print("  λ      c ± σ           c/(4/π²) − 1     z")
    for i, l in enumerate(lam):
        print("  %.2f   %.4f ± %.4f   %+7.2f%%        %+6.1fσ  %s"
              % (l, c[i], sc[i], 100 * (c[i] / C_HIP - 1),
                 (c[i] - C_HIP) / sc[i], "" if sag[i] else "(SAĞLIKSIZ)"))
    out["T5"] = dict(c=[float(x) for x in c], sc=[float(x) for x in sc],
                     C_HIP=C_HIP)

    os.makedirs(SCR173, exist_ok=True)
    json.dump(out, open(SCR173 + "/EK.json", "w"), indent=1, default=float)
    print("\n-> %s/EK.json" % SCR173)


if __name__ == "__main__":
    main()
