"""
173d — HAKEM YÜZLEŞMESİ: 172/G4'ün İKİ YOLU ↔ ÖLÇÜM
====================================================
Yeni ölçüm YOK. `169_k1` AYNEN import edilir (üye tanımı c_WX =
KALİB_u2/W_X — 170/171/172 boyunca aynı); M tanımı 171 §T2b ile birebir
(`M = geo.ort_b [KALİB_b(g)/KALİB_b(Hkeskin)]`, beş bant).

YÜZLEŞİLEN MÜHÜRLER — HİÇBİRİ DEĞİŞTİRİLMEDİ:
  `172/G4.json` (4 Eylül 13:32, gazlar inşa edilmeden):
     λ=0.40  çarpan: M 1.29469 c 0.37905 | doğrudan: M 1.38809 c 0.40516
     λ=1.45  çarpan: M 0.97435 c 0.47471 | doğrudan: M 1.00599 c 0.49098
  `173/ONKAYIT_L040.json`, `173/ONKAYIT_L145.json` (ölçümden önce,
     zaman damgalı): 171'in bütün M ailesi + M9/M9′ + M10.

ÖN-MÜHÜR (bu betiği koşmadan önce, 4 Eylül 2026):
 D1 172 kendi tercihini yazmıştı: "doğrudan yol daha güvenilirdir".
    Bu betik o tercihi SINAR; ıskalarsa da yazılır.
 D2 Bir ölçüm hatası kestirimi: c'nin bant-yayılımı λ ile büyüyor
    (L050 ±0.0072, L130 ±0.0136) ⇒ uçlarda ±0.008 (0.40) ve ±0.018
    (1.45) beklenir; M'de ±0.008 ve ±0.020. Yani λ=0.40'taki %6.7'lik
    ayrışma **~8σ**, λ=1.45'teki %3.2'lik ayrışma **~1.6σ** demektir:
    **HAKEM ESAS OLARAK λ = 0.40'TIR.**
 D3 Çarpan yolunun hatası İKİ parçaya ayrışır: (i) GİRDİ hatası
    (ekstrapole edilen g_E, g_X, ρ_X ile ölçülenler arasındaki fark) ve
    (ii) MODEL hatası (θ = θ₀(ρ_X/ρ_X₀)^0.651 taşıyıcısının kendi
    borcu). 173b'ye göre girdi hatası ≤ %1.5 bekleniyordu.

Kullanım: 173d_hakem.py
Çıktı:    scratchpad/173/HAKEM.json

SONUÇ (yalnız gerçek koşudan; scratchpad/173/log_173d.txt):
 λ = 0.40 (SAĞLIKLI — ASIL HAKEM):
   c = 0.4098 ± 0.0049   M = 1.4086 ± 0.0232   θ = 1.02923  α_g = 0.6474
   **ÇARPAN yolu  M 1.2947 → +8.80% (+4.9σ)  ⇒ ÖLDÜ**
   **DOĞRUDAN yol M 1.3881 → +1.47% (+0.9σ)  ⇒ AYAKTA**
   ⇒ D1 ✓✓✓ 172'nin kendi yazılı tercihi ("doğrudan yol daha
     güvenilirdir") ÖRNEKLEM-DIŞI DOĞRULANDI.
   171'in M ailesi: M0 +17.6σ, M2 +9.3σ, M3 +10.4σ, M8 +7.6σ,
   M1 +11.2σ, M7 +7.2σ ⇒ hepsi ÖLÜ; **M9 (171, 5 gaz) −0.1σ ile TAM
   İSABET**; M9′ +0.9σ.  4/π² c dilinde +1.12% (+0.9σ).
   **D3 ✓✓ HATA AYRIŞTIRMASI:** toplam +8.80% = GİRDİ +1.29% +
   MODEL +7.41% (girdiler g_E +0.35%, g_X +0.47%; θ modeli +8.88%).
   **P5 ÖLDÜ:** θ taşıyıcısı öngörü 0.94506, ölçüm 1.02923 (−%8.18).
   D2 ISKA: hata kestirimim ±0.008 idi, ±0.023 çıktı (hüküm değişmedi).
 λ = 1.40 (SAĞLIKLI; G4 mührü YOK, uydurulmadı):
   c = 0.4338 ± 0.0240   M = 0.8634 ± 0.0310   θ = 0.73954  α_g = 0.1590
   M9 1.0324 → −16.4% (−5.5σ);  M9′ 1.0013 → −13.8% (−4.4σ);
   8-nokta parabolü 1.0046 → −14.1% (−4.6σ).
   ⇒ log-λ parabolü SAĞLIKLI bir gazda, sadakat sınırından ÖNCE ölüyor.
 λ = 1.45 (SAĞLIKSIZ — HÜKÜM YOK, yalnız teşhis):
   R_bant 0.92-0.96 < 0.98; filtre gevşetilmedi.
   Teşhis: c ≈ 0.4098, M ≈ 0.7663; her aday −16%…−26% uzakta;
   G4'ün iki yolu da (−21.4% / −23.8%) ıskalıyor. **P5 orada da öldü**
   (+%27.9).
"""
import importlib
import json
import os
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("169_configs", "167_configs", "166_configs", "165_configs",
           "163_configs", "160_configs", "159_configs", "155_configs",
           "154_configs"):
    sys.path.insert(0, str(QM / _p))
K1 = importlib.import_module("169_k1")

SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
       "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
SCR173 = SCR + "/173"
C_HIP = 4.0 / np.pi ** 2
YENI = {"L040": 0.40, "L140": 1.40, "L145": 1.45}


LOS = [0.52, 0.56, 0.60, 0.64, 0.68]


def pencere(d):
    """Hüküm penceresinin BEŞ bandı (lo ile eşleşerek) + SAĞLIK bayrağı.

    169'un filtresi (R_bant ≥ 0.98, SNR ≥ 3, τ_eff < 0.85, Ms2 > 0)
    DEĞİŞTİRİLMEZ; yalnız düşen gazda bantlar TEŞHİS için de döndürülür
    ve `saglik` False olur. Sağlıksız gazda HÜKÜM VERİLMEZ.
    """
    ban, sag = [], []
    for lo in LOS:
        m = [x for x in d["bant"] if x.get("olculdu")
             and abs(x["lo"] - lo) < 1e-9]
        if not m:
            continue
        b = m[0]
        ban.append(b)
        sag.append(bool(b["R_bant"] >= 0.98 and b["SNR"] >= 3.0
                        and b["tau_eff"] < 0.85 and b["Ms2"] > 0))
    return ban, sag


def gaz_c(D, g):
    """171h.gaz_c ile birebir aynı — c, M girdileri ve hataları."""
    B, sag = pencere(D[g])
    rows = [K1.band_jk(b["cizgi"]) for b in B]
    c = np.array([r["c_WX"] for r in rows])
    s = np.array([r["sc_WX"] for r in rows])
    lg = np.log(c)
    cg = float(np.exp(lg.mean()))
    return dict(c=cg,
                s_jk=float(np.sqrt(np.sum((s / c) ** 2)) / len(c)) * cg,
                s_bant=float(np.std(lg, ddof=1) / np.sqrt(len(lg))) * cg,
                s_tot=float(np.hypot(
                    float(np.sqrt(np.sum((s / c) ** 2)) / len(c)) * cg,
                    float(np.std(lg, ddof=1) / np.sqrt(len(lg))) * cg)),
                c_bant=[float(x) for x in c],
                K=np.array([r["KALIB_u2"] for r in rows]),
                sK=np.array([r["sKALIB_u2"] for r in rows]),
                W=np.array([r["W_X"] for r in rows]),
                t=np.array([r["tau_eff"] for r in rows]),
                R=np.array([b["R_bant"] for b in B]),
                saglik=sag, saglikli=bool(all(sag)))


def M_of(G, H):
    Mo = float(np.exp(np.mean(np.log(G["K"] / H["K"]))))
    sj = np.sqrt(np.sum((G["sK"] / G["K"]) ** 2
                        + (H["sK"] / H["K"]) ** 2)) / len(G["K"]) * Mo
    ss = float(np.std(np.log(G["K"] / H["K"]), ddof=1)
               / np.sqrt(len(G["K"]))) * Mo
    return Mo, float(np.hypot(sj, ss))


def main():
    os.makedirs(SCR173, exist_ok=True)
    D = K1.yukle()
    G4 = json.load(open(SCR + "/172/G4.json"))
    G1 = json.load(open(SCR + "/172/G1.json"))
    G2 = json.load(open(SCR + "/172/G2.json"))
    H = gaz_c(D, "Hkeskin")
    out = {}
    varsa = [g for g in ("L040", "L140", "L145") if g in D]
    print("=" * 104)
    print("173d — HAKEM YÜZLEŞMESİ   (ölçülen yeni gazlar: %s)"
          % ", ".join(varsa))
    print("=" * 104)

    for ad in varsa:
        lam = YENI[ad]
        i = 0 if lam < 1 else 1
        ON = json.load(open(SCR173 + "/ONKAYIT_%s.json" % ad))
        G = gaz_c(D, ad)
        Mo, sMo = M_of(G, H)
        a_ = D[ad]["artik"]
        gcal = a_["gE"] * a_["gX"] ** 2
        orta = [b for b in D[ad]["bant"]
                if b.get("olculdu") and abs(b["lo"] - 0.60) < 1e-9][0]
        th = orta["KALIB_u2"] / gcal          # 172'nin θ tanımı (orta bant)
        th5 = float(np.exp(np.mean(np.log(G["K"])))) / gcal
        alfa = float(-np.polyfit(G["t"] ** 2, np.log(G["K"]), 1)[0])

        print("\n" + "-" * 104)
        print("### %s  (λ = %.2f)   ön kayıt %s" % (ad, lam, ON["zaman"]))
        print("-" * 104)
        if not G["saglikli"]:
            print("  !!! 169'un SAĞLIK FİLTRESİ DÜŞTÜ (R_bant ≥ 0.98): "
                  "bantlar %s" % ["%.4f" % r for r in G["R"]])
            print("  !!! Filtre GEVŞETİLMEDİ. Aşağıdaki sayılar HÜKÜM DEĞİL, "
                  "TEŞHİStir; σ'lar da anlamsızdır.")
        print("  bant τ_eff  KALİB_u2±jk        W_X       c_WX      R_bant")
        for j in range(len(G["t"])):
            print("   %d  %.4f  %.5f±%.5f  %.5f   %.5f   %.4f %s"
                  % (j + 1, G["t"][j], G["K"][j], G["sK"][j], G["W"][j],
                     G["c_bant"][j], G["R"][j],
                     "" if G["saglik"][j] else "  ← SAĞLIKSIZ"))
        print("  **c(%s) = %.4f ± %.4f(jk) ± %.4f(bant) = ±%.4f**   "
              "4/π² farkı %+.2f%%"
              % (ad, G["c"], G["s_jk"], G["s_bant"], G["s_tot"],
                 100 * (G["c"] / C_HIP - 1)))
        print("  **M(%s) = %.4f ± %.4f**   (bant bant: %s)"
              % (ad, Mo, sMo, "  ".join("%.4f" % x
                                        for x in G["K"] / H["K"])))
        print("  θ(orta bant) = %.5f   θ(5 bant geo) = %.5f   "
              "g_cal = %.5f   α_g = %.4f" % (th, th5, gcal, alfa))

        # ---------- HAKEM TABLOSU ----------
        hasG4 = ad in ("L040", "L145")     # G4 yalnız bu ikisini mühürledi
        hk, kaz, ikisi_de = {}, None, False
        if not hasG4:
            print("\n  ** λ = %.2f için 172/G4 MÜHRÜ YOKTUR (G4 yalnız 0.40 "
                  "ve 1.45'i mühürledi) — uydurulmadı. **" % lam)
        else:
            print("\n  ** HAKEM — 172/G4'ün İKİ YOLU (mühür 4 Eyl 13:32) **%s"
                  % ("" if G["saglikli"] else
                     "   [HÜKÜM YOK — sağlık filtresi düştü; yalnız teşhis]"))
            print("  %-26s %9s %9s %8s | %9s %9s %8s"
                  % ("yol", "M(mühür)", "M ölç−ön", "σ_M", "c(mühür)",
                     "c ölç−ön", "σ_c"))
            for nm, mk, ck in (("ÇARPAN yolu", G4["ongoru"]["M"][i],
                                G4["ongoru"]["c"][i]),
                               ("DOĞRUDAN yol (M9′)",
                                G4["ongoru"]["M_dogrudan"][i],
                                G4["ongoru"]["c_dogrudan"][i])):
                dm, zm = 100 * (Mo / mk - 1), (Mo - mk) / sMo
                dc, zc = 100 * (G["c"] / ck - 1), (G["c"] - ck) / G["s_tot"]
                hk[nm] = dict(M=mk, c=ck, dM=dm, zM=float(zm), dc=dc,
                              zc=float(zc))
                print("  %-26s %9.4f %+8.2f%% %+7.1fσ | %9.4f %+8.2f%% "
                      "%+7.1fσ" % (nm, mk, dm, zm, ck, dc, zc))
            kaz = min(hk, key=lambda k: abs(hk[k]["zM"]))
            print("  ⇒ M'de daha yakın olan: **%s**  (|z| = %.1f vs %.1f)"
                  % (kaz, abs(hk[kaz]["zM"]),
                     max(abs(hk[k]["zM"]) for k in hk)))
            ikisi_de = all(abs(hk[k]["zM"]) > 3 for k in hk)
            print("  ⇒ İKİSİ DE 3σ'nın DIŞINDA mı? **%s**"
                  % ("EVET — iki yol da ıskaladı" if ikisi_de else "hayır"))

        # ---------- 171'in M ailesi ----------
        print("\n  ** 171'in M AİLESİ (parametreler 171e'de donduruldu) **")
        print("  %-40s %8s %8s | %9s %8s | %9s %8s"
              % ("aday", "M(ön)", "c(ön)", "M ölç−ön", "σ_M", "c ölç−ön",
                 "σ_c"))
        aile = {}
        for nm, r in ON["ongoru"].items():
            dm, zm = 100 * (Mo / r["M"] - 1), (Mo - r["M"]) / sMo
            dc, zc = 100 * (G["c"] / r["c"] - 1), (G["c"] - r["c"]) / G["s_tot"]
            aile[nm] = dict(M=r["M"], c=r["c"], dM=dm, zM=float(zm), dc=dc,
                            zc=float(zc))
            print("  %-40s %8.4f %8.4f | %+8.2f%% %+7.1fσ | %+8.2f%% %+7.1fσ"
                  % (nm, r["M"], r["c"], dm, zm, dc, zc))
        print("  %-40s %8s %8.4f | %9s %8s | %+8.2f%% %+7.1fσ"
              % ("(çapa) 4/π²", "—", C_HIP, "", "",
                 100 * (G["c"] / C_HIP - 1), (G["c"] - C_HIP) / G["s_tot"]))

        # ---------- D3: girdi hatası vs model hatası ----------
        print("\n  ** D3 — ÇARPAN YOLUNUN HATA AYRIŞTIRMASI **")
        gE0 = G1["Hkeskin"]["gE"]
        gX0 = G1["Hkeskin"]["gX"]
        th0 = G2["theta"]["Hkeskin"]
        thp = ON["theta_tasiyici"]                  # ÖLÇÜLEN ρ_X ile
        M_girdi = ((a_["gE"] / gE0) * (a_["gX"] / gX0) ** 2 * (thp / th0))
        if hasG4:
            print("  %-46s %10s %10s %9s"
                  % ("", "MÜHÜR", "ÖLÇÜLEN", "fark%"))
            for nm, mu, ol in (("g_E (girdi)", G4["ongoru"]["gE"][i],
                                a_["gE"]),
                               ("g_X (girdi)", G4["ongoru"]["gX"][i],
                                a_["gX"]),
                               ("θ  (model: (ρ_X/ρ_X₀)^0.651)",
                                G4["ongoru"]["theta"][i], th)):
                print("  %-46s %10.5f %10.5f %+8.2f%%"
                      % (nm, mu, ol, 100 * (ol / mu - 1)))
            print("  %-46s %10.5f %10.5f %+8.2f%%"
                  % ("M: mühür → ölçülen girdilerle (θ hâlâ modelden)",
                     G4["ongoru"]["M"][i], M_girdi,
                     100 * (M_girdi / G4["ongoru"]["M"][i] - 1)))
            print("  %-46s %10.5f %10.5f %+8.2f%%"
                  % ("M: ölçülen girdiler+model → GERÇEK M", M_girdi, Mo,
                     100 * (Mo / M_girdi - 1)))
            print("  ⇒ toplam %+.2f%% = GİRDİ %+.2f%% + MODEL %+.2f%%"
                  % (100 * (Mo / G4["ongoru"]["M"][i] - 1),
                     100 * (M_girdi / G4["ongoru"]["M"][i] - 1),
                     100 * (Mo / M_girdi - 1)))
        # θ taşıyıcısının kendi artığı (P5)
        print("  ** P5 (bonus) — θ taşıyıcısı: öngörü %.5f, ölçüm %.5f, "
              "artık %+.2f%%  ⇒ %s"
              % (thp, th, 100 * (thp / th - 1),
                 "AYAKTA (≤%3)" if abs(100 * (thp / th - 1)) <= 3
                 else "ÖLDÜ (>%3)"))

        out[ad] = dict(lam=lam, zaman_onkayit=ON["zaman"], c=G["c"],
                       s_tot=G["s_tot"], s_jk=G["s_jk"], s_bant=G["s_bant"],
                       c_bant=G["c_bant"], M=Mo, sM=sMo, theta=float(th),
                       theta5=float(th5), gcal=float(gcal), alfa=alfa,
                       gE=a_["gE"], gX=a_["gX"],
                       tau=[float(x) for x in G["t"]],
                       KALIB=[float(x) for x in G["K"]],
                       W_X=[float(x) for x in G["W"]],
                       M_bant=[float(x) for x in (G["K"] / H["K"])],
                       hakem=hk, aile=aile, kazanan=kaz,
                       saglikli=G["saglikli"], R_bant=[float(x) for x in G["R"]],
                       ikisi_de_iska=bool(ikisi_de),
                       M_girdi_duzeltmeli=float(M_girdi),
                       theta_tasiyici=float(thp))

    json.dump(out, open(SCR173 + "/HAKEM.json", "w"), indent=1, default=float)
    print("\n-> %s/HAKEM.json" % SCR173)


if __name__ == "__main__":
    main()
