# -*- coding: utf-8 -*-
"""
174c — K2: DONDURULMUŞ KURALIN UYGULANMASI ve ÖLÇÜMLE YÜZLEŞME
================================================================
Girdi: scratchpad/174/ONKAYIT_K2.json (174a, zaman damgalı, sha256'lı),
       scratchpad/174/K1_<gaz>.json (174b'nin ölçtükleri),
       scratchpad/172/G1.json + G2.json + G3.json (172'nin defteri),
       scratchpad/167/C_<gaz>.json (167'nin artıkları).
YENİ ÖLÇÜM YOK — bu betik yalnız dondurulmuş kuralı uygular.

DENETİM D1: 172e'nin `lam_es` ters çevirmesi BİREBİR kopyalandı (172e üst
düzey kod koştuğu için import edilemiyor); doğruluğu, 172e'nin yayımlanmış
λ_eş tablosunu G1/G2.json'dan yeniden üreterek sınanır (fark ≤ 1e−9
olmalı).
"""
import json
import sys
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S174, S172, S167 = SCR / "174", SCR / "172", SCR / "167"

LAM7 = ["L050", "L060", "L070", "L085", "Hkeskin", "L115", "L130"]
LV = np.array([.50, .60, .70, .85, 1.00, 1.15, 1.30])
G1 = json.load(open(S172 / "G1.json"))
G2 = json.load(open(S172 / "G2.json"))
G3 = json.load(open(S172 / "G3.json"))
ONK = json.load(open(S174 / "ONKAYIT_K2.json"))
H = ONK["hedef"]


def lam_es(y, v):                      # 172e_gercek_adres.py'den BİREBİR
    """y(λ) eğrisinin v'yi verdiği λ (tekdüze kolda; kol otomatik)."""
    for kol in (slice(0, 6), slice(0, 5), slice(None)):
        yy, ll = y[kol], LV[kol]
        d = np.diff(yy)
        if np.all(d > 0):
            if yy[0] <= v <= yy[-1]:
                return float(np.interp(v, yy, ll)), True
        elif np.all(d < 0):
            if yy[-1] <= v <= yy[0]:
                return float(np.interp(v, yy[::-1], ll[::-1])), True
    return float("nan"), False


def egri(f, gazlar=LAM7):
    return np.array([f(g) for g in gazlar])


# ═════════════════════════════════════════════════════════════════════
print("=" * 74)
print("174c — K2: ÖN-KAYIT %s  (sha256 %s)" % (ONK["zaman"],
                                               ONK["sha256"][:16]))
print("=" * 74)

# --- D1: ters çevirici denetimi --------------------------------------
DEN = {"rE": lambda g: G1[g]["E"]["r"], "Q_E": lambda g: G1[g]["E"]["Q"],
       "ρ_E": lambda g: G1[g]["E"]["rho"], "μ̂²_E": lambda g: G1[g]["E"]["mu2"],
       "g_E": lambda g: G1[g]["gE"], "rX": lambda g: G1[g]["X"]["r"],
       "Q_X": lambda g: G1[g]["X"]["Q"], "ρ_X": lambda g: G1[g]["X"]["rho"],
       "g_X": lambda g: G1[g]["gX"], "θ": lambda g: G2["theta"][g],
       "σ_X̃": lambda g: G1[g]["sigX"], "σ_ds": lambda g: G1[g]["sigds"],
       "σ_Ĉ": lambda g: G1[g]["sigC"]}
print("\n[D1] ters çevirici denetimi (172e'nin G3.json'u ile):")
worst = 0.0
for ad, f in DEN.items():
    le, _ = lam_es(egri(f), f("son"))
    ref = G3["lam_es"].get(ad, float("nan"))
    if ref == ref:
        worst = max(worst, abs(le - ref))
print("     maks |λ_eş(174c) − λ_eş(172e)| = %.2e  %s"
      % (worst, "✓" if worst < 1e-9 else "✗ KIRIK"))

# --- K1 ölçümlerini yükle --------------------------------------------
K1 = {}
for p in sorted(S174.glob("K1_*.json")):
    d = json.load(open(p))
    K1[d["veri"]] = d
mev = [g for g in LAM7 if g in K1]
print("\n[K1] ölçülmüş gazlar: %s   (λ merdiveni: %s)"
      % (", ".join(sorted(K1)), ", ".join(mev)))
if "son" not in K1:
    sys.exit("son ölçülmemiş — 174b önce koşmalı")

# ═════════════════════════════════════════════════════════════════════
# K2-F — KÖPRÜ KİMLİĞİ: 162'nin R'si = 172'nin π'si mi?
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 74)
print("K2-F — KÖPRÜ KİMLİĞİ (162 makinesi τ≤0.86, site m_n  ↔  172 defteri"
      " τ≤0.95, site s_n)")
print("=" * 74)
print("%-10s %10s %10s %9s | %10s %10s %9s | %9s"
      % ("gaz", "R_η(162)", "π_E(172)", "fark%", "R_X(174)", "π_X(172)",
         "fark%", "R_Ĉ(162)"))
KF = {}
for g in sorted(K1, key=lambda z: (z != "son", z)):
    d = K1[g]
    piE = G1[g]["E"]["pi"]          # = K/V_O  (172b'nin defteri, τ≤0.95)
    piX = G1[g]["X"]["pi"]
    fE = 100 * (d["eta"]["R"] / piE - 1)
    fX = 100 * (d["Xtil"]["R"] / piX - 1)
    KF[g] = dict(R_eta=d["eta"]["R"], piE=piE, fE=fE, R_X=d["Xtil"]["R"],
                 piX=piX, fX=fX, R_C=d["Chat"]["R"])
    print("%-10s %10.5f %10.5f %+8.2f%% | %10.5f %10.5f %+8.2f%% | %9.5f"
          % (g, d["eta"]["R"], piE, fE, d["Xtil"]["R"], piX, fX,
             d["Chat"]["R"]))
mx = max(abs(v["fE"]) for v in KF.values())
print("  ön-kayıt K2-F: |R_η/π_E − 1| ≤ %%%.0f ⇒ ölçülen maks %.2f%%  → %s"
      % (ONK["kural"]["K2-F"]["esik_pct"], mx,
         "✓ GEÇTİ" if mx <= ONK["kural"]["K2-F"]["esik_pct"] else "✗ ÖLDÜ"))
print("  (162'nin gerçek gazda yayımladığı R_η = %.3f ; 174b'nin yeniden "
      "ölçümü = %.4f)" % (H["R_eta_162_gercek"], K1["son"]["eta"]["R"]))

# ═════════════════════════════════════════════════════════════════════
# K2-A / K2-G — λ_eş ÖNGÖRÜLERİ
# ═════════════════════════════════════════════════════════════════════
def iki_nokta(y85, y100, v):
    """log-doğrusal iki nokta (λ=0.85, 1.00) tersi."""
    b = (np.log(y100) - np.log(y85)) / 0.15
    return 0.85 + (np.log(v) - np.log(y85)) / b if b else float("nan")


print("\n" + "=" * 74)
print("K2-A / K2-G — λ_eş ÖNGÖRÜLERİ (dondurulmuş kural)")
print("=" * 74)
SON = {"R_η": K1["son"]["eta"]["R"], "R_Ĉ": K1["son"]["Chat"]["R"],
       "R_X": K1["son"]["Xtil"]["R"],
       "m3": K1["son"]["m3"]["olc"], "m3_çizgi": K1["son"]["m3"]["cizgi"],
       "skew_η": K1["son"]["eta"]["skew"],
       "skew_η_çizgi": K1["son"]["eta"]["skew_cizgi"],
       "P_η": K1["son"]["eta"]["P"], "Var_η": K1["son"]["eta"]["Var"]}
ALICI = {"R_η": lambda d: d["eta"]["R"], "R_Ĉ": lambda d: d["Chat"]["R"],
         "R_X": lambda d: d["Xtil"]["R"], "m3": lambda d: d["m3"]["olc"],
         "m3_çizgi": lambda d: d["m3"]["cizgi"],
         "skew_η": lambda d: d["eta"]["skew"],
         "skew_η_çizgi": lambda d: d["eta"]["skew_cizgi"],
         "P_η": lambda d: d["eta"]["P"], "Var_η": lambda d: d["eta"]["Var"]}
LE = {}
print("  λ merdiveni: " + " ".join("%s(%.2f)" % (g, l)
                                   for g, l in zip(LAM7, LV)))
print("%-14s %10s | %-58s | %8s %8s %s"
      % ("nicelik", "son", "λ merdiveni değerleri", "λ_eş(7n)", "λ_eş(2n)",
         "durum"))
for ad, f in ALICI.items():
    v = SON[ad]
    l2 = (iki_nokta(f(K1["L085"]), f(K1["Hkeskin"]), v)
          if ("L085" in K1 and "Hkeskin" in K1 and v > 0
              and f(K1["L085"]) > 0 and f(K1["Hkeskin"]) > 0) else float("nan"))
    l7, y = float("nan"), None
    dur = "merdiven eksik"
    if len(mev) == 7:
        y = egri(lambda g: f(K1[g]))
        l7, _ok = lam_es(y, v)
        d = np.diff(y)
        tek = "tekdüze↑" if np.all(d[:5] > 0) else (
            "tekdüze↓" if np.all(d[:5] < 0) else "TÜMSEK/ÇUKUR")
        ic = y.min() <= v <= y.max()
        dur = "%s, son %s" % (tek, "içeride" if ic else
                              "**ARALIK DIŞI** (%s max %.5f, min %.5f)"
                              % ("üstünde" if v > y.max() else "altında",
                                 y.max(), y.min()))
    LE[ad] = dict(son=v, lam7=l7, lam2=l2, durum=dur,
                  merdiven=[float(x) for x in y] if y is not None else None,
                  L085=f(K1["L085"]) if "L085" in K1 else None,
                  Hk=f(K1["Hkeskin"]) if "Hkeskin" in K1 else None)
    print("%-14s %10.5f | %-58s | %8s %8s %s"
          % (ad, v,
             " ".join("%.4f" % x for x in y) if y is not None else "—",
             ("%.4f" % l7) if l7 == l7 else "  —  ",
             ("%.4f" % l2) if l2 == l2 else "  —  ", dur))

def hukum(deger, hedef, tam, kismi):
    d = abs(deger - hedef)
    return ("TAM İSABET" if d <= tam else
            "KISMİ" if d <= kismi else "ÖLDÜ"), d


lamA = LE["R_η"]["lam7"] if LE["R_η"]["lam7"] == LE["R_η"]["lam7"] \
    else LE["R_η"]["lam2"]
hA, dA = hukum(lamA, H["lam_es_E"], 0.05, 0.10)
if lamA != lamA:
    print("\n  K2-A: λ_eş^ön(E) **YOK** — R_η(son) λ ailesinin ARALIĞI "
          "DIŞINDA.  %s" % LE["R_η"]["durum"])
    print("        → ön-kayıtlı ölçüt (|λ_eş − 0.7768| ≤ 0.05) "
          "SAĞLANAMAZ ⇒ **K2-A ÖLDÜ** (kurtarma yok)")
else:
    print("\n  K2-A: λ_eş^ön(E) = %.4f   HEDEF λ_eş(E kanalı) = %.4f  "
          "(fark %+.4f) → **%s**"
          % (lamA, H["lam_es_E"], lamA - H["lam_es_E"], hA))
    print("        (ikinci hedef λ_eş(g_E) = %.4f, fark %+.4f)"
          % (H["lam_es_gE"], lamA - H["lam_es_gE"]))
# defterin kendi π_E'si de aynı sınavdan geçirilir (özdeşlik denetimi)
yp = egri(lambda g: G1[g]["E"]["pi"])
lp, _ = lam_es(yp, G1["son"]["E"]["pi"])
print("        [defter denetimi] π_E merdiveni %s ; π_E(son)=%.5f ⇒ "
      "λ_eş(π_E) = %s"
      % (" ".join("%.4f" % x for x in yp), G1["son"]["E"]["pi"],
         ("%.4f" % lp) if lp == lp else "**YOK (aralık dışı)**"))
lamG = LE["R_X"]["lam7"] if LE["R_X"]["lam7"] == LE["R_X"]["lam7"] \
    else LE["R_X"]["lam2"]
hG, dG = hukum(lamG, H["lam_es_X"], 0.05, 0.10)
print("  K2-G: λ_eş^ön(X) = %s   HEDEF λ_eş(X) = %.4f  (fark %+.4f) "
      "→ **%s**  [%s]"
      % (("%.4f" % lamG) if lamG == lamG else "YOK", H["lam_es_X"],
         lamG - H["lam_es_X"], hG, LE["R_X"]["durum"]))

# ═════════════════════════════════════════════════════════════════════
# K2-B — g_E PAYI (parametresiz)
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 74)
print("K2-B — g_E PAYI: Δlog g_E^ön = log R_η(son) − log R_η(λ_çapa)")
print("=" * 74)
capa = H["lam_capa"]
if len(mev) == 7:
    Rlam = egri(lambda g: K1[g]["eta"]["R"])
    Rcapa = float(np.interp(capa, LV, Rlam))
    rhoLam = egri(lambda g: G1[g]["E"]["rho"])
    rhocapa = float(np.interp(capa, LV, rhoLam))
else:
    b = (np.log(K1["Hkeskin"]["eta"]["R"]) - np.log(K1["L085"]["eta"]["R"])) / .15
    Rcapa = float(np.exp(np.log(K1["L085"]["eta"]["R"]) + b * (capa - 0.85)))
    rhocapa = float(np.interp(capa, LV, egri(lambda g: G1[g]["E"]["rho"])))
dlog = float(np.log(K1["son"]["eta"]["R"] / Rcapa))
dlog_rho = float(np.log(G1["son"]["E"]["rho"] / rhocapa))
pB = 100 * (np.exp(dlog) - 1)
pB2 = 100 * (np.exp(dlog - dlog_rho) - 1)
hB, dB = hukum(pB, H["d_gE_pct"], 0.5, 1.0)
print("  R_η(son) = %.5f   R_η(λ=%.4f) = %.5f" % (K1["son"]["eta"]["R"],
                                                  capa, Rcapa))
print("  Δlog g_E^ön (ρ düzeltmesiz) = %+.5f  ⇒ **%+.2f%%**   HEDEF %+.2f%%"
      "  (fark %+.2f puan) → **%s**"
      % (dlog, pB, H["d_gE_pct"], pB - H["d_gE_pct"], hB))
print("  [yan bilgi, hükme girmez] ρ_E düzeltmeli = %+.2f%%  "
      "(Δlog ρ_E = %+.5f)" % (pB2, dlog_rho))
# --- defterin kendi ayrıştırması: Δlog g_E ≡ Δlog π_E − Δlog ρ_E ------
piC = float(np.interp(capa, LV, egri(lambda g: G1[g]["E"]["pi"])))
gC = float(np.interp(capa, LV, egri(lambda g: G1[g]["gE"])))
dpi = float(np.log(G1["son"]["E"]["pi"] / piC))
dg = float(np.log(G1["son"]["gE"] / gC))
print("  [ÖZDEŞ AYRIŞTIRMA, 172 defteri] Δlog g_E = Δlog π_E − Δlog ρ_E :"
      " %+.5f = %+.5f − (%+.5f)  (kalıntı %.1e)"
      % (dg, dpi, dlog_rho, abs(dg - dpi + dlog_rho)))
print("      ⇒ g_E artığının **%.0f%%'i GİRİŞİM (π)**, **%.0f%%'i MODEL "
      "GÜCÜ (ρ)**" % (100 * dpi / dg, -100 * dlog_rho / dg))

# ═════════════════════════════════════════════════════════════════════
# K2-C — θ PAYI
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 74)
print("K2-C — θ PAYI: Ĉ kanalı girişimsizse öngörü SIFIR")
print("=" * 74)
RC = {g: K1[g]["Chat"]["R"] for g in K1}
yay = 100 * (max(RC.values()) / min(RC.values()) - 1)
print("  R_Ĉ: " + "  ".join("%s=%.4f" % (g, RC[g]) for g in sorted(RC)))
print("  yayılım %.2f%%  (ön-kayıt: ±%%2 içinde ⇒ %s)"
      % (yay, "✓ tuttu" if yay <= 4.0 else "✗ ön-kayıt YANILDI"))
print("  Δθ^ön = %+.2f%%   HEDEF %+.2f%%  → **%s**"
      % (0.0, H["d_theta_pct"], "ÖLDÜ (öngörülemedi)"))

# ═════════════════════════════════════════════════════════════════════
# K2-E — ÜÇÜNCÜ MOMENT KÖPRÜSÜ
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 74)
print("K2-E — ÜÇÜNCÜ MOMENT: λ_eş(m3) < λ_eş(R_η) mi?")
print("=" * 74)
lamE = LE["m3"]["lam7"] if LE["m3"]["lam7"] == LE["m3"]["lam7"] \
    else LE["m3"]["lam2"]
print("  m3(son) = %+.5f   λ_eş(m3) = %s   λ_eş(R_η) = %.4f"
      % (SON["m3"], ("%.4f" % lamE) if lamE == lamE else "—", lamA))
print("  tek yönlü ön-kayıt (λ_eş(m3) < λ_eş(R_η)) → **%s**"
      % ("✓ TUTTU" if (lamE == lamE and lamE < lamA) else
         "✗ ÖLDÜ" if lamE == lamE else "— ölçülemedi"))
if lamE == lamE:
    print("  nicel hedef λ_eş(θ) = %.4f  (fark %+.4f) → %s"
          % (H["lam_es_theta"], lamE - H["lam_es_theta"],
             "TAM" if abs(lamE - H["lam_es_theta"]) <= 0.10 else "ıska"))

# ═════════════════════════════════════════════════════════════════════
# K2-D — ΔM KAPANIŞI / H-G3
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 74)
print("K2-D — ΔM KAPANIŞI ve H-G3 MÜHRÜ")
print("=" * 74)
dM = np.log(1 + H["dM_pct"] / 100)
kap = 100 * dlog / dM
print("  Δlog M(ölçülen) = %+.5f (+%.2f%%)" % (dM, H["dM_pct"]))
print("  köprünün kapattığı (yalnız g_E payı) = %+.5f  ⇒ **%.1f%%**"
      % (dlog, kap))
print("  ön-kayıtlı beklenti %.0f%%, mühür eşiği %.0f%% → **H-G3 %s**"
      % (ONK["kural"]["K2-D"]["beklenen_pct"],
         ONK["kural"]["K2-D"]["muhur_esigi_pct"],
         "MÜHÜRLENDİ" if kap >= 70 else "MÜHÜRLENMEDİ (ön-kayıt böyle "
         "diyordu)"))

json.dump(dict(onkayit=ONK["zaman"], sha=ONK["sha256"], KF=KF,
               LE={k: v for k, v in LE.items()}, lamA=lamA, lamG=lamG,
               dlog_gE=dlog, pct_gE=pB, pct_gE_rho=pB2, kapanis_pct=kap,
               R_capa=Rcapa, capa=capa, yay_RC=yay,
               lam_m3=lamE if lamE == lamE else None),
          open(S174 / "K2.json", "w"), indent=1)
print("\n-> %s" % (S174 / "K2.json"))
