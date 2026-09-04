# -*- coding: utf-8 -*-
"""
176f — K2 YÜZLEŞME + K3 NİHAİ AYRIŞIM TABLOSU
==============================================
YENİ ÖLÇÜM YOKTUR. Yalnızca `176/ONKAYIT_K2.json`'daki dondurulmuş kural
(F0…F9) ölçülen sayılara uygulanır ve hükümler yazılır.

Girdi: 176/ONKAYIT_K2.json, 176/insa_VF*.json, 176/G_VF*.json,
       176/K1_VF*.json, 176/K3_VF*.json,
       172/G1.json, 174/K1_{son,Hkeskin,HA4}.json,
       174/K3_{son,Hkeskin,HA4}.json
Çıktı: 176/K2.json
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S176, S174, S172 = SCR / "176", SCR / "174", SCR / "172"
VEK = sys.argv[1:] or ["VF1", "VF2"]

ONK = json.load(open(S176 / "ONKAYIT_K2.json"))
KUR = ONK["kural"]
G1 = json.load(open(S172 / "G1.json"))
K1 = {g: json.load(open(S174 / f"K1_{g}.json"))
      for g in ("son", "Hkeskin", "HA4")}
K3 = {g: json.load(open(S174 / f"K3_{g}.json"))
      for g in ("son", "Hkeskin", "HA4")}
GV = {v: json.load(open(S176 / f"G_{v}.json")) for v in VEK}
IV = {v: json.load(open(S176 / f"insa_{v}.json")) for v in VEK}
for v in VEK:
    K1[v] = json.load(open(S176 / f"K1_{v}.json"))
    K3[v] = json.load(open(S176 / f"K3_{v}.json"))

Hk, sn, ha = G1["Hkeskin"], G1["son"], G1["HA4"]
OUT = dict(onkayit=ONK["zaman"], sha=ONK["sha256"], vekiller=VEK)

print("=" * 78)
print("176f — K2 YÜZLEŞME   [ön-kayıt %s   sha %s]"
      % (ONK["zaman"], ONK["sha256"][:16]))
print("=" * 78)


def ort_ve_sacilim(f):
    """İki tohumun ortalaması ve |y1−y2| (F7 hata çubuğu)."""
    y = [float(f(v)) for v in VEK]
    return float(np.mean(y)), (abs(y[0] - y[1]) if len(y) > 1 else 0.0), y


def hukum(ad, deger, sac, esik, yon, olum_metni, yasa_metni):
    """F7 kuralı: |y1−y2| > |y_ort − eşik| ise HÜKÜMSÜZ."""
    mesafe = abs(deger - esik)
    if sac > mesafe:
        return "HÜKÜMSÜZ (tohum gürültüsü)", False
    ok = (deger >= esik) if yon == ">=" else (deger <= esik)
    return (yasa_metni if ok else olum_metni), bool(ok)


# ═══════════════════════════════════════════════════════════════════════
# F0 — İNŞA + ÖLÇÜM KAPILARI
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 78)
print("F0 — İNŞA VE ÖLÇÜM KAPILARI")
print("=" * 78)
print("%-6s %10s %8s %10s %8s %8s %9s %9s"
      % ("gaz", "ilk-kök", "maks|F|", "sıralılık", "zarf", "φ≡0", "R_bant",
         "σ_ds"))
f0 = {}
for v in VEK:
    t = IV[v]
    rb = GV[v]["R_bant_min"]
    kap = t["kapilar"]
    f0[v] = dict(kapilar=kap, R_bant_min=rb, F0_R=GV[v]["F0_R_bant"],
                 maxF=t["maxF"], hucre=t["hucre_benzersiz"],
                 zarf=t["zarf_fark"], sigds=t["sigma_ds"],
                 hepsi=bool(all(kap.values()) and GV[v]["F0_R_bant"]))
    print("%-6s %6d/%d %8.2e %10s %8.1e %8.1e %9.4f %9.5f"
          % (v, t["hucre_benzersiz"], t["n"], t["maxF"],
             "TAM" if t["sirali"] else "BOZUK", t["zarf_fark"], t["g2_dS"],
             rb, t["sigma_ds"]))
print("  F0: %s" % ("✓ HEPSİ GEÇTİ" if all(f0[v]["hepsi"] for v in VEK)
                    else "✗ KAPI TUTMADI"))
OUT["F0"] = f0

# ═══════════════════════════════════════════════════════════════════════
# ÖLÇÜLEN DEFTER
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 78)
print("ÖLÇÜLEN DEFTER — vekiller ↔ çapalar")
print("=" * 78)
print("%-9s %9s %9s %9s %9s %9s %9s %9s"
      % ("gaz", "R_η", "R_Ĉ", "R_X", "g_E", "g_X", "θ", "M"))


def satir(g, R):
    return ("%-9s %9.5f %9.5f %9.5f %9.6f %9.6f %9.6f %9.6f"
            % (g, K1[g]["eta"]["R"], K1[g]["Chat"]["R"], K1[g]["Xtil"]["R"],
               R["gE"], R["gX"], R["th"], R["KAL"]))


for g, R in (("son", sn), ("Hkeskin", Hk), ("HA4", ha)):
    print(satir(g, R))
for v in VEK:
    print(satir(v, GV[v]))
OUT["defter"] = {g: dict(R_eta=K1[g]["eta"]["R"], R_C=K1[g]["Chat"]["R"],
                         R_X=K1[g]["Xtil"]["R"],
                         gE=R["gE"], gX=R["gX"], th=R["th"], KAL=R["KAL"],
                         pi_E=R["E"]["pi"], rho_E=R["E"]["rho"],
                         Q_E=R["E"]["Q"], mu2_E=R["E"]["mu2"],
                         r_E=R["E"]["r"], sigds=R["sigds"], sigC=R["sigC"])
                 for g, R in [("son", sn), ("Hkeskin", Hk), ("HA4", ha)]
                 + [(v, GV[v]) for v in VEK]}

# ═══════════════════════════════════════════════════════════════════════
# F1 — R_η(vekil)
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 78)
print("F1 — R_η(vekil):  merkez %.4f  bant [%.2f, %.2f]  ölüm ≥ 1.15"
      % (KUR["F1"]["merkez"], *KUR["F1"]["bant"]))
print("=" * 78)
Rh, Rs = K1["Hkeskin"]["eta"]["R"], K1["son"]["eta"]["R"]
Rv, Rv_sac, Rv_l = ort_ve_sacilim(lambda v: K1[v]["eta"]["R"])
lamR = Rh - Rv
gerekR = KUR["F1"]["gerek_Lambda_R"]
print("  R_η(vekil) = %s  ⇒ ort %.5f   tohum saçılımı %.5f"
      % (" / ".join("%.5f" % x for x in Rv_l), Rv, Rv_sac))
print("  bant içinde mi? %s      ölüm (≥1.15)? %s"
      % ("✓" if KUR["F1"]["bant"][0] <= Rv <= KUR["F1"]["bant"][1] else "✗",
         "✗ ÖLDÜ" if Rv >= 1.15 else "✓ ölmedi"))
print("  Λ_R = R_η(Hk) − R_η(vekil) = %.5f − %.5f = %+.7f   gerek ≥ %+.7f"
      % (Rh, Rv, lamR, gerekR))
h1, ok1 = hukum("F1", lamR, Rv_sac, gerekR, ">=", "✗ GEREK ŞART TUTMADI",
                "✓ GEREK ŞART TUTTU")
print("  ⇒ F1 %s" % h1)
OUT["F1"] = dict(R_vekil=Rv_l, ort=Rv, sacilim=Rv_sac, R_Hk=Rh, R_son=Rs,
                 Lambda_R=lamR, gerek=gerekR, bant_ici=bool(
                     KUR["F1"]["bant"][0] <= Rv <= KUR["F1"]["bant"][1]),
                 olum=bool(Rv >= 1.15), hukum=h1, gecti=ok1,
                 merkez=KUR["F1"]["merkez"])

# ═══════════════════════════════════════════════════════════════════════
# F2 — g_E(vekil)
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 78)
print("F2 — g_E(vekil):  bant [%.2f, %.2f]  ölüm ≥ %.6f"
      % (*KUR["F2"]["bant"], Hk["gE"]))
print("=" * 78)
gv, gv_sac, gv_l = ort_ve_sacilim(lambda v: GV[v]["gE"])
lamE = math.log(Hk["gE"]) - math.log(gv)
lamE_l = [math.log(Hk["gE"] / x) for x in gv_l]
gerekE = KUR["F2"]["gerek_Lambda_E"]
print("  g_E(vekil) = %s  ⇒ ort %.6f   tohum saçılımı %.6f"
      % (" / ".join("%.6f" % x for x in gv_l), gv, gv_sac))
print("  bant içinde mi? %s      ölüm (≥ g_E(Hk))? %s"
      % ("✓" if KUR["F2"]["bant"][0] <= gv <= KUR["F2"]["bant"][1] else "✗",
         "✗ ÖLDÜ" if gv >= Hk["gE"] else "✓ ölmedi"))
print("  Λ_E = log g_E(Hk) − log g_E(vekil) = %+.7f  (tohumlar %s)   "
      "gerek ≥ %+.7f"
      % (lamE, " / ".join("%+.7f" % x for x in lamE_l), gerekE))
h2, ok2 = hukum("F2", lamE, abs(lamE_l[0] - lamE_l[-1]), gerekE, ">=",
                "✗ GEREK ŞART TUTMADI", "✓ GEREK ŞART TUTTU")
print("  ⇒ F2 %s" % h2)
OUT["F2"] = dict(gE_vekil=gv_l, ort=gv, sacilim=gv_sac, gE_Hk=Hk["gE"],
                 gE_son=sn["gE"], Lambda_E=lamE, Lambda_E_tohum=lamE_l,
                 gerek=gerekE,
                 bant_ici=bool(KUR["F2"]["bant"][0] <= gv
                               <= KUR["F2"]["bant"][1]),
                 olum=bool(gv >= Hk["gE"]), hukum=h2, gecti=ok2)

# ═══════════════════════════════════════════════════════════════════════
# F3 — θ(vekil)
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 78)
print("F3 — θ(vekil):  (a) ≤ %.7f H-F1b YAŞAR | (b) ≥ %.7f H-F1b ÖLÜR   "
      "[ön-kayıtlı beklenti: b]"
      % (KUR["F3"]["dal_a_yasar"], KUR["F3"]["dal_b_olur"]))
print("=" * 78)
tv, tv_sac, tv_l = ort_ve_sacilim(lambda v: GV[v]["th"])
th_a, th_b = KUR["F3"]["dal_a_yasar"], KUR["F3"]["dal_b_olur"]
if tv <= th_a:
    dal, dal_ad = "a", "H-F1b YAŞAR (θ fazlası çöktü)"
elif tv >= th_b:
    dal, dal_ad = "b", "H-F1b ÖLDÜ (θ fazlası çökmedi)"
else:
    dal, dal_ad = "ara", "KISMİ (ne yaşar ne ölür)"
print("  θ(vekil) = %s  ⇒ ort %.6f   tohum saçılımı %.6f"
      % (" / ".join("%.6f" % x for x in tv_l), tv, tv_sac))
print("  θ(Hkeskin) = %.6f   θ(son) = %.6f   θ(HA4) = %.6f"
      % (Hk["th"], sn["th"], ha["th"]))
print("  ⇒ DAL = %s : %s" % (dal.upper(), dal_ad))
# F7 (dondurulmuş): |y1−y2| eşiğe uzaklıktan büyükse madde HÜKÜMSÜZ
mes_a, mes_b = abs(tv - th_a), abs(tv - th_b)
f7a = tv_sac <= mes_a
f7b = tv_sac <= mes_b
print("  F7 (tohum saçılımı %.6f):" % tv_sac)
print("     (a) eşiğine uzaklık %.6f ⇒ %s"
      % (mes_a, "kesin" if f7a else "**HÜKÜMSÜZ (tohum gürültüsü)**"))
print("     (b) eşiğine uzaklık %.6f ⇒ %s"
      % (mes_b, "kesin" if f7b else "**HÜKÜMSÜZ (tohum gürültüsü)**"))
if dal == "ara" and f7b:
    print("  ⇒ (b) KESİN OLARAK DIŞLANDI: ön-kayıtlı BEKLENTİ (b) YANLIŞ "
          "çıktı — θ çöktü, ama (a) eşiğine tohum gürültüsü içinde "
          "ulaşamadı.")
OUT["F3"] = dict(th_vekil=tv_l, ort=tv, sacilim=tv_sac, dal=dal,
                 hukum=dal_ad, a=th_a, b=th_b, th_Hk=Hk["th"],
                 th_son=sn["th"], th_HA4=ha["th"],
                 beklenti=KUR["F3"]["beklenti"],
                 F7_a_kesin=bool(f7a), F7_b_kesin=bool(f7b),
                 mesafe_a=mes_a, mesafe_b=mes_b)

# ═══════════════════════════════════════════════════════════════════════
# F4 — ΔlogM(vekil ← Hkeskin)
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 78)
print("F4 — Δlog M(vekil←Hk):  bant [%.2f, %.2f]  ölüm ≥ 0   gerek Λ_M ≥ %.7f"
      % (*KUR["F4"]["bant"], KUR["F4"]["gerek_Lambda_M"]))
print("=" * 78)
dM_l = [math.log(GV[v]["KAL"] / Hk["KAL"]) for v in VEK]
dM = float(np.mean(dM_l))
dM_sac = abs(dM_l[0] - dM_l[-1]) if len(dM_l) > 1 else 0.0
lamM = -dM
gerekM = KUR["F4"]["gerek_Lambda_M"]
print("  M(vekil) = %s   M(Hk) = %.6f"
      % (" / ".join("%.6f" % GV[v]["KAL"] for v in VEK), Hk["KAL"]))
print("  Δlog M(vekil←Hk) = %s ⇒ ort %+.7f  (saçılım %.7f)"
      % (" / ".join("%+.7f" % x for x in dM_l), dM, dM_sac))
print("  Λ_M = %+.7f   gerek ≥ %+.7f   bant içinde mi? %s"
      % (lamM, gerekM,
         "✓" if KUR["F4"]["bant"][0] <= dM <= KUR["F4"]["bant"][1] else "✗"))
h4, ok4 = hukum("F4", lamM, dM_sac, gerekM, ">=", "✗ GEREK ŞART TUTMADI",
                "✓ GEREK ŞART TUTTU")
print("  ⇒ F4 %s   (ölüm ΔlogM ≥ 0 : %s)"
      % (h4, "✗ ÖLDÜ" if dM >= 0 else "✓ ölmedi"))
OUT["F4"] = dict(dM_tohum=dM_l, ort=dM, sacilim=dM_sac, Lambda_M=lamM,
                 gerek=gerekM, hukum=h4, gecti=ok4, olum=bool(dM >= 0),
                 bant_ici=bool(KUR["F4"]["bant"][0] <= dM
                               <= KUR["F4"]["bant"][1]))

# ═══════════════════════════════════════════════════════════════════════
# F5 — ÜÇÜNCÜ MOMENTLERİN RICE SIFIRI
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 78)
print("F5 — RICE SIFIRI: |m3_çizgi| ≤ %.3f  ve  |skew(η_çizgi)| ≤ %.3f"
      % (KUR["F5"]["m3_cizgi_esik"], KUR["F5"]["skew_eta_ciz_esik"]))
print("=" * 78)
print("%-9s %11s %11s %11s %11s %11s %11s"
      % ("gaz", "m3", "m3_çizgi", "skew(η_çiz)", "skew(e1)", "skew(x1)",
         "skew(ds)"))
for g in ["son", "Hkeskin", "HA4"] + VEK:
    m = K1[g]["m3"]
    print("%-9s %+11.5f %+11.5f %+11.5f %+11.5f %+11.5f %+11.5f"
          % (g, m["olc"], m["cizgi"], m["skew_eta_ciz"], m["skew_e1"],
             m["skew_x1"], m["skew_ds"]))
mc, mc_sac, mc_l = ort_ve_sacilim(lambda v: K1[v]["m3"]["cizgi"])
sc, sc_sac, sc_l = ort_ve_sacilim(lambda v: K1[v]["m3"]["skew_eta_ciz"])
f5a = abs(mc) <= KUR["F5"]["m3_cizgi_esik"]
f5b = abs(sc) <= KUR["F5"]["skew_eta_ciz_esik"]
print("  |m3_çizgi|(vekil ort) = %.5f  ≤ %.3f ? %s   [saçılım %.5f]"
      % (abs(mc), KUR["F5"]["m3_cizgi_esik"], "✓" if f5a else "✗ ÖLDÜ",
         mc_sac))
print("  |skew(η_çizgi)|(vekil ort) = %.5f ≤ %.3f ? %s   [saçılım %.5f]"
      % (abs(sc), KUR["F5"]["skew_eta_ciz_esik"], "✓" if f5b else "✗ ÖLDÜ",
         sc_sac))
print("  (skew(ds) için ÖN-KAYIT YOKTUR — yalnız yazılır: vekil %s)"
      % " / ".join("%+.5f" % K1[v]["m3"]["skew_ds"] for v in VEK))
OUT["F5"] = dict(m3_cizgi=mc_l, m3_cizgi_ort=mc, m3_sac=mc_sac,
                 skew_eta_ciz=sc_l, skew_ort=sc, skew_sac=sc_sac,
                 gecti_m3=bool(f5a), gecti_skew=bool(f5b),
                 skew_ds=[K1[v]["m3"]["skew_ds"] for v in VEK])

# ═══════════════════════════════════════════════════════════════════════
# F6 — DC KAÇAĞI KİLİDE KÖR MÜ (T-3)
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 78)
print("F6 — T-3 (DC kaçağı kilide kördür):  ort(E) > 0 ve "
      "|ort(E)/%.7f − 1| ≤ %.2f ; ⟨cosΔφ⟩ ≥ %.2f (0.50–0.80)"
      % (KUR["F6"]["ortE_capa"], KUR["F6"]["bagil_tolerans"],
         KUR["F6"]["cos_esik"]))
print("=" * 78)
print("%-9s %12s %12s   %s"
      % ("gaz", "ort(E)", "μ̂²_E", "⟨cosΔφ⟩ bant bant (0.4→0.95)"))
for g in ["son", "Hkeskin", "HA4"] + VEK:
    b = [x for x in K3[g]["bant"] if x["lo"] >= 0.40 - 1e-9]
    print("%-9s %+12.6f %12.6f   %s"
          % (g, K3[g]["ortE"], K3[g]["mu2"],
             " ".join("%+.3f" % x["cos"] for x in b)))
ov, ov_sac, ov_l = ort_ve_sacilim(lambda v: K3[v]["ortE"])
capa = KUR["F6"]["ortE_capa"]
bagil = abs(ov / capa - 1)
cosb = []
for v in VEK:
    for x in K3[v]["bant"]:
        if x["lo"] >= 0.50 - 1e-9 and x["hi"] <= 0.80 + 1e-9:
            cosb.append(x["cos"])
f6a = ov > 0
f6b = bagil <= KUR["F6"]["bagil_tolerans"]
f6c = bool(cosb) and min(cosb) >= KUR["F6"]["cos_esik"]
print("  ort(E)(vekil) = %s ⇒ ort %+.6f  (saçılım %.6f)"
      % (" / ".join("%+.6f" % x for x in ov_l), ov, ov_sac))
print("  |ort(E)/çapa − 1| = %.4f  ≤ 0.30 ? %s      ort(E) > 0 ? %s"
      % (bagil, "✓" if f6b else "✗", "✓" if f6a else "✗"))
print("  ⟨cosΔφ⟩ (0.50–0.80 bantları, iki tohum) = %s ⇒ min %.4f ≥ 0.90 ? %s"
      % (" ".join("%+.4f" % c for c in cosb),
         min(cosb) if cosb else float("nan"), "✓" if f6c else "✗"))
f6 = f6a and f6b and f6c
print("  ⇒ F6 / T-3: %s" % ("✓ TUTTU — DC kaçağı faz kilidine KÖRDÜR"
                            if f6 else "✗ T-3 ÖLDÜ"))
OUT["F6"] = dict(ortE=ov_l, ort=ov, sacilim=ov_sac, capa=capa, bagil=bagil,
                 cos=cosb, gecti=bool(f6), a=bool(f6a), b=bool(f6b),
                 c=bool(f6c))

# ═══════════════════════════════════════════════════════════════════════
# F8 — H-F2: KİLİT PAYININ İLK-İLKE (kaba kalem) ÖNGÖRÜSÜ
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 78)
print("F8 — H-F2 kaba kalem:  kilit^HF2 = Δ_faz/(Δ_amp+Δ_res+Δ_faz), τ>0.70")
print("=" * 78)


def uclu(g, lo=0.70, hi=1.01):
    d = K3[g]
    tau = np.array(d["tau"])
    hp = np.array(d["hp_abs"])
    ka = np.array(d["kap_abs"])
    xi = np.array(d["xi"])
    m = (tau > lo) & (tau <= hi)
    amp = float(hp[m].sum())
    mod = float((hp[m] * ka[m]).sum())
    xs = float(xi[m].sum())
    return amp, mod / amp, xs / mod, xs, int(m.sum())


A_s, R_s, P_s, X_s, n_s = uclu("son")
A_h, R_h, P_h, X_h, n_h = uclu("Hkeskin")
d_amp = math.log(A_s / A_h)
d_res = math.log(R_s / R_h)
d_faz = math.log(P_s / P_h)
tot = d_amp + d_res + d_faz
kilit_hf2 = d_faz / tot
print("  τ > 0.70  (n = %d çizgi)" % n_s)
print("    amp = Σ|hp|        : son %.6f   Hk %.6f   Δlog %+.6f"
      % (A_s, A_h, d_amp))
print("    res = Σ|hp||κ|/Σ|hp|: son %.6e   Hk %.6e   Δlog %+.6f"
      % (R_s, R_h, d_res))
print("    faz = ξ/Σ|hp||κ|   : son %.6f   Hk %.6f   Δlog %+.6f"
      % (P_s, P_h, d_faz))
print("    Σξ(τ>0.70)         : son %+.6e   Hk %+.6e   Δlog %+.6f (=toplam)"
      % (X_s, X_h, math.log(X_s / X_h)))
print("  ⇒ kilit^HF2 = %+.6f / %+.6f = **%.4f**   (hedef 0.61 ± 0.10)"
      % (d_faz, tot, kilit_hf2))
f8 = abs(kilit_hf2 - KUR["F8"]["hedef"]) <= KUR["F8"]["tolerans"]
print("  ⇒ H-F2 %s   (|%.4f − 0.61| = %.4f)"
      % ("✓ YAŞADI" if f8 else "✗ ÖLDÜ", kilit_hf2,
         abs(kilit_hf2 - KUR["F8"]["hedef"])))
# bant bant (bilgi)
print("\n  bant bant (bilgi; ön-kayıt tek bölge diyor):")
print("    %-12s %10s %10s %10s %10s" % ("bant", "Δ_amp", "Δ_res", "Δ_faz",
                                         "kilit"))
bb = []
for lo, hi in ((0.70, 0.80), (0.80, 0.95)):
    a1, r1, p1, x1, _ = uclu("son", lo, hi)
    a0, r0, p0, x0, _ = uclu("Hkeskin", lo, hi)
    da, dr, df = (math.log(a1 / a0), math.log(r1 / r0), math.log(p1 / p0))
    kk = df / (da + dr + df)
    bb.append(dict(lo=lo, hi=hi, d_amp=da, d_res=dr, d_faz=df, kilit=kk))
    print("    %.2f–%.2f    %+10.6f %+10.6f %+10.6f %10.4f"
          % (lo, hi, da, dr, df, kk))
OUT["F8"] = dict(d_amp=d_amp, d_res=d_res, d_faz=d_faz, toplam=tot,
                 kilit_HF2=kilit_hf2, hedef=KUR["F8"]["hedef"],
                 tolerans=KUR["F8"]["tolerans"], gecti=bool(f8), bant=bb,
                 son=dict(amp=A_s, res=R_s, faz=P_s, xi=X_s, n=n_s),
                 Hkeskin=dict(amp=A_h, res=R_h, faz=P_h, xi=X_h, n=n_h))

# ═══════════════════════════════════════════════════════════════════════
# F9 — NİHAİ AYRIŞIM TABLOSU (K3)
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 78)
print("F9 / K3 — NİHAİ AYRIŞIM TABLOSU:  ΔM = KESİM + KİLİT  (parametresiz)")
print("=" * 78)
print("  NOT (dürüstlük): F9'un ön-kayıttaki LİTERAL metni")
print("       kesim_j := f_j·Δ_j(erfc←Hk) ,  f_j := Δ_j(son←Hk)/Δ_j(erfc←Hk)")
print("  dejeneredir (kesim_j ≡ Δ_j(son←Hk), kilit_j ≡ 0). Uygulanan okuma,")
print("  F1 ve F2'nin SAYISAL eşiklerinde zaten yazılı olandır:")
print("       kesim_j = f_j·Δ_j(son←Hk) ,  kilit_j = (1−f_j)·Δ_j(son←Hk)")
print("  (0.0149735 = 0.61×0.0245468 ve 0.0059398 = 0.2555×0.0232450).")
print("  Hiçbir ölüm eşiği değiştirilmedi.\n")

carp = [("g_E", lambda R: R["gE"], 1.0),
        ("g_X² (2·g_X)", lambda R: R["gX"], 2.0),
        ("θ", lambda R: R["th"], 1.0),
        ("M = g_E g_X² θ", lambda R: R["KAL"], 1.0)]
print("%-15s %11s %11s %8s %11s %11s %11s %8s"
      % ("çarpan", "Δ(son←Hk)", "Δ(erfc←Hk)", "f_j", "KESİM_j", "KİLİT_j",
         "Δ(Hk←vekil)", "ω_j"))
tab = []
for ad, f, k in carp:
    ds_ = k * math.log(f(sn) / f(Hk))
    de_ = k * math.log(f(ha) / f(Hk))
    fj = ds_ / de_ if de_ else float("nan")
    kes = fj * ds_
    kil = (1.0 - fj) * ds_
    dv_l = [k * math.log(f(Hk) / f(GV[v])) for v in VEK]
    dv = float(np.mean(dv_l))
    om = kil / dv if dv else float("nan")
    om_l = [kil / x if x else float("nan") for x in dv_l]
    sac = abs(dv_l[0] - dv_l[-1]) if len(dv_l) > 1 else 0.0
    tab.append(dict(ad=ad, d_son=ds_, d_erfc=de_, f=fj, kesim=kes, kilit=kil,
                    d_vekil=dv, d_vekil_tohum=dv_l, sacilim=sac,
                    omega=om, omega_tohum=om_l,
                    odenebilir=bool(0.0 < om <= 1.0),
                    F7_kesin=bool(sac <= abs(dv))))
    print("%-15s %+11.6f %+11.6f %+8.4f %+11.6f %+11.6f %+11.6f %+8.4f"
          % (ad, ds_, de_, fj, kes, kil, dv, om))
    print("%-15s %11s %11s %8s %11s %11s %11s %8s"
          % ("", "", "", "", "", "tohumlar:",
             "/".join("%+.5f" % x for x in dv_l),
             "/".join("%+.3f" % x for x in om_l)))
kal = (tab[0]["d_son"] + tab[1]["d_son"] + tab[2]["d_son"] - tab[3]["d_son"])
print("  özdeşlik kalıntısı Δlog M − (Δlog g_E + 2Δlog g_X + Δlog θ) = %.2e"
      % kal)
kalv = [(math.log(Hk["gE"] / GV[v]["gE"])
         + 2 * math.log(Hk["gX"] / GV[v]["gX"])
         + math.log(Hk["th"] / GV[v]["th"])
         - math.log(Hk["KAL"] / GV[v]["KAL"])) for v in VEK]
print("  aynı özdeşlik VEKİL sütununda: %s"
      % " / ".join("%.2e" % x for x in kalv))

# E-payı ve θ-payı ayrı satır (KALEM K3)
print("\n  ΔM = %+.6f  (=%+.2f%%) ayrışımı — E-payı ve θ-payı AYRI:"
      % (tab[3]["d_son"], 100 * (math.exp(tab[3]["d_son"]) - 1)))
E_kes = tab[0]["kesim"] + tab[1]["kesim"]
E_kil = tab[0]["kilit"] + tab[1]["kilit"]
T_kes, T_kil = tab[2]["kesim"], tab[2]["kilit"]
tot_dM = tab[3]["d_son"]
print("    %-22s %11s %11s %9s" % ("", "KESİM", "KİLİT", "toplam"))
for ad, a, b in (("E kanalı (g_E, g_X²)", E_kes, E_kil),
                 ("θ kanalı", T_kes, T_kil),
                 ("TOPLAM (M)", E_kes + T_kes, E_kil + T_kil)):
    print("    %-22s %+11.6f %+11.6f %+9.6f" % (ad, a, b, a + b))
print("    ΔM'nin yüzdesi:        KESİM %+.1f%%   KİLİT %+.1f%%"
      % (100 * (E_kes + T_kes) / tot_dM, 100 * (E_kil + T_kil) / tot_dM))
OUT["F9"] = dict(tablo=tab, kalinti=kal, kalinti_vekil=kalv,
                 E_kesim=E_kes, E_kilit=E_kil, th_kesim=T_kes,
                 th_kilit=T_kil, dM=tot_dM)

# R_η satırı (162 makinesi, DOĞRUSAL)
Ra = K1["HA4"]["eta"]["R"]
fR = (Rs - Rh) / (Ra - Rh)
print("\n  R_η satırı (162 makinesi, doğrusal):")
print("    R_η: Hk %.5f  son %.5f  HA4 %.5f  ⇒ f = %+.5f" % (Rh, Rs, Ra, fR))
print("    KESİM = %+.7f   KİLİT = %+.7f   Λ_R(Hk←vekil) = %+.7f   ω = %+.4f"
      % (fR * (Rs - Rh), (1 - fR) * (Rs - Rh), lamR,
         (1 - fR) * (Rs - Rh) / lamR if lamR else float("nan")))
OUT["F9"]["R_eta"] = dict(f=fR, kesim=fR * (Rs - Rh),
                          kilit=(1 - fR) * (Rs - Rh), d_vekil=lamR,
                          omega=(1 - fR) * (Rs - Rh) / lamR if lamR else None)

# ═══════════════════════════════════════════════════════════════════════
# MÜHÜR
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 78)
print("MÜHÜR KURALI (ön-kayıtta dondurulmuş)")
print("=" * 78)
if not ok2:
    muhur = "ÖLÜM"
    aciklama = "F2 ✗ — kilit faturası ödemiyor."
elif ok1 and ok2 and dal == "a" and ok4:
    muhur = "MÜHÜR-TAM"
    aciklama = ONK["kural"]["MUHUR"]["cumle"]
elif ok1 and ok2 and dal == "b":
    muhur = "MÜHÜR-E"
    aciklama = ("Cümle YALNIZ E kanalı için mühürlenir; θ 'üçüncü eksen' "
                "olarak keskin bir negatifle açık kalır.")
else:
    muhur = "KISMİ"
    aciklama = ("F1=%s F2=%s F3=%s F4=%s — ön-kayıtlı mühür dallarının "
                "hiçbiri tam karşılanmadı." % (ok1, ok2, dal, ok4))
print("  F1 %s | F2 %s | F3 dal=%s | F4 %s | F5 %s/%s | F6 %s | F8 %s"
      % (ok1, ok2, dal, ok4, f5a, f5b, f6, f8))
print("  ⇒ **%s**" % muhur)
print("  %s" % aciklama)
OUT["MUHUR"] = dict(sonuc=muhur, aciklama=aciklama, F1=ok1, F2=ok2,
                    F3=dal, F4=ok4, F5=[bool(f5a), bool(f5b)], F6=bool(f6),
                    F8=bool(f8))

(S176 / "K2.json").write_text(json.dumps(OUT, indent=1, ensure_ascii=False))
print("\n-> %s" % (S176 / "K2.json"))
