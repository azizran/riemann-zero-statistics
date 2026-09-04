"""
173e — EĞRİ UÇLARI: DOKUZ NOKTALI λ MERDİVENİ (0.40 … 1.45)
============================================================
Yeni ölçüm YOK; 173c'nin C_L040/C_L145'i + önbellekteki yedi λ gazı.
`169_k1` AYNEN import edilir; ν tanımı 171k ile birebir
(ν = Δlog KALİB/Δlog W_X, orta bant ve 5-bant ortalaması, jackknife).

YÜZLEŞİLEN ÖN-KAYITLAR — `173b_onkayit.py` docstring'i (§C), ölçümden
ÖNCE yazıldı; burada DEĞİŞTİRİLMEDEN sınanır:
  P1 ν(0.50→0.40) ≈ 1.99 ∈ [1.55, 2.45];  ν(1.30→1.45) ≈ 0.58 ∈ [0.30, 0.90]
     (ikisi de: c(0.40) > c(0.50) ve c(1.45) > c(1.30))
  P2 dokuz noktada kesişim SAYISI 1, yeri λ ∈ [0.63, 0.67];
     fark(0.45) ≈ −0.60 [−0.85, −0.40];  fark(1.375) ≈ +0.24 [+0.08, +0.42]
  P3 α(0.40) ∈ [0.85, 1.15];  α(1.45) ∈ [0.35, 0.75];
     şekil rms: λ=1.45 için > %4, λ=0.40 için ≤ %2.5;  σ*/2(1.45) < 0.18661
  P4 Q_E(0.40) ∈ [0.54, 0.68];  Q_E(1.45) ∈ [0.70, 0.83];
     dokuz noktalı tümsek tepesi λ_c = 1.0109 ± 0.05
  P5 θ taşıyıcısı artığı |·| ≤ %3 (173d'de de basılır)

(λ = 1.40 sadakat-sınırı gazı sonradan eklendiği için merdiven ON
noktalıdır; ön-kayıtlı ADIMLAR — özellikle ν(1.30→1.45) — ara nokta
eklense de AYNEN, uçtan uca hesaplanır.)

Kullanım: 173e_egri_uclari.py
Çıktı:    scratchpad/173/EGRI.json

SONUÇ (yalnız gerçek koşudan; scratchpad/173/log_173e.txt):
 P1 ✓ ν(0.50→0.40) = +2.203 ± 0.056, ön-kayıt [1.55, 2.45] İÇERİDE;
    c(0.40) = 0.4098 > c(0.50) = 0.3800 ✓ (171 §T2c-d ikinci kez).
    ✗ ν(1.30→1.45) = +1.697 (ön-kayıt [0.30, 0.90]) ve c(1.45) < c(1.30).
    YENİ: ν(1.30→1.40) = +1.541 ± 0.103 (+5.3σ) — SAĞLIKLI adım,
    yani c'nin üst-uç DÖNÜŞÜ gerçek, sadakat artefaktı değil.
 P2 yer ✓✓✓ / sayı ✗: kesişimler λ = **0.6486** ve **1.2787**.
    Birincisi 171'in bağımsız λ* = 0.6487'sinden 0.0001 fark (172'nin
    0.6506'sı ‰3 idi; şimdi ‰0.15). İkincisi bir TEPE ve ön-kaydım
    "kesişim sayısı 1" demişti ⇒ ÖLDÜ. fark(0.45) = −0.7552 ✓.
 P3 ✗✗ α(0.40) = 0.6474 (ön-kayıt [0.85, 1.15]), α(1.45) = 0.1866
    ([0.35, 0.75]), şekil rms(0.40) = %3.69 (≤%2.5 dedim) — üçü de ıska.
    AMA 171 §T2c-c'nin [0.50, 1.15] penceresi İKİ YANDAN DOĞRULANDI
    (rms: 0.40'ta %3.69, 1.30'da %3.61, 1.40'ta %8.59, 1.45'te %9.29).
    α tepesi uyum penceresine bağımlı: 5 orta 0.6675 ✓ / 7 nokta 0.7877
    / 9 nokta 0.8058 ⇒ 172'nin "5 orta nokta" şartı GEREKLİLİKMİŞ.
 P4 ✓/✗ Q_E(0.40) = 0.56708 ∈ [0.54, 0.68] ✓ (tümseğin SOL kolu
    doğrulandı); Q_E(1.45) = 0.55489 ∉ [0.70, 0.83] ✗ ve Q_E(1.40) =
    0.66240 de ✗ — sağ kol SAĞLIKLI gazda çöküyor. Tepe: 5 orta 1.0191,
    7 nokta 1.0009 (λ_c = 1.0109 ✓✓), 9 nokta 0.9407 (sağlıksız kaldıraç).
 P6 M9 alt uçta YAŞIYOR (−0.1σ / +0.9σ), λ = 1.40'ta (SAĞLIKLI) ÖLÜYOR
    (−5.5σ / −4.4σ). Geçerlilik penceresi λ ∈ [0.40, 1.30].
 ÖZDEŞLİK DENETİMİ: g = 1 − Q/ρ on gazın hepsinde |ΔK|/K ≤ 5.3e−15.
"""
import importlib
import json
import math
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
LAM_C, LSTAR, P_TH = 1.9147 / 1.894, 0.6486573321625787, 0.651
SIRA = ["L040", "L050", "L060", "L070", "L085", "Hkeskin", "L115", "L130",
        "L140", "L145"]
LAMS = [0.40, 0.50, 0.60, 0.70, 0.85, 1.00, 1.15, 1.30, 1.40, 1.45]


def kusur(V_M, V_O, V_R, kor, g):
    Kv = 0.5 * (V_O + V_M - V_R)
    S = Kv / g
    return dict(mu2=(S - V_M) / V_O, rho=S / V_O, Q=(S - Kv) / V_O,
                r=V_R / V_O, dK=abs(Kv - kor * math.sqrt(V_M * V_O))
                / abs(Kv))


LOS = [0.52, 0.56, 0.60, 0.64, 0.68]


def pencere(d):
    """173d.pencere ile aynı: beş hüküm bandı + 169 filtresinin bayrağı."""
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


def main():
    os.makedirs(SCR173, exist_ok=True)
    D = K1.yukle()
    sira = [g for g in SIRA if g in D]
    lams = [LAMS[SIRA.index(g)] for g in sira]
    B, SAG = {}, {}
    for g in sira + ["son"]:
        ban, sg = pencere(D[g])
        SAG[g] = bool(all(sg))
        rows = [K1.band_jk(b["cizgi"]) for b in ban]
        B[g] = dict(R=np.array([b["R_bant"] for b in ban]),
                    K=np.array([r["KALIB_u2"] for r in rows]),
                    sK=np.array([r["sKALIB_u2"] for r in rows]),
                    W=np.array([r["W_X"] for r in rows]),
                    sW=np.array([r["sW_X"] for r in rows]),
                    t=np.array([r["tau_eff"] for r in rows]),
                    c=np.array([r["c_WX"] for r in rows]),
                    sc=np.array([r["sc_WX"] for r in rows]))
    HK = B["Hkeskin"]
    sX_hk = D["Hkeskin"]["sigX"]
    A2 = np.exp(-2 * np.pi ** 2 * HK["t"] ** 2 * sX_hk ** 2)

    print("=" * 112)
    print("173e — DOKUZ NOKTALI λ MERDİVENİ (%s)" % ", ".join(sira))
    print("=" * 112)

    # ---------------- (0) DEFTER ---------------------------------------
    D9 = {}
    print("  (SAĞLIKSIZ = 169'un R_bant ≥ 0.98 filtresi düştü; filtre "
          "gevşetilmedi, o satır HÜKÜM DEĞİL TEŞHİStir)")
    print("\n  %-9s %5s | %8s %7s | %8s | %7s %7s %7s | %7s %7s | %7s %7s %s"
          % ("gaz", "λ", "c_WX", "±", "M", "g_E", "g_X", "θ", "Q_E", "ρ_X",
             "α_g", "σ_X̃", "R_bant"))
    for g, l in zip(sira, lams):
        G, a_ = B[g], D[g]["artik"]
        lg = np.log(G["c"])
        cg = float(np.exp(lg.mean()))
        s_jk = float(np.sqrt(np.sum((G["sc"] / G["c"]) ** 2)) / len(G["c"])) * cg
        s_bt = float(np.std(lg, ddof=1) / np.sqrt(len(lg))) * cg
        Mv = float(np.exp(np.mean(np.log(G["K"] / HK["K"]))))
        sM = float(np.hypot(
            np.sqrt(np.sum((G["sK"] / G["K"]) ** 2
                           + (HK["sK"] / HK["K"]) ** 2)) / len(G["K"]) * Mv,
            np.std(np.log(G["K"] / HK["K"]), ddof=1) / np.sqrt(len(G["K"]))
            * Mv))
        gcal = a_["gE"] * a_["gX"] ** 2
        orta = [b for b in D[g]["bant"]
                if b.get("olculdu") and abs(b["lo"] - 0.60) < 1e-9][0]
        th = orta["KALIB_u2"] / gcal
        E_ = kusur(a_["varE_mod"], a_["varE_olc"], a_["varE_res"], a_["korE"],
                   a_["gE"])
        X_ = kusur(a_["varX_mod"], a_["varX_olc"], a_["varX_res"], a_["korX"],
                   a_["gX"])
        alfa = float(-np.polyfit(G["t"] ** 2, np.log(G["K"]), 1)[0])
        sekil = 100 * ((G["K"] / G["K"][2]) / (A2 / A2[2]) - 1)
        D9[g] = dict(lam=l, c=cg, s_jk=s_jk, s_bant=s_bt,
                     s_tot=float(np.hypot(s_jk, s_bt)), M=Mv, sM=sM,
                     gE=a_["gE"], gX=a_["gX"], gcal=gcal, theta=float(th),
                     QE=E_["Q"], rhoE=E_["rho"], rE=E_["r"], mu2E=E_["mu2"],
                     QX=X_["Q"], rhoX=X_["rho"], rX=X_["r"],
                     dK_E=E_["dK"], dK_X=X_["dK"],
                     alfa=alfa, sX_eff=float(math.sqrt(alfa / (2 * np.pi ** 2))),
                     sigX=D[g]["sigX"], sigds=D[g]["sigds"],
                     sigC=D[g]["sigC"],
                     KAL=[float(x) for x in G["K"]],
                     WX=[float(x) for x in G["W"]],
                     tau=[float(x) for x in G["t"]],
                     sekil=[float(x) for x in sekil],
                     sekil_rms=float(np.sqrt((sekil ** 2).mean())),
                     saglikli=SAG[g],
                     R_bant=[float(x) for x in G["R"]])
        print("  %-9s %5.2f | %8.4f %7.4f | %8.4f | %7.4f %7.4f %7.4f | "
              "%7.4f %7.4f | %7.4f %7.5f  %s"
              % (g, l, cg, np.hypot(s_jk, s_bt), Mv, a_["gE"], a_["gX"], th,
                 E_["Q"], X_["rho"], alfa, D[g]["sigX"],
                 "%.2f-%.2f" % (G["R"].min(), G["R"].max())
                 + ("" if SAG[g] else "  ← SAĞLIKSIZ")))

    lam = np.array(lams)
    KALm = np.array([np.exp(np.mean(np.log(B[g]["K"]))) for g in sira])
    WXm = np.array([np.exp(np.mean(np.log(B[g]["W"]))) for g in sira])
    cO = np.array([D9[g]["c"] for g in sira])
    MOl = np.array([D9[g]["M"] for g in sira])
    alfa = np.array([D9[g]["alfa"] for g in sira])
    QE = np.array([D9[g]["QE"] for g in sira])

    # ---------------- (1) ν MERDİVENİ (P1) -----------------------------
    print("\n" + "=" * 112)
    print("P1 — ν MERDİVENİ  (ν = Δlog KALİB/Δlog W_X; 171k ile birebir)")
    print("=" * 112)
    def nu_cift(a, b):
        dk = np.log(B[b]["K"] / B[a]["K"])
        dw = np.log(B[b]["W"] / B[a]["W"])
        sdk = np.sqrt((B[a]["sK"] / B[a]["K"]) ** 2
                      + (B[b]["sK"] / B[b]["K"]) ** 2)
        la, lb = D9[a]["lam"], D9[b]["lam"]
        return dict(dan=a, ye=b, lam_a=la, lam_b=lb,
                    lam_orta=0.5 * (la + lb),
                    nu_orta=float(dk[2] / dw[2]),
                    nu_ort=float(dk.mean() / dw.mean()),
                    s_orta=float(sdk[2] / abs(dw[2])),
                    s_ort=float(np.sqrt((sdk ** 2).sum()) / len(sdk)
                                / abs(dw.mean())))

    nu = []
    for a, b in zip(sira[:-1], sira[1:]):
        r = nu_cift(a, b)
        r["z1"] = float((r["nu_ort"] - 1.0) / r["s_ort"])
        nu.append(r)
        print("  %-8s → %-8s (λ %.2f → %.2f, orta %.3f): ν = %+.3f ± %.3f "
              "(orta bant)  %+.3f ± %.3f (5 bant)  [ν−1 = %+.1fσ]"
              % (a, b, r["lam_a"], r["lam_b"], r["lam_orta"], r["nu_orta"],
                 r["s_orta"], r["nu_ort"], r["s_ort"], r["z1"]))
    P1 = {}
    for et, (a, b), band in (("ν(0.50→0.40)", ("L050", "L040"), (1.55, 2.45)),
                             ("ν(1.30→1.45)", ("L130", "L145"), (0.30, 0.90))):
        if a not in B or b not in B:
            continue
        r = [nu_cift(a, b)]      # ara nokta eklense de ön-kayıtlı ADIM aynı
        v = r[0]["nu_ort"]
        P1[et] = dict(olcum=v, band=band,
                      hukum="İÇERİDE" if band[0] <= v <= band[1] else "ISKA")
        print("  ** %s: ölçüm %+.3f ± %.3f, ön-kayıt [%.2f, %.2f] ⇒ %s"
              % (et, v, r[0]["s_ort"], band[0], band[1], P1[et]["hukum"]))
    if "L040" in D9:
        P1["c040_vs_c050"] = dict(c040=D9["L040"]["c"], c050=D9["L050"]["c"],
                                  artti=bool(D9["L040"]["c"] > D9["L050"]["c"]))
        print("  ** c(0.40) = %.4f vs c(0.50) = %.4f ⇒ ön-kayıt 'artar' %s"
              % (D9["L040"]["c"], D9["L050"]["c"],
                 "✓" if D9["L040"]["c"] > D9["L050"]["c"] else "✗"))
    if "L145" in D9:
        P1["c145_vs_c130"] = dict(c145=D9["L145"]["c"], c130=D9["L130"]["c"],
                                  artti=bool(D9["L145"]["c"] > D9["L130"]["c"]))
        print("  ** c(1.45) = %.4f vs c(1.30) = %.4f ⇒ ön-kayıt 'artar' %s"
              % (D9["L145"]["c"], D9["L130"]["c"],
                 "✓" if D9["L145"]["c"] > D9["L130"]["c"] else "✗"))

    # ---------------- (2) VADİ / EĞİM KESİŞMESİ (P2) --------------------
    print("\n" + "=" * 112)
    print("P2 — VADİ: dlogKALİB/dλ = dlogW_X/dλ  (172 §G4.1 ile aynı cebir)")
    print("=" * 112)
    lm = 0.5 * (lam[:-1] + lam[1:])
    dK = np.diff(np.log(KALm)) / np.diff(lam)
    dW = np.diff(np.log(WXm)) / np.diff(lam)
    f = dK - dW
    print("  λ_orta   dlogKALİB/dλ  dlogW_X/dλ   fark (= dlog c/dλ)")
    for i in range(len(lm)):
        print("  %.3f     %+9.4f     %+9.4f    %+9.4f" % (lm[i], dK[i], dW[i],
                                                          f[i]))
    kes = [float(lm[i] + (lm[i + 1] - lm[i]) * (-f[i]) / (f[i + 1] - f[i]))
           for i in range(len(f) - 1) if f[i] * f[i + 1] < 0]
    print("  KESİŞİM(ler): %s   (ön-kayıt: SAYI 1, yer [0.63, 0.67]; "
          "171 λ* = %.4f)" % (["%.4f" % k for k in kes], LSTAR))
    P2 = dict(lam_orta=[float(x) for x in lm], dK=[float(x) for x in dK],
              dW=[float(x) for x in dW], fark=[float(x) for x in f],
              kesisim=kes,
              hukum=("TUTTU" if len(kes) == 1 and 0.63 <= kes[0] <= 0.67
                     else "ISKA"))
    print("  ⇒ P2 %s" % P2["hukum"])
    # uç noktalardaki farkın ön-kaydı
    for et, xk, band in (("fark(0.45)", 0.45, (-0.85, -0.40)),
                         ("fark(1.375)", 1.375, (0.08, 0.42))):
        j = int(np.argmin(np.abs(lm - xk)))
        if abs(lm[j] - xk) < 1e-6:
            P2[et] = dict(olcum=float(f[j]), band=band,
                          hukum="İÇERİDE" if band[0] <= f[j] <= band[1]
                          else "ISKA")
            print("  ** %s: ölçüm %+.4f, ön-kayıt [%+.2f, %+.2f] ⇒ %s"
                  % (et, f[j], band[0], band[1], P2[et]["hukum"]))

    # ---------------- (3) α TÜMSEĞİ ve ŞEKİL PENCERESİ (P3) -------------
    print("\n" + "=" * 112)
    print("P3 — α(λ) TÜMSEĞİ ve A-ÇARPANLAŞMASININ PENCERESİ")
    print("=" * 112)
    print("  gaz      λ      α_g     σ*/2      S/A2Hk−1 (%%, 5 bant)"
          "                 rms%%")
    for g in sira:
        print("  %-8s %.2f  %.4f  %.5f   %s   %.2f"
              % (g, D9[g]["lam"], D9[g]["alfa"], D9[g]["sX_eff"],
                 " ".join("%+6.2f" % x for x in D9[g]["sekil"]),
                 D9[g]["sekil_rms"]))
    P3 = {}
    for et, g, band in (("α(0.40)", "L040", (0.85, 1.15)),
                        ("α(1.45)", "L145", (0.35, 0.75))):
        if g in D9:
            v = D9[g]["alfa"]
            P3[et] = dict(olcum=v, band=band,
                          hukum="İÇERİDE" if band[0] <= v <= band[1] else "ISKA")
            print("  ** %s = %.4f, ön-kayıt [%.2f, %.2f] ⇒ %s"
                  % (et, v, band[0], band[1], P3[et]["hukum"]))
    for et, g, esik, yon in (("şekil rms(0.40) ≤ %2.5", "L040", 2.5, "<="),
                             ("şekil rms(1.45) > %4", "L145", 4.0, ">")):
        if g in D9:
            v = D9[g]["sekil_rms"]
            ok = (v <= esik) if yon == "<=" else (v > esik)
            P3[et] = dict(olcum=v, esik=esik, hukum="TUTTU" if ok else "ISKA")
            print("  ** %s: ölçüm %.2f%% ⇒ %s" % (et, v, P3[et]["hukum"]))
    if "L145" in D9:
        v = D9["L145"]["sX_eff"]
        P3["σ*/2(1.45) < 0.18661"] = dict(
            olcum=v, hukum="TUTTU" if v < 0.18661 else "ISKA")
        print("  ** σ*/2(1.45) = %.5f (L130: 0.18661) ⇒ %s"
              % (v, P3["σ*/2(1.45) < 0.18661"]["hukum"]))
    for et, sl in (("5 orta (0.60-1.15)", slice(2, 7)),
                   ("7 nokta (0.50-1.30)", slice(1, 8)),
                   ("9 nokta (0.40-1.45)", slice(None))):
        if len(lam[sl]) < 3:
            continue
        cc = np.polyfit(lam[sl], np.log(alfa[sl]), 2)
        tp = float(-cc[1] / (2 * cc[0]))
        P3["tepe_" + et] = tp
        print("  α tepesi (%-20s) = %+.4f   (171 λ* = %.4f, fark %+.4f)"
              % (et, tp, LSTAR, tp - LSTAR))

    # ---------------- (4) Q_E TÜMSEĞİ (P4) ------------------------------
    print("\n" + "=" * 112)
    print("P4 — Q_E TÜMSEĞİ (172 §G1.5: tepe = inşa doyum eşiği λ_c = %.4f)"
          % LAM_C)
    print("=" * 112)
    print("  gaz      λ      Q_E      ρ_E      rE       μ̂²_E     g_E      "
          "|ΔK|/K")
    for g in sira:
        print("  %-8s %.2f  %.5f  %.5f  %.5f  %+.5f  %.5f  %.1e"
              % (g, D9[g]["lam"], D9[g]["QE"], D9[g]["rhoE"], D9[g]["rE"],
                 D9[g]["mu2E"], D9[g]["gE"], D9[g]["dK_E"]))
    P4 = {}
    for et, g, band in (("Q_E(0.40)", "L040", (0.54, 0.68)),
                        ("Q_E(1.45)", "L145", (0.70, 0.83))):
        if g in D9:
            v = D9[g]["QE"]
            P4[et] = dict(olcum=v, band=band,
                          hukum="İÇERİDE" if band[0] <= v <= band[1] else "ISKA")
            print("  ** %s = %.5f, ön-kayıt [%.2f, %.2f] ⇒ %s"
                  % (et, v, band[0], band[1], P4[et]["hukum"]))
    for et, sl in (("5 orta", slice(2, 7)), ("7 nokta", slice(1, 8)),
                   ("9 nokta", slice(None))):
        if len(lam[sl]) < 3:
            continue
        cc = np.polyfit(lam[sl], np.log(QE[sl]), 2)
        tp = float(-cc[1] / (2 * cc[0]))
        P4["tepe_" + et] = tp
        print("  Q_E tepesi (%-8s) = %+.4f   (λ_c = %.4f, fark %+.4f)"
              % (et, tp, LAM_C, tp - LAM_C))

    # ---------------- (5) M9 ÖRNEKLEM-DIŞI ------------------------------
    print("\n" + "=" * 112)
    print("P6 — M9 EMPİRİK KUADRATİĞİ İKİ UÇTA ÖRNEKLEM-DIŞI YAŞIYOR MU?")
    print("=" * 112)
    L5 = ["L115", "Hkeskin", "L085", "L070", "L060"]
    l5 = np.array([1.15, 1.00, 0.85, 0.70, 0.60])
    M5 = np.array([D9[g]["M"] for g in L5])
    p9 = np.polyfit(np.log(l5), np.log(M5), 2)
    i7 = [sira.index(g) for g in
          ["L050", "L060", "L070", "L085", "Hkeskin", "L115", "L130"]
          if g in sira]
    p9b = np.polyfit(np.log(lam[i7]), np.log(MOl[i7]), 2)
    print("  %-30s %9s %9s %9s %8s" % ("uyum", "λ=0.40", "λ=1.45", "ölçüm",
                                       "sapma"))
    P6 = {}
    for nm, pp in (("M9  (171, 5 gaz: 0.60-1.15)", p9),
                   ("M9′ (172, 7 gaz: 0.50-1.30)", p9b)):
        v = np.exp(np.polyval(pp, np.log([0.40, 1.45])))
        row = dict(M040=float(v[0]), M145=float(v[1]))
        for j, g in (("040", "L040"), ("145", "L145")):
            if g in D9:
                row["d" + j] = 100 * (D9[g]["M"] / row["M" + j] - 1)
                row["z" + j] = (D9[g]["M"] - row["M" + j]) / D9[g]["sM"]
        P6[nm] = row
        print("  %-30s %9.4f %9.4f | %s"
              % (nm, v[0], v[1],
                 "  ".join("%s: %.4f (%+.2f%%, %+.1fσ)"
                           % (g, D9[g]["M"], row["d" + j], row["z" + j])
                           for j, g in (("040", "L040"), ("145", "L145"))
                           if g in D9)))
    # dokuz noktalı yeni kuadratik ve artıkları
    p9c = np.polyfit(np.log(lam), np.log(MOl), 2)
    art = 100 * (MOl / np.exp(np.polyval(p9c, np.log(lam))) - 1)
    P6["9nokta_artik"] = [float(x) for x in art]
    print("  9 noktalı yeniden uyum artıkları (%%): %s   rms %.2f%%"
          % (" ".join("%+5.2f" % x for x in art),
             float(np.sqrt((art ** 2).mean()))))
    print("  9 noktalı kuadratiğin katsayıları (log M = a logλ² + b logλ + c):"
          " a=%+.4f b=%+.4f c=%+.4f" % tuple(p9c))
    P6["9nokta_kat"] = [float(x) for x in p9c]

    json.dump(dict(defter=D9, nu=nu, P1=P1, P2=P2, P3=P3, P4=P4, P6=P6,
                   lam=[float(x) for x in lam],
                   KALm=[float(x) for x in KALm],
                   WXm=[float(x) for x in WXm]),
              open(SCR173 + "/EGRI.json", "w"), indent=1, default=float)
    print("\n-> %s/EGRI.json" % SCR173)


if __name__ == "__main__":
    main()
