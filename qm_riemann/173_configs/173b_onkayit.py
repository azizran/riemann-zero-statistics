"""
173b — HAKEM ÖN KAYDI: λ = 0.40 / λ = 1.45 (KORELATÖRE BAKILMADAN)
==================================================================
`c`'yi, KALİB'i, θ'yı HESAPLAMAZ. Yalnız yeni gazın MARJİNALLERİ
(σ_ds, σ_X̃, σ_Ĉ, bant-başına W_X) ve MODEL-ALAN büyüklükleri
(g_E, g_X, Q, ρ, μ̂² — `Model165.alanlar`, korelatörden ÖNCE hesaplanır)
kullanılır. Ölçüm `173c_olcum.py` ile BU BETİKTEN SONRA alınır.
Yapı 171f_onkayit_uclar.py ile birebir; hiçbir ölçüm parçası kopyalanmaz.

──────────────────────────────────────────────────────────────────────
A. 172'NİN MÜHÜRLERİ — DEĞİŞTİRİLMEDEN TAŞINIYOR (`172/G4.json`)
──────────────────────────────────────────────────────────────────────
Bunlar 4 Eylül 13:32'de, gazlar İNŞA EDİLMEDEN yazıldı. Bu betik onları
yeniden TÜRETMEZ, dosyadan okur ve aynen basar:

  λ = 0.40  çarpan yolu : g_E 0.68415  g_X 0.71754  θ 0.94525
                          M 1.29469  c 0.37905
            doğrudan yol (171 M9 çıpası, 7 gazlı log-λ kuadratiği):
                          M 1.38809  c 0.40516
  λ = 1.45  çarpan yolu : g_E 0.59391  g_X 0.70852  θ 0.84045
                          M 0.97435  c 0.47471
            doğrudan yol: M 1.00599  c 0.49098

İki yol λ = 0.40'ta M'de **%6.7**, λ = 1.45'te **%3.2** ayrışıyor.
172 kendi tercihini de mühürledi: "çarpan yolu λ = 0.50'yi −%4.4
ıskaladığı için DOĞRUDAN YOL daha güvenilirdir."

──────────────────────────────────────────────────────────────────────
B. 171'İN M AİLESİ — PARAMETRELER 171e'DE DONDURULDU
──────────────────────────────────────────────────────────────────────
  M0 düz            M = 1
  M2 kırık kuvvet   M = max(1, (λ_c/λ)^0.1898),  λ_c = 1.9147/1.894
  M3 varyans-açığı  M = 1 + 0.2316·(u_c − σ_X̃²)_+/u_c
  M8 ölçülen çarpan M = (g_E/g₀)(g_X/g₀)²
  M1 Gram-doyum     M = 1.2728/(1+0.2728λ)              [171e/171h'de ÖLÜ]
  M7 merdiven kesri M = (f_L/f_L⁰)^(−0.1965)            [171e/171h'de ÖLÜ]
  M9 log-λ kuadratik (171'in 5 gazı)                    [EMPİRİK ÇIPA]
  M9′ log-λ kuadratik (172'nin 7 gazı) = G4 "doğrudan yol"

══════════════════════════════════════════════════════════════════════
C. 173'ÜN YENİ ÖN-KAYITLARI (ÖLÇÜMDEN ÖNCE, 4 Eylül 2026 13:50)
══════════════════════════════════════════════════════════════════════
172'nin mührü M, c, g_E, g_X, θ'yı kapsıyor. 173'ün EĞRİ-UÇLARI kapısı
için aşağıdakiler EK ve BAĞIMSIZ ön-kayıtlardır; hepsi yalnız ölçülmüş
yedi λ gazından türetildi, hiçbiri yeni gaza bakmadan yazıldı.

 P1 **ν(λ) UÇLARDA.** ν := Δlog KALİB/Δlog W_X (171 §T2c.3). Ölçülen
    merdiven ν(λ_orta) = 1.4887(0.55) 0.9927(0.65) 0.5627(0.775)
    0.2045(0.925) −0.0337(1.075) +0.2709(1.225). Yerel eğimlerle uzatma:
      **ν(0.50→0.40) ≈ 1.99, ön-kayıt bandı [1.55, 2.45]**  (ν > 1 ⇒
      c artmaya devam eder: **c(0.40) > c(0.50) = 0.3800**)
      **ν(1.30→1.45) ≈ 0.58, ön-kayıt bandı [0.30, 0.90]**  (ν < 1 ⇒
      **c(1.45) > c(1.30) = 0.4540**, İKİNCİ BİR VADİ/TEPE YOK)
    ÖLÜM ÖLÇÜTÜ: bandın dışına çıkarsa ıska yazılır; ν'nün 1'i üst uçta
    KESMESİ (ν(1.30→1.45) > 1) "c'nin ikinci durgun noktası" demektir ve
    bu resmi kökten değiştirir.

 P2 **VADİ/EĞİM-KESİŞME RESMİ UÇLARDA.** dlog c/dλ = dlogKALİB/dλ −
    dlogW_X/dλ. Ölçülen fark: −0.3023(0.55) −0.0010(0.65) +0.2205(0.775)
    +0.3662(0.925) +0.4381(1.075) +0.3397(1.225) — TEK kesişim λ=0.6506.
    Yerel uzatma: **fark(0.45) ≈ −0.60 [−0.85, −0.40]** ve
    **fark(1.375) ≈ +0.24 [+0.08, +0.42]**. ÖN-KAYIT: dokuz noktada da
    **kesişim SAYISI hâlâ 1** ve yeri **λ ∈ [0.63, 0.67]**.

 P3 **α(λ) TEPESİ ve A-ÇARPANLAŞMASININ PENCERESİ.** 5-orta-nokta
    parabolü α(0.40) = **1.144**, α(1.45) = **0.754** der; 7-nokta
    parabolü 0.902 / 0.489 der; L130'dan doğrusal uzatma 0.359 der.
    ÖN-KAYIT (geniş, çünkü üç araç üçe ayrılıyor — dürüstlük):
      **α(0.40) ∈ [0.85, 1.15]** (yani λ=0.50'nin 1.0242'sinden AŞAĞI
      ya da yanında: tümseğin sol kolu düşüyor)
      **α(1.45) ∈ [0.35, 0.75]** (L130'un 0.6874'ünden AŞAĞI)
    ve ayrıca: 171 §T2c-c "çarpanlaşma λ ∈ [0.50, 1.15] penceresinin
    yasasıdır" demişti. ÖN-KAYIT: **λ=1.45'te şekil sapması (S/A2Hk rms)
    L130'un %3.62'sinden BÜYÜK (> %4)**; **λ=0.40'ta ≤ %2.5**
    (L050 %1.50 idi, pencere alt uçta hâlâ geçerli).
    Ek: **σ*/2(1.45) < 0.18661** (L130'unkinin altında).

 P4 **Q_E TÜMSEĞİ.** Tepe λ_c = 1.0109'da (172 §G1.5, ✓✓). Uçlar tepenin
    İKİ YANINDA olduğu için ikisi de düşmeli. Parabolün (5 orta, log Q,
    λ) uzatması **Q_E(0.40) = 0.652**, **Q_E(1.45) = 0.779**; 7-nokta
    parabolü 0.617 / 0.739; log-λ kuadratiği 0.545 / 0.809.
    ÖN-KAYIT (işaret + bant): **Q_E(0.40) < 0.69027** (L050'nin altında)
    bandı **[0.54, 0.68]**; **Q_E(1.45) < 0.83903** (L130'un altında)
    bandı **[0.70, 0.83]**. Tepe dokuz noktayla yeniden uydurulduğunda
    **λ_c'nin ±0.05'inde** kalmalı.

 P5 **θ TAŞIYICISI (172-G2, bonus).** θ/θ₀ = (1−βφ)(ρ_X/ρ_X₀)^0.651;
    λ ailesinde φ = 0. Yedi gazda taşıyıcı artığı %: −2.42 +0.04 +1.13
    +1.10 0.00 −1.55 −1.20 (rms %1.37). ÖN-KAYIT: iki uçta da
    **|artık| ≤ %3 ⇒ taşıyıcı AYAKTA; > %3 ⇒ ÖLDÜ.** Bu betik ρ_X'i
    ÖLÇER (korelatör öncesi) ve θ öngörüsünü ondan üretir; θ'nın kendisi
    173c'de ölçülür.

 P6 **G4'ün ALAN-GİRDİLERİ BU BETİKTE HAKEMLENİR.** g_E, g_X, ρ_X, Q_E
    korelatörden ÖNCE ölçülür. 172'nin çarpan yolu bu dördünü
    ekstrapolasyondan almıştı; bu betik onları GERÇEKTEN ölçer. Yani
    çarpan yolunun M öngörüsünün hatası burada **girdi hatası** ile
    **model hatası**na ayrışır. ÖN-KAYIT: girdi hataları (g_E, g_X)
    **≤ %1.5**, yani λ=0.40'taki %6.7'lik ayrışma girdiden DEĞİL,
    modelden (θ taşıyıcısı + ρ_X üssü) gelmeli.

Kullanım: 173b_onkayit.py <L040|L145|L140>
Çıktı:    scratchpad/173/ONKAYIT_<gaz>.json  (zaman damgalı)

══════════════════════════════════════════════════════════════════════
SONUÇ (yalnız gerçek koşudan; ÖNGÖRÜLER — hüküm 173d/173e'de)
══════════════════════════════════════════════════════════════════════
 L040 (zaman damgası 2026-09-04 14:00:34):
   σ_ds = 0.22728  σ_X̃ = 0.12673  σ_Ĉ = 0.15011
   g_E = 0.68657  g_X = 0.72094  Q_E = 0.56708  ρ_E = 1.80928
   Q_X = 0.47940  ρ_X = 1.71790  μ̂²_E = +0.01337  μ̂²_X = −3.5e−15
   ÖZDEŞLİK |ΔK|/K = 1.6e−15 (E) / 3.3e−15 (X) — 172 §G1.1 mühürleri ✓
   **P6 TUTTU:** G4'ün mührü g_E 0.68415 → ölçüm 0.68657 (+0.35%);
   g_X 0.71754 → 0.72094 (+0.47%);  ρ_X uzatması 1.71842 → 1.71790
   (−0.03%!).  P4 ✓ (Q_E ∈ [0.54, 0.68]).
   M öngörüleri: M0 1.0000 | M2 1.1924 | M3 1.1674 | M8 1.2328 |
   M1 1.1476 | M7 1.2407 | **M9 1.4112** | **M9′ 1.3881** |
   G4 çarpan 1.2947 | M10 1.3113.
 L145 (14:08:29):
   σ_ds = 0.54445  σ_X̃ = 0.31782  σ_Ĉ = 0.37422
   g_E = 0.60919  g_X = 0.69707  Q_E = 0.55489  ρ_E = 1.41984
   Q_X = 0.42012  ρ_X = 1.38687  μ̂²_E = +0.02067
   **P6 ISKA:** g_E +2.57%, g_X −1.62% (ön-kayıt ≤ %1.5);
   ρ_X −3.33%;  **Q_E uzatması −%31.4 ıskaladı** ve ölçüm 0.55489,
   P4 bandı [0.70, 0.83]'ün DIŞINDA ⇒ **P4 ÜST UÇTA ÖLDÜ.**
   M öngörüleri: M0/M2/M3 1.0000 | M8 1.0226 | M1 0.9120 | M7 0.9572 |
   M9 1.0410 | M9′ 1.0060 | G4 çarpan 0.9744 | M10 0.9463.
 L140 (14:27:58, YEDEK):
   σ_ds = 0.52231  σ_X̃ = 0.30150  σ_Ĉ = 0.35174
   g_E = 0.59976  g_X = 0.69817  Q_E = 0.66240  ρ_E = 1.65501
   ρ_X = 1.42625.  **G4 mührü yoktur ve uydurulmadı** (Y4 ✓).
   Q_E ön-kayıt bandı [0.71, 0.84]'ün DIŞINDA (0.66240) — yani Q_E'nin
   çöküşü SAĞLIKLI gazda, sadakat sınırından ÖNCE başlıyor.
"""
import importlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("169_configs", "167_configs", "166_configs", "165_configs",
           "163_configs", "160_configs", "159_configs", "155_configs",
           "154_configs"):
    sys.path.insert(0, str(QM / _p))
K = importlib.import_module("165_cekirdek")
C163 = importlib.import_module("163_cekirdek")
B166 = importlib.import_module("166_bacak")
ORT = importlib.import_module("167_ortak")
K1 = importlib.import_module("169_k1")

TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
SCR173 = SCR / "173"
LO_MIN, LO_MAX = 0.52, 0.68
C_HIP = 4.0 / np.pi ** 2
LAM_C = 1.9147 / 1.894
YENI = {"L040": 0.40, "L145": 1.45, "L140": 1.40}
P_M2, K_M3, C_M1, Q_M7, P_TH = 0.1898, 0.2316, 0.2728, 0.1965, 0.651
LAM5 = ["L115", "Hkeskin", "L085", "L070", "L060"]
LAM7 = ["L050", "L060", "L070", "L085", "Hkeskin", "L115", "L130"]
LV7 = np.array([.50, .60, .70, .85, 1.00, 1.15, 1.30])


def kusur(V_M, V_O, V_R, kor, g):
    """172b'nin (Ö-A..Ö-D) özdeşliği — AYNEN (kopya değil, aynı 5 satır)."""
    Kv = 0.5 * (V_O + V_M - V_R)
    Kk = kor * math.sqrt(V_M * V_O)
    S = Kv / g
    return dict(K=Kv, dK=abs(Kv - Kk) / abs(Kv), S=S, mu2=(S - V_M) / V_O,
                pi=Kv / V_O, rho=S / V_O, Q=(S - Kv) / V_O, r=V_R / V_O)


def main(ad):
    t0 = time.time()
    lam = YENI[ad]
    ORT.KUNYE[ad] = dict(gercek=False, lam=lam, tau_ust=1.00, pen=None,
                         aile="lam", T=1.0)
    K.PENCERE.update(ORT.pencere_dict())
    D = K1.yukle()
    G1 = json.load(open(SCR / "172/G1.json"))
    G2 = json.load(open(SCR / "172/G2.json"))
    G4 = json.load(open(SCR / "172/G4.json"))
    MO = json.load(open(SCR / "171/M_ONKAYIT.json"))
    NU = json.load(open(SCR / "171/NU_MERDIVEN.json"))

    # ---- yeni gazın MARJİNALLERİ ve ALANLARI (korelatör YOK) ----------
    Y = K.gaz(ad, 0.40, 4000)
    Mo = K.Model165(ad, 0.40, 4000, 0.95, "olculen", Y=Y, tau_c=0.95)
    Mo.sec("olculen", 0.95).alanlar(kmax=3)
    A_ = Mo.artik
    sX = float(np.std(Y.Xtil0))
    sds = float(np.std(Y.ds))
    E_ = kusur(A_["varE_mod"], A_["varE_olc"], A_["varE_res"], A_["korE"],
               A_["gE"])
    X_ = kusur(A_["varX_mod"], A_["varX_olc"], A_["varX_res"], A_["korX"],
               A_["gX"])
    print("=" * 100)
    print("173b — HAKEM ÖN KAYDI: %s (λ = %.2f)   [KORELATÖRE BAKILMADI]"
          % (ad, lam))
    print("=" * 100)
    print("  N=%d L=%.5f  σ_ds=%.5f  σ_X̃=%.5f  σ_Ĉ=%.5f  çizgi=%d"
          % (len(Y.m0), Y.L, sds, sX, Mo.sigC, A_["nline"]), flush=True)
    print("  g_E=%.5f  g_X=%.5f  g_cal=%.5f"
          % (A_["gE"], A_["gX"], A_["gE"] * A_["gX"] ** 2))
    print("  ÖZDEŞLİK: |ΔK|/K = %.2e (E) / %.2e (X)   μ̂²_E=%+.5f  "
          "μ̂²_X=%+.2e" % (E_["dK"], X_["dK"], E_["mu2"], X_["mu2"]))
    print("  Q_E=%.5f  ρ_E=%.5f  rE=%.5f | Q_X=%.5f  ρ_X=%.5f  rX=%.5f"
          % (E_["Q"], E_["rho"], E_["r"], X_["Q"], X_["rho"], X_["r"]))

    # ---- P6: G4'ün ALAN-GİRDİLERİNİN HAKEMİ ---------------------------
    i = 0 if lam < 1 else 1
    print("\n  ** P6 — 172/G4 ÇARPAN YOLUNUN ALAN GİRDİLERİ (mühür ↔ ölçüm)")
    print("  %-10s %10s %10s %9s" % ("büyüklük", "G4 mührü", "ÖLÇÜLEN",
                                     "fark%"))
    p6 = {}
    if ad in ("L040", "L145"):
        for nm, mu, ol in (("g_E", G4["ongoru"]["gE"][i], A_["gE"]),
                           ("g_X", G4["ongoru"]["gX"][i], A_["gX"])):
            p6[nm] = dict(muhur=mu, olcum=float(ol),
                          fark=100 * (ol / mu - 1))
            print("  %-10s %10.5f %10.5f %+8.2f%%" % (nm, mu, ol,
                                                      100 * (ol / mu - 1)))
    # ρ_X ve Q_E ekstrapolasyonu (172f'nin kuadratiği ile aynı)
    def kuad(y, x):
        c = np.polyfit(np.log(LV7), np.log(y), 2)
        return float(np.exp(np.polyval(c, np.log(x))))
    rhoX7 = np.array([G1[g]["X"]["rho"] for g in LAM7])
    QE7 = np.array([G1[g]["E"]["Q"] for g in LAM7])
    for nm, ext, ol in (("ρ_X", kuad(rhoX7, lam), X_["rho"]),
                        ("Q_E", kuad(QE7, lam), E_["Q"])):
        p6[nm] = dict(muhur=ext, olcum=float(ol), fark=100 * (ol / ext - 1))
        print("  %-10s %10.5f %10.5f %+8.2f%%   (172f uzatması)"
              % (nm, ext, ol, 100 * (ol / ext - 1)))

    # ---- P4/P5: Q_E bandı ve θ taşıyıcısı ------------------------------
    th0, rhoX0 = G2["theta"]["Hkeskin"], G1["Hkeskin"]["X"]["rho"]
    th_pred = th0 * (X_["rho"] / rhoX0) ** P_TH      # ÖLÇÜLEN ρ_X ile
    print("\n  ** P5 — θ TAŞIYICISI (ölçülen ρ_X ile): θ öngörüsü = %.5f"
          "   [G4 mührü %.5f]" % (th_pred, G4["ongoru"]["theta"][i]
                                  if ad in ("L040", "L145") else float("nan")))
    QEband = {"L040": (0.54, 0.68), "L145": (0.70, 0.83),
              "L140": (0.71, 0.84)}[ad]
    print("  ** P4 — Q_E ölçülen %.5f;  ön-kayıt bandı [%.2f, %.2f] ⇒ %s"
          % (E_["Q"], QEband[0], QEband[1],
             "İÇERİDE" if QEband[0] <= E_["Q"] <= QEband[1] else "DIŞARIDA"))

    # ---- bant-başına W_X (marjinal) -----------------------------------
    ban = C163.bant_adaylari(Y, K.IZGARA_T1)
    hedef = [b for b in ban if LO_MIN - 1e-9 <= b["lo"] <= LO_MAX + 1e-9]
    tauc = np.array([r["tau"] for b in hedef for r in b["cizgi"]])
    WXc = B166.karakteristik(Y.Xtil0, -TWO_PI * tauc).real
    WX, j = [], 0
    for b in hedef:
        nb = len(b["cizgi"])
        gp = np.array([r["gp"] for r in b["cizgi"]])
        WX.append(float((WXc[j:j + nb] * (gp / gp.sum())).sum()))
        j += nb
    WX = np.array(WX)

    # ---- Hkeskin çapası ------------------------------------------------
    RH = [K1.band_jk(b["cizgi"]) for b in K1.saglikli(D["Hkeskin"], 0.68)]
    KH = np.array([r["KALIB_u2"] for r in RH])
    TH = np.array([r["tau_eff"] for r in RH])
    gE0, gX0 = D["Hkeskin"]["artik"]["gE"], D["Hkeskin"]["artik"]["gX"]
    Aan = np.exp(-2 * np.pi ** 2 * TH ** 2 * D["Hkeskin"]["sigX"] ** 2)
    K0 = float(np.exp(np.mean(np.log(KH / Aan))))

    # ---- M adayları ----------------------------------------------------
    lam5 = np.array([D[g]["lam"] for g in LAM5], float)
    sX5 = np.array([D[g]["sigX"] for g in LAM5])
    AX, BX = np.polyfit(lam5 ** 2, sX5 ** 2, 1)
    uc = float(AX * LAM_C ** 2 + BX)
    fL = lambda l: AX * l ** 2 / (AX * l ** 2 + BX)          # noqa: E731
    Mv5 = np.array([MO["M"][g] for g in LAM5])
    p9 = np.polyfit(np.log(lam5), np.log(Mv5), 2)
    MOl7 = np.array(NU["M"])[::-1]
    p9b = np.polyfit(np.log(LV7), np.log(MOl7), 2)

    Mad = {
        "M0 düz": 1.0,
        "M2 kırık kuvvet (λ_c=%.4f)" % LAM_C:
            float(max(1.0, (LAM_C / lam) ** P_M2)),
        "M3 varyans-açığı": float(1 + K_M3 * max(0.0, (uc - sX ** 2) / uc)),
        "M8 ölçülen (g_E/g₀)(g_X/g₀)²":
            float((A_["gE"] / gE0) * (A_["gX"] / gX0) ** 2),
        "M1 Gram-doyum [ÖLÜ]": float((1 + C_M1) / (1 + C_M1 * lam)),
        "M7 merdiven kesri [ÖLÜ]": float((fL(lam) / fL(1.0)) ** (-Q_M7)),
        "M9 log-λ kuad (171, 5 gaz) [EMPİRİK]":
            float(np.exp(np.polyval(p9, np.log(lam)))),
        "M9′ log-λ kuad (172, 7 gaz) = G4 DOĞRUDAN":
            float(np.exp(np.polyval(p9b, np.log(lam)))),
        "G4 ÇARPAN YOLU [172 mührü]":
            float(G4["ongoru"]["M"][i]) if ad in ("L040", "L145")
            else float("nan"),
        "M10 ölçülen çarpan + θ-taşıyıcı":
            float((A_["gE"] / gE0) * (A_["gX"] / gX0) ** 2 * (th_pred / th0)),
    }

    print("\n  bant lo   n     W_X(%s) | KALİB_b(Hk)  A_an(τ)  τ_eff" % ad)
    for j, b in enumerate(hedef):
        print("  %.2f    %4d   %.5f  | %.5f      %.5f  %.5f"
              % (b["lo"], len(b["cizgi"]), WX[j], KH[j], Aan[j], TH[j]))

    rows = {}
    print("\n  " + "*" * 88)
    print("  ** ÖN KAYIT — λ = %.2f (ÖLÇÜMDEN ÖNCE YAZILDI) **" % lam)
    print("  %-38s %8s | %-41s | %8s"
          % ("M adayı", "M", "c_b (5 bant)", "c"))
    for nm, m in Mad.items():
        if m != m:
            continue
        cb = KH * m / WX
        cg = float(np.exp(np.mean(np.log(cb))))
        cb2 = K0 * Aan * m / WX
        rows[nm] = dict(M=m, c_bant=[float(x) for x in cb], c=cg,
                        c_analitikA=float(np.exp(np.mean(np.log(cb2)))))
        print("  %-38s %8.4f | %s | %8.4f"
              % (nm, m, " ".join("%.4f" % x for x in cb), cg))
    if ad in ("L040", "L145"):
        print("  %-38s %8.4f | %-41s | %8.4f"
              % ("G4 ÇARPAN YOLU c [172 mührü]", G4["ongoru"]["M"][i], "",
                 G4["ongoru"]["c"][i]))
        print("  %-38s %8.4f | %-41s | %8.4f"
              % ("G4 DOĞRUDAN YOL c [172 mührü]",
                 G4["ongoru"]["M_dogrudan"][i], "",
                 G4["ongoru"]["c_dogrudan"][i]))
    print("  %-38s %8s | %-41s | %8.4f" % ("(çapa) 4/π²", "—", "", C_HIP))
    print("  " + "*" * 88)

    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    SCR173.mkdir(parents=True, exist_ok=True)
    out = dict(zaman=ts, gaz=ad, lam=lam, sigds=sds, sigX=sX,
               sigC=float(Mo.sigC), gE=float(A_["gE"]), gX=float(A_["gX"]),
               nline=int(A_["nline"]), L=float(Y.L), N=int(len(Y.m0)),
               alan_E={k: float(v) for k, v in E_.items()},
               alan_X={k: float(v) for k, v in X_.items()},
               artik={k: float(v) for k, v in A_.items()
                      if isinstance(v, (int, float))},
               theta_tasiyici=float(th_pred), p6=p6,
               lo=[float(b["lo"]) for b in hedef],
               W_X=[float(x) for x in WX],
               KALIB_Hk=[float(x) for x in KH],
               tau_eff_Hk=[float(x) for x in TH],
               A_analitik=[float(x) for x in Aan], K0=K0, uc=uc,
               C_HIP=C_HIP, ongoru=rows,
               G4_muhur=(dict(gE=G4["ongoru"]["gE"][i],
                              gX=G4["ongoru"]["gX"][i],
                              theta=G4["ongoru"]["theta"][i],
                              M=G4["ongoru"]["M"][i], c=G4["ongoru"]["c"][i],
                              M_dogrudan=G4["ongoru"]["M_dogrudan"][i],
                              c_dogrudan=G4["ongoru"]["c_dogrudan"][i])
                         if ad in ("L040", "L145") else None),
               sure_s=time.time() - t0)
    p = SCR173 / ("ONKAYIT_%s.json" % ad)
    p.write_text(json.dumps(out, indent=1, default=float))
    print("\n-> %s   ZAMAN DAMGASI %s   (%.1f dk)"
          % (p, ts, (time.time() - t0) / 60))


if __name__ == "__main__":
    main(sys.argv[1])
