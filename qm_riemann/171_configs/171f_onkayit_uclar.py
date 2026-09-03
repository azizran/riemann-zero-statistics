"""
171f — T2c ÖN KAYIT: λ = 0.50 ve λ = 1.30 (KORELATÖRE BAKILMADAN)
==================================================================
`c`'yi HESAPLAMAZ, J₂/KALİB_u2'ye BAKMAZ. Yalnız yeni gazın MARJİNALLERİ
(σ_ds, σ_X̃, σ_Ĉ, bant-başına W_X) ve model-alan oranları (g_E, g_X —
`Model165.alanlar`, korelatörden ÖNCE hesaplanır) kullanılır. Ölçüm
`171g_olcum.py` ile BU BETİKTEN SONRA alınır.

ÖNGÖRÜ ZİNCİRİ (171'in çerçevesi):
  KALİB_b(yeni) = KALİB_b(Hkeskin) · M        [çarpanlaşma; A ÖZDEŞ ELENİR]
  c_b(yeni)     = KALİB_b(yeni) / W_X_b(yeni) [W_X marjinalden ÖLÇÜLÜR]
  c(yeni)       = geo.ort_b c_b
ve ayrıca çapraz-kontrol olarak analitik A ile:
  KALİB_b = K₀ · exp(−2π²τ_b²σ_X̃(λ=1)²) · M   (T2a'nın mührü)

M ADAYLARI (hepsi 171e'de beş λ-gazında uydurulmuş, PARAMETRE DONDURULDU):
  M0 düz            M = 1                                   (168 §A3)
  M2 kırık kuvvet   M = max(1, (λ_c/λ)^0.1898), λ_c = 1.0109 (inşa doyumu)
  M3 varyans-açığı  M = 1 + 0.2316·(u_c−σ_X̃²)_+/u_c, u_c = 0.057957
  M8 ölçülen çarpan M = (g_E/g₀)(g_X/g₀)²      (θ λ-değişmez varsayımı)
  M1 Gram-doyum     M = 1.2728/(1+0.2728λ)                  (171e'de ÖLDÜ)
  M7 merdiven kesri M = (f_L/f_L⁰)^{−0.1965}                (171e'de ÖLDÜ)
  M9 log-λ kuadratik                                        (EMPİRİK ÇIPA)

ÖN-MÜHÜR (bu betiği koşmadan önce, 4 Eylül 2026 00:00):
 R1 λ=1.30 HAKEMDİR: M2/M3 **tam 1.0000** der (λ > λ_c ⇒ düz);
    M1 **0.940**, M7 **0.967**, M9 ~**1.005** der. Ayrışma %3–6, ölçüm
    hatası ±%0.6 ⇒ **5–10σ'lık bir ayrım.**
 R2 λ=0.50: M2 **1.143**, M3 **1.133**, M9 **1.19** (kuadratik uzatma),
    M0 1.000. M2/M3 ayrımı yalnız %1 (≈1.4σ) — bu uçta aileler ayrışmaz;
    asıl işlevi M'nin yükselen kolunun DEVAM edip etmediğidir.
 R3 c öngörüleri (W_X marjinalden): c(1.30) ≈ 0.455–0.465 (M=1 ile),
    c(0.50) ≈ 0.360–0.372. Yani `c` λ ile artmaya devam eder ve 4/π²'nin
    ÇOK üstüne çıkar (λ=1.30'da +%12-15) — 170 §K2′'nin hükmünün
    örneklem-dışı uzantısı.
 R4 Marjinaller: σ_X̃(0.50) ≈ 0.156, σ_X̃(1.30) ≈ 0.271 (171b ön-mührü).

Kullanım: 171f_onkayit_uclar.py <L050|L130>
Çıktı:    scratchpad/171/ONKAYIT_<gaz>.json

SONUÇ (gerçek koşudan — ÖNGÖRÜLER; hüküm 171h'de):
 L050 (00:00:05): σ_ds = 0.27903, σ_X̃ = 0.15576 ✓ (R4), σ_Ĉ = 0.18131,
   g_E = 0.6470, g_X = 0.7113.
   M/c öngörüleri: M0 1.0000/0.3098 | M2 1.1430/0.3541 | M3 1.1347/0.3515 |
   M8 1.1309/0.3504 | M1 1.1200/0.3470 | M7 1.1614/0.3598 | M9 1.2238/0.3791.
   R3 ISKA: elle yazılmış c(0.50) ≈ 0.360–0.372 aralığı, çarpanlaşmanın
   kendi öngörüsüyle (M2 → 0.3541) uyuşmuyor.
 L130 (00:03:56): σ_ds = 0.48890, σ_X̃ = 0.27713 (R4 ISKA, +%2.3),
   σ_Ĉ = 0.31923, g_E = 0.5877, g_X = 0.7046.
   M/c: M0/M2/M3 1.0000/0.4618 | M8 1.0078/0.4654 | M1 0.9396/0.4339 |
   M7 0.9674/0.4467 | M9 1.0175/0.4699.  R3 ✓ (0.455–0.465 aralığı 0.4618).
 R1 ✓ ayrışma yapısı (M2/M3 tam 1.0000 ↔ M1 0.940 ↔ M7 0.967 ↔ M9 1.018).
"""
import importlib
import json
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
LO_MIN, LO_MAX = 0.52, 0.68
C_HIP = 4.0 / np.pi ** 2
LAM_C = 1.9147 / 1.894
YENI = {"L050": 0.50, "L130": 1.30}

# --- 171e'de dondurulmuş parametreler --------------------------------
P_M2, K_M3, C_M1, Q_M7 = 0.1898, 0.2316, 0.2728, 0.1965
LAM5 = ["L115", "Hkeskin", "L085", "L070", "L060"]


def main(ad):
    t0 = time.time()
    lam = YENI[ad]
    ORT.KUNYE[ad] = dict(gercek=False, lam=lam, tau_ust=1.00, pen=None,
                         aile="lam", T=1.0)
    K.PENCERE.update(ORT.pencere_dict())
    D = K1.yukle()

    # ---- yeni gazın MARJİNALLERİ --------------------------------------
    Y = K.gaz(ad, 0.40, 4000)
    Mo = K.Model165(ad, 0.40, 4000, 0.95, "olculen", Y=Y, tau_c=0.95)
    Mo.sec("olculen", 0.95).alanlar(kmax=3)
    A_ = Mo.artik
    sX = float(np.std(Y.Xtil0))
    sds = float(np.std(Y.ds))
    print("=== 171 T2c ÖN KAYIT — %s (λ = %.2f) ===" % (ad, lam))
    print("  N=%d L=%.5f  σ_ds=%.5f  σ_X̃=%.5f  σ_Ĉ=%.5f"
          % (len(Y.m0), Y.L, sds, sX, Mo.sigC), flush=True)
    print("  g_E=%.4f  g_X=%.4f  g_cal=%.4f  çizgi=%d"
          % (A_["gE"], A_["gX"], A_["gE"] * A_["gX"] ** 2, A_["nline"]),
          flush=True)
    print("  [171b ön-mührü: σ_ds(0.50)∈[0.279,0.285] σ_X̃∈[0.1555,0.1590]; "
          "σ_ds(1.30)∈[0.477,0.487] σ_X̃∈[0.2685,0.2742]]")

    # ---- bant-başına W_X (marjinal) -----------------------------------
    ban = C163.bant_adaylari(Y, K.IZGARA_T1)
    hedef = [b for b in ban if LO_MIN - 1e-9 <= b["lo"] <= LO_MAX + 1e-9]
    tauc = np.array([r["tau"] for b in hedef for r in b["cizgi"]])
    WXc = B166.karakteristik(Y.Xtil0, -TWO_PI * tauc).real
    WX, i = [], 0
    for b in hedef:
        nb = len(b["cizgi"])
        gp = np.array([r["gp"] for r in b["cizgi"]])
        WX.append(float((WXc[i:i + nb] * (gp / gp.sum())).sum()))
        i += nb
    WX = np.array(WX)

    # ---- Hkeskin çapası (bant bant) -----------------------------------
    RH = [K1.band_jk(b["cizgi"]) for b in K1.saglikli(D["Hkeskin"], 0.68)]
    KH = np.array([r["KALIB_u2"] for r in RH])
    TH = np.array([r["tau_eff"] for r in RH])
    sX_hk = D["Hkeskin"]["sigX"]
    gE0, gX0 = D["Hkeskin"]["artik"]["gE"], D["Hkeskin"]["artik"]["gX"]
    Aan = np.exp(-2 * np.pi ** 2 * TH ** 2 * sX_hk ** 2)
    K0 = float(np.exp(np.mean(np.log(KH / Aan))))

    # ---- M adaylarının öngörüleri -------------------------------------
    lam5 = np.array([D[g]["lam"] for g in LAM5], float)
    sX5 = np.array([D[g]["sigX"] for g in LAM5])
    AX, BX = np.polyfit(lam5 ** 2, sX5 ** 2, 1)
    uc = float(AX * LAM_C ** 2 + BX)
    fL = lambda l: AX * l ** 2 / (AX * l ** 2 + BX)          # noqa: E731
    MO = json.load(open(SCR / "171/M_ONKAYIT.json"))
    Mv = np.array([MO["M"][g] for g in LAM5])
    p9 = np.polyfit(np.log(lam5), np.log(Mv), 2)

    Mad = {
        "M0 düz": 1.0,
        "M2 kırık kuvvet (λ_c=%.4f)" % LAM_C:
            float(max(1.0, (LAM_C / lam) ** P_M2)),
        "M3 varyans-açığı":
            float(1 + K_M3 * max(0.0, (uc - sX ** 2) / uc)),
        "M8 ölçülen (g_E/g₀)(g_X/g₀)²":
            float((A_["gE"] / gE0) * (A_["gX"] / gX0) ** 2),
        "M1 Gram-doyum [171e'de ÖLDÜ]":
            float((1 + C_M1) / (1 + C_M1 * lam)),
        "M7 merdiven kesri [171e'de ÖLDÜ]":
            float((fL(lam) / fL(1.0)) ** (-Q_M7)),
        "M9 log-λ kuadratik [EMPİRİK]":
            float(np.exp(np.polyval(p9, np.log(lam)))),
    }

    print("\n  bant lo  n   W_X(%s) | KALİB_b(Hk)  A_an(τ)  τ_eff" % ad)
    for j, b in enumerate(hedef):
        print("  %.2f   %4d  %.4f     | %.4f      %.4f   %.4f"
              % (b["lo"], len(b["cizgi"]), WX[j], KH[j], Aan[j], TH[j]))

    rows = {}
    print("\n  " + "*" * 76)
    print("  ** ÖN KAYIT — λ = %.2f (ÖLÇÜMDEN ÖNCE YAZILDI) **" % lam)
    print("  %-34s %8s | %-39s | %8s" % ("M adayı", "M", "c_b (5 bant)", "c"))
    for nm, m in Mad.items():
        cb = KH * m / WX                       # çarpanlaşma (A özdeş elenir)
        cg = float(np.exp(np.mean(np.log(cb))))
        cb2 = K0 * Aan * m / WX                # analitik A ile çapraz kontrol
        cg2 = float(np.exp(np.mean(np.log(cb2))))
        rows[nm] = dict(M=m, c_bant=[float(x) for x in cb], c=cg,
                        c_analitikA=cg2)
        print("  %-34s %8.4f | %s | %8.4f" %
              (nm, m, " ".join("%.4f" % x for x in cb), cg))
    print("  %-34s %8s | %-39s | %8.4f"
          % ("(çapa) 4/π²", "—", "", C_HIP))
    print("  %-34s %8s | analitik-A çapraz kontrolü: %s"
          % ("", "", "  ".join("%.4f" % rows[n]["c_analitikA"] for n in rows)))
    print("  " + "*" * 76)

    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    out = dict(zaman=ts, gaz=ad, lam=lam, sigds=sds, sigX=sX,
               sigC=float(Mo.sigC), gE=float(A_["gE"]), gX=float(A_["gX"]),
               nline=int(A_["nline"]), L=float(Y.L), N=int(len(Y.m0)),
               lo=[float(b["lo"]) for b in hedef],
               nline_bant=[len(b["cizgi"]) for b in hedef],
               W_X=[float(x) for x in WX],
               KALIB_Hk=[float(x) for x in KH],
               tau_eff_Hk=[float(x) for x in TH],
               A_analitik=[float(x) for x in Aan], K0=K0,
               uc=uc, AX=float(AX), BX=float(BX), C_HIP=C_HIP,
               ongoru=rows, sure_s=time.time() - t0)
    p = SCR / ("171/ONKAYIT_%s.json" % ad)
    p.write_text(json.dumps(out, indent=1, default=float))
    print("\n-> %s   ZAMAN DAMGASI %s   (%.1f dk)"
          % (p, ts, (time.time() - t0) / 60))


if __name__ == "__main__":
    main(sys.argv[1])
