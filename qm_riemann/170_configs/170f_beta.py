"""
170 — BONUS: (i) T'nin DERİN-SARILMA açılımı = sonlu-doyum düzeltme yasası
             (ii) β = 0.2175 bu açılımdan çıkıyor mu?
             (iii) T/4 artığı (%8–11) aynı çerçevede eriyor mu?
=========================================================================
Yeni ölçüm YOK: `scratchpad/170/K1_T.json` + `K1_T_kesim.json` (170b'nin
bacak ölçümleri) ve `scratchpad/170/K0.json` (170a'nın muhasebesi) okunur.
T yasası `170b_T_yasasi` modülünden import edilir.

(i) DERİN-SARILMA AÇILIMI (kalemle). ρ = 1−ε, ε = e^{−σ²}:
      arcsin(1−ε) = π/2 − √(2ε) + O(ε^{3/2})
      T² = (1−ε)/(π/2−√(2ε)) = (2/π)[1 + (2√2/π)√ε − ε + O(ε^{3/2})]
    ⇒ **T(σ) = √(2/π) · [ 1 + (√2/π) e^{−σ²/2} + O(e^{−σ²}) ]**        (T5)
    Yani SONLU-DOYUM DÜZELTMESİ `e^{−σ²/2}` — tam olarak sarılmış
    Gauss'un KOHERENT genliği — ile doğrusaldır, katsayısı √2/π = 0.45016.
    Dört bacak için:
      c ≈ (4/π²)·[1 + (√2/π) Σ_b e^{−σ_b²/2}]                          (T6)

(ii) β SINAVI. Kesim ekseninde ölçülen yasa θ/θ₀ = 1 − 0.2175φ. (T6)
     kesim gazlarında ne veriyor? (σ_Ĉ ve τ_b kesimle neredeyse
     değişmiyorsa, (T6) φ'ye DUYARSIZ kalır ve β'yı ÜRETEMEZ.)

(iii) T/4 ARTIĞI. Pencere ekseninde σ_Ĉ pencereyle değişmiyorsa (T6)
      yine duyarsızdır; artığın `g_E` ve `g_X` arasındaki paylaşımı
      ayrıca ölçülür (yeni bilgi).

ÖN-MÜHÜR (koşudan ÖNCE):
  * (T5) σ ≥ 1.5'te T'yi %1 içinde vermeli; σ = 1.23'te %2–4 hata.
  * K090'ın ΠT'si Hkeskin'inkinden ‰3'ten az farklı çıkacak
    (σ_Ĉ = 0.27872 ↔ 0.27768), oysa θ/θ₀ = 0.9080 ⇒ **β TÜRETİLEMEZ**,
    (T6)'nın kesim öngörüsü ≈ 1.000 olacak (ölçülen 0.79–0.91).
  * T/4'te σ_Ĉ = 0.2720/0.2724 (Hkeskin 0.2777, −%2) ⇒ (T6) ≈ +%0.5,
    ölçülen artık +%8…+%11 ⇒ **T/4 artığı da erimeyecek.**

Çıktı: scratchpad/170/BONUS.json
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "170_configs"))
TB = importlib.import_module("170b_T_yasasi")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/170")
TWO_PI = 2 * np.pi
T_INF = np.sqrt(2 / np.pi)
C_HIP = 4.0 / np.pi ** 2
BETA = 0.2175


def T_derin(sig):
    """(T5) derin-sarılma açılımı."""
    return T_INF * (1 + (np.sqrt(2) / np.pi) * np.exp(-np.asarray(sig, float) ** 2 / 2))


def main():
    K0 = json.load(open(SCR / "K0.json"))["gaz"]
    KT = json.load(open(SCR / "K1_T.json"))["gaz"]
    p = SCR / "K1_T_kesim.json"
    if p.exists():
        KT.update(json.load(open(p))["gaz"])

    print("=" * 100)
    print("BONUS (i) — (T5): SONLU-DOYUM DÜZELTME YASASININ KAPALI BİÇİMİ")
    print("=" * 100)
    print("    T(σ) = √(2/π)·[1 + (√2/π) e^{−σ²/2} + O(e^{−σ²})],   "
          f"√2/π = {np.sqrt(2)/np.pi:.5f}")
    print("\n     σ      T(σ) TAM     (T5) açılım     fark")
    exp_tab = []
    for s in (0.8, 1.0, 1.25, 1.5, 2.0, 3.0, 4.8):
        tt, ta = float(TB.T_sigma(s)), float(T_derin(s))
        print(f"   {s:5.2f}   {tt:.6f}     {ta:.6f}     {100*(ta/tt-1):+7.3f}%")
        exp_tab.append(dict(sigma=s, T=tt, T5=ta))

    # ---------- (ii) β ----------------------------------------------
    print("\n" + "=" * 100)
    print("BONUS (ii) — β = 0.2175 (T)'DEN TÜRÜYOR MU?  [kesim ekseni]")
    print("=" * 100)
    h = KT.get("Hkeskin")
    piT0 = h["PiT"]
    th0 = K0["Hkeskin"]["theta"]
    print(f"  {'gaz':9s} {'φ':>7s} {'σ_Ĉ':>8s} {'σ_Ĉ/σ₀':>7s} "
          f"{'Π_b T(σ_b)':>11s} {'ΠT/ΠT₀':>8s} | {'θ/θ₀ ÖLÇÜLEN':>13s} "
          f"{'1−0.2175φ':>10s}")
    beta_tab = []
    for ad in ("Hkeskin", "K090", "K070", "HA4", "E060"):
        if ad not in KT or ad not in K0:
            continue
        r, k = KT[ad], K0[ad]
        rat = r["PiT"] / piT0
        tr = k["theta"] / th0
        beta_tab.append(dict(gaz=ad, phi=k["phi"], sigC=r["sigC"],
                             PiT=r["PiT"], PiT_rat=rat, theta_rat=tr,
                             yasa=1 - BETA * k["phi"]))
        print(f"  {ad:9s} {k['phi']:7.4f} {r['sigC']:8.5f} "
              f"{r['sigC']/h['sigC']:7.4f} {r['PiT']:11.5f} {rat:8.4f} | "
              f"{tr:13.4f} {1-BETA*k['phi']:10.4f}")
    if len(beta_tab) > 1:
        ph = np.array([b["phi"] for b in beta_tab])
        pt = np.array([b["PiT_rat"] for b in beta_tab])
        bT = float(-np.polyfit(ph, pt, 1)[0]) if len(ph) > 2 else float("nan")
        print(f"\n  (T6)'nın ürettiği eğim:  β_T = {bT:+.5f}   "
              f"(ölçülen β = {BETA:.4f})   oran = {bT/BETA:+.3f}")

    # ---------- (iii) T/4 artığı --------------------------------------
    print("\n" + "=" * 100)
    print("BONUS (iii) — T/4 ARTIĞI (%8–11): (T6) ERİTİYOR MU? + g_E/g_X "
          "PAYLAŞIMI")
    print("=" * 100)
    g0 = K0["Hkeskin"]
    n0 = g0["nline"] / g0["N"]
    print(f"  {'gaz':7s} {'N':>7s} {'σ_Ĉ':>8s} {'σ_Ĉ/σ₀':>7s} | "
          f"{'g_E/g₀':>7s} {'g_X/g₀':>7s} {'(g_X/g₀)²':>9s} "
          f"{'g_cal/g₀':>9s} {'(1+2n/N)₀/(1+2n/N)':>19s} {'artık':>7s} | "
          f"{'θ/θ₀':>7s}")
    t4 = []
    for ad in ("Hkeskin", "HkT2a", "HkT2b", "HkT4a", "HkT4b"):
        k = K0[ad]
        x = k["nline"] / k["N"]
        lk = (1 + 2 * n0) / (1 + 2 * x)
        gr = k["gcal"] / g0["gcal"]
        t4.append(dict(gaz=ad, N=k["N"], sigC=k["sigC"], gE=k["gE"] / g0["gE"],
                       gX=k["gX"] / g0["gX"], gcal=gr, leak=lk, artik=gr / lk,
                       theta=k["theta"] / th0))
        print(f"  {ad:7s} {k['N']:7d} {k['sigC']:8.5f} "
              f"{k['sigC']/g0['sigC']:7.4f} | {k['gE']/g0['gE']:7.4f} "
              f"{k['gX']/g0['gX']:7.4f} {(k['gX']/g0['gX'])**2:9.4f} "
              f"{gr:9.4f} {lk:19.4f} {gr/lk:7.4f} | "
              f"{k['theta']/th0:7.4f}")
    print("\n  (T6)'nın pencere ekseni öngörüsü: σ_Ĉ pencereyle ‰20 içinde "
          "sabit ⇒ ΠT sabit ⇒ T/4 artığına KATKI YOK.")

    json.dump(dict(T5=exp_tab, beta=beta_tab, T4=t4, BETA=BETA),
              open(SCR / "BONUS.json", "w"), indent=1, default=float)
    print(f"\n-> {SCR}/BONUS.json")


if __name__ == "__main__":
    main()
