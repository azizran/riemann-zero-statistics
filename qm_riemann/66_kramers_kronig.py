"""
66 — KRAMERS-KRONIG TESTİ: w VE v TEK KOMPLEKS YANITIN İKİ YÜZÜ MÜ?
=====================================================================
(18 Ağustos 2026 — Sıçrama 1)

Fikir: sıfır gazı bir "malzeme", asal dalgası ona giren ışık.
Nedensel lineer yanıt χ(τ) = χ'(τ) + i·χ''(τ) için:
  - χ' ÇİFT fonksiyon; χ'' TEK fonksiyon (χ(−τ) = χ*(τ))
  - Bizde: w(0)≈sonlu (çift ✓), v(0)=0, v≈2τ (tek ✓) — parite UYUYOR
  - χ' ile χ'' Kramers-Kronig ile bağlı: biri ötekini belirler

Model: Lorentz osilatör(ler)i — nedensel malzemenin standart formu:
  χ(τ) = ε + Σ_k A_k / (τ_k² − τ² − i·γ_k·τ)
  w = Re χ,  v = λ·|Im χ|  (λ: iki kanalın kalibrasyon oranı — TEK serbestlik)

ASIL TEST (öngörü): Re-formu YALNIZ w verisine fit et (v'yi hiç görmeden),
sonra Im'i λ'yla ölçekleyip v'yle kıyasla. Şekil tutarsa (başlangıç eğimi
~2.01, doyum ~0.55 @ τ≈0.4) → w ve v tek analitik fonksiyonun iki yüzü.

Kontrol: ortak-fit cezası (AIC) — bağ ne kadara mal oluyor?
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent

# ---------------------------------------------------------------
# VERİ (56/47/54/58/62 tablolarından; üretim scriptleri repoda)
# ---------------------------------------------------------------
W_DATA = np.array([  # (tau, w, se)
[0.0283,0.8663,0.0079],[0.0295,0.8628,0.0054],[0.0327,0.8510,0.0046],
[0.0367,0.8321,0.0038],[0.0418,0.8201,0.0035],[0.0449,0.7964,0.0097],
[0.0468,0.7956,0.0067],[0.0519,0.7730,0.0056],[0.0557,0.7745,0.0024],
[0.0579,0.7710,0.0024],[0.0582,0.7527,0.0047],[0.0604,0.7635,0.0023],
[0.0634,0.7552,0.0022],[0.0658,0.7084,0.0126],[0.0663,0.7271,0.0044],
[0.0669,0.7475,0.0022],[0.0685,0.7177,0.0109],[0.0703,0.7350,0.0021],
[0.0760,0.6814,0.0072],[0.0795,0.6643,0.0149],[0.0828,0.6580,0.0129],
[0.0852,0.6509,0.0061],[0.0883,0.6626,0.0030],[0.0917,0.6526,0.0030],
[0.0918,0.6169,0.0086],[0.0958,0.6448,0.0029],[0.0971,0.6193,0.0056],
[0.1005,0.6316,0.0028],[0.1030,0.5870,0.0072],[0.1060,0.6202,0.0027],
[0.1115,0.6037,0.0026],[0.1173,0.5424,0.0066],[0.1293,0.5260,0.0039],
[0.1344,0.5080,0.0039],[0.1404,0.4949,0.0037],[0.1473,0.4778,0.0036],
[0.1552,0.4590,0.0035],[0.1563,0.4321,0.0046],[0.1625,0.4191,0.0045],
[0.1633,0.4404,0.0033],[0.1697,0.4039,0.0044],[0.1781,0.3876,0.0042],
[0.1877,0.3616,0.0041],[0.1974,0.3403,0.0039]])

PRIMES = [2, 3, 5, 7, 11, 13]
V_ROWS = {  # L: ([v'ler], se)
    5.60: ([0.2155,0.3414,0.4662,0.5152,0.5690,0.5601], 0.025),
    6.30: ([0.2058,0.3156,0.4399,0.5099,0.5703,0.5902], 0.0166),
    6.99: ([0.1868,0.2978,0.4200,0.4879,0.5630,0.5853], 0.0111),
    7.69: ([0.1756,0.2784,0.3963,0.4650,0.5450,0.5728], 0.0075),
    8.39: ([0.1632,0.2590,0.3714,0.4382,0.5231,0.5494], 0.0051),
    9.08: ([0.1517,0.2420,0.3478,0.4133,0.4953,0.5233], 0.0034),
    9.86: ([0.1405,0.2245,0.3254,0.3877,0.4663,0.4928], 0.004),
    10.37:([0.1358,0.2159,0.3108,0.3713,0.4473,0.4744], 0.004),
    10.93:([0.1279,0.2034,0.2962,0.3576,0.4299,0.4581], 0.004),
    11.47:([0.1222,0.1950,0.2823,0.3403,0.4172,0.4404], 0.004),
    11.98:([0.1174,0.1847,0.2701,0.3258,0.3988,0.4220], 0.004),
    12.45:([0.1114,0.1787,0.2628,0.3123,0.3807,0.4082], 0.004),
    24.48:([0.054,0.091,0.126,0.163,0.202,0.202], 0.008),
}
V_LIST = [(np.log(p)/L, v, se) for L,(vs,se) in V_ROWS.items()
          for p, v in zip(PRIMES, vs)]
V_LIST += [(0.0155,0.0346,0.0097),(0.0246,0.0485,0.0118),(0.0361,0.0716,0.0153),
           (0.0437,0.0941,0.0181),(0.0538,0.1109,0.0226),(0.0575,0.1226,0.0246),
           (0.0148,0.0297,0.0095),(0.0235,0.0442,0.0117),(0.0344,0.0698,0.0151),
           (0.0416,0.0872,0.0179),(0.0512,0.1034,0.0224),(0.0548,0.1048,0.0243)]
V_DATA = np.array(sorted(V_LIST))
tw, w, sw = W_DATA.T
tv, v, sv = V_DATA.T
print(f"Veri: {len(tw)} w-noktası (τ: {tw.min():.3f}-{tw.max():.3f}), "
      f"{len(tv)} v-noktası (τ: {tv.min():.3f}-{tv.max():.3f})")

# ---------------------------------------------------------------
# MODEL
# ---------------------------------------------------------------
def chi(tau, params, K):
    eps = params[0]
    z = np.full_like(tau, eps, dtype=complex)
    for k in range(K):
        A, t0, g = params[1 + 3*k : 4 + 3*k]
        z += A / (t0**2 - tau**2 - 1j * g * tau)
    return z

def fit_w_only(K, tries=40):
    rng = np.random.default_rng(66)
    best = None
    for _ in range(tries):
        p0 = [rng.uniform(-0.5, 0.5)]
        for _k in range(K):
            p0 += [rng.uniform(0.05, 0.6), rng.uniform(0.2, 1.2), rng.uniform(0.05, 1.0)]
        try:
            r = least_squares(
                lambda p: (w - chi(tw, p, K).real) / sw, p0,
                bounds=([-2] + [0.0, 0.05, 0.005] * K,
                        [2] + [5.0, 3.0, 3.0] * K), max_nfev=20000)
            if best is None or r.cost < best.cost:
                best = r
        except Exception:
            pass
    return best

for K in (1, 2):
    fw = fit_w_only(K)
    chi2_w = 2 * fw.cost
    dof_w = len(tw) - len(fw.x)
    print(f"\n=== K={K} osilatör ===")
    print(f"w-yalnız fit: χ²/dof = {chi2_w/dof_w:.2f}")
    for k in range(K):
        A, t0, g = fw.x[1+3*k:4+3*k]
        print(f"  osilatör {k+1}: A={A:.3f}, τ₀={t0:.3f}, γ={g:.3f}")
    print(f"  ε = {fw.x[0]:+.3f},  Re χ(0) = {chi(np.array([0.0]), fw.x, K).real[0]:.3f}")

    # ÖNGÖRÜ: Im şeklini v'yle kıyasla — tek serbestlik λ
    imv = np.abs(chi(tv, fw.x, K).imag)
    lam = np.sum(v * imv / sv**2) / np.sum(imv**2 / sv**2)
    chi2_v = np.sum(((v - lam * imv) / sv) ** 2)
    print(f"  ÖNGÖRÜ (v'yi görmeden): λ = {lam:.3f} → v-χ²/dof = "
          f"{chi2_v/(len(tv)-1):.2f}")
    tt = np.linspace(1e-3, 0.5, 400)
    imf = np.abs(chi(tt, fw.x, K).imag) * lam
    slope0 = imf[10] / tt[10]
    i04 = np.argmin(np.abs(tt - 0.42))
    print(f"  öngörülen v-başlangıç eğimi: {slope0:.2f}  (ölçülen ~2.01)")
    print(f"  öngörülen v(0.42): {imf[i04]:.3f}  (ölçülen doyma ~0.55-0.60)")
    if K == 2:
        fw2, lam2 = fw, lam

# grafik (K=2)
fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
tt = np.linspace(1e-3, 0.5, 400)
ax = axes[0]
ax.errorbar(tw, w, yerr=sw, fmt="o", ms=3.5, c="firebrick", label="w ölçümleri (44)")
ax.plot(tt, chi(tt, fw2.x, 2).real, "k-", lw=1.3, label="Lorentz fiti (yalnız w'ye)")
ax.set_xlabel("τ"); ax.set_ylabel("w = Re χ")
ax.set_title("Gerçek kısım: fit"); ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax = axes[1]
ax.errorbar(tv, v, yerr=sv, fmt="s", ms=3.5, c="teal", label="v ölçümleri (90)")
ax.plot(tt, lam2 * np.abs(chi(tt, fw2.x, 2).imag), "k-", lw=1.3,
        label=f"ÖNGÖRÜ: λ·|Im χ| (λ={lam2:.2f})")
ax.set_xlabel("τ"); ax.set_ylabel("v")
ax.set_title("Sanal kısım: v'yi görmeden öngörü")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "66_kramers_kronig.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
