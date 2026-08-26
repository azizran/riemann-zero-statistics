"""
121 — TERSİNME ÖZDEŞLİĞİ HAKEMİ: YOĞUNLUK-MİMARİSİ R=+2 Mİ? (26 Ağu)
==========================================================================
"Neden −2" kaleminin dökümü: sayım kısıtı (her gap'te ∫ρ=1) kesin
haritayı verir: ds = (W−ν)/(1−W+ν) → kinematik nefes TÜM momentlerde
+2 (varyans (1−W)^{-4}, m₃ (1−W)^{-6}). İnvolüsyon R_g + R_ρ = 2:
sabit nokta +1 = dengenin adyabatiği (termal aynanın açıklaması);
R_ρ=0 → R_g=+2; ölçülen R_g=−2 ⟺ R_ρ=+4 ("yoğunluk dilinde kinematiğin
iki katı" — yeni soru).

BU HAKEM: yoğunluk-mimarili sentetik — ρ(t) = ρ̄(1 + A_ρcos(ωt) + ν(t)),
ν bant-sınırlı düz Gaussian alan (gap-ölçeği korelasyonlu, FAZ-BAĞIMSIZ);
sıfırlar Φ(t)=∫ρ = n tersinmesiyle. Standart zincirle (D, R_p, R_nn, R₃).
ÖN-MÜHÜR:
  P1  R_p ≈ +2 (özdeşliğin kinematiği; boyalı-mimariden [+1] farklı
      olmalı — mimari ayrımının kanıtı).
  P2  R₃(m₃-norm) ≈ +2.
  P3  Kontrast için aynı ν ile boyalı-mimari satırı: R_p ≈ +1.
Çıkarsa: özdeşlik hakem-onaylı; "neden −2" resmen "yoğunluk dilinde
neden +4" olur ve bütçe türetimi o dile taşınır.

SONUÇ (26 Ağustos) — P1 ✓ P2 ✓ P3 ✓, ÖZDEŞLİK HAKEM-ONAYLI:
  yoğunluk-mimarisi: R_p = +1.90/+1.62, R₃(m₃) = +2.10/+1.61
  (τ=0.15/0.30; +2 yasası, τ-zayıflamasıyla ✓).
  boyalı-mimari (aynı ν): +0.96/+0.71 (+1 yasası ✓) — mimari ayrımı
  tam kalem oranında (2×).
  BONUS: yoğunluk-mimarisi faz-bağımsız gürültüyle bile PERDE üretiyor
  (D 0.970/0.865) — tersinme nonlineerliği + pozitif nefes de koheran
  okumayı yiyor: perde-mekanizması |modülasyon|-sürücülü (işaretten
  bağımsız), 110c yasasıyla tutarlı.
  KALEM DURUMU: (i) tersinme özdeşliği kesin ve hakemli; (ii) denge
  +1'i = involüsyonun (R_g + R_ρ = 2) öz-eş noktası — termal aynanın
  ilk-ilke açıklaması; (iii) soru nihai biçiminde: GERÇEK GAZ YOĞUNLUK
  DİLİNDE R_ρ = +4 — artık-yoğunluk gürültüsü, yoğunluk dalgasına
  kinematik oranın TAM İKİ KATIYLA biner. "Neden −2" = "neden 2×
  kinematik" (iki kuadratür? |Z|²-kanalı bağlantısı [R_w>0] aday).
"""

import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(121)
L0 = 10.37
NZ = 200000
gbar = TWO_PI / L0
rho0 = 1.0 / gbar
A_DS = 0.10

def olc(z, om):
    g = np.diff(z)
    m = 0.5 * (z[:-1] + z[1:])
    ds = g / gbar - 1
    ds = ds - ds.mean()
    X = np.vstack([np.ones_like(m), np.cos(om * m), np.sin(om * m),
                   np.cos(2 * om * m), np.sin(2 * om * m)]).T
    b, *_ = np.linalg.lstsq(X, ds, rcond=None)
    A1 = np.hypot(b[1], b[2]); ph = np.arctan2(b[2], b[1])
    eta = ds - X @ b
    s2 = float((eta**2).mean()); m3 = float((eta**3).mean())
    b2, *_ = np.linalg.lstsq(X, eta**2 - s2, rcond=None)
    Rp = (b2[1] * np.cos(ph) + b2[2] * np.sin(ph)) / (2 * s2 * A1)
    ee = eta[:-1] * eta[1:]
    c1 = float(ee.mean())
    mm = 0.5 * (m[:-1] + m[1:])
    X2 = np.vstack([np.ones_like(mm), np.cos(om * mm), np.sin(om * mm),
                    np.cos(2 * om * mm), np.sin(2 * om * mm)]).T
    b3, *_ = np.linalg.lstsq(X2, ee - c1, rcond=None)
    Rnn = (b3[1] * np.cos(ph) + b3[2] * np.sin(ph)) / (2 * c1 * A1)
    b4, *_ = np.linalg.lstsq(X, eta**3 - m3, rcond=None)
    R3m = ((b4[1] * np.cos(ph) + b4[2] * np.sin(ph)) / (3 * m3 * A1)
           if abs(m3) > 1e-9 else np.nan)
    return A1 / A_DS, Rp, Rnn, R3m, np.sqrt(s2)

def duz_alan(T, dt, lcor, sig):
    """Bant-sınırlı düz Gaussian alan (korelasyon ~lcor)."""
    n = int(T / dt)
    w = rng.normal(0, 1, n)
    k = int(max(3, lcor / dt))
    ker = np.exp(-0.5 * np.linspace(-3, 3, 6 * k + 1)**2)
    ker /= np.sqrt((ker**2).sum())
    f = np.convolve(w, ker, mode="same")
    return f * sig / f.std()

print(f"{'mimari':>10} {'τ':>5} {'D':>6} {'R_p':>7} {'R_nn':>7} "
      f"{'R₃(m₃)':>8} {'σ_η':>6}")
T = NZ * gbar * 1.05
dt = gbar / 8
tt = np.arange(0, T, dt)
for tau in [0.15, 0.30]:
    om = tau * L0
    # yoğunluk dalga-genliği: gap-dalga A_DS'yi verecek şekilde (1. mertebe eş)
    nu = duz_alan(T, dt, gbar, 0.148 * 0.85)   # ds-artığı ~0.148 hedefi (ayar)
    for mimari in ["yoğunluk", "boyalı"]:
        if mimari == "yoğunluk":
            rho = rho0 * (1 + A_DS * np.cos(om * tt) + nu)
            Phi = np.cumsum(rho) * dt
            z = np.interp(np.arange(1, NZ + 2, dtype=float), Phi, tt)
        else:
            rho = rho0 * (1 + nu)
            Phi = np.cumsum(rho) * dt
            y = np.interp(np.arange(1, NZ + 2, dtype=float), Phi, tt)
            UA = A_DS / (2 * np.sin(np.pi * tau)) * gbar
            z = np.sort(y + UA * np.cos(om * y))
        D, Rp, Rnn, R3m, se = olc(z, om)
        print(f"{mimari:>10} {tau:>5.2f} {D:>6.3f} {Rp:>+7.3f} "
              f"{Rnn:>+7.3f} {R3m:>+8.3f} {se:>6.3f}", flush=True)
print("\nÖN-MÜHÜR: yoğunluk-mimarisi R_p≈+2, boyalı ≈+1; "
      "GERÇEK gaz: −2 (ayna).")
