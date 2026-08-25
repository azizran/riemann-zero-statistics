"""
111 — BOND-ZENGİN ÇEKİRDEK: SON %25 BASAMAĞI (25 Ağustos)
==========================================================================
110c kalıntısı: konum-jitter sınıfında bond KİLİTLİ (R_nn ≈ β + kin) —
gerçeğin düşük-τ oranı (R_nn/R_p ≈ 2.75) verilemiyor ve nefes köprüsü
perdenin %59-81'inde kalıyor. Burada çekirdek zenginleştirilir:
PAYLAŞIMLI-ÇİFT bileşeni — w_n = ξ'_n + c_n ζ_n + c_{n-1} ζ_{n-1};
komşu konumlar c_nζ_n'yi paylaşır → gap-bond'a girmeyen, YALNIZ
varyansı oynatan ikinci düğme:
   β: ξ genlik-modülasyonu (R_p ve R_nn birlikte)
   γ: c² modülasyonu (yalnız R_p)
İkisiyle gerçeğin (R_p, R_nn) çifti her τ'da yakalanabilir.
(Not: log'daki "Δ²" etiketi kavramsal hedefti; uygulanan eşdeğer-amaçlı
paylaşımlı-çift sınıfıdır — bond/var tabanı −0.5/(1+c̄²) ≈ −0.45,
gerçeğin −0.52'sinden hafif sapar, dürüstçe kayıtlı.)

ÖN-MÜHÜR:
  E1  γ ayrıştırır: R_p oynar, R_nn (γ'ya) sağır kalır.
  E2  KARAR SORUSU: (R_p,R_nn) gerçeğe eşlenince D, 1−0.36τ'ya iner mi?
      İNERSE köprü kapanır (bond kanalı D'yi besliyor);
      D γ'ya sağırsa bond kanalı D'yi BESLEMİYOR → kalan %20-40 başka
      kaynaktan (o da kesin bir hüküm; basamak iki yönde de biter).
Gerçek hedefler: R_p=−0.85; R_nn(τ)≈−1.95/−1.85/−1.75/−1.3
(109 p-listesi interp.); D_gerçek = 1−0.36τ.

SONUÇ (25 Ağustos):
  E1 ✓ γ AYRIŞTIRIYOR: R_p oynar (Δ≈+0.04/γ-birimi), R_nn sağır.
  E2 — D VARYANS KANALINA bağlı: γ (varyans-modunu söndürünce) D'yi
  geri yükseltir; bond kanalı D'yi doğrudan beslemez. İki-nüfus
  doğrulama noktası (β=−3, γ=17-18; yerel gürültü güçlü ters nefes +
  kolektif bileşen dalgayla birlikte): (R_p, R_nn) çifti gerçeğe
  eşlendi (−0.93..−1.23 / −2.2..−2.9) ve D açığı gerçeğin %65-78'i
  (τ=0.10/0.30 uçlarında %78, ortada %65). HÜKÜM: faz-modülasyonlu
  GAUSSIAN jitter sınıflarının TÜMÜ ~%70 tavanında — kalan ~%25-35
  bu sınıfın DIŞINDA: adaylar (1) alt-Gauss yapı (gerçek kurt
  −0.75!), (2) üçüncü-moment (çarpıklık) dalgası ⟨η³⟩(faz),
  (3) kaynak-verteksi (dalga genliğinin kendisinin renormalizasyonu).
  Sıradaki basamak (112 adayı): alt-Gauss/çarpıklık-dalgalı merdiven.
  Yüksek-τ R_nn kilit-uyumsuzluğu da açık (gerçek −1.3 vs sınıf −2.9;
  gerçek yüksek-τ R zaten zayıf ölçülü).
"""

import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(1110)
L0 = 10.37
NZ = 200000
gbar = TWO_PI / L0

def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8

t0 = TWO_PI * np.exp(L0)
idx = np.arange(NZ + 1, dtype=float)
x = t0 + idx * gbar
for _ in range(8):
    x = x - (rvm_N(x) - rvm_N(t0) - idx) / (np.log(x / TWO_PI) / TWO_PI)

SIG_DS = 0.148
A_DS = 0.10
CB2 = 0.12                                   # paylaşımlı taban c̄²

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
    s2 = float((eta**2).mean())
    b2, *_ = np.linalg.lstsq(X, eta**2 - s2, rcond=None)
    Rp = (b2[1] * np.cos(ph) + b2[2] * np.sin(ph)) / (2 * s2 * A1)
    ee = eta[:-1] * eta[1:]
    c1 = float(ee.mean())
    mm = 0.5 * (m[:-1] + m[1:])
    X2 = np.vstack([np.ones_like(mm), np.cos(om * mm), np.sin(om * mm),
                    np.cos(2 * om * mm), np.sin(2 * om * mm)]).T
    b3, *_ = np.linalg.lstsq(X2, ee - c1, rcond=None)
    Rnn = (b3[1] * np.cos(ph) + b3[2] * np.sin(ph)) / (2 * c1 * A1)
    return A1, Rp, Rnn, c1 / s2

# jitter ölçeği: var(gap-noise) = 2σξ²(1+c̄²) hedef (SIG_DS·ḡ)²
sxi = SIG_DS * gbar / np.sqrt(2 * (1 + CB2))

RNN_GER = {0.10: -1.95, 0.15: -1.85, 0.20: -1.75, 0.30: -1.30}
print(f"{'τ':>5} {'β':>5} {'γ':>5} {'D':>6} {'R_p':>7} {'R_nn':>7} "
      f"{'c₁/σ²':>6}")
SON = {}
for tau in [0.10, 0.15, 0.20, 0.30]:
    om = tau * L0
    UA = A_DS / (2 * np.sin(np.pi * tau)) * gbar
    s_ph = -np.sin(om * x)
    for beta, gam in [(0.0, 0.0), (-2.0, 0.0), (-2.0, 4.0), (-3.0, 0.0),
                      (-3.0, 4.0), (-3.0, 8.0), (-4.0, 8.0), (-4.0, 12.0)]:
        xi = rng.normal(0, sxi, NZ + 1) * (1 + beta * A_DS * s_ph)
        ze = rng.normal(0, sxi, NZ + 1)
        cn = np.sqrt(CB2 * np.clip(1 + gam * A_DS * s_ph, 0.05, None))
        w = xi + cn * ze + np.concatenate([[0], (cn * ze)[:-1]])
        y = np.sort(x + w)
        z = np.sort(y + UA * np.cos(om * y))
        A1, Rp, Rnn, c1r = olc(z, om)
        SON[(tau, beta, gam)] = (A1 / A_DS, Rp, Rnn)
        print(f"{tau:>5.2f} {beta:>5.1f} {gam:>5.1f} {A1/A_DS:>6.3f} "
              f"{Rp:>+7.3f} {Rnn:>+7.3f} {c1r:>+6.2f}", flush=True)
    print(f"      [gerçek hedef: R_p −0.85, R_nn {RNN_GER[tau]:+.2f}, "
          f"D {1-0.36*tau:.3f}]")
