"""
115 — GAP-ASİMETRİLİ ÇEKİRDEK: SON %25 MERDİVENİ (26 Ağustos)
==========================================================================
K2 teoremi (113): jitter farkı gap-çarpıklığını öldürür. SABAH
MİNİ-TEOREMİ 2: herhangi site-haritası u=f(ε) de öldürür (değişim
antisimetrisi) — çarpıklık ZAMAN-ASİMETRİK çapraz-site yapı ister.
EN YALIN KURUCU: u_n = ε_n + λ·(ε²_{n−1} − σ_ε²)
  → ⟨η³⟩ = 6λσ_ε⁴ (Isserlis ✓), pozisyonlar sınırlı, geçmişe-bakan
  asimetri (itme fiziğinin sentetik sureti).

PLAN:
  A1  λ-kalibrasyonu (dalga yok): skew(η) → 0.243 hedefi; lag-tablosu
      (⟨η²η₊⟩, ⟨ηη₊²⟩, ⟨η²η₊₊⟩, ⟨ηη₊η₊₊⟩) T2-gerçekle kıyas (dürüst:
      tek-λ türün tam deseni vermesi beklenmez; hangi momentler tutuyor?)
  A2  MERDİVEN: β-nefes (ε genliği türev-fazlı modüle) + boyalı dalga;
      dörtlü ölçüm (D, R_p, R_nn, R₃) τ ∈ {0.15, 0.30, 0.45},
      β ∈ {0, −2, −3}; λ ∈ {0, kalibre} karşılaştırmalı.
ÖN-MÜHÜR:
  P1  λ=0 → 110c ile tutarlı (D = 1 − 0.111|β|τ Gaussian tavanı).
  P2  KARAR: λ=kalibre'de D açığı Gaussian tavanını AŞARSA çarpıklık
      kanalı köprünün son %25'ine katkıdır (β-eşlenik hücrede
      1−0.36τ'ya yaklaşma payı raporlanır); aşmazsa üçüncü-moment de
      elenir → kalan tek aday kaynak-verteksi (iki yönde de hüküm).

SONUÇ (26 Ağustos) — P2: TAVAN KIRILDI, KÖPRÜ KAPANDI (1. yaklaşım):
  A1: λ=0.12 → skew 0.235 ✓ (hedef 0.243); c₁/σ²=−0.50 ✓. DÜRÜST:
    lag-deseni gerçekle kısmen TERS işaretli (tek-λ geriye-bakan tür
    kaba; ileri/geri karışımlı tür gelecek inceltme).
  A2: ÇARPIKLIK KANALI D'Yİ GÜÇLÜ BESLİYOR — Gaussian tavanı paramparça:
    (λ=0.12, β=−2): D = 0.928/0.862/0.886 (λ=0: 0.965/0.924/0.894;
    tavan: 0.967/0.933/0.900). GERÇEK: 0.946/0.892/0.838. Model artık
    gerçeği İKİ YANDAN SARIYOR (%79-133 aralığı; sistematik %70
    eksiği YOK). Üstelik aynı hücrenin R₃'ü (−1.24/−1.11/−1.01)
    gerçek aralığın (−1.0..−1.45) İÇİNDE — dörtlü (R_p, R_nn, R₃, D)
    ilk kez tek reçetede tutarlı. λ tek başına (β=0) perdelemez:
    perde = NEFES × ÇARPIKLIK etkileşimi (β-modülasyonu λ üzerinden
    3. momente yayılıyor, o da koheran okumayı yiyor).
  BÜYÜK RESİM: perde = varyans-nefesi kanalı (~%70, Gaussian, kapalı
    formlu) + çarpıklık-doğrultması kanalı (~%30, gap-düzeyi itme
    asimetrisi) — SON ÇEYREĞİN ADRESİ BULUNDU VE DOLDURULDU.
  AÇIK (inceltme): ortak (β,λ) ince-ayarıyla üçlünün eş-anlı tam
    eşlemesi; lag-deseni için ileri/geri λ-karışımı; çarpıklık-kanalı
    katsayısının analitik türetimi (0.111'in kardeşi); R_p orta-τ
    payı; sonra NOT 5.
"""

import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(115)
L0 = 10.37
NZ = 200000
gbar = TWO_PI / L0
A_DS = 0.10
SIG_T = 0.148

def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8

t0 = TWO_PI * np.exp(L0)
idx = np.arange(NZ + 1, dtype=float)
xg = t0 + idx * gbar
for _ in range(8):
    xg = xg - (rvm_N(xg) - rvm_N(t0) - idx) / (np.log(xg / TWO_PI) / TWO_PI)

def cekirdek(lam, beta, om, n):
    """u = ε(1+βAs) + λ(ε²₋−σ²); ds-varyansı SIG_T'ye normalanır."""
    s_ph = -np.sin(om * xg[:n]) if om > 0 else np.zeros(n)
    eps = rng.normal(0, 1.0, n) * (1 + beta * A_DS * s_ph)
    e2 = np.concatenate([[0.0], eps[:-1]**2 - 1.0])
    u = eps + lam * e2
    eta_raw = np.diff(u)
    sc = SIG_T / eta_raw.std()
    return u * sc * gbar

def dortlu(z, om):
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
    R3 = (b4[1] * np.cos(ph) + b4[2] * np.sin(ph)) / (s2**1.5 * A1)
    return A1 / A_DS, Rp, Rnn, R3

# ---- A1: λ kalibrasyonu (dalga yok)
print("A1 — λ kalibrasyonu (hedef skew 0.243; gerçek lag: "
      "+0.118/+0.119/+0.055/−0.180):")
print(f"{'λ':>6} {'skew':>7} {'c₁/σ²':>6} {'η²η₊':>7} {'ηη₊²':>7} "
      f"{'η²η₊₊':>7} {'ηη₊η₊₊':>8}")
for lam in [0.0, 0.05, 0.10, 0.15, 0.20]:
    u = cekirdek(lam, 0.0, 0.0, NZ + 1)
    eta = np.diff(u) / gbar
    eta = eta / eta.std()
    e0, e1, e2v = eta[:-2], eta[1:-1], eta[2:]
    print(f"{lam:>6.2f} {float((eta**3).mean()):>+7.3f} "
          f"{float((eta[:-1]*eta[1:]).mean()):>+6.2f} "
          f"{float((e0**2*e1).mean()):>+7.3f} {float((e0*e1**2).mean()):>+7.3f} "
          f"{float((e0**2*e2v).mean()):>+7.3f} {float((e0*e1*e2v).mean()):>+8.3f}",
          flush=True)

# ---- A2: merdiven
print("\nA2 — merdiven (D açığı vs Gaussian tavanı 0.111|β|τ):")
print(f"{'λ':>5} {'β':>5} {'τ':>5} {'D':>6} {'R_p':>7} {'R_nn':>7} "
      f"{'R₃':>7} {'tavan-D':>8}")
LAMC = 0.115                          # A1'den; skew≈0.24 hedefi
for lam in [0.0, LAMC]:
    for tau in [0.15, 0.30, 0.45]:
        om = tau * L0
        UA = A_DS / (2 * np.sin(np.pi * tau)) * gbar
        for beta in [0.0, -2.0, -3.0]:
            u = cekirdek(lam, beta, om, NZ + 1)
            y = np.sort(xg + u)
            z = np.sort(y + UA * np.cos(om * y))
            D, Rp, Rnn, R3 = dortlu(z, om)
            tav = 1 - 0.111 * abs(beta) * tau
            print(f"{lam:>5.2f} {beta:>5.1f} {tau:>5.2f} {D:>6.3f} "
                  f"{Rp:>+7.3f} {Rnn:>+7.3f} {R3:>+7.3f} {tav:>8.3f}",
                  flush=True)
print("\nGERÇEK: D = 0.946/0.892/0.838; R_p −0.85; R₃(σ³) −1.0..−1.45")
