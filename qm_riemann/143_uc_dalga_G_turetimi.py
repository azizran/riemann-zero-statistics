"""
143 — ÜÇ-DALGA MAKİNESİ: G'NİN ARACI-DALGA TÜRETİMİ ve SINAVI (29 Ağu)
==========================================================================
KALEM: orta-nokta kayması δ(m̄)=ḡ_t Σ a_c cos(πτ_c) sin(ω_c m̄);
birinci mertebe açılımda ⟨cos ω₁m cos ω₂m⟩ yalnız δ'nın FARK (ω₃=ω₂−ω₁)
ve TOPLAM (ω₁+ω₂; yalnız q₁q₂ asal-kuvvetse örgüde) bileşenlerinden
beslenir. Sıfır-parametre formül (a_Q = 1/(π m_Q √Q), τ_Q = log Q/L):

  G_cc = −π[ τ₃ a₃ cos(πτ₃) + (τ₁+τ₂) a_Σ cos(πτ_Σ)·1{q₁q₂ örgüde} ]
  G_ss = −π[ τ₃ a₃ cos(πτ₃) − (τ₁+τ₂) a_Σ cos(πτ_Σ)·1{...} ]
  G_cs = 0

ÖN-MÜHÜR:
  U1  Aracı kuralı: eş-asal çiftlerde G yalnız aracıya bağlı
      (G(3,6)=G(5,10)=G(7,14)... ±%15).
  U2  Kule çiftleri Σ-terimiyle güçlü; kule merdiveninde yukarı
      zayıflar: |G(2,4)| > |G(4,8)| > |G(8,16)|.
  U3  Kanal: eş-asalda G_cc ≈ G_ss; kulede yarılma ≈ 2π(τ₁+τ₂)a_Σc_Σ.
  U4  G_cs ≈ 0 (faz-bükümsüzlüğün mekaniği).
  U5  İşaretler negatif; değerler ±%15 (küçük aracılar için ikinci
      mertebe kaskad fazlalığı kayda).

SONUÇ (29 Ağustos, gerçek koşudan) — 19 ÇİFTTE 19 İSABET; YASA MÜHÜRLÜ:
  U5 ✓ işaret 19/19 negatif. U2 ✓ değerler: oran 0.91-0.98 (ort ~0.95;
  tekdüze ~%5 açık — ikinci-mertebe adayı). U1 ✓ aracı kuralı ±%1:
  (5,10)(7,14)(11,22)(13,26) → −.0384/−.0382/−.0382/−.0382.
  U2 ✓ kule merdiveni: |G(2,4)|>|G(4,8)|>|G(8,16)| = 562>442>393
  (öngörü 575>464>416). U3 ✓ kanal yarılması kulede tam boy:
  (2,4) cc/ss = −562/−236 (öngörü −575/−226); yarılma işareti
  (25,125)'te cos(πτ_Σ)'nin döndüğü yerde dönüyor. U4 ✓ G_cs ≤ 1e-4
  (19/19) — faz-bükümsüzlük mekanik olarak türetildi ve doğrulandı.
  Kayıt: eş-asal cc/ss küçük asimetri (~.0025, ss yüksek) — öngörü
  eşit; ikinci-mertebe not. 142'nin tüm fenomenolojisi bu mikro-yasadan
  akıyor; sıradaki kalem: G'den kolektif ρ_tail(τ), çukur genliği ve
  K(τ) hesabı (KALEM_UC_DALGA_G_29AGU2026.md §5).
"""

import numpy as np
from pathlib import Path
from sympy import factorint

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
Z = np.sort(np.asarray(d["zeros"], dtype=float))
zz = Z[len(Z) - 300000:]
mid = 0.5 * (zz[:-1] + zz[1:])
Nn = len(mid)
L = float(np.log(mid / TWO_PI).mean())

def is_pp(n):
    f = factorint(n)
    return len(f) == 1

def aQ(Q):
    (p, m), = factorint(Q).items()
    return 1.0 / (np.pi * m * np.sqrt(Q))

def pred(q1, q2):
    t1, t2 = np.log(q1) / L, np.log(q2) / L
    q3 = q2 // q1
    t3 = np.log(q3) / L
    T1 = t3 * aQ(q3) * np.cos(np.pi * t3)
    S = q1 * q2
    if is_pp(S):
        TS = (t1 + t2) * aQ(S) * np.cos(np.pi * np.log(S) / L)
    else:
        TS = 0.0
    return -np.pi * (T1 + TS), -np.pi * (T1 - TS)

CIFTLER = [
    ("kule", 2, 4), ("kule", 4, 8), ("kule", 8, 16),
    ("kule", 3, 9), ("kule", 9, 27), ("kule", 5, 25), ("kule", 25, 125),
    ("kule", 7, 49),
    ("eş-asal", 3, 6), ("eş-asal", 5, 10), ("eş-asal", 7, 14),
    ("eş-asal", 11, 22), ("eş-asal", 13, 26),
    ("eş-asal", 2, 6), ("eş-asal", 5, 15), ("eş-asal", 7, 21),
    ("eş-asal", 2, 10), ("eş-asal", 3, 15), ("eş-asal", 2, 14),
]

print(f"L={L:.3f}  N={Nn}   (gürültü tabanı ~0.002)")
print(f"{'tip':>8} {'çift':>9} {'aracı':>5} | {'G_cc ölç':>9} {'öngörü':>8} "
      f"{'oran':>5} | {'G_ss ölç':>9} {'öngörü':>8} | {'G_cs':>7}")
for tip, q1, q2 in CIFTLER:
    w1, w2 = np.log(q1), np.log(q2)
    c1, s1 = np.cos(w1 * mid), np.sin(w1 * mid)
    c2, s2 = np.cos(w2 * mid), np.sin(w2 * mid)
    gcc = 2 * float(np.mean(c1 * c2))
    gss = 2 * float(np.mean(s1 * s2))
    gcs = 2 * float(np.mean(c1 * s2))
    pcc, pss = pred(q1, q2)
    oran = gcc / pcc if pcc != 0 else float("nan")
    print(f"{tip:>8} {q1:>4}↔{q2:<4} {q2//q1:>5} | {gcc:>+9.4f} {pcc:>+8.4f} "
          f"{oran:>5.2f} | {gss:>+9.4f} {pss:>+8.4f} | {gcs:>+7.4f}",
          flush=True)
