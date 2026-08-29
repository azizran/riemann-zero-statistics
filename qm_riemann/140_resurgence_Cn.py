"""
140 — RESURGENCE AVI: C(n)'DE ALÇAK SIFIR REZONANSLARI (29 Ağu)
==========================================================================
Kalem (KALEM_BK423): (4.23)'ün rezonansları gap-gecikme dilinde
n*_k = γ_k·L/2π'ye düşer. Köşegen TAM merdiven (asal-kuvvetler, τ≤1)
bu yapıyı taşır; PNT-pürüzsüz merdiven taşımaz.
A) Özdeşlik kontrolü: F = (ζ'/ζ)' − Σ log²p/(p^s−1)², Re s=1.05'te.
B) Öngörü: C_tam(n) = Σ_{q:τ≤1} 2a_q²sin²(πτ)cos(2πnτ) ile
   C_düz(n) (m=1 asalları PNT süreklisiyle değiştir; m≥2 aynı) —
   fark = aritmetik/resurgence taşıyıcısı.
C) Ölçüm: C(n) = ⟨(ds₀−μ)(ds_n−μ)⟩, n=2..60, son-300k ve orta-300k.
   Anomali: A(n) = C − yerel medyan taban; σ, rezonans-dışı laglardan.
ÖN-MÜHÜR:
  R1  son-300k: n*≈27.1'de yerel anomali, işareti tam-merdiven
      öngörüsüyle aynı, taban gürültüsünden ≥3σ.
  R2  orta-300k: anomali n*≈25.8'e KAYAR (γ₁L/2π izleme) — sabit-n
      artefaktından ayrışır.
  R3  Pürüzsüz merdiven aynı yerde ≈0 — anomali aritmetik imza.
  Kayıt: genlik vs öngörü oranı (koşullama/ekran çarpanı serbest, 139).

SONUÇ (29 Ağustos, gerçek koşudan) — R1 ✓ R2 ✓ R3 ✓, KEŞİF:
  ALTI ÇUKUR, ALTISI YERİNDE. son-300k (L=12.030): n=27/40/48'de
  −4.3σ/−4.1σ/−4.4σ (öngörü n*=27.1/40.2/47.9). orta-300k (L=11.464):
  çukurlar n=26/38/46'ya KAYDI (öngörü 25.8/38.4/45.6) — γ_k·L/2π
  izlemesi birebir; sabit-n artefaktı dışlandı. Tam-merdiven öngörüsü
  aynı yer/işaret/şekil (omuz asimetrisi dahil); PÜRÜZSÜZ merdiven
  hiçbir şey göstermiyor (≤3e-3) — imza aritmetik. Genlik oranı
  ölçüm/öngörü altı çukurda tekdüze: 0.73/0.77/0.79/0.73/0.76/0.75
  ≈ 0.76±0.03 (koşullama/ekran çarpanının uzun-gecikme hali — açık iş).
  ⇒ "GAZ KENDİ ALÇAK SIFIRLARINI TANIR": γ₁,γ₂,γ₃ iki milyon tekne
  ötede, gap-kovaryansında rezonans çukuru olarak görünüyor — sıfır
  serbest parametre. Resurgence gap dilinde ölçüldü; (4.23) kapısı
  ölçümle bağlandı.
  A-notu: özdeşlik kontrolü Re s=1.05'te %15 sapmıştı — teşhis: kesme
  kuyruğu (p^{-0.05} yavaş yakınsama), özdeşlik hatası değil; ayrı
  koşuda Re s=1.5→1.7e-3, Re s=2.0→7.3e-6 (kesmeyle ölçekleniyor ✓).
"""

import numpy as np
from pathlib import Path
from sympy import primerange
from mpmath import mp, mpf, mpc, zeta, diff as mpdiff, log as mplog

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
GAMS = [14.134725142, 21.022039639, 25.010857580]

# --- A) özdeşlik kontrolü (Re s = 1.05, mutlak yakınsak)
mp.dps = 30
s0 = mpc(mpf("1.05"), mpf("-14.0"))
lhs = mpf(0)
for p in primerange(2, 200000):
    lp = mplog(p)
    q, m = p, 1
    while m * 1.05 < 40 and q < 10**16:
        lhs += lp * lp * (mpf(q)) ** (-s0)  # q^{-s} = p^{-ms}
        q *= p; m += 1
zz2 = mpdiff(lambda t: zeta(t, derivative=1) / zeta(t), s0)
tail = mpf(0)
for p in primerange(2, 200000):
    lp = mplog(p)
    tail += lp * lp / (mpf(p) ** s0 - 1) ** 2
rhs = zz2 - tail
print(f"A) özdeşlik @s=1.05-14i: |LHS-RHS|/|RHS| = "
      f"{abs(lhs - rhs) / abs(rhs)}", flush=True)

# --- veri pencereleri
d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
Z = np.sort(np.asarray(d["zeros"], dtype=float))
N = len(Z)
PEN = {"son": Z[N - 300000:], "orta": Z[850000:1150000]}
NMAX = 60

def pk_arrays(L):
    lim = int(np.exp(L))
    om, a2, m1 = [], [], []
    for p in primerange(2, lim + 1):
        q, mm = p, 1
        while q <= lim:
            om.append(np.log(q))
            a2.append(1.0 / (np.pi**2 * mm**2 * q))
            m1.append(mm == 1)
            q *= p; mm += 1
    return np.array(om), np.array(a2), np.array(m1)

for ad, zz in PEN.items():
    g = np.diff(zz)
    mid = 0.5 * (zz[:-1] + zz[1:])
    ds = g * np.log(mid / TWO_PI) / TWO_PI - 1
    ds = ds - ds.mean()
    L = float(np.log(mid / TWO_PI).mean())
    nstars = [gk * L / TWO_PI for gk in GAMS]
    om, a2, m1 = pk_arrays(L)
    tau = om / L
    # ölçüm
    Cm = np.array([np.mean(ds[:-n] * ds[n:]) for n in range(1, NMAX + 1)])
    # öngörüler
    def pred(nn, exact_m1):
        base = 2 * a2 * np.sin(np.pi * tau) ** 2
        if exact_m1:
            return float(np.sum(base * np.cos(2 * np.pi * nn * tau)))
        # m=1'i PNT süreklisiyle değiştir, m>=2 aynı
        tt = np.linspace(np.log(2) / L, 1.0, 4000)
        f = 2 / np.pi**2 * np.sin(np.pi * tt) ** 2 * np.cos(2 * np.pi * nn * tt) / tt
        smooth = float(np.trapz(f, tt))
        rest = float(np.sum((base * np.cos(2 * np.pi * nn * tau))[~m1]))
        return smooth + rest
    Ct = np.array([pred(n, True) for n in range(1, NMAX + 1)])
    Cd = np.array([pred(n, False) for n in range(1, NMAX + 1)])
    # anomali: yerel medyan taban (±3..8 komşu, ±2 hariç)
    def anom(C):
        A = np.zeros_like(C)
        for i in range(len(C)):
            js = [j for j in range(len(C))
                  if 3 <= abs(j - i) <= 8]
            A[i] = C[i] - np.median(C[js])
        return A
    Am, At = anom(Cm), anom(Ct)
    off = [i for i in range(12, NMAX - 1)
           if all(abs((i + 1) - ns) > 2.5 for ns in nstars)]
    sig = float(np.std(Am[off]))
    print(f"\n=== {ad}: L={L:.3f}  n* = " +
          ", ".join(f"{x:.1f}" for x in nstars) + f"  σ_taban={sig:.2e}")
    for ns in nstars:
        i = int(round(ns)) - 1
        win = Am[max(0, i - 1):i + 2]
        j = i - 1 + int(np.argmax(np.abs(win)))
        print(f"  n*≈{ns:5.1f}: ölçüm A({j+1})={Am[j]:+.2e} "
              f"({Am[j]/sig:+.1f}σ)  tam-merdiven A={At[j]:+.2e}  "
              f"pürüzsüz A={anom(Cd)[j]:+.2e}", flush=True)
    print("  n :  C_ölçüm   A_ölçüm   A_tam    A_düz")
    for n in range(22, 52):
        print(f"  {n:2d}: {Cm[n-1]:+.2e} {Am[n-1]:+.2e} "
              f"{At[n-1]:+.2e} {anom(Cd)[n-1]:+.2e}")
