"""
37 — CONREY-GHOSH SANITY CHECK (16 Ağustos 2026)
=================================================

Kanıtlı teorem (CG 1985, RH altında):
    (1/N(T)) Σ max_{γ<t≤γ⁺} |ζ(½+it)|² ~ ½(e²−5)·log T ≈ 1.19453·log T

Bizim 36_T100k.npz verisi tam bu objeyi içeriyor (her aralıkta max|Z|).
Test: mean(max²) eğimi log(t/2π)'ye karşı 1.19453'e oturuyor mu?

Oturuyorsa → pipeline kanıtlı teoremle tutarlı, r=0.83 ölçümü güvenilir zemin üstünde.
Oturmuyorsa → ya optimizer bias var (36'da xatol=0.05, maxiter=10 kabaydı) ya da hata.

Ek: optimizer-bias sondası — rastgele aralıklarda max|Z| ince taramayla yeniden
hesaplanır, kayıtlı değerle karşılaştırılır.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import mpmath as mp
import time

mp.mp.dps = 10  # 36 ile aynı precision

HERE = Path(__file__).resolve().parent
CG = 0.5 * (np.e**2 - 5)  # 1.19452804946...

d = np.load(HERE / "36_T100k.npz")
intervals, max_amps, t_mid = d["intervals"], d["max_amps"], d["t_mid"]
# Aralık sınırları t_mid ± interval/2'den birebir geri çıkar (zeros dosyası gerekmez)
gam_lo = t_mid - intervals / 2
gam_hi = t_mid + intervals / 2

print(f"Veri: {len(intervals)} aralık, t ∈ [{gam_lo[0]:.1f}, {gam_hi[-1]:.1f}]")
print(f"CG sabiti ½(e²−5) = {CG:.8f}\n")

# ---------------------------------------------------------------
# 1) PENCERE ORTALAMALARI + EĞİM FİTİ
# ---------------------------------------------------------------
max2 = max_amps**2
logt = np.log(t_mid / (2 * np.pi))

n_w = 40
edges = np.geomspace(t_mid[0], t_mid[-1] * 1.0001, n_w + 1)
w_logt, w_mean, w_n = [], [], []
for i in range(n_w):
    m = (t_mid >= edges[i]) & (t_mid < edges[i + 1])
    if m.sum() < 50:
        continue
    w_logt.append(logt[m].mean())
    w_mean.append(max2[m].mean())
    w_n.append(m.sum())
w_logt, w_mean, w_n = map(np.array, (w_logt, w_mean, w_n))

# Ağırlıklı doğrusal fit: mean(max²) = a·log(t/2π) + b
A = np.vstack([w_logt, np.ones_like(w_logt)]).T
coef, res_, *_ = np.linalg.lstsq(A * np.sqrt(w_n)[:, None],
                                 w_mean * np.sqrt(w_n), rcond=None)
a_fit, b_fit = coef
print("1) EĞİM TESTİ  mean(max²) = a·log(t/2π) + b")
print(f"   a = {a_fit:.4f}   (CG: {CG:.4f},  sapma {100*(a_fit/CG-1):+.2f}%)")
print(f"   b = {b_fit:.4f}   (alt-mertebe terimi, CG leading-order'da öngörmez)\n")

# ---------------------------------------------------------------
# 2) GLOBAL CESÀRO ORANI (CG'nin birebir ifadesi), büyüyen T
# ---------------------------------------------------------------
print("2) CESÀRO ORANI  (1/N)Σmax² / log(·)")
print(f"   {'T':>8} {'N':>7} {'/log(T/2π)':>11} {'/log T':>8}")
csum = np.cumsum(max2)
for T in [1000, 3000, 10000, 30000, gam_hi[-1]]:
    k = np.searchsorted(gam_hi, T)
    if k < 100:
        continue
    mean_ = csum[k - 1] / k
    print(f"   {T:>8.0f} {k:>7} {mean_/np.log(T/(2*np.pi)):>11.4f} "
          f"{mean_/np.log(T):>8.4f}")
print(f"   (CG asimptotik: → {CG:.4f}; log T ile log(T/2π) farkı alt-mertebe)\n")

# ---------------------------------------------------------------
# 3) OPTIMIZER-BIAS SONDASI
#    36'daki kaba arama (xatol=0.05, maxiter=10) maksimumu kaçırıyor mu?
#    250 rastgele aralıkta iki aşamalı ince grid ile yeniden hesapla.
# ---------------------------------------------------------------
rng = np.random.default_rng(37)
idx = np.sort(rng.choice(len(intervals), 250, replace=False))

def fine_max(a, b, n1=60, n2=40):
    ts = np.linspace(a, b, n1)
    vs = [abs(float(mp.siegelz(t))) for t in ts]
    j = int(np.argmax(vs))
    lo = ts[max(j - 1, 0)]
    hi = ts[min(j + 1, n1 - 1)]
    ts2 = np.linspace(lo, hi, n2)
    vs2 = [abs(float(mp.siegelz(t))) for t in ts2]
    return max(max(vs), max(vs2))

print("3) OPTIMIZER-BIAS SONDASI (250 aralık, ince grid)")
t0 = time.time()
probe = np.array([fine_max(gam_lo[i] + 0.005, gam_hi[i] - 0.005) for i in idx])
stored = max_amps[idx]
ratio2 = probe**2 / stored**2
print(f"   süre: {time.time()-t0:.0f} s")
print(f"   max²_ince / max²_kayıtlı:  ortalama = {ratio2.mean():.4f}, "
      f"medyan = {np.median(ratio2):.4f}")
print(f"   kaçırma >%1 olan aralık oranı: {(ratio2 > 1.01).mean()*100:.1f}%")
corr = ratio2.mean()
print(f"   → bias-düzeltmeli eğim: a·corr = {a_fit*corr:.4f} "
      f"(CG'ye sapma {100*(a_fit*corr/CG-1):+.2f}%)\n")

# ---------------------------------------------------------------
# 4) GRAFİK
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5))

ax = axes[0]
ax.scatter(w_logt, w_mean, s=np.sqrt(w_n), c="steelblue", zorder=3,
           label="pencere ortalamaları (40 pencere)")
xs = np.linspace(w_logt.min(), w_logt.max(), 100)
ax.plot(xs, a_fit * xs + b_fit, "r-", lw=1.5,
        label=f"fit: a={a_fit:.3f}, b={b_fit:.3f}")
ax.plot(xs, CG * xs + b_fit, "g--", lw=1.5,
        label=f"CG eğimi ½(e²−5)={CG:.4f} (+aynı b)")
ax.set_xlabel("log(t/2π)")
ax.set_ylabel("mean(max|Z|²)")
ax.set_title("Conrey–Ghosh 1985 sanity check")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

ax = axes[1]
ax.plot(w_logt, w_mean / w_logt, "o-", ms=4, c="steelblue",
        label="mean(max²)/log(t/2π)")
ax.axhline(CG, color="g", ls="--", label=f"CG asimptotik {CG:.4f}")
ax.set_xlabel("log(t/2π)")
ax.set_ylabel("oran")
ax.set_title("Orana yakınsama (alt-mertebe b/log t etkisiyle)")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

plt.tight_layout()
out = HERE / "37_conrey_ghosh_check.png"
plt.savefig(out, dpi=110)
print(f"Grafik: {out.name}")

# ---------------------------------------------------------------
# 5) KARAR
# ---------------------------------------------------------------
dev = abs(a_fit * corr / CG - 1)
print("\n" + "=" * 60)
if dev < 0.05:
    print(f"SONUÇ: eğim CG sabitine %{dev*100:.1f} yakınlıkta → PIPELINE DOĞRULANDI")
    print("r=0.83 ölçümü kanıtlı teoremle tutarlı zeminde.")
else:
    print(f"SONUÇ: eğim CG'den %{dev*100:.1f} sapıyor → incelenmeli!")
    print("Olası sebepler: optimizer bias, alt-mertebe terimler, veri hatası.")
