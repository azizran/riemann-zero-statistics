"""
147 — C₀ KALEMİ: SOĞURMA-NEFESİ KANALI ve VARYANS-KANAL ÇAPRAZ SINAVI
==========================================================================
(29 Ağu) K₀ = κ_B(−1.845) + C₀(−1.03). ADAY: soğurma nefes alır —
ρ = A·e^{−ω g_yerel} ise dalga q'nun gerdiği yerde kuyruk sağkalımı
δlnρ = −2πτ_Q·R·ds_q ile modüle olur (R: tepki katsayısı; tam yerel
yanıt R=1). Kanal katkıları (ρw ağırlıklı, w=a²sin²πτ):
  bond   : κ_C = −2πR·Σρw τ cos2πτ / Σρw cos2πτ
  lag-0  : κ_C0 = −2πR·Σρw τ / Σρw
  B-kanalı lag-0: κ_B0 = Σρw·2πτ·cot(πτ) / Σρw   (genlik-mod; lag yok)
PLAN: (i) VARYANS kanalını İLK KEZ haritala: gerçek çizgilerde
R0 = η²-modülasyon/(2σ²A1) (ekran-gölge düzeltmesiyle R0·r) ve boyalı
çizgi-dışı ad0(τ*); (ii) bond'dan R'yi çöz (K₀ − κ_B = −3.8R·düzelt);
(iii) aynı R ile lag-0 çekirdek öngörüsü κ_B0+κ_C0'ı, ölçülen
(R0·r − ad0) ile karşılaştır.
ÖN-MÜHÜR:
  V1  R0(τ_p) ve ad0(τ*) ölçülür (yeni kanal haritası — kayıt).
  V2  Bond'dan çözülen tek R, lag-0 çekirdeğini ±%20 vurursa
      soğurma-nefesi kanalı MÜHÜR (C₀'ın kaynağı; R yeni nesne).
  V3  Tutmazsa C₀ AÇIK kalır; iki yeni ölçüm yine kârda.

SONUÇ (29 Ağustos, gerçek koşudan) — TEK-R DÜŞTÜ; İKİ MÜHÜR + ADRES:
  V1 ✓ İki yeni harita: R0 gerçek çizgiler (ham −0.70→−0.88; ekran-
  düzeltmeli plato −0.79±0.02), boyalı ad0 ≈ +0.85 DÜZ (bond ad1
  kontrolü 137'nin +0.92→+1.34'ünü birebir üretti — zincir sağlam).
  BÜYÜK MÜHÜR — VARYANS KANALI PARAMETRESİZ KAPANDI: ölçülen lag-0
  çekirdek R0·r − ad0 = −1.63±0.04 (yedi çizgi, düz); öngörü
  κ_B0 = Σρw·2πτcot(πτ)/Σρw = −1.578 — %3! (πτcotπτ ailesi, ρ-ağırlıklı,
  kanal-özgül; soğurma-nefesi katkısı YOK: R_lag0 ≈ 0.01.)
  V2/V3: tek-R modeli RED (bond R=0.286 vs lag-0 ≈0) — ama ayrıştırıcı:
  C₀ = −1.03 yalnız BOND/lag yapısında yaşıyor: gap-aşırı koherans
  cezasının (γ₁) nefes modülasyonu adayı. İpucu: 144'ün c₁ sağlamasındaki
  0.77 ≈ γ₁ (kapanmayan %30). SIRADAKİ: bant-çözünürlüklü BOND
  spektroskopisi → γ₁(τ) haritası → C₀'ın kapalı formu.
"""

import numpy as np
from pathlib import Path
from sympy import primerange, factorint

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
Z = np.sort(np.asarray(d["zeros"], dtype=float))
zz0 = Z[len(Z) - 300000:]

def pk(lim):
    out = []
    for p in primerange(2, lim + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

def chunked_fit(y, tmid, freqs, chunk=40000):
    tt = (tmid - tmid.mean()) / (tmid[-1] - tmid[0])
    C = 3 + 2 * len(freqs)
    XtX = np.zeros((C, C)); Xty = np.zeros(C)
    fr = np.array(freqs)
    for s0 in range(0, len(y), chunk):
        sl = slice(s0, min(s0 + chunk, len(y)))
        arg = np.outer(tmid[sl], fr)
        Xc = np.empty((sl.stop - sl.start, C))
        Xc[:, 0] = 1; Xc[:, 1] = tt[sl]; Xc[:, 2] = tt[sl]**2
        Xc[:, 3::2] = np.cos(arg); Xc[:, 4::2] = np.sin(arg)
        XtX += Xc.T @ Xc; Xty += Xc.T @ y[sl]
        del Xc, arg
    b = np.linalg.solve(XtX, Xty)
    fit = np.empty(len(y))
    for s0 in range(0, len(y), chunk):
        sl = slice(s0, min(s0 + chunk, len(y)))
        arg = np.outer(tmid[sl], fr)
        fit[sl] = (b[0] + b[1]*tt[sl] + b[2]*tt[sl]**2 +
                   np.cos(arg) @ b[3::2] + np.sin(arg) @ b[4::2])
        del arg
    return b, fit

def kanallar(z, hedefler, ekstra=()):
    g = np.diff(z)
    m = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(m / TWO_PI)
    L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    qs = pk(min(int(np.exp(0.52 * L)), 720))
    freqs = [np.log(q) for q in qs] + list(ekstra)
    b1, f1 = chunked_fit(ds, m, freqs)
    eta = ds - f1
    s2 = float(np.var(eta))
    ee0 = eta * eta - s2
    b0, _ = chunked_fit(ee0, m, freqs)
    ee1 = eta[:-1] * eta[1:]
    c1 = float(ee1.mean())
    mm = 0.5 * (m[:-1] + m[1:])
    b3, _ = chunked_fit(ee1 - c1, mm, freqs)
    out = {}
    for f in hedefler:
        i = freqs.index(f)
        tau = f / L
        cg, sg = b1[3 + 2*i], b1[4 + 2*i]
        A1 = float(np.hypot(cg, sg)); ph = np.arctan2(sg, cg)
        def proj(bb):
            return float((bb[3 + 2*i] * np.cos(ph) +
                          bb[4 + 2*i] * np.sin(ph)))
        R0 = proj(b0) / (2 * s2 * A1)
        R1 = proj(b3) / (2 * c1 * A1)
        out[f] = (tau, A1, R0, R1)
    return out, L, s2, c1

# (i-a) gerçek çizgiler
PR = [2, 3, 5, 7, 11, 13, 17]
out, L, s2, c1 = kanallar(zz0, [np.log(p) for p in PR])
print(f"L={L:.3f}  σ²={s2:.4f}  c₁={c1:+.5f}")
print("GERÇEK ÇİZGİLER (R0: varyans kanalı; R1: bond — kontrol):")
RLINE = {2: 1.001, 3: 0.974, 5: 0.943, 7: 0.925, 11: 0.902,
         13: 0.895, 17: 0.882}   # 141/145 ekran (son, 0.52-taban)
for p in PR:
    tau, A1, R0, R1 = out[np.log(p)]
    print(f"  p={p:>2} τ={tau:.3f}  R0={R0:+.3f}  R0·r={R0*RLINE[p]:+.3f}"
          f"   R1={R1:+.3f} (R1/cosπτ={R1/np.cos(np.pi*tau):+.3f})",
          flush=True)

# (i-b) boyalı çizgi-dışı: ad0(τ*)
gbar = TWO_PI / L
LINES = [np.log(q) for q in pk(720)]
print("\nBOYALI (ad0 = varyans-kanal adyabatik; ad1 bond kontrol):")
AD0 = []
for ts in [0.10, 0.16, 0.24]:
    om = ts * L
    while min(abs(om - l) for l in LINES) < 0.012:
        om += 0.013
    U = 0.1 / (2 * np.sin(np.pi * om / L))
    zp = np.sort(zz0 + U * gbar * np.cos(om * zz0))
    outp, _, _, _ = kanallar(zp, [om], ekstra=[om])
    tau, A1, R0, R1 = outp[om]
    AD0.append((tau, R0))
    print(f"  τ*={tau:.3f}  ad0={R0:+.3f}   ad1={R1:+.3f} "
          f"(ad1/cosπτ={R1/np.cos(np.pi*tau):+.3f})", flush=True)

# (ii-iii) tek-R çapraz bağlama (ρ-modeli 144: exp(2.62−6.20τ))
num_b = den_b = n0B = n0C = nbC = d0 = 0.0
for p in primerange(2, int(np.exp(L)) + 1):
    q, mmn = p, 1
    while q <= int(np.exp(L)):
        t = np.log(q) / L
        if t > 0.52:
            w = np.exp(2.62 - 6.20 * t) * np.sin(np.pi * t)**2 / (
                np.pi**2 * mmn**2 * q)
            num_b += w * np.pi * t * np.cos(3*np.pi*t) / np.sin(np.pi*t)
            den_b += w * np.cos(2*np.pi*t)
            nbC += w * (-TWO_PI * t) * np.cos(2*np.pi*t)
            n0B += w * TWO_PI * t / np.tan(np.pi*t)
            n0C += w * (-TWO_PI * t)
            d0 += w
        q *= p; mmn += 1
kB = num_b / den_b
K0 = -2.87
R = (K0 - kB) * den_b / nbC
print(f"\nkB={kB:+.3f}  K₀={K0}  →  R = {R:.3f}")
print(f"lag-0 çekirdek öngörüsü: κ_B0={n0B/d0:+.3f}  κ_C0={R*n0C/d0:+.3f}"
      f"  toplam={n0B/d0 + R*n0C/d0:+.3f}")
print("karşılaştırma: ölçülen lag-0 çekirdek = R0·r − ad0(interp):")
at, av = [t for t, _ in AD0], [v for _, v in AD0]
for p in PR:
    tau, A1, R0, R1 = out[np.log(p)]
    ker0 = R0 * RLINE[p] - float(np.interp(tau, at, av))
    print(f"  p={p:>2}: {ker0:+.3f}", flush=True)
