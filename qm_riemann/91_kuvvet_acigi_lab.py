"""
91 — KUVVET AÇIĞI LABORATUVARI: İKİNCİ-MERTEBE TÜRETİM (21 Ağustos 2026)
==========================================================================
Gözlem (85/90): kuvvet çizgileri tam Λ-ağırlığının altında (q,k ile
büyüyen açık); asallar tam. Aday mekanizma: sıfır koşulu ÖRTÜK —
   ρ̄·u_n = −S(t_n^s + u_n)   (sıfır, kaydığı yerdeki alanı hisseder)
⟹ ikinci mertebe u ≈ −S/ρ̄ + S'S/ρ̄²; S'S çarpımı p^i·p^j = p^k
kombinasyonlarıyla YALNIZ kuvvet çizgilerini besler (asalları p^{-1}
bastırmalı fark terimleri dışında beslemez) — asal-bağışıklığı doğal.

LAB: EF çizgilerinden (tüm p^k ≤ 720, bilinen genlikler, GÜRÜLTÜSÜZ)
sentetik örgü kur, örtük denklemi tam çöz (sönümlü sabit-nokta), gerçek
ölçüm makinesini aynen çalıştır:
  T1  benekler: ε_lab(q) = |Ĝ|_lab / [Λ/(L√q)·cos(πτ)·DW] — gerçek
      ε_ölçüm ile yan yana (asallar kontrol: ~1 kalmalı).
  T2  gap kanalı: R_lab = k·v_q/v_asal(τ_q) — 90-T2'nin lab kopyası.
  T3  ölçekleme: iki L'de (8.5, 10.37) ε'nin L ve (k,p) bağımlılığı;
      analitik aday yarışı: ε−1 ∝ k·H_{k-1}·τ_p ? / k·τ_p ? / τ_q ?

SONUÇ (21 Ağustos — dürüst muhasebe):
  1. LAB DESENİ ÜRETİYOR: kuvvetler asallardan sistematik derin, sıralama
     gerçekle uyumlu; gap-oranı R_lab ≈ 1.01-1.04 (90-T2'nin "oranda açık
     kaybolur" bulgusu lab'da da ✓). AMA lab açıkları abartıyor (gürültüsüz
     → maksimum koherans) VE fazları döndürüyor (−176° → −150°).
  2. ÖLDÜRME TESTİ: gerçek veride kuvvet-beneği fazları SAF 180°
     (0.1-0.6° içinde; asal kontrol 180.0°) — faz dönmesi YOK.
     → Koheran S'S mekanizması gerçek açığın kaynağı DEĞİL (gerçek
     alanın rastgele bileşeni koheran 2. mertebeyi dekohere ediyor).
  3. YENİ GÜÇLÜ KISIT: gerçek mekanizma FAZ-KORUYUCU (saf reel çarpan),
     aritmetik-seçici (yalnız p^k), q ve k ile büyüyen. Bessel/harmonik-
     örnekleme adayı da nicel elendi (J_4(0.38) küçük; 49-25 sıralaması
     ters). AÇIK: dekohere-filtreli harmonik terimler; S(t) çizgi
     genliklerinin 2. mertebe EF düzeltmeleri; kendi-çizgi × harmonik
     girişiminin dikkatli analitiği. T3 ölçeklemesi kararsız (düşük-L
     lab güvenilmez) — nicel yasa henüz yok.
"""

import numpy as np
from pathlib import Path
from sympy import primerange, factorint

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
P11 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
QPOW = [(4, 2, 2), (8, 2, 3), (9, 3, 2), (16, 2, 4),
        (25, 5, 2), (27, 3, 3), (32, 2, 5), (49, 7, 2)]
EPS_REAL = {4: 0.992, 8: 0.964, 9: 0.975, 16: 0.873,
            25: 0.893, 27: 0.806, 32: 0.652, 49: 0.723}

def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q)
            q *= p
    return sorted(set(out))

QS = pk_list(720)
LAM = {}
for q in QS:
    (pp, kk), = factorint(q).items()
    LAM[q] = (float(np.log(pp)), kk)

def S_field(t):
    s = np.zeros_like(t)
    for q in QS:
        lam, kk = LAM[q]
        s -= (lam / (np.pi * np.sqrt(q) * np.log(q))) * np.sin(t * np.log(q))
    return s

def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8

def Ghat(t, omegas, chunk=30000):
    out = np.zeros(len(omegas), dtype=complex)
    for s0 in range(0, len(t), chunk):
        tt = t[s0:s0 + chunk]
        out += np.exp(1j * np.outer(omegas, tt)).sum(axis=1)
    return out / len(t)

def lab_orgusu(L0, n=40000):
    t0 = TWO_PI * np.exp(L0)
    N0 = rvm_N(t0)
    idx = np.arange(n + 1, dtype=float)
    t = t0 + idx * TWO_PI / L0
    for _ in range(8):
        t = t - (rvm_N(t) - N0 - idx) / (np.log(t / TWO_PI) / TWO_PI)
    rho = np.log(t / TWO_PI) / TWO_PI
    u = np.zeros_like(t)
    for _ in range(60):
        u = 0.5 * u + 0.5 * (-S_field(t + u) / rho)
    z = t + u
    gaps = np.diff(z)
    mids = 0.5 * (z[:-1] + z[1:])
    return z, gaps, mids, float(np.log(mids / TWO_PI).mean()), (mids - 0.5*(t[:-1]+t[1:])).std()

def olc(L0):
    z, gaps, mids, L, sig = lab_orgusu(L0)
    sonuc = {}
    # benekler
    for q in P11 + [q for q, _, _ in QPOW]:
        lam, kk = LAM[q]
        om = np.log(q)
        tau = om / L
        if tau > 0.55:
            continue
        dw = np.exp(-om**2 * sig**2 / 2)
        pred = lam / (L * np.sqrt(q)) * np.cos(np.pi * tau) * dw
        G = Ghat(mids, np.array([om]))[0]
        sonuc[q] = (abs(G) / pred, G, tau, kk)
    # gap kanalı
    Lw = np.log(mids / TWO_PI)
    yg = np.log(gaps * Lw / TWO_PI)
    qs = [q for q in QS if np.log(q) / L <= 0.55]
    cols = [np.ones_like(yg)]
    for q in qs:
        arg = mids * np.log(q)
        cols += [np.cos(arg), np.sin(arg)]
    X = np.vstack(cols).T
    b, *_ = np.linalg.lstsq(X, yg, rcond=None)
    vdict = {}
    for q in qs:
        i = qs.index(q)
        vdict[q] = np.hypot(b[1 + 2*i], b[2 + 2*i]) / q**-0.5
    prim = np.array(sorted((np.log(p) / L, vdict[p]) for p in P11 if p in vdict))
    R = {}
    for q, pp, kk in QPOW:
        tau = np.log(q) / L
        if q in vdict and prim[0, 0] <= tau <= prim[-1, 0]:
            R[q] = kk * vdict[q] / np.interp(tau, prim[:, 0], prim[:, 1])
    return sonuc, R, L

print("LAB L≈10.37 (gerçek pencereyle aynı yükseklik):")
sp, Rlab, L = olc(10.37)
print(f"  örgü L = {L:.3f}")
print(f"  {'q':>3} {'ε_lab':>7} {'ε_ölçüm':>8} {'faz(Ĝ)':>8}")
for q in P11:
    if q in sp:
        r, G, tau, kk = sp[q]
        print(f"  {q:>3} {r:>7.3f} {'(asal)':>8} {np.degrees(np.angle(G)):>7.0f}°")
for q, pp, kk in QPOW:
    if q in sp:
        r, G, tau, _ = sp[q]
        print(f"  {q:>3} {r:>7.3f} {EPS_REAL[q]:>8.3f} {np.degrees(np.angle(G)):>7.0f}°")
print("  GAP KANALI R_lab = k·v/v_asal (90-T2 lab kopyası):")
for q, pp, kk in QPOW:
    if q in Rlab:
        print(f"    q={q:>2}: R_lab = {Rlab[q]:.3f}")

# ---- T3: ölçekleme ve analitik yarış
print("\nT3 — ÖLÇEKLEME (L=8.5 vs 10.37) ve analitik aday yarışı:")
sp2, _, L2 = olc(8.5)
H = {1: 0.0, 2: 1.0, 3: 1.5, 4: 11/6, 5: 25/12}
print(f"  {'q':>3} {'1−ε(10.4)':>9} {'1−ε(8.5)':>9} {'oran':>6} "
      f"{'k·H·τp':>7} {'k·τp':>6} {'τ_q':>6}")
for q, pp, kk in QPOW:
    if q in sp and q in sp2:
        d1 = 1 - sp[q][0]
        d2 = 1 - sp2[q][0]
        taup = np.log(pp) / L
        print(f"  {q:>3} {d1:>9.3f} {d2:>9.3f} "
              f"{(d2/d1 if abs(d1) > 1e-4 else float('nan')):>6.2f} "
              f"{kk*H[kk]*taup:>7.3f} {kk*taup:>6.3f} {np.log(q)/L:>6.3f}")
print(f"  (L-oranı beklentisi 1/L ölçeklemesiyle: {L/L2:.2f})")
