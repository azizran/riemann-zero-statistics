"""
108c — 108b ZİNCİRİNİN ENJEKSİYON KALİBRASYONU (24 Ağustos)
==========================================================================
108b bulgusu: sinüs-normda perde 0.955→0.823 nazik iniş, 8 ailede
EVRENSEL, sıcaklıktan bağımsız (P2 RET, P3 isabet). ŞÜPHE: inişin
tamamı boru-hattı transferi (a_v) olabilir; fiziksel perde SABİT ~0.95.

ÖN-MÜHÜR: a_v(τ) bu zincir için ölçülür (bilinen genlikli yapay
yerdeğiştirme dalgası enjeksiyonu, çizgi-dışı frekans, tam-taban
ortak regresyonla geri-okuma; ζ-120k + chi3 evrensellik çaprazı).
  E1  a_v(τ) düşüyorsa ve D_sin/a_v ≈ SABİT (0.94-0.98) çıkarsa:
      perde τ-bağımsız sabit tutmadır — teori hedefi kökten değişir
      (w₀<1'in %2-5 kalıcı tutmasıyla birleşme adayı).
  E2  D_sin/a_v hâlâ τ ile iniyorsa: gerçek τ-bağımlı perde var.

SONUÇ (24 Ağustos) — E2 + BONUS SEÇİCİLİK HÜKMÜ:
  a_v(τ) ≈ 1.00 (0.997-1.033; ζ ve χ₃ özdeş) → tam-taban zinciri
  TEMİZ ALET; iniş boru-hattı değil. D_fiz = 0.958 → 0.797 (τ 0.1→0.5),
  8 ailede evrensel, sıcaklık-bağımsız, yaklaşık τ-lineer (≈1−0.36τ).
  BONUS (kalibrasyonun kendisi bir fizik deneyi çıktı): enjekte edilen
  yapay dalga AYNI frekans ve aletle PERDELENMİYOR (a_v≈1) — perde
  yalnız gazın KENDİ (aritmetik) dalgalarını bastırır. Yarım-gap-jitter
  DW adayı da böylece ölür (yapay dalgayı da bastırırdı). Kalan tek
  mekanizma sınıfı: gap-dalgalanmalarının aritmetik dalganın FAZIYLA
  korelasyonu — gazın dalga çevresinde yeniden-dengelenmesi (92'nin
  ilk sezgisi, artık üç kontrolle izole).
"""

import numpy as np
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
Q_CAP = 720

def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

LINES = [np.log(q) for q in pk_list(300)]

def chunked_v_ekstra(yg, tmid, qs, om_x, chunk=40000):
    """Tam taban + enjeksiyon frekansı sütun çifti; b döner."""
    tt = (tmid - tmid.mean()) / (tmid[-1] - tmid[0])
    C = 3 + 2 * len(qs) + 2
    XtX = np.zeros((C, C)); Xty = np.zeros(C)
    n = len(yg)
    def cols(sl):
        c = [np.ones(sl.stop - sl.start), tt[sl], tt[sl]**2]
        for q in qs:
            arg = tmid[sl] * np.log(q)
            c += [np.cos(arg), np.sin(arg)]
        c += [np.cos(om_x * tmid[sl]), np.sin(om_x * tmid[sl])]
        return np.vstack(c).T
    for s0 in range(0, n, chunk):
        sl = slice(s0, min(s0 + chunk, n))
        Xc = cols(sl)
        XtX += Xc.T @ Xc; Xty += Xc.T @ yg[sl]
    b = np.linalg.solve(XtX, Xty)
    return b[-2], b[-1]

def kalibre(zz, qeff, taustars, U=0.02):
    g0 = np.diff(zz)
    m0 = 0.5 * (zz[:-1] + zz[1:])
    L = float(np.log(qeff * m0 / TWO_PI).mean())
    gbar = TWO_PI / L
    qs = pk_list(min(np.exp(0.52 * L), Q_CAP))
    out = []
    for ts in taustars:
        om = ts * L
        om += 0.013 if min(abs(om - l) for l in LINES) < 0.012 else 0.0
        zp = zz + U * gbar * np.cos(om * zz)
        gp = np.diff(zp)
        mp_ = 0.5 * (zp[:-1] + zp[1:])
        Lw = np.log(qeff * mp_ / TWO_PI)
        yg = np.log(gp * Lw / TWO_PI)
        bc, bs = chunked_v_ekstra(yg, mp_, qs, om)
        # taban (enjeksiyonsuz) aynı frekansta
        Lw0 = np.log(qeff * m0 / TWO_PI)
        yg0 = np.log(g0 * Lw0 / TWO_PI)
        bc0, bs0 = chunked_v_ekstra(yg0, m0, qs, om)
        v_ana = np.hypot(bc - bc0, bs - bs0)
        beklenen = U * gbar * (2 / gbar) * np.sin(np.pi * om / L)
        out.append((om / L, v_ana / beklenen))
    return out

d41 = np.load(HERE / "41_bigT_windows.npz")
gz, tm = d41["gaps_120k"], d41["tmid_120k"]
zz = np.empty(len(gz) + 1)
zz[0] = tm[0] - gz[0] / 2
zz[1:] = zz[0] + np.cumsum(gz)
zc = np.load(HERE / "101f_chi3_zeros.npz")["zeros"]
lo = np.exp(np.log(zc[0] + 1) + 0.25 * (np.log(zc[-1]) - np.log(zc[0] + 1)))
zc = zc[zc >= lo]

TS = [0.08, 0.15, 0.25, 0.35, 0.45, 0.50]
print("a_v(τ) — ζ-120k ve χ₃, tam-taban zinciri:", flush=True)
AZ = kalibre(zz, 1.0, TS)
AC = kalibre(zc, 3.0, TS)
# 108b'nin ζ D_sin bin değerleri (asal-bin ortalamaları, buradan böl)
D_ZETA = [(0.10, 0.955), (0.20, 0.913), (0.30, 0.883), (0.40, 0.852),
          (0.50, 0.823)]
print(f"{'τ':>6} {'a_v(ζ)':>7} {'a_v(χ₃)':>8}")
for (t1, a1), (t2, a2) in zip(AZ, AC):
    print(f"{t1:>6.3f} {a1:>7.3f} {a2:>8.3f}")
print(f"\n{'τ':>6} {'D_sin(ζ,108b)':>13} {'a_v(interp)':>11} {'D_fiz':>7}")
ts_a = [t for t, a in AZ]; av_a = [a for t, a in AZ]
for t, d in D_ZETA:
    a = np.interp(t, ts_a, av_a)
    print(f"{t:>6.2f} {d:>13.3f} {a:>11.3f} {d/a:>7.3f}")
