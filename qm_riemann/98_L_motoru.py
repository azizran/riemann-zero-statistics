"""
98 — VEKTÖRİZE L-MOTORU: DÖKÜM VE SAĞLAMA (20 Ağustos 2026)
==========================================================================
Hurwitz-mpmath 4.600 sıfıra 10 saat istiyordu; kampanya büyüyecekse
motor gerek. Tasarım: tamamlanmış-fonksiyon ana toplamı (AFE, simetrik
kesim N = √(qt/2π)), karakter fazları gömülü:
   Z_χ(t) = 2 Σ_{n≤N} |χ(n)| n^{-1/2} cos(θ(t) − t log n + arg χ(n))
   θ(t) = (t/2)log(q/π) + Im logΓ((s+a̸)/2) − arg(ε)/2
   (Im logΓ vektörize Stirling; a=1 tek, a=0 çift karakter)
Bilinen eksik: RS-düzeltme terimsiz → hata O((qt)^{-1/4}) — kabul
edilebilirliği SAĞLAMAYLA ölçülür, varsayılmaz:
  S1  θ-fazı vs mpmath (örneklem noktalarında)
  S2  motor-Z, mpmath-üretimi 3 adanın (β, χ₃, χ₅) önbellek sıfırlarında
      ~0 vermeli (medyan |Z|_sıfır / RMS|Z|_ara raporu)
  S3  yeniden-üretim: önbellek aralığında motorla sıfır bul, birebir
      eşle: maks |Δγ|/ort-boşluk raporu
  S4  hız: sıfır/saniye ölçümü (mpmath ~0.13 sıfır/sn idi)
"""

import numpy as np
import mpmath as mp
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
mp.mp.dps = 15

def stirling_imlog_gamma(z):
    """Im log Γ(z), z kompleks numpy dizisi, |z| büyük (Stirling)."""
    return ((z - 0.5) * np.log(z) - z + 0.5 * np.log(TWO_PI)
            + 1.0 / (12 * z) - 1.0 / (360 * z**3) + 1.0 / (1260 * z**5)).imag

class Lmotor:
    def __init__(self, q, chi_table, a_par):
        self.q = q
        self.a = a_par           # 0 çift, 1 tek
        nmax = 4096
        ns = np.arange(1, nmax + 1)
        chiv = np.array([chi_table[n % q] for n in ns], dtype=complex)
        self.abschi = np.abs(chiv)
        self.phichi = np.angle(chiv)
        self.logn = np.log(ns)
        self.wn = ns**-0.5
        tau = mp.mpc(0)
        for a2 in range(1, q):
            tau += chi_table[a2] * mp.e**(2j * mp.pi * a2 / q)
        eps = tau / ((1j if a_par == 1 else 1) * mp.sqrt(q))
        self.arg_eps = float(mp.arg(eps))
        self.abs_eps = abs(complex(eps))

    def theta(self, t):
        z = (0.5 + self.a) / 2 + 0.25 + 0j + 0.5j * t   # ((1/2+it)+a)/2
        z = ( (0.5 + 1j*t) + self.a ) / 2
        return (t / 2) * np.log(self.q / np.pi) + stirling_imlog_gamma(z) \
               - self.arg_eps / 2

    def Z(self, t, chunk=4000):
        t = np.atleast_1d(np.asarray(t, dtype=float))
        out = np.empty_like(t)
        N = np.sqrt(self.q * t / TWO_PI).astype(np.int64)
        th = self.theta(t)
        for Nv in np.unique(N):
            m = N == Nv
            tm, thm = t[m], th[m]
            zvals = np.empty(len(tm))
            for s0 in range(0, len(tm), chunk):
                sl = slice(s0, s0 + chunk)
                ph = (thm[sl, None] - tm[sl, None] * self.logn[None, :Nv]
                      + self.phichi[None, :Nv])
                zvals[sl] = 2 * (np.cos(ph) * (self.abschi[:Nv] * self.wn[:Nv])[None, :]).sum(axis=1)
            out[m] = zvals
        return out

    def sifir_bul(self, T0, T1, grid_frac=0.15):
        ts = [T0]
        t = T0
        while t < T1:
            t += grid_frac * TWO_PI / np.log(self.q * t / TWO_PI)
            ts.append(t)
        ts = np.array(ts)
        v = self.Z(ts)
        sc = np.where(np.sign(v[:-1]) * np.sign(v[1:]) < 0)[0]
        a, b = ts[sc].copy(), ts[sc + 1].copy()
        fa = v[sc].copy()
        for _ in range(30):
            mmid = 0.5 * (a + b)
            fm = self.Z(mmid)
            left = fa * fm <= 0
            b = np.where(left, mmid, b)
            a = np.where(left, a, mmid)
            fa = np.where(left, fa, fm)
        return 0.5 * (a + b)

CHI4 = {0: 0, 1: 1, 2: 0, 3: -1}
CHI3 = {0: 0, 1: 1, 2: -1}
CHI5 = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}

print("S1 — θ-fazı vs mpmath (t = 500, 2000, 4000):")
for isim, q, tab in [("β(χ₄)", 4, CHI4), ("χ₃", 3, CHI3), ("χ₅", 5, CHI5)]:
    M = Lmotor(q, tab, 1)
    for t0 in [500.0, 2000.0, 4000.0]:
        s = mp.mpc(0.5, t0)
        ps_mp = float((t0/2) * mp.log(q/mp.pi) + mp.im(mp.loggamma((s+1)/2))) \
                - M.arg_eps / 2
        ps_np = float(M.theta(np.array([t0]))[0])
        print(f"  {isim} t={t0:6.0f}: Δθ = {abs(ps_mp-ps_np):.2e}  |ε|={M.abs_eps:.6f}")

print("\nS2/S3 — önbellek sıfırlarında motor-Z + yeniden-üretim eşlemesi:")
import time
for isim, q, tab, cache in [("β(χ₄)", 4, CHI4, "96_beta_zeros.npz"),
                            ("χ₃", 3, CHI3, "97_chi3_zeros.npz"),
                            ("χ₅", 5, CHI5, "97_chi5_zeros.npz")]:
    zz = np.load(HERE / cache)["zeros"]
    M = Lmotor(q, tab, 1)
    Zat = np.abs(M.Z(zz))
    mids = 0.5 * (zz[:-1] + zz[1:])
    Zmid = np.abs(M.Z(mids))
    t0 = time.time()
    yeniden = M.sifir_bul(float(zz[0]) - 0.3, float(zz[-1]) + 0.3)
    dt = time.time() - t0
    # eşleme: her önbellek sıfırına en yakın motor sıfırı
    idx = np.searchsorted(yeniden, zz)
    idx = np.clip(idx, 1, len(yeniden) - 1)
    yakin = np.where(np.abs(yeniden[idx] - zz) < np.abs(yeniden[idx-1] - zz),
                     yeniden[idx], yeniden[idx-1])
    d = np.abs(yakin - zz)
    gbar = np.diff(zz).mean()
    print(f"  {isim}: |Z|@sıfır medyan {np.median(Zat):.2e} vs @ara {np.median(Zmid):.3f} "
          f"(oran {np.median(Zat)/np.median(Zmid):.1e})")
    print(f"        motor {len(yeniden)} sıfır buldu (önbellek {len(zz)}); "
          f"eşleme: medyan |Δγ|/⟨g⟩ = {np.median(d)/gbar:.2e}, "
          f"maks = {d.max()/gbar:.2e}; hız = {len(yeniden)/dt:.0f} sıfır/sn "
          f"(mpmath ~0.13/sn idi)")
