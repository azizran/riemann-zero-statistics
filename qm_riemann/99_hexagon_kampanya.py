"""
99 — HEXAGON KAMPANYASI: MOTOR İLE DÖRT ADA, BÜYÜK ÖLÇEK (20 Ağustos)
==========================================================================
98 motoruyla takımada kampanyası: dört ada × ~60-75k sıfır (mpmath
seferlerinin ~15 katı) + her adada 3 iç-pencere (ada-içi ısınma eğrisi,
aynı-ayak termometreyle koro-yasası testi).

ÖN-MÜHÜRLÜ ÖNGÖRÜLER (ölçümden önce):
  H1  MOD 7 ALTIGEN KADRANI (χ(3)=e^{iπ/3}, tek; 6. birim kökleri):
      faz = 180° + arg χ (çapa χ=+1→180°):
        n≡1→180° (29) | n≡3→240° (3,17) | n≡2→300° (2,23) |
        n≡6→0° (13) | n≡4→60° (4,11,25) | n≡5→120° (5,19)
      ALTI konum — kompleks spektroskopinin tam kadranı. Ölü: 7, 49.
  H2  Mutlak yasa dört adada, 15× istatistikle (benek hatası ~3× küçülür).
  H3  KORO-YASASI ada-içi: her adanın 3 penceresi kendi ısınma eğrisini
      çizer; koro-düzeltmeli "tam" sıcaklıklar dört adada TEK eğriye
      binmeli (97'nin ±%8'i ada-içi eğimlerle sınanır).
  H4  Kuvvet açığı 4 adada, kesin ölçekleme: ε(q,k) vs τ_q ve k
      (motor istatistiğiyle ilk nicel yasa denemesi).
Motor kusuru notu: ~%0.3 sığ çift kaçar → tarak ~%1 sistematik; benek/
faz ölçümleri duyarsız (98 sağlaması).
"""

import numpy as np
import time
from pathlib import Path

exec(open("98_L_motoru.py").read().split('CHI4 = ')[0])

HERE = Path(".").resolve()
CHI4 = {0: 0, 1: 1, 2: 0, 3: -1}
CHI3 = {0: 0, 1: 1, 2: -1}
CHI5 = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
z6 = np.exp(1j * np.pi / 3)
CHI7 = {0: 0, 1: 1, 3: z6, 2: z6**2, 6: z6**3, 4: z6**4, 5: z6**5}

ISLANDS = [
    ("chi3", 3, CHI3, 55000.0,
     [(3, np.log(3), 1), (9, np.log(3), 2), (2, np.log(2), 1),
      (5, np.log(5), 1), (8, np.log(2), 3), (11, np.log(11), 1),
      (4, np.log(2), 2), (7, np.log(7), 1), (13, np.log(13), 1),
      (25, np.log(5), 2), (27, np.log(3), 3)]),
    ("beta", 4, CHI4, 50000.0,
     [(2, np.log(2), 1), (4, np.log(2), 2), (8, np.log(2), 3),
      (3, np.log(3), 1), (5, np.log(5), 1), (7, np.log(7), 1),
      (9, np.log(3), 2), (11, np.log(11), 1), (13, np.log(13), 1),
      (25, np.log(5), 2), (27, np.log(3), 3)]),
    ("chi5", 5, CHI5, 48000.0,
     [(5, np.log(5), 1), (25, np.log(5), 2), (2, np.log(2), 1),
      (7, np.log(7), 1), (17, np.log(17), 1), (3, np.log(3), 1),
      (8, np.log(2), 3), (13, np.log(13), 1), (23, np.log(23), 1),
      (4, np.log(2), 2), (9, np.log(3), 2), (19, np.log(19), 1),
      (11, np.log(11), 1), (27, np.log(3), 3)]),
    ("chi7", 7, CHI7, 42000.0,
     [(7, np.log(7), 1), (49, np.log(7), 2), (29, np.log(29), 1),
      (3, np.log(3), 1), (17, np.log(17), 1), (2, np.log(2), 1),
      (23, np.log(23), 1), (13, np.log(13), 1), (4, np.log(2), 2),
      (11, np.log(11), 1), (25, np.log(5), 2), (5, np.log(5), 1),
      (19, np.log(19), 1), (9, np.log(3), 2), (8, np.log(2), 3),
      (27, np.log(3), 3)]),
]

def Ghat(t, omegas, chunk=30000):
    out = np.zeros(len(omegas), dtype=complex)
    for s0 in range(0, len(t), chunk):
        tt = t[s0:s0 + chunk]
        out += np.exp(1j * np.outer(omegas, tt)).sum(axis=1)
    return out / len(t)

SUM = {}
for etiket, q, tab, T1, LAMK in ISLANDS:
    cache = HERE / f"99_{etiket}_zeros.npz"
    M = Lmotor(q, tab, 1)
    if cache.exists():
        zz = np.load(cache)["zeros"]
    else:
        t0 = time.time()
        zz = M.sifir_bul(200.0, T1, grid_frac=0.12)
        np.savez(cache, zeros=zz)
        print(f"[{etiket}] {len(zz)} sıfır, {time.time()-t0:.0f} sn")
    gaps = np.diff(zz)
    mids = 0.5 * (zz[:-1] + zz[1:])
    n = len(mids)
    Leff = float(np.log(q * mids / TWO_PI).mean())

    def Nsm(t):
        x = t / TWO_PI
        return x * np.log(q * x / np.e)

    # ada-içi 3 pencere: log-t eşit bölme; her birinde tarak-termometre
    print(f"\n===== {etiket} (q={q}): n={n}, L_eff={Leff:.3f} =====")
    pens = []
    edges = np.exp(np.linspace(np.log(zz[0] + 1), np.log(zz[-1]), 4))
    for i in range(3):
        m = (mids >= edges[i]) & (mids < edges[i + 1])
        mm, zzw = mids[m], mids[m]
        Lw = float(np.log(q * mm / TWO_PI).mean())
        xw = Nsm(mm) - Nsm(mm[0])
        comb = abs(np.exp(2j * np.pi * xw).mean())
        su = np.sqrt(-2 * np.log(comb)) / TWO_PI
        pens.append((Lw, su, int(m.sum())))
        print(f"  pencere {i+1}: L={Lw:.2f}  σ_u = {su:.4f}  (n={m.sum()})")
    x_all = Nsm(mids) - Nsm(zz[0])
    comb_all = abs(np.exp(2j * np.pi * x_all).mean())
    sig_u = np.sqrt(-2 * np.log(comb_all)) / TWO_PI
    sig_t = sig_u * TWO_PI / np.log(q * mids.mean() / TWO_PI)

    print(f"  {'q':>3} {'sınıf':>6} {'ölçüm':>8} {'öngörü':>8} {'oran':>6} "
          f"{'faz':>7} {'ö-faz':>6}")
    rows = []
    for qq, lam, kk in LAMK:
        om = np.log(qq)
        tau = om / Leff
        G = Ghat(mids, np.array([om]))[0]
        chi = tab[qq % q]
        ph = (np.degrees(np.angle(G)) + 360) % 360
        if chi == 0:
            print(f"  {qq:>3} {'ölü':>6} {abs(G):>8.4f} {'—':>8} {'—':>6} "
                  f"{ph:>6.0f}° {'—':>6}")
            rows.append((qq, kk, tau, abs(G), np.nan, ph, np.nan))
        else:
            dw = np.exp(-om**2 * sig_t**2 / 2)
            pred = lam / (Leff * np.sqrt(qq)) * np.cos(np.pi * tau) * dw
            pph = (180 + np.degrees(np.angle(complex(chi)))) % 360
            print(f"  {qq:>3} {qq % q:>6} {abs(G):>8.4f} {pred:>8.4f} "
                  f"{abs(G)/pred:>6.3f} {ph:>6.1f}° {pph:>5.0f}°")
            rows.append((qq, kk, tau, abs(G), abs(G)/pred, ph, pph))
    SUM[etiket] = dict(L=Leff, sig_u=sig_u, pens=pens, rows=rows, n=n)

print("\n===== H3: ADA-İÇİ ISINMA + KORO TABLOSU =====")
UFAM = {"chi3": [3, 9, 27], "beta": [2, 4, 8], "chi5": [5, 25], "chi7": [7, 49]}
for etiket in SUM:
    d = SUM[etiket]
    q = {"chi3": 3, "beta": 4, "chi5": 5, "chi7": 7}[etiket]
    print(f"{etiket}: L={d['L']:.2f}  σ_u(tüm) = {d['sig_u']:.4f}  "
          f"pencereler: " + " ".join(f"({L:.1f},{s:.3f})" for L, s, _ in d['pens']))
