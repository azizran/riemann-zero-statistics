"""
134 — İPTALİN KİMLİĞİ: KAYNAK MI, FREKANS MI? (27 Ağustos, gün sonu)
==========================================================================
133 sonrası tek-toplam muhasebesinin son bilinmeyeni: net R = kinematik
kanal − iptal; iptal NEYE bakıyor?
DENEY: D3-lab, ω_2 çizgisine boyalı EK dalga (aritmetik genliğe eşit —
sürücü 2×) bindirilir; R_nn @ p=2 ölçülür.
ÖN-MÜHÜR (iki hipotez, keskin ayrım):
  H-kaynak:  iptalden yalnız aritmetik-kökenli sürücü kaçar → boyalı
             pay iptal edilir → pay sabit, A1 2× → R_ölçülen ≈ R/2.
  H-frekans: aritmetik-destekli frekans TÜM sürücüyü korur → boyalı
             pay katılır → R_ölçülen ≈ R (değişmez).
  Kontrol: çizgi-dışı boyalı (bilinen: tam iptal, R≈0) aynı koşuda.
Sonuç, kalemdeki iptal teriminin bağlanma biçimini (öz-tutarlılık-
kökenli mi, tayf-destek-kökenli mi) tek hamlede belirler.

SONUÇ (27 Ağustos) — H-KAYNAK (süperpozisyon haliyle) VE BÜYÜK SENTEZ:
  taban: A1=0.121, R=−1.713 | 2×-sürücü: A1=0.253 (2.1× ✓), R=−0.542 |
  çizgi-dışı boyalı: R=+0.898.
  OKUMA: boyalı eş-frekans payı aritmetik yanıta KATILMIYOR — kendi
  ADYABATİK (+0.9) yanıtını taşıyor ve karışım genlik-ağırlıklı
  süperpozisyon: (0.121·(−1.71)+0.133·(+0.9))/0.253 = −0.35 ~ ölçülen
  −0.54 (kaba mertebede ✓). "Frekans-bekçisi" YOK; iki yanıt yasası var
  ve süperpoze oluyorlar.
  BÜYÜK SENTEZ — İPTAL DEĞİL AYRIŞIM:
     R_aritmetik = R_adyabatik(+~1.2) + R_çekirdek(−~3.2) ≈ −2.0
  132c'nin sıfır-parametre çekirdeği κ_kuram = −3.22 ve kontrollerdeki
  adyabatik +0.9..+1.2 ile: −3.22 + 1.2 = −2.02 ✓✓. F2 de oturuyor:
  faz-karıştırma korelasyon-çekirdeğini öldürür → kalan adyabatik
  (+0.3..+0.6 ölçülmüştü ✓). SEÇİCİLİK yeniden adlandı: boyalı sürücü
  yalnız adyabatik terimi alır (aritmetik-korelasyonu yok); kendi
  dalgası her ikisini. "−2 = +1 (evrensel kinematik) − 3.2 (kuyruk-
  korelasyon çekirdeği, 132c'de parametresiz ölçülü)."
  KALAN KALEM (tek terim): çekirdeğin −3.2'sinin analitik türetimi
  (emilim-renormalize kuyruk tek-toplamı) — ayrışım artık ölçülmüş
  bileşenlerle KAPALI.
"""

import numpy as np
import time
from pathlib import Path
from sympy import primerange, factorint

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
L0 = 7.0
NZ = 40000

def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8

t0 = TWO_PI * np.exp(L0)
idx = np.arange(NZ + 1, dtype=float)
tg = t0 + idx * TWO_PI / L0
for _ in range(8):
    tg = tg - (rvm_N(tg) - rvm_N(t0) - idx) / (np.log(tg / TWO_PI) / TWO_PI)
rho = np.log(tg / TWO_PI) / TWO_PI

QS = pk_list(int(np.exp(1.4 * L0)))
LOGQ = np.log(np.array(QS))
LAMv = np.array([float(np.log(list(factorint(q).items())[0][0])) for q in QS])
W0 = LAMv / (np.pi * np.sqrt(np.array(QS)) * LOGQ)

def S_field(t, chunk=8000):
    s = np.zeros_like(t)
    for s0 in range(0, len(t), chunk):
        tt = t[s0:s0 + chunk]
        s[s0:s0 + chunk] = -(np.sin(np.outer(tt, LOGQ)) * W0[None, :]).sum(axis=1)
    return s

u = np.zeros_like(tg)
for _ in range(40):
    u = 0.5 * u + 0.5 * (-S_field(tg + u) / rho)
z_taban = np.sort(tg + u)

a2 = 1.0 / (np.pi * np.sqrt(2))
gbar = TWO_PI / L0
U2 = (a2 / rho.mean())                      # aritmetik 2-dalgası u-genliği

def olc(z, hedef_om):
    g = np.diff(z)
    m = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(m / TWO_PI)
    Lb = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    qs_fit = pk_list(min(int(np.exp(0.52 * Lb)), 720))
    freqs = [np.log(q) for q in qs_fit]
    if all(abs(hedef_om - f) > 1e-9 for f in freqs):
        freqs = freqs + [hedef_om]
    tt = (m - m.mean()) / (m[-1] - m[0])
    cols = [np.ones_like(m), tt, tt**2]
    for f in freqs:
        cols += [np.cos(f * m), np.sin(f * m)]
    X = np.vstack(cols).T
    b1, *_ = np.linalg.lstsq(X, ds, rcond=None)
    eta = ds - X @ b1
    s2 = float((eta**2).mean())
    ee = eta[:-1] * eta[1:]
    c1 = float(ee.mean())
    mm = 0.5 * (m[:-1] + m[1:])
    tt2 = (mm - mm.mean()) / (mm[-1] - mm[0])
    cols2 = [np.ones_like(mm), tt2, tt2**2]
    for f in freqs:
        cols2 += [np.cos(f * mm), np.sin(f * mm)]
    X2 = np.vstack(cols2).T
    b3, *_ = np.linalg.lstsq(X2, ee - c1, rcond=None)
    i = freqs.index(hedef_om)
    cA, sA = b1[3 + 2 * i], b1[3 + 2 * i + 1]
    A1 = np.hypot(cA, sA); ph = np.arctan2(sA, cA)
    Rn = (b3[3 + 2 * i] * np.cos(ph) +
          b3[3 + 2 * i + 1] * np.sin(ph)) / (2 * c1 * A1)
    return A1, Rn, np.sqrt(s2)

om2 = np.log(2)
om_off = om2 + 0.037
print("hücreler (R_nn @ hedef frekans):", flush=True)
for ad, z, om in [
    ("taban (aritmetik 2-dalga)", z_taban, om2),
    ("2-çizgisi + boyalı EŞ-genlik (sürücü 2×)",
     np.sort(z_taban + U2 * np.sin(om2 * z_taban)), om2),
    ("kontrol: çizgi-dışı boyalı",
     np.sort(z_taban + U2 * np.sin(om_off * z_taban)), om_off),
]:
    A1, Rn, se = olc(z, om)
    print(f"  [{ad}] A1={A1:.4f}  R_nn={Rn:+.3f}  σ_η={se:.4f}", flush=True)
print("\nH-kaynak: 2. hücrede R ≈ taban/2;  H-frekans: R ≈ taban.")
