"""
101 — DONMUŞ KOYLARIN HARİTASI: ÇÖZÜLME SINIRI ADALARDA (20 Ağustos)
==========================================================================
93'ün keşfi (ζ): spontane koşullu gap-yanıtı ilk asal çizgisinin
altında ~0 (donmuş), çizgiler açıldıkça CUE'ya termalleşir.
HİPOTEZ: "gazı çözen aritmetiktir — donma ilk SAĞ KALAN çizgide biter."

ÖN-MÜHÜRLÜ ÖNGÖRÜLER (ölçümden önce):
  P1  χ₃, χ₅, χ₇: çözülme ζ gibi ilk çizgide (log2/L ≈ 0.068-0.072)
      başlar; τ ≥ 0.3'te CUE (1.05-1.18).
  P2  AYRIŞTIRICI — β'nın 2-çizgisi YOK: çözülme log3/L ≈ 0.11'e
      GECİKMELİ; τ = 0.08-0.10'da ötekiler ~0.5-0.65 iken β donuk
      (~0.1-0.3) kalmalı; yarı-çözülme noktası ~×1.6 kaymalı.
  P3  Dört ada τ ≥ 0.3'te TEK CUE eğrisinde (spontane dal evrensel).
  P4  ζ aynı-ayak tarak-termometresi (36/41 pencereleri) — H3 nicel
      tablosuna ζ satırları.

Makine: 93'ün D_spont estimatörü (yerel katlama + kübik detrend +
çizgi-dışı bantlar), her adanın EN BÜYÜK penceresi (L≈9.5-10.1).

SONUÇ (20 Ağustos — HAM SONUÇ ARTEFAKT, ölçüm 101d'de):
  Dört ada τ=0.04'te D≈0.96-0.98 verdi (ζ: 0.041) — "hiç donma yok"
  gibi. 101b kusur kapısı bunu ÇÜRÜTTÜ: motorun ~%0.3 kaçık yakın-çifti
  örgü kusurudur; temiz ζ'ya %0.1 silme enjeksiyonu bile D'yi
  0.09→0.55'e fırlatır (%0.3: 0.76). Donma estimatörü kinematik kilidin
  (yapısız süreçte D→1/cos(κ/2); Poisson kontrolü 1.016) İHLALİNİ ölçer
  ve hiperuniform karanlık alanda kusura aşırı duyarlıdır. Benek/faz
  ölçümleri bu kusura bağışık (98-S2), donma ölçümü DEĞİL.
  → Sertifikalı sıfırlar: 101c; gerçek ada-donma ölçümü: 101d.
  P4 (ζ aynı-ayak tarak) TEMİZ VERİYLE GEÇERLİ: σ_u = 0.229 (L=6.99)
  → 0.283 (L=12.45), 10 pencere — H3 tablosunun ζ satırları.
"""

import numpy as np
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

LINES_ALL = [np.log(q) for q in pk_list(200)]

def D_spont(zz, q_cond, tau_grid):
    """En büyük log-üçte-birlik pencerede, çizgi-dışı bant estimatörü."""
    edges = np.exp(np.linspace(np.log(zz[0] + 1), np.log(zz[-1]), 4))
    m = zz >= edges[2]
    z = zz[m]
    gaps = np.diff(z)
    mids = 0.5 * (z[:-1] + z[1:])
    L = float(np.log(q_cond * mids / TWO_PI).mean())
    ds = gaps * np.log(q_cond * mids / TWO_PI) / TWO_PI - 1
    tt = (mids - mids.mean()) / (mids[-1] - mids[0])
    P = np.vstack([np.ones_like(tt), tt, tt**2, tt**3]).T
    ds = ds - P @ np.linalg.lstsq(P, ds, rcond=None)[0]
    out = []
    for tau0 in tau_grid:
        oms = tau0 * L + np.linspace(-0.018 * L, 0.018 * L, 140)
        oms = np.array([o for o in oms
                        if min(abs(o - l) for l in LINES_ALL) > 0.01])
        num = 0.0 + 0j; den = 0.0
        for s0 in range(0, len(oms), 35):
            ob = oms[s0:s0 + 35]
            rr = np.exp(1j * np.outer(ob, z)).sum(axis=1)
            GG = (np.exp(1j * np.outer(ob, mids)) * ds[None, :]).sum(axis=1)
            num += (GG * np.conj(rr)).sum()
            den += (np.abs(rr)**2).sum()
        kap = 2 * np.pi * tau0
        out.append(abs(num) / ((2 * np.sin(kap / 2) / kap) * den))
    return L, out

TAUG = [0.04, 0.06, 0.08, 0.10, 0.125, 0.15, 0.20, 0.30, 0.40, 0.50]
print("ÇÖZÜLME EĞRİLERİ (en büyük pencere; CUE referans: "
      "1.00/1.01/1.04/1.08/1.13/1.18 @ 0.03-0.5):")
print(f"{'τ':>6} " + " ".join(f"{s:>7}" for s in ["chi3", "beta", "chi5", "chi7"]))
CURVES = {}
for etiket, q in [("chi3", 3), ("beta", 4), ("chi5", 5), ("chi7", 7)]:
    zz = np.load(HERE / f"99_{etiket}_zeros.npz")["zeros"]
    L, D = D_spont(zz, q, TAUG)
    CURVES[etiket] = (L, D)
    print(f"  [{etiket}: L={L:.2f}, ilk-çizgi τ = "
          f"{np.log(3 if etiket == 'beta' else 2)/L:.3f}]")
for i, tau0 in enumerate(TAUG):
    print(f"{tau0:>6.3f} " + " ".join(f"{CURVES[e][1][i]:>7.3f}"
                                       for e in ["chi3", "beta", "chi5", "chi7"]))

# ---- P4: ζ aynı-ayak tarak-termometresi
print("\nP4 — ζ AYNI-AYAK TARAK (36+41 pencereleri):")
def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8
d36 = np.load(HERE / "36_T100k.npz")
d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
WND = []
edges = np.geomspace(d36["t_mid"][0], d36["t_mid"][-1] * 1.0001, 13)
for i in range(12):
    m = (d36["t_mid"] >= edges[i]) & (d36["t_mid"] < edges[i + 1])
    if m.sum() >= 3000:
        WND.append((d36["t_mid"][m], "36"))
for k in K41[:6]:
    WND.append((d41[f"tmid_{k}"], "41"))
for tmid, src in WND:
    L = float(np.log(tmid / TWO_PI).mean())
    x = rvm_N(tmid)
    comb = abs(np.exp(2j * np.pi * x).mean())
    su = np.sqrt(-2 * np.log(comb)) / TWO_PI
    print(f"  ζ[{src}] L={L:5.2f}: σ_u(tarak) = {su:.4f}")
