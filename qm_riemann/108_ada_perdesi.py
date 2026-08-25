"""
108 — ADA PERDESİ: SICAKLIK-SIRALAMASI TESTİ (24 Ağustos, kaptan kuyusu)
==========================================================================
Perde-teorisi hipotezi (106 sonrası): D(τ) = fark-kanalı DW'si —
komşu sıfırların GÖRELİ artık-kıpırdanması gap-kanalının koheran
okumasını söndürür: ln D ≈ −κ²σ_rel²/2, σ_rel² = 2σ²(1−ρ₁).
105'in soğukluk yasası kadran verdi: ölü-2 üçlüsü (β, χ₈ᵉ, χ₈ᵒ)
σ_tam² ≈ 0.050, canlı-2 dörtlüsü ≈ 0.058-0.064 (−%20).

ÖN-MÜHÜRLÜ ÖNGÖRÜLER (ölçümden önce):
  P1  Her ada perdeler: D(τ_q) < 1, τ ile düşer (ζ deseni).
  P2  MEKANİZMA TESTİ: aynı τ'da ölü-2 üçlüsü AZ perdeler;
      ln-bastırma oranı ≈ σ² oranı ≈ 0.80±0.1 → τ≈0.35'te D farkı
      ~%5-7, τ≈0.45'te ~%9; SIRALAMA σ_tam²'yi izler.
  P3  Kervan-içi çökme: her sıcaklık sınıfının çizgi-başına D'leri
      tek D(τ) eğrisine biner.
Estimatör notu: lineer ds kanalı, çizgi-frekansında cos/sin regresyonu;
katı payda v_rigid = 2sin(πτ)/(π√p) (Λ/log p = 1, asallar). Boru-hattı
çarpanı a_v mutlak D'yi ~1/0.85'e dek şişirebilir (92-T1a) — adalar-arası
KARŞILAŞTIRMA ortak-estimatörle bundan bağımsız. Parlak-kanal ölçümü:
kusur duyarlılığı düşük (98-S2) → 101f/105b sertifikalı setler yeterli.
"""

import numpy as np
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

def v_olc(z, qeff, om):
    g = np.diff(z)
    m = 0.5 * (z[:-1] + z[1:])
    Ldz = np.log(qeff * m / TWO_PI)
    ds = g * Ldz / TWO_PI - 1
    tt = (m - m.mean()) / (m[-1] - m[0])
    P = np.vstack([np.ones_like(tt), tt, tt**2, tt**3]).T
    ds = ds - P @ np.linalg.lstsq(P, ds, rcond=None)[0]
    X = np.vstack([np.ones_like(m), np.cos(om * m), np.sin(om * m)]).T
    b, *_ = np.linalg.lstsq(X, ds, rcond=None)
    return np.hypot(b[1], b[2]), float(Ldz.mean())

SETS = [
    ("zeta",  1, None, [], "41"),
    ("chi3",  3, "101f_chi3_zeros.npz",  [3], None),
    ("beta",  4, "101f_beta_zeros.npz",  [2], None),
    ("chi5",  5, "101f_chi5_zeros.npz",  [5], None),
    ("chi7",  7, "101f_chi7_zeros.npz",  [7], None),
    ("chi5e", 5, "105b_chi5e_zeros.npz", [5], None),
    ("chi8e", 8, "105b_chi8e_zeros.npz", [2], None),
    ("chi8o", 8, "105b_chi8o_zeros.npz", [2], None),
]
OLU2 = {"beta", "chi8e", "chi8o"}
PRIMES = list(primerange(2, 160))

SONUC = {}
for ad, qeff, dosya, olu, src in SETS:
    if src == "41":
        d41 = np.load(HERE / "41_bigT_windows.npz")
        zlist = []
        for k in ["120k", "200k"]:
            gz, tm = d41[f"gaps_{k}"], d41[f"tmid_{k}"]
            zz = np.empty(len(gz) + 1)
            zz[0] = tm[0] - gz[0] / 2
            zz[1:] = zz[0] + np.cumsum(gz)
            zlist.append(zz)
    else:
        zc = np.load(HERE / dosya)["zeros"]
        lo = np.exp(np.log(zc[0] + 1) + 0.25 * (np.log(zc[-1]) - np.log(zc[0] + 1)))
        zlist = [zc[zc >= lo]]
    rows = []
    for p in PRIMES:
        if p in olu:
            continue
        om = np.log(p)
        vs, Ls = [], []
        for zz in zlist:
            v, L = v_olc(zz, float(qeff), om)
            vs.append(v); Ls.append(L)
        v, L = float(np.mean(vs)), float(np.mean(Ls))
        tau = om / L
        if not (0.06 < tau < 0.52):
            continue
        vr = 2 * np.sin(np.pi * tau) / (np.pi * np.sqrt(p))
        rows.append((p, tau, v / vr))
    SONUC[ad] = rows
    print(f"[{ad}] {len(rows)} çizgi, τ {rows[0][1]:.3f}-{rows[-1][1]:.3f}",
          flush=True)

print(f"\n{'p':>4}", end="")
for ad, *_ in SETS:
    print(f" {ad:>7}", end="")
print("   (D = v_ölç/v_katı; τ ada başına hafif farklı)")
for i, p in enumerate([pp for pp, t, d in SONUC["beta"]]):
    print(f"{p:>4}", end="")
    for ad, *_ in SETS:
        dv = [d for pp, t, d in SONUC[ad] if pp == p]
        print(f" {dv[0]:>7.3f}" if dv else f" {'—':>7}", end="")
    print()

print("\nτ-binli kervan ortalamaları (P2 testi):")
print(f"{'τ-bin':>12} {'canlı-2 ⟨D⟩':>12} {'ölü-2 ⟨D⟩':>11} {'ln-oran':>8}")
for lo, hi in [(0.06, 0.15), (0.15, 0.25), (0.25, 0.35), (0.35, 0.45),
               (0.45, 0.52)]:
    cv, ov = [], []
    for ad, *_ in SETS:
        for pp, t, d in SONUC[ad]:
            if lo <= t < hi:
                (ov if ad in OLU2 else cv).append(d)
    if cv and ov:
        mc, mo = np.mean(cv), np.mean(ov)
        lr = np.log(mo) / np.log(mc) if (0 < mc < 1 and 0 < mo) else np.nan
        print(f"[{lo:.2f},{hi:.2f}) {mc:>12.3f} {mo:>11.3f} {lr:>8.3f} "
              f"(n={len(cv)}/{len(ov)})")
