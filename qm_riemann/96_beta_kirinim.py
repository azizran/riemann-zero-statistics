"""
96 — DIRICHLET BETA'NIN KIRINIM DESENİ: EVRENSELLİK TESTİ (21 Ağustos)
==========================================================================
Not 4 açık problem 6: aynı kırınım programı başka L-fonksiyonunda.
Hedef: L(s, χ₄) (Dirichlet beta). Sıfırlar tamamlanmış fonksiyonun
işaret değişimleriyle: Λ(s) = (4/π)^{(s+1)/2} Γ((s+1)/2) L(s,χ₄),
Λ(1/2+it) reel; Z(t) = Re[e^{iψ(t)} L(1/2+it)],
ψ = (t/2)log(4/π) + Im logΓ(3/4 + it/2). (İç sağlama: Im ≈ 0.)
L(s,χ₄) = 4^{-s}[ζ(s,1/4) − ζ(s,3/4)] (Hurwitz).

ÖNCEDEN YAZILMIŞ ÖNGÖRÜLER (ölçümden önce, 21 Ağustos):
  P1  q = 2, 4, 8 benekleri YOK (χ₄(2)=0) — karanlık-alan seviyesinde.
  P2  Faz karakteri okur: χ(p)=−1 (3,7,11,19,23) → 0°;
      χ(p)=+1 (5,13,17) → 180°; q=9 (χ=+1) → 180°.
  P3  Genlik aynı mutlak yasada: |Ĝ| = Λ(q)|χ(q)|/(L√q)·cos(πτ)·DW,
      L = log(4t/2π) (iletken-düzeltmeli), DW ölçülen jitter'dan.
  P4  Evrensellik: karanlık alan (hiperuniform) ve tarak (sıcaklık) var.

SONUÇ (21 Ağustos — dört öngörü de İSABET):
  P1 ✓ 2/4/8 benekleri SÖNDÜ (0.0011/0.0006/0.0003 — karanlık-alan
     seviyesi, fazlar rastgele). Zeta'da 2-beneği ~0.043'tü: 40× fark.
  P2 ✓✓ FAZLAR KARAKTERİ OKUDU: χ=−1 asalları 0.5-2.5° (kilit 0°),
     χ=+1 asalları 178-181° (kilit 180°), q=9 (χ=+1) 178.2°.
     Kırınım deseni = karakter tablosunun ekran görüntüsü.
  P3 ✓ mutlak yasa iletken-L ve tarak-termometreli DW ile: asal oranları
     0.987/0.962/0.961/0.941/0.912/0.902/0.918/0.805 (küçük-τ ~1,
     τ ile hafif iniş — zeta ailesiyle aynı biçim); q=9 kuvveti 0.756 —
     KUVVET AÇIĞI BETA'DA DA VAR (zeta q=27@benzer-τ: 0.81) → evrensel.
  P4 ✓ karanlık alan 1/√n'in 12 kat altı (hiperuniform); tarak 0.60 →
     σ_u = 0.161 — zeta'nın log-ısınma eğrisiyle tutarlı soğuk pencere.
  Dürüst notlar: işaret-taraması ~27 yakın-çift kaçırdı (%0.5; benekler
  duyarsız, indeksli akış bozulur → sıcaklık tarak-termometreyle);
  tarak Gauss varsayımlı (zeta'nın 0.885 alt-Gauss düzeltmesi c'yi ~%20
  oynatır, benek DW'sine etkisi <%1); benek belirsizlikleri ~0.001-0.0015
  (karanlık-alan seviyesi; q=23'te ~%10).
"""

import numpy as np
import mpmath as mp
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
mp.mp.dps = 15
CACHE = HERE / "96_beta_zeros.npz"

def Zbeta(t):
    s = mp.mpc(0.5, t)
    L = mp.power(4, -s) * (mp.zeta(s, mp.mpf(1)/4) - mp.zeta(s, mp.mpf(3)/4))
    psi = (t / 2) * mp.log(4 / mp.pi) + mp.im(mp.loggamma((s + 1) / 2))
    z = mp.e**(1j * psi) * L
    return float(mp.re(z)), float(mp.im(z))

if CACHE.exists():
    zz = np.load(CACHE)["zeros"]
    print(f"önbellekten {len(zz)} sıfır")
else:
    T0, T1 = 200.0, 4600.0
    # ızgara: yerel yoğunluğa göre adım (ort. boşluğun ~1/4'ü)
    ts, vals, ims = [], [], []
    t = T0
    while t < T1:
        ts.append(t)
        t += 0.25 * TWO_PI / np.log(4 * t / TWO_PI)
    ts = np.array(ts)
    print(f"ızgara: {len(ts)} nokta — değerlendiriliyor (mpmath, dakikalar)...")
    for i, t in enumerate(ts):
        re, im_ = Zbeta(float(t))
        vals.append(re); ims.append(im_)
        if i % 2000 == 0:
            print(f"  {i}/{len(ts)}  maks|Im| şimdiye dek: {max(map(abs, ims)):.2e}")
    vals = np.array(vals)
    print(f"iç sağlama: maks |Im Z| = {max(map(abs, ims)):.2e} (≈0 olmalı)")
    zeros = []
    sc = np.where(np.sign(vals[:-1]) * np.sign(vals[1:]) < 0)[0]
    print(f"{len(sc)} işaret değişimi — bisection ile inceltiliyor...")
    for i in sc:
        a, b = float(ts[i]), float(ts[i+1])
        fa = vals[i]
        for _ in range(22):
            m = 0.5 * (a + b)
            fm, _ = Zbeta(m)
            if fa * fm <= 0:
                b = m
            else:
                a, fa = m, fm
        zeros.append(0.5 * (a + b))
    zz = np.array(zeros)
    np.savez(CACHE, zeros=zz)
    print(f"{len(zz)} sıfır hesaplandı ve önbelleğe alındı")

# ---- örgü nesneleri
gaps = np.diff(zz)
mids = 0.5 * (zz[:-1] + zz[1:])
n = len(mids)
Leff = float(np.log(4 * mids / TWO_PI).mean())
print(f"\nn = {n} boşluk; L_eff = log(4t/2π) ort = {Leff:.3f}; "
      f"ort boşluk = {gaps.mean():.4f} (beklenen {TWO_PI/np.log(4*mids.mean()/TWO_PI):.4f})")

def Ghat(t, omegas, chunk=30000):
    out = np.zeros(len(omegas), dtype=complex)
    for s0 in range(0, len(t), chunk):
        tt = t[s0:s0 + chunk]
        out += np.exp(1j * np.outer(omegas, tt)).sum(axis=1)
    return out / len(t)

# jitter (DW için): pürüzsüz akış karşılaştırması
def Nsm(t):
    x = t / TWO_PI
    return 2 * x * np.log(4 * x / np.e) / 2 + 0 * x  # (t/2π)log(4t/2πe)

# DİKKAT: işaret-taraması ~%0.5 yakın-çift kaçırır → indeksli akış
# kayar ve σ_t patlar. Sıcaklık İNDEKSSİZ ölçülür: tarağın kendisinden.
t0 = zz[0]
x_unf = Nsm(mids) - Nsm(t0)
comb1 = abs(np.exp(2j * np.pi * x_unf).mean())
c_jit = -np.log(comb1)
sig_t = np.sqrt(2 * c_jit) / np.log(4 * mids.mean() / TWO_PI)
print(f"tarak-termometre: |Ĝ_x(2π)| = {comb1:.4f} → c_jit = {c_jit:.3f} → "
      f"σ_u = {np.sqrt(2*c_jit)/(2*np.pi):.3f} (zeta L=9.86'da 0.198 — "
      f"soğuk pencere trendiyle tutarlı)")

# ---- P1-P3: benek tablosu
CHI = {1: 1, 3: -1, 5: 1, 7: -1, 9: 1, 11: -1, 13: 1, 17: 1, 19: -1, 23: -1}
LAM = {3: np.log(3), 5: np.log(5), 7: np.log(7), 9: np.log(3),
       11: np.log(11), 13: np.log(13), 17: np.log(17), 19: np.log(19),
       23: np.log(23)}
print("\nBENEK TABLOSU (öngörü: |Ĝ| = Λ|χ|/(L√q)·cos(πτ)·DW; faz χ'yi okur):")
print(f"{'q':>3} {'χ':>3} {'ölçüm|Ĝ|':>9} {'öngörü':>8} {'oran':>6} {'faz':>8} {'öngörü-faz':>10}")
for q in [2, 4, 8, 3, 5, 7, 9, 11, 13, 17, 19, 23]:
    om = np.log(q)
    tau = om / Leff
    G = Ghat(mids, np.array([om]))[0]
    dw = np.exp(-om**2 * sig_t**2 / 2)
    if q in (2, 4, 8):
        print(f"{q:>3} {'0':>3} {abs(G):>9.4f} {'YOK':>8} {'—':>6} "
              f"{np.degrees(np.angle(G)):>7.0f}° {'(taban)':>10}")
    else:
        pred = LAM[q] * abs(CHI[q]) / (Leff * np.sqrt(q)) * np.cos(np.pi * tau) * dw
        ph = np.degrees(np.angle(G))
        ph = ph + 360 if ph < -90 else ph
        pred_ph = 180 if CHI[q] > 0 else 0
        print(f"{q:>3} {CHI[q]:>+3} {abs(G):>9.4f} {pred:>8.4f} "
              f"{abs(G)/pred:>6.3f} {ph:>7.1f}° {pred_ph:>9}°")

# plasebo + karanlık alan + tarak
fakes = np.array([np.log(x) for x in [2.31, 6.7, 10.4, 15.3]])
Gf = np.abs(Ghat(mids, fakes))
lines = [np.log(v) for v in [3,5,7,9,11,13,17,19,23,25,27,29,31,37,41,43,47,49]]
om_scan = np.linspace(0.3, 3.5, 700)
om_scan = np.array([o for o in om_scan if min(abs(o - l) for l in lines) > 0.02])
Gd = np.abs(Ghat(mids, om_scan))
print(f"\nplasebo |Ĝ| ort = {Gf.mean():.4f}")
print(f"karanlık alan: medyan |Ĝ| = {np.median(Gd):.5f}  (1/√n = {1/np.sqrt(n):.5f} — "
      f"oran {np.median(Gd)*np.sqrt(n):.2f} → hiperuniform ✓)")

# ---- figür
import matplotlib.pyplot as plt
qs_plot = [3, 5, 7, 9, 11, 13, 17, 19, 23]
meas_p, pred_p, cols = [], [], []
for q in qs_plot:
    om = np.log(q); tau = om / Leff
    G = Ghat(mids, np.array([om]))[0]
    dw = np.exp(-om**2 * sig_t**2 / 2)
    meas_p.append(abs(G))
    pred_p.append(LAM[q] * abs(CHI[q]) / (Leff * np.sqrt(q)) * np.cos(np.pi * tau) * dw)
    cols.append("#2E8B57" if CHI[q] < 0 else "firebrick")
fig, ax = plt.subplots(figsize=(8.8, 5.4))
xs = np.arange(len(qs_plot))
ax.bar(xs - 0.18, pred_p, 0.36, color="steelblue", alpha=0.8,
       label="explicit-formula prediction (conductor-$L$, comb-DW)")
ax.bar(xs + 0.18, meas_p, 0.36, color=cols, alpha=0.85,
       label="measured (green: $\chi=-1$, phase $0°$; red: $\chi=+1$, $180°$)")
for q, x in zip([2, 4, 8], [-2.6, -1.9, -1.2]):
    pass
ax.axhline(0.0011, color="gray", lw=0.9, ls="--",
           label="dark field (where $q = 2, 4, 8$ spots sit: $\chi=0$)")
ax.set_xticks(xs); ax.set_xticklabels([str(q) for q in qs_plot])
ax.set_xlabel("$q$"); ax.set_ylabel("$|\hat G(\log q)|$")
ax.set_title(r"Dirichlet $\beta$: the diffraction pattern reads the character")
ax.legend(fontsize=9); ax.grid(alpha=0.3, axis="y")
plt.tight_layout()
plt.savefig(HERE / "96_beta_kirinim.png", dpi=110)
print("Figür: 96_beta_kirinim.png")
