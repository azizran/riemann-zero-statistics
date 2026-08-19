"""
92 — D(τ) PERDESİ → GUE KÖPRÜSÜ (21 Ağustos 2026, gece)
==========================================================================
Fizik: sayım özdeşliği x_n = n − 1 − S(z_n) sıfır KONUMLARINI asal
çizgilerini tam taşımaya zorlar (benekler ✓ %1). Gap'ler FARK ölçer;
çekirdek alan asal geriliminin altında yeniden düzenleniyorsa fark
kanalı bastırılır (öteleme serbest, sıkıştırmaya direnç) → D(τ) < 1.
NULL: çekirdek çizgilerden bağımsız (additif) olsaydı D ≡ 1 olurdu —
ölçülen D < 1, çekirdek-çizgi ETKİLEŞİMİNİN kanıtı.

KÖPRÜ HİPOTEZİ: bastırma oranı, gazın İÇ korelasyonlarının evrensel bir
fonksiyonu — GUE/CUE log-gazının doğrusal yanıtından hesaplanabilir:
  D(κ) = [gap-yanıtı] / [katı-dalga gap-yanıtı]  (κ = 2πτ)
CUE'da MC'siz, kovaryansla kesin (zayıf tek-cisim potansiyele yanıt):
  D_CUE(κ) = |⟨G_m ρ_m*⟩| / [ (2 sin(κ/2)/κ) · ⟨|ρ_m|²⟩ ],  κ = 2πm/N
  (ρ_m = Σe^{imθ}; G_m = Σ(δs)e^{imθ̄}; sağlama: ⟨|ρ_m|²⟩ = min(m,N),
   Diaconis–Shahshahani.)

  T1  ölçülen D(τ): v(83, tam-taban) / [2τ · a_v(τ)] — a_v boru-hattı
      çarpanı 6 noktalı enjeksiyon kalibrasyonuyla (89-T3 genişletmesi).
  T2  null (analitik): additif çekirdek ⟹ D ≡ 1 — ölçüm bunu dışlıyor.
  T3  D_CUE(κ) eğrisi (N=256, M=3000 CUE örneği) vs ölçülen D(2πτ).
  T4  zeta SPONTANE koşullu yanıtı (çizgi-dışı bantlar, pseudo-topluluk).

SONUÇ (21 Ağustos gecesi):
  T1: ölçülen perde temiz: D = 0.908/0.863/0.797/0.726/0.634/0.484
      (τ-binleri 0.05→0.55); a_v = 0.999→0.848.
  T3: NAİF KÖPRÜ ÖLDÜ — CUE koşullu yanıtı PERDELEMEZ, TERS yönde:
      D_CUE = 1.00→1.18 (κ→π; sağlama ⟨|ρ_m|²⟩=min(m,N) ✓).
  T4: ASIL KEŞİF — zeta kendi spontane modlarına CUE GİBİ yanıt veriyor:
      D_spont = 1.075/1.088/1.177 @ τ=0.3/0.4/0.5 vs CUE 1.081/1.132/1.183
      (Nyquist'te binde-beş uyum!). τ=0.2: 0.956 (hafif altta); τ=0.1:
      0.604 — ŞÜPHELİ (yoğunluk-sürüklenmesi düşük-ω kirliliği adayı,
      detrend edilmemiş; bayraklı). SENTEZ: perdeleme gazın değil
      SÜRÜCÜNÜN özelliği — koheran aritmetik dalgalar perdeleniyor,
      spontane gürültü CUE-gibi geçiyor. D_asal/D_spont ≈ 0.63/1.09 @
      τ≈0.45: saf aritmetik-koherans fiziği izole edildi. Açık: neden
      koheran sürücü perdelenir (adyabatik-vs-ani; sayım-özdeşliği
      gömülmesi); düşük-τ spontane detrend tekrarı.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

# ============ T1a: a_v(τ) enjeksiyon kalibrasyonu ============
d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
gaps, amps, tmid = d41[f"gaps_{K41[1]}"], d41[f"amps_{K41[1]}"], d41[f"tmid_{K41[1]}"]
L0 = float(np.log(tmid / TWO_PI).mean())
z = np.empty(len(gaps) + 1)
z[0] = tmid[0] - gaps[0] / 2
z[1:] = z[0] + np.cumsum(gaps)

print("T1a — a_v(τ) enjeksiyon kalibrasyonu (L=10.37):")
tau_cal, av_cal = [], []
for taustar in [0.06, 0.12, 0.20, 0.30, 0.40, 0.50]:
    om = taustar * L0
    U = 0.02
    zp = z + U * np.cos(om * z)
    gp = np.diff(zp)
    mp_ = 0.5 * (zp[:-1] + zp[1:])
    yg = np.log(gp * np.log(mp_ / TWO_PI) / TWO_PI)
    X = np.vstack([np.ones_like(mp_), np.cos(om * mp_), np.sin(om * mp_)]).T
    b, *_ = np.linalg.lstsq(X, yg, rcond=None)
    yg0 = np.log(gaps * np.log(tmid / TWO_PI) / TWO_PI)
    X0 = np.vstack([np.ones_like(tmid), np.cos(om * tmid), np.sin(om * tmid)]).T
    b0, *_ = np.linalg.lstsq(X0, yg0, rcond=None)
    a_v = np.hypot(b[1] - b0[1], b[2] - b0[2]) / (U * om)
    tau_cal.append(taustar); av_cal.append(a_v)
    print(f"  τ* = {taustar:.2f}: a_v = {a_v:.3f}")
tau_cal, av_cal = np.array(tau_cal), np.array(av_cal)

# ============ T1b: ölçülen D(τ) ============
d83 = np.load(HERE / "83_tam_taban_egri.npz")
tauv, vv, Lv = d83["tau"], d83["v"], d83["L"]
a_interp = np.interp(tauv, tau_cal, av_cal)
D_emp = vv / (2 * tauv * a_interp)
print(f"\nT1b — ölçülen D(τ) ({len(vv)} nokta):")
for lo, hi in [(0.0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5), (0.5, 0.62)]:
    m = (tauv >= lo) & (tauv < hi)
    if m.sum():
        print(f"  τ∈[{lo:.2f},{hi:.2f}): ⟨D⟩ = {D_emp[m].mean():.3f} ± "
              f"{D_emp[m].std()/np.sqrt(m.sum()):.3f} (n={m.sum()})")
print("T2 — null: additif çekirdek ⟹ D ≡ 1; ölçüm D<1 → etkileşim gerçek.")

# ============ T3: CUE doğrusal-yanıt eğrisi ============
print("\nT3 — D_CUE(κ) (N=256, M=3000; MC'siz kovaryans yanıtı):")
rng = np.random.default_rng(92)
N, M = 256, 3000
ms = [8, 13, 26, 51, 77, 102, 115, 128]
rho = np.zeros((M, len(ms)), dtype=complex)
Gm = np.zeros((M, len(ms)), dtype=complex)
for s in range(M):
    A = (rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(A)
    Q = Q * (np.diagonal(R) / np.abs(np.diagonal(R)))
    th = np.sort(np.angle(np.linalg.eigvals(Q)))
    dth = np.diff(np.concatenate([th, [th[0] + TWO_PI]]))
    mid = th + dth / 2
    ds = dth * N / TWO_PI - 1
    for j, m_ in enumerate(ms):
        rho[s, j] = np.exp(1j * m_ * th).sum()
        Gm[s, j] = (ds * np.exp(1j * m_ * mid)).sum()

print(f"  {'m':>4} {'κ/2π=τ':>7} {'⟨|ρ|²⟩':>7} {'min(m,N)':>8} {'D_CUE':>7}")
D_cue = []
for j, m_ in enumerate(ms):
    kap = TWO_PI * m_ / N
    var_r = np.mean(np.abs(rho[:, j])**2)
    cov = np.mean(Gm[:, j] * np.conj(rho[:, j]))
    D = abs(cov) / ((2 * np.sin(kap / 2) / kap) * var_r)
    D_cue.append((m_ / N, D))
    print(f"  {m_:>4} {m_/N:>7.3f} {var_r:>7.1f} {min(m_, N):>8} {D:>7.3f}")
D_cue = np.array(D_cue)

# ============ T4: ZETA SPONTANE KOŞULLU YANITI ============
print("\nT4 — zeta spontane (çizgi-dışı bantlar, 6 pencere birleşik):")
lines_all = [np.log(q) for q in
    [2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53,59,
     61,64,67,71,73,79,81,83,89,97,101,103,107,109,113,121,125,127,128]]
D_sp = []
for tau0 in [0.10, 0.20, 0.30, 0.40, 0.50]:
    num = 0.0 + 0j; den = 0.0
    for k in K41[:6]:
        gz, tm = d41[f"gaps_{k}"], d41[f"tmid_{k}"]
        Lz = float(np.log(tm / TWO_PI).mean())
        zz = np.empty(len(gz) + 1)
        zz[0] = tm[0] - gz[0] / 2
        zz[1:] = zz[0] + np.cumsum(gz)
        dsz = gz / gz.mean() - 1
        oms = tau0 * Lz + np.linspace(-0.05 * Lz, 0.05 * Lz, 240)
        oms = np.array([o for o in oms
                        if min(abs(o - l) for l in lines_all) > 0.01])
        for s0 in range(0, len(oms), 60):
            ob = oms[s0:s0+60]
            rr = np.exp(1j * np.outer(ob, zz)).sum(axis=1)
            GG = (np.exp(1j * np.outer(ob, tm)) * dsz[None, :]).sum(axis=1)
            num += (GG * np.conj(rr)).sum()
            den += (np.abs(rr)**2).sum()
    kap = 2 * np.pi * tau0
    D_sp.append((tau0, abs(num) / ((2 * np.sin(kap / 2) / kap) * den)))
    print(f"  τ≈{tau0:.2f}: D_spont = {D_sp[-1][1]:.3f}")
D_sp = np.array(D_sp)

# ============ FİGÜR ============
fig, ax = plt.subplots(figsize=(8.6, 5.6))
ax.plot(tauv, D_emp, "o", ms=4, c="steelblue", alpha=0.55,
        label="ölçülen D = v / (2τ·a_v)  (132 nokta)")
ax.plot(D_cue[:, 0], D_cue[:, 1], "s-", c="firebrick", ms=7, lw=1.4,
        label="CUE doğrusal yanıtı D(κ=2πτ)")
ax.plot(D_sp[:, 0], D_sp[:, 1], "^", c="#2E8B57", ms=8,
        label="zeta spontane modları (çizgi-dışı)")
ax.axhline(1, color="gray", lw=0.8, ls=":", label="null: additif çekirdek (D=1)")
ax.set_xlabel(r"$\tau = \kappa/2\pi$")
ax.set_ylabel("D — fark-kanalı perdelemesi")
ax.set_title("Perde köprüsü: sıfır gazı vs CUE log-gazı")
ax.set_ylim(0.3, 1.28)
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(HERE / "92_perde.png", dpi=110)
print("\nFigür: 92_perde.png")
