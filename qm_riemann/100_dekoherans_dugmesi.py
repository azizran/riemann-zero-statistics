"""
100 — DEKOHERANS DÜĞMESİ: f ~ τ^3.3'ÜN MEKANİK TÜRETİMİ (20 Ağustos)
==========================================================================
Bilanço: 91 tam-koheran lab → açık var AMA eğim ~τ^1.8 ve fazlar dönüyor;
gerçek → τ^3.3 ve faz-saf (≤1°); 99 → mekanizma aritmetik (aynı ω'da
asal bağışık), (1−ε)/k = f(τ) çökmesi.

HİPOTEZ: gerçek yasa = aile çapraz-terimleri × DEKOHERANS SÜZGECİ.
Rastgele alan, koheran ikinci-mertebenin döndüren (tek) kısmını öldürür;
kalan faz-koruyucu kısım daha dik bir τ-yasasına iner.

DENEY: 91 lab'ı + ayarlanabilir rastgele yerdeğiştirme η ~ N(0, σ_n):
  σ_n / σ_ölçülen ∈ {0, 0.5, 1.0, 1.5}   (σ_ölçülen: kampanya, L≈9.5)
Her ayarda: kuvvet-beneği fazları (dönme?) + (1−ε)/k vs τ (eğim fiti).
BAŞARI ÖLÇÜTÜ (ön-mühürlü): σ_n = σ_ölç'te lab, gerçeğin ÜÇLÜSÜNE
oturmalı — faz-saflık ≤~1-2°, eğim ~3.0-3.6, genlikler ×2 içinde.
Otursa → mekanizma türetildi; kapalı form = bu modelin asimptotiği.

SONUÇ (20 Ağustos — hipotez KISMEN yanlış, türetim yine de KAPANDI):
  1. Taban-iptalli oran-testi (90-T2 usulü, lab'ın kendi asal eğrisine
     bölme) uygulanınca SAF-KOHERAN lab (σ_n=0) gerçek yasayı VURDU:
     eğim 3.26 (gerçek 3.3), f(0.23)=0.0306 (gerçek 0.030) —
     SIFIR SERBEST PARAMETRE. → KUVVET-AÇIĞI YASASI MEKANİK OLARAK
     TÜRETİLDİ: örtük sıfır-koşulu (ρ̄u = −S(t+u)) + EF aile-çizgileri
     + orta-nokta örneklemesi, oran-testiyle okununca k·f(τ), f≈4τ^3.3
     üretir. 91'in "τ^1.8" hükmü mutlak-taban yanlılığı artefaktıydı —
     DÜZELTİLDİ.
  2. Dekoherans hipotezi ÖLDÜ: iid gürültü eğimi bozuyor (2.6→2.3) —
     gerçek rastgelelik iid değil (hiperuniform antikorelasyon).
  3. AÇIK KALAN: faz-saflığı — lab ~4-7° döndürüyor, gerçek ≤1°.
     Döndüren (tek) bileşeni gerçekte ~10× bastıran şey (korele gürültü?
     öz-tutarlılıkta S'S–örnekleme iptali?) hâlâ türetilmedi.
  4. Kapalı form adayları (validen modelin asimptotiği): saf τ^3
     yetersiz, τ^3.3 etkin — cos/log düzeltmeli τ³ ailesi; analitik
     açılım sıradaki kâğıt-kalem ödevi.
"""

import numpy as np
from pathlib import Path
from sympy import primerange, factorint

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
L0 = 9.5
NZ = 60000
SIG_MEAS_U = 0.23                      # kampanya, L≈9.5 (unfold)
QPOW = [(4, 2, 2), (8, 2, 3), (9, 3, 2), (16, 2, 4), (25, 5, 2),
        (27, 3, 3), (32, 2, 5), (49, 7, 2)]
PRC = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]   # asal eğrisi (taban-iptali)

def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
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

# pürüzsüz örgü + örtük çizgi-alanı (91 makinesi)
t0 = TWO_PI * np.exp(L0)
idx = np.arange(NZ + 1, dtype=float)
t = t0 + idx * TWO_PI / L0
for _ in range(8):
    t = t - (rvm_N(t) - rvm_N(t0) - idx) / (np.log(t / TWO_PI) / TWO_PI)
rho = np.log(t / TWO_PI) / TWO_PI
u = np.zeros_like(t)
for _ in range(60):
    u = 0.5 * u + 0.5 * (-S_field(t + u) / rho)
gap_mean = TWO_PI / L0
rng = np.random.default_rng(100)

print(f"lab: L={L0}, n={NZ}; σ_ölç(unfold) = {SIG_MEAS_U} → t-birim "
      f"{SIG_MEAS_U*gap_mean:.4f}")
print(f"{'σ_n/σ_ölç':>9} {'⟨|Δfaz|⟩kuvvet':>14} {'eğim α':>7} {'f(0.23)':>8} "
      f"{'asal-oran':>9}")

SEEDS = [100, 101, 102, 103]
for frac in [0.0, 0.5, 1.0, 1.5]:
    accF = {q: [] for q, _, _ in QPOW}
    accPh = {q: [] for q, _, _ in QPOW}
    for sd in ([100] if frac == 0.0 else SEEDS):
        r2 = np.random.default_rng(sd)
        eta = r2.normal(0, frac * SIG_MEAS_U * gap_mean, NZ + 1)
        z = t + u + eta
        mids = 0.5 * (z[:-1] + z[1:])
        L = float(np.log(mids / TWO_PI).mean())
        sig_t = float((mids - 0.5 * (t[:-1] + t[1:])).std())
        # lab'ın KENDİ asal eğrisi (taban-iptali — 90-T2 usulü)
        prT, prR = [], []
        for p in PRC:
            om = np.log(p); tau = om / L
            dw = np.exp(-om**2 * sig_t**2 / 2)
            pred = np.log(p) / (L * np.sqrt(p)) * np.cos(np.pi * tau) * dw
            prT.append(tau)
            prR.append(abs(Ghat(mids, np.array([om]))[0]) / pred)
        prT, prR = np.array(prT), np.array(prR)
        for q, pp, kk in QPOW:
            om = np.log(q); tau = om / L
            dw = np.exp(-om**2 * sig_t**2 / 2)
            pred = np.log(pp) / (L * np.sqrt(q)) * np.cos(np.pi * tau) * dw
            G = Ghat(mids, np.array([om]))[0]
            base = np.interp(tau, prT, prR)      # asal-eğri tabanı
            eps = (abs(G) / pred) / base
            ph = (np.degrees(np.angle(G)) + 360) % 360
            accF[q].append(1 - eps)
            accPh[q].append(min(abs(ph-180), abs(ph-0), abs(ph-360)))
    rows = []
    for q, pp, kk in QPOW:
        d = np.mean(accF[q]); dph = np.mean(accPh[q])
        tau = np.log(q) / L
        if d > 1e-4:
            rows.append((tau, d / kk, dph, q))
    A = np.array([(r[0], r[1], r[2]) for r in rows])
    c = np.polyfit(np.log(A[:, 0]), np.log(A[:, 1]), 1)
    alpha, f023 = c[0], np.exp(c[1] + c[0] * np.log(0.23))
    print(f"{frac:>9.1f} {A[:,2].mean():>13.1f}° {alpha:>7.2f} {f023:>8.4f}")
    if frac == 1.0:
        print("   σ_n = σ_ölç ayrıntı (taban-iptalli):")
        for tau, fv, dph, q in rows:
            print(f"     q={q:>2}: τ={tau:.3f}  f={fv:.4f}  Δfaz={dph:.1f}°")
print("\nGERÇEK HEDEF: eğim ≈ 3.3, f(0.23) ≈ 0.030, Δfaz ≤ ~1-2°")
