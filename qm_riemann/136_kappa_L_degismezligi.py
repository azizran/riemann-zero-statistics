"""
136 — κ'NIN L-DEĞİŞMEZLİĞİ: EVRENSEL −2'NİN KESİM-KARARLILIĞI (28 Ağu)
==========================================================================
135'in kapalı formu κ_çekirdek = ΣA²πτcos3πτ/sinπτ / ΣA²cos2πτ
[kuyruk] lab kesimiyle −3.18 verdi. SORU: gerçek gazda kuyruğu kesen
şey DW (jitter) — κ, DW-ağırlıklı kuyrukta L'den bağımsız mı?
(120/127/129: R_nn = −2cosπτ yedi adada ve L 9.9-12.5'te evrensel;
κ_ad ≈ +1.2 ile κ_çekirdek ≈ −3.2 sabit kalmalı.)

DW ağırlığı: w(Q) = exp(−ω_Q²σ_t²), σ_t = σ_u(L)·2π/L; σ_u(L)
ısınma eğrisinden (101-P4: 0.2287@6.99 → 0.2833@12.45; log-interp).
ÖN-MÜHÜR: κ(L) L=8..13'te −3.2 civarında ~±%10 bandında kalmalı
(evrenselliğin kesim-kararlılık açıklaması); L=24.5/44.6'da da
basılır (öngörü keskin değil — kuyruk bileşimi orada bambaşka).

SONUÇ (28 Ağustos) — L-DEĞİŞMEZLİK ✓, DEĞER SÜRPRİZİ:
  κ_DW(L) = −2.16/−2.26/−2.22/−2.15/−2.14/−2.13/−2.12 (L=8→12.45):
  DÜZ ±%3 — EVRENSELLİĞİN AÇIKLAMASI BU (DW-ağırlıklı çekirdek toplamı
  fiziksel L-aralığında sabit). DEĞER: DW-kesimli κ ≈ −2.15 — lab'ın
  keskin-kesimli −3.18'inden farklı, ve TEK BAŞINA ölçülen −2.00±0.09'a
  %6-10 yakın. YENİDEN-OKUMA: gerçek gazda ayrışım κ ≈ κ_DW(−2.15) +
  küçük adyabatik artık (+~0.15); lab'ın (−3.18 + 1.2) bölünmesi
  keskin-kesimin eseri. "−2" sihirli tamsayı değil, DW-ağırlıklı
  toplamın fiziksel rejimdeki değeri — TAM −2 mi (derin özdeşlik)
  yoksa ≈−2.1 mi: AÇIK (adyabatik artığın hassas ölçümüyle ayrışır).
  Not: L=24.5/44.6 NaN — asal listesi (2e5) o derinliğe yetmiyor,
  kuyruk boş kalıyor (bilinçli sınır; derin-L ayrı iş).
"""

import numpy as np
from sympy import primerange, factorint

def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

LS = [6.99, 7.69, 8.39, 9.08, 9.86, 10.37, 10.93, 11.47, 11.98, 12.45]
SU = [0.2287, 0.2384, 0.2471, 0.2549, 0.2626, 0.2673, 0.2721, 0.2766,
      0.2803, 0.2833]

def sigma_u(L):
    return float(np.interp(L, LS, SU))

QALL = pk_list(200000)
LAM = {}
for q in QALL:
    (pp, kk), = factorint(q).items()
    LAM[q] = float(np.log(pp))

print(f"{'L':>6} {'σ_u':>6} {'κ_çekirdek':>11} {'κ_ad+κ':>8}  (hedef ≈ −3.2 / −2.0)")
for L in [8.0, 9.0, 9.86, 10.37, 11.0, 12.0, 12.45, 24.5, 44.6]:
    su = sigma_u(min(L, 12.45)) if L <= 13 else 0.30
    sig_t = su * 2 * np.pi / L
    num = den = 0.0
    for q in QALL:
        om = np.log(q)
        tau = om / L
        if tau <= 0.52:
            continue
        w = np.exp(-om**2 * sig_t**2)
        if w < 1e-8:
            break
        a = LAM[q] / (np.pi * np.sqrt(q) * om)
        A2 = (2 * a * np.sin(np.pi * tau))**2 * w
        num += A2 * np.pi * tau * np.cos(3 * np.pi * tau) / np.sin(np.pi * tau)
        den += A2 * np.cos(2 * np.pi * tau)
    kap = num / den if den != 0 else float("nan")
    print(f"{L:>6.2f} {su:>6.3f} {kap:>+11.3f} {1.2 + kap:>+8.2f}", flush=True)
