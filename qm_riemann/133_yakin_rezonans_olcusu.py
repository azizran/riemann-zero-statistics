"""
133 — "≈" ÖLÇÜSÜNÜN DARALTILMASI: KÂĞIT FORMÜLÜNÜN DOĞRUDAN DEĞERİ (27 Ağu)
==========================================================================
132c'nin bıraktığı tek açık uç kapatılıyor. Kâğıt formülü (tam açık):
  R_nn(τ_p) = [Σ_{Q'∈kuyruk} Σ_{Q: asal-kuvvet, δ=log(Q/(pQ'))}
               A_Q A_{Q'} cos(π(τ_Q+τ_{Q'})) · Φ(δ)] / (2 c₁ A_p)
  Φ(δ) = cos(δ·t_c)·sinc(δT/2)      [OLS pencere-transferi, KESİN]
  A_Q = 2a_Q sin(πτ_Q);  c₁ = ½Σ_{Q∈kuyruk} A_Q² cos(2πτ_Q)
"Aritmetik artık" = asal-kuvvetlerin pQ' çevresindeki dağılımı
(δ ≈ r/(pQ'), r = en yakın komşu uzaklıkları).

ÖN-MÜHÜRLER:
  N1  Formül, lab parametreleriyle (L̄=8.11, D3-kuyruk τ∈(0.52,1.4],
      pencere T ve t_c lab'ınki) lab ölçümünü (R_nn −1.71/−1.79/−1.68/
      −1.55 @ p=2,3,5,7) ±%25 içinde vurmalı → "≈" ölçüsü KAPANIR.
  N2  δ-kabuk ayrıştırması: katkının çekirdeği |δ| ≲ 2π/T kabuğunda
      olmalı (uzak kabuklar sinc'le sönmeli) — ölçünün anatomisi.
  N3  c₁ kapalı formu, lab'ın ölçtüğü c₁'i (c₁/σ² ≈ −0.57'den) ±%20
      içinde vermeli (payda sağlaması).

SONUÇ (27 Ağustos) — N1 KESİN RET, VE İKİ ÖĞRETİCİ İMZA:
  R_kağıt = +0.05/+0.06/+0.01/−0.02 (lab −1.71..−1.55): STATİK
  ÇİFT-VURUŞ TOPLAMI SIFIR — pürüzsüz ızgarada kuyruk çift-beatleri
  ω_p'ye net katkı vermiyor (130'un kule-retiyle tam tutarlı; 132c'nin
  "≈"-toplamı bu haliyle BOŞ kapıymış — kapandı).
  İMZA 1 (σ-uyuşmazlığı → EMİLİM KANITI): kapalı-form σ_η² = 0.0743;
  koheran lab 0.0266; KARIŞIK-FAZLI lab (130-F2) 0.243² = 0.0590 ≈
  kapalı forma yakın → öz-tutarlı örgü, cos-kilitli kuyruğun ~2/3
  gücünü EMİYOR; kapalı-form genlikler ancak fazsız/karışık durumda
  geçerli. Gerçek nesne, emilmiş/renormalize kuyruk.
  İMZA 2 (mekanizmanın adresi kesinleşti): ω_p-içeriği statik beat'te
  değil, KONUM-FAZ MODÜLASYONUNDA — her kuyruk çizgisi p-dalgasınca
  yerinden oynatılmış noktalarda örnekleniyor: cos(ω_Q(m̄+u_p)) →
  Σ_Q A_Q²·(mod-indeksi ω_Q a_p/ρ̄)·(geometri) TEK-TOPLAMI (132b'nin
  C4 kanalı, %66-92) EKSİ faz-duyarlı KARŞIT-KANAL (iptal; seçicilik
  = iptalin aritmetik bozulması, F2/FZ-C3 ölçüleriyle tutarlı).
  KALEMİN SIRADAKİ NESNESİ: tek-toplam + karşıt-kanal muhasebesi
  (çift enumerasyonu yok; Σ_Q A_Q²... kapalı form yazılabilir) —
  "−2 = kinematik tek-toplamın iptalden artakalanı."
"""

import numpy as np
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
T_win = float(tg[-1] - tg[0])
t_c = float(0.5 * (tg[0] + tg[-1]))
Lb = float(np.log((0.5 * (tg[:-1] + tg[1:])) / TWO_PI).mean())
gbar = TWO_PI / Lb
print(f"lab parametreleri: L̄={Lb:.3f}  T={T_win:.0f}  t_c={t_c:.0f}  "
      f"çözünürlük 2π/T={TWO_PI/T_win:.2e}")

QCAP = int(np.exp(1.4 * L0))
TAIL = [q for q in pk_list(QCAP) if np.log(q) / Lb > 0.52]
TAILset = TAIL
logT = np.array([np.log(q) for q in TAIL])
tauT = logT / Lb

def a_of(q):
    (pp, kk), = factorint(q).items()
    return np.log(pp) / (np.pi * np.sqrt(q) * np.log(q))

aT = np.array([a_of(q) for q in TAIL])
AT = 2 * aT * np.sin(np.pi * tauT)

# N3: c1 kapalı formu
c1_kapali = 0.5 * float((AT**2 * np.cos(2 * np.pi * tauT)).sum())
s2_kapali = 0.5 * float((AT**2).sum())
print(f"N3 — kapalı form: σ_η² = {s2_kapali:.4f} (lab 0.163²={0.163**2:.4f})"
      f"  c₁ = {c1_kapali:+.5f}  c₁/σ² = {c1_kapali/s2_kapali:+.3f} "
      f"(lab −0.57)")

# tüm asal-kuvvetler (Q adayları) 7·QCAP'e kadar
ALLQ = pk_list(7 * QCAP + 100)
logA = np.array([np.log(q) for q in ALLQ])

def sinc(x):
    return np.sinc(x / np.pi)

print(f"\nN1/N2 — kâğıt formülünün değeri (lab ölçümü: "
      f"−1.71/−1.79/−1.68/−1.55):")
print(f"{'p':>3} {'τ_p':>6} {'R_kagit':>8} {'lab':>7} {'−2cosπτ':>8}  "
      f"kabuklar[|δ|T/2π: <0.5 / 0.5-2 / 2-8 / >8]")
LAB = {2: -1.713, 3: -1.789, 5: -1.678, 7: -1.551}
for p in [2, 3, 5, 7]:
    om_p = np.log(p)
    tau_p = om_p / Lb
    A_p = 2 * a_of(p) * np.sin(np.pi * tau_p)
    kab = np.zeros(4)
    top = 0.0
    for iq, Qp in enumerate(TAIL):
        hedef = np.log(Qp) + om_p            # log(p·Q')
        # yakın asal-kuvvetler: |δ| < 40·(2π/T) penceresi
        dmax = 40 * TWO_PI / T_win
        i0 = np.searchsorted(logA, hedef - dmax)
        i1 = np.searchsorted(logA, hedef + dmax)
        for j in range(i0, i1):
            Q = ALLQ[j]
            if Q == p * Qp:                  # kule (tam rezonans) DAHİL
                pass
            dlt = logA[j] - hedef
            tQ = logA[j] / Lb
            if tQ * Lb / Lb <= 0.52:         # Q kuyrukta olmalı (fit-dışı)
                continue
            AQ = 2 * a_of(Q) * np.sin(np.pi * tQ)
            w = np.cos(dlt * t_c) * sinc(dlt * T_win / 2)
            katki = AQ * AT[iq] * np.cos(np.pi * (tQ + tauT[iq])) * w
            top += katki
            s = abs(dlt) * T_win / TWO_PI
            kab[0 if s < 0.5 else 1 if s < 2 else 2 if s < 8 else 3] += katki
    R_kagit = top / (2 * c1_kapali * A_p)
    print(f"{p:>3} {tau_p:>6.3f} {R_kagit:>+8.3f} {LAB[p]:>+7.3f} "
          f"{-2*np.cos(np.pi*tau_p):>+8.3f}  "
          f"[{kab[0]/(2*c1_kapali*A_p):+.2f} {kab[1]/(2*c1_kapali*A_p):+.2f} "
          f"{kab[2]/(2*c1_kapali*A_p):+.2f} {kab[3]/(2*c1_kapali*A_p):+.2f}]",
          flush=True)
