"""
116 — TERMAL AYNA: DENGE-CUE'NUN KENDİ DALGALARINA NEFESİ (26 Ağustos)
==========================================================================
Köprü kapandı (115): perde = nefes(β≈−2) × çarpıklık(λ). SON SORU:
bu davranış EVRENSEL log-gaz dengesi mi, aritmetiğe mi özgü?
Ölçüm: CUE'da dış dalga YOK — gazın KENDİ spontane yoğunluk dalgaları
sürücü rolünde (gerçek ölçümün birebir termal aynası): örneklem-başına
gap-dalgası G_m yön belirler; varyans-tarağı V_m ve küp-tarağı T_m
o yöne izdüşürülür, toplulukta ortalanır (92-T3 kovaryans usulü,
Metropolis'siz).

Tanımlar (örneklem s, mod m; θ sıralı, dθ sarmal, ds = dθN/2π−1,
θ̄ = orta): G = Σds·e^{imθ̄}; V = Σ(ds²−⟨ds²⟩_s)e^{imθ̄};
T = Σ(ds³−⟨ds³⟩_s)e^{imθ̄}; ρ = Σe^{imθ}.
  D_eq  = |⟨G ρ*⟩| / (c·⟨|ρ|²⟩), c = 2sin(κ/2)/κ   [92-T3 çaprazı]
  R_eq  = Re⟨V G*⟩ / (2σ²⟨|G|²⟩)
  R3_eq = Re⟨T G*⟩ / (σ³·... σ³⟨|G|²⟩/σ... birim: Re⟨VG*⟩ normuyla aynı
          yapı: /(σ³⟨|G|²⟩)·σ² → raporda iki norm da basılır)
  Rnn_eq = Re⟨W G*⟩/(2c₁⟨|G|²⟩), W = Σ(ds_nds_{n+1}−c₁)e^{imθ̄}

ÖN-MÜHÜR:
  T1  D_eq bilinen 92-T3 eğrisini (1.00→1.18) üretmeli (çapraz ✓).
  T2  KARAR: R_eq derin negatifse (≈ −1..−3; gerçeğin içsel β'sı
      gibi) → TERS NEFES EVRENSEL DENGE FİZİĞİ (büyük birleşme:
      aritmetik yalnız dalgaları sağlar, taşıma log-gazın malı).
      R_eq ≈ 0/pozitifse → tanıma ARİTMETİĞE ÖZGÜ.
  T3  R3_eq ve Rnn_eq işaretleri gerçekle kıyas (tam dörtlü ayna).
N=256, M=3000; modlar m = 8, 26, 51, 77, 115 (τ ≈ 0.03-0.45).

SONUÇ (26 Ağustos) — AYNA KONUŞTU: TANIMA ARİTMETİĞE ÖZGÜ:
  T1 ✓ D_eq = 1.001→1.160 (92-T3'ün 1.00→1.18'i birebir; makine doğru).
  T2 KARAR — DENGE-CUE KENDİ DALGALARINA ADYABATİK NEFES ALIR:
     R_eq = +1.02/+0.97/+0.89/+0.78/+0.47 (boyalı-kinematik eğrinin
     üstünde/üzerinde; NEGATİF DEĞİL). Rnn_eq de pozitif (+0.98→+0.20);
     R3_eq güçlü POZİTİF (σ³-norm +3.6; m₃-norm +5.7 — termal gazın
     çarpıklığı dalgasıyla BİRLİKTE sürüklenir, ölçek-sürüklemesi gibi).
  HÜKÜM: ζ-gazının ters nefesi (üç momentte) HİÇBİR denge log-gaz
  davranışında yok — ne boyalı, ne spontane-denge; işaretler TERS.
  "Gaz kendi aritmetiğini tanır" artık termal-aynalı, üç-momentli,
  işaret-çevrik bir ifade: denge gazı kendi dalgalarını adyabatik
  taşır; ζ-gazı kendi aritmetik dalgalarına karşı gürültüsünü ters
  örgütler. Aritmetik dalgalar gazın "kendi termal dalgası gibi"
  DEĞİL — özel muamele görüyor. NOT 6 ADAYI ANA CÜMLE.
  (Yan not: CUE tam-gap skew +0.50 / c₁ −0.31 vs gerçek η +0.24 /
  −0.52 — artık-istatistik de yapıca farklı; karşılaştırmalar
  G-yönü izdüşüm konvansiyonunda, birim-notu docstring'de.)
"""

import numpy as np

TWO_PI = 2 * np.pi
rng = np.random.default_rng(116)
N, M = 256, 3000
MS = [8, 26, 51, 77, 115]

acc = {m: np.zeros(6) for m in MS}   # [ReVG*, ReTG*, ReWG*, |G|², |ρ|², ReGρ*/|Gρ|-payı]
accD = {m: [0.0 + 0j, 0.0] for m in MS}
s2s, m3s, c1s = [], [], []
for s in range(M):
    A = (rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))) / np.sqrt(2)
    Q, Rq = np.linalg.qr(A)
    Q = Q * (np.diagonal(Rq) / np.abs(np.diagonal(Rq)))
    th = np.sort(np.angle(np.linalg.eigvals(Q)))
    dth = np.diff(np.concatenate([th, [th[0] + TWO_PI]]))
    mid = th + dth / 2
    ds = dth * N / TWO_PI - 1
    s2 = float((ds**2).mean()); m3 = float((ds**3).mean())
    ee = ds * np.roll(ds, -1)
    c1 = float(ee.mean())
    s2s.append(s2); m3s.append(m3); c1s.append(c1)
    for m in MS:
        e = np.exp(1j * m * mid)
        G = (ds * e).sum()
        V = ((ds**2 - s2) * e).sum()
        T = ((ds**3 - m3) * e).sum()
        W = ((ee - c1) * e).sum()
        rho = np.exp(1j * m * th).sum()
        a = acc[m]
        a[0] += (V * np.conj(G)).real
        a[1] += (T * np.conj(G)).real
        a[2] += (W * np.conj(G)).real
        a[3] += abs(G)**2
        a[4] += abs(rho)**2
        accD[m][0] += G * np.conj(rho)
        accD[m][1] += abs(rho)**2

s2 = float(np.mean(s2s)); m3 = float(np.mean(m3s)); c1 = float(np.mean(c1s))
print(f"CUE: σ²={s2:.4f}  skew={m3/s2**1.5:+.3f}  c₁/σ²={c1/s2:+.3f}  "
      f"(gerçek η: 0.148²=0.022 ölçekli; skew +0.24; c₁/σ² −0.52)")
print(f"\n{'m':>4} {'τ':>6} {'D_eq':>6} {'R_eq':>7} {'Rnn_eq':>8} "
      f"{'R3_eq(σ³)':>10}")
for m in MS:
    tau = m / N
    kap = TWO_PI * tau
    c = 2 * np.sin(kap / 2) / kap
    a = acc[m]
    D = abs(accD[m][0]) / (c * accD[m][1])
    Rp = a[0] / (2 * s2 * a[3])
    Rnn = a[2] / (2 * c1 * a[3])
    R3 = a[1] / (s2**1.5 * a[3]) * np.sqrt(s2)   # = ReTG*/(σ·s2·|G|²)·...
    R3b = a[1] / (3 * m3 * a[3])
    print(f"{m:>4} {tau:>6.3f} {D:>6.3f} {Rp:>+7.3f} {Rnn:>+8.3f} "
          f"{R3:>+10.3f}  [m₃-norm {R3b:+.3f}]", flush=True)
print("\nGERÇEK (ζ, kendi dalgaları): R_p −0.85 düz; R_nn −1.9→−1.2; "
      "R₃(σ³) −1.0..−1.45; içsel β≈−2..−3")
