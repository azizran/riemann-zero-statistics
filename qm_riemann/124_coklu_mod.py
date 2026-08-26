"""
124 — ÇOKLU-MOD MIHLAMASI: TAM ARİTMETİK TAYF + FAZ AYRIŞTIRICISI (26 Ağu)
==========================================================================
123'te tek-mod mıhlama 10× termalde bile adyabatik kaldı. Gerçek gaz
~45 çizgiyi p^k-kule yapısı ve ORTAK-ORİJİNLİ SİNÜS FAZLARIYLA birden
taşır. Burada dairesel log-gaza TAM tayf mıhlanır:
  modlar m_q = round(N·log q/L₀) (kule ilişkileri korunur: m₄=2m₂ vb.)
  genlikler Ã_q = (2/π)sin(πτ_q)·Λ(q)/(√q·log p)  [ζ'nın katı yasası]
  → toplam dalga payı ~%85-87 (gerçek rejim!)
İKİ KOŞU (ön-mühürlü ayrıştırıcı):
  A) ARİTMETİK-FAZLI: hepsi ortak-orijinli sinüs (gerçek desen).
  B) KARIŞIK-FAZLI: aynı genlikler, rastgele fazlar.
ÖN-MÜHÜR:
  Ç1  A'da R → eksiye kayarsa ve B'de +1 kalırsa: ANOMALİNİN AJANI FAZ
      TUTARLILIĞI (tanımanın mekanik karşılığı — büyük kapanış).
  Ç2  İkisi de eksiyse: çoklu-mod/etkileşim yeter, faz gerekmez.
  Ç3  İkisi de +1'se: H-çoklu da ölür → tek aday determinizm/iz-formülü.
N=256, L₀=10.37 eşleniği; Metropolis λ=2; 2000 ısınma + 4000 örnekleme.

SONUÇ (26 Ağustos, İLK GEÇİŞ — SONUÇSUZ, DÜRÜST KAYIT):
  A-aritmetik: R_p = +1.50/+0.89/+0.61/−0.77/+0.79 (karışık işaret);
  B-karışık: hepsi pozitif (+0.57..+1.37). ANCAK kabul %11'e düştü,
  A'da σ_η = 0.575 (B: 0.408) — 46 kısıt + faz-uyumlu üst üste binme
  gazı zorluyor; zincirlerin DENGELENDİĞİ ŞÜPHELİ. Bu istatistikle
  Ç1/Ç2/Ç3 ayrımı YAPILAMAZ. Ayrıca katı-genlik (perde uygulanmamış)
  dalga payını 0.234'e şişiriyor (ζ: 0.142) — perdeli-genlik varyantı
  gerekli. → Düzgün ayarlı versiyon (adaptif adım, yakınsama tanıları,
  çok-tohum, perdeli genlik) 125 olarak tayfaya devredildi.
"""

import numpy as np
import time
from sympy import primerange, factorint

TWO_PI = 2 * np.pi
rng = np.random.default_rng(124)
N = 256
L0 = 10.37
LAM = 2.0
STEP = 0.35 * TWO_PI / N

def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

QS = [q for q in pk_list(int(np.exp(0.5 * L0))) if q >= 2]
_agg = {}
for q in QS:
    tau = np.log(q) / L0
    mq = int(round(N * tau))
    if mq < 2 or mq > N // 2:
        continue
    A = (2 / np.pi) * np.sin(np.pi * tau) / (np.sqrt(q))
    _agg[mq] = _agg.get(mq, 0.0) + A       # çakışan çizgiler: hedefler toplanır
MODS = np.array(sorted(_agg))
AMPS = np.array([_agg[m] for m in MODS])
print(f"{len(MODS)} tekil mod mıhlanıyor ({len(QS)} çizgiden; çakışanlar "
      f"birleşik); dalga payı ≈ {0.5*np.sum(AMPS**2):.4f} (ζ: ~0.142; "
      f"katı-genlik, perde yok — bilinçli)", flush=True)

def kos(fazlar, n_burn=2000, n_samp=4000, her=20):
    C0 = N * AMPS / 2 * np.cos(fazlar)
    S0 = N * AMPS / 2 * np.sin(fazlar)
    th = np.sort(rng.uniform(0, TWO_PI, N))
    C = np.array([np.cos(m * th).sum() for m in MODS])
    S = np.array([np.sin(m * th).sum() for m in MODS])
    kabul = 0; toplam = 0
    ornekler = []
    for sw in range(n_burn + n_samp):
        for _ in range(N):
            i = rng.integers(N)
            yeni = th[i] + rng.normal(0, STEP)
            d_es = th - th[i]; d_ye = th - yeni
            d_es[i] = 1.0; d_ye[i] = 1.0
            le = np.log(np.abs(np.sin(0.5 * d_es)))
            ly = np.log(np.abs(np.sin(0.5 * d_ye)))
            le[i] = 0.0; ly[i] = 0.0
            dE = -2.0 * (ly.sum() - le.sum())
            dC = np.cos(MODS * yeni) - np.cos(MODS * th[i])
            dS = np.sin(MODS * yeni) - np.sin(MODS * th[i])
            dV = LAM * (((C + dC - C0)**2 - (C - C0)**2).sum() +
                        ((S + dS - S0)**2 - (S - S0)**2).sum())
            toplam += 1
            if dE + dV < 0 or rng.random() < np.exp(-(dE + dV)):
                th[i] = yeni % TWO_PI
                C += dC; S += dS
                kabul += 1
        if sw >= n_burn and (sw - n_burn) % her == 0:
            ornekler.append(np.sort(th.copy()))
    return ornekler, kabul / toplam

HEDEF = [17, 27, 40, 48, 63]           # m_2, m_3, m_5, m_7, m_13

def R_olc(ornekler):
    Ms = list(MODS)
    C = 2 + 2 * len(Ms)
    XtX = np.zeros((C, C)); Xty = np.zeros(C)
    DATA = []
    for th in ornekler:
        dth = np.diff(np.concatenate([th, [th[0] + TWO_PI]]))
        mid = th + dth / 2
        ds = dth * N / TWO_PI - 1
        cols = [np.ones(N), mid / TWO_PI - 0.5]
        for m in Ms:
            cols += [np.cos(m * mid), np.sin(m * mid)]
        X = np.vstack(cols).T
        XtX += X.T @ X; Xty += X.T @ ds
        DATA.append((X, ds))
    b1 = np.linalg.solve(XtX, Xty)
    s2acc = 0.0; nt = 0
    for X, ds in DATA:
        eta = ds - X @ b1
        s2acc += (eta**2).sum(); nt += len(eta)
    s2 = s2acc / nt
    Xty2 = np.zeros(C)
    for X, ds in DATA:
        eta = ds - X @ b1
        Xty2 += X.T @ (eta**2 - s2)
    b2 = np.linalg.solve(XtX, Xty2)
    out = {}
    for m in HEDEF:
        i = Ms.index(m)
        c1, s1 = b1[2 + 2 * i], b1[2 + 2 * i + 1]
        A1 = np.hypot(c1, s1); ph = np.arctan2(s1, c1)
        P = b2[2 + 2 * i] * np.cos(ph) + b2[2 + 2 * i + 1] * np.sin(ph)
        out[m] = (A1, P / (2 * s2 * A1))
    return out, np.sqrt(s2)

for ad, fazlar in [("A-aritmetik", np.full(len(MODS), -np.pi / 2)),
                   ("B-karışık", rng.uniform(0, TWO_PI, len(MODS)))]:
    t0 = time.time()
    orn, acc = kos(fazlar)
    out, se = R_olc(orn)
    print(f"\n[{ad}] kabul={acc:.2f} σ_η={se:.3f} ({time.time()-t0:.0f} sn)")
    print(f"{'m':>4} {'A1':>7} {'R_p':>7}")
    for m in HEDEF:
        A1, Rp = out[m]
        print(f"{m:>4} {A1:>7.4f} {Rp:>+7.3f}", flush=True)
print("\nGERÇEK GAZ: R_p ≈ −0.85..−0.94 (içsel −2); 123 tek-mod: +0.95")
