"""
132b — ÇİFT-SINIF DİSEKSİYONU: R HANGİ BLOKTAN GELİYOR? (27 Ağustos)
==========================================================================
132a'nın SERT BULGUSU: δu₂ = −S'(t)u₁/ρ̄ λ=1'de örgüyü YIKIYOR
(3429 kesişme, σ_η 0.30→0.44, R çöküyor) — çünkü açılım parametresi
BİRDEN BÜYÜK: rms(S'/ρ̄) = 1.172 (kuyruk 1.113 / taban 0.419),
rms(δu₂)/rms(u₁) = 1.07. Lagrange serisi yakınsama yarıçapının
DIŞINDA; kalem 2. mertebede KESEMEZ. Öyleyse diseksiyon iki katta:

KAT-1 — TAYLOR BLOKLARI, DOĞRUSAL-TEPKİ REJİMİNDE (λ küçük):
  δu₂ = (1/ρ̄²) Σ_{Q,Q'} a_Q ω_Q a_{Q'} cos(ω_Q t) sin(ω_{Q'} t)
        = −S'(t)·u₁(t)/ρ̄   → blok başına ÇARPIM: −S'_A·u₁_B/ρ̄
  Örgü: z = t + u₁(tam) + λ·δu₂^{A×B}; λ=0.25'te örgü sağlam
  (kesişme yok) ⟹ eğim  m_B = [R_nn(λ) − R_nn(0)]/λ  toplanabilir.
  Bloklar (kuyruk τ₀=logQ/L₀>0.52 → 2109 çizgi; taban ≤0.52 → 19):
    B1 kuyruk×kuyruk | B2a kuyruk(türev)×taban(dalga)
    B2b taban(türev)×kuyruk(dalga) | B3 taban×taban
    B4  ölçülen p'yi içeren çiftler (alt-sınıflar B4a=kuyruk×p,
        B4b=p×kuyruk) — 130'un "dalga×alan" adayı, p başına ayrı örgü.
  λ-merdiveni ayrıca serinin nerede koptuğunu gösterir.

KAT-2 — YERDEĞİŞTİRME SINIFLARI (RESUMLANMIŞ, ıraksama yok):
  Kesin ilişki: z = t + u,  u = −S(z)/ρ̄ ⟹ ds = −x/(1+x), x = S'(z)/ρ̄.
  Taylor yerine ARGÜMANI böl — hangi sınıf YERİNDEN OYNAMIŞ noktada
  örnekleniyorsa etki oradan doğar:
    C0  u = −S(t+u)/ρ̄                              (=A3, referans)
    C1  u = −[S_kuy(t+u) + S_tab(t)]/ρ̄     (yalnız KUYRUK yer değiştirir)
    C2  u = −[S_kuy(t) + S_tab(t+u)]/ρ̄     (yalnız TABAN yer değiştirir)
    C3  u = −[S_kuy(t+u_tab) + S_tab(t)]/ρ̄ (kuyruk YALNIZ taban-dalgasınca
        oynatılır; u_tab = −S_tab(t)/ρ̄ — iterasyonsuz, tam resumlanmış)
    C4p u = −[S_kuy(t+u_p) + S_tab(t)]/ρ̄   (yalnız TEK p çizgisinin dalgası)
  C3/C4 kalemin ıraksamayan nesnesidir: Taylor kesilmez, yerdeğiştirme
  argümanın içinde tutulur.

FAZ ÇAPRAZI (130-F2, 2. katta): baskın sınıfın katkısı kuyruk fazları
karıştırılınca ÖLMELİ (koherans-esaslı ⟹ faz-duyarlı).

ÖN-MÜHÜR:
  P1  B1/C1 baskınsa → kuyruk×kuyruk (126'nın ilk defteri).
  P2  B2/C3 baskınsa → kuyruk×taban-dalgası (130'un tezi); türevin
      hangi sınıfta durduğu (B2a mı B2b mi) kalem için belirleyici.
  P3  B4/C4p, B2'nin büyük kısmını yiyorsa yasa TEK ÇİZGİYE indirgenir
      (kalem için en ucuz sonuç: p-dalgası özel).

SONUÇ (27 Ağustos) — P2+P3 ✓: YASA "KUYRUĞUN p-DALGASINCA OYNATILMIŞ
ÖRNEKLENİŞİ"NDEN DOĞUYOR; TAYLOR KATI YALNIZ SÜRÜCÜYÜ İŞARET EDİYOR.

TANI: rms(S'/ρ̄) = 1.172 (kuyruk 1.113 / taban 0.419, max 4.24);
rms δu₂: B1 0.168 | B2a 0.256 | B2b 0.070 | B3 0.103 | TAM 0.333
(ḡ=0.898, rms u₁=0.310) — "düzeltme" öncü terim boyunda.

KAT-1a λ MERDİVENİ (tam δu₂): R_nn(p=2) −0.33 → −0.31 → −0.26 → −0.26
→ −0.12 (λ=0/0.1/0.25/0.5/1) ve kesişme 916 → 3429. Terim λ büyüdükçe
R'yi SIFIRA sürüyor; seri hedeften UZAKLAŞTIRIYOR. (Uyarı: λ=0 tabanı
zaten 916 kesişmeli — KAT-1 kırık örgü üstünde ölçülüyor, artefakt
şüphesi kayıtlı; bu yüzden yalnız EĞİM/atıf için okunur.)

KAT-1b BLOK EĞİMLERİ m = ΔR_nn/λ (λ=0.25; + = hedeften uzak):
   p    B1      B2a      B2b      B3   | Σblok    TAM
   2  +0.601  −0.891  +0.884  −0.124  | +0.471  +0.287
   3  +0.345  −0.345  +0.892  −0.014  | +0.879  +0.836
   5  +0.261  +0.013  +0.927  +0.240  | +1.442  +1.459
   7  +0.041  −0.270  +0.959  +0.371  | +1.102  +1.208
  TOPLANABİLİRLİK ✓ (Σblok ≈ TAM) — blok ayrıştırması meşru.
  TEK NEGATİF (hedefe DOĞRU süren) blok: B2a = kuyruk(türev)×taban(dalga).
  En büyük pozitif (zarar) blok: B2b = taban(türev)×kuyruk(dalga) —
  düzgün +0.88..+0.96; ıraksamayı taşıyan da bu.
KAT-1c B4: m[B4a = kuyruk×p] = −0.80/−0.66/−0.58/−0.66 — HER p'DE
  NEGATİF ve p=3,5,7'de tüm B2a'dan DAHA negatif (öteki taban çizgileri
  kısmen götürüyor). m[B4b = p×kuyruk] = +0.49..+0.73 (zarar tarafı).
  ⟹ İŞARET-DOĞRU SÜRÜCÜ TEK BİR ALT-SINIF: (Q ∈ kuyruk, Q' = p).

KAT-2 YERDEĞİŞTİRME SINIFLARI (resumlanmış — ASIL ATIF):
  C0 tam öz-tutarlı : −1.811/−1.782/−1.657/−1.536  (σ_η 0.163, kes. 6)
  C1 yalnız KUYRUK yer değiştirir: −1.435/−1.649/−1.728/−1.845
     = C0'ın %79/%93/%104/%120 — YASA YAŞIYOR.
  C2 yalnız TABAN yer değiştirir : +0.452/+0.275/+0.243/+0.273
     — YASA ÖLÜ, İŞARET TERS. ⟹ Etki tamamen "kuyruk yer değiştirmiş
     noktada örnekleniyor" kanalında.
  C3 kuyruk yalnız TABAN-DALGASINCA oynatılır (iterasyonsuz):
     −1.170/−1.388/−1.366/−1.525 = C0'ın %65/%78/%82/%99.
  C4p kuyruk yalnız TEK p ÇİZGİSİNİN dalgasınca oynatılır:
     −1.225/−1.177/−1.245/−1.416 = C0'ın %68/%66/%75/%92.
  ⟹ P3 ✓: TEK ÇİZGİNİN yerdeğiştirmesi yasanın üçte-ikisinden fazlasını
  tek başına veriyor. 130'un tezi (kuyruk × p-dalgası) DOĞRULANDI.

FAZ ÇAPRAZI: FZ-C0 (tam lab, kuyruk fazı karışık) R_nn = −0.40/−0.42/
  −0.29/−0.31, R_p = +0.28/+0.21/+0.22/+0.14 — 130-F2 BİREBİR tekrar
  (etki ölüyor) ✓. FZ-C3 %55-65 koruyor; FZ-B2a eğimi karıştırınca
  −1.35/−1.14/−1.18/−1.26'ya (daha da negatif) çıkıyor. ⟹ FAZ
  DUYARLILIĞI ÇEKİRDEKTE DEĞİL, ÖZ-TUTARLI SABİT NOKTADA: koherans
  koşulu "tüm çizgiler aynı fazda" olmasını, ızgaranın kuyruğu kısmen
  soğurmasıyla birlikte gerektiriyor (σ_η 0.163→0.220). KALEM İÇİN
  AÇIK UÇ — ve 132c'nin çekirdek kuramının en keskin sınavı.
"""

import numpy as np
import time
from sympy import primerange, factorint

TWO_PI = 2 * np.pi
L0 = 7.0
NZ = 40000
PS = [2, 3, 5, 7]
LAM_SCAN = [0.1, 0.25, 0.5, 1.0]
LAM_BLOK = 0.25
rng = np.random.default_rng(132)


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

QCAP = int(np.exp(1.4 * L0))
QS = pk_list(QCAP)
QARR = np.array(QS)
LOGQ = np.log(QARR)
LAMv = np.array([float(np.log(list(factorint(q).items())[0][0])) for q in QS])
W0 = LAMv / (np.pi * np.sqrt(QARR) * LOGQ)
TAU0 = LOGQ / L0
KUY = TAU0 > 0.52
TAB = ~KUY
HEP = np.ones(len(QS), bool)
FAZ0 = np.zeros(len(QS))
FAZ_K = FAZ0.copy()
FAZ_K[KUY] = rng.uniform(0, TWO_PI, int(KUY.sum()))


def S_field(t, mask=HEP, faz=FAZ0, chunk=8000):
    """S_A(t) = −Σ_{Q∈A} a_Q sin(ω_Q t + φ_Q)"""
    w = W0[mask]; lq = LOGQ[mask]; ph = faz[mask]
    s = np.zeros_like(t)
    if w.size == 0:
        return s
    for a in range(0, len(t), chunk):
        tt = t[a:a + chunk]
        s[a:a + chunk] = -(np.sin(np.outer(tt, lq) + ph[None, :])
                           * w[None, :]).sum(axis=1)
    return s


def Sp_field(t, mask=HEP, faz=FAZ0, chunk=8000):
    """S'_A(t) = −Σ_{Q∈A} a_Q ω_Q cos(ω_Q t + φ_Q)"""
    w = (W0 * LOGQ)[mask]; lq = LOGQ[mask]; ph = faz[mask]
    s = np.zeros_like(t)
    if w.size == 0:
        return s
    for a in range(0, len(t), chunk):
        tt = t[a:a + chunk]
        s[a:a + chunk] = -(np.cos(np.outer(tt, lq) + ph[None, :])
                           * w[None, :]).sum(axis=1)
    return s


def olc(z):
    g = np.diff(z)
    m = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(m / TWO_PI)
    Lb = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    qs_fit = pk_list(min(int(np.exp(0.52 * Lb)), 720))
    freqs = [np.log(q) for q in qs_fit]
    tt = (m - m.mean()) / (m[-1] - m[0])
    cols = [np.ones_like(m), tt, tt**2]
    for f in freqs:
        cols += [np.cos(f * m), np.sin(f * m)]
    X = np.vstack(cols).T
    b1, *_ = np.linalg.lstsq(X, ds, rcond=None)
    eta = ds - X @ b1
    s2 = float((eta**2).mean())
    b2, *_ = np.linalg.lstsq(X, eta**2 - s2, rcond=None)
    ee = eta[:-1] * eta[1:]
    c1 = float(ee.mean())
    mm = 0.5 * (m[:-1] + m[1:])
    tt2 = (mm - mm.mean()) / (mm[-1] - mm[0])
    cols2 = [np.ones_like(mm), tt2, tt2**2]
    for f in freqs:
        cols2 += [np.cos(f * mm), np.sin(f * mm)]
    X2 = np.vstack(cols2).T
    b3, *_ = np.linalg.lstsq(X2, ee - c1, rcond=None)
    out = {}
    for p in PS:
        f = np.log(p)
        i = qs_fit.index(p)
        cA, sA = b1[3 + 2 * i], b1[3 + 2 * i + 1]
        A1 = np.hypot(cA, sA); ph = np.arctan2(sA, cA)
        Rp = (b2[3 + 2 * i] * np.cos(ph) +
              b2[3 + 2 * i + 1] * np.sin(ph)) / (2 * s2 * A1)
        Rn = (b3[3 + 2 * i] * np.cos(ph) +
              b3[3 + 2 * i + 1] * np.sin(ph)) / (2 * c1 * A1)
        out[p] = (f / Lb, Rp, Rn, A1)
    return out, np.sqrt(s2), c1 / s2


def hucre(ad, zraw, sessiz=False):
    kes = int((np.diff(zraw) <= 0).sum())
    out, se, c1r = olc(np.sort(zraw))
    if not sessiz:
        print(f"\n[{ad}] kesişme={kes} σ_η={se:.4f} c₁/σ²={c1r:+.2f}")
        print(f"{'p':>3} {'τ':>6} {'R_p':>7} {'R_nn':>7} {'−2cosπτ':>8} {'A₁':>9}")
        for p in PS:
            tau, Rp, Rn, A1 = out[p]
            print(f"{p:>3} {tau:>6.3f} {Rp:>+7.3f} {Rn:>+7.3f} "
                  f"{-2*np.cos(np.pi*tau):>+8.3f} {A1:>9.2e}", flush=True)
    return out, kes


print(f"lab: {len(QS)} çizgi ({int(KUY.sum())} kuyruk / {int(TAB.sum())} taban)"
      f"; taban: {[int(q) for q in QARR[TAB]]}", flush=True)

u1 = -S_field(tg) / rho
u1_k = -S_field(tg, KUY) / rho
u1_t = -S_field(tg, TAB) / rho
Sp_h = Sp_field(tg)
Sp_k = Sp_field(tg, KUY)
Sp_t = Sp_field(tg, TAB)
gbar = TWO_PI / L0
print(f"\nAÇILIM TANISI: ḡ={gbar:.4f}  rms(u₁)={np.sqrt((u1**2).mean()):.4f}")
for ad, sp in [("hepsi", Sp_h), ("kuyruk", Sp_k), ("taban", Sp_t)]:
    print(f"  rms(S'_{ad}/ρ̄) = {np.sqrt(((sp/rho)**2).mean()):.4f}  "
          f"max = {np.abs(sp/rho).max():.3f}", flush=True)

DU2 = {"B1": -Sp_k * u1_k / rho, "B2a": -Sp_k * u1_t / rho,
       "B2b": -Sp_t * u1_k / rho, "B3": -Sp_t * u1_t / rho,
       "TAM": -Sp_h * u1 / rho}
for k, v in DU2.items():
    print(f"  rms(δu₂[{k}]) = {np.sqrt((v**2).mean()):.4f}", flush=True)

# ================= KAT-1a: λ MERDİVENİ (tam δu₂) =====================
print("\n########## KAT-1a: λ MERDİVENİ (tam δu₂) ##########", flush=True)
base, _ = hucre("λ=0  (yalnız u₁ — taban çizgisi)", tg + u1)
for lam in LAM_SCAN:
    hucre(f"λ={lam}  z = t+u₁+λ·δu₂(tam)", tg + u1 + lam * DU2["TAM"])

# ================= KAT-1b: BLOK ATFI (λ=0.25) ========================
print(f"\n########## KAT-1b: BLOK ATFI (doğrusal tepki, λ={LAM_BLOK}) ##########",
      flush=True)
BL = ["B1", "B2a", "B2b", "B3"]
ETI = {"B1": "kuyruk×kuyruk", "B2a": "kuyruk(türev)×taban(dalga)",
       "B2b": "taban(türev)×kuyruk(dalga)", "B3": "taban×taban",
       "TAM": "TAM (kontrol)"}
egim = {}
for b in BL + ["TAM"]:
    o, kes = hucre(f"{b}  {ETI[b]}  λ={LAM_BLOK}",
                   tg + u1 + LAM_BLOK * DU2[b])
    egim[b] = {p: (o[p][2] - base[p][2]) / LAM_BLOK for p in PS}

print("\n=== BLOK EĞİMLERİ  m = ΔR_nn/λ  ve TAM'a payı ===", flush=True)
print(f"{'p':>3} " + " ".join(f"{b:>9}" for b in BL) +
      f" {'Σblok':>9} {'TAM':>9}")
for p in PS:
    v = [egim[b][p] for b in BL]
    print(f"{p:>3} " + " ".join(f"{x:>+9.3f}" for x in v) +
          f" {sum(v):>+9.3f} {egim['TAM'][p]:>+9.3f}", flush=True)
print(f"{'':>3} " + " ".join(f"{b:>9}" for b in BL) + "   (TAM'ın % payı)")
for p in PS:
    tot = egim["TAM"][p]
    print(f"{p:>3} " + " ".join(f"{egim[b][p]/tot:>8.1%}" for b in BL), flush=True)

# ================= KAT-1c: B4 (p'yi içeren çiftler) ==================
print(f"\n########## KAT-1c: B4 'dalga×alan' alt-sınıfı (λ={LAM_BLOK}) ##########",
      flush=True)
print(f"{'p':>3} {'m[B4a kuy×p]':>13} {'m[B4b p×kuy]':>13} {'m[B4 p-dahil]':>14} "
      f"{'m[B2a]':>9} {'m[B2b]':>9} {'m[TAM]':>9}", flush=True)
for p in PS:
    ip = QS.index(p)
    MP = np.zeros(len(QS), bool); MP[ip] = True
    u1_p = -S_field(tg, MP) / rho
    Sp_p = Sp_field(tg, MP)
    d4a = -Sp_k * u1_p / rho
    d4b = -Sp_p * u1_k / rho
    d4 = -(Sp_p * u1 + Sp_h * u1_p - Sp_p * u1_p) / rho
    r = []
    for d in (d4a, d4b, d4):
        o, _ = hucre("", tg + u1 + LAM_BLOK * d, sessiz=True)
        r.append((o[p][2] - base[p][2]) / LAM_BLOK)
    print(f"{p:>3} {r[0]:>+13.3f} {r[1]:>+13.3f} {r[2]:>+14.3f} "
          f"{egim['B2a'][p]:>+9.3f} {egim['B2b'][p]:>+9.3f} "
          f"{egim['TAM'][p]:>+9.3f}", flush=True)

# ================= KAT-2: YERDEĞİŞTİRME SINIFLARI ====================
print("\n########## KAT-2: YERDEĞİŞTİRME SINIFLARI (resumlanmış) ##########",
      flush=True)


def lab_split(mask_disp, n_iter=5, faz=FAZ0):
    """u = −[S_{mask}(t+u) + S_{~mask}(t)]/ρ̄, sönümlü Picard (131 şeması)"""
    sabit = S_field(tg, ~mask_disp, faz)
    u = np.zeros_like(tg)
    for _ in range(n_iter):
        u = 0.5 * u + 0.5 * (-(S_field(tg + u, mask_disp, faz) + sabit) / rho)
    return tg + u


ts = time.time()
hucre("C0  u=−S(t+u)/ρ̄  (=A3 referansı)", lab_split(HEP))
hucre("C1  yalnız KUYRUK yer değiştirir", lab_split(KUY))
hucre("C2  yalnız TABAN yer değiştirir", lab_split(TAB))
zc3 = tg + (-(S_field(tg + u1_t, KUY) + S_field(tg, TAB)) / rho)
hucre("C3  kuyruk YALNIZ taban-dalgasınca oynatılır (iterasyonsuz)", zc3)
print(f"  (KAT-2 ana hücreler: {time.time()-ts:.0f} sn)", flush=True)

print("\n=== C4p: yalnız TEK p çizgisinin dalgası kuyruğu oynatır ===", flush=True)
print(f"{'p':>3} {'τ':>6} {'R_p':>8} {'R_nn':>8} {'−2cosπτ':>9} {'kesişme':>8}",
      flush=True)
for p in PS:
    ip = QS.index(p)
    MP = np.zeros(len(QS), bool); MP[ip] = True
    u1_p = -S_field(tg, MP) / rho
    z = tg + (-(S_field(tg + u1_p, KUY) + S_field(tg, TAB)) / rho)
    o, kes = hucre("", z, sessiz=True)
    tau, Rp, Rn, _ = o[p]
    print(f"{p:>3} {tau:>6.3f} {Rp:>+8.3f} {Rn:>+8.3f} "
          f"{-2*np.cos(np.pi*tau):>+9.3f} {kes:>8d}", flush=True)

# ================= FAZ ÇAPRAZI =======================================
print("\n########## FAZ ÇAPRAZI (kuyruk fazları karışık, taban sabit) ##########",
      flush=True)
u1_kf = -S_field(tg, KUY, FAZ_K) / rho
Sp_kf = Sp_field(tg, KUY, FAZ_K)
u1_f = u1_kf + u1_t
basef, _ = hucre("FZ-λ=0  (yalnız u₁, kuyruk fazı karışık)", tg + u1_f)
for ad, d in [("FZ-B2a kuyruk×taban", -Sp_kf * u1_t / rho),
              ("FZ-B2b taban×kuyruk", -Sp_t * u1_kf / rho),
              ("FZ-B1  kuyruk×kuyruk", -Sp_kf * u1_kf / rho)]:
    o, _ = hucre(f"{ad}  λ={LAM_BLOK}", tg + u1_f + LAM_BLOK * d)
    print("  eğim m = " + " ".join(
        f"p{p}:{(o[p][2]-basef[p][2])/LAM_BLOK:+.3f}" for p in PS), flush=True)
hucre("FZ-C3  kuyruk(faz karışık) taban-dalgasınca oynatılır",
      tg + (-(S_field(tg + u1_t, KUY, FAZ_K) + S_field(tg, TAB)) / rho))
hucre("FZ-C0  tam öz-tutarlı, kuyruk fazı karışık (130-F2 tekrarı)",
      lab_split(HEP, faz=FAZ_K))
