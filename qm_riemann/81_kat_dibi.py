"""
81 — KAT DİBİ SİS BOĞAZI: DİREKT + GÖRÜNTÜ KARIŞIMI (20 Ağustos 2026)
==========================================================================
79 + 80'in bütün ipleri buraya işaret ediyordu: kat dibinde (τ ≈ 0.5)
çıplak okuma B = −0.344, giydirilmiş −0.124 (×2.77), ve görüntü
frekansları tam burada direkt frekanslarla çakışıyor.

Üç test:
  T1  GİYDİRME MERDİVENİ: band-asal B̂'si taban katmanlarıyla —
      L0 yalnız band | L1 +P11 | L2 +yalnız p^k kuvvetleri | L3 +ikisi |
      L4 +tam direkt taban (tüm p^k ≤ e^{0.45L}); kontrollü VE kontrolsüz
      (g_u, g_u² kovaryatları mekanizmanın parçası mı?).
  T2  GÖRÜNTÜ REGRESÖRLERİ: her band asalı p için durağan-faz bandı
      m ∈ [0.8, 1.25]·m*, m* = (t̄/2π)/p üzerinden ayna toplamı
      R_p = Σ d(m)/√m · cos/sin(2θ(t) − t·log m)  (74'ün makinesi,
      İLK KEZ max kanalında). Kontroller: kaydırılmış-bant (2.5·m* —
      durağan faz aritmetik-olmayan frekansa kilitli) ve yarım-tamsayı
      ızgara plasebosu. v kanalında aynı kolonlar (sessizlik öngörüsü).
      Soru: R^stat içerik taşıyor mu, ve varlığı band-asal B̂'sini
      oynatıyor mu (karışım)?
  T3  MEKANİZMA: giydirme örgülü-örnekleme aracılığıyla mı? Küçük-asal
      kolonları ile band-asal kolonlarının Gram-R²'si gerçek gap-ortaları
      vs RvM-pürüzsüz ızgara. Gerçekte var, pürüzsüzde yoksa → örnekleme
      desenine işlenmiş asal dalgaları kanalları çiftliyor.
  T4  SAHİCİ ARKA ODA: L5 (band-altı p^k TAM taban) profili, 36+41+55.

SONUÇLAR (20 Ağustos sabahı):
  T1: kuvvetler MASUM (L2≈L0), küçük asallar suçlu; merdiven L4→L5
      yakınsıyor: çıplak −0.354 → tam-taban −0.042±0.002. Kontrolsüz
      sütun işaret değiştirir (+1.56→+0.90): kontroller yük taşıyıcı.
  T3: Gram R² gerçek ortalarda 7e-4, pürüzsüz ızgarada 0 → örnekleme
      ızgarası (v-alanı) asal dalgalarını taşıyor; p^{-1/2} normalizasyonu
      sızıntıyı √p ile büyütüyor → giydirme = atlanmış-değişken transferi.
      ÇIPLAK OKUMANIN ~%88'İ HAYALET; sahici oda küçük ama 18σ gerçek.
  T2: durağan-bant görüntü içeriği zayıf (cos +0.0026±0.0005, 5.5σ);
      kaydırılmış-bant ve yarım-ızgara kontrolleri NULL DEĞİL (30σ/13σ,
      hepsi saf-cos) → 2θ-toplam ailesi ortak bir S-kaynaklı bileşen
      paylaşıyor; T2 tasarımı v2 ister (karıştırılmış-gap vekil ızgara?).
      v kanalı öngörüldüğü gibi sessiz (⟨z²⟩=3.4).
  T4: B_tam(L) ≈ −0.05'te yatay (−0.06…−0.026, hafif iniş olabilir),
      36/41 sınırında sıçrama yok.
"""

import numpy as np
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A_CG = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
P11 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]

def theta(t):
    return t / 2 * np.log(t / TWO_PI) - t / 2 - np.pi / 8 + 1 / (48 * t)

def unfold(gaps, amps, tmid):
    Lw = np.log(tmid / TWO_PI)
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A_CG * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    return np.log(a_u), np.log(g_u), g_u, float(Lw.mean())

def reg_bhat(y, tmid, band, extra_freqs, g_u=None, extra_cols=None):
    cols = [np.ones_like(y)]
    if g_u is not None:
        cols += [g_u, g_u**2]
    for q in extra_freqs:
        arg = tmid * np.log(q)
        cols += [np.cos(arg), np.sin(arg)]
    i0 = len(cols)
    for p in band:
        arg = tmid * np.log(p)
        cols += [np.cos(arg), np.sin(arg)]
    if extra_cols:
        cols += extra_cols
    X = np.vstack(cols).T
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    se = np.sqrt((y - X @ b).var() * np.diag(np.linalg.inv(X.T @ X)))
    ws = np.array([b[i0 + 2*i] / p**-0.5 for i, p in enumerate(band)])
    ss = np.array([se[i0 + 2*i] / p**-0.5 for i, p in enumerate(band)])
    Bh = np.sum(ws / ss**2) / np.sum(1 / ss**2)
    Be = 1 / np.sqrt(np.sum(1 / ss**2))
    return Bh, Be, b, se, X

def pk_list(lim, primes_too=False):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p if primes_too else p * p
        while q <= lim:
            out.append(q)
            q *= p
    return sorted(out)

WNDS = []
d41 = np.load(HERE / "41_bigT_windows.npz")
for k in sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1])):
    WNDS.append((d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"]))
d55 = np.load(HERE / "55_win_1e+08.npz")
WNDS.append((d55["gaps"], d55["amps"], d55["tmid"]))

# ============ T1: GİYDİRME MERDİVENİ ============
print("T1 — GİYDİRME MERDİVENİ (band [0.505,0.55]; hücre: kombine B̂)")
print(f"{'katman':>34} {'kontrollü':>12} {'kontrolsüz':>12}")
merdiven = {}
for isim, kucuk_al in [
    ("L0: yalnız band", lambda L: []),
    ("L1: + P11 küçük asallar", lambda L: [p for p in P11 if np.log(p)/L < 0.45]),
    ("L2: + yalnız p^k kuvvetleri", lambda L: [q for q in pk_list(np.exp(0.45*L)) ]),
    ("L3: + P11 + kuvvetler", lambda L: sorted(set(
        [p for p in P11 if np.log(p)/L < 0.45] + pk_list(np.exp(0.45*L))))),
    ("L4: + TAM direkt taban (p^k tam)", lambda L: pk_list(np.exp(0.45*L), primes_too=True)),
]:
    res = {}
    for ctl in [True, False]:
        num = den = 0.0
        for gaps, amps, tmid in WNDS:
            ya, yg, g_u, L = unfold(gaps, amps, tmid)
            band = [p for p in primerange(int(np.exp(0.505*L)), int(np.exp(0.55*L)))][:10]
            Bh, Be, *_ = reg_bhat(ya, tmid, band, kucuk_al(L),
                                  g_u=(g_u if ctl else None))
            num += Bh / Be**2; den += 1 / Be**2
        res[ctl] = (num / den, 1 / np.sqrt(den))
    merdiven[isim] = res
    print(f"{isim:>34} {res[True][0]:>+9.4f}±{res[True][1]:.3f} "
          f"{res[False][0]:>+9.4f}±{res[False][1]:.3f}")

# ============ T2: GÖRÜNTÜ REGRESÖRLERİ ============
print("\nT2 — GÖRÜNTÜ REGRESÖRLERİ (durağan-faz ayna toplamları, z-skorlu)")
M_MAX = 12000
dcount = np.zeros(M_MAX + 1, dtype=np.int64)
for i in range(1, M_MAX + 1):
    dcount[i::i] += 1

def ayna_toplami(tmid, th2, m_lo, m_hi):
    Rc = np.zeros_like(tmid); Rs = np.zeros_like(tmid)
    for m in range(int(m_lo), int(m_hi) + 1):
        wgt = dcount[m] / np.sqrt(m)
        arg = th2 - tmid * np.log(m)
        Rc += wgt * np.cos(arg); Rs += wgt * np.sin(arg)
    return Rc, Rs

def yarim_toplam(tmid, th2, m_lo, m_hi):
    Rc = np.zeros_like(tmid); Rs = np.zeros_like(tmid)
    for m in range(int(m_lo), int(m_hi) + 1):
        wgt = dcount[m] / np.sqrt(m + 0.5)
        arg = th2 - tmid * np.log(m + 0.5)
        Rc += wgt * np.cos(arg); Rs += wgt * np.sin(arg)
    return Rc, Rs

def z(x):
    sd = x.std()
    return (x - x.mean()) / (sd if sd > 0 else 1.0)

acc = {"stat_c": [], "stat_s": [], "off_c": [], "off_s": [],
       "half_c": [], "half_s": [], "v_stat_c": [], "v_stat_s": []}
B_degisim = []
for gaps, amps, tmid in WNDS:
    ya, yg, g_u, L = unfold(gaps, amps, tmid)
    th2 = 2 * theta(tmid)
    tbar = float(tmid.mean())
    band = [p for p in primerange(int(np.exp(0.505*L)), int(np.exp(0.55*L)))][:10]
    tam = pk_list(np.exp(0.45*L), primes_too=True)
    ex_stat, keys = [], []
    for p in band:
        mstar = (tbar / TWO_PI) / p
        if 1.25 * mstar > M_MAX:
            continue
        Rc, Rs = ayna_toplami(tmid, th2, max(1, 0.8*mstar), 1.25*mstar)
        ex_stat += [z(Rc), z(Rs)]; keys.append(p)
    mmed = (tbar / TWO_PI) / np.median(band)
    off_ok = 2.9 * mmed <= M_MAX
    extra = list(ex_stat)
    if off_ok:
        Oc, Os = ayna_toplami(tmid, th2, max(1, 2.2*mmed), 2.9*mmed)
        extra += [z(Oc), z(Os)]
    Hc, Hs = yarim_toplam(tmid, th2, max(1, 0.8*mmed), min(1.25*mmed, M_MAX))
    extra += [z(Hc), z(Hs)]
    # R'siz ve R'li band-B̂ (karışım testi)
    B0_, e0_, *_ = reg_bhat(ya, tmid, band, tam, g_u=g_u)
    B1_, e1_, b, se, X = reg_bhat(ya, tmid, band, tam, g_u=g_u, extra_cols=extra)
    B_degisim.append((L, B0_, e0_, B1_, e1_))
    i0 = X.shape[1] - len(extra)
    for j, p in enumerate(keys):
        acc["stat_c"].append((b[i0+2*j], se[i0+2*j]))
        acc["stat_s"].append((b[i0+2*j+1], se[i0+2*j+1]))
    if off_ok:
        acc["off_c"].append((b[-4], se[-4])); acc["off_s"].append((b[-3], se[-3]))
    acc["half_c"].append((b[-2], se[-2])); acc["half_s"].append((b[-1], se[-1]))
    # v kanalı: aynı kolonlar
    _, _, bv, sev, Xv = reg_bhat(yg, tmid, band, tam, extra_cols=extra)
    i0v = Xv.shape[1] - len(extra)
    for j, p in enumerate(keys):
        acc["v_stat_c"].append((bv[i0v+2*j], sev[i0v+2*j]))
        acc["v_stat_s"].append((bv[i0v+2*j+1], sev[i0v+2*j+1]))

def komb(key):
    a = np.array(acc[key])
    m = np.sum(a[:,0]/a[:,1]**2)/np.sum(1/a[:,1]**2)
    e = 1/np.sqrt(np.sum(1/a[:,1]**2))
    x2 = float(np.mean((a[:,0]/a[:,1])**2))
    return m, e, x2

for isim, kc, ks in [("R^stat (w)", "stat_c", "stat_s"),
                     ("R^off  (w)", "off_c", "off_s"),
                     ("R^half (w)", "half_c", "half_s"),
                     ("R^stat (v)", "v_stat_c", "v_stat_s")]:
    mc, ec, xc = komb(kc); ms, es, xs = komb(ks)
    print(f"  {isim}: cos = {mc:+.5f}±{ec:.5f} (⟨z²⟩={xc:.1f}) | "
          f"sin = {ms:+.5f}±{es:.5f} (⟨z²⟩={xs:.1f})")

print("\n  KARIŞIM: R^stat kolonları band-asal B̂'sini oynatıyor mu?")
for L, B0_, e0_, B1_, e1_ in B_degisim:
    print(f"    L={L:5.2f}: R'siz {B0_:+.4f}±{e0_:.4f} → R'li {B1_:+.4f}±{e1_:.4f}")

# ============ T3: MEKANİZMA — ÖRNEKLEME GRAM TESTİ ============
print("\nT3 — GRAM R² (band-asal kolonu ~ küçük-asal kolonları):")
def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8
for wi in [1, 6]:
    gaps, amps, tmid = WNDS[wi]
    ya, yg, g_u, L = unfold(gaps, amps, tmid)
    t0 = tmid[0] - gaps[0] / 2
    kk = np.arange(len(tmid)) + 0.5
    ts = t0 + kk * TWO_PI / np.log(t0 / TWO_PI)
    for _ in range(6):
        fdel = rvm_N(ts) - rvm_N(t0) - kk
        ts = ts - fdel / (np.log(ts / TWO_PI) / TWO_PI)
    band = [p for p in primerange(int(np.exp(0.505*L)), int(np.exp(0.55*L)))][:10]
    kucuk = [p for p in P11 if np.log(p)/L < 0.45]
    for etiket, tg in [("gerçek ortalar", tmid), ("RvM pürüzsüz", ts)]:
        Xs = np.vstack([f(tg * np.log(q)) for q in kucuk
                        for f in (np.cos, np.sin)]).T
        r2s = []
        for p in band:
            c = np.cos(tg * np.log(p))
            bb, *_ = np.linalg.lstsq(Xs, c - c.mean(), rcond=None)
            r2s.append(1 - ((c - c.mean() - Xs @ bb)**2).sum()
                       / ((c - c.mean())**2).sum())
        print(f"  L={L:5.2f} {etiket:>14}: ⟨R²⟩ = {np.mean(r2s):.5f}")

# ============ T4: SAHİCİ ARKA ODA PROFİLİ (L5 tam-taban) ============
print("\nT4 — SAHİCİ ARKA ODA (L5: band-altı p^k tam taban):")
d36 = np.load(HERE / "36_T100k.npz")
edges = np.geomspace(d36["t_mid"][0], d36["t_mid"][-1] * 1.0001, 13)
LOW = []
for i in range(12):
    m = (d36["t_mid"] >= edges[i]) & (d36["t_mid"] < edges[i + 1])
    if m.sum() >= 3000:
        LOW.append((d36["intervals"][m], d36["max_amps"][m], d36["t_mid"][m]))
for src, WL in [("36", LOW), ("41/55", WNDS)]:
    for gaps, amps, tmid in WL:
        ya, yg, g_u, L = unfold(gaps, amps, tmid)
        band = [p for p in primerange(int(np.exp(0.505*L)), int(np.exp(0.55*L)))][:10]
        if len(band) < 2:
            continue
        t5 = pk_list(np.exp(0.505*L), primes_too=True)
        B5, e5, *_ = reg_bhat(ya, tmid, band, t5, g_u=g_u)
        print(f"  {src:>5} L={L:5.2f}: B_tam = {B5:+.4f} ± {e5:.4f}")
