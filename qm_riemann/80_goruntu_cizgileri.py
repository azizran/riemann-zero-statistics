"""
80 — GÖRÜNTÜ ÇİZGİLERİ: 2θ SPEKTROSKOPİSİ (19 Ağustos 2026)
==========================================================================
Fısıltı (çizim seansı): aynadan yansıyan dalganın frekansı da aynalanır —
görüntü kaynak τ′ = 1 − τ'da, frekansı L − log m (KAYAR: chirp).
Düz cos(t·ω) kulağı bunu duyamaz; 2θ gözlüğü gerekir:
    X_m(t) = cos(2θ(t) − t·log m),  θ = RS faz fonksiyonu.
Bunlar AFE ayna-toplamının çapraz terimleri (74'ün noktasal |Z|²'de
bulduğu chirp'lerin ta kendisi) — İLK KEZ max/gap KANALLARINDA aranıyor.

Öldürme testi (tasarımda):
  P1: w kanalı (genlik = RS toplamı, 2θ taşır) görüntüyü DUYMALI.
  P2: v kanalı (explicit formül, 2θ terimi yok) SAĞIR kalmalı.
  P3: tamsayı-olmayan m plaseboları (aynı alias/örtüşme, aritmetik yok)
      taban çizgisini vermeli.

Konvansiyon: 70/79 makinesi + chirp kolonları; anchored (Odlyzko) pencere
2θ float64 hassasiyeti yüzünden DIŞARIDA (t ≤ 1e10).

SONUÇ (gecenin dürüst muhasebesi):
  1. İlk geçiş: tamsayı-m çizgileri gür (χ²'ler 2700'e dek, saf-cos kilitli),
     plasebo tertemiz (⟨χ²⟩=1.4), v neredeyse sağır (18 vs plasebo 8.6). AMA:
  2. ÖRNEKLEM KİMLİĞİ: gap-ortalarında θ(γ)≈πn olduğundan 2θ gözlüğü kısmen
     S-giydirilmiş DİREKT çizgiye katlanır; ve P11 tabanı asal-kuvvet
     içermiyordu → m=4,8 sesinin çoğu 2²,2³ DİREKT çizgileriymiş
     (düz kolon +0.143/+0.048 alıyor — metodolojik ders!).
  3. AYRIŞTIRMA SONRASI sahici 2θ-özel içerik: m=4: +0.0102±0.0005 (20σ),
     m=8: +0.0024±0.0004 (6σ), m=12: 2.5σ — direkt çizginin ~%5-7'si.
     KISMİ YANSIMA ADAYI. m=1 (saf 2θ/S çizgisi) kimliği açık: düz-kolon
     karşılığı DC'ye çöküyor, temiz test S-vekili tabanla kurulmalı.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A_CG = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
P11 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
M_INT = [1, 2, 3, 4, 5, 6, 8, 10, 12]
M_FAKE = [1.53, 3.7, 8.45]

def theta(t):
    return t / 2 * np.log(t / TWO_PI) - t / 2 - np.pi / 8 + 1 / (48 * t)

def spektro(gaps, amps, tmid):
    Lw = np.log(tmid / TWO_PI)
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A_CG * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    ya, yg = np.log(a_u), np.log(g_u)
    th2 = 2 * theta(tmid)
    cw = [np.ones_like(ya), g_u, g_u**2]
    cv = [np.ones_like(yg)]
    for p in P11:
        arg = tmid * np.log(p)
        cw += [np.cos(arg), np.sin(arg)]
        cv += [np.cos(arg), np.sin(arg)]
    ms = M_INT + M_FAKE
    for mq in ms:
        arg = th2 - tmid * np.log(mq)
        cw += [np.cos(arg), np.sin(arg)]
        cv += [np.cos(arg), np.sin(arg)]
    def reg(y, cols, o):
        X = np.vstack(cols).T
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        se = np.sqrt((y - X @ b).var() * np.diag(np.linalg.inv(X.T @ X)))
        out = []
        i0 = o + 2 * len(P11)
        for i, mq in enumerate(ms):
            out.append((mq, b[i0 + 2*i], se[i0 + 2*i],
                        b[i0 + 2*i + 1], se[i0 + 2*i + 1]))
        return out
    L = float(Lw.mean())
    return L, reg(ya, cw, 3), reg(yg, cv, 1)

WNDS = []
d36 = np.load(HERE / "36_T100k.npz")
edges = np.geomspace(d36["t_mid"][0], d36["t_mid"][-1] * 1.0001, 13)
for i in range(12):
    m = (d36["t_mid"] >= edges[i]) & (d36["t_mid"] < edges[i + 1])
    if m.sum() >= 3000:
        WNDS.append((d36["intervals"][m], d36["max_amps"][m], d36["t_mid"][m]))
d41 = np.load(HERE / "41_bigT_windows.npz")
for k in sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1])):
    WNDS.append((d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"]))
for f in ["55_win_1e+08.npz", "55_win_1e+09.npz", "55_win_1e+10.npz"]:
    d = np.load(HERE / f)
    WNDS.append((d["gaps"], d["amps"], d["tmid"]))

# ---- pencere pencere spektrum; kombinasyon τ′ = 1 − log m / L üzerinden
rows_w, rows_v = [], []
for gaps, amps, tmid in WNDS:
    L, rw, rv = spektro(gaps, amps, tmid)
    for (mq, c, sc, s, ss), (_, cv_, scv, sv_, ssv) in zip(rw, rv):
        tp = 1 - (np.log(mq) / L if mq >= 1 else 0)
        fake = mq not in M_INT
        rows_w.append((tp, mq, L, c, sc, s, ss, fake))
        rows_v.append((tp, mq, L, cv_, scv, sv_, ssv, fake))

W = np.array([(t, c, sc, s, ss) for t, mq, L, c, sc, s, ss, f in rows_w if not f])
Wf = np.array([(t, c, sc, s, ss) for t, mq, L, c, sc, s, ss, f in rows_w if f])
V = np.array([(t, c, sc, s, ss) for t, mq, L, c, sc, s, ss, f in rows_v if not f])
Vf = np.array([(t, c, sc, s, ss) for t, mq, L, c, sc, s, ss, f in rows_v if f])
n_w = len(WNDS)
print(f"{n_w} pencere × {len(M_INT)} tamsayı-m + {len(M_FAKE)} plasebo")

def z2(P):
    return ((P[:, 1] / P[:, 2])**2 + (P[:, 3] / P[:, 4])**2)

print("\nKANAL SESSİZLİK KARŞILAŞTIRMASI (χ² per çizgi; sıfır hipotezi dof=2):")
for isim, P, Pf in [("w (genlik)", W, Wf), ("v (boşluk)", V, Vf)]:
    print(f"  {isim}: tamsayı-m ⟨χ²⟩ = {z2(P).mean():6.2f} (n={len(P)}) | "
          f"plasebo ⟨χ²⟩ = {z2(Pf).mean():6.2f} (n={len(Pf)})")

print("\nEN GÜÇLÜ w-ÇİZGİLERİ (χ² > 12):")
sel = np.argsort(-z2(W))
for i in sel[:12]:
    t, c, sc, s, ss = W[i]
    x2 = z2(W)[i:i+1][0]
    if x2 > 12:
        mq, L = [(mq, L) for tp, mq, L, *_ , f in rows_w
                 if not f and abs(tp - t) < 1e-9][0][:2]
        print(f"  m={int(mq):>3} L={L:5.2f} τ′={t:.3f}  cos={c:+.4f}±{sc:.4f} "
              f"sin={s:+.4f}±{ss:.4f}  χ²={x2:.1f}")

# τ′ binlerinde kombine genlik (kareler toplamı, gürültü çıkarılmış)
print("\nτ′-BİN ÖZETİ (w kanalı; genlik² = cos²+sin²−2σ², plasebo aynı bin):")
print(f"{'τ′-bin':>13} {'n':>3} {'√⟨A²⟩':>8} {'plasebo':>8}")
for lo, hi in [(0.45,0.55),(0.55,0.65),(0.65,0.75),(0.75,0.85),(0.85,0.95),(0.95,1.01)]:
    m = (W[:, 0] >= lo) & (W[:, 0] < hi)
    mf = (Wf[:, 0] >= lo) & (Wf[:, 0] < hi)
    if m.sum():
        a2 = (W[m,1]**2 + W[m,3]**2 - W[m,2]**2 - W[m,4]**2).mean()
        a2f = ((Wf[mf,1]**2 + Wf[mf,3]**2 - Wf[mf,2]**2 - Wf[mf,4]**2).mean()
               if mf.sum() else np.nan)
        print(f"[{lo:.2f},{hi:.2f}) {int(m.sum()):>3} "
              f"{np.sqrt(max(a2,0)):>8.4f} {np.sqrt(max(a2f,0)) if a2f==a2f else float('nan'):>8.4f}")

# ---- ağırlık-normalize edilmiş içerik: u(m) = cos-katsayısı · √m / d(m)
# (|Z|²'de 2θ−t·log m çizgisinin toplam ağırlığı Σ_{nn'=m}(nn')^{-1/2} = d(m)/√m)
DIV = {1: 1, 2: 2, 3: 2, 4: 3, 5: 2, 6: 4, 8: 4, 10: 4, 12: 6}
print("\nNORMALİZE İÇERİK u(m) (pencereler üzerinden ağırlıklı ortalama, yalnız cos):")
print(f"{'m':>3} {'d(m)':>4} {'⟨cos⟩':>8} {'u(m)':>8} {'±':>7}")
for mq in M_INT:
    rows = [(c, sc) for tp, m_, L, c, sc, s, ss, f in rows_w if m_ == mq]
    c = np.array([r[0] for r in rows]); sc = np.array([r[1] for r in rows])
    cbar = np.sum(c / sc**2) / np.sum(1 / sc**2)
    cse = 1 / np.sqrt(np.sum(1 / sc**2))
    norm = np.sqrt(mq) / DIV[mq]
    print(f"{mq:>3} {DIV[mq]:>4} {cbar:>+8.4f} {cbar*norm:>+8.4f} {cse*norm:>7.4f}")

# ---- KİMLİK DENETİMİ: düz gözlük (cos t·log m) ve 2θ gözlüğü AYNI tabanda
print("\nKİMLİK DENETİMİ (m=4,8,12; iki gözlük yan yana):")
print(f"{'m':>3} {'düz (direkt p^k)':>18} {'2θ-özel artık':>15}")
MS_ID = [4, 8, 12]
acc = {m_: [[], []] for m_ in MS_ID}
for gaps, amps, tmid in WNDS:
    Lw = np.log(tmid / TWO_PI)
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A_CG * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    ya = np.log(a_u)
    th2 = 2 * theta(tmid)
    cols = [np.ones_like(ya), g_u, g_u**2]
    for pp in P11:
        arg = tmid * np.log(pp)
        cols += [np.cos(arg), np.sin(arg)]
    for mq in MS_ID:
        arg = tmid * np.log(mq)
        cols += [np.cos(arg), np.sin(arg)]
    for mq in MS_ID:
        arg = th2 - tmid * np.log(mq)
        cols += [np.cos(arg), np.sin(arg)]
    X = np.vstack(cols).T
    b, *_ = np.linalg.lstsq(X, ya, rcond=None)
    se = np.sqrt((ya - X @ b).var() * np.diag(np.linalg.inv(X.T @ X)))
    o = 3 + 2 * len(P11)
    for i, mq in enumerate(MS_ID):
        acc[mq][0].append((b[o + 2*i], se[o + 2*i]))
        acc[mq][1].append((b[o + 6 + 2*i], se[o + 6 + 2*i]))
for mq in MS_ID:
    d = np.array(acc[mq][0]); c = np.array(acc[mq][1])
    db = np.sum(d[:,0]/d[:,1]**2)/np.sum(1/d[:,1]**2); de = 1/np.sqrt(np.sum(1/d[:,1]**2))
    cb = np.sum(c[:,0]/c[:,1]**2)/np.sum(1/c[:,1]**2); ce = 1/np.sqrt(np.sum(1/c[:,1]**2))
    print(f"{mq:>3} {db:>+11.4f}±{de:.4f} {cb:>+9.4f}±{ce:.4f}  (2θ-özel/direkt = {cb/db:+.3f})")

# ---- grafik
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))
for ax, P, Pf, isim in [(axes[0], W, Wf, "w kanalı (genlik)"),
                        (axes[1], V, Vf, "v kanalı (boşluk)")]:
    amp = np.sqrt(np.clip(P[:,1]**2 + P[:,3]**2 - P[:,2]**2 - P[:,4]**2, 0, None))
    ampf = np.sqrt(np.clip(Pf[:,1]**2 + Pf[:,3]**2 - Pf[:,2]**2 - Pf[:,4]**2, 0, None))
    ax.plot(P[:, 0], amp, "o", ms=4, c="firebrick", alpha=0.65, label="tamsayı m")
    ax.plot(Pf[:, 0], ampf, "x", ms=6, c="gray", label="plasebo m")
    ax.axvline(0.5, color="gray", ls="--", lw=1)
    ax.set_xlabel(r"τ′ = 1 − log m / L (görüntü konumu)")
    ax.set_ylabel("chirp genliği (gürültü-düzeltilmiş)")
    ax.set_title(isim); ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(HERE / "80_goruntu_cizgileri.png", dpi=110)
print("\nGrafik: 80_goruntu_cizgileri.png")
