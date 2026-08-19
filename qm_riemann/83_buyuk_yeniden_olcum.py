"""
83 — BÜYÜK YENİDEN ÖLÇÜM: TAM-TABAN KAMPANYASI (20 Ağustos 2026)
==========================================================================
82'nin yer-gerçeği kanıtı üzerine: üçlemenin bütün manşet sayıları tam
tabanla (tüm p^k, pencere başına τ ≤ 0.55 veya q ≤ 720 sınırıyla)
yeniden ölçülüyor. Kalan sistematik: 360→720 merdiven adımı b'yi ~%1
oynatıyor → taban-kesimi belirsizliği olarak not edilir.

  T1  Korunum yasası (44 çift) + öldürme testleri + kuvvet-u'ları + w₀ izi
  T2  132-nokta w(τ) eğrisi + kesiş τ₀ + faz (kayıpsızlık) kontrolü
  T3  Termal Bragg yeniden fiti (ölçülen c_jit, YENİ B, yeni eğri)
  T5  Soyulmuş çekirdek r* (12 pencere, eski vs tam taban)
  (s_eff = 4.11 MUAF: u-serisinin karakteristik fonksiyonu, regresyonsuz)

Çıktılar: ESKİ→YENİ tablosu, 83_tam_taban_egri.npz, 2 karşılaştırma figürü.
"""

import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A_CG = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
T0_ODL = 267653395647
P4 = [2, 3, 5, 7]
P2 = [11, 13]
P11 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
Q_CAP = 720
mp.mp.dps = 30

def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q)
            q *= p
    return sorted(set(out))

def tam_taban(L, tau_max=0.45):
    return pk_list(min(np.exp(tau_max * L), Q_CAP))

def cols_of(tmid, qs, ph, sl, g_u=None):
    cols = [np.ones(sl.stop - sl.start)]
    if g_u is not None:
        cols += [g_u[sl], g_u[sl]**2]
    for q in qs:
        arg = ph[q] + tmid[sl] * np.log(q)
        cols += [np.cos(arg), np.sin(arg)]
    return np.vstack(cols).T

def chunked_reg(y, tmid, qs, anchored, g_u=None, chunk=40000, resid_with=None):
    C = (3 if g_u is not None else 1) + 2 * len(qs)
    XtX = np.zeros((C, C)); Xty = np.zeros(C)
    ph = {q: (float(mp.fmod(T0_ODL * mp.log(q), 2 * mp.pi)) if anchored else 0.0)
          for q in qs}
    n = len(y); ss = 0.0
    for s0 in range(0, n, chunk):
        sl = slice(s0, min(s0 + chunk, n))
        Xc = cols_of(tmid, qs, ph, sl, g_u)
        XtX += Xc.T @ Xc; Xty += Xc.T @ y[sl]
    b = np.linalg.solve(XtX, Xty)
    res = np.empty(n) if resid_with is not None or True else None
    for s0 in range(0, n, chunk):
        sl = slice(s0, min(s0 + chunk, n))
        Xc = cols_of(tmid, qs, ph, sl, g_u)
        r = y[sl] - Xc @ b
        res[sl] = r
        ss += (r**2).sum()
    se = np.sqrt((ss / n) * np.diag(np.linalg.inv(XtX)))
    return b, se, res

def unfold(gaps, amps, tmid, anchored=False):
    Lw = np.log(((T0_ODL + tmid) if anchored else tmid) / TWO_PI)
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A_CG * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    return np.log(a_u), np.log(g_u), g_u, float(Lw.mean())

def kanal_tam(gaps, amps, tmid, targets, anchored, tau_max=0.45, powers_too=True):
    ya, yg, g_u, L = unfold(gaps, amps, tmid, anchored)
    qs = sorted(set(targets) | set(tam_taban(L, tau_max)))
    bw, sew, _ = chunked_reg(ya, tmid, qs, anchored, g_u=g_u)
    bv, sev, _ = chunked_reg(yg, tmid, qs, anchored, g_u=None)
    out = {}
    for p in targets:
        i = qs.index(p)
        iw, iv = 3 + 2 * i, 1 + 2 * i
        out[p] = (np.log(p) / L, bw[iw] / p**-0.5, sew[iw] / p**-0.5,
                  bw[iw + 1] / p**-0.5,
                  float(np.hypot(bv[iv], bv[iv + 1]) / p**-0.5))
    guc = {}
    if powers_too:
        for q, pp, kk in [(4, 2, 2), (8, 2, 3), (9, 3, 2), (25, 5, 2), (27, 3, 3)]:
            if q in qs:
                i = qs.index(q)
                w_efa = pp**(-kk/2) / kk
                guc[q] = (bw[3 + 2*i] / w_efa, sew[3 + 2*i] / w_efa)
    return out, guc, L

# ---- pencereler
LAW = []
d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
for k in K41[:6]:
    LAW.append((d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"], False))
for f in ["55_win_1e+08.npz", "55_win_1e+09.npz", "55_win_1e+10.npz", "55_win_1e+11.npz"]:
    d = np.load(HERE / f)
    LAW.append((d["gaps"], d["amps"], d["tmid"], False))
d53 = np.load(HERE / "53_odlyzko_amps.npz")
LAW.append((d53["gaps"], d53["max_amps"], d53["t_mid"], True))

EXT = []
d36 = np.load(HERE / "36_T100k.npz")
edges = np.geomspace(d36["t_mid"][0], d36["t_mid"][-1] * 1.0001, 13)
for i in range(12):
    m = (d36["t_mid"] >= edges[i]) & (d36["t_mid"] < edges[i + 1])
    if m.sum() >= 500:
        EXT.append((d36["intervals"][m], d36["max_amps"][m], d36["t_mid"][m], False))
for k in K41:
    EXT.append((d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"], False))

# ============ T1: KORUNUM YASASI + ÖLDÜRME TESTLERİ + GÜÇLER ============
print("T1 — KORUNUM YASASI (tam taban):")
tr, oos, widx = [], [], []
u_guc = {4: [], 8: [], 9: [], 25: [], 27: []}
w2_iz = []
for wi, (gaps, amps, tmid, anch) in enumerate(LAW):
    out, guc, L = kanal_tam(gaps, amps, tmid, P4 + P2, anch)
    for p in P4 + P2:
        tau, w_, s_, wsin, v_ = out[p]
        if w_ > 0:
            (tr if p in P4 else oos).append((v_, np.sqrt(w_), w_, tau))
            if p in P4:
                widx.append(wi)
    for q, (u, su) in guc.items():
        u_guc[q].append((u, su, L))
    w2_iz.append((L, out[2][0], out[2][1], out[2][2]))
tr = np.array(tr); oos = np.array(oos); widx = np.array(widx)
vv, sq, wv, _ = tr.T
X = np.vstack([np.ones_like(vv), vv]).T
c, *_ = np.linalg.lstsq(X, sq, rcond=None)
rms = np.sqrt(np.mean((sq - X @ c)**2))
b1 = np.sum((1 - sq) * vv) / np.sum(vv * vv)
rms1 = np.sqrt(np.mean((sq - 1 + b1 * vv)**2))
rng = np.random.default_rng(83)
uw = np.unique(widx); slopes = []
for _ in range(400):
    sel = rng.choice(uw, len(uw), replace=True)
    ii = np.concatenate([np.where(widx == s)[0] for s in sel])
    cc, *_ = np.linalg.lstsq(X[ii], sq[ii], rcond=None)
    slopes.append(-cc[1])
print(f"  SERBEST: √w = {c[0]:.4f} − {-c[1]:.4f}·v  RMS {rms:.4f}  "
      f"(blok-bootstrap eğim ± {np.std(slopes):.4f})")
print(f"  a≡1   : b = {b1:.4f}  RMS {rms1:.4f}")
vo, so, wo, to = oos.T
pred = 1 - b1 * vo
pred_f = c[0] + c[1] * vo
print(f"  OUT-OF-SAMPLE p=11,13 ({len(oos)} çift): a≡1 doğrusuna RMS = "
      f"{np.sqrt(np.mean((so - pred)**2)):.4f} (ofset {np.mean(so - pred):+.4f}) | "
      f"serbest doğruya RMS = {np.sqrt(np.mean((so - pred_f)**2)):.4f} "
      f"(ofset {np.mean(so - pred_f):+.4f})")
cq = np.sum((1 - wv) * vv**2) / np.sum(vv**4)
c2, *_ = np.linalg.lstsq(X, wv, rcond=None)
r_sqrt = np.sqrt(np.mean(((X @ c)**2 - wv)**2))
r_lin = np.sqrt(np.mean((X @ c2 - wv)**2))
r_uni = np.sqrt(np.mean(((1 - cq * vv**2) - wv)**2))
print(f"  FORM YARIŞI (w-ölçeği): √w-lineer {r_sqrt:.4f} | w-lineer {r_lin:.4f} | "
      f"şiddet-üni {r_uni:.4f} (×{r_uni/r_sqrt:.1f})")
print("  KUVVET-u'LARI (u = katsayı·k·p^{k/2}; sum-rule tahmini |u|≈w(τ_q)):")
for q in [4, 8, 9, 25, 27]:
    a = np.array([(u, su) for u, su, L in u_guc[q]])
    if len(a):
        um = np.sum(a[:,0]/a[:,1]**2)/np.sum(1/a[:,1]**2)
        ue = 1/np.sqrt(np.sum(1/a[:,1]**2))
        print(f"    u({q:>2}) = {um:+.3f} ± {ue:.3f}")
print("  w(p=2) İZİ (w₀ sorusu):")
for L, tau, w_, s_ in w2_iz:
    print(f"    L={L:5.2f} τ={tau:.3f}: w = {w_:.4f} ± {s_:.4f}")

# ============ T2: 132-NOKTA EĞRİ + KESİŞ + FAZLAR ============
print("\nT2 — TAM-TABAN w(τ) EĞRİSİ (P11 × 12 pencere):")
pts = []
for gaps, amps, tmid, anch in EXT:
    out, _, L = kanal_tam(gaps, amps, tmid, P11, anch, tau_max=0.55, powers_too=False)
    for p in P11:
        tau, w_, s_, wsin, v_ = out[p]
        pts.append((tau, w_, s_, wsin, v_, L))
pts = np.array(sorted(pts, key=lambda r: r[0]))
np.savez(HERE / "83_tam_taban_egri.npz", tau=pts[:,0], w=pts[:,1], sw=pts[:,2],
         wsin=pts[:,3], v=pts[:,4], L=pts[:,5])
m = (pts[:, 0] >= 0.30) & (pts[:, 0] <= 0.58)
t_, w_, s_ = pts[m, 0], pts[m, 1], pts[m, 2]
Xc_ = np.vstack([np.ones_like(t_), t_]).T / s_[:, None]
cc, *_ = np.linalg.lstsq(Xc_, w_ / s_, rcond=None)
tau0 = -cc[0] / cc[1]
boots = []
Ls = pts[m, 5]
for _ in range(400):
    sel = rng.choice(np.unique(Ls), len(np.unique(Ls)), replace=True)
    ii = np.concatenate([np.where(Ls == s)[0] for s in sel])
    ccb, *_ = np.linalg.lstsq(Xc_[ii], (w_ / s_)[ii], rcond=None)
    boots.append(-ccb[0] / ccb[1])
print(f"  KESİŞ: τ₀ = {tau0:.3f} ± {np.std(boots):.3f} (eski 0.40 ± 0.02)")
for lo, hi in [(0.0, 0.2), (0.2, 0.35), (0.35, 0.45), (0.45, 0.62)]:
    mm = (pts[:, 0] >= lo) & (pts[:, 0] < hi)
    if mm.sum():
        sinr = np.abs(pts[mm, 3]).mean() / max(np.abs(pts[mm, 1]).mean(), 1e-9)
        print(f"  τ∈[{lo:.2f},{hi:.2f}): ⟨w⟩ = {pts[mm,1].mean():+.4f}  "
              f"⟨|sin|⟩/⟨|cos|⟩ = {sinr:.3f}")

# ============ T3: TERMAL BRAGG YENİDEN FİT ============
def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8
cs = []
for gaps, amps, tmid, anch in EXT:
    t0 = tmid[0] - gaps[0] / 2
    kk = np.arange(len(tmid)) + 0.5
    t = t0 + kk * TWO_PI / np.log(t0 / TWO_PI)
    for _ in range(6):
        fdel = rvm_N(t) - rvm_N(t0) - kk
        t = t - fdel / (np.log(t / TWO_PI) / TWO_PI)
    L = float(np.log(tmid / TWO_PI).mean())
    cs.append((L * (tmid - t).std()) ** 2 / 2)
c_jit = float(np.mean(cs))
tauA, wA, swA = pts[:, 0], pts[:, 1], pts[:, 2]
mB = tauA > 0.5
if mB.sum() == 0:
    raise SystemExit("kat-ötesi nokta yok — pencere kesimini kontrol et")
B_new = float(np.sum(wA[mB] / swA[mB]**2) / np.sum(1 / swA[mB]**2))
from scipy.optimize import least_squares
fT = least_squares(lambda p: ((p[0]*(1-2*tauA)*np.exp(-c_jit*tauA**2)*(tauA<0.5)
                               + B_new) - wA) / swA, [1.05])
tt = np.linspace(0.01, 0.499, 4000)
fmod = fT.x[0]*(1-2*tt)*np.exp(-c_jit*tt**2) + B_new
t_model = tt[np.argmax(fmod < 0)] if fmod.min() < 0 else np.nan
print(f"\nT3 — TERMAL BRAGG (tam eğri): c_jit = {c_jit:.2f} (ölçüm, değişmez), "
      f"B_yeni = {B_new:+.4f}, A = {fT.x[0]:.3f}, model kesişi = {t_model:.3f}")

# ============ T5: SOYULMUŞ ÇEKİRDEK ============
print("\nT5 — SOYULMUŞ ÇEKİRDEK r* (12 pencere):")
r_raw, r_old, r_new = [], [], []
for gaps, amps, tmid, anch in EXT:
    ya, yg, g_u, L = unfold(gaps, amps, tmid, anch)
    r_raw.append(np.corrcoef(ya, yg)[0, 1])
    for qs, dst in [(list(P11), r_old),
                    (sorted(set(P11) | set(tam_taban(L, 0.55))), r_new)]:
        _, _, ra = chunked_reg(ya, tmid, qs, anch, g_u=None)
        _, _, rg = chunked_reg(yg, tmid, qs, anch, g_u=None)
        dst.append(np.corrcoef(ra, rg)[0, 1])
print(f"  ham ⟨r⟩ = {np.mean(r_raw):.4f} | P11-soyulmuş ⟨r*⟩ = {np.mean(r_old):.4f} "
      f"| TAM-soyulmuş ⟨r*⟩ = {np.mean(r_new):.4f}")

# ============ FİGÜRLER ============
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))
ax = axes[0]
ax.errorbar(pts[:, 0], pts[:, 1], yerr=pts[:, 2], fmt="o", ms=3.5, c="firebrick",
            alpha=0.6, label="tam taban (yeni)")
d_old = None
try:
    ax.plot(tt, fT.x[0]*(1-2*tt)*np.exp(-c_jit*tt**2) + B_new, "k-", lw=1.3,
            label=f"termal model (B={B_new:+.3f}, kesiş {t_model:.2f})")
except Exception:
    pass
ax.axvline(0.5, color="gray", ls="--", lw=1)
ax.axhline(0, color="gray", lw=0.6)
ax.axvline(tau0, color="steelblue", ls=":", lw=1.2, label=f"ölçülen kesiş {tau0:.3f}")
ax.set_xlabel("τ"); ax.set_ylabel("w"); ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax.set_title("Tam-taban w(τ): eğri, kesiş, plato")
ax = axes[1]
ax.plot(tr[:, 0], tr[:, 1], "o", ms=4, c="steelblue", alpha=0.6, label="P4 (44 çift)")
ax.plot(oos[:, 0], oos[:, 1], "s", ms=6, c="firebrick", zorder=5, label="p=11,13 (oos)")
xx = np.linspace(0, max(vv.max(), vo.max()) * 1.05, 50)
ax.plot(xx, c[0] + c[1] * xx, "k-", lw=1.2, label=f"√w = {c[0]:.3f} − {-c[1]:.3f}v")
ax.plot(xx, 1 - 1.0683 * xx, "k--", lw=1, alpha=0.5, label="eski: 1 − 1.068v")
ax.set_xlabel("v"); ax.set_ylabel("√w"); ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax.set_title("Korunum yasası: tam taban")
plt.tight_layout()
plt.savefig(HERE / "83_kampanya.png", dpi=110)
print("\nFigür: 83_kampanya.png | Eğri: 83_tam_taban_egri.npz")
