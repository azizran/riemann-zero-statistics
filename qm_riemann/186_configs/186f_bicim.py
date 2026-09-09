# -*- coding: utf-8 -*-
"""
186f — G3: BİÇİM DEFTERİ (M(τ)'ye üç aday biçim; türetim değil, kısıt)
======================================================================
ONKAYIT_186 dondurdu: (i) M=a+b·τ, (ii) M=1−k·τ^β (β-tarama, kapalı-form k),
(iii) M=a+b·V(τ) (V=ölçülü v_g/v_Hk). Ağırlık 1/σ² (birincil defter-seM);
dof=6; parametre jk-se'leri 8 loo-replikadan.
Çıktı: 186/G3_bicim.json + ekran.
"""
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S184, S185, S186 = SCR / "184", SCR / "185", SCR / "186"
NJACK = 8

spec = importlib.util.spec_from_file_location(
    "b185", QM / "185_configs" / "185b_oz_muhasebe.py")
b185 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b185)

ONK = json.load(open(S186 / "ONKAYIT_186.json"))
sha = hashlib.sha256(
    (QM / "186_configs" / "186a_onkayit.py").read_bytes()).hexdigest()
assert sha == ONK["sha256"], "ON-KAYIT SHA UYUMSUZ"

G = np.load(S184 / "K1_gercek.npz")
H = np.load(S184 / "K1_Hkeskin.npz")
OZg = np.load(S185 / "OZ_gercek.npz")
OZh = np.load(S185 / "OZ_Hkeskin.npz")
K1f = np.load(S185 / "K1_faktorler.npz", allow_pickle=True)
G2S = json.load(open(S186 / "G2_sonuc.json"))

tau = np.asarray(G["tau"], float)
aq = np.asarray(G["aq"], float)
ae = np.asarray(G["aq_eff"], float)
kenar = ONK["kenar"]
maskeler = [(tau >= kenar[b]) & (tau < kenar[b + 1]) for b in range(8)]
tau_bar = K1f["tau_bar"][:8]
V_tam = np.array(G2S["V_oran"])

# ölçülü M (tam + loo) ve birincil σ (defter)
se185 = K1f["se_tablo"]
wg, wh = K1f["wg"][:8], K1f["wh"][:8]
wog, woh = K1f["wog"][:8], K1f["woh"][:8]
mg, mh = wg - wog, wh - woh
M_b = mg / mh
se_mg = np.sqrt(se185[:8, 7] ** 2 + se185[:8, 9] ** 2)
se_mh = np.sqrt(se185[:8, 8] ** 2 + se185[:8, 10] ** 2)
sigM = M_b * np.sqrt((se_mg / mg) ** 2 + (se_mh / mh) ** 2)
wgt = 1.0 / sigM ** 2


def M_loo(disari):
    ag = np.abs(b185.c_olculu(G, disari))
    ah = np.abs(b185.c_olculu(H, disari))
    og = np.abs(b185.c_oz(OZg, aq, disari))
    oh = np.abs(b185.c_oz(OZh, aq, disari))
    out = np.zeros(8)
    for k, m in enumerate(maskeler):
        sae = ae[m].sum()
        out[k] = ((ag[m].sum() - og[m].sum()) / sae) / \
                 ((ah[m].sum() - oh[m].sum()) / sae)
    return out


def V_loo(disari):
    def vb(D):
        N, nb = int(D["N"]), D["nb"]
        re = D["Gre"].sum(0) - (D["Gre"][disari] if disari >= 0 else 0)
        im = D["Gim"].sum(0) - (D["Gim"][disari] if disari >= 0 else 0)
        n = N - (nb[disari] if disari >= 0 else 0)
        aG = np.abs(re + 1j * im) / n
        return np.array([aG[m].sum() / (float(D["gbar"]) * ae[m].sum())
                         for m in maskeler])
    return vb(OZg) / vb(OZh)


def fit_dogrusal(y, x):
    """ağırlıklı LS y = a + b·x → (a, b, χ²)."""
    W, X, Y = wgt, x, y
    Sw, Sx, Sy = W.sum(), (W * X).sum(), (W * Y).sum()
    Sxx, Sxy = (W * X * X).sum(), (W * X * Y).sum()
    d = Sw * Sxx - Sx ** 2
    a = (Sxx * Sy - Sx * Sxy) / d
    b = (Sw * Sxy - Sx * Sy) / d
    chi = float((W * (Y - a - b * X) ** 2).sum())
    return float(a), float(b), chi


def fit_guc(y):
    """M = 1 − k·τ^β; β-tarama [0.05,6.0) adım 0.001, kapalı-form k."""
    best = None
    for be in np.arange(0.05, 6.0, 0.001):
        tb = tau_bar ** be
        k = np.sum(wgt * (1.0 - y) * tb) / np.sum(wgt * tb * tb)
        chi = float(np.sum(wgt * (y - (1.0 - k * tb)) ** 2))
        if best is None or chi < best[2]:
            best = (float(k), float(be), chi)
    return best


print("=" * 78)
print(f"186f / G3 BİÇİM DEFTERİ  [on-kayit sha {ONK['sha256'][:12]}]")
print("=" * 78)

# tam-örneklem fitler
a1, b1, x1 = fit_dogrusal(M_b, tau_bar)
k2, be2, x2 = fit_guc(M_b)
a3, b3, x3 = fit_dogrusal(M_b, V_tam)

# jackknife (M ve V loo; σ-ağırlıkları donmuş)
p1, p2, p3 = [], [], []
for i in range(NJACK):
    Mi = M_loo(i)
    Vi = V_loo(i)
    p1.append(fit_dogrusal(Mi, tau_bar)[:2])
    p2.append(fit_guc(Mi)[:2])
    p3.append(fit_dogrusal(Mi, Vi)[:2])
se1 = b185.jk_se(p1)
se2 = b185.jk_se(p2)
se3 = b185.jk_se(p3)

print("\nBİÇİM DEFTERİ (dof=6; ±jk se):")
print(f"  doğrusal M=a+b·τ    : a = {a1:.4f}±{se1[0]:.4f}  "
      f"b = {b1:+.4f}±{se1[1]:.4f}  χ²/dof = {x1/6:.2f}")
print(f"  güç     M=1−k·τ^β   : k = {k2:.4f}±{se2[0]:.4f}  "
      f"β = {be2:.3f}±{se2[1]:.3f}  χ²/dof = {x2/6:.2f}")
print(f"  v-afin  M=a+b·V(τ)  : a = {a3:+.4f}±{se3[0]:.4f}  "
      f"b = {b3:+.4f}±{se3[1]:.4f}  χ²/dof = {x3/6:.2f}")
print(f"\n  (KALEM ön-bakışı: 1−M ~ τ^0.56-ish '√τ kokusu' — ölçülen β = "
      f"{be2:.3f})")
print(f"  (184 zarf kıyası: r = 1−0.149·τ^1.30 — M'nin (k,β)'sı = "
      f"({k2:.3f}, {be2:.3f}): zarftan hem ölçek hem üs olarak AYRI nesne)")

json.dump({
    "sha_onkayit": ONK["sha256"], "tau_bar": tau_bar.tolist(),
    "M_olc": M_b.tolist(), "sigM_defter": sigM.tolist(),
    "V_oran": V_tam.tolist(),
    "dogrusal": {"a": a1, "b": b1, "se": se1.tolist(), "chi2dof": x1 / 6},
    "guc": {"k": k2, "beta": be2, "se": se2.tolist(), "chi2dof": x2 / 6},
    "v_afin": {"a": a3, "b": b3, "se": se3.tolist(), "chi2dof": x3 / 6},
}, open(S186 / "G3_bicim.json", "w"), indent=1, ensure_ascii=False)
print(f"\n-> {S186/'G3_bicim.json'}  BİTTİ")
