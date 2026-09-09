# -*- coding: utf-8 -*-
"""
185e — K4: HA4 SINAVI (örneklem-dışı; sıfır yeni ayar)
======================================================
ONKAYIT_185 dondurdu:
  env(τ) = 0.5·erfc((τ−0.68)/0.125)          (164/184 defterinden)
  c_HA4^öz,erfc = env(τ)·c_HA4^öz            (nominal formül, HA4'ün g,m'si)
  BİRİNCİL (K4-b): c_HA4^pred = (c_Hk^olc − c_Hk^öz) + c_HA4^öz,erfc
    w_pred_b = Σ|c^pred|/Σae  vs  ölçülü w_HA4 (K1_HA4.npz)
    σ² = se_olc² + se_pred² (bağımsız yayılım); bant-χ²/dof, dof=8
  YAN (K4-a): r_pred = Σ|c_HA4^öz,erfc|/Σ|c_Hk^öz|  vs  r_olc = Σâ_HA4/Σâ_Hk

Çıktı: 185/OZ_HA4.npz + 185/K4_HA4.json + ekran.
"""
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S184, S185 = SCR / "184", SCR / "185"

spec = importlib.util.spec_from_file_location(
    "b185", QM / "185_configs" / "185b_oz_muhasebe.py")
b185 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b185)

ONK = json.load(open(S185 / "ONKAYIT_185.json"))
G = np.load(S184 / "K1_gercek.npz")
H = np.load(S184 / "K1_Hkeskin.npz")
A4 = np.load(S184 / "K1_HA4.npz")
tau, ae, aq = G["tau"], G["aq_eff"], G["aq"]
w = np.asarray(G["w"], float)
kenar = ONK["kenar"]
maskeler = [(tau >= kenar[b]) & (tau < kenar[b + 1]) for b in range(8)]
et = [f"{kenar[b]:.2f}-{kenar[b+1]:.2f}" for b in range(8)]

print("=" * 78)
print(f"185e / K4 HA4 SINAVI  [on-kayit sha {ONK['sha256'][:12]}]")
print("=" * 78, flush=True)

if not (S185 / "OZ_HA4.npz").exists():
    b185.oz_hesapla("HA4", w)
OZ4 = np.load(S185 / "OZ_HA4.npz")
OZh = np.load(S185 / "OZ_Hkeskin.npz")

try:
    from scipy.special import erfc as _erfc
    env = 0.5 * _erfc((tau - 0.68) / 0.125)
except ImportError:
    env = 0.5 * np.array([math.erfc(x) for x in (tau - 0.68) / 0.125])


def defter(disari):
    ch_o = b185.c_olculu(H, disari)
    c4_o = b185.c_olculu(A4, disari)
    ch_z = b185.c_oz(OZh, aq, disari)
    c4_z = env * b185.c_oz(OZ4, aq, disari)
    w_olc = np.zeros(8)
    w_pred = np.zeros(8)
    r_olc = np.zeros(8)
    r_pred = np.zeros(8)
    for k, m in enumerate(maskeler):
        w_olc[k] = np.abs(c4_o[m]).sum() / ae[m].sum()
        w_pred[k] = np.abs(ch_o[m] - ch_z[m] + c4_z[m]).sum() / ae[m].sum()
        r_olc[k] = np.abs(c4_o[m]).sum() / np.abs(ch_o[m]).sum()
        r_pred[k] = np.abs(c4_z[m]).sum() / np.abs(ch_z[m]).sum()
    return w_olc, w_pred, r_olc, r_pred


tam = defter(-1)
reps = [defter(i) for i in range(8)]
se = [b185.jk_se([x[j] for x in reps]) for j in range(4)]
w_olc, w_pred, r_olc, r_pred = tam
sig_b = np.sqrt(se[0] ** 2 + se[1] ** 2)
z_b = (w_pred - w_olc) / sig_b
chi_b = float(np.sum(z_b ** 2) / 8)

print("\nK4-b BİRİNCİL — TOPLAMSAL ÖZ-DEĞİŞİM (w_HA4 defteri):")
print(f"{'bant':>12} {'w_HA4 ölçülü':>15} {'w_HA4 pred':>15} {'z':>6}")
for k in range(8):
    print(f"{et[k]:>12} {w_olc[k]:8.4f}±{se[0][k]:.4f} "
          f"{w_pred[k]:8.4f}±{se[1][k]:.4f} {z_b[k]:6.1f}")
print(f"bant-χ²/dof = {chi_b:.2f}  (dof=8; bonus-sınav eşiği 2)")

zr = (r_pred - r_olc) / np.sqrt(se[2] ** 2 + se[3] ** 2)
chi_a = float(np.sum(zr ** 2) / 8)
print("\nK4-a YAN — ÇARPIMSAL (teşhis; r_HA4 = â_HA4/â_Hk):")
print(f"{'bant':>12} {'r ölçülü':>15} {'r_pred':>15} {'z':>8}")
for k in range(8):
    print(f"{et[k]:>12} {r_olc[k]:8.4f}±{se[2][k]:.4f} "
          f"{r_pred[k]:8.4f}±{se[3][k]:.4f} {zr[k]:8.1f}")
print(f"bant-χ²/dof = {chi_a:.2f}")

json.dump({
    "sha_onkayit": ONK["sha256"],
    "w_olc": w_olc.tolist(), "se_olc": se[0].tolist(),
    "w_pred": w_pred.tolist(), "se_pred": se[1].tolist(),
    "z": z_b.tolist(), "chi2dof_b": chi_b,
    "r_olc": r_olc.tolist(), "r_pred_carpimsal": r_pred.tolist(),
    "chi2dof_a": chi_a,
}, open(S185 / "K4_HA4.json", "w"), indent=1, ensure_ascii=False)
print("\n-> 185/K4_HA4.json  BİTTİ", flush=True)
