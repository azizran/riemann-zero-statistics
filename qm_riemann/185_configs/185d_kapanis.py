# -*- coding: utf-8 -*-
"""
185d — K3: KAPANIŞ (baskın faktör ileri-hesabı; türetilmiş (c,α); H-W2/H-W3)
============================================================================
ONKAYIT_185 dondurdu:
  r_pred^A = Σâ_g^öz/Σâ_Hk^öz          (çarpımsal-kinematik; ölçülü â girmez)
  r_pred^B = Σ|c_Hk^olc − c_Hk^öz + c_g^öz| / Σâ_Hk
             (toplamsal ortak-karışım; ölçülü â_g girmez)
Mekanizma-eşleme (K2): F1/F3 çifti baskın → birincil = B.
χ²/dof daima 184 defter-σ'sıyla (bağımsız yayılım). Fit: 1−c·τ^α,
α-taramasında kapalı-form c (serbest sabit yok; fit yalnız ÖZETLEME).
H-W3: α=3.3 sabit, en-iyi-ölçek şekil sınavı (ölçek bilgisi erişilemedi).

Çıktı: 185/K3_kapanis.json + ekran.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S184, S185 = SCR / "184", SCR / "185"

spec = importlib.util.spec_from_file_location(
    "b185", QM / "185_configs" / "185b_oz_muhasebe.py")
b185 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b185)  # __main__ değil; yalnız fonksiyonlar

ONK = json.load(open(S185 / "ONKAYIT_185.json"))
K1 = np.load(S185 / "K1_faktorler.npz", allow_pickle=True)
G = np.load(S184 / "K1_gercek.npz")
H = np.load(S184 / "K1_Hkeskin.npz")
OZg = np.load(S185 / "OZ_gercek.npz")
OZh = np.load(S185 / "OZ_Hkeskin.npz")
tau, ae, aq = G["tau"], G["aq_eff"], G["aq"]
kenar = ONK["kenar"]
maskeler = [(tau >= kenar[b]) & (tau < kenar[b + 1]) for b in range(8)]
et = [f"{kenar[b]:.2f}-{kenar[b+1]:.2f}" for b in range(8)]
tau_bar = K1["tau_bar"][:8]
r_b = K1["r"][:8]
sig_r = K1["sig_r"][:8]

print("=" * 78)
print(f"185d / K3 KAPANIŞ  [on-kayit sha {ONK['sha256'][:12]}]")
print("=" * 78)


def bant_pred(disari):
    cg_o = b185.c_olculu(G, disari)
    ch_o = b185.c_olculu(H, disari)
    cg_z = b185.c_oz(OZg, aq, disari)
    ch_z = b185.c_oz(OZh, aq, disari)
    A = np.zeros(8)
    B = np.zeros(8)
    for k, m in enumerate(maskeler):
        A[k] = np.abs(cg_z[m]).sum() / np.abs(ch_z[m]).sum()
        B[k] = np.abs(ch_o[m] - ch_z[m] + cg_z[m]).sum() / np.abs(ch_o[m]).sum()
    return A, B


A_b, B_b = bant_pred(-1)
reps = [bant_pred(i) for i in range(8)]
seA = b185.jk_se([x[0] for x in reps])
seB = b185.jk_se([x[1] for x in reps])

zA = (A_b - r_b) / sig_r
zB = (B_b - r_b) / sig_r
chiA = float(np.sum(zA ** 2) / 8)
chiB = float(np.sum(zB ** 2) / 8)

print("\nİLERİ-HESAP DEFTERİ (birincil σ = 184 defter-σ'sı):")
print(f"{'bant':>12} {'r ölçülü':>15} {'r_pred^A':>15} {'z_A':>6} "
      f"{'r_pred^B':>15} {'z_B':>6}")
for k in range(8):
    print(f"{et[k]:>12} {r_b[k]:8.4f}±{sig_r[k]:.4f} "
          f"{A_b[k]:8.4f}±{seA[k]:.4f} {zA[k]:6.1f} "
          f"{B_b[k]:8.4f}±{seB[k]:.4f} {zB[k]:6.1f}")
print(f"\nχ²/dof (dof=8, parametresiz):  A = {chiA:.2f}   B = {chiB:.2f}")

# ---------- H-W2 hükmü (donmuş eşik; birincil = B, K2 çift-baskınlığı) ----------
birincil = "B"
z_bir, chi_bir = zB, chiB
kuyruk_z = max(abs(zB[5]), abs(zB[6]), abs(zB[7]))  # τ>0.70 bantları
tum_2s = bool(np.all(np.abs(z_bir) <= 2.0))
if tum_2s and chi_bir <= 2.0:
    hw2 = "MÜHÜR"
elif (abs(zB[7]) > 2.0 or kuyruk_z > 2.0) or chi_bir > 2.0:
    hw2 = "ÖLDÜ"
else:
    hw2 = "KISMİ"
print(f"H-W2 (birincil={birincil}; mekanizma-eşleme K2: F1/F3 çifti baskın): {hw2}"
      f"  [kuyruk maks|z|={kuyruk_z:.1f}, χ²/dof={chi_bir:.2f}]")

# ---------- güç-yasası fiti: 1 − c·τ^α (α-tarama + kapalı-form c) ----------
def guc_fit(y, sig):
    wgt = 1.0 / sig ** 2
    alfalar = np.arange(0.05, 6.0, 0.001)
    best = None
    for a in alfalar:
        ta = tau_bar ** a
        c = np.sum(wgt * (1.0 - y) * ta) / np.sum(wgt * ta * ta)
        chi = np.sum(wgt * (y - (1.0 - c * ta)) ** 2)
        if best is None or chi < best[2]:
            best = (c, a, chi)
    c, a, chi = best
    return c, a, chi / (8 - 2)


c_olc, a_olc, x_olc = guc_fit(r_b, sig_r)
c_A, a_A, x_A = guc_fit(A_b, sig_r)
c_B, a_B, x_B = guc_fit(B_b, sig_r)
print("\nGÜÇ-YASASI FİTLERİ  1 − c·τ^α  (fit-χ²/dof, dof=6):")
print(f"  ölçülü (kontrol) : c = {c_olc:+.4f}  α = {a_olc:.3f}  χ²/dof = {x_olc:.2f}"
      f"   [184: c=0.149 α=1.30 χ²/dof=0.34]")
print(f"  r_pred^A         : c = {c_A:+.4f}  α = {a_A:.3f}  χ²/dof = {x_A:.2f}")
print(f"  r_pred^B         : c = {c_B:+.4f}  α = {a_B:.3f}  χ²/dof = {x_B:.2f}")

# ---------- H-W3: script-100 şekil sınavı (α=3.3 sabit; ölçek erişilemedi) ----------
ta = tau_bar ** 3.3
wgt = 1.0 / sig_r ** 2
c3 = float(np.sum(wgt * (1.0 - r_b) * ta) / np.sum(wgt * ta * ta))
chi3 = float(np.sum(wgt * (r_b - (1.0 - c3 * ta)) ** 2) / (8 - 1))
print(f"\nH-W3 — script-100 şekli (1 − c₃·τ^3.3; ölçek bilgisi ERİŞİLEMEDİ,"
      f" yalnız en-iyi-ölçek):")
print(f"  en-iyi c₃ = {c3:.4f}  →  χ²/dof = {chi3:.2f}  (dof=7)"
      f"   [α=3.3 ≠ türetilen α; kıyas raporda]")

json.dump({
    "sha_onkayit": ONK["sha256"],
    "birincil": birincil,
    "r_olc": r_b.tolist(), "sig_r": sig_r.tolist(), "tau_bar": tau_bar.tolist(),
    "r_pred_A": A_b.tolist(), "se_A": seA.tolist(), "z_A": zA.tolist(),
    "chi2dof_A": chiA,
    "r_pred_B": B_b.tolist(), "se_B": seB.tolist(), "z_B": zB.tolist(),
    "chi2dof_B": chiB,
    "H-W2": hw2,
    "fit": {"olculu": [c_olc, a_olc, x_olc], "A": [c_A, a_A, x_A],
            "B": [c_B, a_B, x_B], "hedef_184": [0.149, 1.30, 0.34]},
    "H-W3": {"c3_en_iyi": c3, "chi2dof": chi3, "alfa_sabit": 3.3,
             "olcek": "erişilemedi (okuma listesinde yok)"},
}, open(S185 / "K3_kapanis.json", "w"), indent=1, ensure_ascii=False)
print("\n-> 185/K3_kapanis.json  BİTTİ")
