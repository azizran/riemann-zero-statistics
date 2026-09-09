# -*- coding: utf-8 -*-
"""
185c — K2: ANATOMİ (faktörlerin τ-profilleri, baskın terim, Taylor-teşhisi, H-W0)
=================================================================================
Yalnız-yorum Taylor teşhisi (hesap DAİMA kesin; ONKAYIT_185):
  T1 = (πτ)·cot(πτ)·Re[Ĝ_q]/ḡ      (uyumlu geri-besleme)
  T2 = −(πτ)²·Δσ_ε²/2               (titreşim; Δ = gerçek − Hkeskin)
Baskın faktör: kuyruk bandında log-pay pay_i = log F_i / log r.
H-W0: pencere içi L kayması sınırı (ilk/son sıfırdan L_min/L_maks).

Çıktı: 185/K2_anatomi.json + 185/K2_G_profil.npz + ekran.
"""
import json
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S155, S184, S185 = SCR / "155", SCR / "184", SCR / "185"
TWO_PI = 2 * np.pi

ONK = json.load(open(S185 / "ONKAYIT_185.json"))
K1 = np.load(S185 / "K1_faktorler.npz", allow_pickle=True)
G = np.load(S184 / "K1_gercek.npz")
OZg = np.load(S185 / "OZ_gercek.npz")
OZh = np.load(S185 / "OZ_Hkeskin.npz")

tau, ae = G["tau"], G["aq_eff"]
kenar = ONK["kenar"]
maskeler = [(tau >= kenar[b]) & (tau < kenar[b + 1]) for b in range(8)]
maskeler.append(tau > ONK["kuyruk_tau"])
et = [f"{kenar[b]:.2f}-{kenar[b+1]:.2f}" for b in range(8)] + ["KUYRUK>0.70"]
tau_bar = K1["tau_bar"]

r_b, F1_b, F2_b, F3i_b = K1["r"], K1["F1"], K1["F2"], K1["F3i"]

print("=" * 78)
print(f"185c / K2 ANATOMİ  [on-kayit sha {ONK['sha256'][:12]}]")
print("=" * 78)

# ---------- baskın faktör: log-paylar ----------
print("\nLOG-PAYLAR  pay_i = log F_i / log r  (toplam-oranı lehçesi; Σ=1 kesin):")
print(f"{'bant':>12} {'log r':>8} {'pay(ANOM)':>10} {'pay(KİN)':>10} {'pay(KOMŞU⁻¹)':>12}")
paylar = np.zeros((9, 3))
for k in range(9):
    lr = np.log(r_b[k])
    p = np.array([np.log(F1_b[k]), np.log(F2_b[k]), np.log(F3i_b[k])]) / lr
    paylar[k] = p
    print(f"{et[k]:>12} {lr:8.4f} {p[0]:10.3f} {p[1]:10.3f} {p[2]:12.3f}")
kuyruk_pay = paylar[8]
bask_i = int(np.argmax(np.abs(kuyruk_pay)))
bask_ad = ["F_ANOMALİ", "F_KİNEMATİK", "F_KOMŞU⁻¹"][bask_i]
# F1/F3 çifti mi tek faktör mü: çiftin ortak payı
cift_pay = kuyruk_pay[0] + kuyruk_pay[2]
print(f"\nkuyruk baskın tek faktör: {bask_ad} (pay {kuyruk_pay[bask_i]:.3f});"
      f"  F_ANOM×F_KOMŞU⁻¹ çift payı = {cift_pay:.3f}")

# ---------- Ĝ(τ) profili ve Taylor-teşhisi ----------
def Ghat(OZ):
    N = int(OZ["N"])
    return (OZ["Gre"].sum(0) + 1j * OZ["Gim"].sum(0)) / N

def Ghat2(OZ):
    N = int(OZ["N"])
    return (OZ["G2re"].sum(0) + 1j * OZ["G2im"].sum(0)) / N

Gg, Gh = Ghat(OZg), Ghat(OZh)
G2g, G2h = Ghat2(OZg), Ghat2(OZh)
gbar_g, gbar_h = float(OZg["gbar"]), float(OZh["gbar"])
sig2_g = float(OZg["sdg2"].sum()) / int(OZg["N"]) / gbar_g ** 2
sig2_h = float(OZh["sdg2"].sum()) / int(OZh["N"]) / gbar_h ** 2
d_sig2 = sig2_g - sig2_h

pt = np.pi * tau
cot = np.cos(pt) / np.sin(pt)
T1g = pt * cot * np.real(Gg) / gbar_g
T1h = pt * cot * np.real(Gh) / gbar_h
T2 = -pt ** 2 * d_sig2 / 2.0

def bantla(x):
    return np.array([np.sum(ae[m] * x[m]) / ae[m].sum() for m in maskeler])

bT1g, bT1h, bT2 = bantla(T1g), bantla(T1h), bantla(T2)
bGg_re, bGh_re = bantla(np.real(Gg)), bantla(np.real(Gh))
bGg_abs, bGh_abs = bantla(np.abs(Gg)), bantla(np.abs(Gh))
bG2g_abs, bG2h_abs = bantla(np.abs(G2g)), bantla(np.abs(G2h))

print(f"\nσ_ε²: gerçek {sig2_g:.5f}  Hkeskin {sig2_h:.5f}  Δ = {d_sig2:+.5f}")
print("\nTAYLOR-TEŞHİSİ (YALNIZ YORUM; kesin hesap = F_KİN defteri):")
print(f"{'bant':>12} {'ΔT1=T1g−T1Hk':>13} {'T2(Δσ²)':>9} {'ΔT1+T2':>8} {'log F_KİN':>10}")
for k in range(9):
    print(f"{et[k]:>12} {bT1g[k]-bT1h[k]:13.4f} {bT2[k]:9.4f} "
          f"{bT1g[k]-bT1h[k]+bT2[k]:8.4f} {np.log(F2_b[k]):10.4f}")

print("\nĜ(τ) BANT PROFİLİ (H-W2 köprüsü; ae-ağırlıklı):")
print(f"{'bant':>12} {'Re Ĝ_g':>9} {'Re Ĝ_Hk':>9} {'|Ĝ_g|':>8} {'|Ĝ_Hk|':>8} "
      f"{'|Ĝ2_g|':>8} {'|Ĝ2_Hk|':>8}")
for k in range(9):
    print(f"{et[k]:>12} {bGg_re[k]:9.5f} {bGh_re[k]:9.5f} {bGg_abs[k]:8.5f} "
          f"{bGh_abs[k]:8.5f} {bG2g_abs[k]:8.5f} {bG2h_abs[k]:8.5f}")

# ---------- öz'ün DC/2ω anatomisi ----------
def oz_parcalar(OZ):
    N = int(OZ["N"])
    Sdc = OZ["Sdc"].sum(0) / N
    S2 = (OZ["S2re"].sum(0) + 1j * OZ["S2im"].sum(0)) / N
    return Sdc, S2

Sdc_g, S2_g = oz_parcalar(OZg)
Sdc_h, S2_h = oz_parcalar(OZh)
bdc_g = bantla(Sdc_g / np.sin(pt))       # ⟨sin(ωg/2)⟩/sin(πτ)  (DC/nominal)
bdc_h = bantla(Sdc_h / np.sin(pt))
b2_g = bantla(np.abs(S2_g) / np.sin(pt))
b2_h = bantla(np.abs(S2_h) / np.sin(pt))
print("\nÖZ ANATOMİSİ — DC ve 2ω parçaları (nominal sin(πτ)'ye oran):")
print(f"{'bant':>12} {'DC_g':>8} {'DC_Hk':>8} {'|2ω|_g':>8} {'|2ω|_Hk':>8}")
for k in range(9):
    print(f"{et[k]:>12} {bdc_g[k]:8.4f} {bdc_h[k]:8.4f} {b2_g[k]:8.4f} {b2_h[k]:8.4f}")

# ---------- H-W0: pencere içi L kayması ----------
d = np.load(S155 / "eta_son_t0.4_c4000.npz")
mid = np.asarray(d["mid"], float)
Lw = np.log(mid / TWO_PI)
L_min, L_maks, L_bar = float(Lw.min()), float(Lw.max()), float(Lw.mean())
dL_rms = float(Lw.std())
# bant profili eğimi d(log F2)/dτ (bitişik bant farkları, maks mutlak)
lf2 = np.log(F2_b[:8])
egim_f2 = np.max(np.abs(np.diff(lf2) / np.diff(tau_bar[:8])))
# saf zarf eğimi d r/dτ
egim_r = np.max(np.abs(np.diff(r_b[:8]) / np.diff(tau_bar[:8])))
tau_k = float(tau_bar[8])
dtau_rms = tau_k * dL_rms / L_bar
katki_dif = egim_f2 * dtau_rms                # gerçek−ikiz FARKINDA (ortak-mod düştü)
katki_kaba = egim_r * dtau_rms                # ortak-mod düşmese bile kaba üst sınır
bir_r = 1.0 - float(r_b[8])
print(f"\nH-W0 — PENCERE İÇİ L KAYMASI:")
print(f"  L_min={L_min:.5f}  L_maks={L_maks:.5f}  L̄={L_bar:.5f}  "
      f"ΔL_rms={dL_rms:.5f} ({100*dL_rms/L_bar:.2f}%)")
print(f"  δτ_rms(kuyruk) = τ̄·ΔL_rms/L̄ = {dtau_rms:.5f}")
print(f"  diferansiyel sınır (ortak-mod düşer): |dlogF_KİN/dτ|·δτ = "
      f"{katki_dif:.2e}  → (1−r)'nin %{100*katki_dif/bir_r:.2f}'i")
print(f"  kaba üst sınır (ortak-mod düşmese): |dr/dτ|·δτ = "
      f"{katki_kaba:.2e}  → (1−r)'nin %{100*katki_kaba/bir_r:.2f}'i")
hw0_yuzde = 100 * katki_kaba / bir_r
hw0 = "KAPANDI (<%1)" if hw0_yuzde < 1.0 else f"AÇIK (%{hw0_yuzde:.2f})"
print(f"  H-W0 hükmü: {hw0}")

np.savez_compressed(
    S185 / "K2_G_profil.npz", tau_bar=tau_bar, etiket=np.array(et),
    Gg_re=bGg_re, Gh_re=bGh_re, Gg_abs=bGg_abs, Gh_abs=bGh_abs,
    G2g_abs=bG2g_abs, G2h_abs=bG2h_abs, T1g=bT1g, T1h=bT1h, T2=bT2,
    dc_g=bdc_g, dc_h=bdc_h, iki_g=b2_g, iki_h=b2_h,
    Gg_cizgi=Gg, Gh_cizgi=Gh)
json.dump({
    "sha_onkayit": ONK["sha256"],
    "paylar_kuyruk": {"F_ANOMALI": kuyruk_pay[0], "F_KINEMATIK": kuyruk_pay[1],
                      "F_KOMSU_inv": kuyruk_pay[2]},
    "baskin": bask_ad, "cift_pay_ANOMxKOMSU": cift_pay,
    "sig2_eps": {"gercek": sig2_g, "Hkeskin": sig2_h, "delta": d_sig2},
    "H-W0": {"L_min": L_min, "L_maks": L_maks, "L_bar": L_bar,
             "dL_rms": dL_rms, "dtau_rms_kuyruk": dtau_rms,
             "katki_diferansiyel": katki_dif, "katki_kaba_ust": katki_kaba,
             "yuzde_1_minus_r": hw0_yuzde, "hukum": hw0},
}, open(S185 / "K2_anatomi.json", "w"), indent=1, ensure_ascii=False)
print("\n-> 185/K2_anatomi.json + K2_G_profil.npz  BİTTİ")
