# -*- coding: utf-8 -*-
"""
186d — K3: KAPANIŞ (zincir r_pred; türetilmiş (c,α); sabit-nokta kalem-sayısalı)
================================================================================
ONKAYIT_186 kazanan kuralı: eşiği geçen güzergâh; İKİSİ DE GEÇEMEDİYSE zincir
χ²'si küçük olana YALNIZ TEŞHİS olarak takılır (hüküm ÖLDÜ, kurtarma yok).
Zincir: r_pred_b = (w_g^öz + M^pred·m_Hk)/(w_Hk^öz + m_Hk); 184 defterine χ²/dof;
1−c·τ^α fiti (α-tarama [0.05,6.0) adım 0.001, kapalı-form c — 185d AYNEN).
SABİT-NOKTA: G1 ölse de iki koşu (ρ=yasa, ρ≡1) doğrusallaştırılmış haritayı
verir — kalem-sayısalı burada (bölüm rapora; yeni koşu YOK, mevcut sayılardan).
Çıktı: 186/K3_sonuc.json + ekran.
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
    "b186", QM / "186_configs" / "186b_g1_oztutarlilik.py")
b186 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b186)
b185 = b186.b185

ONK = json.load(open(S186 / "ONKAYIT_186.json"))
sha = hashlib.sha256(
    (QM / "186_configs" / "186a_onkayit.py").read_bytes()).hexdigest()
assert sha == ONK["sha256"], "ON-KAYIT SHA UYUMSUZ"

G1S = json.load(open(S186 / "G1_sonuc.json"))
G2S = json.load(open(S186 / "G2_sonuc.json"))
G = np.load(S184 / "K1_gercek.npz")
H = np.load(S184 / "K1_Hkeskin.npz")
OZg = np.load(S185 / "OZ_gercek.npz")
OZh = np.load(S185 / "OZ_Hkeskin.npz")
Pg = np.load(S186 / "G1_proj_gercek.npz")
K1f = np.load(S185 / "K1_faktorler.npz", allow_pickle=True)

tau = np.asarray(G["tau"], float)
aq = np.asarray(G["aq"], float)
ae = np.asarray(G["aq_eff"], float)
kenar = ONK["kenar"]
maskeler = [(tau >= kenar[b]) & (tau < kenar[b + 1]) for b in range(8)]
rhos = b186.rho_vektorleri(tau)
tau_bar = K1f["tau_bar"][:8]
r_b, sig_r = K1f["r"][:8], K1f["sig_r"][:8]

print("=" * 78)
print(f"186d / K3 KAPANIŞ  [on-kayit sha {ONK['sha256'][:12]}]")
print("=" * 78)

chi_G1 = G1S["G1_chi2dof"]
chi_G2 = G2S["chi2dof_M"]
gecen = [ad for ad, c in [("G1", chi_G1), ("G2", chi_G2)] if c <= 2.0]
if gecen:
    kazanan = min(gecen, key=lambda a: {"G1": chi_G1, "G2": chi_G2}[a])
    statu = "KAZANAN"
else:
    kazanan = "G1" if chi_G1 <= chi_G2 else "G2"
    statu = "TEŞHİS (hiçbiri eşiği geçemedi; kurtarma yok)"
print(f"\nGüzergâh χ²/dof: G1 = {chi_G1:.2f}, G2 = {chi_G2:.2f}  →  "
      f"zincir {kazanan}'e {statu} olarak takılır.")


def zincir_defteri(Mp_fn):
    """Mp_fn(disari) -> (8,) M^pred; zincir r_pred tam + loo, jk se."""
    def tek(disari):
        og = np.abs(b185.c_oz(OZg, aq, disari))
        oh = np.abs(b185.c_oz(OZh, aq, disari))
        ah = np.abs(b185.c_olculu(H, disari))
        wog = np.array([og[m].sum() / ae[m].sum() for m in maskeler])
        woh = np.array([oh[m].sum() / ae[m].sum() for m in maskeler])
        wh = np.array([ah[m].sum() / ae[m].sum() for m in maskeler])
        mh = wh - woh
        Mp = Mp_fn(disari)
        return (wog + Mp * mh) / (woh + mh)
    tam = tek(-1)
    reps = np.array([tek(i) for i in range(NJACK)])
    return tam, b185.jk_se(reps)


def Mp_g1(ad):
    def fn(disari):
        og = np.abs(b185.c_oz(OZg, aq, disari))
        oh = np.abs(b185.c_oz(OZh, aq, disari))
        ah = np.abs(b185.c_olculu(H, disari))
        ya = np.abs(b186.c_ya(Pg, ad, disari))
        out = np.zeros(8)
        for k, m in enumerate(maskeler):
            sae = ae[m].sum()
            mh = (ah[m].sum() - oh[m].sum()) / sae
            mp = (ya[m].sum() - (rhos[ad][m] * og[m]).sum()) / sae
            out[k] = mp / mh
        return out
    return fn


if kazanan == "G1":
    Mp_fn = Mp_g1("yasa")
else:  # G2 köprüsü (loo: v bantları OZ bloklarından)
    def Mp_fn(disari):
        def vb(D):
            N, nb = int(D["N"]), D["nb"]
            if disari < 0:
                re, im, n = D["Gre"].sum(0), D["Gim"].sum(0), N
            else:
                re = D["Gre"].sum(0) - D["Gre"][disari]
                im = D["Gim"].sum(0) - D["Gim"][disari]
                n = N - nb[disari]
            aG = np.abs(re + 1j * im) / n
            return np.array([aG[m].sum() / (float(D["gbar"]) * ae[m].sum())
                             for m in maskeler])
        vg, vh = vb(OZg), vb(OZh)
        return ((1.017 - 0.884 * vg) / (1.017 - 0.884 * vh)) ** 2

r_pred, se_rp = zincir_defteri(Mp_fn)
z_r = (r_pred - r_b) / sig_r
chi_zincir = float(np.sum(z_r ** 2) / 8)

print(f"\nZİNCİR DEFTERİ (r_pred = (w_g^öz + M^pred·m_Hk)/(w_Hk^öz + m_Hk); "
      f"{kazanan} M^pred'i; σ = 184 defter-σ):")
for k in range(8):
    print(f"  τ̄={tau_bar[k]:.3f}  r_ölç={r_b[k]:.4f}±{sig_r[k]:.4f}  "
          f"r_pred={r_pred[k]:.4f}±{se_rp[k]:.4f}  z={z_r[k]:+.1f}")
print(f"  zincir bant-χ²/dof = {chi_zincir:.2f}  (eşik 2)")
zincir_hukum = "KAPANDI" if chi_zincir <= 2.0 else "KAPANMADI"
print(f"  K3 ZİNCİR HÜKMÜ: {zincir_hukum}")


# ---------- güç-yasası fitleri (185d makinesi AYNEN) ----------
def guc_fit(y, sig):
    wgt = 1.0 / sig ** 2
    best = None
    for a in np.arange(0.05, 6.0, 0.001):
        ta = tau_bar ** a
        c = np.sum(wgt * (1.0 - y) * ta) / np.sum(wgt * ta * ta)
        chi = np.sum(wgt * (y - (1.0 - c * ta)) ** 2)
        if best is None or chi < best[2]:
            best = (float(c), float(a), float(chi))
    c, a, chi = best
    return c, a, chi / (8 - 2)


c_olc, a_olc, x_olc = guc_fit(r_b, sig_r)
c_p, a_p, x_p = guc_fit(r_pred, sig_r)
print("\nGÜÇ-YASASI FİTLERİ  1 − c·τ^α  (fit-χ²/dof, dof=6):")
print(f"  ölçülü (kontrol) : c = {c_olc:+.4f}  α = {a_olc:.3f}  "
      f"χ²/dof = {x_olc:.2f}   [184: c=0.149 α=1.30 χ²/dof=0.34]")
print(f"  zincir r_pred    : c = {c_p:+.4f}  α = {a_p:.3f}  χ²/dof = {x_p:.2f}")

# yan: h_yasa (ikiz-çekirdek yan-lehçesi) zinciri — bilgi
rp_h = np.array(G1S["h_yasa"]["r_pred"])
z_h = (rp_h - r_b) / sig_r
chi_h = float(np.sum(z_h ** 2) / 8)
c_h, a_h, x_h = guc_fit(rp_h, sig_r)
print(f"  yan (Hk-çekird.) : c = {c_h:+.4f}  α = {a_h:.3f}  χ²/dof = {x_h:.2f}"
      f"   [zincir χ²/dof = {chi_h:.2f}]")

# ---------- SABİT-NOKTA KALEM-SAYISALI (yeni koşu yok) ----------
# Aile ρ_c(τ) = 1 − c·τ^1.30 (α donmuş). İki koşu: c=0 (g_bir), c=0.149 (g_yasa)
# → bant başına doğrusallaştırılmış harita r_pred,b(c) = r0_b + s_b·c.
# Sabit nokta (bant b): 1 − c*·τ̄_b^1.30 = r0_b + s_b·c*  →
#   c*_b = (1 − r0_b)/(s_b + τ̄_b^1.30)
r0 = np.array(G1S["g_bir"]["r_pred"])
r1 = np.array(G1S["g_yasa"]["r_pred"])
s_b = (r1 - r0) / 0.149
t13 = tau_bar ** 1.30
c_star = (1.0 - r0) / (s_b + t13)
rho_star = 1.0 - c_star * t13
kontraksiyon = np.abs(s_b) / t13   # |dr_pred/dc| / |dρ/dc| — haritanın yerel eğimi
print("\nSABİT-NOKTA KALEM-SAYISALI (ρ_c = 1−c·τ^1.30 ailesi; iki koşudan "
      "doğrusallaştırma):")
print(f"{'τ̄':>6} {'r_pred(c=0)':>12} {'s_b=dr/dc':>10} {'c*_b':>8} "
      f"{'ρ*=r*':>8} {'|eğim|':>7}")
for k in range(8):
    print(f"{tau_bar[k]:6.3f} {r0[k]:12.4f} {s_b[k]:10.3f} {c_star[k]:8.3f} "
          f"{rho_star[k]:8.4f} {kontraksiyon[k]:7.3f}")
print("  → sabit nokta HER bantta r*>1 tarafında (anti-zarf); harita güçlü "
      "kontraktif (|eğim|≪1) ama yanlış merkeze çekiyor: güç-yasası zarfı "
      "kendi karışımından ÖZ-ÜRETMİYOR.")

json.dump({
    "sha_onkayit": ONK["sha256"],
    "kazanan": kazanan, "statu": statu,
    "chi2dof": {"G1": chi_G1, "G2": chi_G2},
    "tau_bar": tau_bar.tolist(),
    "r_olc": r_b.tolist(), "sig_r": sig_r.tolist(),
    "r_pred": r_pred.tolist(), "se_r_pred": se_rp.tolist(),
    "z": z_r.tolist(), "zincir_chi2dof": chi_zincir,
    "zincir_hukum": zincir_hukum,
    "fit": {"olculu": [c_olc, a_olc, x_olc], "zincir": [c_p, a_p, x_p],
            "yan_Hk": [c_h, a_h, x_h], "hedef_184": [0.149, 1.30, 0.34]},
    "yan_Hk_zincir_chi2dof": chi_h,
    "sabit_nokta": {"r0": r0.tolist(), "s_b": s_b.tolist(),
                    "c_star": c_star.tolist(), "rho_star": rho_star.tolist(),
                    "egim": kontraksiyon.tolist()},
}, open(S186 / "K3_sonuc.json", "w"), indent=1, ensure_ascii=False)
print(f"\n-> {S186/'K3_sonuc.json'}  BİTTİ")
