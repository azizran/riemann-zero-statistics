# -*- coding: utf-8 -*-
"""
188f — ÖN-KAYITSIZ KEŞİF (hüküm DIŞI; yalnız KAYIT)
===================================================
Ön-kayıt SONRASI, harita görüldükten sonra sorulan sorular — hiçbir hükmü
değiştirmez, eşik kurtarmaz; bir sonraki kalemin adresini yazmak için:
 (a) 187 katman ızgarasında harita-ζ_i (HAVUZ, −Re Σ_s K) + jk se
     (187'nin "küçük se" ipucunun sınaması)
 (b) derinlik-eşli ζ toplamları: gerçek vs ikiz, τ' ≤ 0.90/1.00/1.10/1.20
 (c) blok-blok tepe konumları vs yerel L_b (L-kayması: Bragg, 1+τ_2, 1.15)
 (d) 8-bant: corr(z_r1, |ζ|_tam), yükseliş ayrışımı, τ'>1.30 kuyruk payı
 (e) naif birinci-mertebe bant çarpanı U_T = Σ c^ölç conj(mix)/Σ|mix|² vs u
 (f) ön-kayıtlı pencerelerde (Bragg + 8 uydu, ±0.006) κ toplamları ± jk
Çıktı: 188/K_kesif.json + ekran.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S = {d: SCR / d for d in ["155", "184", "185", "186", "187", "188"]}
L = 12.029593241726252
spec = importlib.util.spec_from_file_location(
    "b185", QM / "185_configs" / "185b_oz_muhasebe.py")
b185 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b185)


def jk(r):
    r = np.asarray(r)
    return np.sqrt(7 / 8 * ((r - r.mean(0)) ** 2).sum(0))


H = np.load(S["188"] / "harita_K_gercek.npz")
Hk = np.load(S["188"] / "harita_K_Hkeskin.npz")
A = json.load(open(S["188"] / "K3_analiz.json"))
k, n = H["kenar"], H["ncizgi"]
K, Kr = H["K"][-1], H["K_reps"][:, -1]
out = {}

# (a)
iz = [0.86, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30]
print("(a) 187 katman ızgarasında harita-ζ_i (HAVUZ; −Re Σ_s K ± jk):")
a = []
for lo, hi in zip(iz[:-1], iz[1:]):
    m = (k[:-1] >= lo - 1e-9) & (k[1:] <= hi + 1e-9)
    z, zr = K[m].sum(), Kr[:, m].sum(1)
    a.append({"katman": f"({lo:.2f},{hi:.2f}]", "zeta_i": float(-z.real),
              "se": float(jk(-zr.real)), "aci": float(np.degrees(np.angle(z))),
              "yogunluk": float(-z.real / (hi - lo)), "n": int(n[m].sum())})
    print(f"   ({lo:.2f},{hi:.2f}] ζ_i={-z.real:.4f}±{jk(-zr.real):.4f} "
          f"∠{np.degrees(np.angle(z)):+.1f}°  /Δτ'={-z.real/(hi-lo):.3f}  "
          f"n={n[m].sum()}")
out["a_katman"] = a

# (b)
kk, KH, KHr = Hk["kenar"], Hk["K"][-1], Hk["K_reps"][:, -1]
print("\n(b) derinlik-eşli ζ (HAVUZ |Σ K| ± jk): gerçek | ikiz")
b = []
for c in [0.90, 1.00, 1.10, 1.20]:
    m, mh = k[1:] <= c + 1e-9, kk[1:] <= c + 1e-9
    zg, zgr = K[m].sum(), Kr[:, m].sum(1)
    zh, zhr = KH[mh].sum(), KHr[:, mh].sum(1)
    b.append({"derinlik": c, "gercek": float(abs(zg)),
              "se_g": float(jk(np.abs(zgr))), "ikiz": float(abs(zh)),
              "se_h": float(jk(np.abs(zhr)))})
    print(f"   τ'≤{c:.2f}: {abs(zg):.4f}±{jk(np.abs(zgr)):.4f} | "
          f"{abs(zh):.4f}±{jk(np.abs(zhr)):.4f}")
out["b_derinlik_esli"] = b

# (c) blok-blok tepe konumları
PR = np.load(S["188"] / "harita_proj_gercek.npz")
G = np.load(S["184"] / "K1_gercek.npz")
OZ = np.load(S["185"] / "OZ_gercek.npz")
P = np.load(S["186"] / "G1_proj_gercek.npz")
ck = 2 * (P["re_bir"].sum(0) + 1j * P["im_bir"].sum(0)) / int(P["N"])
coz = b185.c_oz(OZ, G["aq"])
mix_all = ck - coz
tau = G["tau"]
win = (tau >= 0.45) & (tau < 0.86)
mix = mix_all[win]
mid = np.load(S["155"] / "eta_son_t0.4_c4000.npz")["mid"]
kj = np.linspace(0, len(mid), 9).astype(int)
o = 0.5 * (k[:-1] + k[1:])
print("\n(c) blok-blok tepe konumu (yalnız-blok K_HAVUZ, −Re) vs yerel L_b:")
c_ = []
for bb in range(8):
    Lb = float(np.log(mid[kj[bb]:kj[bb + 1]] / (2 * np.pi)).mean())
    C = 2 * (PR["re"][bb] + 1j * PR["im"][bb]) / PR["nb"][bb]
    kap = -(C @ np.conj(mix)).real / np.sum(np.abs(mix) ** 2)

    def pk(lo, hi):
        m = (o > lo) & (o < hi)
        return float(o[m][np.argmax(kap[m])])
    r = {"blok": bb, "L_b": Lb,
         "bragg_ong": Lb / L, "bragg_tepe": pk(0.985, 1.015),
         "t2_ong": (Lb + np.log(2)) / L, "t2_tepe": pk(1.04, 1.075),
         "t6_ong": (Lb + np.log(6)) / L, "t115_tepe": pk(1.13, 1.17)}
    c_.append(r)
    print(f"   blok {bb}: L_b={Lb:.4f}  Bragg öng {r['bragg_ong']:.4f} tepe "
          f"{r['bragg_tepe']:.4f} | 1+τ2 öng {r['t2_ong']:.4f} tepe "
          f"{r['t2_tepe']:.4f} | L_b+log6 öng {r['t6_ong']:.4f} tepe "
          f"{r['t115_tepe']:.4f}")
for ad, x, y in [("bragg", "bragg_ong", "bragg_tepe"),
                 ("t2", "t2_ong", "t2_tepe"), ("t6", "t6_ong", "t115_tepe")]:
    xs = np.array([r[x] for r in c_])
    ys = np.array([r[y] for r in c_])
    print(f"   {ad}: corr(öng, tepe) = {np.corrcoef(xs, ys)[0,1]:.3f}  "
          f"maks|tepe−öng| = {np.abs(ys-xs).max():.4f} (dilim 0.005)")
    out[f"c_{ad}_corr"] = float(np.corrcoef(xs, ys)[0, 1])
    out[f"c_{ad}_maksfark"] = float(np.abs(ys - xs).max())
out["c_bloklar"] = c_

# (d) 8-bant
v6 = A["vi"]
zf = np.array([r["zeta_tam_187c"] for r in v6])
z = np.array([r["z"] for r in v6])
zr1 = np.array([r["z_r1"] for r in v6])
za = np.array([r["z_art"] for r in v6])
d = {"corr_zr1_zetatam": float(np.corrcoef(zr1, zf)[0, 1]),
     "corr_z_zetatam": float(np.corrcoef(z, zf)[0, 1]),
     "oran_zr1_zetatam": (zr1 / zf).tolist(),
     "kuyruk_130_ustu": (zf - z).tolist(),
     "kuyruk_pay": ((zf - z) / zf).tolist(),
     "yukselis_tam_min_B7": float(zf[7] - zf.min()),
     "yukselis_z_B4_B7": float(z[7] - z[4]),
     "yukselis_zr1_B4_B7": float(zr1[7] - zr1[4]),
     "yukselis_zart_B4_B7": float(za[7] - za[4])}
out["d_8bant"] = d
print(f"\n(d) corr(z_r1, |ζ|tam) = {d['corr_zr1_zetatam']:.3f}; corr(z, |ζ|tam)"
      f" = {d['corr_z_zetatam']:.3f}; z_r1/|ζ|tam = "
      f"{np.round(zr1/zf, 3).tolist()}")
print(f"    τ'>1.30 kuyruk payı (|ζ|tam − z)/|ζ|tam = "
      f"{np.round((zf-z)/zf, 3).tolist()}")
print(f"    yükseliş: |ζ|tam (min→0.80-0.86) {d['yukselis_tam_min_B7']:.4f}; "
      f"B(0.65-0.70)→B(0.80-0.86): z {d['yukselis_z_B4_B7']:+.4f} = z_r1 "
      f"{d['yukselis_zr1_B4_B7']:+.4f} + z_art {d['yukselis_zart_B4_B7']:+.4f}")

# (e) naif U_T
col = b185.c_olculu(G)
ik = np.round(0.45 + 0.01 * np.arange(42), 2)
UT = []
for bb in range(41):
    m = (tau >= ik[bb]) & (tau < ik[bb + 1])
    UT.append(np.sum(col[m] * np.conj(mix_all[m])) /
              np.sum(np.abs(mix_all[m]) ** 2))
UT = np.array(UT)
u = np.array(A["u_mod"]) * np.exp(1j * np.radians(A["u_aci"]))
out["e_corr_u_UT"] = float(np.corrcoef(np.abs(u), np.abs(UT))[0, 1])
print(f"\n(e) naif bant çarpanı: corr(|u|, |U_T|) = {out['e_corr_u_UT']:.3f} "
      f"(41 ince bant)")

# (f) ön-kayıtlı pencere toplamları
ONK = json.load(open(S["188"] / "ONKAYIT_188.json"))
pen = ONK["uydular"]["pencere"]
kap, kapr = -K.real, -Kr.real
pencereler = {"Bragg": 1.0}
pencereler.update(ONK["uydular"]["konumlar"])
print("\n(f) ön-kayıtlı pencerelerde Σκ (HAVUZ; dilim merkezi ±0.006):")
f_, hepsi = {}, np.zeros(len(o), bool)
for ad, t in pencereler.items():
    m = np.abs(o - t) <= pen
    hepsi |= m
    f_[ad] = {"konum": float(t), "dilim": int(m.sum()),
              "kappa": float(kap[m].sum()), "se": float(jk(kapr[:, m].sum(1)))}
    print(f"   {ad:>9} τ'={t:.4f} ({m.sum()} dilim): Σκ = {kap[m].sum():.4f} "
          f"± {jk(kapr[:, m].sum(1)):.4f}")
top = kap.sum()
print(f"   tüm pencereler: {kap[hepsi].sum():.4f} / toplam {top:.4f} = "
      f"{kap[hepsi].sum()/top:.3f}  (dilim payı {hepsi.mean():.3f})")
f_["tum_pencere_pay"] = float(kap[hepsi].sum() / top)
f_["tum_pencere_dilim_pay"] = float(hepsi.mean())
m115 = (o > 1.14) & (o < 1.155)
print(f"   1.15-tepesi (1.1425–1.1525, pencere dışı): Σκ = {kap[m115].sum():.4f}"
      f" ± {jk(kapr[:, m115].sum(1)):.4f}")
f_["tepe_115"] = float(kap[m115].sum())
out["f_pencereler"] = f_

json.dump(out, open(S["188"] / "K_kesif.json", "w"), indent=1,
          ensure_ascii=False)
print(f"\n-> {S['188']/'K_kesif.json'}  BİTTİ")
