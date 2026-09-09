# -*- coding: utf-8 -*-
"""
186e — K4 (bonus): HA4 SINAVI (karışım-eşitsizliği açıklanıyor mu?)
===================================================================
Kazanan yok (G1/G2 ikisi de öldü) → K0 kuralının K3-analojisi: χ²'si küçük
güzergâhın (G1) makinesi TEŞHİS olarak; G2 köprü varyantı yan sütun.
  Ölçülü hedef: M_HA4,b = m_HA4,b/m_Hk,b;  m_HA4 = w_HA4 − w_HA4^öz,erfc
                (env=0.5·erfc((τ−0.68)/0.125), 185-K4 DONMUŞ)
  G1-env (birincil teşhis): ds^ya(ρ=env) Hkeskin kinematiğinde →
                m^pred = w^ya − Σenv·â^öz_Hk/Σae;  yan: aynı koşu gerçek kin.
  G2-köprü (yan): M^pred = [(1.017−0.884·v_HA4)/(1.017−0.884·v_Hk)]²
Çıktı: 186/K4_HA4.json + ekran.
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

A = np.load(S184 / "K1_HA4.npz")
H = np.load(S184 / "K1_Hkeskin.npz")
OZa = np.load(S185 / "OZ_HA4.npz")
OZh = np.load(S185 / "OZ_Hkeskin.npz")
Pg = np.load(S186 / "G1_proj_gercek.npz")
Ph = np.load(S186 / "G1_proj_Hkeskin.npz")
G2S = json.load(open(S186 / "G2_sonuc.json"))
K1f = np.load(S185 / "K1_faktorler.npz", allow_pickle=True)

G = np.load(S184 / "K1_gercek.npz")
tau = np.asarray(G["tau"], float)
aq = np.asarray(G["aq"], float)
ae = np.asarray(G["aq_eff"], float)
kenar = ONK["kenar"]
maskeler = [(tau >= kenar[b]) & (tau < kenar[b + 1]) for b in range(8)]
rhos = b186.rho_vektorleri(tau)
env = rhos["env"]
tau_bar = K1f["tau_bar"][:8]

print("=" * 78)
print(f"186e / K4 HA4 SINAVI  [on-kayit sha {ONK['sha256'][:12]}]")
print("=" * 78)


OZg = np.load(S185 / "OZ_gercek.npz")


def defter2(disari):
    aa = np.abs(b185.c_olculu(A, disari))
    ah = np.abs(b185.c_olculu(H, disari))
    oa = np.abs(b185.c_oz(OZa, aq, disari))
    oh = np.abs(b185.c_oz(OZh, aq, disari))
    og = np.abs(b185.c_oz(OZg, aq, disari))
    ya_h = np.abs(b186.c_ya(Ph, "env", disari))
    ya_g = np.abs(b186.c_ya(Pg, "env", disari))
    out = []
    for m in maskeler:
        sae = ae[m].sum()
        w_a = aa[m].sum() / sae
        woz_a = (env[m] * oa[m]).sum() / sae
        m_a = w_a - woz_a
        m_h = (ah[m].sum() - oh[m].sum()) / sae
        mp_h = (ya_h[m].sum() - (env[m] * oh[m]).sum()) / sae
        mp_g = (ya_g[m].sum() - (env[m] * og[m]).sum()) / sae
        out.append([w_a, woz_a, m_a, m_h, m_a / m_h, mp_h, mp_h / m_h,
                    mp_g, mp_g / m_h])
    return np.array(out)


tam = defter2(-1)
reps = np.array([defter2(i) for i in range(NJACK)])
se = b185.jk_se(reps)

w_a, woz_a, m_a, m_h, M_a = [tam[:, j] for j in range(5)]
mp_h, Mp_h = tam[:, 5], tam[:, 6]
mp_g, Mp_g = tam[:, 7], tam[:, 8]

# birincil σ(M_HA4): bağımsız yayılım (defter konvansiyonu)
se_ma = np.sqrt(se[:, 0] ** 2 + se[:, 1] ** 2)
se_mh = se[:, 3]
sigMa = M_a * np.sqrt((se_ma / m_a) ** 2 + (se_mh / m_h) ** 2)

print("\nÖLÇÜLÜ HA4 KARIŞIM DEFTERİ (m_HA4 = w_HA4 − w^öz,erfc; "
      "M_HA4 = m_HA4/m_Hk ±defter-σ):")
for k in range(8):
    print(f"  τ̄={tau_bar[k]:.3f}  w_HA4={w_a[k]:.4f}  w^öz,erfc={woz_a[k]:.4f}  "
          f"m_HA4={m_a[k]:.4f}  m_Hk={m_h[k]:.4f}  "
          f"M_HA4={M_a[k]:.4f}±{sigMa[k]:.4f}")

z_h = (Mp_h - M_a) / sigMa
chi_h = float(np.sum(z_h ** 2) / 8)
z_g = (Mp_g - M_a) / sigMa
chi_g = float(np.sum(z_g ** 2) / 8)
print("\nG1-env TEŞHİSİ (birincil: Hk kinematiği; yan: gerçek kinematiği):")
for k in range(8):
    print(f"  τ̄={tau_bar[k]:.3f}  M^pred(Hk)={Mp_h[k]:.4f}±{se[k,6]:.4f} "
          f"z={z_h[k]:+.1f}   M^pred(g)={Mp_g[k]:.4f}±{se[k,8]:.4f} "
          f"z={z_g[k]:+.1f}")
print(f"  bant-χ²/dof: Hk-kin = {chi_h:.2f}  |  gerçek-kin = {chi_g:.2f}  (eşik 2)")

# G2-köprü yan sütunu
va = np.array(G2S["v"]["HA4"])
vh = np.array(G2S["v"]["Hkeskin"])
Mp_v = ((1.017 - 0.884 * va) / (1.017 - 0.884 * vh)) ** 2
z_v = (Mp_v - M_a) / sigMa
chi_v = float(np.sum(z_v ** 2) / 8)
print(f"\nG2-köprü yan sütunu: M^pred = {np.round(Mp_v, 3).tolist()}")
print(f"  z = {np.round(z_v, 1).tolist()};  bant-χ²/dof = {chi_v:.2f}")

# karışım-eşitsizliği sorusu (orta-τ ~%8-10)
print("\nKARIŞIM-EŞİTSİZLİĞİ (185-K4 çatlağının adresi):")
print(f"  ölçülü M_HA4 (orta-τ, ilk 4 bant): {np.round(M_a[:4], 3).tolist()} "
      f"(>1: HA4 karışımı Hkeskin'den BÜYÜK)")
print(f"  G1-env öngörüsü (Hk-kin):          {np.round(Mp_h[:4], 3).tolist()}")
esitsizlik_aciklandi = bool(np.all(Mp_h[:4] > 1.0))
print(f"  → eşitsizliğin YÖNÜ (M>1) yakalandı mı? "
      f"{'EVET' if esitsizlik_aciklandi else 'HAYIR'}")

hukum = "MÜHÜR" if chi_h <= 2.0 else "ÇATLAK"
print(f"\nK4 HÜKMÜ (bonus; birincil teşhis G1-env/Hk-kin): {hukum}  "
      f"[χ²/dof = {chi_h:.2f}]")

json.dump({
    "sha_onkayit": ONK["sha256"], "tau_bar": tau_bar.tolist(),
    "M_HA4_olc": M_a.tolist(), "sigM_defter": sigMa.tolist(),
    "m_HA4": m_a.tolist(), "m_Hk": m_h.tolist(),
    "G1env_Hkkin": {"M_pred": Mp_h.tolist(), "se": se[:, 6].tolist(),
                    "z": z_h.tolist(), "chi2dof": chi_h},
    "G1env_gkin": {"M_pred": Mp_g.tolist(), "se": se[:, 8].tolist(),
                   "z": z_g.tolist(), "chi2dof": chi_g},
    "G2kopru": {"M_pred": Mp_v.tolist(), "z": z_v.tolist(), "chi2dof": chi_v},
    "esitsizlik_yonu_yakalandi": esitsizlik_aciklandi,
    "K4_hukum": hukum,
}, open(S186 / "K4_HA4.json", "w"), indent=1, ensure_ascii=False)
print(f"\n-> {S186/'K4_HA4.json'}  BİTTİ")
