# -*- coding: utf-8 -*-
"""
186c — K2: G2 v-KANALI (yüksek-τ İLK ölçüm) + KISIT KÖPRÜSÜ
===========================================================
ONKAYIT_186 dondurdu:
  v̂_q = |Ĝ_q|/(ḡ·ae_q);  Ĝ_q = ⟨(g−ḡ)e^{−iω_q m}⟩  (185 OZ blokları AYNEN)
  bant: v_b = Σ_b|Ĝ_q|/(ḡ·Σ_b ae_q)  (toplam-oranı); iki denizde ÖZDEŞ tanım
  köprü: M^pred_b = [(1.017−0.884·v_g,b)/(1.017−0.884·v_Hk,b)]²
  yan: her denizde √w_b vs 1.017−0.884·v_b (kısıt-seviye);  köprü² vs r_b.
Çıktı: 186/G2_sonuc.json + ekran defteri.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S184, S185, S186 = SCR / "184", SCR / "185", SCR / "186"
NJACK = 8

ONK = json.load(open(S186 / "ONKAYIT_186.json"))
sha = hashlib.sha256(
    (QM / "186_configs" / "186a_onkayit.py").read_bytes()).hexdigest()
if sha != ONK["sha256"]:
    raise SystemExit("ON-KAYIT SHA UYUMSUZ")

G = np.load(S184 / "K1_gercek.npz")
tau = np.asarray(G["tau"], float)
ae = np.asarray(G["aq_eff"], float)
kenar = ONK["kenar"]
maskeler = [(tau >= kenar[b]) & (tau < kenar[b + 1]) for b in range(8)]
et = [f"{kenar[b]:.2f}-{kenar[b+1]:.2f}" for b in range(8)]

OZ = {gaz: np.load(S185 / f"OZ_{gaz}.npz")
      for gaz in ["gercek", "Hkeskin", "HA4"]}


def G_hat(D, disari=-1):
    N = int(D["N"])
    nb = D["nb"]
    if disari < 0:
        re, im, n = D["Gre"].sum(0), D["Gim"].sum(0), N
    else:
        re = D["Gre"].sum(0) - D["Gre"][disari]
        im = D["Gim"].sum(0) - D["Gim"][disari]
        n = N - nb[disari]
    return (re + 1j * im) / n


def v_bant(D, disari=-1):
    gbar = float(D["gbar"])
    aG = np.abs(G_hat(D, disari))
    return np.array([aG[m].sum() / (gbar * ae[m].sum()) for m in maskeler])


def jk_se(reps):
    reps = np.asarray(reps)
    return np.sqrt((NJACK - 1) / NJACK * np.sum((reps - reps.mean(0)) ** 2, 0))


def kopru(vg, vh):
    return ((1.017 - 0.884 * vg) / (1.017 - 0.884 * vh)) ** 2


print("=" * 78)
print(f"186c / K2 G2 v-KANALI  [on-kayit {ONK['zaman']} sha {ONK['sha256'][:12]}]")
print("=" * 78)

v = {gaz: v_bant(OZ[gaz]) for gaz in OZ}
v_reps = {gaz: np.array([v_bant(OZ[gaz], i) for i in range(NJACK)]) for gaz in OZ}
v_se = {gaz: jk_se(v_reps[gaz]) for gaz in OZ}

# ölçülü M defteri + birincil σ (185 kayıtlarından; 186b ile aynı formül)
K1f = np.load(S185 / "K1_faktorler.npz", allow_pickle=True)
tau_bar = K1f["tau_bar"][:8]
wg, wh = K1f["wg"][:8], K1f["wh"][:8]
wog, woh = K1f["wog"][:8], K1f["woh"][:8]
r_b, sig_r = K1f["r"][:8], K1f["sig_r"][:8]
se185 = K1f["se_tablo"]
mg, mh = wg - wog, wh - woh
M_b = mg / mh
se_mg = np.sqrt(se185[:8, 7] ** 2 + se185[:8, 9] ** 2)
se_mh = np.sqrt(se185[:8, 8] ** 2 + se185[:8, 10] ** 2)
sigM = M_b * np.sqrt((se_mg / mg) ** 2 + (se_mh / mh) ** 2)

# --- v profilleri (bağımsız bulgu) ---
print("\nv(τ) PROFİLLERİ — yüksek-τ İLK ölçüm (v̂=|Ĝ|/(ḡ·ae); ±jk se):")
print(f"{'bant':>12} {'τ̄':>6} {'v_g':>14} {'v_Hk':>14} {'v_HA4':>14} "
      f"{'V=v_g/v_Hk':>12}")
V_b = v["gercek"] / v["Hkeskin"]
V_reps = v_reps["gercek"] / v_reps["Hkeskin"]
V_se = jk_se(V_reps)
for k in range(8):
    print(f"{et[k]:>12} {tau_bar[k]:6.3f} "
          f"{v['gercek'][k]:7.4f}±{v_se['gercek'][k]:.4f} "
          f"{v['Hkeskin'][k]:7.4f}±{v_se['Hkeskin'][k]:.4f} "
          f"{v['HA4'][k]:7.4f}±{v_se['HA4'][k]:.4f} "
          f"{V_b[k]:7.4f}±{V_se[k]:.4f}")

# --- kısıt köprüsü ---
Mp = kopru(v["gercek"], v["Hkeskin"])
Mp_reps = kopru(v_reps["gercek"], v_reps["Hkeskin"])
Mp_se = jk_se(Mp_reps)
z = (Mp - M_b) / sigM
chi = float(np.sum(z ** 2) / 8)
print("\nKISIT KÖPRÜSÜ M^pred vs ölçülü M (σ birincil = defter-seM):")
for k in range(8):
    print(f"{et[k]:>12} τ̄={tau_bar[k]:.3f}  M^pred={Mp[k]:.4f}±{Mp_se[k]:.4f}  "
          f"M_ölç={M_b[k]:.4f}±{sigM[k]:.4f}  z={z[k]:+.1f}")
print(f"  bant-χ²/dof = {chi:.2f}  (eşik 2)")

# --- yan: köprü vs r (köprü cebirsel olarak w_g/w_Hk öngörür) ---
z_r = (Mp - r_b) / sig_r
chi_r = float(np.sum(z_r ** 2) / 8)
print(f"\nYAN-OKUMA köprü² vs r_b (bilgi): χ²/dof = {chi_r:.2f}; "
      f"z = {np.round(z_r, 1).tolist()}")

# --- yan: kısıt-seviye kontrolü her denizde ---
print("\nKISIT-SEVİYE KONTROLÜ √w_b vs 1.017−0.884·v_b (bilgi; konvansiyon farkı "
      "burada görünür):")
seviye = {}
for gaz, wb in [("gercek", wg), ("Hkeskin", wh)]:
    sw = np.sqrt(wb)
    pred = 1.017 - 0.884 * v[gaz]
    fark = sw - pred
    seviye[gaz] = {"sqrt_w": sw.tolist(), "kisit_pred": pred.tolist(),
                   "fark": fark.tolist()}
    print(f"  [{gaz}] √w = {np.round(sw, 3).tolist()}")
    print(f"          k. = {np.round(pred, 3).tolist()}  "
          f"(fark rms {np.sqrt((fark**2).mean()):.3f})")

# --- Not 2 doyum kıyası (düşük-τ ucu) ---
print(f"\nNOT-2 KIYAS NOTU: v_g(τ̄=0.475) = {v['gercek'][0]:.4f}±"
      f"{v_se['gercek'][0]:.4f}; Not 2 doyumu ~0.55-0.60 (@τ≈0.4, farklı "
      f"konvansiyon olabilir). ORAN lehçesi birincil: V(0.475) = "
      f"{V_b[0]:.4f}±{V_se[0]:.4f}.")

hukum = "MÜHÜR" if chi <= 2.0 else "ÖLDÜ"
print(f"\nG2 HÜKMÜ (köprü M^pred, birincil): {hukum}  [χ²/dof = {chi:.2f}]")

json.dump({
    "sha_onkayit": ONK["sha256"], "tau_bar": tau_bar.tolist(),
    "v": {g: v[g].tolist() for g in v},
    "v_se": {g: v_se[g].tolist() for g in v},
    "V_oran": V_b.tolist(), "V_se": V_se.tolist(),
    "M_olc": M_b.tolist(), "sigM_defter": sigM.tolist(),
    "M_pred": Mp.tolist(), "M_pred_se": Mp_se.tolist(),
    "z": z.tolist(), "chi2dof_M": chi,
    "yan_r": {"z": z_r.tolist(), "chi2dof": chi_r},
    "kisit_seviye": seviye,
    "G2_hukum": hukum,
}, open(S186 / "G2_sonuc.json", "w"), indent=1, ensure_ascii=False)
print(f"\n-> {S186/'G2_sonuc.json'}  BİTTİ")
