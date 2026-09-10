# -*- coding: utf-8 -*-
"""
187f — K4: HA4 %26 KAPANIŞI (defter aritmetiği; H-F4)
=====================================================
ONKAYIT_187 dondurdu:
  BİRİNCİL: M^pred_HA4,b = m^env_HA4,b(τ_c=0.86; HA4-KİNEMATİĞİ) /
            m_Hk,b(ölçülü)   — "erfc-ağırlıklı pencere-içi karışım,
            İPTALSİZ; payda = Hkeskin'in tam karışımı" (KALEM AYNEN).
  YAN SÜTUN: aynı defter, env-katmanlar τ_c=1.10'a eklenmiş (tam-env).
  Hedef: 186e ölçülü M_HA4 = 1.2640→1.0190 (aynı makinede yeniden üretilir).
  z = (M^pred − M_ölç)/σ_defter(M_ölç);  bant-χ²/dof ≤ 2 → KAPANDI.

Çıktı: 187/K4_HA4.json + ekran defteri.
"""
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S184, S185, S186, S187 = (SCR / d for d in ["184", "185", "186", "187"])
NJACK = 8

spec = importlib.util.spec_from_file_location(
    "b185", QM / "185_configs" / "185b_oz_muhasebe.py")
b185 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b185)
spec7 = importlib.util.spec_from_file_location(
    "b187", QM / "187_configs" / "187b_katman_defteri.py")
b187 = importlib.util.module_from_spec(spec7)
spec7.loader.exec_module(b187)


def env_f(t):
    return 0.5 * math.erfc((t - 0.68) / 0.125)


if __name__ == "__main__":
    ONK = b187.onkayit()
    print("=" * 78)
    print(f"187f / K4 HA4 KAPANIŞI  [on-kayit sha {ONK['sha256'][:12]}]")
    print("=" * 78, flush=True)

    G = np.load(S184 / "K1_gercek.npz")
    tau = np.asarray(G["tau"], float)
    aq = np.asarray(G["aq"], float)
    ae = np.asarray(G["aq_eff"], float)
    env_win = np.array([env_f(t) for t in tau])
    A184 = np.load(S184 / "K1_HA4.npz")
    H184 = np.load(S184 / "K1_Hkeskin.npz")
    OZa = np.load(S185 / "OZ_HA4.npz")
    OZh = np.load(S185 / "OZ_Hkeskin.npz")
    KPa = np.load(S187 / "katmanproj_HA4.npz")
    NKAT = KPa["re"].shape[0]

    kenar_b = json.load(open(S186 / "ONKAYIT_186.json"))["kenar"]
    maskeler = [(tau >= kenar_b[b]) & (tau < kenar_b[b + 1]) for b in range(8)]
    maskeler.append((tau >= kenar_b[0]) & (tau < kenar_b[-1]))
    et = [f"{kenar_b[b]:.2f}-{kenar_b[b+1]:.2f}" for b in range(8)] + ["HAVUZ"]

    def cplx(reb, imb, disari, N, nb):
        if disari < 0:
            return 2.0 * (reb.sum(0) + 1j * imb.sum(0)) / N
        return 2.0 * ((reb.sum(0) - reb[disari]) +
                      1j * (imb.sum(0) - imb[disari])) / (N - nb[disari])

    def defter(disari):
        """sütunlar: w_HA4, wöz_erfc, w_Hk, wöz_Hk, m_env_kesik, m_env_tam"""
        cA = b185.c_olculu(A184, disari)
        cH = b185.c_olculu(H184, disari)
        ozA = env_win * np.abs(b185.c_oz(OZa, aq, disari))
        ozH = np.abs(b185.c_oz(OZh, aq, disari))
        ck = cplx(KPa["reT"], KPa["imT"], disari, int(KPa["N"]), KPa["nb"])
        ckat = [cplx(KPa["re"][i], KPa["im"][i], disari, int(KPa["N"]),
                     KPa["nb"]) for i in range(NKAT)]
        out = np.zeros((9, 6))
        for k, m in enumerate(maskeler):
            sae = ae[m].sum()
            out[k, 0] = np.abs(cA[m]).sum() / sae
            out[k, 1] = ozA[m].sum() / sae
            out[k, 2] = np.abs(cH[m]).sum() / sae
            out[k, 3] = ozH[m].sum() / sae
            cya = ck[m].copy()
            out[k, 4] = (np.abs(cya).sum() - ozA[m].sum()) / sae
            for i in range(NKAT):
                cya = cya + ckat[i][m]
            out[k, 5] = (np.abs(cya).sum() - ozA[m].sum()) / sae
        return out

    tam = defter(-1)
    reps = np.array([defter(i) for i in range(NJACK)])
    se = b185.jk_se(reps)

    wA, ozA_b, wH, ozH_b, m_env_kes, m_env_tam = [tam[:, j] for j in range(6)]
    m_HA4 = wA - ozA_b
    m_Hk = wH - ozH_b
    M_olc = m_HA4 / m_Hk
    # defter-σ (bağımsız yayılım; 186e AYNEN)
    se_mA = np.sqrt(se[:, 0] ** 2 + se[:, 1] ** 2)
    se_mH = np.sqrt(se[:, 2] ** 2 + se[:, 3] ** 2)
    sigM = M_olc * np.sqrt((se_mA / m_HA4) ** 2 + (se_mH / m_Hk) ** 2)

    M_pred = m_env_kes / m_Hk
    M_pred_tam = m_env_tam / m_Hk
    se_Mp = b185.jk_se(reps[:, :, 4] / (reps[:, :, 2] - reps[:, :, 3]))
    se_Mpt = b185.jk_se(reps[:, :, 5] / (reps[:, :, 2] - reps[:, :, 3]))
    z = (M_pred - M_olc) / sigM
    z_tam = (M_pred_tam - M_olc) / sigM
    chi = float(np.sum(z[:8] ** 2) / 8)
    chi_tam = float(np.sum(z_tam[:8] ** 2) / 8)

    g186 = json.load(open(S186 / "K4_HA4.json"))
    M186 = np.array(g186["M_HA4_olc"][:8])
    print("\nÖLÇÜLÜ M_HA4 DEFTERİ (186e yeniden üretimi — kontrol):")
    for k in range(8):
        print(f"  {et[k]:>10} w_HA4={wA[k]:.4f} w^öz,erfc={ozA_b[k]:.4f} "
              f"m_HA4={m_HA4[k]:.4f} m_Hk={m_Hk[k]:.4f} "
              f"M_HA4={M_olc[k]:.4f}±{sigM[k]:.4f}")
    muhur_186e = float(np.abs(M_olc[:8] - M186).max())
    print(f"  MÜHÜR: 186e M_HA4 defterine maks|Δ| = {muhur_186e:.2e}")

    print("\nK4 DEFTER ARİTMETİĞİ (birincil: iptalsiz pencere-içi erfc "
          "karışımı, HA4-kinematiği):")
    print(f"{'bant':>10} {'M_ölç':>16} {'M^pred(0.86)':>16} {'z':>6}  "
          f"{'M^pred(1.10 tam)':>17} {'z_tam':>6}")
    for k in range(9):
        print(f"{et[k]:>10} {M_olc[k]:7.4f}±{sigM[k]:.4f} "
              f"{M_pred[k]:7.4f}±{se_Mp[k]:.4f} {z[k]:+6.1f}  "
              f"{M_pred_tam[k]:8.4f}±{se_Mpt[k]:.4f} {z_tam[k]:+6.1f}")
    print(f"\n  bant-χ²/dof (birincil, 0.86) = {chi:.2f}   "
          f"(yan, 1.10 tam-env) = {chi_tam:.2f}   (eşik 2)")

    hukum = "KAPANDI" if chi <= 2.0 else "GEÇEMEDİ (ölüm dürüst)"
    print(f"\nH-F4 HÜKMÜ (birincil): {hukum}  [χ²/dof = {chi:.2f}]")

    sonuc = {
        "sha_onkayit": ONK["sha256"], "etiket": et,
        "M_olc": M_olc.tolist(), "sigM": sigM.tolist(),
        "M_pred_086": M_pred.tolist(), "se_M_pred": se_Mp.tolist(),
        "M_pred_110tam": M_pred_tam.tolist(), "se_M_pred_tam": se_Mpt.tolist(),
        "z": z.tolist(), "z_tam": z_tam.tolist(),
        "chi2dof_birincil": chi, "chi2dof_tam": chi_tam,
        "m_HA4": m_HA4.tolist(), "m_Hk": m_Hk.tolist(),
        "m_env_kesik": m_env_kes.tolist(), "m_env_tam": m_env_tam.tolist(),
        "HF4_hukum": hukum,
        "muhur_186e_maksfark": muhur_186e,
    }
    json.dump(sonuc, open(S187 / "K4_HA4.json", "w"), indent=1,
              ensure_ascii=False)
    print(f"\n-> {S187/'K4_HA4.json'}  BİTTİ", flush=True)
