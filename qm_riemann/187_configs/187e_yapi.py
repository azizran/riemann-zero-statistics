# -*- coding: utf-8 -*-
"""
187e — K3: YAPI (ι profili, C(τ) rekonstrüksiyonu, biçim KAYITLARI)
===================================================================
ONKAYIT_187 dondurdu:
  ι(τ'_i; b) = −Δm_i,b / m^kesik_b  (BİRİNCİL; katman defterinden)
  yan sütun: koherent katman-ζ_i,b (Σ_i ζ_i = ζ teleskopik)
  C(τ̄_b) = m^kesik_b/m_ölç_b − 1  vs  C_rekon = m^kesik_b/m^ya(1.20)_b − 1
  biçimler: C'ye k·τ^β (185d fit makinesi); ι(HAVUZ)'a üstel ve güç.
  TÜRETİM DEĞİL, KAYIT.

Çıktı: 187/K3_yapi.json + ekran defteri.
"""
import importlib.util
import json
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


def guc_fit(x, y, w):
    """y = k·x^β; β-tarama [0.05,6.0) adım 0.001, kapalı-form k (185d AYNEN)."""
    best = None
    for beta in np.arange(0.05, 6.0, 0.001):
        xb = x ** beta
        k = np.sum(w * y * xb) / np.sum(w * xb * xb)
        chi = np.sum(w * (y - k * xb) ** 2)
        if best is None or chi < best[2]:
            best = (k, beta, chi)
    return best


def ustel_fit(x0, y, w):
    """y = A·exp(−x0/λ); λ-tarama, kapalı-form A."""
    best = None
    for lam in np.arange(0.005, 2.0, 0.0005):
        xb = np.exp(-x0 / lam)
        A = np.sum(w * y * xb) / np.sum(w * xb * xb)
        chi = np.sum(w * (y - A * xb) ** 2)
        if best is None or chi < best[2]:
            best = (A, lam, chi)
    return best


def guc_neg_fit(x0, y, w):
    """y = A·x0^{−p}; p-tarama [0.05,6.0), kapalı-form A."""
    best = None
    for p in np.arange(0.05, 6.0, 0.001):
        xb = x0 ** (-p)
        A = np.sum(w * y * xb) / np.sum(w * xb * xb)
        chi = np.sum(w * (y - A * xb) ** 2)
        if best is None or chi < best[2]:
            best = (A, p, chi)
    return best


if __name__ == "__main__":
    ONK = b187.onkayit()
    print("=" * 78)
    print(f"187e / K3 YAPI  [on-kayit sha {ONK['sha256'][:12]}]")
    print("=" * 78, flush=True)
    izg = ONK["katman_izgara_gercek"]
    NKAT = len(izg) - 1

    G = np.load(S184 / "K1_gercek.npz")
    OZg = np.load(S185 / "OZ_gercek.npz")
    Pg = np.load(S186 / "G1_proj_gercek.npz")
    KP = np.load(S187 / "katmanproj_gercek.npz")
    tau = np.asarray(G["tau"], float)
    aq = np.asarray(G["aq"], float)
    ae = np.asarray(G["aq_eff"], float)
    kenar_b = json.load(open(S186 / "ONKAYIT_186.json"))["kenar"]
    maskeler = [(tau >= kenar_b[b]) & (tau < kenar_b[b + 1]) for b in range(8)]
    maskeler.append((tau >= kenar_b[0]) & (tau < kenar_b[-1]))
    et = [f"{kenar_b[b]:.2f}-{kenar_b[b+1]:.2f}" for b in range(8)] + ["HAVUZ"]
    tau_bar = np.array([np.sum(ae[m] * tau[m]) / ae[m].sum()
                        for m in maskeler])

    def defter(disari):
        """[m_olc, m_kesik, m_ya(k1..k7), iota(1..7), Re zeta_i(1..7),
            Im zeta_i(1..7), C, C_rekon] her maske için."""
        c_olc = b185.c_olculu(G, disari)
        c_oz = b185.c_oz(OZg, aq, disari)
        a_oz = np.abs(c_oz)
        ck = b187.c_kesik(Pg, disari)
        ckat = [b187.c_katman(KP, i, disari) for i in range(NKAT)]
        mix = ck - c_oz
        out = np.zeros((9, 2 + NKAT + NKAT + 2 * NKAT + 2))
        for k, m in enumerate(maskeler):
            sae = ae[m].sum()
            soz = a_oz[m].sum()
            m_olc = (np.abs(c_olc[m]).sum() - soz) / sae
            cya = ck[m].copy()
            m_kes = (np.abs(cya).sum() - soz) / sae
            mya = []
            for i in range(NKAT):
                cya = cya + ckat[i][m]
                mya.append((np.abs(cya).sum() - soz) / sae)
            mm = [m_kes] + mya
            iota = [-(mm[i + 1] - mm[i]) / m_kes for i in range(NKAT)]
            den = np.sum(np.abs(mix[m]) ** 2)
            zi = [np.sum(ckat[i][m] * np.conj(mix[m])) / den
                  for i in range(NKAT)]
            C = m_kes / m_olc - 1.0
            C_rek = m_kes / mya[-1] - 1.0
            out[k] = ([m_olc, m_kes] + mya + iota +
                      [z.real for z in zi] + [z.imag for z in zi] +
                      [C, C_rek])
        return out

    tam = defter(-1)
    reps = np.array([defter(i) for i in range(NJACK)])
    se = b185.jk_se(reps)

    i_mya0, i_iota = 2, 2 + NKAT
    i_zr, i_zi = i_iota + NKAT, i_iota + 2 * NKAT
    i_C, i_Crek = i_zi + NKAT, i_zi + NKAT + 1
    orta = np.array([(izg[i] + izg[i + 1]) / 2 for i in range(NKAT)])
    dtau = np.diff(np.array(izg))

    print("\nι(τ'; bant) İPTAL-YOĞUNLUK PROFİLİ (birincil: −Δm/m^kesik):")
    bas = "  ".join(f"({izg[i]:.2f},{izg[i+1]:.2f}]" for i in range(NKAT))
    print(f"{'bant':>10} | {bas}")
    for k in range(9):
        s = " ".join(f"{tam[k, i_iota+i]:+.4f}" for i in range(NKAT))
        print(f"{et[k]:>10} | {s}")
    print("\nι HAVUZ (±jk se) ve yoğunluk ι/Δτ':")
    for i in range(NKAT):
        v, s_ = tam[8, i_iota + i], se[8, i_iota + i]
        print(f"  ({izg[i]:.2f},{izg[i+1]:.2f}]: ι = {v:+.4f}±{s_:.4f}   "
              f"ι/Δτ' = {v/dtau[i]:+.3f}")

    # koherent katman-ζ (yan sütun) + kapsama
    K2 = json.load(open(S187 / "K2_zeta.json"))
    zg_mod = np.array(K2["gercek"]["zeta_mod"])
    zg_aci = np.array(K2["gercek"]["zeta_aci"])
    zeta_tam = zg_mod * np.exp(1j * np.radians(zg_aci))
    print("\nKOHERENT KATMAN-ζ_i (HAVUZ; yan sütun) ve kapsama:")
    zsum = 0
    for i in range(NKAT):
        z = tam[8, i_zr + i] + 1j * tam[8, i_zi + i]
        zsum += z
        print(f"  ({izg[i]:.2f},{izg[i+1]:.2f}]: |ζ_i| = {abs(z):.4f}  "
              f"açı = {np.degrees(np.angle(z)):+.1f}°")
    kapsama = abs(zsum) / abs(zeta_tam[8])
    print(f"  Σζ_i(≤1.20) = {abs(zsum):.4f} ∠{np.degrees(np.angle(zsum)):+.1f}° "
          f" /  ζ(HAVUZ) = {abs(zeta_tam[8]):.4f} → koherent kapsama = "
          f"{kapsama:.3f}")

    # C(τ) rekonstrüksiyonu
    C, C_rek = tam[:, i_C], tam[:, i_Crek]
    seC, seCr = se[:, i_C], se[:, i_Crek]
    fark_reps = reps[:, :, i_Crek] - reps[:, :, i_C]
    se_fark = b185.jk_se(fark_reps)
    z_fark = (C_rek - C) / se_fark
    print("\nC(τ) REKONSTRÜKSİYONU (C = m^kesik/m_ölç − 1; 186 hedefi "
          "0.46→0.58):")
    print(f"{'bant':>10} {'C ölçülü':>16} {'C_rekon (1.20)':>16} {'z_fark':>7}")
    for k in range(9):
        print(f"{et[k]:>10} {C[k]:7.4f}±{seC[k]:.4f} "
              f"{C_rek[k]:7.4f}±{seCr[k]:.4f} {z_fark[k]:+7.2f}")

    # biçim KAYITLARI (türetim değil)
    x = tau_bar[:8]
    wC = 1.0 / seC[:8] ** 2
    kC, bC, chiC = guc_fit(x, C[:8], wC)
    # jk-se parametreler
    kb_reps = np.array([guc_fit(x, reps[i, :8, i_C],
                                wC)[:2] for i in range(NJACK)])
    se_kC, se_bC = b185.jk_se(kb_reps)
    print(f"\nBİÇİM KAYDI — C(τ) = k·τ^β  (KAYIT; türetim değil):")
    print(f"  k = {kC:.4f}±{se_kC:.4f}   β = {bC:.3f}±{se_bC:.3f}   "
          f"fit-χ²/dof = {chiC/6:.2f}  (dof=6)")

    y_i = tam[8, i_iota:i_iota + NKAT]
    w_i = 1.0 / se[8, i_iota:i_iota + NKAT] ** 2
    x0 = orta - 0.86
    Au, lam, chiu = ustel_fit(x0, y_i, w_i)
    Ap, pw, chip = guc_neg_fit(x0, y_i, w_i)
    print(f"BİÇİM KAYDI — ι(τ') HAVUZ (KAYIT):")
    print(f"  üstel A·exp(−(τ'−0.86)/λ): A = {Au:.4f}  λ = {lam:.4f}  "
          f"fit-χ²/dof = {chiu/5:.2f}  (dof=5)")
    print(f"  güç   A·(τ'−0.86)^(−p)  : A = {Ap:.5f}  p = {pw:.3f}  "
          f"fit-χ²/dof = {chip/5:.2f}  (dof=5)")

    sonuc = {
        "sha_onkayit": ONK["sha256"], "etiket": et, "izgara": izg,
        "tau_bar": tau_bar.tolist(),
        "iota": tam[:, i_iota:i_iota + NKAT].tolist(),
        "se_iota": se[:, i_iota:i_iota + NKAT].tolist(),
        "iota_yogunluk_havuz": (tam[8, i_iota:i_iota + NKAT] /
                                dtau).tolist(),
        "zeta_katman_havuz": [[tam[8, i_zr + i], tam[8, i_zi + i]]
                              for i in range(NKAT)],
        "koherent_kapsama_havuz": float(kapsama),
        "C": C.tolist(), "se_C": seC.tolist(),
        "C_rekon": C_rek.tolist(), "se_C_rekon": seCr.tolist(),
        "z_fark": z_fark.tolist(),
        "bicim_C": {"k": float(kC), "se_k": float(se_kC),
                    "beta": float(bC), "se_beta": float(se_bC),
                    "chi2dof": float(chiC / 6)},
        "bicim_iota_ustel": {"A": float(Au), "lambda": float(lam),
                             "chi2dof": float(chiu / 5)},
        "bicim_iota_guc": {"A": float(Ap), "p": float(pw),
                           "chi2dof": float(chip / 5)},
    }
    json.dump(sonuc, open(S187 / "K3_yapi.json", "w"), indent=1,
              ensure_ascii=False)
    print(f"\n-> {S187/'K3_yapi.json'}  BİTTİ", flush=True)
