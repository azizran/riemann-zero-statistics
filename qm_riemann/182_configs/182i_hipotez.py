# -*- coding: utf-8 -*-
"""
182i — K2 (HÜKÜM): H-G1 / H-G2 / H-G3 YARIŞI
=============================================
ÖN-MÜHÜR: ölüm ve yaşama koşulları `182/ONKAYIT_182.json` → `hipotezler`
içinde, 182g koşmadan ÖNCE donduruldu. Bu betik SADECE o koşulları
diskteki `182/ANATOMI_VF*.json` çıktılarına uygular; eşik, gaz kümesi
ya da istatistik SEÇMEZ. Ölümler kurtarmasızdır.

Kullanım: 182i_hipotez.py
"""
import json
import time
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S182 = SCR / "182"


def main():
    t0 = time.time()
    OK = json.load(open(S182 / "ONKAYIT_182.json"))
    G1 = OK["capalar"]["GAUSS1"]
    G3 = OK["capalar"]["GAUSS3"]
    gz = [f"VF{i}" for i in range(1, 9)]
    A = {}
    for g in gz:
        p = S182 / f"ANATOMI_{g}.json"
        if not p.exists():
            raise SystemExit(f"{p} yok — hüküm verilmez.")
        A[g] = json.load(open(p))
    print("=" * 78, flush=True)
    print("182i — H-G1 / H-G2 / H-G3 YARIŞI (ön-kayıtlı ölüm koşulları)",
          flush=True)
    print(f"    182 ön-kayıt sha = {OK['sha256'][:16]}…", flush=True)
    print("=" * 78, flush=True)

    def col(f):
        return np.array([f(A[g]) for g in gz], float)

    kE = col(lambda d: d["kappa"]["E"])
    kXa = col(lambda d: d["kappa"]["Xa"])
    kXb = col(lambda d: d["kappa"]["Xb"])
    uE = col(lambda d: d["kurt"]["E"])
    rho = col(lambda d: d["rho_disk"]["100_E"])
    rho3 = col(lambda d: d["rho_disk"]["111_hepsi"])
    rJ = col(lambda d: d["rho_J"]["E"])
    rJ3 = col(lambda d: d["rho_J"]["hepsi"])
    ro = col(lambda d: d["rho_sur"]["E"])
    ro3 = col(lambda d: d["rho_sur"]["hepsi"])
    rm = col(lambda d: d["rho_Esur"])
    Rec = col(lambda d: d["Rec"])
    rXX = col(lambda d: d["korr_XaXb"])
    rXXs = col(lambda d: d["korr_XaXb_sur"])
    bE = col(lambda d: d["beta_E"])
    pC = col(lambda d: d["pay_cizgi"])
    sE = col(lambda d: d["egim"]["E"])
    sXa = col(lambda d: d["egim"]["Xa"])
    sXb = col(lambda d: d["egim"]["Xb"])

    def oz(v):
        return float(v.mean()), float(np.std(v, ddof=1)), \
            float(np.std(v, ddof=1) / np.sqrt(len(v)))

    print("\n  T-A1 — ÇİZGİ / BANT AYRIŞIMI (n=8; ρ(E) ≡ κ_E · β_E)")
    print("  gaz     κ_E       κ_Xa      κ_Xb      kurt_E    ρ(E)      "
          "β_E       PAY_ÇİZGİ")
    for i, g in enumerate(gz):
        print(f"  {g:6s} {kE[i]:.5f}  {kXa[i]:.5f}  {kXb[i]:.5f}  "
              f"{uE[i]:+.5f}  {rho[i]:.5f}  {bE[i]:.5f}  {pC[i]:+.4f}",
              flush=True)
    for nm, v in (("κ_E", kE), ("ρ(E)", rho), ("β_E", bE),
                  ("PAY_ÇİZGİ", pC), ("kurt_E", uE)):
        m, s, e = oz(v)
        print(f"    ⟨{nm}⟩ = {m:+.6f} ± {s:.6f} (sd)  se {e:.6f}")
    print(f"    Gauss: κ = {G1:.5f}, β = 1, PAY_ÇİZGİ = 1 olurdu "
          f"(açığın tamamı çizgide)")

    print("\n  T-A2 — τ EĞİMLERİ (5 bant, doğrusal)")
    for nm, v in (("s_E", sE), ("s_Xa", sXa), ("s_Xb", sXb)):
        m, s, e = oz(v)
        print(f"    ⟨{nm}⟩ = {m:+.5f} ± {s:.5f}  se {e:.5f}")
    sX = np.concatenate([sXa, sXb])
    m_sE, _, se_sE = oz(sE)
    egim_bant = bool(abs(m_sE) >= float(np.mean(np.abs(sX))) + 3 * se_sE)
    _eb = "✓ BANT karakterli" if egim_bant else "✗ gösterilemedi"
    print(f"    ÖN-KAYIT: |⟨s_E⟩| ≥ ⟨|s_X|⟩ + 3·se(s_E) ?  "
          f"{abs(m_sE):.5f} ≥ {float(np.mean(np.abs(sX))):.5f} + "
          f"{3*se_sE:.5f} = {float(np.mean(np.abs(sX)))+3*se_sE:.5f}  ⇒ {_eb}")

    print("\n  T-A3 — SURROGATE'LAR ve BAĞIMSIZ KESTİRİMCİ (n=8)")
    print("  gaz     ρ(E)disk  ρ_J(E)    |ρ_J/ρ−1|  ρ°(E)     ρ*(E)     "
          "Rec       ρ°(3)     ρ_J(3)")
    for i, g in enumerate(gz):
        print(f"  {g:6s} {rho[i]:.5f}  {rJ[i]:.5f}  {abs(rJ[i]/rho[i]-1):.5f}"
              f"    {ro[i]:.5f}  {rm[i]:.5f}  {Rec[i]:+.4f}  {ro3[i]:.5f}  "
              f"{rJ3[i]:.5f}", flush=True)

    print("\n  T-A4 — İNŞA ÖZ-TUTARLILIĞI (korr(Xa,Xb), R_η)")
    print("  gaz     korr(Xa,Xb)  korr(Xa°,Xb°)   R_η       1−R_η     "
          "Kov(artık,η)/Var   Var(çiz)/P")
    Rl, kv = [], []
    for i, g in enumerate(gz):
        d = A[g]["R_eta"]
        if d is None:
            print(f"  {g:6s} {rXX[i]:+.5f}     {rXXs[i]:+.5f}        "
                  f"(R_η ölçülmedi)")
            continue
        Rl.append(d["R"])
        kv.append(d["kov_art_x"] / d["Var"])
        print(f"  {g:6s} {rXX[i]:+.5f}     {rXXs[i]:+.5f}      "
              f"{d['R']:.5f}  {1-d['R']:+.5f}   {d['kov_art_x']/d['Var']:+.5f}"
              f"          {d['Var_cizgi']/d['P']:.5f}", flush=True)
    Rl, kv = np.array(Rl), np.array(kv)

    # =================================================================
    # ÖN-KAYITLI HÜKÜMLER
    # =================================================================
    print("\n" + "=" * 78)
    print("  ÖN-KAYITLI HÜKÜMLER (182a; kurtarma yok)")
    print("=" * 78)

    # ---- H-G1 -------------------------------------------------------
    kE_sap = np.abs(kE / G1 - 1)
    acik = (G1 - rho) / G1
    c1 = bool(np.all(kE_sap <= 0.05) and np.all(acik >= 0.10))
    ro_sap = abs(float(ro.mean()) / G1 - 1)
    c2 = bool(ro_sap <= 0.02)
    kE_acik = np.abs(kE - G1)
    yarim = 0.5 * (G1 - rho)
    yasa1 = bool(ro_sap > 0.05 or np.any(kE_acik >= yarim))
    hg1 = ("ÖLDÜ" if (c1 and c2) else
           ("YAŞIYOR" if yasa1 else "HÜKÜMSÜZ"))
    print(f"\n  H-G1 (Gauss referansının türetimi)")
    print(f"    (i)  maks|κ_E/√(2/π)−1| = {kE_sap.max():.5f} ≤ 0.05 ? "
          f"{'✓' if np.all(kE_sap<=0.05) else '✗'}   "
          f"min açık (√(2/π)−ρ(E))/√(2/π) = {acik.min():.5f} ≥ 0.10 ? "
          f"{'✓' if np.all(acik>=0.10) else '✗'}")
    print(f"    (ii) ⟨ρ°(E)⟩ = {ro.mean():.5f}   |⟨ρ°(E)⟩/√(2/π) − 1| = "
          f"{ro_sap:.5f} ≤ 0.02 ? {'✓' if c2 else '✗'}")
    print(f"    κ_E'nin açığın yarısını taşıması: maks|κ_E−√(2/π)| = "
          f"{kE_acik.max():.5f} vs min yarım açık {yarim.min():.5f} ⇒ "
          f"{'taşıyor' if np.any(kE_acik>=yarim) else 'TAŞIMIYOR'}")
    print(f"    >>> H-G1: {hg1}")

    # ---- H-G3 -------------------------------------------------------
    d3 = np.abs(rJ / rho - 1)
    acikJ = (G1 - rJ) / G1
    o3 = bool(np.all(d3 <= 0.05) and np.all(acikJ >= 0.10))
    y3 = bool(float(d3.mean()) > 0.10)
    hg3 = "ÖLDÜ" if o3 else ("YAŞIYOR" if y3 else "HÜKÜMSÜZ")
    print(f"\n  H-G3 (kestirimci artefaktı)")
    print(f"    maks|ρ_J(E)/ρ(E) − 1| = {d3.max():.5f} ≤ 0.05 ? "
          f"{'✓' if np.all(d3<=0.05) else '✗'}   (ort {d3.mean():.5f})")
    print(f"    bağımsız kestirimcinin Gauss açığı: min = {acikJ.min():.5f} "
          f"≥ 0.10 ? {'✓' if np.all(acikJ>=0.10) else '✗'}")
    print(f"    >>> H-G3: {hg3}")

    # ---- H-G2 -------------------------------------------------------
    m_r, sd_r, se_r = oz(rXX)
    isaret_tek = bool(np.all(np.sign(rXX) == np.sign(rXX[0])))
    sifir = bool(abs(m_r) <= 2 * se_r)
    Rec_dus = bool(float(Rec.mean()) <= 0.30)
    Rec_yuk = int((Rec >= 0.70).sum())
    o2 = bool(sifir or (not isaret_tek) or Rec_dus)
    y2 = bool((abs(m_r) >= 5 * se_r) and isaret_tek and Rec_yuk >= 7)
    hg2 = "ÖLDÜ" if o2 else ("YAŞIYOR" if y2 else "HÜKÜMSÜZ")
    print(f"\n  H-G2 (inşa öz-tutarlılığı)")
    print(f"    ⟨korr(Xa,Xb)⟩ = {m_r:+.5f} ± {sd_r:.5f} (sd), se = "
          f"{se_r:.5f}  ⇒ |⟨r⟩|/se = {abs(m_r)/se_r:.1f}")
    print(f"    işaret 8/8 tekdüze ? {'✓' if isaret_tek else '✗'}   "
          f"sıfırla uyumlu (|⟨r⟩| ≤ 2se) ? {'EVET' if sifir else 'HAYIR'}")
    print(f"    surrogate kontrolü ⟨korr(Xa°,Xb°)⟩ = {rXXs.mean():+.5f} ± "
          f"{np.std(rXXs, ddof=1):.5f}")
    print(f"    ⟨Rec⟩ = {Rec.mean():+.4f}   Rec ≥ 0.70 olan tohum: "
          f"{Rec_yuk}/8   Rec ≤ 0.30 ? {'EVET' if Rec_dus else 'HAYIR'}")
    if len(Rl):
        print(f"    R_η: ⟨R_η⟩ = {Rl.mean():.5f} ± {np.std(Rl, ddof=1):.5f}"
              f"  (n={len(Rl)});  R_η < 1 olan: {int((Rl<1).sum())}/{len(Rl)}"
              f"  ⇒ Kov(artık,η) > 0 ADRESİ: "
              f"⟨Kov/Var⟩ = {kv.mean():+.5f}")
    print(f"    >>> H-G2: {hg2}")

    rec = dict(
        gazlar=gz, GAUSS1=G1, GAUSS3=G3,
        kappa_E=kE.tolist(), kappa_Xa=kXa.tolist(), kappa_Xb=kXb.tolist(),
        kurt_E=uE.tolist(), rho_E=rho.tolist(), rho3=rho3.tolist(),
        beta_E=bE.tolist(), pay_cizgi=pC.tolist(),
        rho_J_E=rJ.tolist(), rho_J_3=rJ3.tolist(),
        rho_sur_E=ro.tolist(), rho_sur_3=ro3.tolist(),
        rho_Esur=rm.tolist(), Rec=Rec.tolist(),
        korr_XaXb=rXX.tolist(), korr_XaXb_sur=rXXs.tolist(),
        egim_E=sE.tolist(), egim_Xa=sXa.tolist(), egim_Xb=sXb.tolist(),
        egim_bant_karakteri=egim_bant,
        R_eta=Rl.tolist(), kov_art_over_var=kv.tolist(),
        hukum=dict(H_G1=hg1, H_G2=hg2, H_G3=hg3),
        H_G1_kanit=dict(kE_sapma_max=float(kE_sap.max()),
                        acik_min=float(acik.min()),
                        rho_sur_ort=float(ro.mean()), rho_sur_sapma=ro_sap),
        H_G2_kanit=dict(korr_ort=m_r, korr_sd=sd_r, korr_se=se_r,
                        isaret_tekduze=isaret_tek, Rec_ort=float(Rec.mean()),
                        Rec_gecen=Rec_yuk),
        H_G3_kanit=dict(fark_max=float(d3.max()), fark_ort=float(d3.mean()),
                        acikJ_min=float(acikJ.min())),
        onkayit_182_sha=OK["sha256"], zaman=time.ctime(),
        sure_s=time.time() - t0)
    p = S182 / "K2_HIPOTEZ.json"
    p.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=float))
    print(f"\n-> {p}   ({time.time()-t0:.1f} s)", flush=True)


if __name__ == "__main__":
    main()
