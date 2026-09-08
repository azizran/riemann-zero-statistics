# -*- coding: utf-8 -*-
"""
182e — K1(a): KIRPMA HAKEMİNİN n=8 HÂLİ + GAUSS-ALTI HÜKMÜ
===========================================================
ÖN-MÜHÜR: karar kuralı 180a'da (H180a) ve Gauss-altı hükmü 182a'da,
ölçümden ÖNCE donduruldu. Bu betik eşik SEÇMEZ, gaz kümesi SEÇMEZ.

ÜREME KAPISI (R-KAPISI): aşağıdaki `hakem()` VF1..VF4 ile koşulduğunda
`180/K2_HAKEM.json`'un dal, ⟨ρ₃⟩, sd, se, C, z_G, saçılım, eşik
alanlarını TAM SIFIR farkla üretmeli. Üretmezse betik ÖLÜR — yani n=8
tablosu, 180e'nin makinesinin bit-bit aynısıyla basılmıştır.

Kullanım: 182e_hakem8.py
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S182, S180, S169 = SCR / "182", SCR / "180", SCR / "169"
AD = ["100_E", "010_Xa", "001_Xb", "011_XaXb", "110_EXa", "111_hepsi"]


def hakem(gazlar, RHK, RSON, G3, OKH):
    """180e'nin dondurulmuş kuralı — satır satır aynısı."""
    T = {}
    for g in gazlar:
        p = S169 / f"K2b_{g}.json"
        if not p.exists():
            raise SystemExit(f"{p} yok — hüküm verilmez.")
        T[g] = json.load(open(p))["ortalama"]
    r3 = np.array([T[g]["111_hepsi"] for g in gazlar])
    n = len(r3)
    bar = float(r3.mean())
    sd = float(np.std(r3, ddof=1))
    se = sd / np.sqrt(n)
    C = (RHK - bar) / (RHK - G3)
    zG = (bar - G3) / se if se > 0 else float("inf")
    sacilim = float(r3.max() - r3.min())
    esik = 0.5 * abs(RHK - G3)
    if sacilim > esik:
        dal, hukum = "HUKUMSUZ", ("tohum saçılımı çözünürlüğü aşıyor — "
                                  + OKH["HUKUMSUZ"]["hukum"])
    elif C >= 0.70 and abs(zG) <= 3:
        dal, hukum = "A", OKH["DAL_A"]["hukum"]
    elif C <= 0.30:
        dal, hukum = "B", OKH["DAL_B"]["hukum"]
    else:
        dal, hukum = "HUKUMSUZ", OKH["HUKUMSUZ"]["hukum"]
    return dict(gazlar=list(gazlar), n=n, rho3=dict(zip(gazlar, r3.tolist())),
                rho3_bar=bar, sd=sd, se=se, C=C, z_G=zG, sacilim=sacilim,
                esik_sacilim=esik, dal=dal, hukum=hukum, tablo=T)


def main():
    t0 = time.time()
    OK180 = json.load(open(S180 / "ONKAYIT_K0.json"))
    OK182 = json.load(open(S182 / "ONKAYIT_182.json"))
    H = OK180["H180a"]
    RHK, RSON = H["capalar"]["rho3_Hkeskin"], H["capalar"]["rho3_son"]
    G3 = H["capalar"]["GAUSS3"]
    G1 = H["capalar"]["GAUSS1"]
    K180 = json.load(open(S180 / "K2_HAKEM.json"))

    print("=" * 78, flush=True)
    print("182e — KIRPMA HAKEMİ (n=8) + GAUSS-ALTI HÜKMÜ", flush=True)
    print(f"    180 ön-kayıt sha={OK180['sha256'][:16]}…   "
          f"182 ön-kayıt sha={OK182['sha256'][:16]}…", flush=True)
    print("=" * 78, flush=True)

    # ----------------------------------------------------- R-KAPISI
    R4 = hakem([f"VF{i}" for i in (1, 2, 3, 4)], RHK, RSON, G3, H)
    kotu = []
    for a in ("dal", "rho3_bar", "sd", "se", "C", "z_G", "sacilim",
              "esik_sacilim"):
        if R4[a] != K180[a]:
            kotu.append((a, R4[a], K180[a]))
    print(f"  R-KAPISI (VF1..VF4 ile 180/K2_HAKEM.json'u üretme): "
          f"8 alan, {len(kotu)} fark", flush=True)
    if kotu:
        for c in kotu:
            print("   ✗", c, flush=True)
        raise SystemExit("R-KAPISI TUTMADI — n=8 hükmü basılmaz.")
    print("  R-KAPISI ✓ TAM SIFIR FARK (makine 180e'ninkiyle birebir)",
          flush=True)

    # -------------------------------------------------------- n = 8
    gz = [f"VF{i}" for i in range(1, 9)]
    R8 = hakem(gz, RHK, RSON, G3, H)
    T = R8["tablo"]
    print("\n  T5 — 3-BACAK KIRPMA AKTARIMI (169_k2b, lo∈[0.52,0.68] ort.)")
    print("  gaz       " + " ".join(f"{a.split('_')[1]:>8s}" for a in AD))
    ek = {}
    for g in ("Hkeskin", "son") + tuple(gz):
        p = S169 / f"K2b_{g}.json"
        if not p.exists():
            continue
        v = json.load(open(p))["ortalama"]
        ek[g] = v
        print(f"  {g:9s} " + " ".join(f"{v[a]:8.4f}" for a in AD), flush=True)
    print("  GAUSS     " + " ".join(
        f"{G1 ** a.split('_')[0].count('1'):8.4f}" for a in AD), flush=True)

    print(f"\n  ⟨ρ₃⟩_VF = {R8['rho3_bar']:.6f}  sd = {R8['sd']:.6f}  "
          f"se = {R8['se']:.6f}  (n = {R8['n']})")
    print("   tek tek: " + "  ".join(f"{g}={v:.4f}"
                                     for g, v in R8["rho3"].items()))
    print(f"  C = {R8['C']:+.4f}    z_G = {R8['z_G']:+.3f}")
    print(f"  saçılım = {R8['sacilim']:.6f}   (ön-kayıtlı HÜKÜMSÜZ eşiği "
          f"{R8['esik_sacilim']:.6f})")
    print(f"\n  >>> ÖN-KAYITLI DAL (n=8): {R8['dal']}")
    print(f"  >>> {R8['hukum']}", flush=True)

    # -------- ÖNGÖRÜ P1'in denetimi (range monotonluğu) --------------
    p1 = bool(R8["sacilim"] >= K180["sacilim"])
    print(f"\n  ÖNGÖRÜ P1 denetimi: saçılım(n=8) = {R8['sacilim']:.6f} "
          f"≥ saçılım(n=4) = {K180['sacilim']:.6f} ?  "
          f"{'✓ TUTTU' if p1 else '✗ TUTMADI (imkânsız — makine bozuk)'}",
          flush=True)

    # ------------------------ GAUSS-ALTI HÜKMÜ (182a'da dondurulmuş) --
    r3 = np.array(list(R8["rho3"].values()))
    hepsi_alt = bool(np.all(r3 < G3))
    z_ok = bool(R8["z_G"] <= -5.0)
    ga = "GAUSS-ALTI (n=8)" if (hepsi_alt and z_ok) else "HÜKÜM YOK"
    print("\n" + "=" * 78)
    print("  GAUSS-ALTI HÜKMÜ (182a/KARAR.A — z_G ≤ −5 VE 8/8 tohum "
          "ρ₃ < GAUSS3)")
    print("=" * 78)
    print(f"    GAUSS3 = {G3:.6f}   ⟨ρ₃⟩ = {R8['rho3_bar']:.6f}   "
          f"açık = {100*(R8['rho3_bar']/G3-1):+.2f}%")
    print(f"    z_G = {R8['z_G']:+.3f}  (eşik ≤ −5)  ⇒ "
          f"{'✓' if z_ok else '✗'}")
    print(f"    tek tek ρ₃ < GAUSS3: {int((r3<G3).sum())}/8  ⇒ "
          f"{'✓' if hepsi_alt else '✗'}")
    print(f"    >>> {ga}", flush=True)

    # -------- bacak-sayısı ortalamaları (kayıt; karar dışı) -----------
    NB = {"100_E": 1, "010_Xa": 1, "001_Xb": 1, "011_XaXb": 2,
          "110_EXa": 2, "111_hepsi": 3}
    bo = {g: {k: float(np.mean([v[a] for a in AD if NB[a] == k]))
              for k in (1, 2, 3)} for g, v in ek.items()}
    print("\n  bacak-sayısı ortalamaları (kayıt; karar YALNIZ ρ₃'e bağlı)")
    for k in (1, 2, 3):
        gv = G1 ** k
        print(f"   {k} bacak (Gauss {gv:.4f}): " + "  ".join(
            f"{g}={bo[g][k]:.4f}({100*(bo[g][k]/gv-1):+.1f}%)"
            for g in bo), flush=True)

    # -------- E ve X bacaklarının n=8 istatistiği (kayıt) -------------
    print("\n  bacak bacak, VF ailesi n=8 (ort ± sd; Gauss 0.7979)")
    for a in AD:
        v = np.array([ek[g][a] for g in gz])
        print(f"   ρ({a.split('_')[1]:7s}) = {v.mean():.4f} ± "
              f"{v.std(ddof=1):.4f}   Gauss {G1**a.split('_')[0].count('1'):.4f}"
              f"   fark {100*(v.mean()/G1**a.split('_')[0].count('1')-1):+6.2f}%",
              flush=True)

    rec = dict(R4=R4, R8=R8, R_kapisi_fark=len(kotu), P1_tuttu=p1,
               gauss_alti=ga, hepsi_alt=hepsi_alt, z_ok=z_ok,
               bacak_ort=bo,
               bacak_n8={a: dict(
                   ort=float(np.mean([ek[g][a] for g in gz])),
                   sd=float(np.std([ek[g][a] for g in gz], ddof=1)))
                   for a in AD},
               onkayit_182_sha=OK182["sha256"],
               onkayit_180_sha=OK180["sha256"], zaman=time.ctime(),
               sure_s=time.time() - t0)
    p = S182 / "K2_HAKEM_n8.json"
    p.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=float))
    print(f"\n-> {p}", flush=True)


if __name__ == "__main__":
    main()
