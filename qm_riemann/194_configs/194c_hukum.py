# -*- coding: utf-8 -*-
"""
194c — K2: HÜKÜM (H-194a..e; ONKAYIT_194 AYNEN, kurtarma yok)
================================================================
K1_bos_194.json'daki boş/μ=0 pencere ölçümlerini (194b) havuzlar (pooled
8-blok jackknife — her pencere AYNI 8 fiziksel bloktan geldiği için leave-
one-block-out replikaları önce pencereler arasında ORTALANIR, sonra jk()
uygulanır) ve ön-kayıtlı eşiklerle (ONKAYIT_194) karşılaştırır:
 H-194a: κ̄_boş bandı + anlamlılık (ölüm: |κ̄_boş|<2σ ya da pozitif).
 H-194b: Δω eğilimi (KAYIT; havuzlanmış eğim, jackknife se).
 H-194c: boş pencerelerde mod 3,5,10 sınıf payları — havuzlanmış s_r,
         1/φ(a) ± [0.10·(1/φ(a)) + 2σ] bandı (eşiksiz, ama rapor edilir).
 H-194d: κ̄_μ0 vs κ̄_boş (ölüm: μ0 ortalaması ≥3σ daha derin/negatif).
 H-194e: 193'ün +log10/+log7/+log5 paylarının κ̄_boş ile düzeltilmesi;
         s_r' = (κ_r − κ̄_boş/φ(a)) / (κ_top − κ̄_boş), delta-yöntemi se
         yayılımı (κ_r,κ_top se'leri 193'ten, κ̄_boş se'si 194'ten;
         bağımsız varsayılır — kovaryans hesaba katılmaz, KALEM bunu istemiyor).
         Ölüm: +log10'da herhangi bir düzeltilmiş pay ±0.15 dışı.
Çıktı: 194/HUKUM_194.json + ekran.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S193, S194 = SCR / "193", SCR / "194"


def onkayit():
    o = json.load(open(S194 / "ONKAYIT_194.json"))
    s = hashlib.sha256((QM / "194_configs" / "194a_onkayit.py").read_bytes()).hexdigest()
    if s != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {s} != {o['sha256']}")
    return o


def jk(r):
    r = np.asarray(r, float)
    n = len(r)
    return float(np.sqrt((n - 1) / n * np.sum((r - r.mean()) ** 2)))


def havuzla(reps_list, point_list):
    """Pencereler arası havuzlama: pooled_reps[j] = mean_i(reps_i[j]); se = jk(pooled_reps);
    point = mean_i(point_i)."""
    reps = np.array(reps_list, float)          # (n_pencere, 8)
    pooled_reps = reps.mean(axis=0)             # (8,)
    se = jk(pooled_reps)
    point = float(np.mean(point_list))
    return point, se, pooled_reps


if __name__ == "__main__":
    ONK = onkayit()
    K1 = json.load(open(S194 / "K1_bos_194.json"))
    assert K1["sha_onkayit"] == ONK["sha256"]
    BOS = K1["bos_pencereler"]
    MU0 = K1["mu0_pencereler"]
    PHI = ONK["phi_tablosu"]
    print("=" * 78)
    print(f"194c / K2 HÜKÜM  [on-kayit {ONK['zaman']} sha {ONK['sha256'][:12]}]")
    print("=" * 78)
    out = {"sha_onkayit": ONK["sha256"]}

    # ================= κ̄_boş, κ̄_μ0 (havuzlanmış) =================
    kappa_bar_bos, se_bos, reps_bos = havuzla(
        [w["kappa_reps"] for w in BOS.values()], [w["kappa"] for w in BOS.values()])
    kappa_bar_mu0, se_mu0, reps_mu0 = havuzla(
        [w["kappa_reps"] for w in MU0.values()], [w["kappa"] for w in MU0.values()])
    out["kappa_bar_bos"] = kappa_bar_bos
    out["se_bos"] = se_bos
    out["n_bos"] = len(BOS)
    out["kappa_bar_mu0"] = kappa_bar_mu0
    out["se_mu0"] = se_mu0
    out["n_mu0"] = len(MU0)
    print(f"\nκ̄_boş = {kappa_bar_bos:+.5f} ± {se_bos:.5f}  (N={len(BOS)} pencere)")
    for ad, w in BOS.items():
        print(f"   {ad}: merkez={w['merkez']:.5f}  κ={w['kappa']:+.5f}±{w['se_kappa']:.5f}  "
              f"n_cizgi={w['n_cizgi']}")
    print(f"\nκ̄_μ0  = {kappa_bar_mu0:+.5f} ± {se_mu0:.5f}  (N={len(MU0)} pencere)")
    for ad, w in MU0.items():
        print(f"   {ad}: merkez={w['merkez']:.5f}  κ={w['kappa']:+.5f}±{w['se_kappa']:.5f}  "
              f"n_cizgi={w['n_cizgi']}")

    # ================= H-194a =================
    bant_lo, bant_hi = ONK["H_194a"]["bant"]
    bantta = bant_lo <= kappa_bar_bos <= bant_hi
    anlamli_3s = kappa_bar_bos <= -3 * se_bos
    olum_a = (abs(kappa_bar_bos) < 2 * se_bos) or (kappa_bar_bos > 0)
    if olum_a:
        hA = "ÖLDÜ"
    elif anlamli_3s and bantta:
        hA = "MÜHÜR"
    else:
        eksik = []
        if not bantta:
            eksik.append(f"bant [{bant_lo},{bant_hi}] dışı")
        if not anlamli_3s:
            eksik.append("<3σ")
        hA = "KAYIT (kısmi: " + ", ".join(eksik) + ")"
    out["H_194a"] = {"hukum": hA, "kappa_bar_bos": kappa_bar_bos, "se_bos": se_bos,
                     "sigma_negatif": kappa_bar_bos / se_bos if se_bos else None,
                     "bantta": bantta, "anlamli_3sigma": anlamli_3s, "olum": olum_a}
    print(f"\nH-194a (EKSİ TABAN VAR): κ̄_boş={kappa_bar_bos:+.5f}±{se_bos:.5f}  "
          f"({kappa_bar_bos/se_bos:+.2f}σ)  bant[{bant_lo},{bant_hi}]={'EVET' if bantta else 'HAYIR'}  "
          f"->  {hA}")

    # ================= H-194b (eğilim; KAYIT) =================
    x = np.array([w["merkez"] for w in BOS.values()])
    y = np.array([w["kappa"] for w in BOS.values()])
    reps_mat = np.array([w["kappa_reps"] for w in BOS.values()])   # (n,8)
    slope_point = float(np.polyfit(x, y, 1)[0])
    slopes_j = [float(np.polyfit(x, reps_mat[:, j], 1)[0]) for j in range(reps_mat.shape[1])]
    se_slope = jk(slopes_j)
    sabit = abs(slope_point) < 2 * se_slope
    hB = ("SABİT (|eğim|<2σ)" if sabit else
          f"EĞİM VAR (|eğim|>=2σ; biçim: {'artan' if slope_point>0 else 'azalan'})")
    out["H_194b"] = {"hukum": hB, "egim": slope_point, "se_egim": se_slope, "sabit": bool(sabit)}
    print(f"\nH-194b (SABİTLİK, KAYIT): eğim={slope_point:+.6f}±{se_slope:.6f}  "
          f"({slope_point/se_slope:+.2f}σ)  ->  {hB}")

    # ================= H-194c (sınıftan bağımsızlık; boş pencerelerde mod 3,5,10) ==
    MODULER = {3: [1, 2], 5: [1, 2, 3, 4], 10: [1, 3, 7, 9]}
    h194c_out = {}
    hepsi_bantta_genel = True
    for a, rlist in MODULER.items():
        phi_a = PHI[str(a)]
        pred = 1.0 / phi_a
        tol = 0.10 * pred
        satirlar = {}
        for r in rlist:
            reps_r = np.array([w[f"mod{a}"]["siniflar"][f"r{r}"]["s_reps"] for w in BOS.values()])
            s_point = float(np.mean([w[f"mod{a}"]["siniflar"][f"r{r}"]["s"] for w in BOS.values()]))
            pooled = reps_r.mean(axis=0)
            se_pooled = jk(pooled)
            bant_toplam = tol + 2 * se_pooled
            b_ = abs(s_point - pred) <= bant_toplam
            hepsi_bantta_genel &= b_
            satirlar[f"r{r}"] = {"s_havuz": s_point, "se_s_havuz": se_pooled,
                                 "pred_esit": pred, "bant_toplam": bant_toplam, "bantta": bool(b_)}
        h194c_out[f"mod{a}"] = {"phi": phi_a, "pred_1_phi": pred, "siniflar": satirlar,
                                "hepsi_bantta": all(v["bantta"] for v in satirlar.values())}
    out["H_194c"] = {"hukum": ("SINIFTAN BAĞIMSIZ (cos deseni YOK)" if hepsi_bantta_genel
                               else "KISMEN BAĞIMSIZ (bazı sınıflar bant dışı)"),
                     "moduller": h194c_out}
    print(f"\nH-194c (SINIFTAN BAĞIMSIZ, eşiksiz KAYIT): {out['H_194c']['hukum']}")
    for modad, mv in h194c_out.items():
        print(f"   {modad} (φ={mv['phi']}, öngörü 1/φ={mv['pred_1_phi']:.4f}):")
        for r, rv in mv["siniflar"].items():
            print(f"      {r}: s_havuz={rv['s_havuz']:+.4f}±{rv['se_s_havuz']:.4f}  "
                  f"bant±{rv['bant_toplam']:.4f}  {'EVET' if rv['bantta'] else 'HAYIR'}")

    # ================= H-194d =================
    sigma_comb = float(np.sqrt(se_bos ** 2 + se_mu0 ** 2))
    fark_d = kappa_bar_mu0 - kappa_bar_bos
    uyumlu_2s = abs(fark_d) <= 2 * sigma_comb
    olum_d = fark_d <= -3 * sigma_comb
    if olum_d:
        hD = "ÖLDÜ"
    elif uyumlu_2s:
        hD = "MÜHÜR"
    else:
        hD = "KAYIT (kısmi: 2σ dışı, <3σ ölüm eşiği)"
    out["H_194d"] = {"hukum": hD, "fark": fark_d, "sigma_comb": sigma_comb,
                     "fark_sigma": fark_d / sigma_comb if sigma_comb else None,
                     "uyumlu_2sigma": uyumlu_2s, "olum": olum_d}
    print(f"\nH-194d (ÇUKUR=TABAN): κ̄_μ0−κ̄_boş={fark_d:+.5f}  σ_comb={sigma_comb:.5f}  "
          f"({fark_d/sigma_comb:+.2f}σ)  ->  {hD}")

    # ================= H-194e (193 paylarının κ̄_boş ile düzeltilmesi) =============
    K193 = json.load(open(S193 / "K1_sinif.json"))
    HEDEF_E = {"+log10": {"a": 10, "olum_bant": 0.15}, "+log7": {"a": 7, "olum_bant": None},
              "+log5": {"a": 5, "olum_bant": None}}
    e_out = {}
    olum_e_liste = []
    hepsi_010 = True
    for ad, meta in HEDEF_E.items():
        h = K193["hedefler"][ad]
        a = meta["a"]
        phi_a = PHI[str(a)]
        kappa_top = h["kappa_top"]
        se_top = h["se_kappa_top"]
        satirlar = {}
        for r_key, c in h["siniflar"].items():
            kappa_r = c["kappa"]
            se_r = c["se_kappa"]
            s_pred = c["s_pred"]
            D = kappa_top - kappa_bar_bos
            Nn = kappa_r - kappa_bar_bos / phi_a
            s_corr = Nn / D
            d_dkr = 1.0 / D
            d_dktop = -s_corr / D
            d_dkbos = (-1.0 / phi_a) / D + s_corr / D
            se_corr = float(np.sqrt((d_dkr * se_r) ** 2 + (d_dktop * se_top) ** 2 +
                                    (d_dkbos * se_bos) ** 2))
            fark = s_corr - s_pred
            bantta10 = abs(fark) <= 0.10
            hepsi_010 &= bantta10
            if ad == "+log10" and abs(fark) > meta["olum_bant"]:
                olum_e_liste.append(f"{ad} {r_key}: |Δ|={abs(fark):.4f} > 0.15")
            satirlar[r_key] = {"s_ham": c["s"], "se_s_ham": c["se_s"], "s_pred": s_pred,
                               "s_duzeltilmis": float(s_corr), "se_duzeltilmis": se_corr,
                               "fark_pred": float(fark), "bantta_010": bool(bantta10)}
        e_out[ad] = {"kappa_top": kappa_top, "se_kappa_top": se_top, "siniflar": satirlar}
    if olum_e_liste:
        hE = "ÖLDÜ"
    elif hepsi_010:
        hE = "MÜHÜR"
    else:
        hE = "KAYIT (kısmi: bazı paylar ±0.10 dışı ama ölüm eşiği ±0.15 aşılmadı)"
    out["H_194e"] = {"hukum": hE, "olum_liste": olum_e_liste, "hedefler": e_out}
    print(f"\nH-194e (GERME KAPANIŞI): κ̄_boş={kappa_bar_bos:+.5f}±{se_bos:.5f} ile düzeltme")
    for ad, ev in e_out.items():
        print(f"   {ad}: κ_top(193)={ev['kappa_top']:+.5f}±{ev['se_kappa_top']:.5f}")
        for r, rv in ev["siniflar"].items():
            print(f"      {r}: ham={rv['s_ham']:+.4f}  düzeltilmiş={rv['s_duzeltilmis']:+.4f}"
                  f"±{rv['se_duzeltilmis']:.4f}  öngörü={rv['s_pred']:+.4f}  "
                  f"Δ={rv['fark_pred']:+.4f}  ±0.10 {'EVET' if rv['bantta_010'] else 'HAYIR'}")
    print(f"   ölüm listesi (+log10, ±0.15): {olum_e_liste if olum_e_liste else 'yok'}")
    print(f"   ->  {hE}")

    out["hukum_ozet"] = {"H-194a": hA, "H-194b": hB, "H-194c": out["H_194c"]["hukum"],
                         "H-194d": hD, "H-194e": hE}
    print("\n" + "=" * 78)
    print("HÜKÜMLER:", out["hukum_ozet"])
    print("=" * 78)
    json.dump(out, open(S194 / "HUKUM_194.json", "w"), indent=1, ensure_ascii=False)
    print(f"-> {S194/'HUKUM_194.json'}  BİTTİ")
