# -*- coding: utf-8 -*-
"""
192d — K3: A2 PARLAKLIK ORANLARI + HÜKÜM (H-192a/b/c; ONKAYIT_192 AYNEN, kurtarma yok)
=====================================================================================
 A2  düşük havuz profili (192b profiller_192.npz): P(h) = Σ κ(j), |c_j − h| ≤ 0.0625
     (5 dilim, |B_j| ≥ 4); oranlar + jk se (replika oranı); bant [öngörü/2, öngörü·2].
     KAYIT: taban-düzeltmeli P' = Σ(κ − taban_medyan) (192b tabanı), tepe oranı
     κ_tepe/κ_tepe(+log2), SON penceresi aynı oranlar (yoğunluk × 0.025).
 H-192a  192b K1_katalog.json (düşük) → MÜHÜR / ÖLDÜ / KAYIT.
 H-192b  192c K2_sinif.json → MÜHÜR / ÖLDÜ / KAYIT (+ yan öngörüler KAYIT).
 H-192c  KAYIT (bant sayımı).
Çıktı: 192/HUKUM_192.json + ekran.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S192 = SCR / "192"
NJACK = 8
EPS = 1e-9


def jk(r):
    r = np.asarray(r, float)
    return float(np.sqrt((NJACK - 1) / NJACK * np.sum((r - r.mean()) ** 2)))


def onkayit():
    o = json.load(open(S192 / "ONKAYIT_192.json"))
    s = hashlib.sha256((QM / "192_configs" / "192a_onkayit.py").read_bytes()).hexdigest()
    if s != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {s} != {o['sha256']}")
    return o


if __name__ == "__main__":
    ONK = onkayit()
    KAT = ONK["B"]["katalog"]
    B = json.load(open(S192 / "K1_katalog.json"))
    A = json.load(open(S192 / "K2_sinif.json"))
    PF = np.load(S192 / "profiller_192.npz")
    J = PF["J"]
    jmin = int(J[0])
    print("=" * 78)
    print(f"192d / K3 A2 + HÜKÜM  [on-kayit {ONK['zaman']} sha {ONK['sha256'][:12]}]")
    print("=" * 78)
    out = {"sha_onkayit": ONK["sha256"]}

    # ======================= A2 =======================
    def parlaklik(ad, kap, reps, nk, taban=None):
        js = KAT[ad]["A2_j"]
        if any(nk[j - jmin] < 4 for j in js):
            return None, None
        ix = np.array(js) - jmin
        p = float(kap[ix].sum())
        pr = reps[:, ix].sum(1)
        if taban is not None:
            p -= len(ix) * taban
            pr = pr - len(ix) * taban
        return p, pr

    ORAN = ONK["A2"]["oranlar"]
    a2 = {}
    say = 0
    print("\n  A2 — parlaklık oranları (DÜŞÜK; P = Σκ, 5 dilim; bant ×/÷2):")
    print(f"   {'oran':>12} {'öngörü':>7} {'bant':>15} | {'ölçülen':>16} {'bantta':>6} | "
          f"{'taban-düz. (KAYIT)':>18} | {'tepe oranı':>16} | {'SON (KAYIT)':>16}")
    for ad_or, ong in ORAN.items():
        pay, payda = ad_or.split("/")
        sonuc = {"ongoru": ong, "bant": [ong / 2, ong * 2]}
        # birincil
        P1, P1r = parlaklik(pay, PF["d_kap"], PF["d_kap_reps"], PF["d_nkap"])
        P2, P2r = parlaklik(payda, PF["d_kap"], PF["d_kap_reps"], PF["d_nkap"])
        if P1 is None or P2 is None:
            sonuc["olculen"] = "erişilemedi"
            a2[ad_or] = sonuc
            print(f"   {ad_or:>12} ERİŞİLEMEDİ")
            continue
        R = P1 / P2
        Rr = P1r / P2r
        sonuc.update({"P_pay": P1, "se_P_pay": jk(P1r), "P_payda": P2,
                      "se_P_payda": jk(P2r), "oran": R, "se_oran": jk(Rr),
                      "bantta": bool(ong / 2 - EPS <= R <= ong * 2 + EPS)})
        say += int(sonuc["bantta"])
        # KAYIT taban-düzeltmeli
        t1 = B["dusuk"][pay]["taban_medyan"]
        t2 = B["dusuk"][payda]["taban_medyan"]
        Q1, Q1r = parlaklik(pay, PF["d_kap"], PF["d_kap_reps"], PF["d_nkap"], t1)
        Q2, Q2r = parlaklik(payda, PF["d_kap"], PF["d_kap_reps"], PF["d_nkap"], t2)
        sonuc["KAYIT_taban_duz"] = {"oran": Q1 / Q2, "se": jk(Q1r / Q2r),
                                    "bantta": bool(ong / 2 <= Q1 / Q2 <= ong * 2)}
        # KAYIT tepe oranı
        j1 = B["dusuk"][pay]["tepe_j"] - jmin
        j2 = B["dusuk"][payda]["tepe_j"] - jmin
        tr = PF["d_kap"][j1] / PF["d_kap"][j2]
        trr = PF["d_kap_reps"][:, j1] / PF["d_kap_reps"][:, j2]
        sonuc["KAYIT_tepe_orani"] = {"oran": float(tr), "se": jk(trr),
                                     "bantta": bool(ong / 2 <= tr <= ong * 2)}
        # KAYIT son
        S1, S1r = parlaklik(pay, PF["s_kap"], PF["s_kap_reps"], PF["s_nkap"])
        S2, S2r = parlaklik(payda, PF["s_kap"], PF["s_kap_reps"], PF["s_nkap"])
        if S1 is not None and S2 is not None:
            sonuc["KAYIT_son"] = {"P_pay": S1 * 0.025, "P_payda": S2 * 0.025,
                                  "oran": S1 / S2, "se": jk(S1r / S2r),
                                  "bantta": bool(ong / 2 <= S1 / S2 <= ong * 2)}
            ss = f"{S1/S2:.3f}±{jk(S1r/S2r):.3f}"
        else:
            sonuc["KAYIT_son"] = "erişilemedi"
            ss = "erişilemedi"
        a2[ad_or] = sonuc
        print(f"   {ad_or:>12} {ong:7.3f} [{ong/2:.3f},{ong*2:.3f}] | {R:7.3f}±{jk(Rr):.3f}"
              f"  {'EVET' if sonuc['bantta'] else 'HAYIR':>6} | "
              f"{Q1/Q2:7.3f}±{jk(Q1r/Q2r):.3f} | {tr:7.3f}±{jk(trr):.3f} | {ss}")
    out["A2"] = a2
    print(f"   bantta: {say}/{len(ORAN)}")
    for ad in ("+log2", "+log3", "+log5", "+log6", "+log7", "+log10", "-log2", "-log3"):
        P1, P1r = parlaklik(ad, PF["d_kap"], PF["d_kap_reps"], PF["d_nkap"])
        out.setdefault("A2_P", {})[ad] = {"P": P1, "se": jk(P1r)}
        print(f"     P({ad}) = {P1:.4e} ± {jk(P1r):.1e}")

    # ======================= H-192a =======================
    D = B["dusuk"]
    son_list = [a for a in KAT if KAT[a]["liste"] == "sonmeli"]
    ana6 = ONK["B"]["ana6"]
    olum = [a for a in son_list if D[a].get("yanar_4se")]
    hepsi_soner = all(D[a]["hukum"] == "söner" for a in son_list)
    ana_yanar = all(D[a]["hukum"] == "yanar" for a in ana6)
    if olum:
        hA = "ÖLDÜ"
    elif hepsi_soner and ana_yanar:
        hA = "MÜHÜR"
    else:
        eks = [a for a in son_list if D[a]["hukum"] != "söner"]
        eky = [a for a in ana6 if D[a]["hukum"] != "yanar"]
        parc = []
        if eks:
            parc.append("sönmeyen sönmeli: " + ", ".join(
                a + " (" + D[a]["hukum"] + (", çukur" if D[a]["cukur"] else "") + ")"
                for a in eks))
        if eky:
            parc.append("yanmayan ana-6: " + ", ".join(eky))
        hA = "KAYIT (kısmi: " + "; ".join(parc) + ")"
    Sn = B["son"]
    olum_s = [a for a in son_list if Sn[a].get("yanar_4se")]
    out["H192a"] = {
        "hukum": hA, "olum_adaylari": olum,
        "sonmeli": {a: {"hukum": D[a]["hukum"], "cukur": D[a]["cukur"],
                        "kappa_bolu_se": D[a]["kappa_bolu_se"],
                        "soner_nedeni": ("|κ|<2se" if abs(D[a]["kappa_tepe"]) < 2 * D[a]["se"]
                                         else ("taban içinde" if D[a]["hukum"] == "söner" else "—"))}
                    for a in son_list},
        "ana6": {a: {"hukum": D[a]["hukum"], "kappa_bolu_se": D[a]["kappa_bolu_se"]}
                 for a in ana6},
        "yanmali_diger": {a: D[a]["hukum"] for a in KAT
                          if KAT[a]["liste"] == "yanmali" and a not in ana6},
        "son_KAYIT": {"sonmeli": {a: Sn[a]["hukum"] for a in son_list},
                      "ana6": {a: Sn[a]["hukum"] for a in ana6},
                      "olum_adaylari_4se": olum_s}}

    # ======================= H-192b =======================
    U = A["uydular"]

    def oran_ok(ad, lo, hi, sin):
        return all(lo <= U[ad]["siniflar"][s]["oran_top"] <= hi for s in sin)

    d3p = U["+log3"]["farklar"]["1->2"]["fark"]
    d3m = U["-log3"]["farklar"]["1->2"]["fark"]
    k3p = (105 <= abs(d3p) <= 135) and oran_ok("+log3", 0.7, 1.4, ["r1", "r2"])
    k3m = (abs(d3m) < 15) and oran_ok("-log3", 0.35, 0.65, ["r1", "r2"])
    olB = (abs(d3p) < 45) or (abs(d3m) > 45)
    hB = "ÖLDÜ" if olB else ("MÜHÜR" if (k3p and k3m) else "KAYIT")
    # yan öngörüler (KAYIT)
    f5 = U["+log5"]["farklar"]
    k5 = all(57 <= abs(f5[f"{r}->{r+1}"]["fark"]) <= 87 for r in (1, 2, 3)) and \
        oran_ok("+log5", 0.6, 1.5, ["r1", "r2", "r3", "r4"])
    d6 = U["+log6"]["farklar"]["1->5"]["fark"]
    k6 = 105 <= abs(d6) <= 135
    r61 = U["+log6"]["siniflar"]["r1"]["aci_goreli_top"]
    r65 = U["+log6"]["siniflar"]["r5"]["aci_goreli_top"]
    k6yan = (45 <= abs(r61) <= 75) and (45 <= abs(r65) <= 75) and (np.sign(r61) != np.sign(r65))
    t2 = U["+log2"]["tek_sinif_mod2"]
    k2 = abs(t2["aci_goreli_top"]) < 15
    d24 = U["+log2"]["farklar"]["r1mod4->r3mod4"]["fark"]
    o24 = [U["+log2"]["siniflar"][s]["oran_top"] for s in ("r1mod4", "r3mod4")]
    s3 = np.sign(d3p)
    s5 = [np.sign(f5[f"{r}->{r+1}"]["fark"]) for r in (1, 2, 3)]
    s6 = -np.sign(d6)
    yon = bool(all(x == s3 for x in s5) and s6 == s3)
    out["H192b"] = {
        "hukum": hB,
        "+log3": {"fark_1_2": d3p, "se": U["+log3"]["farklar"]["1->2"]["se"],
                  "oranlar": {s: U["+log3"]["siniflar"][s]["oran_top"] for s in ("r1", "r2")},
                  "kosul": bool(k3p)},
        "-log3": {"fark_1_2": d3m, "se": U["-log3"]["farklar"]["1->2"]["se"],
                  "oranlar": {s: U["-log3"]["siniflar"][s]["oran_top"] for s in ("r1", "r2")},
                  "kosul": bool(k3m)},
        "olum_kosulu": bool(olB),
        "yan_KAYIT": {
            "+log5": {"farklar": {k: v["fark"] for k, v in f5.items()},
                      "oranlar": {s: U["+log5"]["siniflar"][s]["oran_top"]
                                  for s in ("r1", "r2", "r3", "r4")}, "kosul": bool(k5)},
            "+log6": {"fark_1_5": d6, "goreli": [r61, r65], "kosul": bool(k6),
                      "yan_pm60": bool(k6yan)},
            "+log2": {"tek_sinif_goreli": t2["aci_goreli_top"], "tek_sinif_oran": t2["oran_top"],
                      "kosul": bool(k2), "mod4_fark": d24, "mod4_oranlar": o24},
            "yon_tutarliligi": {"s_log3": float(s3), "s_log5": [float(x) for x in s5],
                                "s_log6_eksi": float(s6), "tutarli": yon}}}

    # ======================= H-192c =======================
    out["H192c"] = {"hukum": f"KAYIT ({say}/{len(ORAN)} oran bantta)",
                    "bantta": {k: v.get("bantta") for k, v in a2.items()}}

    out["hukum"] = {"H-192a": hA, "H-192b": hB, "H-192c": out["H192c"]["hukum"]}
    print("\nHÜKÜMLER (ön-kayıt eşikleri; kurtarma yok):")
    print(f"  H-192a: ana-6 = { {a: D[a]['hukum'] for a in ana6} }")
    print(f"          sönmeli = { {a: D[a]['hukum'] + ('/çukur' if D[a]['cukur'] else '') for a in son_list} }")
    print(f"          ölüm adayları (sönmeli yanar_4se) = {olum}  → {hA}")
    print(f"  H-192b: +log3 Δ12 = {d3p:+.1f}° (koşul {k3p}); −log3 Δ12 = {d3m:+.1f}° "
          f"(koşul {k3m}); ölüm {olB} → {hB}")
    print(f"          yan: +log5 {k5} ({[round(f5[k]['fark'], 1) for k in f5]}); +log6 {k6} "
          f"(Δ15 {d6:+.1f}°, göreli {r61:+.1f}/{r65:+.1f}); +log2 {k2}; "
          f"mod4 Δ {d24:+.1f}° oranlar {np.round(o24, 3).tolist()}; yön tutarlı {yon}")
    print(f"  H-192c: {out['H192c']['hukum']}")
    json.dump(out, open(S192 / "HUKUM_192.json", "w"), indent=1, ensure_ascii=False)
    print(f"\n-> {S192/'HUKUM_192.json'}  BİTTİ")
