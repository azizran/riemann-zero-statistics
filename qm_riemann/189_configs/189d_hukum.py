# -*- coding: utf-8 -*-
"""
189d — K3: HÜKÜM (H-189a..d) + KİMLİK KAYDI
============================================
ONKAYIT_189'da donan tanımlar ve eşiklerle (kurtarma yok):
  H-189a  zarf = merdiven derinliği (yön + f_1.10 ≥ 0.30; ölüm f < 0.15 / azalma)
  H-189b  sabit-kinematik aritmetiği (bant başına ±max(2σ, 0.01)) + ayrışım
  H-189c  σ_ε monotonluğu (geri-besleme kinematiği)
  H-189d  harita geçerliliği |ζ_HD| vs z_D (HAVUZ ≤ %10)
  K3      r_D(τ) → 1 − c·τ^α (184c makinesi AYNEN), D→∞ (model A birincil,
          B ikincil), derinlik/kinematik payları.
Girdi: 189/zincir_{Hkeskin,Hderin110,Hderin120}.{json,npz}, zincir_gercek.npz,
       muhur_Hkeskin.json, izgara_*.json, insa_kapi_*.json, ONKAYIT_189.json,
       188/harita_K_Hkeskin.npz, 184/K1_gercek.npz, 189/K1_<ad>.npz.
Çıktı: 189/HUKUM_189.json (+ ekran defteri).
"""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S184, S188, S189 = SCR / "184", SCR / "188", SCR / "189"
NJACK = 8
DER = {"1.00": "Hkeskin", "1.10": "Hderin110", "1.20": "Hderin120"}
BANT = ["0.45-0.50", "0.50-0.55", "0.55-0.60", "0.60-0.65",
        "0.65-0.70", "0.70-0.75", "0.75-0.80", "0.80-0.86"]
IH = 9          # 185-defterinde HAVUZ satırı (0..7 bant, 8 KUYRUK, 9 HAVUZ)
IH7 = 8         # 187c defterinde HAVUZ satırı

spec = importlib.util.spec_from_file_location(
    "b184c", QM / "184_configs" / "184c_K2_yarisma.py")
b184c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b184c)


def jk(reps):
    reps = np.asarray(reps, float)
    return np.sqrt((NJACK - 1) / NJACK * np.sum((reps - reps.mean(0)) ** 2, 0))


def onkayit():
    o = json.load(open(S189 / "ONKAYIT_189.json"))
    sha = hashlib.sha256(
        (QM / "189_configs" / "189a_onkayit.py").read_bytes()).hexdigest()
    if sha != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {sha} != {o['sha256']}")
    return o


def yukle_zincir(ad):
    js = json.load(open(S189 / f"zincir_{ad}.json"))
    nz = np.load(S189 / f"zincir_{ad}.npz", allow_pickle=True)
    return js, nz


def fit184c(g, tw):
    """184c makinesi AYNEN: bant_defteri + H_Z2 fit; ek: loo-jk yeniden fit."""
    sat = b184c.bant_defteri(g, tw, "x")[:8]
    tau = np.array([s["tau"] for s in sat])
    r = np.array([s["r"] for s in sat])
    se = np.array([s["se_r"] for s in sat])
    p0, bnd = [0.3, 3.0], ([0.0, 0.1], [10.0, 20.0])
    popt, pcov, chi2, err = b184c.fit(b184c.H_Z2, tau, r, se, p0, bnd)
    perr = np.sqrt(np.diag(pcov)) if pcov is not None else [np.nan, np.nan]
    # loo jackknife: r_k bant başına (184c maskeleri), aynı se ile yeniden fit
    ah_g, jk_g = b184c.band_ahat(g)
    ah_t, jk_t = b184c.band_ahat(tw)
    KEN = b184c.KENAR
    reps = []
    for k in range(NJACK):
        rk = []
        for i in range(8):
            msk = (g["tau"] > KEN[i]) & (g["tau"] <= KEN[i + 1])
            rk.append(jk_g[k][msk].sum() / jk_t[k][msk].sum())
        pk, _, _, _ = b184c.fit(b184c.H_Z2, tau, np.array(rk), se, p0, bnd)
        reps.append(pk if pk is not None else [np.nan, np.nan])
    reps = np.array(reps)
    sinir = bool(popt is not None and (popt[0] <= 1e-9 or popt[0] >= 10 - 1e-9
                                       or popt[1] <= 0.1 + 1e-9
                                       or popt[1] >= 20 - 1e-9))
    return dict(tau=tau.tolist(), r=r.tolist(), se_r=se.tolist(),
                c=float(popt[0]), alfa=float(popt[1]),
                perr=[float(x) for x in perr], jk_se=jk(reps).tolist(),
                chi2_dof=float(chi2 / 6), sinirda=sinir, hata=err)


def fit_serbest(tau, r, se):
    """ÖN-KAYITSIZ EK (hüküm dışı): işaret-serbest 1 − c·τ^α (c < 0 izinli)."""
    from scipy.optimize import curve_fit
    f = lambda t, c, a: 1.0 - c * t ** a
    try:
        popt, pcov = curve_fit(f, tau, r, p0=[0.1, 1.3], sigma=se,
                               absolute_sigma=True,
                               bounds=([-10.0, 0.1], [10.0, 20.0]), maxfev=20000)
        chi2 = float(np.sum(((r - f(tau, *popt)) / se) ** 2))
        return dict(c=float(popt[0]), alfa=float(popt[1]),
                    perr=np.sqrt(np.diag(pcov)).tolist(), chi2_dof=chi2 / 6)
    except Exception as e:
        return dict(hata=str(e))


if __name__ == "__main__":
    ONK = onkayit()
    print("=" * 78)
    print(f"189d / HÜKÜM + KİMLİK KAYDI  [on-kayit {ONK['zaman']} "
          f"sha {ONK['sha256'][:12]}]")
    print("=" * 78, flush=True)
    out = {"sha_onkayit": ONK["sha256"], "damga_onkayit": ONK["zaman"]}

    # ---------------- kapılar ----------------
    muhur = json.load(open(S189 / "muhur_Hkeskin.json"))
    out["makine_muhru"] = muhur
    print(f"makine mührü (Hkeskin bit-bit): {'TUTTU' if muhur['HEPSI'] else 'TUTMADI'}")
    kapilar = {}
    for D in ("1.10", "1.20"):
        ad = DER[D]
        iz = json.load(open(S189 / f"izgara_{ad}.json"))
        kp = S189 / f"insa_kapi_{ad}.json"
        kapilar[ad] = {"izgara": iz,
                       "insa": json.load(open(kp)) if kp.exists() else None}
        print(f"  {ad}: ızgara {iz['n_fark']}/{iz['n_sinav']} → h={iz['h_secilen']}; "
              f"inşa kapıları "
              f"{'GEÇTİ' if kapilar[ad]['insa'] and kapilar[ad]['insa']['HEPSI'] else 'YOK/KALDI'}")
    out["kapilar"] = kapilar
    mevcut = {D: (D == "1.00" or (kapilar[DER[D]]["insa"] or {}).get("HEPSI", False))
              and (S189 / f"zincir_{DER[D]}.json").exists() for D in DER}
    if not muhur["HEPSI"]:
        print("MAKİNE MÜHRÜ TUTMADI — zincir ŞÜPHELİ, hüküm verilmez")
    Z = {D: yukle_zincir(DER[D]) for D in DER if mevcut[D]}
    ZG = np.load(S189 / "zincir_gercek.npz")
    js = {D: Z[D][0] for D in Z}
    nz = {D: Z[D][1] for D in Z}
    ong = ONK["ongoru"]["tablo_kalem_AYNEN"]

    # ---------------- bant tablosu ----------------
    print("\nBANT TABLOSU (185-defter; r = w_g/w_HD ± σ_r [184 konv.]):")
    for D in Z:
        j = js[D]
        print(f"  D={D} ({DER[D]}): L={j['L']:.9f}")
        for k in list(range(8)) + [IH]:
            et = j["etiket185"][k]
            o = ""
            if et in ong and D != "1.00":
                o = f"  öng={ong[et][4 if D == '1.10' else 5]:.4f}"
            print(f"   {et:>11} w_HD={j['w_HD'][k]:.4f}±{j['se_w_HD'][k]:.4f} "
                  f"w^öz={j['w_oz_HD'][k]:.4f}±{j['se_w_oz_HD'][k]:.4f} "
                  f"m={j['m_HD'][k]:.4f}±{j['se_m_HD'][k]:.4f} "
                  f"r={j['r'][k]:.4f}±{j['sig_r'][k]:.4f}{o}")

    # ---------------- H-189a ----------------
    ha = {"erisilemedi": not all(mevcut.values())}
    if all(mevcut.values()):
        r = {D: np.array(js[D]["r"]) for D in Z}
        s = {D: np.array(js[D]["sig_r"]) for D in Z}
        sD = np.sqrt(s["1.10"] ** 2 + s["1.20"] ** 2)
        c1 = r["1.10"] > r["1.00"]
        c2 = r["1.10"] <= r["1.20"] + 2 * sD
        satir = list(range(8)) + [IH]
        bant_ok = [bool(c1[k] and c2[k]) for k in range(8)]
        f110 = (r["1.10"] - r["1.00"]) / (1 - r["1.00"])
        f120 = (r["1.20"] - r["1.00"]) / (1 - r["1.00"])
        # f jk (ortak blok loo)
        rr = {D: nz[D]["d185_reps"][:, :, 0] for D in Z}
        f110_se = jk((rr["1.10"] - rr["1.00"]) / (1 - rr["1.00"]))
        f120_se = jk((rr["1.20"] - rr["1.00"]) / (1 - rr["1.00"]))
        havuz_ok = bool(c1[IH] and c2[IH])
        kosul = havuz_ok and sum(bant_ok) >= 7 and f110[IH] >= 0.30
        olum = bool(f110[IH] < 0.15 or r["1.10"][IH] < r["1.00"][IH]
                    or r["1.20"][IH] < r["1.10"][IH] - 2 * sD[IH])
        hk = "ÖLDÜ" if olum else ("MÜHÜR" if kosul else "KAYIT")
        ha.update(dict(
            bant_c1=c1[:8].tolist(), bant_c2=c2[:8].tolist(),
            bant_gecen=int(sum(bant_ok)), havuz_c1=bool(c1[IH]),
            havuz_c2=bool(c2[IH]), sigma_D=sD[satir].tolist(),
            f110=f110[satir].tolist(), f110_se=f110_se[satir].tolist(),
            f120=f120[satir].tolist(), f120_se=f120_se[satir].tolist(),
            f110_havuz=float(f110[IH]), f120_havuz=float(f120[IH]),
            r_havuz={D: float(r[D][IH]) for D in r},
            sig_r_havuz={D: float(s[D][IH]) for D in s},
            hukum=hk))
        print(f"\nH-189a: HAVUZ r = {r['1.00'][IH]:.4f}±{s['1.00'][IH]:.4f} → "
              f"{r['1.10'][IH]:.4f}±{s['1.10'][IH]:.4f} → "
              f"{r['1.20'][IH]:.4f}±{s['1.20'][IH]:.4f}; "
              f"f_1.10 = {f110[IH]:.4f}±{f110_se[IH]:.4f}, f_1.20 = "
              f"{f120[IH]:.4f}±{f120_se[IH]:.4f}; bant {sum(bant_ok)}/8 → {hk}")
        for k in range(8):
            print(f"   {BANT[k]}: r {r['1.00'][k]:.4f}→{r['1.10'][k]:.4f}→"
                  f"{r['1.20'][k]:.4f}  f110={f110[k]:.3f}±{f110_se[k]:.3f} "
                  f"f120={f120[k]:.3f}±{f120_se[k]:.3f}  "
                  f"c1={'E' if c1[k] else 'H'} c2={'E' if c2[k] else 'H'}")
    out["H-189a"] = ha

    # ---------------- H-189b ----------------
    hb = {}
    F = np.load(SCR / "185" / "K1_faktorler.npz", allow_pickle=True)
    for D in ("1.10", "1.20"):
        if D not in Z:
            hb[D] = {"erisilemedi": True}
            continue
        j = js[D]
        col = 4 if D == "1.10" else 5
        zc = 2 if D == "1.10" else 3
        rows = []
        tut = 0
        for k, et in enumerate(BANT):
            rp = ong[et][col]
            rm, sr = j["r"][k], j["sig_r"][k]
            tol = max(2 * sr, 0.01)
            ok = abs(rm - rp) <= tol
            tut += ok
            # ayrışım: öngörü w_HD^öng = w^öz_Hk + m^kesik_Hk (1 − z_D)
            woh = float(F["woh"][k])
            mhk = float(F["wh"][k] - F["woh"][k])
            z100, zD = ong[et][1], ong[et][zc]
            # (tablo z'leri 4 haneli; aritmetik tam-hassas kaynaktan)
            yh = ONK["ongoru"]["yeniden_hesap"][et]
            z100, zD = yh[1], yh[zc]
            m_ong = mhk / (1 - z100) * (1 - zD)
            w_ong = woh + m_ong
            dw = j["w_HD"][k] - w_ong
            dwoz = j["w_oz_HD"][k] - woh
            dm = j["m_HD"][k] - m_ong
            rkin = j["w_g"][k] / (j["w_oz_HD"][k] + m_ong)
            rows.append(dict(bant=et, r_olc=rm, sig_r=sr, r_ong=rp,
                             sapma=rm - rp, tol=tol, tuttu=bool(ok),
                             w_HD=j["w_HD"][k], w_HD_ong=w_ong, dw=dw,
                             dw_oz_kin=dwoz, dm_karisim=dm,
                             kin_pay=dwoz / dw if dw != 0 else float("nan"),
                             r_kin=rkin, m_HD=j["m_HD"][k], m_HD_ong=m_ong,
                             F_KIN_HD=j["w_oz_g"][k] / j["w_oz_HD"][k]))
        hb[D] = dict(satirlar=rows, tutan=int(tut), tuttu=bool(tut == 8))
    if all(D in Z for D in ("1.10", "1.20")):
        hb["hukum"] = "MÜHÜR" if (hb["1.10"]["tuttu"] and hb["1.20"]["tuttu"]) else "ÖLDÜ"
    else:
        hb["hukum"] = "ERİŞİLEMEDİ"
    print(f"\nH-189b: {hb['hukum']}")
    for D in ("1.10", "1.20"):
        if "satirlar" not in hb[D]:
            continue
        print(f"  D={D}: {hb[D]['tutan']}/8 bant tolerans içinde")
        print(f"   {'bant':>10} {'r_ölç':>14} {'r_öng':>7} {'sapma':>8} {'tol':>6} "
              f"{'Δw':>7} {'Δw^öz(kin)':>10} {'Δm(karış)':>9} {'kin pay':>7} "
              f"{'r^kin':>7} {'F_KİN':>6}")
        for x in hb[D]["satirlar"]:
            print(f"   {x['bant']:>10} {x['r_olc']:.4f}±{x['sig_r']:.4f} "
                  f"{x['r_ong']:.4f} {x['sapma']:+.4f} {x['tol']:.4f} "
                  f"{x['dw']:+.4f} {x['dw_oz_kin']:+10.4f} {x['dm_karisim']:+9.4f} "
                  f"{x['kin_pay']:7.2f} {x['r_kin']:.4f} {x['F_KIN_HD']:.4f}")
    out["H-189b"] = hb

    # ---------------- H-189c ----------------
    hc = {}
    sg = {D: nz[D]["sigma_eps"].item() for D in Z}
    sg_se = {D: nz[D]["sigma_eps_se"].item() for D in Z}
    sg_reps = {D: nz[D]["sigma_eps_reps"] for D in Z}
    sgg = float(ZG["sigma_eps"])
    sgg_se = float(ZG["sigma_eps_se"])
    hc["sigma_eps"] = sg
    hc["se"] = sg_se
    hc["gercek"] = [sgg, sgg_se]
    if all(D in Z for D in DER):
        d1 = sg["1.10"] - sg["1.00"]
        d2 = sg["1.20"] - sg["1.10"]
        d1se = float(jk(sg_reps["1.10"] - sg_reps["1.00"]))
        d2se = float(jk(sg_reps["1.20"] - sg_reps["1.10"]))
        kosul = d1 < 0 and d2 < 0
        olum = d1 > 0 or d2 > 0
        hc.update(dict(
            fark_110_100=d1, se_fark_110_100=d1se,
            fark_120_110=d2, se_fark_120_110=d2se,
            kapanan_110=(sg["1.00"] - sg["1.10"]) / (sg["1.00"] - sgg),
            kapanan_120=(sg["1.00"] - sg["1.20"]) / (sg["1.00"] - sgg),
            gercegi_asti={D: bool(sg[D] < sgg) for D in ("1.10", "1.20")},
            hukum="ÖLDÜ" if olum else ("MÜHÜR" if kosul else "KAYIT")))
        print(f"\nH-189c: σ_ε Hk={sg['1.00']:.5f}±{sg_se['1.00']:.5f} → "
              f"110={sg['1.10']:.5f}±{sg_se['1.10']:.5f} → "
              f"120={sg['1.20']:.5f}±{sg_se['1.20']:.5f}  (gerçek "
              f"{sgg:.5f}±{sgg_se:.5f}); Δ1={d1:+.5f}±{d1se:.5f}, "
              f"Δ2={d2:+.5f}±{d2se:.5f}; kapanan pay 110={hc['kapanan_110']:.3f}, "
              f"120={hc['kapanan_120']:.3f} → {hc['hukum']}")
    else:
        hc["hukum"] = "ERİŞİLEMEDİ"
    out["H-189c"] = hc

    # ---------------- H-189d ----------------
    hd = {}
    hk = np.load(S188 / "harita_K_Hkeskin.npz", allow_pickle=True)
    et88 = list(hk["etiket"])
    orta = hk["orta"]

    def zD_harita(etiket, D):
        i = et88.index(etiket)
        sel = orta <= D + 1e-12
        z = float(-np.real(hk["K"][i, sel].sum()))
        zr = -np.real(hk["K_reps"][:, i, sel].sum(1))
        return z, float(jk(zr))

    for D in ("1.10", "1.20"):
        if D not in Z:
            hd[D] = {"erisilemedi": True}
            continue
        j = js[D]
        Dv = float(D)
        rows = []
        for k, et in enumerate(BANT + ["HAVUZ"]):
            z, zse = zD_harita("B" + et if et != "HAVUZ" else "HAVUZ", Dv)
            mod, mse = j["zeta_mod"][k], j["se_mod"][k]
            nre, nrse = j["neg_re_zeta"][k], j["se_neg_re_zeta"][k]
            rows.append(dict(bant=et, z_D=z, z_D_se=zse, zeta_mod=mod,
                             zeta_mod_se=mse, zeta_aci=j["zeta_aci"][k],
                             zeta_aci_se=j["se_aci"][k],
                             aci_180pm15=j["pencere_180pm15"][k],
                             delta=abs(mod - z) / z, delta_negre=abs(nre - z) / z,
                             neg_re=nre, neg_re_se=nrse,
                             genlik=j["genlik_okuma"][k]))
        hv = rows[-1]
        hd[D] = dict(satirlar=rows, havuz_delta=hv["delta"],
                     tuttu=bool(hv["delta"] <= 0.10))
    if all(D in Z for D in ("1.10", "1.20")):
        hd["hukum"] = "MÜHÜR" if (hd["1.10"]["tuttu"] and hd["1.20"]["tuttu"]) else "ÖLDÜ"
    else:
        hd["hukum"] = "ERİŞİLEMEDİ"
    print(f"\nH-189d: {hd['hukum']}")
    for D in ("1.10", "1.20"):
        if "satirlar" not in hd[D]:
            continue
        for x in hd[D]["satirlar"]:
            print(f"   D={D} {x['bant']:>10} |ζ_HD|={x['zeta_mod']:.4f}±"
                  f"{x['zeta_mod_se']:.4f} ∠{x['zeta_aci']:+7.2f}±{x['zeta_aci_se']:.2f} "
                  f"z_D={x['z_D']:.4f}±{x['z_D_se']:.4f} δ={100*x['delta']:.1f}% "
                  f"(−Reζ: {100*x['delta_negre']:.1f}%)")
    out["H-189d"] = hd

    # ---------------- K3 — KİMLİK KAYDI ----------------
    k3 = {"fit": {}, "fit_serbest_EK": {}}
    g = b184c.yukle("gercek")
    for D in Z:
        d = np.load(S189 / f"K1_{DER[D]}.npz")
        tw = {k: d[k] for k in d.files}
        fr = fit184c(g, tw)
        k3["fit"][D] = fr
        k3["fit_serbest_EK"][D] = fit_serbest(np.array(fr["tau"]),
                                             np.array(fr["r"]),
                                             np.array(fr["se_r"]))
    print("\nK3 — r_D(τ) = 1 − c·τ^α (184c makinesi AYNEN; ± perr / ± jk):")
    for D, fr in k3["fit"].items():
        ek = k3["fit_serbest_EK"][D]
        print(f"   D={D}: c={fr['c']:.4f}±{fr['perr'][0]:.4f}/±{fr['jk_se'][0]:.4f} "
              f"α={fr['alfa']:.3f}±{fr['perr'][1]:.3f}/±{fr['jk_se'][1]:.3f} "
              f"χ²/dof={fr['chi2_dof']:.3f}{'  [SINIRDA]' if fr['sinirda'] else ''}"
              f"   | EK işaret-serbest: c={ek.get('c', float('nan')):+.4f} "
              f"α={ek.get('alfa', float('nan')):.3f} χ²/dof={ek.get('chi2_dof', float('nan')):.3f}")

    # D→∞ — model A (birincil) ve B (ikincil)
    if all(D in Z for D in DER):
        rows7 = list(range(9))            # 187c satırları: 8 bant + HAVUZ
        rows5 = list(range(8)) + [IH]     # aynı bantların 185-defter satırları
        N120 = nz["1.20"]
        wg = N120["d185_tam"][rows5, 7]
        wg_r = N120["d185_reps"][:, rows5, 7]
        wog = N120["d185_tam"][rows5, 9]
        wog_r = N120["d185_reps"][:, rows5, 9]
        Gg = ZG["genlik"][rows7]
        Gg_r = ZG["g_reps"][:, rows7]
        zg = ZG["z_tam"][rows7, 0]
        zg_r = ZG["z_reps"][:, rows7, 0]
        mkg = ZG["z_tam"][rows7, 3]
        mkg_r = ZG["z_reps"][:, rows7, 3]
        r100 = nz["1.00"]["d185_tam"][rows5, 0]
        r100_r = nz["1.00"]["d185_reps"][:, rows5, 0]

        def model_A(woh, mk, Gam, wgv):
            return wgv / (woh + mk * (1 - Gam))

        woz = {D: nz[D]["d185_tam"][rows5, 10] for D in DER}
        woz_r = {D: nz[D]["d185_reps"][:, rows5, 10] for D in DER}
        mk = {D: nz[D]["z_tam"][rows7, 3] for D in DER}
        mk_r = {D: nz[D]["z_reps"][:, rows7, 3] for D in DER}
        Gm = {D: nz[D]["genlik"][rows7] for D in DER}
        Gm_r = {D: nz[D]["g_reps"][:, rows7] for D in DER}

        rA = model_A(woz["1.20"], mk["1.20"], Gg, wg)
        rA_r = model_A(woz_r["1.20"], mk_r["1.20"], Gg_r, wg_r)
        rA_z = model_A(woz["1.20"], mk["1.20"], zg, wg)
        rA_z_r = model_A(woz_r["1.20"], mk_r["1.20"], zg_r, wg_r)
        # kimlik sınaması: model A Γ=Γ_HD'de ölçülen r_HD'yi vermeli
        kimlik = {D: float(np.max(np.abs(model_A(woz[D], mk[D], Gm[D], wg)
                                          - nz[D]["d185_tam"][rows5, 0])))
                  for D in DER}
        Pder = (rA - r100) / (1 - r100)
        Pder_r = (rA_r - r100_r) / (1 - r100_r)
        Pder_z = (rA_z - r100) / (1 - r100)
        # kinematik fark ayrışımı: w_H∞ − w_g = Δöz + Δkesik·(1 − Γ_g)
        d_oz = woz["1.20"] - wog
        d_kes = (mk["1.20"] - mkg) * (1 - Gg)
        # model B: w^öz ve m^kesik Γ_HD'ye karşı doğrusal (3 derinlik)
        Gs = np.array([Gm[D] for D in DER])            # (3, 9)
        rB = np.zeros(9)
        rB_r = np.zeros((NJACK, 9))
        for k in range(9):
            x = Gs[:, k]
            a1 = np.polyfit(x, [woz[D][k] for D in DER], 1)
            a2 = np.polyfit(x, [mk[D][k] for D in DER], 1)
            rB[k] = wg[k] / (np.polyval(a1, Gg[k]) + np.polyval(a2, Gg[k]) * (1 - Gg[k]))
            for b in range(NJACK):
                xb = np.array([Gm_r[D][b, k] for D in DER])
                a1b = np.polyfit(xb, [woz_r[D][b, k] for D in DER], 1)
                a2b = np.polyfit(xb, [mk_r[D][b, k] for D in DER], 1)
                rB_r[b, k] = wg_r[b, k] / (np.polyval(a1b, Gg_r[b, k])
                                           + np.polyval(a2b, Gg_r[b, k]) * (1 - Gg_r[b, k]))
        PderB = (rB - r100) / (1 - r100)
        PderB_r = (rB_r - r100_r) / (1 - r100_r)
        et9 = BANT + ["HAVUZ"]
        k3["D_sonsuz"] = dict(
            etiket=et9, r100=r100.tolist(),
            A_r_inf=rA.tolist(), A_r_inf_se=jk(rA_r).tolist(),
            A_P_der=Pder.tolist(), A_P_der_se=jk(Pder_r).tolist(),
            A_P_kin=(1 - Pder).tolist(),
            A_zeta_duyarlilik_r_inf=rA_z.tolist(),
            A_zeta_duyarlilik_r_inf_se=jk(rA_z_r).tolist(),
            A_zeta_duyarlilik_P_der=Pder_z.tolist(),
            A_kimlik_sinamasi_maks=kimlik,
            A_ayrisim_w_oz=d_oz.tolist(), A_ayrisim_m_kesik=d_kes.tolist(),
            A_w_g=wg.tolist(), A_w_oz_g=wog.tolist(), A_m_kesik_g=mkg.tolist(),
            A_Gamma_g=Gg.tolist(), A_zeta_g=zg.tolist(),
            A_w_oz_H120=woz["1.20"].tolist(), A_m_kesik_H120=mk["1.20"].tolist(),
            B_r_inf=rB.tolist(), B_r_inf_se=jk(rB_r).tolist(),
            B_P_der=PderB.tolist(), B_P_der_se=jk(PderB_r).tolist(),
            Gamma_HD={D: Gm[D].tolist() for D in DER},
            w_oz_HD={D: woz[D].tolist() for D in DER},
            m_kesik_HD={D: mk[D].tolist() for D in DER})
        print("\nK3 — D→∞ (A: kinematik 1.20'de donmuş + Γ_g; B: kinematik "
              "eğilimi Γ'ye doğrusal):")
        print(f"   kimlik sınaması (A, Γ=Γ_HD ⇒ r_HD): maks|Δ| = "
              f"{max(kimlik.values()):.1e}")
        for k in range(9):
            print(f"   {et9[k]:>10} r100={r100[k]:.4f}  A: r_∞={rA[k]:.4f}±"
                  f"{jk(rA_r)[k]:.4f} P_der={Pder[k]:.3f}±{jk(Pder_r)[k]:.3f} "
                  f"(|ζ|-duyarlılık r_∞={rA_z[k]:.4f})  B: r_∞={rB[k]:.4f}±"
                  f"{jk(rB_r)[k]:.4f} P_der={PderB[k]:.3f}  | Δw^öz="
                  f"{d_oz[k]:+.4f} Δkesik(1−Γ)={d_kes[k]:+.4f}")
    out["K3"] = k3
    json.dump(out, open(S189 / "HUKUM_189.json", "w"), indent=1,
              ensure_ascii=False, default=float)
    print(f"\n-> {S189/'HUKUM_189.json'}  BİTTİ", flush=True)
