# -*- coding: utf-8 -*-
"""
189c — K2: ZİNCİR SÜRÜCÜSÜ (184b2 → 185b → 186b → 187c; içerik değişikliği YOK)
===============================================================================
ONKAYIT_189 dondurdu. Modüller importlib ile yüklenir; DOSYALAR DÜZENLENMEZ.
Yalnız çalışma-anı yamaları:
  • b185.kinematik (185b'nin ve 186b'nin içindeki iki b185 örneğinde):
    bilinen adlar (gercek/Hkeskin/HA4) ORİJİNAL fonksiyona gider; yeni adlar
    (Hderin110/Hderin120) 155/z_<ad>.npy'den 185b'nin Hkeskin koluyla AYNI
    satırlarla (g = diff z, m = orta nokta, L = mean log(m/2π)).
  • çıktı dizinleri: b184b2.S184 → 189, b185.S185 → 189 (ön-kayıt SHA
    kontrolünden SONRA). 186b izdüşümü yol argümanıyla 189'a yazar.
  • __main__ blokları (185b defteri, 187c ζ defteri) burada aynı işlem
    sırasıyla yeniden kurulur; 185b defterine HAVUZ satırı EKLENİR (ilk 9 satır
    185b ile bit-bit aynı olmalı — makine mührü).
Mevcut 184-188 önbellekleri yalnız OKUNUR.

Kullanım: 189c_zincir.py <ad> [<ad> …]    (ör. Hkeskin  Hderin110  Hderin120)
Çıktı: 189/K1_<ad>.npz, OZ_<ad>.npz, G1_proj_<ad>.npz, zincir_<ad>.{json,npz},
       (Hkeskin için) muhur_Hkeskin.json
"""
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S155, S184, S185, S186, S187, S189 = (SCR / d for d in
                                      ["155", "184", "185", "186", "187", "189"])
TWO_PI = 2 * np.pi
NJACK = 8
BILINEN = ("gercek", "Hkeskin", "HA4")


def yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec)
    sys.modules[ad] = m
    spec.loader.exec_module(m)
    return m


b184b2 = yukle("b184b2", QM / "184_configs" / "184b2_zarf_ikiz.py")
b185 = yukle("b185", QM / "185_configs" / "185b_oz_muhasebe.py")
b186 = yukle("b186", QM / "186_configs" / "186b_g1_oztutarlilik.py")
b187 = yukle("b187", QM / "187_configs" / "187c_zeta_defteri.py")


def kinematik_yamasi(modul):
    orijinal = modul.kinematik

    def kinematik(gaz):
        if gaz in BILINEN:
            return orijinal(gaz)
        # 185b Hkeskin kolu AYNEN (yalnız dosya adı gaz adından)
        z = np.sort(np.load(S155 / f"z_{gaz}.npy").astype(float))
        g = np.diff(z)
        mid = 0.5 * (z[:-1] + z[1:])
        L = float(np.log(mid / TWO_PI).mean())
        return g, mid, L
    modul.kinematik = kinematik


kinematik_yamasi(b185)
kinematik_yamasi(b186.b185)

ONK185 = b185.onkayit()
ONK186 = b186.onkayit()
ONK187 = b187.onkayit()
ONK189 = json.load(open(S189 / "ONKAYIT_189.json"))
b184b2.S184 = S189
b185.S185 = S189

G = np.load(S184 / "K1_gercek.npz")
q, w, mq = G["q"], G["w"], G["mq"]
tau = np.asarray(G["tau"], float)
aq = np.asarray(G["aq"], float)
ae = np.asarray(G["aq_eff"], float)
OZg = np.load(S185 / "OZ_gercek.npz")
Pg = np.load(S186 / "G1_proj_gercek.npz")

# 185b maskeleri AYNEN (+ HAVUZ satırı eklenir)
kenar = ONK185["kenar"]
MASK185 = [(tau >= kenar[b]) & (tau < kenar[b + 1]) for b in range(8)]
MASK185.append(tau > ONK185["kuyruk_tau"])
MASK185.append((tau >= kenar[0]) & (tau < kenar[-1]))
ET185 = [f"{kenar[b]:.2f}-{kenar[b+1]:.2f}" for b in range(8)] + \
    ["KUYRUK>0.70", "HAVUZ"]
# 187c maskeleri AYNEN
kenar_b = ONK186["kenar"]
MASK187 = [(tau >= kenar_b[b]) & (tau < kenar_b[b + 1]) for b in range(8)]
MASK187.append((tau >= kenar_b[0]) & (tau < kenar_b[-1]))
ET187 = [f"{kenar_b[b]:.2f}-{kenar_b[b+1]:.2f}" for b in range(8)] + ["HAVUZ"]


# ------------------------------------------------ 185b defteri (AYNEN + HAVUZ)
def defter185(H, OZh, disari):
    ag = np.abs(b185.c_olculu(G, disari))
    ah = np.abs(b185.c_olculu(H, disari))
    og = np.abs(b185.c_oz(OZg, aq, disari))
    oh = np.abs(b185.c_oz(OZh, aq, disari))
    satir = {}
    for k, m in enumerate(MASK185):
        r = ag[m].sum() / ah[m].sum()
        F1 = ag[m].sum() / og[m].sum()
        F2 = og[m].sum() / oh[m].sum()
        F3i = oh[m].sum() / ah[m].sum()
        F1c = np.sum(ae[m] * (ag[m] / og[m])) / ae[m].sum()
        F2c = np.sum(ae[m] * (og[m] / oh[m])) / ae[m].sum()
        F3ic = np.sum(ae[m] * (oh[m] / ah[m])) / ae[m].sum()
        wg = ag[m].sum() / ae[m].sum()
        wh = ah[m].sum() / ae[m].sum()
        wog = og[m].sum() / ae[m].sum()
        woh = oh[m].sum() / ae[m].sum()
        satir[k] = (r, F1, F2, F3i, F1c, F2c, F3ic, wg, wh, wog, woh)
    return np.array([satir[k] for k in range(len(MASK185))])   # (10, 11)


# ------------------------------------------------ 187c ζ defteri (AYNEN)
def zeta187(D, OZ, P):
    def zeta_defteri(disari):
        c_olc = b185.c_olculu(D, disari)
        c_oz = b185.c_oz(OZ, aq, disari)
        ck = b187.c_kesik(P, disari)
        delta = c_olc - ck
        mix = ck - c_oz
        out = np.zeros((9, 4))
        zc = np.zeros(9, complex)
        for k, m in enumerate(MASK187):
            z = np.sum(delta[m] * np.conj(mix[m])) / \
                np.sum(np.abs(mix[m]) ** 2)
            sae = ae[m].sum()
            soz = np.abs(c_oz[m]).sum()
            m_olc = (np.abs(c_olc[m]).sum() - soz) / sae
            m_kes = (np.abs(ck[m]).sum() - soz) / sae
            out[k] = [np.abs(z), np.degrees(np.angle(z)), m_olc, m_kes]
            zc[k] = z
        return out, zc

    tam, zc = zeta_defteri(-1)
    rr = [zeta_defteri(i) for i in range(NJACK)]
    reps = np.array([x[0] for x in rr])
    zreps = np.array([x[1] for x in rr])
    se_mod = b185.jk_se(reps[:, :, 0])
    dfark = b187.sar(reps[:, :, 1] - tam[None, :, 1])
    se_aci = np.sqrt((NJACK - 1) / NJACK *
                     np.sum((dfark - dfark.mean(0)) ** 2, 0))
    genlik = 1.0 - tam[:, 2] / tam[:, 3]
    g_reps = 1.0 - reps[:, :, 2] / reps[:, :, 3]
    se_genlik = b185.jk_se(g_reps)
    icinde = np.abs(b187.sar(tam[:, 1] - 180.0)) <= 15.0
    return dict(tam=tam, reps=reps, zc=zc, zreps=zreps, se_mod=se_mod,
                se_aci=se_aci, genlik=genlik, g_reps=g_reps,
                se_genlik=se_genlik, icinde=icinde)


def sigma_eps(OZ):
    N = int(OZ["N"])
    nb = OZ["nb"]
    gbar = float(OZ["gbar"])
    sdg, sdg2 = OZ["sdg"], OZ["sdg2"]
    tam = np.sqrt(sdg2.sum() / N) / gbar          # 185b oz_hesapla AYNEN
    reps = []
    for b in range(NJACK):
        n_ = N - nb[b]
        s1 = sdg.sum() - sdg[b]
        s2 = sdg2.sum() - sdg2[b]
        mu = s1 / n_
        reps.append(np.sqrt(s2 / n_ - mu ** 2) / (gbar + mu))
    reps = np.array(reps)
    return float(tam), float(b185.jk_se(reps[:, None])[0]), reps


def zincir(ad):
    t0 = time.time()
    print(f"\n{'='*78}\n189c ZİNCİR [{ad}]\n{'='*78}", flush=True)
    zf = S155 / f"z_{ad}.npy"
    k1 = S189 / f"K1_{ad}.npz"
    if not k1.exists():
        print("  184b2 …", flush=True)
        b184b2.kos(ad, str(zf), q, w, mq)
    oz = S189 / f"OZ_{ad}.npz"
    if not oz.exists():
        print("  185b oz_hesapla …", flush=True)
        b185.oz_hesapla(ad, np.asarray(w, float))
    pj = S189 / f"G1_proj_{ad}.npz"
    if not pj.exists():
        print("  186b merdiven_izdusum (ρ≡1) …", flush=True)
        rhos = b186.rho_vektorleri(tau)
        b186.merdiven_izdusum(ad, np.asarray(w, float), aq,
                              {"bir": rhos["bir"]}, pj)
    H = np.load(k1)
    OZh = np.load(oz)
    P = np.load(pj)
    if not (np.array_equal(OZh["nb"], G["nb"]) and np.array_equal(H["nb"], G["nb"])
            and np.array_equal(P["nb"], G["nb"])):
        raise SystemExit("BLOK KENARLARI 184 İLE UYUMSUZ")
    if not (np.allclose(H["aq_eff"], ae) and np.allclose(H["tau"], tau)):
        raise SystemExit("EVREN UYUMSUZ (186b denetimi)")

    tam = defter185(H, OZh, -1)
    reps = np.array([defter185(H, OZh, i) for i in range(NJACK)])
    se = b185.jk_se(reps)
    r_b, wg_b, wh_b, wog_b, woh_b = (tam[:, j] for j in (0, 7, 8, 9, 10))
    sig_r = r_b * np.sqrt((se[:, 7] / wg_b) ** 2 + (se[:, 8] / wh_b) ** 2)
    sig_r_ortak = se[:, 0]
    tau_bar = np.array([np.sum(ae[m] * tau[m]) / ae[m].sum() for m in MASK185])
    m_h = wh_b - woh_b
    se_m_h = b185.jk_se(reps[:, :, 8] - reps[:, :, 10])

    Z = zeta187(H, OZh, P)
    sg, sg_se, sg_reps = sigma_eps(OZh)
    negre = -np.real(Z["zc"])
    se_negre = b185.jk_se(-np.real(Z["zreps"]))

    print(f"  L={float(H['L']):.9f}  σ_ε={sg:.5f}±{sg_se:.5f}", flush=True)
    print(f"{'bant':>12} {'τ̄':>6} {'w_g':>7} {'w_HD':>14} {'w^öz_HD':>14} "
          f"{'m_HD':>14} {'r=g/HD':>15}")
    for k in range(10):
        print(f"{ET185[k]:>12} {tau_bar[k]:6.3f} {wg_b[k]:7.4f} "
              f"{wh_b[k]:7.4f}±{se[k,8]:.4f} {woh_b[k]:7.4f}±{se[k,10]:.4f} "
              f"{m_h[k]:7.4f}±{se_m_h[k]:.4f} {r_b[k]:7.4f}±{sig_r[k]:.4f}")
    print(f"\n  ζ DEFTERİ (187c AYNEN) [{ad}]:")
    for k in range(9):
        print(f"{ET187[k]:>12} |ζ|={Z['tam'][k,0]:.4f}±{Z['se_mod'][k]:.4f} "
              f"∠{Z['tam'][k,1]:+8.2f}±{Z['se_aci'][k]:5.2f} "
              f"−Reζ={negre[k]:.4f}±{se_negre[k]:.4f} "
              f"genlik={Z['genlik'][k]:.4f}±{Z['se_genlik'][k]:.4f} "
              f"m_ölç={Z['tam'][k,2]:.4f} m^kesik={Z['tam'][k,3]:.4f}")

    np.savez_compressed(
        S189 / f"zincir_{ad}.npz", et185=np.array(ET185), et187=np.array(ET187),
        d185_tam=tam, d185_reps=reps, d185_se=se, sig_r=sig_r,
        sig_r_ortak=sig_r_ortak, tau_bar=tau_bar,
        z_tam=Z["tam"], z_reps=Z["reps"], zc=Z["zc"], zreps=Z["zreps"],
        z_se_mod=Z["se_mod"], z_se_aci=Z["se_aci"], genlik=Z["genlik"],
        g_reps=Z["g_reps"], se_genlik=Z["se_genlik"],
        sigma_eps=sg, sigma_eps_se=sg_se, sigma_eps_reps=sg_reps)
    js = {"ad": ad, "L": float(H["L"]), "N": int(H["N"]),
          "etiket185": ET185, "tau_bar": tau_bar.tolist(),
          "r": r_b.tolist(), "sig_r": sig_r.tolist(),
          "sig_r_ortak": sig_r_ortak.tolist(),
          "w_g": wg_b.tolist(), "se_w_g": se[:, 7].tolist(),
          "w_HD": wh_b.tolist(), "se_w_HD": se[:, 8].tolist(),
          "w_oz_g": wog_b.tolist(), "se_w_oz_g": se[:, 9].tolist(),
          "w_oz_HD": woh_b.tolist(), "se_w_oz_HD": se[:, 10].tolist(),
          "m_HD": m_h.tolist(), "se_m_HD": se_m_h.tolist(),
          "sigma_eps": sg, "se_sigma_eps": sg_se,
          "etiket187": ET187,
          "zeta_mod": Z["tam"][:, 0].tolist(), "se_mod": Z["se_mod"].tolist(),
          "zeta_aci": Z["tam"][:, 1].tolist(), "se_aci": Z["se_aci"].tolist(),
          "pencere_180pm15": Z["icinde"].tolist(),
          "neg_re_zeta": negre.tolist(), "se_neg_re_zeta": se_negre.tolist(),
          "genlik_okuma": Z["genlik"].tolist(),
          "se_genlik": Z["se_genlik"].tolist(),
          "m_olc": Z["tam"][:, 2].tolist(), "m_kesik": Z["tam"][:, 3].tolist(),
          "sure_s": time.time() - t0}
    json.dump(js, open(S189 / f"zincir_{ad}.json", "w"), indent=1,
              ensure_ascii=False)
    print(f"  -> zincir_{ad}.json ({time.time()-t0:.0f}s)", flush=True)
    return tam, se, sig_r, tau_bar, Z


def gercek_zeta():
    """Gerçeğin 187c defteri (184/185/186 dosyalarından; K3 model A için)."""
    Z = zeta187(G, OZg, Pg)
    sg, sg_se, sg_reps = sigma_eps(OZg)
    np.savez_compressed(
        S189 / "zincir_gercek.npz", z_tam=Z["tam"], z_reps=Z["reps"],
        zc=Z["zc"], zreps=Z["zreps"], z_se_mod=Z["se_mod"],
        z_se_aci=Z["se_aci"], genlik=Z["genlik"], g_reps=Z["g_reps"],
        se_genlik=Z["se_genlik"], sigma_eps=sg, sigma_eps_se=sg_se,
        sigma_eps_reps=sg_reps)
    return Z, sg, sg_se


def esit(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    return bool(a.shape == b.shape and np.array_equal(a, b))


def muhur_hkeskin(tam, se, sig_r, tau_bar, Z, Zg):
    """AYNI sarmalayıcı zinciri Hkeskin'de 184/185/186/187'yi bit-bit vermeli."""
    s = {}
    A, B = np.load(S189 / "K1_Hkeskin.npz"), np.load(S184 / "K1_Hkeskin.npz")
    s["K1_Hkeskin_(184b2)"] = all(esit(A[k], B[k]) for k in B.files)
    A, B = np.load(S189 / "OZ_Hkeskin.npz"), np.load(S185 / "OZ_Hkeskin.npz")
    s["OZ_Hkeskin_(185b)"] = all(esit(A[k], B[k]) for k in B.files)
    A, B = np.load(S189 / "G1_proj_Hkeskin.npz"), np.load(S186 / "G1_proj_Hkeskin.npz")
    s["G1_proj_bir_(186b)"] = all(esit(A[k], B[k]) for k in
                                  ("re_bir", "im_bir", "d_ss", "dyy_bir",
                                   "dxy_bir", "drr_bir", "nb", "w"))
    F = np.load(S185 / "K1_faktorler.npz", allow_pickle=True)
    ad = ["r", "F1", "F2", "F3i", "F1c", "F2c", "F3ic", "wg", "wh", "wog", "woh"]
    s["K1_faktorler_9satir_(185b)"] = (
        all(esit(tam[:9, j], F[k]) for j, k in enumerate(ad))
        and esit(se[:9], F["se_tablo"]) and esit(sig_r[:9], F["sig_r"])
        and esit(tau_bar[:9], F["tau_bar"]))
    K2 = json.load(open(S187 / "K2_zeta.json"))
    for gaz, ZZ in (("Hkeskin", Z), ("gercek", Zg)):
        ref = K2[gaz]
        s[f"K2_zeta_{gaz}_(187c)"] = (
            ZZ["tam"][:, 0].tolist() == ref["zeta_mod"]
            and ZZ["se_mod"].tolist() == ref["se_mod"]
            and ZZ["tam"][:, 1].tolist() == ref["zeta_aci"]
            and ZZ["se_aci"].tolist() == ref["se_aci"]
            and ZZ["genlik"].tolist() == ref["genlik_okuma"]
            and ZZ["se_genlik"].tolist() == ref["se_genlik"]
            and ZZ["tam"][:, 2].tolist() == ref["m_olc"]
            and ZZ["tam"][:, 3].tolist() == ref["m_kesik"])
    s["HEPSI"] = all(s.values())
    json.dump(s, open(S189 / "muhur_Hkeskin.json", "w"), indent=1,
              ensure_ascii=False)
    print("\nMAKİNE MÜHRÜ (Hkeskin, aynı sarmalayıcı zinciri; bit-bit):")
    for k, v in s.items():
        print(f"  {k:32s} {'TUTTU' if v else 'TUTMADI'}")
    return s


if __name__ == "__main__":
    print(f"189c ZİNCİR  [ONKAYIT_189 {ONK189['zaman']} sha {ONK189['sha256'][:12]}; "
          f"185 {ONK185['sha256'][:12]} 186 {ONK186['sha256'][:12]} "
          f"187 {ONK187['sha256'][:12]}]", flush=True)
    Zg, sgg, sgg_se = gercek_zeta()
    print(f"  gerçek: σ_ε={sgg:.5f}±{sgg_se:.5f}  |ζ_g|(HAVUZ)="
          f"{Zg['tam'][8,0]:.4f}  genlik={Zg['genlik'][8]:.4f}", flush=True)
    for ad in sys.argv[1:]:
        tam, se, sig_r, tau_bar, Z = zincir(ad)
        if ad == "Hkeskin":
            muhur_hkeskin(tam, se, sig_r, tau_bar, Z, Zg)
    print("\nBİTTİ", flush=True)
