# -*- coding: utf-8 -*-
"""
190k0 — K0: ZİNCİR SÜRÜCÜSÜ (155 η → 184b → 185b → 186b → 187c; içerik değişikliği YOK)
=====================================================================================
KALEM_TARAK_EVRENSELLIK_23EYL2026 K0. Modüller importlib ile yüklenir; DOSYALAR
DÜZENLENMEZ; mevcut 155-188 önbellekleri yalnız OKUNUR. Çalışma-anı yamaları:
  • 155: k155.K.SCR → hedef dizin (η önbelleği 155'in KENDİ üreticisiyle:
    k155.veri_yukle(<pencere>) + K.eta_onbellek(z, <pencere>, 0.4, 4000)).
  • 184b: ONK = onkayit() (SHA denetimi ÖNCE), sonra b184.S155/S184 → hedef;
    b184.kos(<etiket>, <pencere>) AYNEN (pk_m(e^{0.86 L}) evreni, L = pencerenin L'si).
  • 185b: kinematik yaması (185b'nin ve 186b'nin içindeki iki b185 örneğinde):
    bilinen adlar ORİJİNAL fonksiyona; yeni etiket → 185b 'gercek' kolunun AYNI
    satırları, yalnız önbellek dosyası pencerenin η'sı. b185.S185 → hedef;
    oz_hesapla(<etiket>, w) AYNEN.
  • 186b: merdiven_izdusum(<etiket>, w, aq, {"bir": ρ≡1}, yol) AYNEN (189c emsali:
    ρ≡1 kolu diğer ρ'lardan bağımsız → re_bir/im_bir bit-bit aynı).
  • 187c: ζ defteri __main__ bloğunun işlem sırası AYNEN (189c.zeta187 emsali).
Pencereler:
  dusuk → scratchpad/190/zincir_dusuk/  (etiket 'gercek_dusuk')
  son   → scratchpad/190/muhur_son/     (GEÇİCİ etiket 'gercek_sonM'; MAKİNE MÜHRÜ:
          155 η, 184 K1, 185 OZ, 186 G1_proj(bir), 187 K2_zeta['gercek'] bit-bit)
Kullanım: 190k0_zincir.py <dusuk|son>
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
S155, S184, S185, S186, S187, S190 = (SCR / d for d in
                                      ["155", "184", "185", "186", "187", "190"])
TWO_PI = 2 * np.pi
NJACK = 8
AYAR = {"dusuk": ("zincir_dusuk", "gercek_dusuk"),
        "son": ("muhur_son", "gercek_sonM")}


def yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec)
    sys.modules[ad] = m
    spec.loader.exec_module(m)
    return m


def esit(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return bool(a.shape == b.shape and np.array_equal(a, b))


def zincir(pencere):
    t00 = time.time()
    dizin, etiket = AYAR[pencere]
    OUT = S190 / dizin
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"{'='*78}\n190k0 ZİNCİR [{pencere}] etiket={etiket} -> {OUT}\n{'='*78}",
          flush=True)

    # ---------------- 155 η önbelleği (155'in kendi üreticisi) ----------------
    k155 = yukle("k155", QM / "155_configs" / "155_kos.py")
    k155.K.SCR = OUT
    z = k155.veri_yukle(pencere)
    print(f"  [155] {pencere}: n={len(z)} t∈[{z[0]:.3f},{z[-1]:.3f}]", flush=True)
    t0 = time.time()
    C = k155.K.eta_onbellek(z, pencere, 0.4, 4000)
    mid, ds, L = C["mid"], C["ds"], C["L"]
    g = np.diff(z)
    g2 = (ds + 1.0) * TWO_PI / np.log(mid / TWO_PI)
    tut = float(np.max(np.abs(g - g2)))
    print(f"  [155] L={L:.9f} N={len(mid)} nq={C['nq']} var(ds)={np.var(ds):.5f} "
          f"tutarlılık maks|g−(ds+1)2π/log(m/2π)|={tut:.2e} ({time.time()-t0:.0f}s)",
          flush=True)
    eta_yol = OUT / f"eta_{pencere}_t0.4_c4000.npz"

    # ---------------- 184b ----------------
    b184 = yukle("b184", QM / "184_configs" / "184b_K1_zarf.py")
    b184.ONK = b184.onkayit()                 # SHA denetimi orijinal yerde
    b184.S155 = OUT
    b184.S184 = OUT
    k1 = OUT / f"K1_{etiket}.npz"
    if not k1.exists():
        print(f"  [184b] kos({etiket}, {pencere}) …", flush=True)
        print("   ", b184.kos(etiket, pencere), flush=True)
    G = np.load(k1)
    w = np.asarray(G["w"], float)
    tau = np.asarray(G["tau"], float)
    aq = np.asarray(G["aq"], float)
    ae = np.asarray(G["aq_eff"], float)
    win = (tau >= 0.45) & (tau < 0.86)
    print(f"  [184b] çizgi evreni τ≤0.86: {len(w)}; pencere [0.45,0.86): "
          f"{int(win.sum())}; medyan w={float(np.median(G['wq'])):.4f}", flush=True)

    # ---------------- 185b / 186b / 187c modülleri ----------------
    b185 = yukle("b185", QM / "185_configs" / "185b_oz_muhasebe.py")
    b186 = yukle("b186", QM / "186_configs" / "186b_g1_oztutarlilik.py")
    b187 = yukle("b187", QM / "187_configs" / "187c_zeta_defteri.py")
    ONK185 = b185.onkayit()
    ONK186 = b186.onkayit()
    ONK187 = b187.onkayit()

    def kinematik_yamasi(modul):
        orijinal = modul.kinematik

        def kinematik(gaz):
            if gaz != etiket:
                return orijinal(gaz)
            # 185b 'gercek' kolu AYNEN (yalnız önbellek dosyası pencereninki)
            d = np.load(eta_yol)
            mid_ = np.asarray(d["mid"], float)
            ds_ = np.asarray(d["ds"], float)
            L_ = float(d["L"])
            g_ = (ds_ + 1.0) * TWO_PI / np.log(mid_ / TWO_PI)
            return g_, mid_, L_
        modul.kinematik = kinematik

    kinematik_yamasi(b185)
    kinematik_yamasi(b186.b185)
    b185.S185 = OUT

    oz = OUT / f"OZ_{etiket}.npz"
    if not oz.exists():
        print("  [185b] oz_hesapla …", flush=True)
        b185.oz_hesapla(etiket, w)
    pj = OUT / f"G1_proj_{etiket}.npz"
    if not pj.exists():
        print("  [186b] merdiven_izdusum (ρ≡1) …", flush=True)
        rhos = b186.rho_vektorleri(tau)
        b186.merdiven_izdusum(etiket, w, aq, {"bir": rhos["bir"]}, pj)
    OZ = np.load(oz)
    P = np.load(pj)
    if not (np.array_equal(OZ["nb"], G["nb"]) and np.array_equal(P["nb"], G["nb"])):
        raise SystemExit("BLOK KENARLARI 184 İLE UYUMSUZ")

    # ---------------- 187c ζ defteri (AYNEN) ----------------
    kenar_b = ONK186["kenar"]
    maskeler = [(tau >= kenar_b[b]) & (tau < kenar_b[b + 1]) for b in range(8)]
    maskeler.append((tau >= kenar_b[0]) & (tau < kenar_b[-1]))   # HAVUZ
    et = [f"{kenar_b[b]:.2f}-{kenar_b[b+1]:.2f}" for b in range(8)] + ["HAVUZ"]

    def zeta_defteri(disari):
        c_olc = b185.c_olculu(G, disari)
        c_oz = b185.c_oz(OZ, aq, disari)
        ck = b187.c_kesik(P, disari)
        delta = c_olc - ck
        mix = ck - c_oz
        out = np.zeros((9, 4))
        zc = np.zeros(9, complex)
        for k, m in enumerate(maskeler):
            zz = np.sum(delta[m] * np.conj(mix[m])) / \
                np.sum(np.abs(mix[m]) ** 2)
            sae = ae[m].sum()
            soz = np.abs(c_oz[m]).sum()
            m_olc = (np.abs(c_olc[m]).sum() - soz) / sae
            m_kes = (np.abs(ck[m]).sum() - soz) / sae
            out[k] = [np.abs(zz), np.degrees(np.angle(zz)), m_olc, m_kes]
            zc[k] = zz
        return out, zc

    tam, zc = zeta_defteri(-1)
    rr = [zeta_defteri(i) for i in range(NJACK)]
    reps = np.array([x[0] for x in rr])
    zreps = np.array([x[1] for x in rr])
    se_mod = b185.jk_se(reps[:, :, 0])
    dfark = b187.sar(reps[:, :, 1] - tam[None, :, 1])
    se_aci = np.sqrt((NJACK - 1) / NJACK * np.sum((dfark - dfark.mean(0)) ** 2, 0))
    genlik = 1.0 - tam[:, 2] / tam[:, 3]
    g_reps = 1.0 - reps[:, :, 2] / reps[:, :, 3]
    se_genlik = b185.jk_se(g_reps)
    icinde = np.abs(b187.sar(tam[:, 1] - 180.0)) <= 15.0

    # w-bant defteri (185b konvansiyonu; bilgi): w_g = Σâ/Σae, w^öz = Σ|c^öz|/Σae
    def wdef(disari):
        ag = np.abs(b185.c_olculu(G, disari))
        og = np.abs(b185.c_oz(OZ, aq, disari))
        return np.array([[ag[m].sum() / ae[m].sum(), og[m].sum() / ae[m].sum()]
                         for m in maskeler])
    wt = wdef(-1)
    wse = b185.jk_se(np.array([wdef(i) for i in range(NJACK)]))

    print(f"\n  ζ DEFTERİ (187c AYNEN) [{etiket}]  (±jk se):")
    for k in range(9):
        print(f"{et[k]:>10} |ζ|={tam[k,0]:.4f}±{se_mod[k]:.4f} "
              f"∠{tam[k,1]:+8.2f}±{se_aci[k]:5.2f} "
              f"{'EVET' if icinde[k] else 'HAYIR':>5} "
              f"genlik={genlik[k]:.4f}±{se_genlik[k]:.4f} "
              f"m_ölç={tam[k,2]:.4f} m^kesik={tam[k,3]:.4f}  "
              f"w_g={wt[k,0]:.4f}±{wse[k,0]:.4f} w^öz={wt[k,1]:.4f}", flush=True)

    sonuc = {"pencere": pencere, "etiket": etiket, "L": float(L),
             "N": int(G["N"]), "nline_tau086": int(len(w)),
             "pencere_cizgi_045_086": int(win.sum()),
             "t_lo": float(z[0]), "t_hi": float(z[-1]),
             "var_ds": float(np.var(ds)), "tutarlilik_g": tut,
             "medyan_w": float(np.median(G["wq"])),
             "sha_onkayit": {"184": b184.ONK["sha256"], "185": ONK185["sha256"],
                             "186": ONK186["sha256"], "187": ONK187["sha256"]},
             "etiket187": et,
             "zeta_mod": tam[:, 0].tolist(), "se_mod": se_mod.tolist(),
             "zeta_aci": tam[:, 1].tolist(), "se_aci": se_aci.tolist(),
             "zeta_re": zc.real.tolist(), "zeta_im": zc.imag.tolist(),
             "pencere_180pm15": icinde.tolist(),
             "genlik_okuma": genlik.tolist(), "se_genlik": se_genlik.tolist(),
             "m_olc": tam[:, 2].tolist(), "m_kesik": tam[:, 3].tolist(),
             "w_g": wt[:, 0].tolist(), "se_w_g": wse[:, 0].tolist(),
             "w_oz_g": wt[:, 1].tolist(), "se_w_oz_g": wse[:, 1].tolist()}

    # ---------------- MAKİNE MÜHRÜ (yalnız son) ----------------
    if pencere == "son":
        s = {}
        A, B = np.load(eta_yol), np.load(S155 / "eta_son_t0.4_c4000.npz")
        s["eta_son_(155)"] = all(esit(A[k], B[k]) for k in B.files)
        A, B = G, np.load(S184 / "K1_gercek.npz")
        s["K1_gercek_(184b)"] = all(esit(A[k], B[k]) for k in B.files)
        A, B = OZ, np.load(S185 / "OZ_gercek.npz")
        s["OZ_gercek_(185b)"] = all(esit(A[k], B[k]) for k in B.files)
        A, B = P, np.load(S186 / "G1_proj_gercek.npz")
        s["G1_proj_bir_(186b)"] = all(esit(A[k], B[k]) for k in
                                      ("re_bir", "im_bir", "d_ss", "dyy_bir",
                                       "dxy_bir", "drr_bir", "nb", "w"))
        ref = json.load(open(S187 / "K2_zeta.json"))["gercek"]
        s["K2_zeta_gercek_(187c)"] = (
            tam[:, 0].tolist() == ref["zeta_mod"]
            and se_mod.tolist() == ref["se_mod"]
            and tam[:, 1].tolist() == ref["zeta_aci"]
            and se_aci.tolist() == ref["se_aci"]
            and genlik.tolist() == ref["genlik_okuma"]
            and se_genlik.tolist() == ref["se_genlik"]
            and tam[:, 2].tolist() == ref["m_olc"]
            and tam[:, 3].tolist() == ref["m_kesik"])
        s["HEPSI"] = all(s.values())
        sonuc["makine_muhru"] = s
        print(f"\n  MAKİNE MÜHRÜ (son, aynı sarmalayıcı zinciri; bit-bit):")
        for k, v in s.items():
            print(f"    {k:28s} {'TUTTU' if v else 'TUTMADI'}", flush=True)
        print(f"    ζ_g(HAVUZ) = {tam[8,0]:.4f}±{se_mod[8]:.4f} "
              f"∠{tam[8,1]:+.2f}°±{se_aci[8]:.2f}  (187: 0.3287±0.0023 ∠179.96°)",
              flush=True)

    sonuc["sure_s"] = time.time() - t00
    json.dump(sonuc, open(OUT / f"zincir_{etiket}.json", "w"), indent=1,
              ensure_ascii=False)
    np.savez_compressed(OUT / f"zincir_{etiket}.npz", z_tam=tam, z_reps=reps,
                        zc=zc, zreps=zreps, se_mod=se_mod, se_aci=se_aci)
    print(f"\n  -> {OUT/('zincir_'+etiket+'.json')}  ({time.time()-t00:.0f}s)  BİTTİ",
          flush=True)


if __name__ == "__main__":
    zincir(sys.argv[1])
