# -*- coding: utf-8 -*-
"""
187d — K1-kontrol: İKİZ SAĞIRLIĞI (Hkeskin) + HA4 defteri (K4 hazırlığı)
========================================================================
ONKAYIT_187 dondurdu — kontrol ızgarası [0.86, 0.90, 1.00, 1.10]:
  Hkeskin kinematiği: NOMİNAL katmanlar; c^kesik = G1_proj_Hkeskin
  re_bir/im_bir AYNEN. H-F3: (1.00,1.10] katmanında |Δm(HAVUZ)| ≤ 2se.
  Makine mührü: m^ya_Hk(1.00) ölçülü m_Hk'yi yeniden vermeli
  (ikizin kendi merdiveni τ≤1.00).
  HA4 kinematiği: env-ağırlıklı (0.5·erfc((τ−0.68)/0.125)) taban (pencere
  içi, 3425 çizgi) İLK KEZ HA4'ün kendi kinematiğinde + env-katmanlar;
  makine mührü: m^ya,env(1.10) ölçülü m_HA4'ü yeniden vermeli.

Çıktı: 187/katman_<i>_{Hkeskin,HA4env}.npy, katmanproj_{Hkeskin,HA4}.npz,
       K1k_ikiz_sagirlik.json
"""
import hashlib
import importlib.util
import json
import math
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S184, S185, S186, S187 = (SCR / d for d in ["184", "185", "186", "187"])
NJACK = 8
PBLOK = 2048
LBLOK = 4096
PROJ_BLOK = 1500

spec = importlib.util.spec_from_file_location(
    "b185", QM / "185_configs" / "185b_oz_muhasebe.py")
b185 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b185)

spec7 = importlib.util.spec_from_file_location(
    "b187", QM / "187_configs" / "187b_katman_defteri.py")
b187 = importlib.util.module_from_spec(spec7)
spec7.loader.exec_module(b187)   # asal_kuvvetler, katman_serisi, c_kesik


def env_f(t):
    return 0.5 * math.erfc((t - 0.68) / 0.125)


if __name__ == "__main__":
    ONK = b187.onkayit()
    print("=" * 78)
    print(f"187d / K1-KONTROL İKİZ SAĞIRLIĞI + HA4  [on-kayit sha "
          f"{ONK['sha256'][:12]}]")
    print("=" * 78, flush=True)
    L = ONK["L"]
    izg = ONK["katman_izgara_kontrol"]        # [0.86, 0.90, 1.00, 1.10]
    NKAT = len(izg) - 1                       # 3 katman

    G = np.load(S184 / "K1_gercek.npz")
    w_win = np.asarray(G["w"], float)
    tau = np.asarray(G["tau"], float)
    aq = np.asarray(G["aq"], float)
    ae = np.asarray(G["aq_eff"], float)
    env_win = np.array([env_f(t) for t in tau])

    qmax = np.exp(izg[-1] * L)
    qa, lam, taua, aqa = b187.asal_kuvvetler(qmax, L)
    katman_idx = [np.where((taua > izg[i]) & (taua <= izg[i + 1]))[0]
                  for i in range(NKAT)]
    print(f"\nKatman çizgi sayıları (q≤{qmax:.0f}):")
    for i in range(NKAT):
        print(f"  katman {i+1} ({izg[i]:.2f},{izg[i+1]:.2f}]: "
              f"{len(katman_idx[i])} çizgi", flush=True)

    kenar_b = json.load(open(S186 / "ONKAYIT_186.json"))["kenar"]
    maskeler = [(tau >= kenar_b[b]) & (tau < kenar_b[b + 1]) for b in range(8)]
    maskeler.append((tau >= kenar_b[0]) & (tau < kenar_b[-1]))   # HAVUZ
    et = [f"{kenar_b[b]:.2f}-{kenar_b[b+1]:.2f}" for b in range(8)] + ["HAVUZ"]
    sonuc = {"sha_onkayit": ONK["sha256"], "izgara": izg, "etiket": et,
             "katman_cizgi": [int(len(k)) for k in katman_idx]}

    for gaz, agirlik in [("Hkeskin", "bir"), ("HA4", "env")]:
        g, mid, Lg = b185.kinematik(gaz)
        N = len(mid)
        kenar_jk = np.linspace(0, N, NJACK + 1).astype(int)
        ds_ref = g * np.log(mid / (2 * np.pi)) / (2 * np.pi) - 1.0

        # katman ağırlıkları: nominal (bir) veya env·nominal (HA4)
        if agirlik == "bir":
            kat_a = [aqa[idx] for idx in katman_idx]
        else:
            kat_a = [aqa[idx] * np.array([env_f(t) for t in taua[idx]])
                     for idx in katman_idx]

        # --- AŞAMA A: katman serileri ---
        son_ek = "Hkeskin" if gaz == "Hkeskin" else "HA4env"
        for i in range(NKAT):
            yol = S187 / f"katman_{i+1}_{son_ek}.npy"
            if yol.exists():
                print(f"  [{gaz} katman {i+1}] mevcut, atlanıyor", flush=True)
                continue
            idx = katman_idx[i]
            dds = b187.katman_serisi(g, mid, np.log(qa[idx]), kat_a[i],
                                     f"{gaz} katman {i+1}")
            np.save(yol, dds)
            print(f"  -> {yol}", flush=True)

        # --- AŞAMA B: taban + izdüşümler ---
        yolB = S187 / f"katmanproj_{gaz}.npz"
        if not yolB.exists():
            t0 = time.time()
            katmanlar = [np.load(S187 / f"katman_{i+1}_{son_ek}.npy")
                         for i in range(NKAT)]
            taban = np.zeros(N)
            re = np.zeros((NKAT, NJACK, len(w_win)))
            im = np.zeros((NKAT, NJACK, len(w_win)))
            reT = np.zeros((NJACK, len(w_win)))   # taban izdüşümü (HA4 için)
            imT = np.zeros((NJACK, len(w_win)))
            kats_w = 2.0 * aq * (env_win if agirlik == "env" else 1.0)
            dyy = dxy = 0.0
            for b in range(NJACK):
                lo, hi = kenar_jk[b], kenar_jk[b + 1]
                for s0 in range(lo, hi, PROJ_BLOK):
                    sl = slice(s0, min(s0 + PROJ_BLOK, hi))
                    P = np.outer(mid[sl], w_win)
                    cP = np.cos(P)
                    sP = np.sin(P)
                    del P
                    Kw = np.sin(0.5 * np.outer(g[sl], w_win))
                    Kw *= cP
                    tb = Kw @ kats_w
                    taban[sl] = tb
                    dyy += (tb ** 2).sum()
                    dxy += (tb * ds_ref[sl]).sum()
                    reT[b] += tb @ cP
                    imT[b] -= tb @ sP
                    for i in range(NKAT):
                        v = katmanlar[i][sl]
                        re[i, b] += v @ cP
                        im[i, b] -= v @ sP
                    del cP, sP, Kw
                print(f"  [{gaz} izdüşüm] jk-blok {b+1}/8  "
                      f"({time.time()-t0:.0f}s)", flush=True)
            np.save(S187 / f"katman_0_{son_ek}.npy", taban)
            np.savez_compressed(yolB, re=re, im=im, reT=reT, imT=imT,
                                N=N, nb=np.diff(kenar_jk),
                                izgara=np.array(izg), dyy_taban=dyy,
                                dxy_taban=dxy, agirlik=agirlik)
            print(f"  -> {yolB}  ({time.time()-t0:.0f}s)", flush=True)
        else:
            print(f"  [{gaz} izdüşüm] mevcut, atlanıyor", flush=True)

        KP = np.load(yolB)
        D = np.load(S184 / f"K1_{gaz}.npz")
        OZ = np.load(S185 / f"OZ_{gaz}.npz")
        if gaz == "Hkeskin":
            Pk = np.load(S186 / "G1_proj_Hkeskin.npz")
            # MÜHÜR: yeniden hesaplanan taban momentleri 186 ile örtüşmeli
            f1 = abs(float(KP["dyy_taban"]) - float(Pk["dyy_bir"].sum())) / \
                float(Pk["dyy_bir"].sum())
            print(f"  MÜHÜR taban-momenti (Hk, ρ≡1) vs 186: "
                  f"|Δdyy|/dyy = {f1:.2e}", flush=True)
            sonuc["muhur_taban_Hk"] = f1

        def cplx(reb, imb, disari, N, nb):
            if disari < 0:
                return 2.0 * (reb.sum(0) + 1j * imb.sum(0)) / N
            return 2.0 * ((reb.sum(0) - reb[disari]) +
                          1j * (imb.sum(0) - imb[disari])) / (N - nb[disari])

        oz_w = env_win if agirlik == "env" else np.ones_like(tau)

        def defter(disari):
            c_olc = b185.c_olculu(D, disari)
            a_oz = oz_w * np.abs(b185.c_oz(OZ, aq, disari))
            if gaz == "Hkeskin":
                ck = b187.c_kesik(Pk, disari)
            else:
                ck = cplx(KP["reT"], KP["imT"], disari, int(KP["N"]),
                          KP["nb"])
            ckat = [cplx(KP["re"][i], KP["im"][i], disari, int(KP["N"]),
                         KP["nb"]) for i in range(NKAT)]
            out = np.zeros((9, 2 + NKAT))
            for k, m in enumerate(maskeler):
                sae = ae[m].sum()
                soz = a_oz[m].sum()
                out[k, 0] = (np.abs(c_olc[m]).sum() - soz) / sae
                cya = ck[m].copy()
                out[k, 1] = (np.abs(cya).sum() - soz) / sae
                for i in range(NKAT):
                    cya = cya + ckat[i][m]
                    out[k, 2 + i] = (np.abs(cya).sum() - soz) / sae
            return out

        tam = defter(-1)
        reps = np.array([defter(i) for i in range(NJACK)])
        se = b185.jk_se(reps)
        dm = np.diff(tam[:, 1:], axis=1)
        se_dm = b185.jk_se(np.diff(reps[:, :, 1:], axis=2))
        m_olc = tam[:, 0]

        print(f"\n[{gaz}] KATMAN DEFTERİ (ağırlık={agirlik}; ±jk se):")
        print(f"{'bant':>10} {'m_ölç':>8} | m^ya(0.86, 0.90, 1.00, 1.10)")
        for k in range(9):
            s = " ".join(f"{tam[k, 1+j]:.4f}" for j in range(NKAT + 1))
            print(f"{et[k]:>10} {m_olc[k]:8.4f} | {s}")
        print(f"[{gaz}] KATMAN ARTIMLARI (HAVUZ):")
        for i in range(NKAT):
            print(f"  ({izg[i]:.2f},{izg[i+1]:.2f}]: Δm = {dm[8, i]:+.5f} "
                  f"± {se_dm[8, i]:.5f}  |Δm|/2se = "
                  f"{abs(dm[8, i])/(2*se_dm[8, i]):.2f}")

        sonuc[gaz] = {
            "agirlik": agirlik, "m_olc": m_olc.tolist(),
            "se_m_olc": se[:, 0].tolist(),
            "m_ya": tam[:, 1:].tolist(), "se_m_ya": se[:, 1:].tolist(),
            "dm": dm.tolist(), "se_dm": se_dm.tolist()}

        if gaz == "Hkeskin":
            # makine mührü: m^ya(1.00) ölçülü m_Hk'yi vermeli
            rel = (tam[:, 3] - m_olc) / m_olc
            z_seal = (tam[:, 3] - m_olc) / np.sqrt(se[:, 3] ** 2 +
                                                   se[:, 0] ** 2)
            print(f"[{gaz}] MAKİNE MÜHRÜ m^ya(1.00) vs m_ölç: bağıl fark "
                  f"{np.round(rel[:8], 4).tolist()}  HAVUZ={rel[8]:+.4f}")
            sonuc["muhur_Hk_tau100"] = {"bagil_fark": rel.tolist(),
                                        "z": z_seal.tolist()}
            # H-F3: son katman (1.00,1.10]
            hf3_havuz = bool(abs(dm[8, -1]) <= 2 * se_dm[8, -1])
            hf3_bant = [bool(abs(dm[k, -1]) <= 2 * se_dm[k, -1])
                        for k in range(8)]
            print(f"[{gaz}] H-F3 (τ'>1.00 sağırlık, HAVUZ): "
                  f"|Δm|={abs(dm[8, -1]):.5f} ≤ 2se={2*se_dm[8, -1]:.5f} ? "
                  f"{'EVET — MÜHÜR' if hf3_havuz else 'HAYIR — ÖLDÜ'}   "
                  f"bant: {sum(hf3_bant)}/8")
            sonuc["HF3"] = {"havuz": hf3_havuz, "bant": hf3_bant,
                            "dm_son": dm[:, -1].tolist(),
                            "se_dm_son": se_dm[:, -1].tolist()}
        else:
            # makine mührü: m^ya,env(1.10) ölçülü m_HA4'ü vermeli
            rel = (tam[:, -1] - m_olc) / m_olc
            print(f"[{gaz}] MAKİNE MÜHRÜ m^ya,env(1.10) vs m_ölç: bağıl fark "
                  f"{np.round(rel[:8], 4).tolist()}  HAVUZ={rel[8]:+.4f}")
            sonuc["muhur_HA4_tau110"] = {"bagil_fark": rel.tolist()}

    json.dump(sonuc, open(S187 / "K1k_ikiz_sagirlik.json", "w"), indent=1,
              ensure_ascii=False)
    print(f"\n-> {S187/'K1k_ikiz_sagirlik.json'}  BİTTİ", flush=True)
