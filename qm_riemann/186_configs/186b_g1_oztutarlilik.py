# -*- coding: utf-8 -*-
"""
186b — K1: G1 ÖZ-TUTARLILIK (yeniden-ağırlıklı merdiven; ρ=184 yasası vs ρ≡1)
=============================================================================
ONKAYIT_186 dondurdu:
  ds^ya_n = Σ_q' 2 a_q' ρ(τ_q') sin(ω_q' g_n/2) cos(ω_q' m_n)
  c^ya_q  = 2⟨ds^ya_n e^{−iω_q m_n}⟩ ;  öz-terim = ρ_q·c_q^öz (KESİN)
  m^pred_b = [Σ|c^ya| − Σρ·â^öz]/Σae ;  M^pred_b = m^pred_b/m_Hk,b
BİRİNCİL kinematik GERÇEK; DONMUŞ YAN-LEHÇE: aynı koşular Hkeskin kinematiğinde.
Kontrol ρ≡1: zincir r_pred 8 bantta > 1 vermeli (makine mührü; vermezse ŞÜPHELİ).
ρ-varyantları: yasa (1−0.149τ^1.30), bir (≡1), env (0.5·erfc((τ−0.68)/0.125), K4 için).

Çıktı: 186/G1_proj_<gaz>.npz + 186/G1_sonuc.json + ekran defteri.
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
S184, S185, S186 = SCR / "184", SCR / "185", SCR / "186"
TWO_PI = 2 * np.pi
BLOK = 1500
NJACK = 8

spec = importlib.util.spec_from_file_location(
    "b185", QM / "185_configs" / "185b_oz_muhasebe.py")
b185 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b185)


def onkayit():
    o = json.load(open(S186 / "ONKAYIT_186.json"))
    sha = hashlib.sha256(
        (QM / "186_configs" / "186a_onkayit.py").read_bytes()).hexdigest()
    if sha != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {sha} != {o['sha256']}")
    return o


def rho_vektorleri(tau):
    yasa = 1.0 - 0.149 * tau ** 1.30
    bir = np.ones_like(tau)
    env = np.array([0.5 * math.erfc((t - 0.68) / 0.125) for t in tau])
    return {"yasa": yasa, "bir": bir, "env": env}


def merdiven_izdusum(gaz, w, aq, rhos, yol):
    """Tek geçiş: ds^ya (her ρ için) kur + tüm çizgilere izdüşür; jk-blok toplamları.
    Ayrıca ds^ya vs referans ds (g·L_n/2π − 1) sızıntı/kesme teşhis toplamları."""
    t0 = time.time()
    g, mid, L = b185.kinematik(gaz)
    N = len(mid)
    Ln = np.log(mid / TWO_PI)
    ds_ref = g * Ln / TWO_PI - 1.0
    kenar = np.linspace(0, N, NJACK + 1).astype(int)
    nb = np.diff(kenar)
    adlar = list(rhos)
    Q = len(w)
    re = {a: np.zeros((NJACK, Q)) for a in adlar}
    im = {a: np.zeros((NJACK, Q)) for a in adlar}
    d_ss = np.zeros(NJACK)                       # Σ ds_ref²
    d_yy = {a: np.zeros(NJACK) for a in adlar}   # Σ (ds^ya)²
    d_xy = {a: np.zeros(NJACK) for a in adlar}   # Σ ds^ya·ds_ref
    d_rr = {a: np.zeros(NJACK) for a in adlar}   # Σ (ds^ya−ds_ref)²
    katsayi = {a: 2.0 * aq * rhos[a] for a in adlar}
    for b in range(NJACK):
        lo, hi = kenar[b], kenar[b + 1]
        d_ss[b] = (ds_ref[lo:hi] ** 2).sum()
        for s0 in range(lo, hi, BLOK):
            sl = slice(s0, min(s0 + BLOK, hi))
            P = np.outer(mid[sl], w)
            cP = np.cos(P)
            sP = np.sin(P)
            del P
            K = np.sin(0.5 * np.outer(g[sl], w))
            K *= cP
            dref = ds_ref[sl]
            for a in adlar:
                dsv = K @ katsayi[a]
                re[a][b] += dsv @ cP
                im[a][b] -= dsv @ sP
                d_yy[a][b] += (dsv ** 2).sum()
                d_xy[a][b] += (dsv * dref).sum()
                d_rr[a][b] += ((dsv - dref) ** 2).sum()
            del cP, sP, K
        print(f"  [{gaz}] jk-blok {b+1}/8  ({time.time()-t0:.0f}s)", flush=True)
    kayit = dict(w=w, N=N, L=L, nb=nb, adlar=np.array(adlar), d_ss=d_ss)
    for a in adlar:
        kayit[f"re_{a}"] = re[a]
        kayit[f"im_{a}"] = im[a]
        kayit[f"dyy_{a}"] = d_yy[a]
        kayit[f"dxy_{a}"] = d_xy[a]
        kayit[f"drr_{a}"] = d_rr[a]
    np.savez_compressed(yol, **kayit)
    print(f"  [{gaz}] -> {yol}  ({time.time()-t0:.0f}s)", flush=True)


def c_ya(D, ad, disari=-1):
    N = int(D["N"])
    nb = D["nb"]
    if disari < 0:
        re, imv, n = D[f"re_{ad}"].sum(0), D[f"im_{ad}"].sum(0), N
    else:
        re = D[f"re_{ad}"].sum(0) - D[f"re_{ad}"][disari]
        imv = D[f"im_{ad}"].sum(0) - D[f"im_{ad}"][disari]
        n = N - nb[disari]
    return 2.0 * (re + 1j * imv) / n


if __name__ == "__main__":
    ONK = onkayit()
    print("=" * 78)
    print(f"186b / K1 G1 ÖZ-TUTARLILIK  [on-kayit {ONK['zaman']} "
          f"sha {ONK['sha256'][:12]}]")
    print("=" * 78, flush=True)

    G = np.load(S184 / "K1_gercek.npz")
    H = np.load(S184 / "K1_Hkeskin.npz")
    OZg = np.load(S185 / "OZ_gercek.npz")
    OZh = np.load(S185 / "OZ_Hkeskin.npz")
    w = np.asarray(G["w"], float)
    tau = np.asarray(G["tau"], float)
    aq = np.asarray(G["aq"], float)
    ae = np.asarray(G["aq_eff"], float)
    if not (np.allclose(H["aq_eff"], ae) and np.allclose(H["tau"], tau)):
        raise SystemExit("EVREN UYUMSUZ: Hkeskin çizgi evreni gerçekle özdeş değil")
    kenar = ONK["kenar"]
    maskeler = [(tau >= kenar[b]) & (tau < kenar[b + 1]) for b in range(8)]
    et = [f"{kenar[b]:.2f}-{kenar[b+1]:.2f}" for b in range(8)]
    rhos = rho_vektorleri(tau)

    # --- ağır geçişler (mevcutsa atla) ---
    for gaz in ["gercek", "Hkeskin"]:
        yol = S186 / f"G1_proj_{gaz}.npz"
        if not yol.exists():
            merdiven_izdusum(gaz, w, aq, rhos, yol)
        else:
            print(f"  [{gaz}] izdüşüm mevcut, atlanıyor", flush=True)
    Pg = np.load(S186 / "G1_proj_gercek.npz")
    Ph = np.load(S186 / "G1_proj_Hkeskin.npz")
    if not (np.array_equal(Pg["nb"], G["nb"]) and np.array_equal(Ph["nb"], H["nb"])):
        raise SystemExit("JK BLOK KENARLARI 184 İLE UYUMSUZ")

    # --- defter: tam + loo ---
    def defter(disari):
        ag = np.abs(b185.c_olculu(G, disari))
        ah = np.abs(b185.c_olculu(H, disari))
        og = np.abs(b185.c_oz(OZg, aq, disari))
        oh = np.abs(b185.c_oz(OZh, aq, disari))
        ya = {(gaz, a): np.abs(c_ya(P, a, disari))
              for gaz, P in [("g", Pg), ("h", Ph)] for a in ["yasa", "bir"]}
        out = []
        for m in maskeler:
            sae = ae[m].sum()
            wg, wh = ag[m].sum() / sae, ah[m].sum() / sae
            wog, woh = og[m].sum() / sae, oh[m].sum() / sae
            mg, mh = wg - wog, wh - woh
            satir = [wg, wh, wog, woh, mg, mh, mg / mh]
            # birincil (gerçek kinematik): m^pred = w^ya − Σρ â^öz_g/Σae
            for a in ["yasa", "bir"]:
                wya = ya[("g", a)][m].sum() / sae
                woz_ya = (rhos[a][m] * og[m]).sum() / sae
                mp = wya - woz_ya
                satir += [wya, mp, mp / mh]
            # yan-lehçe (Hkeskin kinematiği): öz-terim ρ·â^öz_Hk
            for a in ["yasa", "bir"]:
                wya = ya[("h", a)][m].sum() / sae
                woz_ya = (rhos[a][m] * oh[m]).sum() / sae
                mp = wya - woz_ya
                satir += [wya, mp, mp / mh]
            out.append(satir)
        return np.array(out)   # (8, 19)

    tam = defter(-1)
    reps = np.array([defter(i) for i in range(NJACK)])
    se = b185.jk_se(reps)

    K1f = np.load(S185 / "K1_faktorler.npz", allow_pickle=True)
    tau_bar = K1f["tau_bar"][:8]
    se185 = K1f["se_tablo"]  # (9,11): 7..10 = wg, wh, wog, woh jk-se
    wg_b, wh_b, wog_b, woh_b, mg_b, mh_b, M_b = [tam[:, j] for j in range(7)]

    # BİRİNCİL σ(M_ölç) — 185f defter konvansiyonu (bağımsız yayılım):
    se_mg = np.sqrt(se185[:8, 7] ** 2 + se185[:8, 9] ** 2)
    se_mh = np.sqrt(se185[:8, 8] ** 2 + se185[:8, 10] ** 2)
    sigM = M_b * np.sqrt((se_mg / mg_b) ** 2 + (se_mh / mh_b) ** 2)

    print("\nÖLÇÜLÜ M DEFTERİ (185f yeniden üretimi — kontrol):")
    for k in range(8):
        print(f"  τ̄={tau_bar[k]:.3f}  m_g={mg_b[k]:.4f} m_Hk={mh_b[k]:.4f} "
              f"M={M_b[k]:.4f}±{sigM[k]:.4f}")

    sut = {"g_yasa": (7, 8, 9), "g_bir": (10, 11, 12),
           "h_yasa": (13, 14, 15), "h_bir": (16, 17, 18)}
    sonuc = {"sha_onkayit": ONK["sha256"], "tau_bar": tau_bar.tolist(),
             "M_olc": M_b.tolist(), "sigM_defter": sigM.tolist(),
             "m_g": mg_b.tolist(), "m_Hk": mh_b.tolist(),
             "wog": wog_b.tolist(), "woh": woh_b.tolist()}

    def zincir(Mp):
        return (wog_b + Mp * mh_b) / (woh_b + mh_b)

    for ad, (jw, jm, jM) in sut.items():
        Mp = tam[:, jM]
        seMp = se[:, jM]
        z = (Mp - M_b) / sigM
        chi = float(np.sum(z ** 2) / 8)
        rp = zincir(Mp)
        sonuc[ad] = {"w_ya": tam[:, jw].tolist(), "m_pred": tam[:, jm].tolist(),
                     "M_pred": Mp.tolist(), "se_M_pred": seMp.tolist(),
                     "z": z.tolist(), "chi2dof_M": chi, "r_pred": rp.tolist()}
        print(f"\n[{ad}] M^pred defteri (σ birincil = defter-seM):")
        for k in range(8):
            print(f"  τ̄={tau_bar[k]:.3f}  M^pred={Mp[k]:.4f}±{seMp[k]:.4f}  "
                  f"M_ölç={M_b[k]:.4f}±{sigM[k]:.4f}  z={z[k]:+.1f}  "
                  f"r_pred={rp[k]:.4f}")
        print(f"  bant-χ²/dof = {chi:.2f}")

    # --- makine mührü: birincil kontrol (g_bir) zinciri 8 bantta > 1 mi? ---
    rp_bir = np.array(sonuc["g_bir"]["r_pred"])
    muhur_g = bool(np.all(rp_bir > 1.0))
    rp_bir_h = np.array(sonuc["h_bir"]["r_pred"])
    muhur_h = bool(np.all(rp_bir_h > 1.0))
    print(f"\nMAKİNE MÜHRÜ (birincil, gerçek-kinematik ρ≡1): r_pred>1 tüm bantlar? "
          f"{'EVET — 185-B ters-yönü yeniden üretildi' if muhur_g else 'HAYIR — makine ŞÜPHELİ'}")
    print(f"  r_pred(ρ≡1) = {np.round(rp_bir, 4).tolist()}")
    print(f"MAKİNE MÜHRÜ (yan-lehçe, Hk-kinematik ρ≡1): r_pred>1 tüm bantlar? "
          f"{'EVET' if muhur_h else 'HAYIR'}")
    print(f"  r_pred(ρ≡1) = {np.round(rp_bir_h, 4).tolist()}")

    # --- sızıntı/kesme teşhisi: ds^ya vs referans ds ---
    print("\nSIZINTI/KESME TEŞHİSİ (tam-örneklem):")
    tes = {}
    for gaz, P in [("gercek", Pg), ("Hkeskin", Ph)]:
        ss = P["d_ss"].sum()
        tes[gaz] = {}
        for a in ["yasa", "bir", "env"]:
            yy = P[f"dyy_{a}"].sum()
            xy = P[f"dxy_{a}"].sum()
            rr = P[f"drr_{a}"].sum()
            kor = xy / np.sqrt(ss * yy)
            tes[gaz][a] = {"korelasyon": float(kor),
                           "var_orani": float(yy / ss),
                           "artik_var_pay": float(rr / ss)}
            print(f"  [{gaz} ρ={a}]  korr(ds^ya,ds)={kor:.4f}  "
                  f"Var(ds^ya)/Var(ds)={yy/ss:.4f}  "
                  f"Var(ds^ya−ds)/Var(ds)={rr/ss:.4f}")
    sonuc["makine_muhru"] = {"g_bir_rpred_hepsi_bir_ustu": muhur_g,
                             "h_bir_rpred_hepsi_bir_ustu": muhur_h}
    sonuc["sizinti_teshisi"] = tes

    # pencere-içi karışım payı (bilgi): m^trunc(ρ≡1)/m_ölçülü
    mtr_g = np.array(sonuc["g_bir"]["m_pred"]) / mg_b
    mtr_h = np.array(sonuc["h_bir"]["m_pred"]) / mh_b
    print("\nPENCERE-İÇİ KARIŞIM PAYI m^ya(ρ≡1)/m_ölç (bilgi):")
    print(f"  gerçek : {np.round(mtr_g, 3).tolist()}")
    print(f"  Hkeskin: {np.round(mtr_h, 3).tolist()}")
    sonuc["pencere_ici_pay"] = {"gercek": mtr_g.tolist(), "Hkeskin": mtr_h.tolist()}

    # --- G1 hükmü (donmuş eşik: birincil lehçe χ²/dof ≤ 2) ---
    chi_g1 = sonuc["g_yasa"]["chi2dof_M"]
    hukum = "MÜHÜR" if chi_g1 <= 2.0 else "ÖLDÜ"
    if not muhur_g:
        hukum += " (makine ŞÜPHELİ: ρ≡1 mührü vermedi)"
    sonuc["G1_hukum"] = hukum
    sonuc["G1_chi2dof"] = chi_g1
    print(f"\nG1 HÜKMÜ (birincil = gerçek-kinematik, ρ=yasa): {hukum}  "
          f"[χ²/dof = {chi_g1:.2f}, eşik 2]")

    json.dump(sonuc, open(S186 / "G1_sonuc.json", "w"), indent=1,
              ensure_ascii=False)
    print(f"\n-> {S186/'G1_sonuc.json'}  BİTTİ", flush=True)
