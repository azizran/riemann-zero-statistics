# -*- coding: utf-8 -*-
"""
185b — K1: ÖZ-MUHASEBE (â^öz kesin beklenen-değer; üç faktörün bant defteri)
============================================================================
ONKAYIT_185 dondurdu:
  c_q^öz = (4 a_q / N) Σ_n sin(w_q g_n/2) cos(w_q m_n) e^{-i w_q m_n}
         = (2 a_q / N) (S_dc + S_2w),   S_dc = Σ sin(w_q g_n/2),
           S_2w = Σ sin(w_q g_n/2) e^{-2 i w_q m_n}
  Ĝ_q   = <(g_n - ḡ) e^{-i w_q m_n}>       (+ 2w kopyası anatomi için)

Kinematik kaynaklar (donmuş):
  gercek : mid,ds eta_son önbelleğinden; g = (ds+1)·2π/log(mid/2π)  (kesin inv.)
  Hkeskin: z_Hkeskin.npy → g=diff(z), m=(z_i+z_{i+1})/2
Ölçülü â'lar 184 K1 npz'lerinden AYNEN (yeniden ölçüm YOK).

Çıktı: 185/OZ_<gaz>.npz (blok toplamları) + 185/K1_faktorler.npz + ekran defteri.
"""
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S155, S184, S185 = SCR / "155", SCR / "184", SCR / "185"
TWO_PI = 2 * np.pi
BLOK = 2000
NJACK = 8


def onkayit():
    o = json.load(open(S185 / "ONKAYIT_185.json"))
    sha = hashlib.sha256((QM / "185_configs" / "185a_onkayit.py").read_bytes()).hexdigest()
    if sha != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {sha} != {o['sha256']}")
    return o


def kinematik(gaz):
    """Donmuş kaynaklardan (g, m, L) üçlüsü."""
    if gaz == "gercek":
        d = np.load(S155 / "eta_son_t0.4_c4000.npz")
        mid = np.asarray(d["mid"], float)
        ds = np.asarray(d["ds"], float)
        L = float(d["L"])
        g = (ds + 1.0) * TWO_PI / np.log(mid / TWO_PI)
        return g, mid, L
    zdos = {"Hkeskin": "z_Hkeskin.npy", "HA4": "z_HA4.npy"}[gaz]
    z = np.sort(np.load(S155 / zdos).astype(float))
    g = np.diff(z)
    mid = 0.5 * (z[:-1] + z[1:])
    L = float(np.log(mid / TWO_PI).mean())
    return g, mid, L


def oz_hesapla(gaz, w, njack=NJACK, blok=BLOK):
    """Tek geçişte blok toplamları: S_dc, S_2w, Ĝ (w ve 2w). KESİN — Taylor yok."""
    t0 = time.time()
    g, mid, L = kinematik(gaz)
    N = len(mid)
    gbar = g.mean()
    dg = g - gbar
    nline = len(w)
    kenar = np.linspace(0, N, njack + 1).astype(int)
    nb = np.diff(kenar)
    Sdc = np.zeros((njack, nline))
    S2re = np.zeros((njack, nline))
    S2im = np.zeros((njack, nline))
    Gre = np.zeros((njack, nline))
    Gim = np.zeros((njack, nline))
    G2re = np.zeros((njack, nline))
    G2im = np.zeros((njack, nline))
    sdg = np.zeros(njack)
    sdg2 = np.zeros(njack)
    for b in range(njack):
        lo, hi = kenar[b], kenar[b + 1]
        sdg[b] = dg[lo:hi].sum()
        sdg2[b] = (dg[lo:hi] ** 2).sum()
        for s0 in range(lo, hi, blok):
            sl = slice(s0, min(s0 + blok, hi))
            P = np.outer(mid[sl], w)
            cP = np.cos(P)
            sP = np.sin(P)
            del P
            Gm = np.sin(0.5 * np.outer(g[sl], w))
            c2 = 2.0 * cP * cP - 1.0
            s2 = 2.0 * sP * cP
            Sdc[b] += Gm.sum(0)
            S2re[b] += (Gm * c2).sum(0)
            S2im[b] -= (Gm * s2).sum(0)
            d = dg[sl]
            Gre[b] += d @ cP
            Gim[b] -= d @ sP
            G2re[b] += d @ c2
            G2im[b] -= d @ s2
            del cP, sP, Gm, c2, s2
    np.savez_compressed(
        S185 / f"OZ_{gaz}.npz", w=w, L=L, N=N, gbar=gbar, nb=nb,
        Sdc=Sdc, S2re=S2re, S2im=S2im, Gre=Gre, Gim=Gim,
        G2re=G2re, G2im=G2im, sdg=sdg, sdg2=sdg2)
    print(f"  [OZ {gaz}] N={N} L={L:.9f} ḡ={gbar:.7f} "
          f"σ_ε={np.sqrt(sdg2.sum()/N)/gbar:.5f}  ({time.time()-t0:.0f}s)", flush=True)
    return dict(gaz=gaz, N=N, L=L, gbar=gbar)


def c_oz(OZ, aq, disari=-1):
    """Kompleks c^öz (çizgi başına); disari>=0 ise o blok dışarıda (loo)."""
    N = int(OZ["N"])
    nb = OZ["nb"]
    if disari < 0:
        S = (OZ["Sdc"].sum(0) + OZ["S2re"].sum(0) + 1j * OZ["S2im"].sum(0))
        n = N
    else:
        S = ((OZ["Sdc"].sum(0) - OZ["Sdc"][disari])
             + (OZ["S2re"].sum(0) - OZ["S2re"][disari])
             + 1j * (OZ["S2im"].sum(0) - OZ["S2im"][disari]))
        n = N - nb[disari]
    return 2.0 * aq * S / n


def c_olculu(D, disari=-1):
    """184 npz bloklarından kompleks ölçülü c (loo destekli)."""
    N = int(D["N"])
    nb = D["nb"]
    if disari < 0:
        re, im, n = D["re_blok"].sum(0), D["im_blok"].sum(0), N
    else:
        re = D["re_blok"].sum(0) - D["re_blok"][disari]
        im = D["im_blok"].sum(0) - D["im_blok"][disari]
        n = N - nb[disari]
    return 2.0 * (re + 1j * im) / n


def jk_se(reps):
    reps = np.asarray(reps)
    return np.sqrt((NJACK - 1) / NJACK * np.sum((reps - reps.mean(0)) ** 2, 0))


if __name__ == "__main__":
    ONK = onkayit()
    print("=" * 78)
    print(f"185b / K1 ÖZ-MUHASEBE  [on-kayit {ONK['zaman']} sha {ONK['sha256'][:12]}]")
    print("=" * 78, flush=True)

    G = np.load(S184 / "K1_gercek.npz")
    H = np.load(S184 / "K1_Hkeskin.npz")
    q, w, tau, ae, aq = G["q"], G["w"], G["tau"], G["aq_eff"], G["aq"]
    kenar = ONK["kenar"]
    maskeler = [(tau >= kenar[b]) & (tau < kenar[b + 1]) for b in range(8)]
    maskeler.append(tau > ONK["kuyruk_tau"])
    et = [f"{kenar[b]:.2f}-{kenar[b+1]:.2f}" for b in range(8)] + ["KUYRUK>0.70"]

    # --- öz hesapları (kesin beklenen-değer) ---
    for gaz in ["gercek", "Hkeskin"]:
        if not (S185 / f"OZ_{gaz}.npz").exists():
            oz_hesapla(gaz, np.asarray(w, float))
        else:
            print(f"  [OZ {gaz}] mevcut, atlanıyor", flush=True)
    OZg = np.load(S185 / "OZ_gercek.npz")
    OZh = np.load(S185 / "OZ_Hkeskin.npz")
    if not (np.array_equal(OZg["nb"], G["nb"]) and np.array_equal(OZh["nb"], H["nb"])):
        raise SystemExit("BLOK KENARLARI 184 İLE UYUMSUZ")

    # --- tam-örneklem ve loo replikaları ---
    def defter(disari):
        ag = np.abs(c_olculu(G, disari))
        ah = np.abs(c_olculu(H, disari))
        og = np.abs(c_oz(OZg, aq, disari))
        oh = np.abs(c_oz(OZh, aq, disari))
        satir = {}
        for k, m in enumerate(maskeler):
            r = ag[m].sum() / ah[m].sum()
            F1 = ag[m].sum() / og[m].sum()
            F2 = og[m].sum() / oh[m].sum()
            F3i = oh[m].sum() / ah[m].sum()
            # çizgi lehçesi (H-W1b)
            F1c = np.sum(ae[m] * (ag[m] / og[m])) / ae[m].sum()
            F2c = np.sum(ae[m] * (og[m] / oh[m])) / ae[m].sum()
            F3ic = np.sum(ae[m] * (oh[m] / ah[m])) / ae[m].sum()
            wg = ag[m].sum() / ae[m].sum()
            wh = ah[m].sum() / ae[m].sum()
            wog = og[m].sum() / ae[m].sum()
            woh = oh[m].sum() / ae[m].sum()
            satir[k] = (r, F1, F2, F3i, F1c, F2c, F3ic, wg, wh, wog, woh)
        return np.array([satir[k] for k in range(len(maskeler))])  # (9, 11)

    tam = defter(-1)
    reps = np.array([defter(i) for i in range(NJACK)])   # (8, 9, 11)
    se = jk_se(reps)

    (r_b, F1_b, F2_b, F3i_b, F1c_b, F2c_b, F3ic_b,
     wg_b, wh_b, wog_b, woh_b) = [tam[:, j] for j in range(11)]

    # birincil sigma (184 defter konvansiyonu: bağımsız yayılım w_g, w_Hk)
    sig_r = r_b * np.sqrt((se[:, 7] / wg_b) ** 2 + (se[:, 8] / wh_b) ** 2)
    sig_r_ortak = se[:, 0]  # yan sütun: korelasyonlu ortak-jk

    # H-W1a: yapısal kapanış (teleskopik)
    kap_a = np.abs(F1_b * F2_b * F3i_b - r_b)
    # H-W1b: çizgi-lehçesi rekonstrüksiyonu
    r_rekon = F1c_b * F2c_b * F3ic_b
    z_b = (r_rekon[:8] - r_b[:8]) / sig_r[:8]
    chi2_dof_w1b = float(np.sum(z_b ** 2) / 8)

    tau_bar = np.array([np.sum(ae[m] * tau[m]) / ae[m].sum() for m in maskeler])

    print("\nK1 BANT DEFTERİ (toplam-oranı lehçesi; ±jk se):")
    print(f"{'bant':>12} {'τ̄':>6} {'r=g/Hk':>14} {'F_ANOM':>14} "
          f"{'F_KİN':>14} {'F_KOMŞU⁻¹':>14}")
    for k in range(9):
        print(f"{et[k]:>12} {tau_bar[k]:6.3f} "
              f"{r_b[k]:7.4f}±{sig_r[k]:.4f} {F1_b[k]:7.4f}±{se[k,1]:.4f} "
              f"{F2_b[k]:7.4f}±{se[k,2]:.4f} {F3i_b[k]:7.4f}±{se[k,3]:.4f}")

    print("\nÇİZGİ LEHÇESİ (H-W1b; ae-ağırlıklı çizgi-faktör ortalamaları):")
    print(f"{'bant':>12} {'F̄1':>14} {'F̄2':>14} {'F̄3⁻¹':>14} "
          f"{'r̂=ΠF̄':>8} {'r':>8} {'z':>6}")
    for k in range(8):
        print(f"{et[k]:>12} {F1c_b[k]:7.4f}±{se[k,4]:.4f} "
              f"{F2c_b[k]:7.4f}±{se[k,5]:.4f} {F3ic_b[k]:7.4f}±{se[k,6]:.4f} "
              f"{r_rekon[k]:8.4f} {r_b[k]:8.4f} {z_b[k]:6.2f}")
    print(f"\nH-W1a yapısal kapanış: maks|ΠF−r| = {kap_a.max():.2e}  (beklenti ~1e-15)")
    print(f"H-W1b bant-χ²/dof = {chi2_dof_w1b:.3f}  (ölüm eşiği > 2)")

    # TUTARLILIK: F_KOMŞU (= â_Hk/â_Hk^öz) 1'e ne kadar yakın?
    FK = 1.0 / F3i_b
    print("\nTUTARLILIK — F_KOMŞU = â_Hk/â_Hk^öz (ikiz ds'i TAM nominal merdiven):")
    for k in range(9):
        print(f"{et[k]:>12}  F_KOMŞU = {FK[k]:7.4f}   (1'den sapma {FK[k]-1:+.4f})")

    # w^öz profilleri (bilgi)
    print("\nw^öz = â^öz/ae profilleri (bilgi):")
    for k in range(9):
        print(f"{et[k]:>12}  w_g={wg_b[k]:6.4f} w_g^öz={wog_b[k]:6.4f}  "
              f"w_Hk={wh_b[k]:6.4f} w_Hk^öz={woh_b[k]:6.4f}")

    np.savez_compressed(
        S185 / "K1_faktorler.npz",
        etiket=np.array(et), tau_bar=tau_bar,
        r=r_b, sig_r=sig_r, sig_r_ortak=sig_r_ortak,
        F1=F1_b, F2=F2_b, F3i=F3i_b,
        F1c=F1c_b, F2c=F2c_b, F3ic=F3ic_b,
        se_tablo=se, r_rekon=r_rekon, z_w1b=z_b,
        chi2_dof_w1b=chi2_dof_w1b, kap_a=kap_a,
        wg=wg_b, wh=wh_b, wog=wog_b, woh=woh_b,
        kenar=np.array(kenar), kuyruk_tau=ONK["kuyruk_tau"])
    print("\n-> 185/K1_faktorler.npz  BİTTİ", flush=True)
