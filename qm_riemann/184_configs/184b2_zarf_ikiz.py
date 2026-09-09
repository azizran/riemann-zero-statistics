# -*- coding: utf-8 -*-
"""
184b2 — K1 (ikizler): â_q PROFİLİ SADAKATLİ İKİZLERDE (Hkeskin, HA4, …)
======================================================================
164/180'in SADAKATLİ çözücüsüyle kurulan ikizler z-dosyadan okunur; mid,ds
DOĞRUDAN z'den hesaplanır (â_q ham `ds` izdüşümü olduğundan η-regresyonu
GEREKMEZ — 155/159 ile birebir: mid=½(z_i+z_{i+1}), ds=g·Lw/2π−1).

Çizgi kümesi ORTAK: son'un q,w,mq'su (K1_gercek.npz) — böylece r(τ) line-
by-line hizalı. Nominal a_q^eff = 2 a_q sin(πτ), τ = w/L (gazın kendi L'si).

Çıktı: 184/K1_<ad>.npz (184b ile aynı şema).
Kullanım: 184b2_zarf_ikiz.py <ad>=<zfile> [...]
   ör: 184b2_zarf_ikiz.py Hkeskin=.../164/z_Hkeskin.npy HA4=.../164/z_HA4.npy
"""
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S184 = SCR / "184"
TWO_PI = 2 * np.pi
BLOK = 2000
NJACK = 8


def izdusum_bloklu(ds, mid, w, njack=NJACK, blok=BLOK):
    N = len(mid)
    nline = len(w)
    kenar = np.linspace(0, N, njack + 1).astype(int)
    nb = np.diff(kenar)
    re_blok = np.zeros((njack, nline))
    im_blok = np.zeros((njack, nline))
    for b in range(njack):
        lo, hi = kenar[b], kenar[b + 1]
        for s0 in range(lo, hi, blok):
            sl = slice(s0, min(s0 + blok, hi))
            A = np.outer(mid[sl], w)
            re_blok[b] += np.cos(A).T @ ds[sl]
            im_blok[b] -= np.sin(A).T @ ds[sl]
            del A
    return re_blok, im_blok, nb, N


def kos(ad, zfile, q, w, mq):
    t0 = time.time()
    z = np.sort(np.load(zfile).astype(float))
    g = np.diff(z)
    mid = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(mid / TWO_PI)
    L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1.0
    tau = w / L
    aq = 1.0 / (np.pi * mq * np.sqrt(q))
    aq_eff = 2.0 * aq * np.sin(np.pi * tau)
    print(f"  [{ad}] z={Path(zfile).name} N={len(mid)} L={L:.9f} "
          f"σ_ds={np.std(ds):.5f}", flush=True)
    re_blok, im_blok, nb, N = izdusum_bloklu(ds, mid, w)
    c = 2.0 * (re_blok.sum(0) + 1j * im_blok.sum(0)) / N
    ahat = np.abs(c)
    wq = ahat / aq_eff
    np.savez_compressed(
        S184 / f"K1_{ad}.npz", q=q, w=w, tau=tau, mq=mq, aq=aq,
        aq_eff=aq_eff, ahat=ahat, wq=wq, re_blok=re_blok, im_blok=im_blok,
        nb=nb, N=N, L=L)
    print(f"      -> K1_{ad}.npz ({time.time()-t0:.0f}s) medyan w={np.median(wq):.4f}",
          flush=True)


if __name__ == "__main__":
    G = np.load(S184 / "K1_gercek.npz")
    q, w, mq = G["q"], G["w"], G["mq"]
    print(f"184b2 — ortak çizgi kümesi: son (nline={len(q)})", flush=True)
    for arg in sys.argv[1:]:
        ad, zf = arg.split("=", 1)
        kos(ad, zf, q, w, mq)
