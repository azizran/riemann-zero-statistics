# -*- coding: utf-8 -*-
"""
184b — K1: ZARF PROFİLİ w(τ) = â_q/a_q^eff  (gerçek 'son' + ikiz 'keskin')
==========================================================================
K0/ONKAYIT_184 dondurdu:  â_q = |c_q(ds)|, c_q(ds)=2⟨ds_n e^{−iω_q m_n}⟩,
a_q^eff = 2 a_q sin(πτ_q), w_q = â_q/a_q^eff.

Makine 155/159/162 ile BİT-BİT aynı ham izdüşüm: `ds`,`mid` doğrudan
eta_onbellek'ten (yeniden regresyon YOK); çizgi kümesi 162'nin evreni
(pk_m(e^{0.86L})). Tek geçişte 8 bitişik JACKKNIFE bloğunun kısmi
toplamları toplanır (leave-one-out sonradan aritmetik).

Çıktı: 184/K1_<etiket>.npz  (q,w,tau,mq,aq_eff,ahat,re_blok,im_blok,N,nb).
Kullanım: 184b_K1_zarf.py
"""
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "154_configs"))
import importlib
C154 = importlib.import_module("154_cekirdek")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S184 = SCR / "184"
S155 = SCR / "155"
TAU_CIZGI = 0.86
BLOK = 2000


def onkayit():
    p = S184 / "ONKAYIT_184.json"
    o = json.load(open(p))
    sha = hashlib.sha256((QM / "184_configs" / "184a_onkayit.py").read_bytes()).hexdigest()
    if sha != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {sha} != {o['sha256']}")
    return o


def izdusum_bloklu(ds, mid, w, njack, blok=BLOK):
    """Her çizgi için 8 bitişik blokta Σ ds cos(w·mid), Σ ds sin(w·mid)."""
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


def kos(etiket, gaz):
    t0 = time.time()
    d = np.load(S155 / f"eta_{gaz}_t0.4_c4000.npz")
    mid = np.asarray(d["mid"], float)
    ds = np.asarray(d["ds"], float)
    L = float(d["L"])
    qm = C154.pk_m(int(np.exp(TAU_CIZGI * L)))
    q = np.array(sorted(qm), dtype=float)
    mq = np.array([qm[int(x)] for x in q], dtype=float)
    w = np.log(q)
    tau = w / L
    aq = 1.0 / (np.pi * mq * np.sqrt(q))
    aq_eff = 2.0 * aq * np.sin(np.pi * tau)
    print(f"  [{etiket}={gaz}] L={L:.9f} N={len(mid)} nline={len(q)} "
          f"tau∈[{tau.min():.3f},{tau.max():.3f}]", flush=True)
    re_blok, im_blok, nb, N = izdusum_bloklu(ds, mid, w, ONK["njack"])
    re_tot = re_blok.sum(0)
    im_tot = im_blok.sum(0)
    c = 2.0 * (re_tot + 1j * im_tot) / N
    ahat = np.abs(c)
    wq = ahat / aq_eff
    np.savez_compressed(
        S184 / f"K1_{etiket}.npz", q=q, w=w, tau=tau, mq=mq, aq=aq,
        aq_eff=aq_eff, ahat=ahat, wq=wq, re_blok=re_blok, im_blok=im_blok,
        nb=nb, N=N, L=L)
    print(f"      -> K1_{etiket}.npz  ({time.time()-t0:.0f}s)  "
          f"medyan w={np.median(wq):.4f}", flush=True)
    return dict(etiket=etiket, gaz=gaz, L=L, N=int(N), nline=int(len(q)))


if __name__ == "__main__":
    S184.mkdir(parents=True, exist_ok=True)
    ONK = onkayit()
    print("=" * 70)
    print(f"184b / K1 ZARF PROFİLİ  [on-kayit {ONK['zaman']} sha {ONK['sha256'][:12]}]")
    print("=" * 70)
    ozet = []
    for et, gaz in ONK["gazlar"].items():
        ozet.append(kos(et, gaz))
    # ozet + L-tutarlılık
    Ls = {o["etiket"]: o["L"] for o in ozet}
    print(f"  L(gercek)={Ls.get('gercek'):.9f}  L(ikiz)={Ls.get('ikiz'):.9f}  "
          f"fark={abs(Ls['gercek']-Ls['ikiz']):.2e}", flush=True)
    (S184 / "K1_ozet.json").write_text(json.dumps(ozet, indent=1))
    print("  BİTTİ", flush=True)
