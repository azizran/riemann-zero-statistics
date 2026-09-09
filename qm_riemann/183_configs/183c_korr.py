# -*- coding: utf-8 -*-
"""
183c — korr(Xa,Xb) ÖLÇER (182g'nin İLK YARISI, AYNEN)
======================================================
ÖN-MÜHÜR: 183/ONKAYIT_183.json.

`182g_anatomi.py`'nin bölünmüş-merdiven zincirinin BİREBİR aynısı
(K.gaz → Model165(...,'olculen',τ_c=0.95) → .alanlar(kmax=3) →
τ-sıralı bir-atlamalı bölünme → K.sentez), yalnız pahalı tayf/ρ_J
bölümü YOKTUR. Ölçülenler: korr(Xa,Xb), κ_E/κ_Xa/κ_Xb, kurt(E),
ve aynı tohumla SUR-A surrogate kontrolü korr(Xa°,Xb°).

ÜREME KAPISI: `VF1` ile koşulduğunda 182/ANATOMI_VF1.json'un
korr_XaXb / kappa / kurt alanlarını 0 farkla üretmelidir.

Kullanım: 183c_korr.py VD0 VD1 VD2 VD4 [VF1 ...]
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("167_configs", "166_configs", "165_configs", "163_configs",
           "160_configs", "159_configs", "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
K = importlib.import_module("165_cekirdek")
ORT = importlib.import_module("167_ortak")
K.PENCERE.update(ORT.pencere_dict())

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S183 = SCR / "183"
KMAX = 3
G1 = np.sqrt(2 / np.pi)


def kap(F):
    Fc = float(np.std(F)) * np.sign(F)
    return float(np.dot(Fc, F) / np.dot(F, F))


def kurt(F):
    d = F - F.mean()
    s = d.std()
    return float((d ** 4).mean() / s ** 4 - 3.0)


def kos(ad, sur_tohum):
    t0 = time.time()
    print("=" * 78, flush=True)
    print(f"183c — korr(Xa,Xb): {ad}  (182g zinciri aynen)", flush=True)
    print("=" * 78, flush=True)
    Y = K.gaz(ad, 0.40, 4000)
    Mo = K.Model165(ad, 0.40, 4000, 0.95, "olculen", Y=Y, tau_c=0.95)
    Mo.sec("olculen", 0.95).alanlar(kmax=KMAX)
    idx = np.nonzero(Mo.msk)[0]
    o = np.argsort(Mo.M["tau"][idx])
    ia, ib = idx[o[0::2]], idx[o[1::2]]
    w = Mo.M["w"]
    Xa = K.sentez(Mo.s, w[ia], Mo.y[ia]); Xa -= Xa.mean()
    Xb = K.sentez(Mo.s, w[ib], Mo.y[ib]); Xb -= Xb.mean()
    E = Mo.E
    r = float(np.corrcoef(Xa, Xb)[0, 1])
    kE, kXa, kXb = kap(E), kap(Xa), kap(Xb)
    uE = kurt(E)

    rng = np.random.default_rng(1000 + int(sur_tohum))

    def sur(W, amp):
        ps = rng.uniform(0, 2 * np.pi, len(amp))
        f = K.sentez(Mo.s, W, np.abs(amp) * np.exp(1j * ps))
        return f - f.mean()

    Es = sur(w[Mo.msk], Mo.hp[Mo.msk])
    Xas = sur(w[ia], Mo.y[ia])
    Xbs = sur(w[ib], Mo.y[ib])
    rs = float(np.corrcoef(Xas, Xbs)[0, 1])
    # SUR-C sınavı: surrogate model alanı ölçülen alandan bağımsız mı?
    rEEs = float(np.corrcoef(Es, E)[0, 1])

    print(f"  |A|={len(ia)} |B|={len(ib)}  N={len(Mo.s)}", flush=True)
    print(f"  korr(Xa,Xb)  = {r:+.5f}", flush=True)
    print(f"  korr(Xa°,Xb°)= {rs:+.5f}   korr(E°,E) = {rEEs:+.5f}", flush=True)
    print(f"  κ_E={kE:.5f} κ_Xa={kXa:.5f} κ_Xb={kXb:.5f} (Gauss {G1:.5f})  "
          f"kurt(E)={uE:+.5f}", flush=True)

    rec = dict(ad=ad, korr_XaXb=r, korr_XaXb_sur=rs, korr_E_Esur=rEEs,
               kappa=dict(E=kE, Xa=kXa, Xb=kXb), kurt_E=uE,
               nA=int(len(ia)), nB=int(len(ib)), N=int(len(Mo.s)),
               sur_tohum=1000 + int(sur_tohum), sure_s=time.time() - t0)
    S183.mkdir(parents=True, exist_ok=True)
    p = S183 / f"KORR_{ad}.json"
    p.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=float))
    print(f"  -> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)
    return rec


if __name__ == "__main__":
    args = sys.argv[1:] or ["VD0", "VD1", "VD2", "VD4"]
    for a in args:
        # derinlik merdiveninin tamamı tohum 1'dir; VF<n>/VS<n> kendi tohumu
        th = 1 if a.startswith("VD") else int("".join(c for c in a
                                                      if c.isdigit()) or 1)
        kos(a, th)
