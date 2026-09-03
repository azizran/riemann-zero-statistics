"""
168 — A2(iii): κ'nın TEPE PROFİLİ, MERKEZ FAZI ÇIKARILMIŞ
=========================================================
168_rezonans'ın (P1)'i ham Re κ ile ölçüldüğü için 165 §2b'nin MERKEZ FAZINA
(arg κ ≈ ν·s̄, s̄ ≈ 1.05e6 ⇒ ν adımı dres/6 = 6.7e−6 iken faz > 2π dönüyor)
takıldı. Burada 165 §4c'nin ayrıştırması kullanılır:

    κ(ν) = e^{iν s̄} · κ̃(ν)      (κ̃ = zarf)

ve profil |κ̃| ile ölçülür. Ölçülenler, her seçili merdiven çizgisi r için:
    w_eff = ∫|κ̃(ω_r+x)| dx / |κ̃(ω_r)|        (tepe-integrali / tepe-yüksekliği)
    ve karşılaştırma:  sinc ailesi ⇒ w_eff = 2π/T = dres  (∫sinc = ∫sinc² = π)

Ayrıca Δ ≈ 0 çevresindeki Dirichlet zarfı (165 §4d'nin tablosu) yeniden
ölçülür — normalizasyon çapası.

Kullanım: 168_profil.py <gaz> [taban] [tau_c]
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("167_configs", "165_configs", "163_configs", "160_configs",
           "159_configs", "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
K = importlib.import_module("165_cekirdek")
ORT = importlib.import_module("167_ortak")
K.PENCERE.update(ORT.pencere_dict())
TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/168")


def kos(veri, taban=0.40, tau_c=0.95):
    t0 = time.time()
    SCR.mkdir(parents=True, exist_ok=True)
    Y = K.gaz(veri, taban, 4000)
    Mo = K.Model165(veri, taban, 4000, 0.95, "olculen", Y=Y, tau_c=tau_c)
    Mo.sec("olculen", tau_c)
    M = Mo.M
    dres, sbar = Mo.dres, Mo.sbar
    w = M["w"][Mo.msk]
    tau_r = M["tau"][Mo.msk]
    a_r = M["a"][Mo.msk]
    q_r = M["q"][Mo.msk]
    G = -np.pi * tau_r * a_r * np.cos(np.pi * tau_r)
    print(f"=== 168-PROFİL {veri} === N={len(Mo.s)} T={Mo.T:.1f} "
          f"dres={dres:.4e} s̄={sbar:.6g} çizgi={len(w)} "
          f"({time.time()-t0:.0f}s)", flush=True)

    # --- Δ ≈ 0 çapası (165 §4d) ---------------------------------------
    x0 = np.array([0., .25, .5, .75, 1., 1.5, 2., 5., 10.]) * dres
    k0 = K.kappa(Mo.s, x0)
    print("\n  Δ/dres :  " + "  ".join(f"{x/dres:6.2f}" for x in x0))
    print("  |κ̃|    :  " + "  ".join(f"{abs(v):6.4f}" for v in k0))
    print("  Dirichlet(düzgün) : 1.0000 0.9003 0.6366 0.3001 0.0000 ...")

    # --- ω_r çevresinde ince profil, merkez fazı çıkarılmış ------------
    sel = [0, 1, 2, 4, 8, 20, 60, 300, 1500, 3000, 6000, 8000]
    sel = [i for i in sel if i < len(w)]
    JJ = np.arange(-40, 41)
    step = dres / 8.0
    nu = np.concatenate([w[i] + JJ * step for i in sel])
    tm = time.time()
    kp = K.kappa(Mo.s, nu)
    print(f"\n  κ ölçüldü: {len(nu)} frekans ({time.time()-tm:.0f}s)",
          flush=True)
    out = []
    print("\n      q      τ_r      𝒢_r       |κ̃(ω_r)|   Reκ̃(ω_r)   "
          "∫|κ̃|dx/|κ̃(0)|   ∫κ̃dx/κ̃(0)   |κ̃|(±dres)/|κ̃|(0)   "
          "|κ̃|(±3dres)/(0)")
    for k, i in enumerate(sel):
        z = kp[k * len(JJ):(k + 1) * len(JJ)]
        zt = z * np.exp(-1j * nu[k * len(JJ):(k + 1) * len(JJ)] * sbar)
        c = len(JJ) // 2
        A = np.abs(z)
        w_abs = float(A.sum() * step / A[c]) if A[c] else float("nan")
        # faz-referanslı (kompleks) integral: κ̃'nın kendi merkez fazına göre
        ph = zt[c] / abs(zt[c]) if abs(zt[c]) else 1.0
        w_coh = float(np.real(zt / ph).sum() * step / abs(zt[c])) \
            if abs(zt[c]) else float("nan")
        r1 = float(0.5 * (A[c + 8] + A[c - 8]) / A[c]) if A[c] else np.nan
        r3 = float(0.5 * (A[c + 24] + A[c - 24]) / A[c]) if A[c] else np.nan
        out.append(dict(q=int(q_r[i]), tau=float(tau_r[i]), G=float(G[i]),
                        kabs=float(A[c]), kre=float(zt[c].real),
                        w_abs=w_abs / dres, w_coh=w_coh / dres, r1=r1, r3=r3,
                        prof=[float(x) for x in A / A[c]]))
        print(f"  {int(q_r[i]):8d} {tau_r[i]:.4f} {G[i]:+.6f}  {A[c]:.6f}  "
              f"{zt[c].real:+.6f}   {w_abs/dres:9.4f}      "
              f"{w_coh/dres:9.4f}    {r1:9.4f}          {r3:9.4f}",
              flush=True)

    # --- MERDİVEN ÇÖZÜNÜRLÜĞÜ: komşu çizgi aralığı / dres --------------
    ws = np.sort(w)
    d = np.diff(ws)
    print("\n  MERDİVEN ARALIĞI (δω) / dres, τ kutularında  "
          "[<1 ⇒ çizgiler pencerede AYRIŞMIYOR]")
    print("   τ kutu    n_çizgi   ⟨δω⟩/dres   medyan   %(δω<dres)")
    tm2 = tau_r[np.argsort(w)][1:]
    kut = []
    for lo in np.arange(0.1, 0.95, 0.1):
        s = (tm2 >= lo) & (tm2 < lo + 0.1)
        if s.sum() < 2:
            continue
        v = d[s] / dres
        kut.append(dict(lo=float(lo), n=int(s.sum()), ort=float(v.mean()),
                        med=float(np.median(v)), frac=float((v < 1).mean())))
        print(f"   {lo:.1f}-{lo+0.1:.1f}  {s.sum():7d}   {v.mean():9.2f}  "
              f"{np.median(v):8.2f}   {100*(v<1).mean():7.1f}%")

    p = SCR / f"PROF_{veri}.json"
    p.write_text(json.dumps(dict(veri=veri, T=Mo.T, dres=dres, sbar=sbar,
                                 N=len(Mo.s), anchor=[abs(v) for v in k0],
                                 lines=out, kutu=kut), indent=1))
    print(f"\n-> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)


if __name__ == "__main__":
    kos(sys.argv[1],
        float(sys.argv[2]) if len(sys.argv) > 2 else 0.40,
        float(sys.argv[3]) if len(sys.argv) > 3 else 0.95)
