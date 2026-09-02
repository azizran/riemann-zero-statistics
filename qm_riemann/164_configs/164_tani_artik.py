"""
164 — TANI: 152'nin merdiveni NEDEN eksik kuruyor? (ARTIK AYRIŞIMI)
===================================================================
Bir gazın ds'i, çözümün artığıyla ÖZDEŞ olarak şöyle ayrışır:

    N̄(z_{n+1}) + S(z_{n+1}) = (n+1) + F_{n+1}
    N̄(z_n)     + S(z_n)     =  n    + F_n
    ⇒ [N̄(z_{n+1}) − N̄(z_n)] = 1 − ΔS + (F_{n+1} − F_n)
    ⇒ ds_n ≡ g_n·log(m_n/2π)/2π − 1 = −ΔS_n + ΔF_n + O(g³N̄‴)

    (N̄‴ = −1/(2πt²) ≈ 2.6e−13 ⇒ üçüncü terim ≤ 1e−14, ihmal edilir)

Yani ölçülen çizgi genliği İKİ parçadan gelir:

    c_q(ds) = c_q(−ΔS)  +  c_q(ΔF)
              ⟵ merdiven      ⟵ ÇÖZÜM ARTIĞI

Sadakatli inşada |F| ≤ 2e−9 olduğundan ikinci terim ÖLÜDÜR. 152'nin
sönümlü/kelepçeli Newton'unda ise |F|'nin %99'luk dilimi 6.7e−03,
maksimumu 2.16'dır (152'nin kendi logu) — bu tanı c_q(ΔF)'yi ÖLÇER ve
"gerçekleşen/nominal = 0.54"ün tam olarak buradan geldiğini gösterir.

Kullanım:  164_tani_artik.py <veri> [<veri> ...]
   veri ∈ {keskin, A4, Skeskin, SA4, NKkeskin, eski_keskin, eski_A4}
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "155_configs"))
sys.path.insert(0, str(QM / "154_configs"))
sys.path.insert(0, str(QM / "164_configs"))
KOS155 = importlib.import_module("155_kos")
C154 = importlib.import_module("154_cekirdek")
INSA = importlib.import_module("164_insa")

SCR164 = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
              "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/164")
TWO_PI = 2 * np.pi
PEN = {"A4": (0.68, 0.125), "SA4": (0.68, 0.125), "eski_A4": (0.68, 0.125),
       "SA1": (0.75, 0.10)}
Q163 = (2, 3, 5, 7, 11, 101, 1009)


def olc(veri):
    tb = time.time()
    z = np.sort(KOS155.veri_yukle(veri))
    d = np.load(QM / "128_odl_zeros6_2e6_zeros.npz")
    Z = np.sort(np.asarray(d["zeros"], float))
    zr = Z[len(Z) - 300000:]
    t0 = float(zr[0])
    Lh = float(np.log(0.5 * (zr[0] + zr[-1]) / TWO_PI))
    p = PEN.get(veri)
    om, a, a_ham, _ = INSA.merdiven(Lh, p[0] if p else None, p[1] if p else None)
    n0 = int(np.ceil(INSA.rvm_N(t0)))
    ns = np.arange(n0, n0 + len(z), dtype=float)

    S, _ = INSA.S_ve_Sp(z, om, a, deriv=False)
    F = INSA.rvm_N(z) + S - ns
    g = np.diff(z)
    mid = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(mid / TWO_PI)
    L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    dS = -(S[1:] - S[:-1])              # −ΔS
    dF = F[1:] - F[:-1]                 # ΔF
    kap = float(np.max(np.abs(ds - (dS + dF))))
    aF = np.abs(F)
    print(f"=== ARTIK TANISI {veri}: n={len(z)}  L={L:.4f}", flush=True)
    print(f"  |F|: medyan={np.median(aF):.3e}  %99={np.percentile(aF,99):.3e} "
          f" maks={aF.max():.3e}   (|F|>1e-8 tekne: {int((aF>1e-8).sum())}, "
          f"%{100*np.mean(aF>1e-8):.2f})", flush=True)
    print(f"  özdeşlik kapanışı  maks|ds − (−ΔS + ΔF)| = {kap:.2e}", flush=True)
    print(f"  Var(ds)={np.var(ds):.5f}  Var(−ΔS)={np.var(dS):.5f}  "
          f"Var(ΔF)={np.var(dF):.5f}  kor(−ΔS,ΔF)="
          f"{np.corrcoef(dS,dF)[0,1]:+.4f}", flush=True)

    dsm = ds - ds.mean()
    dSm = dS - dS.mean()
    dFm = dF - dF.mean()
    qm = C154.pk_m(int(np.exp(1.0 * L)))
    print("\n   q   τ      b_nom     c(ds)     c(−ΔS)    c(ΔF)   "
          "  ds/b   −ΔS/b   ΔF/b", flush=True)
    sat = []
    for q in Q163:
        if q not in qm:
            continue
        w = np.log(q)
        tau = w / L
        aq = 1.0 / (np.pi * qm[q] * np.sqrt(q))
        if p is not None:
            import math
            aq *= 0.5 * math.erfc((tau - p[0]) / p[1])
        b = 2 * aq * np.sin(np.pi * tau)
        e = np.exp(-1j * w * mid)
        c1 = 2 * np.mean(dsm * e)
        c2 = 2 * np.mean(dSm * e)
        c3 = 2 * np.mean(dFm * e)
        print(f"{q:5d} {tau:.4f} {b:9.5f} {abs(c1):9.5f} {abs(c2):9.5f} "
              f"{abs(c3):9.5f}  {abs(c1)/b:6.3f}  {abs(c2)/b:6.3f}  "
              f"{abs(c3)/b:6.3f}   (Re c(ΔF)/b = {(c3/b).real:+.3f})",
              flush=True)
        sat.append(dict(q=int(q), tau=float(tau), b=float(b),
                        c_ds=[c1.real, c1.imag], c_dS=[c2.real, c2.imag],
                        c_dF=[c3.real, c3.imag]))
    out = dict(veri=veri, L=L, medF=float(np.median(aF)),
               p99F=float(np.percentile(aF, 99)), maxF=float(aF.max()),
               nF=int((aF > 1e-8).sum()), kapanis=kap,
               var_ds=float(np.var(ds)), var_dS=float(np.var(dS)),
               var_dF=float(np.var(dF)),
               kor=float(np.corrcoef(dS, dF)[0, 1]), cizgi=sat,
               sure_s=time.time() - tb)
    (SCR164 / f"artik_{veri}.json").write_text(json.dumps(out, indent=1))
    print(f"-> artik_{veri}.json ({time.time()-tb:.0f}s)\n", flush=True)


if __name__ == "__main__":
    for v in sys.argv[1:]:
        olc(v)
