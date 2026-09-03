"""
167 — ABLASYON GAZLARININ İNŞASI (164'ün SADAKATLİ çözücüsüyle)
==============================================================
164_insa'nın `coz_sadakatli`'si (ızgara braketi + korumalı Newton, sıralı
İLK-KÖK) AYNEN kullanılır; hiçbir çözüm parçası kopyalanmaz. Değişen
YALNIZ merdivendir:

    S(t) = −Σ_q (λ·a_q·w_q) sin(ω_q t) ,  a_q = 1/(π m √q)

  (a) GENLİK ÖLÇEĞİ  λ ∈ {0.70, 0.85, 1.00, 1.15}      → σ'lar λ ile ölçeklenir
  (b) KESİM          keskin τ ≤ τ_ust  ya da  erfc(τ_c, Δ)

λ = 1.00 + keskin τ≤1.00 = **Hkeskin** (164'ün gazı, yeniden üretilmez)
λ = 1.00 + erfc 0.68/0.125 = **HA4**   (164'ün gazı)

Seviye konvansiyonu her gazda c = −½ (164 §3d): N̄ + S = n − ½.

Kullanım:  167_insa.py <ad> [h] [nz] [nwork]
"""
import importlib
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "164_configs"))
I164 = importlib.import_module("164_insa")

SCRR = I164.SCRR
SCR167 = SCRR / "167"
SCR155 = SCRR / "155"
TWO_PI = 2 * np.pi

# ad -> (lam, tau_ust, tau_c, delta)  ; hepsi c = −½
KONFIG = {
    # (a) GENLİK ÖLÇEĞİ ekseni (keskin τ≤1.00; λ=1.00 → Hkeskin)
    "L070": dict(lam=0.70),
    "L085": dict(lam=0.85),
    "L115": dict(lam=1.15),
    # (c) KESİM ekseni (λ=1.00)
    "K090": dict(tau_ust=0.90),
    "K070": dict(tau_ust=0.70),
    "E060": dict(tau_c=0.60, delta=0.125),
}


def main(ad, h=I164.HIZGARA, nz=I164.NZERO, nwork=I164.NWORK):
    par = KONFIG[ad]
    lam = float(par.get("lam", 1.0))
    tau_ust = float(par.get("tau_ust", 1.00))
    tau_c = par.get("tau_c")
    delta = par.get("delta")
    tbas = time.time()
    SCR167.mkdir(parents=True, exist_ok=True)
    SCR155.mkdir(parents=True, exist_ok=True)
    print(f"=== 167 İNŞA {ad}: λ={lam} τ_ust={tau_ust} erfc={tau_c}/{delta} "
          f"h={h} nz={nz} ===", flush=True)

    d = np.load(QM / "128_odl_zeros6_2e6_zeros.npz")
    Z = np.sort(np.asarray(d["zeros"], dtype=float))
    zr = Z[len(Z) - 300000:]
    t0, t1 = float(zr[0]), float(zr[-1])
    L_hedef = float(np.log(0.5 * (t0 + t1) / TWO_PI))

    om, a, a_ham, wgt = I164.merdiven(L_hedef, tau_c, delta, tau_ust=tau_ust)
    a = lam * a                                   # GENLİK ÖLÇEĞİ
    print(f"merdiven: {len(om)} çizgi (τ≤{tau_ust}, L={L_hedef:.4f})  "
          f"Σa/Σa_ham={a.sum()/a_ham.sum():.4f}", flush=True)
    print(f"  N̄'={I164.rvm_d(0.5*(t0+t1)):.4f}  rms S'="
          f"{np.sqrt(0.5*np.sum((a*om)**2)):.4f}  Σa_qω_q={np.sum(a*om):.1f}"
          f"  rms S={np.sqrt(0.5*np.sum(a**2)):.4f}", flush=True)

    n0 = int(np.ceil(I164.rvm_N(t0)))
    ns = np.arange(n0, n0 + nz, dtype=float) - 0.5     # c = −½

    ctx = mp.get_context("spawn")
    pool = ctx.Pool(nwork, initializer=I164._init, initargs=(om, a))
    try:
        I164.NWORK = nwork
        z, F, tani = I164.coz_sadakatli(om, a, ns, t0, pool, h=h)
    finally:
        pool.close()
        pool.join()

    dz = np.diff(z)
    sirali = bool(np.all(dz > 0))
    g = dz
    mid = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(mid / TWO_PI)
    L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    print(f"  sıralılık: {'TAM' if sirali else 'BOZUK'}  min Δz={dz.min():.6f}",
          flush=True)
    print(f"  σ_ds = {np.std(ds):.4f}  (σ_ds² = {np.var(ds):.4f})   L={L:.4f}",
          flush=True)

    np.save(SCR167 / f"z_{ad}.npy", z)
    np.save(SCR155 / f"z_{ad}.npy", z)     # 155_kos.veri_yukle buradan okur
    tani.update(ad=ad, lam=lam, tau_ust=tau_ust, tau_c=tau_c, delta=delta,
                L=L, L_hedef=L_hedef, sirali=sirali,
                min_dz=float(dz.min()), sigma_ds=float(np.std(ds)),
                sigma_ds2=float(np.var(ds)), n=int(len(z)),
                nline=int(len(om)), sure_s=time.time() - tbas)
    (SCR167 / f"insa_{ad}.json").write_text(json.dumps(tani, indent=1))
    print(f"-> z_{ad}.npy  ({(time.time()-tbas)/60:.1f} dk)", flush=True)


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in KONFIG:
        raise SystemExit(f"kullanım: 167_insa.py <{'|'.join(KONFIG)}>")
    main(sys.argv[1],
         float(sys.argv[2]) if len(sys.argv) > 2 else I164.HIZGARA,
         int(sys.argv[3]) if len(sys.argv) > 3 else I164.NZERO,
         int(sys.argv[4]) if len(sys.argv) > 4 else I164.NWORK)
