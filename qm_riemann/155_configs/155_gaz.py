"""
155 — EK SENTETİK TANIKLAR (153'ün gazları, τ₀ için)
====================================================
153_gaz.py'nin gaz üreticilerini KOPYALAMAZ, import eder. Üretilen z
dizileri scratchpad/155/z_<ad>.npy'ye yazılır; 155_kos.py oradan okur.

NİÇİN: A4 (erfc-0.68) ve keskin gazlarında τ₀ ≈ 0.489, gerçekte ≈ 0.509.
Fark hangi istatistikten geliyor? 153'ün gaz merdiveni tam bu soruyu
ayırmak için kurulmuştu:

  N5 / N5z : kısa-menzilli YAPI değişikliği (itme süpürmesi, iki doz).
             σ_ds² büyür VE kısa-menzil dağılımı gerçeğe yaklaşır.
  J14/ J26 : YAPISIZ titreşim, N5/N5x ile AYNI σ_ds²'yi üretecek dozda.
             ⇒ τ₀ N-serisinde kayıp J-serisinde kaymıyorsa, τ₀ momentin
             değil KISA-MENZİL YAPISININ fonksiyonudur (H-N'in nicel
             sürümüne karşı doğrudan kanıt).
  P1 / P0  : bağımsız gap örnekleyicili taban (GUE / Poisson) + tek
             geçiş merdiven boyası. P1, gerçeğin ikili itmesine en yakın
             sentetik.

Kullanım: 155_gaz.py <N5|N5z|J14|J26|P1|P0|Pd>
"""
import importlib
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "153_configs"))
G153 = importlib.import_module("153_gaz")

HERE = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR153 = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
              "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/153")
OUT = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/155")
TWO_PI = 2 * np.pi


def merdiven():
    """153/152 ile birebir aynı erfc-0.68/0.125 kesimli merdiven."""
    d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
    Z = np.sort(np.asarray(d["zeros"], dtype=float))
    zr = Z[len(Z) - 300000:]
    t0, t1 = float(zr[0]), float(zr[-1])
    Lh = float(np.log(0.5 * (t0 + t1) / TWO_PI))
    from sympy import primerange
    lim = int(np.exp(1.00 * Lh))
    lad = []
    for p in primerange(2, lim + 1):
        q, m = p, 1
        while q <= lim:
            lad.append((np.log(q), 1.0 / (np.pi * m * np.sqrt(q))))
            q *= p
            m += 1
    om = np.array([w for w, _ in lad])
    a = np.array([x for _, x in lad])
    w = np.array([0.5 * math.erfc((x - G153.TAU_C) / G153.DELTA)
                  for x in om / Lh])
    return om, a * w, t0, t1, Lh


def main(ad):
    tb = time.time()
    tip, par = G153.KONFIG[ad]
    print(f"=== 155 gaz {ad}  ({tip} {par}) ===", flush=True)
    zt = SCR153 / "z_taban.npy"
    if tip == "itme":
        z = np.load(zt)
        print(f"  153 taban (erfc-0.68) okundu, {len(z)} nokta", flush=True)
        z = G153.itme_supur(z, par["eps"], par["nsup"])
    elif tip == "titresim":
        z = np.load(zt)
        sj = par["sigma_j"]
        gbar = 1.0 / G153.rvm_d(z.mean())
        rj = np.random.default_rng(7)          # 153 ile AYNI tohum
        z = z + rj.normal(0.0, sj * gbar, size=len(z))
        capraz = int(np.sum(np.diff(z) < 0))
        z = np.sort(z)
        print(f"  titreşim σ_j={sj} (σ={sj*gbar:.4f}); sıra bozan {capraz}",
              flush=True)
    elif tip == "boya":
        import multiprocessing as mp
        om, a, t0, t1, Lh = merdiven()
        n0 = int(np.ceil(G153.rvm_N(t0)))
        rng = np.random.default_rng(11)        # 153 ile AYNI tohum
        s = G153.gap_ornekle(par["taban"], 299999, rng)
        s = s / s.mean()
        u = np.concatenate(([0.0], np.cumsum(s))) + n0
        print(f"  taban gap: {par['taban']} var={s.var():.4f} "
              f"P(s<0.3)={np.mean(s<0.3):.5f}", flush=True)
        x = np.full_like(u, 0.5 * (t0 + t1))
        for _ in range(60):
            F = G153.rvm_N(x) - u
            x = np.clip(x - F / G153.rvm_d(x), 100.0, None)
            if np.max(np.abs(F)) < 1e-9:
                break
        ctx = mp.get_context("spawn")
        pool = ctx.Pool(G153.NWORK, initializer=G153._init,
                        initargs=(om, a))
        try:
            S = np.concatenate([r[0] for r in
                                pool.map(G153._work, np.array_split(x,
                                                                    G153.NWORK))])
        finally:
            pool.close(); pool.join()
        z = x - S / G153.rvm_d(x)
        print(f"  merdiven boyası: |Δz|ort={np.mean(np.abs(S/G153.rvm_d(x))):.5f}",
              flush=True)
    else:
        raise SystemExit(f"{ad}: bu sürücü yalnız itme/titresim/boya üretir")
    z = np.sort(z)
    g = np.diff(z)
    mid = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(mid / TWO_PI)
    ds = g * Lw / TWO_PI - 1
    print(f"  σ_ds² = {np.var(ds):.4f}   P(s<0.3) = "
          f"{np.mean((g/(TWO_PI/Lw)) < 0.3):.5f}   "
          f"(gerçek: 0.1674 / 0.02420)", flush=True)
    p = OUT / f"z_{ad}.npy"
    np.save(p, z)
    print(f"-> {p}  ({time.time()-tb:.0f}s)", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
