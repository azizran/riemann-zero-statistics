"""
59 — YÜKSEK-HASSASİYETLİ CUE KAMPANYASI (K1 düzeltmesi) (17 Ağustos, gece)
===========================================================================

Denetim K1: 42b'nin kuadratik r_CUE(N) fiti kendi MC hatalarına karşı
χ²/dof≈29 ile reddediliyor; N_eff sabitliği o fitin artefaktı olabilir.

Çözüm: fit AİLESİNİ tamamen terk et —
  - her N için 1.5M aralık (MC se ~4e-4; eskinin ~7.5 katı istatistik)
  - N = 5..20 tamsayı noktaları
  - hata: 10-parça batch dağılımından (formül değil, ampirik)
  - interpolasyon: PCHIP (monoton, form varsayımı yok) → 60'ta

Çıktı: 59_cue_hassas.npz (N, rP, seP, rS, seS)
"""

import numpy as np
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
import time

rng = np.random.default_rng(59)
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

def haar_unitary_batch(M, N):
    G = (rng.standard_normal((M, N, N)) + 1j * rng.standard_normal((M, N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(G)
    diag = np.einsum("mii->mi", R)
    return Q * (diag / np.abs(diag))[:, None, :]

def cue_batch(N, n_gaps, grid=48, chunk=1200):
    M = int(np.ceil(n_gaps / N))
    gs, ms = [], []
    for s in range(0, M, chunk):
        m = min(chunk, M - s)
        U = haar_unitary_batch(m, N)
        ph = np.sort(np.angle(np.linalg.eigvals(U)), axis=1)
        gp = np.diff(np.concatenate([ph, ph[:, :1] + TWO_PI], axis=1), axis=1)
        u = np.arange(1, grid + 1) / (grid + 1)
        th = ph[:, :, None] + gp[:, :, None] * u[None, None, :]
        df = th[:, :, :, None] - ph[:, None, None, :]
        amp = np.prod(2 * np.abs(np.sin(df / 2)), axis=-1)
        gs.append(gp.ravel()); ms.append(amp.max(axis=-1).ravel())
    return np.concatenate(gs), np.concatenate(ms)

NS = np.arange(5, 21)
N_GAPS = 1_500_000
N_BATCH = 10

out_N, out_rP, out_seP, out_rS, out_seS = [], [], [], [], []
for N in NS:
    t1 = time.time()
    rPs, rSs = [], []
    for b in range(N_BATCH):
        g, a = cue_batch(N, N_GAPS // N_BATCH)
        rPs.append(pearsonr(g, a)[0])
        rSs.append(spearmanr(g, a)[0])
    rPs, rSs = np.array(rPs), np.array(rSs)
    rP, seP = rPs.mean(), rPs.std(ddof=1) / np.sqrt(N_BATCH)
    rS, seS = rSs.mean(), rSs.std(ddof=1) / np.sqrt(N_BATCH)
    out_N.append(N); out_rP.append(rP); out_seP.append(seP)
    out_rS.append(rS); out_seS.append(seS)
    print(f"N={N:>2}: P {rP:.5f}±{seP:.5f}  S {rS:.5f}±{seS:.5f}  "
          f"({(time.time()-t1)/60:.1f} dk)", flush=True)
    # ara kayıt (ders alındı: pencere-başına kaydet)
    np.savez(HERE / "59_cue_hassas.npz", N=np.array(out_N),
             rP=np.array(out_rP), seP=np.array(out_seP),
             rS=np.array(out_rS), seS=np.array(out_seS))

print("Tamam: 59_cue_hassas.npz", flush=True)
