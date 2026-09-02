"""
164 — TANI: İLK-KÖK KURALININ FAZ YANLILIĞI ve SEVİYE KONVANSİYONU
==================================================================
Sadakat denetimi (§3) beklenmedik bir şey buldu: sadakatli gazın çizgileri
GENLİKÇE 1.00'i tutuyor ama bir FAZ taşıyor —

    arg c_q ≈ −0.18·ω_q   (Skeskin),   ≈ −0.12·ω_q   (SA4)

gerçek gazda ve 152'nin gazlarında bu faz YOK (|arg| ≤ 0.02 rad).
`ds ⊃ A cos(ω_q(m − Δ))` demek olduğundan bu, gerçekleşen yapının
merdivene göre `Δ ≈ 0.18` (≈ 0.34 ḡ) kaymış olması demektir — İLK-KÖK
kuralının "yükselen kenar" yanlılığının beklenen imzası.

SINAV: denklemi `N̄(z) + S(z) = n + c` olarak koş ve `c`'yi tara.
`c` bir KONVANSİYONdur: gerçek sayma fonksiyonu sıfırda +1 atlar,
asal-toplam onun DÜZGÜN sürümüdür, yani `N̄(γ_n) + S_düz(γ_n) ≈ n − ½`.
152 (ve §2) `c = 0` kullandı. Seviye kayması KATI BİR ÖTELEME DEĞİLDİR
(`δz ≈ c/F'` ve `F'` merdivene kilitli dalgalanır), bu yüzden `c`
fazı gerçekten değiştirebilir.

Ucuz: ızgara `c`'den bağımsızdır — bir kez hesaplanır, her `c` için
yalnız `searchsorted` + korumalı Newton tekrarlanır.

Kullanım:  164_tani_seviye.py [nz]     (varsayılan 40000 tekne)
"""
import json
import sys
import time
import multiprocessing as mp
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "164_configs"))
sys.path.insert(0, str(QM / "154_configs"))
import importlib
I = importlib.import_module("164_insa")
C154 = importlib.import_module("154_cekirdek")

SCR164 = I.SCR164
TWO_PI = 2 * np.pi
CLER = (-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5)
Q = (2, 3, 5, 7, 11, 101)


def main(nz=40000):
    tb = time.time()
    d = np.load(QM / "128_odl_zeros6_2e6_zeros.npz")
    Z = np.sort(np.asarray(d["zeros"], float))
    zr = Z[len(Z) - 300000:]
    t0 = float(zr[0])
    Lh = float(np.log(0.5 * (zr[0] + zr[-1]) / TWO_PI))
    om, a, _, _ = I.merdiven(Lh)
    n0 = int(np.ceil(I.rvm_N(t0)))
    ns = np.arange(n0, n0 + nz, dtype=float)
    print(f"=== 164 SEVİYE TARAMASI  nz={nz}  çizgi={len(om)} ===", flush=True)
    ctx = mp.get_context("spawn")
    pool = ctx.Pool(I.NWORK, initializer=I._init, initargs=(om, a))
    out = []
    try:
        for c in CLER:
            z, F, tani = I.coz_sadakatli(om, a, ns + c, t0, pool,
                                         h=I.HIZGARA, log=lambda *x: None)
            g = np.diff(z)
            mid = 0.5 * (z[:-1] + z[1:])
            Lw = np.log(mid / TWO_PI)
            L = float(Lw.mean())
            ds = g * Lw / TWO_PI - 1
            dsm = ds - ds.mean()
            qm = C154.pk_m(int(np.exp(1.0 * L)))
            sat = {}
            for q in Q:
                w = np.log(q)
                tau = w / L
                aq = 1.0 / (np.pi * qm[q] * np.sqrt(q))
                b = 2 * aq * np.sin(np.pi * tau)
                cc = 2 * np.mean(dsm * np.exp(-1j * w * mid))
                sat[q] = (abs(cc) / b, float(np.angle(cc)),
                          float(np.angle(cc)) / w)
            egim = np.mean([sat[q][2] for q in Q])
            print(f"  c={c:+.2f}  maks|F|={tani['maxF']:.2e}  "
                  f"σ_ds²={np.var(ds):.4f}  sıra bozan="
                  f"{int((g<=0).sum())}  |c_q|/b: "
                  + " ".join(f"{sat[q][0]:.3f}" for q in Q)
                  + f"   arg/ω ort = {egim:+.4f}", flush=True)
            print("        arg c_q: "
                  + " ".join(f"{sat[q][1]:+.4f}" for q in Q), flush=True)
            out.append(dict(c=c, maxF=tani["maxF"], sigma_ds2=float(np.var(ds)),
                            sat={str(q): sat[q] for q in Q}, egim=float(egim)))
    finally:
        pool.close()
        pool.join()
    (SCR164 / "seviye_taramasi.json").write_text(json.dumps(out, indent=1))
    print(f"-> seviye_taramasi.json  ({(time.time()-tb)/60:.1f} dk)",
          flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 40000)
