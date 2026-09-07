# -*- coding: utf-8 -*-
"""
180b — K1 (İNŞA): GERÇEK-ZARFLI KARIŞIK VEKİL (VS)
===================================================
ÖN-MÜHÜR: 180a'nın ONKAYIT_K0.json'u koşudan ÖNCE diske yazıldı; bu
betik onu okur, sha'sını doğrular ve YALNIZ orada dondurulmuş zarfı
kullanır. Hiçbir yüzleşme niceliği burada hesaplanmaz.

Makine 176'nınkidir, AYNEN:
  * alan değerlendiricisi `176_vekil_cekirdek` (faz taşır, φ≡0'da
    `164_insa.S_ve_Sp` ile bit-bit aynı),
  * çözücü `164_insa.coz_sadakatli` (ızgara braketi + SIRALI İLK-KÖK +
    korumalı Newton, h = 0.015, nz = 300000, c = −½) — kopyalanan satır
    yok, yalnız `S_par`/`SSp_par` isimleri değiştiriliyor.
Tek fark: genlik dizisi Hkeskin'in a_q'su değil, 180a'nın dondurduğu
    A_q(VS) = a_q · r(τ_q),  r(τ) = R_bant(son)/R_bant(Hkeskin)
dizisidir; her iki tohumda BİT-BİT aynıdır (kapı V1/V2, sha ile).

KAPILAR (180a'da dondurulmuş):
  V1 sha256(A_VS) == ön-kayıttaki sha
  V2 maks|A_VS − a_Hk·r(τ)| == 0.0  (diskteki diziyle tarif birebir)
  V3 φ≡0 sağlaması 0.0 / 0.0
  V4 ilk-kök hücre 300000 / 300000
  V5 maks|F| ≤ 1e−8, aşan tekne 0
  V6 sıralılık TAM
Biri tutmazsa gaz ÖLÇÜLMEZ.

Kullanım: 180b_vs_insa.py <ad> <tohum>       ör. 180b_vs_insa.py VS1 1
"""
import hashlib
import importlib
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "176_configs"))
sys.path.insert(0, str(QM / "164_configs"))
I164 = importlib.import_module("164_insa")
VC = importlib.import_module("176_vekil_cekirdek")

SCRR = I164.SCRR
S180 = SCRR / "180"
S167 = SCRR / "167"
S155 = SCRR / "155"
TWO_PI = 2 * np.pi


def main(ad, tohum, h=I164.HIZGARA, nz=I164.NZERO, nwork=VC.NWORK):
    tbas = time.time()
    for d in (S180, S167, S155):
        d.mkdir(parents=True, exist_ok=True)
    OK = json.load(open(S180 / "ONKAYIT_K0.json"))
    print("=" * 74, flush=True)
    print(f"=== 180b VS İNŞA {ad}: tohum={tohum} h={h} nz={nz} "
          f"nwork={nwork} ===", flush=True)
    print(f"    ön-kayıt: {OK['zaman']}  sha={OK['sha256'][:16]}…", flush=True)
    print("=" * 74, flush=True)

    # ---------------------------------------------------------- merdiven
    d = np.load(QM / "128_odl_zeros6_2e6_zeros.npz")
    Z = np.sort(np.asarray(d["zeros"], dtype=float))
    zr = Z[len(Z) - 300000:]
    t0, t1 = float(zr[0]), float(zr[-1])
    L_hedef = float(np.log(0.5 * (t0 + t1) / TWO_PI))
    om, a_h, _, _ = I164.merdiven(L_hedef, None, None)

    A = np.load(S180 / "A_sonzarf.npy")           # 180a'nın dondurduğu dizi
    tg, rr = np.load(S180 / "r_tau_profil.npy")
    A_tarif = a_h * np.interp(om / L_hedef, tg, rr)

    sha_A = hashlib.sha256(np.ascontiguousarray(A).tobytes()).hexdigest()
    v2 = float(np.max(np.abs(A - A_tarif)))
    V1 = (sha_A == OK["zarf"]["sha_A_VS"])
    V2 = (v2 == 0.0) and (len(A) == len(om))
    print(f"  merdiven: {len(om)} çizgi (τ≤1.00, L_hedef={L_hedef:.9f})",
          flush=True)
    print(f"  V1 sha256(A_VS) = {sha_A}  {'✓' if V1 else '✗'}", flush=True)
    print(f"     ön-kayıt      = {OK['zarf']['sha_A_VS']}", flush=True)
    print(f"  V2 maks|A_VS − a_Hk·r(τ)| = {v2:.1e}  {'✓' if V2 else '✗'}",
          flush=True)
    print(f"     Σ|A| = {np.abs(A).sum():.9f}  ΣA² = {np.sum(A**2):.9f}  "
          f"(Hk: {np.abs(a_h).sum():.9f} / {np.sum(a_h**2):.9f})", flush=True)

    rng = np.random.default_rng(int(tohum))
    phi = rng.uniform(0.0, TWO_PI, len(om))
    print(f"     φ: n={len(phi)} ort={phi.mean():.6f} "
          f"⟨cosφ⟩={np.cos(phi).mean():+.6f} "
          f"⟨sinφ⟩={np.sin(phi).mean():+.6f}", flush=True)

    zt = np.linspace(t0, t0 + 50.0, 401)
    s_ref, sp_ref = I164.S_ve_Sp(zt, om, A, deriv=True)
    s_v, sp_v = VC.S_ve_Sp_fazli(zt, om, A, np.zeros_like(A), deriv=True)
    g3_S = float(np.max(np.abs(s_ref - s_v)))
    g3_Sp = float(np.max(np.abs(sp_ref - sp_v)))
    V3 = (g3_S == 0.0) and (g3_Sp == 0.0)
    print(f"  V3 φ≡0 SAĞLAMASI: maks|ΔS| = {g3_S:.1e}  maks|ΔS'| = "
          f"{g3_Sp:.1e}  {'✓' if V3 else '✗'}", flush=True)

    if not (V1 and V2 and V3):
        raise SystemExit("KAPI V1/V2/V3 TUTMADI — vekil kurulmaz.")

    # ------------------------------------------------------------ çözüm
    n0 = int(np.ceil(I164.rvm_N(t0)))
    ns = np.arange(n0, n0 + nz, dtype=float) - 0.5
    I164.S_par = VC.S_par_v
    I164.SSp_par = VC.SSp_par_v
    ctx = mp.get_context("spawn")
    pool = ctx.Pool(nwork, initializer=VC._init_v, initargs=(om, A, phi))
    try:
        z, F, tani = I164.coz_sadakatli(om, A, ns, t0, pool, h=h)
    finally:
        pool.close()
        pool.join()

    dz = np.diff(z)
    sirali = bool(np.all(dz > 0))
    mid = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(mid / TWO_PI)
    L = float(Lw.mean())
    ds = dz * Lw / TWO_PI - 1
    print(f"  V4 ilk-kök hücre: {tani['hucre_benzersiz']} / {nz}", flush=True)
    print(f"  V5 maks|F| = {tani['maxF']:.3e}  (>1e−8 tekne: "
          f"{tani['nF_asan']})", flush=True)
    print(f"  V6 sıralılık: {'TAM' if sirali else 'BOZUK'}  "
          f"min Δz = {dz.min():.6f}", flush=True)
    print(f"  σ_ds = {np.std(ds):.5f}   L = {L:.9f}", flush=True)

    np.save(S180 / f"z_{ad}.npy", z)
    np.save(S155 / f"z_{ad}.npy", z)
    np.save(S167 / f"z_{ad}.npy", z)

    tani.update(ad=ad, tip="vekil_faz_sonzarf", tohum=int(tohum), lam=1.0,
                tau_ust=1.00, tau_c=None, delta=None, c=-0.5, L=L,
                L_hedef=L_hedef, sirali=sirali, min_dz=float(dz.min()),
                sigma_ds=float(np.std(ds)), sigma_ds2=float(np.var(ds)),
                n=int(len(z)), nline=int(len(om)), sha_A=sha_A,
                sha_A_onkayit=OK["zarf"]["sha_A_VS"], V2_fark=v2,
                g3_dS=g3_S, g3_dSp=g3_Sp,
                sum_abs_A=float(np.abs(A).sum()),
                sum_A2=float(np.sum(A ** 2)),
                kapilar=dict(V1=bool(V1), V2=bool(V2), V3=bool(V3),
                             V4=bool(tani["hucre_benzersiz"] == nz),
                             V5=bool(tani["nF_asan"] == 0
                                     and tani["maxF"] <= 1e-8),
                             V6=bool(sirali)),
                sure_s=time.time() - tbas)
    (S180 / f"insa_{ad}.json").write_text(json.dumps(tani, indent=1,
                                                     default=float))
    print(f"-> z_{ad}.npy  ({(time.time()-tbas)/60:.1f} dk)", flush=True)
    print(f"   KAPILAR: {tani['kapilar']}", flush=True)
    return tani


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise SystemExit("kullanım: 180b_vs_insa.py <ad> <tohum>")
    main(sys.argv[1], int(sys.argv[2]))
