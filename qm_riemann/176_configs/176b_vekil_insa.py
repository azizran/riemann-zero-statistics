# -*- coding: utf-8 -*-
"""
176b — K1 (İNŞA): FAZ-KARIŞTIRILMIŞ VEKİL GAZ
==============================================
ÖN-MÜHÜR (koşudan ÖNCE yazıldı; 176a'nın ön-kaydıyla tutarlı).

Gerçek ζ gazının NOMİNAL merdiveni `Hkeskin`'inkidir:
    S(t) = − Σ_q a_q sin(ω_q t) ,  a_q = 1/(π k √q) ,  q = p^k , τ_q ≤ 1.00
Vekil, bu GENLİK dizisini (zarfı) BİT-BİT korur ve yalnız her çizgiye
bağımsız bir φ_q ~ U(0,2π) takar:
    S_v(t) = − Σ_q a_q sin(ω_q t + φ_q)
Sonra merdiven, 164'ün SADAKATLİ zinciriyle (ızgara braketi + SIRALI
İLK-KÖK + korumalı Newton, h = 0.015, nz = 300000, c = −½) SIFIRDAN
çözülür. Çözücü `164_insa.coz_sadakatli`'den AYNEN çağrılır; kopyalanan
tek satır yoktur. Değiştirilen tek şey, çözücünün çağırdığı ALAN
DEĞERLENDİRİCİSİDİR (`176_vekil_cekirdek`, faz taşır).

İNŞA KAPILARI (176a/F0 — ön-kayıtta dondurulmuş):
  G1  ZARF ÖZDEŞLİĞİ  maks|A_q(vekil) − A_q(Hkeskin)| = 0.0 (TAM sıfır)
      ve türev nicelikleri (Σ|A|, ΣA², rms S, rms S', Σ A_qω_q) özdeş.
  G2  φ ≡ 0 SAĞLAMASI: faz taşıyan değerlendirici, φ = 0'da
      `164_insa.S_ve_Sp` ile maks fark 0.0 vermeli (aksi hâlde
      değerlendiricinin kendisi bozuktur).
  G3  ilk-kök hücre 300000 / 300000 (benzersiz)
  G4  maks|F| ≤ 1e−8, aşan tekne 0
  G5  sıralılık TAM
Biri tutmazsa gaz ÖLÇÜLMEZ.

Kullanım: 176b_vekil_insa.py <ad> <tohum> [h] [nz]
   ör.    176b_vekil_insa.py VF1 1
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
S176 = SCRR / "176"
S167 = SCRR / "167"
S155 = SCRR / "155"
TWO_PI = 2 * np.pi


def main(ad, tohum, h=I164.HIZGARA, nz=I164.NZERO, nwork=VC.NWORK):
    tbas = time.time()
    for d in (S176, S167, S155):
        d.mkdir(parents=True, exist_ok=True)
    print("=" * 74, flush=True)
    print(f"=== 176b VEKİL İNŞA {ad}: tohum={tohum} h={h} nz={nz} "
          f"nwork={nwork} ===", flush=True)
    print("=" * 74, flush=True)

    # ---------------------------------------------------------- merdiven
    d = np.load(QM / "128_odl_zeros6_2e6_zeros.npz")
    Z = np.sort(np.asarray(d["zeros"], dtype=float))
    zr = Z[len(Z) - 300000:]
    t0, t1 = float(zr[0]), float(zr[-1])
    L_hedef = float(np.log(0.5 * (t0 + t1) / TWO_PI))

    # Hkeskin'in merdiveni (164_insa.main("Hkeskin") ile AYNI çağrı)
    om_h, a_h, a_ham_h, w_h = I164.merdiven(L_hedef, None, None)
    # vekilin merdiveni (aynı çağrı — zarf birebir)
    om, A, a_ham, wgt = I164.merdiven(L_hedef, None, None)

    rng = np.random.default_rng(int(tohum))
    phi = rng.uniform(0.0, TWO_PI, len(om))

    # --- G1: ZARF ÖZDEŞLİĞİ ------------------------------------------
    zarf_fark = float(np.max(np.abs(A - a_h)))
    om_fark = float(np.max(np.abs(om - om_h)))
    sha_A = hashlib.sha256(np.ascontiguousarray(A).tobytes()).hexdigest()
    sha_Ah = hashlib.sha256(np.ascontiguousarray(a_h).tobytes()).hexdigest()
    rmsS = float(np.sqrt(0.5 * np.sum(A ** 2)))
    rmsSp = float(np.sqrt(0.5 * np.sum((A * om) ** 2)))
    sum_aw = float(np.sum(A * om))
    sum_abs = float(np.abs(A).sum())
    print(f"  merdiven: {len(om)} çizgi (τ≤1.00, L_hedef={L_hedef:.6f})",
          flush=True)
    print(f"  G1 ZARF: maks|A(vekil)−A(Hkeskin)| = {zarf_fark:.1e}   "
          f"maks|ω−ω| = {om_fark:.1e}", flush=True)
    print(f"     sha256(A)  vekil  = {sha_A}", flush=True)
    print(f"     sha256(A) Hkeskin = {sha_Ah}", flush=True)
    print(f"     Σ|A| = {sum_abs:.9f}  ΣA² = {np.sum(A**2):.9f}  "
          f"rms S = {rmsS:.6f}  rms S' = {rmsSp:.6f}  ΣA_qω_q = {sum_aw:.4f}",
          flush=True)
    print(f"     φ: n={len(phi)}  ort={phi.mean():.6f}  "
          f"⟨cos φ⟩={np.cos(phi).mean():+.6f}  ⟨sin φ⟩={np.sin(phi).mean():+.6f}",
          flush=True)
    G1_ok = (zarf_fark == 0.0) and (om_fark == 0.0) and (sha_A == sha_Ah)

    # --- G2: φ ≡ 0 SAĞLAMASI (değerlendiricinin kendisi) ---------------
    zt = np.linspace(t0, t0 + 50.0, 401)
    s_ref, sp_ref = I164.S_ve_Sp(zt, om, A, deriv=True)
    s_v, sp_v = VC.S_ve_Sp_fazli(zt, om, A, np.zeros_like(A), deriv=True)
    g2_S = float(np.max(np.abs(s_ref - s_v)))
    g2_Sp = float(np.max(np.abs(sp_ref - sp_v)))
    print(f"  G2 φ≡0 SAĞLAMASI: maks|ΔS| = {g2_S:.1e}   "
          f"maks|ΔS'| = {g2_Sp:.1e}   (401 nokta)", flush=True)
    G2_ok = (g2_S == 0.0) and (g2_Sp == 0.0)

    if not (G1_ok and G2_ok):
        raise SystemExit("KAPI G1/G2 TUTMADI — vekil kurulmaz.")

    # ------------------------------------------------------------ çözüm
    n0 = int(np.ceil(I164.rvm_N(t0)))
    ns = np.arange(n0, n0 + nz, dtype=float) - 0.5      # c = −½

    # ÇÖZÜCÜ AYNEN 164'ünkidir; yalnız alan değerlendiricisi faz taşır.
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
    g = dz
    mid = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(mid / TWO_PI)
    L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    print(f"  G3 ilk-kök hücre: {tani['hucre_benzersiz']} / {nz}", flush=True)
    print(f"  G4 maks|F| = {tani['maxF']:.3e}   (>1e−8 tekne: "
          f"{tani['nF_asan']})", flush=True)
    print(f"  G5 sıralılık: {'TAM' if sirali else 'BOZUK'}  "
          f"min Δz = {dz.min():.6f}", flush=True)
    print(f"  σ_ds = {np.std(ds):.5f}  (σ_ds² = {np.var(ds):.5f})   "
          f"L = {L:.9f}", flush=True)

    np.save(S176 / f"z_{ad}.npy", z)
    np.save(S155 / f"z_{ad}.npy", z)       # 155_kos.veri_yukle buradan okur
    np.save(S167 / f"z_{ad}.npy", z)

    tani.update(ad=ad, tip="vekil_faz", tohum=int(tohum), lam=1.0,
                tau_ust=1.00, tau_c=None, delta=None, c=-0.5,
                L=L, L_hedef=L_hedef, sirali=sirali,
                min_dz=float(dz.min()), sigma_ds=float(np.std(ds)),
                sigma_ds2=float(np.var(ds)), n=int(len(z)),
                nline=int(len(om)),
                zarf_fark=zarf_fark, om_fark=om_fark,
                sha_A=sha_A, sha_A_Hkeskin=sha_Ah,
                g2_dS=g2_S, g2_dSp=g2_Sp,
                sum_abs_A=sum_abs, sum_A2=float(np.sum(A ** 2)),
                rms_S=rmsS, rms_Sp=rmsSp, sum_a_om=sum_aw,
                kapilar=dict(G1=bool(G1_ok), G2=bool(G2_ok),
                             G3=bool(tani["hucre_benzersiz"] == nz),
                             G4=bool(tani["nF_asan"] == 0
                                     and tani["maxF"] <= 1e-8),
                             G5=bool(sirali)),
                sure_s=time.time() - tbas)
    (S176 / f"insa_{ad}.json").write_text(json.dumps(tani, indent=1))
    print(f"-> z_{ad}.npy  ({(time.time()-tbas)/60:.1f} dk)", flush=True)
    print(f"   KAPILAR: {tani['kapilar']}", flush=True)
    return tani


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise SystemExit("kullanım: 176b_vekil_insa.py <ad> <tohum> [h] [nz]")
    main(sys.argv[1], int(sys.argv[2]),
         float(sys.argv[3]) if len(sys.argv) > 3 else I164.HIZGARA,
         int(sys.argv[4]) if len(sys.argv) > 4 else I164.NZERO)
