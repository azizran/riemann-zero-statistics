# -*- coding: utf-8 -*-
"""
183b — K1 (İNŞA): DERİNLİK-DÜĞMELİ VEKİL GAZ (TEK KOŞU, BEŞ BASAMAK)
=====================================================================
ÖN-MÜHÜR: 183/ONKAYIT_183.json (183a, 19:47:00, veriden önce).

ÇEKİRDEK DEĞİŞTİRİLMEZ. `164_insa.coz_sadakatli` AYNEN çağrılır; tek
satırı düzenlenmez. Derinlik düğmesi, çekirdeğin çağırdığı
`SSp_par` adına takılan bir SAYAÇ-SARMALAYICIDIR:

  * sarmalayıcı gerçek değerlendiriciyi (176_vekil_cekirdek.SSp_par_v)
    AYNEN çağırır ve DEĞİŞTİRMEDEN döndürür — sayısal hiçbir şey
    bozulmaz, d=∞ çıktısı 176/z_VF1.npy ile bit-bit aynı olmalıdır;
  * çekirdeğin aktif-küme defterini (akt) DIŞARIDAN yeniden kurar
    (aynı TOLIN, aynı F tanımı) ve her Newton yinelemesinin GİRİŞİNDE
    z'nin tam kopyasını saklar.

Yineleme d'nin girişindeki z = tam olarak d korumalı-Newton güncellemesi
almış z. Böylece:
    d = 0 : yalnız ızgara braketi + doğrusal ara değer (GERİ-BESLEME YOK)
    d = 1, 2, 4 : sığ öz-tutarlılık
    d = ∞ : çekirdeğin kendi yakınsaması (= VF1)
hepsi TEK koşudan, aynı tohum/ızgara/braket/faz dizisiyle çıkar.

R-KAPISI: d=∞ ile 176/z_VF1.npy arasında maks|Δz| = 0.0 olmalı.
G4/G5 sığ derinliklerde UYGULANMAZ (ön-kayıt): sığ gazlar TEŞHİS gazıdır.

Kullanım: 183b_derinlik_insa.py [tohum]
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
S183, S176, S167, S155 = SCRR / "183", SCRR / "176", SCRR / "167", SCRR / "155"
TWO_PI = 2 * np.pi
DERINLIKLER = [0, 1, 2, 4]           # ön-kayıt; ∞ ayrıca
ADLAR = {0: "VD0", 1: "VD1", 2: "VD2", 4: "VD4"}


class DerinlikSarmalayici:
    """SSp_par'ı sarar: değeri DEĞİŞTİRMEZ, yalnız z'yi anlık-görüntüler."""

    def __init__(self, ns, tolin, gercek):
        self.ns = ns
        self.tolin = tolin
        self.gercek = gercek
        self.akt = None
        self.it = 0
        self.z = None
        self.anlik = {}
        self.tutarli = True
        self.nakt = []

    def __call__(self, pool, za):
        akt = np.arange(len(self.ns)) if self.akt is None else self.akt
        if len(akt) != len(za):
            self.tutarli = False          # defter tutmadı — R-KAPISI düşer
            akt = akt[:len(za)]
        if self.z is None:
            self.z = np.array(za, dtype=float, copy=True)
        else:
            self.z[akt] = za
        if self.it in DERINLIKLER:
            self.anlik[self.it] = self.z.copy()
        self.nakt.append(int(len(akt)))
        Sa, Spa = self.gercek(pool, za)
        F = I164.rvm_N(za) + Sa - self.ns[akt]
        kal = np.abs(F) > self.tolin
        self.akt = akt[kal]
        self.it += 1
        return Sa, Spa                    # DEĞİŞTİRİLMEDEN


def main(tohum=1, h=I164.HIZGARA, nz=I164.NZERO, nwork=VC.NWORK):
    tb = time.time()
    for d in (S183, S176, S167, S155):
        d.mkdir(parents=True, exist_ok=True)
    print("=" * 78, flush=True)
    print(f"183b — DERİNLİK MERDİVENİ (tohum={tohum}, h={h}, nz={nz})",
          flush=True)
    print("=" * 78, flush=True)

    d = np.load(QM / "128_odl_zeros6_2e6_zeros.npz")
    Z = np.sort(np.asarray(d["zeros"], dtype=float))
    zr = Z[len(Z) - 300000:]
    t0, t1 = float(zr[0]), float(zr[-1])
    L_hedef = float(np.log(0.5 * (t0 + t1) / TWO_PI))

    om_h, a_h, _, _ = I164.merdiven(L_hedef, None, None)
    om, A, _, _ = I164.merdiven(L_hedef, None, None)
    rng = np.random.default_rng(int(tohum))
    phi = rng.uniform(0.0, TWO_PI, len(om))

    # --- G1 / G2 (176b'nin kapıları, AYNEN) ---------------------------
    zarf_fark = float(np.max(np.abs(A - a_h)))
    om_fark = float(np.max(np.abs(om - om_h)))
    sha_A = hashlib.sha256(np.ascontiguousarray(A).tobytes()).hexdigest()
    sha_Ah = hashlib.sha256(np.ascontiguousarray(a_h).tobytes()).hexdigest()
    G1 = (zarf_fark == 0.0) and (om_fark == 0.0) and (sha_A == sha_Ah)
    zt = np.linspace(t0, t0 + 50.0, 401)
    s_ref, sp_ref = I164.S_ve_Sp(zt, om, A, deriv=True)
    s_v, sp_v = VC.S_ve_Sp_fazli(zt, om, A, np.zeros_like(A), deriv=True)
    g2S, g2Sp = float(np.max(np.abs(s_ref - s_v))), float(np.max(np.abs(sp_ref - sp_v)))
    G2 = (g2S == 0.0) and (g2Sp == 0.0)
    print(f"  G1 ZARF: maks|ΔA| = {zarf_fark:.1e}  maks|Δω| = {om_fark:.1e}",
          flush=True)
    print(f"     sha256(A) vekil = {sha_A}", flush=True)
    print(f"     sha256(A) Hkesk = {sha_Ah}", flush=True)
    print(f"  G2 φ≡0: maks|ΔS| = {g2S:.1e}  maks|ΔS'| = {g2Sp:.1e}", flush=True)
    if not (G1 and G2):
        raise SystemExit("KAPI G1/G2 TUTMADI — kurulmaz.")

    n0 = int(np.ceil(I164.rvm_N(t0)))
    ns = np.arange(n0, n0 + nz, dtype=float) - 0.5
    TOLIN = max(2e-9, 8.0 * np.spacing(float(ns[-1])))
    print(f"  TOLIN (çekirdekle aynı formül) = {TOLIN:.3e}", flush=True)

    I164.S_par = VC.S_par_v
    sar = DerinlikSarmalayici(ns, TOLIN, VC.SSp_par_v)
    I164.SSp_par = sar

    ctx = mp.get_context("spawn")
    pool = ctx.Pool(nwork, initializer=VC._init_v, initargs=(om, A, phi))
    try:
        z, F, tani = I164.coz_sadakatli(om, A, ns, t0, pool, h=h)
    finally:
        pool.close()
        pool.join()

    print(f"\n  sarmalayıcı: {sar.it} yineleme; aktif küme dizisi = "
          f"{sar.nakt}", flush=True)
    print(f"  aktif-küme defteri tutarlı mı: {sar.tutarli}", flush=True)

    # ---------------- R-KAPISI: d=∞ ↔ 176/z_VF1.npy -------------------
    ref_p = S176 / f"z_VF{tohum}.npy"
    R_KAPI = None
    dz_ref = None
    if ref_p.exists():
        zref = np.load(ref_p)
        dz_ref = float(np.max(np.abs(z - zref))) if len(zref) == len(z) else np.inf
        R_KAPI = bool(dz_ref == 0.0)
        print(f"  R-KAPISI: maks|z(∞) − z_VF{tohum}| = {dz_ref:.3e}   "
              f"{'✓ BİT-BİT' if R_KAPI else '✗ TUTMADI'}", flush=True)
    else:
        print(f"  R-KAPISI: {ref_p} yok — sınav yapılamadı", flush=True)

    # ---------------- basamakları yaz ---------------------------------
    out = dict(tohum=int(tohum), h=h, nz=nz, L_hedef=L_hedef,
               sha_A=sha_A, sha_A_Hkeskin=sha_Ah, G1=bool(G1), G2=bool(G2),
               tolin=float(TOLIN), n_yineleme=int(sar.it),
               aktif_dizi=sar.nakt, defter_tutarli=bool(sar.tutarli),
               R_KAPISI=R_KAPI, dz_ref=dz_ref, tani_inf=tani, basamak={})

    def kayit(ad, zz, etiket):
        dz = np.diff(zz)
        sirali = bool(np.all(dz > 0))
        n_ihlal = int(np.sum(dz <= 0))
        zs = np.sort(zz)                       # ölçüm zinciri böyle okur
        dzs = np.diff(zs)
        mid = 0.5 * (zs[:-1] + zs[1:])
        Lw = np.log(mid / TWO_PI)
        L = float(Lw.mean())
        ds = dzs * Lw / TWO_PI - 1
        Sf = None
        r = dict(ad=ad, etiket=etiket, sirali=sirali, n_ihlal=n_ihlal,
                 min_dz=float(dz.min()), L=L,
                 sigma_ds=float(np.std(ds)), sigma_ds2=float(np.var(ds)),
                 n=int(len(zz)),
                 maks_kayma=float(np.max(np.abs(zz - z))),
                 rms_kayma=float(np.sqrt(np.mean((zz - z) ** 2))))
        for dd in (S183, S176, S167, S155):
            np.save(dd / f"z_{ad}.npy", zz)
        print(f"  [{etiket}] {ad}: sıralı={sirali} (ihlal {n_ihlal})  "
              f"min Δz={dz.min():+.6f}  σ_ds={r['sigma_ds']:.5f}  "
              f"L={L:.9f}  maks|z−z(∞)|={r['maks_kayma']:.3e}", flush=True)
        return r

    for dd in sorted(sar.anlik):
        out["basamak"][str(dd)] = kayit(ADLAR[dd], sar.anlik[dd], "TEŞHİS")
    # d=∞ ayrıca kaydedilmez: 176/z_VF1.npy zaten odur (R-KAPISI kanıtı)
    dzi = np.diff(z)
    out["basamak"]["inf"] = dict(ad=f"VF{tohum}", etiket="İNŞA (mevcut)",
                                 sirali=bool(np.all(dzi > 0)),
                                 n_ihlal=int(np.sum(dzi <= 0)),
                                 min_dz=float(dzi.min()),
                                 maks_F=float(np.abs(F).max()),
                                 nF_asan=int(tani["nF_asan"]),
                                 hucre_benzersiz=int(tani["hucre_benzersiz"]),
                                 n=int(len(z)), maks_kayma=0.0, rms_kayma=0.0)
    out["sure_s"] = time.time() - tb
    (S183 / f"DERINLIK_t{tohum}.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False, default=float))
    print(f"\n  -> {S183}/DERINLIK_t{tohum}.json  "
          f"({(time.time()-tb)/60:.1f} dk)", flush=True)
    return out


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
