# -*- coding: utf-8 -*-
"""
189b — K1: DERİN İKİZ İNŞASI (164_insa.py DÜZENLENMEDEN, sarmalayıcıyla)
=======================================================================
ONKAYIT_189 dondurdu: 164_insa modülü importlib ile yüklenir; `merdiven`
modül özniteliği functools.partial(merdiven, tau_ust=D) ile değiştirilir
(main onu çağrı anında modül globalinden okur); KONFIG'e
Hderin110/Hderin120 = ("sadakatli", {"c": -0.5}) eklenir; main(ad, h, nz)
AYNEN çağrılır (çıktı: 164/z_<ad>.npy, 155/z_<ad>.npy, 164/insa_<ad>.json).

Süreç havuzu 'spawn': işçiler modülü '164_insa' adıyla sys.path üzerinden
içe aktarır (işçiler yalnız _init/_work_S/_work_SP koşar; om, a initargs ile
gelir — yamaya ihtiyaçları yoktur).

Kullanım:
  189b_insa.py izgara 1.10     # ızgara sınavı (ilk 5000 seviye, h .015 vs .0075)
  189b_insa.py insa 1.10       # inşa (h sınav JSON'undan otomatik) + kapılar
"""
import functools
import importlib
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S189 = SCR / "189"
sys.path.insert(0, str(QM / "164_configs"))
insa = importlib.import_module("164_insa")      # dosya DEĞİŞMEZ

AD = {"1.10": "Hderin110", "1.20": "Hderin120"}
N_SINAV = 5000
H_KABA, H_INCE = 0.015, 0.0075
FARK_TOL = 1e-6
FARK_ESIK = 5
L_KAPI = 12.02959324


def yamala(tau_ust):
    """Çalışma-anı yaması: merdiven(…, tau_ust=D) + KONFIG girişi."""
    ad = AD[f"{tau_ust:.2f}"]
    if not hasattr(insa, "_merdiven_orijinal"):
        insa._merdiven_orijinal = insa.merdiven
    insa.merdiven = functools.partial(insa._merdiven_orijinal, tau_ust=tau_ust)
    insa.KONFIG[ad] = ("sadakatli", {"c": -0.5})
    return ad


def pencere():
    d = np.load(QM / "128_odl_zeros6_2e6_zeros.npz")
    Z = np.sort(np.asarray(d["zeros"], dtype=float))
    zr = Z[len(Z) - 300000:]
    t0, t1 = float(zr[0]), float(zr[-1])
    L_hedef = float(np.log(0.5 * (t0 + t1) / insa.TWO_PI))
    return t0, t1, L_hedef


def izgara_sinavi(tau_ust):
    ad = yamala(tau_ust)
    tb = time.time()
    t0, t1, L_hedef = pencere()
    om, a, _, _ = insa.merdiven(L_hedef, None, None)
    print(f"=== 189b IZGARA SINAVI {ad} (τ_üst={tau_ust}) ===", flush=True)
    print(f"merdiven: {len(om)} çizgi (τ≤{tau_ust:.2f}, L_hedef={L_hedef:.4f}); "
          f"ω_max={om.max():.4f} → en kısa dalga boyu 2π/ω_max="
          f"{2*np.pi/om.max():.4f}; Σa={a.sum():.3f}", flush=True)
    n0 = int(np.ceil(insa.rvm_N(t0)))
    ns = np.arange(n0, n0 + N_SINAV, dtype=float) - 0.5
    ctx = mp.get_context("spawn")
    pool = ctx.Pool(insa.NWORK, initializer=insa._init, initargs=(om, a))
    try:
        zk, Fk, tk = insa.coz_sadakatli(om, a, ns, t0, pool, h=H_KABA)
        zi, Fi, ti = insa.coz_sadakatli(om, a, ns, t0, pool, h=H_INCE)
    finally:
        pool.close()
        pool.join()
    dz = np.abs(zk - zi)
    nfark = int((dz > FARK_TOL).sum())
    h_sec = H_KABA if nfark <= FARK_ESIK else H_INCE
    sonuc = dict(ad=ad, tau_ust=tau_ust, nline=int(len(om)), n_sinav=N_SINAV,
                 n0=n0, h_kaba=H_KABA, h_ince=H_INCE, fark_tol=FARK_TOL,
                 n_fark=nfark, maks_dz=float(dz.max()),
                 medyan_dz=float(np.median(dz)),
                 maksF_kaba=tk["maxF"], maksF_ince=ti["maxF"],
                 sirali_kaba=bool(np.all(np.diff(zk) > 0)),
                 sirali_ince=bool(np.all(np.diff(zi) > 0)),
                 esik=FARK_ESIK, h_secilen=h_sec, sure_s=time.time() - tb)
    (S189 / f"izgara_{ad}.json").write_text(json.dumps(sonuc, indent=1))
    print(f"IZGARA SINAVI {ad}: ilk-kök farkı {nfark}/{N_SINAV} "
          f"(|Δz|>{FARK_TOL:g}); maks|Δz|={dz.max():.3e}  → h={h_sec} "
          f"({time.time()-tb:.0f}s)", flush=True)


def insa_et(tau_ust):
    ad = yamala(tau_ust)
    sj = S189 / f"izgara_{ad}.json"
    if not sj.exists():
        raise SystemExit(f"ızgara sınavı yok: {sj}")
    h = json.loads(sj.read_text())["h_secilen"]
    print(f"=== 189b İNŞA {ad}: merdiven tau_ust={tau_ust} (164 main'in "
          f"'τ≤1.00' etiketi sabit metindir; gerçek kesim τ≤{tau_ust:.2f}), "
          f"h={h} (ızgara sınavından) ===", flush=True)
    insa.main(ad, h, insa.NZERO)
    tani = json.loads((SCR / "164" / f"insa_{ad}.json").read_text())
    kapilar = {
        "maks_F": tani["maxF"], "maks_F_gecti": tani["maxF"] <= 1e-8,
        "sirali": tani["sirali"], "min_dz": tani["min_dz"],
        "L": tani["L"], "L_8hane": round(tani["L"], 8),
        "L_gecti": round(tani["L"], 8) == L_KAPI,
        "nline": tani["nline"], "h": tani["h"], "sure_s": tani["sure_s"],
        "sigma_ds2": tani["sigma_ds2"], "n": tani["n"],
    }
    kapilar["HEPSI"] = bool(kapilar["maks_F_gecti"] and kapilar["sirali"]
                            and kapilar["L_gecti"])
    (S189 / f"insa_kapi_{ad}.json").write_text(json.dumps(kapilar, indent=1))
    print(f"İNŞA KAPILARI {ad}: maks|F|={tani['maxF']:.3e} "
          f"({'≤1e-8' if kapilar['maks_F_gecti'] else '>1e-8 KALDI'}), "
          f"sıralılık {'TAM' if tani['sirali'] else 'BOZUK'}, "
          f"L={tani['L']:.9f} ({'=' if kapilar['L_gecti'] else '≠'}12.02959324) "
          f"→ {'GEÇTİ' if kapilar['HEPSI'] else 'KALDI'}", flush=True)


if __name__ == "__main__":
    mod, D = sys.argv[1], float(sys.argv[2])
    if f"{D:.2f}" not in AD:
        raise SystemExit("derinlik 1.10 ya da 1.20")
    if mod == "izgara":
        izgara_sinavi(D)
    elif mod == "insa":
        insa_et(D)
    else:
        raise SystemExit("kullanım: 189b_insa.py <izgara|insa> <1.10|1.20>")
