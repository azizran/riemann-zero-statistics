# -*- coding: utf-8 -*-
"""
198k0 — W_alt ZİNCİRİ (190k0 zincir sırası AYNEN: 155 η → 184b → 185b → 186b; içerik
değişikliği YOK). KALEM_KAPPA_YUKSEKLIK_25EYL2026 "Ölçüm notları": W_alt için 190k0 muadili.
Modüller importlib ile yüklenir; DOSYALAR DÜZENLENMEZ. Çalışma-anı yamaları (190k0 emsali):
  • 155_kos: PENCERE['alt'] = (20000, 200000)  (zeros6 Z[20000:200000]); K.SCR → çıktı dizini
    (η önbelleği 155'in KENDİ üreticisiyle: veri_yukle('alt') + K.eta_onbellek(z,'alt',0.4,4000)).
  • 184b: ONK = onkayit() (SHA ÖNCE), S155/S184 → çıktı; kos('gercek_alt', 'alt') AYNEN
    (pk_m(e^{0.86 L}) evreni, L = pencerenin L'si).
  • 185b/186b: kinematik yaması 190k0 AYNEN (yeni etiket → 185b 'gercek' kolu, pencerenin η'sı);
    oz_hesapla, merdiven_izdusum(ρ≡1) AYNEN.
KÖRLÜK (KALEM 198: W_alt'ta ölçüm-öncesi κ/K/P/A görülmez): 190k0'ın 187c ζ DEFTERİ (çekirdek
niceliği) KOŞULMAZ; 184b'nin ekrana bastığı ölçülen genlik özeti (medyan w) GİZLENİR — yalnız
sayaçlar/kinematik yazılır.
Çıktı: scratchpad/198/zincir_alt/{eta_alt_t0.4_c4000, K1_gercek_alt, OZ_gercek_alt,
       G1_proj_gercek_alt}.npz + zincir_gercek_alt.json
Kullanım: 198k0_zincir.py
"""
import contextlib
import importlib.util
import io
import json
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
S198 = QM / "scratchpad" / "198"
OUT = S198 / "zincir_alt"
PENCERE, ARALIK, ETIKET = "alt", (20000, 200000), "gercek_alt"
TWO_PI = 2 * np.pi


def yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec)
    sys.modules[ad] = m
    spec.loader.exec_module(m)
    return m


def zincir():
    t00 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"{'='*78}\n198k0 ZİNCİR [{PENCERE}] Z[{ARALIK[0]}:{ARALIK[1]}] etiket={ETIKET} "
          f"-> {OUT}\n{'='*78}", flush=True)
    # ---------------- 155 η önbelleği (155'in kendi üreticisi) ----------------
    k155 = yukle("k155", QM / "155_configs" / "155_kos.py")
    assert PENCERE not in k155.PENCERE
    k155.PENCERE[PENCERE] = ARALIK
    k155.K.SCR = OUT
    z = k155.veri_yukle(PENCERE)
    print(f"  [155] {PENCERE}: n={len(z)} t∈[{z[0]:.3f},{z[-1]:.3f}]", flush=True)
    t0 = time.time()
    C = k155.K.eta_onbellek(z, PENCERE, 0.4, 4000)
    mid, ds, L = C["mid"], C["ds"], C["L"]
    g = np.diff(z)
    g2 = (ds + 1.0) * TWO_PI / np.log(mid / TWO_PI)
    tut = float(np.max(np.abs(g - g2)))
    Lw = np.log(mid / TWO_PI)
    print(f"  [155] L={L:.9f} N={len(mid)} L∈[{Lw.min():.4f},{Lw.max():.4f}] nq={C['nq']} "
          f"tutarlılık maks|g−(ds+1)2π/log(m/2π)|={tut:.2e} ({time.time()-t0:.0f}s)", flush=True)
    eta_yol = OUT / f"eta_{PENCERE}_t0.4_c4000.npz"

    # ---------------- 184b (ölçülen genlik özeti ekrana BASILMAZ) ----------------
    b184 = yukle("b184", QM / "184_configs" / "184b_K1_zarf.py")
    b184.ONK = b184.onkayit()
    b184.S155 = OUT
    b184.S184 = OUT
    k1 = OUT / f"K1_{ETIKET}.npz"
    if not k1.exists():
        t0 = time.time()
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            ozet = b184.kos(ETIKET, PENCERE)
        print(f"  [184b] kos({ETIKET}, {PENCERE}) AYNEN: L={ozet['L']:.9f} N={ozet['N']} "
              f"nline={ozet['nline']} ({time.time()-t0:.0f}s; K1 genlik özeti gizlendi)",
              flush=True)
    G = np.load(k1)
    w = np.asarray(G["w"], float)
    tau = np.asarray(G["tau"], float)
    aq = np.asarray(G["aq"], float)
    print(f"  [184b] çizgi evreni τ≤0.86: {len(w)}; HAVUZ' [0.45,0.74): "
          f"{int(((tau >= 0.45) & (tau < 0.74)).sum())}", flush=True)

    # ---------------- 185b / 186b (190k0 kinematik yaması AYNEN) ----------------
    b185 = yukle("b185", QM / "185_configs" / "185b_oz_muhasebe.py")
    b186 = yukle("b186", QM / "186_configs" / "186b_g1_oztutarlilik.py")
    b185.onkayit()
    b186.onkayit()

    def kinematik_yamasi(modul):
        orijinal = modul.kinematik

        def kinematik(gaz):
            if gaz != ETIKET:
                return orijinal(gaz)
            d = np.load(eta_yol)
            mid_ = np.asarray(d["mid"], float)
            ds_ = np.asarray(d["ds"], float)
            L_ = float(d["L"])
            g_ = (ds_ + 1.0) * TWO_PI / np.log(mid_ / TWO_PI)
            return g_, mid_, L_
        modul.kinematik = kinematik

    kinematik_yamasi(b185)
    kinematik_yamasi(b186.b185)
    b185.S185 = OUT
    oz = OUT / f"OZ_{ETIKET}.npz"
    if not oz.exists():
        print("  [185b] oz_hesapla …", flush=True)
        b185.oz_hesapla(ETIKET, w)
    pj = OUT / f"G1_proj_{ETIKET}.npz"
    if not pj.exists():
        print("  [186b] merdiven_izdusum (ρ≡1) …", flush=True)
        rhos = b186.rho_vektorleri(tau)
        b186.merdiven_izdusum(ETIKET, w, aq, {"bir": rhos["bir"]}, pj)
    OZ = np.load(oz)
    P = np.load(pj)
    if not (np.array_equal(OZ["nb"], G["nb"]) and np.array_equal(P["nb"], G["nb"])):
        raise SystemExit("BLOK KENARLARI 184 İLE UYUMSUZ")
    sonuc = {"pencere": PENCERE, "aralik": list(ARALIK), "etiket": ETIKET, "L": float(L),
             "N": int(G["N"]), "L_min": float(Lw.min()), "L_max": float(Lw.max()),
             "t_lo": float(z[0]), "t_hi": float(z[-1]), "nline_tau086": int(len(w)),
             "havuz_ussu_045_074": int(((tau >= 0.45) & (tau < 0.74)).sum()),
             "tutarlilik_g": tut, "nb": G["nb"].tolist(),
             "korluk": "187c ζ defteri KOŞULMADI; 184b genlik özeti gizlendi",
             "sure_s": time.time() - t00}
    json.dump(sonuc, open(OUT / f"zincir_{ETIKET}.json", "w"), indent=1, ensure_ascii=False)
    print(f"\n  -> {OUT/('zincir_'+ETIKET+'.json')}  ({time.time()-t00:.0f}s)  BİTTİ", flush=True)


if __name__ == "__main__":
    zincir()
