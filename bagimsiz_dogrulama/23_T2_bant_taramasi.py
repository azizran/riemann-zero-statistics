# -*- coding: utf-8 -*-
"""
T2 (kalem 188) — KAYAN BANT TARAMASI
=====================================
Sabit 8 bant yerine τ boyunca kayan pencere: |ζ|(τ) gerçekten DÜZ mü?
  * Düz ve 1/π'ye sabit  → plato KİNEMATİK bir sabit olabilir (tehlike)
  * Düzgün değişiyor     → ζ bir izdüşüm katsayısı (dinamik)
Kontrol: aynı tarama faz-rastgele vekilde (VF1) — "gürültü nasıl görünür".
"""
import sys, json, importlib.util
import numpy as np
from pathlib import Path

QM = Path("/Users/ugur/Desktop/Deney/qm_riemann"); SCR = QM / "scratchpad"
S184, S185, S186, S187 = (SCR / d for d in ["184", "185", "186", "187"])
TWO_PI = 2 * np.pi
DISI_QM = "/Users/ugursezen/Desktop/arin/deney/qm_riemann"
DISI_SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
            "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")

def yerel_yukle(ad, yol, ek=()):
    k = Path(yol).read_text(encoding="utf-8").replace(DISI_QM, str(QM)).replace(DISI_SCR, str(SCR))
    for a, b in ek:
        if a not in k: raise SystemExit(f"yama hedefi yok: {a[:50]}")
        k = k.replace(a, b)
    spec = importlib.util.spec_from_loader(ad, loader=None)
    m = importlib.util.module_from_spec(spec); m.__file__ = str(yol); sys.modules[ad] = m
    exec(compile(k, str(yol), "exec"), m.__dict__); return m

_kin_yama = [("    zdos = {\"Hkeskin\": \"z_Hkeskin.npy\", \"HA4\": \"z_HA4.npy\"}[gaz]",
  "    if gaz not in (\"Hkeskin\", \"HA4\"):\n"
  "        d = np.load(S155 / f\"eta_{gaz}_t0.4_c4000.npz\")\n"
  "        mid = np.asarray(d[\"mid\"], float); ds = np.asarray(d[\"ds\"], float)\n"
  "        return (ds + 1.0) * TWO_PI / np.log(mid / TWO_PI), mid, float(d[\"L\"])\n"
  "    zdos = {\"Hkeskin\": \"z_Hkeskin.npy\", \"HA4\": \"z_HA4.npy\"}[gaz]")]
b185 = yerel_yukle("b185", QM / "185_configs" / "185b_oz_muhasebe.py", _kin_yama)
b187 = yerel_yukle("b187", QM / "187_configs" / "187c_zeta_defteri.py")

NJ = 8
def yukle_uclu(gaz):
    D = np.load(S184 / f"K1_{gaz}.npz"); OZ = np.load(S185 / f"OZ_{gaz}.npz"); P = np.load(S186 / f"G1_proj_{gaz}.npz")
    tau = np.asarray(D["tau"], float); aq = np.asarray(D["aq"], float)
    co = [b185.c_olculu(D, i) for i in range(-1, NJ)]
    cz = [b185.c_oz(OZ, aq, i) for i in range(-1, NJ)]
    ck = [b187.c_kesik(P, i) for i in range(-1, NJ)]
    return tau, co, cz, ck

def zeta_hesap(co, cz, ck, m):
    delta = co - ck; mix = ck - cz
    z = np.sum(delta[m] * np.conj(mix[m])) / np.sum(np.abs(mix[m]) ** 2)
    return abs(z), float(b187.sar(np.degrees(np.angle(z))))

def tara(gaz, genislik=0.05, adim=0.01):
    tau, co, cz, ck = yukle_uclu(gaz)
    sat = []
    for c0 in np.arange(0.45, 0.861, adim):
        m = (tau >= c0 - genislik/2) & (tau < c0 + genislik/2)
        if m.sum() < 5: continue
        tam = zeta_hesap(co[0], cz[0], ck[0], m)
        reps = np.array([zeta_hesap(co[i+1], cz[i+1], ck[i+1], m) for i in range(NJ)])
        se_m = np.sqrt((NJ-1)/NJ * np.sum((reps[:,0]-reps[:,0].mean())**2))
        df = b187.sar(reps[:,1]-tam[1]); se_a = np.sqrt((NJ-1)/NJ*np.sum((df-df.mean())**2))
        sat.append((float(c0), tam[0], se_m, tam[1], se_a, int(m.sum())))
    return sat

if __name__ == "__main__":
    import math
    print(f"1/π = {1/math.pi:.5f}    (gerçek deniz havuzu |ζ| = 0.3287)")
    for gaz in ["gercek", "VF1"]:
        for gen in [0.05, 0.10]:
            sat = tara(gaz, gen)
            mods = np.array([s[1] for s in sat])
            print(f"\n[{gaz}] genişlik {gen:.2f} — {len(sat)} pencere")
            print(f"  |ζ|: min {mods.min():.4f}  medyan {np.median(mods):.4f}  maks {mods.max():.4f}  "
                  f"std {mods.std():.4f}  (ort. hata {np.mean([s[2] for s in sat]):.4f})")
            print(f"  {'τ merkez':>9} {'|ζ|':>8} {'±se':>7} {'açı°':>8} {'±se':>6} {'çizgi':>6}")
            for (c, m, sm, a, sa, n) in sat[::4]:
                print(f"  {c:9.3f} {m:8.4f} {sm:7.4f} {a:8.2f} {sa:6.2f} {n:6d}")
            json.dump(sat, open(f"23_T2_{gaz}_{gen:.2f}.json", "w"), indent=1)
