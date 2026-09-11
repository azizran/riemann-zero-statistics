# -*- coding: utf-8 -*-
"""
T3 (kalem 188) — KESİM TARAMASI: |ζ|(τ) profilinin minimumu kesimle kayıyor mu?
================================================================================
Üç kesim: τ_c = 0.70, 0.86, 1.00.
  * Minimum τ_c ile KAYIYORSA  → profil KESİM kaynaklı (yapısal, deniz değil)
  * Minimum yerinde (≈0.65) KALIYORSA → profil DENİZ kaynaklı (aritmetik içerik)
0.86 koşusu, saklı G1_proj_gercek ile karşılaştırılarak doğrulanır (yeniden inşa kontrolü).
"""
import sys, json, importlib.util, math
import numpy as np
from pathlib import Path

QM = Path("/Users/ugur/Desktop/Deney/qm_riemann"); SCR = QM / "scratchpad"
S155, S184, S185, S186 = (SCR / d for d in ["155", "184", "185", "186"])
TWO_PI = 2 * np.pi
DISI_QM = "/Users/ugursezen/Desktop/arin/deney/qm_riemann"
DISI_SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
            "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")

def yerel_yukle(ad, yol, ek=()):
    k = Path(yol).read_text(encoding="utf-8").replace(DISI_QM, str(QM)).replace(DISI_SCR, str(SCR))
    for a, b in ek:
        if a not in k: raise SystemExit(f"yama hedefi yok ({ad}): {a[:60]}")
        k = k.replace(a, b)
    spec = importlib.util.spec_from_loader(ad, loader=None)
    m = importlib.util.module_from_spec(spec); m.__file__ = str(yol); sys.modules[ad] = m
    exec(compile(k, str(yol), "exec"), m.__dict__); return m

_kin = [("    zdos = {\"Hkeskin\": \"z_Hkeskin.npy\", \"HA4\": \"z_HA4.npy\"}[gaz]",
  "    if gaz not in (\"Hkeskin\", \"HA4\"):\n"
  "        d = np.load(S155 / f\"eta_{gaz}_t0.4_c4000.npz\")\n"
  "        mid = np.asarray(d[\"mid\"], float); ds = np.asarray(d[\"ds\"], float)\n"
  "        return (ds + 1.0) * TWO_PI / np.log(mid / TWO_PI), mid, float(d[\"L\"])\n"
  "    zdos = {\"Hkeskin\": \"z_Hkeskin.npy\", \"HA4\": \"z_HA4.npy\"}[gaz]")]

b184 = yerel_yukle("b184", QM / "184_configs" / "184b_K1_zarf.py")
b185 = yerel_yukle("b185", QM / "185_configs" / "185b_oz_muhasebe.py", _kin)
b186 = yerel_yukle("b186", QM / "186_configs" / "186b_g1_oztutarlilik.py", [(
    "b185 = importlib.util.module_from_spec(spec)\nspec.loader.exec_module(b185)",
    "b185 = importlib.util.module_from_spec(spec)\nspec.loader.exec_module(b185)\n"
    "_kin_orijinal = b185.kinematik\n"
    "def _kin_genel(gaz):\n"
    "    if gaz in (\"Hkeskin\", \"HA4\"): return _kin_orijinal(gaz)\n"
    "    from pathlib import Path as _P\n"
    "    _S155 = _P(\"/Users/ugur/Desktop/Deney/qm_riemann/scratchpad/155\")\n"
    "    d = np.load(_S155 / f\"eta_{gaz}_t0.4_c4000.npz\")\n"
    "    mid = np.asarray(d[\"mid\"], float); ds = np.asarray(d[\"ds\"], float)\n"
    "    return (ds + 1.0) * TWO_PI / np.log(mid / TWO_PI), mid, float(d[\"L\"])\n"
    "b185.kinematik = _kin_genel")])
b187 = yerel_yukle("b187", QM / "187_configs" / "187c_zeta_defteri.py")
b184.ONK = b184.onkayit()

NJ = 8
def kur(tc):
    """Verilen kesim için (tau, c_olc, c_oz, c_kesik, N, nb) döndürür."""
    T3 = SCR / "T3"; T3.mkdir(parents=True, exist_ok=True)
    if tc <= 0.865:                          # saklı 0.86 evreninden maskele
        D = np.load(S184 / "K1_gercek.npz"); OZ = np.load(S185 / "OZ_gercek.npz")
        tau = np.asarray(D["tau"], float); m = tau <= tc
        w = np.asarray(D["w"], float)[m]; aq = np.asarray(D["aq"], float)[m]; taum = tau[m]
        co = b185.c_olculu(D, -1)[m]; cz = b185.c_oz(OZ, np.asarray(D["aq"], float), -1)[m]
        rh = {k: v[m] for k, v in b186.rho_vektorleri(tau).items()}
        gaz = "son"
    else:                                    # 1.00: K1'i yeniden inşa et
        K1y = T3 / "K1_gercek_tc100.npz"
        if not K1y.exists():
            b184.S184 = T3
            kaynak = (QM / "184_configs" / "184b_K1_zarf.py").read_text(encoding="utf-8")
            # TAU_CIZGI'yı 1.00'a çek (bellekte)
            b184_100 = yerel_yukle("b184_100", QM / "184_configs" / "184b_K1_zarf.py",
                                   [("TAU_CIZGI = 0.86", "TAU_CIZGI = 1.00")])
            b184_100.ONK = b184_100.onkayit(); b184_100.S184 = T3
            b184_100.kos("gercek_tc100", "son")
        D = np.load(T3 / "K1_gercek_tc100.npz")
        tau = np.asarray(D["tau"], float); m = np.ones_like(tau, bool)
        w = np.asarray(D["w"], float); aq = np.asarray(D["aq"], float); taum = tau
        OZy = T3 / "OZ_son.npz"
        if not OZy.exists():
            b185.S185 = T3; b185.oz_hesapla("son", w)   # eta_son üzerinden
        OZ = np.load(OZy)
        co = b185.c_olculu(D, -1); cz = b185.c_oz(OZ, aq, -1)
        rh = b186.rho_vektorleri(tau)
    P = SCR / "T3" / f"proj_tc{int(round(tc*100))}.npz"
    if not P.exists():
        P.parent.mkdir(parents=True, exist_ok=True)
        b186.merdiven_izdusum("son", w, aq, rh, P)
    Pz = np.load(P)
    ck = b187.c_kesik(Pz, -1)
    return taum, co, cz, ck

def tara(taum, co, cz, ck, gen=0.05):
    sat = []
    for c0 in np.arange(0.45, taum.max() + 1e-9, 0.01):
        m = (taum >= c0 - gen/2) & (taum < c0 + gen/2)
        if m.sum() < 5: continue
        delta = co - ck; mix = ck - cz
        z = np.sum(delta[m]*np.conj(mix[m]))/np.sum(np.abs(mix[m])**2)
        sat.append((float(c0), float(abs(z)), float(b187.sar(np.degrees(np.angle(z)))), int(m.sum())))
    return sat

if __name__ == "__main__":
    print(f"1/π = {1/math.pi:.5f}")
    sonuc = {}
    for tc in [0.70, 0.86, 1.00]:
        taum, co, cz, ck = kur(tc)
        sat = tara(taum, co, cz, ck)
        mods = np.array([s[1] for s in sat])
        imin = int(np.argmin(mods))
        sonuc[str(tc)] = sat
        print(f"\nkesim τ_c = {tc:.2f}   ({len(taum)} çizgi, {len(sat)} pencere)")
        print(f"  minimum |ζ| = {mods.min():.4f} @ τ = {sat[imin][0]:.3f}   "
              f"uçlar: {sat[0][1]:.4f} (τ={sat[0][0]:.2f}) … {sat[-1][1]:.4f} (τ={sat[-1][0]:.2f})")
        print(f"  açı aralığı: {min(s[2] for s in sat):.2f}° … {max(s[2] for s in sat):.2f}°")
    json.dump(sonuc, open("24_T3_sonuc.json", "w"), indent=1)
    print("\n-> 24_T3_sonuc.json")
