# -*- coding: utf-8 -*-
"""
T1 (kalem 188) — FAZ-RASTGELE VEKİL DENİZDE ζ
==============================================
176b'nin "vekil" gazları: S_v(t) = -Σ_q a_q sin(ω_q t + φ_q) — genlik zarfı
Hkeskin'le BİREBİR, yalnız fazlar bağımsız rastgele. Eğer ζ ≈ 0.33∠180°
bu denizlerde de çıkıyorsa iptal KİNEMATİK; plasebo tabanına düşüyorsa ARİTMETİK.

Zincir: 184b (K1) → 185b (OZ) → 186b (G1_proj) → 187c (ζ defteri).
Orijinal dosyalar DEĞİŞMEZ; yeni gazlar K1_<gaz>/OZ_<gaz>/G1_proj_<gaz> olarak yazılır.

Kullanım: python3 22_T1_vekil_zeta.py gercek VF5 VS1 VF1
"""
import sys, json, time, importlib.util
import numpy as np
from pathlib import Path

QM = Path("/Users/ugur/Desktop/Deney/qm_riemann"); SCR = QM / "scratchpad"
S155, S184, S185, S186, S187 = (SCR / d for d in ["155", "184", "185", "186", "187"])
TWO_PI = 2 * np.pi
sys.path.insert(0, str(QM / "154_configs")); sys.path.insert(0, str(QM / "185_configs"))

DISI_QM  = "/Users/ugursezen/Desktop/arin/deney/qm_riemann"
DISI_SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
            "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")

# 185b'nin kinematik()'ine, TANIMADIĞI gazlar için eta önbelleğinden okuyan bir dal ekle.
KINEMATIK_YAMA = [
    ('    zdos = {"Hkeskin": "z_Hkeskin.npy", "HA4": "z_HA4.npy"}[gaz]',
     '    if gaz not in ("Hkeskin", "HA4"):\n'
     '        d = np.load(S155 / f"eta_{gaz}_t0.4_c4000.npz")\n'
     '        mid = np.asarray(d["mid"], float); ds = np.asarray(d["ds"], float)\n'
     '        return (ds + 1.0) * TWO_PI / np.log(mid / TWO_PI), mid, float(d["L"])\n'
     '    zdos = {"Hkeskin": "z_Hkeskin.npy", "HA4": "z_HA4.npy"}[gaz]'),
]

def yerel_yukle(ad, yol, ek_yama=()):
    """Kaynağı BELLEKTE yamalar (mutlak yollar + ek yama); diskteki dosya değişmez."""
    kaynak = Path(yol).read_text(encoding="utf-8")
    kaynak = kaynak.replace(DISI_QM, str(QM)).replace(DISI_SCR, str(SCR))
    for a, b in ek_yama:
        if a not in kaynak:
            raise SystemExit(f"YAMA HEDEFİ BULUNAMADI ({ad}): {a[:60]}")
        kaynak = kaynak.replace(a, b)
    spec = importlib.util.spec_from_loader(ad, loader=None)
    m = importlib.util.module_from_spec(spec)
    m.__file__ = str(yol)
    sys.modules[ad] = m
    exec(compile(kaynak, str(yol), "exec"), m.__dict__)
    return m

b184 = yerel_yukle("b184", QM / "184_configs" / "184b_K1_zarf.py")
b185 = yerel_yukle("b185", QM / "185_configs" / "185b_oz_muhasebe.py", KINEMATIK_YAMA)
_186_YAMA = [(
    "b185 = importlib.util.module_from_spec(spec)\nspec.loader.exec_module(b185)",
    "b185 = importlib.util.module_from_spec(spec)\nspec.loader.exec_module(b185)\n"
    "_kin_orijinal = b185.kinematik\n"
    "def _kin_genel(gaz):\n"
    "    if gaz in (\"Hkeskin\", \"HA4\"):\n"
    "        return _kin_orijinal(gaz)\n"
    "    from pathlib import Path as _P\n"
    "    _S155 = _P(\"/Users/ugur/Desktop/Deney/qm_riemann/scratchpad/155\")\n"
    "    d = np.load(_S155 / f\"eta_{gaz}_t0.4_c4000.npz\")\n"
    "    mid = np.asarray(d[\"mid\"], float); ds = np.asarray(d[\"ds\"], float)\n"
    "    return (ds + 1.0) * TWO_PI / np.log(mid / TWO_PI), mid, float(d[\"L\"])\n"
    "b185.kinematik = _kin_genel"
)]
b186 = yerel_yukle("b186", QM / "186_configs" / "186b_g1_oztutarlilik.py", _186_YAMA)
b187 = yerel_yukle("b187", QM / "187_configs" / "187c_zeta_defteri.py")
b184.ONK = b184.onkayit()          # kos() modül-globali ONK'u bekliyor

def zeta_hesapla(gaz, hesapla=False, kuyruk=False):
    K1 = S184 / f"K1_{gaz}.npz"
    if hesapla and not K1.exists():
        t = time.time(); b184.kos(gaz, gaz); print(f"    [K1 {gaz}] {time.time()-t:.0f}s", flush=True)
    D = np.load(K1)
    w = np.asarray(D["w"], float); tau = np.asarray(D["tau"], float)
    aq = np.asarray(D["aq"], float); ae = np.asarray(D["aq_eff"], float)
    OZp = S185 / f"OZ_{gaz}.npz"
    if not OZp.exists():
        t = time.time(); b185.oz_hesapla(gaz, w); print(f"    [OZ {gaz}] {time.time()-t:.0f}s", flush=True)
    Pp = S186 / f"G1_proj_{gaz}.npz"
    if not Pp.exists():
        rhos = b186.rho_vektorleri(tau); b186.merdiven_izdusum(gaz, w, aq, rhos, Pp)
    OZ = np.load(OZp); P = np.load(Pp)
    kenar = json.load(open(S186 / "ONKAYIT_186.json"))["kenar"]
    maskeler = [(tau >= kenar[b]) & (tau < kenar[b+1]) for b in range(8)]
    maskeler.append((tau >= kenar[0]) & (tau < kenar[-1]))
    c_olc = b185.c_olculu(D, -1); c_oz = b185.c_oz(OZ, aq, -1); ck = b187.c_kesik(P, -1)
    delta = c_olc - ck; mix = ck - c_oz
    sat = []
    for m in maskeler:
        z = np.sum(delta[m] * np.conj(mix[m])) / np.sum(np.abs(mix[m]) ** 2)
        sae = ae[m].sum(); soz = np.abs(c_oz[m]).sum()
        m_olc = (np.abs(c_olc[m]).sum() - soz) / sae
        m_kes = (np.abs(ck[m]).sum() - soz) / sae
        sat.append((abs(z), float(b187.sar(np.degrees(np.angle(z)))), 1 - m_olc / m_kes))
    return sat

if __name__ == "__main__":
    gazlar = sys.argv[1:] or ["gercek", "Hkeskin"]
    sonuc = {}
    for gaz in gazlar:
        t = time.time()
        try:
            s = zeta_hesapla(gaz, hesapla=True)
        except Exception as e:
            print(f"  [{gaz}] HATA: {type(e).__name__}: {e}", flush=True); continue
        hav = s[-1]
        sonuc[gaz] = s
        print(f"\n[{gaz}]  ({time.time()-t:.0f}s)", flush=True)
        print(f"  {'bant':>10} {'|ζ|':>8} {'açı°':>9} {'genlik-okuma':>13}")
        for et, (m, a, g) in zip(["0.45-0.50","0.50-0.55","0.55-0.60","0.60-0.65","0.65-0.70",
                                  "0.70-0.75","0.75-0.80","0.80-0.86","HAVUZ"], s):
            print(f"  {et:>10} {m:8.4f} {a:9.2f} {g:13.4f}")
    json.dump({k: [[float(x) for x in row] for row in v] for k, v in sonuc.items()},
              open("22_T1_sonuc.json", "w"), indent=1, ensure_ascii=False)
    print("\n-> 22_T1_sonuc.json")
