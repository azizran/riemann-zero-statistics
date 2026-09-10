# -*- coding: utf-8 -*-
"""
187c — K2: ζ DEFTERİ (iptal katsayısı + açısı; gerçek VE ikiz)
==============================================================
ONKAYIT_187 dondurdu:
  Δ_q = c^ölç_q − c^kesik_q   (c^ölç 184 K1 blokları; c^kesik 186 G1_proj
  re_bir/im_bir AYNEN).  karışım^kesik_q = c^kesik_q − c^öz_q (ρ≡1 öz-terim,
  185 OZ blokları AYNEN).
  ζ_b = Σ_{q∈b} Δ_q·conj(karışım_q) / Σ_{q∈b} |karışım_q|²
  modül = iptal payı; açı (derece) faz sınavı; jackknife açı-se loo
  replika açılarının tam-örneklem açısına sarılmasıyla.
  H-F2: açı ∈ 180°±15° (bant başına) → 'saf yıkıcı'; dışı → döndürülmüş
  faz-örgüsü KAYDI. Yan sütun: genlik-okuması 1 − m_ölç/m^kesik.

Çıktı: 187/K2_zeta.json + ekran defteri.
"""
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S184, S185, S186, S187 = (SCR / d for d in ["184", "185", "186", "187"])
NJACK = 8

spec = importlib.util.spec_from_file_location(
    "b185", QM / "185_configs" / "185b_oz_muhasebe.py")
b185 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b185)


def onkayit():
    o = json.load(open(S187 / "ONKAYIT_187.json"))
    sha = hashlib.sha256(
        (QM / "187_configs" / "187a_onkayit.py").read_bytes()).hexdigest()
    if sha != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {sha} != {o['sha256']}")
    return o


def c_kesik(P, disari=-1):
    N = int(P["N"])
    nb = P["nb"]
    if disari < 0:
        re, im, n = P["re_bir"].sum(0), P["im_bir"].sum(0), N
    else:
        re = P["re_bir"].sum(0) - P["re_bir"][disari]
        im = P["im_bir"].sum(0) - P["im_bir"][disari]
        n = N - nb[disari]
    return 2.0 * (re + 1j * im) / n


def sar(a):
    """Açıyı (−180, 180]'e sar."""
    return (a + 180.0) % 360.0 - 180.0


if __name__ == "__main__":
    ONK = onkayit()
    print("=" * 78)
    print(f"187c / K2 ζ DEFTERİ  [on-kayit {ONK['zaman']} "
          f"sha {ONK['sha256'][:12]}]")
    print("=" * 78, flush=True)

    G = np.load(S184 / "K1_gercek.npz")
    tau = np.asarray(G["tau"], float)
    aq = np.asarray(G["aq"], float)
    ae = np.asarray(G["aq_eff"], float)
    kenar_b = json.load(open(S186 / "ONKAYIT_186.json"))["kenar"]
    maskeler = [(tau >= kenar_b[b]) & (tau < kenar_b[b + 1]) for b in range(8)]
    maskeler.append((tau >= kenar_b[0]) & (tau < kenar_b[-1]))   # HAVUZ
    et = [f"{kenar_b[b]:.2f}-{kenar_b[b+1]:.2f}" for b in range(8)] + ["HAVUZ"]

    sonuc = {"sha_onkayit": ONK["sha256"], "etiket": et}
    for gaz, Kad, OZad, Pad in [
            ("gercek", "K1_gercek.npz", "OZ_gercek.npz",
             "G1_proj_gercek.npz"),
            ("Hkeskin", "K1_Hkeskin.npz", "OZ_Hkeskin.npz",
             "G1_proj_Hkeskin.npz")]:
        D = np.load(S184 / Kad)
        OZ = np.load(S185 / OZad)
        P = np.load(S186 / Pad)

        def zeta_defteri(disari):
            c_olc = b185.c_olculu(D, disari)
            c_oz = b185.c_oz(OZ, aq, disari)
            ck = c_kesik(P, disari)
            delta = c_olc - ck
            mix = ck - c_oz
            out = np.zeros((9, 4))
            for k, m in enumerate(maskeler):
                z = np.sum(delta[m] * np.conj(mix[m])) / \
                    np.sum(np.abs(mix[m]) ** 2)
                sae = ae[m].sum()
                soz = np.abs(c_oz[m]).sum()
                m_olc = (np.abs(c_olc[m]).sum() - soz) / sae
                m_kes = (np.abs(ck[m]).sum() - soz) / sae
                out[k] = [np.abs(z), np.degrees(np.angle(z)),
                          m_olc, m_kes]
            return out

        tam = zeta_defteri(-1)
        reps = np.array([zeta_defteri(i) for i in range(NJACK)])
        se_mod = b185.jk_se(reps[:, :, 0])
        # açı-se: replika açıları tam-örneklem açısına sarılır
        dfark = sar(reps[:, :, 1] - tam[None, :, 1])
        se_aci = np.sqrt((NJACK - 1) / NJACK *
                         np.sum((dfark - dfark.mean(0)) ** 2, 0))
        genlik = 1.0 - tam[:, 2] / tam[:, 3]     # 1 − m_ölç/m^kesik
        g_reps = 1.0 - reps[:, :, 2] / reps[:, :, 3]
        se_genlik = b185.jk_se(g_reps)
        icinde = np.abs(sar(tam[:, 1] - 180.0)) <= 15.0

        print(f"\n[{gaz}] ζ DEFTERİ (±jk se):")
        print(f"{'bant':>10} {'|ζ|':>16} {'açı°':>16} "
              f"{'180°±15?':>9} {'genlik-okuma':>15}")
        for k in range(9):
            print(f"{et[k]:>10} {tam[k,0]:7.4f}±{se_mod[k]:.4f} "
                  f"{tam[k,1]:+8.2f}±{se_aci[k]:5.2f} "
                  f"{'EVET' if icinde[k] else 'HAYIR':>9} "
                  f"{genlik[k]:7.4f}±{se_genlik[k]:.4f}")
        sonuc[gaz] = {
            "zeta_mod": tam[:, 0].tolist(), "se_mod": se_mod.tolist(),
            "zeta_aci": tam[:, 1].tolist(), "se_aci": se_aci.tolist(),
            "pencere_180pm15": icinde.tolist(),
            "genlik_okuma": genlik.tolist(),
            "se_genlik": se_genlik.tolist(),
            "m_olc": tam[:, 2].tolist(), "m_kesik": tam[:, 3].tolist()}

    # H-F2 hükmü (gerçek; bant bant)
    ic = sonuc["gercek"]["pencere_180pm15"][:8]
    print(f"\nH-F2 (gerçek, 8 bant): 180°±15° içinde {sum(ic)}/8 bant → "
          f"{'SAF YIKICI mührü' if all(ic) else 'döndürülmüş faz-örgüsü KAYDI (pencere dışı bantlar var)'}")
    sonuc["HF2"] = {"bant_icinde": ic, "hepsi_icinde": bool(all(ic))}
    json.dump(sonuc, open(S187 / "K2_zeta.json", "w"), indent=1,
              ensure_ascii=False)
    print(f"\n-> {S187/'K2_zeta.json'}  BİTTİ", flush=True)
