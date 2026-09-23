# -*- coding: utf-8 -*-
"""
189f — ÖN-KAYITSIZ EK (hüküm DIŞI, KAYIT): H-189b sapmasının karışım payı ve
H-189d açığının derinlik-eşli kıyası. Veriler görüldükten SONRA yazıldı;
hiçbir hükmü değiştirmez.

(a) Karışım payının iki parçası: m_HD = m^kesik_HD(1 − Γ_HD) [187c özdeşliği],
    m_öng = m^kesik_öng(1 − z_D), m^kesik_öng = m_Hk/(1 − z_1.00) ⇒
    Δm = (m^kesik_HD − m^kesik_öng)(1 − z_D)   [pencere-içi karışım kinematiği]
       + m^kesik_HD·(z_D − Γ_HD)              [iptal açığı: ölçülen iptal < harita]
    D = 1.00 satırı lehçe tabanını gösterir (iki terim tam iptal eder).
(b) |ζ_HD| ile derinlik-eşli harita toplamları: ikiz haritası (Hkeskin
    kinematiği, 188c) ve gerçek haritası (188b), Σ_{orta ≤ D} K(HAVUZ).
Çıktı: 189/EK_189f.json
"""
import json
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S188, S189 = SCR / "188", SCR / "189"
DER = {"1.00": "Hkeskin", "1.10": "Hderin110", "1.20": "Hderin120"}
BANT = ["0.45-0.50", "0.50-0.55", "0.55-0.60", "0.60-0.65",
        "0.65-0.70", "0.70-0.75", "0.75-0.80", "0.80-0.86"]
IH, IH7 = 9, 8
NJACK = 8


def jk(r):
    r = np.asarray(r, float)
    return float(np.sqrt((NJACK - 1) / NJACK * np.sum((r - r.mean(0)) ** 2, 0)))


Z = {D: json.load(open(S189 / f"zincir_{ad}.json")) for D, ad in DER.items()}
hk = np.load(S188 / "harita_K_Hkeskin.npz", allow_pickle=True)
hg = np.load(S188 / "harita_K_gercek.npz", allow_pickle=True)


def zmap(h, et, D):
    i = list(h["etiket"]).index(et)
    sel = h["orta"] <= D + 1e-12
    return (float(-np.real(h["K"][i, sel].sum())),
            jk(-np.real(h["K_reps"][:, i, sel].sum(1))),
            float(np.abs(h["K"][i, sel].sum())))


out = {"a_karisim": {}, "b_derinlik_esli": {}}
print("(a) KARIŞIM PAYI AYRIŞIMI Δm = Δm^kesik·(1−z_D) + m^kesik_HD·(z_D − Γ_HD)")
for k, et in enumerate(BANT + ["HAVUZ"]):
    r5 = k if k < 8 else IH
    r7 = k if k < 8 else IH7
    lab = "B" + et if et != "HAVUZ" else "HAVUZ"
    z100 = zmap(hk, lab, 1.00)[0]
    mHk = Z["1.00"]["m_HD"][r5]
    mk_ong = mHk / (1 - z100)
    satir = {}
    for D in DER:
        zD = zmap(hk, lab, float(D))[0]
        mk = Z[D]["m_kesik"][r7]
        G = Z[D]["genlik_okuma"][r7]
        m = Z[D]["m_HD"][r5]
        m_ong = mk_ong * (1 - zD)
        t1 = (mk - mk_ong) * (1 - zD)
        t2 = mk * (zD - G)
        satir[D] = dict(m_HD=m, m_ong=m_ong, dm=m - m_ong, kesik_terim=t1,
                        iptal_acigi_terim=t2, kapanis=m - m_ong - t1 - t2,
                        m_kesik_HD=mk, m_kesik_ong=mk_ong, Gamma_HD=G, z_D=zD)
    out["a_karisim"][et] = satir
    if et in ("0.45-0.50", "0.65-0.70", "0.80-0.86", "HAVUZ"):
        for D, s in satir.items():
            print(f"  {et:>10} D={D}: Δm={s['dm']:+.4f} = kesik {s['kesik_terim']:+.4f}"
                  f" + iptal-açığı {s['iptal_acigi_terim']:+.4f}  "
                  f"(m^kesik {s['m_kesik_HD']:.4f} vs öng {s['m_kesik_ong']:.4f}; "
                  f"Γ_HD {s['Gamma_HD']:.4f} vs z_D {s['z_D']:.4f}; kapanış {s['kapanis']:.1e})")

print("\n(b) |ζ_HD| vs derinlik-eşli harita toplamları (HAVUZ):")
for D in DER:
    zi, zise, zim = zmap(hk, "HAVUZ", float(D))
    zg, zgse, zgm = zmap(hg, "HAVUZ", float(D))
    zeta = Z[D]["zeta_mod"][IH7]
    out["b_derinlik_esli"][D] = dict(
        zeta_HD=zeta, se=Z[D]["se_mod"][IH7], ikiz_harita_negRe=zi,
        ikiz_harita_se=zise, ikiz_harita_mod=zim, gercek_harita_negRe=zg,
        gercek_harita_se=zgse, gercek_harita_mod=zgm,
        oran_ikiz=zeta / zi, oran_gercek=zeta / zg)
    print(f"  D={D}: |ζ_HD|={zeta:.4f}±{Z[D]['se_mod'][IH7]:.4f}  ikiz-harita "
          f"{zi:.4f}±{zise:.4f} (oran {zeta/zi:.3f})  gerçek-harita "
          f"{zg:.4f}±{zgse:.4f} (oran {zeta/zg:.3f})")
json.dump(out, open(S189 / "EK_189f.json", "w"), indent=1, ensure_ascii=False)
print(f"\n-> {S189/'EK_189f.json'}")
