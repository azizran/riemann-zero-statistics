# -*- coding: utf-8 -*-
"""
191c — K2: ZİNCİR (Hderin130) — 189c makinesi AYNEN, ince sarmalayıcı
=======================================================================
189c_zincir.py DOSYASI DEĞİŞTİRİLMEDEN importlib ile yüklenir. Tek yeni
şey: b189c.zincir("Hderin130") çağrısı — 184b2 → 185b → 186b → 187c
zincirini YENİ gaz için koşturur (K1/OZ/G1_proj/zincir_Hderin130.* YENİ
dosyalar; mevcut Hkeskin/Hderin110/Hderin120/gercek dosyalarının hiçbiri
ÜZERİNE YAZILMAZ — bu yüzden b189c.zincir("Hkeskin") ve
b189c.gercek_zeta() BURADA TEKRAR ÇAĞRILMAZ; 189'un zaten geçmiş 6/6
makine mührü (scratchpad/189/muhur_Hkeskin.json) ve zincir_gercek.npz
OLDUĞU GİBİ okunur/referans alınır — aynı DEĞİŞTİRİLMEMİŞ kod (sha256
denetimi) bit-bit aynı sonucu vereceği için yeniden koşmak yalnız
üzerine-yazma riski katardı.

Kullanım: 191c_zincir.py
Çıktı: scratchpad/189/{K1,OZ,G1_proj,zincir}_Hderin130.{npz,json} (YENİ),
       scratchpad/191/ZINCIR_191_ozet.json
"""
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S189, S191 = SCR / "189", SCR / "191"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec)
    sys.modules[ad] = m
    spec.loader.exec_module(m)
    return m


if __name__ == "__main__":
    t0 = time.time()
    print("=" * 78)
    print("191c — K2: ZİNCİR (Hderin130; 189c makinesi AYNEN)")
    print("=" * 78, flush=True)

    ONK191 = json.load(open(S191 / "ONKAYIT_191.json"))
    sha_189c_now = sha(QM / "189_configs" / "189c_zincir.py")
    sha_189d_now = sha(QM / "189_configs" / "189d_hukum.py")
    ok_189c = sha_189c_now == ONK191["zincir"]["sha_189c"]
    ok_189d = sha_189d_now == ONK191["zincir"]["sha_189d"]
    print(f"  makine değişmedi mi (K0'daki sha256 ile bit-bit): "
          f"189c={'EVET' if ok_189c else 'HAYIR — DUR'}  "
          f"189d={'EVET' if ok_189d else 'HAYIR — DUR'}")
    if not (ok_189c and ok_189d):
        raise SystemExit("189c/189d K0'dan beri DEĞİŞMİŞ — zincir ŞÜPHELİ, "
                         "koşulmaz")

    muhur = json.load(open(S189 / "muhur_Hkeskin.json"))
    print(f"  189'un mevcut makine mührü (Hkeskin, bit-bit; TEKRAR "
          f"KOŞULMADI — aynı değişmemiş kod): "
          f"{'6/6 TUTTU' if muhur['HEPSI'] else 'TUTMADI'}")
    for k, v in muhur.items():
        if k != "HEPSI":
            print(f"    {k:32s} {'TUTTU' if v else 'TUTMADI'}")
    if not muhur["HEPSI"]:
        raise SystemExit("MEVCUT Hkeskin MÜHRÜ TUTMAMIŞ — zincir koşulmaz")

    zg_var = (S189 / "zincir_gercek.npz").exists()
    print(f"  scratchpad/189/zincir_gercek.npz mevcut (TEKRAR KOŞULMADI): "
          f"{zg_var}")
    if not zg_var:
        raise SystemExit("zincir_gercek.npz yok — K3 modelleri için gerekli")

    print("\n  189c_zincir.py importlib ile yükleniyor (dosya AYNEN)…",
          flush=True)
    b189c = yukle("b189c_191", QM / "189_configs" / "189c_zincir.py")

    hedef_k1 = S189 / "K1_Hderin130.npz"
    if hedef_k1.exists():
        raise SystemExit(f"{hedef_k1} ZATEN VAR — üzerine yazılmaz, betik "
                         "yalnız YENİ gaz için koşar")

    print("\n  b189c.zincir('Hderin130') çağrılıyor "
          "(184b2 → 185b → 186b → 187c) …", flush=True)
    tam, se, sig_r, tau_bar, Z = b189c.zincir("Hderin130")

    js = json.load(open(S189 / "zincir_Hderin130.json"))
    ozet = {
        "zaman": time.strftime("%a %b %d %H:%M:%S %z %Y"),
        "sha_189c_dogrulandi": ok_189c, "sha_189d_dogrulandi": ok_189d,
        "muhur_Hkeskin_referans": muhur["HEPSI"],
        "L_Hderin130": js["L"], "N_Hderin130": js["N"],
        "HAVUZ_index_185": 9, "HAVUZ_index_187": 8,
        "r_HAVUZ": js["r"][9], "sig_r_HAVUZ": js["sig_r"][9],
        "w_g_HAVUZ": js["w_g"][9], "w_HD_HAVUZ": js["w_HD"][9],
        "w_oz_HD_HAVUZ": js["w_oz_HD"][9], "m_HD_HAVUZ": js["m_HD"][9],
        "sigma_eps": js["sigma_eps"], "se_sigma_eps": js["se_sigma_eps"],
        "zeta_mod_HAVUZ": js["zeta_mod"][8], "se_mod_HAVUZ": js["se_mod"][8],
        "zeta_aci_HAVUZ": js["zeta_aci"][8],
        "genlik_okuma_HAVUZ": js["genlik_okuma"][8],
        "sure_zincir_s": js["sure_s"], "sure_toplam_s": time.time() - t0,
    }
    json.dump(ozet, open(S191 / "ZINCIR_191_ozet.json", "w"), indent=1,
              ensure_ascii=False, default=float)
    print(f"\n  HAVUZ: r_1.30={ozet['r_HAVUZ']:.4f}±{ozet['sig_r_HAVUZ']:.4f}"
          f"  σ_ε(1.30)={ozet['sigma_eps']:.5f}±{ozet['se_sigma_eps']:.5f}"
          f"  |ζ_H130|={ozet['zeta_mod_HAVUZ']:.4f}±"
          f"{ozet['se_mod_HAVUZ']:.4f}")
    print(f"\n-> {S191/'ZINCIR_191_ozet.json'}")
    print(f"BİTTİ ({time.time()-t0:.0f}s)", flush=True)
