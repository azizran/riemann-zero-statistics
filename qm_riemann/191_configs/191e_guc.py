# -*- coding: utf-8 -*-
"""
191e — K4: KİLİT GÜÇ ANALİZİ (deep-twin KİLİT ÖLÇÜMÜNDEN ÖNCE)
================================================================
ONKAYIT_191["K4_guc_analizi"]'nde donan kural: 180 raporunun HAM
gerçek↔Hkeskin farkları (Δlog M = +0.057542, Δlog Q_E = −0.044475) için,
D=1.30'da iki rakip okumanın ("kilit de derinlik" vs "kilit tabanı")
AYRIMINI, f_1.30 (191d) ve 180'in kilit payı s_kilit (± kendi se'si)
kullanarak hesaplar — ÖLÇÜMDEN (yeni gaz inşası) ÖNCE.

Cebirsel özdeşlik (iki okumanın farkı):
  oku_i  (kilit de derinlik) : Δ(D) = Δ(1.00)·(1 − f_D)
  oku_ii (kilit tabanı)      : Δ(D) = Δ(1.00)·[(1−s_kilit)(1−f_D) + s_kilit]
  Ayrım(D) := oku_ii − oku_i = Δ(1.00)·s_kilit·f_D

σ_birleşik: Δ(1.00) NOKTA DEĞER kabul edilir (180'in kendi formülü: ΔN
"tohumsuz", ayrı bir se YOK; scratchpad/180 bu oturumda TEMİZLENMİŞ,
yeniden hesaplanamaz — AÇIKÇA yazılan sınırlama). s_kilit VE f_D'nin
KENDİ belirsizlikleri bağımsız yayılır:
  σ_birleşik = |Ayrım| · √((se_s_kilit/s_kilit)² + (se_f_D/f_D)²)

Karar: Ayrım/σ_birleşik ≥ 3 (HEM M HEM Q_E'de) VE 180'in para birimi
hattı (167_olcum→172b) yeni bir gaz için 60 dakikadan kısa sürede
koşulabiliyorsa → ölç ve hükmet; aksi hâlde ÖLÇME YAPILMAZ (kurtarma yok).

Girdi: 191/ONKAYIT_191.json, 191/HUKUM_191.json.
Çıktı: 191/GUC_191.json (+ ekran defteri).
"""
import hashlib
import json
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S189, S191 = SCR / "189", SCR / "191"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


if __name__ == "__main__":
    print("=" * 78)
    print("191e / K4 — KİLİT GÜÇ ANALİZİ (ölçümden ÖNCE)")
    print("=" * 78, flush=True)

    ONK = json.load(open(S191 / "ONKAYIT_191.json"))
    if sha(QM / "191_configs" / "191a_onkayit.py") != ONK["sha256"]:
        raise SystemExit("ON-KAYIT SHA UYUMSUZ")
    K4 = ONK["K4_guc_analizi"]
    HUK = json.load(open(S191 / "HUKUM_191.json"))
    f_D = HUK["H-191a"]["f130_havuz"]
    f_D_se = HUK["H-191a"]["f130_se"]
    print(f"  girdi f_1.30 (191d, HAVUZ, ortak-blok loo se): "
          f"{f_D:.4f}±{f_D_se:.4f}")

    girdi = K4["girdi_180"]
    sonuc = {}
    for ad, dlog, dlog_ad in (
        ("M", girdi["Delta_log_M_son_Hk"], "s_kilit_M_log"),
        ("Q_E", girdi["Delta_log_Q_E_son_Hk"], "s_kilit_QE_log"),
    ):
        s_kilit = girdi[dlog_ad]["merkez"]
        se_s = girdi[dlog_ad]["se"]
        oku_i = dlog * (1 - f_D)
        oku_ii = dlog * ((1 - s_kilit) * (1 - f_D) + s_kilit)
        ayrim = dlog * s_kilit * f_D
        ayrim_check = oku_ii - oku_i
        rel = np.sqrt((se_s / s_kilit) ** 2 + (f_D_se / f_D) ** 2)
        sig_birlesik = abs(ayrim) * rel
        oran = abs(ayrim) / sig_birlesik if sig_birlesik > 0 else float("inf")
        sonuc[ad] = dict(
            Delta_log_1_00=dlog, s_kilit=s_kilit, se_s_kilit=se_s,
            f_D=f_D, se_f_D=f_D_se,
            oku_kilit_de_derinlik=oku_i, oku_kilit_tabani=oku_ii,
            ayrim=ayrim, ayrim_ozdeslik_kontrolu=ayrim_check,
            rel_belirsizlik=rel, sigma_birlesik=sig_birlesik,
            ayrim_sigma_orani=oran, yeterli_guc=bool(oran >= 3.0))
        print(f"\n  [{ad}] Δ(1.00)={dlog:+.6f}  s_kilit={s_kilit:.3f}"
              f"±{se_s:.3f}")
        print(f"    oku_i  (kilit de derinlik) Δ(1.30) = {oku_i:+.6f}")
        print(f"    oku_ii (kilit tabanı)      Δ(1.30) = {oku_ii:+.6f}")
        print(f"    Ayrım = Δ(1.00)·s_kilit·f_D = {ayrim:+.6f}  "
              f"(özdeşlik kontrolü: oku_ii−oku_i = {ayrim_check:+.6f}, "
              f"fark {abs(ayrim-ayrim_check):.2e})")
        print(f"    σ_birleşik = |Ayrım|·{rel:.4f} = {sig_birlesik:.6f}")
        print(f"    Ayrım/σ_birleşik = {oran:.3f}  "
              f"({'YETERLİ (≥3)' if oran >= 3.0 else 'YETERSİZ (<3)'})")

    guc_yeterli_ikisi_de = all(sonuc[a]["yeterli_guc"] for a in sonuc)
    print(f"\n  GÜÇ HÜKMÜ: M {sonuc['M']['ayrim_sigma_orani']:.2f}σ, "
          f"Q_E {sonuc['Q_E']['ayrim_sigma_orani']:.2f}σ — ikisi de ≥3σ mü: "
          f"{guc_yeterli_ikisi_de}")

    # ---- ikincil kapı: 180 hattının altyapısı hâlâ erişilebilir mi ----
    S166_183 = [d for d in range(166, 184) if (SCR / str(d)).exists()]
    hat_erisilebilir = len(S166_183) > 0
    if hat_erisilebilir:
        durum_metni = "MEVCUT " + str(S166_183)
    else:
        durum_metni = ("YOK (bu oturumda temizlenmiş — 167/C_*.json, "
                       "172b girdileri, 176 VF vekilleri vb. hiçbiri yok)")
    print(f"\n  180 para-birimi hattının altyapısı (scratchpad/166-183): "
          f"{durum_metni}")
    print("  (Bu ikincil kapı yalnız BİLGİ amaçlı kontrol edildi: güç "
         "hükmü zaten YETERSİZ çıktığı için '60 dakika' sorusu pratikte "
         "ele alınmaz — kurala göre HER İKİ kapı da açık olmalı.)")

    if guc_yeterli_ikisi_de:
        karar = ("(varsayımsal — bu dalda tetiklenmedi) güç yeterli, "
                "60-dk hat sınavı gerekirdi")
    else:
        karar = "ÖLÇME YAPILMAZ: güç yetersiz (Ayrım/σ_birleşik < 3, M ve/veya Q_E)"
    print(f"\n  KARAR: {karar}")
    print("  Kurtarma yok, eşik gevşetme yok.")

    out = dict(f_1_30=f_D, f_1_30_se=f_D_se, sonuc=sonuc,
              guc_yeterli_ikisi_de=guc_yeterli_ikisi_de,
              hat_altyapisi_mevcut_dizinler=S166_183,
              hat_erisilebilir_bilgi_amacli=hat_erisilebilir,
              karar=karar,
              sinirlama=("Δ(1.00) (ΔlogM, ΔlogQ_E) 180 raporunda NOKTA "
                        "değer; kendi ayrı jackknife se'si YOK (180a "
                        "formülü: ΔN 'tohumsuz') ve kaynak scratchpad "
                        "(166-183) bu oturumda temizlenmiş olduğu için "
                        "yeniden hesaplanamadı — σ_Δ(1.00)=0 varsayıldı; "
                        "bu YETERSİZ güç sonucunu SIKILAŞTIRAN değil "
                        "GEVŞETEN bir varsayımdır (gerçek σ_birleşik "
                        "muhtemelen daha büyük, oran muhtemelen daha "
                        "küçük olurdu) — yani 'yetersiz güç' hükmü bu "
                        "sınırlamadan ETKİLENMEZ / tersine daha da "
                        "sağlamlaşır."))
    json.dump(out, open(S191 / "GUC_191.json", "w"), indent=1,
              ensure_ascii=False, default=float)
    print(f"\n-> {S191/'GUC_191.json'}  BİTTİ", flush=True)
