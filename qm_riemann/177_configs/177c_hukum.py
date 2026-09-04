# -*- coding: utf-8 -*-
"""
177c — K3 HÜKÜM: DONDURULMUŞ KURALIN YÜZLEŞMESİ
================================================
ÖN-MÜHÜR: bu betik, `tohum = 3` vekilinin HİÇBİR ölçümü var olmadan önce
yazılmıştır (176b VF3 inşası koşarken; 176c hiç koşmamışken).

YENİ ÖLÇÜM YOKTUR. `177/ONKAYIT_177.json`'daki dondurulmuş kural
(SAÇ_n, H1, H2, merdiven) ölçülen θ'lara uygulanır ve hüküm yazılır.

Girdi : 177/ONKAYIT_177.json, 176/G_VF*.json, 176/insa_VF*.json
Çıktı : 177/HUKUM_177.json

Kullanım: 177c_hukum.py [VF1 VF2 VF3 ...]     (öntanımlı: VF1 VF2 VF3)
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S176, S177 = SCR / "176", SCR / "177"

VEK = sys.argv[1:] or ["VF1", "VF2", "VF3"]
ONK = json.load(open(S177 / "ONKAYIT_177.json"))
K = ONK["kural"]
TH_A = K["F3"]["dal_a_yasar"]
TH_B = K["F3"]["dal_b_olur"]
TH_HK = K["F3"]["capa"]["Hkeskin"]
D_SON = K["F3"]["capa"]["dlog_th_son_Hk"]
D_ERFC = K["F3"]["capa"]["dlog_th_erfc_Hk"]
KILIT_TH = K["H2"]["kilit_th"]
KESIM_TH = K["H2"]["kesim_th"]
F_TH = K["H2"]["f_th"]


def sac_n(y):
    y = np.asarray(y, dtype=float)
    return float(2.0 * np.std(y, ddof=1) / math.sqrt(len(y))) if len(y) > 1 \
        else float("nan")


def sac_kat(y):
    y = np.asarray(y, dtype=float)
    return float(y.max() - y.min()) if len(y) > 1 else float("nan")


def karar(ort, sac, esik, yon):
    """(kesin?, marj). yon '<=' dal (a), '>=' dal (b)."""
    marj = abs(ort - esik)
    tarafta = (ort <= esik) if yon == "<=" else (ort >= esik)
    return bool(tarafta and marj >= sac), marj, bool(tarafta)


print("=" * 78)
print("177c — K3 HÜKÜM   [ön-kayıt %s   sha %s]"
      % (ONK["zaman"], ONK["sha256"][:16]))
print("=" * 78)
print("  vekiller: %s   (n = %d)" % (" ".join(VEK), len(VEK)))

GV = {v: json.load(open(S176 / f"G_{v}.json")) for v in VEK}
IV = {v: json.load(open(S176 / f"insa_{v}.json")) for v in VEK}

# ── F0: inşa kapıları (raporlanır; ölüm maddesi ön-kayıtta) ────────────
print("\n" + "=" * 78)
print("F0 — İNŞA KAPILARI (G1–G5) ve R_bant (BAĞLAYICI DEĞİL)")
print("=" * 78)
print("%-5s %-9s %10s %10s %9s %9s %9s %8s"
      % ("gaz", "tohum", "ilk-kök", "maks|F|", "sıralı", "zarf", "φ≡0",
         "R_bant"))
kapi_hepsi = True
for v in VEK:
    t = IV[v]
    kp = t["kapilar"]
    ok = all(kp.values())
    kapi_hepsi &= ok
    print("%-5s %-9s %10s %10.3e %9s %9.1e %9.1e %8.4f  %s"
          % (v, t.get("tohum"), "%d/300000" % t["hucre_benzersiz"],
             t["maxF"], "TAM" if t["sirali"] else "BOZUK",
             t["zarf_fark"], t["g2_dS"], GV[v]["R_bant_min"],
             "✓" if ok else "✗ KAPI ÖLÜMÜ"))
print("  sha256(A): %s"
      % (" | ".join("%s %s" % (v, IV[v]["sha_A"][:16]) for v in VEK)))
print("  hepsi Hkeskin'in sha'sıyla aynı mı: %s"
      % all(IV[v]["sha_A"] == IV[v]["sha_A_Hkeskin"] for v in VEK))
print("  NOT: R_bant kapısı 176'da ıskalamıştı ve kurtarılmadı; ön-kayıt "
      "onu 177'de de\n       HİÇBİR HÜKME DAYANAK YAPMAZ diye yazdı.")

# ── H1: F3 / H-F1b ─────────────────────────────────────────────────────
print("\n" + "=" * 78)
print("H1 — F3 / H-F1b:  (a) ≤ %.7f YAŞAR | (b) ≥ %.7f ÖLÜR"
      % (TH_A, TH_B))
print("=" * 78)
th = [GV[v]["th"] for v in VEK]
ort = float(np.mean(th))
s_n, s_k = sac_n(th), sac_kat(th)
print("  θ = %s" % "  ".join("%s %.7f" % (v, x) for v, x in zip(VEK, th)))
print("  ȳ = %.7f    SAÇ_n = 2s/√n = %.7f    [katı max−min = %.7f]"
      % (ort, s_n, s_k))
ka, ma, ta = karar(ort, s_n, TH_A, "<=")
kb, mb, tb = karar(ort, s_n, TH_B, ">=")
ka_k = bool(ta and ma >= s_k)
print("  dal (a): taraf %s  marj %.7f  vs SAÇ_n %.7f  ⇒ %s"
      % ("✓" if ta else "✗", ma, s_n, "KESİN" if ka else "HÜKÜMSÜZ"))
print("           katı okumada: marj %.7f vs %.7f ⇒ %s"
      % (ma, s_k, "KESİN" if ka_k else "hükümsüz"))
print("  dal (b): taraf %s  marj %.7f  vs SAÇ_n %.7f  ⇒ %s"
      % ("✓" if tb else "✗", mb, s_n, "KESİN" if kb else "HÜKÜMSÜZ"))
if ka:
    H1 = "YAŞADI"
elif kb:
    H1 = "ÖLDÜ"
else:
    H1 = "HÜKÜMSÜZ"
print("  ⇒ H1 = **H-F1b %s**" % H1)

# ── H2: F9 / ω_θ bandı ─────────────────────────────────────────────────
print("\n" + "=" * 78)
print("H2 — F9 / ω_θ bandı:  ω_θ ∈ (0,1]  ⟺  Δ̄ ≥ kilit_θ = %.7f"
      % KILIT_TH)
print("=" * 78)
dl = [math.log(TH_HK) - math.log(x) for x in th]
dbar = float(np.mean(dl))
ds_n, ds_k = sac_n(dl), sac_kat(dl)
om = [KILIT_TH / d if d != 0 else float("inf") for d in dl]
om_bar = KILIT_TH / dbar
print("  Δ_i = log θ(Hk) − log θ(VF_i) = %s"
      % "  ".join("%+.7f" % x for x in dl))
print("  Δ̄ = %+.7f   SAÇ_n(Δ) = %.7f   [katı %.7f]" % (dbar, ds_n, ds_k))
print("  ω_θ (tohum başına) = %s" % "  ".join("%.4f" % x for x in om))
print("  ω_θ (Δ̄ üzerinden) = %+.4f" % om_bar)
ust_k, ust_m, ust_t = karar(dbar, ds_n, KILIT_TH, ">=")     # ω ≤ 1
alt_k, alt_m, alt_t = karar(dbar, ds_n, 0.0, ">=")          # ω > 0
ust_kk = bool(ust_t and ust_m >= ds_k)
print("  üst kenar (ω ≤ 1 ⟺ Δ̄ ≥ %.7f): taraf %s marj %.7f vs %.7f ⇒ %s"
      % (KILIT_TH, "✓" if ust_t else "✗", ust_m, ds_n,
         "KESİN" if ust_k else "HÜKÜMSÜZ"))
print("            katı okumada: marj %.7f vs %.7f ⇒ %s"
      % (ust_m, ds_k, "KESİN" if ust_kk else "hükümsüz"))
print("  alt kenar (ω > 0 ⟺ Δ̄ > 0):        taraf %s marj %.7f vs %.7f ⇒ %s"
      % ("✓" if alt_t else "✗", alt_m, ds_n,
         "KESİN" if alt_k else "HÜKÜMSÜZ"))
H2 = ("KESİN: ω_θ ∈ (0,1]" if (ust_k and alt_k)
      else ("BANDIN İÇİNDE ama kesin değil" if (ust_t and alt_t)
            else "BANDIN DIŞINDA"))
print("  ⇒ H2 = **%s**" % H2)

# ── F9'un θ satırı (parametresiz; kesim/kilit tohumdan bağımsız) ───────
print("\n" + "=" * 78)
print("F9 — NİHAİ AYRIŞIM TABLOSUNUN θ SATIRI (güncel)")
print("=" * 78)
print("  Δ(son←Hk) = %+.6f   Δ(erfc←Hk) = %+.6f   f_θ = %+.4f"
      % (D_SON, D_ERFC, F_TH))
print("  KESİM_θ = %+.6f   KİLİT_θ = %+.6f   (tohumdan BAĞIMSIZ)"
      % (KESIM_TH, KILIT_TH))
print("  Δ(Hk←vekil) = %+.6f  (tohum başına %s)"
      % (dbar, " / ".join("%+.5f" % x for x in dl)))
print("  ω_θ = KİLİT_θ / Δ(Hk←vekil) = %+.4f" % om_bar)
print("  gerçeğin fazlasının yüzdesi: Δ(Hk←vekil)/Δ(son←Hk) = %%%.1f"
      % (100.0 * dbar / D_SON))

# ── merdiven ───────────────────────────────────────────────────────────
n = len(VEK)
if H1 != "HÜKÜMSÜZ" and (ust_k and alt_k):
    sonraki = "TAMAM — merdiven durur."
elif n < K["merdiven"]["ust_sinir"]:
    sonraki = "n ← %d (VF%d, tohum %d) kur." % (n + 1, n + 1, n + 1)
else:
    sonraki = "n = %d üst sınırdır ⇒ KALICI HÜKÜMSÜZ." % n
print("\n  MERDİVEN: %s" % sonraki)

OUT = dict(onkayit=ONK["zaman"], sha=ONK["sha256"], vekiller=VEK, n=n,
           F0=dict(kapilar={v: IV[v]["kapilar"] for v in VEK},
                   hepsi_gecti=bool(kapi_hepsi),
                   R_bant_min={v: GV[v]["R_bant_min"] for v in VEK},
                   R_bant_baglayici=False),
           H1=dict(theta={v: GV[v]["th"] for v in VEK}, ort=ort,
                   sac_n=s_n, sac_kat=s_k, esik_a=TH_A, esik_b=TH_B,
                   marj_a=ma, marj_b=mb, taraf_a=ta, taraf_b=tb,
                   kesin_a=ka, kesin_a_kati=ka_k, kesin_b=kb, hukum=H1),
           H2=dict(delta={v: d for v, d in zip(VEK, dl)}, dbar=dbar,
                   sac_n=ds_n, sac_kat=ds_k, kilit_th=KILIT_TH,
                   omega_tohum={v: o for v, o in zip(VEK, om)},
                   omega=om_bar, ust_kesin=ust_k, ust_kesin_kati=ust_kk,
                   alt_kesin=alt_k, hukum=H2),
           F9_theta=dict(d_son_Hk=D_SON, d_erfc_Hk=D_ERFC, f=F_TH,
                         kesim=KESIM_TH, kilit=KILIT_TH,
                         d_Hk_vekil=dbar, omega=om_bar,
                         yuzde_fazla=100.0 * dbar / D_SON),
           merdiven=sonraki)
p = S177 / ("HUKUM_177_n%d.json" % n)
p.write_text(json.dumps(OUT, indent=1, ensure_ascii=False))
print("  -> %s" % p)
