# -*- coding: utf-8 -*-
"""
190e — ÖN-KAYITSIZ KEŞİF (hüküm DIŞI; yalnız KAYIT)
===================================================
Hükümler (190c) görüldükten SONRA sorulan sorular; hiçbir hükmü değiştirmez:
 (a) Profil bütünü: düşük ω-dilim profili ile son profili arasındaki Pearson
     korelasyonu iki koordinatta — Δω (evrensel okuma) ve τ'−1 ölçekli
     (τ'-sabit okuma: son profili Δω·L_son/L_düşük'te okunur). Aralık [−1.3, 2.3].
 (b) ω-dilim tepelerindeki sistematik pozitif kayma: blok başına Bragg kayması
     vs blok-içi L yayılımı ve blok-içi L dağılımının çarpıklığı
     (medyan L − L_b; ağırlık: noktalar indekse göre düzgün).
 (c) Diğer belirgin tepeler (±0.4 çifti): iki pencerede Δω konumu (aynı kural,
     merkez ±0.15; aday etiket ±log(3/2) = ±0.405 — SINANMADI).
Çıktı: 190/K_kesif_190.json + ekran.
"""
import json
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S190, S155 = SCR / "190", SCR / "155"
TWO_PI = 2 * np.pi
P = np.load(S190 / "profiller_190.npz")
H = json.load(open(S190 / "HUKUM_190.json"))
ONK = json.load(open(S190 / "ONKAYIT_190.json"))
L_D = ONK["pencere"]["L"]
L_S = ONK["L_son"]
out = {}

m = P["om_merkez"]
nb = P["om_nb"].astype(float)
yog_d = (P["om_kap"] * (nb / nb.sum())[:, None]).sum(0) / 0.025
nbs = P["nb_s"].astype(float)


def son_yog(x):
    y = np.zeros_like(x)
    for b in range(8):
        y += (nbs[b] / nbs.sum()) * np.interp(x, P["ts_dw"][b],
                                              P["ts_kap"][b] / (0.005 * L_S),
                                              left=0.0, right=0.0)
    return y


# (a)
msk = (m >= -1.3) & (m <= 2.3)
c_dw = float(np.corrcoef(yog_d[msk], son_yog(m[msk]))[0, 1])
c_tau = float(np.corrcoef(yog_d[msk], son_yog(m[msk] * L_S / L_D))[0, 1])
out["a_profil_corr_Δω"] = c_dw
out["a_profil_corr_tau_olcekli"] = c_tau
print(f"(a) profil korelasyonu (düşük ω-dilim vs son, Δω∈[−1.3,2.3], "
      f"{msk.sum()} dilim): Δω ekseninde {c_dw:.3f} | τ'-ölçekli {c_tau:.3f}")

# (b)
mid = np.load(S190 / "zincir_dusuk" / "eta_dusuk_t0.4_c4000.npz")["mid"]
kj = np.linspace(0, len(mid), 9).astype(int)
Lb = np.array(ONK["bloklar"]["L_b"])
yay = np.array(ONK["bloklar"]["blok_ici_L_yayilimi"])
carp = np.array([float(np.median(np.log(mid[kj[b]:kj[b + 1]] / TWO_PI)) - Lb[b])
                 for b in range(8)])
tab = H["H190a_tablo"]
print("(b) blok-içi L dağılımı ve ω-dilim tepe kaymaları (tepe_b − hedef):")
b_ = {"yayilim": yay.tolist(), "medyanL_eksi_Lb": carp.tolist()}
for ad in ("0", "+log2", "+log3", "+log6", "-log2", "-log3"):
    kay = np.array(tab[ad]["tepe_blok"]) - tab[ad]["hedef"]
    r = float(np.corrcoef(kay, yay)[0, 1])
    b_[ad] = {"kayma": kay.tolist(), "ort_kayma": float(kay.mean()),
              "corr_kayma_yayilim": r}
    print(f"   {ad:>6}: kayma_b = {' '.join('%+.3f' % k for k in kay)} "
          f"ort {kay.mean():+.4f}  corr(kayma, yayılım) = {r:+.2f}")
print(f"   yayılım = {np.round(yay, 3).tolist()}")
print(f"   medyan L − L_b = {np.round(carp, 4).tolist()}")
tum = np.concatenate([np.array(tab[a]["tepe_blok"]) - tab[a]["hedef"]
                      for a in ("0", "+log2", "+log3", "+log6", "-log2", "-log3")])
b_["tum_kayma_ort"] = float(tum.mean())
b_["tum_kayma_se"] = float(tum.std(ddof=1) / np.sqrt(len(tum)))
print(f"   6 hedef × 8 blok ortalama kayma = {tum.mean():+.4f} ± "
      f"{tum.std(ddof=1)/np.sqrt(len(tum)):.4f} (naif se)")
out["b_kayma"] = b_

# (c)
c_ = {}
for ad, h in (("+0.405", np.log(1.5)), ("-0.405", -np.log(1.5))):
    td, ts = [], []
    for b in range(8):
        ok = (np.abs(m - h) <= 0.15 + 1e-9) & (P["om_ncz"][b] > 0)
        td.append(float(m[ok][np.argmax(P["om_kap"][b][ok])]))
        ok = np.abs(P["ts_dw"][b] - h) <= 0.15 + 1e-9
        ts.append(float(P["ts_dw"][b][ok][np.argmax(P["ts_kap"][b][ok])]))
    c_[ad] = {"dusuk_om": td, "dusuk_medyan": float(np.median(td)),
              "son": ts, "son_medyan": float(np.median(ts)),
              "tau_sabit_ong_dusuk": float(np.median(ts)) * L_D / L_S}
    print(f"(c) {ad} çevresi: düşük medyan {np.median(td):+.4f} | son medyan "
          f"{np.median(ts):+.4f} | τ'-sabit okuma düşükte "
          f"{np.median(ts)*L_D/L_S:+.4f}")
out["c_cift"] = c_
json.dump(out, open(S190 / "K_kesif_190.json", "w"), indent=1, ensure_ascii=False)
print(f"-> {S190/'K_kesif_190.json'}  BİTTİ")
