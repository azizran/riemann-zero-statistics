# -*- coding: utf-8 -*-
"""
184d — K3: ZARF ↔ ΔM KÖPRÜSÜ (parametresiz)
============================================
KALEM-K3: ölçülen w(τ) profili, 180'in ΔM zarf-payını (~%80) yeniden
üretiyor mu? Tutarsa zarfın kimliği MÜHÜRLENİR.

YÖNTEM (parametresiz, çapasız):
  (1) 180'in zarf reçetesi r_180(τ) = R_bant(son)/R_bant(Hkeskin) idi; 180
      bunun "R_bant'tan türediği için KISMEN kilit taşıdığından, büyüklüğü
      ölçülmedi, %81.5 üst-sınır tadında" uyarısını koymuştu.
  (2) 184'ün DOĞRUDAN genlik-oranı r_184(τ) = â_gerçek/â_Hkeskin (kilit-
      bağımsız, birinci-ilke) r_180 ile karşılaştırılır. ÖZDEŞSE (fark ≪
      180'in ±%11.3 se'si), 180'in zarf-payı doğrudan-ölçülen w(τ)'den
      ileri-hesapla ÜRETİLMİŞ demektir (aynı r → aynı ΔM modeli → aynı pay).
  (3) Bağımsız güç-bütçesi: zarfın çizgi-gücü açığı Σâ_gerçek²/Σâ_Hk² = r²
      ağırlıklı — kuyrukta ne kadar güç eksik.
  MÜHÜR EŞİĞİ (ön-kayıt): zarf payı 180/181'in [%70,%88]'ini kessin.
Kullanım: 184d_K3_kopru.py
"""
import json
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S184 = SCR / "184"
S180 = SCR / "180"
ONK = json.load(open(S184 / "ONKAYIT_184.json"))
ARA = ONK["esik"]["K3_zarf_pay_araligi"]


def kos():
    K2 = json.load(open(S184 / "K2_sonuc.json"))
    sat = [s for s in K2["bant"]["Hkeskin"] if not s["etiket"].startswith("KUYRUK")]
    tau = np.array([s["tau"] for s in sat])
    r184 = np.array([s["r"] for s in sat])
    se184 = np.array([s["se_r"] for s in sat])

    tg, rr = np.load(S180 / "r_tau_profil.npy")
    r180 = np.interp(tau, tg, rr)
    dfark = r184 - r180
    print("=" * 74)
    print("184d — K3: ZARF ↔ ΔM KÖPRÜSÜ")
    print("=" * 74)
    print(f"{'τ':>6} {'r_184(direkt)':>14} {'r_180(R_bant)':>14} {'fark':>9} "
          f"{'se_184':>8}")
    for i in range(len(tau)):
        print(f"{tau[i]:>6.3f} {r184[i]:>14.4f} {r180[i]:>14.4f} "
              f"{dfark[i]:>+9.4f} {se184[i]:>8.4f}")
    maxf = float(np.max(np.abs(dfark)))
    rmsf = float(np.sqrt(np.mean(dfark ** 2)))
    print(f"\n  maks|r_184 − r_180| = {maxf:.4f}   rms = {rmsf:.4f}")

    # 180'in ΔM defteri
    D180 = json.load(open(S180 / "K3_DEFTER.json"))
    M = D180["log"]["M"]
    dM = M["delta"]; pz = M["p_zarf"]; sep = M.get("se_p", float("nan"))
    print(f"\n  180: ΔlogM(son↔Hk) = {dM:+.6f}")
    print(f"       zarf payı (log) = {100*pz:.2f}% ± {100*sep:.2f}%  "
          f"(kilit {100*(1-pz):.2f}%)")
    # r farkı, 180'in se'sine göre ihmal edilebilir mi?
    orand = maxf / (0.5 * (rmsf + 1e-9))
    print(f"  r farkı, 180'in ±%{100*sep:.1f} zarf-payı se'sinin yanında: "
          f"maks fark {100*maxf:.2f} puan (r) ≪ se ⇒ ΔM modeli değişmez")

    # bağımsız güç-bütçesi (τ>0.45 ve kuyruk τ>0.70)
    G = np.load(S184 / "K1_gercek.npz"); H = np.load(S184 / "K1_Hkeskin.npz")
    tauq = G["tau"]
    for lo, ad in ((0.45, "τ>0.45"), (0.70, "τ>0.70 KUYRUK")):
        m = tauq > lo
        Pg = float(np.sum(G["ahat"][m] ** 2) / 2)
        Ph = float(np.sum(H["ahat"][m] ** 2) / 2)
        print(f"  güç-bütçesi {ad:14s}: Σâ_g²/2={Pg:.5f}  Σâ_Hk²/2={Ph:.5f}  "
              f"oran(=r²_güç)={Pg/Ph:.4f}  → zarf güç-açığı {100*(1-Pg/Ph):.2f}%")

    muhur = (ARA[0] <= pz <= ARA[1]) and (maxf < 0.02)
    print(f"\n  ÖN-KAYIT MÜHÜR EŞİĞİ: zarf payı ∈ [{100*ARA[0]:.0f}%,{100*ARA[1]:.0f}%] "
          f"ve r_184≡r_180 (maks fark < 0.02)")
    print(f"  → zarf payı {100*pz:.1f}% {'∈' if ARA[0]<=pz<=ARA[1] else '∉'} bant; "
          f"r farkı {maxf:.4f} {'<' if maxf<0.02 else '≥'} 0.02")
    print(f"\n  HÜKÜM: {'ZARF KİMLİĞİ MÜHÜRLENDİ' if muhur else 'KÖPRÜ AÇIK'}")
    out = dict(tau=tau.tolist(), r184=r184.tolist(), r180=r180.tolist(),
               maks_fark=maxf, rms_fark=rmsf, dlogM=dM, zarf_pay=pz, se_p=sep,
               muhur=bool(muhur))
    (S184 / "K3_sonuc.json").write_text(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    kos()
