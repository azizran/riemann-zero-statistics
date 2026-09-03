"""
170 — K3: BACAK-BACAK KORELATÖR AKTARIMI, λ EKSENİ ve GERÇEK GAZ
================================================================
169 §K2.3'ün BÖLÜNMÜŞ-MERDİVEN kurulumu (`169_k2b.main`) AYNEN çağrılır
(kod kopyalanmaz); tek yaptığımız onu Hkeskin'in DIŞINDAKİ gazlarda da
koşturup sonuçları yan yana koymaktır:

    G = E · X_a · X_b   (üç BAĞIMSIZ alan bacağı + taşıyıcı)
    ρ(varyant) = Pu2(kırpılmış)/Pu2(tam),  clip(F) = σ_F sgn(F)

H-D1'in korelatör-düzeyi sınavı: eğer `c(λ)` bacak-başına aktarım
yasasıysa, ÖLÇÜLEN bacak-başı aktarım λ ile `c(λ)^{1/4}` gibi düşmelidir:
    c^{1/4}: 0.7970 (λ=1.00) → 0.7861 (0.85) → 0.7795 (0.70)   [−%2.2]
K3'ün sorusu (H-D2): gerçek gazın (`son`) bacak aktarımları Hkeskin'den
BÜYÜK mü (162: kilitli fazlar η'yı %29 sessizleştirir ⇒ daha az doyum
⇒ daha büyük aktarım ⇒ daha büyük c)?

ÖN-MÜHÜR (koşudan ÖNCE):
  * Hkeskin (169'da ölçüldü): 1 bacak 0.8293, 2 bacak 0.6586,
    3 bacak 0.5284; alan düzeyi E 0.82775, Xa 0.80854, Xb 0.80915.
  * H-D1 tutuyorsa L070'in 3-bacak aktarımı
    0.5284·(0.3690/0.4035)^{3/4} = **0.4945** olmalı (−%6.4).
  * 170b §K1.5 alan düzeyinde λ-değişmezlik ölçtü (‰1–5) ⇒ BEKLENTİM
    korelatör düzeyinde de aktarımın λ ile ~SABİT kalması (|Δ| < %2),
    yani H-D1'in −%6.4'ünün ÇIKMAMASI.
  * `son` için: alan-düzeyi aktarım Hkeskin'e çok yakın çıkmalı
    (170b: E 0.8268 ↔ 0.8296, X 0.8225 ↔ 0.8214) ⇒ K3'ün +3.06σ'sı
    bacak aktarımıyla KAPANMAMALI.

SONUÇ bloğu YALNIZ gerçek koşu çıktısındandır.
Kullanım: 170e_bacak_aktarim.py [gaz1,gaz2,...]
Çıktı:    scratchpad/169/K2b_<gaz>.json (169'un düzeni) + ekrana özet
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("170_configs", "169_configs", "167_configs", "166_configs",
           "165_configs", "163_configs", "160_configs", "159_configs",
           "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
ORT = importlib.import_module("167_ortak")
ORT.KUNYE["L060"] = dict(gercek=False, lam=0.60, tau_ust=1.00, pen=None,
                         aile="lam", T=1.0)
K2B = importlib.import_module("169_k2b")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
G1 = np.sqrt(2 / np.pi)
NB = {"100_E": 1, "010_Xa": 1, "001_Xb": 1,
      "011_XaXb": 2, "110_EXa": 2, "111_hepsi": 3}


def ozet(gazlar):
    print("\n" + "=" * 104)
    print("K3 — BACAK-BACAK KORELATÖR AKTARIMI (bütün gazlar yan yana)")
    print("=" * 104)
    print(f"  {'gaz':9s} {'1 bacak':>9s} {'2 bacak':>9s} {'3 bacak':>9s} "
          f"{'(3b)^{4/3}':>10s} {'ölçülen c':>10s} {'oran':>7s}")
    C = {a: json.load(open(SCR / f"167/C_{a}.json")) for a in gazlar
         if (SCR / f"167/C_{a}.json").exists()}
    out = {}
    for g in gazlar:
        p = SCR / f"169/K2b_{g}.json"
        if not p.exists():
            continue
        d = json.loads(p.read_text())
        o = d["ortalama"]
        v = {n: float(np.mean([o[a] for a in o if NB[a] == n]))
             for n in (1, 2, 3)}
        cc = float(np.exp(np.mean(np.log(
            [b["KALIB_u2"] / b["W_X"] for b in C[g]["bant"]
             if b.get("olculdu") and 0.52 - 1e-9 <= b["lo"] <= 0.68 + 1e-9
             and b["tau_eff"] < 0.85 and b["Ms2"] > 0
             and b["R_bant"] >= 0.98 and b["SNR"] >= 3.0]))))
        ext = v[3] ** (4.0 / 3.0)
        out[g] = dict(bacak=v, c=cc, ext4=ext, oran=ext / cc)
        print(f"  {g:9s} {v[1]:9.4f} {v[2]:9.4f} {v[3]:9.4f} {ext:10.4f} "
              f"{cc:10.4f} {ext/cc:7.4f}")
    print(f"  {'Gauss':9s} {G1:9.4f} {G1**2:9.4f} {G1**3:9.4f} "
          f"{G1**4:10.4f}   (4/π²)")
    return out


if __name__ == "__main__":
    gz = (sys.argv[1].split(",") if len(sys.argv) > 1
          else ["L085", "L070", "son"])
    for g in gz:
        if not (SCR / f"169/K2b_{g}.json").exists():
            K2B.main(g, 0.40, 0.95)
    ozet(["Hkeskin"] + [g for g in gz if g != "Hkeskin"])
