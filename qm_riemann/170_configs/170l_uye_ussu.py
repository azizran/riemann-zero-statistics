"""
170 — K0(d): ÜYE ÜSSÜ ν'nün İKİ YÖNDE AYRIŞTIRILMASI
====================================================
Yeni ölçüm YOK; `scratchpad/167/C_<gaz>.json` okunur (bant seçimi
`170c_onkayit60.bant_c` = 169'un ortak penceresi, import).

`c = KALİB_u2/W_X` üyesinin doğruluğu şu soruya indirgenir:
    **KALİB_u2 ∝ W_X^ν'de ν kaçtır?**   ν = 1 ⇒ c gerçekten değişmez.
İki BAĞIMSIZ yön vardır ve karıştırılırsa yanıltır:
  (a) BANT yönü — aynı gazda τ arttıkça W_X düşer (DW eğrisi boyunca);
  (b) λ yönü    — aynı bantta gaz değiştikçe W_X değişir (σ_X̃ ile).
169 §K1.5'in üye sınavı (a)+(b)'yi karıştırıyordu; burada ayrılıyor.

ÖN-MÜHÜR (koşudan ÖNCE): 170c'nin havuzlanmış (a+b) uyumu ν = 0.823
verdi. Beklentim: (a) bant yönü ν'yü 1'e YAKIN verecek (KALİB bant
boyunca W_X ile birlikte düşüyor), (b) λ yönü ise 1'den KÜÇÜK ve λ ile
değişken çıkacak (λ→1'de ≈0.2, λ→0.6'da ≈1.0), yani P4'ün başarısı
büyük ölçüde BANT yönünden geliyor olacak.

Çıktı: scratchpad/170/UYE_USSU.json
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
OK60 = importlib.import_module("170c_onkayit60")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/170")
LAM = [("Hkeskin", 1.00), ("L085", 0.85), ("L070", 0.70), ("L060", 0.60)]
if (Path(str(SCR).replace("/170", "/167")) / "C_L115.json").exists():
    LAM = [("L115", 1.15)] + LAM


def main():
    REF = {a: OK60.bant_c(a)[0] for a, _ in LAM}
    los = sorted(set.intersection(*[set(REF[a]) for a, _ in LAM]))
    print("=" * 92)
    print("K0(d) — ÜYE ÜSSÜ  ν = d log KALİB_u2 / d log W_X  "
          "(BANT yönü ↔ λ yönü)")
    print("=" * 92)

    print("\n  (a) BANT YÖNÜ (her gazda ayrı ayrı, beş bant boyunca):")
    print(f"  {'gaz':9s} {'λ':>5s} {'ν_bant':>8s} {'artık rms':>10s}")
    A = {}
    for a, l in LAM:
        x = np.log([REF[a][lo][1] for lo in los])
        y = np.log([REF[a][lo][0] for lo in los])
        p = np.polyfit(x, y, 1)
        A[a] = dict(lam=l, nu=float(p[0]),
                    rms=float(np.std(y - np.polyval(p, x))))
        print(f"  {a:9s} {l:5.2f} {p[0]:8.4f} {A[a]['rms']:10.5f}")

    print("\n  (b) λ YÖNÜ (her bantta ayrı ayrı, λ gazları boyunca):")
    print(f"  {'bant lo':>8s} " + " ".join(f"{a:>9s}" for a, _ in LAM) +
          f" {'ν_λ (hepsi)':>12s}")
    B = {}
    for lo in los:
        wx = np.array([REF[a][lo][1] for a, _ in LAM])
        kb = np.array([REF[a][lo][0] for a, _ in LAM])
        p = np.polyfit(np.log(wx), np.log(kb), 1)
        B[lo] = float(p[0])
        print(f"  {lo:8.2f} " + " ".join(f"{v:9.4f}" for v in kb) +
              f" {p[0]:12.4f}")
    print(f"\n  ν_λ ortalaması = {np.mean(list(B.values())):.4f} "
          f"± {np.std(list(B.values())):.4f}")

    print("\n  (b') λ YÖNÜ, KOMŞU ÇİFTLER (yerel eğim; ν=1 ⇒ c sabit):")
    print(f"  {'λ çifti':>14s} " + " ".join(f"{lo:9.2f}" for lo in los) +
          f" {'ort':>8s}")
    C = []
    for i in range(len(LAM) - 1):
        a1, l1 = LAM[i]
        a2, l2 = LAM[i + 1]
        v = [np.log(REF[a1][lo][0] / REF[a2][lo][0])
             / np.log(REF[a1][lo][1] / REF[a2][lo][1]) for lo in los]
        C.append(dict(cift=f"{l1:.2f}→{l2:.2f}", nu=[float(z) for z in v],
                      ort=float(np.mean(v))))
        print(f"  {l1:.2f} → {l2:.2f}   " +
              " ".join(f"{z:9.4f}" for z in v) + f" {np.mean(v):8.4f}")

    print("\n  YORUM: ν=1 olduğu yerde c(λ) DÜZ olur; ν<1 olduğu yerde "
          "c λ ile düşer.")
    json.dump(dict(bant=A, lam_bant=B, cift=C),
              open(SCR / "UYE_USSU.json", "w"), indent=1, default=float)
    print(f"\n-> {SCR}/UYE_USSU.json")


if __name__ == "__main__":
    main()
