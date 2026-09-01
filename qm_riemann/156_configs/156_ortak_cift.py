"""
156 — ORTAK ÇİFT ARAMASI: tek bir (R_l, R_e) HER İKİ gazın fazını da
açıklıyor mu?

NİÇİN AYRI DOSYA: 156_nakil.py'nin ilk sürümü iki gazın faz-eşlenmiş
sırtlarını R_e(R_l) eğrileri gibi doğrusal ara değerle kesiştiriyordu.
Sırt ÇOK DALLI çıktı (kök bulucu "|R_e|'si en küçük kök" kuralıyla
dallar arasında sıçrıyor: A4'te τ=0.66'da R_l≈0.5'te R_e +8.4'ten
−3.5'e atlıyor). Doğrusal ara değer bu sıçramanın üstünden köprü kurup
OLMAYAN bir kesişim uyduruyordu (nakil çıktısında A4 artığı 1.48 rad
çıkması bunun işaretidir). Burada kesişim, ara değersiz, doğrudan
2-B ızgarada aranır:

  J(R_l,R_e) = |Δφ_gerçek| + |Δφ_A4|,  Δφ = sarılı (arg M − φ_Γ)
  önce kaba ızgara, sonra en iyi hücrenin çevresinde ince ızgara.
  min J < 0.04 rad ise ORTAK ÇÖZÜM VAR; değilse yok.

Kabul edilen noktada ayrıca her iki gazın |ΔRe|'si ve n_eff'i basılır —
faz kesişimi bulunsa bile aşırı-belirleme ve ağırlık çöküşü sınavını
geçmesi gerekir.

Kullanım: 156_ortak_cift.py
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
NAK = importlib.import_module("156_nakil")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/156")
TWO_PI = 2 * np.pi
RL = np.linspace(-1.0, 4.0, 81)
RE = np.linspace(-8.0, 14.0, 177)

# İKİNCİ GEÇİŞ (bu dosyanın ilk koşusu yalnız GLOBAL en iyiyi arıyordu;
# çıktısı log_ortak_1gecis.txt'de duruyor). Eklenen: J'nin BÜTÜN yerel
# minimumları taranır. Gerekçe — iki bilinmeyen (R_l,R_e) ve iki kısıt
# (iki gazın fazı) olduğu için ortak bir çözümün VAR OLMASI sayım gereği,
# kanıt değildir; ayrıca çözüm TEK olmayabilir. Asıl sınav, ortak çözümün
# hiçbir ek serbestlik olmadan İKİ gazın Re'sini de tutup tutmadığıdır.
# Bu yüzden tüm çözümler listelenir ve aralarında Re'yi en iyi tutan
# raporlanır.


def sarili(x):
    return (x + np.pi) % TWO_PI - np.pi


def tara(A, X, ex, phi, rl, re):
    F = np.empty((len(rl), len(re)))
    for i, a in enumerate(rl):
        base = a * X["lad"]
        for j, b in enumerate(re):
            M, _ = NAK.M_of(A, base + b * X["eta"], ex)
            F[i, j] = sarili(np.angle(M) - phi)
    return F


def main():
    D = {a: json.load(open(SCR / f"k3_{a}_std.json"))
         for a in ("gercek", "A4")}
    X = {a: NAK.kanallar(a) for a in ("gercek", "A4")}
    kayit = []
    bg = {b["tau"]: b for b in D["gercek"]["bantlar"] if b.get("olculdu")}
    bs = {b["tau"]: b for b in D["A4"]["bantlar"] if b.get("olculdu")}
    for tau in sorted(set(bg) & set(bs)):
        A = TWO_PI * tau
        ex = {a: np.exp(-1j * A * X[a]["tam"]) for a in X}
        Fg = tara(A, X["gercek"], ex["gercek"], bg[tau]["phi"], RL, RE)
        Fs = tara(A, X["A4"], ex["A4"], bs[tau]["phi"], RL, RE)
        J = np.abs(Fg) + np.abs(Fs)
        dl = RL[1] - RL[0]; de = RE[1] - RE[0]
        # tüm yerel minimumlar (8-komşuluk), J küçük olanlar
        aday = []
        for i in range(1, len(RL) - 1):
            for j in range(1, len(RE) - 1):
                blok = J[i - 1:i + 2, j - 1:j + 2]
                if J[i, j] <= blok.min() and J[i, j] < 0.25:
                    aday.append((i, j))
        pat = ("  [bant patlak — hüküm yok]"
               if (bs[tau]["patlak"] or bg[tau]["patlak"]) else "")
        print(f"=== τ={tau:.4f}   φ_gerçek={bg[tau]['phi']:+.4f}  "
              f"φ_A4={bs[tau]['phi']:+.4f}   kaba aday {len(aday)}{pat}")
        cozum = []
        for i, j in aday:
            rl2 = np.linspace(RL[i] - dl, RL[i] + dl, 25)
            re2 = np.linspace(RE[j] - de, RE[j] + de, 25)
            Fg2 = tara(A, X["gercek"], ex["gercek"], bg[tau]["phi"], rl2, re2)
            Fs2 = tara(A, X["A4"], ex["A4"], bs[tau]["phi"], rl2, re2)
            J2 = np.abs(Fg2) + np.abs(Fs2)
            i2, j2 = np.unravel_index(int(np.argmin(J2)), J2.shape)
            if J2[i2, j2] >= 0.04:
                continue
            Rl, Re = float(rl2[i2]), float(re2[j2])
            c = dict(Rl=Rl, Re=Re, J=float(J2[i2, j2]))
            for a, b in (("gercek", bg[tau]), ("A4", bs[tau])):
                M, ne = NAK.M_of(A, Rl * X[a]["lad"] + Re * X[a]["eta"], ex[a])
                c[a] = dict(dphi=float(abs(sarili(np.angle(M) - b["phi"]))),
                            dRe=float(abs(M.real - b["Gre"])), neff=ne,
                            ReM=float(M.real))
            c["dRe_max"] = max(c["gercek"]["dRe"], c["A4"]["dRe"])
            c["neff_min"] = min(c["gercek"]["neff"], c["A4"]["neff"])
            if any(abs(Rl - q["Rl"]) < 1e-3 and abs(Re - q["Re"]) < 1e-3
                   for q in cozum):
                continue
            cozum.append(c)
        if not cozum:
            print("    ORTAK FAZ ÇÖZÜMÜ YOK (J<0.04 rad hiçbir yerde)\n")
            kayit.append(dict(tau=tau, cozum=[], patlak=bool(pat)))
            continue
        cozum.sort(key=lambda c: c["dRe_max"])
        for c in cozum:
            print(f"    (R_l,R_e)=({c['Rl']:+.4f},{c['Re']:+.4f})  J="
                  f"{c['J']:.4f}   gerçek: ReM={c['gercek']['ReM']:+.4f} "
                  f"|ΔRe|={c['gercek']['dRe']:.4f} n_eff="
                  f"{c['gercek']['neff']:.0f}   A4: ReM="
                  f"{c['A4']['ReM']:+.4f} |ΔRe|={c['A4']['dRe']:.4f} n_eff="
                  f"{c['A4']['neff']:.0f}   maks|ΔRe|={c['dRe_max']:.4f}")
        en = cozum[0]
        hk = ("ORTAK ÇİFT AŞIRI-BELİRLEMEYİ GEÇTİ" if en["dRe_max"] < 0.06
              else "ortak çift Re'yi TUTMUYOR")
        print(f"    hüküm: {len(cozum)} ortak faz çözümü; en iyisinde "
              f"maks|ΔRe|={en['dRe_max']:.4f} → {hk}\n")
        kayit.append(dict(tau=tau, cozum=cozum, en_iyi=en, hukum=hk,
                          patlak=bool(pat)))
    (SCR / "ortak_cift.json").write_text(json.dumps(kayit, indent=1))
    print(f"-> {SCR / 'ortak_cift.json'}")


if __name__ == "__main__":
    main()
