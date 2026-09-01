"""
156 — FAZ KAPASİTESİ (sarma-açılmış): kanal TEK BAŞINA hangi fazı
üretebiliyor?

156_cekirdek.olc3 içindeki `kapasite` alanı arg M'nin SARILI (principal
value) menzilini kaydeder; bu, |arg| > π'ye giden eğrilerde yanıltıcıdır
(sarma sıçraması sahte bir "menzil genişliği" üretir). Burada aynı tarama
np.unwrap ile SARMA-AÇILMIŞ olarak yeniden yapılır ve rapor tablosu
bundan üretilir. k3_*.json'daki ham alan koşulduğu gibi bırakıldı;
raporun kapasite tablosu YALNIZ bu dosyadan gelir.

Kullanım: 156_kapasite.py
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
C3 = importlib.import_module("156_cekirdek")
NAK = importlib.import_module("156_nakil")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/156")
TWO_PI = 2 * np.pi
GRID = np.arange(-4.0, 8.0001, 0.05)


def main():
    D = {a: json.load(open(SCR / f"k3_{a}_std.json"))
         for a in ("gercek", "A4")}
    X = {a: NAK.kanallar(a) for a in ("gercek", "A4")}
    kayit = []
    print("gaz     τ̄      φ_Γ      kanal  sarma-açılmış arg menzili        "
          "kök var mı  (R̃ @ menzil ucu, n_eff orada)")
    for a in ("gercek", "A4"):
        s = D[a]["olcek"]
        for b in D[a]["bantlar"]:
            if not b.get("olculdu"):
                continue
            tau = b["tau"]; A = TWO_PI * tau
            ex = np.exp(-1j * A * X[a]["tam"])
            phm = b["phi"]
            for k in ("tam", "lad", "eta"):
                Y = X[a][k] * s[k]
                ar = np.empty(len(GRID)); ne = np.empty(len(GRID))
                for i, R in enumerate(GRID):
                    M, n = NAK.M_of(A, R * Y, ex)
                    ar[i] = np.angle(M); ne[i] = n
                au = np.unwrap(ar)
                # R̃=0 noktasını sarılı değere sabitle (referans)
                i0 = int(np.argmin(np.abs(GRID)))
                au = au - au[i0] + ar[i0]
                lo, hi = float(au.min()), float(au.max())
                # hedef faz (2π katları dahil) menzilde mi?
                kmin = int(np.floor((lo - phm) / TWO_PI))
                kmax = int(np.ceil((hi - phm) / TWO_PI))
                var = any(lo <= phm + kk * TWO_PI <= hi
                          for kk in range(kmin, kmax + 1))
                imx = int(np.argmax(au)); imn = int(np.argmin(au))
                # n_eff kısıtlı menzil: ağırlık çökmeden erişilebilen faz
                m = ne >= 3000
                if m.any():
                    lo_s, hi_s = float(au[m].min()), float(au[m].max())
                    var_s = any(lo_s <= phm + kk * TWO_PI <= hi_s
                                for kk in range(kmin, kmax + 1))
                else:
                    lo_s = hi_s = float("nan"); var_s = False
                print(f"{a:7s} {tau:.4f} {phm:+.4f}  {k:4s}  "
                      f"[{lo:+.3f} … {hi:+.3f}]  span={hi-lo:5.2f}   "
                      f"{'kök VAR' if var else 'KÖK YOK'}   | n_eff≥3000 "
                      f"kısıtlı: [{lo_s:+.3f} … {hi_s:+.3f}] "
                      f"{'kök VAR' if var_s else 'KÖK YOK'}"
                      f"   (max@R̃={GRID[imx]:+.2f} n_eff={ne[imx]:.0f})")
                kayit.append(dict(gaz=a, tau=tau, kanal=k, phi=phm,
                                  arg_lo=lo, arg_hi=hi, kok_var=bool(var),
                                  arg_lo_s=lo_s, arg_hi_s=hi_s,
                                  kok_var_s=bool(var_s),
                                  R_at_max=float(GRID[imx]),
                                  neff_at_max=float(ne[imx]),
                                  R_at_min=float(GRID[imn]),
                                  fit_artik=b["kanal"][k]["artik"]))
    (SCR / "kapasite.json").write_text(json.dumps(kayit, indent=1))
    print(f"\n-> {SCR / 'kapasite.json'}")


if __name__ == "__main__":
    main()
