"""
156 — koşucu: üç-kanallı ağırlık ayrışımı, gerçek ve sentetik-A4.

Kullanım:  156_kos.py <gercek|A4|keskin> [std]
  gercek : 128_odl_zeros6_2e6_zeros.npz, son 300k sıfır (150/154 penceresi)
  A4     : 152'nin erfc-0.68/0.125 gazı — 154'ün önbelleğe aldığı Newton
           çözümü (z_A4.npy) OKUNUR, yeniden çözülmez.
  keskin : 152'nin keskin kontrolü (z_keskin.npy) — ek referans.

KOPYA DENETİMİ: koşu sonunda 154'ün aynı pencere için kaydettiği
|Γ|, φ, Re Γ ve R değerleriyle karşılaştırma basılır. Bu, 156'nın Γ
döngüsünün 154'ten sapmadığını gösterir (fark 4. basamakta olmalı).
"""
import json
import sys
import time
from pathlib import Path

import importlib
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
C3 = importlib.import_module("156_cekirdek")

HERE = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
OUT = SCR / "156"
REF = {"gercek": SCR / "154/R_gercek_son_std.json",
       "A4": SCR / "154/R_sentetik_A4_std.json",
       "keskin": SCR / "154/R_sentetik_keskin_std.json"}


def yukle(ad):
    if ad == "gercek":
        d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
        Z = np.sort(np.asarray(d["zeros"], dtype=float))
        return Z[len(Z) - 300000:]
    if ad == "A4":
        return np.sort(np.load(SCR / "154/z_A4.npy"))
    if ad == "keskin":
        return np.sort(np.load(SCR / "152/z_keskin.npy"))
    raise SystemExit(f"bilinmeyen: {ad}")


def main(ad, bset="std", taban=0.52, cap=720):
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    z = yukle(ad)
    print(f"=== 156 / {ad} / bant {bset} / taban {taban} cap {cap} — "
          f"{len(z)} nokta ===", flush=True)
    bant = C3.BANTLAR_STD if bset == "std" else C3.GOREV_BANTLARI
    r = C3.olc3(z, f"156-{ad}-t{taban}", bantlar=bant, taban=taban, cap=cap)
    r["ad"] = ad; r["bantset"] = bset; r["taban"] = taban; r["cap"] = cap
    r["sure_s"] = time.time() - t0

    # ---- 154 ile kopya denetimi (yalnız standart 0.52 tabanında anlamlı) ----
    p = REF.get(ad) if abs(taban - 0.52) < 1e-9 else None
    if p and p.exists():
        ref = {b["tau"]: b for b in json.load(open(p))["bantlar"]
               if b.get("olculdu")}
        print("\nKOPYA DENETİMİ (154 vs 156)  τ    Δ|Γ|      Δφ        ΔReΓ"
              "      R_154   R̃_tam(156)", flush=True)
        for b in r["bantlar"]:
            if not b.get("olculdu") or b["tau"] not in ref:
                continue
            q = ref[b["tau"]]
            print(f"   {b['tau']:.4f}  {abs(b['absG']-q['absG']):.2e}  "
                  f"{abs(b['phi']-q['phi']):.2e}  "
                  f"{abs(b['Gre']-q['Gre']):.2e}   {q['R']:+7.3f}  "
                  f"{b['kanal']['tam']['Rn']:+7.3f}", flush=True)

    ek = "" if abs(taban - 0.52) < 1e-9 else f"_t{taban}"
    (OUT / f"k3_{ad}_{bset}{ek}.json").write_text(json.dumps(r, indent=1))
    print(f"\nbitti — {(time.time()-t0)/60:.1f} dk  -> "
          f"{OUT / f'k3_{ad}_{bset}{ek}.json'}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1],
         sys.argv[2] if len(sys.argv) > 2 else "std",
         float(sys.argv[3]) if len(sys.argv) > 3 else 0.52,
         int(sys.argv[4]) if len(sys.argv) > 4 else 720)
