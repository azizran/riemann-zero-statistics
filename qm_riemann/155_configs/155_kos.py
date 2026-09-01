"""
155 — KOŞU SÜRÜCÜSÜ
===================
Kullanım:  155_kos.py <veri> <taban> <bantset>

veri     : son | orta | dusuk | keskin | A4
             son    = zeros6 son 300k          (L=12.030)
             orta   = zeros6 Z[850000:1150000] (L=11.464)  [154'ün 'orta'sı]
             dusuk  = zeros6 Z[200000:500000]  (L=10.484)  ← H-L AYIRICISI
             keskin = 152'nin saf asal merdiveni (kesim τ≤1.00, keskin)
             A4     = 152'nin erfc-0.68/0.125 gazı (154'ün önbelleği)
taban    : regresyon tabanı, ör. 0.40 / 0.46 / 0.52
bantset  : ince (0.02 ızgara, τ 0.42–0.56) | kaba (0.03, τ 0.42–0.60)
           | std154 (154'ün 0.03'lük ORTAK ızgarası — doğrulama)

NOT (taban kısıtı): eta_zinciri τ ≤ taban çizgilerini regresyonda
ÇIKARIR. Bu yüzden bir bant ancak lo ≥ taban ise ölçülebilir; sürücü
otomatik olarak alt bantları eler ve kaç bant elendiğini yazar.
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
K = importlib.import_module("155_cekirdek")

HERE = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR154 = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
              "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/154")
SCR152 = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
              "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/152")
OUT = K.SCR

PENCERE = {"son": (-300000, None), "orta": (850000, 1150000),
           "dusuk": (200000, 500000)}
STD154 = [(round(0.465 + 0.03 * i, 3), round(0.495 + 0.03 * i, 3))
          for i in range(11)]


def veri_yukle(ad):
    if ad in PENCERE:
        d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
        Z = np.sort(np.asarray(d["zeros"], dtype=float))
        a, b = PENCERE[ad]
        return Z[len(Z) + a:] if b is None else Z[a:b]
    if ad == "keskin":
        return np.sort(np.load(SCR152 / "z_keskin.npy"))
    if ad == "A4":
        return np.sort(np.load(SCR154 / "z_A4.npy"))
    p = OUT / f"z_{ad}.npy"                 # 155_gaz.py'nin ürettikleri
    if p.exists():
        return np.sort(np.load(p))
    raise SystemExit(f"bilinmeyen veri: {ad}")


if __name__ == "__main__":
    veri, taban, bset = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    # 4. argüman "sahte:<gercek>-<ust>" → SAHTE-TABAN kontrolü:
    #    τ≤gercek GERÇEK asal frekansları + (gercek,ust] aralığında AYNI
    #    SAYIDA ama asal OLMAYAN frekans. Bant elemesi 'gercek' tabanına
    #    göre yapılır (asal çizgiler orada hâlâ η'nın içinde).
    sahte = None
    if len(sys.argv) > 4 and sys.argv[4].startswith("sahte:"):
        a, b = sys.argv[4][6:].split("-")
        sahte = (float(a), float(b))
    t0 = time.time()
    z = veri_yukle(veri)
    bant = {"ince": K.INCE_NEG, "kaba": K.KABA_NEG, "std154": STD154,
            "genis": K.GENIS, "yuksek": K.YUKSEK}[bset]
    elenen = [b for b in bant if b[0] < taban - 1e-9]
    bant = [b for b in bant if b[0] >= taban - 1e-9]
    print(f"=== 155 / veri={veri} / taban={taban} / bant={bset} ===")
    print(f"    n={len(z)}  t∈[{z[0]:.1f},{z[-1]:.1f}]  "
          f"{len(bant)} bant ({bant[0][0]}–{bant[-1][1]}); "
          f"{len(elenen)} bant taban altında elendi", flush=True)
    r = K.olc155(z, f"{veri}-t{taban}-{bset}", bant,
                 anahtar=veri, taban=taban, cap=4000, sahte=sahte)
    r["veri"] = veri
    r["bantset"] = bset
    r["sahte"] = sahte
    r["t_lo"] = float(z[0])
    r["t_hi"] = float(z[-1])
    r["sure_s"] = time.time() - t0
    ekad = "" if sahte is None else f"S{sahte[0]}-{sahte[1]}"
    p = OUT / f"tau0_{veri}{ekad}_t{taban}_{bset}.json"
    p.write_text(json.dumps(r, indent=1))
    print(f"\n-> {p}   ({(time.time()-t0)/60:.1f} dk)", flush=True)
