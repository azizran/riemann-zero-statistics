"""
158 — a'nın KİMLİK SINAVI: KOŞU SÜRÜCÜSÜ
=======================================
Kullanım:  158_kos.py <veri> <taban>

157, ilkel gözlenebilir φ_Γ(τ_eff)'in sıfırdaki EĞİMİNİ ölçtü:
    a = dφ_Γ/dτ|_{τ₀} = 10.759 ± 0.114   (2 pencere × 5 taban, yayılım %3.7)
ve bunun taban konvansiyonundan bağımsız olduğunu gösterdi (taban φ'yi
salt ÖTELİYOR). 157'nin kendi dürüstlük notu: "Sentetik gazlarda a'nın ne
olduğu bu koşuda ölçülmedi." Bu koşu onu ölçüyor.

SORU: a MAKİNE-EVRENSELİ mi (her gazda aynı → tahminci kinematiği),
yoksa KORELASYON-DUYARLI mı (τ₀ gibi merdiveni izliyor → yeni fizik)?

ÖLÇÜM KOPYALANMADI, IMPORT EDİLDİ:
    155_cekirdek.olc155   ← bant döngüsü (çizgi-bazında kayıt + τ_eff)
    155_kos.veri_yukle    ← gaz/pencere yükleyici (z_*.npy önbellekleri)
155_cekirdek zaten 154_cekirdek'in ölçüm zincirini import ediyor
(eta_zinciri, pk_m, _M, _R_ile_faz_es); 155'in kopya-kayması denetimi
(154 ↔ 155, maks fark 2.2e−16) bu koşu için de geçerlidir ve
158_dogrulama.py ile yeniden sınanır.

TEK EKLEME — bant ızgarası:
    IZGARA158 = izgara(0.28, 0.64, 0.02)   → 18 bant, merkezler 0.29…0.63
Bu ızgara 0.42'ye HİZALIDIR (0.42−0.28 = 7×0.02), yani merkezleri
0.43'ten itibaren 157'nin `ince` ızgarasıyla BİREBİR çakışır. Böylece
  * gerçek gazda 157'nin a/b tablosu bant bant yeniden üretilebilir,
  * sentetik gazlarda (τ₀ ≈ 0.489) negatif dal 0.29'a kadar açılır.

Taban altı bantlar elenir (eta_zinciri τ ≤ taban çizgilerini regresyondan
ÇIKARIR, o bantlarda ölçülecek güç kalmaz) — 155_kos'un kuralı aynen.
"""
import importlib
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, "/Users/ugursezen/Desktop/arin/deney/qm_riemann/155_configs")
K = importlib.import_module("155_cekirdek")
KOS155 = importlib.import_module("155_kos")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/158")

# 0.42'ye hizalı ızgara: merkezler 0.29, 0.31, …, 0.63
IZGARA158 = K.izgara(0.28, 0.64, 0.02)


if __name__ == "__main__":
    veri, taban = sys.argv[1], float(sys.argv[2])
    t0 = time.time()
    SCR.mkdir(parents=True, exist_ok=True)
    z = KOS155.veri_yukle(veri)
    bant = [b for b in IZGARA158 if b[0] >= taban - 1e-9]
    elenen = len(IZGARA158) - len(bant)
    print(f"=== 158 / veri={veri} / taban={taban} / ızgara 0.02 ===")
    print(f"    n={len(z)}  t∈[{z[0]:.1f},{z[-1]:.1f}]  "
          f"{len(bant)} bant ({bant[0][0]}–{bant[-1][1]}); "
          f"{elenen} bant taban altında elendi", flush=True)
    r = K.olc155(z, f"{veri}-t{taban}-158", bant,
                 anahtar=veri, taban=taban, cap=4000)
    r["veri"] = veri
    r["bantset"] = "158"
    r["t_lo"] = float(z[0])
    r["t_hi"] = float(z[-1])
    r["sure_s"] = time.time() - t0
    p = SCR / f"F_{veri}_t{taban}.json"
    p.write_text(json.dumps(r, indent=1))
    print(f"\n-> {p}   ({(time.time()-t0)/60:.1f} dk)", flush=True)
