"""
173c — HAKEM ÖLÇÜMÜ: λ = 0.40 (L040) ve λ = 1.45 (L145)
========================================================
ÖLÇÜM ZİNCİRİ KOPYALANMAZ: `167_olcum.kos` AYNEN çağrılır — 170'in
L060/L115 ve 171'in L050/L130 ölçümleriyle birebir aynı yol
(`kos(<gaz>, 0.40, 0.95, 0, 0)`). Tek yapılan `167_ortak.KUNYE`'ye
künyeyi eklemektir.

ÖN KAYITLAR bu koşulardan ÖNCE yazıldı:
   scratchpad/173/ONKAYIT_L040.json
   scratchpad/173/ONKAYIT_L145.json
ve 172'nin mühürleri (scratchpad/172/G4.json, 4 Eylül 13:32) hiç
değiştirilmedi.

Kullanım: 173c_olcum.py <L040|L145|L140>
Çıktı:    scratchpad/167/C_<gaz>.json  (167'nin kendi dosya düzeni)

SONUÇ (yalnız gerçek koşudan; 3.7-3.8 dk her biri):
 L040: 169'un sağlık filtresi BEŞ HÜKÜM BANDINDA DA GEÇTİ
   (R_bant 1.5289-1.6835, SNR 245-29). Hüküm üyeleri 173d'de:
   c_WX = 0.4098 ± 0.0049, M = 1.4086 ± 0.0232, θ = 1.02923.
 L140: sağlık filtresi BEŞ BANTTA DA GEÇTİ (R_bant 0.9829-1.0370;
   pencerenin ALTINDAKİ iki bant, lo = 0.44/0.48, 0.9682/0.9773 ile
   düşüyor ama onlar hüküm penceresinde değil).
   c_WX = 0.4338 ± 0.0240, M = 0.8634 ± 0.0310, θ = 0.73954.
 L145: **169'un SAĞLIK FİLTRESİ BEŞ BANDIN HEPSİNDE DÜŞTÜ**
   (R_bant = 0.9209 0.9303 0.9373 0.9505 0.9598 < 0.98).
   **FİLTRE GEVŞETİLMEDİ.** λ = 1.45 için c ve M ailenin kendi
   konvansiyonuyla HÜKME BAĞLANAMAZ; sayılar yalnız TEŞHİS olarak
   taşındı (c ≈ 0.4098, M ≈ 0.7663, θ ≈ 0.6426).
   Ayrıca λ = 1.45'te KALİB_u2 τ ile tekdüze inmiyor: 0.2578 0.2354
   0.2104 0.1933 0.1902 0.1909 0.2013 0.2212 0.2441 — τ ≈ 0.62'de
   MİNİMUM yapıp yükseliyor, yani A(τ) = exp(−ατ²) biçimi kırılmış.
"""
import importlib
import sys
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("167_configs", "166_configs", "165_configs", "163_configs",
           "160_configs", "159_configs", "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
ORT = importlib.import_module("167_ortak")
for _ad, _l in (("L040", 0.40), ("L145", 1.45), ("L140", 1.40)):
    ORT.KUNYE[_ad] = dict(gercek=False, lam=_l, tau_ust=1.00, pen=None,
                          aile="lam", T=1.0)
O167 = importlib.import_module("167_olcum")

if __name__ == "__main__":
    ad = sys.argv[1]
    if ad not in ("L040", "L145", "L140"):
        raise SystemExit("kullanım: 173c_olcum.py <L040|L145|L140>")
    O167.kos(ad, 0.40, 0.95, 0, 0)
