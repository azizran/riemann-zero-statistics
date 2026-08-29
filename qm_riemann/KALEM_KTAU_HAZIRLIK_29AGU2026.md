# KALEM HAZIRLIK — K(τ) türetimi: hedef eğri ve kısıtlar (29 Ağu 2026)

Sefer dönüşü ortak kalem oturumunun masası. Amaç: çekirdek eğrisi
K(τ)'nin parametresiz türetimi.

## Hedef eğri (138 ölçümü, üç zeros6 penceresi)

Çekirdek K(τ) = κ_arit − κ_ad (τ_p noktalarında):

| p  | erken (L=10.51) | orta (L=11.47) | son (L=12.03) |
|----|-----------------|----------------|----------------|
| 2  | −2.929 @0.066   | −2.867 @0.060  | −2.867 @0.058  |
| 3  | −3.035 @0.105   | −2.957 @0.096  | −2.939 @0.091  |
| 5  | −3.143 @0.154   | −3.034 @0.140  | −3.017 @0.134  |
| 7  | −3.179 @0.186   | −3.061 @0.170  | −3.107 @0.162  |
| 11 | −3.250 @0.229   | −3.196 @0.209  | −3.151 @0.199  |
| 13 | −3.289 @0.245   | −3.245 @0.224  | −3.181 @0.213  |
| 17 | −3.169 @0.270   | −3.265 @0.247  | −3.247 @0.236  |

Adyabatik eğri κ_ad(τ*) (boyalı, üç pencere ort.):
τ*=0.10 → +0.94 ; 0.16 → +1.04 ; 0.24 → +1.35.

## KISAYOL GÖZLEMİ (K4 — türetim stratejisini değiştirir)

Toplam κ(τ) = −2 DÜZ olduğu için çekirdek eğrisi özdeş olarak
**K(τ) = −2 − κ_ad(τ)**. Sağlama: −2−0.94=−2.94 (ölç. −2.94/−2.96 ✓);
−2−1.04=−3.04 (ölç. −3.02..−3.11 ✓); −2−1.35=−3.35 (ölç. −3.25, %3).
⇒ K eğrisinin şekli, düz −2'nin GÖLGESİ. İki bağımsız türetim yolu:
  (Y1) K'yi iz-formülünden doğrudan türet (çekirdek mekaniği);
  (Y2) κ_ad(τ)'yi kinematikten türet (adyabatik yanıt; 109b/110c
       makineleri) + "toplam=−2" ilkesini ayrıca kanıtla → K bedava.
Y2'nin avantajı: κ_ad denge-fiziği, tek gaz gerektirmez; en sert
çapraz test: kinematik κ_ad öngörüsü → K parametresiz.

## Kısıtlar (her türetim adayının vurması gerekenler)

- **K1 (L-değişmezlik):** ölçü dτ/τ (PNT) + yalnız-τ regülatör
  yapısını korumalı (138).
- **K2 (değer/şekil):** ort ≈ −3.05; τ→0.06'da ≈ −2.87; τ=0.25'te
  ≈ −3.25; tam cosπτ-çarpanlı düz DEĞİL (δ(τ) var).
- **K3 (payda sorunu):** 135 oran formunun paydası ΣA²cos2πτ derin
  kuyrukta near-cancel → oran kesime aşırı duyarlı (137 ek taraması:
  σ_u ile −0.6..−2.65 savruluyor). Türetim ya bu paydayı değiştirmeli
  ya fiziksel regülatörle doğal kesmeli. Adaylar:
  (a) yalnız-inkoherent DW (σ≈0.148; tam-σ_u DW değerce RED, 137);
  (b) soğurma (kuyruk gücünün ~2/3'ü örgüce soğuruluyor);
  (c) 135'te atılan τ_p-bağımlı terimlerin geri alınması — 135 sabit
      verir, ölçüm eğri diyor; şekil o terimlerde olabilir.
- **K4 (toplam ilkesi):** K(τ) + κ_ad(τ) = −2 her τ'da (ölçüm ±%3) —
  türetim bu özdeşliği ya varsaymalı ya ÜRETMELİ; üretirse "neden tam
  −2" sorusunun kendisi çözülür.

## Sefer girdisi beklenenler
Kuantum-kaos raporundan (KESIF_SEFERI_KUANTUM_KAOS_29AGU2026.md):
Argaman eylem-korelasyonları, Bogomolny–Keating asal-tarafı,
Sieber–Richter/encounter aletleri — Y1 için hazır alet var mı?
