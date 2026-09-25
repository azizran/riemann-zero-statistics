# 198 — κ_p'NİN YÜKSEKLİK BAĞIMLILIĞI: paydaki asal sönmesi sonlu-yükseklik etkisi

**Soru (KALEM_KAPPA_YUKSEKLIK_25EYL2026):** 197'de (veri-sonrası) görülen, paydaki (k_p = +1)
asallara bağlı sönme κ_p (≈ p^{−1/3} gibi görünüyordu) yükseklikle BK'ya (κ → 1) mı
yaklaşıyor (H_S, γ = 1), yoksa sabit bir yapı mı (H_C, γ = 0)? Model 1 − κ_p(L) =
B_p (L/L_d)^{−γ}.

**Ön-kayıt:** KALEM + 198_configs (198k0 zincir, 198a/b/c) + ONKAYIT_198.json (sha256 d93cac2c…)
commit fc079c7, ölçümden ÖNCE push. Süreç (kör pencereler görülmeden): Sonnet türetim teftişi
(2 KRİTİK: W_son yarı-körlüğü → nicel gerekçe; W_alt'ta blok-içi L yayılımı 0.026-0.217 → eşit-ΔL
bloklar + M7/M8 kapıları), iki makine raporu (W_alt kısmi kapsama; M8 karar niceliği γ ile).

## Kapılar (hepsi GEÇTİ)

| kapı | sonuç |
|---|---|
| M1 kod yolu | 197 İKİNCİL genlikleri ve s bağıl fark 0 |
| M2 (üç pencere) | ana-6 yanar: alt 15-22σ, düşük 23-60σ, son 29-165σ; ayna 2A₂/3A₃ = 1.04 / 0.95 / 0.94 |
| M5' uyum kalitesi | artık/std 0.04-0.07; −log2/+log2 profil korelasyonu ≥ 0.997; ayna R_2, R_3 ∈ [0.97, 1.04] |
| M7 havuz değişmezliği | χ²₄ = 5.79, p = 0.22 |
| M8 geometri yanlılığı | γ̂ yanlılığı +0.013 (0.054σ); Λ yanlılıkları ≈ −0.01 |
| M6 güç | H-198a doğru 0.930 / 0.935, yanlış 0.01; H-198b güç 0.33 / 0.40 ⇒ KAYIT |
| M4 | σ_eff(γ) = 0.214 ≤ 0.35 |

## Ölçülen κ_p (8-grup jk se)

| pencere | L_W | κ_2 | κ_3 | κ_5 | κ_7 |
|---|---|---|---|---|---|
| W_alt (KÖR) | 9.343 | 0.784 ± 0.014 | 0.663 ± 0.013 | 0.583 ± 0.019 | 0.521 ± 0.028 |
| W_düşük (ref.) | 10.484 | 0.789 ± 0.011 | 0.713 ± 0.015 | 0.604 ± 0.021 | 0.538 ± 0.020 |
| W_son (yarı-kör) | 12.030 | 0.837 ± 0.009 | 0.770 ± 0.007 | 0.712 ± 0.014 | 0.616 ± 0.026 |
| öngörü H_S, W_son | | 0.821 | 0.756 | 0.659 | 0.591 |
| öngörü H_C | | 0.795 | 0.720 | 0.608 | 0.530 |

## HÜKÜM (eşikler donmuş; kurtarma yok)

| hipotez | ölçülen | hüküm |
|---|---|---|
| **H-198a** (γ, üç pencere) | **γ̂ = 1.36 ± 0.21** (σ_eff); \|γ̂ − 1\| = 1.7σ ≤ 2σ, \|γ̂\| = 6.4σ | **H_S TUTAR** — sabit yapı (H_C) 6.4σ dışlandı |
| H-198b (yalnız kör W_alt) | Λ_alt = 1.074 ± 0.055 (H_S 1.122, H_C 1.000) | KAYIT (güç < 0.80); hesaplanan: belirsiz |

## KAYIT

- **Çapasız oran** (desenin fark edildiği düşük pencere HİÇ kullanılmadan):
  Λ_son/Λ_alt = 0.707 ± 0.037 (H_S 0.777, H_C 1.000) ⇒ H_C 7.9σ dışı; H_S'ten 1.9σ dik.
  Seçim etkisi (ortalamaya dönüş) bu sonucu açıklamaz.
- Λ_son = 0.779 ± 0.046 (H_S 0.872, H_C 1).
- Diğer pozitif-üs sönmeleri de L ile BK'ya gevşiyor: buçuklu aile R(5/2) 0.79 → 0.83 → 0.92,
  R(7/2) 0.75 → 0.81 → 0.88; çeyrek aile q_f/s 0.65 → 0.69 → 0.75. Ayna yanı tam sayıları (negatif
  üsler) her pencerede 1 (R_2, R_3 ∈ [0.97, 1.04]).
- Bileşik çarpımsallık L ile iyileşiyor: +log6 κ_r/(κ_2κ_3) 1.20 → 1.15 → 1.10; +log10 1.16 →
  1.22 → 1.07. +log(3/2), +log(5/2) her pencerede κ-çarpımıyla uyumlu (0.98-1.07); +log(5/4),
  +log(5/3) 0.72-0.97 (tutarlı biçimde düşük).
- μ=0 konumları her üç pencerede küçük negatif, L'den bağımsız: +log4 −0.023 / −0.023 / −0.032 s,
  +log8 −0.021 / −0.018 / −0.017 s, +log9 −0.020 / −0.016 / −0.020 s (BK: 0).
- M8 düzeltmeli κ (model yanlılığı çıkarılmış): alt 0.772 / 0.653 / 0.600 / 0.536; düşük 0.778 /
  0.707 / 0.626 / 0.553; son 0.828 / 0.770 / 0.723 / 0.633 — aynı resim.

## Dürüstlük notları

1. Kararı büyük ölçüde YARI-KÖR W_son taşıyor; TAMAMEN KÖR W_alt tek başına belirsiz (M6'nın
   önceden söylediği gibi). Çapasız oran da W_son'a dayanır.
2. Biçim temiz bir 1/L değil: alt → düşük arasında değişim küçük (κ_2 0.784 → 0.789), düşük → son
   arasında büyük (κ_5 0.604 → 0.712). γ̂ = 1.36 bu eğriliğin ortalaması; üç noktayla fonksiyon
   biçimi (1/L, 1/L², log p/L…) ayrıştırılamaz.
3. 197'nin "κ_p ≈ p^{−1/3}" okuması bir YÜKSEKLİK TESADÜFÜYDÜ (L ≈ 10.5'te); sabit bir yasa değil.

## MANŞET (aday)

> 197'de paydaki asallara bağlı görülen BK-altı sönme (κ_p) sabit bir yapı değil, sonlu-yükseklik
> etkisi: üç pencerede (L = 9.3, 10.5, 12.0) κ_p BK'ya doğru gevşiyor, 1 − κ_p ∝ L^{−γ} ile
> γ̂ = 1.36 ± 0.21; sabit-yapı hipotezi 6.4σ (çapasız oranla 7.9σ) dışlandı. Ayna yanı (negatif
> üsler) her yükseklikte tam BK. Tablo, BK sanısının pozitif-üs katsayılarına asimptotik olarak,
> ölçülebilir bir 1/L-mertebesi düzeltmeyle yaklaşıldığını söylüyor.

## Açık kalemler

- Düzeltmenin biçimi: daha fazla pencere (L ∈ [8, 12.5], zeros6 içinde) ile 1/L mi, log p/L mi,
  1/L² mi? Ya da BK'nın bilinen alt-mertebe terimlerinden (Bogomolny 2007 / Conrey–Snaith) κ_p(L)
  türetimi — ön-kayıtlı öngörüyle.
- μ=0 negatif artıkları (s'nin ~%2'si, L'den bağımsız): BK-ötesi (ikinci mertebe) terim adayı.

## Teslim

- 198_configs/: 198k0_zincir.py, 198a_onkayit.py, 198b_harita.py, 198c_hukum.py, ONKAYIT_198.json,
  F_198.json, M1/M6/M7/M8 json'ları, HUKUM_198.json
- 198_kappa_yukseklik.png
- Veri: scratchpad/198 (gitignore)
