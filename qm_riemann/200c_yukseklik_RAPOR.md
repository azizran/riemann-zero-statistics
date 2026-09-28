# 200-C RAPOR — Merdivenin yükseklik bağımlılığı (LMFDB/Platt, L = 14.2–22.3)

**Ön kayıt:** ONKAYIT_200C.json (sha256 391257a4…), commit 567af03 — ölçümden ÖNCE push'landı.
**Ölçüm:** 27 Eyl 2026, `200c_olcum.py` (16 533 s ≈ 4.6 sa; C4'te b=1 için 10⁶ rastgele sıfır),
`200c_analiz.py` ⇒ `HUKUM_200C.json`.
**Durum:** MÜHÜRLÜ — ortak teftiş onayı 27 Eyl 2026 (kaptan: "mühürleyelim").

## HÜKÜM (önceden kayıtlı kurallar; manşet H-200C-2)

| sınav | sonuç |
|---|---|
| **H-200C-2 (MANŞET, sınıf)** | **KESİN: H_Sγ** — χ²: H_Sγ 10.3 / 12 nokta, H_S1 30.1, H_C 166.2; Δ = 19.8 |
| H-200C-1 (eğim) | **SONLU-YÜKSEKLİK** — γ̂ = 0.749 ± 0.170 (z = 4.40 ≥ 3.5) |
| H-200C-4 (asimptot, 7 pencere) | **TUTAR** — Â = 0.277 ± 0.041 (a_k = 0.2337, 1.0σ), γ̂ = 0.447 ± 0.144 (3.1σ) |

## κ₃ asal izi (kayma = gözlenen − eğik CUE)

| pencere | L | b=0 | b=1 | b=2 | r₁ = s₁/s₀ | r₂ = s₂/s₁ |
|---|---|---|---|---|---|---|
| (A/B havuz) | 9.3–11.7 | +0.283 | +0.150 | +0.086 | 0.53 | 0.57 |
| C1 | 14.19 | +0.281 ± 0.016 | +0.1633 ± 0.0013 | +0.1053 ± 0.0077 | 0.58 | 0.64 |
| C2 | 16.58 | +0.289 ± 0.020 | +0.1744 ± 0.0012 | +0.1183 ± 0.0060 | 0.60 | 0.68 |
| C3 | 18.88 | +0.275 ± 0.017 | +0.1787 ± 0.0015 | +0.1242 ± 0.0069 | 0.65 | 0.70 |
| C4 | 22.31 | +0.283 ± 0.019 | +0.1839 ± 0.0019 | +0.1466 ± 0.0091 | 0.65 | 0.80 |

Kör tahminler (H_Sγ, 200-A/B'den): b=1 +0.159 / +0.164 / +0.169 / +0.174; b=2 +0.105 / +0.114 /
+0.121 / +0.130. Gözlenen b=1 değerleri tahminin ~0.005–0.01 üstünde (yakınsama biraz daha hızlı).

κ₂ kaymaları: b0 −0.082…−0.094, b1 −0.075…−0.079, b2 −0.095…−0.115 — yükseklikten bağımsız,
a_k (−0.088) çevresinde (200-A/B ile aynı tablo).

## Okuma

 - **"Her sıfırda yarıya iner" bir SONLU-YÜKSEKLİK etkisiymiş.** Merdiven yükseldikçe
   düzleşiyor: oranlar 0.53/0.57'den (L≈11) 0.65/0.80'e (L≈22) çıkıyor. Basamaklar arası fark
   L^{−γ} ile (γ ≈ 0.45–0.75) kapanıyor. 198'deki κ_p ≈ p^{−1/3} gibi, 200-A/B'nin "0.55"i de bir
   yükseklik tesadüfü; bu kez hız ve biçim de ölçüldü.
 - **Keating–Snaith/HKO'nun asimptotik resmi ayakta:** üç basamak ortak bir değere yakınsıyor
   (serbest asimptot 0.277 ± 0.041, a_k ile 1σ uyumlu). Açık nokta: b=0 L = 9–22 boyunca ~0.28'de
   DÜZ duruyor (a_k'nın 0.05 üstünde); ortak asimptotun a_k mı (0.234) yoksa ~0.28 mi olduğu bu
   veriyle ayırt edilemiyor.
 - Varyans (κ₂) izi her yükseklikte ve basamakta aynı: a_k ile uyumlu.

## Dürüst notlar

 1. H_Sγ bir kuram değil, 200-A/B'den dışdeğerlemeydi; yeni ve iki kat yüksek veride tuttu.
 2. Faz-kilidi fikri (Landau–Gonek, 1/L ile zayıflar) bu tabloyla UYUMLU: sönmenin yükseklikle
    azalması bekleniyordu. Nicel bağ kurulmadı.
 3. İkincil kontroller temiz: aralık–tepe Pearson'ı CUE MC ile uyumlu (|ρ| ≤ 0.01); M8c ≤ %11.
