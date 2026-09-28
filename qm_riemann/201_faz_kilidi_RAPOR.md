# 201 RAPOR — Faz kilidi: sıfırlar asal dalgalarını hangi fazda yakalıyor?

**Ön kayıt:** ONKAYIT_201.json, commit 16e7c3b (ölçümden önce). **Ölçüm:** 28 Eyl, 6 saniye (yalnız
sıfır konumları; 7 pencere, L 9.3–22.3). **Durum:** MÜHÜRLÜ — ortak teftiş onayı 28 Eyl 2026 (kaptan: "mühürleyelim").

## HÜKÜM

| sınav | sonuç |
|---|---|
| **H-201-1** (tek sıfırlar, Landau–Gonek) | **TUTAR** — R₁ = 1.000 (pencereler: 1.003 ± 0.007, 1.001, 1.001, 1.000, 1.000, 1.000, 1.000); m=2 de 1.000; 1/L ölçeklemesi TUTAR |
| **H-201-2** (yakın çiftlerin orta noktası, δ̃ < 0.2) | **KESİN: 4** — R₂ = 3.874 ± 0.023 |
| b=0 denetimi (rastgele t) | R ≈ 0 (−0.010…+0.005) ✓ |

Ham örnek (C4, L = 22.31, 1.5·10⁶ sıfır): E cos(γ log 2) = −0.021975, Landau–Gonek −0.021973;
E sin ≈ 0 (±4·10⁻⁶). Asal başına çift oranları: W2 3.55–4.08, C4 3.54–4.05.

## Okuma

 1. **Tek sıfırlar:** Landau'nun formülü (1911) dört ondalığa kadar tutuyor — bu bir teoremin
    doğrulaması, keşif değil. Ama yan bulgu çarpıcı: blok başına dalgalanma, bağımsız noktalar için
    beklenenden 7–19 kat küçük (sıfırların "spektral sertliği": doğrusal istatistiklerin varyansı
    rastgele noktalardaki gibi n ile değil, çok yavaş büyür).
 2. **Yakın çiftler (YENİ):** orta noktada asal fazı kilidi tek sıfırınkinin **~4 katı**. Bu,
    teftiş sonrası veriden ÖNCE yazılan "yerel GUE itmesi yerel yoğunlukla ölçeklenir" öngörüsü:
    küçük aralıklı çift yoğunluğu ∝ n(t)⁴ (n(t) asal dalgalarıyla modüle olan yerel sıfır
    yoğunluğu). "Yoğunluklar çarpılır, itme duyarsız" (2×) ve "tek sıfır gibi" (1×) elendi.
 3. **ε bağımlılığı kuramla uyumlu:** R₂ = 4.07 (ε 0.1) → 3.87 (0.2) → 3.74 (0.3) — aralık büyüdükçe
    n⁴δ² yaklaşımı zayıflıyor ve oran 2'ye doğru iniyor. m = 2 harmoniğinde R₂ ≈ 3.05 (ikinci
    mertebe n⁴ açılımı ve M6'daki doğrusal olmama ile uyumlu: sentetik k=4'te 2.90).
 4. **Kaptanın sezgisi ("tepede mi çukurda mı?"):** sıfırlar asal dalgasının tepesinden kaçıyor;
    birbirine çok yakın iki sıfır bu kaçışı dört kat güçlü yapıyor.

## Dürüst notlar

 - "Kilit asal izinin sönmesini açıklıyor mu?" sorusu bu kalemde karara bağlanmadı (bağımsız-Euler
   çevirisi kararsız; teftiş). Ölçülen kilit kuramın öngördüğü büyüklükte — yani 200-A/B/C'deki
   κ₃ sönmesini açıklayacak "fazladan" bir kilit YOK; sönmenin mekanizması başka bir yerde.
 - 4× sonucu literatürde var mı (yerel çift yoğunluğunun asal fazlarına bağlılığı;
   Bogomolny–Keating'in aritmetik düzeltmeleri) — kontrol edilmeli.
