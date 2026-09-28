# 200-B RAPOR — Yakın çift basamağı (b=2): asalların çarpıklık izi

**Ön kayıt:** ONKAYIT_200B.json (sha256 b37dc244…), commit 85f0489 — ölçümden ÖNCE push'landı.
**Ölçüm:** 26 Eyl 2026, `200b_olcum.py --n-grid 129` (30 sn; 47 074 olay, δ̃ < 0.3),
`200b_analiz.py` ⇒ `HUKUM_200B.json`, olaylar `200b_olaylar.npz`.
**Durum:** MÜHÜRLÜ — ortak teftiş onayı 27 Eyl 2026 (kaptan: "onaylıyorum, 200-B'yi mühürleyelim").

## HÜKÜM (önceden kayıtlı kurallarla, birincil ε = 0.2, 14 168 olay)

| sınav | sonuç |
|---|---|
| **H-200B-1** (b=2'de aritmetik imza) | **VAR** — D̄₃ = **+0.0859 ± 0.0030** (29σ) |
| **H-200B-1b** (manşet: çarpıklığın işareti) | **DÖNDÜ** — düzeltilmiş havuz κ₃ = +0.0143 ± 0.0030 (4.8σ) |
| **H-200B-2** (A'nın taşıdığı a_k) | **KISMİ** — κ₂: D̄₂ = −0.0988 (a_k −0.0881; fark −0.011, tol 0.024) TUTAR; κ₃: fark −0.148 (tol 0.142) TUTMAZ |
| **H-200B-3** (merdivenin biçimi) | **KESİN: GEOMETRİK** — χ²: GEOMETRİK 0.36, SÖNMÜŞ 24.0 (lin 24.0, RMT 67.8), SABİT 200.7; Δ = 23.6 |

Öncelik kuralı gereği manşet H-200B-3: **GEOMETRİK**. Analiz betiğinin çelişki notu
("D̄₃ bandında") yanlış ifadeli: D̄₃ a_k'nın κ₃ bandının DIŞINDA (0.148 > 0.142); KISMİ
sonucunu κ₂ getiriyor. Yani κ₃ düzeyinde iki kural çelişmiyor: a_k κ₃'te düşüyor, biçim
geometrik.

## Asıl bulgu: asal izi her koşullanan sıfırla ~0.55 katına iniyor

| basamak | ne ölçüldü | κ₃ kayması (asal izi) | bir öncekine oran |
|---|---|---|---|
| b=0 | rastgele t, log\|Z\| | +0.2828 ± 0.0093 | — |
| b=1 | sıfırlarda log\|Z′\| | +0.1499 ± 0.0009 | 0.530 ± 0.018 |
| b=2 | yakın çiftte log(M/δ̃²) | **+0.0859 ± 0.0030** | 0.573 ± 0.020 |

 - **Kör tahmin:** 200-A'nın iki noktasından geometrik uzantı, b=2 verisi görülmeden
   **+0.0794 ± 0.0028** öngörmüştü; ölçülen +0.0859 ± 0.0030 (fark 0.0065, birleşik ~1.6σ).
   Doğrusal uzantı (+0.017) ve salt RMT (0) kesin reddedildi; a_k'nın sabit merdiveni
   (+0.234) χ² = 201 ile reddedildi.
 - İki oran (0.530, 0.573) birbiriyle ~1.6σ uyumlu.
 - **Çarpıklığın işareti döndü ama kıl payı:** yakın çiftte log-tepe çarpıklığı ≈ +0.15
   (κ₃/κ₂^{3/2}); iskeletin a_k tahmini +1.7 idi. Manşet tahmin işaret düzeyinde tuttu,
   büyüklükte tutmadı.

## κ₂ (varyans) kayması

| basamak | κ₂ kayması |
|---|---|
| b=0 | −0.0820 ± 0.0024 |
| b=1 | −0.0726 ± 0.0008 |
| b=2 | −0.0988 ± 0.0036 |
| a_k öngörüsü | −0.0881 |

κ₂ kayması sönmüyor; üç basamakta −0.07…−0.10, a_k sabitinin çevresinde. Geometrik/doğrusal
κ₂ uzantıları (−0.064) b=2'de ~9σ dışarıda. **Tablo:** varyanstaki asal izi basamaktan
bağımsız (a_k), çarpıklıktaki asal izi geometrik sönüyor.

## Sağlamlık (ikincil, KAYIT)

 - ε = 0.1: D̄₃ = +0.0867 ± 0.0066; ε = 0.3: +0.0828 ± 0.0044 ⇒ ε'den bağımsız.
 - Pencere başına D₃ (ε = 0.2): W1 +0.076 ± 0.006, W2 +0.084 ± 0.005, W3 +0.093 ± 0.005 —
   hafif yükselen eğilim (~2.5σ; b=0/b=1'de eğilim yoktu). Kayıt.
 - Bağımsızlık (Teorem 1(b)): Pearson(δ̃, x) gerçek +0.021 / −0.004 / −0.005, CUE MC
   +0.005…+0.008 ⇒ uyumlu (ε = 0.1'de W2 −0.106, n = 356, ~2σ).
 - Aralık yasası: δ̃ < ε/2 oranı 0.124–0.129 (U^{1/3} ⇒ 0.125; CUE 0.128) ✔.
 - M8b: ε = 0.2'de tetiklenmedi.

## Dürüst notlar

 1. Geometrik sönme bir KURAM değil: A'nın iki noktasından türetilmiş bir dışdeğerlemeydi.
    Üçüncü nokta onu kör sınavda doğruladı; ama "neden 0.55?" sorusunun cevabı yok.
 2. Hepsi L ≈ 9–12 aralığında; yükseklikle değişip değişmediği açık (b=2'de hafif eğilim).
 3. Sonlu-ε düzeltmesinin zeta'ya aktarımı varsayım; belirsizliği hata payına %100 eklendi.

## Kapanış cümlesi (teftiş için taslak)

Sıfırlara koşulladıkça asalların varyanstaki izi aynı kalıyor, çarpıklıktaki izi her
koşullanan sıfırla kabaca yarıya iniyor (+0.28 → +0.15 → +0.086). Üçüncü basamağın değeri
ikinci basamaktan önceden, kör olarak öngörüldü ve tuttu.
