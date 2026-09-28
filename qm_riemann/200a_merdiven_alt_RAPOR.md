# 200-A RAPOR — Eğim merdiveninin alt basamakları (b=0 rastgele t, b=1 sıfırlarda Z′)

**Ön kayıt:** ONKAYIT_200A.json (sha256 f0f0ce57…), commit aa861b5 — ölçümden ÖNCE push'landı.
**Ölçüm:** 26 Eyl 2026, `200a_olcum.py` (35 sn; 3×10⁶ rastgele t + 1 980 184 sıfırda Z′),
`200a_analiz.py` ⇒ `HUKUM_200A.json`. Diske yalnız blok güç toplamları yazıldı (M9).
**Durum:** MÜHÜRLÜ — ortak teftiş onayı 26 Eyl 2026 (kaptan: "200-A'yı mühürleyelim").

## HÜKÜM (önceden kayıtlı kurallarla)

| sınav | sonuç |
|---|---|
| **H-200A-1** (alt basamaklarda aritmetik imza) | **TUTAR** — Δ₁ = χ²_CUE − χ²_{a_k} = **287.1** (eşik 25) |
| Gauss rakibi | χ²_Gauss − χ²_{a_k} = 3661.6 — reddedildi |
| **H-200A-2** (model seçimi) | **M* = a_k, KESİN** — ikinci en iyi M_CUE, Δχ² = 287.1 (eşik 9); S* = {a_k} |
| **Uyum kalitesi G** | χ²_min/12 = 28.2/12 = **2.35 ⇒ İYİ** (teori payı dahil) |
| M8 (blok sayısı) | 64→128'de en büyük SE değişimi %27 > %20 ⇒ kural gereği 128 blok |
| İkincil: basamak başına seçim | b=0: a_k (χ² 6.1); b=1: a_k (χ² 22.1) ⇒ **iki basamak aynı modeli seçti** |
| İkincil: tam korelasyonlu σ_teori | M* = a_k (değişmedi) |
| İkincil: sürekli-X hibrit | X* = 2.99, χ² = 129.3 (a_k'dan 101 kötü) ⇒ hibrit aile toptan kötü |
| 200-B toleransı (kural) | κ₂: 0.0238, κ₃: 0.1422 |

χ² (12 bileşen, kovaryanslı + σ_teori): CUE 315.2 · **a_k 28.2** · hyb2 1243.0 · hyb3 1133.2 ·
hyb5 872.5 · hyb7 991.1 · hyb11 959.4 · Gauss 3689.7.

## Sayılar

Gözlenen k-istatistikleri (± f-kalibre jackknife SE, 128 blok) ve tahminler:

| | gözlenen | M_CUE | M_ak | a_k artığı |
|---|---|---|---|---|
| W1 b0 κ₂ | 1.8191 ± 0.0045 | 1.905 | 1.817 | +0.002 |
| W1 b0 κ₃ | −2.0887 ± 0.0181 | −2.387 | −2.153 | +0.065 |
| W1 b1 κ₂ | 0.4672 ± 0.0022 | 0.540 | 0.452 | +0.015 |
| W1 b1 κ₃ | −0.0465 ± 0.0027 | −0.190 | +0.043 | −0.090 |
| W2 b0 κ₂ | 1.8899 ± 0.0044 | 1.969 | 1.880 | +0.010 |
| W2 b0 κ₃ | −2.1237 ± 0.0195 | −2.397 | −2.163 | +0.039 |
| W2 b1 κ₂ | 0.5200 ± 0.0014 | 0.592 | 0.503 | +0.017 |
| W2 b1 κ₃ | −0.0519 ± 0.0021 | −0.198 | +0.036 | −0.088 |
| W3 b0 κ₂ | 1.9348 ± 0.0039 | 2.017 | 1.929 | +0.006 |
| W3 b0 κ₃ | −2.1240 ± 0.0130 | −2.403 | −2.170 | +0.046 |
| W3 b1 κ₂ | 0.5586 ± 0.0011 | 0.632 | 0.544 | +0.015 |
| W3 b1 κ₃ | −0.0509 ± 0.0010 | −0.203 | +0.031 | −0.082 |

Tanı: Var_jk/Var_iid = 1.2–8.9 (b=0 κ₃'te en büyük: rastgele t'ler açılmış birimde yoğun
örneklendi ve sol kuyruk ağır); beklendiği gibi ≫ 1.

## Aritmetik kayma = gözlenen − M_CUE (veri-sonrası, KAYIT)

| basamak | κ₂ kayması | κ₃ kayması |
|---|---|---|
| b=0 (rastgele t) | −0.0820 ± 0.0024 (W: −0.086, −0.079, −0.082) | **+0.283 ± 0.009** (+0.298, +0.273, +0.279) |
| b=1 (sıfırlarda Z′) | −0.0726 ± 0.0008 (−0.073, −0.072, −0.073) | **+0.150 ± 0.001** (+0.144, +0.146, +0.152) |
| a_k öngörüsü (her basamak) | −0.0881 | +0.2337 |

Gözlemler:
1. **κ₂ kayması neredeyse basamaktan bağımsız** ve a_k sabitine yakın (−0.082 / −0.073 vs
   −0.088). Keating–Snaith normalizasyonu (N = L) L ≈ 9–12'de bile κ₂'yi ~%10–20 içinde
   veriyor. Hibrit modeller (küçük X) κ₂'yi çok büyük öngördüğü için eleniyor.
2. **κ₃ kayması basamağa bağlı:** b=0'da +0.283, b=1'de +0.150 (a_k: +0.234, ikisinin arası).
   Pencereler arasında sabit (W1→W3 eğilim yok) ⇒ sonlu-L gürültüsü değil, basamağın özelliği.
   Yani "aritmetik katkı her basamakta aynı" merdiven sanısı **κ₃ düzeyinde tam tutmuyor**:
   koşullanan her sıfırla sağa çarpıtan asal katkısı azalıyor.
3. Aday açıklama (spekülatif, sınanmadı): sıfıra koşullamak asal fazlarını Landau–Gonek
   yanlılığıyla kaydırıyor; teftişin kestirimi asal başına Δκ₃ ≈ −0.002…−0.005 idi — gözlenen
   fark (−0.13) bundan çok büyük. Açıklama açık.

## 200-B için sonuç (karar kaptanla)

 - **Önceden kayıtlı taşıma:** M* = a_k ⇒ b=2 tahmini κ₂ = 0.188 / 0.230 / 0.263 (± 0.024),
   κ₃ = +0.168 / +0.162 / +0.158 (± 0.142); D₃ = +0.234.
 - **A'nın açtığı soru:** κ₃ kayması b ile düşüyorsa yakın çiftte (b=2) ne olur?
   - merdiven sabit (a_k): D₃ = +0.234,
   - geometrik azalma (oran 0.150/0.283 = 0.53): D₃ ≈ +0.080,
   - doğrusal azalma (fark −0.133): D₃ ≈ +0.017,
   - salt RMT: 0; Gauss: ≈ −0.03.
   İşaret dönmesi tahmini artık kesin değil: doğrusal azalma sürerse yakın çiftte asal imzası
   neredeyse sıfırlanır. Bu rakipler B verisi görülmeden, A'dan bilgilenmiş hipotezler olarak
   200-B KALEM'ine girebilir.
 - Güç: b=2'de δ̃ < 0.2 için ~14 000 olay; SE(κ₃) ≈ 0.01–0.02 (dar formülün 2.7 katı) ⇒ a_k,
   geometrik ve doğrusal ayrışır (fark ≥ 0.06); doğrusal ile salt RMT (fark 0.017) ayrışmaz.
