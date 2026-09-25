# 197 — AYNA YASASI: Bogomolny–Keating'in ayna-yanı uyduları kör bantta

**Soru (KALEM_AYNA_YASASI_25EYL2026):** BK çift-korelasyon sanısının Euler çarpımı, tarak
uydularının ayna yanında (Δω = −log x) zarfın TAM 1/x olduğunu, tam sayı ve buçuklu
ailelerin aynı zarfta durduğunu ve asal kuvvetlerinin aynada bastırılmadığını öngörür
(f_p(k<0) = p^k, f_p(1) = p/(p−1)², f_p(k≥2) = 0). Bu, hiç görüntülenmemiş
Δω ∈ [−2.12, −1.58] bandında tutuyor mu?

**Ön-kayıt:** KALEM + 197a/b/c + ONKAYIT_197.json (sha256 513327a0…) commit ed9d401,
ölçümden ÖNCE push. Süreç (hepsi kör veri görülmeden): Sonnet türetim teftişi (KRİTİK:
ilk eşikler BK'yı R-g'den ayırmıyordu → profil uyumu + öz-tutarlılık kuralı); dört makine
raporu (K_düz; HAVUZ' τ<0.74; sıkışık tarakta serbest q_f dejenere → çeyrek aile
kalibrasyondan taşınır; jk σ düşük → f_Z = 1.565, f_ρ = 1.374; β KAYIT'a; M6 (ii) karara
giren niceliklere). Ölçüm 82 sn, hüküm betiği ön-kayıt sha'larını denetleyerek koştu.

## Kapılar (hepsi GEÇTİ)

| kapı | sonuç |
|---|---|
| M1 kod yolu | 190 haritası bit-bit (K_ham, K_düz; 70 650 karşılaştırma) |
| M2 kontrol bandı | ana-6 yanar (22.8-57.5σ); +log3/+log2 0.299 (192: 0.310), +log6/+log2 0.614 (0.602), −log2/+log2 0.272 (0.240), −log3/+log3 0.510 (0.466) — hepsi ±%30 içinde |
| M3 öz-terim | tarama ∩ HAVUZ' = ∅ |
| M4 güç | σ_eff^BK = 0.023, σ_eff^Rg = 0.022 (eşik 0.10) |
| M5 çizgi biçimi | artık rms / pencere std = 0.048 (≤ 0.35); −log2 profil korelasyonu 0.9997 |
| M6 (hüküm anı, ölçülen gürültü) | doğru kazanır 0.955 / 0.930, yanlış 0.000 / 0.000 |

## HÜKÜM (eşikler donmuş; kurtarma yok)

| hipotez | ölçülen | hüküm |
|---|---|---|
| **H-197a** zarf, tam sayı ailesi | R_5 = 1.016±0.018, R_6 = 1.014±0.033, R_7 = 1.046±0.031, R_8 = 1.048±0.060 → **R̄_Z^BK = 1.023 ± 0.023** (σ_eff); M_Rg: R̄_Z^Rg = 0.972 vs öngörü R̄^g = 0.706 (> 12σ_eff uzak) | **BK TUTAR** |
| **H-197b** buçuklu aile aynı zarfta | R_{11/2} = 0.813, R_{13/2} = 0.820, R_{15/2} = 0.754 → **ρ = 0.789 ± 0.037** (σ_eff); \|ρ−1\| = 0.211 ≤ 0.25 | **TUTAR** (bant içinde — AMA bkz. dürüstlük notu) |
| **H-197c** −log8 aynada bastırılmaz | **ψ = 0.993 ± 0.039** | **TUTAR** |

KAYIT: β^BK = 0.92 ± 0.11 (BK 1); χ²_Z (4 tam sayı hedefi) BK 3.8 / R-g 551; χ²_ana7 BK 82.9
(buçuklu açığı) / R-g 584. K_ham ile hükümler ve sayılar aynı (K_düz − K_ham ≤ %0.6).
Kalibrasyon içi: A(−log2)·2 = 0.0479±0.0013, A(−log3)·3 = 0.0502±0.0017 (1/x yasası
x = 2-3'te de tutuyor). M_Rg uyumuyla H-197b/c aynı.

**Dürüstlük notu (H-197b):** ön-kayıtlı bant geniş (±0.25) seçilmişti; ρ = 0.789 istatistik
olarak 1'den ~5.7σ_eff AŞAĞIDA. Buçuklu aile BK'nın 2/m'sinin ~%80'inde — görülmüş
kalibrasyon penceresinde de aynı (x = 5/2 ve 7/2: 0.82). Hüküm TUTAR, ama doğru okuma:
"buçuklu aile yanıyor ve aynı 1/x biçimini izliyor; genliği BK'nın ~0.8'i".

## İKİNCİL (görülmüş bant, HÜKÜM DIŞI) — 32 bloklu dar çizgi biçimiyle 192 konumları

Normalizasyon ayna-yanı s (BK'nın çıplak ölçeği), A/(s·c_BK):

| konum | 192 (8 blok, +log2 norm.) | 197 (32 blok, profil uyumu) |
|---|---|---|
| +log2 / +log3 / +log5 / +log7 | (norm.) / 0.84 / 0.64 / 0.61 | 0.795 / 0.720 / 0.608 / 0.530 |
| +log6 / +log10 | 0.82 / 0.75 | 0.654 / 0.584 |
| +log(3/2) | 0.72 | 0.721 |
| **+log(5/3)** | **−0.44 (−9.6σ)** | **+0.424 (7.3σ yanıyor)** |
| **+log(5/2)** | **0.17 (−7.8σ)** | **0.571** |
| **+log(5/4)** | **2.33 (+3.4σ)** | **0.495** |
| −log(3/2) / **−log(4/3)** | 0.82 / **0.31 (−4.5σ)** | 0.849 / **0.788 (21σ)** |
| μ=0: +log4, 8, 9, 4/3, 8/5, 9/2 | −5…−9σ çukur | küçük negatif: −0.0006…−0.0012 (≈ s'nin −%1…−2.5'i; 3.6-5.5σ); +log(8/3), 9/4, 9/5 sıfırla uyumlu |

⇒ 196 ŞERHİNİN (ii) okuması büyük ölçüde DOĞRULANDI: karışık oranlardaki kaba ıskalar
(işaret dönmesi dahil) 8-blok çizgi biçiminin komşu kirliliğiydi.

## ÖN-KAYITSIZ GÖZLEM (veri-SONRASI; ipucu, sınav değil)

Sapmalar yalnız PAYDAKİ (k_p = +1) asallara bağlı, asal-başına bir çarpan gibi:
κ_2 = 0.795±0.015, κ_3 = 0.720±0.012, κ_5 = 0.608±0.018, κ_7 = 0.530±0.020 —
**κ_p ≈ p^{−1/3}** (0.794, 0.693, 0.585, 0.523; %4 içinde). Aynı κ'lar buçuklu aileyi
(κ_2 = 0.795 vs 0.807), çeyrek aileyi (q_f/s = 0.677; κ_3 = 0.72 ile κ_2κ_3 = 0.57 arası),
−log(3/2) (κ_2; 1.07×) ve −log(4/3) (κ_3; 1.09×) açıklıyor. Çarpımsallık kusursuz DEĞİL:
+log6 κ_2κ_3'ün 1.14×'ü, +log10 1.21×, +log(5/4) 0.81×, +log(5/3) 0.70×. Negatif üsler
(f_p(k<0) = p^k) ise TAM tutuyor (R̄_Z = 1.02).
Olası okumalar: (a) sonlu yükseklik düzeltmesi (L = 10.48'de κ_p = 1 − c·log p/L — son
pencerede ölçülebilir), (b) çekirdeğin k = +1 terimlerini farklı ağırlıklandırması (ölçüm
geometrisi), (c) gerçek bir p^{−1/3} yapısı. (a) ile (c) L-bağımlılığıyla ayrılır.

## MANŞET (aday)

> Bogomolny–Keating'in Euler çarpımından türetilen ayna yasası — tarak uydularının ayna
> yanında genlik tam 1/x — hiç görüntülenmemiş bir bantta, parametresiz ve kör olarak
> tuttu: x = 5-8 için R̄ = 1.023 ± 0.023 (rakip çizgi-yoğunluğu modeli > 12σ dışlandı);
> asal kuvveti −log 8 zarfta (ψ = 0.99). Buçuklu aile yanıyor ama genliği BK'nın ~0.8'i.
> Dar çizgi biçimi 196'nın karışık-oran ıskalarını büyük ölçüde sildi; kalan sapma
> paydaki asallara bağlı asal-başına bir çarpan (κ_p ≈ p^{−1/3}, veri-sonrası).

## Açık kalemler

- 198 adayı: κ_p'nin L-bağımlılığı (son pencere, L = 12.45): sonlu-yükseklik mi (κ → 1),
  sabit mi (p^{−1/3})? Ön-kayıtlı, 197 makinesiyle.
- μ=0 konumlarındaki küçük negatif artıklar (s'nin %1-2.5'i): BK'nın ötesi (ikinci mertebe)
  mi, taban modeli mi?

## Teslim

- 197_configs/: 197a_onkayit.py, 197b_harita.py, 197c_hukum.py, ONKAYIT_197.json,
  F_KALIBRASYON.json, M6_on_deneme.json, HUKUM_197.json
- 197_ayna_yasasi.png (kör bant + iki model uyumu; kalibrasyon; ikincil + yan; R_x paneli)
- Veri: scratchpad/197 (harita_omega32.npz, harita_K32.npz, …; gitignore)
