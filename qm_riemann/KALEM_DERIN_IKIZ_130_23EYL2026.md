# KALEM — 191: DERİN İKİZ τ≤1.30 — zarf tam mı kapanır, yoksa gerçek bir kalıntı mı kalır?
(23 Eylül 2026 akşamı — gece inşası; kullanıcı onaylı akşam planı)

## Durum

189: derin ikiz (τ≤1.10, 1.20) zarfı kapatıyor: HAVUZ f = 0.324, 0.509;
σ_ε 0.43135 → 0.42248 → 0.41788 (gerçek 0.40921). D→∞ modelleri (A, B)
zarfın "tamamının" derinlik olduğunu söylüyordu — ama bunlar dış-değerleme.

## ÖN-MÜHÜR — veri-öncesi keşfedilen yasa ve öngörüler (kaptan, 23 Eyl akşamı)

188'in GERÇEK-kinematikli harita toplamları (derinlik-eşli ζ): τ'≤1.00: 0.0758,
≤1.10: 0.1717, ≤1.20: 0.2261, ≤1.30: 0.2566; tam derinlik (187c): 0.3287.
Δz_D := z_D − z_1.00 = 0.0959 / 0.1503 / 0.1808 / (∞) 0.2529.

**Harita-doğrusal yasa (189 verisinden, veri-öncesi fark edildi):**
- Zarf: f_D / Δz_D = 0.324/0.0959 = **3.38**, 0.509/0.1503 = **3.39** — sabit.
- σ_ε kapanış payı g_D = (σ_Hk − σ_D)/(σ_Hk − σ_g): 0.401/0.0959 = 4.18,
  0.608/0.1503 = 4.05 — yaklaşık sabit (hafif azalan).

**Öngörüler (τ≤1.30):**
- f_1.30 = 3.385 × 0.1808 = **0.612** (bant [0.58, 0.64]).
  (Bağımsız kontrol: geometrik artım kestirimi 0.615 — örtüşüyor.)
- σ_ε(1.30): g = (4.05…4.18) × 0.1808 = 0.732…0.756 ⇒ **σ_ε ∈ [0.4140, 0.4160]**
  (merkez 0.4150).

**Yasanın sonsuz-derinlik imaları (sınanan şey bu):**
- Zarf: f_∞ = 3.385 × 0.2529 ≈ **0.86** ⇒ zarfın ~%14'ü derinlikle KAPANMAZ —
  gerçek-ikiz farkının derinlik-dışı bir kalıntısı (189 modelleri A/B ~%100
  diyordu; çelişki).
- σ_ε: g_∞ ≈ 1.02-1.06 ⇒ aralık saçılımı derinlikle gerçeğe TAM yakınsar.

## Hipotezler (ön-kayıtla donar; ölümler kurtarmasız)

- **H-191a — ZARF SÜRER (birincil yön):** f_1.30 > f_1.20 + 2se (HAVUZ) ve
  8 bantta r_1.20 ≤ r_1.30 + 2σ. ÖLÜM: f_1.30 ≤ f_1.20 + 2se (duraklama/geri dönüş).
- **H-191b — HARİTA-DOĞRUSAL YASA (nicel):** f_1.30 ∈ [0.58, 0.64]. Tutarsa
  yasa üç derinlikte geçerli ⇒ f_∞ ≈ 0.86 okuması (zarfın ~%14'ü kalıntı)
  güç kazanır. f_1.30 > 0.64 ⇒ yasa bozulur, tam kapanış yolunda (189 A/B lehine).
  f_1.30 < 0.58 ⇒ yasa bozulur, erken doyum.
- **H-191c — σ_ε YAKINSAMASI:** σ_ε(1.30) < σ_ε(1.20) − 2se ve ≥ 0.40921 − 2se
  (aşma yok); nicel bant [0.4140, 0.4160].
- **K3 — KİMLİK (KAYIT):** dört derinlikte (c, α); f ve g'nin Δz'ye karşı
  doğrusallığı (dört nokta, eğim ± se); yasa tutarsa zarfın derinlik payı
  f_∞ ± se ve kalıntı payı; 189 modelleriyle yan yana.
- **K4 — KİLİT GÜÇ ANALİZİ (veri-öncesi; ölçümden ÖNCE hesaplanır):** 180'in ham
  para birimleri (M, Q_E) için gerçek↔Hkeskin farkının jackknife se'leri
  kullanılarak, "kilit de derinlik" (fark derinlikle zarfla orantılı küçülür) ve
  "kilit tabanı" (fark ~%20 tabanda durur) okumalarının 1.30'daki ayrımı
  hesaplanır. Ayrım ≥ 3σ ise derin ikizlerde (1.10/1.20/1.30) ham para birimleri
  ölçülür ve hüküm verilir; < 3σ ise ölçüm YAPILMAZ, "güç yetersiz" diye kaydedilir
  (kurtarma yok, eşik gevşetme yok).

## Kapılar

- **K0 — ÖN-KAYIT:** bu kalemin commit'i (GitHub'a push'lanmış hash + zaman
  damgası) inşadan ÖNCE. Ek olarak ajan zincirden önce ONKAYIT_191.json yazar.
- **K1 — İNŞA (kaptan, gece):** 189b sarmalayıcısı üzerinden (dosyalar
  DEĞİŞMEDEN; 191b ince sarmalayıcı AD["1.30"]="Hderin130" ekler): önce ızgara
  sınavı (h 0.015 vs 0.0075, ≤5/5000), sonra inşa (~425k çizgi, ~5-6 saat).
  Kapılar: maks|F| ≤ 1e-8, sıralılık TAM, L = 12.02959324.
- **K2 — ZİNCİR + HÜKÜM (ajan, inşa bitince):** 189c/189d makinesi AYNEN
  (Hderin130 eklenerek); makine mührü (Hkeskin 6/6); H-191a/b/c, K3, K4.

Ölümler kurtarmasız; ajan git'e DOKUNMAZ; sonuç ORTAK TEFTİŞE (sabah).
