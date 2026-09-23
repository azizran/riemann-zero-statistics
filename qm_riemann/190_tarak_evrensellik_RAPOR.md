# 190 — TARAĞIN EVRENSELLİĞİ: örneklem-dışı sınav (düşük pencere)

**Soru (KALEM_TARAK_EVRENSELLIK_23EYL2026):** 188, pencere-ötesi iptalin
bant-bağımsız bir τ'-tarağı olduğunu ölçtü. Tepeler Bragg'de ve asal uydularındaydı
ve blok-blok yerel L ile kayıyordu. Uydu koşulu ω' = L_yerel ± log n ise tepeler
**Δω := ω' − L_yerel** ekseninde yükseklikten bağımsız olmalı: Δω = 0, ±log 2
(0.693), ±log 3 (1.099), +log 6 (1.792). Rakip okuma "τ'-sabit"tir: tepeler son
penceredeki τ' konumlarında kalır. O zaman düşük pencerede (L = 10.484) +log 2
tepesi 0.604'e, +log 3 tepesi 0.957'ye, +log 6 tepesi 1.562'ye düşer. İkinci
soru: imzalı öngörücü S_Re, rang-1 profil Re v'yi yeni pencerede de izliyor mu?

Her sayının yanında onu üreten betik adı vardır. Tek dalga koşuldu. Ölümler
kurtarılmadı. Git'e dokunulmadı. Sonuç ORTAK TEFTİŞE sunulur, commit sonra.

---

## Ölçüm tanımı (K0b'de donmuş; 188b/188d/188f makinesi AYNEN)

- **Pencere "dusuk":** zeros6 Z[200000:500000] (155_kos.PENCERE).
  t ∈ [139 502.6, 319 387.2]. N = 299 999 orta nokta. L = 10.483929541.
  L ∈ [10.008, 10.836].
- **Zincir:** 155 η önbelleği → 184b K1 → 185b öz → 186b G1_proj (ρ≡1) →
  187c ζ defteri. Modüller importlib ile yüklendi. Hiçbir dosya düzenlenmedi.
  Yeni etiket 'gercek_dusuk', çıktılar `scratchpad/190/zincir_dusuk/`'a gitti.
  Çalışma-anı yamaları şunlardı: `K.SCR`, `b184.S155/S184`, `b185.S185`
  (SHA denetiminden SONRA) ve `b185.kinematik` (185b ile 186b içindeki iki
  örnekte). Kinematik yaması yeni etiketi 185b'nin 'gercek' koluyla aynı
  satırlarla okur; tek fark önbellek dosyasıdır.
- **Harita (190b):** 188b modülünün fonksiyonları AYNEN çağrıldı: seri_ve_G,
  izdusum, c_proj, K_matris, G_loo, bant_maskeleri, karisim, dilimleri_kur.
  Pencere-ötesi τ' ∈ (0.86, 1.30], τ' = log q'/L_düşük. Dilimler 88 × 0.005
  (188 ızgarası AYNEN). Pencere HAVUZ τ ∈ [0.45, 0.86). K(b,s), Ĝ_mid (MUTLAK
  m_n), S, W ve S_Re = Σ a'·πτ'·cos(πτ')·Re Ĝ_mid 188'deki gibi tanımlı.
  **188'den tek sapma TAM örneklemdir** (ALT = 1). Görev tanımı alt-örneklem
  gerektirmiyordu; 188 her 3. noktayı kullanmıştı.
- **Δω profili (190b-D):** Her blok b için çizgi q' şu dilime atanır:
  j = floor((log q' − L_b)/0.025 + ½), dilim merkezi j·0.025. Seri yalnız blok b
  noktalarında kurulur (seri_ve_G AYNEN). İzdüşüm de blok içinde yapılır: 188f(c)
  tek-blok deseni, mix TAM örneklem. κ_b(j) = −Re K_HAVUZ(b,j). L_b, blok içi
  mean log(m/2π)'dir (188f(c) AYNEN).
- **Jackknife:** 8-blok loo, se = √(7/8·Σ(θ_i−θ̄)²).

---

## K0 — ZİNCİR + MAKİNE MÜHRÜ  [190k0_zincir.py]

**Makine mührü.** Aynı sarmalayıcı 'son' için geçici etiketle (`gercek_sonM`,
`scratchpad/190/muhur_son/`) koşuldu. Mevcut hiçbir dosyanın üzerine yazılmadı.
Beş dosya da **BİT-BİT** aynı çıktı [190k0]:

| mühür | sonuç |
|---|---|
| 155 eta_son_t0.4_c4000.npz (tüm diziler) | TUTTU |
| 184 K1_gercek.npz (tüm diziler) | TUTTU |
| 185 OZ_gercek.npz (tüm diziler) | TUTTU |
| 186 G1_proj_gercek.npz (re/im/moment, ρ≡1) | TUTTU |
| 187 K2_zeta['gercek'] (mod, açı, se'ler, genlik, m'ler; liste eşitliği) | TUTTU |

**ζ_g(HAVUZ) = 0.3287±0.0023 ∠+179.96°±0.07**. Bu 187 ile hane hane aynıdır.
Yan ürün olarak w_g 8 bandı da 188 K0 ile aynı çıktı: 1.1698±0.0010 … 1.2560
[190k0].

**Düşük pencere zinciri** [190k0]:
- η önbelleği 155'in kendi üreticisiyle kuruldu. nq = 27, var(ds) = 0.16506.
  Tutarlılık için maks|g − (ds+1)·2π/log(m/2π)| = 2.2e-16.
- 184b çizgi evreni τ ≤ 0.86 (L_düşük ile): 1083 çizgi, τ ∈ [0.066, 0.860].
  Bunların **1044'ü** pencere [0.45, 0.86) içinde. Medyan w = 1.2525; son
  pencerede 1.2516'ydı.
- 185b öz: ḡ = 0.5996173, σ_ε = 0.40698. 186b ρ≡1 izdüşümü: 11 s.

---

## K0b — DONMUŞ ÖN-KAYIT  [190a_onkayit.py]

`scratchpad/190/ONKAYIT_190.json`:
- **sha256 = 38d122a260bf5f78…** (tam:
  38d122a260bf5f78c6fbfa9b85fa45a88e5843b881d781cc48c79b94d5beb1f4)
- **damga Wed Sep 23 17:42:39 +03 2026**

Sıralama dosya zamanlarından okunur. K0 bitişleri 17:38:48 (düşük) ve 17:39:58
(son mührü). Ön-kayıt 17:42:39'da yazıldı, ilk harita dilimi 17:43:36'da. Yani
ön-kayıt haritadan ÖNCE yazıldı. Betik ONKAYIT zaten varsa yazmayı reddeder.

Donanlar:
- **Hedefler (hüküm):** 0, +0.6931, +1.0986, +1.7918. **Kayıt:** −0.6931 ve
  −1.0986.
- **τ'-sabit rakipler:** ± log n · L_düşük/L_son = **0.6041 / 0.9575 / 1.5615**
  (KALEM 0.604/0.957/1.562). Kayıt için −0.6041 ve −0.9575. Bragg için sayısal
  rakip yoktur; orada yalnız ±0.05 koşulu uygulanır.
- **ω-dilim ızgarası:** genişlik 0.025, merkezler j·0.025. Δω aralığı blok b için
  [0.86L − L_b, 1.30L − L_b].
- **8 eşit blok** (linspace(0,N,9)). **L_b** = 10.0876, 10.2305, 10.3540,
  10.4627, 10.5599, 10.6477, 10.7277, 10.8014. Blok-içi L yayılımı 0.155 →
  0.071. Bu kinematiktir, sonuç değildir.
- **Tepe tanımı:** Blok başına HAVUZ κ_b profilinin hedef ±0.15 penceresindeki
  maksimumu, dilim merkezi olarak alınır; ara değerleme yoktur. Ardından
  blok-medyan alınır. se, blok-loo jackknife'tır.
- **Eşikler:** H-190a, H-190b ve H-190c aşağıdaki HÜKÜM tablosundaki gibidir.
  Ek kural: H-190c tutarsa H-190a ve H-190b ölür.
- **Ön-kayıtta yazılı tasarım notu:** +log 6 için rakip (1.562) arama
  penceresinin [1.642, 1.942] dışında kalır. Bu yüzden rakip-yakınlık koşulu
  "tepe > 1.677" ile eşdeğerdir. Ayırıcı sınav, ölüm kolu olan +log 2'dir. Onun
  penceresi [0.543, 0.843] rakibi içerir; hedefle rakibin orta noktası 0.6485.
- **KAYIT tanımları:** tepe_b'nin L_b'ye eğimi (evrensel okuma 0, τ'-sabit okuma
  −1). Son profili 188 dosyalarından (yeniden koşu yok). Düşük pencerenin τ'-dilim
  çapraz kontrolü.
- **Makine mühürleri M1-M3** (sonuçları K1'de).
- **Girdi dosyalarının sha'ları:** η, K1, OZ, G1_proj, 188b, 188d, 188f.

---

## K1 — HARİTA (düşük pencere)  [190b_harita.py]

q ≤ e^{1.30·L} = 829 940 aralığında 66 378 asal-kuvvet var. Bunların **65 295'i
pencere-ötesi**. 88 τ'-dilim 6 süreçle 145 s sürdü. Δω profili için blok başına
185-186 ω-dilimi ve 64 görev vardı; 212 s sürdü.

**Makine mühürleri** [190b]:

| mühür | sonuç |
|---|---|
| M1: seri_ve_G vs 187b.katman_serisi (düşük, dilim 0, 47 çizgi, tam örneklem) | maks\|Δ\| = **0.0** |
| M2: 188'in `dilim_0.npz`'i (son, alt-örneklem, 195 çizgi) bu kod yoluyla | dds, Gre, Gim, q, a, τ **bit-bit** |
| M3: her blokta Σ_j K_HAVUZ(b,j) (ω) = Σ_s K_HAVUZ(b,s) (τ') | göreli fark maks **2.2e-16** |

**Harita toplamları (HAVUZ)** [190b; son sütunu 188b]:

| derinlik | düşük | ζ_g'ye oranı | son (188) | ζ_g'ye oranı |
|---|---|---|---|---|
| Σ_{τ'≤1.20} K | 0.2345±0.0046 ∠−179.68°±0.13 | %70.4 | 0.2261±0.0026 ∠−179.82° | %68.8 |
| Σ_{τ'≤1.30} K | **0.2638±0.0019 ∠−179.62°±0.13** | **%79.2** | 0.2566±0.0032 ∠+179.92° | %78.1 |

### Blok-blok tepe konumları (ω-dilim 0.025; HAVUZ κ; hedef ±0.15)  [190c]

| hedef | rakip | tepe_b (blok 0…7) | blok-medyan ± jk | medyan − hedef | medyan − rakip | koşul |
|---|---|---|---|---|---|---|
| **0** | — | +0.050 +0.050 +0.025 +0.025 +0.025 0.000 −0.025 0.000 | **+0.0250** ±0.0000 | +0.0250 | — | GEÇTİ |
| **+log 2 = 0.6931** | 0.6041 | 0.725 ×5, 0.700 ×3 | **+0.7250** ±0.0000 | +0.0319 | +0.1209 | GEÇTİ |
| **+log 3 = 1.0986** | 0.9575 | 1.150 1.150 1.050 1.125 1.100 1.075 1.100 1.125 | **+1.1125** ±0.0331 | +0.0139 | +0.1550 | GEÇTİ |
| **+log 6 = 1.7918** | 1.5615 | 1.850 1.750 1.800 1.825 1.800 1.800 1.800 1.800 | **+1.8000** ±0.0000 | +0.0082 | +0.2385 | GEÇTİ |
| −log 2 (KAYIT) | −0.6041 | −0.650 −0.725 −0.650 −0.700 −0.725 −0.675 −0.675 −0.675 | −0.6750 ±0.0000 | +0.0181 | −0.0709 | (geçerdi) |
| −log 3 (KAYIT) | −0.9575 | −1.050 −1.050 −1.075 −1.075 −1.125 −1.075 −1.075 −1.075 | −1.0750 ±0.0000 | +0.0236 | −0.1175 | (geçerdi) |

**jk se uyarısı:** Birçok satırda medyanın jk se'si 0.0000 çıkıyor. Bunun
sebebi kesinlik değil, medyanın dilim kuantizasyonudur: 7-blok medyanları aynı
dilim merkezine düşüyor. Gerçek konum belirsizliği en az dilimin yarı genişliği
kadardır (±0.0125). Blokların kendi yayılımı tabloda görülüyor. Ön-kayıt bu se'yi
"kaba, hükme girmez" diye tanımlamıştı.

**Ayırıcı güç** [190c]: +log 2'nin sekiz bloktan HİÇBİRİ rakibe doğru düşmüyor.
Tepeler 0.700-0.725'te, rakip 0.604'te, orta nokta 0.6485. +log 3 ve +log 6'da
da bütün bloklar hedefe rakipten yakın.

**Eğim sınavı (KAYIT)** — tepe_b = a + s·L_b, blok-loo jk [190c]. Evrensel okuma
s = 0, τ'-sabit okuma s = −1 verir:

| hedef | s (düşük, ω-dilim) | s (son, 188 τ'-dilim) |
|---|---|---|
| 0 | −0.094±0.022 | +0.155±0.206 |
| +log 2 | −0.042±0.014 | −0.076±0.129 |
| +log 3 | −0.053±0.058 | +0.037±0.164 |
| +log 6 | −0.018±0.080 | −0.152±0.194 |
| −log 2 | −0.009±0.065 | −0.076±0.129 |
| −log 3 | −0.046±0.025 | +0.041±0.143 |

Bütün eğimler 0 çevresinde. −1'den düşük pencerede 12-68 σ, son pencerede
4.4-7.2 σ uzaktalar. τ'-sabit okuma blok düzeyinde de dışlanıyor. Bragg'in küçük negatif eğimi (4σ) KEŞİF (b)'de
tartışılıyor.

**Çapraz kontrol + son profili (KAYIT; τ'-dilim 0.005 çözünürlüğü, aynı tepe
kuralı)** [190c]:

| hedef | düşük, τ'-dilim (0.052 ω) | SON, 188 dosyaları (0.060 ω) |
|---|---|---|
| 0 | +0.0008±0.0050 | −0.0006±0.0058 |
| +log 2 | 0.6916±0.0196 | 0.6930±0.0132 |
| +log 3 | 1.1016±0.0050 | 1.1058±0.0087 |
| +log 6 | 1.7924±0.0196 | 1.7855±0.0073 |
| −log 2 | −0.6879±0.0038 | −0.6904±0.0132 |
| −log 3 | −1.0907±0.0196 | −1.1054±0.0027 |

Daha kaba τ'-dilim profilinde iki pencerenin blok-medyanları hedefleri
≤ 0.008 içinde veriyor (en büyük sapma: düşük −log 3, +0.0079). Ön-kayıtlı ω-dilim profilinde
sistematik küçük bir pozitif kayma var: ortalama +0.016. Bkz. KEŞİF (b).

---

## K2 — İMZALI ÖNGÖRÜCÜ VE YAPISIZ SIFIR  [190c_hukum.py → 188d.analiz AYNEN]

| nicelik (88 τ'-dilim) | düşük | son (188 yeniden koşu) |
|---|---|---|
| **corr(Re v, S_Re)** | **−0.8542 ± 0.0137** | −0.9481 ± 0.0062 |
| corr(Re v, W) | −0.3047 ± 0.0885 | −0.1007 ± 0.0154 |
| \|corr(Re v,W)\| − \|corr(Re v,S_Re)\| | −0.5496 ± 0.0962 | −0.8474 ± 0.0199 |
| rang-1 payı (KAYIT) | 0.8205 ± 0.0098 | 0.7730 ± 0.0239 |
| corr(\|v\|, S) (KAYIT) | +0.6772 ± 0.0980 | +0.1269 ± 0.0398 |
| corr(\|v\|, W) (KAYIT) | −0.3045 ± 0.0855 | −0.1084 ± 0.0175 |

Son sütun, analiz kod yolunun 188 dosyalarından yeniden koşusudur. 188d'nin
sayılarını hane hane veriyor (0.773±0.024; 0.127±0.040; −0.108±0.018;
−0.948±0.006). Böylece analiz makinesi de mühürlü.

S_Re yeni pencerede de Re v'yi imzalı olarak izliyor. Ama izleme daha zayıf:
son penceredeki −0.948'den −0.854'e düşüyor. Bu 6.2σ'lık bir fark (se'ler
bağımsız sayıldı). Yapısız ağırlık W da düşük pencerede daha çok iz bırakıyor
(−0.30'a karşı −0.10), yine de S_Re'nin çok gerisinde kalıyor.

---

## KAYIT — ζ_g, düşük vs son  [190k0 (düşük) · 187c/K2_zeta.json (son)]

| bant | düşük \|ζ\| | açı° | son \|ζ\| | açı° |
|---|---|---|---|---|
| 0.45-0.50 | 0.3370±0.0049 | −179.92±0.41 | 0.3284±0.0032 | +179.55±0.48 |
| 0.50-0.55 | 0.3333±0.0037 | −179.68±0.33 | 0.3230±0.0032 | −179.46±0.65 |
| 0.55-0.60 | 0.3236±0.0023 | −179.41±0.26 | 0.3184±0.0060 | −179.75±0.48 |
| 0.60-0.65 | 0.3247±0.0041 | +179.69±0.47 | 0.3182±0.0030 | +179.33±0.34 |
| 0.65-0.70 | 0.3205±0.0027 | −179.69±0.52 | 0.3184±0.0039 | +179.51±0.55 |
| 0.70-0.75 | 0.3310±0.0033 | −179.83±0.37 | 0.3274±0.0041 | −179.94±0.45 |
| 0.75-0.80 | 0.3454±0.0042 | +179.45±0.40 | 0.3456±0.0052 | −179.31±0.33 |
| 0.80-0.86 | 0.3742±0.0067 | +179.65±0.35 | 0.3811±0.0072 | +179.93±0.38 |
| **HAVUZ** | **0.3330±0.0020** | **−179.93±0.08** | **0.3287±0.0023** | **+179.96±0.07** |

İptalin ~1/3 büyüklüğü ve 180°'si daha düşük yükseklikte de sürüyor. HAVUZ
farkı 0.0043±0.0031'dir (1.4σ). 9 satırın dokuzu da 180°±1° içinde. Bant
profilinin şekli aynı: orta taban 0.65-0.70'te (0.3205 / 0.3184), kesime doğru
yükseliş var (0.3742 / 0.3811). Alt bantlar düşük pencerede ~0.005-0.010 daha
yüksek. Genlik-okuma Γ(HAVUZ) = 0.3377±0.0023 (son 0.3335±0.0026) [190k0].

---

## ÖN-KAYITSIZ KEŞİF (hükümler görüldükten SONRA; hüküm DIŞI)  [190e_kesif.py]

**(a) Profil bütünü.** İki penceredeki HAVUZ κ yoğunluk profillerinin Pearson
korelasyonu Δω ∈ [−1.3, 2.3] aralığında 144 ω-dilimi üzerinden ölçüldü. **Δω
ekseninde 0.980.** τ'-ölçekli eksende, yani son profili Δω·L_son/L_düşük'te
okunduğunda **0.002**. Figürün sol paneli bunu gösteriyor: yalnız dört hedef
değil, BÜTÜN profil Δω ekseninde üst üste biniyor. Buna ±0.4 çifti ve öteki
küçük tepeler de dahil.

**(b) ω-dilim tepelerindeki pozitif kayma.** 6 hedef × 8 blokta ortalama kayma
+0.016±0.004 (naif se). Kayma blok-içi L yayılımıyla korele: Bragg'de +0.88,
+log 2'de +0.75, öteki hedeflerde +0.12…+0.55. Blok 0'da yayılım 0.155 ve Bragg
+0.050 kayıyor; blok 7'de yayılım 0.071 ve kayma 0.000. Blok-içi L dağılımının
çarpıklığı bunu AÇIKLAMIYOR: medyan L − L_b ≤ 0.0011. Aday açıklama şöyle:
blok içinde tepe, yayılım genişliğinde bir plato olur. W = Σa' ω ile arttığı
için (~e^{ω/2}/ω) platonun üst ucu hafif ağır basar ve argmax oraya kayar.
**SINANMADI.** Hükmü etkilemiyor: bütün medyan kaymaları ≤ 0.032 < 0.05 ve
+log 2/3/6'da rakibe ters yönde. Daha kaba τ'-dilim profilinde ortalama kayma
yalnız +0.003.

**(c) ±0.4 çifti.** 188'in "1∓0.035 simetrik çifti" Δω ekseninde iki pencerede
de aynı yerde duruyor: düşük +0.4125/−0.3750, son +0.3983/−0.3979. τ'-sabit
okuma düşük pencerede ±0.347 verirdi. Aday etiket ±log(3/2) = ±0.405 olabilir.
**SINANMADI.** 0.025 çözünürlükte küçük rasyonellerin yoğunluğu yüzünden etiket
tekil değil (188 şerhi).

---

## HÜKÜM (eşikler K0b'de donmuş; kurtarma yok)

| hipotez/kapı | hüküm | dayanak |
|---|---|---|
| **K0** zincir + makine mührü | **GEÇTİ** | 'son' sarmalayıcı zinciri 155/184/185/186/187'yi bit-bit verdi; ζ_g = 0.3287±0.0023 ∠+179.96° [190k0] |
| M1 / M2 / M3 | **TUTTU** | 0.0 / bit-bit / 2.2e-16 [190b] |
| **H-190a** Δω evrenselliği | **MÜHÜR** | blok-medyan tepeleri 0 → +0.025, +log 2 → 0.725, +log 3 → 1.1125, +log 6 → 1.800. Dördü de ±0.05 içinde ve rakipten yakın. +log 2: \|Δh\| = 0.032, \|Δr\| = 0.121 [190c] |
| **H-190b** imzalı öngörücü | **MÜHÜR** | corr(Re v, S_Re) = −0.854±0.014 ≤ −0.80 [190c] |
| **H-190c** yapısız sıfır | **ÖLDÜ** | \|corr(Re v,W)\| = 0.305 < \|corr(Re v,S_Re)\| = 0.854 (fark −0.550±0.096) [190c] |
| KAYIT ζ_g düşük | KAYIT | 0.3330±0.0020 ∠−179.93°±0.08; bant profili son ile aynı şekil [190k0] |
| KAYIT eğim / son profili | KAYIT | eğimler 0 çevresinde (−1 dışlanıyor); son τ'-dilim medyanları hedeflerin ≤0.008 içinde [190c] |

Dürüst okuma için iki sınır var:
- H-190a'nın geçişi 0.025'lik dilim ölçeğindedir. ω-dilim profilinde sistematik
  +0.016 kayma var ve kaynağı açık (KEŞİF b).
- H-190b eşiği geçti ama izleme sondakinden belirgin biçimde zayıf: −0.854'e
  karşı −0.948.

---

## TÜRETİLEN vs ÖLÇÜLEN

**Türetilen (ön-kayıtta, ölçümden önce):**
- Uydu konum koşulu ω' = L_yerel + log n, yani Δω = 0, ±log 2, ±log 3, +log 6.
- τ'-sabit rakip konumları log n · L_düşük/L_son.

İkisi de kinematik öngörüdür. Parametre içermezler ve düşük pencerenin hiçbir
verisine bakılmadan yazıldılar (L_düşük ve L_b yalnız mid'den).

**Ölçülen:** Bütün tepe konumları, bütün genlikler, korelasyonlar, ζ_g ve Σ K.
Tarağın GENLİKLERİ, ζ'nin ~1/3 büyüklüğü ve S_Re izlemesinin gücü türetilmedi.
Bunlar açık. α=1.30 bağı da açık; bu kalemde dokunulmadı.

---

## MANŞET (aday cümle)

> **Tarak örneklem-dışı sınavı geçti: iptal çekirdeğinin tepeleri yükseklikten
> bağımsız bir Δω = ω' − L_yerel ızgarasında oturuyor. Yeni pencerede (L = 10.48;
> son 12.03) ön-kayıtlı blok-medyan tepeler 0.025, 0.725, 1.113, 1.800 çıktı;
> hedefler 0, log 2, log 3, log 6 idi. Dördü de ±0.05 içinde ve τ'-sabit rakipten
> (0.604/0.957/1.562) uzakta. +log 2'nin 8 bloğunun hiçbiri rakibe düşmüyor.
> Blok eğimleri 0 (τ'-sabit okuma −1 isterdi). H-190a MÜHÜR.** İmzalı öngörücü
> S_Re yeni pencerede de Re v'yi izliyor (corr −0.854±0.014 ≤ −0.80; H-190b
> MÜHÜR), ama sondakinden (−0.948) zayıf. Yapısız sıfır öldü (H-190c). İptalin
> kendisi de yükseklikten bağımsız: ζ_g = 0.333∠180° (son 0.329∠180°), bant
> profili aynı şekilde. Ön-kayıtsız keşif: iki pencerenin profilleri Δω ekseninde
> bütünüyle üst üste biniyor (corr 0.980; τ'-ölçekli eksende 0.002).

---

Teslim:
- Bu rapor.
- `190_configs/`: 190k0_zincir, 190a_onkayit, 190b_harita, 190c_hukum,
  190d_figur, 190e_kesif.
- `190_tarak_evrensellik.png`. Sol panel: Δω profili, düşük (ω-dilim 0.025) ile
  son (188 τ'-dilim) normalize; hedefler düz, rakipler kesikli, blok-medyan
  tepeler ▼. Sağ panel: düşük pencere Re v, −S_Re ve W normalize; corr'lar
  başlıkta.
- `scratchpad/190/`:
  - ONKAYIT_190.json
  - zincir_dusuk/: eta_dusuk, K1/OZ/G1_proj_gercek_dusuk, zincir_gercek_dusuk.{json,npz}
  - muhur_son/: geçici 'gercek_sonM' zinciri ve mühür sonucu
  - dilim_0..87.npz
  - harita_proj_dusuk.npz, harita_K_dusuk.npz
  - omega/b*_c*.npz, harita_omega_dusuk.npz
  - K1_harita_dusuk.json, HUKUM_190.json, profiller_190.npz, K_kesif_190.json
  - tüm loglar (log_k0_*, log_190a…e)
