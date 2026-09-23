# 192 — UYDULARIN TEORİSİ: tarak uydusu kalıntı sınıfını okur mu? (μ(a)/φ(a) seçim kuralı)

**Soru (KALEM_UYDU_TEORISI_23EYL2026):** 188-190, pencere-ötesi iptalin orta-nokta
örgüsünün Bragg tarağı ve uydularıyla taşındığını gösterdi. Uydular Δω = ω' − L_yerel
ekseninde log 2, log 3, log 6'da duruyor ve yükseklikten bağımsız (190 MÜHÜR).
Uyduların KONUMLARI anlaşıldı, PARLAKLIKLARI anlaşılmadı. Kalem bunun için veri-öncesi
bir durağan-faz teorisi öneriyor ve üç sınav koyuyor:
- **B:** hangi uydular yanar, hangileri söner (seçim kuralı).
- **A1:** kalıntı-sınıfı dönmesi (en keskin, kör sınav).
- **A2:** parlaklık oranları.

Her sayının yanında onu üreten betik adı vardır. Tek dalga koşuldu. Ölümler
kurtarılmadı. Git'e dokunulmadı. Derin ikiz 1.30 inşasına (191) dokunulmadı. Bütün
koşular `nice -n 19` ile, tek süreçle ve OMP/BLAS = 1 ile yapıldı. Sonuç ORTAK TEFTİŞE
sunulur, commit sonra.

---

## Türetimin kısa özeti (kalemden; kaptanın veri-öncesi cebiri)

1. Seviye koşulu N(t_n) = n − ½, N = N̄ + S, S(t) = −Σ_q a_q sin(ω_q t),
   a_q = Λ(q)/(π√q log q).
2. Jacobi–Anger açılımı: e^{−2πiS} = Π_q Σ_k J_k(2πa_q) e^{ikω_q t}.
3. {k_q} teriminin durağan noktası L(t*) = ω' − log n'dir; burada
   n := Π q^{−k_q} ve n rasyonel olabilir. Uydular bu yüzden **Δω = log n**'de
   durur. 190'ın evrenselliği buradan gelir.
4. Durağan fazın değeri, n = a/b sade kesir için
   **φ* ≡ 2π q' b/a − 7π/4 (mod 2π)**. Faz q' mod a'ya bağlıdır.
5. Çizgiler a'ya asal kalıntı sınıflarına eşit dağılır. Bu yüzden uydunun koherens
   çarpanı Ramanujan toplamıdır: c_a(1)/φ(a) = **μ(a)/φ(a)**. μ(a) = 0 olan uydular
   TÜM mertebelerde söner.
6. Birinci mertebede artı-taraf parlaklığı ∝ |μ(n)|/φ(n), eksi-taraf parlaklığı
   ∝ 1/m. İşaret hepsinde aynıdır (iptal).

Ölçülebilir öngörüler:
- **B:** yanmalı ve sönmeli konum kataloğu.
- **A1:** sınıflar arası açılar. +log 3'te sınıflar 120°±15° ayrık olmalı; kontrol
  −log 3'te aynı yönde olmalı. +log 5'te 72° adımlar, +log 6'da 120° beklenir.
- **A2:** oranlar, bant ×/÷ 2 ile.

---

## Ölçüm tanımı (K0'da donmuş; 188b/190b makinesi AYNEN)

**B, birincil (DÜŞÜK pencere, L = 10.484).** 190b'nin ω-dilim (0.025) blok
dosyaları kullanıldı (`190/omega/b*_c*.npz`, yeni harita yok). Hesap şöyle:
- **Havuz profili.** Her ω-dilimi j için yalnız o dilimi TAM kapsayan bloklar B_j
  alınır:
  - C_j = 2·Σ_{b∈B_j}(re+i·im)_{b,j} / Σ_{b∈B_j} n_b
  - K_j = Σ_{q∈HAVUZ} C_j·conj(mix_q) / Σ|mix_q|²
  - κ = −Re K
  - mix = 188b.karisim, 'dusuk' zinciri.
- **Jackknife.** 8-blok loo; C'de blok çıkarılır, mix de loo'dur (188b K_hesap
  AYNEN).
- **Tepe.** Merkezi |c − h| ≤ 0.03 olan dilimlerde κ'nın maksimumu. se o dilimin
  jk se'sidir.
- **Yerel taban.** 0.05 < |c − h| ≤ 0.15 halkasındaki dilimlerin medyanı. Yayılım
  1.4826·MAD'dir. Tabanın üst sınırı medyan + 2·yayılım.
- **Hüküm etiketleri:**
  - "yanar": κ_tepe > 3 se VE üstte.
  - "söner": |κ_tepe| < 2 se YA DA |κ_tepe − medyan| ≤ 2·yayılım.
  - "belirsiz": ikisi de değil. Tabanın altındaki anlamlı değer ayrıca "çukur"
    olarak işaretlenir.

**B, ikincil (SON pencere, KAYIT).** 188 τ'-dilim blok profili kullanıldı (190c
tau_blok_profili AYNEN, mix loo eklendi). Yoğunluk 0.025 ızgarasına np.interp ile
taşındı (190e son_yog AYNEN). Kurallar birincille aynı.

**A1 (SON pencere).**
- Çizgiler τ' = log q'/L_son değerine göre beş pencerede seçildi. Sınıf
  r = q' mod a'dır. gcd(q', a) > 1 olanlar "bölünen" sınıfına gider.
- Seri 188b.seri_ve_G AYNEN hesaplandı. Alt-örneklem her 3. noktadır (100 000),
  188 ile aynı.
- İzdüşüm 188b.izdusum ile yapıldı. K_r = 188b.K_matris[HAVUZ] olarak, 8 loo
  ile hesaplandı.
- Açılar toplam K_top'a göre verilir; sınıf farkları sarılı olarak alınır.

**A2 (DÜŞÜK pencere).** Parlaklık P(h) = Σκ, |c − h| ≤ 0.0625 içindeki 5 dilim
üzerinden alınır. Bu, 188'in ±0.006 τ' uydu penceresinin düşük penceredeki ω
karşılığıdır. Oranların se'si replika oranlarının jk se'sidir.

---

## K0 — DONMUŞ ÖN-KAYIT  [192a_onkayit.py]

`scratchpad/192/ONKAYIT_192.json`:
- **sha256 = bdbcf088e55eb2eb…** (tam:
  bdbcf088e55eb2eb88e58608ddc5a80b6238d7567ec010ca088ade81f47031e6)
- **damga Wed Sep 23 19:02:31 +03 2026**

**Sıralama (dosya zamanlarından):** ön-kayıt 19:02:31'de yazıldı. İlk makine mührü
(192c) 19:04:43'te, ilk profil okuması (192b) 19:05:24'te geldi. Yani ön-kayıt,
PROFİLLERE ve HARİTALARA BAKILMADAN yazıldı. Betik yalnız kinematik bilgi okur:
L, L_b ve asal sayma. Aynı dosya varsa yazmayı reddeder. Bir kuru koşu (`--kuru`)
yazmadan yalnız ızgara ve kapsama listelerini basmıştır.

Donanlar:
- 23 konumluk B kataloğu: Δω, n = a/b, μ(a), φ(a).
- Her konum için arama, halka ve A2 dilim listeleri; kapsama listeleri (kinematik).
- Yanar, söner ve belirsiz ölçütleri.
- A1 pencereleri, sınıflar, sınıf başına çizgi sayıları ve eşikler.
- Yön-tutarlılık KAYDI.
- Makine mühürleri ve tam-örneklem kontrolü.
- A2 oranları ve bantları.
- H-192a/b/c hüküm kuralları.
- 19 girdi dosyasının sha256'sı; kalem dahil.

**Veri-öncesi ızgara notu (ön-kayıtta yazılı):** +log(8/5)'in ±0.03 araması, 0.500
dilimini +log(5/3) ile paylaşıyor. Hüküm literal kuralla verildi. Çakışmasız arama
KAYIT olarak yan yana raporlandı; sonuç aynı çıktı, çünkü tepe 0.450'de.

A1 çizgi evreni (asal sayma) [192a]:

| uydu | τ' penceresi | çizgi | sınıflar | bölünen |
|---|---|---|---|---|
| +log 3 | [1.078, 1.104] | 12 009 | r1 6 007 / r2 6 001 | 3¹² |
| −log 3 | [0.896, 0.922] | 1 616 | r1 815 / r2 800 | 3¹⁰ |
| +log 5 | [1.121, 1.146] | 18 498 | r1 4 653 / r2 4 590 / r3 4 646 / r4 4 609 | — |
| +log 6 | [1.137, 1.161] | 21 086 | r1 10 544 / r5 10 541 | 2²⁰ |
| +log 2 | [1.045, 1.070] | 7 948 | tek (ek KAYIT: mod-4 r1 3 986 / r3 3 962) | — |

---

## MAKİNE MÜHÜRLERİ  [192c_sinif.py, 192b_katalog.py]

| mühür | sonuç |
|---|---|
| 188 `dilim_44.npz` (τ' ∈ (1.080,1.085], 2 070 çizgi, +log 3 penceresi içinde) aynı liste + seri_ve_G ile | dds, Gre, Gim, q, a, τ **BİT-BİT** [192c] |
| tek-dilim K (50 bant) vs 188 `harita_K_gercek.npz` K[:,44] | maks\|Δ\| = 1.0e-17 (tam), 1.3e-17 (loo) [192c] |
| doğrusallık: dilim-44 çizgilerinin mod-3 sınıf serileri toplamı = dilim serisi | maks\|Δ\| = 9.7e-17 (maks\|seri\| 0.099) [192c] |
| düşük havuz profili: T_{b,j}/n_b ≡ 190 `harita_omega_dusuk` K[b, HAVUZ, j] | maks\|Δ\| = 4.3e-17 [192b] |
| son blok profili ≡ 190 `profiller_190` ts_kap / ts_dw | 0.0 / 0.0 [192b] |
| kapsama listeleri ≡ ön-kayıt (23 konum, düşük + son) | assert geçti [192b] |

**Tam-örneklem kontrolü (KAYIT)** [192c]. Alt-örneklem yansız çıktı:

| sınıf | K_alt | K_tam | \|Δ\|/se_alt |
|---|---|---|---|
| +log 3 r1 | −0.00899−0.00077i (∠−175.1°) | −0.00967−0.00070i (∠−175.9°) | 0.98 |
| −log 3 r1 | −0.00578+0.00011i | −0.00575+0.00001i | 0.07 |

---

## K1 — B KATALOĞU  [192b_katalog.py]

DÜŞÜK pencere (birincil) tablosu. κ değerleri ×10⁻³, ± 8-blok jk. Son sütun SON
penceresinin hükmüdür (KAYIT, aynı kurallar).

**Yanmalı (μ(a) ≠ 0):**

| hedef | Δω | a | tepe | κ_tepe | κ/se | taban medyan ± yayılım | hüküm (düşük) | son (KAYIT) |
|---|---|---|---|---|---|---|---|---|
| **+log 2** | 0.693 | 2 | 0.700 | **+19.19 ± 1.55** | +12.4 | −0.73 ± 0.44 | **yanar** | yanar (+10.5σ) |
| **+log 3** | 1.099 | 3 | 1.125 | **+6.02 ± 0.54** | +11.2 | +0.34 ± 0.72 | **yanar** | yanar (+10.2σ) |
| +log 5 | 1.609 | 5 | 1.625 | +1.91 ± 0.22 | +8.8 | −0.07 ± 0.36 | yanar | yanar |
| **+log 6** | 1.792 | 6 | 1.800 | **+11.85 ± 1.02** | +11.6 | +0.39 ± 0.90 | **yanar** | yanar (+10.0σ) |
| +log 7 | 1.946 | 7 | 1.975 | +1.14 ± 0.24 | +4.8 | +0.56 ± 1.72 | söner (taban içinde) | söner |
| +log 10 | 2.303 | 10 | 2.325 | +4.47 ± 0.37 | +12.0 | −0.06 ± 0.64 | yanar | yanar |
| **+log(3/2)** | 0.405 | 3 | 0.425 | **+2.58 ± 0.34** | +7.7 | −0.68 ± 0.49 | **yanar** | yanar (+13.5σ) |
| +log(5/3) | 0.511 | 5 | 0.525 | −0.44 ± 0.15 | −2.9 | +1.88 ± 1.51 | söner (taban içinde) | söner |
| +log(5/2) | 0.916 | 5 | 0.900 | +0.25 ± 0.16 | +1.6 | −0.53 ± 0.30 | söner (\|κ\|<2se) | belirsiz |
| +log(5/4) | 0.223 | 5 | 0.200 | +1.75 ± 0.29 | +6.1 | −0.72 ± 0.75 | yanar | yanar |
| **−log 2** | −0.693 | 1 | −0.675 | **+4.60 ± 0.61** | +7.6 | −0.30 ± 0.34 | **yanar** | yanar (+7.7σ) |
| **−log 3** | −1.099 | 1 | −1.075 | **+3.37 ± 0.36** | +9.3 | −0.26 ± 0.55 | **yanar** | yanar (+14.9σ) |
| −log(3/2) | −0.405 | 2 | −0.425 | +5.23 ± 0.52 | +10.0 | +0.15 ± 0.33 | yanar | yanar |
| −log(4/3) | −0.288 | 3 | −0.275 | +0.55 ± 0.28 | +2.0 | +0.74 ± 2.08 | söner (\|κ\|<2se) | söner |

**Sönmeli (μ(a) = 0):**

| hedef | Δω | a | tepe | κ_tepe | κ/se | taban medyan ± yayılım | hüküm (düşük) | son (KAYIT) |
|---|---|---|---|---|---|---|---|---|
| +log 4 | 1.386 | 4 | 1.400 | **−1.08 ± 0.12** | **−8.7** | −0.18 ± 0.22 | **belirsiz — ÇUKUR** | belirsiz — çukur (−7.2σ) |
| +log 8 | 2.079 | 8 | 2.050 | −0.10 ± 0.15 | −0.7 | −0.02 ± 0.86 | söner (\|κ\|<2se) | söner |
| +log 9 | 2.197 | 9 | 2.175 | −0.39 ± 0.18 | −2.2 | +0.46 ± 1.48 | söner (taban içinde) | söner |
| +log(4/3) | 0.288 | 4 | 0.300 | **−1.20 ± 0.21** | **−5.8** | +1.52 ± 1.25 | **belirsiz — ÇUKUR** | belirsiz — çukur (−4.3σ) |
| +log(8/3) | 0.981 | 8 | 0.975 | −0.50 ± 0.12 | −4.3 | +1.36 ± 2.66 | söner (taban içinde) | söner |
| +log(8/5) | 0.470 | 8 | 0.450 | +0.61 ± 0.39 | +1.6 | −0.49 ± 0.88 | söner (\|κ\|<2se) | söner |
| +log(9/2) | 1.504 | 9 | 1.475 | −0.01 ± 0.12 | −0.1 | +0.27 ± 1.94 | söner (\|κ\|<2se) | söner |
| +log(9/4) | 0.811 | 9 | 0.825 | −0.87 ± 0.11 | −7.9 | +2.09 ± 3.53 | söner (taban içinde) | söner |
| +log(9/5) | 0.588 | 9 | 0.575 | −0.92 ± 0.15 | −6.3 | +5.92 ± 9.74 | söner (taban içinde) | söner |

Ana-6 (kalın satırlar) düşük pencerede de son pencerede de yanıyor: 7.6-12.4σ [192b].

**Hiçbir μ(a) = 0 konumu yanmıyor.** Sönmeli aramalardaki en yüksek κ/se değeri
+1.6'dır (+log(8/5)). Hiçbir sönmeli "4se-yanar" ölçütüne ulaşmıyor. Yalın
"κ > 4 se" sayımı da 0/9 [192b].

Sönmelilerin sönme nedenleri iki grupta toplanıyor:
- **3/9 temiz sıfır:** |κ| < 2 se. Bunlar log 8, log(8/5), log(9/2).
- **4/9 geniş taban:** "söner" hükmü yalnız geniş taban bandından geliyor. Bunlar
  log 9, log(8/3), log(9/4), log(9/5). Bu dördünde κ_tepe −2.2σ ile −7.9σ arasında,
  ANLAMLI NEGATİF. Halka komşu parlak uyduları içeriyor; örneğin log(9/5)'in
  halkası +log 2 omzunu içeriyor ve yayılım 9.7·10⁻³ oluyor.
- **Kalan 2/9 çukur:** log 4 ve log(4/3) taban bandının da altında. Ön-kayıtlı
  kurala göre belirsiz (çukur).

Yanmalı listenin ana-6 dışındaki sekiz konumunun dördü yanmıyor:
- +log 7: +4.8σ ama taban içinde, çünkü halkası +log 6 tepesinin 1.800 dilimini
  içeriyor.
- +log(5/3): −2.9σ.
- +log(5/2): +1.6σ.
- −log(4/3): +2.0σ.

Bu dört konum hükme girmez.

---

## K2 — A1 SINIF AYRIŞTIRMASI (SON penceresi)  [192c_sinif.py]

Açılar derece cinsinden, ± jk. "Göreli" sütunu sınıfın açısı eksi toplamın açısıdır.

| uydu | K_top (\|K\|, ∠) | κ_top | sınıf | \|K_r\| | ∠K_r | göreli | \|K_r\|/\|K_top\| |
|---|---|---|---|---|---|---|---|
| **+log 3** | 0.01962±0.00076 ∠−179.7°±0.4 | 0.01962±0.00076 | r1 | 0.00902±0.00060 | −175.1±2.2 | **+4.6±2.2** | **0.460±0.024** |
| | | | r2 | 0.01066±0.00061 | +176.3±1.7 | **−3.9±1.7** | **0.543±0.025** |
| | | | bölünen (3¹²) | 1.7e-6 | — | — | 0.000 |
| **−log 3** | 0.01173±0.00203 ∠+178.5°±2.7 | 0.01173±0.00203 | r1 | 0.00578±0.00129 | +178.9±4.1 | +0.4±2.8 | 0.493±0.041 |
| | | | r2 | 0.00595±0.00092 | +178.0±3.5 | −0.5±2.6 | 0.507±0.042 |
| | | | bölünen (3¹⁰) | 1.7e-6 | — | — | 0.000 |
| +log 5 | 0.01454±0.00623 ∠+177.1°±1.0 | 0.01452±0.00623 | r1 | 0.00060±0.00149 | −3.7±139 | +179.1±139 | 0.041±0.149 |
| | | | r2 | 0.00821±0.00175 | −178.9±1.6 | +4.0±2.0 | 0.565±0.160 |
| | | | r3 | 0.00730±0.00187 | +175.4±2.7 | −1.7±2.2 | 0.502±0.120 |
| | | | r4 | 0.00050±0.00095 | +43.9±92 | −133.2±91 | 0.034±0.104 |
| +log 6 | 0.03683±0.00121 ∠−179.8°±0.4 | 0.03683±0.00121 | r1 | 0.01855±0.00068 | +177.9±0.9 | −2.2±0.6 | 0.504±0.013 |
| | | | r5 | 0.01831±0.00083 | −177.5±0.6 | +2.3±0.7 | 0.497±0.012 |
| | | | bölünen (2²⁰) | 4.3e-7 | — | — | 0.000 |
| +log 2 | 0.05295±0.00103 ∠−179.4°±0.4 | 0.05295±0.00103 | tek (mod 2) | = K_top | — | 0.00 | 1.000 |
| | | | ek: r1 mod 4 | 0.02592±0.00058 | −179.3±0.4 | +0.1±0.2 | 0.490±0.006 |
| | | | ek: r3 mod 4 | 0.02703±0.00062 | −179.4±0.4 | −0.1±0.2 | 0.510±0.006 |

**Sınıflar arası açı farkları ve ön-kayıtlı koşullar** [192c, 192d]:

| uydu | fark | ölçülen | öngörü (kalem) | koşul |
|---|---|---|---|---|
| **+log 3** | Δ_{1→2} | **−8.6° ± 3.9** | ±120° ± 15° | **TUTMADI** (ölüm eşiği < 45°: ÖLÜM) |
| | oranlar | 0.460 / 0.543 | [0.7, 1.4] | TUTMADI |
| **−log 3** | Δ_{1→2} | **−0.9° ± 5.5** | < 15° | **TUTTU** |
| | oranlar | 0.493 / 0.507 | [0.35, 0.65] | TUTTU |
| +log 5 (KAYIT) | Δ_{1→2}, Δ_{2→3}, Δ_{3→4} | −175±139, −5.7±3.7, −131.5±90 | ±72° ± 15° | TUTMADI |
| | oranlar | 0.041 / 0.565 / 0.502 / 0.034 | [0.6, 1.5] | TUTMADI |
| +log 6 (KAYIT) | Δ_{1→5} | +4.5° ± 1.2 | ±120° ± 15° | TUTMADI |
| | göreli ∓60° yan kaydı | −2.2 / +2.3 | ∓60 ± 15 | TUTMADI |
| +log 2 (kontrol) | tek sınıf göreli | 0.00° | < 15° | TUTTU (özdeşlik) |
| | ek mod-4 | Δ −0.1°±0.4; 0.490 / 0.510 | aynı yön, ≈0.5 | TUTTU |

**Yön-tutarlılık KAYDI anlamsız çıktı.** Kayda göre işaretler "tutarlı" [192d]. Ama
dönme olmadığı için işaretler gürültü düzeyindeki küçük farkların işaretidir.

**Okuma:** Tahmin edilen dönme YOK. +log 3'ün iki sınıfı, kontrol −log 3'ünkiler gibi
toplamın yönünde duruyor ve her biri yaklaşık yarı ağırlık taşıyor. Bu desen +log 6
ve +log 2 mod-4 bölmesinde de aynı. +log 5 ise dönmüyor ama sınıflar EŞİT de değil:
toplamı r2 ve r3 taşıyor (0.57, 0.50); r1 ve r4 sıfırla uyumlu (0.04±0.15,
0.03±0.10). +log 5'in toplam se'si büyük (|K_top| ±%43). Sınıf se'leri daha küçük;
bu da bloklar arası ortak bir dalgalanmaya işaret ediyor.

---

## K3 — A2 PARLAKLIK ORANLARI (DÜŞÜK pencere)  [192d_hukum.py]

P(h) = Σκ, 5 dilim, ± jk [192d]:

| uydu | +log 2 | +log 3 | +log 5 | +log 6 | +log 7 | +log 10 | −log 2 | −log 3 |
|---|---|---|---|---|---|---|---|---|
| P ×10⁻³ | 70.66±1.2 | 21.89±0.50 | 6.32±0.25 | 42.57±0.88 | 3.44±0.28 | 15.32±0.51 | 16.94±0.37 | 10.20±0.62 |

| oran | öngörü | bant | **ölçülen (düşük)** | bantta? | taban-düz. (KAYIT) | tepe oranı (KAYIT) | SON (KAYIT) |
|---|---|---|---|---|---|---|---|
| +log3/+log2 | 0.500 | [0.25, 1.00] | **0.310 ± 0.006** | EVET | 0.272±0.005 | 0.314±0.008 | 0.327±0.014 |
| +log5/+log2 | 0.250 | [0.125, 0.50] | **0.089 ± 0.004** | HAYIR | 0.090±0.003 | 0.100±0.006 | 0.106±0.007 |
| +log6/+log2 | 0.500 | [0.25, 1.00] | **0.602 ± 0.007** | EVET | 0.547±0.007 | 0.618±0.010 | 0.619±0.036 |
| +log7/+log2 | 0.167 | [0.083, 0.333] | **0.049 ± 0.004** | HAYIR | 0.009±0.004 | 0.060±0.010 | 0.050±0.009 |
| +log10/+log2 | 0.250 | [0.125, 0.50] | **0.217 ± 0.005** | EVET | 0.210±0.004 | 0.233±0.007 | 0.245±0.011 |
| −log2/+log2 | 0.500 | [0.25, 1.00] | **0.240 ± 0.004** | HAYIR (sınırda) | 0.248±0.004 | 0.240±0.014 | 0.247±0.021 |
| −log3/+log3 | 0.667 | [0.333, 1.333] | **0.466 ± 0.035** | EVET | 0.570±0.040 | 0.560±0.039 | 0.434±0.021 |

**Bantta: 4/7** [192d].

"Bantta" ifadesi uyum anlamına gelmiyor. Örneğin +log3/+log2 = 0.310±0.006, öngörü
0.5'in 0.62 katıdır. Bant ×/÷2 olduğu için geçiyor.

**Tek asallar |μ|/φ'den HIZLI düşüyor:** log 3, 5, 7 için 0.31 / 0.089 / 0.049,
öngörü 0.5 / 0.25 / 0.17. Kalemin son penceresi için andığı gerilim düşük pencerede
de aynen duruyor:
- Kalemin son değerleri (0.27, 0.08, 0.62, 0.03) 188'in ±0.006 τ' pencere
  dökümündendir.
- Bu kalemin tanımıyla son değerleri 0.327, 0.106, 0.619, 0.050 [192d].

---

## HÜKÜM (eşikler K0'da donmuş; kurtarma yok)

| hipotez/kapı | hüküm | dayanak |
|---|---|---|
| **K0** ön-kayıt | **GEÇTİ** | sha bdbcf088…, 19:02:31; ilk profil okuması 19:05:24 [192a] |
| makine mühürleri | **TUTTU** | dilim_44 bit-bit; K 1.0e-17; doğrusallık 9.7e-17; ω-profili 4.3e-17; son profil 0.0 [192c, 192b] |
| tam-örneklem kontrolü | **TUTTU** (KAYIT) | 0.98 / 0.07 se_alt [192c] |
| **H-192a** seçim kuralı (B) | **KAYIT (kısmi)** — ÖLÜM YOK | Ana-6'nın hepsi yanıyor (7.6-12.4σ). Sönmelilerin hiçbiri yanmıyor; en yüksek +1.6σ, 4se-yanar 0/9. Ancak 9 sönmeliden 2'si "söner" değil, ÇUKUR: +log 4 −8.7σ, +log(4/3) −5.8σ. Son penceresi aynı: 2 çukur, ölüm adayı yok [192b, 192d] |
| **H-192b** kalıntı dönmesi (A1) | **ÖLDÜ** | +log 3: Δ_{1→2} = −8.6°±3.9 < 45° (ölüm koşulu). −log 3 koşulu tuttu (−0.9°±5.5; 0.49/0.51). Yan öngörüler +log 5 ve +log 6 de tutmadı [192c, 192d] |
| **H-192c** parlaklık (A2) | **KAYIT** | 4/7 oran bantta; tek asallar öngörüden hızlı düşüyor; −log2/+log2 = 0.240 bandın hemen altında [192d] |

**Dürüst okuma:**
- **A1 ölümü temiz ve kör.** Hiçbir sınıf dönmüyor. +log 3 ile −log 3 aynı deseni
  veriyor: iki sınıf toplam yönünde, her biri ~½. Kalemin 4. adımındaki faz
  "2π q' b/a" gözlenebilir K_r'de GÖRÜNMÜYOR.
- **B'nin "yanma" yüzü güçlü tuttu.** μ(a) = 0 konumlarının hiçbiri yanmıyor. Ama
  "söner" etiketi 4 konumda yalnız geniş taban bandı sayesinde kazanıldı ve bu
  konumlar da anlamlı negatif. Yani μ(a) = 0 uyduları "sıfır" değil, çoğunlukla
  NEGATİF (bkz. KEŞİF b). Kalemin "tüm mertebelerde koherens sıfır" iddiası bu
  ölçüde desteklenmiyor. Ön-kayıtlı ölüm yalnız pozitif yanma içindi ve o
  gerçekleşmedi.
- **Ön-kayıtlı taban tanımı kalabalık bölgelerde geniş.** Halka komşu parlak uyduları
  içeriyor. Bu, sönmeli hükmünü kolaylaştırdı (4 konum) ve yanmalı +log 7'yi
  söndürdü. Bu bir tasarım sınırıdır; kural değiştirilmedi.

---

## ÖN-KAYITSIZ KEŞİF (hükümler görüldükten SONRA; hüküm DIŞI)  [192f_kesif.py]

**(a) Sınıflar dönmüyor, gerçek eksende "cos" deseniyle diziliyor. SINANMADI.**

Her sınıfın toplam yönündeki imzalı payı f_r = Re(K_r·conj K_top)/|K_top|² olarak
ölçüldü (Σf_r = 1). Aday desen f_r ∝ cos(2π r b/a); Σ_r cos(2πr b/a) = μ(a), yani
koherens çarpanı kalemdekiyle aynıdır:

| uydu | ölçülen f_r (± jk) | cos deseni |
|---|---|---|
| +log 3 | r1 +0.458±0.025, r2 +0.542±0.025 | +0.5, +0.5 |
| −log 3 | r1 +0.493±0.042, r2 +0.507±0.042 | +0.5, +0.5 |
| +log 6 | r1 +0.503±0.012, r5 +0.497±0.012 | +0.5, +0.5 |
| +log 5 | r1 −0.041±0.155, r2 +0.563±0.158, r3 +0.502±0.120, r4 −0.024±0.122 | −0.309, +0.809, +0.809, −0.309 |

Dik paylar küçük: |g_r| ≤ 0.04.
- +log 3, −log 3 ve +log 6'da desen tutuyor.
- +log 5'te işaret deseni tutuyor: r2 ve r3 taşıyor; r1 ve r4 negatif yönde ama
  küçük. Büyüklükler 1.5-2.6σ sapıyor.

Aday açıklama: seri gerçel bir fonksiyon, cos(ω'm) = (e^{iω'm} + e^{−iω'm})/2. Bu
yüzden her sınıf r, eşleniği a−r ile birlikte görünür ve e^{2πirb/a} yerine
cos(2πrb/a) taşır. **SINANMADI.** Bu yeni bir ön-kayıt gerektirir. Keskin bir aday:
+log 7'nin r → cos(2πr/7) deseni, yani (r1,r6) : (r2,r5) : (r3,r4) =
−0.62 : +0.22 : +0.90 (toplama göre). Bir diğeri: +log 5'te r1 ve r4'ün toplama
TERS işareti.

**(b) μ(a) = 0 konumlarında negatif çukurlar. SINANMADI; kaba taban
karşılaştırması.**

Sönmeli konumların ±0.05 içindeki en düşük κ/se değerleri:

| pencere | +log 4 | +log 8 | +log 9 | +log(4/3) | +log(8/3) | +log(8/5) | +log(9/2) | +log(9/4) | +log(9/5) |
|---|---|---|---|---|---|---|---|---|---|
| düşük | −5.4σ | −8.2σ | −8.3σ | −13.6σ | −3.1σ | −7.6σ | −1.1σ | −7.9σ | −9.8σ |
| son | −4.6σ | −4.6σ | −4.6σ | −6.4σ | −1.7σ | −0.3σ | +0.2σ | −7.8σ | −9.2σ |

Taban karşılaştırması: Δω ∈ [−1.3, 2.3] içinde rastgele bir h için aynı büyüklüğün
≤ −5σ olma olasılığı düşükte 0.31, sonda 0.15. Sönmelilerde bu sayı düşükte 7/9,
sonda 3/9. ≤ −3σ için oranlar: düşük 8/9 (taban 0.46), son 6/9 (taban 0.31).

Yani μ(a) = 0 konumları rastgeleden sık negatif çukur taşıyor. En belirgin olanlar
log(9/5), log(9/4) ve log(4/3); bunlar iki pencerede de −6σ'dan derin. Konumlar
bağımsız değil ve ±0.05 penceresi komşu tepelere dokunabiliyor. Bu yüzden bu bir
ipucudur, sınav değildir.

**(c) A2 ek gözlemleri** [192f]:
- 2m/m oranları: P(+log 6)/P(+log 3) = 1.944±0.038 (son 1.894±0.089) ve
  P(+log 10)/P(+log 5) = 2.424±0.133 (son 2.308±0.162). 2 çarpanı parlaklığı
  yaklaşık İKİYE katlıyor. |μ|/φ bu oran için 1 öngörür (φ(2m) = φ(m)).
- Düşük ile son pencere arasındaki A2 oran farkları −2.3σ ile +0.8σ arasında. Yedi
  oranın beşi 1.1σ içinde; en büyük sapmalar +log10/+log2 −2.3σ ve +log5/+log2
  −2.1σ. Uyduların GÖRELİ parlaklıkları kabaca yükseklikten bağımsız; bu, 190'ın
  konum evrenselliğinin genlik karşılığıdır.

---

## TÜRETİLEN vs ÖLÇÜLEN

**Türetilen (veri-öncesi, kaptan; kalem):**
- Seçim kuralı μ(a)/φ(a).
- Faz φ* ≡ 2πq'b/a − 7π/4 ve kalıntı-sınıfı dönmesi.
- Birinci mertebe parlaklıklar |μ(n)|/φ(n) ve 1/m.

A1, hiçbir verisine bakılmadan yazılmış en keskin öngörüdür: sınıf bölmesi daha önce
hiç yapılmamıştı.

**Körlük sınırları:**
- **B kataloğu tam kör değil.** Son penceresi figürü kaptan tarafından kaba
  çözünürlükte görülmüştü (188 τ'-dilim, 0.060 ω). 190 figürünün sol panelinde düşük
  pencerenin ω-dilim profili de normalize olarak yer alıyordu (190 raporu, Teslim).
  Hangi konumların parlak göründüğü kalem yazılmadan önce kısmen biliniyordu.
- **A2'nin son pencere oranları görülmüştü** (kalemde yazılı). Sınav bu yüzden düşük
  pencereye taşındı. Ama düşük pencere oranları sonunkilerle aynı çıktı (KEŞİF c).
  Bu yüzden A2 fiilen yarı-kör bir sınavdır.
- **Tek tam kör sınav A1'dir, ve öldü.**

**Ölçülen:** κ profilleri, bütün K_r ve açılar, oranlar, çukurlar ve cos deseni.

**Türetilmemiş ve açık kalanlar:**
- Tek asalların hızlı düşüşü (log 3, 5, 7: 0.31, 0.089, 0.049).
- 2 çarpanının yaklaşık ×2 etkisi.
- μ(a) = 0 çukurlarının işareti.
- Sınıfların cos-tipi dizilişi.
- α = 1.30 bağı; bu kalemde dokunulmadı.

Landau–Gonek akrabalığı ve literatürdeki yeri kalemdeki gibi KONTROL EDİLMEDİ.
Yenilik iddiası yok.

---

## MANŞET (aday cümle)

> **Uydular kalıntı sınıfına göre DÖNMÜYOR: en keskin kör öngörü öldü.** +log 3
> uydusunun q' mod 3 sınıfları 120° değil, −8.6°±3.9 ayrık çıktı. Her sınıf toplamın
> yarısını taşıyor (0.46/0.54), yani kontrol −log 3 ile aynı desen (−0.9°, 0.49/0.51).
> +log 5 ve +log 6 de dönmüyor. H-192b ÖLDÜ. Seçim kuralının "yanma" yüzü ise duruyor:
> μ(a) = 0 olan 9 konumun hiçbiri yanmıyor (en yüksek +1.6σ), ana-6'nın altısı da
> 7.6-12.4σ yanıyor. Ama μ(a) = 0 konumları sıfır değil, çoğunlukla NEGATİF; log 4 ve
> log(4/3) ön-kayıtlı tabanın da altında (−8.7σ, −5.8σ). H-192a KAYIT (kısmi), ölüm
> yok. Parlaklıklar |μ|/φ'den hızlı düşüyor (log3/log2 0.31, log5 0.089, log7 0.049;
> 4/7 bantta; H-192c KAYIT) ve iki pencerede aynı. Ön-kayıtsız: sınıf katkıları
> dönmek yerine gerçek eksende cos(2πrb/a) işaretleriyle diziliyor, ki bu μ(a)
> koherensini korur. Kalemin fazı eşlenik çiftlerle ortalanıyor olabilir; SINANMADI.

---

Teslim:
- Bu rapor.
- `192_configs/`:
  - 192a_onkayit
  - 192b_katalog
  - 192c_sinif
  - 192d_hukum
  - 192e_figur
  - 192f_kesif
- `192_uydu_teorisi.png`:
  - Sol üst: DÜŞÜK Δω profili (± jk). Yanmalı konumlar yeşil düz çizgi ve ▲,
    sönmeli konumlar kırmızı kesikli çizgi ve ▼. Boş işaret, hükmün öngörüyle
    uyuşmadığını gösterir; çukurlar etiketli.
  - Sol alt: SON profili ile düşük profilin normalize karşılaştırması.
  - Sağ: beş A1 kutupsal paneli. Siyah ok K_top'tur; renkli oklar sınıf
    K_r/|K_top| değerleridir, açı toplama göredir. Gri kesikli çizgiler KALEM
    öngörüsüdür.
- `scratchpad/192/`:
  - ONKAYIT_192.json
  - muhur_192.json
  - K1_katalog.json, profiller_192.npz
  - seri_*.npz (15 sınıf serisi), sinif_proj.npz
  - tamseri_*/tam_*.npz (tam-örneklem kontrolü)
  - K2_sinif.{json,npz}
  - HUKUM_192.json, K_kesif_192.json
  - loglar (log_192a…f)
