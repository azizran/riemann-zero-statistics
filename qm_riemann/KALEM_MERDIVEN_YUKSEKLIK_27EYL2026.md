# KALEM — 200-C: MERDİVENİN YÜKSEKLİK BAĞIMLILIĞI — asal izinin sönmesi sabit yapı mı, sonlu-yükseklik etkisi mi?
(27 Eylül 2026 — TÜRETİM TEFTİŞİNDEN GEÇTİ (Sonnet): tahminler bit bit doğrulandı; işlenen
düzeltmeler: (1) çapa aralığı kuralı + M1c'de çapaya en uzak noktalar; (2) H_C'nin hata payı
gerçek A/B SE'si (düz 0.01 değil), bütün hipotezlere ortak σ_sys = 0.005 (N vs N_eff taban
sistematiği, teftişte 0.002–0.006 ölçüldü); (3) H-200C-1 eşiği 5σ → 3.5σ (teftiş MC'si: 5σ
kötümser SE'de %11 güç), M6c ile yeniden kalibre; (4) basamaklar arası kovaryans (ortak t-blokları)
γ uyumuna girer; (5) serbest-asimptot 7-pencere uyumu birincil çapraz denetim (H-200C-4);
(6) C1'in L yayılımı (0.072) M7c'de; (7) pilot sayımlar yapıldı. Veri: LMFDB (Platt) dosyaları
kullanıcı tarafından tarayıcıyla indirildi (27 Eyl). Eşikler M6c'den sonra donar.
MAKİNE RAPORU (ölçümden önce; gerçek veride hiçbir dağılım özeti hesaplanmadı;
`200_configs/200c_MAKINE_RAPORU.md`): Motor `200c_motor.py` — çapalı RS (C0+C1); çapa aralığı
C1 90 / C2 455 / C3 2122 / C4 20 792 (t biriminde), 7381 / 1253 / 239 / 24 çapa; RS terim sayısı
1238 / 3991 / 12 613 / 69 794. Ajan mühür öncesi bir hata yakaladı ve düzeltti (2π modül
yüklenirken düşük hassasiyette donmuştu ⇒ C3/C4'te 10⁻⁴ rad faz hatası). M0c GEÇTİ (sha256,
artanlık, ortalama aralık 1.000000, N(T) farkı ≤ 0.84). M1c/M2c GEÇTİ (1600 mpmath noktası:
|ΔZ| ~ 10⁻⁹, |ΔZ′| ~ 2·10⁻⁸, k₃ farkı ≤ 2.6·10⁻⁴). M3c GEÇTİ (129 vs 513: |ΔM|/M ~ 2·10⁻⁷).
b=1: C1–C3'te bütün 1.5·10⁶ sıfır, C4'te 10⁶ rastgele sıfır (tohum 301; 30 dk kuralı).
M6c (50 tekrar; 10 tekrarlık ilk sürüm `M6c_200C_v1_10tekrar.json` gürültülüydü): f = 0.94–1.25;
Δκ^{CUE}(0.2): N=14 −0.0064/+0.0066, 17 −0.0062/+0.0028, 19 −0.0031/+0.0034, 22 −0.0056/+0.0045.
H-200C-1 güç tablosu (4000 sim.): eşik 3.5σ'da H_C altında yanlış-pozitif 0.000, güç H_Sγ 0.25,
H_S1 0.86 ⇒ ikili kural sağlanamadı, eşik KALEM kuralıyla 3.5σ'da DONDU. H-200C-2 karışıklık
matrisi: üç sınıf da %100 doğru. M7c: düzeltmeler ~0 ⇒ uygulanmadı.
**ÖNCELİK KURALI (veriden ÖNCE, M6c güç tablosuna dayanarak):** H-200C-1 γ≈0.5'e karşı zayıf
(%25 güç) olduğundan **manşet H-200C-2'nin sınıfıdır**; H-200C-1 (eğim) ve H-200C-4 (asimptot)
destekleyici olarak raporlanır; H-200C-1 "BELİRSİZ" derken H-200C-2 KESİN bir sınıf verirse bu
çelişki değil, beklenen güç farkıdır. EŞİKLER DONDU.)

## Durum

200-A/B MÜHÜRLÜ (509d7a5, f7462fe). Aritmetik kayma = gözlenen − eğik CUE (N = L):

| basamak | κ₃ kayması (havuz, L ≈ 9.3–11.7) | κ₂ kayması |
|---|---|---|
| b=0 rastgele t, log\|Z\| | +0.2828 ± 0.0093 | −0.0820 ± 0.0024 |
| b=1 sıfırlarda log\|Z′\| | +0.1499 ± 0.0009 | −0.0726 ± 0.0008 |
| b=2 yakın çift log(M/δ̃²) | +0.0859 ± 0.0030 | −0.0988 ± 0.0036 |
| a_k (Keating–Snaith/HKO, önde gelen mertebe, her b) | +0.2337 | −0.0881 |

κ₃'te sönme (+0.28 → +0.15 → +0.086; oran ~0.55) kör geometrik tahminle doğrulandı (200-B).

**Veri-sonrası ipucu (bulgu DEĞİL):** A/B pencere değerleri (b0: +0.298/+0.273/+0.279; b1:
+0.144/+0.146/+0.152; b2: +0.076/+0.084/+0.093) L arttıkça a_k'ya doğru kayıyor gibi.
s_b(L) = a_k + c_b·L^{−γ} ortak-γ uyumu: **γ̂ = 0.48 ± 0.12** (χ² 2.07/5; `200c_tahmin.py`).
Bu doğruysa sönme bir SONLU-YÜKSEKLİK etkisi, L → ∞'da merdiven a_k'da birleşir (198'deki
κ_p ≈ p^{−1/3} gibi bir yükseklik tesadüfü). 200-A/B'nin L aralığı dar (ΔL ≈ 2.3); karar yeni,
hiç görülmemiş ve iki kat yüksek veride verilecek.

## Soru

L = 14–22'de merdivenin κ₃ kaymaları (i) L = 11'deki değerlerinde mi kalıyor (sabit yapı),
(ii) a_k'ya doğru mu kayıyor (sonlu-yükseklik)? İkincisiyse hangi hızla?

## Veri ve pencereler (LMFDB, Platt; ±2⁻¹⁰² kesinlik; hiç görülmedi)

| pencere | dosya | nominal L | içerik |
|---|---|---|---|
| C1 | zeros_8846000.dat | L 14.158–14.230, L̄ 14.194 | t 8 846 000 – 9 509 987 | δ̃<0.1/0.2/0.3: 1480 / **11 515** / 37 751 |
| C2 | zeros_99146000.dat | 16.574–16.580, 16.577 | 99 146 000 – 99 714 542 | 1410 / **11 546** / 38 247 |
| C3 | zeros_997946000.dat | 18.883–18.884, 18.884 | 997 946 000 – 998 445 099 | 1489 / **12 051** / 39 305 |
| C4 | zeros_30599546000.dat | 22.306, 22.306 | 30 599 546 000 – 30 599 968 515 | 1534 / **12 143** / 39 612 |

Her dosyada ilk bloktan başlayan **1.5·10⁶ ardışık sıfır** (`200_configs/200c_veri.py`; sıfırlar
TAMSAYI taban + float64 ofset olarak tutulur — t ≈ 3·10¹⁰'da mutlak float64 çözünürlüğü 7.6·10⁻⁶
olduğundan). Pilot bütünlük (yalnız konumlar, 27 Eyl): ilk örnek sıfırlar bilinen değerlerle
aynı; kesin artan; açılmış ortalama aralık 1.00000; Riemann–von Mangoldt ana terimiyle sayım
farkı +0.32 / +0.84 / +0.45 / +0.24 (S(T) mertebesi). İndirilen dosyaların sha256'ları
ONKAYIT'a girer.
Atıf: LMFDB + Platt (Math. Comp. 2015, doi:10.1090/S0025-5718-2014-02884-6).

## Gözlenebilirler (200-A/B ile AYNI tanımlar)

 - b=0: blok t-aralığında 10⁶ düzgün rastgele t (tohum 300), x = log|Z(t)|.
 - b=1: bloktaki sıfırlar (çalışma süresi izin verirse hepsi, yoksa tohum 301 ile 10⁶
   rastgele), x = log(|Z′(γ_n)|·2π/L_n).
 - b=2: bloktaki bütün δ̃ < 0.2 çiftleri, x = log(M_n/δ̃_n²).
 - Kayma: s_b = k₃ − κ₃^{CUE,b}(L̄) (b=2'de ayrıca − Δκ₃^{CUE}(ε, L̄), M6c); aynı κ₂ için.
 - SE: 64 bitişik t-bloklu jackknife × f (M6c); b=2'de + Δκ² (%100 aktarım belirsizliği).

## Hipotezler — κ₃ kaymaları (birincil; `200c_tahmin.json`, kesin L̄ ile; ± tahmin yayılımı σ_h)

| b | L̄ = 14.194 / 16.577 / 18.884 / 22.306 | H_C (sabit) | H_Sγ (γ = 0.48 ± 0.12) | H_S1 (γ = 1, a_k'ya) |
|---|---|---|---|---|
| 0 | | +0.2827 ± 0.0093 | +0.277 / +0.274 / +0.271 / +0.268 (± 0.007–0.008) | +0.271 / +0.266 / +0.262 / +0.257 (± 0.004–0.007) |
| 1 | | +0.1500 ± 0.0009 | +0.159 / +0.164 / +0.169 / +0.174 (± 0.002–0.005) | +0.168 / +0.177 / +0.184 / +0.192 (± 0.0004–0.0007) |
| 2 | | +0.0859 ± 0.0030 | +0.105 / +0.114 / +0.121 / +0.130 (± 0.005–0.009) | +0.123 / +0.139 / +0.151 / +0.163 (± 0.0014–0.0022) |

H_Sγ ve H_S1 KURAM değil: 200-A/B'nin pencere değerlerinden türetilmiş dışdeğerlemeler
(H_S1'de üs 1 sabit, asimptot a_k; H_Sγ'da ortak üs uydurulmuş). Ayırt edici güç en çok
b=1 (küçük SE) ve b=2'de (büyük ayrılık).

κ₂ (ikincil): üç hipotez de −0.08…−0.10 arası, farklar ≤ 0.007 ⇒ yalnız KAYIT.

## Karar kuralları (taslak — teftiş sonrası donacak)

 - **H-200C-1 (yükseklik bağımlılığı):** yalnız YENİ pencerelerin 12 κ₃ kaymasına
   s_b(L) = a_k + c_b·L^{−γ} (c_b serbest, ortak γ) genelleştirilmiş en küçük kareler ile uydurulur
   (kovaryans: pencere başına ortak t-bloklu jackknife'tan 3×3 basamaklar-arası blok, f-kalibre,
   + köşegene σ_sys²) ⇒ γ̂ ± σ_γ. **γ̂/σ_γ ≥ 3.5 ⇒ SONLU-YÜKSEKLİK**; |γ̂| ≤ 2σ_γ ⇒ SABİT YAPI;
   arası BELİRSİZ. (Eşik M6c'de gerçek SE'lerle sınanır: H_C altında yanlış-pozitif < %1 ve
   H_Sγ altında güç ≥ 0.80 olmalı; olmazsa eşik [3, 5] içinde M6c kuralıyla ayarlanıp DONAR.)
 - **H-200C-2 (sınıf):** 12 kaymanın her hipotezin tablosuna χ²'si,
   σ² = SE² + σ_h² + σ_sys² (σ_h: tablodaki tahmin yayılımı — H_C için A/B ağırlıklı ortalama SE'si;
   σ_sys = 0.005 hepsine ortak: N yerine N_eff = L + c, c ∈ [0, 2] tabanının teftişte ölçülen
   0.002–0.006'lık kayma etkisi); en iyi sınıf, bir sonrakine Δχ² ≥ 9 ⇒ KESİN; en iyinin χ²'si
   > 36 ⇒ HİÇBİRİ.
 - **H-200C-4 (asimptot, birincil çapraz denetim):** 7 pencerenin (A/B 3 + C 4) κ₃ kaymalarına
   s_b(L) = A + c_b·L^{−γ} (ORTAK serbest A ve γ) uyumu ⇒ Â ± σ_A. |Â − a_k| ≤ 2σ_A ve
   γ̂ > 0 (≥ 3σ) ⇒ "a_k'ya yakınsama" TUTAR. H-200C-1 SONLU-YÜKSEKLİK derken bu uyum γ ≤ 0 (2σ)
   verirse ÇELİŞKİ yazılır ve manşet BELİRSİZ'e iner.
 - **H-200C-3 (manşet: "yarıya iner" oranı):** pencere başına r₁ = s₁/s₀, r₂ = s₂/s₁;
   L = 22.3'te r₂ H_C'de 0.57, H_S1'de ~0.85. r₂(22.3) − r₂(11.7)'nin işareti ve büyüklüğü
   raporlanır (karar H-200C-1'e bağlı).

**İKİNCİL (KAYIT):** κ₂ kaymaları; bağımsızlık ve aralık yasası (200-B ikincilleri); Not 1'in N_eff'iyle karşılaştırma.

## Kapılar (K0 — ölçümden ÖNCE)

 - **M0c (veri bütünlüğü):** 13-baytlık kayıtların çözümü; ilk sıfırlar bilinen değerlerle;
   blok sayımları Nt0/Nt1 ile; kesin artanlık; sayım N(T) (Riemann–von Mangoldt) ile ±5
   içinde; ortalama aralık 2π/L ile %1 içinde.
 - **M1c (motor):** çapalı Riemann–Siegel (53 tipi: faz = [θ(t₀) − t₀ log n mod 2π] mpmath ile
   + (θ′(t₀) − log n)·dt + θ″(t₀)dt²/2; C0 + C1). **Çapa aralığı kuralı:** atılan kübik terim
   |θ‴|·dt³/6 ≈ dt³/(12 t₀²) her değerlendirme noktasında ≤ 10⁻¹⁰ rad olacak sıklıkta yeni çapa
   (teftiş: tek çapa C1'de tolerayı 7 mertebe aşar); float64 fazın mutlak hatası da raporlanır.
   Pencere başına 200 rastgele t (yarısı çapaya EN UZAK noktalardan) ve 200 sıfırda mpmath'e
   karşı: |ΔZ| ≤ 2·10⁻⁵, |ΔZ′| ≤ 2·10⁻⁴; iki motorun k₂/k₃ FARKI ≤ 0.002/0.005 (değer raporlanmaz).
 - **M2c (sıfır sağlaması):** |Z(γ_n)| ≤ 10⁻⁴·|Z′(γ_n)|.
 - **M3c (tepe):** 129 vs 513 ızgara, 200 olay/pencere: |ΔM|/M ≤ 10⁻⁶.
 - **M6c (sentetik):** Δκ^{CUE}(ε) tablosu N ∈ {14, 17, 19, 22}; f-kalibrasyonu (gerçek
   olay sayılarıyla); güç: her hipotez altında γ̂ ve sınıf karışıklık matrisi.
 - **M7c (karışım):** pencere içi L yayılımının (C1'de 0.072; diğerlerinde ≤ 0.006) taban ve
   Δκ düzeltmelerine etkisi; > 0.002 ise uygulanır. **M8c** (blok sayısı; 200-B kuralı).
   **M9c** (körlük: ölçümden önce gerçek veride dağılım özeti yok).

## Süreç

KALEM → Sonnet türetim teftişi → kullanıcı dosyaları indirir → makine (M0c–M9c) → eşikler
donar → ONKAYIT_200C push → ölçüm → ortak teftiş → mühür.
