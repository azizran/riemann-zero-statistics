# KALEM — 200-B: YAKIN ÇİFT BASAMAĞI (b=2) — asalların imzası iki sıfır birleşirken ne olur?
(26 Eylül 2026 — kullanıcı onaylı iki aşamalı tasarımın ikinci aşaması; 200-A MÜHÜRLÜ (509d7a5).
TÜRETİM TEFTİŞİNDEN GEÇTİ (Sonnet, 26 Eyl): tahmin aritmetiği bağımsız olarak doğrulandı;
işlenen düzeltmeler: (1) M8b blok-sayısı kapısı eklendi; (2) H-200B-3'ün manşette önceliği
açıkça yazıldı (H-200B-2 ile çelişki riski); (3) sonlu-ε düzeltmesi kestirimi teftişin MC'siyle
güncellendi ve aktarım belirsizliği hata payına eklendi; (4) H_geo/H_lin hata payları 200-A'nın
b0–b1 kovaryansıyla yeniden hesaplandı (`200b_A_kaymalari.json`); (5) nicelik-çarpıklığı ortak
birincili düşürüldü (a_k biçimi için olasılık yasası yok, §6b); (6) yüksek-L ikincili çıkarıldı
(eski motor C0-only; ayrı bir kalemde yeniden doğrulanmadan kullanılmaz). Eşikler M6'dan sonra
donar.
MAKİNE RAPORU (ölçümden önce; gerçek yakın çiftlerde hiçbir dağılım özeti hesaplanmadı;
`200_configs/200b_MAKINE_RAPORU.md`): M1b GEÇTİ (motor+ızgara vs mpmath altın-oran, 900 olay:
|ΔM|/M ≤ 5·10⁻⁵; k₂/k₃ farkı ~10⁻⁷). MÜHÜR ÖNCESİ SAPMA: M2b 33 noktalı ızgarayla KALDI
(33 vs 129: |ΔM|/M 6.9·10⁻⁶…1.4·10⁻⁵ > 10⁻⁶; k₂/k₃ etkisi ~10⁻⁷, karar eşiklerinin 4–5 mertebe
altında) ⇒ kapı gevşetilmedi, IZGARA 129 NOKTAYA ÇIKARILDI; aynı 900 olayda 129 vs 513:
|ΔM|/M ≤ 3.7·10⁻⁷ ⇒ M2b GEÇTİ; M1b 129 ile yeniden: ≤ 5.0·10⁻⁵ GEÇTİ (`M1b_200B_g129.json`,
`M2b_200B_g129.json`). M6b(i): Δκ₂(0.2) −0.0043…−0.0038, Δκ₃(0.2) +0.0032…+0.0039 (N = 9–12,
≥ 4·10⁵ olay/N); bağımsız denetçinin değerleriyle 8 karşılaştırmanın hepsi ≤ 1.6σ uyumlu.
M6b(ii): f = 0.94–1.13. M6b(iii): havuz SE(D̄₃) ≈ 0.0042, SE(D̄₂) ≈ 0.0041 (ön kestirimden iyi);
H-200B-3 karışıklık matrisinde P(doğru sınıf) ≥ 0.9998; P < 0.80 çift yok. M7b: bütün
düzeltmeler ≤ 0.0003 ⇒ uygulanmadı. EŞİKLER DONDU.)

## Durum

 - Teorem (CUE, `TEOREM_KUCUK_ARALIK_ISKELET_26EYL2026.md`): iki özdeğer birleşirken
   M/s² → |Λ_{N−2}|/4, geri kalanlar |Λ|⁴ ile eğik ⇒ log-tepenin kümülantları
   κ_r^{CUE,2}(N) kapalı formda.
 - 200-A (mühürlü): alt basamaklarda (b=0 rastgele t, b=1 sıfırlarda Z′) aritmetik imza
   TUTAR; seçilen sonlu-L modeli **M* = a_k (KESİN)**. Veri-sonrası: aritmetik KAYMA
   (gözlenen − CUE) κ₂'de basamaktan ~bağımsız (b0 −0.0820 ± 0.0024, b1 −0.0726 ± 0.0008),
   κ₃'te basamakla sönüyor (b0 +0.2828 ± 0.0093, b1 +0.1499 ± 0.0009; a_k: +0.2337).
 - SORU: yakın çift basamağında (b=2) asalların çarpıklık imzası ne kadar? İskeletteki manşet
   tahmin (a_k, "çarpıklığın işareti döner") mi, yoksa 200-A'daki sönme sürüyor mu?

## Gözlenebilir

Zeros6'da ardışık γ_n < γ_{n+1}; m_n = (γ_n+γ_{n+1})/2, L_n = log(m_n/2π),
δ̃_n = (γ_{n+1}−γ_n)·L_n/2π. Olay: δ̃_n < ε (**birincil ε = 0.2**; ikincil 0.1, 0.3).
M_n = max_{[γ_n, γ_{n+1}]} |Z(t)| (200-A motoru, C0+C1). **x_n = log(M_n/δ̃_n²).**
k-istatistikleri k₂, k₃ (yansız). Pencereler 200-A ile aynı (olay, m_n'ye göre pencereye düşer):

| pencere | olay sayısı ε=0.1 / **0.2** / 0.3 | L̄ (ε=0.2) | κ₂^{CUE,2}(L̄) | κ₃^{CUE,2}(L̄) |
|---|---|---|---|---|
| W1 | 147 / **1 143** / 3 916 | 9.355 | 0.2768 | −0.0658 |
| W2 | 356 / **2 860** / 9 445 | 10.594 | 0.3181 | −0.0713 |
| W3 | 1 308 / **10 165** / 33 713 | 11.662 | 0.3515 | −0.0754 |

(Sayımlar yalnız sıfır konumlarından: `200_configs/200b_sayimlar.json`.)

**Kaymalar:** D_r := k_r − κ_r^{CUE,2}(L̄_W) − Δκ_r^{CUE}(ε, L̄_W), burada Δκ_r^{CUE}(ε) sonlu-ε
düzeltmesidir (CUE Monte Carlo, aynı ε'de; K0'da kesinleşir; ön kestirim ε=0.2'de
Δκ₃ ≈ +0.004, Δκ₂ ≈ −0.005 — teftiş MC'si, N = 9–10). **Aktarım belirsizliği:** düzeltmenin
zeta'ya aynen geçtiği varsayılır (iskelet W3, sınanmadı); bu yüzden her pencerede
σ_ε,r = |Δκ_r^{CUE}(ε)| (düzeltmenin %100'ü) SE'ye karesel eklenir. Havuz: D̄_r = pencerelerin ters-varyans ağırlıklı ortalaması
(200-A'da kaymalar pencereden bağımsızdı).

## Hipotezler (havuzlanmış kaymalar)

| hipotez | kaynak | D̄₃ | D̄₂ | κ₃ (W2 için) |
|---|---|---|---|---|
| **H_sabit** (a_k; 200-A'nın taşıdığı M*) | önceden kayıtlı taşıma | +0.2337 | −0.0881 | +0.162 |
| **H_geo** (geometrik sönme: s₂ = s₁²/s₀) | 200-A'dan bilgilenmiş | +0.0794 ± 0.0028 | −0.0644 ± 0.0020 | +0.008 |
| **H_lin** (doğrusal sönme: s₂ = 2s₁ − s₀) | 200-A'dan bilgilenmiş | +0.0170 ± 0.0095 | −0.0633 ± 0.0024 | −0.054 |
| **H_RMT** (salt eğik CUE) | iskelet | 0 | 0 | −0.071 |
| H_Gauss (rapor için) | iskelet | ≈ −0.03 (W1 −0.038, W2 −0.032, W3 −0.028) | W1 −0.043, W2 −0.084, W3 −0.118 | −0.104 |

H_geo ve H_lin, 200-B verisi görülmeden, yalnız 200-A'nın mühürlü sayılarından türetildi.
Hata payları A'nın ortak-blok jackknife kovaryansıyla yayıldı: κ₃'te b0–b1 korelasyonu
ρ = −0.003 (etkisiz), κ₂'de ρ = +0.32. **Çerçeve uyarısı:** bunlar KURAM DEĞİL, iki noktadan
iki parametreli dışdeğerlemedir; A verisini tanım gereği tutarlar. Sınadıkları şey yalnız
"sönme hangi hızla sürüyor" sorusudur. "κ₃ (W2)" sütunu κ₃^{CUE,2}(L̄) + D̄₃'tür (sonlu-ε
düzeltmesi dahil değil; ~0.004). H_lin ile H_RMT'nin farkı (0.017) beklenen hata payından küçük
⇒ önceden **SÖNMÜŞ** sınıfında birleştirilir.

## Karar kuralları

SE_D: pencere başına t-bloklu jackknife (64 bitişik t-bloğu; olaylar m_n'ye göre), × f (M6);
havuzlanmış SE ters-varyans birleşimi. Ön kestirim (ürün yasası simülasyonu): salt eğik CUE
altında havuz SE(D̄₃) ≈ 0.005, sağa çarpık (Euler p=2,3 eklenmiş) seçenekte ≈ 0.014.

 - **H-200B-1 (b=2'de aritmetik imza):** D̄₃ ≥ 5·SE ve D̄₃ ≥ 0.05 ⇒ VAR; D̄₃ ≤ 2·SE ⇒ YOK;
   arası ZAYIF.
 - **H-200B-1b (manşet: çarpıklığın işareti):** havuzlanmış düzeltilmiş κ₃
   (k₃ − Δκ₃^{CUE}(ε)) > 0 ve ≥ 3σ ⇒ DÖNDÜ; < 0 ve ≥ 3σ ⇒ DÖNMEDİ; arası BELİRSİZ.
 - **H-200B-2 (200-A'nın önceden kayıtlı taşıması, a_k):** |D̄₃ − 0.2337| ≤ 0.1422 ⇒ κ₃'te
   a_k TUTAR; |D̄₂ + 0.0881| ≤ 0.0238 ⇒ κ₂'de a_k TUTAR. İkisi ⇒ a_k TUTAR; biri ⇒ KISMİ;
   hiçbiri ⇒ ÖLÜ. (Toleranslar 200-A kuralıyla HUKUM_200A.json'da sabitlendi.)
 - **Öncelik (çelişki kuralı):** H_geo (+0.0794) a_k'nın tolerans bandının [0.0915, 0.3759]
   hemen dışında; D̄₃ ≈ 0.085–0.09 çıkarsa H-200B-2 "a_k TUTAR" derken H-200B-3 GEOMETRİK
   diyebilir. Kural: **merdivenin biçimi ve manşet H-200B-3'e göre verilir**; H-200B-2 A'nın
   önceden kayıtlı taşıma sınavı olarak AYNEN raporlanır, çelişki varsa rapor bunu açıkça yazar
   ("a_k toleransı içinde, ama biçim GEOMETRİK"). A'nın toleransı a_k'nın kendi b0/b1 artıklarından
   gelen kaba bir bütçedir; H-200B-3'ün χ²'si b=2'nin gerçek SE'sini kullanır.
 - **H-200B-3 (merdivenin biçimi):** sınıflar SABİT (+0.2337), GEOMETRİK (+0.0795), SÖNMÜŞ
   (H_lin +0.0170 ∪ H_RMT 0; sınıfın χ²'si ikisinin küçüğü). χ²_h = (D̄₃ − tahmin_h)² /
   (SE² + σ_h² + 0.01²), σ_h tahminin A'dan gelen yayılımı (sabit için 0). En küçük χ²'li sınıf;
   bir sonrakine Δχ² ≥ 9 ⇒ KESİN, aksi hâlde BELİRSİZ; en iyi sınıfın χ²'si > 9 ⇒ HİÇBİRİ
   (merdiven başka bir biçimde).

**İKİNCİL (KAYIT):** D̄₂'nin hipotezlere uzaklığı; pencere başına D₃ (L eğilimi — 198'de bu
L aralığında gerçek sonlu-yükseklik kayması görüldü, ⇒ havuzlanmış D̄'nin yükseklik sistematiği
taşıyabileceği raporda tartışılır); ε = 0.1 ve 0.3'te aynı çözümleme; bağımsızlık:
Pearson(δ̃_n, x_n) olaylar içinde, CUE MC değerine göre (teorem 1(b): ε → 0'da 0; ε=0.2'de
CUE ≈ 0.01); aralık yasası: ε/2 altındaki olay oranı (U^{1/3} ⇒ 1/8). Yüksek-L pencereleri
bu kalemden ÇIKARILDI (eski motor C0-only; ileride yeniden doğrulanarak ayrı kalem).
Nicelik-çarpıklığı gibi dayanıklı bir ortak birincil ÖNCEDEN kaydedilemiyor: a_k biçiminin
sonlu L'de olasılık yasası yok (§6b), yalnız kümülantları tanımlı; bu yüzden birincil κ₃
kalır ve kuyruk duyarlılığı M8b + jackknife ile izlenir.

## Görülmüşlük beyanı

 - Yakın çiftlerde M_n, Z″ ya da log(M/δ̃²) programda HİÇ hesaplanmadı (25 Eyl envanteri;
   iskelet ve 200-A körlük bariyeri M9: sıfır başına Z′ saklanmadı). K0'da grep ile
   doğrulanır.
 - 55_win ve 53 dosyalarında M bütün aralıklar için mevcut ama küçük-aralık koşullu hiçbir
   istatistik hesaplanmadı (Not 1 bütün aralıkları kullandı: marjinal ortalama ve korelasyon).
 - Sayımlar (yukarıdaki tablo) yalnız sıfır konumlarından.

## Ölçüm tanımı

 - Motor: `200a_motor.py` (C0+C1; 200-A M1'de |ΔZ| ≤ 2.4·10⁻⁷). Tepe: aralık içinde 33 noktalı
   ızgara + parabolik rafinasyon (M2b ile doğrulanır). [Sapma: 129 noktalı ızgara, yukarıya bakın.]
 - Olay başına (t = m_n, δ̃_n, M_n) ölçümden SONRA saklanabilir (Not 8 figürleri için);
   ölçümden önce gerçek veride hiçbir M hesaplanmaz.
 - Sonlu-ε düzeltmesi Δκ_r^{CUE}(ε, L̄_W): N ∈ {9, 10, 11, 12} CUE MC (her N'de ε=0.2 için
   ≥ 4·10⁵ olay), L̄_W'ye doğrusal ara değer.

## Kapılar (K0 — ölçümden ÖNCE)

 - **M1b (tepe motoru):** pencere başına 300 rastgele yakın çiftte M_n: motor+ızgara vs mpmath
   (siegelz, dps ≥ 20, altın-oran maks. arama). Kapı: |ΔM|/M ≤ 10⁻⁴; iki motorla hesaplanan
   k₂/k₃ FARKI ≤ 0.002 / 0.005 (değerler raporlanmaz).
 - **M2b (ızgara çözünürlüğü):** 33 vs 129 nokta: |ΔM|/M ≤ 10⁻⁶ (yalnız fark).
 - **M6b (sentetik):** (i) Δκ_r^{CUE}(ε) tablosu; (ii) f-kalibrasyonu: gerçek olay sayılarıyla
   sentetik örnekler (eğik CUE ürün yasası, ve sağa çarpık bir örnekleyici: + bağımsız Euler
   p=2,3) aynı blok-jackknife hattından; (iii) güç: her hipotez altında havuzlanmış D̄₃ ve κ₃
   için Gauss yaklaşımıyla karar kuralları ⇒ karışıklık matrisi; P(doğru sınıf) < 0.80 olan
   çiftler raporlanır (kural gereği birleştirilir).
 - **M8b (blok sayısı):** 32 / 64 / 128 t-bloğu (W1'de 64 blokta ~18 olay/blok). Herhangi bir
   pencerede SE_{32} ile SE_{64} farkı > %20 ise o pencerede İKİSİNİN BÜYÜĞÜ kullanılır
   (tutucu); 128 blok yalnız tanı.
 - **M7b (karışım):** pencere içi L yayılımının κ tahminlerine etkisi < 0.002 mi.
 - **M9b:** makine gerçek veride M hesaplamadığını beyan eder; ölçüm betiği yalnız ONKAYIT
   push'landıktan sonra çalışır.

## Süreç

KALEM → Sonnet türetim teftişi → makine ajanı (M1b, M2b, M6b, M7b) → eşikler donar →
ONKAYIT_200B (sha256) push → ölçüm → ortak teftiş → mühür.
