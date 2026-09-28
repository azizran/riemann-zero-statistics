# KALEM — 200-A: EĞİM MERDİVENİNİN ALT BASAMAKLARI — sonlu yükseklikte aritmetik modelin seçimi
(26 Eylül 2026 — kullanıcı onaylı İKİ AŞAMALI tasarım: 200-A (bu kalem) alt basamaklarda
(b=0, b=1) sonlu-L modelini SEÇER ve mühürler; 200-B yakın çift basamağını (b=2) seçilen
modelle KÖR sınar. TÜRETİM TEFTİŞİNDEN GEÇTİ (Sonnet, 26 Eyl): model tahminleri bağımsız
olarak doğrulandı; altı tasarım düzeltmesi bu sürümde işlendi: (1) χ² kovaryanslı (k₂–k₃
korelasyonu ≈ −0.77, ürün yasasıyla ölçüldü); (2) önceden ilan edilmiş teori payı σ_teori
(yoksa mutlak uyum kuralı her modeli reddederdi); (3) SE sezgileri jackknife ile (dar formül
~2.7× iyimserdi); (4) 200-B körlüğü teknik olarak: sıfır başına Z′ SAKLANMAZ; (5) M_ak ve
M_hyb'nin asal başına formüllerinin farkı açık yazıldı; (6) Var_jk/Var_iid tanı raporu,
200-B toleransı kümülant başına. Eşikler M6'dan sonra donar.
MAKİNE RAPORU (ölçümden önce; gerçek veride hiçbir kümülant hesaplanmadı; `200_configs/
200a_MAKINE_RAPORU.md`): Motor `200a_motor.py` (vektörize RS float64). MÜHÜR ÖNCESİ SAPMA:
yalnız C0 kalanıyla M1 W1'de KALDI (|ΔZ| 6.2·10⁻⁵ > 2·10⁻⁵, düşük-t ucu, t^{−3/4}) ⇒ C1
terimi eklendi (C1 = −C0‴/96π²; mpmath dps 40, 28. derece Chebyshev, hata ≤ 2·10⁻¹⁵) ⇒ M1
GEÇTİ: |ΔZ| max W1 2.4·10⁻⁷, W2 2.0·10⁻⁸, W3 1.9·10⁻⁸; |ΔZ′| ≤ 8·10⁻⁸; iki motorun k₂/k₃
FARKI ≤ 3·10⁻⁶ / 4·10⁻⁵. M2 GEÇTİ (%100; en kötü |Z|/|Z′| 1.9·10⁻⁷). M7: bütün düzeltmeler
< 0.002 ⇒ ihmal. M6: f = 0.82–1.22 (12 istatistik; CUE ve hyb3 örnekleyicilerinden elemanca
maksimum); seçim gücü: gürültü-yalnız P(M* = doğru) ≥ 0.9986, gürültü + model hatası ≥ 0.857
(en kötü hyb5), P < 0.80 çift YOK ⇒ S*'ta önceden birleştirme gerekmedi. M9: yalnız blok güç
toplamları (3×2×128×5) yazılıyor, sahte uçtan uca koşuda assert geçti. Analiz betiğinde
Gauss için Δ ölçütü KALEM tanımına çekildi (χ²_Gauss − en iyi aritmetik). EŞİKLER DONDU.)

## Durum

Teorem iskeleti (`TEOREM_KUCUK_ARALIK_ISKELET_26EYL2026.md`, 38d8efa + a9f8fe0; Sonnet
teftişinden geçti): CUE'de iki özdeğer birleşirken geri kalanlar |Λ|⁴ ile eğilir;
koşullanan sıfır sayısı b = 0, 1, 2 ⇒ |Λ|^{2b} eğimi ("eğim merdiveni"). Zeta sanısı:
κ_r(zeta, b) = κ_r(model, b), aritmetik katkı her basamakta aynı.

Sorun (§6b): sonlu L'de "aritmetik model" tek değil.
 - a_k biçimi (κ^{arit}: κ₂ −0.088124, κ₃ +0.233653) L≈11'de olasılık yasası değil;
 - hibrit(X) (rastgele Euler çarpımı p ≤ X ⊗ eğik CUE(N_X), N_X = L/(e^γ log X)) gerçek bir
   olasılık modeli ama X'e bağlı;
 - b=2'de κ₃'ün İŞARETİ bütün aritmetik modellerde pozitif (sağlam), κ₂ ise 0.19–0.72
   arasında modele bağlı.
Karar: modeli b=2'ye bakmadan, bilinen bölgede (b=0, b=1) seç; b=2'ye kör taşı.

## Soru

L ≈ 9.3–11.7'de (i) rastgele t'de log|Z(t)| (b=0) ve (ii) sıfırlarda log|Z′(γ_n)| (b=1)
dağılımlarının κ₂ ve κ₃'ünü hangi sonlu-L modeli betimliyor?

## Model ailesi F (hepsi PARAMETRESİZ; N = L)

 - **M_CUE:** κ_r = κ_r^{CUE,b}(L) — eğik CUE, eğim |Λ|^{2b}, kalan n = L − b.
 - **M_ak:** κ_r^{CUE,b}(L) + κ_r^{arit} (Keating–Snaith / HKO biçimi).
 - **M_hyb(X), X ∈ {2, 3, 5, 7, 11}:** κ_r^{CUE,b}(N_X) + Σ_{p≤X} κ_r(−log|1 − p^{−1/2}e^{iα_p}|),
   N_X = L/(e^γ log X); asal başına κ₂ = ½Li₂(1/p), κ₃ = (3/2)S_{1,2}(1/p).
   **DİKKAT (yazım tuzağı):** M_ak'nın asal başına κ₂ katkısı ½[log(1−1/p) + Li₂(1/p)]'dir
   (log terimi CUE(L)'ye göre normalizasyondan gelir ve toplamı yakınsatır); M_hyb'de log
   terimi YOKTUR, çünkü CUE kısmı N_X ile küçültülmüştür (sonlu X toplamı). İki formül farklı
   nesnelerdir; birbirine "düzeltilmemelidir".
 - **M_Gauss** (Not 1'in Null 1'i, aynı spektrumlu Gauss yüzeyi): b=0 log|N(0,σ²)|:
   κ₂ = π²/8 = 1.2337, κ₃ = −7ζ(3)/4 = −2.1036; b=1 (Kac–Rice, |X′| ile eğik ⇒ χ₂):
   κ₂ = π²/24 = 0.4112, κ₃ = −ζ(3)/4 = −0.3005. L'den bağımsız.

κ_r^{CUE,b}(N) = Σ_{j=1}^{N−b} [ψ^{(r−1)}(j+2b) − 2^{1−r}ψ^{(r−1)}(j+b)], tam-sayı olmayan N
için analitik devam (`200t_merdiven.py`). Tahminler `200_configs/200a_tahmin.json`
(`200a_pencereler_tahmin.py`).

## Pencereler (yalnız sıfır konumlarından; zeros6)

| pencere | L aralığı | sıfır indeksi | n sıfır | t aralığı | L̄ (b=0, t-düzgün) | L̄ (b=1, sıfır ağırlıklı) |
|---|---|---|---|---|---|---|
| W1 | [8, 10) | 20 868 – 198 239 | 177 371 | 18 731 – 138 396 | 9.313 | 9.343 |
| W2 | [10, 11) | 198 239 – 598 742 | 400 503 | 138 396 – 376 200 | 10.582 | 10.589 |
| W3 | [11, 12.10] | 598 742 – 2 001 052 | 1 402 310 | 376 200 – 1 132 491 | 11.650 | 11.658 |

Aynı pencereler 200-B'de b=2 için kullanılacak. Karışım düzeltmesi (M7) teftişte W1 için
kestirildi: b=0 k₂ −0.0008, b=1 k₂ −0.0005 (κ₃ +0.0003/+0.0002) ⇒ 0.002 eşiğinin altında.

## Tahmin tablosu (κ₂ / κ₃)

| pencere | b | M_CUE | M_ak | hyb2 | hyb3 | hyb5 | hyb7 | hyb11 | Gauss |
|---|---|---|---|---|---|---|---|---|---|
| W1 | 0 | 1.905/−2.387 | 1.817/−2.153 | 2.091/−2.226 | 2.045/−2.115 | 1.961/−2.027 | 1.942/−1.973 | 1.887/−1.909 | 1.234/−2.104 |
| W1 | 1 | 0.540/−0.190 | 0.452/+0.043 | 0.749/−0.035 | 0.768/+0.056 | 0.760/+0.110 | 0.787/+0.138 | 0.789/+0.163 | 0.411/−0.300 |
| W2 | 0 | 1.969/−2.397 | 1.880/−2.163 | 2.155/−2.238 | 2.108/−2.134 | 2.024/−2.053 | 2.005/−2.004 | 1.949/−1.947 | 1.234/−2.104 |
| W2 | 1 | 0.592/−0.198 | 0.503/+0.036 | 0.797/−0.043 | 0.810/+0.045 | 0.795/+0.098 | 0.817/+0.125 | 0.815/+0.150 | 0.411/−0.300 |
| W3 | 0 | 2.017/−2.403 | 1.929/−2.170 | 2.203/−2.246 | 2.156/−2.146 | 2.072/−2.071 | 2.052/−2.026 | 1.996/−1.973 | 1.234/−2.104 |
| W3 | 1 | 0.632/−0.203 | 0.544/+0.031 | 0.836/−0.049 | 0.844/+0.037 | 0.823/+0.088 | 0.843/+0.115 | 0.836/+0.140 | 0.411/−0.300 |

Başvuru için (200-B'de mühürlenecek, burada karar yok) b=2: W2'de M_CUE 0.318/−0.071,
M_ak 0.230/+0.162, hyb2 0.540/+0.080, hyb3 0.598/+0.157, hyb5 0.630/+0.194, hyb7
0.680/+0.211, hyb11 0.708/+0.221, Gauss 0.234/−0.104.

## Görülmüşlük beyanı

 - b=0: `22_selberg_test.py` (programın ilk günleri, düşük yükseklik) log|Z(t)|'yi Selberg
   normaline karşı KS testiyle sınadı. κ₂/κ₃ bu model ailesine karşı HİÇ hesaplanmadı;
   W1–W3 yüksekliklerinde rastgele-t değer dağılımı hiç hesaplanmadı (K0'da grep ile
   doğrulanır).
 - b=1: Z′ programda HİÇ hesaplanmadı (25 Eyl envanteri §3.4; K0'da grep).
 - **200-B'nin körlüğü:** b=1 örneği yakın çiftlere ait sıfırları da içerir (δ̃<0.2 için
   sıfırların ~%1.4'ü). 200-A'da **hiçbir istatistik δ̃'ye göre alt kümelenmez**, M_n ya da
   Z″ hesaplanmaz. Sıfır başına Z′ değerleri `200a_b1_*.npz`'ye yazılır ve dosya
   **Teknik bariyer:** sıfır başına Z′ değerleri SAKLANMAZ. Makine yalnız jackknife blokları
   başına güç toplamlarını (n, Σx, Σx², Σx³, Σx⁴) yazar (`200a_bloklar.npz`); blok başına
   binlerce sıfır olduğundan yakın çift bilgisi geri çıkarılamaz. Ara bellekteki değerler
   betik sonunda atılır (M9).

## Ölçüm tanımı

 - **Motor:** Motor B (vektörize Riemann–Siegel, float64; `41_bigT_scan.py` /
   `55_kucuk_tau_scan.py` ailesi) + ana toplamın terim terim analitik türevi ile Z′
   (Z′ = −2Σ n^{−1/2}(θ′(t) − log n) sin(θ(t) − t log n) + kalan terimin türevi; kalanın
   türevi M1'de ölçülüp ya dahil edilir ya da ihmal gerekçelendirilir).
 - **b=0:** her pencerede n₀ = 10⁶ nokta, t ~ Düzgün[t_a, t_b] (tohum 200), x = log|Z(t)|.
 - **b=1:** penceredeki bütün sıfırlar, x = log(|Z′(γ_n)|·2π/L(γ_n)) (açılmış; κ₂, κ₃
   ölçekten bağımsız).
 - **İstatistikler:** k-istatistikleri k₂, k₃ (yansız). Hata: pencere başına ORTAK t-blokları
   (64 bitişik blok, sınırlar t'de; b=0 örnekleri ve b=1 sıfırları aynı bloklara düşer) ile
   bir-blok-dışarıda jackknife ⇒ pencere başına 4×4 kovaryans C_W (k₂^{b0}, k₃^{b0}, k₂^{b1},
   k₃^{b1}); pencereler bağımsız ⇒ 12×12 blok-köşegen C. M6 kalibrasyonu: C_ij → f_i f_j C_ij.
   Tanı: Var_jk / Var_iid oranı her istatistik için raporlanır (≫1 beklenir; ≈1 hata işareti).
 - **Teori payı (önceden ilan):** σ_teori,r = 0.1·|κ_r^{arit}| + 0.01 ⇒ κ₂: 0.0188,
   κ₃: 0.0334 (her basamak ve pencerede). Gerekçe: teftişin Landau–Gonek faz yanlılığı
   kestirimi (aritmetik sinyalin %5–15'i; p = 2, 3, 5 için asal başına Δκ₂ ≈ +0.005…+0.0015,
   Δκ₃ ≈ −0.005…−0.002, aynı işaretli) + O(1/L²) ve kesikli-X ızgarası için 0.01. Pencereler
   arası ortak bir sistematik olduğu için σ_teori tam korelasyonlu değil, BAĞIMSIZ eklenir
   (tutucu olmayan seçim; ikincil olarak tam korelasyonlu hâli de raporlanır).
 - **Karışım (M7):** tahmin = örneklerin L_i'si üzerinden model kümülantlarının ortalaması
   + karışım terimi (κ₁(L)'nin pencere içi varyansı vb.). |düzeltme| < 0.002 ise ihmal,
   değilse uygulanır.

## Hipotezler ve karar kuralları

Model χ²'si (12 bileşenli artık vektörü r_m = k^{göz} − κ^{m}; 3 pencere × 2 basamak ×
{k₂, k₃}):
 χ²_m = r_mᵀ (C_f + diag(σ_teori²))⁻¹ r_m,  C_f = f-kalibre jackknife kovaryansı.
Not: teftiş k₂ yalnız başına komşu hibritleri (ör. hyb7–hyb11, Δκ₂ ≈ 0.002) HİÇBİR
basamakta ayıramaz; ayrımı κ₃ yapar.

 - **H-200A-1 (alt basamaklarda aritmetik imza):** Δ₁ := χ²_{M_CUE} − min_{m∈arit} χ²_m.
   Δ₁ ≥ 25 ⇒ TUTAR; Δ₁ < 4 ⇒ ÖLÜ; arası BELİRSİZ. (Ön beklenti TUTAR: HKO 2000 ve
   Keating–Snaith b=0/b=1'de aritmetik çarpanı öngörüyor; bu bilinen bölgenin sınavı.)
   Aynı ölçüt M_Gauss için de raporlanır.
 - **H-200A-2 (model seçimi):** M* = argmin_{F} χ²_m. Δχ²(ikinci en iyi) ≥ 9 ⇒ KESİN; aksi
   hâlde S* = {m : χ²_m − χ²_min < 9} kümesi 200-B'ye bant olarak taşınır. "KESİN" yalnız
   "ızgaradaki en yakın üye" anlamındadır (gerçek iki X arasında olabilir; ikincil sürekli-X
   uyumu bunu raporlar).
 - **Uyum kalitesi:** G := χ²_min / 12.
   G ≤ 3 İYİ; 3 < G ≤ 10 KABA (200-B'de κ₂ sınavı KAYIT'a iner, κ₃ tahmini S* bandı ±
   A'daki κ₃ artıklarının RMS'i ile taşınır); G > 10 AİLE YETERSİZ (200-B yalnız modelden
   bağımsız D₃ işaret sınavını yapar). (G teori payı dahil hesaplanır; payın kendisi M6'dan
   ÖNCE sabittir ve veriye göre ayarlanmaz.)
 - **İKİNCİL (KAYIT):** basamak başına ayrı seçim (M*_0 = M*_1 mi? merdiven sanısının ilk
   işareti); sürekli X uyumu; model başına ortak N_eff = L + c uyumu (Not 1 ile karşılaştırma);
   artıkların W1→W3 eğilimi.

## 200-B'ye taşınan (200-A mühürlenince)

 - M* (ya da S*) ve onun b=2 tahminleri (pencere başına κ₂, κ₃).
 - Tolerans kuralı (A'nın sonucu görülmeden önce burada sabit): KÜMÜLANT BAŞINA
   tolerans_r = max(σ_teori,r, 2 × A'daki M* artıklarının aynı r için RMS'i) (b=0 ve b=1
   artıkları birlikte, pencereler birlikte; r = 2 ve r = 3 ayrı).
 - 200-B'nin BİRİNCİL sınavı model seçiminden bağımsız kalır: D₃ := κ₃(b=2) − κ₃^{CUE,2}(L);
   aritmetik aile [+0.15, +0.28], M_CUE 0, Gauss ≈ −0.03.

## Kapılar (K0 — ölçümden ÖNCE, hepsi zorunlu)

 - **M1 (motor doğruluğu):** her pencerede 2000 rastgele t ve 2000 rastgele sıfırda Motor B
   ile mpmath (`siegelz(t)`, `siegelz(t, derivative=1)`, dps ≥ 20) karşılaştırması:
   |ΔZ| ≤ 2·10⁻⁵, |ΔZ′| ≤ 2·10⁻⁴; bu alt örneklerde iki motorla hesaplanan k₂, k₃ farkı
   ≤ 0.002 / 0.005.
 - **M2 (sıfır sağlaması):** her pencerede 2000 sıfırda |Z(γ_n)| ≤ 10⁻⁴·|Z′(γ_n)| ⇒ tablo
   ile motor uyumlu.
 - **M6 (sentetik güç + SE kalibrasyonu):** F'nin her modelinden, pencerelerin gerçek n'siyle
   sentetik örnekler (CUE kısımları çarpım yasasıyla: b=0 γ_j, b=1 |Λ|²-eğik ξ_j; tam-sayı
   olmayan N için komşu tam sayıların karışımı; Euler çarpanları bağımsız düzgün fazlarla) ⇒
   seçim hattı 200 tekrar: doğru modeli seçme olasılığı ≥ 0.80 olmalı; olmayan çiftler
   önceden S*'ta birleştirilir (kural). f_b = tekrarlar arası SD / ortalama SE_jk.
 - **M7 (karışım):** her pencere/basamak için düzeltme hesaplanır ve raporlanır.
 - **M8 (etkin örnek):** blok boyu duyarlılığı (32/64/128 blok); SE değişimi > %20 ise
   128 blok kullanılır. Teftiş: b=0 örnek aralığı açılmış birimde W1 ≈ 0.18, W2 ≈ 0.40,
   W3 ≈ 1.4 (yerel aşırı örnekleme); 64 blok ≈ 2 800–21 900 açılmış birim ⇒ yeterli.
 - **M9 (körlük bariyeri):** makine betiği sıfır başına Z′ dizisini diske yazmadığını ve
   yalnız blok güç toplamlarını sakladığını assert eder; K0 raporu çıktı dosyalarının
   listesini ve boyutlarını verir.

## Süreç

KALEM → Sonnet türetim teftişi → makine ajanı (Motor B + Z′; M1, M2, M7, M8; M6) → eşikler
donar → ONKAYIT_200A.json (sha256) push → ölçüm → ortak teftiş → mühür → KALEM 200-B
(b=2 tahminleri M*'dan) → teftiş → makine → ONKAYIT_200B push → yakın çift ölçümü.
