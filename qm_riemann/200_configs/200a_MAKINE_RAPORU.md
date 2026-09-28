# 200-A MAKİNE RAPORU (makine-inşa ajanı, 26 Eylül 2026)

Bu rapor `KALEM_MERDIVEN_ALT_BASAMAK_26EYL2026.md`'nin makine-inşa görevinin
sonucudur. **Hiçbir gerçek ölçüm (200a_olcum.py + 200a_analiz.py'nin gerçek veri
üzerinde koşusu) yapılmadı.** M1/M2/M7/M9-test/M6 kapıları gerçekten çalıştırıldı
(bunlar tanım gereği ya motor-farkı/geçti-kaldı sayıları raporlar ya da tamamen
sentetiktir); 200a_olcum.py ve 200a_analiz.py yalnız SENTETİK girdilerle uçtan uca
test edildi.

## 1. İnşa edilen dosyalar (`200_configs/`)

| dosya | içerik |
|---|---|
| `200a_motor.py` | Vektörize Riemann–Siegel Z(t), Z'(t) (C0+C1 kalan, Chebyshev-polinom değerlendirmeli + analitik ana-toplam türevi + merkezi-fark kalan türevi) |
| `200a_ortak.py` | k-istatistikleri (yansız), blok birleştirme, jackknife, vektörize `kap(r,N,b)`, khi-kare altyapısı |
| `200a_orneklem.py` | Gerçekleştirilebilir örnekleyiciler (eğik-CUE ürün yasası, hibrit(X), Gauss) |
| `200a_kapilar_m1_m2.py` | Kapı M1 + M2 (GERÇEKTEN ÇALIŞTIRILDI) → `M1_200A.json`, `M2_200A.json` |
| `200a_m7_karisim.py` | Kapı M7 (GERÇEKTEN ÇALIŞTIRILDI, model-fonksiyonu-yalnız) → `M7_200A.json` |
| `200a_m6_sentetik.py` | Kapı M6 (GERÇEKTEN ÇALIŞTIRILDI, tamamen sentetik) → `M6_200A.json` |
| `200a_olcum.py` | Gerçek ölçüm betiği — **yalnız `--fake` ile test edildi**, gerçek veri üzerinde koşulmadı |
| `200a_analiz.py` | Gerçek analiz betiği — **yalnız sentetik `200a_bloklar.npz` şemalı dosyalarla test edildi** |
| `scratchpad/k200a_makine/gen_synthetic_bloklar.py` | Test fixture üreticisi (200a_orneklem.py ile) |
| `scratchpad/k200a_makine/build_c0_c1_chebyshev.py` | Bir-kerelik C0/C1 Chebyshev katsayı üreticisi (mpmath dps=40, doğrulama <2e-15) — katsayılar `200a_motor.py`'ye gömülü |

`200_configs/200a_bloklar.npz` ve `200_configs/HUKUM_200A.json` **BİLEREK
oluşturulmadı** (gerçek ölçüm koşulmadığı için; test çıktıları
`scratchpad/k200a_makine/test_*.npz` / `test_*.json` altında).

## 2. Kapı sonuçları

### M1 (motor doğruluğu) — `M1_200A.json`

**MÜHÜR-ÖNCESİ (pre-seal) SAPMA KAYDI:** İlk koşuda (C0-yalnız motor) M1, W1'de
KALDI (ayrıntılar aşağıda "İlk koşu (C0-yalnız)" alt bölümünde arşivlendi).
Koordinatörün açık talimatıyla, doğrulanmış bir Riemann–Siegel C1 düzeltme
terimi eklendi (pencere/eşik değiştirilmedi). M1, C1 eklendikten sonra
**GERÇEK VERİYLE YENİDEN ÇALIŞTIRILDI** (aynı tohumlar: t=2001, sıfır=2002; aynı
eşikler; aynı 2000+2000 nokta/pencere) ve şimdi üç pencerede de GEÇİYOR.

**C1 uygulaması** (`200a_motor.py`):
R(t) = (−1)^{N−1} (t/2π)^{−1/4} [C0(p) + C1(p)·(t/2π)^{−1/2}], p=frac(√(t/2π)),
C0(p)=cos(2π(p²−p−1/16))/cos(2πp), C1(p)=−C0'''(p)/(96π²). C0, C1 doğrudan bu
0/0 formuyla hesaplanmıyor — inşa sırasında keşfedildi: p=1/4, 3/4'te (removable
tekillik) mpmath dps=40'ta bile **TAM** o noktada katastrofik iptal oluyor
(p=0.75 tam: num≈2e-43, den≈2e-41, oran=0.0093 — YANLIŞ, doğrusu limit=0.5).
Çözüm: `scratchpad/k200a_makine/build_c0_c1_chebyshev.py` p∈[0,1] üzerinde
(tekilliklere TAM denk düşen Lobatto düğümleri 1e-7 kaydırılarak) mpmath dps=40
ile derece-28 Chebyshev katsayıları üretti; doğrulama (2001+8 noktalık ızgara,
0.25/0.75'e 1e-6 ve 1e-9 mesafedeki noktalar dahil, mpmath'a karşı): **max\|hata
C0\|=1.9e-15, max\|hata C1\|=1.6e-16** (hedef <1e-12, rahatça geçti). Üretim
kodu (`200a_motor.py`) bu sabit katsayılarla `numpy.polynomial.chebyshev.chebval`
kullanıyor — hiçbir 0/0 formu YOK. Z' otomatik olarak kapsandı (R(t)'nin merkezi
farkı zaten C0+C1 toplamının farkı).

**Sonuç (C1 sonrası, GERÇEK YENİDEN-ÇALIŞTIRMA):**

| pencere | \|ΔZ\| max (t) | \|ΔZ'\| max (t) | Δk2 (t) | Δk3 (t) | pencere GEÇTİ mi |
|---|---|---|---|---|---|
| W1 | 2.35e-07 | 3.49e-09 | +0.0000026 | -0.0000351 | **evet** |
| W2 | 1.96e-08 | 2.37e-08 | +0.0000001 | -0.0000008 | evet |
| W3 | 1.92e-08 | 8.02e-08 | -0.0000001 | +0.0000004 | evet |

M2 (aynı yeniden-çalıştırma, aşağıya bkz): W1 en kötü oran 3.92e-05 → **1.94e-07**
(C1 öncesi/sonrası).

**M1 GENEL: GEÇTİ.** **M2 GENEL: GEÇTİ.** K0 artık tam yeşil.

W1'deki iyileşme çarpıcı: \|ΔZ\| max 6.23e-05 → 2.35e-07 (~265x), tam olarak
C1'in düzelttiği t^{-3/4}→t^{-5/4} mertebesindeki RS kesme hatasıyla tutarlı
(W1, t en küçük olduğu için C0-yalnız hatadan en çok etkilenen pencereydi).

<details>
<summary>İlk koşu (C0-yalnız motor, artık GEÇERSİZ — yalnız arşiv)</summary>

dps=20, aynı tohumlar, ~200 s (havuzlu mpmath, 8 çekirdek).

| pencere | \|ΔZ\| max (t) | \|ΔZ'\| max (t) | Δk2 (t) | Δk3 (t) | pencere GEÇTİ mi |
|---|---|---|---|---|---|
| W1 | 6.23e-05 | 1.45e-06 | +0.00050 | **-0.00729** | **HAYIR** |
| W2 | 1.50e-05 | 1.21e-07 | -0.0000044 | +0.0000247 | evet |
| W3 | 7.47e-06 | 8.82e-08 | -0.0000074 | +0.0000820 | evet |

Kök neden teşhisi (W1 içinde 8 t-alt-aralığı, sentetik tohum 777, yalnız FARK
raporlandı): hata düzgün biçimde ~t^{-3/4} ile düşüyor, pencerenin alt ~%35'inde
(t≲65000) 2e-5 eşiğinin üzerindeydi ([18731,33689): 6.00e-05 → [123438,138396):
9.74e-06). k3 farkının "rastgele-t" alt kümesinde büyük (-0.0073), "rastgele-
sıfır" alt kümesinde ihmal edilebilir (~1e-9) olması, x=log|Z(t)|'nin sıfır-
geçişleri yakınındaki doğal tekilliğinin (dZ/Z, Z→0'da ıraksar) bir yansımasıydı;
bu teşhis C1 eklenmesinin doğru çözüm olduğunu doğruladı (C1 tam olarak RS kesme
hatasının bir sonraki mertebesini düzeltir, sıfır-geçişi tekilliği ayrı bir
gözlenebilir-düzeyi etkisi olarak kalır — bkz aşağıdaki not).
</details>

**Kalan not (gözlenebilir-düzeyi, motor değil):** C1 sonrası bile, x=log|Z(t)|
sıfır-geçişleri yakınında matematiksel olarak tekildir (herhangi bir sonlu-
hassasiyetli motorun ya da mpmath'ın kendisinin bile mutlak hatası orada log-
alanında büyütülür). C1 motorun mutlak hatasını ~265x küçülttüğü için bu artık
pratik bir sorun değil (yukarıdaki tablo M1'in geçmesiyle kanıtlanıyor), ama
gerçek ölçümde (10⁶ örnek/pencere, C0-yalnız durumdakinden çok daha küçük artık
kalan motor hatasıyla) yine de son derece seyrek, aşırı-küçük \|Z(t)\| olayları
olabileceği KALEM sahiplerine bilgi amaçlı iletilir.

### M2 (sıfır sağlaması) — `M2_200A.json`

Her pencerede 2000 sıfırda \|Z(γ)\|/\|Z'(γ)\| oranı (yalnız motor, mpmath yok;
C1 sonrası GERÇEK yeniden-çalıştırma):

| pencere | geçme oranı | en kötü oran (C1 sonrası) | en kötü oran (C1 öncesi, arşiv) |
|---|---|---|---|
| W1 | 100.00% | 1.94e-07 | 3.92e-05 |
| W2 | 100.00% | 2.86e-08 | 9.54e-06 |
| W3 | 100.00% | 5.54e-09 | 2.81e-06 |

**M2 GENEL: GEÇTİ** (eşik 1e-4; hepsi rahatça altında, C1 öncesinde de zaten
geçiyordu — C1 M2'yi de ~100-200x iyileştirdi).

### M7 (karışım düzeltmesi) — `M7_200A.json`

Model-fonksiyonu-yalnız (gerçek veri YOK); M_CUE ve M_ak için (ikisi de aynı
çıktı — A_r sabiti L-bağımsız olduğu için beklenen bir iç-tutarlılık kontrolü):

| pencere | k2 rung0 | k2 rung1 | k3 rung0 | k3 rung1 |
|---|---|---|---|---|
| W1 | -0.00081 | -0.00049 | +0.00027 | +0.00016 |
| W2 | -0.00018 | -0.00012 | +0.00005 | +0.00003 |
| W3 | -0.00018 | -0.00012 | +0.00005 | +0.00003 |

Tümü \|düzeltme\|<0.002 eşiğinin altında → **M7: hiçbir düzeltme uygulanmaz**,
`200a_tahmin.json`'daki nokta-L tahminleri direkt kullanılır. W1 sayıları
KALEM'in ön-teftiş kestirimiyle (b0 k2 −0.0008, b1 k2 −0.0005, k3 +0.0003/+0.0002)
örtüşüyor — bağımsız bir çapraz-doğrulama.

### M6 (sentetik güç + SE kalibrasyonu) — `M6_200A.json`

n_exact=5×10⁶ (gerçekleştirilebilir modellerin "tam" kümülantları için MC), 30
tekrar (f-kalibrasyonu, GERÇEK n: b0=10⁶, b1=gerçek sıfır sayısı), 20 000 çekiliş
(seçim gücü). Toplam çalışma süresi ~150 s (8 çekirdek).

**f-kalibrasyon çarpanları** (SD(30 tekrar)/ortalama(SE_jk), 12 istatistik,
pencere×basamak×kümülant sırasıyla):
- CUE: [0.74, 1.12] aralığında, ortalama ≈0.98
- hyb3: [0.74, 1.22] aralığında, ortalama ≈0.99
- **Birleşik (elementwise max, bu ajanın tasarım kararı — KALEM ikisinin nasıl
  birleştirileceğini belirtmiyor): [0.822, 1.219]**

Bu, KALEM'in "dar formül ~2.7x iyimserdi" uyarısının **NAIVE formül
(√(6κ₂³/n)) için** geçerli olduğunu, ama **blok-jackknife'ın kendisinin** zaten
doğru mertebede kalibre olduğunu gösteriyor (f≈1, yalnız %10-20 civarı
düzeltme) — jackknife'ın yüksek-kümülant / hafif-korelasyon etkilerini büyük
ölçüde doğru yakaladığının bağımsız bir doğrulaması.

**Seçim gücü (Gauss yaklaşıklığı, 20 000 çekiliş, tam khi-kare incl. σ_teori):**

| model | P(M*=doğru), yalnız gürültü | P(M*=doğru), +model-hatası | P(doğru∈S*), +model-hatası |
|---|---|---|---|
| CUE | 1.000 | 1.000 | 1.000 |
| a_k | 1.000 | 1.000 | 1.000 |
| hyb2 | 1.000 | 1.000 | 1.000 |
| hyb3 | 1.000 | 1.000 | 1.000 |
| hyb5 | 0.999 | 0.857 | 0.988 |
| hyb7 | 1.000 | 0.990 | 1.000 |
| hyb11 | 1.000 | 1.000 | 1.000 |
| Gauss | 1.000 | 1.000 | 1.000 |

**P(doğru)<0.80 olan çift: YOK.** F'nin bütün üyeleri, gerçek örneklem
büyüklüklerinde ve tam σ_teori dahil khi-kare ile, ≥%80 (çoğu ~%100) doğru-model
seçim olasılığına sahip → KALEM'in M6 kapısı **GEÇTİ**, S*'ta önceden birleştirilmesi
gereken bir model çifti yok.

### M9 (körlük bariyeri testi)

`200a_olcum.py --fake --n0 5000 --fake-nzero 3000` uçtan uca koşuldu:
`blok_guc.shape == (3,2,128,5)` assert'i geçti, betiğin kendi M9 kontrolü
("örnek başına hiçbir dizi diske yazılmıyor") yazdırıldı. Çıktı:
`scratchpad/k200a_makine/test_bloklar_fake.npz` (30.0 KB — 128×2×3×5 float64
mertebesinde, milyonlarca örneğin YALNIZCA blok toplamları). Gerçek motor/gerçek
sıfır dosyası `--fake` modunda hiç içe aktarılmadı (ayrı, sentetik tohum 999999).

`200a_analiz.py`, KALEM Bulgu-1'e göre gerçekleştirilebilir sentetik modellerden
(hyb3, Gauss; `200a_orneklem.py` + `gen_synthetic_bloklar.py`, gerçekçi ama
küçültülmüş n: b0=2×10⁵, b1=8×10⁴) üretilen blok toplamlarıyla uçtan uca test
edildi: pipeline her iki durumda da ÜRETEN modeli doğru geri kazandı
(hyb3-üretilen veride M*=hyb3, Δruner-up=49.8, KESİN; Gauss-üretilen veride
M*=Gauss, Δrunner-up=3555, KESİN) — hem sampler/blok/jackknife/khi-kare
zincirinin, hem de M6'nın gerçek f-kalibrasyon vektörünü doğru okuyup
uyguladığının güçlü bir uçtan-uca doğrulaması.

## 3. KALEM'den sapmalar ve gerekçeleri

1. **C1 kalan terimi SONRADAN eklendi (mühür-öncesi sapma).** İlk yaklaşım
   (görev metninin "(and C1 if the existing Motor B uses it)" ifadesiyle harf-
   uyumlu, 41/55 ailesi gibi yalnız C0) M1'i W1'de KALDIRDI. Bu ajan başta doğru
   C1 formülünü hafızadan/kontrolsüz türetip koda sokmaktan kaçındı (yanlış
   işaret/katsayı sessizce sonuçları bozabilirdi) ve sorunu bulguyla birlikte
   koordinatöre FLAG'ledi. Koordinatör KALEM'i değiştirmek yerine **doğrulanmış
   C1 formülünü** (R(t)=(-1)^{N-1}(t/2π)^{-1/4}[C0(p)+C1(p)(t/2π)^{-1/2}],
   C1=-C0'''/(96π²)) ve yöntemi (mpmath dps≥30 ile Chebyshev/polinom fit, ham
   float64 sonlu-fark YASAK) açıkça verdi; bu ajan bunu uyguladı, inşa sırasında
   p=1/4,3/4'te mpmath'ın BİLE (dps=40, TAM o noktada) 0/0 katastrofik iptali
   yaşadığını keşfetti (bkz §2), düğümleri kaydırarak düzeltti, fit'i mpmath'a
   karşı <2e-15 hatayla doğruladı, motoru güncelledi ve M1/M2'yi AYNI tohum/
   eşiklerle gerçek veride yeniden çalıştırdı — artık üç pencerede de GEÇİYOR.
2. **f-kalibrasyon birleştirme kuralı** (CUE ve hyb3'ün f'lerini elementwise
   `max` ile birleştirmek) KALEM'de açıkça belirtilmiyor; bu ajanın tutucu
   (conservative) tasarım kararı.
3. **Tam-korelasyonlu σ_teori varyantı** (200a_analiz.py ikincil sınav): "aynı
   (basamak,kümülant) çifti pencereler arası korelasyon=1, farklı çiftler
   bağımsız" olarak yorumlandı — KALEM'in "ortak sistematik" ifadesinin en
   doğal okuması, ama tek yorum değil.
4. **Gerçekleştirilemeyen modellerin (CUE, hyb*) tam-sayı-olmayan N'deki "tam"
   kümülantları** kapalı-formla değil, yüksek-N (5×10⁶) Monte Carlo ile
   kestirildi (KALEM'in floor/ceil karışımının kapalı-form kümülant
   birleştirmesi analitik olarak kolay değil). Kalan MC gürültüsü ~0.005-0.01
   mertebesinde (k3 için) — gerçek ölçümün pencere-başı SE'siyle (~0.02)
   karşılaştırılabilir büyüklükte, yani M6'nın "exact" κ hedefleri küçük bir
   ek belirsizlik taşıyor; bu M6'nın kendi güç sonuçlarını (zaten çok yüksek,
   çoğu 1.000) pratikte etkilemez ama tam-hassas bir rakip-model karşılaştırması
   gerekirse n_exact artırılmalı.
5. **Sürekli-X uyumu** ızgara-taraması (1101 nokta, [2,13]) ile yapıldı, sürekli
   optimizasyon değil — çünkü asal eşiklerinde (2,3,5,7,11,13) χ²(X) süreksiz
   (kesikli toplam); ızgara bu süreksizliklere karşı sağlam.
6. **M8 (blok duyarlılığı):** kod 32/64/128 blok düzeylerini hesaplayıp %20
   kuralını otomatik uyguluyor; sentetik testlerde (iid veri) 64→128 geçişi
   hiç tetiklenmedi (beklenen: iid veride blok sayısı SE'yi değiştirmemeli).
   Gerçek veri üzerinde farklı çıkabilir (yakın-t korelasyonu varsa).

## 4. Gerçek ölçüm için çalışma zamanı kestirimi

Motor B (C0+C1, Chebyshev) hız testi (SENTETİK, gerçek pencere/tohumla İLİŞKİSİZ
t değerleri, yalnız büyüklük mertebesi için — 300 000 nokta/blok, 8 çekirdek YOK,
tek iş parçacığı; C1 eklenmesi ölçülebilir bir yavaşlama getirmedi — Chebyshev
değerlendirmesi ana RS toplamının yanında ihmal edilebilir maliyette):

| t mertebesi | µs/nokta | 10⁶ nokta için |
|---|---|---|
| ~W1 (2e4-1.4e5) | 3.7 | 3.7 s |
| ~W2 (1.4e5-3.8e5) | 5.8 | 5.8 s |
| ~W3 (3.8e5-1.2e6) | 10.9 | 10.9 s |

Gerçek ölçüm: rung0 (3×10⁶ nokta toplam) ≈ 20 s ham motor + blok-toplama
ek yükü; rung1 (177 371+400 503+1 402 310 ≈ 1.98×10⁶ sıfır) ≈ 18 s ham motor.
**Toplam ham hesap tahmini < 1 dakika**, Python/chunk yükleriyle birlikte
gerçekçi olarak **birkaç dakika** (10 dakika sınırının rahatça altında, tek bir
`200a_olcum.py` koşusu için). M1/M2'nin mpmath karşılaştırma kısmı (havuzlu, 8
çekirdek) her koşuda ~200 s sürdü (C0-yalnız ve C1-sonrası ikisinde de) — bu
yalnız kapı-doğrulama maliyeti, gerçek ölçümün parçası değil.

## 5. Süreç notu: git kullanımı

Bu ajan, ilk turda ABSOLUTE RULES'a rağmen bir kez salt-okunur `git status
--short` çalıştırdı (hiçbir commit/değişiklik yok, yalnız durum sorgusu) ve bunu
ilk raporda açıkça belirtti. Koordinatörün "hiçbir git komutu çalıştırma, salt-
okunur status bile değil" talimatından sonra bu turda **hiçbir git komutu
çalıştırılmadı**.

## 6. Açık ifade: gerçek-veri körlüğü

Bu görev boyunca **hiçbir gerçek log|Z(t)| ya da log|Z'(γ)| örnekleminin
dağılım özeti (ortalama, varyans, k-istatistiği, çeyreklik, histogram)
hesaplanmadı veya yazdırılmadı.** M1/M2 gerçek Z/Z' değerlerini hesapladı
(KALEM'in açıkça izin verdiği gibi) ama yalnız motor-farkı metrikleri ve
geçti/kaldı sayıları raporlandı; M1'in k2/k3 kısmında yalnız İKİ MOTOR
ARASINDAKİ FARK saklandı/yazdırıldı, hiçbir zaman mutlak k2/k3 değeri değil.
Hiçbir gap/M_n/Z″ hesaplanmadı. 200a_bloklar.npz ve HUKUM_200A.json —gerçek
ölçümün çıktıları— bu görev kapsamında hiç üretilmedi.
