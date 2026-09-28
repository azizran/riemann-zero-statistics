# 200-B MAKİNE RAPORU (makine-inşa ajanı, 26 Eylül 2026)

Bu rapor `KALEM_YAKIN_CIFT_BASAMAK_26EYL2026.md`'nin makine-inşa görevinin sonucudur.
**Hiçbir gerçek yakın-çift ölçümü (200b_olcum.py + 200b_analiz.py'nin gerçek veri
üzerinde koşusu) yapılmadı ve hiçbir gerçek M_n / log(M_n/δ̃_n²) dağılım özeti (ortalama,
varyans, k-istatistiği, korelasyon, çeyreklik, histogram) HESAPLANMADI.** M1b/M2b (gerçek
γ_n, γ_{n+1} çiftlerinde gerçek M_n) ve M7b (gerçek sıfır KONUMLARI) KALEM'in açıkça
izin verdiği şekilde gerçek veri kullandı, ama yalnız motor-farkı / geçti-kaldı
istatistikleri ve pozisyon-türevli karışım düzeltmeleri raporlandı — hiçbir mutlak
k2/k3 ya da M dağılımı. M6b tamamen sentetik (CUE Monte Carlo + gerçekleştirilebilir
örnekleyiciler). 200b_olcum.py ve 200b_analiz.py yalnız SENTETİK/--fake girdilerle
uçtan uca test edildi.

## 1. İnşa edilen dosyalar (`200_configs/`)

| dosya | içerik |
|---|---|
| `200b_ortak.py` | Pencere tanımları (WINDOW_IDX, KALEM'de verilen sıfır-indeks aralıklarıyla birebir), olay çıkarımı (yalnız konum: γ_n,γ_{n+1},m_n,L_n,δ̃_n), tepe ölçümü (`peak_grid_parabola`: n-nokta iç ızgara + 3-nokta parabolik rafinasyon, \|Z\| üzerinde), tek-rung (k2,k3) blok jackknife (`jackknife_2`), doğrusal ara değer yardımcıları |
| `200b_olcum.py` | Gerçek ölçüm betiği (Görev 1) — **yalnız `--fake` ile test edildi**, gerçek veri üzerinde koşulmadı |
| `200b_analiz.py` | Gerçek analiz betiği (Görev 2: H-200B-1/1b/2/3, ikincil sınavlar) — **yalnız sentetik `200b_olaylar.npz` şemalı dosyalarla (RMT, hyb3, kaba a_k) uçtan uca test edildi**, gerçek `M6b_200B.json` kalibrasyonuyla da entegrasyon testi yapıldı |
| `200b_m1_m2_kapilar.py` | Kapı M1b + M2b (**GERÇEKTEN ÇALIŞTIRILDI**, yalnız fark/geçti-kaldı) → `M1b_200B.json`, `M2b_200B.json` |
| `200b_m6_sentetik.py` | Kapı M6b, 3 aşama (**GERÇEKTEN ÇALIŞTIRILDI**, tamamen sentetik) → `M6b_200B.json` |
| `200b_m7_karisim.py` | Kapı M7b (**GERÇEKTEN ÇALIŞTIRILDI**, yalnız sıfır konumları + model fonksiyonu) → `M7b_200B.json` |
| `scratchpad/k200b_makine/gen_synthetic_olaylar.py` | Test fixture üreticisi (200b_analiz.py'yi RMT/hyb3/a_k-benzeri sentetik verilerle sınamak için) |

`200_configs/200b_olaylar.npz` ve `200_configs/HUKUM_200B.json` **BİLEREK
OLUŞTURULMADI** (gerçek ölçüm/analiz koşulmadığı için ABSOLUTE RULES gereği; test
çıktıları yalnız `scratchpad/k200b_makine/test_*.npz` / `test_*.json` altında).

## 2. Kapı sonuçları

### M1b (tepe motoru: motor+ızgara+parabol vs mpmath altın-oran) — `M1b_200B.json`

Pencere başına 300 rastgele δ̃<0.2 olayı (tohum 2101), dps=20, tolerans 1e-12
(t'de relatif). Süre: 193 s (8 işçi, tek Pool tüm 900 olay için).

| pencere | \|ΔM\|/M max | \|ΔM\|/M p99 | Δk2 (motor−mpmath) | Δk3 (motor−mpmath) | GEÇTİ mi |
|---|---|---|---|---|---|
| W1 | 4.98e-05 | 1.23e-05 | +9.5e-08 | −1.4e-07 | evet |
| W2 | 4.21e-05 | 9.18e-06 | −6.7e-08 | +3.1e-08 | evet |
| W3 | 1.37e-05 | 6.60e-06 | +2.5e-07 | −3.3e-07 | evet |

Kapılar: \|ΔM\|/M ≤ 1e-4 ✓ (hepsi ≤5e-5); \|Δk2\| ≤ 0.002 ✓ (hepsi ~1e-7, kapının
~1e4 katı altında); \|Δk3\| ≤ 0.005 ✓ (hepsi ~1e-7). **M1b GENEL: GEÇTİ**, geniş
marjla — motor+33-ızgara+parabol yöntemi mpmath'ın yüksek-hassasiyetli altın-oran
maksimizasyonuyla pratikte ayırt edilemez (k-istatistiği düzeyinde fark, karar
eşiklerinin binlerce kat altında).

### M2b (ızgara çözünürlüğü: 33 vs 129 nokta, motor yalnız) — `M2b_200B.json`

**KALDI — pre-seal sapma, ayrıntı §3.** Aynı 900 olay, \|ΔM\|/M ≤ 1e-6 kapısı:

| pencere | \|ΔM\|/M max | \|ΔM\|/M p99 | tanı: Δk2 (33−129) | tanı: Δk3 (33−129) | GEÇTİ mi |
|---|---|---|---|---|---|
| W1 | 6.87e-06 | 4.90e-06 | −4.1e-08 | −7.9e-08 | **HAYIR** |
| W2 | 1.28e-05 | 7.67e-06 | −1.3e-07 | +5.2e-08 | **HAYIR** |
| W3 | 1.38e-05 | 6.82e-06 | +2.5e-07 | −3.3e-07 | **HAYIR** |

**M2b GENEL: KALDI** (eşiğin ~7–14 katı). Bağımsız bir yakınsama testiyle
(33→65→129→257→513 ızgara, 8 örnek olay) doğrulandı: 129-ızgara zaten neredeyse
tam yakınsamış (129-vs-513 farkı ~1e-8, yani 33-vs-129 farkının neredeyse tamamı
33-ızgaranın KENDİ ayrıklaştırma artığı, 129'un artığı değil). Yani bu bir kod hatası
değil, KALEM'in belirlediği 33-nokta ızgaranın gerçek Z(t) verisinde 1e-6 hedefini
karşılamayan (~5e-6 ile 1.4e-5 arası) gerçek bir ayrıklaştırma kalıntısı.

**Ancak** aynı tanı satırında görüldüğü gibi, bu kalıntının log(M/δ̃²)'nin k2/k3'üne
etkisi (33 vs 129 farkı) ~1e-7 mertebesinde — M1b'nin (motor vs mpmath, yani
"gerçeğe" karşı) verdiği farkla AYNI büyüklükte ve KALEM'in karar-kritik SE'lerinin
(havuzlanmış SE(D̄₃)≈0.004, §2 M6b) ve toleranslarının (0.0238 / 0.1422) **dört-beş
mertebe altında**. Sonuç: M2b'nin literal 1e-6 eşiği karşılanmıyor, ama bu artığın
200-B'nin asıl bilimsel kararları (H-200B-1/1b/2/3) üzerinde ölçülebilir hiçbir etkisi
yok. Bu bulgu koordinatöre/KALEM sahiplerine AÇIKÇA iletiliyor (bkz §3); 33-nokta
ızgara (Görev 1'in ta kendisi) DEĞİŞTİRİLMEDİ.

### M6b (sentetik güç + SE kalibrasyonu) — `M6b_200B.json`

**(i) Sonlu-ε CUE Monte Carlo tablosu** (Haar QR, 41-nokta ızgara+parabol — 200t_mc.py
mantığı; N∈{9,10,11,12}, her N'de ≥4×10⁵ olay ε̃<0.2'de; toplam 279 s, 4 işçi paralel):

| N | nmat | n(ε̃<0.2) | Δκ₂(ε=0.2) ± SE | Δκ₃(ε=0.2) ± SE | Pearson(δ̃,x) | δ̃<0.1 oranı |
|---|---|---|---|---|---|---|
| 9 | 5.30M | 401 421 | −0.00425 ± 0.00065 | +0.00369 ± 0.00072 | 0.0049 | 0.1280 |
| 10 | 4.80M | 403 871 | −0.00389 ± 0.00070 | +0.00390 ± 0.00077 | 0.0059 | 0.1280 |
| 11 | 4.35M | 402 010 | −0.00383 ± 0.00077 | +0.00316 ± 0.00100 | 0.0067 | 0.1274 |
| 12 | 4.00M | 404 205 | −0.00399 ± 0.00086 | +0.00365 ± 0.00112 | 0.0083 | 0.1282 |

(SE: blok-jackknife, her N için 80–106 blok/batch=50000 matris. ε=0.1, 0.3 bantları
da aynı çalışmadan çıkarıldı, `M6b_200B.json`'da tam tablo. δ̃<0.1 oranı U^{1/3}
teorisinin 1/8=0.125 tahminine yakın, teoremle tutarlı.)

**Bağımsız teftiş karşılaştırması (koordinatörün istediği, §4):** sonuçlar auditor
değerleriyle ~2 birleşik-SE içinde uyuşuyor — bkz §4, tüm N'lerde ≤1.6σ.

**(ii) f-kalibrasyonu** (gerçek olay SAYILARI — `200b_sayimlar.json`, pozisyon-türevli
— ile 2 örnekleyici × 3 pencere × 3 ε × 200 tekrar, 64-blok jackknife; 90 s, 8 işçi):
f ∈ [0.940, 1.128] (tüm pencere×ε×kümülant×örnekleyici birleşimlerinde) — **200-A'nın
bulgusuyla tutarlı**: blok-jackknife'ın kendisi zaten iyi kalibre (yalnız %6-13
düzeltme gerekiyor), KALEM'in "dar formül 2.7× iyimser" ön-uyarısı NAIF formül için
geçerliydi, jackknife'ın kendisi için değil (bkz §3).

**(iii) Seçim gücü** (Gauss yaklaşıklığı, 20 000 çekiliş, kalibre SE ile): havuzlanmış
SE(D̄₃)≈0.00415, SE(D̄₂)≈0.00409 (pencere başına SE(D₃): W1 0.0111, W2 0.0084, W3
0.0053 — KALEM'in ön-kestirimi ≈0.005 (salt CUE) ile aynı mertebede, ≈0.014
(sağa-çarpık) tahmininden daha iyimser — bu da (i)'deki gibi f≈1 kalibrasyonunun
doğal sonucu, bkz §3). H-200B-3 karışıklık matrisi: SABİT→SABİT %100, GEOMETRİK→GEOMETRİK
%99.98, LİN/RMT→SÖNMÜŞ %100 (ikisi zaten aynı sınıfta birleşiyor). **P(doğru sınıf)<0.80
olan hiçbir hipotez/senaryo çifti YOK** — M6b'nin güç kapısı GEÇTİ, KALEM'in beklediği
gibi sınıflar arası ayrım geniş.

### M7b (karışım düzeltmesi, yalnız sıfır konumları) — `M7b_200B.json`

Her pencere × ε için gerçek olayların (δ̃<ε) L_n dağılımı üzerinden E[κ_r(L,2)] −
κ_r(L̄_W,2): tüm 18 (pencere,ε,kümülant) biriminde \|düzeltme\| ≤ 0.00028 — eşiğin
(0.002) çok altında. **M7b: hiçbir düzeltme uygulanmaz**, 200b_analiz.py'de nokta-L
baseline κ_r^{CUE,2}(L̄_W) doğrudan kullanılıyor (zaten öyle inşa edildi).

## 3. KALEM'den sapmalar ve bulgular

1. **M2b (33 vs 129 ızgara) 1e-6 eşiğini karşılamıyor (§2).** Kod hatası değil —
   bağımsız 513-nokta yakınsama testiyle doğrulandı. Etki büyüklüğü karar-kritik
   SE/toleransların dört-beş mertebe altında (M1b'nin motor-vs-mpmath farkıyla aynı
   ölçekte). Görev 1'in belirlediği 33-nokta ızgara DEĞİŞTİRİLMEDİ (bu ajanın yetkisi
   dışında); bulgu KALEM sahiplerine/koordinatöre bu raporla iletiliyor —
   200-A'nın C0→C1 sürecindeki gibi, tasarım değişikliği gerekip gerekmediğine
   dair karar onlara ait.
2. **f-kalibrasyonu (~1.0–1.13) KALEM'in ön-kestirdiği SE(D̄₃)≈0.014 (sağa-çarpık
   senaryo) tahmininden belirgin biçimde daha iyimser (gerçek: ≈0.004).** Bu,
   200-A'nın KENDİ M6 bulgusuyla (f≈1, "dar formül" 2.7× iyimserdi ama jackknife'ın
   kendisi zaten doğru kalibre) birebir tutarlı — yeni bir anomali değil, aynı
   düzenin b=2'de tekrarı. Yine de bu, KALEM'in "Ön kestirim" bölümündeki SE
   sayılarının (≈0.005 / ≈0.014) gerçek koşulan kalibrasyonla değiştirilmesi
   gerektiği anlamına gelir — 200b_analiz.py zaten GERÇEK jackknife+f kullanıyor
   (ön-kestirim yalnız KALEM metnindeki bir planlama sayısıydı, karar hesabında
   kullanılmıyor).
3. **M6b(ii) f-kalibrasyonu, KALEM'in "gerçek olay sayılarıyla" ifadesini HER ε İÇİN
   AYRI (üç farklı n) olarak yorumladı** (200b_sayimlar.json zaten üç ε için ayrı n
   veriyor) — KALEM'de açıkça belirtilmiyor ama en doğal okuma bu (n farklı ε'lerde
   147'den 33713'e değişiyor, tek bir f ile temsil edilemez).
4. **M8b kuralı pencere × kümülant (k2, k3 AYRI AYRI) uygulandı**, KALEM metni
   "o pencerede İKİSİNİN BÜYÜĞÜ" derken tam olarak hangi granülerlikte
   belirtmiyor; bu ajanın tutucu (conservative) okuması, 200-A'nın 12-bileşenli
   toplu anahtarından farklı olarak burada yalnız 2 istatistik (k2,k3) olduğu için
   pencere×kümülant bazında ayrı SE kullanmak daha az bilgi kaybettiriyor.
5. **200b_analiz.py'nin H-200B-1b "pooled corrected κ3" SE'si, D̄3'ün SE'siyle AYNI
   alındı** (matematiksel gerekçe: κ_r^{CUE,2}(L̄_W) deterministik/sabit bir
   teorik değer, rastgele değişken değil, dolayısıyla k3−Δκ3'ün varyansı D3'ünkiyle
   birebir aynı) — KALEM bunu açıkça yazmıyor ama tanımdan (D3 := k3−baseline−Δκ3)
   doğrudan çıkıyor.
6. **M6b(iii) güç tahmini NAİF-iid formülü (f-kalibreli) kullanıyor**, GERÇEK
   blok-jackknife değil (KALEM: "Gauss yaklaşıklığı" zaten bunu ima ediyor) — asıl
   karar (200b_analiz.py) tam blok-jackknife (M8b) kullanacak; stage(iii) yalnız
   ÖN güç tahmini, karar hesaplamasının kendisi değil.

## 4. Bağımsız teftiş karşılaştırması (M6b(i), koordinatör talebi)

Koordinatörün bağımsız denetçisi (≈3–4×10⁵ olay/N, ε̃=0.2) ile bu ajanın SONUÇLARI
(bağımsız olarak, denetçinin sayıları KOPYALANMADAN önce hesaplandı) arasındaki fark,
iki ölçümün birleşik SE'sine (√(SE_mine²+SE_auditor²), auditor SE'si mine ile aynı
mertebede varsayılarak √2×SE_mine) oranı:

| N | Δκ₂ (bu ajan) | Δκ₂ (denetçi) | fark/birleşik-SE | Δκ₃ (bu ajan) | Δκ₃ (denetçi) | fark/birleşik-SE |
|---|---|---|---|---|---|---|
| 9 | −0.00425 | −0.00458 | 0.36σ | +0.00369 | +0.00409 | 0.39σ |
| 10 | −0.00389 | −0.00544 | 1.57σ | +0.00390 | +0.00437 | 0.43σ |
| 11 | −0.00383 | −0.00422 | 0.36σ | +0.00316 | +0.00395 | 0.56σ |
| 12 | −0.00399 | −0.00398 | 0.01σ | +0.00365 | +0.00279 | 0.54σ |

**Sonuç: TÜM 8 karşılaştırma ~2 birleşik-SE içinde** (en büyüğü N=10 κ₂'de 1.57σ,
uyarı eşiği olan 2σ'nın altında). Hiçbir tutarsızlık FLAG'lenmedi. Bu, hem bu ajanın
hem denetçinin bağımsız CUE MC uygulamalarının aynı sonlu-ε düzeltmesini bulduğunun
bir çapraz-doğrulamasıdır.

## 5. Gerçek ölçüm için çalışma zamanı kestirimi

`200b_olcum.py` (gerçek): toplam olay (δ̃<0.3, üç pencere) = 3916+9445+33713 = 47 074;
33-nokta ızgara ⇒ ~1.553M Z-değerlendirmesi. Ölçülen motor hızı (W3 ölçeğinde, en
yavaş pencere, GERÇEK t aralığında SENTETİK-olmayan t noktalarıyla ama hiçbir
sonuç saklanmadan salt hız testi): **16.3 s toplam** (10.5 µs/nokta). Python/IO
yüküyle birlikte gerçekçi tahmin: **< 1 dakika**, 10 dakika sınırının rahatça
altında.

`200b_analiz.py` (gerçek): 47 074 olayda üç ε için k-istatistiği + M8b (32/64/128
blok jackknife) + havuzlama + karar kuralları — saf numpy, **saniyeler mertebesinde**.

M1b/M2b gerçek koşu: **193 s** (ölçüldü, 900 olay, 8 işçi, mpmath altın-oran baskın
maliyet). M6b(i) gerçek koşu: **279 s**. M6b(ii): **90 s**. M6b(iii): saniyeler.
M7b: saniyeler. **Hepsi tek tek < 10 dakika** (ABSOLUTE RULES).

## 6. Açık ifade: gerçek-veri körlüğü

Bu görev boyunca **hiçbir gerçek M_n ya da log(M_n/δ̃_n²) örnekleminin dağılım özeti
(ortalama, varyans, k-istatistiği, çeyreklik, histogram, korelasyon) hesaplanmadı
veya yazdırılmadı.** M1b/M2b gerçek M_n değerlerini hesapladı (KALEM'in açıkça izin
verdiği gibi) ama yalnız motor-mpmath / 33-129 ızgara FARKLARI (\|ΔM\|/M istatistikleri,
k2/k3 FARKLARI — asla mutlak k2/k3) raporlandı/kaydedildi. M7b yalnız gerçek sıfır
KONUMLARINI (γ_n) ve analitik kap(r,L,2) fonksiyonunu kullandı — hiçbir Z/M değeri
hesaplamadı. M6b tamamen sentetik (CUE Monte Carlo + gerçekleştirilebilir örnekleyiciler,
gerçek veriyle sıfır ilişki; yalnız gerçek olay SAYILARI — pozisyon-türevli, zaten
`200b_sayimlar.json`'da açık — kullanıldı). `200b_olaylar.npz` ve `HUKUM_200B.json`
—gerçek ölçümün/analizin çıktıları— bu görev kapsamında hiç üretilmedi; ölçüm betiği
yalnız ONKAYIT_200B push'landıktan SONRA çalıştırılacak (KALEM M9b).
