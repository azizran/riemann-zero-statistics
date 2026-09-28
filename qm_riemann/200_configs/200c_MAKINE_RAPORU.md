# 200-C MAKİNE RAPORU (makine-inşa ajanı, 27 Eylül 2026)

Bu rapor `KALEM_MERDIVEN_YUKSEKLIK_27EYL2026.md`'nin makine-inşa görevinin sonucudur.
**Hiçbir gerçek ölçüm (200c_olcum.py + 200c_analiz.py'nin gerçek veri üzerinde koşusu)
yapılmadı ve hiçbir gerçek log|Z(t)|, log|Z'(γ_n)|, M_n ya da log(M_n/δ̃_n²) dağılım
özeti (ortalama, varyans, k-istatistiği, çeyreklik, histogram, korelasyon) HİÇBİR
YERDE hesaplanmadı/yazdırılmadı.** M0c (yalnız konum/sha256/sayım), M1c/M2c (motor
doğruluğu — KALEM'in açıkça izin verdiği gibi gerçek Z/Z' hesaplandı ama yalnız
motor-mpmath FARKLARI ve k2/k3 FARKLARI raporlandı), M3c (tepe bulucu — yalnız
129-vs-513 FARKI), M6c (tamamen sentetik CUE Monte Carlo + gerçekleştirilebilir
örnekleyiciler; yalnız gerçek olay SAYILARI — pozisyon-türevli — kullanıldı), M7c
(yalnız sıfır KONUMLARI + analitik kap() fonksiyonu) GERÇEKTEN ÇALIŞTIRILDI (kapı
tanımı gereği zaten öyle). 200c_olcum.py ve 200c_analiz.py yalnız `--fake` / sentetik
şema-uyumlu dosyalarla uçtan uca test edildi.

**Özet (tüm kapılar GEÇTİ):** M0c GEÇTİ (0.3s) · M1c GEÇTİ (135.9s) · M2c GEÇTİ ·
M3c GEÇTİ (310.8s) · M6c GERÇEKTEN ÇALIŞTIRILDI (i: 441s, ii: 181s @10 tekrar,
iii: saniyeler) · M7c GEÇTİ (flag yok). **H-200C-1 eşiği 3.5'te DONDU** (KALEM'in
ikili kuralı hiçbir aday eşikte karşılanmadı — bkz §5.3). Çapa sayıları: C1 7381,
C2 1253, C3 239, C4 24. b1 kararı: C1/C2/C3 tümü, **C4 alt-örnek** (1e6, tohum 301).

## 1. İnşa edilen dosyalar (`200_configs/`)

| dosya | içerik |
|---|---|
| `200c_motor.py` | Çapalı Riemann–Siegel Z(t), Z'(t) (tamsayı taban + float64 ofset; mpmath çapa fazı + float64 Taylor açılımı; C0+C1 kalan, 200a_motor'un Chebyshev katsayılarını YENİDEN KULLANIR) |
| `200c_veri.py` | (Görev tarafından verilen, DEĞİŞTİRİLMEDİ) LMFDB/Platt sıfır dosyası okuyucu |
| `200c_ortak.py` | Pencere/dosya tanımları, sıfır yükleme (T0+ofset), olay çıkarımı (200b_ortak.events_from_zeros'un ofset-uyarlanmış sürümü), çapa tablosu önbelleği |
| `200c_m0_veri_butunlugu.py` | Kapı M0c (**GERÇEKTEN ÇALIŞTIRILDI**) → `M0c_200C.json` |
| `200c_m1_m2_kapilar.py` | Kapı M1c + M2c (**GERÇEKTEN ÇALIŞTIRILDI**, yalnız motor-farkı/geçti-kaldı) → `M1c_200C.json`, `M2c_200C.json` |
| `200c_m3_tepe.py` | Kapı M3c (**GERÇEKTEN ÇALIŞTIRILDI**, yalnız 129-vs-513 farkı) → `M3c_200C.json` |
| `200c_m6_sentetik.py` | Kapı M6c, 3 aşama (**GERÇEKTEN ÇALIŞTIRILDI**, tamamen sentetik + gerçek olay SAYILARI) → `M6c_200C.json` |
| `200c_m7_karisim.py` | Kapı M7c (**GERÇEKTEN ÇALIŞTIRILDI**, yalnız sıfır konumları + analitik kap()) → `M7c_200C.json` |
| `200c_olcum.py` | Gerçek ölçüm betiği (Görev 6) — **yalnız `--fake` ile test edildi** |
| `200c_analiz.py` | Gerçek analiz betiği (Görev 7, H-200C-1..4 + M8c blok-duyarlılık tanısı dahil) — **yalnız 200c_olcum.npz şemasıyla uyumlu SENTETİK dosyalarla (ve gerçek pencere konumlarıyla karışık bir ek testle) test edildi** |
| `200c_tahmin.json`, `200c_tahmin.py` | (Görev tarafından verilen, DEĞİŞTİRİLMEDİ) hipotez tabloları |
| `scratchpad/k200c_makine/*` | Test fixture üreticileri, hız kıyaslamaları, hata-teşhis betikleri |

`200_configs/200c_olcum.npz` ve `HUKUM_200C.json` (gerçek ölçümün/analizin çıktıları)
**BİLEREK OLUŞTURULMADI** — test çıktıları yalnız `scratchpad/k200c_makine/test_*`
altında.

## 2. Motor (200c_motor.py) — tasarım, bulunan hata, doğrulama

### 2.1 Çapa aralığı kuralı ve sayısı (pencere başına)

Kübik terim sınırı: |θ'''(t_a)|·h³/6 ≤ 1e-10 rad, θ'''(t)≈−1/(2t²) baskın terimiyle,
h=(12·1e-10·t_min²)^(1/3), aralık=2h (aşağı yuvarlanmış tam sayı), pencerenin ALT
UCUNDAKİ t_min kullanılarak MUHAFAZAKAR/worst-case:

| pencere | t_min | h (float) | çapa aralığı (tam sayı) | çapa sayısı (tüm pencere) | N_max (RS terim sayısı) |
|---|---|---|---|---|---|
| C1 | 8 846 000 | 45.45 | 90 | **7 381** | 1 238 |
| C2 | 99 146 000 | 227.64 | 455 | **1 253** | 3 991 |
| C3 | 997 946 000 | 1061.20 | 2 122 | **239** | 12 613 |
| C4 | 30 599 546 000 | 10396.10 | 20 792 | **24** | 69 794 |

Çapa kurma süresi (mpmath, dps=30, 8 işçi paralel, tüm pencere): C1 6.30s, C2 3.57s,
C3 2.54s, C4 2.27s — **toplam ~15s** (4 pencere). Bellek: `phase` tablosu (çapa ×
N_max float64) C1 69.7MB, C2 38.2MB, C3 23.0MB, C4 12.8MB.

### 2.2 Mühür-öncesi bulunan hata (pre-seal deviation) — `MP_TWO_PI` donması

İlk sürümde `MP_TWO_PI = 2 * mp.pi` MODÜL YÜKLENİRKEN (mpmath varsayılan dps=15'te)
bir kez hesaplanıyordu; sonraki `mp.mp.dps = 30` ayarları bu DONMUŞ sabiti
etkilemiyordu. Sonuç: C3/C4 gibi büyük t'lerde `fmod(θ(t_a), 2π)` ve
`fmod(t_a·log n, 2π)` indirgemelerinde sistematik faz hatası (gözlenen: C4'te
|ΔZ|≈1.3e-4, |ΔZ'|≈1.5e-3 — M1c eşiklerini (2e-5/2e-4) aşıyordu). Teşhis
(`scratchpad/k200c_makine/debug_c4_phase*.py`): n=1'de (log(1)=0, yalnız θ(t_a)
etkili) bile dps=30 tablo değeri dps=50 taze hesaptan 1.33e-5 rad sapıyordu;
n arttıkça sapma küçülüyordu (θ(t_a) payının fazdaki payı sabit ama t_a·log(n)
payı n'e göre değişiyor) — bu desen "yalnız θ(t_a)'nın kendisi yanlış" ipucunu
verdi, izole test 2π'nin donduğunu doğruladı. **Düzeltme:** `2*pi` artık HER
zaman `mp.mp.dps` ayarlandıktan HEMEN SONRA, fonksiyon içinde taze hesaplanır
(`_anchor_worker`); hiçbir modül-seviyesi dondurulmuş mpmath sabiti yok.
Düzeltme sonrası tüm pencerelerde motor-mpmath farkı ~1e-9–1e-10'a düştü (bkz M1c).

### 2.3 M1c/M2c sonuçları (GERÇEK mpmath karşılaştırması, dps=25)

200 rastgele t (tohum 3101; yarısı geniş bir adaydan [2000 nokta] çapaya EN UZAK
|dt|'ye sahip 100 nokta + 100 sıradan rastgele) ve 200 rastgele sıfır (tohum 3102)
/ pencere, havuzlu mpmath (8 işçi), toplam 1600 nokta, **118s mpmath + 135.9s toplam**:

| pencere | \|ΔZ\| max (t) | \|ΔZ'\| max (t) | Δk2 (t) | Δk3 (t) | \|ΔZ\| max (sıfır) | M2c en kötü oran |
|---|---|---|---|---|---|---|
| C1 | 1.19e-09 | 5.12e-09 | +1.3e-10 | −1.3e-09 | 1.40e-09 | 9.79e-11 |
| C2 | 1.20e-09 | 1.83e-08 | +1.5e-05 | −2.6e-04 | 1.54e-09 | 7.28e-11 |
| C3 | 7.98e-10 | 1.33e-08 | +7.8e-10 | −3.5e-09 | 7.21e-10 | 6.44e-11 |
| C4 | 2.59e-09 | 2.23e-08 | −1.1e-09 | +6.1e-09 | 7.21e-10 | 7.11e-11 |

Gates: |ΔZ|≤2e-5 ✓ (4 mertebe marjla), |ΔZ'|≤2e-4 ✓, |Δk2|≤0.002 ✓, |Δk3|≤0.005 ✓
(C2 t-örneğinde en büyük, 2.6e-4, yine de eşiğin ~20 katı altında). **M1c GENEL:
GEÇTİ. M2c GENEL: GEÇTİ** (tüm sıfırlarda |Z|/|Z'|≤1e-4, en kötü oran ~1e-10).

### 2.4 M3c sonuçları (tepe bulucu, 129 vs 513 ızgara, GERÇEK)

200 rastgele δ̃<0.2 olayı / pencere (tohum 3103):

| pencere | olay sayısı (δ̃<0.2) | \|ΔM\|/M max | süre |
|---|---|---|---|
| C1 | 11 515 | 1.50e-07 | 2.5s |
| C2 | 11 546 | 2.25e-07 | 10.4s |
| C3 | 12 051 | 2.22e-07 | 37.2s |
| C4 | 12 143 | 1.83e-07 | 243.9s |

Gate 1e-6 — hepsi ~4 mertebe altında. **M3c GENEL: GEÇTİ** (toplam 310.8s). Olay
sayıları KALEM'in tablosuyla (11515/11546/12051/12143 @ δ̃<0.2) BİREBİR eşleşiyor —
bağımsız çapraz doğrulama.

## 3. M0c (veri bütünlüğü, GERÇEK, GEÇTİ)

| pencere | sha256 (ilk 16) | blok | sıfır | t aralığı | açılmış ort. aralık | L̄ (pilot) | sayım farkı |
|---|---|---|---|---|---|---|---|
| C1 | a3f3dd838cc50b45 | 317 | 1 500 000 | [8846000.35, 9509986.68] | 1.000000 | 14.1943 | +0.324 |
| C2 | f181f573866fb1cd | 271 | 1 500 000 | [99146000.11, 99714541.97] | 1.000000 | 16.5771 | +0.842 |
| C3 | fb696127ea715c34 | 238 | 1 500 000 | [997946000.12, 998445098.60] | 1.000000 | 18.8836 | +0.445 |
| C4 | 2c5e764cedb2a0d0 | 202 | 1 500 000 | [30599546000.12, 30599968514.72] | 1.000000 | 22.3064 | +0.242 |

Sayım farkları (S(T) mertebesi) KALEM'in pilot değerleriyle (+0.32/+0.84/+0.45/+0.24)
BİREBİR eşleşiyor. Kesin artanlık: hepsi TRUE. **M0c GENEL: GEÇTİ** (0.3s).

## 4. b1 (all-vs-subsample) kararı

`200c_olcum.py`'nin M3c'nin GERÇEK 200-olay hız ölçümünden (129+513 ızgara)
çıkarılan per-nokta maliyeti: C1=0.0195ms, C2=0.081ms, C3=0.29ms, C4=1.9ms.
1e6/1.5e6 nokta için projeksiyon:

| pencere | b0 (1e6) projeksiyon | b1-tümü (1.5e6) projeksiyon | KALEM eşiği (30dk) | karar |
|---|---|---|---|---|
| C1 | 19.5s | 29s | altında | **tümü** |
| C2 | 81s | 122s | altında | **tümü** |
| C3 | 290s (4.8dk) | 435s (7.3dk) | altında | **tümü** |
| C4 | 1900s (31.7dk) | 2850s (47.5dk) | **AŞIYOR** | **alt-örnek (1e6, tohum 301)** |

`200c_olcum.py`'de `B1_DECISION = {"C1":"all","C2":"all","C3":"all","C4":"subsample"}`
olarak kayıtlı; `200c_m6_sentetik.py`'nin f-kalibrasyonu (stage ii) AYNI kararı
kullanır (bir assert ile tutarlılık zorlanır).

## 5. M6c (sentetik güç + SE kalibrasyonu + eşik donma)

### 5.1 Aşama (i): sonlu-ε CUE Monte Carlo tablosu, N ∈ {14,17,19,22}

200b_m6_sentetik.py'nin AYNI mantığı (Haar QR, 41-nokta ızgara+parabol) YENİDEN
KULLANILDI (N_list ve hedef farklı). Hedef ≥4×10⁵ olay/N @ ε̃<0.2, 4 N paralel
(8 işçi havuzu, yalnız 4 iş aktif kullanıldı):

| N | matris sayısı | olay (ε̃<0.2) | süre | Δκ2(ε=0.2) | Δκ3(ε=0.2) |
|---|---|---|---|---|---|
| 14 | 3 400 000 | 403 974 | 276s | −0.00637 | +0.00662 |
| 17 | 2 800 000 | 404 463 | 336s | −0.00618 | +0.00278 |
| 19 | 2 500 000 | 403 977 | 378s | −0.00309 | +0.00341 |
| 22 | 2 150 000 | 402 191 | 440s | −0.00556 | +0.00452 |

Toplam duvar-saati **441s (~7.35dk)** — 4 N paralel koştuğu için en yavaş (N=22)
belirleyici. KALEM'in kendi notu ("N=22 için uzun sürebilir; bütçele/paralelleştir")
ile tutarlı; tek bir hesaplama olarak ~15dk sınırının altında.

### 5.2 Aşama (ii): f-kalibrasyonu (3 basamak × 4 pencere)

3 basamak (b=0,1,2) × 2 örnekleyici (saf eğik-CUE, sağa-çarpık +Euler p=2,3) ×
4 pencere × tekrar, GERÇEK n (b0=1e6, b1 `200c_olcum.B1_DECISION`'a göre
1.5e6/1e6, b2≈11.5–12.1k @ ε=0.2) ve GERÇEK t-aralıklarıyla, 64-blok jackknife.

**Mühür-öncesi sapma (tekrar sayısı 200→40→10):** İlk koşu (KALEM/200-A/B
geleneğine uyan 200 tekrar) rung=1'de (sıfırlar, n≈1–1.5×10⁶) b=1 için
`200a_orneklem.sample_log_tilted_product`'ın disk-çarpanı reddetme-örneklemesi
(`_oversample_mult[1]=6.0`, 200-A/B'nin L≈9–12 aralığı için kalibre edilmiş)
200-C'nin ÇOK DAHA GENİŞ L aralığında (14–22.3, j çarpan sayısı ~L'ye kadar)
beklenenden çok daha yavaş yakınsadı (200 tekrar 20+ dk'da bitmedi, 40 tekrar
~11.5 dk'da bitmedi — kesin kanıt yok ama örüntü, büyük j'de düşük kabul
oranına işaret ediyor). **Tekrar sayısı 10'a düşürüldü** (n GERÇEK sayılarla
AYNI kaldı — yalnız kaç kez tekrarlandığı azaldı); bu koşu **181s (3.0dk)**'de
GERÇEKTEN tamamlandı (240 iş: 3×2×4×10).

f çarpanları (SD(10 tekrar)/ortalama(SE_jk), elementwise max(cue,skew)) aralığı:
**[0.92, 1.84]** — 200-A/B'nin bulduğu ~[0.9,1.3] aralığından biraz daha geniş
(10 tekrarla SD tahmininin kendisi gürültülü — bkz §7 sapma notu), ama aynı
mertebede ve f≈1 civarında (jackknife'ın temelde iyi kalibre olduğu bulgusunu
DOĞRULUYOR, "dar formül iyimserdi" uyarısı NAİF formül için geçerliydi).

### 5.3 Aşama (iii): H-200C-1 gücü + eşik donma + H-200C-2 karışıklık

4000 çekiliş/senaryo, kalibre SE + σ_sys ile pencere-başına 3×3 kappa3 alt-blok
kovaryans (basit sabit cross-rung korelasyon varsayımı corr=0.3 — GERÇEK analiz
gerçek ortak-blok jackknife kullanacak):

| eşik | FP oranı (H_C) | güç (H_Sγ, γ=0.48) | güç (H_S1, γ=1) |
|---|---|---|---|
| 3.0 | 0.002 | 0.377 | 0.922 |
| **3.5** | **0.000** | **0.192** | **0.806** |
| 4.0 | 0.000 | 0.083 | 0.629 |
| 5.0 | 0.000 | 0.010 | 0.222 |

**Hiçbir eşik KALEM'in ikili kuralını (FP<%1 VE güç(H_Sγ)≥0.80) karşılamıyor**
(FP her zaman rahatça düşük ama güç(H_Sγ) hiçbir zaman %80'e ulaşmıyor — en
iyisi eşik=3.0'da %37.7). KALEM'in kendi yedek kuralı gereği ("olmazsa eşik
[3,5] içinde M6c kuralıyla ayarlanıp 3.5 DONAR") **eşik 3.5'te DONDU**
(kod: `frozen_threshold` alanı, varsayılan 3.5'e düşüyor). Not: γ=1 (H_S1,
a_k'ya tam yakınsama) senaryosunda güç zaten iyi (%80.6 @ 3.5); asıl güçsüzlük
KALEM'in kendi kör tahmini γ̂=0.48±0.12'nin (H_Sγ) MERKEZ değerine karşı —
bu, 200-C ölçümünün H-200C-1'i SONLU-YÜKSEKLİK olarak KESİN damgalamak için
yeterli güce sahip OLMAYABİLECEĞİ, ama SABİT YAPI'yı ayırt etmek için (FP≈0)
güvenilir olduğu anlamına geliyor — H-200C-4 (asimptot, 7-pencere) bu yüzden
birincil çapraz denetim olarak KALEM'de doğru şekilde öne çıkarılmış.

H-200C-2 karışıklık matrisi (χ²-en-iyi sınıf seçimi, aynı kovaryans varsayımıyla):
H_C→H_C %100, H_Sγ→H_Sγ %100, H_S1→H_S1 %100 — üç sınıf ayrımı net (KALEM'in
"σ_h" terimiyle üç hipotez tablosu birbirinden yeterince ayrışıyor).

f-kalibrasyonu ve güç aşamaları toplam **~3.5 dakika** sürdü (stage ii 181s +
stage iii saniyeler).

## 6. M7c (karışım düzeltmesi)

M6c'nin sonlu-ε tablosu hazır olduktan sonra yeniden çalıştırıldı (`m6c_available:
true`). Hem taban (kap(r,L,b)) hem de Δκ_r^CUE(ε=0.2,L) düzeltmeleri hesaplandı:

| pencere | b0 k2 | b0 k3 | b1 k2 | b1 k3 | b2 Δκ2 (ε=0.2) | b2 Δκ3 (ε=0.2) |
|---|---|---|---|---|---|---|
| C1 | −2.9e-07 | +6.5e-08 | −4.0e-07 | +8.6e-08 | ~9e-19 | ~9e-19 |
| C4 | (benzer, ~1e-15) | | | | 0.0 | ~9e-19 |

Tüm düzeltmeler makine-hassasiyeti düzeyinde sıfır (pencere-içi L yayılımı çok
küçük olduğundan — C1'de bile ΔL≈0.001 mertebesinde ortalama sapma, KALEM'in
belirttiği 0.072'nin kendisi δ̃<0.3 olaylarının L yayılımı, ama ORTALAMA'nın
ETRAFINDAKİ simetrik dağılım nedeniyle E[κ(L)]−κ(L̄) ≈ 0 — kap() fonksiyonu bu
L aralığında neredeyse doğrusal). **Eşik (0.002) aşan HİÇBİR (pencere,basamak,
kümülant[,ε]) yok** → **M7c: hiçbir düzeltme uygulanmaz**, 200c_analiz.py
zaten nokta-L̄ taban + Δκ(ε,L̄) kullanıyor (tasarım gereği).

## 7. KALEM'den sapmalar ve tasarım kararları

1. **`MP_TWO_PI` donma hatası** (§2.2) — mühür-öncesi bulundu ve düzeltildi; görev
   metnini DEĞİŞTİRMEDİ (yalnız uygulama hatasıydı), motoru gerçek veriyle yeniden
   doğruladı (M1c/M2c/M3c hepsi GEÇTİ).
2. **Çapa aralığı worst-case (pencere t_min) kullanılarak SABİT tutuldu** — KALEM
   "anchors spaced so that... report it and the number of anchors" diyor (tekil bir
   sayı ima ediyor); her çapaya kendi t_a'sına göre DEĞİŞKEN aralık vermek (daha az
   çapa) yerine, basit/raporlanabilir tek bir aralık tercih edildi (muhafazakar,
   fazladan çapa maliyeti ihmal edilebilir — toplam kurma süresi zaten ~15s).
3. **f-kalibrasyonu 3 basamak (b0,b1,b2) için AYRI AYRI hesaplandı** (200b'nin
   yalnız b2 için yaptığının genellemesi) — KALEM'in "GLS covariance = f-kalibre
   jackknife 3x3 cross-rung blok" ifadesi zaten 3 basamağın ortak kalibrasyonunu
   gerektiriyor.
4. **H-200C-1/H-200C-4 GLS uyumu** `scipy.optimize.least_squares` + Cholesky-
   beyazlatma ile yapıldı (200c_tahmin.py'nin kendi GLS yönteminin genellemesi).
5. **M6c(iii) güç simülasyonunun kovaryansı** (cross-rung, pencere başına 3x3)
   basit bir SABİT korelasyon katsayısı (corr=0.3) varsayımıyla modellendi — GERÇEK
   analiz (200c_analiz.py) GERÇEK ortak-blok jackknife kovaryansını kullanacak; bu
   yalnız ÖN güç/eşik-donma tahmini (200b_m6_sentetik'in stage(iii)'ü ile aynı
   ruhta bir basitleştirme, orada da "GERÇEK karar SE'leri analiz.py'de blok-
   jackknife'tan gelecek" notu var).
6. **`200c_olcum.py`'nin b1 kararı pencere-başına SABİT kodlandı** (dinamik hız
   ölçümü yerine) — M3c'nin GERÇEK hız verisinden türetildi, KALEM'in "decide from
   your benchmark and record the decision" talimatına birebir uyuyor.
7. **M6c(ii) tekrar sayısı 200 → 10'a düşürüldü (mühür-öncesi sapma, gerekçeli).**
   200-A/B geleneği 200 tekrardı; bu, 200-C'nin GENİŞ L aralığında (14–22.3, disk-
   çarpanı reddetme-örneklemesinin `_oversample_mult` sabitleri 200-A/B'nin dar
   L≈9–12 aralığı için kalibre edilmişti) rung=1 (b1, n≈1–1.5×10⁶) örneklemesini
   beklenenden ÇOK yavaşlattı: 200 tekrar >20 dk'da, 40 tekrar >11.5 dk'da
   bitmedi (ikisi de iptal edildi — bkz §5.2 ve `scratchpad/k200c_makine/
   m6c_stage_ii_v*.log`). **10 tekrarla GERÇEKTEN 181 saniyede tamamlandı**
   (n'ler GERÇEK kaldı, yalnız kaç kez tekrarlandığı azaldı). Sonuç: f
   çarpanları (SD(10)/ortalama(SE_jk)) 200-A/B'nin 200-tekrarlı tahminlerinden
   DAHA GÜRÜLTÜLÜ (10 noktalı bir SD tahmini doğası gereği belirsiz — göreli
   hata ~1/√(2·9)≈24%) ama aynı mertebede (~[0.9,1.8]) ve f≈1 civarında olma
   eğilimini DOĞRULUYOR. Bu, KALEM sahiplerine/koordinatöre AÇIKÇA iletilir:
   gerçek ölçümden ÖNCE, zaman kısıtı olmayan bir oturumda M6c(ii) daha yüksek
   tekrar sayısıyla (örn. 50-100, `--n-reps-ii`) YENİDEN çalıştırılıp
   `M6c_200C.json` güncellenebilir — 200c_analiz.py bunu otomatik okur, hiçbir
   kod değişikliği gerekmez.
8. **M8c (200-B kuralı, blok sayısı duyarlılığı)** `200c_analiz.py`'ye TANI
   olarak eklendi (32/64/128 blok, 64→128 göreli SE değişimi >%20 ise flag) —
   ancak KALEM'in kendisi `200c_olcum.py` için nblocks=128'i ZATEN SABİTLEDİĞİNDEN
   (Görev 6: "common 128 equal-t-length blocks per window"), 200-A/B'nin aksine
   burada birincil blok sayısını DEĞİŞTİRMEZ — yalnız 128'in 64'e göre yakınsadığını
   doğrular ve raporlar (sentetik testte C1 için 64→128 değişimi ~%10, eşiğin
   altında).

## 8. Açık ifade: körlük

Bu görev boyunca **hiçbir gerçek log|Z(t)|, log(|Z'(γ_n)|·2π/L_n), M_n ya da
log(M_n/δ̃_n²) örnekleminin dağılım özeti (ortalama, varyans, k-istatistiği,
çeyreklik, histogram, korelasyon) hesaplanmadı veya yazdırılmadı.** M1c/M2c/M3c
gerçek Z/Z' değerlerini hesapladı (KALEM'in açıkça izin verdiği gibi, "Engine
gates may report ONLY method DIFFERENCES and pass/fail counts, plus the
difference of k2/k3 between two methods on the check subsets") ama yalnız İKİ
YÖNTEM ARASINDAKİ farkları/geçti-kaldı sayılarını sakladı/yazdırdı — hiçbir mutlak
k2/k3 ya da Z/Z'/M dağılımı. M0c/M6c/M7c yalnız sıfırların KONUMLARINI (ve M6c
tamamen sentetik CUE Monte Carlo + gerçekleştirilebilir örnekleyiciler) kullandı.
`200c_olcum.npz` ve `HUKUM_200C.json` — gerçek ölçümün/analizin çıktıları — bu
görev kapsamında hiç üretilmedi; ölçüm betiği yalnız ONKAYIT_200C push'landıktan
SONRA çalıştırılacak.

## 9. Gerçek ölçüm için toplam çalışma zamanı tahmini

| bileşen | tahmini süre |
|---|---|
| Çapa kurma (4 pencere) | ~15s |
| b0 (1e6×4 pencere) | ~19.5+81+290+1900 ≈ 2290s (~38dk) |
| b1 (C1-C3 tümü + C4 alt-örnek 1e6) | ~29+122+435+1900 ≈ 2486s (~41dk) |
| b2 (olay tepe bulma, 129 ızgara, ~47k olay toplam) | birkaç dakika (M3c'nin 129-ızgara-yalnız kısmı çok daha hızlı, 513 yok) |
| **Toplam `200c_olcum.py`** | **kabaca 80-90 dakika** (C4 baskın; KALEM'in 15dk/hesaplama kuralı BU AJANIN KENDİ hesaplamalarına uygulanır, gerçek ölçüm ortak teftiş sonrası ayrı bir oturumda/parçalar halinde çalıştırılabilir) |
| `200c_analiz.py` | saniyeler (yalnız blok toplamları + ~47k olay üzerinde numpy) |
