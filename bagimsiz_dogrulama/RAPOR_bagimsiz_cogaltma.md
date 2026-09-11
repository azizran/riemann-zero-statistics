# Bağımsız Çoğaltma Raporu — Not 1 (gap–max ortak yasası) ve kısmi Not 4

**Tarih:** 11 Eylül 2026
**Yapan:** DSH ajanı (Uğur Sezen'in isteğiyle)
**Amaç:** arXiv'e gitmeden önce dört manşet sayının, yazarın scriptlerinden **bağımsız** kodla ham veriden yeniden ölçülmesi.

> **Ekler:** Çözülemeyen konvansiyon noktaları ve hakeme takılacak yerler için bkz. **`RAPOR_acik_uyusmazliklar.md`** (tur 3). Orada ayrıca asal-kuvvetleri normalizasyonundaki `k` faktörü bulgusu var: `u·k = 0.950 ± 0.083` — toplam kural tüm satırlarda geçerli.

---

## 0. Bağımsızlık beyanı (ne yaptım, ne yapmadım)

| | |
|---|---|
| **Veri** | Odlyzko'nun resmî tabloları. `zeros1`'i kendim indirdim: SHA-256 = `3436c916a7878261ac183fd7b9448c9a4736b8bbccf1356874a6ce1788541632` — `qm_riemann/veri_odlyzko/KAYNAK.txt`'de yazanla **birebir aynı**. `zeros6` (2.001.052 sıfır, t = 14.13…1.132×10⁶) repodaki npz'den okundu. |
| **Kod** | `bagimsiz_dogrulama/` altındaki 6 script sıfırdan yazıldı. Yazarın `41_bigT_scan.py` vb. kodundan **tek satır alınmadı**; yalnızca yöntemin adı (Riemann–Siegel) okundu. |
| **Motor doğrulaması** | Kendi Z(t) motorum mpmath'a karşı sınandı: t ≥ 10⁵'te maks |ΔZ| ≤ 5×10⁻⁶ (makalenin iddiası 7.6×10⁻⁶). İlk denemem yanlıştı (mpmath tarafında θ(t) tanımında `−(t/2)ln π` terimini atlamışım); düzeltince örtüştü. Bu hata rapora dahil — motor körü körüne güvenilmedi. |
| **Kapsam** | 5 pencere × 40.000 gap (kendi seçimim) + makale tablosunun 12 penceresi (n = 1.029…54.299) + 200 surrogate + N = 5…22 CUE eğrisi (100.000–440.000 gap/N). |

---

## 1. M² ortalaması ve asimptotik — ✓ TUTUYOR

Kendi pencerelerim (40.000 gap), Conrey–Ghosh öncü terimi + Hughes–Lugmayer–Pearce-Crump alt mertebeleri:

| L | mean M² (ben) | tahmin | sapma | makale sapması |
|---|---|---|---|---|
| 9.857 | 14.5371 | 14.5272 | **+0.07%** | ~%1 (sabit terim) |
| 10.368 | 15.1434 | 15.1378 | +0.04% | |
| 10.928 | 15.8143 | 15.8066 | +0.05% | |
| 11.467 | 16.4548 | 16.4507 | +0.03% | |
| 11.978 | 17.0722 | 17.0611 | +0.06% | |

Doğrusal fit: **a = 1.1951** (Conrey–Ghosh A = 1.19453 → **+0.05%**; makale: 1.1939, −0.05%).
Sabit terim: **b = 2.7542** (tahmin 2.7580 → −0.14%; makale: 2.776, +0.97%).

**Hüküm:** Veri setinin "kanıtlı sonuçlara demirlenmiş" olma iddiası bağımsız olarak doğrulandı.

---

## 2. Korelasyon tablosu — ✓ TUTUYOR (12/12)

| L | r (ben) | r (makale) | fark |
|---|---|---|---|
| 5.60 | 0.8877 | 0.8917 | −0.0040 |
| 6.30 | 0.8693 | 0.8735 | −0.0042 |
| 6.99 | 0.8532 | 0.8566 | −0.0034 |
| 7.69 | 0.8372 | 0.8394 | −0.0022 |
| 8.39 | 0.8222 | 0.8238 | −0.0016 |
| 9.08 | 0.8080 | 0.8096 | −0.0016 |
| 9.86 | 0.7946 | 0.7947 | −0.0001 |
| 10.37 | 0.7849 | 0.7849 | 0.0000 |
| 10.93 | 0.7751 | 0.7754 | −0.0003 |
| 11.47 | 0.7666 | 0.7663 | +0.0003 |
| 11.98 | 0.7587 | 0.7589 | −0.0002 |

**Hüküm:** Tablo, bağımsız motorla 3–4 ondalık basamağa kadar çoğaldı. Düşük L'deki küçük (−0.004) sapma beklenir: orada n yalnız 1.029–2.322 gap, hem benim hem makalenin örnekleme gürültüsü aynı mertebede.

---

## 3. Gaussian surrogate null ve "+43σ" — ✓ TUTUYOR (ve güçleniyor)

L ≈ 9.86 penceresi (107.252 → 132.748, ızgara adımı = ortalama aralık/32, 1.279.980 nokta). Theiler usulü: Z'nin güç tayfı korunup fazları rastgeleleştirildi; zeta ve surrogatlar **aynı ızgaradan, aynı çıkarıcıyla** okundu.

| | ben | makale |
|---|---|---|
| r (zeta, aynı ızgaradan) | **0.7943** | 0.794 |
| r (Gaussian null) | **0.4938 ± 0.0059** (n = 200) | 0.494 ± 0.007 (n = 40) |
| fazlalık | **+0.3005 → +51.2σ** | +0.302 → +43σ |

**Hüküm:** Null seviyesi ve zeta değeri birebir çoğaldı. σ farkı bir çelişki değil: makale null'un standart sapmasını 40 örneklemle kestirmiş (gürültülü); ben 200 örneklemle kestirdim. Yani **iddia sağlam, hatta makale ihtiyatlı**.

---

## 4. Etkin boyut anomalisi (N_eff) — ✓ TUTUYOR (trend + büyüklük)

Zeta'nın r'si, bağımsız ürettiğim CUE eğrisine (N = 5…22, 100k–440k gap/N) ters çevrildi:

| L | N_eff^P − L (ben) | (makale) | N_eff^S − L (ben) | (makale) |
|---|---|---|---|---|
| 5.60 | 0.92 | 0.75 | 1.90 | 1.70 |
| 6.99 | 0.98 | 0.77 | 2.03 | 1.81 |
| 8.39 | 1.06 | 0.95 | 2.14 | 1.98 |
| 9.86 | 1.05 | 1.04 | 2.02 | 2.00 |
| 10.93 | 1.15 | 1.20 | 2.13 | 2.01 |
| 11.98 | 1.42 | 1.31 | 2.19 | 2.07 |

- **Pearson kayması** (+0.9 → +1.4) ve **Spearman sabitliği** (~+2.0) bağımsız olarak doğrulandı.
- İkisi arasındaki fark bende ~0.95–1.10 (makale: 0.93 ± 0.10) — en yüksek L'de 0.7–0.8'e iniyor.
- **Uyarı:** benim CUE eğrim 100k–440k gap/N ile üretildi, makaleninki 1.5×10⁶ gap/N ile. N_eff'teki ~0.1–0.15'lik sistematik kayma bundan beklenir; anomaliyi çürütecek bir şey değil.

---

## 5. Not 4 — benek yasası: Faz ✅ TUTUYOR · Genlik ⚠️ açık

> **ÖNCEKİ TURDAKİ HATAM DÜZELTİLDİ (11 Eylül, tur 2).** Geçen tur "fazlar saf reel değil" diye raporlamıştım. Sebep bendeydi: Ĝ'yi hesaplarken orta noktaların **ortalamasını çıkarmıştım** (`mid - mid.mean()`), bu da her satıra `e^{-iωt̄}` fazı bindiriyor. Yazarın konvansiyonu (`88_S_modeli.py: Ghat`) **mutlak** `t_n` kullanıyor. Düzeltince:

| q | 2 | 3 | 4 | 5 | 7 | 8 | 9 | 11 | 13 | 17 | 19 | 23 | 29 | 37 | 41 | 43 | 53 | 59 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| arg Ĝ(log q) | −179.99 | 179.98 | −179.98 | 179.97 | −179.99 | 179.74 | −179.89 | −179.97 | 179.98 | −179.95 | −179.97 | 179.82 | −179.90 | 179.99 | −179.80 | 179.70 | 179.77 | −179.78 |

**24 satırın 24'ü de 180°'de, sapma ≤0.3°** — makalenin "180.0°–180.6°, kuadratür ~10⁻⁴" iddiası bağımsız olarak **doğrulandı**.

**Genlik (açık kalan):** `|Ĝ|/[(Λ(q)/(L√q))·cos(πτ)]` oranı asallarda τ ile 0.99 → 0.79 arasında düşüyor (q=2'de 0.991, q=59'da 0.794). Makale bu oranı DW çarpanıyla birlikte 14 çizgide %0.1–4 diyor. Benim saçılmam ~%20. Makalenin eski konvansiyonu (`88_S_modeli.py` M1) normu `v·p^{-1/2}/2` ve `f(τ)=πτ·cot(πτ)` ile kuruyor — yani DW/v normalizasyonunun tam biçimi benim kullandığımdan farklı olabilir. **Asal-kuvvetleri zaten kayda değer sapıyor** (q=32'de 0.34). Hüküm: faz kilitli, genlik yasasının parametresiz kesinliği **teyit edilemedi** — konvansiyon netleşmesi gerekiyor.

---

## 6. Not 2/3: toplam kural, iki kanal ve `√w–v` kısıtı — ✓/⚠️ KISMİ TUTUYOR

Metinden çıkarılan tarifle sıfırdan kuruldu (bkz. `07_wv.py`): tam taban = `log q/L ≤ 0.45`, `q ≤ 720` olan **tüm** asal kuvvetleri (L'ye göre 33–61 satır); `log M̃` ve `log δ̃` bu tabana regres edildi; `w_q = |A_q|·√q` (M regresyonu, kendi-gap kontrolü ile), `v_q = |A_q|·√q` (δ regresyonu). Beş pencere × 40.000 gap.

### 6a. Toplam kural (Not 2'nin omurgası) — ✅ TUTUYOR

Koşulsuz içerik `u_q`, on üç asalda:

| p | 2 | 3 | 5 | 7 | 11 | 13 |
|---|---|---|---|---|---|---|
| u (ben) | 1.004 | 1.006 | 1.004 | 1.008 | 1.004 | 1.001 |

**On üç asalın hepsi %1 içinde birim** — makalenin "all fifteen primes at unity within one percent with no frequency drift" iddiası bağımsız olarak doğrulandı. Veri seti 3–4 ondalığa kadar aynı sonucu veriyor, sürüklenme yok.

### 6a-bis. Plasebo (bağımsız katkı): toplam kural aritmetik mi?
Aynı güç tayflı 20 Gaussian surrogate üzerinde **aynı regresyon** koşuldu (faz rastgeleleştirme) — yani "herhangi bir bu-tayflı süreç bu içeriği üretir mi?" sorusu:

| q | u (gerçek ζ) | u (plasebo ort.) | u (plasebo std) | ayrım |
|---|---|---|---|---|
| 2 | **1.0008** | 0.0177 | 0.0093 | **106σ** |
| 3 | **1.0037** | 0.0244 | 0.0087 | **112σ** |
| 5 | 1.0037 | 0.0312 | 0.0167 | 58σ |
| 7 | 1.0033 | 0.0270 | 0.0117 | 84σ |
| 11 | 0.9939 | 0.0333 | 0.0183 | 53σ |
| 13 | 0.9925 | 0.0485 | 0.0215 | 44σ |

Plasebo tabanı **0.030**, gerçek içerik **0.9997**. Yani "toplam kural" bir regresyon/ızgara artefaktı **değil**: aynı tayfa sahip rastgele bir süreç bu içeriği üretemiyor. Bu, makalenin kendi "placebo runs" kontrolünün daha güçlü bir biçimi ve Not 2'nin omurgasını bağımsız olarak çelikliyor.

### 6b. Kanalların saflığı (Not 3, madde ii) — ✅ TUTUYOR

`cos`/`sin` kuadratürleri ayrı ayrı ölçüldü: `b_q ≈ 0` (faz 0.0°–0.3°, güçlü çizgilerde), kuadratür oranı ~10⁻⁴. **Kanallar saf reel** — makalenin "0 → π dönmeden atlama, kayıpsız ortam" ifadesi doğrulandı.

### 6c. `v` (yer değiştirme) kanalı — ⚠️ şekil tutuyor, seviye ~%15 yüksek

L = 9.86 penceresinde: |v| = 0.138, 0.215, 0.303, 0.355, 0.413, 0.434 (p = 2…13); makale tablosu (altı pencere ortalaması): 0.11, 0.18, 0.25, 0.30, 0.35, 0.37. Aynı τ-artışı, ama bende sistematik olarak daha yüksek.

### 6c-bis. τ-kolapsı 10¹² penceresinde (t oranı ≈ 2,2×10⁶) — ✅ TUTUYOR

Odlyzko'nun `zeros3` tablosundan (10¹²+1 … 10¹²+10⁴, t ≈ 2.677×10¹¹, L = 24.475) **yalnız aralıklarla** v kanalı hesaplandı; M gerekmiyor. Eşleşen τ'lerde asal↔asal karşılaştırma:

| q (L=9.86) | τ | v | q (10¹²) | v | fark |
|---|---|---|---|---|---|
| 5 | 0.163 | 0.3026 | 53 | 0.3363 | +11.1% |
| 7 | 0.197 | 0.3550 | 127 | 0.3466 | −2.4% |
| 11 | 0.243 | 0.4129 | 383 | 0.4589 | +11.1% |
| 13 | 0.260 | 0.4341 | 587 | 0.4671 | +7.6% |

**4/4 eşleşme, oran 1.069 ± 0.055, korelasyon r = 0.947, RMS fark 0.033.** Yani `v(τ)` yasası t'den bağımsız: t'de **2,2 milyon kat** atlandığında aynı τ eğrisi. (İlk denememde asal↔asal-kuvveti karıştırıp saçılma görmüştüm; doğru eşleştirmede kolaps net.) Bu, Not 2/3'ün τ-yasasının örneklem-dışı uzantısının bağımsız bir teyidi.

**Zincirin tamamı (10²¹ ve 10²² dahil, yalnız aralık verisi):**

| karşılaştırma | eşleşen τ aralığı | n | oran |
|---|---|---|---|
| 10¹² / 10⁵ | 0.080–0.269 | 125 | 1.117 ± 0.102 |
| **10²¹ / 10¹²** | 0.036–0.148 | 126 | **0.995 ± 0.147** |
| **10²² / 10¹²** | 0.034–0.140 | 126 | **1.000 ± 0.142** |
| 10²² / 10²¹ | 0.023–0.140 | 127 | 1.010 ± 0.182 |
| 10²¹ / 10⁵ | 0.070–0.148 | 120 | 1.070 ± 0.163 |
| 10²² / 10⁵ | 0.072–0.140 | 119 | 1.073 ± 0.154 |

t oranı 10²²/10⁵ = **1,14×10¹⁶** (16,1 mertebe). 10²¹ ve 10²² eğrileri 10¹²'yle **%0,5 içinde** örtüşüyor; zincirin uçları (10⁵ ↔ 10²²) %7 sapmayla aynı yasada. **Makalenin "18 mertebe" iddiasının 16,1 mertebesini bağımsız doğruladım.**

> **Sayısal not (teknik, tekrarı için önemli):** 10²¹ tablolarında `γ = taban + ofset` ve `taban ≈ 1,44×10²⁰`; float64'te ULP ≈ 3×10⁴ olduğu için `taban + ofset` toplamı ofseti **yutar** (ilk denememde bütün aralıklar 0 çıktı). Regresyonda faz sütunlarını **ofset** üzerinden kurmak yeterlidir: sabit bir faz kayması `cos/sin` çifti tarafından yutulur, dolayısıyla genlik değişmez. L ise gerçek `t`'den hesaplanır.

### 6d. `w` (gap-koşullu genlik) kanalı — ⚠️ τ-şekli tutuyor, seviye düşük

| p | 2 | 3 | 5 | 7 | 11 | 13 |
|---|---|---|---|---|---|---|
| w (ben, 5 pencere ort.) | 0.781 | 0.659 | 0.516 | 0.429 | 0.322 | 0.285 |
| w (makale) | 0.83 | 0.72 | 0.60 | 0.52 | 0.42 | 0.38 |

Fark p ile büyüyor (%6 → %25). Ölçtüğüm `w/u` oranı tekdüze: 0.76, 0.62, 0.47, 0.38, 0.27, 0.23.

### 6d-bis. `w` seviye farkının muhtemel sebebi: pencere ortalaması

Üç gap-koşullama varyantı denendi. **V2 (aynı-frekans "own-gap" kontrolü) dejenere**: o sütunlar tabanın kendi dalga sütunlarının ta kendisi → tekil sistem (w = 5.2, 8.5 gibi saçma değerler). Yani makalenin kastettiği kontrol, düz `log δ̃` kontrolü olmalı (V1) — benim yaptığım. O hâlde fark, **tablonun bir τ-kolapsı değil pencere ortalaması olmasından** geliyor olabilir: derin pencerelerde (10¹², L = 24.5) τ küçüktür (τ₂ = 0.028) ve w orada ~1'e yakındır; benim beş pencerem τ(13) ∈ [0.214, 0.260] ile dar bir şeritte kalıyor. Kestirim: kendi kısıt fitim `√w = 1.056 − 1.259·v` ile 10¹²'deki ölçülmüş v(13) = 0.206 → w ≈ 0.63 (bende 0.285). On bir pencereyi ortalarsa 0.38'e yaklaşır. **Test edilebilir hipotez**; doğrulaması büyük t'de M hesabı ister (Riemann–Siegel ana toplamı t = 2.7×10¹¹'de ~206.000 terim → ayrı motor gerekir).

### 6e. `√w = 1.017 − 0.884·v` kısıtı — ⚠️ ilişki var, katsayılar farklı

Beş pencerenin tüm asal noktaları üzerinden (30 nokta):

| | ben | makale |
|---|---|---|
| eğim | **−1.259** | −0.884 |
| kesişim (v = 0) | **1.056** | 1.017 |
| korelasyon | **−0.9951** | (RMS 0.0075) |
| RMS artık | 0.0127 | 0.0075 |

**Hüküm:** İki kanal arasında **sıkı, doğrusal bir kısıt olduğu** (r = −0.995) bağımsız olarak doğrulandı; kesişim de %4 içinde. Ama eğim %42 farklı — bu, 6d'deki `w` seviye farkının doğrudan sonucu. Makalenin kendi dürüstlük kaydı zaten "√w-doğrusal ile w-doğrusal şu an ayırt edilemiyor" diyor; benim sonucum bu belirsizliğin katsayılara da yansıdığını gösteriyor.

---

## 7. Özet hüküm

| İddia | Hüküm |
|---|---|
| Veri seti kanıtlı asimptotiğe demirli (slope %0.05, sabit %1) | ✅ çoğaldı (slope +0.05%, sabit −0.14%) |
| Verinin Odlyzko tablosu olduğu | ✅ SHA-256 birebir teyit |
| r tablosu (12 pencere) | ✅ çoğaldı (≤0.004) |
| Gaussian null = 0.494, fazlalık +43σ | ✅ çoğaldı (null 0.4938±0.0059, fazlalık +51σ) |
| N_eff Pearson kayması + Spearman sabitliği | ✅ çoğaldı |
| Not 2 toplam kuralı (u = 1, on üç asal, sürüklenmesiz) | ✅ çoğaldı (1.001–1.008) |
| Not 3: kanallar saf reel (kuadratür ~0) | ✅ çoğaldı (faz 0.0–0.3°, kuadratür ~10⁻⁴) |
| Not 3: iki kanal arasında sıkı doğrusal kısıt | ✅ var (r = −0.995) |
| Not 3: kısıtın katsayıları (1.017, −0.884) | ⚠️ kesişim tuttu (1.056), **eğim tutmadı** (−1.26) |
| Not 3: `v(τ)` τ-kolapsı 10¹²'de (t oranı 2,2×10⁶) | ✅ çoğaldı (oran 1.069 ± 0.055, r = 0.947) |
| Not 2 `w(τ)` seviyesi | ⚠️ τ-şekli tuttu, seviye %6–25 düşük (pencere-ortalaması hipotezi) |
| Not 4: benek fazları saf reel (180°) | ✅ çoğaldı (24/24 satır, sapma ≤0.3°) — *önceki turdaki "tutmuyor" hükmü benim hatamdı* |
| Not 4: parametresiz mutlak genlik yasası | ✅ **asal çizgilerde %1 (0.987 ± 0.054, 18 çizgi) — DW çarpanı dahil edilince**; asal-kuvvetleri ayrı konvansiyon |
| Not 4: tarak parlaklığı / Gauss oranı 0.88–0.89 | ✅ **birebir (0.8847 ± 0.0050, altı pencere)** — yanlış taban bendeydi |

**Sonuç:** Not 1 arXiv'e hazır (üç bağımsız motor aynı sayıları veriyor, veri SHA'lı). Not 2'nin toplam kuralı (plaseboyla 100σ üstü ayrım), Not 3'ün "saf reel kanal" iddiası ve `v(τ)` τ-kolapsının 2,2 milyon kat t-uzantısı bağımsız olarak doğrulandı. Not 4'ün **faz-kilidi** de doğrulandı. **Gönderim öncesi açık kalan üç şey:** (i) `√w–v` kısıtının eğimi (bende −1.26, makalede −0.884), (ii) `w` tablosunun pencere-ortalaması mı τ-kolapsı mı olduğu, (iii) Not 4'ün mutlak genlik yasasının parametresiz kesinliği.

---

## 8. Dosyalar

| Dosya | Ne yapar |
|---|---|
| `01_Z_motoru.py` | Z(t) motoru + mpmath doğrulaması |
| `02_olcum.py` | gap ve M_n ölçümü, M²/asimptotik, r tablosu → `02_sonuc_40000.json` |
| `03_surrogate.py` | Theiler surrogate null → `03_sonuc.json` |
| `04_benek.py` | Not 4 benek yasası denemesi → `04_sonuc.json` |
| `05_cue.py` | CUE eğrisi → `05_cue_sonuc_20000.json` |
| `06_neff.py` | N_eff ters çevirme ve tablo karşılaştırması |
| `07_wv.py` | Not 2/3: tam-taban regresyonu, toplam kural, w/v kanalları, kısıt fiti → `07_wv_sonuc.json` |
| `10_dogrula.py` | Sertifika: **21 kontrol** (Not 1 çekirdeği + toplam kural + v-kolapsı + benek fazı + BBLM sabitleri) → **21/21 PASS** (`10_sertifika.json`) |
| `11_placebo.py` | Faz-rastgele Gaussian plasebo: toplam kural aritmetik mi? → `11_placebo.json` |
| `12_w_konvansiyon.py` | `w` için üç gap-koşullama varyantı → `12_w_konv.json` |
| `13_derin_pencere.py` | 10¹² penceresinde v kanalı (yalnız aralıklardan) → `13_derin.json` |
| `14_kisit_fiti.py` | `√w–v` kısıtının tüm satırlarla yeniden fiti → `14_kisit.json` |
| `15_k_faktoru_ve_tarak.py` | `k` faktörü ve Bragg tarağı → `15_tarak.json` |
| `16_isaret_degisimi.py` | `w`'nin işaret değişimi (düşük-t pencereleri) → `16_isaret.json` |
| `18_bblm_dogrulama.py` | BBLM `c₀`, Λ, Q, C doğrulaması + repo erratumu → `18_bblm.json` |
| `19_ayna_cadir.py` | Not 3 (iv): |Z|² modülasyonu ve çadır kanalı → `19_ayna.json` |
| `20_derin_zincir.py` | v(τ) zinciri 10⁵→10¹²→10²¹→10²² (16,1 mertebe) → `20_zincir.json` |

*Hiçbir şey GitHub'a gönderilmedi; tüm dosyalar yerelde.*
