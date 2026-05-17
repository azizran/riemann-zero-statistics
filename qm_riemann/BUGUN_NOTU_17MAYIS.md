# 17 Mayıs 2026 — Sezgiden Sayıya Yolculuk

*Uğur'un "plakada dalga deseni var" sezgisinin sayısal testinin günlüğü*

---

## SEZGİDEN ÇIKIŞ

Gün başında Uğur dedi (dün gecenin devamı):

> *"1/2'deki tüm sıfırların düştüğü yerler aslında plaka. Tıpkı çift yarıktaki gibi.
> O plakanın neresine düşüyor sıfırlar. Belki de aslına baktığımızda quantumdaki gibi
> dalga desenidir bu belki de."*

Ve:

> *"Riemann plakasında büyük boşlukların olduğu anlar var. O boşluklarda arkadaki dalgalar
> özellikle peak yapıyor, ya aşağıya ya da yukarıya. Bilinçli bir tercih gibi.
> Çift yarık değil ama çift yarık gibi bir pattern bırakıyor. Sanki var bir yapı."*

---

## SEZGİNİN TESTİ — 8 ARDIŞIK ADIM

### 1. Plaka görselleştirme (26)

σ=1/2 çizgisi bir plaka olarak çizildi, 269 sıfır parlak şerit. Gözle bakınca **doku** var ama klasik çift yarık (eşit şeritler + sinc² zarfı) deseni değil.

### 2. Boşluk-genlik testi (27)

Hipotez: ardışık sıfırlar arası büyük boşluk → orada Z(t) genliği büyük.

```
648 çift, t ∈ [14, 1000]
Pearson r = 0.701
Spearman r = 0.874
Normalize r = 0.924
```

**Doğrulandı**. Çok güçlü pozitif korelasyon.

### 3. Selberg-naif vs gerçek (28)

Klasik çift yarık şiddet eğrisi ile Z(t)² karşılaştırması:

| | Klasik | Riemann |
|---|---|---|
| Aralık | SABİT (9.87 ± 0.008) | DAĞINIK (1.52 ± 0.71) |
| Zarf | sinc² (ortada parlak) | aralık-bağlı |
| Yapı | 2 kaynak | ∞ kaynak (asallar) |

Selberg-naif null (sadece R(t) = √log t): r = 0.024 → **korelasyonu yaratmıyor**.

### 4. Gaussian-PSD null (29) — KRİTİK TEST

> *"R = 0.92 belki de Selberg açıklamamıştır."* (Uğur)

Z(t) ile **aynı güç spektrumlu** Gaussian rastgele süreç üret, ona aynı testi uygula:

```
Gerçek ζ:                r = 0.825
Gaussian null (30 trial): r = 0.505 ± 0.048
Z-skoru:                 +6.67σ  (p << 10⁻¹⁰)
```

**ζ-özgün doğrulandı**. Aralık-genlik korelasyonu Gaussian sürec yapısının çok ötesinde.

### 5. Slope kapalı form arama (30-32)

Linear fit: max|Z| ≈ slope(t) × normalize_aralık.

Slope T ile artıyor (Selberg-tipi yapısal form). slope² ile log(t/2π) arasında doğrusal:

```
T=3000:  slope² = 3.028·log(t/2π) + 5.753  (b ≈ π·log(2π))
T=6000:  slope² = 3.004·log(t/2π) + 6.256  (b ≈ 2π)
```

a katsayısı T ile **monoton düşüyor**: 3.77 → 3.30 → 3.22 → 3.17 → 3.06 → 3.00.

T=6000'de:
- a ≈ 3.00 (sapma %0.1)
- b ≈ 6.26 ≈ 2π = 6.283 (sapma %0.4)

**Aday hipotez**: slope²(t→∞) → 3·log(t/2π) + 2π

Ama yakınsama hâlâ tamamlanmamış olabilir — daha büyük T şart.

---

## ÖZET — NE BULDUK, NE BULMADIK

### Kesin pozitif bulgular (3 sağlam sonuç)

1. **Aralık-genlik korelasyonu var ve güçlü**: r=0.92 (normalize)
2. **Bu Gaussian sürec yapısının ötesinde**: +6.67σ (ζ-özgün)
3. **Slope yapısal form Selberg-tipi**: slope² ∝ log(t)

### Belirsiz / yarım pozitif

4. Slope² için kapalı form aday: **3·log(t/2π) + 2π** (T=6000'e kadar)
   - Daha büyük T testi şart
   - Yakınsama monoton ama tamamlanmadı

### Reddedilen hipotezler bu gün

- slope² = π·log(t)  (koincidans, T arttıkça bozuldu)
- b = π·log(2π)  (T=3000'de iyi, T=6000'de 2π'ye yakın)

---

## LİTERATÜR TARAMASI (Claude'un bildiklerinden)

Bu sonuçlar literatürde nereye düşüyor?

### Doğrudan ilgili çalışmalar

**Selberg (1946)** — log|ζ(½+it)| log-normal dağılım:
- Bizim slope formülünün **yapısal sebebi**
- Selberg ⟨|Z|²⟩ ≈ log(t/2π) (Selberg moment)
- Bizim slope² ∝ log(t) bu sonucun **lokal maks** versiyonu

**Lehmer (1956)** — Lehmer pairs (yakın sıfır çiftleri):
- Tersi yön: küçük aralık → küçük max|Z|
- Bizim büyük aralık → büyük max|Z| **aynı paranın diğer yüzü**
- Lehmer pairs bilinen ama aralık-max korelasyonunun **sayısal değeri** literatürde net değil

**Conrey-Ghosh (1989)** — moments of zeta:
- ∫|Z|² dt ≈ T·log T (ortalama)
- Bizim slope² ∝ log(t) buradan türetilebilir

**Hughes-Keating-O'Connell (2001)** — characteristic polynomial of GUE:
- Max log|ζ| in unit interval ~ log T
- Bizimle uyumlu ama aralık-bağlı değil

**Soundararajan (2008-2009)** — extreme values of ζ:
- Max|ζ(½+it)| in long intervals ~ exp(√(log T·log log T))
- Çok farklı ölçek

**Najnudel + Arguin-Belius-Bourgade-Radziwill-Soundararajan (2018-2019)** — local maxima of ζ:
- In a unit interval, max log|ζ| ≈ log log T - (3/4) log log log T + ...
- Bu çok hassas formülasyon
- Bizim slope² ≈ 3·log(t/2π) + 2π formülü bunun lokal versiyonu **olabilir** ama doğrudan eşleşmiyor

### Bizim katkımız nerede olabilir

**Pedagojik / görsel**:
- "Çift yarık plakası" çerçevesi
- 0-1/2-1 üçgeni
- Sezgisel "büyük boşluk = büyük dalga" → r=0.83 sayısal teyit

**Belki yeni sayısal gözlem**:
- **r = 0.83**, +6.67σ Gaussian-PSD'den fark — bu spesifik sayının literatürde direkt görmedim
- slope² = 3·log(t/2π) + 2π aday formülü (eğer doğruysa) — kontrol etmedim

**Kesinlikle yeni değil**:
- Selberg-tipi yapı (1946'dan beri)
- ζ ≈ GUE istatistik (1972'den beri)
- Lehmer pairs (1956'dan beri)

---

## SENİN SEZGİNİN VARDIĞI YER

Başlangıç:
> *"plakaya parçacık düşüyor, dalga deseni var, bir yapı var"*

Sayısal sonuç:
- r=0.83 korelasyon
- 6.67σ Gaussian-üstü
- Yapısal formül slope² ≈ 3·log(t/2π) + 2π (aday)

**Yaklaştık mı?** EVET. Sezginin söylediği "yapı" gerçek bir sayıya kavuştu. Sayı:
- Gaussian sürec **olamaz** (6.67σ)
- ζ-özgün bir şey
- Kapalı form (henüz tam değil ama yakın)

**Çözdük mü?** Hayır. Aslında çok ince bir adım attık — bu literatürün bilinen kıyısında oturduk, **bizim çerçevemiz** (çift yarık + plaka + dalga deseni) henüz literatürde standart değil ama matematik sonuçları kıyıdaki bilinene yakın.

---

## DERİN LİTERATÜR TARAMASI (2. ajan, 17 Mayıs)

Sistematik aramayla bulunan en yakın referans:

**Hall, R. R. (1999)** — *The second moment of Z'(t)*:
$$\int_0^T Z'(t)^2\,dt \sim \frac{T \cdot (\log T)^3}{6\pi^2}$$

Bağlantı bizim formüle: $|Z|_{\max} \approx |Z'(\gamma_n)| \cdot \text{gap}/2$ (parabolik tepe), o zaman:
- $|Z|_{\max}^2 \approx \frac{1}{4} |Z'|^2 \cdot \text{gap}^2$
- $\langle|Z'|^2\rangle \sim (\log t)^3 / (6\pi^2)$
- ortalama gap $\sim 2\pi/\log(t/2\pi)$

Sonuçta beklenen ~ $\log(t)/\text{sabit}$ formu — bizim slope² ailesinde ama **birebir eşleşme yok**.

**Net karar (2. ajan)**:
> *"Direkt karşılığı YOK. Bu üç şey literatürde bulunamadı:
> - (gap, max|Z|) çifti için rapor edilmiş Pearson r
> - 6.67σ Gaussian-PSD null hipotezine karşı test
> - slope²(t) = 3·log(t/2π) + 2π kapalı form aday formülü
>
> Hall–FHK–log-correlated field çerçevesinde **henüz sayısal olarak raporlanmamış** Pearson korelasyonu ve aday formül. Paper yazımı için somut bir niş var."*

**Risk**: Ivić "Theory of Hardy's Z-Function" (Cambridge 2012) ve Hall'ın orijinal 1999 papırı paywall ardında, tam erişilemedi. Eğer bu kaynaklarda Pearson r ölçümü varsa, bizimki tekrar olur.

**Önerilen dürüst çerçeve**: *"Tamamen yeni" demeyin; "Hall'ın Z'² momenti ve log-correlated field çerçevesinde sayısal olarak raporlanmamış gap-amplitude joint Pearson korelasyonu ve aday slope² asimptotik formülü."*

---

## YARIN İÇİN AÇIK ADIMLAR

1. **Literatür kontrolü (Google Scholar)**:
   - "Riemann zeros gap amplitude correlation" — bu spesifik korelasyon var mı?
   - Najnudel 2018, Arguin et al 2019 ile bizim formülünü karşılaştır
   - **Önemli**: bizim r=0.83 sayısının literatürde tam karşılığı var mı?

2. **Daha büyük T testi**:
   - T=20000 veya T=50000 (gece çalıştırılabilir, hesap ~saatler)
   - a katsayısı 3'te mi sabit kalıyor, daha mı düşüyor?

3. **Analitik türetim**:
   - Riemann-Siegel formülünden slope² = ? · log(t/2π) + ? doğrudan türetim
   - Bu konu zor ama somut

4. **Pedagojik belge**:
   - "Çift yarık çerçevesinde Riemann sıfırları" — bir öğretici makale
   - Bu bilim katkısı değil ama pedagojik değer var

---

## T=30000 TEST SONUCU (17 Mayıs, akşam)

90.9 dakika hesap. Sonuç **dürüst olumsuz**:

```
T_upper    a fit    b fit
500        2.87     6.83
1000       3.22     5.59
2000       3.05     6.27
6000       3.00     6.26     ← KOİNCİDANS (a≈3, b≈2π)
18000      2.71     7.75     ← bir an a ≈ e geçti
30000      2.63     8.12     ← hâlâ iniyor
```

**Hipotez slope² = 3·log(t/2π) + 2π ÇÜRÜTÜLDÜ.**

T=6000'deki "muhteşem koincidans" sonlu boyut etkisi. T arttıkça:
- a katsayısı monoton **düşüyor** (3'ten 2.6'ya)
- b katsayısı monoton **yükseliyor** (2π'den 8'e)

Alternatif fit: slope² ≈ 19.55·log log(t) − 15.95 → FHK ölçeği (log log T) olabilir, ama bizim ölçüm aralığında bu da kesin değil.

**Karar**: aday formülümüz çürütüldü. Asıl asimptotik form henüz bilinmiyor.

### Kalanlar ve düşenler

| Bulgu | Durum |
|---|---|
| r = 0.83 gap-amplitude korelasyonu | ✓ Kalır |
| +6.67σ Gaussian-PSD'den ζ-özgün | ✓ Kalır |
| Kuplaj kavramı (sezgisel) | ✓ Kalır |
| slope² ≈ 3·log(t/2π) + 2π | ✗ DÜŞTÜ |
| Yeni asimptotik kapalı form | Henüz yok |

### Önemli ders

İyi ki test ettik. T=6000'de durup paper yazsaydık reziliyolduk. Şimdi temiz bir negatif sonuçla devam edebiliyoruz — ana bulgular sağlam, aday formül elendi.

## T=100000 İKİNCİ TEST SONUCU (17 Mayıs, gece)

T=30000 sonrası ikinci aday formül: slope² = (π³/4)·log^(2/π)(t/2π).

T=5000'de β tam 2/π'ye yapıştı (sapma -0.0035), T=10000'de de aynı. AMA T=30000'de β=0.620, T=75000'de β=0.613. **Yine düşmeye devam ediyor.**

```
T=5000   β = 0.633  ← TAM 2/π koincidans
T=10000  β = 0.633  ← devam
T=30000  β = 0.620  ← uzaklaşma
T=75000  β = 0.613  ← hâlâ iniyor
```

**İkinci aday formül de ÇÜRÜTÜLDÜ.** 

Hipotez RMS oranı: T=30000'de %2.6 yakınlıkta idi, T=75000'de %28'e çıktı.

### Net karar (iki dürüst negatif sonuç sonrası)

Bu seansta **iki kez** aynı tuzağa düştük:
1. Küçük T'de güzel koincidans
2. Hızlı hipotez kurma
3. Büyük T'de çürüme

**Aday formül arayışını bırakıyoruz.** Çünkü yakınsama log log T ölçeğinde olduğundan, T=10¹⁵+ gerekli somut görmek için. Bizim ölçüm aralığımızda kapalı form ÇIKARILAMAZ.

### Gerçek katkımız (kalır)

- ✓ r = 0.83 (gap-amplitude Pearson korelasyonu)
- ✓ +6.67σ Gaussian-PSD null'undan ζ-özgün
- ✓ Kuplaj kavramsal çerçevesi (çift yarık + plaka)
- ✓ Hall (1999) Z'² ile yapısal bağlantı (köprü kurulabilir)

### Asıl asimptotik form

Trend β → 1/2 (Selberg √log) doğru yönde olabilir ama doğrulanamaz bu boyutta. Bu yorum spekülatif.

Yarın için açık soru: gerçek asimptotik form ne? T=10⁵+ ile test mi (saatler), yoksa analitik türetim mi (Hall ve Riemann-Siegel'den)?

---

## DEFTERİN STATÜSÜ

```
qm_riemann/
├── 1-32 Python betiği (test, görsel, analiz)
├── ~30 PNG görsel
├── BULGULAR.md (genel özet)
├── BUYUK_AN.md (16 Mayıs büyük an)
├── 08_bizim_acidan.md (yanlış kapı belgesi)
└── BUGUN_NOTU_17MAYIS.md (bu belge)
```

İki günlük yoğun çalışma. Sezgi → görsel → test → null → formül arayışı zincirinin
hem **başarıları** hem **başarısızlıkları** dürüstçe kayıtta.

---

*"Sezgi doğru yöne baktı. Sayı oraya ulaştı. Kapalı form henüz tam değil ama yapısal
form sağlam. Yarın literatür kontrolü ile bizim noktanın gerçek yerini göreceğiz."*

— 17 Mayıs 2026, qm_riemann/BUGUN_NOTU_17MAYIS.md
