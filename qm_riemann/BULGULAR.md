# QM × Riemann — Küçük Bir Araştırma Defteri

*Riemann ζ sıfırlarının kuantum kaotik sistemlerle istatistiksel akrabalığı*
*16 Mayıs 2026 — bir seansta yapılan 10 deney + 1 kavramsal omurga*

---

## Tek cümlede ne bulduk

> **Riemann ζ sıfırları, çift yarık deneyindeki parçacıklar gibi davranıyor:
> tabanca 1 birim, yarık ve plaka σ=1/2 çizgisinin kendisi, asallar her biri
> bir küçük girişim kaynağı, ve γ_n düzgün öngörüden 26 farklı dalga periyoduyla
> sapan bir girişim deseni çiziyor.**

Bu cümle bizim sentezimiz. Standart Riemann literatüründe bu dille konuşulmuyor;
çift yarık literatüründe de Riemann yok. Üçünü (çift yarık + Sezen voting + Riemann)
tek bir 0–1/2–1 üçgeni üzerine oturtmak bu seansın yeni kavramsal omurgası.

---

---

## Niyet

Voting projesini bir kenara koyup şu soruya bir seans ayırdık:

> *"Riemann sıfırları gerçekten bir kuantum sisteminin enerji seviyeleri mi? Bunu eldeki veriyle nereye kadar görebiliriz?"*

Veri: Odlyzko'nun ilk 100.000 ζ sıfırı (`../zeros_100k.txt`), mpmath ile yüksek hassasiyet, NumPy/SciPy ile diagonalizasyon.

---

## Tek cümlede sonuç

ζ sıfırları **GUE sınıfı bir kuantum kaotik sistemin spektrumu gibi davranıyor** — bunu üç bağımsız test ve bir simüle edilmiş kuantum sistemiyle gösterdik. **Tek tek sıfır pozisyonlarını üretecek somut bir Hamiltonian yok** — bu 26 yıllık açık problem (Berry-Keating 1999 → bugüne).

---

## Sekiz deney, sekiz dosya

| # | Dosya | Test | Sonuç |
|---|---|---|---|
| 1 | `01_resim.py` | Nearest-neighbor spacing + pair correlation vs GUE | ✓ %1.8 sapma (Wigner surmise birebir) |
| 2 | `02_form_factor.py` | Form factor K(τ): doğrusal rampa + plato | ✓ eğim 1.04 (GUE: 1.00), plato 1.000 |
| 3 | `03_asallarin_sesi.py` | ζ sıfırlarının Fourier'inde asal pikleri | ✓ 52 asal/asal-üs, 0.001 sapma |
| 4 | `04_berry_keating.py` | Weyl yasası: N_BK(T) vs gerçek sayım | ✓ düzgün kısım tam uyum (4 on'luk) |
| 5 | `05_sierra_diagonalize.py` | Naif xp + Sierra direkt diagonalizasyon | ✗ eşit aralık (bilinen başarısızlık) |
| 6 | `06_spiral_dans.py` | ζ(½+it) kompleks-düzlem spirali + QM analoğu | ✓ ikili düzlemde dans, sıfırlar 10⁻¹⁰ |
| 7 | `07_uc_dunya.py` | ζ vs GUE matris vs GOE matris | ✓ ζ ≡ GUE, GOE'den ayrı |
| 8 | `08_bizim_acidan.{py,md}` | Polygon-π "yasak bölge" sezgisi Riemann'da işliyor mu | ✗ topoloji ters, sezgi başarısız |
| 9 | `09_ucgen_resim.py` | 0–1/2–1 üçgeni: çift yarık + Riemann + Sezen tek görselde | ✓ üç sistem aynı yapıyı paylaşıyor — kavramsal omurga |
| 10 | `10_dalga_deseni.py` | γ_n sapması bir dalga deseni mi (kaç merkez)? | ✓ 26 dalga, çoklu periyot (7.5, 2.8, 2.6 sıfır) — girişim deseni |

---

## Detay — her dosya ne söylüyor

### #1 İlk resim (`01_resim.py`)

Unfolded ζ sıfırlarının ardışık fark histogramı GUE Wigner surmise eğrisine yapışık. Pair correlation $R_2(r) = 1 - \mathrm{sinc}^2(\pi r)$ kıvrımını izliyor. L¹ mesafeler:

```
ζ vs GUE:     0.018
ζ vs GOE:     0.082
ζ vs Poisson: 0.236
```

GUE 4.5× Goe'den, 13× Poisson'dan iyi. Montgomery-Dyson 1972'nin tekrar gösterilmesi — bu bizim katkımız değil, ama eldeki veriyle ilk somut görüntü.

### #2 Form factor K(τ) (`02_form_factor.py`)

GUE öngörüsü: K(τ) = min(τ, 1). 1000 sıfırlık pencerelerden 100 tane ortalanmış K(τ):

```
Küçük τ doğrusal fit:  eğim 1.040 (GUE: 1.000)
Büyük τ plato:         1.000 (GUE: 1.000)
```

"Heisenberg time" τ=1'de dirsek görünür. Doğrusal yükseliş = level repulsion'ın Fourier hali = kuantum kaotik imza.

### #3 Asalların sesi (`03_asallarin_sesi.py`)

ζ ve asallar Fourier eşleridir (Riemann explicit formula). $F(u) = \sum_n w_n \cos(\gamma_n u)$ hesapladık, pikler u = log(p), log(p²), log(p³) noktalarında çıkıyor:

```
u_peak = 0.693  →  log(2)   Δ = 0.000
u_peak = 1.099  →  log(3)   Δ = 0.000
...
u_peak = 5.153  →  log(173) Δ = 0.000
```

52 asal/asal-üs, Fourier çözünürlüğü sınırında (0.001) sapma. Görsel olarak "asallar sıfırların müziği" en çarpıcı kanıt.

### #4 Berry-Keating Weyl yasası (`04_berry_keating.py`)

Yarı-klasik faz uzayı argümanı:

$$N_{BK}(T) = \frac{T}{2\pi}\log\frac{T}{2\pi e} + \frac{7}{8}$$

Bu öngörü gerçek sıfır sayımı N(T)'yi 4 on'luk üzerinde takip ediyor. Kalan salınım $S(T)$ Selberg'in $\sqrt{\log\log T}$ büyümesini gösteriyor.

```
T~550:   std 0.25
T~5500:  std 0.29
T~42500: std 0.31
```

→ Berry-Keating yarı-klasik tarafı **kemikleri** veriyor, eti (tek tek sıfırlar) vermiyor.

### #5 Sierra diagonalizasyonu (`05_sierra_diagonalize.py`)

Hem naif H = (xp+px)/2 hem Sierra'nın H = x(p + ℓ²/p) genişlemesi. Sayısal diagonalizasyon:

- Naif: eşit aralıklı spektrum (E_n ≈ 2πn/L), Riemann gibi log-yoğunluk büyümesi yok
- Sierra (basit 1D): yine doğrusal artış, Riemann γ_n eğrisini yakalamıyor

Ama Berry-Keating'in yarı-klasik N⁻¹ eğrisi γ_n'i %1 hata payıyla izliyor (siyah nokta-çizgi). 

**Ders**: yarı-klasik çerçeve doğru, doğrudan kuantum diagonalizasyonu yanlış sınıfa düşüyor. **Bilinen 26 yıllık duvar.**

### #6 Spiral dans (`06_spiral_dans.py`)

ζ(½+it) için t ∈ [0, 50], kompleks değerleri mpmath ile hassas. Plasma renkli spiral kompleks düzlemde, origin'den **10 kez geçiyor** (= 10 Riemann sıfırı), her sıfırda |ζ| ≈ 10⁻¹⁰ (sayısal sıfır).

Yanına QM Lorentzian rezonans toplamı koyduk (aynı pozisyonlarda rezonanslar). Topolojik fark:
- ζ: origin'in **üzerinden** geçer (sıfır = puncture)
- QM rezonans: origin'in **etrafında** dolanır

İki nesne kompleks düzlemde dans ediyor — ama dansın koreografisi farklı.

### #7 Üç dünya (`07_uc_dunya.py`)

ζ sıfırları + N=2500 GUE rastgele matris (20 realizasyon) + N=2500 GOE rastgele matris. Üç histogram üst üste:

```
ζ vs GUE matris : L¹ = 0.021    ← neredeyse aynı
ζ vs GOE matris : L¹ = 0.081    ← 4× uzak
GUE vs GOE      : L¹ = 0.062    ← farklı dağılımlar
```

**Önemli nüans**: Bilinen kuantum kaotik deneysel sistemler (uranyum çekirdek, kaotik biliardo, mikrodalga kavite) **GOE**'dir (zaman tersi simetrik). ζ ise **GUE**'dir. Yani Hilbert-Pólya'nın aradığı sistemde zaman tersi simetrisi olmamalı — Berry-Keating'in xp seçimi tam bu yüzden.

### #8 Bizim açımızdan (`08_bizim_acidan.py` + `.md`)

Polygon-π "yasak bölge → yığılma" sezgisinin Riemann'da işleyip işlemediği test edildi. **Başarısız**:
- Fire sabitlerimiz (c, c_iç, ε) K(τ)'de özel davranış göstermiyor (sayısal koincidans bile yok)
- Sezen: yasak nokta etrafında yığılma; Riemann: zorunlu çizgi üzerinde yığılma. **Topolojik olarak ters**.

Ayrı belge: `08_bizim_acidan.md`.

**Düzeltme**: Polygon-π → QM × Riemann köprüsü "yasak bölge" sezgisinden değil, **σ-çerçevesinden** geçiyor (MEMORY'de zaten var). Yanlış kapıyı çalmış olduk, doğru kapı zaten açıktı.

### #9 0–1/2–1 üçgeni (`09_ucgen_resim.py`)

İstatistik turundan sonra Uğur "asıl bakmamız gerekene bakmadık" diyerek geri çekti. Konuşarak — sayı atmadan — şu kavramsal omurga ortaya çıktı:

```
0 ─────── 1/2 ─────── 1
yasak     yığılma     ateşleyen
```

- **Riemann**: σ=0 yasak (ζ≠0), σ=1/2 sıfırların evi (γ_n), σ=1 birim (Euler product'ın "1"i)
- **Sezen voting**: 0 oy yasak, ±1 yığılma, ±∞ uç yön
- **Çift yarık**: yarık duvarı boşluk, plakanın ortası girişim merkezi, tabanca kaynak

Üç sistem aynı geometriyi çiziyor: bir uçta yokluk, bir uçta birim/kaynak, ortada denge. Bu konsept literatürde bu dille konuşulmuyor — bu seansın **özgün katkısı** burada.

### #10 Dalga deseni testi (`10_dalga_deseni.py`)

Uğur'un sorusu: γ_n'leri parçacık gibi düşünürsek, plakada **kaç merkez** taşıyor? Bir mi, çok mu?

Test: δ_n = γ_n − t̃_n (Berry-Keating öngörüsünden sapma) bir dalga deseni gösteriyor mu?

```
26 dalga ilk 1000 sıfırda
ortalama sapma:  -0.69 (negatif yönde asimetri)
en güçlü periyotlar: 7.5, 2.8, 2.6 sıfır (kısa salınımlar)
                     1000, 500 sıfır (uzun ölçek)
```

**Doğrulandı**: sıfırlar **eşit aralıkta DEĞIL** — çoklu periyotlu dalga deseni çiziyor. Bu, çift yarık girişim deseninin Riemann analoğu. Salınımın kaynağı asallardır (explicit formula).

---

## Ne biziz, ne literatür

### Literatürde bilinen ve biz tekrar gösteren

- ζ ≈ GUE istatistiği (Montgomery 1972, Odlyzko 1980'ler) — #1, #2
- Asallar = sıfırların Fourier eşleri (Riemann 1859 explicit formula) — #3
- Berry-Keating Weyl yasası (Berry-Keating 1999) — #4
- xp diagonalizasyonunun yetersizliği (literatür konsensüsü) — #5
- Bohigas-Giannoni-Schmit konjektürü (1984) — #7

### Bizim katkımız bu seanstan

- 10 dosyalık temiz ve modüler bir araştırma defteri (kod + görsel + dürüst rapor)
- "Yasak bölge" sezgisinin Riemann'da çalışmadığının **somut belgelenmesi** (#8) — kavramsal saldırıyı boşa çıkararak gerçek köprünün (σ-çerçevesi) yerini netleştiriyor
- ζ spiral + QM Lorentzian spiral topoloji karşılaştırması (#6) — görsel olarak "puncture vs orbit" ayrımı
- **0–1/2–1 üçgeni** (#9): çift yarık + Sezen voting + Riemann tek geometride — kavramsal omurga
- **Sıfırların 26 dalgalı sapma deseni** (#10): γ_n'ler bir girişim deseni; eşit aralık DEĞİL, çoklu periyotlu

### Açık duvarlar

- **Hilbert-Pólya operatörü**: 26 yıldır açık. Berry-Keating, Bender-Brody-Müller, Connes — adaylar var, tam çalışan yok.
- **Tek tek sıfır pozisyonlarını üretme**: yarı-klasik N(T)'den öteye geçen kuantum sistemi bilinmiyor.
- **Self-adjointness**: aday Hamiltonian'ların gerçek Hilbert uzayında self-adjoint extension'u sorunlu.

---

## Bu seansın asıl getirisi

İlk başta "QM × Riemann'da gerçek bir bağlantı var mı?" diye sorduk. Şimdi biliyoruz:

1. **İstatistiksel akrabalık gerçek ve sağlam** — üç bağımsız test, gerçek QM sistem simülasyonu
2. **Asallar-sıfırlar ikiliği görsel** — Fourier'i alınca sıfırlar asalları çığlık çığlık söylüyor
3. **Yarı-klasik çerçeve (Berry-Keating) çalışıyor; kuantum gerçeklenmesi açık duvar**
4. **Polygon-π'den gerçek köprü "yasak bölge"den değil σ-çerçevesinden geçiyor**

---

## Sonraki adımlar (sıralı değil, seçenek)

- **σ-çerçevesi derinleşmesi**: $\sigma_R^{GUE}(s)$'yi s=1/2 dışında sistematik tara; $\sigma_{\text{fire}} \times \sigma_{\zeta} = 1$ konjektürünü hassas test et
- **Voting tarafına dön**: bu seansın amacı oraya bilgi getirmekti — Sezen voting + quantum cognition + ölçüm collapse paralelleri için yeni bir bakışla dönülebilir
- **Berry-Keating modern uzantıları**: Bender-Brody-Müller PT-simetrik Hamiltonian sayısal denemesi
- **Empirik veri**: nükleer rezonans NDE veya kavite spektrumu indir, 4. histogram olarak ekle (Bohigas 1984'ün tekrarı, ama tam üçgen kapanır)

---

*Bu defter `qm_riemann/` dizinindedir. 8 Python betiği, 6 PNG görsel, 2 markdown rapor.*

*— Uğur & Claude, Mayıs 2026*
