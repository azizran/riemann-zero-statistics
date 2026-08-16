# Literatür Taraması — 16 Ağustos 2026

*17 Mayıs notundaki 1. açık adımın kapatılması: r=0.83 gap-amplitude korelasyonu ve
+6.67σ Gaussian-PSD testi literatürde var mı?*

---

## ANA SONUÇ

**Sayısal bulgularımızın (Pearson r, surrogate null testi) doğrudan karşılığı bulunamadı.
Ama nitel bağlantı sandığımızdan daha derin biçimde biliniyor — ve tam bu obje üzerinde
2024-25'te aktif çalışma var.**

---

## 1. KRİTİK KEŞİF: Conrey–Ghosh 1985 (tam metin okundu)

**"A mean value theorem for the Riemann zeta-function at its relative extrema on the
critical line"**, J. London Math. Soc. (2) 32 (1985) 193–202.

RH varsayımı altında kanıtlanmış teorem:

$$\frac{1}{N(T)} \sum_{0<\gamma\le T} \max_{\gamma < t \le \gamma^+} |\zeta(\tfrac12+it)|^2 \sim \tfrac12(e^2-5)\log T \approx 1.1945\,\log T$$

Yani: **her sıfır aralığındaki max|Z|²'nin ortalaması** kanıtlı olarak biliniyor
(ortalama |ζ|² olan log T'den ½(e²−7) ≈ 0.1945 kat büyük).

Ayrıca aynı makalede:
- Kaydırılmış moment: Σ|ζ(½+it_γ+iα·2π/L)|² ~ C(α)·(T/2π)L², C(α) açık formüllü
- Ardışık maksimumlar arası boşluk bazen > 1.4·(2π/L)

**Bizim için anlamı:** Literatür max|Z|'nin *koşulsuz ortalamasını* biliyor.
Bizim ölçtüğümüz şey — max|Z|'nin *aralığa koşullu* davranışı (regresyon eğimi,
Pearson r) — Conrey–Ghosh'un **gap-koşullu rafinesidir**. Makalede joint (gap, max)
istatistiği YOK, korelasyon katsayısı YOK.

**Bonus — pipeline sanity check:** Verimizde max² ortalamasının 1.1945·log(t/2π)'ye
oturması gerekir. Kanıtlı teoremle kendi verimizi test edebiliriz. (Yapılacak!)

## 2. Aynı damarın 2024–25 devamı (taze literatür!)

- **Hughes–Lugmayer–Pearce-Crump** (arXiv:2411.05573, JLMS 2025): aynı ikinci momentin
  alt-mertebe terimleri. Yine sadece ortalama — gap-koşullu analiz yok (tam metin tarandı).
- **Pearce-Crump** (arXiv:2411.05568, Mathematika 2025): ekstremumlarda birinci moment
  ve türevler. Yine joint istatistik yok.

**Anlam:** Camia şu anda tam bu objeyle (sıfırlar arası ekstremumlar) ilgileniyor.
Niş güncel ve canlı — hem fırsat hem risk (joint istatistiği yakında biri yapabilir).

## 3. Nitel bağlantının bilinen kökleri (dürüstlük notu)

r=0.83'ün **işareti** uzmanlar için sürpriz değil:

- **Lehmer çiftleri** (1956): küçük aralık ↔ küçük |Z| tepesi — folklor düzeyinde bilinir.
- **Gonek–Hejhal konjektürü**: |ζ'(ρ)| dağılımı ↔ komşu sıfır mesafesi ilişkili.
  Bizim (gap, max|Z|) ölçümümüz bunun integre versiyonu.
- **Conrey–Ghosh–Gonek 1986** büyük-boşluk yöntemi: "büyük boşluk ↔ yakınında büyük
  değer" prensibi yöntemin içine gömülü (Mueller 1981 hattı).

**Yeni olan işaret değil, SAYI**: Pearson r'nin kendisi, +6.67σ surrogate testi,
ve gap-koşullu eğim ölçümü.

## 4. Gaussian-PSD (surrogate) testi

Faz-randomizasyon / surrogate data yöntemi (Theiler et al. 1992) nonlineer zaman serisi
literatüründe standart, ama **Z(t)'ye uygulanmış hali bulunamadı**. "Aynı güç spektrumlu
Gaussian süreçten +6.67σ sapma" ölçümü görünüşe göre raporlanmamış. Bu, en savunulabilir
özgün katkımız.

## 5. Yakın ama farklı işler (karışmasın)

- arXiv:2511.18275 (Kas 2025): Hardy Z varyasyonel yaklaşımla Montgomery pair correlation —
  sıfır ARALIKLARı hakkında, genlik yok. (İddialı preprint, hakemli mi belirsiz.)
- arXiv:2507.10193 (2025): CUE spacing ratios finite-size düzeltmeleri vs ζ sıfırları —
  yine sadece spacing.
- FHK / Najnudel / Arguin et al.: uzunluğu O(1) aralıklarda max (çok sıfır içerir) —
  bizimki tek-gap ölçeği, farklı obje (17 Mayıs notundaki değerlendirme doğru).

## 6. Kalan risk

- Ivić, *The Theory of Hardy's Z-Function* (Cambridge 2012) — içindekiler çekilemedi
  (paywall/500). Ama risk artık düşük: gap-maksimum damarının otoritatif hattı
  CG85 → HLP-C 2024-25, ve bu taze makaleler joint istatistik içermiyor; literatürde
  olsaydı onların girişinde atıf olurdu.
- Hall (1999) Z'² momenti: yapısal köprü duruyor, Pearson r içermiyor (2. ajanın
  Mayıs tespiti bu taramayla tutarlı).

---

## SONUÇ VE KONUMLANDIRMA

Makale çerçevesi (dürüst versiyon):

> *"Conrey–Ghosh (1985) sıfır aralıklarındaki |Z| maksimumlarının ortalamasını kanıtladı;
> Hughes–Lugmayer–Pearce-Crump (2025) alt-mertebeleri hesapladı. Biz bu maksimumların
> aralık uzunluğuna KOŞULLU dağılımını sayısal olarak inceliyoruz: Pearson r = 0.83
> (normalize 0.92), aynı güç spektrumlu Gaussian surrogate'lerden +6.67σ sapma —
> korelasyonun ζ-özgün, spektrum-ötesi olduğunun göstergesi."*

### Sonraki adımlar

1. ~~**CG85 sanity check**~~ → **YAPILDI (aynı gün, `37_conrey_ghosh_check.py`):**
   - Eğim fiti (40 pencere, T=100k verisi): mean(max²) = **1.1926**·log(t/2π) + 2.776
   - Kanıtlı CG sabiti 1.19453'e sapma: **−0.16%**; optimizer-bias düzeltmesiyle **−0.05%**
   - Bias sondası (250 aralık, ince grid): 36'nın kaba optimizeri max²'yi sadece
     ‰1.1 eksik tahmin ediyor — ihmal edilebilir
   - Ders: Cesàro oranı T=75k'da hâlâ 1.38/1.15 (yavaş yakınsama, O(1) terimi) ama
     eğim sabiti hemen yakalıyor — önceki slope² fiyaskolarının tersine burada
     teorem gerçekten a·log t + b formunda
   - **PIPELINE DOĞRULANDI**: r=0.83 ölçümü kanıtlı teoremle tutarlı zeminde
   - ~~Açık uç: b karşılaştırması~~ → **YAPILDI (`38_hlpc_alt_mertebe.py`):**
     HLP-C Teorem 1'den yerel ortalama türetildi:
     mean(max²|t) = A·L + (α₋₁+2A) + (α₀+α₋₁)/L,  L = log(t/2π)
     - α₋₁ = 5−e²−10γ₀+2e²γ₀ = 0.368945, α₀ = −0.423285 (Stieltjes γ₁ ile)
     - Teorik sabit terim: 2.7494 | Ölçülen (37): 2.7760 → sapma **+0.97%**
     - Sıfır serbest parametreli teori eğrisi: 26 pencerenin **26'sı** 2σ içinde,
       ortalama |z| = 0.06 — veri teoremin üstüne yapışık
     - Dürüst not: HLP-C'nin kendi makalesinde de numerik bölüm var (Section 6),
       "ilk sayısal doğrulama" iddiası YOK; bizimki bağımsız yöntem ve aralıkta
       (per-gap optimizer, T≤75k) bir çapraz doğrulama
     - **ÇİFTE DOĞRULAMA TAMAM**: eğim → CG 1985 (−0.05%), sabit terim → HLP-C
       2024 (+0.97%). Pipeline artık iki kanıtlı/yayınlanmış sonuçla çapalı.
2. ~~r'nin beklenen kısmının ayrıştırılması~~ → **YAPILDI (`39_cue_null.py` + `40_unfold_recheck.py`):**

   **CUE null** (Keating–Snaith modeli, N ↔ L yoğunluk eşlemesi):
   - r_CUE(N) ilk kez ölçüldü (N=5..16, 100k aralık/N): 0.9295 (N=5) → 0.7222 (N=16)
   - Bu tablo literatürde yok gibi — kendi başına küçük bir katkı
   - ζ tarafı unfold edildi (g̃ = gap·L/2π, ã = amp/√(AL+b), 38'in doğrulanmış ölçeğiyle)
   - Pencere-içi kayma düzeltmesi farkın sadece ~0.003'ünü götürdü

   **ANA BULGU:** r_ζ − r_CUE(N=L) = **−0.019**, altı pencerede DÜZ (−0.0188…−0.0200),
   anlamlılık pencere başına 3–13σ. Eşdeğer ifade: **N_eff = L + 0.84 ± 0.08**
   (ζ, kendi boyutundan ~0.84 büyük CUE gibi korele).

   **Üç null'un hikayesi (makalenin omurgası olabilir):**
   | Null | r | ζ'ya göre |
   |---|---|---|
   | Gaussian-PSD (spektrum var, sıfır yapısı yok) | 0.505 | çok DÜŞÜK (+6.67σ) |
   | CUE N=L (tam RMT yapısı) | 0.841 | hafif YÜKSEK (−0.019) |
   | **ζ gerçek** | **0.822** | — |

   ζ iki null'un ARASINDA: korelasyon Gaussian'dan çok güçlü, saf RMT'den ölçülebilir
   zayıf. Aradaki −0.019, asal-salınım (explicit formula) dalgalarının gap'ten bağımsız
   genlik varyansı eklemesiyle tutarlı — dekorele edici aritmetik katkı adayı.

   **Literatür bağlamı:** BBLM 2006 (arXiv math/0602270) "sonlu-yükseklik ζ ↔ sonlu-N
   CUE, aritmetik sabitli N_eff" paradigmasını SPACING DAĞILIMI için kurdu
   (N_eff = L/√(12Λ), Λ=1.57314 — farklı gözlemlenebilir, farklı eşleme).
   Bizim gap-genlik r'si bu paradigmaya YENİ bir gözlemlenebilir ekliyor; +0.84
   kaymasının aritmetik türetimi açık problem.

   **Dürüst kayıtlar:**
   - N_eff−L'de hafif trend (+0.06/L) var ama 6 noktayla gürültü içinde → 41-42b'de çözüldü
   - ζ genlikleri kaba optimizerden; 37 sondası bias'ın ihmal edilebilir olduğunu gösterdi
   - N↔L eşlemesindeki O(1) konvansiyon belirsizliği tek başına −0.019'u açıklayabilirdi;
     bunu N_eff−L'nin L-boyunca sabitliğiyle karşıladık ama aralığımız dar (L: 5.6–9.1)

## PLAN A SONUCU: BÜYÜK T (41 + 42 + 42b, aynı gün akşam)

**Altyapı (`41_bigT_scan.py`):** Vektörize Riemann-Siegel (C0 düzeltmeli) motoru.
- mpmath'e karşı max hata 7.6e-6; 36'nın çifte-doğrulanmış verisiyle örtüşme
  bölgesinde aralıklar 1e-4, genlikler ‰0.8 içinde eşleşti (iki FARKLI motor!)
- 6 yeni pencere, t=1.6M'a kadar (L=12.45), 240k yeni aralık, 6.2 dakika
- Yakın-çift kurtarma: şüpheli aralıklarda (iç min|Z|<0.1) 10× ince tarama

**42'deki tuzak:** ham CUE tablosuyla interpolasyon zayıf trend gösterdi (Δχ²=7.4).
42b pürüzsüz fit (r_CUE(N) = r∞ + c₁/N + c₂/N², RMS 0.0037/0.0015) bunu ARTEFAKT
olarak teşhis etti — tablo zikzağı N_eff'e ±0.1-0.2 sahte dalga katıyormuş.

**NİHAİ TABLO (12 pencere, L = 5.6 → 12.45, t = 200 → 1.6M):**

| Metrik | N_eff − L | χ²/dof (sabit) | Trend Δχ² |
|---|---|---|---|
| Pearson | **+0.93** | 0.09 | 0.5 (yok) |
| Spearman | **+1.86** | 0.06 | 0.1 (yok) |

1. **Kayma iki metrikte de L'den BAĞIMSIZ** — t iki büyüklük mertebesi değişirken
   sabit. "ζ, N = L + c boyutlu CUE gibi korele" gözlemi güçlü zeminde.
2. **AMA c metriğe bağlı: Pearson +0.93 ≠ Spearman +1.86.** Tek bir N_eff iki
   projeksiyonu birden açıklayamıyor → ζ'nın (gap, max) ortak dağılımı hiçbir
   sonlu-CUE ile TAM eşleşmiyor; fark mutlak-CUE-ötesi (aritmetik) yapı imzası.
   Bu, −0.019'dan daha keskin bir ifade.
3. BBLM paralelliği tamamlandı: onların spacing-şekil N_eff'i, bizim korelasyon
   N_eff'lerimiz — her gözlemlenebilir kendi efektif boyutunu seçiyor.

## AÇIĞIN ANATOMİSİ (43, aynı gün gece)

−0.019 nerede yaşıyor? 3 pencere (L=7.0/9.9/12.4) vs CUE(N=L) koşullu eğriler:

1. **Koşullu ORTALAMA E[ã|g̃] eşleşiyor** (fark ~0.000-0.007, bin RMS ~0.03).
   Yani boşluk-tepe REGRESYON EĞRİSİ RMT-evrensel — bağın şekli aynı.
2. **Açık tamamen koşullu VARYANSTA ve büyük boşluklara yığılmış:**
   g̃ < 1'de fark ≈ 0; g̃ = 1.3'ten sonra hızla büyüyor (L=12.45'te +0.05 →
   +0.32). ζ'nın GENİŞ boşluklardaki tepe yüksekliği CUE'dan daha oynak.
3. Beklediğim "boşluktan bağımsız eklenen varyans" (additive) modeli YANLIŞ
   çıktı; çarpımsal modelden de dik. Büyük boşluk ↔ asal-rezonans bölgeleri
   birlikteliği adayı — test edilmedi henüz.
4. UYARI: en uç bin (g̃∈[2.2,3.5]) işaret değiştiriyor (−0.9…−1.1) —
   az örnek / kuyruk etkisi olabilir, üstüne iddia kurulmamalı.

Sonraki keskin test (44 adayı): |Z|'yi küçük asalların kısmi Euler çarpanına
bölüp (p=2,3,5,7; uzun periyotlu modülasyonu taşıyanlar bunlar) büyük-boşluk
varyans fazlasının çöküp çökmediğine bakmak. Çökerse mekanizma = asal dalgaları,
kanıtlanmış olur.

**Dürüst kayıtlar (güncel):**
- r∞ ekstrapolasyonları (Pearson 0.52, Spearman 0.77) N≤19 fitinden — güvenilmez,
  satılmamalı; sadece interpolasyon aracı
- Pürüzsüz fit RMS'i MC hatasının ~3 katı → c değerlerinde ±0.15 sistematik pay
- χ²/dof'ların çok küçük olması (0.06-0.09) hata çubuklarının muhafazakâr
  olduğunu gösteriyor; DÜZLÜK sonucu bundan etkilenmez
- 36 (mpmath) ve 41 (RS) motorlarının pencereleri L=9.1/9.9 sınırında pürüzsüz
  birleşiyor — motorlar-arası tutarlılık görünür durumda

3. ~~Kısa not taslağı~~ → **YAZILDI: `arxiv_gap_amplitude.tex`** (aynı gün gece).
   Omurga: (i) veri + CG85/HLP-C çapa, (ii) iki null (Gaussian-PSD, CUE) +
   r_CUE(N) tablosu, (iii) N_eff = L + c gözlemi (Pearson +0.93, Spearman +1.86,
   iki Observation olarak), (iv) 4 açık problem. TeX makinede derlenemedi
   (pdflatex yok) — arXiv/Overleaf'te derlenecek. Gözden geçirilecekler:
   e-posta adresi, Milinovich/Lehmer künye teyidi, +6.67σ aralık ifadesi.

*— 16 Ağustos 2026, tarama: Claude (WebSearch + tam metin: CG85 PDF, arXiv 2411.05573/68)*
