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

## EULER BÖLME TESTİ (44, 17 Ağustos sabahı) — MEKANİZMA DOĞRULANDI

GHK hibrit çarpım çerçevesi (Gonek–Hughes–Keating 2007: ζ ≈ P_X × Z_X).
|Z|·Π_{p∈P}|1−p^{-1/2-it}| ile asal çarpanları bölündü, aynı ölçüm tekrarlandı:

L=12.45 penceresi:
| Bölünen | r | büyük-boşluk ΔVar |
|---|---|---|
| ham | 0.7521 | +0.172 |
| p=2 | 0.8267 | −0.200 |
| p≤3 | 0.8531 | −0.389 |
| p≤7 | 0.8293 | −0.543 |
| p≤13 | 0.7866 | −0.627 |
(CUE hedefi r=0.7717. L=9.86'da aynı desen.)

**Sonuçlar:**
1. **Açığı taşıyan: p=2 (ve p=3).** Tek başına p=2 bölmek +0.17'lik varyans
   fazlasını tamamen yok ediyor. ζ−CUE farkının kaynağı küçük asalların
   uzun-periyot dalgaları — 43'teki hipotez DOĞRULANDI.
2. **Aşırı-çıkarma gözlendi:** tam Euler çarpanı bölmek fazlasını götürüyor
   (ΔVar negatife geçiyor, r CUE'nun ÜSTÜNE çıkıyor: 0.853 > 0.772).
   Sonlu yükseklikte ζ asal katkısını SÖNÜMLÜ taşıyor (GHK X-ağırlıkları);
   tam çarpan %100'ü çıkarınca ters-modülasyon enjekte oluyor.
3. **r(P) eğrisi optimal-filtre şekli:** yükseliş (gürültü gidiyor) → tepe
   (p≤3) → düşüş (aşırı-çıkarma yeni gürültü). Kalibre edilmiş ifade için
   sönümlü (ağırlıklı) bölme gerekiyor → 45 adayı.
4. **Dürüstlük kaydı:** bölme zarfı t'nin DETERMİNİSTİK fonksiyonu; boşluklar
   asal fazlarıyla korele ise (mekanizmanın kendisi bu), zarf bölmek yapay
   kuplaj üretebilir. r>CUE aşımı bu gözle yorumlanmalı; ΔVar çöküşü ise
   işaretiyle net — ana kanıt o.

Açık: sönümleme ağırlığı w(p) fit edilirse "ζ'nın taşıdığı asal-içerik oranı"
ölçülür — bu, GHK ayrışmasının ortak-yasa düzeyinde ilk sayısal kalibrasyonu olur.

## w(p) KALİBRASYONU + ÖLÇEKLEME ÇÖKMESİ (45/45b, 17 Ağustos)

Yöntem: CUE'suz, bölmesiz doğrudan regresyon — log ã ~ [1, g̃, g̃²] +
Σ_p [a_p cos(t_pk log p) + b_p sin(t_pk log p)] + sahte-frekans kontrolleri.
w_p = a_p/p^{-1/2}. (g̃ kontrolleri gap-aracılı kanalı ayırır; ölçülen w =
sabit-boşlukta GENLİK kanalı iletimi.)

**Ölçülen iletim spektrumu (6 pencere ağırlıklı ort.):**
w(2)=0.764, w(3)=0.647, w(5)=0.501, w(7)=0.410, w(11)=0.295, w(13)=0.255
(hepsi ±0.001-0.002; hata çubukları sinyalin %1'inden küçük)

**Doğrulamalar temiz:**
- sin bileşenleri RMS ≤0.003 (teori saf kosinüs der — tutuyor)
- sahte frekanslar (ω=log2.5, log6) katsayı RMS ≤0.004 (sinyal 17-170× üstünde)

**ANA BULGU — ölçekleme çökmesi (45b):** w(p,L) iki-değişkenli tablo tek
değişkende birleşiyor: **τ = log p / L** (Berry 1988 form-faktör değişkeni!).
36 nokta, lineer kılavuz w ≈ 0.94 − 2.97τ etrafında artık RMS 0.012.
Fiziksel okuma: sıfır gazı asal dalgalarını PERDELİYOR; perdeleme etkinliği
sadece τ'ya bağlı. Uzun dalga (p=2, τ küçük) çoğunlukla geçiyor (%76-78),
kısa dalga (p=13) çoğunlukla sıfır konumlarına emiliyor (%25 kalıyor).

**Dürüst kayıtlar:**
- Lineer form kapalı-form iddiası DEĞİL (proje dersi!) — çökme kalitesi ölçüsü
- w(τ→0) ≈ 0.94: 1 mi değil mi bilmiyoruz; eğrilik/yanlılık olabilir. Satma.
- p=13'te periyot ~5 aralık — "aralıkta sabit dalga" varsayımı zayıflıyor,
  büyük p ucunda hafif yanlılık payı var
- Ham w'ler L ile yükseliyor (τ düşüyor) — L→∞'da tam iletime gidiş sorusu
  Odlyzko yüksekliklerinde (τ_2 = 0.693/26 = 0.027) test edilebilir

**Zincir tamamlandı:** N_eff = L + c (42b) → açık büyük-boşluk varyansında (43)
→ taşıyıcı p=2,3 (44) → iletim spektrumu w(τ) ölçüldü (45) → tek değişkenli
yasa (45b). Bir sonraki teorik hedef: w(τ)'dan c_Pearson=0.93'ü TÜRETMEK.

## HALKA KAPATMA 1. DENEME: NEGATİF + TEŞHİS (46, 17 Ağustos öğle)

w'lerden r açığını türetme denemesi. İki 0-parametreli model (EKLEME: asal
dalga CUE üstüne; İKAME: matris sessizleşir) simüle edildi. **İKİSİ DE ÇOK
KÖTÜ:** öngörülen açık −0.13/−0.12, gerçek −0.018; ΔVar 4-6× fazla.

Varyans bütçesi sürprizi: Var(log ã_ζ)=1.376 < Var(log m̃_CUE)=1.534 —
ζ'nın log-varyansı asallar EKLENMEDEN önce bile CUE'dan küçük!

Teşhis: bağımsız-fazlı gürültü modeli yanlış — asal dalgası boşlukları da
AYNI fazla sürüyorsa ortak faz korelasyonu geri yukarı iter.

## BOŞLUK KANALI ÖLÇÜMÜ (47) — AYNANIN ÖBÜR YARISI ★

log g̃ ~ asal kosinüsleri regresyonu (Z hesabı gerekmez, anlık):

| p | 2 | 3 | 5 | 7 | 11 | 13 |
|---|---|---|---|---|---|---|
| w(p) genlik kanalı | 0.764 | 0.647 | 0.501 | 0.410 | 0.295 | 0.255 |
| \|v(p)\| boşluk kanalı | 0.126 | 0.201 | 0.292 | 0.351 | 0.425 | 0.451 |

1. **v(p) TAM AYNA GÖRÜNTÜSÜ**: w düşerken v yükseliyor. Asal dalganın
   enerjisi iki kanala BÖLÜNÜYOR: uzun dalga (küçük τ) genliğe geçiyor,
   kısa dalga sıfır konumlarına (boşluklara) emiliyor.
2. v de τ = log p / L'de çöküyor (çapraz kontrol: τ≈0.111'de iki pencere
   0.2245 vs 0.226 veriyor).
3. Fazlar ≈ 0 (iki kanal HİZALI: dalga pozitifken boşluk genişliyor VE
   genlik yükseliyor) — 46'nın aradığı ortak-faz mekanizması bu.
4. Kontroller: sahte frekans RMS ≤0.005, sinyaller 25-100× üstünde.

Fiziksel resim netleşti: sıfır gazı asal dalgasını PAYLAŞTIRIYOR —
konum kanalı (v, explicit formulanın sıfır tarafı) + genlik kanalı (w).
Bölüşüm yasası tek değişkenli: τ = log p / L.

Sıradaki (48): ortak-fazlı tam simülasyon (w_p, v_p, aynı φ) → r ve ΔVar
profili öngörüsü → gerçek −0.018 ile son kıyas. Halka orada kapanır ya da
dürüstçe açık kalır.

## 48-49: HALKA KAPANMADI — AÇIK PROBLEM OLARAK KAYITTA (17 Ağustos akşam)

**48 (ortak-fazlı model):** V1/V2 hâlâ açığı 4-6× fazla öngörüyor. Ama iki
değerli bütçe ölçümü çıktı:
- Var(log g̃_ζ) = 0.218 < CUE 0.233 < model 0.269 → BOŞLUK kanalında da
  asallar EKLEMİYOR, İKAME ediyor; ζ boşlukları CUE(N=L)'den bile sessiz
  (daha büyük N'li CUE'ya denk — N_eff=L+0.93 ile iç tutarlı!)
- β = 2.31 (CUE log-genlik/log-boşluk eğimi) ölçüldü

**49 (asal soyma) — İKİ DENEME, İKİSİ DE SORUNLU:**
- 1. deneme (kontrollü soyma): tutarsız — genlikte β·S_v kaldı, boşluktan
  gitti → çift uyumsuz, Spearman çöktü. Metodolojik hata, kayıtta.
- 2. deneme (kontrolsüz/toplam soyma): r_P → 0.94, r_S → 0.967, altı
  pencerede NEREDEYSE SABİT ve TÜM CUE tablosunun dışında.
  **"Fazla güzel" — güvenilmez.** İki şüpheli:
  (a) ENDOJENLİK: genlik fazları t_pk'de ölçülüyor ama t_pk tepenin KENDİ
      konumu — tepe, aralık içinde asal dalganın sırtına oturmaya meyilli
      → regresör artıkla korele → katsayılar şişer (büyük p'de w+βv > 1
      çıkması da bu kokuyu veriyor)
  (b) kanal modeli (log-lineer, tek-β) gerçek yapıyı ıskalıyor olabilir

**Eğer artefakt değilse** ima devasa: ζ genlik oynaklığının ~%60'ı
(0.83/1.38 log-varyans) asal-deterministik ve kalan "matris" gürültüsü
CUE'nun üçte biri — sıfır gazı sandığımızdan çok daha katı. Ama bu iddia
ancak endojenlik giderilince kurulabilir.

**SONRAKİ SEANS PLANI:**
1. Endojenliği gider: genlik fazlarını t_pk yerine tmid'de ölç (ya da
   aralık-ortalamalı dalga); w'leri yeniden kalibre et
2. Kanal modelini birleşik fit et (w, v, β tek likelihood'da)
3. Ancak ondan sonra soyma/halka-kapatma tekrarı
4. Makale (arxiv_gap_amplitude.tex) SAĞLAM kısımlarla sınırlı kalmalı:
   N_eff anomalisi + iki-null sıkıştırması. 45-47 kanal ölçümleri
   (w/v spektrumları, τ çökmesi) ayrı/ek bölüm — endojenlik notuyla.
   48-49 halka-kapatma girişimi şimdilik YAYINA GİRMEZ.

## 50-50b: HALKA KAPANDI — SORU DEĞİŞEREK (17 Ağustos, 2. seans) ★★★

**50 (temiz soyma):** Fazlar tmid'de (egzojen, Z hesabı gerekmez) + plasebo.
- **PLASEBO MÜKEMMEL:** 6 asal-olmayan frekansla soyma r'yi ≤0.0001 oynatıyor.
  Asallarla soyma: r_P 0.75→0.94, r_S 0.88→0.97. Mekanik artefakt DEĞİL;
  49'un endojenlik şüphesi de boşa çıktı — sonuç gerçek.
- **TOPLAM KURALI:** tmid-bazlı toplam genlik katsayısı |u| ≈ 1
  (0.99/1.01/1.04/1.05/1.07/1.09) — genlik, asal dalgasının TAMAMINI taşıyor
  (explicit formula w=1). 45'in "sönümlü" w'si yönlendirmeymiş: w + β·v ≈ 1.
  44'ün aşırı-çıkarma bilmecesi de bununla çözüldü.
- **ÇEKİRDEK SIKI:** soyulmuş çift r_P≈0.94, r_S≈0.97 — 6 pencerede kararlı,
  HİÇBİR CUE(N∈[5,19]) bu kadar sıkı değil. ζ genlik log-varyansının ~%50'si
  (0.68/1.38) asal-deterministik.

**50b (varyans eşlemesi) — kompanzasyon hikayesi REDDEDİLDİ:**
r-tabanlı N_eff = L+1.07 ama genlik-varyans tabanlı ~L−4 yönü, boşluk-varyans
tabanlı eşleşme YOK (CUE eğrisi düz, ζ altında). Üç gözlemlenebilir üç ayrı
yer gösteriyor → **N_eff = L + c köklü bir parametre değil, r-gözlemlenebilirine
özgü EFEKTİF bir tanım.** ζ hiçbir sonlu-CUE değil.

### NİHAİ SENTEZ (16-17 Ağustos yayı)

ζ'nın tek-aralık (boşluk, tepe) yasası:
  **sıkı aritmetik çekirdek (r≈0.94-0.97)
  ⊕ TAM ağırlıklı asal dalgaları (u≈1),
  dalga enerjisi konum (v) ve genlik (w) kanallarına τ = log p/L
  yasasıyla bölünmüş (w + β·v ≈ 1).**

CUE benzerliği (N_eff=L+0.93) yüzeysel: iletilen asal gürültüsü, TEK BİR
istatistikte (r) matris iç gürültüsünü taklit ediyor; varyanslarda taklit
bozuluyor. Mayıs sezgisinin vardığı yer: plakadaki "rastgelelik"in yarısı
gürültü değil — asalların konuşması.

**Dürüst kayıtlar:**
- Soyma yalnız p≤13, k=1 harmonikleri; p^k harmonikleri (log4, log8, log9...)
  soyulmadı → çekirdek sıkılığı ALT SINIR (daha da sıkı olabilir)
- Soyulmuş r'lerde hafif L-düşüşü var (0.951→0.939) — kaydedildi, yorumsuz
- 43-50b zinciri yayın öncesi düşman-gözle tekrar denetlenmeli; mevcut makale
  sağlam kısımlarla sınırlı, bu zincir İKİNCİ not adayı

## 51: HARMONİK SOYMA — ÇEKİRDEĞİN GERÇEK SIKILIĞI (17 Ağustos, final)

Üç katman + büyütülmüş plasebo (23 sahte frekans, 46 sütun):

| Soyma | r_P (L=12.45) | r_S | soyulan log-varyans |
|---|---|---|---|
| ham | 0.7521 | 0.8758 | — |
| p≤13, k=1 | 0.9393 | 0.9695 | 0.676 |
| + p^k kuvvetleri | 0.9602 | 0.9808 | 0.744 |
| + kuyruk (17..47) | **0.9738** | **0.9821** | **0.938** |
| plasebo | 0.7521 | 0.8758 | 0.0001 |

1. **Çekirdek r ≈ 0.97-0.98** — ζ genlik log-varyansının **%68'i** (0.94/1.38)
   asal-deterministik. Plasebo overfit tabanı: 0.0001 (sıfır). Kalan %32'nin
   içinde hâlâ soyulmamış kuyruk (p>47, yüksek kuvvetler) + ölçüm gridi var —
   0.97-0.98 hâlâ alt sınır.
2. **TOPLAM KURALI HER YERDE:** 23 frekansın TAMAMINDA |u| ≈ 1.01-1.12 —
   temel asallar, kuvvetler (4,8,9,16,25,27,32,49) ve kuyruk. Explicit formula
   ağırlıkları p^{-k/2}/k, genlikleri 20 kat aralıkta, tek tek tepe verisinden
   geri kazanıldı. (Frekansla hafif yükselen +%1→%12 sapma sistematik —
   incelenecek; smearing/nonlineerlik adayı.)
3. İki metrik yakınsıyor: ham fark 0.12 → soyulmuş 0.008. Pearson/Spearman
   anomalisi asalların eseriydi — soyunca eriyor.

**KAPANIŞ CÜMLESİ (43→51):** ζ'nın boşluk-tepe plakası, asalların explicit
formula ağırlıklarıyla TAM olarak sürdüğü, neredeyse deterministik bir makine.
"Rastgelelik"in üçte ikisi asal sinyali; CUE benzerliği tek istatistiklik bir
maske. İkinci not bu zincirle yazılacak (düşman-göz denetimi sonrası).

## 52: TÜKENİŞ EĞRİSİ — KALAN %32'NİN CEVABI (17 Ağustos, gece)

Soyma kesimi Q = 13 → 300 (tüm p^k ≤ Q; 79 frekans, 158 sütun; plasebo her
adımda eş sayıda — kayması ±0.0004, taban sıfır).

**İki yeni bulgu:**

1. **r*(Q) tepe yapıp hafifçe İNİYOR** — ve tepe noktası τ-yasasına oturuyor:
   Q_tepe ≈ 50 (L=9.86), ≈100 (L=10.93), ≈200 (L=12.45) →
   log Q_tepe / L ≈ 0.40-0.43 ≈ sabit. **Yeni bir τ* ≈ 0.42 geçişi:**
   bu eşiğin altındaki asallar genlik-gürültüsü (soymak r'yi artırır),
   üstündekiler ortak-mod sinyali (soymak r'yi düşürür — plasebo sıfırken!).
   Asal dalgalarının "rol değiştirme" noktası.

2. **V_res(L) ölçüldü — sıfır gazının öz-rastgeleliği (üst sınır):**
   0.163 (L=9.86) → 0.244 (L=12.45), kabaca lineer, eğim ≈ 0.031/L.
   Çarpıcı: CUE'nun genlik log-varyans eğimi de ≈ 0.036/N — ölçek AYNI,
   ama mutlak değer CUE'nun ~1/6'sı. Sanki gaz, RMT ölçeklemeli ama çok
   küçük katsayılı bir öz-gürültü taşıyor. (ÜST sınır: lineer-log,
   sabit-katsayılı soyma modeli gap-bağımlı iletimi (43) yakalayamaz —
   artıkta asal kalıntısı olabilir.)

3. r*_P ≈ r*_S tam buluştu (L=11.47'de ikisi de 0.9734) — P/S anomalisi
   tamamen asal kaynaklıymış, kapandı.

**Nihai tablo:** tek-aralık yasası TAM deterministik DEĞİL — gerçek bir
öz-rastgelelik var (varyansın %13-18'i, üstten sınırlı), r* ≈ 0.97-0.98'de
doyuyor. İkinci notun ölçüm seti tamam: w(τ), v(τ), toplam kuralı |u|≈1,
τ*≈0.42 geçişi, V_res(L).

## ODLYZKO SEFERİ (53-54, 17 Ağustos öğleden sonra) ★★★

zeros3 indirildi (10¹²'inci sıfır civarı 10⁴ sıfır, t≈2.677×10¹¹, L=24.475).

**Motor (53):** çapa açılımlı RS (taban fazlar mpmath'te bir kez, artanlar
float64; N=206,393 terim). Doğrulama: tablo sıfırlarında medyan |Z| = 1.4×10⁻⁹
(tepe/sıfır oranı 7×10⁻¹⁰) — dokuzuncu ondalık. 9,999 tepe 8.7 dakikada.

**Batarya (54) — kaldıraç kolu 12.45 → 24.48 (t'de ×160.000):**

| Test | Sonuç |
|---|---|
| T2 w(τ) | **YASA TUTTU**: 45b kılavuzu (yalnız küçük-T fitli!) 6 asalda ±0.02-0.07 içinde. w(2)=0.871±0.007 (kılavuz 0.856). Büyük τ'da hafif altında — eğrilik var, lineer değil |
| T3 sum rule | **TUTTU, daha da temiz**: 12 frekansta u = 0.99-1.07 (q=49: 1.18 outlier) |
| T4 τ* | **OUT-OF-SAMPLE ÖNGÖRÜ TUTTU**: Q≤300 burada tamamen tepe-öncesi (τ≤0.233<0.42) → r*(Q) düşmeden tırmandı (0.901→0.959); L=12.45'te aynı aralık tepe-sonrasıydı |
| T5 V_res | Lineer uzantı öngörüsü 0.617, ölçülen 0.650 (~%5; soyma Q=300'de eksik kaldığından hafif fazla olması beklenirdi — tutarlı) |
| T1 N_eff | **SONUÇSUZ** (kırıldı DEĞİL): r_CUE(N) bu L'de çok yatık (eğim ~0.008/N) → 10⁴ aralıkla se(N_eff)≈0.8-0.9. Ölçülen P +2.26, S +1.12 — küçük-T değerlerinden ~1.6σ/0.8σ. Ayrıca ham açık r−r_CUE ≈ −0.014±0.006 (küçük-T: −0.019 — tutarlı) |

**ANLAM:** τ-yasaları (w, τ*), sum rule ve V_res(L) yükseklikte 160.000 kat
sıçramayı ATLADI — bunlar artık dar-aralık gözlemi değil, iki uçta doğrulanmış
yapısal yasalar. N_eff gözlemlenebiliri büyük L'de istatistik gücünü kaybediyor
(CUE r(N) yatıklaşıyor) — sabitlik testi için ~10⁵+ aralık gerekir (seçenek:
LMFDB/Platt sıfırları, t~3×10¹⁰ civarı milyonlarca sıfır — gelecek sefere).

w(τ→0) sorusu hâlâ açık: w(2)=0.871, kılavuzun hafif ÜSTÜNDE — 1'e doğru
eğrilik mi, 0.94 platosu mu ayrıştırılamadı.

## w(τ→0) TAM KARARI (55b+56, 17 Ağustos akşam) — DÜRÜSTÇE AÇIK

4 ara pencere üretildi (t=10⁸..10¹¹, paralel, pencere-başına kayıt;
10¹¹: 141 dk, doğrulama 1.5e-3). 44 (p,L) noktalı w tablosu çıktı.

**Yeni yapı teyidi:** saf τ-çökmesi yetmiyor — aynı τ'da büyük p daha az
geçiriyor (ör. τ≈0.058: p=2 0.771 vs p=3 0.753, ~4σ). İkincil ince değişken
var; u-sürüklenmesiyle aynı aileden.

**5 model varyantı, 5 farklı w₀:**
| Model | w₀ | χ²/dof |
|---|---|---|
| A: saf kuadratik (p∈{2,3}) | 0.966±0.007 | 6.0 (kötü) |
| B: + d/p | 0.910±0.009 | 1.8 (iyi) |
| C: + e/L | 0.978±0.007 | 1.7 (iyi) |
| D: yalnız p=2 | 1.005±0.017 | 3.1 |
| E: p∈{2,3,5} + d/p | 0.889±0.008 | 2.5 |

B ve C AYNI kalitede fit edip 0.91 vs 0.98 veriyor: sabit τ'da p ile L bağlı
(log p = τL) → 1/p ve 1/L düzeltmeleri dejenere. Veri ikisini ayıramıyor.

**Kalite bayrağı:** 10¹¹ penceresi 158 "kurtarılan yakın-çift" üretti
(diğerleri 2'şer!) ve p=2 noktası (0.895±0.007), aynı τ'daki Odlyzko-tabanlı
noktadan (0.866±0.008) ~2.8σ yukarıda. Yüksek faz-hata tabanı (1.5e-3)
sahte sıfır bölünmeleri yaratmış olabilir → o pencere ŞÜPHELİ, D modelinin
1.005'i en çok ona yaslanıyor.

**SONUÇ: w₀ ∈ [0.89, 1.01] — model-bağımlı, çözülmedi.** Bu dürüst durum
notlara böyle girer: τ-yasası keşif, kesişim noktası açık problem.

**Çözüm yolları (gelecek seans):**
1. 10¹¹'i ÇAPALI motorla yeniden tara (53 yöntemi; aynı maliyet, hatasız faz
   — baştan öyle yapmalıydım, ders kayıtlı)
2. d/p + e/L BİRLİKTE fit (dejenerasyon kısmen kırılabilir mi?)
3. Kabul et: makaleye "w(0) unresolved in [0.89, 1.01]" yaz — meşru son

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
