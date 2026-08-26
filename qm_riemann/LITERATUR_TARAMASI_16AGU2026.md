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

## 57: ÇAPALI 10¹¹ + NİHAİ w₀ KARARI (17 Ağustos, gece) ★

Çapalı motor 10¹¹'de: doğrulama 3.1e-5 (f64'ün 50'de biri), 20.000 aralık.
**İki kayıtlı öngörü de TUTTU:**
1. Buharlaşma: 158 "kurtarılan yakın çift" → **2** (156'sı faz-gürültüsü
   sahtesiymiş — bayrak haklıydı)
2. Yakınsama: p=2@L=23.49 noktası 0.895 → **0.863±0.005**, Odlyzko-tabanlı
   komşusuyla (0.866±0.008) artık 0.4σ uyumda; p=3 çifti de 0.1σ.
   İki bağımsız motor (bizim çapalı tarama / Odlyzko tablosu) bitişik
   yüksekliklerde aynı fiziği veriyor.

**NİHAİ FIT (temiz veri):**
| Model | w₀ | χ²/dof | 1'e uzaklık |
|---|---|---|---|
| A: saf kuadratik | 0.955±0.007 | 4.6 | 6.6σ |
| B: +d/p | 0.903±0.009 | 0.59 | 10.7σ |
| C: +e/L | 0.967±0.007 | 0.34 | 4.7σ |
| D: yalnız p=2 | 0.976±0.016 | 0.45 | 1.5σ |
| E: geniş +d/p | 0.884±0.008 | 1.4 | 15.3σ |

Temiz pencereyle B ve C artık MÜKEMMEL fit (0.59/0.34) ve D'nin 1.005'i
0.976'ya indi (sahte pencereye yaslanıyormuş).

## DÜŞMAN-GÖZ DENETİM RAPORU ÖZETİ (17 Ağustos, gece — bağımsız ajan)

~60 sayı ham veriden yeniden üretildi; mekanik şüpheler (soyma şişirmesi,
grid asimetrisi, unfolding döngüselliği) ÇÜRÜDÜ. Ama 4 KIRMIZI bulgu:

**K1 (Not 1, ölümcül):** "N_eff − L sabit" büyük olasılıkla interpolant
artefaktı — 42b'nin kuadratik CUE fiti kendi MC hatalarına karşı χ²/dof≈29
ile reddediliyor; yeterli fitlerle kayma L ile YÜKSELİYOR (+0.74→+1.24).
(42'deki trendi 42b'de "zikzak" diye silmiştik — yanlış teşhismiş; asıl
sorun kuadratik ailenin sistematik misfit'iymiş.) KURTULAN: iki-metrik
yarılması (~0.9, interpolant-bağımsız) ve iki-null sıkıştırması.
**K2:** "plasebo ≤1e-4" cümleleri abartılı (52/54'te 1e-3 mertebesi) +
plasebo çekim aralığı kaçınma listesini aşıyor (gerçek asallar karışıyor).
**K3:** "eighteen orders" yanlış — v ölçümleri 16.1 mertebe (36 pencerelerinde
v hiç ölçülmedi; ölçülürse ~18 gerçek olur — anlık iş).
**K4:** sum rule seviyesi ve u-drift'i faz-çerçevesine ±%5-15 duyarlı;
pürüzsüz çerçevede u≈0.96, drift kayboluyor (w₀<1 ile tutarlı!).
"tmid egzojen" ifadesi daraltılmalı (genlik-ölçümüne-göre egzojen).

SARI'lardan önemliler: w+βv≈1 sayısal tutmuyor (%30'a dek); kombinasyon
frekanslarında (log p±log q) 11-13σ faz-kilitli içerik (YENİ BULGU —
u-drift/V_res yorumlarını etkiler); Gaussian null künyesi yanlış
(t∈[50,500]); τ*=0.40±0.03 (0.42 değil); "incomplete transparency"
yumuşatılmalı (D modeli 1'le 1.5σ uyumlu); Not1/Not2 çapraz pürüzler.

HÜKÜM: bu haliyle arXiv'e gitmez; veri altyapısı sağlam, iddia katmanı
revizyon ister. Not 1 Observation-2 omurgasına yeniden kurulmalı (CUE MC
büyütülerek); Not 2'de K2-K4 + S düzeltmeleri.

## 59-60: N_eff NİHAİ KARARI (K1 çözüldü) (18 Ağustos gecesi) ★★

Hassas CUE kampanyası (16 N × 1.5M aralık, ampirik batch-hataları) + PCHIP
(form varsayımsız). NİHAİ TABLO:

| Metrik | Sabit model | Lineer eğim | Anlamlılık |
|---|---|---|---|
| Pearson | +1.04 (χ²/dof 1.89, kötü) | **+0.098±0.022 /L** | **4.5σ — YÜKSELIYOR** |
| Spearman | **+1.99 (χ²/dof 0.34, iyi)** | +0.046±0.027 | 1.7σ — sabitle uyumlu |

1. **Denetçi haklıydı:** Pearson kayması sabit değil (+0.75 → +1.30).
   42b'nin "sabit +0.93"ü kuadratik-fit artefaktıydı. Ders: fit ailesi
   yerine hassas noktalar + monoton interpolasyon.
2. **Spearman kayması ≈ +2.0 SABİT** — sabitlik hikâyesi rank-metriğinde
   yaşıyor.
3. **İki-metrik yarılması hayatta: ort +0.93 ± 0.10** (interpolant-bağımsız).
4. **Yeni resim Not 2 ile daha tutarlı:** P/S anomalisi asal-kaynaklı
   (soymayla eriyor); Pearson'ın L-sürüklenmesi de doğal olarak iletilen
   asal içeriğin L-bağımlılığı — rank-tabanlı N_eff sabitken genlik-duyarlı
   N_eff asal dalgalarıyla sürükleniyor. Not 1'in yeni omurgası:
   "rank-based effective dimension is constant (+2.0); the Pearson one
   drifts in the manner expected from transmitted prime content."

Not 1 yeniden yazımı bekliyor (sonraki iş).

## BÜYÜK REVİZYON UYGULANDI (18 Ağustos gecesi) — commit 5d912d8

**Not 1 yeni omurga:** abstract (iii)/(iv) yeniden; PCHIP metodolojisi
(kuadratik-fit hikâyesi dürüstçe anlatıldı — "an earlier draft of this very
note fell the same way"); Observation 1 = Spearman sabit +1.99 & Pearson
drift 4.5σ; Observation 2 = yarılma 0.93±0.10; tablo/figür/CUE tablosu
(1.5M kampanya değerleri) yenilendi; Gaussian null künyesi düzeltildi
(t∈[50,500]); companion2 atıfı eklendi.

**Not 2 denetim düzeltmeleri:** giriş yeni resmi anlatıyor; "sixteen orders"
(K3); plasebo cümleleri dürüstleştirildi (K2 metin kısmı); τ*=0.40±0.03
(S4); w+βv %5-30 nitelemesi (S1); transparency yumuşatıldı (S5); Odlyzko
u hataları (S6); scripts 27-60; factor 13.

**61 (çerçeve kontrolü) — K4 DOĞRULANDI VE ÇÖZÜLDÜ:** pürüzsüz (RvM-akışı)
fazlarla dokuz asalın hepsi u = 0.95-0.97, drift YOK (tmid: 1.01→1.12
driftli). u-drift bilmecesi = çerçeve artefaktı. Pürüzsüz seviye 0.96 ≈
w₀ bulgusuyla tutarlı — nota işlendi.

**Denetimden kalan işler:** K2 tasarım düzeltmesi + yeniden koşum (52/54
plasebo çekim-aralığı); S2 kombinasyon-çizgileri bulgusunun nota eklenmesi;
S3 ideali (Gaussian null'u ana pencerelerde tekrar); S11 novelty literatür
turu (gönderim öncesi); K3 opsiyonu (36'da v ölç → gerçek ~18 mertebe).

## 58: v-YASASI 18 BÜYÜKLÜK MERTEBESİNDE (17 Ağustos, gece geç) ★★★
(NOT: denetim K3 — doğrusu 16.1 mertebe; 36 pencerelerinde v ölçülürse ~18)

zeros4/zeros5 indirildi (10²¹ ve 10²²'inci sıfırlar; L=44.58 ve 46.83).
Bu yüksekliklerde genlik hesabı imkânsız ama v-kanalı yalnız sıfır
konumlarından ölçülür — tablolar yetti (Z hesabı YOK, anında).

**SONUÇ: 12/12 nokta kılavuza oturdu.** Küçük-τ kılavuzu v ≈ 2.014·τ
(yalnız t ≤ 2.7×10¹¹ verisine fit) yeni noktaları öngördü:
- 10²¹: v(2)=0.035±0.010 (kılavuz 0.031), ... hepsi ≤0.5σ içinde
- 10²²: v(2)=0.0297±0.0095 (kılavuz 0.0298 — 0.01σ!)
- Sinyaller 3-5σ (12/12 tespit); sahte frekanslar gürültü tabanında

**Anlam:** boşluk-kanalı yasası artık t ~ 10³'ten 1.4×10²¹'e — 18 büyüklük
mertebesi — doğrulanmış durumda. 10²²'inci sıfırın komşuları bile asalların
konumunu tam öngörülen genlikte fısıldıyor. (Explicit formulanın nitel
içeriği bilinir; ölçülen İLETİM KATSAYISININ tek-değişkenli τ-yasasıyla
bu menzilde doğrulanması bizim katkı.)

Not 2'ye eklenecek; v ≈ 2τ küçük-τ davranışı teorik türetim için ekstra
ipucu (lineer başlangıç → perdeleme teorisinin ilk katsayısı).

**KARAR: w₀ < 1 BEŞ MODELDE BİRDEN — perde τ→0'da tam saydamlaşMIyor.**
Kalıcı tutma en az ~%2 (D), muhtemelen %3-10 (B-C aralığı). Kesin değer
1/p-vs-1/L dejenerasyonu yüzünden [0.88, 0.98] içinde açık; ama nitel soru
kapandı: sonsuz uzun dalga bile sıfır gazından kayıpsız geçemiyor.
Notlara işlenecek ifade: "the transmission does not reach unity as τ→0;
a persistent absorption of 2-12% survives, its exact value blocked by a
1/p-vs-1/L degeneracy."

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

## KALAN DENETİM İŞLERİ TAMAMLANDI (18 Ağustos) — 62-65

| İş | Sonuç |
|---|---|
| K3 (62) | v-kanalı 36 pencerelerinde ölçüldü → menzil GERÇEK 18.1 mertebe (t=1148'den); BONUS: v doyumu keşfi (~0.55-0.60, τ≈0.4 civarı — τ* ile aynı bölge) |
| S2 (63) | Kombinasyon çizgileri sistematik: çarpım (log pq) çizgileri ~0.005, 13/13 cos-negatif kilitli (P~2⁻¹³); fark çizgileri YOK → karesel-doğrultma çürüdü, ANALİTİK-yanıt ipucu (e^{i(θp+θq)}, eşleniksiz); V_res etkisi ~1e-4 (önemsiz). Not: denetçinin 11-13σ/0.04-0.05 ön-sondajı bizim sistematikle uyuşmadı (muhtemelen farklı konvansiyon); nota BİZİM ölçüm girdi. Not 2'ye bölüm + açık problem 7 |
| K2 (64) | Plasebo tasarımı düzeltildi (tam kaçınma listesi ≤360); tabanlar: 9-frek 6e-4, 79-frek 5e-3 (r-kayması). Sonuç: L=12.45 tepe-sonrası düşüş (−0.0013) taban içinde → τ* kanıtı alçak pencerelerde (−0.013, tabanın katbekat üstü). 52/54 scriptleri + Not 2 güncellendi |
| S3 (65) | Modern Gaussian null: L=9.86, 40k aralık, 40 surrogate → r_null=0.494±0.007 vs ζ 0.794 = **+43σ** (eski: +6.7σ). Grid-tabanlı ζ r'si optimizer-tabanlıyla 3 ondalıkta aynı (çapraz doğrulama). Not 1 abstract+metin güncellendi |

Gönderim öncesi kalanlar: S11 novelty literatür turu; kullanıcının son okuması;
arXiv hesap/endorsement mekaniği.

## S11 NOVELTY TURU (18 Ağustos) — TAMAMLANDI

1. **Lehmer literatürü**: nitel kalıyor (Odlyzko'nun min|Z| sayımları var,
   Pearson r yok) → Not 1 iddiası sağlam.
2. **YENİ ATIF — FGK 2018** (Fyodorov–Gnutzmann–Keating, J. Phys. A 51,
   464001): CUE'da GLOBAL maksimum ↔ büyük boşluk korelasyonunun ilk nicel
   kanıtı. Bizim aralık-başına tablo değil ama en yakın akraba → Not 1'e
   atıf + "closest published relative" nitelemesi eklendi.
3. **YENİ ATIF — Ford–Zaharescu 2005** (J. reine angew. Math. 579) +
   FSZ 2009 (Math. Ann. 343): {γ·log p/2π} kesirli kısımlarının asal-bağlı
   non-uniform limit dağılımı — v-KANALININ klasik atası! Bizim katkı
   "aralık-başına iletim katsayısı + τ-yasası" olarak konumlandı, Not 2'ye
   atıflar eklendi.
4. w/v kanal-ayrışması + τ-iletim ölçümü türünde başka iş bulunamadı.

İki not da derlendi. Gönderim öncesi kalan: kullanıcının son okuması +
arXiv mekaniği.

## SIÇRAMA 1: KRAMERS-KRONIG ARAYIŞI — 1. TUR (66/66b, 18 Ağustos)

**Soru:** w(τ) ve v(τ) tek kompleks nedensel yanıtın iki yüzü mü?

**Yol açan iki yapısal bulgu (kalıcı değerli):**
1. Parite uyumu: w(0) sonlu (çift), v(0)=0 lineer başlangıç (tek) — nedensel
   χ'nin Re/Im parite yapısıyla uyumlu.
2. FAZ MUHASEBESİ: boşluk = sayımın türevi → türev cos'u sin'e çevirir →
   w EŞ-FAZLI, yer-değiştirme U := v/τ ÇEYREK-FAZLI. Doğal çift (w, U).
   Ve U(0)/w(0) = 2.01/0.94 ≈ 2.14 ≈ β (bağımsız ölçülen log-log eğim!) —
   sum rule'un β'sı faz muhasebesinden kendiliğinden düşüyor.

**Testler ve sonuçlar:**
- 66 (Lorentz osilatörleri, w-fit → v öngörüsü): RED. Fit γ'yı sınıra itti
  (aşırı-sönümlü rejim istiyor), öngörülen v-eğimi 4.5 (ölçülen 2.01).
- 66b (Debye sürekliliği, pozitiflik=nedensellik, ortak ρ): RED.
  Ama ASİMETRİK red: her kanal AYRI AYRI pozitif-Debye ile güzel fit
  oluyor (w kendi Re-formuyla RMS 0.016 ≈ per-asal taban; U kendi
  Im-formuyla χ²/dof ≈ 1.1) — fakat İKİ KANAL FARKLI SPEKTRUM İSTİYOR
  (ortak ρ: Δχ²=1034, μ sınıra kaçıyor).

**1. tur hükmü:** perde TEK skaler dielektrik değil. Eğer KK burada
yaşıyorsa ya TENSÖREL (çift-kırılımlı malzeme gibi: iki kanal iki ayrı
duyarlılık) ya da eşlenik çift başka. w+βv≈1'in p ile büyüyen sapması da
aynı semptom (τ-bağımlı kuplaj).

**2. tur adayları:** (i) μ(τ) serbest bırakılıp gereken kuplajın şekli
çıkarılır — sum-rule sapmasıyla aynı fonksiyon mu? (ii) numerik Hilbert
dönüşümü (çıkarmalı) doğrudan U ölçümünden w'yi kestirmeyi dener;
(iii) paralelde Sıçrama 3 (Spearman +2.0 türetimi) hızlı zafer olarak.

## SIÇRAMA 1, 2. TUR (67 + hızlı testler): ADAY YASA BULUNDU ★★★

Zincir: iki-akışkan modeli de RED (χ²/dof 10.8) → AMA Test B mücevheri:
sum-rule ihlali g(τ) = w + 2.14v − 1 küçük-τ'da SIFIR (−0.003), sonra
büyüyor; iki-akışkan artığıyla corr −0.85 (iki anomali = tek fonksiyon ✓).

g'nin şekli: g ≈ 1.167·v² (corr 0.96, artık tabanda). Ve 1.167 ≈ β²/4 =
1.145 → cebir katlanıyor: w + βv − 1 = (βv/2)² ⟺ **w = (1 − βv/2)²**

**ADAY YASA (genlik korunumu):  √w + (β/2)·v = 1,  β/2 ≈ 1.07**
- Serbest fit: √w = 0.9913 − 1.0313·v, artık RMS 0.0091 (per-asal taban!)
- a≡1 sabitli: b = 1.0683 ≈ β/2 = 1.07 (faz muhasebesinden bağımsız türetilen β!)
- corr(√w, v) = −0.9955 (44 eşleşmiş çift, 11 pencere × 4 asal)
- Aynı parametre sayılı düz-lineer w(v) fitinden 1.6× iyi
- Fiziksel okuma: ışın-bölücü ünitaritesi gibi — GENLİK-iletimi (√w) ile
  konum-soğurması (v) lineer bölüşüyor; şiddet değil genlik korunuyor.

**DİSİPLİN NOTU (Mayıs dersleri):** bu kapalı-form-komşusu bir iddia —
ama numeroloji değil: 2-değişkenli fonksiyonel ilişki, 44 noktada taban
RMS'de, katsayısı bağımsız yoldan (faz muhasebesi β'sı) türetilmiş.
ÖLDÜRME TESTLERİ (3. tur): (1) p=11,13 için w ölç (yasada HİÇ kullanılmadı)
→ out-of-sample; (2) Odlyzko penceresinde p=11,13 çiftleri; (3) düşman-göz.
Testlerden geçerse Not 2'ye "conservation law" bölümü + belki Not 3.

## ÖLDÜRME TESTİ 1 (68): YASA SAĞ ÇIKTI (18 Ağustos, gece)

p=11,13 için w İLK KEZ ölçüldü (22 yeni çift; yasa türetiminde yok;
τ=0.10-0.26 — eğitim menzilinin ötesi; √w 0.45-0.78 aralığı):

- Out-of-sample RMS = 0.0148 (eğitim tabanı 0.0091) — yasa çizgisi izleniyor
- 22 çiftin 16'sı ≤1.5σ; en büyük sapmalar p=13'ün büyük-L noktaları
  (−0.037@3σ, −0.034@4σ) — bilinen sabit-τ per-asal yarılma deseniyle uyumlu
- Hafif sistematik: ortalama −0.009 (ölçüm çizginin hafif altında) —
  ya yasanın 3. mertebe düzeltmesi ya ikincil değişken

**DURUM: √w + (β/2)·v = 1 genlik-korunumu adayı 66 çift, 11 pencere,
6 asal, t=10⁵→2.7×10¹¹ boyunca ~%1.5 içinde tutuyor.** İlk-mertebe gölgesi
eski sum rule (u≈1); √ formu onu tüm ölçülen mertebelere taşıyor.

3. tur (yeni oturum): düşman-göz turu (dolaşıklık: w ve v aynı regresyon
verisinden — ortak-mod hata var mı?); teorik soru: NEDEN genlik-lineer?
(ünitarite şiddet-lineer verirdi; bu koherent genlik-bölüşümü) ve neden β/2?
Karar: Not 2'ye bölüm mü, Not 3 mü?

## YASANIN DÜŞMAN-GÖZ TURU (69) — ÜÇ SALDIRI, ÜÇ SAĞ ÇIKIŞ (18 Ağustos, gece)

S1 ORTAK-MOD (yarı-bölme, gürültüler bağımsız):
  aynı-yarı RMS 0.0175 (a=0.999, b=1.053) | ÇAPRAZ RMS 0.0193 (a=0.998, b=1.051)
  → neredeyse özdeş; ortak-mod hata açıklaması ÖLDÜ. Ve kesişim a = 1'e
  üç ondalıkta oturuyor — iki bağımsız yarıda da.

S2 ÜS ÖLÇÜMÜ: w = (1−bv)^k serbest fit → k = 1.938 ± 0.184 (2'den 0.3σ).
  Genlik-korunumu üssü VARSAYILMADI, ÖLÇÜLDÜ.

S3 FORM YARIŞI (w-ölçeği artık RMS): √w-lineer 0.0217 ≈ (1−bv)^k 0.0217
  < w-lineer 0.0230 < üstel 0.0372 << şiddet-ünitaritesi 0.1656.
  (w-lineer'in yakınlığı yanıltmasın: g-analizi v² teriminin gerçekliğini
  ayrıca kanıtlıyor — lineer form onu üretemez.)

**SIÇRAMA 1 DURUM: √w + (β/2)v = 1 üç ayrı düşman testinden geçti.**
Kalan: teori (neden genlik-lineer? neden β/2?) + yayın kararı (Not 2
bölümü / Not 3) — yeni oturum.

## 70 + FAZ KONTROLÜ: PERDE KAPANIYOR — VE AYNA-FAZINDA YENİDEN AÇILIYOR ★★★
(18 Ağustos, gündüz — "neden tam kapanmıyor?" sorusunun ölçülmüş cevabı)

Genişletilmiş ölçüm: p ∈ {2..31} × 12 pencere → 132 üçlü, τ = 0.056-0.613.

**BULGU 1 — w İŞARET DEĞİŞTİRİYOR:** w(τ) sıfırı τ₀ ≈ 0.395-0.40'ta kesiyor
(+0.018@0.378 → −0.004@0.402) ve negatifleşip ≈ −0.17'de platoluyor
(τ=0.5-0.61; 23 ufuk-bölgesi noktasının HEPSİ negatif — işaret deseni 2⁻²³).
**τ₀ = τ*: tükeniş analizinin rol-değişim noktası, w'nin sıfır geçişiymiş.**

**BULGU 2 — FAZLAR DÖNMÜYOR, ATLIYOR:** w fazı 0° → (kesişte) → ±180°;
sin bileşenleri her yerde ≤%3. v fazı hep ~0° (±3°), |v| pürüzsüz
(tepe ~0.60 tam kesişme civarı, sonra 0.43'e iniş). Sönümlü rezonans ±90°'den
geçerdi — geçmiyor. **Her iki kanal da REEL: sistem KAYIPSIZ.**
Kayıpsız (Hermitsel) sistemin reaktif yanıtı: iç özfrekans geçilirken
işaret değişir, soğurma yok. Hilbert-Pólya kokusu: perde yutmuyor, YANSITIYOR.

**BULGU 3 — RS-ufku hipotezi (P1) YANLIŞLANDI:** τ=1/2'de özel bir şey yok;
kapanış 0.40'ta, ufuktan önce. (P2 doğru: v 1/2'yi pürüzsüz geçiyor.)

**SORUNUN CEVABI:** Perde tam kapanıyor — τ* ≈ 0.40'ta — ve ayna-fazında
yeniden açılıyor. √w + (β/2)v = 1 yasası kesiş-öncesi dalın denklemi;
v'nin 0.935'e ulaşması hiç gerekmiyordu. τ* böylece yasaya bağlandı ✓.
V_res bağlantısı hâlâ açık (negatif dalın varyans katkısı w² ile girer —
sonraki hesap).

Dürüstlük: negatif-w bölgesi alçak pencerelerden (n=1k-25k); işaret deseni
ve 5-8σ genlikler sağlam ama per-asal yarılma payı hesaba katılmalı;
τ₀'ın tam yeri ±0.02.

## 71: V_res'İN KİMLİĞİ ÇÖZÜLDÜ — ÜÇ BİLEŞEN (18 Ağustos, öğle)

(A) KUYRUK ASALLARI (p>300, ölçülmüş w(τ) ile hesap): V_res'in yalnız
    ~%4.5'i (0.008/0.163, 0.011/0.244). "Soyulmamış kuyruk" açıklaması ÖLDÜ.
(B) BOŞLUK-MODÜLASYONLU ASAL KANALI (g̃-etkileşimli soyma, plasebo-düzeltmeli):
    NET 0.074→0.062 (L=9.86→12.45) — V_res'in %45→%25'i. 43'ün w(g̃)
    bulgusunun varyans bütçesindeki karşılığı. r* bu soymayla 0.984'e çıkıyor.
(C) ÇEKİRDEK (artakalan): 0.081 → 0.171, eğim ≈ 0.035/L —
    **CUE'nun genlik-varyans hızı 0.036/N ile NEREDEYSE AYNI, genlik ~1/10.**
    Ham V_res'in 0.031/L'si bileşenler ayrılınca 0.035/L'ye oturdu:
    gazın öz-gürültüsü RMT hızında büyüyor, RMT genliğinin onda birinde.

**BÜYÜK RESİM TAMAM (bugünkü soru kapandı):**
  yasa (√w+βv/2=1, kesiş-öncesi dal) + τ* (w'nin sıfır geçişi; ayna-açılım;
  kayıpsız) + V_res (= %5 kuyruk + %25-45 modüle asal + CUE-hızlı öz-çekirdek)
  — üçü tek yanıt-teorisi resminde. NOT 3 İSKELETİ HAZIR:
  "The response theory of the zero gas".

Dürüstlük: kuyruk hesabı w-eğrisi modeline bağlı (±%50 oynasa pay ≤%7);
çekirdek eğimi 6 noktadan; modüle-kanal payının L ile düşüşü ayrıca ilginç
(mutlak sabit ~0.06-0.07 → payı düşüyor) — teorik iş.

## 72: SICAK KRİSTAL (DEBYE-WALLER) — AYNA KAYMASI HESAPLANDI ★★★
(18 Ağustos — "çocuk resmi" seansının kapısı)

Bragg resmi: sıfır örgüsü = kristal; ideal ayna τ=1/2 (RS ufku = "2 salıncak");
ısıl titreşim aynayı kaydırır. TÜM girdiler ölçülü: jitter σ_u (T1),
plato B=−0.161, ideal kesiş 1/2.

**T1 — Kristal ISINIYOR:** unfold jitter sabit değil: 0.198 (L=5.6) →
0.272 (L=12.45) — S(t)'nin bilinen log-ısınması. c_jitter: 0.77→1.46.

**T2 — SIFIR-AYAR VURUŞU:** c=jitter, B=plato, tek serbest A →
A = 1.03 ≈ 1 (TAM koherent genlik) ve modelin sıfır-geçişi **τ = 0.406**
(ölçülen τ* = 0.395-0.40). AYNANIN 0.50→0.40 KAYMASI, ölçülmüş sıcaklık +
dolaylı kanal ile HESAPLANDI.

**V2 — İKİNCİL DEĞİŞKENİN KİMLİĞİ:** pencere-bazlı c(L)·s (tek ortak s):
RMS 0.052→0.023 ve soğuk/sıcak pencere artık-yarılması +0.034/+0.040 →
±0.003'e ÇÖKTÜ. Günlerdir izlediğimiz "sabit-τ'da L-yarılması" ikincil
değişkeni = PENCERE SICAKLIĞI. (Üçüncü anomali de resme katıldı.)

**Açık kalanlar (dürüst):** etkin sönüm ölçeği s ≈ 3.3 (naif bağımsız-DW'nin
üstü — korele düzensizlik/öz-tutarlı alan teorisi ister); yapı çarpanı
(1−2τ) kaba tahmin (şekil artığı 0.023 > taban 0.012-0.016); T3'te serbest-c
fit 4.2 istiyor (s·c_ort ≈ 3.9 ile tutarlı).

**GÜNÜN CÜMLESI:** "Sıfırlar, asal dalgalarını Bragg-yansıtan, logaritmik
ısınan bir kristaldir; aynanın yeri kristalin sıcaklığından hesaplanır;
teli pencereden pencereye farklı gösteren şey sıcaklık farkıymış."
Not 3 iskeleti artık: yasa + kayıpsız Bragg aynası + sıcak-kristal kayması +
V_res ayrışımı + sıcaklık-ikincil-değişken.

## 73: s'NİN TEORİSİ ÇÖZÜLDÜ — DURAN DALGA SÖNÜMÜ ★★★

Hipotez: Bragg yansımasında kristal içi alan 2k periyotlu DURAN DALGA →
sönüm |⟨e^{i2ku}⟩| (Gauss'ta s=4). AMPİRİK test (Gauss varsayımı yok,
karakteristik fonksiyon doğrudan u-serisinden):

**s_eff = ln φ₂ / ln φ₁ = 4.11 ± 0.04** — 72'nin fitlerinin istediği 3.3-4.2
aralığının tam içi. Gizemli çarpan = duran-dalga (2k) faktörü; Gauss 4'ü,
alt-Gauss u-dağılımı (kurtosis −0.3…−0.5 — seviye itmesinin kuyruk bastırması!)
4.11'e giydiriyor.

MODEL YARIŞI (yalnız A serbest, B=plato): sönümsüz 0.074 → naif DW 0.052 →
**duran dalga 0.029**. Kesişler: naif 0.408 / duran 0.368 (ölçülen 0.395-0.40
arada) → kısmi-yansıma karışımı hipotezi (Bragg altında yansıma kesirli;
Φ = (1−r)φ₁ + rφ₂) doğal sonraki adım.

Ek ölçümler: ρ₁(komşu-u) L ile 0.40→0.63 yükseliyor; u alt-Gauss.

**ÇOCUK-RESMİ SEANSININ TOPLAM HASADI (70-73):** ayna = Bragg; kayma =
sıcaklık (hesaplandı, 0.406); ikincil değişken = pencere sıcaklığı;
s = duran-dalga çarpanı (ölçüldü, 4.11). Not 3 omurgası 6 omur:
yasa + kayıpsız ayna + sıcak-kristal + duran-dalga + V_res ayrışımı +
sıcaklık-değişkeni. Açık: yapı çarpanı (1−2τ) kaba (RMS 0.029 > taban),
kısmi-yansıma karışımı, s'nin 4'ten 4.11'e mikro-sapması.

## 74: KAT KAPISI — TÜRETİM + BÜYÜK SÜRPRİZ: AYNA-TOPLAMI KONUŞUYOR ★★★
(18 Ağustos, akşam)

**TÜRETİM (3 satır, kalıcı):** |Z|²'de cos(t log p) frekansını (m, mp)
çiftleri üretir (mp ≤ N) → bağıl genlik 2p^{-1/2}·S₁(N/p)/S₁(N) →
sürekli limitte 2(1−2τ), p>N'de TAM SIFIR. (1−2τ) yapı çarpanının kimliği:
katlanmış toplamın çift-sayma aritmetiği. ✓

**ÖLÇÜM (noktasal |Z|², 1.27M grid, L=9.86, parametresiz kıyas):**
1. Fold altında ölçülen R tahminden SİSTEMATİK YÜKSEK (1.32 vs 0.67 @τ=0.38)
2. **FOLD ÖTESİNDE R ≈ 1.0 PLATOSU** (p=139→499, τ=0.50→0.63; tahmin 0!)
   — cos-kilitli (sin ≤ 0.04), 20σ+ sıfırdan uzak.

**TEŞHİS:** "hızlı" diye ihmal ettiğim 2θ-chirp terimleri (cos(2θ − t log nm))
tam nm ≈ N²/p'de log-p frekansına KİLİTLENİR (durağan-faz geçişi) — bunlar
AFE'nin AYNA-TOPLAMI ile ana toplamın çapraz terimleri. Fold ötesinde direkt
çift yok → görülen HER ŞEY ayna-kanalı: **χ-girişimi ilk kez noktasal olarak
ölçüldü, büyüklüğü ~1·p^{-1/2}, düz plato.** Max-kanalındaki B = −0.16
platosunun kimlik adayı da bu (işaret/büyüklük eşlemesi açık iş).

**Resim derinleşti:** "ayna = fonksiyonel denklem" çürümedi — tersine:
katın ötesinde görünen her şey aynadaki görüntünün kendisi çıktı.

SONRAKİ (yeni oturum): (a) ayna-kanalı tahmini: Σ_{nm≈N²/p}(nm)^{-1/2}
bölen-toplamları + durağan-faz ağırlığı → R_öte ≈ 1 hesaplanabilir mi?
(b) B'nin işaret/büyüklüğünü ayna-kanalından türet; (c) fold-altı fazlalık =
direkt + ayna toplamı ayrışımı; (d) Not 3'e "the two-sum interference" bölümü.

## 75: AYNANIN SESİ HESAPLANDI — ÇADIR ŞARKISI ★★★ (18 Ağustos, gece)

Tahmin: noktasal |Z|² modülasyonu kaydırılmış-ikinci-moment (BCHB-tipi)
ana terimini izler: R = 2(L−log p + c₀)/(L+c₀) ≈ 2(1−τ) — katta KESİLMEZ,
tam ufukta (p≈t/2π) ölür. Ayrışım: ayna = ölçüm − direkt = ÇADIR
(2τ ↑ katta 1.0 ↓ 2(1−τ)).

SONUÇ (2 pencere, 44 nokta, parametresiz):
- Şekil DOĞRULANDI: iki pencerede de ölçüm 2(1−τ)-ailesini kat boyunca
  izliyor (p=2'de %0.4-0.7 uyum; fold-ötesi 0.92/0.77 vs tahmin ~0.77);
  çadır deseni ayrışımda net görünüyor.
- Sistematik: ölçüm tahminden ~%8-15 yüksek; c₀ serbest fit 1.12 (tahmin
  2γ−1=0.154) — kesin sabitler GERÇEK BCHB ana terimlerini ister (χ-tarafı
  ikinci ana terim dahil). Teorik ödev: sabitin türetimi.
- UYARI: τ ≳ 0.7'de sıfır-tarağı bandı (ω ~ 2π/gap dağılımı) regresyonu
  kirletiyor — p=21169@L=12.45 (τ=0.80) ölçümü 5.2 (tahmin 0.42): tarak
  kontaminasyonu, formül testi τ ≲ 0.7 ile sınırlı.

**KAT KAPISININ TOPLAM HASADI (74-75):** (1−2τ) türetildi (çift-sayma);
ayna kanalı keşfedildi ve şarkısı hesaplandı (çadır, kaydırılmış-moment
ailesi); fonksiyonel denklem yanıt-teorisine ÖLÇÜLEBİLİR biçimde bağlandı.
Not 3 omurgası 8 omur oldu: + iki-toplam girişimi + çadır şarkısı.
Açık: kesin sabit (BCHB türetimi), tarak-bandı ayrıştırması, B'nin
max-kanalındaki işaret/büyüklük eşlemesi.

## NOT 3 YAZILDI: "A response theory for the Riemann zero gas" (19 Ağustos)

`arxiv_response_theory.tex` — 8 omur tek çatıda: (1) genlik-korunumu yasası
(3 kill-test metinde), (2) işaret geçişi + kayıpsızlık (fazlar atlar),
(3) termal Bragg (kesiş 0.406 hesaplı; ikincil değişken = sıcaklık),
(4) duran-dalga s=4.11, (5) (1−2τ) türetimi (Derivation ortamında),
(6) ayna kanalı = çadır (2 pencere, BCHB-ailesi), (7) V_res ayrışımı
(CUE-hızlı çekirdek), (8) sentez + 6 açık problem. 3 İngilizce figür (76).
Derlendi. GÖZDEN GEÇİRİLECEK: BCHB 1985 künyesi (J. reine angew. Math. 357,
161-181) hafızadan — gönderim öncesi teyit şart; c₀ ~1.1 vs 2γ−1 dürüstçe
metinde. Üçleme tamam: Not 1 (gözlem) + Not 2 (anatomi) + Not 3 (teori).

## NOT 3 DÜŞMAN-GÖZ DENETİMİ + DÜZELTMELER (19 Ağustos)

Bağımsız ajan denetimi (Not 1-2 protokolüyle): ~50 sayı yeniden üretildi,
BCHB 1985 künyesi web'den TEYİT (J. reine angew. Math. 357 (1985), 161-181).
Karar: "omurga sağlam, engeller metin düzeyinde — revizyonlarla gönderilebilir."

3 KIRMIZI (hepsi işlendi):
- K1: "31 asal / 168 nokta" → gerçekte 11 asal / 132 nokta (metin + figür
  etiketi düzeltildi, fig_crossing_en yeniden üretildi).
- K2: "kesiş HESAPLANDI (0.406)" abartıydı — naif-c modeli fit olarak
  reddediliyor (serbest c=4.17, jitter'dan 28.9σ); duran-dalga en iyi fit ama
  0.368'de kesiyor; ara ölçekler 0.39. Yeni çerçeve: TERMAL AİLE KESİŞİ
  [0.37, 0.41]'DE KISKAÇLIYOR (abstract + Observation + giriş yeniden yazıldı).
- K3: "β faz muhasebesinden TÜRETİLDİ" savunulamaz — w(0)∈[0.88,0.98]
  çözülmemiş (Not 2), serbest eğim 1.029±0.012 ile β/2=1.07 arası ~3σ gerilim.
  Yeni dil: "motive edilmiş, türetilmemiş"; Not 2'nin β≈2.3'üyle köprü kuruldu.

8 SARI (hepsi işlendi): B taban-bağımlılığı + "sabit boşlukta" niteleyicisi;
serbest-üs dürüstlüğü (aşağıda); s_eff=4.11 "Gauss'luk ölçüyor" itirafı +
1.5k ara-ölçek eşdeğerliği; çadır "%6-20 sistematik altta" + tarak sınırı
τ≲0.62 + noktasal-regresyon plasebo eksikliği bayrağı; out-of-sample −0.009
tek-yanlı ofset; "one sixth/one tenth" köprüsü; abstract öncül düzeltmesi
(ayna kanalı = çadır); "lossless" → "consistent with lossless";
"independently" → "by a different method"; Berry88 + GHK07 metinde atıflandı;
Reproducibility 62-78.

YENİ SCRIPTLER (denetçinin S5'i — manşet sayıların numaralı scripti yoktu):
- 77_yasa_fiti.py: 44-çift yasa fiti. Serbest: √w=0.990−1.025v RMS 0.0096;
  a≡1: b=1.069; pencere-blok bootstrap eğim 1.029±0.012; form yarışı ×9.2.
  ÖNEMLİ TEŞHİS: ortak (b,k) fiti DEJENERE VADİDE (b=0.75,k=3.0 da aynı RMS;
  bootstrap medyan 2.9, %16-84 [2.3,3.9]) — ama b'yi doğrusal yasadan
  koşullayınca k=2.01 taş gibi. Metin buna göre: "üs tek başına O(1);
  koşullu k=2.01; ayırt edici test form yarışı" (denetçinin 1.94±0.29'undan
  da dürüst — o da kendi boru hattının vadi konumuydu).
- 78_sicaklik_cokmesi.py: pencere-c × ortak-s modeli: RMS 0.052→0.023,
  soğuk/sıcak artık yarılması +0.030/+0.044 → ±0.003. Not 3 V2 iddiası
  scriptle mühürlendi.

Not 3 yeniden derlendi (240 KB). Üçleme REVİZE HALİYLE TAMAM.
Kalan: kullanıcı okuması; arXiv mekaniği (hesap, endorsement, çapraz-ID);
teori ödevleri (çadır sabiti/BCHB ana terimleri, B eşlemesi, β türetimi,
r(τ) kısmi yansıma, çekirdek dağılımı GUE?, L-fonksiyon evrenselliği).

## 79: B'NİN EĞİMİ — ÇADIR YAMACI SORUSU ARKA ODANIN İKİ OKUMASINA VARDI (19 Ağustos)

Fısıltı (çizim seansı): B = −0.16 ayna-sızıntısıysa düz olmamalı, çadırın
yamacını (−k·2(1−τ)) izlemeli. Gece dört perdede ilerledi:

1. SIMPSON TUZAĞI: havuzlanmış kat-ötesi eğim negatif (−0.24±0.09) çıktı ama
   yüksek-τ noktaları yalnız küçük-L pencerelerinden geliyor; L-kontrollü eğim
   +0.18±0.11 (pencere-içi +0.22, +0.18) — çadır yönünde ama düzle de uyumlu.
   BİÇİM ÇÖZÜLMEDİ. Termal-çadır (DW×tent) fiti anti-sönüm istiyor (s<0) →
   basit termal-çadır ÖLÜ.

2. PLASEBO DERSİ: p+0.5 sahteleri asal çizgisine 3-16 çözünürlük-birimi
   mesafede — SIZINTI (+0.07 yalancı taban, işaret bile ters). Geometrik-orta
   sahtelerle taban ±0.02-0.05; sinyalin ~%5-15'i. Sahte tasarımında kural:
   en yakın gerçek çizgiye ≥20 çözünürlük-birimi.

3. "ERİYEN PLATO" ARTEFAKTI: pencere-ortalaması B̂(L) monoton görünüyordu
   (−0.224→−0.143, dB/dL=+0.039, 6σ!) ve L≈11.7'de sıfır öngörüyordu —
   büyük-pencere testi bunu ÖLDÜRDÜ (B̂≈−0.35, L-bağımsız). Sebep: pencere
   ortalamalarının τ-bileşimi L ile kayıyordu. AYNI bandda ([0.505,0.55])
   ölçünce süreksizlik yok.

4. ANA BULGU — ARKA ODANIN İKİ OKUMASI, İKİSİ DE EVRENSEL: aynı τ-bandında
   ÇIPLAK okuma (yalnız band asalları tabanda) B = −0.344; GİYDİRİLMİŞ okuma
   (küçük asallar kovaryat) B = −0.124. Oran 2.77. İKİSİ DE L=5.6→16.6
   boyunca SABİT (36+41+55 veri setleri, 12+1 pencere; 55@1e8 −0.32±0.11
   aynı-motor kontrol). Not 3'ün −0.16'sı kendi konvansiyonunda (P11 karışık
   band) tutarlı, arada. Giydirme mekanizması AÇIK — küçük-asal kovaryatları
   kat-dibi yanıtının ~%64'ünü emiyor; 71'in gap-modüle etkileşim kanalıyla
   akrabalık ilk şüpheli. Not 3 metnine dokunmadı (taban-bağımlılığı uyarısı
   zaten girmişti; şimdi UÇLARI SAYISALLANDI).

Açık: (a) biçim sorusu için büyük pencerelerde τ∈[0.55,0.75] gerekiyor
(L=16.6'da p≈9e3..2.6e5 — ağır ama olanaklı); (b) giydirme mekanizması
(çapraz-frekans etkileşim modeli); (c) Odlyzko bağımsız-motor teyidi
(sinyal/gürültü ~6σ tahmini). Script: 79_B_egimi.py (+ para-grafiği
79_B_egimi.png: iki düz çizgi, −0.34 ve −0.12).

## 80: GÖRÜNTÜ ÇİZGİLERİ — 2θ SPEKTROSKOPİSİ İLK ÖLÇÜM (19 Ağustos)

Fısıltı: yansıyan dalganın frekansı aynalanır (görüntü kaynak, kayan frekans
L−log m) — düz kulak duyamaz, 2θ gözlüğü gerekir: cos(2θ(t)−t·log m).
Bunlar AFE ayna-toplamı çapraz terimleri (74'ün chirp'leri) — ilk kez
max/gap KANALLARINDA arandı (13 pencere, m=1..12 + tamsayı-olmayan plasebo).

ÜÇ PERDE:
1. İlk geçiş GÜR: tamsayı-m çizgileri χ²=2700'e dek (52σ), HEPSİ saf-cos
   kilitli (sin≈0); plasebo tertemiz (⟨χ²⟩=1.44); v kanalı neredeyse sağır
   (18 vs plasebo 8.6 — explicit formülde 2θ yok öngörüsüyle uyumlu).
   u(m) deseni çarpıcı: 2-kuvvetleri (4,8) gür, 2×tek (2,6,10) SESSİZ.
2. ÖRNEKLEM KİMLİĞİ UYARISI: gap-ortalarında θ(γ)≈πn → 2θ gözlüğü kısmen
   S-giydirilmiş DİREKT çizgiye katlanır. Ve P11 tabanında asal-kuvvet YOKTU:
   düz kolonlar eklenince m=4,8 sesinin çoğu 2²,2³ direkt çizgisi çıktı
   (+0.148/+0.051 — sum-rule'daki |u|≈1 ile tutarlı yönde). DERS: taban
   daima p^k içermeli; 79'un "giydirme" bilmecesine de şüpheli oldu.
3. AYRIŞTIRMA SONRASI SAHİCİ 2θ-ÖZEL İÇERİK HAYATTA: m=4: +0.0092±0.0004
   (20σ+), m=8: +0.0018±0.0004 (5σ), oranlar direkt çizginin %3.6-6.2'si —
   KISMİ YANSIMA ADAYI r ~ %4-6. m=1 (saf 2θ/S çizgisi, ilk geçişte +0.041)
   kimliği AÇIK: düz karşılığı DC'ye çöker, temiz test S-vekili ister.

Açık: (a) m=1 S-çizgisinin temiz tasarımı; (b) kat-dibi bölgesinde
direkt+görüntü karışımı (inst. frekans çakışması — 79'un çıplak −0.35'iyle
bağ); (c) r(τ′) yasası için m-taraması p^k-tam tabanla; (d) v'nin zayıf
2θ artığı (plasebo-üstü ~2×) ne? Script: 80_goruntu_cizgileri.py + png.

## HAREKAT PLANI (19 Ağustos gecesi yazıldı; 20 Ağustos için)

BUGÜNÜN HASADI (dört commit): Not 3 denetim düzeltmeleri işlendi + 77/78
(4d6ce8b); üç notun okuma rehberi + üç yeni fısıltı (çizim seansı);
79: arka odanın iki evrensel okuması (b17d319); 80: görüntü çizgileri,
sahici 2θ-özel içerik %4-6 (34f08fb).

İKİ KOL + BİR KAPI:

KOL 1 — YAYIN (üçleme arXiv'e):
  1a. Kullanıcının son okuması (devam ediyor; rehber verildi).
  1b. arXiv hesabı + math.NT endorsement DURUMU — DIŞ BEKLEME SÜRESİ VAR,
      erkenden başlatılmalı (kullanıcı işi; metadata/abstract'ları ben
      hazırlarım). Gönderim sırası önerisi: Not 1 → ID al → Not 2-3'e
      çapraz-ID işle → aynı oturumda gönder.
  1c. KAPI (ön-gönderim sigortası): TABAN TAMLIĞI DENETİMİ (aşağıda 2a).
      Not 1 kanal-regresyonu içermiyor → muaf. Not 2-3 sayıları p^k-eksik
      tabanla ölçüldü → denetim geçilmeden gönderilmez.

KOL 2 — ARAŞTIRMA (79-80'in açtığı damar; Not 4 adayı):
  2a. ÖNCE: 81 = taban tamlığı denetimi. 80'in dersi: P-tabanlarında asal
      kuvvet yoktu (log4, log8 direkt +0.148/+0.051 taşıyor!). Standart
      boru hattına p^k (4,8,9,16,25,27,32) kolonları ekle; üçlemenin manşet
      sayıları (τ-yasaları, korunum yasası fiti, sum rule, soyulmuş r*,
      B'nin 2×2'si) OYNUYOR MU? Oynamazsa → üçleme mühürlü, gönder.
      Oynarsa → düzelt, yeniden derle, sonra gönder. AYRICA 79'un
      "giydirme ×2.77" bilmecesinin 1 numaralı şüphelisi bu.
  2b. Kat-dibi karışımı: fold civarında direkt(p^k-tam) + 2θ ortak taban —
      çıplak −0.35 ayrışıyor mu (direkt kuvvetler + görüntü)? 79↔80 köprüsü.
  2c. r(τ′) taraması: p^k-tam tabanla m-süpürmesi; yansıma yasası biçimi,
      pencereler arası çökme.
  2d. m=1 S-çizgisi temiz tasarım: 61'in RvM-pürüzsüz konum makinesiyle
      S(t_mid) vekil kolonu kur; genlik-S çiftlenimini doğrudan ölç.

BEKLEME LİSTESİ (sıra sonrası): geçişin bulanıklığı (termal 2. öngörü);
tepe-konumu kanalı (kemerin biçimi); Odlyzko bağımsız-motor B teyidi
(çapalı-θ mpmath tasarımı gerek); v'nin zayıf 2θ artığı; teori ödevleri
(çadır sabiti/BCHB, β türetimi, çekirdek GUE?).

YARIN SABAH BAŞLANGIÇ SIRASI: (1b'yi kullanıcı başlatır — bekleme süresi
yüzünden ilk iş) → 2a (kapı; yarım gün) → sonuca göre ya gönderim
mekaniği ya 2b. Gerekçe: 2a hem sigorta hem bilmece çözümü — çift verim.

## GECE UÇUŞU HARİTASI — DURUM DEĞERLENDİRMESİ (19 Ağustos, gece sonu)

Ülkenin haritası çizildi (sohbette; iki resim: "son gece uçuşu" + "Riemann
ülkesi haritası"). Dürüst kartografya:

AYDINLIK (fenerli) BÖLGE — ölçüldü, savunuldu, denetimden geçti:
Kemerler Vadisi (Not 1: r, +43σ null, N_eff kayması/yarılması), İki Kanal
Ovası (Not 2: w/v, τ-yasaları, sum rule; v 18 mertebe), Korunum Tepesi
(√w+βv/2=1, 3 kill-test; ~3σ eğim gerilimi bilinen çatlağı), Sıcak Ayna
(kesiş [0.37,0.41] kıskacı), Çadır Sahili (noktasal ayna kanalı).

SİS KUŞAĞI — ölçüldü ama anlaşılmadı:
Arka Oda (iki evrensel kapı −0.34/−0.12; giydirme mekanizması ?),
Görüntü Kıyısı (%5 sahici yankı; yasası ?), S-feneri (m=1, kimliği ?),
w₀ ([0.88,0.98] açık), ve KAT DİBİ SİS BOĞAZI — iki gecenin bütün
iplerinin işaret ettiği yer; yarınki seferin rotası oraya.

KARANLIK KUZEY: Çekirdek Denizi (V_res özü: GUE mi?), Neden Dağları
(türetimler: β, çadır sabiti, B kimliği, kısmi yansıma teorisi).
En uzakta RH Zirvesi: patikamız oraya ÇIKMIYOR ve bunu iddia etmiyoruz;
ama ölçtüğümüz her şey (reel yanıtlar, dönmeden atlayan fazlar, kayıpsızlık)
Hermitsel bir üreteçle TUTARLI — bulut arasından silüet.

MESAFELER: yayına günler; fenomenoloji→teori orta mesafe (aylar, patika
belli); RH'ye çok uzak ve rota o değil. HAVA DURUMU: iyi ve güzel —
yasalar tek değişkende çöküyor, plasebolar temiz dönüyor, ve bugüne dek
her tuzak doğanın değil BİZİM artefaktımız çıktı; düzeltilince resim hep
netleşti. Epsikl bataklığı böyle davranmaz; gerçek ülke böyle davranır.

## 81: KAT DİBİ SİS BOĞAZI GEÇİLDİ — GİYDİRME MEKANİZMASI ÇÖZÜLDÜ (20 Ağustos)

Boğaza girdik; sis beklediğimizden başka bir şey saklıyormuş:

1. GİYDİRME MERDİVENİ (T1): band-asal B̂'si taban katmanlarıyla —
   yalnız band −0.354 | +P11 −0.126 | +YALNIZ p^k kuvvetleri −0.358
   (KUVVETLER MASUM — 80'in şüphelisi aklandı) | +tam direkt taban
   (tüm p^k ≤ e^{0.45L}) −0.039 | band-altı TAM (≤ e^{0.505L}) −0.042
   → YAKINSADI. Kontrolsüz sütun: +1.56→+0.90 (işaret bile kontrollerin
   eseri; Not 3'ün "sabit boşlukta" niteleyicisi kat dibinde ×4 önemli).

2. MEKANİZMA KANITI (T3): band-kolonu ~ küçük-asal kolonları Gram-R²:
   gerçek gap-ortalarında 7×10⁻⁴, RvM-pürüzsüz ızgarada 0.00000.
   → ÖRNEKLEME IZGARASI ASAL DALGALARINI TAŞIYOR (v-alanı deseni
   örgüde). Atlanmış-değişken transferi + p^{-1/2} normalizasyonunun
   √p büyütmesi = giydirme. Büyüklük mertebesi tutuyor (ρ~0.006 ×
   büyük küçük-asal katsayıları × √p ≈ 0.2). ÇIPLAK −0.35'İN ~%88'İ
   HAYALETMİŞ. 79'un "iki okuması" tek resimde birleşti: her taban
   kendi atlanmış-değişken yüküyle okur; fizik tam-taban okumasında.

3. SAHİCİ ARKA ODA (T4): B_tam(L) ≈ −0.05 (−0.06…−0.026), 18σ gerçek,
   hafif iniş eğilimli olabilir, 36/41 veri sınırında sıçrama YOK.
   Not 3 iması: metnin taban-bağımlılığı uyarısı artık MEKANİZMALI;
   gönderim öncesi uyarıya tek cümle güçlendirme önerilecek (platonun
   fiziksel değeri konvansiyon değerlerinden hayli küçük).

4. GÖRÜNTÜ REGRESÖRLERİ (T2, ilk geçiş): durağan-bant ayna toplamı
   zayıf ama var (cos +0.0026±0.0005, 5.5σ, saf-cos); kaydırılmış-bant
   ve yarım-ızgara kontrolleri null ÇIKMADI (30σ/13σ) — 2θ-toplam
   ailesi ortak S-kaynaklı bileşen paylaşıyor; kontroller aslında ayna
   sürekliliğinin ölçümü. v kanalı sessiz (öngörü ✓). T2 v2 tasarımı:
   karıştırılmış-gap vekil ızgara + S-vekili kolon gerekiyor.

Ders defteri: (a) kat dibinde hiçbir tek-taban okuması ham haliyle
fiziksel değil — merdiven yakınsaması şart; (b) örnekleme deseni bir
kovaryat gibi davranır; (c) "kontrol" diye eklenen her kolon önce
kendisi ölçülmeli. Script: 81_kat_dibi.py (T1-T4).

## 82: KAPI TESTİ ALARM VERDİ — ÜÇLEME BEKLEMEYE ALINDI (20 Ağustos)

81'in mekanizması "büyük sonuç" olarak kovalandı ve kapı görevini yaptı:

1. YER-GERÇEĞİ KANITI (T4b): gerçek gap-ortası ızgarası + BİLİNEN sentetik
   sinyal (79 çizgi + kontroller + gürültü): P4-tabanlı regresyon w'yi
   %15-25 SİSTEMATİK DÜŞÜK ölçüyor (21σ); tam taban gerçeği ~%1 içinde
   buluyor. Kısıtlı-taban kanal ölçümü ızgara-çiftlenimi yüzünden yanlı;
   FİZİKSEL OKUMA = TAM-TABAN OKUMASI. (Öngörümün aksine küçük τ bağışık
   DEĞİL: R² küçük τ'da yeterince düşmüyor, T3.)

2. GERÇEK VERİDE ETKİ: korunum yasası fiti P4→TAM: a 0.990→1.017,
   b 1.025→0.884 (~12σ), a≡1 eğimi 1.069→0.806; RMS 0.0096→0.0075 —
   YASANIN FORMU GÜÇLENIYOR, KATSAYILARI DEĞİŞİYOR. Merdiven TAM/2→TAM
   arasında yakınsıyor. 132 τ-yasası noktasının hepsi oynuyor (ort ~15σ,
   %2-6 mutlak); Δv/v %13-20.

3. KARAR: Not 2-3 manşet sayıları kısıtlı-taban konvansiyonuyla ölçülmüş
   → GÖNDERİM DONDURULDU; "büyük yeniden ölçüm" kampanyası gerek (tüm
   boru hattı tam tabanla: τ-yasaları, sum rule, soyulmuş çekirdek,
   kesiş, s_eff, w₀). Not 1 MUAF (kanal regresyonu yok). Not: b~0.88,
   a~1.02 — β/2=1.07 hikâyesi ve K3 gerilimi yeniden yazılacak; w₀ tam
   tabanda 1'e yaklaşabilir (perde tam saydam?) — heyecan verici ama
   önce ölçüm. Süreç kültürün zaferi: bu, arXiv'den SONRA yakalansaydı
   erratum olurdu; kapıda yakalandı.

Scriptler: 82_izgara_denetimi.py (T1-T4). Sıradaki: 83 = büyük yeniden
ölçüm (tam-taban boru hattı, tüm manşetler) + Not 2-3 revizyon turu.

## 83: BÜYÜK YENİDEN ÖLÇÜM TAMAM — ESKİ→YENİ TABLOSU (20 Ağustos)

Tam-taban kampanyası (tüm p^k, τ≤0.45-0.55 / q≤720; Odlyzko dahil 11+12
pencere). Sonuçlar üç sınıf:

GÜÇLENENLER:
- TERMAL BRAGG ÇAKIŞTI: yeni kesiş τ₀ = 0.447±0.005 (eski 0.40); naif
  termal model (ölçülen c_jit=1.17, yeni B, A=1.020) kesişi 0.443'te
  veriyor — ÖLÇÜMLE 1σ İÇİNDE. Eski "aile kıskacı [0.37,0.41]" gerilimi
  ERİDİ: düzeltilmiş veriyle model basitleşti ve isabet etti.
- τ-YASASI EVRENSELLEŞTİ: kuvvet çizgileri u(4)=0.700, u(8)=0.560,
  u(9)=0.530, u(25)=0.338, u(27)=0.335 — her biri kendi τ'sundaki w
  değerinde: yasa asal DEĞİL, asal-kuvveti yasası (sum rule p^k dahil).
- Kayıpsızlık tutuyor: sin/cos ≤ 0.08 (geçiş altı), 0.14 (ötesi).
- w₀ izi: w(2) τ=0.028'de 0.9097±0.0035; kaba uzatma w₀≈0.97-0.98 —
  perde eski tahminden çok daha saydam (eski aralığın üst ucu).

DEĞİŞENLER (revizyon gerektirir):
- Yasa fiti: √w = 1.0171 − 0.8835·v, RMS 0.0075 (eski 0.990−1.025v,
  0.0096). Blok-bootstrap eğim ±0.021. a≡1 artık ZAYIF (RMS 0.0102 >
  0.0075): kesme 1.017, 1'den ayrılıyor. β/2 anlatısı yeniden yazılacak.
- FORM YARIŞI BULANIKLAŞTI: w-lineer 0.0099 ← √w-lineer 0.0127 (hafif
  öne geçti); şiddet-ünitaritesi hâlâ ölü (×6.3). "Genlik korunumu"
  formu artık tek galip değil — dürüstçe iki-form belirsizliği yazılacak.
- Out-of-sample p=11,13: serbest doğruya RMS 0.0269, ofset −0.018 (eski
  0.0148) — yasa out-of-sample hayatta ama daha gevşek; τ-menzil
  büyümesi/eğrilik şüphesi not edildi.
- Plato: B = −0.05..−0.09 (konvansiyona göre; eski −0.16'dan çok sığ;
  81-83 konvansiyon farkı AÇIK madde).

DEĞİŞMEYENLER:
- Soyulmuş çekirdek: ham 0.920 → P11-soyulmuş 0.9783 → TAM-soyulmuş
  0.9774 — r* tabana DUYARSIZ ✓ (Not 2'nin en sağlam sonucu).
- s_eff = 4.11 muaf (regresyonsuz); Not 1 komple muaf.

Çıktılar: 83_buyuk_yeniden_olcum.py, 83_tam_taban_egri.npz, 83_kampanya.png.
SIRADAKİ: Not 2-3 revizyon turu (yeni sayılar + iki-form dürüstlüğü +
basitleşen termal bölüm), ardından taze düşman-göz denetimi, SONRA arXiv.

## NOT 2-3 REVİZYONU TAMAM — TAM-TABAN SAYILARIYLA (20 Ağustos, akşam)

84: revizyon figürleri + Odlyzko tam-taban bataryası (fig_law/fig_crossing/
fig_tau yeniden; fig_tent ve fig_v değişmedi — çadır muaf, v yeniden-ölçüm
bekliyor). Odlyzko: w₂=0.9097±0.0035, yeni kılavuz 0.982−2.598τ farkı
+0.0011 (!); farklar τ ile −0.055'e büyüyor (eğrilik).

NOT 3 revizyonu: yeni abstract (sistematik + tam-taban konvansiyonu başa);
YENİ BÖLÜM "The sampling-grid systematic" (mekanizma, merdiven, yer-gerçeği,
konvansiyon, muafiyet sınıfları); §yasa → "The two-channel constraint"
(1.017−0.884v; kuvvet-u'ları; form-yarışı dürüstlüğü; β-koincidansı
"cautionary record" olarak kayda geçti); §kesiş 0.447±0.005 + plato
−0.05..−0.09 + faz oranları; §termal "The crossing is REPRODUCED"
(model 0.443 vs ölçüm 0.447; sıcaklık-çökmesi ve duran-dalga yarışı eski
konvansiyon bayrağıyla; s_eff duruyor); çadır bölümüne muafiyet cümlesi;
V_res bölümü r*-duyarsızlığıyla yeniden; sentez + açık problemler (7 oldu:
+ızgara-çiftlenimi teorisi); repro 62-84.

NOT 2 revizyonu: abstract'a sistematik + taban/çerçeve ayrımı; konvansiyon
bölümüne revizyon alt-bölümü; kanal tablosu yeni (w 0.83..0.38, v 0.11..0.37);
β≈2.3 cümlesi → kısıt yasası göndermesi; sum rule'a ortak-geçirim paragrafı;
r* taban-duyarsızlığı; τ* ↔ 0.447 bağlantısı; V_res bayrağı; τ→0 alt-bölümü
yeniden ("Near-complete transparency": eski %2-12 emilimin çoğu ızgara
yanlılığıymış, yeni iz w(0)≈0.97-0.98, beş-model yeniden-uzatması AÇIK);
Odlyzko maddesi yeni sayılarla; 10²¹-²² bölümü iç-tutarlı-eski-konvansiyon
notuyla; ikinci-derece çizgilere aday-mekanizma; açık problemler güncel;
repro 27-84. Her iki PDF derli.

KALAN yeniden-ölçüm listesi (dürüst bayraklar metinde): beş-model w₀;
exhaustion eğrisi; V_res oranları + etkileşim payı; sıcaklık-çökmesi;
duran-dalga yarışı; v-yasası 10²¹-²² + onset; ikinci-derece çizgiler.
SONRA: taze düşman-göz denetimi → kullanıcı okuması → arXiv.

## REVİZYON DENETİMİ + DÜZELTMELER (20 Ağustos, gece)

Taze düşman-göz denetimi (bağımsız ajan; 81-84 + 54 sıfırdan yeniden koştu,
npz'ler BİT-BİT AYNI çıktı — deterministik yeniden-üretim onayı). Karar:
"bu haliyle geçmez ama az farkla, tamamen onarılabilir." 5 KIRMIZI +
7 SARI — HEPSİ İŞLENDİ (28/28 yama):

- K1 (en değerlisi): u(4)=0.70 vs sum-rule u≈1 farkı taban değil KOŞULLAMA
  farkıymış (83'ün w-regresyonu g_u-kovaryatlı). 85_sum_rule_tam.py yazıldı:
  KOŞULSUZ tam-taban sum rule → ON BEŞ ASALIN HEPSİ u=0.988–1.007, DRİFTSİZ
  (eski +%1..12 drift taban transferiymiş); kuvvetler birimin hafif altında
  (0.83–0.99, ağırlıkla azalan) — YENİ ince yapı, açık bırakıldı. Not 2
  sum-rule bölümü koşulsuz/koşullu ayrımıyla yeniden yazıldı — SUM RULE
  ESKİSİNDEN GÜÇLÜ ÇIKTI.
- K2: silinen q=49 outlier'ı (1.18) geri kondu — dürüstlük gerilemesiydi.
- K3: termal formül birim hatası: c=(Lσ_u)²/2 → c=(2πσ_u)²/2 (sayılar
  zaten doğruydu, yazım yanlıştı).
- K4: tablo parantezi %5-20 → %8-30 (p ile büyüyen).
- K5: yer-gerçeği "%15-25" → "0.14-0.24 mutlak (%16-40 göreli)", "~%1" →
  "≤%1.5" (dört yerde).
- S1: kesiş tanım-sistematiği eklendi (±0.009 aralık; kuadratik-sıfır
  0.447-0.448 sabit). S2: termal "reproduced" dürüstlüğü (model sıfırı
  konvansiyon taramasında 0.44-0.47; B aynı eğriden → kısmen iç-tutarlılık);
  plato −0.04..−0.09. S3: 1.5σ → 1.6σ. S4: w(0) 0.94–1.01 model-bağımlı;
  açık problem (5) "0–0.06". S5: u(p^k)≈w(τ) "within errors" → "%1-2
  düzeyinde (q=4'te 2-4σ artık açık)"; koşullama dili netleşti. S6:
  yer-gerçeği sınırlaması (tek pencere/iid gürültü; genişbant bileşen
  dışlanmadı — 81-T2 bağlantısı). S7: N∈[5,19]; ≤0.011; onset 1.6-1.7
  (tanıma bağlı); kesme +%1.7 vs cap +%0.3 (statü açık); abstract V_res
  ve 10²¹⁻²² bayrakları; her iki not derli.

YEŞİL özeti: manşetlerin tamamı (yasa 1.017−0.884v ±0.021; kesiş; termal;
r*; Odlyzko; tablo; merdiven; Gram; cautionary record sayıları) birebir
doğrulandı. Denetçinin sözü: "bunlar işlendikten sonra iki not kendi
denetim kültürünün standardını karşılar." Script: 85_sum_rule_tam.py.

## 86: SIFIR ÖRGÜSÜNÜN KIRINIM DESENİ — CETVEL TEORİSİ İLK ÖLÇÜM (20 Ağustos, gece)

Fikir: bükülmenin tamamı tek nesnede — ızgaranın kendi Fourier'i
Ĝ(ω) = ⟨e^{iωt_n}⟩ (yapı çarpanı / kırınım deseni). Birinci-ilke teori
(t = t^s + u, birinci mertebe): üç parametresiz öngörü. ÜÇÜ DE İSABET:

P1 ASAL BENEKLERİ: |Ĝ(log p)| = v_p·p^{-1/2}/2 (83'ün v'siyle,
   parametresiz). Ölçüm/öngörü: p=2..7'de 1.07/1.05/1.01/0.98 (—
   küçük ω'da teori TAM), p büyüdükçe sistematik açık: 31'de 0.70.
   Ortalama 0.887±0.123. Açık ince yapı: benek-sönümü yasası (çizginin
   kendi DW'si ~0.90 veriyor, yetmiyor — v-doyumu/2. mertebe adayı).
P2 FAZ KİLİDİ: Ĝ(log p) SAF REEL, hepsi NEGATİF işaretli; Im ≈ 0.0000
   (Re ~0.03-0.06'ya karşı <1e-4!). v cos-kilitli ⟹ u sin ⟹ iω ⟹ reel:
   teorinin faz zinciri kusursuz doğrulandı. (Negatif işaretin fiziği —
   itme fazı — yorumlanacak.)
P3 TARAK/SICAKLIK: katlanmamış çerçevede |Ĝ_x(2π)| vs e^{-c_jit}:
   altı pencerede oran 0.88-0.89 SABİT — tarak parlaklığı pencere
   sıcaklığını birebir izliyor (sıcaklığın 3. bağımsız ölçümü ✓);
   sabit 0.885 çarpanı yeni açık sabit (alt-Gauss düzeltmesi adayı).
T4 KARANLIK ALAN: çizgi-dışı taban gerçek ızgarada 0.0002 — beyaz
   gürültü 1/√n=0.005'in 25 KAT ALTI ve karıştırılmış-gap vekilinin
   (0.0014) altında: HİPERUNİFORMLUK doğrudan kırınımda görüldü. S-ışıltısı
   parıltı değil KARANLIKMIŞ: sıfır gazının katılığı ekranda. Aritmetik
   içerik tamamen ayrık beneklerde.

SENTEZ: "cetvel bükülmesi" artık ne hata ne düzeltme — ölçülmüş, öngörülü
bir nesne. Ĝ(ω) tek başına birleştiriyor: giydirme (Gram çiftlenimleri =
benekler), v-alanı (benek parlaklığı = v'nin bağımsız ikinci ölçümü!),
sıcaklık (tarak), hiperuniformluk (karanlık alan). Not 4'ün omurga adayı.
Açık: benek-sönümü yasası; 0.885 sabiti; Re<0 işaretinin fiziği;
karanlık-alan tayfının ince haritası (S-denizi); kuvvet benekleri.
Script: 86_kirinim_deseni.py + 86_kirinim.png.

## 87: KAPALI DEVRE — BÜKÜLME, KRİSTALİN TERMAL DİFÜZ SAÇILMASINDA YAŞIYOR (20 Ağustos, gece-2)

Soru: 81'in ölçtüğü giydirme, 86'nın ölçtüğü Ĝ'den türer mi? ÜÇ ADIM:

1. KİMLİK (T1): iki dalga kolonunun Gram'ı = ½[Ĝ(Δω) ± Ĝ(Σω)] bileşimleri
   — sayısal doğrulama: maks hata ~5e-8 (Gram ölçeği 0.006). Cetvel
   teorisinin cebiri KESİN: bütün çiftlenim Ĝ'de yaşıyor.
2. AYRIŞTIRMA (T2): atlanmış-değişken transferi kesin cebirle kapandı
   (çıplak−tam = bias, iki pencerede +0.429/+0.450 birebir). Katkı
   dağılımı: TEK taşıyıcı YOK — 11 atlanan çizginin hepsi pozitif katkı,
   q ile yavaş azalan (+0.06'dan +0.02'ye): transfer KOLEKTİF.
3. TAŞIYICI KİMLİĞİ (T3): en büyük çiftlerin frekanslarında |Ĝ| =
   0.002-0.009 — p^k çizgilerinde DEĞİL (uzaklıklar 0.03-0.44), ama
   yerel karanlık tabanın (5e-4) 4-15 KATI, ve Σω büyüdükçe tarağa
   (ω=L) doğru YÜKSELİYOR. Kimlik: TERMAL DİFÜZ SAÇILMA (TDS) —
   hiperuniform karanlıkla Bragg tarağı arasındaki termal omuz.
   Kaba 1B kristal kestirimi √((1−e^{−ω²σ²})/2n) ≈ 0.0035 @ω=6.5;
   ölçülen 0.005-0.009 — aynı mertebe.

BÜYÜK RESİM: üç günün draması tek cümleye indi — "regresyonu büken şey,
sıfır kristalinin termal difüz saçılmasıdır." Halka: v-alanı → Ĝ
(benekler+tarak+TDS) → Gram → giydirme → düzeltilmiş fizik. Not 4'ün
merkezi bölümü hazır. Açık: TDS'nin sıcaklık ölçeklemesi (pencereler
arası test: kaba işaret karışık — n ve √p farkları ayıklanmalı);
kontrollü-rung kontrol-aracılı pay; tam difüz-model fiti (DW²·tarak +
(1−DW²)·difüz + benekler = tam S(ω) modeli). Script: 87_kapali_devre.py.

## 88: TAM S(ω) MODELİ — BENEK YASASI ÇÖZÜLDÜ, RAMPA GÖRÜLDÜ (21 Ağustos)

M1 — BENEK-SÖNÜMÜ YASASI ÇÖZÜLDÜ: gap regresyonu dalganın FARKINI ölçer
(2sin(πτ)), ortalar ORTALAMASINI taşır (cos(πτ)) → saf geometri düzeltmesi
f(τ) = πτ·cot(πτ). Ölçülen oran/f = 1.083 ± 0.007 — ON BİR ASALDA DÜMDÜZ
(86'nın 1.07→0.70 düşüşü tamamen geometriymiş). f(1/2)=0: beneklerin
katta sönmesi = Nyquist. Kalan tek sabit k₀=1.083 (v-konvansiyonu ya da
çizgi-üstü korelasyon adayı) — yeni açık sabit.

M2 — TARAK YARISI ÇÖZÜLDÜ: k=1'de Gauss 0.32 vs ölçüm 0.256; ölçülen
kurtosis (−0.75) ile Edgeworth 0.27 — 0.885 açığının ÇOĞU alt-Gauss
düzeltmesi, ~%6 artık (yüksek kümülantlar). k=2 tarağı fiilen sönmüş
(ölçüm ~gürültü; Edgeworth k=2'de çöker — beklenen).

M3 — BERRY AYRIŞIMI DOĞRUDAN GÖRÜLDÜ (ince-ızgara F(α)=n|Ĝ|², 36k nokta):
İlk çizgiden önce (α<0.65/L·L) F̄=0.002 — GUE rampasının 18 KAT ALTINDA
(sonlu-yükseklik hiperuniform karanlığı); log2-log3 binine gelince
çizgiler DAHİL F̄=0.102 vs rampa 0.089 — RAMPA BENEKLERDEN GERİ GELİYOR
(Berry 1988'in resmi, ortalar üzerinde ilk doğrudan ölçüm); çizgi-arası
süreklilik 10-20 kat altta kalıyor. Yüksek α'da F̄ rampanın gerisine
düşüyor — v-DOYUMUYLA AYNI BÖLGE: yeni teori köprüsü adayı —
MONTGOMERY RAMPASI (teorem!) + kırınım kimliği ⟹ v(τ) üzerinde kesin
toplam kuralı ⟹ v-yasası normalizasyonunun (ve belki k₀'ın) TÜRETİM YOLU.

Model artık tam: S(ω) = benekler[(v·q^{-1/2}/2)·πτcot(πτ)·k₀] +
tarak[char-fonksiyon, alt-Gauss] + karanlık alan[sonlu-L rampa inşası].
Not 4 omurgası üç bölümüyle hazır. Script: 88_S_modeli.py + png.

## 89: MONTGOMERY KÖPRÜSÜ — MUTLAK BENEK YASASI DOĞRULANDI, v-KURALI FİZİĞE İNCELDİ (21 Ağustos)

TÜRETİM: S(t) çizgisi Λ(q)/(π√q log q) ⟹ u-çizgisi U_q = 2Λ/(L√q log q);
ortalar ortalama (cos πτ), gap'ler fark alır ⟹
  BENEK: |Ĝ(log q)| = Λ(q)/(L√q)·cos(πτ)·DW   [kanal ölçümsüz, mutlak]
  GAP:   v(τ) = (2/π)·sin(πτ)  [naif fark-geometrisiyle]
  → onset 2τ TÜREDİ (ölçülen 2.014); doyum ölçeği 2/π = 0.637;
  → kuvvetler için v(p^k) = (2/π)sin(πτ)/k (yeni öngörü, test edilecek);
  → Pythagoras: benek ∝ cos²(πτ), gap ∝ sin²(πτ) — rampa bölüşümü.

T2 — MUTLAK BENEK YASASI: 1-4% İÇİNDE DOĞRULANDI. p=2: 0.999, p=3: 1.002
(!!), p=31: 1.040; ortalama 1.022±0.014. Explicit formül, sıfır örgüsünün
kırınım desenini HİÇBİR kanal ölçümü girmeden öngörüyor — programın en
temiz ilk-ilke temas noktası. (+%2-4 hafif eğim: Gauss-DW'nin çizgiler
için fazla sönümü — alt-Gauss düzeltme adayı, açık.)

T3 — ENJEKSİYON KALİBRASYONU (gerçek diziye bilinen dalga): cos-tarafı
KUSURSUZ (a_G = 0.943/0.674/0.281 vs cos(πτ) = 0.951/0.707/0.309) —
86'nın benek geometrisi bağımsız doğrulandı. SÜRPRİZ: sin-tarafında
boru hattı naif U·ω'yu neredeyse aynen döndürüyor (a_v = 1.00/0.98/0.92),
sabit-gap sin-formülünden ÇOK daha az zayıflatıyor.

T1 + T3 BİRLEŞİMİ — YENİ FİZİK NESNESİ: benekler U'nun EF değerini
kanıtladığına ve boru hattı sadık olduğuna göre, ölçülen v'nin
2τ·a_v'den kalan açığı (0.92→0.75, τ ile büyüyen) GERÇEK bir bastırma:
FARK-KANALI PERDELEMESİ D(τ) — sıfırlar asal dalgaya katı yerdeğiştirme
olarak değil KISMEN KORELE hareketle yanıt veriyor (komşular birlikte
kayınca gap'ler daha az nefes alır; ortalama kanalı etkilenmez —
benekler tam, v eksik: tutarlı!). D(τ) ≈ sinc(πτ)·κ, κ ≈ 0.86-0.92:
sinc = "bir sıfırın yanıtı dalgayı ~bir örgü aralığı üzerinden ortalar"
(fiziksel çekirdek) → v(τ) = (2/π)sin(πτ)·κ formu fizikle geri geliyor;
κ<1 artığı açık. T4 rampa: orta binlerde uyum (0.076 vs 0.081; 0.061 vs
0.050), düşük binlerde Λ²-topaklanması — nitel ✓, nicel pürüzlü.

Not 4 hiyerarşisi netleşti: (1) mutlak benek yasası [teorem-komşusu],
(2) tarak/sıcaklık [char-fonksiyon], (3) fark-kanalı perdelemesi D(τ)
[yeni gözlemlenebilir], (4) rampa bölüşümü [Berry-Montgomery bağı].
Script: 89_montgomery_koprusu.py + 89_kopru.png.

## 90: 1/k TESTİ — Λ-AĞIRLIĞI İKİ KANALDA KESİN DOĞRULANDI (21 Ağustos)

89'un öngörüsü v(p^k) = (2/π)sin(πτ)·D/k iki bağımsız kanalda test edildi:

T2 GAP KANALI (keskin oran testi, κ/D-bağımsız): R = k·v_q/v_asal(τ_q):
  q=4: 1.013±0.009 | 8: 1.005±0.014 | 9: 1.005±0.009 | 16: 0.953±0.021 |
  25: 0.984±0.012 | 27: 0.953±0.019 — havuz ⟨R⟩ = 0.986.
  1/k HİPOTEZİ KAZANDI (karşı hipotez R≈k=2-5 derdi; ölçülen R/k
  0.24-0.51'de ezildi). k·v(p^k) asal eğrisinin ÜSTÜNE düşüyor:
  v-yasası gerçekten Λ-ağırlıklı. Plasebo 0.013 (sinyal ~0.1-0.4).

T1 KUVVET BENEKLERİ (mutlak, regresyonsuz): 1/k formu her yerde kazanıyor
  (4: 0.992, 8: 0.964, 9: 0.975 — asallarla aynı kalite!); zayıf/yüksek-k
  çizgilerde ek açık: 16: 0.87, 25: 0.89, 27: 0.81, 32: 0.65, 49: 0.72.

YENİ TUTARLI İNCE YAPI — "KUVVET AÇIĞI" ε(q,k): kuvvet çizgileri tam
Λ-ağırlığının birkaç %-%35 altında, ÜÇ bağımsız nesnede aynı desen:
85'in koşulsuz u'ları (0.99/0.96/0.90/0.93/0.90/0.83/0.90), 90'ın
benekleri (0.99/0.96/0.87/0.89/0.81/0.65/0.72), T2'nin hafif 0.95'leri.
Asallar tam, kuvvetler eksik — q ve k ile büyüyen açık. Aday teori:
S(t) açılımının ikinci-mertebe (q^{-3/2}) düzeltmeleri. AÇIK NESNE.

Zincirin durumu: EF → u-çizgileri (Λ-ağırlık ✓ 1/k ✓) → benekler
(mutlak %1-4 ✓) → v-yasası (onset 2τ ✓, doyum 2/π ✓, perde D(τ)) →
rampa bölüşümü. Not 4 malzemesi: 86-90, beş script, tek anlatı.
Script: 90_kuvvet_v_testi.py.

## 91: KUVVET AÇIĞI LABORATUVARI — İKİ ADAY ÖLDÜ, GÜÇLÜ KISIT DOĞDU (21 Ağustos)

Türetim denemesi laboratuvarla: EF çizgilerinden (tüm p^k ≤ 720, bilinen
genlik, gürültüsüz) sentetik örgü, örtük sıfır koşulu ρ̄u = −S(t+u) tam
çözüldü, gerçek ölçüm makinesi üstünde koşuldu.

1. LAB DESENİ ÜRETTİ: kuvvetler asallardan sistematik derin (0.92→0.46,
   sıralama gerçekle uyumlu); gap-oranının ≈1 kalması bile kopyalandı
   (R_lab 1.01-1.04). Ama gürültüsüz lab açıkları abartıyor VE kuvvet
   fazlarını döndürüyor (−176°→−150° — 2. mertebenin parmak izi).
2. ANINDA ÖLDÜRME TESTİ (lab'ın öngörüsüyle): gerçek veride kuvvet-beneği
   fazları ölçüldü — SAF 180°, 0.1-0.6° içinde (32 ve 49 dahil!).
   FAZ DÖNMESİ YOK → koheran S'S mekanizması gerçek açığın kaynağı değil;
   gerçek alanın rastgele bileşeni koheran 2. mertebeyi dekohere ediyor.
3. Bessel/harmonik-örnekleme adayı da nicel elendi: J_4(0.38) 16'nın
   %13'ünü veremez; ε(49) < ε(25) sıralaması U_p-hiyerarşisine ters.
4. KALAN KISIT SETİ (yeni ölçüm): mekanizma faz-koruyucu (saf reel),
   yalnız p^k'ları seçiyor, q VE k ile büyüyor, asalları es geçiyor.
   Açık adaylar: dekohere-filtreli harmonik terimler; S(t) çizgi
   genliklerinin kendi 2. mertebe düzeltmeleri; kendi-çizgi × harmonik
   girişim analitiği. Dürüst durum: kuvvet açığı TÜRETİLEMEDİ — ama iki
   yanlış kapı kapandı ve faz-saflığı kısıtı (yeni, keskin) kayda geçti.

Script: 91_kuvvet_acigi_lab.py (+faz testi çıktısı logda).

## 92: PERDE KÖPRÜSÜ — NAİF KÖPRÜ ÖLDÜ, DAHA DERİN BİR UYUM DOĞDU (21 Ağustos, gece)

Hipotez: D(τ) perdelemesi CUE log-gazının doğrusal yanıtından türer.
Dört katlı test (a_v 6-nokta enjeksiyon kalibrasyonlu):

T1 ÖLÇÜLEN PERDE (temiz eğri): D = 0.908 → 0.863 → 0.797 → 0.726 →
   0.634 → 0.484 (τ 0.05→0.55). T2 null (analitik): additif çekirdek
   D≡1 verirdi → ölçüm çekirdek-çizgi ETKİLEŞİMİNİ kanıtlıyor.

T3 NAİF KÖPRÜ ÖLDÜ: CUE koşullu yanıtı perdelemiyor, TERS yönde
   (D_CUE = 1.00→1.18, κ→π; sağlama ⟨|ρ_m|²⟩ = min(m,N) ✓ kusursuz —
   Diaconis-Shahshahani). Zeta perdesi hem null'un hem CUE'nun altında.

T4 ASIL KEŞİF: zeta KENDİ spontane modlarına (çizgi-dışı bantlar,
   pseudo-topluluk) CUE GİBİ yanıt veriyor: D_spont = 1.075 / 1.088 /
   1.177 @ τ = 0.3 / 0.4 / 0.5 — CUE'nun 1.081 / 1.132 / 1.183'üne
   karşı (Nyquist'te binde-beş!). τ=0.2: 0.956; τ=0.1: 0.604 ŞÜPHELİ
   (pencere-içi yoğunluk sürüklenmesi düşük-ω kirliliği; detrend
   tekrarı gerek — bayraklı).

SENTEZ: perdeleme GAZIN değil SÜRÜCÜNÜN özelliği. Sıfır gazı kendi
gürültüsüne CUE gibi yanıt veriyor (yeni, keskin bir CUE-uyum ölçümü:
koşullu gap-yanıt fonksiyoneli — bilinen testlerden bağımsız); asal
korosuna ise perdeyle. D_asal/D_spont ≈ 0.63/1.09 @ τ≈0.45 → saf
aritmetik-koherans fiziği İZOLE edildi: perde, koheran deterministik
sürücülere özgü. Açık: neden koheran sürücü perdelenir (adyabatik-vs-
ani yanıt; sayım-özdeşliği gömülmesi); düşük-τ spontane detrend.
Script: 92_perde_koprusu.py + 92_perde.png.

## 93: DÜŞÜK-τ SORGUSU — ARTEFAKT DEĞİL, DONMA→CUE GEÇİŞİ (21 Ağustos, kapanış)

92'nin şüpheli τ=0.1 noktası sorgulandı: yerel katlama + kübik detrend
HİÇBİR ŞEYİ DEĞİŞTİRMEDİ (üç ondalık aynı); faz-karıştırmalı vekil taban
0.01-0.05 (estimatör temiz). Düşüklük GERÇEK — ve ince ızgara altında
yepyeni bir yapıya çözüldü:

D_spont(τ): 0.041 (τ=0.04) → 0.240 (0.06) → 0.492 (0.08) → 0.662 (0.10)
→ 0.788 (0.125) → 0.925 (0.15) → 0.964 (0.20) → 1.05/1.11/1.18 (CUE ✓).

OKUMA: ilk asal çizgisi eşiğinin (τ₂ = log2/L ≈ 0.067) altında sıfır
gazı fiilen DONMUŞ — konum modülasyonu var, gap yanıtı ~sıfır: rampa-altı
karanlık bölgenin (86'da 18× ölçülmüştü) dalgalanmaları termal değil,
katı kolektif öteleme. Çizgiler devreye girdikçe süreklilik CUE-termal
davranışa geçiyor; Nyquist'te binde-beş CUE uyumu. DONMA→TERMALLEŞME
geçişi tam asal çizgilerinin açıldığı bantta.

İKİ-DAL YANIT RESMİ (Not 4'ün taç figürü adayı): asal dalı 0.91→0.48
İNER, spontane dal ~0→1.18 ÇIKAR, kesişme τ≈0.15-0.2. Uzun dalgada
asallar gap'leri spontane modlardan çok daha etkin sıkıştırıyor
(adyabatik denge-deformasyonu); kısa dalgada perdeleniyor. Açık soru
inceldi ve ikiye bölündü: (a) asal dalının perde çekirdeği; (b) spontane
dalın fonon-katılığı — ikisinin kesişme noktası yeni bir ölçek (~0.17).

Script: 93_spontane_detrend.py. 86→93: cetvel teorisi seferi TAMAM —
sekiz script, iki teorem-teması, üç yeni nesne (D-dalları, kuvvet açığı,
kırınım deseni), dört dürüst ölüm (KK-köprüsü adayları). LİMANA DÖNÜŞ:
sıradaki oturum = Not 4 taslağı + üçleme okuması + arXiv.

## NOT 4 TASLAK YAZILDI — TERSANEDEN DENİZE (21 Ağustos, akşam)

`arxiv_warm_crystal.tex`: "The Riemann zero lattice as a warm crystal:
diffraction, screening, and the two-branch response". 86-93 seferinin
tamamı tek çatıda, dokuz bölüm:
  §2 kırınım deseni (benekler saf-180°/tarak-termometre/karanlık alan 25×)
  §3 mutlak benek yasası (parametresiz, %0.1-4; kuvvetlerde Λ ve 1/k)
     + benek geometrisi πτ·cot(πτ), k₀=1.083 = 1/D(0) kapanışı
  §4 cetvel kimliği (Gram=Ĝ, 5e-8) + giydirme=TDS + Berry-rampa inşası
  §5 gap toplam-kuralı (2τ onset, 2/π doyum — türetilmiş) + enjeksiyon
     kalibrasyonu + perdeleme D(τ) 0.91→0.48 + Λ/1k oran testleri
  §6 iki-dal yanıtı: CUE referansı (min(m,N) sağlamalı), spontane dal
     (donma→CUE, τ₂ eşiğinde çözülme, Nyquist %0.5), asal dal; kesişme
     τ≈0.17
  §7 kuvvet açığı: iki ölü mekanizma (lab + faz-saflık testi; Bessel)
     + üç-kısıt seti
  §8 sentez + 6 açık problem; repro 79-94.
4 İngilizce figür (94): fig_diffraction/geometry/bridge/twobranch_en.
Derlendi (409 KB). KÜNYE TEYİDİ GEREK (hafızadan): Montgomery 1973
(PSPM 24, 181-193), Diaconis-Shahshahani 1994 (J. Appl. Probab. 31A,
49-62) — gönderim öncesi web teyidi şart (BCHB geleneği). SIRADAKİ:
taze düşman-göz denetimi (Not 4) → kullanıcı okuması → dörtleme arXiv.

## NOT 4 DENETİMİ + DÜZELTMELER (21 Ağustos, gece — 2 KIRMIZI + 12 SARI işlendi)

Taze düşman-göz denetimi: 86-93 + AĞIR 92 (8 dk CUE) sıfırdan bağımsız
koşuldu; üç künye web-teyitli (Montgomery 1973 PSPM XXIV 181-193 ✓,
Diaconis-Shahshahani JAP 31A 49-62 ✓, Berry Nonlinearity 1 ✓); türetim
cebiri elle doğrulandı; scriptlerde olmayan iki ölçüm bağımsız üretildi.
Karar: "sayısal omurga olağanüstü sağlam; iki kırmızı bir akşamlık iş."

K1: τ≈0.17 kesişmesi 92'nin KABA ızgarasından kalmaydı (0.604@0.1
    aradeğerlemesi) — 93'ün düzelttiği veriyle τ* = 0.140 (0.13-0.15;
    çapraz-estimatör sistematiği bütçelenmedi notuyla). Dört yerde
    düzeltildi. 95-T2 kesişmeyi numaralı script yaptı: 0.140 ✓.
K2: §5'in "onset 2.014τ / doyum 0.55-0.60" doğrulaması eski-konvansiyon
    sayılarıydı — tam-taban gerçeği: onset (1.7-1.8)τ = 2τ·D(0), doyum
    ≈0.47 = (2/π)·D(kat); 2.014≈2 uyumu iki sistematiğin tesadüfi iptali
    olarak KAYDA geçti. Rijit-öngörü/perdeli-ölçüm ayrımı netleşti.
SARI'lar: 15→14 çizgi; k₀=1/D kapanışı "%1.5'e kadar" (1.101±0.003 vs
1.083±0.007) + kapanış cebiri Açık Problem 3'e; abstract'ta iki-mekanizma
özeti düzeltildi (biri fazla, biri nicel öldü); faz-saflığı belirsizlik-
dürüstlüğü (güçlü beneklerde 3-5°, zayıflarda ±24°; kill güçlülerden);
"exactly like CUE" → "%3 içinde, Nyquist'te %0.34"; tarak çift-taban
açıklaması; 978/1200 frekans; enjeksiyon %1-9 + doğrusallık (95-T3:
5× genlik → <%0.7); bant-karışımı sistematiği (donma yönsel sağlam —
karışım yukarı çeker); "conditional (regression) response" adlandırması;
TDS 4-18×; 630-950 frekans; TDS cümlesi "consistent with" tonuna.

95_denetim_ekleri.py: faz kill-testi (KOD olarak; 180.00-180.65°,
ist. belirsizlikler), dal kesişmesi (0.140), enjeksiyon doğrusallığı.
Yeniden derlendi (422 KB). Repro 79-95. DÖRTLEME DURUMU: dördü de
yazılmış + denetimli + düzeltilmiş. SIRADAKİ: kullanıcı okuması (4 not)
→ arXiv mekaniği (endorsement erken başlatılabilir).

## 96: DIRICHLET BETA'NIN KIRINIMI — EVRENSELLİK + KARAKTER OKUMA (21 Ağustos, akşam)

Not 4 açık problem 6 aynı gün saldırıya uğradı: programın İLK zeta-dışı
veri seti üretildi — L(s,χ₄) (Dirichlet beta), 4.962 sıfır, kendi
tamamlanmış-fonksiyon makinemizle (Hurwitz-zeta + faz; iç sağlama
maks|Im Z| = 1.7e-11 ✓). Dört öngörü ÖLÇÜMDEN ÖNCE script başına
mühürlendi; DÖRDÜ DE İSABET:

P1 ✓ 2/4/8 BENEKLERİ SÖNDÜ (χ₄(2)=0): 0.0011/0.0006/0.0003 —
   karanlık-alan seviyesinde, fazlar rastgele. Zeta'da 2-beneği 0.043.
P2 ✓✓ FAZLAR KARAKTERİ OKUDU: χ=−1 asalları (3,7,11,19,23) 0.5-2.5°
   (kilit 0°); χ=+1 (5,13,17) ve q=9 → 178-181° (kilit 180°).
   "Kırınım deseniyle aritmetik spektroskopi" İLK KEZ gösterildi.
P3 ✓ mutlak yasa iletken-L (log 4t/2π) + tarak-termometreli DW ile:
   asallar 0.99/0.96/0.96/0.94/0.91/0.90/0.92/0.81 (zeta ailesiyle aynı
   biçim); q=9: 0.76 — KUVVET AÇIĞI BETA'DA DA (evrensel ince yapı!).
P4 ✓ karanlık alan 1/√n'in 12 kat altı (hiperuniformluk L'de de);
   tarak 0.60 → σ_u = 0.161 — zeta'nın log-ısınma eğrisinin soğuk ucu
   (L=7.2), SICAKLIK YASASI L-FONKSİYONLAR ARASI uzuyor gibi.

Dürüst notlar: tarama ~27 yakın-çift kaçırdı (%0.5; benekler duyarsız;
indeksli akış bozulduğundan sıcaklık tarak-termometreyle — indekssiz);
tarak Gauss-varsayımlı (0.885 düzeltmesi benek-DW'ye <%1 etkir); benek
belirsizliği ~0.001-0.0015. NOT 5 TOHUMU: "Arithmetic spectroscopy:
reading Dirichlet characters off the diffraction pattern of L-function
zeros" — evrensellik-taşınan sabitler (cos-geometri, DW-sıcaklık,
karanlık alan, kuvvet açığı) vs karakter-taşıyan yapılar (benek
varlığı/fazı). Script: 96_beta_kirinim.py + npz önbellek + figür.

## 97 (KISIM 1): χ₅ KADRANI OKUNDU — KOMPLEKS KARAKTER FAZLARDAN SÖKÜLDÜ (20 Ağustos, sabah)

4601 χ₅ sıfırı (10 saatlik mpmath seferi; Hurwitz büyük-t'de tahminden
~50× yavaş — ders). ÖN-MÜHÜRLÜ ÖNGÖRÜLER:

P2 KADRAN ✓✓✓ — dört konum, derece hassasiyeti:
  q≡2 (χ=i): 269.3/270.0/271.5° | q≡3 (χ=−i): 89.9/89.4/89.7/89.5°
  q≡4 (−1): 1.6/359.6/0.3° | q≡1 (+1): 180.0°.
  SANAL BİRİM i KIRINIM FAZINDAN OKUNDU — aritmetik spektroskopinin
  kompleks hali ilk kez. P1 ✓ mezarlık taşındı: 5, 25 ölü (plasebo
  dibi), 2 dirildi. P3 ✓ genlikler mutlak yasada (asallar 0.94-0.99);
  kuvvet açığı ÜÇÜNCÜ örgüde de aynı bant (8: 0.77, 9: 0.81).

P5 SÜRPRİZİ → İNCELMİŞ YASA: χ₅ σ_u=0.220 vs β 0.161 (aynı L!) — naif
tek-eğri YANLIŞ; ama fark tam 2-ailesinin varyansı: Δσ² = 0.0225 vs
hesap Σ U²(2,4,8)/2 = 0.0216 (%4!). "HER ADA KENDİ KOROSU KADAR
ISINIR." ÖN-MÜHÜR (χ₃ fırındayken): χ₃ tarak-termometresi σ_u ≈
0.20-0.21 vermeli (3-korosu susmuş). ζ-adalar tam kıyası: termometre-
konvansiyonu inceliğiyle açık (χ₃ inince aynı ayakla).

## 97 (KISIM 2): χ₃ İNDİ — KADRAN YİNE OKUNDU, KORO-YASASI ÜÇ NOKTADA (20 Ağustos)

4755 χ₃ sıfırı (L_eff=6.92). P4 KADRAN ✓: 3 ve 9 ölü (0.0008/0.0002,
plasebo dibi); χ=−1 sınıfı (2,5,8,11) → 0.4/0.3/359.5/0.1°; χ=+1
(4,7,13) → 180.2/179.3/180.3°. Genlikler: asallar 0.97-1.00 (q=2:
1.003 — mutlak yasa DÖRDÜNCÜ örgüde de tam); kuvvet açığı yine aynı
bant (8: 0.74, 4: 0.93) — dört örgüde evrensel ince yapı.

ÖN-MÜHÜR KARNESİ: öngörü σ_u ≈ 0.20-0.21 (L'yi 7.6 sanarak); ölçüm
0.196 @ L=6.92 — bandın dibinde, İSABET sayılır (kaba öngörü hassasiyeti
içinde). KORO-YASASININ ÜÇ-ADA TESTİ: her adaya susturulmuş ailesinin
varyansı geri eklenince "tam" sıcaklıklar: χ₃ 0.0535, β 0.0481,
χ₅ 0.0562 (L 6.9-7.3) — ham yayılım ±%32'den ±%8'e ÇÖKÜYOR.
"Her ada kendi korosu kadar ısınır" ilk denetimden geçti; kalan ±%8
(pencere-içeriği, alt-Gauss, kaçırılan-çift oranı) açık. ζ-adalar tam
kıyası hâlâ konvansiyon-inceliği bekliyor (aynı-ayak termometre).

NOT 5 İSKELETİ TAMAM: dört örgü (ζ, β, χ₃, χ₅) — mutlak benek yasası
4/4; kadran iki karakterde (biri KOMPLEKS) derece-hassas; mezarlık üç
kez taşındı; kuvvet açığı 4/4; koro-sıcaklık yasası. Başlık adayı:
"Arithmetic spectroscopy: reading Dirichlet characters off the
diffraction patterns of L-function zero lattices". Ders: Hurwitz
büyük-t'de yavaş — Not 5 kampanyası için vektörize L-motoru gerek.

## 99: HEXAGON KAMPANYASI — ALTIGEN OKUNDU, KUVVET-AÇIĞI YASASI DOĞDU (20 Ağustos)

Motorla dört ada, ~291.000 sıfır, TOPLAM 19 SANİYE (χ₃ 79.9k, β 74.0k,
χ₅ 72.4k, χ₇ 64.4k; önbellekler 99_*_zeros.npz).

H1 ✓✓✓ ALTIGEN KADRAN (mod 7, 6. birim kökleri) YARIM-DERECE İÇİNDE:
  180°: 29(179.8), 8(180.5) | 240°: 3(240.4), 17(240.3) |
  300°: 2(300.6), 23(299.9), 9(299.8) | 0°: 13(359.5), 27(357.2) |
  60°: 4(60.0), 11(60.2), 25(59.8) | 120°: 5(120.4), 19(119.7).
  Ölü: 7(0.0003), 49(0.0002). e^{iπ/3} bir nokta kümesinden ölçüldü.
  χ₅ kadranı da 15× istatistikle keskinleşti (270.1/90.1/0.3/180.0).

H2 ✓ MUTLAK YASA BİNDE-DÜZEYİNE İNDİ: asal oranları dört adada
  0.96-1.00; χ₅'te 17: 1.001, 13: 1.002, 23: 1.001, 11: 0.998 —
  parametresiz yasa, binde-birkaç isabetle.

H4 ✓ → YENİ NİCEL YASA — KUVVET-AÇIĞI ÇÖKMESİ: (1−ε)/k = f(τ) tek
  eğriye biniyor: f(0.146)=0.007 [4@χ₃], f(0.226)=0.0285 [9:0.0285 ve
  8/3:0.029 — k=2 ve k=3 AYNI f!], f(0.333)=0.105-0.112 [25/2 ve 27/3
  yine aynı!]. Eğim log-log ~3.2-3.4: f ≈ c·τ^{~3.3}. KUVVET AÇIĞI
  ARTIK YASALI: 1−ε ≈ k·f(τ), f ~ τ^3.3 — üçüncü-mertebe koku
  (teori adayı: S(t) genliklerinin τ³-düzeltmesi?). 91'in üç kısıtı
  (faz-koruyucu ✓ burada da, p^k-seçici, q&k-büyüyen) yasayla uyumlu.

H3 ✓/△ ADA-İÇİ ISINMA: dört adada da log-ısınma NET (β: 0.142→0.194,
  χ₃: 0.181→0.226, χ₅: 0.207→0.248, χ₇: 0.216→0.255; L~6→10) —
  ısınma evrensel. Koro-sıralaması ✓ (β en soğuk < χ₃ < χ₅ ≲ χ₇).
  Nicel koro-düzeltmesi: küçük-L üçlüsünde ±%8 idi; büyük-L dörtlüsünde
  ±%15 ve β düzeltme-sonrası bile sistematik soğuk — açık inceltme
  (pencere-içeriği/alt-Gauss/kaçırılan-çift; iletken-4 imprimitiflik?).

NOT 5 MALZEMESİ TAŞTI: 5 örgü (ζ dahil), kadran 3 karakterde (kare +
altıgen), mutlak yasa binde-düzeyi, kuvvet-açığı YASASI, ada-içi ısınma
+ koro. Scriptler: 98 motor + 99 kampanya.

## 100: KUVVET-AÇIĞI YASASI MEKANİK OLARAK TÜRETİLDİ (20 Ağustos, gece)

f ~ τ^3.3'ün türetim oturumu, iki perdelik:

PERDE 1 — kaba düğme: lab tabanı %5 yanlıydı (asal kontrol 1.05),
eğim fiti savruldu (4.9/2.4/3.5). Teşhis: mutlak-taban sistematiği
açıklarla aynı mertebede. Çare gerçek-veri metodolojisinin aynısı:
ORAN-TESTİ (lab'ın kendi asal eğrisine bölme, taban iptali) + çok-tohum.

PERDE 2 — SONUÇ: saf-koheran lab (gürültüsüz!), taban-iptalli okumayla
gerçek yasayı VURDU: eğim 3.26 (gerçek 3.3), f(0.23) = 0.0306 (gerçek
0.030) — sıfır serbest parametre. TÜRETİM KAPANDI (mekanik düzeyde):
(1−ε)/k = f(τ) ≈ 4·τ^{3.3} yasası, örtük sıfır-koşulu + EF aile-
çizgileri + orta-nokta örneklemesinin oran-testi görüntüsüdür.
DÜZELTME: 91'in "koheran mekanizma τ^1.8 verir" hükmü mutlak-taban
artefaktıydı; koheran mekanizma DOĞRU yasayı veriyor.
DEKOHERANS HİPOTEZİ ÖLDÜ: iid gürültü eğimi bozuyor (2.6-2.3) —
gerçek rastgelelik iid değil (hiperuniformlukla tutarlı).

AÇIK: (a) faz-saflığı — lab ~4-7° döndürür, gerçek ≤1°: döndüren
bileşeni gerçekte bastıran şey (korele gürültü / S'S-örnekleme iptali)
türetilmedi; (b) kapalı form — validen modelin asimptotiği (τ³-ailesi
cos/log düzeltmeli), kâğıt-kalem ödevi. Script: 100_dekoherans_dugmesi.py.
Not 5'e girecek hali: "the power-deficit law is reproduced, slope and
amplitude, by the coherent family mechanism read through the
ratio-test — with its phase purity still stricter in nature than in
the model."

## 101: DONMUŞ KOYLAR — ÇÖZÜLME SINIRI KARAKTERLE PROGRAMLANIYOR (20 Ağustos)

Tatil modu kararı: gönderim/Not-5 yazımı kaptan dönene dek ertelendi;
ölçüm sürüyor. Kuyu seçimi: 93'ün donma→CUE geçişinin takımada testi —
"gazı çözen aritmetiktir" iddiasının KONTROLLÜ deneyi, çünkü doğa
mükemmel bir kontrol hediye etmiş: β'nın 2-çizgisi ÖLÜ (χ₄(2)=0).
Ön-mühür: P1 (χ₃/χ₅/χ₇ ζ gibi log2/L'de çözülür), P2 (β log3/L'e
GECİKİR, ~×1.585), P3 (τ≥0.3'te hepsi tek CUE), P4 (ζ aynı-ayak tarak).

PERDE 1 — ARTEFAKT VE KAPI (101, 101b): ham ölçüm dört adada "donma
yok" dedi (τ=0.04'te D≈0.97). Kusur kapısı bunu ÇÜRÜTTÜ: donma
estimatörü kinematik kilidin (yapısız süreç → D=1/cos(κ/2); Poisson
kontrolü 1.016 ✓) İHLALİNİ ölçer; temiz ζ'ya %0.1 sıfır-silme
enjeksiyonu D'yi 0.09→0.55'e fırlatır. Benek/faz ölçümleri kusura
bağışık; donma ölçümü DEĞİL.

PERDE 2 — MOTORUN KÖR NOKTASI (101c, 101e, 101f): ince ızgara (4×)
sadece +54-72 buldu ama sayım sertifikası (d_i = i − Δθ/π sürüklenmesi)
ada başına ~285-620 kayıp gösterdi. Teşhis: KALDIRILMIŞ DİPLER —
motorun RS-düzeltmesiz O((qt)^{-1/4}) hatası dar çiftin |Z| çukurunu
sıfır üstüne kaldırır; işaret-taraması HİÇBİR ızgarada göremez.
Hibrit kurtarma: dip adayı (işaretsiz |Z| minimumu) + 3-5 mpmath
(Hurwitz, gerçek-Z) + parabol kökleri → ada başına ~290-790 sıfır
kurtarıldı. Kalan pürüzler (tekrar-ekleme kopyaları, χ₇'de ~10 inatçı)
için nihai disiplin: DÜZLÜK SEGMENTASYONU — sayım-sürüklenmesi düz
(|Δmedyan|<0.5) parçalar sertifikalıdır, her basamakta kes, çevresini
at. Kusur analize giremez. ζ'ya da aynı disiplin (0 kesim çıktı ✓).
YENİ STANDART: motorlu kampanyalarda hiperuniformluk-hassas ölçümler
(donma, karanlık alan, sıkıştırılabilirlik) sayım-sertifikası +
düzlük-segmentasyonu İSTER; benek/faz ölçümleri istemez (98-S2).

PERDE 3 — SONUÇ (101d, figür 101_donmus_koylar.png): P1 ✓ P2 ✓✓✓ P3 ✓
  τ=0.068 (log2/L): ζ 0.34, χ₃ 0.39, χ₅ 0.35, χ₇ 0.47 — çözülmede;
  β 0.040 (vekil taban 0.028) — DONUK. β 0.085→0.113'te dik çözülür
  (0.116 → 0.460 → 0.705 tam log3/L'de), 0.14'te ortak eğriye biner.
  τ=0.3: beşi 1.01-1.06 (tek CUE). Yarı-çözülme kayması ×1.4-1.5
  (öngörü ×1.585). HÜKÜM: DONMA SINIRI ADANIN İLK SAĞ KALAN
  ÇİZGİSİDİR. Ölü kırınım beneği dinamikte de ölü: benek yoksa çözülme
  yok. Karakter seçimiyle donma sınırı KAYDIRILABİLİYOR — aritmetik
  gazı çözer iddiası artık müdahaleli-deney statüsünde.
P4 ✓: ζ aynı-ayak tarak-termometresi: σ_u = 0.2287 (L=6.99) → 0.2833
  (L=12.45), 10 pencere — H3 nicel tablosunun ζ satırları hazır.

Scriptler: 101 (ham+kapı notu), 101b (kusur kapısı), 101c (sayım-
denetimli yeniden türetim), 101e (hibrit dip-kurtarma), 101f (temizlik),
101d (segmentli ölçüm), 101g (figür). Not 5 malzemesi: "the thaw
boundary is programmable: kill a character line and the freeze extends
to the next surviving one."

## 101h-j: OPUS TAYFASI İLK SEFERİ — İKİ HİPOTEZ ÖLDÜ, KİLİT KESKİNLEŞTİ (20 Ağustos)

Bütçe rejimi: Fable %91'de → mekanik ölçümler tarifli+ön-mühürlü olarak
Opus alt-ajanına verildi (ilk deneme; ajan kendi kapılarını ekledi ve
kaptanın iki hipotezini öldürdü — düzen İŞLİYOR). Denetim+commit Fable.

Kuram sohbetinin üç-adımlı iskeleti (mod envanteri / moiré / TDS-omzu)
teste girdi; SKOR: 1 güçlendi, 2 öldü:

1) 101h (dar bant): çözülme yamacı bandla dikleşiyor (±0.02L→±0.005L:
   1.8-4.3×); "yumuşak öncü" büyük ölçüde BAND KAÇAĞIYMIŞ (band ilk
   çizgiyi kapsayınca sızıyor; kaçak kapısı eklendi). β'nın 101d'deki
   0.095/0.460 noktası da kaçaktı: kaçaksız β 0.156'da sürünüp 0.847'ye
   TAM τ_ilk=0.1143'te sıçrar → ÇİZGİ KİLİDİ DAHA NET. x=τ/τ_ilk
   çökmesi eşik civarı %38 iyileşir. AÇIK: χ₇ τ=0.04 anomalisi (band
   daralınca yükseliyor, kaçak yok; segmentasyon şüphesi, test edilmedi).

2) 101i (vuruş/moiré): HİPOTEZ RET — donuk bölge artığı vuruş moiré'si
   DEĞİL. Çözünürlük kapısı: ω=log2'de ζ keskin tepe (D 0.35→1.005→
   0.43; kontrast 1.93×), β aynı ω'da DÜZ (0.97×) — ölü-çizgi hükmü
   ilk kez TEK FREKANSTA doğrudan. Vuruşlarda %1 tümsek bile yok;
   D(ω) çizgiler arasında pürüzsüz vadi (0.36→0.19→0.23, log3'e
   yaklaşırken yükseliş). Artık D'nin kaynağı: çizgi KOMŞULUĞU.

3) 101j (TDS örtüşmesi): HİPOTEZ RET + TEHLİKELİ BULGU — Hann taperli
   I(ω) log2 çevresinde 15 kademe Δω^{-6} (pencere yan-lobu) iner,
   kusursuz-örgüden ayırt edilemez → ÇİZGİ ÇEVRESİNDE FİZİKSEL DİFÜZ
   BİLEŞEN YOK (açık formülün atomik tayfıyla tutarlı: log2 altı ve
   log2-log3 arası atom YOK, boşluk). Tapersiz (dikdörtgen) çekirdek
   1/Δω² verir = TAM TDS'NİN BEKLEDİĞİ YASA → sahte doğrulama tuzağı.
   87'NİN TDS HÜKMÜ ŞÜPHEDE (tapersiz ölçülmüştü; "karanlık taban"ı
   sızıntı düzeyinde). 86'nın karanlık-alan sayıları da alet-sınırlı
   ÜST SINIR olur (hiperuniformluk yönü sağlam, "25×" sayısı değil).
   Taperli-D kapısı 0/0'a çöktü (E) → "eşik-altı D artığı sızıntıdır"
   iddiası KANITSIZ kaldı (yalnız I(ω) düzeyinde kesin).

YENİ RESİM ADAYI (bir sonraki derin oturuma): "her şey çizgidir" —
tayf saf atomik (log q); karanlık bölgeler gerçek BOŞLUK; çözülme =
çizgi yoğunlaşması; 92'nin iki-dal ayrımı (asal 0.48 vs spontane 1.18)
belki atom-üstü vs atom-kanadı okumasının pencere optiği. Sıradaki
mekanik iş (102, Opus): 86-T4/87-T3/88-M3'ün taperli yeniden denetimi +
kanat-öngörüsü testi (Ĝ_win(ω) = Σ A_q·K(ω−log q) parametresiz fit) +
Not 4 etki listesi. Commit: 87cd917.

## 102: KARANLIK ALANIN VE "TDS"NİN SAHİBİ BULUNDU (20 Ağustos, Opus tayfası 2. sefer)

86/87/88'in taperli yeniden denetimi + kanat-öngörüsü + Not 4 etki
listesi. 101j'nin şüphesi kısmen doğru, kısmen yanlış çıktı; asıl bulgu
üçüncü bir şey: MERKEZİ KİMLİK Ĝ_orta(ω) = (1/n)Σ e^{iωz}·e^{iωg/2} —
orta noktalar örgüyü DALGALANAN YARIM-GAP kaymasıyla örnekler; bu
kaymanın dalgalanması difüz saçılmadır. Aynı ω'da SIFIRLAR 8-10 kademe
daha karanlık: sıfır tayfı ÇÖZÜLMÜŞ ÇİZGİ ORMANI + GERÇEK BOŞLUK
("her şey çizgidir" resmi ölçümle doğrulandı); "karanlık alan" ve
"TDS" ORTA-NOKTA IZGARASININ malı. Örnekleme-fazı eğrisi (yeni
gözlemlenebilir): c=0→karanlık, c=1(orta)→tepe, c=2→yine karanlık.

102a (86-T4): 86 tabanı ω<1'de sızıntı düzeyindeydi (kusursuz örgüyle
1.20×), tüm bantta değil (28×). Taperli tabanda orta-noktalarda GERÇEK
bileşen var; sıfırlarda taban 7.3e-16 = VERİ-SINIRLI üst sınır.
Dürüst hiperuniformluk: bağımsız-gap vekilinin 40-1150× altı; rampaya
göre ortalar ≥78-1.3e4×, sıfırlar ≥5.6e11-7.5e13×.

102b (87-T3, kalp): KANAT HİPOTEZİ RET (1. mertebe kanatlar ölçümün
1/2300'ü) AMA içerik GERÇEK (taperi geçer ×0.62-0.67, sızıntının
10⁴-10²⁰ katı) ve TAMAMEN ARİTMETİK: gürültüsüz açık-formül lab örgüsü
düzeyi medyanda %6 içinde veriyor. KİMLİK DEĞİŞTİ: sıfır kristalinin
TDS'si değil, ASAL ÇİZGİLERİNİN ORTA-NOKTA ÖRNEKLEMESİYLE ÜRETİLEN
DOĞRUSAL-OLMAYAN YAN BANTLARI. 87'nin giydirme mekanizması AYAKTA
(taban gerçek, sızıntı değil — 101j'nin o şüphesi ω≈4-7'de yanlıştı);
87'nin "çizgide değil" mesafeleri tam p^k listesiyle düzeltildi
(çözünürlüğün 7-199 katı).

102c (88-M3): "18× karanlık" %100 ALETMİŞ (aritmetiksiz RvM kontrol
aynı 17.2×'i veriyor!); gerçek sınırlar ≥6.5e3× (orta) / ≥2.1e12×
(sıfır). RAMPA DÖNÜŞÜ GERÇEK VE GÜÇLENDİ: log2-log3 bini taper-
değişmez 0.1013, %99.9 çizgi-üstü, MUTLAK YASA %2-4 içinde parametresiz
(1.016/1.038). "Çizgi-arası 10-20×" sinc kuyruklarıyla ~90× şişikmiş
(gerçek 924×/2.5e5×). BONUS: sıfır ızgarası Montgomery rampasını
α≳0.2'de ±%25 izler; 88/89'un "yüksek-α açığı" fizik değil orta-nokta
cos(πτ)+DW sönümü.

NOT 4 ETKİ LİSTESİ → NOT4_ETKI_LISTESI_102.md (11 değişiklik + 3 kazanç
+ ayakta kalanlar + dürüst kayıtlar). Çoğu iddia ayakta; Berry bölümü
GÜÇLENDİ; TDS kimliği ve karanlık sayıları değişecek. Tex'e dokunulmadı
(revizyon kaptanın dönüşüne). Scriptler: 102a/b/c (+png/npz).
Dürüst kayıtlar: sıfır tabanı veri-sınırlı; lab uyumu düzey uyumu;
I üstel (medyan/ort ln2); DW konvansiyonu ~1.3×; χ₇ 0.04 anomalisi
hâlâ açık. SENTEZ CÜMLESİ: "regresyonları büken şey aritmetiktir ve
örnekleme ızgarasının difüz saçılmasında yaşar."

## 87b: 87'NİN ÜÇ TESTİNİN TAPERLİ KAPANIŞ DENETİMİ (20 Ağustos, tayfa 3. sefer)

Kaptanın görevi (101j şüphesi üzerine) 87'yi yeniden denetlemekti;
102b'nin çoğunu kapattığı görüldü, 87b kalan boşlukları mühürledi:
T1 hiç yeniden koşulmamıştı, T2 kapanış artığı ölçülmemişti,
frekans-eşlenik taperli taban tablosu yoktu.

87b (87b_taperli_denetim.py): T1 Gram=Ĝ kimliği tapersiz VE
Hann-ağırlıklı biçimde ~5e-8 ile geçiyor; T2 kapanış artığı ~1e-14,
B̂ +0.429/+0.450 birebir — İKİSİ DE TABANDAN BAĞIMSIZ, AYAKTA.
T3: taşıyıcılar taperi geçiyor (|Ĝ|² oranı 0.67/0.62, 102b ✓);
taperli çizgi-dışlamalı taban 4e-4–1.6e-3, M-S sızıntısının
~3e9 katı (101j'nin "taban alet" şüphesi ω≈4–7'de kesin yanlış);
AMA frekans-eşlenik oran 3–7× (medyan 4.0/5.0) — 87'nin "15×" ucu
Σω taşıyıcılarını Δω tabanına bölme eseriymiş. Sıfır ızgarası aynı
frekanslarda 9.3e3–1.6e4× (|Ĝ|) karanlık → kimlik düzeltmesi kesinleşti:
87-T3 "sıfır kristalinin TDS'si" → "orta-nokta ızgarasının asal-çizgi
yan bant difüz alanı, tabanın 3–7× üstünde". Script: 87b (+npz).

## 107: YEDİ ADA + ζ, İKİ YENİ KAPIYLA DENETLENDİ — LOG2/LOG3 KERVANLARI
       (25 Ağustos)

105'in bıraktığı iş ("101f verilerinin kısa-çukur kapısıyla yeniden
denetimi") yapıldı ve iki yeni kapı (105e KISA-ÇUKUR + 104e SIÇRAMA)
TEK BORU HATTINDA sekiz veri kümesine uygulandı. Ölçüm 101d estimatörü,
τ ızgarası 101d'ninki, band ±0.01L, band-kaçağı kapılı, havuz DAİMA
pencere-ayrıştırmalı (104'ün pratik kuralı 4).

POZİTİF KONTROL ✓✓ (107a/107c): kapılar 104d'nin ve 105e'nin ELLE
bulduğu kusurların TAMAMINI otomatik buldu ve analiz dışında tuttu —
χ₃ t∈[3872.54,3892.79], χ₇ t∈[3081.89,3109.55] (104d'nin kopya çiftleri),
χ₇ t∈[3664.60,3680.83] (104d'nin "i₀≥1000, t≥3680" yumuşak kuyruğu, bu
kez K3 otomatik buldu), chi8e t∈[19353.5,19359.3] (105e'nin çukuru).
ADLI KUSUR DENETİMİ 8/8 ✓.

ÜÇ YENİ KUSUR (hiç raporlanmamıştı), üçü de KISA-ÇUKUR kapısından:
  χ₅  t ∈ [5810.70, 5828.48]  (23 sıfır) ← χ₅ "TERTEMİZ" SANILIYORDU
      (104e'nin sıçrama kapısı ona 0 işaret vermişti). "χ₅ temiz"
      hükmü artık geçersiz; etkisi küçük (D'ler ≤0.006 oynadı).
  χ₇  t ∈ [767.62, 831.45] ve [1903.73, 1915.34] — ikisi de χ₇'nin
      n=2578 < 3000 olduğu için zaten analize girmeyen en alt
      log-penceresinde (sertifika kaydına girer, ölçümü etkilemez).
KOPYA ENVANTERİ genişledi: χ₇'de u<0.005 olan ALTI çift var (104d ikisini
bulmuştu), chi8e'de bir. YEDİSİ DE nihai pencerelerin dışında — bu yüzden
K3'ün u-ölçütü hiç tetiklenmedi; SON EMNİYET KİLİDİ olarak duruyor.
TERTEMİZ KÜMELER: ζ (maskelenen 0), chi8o (hiçbir kapı), β (yalnız K1).

χ₃'ÜN ORTA-τ ÇUKURU — AÇIK SORU OLARAK KALDI (ve teşhis DEĞİŞTİ):
  τ         0.095  0.105  0.113  0.125  0.140
  ÖNCE      0.590  0.565  0.512  0.320  0.766   (yalnız 101d kapısı)
  SONRA     0.426  0.199  0.161  0.134  0.531   (üç kapı)
  Yalnız 816 sıfır (%1.1) maskelendi, D(0.105) 0.565 → 0.199. Çukur
  DÜZELMEDİ, DERİNLEŞTİ. 104'ün "havuzlama girişimi" açıklaması bu
  τ'larda ÇÜRÜDÜ: pencere ayrıştırmasında ÜÇ PENCERENİN ÜÇÜ DE çukurda
  (0.282/0.134/0.137, havuz 0.199 — havuz onların ORTASINDA, ALTINDA
  değil; oran 0.71-0.77, beş kümenin genel bandı 0.59-0.92 içinde).
  Buna karşılık τ=0.068'de yıkıcı girişim GERÇEKTEN var (0.320/0.145/
  0.544, havuz 0.194, oran 0.36) — 104'ün teşhisi O NOKTADA doğru.
  ⟹ χ₃'ün orta-τ çukuru ne kusur ne havuzlama; ya fiziksel, ya da henüz
  görülmemiş bir ÜÇÜNCÜ kusur sınıfı. SIRADAKİ İŞ.

ANA HÜKÜM 5'e 3 İLE GÜÇLENDİ (107b, band ±0.01L, sekiz küme):
  EŞİKSİZ AYRIŞTIRICI — D(log2/L), yani ORTAK çizgide:
    LOG2 KERVANI  ζ .461  χ₃ .508  χ₅ .482  χ₇ .512  χ₅ᵉ .492
    LOG3 KERVANI  β .046  χ₈ᵉ .009  χ₈ᵒ .048          ← 10-50× düşük
  β artık tek ada değil ÜÇ ÜYELİ BİR SINIF. "Donma sınırı adanın ilk
  sağ kalan çizgisidir" beş çözülmüş + üç donmuş ada ile duruyor.

YENİ NİCEL SONUÇ — YARI-ÇÖZÜLME MUTLAK FREKANSTA log 2'YE OTURUYOR:
  ω_yarı = τ_yarı·⟨L⟩ ;  LOG2 kervanının dört Dirichlet adası:
  χ₃ 0.6918, χ₅ 0.7054, χ₇ 0.6792, χ₅ᵉ 0.6966 → ORT 0.6932
  ve log 2 = 0.69315  (%0.01!). LOG3 kervanı: β 0.9790, χ₈ᵉ 0.9668,
  χ₈ᵒ 0.9555 → ort 0.9671, log 3 = 1.0986'nın %12 ALTINDA.
  Kervan kayması ×1.395 (101d ×1.4-1.5, 105d ×1.376-1.410 ile aynı
  bant; saf öngörü log3/log2 = 1.585'in altında). Bu %12'lik ERKEN
  AÇILMA üç adada da sistematik (yarı/τ_ilk 0.870-0.891) ve AÇIKTIR.
  DÜRÜST KAYIT: ζ'nın yarı/τ_ilk = 1.338'i eşik-kırılganlığıdır
  (D(0.068) = 0.489, 0.5'in 0.011 altında; τ=0.075 çukuru geçişi
  0.085'e atlatıyor). Eşiksiz D(τ_ilk) ζ'da 0.461 — bandın içinde.
  YARI-ÇÖZÜLME ÖLÇÜTÜ EŞİĞE KIRILGANDIR; D(τ_ilk) tercih edilmeli.

YENİDEN ÖLÇEKLEME ÇÖKMESİ (x = τ/τ_ilk): eşik bölgesinde sekiz eğrinin
saçılması ham τ'da std/ort 0.5125 iken ölçekli x'te 0.2559 — YARIYA
iniyor (χ₃ hariç 0.2296, ζ hariç 0.2503, ikisi de hariç 0.2141).
İki kervan tek eşiğe kilitleniyor.

SERTİFİKA DEFTERİ (107c, Not 5 §reproducibility'ye hazır): sekiz küme
için n_ham / n_analiz / kalan bölge / K1-K2-K3 ihlalleri / kopya /
maskelenen sıfır / pencere sayısı / ⟨L⟩ / τ_ilk tek tabloda; ayrıca
nihai analiz penceresi t-aralıkları birebir listeli. Maskelenen oran
%0.0 (ζ) - %10.0 (χ₇). χ₇ tek başına K1'in %54'ünü, K2'nin %61'ini,
K3'ün TAMAMINI taşıyor.

Scriptler: 107a (eski dörtlünün yeniden denetimi, 81 s), 107b (yedi-ada
tablosu + figür, 16 s), 107c (sertifika defteri, 9 s).
npz: 107a_temiz, 107b_yedi_ada, 107c_sertifika (+107a_kusurlar.json).
Figür: 101_donmus_koylar_v2.png (panel a: ham τ, "donma penceresi";
panel b: x = τ/τ_ilk çökmesi).

## 104: χ₇'NİN τ=0.04 ANOMALİSİ KAPANDI — İKİ KOPYA SIFIR (25 Ağustos)

101h'den beri açık duran ve 102'nin "dürüst kayıtlar"ında da anılan
tek soru: χ₇'nin D_spont(τ=0.04)'ü band daraldıkça neden YÜKSELİYOR
(0.078 → 0.106 → 0.125, kaçak yok)? Üç hipotez ön-mühürlü ayrıştırıldı;
ikisi öldü, biri kanıtlandı ve kusurun YERİ tek tek sıfırlara indi.

H-seg (104a) KESİN RET. χ₃/χ₅'i χ₇'nin segment profiline, 2 uzun
  parçaya, 6 eşit parçaya bölmek (eş-n) ve ortak n*=24000 üstünde
  m = 1,2,4,6,8 eşit-parça merdiveni: parça sayısıyla HİÇBİR eğilim
  yok, m=1'de oran en büyük. Ortak ω₀ + ortak ω-bandı da yükselişi
  öldürmüyor (χ₇ 1.39×) → τ↔ω hizalaması da suçsuz.
  YAN KAZANÇ 1: "band daraldıkça yükselme" 20+ düzenekte GENEL
  DAVRANIŞ (1.1-1.3×); anormal olan χ₇'nin 1.60'ı kadar χ₃'ün
  0.26'sıydı. YAN KAZANÇ 2: vekil taban V band'dan ve (tek pencerede)
  n'den pratikte BAĞIMSIZ, AİLEYE bağlı.
  YAN KAZANÇ 3 (asıl ipucu): χ₇'nin havuz D'sinin tamamı TEK
  PENCEREDEN (w1, n=11539) geliyor — w1'de D 0.30-0.45, öteki beşinde
  0.02-0.14. χ₃'ün ters anomalisi de kendi w1'inden.

H-kusur (104b, 104d, 104e) DOĞRULANDI — KESİN.
  101b'nin kusur duyarlılığı ilk kez BAND'A genişletildi: temiz χ₃'e
  SADECE 10 SİLME (0.13‰) D'yi 0.018 → 0.134'e çıkarıyor ve ΔD dar
  bantta 1.8× büyük (+0.117 vs +0.063). "Band daralınca yükselme"
  KUSURUN İMZASIDIR.
  Kusurun yeri (104d): χ₇ w1, t = 3081.6688 ve 3081.6695 (fark 7e-4),
  3081.8885 ve 3081.8891 (fark 6e-4) — İKİ KOPYA SIFIR (normalize
  boşluk 0.0008; pencere medyanı 0.974) + ilerideki iki telafi edici
  kayıp. NET sürüklenme sıfır olduğu için 101f'nin sayım sertifikası
  DA 101d'nin düzlük segmentasyonu DA göremedi. Aynı sınıf ikinci
  kusur χ₃ w1'de, t ∈ [3872.5, 3890.2].
  101f'nin listelediği 8 şüpheli χ₇ bölgesinin analiz pencereleriyle
  örtüşmesi 0 — suçlu listede değildi.

H-fizik (104c) RET. Sabit-ω taraması ω ∈ [0.20,0.65], band ±0.004:
  χ₇'nin en güçlü artığı z = +2.00 (ω=0.380), ön-mühürdeki z ≥ 3'ün
  altında ve ince taramada tekrarlanmıyor. Çözünürlük kapısı χ₇'de de
  geçildi (log2 kontrastı 1.55×), yani null gerçek.
  BAĞIMSIZ DOĞRULAMA: D(ω) artık saçılması ζ 0.070 / χ₅ 0.056 (pürüzsüz,
  101i'yle uyumlu) ama χ₃ 0.464 / χ₇ 0.303 — gürültülü iki aile, kusur
  taşıyan iki aileyle birebir aynı.
  (Dürüst kayıt: ω=log7 ölü/canlı çaprazı BİLGİSİZ çıktı — orada
  τ ≈ 0.19 ve beş aile de çoktan çözülmüş, D ≈ 1; doymuş zeminde
  kontrast imkânsız. Ölü-çizgi hükmünün nulli değildir.)

YENİ KAPI — SIÇRAMA KAPISI (104e), sayım sertifikası + düzlük
segmentasyonunun ÜÇÜNCÜ basamağı: r = (d'nin 21'lik kayan ortalaması)
− (801'lik kayan medyanı); temiz 17 pencerede maks|r| = 0.06-0.12,
kusurlu ikisinde 1.99/2.01 → eşik 0.30, pad 600 (sıçramanın yumuşak
kuyruğu ~1000 sıfır sürüyor). Ek: normalize boşluk u < 0.005 (kopya).
Beş ailede toplam İKİ küme işaretlendi; ζ, β, χ₅ tertemiz.

TEMİZLENMİŞ SONUÇ: χ₇ 0.078/0.106/0.125 (1.60×) → 0.033/0.037/0.038
(1.15×); χ₃ 0.069/0.038/0.018 (0.26×) → 0.004/0.003/0.006 (1.52×, ama
vekil taban 0.0045 ile aynı mertebede = ölçüm sınırında donmuş).
ANOMALİ KAPANDI. 101d/101h'nin ANA HÜKMÜ ETKİLENMEDİ: β log3/⟨L⟩ =
0.1143'e dek donuk (0.013-0.091) ve tam orada 0.847'ye sıçrıyor;
ζ/χ₅/χ₇ log2/⟨L⟩'de çözülüyor. "Donma sınırı adanın ilk sağ kalan
çizgisidir" ayakta.

DÜRÜST KAYIT — YENİ ALET SINIRI: D = |Σ_w num_w| / Σ_w den_w
pencereler-arası KOHERENT toplamdır; pencere katkıları ters fazda
gelince YIKICI GİRİŞİM olur. Ölçüldü: χ₃ τ=0.068'de pencere-başına
0.301/0.033/0.384 iken havuz 0.109 — üçünün de altında. 101h'nin χ₃
satırındaki dalgalanmalar (0.583→0.320→0.232) bunun eseriymiş.
Eşik-üstü tek tek τ noktaları ±0.2 oynayabilir; EŞİĞİN YERİ oynamaz.

PRATİK KURAL (donuk-taban ölçümleri için, 98-S2 ve 101'in standardına
ek): (1) sıçrama kapısı zorunlu — net sürüklenme sıfır olsa bile
pencere-içi |r| > 0.30 varsa kes (pad ≥ 600); (2) donuk D ancak
D > 3·V ise ölçümdür, altı üst sınırdır; (3) band-daraltma bir
KUSUR TESTİDİR: D(dar)/D(geniş) genel bandı (1.0-1.4×) belirgin
aşıyorsa veriden şüphelen, fizikten değil; (4) havuzlanmış D'yi
daima pencere-başına ayrıştırarak rapor et.

Scriptler: 104a (segment aleti), 104b (kusur/sıkılaştırma), 104c
(sabit-ω taraması), 104d (kusurun yeri), 104e (sıçrama kapısı +
temizlenmiş taban). Toplam koşu ~8 dk. npz: 104a/104b/104c/104e.

## 103: NOT 4 REVİZYONU UYGULANDI (20 Ağustos, gece — Opus tayfası 3. sefer)

Etki listesi (NOT4_ETKI_LISTESI_102.md) tex'e işlendi: 11 değişiklik +
3 kazanç (örnekleme-fazı paragrafı "Whose darkness?", bin-düzeyi Berry
güçlendirmesi, atomik/sürekli tayf ayrımı Synthesis'te) + dürüst
sınır-dili. TDS artık yalnız reddedildiği cümlede geçiyor; tutarlılık
taraması temiz. fig_diffraction_en 3 panele çıktı (103 scripti; orta
panel 5-ızgara karşılaştırması, sağ panel örnekleme-fazı eğrisi —
102a değerleri her basamağa dek bağımsız yeniden üretildi). Derleme
temiz, 10 sayfa. TEK ESASLI SAPMA (onaylandı): 4-18× → 3-7× (medyan
4-5) — kaynak 87b (kaptanın başlattığı bağımsız oturum; 87'nin 15×
ucu Σω taşıyıcısını Δω tabanına bölme artefaktıydı; taban gerçek,
sızıntının ~10⁹ katı — 101j'nin o şüphesi ω≈4-7'de kesin yanlış).
Kalan: kaptanın son okuması; abstract'taki "hyperuniformity seen
directly" ibaresi gövdeden hafif yumuşak (son okumada karar).

## 106: İKİ-DAL GİZEMİ ÇÖZÜLDÜ — BİRİ ÖZDEŞLİK, BİRİ OPTİK, FİZİK PERDEDE (24 Ağustos)

Kota tazelendi; 104 (χ₇) ve 105 (çift karakter + mod 8) Opus tayfalarında,
kaptan köşkü 106'da. Belirleyici deney: öz-tutarlı (ρ̄u=−S(t+u)) vs naif
(u=−S(t)/ρ̄) lab örgüleri + 101d spontane estimatörü + 101i vuruş tarama.

ÜÇ HÜKÜM (106 + 106b):
1) PERDE = GERÇEK ÇOK-CİSİM FİZİĞİ: öz-tutarlı lab v-genlikleri katı-
   analitiğin %±5'inde DÜZ — perde (0.91→0.48) lab'da YOK. Öz-tutarlılık
   hipotezi RET (naif-oran kanalı 671 kesişmeyle kirli; yükselen oran
   artefakt olarak kayıtlı). 92'nin açık sorusu temiz izole edildi:
   perdeyi üreten mekanizma EF-yerdeğiştirme alanının ötesinde.
2) SPONTANE "CUE-UYUMU" = ÖRNEKLEME OPTİĞİ: lab (CUE'suz!) spontane
   eğriyi τ≥0.15'te %2-6 içinde üretiyor (1.041/1.112/1.192 vs
   1.075/1.088/1.177). 92-T4'ün "ζ spontane modlara CUE gibi yanıt
   verir" cümlesi yeniden yorumlandı; Not 4 §branches'e son-okuma
   kaydı (etki listesine EK yazıldı).
3) DONMA = ATOMİK TAYFIN KESİNLİĞİ: lab düşük-ω'da gerçeğin 109 katı
   parlak; lab SIFIR-ızgarası bile parlak (4.7e-4 vs gerçek ~7e-16).
   Gerçek boşluklar TAM boş çünkü açık formül kesin özdeşlik; lab'ın
   t+u inşası u² vuruşlarını sahte üretir (lab boşluk modeli değil;
   orta-ω yan-bant geçerliliği [102b %6] sürer). Donmanın "çizgi
   envanteri" okuması (101) bu temele oturdu: boşluğun boşluğu
   özdeşliktir, çözülme çizginin girişidir.
Belirsiz: lab vuruş-tepe kontrastı yok (P3'). Scriptler: 106, 106b.

## 105: ÇİFT KARAKTER + MOD 8 — DÖRT ÖN-MÜHÜR DE İSABET (24 Ağustos, Opus tayfası)

Üç yeni ada: chi5e (mod 5 Legendre, İLK ÇİFT karakter a=0), chi8e
(mod 8 çift), chi8o (mod 8 tek) — motor a=0 sağlaması geçti (yanlış
parite 38 mertebe patlatıyor → kapı duyarlı; a=0 dalı konumda a=1'den
iyi). Kampanya + tam sertifika: 72.8k/65.9k/65.9k sıfır; kaldırılmış-dip
dersi birebir tekrarlandı (+376/+130/+432 mpmath kurtarması).

P2 ✓✓✓ SPEKTROSKOPİ: fazlar 0.04-0.11° içinde; mezarlıklar (chi5e:
  5,25; chi8'ler: 2,4,8) binde-1 düzeyinde ölü.
P1 ✓✓ GAMMA EVRENSELLİĞİ: mutlak yasa oranı chi5e (a=0) 1.009 = chi5
  (a=1) 1.009 (aynı q, aynı L!); chi8e/chi8o |Ĝ| ve kuvvet-açığı
  birebir, yalnız fazlar arg χ kadar ayrık; hiperuniformlukta a-farkı
  ≤0.01 (kaba Σ²(n) ölçümü). Gamma faktörü HİÇBİR yasaya girmiyor.
P4 ✓✓✓ β-SOĞUKLUĞU ÇÖZÜLDÜ: σ_tam²@L=9 — canlı-2 dörtlüsü 0.058-0.064,
  ölü-2 üçlüsü (β, chi8e, chi8o) 0.0496/0.0500/0.0500: iki iletken
  (4 ve 8), iki parite, AYNI −%20 artık. İletken-4/imprimitiflik
  hipotezi ÖLDÜ; soğukluk = 2-ailesinin ölülüğü. (Koro-varyansı artığı
  tam açıklamıyor — nicel artık-yasası açık; ama artık yalnız çizgi
  envanterine bağlı.)
P3 ✓✓✓ DONMA İKİ BAĞIMSIZ ADADA REPLİKE: chi5e log2/L'de çözülür
  (yarı-çözülme/τ_ilk = 1.006, χ₃ ile birebir); chi8e VE chi8o log3'e
  gecikir (chi8e log2'de D=0.009 — vekil tabanında: ölü çizgide çözülme
  YOK); parite farkı %1.3. HÜKÜM: donma sınırı YALNIZ çizgi
  envanterinin fonksiyonu — parite/gamma/iletken değil.

105e DÜRÜST KOVALAMACA: chi8e eşik-altı anomalisi kovalandı → YENİ
KUSUR SINIFI: kısa çukur (2 kayıp + anında telafi) — 80'lik medyan
dedektörü bunları YIKAR, donma estimatörü aşırı duyarlı (0.02→0.48).
Yeni kapı: kısa-20/uzun-200 medyan farkı > 0.7 (yanlış-pozitif temiz).
Anomali aletseldi (0.188→0.045). AÇIK: χ₃'ün orta-τ çukuru yeni kapıyla
derinleşti → 101f verilerinin (4 eski ada) kısa-çukur kapısıyla yeniden
denetimi sıradaki iş; hiperuniformluk ölçümü kaba (102a titizliğinde
değil). Scriptler: 105a-e + npz'ler.

## 108: ADA PERDESİ — EVRENSEL, SICAKLIK-BAĞIMSIZ, SEÇİCİ (24 Ağustos, kaptan)

Perde-teorisi kuyusunun ölçüm ayağı, üç script:
108 (tek-çizgi) TUZAĞA DÜŞTÜ ve ζ kontrol sütunu yakaladı: D≈1.05
  her yerde — 81-83'ün örnekleme-ızgarası şişmesi. DERS TAZELENDİ:
  perde YALNIZ tam-taban ortak regresyonla ölçülür.
108b (tam-taban, 8 küme): alet doğrulandı (ζ 2τ-normda 0.94→0.54 =
  kampanya). KAVRAM DENETİMİ: sinüs-normda perde 0.955→0.823 —
  "0.48'e düşüş"ün çoğu 2τ/sinüs geometrisi. P3 ✓✓✓ SEKİZ SÜTUN TEK
  EĞRİ (±0.01-0.02): PERDE EVRENSEL. P2 RET: ölü-2 üçlüsü canlılarla
  özdeş (ln-oran ~1.0) → sıcaklık-DW hipotezi ÖLDÜ; perde aritmetik
  sıcaklığa değil yerel/evrensel istatistiğe bağlı.
108c (enjeksiyon kalibrasyonu): a_v(τ) ≈ 1.00 (ζ ve χ₃ özdeş) →
  zincir temiz; D_fiz = 0.958→0.797, ~1−0.36τ. BONUS — SEÇİCİLİK:
  enjekte yapay dalga AYNI frekansta PERDELENMİYOR → perde yalnız
  gazın KENDİ aritmetik dalgalarını bastırır; yarım-gap-jitter DW
  adayı da öldü (yapayı da bastırırdı).
GÜNÜN PERDE BİLANÇOSU (106+108): perde (a) EF-ötesi gerçek fizik,
  (b) evrensel, (c) sıcaklık-bağımsız, (d) alet-temiz, (e) SEÇİCİ.
  Kalan tek mekanizma sınıfı: gap-dalgalanmalarının aritmetik dalga
  FAZIYLA korelasyonu — gazın dalga çevresinde yeniden-dengelenmesi.
  Sıradaki teori ödevi: bu korelasyonu doğrudan ölçmek (faz-çözümlü
  ⟨δg·cos(ω m+φ)⟩ üçlü-korelatörü) ve GUE'de karşılığını hesaplamak.

## 109: GAZIN NEFESİ TERS — ANTİ-ADYABATİK ÜÇLÜ-KORELATÖR (24 Ağustos, kaptan)

Perde mekanizması avının faz-çözümlü ölçümü: dalga-ortalaması tam-tabanla
söküldü, artık η'nın ikinci momentleri çizgi fazında okundu (ζ, iki
büyük pencere; plasebolu; asallar 2-17).

YAN MANŞET (M0): var(ds)=0.164 (GUE-yakın; 41/36/101f üç kaynakta
teşhisle doğrulandı) ama σ_η²=0.022 → GAP VARYANSININ ~%87'Sİ
DETERMİNİSTİK açık-formül dalgaları (Σ A²/2 ≈ 0.14 muhasebesi tutuyor).
GUE-biçimli spacing dağılımı ~60 deterministik dalganın süperpozisyonu
olarak doğuyor. Artık komşu-bağı c₁/σ² = −0.52.

ANA MANŞET (Q2): NEFES DALGASI VAR (Q1 ✓, plasebo temiz) ve TERS —
R_p ≈ −0.85, yedi asalda düz. Adyabatik yerel-ölçekleme R=+1 öngörür;
gaz tam tersini yapıyor: aritmetik dalganın GERDİĞİ bölgeler sessizleşir,
SIKIŞTIRDIĞI bölgeler gürültülenir (%85 güçle). Q3: komşu-bağı dalgası
R_nn −1.9→−1.2. TUZAK KONTROLÜ: dalga-çifti (beat) özdeşliği p=11,13,17
için mevcut değilken R düz → özdeşlik kaynak olamaz (açık hesap
doğrulaması yapılacak).

SIRADAKİ: (a) CUE-referans R (dış dalga enjekteli CUE topluluğu — termal
gaz ne yapar?); (b) adalarda R evrenselliği; (c) ters-nefes ↔ perde
(D≈1−0.36τ) nicel köprüsü; (d) beat-özdeşliğinin açık hesabı.
Script: 109_uclu_korelator.py.

## 109b: HAKİKAT MATRİSİ — "TAM TERS NEFES" (24 Ağustos, kaptan)

Dört hücre: CUE+boyalı / ζ+boyalı / ζ+kendi / lab+kendi. Kinematik
teorem (boyalı dalga → R=+1) küçük τ'da doğrulandı; 2.-moment kanalının
kendi τ-transferi çıktı (S1-S2 referans eğrisi: +0.94→+0.15; 1. momentin
a_v≈1'inden farklı — yeni alet bilgisi). ζ boyalı dalgaya CUE ile
örtüşerek termal davranıyor (+0.95/+0.77/+0.27); KENDİ dalgalarına
−0.85 DÜZ; lab kendi dalgalarına POZİTİF (+0.2..+0.6) → ters nefes
EF-ötesi ve YALNIZ gerçek gazda. TRANSFER-DÜZELTMELİ: R_true ≈ −1 —
gaz kendi aritmetiğine adyabatiğin TAM TERSİ, EŞİT GÜÇTE nefesle yanıt
veriyor. Seçicilik artık iki momentte kanıtlı (108c: 1. moment; 109b:
2. moment). Perde-köprüsü (D≈1−0.36τ ↔ R_true≈−1) sıradaki teori ödevi.
Script: 109b_cue_referans_nefes.py.

## 110: KÖPRÜ DENKLEMİ — AYAKLAR DİKİLDİ (24 Ağustos gecesi, kaptan)

Soru: tam-ters-nefes (R≈−1) perdeyi (D≈1−0.36τ) üretir mi? Kalem
(tek-gap teleskop, 2. mertebe): D = 1/(1+(κ/2)tan(κ/2)σ_η²|R|) —
düzey (0.97) ve biçim (κ²) YETERSİZ → sentetik merdiven.

PERDE 1 — İKİ ÖĞRETİCİ NULL: (a) eklemeli mimari (dalga + jitter)
nefes üretmez; dilatasyon teoremi gürültünün dalganın İÇİNDEN geçmesini
ister (CUE kontrolü +0.94 vs eklemeli 0 çelişkisi yakalattı).
(b) FAZ DERSİ: nefes GAP-dalgası (türev, −sin) fazında yaşar;
yerdeğiştirme (cos) fazına boyanan reçete R'ye DİK kalır — kuadratür
ayrışımı (K1-K6 suç merdiveni) sırrı çözdü. İkisi de mekanizma bilgisi.

PERDE 2 — DÜZELTİLMİŞ MERDİVEN (dilatasyon + türev-fazlı jitter-mod):
β=0 kinematik sağlama ✓ (R +0.93→+0.26, D≈1). β=−2: R −0.86/−0.49/
−0.11, D 0.982/0.961/0.913. β=−3: R −1.78/−1.19/−0.35, D 0.971/0.917/
0.872. HÜKÜM: ANTİ-NEFES PERDEYİ MEKANİK ÜRETİR ve biçim LİNEER-τ
(gerçeğin 1−0.36τ biçimi ✓ — kalemin κ²'si değil). Büyüklük: gerçek
açığın ~%70-75'i (orta/yüksek τ, R eşlenik); düşük-τ ucu eksik.
KALAN AYAK: komşu-bağı kanalı (gerçek R_nn≈−1.5, reçetesiz kaldı) +
iid-dışı gürültü yapısı → M-c basamağı ve kapalı form sıradaki oturum.
Scriptler: 110, 110b. GÜNÜN KAPANIŞI: perde artık karakterize (108) +
faz-uzayı yüzü ölçülü (109-109b) + mekanik köprüsü yarı-kurulu (110).

## 110c: M-c BASAMAĞI — KÖPRÜ %60-80 KURULDU, İKİ KAPALI FORM (25 Ağustos)

25 hücrelik (β,τ) ızgarası + gerçekle hizalama. İKİ KAPALI FORM:
  R_p(β,τ) = πτ·cot(πτ) + β·cos(πτ)   [25 hücrede ±0.05 ✓]
  D(β,τ) = 1 − 0.111·|β|·τ            [artık std 0.008]
πτcot(πτ) ÜÇÜNCÜ kez sahnede (88 benek sönümü, 102c bin-rampa, şimdi
nefes kinematiği) — örnekleme geometrisinin ana fonksiyonu olduğu
kesinleşti. HİZALAMA (gerçek tarafta sıfır serbest parametre: β_eff
ölçülen R'den, c sentetik merdivenden): τ=0.10→0.30 aralığında nefes
köprüsü gerçek perdenin %59→%81'ini taşıyor. Kalan %20-40 (düşük-τ
ağırlıklı): konum-jitter sınıfının bond kilidi gerçeğin düşük-τ
R_nn/R_p oranını (2.75) vermiyor — daha zengin gürültü çekirdeği
(Δ² bileşenleri) ister. τ=0.45 hizalaması geçersiz (cosπτ→0 ıraksama;
gerçek R yüksek-τ'da ölçülmedi — kompozit asallarla ölçüm adayı).
AÇIK: c=0.111'in analitik türetimi + σ-ölçeklemesi.
BÜYÜK RESİM: perde hikâyesi üç kapalı-formlu bir mekanizma zincirine
indi: gaz kendi gap-dalgasının fazında ters nefes alır (R_true≈−1,
yalnız kendi dalgalarında) → nefes, örnekleme geometrisi (πτcotπτ)
üzerinden koheran okumayı lineer-τ yasasıyla söndürür (D=1−c|β|τ) →
ölçülen perdenin ~3/4'ü budur. Not 5'in mekanizma bölümü hazır.
Scriptler: 110c. Commit ile: 110 serisi tamam.

## 111: BOND-ZENGİN ÇEKİRDEK — GAUSSIAN SINIFIN TAVANI %70 (25 Ağustos)

İki-düğmeli (β: yerel genlik-mod; γ: paylaşımlı-çift, yalnız-varyans)
çekirdek, 32 hücre + iki-nüfus doğrulama noktası. E1 ✓ (γ ayrıştırır:
R_p oynar, R_nn sağır). E2: D VARYANS kanalına bağlı — bond D'yi
doğrudan beslemez; bond, zayıflamasız dürüst GÖSTERGE (içsel β'yı düşük
τ'da o verir: gerçek R_nn −1.95 → β≈−3 → D-yasası 1−0.333τ = gerçeğin
%92 eğimi!). AMA tek-β R_p'yi aşırı yapar; iki-nüfus noktası (β=−3,
γ=17: yerel ters nefes + kolektif birlikte-nefes) çifti eşler ve D
açığının %65-78'ini verir. HÜKÜM: faz-modülasyonlu GAUSSIAN jitter
sınıfları ~%70 TAVANINDA — kalan %25-35 sınıf dışı: alt-Gauss yapı
(gerçek kurt −0.75), çarpıklık dalgası ⟨η³⟩(faz), veya kaynak-verteksi.
Sıradaki basamak: 112 (alt-Gauss/çarpıklık merdiveni). Yüksek-τ R_nn
kilit-uyumsuzluğu açık. KÖPRÜ DURUMU: ~3/4 kapalı-formlu kuruldu
(R_p = πτcotπτ + βcosπτ; D = 1−0.111|β|τ; iki-nüfus resmi), son
çeyreğin adresi Gaussian-ötesi yapıya daraltıldı. Script: 111.

## 112: ÇARPIKLIK DALGASI — YENİ ODA AÇILDI, İKİ SÜRPRİZ (25 Ağustos)

Üçüncü momentin ilk faz-çözümlü ölçümü (ζ, iki pencere, plasebolu,
boyalı-kontrollü). S1: η AĞIR-KUYRUKLU — skew +0.24, kurt +0.86 →
111'in "alt-Gauss yapı" adayı η düzeyinde RET (73'ün −0.75'i u içindi,
η'ya taşınmıyor). S2: çarpıklık dalgası VAR ve anti-fazlı (R₃ −1.4..
−2.0; plasebo tabanının 3-6 katı). S3 SÜRPRİZ: boyalı referans da
NEGATİF ve daha güçlü (−4.5/−4.1) — naif adyabatik +1 referansı bu
kanalda YANLIŞ (ön-mühür çerçevesi kısmen geçersiz; dürüst kayıt).
Geçici okuma: gerçek gaz, kinematik çarpıklık yanıtını ~0.4'e BASTIRIYOR
— tanıma 3. momentte işaret-çevirme değil güçlü bastırma. Köprüye
sayılması R₃-kinematiği teorisi yazılana dek beklemede. Kalan-%25
adayları: alt-Gauss RET; 3.-moment yorumu açık; kaynak-verteksi masada.
Script: 112. — Bugünün toplamı (110c+111+112): köprü ~3/4 kapalı-form,
son çeyrek iki daralmış adaya inmiş durumda; üç yeni gözlemlenebilir
(R_p, R_nn, R₃) ve iki kapalı form Not 5'e hazır.

## 113: R₃ KİNEMATİĞİ HAKEM EĞRİLERİ — ÜÇ YASA, BİR YAPISAL KEŞİF (25 Ağu gecesi)

K1: Gaussian kinematiği sıfır değil — R₃(σ³) ≈ −0.3·κ³·σ (σ-lineer,
τ-kübik; katlanma değil komşu-lag [⟨η³δm⟩=3σ⁴] adayı). K2 YAPISAL
KEŞİF: konum-jitter farkı gap-çarpıklığını HER dağılımda öldürür
(⟨(Δw)³⟩≡0) → gerçek η'nın +0.24 çarpıklığı GAP-DÜZEYİ (itme fiziği)
asimetrisi — tüm jitter sınıflarının dışında; kalan-%25 için en somut
iz. K3: ζ+boyalı eğrisi ~düz (−2.7..−3.2, σ³-norm), Gaussian-kinematiğin
çok üstünde; işaret ANTİ-ölçek-sürükleme (naif +2.1'e karşı −2.8).
Gerçek/boyalı oranı ≈0.4 eğri düzeyinde doğrulandı (112 okuması ✓).
KALEM OTURUMU HEDEFLERİ: −0.3κ³σ Isserlis türetimi; anti-sürükleme
kübik muhasebesi; gap-asimetrili gürültü sınıfıyla D-merdiveni.
Script: 113. Kübik gece bereketli kapandı.

## 114-115: KALEM OTURUMU — KÖPRÜ KAPANDI (26 Ağustos sabahı)

114 (kalem+hakem): Gaussian kinematik çarpıklık yasası TÜRETİLDİ ve
DOĞRULANDI — R₃_kin = −(3/4)·σ_ds·κ² (eğrilik-doğrultması; fit'in
⟨N⟩-merkezlemesi katsayıyı 12→6 yapar; yüksek-sinyal hücrelerde %2-7).
MİNİ-TEOREM 1: yasa durağan-Gaussian'da korelasyon-bağımsız (cov(η²,u²)
=2⟨ηu⟩²=0 iptali) → ζ-boyalı düz −2.8 KESİN non-Gaussian. Gerçek η'nın
lag-momentleri ilk kez ölçüldü: skew +0.243, ⟨η²η₊⟩ +0.118,
⟨ηη₊²⟩ +0.119, ⟨η²η₊₊⟩ +0.055, ⟨ηη₊η₊₊⟩ −0.180 (üç-gap itme imzası).

115 (hedef 3): MİNİ-TEOREM 2 — her site-haritası u=f(ε) gap-çarpıklığını
öldürür; çarpıklık ZAMAN-ASİMETRİK çapraz-site yapı ister. En yalın
kurucu: u_n = ε_n + λ·(ε²_{n−1}−σ²) → ⟨η³⟩=6λσ⁴. λ=0.12 skew'i
kalibre etti (lag-deseni kısmen ters — kaba tür, dürüstçe kayıtlı).
MERDİVEN KARARI — GAUSSIAN TAVAN KIRILDI: (λ=0.12, β=−2) hücresi
D = 0.928/0.862/0.886 (gerçek 0.946/0.892/0.838; tavan 0.967/0.933/
0.900) — model gerçeği iki yandan sarıyor (%79-133), sistematik %70
eksiği YOK; aynı hücrenin R₃'ü (−1.24/−1.11/−1.01) gerçek aralıkta.
λ tek başına perdelemez: PERDE = NEFES × ÇARPIKLIK etkileşimi.

HÜKÜM — KÖPRÜ (1. yaklaşımda) KAPANDI: perde iki kanalın toplamı:
  (1) varyans-nefesi (~%70; Gaussian; R_p = πτcotπτ + βcosπτ,
      D = 1−0.111|β|τ kapalı formları),
  (2) çarpıklık-doğrultması (~%30; gap-düzeyi itme asimetrisi;
      β-modülasyonunun λ üzerinden 3. momente yayılması).
Mekanizma cümlesi: "gaz kendi aritmetik dalgasının fazında hem
varyansını hem çarpıklığını ters-modüle eder; bu iki doğrultma birlikte
koheran gap-okumasının 1−0.36τ'sunu üretir." AÇIK inceltmeler: ortak
(β,λ) ince-ayarı; lag-deseni; çarpıklık-katsayısının analitik türetimi;
orta-τ R_p payı. SIRADAKİ BÜYÜK İŞ: NOT 5. Scriptler: 114, 115.

## 116: TERMAL AYNA — TANIMA ARİTMETİĞE ÖZGÜ (26 Ağustos)

Büyük hamle A: denge-CUE'nun KENDİ spontane dalgalarına dörtlüsü
(kovaryans koşullaması, Metropolis'siz; 3000 örnek). T1 ✓ (D_eq
92-T3'ü birebir üretti — makine doğru). KARAR: DENGE GAZI KENDİ
DALGALARINA ADYABATİK NEFES ALIR — R_eq +1.02→+0.47, Rnn_eq pozitif,
R3_eq güçlü pozitif (+3.6 σ³) — ÜÇ MOMENTTE DE GERÇEĞİN TERSİ.
HÜKÜM: ζ-gazının ters nefesi hiçbir denge log-gaz davranışında yok;
aritmetik dalgalar gazın gözünde "kendi termal dalgası" değil — özel
muamele görüyor. Program cümlesi artık termal-aynalı: "denge gazı
kendi dalgalarını adyabatik taşır; ζ-gazı kendi aritmetiğine karşı
gürültüsünü ters örgütler." NOT 6'NIN ANA CÜMLESİ ADAYI. Mekanizma
teorisinin kalan sorusu keskinleşti: statik-adyabatik denge de +1
verir (yerel-evrensellik) — ters işaret için denge-dışı/kural-bağlı
bir ilke gerek (varyans bütçesi? sayım-özdeşliği kısıtı?). Script: 116.
Büyük hamle B (GL(2)/Δ, 117) tayfada sürüyor.

## 117: RAMANUJAN Δ — PROGRAMIN İLK GL(2) ÖRGÜSÜ; HECKE İŞARETLERİ
## KIRINIMDAN OKUNDU (26 Ağustos, Opus tayfası)

Bugüne kadarki sekiz ada GL(1)'di (ζ + yedi Dirichlet). Δ = q·Π(1−q^n)^24
(seviye 1, ağırlık 12) ile derece 2'ye geçildi. Λ(s)=(2π)^{-s}Γ(s)L(s,Δ),
FE s→12−s, ε=+1, kritik doğru Re s=6; θ(t)=−t log2π+Im logΓ(6+it);
ana toplam X(t)=t/2π (analitik iletken (t/2π)²).
KONVANSİYON: yoğunluk θ'/π=(1/π)log(t/2π) ⟹ L_eff = 2·log(t/2π).
(Harekât mühründeki (2/π)·log yoğunluğu 2 kat fazlaydı; N(36000)=87 700
öngörüsü ölçülen 87 543 ile %0.2 uyuştu — konvansiyon böyle sabitlendi.)

MOTOR (117a) — ve seferin ilk sürprizi: TAPER.
τ(n) kesin (Jacobi E³ seyrek çekirdeği, 7 çarpım, 0.3 sn; Deligne 783
asalda ✓). BAĞIMSIZ ORACLE bulundu: Hecke'nin kesin integral formülü
Λ = Σ τ(n)[Γ(S,2πn)/(2πn)^S + Γ(12−S,2πn)/(2πn)^{12−S}] — e^{πt/2}
iptali yüzünden t≲250'de ama TAMAMEN bağımsız. Oracle: (a) Λ tam
gerçek ⟹ ε=+1 ve FE doğrulandı; (b) γ₁=9.222379399 (LMFDB ile 4.6e-10);
(c) KESKİN kesim %5.88 hata, KOSİNÜS-TAPER (W=0.5√X) %0.51 — 11-14 KAT.
(d) kesim uzunluğu ölçüldü: u=0.7→0.396, u=1.0→0.008, u=1.3→0.289 —
derece-2 uzunluğu t/2π'de keskin optimum (iletkenin doğrudan ölçümü).
SONUÇ: "kaldırılmış dip" GL(1)'de de ALETSELMİŞ — derecenin değil
KESİM PROFİLİNİN eseri. Kanıt S2'de canlı: t∈[10000,10400]'de keskin
motor 4 sıfır kaçırıyor (93 basamak bayrağı), taper 0 basamak.
⟹ AÇIK İŞ: 101/105'in yedi GL(1) adasının taper'lı yeniden dökümü.

KAMPANYA (117a2): T∈[200, 36000], 87 543 sıfır, 3.7 dk. Sayım
sertifikası +2, dip-kurtarma +0, temizlik +0, KALAN BÖLGE 1.
GL(1)'de ada başına 128-435 kayıp dip vardı; burada SIFIR. Öngörü
mühürlüydü ve isabet etti.

KIRINIM (117b) — TAÇ SORU CEVAPLANDI.
Kod kapısı: aynı boru hattı ζ'da 89'un yasasını birebir veriyor
(faz 180.0°±0.4°, oran 0.99-1.03) ⟹ okumalar kalibre.
P1 ✓✓✓ HECKE İŞARETLERİ 14/14, iki ızgarada da. τ(p) işareti
  (p=2..43: −++−+−−+++−−+−) kırınım FAZINDAN hatasız okundu;
  maks |Δfaz| 0.24° (p=43 hariç — |a_43|=0.018, neredeyse ölü çizgi,
  5.8° sapıyor ama işareti yine doğru).
P2 ✓✓✓ SATAKE p^k 9/9. b(p^k)=α^k+β^k özyinelemesi hem ÇEVRİK hem
  ÇEVRİLMEMİŞ öngördü: 4,9,25,27,49,121,169 → 0° (ζ'da 180°'ydi);
  8 (b=+1.442) ve 16 (b=+0.954) → 180°, ÇEVRİLMEDİ. Faz b'nin
  İŞARETİNİ izliyor, "hepsi çevrik" değil.
MUTLAK BENEK YASASI DERECEDEN BAĞIMSIZ ✓✓: |Ĝ|=|b(q)|Λ(q)/(L_eff√q)
  ·cos(πτ)·DW, üst pencerede oranlar 1.003-1.035 (20 çizgi) — GL(1)
  kalitesi. Tek değişen: açık formül katsayısı 1 → b(q).
MEZARLIK ✓✓✓: bileşik ω (log 6,10,14,15,21,22,33,35) plasebo
  tabanında (3e-5…6e-4 vs canlı 0.008-0.047) — "her şey çizgidir"
  GL(2)'de doğrudan.
P3 — KORO YASASI KIRILDI VE YENİDEN KURULDU. Ön-mühür "Δ eşleştirilmiş
  L'de ζ'dan soğuk" dedi: İSABET (σ_u 0.2129 vs 0.2858 @L=12.66;
  0.2342 vs 0.3055 @L=16.10). Ama 97/105'in TOPLAMSAL biçimi
  (Δσ² = (2/L²)(Σa_p²/p − Σ1/p)) büyüklüğü 3.1-5.3 kat ıskaladı ve
  L ile yanlış yöne gitti — RET. Veriden ÇARPIMSAL biçim doğdu:
     σ_u²(Δ)/σ_u²(ζ) = [Σ a_p² w_p/p] / [Σ w_p/p]   (w = DW ağırlığı)
  %5.2 ve %2.2 sapma, L ile İYİLEŞİYOR. GL(1)'de iki biçim ayırt
  edilemiyordu (envanter yalnız birkaç ölü çizgi kadar değişiyordu);
  Δ'da Sato-Tate BÜTÜN çizgileri yeniden ağırlıklandırdığı için ilk kez
  AYRIŞTILAR. Fizik: Δ soğuk çünkü küçük asallarda a_p²≪1 (0.281,
  0.358, 0.142) ve koro ağırlığı 1/p; ⟨a_p²⟩=0.9775 (Sato-Tate ✓).
  β'nın "ölü 2-ailesi soğukluğu"nun SÜREKLİ analoğu.
Σ²(n) = 0.504/0.460/0.443/0.470 — GUE'nin 0.76-1.31 katı, GL(1)
  bandıyla aynı.

ARTEFAKT AVI (gizleme yok — beşi de kayıtlı):
 1) G3'ün ilk parantezi ±0.3⟨g⟩ idi; t≈146'daki YAKIN ÇİFT
    (γ=146.1487/146.4091, s=0.26) parantezde iki kesişim bıraktı,
    sahte 2.6e-1 "hata" üretti. Oracle o bölgede motorla 1e-3 içinde.
 2) S1 ölçütü "1.3X ile fark<%10" idi — YANLIŞ KURULMUŞ ölçüt: AFE
    uzunluğu iletkenle sabittir, uzatmak bozar. (b2) bunu ölçtü.
 3) S2 ölçütü "sürüklenme aralığı<1.2" idi — yanlış istatistik;
    d_i sabit değil S(t) dalgalanmasıdır (1000 sıfırda ~4σ≈1.8 normal).
    Ölçüt basamak dedektörüne çevrildi.
 4) Σ²(n) ilk koşuda 0.70/1.13/3.02/16.14 çıktı (GUE'nin 26 katı!) —
    kesilen kusurlu bölgenin BOŞLUĞU kutuları boş sayıyordu.
    Segment-bazlı ölçümle düzeldi.
 5) ζ kontrol örgüsü keskin kesimle üretilmişti, 40 birim sürüklenme
    (kayıp sıfır ζ'yı SAHTE ISITIR ⟹ "Δ soğuk"u sahte güçlendirirdi).
    Taper'lı yeniden döküm: sürüklenme 40.0→2.4, 38 sıfır geri geldi,
    σ_u 0.2856→0.2858. Etki 2e-4, fark 0.073 — HÜKÜM AYAKTA.

AÇIK İŞLER: (i) GL(1) adalarının taper'lı yeniden dökümü; (ii) SIFIR
ızgarasında (orta-noktasız) mutlak yasa oranları ω ile büyüyor
(1.03→2.42) — o ızgaranın sönümü Gauss DEĞİL, yeni gözlemlenebilir;
(iii) çarpımsal koro yasasının türetimi; (iv) Δ'da donma/çözülme
ölçümü (117 bu seferde YAPILMADI — ilk sağ kalan çizgi log 2, ölü
çizgi yok, o yüzden 101'in müdahaleli deneyi Δ'da kurulamıyor;
GL(2)'de ölü çizgi için CM formu veya ölü Euler çarpanı gerek).
Scriptler: 117a_delta_motoru.py (motor+kapılar), 117a2_delta_kampanya.py,
117b_delta_kirinim.py; veriler 117a_delta_zeros.npz (87 543),
117b_delta_kirinim.npz, 117a_oracle_ornek.npz.

## KALEM DEFTERİ — DENGE-DIŞI İLKE: ÇERÇEVE (26 Ağustos, kaptan)

Termal aynanın (116) bıraktığı soru: her denge referansı +1 verirken
ζ'nın üç-moment ters nefesi hangi ilkeden doğar? Kalem dökümü üç parça:

1) NO-GO LEMMASI (yerel-evrensellik ⟹ adyabatik): dalga-boyundan kısa
pencerelerde gap istatistiği "yerel ortalama aralıkla ölçeklenmiş
evrensel yasa" ise, ⟨η²⟩(m) = σ²(1+W)² zorunlu → R_p = R_nn = R₃ = +1
(TÜM momentler adyabatik). Ölçülen (−0.85/−1.5/−1.4) ⟹ YEREL
EVRENSELLİK, ARİTMETİK ALANIN FAZINDA O(1) GÜÇLE KIRILIYOR: gazın
yerel istatistiği yalnız yerel yoğunluğun fonksiyonu DEĞİL — aritmetik
fazı taşıyor. "Gaz kendi aritmetiğini tanır"ın mikroskobik içeriği bu.
(Statik-adyabatik, boyalı-kinematik ve spontane-denge referanslarının
üçünün de +1 vermesi lemmanın üç ampirik köşesi: 109b-S1/S2, 116.)

2) YENİDEN-ÇERÇEVELEME — AŞIRI-TERMAL SABİTLENME: aritmetik modlar
denge dalgalanması değil; kıyas: koheran çizgi genliği A≈0.1 (ds),
aynı modun termal RMS'i √(σ_η²/n)·√S ≈ 7e-4 → çizgiler ~100× termal
GENLİKTE (~10⁴× eşbölüşüm enerjisi) SABİTLENMİŞ. Soru "gaz dengeyi
neden ihlal ediyor" değil; "katı log-gaz, aşırı-termal sabitlenmiş bir
mod çevresinde nasıl düzenlenir" — parametrik/nonlineer rejim,
pertürbatif denge değil. (Açık formülün KESİNLİĞİ = modların genlik
dalgalanması SIFIR — denge topluluklarında modlar nefes alır, ζ'da
alamaz; artık serbestlik yalnız gürültüde.)

3) BAĞ ADAYI — KAYIPSIZLIK/KORUNUM: 68'in korunum yasası
(√w+(β/2)v=1; g≈−(β²/4)v²) genlik-kanalında "birinci moment büyükken
ikinciden eksilt" yapısındaydı; 70-73'ün kayıpsızlık/Hermitsel koku
gözlemleriyle birlikte, ters nefes bu korunumun gap-kanalı sureti
olabilir. Türetim açık.

MÜHÜRLÜ AYRIŞTIRICI ÖNGÖRÜLER (sıradaki ölçümler):
  Ö1 ADALARDA R: sabitlenme-çerçevesi R'nin ada sıcaklığıyla (ölü-2
     üçlüsü −%20) oynayabileceğini söyler; katı-evrensellik aynı kalır
     der. (D evrenseldi [108b]; R henüz adalarda ölçülmedi — ayrıştırıcı.)
  Ö2 GENLİK-KANALI NEFESİ R_w: 41 verisinde amps var — |Z| gürültüsünün
     varyansı gap-dalgası fazında nasıl salınıyor? Korunum-çerçevesi
     belirli bir çapraz-kanal deseni öngörür (gap-gürültüsü kısılırken
     genlik-gürültüsü?); tamamen yeni gözlemlenebilir.

## 119: ÇAPRAZ-KANAL NEFES TAKASI — H-KORUNUM KAZANDI (26 Ağustos)

Kalem defterinin Ö2 hakem ölçümü: |Z| genlik-gürültüsünün varyansı,
gap-dalgası fazında İLK KEZ dinlendi. SONUÇ: R_w GÜÇLÜ POZİTİF
(+6.26 → +0.94, p=2..17; plasebo 0.05-0.18) — R_p aynı koşuda
−0.69..−0.92. Aritmetik dalganın gerdiği bölgede gap-gürültüsü
kısılırken GENLİK-GÜRÜLTÜSÜ KÜKRÜYOR: iki kanal, gazın kendi
dalgasının fazına kilitli DALGALANMA TAKASI yapıyor. 68'in korunum
yasası ve Not 1'in gap-genlik korelasyonu gürültü düzeyinde doğrulandı.
Denge-dışı ilkenin adı: SABİT DALGALANMA BÜTÇESİ (faz-çözümlü korunum,
kayıpsızlığın 2.-moment sureti). Yapı notu: R_w hiperbolik-yakın düşer
(R_w·τ 0.44→0.26). AÇIK: hangi büyüklüğün korunduğunun türetimi;
R_w biçim yasası; Ö1 (adalarda R) mühürlü ayrıştırıcı duruyor.
Script: 119. — Not 6 malzemesi (bugün): termal ayna + GL(2) Hecke +
çarpımsal koro + no-go lemması + NEFES TAKASI.

## 120: Ö1 — NEFES YEDİ ADADA EVRENSEL; R_nn ≈ −2cos(πτ) (26 Ağustos)

Mühürlü ayrıştırıcı Ö1 ölçüldü: ⟨R_p⟩ = −0.880..−0.941, ölü-2 üçlüsü
canlılardan AYRIŞMIYOR; hepsi ζ ile uyumlu — NEFES DE (D gibi)
EVRENSEL → sabit dalgalanma bütçesi yalnız evrensel büyüklüklerden
kurulu; sıcaklık/çizgi-envanteri katsayıya girmiyor. YENİ KAPALI FORM:
R_nn ≈ −2·cos(πτ) yedi adada ve ζ'da (bond kanalı kinematik-temiz;
İÇSEL β = −2.0 EVRENSEL SABİT). YAN BULGU: σ_η (0.149-0.157) ve
c₁/σ² (−0.58..−0.60) yedi adada özdeş — artık-gürültü evrensel;
−%20 sıcaklık farkı tamamen deterministik dalga payında (tutarlılık).
Türetim hedefi: "β=−2 neden evrensel?" (iki kanal/iki kuadratür
bütçesi −1+−1?). Script: 120. Günün beşinci büyük bulgusu.

## 118: YEDİ GL(1) ADASININ TAPER'LI YENİDEN DÖKÜMÜ — 104/105/107'NİN
##      BÜTÜN KUSURLARI TEK BİR ALET DEĞİŞİKLİĞİYLE YOK OLDU (26 Ağustos)

117'nin açık işi ("101/105'in yedi GL(1) adasının taper'lı yeniden
dökümü") yapıldı. Sonuç 117'nin öngördüğünden geniş çıktı: 104, 105e ve
107a'nın ÜÇ AYRI SEFERDE ELLE kovaladığı yedi adlı kusurun YEDİSİ DE
ALETSELMİŞ ve keskin AFE kesimi düzeltilince ortadan kalktı.

118a — TAPER'LI GL(1) MOTORU + YER GERÇEĞİ (645 s, dört kapı da ✓)
  Taper 117a kalıbı; tek değişen X = √(qt/2π) (GL(1) uzunluğu).
  G1 ✓ yedi adada |Δθ|=0, |Δε|=0, c=0 ⇔ 98 BİT-BİT aynı.
  G2 ✓✓ ORACLE GL(1)'DE 117'DEN GÜÇLÜ: mpmath Hurwitz her yükseklikte
     kesin Z veriyor (117'de oracle t≲250'de ölüyordu). Keskin hata
     %1.34-4.34 → taper(c=0.5) %0.14-0.46; KAZANÇ 4.9-12.5 kat.
     Optimum c GL(1)'de 0.7-1.0'a kayıyor (117'de 0.5) — c=0.5 tek
     konvansiyon için korundu; c=0.7 ~1.5× daha kazandırır (AÇIK İŞ).
  G3 ✓✓ konum doğruluğu, mpmath köküne karşı: medyan |Δγ|/⟨g⟩
     3.5e-3…1.6e-2 → 2.1e-4…6.6e-4 (11-50 kat). 105a-G5'in ~1e-2'si
     ~1e-4'e indi.
  G4 ✓✓✓ ASIL SINAV: 101e/101f/105b'nin mpmath dip-kurtarmasıyla ELLE
     çıkardığı 2499 sıfırın %98.7'si (2467) taper motorunda SALT
     İŞARET-TARAMASIYLA geri geldi. Keskin motor aynı kümede 1/2499.
     Ada başına: chi3 98.6, beta 97.5, chi5 98.1, chi7 98.5,
     chi5e 98.9, chi8e 99.2, chi8o 100.0.
  BULUNAMAYAN 32'NİN TEŞHİSİ (ve bir hipotez ÖLDÜ): "artık sorun
     ızgara adımı" dedim; ÇÜRÜDÜ. 18 kümenin yalnız 5'inde taper-Z
     işaret değiştiriyor; 13'ünde DİP HÂLÂ KALDIRILMIŞ (min|Z|
     2e-4…1.5e-2). grid_frac 0.03→0.01 yalnız χ₃'te +2 getirdi,
     β/χ₇'de HİÇ. ⟹ TAPER MEKANİZMAYI YOK ETMİYOR, NÜFUSUNU ~60 KAT
     AZALTIYOR. 117'nin Δ'da 0 çıkması GL(2)'nin şansıymış.
     Ayrıca 5'i 107c'nin KOPYA envanterindeki sahte çiftler — taper'ın
     onları üretmemesi kaçırma değil düzeltme.
  DÜRÜST KAYIT: G3'ün ilk ölçütü ("önbellek sıfırlarında |Z|")
     DAİRESELDİ — o sıfırlar zaten keskin motorun kökleri, keskin
     orada tanım gereği 1.7e-11 veriyor. Bağımsız yer gerçeğine
     çevrildi. (117'nin S1/S2 düzeltmeleriyle aynı sınıf hata.)

118b — YEDİ ADANIN DÖKÜMÜ (59 s; boru hattı 105b'nin A/B/C/D'si aynen)
  DİP-KURTARMA İHTİYACI ÇÖKTÜ: C adımı ada başına 128-435 → 2-8 sıfır;
  mpmath eval 217-741 → 3-12 (60 kat). Toplam sayılar eski dökümle
  0-6 sıfır farkla örtüşüyor (chi8o birebir 65943) — iki BAĞIMSIZ yol
  aynı sayıya varıyor.
  KALAN ŞÜPHELİ BÖLGE 18 → 0 (yedi adanın yedisinde).
  ÜÇ KAPI DA SUSTU: K1 ham 85/254/789/169/84/85/0 → hepsi 0;
  K2 24/23/121/20/0/12/0 → hepsi 0; K3 yalnız χ₅'te 2 (yeni bir kopya
  çifti, t=34633.32; kapı yakaladı, analiz dışında).
  maks|med₂₀−med₂₀₀| 1.85/2.01/1.97/1.90/0.37/1.40/0.20 →
  0.20-0.29 (eşik 0.7). Maskelenen: χ₇ %10.0→%4.0, χ₃ %5.2→%3.4,
  χ₈ᵉ %5.4→%4.0, β %4.1→%3.6; kalan maskeleme artık kusurdan değil,
  n<3000 düşen en alt log-penceresinden.
  ★ ADLI KUSUR DENETİMİ — SEFERİN EN ÇARPICI SATIRI: 104d/105e/107a'nın
  elle bulduğu YEDİ kusurun yedisi de artık analiz İÇİNDE, çünkü
  KUSURUN KENDİSİ YOK. Ölçüldü (±400 sıfır komşuluğu, maks|Δmed₈₀| ve
  maks|m₂₀−m₂₀₀|): χ₇ 3081.67: 0.49→0.05 / 1.97→0.13; χ₅ 5820:
  2.00→0.10 / 2.01→0.19; χ₃ 3872.5: 0.18→0.10 / 1.85→0.15;
  χ₈ᵉ 19352: 0.13→0.09 / 1.40→0.16. χ₇'nin altı "kopya sıfırı" ve
  chi8e'nin biri de gitti.

118c — TAÇ SONUÇLAR (96 s) — DÖRT SORUNUN HÜKMÜ
  C1 ✓✓ KERVAN AYRIMI AYAKTA. D(log2/L): LOG2 ζ .461 χ₃ .533 χ₅ .516
     χ₇ .512 χ₅ᵉ .492 | LOG3 β .032 χ₈ᵉ .025 χ₈ᵒ .061. Ortalama oran
     14.3× → 12.8×. D(τ_ilk) LOG3: .742/.800/.722 → .756/.757/.722
     (DARALDI). "Donma sınırı adanın ilk sağ kalan çizgisidir" 5'e 3.
     ÇÖKME İYİLEŞTİ: x=τ/τ_ilk std/ort 0.2559 → 0.2184 (ham τ 0.509
     sabit) — temizlik arttıkça kilitlenme SIKILAŞIYOR.
  C2 ✓/DÜZELTME — ω_yarı ≈ log 2 SAĞ ama "%0.01" TESADÜFMÜŞ.
     χ₃ .6918→.6857, χ₅ .7054→.6834, χ₇ .6792→.6840, χ₅ᵉ .6966→.6966;
     ort .6932→.6874 (log 2 = .69315). Sapma %0.01 → %0.83. AMA
     SAÇILMA YARIYA İNDİ (std .0110→.0062, sem .0055→.0031) ve log 2
     hâlâ 1.9 sem mesafede. Doğru ifade: "~%1 doğrulukla log 2".
     Not 5'te %0.01 YAZILMAMALI.
  C3 ✓ LOG3'ÜN %12 ERKEN AÇILMASI GERÇEK, ALETSEL DEĞİL: %12.0→%11.0
     (log 3'ten 9.9 sem). yarı/τ_ilk .891/.880/.870 → .890/.909/.870.
     Kervan kayması ×1.395 → ×1.422 (öngörü 1.585). "Kaçırılan
     sıfırların kalıntısı" açıklaması ÖLDÜ; mekanizma AÇIK SORU.
  C4 ✓ χ₃'ÜN ORTA-τ ÇUKURU DURUYOR; "ÜÇÜNCÜ KUSUR SINIFI" ELENDİ.
     τ=0.105/0.113/0.125: .199/.161/.134 → .237/.210/.191 (χ₅ .734/
     .754/.741 ve χ₇ .716/.827/.767 ile kıyasla hâlâ 3-4 kat düşük).
     Band ±0.005/0.01/0.02L'de pencere-başına .31/.16, .33/.13,
     .58/.16 — İKİ PENCEREDE DE, ÜÇ BANTTA DA düşük; çukur χ₃'ün
     YÜKSEK-t yarısında (L=9.65, n=63549) yaşıyor. Kusur değil,
     havuzlama değil ⟹ geriye FİZİKSEL şık kaldı. AÇIK SORU.
  ARTEFAKT (gizleme yok): χ₇'nin D(0.095)'i 0.527→0.036'ya düştü.
     C5 teşhis etti — band ±0.01L'de havuz 0.036 iken pencereler
     0.42/0.42: 104'ün YIKICI GİRİŞİM patolojisi. ±0.005L'de 0.361,
     ±0.02L'de 0.602 ⟹ ALETSEL. Taç sonuçları etkilemez (χ₇'nin
     τ_ilk 0.0693, yarı-τ 0.0684). İkinci küçük kayıt: χ₈ᵉ D(0.140)
     .638→.474, D(0.200) .953→.778; pencere oranı 0.90, girişim yok,
     estimatörün bilinen duyarlılık bandında.

DERS (boru hattı için kalıcı): motorlu kampanyalarda AFE'nin KESİM
PROFİLİ, RS-düzeltmesi kadar önemli bir sistematiktir. "Kaldırılmış
dip" bir doğa olgusu değil bir pencere olgusuydu. Kosinüs-taper
(W=0.5√X) GL(1)+GL(2)'de tek satırlık bir değişiklikle sekiz veri
kümesinin kusur envanterini sıfırladı. YENİ STANDART: yeni ada
dökümlerinde taper VARSAYILAN; sayım sertifikası + dip-kurtarma
KALDIRILMIYOR (13/18 kalıntı dip hâlâ onu gerektiriyor) ama artık
seyrek bir emniyet kilidi.

AÇIK İŞLER: (i) c = 0.7 taper'ı (~1.5× daha) ve kalan dar-çift
diplerinin kapatılması; (ii) χ₃'ün orta-τ çukurunun fiziği;
(iii) log3 kervanının %11 erken açılmasının mekanizması; (iv) Not 5'in
sertifika defteri ve ω_yarı ifadesi 118b/118c sayılarıyla güncellenmeli.
Scriptler: 118a_taperli_L_motoru.py (motor+dört kapı), 118b_yedi_ada_
taperli_dokum.py (döküm+sertifika defteri), 118c_tac_sonuc_dogrulamasi.py.
npz: 118a_yer_gercegi, 118b_{ada}_ham, 118b_{ada}_zeros, 118b_sertifika,
118c_yedi_ada.

## KALEM DEFTERİ 2 + 121: "NEDEN −2" → "NEDEN 2× KİNEMATİK" (26 Ağu)

Kalem dökümü üç parça, üçü de 121 hakemiyle onaylı:
(1) TERSİNME ÖZDEŞLİĞİ (kesin): sayım kısıtı ∫ρ=1/gap →
    ds = (W−ν)/(1−W+ν); kinematik nefes TÜM momentlerde +2
    (varyans (1−W)^{-4}, m₃ (1−W)^{-6}). Hakem: yoğunluk-mimarili
    sentetik R_p +1.90/+1.62, R₃ +2.10/+1.61 ✓; boyalı +0.96/+0.71 ✓
    (mimari ayrımı tam 2×).
(2) İNVOLÜSYON R_g + R_ρ = 2 (faz referansı dillerde ters döner):
    sabit nokta +1 = DENGENİN ADYABATİĞİ — termal aynanın (116)
    ilk-ilke açıklaması: denge gazı gap↔yoğunluk ikiliğinin öz-eş
    noktasında oturur.
(3) SORU NİHAİ BİÇİMDE: ölçülen R_g = −2 ⟺ R_ρ = +4 — gerçek gazın
    artık-yoğunluk gürültüsü, yoğunluk dalgasına KİNEMATİK ORANIN TAM
    İKİ KATIYLA biner. "Neden −2" = "neden 2× kinematik"; adaylar:
    iki kuadratür (+2+2), |Z|²-kanalı bağlantısı (R_w>0 ile tutarlı).
BONUS (121): yoğunluk-mimarisi faz-bağımsız gürültüyle bile perde
üretir (D 0.87 @ τ0.3) — perde |modülasyon|-sürücülü, işaret-bağımsız.
Script: 121. Denge-dışı ilkenin kalem arkı bugünlük burada: soru
en yalın biçimine indirildi ve tüm referans noktaları hakemli.

## 122: ZARF SÜRÜCÜSÜ HAKEMİ — ANOMALİ SEKTÖR-SEÇİCİ (26 Ağustos)

Z1 yön isabet üs RET: Var(η|a_u) 6× düşüyor ama eğim −0.60 (kare-ters
zarf zincirinin −2..−4 bandı değil) — zincir üs düzeyinde çürüdü.
Z2 BÜYÜK İSABET: R_w'nun vahşi +6.3→+0.9 düşüşü payda hatasıymış;
doğru paydayla (genlik-kanalının kendi dalga genliği A1_w)
R_w' = +1.21→+0.72 ≈ +1 — GENLİK KANALI KENDİ DALGASINA ADYABATİK.
YENİ KESKİN YAPI: anomali GAP/YOĞUNLUK SEKTÖRÜNE HAS — |Z| sektörü
öz-eş (+1) noktada, gap sektörü −2 ayna-noktasında; aynı dalgalar,
iki sektör, iki rejim. Z3: genlik gürültüsü çarpımsal (doğal).
"Neden −2" iki kez darald: evrensel sabit (120) + sektör-seçici (122).
SIRADAKİ KALEM ADAYI: kısıt-hipotezi — gap'ler sayım-özdeşliğini
(∫ρ=1/gap) TAM taşır, |Z| taşımaz → kısıtlı sektör anomalik, serbest
sektör adyabatik. Script: 122.

## KALEM DEFTERİ 3 + 123: H-KISIT ÖLDÜ — SABİTLENMİŞ GAZ ADYABATİK (26 Ağu)

Kısıt-hipotezi kalemi iki mini-no-go üretti (Gaussian koşullama faz-
kilitli nefes veremez; lineer-yanıt R_eq=+1'i kopyalar) → tek umut
aşırı-termal eyer rejimiydi. 123 DENEYİ (Metropolis, dairesel β=2 +
mod-mıhlama): C₀/σ = 1.9/5/10'da R_p = +1.07/+0.96/+0.95 — ON KAT
termal sabitlemede, %29 dev dalgada bile İNATLA ADYABATİK. H-KISIT RET.
DÖRDÜNCÜ DARALMA: −2 (1) evrensel, (2) denge-dışı, (3) sektör-seçici,
(4) sabitlenmiş-gazla üretilemez. KALAN İKİ ADAY: (a) çoklu-mod
sabitleme (~60 çizgi + p^k kule yapısı birden — modlar-arası yapı;
sıradaki basamak), (b) determinizm/spektral yapı (ζ örgüsü Gibbs
ölçüsü değil deterministik tayfın ergodik izi; −2 iz-formülünün
2.-mertebe yapısı olabilir — derin teori). Script: 123.
Kalem arkının bugünkü bilançosu: dört no-go/RET, iki kapalı form,
bir involüsyon, bir sektör haritası — soru tarihinin en dar hâlinde.
