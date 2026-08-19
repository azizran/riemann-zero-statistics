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
