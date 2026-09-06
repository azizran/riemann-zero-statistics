# 179 — BAĞIMSIZ TEFTİŞ RAPORU

**Beş mührün dış denetimi · Sentez ve nihai hükümler**
Tarih: 2026-09-06 · Teftiş görevi 179 · Dil: Türkçe

---

## 0. BU RAPOR NEDİR

Bu rapor, 169–177 seferlerinde mühürlenmiş **beş büyük iddianın** (S1–S5)
bağımsız teftişinin sonucudur. Teftişin amacı mührü **övmek değil, kırmaya
çalışmaktır**. Rapor tek başına okunacak biçimde yazıldı: her mühür için
iddia, sayıların yeniden-üretim tablosu, iki bağımsız hakemin gerekçeleri,
ve nihai hüküm ayrı ayrı verildi. **Her sayının yanında kaynağı yazılıdır.**

### 0.1 Teftiş düzeneği

Her mühür üç bağımsız mercekten geçti:

| Rol | Ne yaptı |
|---|---|
| **Yeniden-üretici** | Manşet sayıları **birincil veriden** (128_odl_zeros6_2e6_zeros.npz, 165/tayf önbellekleri, 167/C_\*.json) **sıfırdan yazdığı kodla** yeniden üretti. Proje modülü import edilmedi. |
| **Hakem A (istatistik / ön-kayıt merceği)** | Hata çubukları, ön-kayıt bütünlüğü (sha256, zaman damgaları), seçim etkisi, şans tabanı, çoklu-deneme muhasebesi. |
| **Hakem B (sistematik / alet merceği)** | Konvansiyon sızıntısı, pencere seçimi, kestirimci artefaktı, döngüsellik, taze-veri sınavları. |

Ayrıca bir **eksiklik eleştirmeni** beş mührün *dışında* kalan denetimsiz
iddiaları, riskli ölümleri ve seferler-arası tutarsızlıkları taradı (§6).

**Sentezcinin kendi doğrulamaları** §8'de listelendi — hakemlerin en ağır
üç bulgusunu (S1 totolojisi, S2 ön-kayıt eşiği, S4 ön-kayıt formülü) kendi
kodumla/kaynak dosyadan bağımsız olarak sınadım.

### 0.2 HÜKÜMLERİN ÖZETİ

| # | Mühür | Yeniden-üretim | Hakem A | Hakem B | **NİHAİ HÜKÜM** |
|---|---|---|---|---|---|
| **S1** | Vadi türetimi | Aritmetik **birebir**, iddia **yıkık** | ÇÜRÜTÜLDÜ | ÇÜRÜTÜLDÜ | **ÇÜRÜTÜLDÜ** |
| **S2** | A(τ) kimliği | Sayılar birebir, **uyum yok** | ŞÜPHELİ | ŞÜPHELİ | **ŞÜPHELİ** |
| **S3** | İki arkın özdeşliği | **Tam uyum** | ŞÜPHELİ | ŞÜPHELİ | **ŞÜPHELİ** |
| **S4** | Fazlanın ayrışımı | **Tam uyum** | ŞÜPHELİ | ŞÜPHELİ | **ŞÜPHELİ** |
| **S5** | κ çizgi yasası | Sayılar birebir, **uyum yok** | ŞÜPHELİ | ŞÜPHELİ | **ŞÜPHELİ** |

> **Teftişin genel bulgusu.** Beş mührün **aritmetiği kusursuzdur**.
> Yeniden-üreticiler, birincil veriden sıfırdan kurdukları kodla manşet
> sayıların neredeyse tamamını **basılan haneye kadar** yeniden ürettiler;
> uydurulmuş tek bir sayı, oynanmış tek bir zaman damgası, gizlenmiş tek
> bir ölüm bulunmadı. Ön-kayıt betiklerinin sha256'ları tutuyor
> (175a `13a5d195…`, 176a `6173943f…`, 177a `03a9b9ad…`, 177c `1e5cdd1a…`;
> 179/S4 hakem doğrulaması).
>
> **Çöken şey sayılar değil, sayıların üstüne kurulan epistemik
> cümlelerdir.** Beş mührün beşinde de aynı desen tekrarlıyor: bir
> **cebirsel özdeşlik** ampirik bir buluşma gibi sunulmuş; bir **hata
> çubuğu yazılmamış** ve dört hane basılmış; bir **pencere veya
> konvansiyon seçimi** olgunun kendisi sayılmış. Bu rapor bu üç hatayı
> mühür mühür adresliyor.

---

# 1. S1 — VADİ TÜRETİMİ

## 1.1 Mühür cümlesi (iddia)

> c(λ)'nın minimumu **eğim-kesişmesinden türetildi**:
> `dlogKALİB/dλ = dlogW_X/dλ` kesişimi λ = 0.6506 (172) / 0.6486 (173,
> on-nokta) ↔ **bağımsız ölçülen** vadi λ\* = 0.6487; ayrıca Q_E
> tümseğinin tepesi = inşa doyum eşiği λ_c = 1.0109; c(λ) = U + tepe
> (tepe λ ≈ 1.279).

**Kaynak:** `172_carpan_RAPOR.md` §G4, `173_hakem_RAPOR.md` §H4,
`172_configs/172f_sentez.py`, `173_configs/173e_egri_uclari.py`;
veriler `167/C_L*.json`.

## 1.2 Yeniden-üretim tablosu

Yeniden-üretici çizgi-düzeyi kayıtlardan (167/C_\*.json) bant
birleştirmesini, gp-ağırlıklı W_X'i, bant-içi 8-grup jackknife'ı,
geo-ortalamaları ve bütün eğim/kesişim/parabol hesaplarını kendi eliyle
kurdu. **Hiçbir proje modülü import edilmedi.**

| İddia | Kaynak değer | Yeniden üretilen | Fark |
|---|---|---|---|
| 172 §G4.1 eğim-kesişmesi (7 gaz) | λ = 0.6506 | 0.650553 | **0** (birebir) |
| 173 §H4.2 alt kesişim (on nokta) | λ = 0.6486 | 0.648638 | **0** |
| 173 §H4.2 üst kesişim | λ = 1.2787 | 1.278666 | **0** |
| 171 λ\* — *"bağımsız ölçülen vadi"* | 0.6487 (NU kaydı 0.6486573) | 0.648657 | **Sayı birebir — ama bağımsız değil** (§1.3) |
| 173 *"fark 0.0001, ‰0.15"* | 0.0001 / ‰0.15 | **0.000019 / ‰0.03** | Rapor yuvarlanmış hanelerden çıkarma yapmış; **5× abartılı** |
| 172 *"fark +0.0019, ‰3"* | +0.0019 / ‰3 | +0.001896 / ‰2.92 | 0 |
| dlogKALİB/dλ tablosu (172, 6 satır) | −0.9046 … −0.1274 | aynı altı sayı, 4 hanede | 0 |
| dlogW_X/dλ (172, 6 satır, ORTA bant) | −0.6023 … −0.4671 | aynı altı sayı | 0 — **ama tutarsız eşleşme** (§1.3) |
| dlogW_X/dλ (173, 9 satır, 5-bant geo-ort) | −0.6277 … −1.2506 | dokuz satır birebir | 0 |
| Aynı 7 gaz, **tutarlı** tanım (W_X de 5-bant) | *raporlarda yok* | kesişim = **0.648638** | 173'ün "keskinleşme"si yeni gazdan değil, **tanım düzeltmesinden** |
| Alt kesişimin L040'a bağlılığı | 173: "alt uca yeni nokta eklendiğinde" | L040 KALİB'i **%20 bozuldu → kesişim altı hanede değişmedi** (L085/Hkeskin/L115/L130/L140/L145 için de aynı) | **Atıf yanlış**; kesişim yalnız L050/L060/L070'e bağlı |
| Q_E tepesi, 5 orta pencere (0.60–1.15) | 1.0191 | 1.0191 | 0; birini-dışarıda yayılımı **0.0627** |
| Q_E tepesi, 7 nokta (0.50–1.30) | 1.0009 | 1.0009 | 0; yayılım **0.0075** |
| Q_E tepesi, *"9 nokta (0.40–1.45)"* | 0.9407 | 0.9407 ancak **ON** noktayla; dokuz noktayla 0.9391 | **Etiket hatası** (173e `slice(None)` on gazı kapsıyor) |
| λ_c = inşa doyum eşiği | 1.0109 (170 §K0.2) | 1.9147/1.894 = 1.010929; 1.894 iki bağımsız inşa kaydıyla tutarlı | 0 — bu büyüklük **gerçekten bağımsız** |
| c(λ) üst tepesi (1.15/1.30/1.40 parabolü) | λ = 1.2786, c_max ≈ 0.4547 | 1.278577 / 0.454660 | 0; MC 95% aralığı **[1.2253, 1.6144]** |
| Üst dönüşün anlamlılığı | 173 §H4.1: ν−1 = **+5.3σ** | c'nin kendi σ'suyla c(1.30)−c(1.40) = +0.0202 ± 0.0276 = **0.73σ** | **Hata modeli eksik** |
| c defteri (10 gaz) | 0.4098 … 0.4098 | aynı on sayı, maks \|Δ\| = 5e−5 | 0 |
| σ_c defteri (10 gaz) | 0.0049 … 0.0246 | kendi jackknife + bant saçılımıyla aynı on sayı | 0 |
| Q_E defteri (10 gaz) | 0.56708 … 0.55489 | aynı on sayı, maks \|Δ\| = 4.8e−6 | 0 |
| g = 1 − Q/ρ özdeşliği | \|ΔK\|/K ≤ 5.3e−15 | maks 5.3e−15 | 0 |
| α tepesi (172 §G4.2 / 173 §H4.3) | 0.6675 / 0.6594 / 0.7877 / 0.8058 | aynı dördü | 0; yine "9 nokta" etiketi **on** noktalı |

**Yeniden-üreticinin özeti:** *"Bu tarafta tek bir tutarsızlık yok."*
Ama `uyum = false` — çünkü **mühür iddiası** yıkıldı.

## 1.3 Hakem A (istatistik / ön-kayıt merceği) — **ÇÜRÜTÜLDÜ**

1. **TOTOLOJİ.** `c ≡ KALİB/W_X` olduğu için `dlogKALİB/dλ = dlogW_X/dλ`,
   harfi harfine `dlog c/dλ = 0`'dır. Dahası: **herhangi üç nokta için**,
   iki sekant eğimini orta noktalarına koyup doğrusal ara-değerlemenin
   sıfırı, o üç noktadan geçen **parabolün tepesine tam olarak eşittir**
   (20 rastgele üçlüde \|fark\| ≤ 5.3e−14; eşit aralık şartı yok).
   171'in λ\* = 0.648657'si `171k_figur.py:102-106`'da
   `np.polyfit(x, c, 2)` ile alınan **doğrusal** parabol tepesi;
   173'ün kesişimi 0.648638 **aynı üç c değerinin log** parabol tepesi.
   Kutlanan "0.0001 fark" = **0.000020'lik log/doğrusal ölçek farkı**.
   Üst uçta aynı: 1.278666 (log) ↔ 1.278577 (doğrusal).
2. **AŞIRI KESİNLİK.** Raporun kendi σ_tot'larıyla MC (100k çekiliş):
   çekilişlerin **%25.5'i minimum bile üretmiyor**; tepe 68% [0.606, 0.700].
   En iyimser modelde (bant saçılımı tamamen ortak) bile 68% [0.638, 0.665]
   = **±0.014**. `c(0.70) − c(0.60) = +0.02σ` ⇒ **vadi tabanı ölçüm içinde
   düzdür.** "0.0001 isabet" her iki hata modelinde de **100–600 kat**
   aşırı-kesinlik.
3. **ÖN-KAYIT BÜTÜNLÜĞÜ — damgalar tutarlı, sahtecilik yok**
   (ONKAYIT_L040 14:00:34 < C_L040 14:08:47; G4.json 13:32, inşadan önce).
   **Ama** mührün "geçen" sınavlarının hepsi **eldeki eski veriyle
   kararlaşan** sınavlardır. Gerçekten öngörüsel ön-kayıtların **dördü de
   öldü**: "kesişim sayısı = 1" ✗ (2 çıktı), fark(1.375) bandı ✗,
   ν(1.30→1.45) ∈ [0.30, 0.90] ✗ (1.697), 173b-P4'ün dokuz-noktalı Q_E
   tepesi λ_c ± 0.05 ✗ (0.9407, Δ = −0.070).
4. **Q_E ↔ λ_c: SEÇİM ETKİSİ VE AYRIM GÜCÜ YOK.** 36 bitişik pencerenin
   **16'sı** (log; doğrusalda 18) λ_c ± 0.05 içinde tepe veriyor ⇒
   "isabet"in **şans tabanı ~%44–50** ve raporda hiç verilmemiş. Kararlı
   pencere (7 nokta, yayılım 0.0075) **1.0009** veriyor: boş hipoteze
   (tepe referans gazda, λ = 1.0000) 0.0009, λ_c'ye 0.0100 uzak.
5. **ÜST TEPE σ'SU ŞİŞİRİLMİŞ.** `173e:196-205`: ν'nün hatası **yalnız**
   KALİB jackknife'ını taşıyor (±0.1028 birebir üretildi). Bant-başına
   ν = +2.53, +1.92, +1.82, +1.37, +0.52 ⇒ **tam hata σ = 0.349 ⇒
   ν−1 = +1.6σ** (rapor: +5.3σ).

## 1.4 Hakem B (sistematik / alet merceği) — **ÇÜRÜTÜLDÜ**

1. **TOTOLOJİ bağımsız doğrulandı** (kendi kodu, proje modülü yok):
   10-nokta sekant kesişimi 0.648638, aynı üç gazın log-c parabol tepesine
   **8e-16 içinde eşit**.
2. **KONVANSİYON SIZINTISI (yeni bulgu).** Vadi yeri makul konvansiyonlar
   arasında **[0.6209, 0.7002]** geziyor (tek bantlar 0.62–0.70,
   birini-dışarıda-bırak 0.6413–0.6602, λ→logλ apsisi tek başına −0.0024).
   İddia edilen "‰0.15 isabet" (0.0001), **en küçük konvansiyon
   oynamasının 1/24'ü**.
3. **KESTİRİMCİ ARTEFAKTI KOD DÜZEYİNDE.** `172b_gE_yasasi.py:131`
   `WX=orta['W_X']` (**ORTA** bant) ↔ `172d_theta_yasasi.py:138`
   `KAL=exp(mean(log k))` (**5-bant geo**). 172'nin 0.6506'sı bu **tutarsız
   eşleşmeden** geliyor; tutarlı tanım aynı 7 gazda 0.6486 veriyor.
4. **VERİ-KOKLAMA.** 173b'nin "kesişim ∈ [0.63, 0.67]" ön-kaydı, kayıt
   anında **matematiksel olarak garanti** bir "sınav"dı. `172f_sentez.py`
   kendi docstring'inde Ö2'de dürüstçe *"özdeşlik"* diyor; **rapor dili
   bunu bağımsız doğrulamaya çevirmiş.**
5. **YENİDEN-ÜRETİCİYE İKİ DÜZELTME (mührü kurtarmıyor):**
   (a) **Vadinin VARLIĞI sağlamdır** — eşleştirilmiş testte L050→L060
   duvarı **−7.7σ**, L040→L050 **−9.2σ**; jk-yalnız MC'de çekilişlerin
   %100'ü minimum üretiyor. "Hiçbir şey sınamaz" yalnız **konum**
   hassasiyeti için doğrudur.
   (b) Üst tepede bant-çözünürlüklü hesap **+1.5σ ile −2.8σ** arası
   veriyor; yeniden-üreticinin 0.73σ'sı **en tutucu uç**.

## 1.5 Sentezcinin tartısı ve NİHAİ HÜKÜM

**Hakemler çelişmiyor** — ikisi de ÇÜRÜTÜLDÜ dedi ve birbirini bağımsız
kodla doğruladı. Tek nüans farkı **üst tepenin anlamlılığında**:
yeniden-üretici 0.73σ, Hakem A 1.6σ, Hakem B "+1.5σ ile −2.8σ arası".
Bu üçünü tartınca: **hangi hata modeli seçilirse seçilsin +5.3σ düşer**;
gerçek anlamlılık **1–2σ mertebesindedir**, yani işaret düzeyi, mühür
düzeyi değil.

**Sentezcinin kendi doğrulaması:** Totoloji iddiasını kendim sınadım.
`Fraction` ile **kesin rasyonel aritmetikte**, 300 rastgele üçlüde
`max |parabol tepesi − sekant kesişmesi| = 0` (tam sıfır).
Yani bu bir yaklaşıklık değil, **cebirsel özdeşliktir**
(`179/sentez/179_sentez_dogrula2.py`). Ayrıca 173 raporunun literal
ifadesini kaynaktan doğruladım: satır 39 *"171'in bağımsız ölçtüğü
0.6487"*, satır 468 aynısı, satır 799-800 *"171'in bağımsız λ\* = 0.6487'si
ile fark 0.0001"*. **Bu ifade yanlıştır.**

> ## ⛔ NİHAİ HÜKÜM S1: **ÇÜRÜTÜLDÜ**
>
> **Düşen iddia:** *"Minimum eğim-kesişmesinden TÜRETİLDİ ↔ BAĞIMSIZ
> ölçülen vadi λ\* = 0.6487"* cümlesi. İki sayı **aynı üç gazın
> (L050/L060/L070) aynı üç c değerinin** log ve doğrusal parabol
> tepeleridir; aralarındaki 0.000019 bir ölçüm uyuşması değil,
> **ölçek seçimidir**. Sekant-eğim kesişmesi = parabol tepesi,
> herhangi üç nokta için **kesin özdeşliktir** (sentezci doğrulaması:
> rasyonel aritmetikte fark tam 0).
>
> **Düşen sayılar:** dört haneli λ\* = 0.6487 (gerçek çözünürlük
> ±0.014–0.10, yani 100–600× aşırı-kesinlik); "0.0001 fark / ‰0.15"
> (gerçek 0.000019 / ‰0.03); "+5.3σ" üst dönüş (gerçek ~1–2σ);
> "Q_E tepesi = λ_c" isabeti (şans tabanı ~%44–50, kararlı pencere
> boş hipotezi ~11 kat tercih ediyor); 173 §H4b'nin "alt uca yeni nokta
> eklenince keskinleşti" atfı (bozma sınavı: alt kesişim L040'a **hiç**
> bağlı değil).
>
> **Ayakta kalan (gerçek bulgular):** ① Vadinin **varlığı** —
> c(0.40) > c(0.50) örneklem-dışı **+3.4σ**, sol duvar eşleştirilmiş
> testte **−7.7σ**, ν(0.50→0.40) = 2.203 öngörüsel bandın [1.55, 2.45]
> içinde. Minimum [0.5, 0.85] civarında **gerçektir**. ② `dlogW_X/dλ`'nın
> λ ≤ 1.30'da yaklaşık sabitliği (−0.42…−0.63) — **tek gerçek ampirik
> içerik**. ③ λ_c = 1.9147/1.894 = 1.010929 büyüklüğünün kendisi
> Q_E'den gerçekten bağımsızdır (iki inşa kaydıyla tutarlı); sorun
> köken değil, **eşleştirmenin ayrım gücüdür**. ④ On gazın c/σ_c/Q_E
> defterleri ve bütün eğim tabloları **hanesine kadar doğrudur**.

---

# 2. S2 — A(τ) KİMLİĞİ

## 2.1 Mühür cümlesi (iddia)

> Kalibrasyonun λ-değişmez bant şekli `A(τ) = exp(−2π²τ²·σ_X̃(λ=1)²)` —
> **SIFIR serbest parametre**; 5-bant rms %0.62, **ÖRNEKLEM-DIŞI** 9-bant
> rms %1.49 (**1-parametreli en iyi uyumu yener**); geçerlilik penceresi
> λ ∈ [0.50, 1.15]; rakipler (W_amp·W_X, 165-F, tarak-sayımları)
> ön-kayıtlı öldü.

**Kaynak:** `171_nu_RAPOR.md` §T2a, `173_hakem_RAPOR.md` §H4c,
`171_configs/171c_A_kimlik.py` + `171d_A_yuzlesme.py`;
veriler `167/C_*.json` ve `170/K0.json`.

## 2.2 Yeniden-üretim tablosu

| İddia | Kaynak değer | Yeniden üretilen | Fark |
|---|---|---|---|
| A(τ) 5-bant vektörü | 1.1058 / 1.0454 / 1.0000 / 0.9405 / 0.8810 | aynı beşi | **BİREBİR** |
| τ_ort (5 bant) | 0.5393 … 0.6982 | aynı beşi | BİREBİR |
| Kendi KALİB_u2 ↔ JSON bant alanı | — | en büyük bağıl fark **0.000e+00** | Hesap yolu doğrulandı |
| Bant-içi jackknife (5 bant) | 0.60 / 0.68 / 0.41 / 0.77 / 0.81 % | aynı beşi | BİREBİR |
| **MÜHÜR: A2Hk 5-bant rms** | **%0.62** | %0.62 | BİREBİR (dof düzeltmesiyle %0.69) |
| A2s serbest uyum α | 1.146 | 1.14573 | BİREBİR |
| σ\*_eff = √(2α)/π | **0.48190** (rapor §T2a.2) | **0.481845** (A_ONKAYIT.json: 0.4818447623) | ⚠ **RAPORUN YUVARLAMA HATASI** |
| σ_X̃-eşdeğeri = σ\*_eff/2 | **0.24095** (rapor) | **0.240922** (A_ONKAYIT.json: 0.2409223811) | ⚠ **RAPORUN YUVARLAMA HATASI** |
| KİMLİK: σ_X̃-eş ↔ σ_X̃(Hkeskin) = 0.24204 | −%0.46 | −%0.46 | BİREBİR |
| **σ_X̃-eşdeğerinin HATA ÇUBUĞU** | **RAPOR VERMİYOR** | **0.24092 ± 0.00410 (±%1.70)**; σ(α) = 0.0390 (gaz-jk) / 0.0211 (bant-jk) / 0.0350 (gaz s.h.) | ⚠ **EKSİK** — −%0.46 yalnız **−0.27σ** |
| Rakip σ_X̃(son, GERÇEK ζ) elenmesi | +%3.03, "kaybetti" | +%3.03 ama = **+1.73σ** | ⚠ **HÜKÜM AŞIRI** — 2σ içinde, **elenmiyor** |
| λ_eff | 0.9886 | 0.98861 | BİREBİR — ama doğrusal ara değerle 0.9972; ±1σ ⇒ **[0.9549, 1.0238]** |
| ÖRNEKLEM-DIŞI A(0.7386) | 0.8265 ± 0.0122 (5 gaz) | aynı | BİREBİR |
| ÖRNEKLEM-DIŞI A(0.7774) | 0.7456 ± 0.0223 (**2 gaz**) | aynı | BİREBİR (n=2 ⇒ pratikte bilgisiz) |
| **MÜHÜR: A2Hk 9-bant rms** | **%1.49** | %1.49 | BİREBİR — ama **7 NOKTALI** bir rms |
| A2s (1 param) 9-bant rms | %1.58 | %1.58 | BİREBİR — fark 0.09 pp, ölçüm hatası 1.47/2.99 pp ⇒ **AYIRT EDİLEMEZ** |
| Rakip A1 (W_amp·W_X) | ±%9.4, rms9 %10.58 | aynı | BİREBİR — **W_amp'ın yokluğu doğrulandı** |
| Rakip A2iv (λ-değişmez taban) | rms9 %9.07 | aynı | BİREBİR |
| Rakip A3f/A3s (165-F) | rms9 %13.82 / %10.25 | aynı | BİREBİR |
| Rakip A4 (tarak-sayım) | rms5 %79.15 | aynı | BİREBİR |
| Tarak-sayımı İŞARET sınavı: **öngörülen artış** | **×6.0 ARTIŞ** | **×5.224** (betiğin kendi A4 fonksiyonundan; artıklardan çapraz 5.223) | ⚠ **RAPORUN SAYISI −%13 YANLIŞ** (hüküm değişmiyor: işaret hâlâ zıt) |
| Tarak-sayımı: ölçülen oran | ×0.797 DÜŞÜŞ | ×0.7967 | BİREBİR |
| 9-bant SIRALAMA (ilk beş) | A2Hk 1.49 / A2s 1.58 / A2sn 2.19 / A5 3.97 / A1λ 5.64 | aynı beşi | BİREBİR |
| 173 §H4.3 α_g (on gaz) | 0.6474 … 0.1866 | aynı on değer | BİREBİR (≤%0.01) |
| 173 §H4c pencere: rms(0.40) / rms(1.30) | %3.69 / %3.61 | %3.68 / %3.62 | ±0.01 pp — **PENCERE İDDİASI DOĞRULANDI** |
| L145'in sağlık durumu | rapor ✓/✗ hükmü basıyor | **BEŞ BANDIN HEPSİ DÜŞÜYOR** (R_bant 0.9209–0.9598 < 0.98) | ⚠ **USUL SORUNU** |
| İÇ TUTARLILIK: λ=1 gazının **kendi** şeklinin genişliği | RAPOR SORMUYOR | Hkeskin öz-şekli σ\*/2 = 0.23596 = σ_X̃(λ=1)'den **−%2.51** (mührün dayandığı −%0.46'nın **5 katı**) | ⚠ **EKSİK SINAV** |
| Gaz-başına σ\*/2 saçılımı | RAPOR VERMİYOR | L115 −6.31 / Hk −2.51 / L085 +1.12 / L070 +3.21 / L060 +1.88 %; sd %3.89 | ⚠ **EKSİK** |
| SINANMAMIŞ RAKİP: **ölçülen** W_X | RAPOR HİÇ SINAMAMIŞ | rms %0.61 ↔ Gauss %0.62 | Boşluk — **iddiayı destekliyor** |
| Ön-kayıt zaman damgası | A_ONKAYIT 23:44:17 | 23:44:17 → A_YUZLESME 23:45:08 (51 s) | TUTARLI — ama TAU9 **koda gömülü** |

## 2.3 Hakem A (istatistik / ön-kayıt merceği) — **ŞÜPHELİ**

1. **YENİ BULGU (en ağır): MÜHÜR KENDİ ÖN-KAYITLI ÖLÇÜTÜNÜ İSKALADI.**
   `171d_A_yuzlesme.py`'nin ön-mühür ölçütü **Q4: "dokuz bantta rms ≤ %1.2
   tutarsa mühürlenir"** idi. Ölçülen **%1.49**. Mühür **yine de basıldı**;
   gerekçe sonradan seçilen *"A2s'yi yener"* ölçütüne kaydırıldı — ve o
   ölçüt istatistiksel olarak **boştur** (α'lar yalnız %0.93 ayrık;
   örneklem-dışı fark 0.17/0.24 pp iken ölçüm hatası 1.47/2.99 pp).
   **Kale direği taşınmış.**
2. **HAFİFLETİCİ (aynı bulgunun öbür yüzü):** H0 altında ölçülen gürültü
   bütçesiyle **beklenen rms9 = %1.32** (χ²/dof ≈ 1.27). Yani %1.49
   kimliğe **karşı kanıt değil** — **%1.2 eşiği gürültü bütçesi hiç
   hesaplanmadan konmuş.** Q4 ıskası mührü çürütmez, **usul hatasını
   belgeler.**
3. **HATA ÇUBUĞU DAHA DA GENİŞ OLABİLİR:** bant-bırak-1 jackknife
   σ(α) = 0.0574 (σ_X̃-eş için **±%2.5**). "%0.5 içinde" iddiası her hata
   modelinde çözünürlüğün **3–5 kat altında**; σ_X̃(son) rakibi
   muhafazakâr hatayla +1.7σ, bant-LOO ile **+1.2σ — ELENMİYOR**.
4. **SEÇİM ETKİSİ KISMEN TELAFİLİ:** P3 kazanan kimliği (Hkeskin)
   koşudan **ÖNCE adlandırmış** — 8 adaylık tablodan sonradan seçilmemiş;
   bu bak-başka-yere etkisini büyük ölçüde giderir. **Ama** ön-kayıtlı
   "%1 içinde" testinin gücü yalnız **%44** ve −%0.46'lık isabetin şans
   p'si **≈ 0.086** — tek başına "mühür" değil **"işaret" düzeyi**.
5. **DÖNGÜSELLİK YOK:** A(τ) bant şekli ile σ_X̃ marjinal genişliği farklı
   gözlenebilirler; döngü olsa Hkeskin'in kendi şekli kendi σ_X̃'ini
   birebir verirdi — vermiyor (−%2.51).
6. **NİTEL ÇEKİRDEK ÇÜRÜTÜLEMEDİ:** A1 örneklem-içi **10–29σ**, tarak
   sayımı işaret sınavında devasa marjla, 165-F ve λ-değişmez taban
   örneklem-dışı **3–30σ** ile ölüyor. **Hiçbir hata modeli bunları
   diriltmiyor.**

## 2.4 Hakem B (sistematik / alet merceği) — **ŞÜPHELİ**

1. **PENCERE SEÇİM ARTEFAKTI (en ağır bulgu).** `LO_MIN = 0.52` tabanı bir
   sağlık ölçütü **değil**, 159/160'ın karşılaştırılabilirlik
   konvansiyonunun **mirasıdır** (`166_yaris.py`: "görevin τ kısıtı").
   167'nin beş λ-gazının **hepsi taban = 0.40 ile inşa edilmiş** ve
   hepsinde **lo = 0.44 ve 0.48 bantları SAĞLIKLI** (R_bant 1.13–1.51,
   SNR 670–2717) — ama A(τ) analizi bu bantları **hiç kullanmamış** ve
   raporlar varlıklarını **hiç anmamış**. Yüzleştirildiğinde: mühürlü
   A2Hk sapması **−%2.01 (−3.9σ)** ve **−%1.39 (−2.9σ)**; raporun "yendi"
   dediği rakip **A2sn (σ_X̃(son), GERÇEK ζ gazı) −%0.72 / −%0.37 ile
   TUTUYOR.** Sapma 10/10 gaz-bant kombinasyonunda aynı işaretli.
2. **KİMLİK PENCEREYE GÖRE EL DEĞİŞTİRİYOR.** 7-bant (0.44–0.68)
   penceresinde sıralama **tersine dönüyor**: A2sn %0.71 < A2s %0.99 <
   A2Hk %1.06 (4 normalizasyon seçeneğinin 3'ünde A2sn birinci, A2Hk
   **hiçbirinde** birinci). 7-bant uyumundan σ_X̃-eşdeğeri **0.23367**:
   σ_X̃(son) = 0.23384'ten yalnız **−%0.07**, σ_X̃(λ=1)'den −%3.46.
3. **4 ÖRNEKLEM-DIŞI NOKTALI χ²:** A2Hk χ²/4 = **6.21** (p ≈ 5e−5 →
   dışlanır), A2sn χ²/4 = **1.46** (uyumlu).
   **ŞERH:** lo = 0.44 bandı taban 0.40'a 1 bant genişliği mesafede; beş
   gazın paylaştığı **~%1.5'lik ortak bir kenar sistematiği** gaz-arası
   saçılımla dışlanamaz. **Ama bu şerh mührü kurtarmaz:** bant düzeyinde
   %1.5+ ortak sistematik mümkünse "%0.5 içinde" hassasiyet çerçevesi
   zaten çöker.
4. **KONVANSİYON DUYARLILIĞI:** σ_X̃-eşdeğeri makul varyantlar arasında
   **0.23721–0.24433** (%3.0 genişlik) oynuyor; λ_eff aynı varyantlarla
   **0.9580–1.0177**.
5. **SEFERLER-ARASI VERİ-KOKLAMA:** 171c'nin docstring'i ÖN-MÜHÜR'ün
   *"171a'nın mühürlü A vektöründen EL HESABI"* olduğunu **kendisi
   söylüyor**. Dolayısıyla 5-bant rms %0.62 **tamamen örneklem-içidir**;
   gerçek örneklem-dışı içerik yalnız **2 nokta** ve o iki nokta A2Hk ile
   A2sn'yi **ayırt edemiyor** (ikisi de < 2σ).
6. **DÖNGÜSELLİK BULUNAMADI (mührün lehine):** `KALIB_u2 = Mu2/Pu2` tanımı
   W_X'e bölünmüyor; gaz-düzeyi çarpanlar S = K/K_MID oranında sadeleşiyor;
   tau_eff beş gazda dört haneye kadar özdeş.

## 2.5 Sentezcinin tartısı ve NİHAİ HÜKÜM

**Hakemler aynı hükümde ama farklı yerden vuruyorlar** ve iki bulgu
birbirini **güçlendiriyor**:

- Hakem A: mühür **kendi eşiğini ıskaladı**, sonra ölçüt değişti.
- Hakem B: mührün **ayırt edici sayısı** (λ=1 kimliği) genişletilmiş
  pencerede **rakibe geçiyor**.

Bu ikisi bağımsız yollardan aynı sonuca varıyor: **"kimlik λ=1
gazınındır" ifadesi verinin taşıyabileceğinden fazlasıdır.**

**Sentezcinin kendi doğrulaması:** `171d_A_yuzlesme.py` satır 18-19'u ve
33'ü kaynaktan okudum. Literal metin: *"Q4 Kimlik hükmü: … dokuz bantta
rms ≤ %1.2 tutarsa mühürlenir"* ve ölçüm satırı: *"Q4 ~ A2Hk'nın 9-bant
rms'i %1.49 (öngörü ≤%1.2) — yine de 1-parametreli…"*. **Hakem A'nın
bulgusu kod düzeyinde doğrudur.** Ayrıca `A_ONKAYIT.json`'dan okudum:
`sig_eff = 0.4818447622796531`, `sX_eff = 0.24092238113982656`,
`TAU9 = [0.7386, 0.7774]` (**koda gömülü** ⇒ veri-kör değil),
`sapma(σ_X̃(son)) = +3.0302%`. `S9` alanı bant sayılarını da veriyor:
L115 **6**, Hkeskin **7**, L085 **7**, L070 **6**, L060 **6** ⇒
**"dokuz bant" etiketi 7 noktalı bir rms'i adlandırıyor.**

**Neden ÇÜRÜTÜLDÜ değil?** Hakem B'nin kendi şerhi ciddidir: taban
kenarına 1 bant mesafedeki lo = 0.44/0.48 bantlarında beş gazın paylaştığı
ortak bir alet sistematiği dışlanamıyor — ve bu sınav (bilinen-doğrulu
sentetik kontrol koşusu) hiç yapılmadı. Ortak sistematik gerçekse
düşük-bant reddi geçersizdir. **Bu yüzden kimlik adresi "reddedildi" değil,
"tek bir sınava kadar askıda"dır.**

> ## ⚠ NİHAİ HÜKÜM S2: **ŞÜPHELİ**
>
> **Çürüyen alt-iddialar (bunlar Not 5'e girmemeli):**
> ① *"Genişlik λ=1 gazının σ_X̃'idir, %0.5 içinde"* — ölçülen çözünürlük
> **±%1.7–2.5**, yani −%0.46 sadece **−0.27σ**; rakip σ_X̃(son)
> **+1.2–1.7σ ile ELENMİYOR** ve 7-bant penceresinde **kazanıyor**.
> ② *"1-parametreli en iyi uyumu YENER"* — fark ölçüm hatasının
> **1/9–1/13'ü**; istatistiksel olarak **boş**.
> ③ *"Ön-kayıtlı mühürlendi"* — ön-kayıtlı eşik **rms9 ≤ %1.2 ıskalandı
> (%1.49)** ve mühür sonradan seçilen ölçütle basıldı (kod düzeyinde
> sentezci tarafından doğrulandı: `171d:18-19, 33`).
> ④ İki basılı sayı yanlış: σ\*_eff = 0.48190 (gerçek **0.481845**),
> σ_X̃-eş = 0.24095 (gerçek **0.240922**); tarak-sayımı "×6.0 ARTIŞ"
> (gerçek **×5.224**).
> ⑤ *"9 bant"* etiketi **7 noktalı** bir rms'i adlandırıyor.
>
> **Ayakta kalan (gerçek bulgular):** ① `A(τ)`'nun **λ-değişmez,
> sıfır-parametreli ~Gauss biçimi** — hiçbir hata modelinde kırılmadı.
> ② **W_amp yoktur** (örneklem-içi 10–29σ); **tarak-sayımı işaret
> sınavında ölür** (öngörü ×5.22 artış ↔ ölçüm ×0.797 düşüş);
> **165-F ve λ-değişmez taban örneklem-dışı 3–30σ ile ölür.**
> ③ **Geçerlilik penceresi λ ∈ [0.50, 1.15]** iki yandan doğrulandı
> (rms %3.68 @ 0.40, %3.62 @ 1.30; L040 ve L130 sağlıklı).
> ④ **Genişlik ölçeği:** σ_X̃-eşdeğeri = **0.24092 ± 0.00410**, yani
> **λ_eff = 0.99 ± 0.035** — savunulabilir ifade budur.
>
> **Mührü kapatacak TEK sınav:** Beş gazın **lo = 0.44 ve 0.48**
> bantlarında paylaştığı ortak taban-kenarı sistematiğini ölç
> (taban 0.40 → 0.36 kaydırılıp sapmanın hareket edip etmediğine bak,
> ya da A'sı a priori bilinen sentetik bir kontrol koşusu yap).
> **Sistematik yoksa mühür ÇÜRÜTÜLDÜ'ye düşer; varsa hassasiyet iddiası
> zaten düşer ve mühür "λ_eff = 0.99 ± 0.035" biçiminde SAĞLAM'a çıkar.**

---

# 3. S3 — İKİ ARKIN ÖZDEŞLİĞİ

## 3.1 Mühür cümlesi (iddia)

> `K = ⟨E·e1⟩ ≡ Σ|hp_q|²/2` özdeşliği (**ölçümden önce türetildi**) ⇒
> 162'nin girişim oranı = 172'nin projeksiyon nesnesi ⇒ `g_E ≡ R_η/ρ_E`;
> 8 gazda ≤%1.83 (E) / ≤%4.06 (X); λ_eş(R_X) = 0.9026 ↔ 0.9089.

**Kaynak:** `174_gercek_imza_RAPOR.md` (K2-F, K2-G),
`174_configs/174b_K1_girisim.py` + `174c_K2_yuzlesme.py`.

## 3.2 Yeniden-üretim tablosu

Yeniden-üretici 162 makinesini **kendi ekstraktörüyle** yeniden kurdu
(160.Yerel160 + 154.pk_m birincil; `c_q = 2⟨x e^{−iωm}⟩` **q-bloklu
complex128**, 174b/162'nin n-bloklu cos/sin kurgusundan bağımsız;
8 rastgele çizgide bloksuz doğrudan toplamla sağlandı, sapma ≤ 4e−18).

| İddia | Kaynak değer | Yeniden üretilen | Fark |
|---|---|---|---|
| **MÜHÜR (E):** K = Σ\|hp_q\|²/2, *"hiçbir varsayım olmadan"* | tam özdeşlik | **TAM biçim: K = Σ\|hp\|²/2 − ē·ort(E)**, ē = mean(Y.e1); ölçülen/öngörülen oran = **1.0000000**, bağıl kalıntı 5.1e−7 … 8.6e−7 | ⚠ **özdeşlik değil**; `165_cekirdek.py:222` hp'yi **ortalanmamış** Y.e1'den alıyor, `alanlar()` ise e1'i ortalıyor |
| Aynı mühür, **X kanalı** | — | maks bağıl fark **5.44e−16** | **X'te özdeşlik TAM** |
| R_η(gerçek) | 1.28829 | 1.2882882 | **0.0** |
| R_Ĉ(gerçek) | 1.02599 | 1.0259898 | 0.0 |
| R_X(gerçek) | 1.09022 | 1.0902194 | 0.0 |
| R_η λ merdiveni (7 gaz) | 1.2541 … 1.1739 | aynı yedisi | 0.0 |
| R_X λ merdiveni | 1.1633 … 1.0124 | aynı yedisi | 0.0 |
| **K2-F:** \|R_η/π_E − 1\| maks (8 gaz) | **%1.83** (L130) | %1.83; sekiz satır birebir | 0.0 |
| **K2-F:** \|R_X/π_X − 1\| maks | **%4.06** (L050) | %4.06; aynı sekiz değer | 0.0 |
| **K2-F HÜKMÜ:** *"farkı çizgi kümesi İLE SİTE açıklar"* | iki kaynak | **Çizgi kümesi eşlenince: R_X ≡ π_X^{≤0.86} bağıl fark 0.0e+00…1.0e−15** (bit düzeyinde aynı float); **R_η ≡ π_E^{≤0.86} 1.50e−6…3.03e−6**. Kuyruk payı π_E/π_E^{≤0.86}−1 = %0.87…%1.87 ⇒ **fark sütununu birebir veriyor** | ⚠ **SAYI aynı, YORUM farklı:** sitenin payı ölçülen olarak **SIFIR** |
| **`g_E ≡ R_η/ρ_E`** | özdeşlik (≡) | R_η/ρ_E ↔ g_E: **−%0.86 … −%1.83** (gerçek gazda −%1.19). Gerçek özdeşlik **g_E = π_E/ρ_E**: kalıntı **0.0e+00**, 8/8 gaz | ⚠ **≡ HAK EDİLMEMİŞ** |
| λ_eş(R_X) | 0.9026 | 0.902580 | 0.0 |
| Hedef λ_eş(X) = 0.9089 ⇒ *"TAM İSABET"* | fark −0.0063 | hedef = ort(rX 0.889198, ρ_X 0.928512) = 0.908855. **AYNI nesne τ≤0.95 ile: λ_eş(π_X) = 0.919186** | ⚠ hedefin **iç yayılımı 0.0393**; τ-kesim seçimi λ_eş'i **0.0166** kaydırıyor |
| K2-B: Δlog g_E^ön = +%1.53 | +1.53% | +1.53% | 0.0 |
| K2-D: ΔM kapanışı | %31.4 | %31.4 | 0.0 |
| K1-b: R_η tavanı, gerçek gaz +%1.1 üstünde | 1.2745 @ λ=0.70 | 1.27449; +%1.08 | 0.0 |
| *"Gerçek gaz BÜTÜN ailenin üstünde"* | R_η(son) = 1.28829 > tavan | λ ailesi için doğru. **AMA 174'ün KENDİ K1 dosyalarında R_η(HA4) = 1.29627, R_η(K070) = 1.29543 > 1.28829** (bağımsız doğrulandı: HA4 = 1.29627) | ⚠ İki kesim gazı gerçek gazın **ÜSTÜNDE**; raporun hiçbir tablosunda yok |
| Köprü λ merdiveni **dışında** (174'te sınanmadı) | yok | K090 −%1.35/−%2.52, K070 −%0.87/−%2.13, HA4 −%0.84/−%2.13, E060 −%1.00/−%2.22 ⇒ **12 gazda maks DEĞİŞMİYOR** (%1.83 / %4.06) | Köprü 12 gazda da geçiyor |

**Yeniden-üreticinin kutusu: `uyum = true`** — sayısal uyuşmazlık yok.

## 3.3 Hakem A (istatistik / ön-kayıt merceği) — **ŞÜPHELİ**

1. **ÖN-KAYIT BİÇİMSEL OLARAK SAĞLAM:** `174a_onkayit.py`'nin bugünkü
   sha256'sı dondurulan kayıtla **birebir** (`b643e6c7…`); zaman zinciri
   tutarlı (betik 14:59:12 → mühür 14:59:18 → 174b 15:01:13 → ilk ölçüm
   15:03:33 → K2 15:17 → rapor 15:26); üzerine-yazma reddi kodda var;
   ölen maddeler (K2-A, C, E, D) **kurtarılmadan** rapora geçmiş.
   **Damga sahteciliği veya formül değiştirme izi YOK.**
2. **AMA AYAKTA KALAN İKİ AYAK DÜŞÜK ŞİDDETLİ SINAVLAR:** K2-F'nin manşet
   satırı **ön-kayıt anında hesaplanabilirdi** — hem R_η(son) = 1.288 hem
   π_E(son) = 1.30385 ön-kaydın **kendi "BİLİNEN" listesinde**;
   1.288/1.30385 − 1 = −%1.22 ≪ %5 eşiği. K2-G'nin de "bilinen ikizi"
   vardı: yalnız G1'deki π_X merdivenini dondurulmuş kuralla ters çevirmek
   **λ_eş = 0.9192** verir — hedef 0.9089'un ±0.05 toleransının **zaten
   içinde**. *"TAM İSABET"in gerçekten yeni ampirik içeriği yalnız
   τ-kesim kaymasıdır (−0.0166).*
   **Ölen maddeler (A, C, E, D) ise tam da gerçekten riskli olanlardı.**
3. **DÖNGÜSELLİK DOĞRULANDI (kendi kodu):** 174b'nin c_X'i 165
   önbelleğinin y'siyle **bit-özdeş** (8 gazda maks \|Δ\| ≤ 2e−16).
   V sabit, P eklenen her çizgiyle büyüdüğünden **R ≤ π işareti bir
   TEOREMDİR** — raporun "işareti hep aynı" bulgusu ampirik değil zorunlu.
4. **HEDEF TANIMLARINDA KAYITSIZ SEÇİM:** K2-A hedefi 0.7768,
   `log_172e.txt`'nin kendi bastığı *"E (η) λ_eş ortalaması = 0.7712
   [5 büyüklük]"* ile çelişir — 0.7768, μ̂²_E atılarak alınan **4'lü
   ortalamadır** ve atma kuralı **hiçbir yerde yazılı değil**.
   Hiçbir R, π veya λ_eş için **jackknife/bootstrap σ yok**.
5. **ADİL NOT:** Bütün hedef/ölçüm varyant kombinasyonları denendi,
   **hepsi ±0.05 ölçütünü geçiyor** (maks sapma 0.030) — **GEÇTİ hükmü
   sağlam**, *"TAM İSABET"* dili şişkin.
6. **HÜKÜM (ii) PROGRAMCA DÜZELTİLDİ:** aynı günün 175 raporu (16:20)
   K1_HA4/K1_K070 ile **açıkça yüzleşiyor** (*"gerçek gaz kesim ailesinin
   İÇİNDE"*) — **öz-düzeltme var**, fakat 174'ün retoriği kendi raporunda
   düzeltilmemiş duruyor.

## 3.4 Hakem B (sistematik / alet merceği) — **ŞÜPHELİ**

1. **K2-F BİR SINAV DEĞİL, KESİM ARİTMETİĞİ — kendi hesabıyla teyit:**
   165/tayf önbelleğinden üç satırlık bağımsız hesapta
   **R_X(174b) = π_X^{≤0.86} bağıl fark 0.0e+00** (X, 3 gazda ≤ 3e−16),
   R_η = π_E^{≤0.86} 1.9–3.9e−6.
2. **GEÇME RİSKİ SIFIRDI.** K2-F/K2-G'yi belirleyen **bütün girdiler**
   (8 gazın tayf_\*.npz'si, G1/G3.json, 167/C_\*) ön-kayıttan
   (4 Eyl 14:59:18) **ÖNCE** diskteydi ve sonuç deterministik fonksiyondu.
   Kuyruk payı kasadan hesaplanabilirdi: ≤%5 eşiğinin aşılması için
   162'nin kesimi τ_c ≲ 0.75 olmalıydı; dondurulmuş 0.86'da maks %1.87 —
   **sınavın ölme ihtimali yoktu.**
3. **K2-F EŞİĞİ KODDA YALNIZ E KANALINA UYGULANMIŞ:**
   `174c_K2_yuzlesme.py:110` `mx = max(abs(v["fE"]))` — **%4.06'lık X
   sütunu geçme/ölme kararına hiç girmiyor.** Mühür cümlesindeki
   "≤%4.06 (X)" ön-kayıtlı bir sınavın sonucu değil, **sonradan
   raporlamadır**.
4. **AMA İSABET BOŞ DA DEĞİL (çürütme denemesi başarısız):** ters çevirme
   kuralı duyarlılığı sıfıra yakın (λ-doğrusal/logR/sekant/parabol:
   **0.9021–0.9026**); hedefin **her tek bileşeni** ±0.05'i geçirir;
   τ_c ≥ 0.80 **her kesimde** ±0.05 içinde. R_X(son) ancak
   **−%0.94…+%0.73**'lük pencereye düşseydi geçerdi **ve düştü**.
   X kanalının λ-ailesiyle tutarlılığı **±0.03–0.05 düzeyinde gerçek bir
   bulgudur**.
5. **KONVANSİYON/TABAN SIZINTISI BULUNMADI:** R ve π boyutsuz oranlar;
   iki ark aynı pencere/taban/taper'ı paylaşıyor; yalnız-rX hedefi de
   geçiyor ⇒ kestirimci (g_cal/W_X) sızıntısı **mühre işlemiyor**.
   **Sonuç bir alet artefaktı DEĞİLDİR.**

## 3.5 Sentezcinin tartısı ve NİHAİ HÜKÜM

**Hakemler tam uyumlu** ve ikisi de aynı yerde duruyor: *"aritmetik
kusursuz, kanıtsal ambalaj şişkin."* İkisi de **çürütme aradı, bulamadı** —
uydurulmuş sayı, damga tutarsızlığı, sonradan gevşetilmiş eşik, gizlenmiş
ölüm yok.

**Sentezcinin kendi doğrulaması:** `174_gercek_imza_RAPOR.md`'nin literal
ifadelerini kaynaktan okudum: satır 32 ve 584 *"`g_E ≡ R_η/ρ_E`"*
(**≡ işaretiyle**), satır 607 *"tamamen bağımsız bir ölçüyle (girişim
oranı)"*. Her iki ifade de yeniden-üretici ve iki hakem tarafından
**ölçülerek** yanlışlandı: gerçek özdeşlik `g_E = π_E/ρ_E`'dir (kalıntı 0),
R_η sürümü %0.86–1.83 sistematik taşır; ve R_X, 172'nin π_X'ini üreten
**aynı y_q tayfının kırpılmışıdır** (bit düzeyinde ispatlı).

**Neden ÇÜRÜTÜLDÜ değil?** Çünkü mührün **sayısal** iddialarının hiçbiri
düşmedi ve iki gerçek bulgu ayakta: (a) muhasebe kapanışı **12 gaza
genişletildiğinde de tutuyor** (maks %1.83/%4.06 değişmiyor);
(b) λ_eş(X) tutarlılığı dört ters çevirme kuralına, her hedef bileşenine
ve τ_c ≥ 0.80 kesimlerine **dayanıklı** ve geçiş penceresi dardı
(−%0.94…+%0.73). **Boş bir sınav değil — sadece "TAM İSABET" değil,
"tutarlı".**

> ## ⚠ NİHAİ HÜKÜM S3: **ŞÜPHELİ**
>
> **Çürüyen alt-iddialar:**
> ① *"İki ark"* — **tek arktır**. `c_X` 165'in `y`'siyle **bit-özdeş**;
> K2-F'nin sekiz "fark" değerinin **tamamı** (0.86, 0.95] kuyruk
> aritmetiğidir; sitenin payı ölçülen olarak **sıfır** (E'de 2e−6, X'te 0);
> `R ≤ π` işareti bir **teoremdir**, ampirik bulgu değil.
> ② *"tamamen bağımsız bir ölçüyle doğrulandı"* (HÜKÜM iv) — **yanlış**.
> ③ *"`g_E ≡ R_η/ρ_E`"* — **≡ değil**; gerçek sapma −%0.86…−%1.83.
> Cebirsel köprü **π ile** kurulur (kalıntı 0.0e+00), R ile değil.
> ④ *"Mühür özdeşliği hiçbir varsayım olmadan"* — E kanalında bir
> **düzeltme terimi taşır**: K = Σ\|hp\|²/2 **− ē·ort(E)** (bağıl 5–9e−7;
> sayısal olarak zararsız, retorik olarak yanlış). X'te özdeşlik **tamdır**
> (≤5.4e−16).
> ⑤ *"λ_eş = 0.9026, fark −0.0063, TAM İSABET"* — alıntılanan 0.0063
> hassasiyeti, hedefin **iç yayılımından (0.0393)** ve **τ-kesim
> kaymasından (0.0166)** küçüktür. Doğru ifade:
> **λ_eş(X) = 0.90–0.92 ↔ hedef 0.909 ± 0.02, tutarlı.**
> ⑥ *"Gerçek gaz bütün ailenin üstünde"* — yalnız **λ merdiveni** için
> doğru; 174'ün **kendi** K1 dosyalarında R_η(HA4) = 1.29627 ve
> R_η(K070) = 1.29543 gerçek gazın **üstünde** (175 bunu düzeltti,
> 174 metni düzeltmedi).
> ⑦ K2-F eşiği kodda **yalnız E kanalına** uygulanmış — "≤%4.06 (X)"
> ön-kayıtlı bir sınav sonucu değil.
>
> **Ayakta kalan:** ① **Aritmetik kusursuz** — bütün manşet sayılar iki
> bağımsız kodla birebir; ön-kayıt sha256'sı ve zaman zinciri **gerçek**.
> ② **Muhasebe kapanışı sağlam ve genişletilebilir:** 12 gazda maks fark
> **değişmiyor** (%1.83 E / %4.06 X). ③ **λ_eş(X) = 0.90–0.92**, hedef
> 0.909 ± 0.02 ile tutarlı — **boş olmayan, dayanıklı bir bulgu**.
> ④ K2-A, C, E, D ölümleri **dürüstçe** yazılmış.
>
> **Mührü kapatacak iş:** Mühür cümlesini *"iki ark yapım gereği aynı
> nesnedir (muhasebe kapanışı, keşif değil); X kanalının λ-ailesi
> içindeki yeri ±0.03–0.05 düzeyinde gerçektir"* biçiminde yeniden yaz.
> Hedef λ_eş için **jackknife σ** hesapla ve hedef bileşen atma kuralını
> **yaz**.

---

# 4. S4 — FAZLANIN AYRIŞIMI

## 4.1 Mühür cümlesi (iddia)

> Gerçek gazın E-kanal fazlası = **%39 kesim + %61 kilit** (**üç bağımsız
> kestirici** — f_b, g_E, μ̂²_E — **%2.7 içinde aynı %39**); toplam
> ΔM = **%8.6 kesim + %91.4 kilit**; karıştırma sınavında kilit nedensel
> ölçüldü (R_η 1.265→0.809, **gereğin 76.9 katı**; μ̂²_E −%86) **AMA
> %91.4 atfının nedensel dayanağı 177'de askıya düştü** (H-F1b hükümsüz)
> — resmi statü *"desteksiz, çürütülmemiş"*.

**Kaynak:** `175_kesim_kilit_RAPOR.md`, `176_kilit_faturasi_RAPOR.md`,
`177_ucuncu_tohum_RAPOR.md`; veriler `167/` ve `176/`.

## 4.2 Yeniden-üretim tablosu

Yeniden-üretici **BLOK 2b'yi tam sıfırdan** kurdu:
`128_odl_zeros6_2e6_zeros.npz`'nin son 300000 sıfırı → kendi η zinciri →
kendi merdiveni (asal kuvvetler, τ≤0.95, 8981 çizgi) → kendi
`hp_q`/`κ_q`/`ξ_q`. Sonuç: **η önbelleğiyle bit-bit aynı**
(maks\|Δη\| ≤ 8.9e−16), ξ dizileriyle maks\|Δξ\| ≤ 2.6e−18.

| İddia | Kaynak değer | Yeniden üretilen | Fark |
|---|---|---|---|
| f_b(τ 0.70–0.80) | +0.3975 | +0.3975 (hem 174/K3'ten hem **SIFIRLARDAN**) | **0** |
| f_b(τ 0.80–0.95) | +0.3871 | +0.3871 | 0 |
| iki f_b *"%2.7 farkla aynı"* | %2.7 | %2.69 | Yok — ama bu %2.7 **YALNIZ iki f_b'ye** aittir |
| f(g_E) *"üçüncü kestirici"* | +0.390 | +0.3903 (**log**) — ama 175 §3b tablosunun **kendi konvansiyonu DOĞRUSAL: +0.3828** | ⚠ Konvansiyon seçimi %2.61 ↔ %3.68 yayılım farkı yaratıyor |
| f(μ̂²_E) | +0.397 | +0.3968 | 0.0002 |
| **E kanalı %39 / %61** | %39 / %61 | g_E tek başına %39.0/%61.0 | 0 |
| **ΔM = %8.6 kesim + %91.4 kilit** | +0.004929 / +0.052612 | +0.0049292 (%8.57) / +0.0526124 (%91.43) | toplam kalıntı **9.0e−17** |
| Δlog M özdeşliği | +0.057542 = +0.024547 + 0.001477 + 0.031518 | kalıntı **−9.7e−17** | ⚠ **Bu bir sınav değil**: θ := M/(g_E g_X²) **tanım gereği** |
| R_η: gerçek / ikiz / erfc | 1.28829 / 1.26504 / 1.29627 | 1.288288 / 1.265043 / 1.296267 | **0.0** |
| Karıştırma: R_η 1.265 → 0.809 | 0.80850 | 0.808495 | 0.0 |
| Λ_R, *"gereğin 76.9 katı"* | 76.9 kat | Λ_R = +0.4565467 ⇒ **76.86 kat** | 0.04 — **ama bu oran sınavın GÜCÜNÜ değil ÇÖZÜNÜRLÜK EKSİKLİĞİNİ ölçüyor** |
| μ̂²_E vekilde | −%86 | −%86.1 (n=2), −%86.0 (n=4) | 0.1 puan |
| Q_E / ρ_E / π_E / g_E / M | −68 / −49 / −35 / +26 / +64 | −67.5 / −48.6 / −35.2 / +26.2 / +64.0 | ≤0.5 puan |
| Σξ(τ>0.70): gerçek / keskin / erfc | +2.512385e−02 / +3.782353e−02 / +5.368003e−03 | aynı üçü (**SIFIRLARDAN**) | 0 |
| ω_θ: n=2 / n=3 / n=4 | +0.716 / +0.932 / +1.2235 | aynı üçü; θ̄₄ = 0.8675579, SAÇ_4 = 0.0257336; dal(b) marjı **0.0257759 (kıl payı 4.2e−05)** | 0 |
| 175 §3b: f(θ) | **−0.1430** | G1 defterinde **doğrusal −0.1685, log −0.1495**; −0.1430 ancak **172/G2 defterinin** θ'sıyla çıkıyor | ⚠ **ÇELİŞKİ:** tablonun öteki satırları G1-doğrusal; f(θ) satırı **G2'den** |
| 175 girişi: *"ΔM = +%4.97"* | +%4.97 | Aynı defterde **+%5.92**. +%4.97, 172'nin **λ-EĞRİSİ çapasına** göredir ve `174a_onkayit.py`'de **sabit kodlu** (dM_pct=4.97) | ⚠ **0.95 puan** — iki çapa aynı isimle |
| *"ΔM'nin %51.5'lik θ payı"* (175 §5, 176 §6.3, 177) | %51.5 | Bu raporların **kendi defterinde** 0.031518/0.057542 = **%54.8** | ⚠ **3.3 puan — bayat sayı üç rapor boyunca taşınmış** |
| 177 §5 Var(log θ) payları | 92.4 / 7.5 / 0.1 | **Ham paylar** 85.1 / 6.9 / 0.1; toplam %92.1 (çapraz terimler %7.9). Raporun yüzdeleri bu toplama **normalize edilmiş** | Sonuç (g_X baskın) değişmiyor; **ifade yanlış** |
| *"kesim kesri f — PARAMETRE YOK"* (P5) | parametresiz | **REFERANS GAZ parametredir.** E-kanalı kesim payı: HA4 **%37.1**, K070 **%50.1**, E060 **%20.7**, K090 **−%544** (dejenere). ΔM ayrışımı: +%8.6 / +%14.7 / −%0.3 / **−%264.3** | ⚠ **İDDİANIN EN ZAYIF YERİ** |
| *"üç BAĞIMSIZ kestirici"* | bağımsız | **μ̂²_E ≡ π_E/g_E − Var(E_mod)/Var(e1)** — g_E'nin cebirsel fonksiyonu; kaldıraç **dlnμ̂²/dlng = −53.8** (son), **−134.3** (HA4). g_E %0.1 oynatılırsa f(μ̂²_E) 0.3968→0.4497 (**+%13.3**). f_b ise aynı Σξ'nin **kısmi toplamı** | ⚠ **"bağımsız" nitelemesi desteklenmiyor** |
| DC kaçağının kesim kesri | 0.397 | **AYNI ölçüm, üç okuma:** ort(E) ⇒ **0.4957**; μ̂²_E ⇒ 0.3968; Σξ(τ>0.70) ⇒ 0.3913 | ⚠ **%20 fark** — f yeniden parametrelemeye **değişmez değil**; 0.4957 raporun **kendi §3b tablosunda basılı** |
| %8.6 kesim payı | ölçülmüş nicelik | **İki ters işaretli terimin artığı:** E kanalı **+0.009643**, θ kanalı **−0.004713**. KİLİT_θ/Δθ = **1.1495** ⇒ θ'nın "kilit payı" kendi fazlasının **%115'i** | ⚠ Duyarlılık: θ'ya E'nin f'si verilseydi %38.1; f'ler doğrusal alınsaydı %7.2; ortak α=0.4957 ile %49.6 |
| **LEHTE:** f_b'nin τ>0.70 kararlılığı | 0.3975 ve 0.3871 | 0.05'lik ince bantlarda **0.4287 / 0.3726 / 0.3814 / 0.3851 / 0.3977** (aralık [0.373, 0.429]); 0.85 bölmesinde iki yarı **%0.3 farkla aynı** | **Kararlılık raporun iddia ettiğinden GENİŞ**; olgu bölmeden bağımsız |

**Yeniden-üreticinin kutusu: `uyum = true`.**

## 4.3 Hakem A (istatistik / ön-kayıt merceği) — **ŞÜPHELİ**

1. **KURAL MÜHÜRLÜ, MANŞET ARİTMETİĞİ DEĞİL.** Sha256'lar birebir
   (175a `13a5d195…`, 176a `6173943f…`, 177a `03a9b9ad…`,
   177c `1e5cdd1a…`); eşikler gerçekten dondurulmuş; ölümler (P4, F2, F6,
   F8, H-F1b) **kurtarmasız** yazılmış. **AMA** %8.6/%91.4'ü üreten
   ayrışım formülü **hiçbir ön-kayıtta tutarlı dondurulmamış**:
   175/ONKAYIT'ın P6'sı **başka bir formüldür** ve raporun kendi
   tablosunda −%145.6/−%81.9/+%327.5 verir; 176/ONKAYIT'ın **F9 literal
   metni dejeneredir** (kesim_j = f_j·Δ_j(erfc←Hk) ⇒ **kilit_j ≡ 0**).
   Ayrıca yüzleşme betikleri **ölçümden sonra değişmiştir**
   (175d mtime 16:08 > ilk ölçüm 15:46; 176f 19:26 > 19:05) — sha mührü
   yalnız **onkayit** betiklerini kapsıyor.
2. **HATA ÇUBUĞU YOK — VE JACKKNIFE "AYNI" İFADESİNİ REDDEDİYOR.**
   Delete-1 jackknife (8981 eşleşmiş çizgi): f_b(0.70–0.80) =
   **+0.3975 ± 0.0031**, f_b(0.80–0.95) = **+0.3871 ± 0.0028**;
   fark **+0.0104 ± 0.0042 ⇒ z ≈ 2.5**. Mührün *"%2.7 içinde AYNI"*
   cümlesi, nominal düzeyde **istatistiksel olarak anlamlı** bir farkı
   "aynı" diye sunuyor. Daha ağırı: f(toplam DC) = **0.4957 ± 0.0098** ile
   f(τ>0.70) = 0.3913 ± 0.0021 arasında **>10σ** var.
3. **SEÇİM ETKİSİ + ŞANS TABANI RAPORLANMAMIŞ.** 175 §3b'nin 18 satırlık
   tablosunda **ikinci bir sıkı çakışma** var: f(ort E) = 0.4957 ↔
   f(çizgi payı η) = 0.4845 (**%2.3** — manşetin %2.7'sinden **sıkı**).
   Hangi kestiricilerin "E-kanalı konsensüsü" sayılacağı **hiçbir
   ön-kayıtta tanımlı değil.**
4. **KÜÇÜK-n "KESİN" İLANLARI t-DÜZELTMESİZ:** 176 F3 dal(b) "KESİN
   dışlandı" (n=2, t=2.54, dof=1) gerçekte **tek-yanlı p ≈ 0.12**;
   177 n=3'ün "KESİN" işaret sınavı **p ≈ 0.07** ve n=4'te **kendisi
   düştü**; n=4 dal(b) dışlaması **p ≈ 0.069**. θ-sınavı tasarım gereği
   güçsüzdü: kendi §5a'ları ω_θ ≤ 1 için **n ≈ 82 tohum** gerektiğini
   ölçtü, **tavan 4'tü**.
5. **BOŞ DENETİMLER:** P3'ün "altı gazda tam isabet" satırı 4.1e−13
   düzeyinde tutar — **iki gerçekleştirimin aynı özdeşliği hesapladığının
   kanıtıdır**, fizik öngörüsü değil. D6/D8 denetimleri
   θ := M/(g_E g_X²) tanımının yeniden yazılmasıdır.
6. **LEHTE:** Kilidin **varlığı ve büyüklüğü sağlam** — R_η 1.265→0.8085
   çöküşü tohum saçılımının (0.0028) **~160 katı**; μ̂²_E −%86 dört tohumda
   **%3.3 içinde** aynı. Tayfa **ölümleri saklamıyor**.

## 4.4 Hakem B (sistematik / alet merceği) — **ŞÜPHELİ**

1. **REFERANS GAZ GİZLİ PARAMETREDİR** (kendi bağımsız koşumuyla teyit):
   E-kanalı kesim payı E060 **%20.7**, HA4 **%37.1**, K070 **%50.1**;
   K090 ile **dejenere** (g_E(K090) = 0.5812 < Hkeskin 0.5837, payda işaret
   değiştiriyor). **ADİL NOT:** üçlü 164/172'de **önceden** dondurulmuştu,
   seçim 175'te yapılmadı — bu yüzden **çürütme değil**, "parametresiz"
   nitelemesinin düşmesi.
2. **"%2.7 İÇİNDE AYNI" KONVANSİYON KARMASININ ESERİ:** raporun karması
   log-g_E (0.3903) + doğrusal-μ̂² (0.3968) + iki f_b ⇒ %2.7;
   **tutarlı doğrusal** alınırsa yayılım **%3.8**; **tutarlı log** alınırsa
   **%51** (f(μ̂²)_log = 0.2629). Kümelenme olgusu tutarlı doğrusalda da
   yaşıyor (~%4) — **sahte değil, ama "%2.7" seçilmiş bir montajdır**.
3. **KALIB BANT KONVANSİYONUNA SIZINTI (yeni bulgu):** lo = 0.52 → **%2.0**,
   0.56 → %11.2, **0.60 → %8.6**, 0.64 → %11.4, 0.68 → **%20.3** kesim;
   **lo ≥ 0.72'de Δlogθ işaret değiştirip ayrışım çöküyor (%94.3)**.
   E-kanalı içi %37.1 ise lo'dan **bağımsız**. Yani manşet ayrışımın
   **θ-yarısı bir bant konvansiyonunun eseridir** — üstelik 175'in kendi
   notu W_pos düzeltmesinin 166/167 kalibrasyon zincirinde (θ'nın paydası)
   **yanlış olabileceğini** yazıyor.
4. **"76.9 KAT" LEHTE KANIT DEĞİL:** Λ_R = 0.4565 **tam karıştırma (α=1)**
   bütçesi; gerçeğin ikizden fazlası bütçenin **%5.1'i**, atfedilen kilit
   dilimi **%1.3'ü**. **Ara doz hiç örneklenmedi** (176 §7.5'in kendi
   önerisi bu boşluğu itiraf ediyor).
5. **177'NİN İSTATİSTİĞİ KIRILGAN:** ω_θ tohum başına **0.51–19.93**
   (ortalamayı tek tohum sürüklüyor; medyan 1.78, o da bandın dışında);
   θ tohum varyansının **%92'si g_X böleninden** geliyor ve g_X vekilde
   tek bantlık kestirim (**R_bant filtresi vekili eliyor** — sınamak için
   üretilen nesneyi eleyen alet).
6. **LEHTE:** Bütün manşet sayıları **sıfırlardan bit düzeyinde**
   yeniden üretildi; f_b'nin τ>0.70 kararlılığı **raporun iddiasından
   geniş**; kilit çökertme ölçerleri dört tohumda **%0.8–3.3 içinde**.

## 4.5 Sentezcinin tartısı ve NİHAİ HÜKÜM

**Hakemler uyumlu.** İkisi de aynı üç nitelemeyi hedefliyor
("parametresiz", "bağımsız", "nedensel") ve ikisi de olgunun kendisini
(kilidin varlığı) **sağlam buluyor**.

**Sentezcinin kendi doğrulaması:** Hakem A'nın en ağır bulgusunu
(ön-kayıt formülü) kendim sınadım. `175/ONKAYIT_K2.json` ve
`176/ONKAYIT_K2.json`'u okudum:

- **175-P5 defteri f'yi DOĞRUSAL kaydediyor:** g_E → **0.3828**,
  θ(G1) → −0.1685, θ(G2) → **−0.1430**, ort(E) → **0.4957**.
  Ben aynı defterden log sürümlerini hesapladım: g_E → **0.3903**,
  θ(G1) → **−0.1495**, ort(E) → 0.4122.
  ⇒ **Manşet %39, dondurulmuş defterin doğrusal f'siyle değil, LOG
  f'siyle üretilmiştir**; ve 175 §3b'nin f(θ) = −0.1430 satırı gerçekten
  **G2 defterinden** gelmektedir (birebir tuttu).
- **175-P6 formülü:** `"kesim=f*·ΔlogM(erfc); kilit=(E,X artığı);
  θ=θ artığı"` — uygulanan okuma bu değildir.
- **176-F9 formülü:** `"kesim_j = f_j·Δlog_j(erfc←Hk);
  kilit_j = Δlog_j(son←Hk) − kesim_j"`. Bunu manşeti üreten **log f** ile
  hesapladım: g_E için kesim = +0.0245468 = Δ_j(son←Hk) tam olarak,
  **kilit = +0.000e+00**; θ(G1) için de **kilit = 0**.
  ⇒ **F9'un literal metni dejeneredir — Hakem A haklıdır.**
- **Uygulanan okuma** (kesim_j = f_j^log · Δ_j(son←Hk)) ise
  **%8.57 kesim / %91.43 kilit** veriyor — manşetle birebir.

**Yani manşet ayrışım, dondurulmuş iki formülün hiçbiri değil, üçüncü bir
okumadır.** Bu tek başına bir sahtecilik değil (rapor okumayı açıkça
yazıyor), ama *"ön-kayıtlı ayrışım"* havasını hak etmiyor.

> ## ⚠ NİHAİ HÜKÜM S4: **ŞÜPHELİ**
>
> **Çürüyen alt-iddialar:**
> ① *"kesim kesri f — koordinat yok, ara değer yok, PARAMETRE YOK"* —
> **REFERANS GAZ parametredir**: %20.7 (E060) … %50.1 (K070), K090'la
> dejenere. "%39" gerçek gazın değil, **(son, Hkeskin, HA4) üçlüsünün**
> sayısıdır.
> ② *"üç BAĞIMSIZ kestirici %2.7 içinde aynı"* — üç katmanda çürüyor:
> **cebirsel akrabalık** (μ̂²_E = π_E/g_E − Var oranı, kaldıraç −54…−134;
> f_b aynı toplamın parçası), **konvansiyon karması** (tutarlı doğrusalda
> %3.8, tutarlı logda %51), ve **jackknife** (iki f_b farkı **z ≈ 2.5**,
> yani "aynı" değil "yakın"; f(toplam DC) = 0.4957 ise **>10σ** uzakta).
> ③ *"%8.6 kesim + %91.4 kilit"* — **iki ters işaretli terimin artığıdır**
> (E +0.009643, θ −0.004713; KİLİT_θ = 1.1495·Δθ) ve **KALIB bant
> konvansiyonuna sızar** (lo 0.52 → %2.0, 0.68 → %20.3, ≥0.72'de çöküyor).
> Üreten formül **hiçbir ön-kayıtta tutarlı dondurulmamış** (sentezci
> doğrulaması: 176-F9 literal metni **dejenere**, kilit ≡ 0).
> ④ *"gereğin 76.9 katı"* — sınavın **gücünü değil çözünürlük
> eksikliğini** ölçüyor (atfedilen dilim tam-karıştırma bütçesinin
> **%1.3'ü**; ara doz **hiç örneklenmedi**).
> ⑤ Bayat/karışık sayılar: *"%51.5 θ payı"* (kendi defterinde **%54.8**);
> iki farklı ΔM (+%4.97 λ-eğrisi çapası ↔ +%5.92 Hkeskin çapası) **aynı
> isimle**; f(θ) = −0.1430 satırı **başka defterden** (ayrışımda kullanılan
> −0.1495).
>
> **Ayakta kalan (gerçek bulgular):** ① **Kilidin VARLIĞI nedensel ve
> tohuma duyarsız** — R_η 1.265→0.809 (tohum saçılımının **~160 katı**),
> μ̂²_E **−%86** dört tohumda **%3.3 içinde**, Q_E −67.5, ρ_E −48.6.
> Bu **şüphe götürmez bir olgudur.** ② **f_b'nin τ>0.70 kararlılığı
> raporun iddiasından da geniştir** ([0.373, 0.429] ince bantlarda;
> bölme noktasına duyarsız) — **sabit üçlü koşuluyla** "%39 civarı kararlı
> bir ara-konum" gerçektir. ③ Ön-kayıt **eşikleri** gerçekten dondurulmuş
> ve **ölümler kurtarılmadan** yazılmış. ④ Bütün manşet sayılar
> **sıfırlardan bit düzeyinde** yeniden üretildi.
>
> **Savunulabilir daraltılmış cümle:** *"Sabit (Hkeskin, HA4) çapasına ve
> doğrusal konvansiyona GÖRE E-kanal ara-konumu %39 ± 4'tür; toplam ΔM
> ayrışımı ve %91.4'ün nedensel atfı desteksizdir."*
>
> **Mührü kapatacak iş (eleştirmenin ŞART'ıyla aynı, bkz. §6.4):**
> ara-doz karıştırma (α taraması), referans-aile taraması, jackknife
> σ'ların deftere girmesi, kestirici kümesinin **ön-kayıtlı tanımı**.

---

# 5. S5 — κ ÇİZGİ YASASI

## 5.1 Mühür cümlesi (iddia)

> `κ = −πA_Qτ_Q·cos(πτ_Q)` **SÖNÜMSÜZ** — **7000+ asal çizgide ±%10**
> (τ = 0.60–0.95, medyan 0.965); kulelerin ×2 aşımı seviye denkleminin
> 2. mertebesinden: `A^eff/A = 1+τ` (aşım 1.97→1.07); aynı köşe formu
> 143'ün G-yasası ve 165'in üçüncü-moment tepeleriyle uyumlu.

**Kaynak:** `175_kesim_kilit_RAPOR.md` §K3 (satır 506, 516, 729),
`175_configs/175e_K3_carpimsal.py` + `175g_kappa_yasasi.py`,
`KALEM_UC_DALGA_G_29AGU2026.md`.

## 5.2 Yeniden-üretim tablosu

Yeniden-üretici sıfırları kendi yükledi, mid/L/s zincirini, kendi
asal-kuvvet eleğini, kendi bloklu `κ(ν)=⟨e^{iνs}⟩` toplayıcısını, kendi
W_pos'unu ve kendi 2.-mertebe kanallarını kurdu. **Hiçbir 154–175 modülü
import edilmedi.**

| İddia | Kaynak değer | Yeniden üretilen | Fark |
|---|---|---|---|
| L ve merdiven boyu | L=12.030; 8981 çizgi (8876 asal, 105 kule) | L=12.029593; 8981 (8876+105) | yok |
| ASAL bant medyanları (0.60→0.95) | 1.044, 1.000, 0.965, 0.945, 0.924, **0.899**, 0.998 | 1.0445, 0.9996, 0.9648, 0.9450, 0.9240, **0.8989**, 0.9975 | **yok (4 hane birebir)** |
| τ>0.60 ASAL medyanı ("medyan 0.965") | 0.9648 | 0.9648 | yok — ama bu **ÇİZGİ medyanı değil**, 7 bant medyanının medyanı (**çizgi medyanı 0.9538**) |
| **"τ=0.60-0.95, 7000+ asal çizgide ±%10"** | ±%10 | çizgi sayısı **8658** ✓; ama **çizgi çizgi ±%10 içinde kalan pay %50.9** (4403/8658); std 0.184; %5–%95: **0.600–1.274**. Bant medyanlarından **biri (0.8989 = −%10.1) dışarıda** | ⚠ **İDDİA LAFZEN TUTMUYOR** |
| Bant içi ±%10 yakalama oranı | — | 0.60-0.65 **%98.7** → 0.75-0.80 %81.4 → 0.85-0.90 **%39.9** → 0.90-0.95 **%35.9** | ⚠ Yasanın "sıkı" olduğu yer **τ ≲ 0.75** |
| K3-c ön-mühürlü sınav | gerçek 0.982 / ikiz 1.019 / HA4 4.798 | 0.9823 / 1.0195 / 4.7984 | yok |
| K3-a: med(ASAL, τ>0.60); σ_Ĉ | 2.380; 0.27303 | 2.380; 0.27303 | yok |
| KULE aşımı ve κ̂₂ kapanışı | 1.966 → **1.067** | 1.9658 → 1.0666 | yok |
| KULE bant bant \|κ\|/\|κ̂₂\| | 1.192 … 1.067 | aynı altısı | yok — **0.985–1.192 arası saçılıyor** (0.65-0.70'te **+%19**) |
| **"A^eff/A (ölçülen)" = "öngörü 1+c_kτ", altı bantta birebir** | 1.691 … 1.935 (iki satır **özdeş**) | aynı altısı (iki satır özdeş) | ⚠ **SAYI AYNI, ANLAM FARKLI** — Dfar=0 kulelerde `A^eff/A − (1+H_{k−1}τ)` = **8.9e−16** ⇒ **tanımın cebirsel özdeşliği, ölçüm değil** |
| Çarpımsal kanal katsayısı c_k | k=2 → 1.0; **aksi → 1.5** (175g kodu) | Doğrusu **H_{k−1}**: 1.0, 1.5, **1.8333**, **2.0833**, … ; **105 kulenin 29'u k≥4** | ⚠ **175g'nin c_k tablosu k≥4'te YANLIŞ** |
| *"Kule çizgileri asallara göre tam ×2.03 – 2.33 fazla"* | ×2.03 – 2.33 | med(KULE)/med(ASAL), τ>0.60: 2.390, 2.049, 2.037, **1.946**, 2.077, 2.128, 2.100 ⇒ **1.946 – 2.390** | ⚠ **ARALIK TUTMUYOR (iki uçta da)** |
| *"Kulelerin medyan çokluğu tam 2.0"* | 2.0 | 2.0 | yok |
| §4d — HA4'te asal oranın patlaması | 1.64 → 35.35 | 1.6432 → 35.3544 | yok |
| **GÜRÜLTÜ TABANI** (raporda **yok**) | — | Merdiven-**dışı** ara-nokta frekanslarda: **6.4e−4** (0.75-0.80), **1.02e−3** (0.85-0.90), **1.61e−3** (0.90-0.95); κ̂₀ medyanı 5.53e−3 / 4.14e−3 / 3.40e−3 ⇒ **SNR = 8.6 / 4.0 / 2.1**. Taban düzeltilince: 0.9450→0.9381, 0.8989→**0.8630**, 0.9975→**0.8888** | ⚠ **YENİ BULGU** — tepe banttaki "0.998 toparlanması" **taban kirlenmesi**; düzeltilmiş oran **tekdüze DÜŞÜYOR** (1.034→0.889) |
| *"143'ün G-yasası ve 165'in üçüncü-moment tepeleriyle uyumlu"* | — | **SINANMADI** (kapsam dışı) | ⚠ **DOĞRULANMAMIŞ sayılmalı** |

## 5.3 Hakem A (istatistik / ön-kayıt merceği) — **ŞÜPHELİ**

1. **ÖN-KAYIT MEKANİĞİ SAĞLAM (lehte):** `175a_onkayit.py`'nin sha256'sı
   ONKAYIT_K2.json'daki değerle **birebir** (`13a5d195…`), dosya
   04.09.2026 15:43:46'da **bir kez** yazılmış, tüm ölçüm çıktıları
   (K2 16:08, K3M 16:08, K3G 16:13, K4 16:20) **sonrasında**; eşikler
   betik/JSON/rapor'da özdeş. **Damga oynanmamış.**
2. **AMA "ÖN-MÜHÜRLÜ" SINAVLAR RİSKSİZDİ:** sınanan κ verisi
   (`174/K3_son.json`, mtime 15:10:45) mühürden **33 dakika ÖNCE**
   diskteydi ve 174 raporu bant tablosunu **yayımlamıştı**. 174'ün yayımlı
   değerlerini ön-kayıt betiğinin **kendi Gauss W_pos'uyla** (σ=0.2730,
   `175a:323`) çarpmak yetiyor: K3-c öngörüsü **≈1.014** (ölçülen 0.982),
   K3-a öngörüsü **2.383** (ölçülen 2.380). İkisi de kayıt anında
   **~3 haneye kadar hesaplanabilirdi**.
3. **PENCERE SEÇİMİ SONUCU İZLİYOR:** K3-c'nin τ>0.55 penceresi,
   aynı aritmetikle öngörülebilir değeri **1.563** olan (ölçülen 1.508)
   0.50–0.55 bandını **tam sınırından dışlıyor**.
4. **"SÖNÜMSÜZ" HÜKMÜ PENCERE-KOŞULLU:** aynı gazda τ 0.40–0.45'te
   med \|κ\|/\|κ̂₀\|(asal) = **0.742** — κ̂₀ orada 7.5e−3, tabanın ~5 katı,
   **yani taban açıklamaz** — sönümsüz yasa **−%26 ıskalıyor**;
   174e'nin W_pos'lu yasası ise aynı bantta **0.966** tutuyordu.
5. **HATA ÇUBUĞU EKLENİNCE ±%10 BANT DÜZEYİNDE DE REDDEDİLİYOR:**
   bootstrap SE ≈ 0.002–0.003; 0.85–0.90 bandı med = 0.8989,
   **%95 CI [0.8925, 0.9054] — aralığın TAMAMI ±%10'un dışında** (z ≈ −30).
   Çizgi medyanı 0.9538, CI [0.9508, 0.9568], **1'den 31σ uzak**.
6. **RAPOR TABLOSU SEÇİLMİŞ:** §4c tam da sütunların **özdeş olduğu**
   6 bandı gösteriyor; ham logda (`175/log_175g.txt`) duran ve sütunların
   **gerçekten ayrıştığı** beş alt bant (0.50–0.55: 1.9769↔1.7903;
   0.60–0.65: 1.7813↔1.7954, kule kapanışı **1.3835**) **tablodan
   çıkarılmış** — karşılaştırmanın ampirik içerik taşıdığı tek bantlar
   gösterilmeyenler.
7. **KULE KAPANIŞININ CI'SI 1'İ DIŞLIYOR:** 78 kulede med = **1.047**,
   %95 CI **[1.026, 1.110]**; bant bant 0.985–1.38.
8. **ÇOKLU-DENEME MUHASEBESİ YOK:** aynı 300k sıfır / 8981 çizgi üzerinde
   **en az dört kapalı form** denendi ve nihai yasa **ön-mühürsüz turda**
   seçildi.
9. **LEHTE:** rapor §7a dürüstlük notları **alışılmadık ölçüde açık**
   (175g'nin post-hoc etiketi, Ö-c3 ıskası, P4 ölümünün önceden yazılması);
   olgunun çekirdeği **gerçek**.

## 5.4 Hakem B (sistematik / alet merceği) — **ŞÜPHELİ**

1. **TABAN SIZINTISI DOĞRULANDI:** bağımsız örneklemeyle 0.90-0.95'te
   **1.556e−3** (yeniden-üreticininki 1.614e−3) ⇒ **SNR ≈ 2**;
   taban-düzeltilmiş bant medyanı **0.897**, Hann-pencereli kestirimciyle
   **0.851**. Tepe bandın sayıları **kestirimciye duyarlı**.
2. **AMA VERİ-KOKLAMA İTİRAZI MEKANİZMA İÇİN ÖLDÜ — TAZE PENCERE SINAVI
   (yeni, güçlü kanıt):** 174/175 zincirinin **hiç dokunmadığı** ayrık
   pencerede (sıfır indeksi **1.1M–1.4M**, L = 11.6689, merdiven/Smul/Dfar
   sıfırdan) yasa **aynen tutuyor**: asal bant medyanları
   **1.049 / 0.997 / 0.959 / 0.944 / 0.930 / 0.944 / 0.913**
   (**yedisi de ±%10 içinde**), kule aşımı **×1.9685 → κ̂₂ ile 1.0221**
   (66 kule). **Yasa örneklem-dışında yaşıyor.**
3. **KONVANSİYON SIZINTISI YOK — SİTE SINAVI (yeni, en güçlü kanıt):**
   κ **orta-nokta yerine SIFIRLARIN kendisinde** ölçülüp
   site-dönüştürülmüş kehanetle (κ_z = −πAτ, **cos'suz**)
   karşılaştırıldığında oran **YEDİ BANTTA DA 1.000 ± 0.011**, taban
   **5–10 kat düşük (SNR > 20)**; taze pencerede de 0.993–1.001.
   ⇒ **cos(πτ) çarpanı tam olarak orta-nokta yarım-gap kinematiğidir ve
   temel tepki gerçekten SÖNÜMSÜZDÜR.** Orta-noktadaki taban-düzeltilmiş
   tekdüze düşüş sönüm değil, **siteye özgü 2. mertebe etkidir**.
4. **BİRİM/ÖLÇEK SIZINTISI YOK:** λ = 1.0 tam; L konvansiyonu ±%0.5
   oynatıldığında bant medyanları ≤%3 kayıyor; bant genişliği/kaydırma
   varyantlarında yasa yaşıyor (min 0.889, maks 1.070).
5. **YASANIN PENCERESİ CÜMLEYE TAŞINMAMIŞ:** κ̂₀ ∝ cos(πτ), τ→0.5'te
   sıfırlanır; oran 0.45–0.55'te anlamsız (0.52 / 1.51). *(Sıfır-site
   ölçümünde bu koşulluluk sorunu **yoktur** ve oran yine 1.00'dır —
   sınırlama yasanın değil, **orta-nokta ORANININ** özelliğidir.)*
6. **SAYISAL DÜRÜSTLÜK TAM:** bütün manşet sayılar bağımsız zincirde
   4 haneye birebir; 175g **"ön-mühürsüz" diye dürüstçe etiketli.**

## 5.5 Sentezcinin tartısı ve NİHAİ HÜKÜM

**Bu, beş mührün en ilginç durumu.** Hakemler **hükümde uyumlu** ama
**zıt yönlerde vuruyorlar**:

- **Hakem A** mührün *kanıtsal ambalajını* yıktı: ön-mühürlü sınavlar
  risksizdi, pencere sonucu izliyor, hata çubuğu yok, ±%10 bant düzeyinde
  de reddediliyor, §4c tablosu ayrışan bantları gizliyor.
- **Hakem B** mührün *mekanizmasını güçlendirdi*: taze pencere ve
  sıfır-site sınavları, hakem A'nın veri-koklama itirazını **mekanizma
  için öldürüyor** ve yasanın **gerçekten sönümsüz** olduğunu daha temiz
  bir kestirimciyle (SNR > 20, 1.000 ± 0.011) gösteriyor.

**Bu bir çelişki değil, bir iş bölümüdür:** *fizik sağlam, cümle yanlış.*
Hakem B'nin taze-pencere sınavı (1.1M–1.4M sıfırları) programın hiçbir
seferinde yapılmamış türden bir sınavdır ve **mührün lehine ciddi yeni
kanıttır** — ama o sınav bile bant medyanlarını **±%10'un kenarında**
buluyor (0.913 alt uçta), yani **±%10 rakamını değil, yasanın kendisini**
doğruluyor.

**Sentezcinin kendi doğrulaması:** 175 raporunun literal ifadelerini
kaynaktan okudum: satır 506 *"7000+ asal çizgide **±%10 içinde** tutuyor"*,
satır 729 aynısı, satır 516 *"Kule çizgileri asallara göre tam
**×2.03 – 2.33** fazla veriyor"*. Her iki ifade de üç bağımsız teftişçi
tarafından **ölçülerek** yanlışlandı (çizgi düzeyinde %50.9; gerçek aralık
1.946–2.390).

> ## ⚠ NİHAİ HÜKÜM S5: **ŞÜPHELİ**
>
> **Çürüyen alt-iddialar:**
> ① *"7000+ asal çizgide ±%10"* — **çizgi düzeyinde %50.9** (±%20'de
> %77.2, ±%50'de %97.7); üstelik **bant medyanlarından biri de dışarıda**
> (0.8989 = −%10.1, bootstrap CI'sının **tamamı** ±%10'un dışında).
> Bu bir **bant-medyanı** ifadesidir, çizgi ifadesi değil.
> ② *"medyan 0.965"* — **çizgi medyanı 0.9538**; 0.9648 bant
> medyanlarının medyanıdır.
> ③ *"A^eff/A (ölçülen) = öngörü, altı bantta birebir"* — **cebirsel
> özdeşlik** (Dfar=0 kulelerde fark **8.9e−16**); iki sütun aynı tanımın
> iki yazımı. Rapor tablosu, sütunların **gerçekten ayrıştığı** beş alt
> bandı (kapanışın **1.3835**'e patladığı 0.60–0.65 dahil) ham logdan
> **çıkarmış**.
> ④ *"×2.03 – 2.33"* — gerçek aralık **1.946 – 2.390** (iki uçta da tutmuyor).
> ⑤ `c_k` katsayısı **k≥4'te yanlış** (kodda 1.5; doğrusu **H_{k−1}**:
> k=4 → 1.8333, k=5 → 2.0833). **105 kulenin 29'u k≥4.**
> ⑥ Tepe banttaki *"0.998 toparlanması"* büyük ölçüde **gürültü tabanı
> kirlenmesidir** (SNR ≈ 2; düzeltilmiş 0.86–0.89) — rapor bu mekanizmayı
> yalnız HA4'e (§4d) yazmış, gerçek gazda **hiç sınamamış**.
> ⑦ Kule kapanışı ±%10 değil **±%20** mertebesinde (78 kule; bant bant
> 0.985–1.38; %95 CI [1.026, 1.110] **1'i dışlıyor**).
> ⑧ *"143'ün G-yasası ve 165'in üçüncü-moment tepeleriyle uyumlu"* —
> **bu teftişte sınanmadı; doğrulanmamış sayılmalıdır.**
>
> **Ayakta kalan — ve GÜÇLENEN:** ① **Yasanın kendisi gerçektir.**
> Asal çizgilerde W_pos'lu yasanın **×1.9–3.8 sistematik ıskasına** karşı
> κ̂₀ bant medyanları **0.90–1.05**. ② **Kulelerin ×2 fazlası ikinci
> mertebe çarpımsal kanalla ~1.05'e iniyor** ve `|κ|/|κ̂₂|` sütunu
> **gerçekten ampiriktir**. ③ **TAZE PENCERE:** zincirin hiç görmediği
> 1.1M–1.4M sıfırlarında yasa **ve** kule kapanışı (×1.9685 → **1.0221**)
> örneklem-dışı yeniden üretiliyor ⇒ **veri-koklama itirazı mekanizma
> için ölüyor.** ④ **SIFIR-SİTE:** κ sıfırların kendisinde ölçülünce
> cos'suz kehanetle oran **yedi bantta 1.000 ± 0.011, SNR > 20** ⇒
> **cos(πτ) tam olarak orta-nokta kinematiğidir ve temel tepki gerçekten
> SÖNÜMSÜZDÜR.** ⑤ Birim/ölçek/konvansiyon sızıntısı **yok** (λ = 1.0 tam).
>
> **S5, mühür cümlesi yeniden yazılırsa SAĞLAM'a çıkacak tek mühürdür.**
> Yeniden yazım adresi: **sıfır-site formülasyonu** (orada oran
> 1.000 ± 0.011, cos koşulluluk sorunu yok, taban 5–10 kat düşük) +
> **taze pencere doğrulaması**. Orta-nokta dilinde kalınacaksa:
> *"bant medyanları τ ∈ [0.60, 0.95]'te ±%10–14 (taban-düzeltmeli),
> çizgi payı %51; kule kapanışı ±%20"*.

---

# 6. EKSİKLİK ELEŞTİRMENİNİN BULGULARI

Beş mühür arkın omurgasını kapsıyor, ama manşet katındaki bazı büyük
iddialar **hiçbir mührün içinde değil** — yani bu teftişten **denetimsiz**
çıktılar.

## 6.1 Denetimsiz kalan büyük iddialar

| # | İddia | Kaynak | Neden önemli |
|---|---|---|---|
| E1 | **"4/π² limit değil, λ=1.00 tesadüfü"** — 169'un c₀ = **0.4057 ± 0.0016 = 4/π² + 0.23σ** φ→0 kesişimi hâlâ **açıklamasız** | 170 §K2′, 171 §T2c″; açık borç 170 §6-7 | "Tesadüf" bir **ölçüm değil yorumdur** ve Not 5'in c-anlatısının **temel taşıdır** |
| E2 | **Tek "pozitif taç":** gerçeğin +3.00σ fazlasının **%86/%93'ünün** ölçülen **bacak kırpma aktarımıyla** kapanması | 170 §K3, 171 §T3.3 | Mühür listesi **dışında** ve 174–177'nin kilit anlatısıyla **hiç yüzleştirilmedi** (bkz. T2) |
| E3 | `T(σ) = √(ρ/arcsin ρ)` yasası | 170 §K1 | Yalnız **MC ile** mühürlü; bağımsız denetimi yok |
| E4 | **M(λ) log-λ parabolü** (log M = 0.328(logλ)² − 0.077logλ, rms %0.73) "jilet" olarak kullanılıyor | 173 §H4 | Empirik; **fit penceresi [0.40, 1.30]**, A·M çarpanlaşmasının doğrulanmış penceresi **[0.50, 1.15]**'in **dışına taşıyor**; uç noktalardaki şekil-kirliliği payı **hiç hesaplanmadı** |
| E5 | **λ_sad ≈ 1.38–1.42** "İNŞA SADAKAT SINIRI" ve λ≈1.279 tepesi | 173 §H5 | Tamamen `R_bant`'ın "sadakat" yorumuna dayanıyor; **176 §2c bu yorumu çürüttü** (bkz. T1) — **denetimsiz ve muhtemelen yanlış adlandırılmış** |
| E6 | **β = 0.2175** kesim-seviye yasası θ/θ₀ = 1 − βφ | 168'den devir; 171 §T4 | **Türetimsiz** ama ayrışım anlatısının **girdisi** |
| E7 | 172'nin öteki iki sıfır-parametre isabeti (Q_E tümseği = λ_c; τ-eğimi = −2Δα·τ̄) | 172 Ö5, Ö6/7 | Büyük pozitif iddialar, **bağımsız denetimsiz** *(not: Q_E ↔ λ_c ayağı S1 teftişinde **ayrım gücü yok** çıktı)* |
| E8 | **"X kanalı λ-gazı, η kanalı KESİM-gazı"** tamamlayıcılığı | 175 HÜKÜM (i) | Yalnız **5 noktalı** bir kesim ailesine dayanıyor (175 bunu kendisi söylüyor); 174'ün "aile-dışılık" iddiaları **hata çubuksuz** tavan karşılaştırmalarıyla verildi |
| E9 | **W_pos zincirinin sağlığı** — `W_pos = ⟨e^{2πiτĈ}⟩` çarpanının κ'ya takılmasının **YANLIŞ** olduğu ölçüldü (175 §4b), **aynı çarpan ailesi 166/167 ölçüm zincirinde hesaplanıp taşınıyor** (`166_bacak.py:104-205`, `167_olcum.py:144-270`) | 175 "Sıradaki adım 2" | **KALİB/θ zincirinde nerede doğru nerede yanlış olduğu sınanmadı**; 175 bunun θ'nın ΔM payına dokunabileceğini **açıkça yazdı**, 176–177 **bakmadı**. *(S4 Hakem B'nin lo-konvansiyonu bulgusu bunu bağımsız olarak destekliyor.)* |

## 6.2 Riskli / koşullu ölümler

| # | Ölen madde | Neden riskli |
|---|---|---|
| Ö1 | **176 (iv)** "T-4'ün θ kilide kör savı ÖLDÜ" | Dayanağı ω_θ = 0.716 idi; **177 bunun iki tohumun şansı olduğunu gösterdi** (tohum başına 0.51→19.93) ve dal(b) dışlaması n=4'te **kıl payı** (4.2e−05; 177 §8a kendi raporu "ağırlık taşımayan satır" diyor). Geriye kalan tek kanıt **4/4 işaret tutarlılığı** — ölüm fiilen **işaret sınavına indirgenmiş** |
| Ö2 | **H-F2 ölümü** ve "sahte isabet" etiketi (176 §3g) | Gerekçe τ>0.70 açığının **%61.5'inin GENLİK** olması. Ama aynı raporun (vi) bulgusu **genlik fazlasının KENDİSİNİN kilit ürünü** olduğu — genlik/κ/faz ayrışımı kilit eksenine **dik değilse**, %61.5'in bir kısmı kilittir ve **H-F2'nin 0.61'i rastlantı olmayabilir**. **Ölüm koşullu** |
| Ö3 | **H-K3 ölümü** (175 §4a) | Ölen yalnız *"safça q₁q₂"* **literal biçimdir**; aynı raporun §4c'si çarpımsal kanalın **GERÇEK** olduğunu bulmuştur. Özet satırından "çarpımsal kanal öldü" okunursa **yanlış** — **ölüm etiketi fazla geniş**; diriltilen biçim **ön-kayıtsız** |
| Ö4 | **λ=1.45'in hükümsüzleştirilmesi** (173 §H2, R_bant 0.92–0.96 < 0.98) | 176 §2c aynı kapının **"sınamak için üretilen nesnenin ta kendisini elediğini"** ölçtü. R_bant sadakat değil **kilit-metre** ise, λ=1.45 **çizgileri nominalde, kilidi zayıf SAĞLIKLI bir gaz** olabilir; çöpe atılan hakem verisi (çarpan-yolu 0.4747 ↔ doğrudan 0.4910) **diri olabilir** |
| Ö5 | **172d-Ö4'ün ölümü** ("hiçbir tek değişken θ'yı taşımaz" → ρ_X^0.651 "TAŞIYICI") | **173 §H3.3 taşıyıcıyı iki uçta öldürdü** (−%8.2/+%27.9) ve ıskanın %84'ünün bu model olduğunu ölçtü — **ölen ön-kayıt aslında haklıymış**; "TAŞIYICI" hükmü **geri alınmış sayılmalı** ama survey'de hâlâ öyle duruyor |
| Ö6 | **170 §K2 "taban" okumasının ölümü** (171 §T2c′, ν = 1.489 ± 0.063, **+7.7σ**) | 171'in **kendi dürüstlük notu 4**, ham c(0.50)−c(0.60) farkının yalnız **+1.1σ** olduğunu, +7.7σ'nın **bantlar-arası saçılımın oranda sadeleştiği varsayımına** dayandığını söylüyor. **S1 teftişi bu hata modelini tam olarak çürüttü** — bu ölüm de aynı zeminde duruyor |

## 6.3 Seferler-arası tutarsızlıklar

| # | Tutarsızlık | Durum |
|---|---|---|
| T1 | **`R_bant`'ın iki uzlaşmaz yorumu.** 173 §H5 onu **İNŞA SADAKATİ** sayıp λ_sad ≈ 1.40'ı "yeni sabit" ilan etti ve λ=1.45'i hükümsüz kıldı; **176 §2c/(vi)** `R_bant > 1`'in **KİLİDİN ÜRÜNÜ** olduğunu, faz karıştırılınca 1'e oturduğunu **ölçtü** | **Doğrudan çatışma.** *Sentezci kaynaktan doğruladı:* 173 satır 42-43 *"inşa sadakat sınırı λ_sad ≈ 1.38-1.42"* ↔ 176 satır 562 *"YENİ, ÖLÇÜLMÜŞ BİR OLGU: R_bant fazlası kilidin malıdır"*, satır 567 *"R_bant ≥ 0.98 filtresinin ne ölçtüğünü yeniden tanımlar"*. **Ne 176 ne 177 geri dönüp 173'ü güncelledi** |
| T2 | **Aynı ΔM'nin iki ~%90 açıklaması.** 171 §T3.3: *"ölçülen kırpma aktarımı ΔM'nin **%93**'ünü kapatıyor (+0.7σ artık)"* ↔ 176 K3: *"ΔM'nin **%91.4**'ü KİLİT"* | **İki mekanizma toplamda %184 kapatıyor.** Ya aynı paranın iki para birimi (**gösterilmedi**) ya **çifte sayım**. **Hiçbir sefer ikisini tek defterde yüzleştirmedi.** ⚠ **En sert yapısal sorun** |
| T3 | 172 §G3 *"dar adres λ_eş(E) = 0.777"* ↔ 174 K2-A *"η kanalında λ_eş YOK; 0.777 tekdüze güçlerin eseriymiş"* | 174 **düzeltti**, ama survey'in 172 girdisi hâlâ 0.777'yi *"YÖN ÖN-KAYITLIYDI, TUTTU"* ile taşıyor — **Not 5'e ölü bir sayının kopyalanma riski** |
| T4 | 174 ANA BULGU: *"fazla, sentetik ailenin **HİÇ** taşımadığı bir girişim fazlasıdır"* ↔ 175 HÜKÜM (i): aynı dört nicelik **kesim ailesinin İÇİNDE** | *"Aile-dışılık" aileye göreliymiş.* **S3 teftişi bunu bağımsız olarak pekiştirdi:** 174'ün **kendi** K1 dosyalarında R_η(HA4) = 1.29627 > R_η(son) = 1.28829. İkisi de survey'de **düzeltilmeden yan yana** |
| T5 | 176 HÜKÜM (iv): *"TEK ÖDENEN FATURA: θ, ω_θ = +0.716 ∈ (0,1]"* ↔ 177 §3b: ω_θ = **1.2235**, bandın dışında, **KALICI HÜKÜMSÜZ** | Kaskad **dürüstçe işledi** ama 176 raporunun HÜKÜM bölümü **düzeltme şerhi almadı** — 176'yı tek başına okuyan **yanılır** |
| T6 | 171 §T4: *"HA4 aynı A(τ) üstünde DEĞİL (%4.2)"* ↔ 175'in HA4'ü **erfc-ikiz/kesim çapası** yapması ve %39'u büyük ölçüde τ>0.70 bantlarından ölçmesi | Üstelik 175 §4d aynı gazın en üst bantlarında *"ölçülen şey artık merdivenin tepkisi değildir"* diyor. **Şekil-uyumsuzluğun ve taban-kirliliğinin %39'a katkısı hiç tartışılmadı** |

## 6.4 Eleştirmenin "ŞART OLAN" tek işi

> **ΔM'nin (Q_E, ρ_E) para biriminde nedensel ayrışımı + kırpma/kilit
> çifte-sayım yüzleşmesi.**
>
> Gerçeğin fazlasının (ΔM = +%5.6) **nedensel defterini tek, tohum-duyarsız
> para biriminde** kapatmak: ΔM'yi `g_E/θ` yerine **(Q_E, ρ_E)** çiftinde
> ayrıştırıp (Q_E dört tohumda **%0.8**, μ̂²_E **%3.3** tutarlı — yani
> **ölçülebilir**) "%91.4 kilit" atfını desteksizlikten **kurtarmak YA DA
> öldürmek**; ve **aynı defterde** 170 §K3.2 / 171 §T3.3'ün **%86–93
> kırpma-kapanışıyla çifte-sayım sorusunu yanıtlamak**.
>
> **Gerekçe:** Not 5'in manşet cümlesi kaçınılmaz olarak *"gerçeğin fazlası
> aritmetik faz kilididir"* olacak; 177 §7 bu atfın bugün **"nedensel
> karşılığı olmayan bir defter satırı"** olduğunu kayda geçirdi ve rakip
> mekanizma (kırpma aktarımı) **aynı sayının %93'ünü** kapatıyor.
> **Bu kapanmadan yazılan Not 5, merkez iddiasını desteksiz ve
> iç-çelişkili basar.**

**Sentezcinin değerlendirmesi:** Eleştirmenin ŞART'ına **katılıyorum** ve
**öncelik sırasını onaylıyorum**. Ek olarak: S4 teftişi bu ŞART'ı
**bağımsız olarak doğruladı** — iki hakem de "%91.4"ün ara-doz taraması
olmadan taşınamayacağını ölçtü, ve Hakem B'nin lo-konvansiyonu bulgusu
(ayrışımın θ-yarısı lo 0.52→0.68 arasında %2 ↔ %20 geziyor, ≥0.72'de
çöküyor) eleştirmenin **E9 (W_pos zinciri)** eksikliğiyle **aynı yeri
gösteriyor**. Yani ŞART aslında **iki değil üç ayaklıdır**:
① (Q_E, ρ_E) para biriminde nedensel ayrışım,
② kırpma ↔ kilit çifte-sayım yüzleşmesi,
③ **W_pos'un KALİB/θ zincirindeki doğruluğunun sınanması**.

---

# 7. NOT 5 İÇİN YEŞİL / SARI / KIRMIZI LİSTE

> **Kullanım:** 🟢 **YEŞİL** = olduğu gibi basılabilir (yanındaki
> sayıyla ve kaynağıyla). 🟡 **SARI** = ancak **yazılı hata çubuğuyla** ve
> **daraltılmış cümleyle** basılabilir. 🔴 **KIRMIZI** = **basılmamalı**;
> yanlış, ölü, ya da desteksiz.

## 🟢 YEŞİL — basılabilir

| # | İfade | Sayı ve kaynak |
|---|---|---|
| Y1 | **c(λ) vadisinin VARLIĞI** | c(0.40) > c(0.50) **örneklem-dışı +3.4σ**; sol duvar eşleştirilmiş testte **−7.7σ**, L040→L050 **−9.2σ**; ν(0.50→0.40) = **2.203**, öngörüsel bandın [1.55, 2.45] **içinde**. Minimum **[0.5, 0.85]** civarında. *(179/S1, iki hakem bağımsız)* |
| Y2 | **`dlogW_X/dλ` λ ≤ 1.30'da yaklaşık sabittir** | **−0.42 … −0.63**; üç bağımsız kod birebir. **S1'in tek gerçek ampirik içeriği** *(172 §G4, 179/S1)* |
| Y3 | **A(τ) λ-değişmez ve ~Gauss'tur** | 5-bant rms **%0.62**, 9(=7 nokta)-bant rms **%1.49**; ölçülen W_X ile de rms **%0.61** (Gauss/ölçülen sapması %0.13–0.31) *(171 §T2a, 179/S2)* |
| Y4 | **W_amp YOKTUR** | A1 (W_amp·W_X) örneklem-içi **10–29σ**, 5 bantta ±%9.4, rms9 %10.58 *(171 §T2a, 179/S2)* |
| Y5 | **Tarak-sayımı İŞARET sınavında ölür** | Öngörü **×5.224 ARTIŞ** ↔ ölçüm **×0.7967 DÜŞÜŞ** (τ 0.539→0.698). *(⚠ raporun "×6.0"ı yanlış — **×5.224** bas)* |
| Y6 | **165-F ve λ-değişmez taban örneklem-dışı ölür** | rms9 %13.82 / %10.25 / %9.07; sapmalar 3–30σ *(171 §T2a, 179/S2)* |
| Y7 | **A(τ)'nun geçerlilik penceresi λ ∈ [0.50, 1.15]** | İki yandan doğrulandı: şekil rms **%3.68** @ λ=0.40, **%3.62** @ λ=1.30; pencere içi %0.53–1.50. Dayandığı gazlar (L040, L130) **sağlıklı** *(173 §H4c, 179/S2)* |
| Y8 | **Mühür özdeşliği `K = Σ|hp_q|²/2`** | **X kanalında TAM** (≤**5.44e−16**); E kanalında ortalama-merkezi düzeltme terimiyle tam (bağıl kalıntı **5–9e−7**) *(179/S3)* |
| Y9 | **Muhasebe kapanışı `g_E = π_E/ρ_E`** | Kalıntı **0.0e+00**, 8/8 gaz. *(⚠ `π_E` ile — `R_η` ile **değil**)* |
| Y10 | **Girişim oranı köprüsü 12 gaza genişletildiğinde de tutuyor** | Maks fark **değişmiyor**: %1.83 (E) / %4.06 (X); eklenen dört gaz −%0.84…−%1.35 (E), −%2.13…−%2.52 (X) *(179/S3)* |
| Y11 | **Kilidin VARLIĞI nedenseldir ve tohuma duyarsızdır** | R_η **1.265 → 0.809** (tohum saçılımının **~160 katı**); μ̂²_E **−%86**, dört tohumda **%3.3 içinde**; Q_E −%67.5, ρ_E −%48.6, π_E −%35.2 *(176, 177, 179/S4)* |
| Y12 | **f_b, τ > 0.70'te kararlıdır** | 0.05'lik ince bantlarda **[0.373, 0.429]**; bölme noktası 0.85'te iki yarı **%0.3 farkla** aynı; delete-1 jackknife **±0.003**. **Bölme seçiminden bağımsız** *(179/S4)* |
| Y13 | **κ̂₀ = −πA_Qτ_Q cos(πτ_Q) yasası gerçektir** | Asal bant medyanları τ∈[0.60,0.95]'te **0.90–1.05**; rakip W_pos'lu yasanın aynı bantlardaki sistematik ıskası **×1.9–3.8** *(175 §K3, 179/S5)* |
| Y14 | **Kulelerin ×2 aşımı 2. mertebe çarpımsal kanalla kapanıyor** | **1.9658 → 1.0666** (τ>0.60); `|κ|/|κ̂₂|` sütunu **gerçekten ampiriktir** *(179/S5)* |
| Y15 | **⭐ TAZE PENCERE: yasa örneklem-dışında yaşıyor** | Zincirin **hiç görmediği** 1.1M–1.4M sıfırlarında (L = 11.6689) bant medyanları **1.049 / 0.997 / 0.959 / 0.944 / 0.930 / 0.944 / 0.913** (**yedisi de ±%10 içinde**); kule **×1.9685 → 1.0221** (66 kule) *(179/S5 Hakem B)* |
| Y16 | **⭐ SIFIR-SİTE: temel tepki gerçekten SÖNÜMSÜZDÜR** | κ sıfırların kendisinde ölçülüp cos'suz kehanetle (κ_z = −πAτ) karşılaştırılınca oran **yedi bantta 1.000 ± 0.011**, taban **5–10 kat düşük (SNR > 20)**; taze pencerede 0.993–1.001. ⇒ **cos(πτ) tam olarak orta-nokta yarım-gap kinematiğidir** *(179/S5 Hakem B)* |
| Y17 | **Ön-kayıt disiplini gerçektir** | 174a sha256 `b643e6c7…`, 175a `13a5d195…`, 176a `6173943f…`, 177a `03a9b9ad…`, 177c `1e5cdd1a…` — **hepsi birebir**; zaman zincirleri tutarlı; **ölümler kurtarılmadan yazılmış** *(179/S3, S4, S5)* |
| Y18 | **Aritmetik defterler hanesine kadar doğrudur** | On gazın c/σ_c/Q_E defterleri, bütün eğim tabloları, A(τ) vektörü, R/π defterleri, κ bant medyanları — **birincil veriden sıfırdan yazılan kodla birebir**; η zinciri **bit-bit** aynı (≤8.9e−16) *(179/S1–S5)* |

## 🟡 SARI — yalnız hata çubuğuyla ve daraltılmış cümleyle

| # | Basılabilir daraltılmış biçim | Basılmaması gereken biçim |
|---|---|---|
| S1 | **Vadi konumu: λ\* ≈ 0.65 ± 0.05** (en iyimser hata modelinde ±0.014; konvansiyonlar arası **[0.62, 0.70]**) | ~~λ\* = 0.6487~~ (dört hane = **100–600× aşırı-kesinlik**) |
| S2 | **c(λ)'nın üst durgun noktası λ ≈ 1.28** — beş bandın **beşi de** aşağı-bükümlü parabol veriyor; **anlamlılık ~1–2σ**; MC 95% aralığı **[1.22, 1.62]** | ~~λ = 1.2786, +5.3σ~~ |
| S3 | **A(τ)'nun genişliği: λ_eff = 0.99 ± 0.035** (σ_X̃-eşdeğeri **0.24092 ± 0.00410**) | ~~"λ=1 gazının σ_X̃'i, %0.5 içinde"~~ |
| S4 | **λ_eş(X girişim oranı) = 0.90–0.92 ↔ hedef 0.909 ± 0.02 — tutarlı**; geçiş penceresi R_X biriminde dardı (−%0.94…+%0.73) | ~~"λ_eş = 0.9026, fark −0.0063, TAM İSABET"~~ |
| S5 | **E-kanal ara-konumu: sabit (Hkeskin, HA4) çapasına ve doğrusal konvansiyona GÖRE %39 ± 4** | ~~"%39 kesim + %61 kilit, üç bağımsız kestirici %2.7 içinde"~~ |
| S6 | **ΔM'nin θ payı: %54.8** (bu raporların **kendi** defterinde: 0.031518/0.057542) | ~~%51.5~~ (bayat; λ-eğrisi çapalı eski defterden) |
| S7 | **ΔM (Hkeskin çapası) = +%5.92** — çapayı **her seferinde yaz** | ~~"ΔM = +%4.97"~~ ile ~~"+%5.92"~~ aynı isimle yan yana |
| S8 | **κ yasası (orta-nokta dilinde): bant medyanları ±%10–14 (taban-düzeltmeli), çizgi payı %51; kule kapanışı ±%20** — *veya daha iyisi:* **sıfır-site dilinde 1.000 ± 0.011 (Y16)** | ~~"7000+ asal çizgide ±%10"~~ |
| S9 | **λ_c = 1.9147/1.894 = 1.0109** büyüklüğünün **kendisi** (iki bağımsız inşa kaydıyla tutarlı) — **ama Q_E tepesiyle eşleştirilmeden** | ~~"Q_E tümseğinin tepesi = λ_c"~~ |

## 🔴 KIRMIZI — basılmamalı

**A. Yanlış/çürütülmüş epistemik cümleler**

| # | Basılmaması gereken | Neden |
|---|---|---|
| K1 | *"Vadi bağımsız türetimle dört hanede bulundu"* / *"171'in **bağımsız ölçtüğü** 0.6487"* | **Totoloji.** Sekant kesişmesi = parabol tepesi (**kesin özdeşlik**; sentezci: rasyonel aritmetikte fark **tam 0**). İki sayı aynı üç c değerinin log ve doğrusal okumasıdır *(173 satır 39, 468, 799-800)* |
| K2 | *"fark 0.0001, ‰0.15"* | Gerçek fark **0.000019 / ‰0.03**; rapor **yuvarlanmış hanelerden** çıkarma yapmış |
| K3 | *"alt uca yeni nokta eklenince ‰3 → ‰0.15 keskinleşti"* | **Yanlış adres.** Bozma sınavı: alt kesişim L040'a (ve öteki 6 gaza) **hiç bağlı değil**; keskinleşme **tanım düzeltmesinden** |
| K4 | *"Q_E tepesi = inşa doyum eşiği λ_c"* isabeti | **Ayrım gücü yok.** 36 pencerenin 16'sı λ_c ± 0.05'te (şans tabanı **~%44–50**); kararlı pencere **1.0009** — boş hipotezi (λ=1.0000) **~11 kat** tercih ediyor |
| K5 | *"üst dönüş +5.3σ, sadakat artefaktı değil GERÇEK"* | Hata **yalnız KALİB jackknife'ını** taşıyor; bant saçılımıyla **+1.6σ** (bant-başına ν 0.52–2.53) |
| K6 | *"A(τ)'nun genişliği λ=1 gazının σ_X̃'idir (%0.5 içinde)"* | Çözünürlük **±%1.7–2.5** ⇒ −%0.46 sadece **−0.27σ**; rakip σ_X̃(son) **+1.2–1.7σ ile elenmiyor** ve 7-bant penceresinde **kazanıyor** |
| K7 | *"1-parametreli en iyi uyumu YENER"* | Fark ölçüm hatasının **1/9–1/13'ü**; α'lar yalnız %0.93 ayrık — **istatistiksel olarak boş** |
| K8 | *"ön-kayıtla mühürlendi"* (S2 için) | Ön-kayıtlı eşik **rms9 ≤ %1.2 ıskalandı (%1.49)**; mühür **sonradan seçilen** ölçütle basıldı *(sentezci doğrulaması: `171d:18-19, 33`)* |
| K9 | *"`g_E ≡ R_η/ρ_E`"* (≡ işaretiyle) | **Özdeşlik değil** — sapma **−%0.86…−%1.83**. Gerçek özdeşlik **`g_E = π_E/ρ_E`** (kalıntı 0) |
| K10 | *"tamamen bağımsız bir ölçüyle doğrulandı"* (174 HÜKÜM iv) | `c_X`, 165'in `y`'siyle **bit-özdeş**; `R_X ≡ π_X^{≤0.86}` (fark **0.0e+00**). **İki ark tek arktır** |
| K11 | *"162 ile 172 aynı sayıyı ölçüyor; farkı çizgi kümesi İLE SİTE açıklar"* | **Sitenin payı ölçülen olarak SIFIR** (E'de 2e−6, X'te 0); farkın **tamamı** (0.86, 0.95] kuyruğu. `R ≤ π` işareti bir **teoremdir** |
| K12 | *"kesim kesri f — koordinat yok, ara değer yok, **PARAMETRE YOK**"* | **Referans gaz parametredir**: %20.7 (E060) … %50.1 (K070); K090 ile **dejenere (−%544)** |
| K13 | *"üç BAĞIMSIZ kestirici %2.7 içinde aynı"* | Kestiriciler **cebirsel akraba** (μ̂²_E kaldıracı −54…−134; f_b aynı toplamın parçası); **%2.7 konvansiyon karmasıdır** (tutarlı doğrusalda %3.8, logda %51); jackknife **z ≈ 2.5** |
| K14 | *"ΔM = %8.6 kesim + %91.4 kilit"* — **ölçülmüş nedensel ayrışım** olarak | **İki ters işaretli terimin artığı**; **KALİB bant konvansiyonuna sızıyor** (lo 0.52 → %2.0, 0.68 → %20.3, ≥0.72'de çöküyor); üreten formül **hiçbir ön-kayıtta tutarlı dondurulmamış** (176-F9 literal metni **dejenere**) |
| K15 | *"gereğin 76.9 katı"* — **sınav gücü** kanıtı olarak | Atfedilen dilim tam-karıştırma bütçesinin **%1.3'ü**; **ara doz hiç örneklenmedi** |
| K16 | *"7000+ asal çizgide ±%10"* | Çizgi düzeyinde **%50.9**; bant medyanlarından biri (**0.8989**) de dışarıda (CI'sının **tamamı**) |
| K17 | *"A^eff/A (ölçülen) = öngörü, altı bantta birebir"* | **Cebirsel özdeşlik** (**8.9e−16**); rapor tablosu sütunların **ayrıştığı** beş alt bandı ham logdan **çıkarmış** |

**B. Yanlış basılmış sayılar (doğrusuyla değiştir)**

| # | Basılı | Doğrusu | Kaynak |
|---|---|---|---|
| K18 | σ\*_eff = ~~0.48190~~ | **0.481845** | `A_ONKAYIT.json: 0.4818447623` *(sentezci okudu)* |
| K19 | σ_X̃-eşdeğeri = ~~0.24095~~ | **0.240922** | `A_ONKAYIT.json: 0.2409223811` |
| K20 | Tarak-sayımı ~~"×6.0 ARTIŞ"~~ | **×5.224** | Betiğin **kendi** A4 fonksiyonu; artıklardan çapraz 5.223 |
| K21 | ~~"×2.03 – 2.33"~~ (kule/asal) | **1.946 – 2.390** | *(179/S5)* |
| K22 | `c_k` = ~~1.5 (k≥4)~~ | **H_{k−1}** (k=4 → 1.8333, k=5 → 2.0833) | **105 kulenin 29'u k≥4** |
| K23 | ~~"%51.5 θ payı"~~ | **%54.8** | 0.031518/0.057542, **kendi defterinde** |
| K24 | ~~"9 bant"~~ (A(τ) rms'i) | **7 nokta** (5 örneklem-içi + 2 dışı) | `A_YUZLESME.json/S9`: L115 6, Hk 7, L085 7, L070 6, L060 6 |
| K25 | ~~"9 nokta (0.40–1.45)"~~ (Q_E/α tepesi) | **10 gaz** (173e `slice(None)` L140'ı kapsıyor); gerçek 9-nokta Q_E tepesi **0.9391** | *(179/S1)* |

**C. Ölü sayılar — kopyalanmaya hazır bekliyorlar**

| # | Ölü sayı | Nerede ölmüş | Nerede hâlâ diri duruyor |
|---|---|---|---|
| K26 | **λ_eş(E) = 0.777** *("dar adres")* | **174 K2-A**: η kanalında λ_eş **YOK** | survey 172 girdisi, *"YÖN ÖN-KAYITLIYDI, TUTTU"* ile |
| K27 | **ρ_X^0.651 "TAŞIYICI"** | **173 §H3.3**: iki uçta öldü (−%8.2/+%27.9), ıskanın %84'ü bu model | survey satır ~3295, hâlâ "TAŞIYICI" |
| K28 | **"TEK ÖDENEN FATURA: θ, ω_θ = +0.716"** | **177 §3b**: ω_θ = **1.2235**, **KALICI HÜKÜMSÜZ** | 176 HÜKÜM (iv) bölümü, **şerhsiz** |
| K29 | **"fazla, sentetik ailenin HİÇ taşımadığı"** | **175 HÜKÜM (i)** + 174'ün **kendi** K1 dosyaları (R_η(HA4) = 1.29627 > 1.28829) | 174 ANA BULGU başlık cümlesi, survey'de düzeltmesiz |
| K30 | **λ_sad ≈ 1.38–1.42 "İNŞA SADAKAT SINIRI"** | **176 §2c/(vi)**: `R_bant > 1` **kilidin ürünüdür**, sadakat metresi değil | 173 §H5, survey satır ~3305 |
| K31 | **λ=1.45'in "hükümsüz"lüğü** | **176 §2c**: aynı kapı "sınamak için üretilen nesneyi eliyor" | 173 §H2; çöpe atılan hakem verisi (0.4747 ↔ 0.4910) **diri olabilir** |

**D. Denetimsiz — Not 5'te ancak "denetlenmedi" şerhiyle anılabilir**

| # | İddia | Şerh |
|---|---|---|
| K32 | **"4/π² limit değil, λ=1.00 tesadüfü"** | 169'un c₀ = 0.4057 ± 0.0016 = 4/π² + **0.23σ** kesişimi **açıklamasız**; "tesadüf" bir **yorumdur** |
| K33 | **"Kırpma aktarımı ΔM'nin %86–93'ünü kapatıyor"** | **Kilit anlatısıyla hiç yüzleştirilmedi**; iki mekanizma toplamda **%184** kapatıyor ⇒ **çifte sayım riski** |
| K34 | `T(σ) = √(ρ/arcsin ρ)` | Yalnız MC ile mühürlü |
| K35 | **β = 0.2175** kesim-seviye yasası | **Türetimsiz** |
| K36 | **M(λ) log-λ parabolü** "jilet" olarak | Fit penceresi [0.40, 1.30], doğrulanmış çarpanlaşma penceresi **[0.50, 1.15]** — uçlarda şekil-kirliliği payı **hesaplanmadı** |
| K37 | *"143'ün G-yasası ve 165'in üçüncü-moment tepeleriyle uyumlu"* | **Bu teftişte sınanmadı** |

---

# 8. SENTEZCİNİN KENDİ DOĞRULAMALARI

Hakemlerin en ağır üç bulgusunu, onlardan bağımsız olarak kendim sınadım.

| # | Ne sınadım | Sonuç |
|---|---|---|
| D1 | **S1 totolojisi** — "sekant-eğim kesişmesi = parabol tepesi" gerçekten özdeşlik mi? | `Fraction` ile **kesin rasyonel aritmetikte**, 300 rastgele üçlüde `max\|tepe − kesişme\| = **0**` (tam sıfır; eşit aralık şartı yok). Float ile 2000 üçlüde 5.3e−10 (yalnız polyfit koşullanması). ⇒ **Özdeşlik doğrulandı** |
| D2 | **173'ün "bağımsız" ifadesi** gerçekten yazılı mı? | `173_hakem_RAPOR.md` satır **39**: *"171'in bağımsız ölçtüğü 0.6487"*; satır **468** aynısı; satır **799-800**: *"171'in bağımsız λ\* = 0.6487'si ile fark 0.0001"*. ⇒ **Yazılı ve yanlış** |
| D3 | **S2 ön-kayıt eşiği** ıskalandı mı? | `171d_A_yuzlesme.py` satır 18-19: *"Q4 … dokuz bantta rms ≤ %1.2 tutarsa mühürlenir"*; satır 33: *"Q4 ~ A2Hk'nın 9-bant rms'i %1.49 (öngörü ≤%1.2) — yine de…"*. ⇒ **Kod düzeyinde doğrulandı** |
| D4 | **S2'nin basılı sayıları** | `A_ONKAYIT.json`: `sig_eff = 0.4818447622796531`, `sX_eff = 0.24092238113982656`, `lam_eff = 0.9886089723153767`, `TAU9 = [0.7386, 0.7774]` (**koda gömülü**), σ_X̃(son) sapması **+3.0302%**. ⇒ **Rapor 0.48190/0.24095 basmış — yuvarlama hatası doğrulandı** |
| D5 | **S2'nin "9 bant" etiketi** | `A_YUZLESME.json/S9`: L115 **6**, Hkeskin **7**, L085 **7**, L070 **6**, L060 **6** bant. ⇒ **"9 bant" = 7 nokta; etiket yanıltıcı** |
| D6 | **S4 ön-kayıt formülü** dejenere mi? | `176/ONKAYIT_K2.json/kural/F9/formul`: `kesim_j = f_j·Δlog_j(erfc←Hk); kilit_j = Δlog_j(son←Hk) − kesim_j`. Manşeti üreten **log f** ile hesapladım: g_E → kesim = +0.0245468 = Δ_j(son←Hk) **tam olarak**, **kilit = +0.000e+00**; θ(G1) için de **kilit = 0**. ⇒ **F9 literal metni DEJENERE — doğrulandı** |
| D7 | **S4'ün uygulanan okuması** manşeti veriyor mu? | kesim_j = f_j^log · Δ_j(son←Hk) ⇒ **%8.57 kesim / %91.43 kilit** — **manşetle birebir**. ⇒ Manşet, **dondurulmuş iki formülün hiçbiri değil, üçüncü bir okumadır** |
| D8 | **S4'ün konvansiyon karması** | `175/ONKAYIT_K2.json/kural/P5/defter` **doğrusal** f kaydediyor: g_E **0.3828**, θ(G1) −0.1685, θ(G2) **−0.1430**, ort(E) **0.4957**. Log sürümleri: g_E **0.3903**, θ(G1) **−0.1495**, ort(E) 0.4122. ⇒ **Manşet %39 LOG f'yle; §3b'nin f(θ) = −0.1430 satırı gerçekten G2 defterinden** — doğrulandı |
| D9 | **S5'in mühür cümlesi** gerçekten öyle mi yazılı? | `175_kesim_kilit_RAPOR.md` satır **506** ve **729**: *"7000+ asal çizgide ±%10 içinde tutuyor"*; satır **516**: *"Kule çizgileri asallara göre tam ×2.03 – 2.33 fazla veriyor"*. ⇒ **Yazılı; ikisi de ölçümle yanlışlandı** |
| D10 | **R_bant tutarsızlığı** kaynakta var mı? | 173 satır **42-43**: *"inşa sadakat sınırı λ_sad ≈ 1.38-1.42"*; 176 satır **562**: *"YENİ, ÖLÇÜLMÜŞ BİR OLGU: R_bant fazlası kilidin malıdır"*, satır **567**: *"R_bant ≥ 0.98 filtresinin ne ölçtüğünü yeniden tanımlar"*, satır **289**: *"kilitli gazların nominali AŞTIĞINI"*. ⇒ **Doğrudan çatışma doğrulandı** |

**Sentezcinin betikleri:**
`/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/179/sentez/179_sentez_dogrula.py`
ve `179_sentez_dogrula2.py`.

**Teftiş kuralları:** `qm_riemann` altında **yalnız bu rapor yazıldı**;
başka hiçbir dosya değiştirilmedi. **Git'e dokunulmadı.**
`178_*` dosyalarına ve `scratchpad/178`'e **girilmedi**.

---

# 9. KAPANIŞ — TEFTİŞİN TEK CÜMLESİ

> **Beş mührün aritmetiği kusursuz, dürüstlüğü yüksek, epistemik ambalajı
> şişkindir.** Birincil veriden sıfırdan yazılan kodlarla neredeyse her
> manşet sayı basılan haneye kadar yeniden üretildi; uydurulmuş bir sayı,
> oynanmış bir damga, gizlenmiş bir ölüm bulunmadı — ön-kayıt sha256'ları
> tutuyor ve **ölümler kurtarılmadan yazılmış**. Buna karşılık beş mührün
> **beşinde de** aynı üç hata tekrarlıyor: **(1)** bir **cebirsel özdeşlik**
> ampirik bir buluşma gibi sunulmuş (S1'in kesişim ↔ λ\*'ı, S3'ün "iki
> ark"ı, S5'in `A^eff/A` tablosu, S4'ün D6/D8 denetimleri); **(2)** hiçbir
> **hata çubuğu yazılmamış** ve ölçüm çözünürlüğünün 100–600 katı hassasiyet
> basılmış; **(3)** bir **pencere veya konvansiyon seçimi** olgunun kendisi
> sayılmış (S2'nin LO_MIN = 0.52'si, S4'ün lo bandı, S5'in τ > 0.55'i).
>
> **S1 çürütüldü** — mühür cümlesi bir totolojidir; ama **vadinin varlığı
> gerçektir**. **S2, S3, S4 şüphelidir** ve üçünün de kurtarma yolu ucuz ve
> adreslidir. **S5 şüphelidir ama tersine güçlenmiştir:** hakemin taze
> pencere ve sıfır-site sınavları, cümlesi yeniden yazılırsa **SAĞLAM'a
> çıkacak tek mühür** olduğunu gösteriyor.
>
> **Not 5 yazılmadan önce kapatılması ŞART olan tek iş** beş mührün hiçbiri
> değildir: **ΔM'nin (Q_E, ρ_E) para biriminde nedensel ayrışımı, kırpma
> (%93) ↔ kilit (%91.4) çifte-sayım yüzleşmesi, ve W_pos'un KALİB/θ
> zincirindeki doğruluğunun sınanması.** Bu üçü kapanmadan yazılan Not 5,
> merkez cümlesini **desteksiz ve iç-çelişkili** basar.

---

*179 — Bağımsız Teftiş Tayfası · Yeniden-üreticiler, istatistik hakemleri,
sistematik hakemleri, eksiklik eleştirmeni ve sentezci.*
*Teftişin değeri dürüstlüğündedir — övmek için değil, kırmak için bakıldı.*
