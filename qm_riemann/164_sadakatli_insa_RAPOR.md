# 164 — SENTETİK GAZ İNŞASININ GENLİK-SADAKATLİ ONARIMI ve ×1.4'ÜN YENİDEN YARGILANMASI

**Hedef.** 163 §1c bir yan-kayıt bıraktı: 152'nin sönümlü/kelepçeli
Newton'u `N̄(z) + S(z) = n`'i **eksik gerçekliyor** — gerçekleşen/nominal
merdiven genliği keskin'de **0.54**, A4'te **0.85**, gerçek veride
**1.00–1.04**. Bu, 152'den bu yana kurulan bütün sentetik kıyasları
(×1.4 faz açığı, τ₀ merdiveni, δ parmak izi) etkileyen bir sistematiktir.
Bu rapor inşayı **sadakatle onardı**, açığın **mekanizmasını ölçtü** ve
faz açığını yeniden yargıladı.

---

> ## Kısa hüküm
>
> 1. **İNŞA ONARILDI.** `N̄ + S = n + c` **sıralı ilk-kök** kuralıyla
>    çözüldü: `|F| ≤ 1.863e−09` (= 8 ulp) **HER 300 000 teknede**,
>    sıralılık TAM, ızgara yakınsamış (h = 0.015, 10.4M nokta).
>    Gerçekleşen/nominal merdiven genliği **0.54 → 1.00** (keskin),
>    **0.85 → 1.00** (A4).
> 2. **AÇIĞIN MEKANİZMASI ÖLÇÜLDÜ ve BEKLENENDEN BAŞKA ÇIKTI.**
>    `rms S' = 1.894 ≈ N̄' = 1.915` olduğundan `F` **monoton değil**
>    (ince ızgarada `ΔG < 0` kesri **%19**); `F = n`'in çok kökü var ve
>    152'nin Newton'u seviyelerin **%22'sinde BAŞKA bir köke**
>    yakınsamış. Yakınsama artığı suçlu **değil**: kontrol gazı
>    `NKkeskin` (tam yakınsamış "en yakın kök") genliği 0.538 → **0.553**
>    yapıyor, yani hiçbir şey. Bütün etkiyi **kök seçimi** taşıyor.
> 3. **İKİNCİ VE DAHA DERİN KUSUR BULUNDU: SEVİYE KONVANSİYONU YARIM
>    SEVİYE KAYMIŞ.** Sadakatli gazın çizgileri genlikçe 1.00'i tutuyor
>    ama bir FAZ taşıyor (`arg c_q ≈ −0.18·ω_q`). Seviye taraması fazın
>    **tam olarak `c = −½`'de sıfırlandığını** ölçtü (−0.0003). Kalem
>    karşılığı kesindir: `S = π^{-1} arg ζ` her sıfırda **+1 atlar**,
>    asal-toplam onun simetrik sürümüdür ⇒
>    **`N̄(γ_n) + S_asal(γ_n) = n − ½`**. **152 (ve 151b/153/155'in
>    Newton-tabanlı bütün gazları) `c = 0` kullandı** — `P1`/`P0` gibi
>    gap-örnekleyicili gazlar bu kusurdan etkilenmez.
> 4. **×1.4 FAZ AÇIĞI KAPANDI (senaryo H-İ).** 152'nin ÖZGÜN ölçüsüyle
>    (`arg Γ_sent / arg Γ_gerçek`; taban 0.52, cap 720) bant ortalaması:
>    **152-keskin 1.571 → kök seçimi onarılınca 1.129 → seviye
>    konvansiyonu da onarılınca 0.948.** 152'nin en dayanıklı saydığı
>    τ = 0.66 bandı **1.319 → 1.058 → 0.971**. A4: **1.717 → 1.351 →
>    1.121**. Gerçek gaz, `Hkeskin` (0.948) ile `HA4` (1.121)
>    **ARASINDA**.
> 5. **EŞLEŞEN ŞEY YALNIZ FAZ DEĞİL, TAM KOMPLEKS Γ.** 158 ızgarasının
>    on iki bandında `|Γ|` için gerçek ↔ `Hkeskin`
>    **maks fark 0.0042** (0.922→0.787 ↔ 0.926→0.784). 152-keskin aynı
>    bantlarda 0.81 → 0.80: **τ ile sönmüyor bile.**
> 6. **156'NIN İKİ HÜKMÜ DE DÜŞTÜ.** Bağlaşım `R_sent/R_gerçek`
>    **1.9–4.0 → 0.79–0.97** (Hkeskin), **0.99–1.25** (HA4); rejim
>    metriği `n_eff(τ=0.74)`: gerçek 17 361, 152-keskin **58**,
>    sadakatli **19 337 / 19 798**.
>    "Sentetik uç-değer rejiminde, gerçek yumuşak rejimde" farkı YOK.
> 7. **161'İN PARMAK İZİ MÜHRÜ DÜŞTÜ.** Üç-eksenli (δ(½), dδ/dτ, b)
>    uzaklık (aynı metrik): 152-keskin **16.1**, 152-A4 **47.9**;
>    sadakatli+doğru-konvansiyon **Hkeskin 2.3**, **HA4 1.4** — yani
>    gerçeğin konvansiyon bulutunun **İÇİNDE** (161'in en yakın
>    sentetiği N5z 2.9σ idi). "δ(½)<0 olan tek gaz gerçek" ve "(−,−,−)
>    yalnız gerçeğin" sınıflandırmaları da geçersiz.
> 8. **S3 ÜÇLÜSÜ:** `HA4` σ_η² = **0.0229** (gerçek **0.0227**, %+0.9),
>    c₁ = −0.01041 (gerçek −0.01158), σΔ² = 0.0574 (0.0547);
>    tek belirgin sapma σ_ds² **0.1886 vs 0.1674 (%+13)**.
>    152'nin "iki gözlenebilir ters yönde hareket ediyor, ortak sebeple
>    açıklanamaz" hükmü düştü: **ortak sebep inşaymış.**
> 9. **KALAN AÇIK ve DÜRÜSTLÜK.** Gerçek gaz iki sadakatli gazın
>    ARASINDA: `Hkeskin` δ'yı %8–15 aşıyor, `HA4` %15–25 altında
>    kalıyor; σ_ds² her ikisinde %11–13 yüksek. Yani kalan fark bir
>    "eksik malzeme" değil, **merdiven penceresinin ayarı**.
> 10. **AÇIK BORÇ.** 161'in dokuz-gaz merdiveninin yedi basamağı
>    (J14, J26, N5, N5z, P1, orta, dusuk) bu koşuda YENİDEN İNŞA
>    EDİLMEDİ. 158/161'in sıralama hükümleri (N5z ayrımı, "üçlü
>    gerekli", 2.9σ) **yeniden ölçülmeden kullanılmamalıdır**.

## 1. TEŞHİS — 152 neden eksik kuruyor

### 1a. F MONOTON DEĞİL (ölçüldü)

`F(z) ≡ N̄(z) + S(z) − n`, `S(t) = −Σ_q a_q sin(ω_q t)`, `a_q = 1/(πm√q)`.
Bu pencerede (son 300k, `L = 12.0304`, τ ≤ 1.00 ⇒ **15450 çizgi**):

| nicelik | keskin merdiveni | A4 merdiveni (erfc 0.68/0.125) |
|---|---|---|
| `N̄' = log(t/2π)/2π` | **1.9147** | 1.9147 |
| `rms S' = √(½Σ(a_qω_q)²)` | **1.8940** ← N̄' ile aynı mertebe | 1.1844 |
| `max\|S'\| ≤ Σ a_q ω_q` | 259.9 | 43.0 |
| `rms S` | 0.3825 | 0.3503 |
| ince ızgarada `ΔG < 0` kesri | **0.1915** | **0.0673** |

> **`F' = N̄' + S'` işaret değiştiriyor: ince ızgarada `G = N̄+S`'nin
> ardışık farkı keskin'de noktaların %19.2'sinde NEGATİF.** Yani
> `F(z) = n` denkleminin ÇOK KÖKÜ vardır ve "kök bul" tanımı tek başına
> yetmez.

### 1b. 152'nin Newton'u — ve neden yakınsamıyor

```
payda = max(N̄' + S', 0.3·N̄')          ← negatif türevi 0.3·N̄'e KIRPAR
adım  = clip(0.8·F/payda, ±ḡ)          ← sönüm 0.8 + ±ḡ kelepçesi
20 iterasyon; |F|max < 1e-3 olunca dur
```

`164_insa.py eski_keskin` bu Newton'u satır satır yeniden koştu ve
**152'nin kayıtlı `z_keskin.npy`'sini `maks|Δz| = 0.0e+00` ile bit
düzeyinde yeniden üretti** (aynısı `eski_A4` ↔ `z_A4.npy`: **0.0e+00**).
Karşılaştırmanın meşruiyeti böylece kurulmuştur (§7 V1).

Durma ölçütü **hiç sağlanmadı**: 20 iterasyonun sonunda

| gaz | medyan\|F\| | %99\|F\| | **maks\|F\|** | \|F\|>1e−8 tekne |
|---|---|---|---|---|
| 152-keskin | 2.33e−10 | 6.15e−03 | **1.300** | 8795 (%2.93) |
| 152-A4 | 2.33e−10 | 3.55e−05 | **7.07e−02** | 11134 (%3.71) |

### 1c. AMA ASIL SEBEP YAKINSAMA DEĞİL — **KÖK SEÇİMİ**

Özdeşlik (`164_tani_artik.py`, kapanış ≤ 9.9e−10):

```
ds_n ≡ g_n·log(m_n/2π)/2π − 1 = −ΔS_n + ΔF_n + O(g³N̄‴)
                                              (N̄‴ ≈ 2.6e−13 ⇒ ≤ 1e−14)
ΔF_n = F_{n+1} − F_n      ← ÇÖZÜM ARTIĞININ ds'e sızması
```

Artık kanalı ölçüldü ve **küçük** çıktı:

| gaz | Var(ds) | Var(−ΔS) | **Var(ΔF)** | kor(−ΔS,ΔF) | `c_2(ΔF)/b_2` |
|---|---|---|---|---|---|
| 152-keskin | 0.12795 | 0.14406 | **0.02034** | −0.337 | −0.057 |
| 152-A4 | 0.11284 | 0.11285 | **0.0000** | −0.023 | −0.000 |
| **164-Skeskin** | 0.16962 | 0.16962 | **0.0000** | +0.031 | +0.000 |
| **164-SA4** | 0.12733 | 0.12733 | **0.0000** | +0.038 | +0.000 |

152-A4'te artık kanalı **tamamen ölü** olmasına rağmen genlik oranı yine
0.846. Ve 152-keskin'de yalnız `|F| ≤ 1e−8` olan **iyi** teknelere
kısıtlansa bile oran **0.523** (tümü: 0.538). Yani eksiklik "yakınsamamış
tekneler"den GELMİYOR.

**Ölçülen sebep: 152'nin Newton'u BAŞKA BİR KÖKE yakınsıyor.**
Aynı `n` indeksinde iki çözüm karşılaştırıldığında:

| çift | farklı tekne | ort Δz | sd Δz | maks \|Δz\| |
|---|---|---|---|---|
| 152-keskin ↔ 164-Skeskin | **%22.06** | +0.0837 (+0.160 ḡ) | 0.1745 (0.334 ḡ) | 1.027 (1.97 ḡ) |
| 152-A4 ↔ 164-SA4 | **%5.74** | +0.0263 | 0.1142 | 1.053 |

Yanlış kök de `F = n`'i sağlar (bu yüzden `|F|` küçüktür!) ama merdivenin
FAZINA KİLİTLİ bir konumdadır; seviyelerin beşte birinde yanlış kök almak
çizgi genliğinin **%46'sını** siliyor. A4'te oran %5.7 yanlış kök ↔ %15
genlik kaybı. **Genlik açığı bir yakınsama kusuru değil, bir KÖK SEÇİMİ
kusurudur.**

---

## 2. SADAKATLİ İNŞA — SIRALI İLK-KÖK

### 2a. Tanım

```
z_n = F(z) = n denkleminin, alt sınırdan itibaren İLK kökü
```

Bu tanım (i) **sıralılığı özdeş olarak garanti eder** (seviye n+1 > n ⇒
ilk kökü de sağdadır), (ii) **sayımı korur** (z_n'in altında tam n tekne),
(iii) F'nin monoton olduğu her yerde tek köke indirgenir.

### 2b. Algoritma (`164_configs/164_insa.py`)

1. **İnce ızgara.** `S` ızgarada BİR KEZ, TAM (yaklaşıksız), 7 süreçli
   havuzla hesaplanır — `10 449 280` nokta, `h = 0.015`.
2. `G_k = N̄(z_k) + S(z_k)`; **koşan maksimum** `M_k = max(G_0..G_k)`
   monotondur ve `M_k ≥ n`'i sağlayan İLK k'da `G_k = M_k`'dir. Yani
   `k_n = searchsorted(M, n)` **ilk ızgara-geçişini** verir;
   braket `[z_{k−1}, z_k]`, `G_{k−1} < n ≤ G_k`.
3. Braket içinde **korumalı Newton**: doğrusal ara değerden başlanır;
   her adımda TAM `S, S'`; Newton adımı braketten çıkarsa ya da `F' ≤ 0`
   ise ikiye bölme. Yakınsama **garantilidir**. Yalnız yakınsamamış
   noktalar üzerinde çalışılır.
4. **Kübik interpolasyon kullanılmadı** ve nedeni ölçüldü:
   `rms S⁗ = 1665` ⇒ kübik spline hatası `h⁴·1665/384`, `h = 0.0522`'de
   `3.2e−05`. Yani interpolant zaten `|F| ≤ 1e−8` ölçütünü sağlayamaz;
   **ölçüt İNTERPOLANTTA değil GERÇEK F'de** sağlanmalıdır. Bu yüzden
   ızgara yalnız **braket** için kullanıldı, kök TAM `S` ile bulundu.

### 2c. Izgara adımı ÖLÇÜLDÜ (görevin 1/10-dalga-boyu ölçütü YETMİYOR)

En kısa dalga boyu `2π/ω_max = 0.5223` ⇒ görevin ölçütü `h ≤ 0.0522`.
5000-tekne smoke koşusunda ızgara-yakınsama sınavı:

| h | `h = 0.00375`'e göre farklı tekne | maks \|Δz\| |
|---|---|---|
| 0.0400 (≈1/13 dalga boyu) | **4 / 5000** | 0.276 |
| **0.0150** (≈1/35) | **0 / 5000** | 0.0 |
| 0.0075 (≈1/70) | 0 / 5000 | 0.0 |

> **1/10 dalga boyu ölçütü tek başına yetmiyor**: alt-ızgara çift
> geçişleri (F'nin bir hücre içinde n'i aşıp geri düşmesi) kaçırılıyor.
> Üretim koşuları `h = 0.015` ile yapıldı — **ızgara yakınsamış**.

### 2d. Yakınsama — HER teknede

| gaz | ızgara | ΔG<0 | maks\|F\| | **\|F\|>1e−8 tekne** | ikiye-bölme adımı | sıralı | min Δz |
|---|---|---|---|---|---|---|---|
| **Skeskin** | 10 449 280 | 0.1915 | **1.863e−09** | **0** | 46 | TAM | 0.1045 (0.200 ḡ) |
| **SA4** | 10 447 842 | 0.0673 | **1.863e−09** | **0** | 10 | TAM | 0.1422 (0.272 ḡ) |

**Kayan-nokta tabanı raporlanıyor:** `F = (N̄(z) − n) + S(z)` ve
`N̄ ≈ 1.7e6` olduğundan çıkarma iptali `|F|`'yi `ulp(n) = 2.33e−10`'un
altına indiremez. Ölçülen maksimum **1.86e−09 = 8 ulp**; görev ölçütü
`1e−08` bu tabanın **4 katı** üstündedir. Yani ölçüt SAĞLANDI ve
"daha iyisi çift-duyarlıkta mümkün değildir".

---

## 3. SADAKAT DENETİMİ (bant bant)

Yöntem 163 §1c'nin aynısı (`164_configs/164_sadakat.py`):
`c_q = 2⟨(ds−⟨ds⟩)e^{−iω_q m}⟩` ölçülür, `b_q = 2a_q^{(pencereli)}sin(πτ_q)`
ile karşılaştırılır. **163'ün çizgileri bit düzeyinde yeniden üretildi**
(152-keskin: 0.538/0.540/0.556/0.547/0.544/0.571/0.491 — 163 §1c ile aynı).

Bant oranı **gürültü çıkarımlıdır**: ölçüm zincirinin kendi ara-nokta
referansı (`W' = w + gap/2`, aynı `2.5·dres` eleği) her çizgi için
ölçülür ve `R = √((Σ|c_on|² − Σ|c_off|²)/Σb²)`.

### 3a. Çizgi bazında (163 §1c tablosunun devamı)

| q | τ | **gerçek (son)** | 152-keskin | **164-Skeskin** | 152-A4 | **164-SA4** |
|---|---|---|---|---|---|---|
| 2 | 0.0576 | 1.004 | 0.538 | **1.003** | 0.846 | **0.997** |
| 3 | 0.0913 | 1.008 | 0.540 | **1.000** | 0.852 | **0.991** |
| 5 | 0.1338 | 1.017 | 0.556 | **1.011** | 0.843 | **0.979** |
| 7 | 0.1618 | 1.024 | 0.547 | **1.008** | 0.832 | **0.968** |
| 11 | 0.1993 | 1.035 | 0.544 | **1.010** | 0.824 | **0.957** |
| 101 | 0.3836 | 1.116 | 0.571 | **1.033** | 0.732 | **0.833** |
| 1009 | 0.5750 | 1.215 | 0.491 | **0.979** | 0.491 | **0.614** |

### 3b. Bant bazında (gürültü çıkarımlı)

| τ-bant | gerçek (son) | 152-keskin | NK-keskin | **164-Skeskin** | 152-A4 | **164-SA4** |
|---|---|---|---|---|---|---|
| 0.05–0.10 | 1.007 | 0.539 | 0.552 | **1.001** | 0.850 | 0.993 |
| 0.10–0.15 | 1.017 | 0.532 | 0.540 | **1.012** | 0.834 | 0.980 |
| 0.15–0.20 | 1.031 | 0.533 | 0.534 | **1.016** | 0.812 | 0.963 |
| 0.20–0.25 | 1.047 | 0.548 | 0.547 | **1.021** | 0.803 | 0.938 |
| 0.25–0.30 | 1.067 | 0.547 | 0.543 | **1.031** | 0.762 | 0.915 |
| 0.30–0.35 | 1.087 | 0.554 | 0.548 | **1.030** | 0.744 | 0.880 |
| 0.35–0.40 | 1.113 | 0.552 | 0.548 | **1.033** | 0.721 | 0.844 |
| 0.40–0.45 | 1.139 | 0.546 | 0.545 | **1.036** | 0.694 | 0.808 |
| 0.45–0.50 | 1.165 | 0.536 | 0.536 | **1.032** | 0.673 | 0.761 |
| 0.50–0.55 | 1.190 | 0.521 | 0.522 | **1.023** | 0.613 | 0.695 |
| 0.55–0.60 | 1.213 | 0.488 | 0.490 | **0.999** | 0.472 | 0.572 |
| 0.60–0.65 | 1.235 | 0.460 | 0.465 | **0.974** | 0.112 ⚠ | ⚠ |
| 0.65–0.70 | 1.249 | 0.418 | 0.421 | **0.948** | ⚠ | ⚠ |
| 0.70–0.75 | 1.259 | 0.370 | 0.371 | 0.904 | ⚠ | ⚠ |
| 0.75–0.80 | 1.245 | 0.305 ⚠ | 0.302 ⚠ | 0.845 | ⚠ | ⚠ |
| 0.80–0.85 | 1.203 | 0.050 ⚠ | 0.024 ⚠ | 0.738 | ⚠ | ⚠ |

(⚠ = SNR = Σ|c_on|²/Σ|c_off|² < 3, yani bant ölçüm gürültüsünün içinde;
A4'ün τ > 0.60 bantlarında erfc penceresi nominali sıfıra sürdüğü için
oran **tanımsızdır** — 0.112 / 2.27 / 6.71 / 20.3 gibi sayılar bir
gürültü/sıfır bölümüdür, hükme girmez.)

> **SADAKAT HÜKMÜ.** Görevin ölçütü "τ-bantlarında ≥ 0.98".
> `Skeskin` (c = 0) **τ ≤ 0.60'ta 0.999 – 1.036** ile ölçütü sağlıyor
> (`0.60–0.65`'te 0.974, ölçütün %0.6 altı); τ > 0.65'te 0.95 → 0.74'e
> iniyor. Eski inşa aynı bantlarda **0.49 – 0.56**.
> **Doğru seviye konvansiyonunda (`Hkeskin`, §3g) ölçüt τ ≤ 0.90'ın
> BÜTÜN bantlarında sağlanıyor** (1.010 … 1.397) — ama orada da
> yüksek τ'da gerçeği %11 aşıyor. **Genlik açığı KAPANDI**; yüksek-τ'da
> kalan fark artık bir çözüm kusuru değil, ölçünün rezonant beslemesi
> ve gaz farkıdır (§3c).

### 3c. Gerçek gazda oran 1.00'de KALMIYOR — ve bu yöntemin sınırı

Gerçek veri kontrolü, oranın düşük τ'da **1.007** (yöntemin doğrulaması)
ama yüksek τ'da **1.26**'ya çıktığını gösteriyor. Bu bir kusur değil,
ölçünün kendisinin bir özelliğidir: `c_q`, kendi çizgisinin yanında
`ω_q`'ya rezonant ÇİFTLERDEN de beslenir (163 §2c'nin çarpımsal
üçlüleri) ve o besleme gazın kendi `ds` istatistiğine bağlıdır.
`⟨sin(ω_q g/2)⟩/sin(πτ_q)` doğrudan ölçüldüğünde beş gazda da
**0.997 (q=2) / 0.96–0.98 (q=11) / 0.77–0.84 (q=1009)** çıkıyor —
yani "kendi çizgi" terimi bütün gazlarda aynı; oranların farkı
**tamamen rezonant beslemeden** geliyor. **Bu yüzden 0.98 ölçütü ancak
düşük-τ'da (τ ≲ 0.6) bir SADAKAT sınavıdır; yüksek τ'da bir gaz-farkı
ölçüsüdür.**

### 3d. YENİ VE BEKLENMEDİK: sadakatli gaz çizgilerinde bir FAZ var

163 §1c "bütün çizgilerin fazı `arg ≈ 0` (maks 0.05 rad): belirlenimli
merdiven varsayımı üç gazda da doğrudan doğrulanmıştır" demişti. Sadakatli
inşada bu **bozuluyor**:

| q | τ | **son** | 152-keskin | NK-keskin | **164-Skeskin** | 152-A4 | **164-SA4** |
|---|---|---|---|---|---|---|---|
| 2 | 0.0576 | +0.0000 | +0.0020 | +0.0024 | **−0.1282** | +0.0012 | **−0.0872** |
| 3 | 0.0913 | −0.0000 | −0.0039 | −0.0042 | **−0.2087** | −0.0009 | **−0.1337** |
| 5 | 0.1338 | +0.0000 | +0.0023 | +0.0019 | **−0.2826** | −0.0043 | **−0.1876** |
| 7 | 0.1618 | −0.0000 | +0.0060 | +0.0047 | **−0.3457** | +0.0022 | **−0.2196** |
| 11 | 0.1993 | +0.0000 | +0.0053 | +0.0038 | **−0.4138** | +0.0066 | **−0.2666** |
| 101 | 0.3836 | −0.0004 | +0.0161 | +0.0216 | **−0.7038** | +0.0078 | **−0.4647** |
| 1009 | 0.5750 | −0.0083 | −0.0486 | −0.0616 | **−0.9758** | −0.0358 | **−0.8675** |

`arg c_q / ω_q` neredeyse sabittir: Skeskin **−0.185 … −0.141**,
SA4 **−0.126 … −0.111**. `ds ⊃ A·cos(ω_q(m − Δ))` demek olduğundan bu,
gerçekleşen yapının merdivene göre **Δ ≈ 0.18 (0.34 ḡ)** / **0.12
(0.23 ḡ)** kaymış olması demektir.

> **İlk okuma (ve neden yanlış çıktı):** faz ilk bakışta İLK-KÖK
> kuralının "yükselen kenar" yanlılığı gibi duruyor — `NKkeskin`
> (en yakın kök) fazsız (≤ 0.02 rad) ama genliksiz (0.55); `Skeskin`
> (ilk kök) genlikli (1.00) ama fazlı (−0.18ω). §3e bunu **sınadı ve
> çürüttü**: faz kök seçiminin değil **SEVİYE KONVANSİYONUNUN**
> işidir ve `n → n − ½`'de tam olarak sıfırlanır.


### 3e. FAZIN KAYNAĞI ÇÖZÜLDÜ: SEVİYE KONVANSİYONU `n → n − ½`

§3d'nin fazı, kök seçiminin değil **seviye konvansiyonunun** işidir ve
bu ölçülerek gösterildi (`164_tani_seviye.py`, 40 000 tekne; ızgara
`c`'den bağımsız olduğu için tek pasoyla yedi `c` değeri):

```
denklem:  N̄(z) + S(z) = n + c
```

| c | maks\|F\| | σ_ds² | \|c_q\|/b (q = 2,3,5,7,11,101) | **arg c_q / ω_q** |
|---|---|---|---|---|
| −1.00 | 1.86e−09 | 0.1715 | 1.011 1.005 1.016 1.016 1.020 1.047 | **−0.1735** |
| −0.75 | 1.86e−09 | 0.1811 | 1.001 1.012 1.022 1.035 1.047 1.145 | +0.0702 |
| **−0.50** | 1.86e−09 | 0.1847 | **1.000 1.009 1.020 1.032 1.048 1.151** | **−0.0003** |
| −0.25 | 1.86e−09 | 0.1805 | 1.000 1.009 1.017 1.030 1.045 1.141 | −0.0719 |
| **0.00** | 1.86e−09 | 0.1714 | 1.011 1.005 1.016 1.016 1.019 1.048 | **−0.1735** |
| *(gerçek gaz referansı)* | | 0.1674 | *1.004 1.008 1.017 1.024 1.035 1.116* | *≈ 0.000* |

`c` ve `c+1` **tam olarak bir indis ötelemesidir** (`z_n(c+1) = z_{n+1}(c)`),
ve ölçüm bunu doğruluyor: `c = 0` ile `c = −1` satırları özdeş
(σ_ds² 0.1714 / 0.1715, arg/ω −0.1735 / −0.1735). Yani `c` fiziksel
olarak yalnız **kesirli kısmıyla** anlamlıdır ve fazı sıfırlayan değer
**tam olarak `c = −½`**'dir.

> **KALEM KARŞILIĞI (ve bu yüzden bir tesadüf değil).** Riemann–von
> Mangoldt: `N(T) = N̄(T) + S(T)`, `S = π^{-1} arg ζ(½+iT)` her sıfırda
> **+1 atlar**. Standart (Titchmarsh) konvansiyonunda sıfırın kendisinde
> `N(γ_n) = ½[N(γ_n^+) + N(γ_n^-)] = n − ½`, ve asal-toplam tam olarak
> `S`'nin bu **simetrik/düzgün** sürümüdür. Dolayısıyla
> **`N̄(γ_n) + S_asal(γ_n) = n − ½`** — yani doğru denklem `c = −½`'dir.
> **152 (ve bu raporun §2 gazları) `c = 0` kullandı: yarım seviye
> kaymış.** Ölçülen faz (`−0.18·ω`, yani `Δ ≈ 0.34 ḡ`) bunun imzasıdır
> ve `c = −½`'de **0.000**'a düşüyor.
>
> Üstelik `c = −½`'de genlik profili de gerçeğin profiline oturuyor
> (τ ile YÜKSELEN 1.00 → 1.15; `c = 0` düz 1.01 → 1.05, gerçek
> 1.00 → 1.12). **Yani doğru konvansiyonda sadakatli gaz merdivenin
> hem genliğini hem fazını gerçek gibi kuruyor.**

Bu bulgu üzerine **`Hkeskin`** (`c = −½`, keskin) ve **`HA4`**
(`c = −½`, erfc 0.68/0.125) gazları da tam ölçekte (300 000 tekne)
inşa edildi. Rapor **iki konvansiyonu da ayrı ayrı taşıyor**:
`Skeskin`/`SA4` (`c = 0`) satırları **152'nin kendi konvansiyonundadır**
— ×1.4'ün yeniden yargılanması orada yapılmalıydı (§4a) —, `Hkeskin`/`HA4`
(`c = −½`) satırları **doğru konvansiyondadır** (§3g, §4d–4f, §5).

### 3f. GAP DAĞILIMI — kısa-mesafe itmesi

| gaz | σ_ds² | **P(s<0.3)** | P(s<0.5) | çarpıklık | basıklık |
|---|---|---|---|---|---|
| **gerçek (son)** | **0.1674** | **0.02420** | **0.10385** | +0.468 | +0.134 |
| 152-keskin | 0.1280 | 0.00110 | 0.04240 | +0.664 | −0.145 |
| NK-keskin | 0.1317 | 0.00312 | 0.04516 | +0.665 | −0.003 |
| **164-Skeskin** (c=0) | **0.1696** | 0.00492 | 0.08452 | +0.737 | +0.170 |
| **164-Hkeskin** (c=−½) | 0.1861 | **0.02838** | 0.16247 | +0.318 | −0.098 |
| 152-A4 | 0.1128 | 0.00013 | 0.03550 | +0.698 | +0.685 |
| 164-SA4 | 0.1273 | 0.00013 | 0.03802 | +0.885 | +1.153 |

> 152'nin en eski kusurlarından biri, sentetik gazın **kısa gapları
> hiç üretmemesiydi** (`P(s<0.3)`: gerçek 0.0242, keskin **0.0011** —
> yirmi kat az). Sadakatli inşa bunu 0.0049'a, doğru seviye
> konvansiyonunda (`c = −½`) **0.0284'e** taşıyor — gerçeğin
> **%17 üstü**. `P(s<0.5)`'te ise `c = −½` aşırı düzeltiyor
> (0.162 vs 0.104) ve `c = 0` daha iyi (0.085). **İki konvansiyon
> gap dağılımının farklı yerlerini tutturuyor; hiçbiri hepsini
> tutturmuyor ve bu böyle yazılıyor.**


### 3g. DOĞRU KONVANSİYONDA (`c = −½`) İNŞA EDİLEN GAZLAR

`Hkeskin` ve `HA4` tam olarak §2'nin algoritmasıyla, yalnız seviye
`n − ½` alınarak kuruldu. Yakınsama aynı: **maks|F| = 1.863e−09,
aşan tekne 0, sıralılık TAM** (min Δz = 0.1067 / 0.1499).

| q | τ | **gerçek (son)** | **Hkeskin** (c=−½) | Skeskin (c=0) | **HA4** (c=−½) | SA4 (c=0) |
|---|---|---|---|---|---|---|
| 2 | 0.0576 | 1.004 / +0.0000 | **1.005 / −0.0002** | 1.003 / −0.1282 | **1.006 / −0.0000** | 0.997 / −0.0872 |
| 3 | 0.0913 | 1.008 / −0.0000 | **1.014 / −0.0004** | 1.000 / −0.2087 | **1.015 / +0.0002** | 0.991 / −0.1337 |
| 5 | 0.1338 | 1.017 / +0.0000 | **1.023 / −0.0003** | 1.011 / −0.2826 | **1.031 / +0.0001** | 0.979 / −0.1876 |
| 7 | 0.1618 | 1.024 / −0.0000 | **1.036 / −0.0004** | 1.008 / −0.3457 | **1.044 / +0.0002** | 0.968 / −0.2196 |
| 11 | 0.1993 | 1.035 / +0.0000 | **1.052 / +0.0008** | 1.010 / −0.4138 | **1.066 / −0.0006** | 0.957 / −0.2666 |
| 101 | 0.3836 | 1.116 / −0.0004 | **1.162 / −0.0004** | 1.033 / −0.7038 | **1.232 / −0.0006** | 0.833 / −0.4647 |
| 1009 | 0.5750 | 1.215 / −0.0083 | **1.321 / +0.0001** | 0.979 / −0.9758 | **1.536 / −0.0004** | 0.614 / −0.8675 |

*(hücre = `|c_q|/b_q` / `arg c_q`)*

> **`c = −½`'de FAZ TAM OLARAK SIFIR** (|arg| ≤ 8e−04, gerçek gazın
> kendisiyle aynı mertebe) **ve genlik profili gerçeğin τ-ile-YÜKSELEN
> profilini izliyor** (1.005 → 1.32; gerçek 1.004 → 1.22; `c = 0` düz
> 1.00 → 0.98). Sadakat ölçütü (≥0.98) `Hkeskin`'de τ ≤ 0.90'ın
> **bütün bantlarında** sağlanıyor (1.010 … 1.397). Aşırı-düzeltme
> vardır ve saklanmıyor: yüksek τ'da Hkeskin gerçeğin **%11'i
> üstünde** (1.397 vs 1.259).



---

## 4. ×1.4'ÜN YENİDEN YARGILANMASI — 152'nin ÖZGÜN ÖLÇÜSÜYLE

152'nin "faz oranı" tam olarak `arg Γ_rot(sentetik)/arg Γ_rot(gerçek)`tir,
yani 155/158'in `phi` alanı. `164_configs/164_faz.py` ölçümü kopyalamaz;
`155_cekirdek.olc155`'i **152'nin konvansiyonuyla** çağırır (bantlar
0.525–0.85, η zinciri taban 0.52 / cap 720, aday çizgiler `pk(e^{0.86L})`,
220 örnek, tohum 21, `gap ≥ 2.5·dres`).

**MEŞRUİYET (V8).** Bu koşu 152'nin kayıtlı Γ'sını `keskin`'de
**3.3e−16**, `A4`'te **1.1e−16** ile yeniden üretiyor; `son`'da 152'nin
üç haneli GERCEK satırıyla fark ≤ **4.8e−04** ve S3 üçlüsü
**σ_ds² = 0.1674 / σ_η² = 0.0227 / c₁ = −0.01158** — 152'nin "gerçek"
referansının aynısı.

### 4a. FAZ ORANI — 152'NİN KONVANSİYONUNDA (`c = 0`)

*(§4d aynı tabloyu doğru seviye konvansiyonunda, `c = −½`, tekrarlar.)*

| gaz | τ=0.5375 | τ=0.585 | τ=0.66 | τ=0.74 | τ=0.815 | **ort (b1–4)** |
|---|---|---|---|---|---|---|
| **gerçek (son)** | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | **1.000** |
| 152-keskin | 2.110 | 1.437 | 1.319 | 1.418 | (sarmış) | **1.571** |
| NK-keskin (kontrol) | 2.064 | 1.423 | 1.311 | 1.415 | (sarmış) | **1.553** |
| **164-Skeskin** | **1.209** | **1.062** | **1.058** | **1.187** | (sarmış) | **1.129** |
| 152-A4 | 2.327 | 1.505 | 1.319 | ‡ | (sarmış) | **1.717** |
| **164-SA4** | **1.547** | **1.176** | **1.113** | 1.569 | ‡ | **1.351** |

(‡ = `|Γ| > 1.5`, 152'nin patlak ölçütü — ortalamaya girmiyor. τ = 0.815
bandı 152'de de sekiz koşunun altısında fazını sarıyordu; hükme girmez.)

Ham fazlar (`φ_Γ`):

| gaz | 0.5375 | 0.585 | 0.66 | 0.74 | 0.815 |
|---|---|---|---|---|---|
| **son** | +0.2556 | +0.7251 | +1.3810 | +1.8595 | +1.7274 |
| 152-keskin | +0.5395 | +1.0417 | +1.8216 | +2.6364 | −2.9446 |
| NK-keskin | +0.5276 | +1.0321 | +1.8099 | +2.6309 | −2.9722 |
| **164-Skeskin** | **+0.3092** | **+0.7699** | **+1.4606** | +2.2072 | −2.5225 |
| 152-A4 | +0.5947 | +1.0912 | +1.8214 | −2.4989 | −2.7252 |
| **164-SA4** | **+0.3954** | **+0.8525** | **+1.5376** | +2.9168 | −2.7365 |

> **HÜKÜM (×1.4, birinci aşama).** 152'nin en dayanıklı saydığı sayı —
> *"τ = 0.66'da yedi koşuda 1.30–1.37, %5'ten dar; hiçbir knob
> kıpırdatmıyor"* — **yalnız kök seçimi onarılınca 1.058'e düştü.**
> Bant ortalaması **1.571 → 1.129** (fazlalığın **%77'si**); A4'te
> 1.717 → 1.351 (**%51**). §4d bunun geri kalanını da kapatıyor.

### 4b. KONTROL: açığı taşıyan şey KÖK SEÇİMİ, yakınsama DEĞİL

`NKkeskin` = aynı denklem, 152'nin çözümünden başlatılan **sönümsüz**
Newton (yani "en yakın kök", tam yakınsama hedefiyle):

| gaz | kök kuralı | yakınsama | **genlik (q=2)** | **faz oranı** | σ_ds² |
|---|---|---|---|---|---|
| 152-keskin | en yakın (sönümlü) | medyan 2.3e−10, **maks 1.30** | **0.538** | **1.571** | 0.1280 |
| NK-keskin | en yakın (sönümsüz) | medyan 2.3e−10, maks 4.19 | **0.553** | **1.553** | 0.1317 |
| **164-Skeskin** | **İLK KÖK** | maks **1.9e−09** | **1.003** | **1.129** | **0.1696** |
| gerçek (son) | — | — | 1.004 | 1.000 | 0.1674 |

> Yakınsamayı düzeltmek tek başına **hiçbir şey yapmıyor** (0.538 → 0.553;
> 1.571 → 1.553). Bütün etkiyi **ilk-kök seçimi** taşıyor. (NK'nin
> sıralılığı garanti değil ve bozuldu: sıralamadan önce **872** sıra bozan
> çift, sonra 2 özdeş tekne — ilk-kök kuralının gerekliliğinin ayrıca
> kanıtı.)

### 4c. S3, bağlaşım R ve 156'nın n_eff REJİMİ — hepsi birlikte kapanıyor

| gaz | σ_ds² | σ_η² | c₁ | σΔ² | R(0.5375) | R(0.585) | R(0.66) | n_eff(0.585) | n_eff(0.66) |
|---|---|---|---|---|---|---|---|---|---|
| **gerçek (son)** | **0.1674** | **0.0227** | **−0.01158** | 0.0547 | +0.391 | +1.086 | +1.852 | 153 942 | 41 855 |
| 152-keskin | 0.1280 | 0.0767 | −0.03535 | 0.0350 | +1.574 | +2.721 | +3.626 | 18 590 | **1 824** |
| NK-keskin | 0.1317 | 0.0801 | −0.03494 | 0.0372 | +1.397 | +2.256 | +2.756 | 20 213 | 3 080 |
| **164-Skeskin** | **0.1696** | 0.0386 | **−0.01384** | 0.0572 | **+0.463** | **+1.133** | **+1.842** | **139 641** | **35 825** |
| 152-A4 | 0.1128 | 0.0230 | −0.00918 | 0.0429 | +1.215 | +2.034 | +3.867 | 59 808 | 9 029 |
| **164-SA4** | 0.1273 | 0.0212 | −0.00649 | 0.0518 | +0.677 | +1.308 | +2.218 | 121 553 | 35 773 |

Üç ayrı 152–156 hükmü aynı anda düşüyor:

1. **152'nin "sentetik aynı anda hem az dağınık hem fazla dik" ikilemi
   YOK.** Sadakatli keskin σ_ds² = **0.1696** (gerçek 0.1674, **%+1.3**);
   152-keskin 0.1280 (**%−24**). Aynı gazda faz oranı da 1.13.
2. **156'nın `R_lad(sent)/R_lad(gerçek) = 1.9–3.3` açığı KAPANDI:**
   Skeskin/gerçek = **1.18 / 1.04 / 0.99**; 152-keskin/gerçek =
   4.03 / 2.51 / 1.96.
3. **156'nın "gerçek YUMUŞAK rejimde, sentetik UÇ-DEĞER rejiminde"
   hükmü DÜŞTÜ.** τ = 0.66'da `n_eff`: gerçek 41 855, 152-keskin
   **1 824**, **sadakatli keskin 35 825**. Aynı rejim.
4. Kalan tek belirgin S3 farkı **σ_η²**: Skeskin 0.0386, gerçek 0.0227
   (152-keskin 0.0767). Yani sadakatli keskin gerçeğe **iki kat** daha
   yakın ama hâlâ %70 yüksek; erfc penceresi (SA4: 0.0212) onu gerçeğin
   biraz altına indiriyor — 152'nin A4'ü seçme gerekçesi (S3'ü oturtmak)
   sadakatli inşada da geçerli.


### 4d. DOĞRU KONVANSİYONDA (`c = −½`) FAZ ORANI — AÇIK KAPANIYOR

| gaz | τ=0.5375 | τ=0.585 | τ=0.66 | τ=0.74 | τ=0.815 | **ort (b1–4)** |
|---|---|---|---|---|---|---|
| **gerçek (son)** | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | **1.000** |
| 152-keskin | 2.110 | 1.437 | 1.319 | 1.418 | (sarmış) | **1.571** |
| 164-Skeskin (c=0) | 1.209 | 1.062 | 1.058 | 1.187 | (sarmış) | **1.129** |
| **164-Hkeskin (c=−½)** | **0.871** | **0.945** | **0.971** | **1.003** | 1.158 | **0.948** |
| 152-A4 | 2.327 | 1.505 | 1.319 | ‡ | (sarmış) | **1.717** |
| 164-SA4 (c=0) | 1.547 | 1.176 | 1.113 | 1.569 | ‡ | **1.351** |
| **164-HA4 (c=−½)** | **1.259** | **1.096** | **1.065** | 1.064 | ‡ | **1.121** |

Ham fazlar ve ölçüm sağlığı:

| gaz | φ(0.5375) | φ(0.585) | φ(0.66) | φ(0.74) | \|Γ\| (aynı sırayla) |
|---|---|---|---|---|---|
| **son** | +0.2556 | +0.7251 | +1.3810 | +1.8595 | 0.813 / 0.735 / 0.625 / 0.520 |
| **Hkeskin** | **+0.2226** | **+0.6855** | **+1.3409** | **+1.8656** | **0.811 / 0.732 / 0.630 / 0.549** |
| **HA4** | +0.3219 | +0.7948 | +1.4704 | +1.9788 | 0.806 / 0.737 / 0.616 / **0.215** |
| 152-keskin | +0.5395 | +1.0417 | +1.8216 | +2.6364 | 0.807 / 0.794 / 0.790 / 0.836 |

> **HÜKÜM (kesin biçimi).** 152'nin ×1.4'ü, **iki ayrı inşa kusurunun
> toplamıydı**: (i) kök seçimi (1.571 → 1.129) ve (ii) yarım-seviye
> konvansiyon kayması (1.129 → **0.948**). **Doğru kurulmuş saf asal
> merdiven gazı, gerçek zeta sıfırlarının dispersiyon fazını %5 içinde
> — hatta hafifçe ALTINDA — üretiyor.** `|Γ|` de aynı anda oturuyor
> (0.811/0.732/0.630 vs gerçek 0.813/0.735/0.625, **≤%1**), yani bu bir
> faz-eşleme tesadüfü değil, tam kompleks Γ'nın eşleşmesidir.
>
> (HA4'ün τ = 0.74 bandında `|Γ| = 0.215` ve `±φ = 0.178`; o hücrenin
> hata payı ±%10'dur ve ortalamaya bu belirsizlikle giriyor.)

### 4e. S3 ÜÇLÜSÜ — `HA4` GERÇEĞE OTURUYOR

| gaz | σ_ds² | σ_η² | c₁ | σΔ² |
|---|---|---|---|---|
| **gerçek (son)** | **0.1674** | **0.0227** | **−0.01158** | **0.0547** |
| 152-keskin | 0.1280 | 0.0767 | −0.03535 | 0.0350 |
| 152-A4 | 0.1128 | 0.0230 | −0.00918 | 0.0429 |
| 164-Skeskin (c=0) | 0.1696 | 0.0386 | −0.01384 | 0.0572 |
| 164-SA4 (c=0) | 0.1273 | 0.0212 | −0.00649 | 0.0518 |
| **164-Hkeskin (c=−½)** | 0.1861 | 0.0319 | −0.01779 | 0.0586 |
| **164-HA4 (c=−½)** | 0.1886 | **0.0229** | **−0.01041** | 0.0574 |

> `HA4`: **σ_η² %0.9, c₁ %10, σΔ² %5** farkla gerçeğin üstünde/altında;
> tek belirgin sapma σ_ds² (**%+13**). 152'nin A4'ü S3'ü zaten
> oturtuyordu (σ_η² 0.0230) **ama faz oranı 1.72 idi**; sadakatli inşada
> aynı gaz S3'ü koruyor **ve** faz oranını **1.121**'e getiriyor.
> **152'nin "iki gözlenebilir ters yönde hareket ediyor, ortak bir
> sebeple açıklanamazlar" hükmü DÜŞTÜ: ortak sebep inşaymış.**

### 4f. BAĞLAŞIM R ve n_eff — `Hkeskin` gerçeğin üstüne oturuyor

| gaz | R(0.5375) | R(0.585) | R(0.66) | R(0.74) | n_eff(0.585) | n_eff(0.66) | n_eff(0.74) |
|---|---|---|---|---|---|---|---|
| **gerçek (son)** | +0.391 | +1.086 | +1.852 | +2.120 | 153 942 | 41 855 | 17 361 |
| 152-keskin | +1.574 | +2.721 | +3.626 | (sarmış) | 18 590 | **1 824** | **58** |
| **164-Hkeskin** | **+0.308** | **+0.917** | **+1.661** | **+2.055** | **177 908** | **56 543** | **19 337** |
| **164-HA4** | +0.487 | +1.144 | +1.830 | +2.096 | 137 169 | 41 031 | 19 798 |
| 164-Skeskin | +0.463 | +1.133 | +1.842 | (sarmış) | 139 641 | 35 825 | 2 254 |

> 156'nın **"gerçek gaz fazını YUMUŞAK rejimde üretiyor, sentetik
> uç-değer rejiminde"** hükmü: `n_eff(τ=0.74)` gerçek 17 361,
> 152-keskin **58**, sadakatli-doğru-konvansiyon **19 337 / 19 798**.
> **Rejim farkı yok.** `R` de τ = 0.74'e kadar gerçeğin ±%20'si içinde
> (152-keskin'de faz sarıyordu).


---

## 5. HÜKÜM TABLOSU — eski inşa ↔ sadakatli inşa ↔ gerçek

### 5a. δ PARMAK İZİ (158 ızgarası, taban 0.40, W-A, `w = 1/σ_φ`)

| gaz | **δ(½)** | ±jk | **dδ/dτ\|½** | ±jk | **b** | ±jk | a_τ₀ | **τ₀\*** |
|---|---|---|---|---|---|---|---|---|
| **gerçek (son)** | **−0.0936** | 0.0012 | **−1.815** | 0.015 | **−7.759** | 0.372 | 10.615 | **0.5087** |
| 152-keskin | +0.1329 | 0.0048 | −0.993 | 0.045 | −3.811 | 0.853 | 11.660 | 0.4885 |
| NK-keskin | +0.1209 | 0.0045 | −1.031 | 0.050 | −2.844 | 0.820 | 11.595 | 0.4895 |
| **164-Skeskin** | **−0.0642** | 0.0021 | **−2.113** | 0.043 | **−6.647** | 0.734 | **10.372** | **0.5061** |
| 152-A4 | +0.1347 | 0.0033 | +0.602 | 0.039 | +5.541 | 1.100 | 13.054 | 0.4898 |
| **164-SA4** | **−0.0996** | 0.0020 | **−1.856** | 0.050 | **−15.322** | 0.856 | 10.421 | **0.5093** |
| **164-Hkeskin** (c=−½) | **−0.1189** | 0.0013 | **−1.917** | 0.021 | **−7.897** | 0.449 | **10.472** | **0.5112** |
| **164-HA4** (c=−½) | **−0.0528** | 0.0014 | **−1.777** | 0.018 | **−6.384** | 0.465 | **10.727** | **0.5049** |

> **161'İN İKİ ANA SINIFLANDIRMASI DÜŞTÜ.**
> * 161 §9: *"δ(½) < 0 olan tek gaz gerçektir (τ₀ > ½ olan tek gaz
>   olmasının aynısı)"* — **sadakatli her iki gaz da δ(½) < 0 ve
>   τ₀\* > ½** (0.5061 / 0.5093, gerçek 0.5087).
> * 161 §9'un dörtlü sınıfı: **(−,−,−) yalnız gerçeğin** idi;
>   sadakatli keskin de **(−,−,−)**, sadakatli A4 de **(−,−,−)**.

### 5a-2. KONVANSİYON BULUTU (32 fit: 4 taban × 4 pencere × 2 ağırlık) ve AYIRMA GÜCÜ

| gaz | **δ(½) ± σ_k** | **dδ/dτ\|½ ± σ_k** | **b ± σ_k** | a_τ₀ ± σ_k | **τ₀\* ± σ_k** |
|---|---|---|---|---|---|
| **gerçek (son)** | **−0.0685 ± 0.0510** | **−1.801 ± 0.050** | **−5.980 ± 1.295** | 10.691 ± 0.071 | **0.5064 ± 0.0047** |
| 152-keskin | +0.1344 ± 0.0032 | −1.039 ± 0.094 | −3.984 ± 1.382 | 11.621 ± 0.120 | 0.4883 ± 0.0002 |
| NK-keskin | +0.1238 ± 0.0042 | −1.058 ± 0.098 | −3.472 ± 1.698 | 11.584 ± 0.126 | 0.4892 ± 0.0003 |
| **164-Skeskin** | **−0.0361 ± 0.0441** | **−2.126 ± 0.038** | **−6.760 ± 1.606** | **10.401 ± 0.068** | **0.5035 ± 0.0042** |
| 152-A4 | +0.1374 ± 0.0046 | +0.473 ± 0.231 | +4.050 ± 2.633 | 12.961 ± 0.218 | 0.4895 ± 0.0004 |
| **164-SA4** | **−0.0940 ± 0.0136** | **−1.955 ± 0.716** | **−13.256 ± 5.092** | 10.361 ± 0.799 | **0.5089 ± 0.0015** |
| **164-Hkeskin** (c=−½) | **−0.0958 ± 0.0499** | **−1.909 ± 0.049** | **−6.181 ± 1.247** | — | — |
| **164-HA4** (c=−½) | **−0.0256 ± 0.0477** | **−1.753 ± 0.050** | **−5.543 ± 0.543** | — | — |

Gerçek havuz (son + orta, **64 fit**): ort **(−0.0692, −1.7962, −6.0187)**,
havuz sd **(0.0508, 0.0564, 1.2263)** — 161 §3b'nin kayıtlı satırıyla
**birebir** (V3).

**Gerçekten uzaklık** (gerçek havuzun kendi sd'si birimiyle; 3B
Mahalanobis gerçek havuzun 64-fit kovaryansıyla — *161'in normalizasyonu
farklı olduğundan mutlak sayılar 161'inkiyle karşılaştırılamaz, aynı
metrikle GAZDAN GAZA karşılaştırma geçerlidir*):

| gaz | δ(½) | dδ/dτ | b | **3B Mahalanobis** |
|---|---|---|---|---|
| 152-keskin | 4.0σ | 13.4σ | 1.7σ | **16.1** |
| NK-keskin | 3.8σ | 13.1σ | 2.1σ | **15.9** |
| **164-Skeskin** | **0.7σ** | 5.9σ | **0.6σ** | **6.5** |
| 152-A4 | 4.1σ | 40.3σ | 8.2σ | **47.9** |
| **164-SA4** (c=0) | **0.5σ** | **2.8σ** | 5.9σ | **8.1** |
| **164-Hkeskin** (c=−½) | **0.5σ** | **2.0σ** | **0.1σ** | **2.3** |
| **164-HA4** (c=−½) | **0.9σ** | **0.8σ** | **0.4σ** | **1.4** |

> Aynı metrikte sadakatli inşa uzaklığı **16.1 → 6.5** (keskin, c=0),
> **47.9 → 8.1** (A4, c=0) ve doğru seviye konvansiyonunda
> **16.1 → 2.3** (Hkeskin), **47.9 → 1.4** (HA4) düşürüyor.
> `c = 0` gazlarında geriye kalan tek ayırıcı eksen `dδ/dτ`'dur ve
> orada sadakatli gaz gerçeği **AŞAR** (−2.126 vs −1.801); `c = −½`'de
> o eksen de kapanıyor (Hkeskin 2.0σ, HA4 0.8σ).
>
> > **HÜKÜM: `c = −½` ile kurulmuş SAF ASAL MERDİVEN gazı, 161'in
> > üç-eksenli parmak izinde gerçek zeta sıfırlarının konvansiyon
> > bulutunun İÇİNDEDİR** (HA4 1.4σ, Hkeskin 2.3σ). 161'in
> > "parmak izi üçlüsü gerçeği yedi sentetiğin hepsinden ayırıyor"
> > hükmü, **inşa onarılınca ayakta kalmıyor.**


### 5b. δ BANT BANT (taban 0.40; δ = φ_Γ − (4π·τ_eff − 2π))

| τ̄ | **gerçek** | 152-keskin | **164-Skeskin** (c=0) | **164-Hkeskin** (c=−½) | 152-A4 | **164-SA4** (c=0) | **164-HA4** (c=−½) |
|---|---|---|---|---|---|---|---|
| 0.41 | **+0.0263** | +0.1919 | **+0.0666** | **+0.0070** | +0.0984 | +0.0033 | **+0.0626** |
| 0.43 | **+0.0048** | +0.1810 | **+0.0462** | **−0.0145** | +0.1102 | −0.0152 | **+0.0431** |
| 0.45 | **−0.0200** | +0.1676 | **+0.0224** | **−0.0426** | +0.1189 | −0.0473 | **+0.0200** |
| 0.47 | **−0.0465** | +0.1645 | **−0.0037** | **−0.0693** | +0.1346 | −0.0695 | **−0.0060** |
| 0.49 | **−0.0818** | +0.1555 | **−0.0400** | **−0.1084** | +0.1321 | −0.1029 | **−0.0395** |
| 0.51 | **−0.1152** | +0.1126 | **−0.0883** | **−0.1402** | +0.1453 | −0.1249 | **−0.0698** |
| 0.53 | **−0.1593** | +0.1018 | **−0.1371** | **−0.1875** | +0.1593 | −0.1723 | **−0.1130** |
| 0.55 | **−0.2089** | +0.0715 | **−0.1855** | **−0.2397** | +0.1741 | −0.2138 | **−0.1602** |
| 0.57 | **−0.2568** | +0.0528 | **−0.2532** | **−0.2913** | +0.1960 | −0.2912 | **−0.2097** |
| 0.59 | **−0.3123** | −0.0006 | **−0.3020** | **−0.3467** | +0.2365 | −0.3864 | **−0.2577** |
| 0.61 | **−0.3895** | −0.0195 | **−0.3759** | **−0.4274** | +0.2756 | −0.5286 | **−0.3261** |
| 0.63 | **−0.4545** | −0.0495 | **−0.4671** | **−0.4905** | +0.4527 | −0.9973 | **−0.3875** |

> **Gerçek gaz, iki sadakatli gazın ARASINDA kalıyor.** `Hkeskin`
> (c=−½, penceresiz) her bantta gerçeği **%8–15 AŞIYOR**; `HA4`
> (c=−½, erfc) her bantta **%15–25 ALTINDA**; ikisi gerçeği
> **kuşatıyor** (ör. τ̄=0.57: −0.2913 / **−0.2568** / −0.2097).
> `Skeskin` (c=0) τ̄ ≥ 0.57'de gerçeğin %1.4–3.5'i içinde ama düşük
> τ'da işaret ayrılıyor. Eski inşa ise **her bantta** işaretçe ayrı
> (+0.19 … −0.05).
>
> **SÖNÜM DE EŞLEŞİYOR — ve bu koşunun en çarpıcı tek sayısı budur.**
> `|Γ|` on iki bantta:
>
> | | 0.41 | 0.43 | 0.45 | 0.47 | 0.49 | 0.51 | 0.53 | 0.55 | 0.57 | 0.59 | 0.61 | 0.63 |
> |---|---|---|---|---|---|---|---|---|---|---|---|---|
> | **gerçek** | 0.922 | 0.916 | 0.912 | 0.903 | 0.892 | 0.885 | 0.870 | 0.857 | 0.844 | 0.831 | 0.797 | 0.787 |
> | **Hkeskin** | 0.926 | 0.918 | 0.912 | 0.904 | 0.894 | 0.883 | 0.870 | 0.856 | 0.840 | 0.827 | 0.796 | 0.784 |
> | **HA4** | 0.909 | 0.905 | 0.900 | 0.893 | 0.886 | 0.881 | 0.866 | 0.855 | 0.845 | 0.835 | 0.805 | 0.796 |
> | 152-keskin | 0.807 | 0.801 | 0.818 | 0.801 | 0.785 | 0.786 | 0.786 | 0.782 | 0.775 | 0.784 | 0.806 | 0.795 |
>
> On iki bantta gerçekten **maks \|Δ\|Γ\|\|**:
> **`Hkeskin` 0.0042**, `HA4` 0.0126, `Skeskin` 0.0954,
> `SA4` 0.1319, 152-keskin 0.1153.
> Sönüm profili (`|Γ|` ilk → son bant): gerçek **0.922 → 0.787**,
> `Hkeskin` **0.926 → 0.784**, `HA4` 0.909 → 0.796; buna karşılık
> `Skeskin` 0.883 → **0.882** ve 152-keskin 0.807 → **0.795** —
> **`c = 0` gazları τ ile SÖNMÜYOR bile.**
> Yani eşleşen şey yalnız faz değil, **TAM KOMPLEKS Γ**'dır ve bunu
> sağlayan şey **seviye konvansiyonudur**, σ_ds² değil (`Skeskin`
> σ_ds²'de gerçeğe daha yakın — 0.1696 vs 0.1861 — ama `|Γ|`'de
> yirmi kat uzak).

### 5c. SENARYO HÜKMÜ

| senaryo | öngörü | ölçülen |
|---|---|---|
| **(H-İ)** açık büyük ölçüde İNŞA açığı | sadakatli gaz gerçeğe yaklaşır | faz oranı **1.571 → 1.129** (kök seçimi) **→ 0.948** (seviye konvansiyonu); δ parmak-izi uzaklığı **16.1 → 2.3**; τ₀\* **0.4885 → 0.5112** (gerçek 0.5087); `\|Γ\|` bant bant **≤%1**; R oranı **1.9–4.0 → 0.79–0.97** (Hkeskin) / **0.99–1.25** (HA4); n_eff(0.74) **58 → 19 337** (gerçek 17 361); σ_η² **0.0767 → 0.0229** (HA4; gerçek 0.0227); **maks \|Δ\|Γ\|\| = 0.0042** |
| **(H-F)** açık gerçek fizik, inşa aklanır | fark kalır | **KALMADI.** Bu senaryo düştü |
| **(kısmi)** | payı ölç | kalan pay: `Hkeskin` fazı gerçeğin **%5 ALTINDA** (0.948) ve δ'yı %8–15 aşıyor; `HA4` fazı %12 üstünde ve δ'yı %15–25 altında; **gerçek ikisinin ARASINDA.** σ_ds² her iki c=−½ gazında **%11–13 yüksek** |

> **HÜKÜM: senaryo H-İ, tam biçimiyle.** 152'nin ×1.4'ü **iki ayrı inşa
> kusurunun toplamıydı**:
> 1. **kök seçimi** — `F = n`'in çok kökü var, 152 seviyelerin %22'sinde
>    yanlışını seçmiş (1.571 → 1.129);
> 2. **yarım-seviye konvansiyon kayması** — doğru denklem `N̄+S = n−½`
>    (1.129 → 0.948).
>
> Kalan fark **fizik değil, gazın penceresi**: penceresiz (`Hkeskin`)
> gerçeği aşıyor, erfc-pencereli (`HA4`) altında kalıyor, gerçek
> **ikisinin arasında**. Yani "eksik malzeme" arayışı (GUE itmesi,
> titreşim, kısa-menzil yapı — 153/155'in bütün knob'ları) artık
> **gereksizdir**: saf asal merdiven, doğru çözülünce yetiyor.

---

## 6. 155/161'İN KORELASYON-PROBU OKUMASI — REVİZYON ZORUNLU

155/157/158/161 zinciri şu merdiveni kurmuştu (τ₀ artan):
`keskin < A4 < J14 < J26 < N5 < N5z < P1 < son/orta` ve okuması
*"τ₀, gazın kısa-menzil korelasyon yapısını izliyor; katı/yapısız
gazlar ½'nin altında, gerçek üstünde"* idi.

| gaz | τ₀\* (161, eski inşa) | τ₀\* (164, sadakatli c=0) | τ₀\* (164, sadakatli c=−½) |
|---|---|---|---|
| keskin | 0.4883 | 0.5061 | **0.5112** |
| A4 | 0.4895 | 0.5093 | **0.5049** |
| **gerçek (son)** | **0.5064** | — | — |

> **Merdivenin İKİ EN ALT BASAMAĞI, inşa onarılınca gerçeğin
> ETRAFINA yerleşiyor.** "τ₀ < ½ ⇒ sentetik" kuralı sentetik gazın
> fiziğinden değil, **çözücünün kusurundan** geliyordu.
>
> Düşen hükümler (hepsi 161'in kendi cümleleriyle):
> * §9 *"δ(½) < 0 olan tek gaz gerçektir"* — sadakatli **dört gazın
>   dördü de** δ(½) < 0.
> * §9 *"(−,−,−) sınıfı yalnız son/orta'nın"* — `Skeskin`, `Hkeskin`,
>   `SA4`, `HA4` **hepsi (−,−,−)**.
> * §7 *"parmak izi üçlüsü gerçeği yedi sentetiğin hepsinden
>   ayırıyor"* — `HA4` **1.4σ**, `Hkeskin` **2.3σ**; 161'in en yakın
>   sentetiği (N5z) 2.9σ idi. **Artık ayırmıyor.**
>
> Bunun anlamı iki yönlüdür ve ikisi de kayda geçiyor:
> 1. **KAYIP:** korelasyon probunun sentetiği gerçekten ayırma gücü
>    büyük ölçüde bir ARTEFAKTMIŞ. 161'in 2.9σ'lık N5z ayrımı,
>    "üçlü gerekli" hükmü ve bütün merdiven sıralaması **MERDİVENİ
>    EKSİK KURAN** gazlar üzerinde ölçülmüştü.
> 2. **KAZANÇ — ve bu koşunun asıl bulgusu:** `N̄ + S = n − ½`'in
>    SADAKATLİ çözümü — hiçbir GUE itmesi, hiçbir titreşim, hiçbir
>    ek parametre olmadan, **yalnız asal merdiven** — gerçek zeta
>    sıfırlarının `|Γ|`'sını (%1), fazını (%5), δ tayfını,
>    parmak izini (1.4–2.3σ), τ₀'ını, a'sını, bağlaşım rejimini
>    (n_eff) ve S3 üçlüsünü (σ_η² %1) üretiyor.
>    **Taç sınavın nicel kapanışına açılan kapı budur.**
>
> **J14/J26/N5/N5z/P1/orta/dusuk bu koşuda YENİDEN İNŞA EDİLMEDİ**
> (N5/N5z/J14/J26 doğrudan 153'ün kusurlu `z_taban.npy`'sinden türüyor;
> `P1` kendi gap örnekleyicisini kullandığı için bu kusurdan
> ETKİLENMEZ). 161'in merdiveninin tamamı yeniden ölçülmeden
> 158/161'in sıralama hükümleri geçerli sayılmamalıdır — bu bir
> **açık borçtur**.

---

## 7. DENETİM

| | sınav | sonuç |
|---|---|---|
| **V1** | 152'nin inşası bit düzeyinde yeniden üretiliyor mu? | `eski_keskin` ↔ `scratchpad/152/z_keskin.npy`: **maks\|Δz\| = 0.0e+00**; `eski_A4` ↔ `scratchpad/154/z_A4.npy`: **0.0e+00** |
| **V2** | Sadakat ölçümü 163 §1c'yi yeniden üretiyor mu? | 152-keskin çizgileri **0.538 / 0.540 / 0.556 / 0.547 / 0.544 / 0.571 / 0.491**, 163 §1c ile aynı; A4 **0.846 / 0.852 / 0.843 / 0.832 / 0.824 / 0.732 / 0.433** — aynı |
| **V3** | δ fit zinciri 161'i yeniden üretiyor mu? | 161 §3a birincil satır: son **−0.0936 / −1.815 / −7.759**, keskin **+0.1329 / −0.993 / −3.811**, A4 **+0.1347 / +0.602 / +5.541** — hepsi bu koşuda AYNI. 161 §3b konvansiyon bulutu (32 fit): son **−0.0685 / −1.801 / −5.980**, keskin **+0.1344 / −1.039 / −3.984**, A4 **+0.1374 / +0.473 / +4.050** — hepsi AYNI |
| **V4** | `\|F\| ≤ 1e−8` her teknede | Skeskin **maks 1.863e−09**, aşan tekne **0**; SA4 aynı. Kayan-nokta tabanı `ulp(n) = 2.33e−10` (maks = 8 ulp) |
| **V5** | sıralılık | Skeskin ve SA4: **sıra bozan çift 0**, min Δz = 0.200 ḡ / 0.272 ḡ |
| **V6** | ızgara yakınsaması | `h = 0.015` ↔ `h = 0.00375`: 5000 teknede **0 fark**; `h = 0.04` 4 tekne kaçırıyor |
| **V7** | `ds = −ΔS + ΔF` özdeşliği | beş gazda **maks\|fark\| ≤ 9.9e−10** (`Hkeskin` için `164_tani_artik.py` `c = 0` seviyesini kullanıyor ⇒ `F ≡ +0.5` sabiti çıkıyor; `ΔF` yine 0 ve ayrışım geçerli) |
| **V8** | 152'nin faz ölçüsü bit düzeyinde yeniden üretiliyor mu? | `keskin` Γ farkı **3.3e−16**, `A4` **1.1e−16**; `son` 152'nin üç haneli GERCEK satırıyla **≤ 4.8e−04** ve S3 üçlüsü **0.1674 / 0.0227 / −0.01158** (152'nin referansı) |
| **V9** | seviye periyodikliği (`c` ↔ `c+1` indis ötelemesi) | `c = 0` ↔ `c = −1`: σ_ds² 0.1714 / 0.1715, `arg/ω` −0.1735 / −0.1735 — **özdeş** |
| **V10** | `c = −½` gazlarında da tam yakınsama | `Hkeskin` / `HA4`: maks\|F\| **1.863e−09**, aşan tekne **0**, sıralılık TAM (min Δz 0.1067 / 0.1499) |
| **V11** | `\|Γ\|` eşleşmesi (bağımsız sınav) | gerçek ↔ `Hkeskin`, taban 0.40, 12 bant: **maks \|Δ\|Γ\|\| = 0.0042** |

---

## 8. DÜRÜSTLÜK NOTLARI

* **Uydurma yok; her sayı bu koşuda üretildi.** Tablolar
  `scratchpad/164/`'teki `insa_*.json` (5), `sadakat_*.json` (7),
  `artik_*.json` (4), `faz152_*.json` (6), `F_<gaz>_t<taban>.json`
  (3 gaz × 5 taban) ve `log_*.txt` çıktılarından alındı. 152/161/163'ten
  ALINTILANAN tek şey referans satırlarıdır ve her biri bu koşunun kendi
  ölçümüyle V1–V3, V8'de karşılaştırılmıştır (fark 0.0e+00 … 3.3e−16).
* **Karşılaştırmanın meşruiyeti kurulmuştur.** Eski inşa `164_insa.py`
  içinde satır satır yeniden koşuldu ve `z_keskin.npy` / `z_A4.npy`
  ile **maks|Δz| = 0.0e+00** verdi. Yani "eski ↔ yeni" farkı bir kod
  farkı değil, yalnız KÖK KURALI farkıdır.
* **`|F| ≤ 1e−8` ölçütü SAĞLANDI ama tabanı raporlanıyor.** Ölçülen
  maksimum 1.863e−09'dur ve bu `8 × ulp(n)`'dir; çift duyarlıkta
  `F = (N̄ − n) + S` bundan daha küçük olamaz. Bu bir başarı değil,
  bir SINIRDIR ve saklanmıyor.
* **Görevin "kübik interpolasyon" önerisi UYGULANMADI ve nedeni ölçüldü**
  (§2b.4): `rms S⁗ = 1665` ile kübik spline hatası ölçütün 3000 katı
  olurdu. Izgara yalnız braket için kullanıldı; kök TAM `S` ile bulundu.
  Görevin "adım ≤ dalga boyunun 1/10" ölçütü de **yetmedi** ve bu
  ölçülerek gösterildi (§2c): `h = 0.04`'te 5000 teknenin 4'ü kaçıyor.
* **§2–§5'in `Skeskin`/`SA4` gazları YARIM SEVİYE KAYMIŞ konvansiyonda
  (`c = 0`) kurulmuştur** ve bu **kasıtlıdır**: ×1.4'ün yeniden
  yargılanması 152'nin kendi konvansiyonunda yapılmalıydı. Doğru
  konvansiyon (`c = −½`) §3e'de ölçülerek bulundu ve `Hkeskin`/`HA4`
  ile ayrıca koşuldu (§3g, §4d–4f, §5). **İki konvansiyon da tabloda
  ayrı ayrı duruyor; hiçbiri diğerinin yerine geçirilmedi.**
* **`c = −½` bulgusu bu koşunun İÇİNDE, planlanmamış biçimde çıktı.**
  Sıra şuydu: (i) sadakat denetimi genliği 1.00'a getirdi ama beklenmedik
  bir faz gösterdi (§3d); (ii) ilk hipotez "ilk-kök yanlılığı" idi;
  (iii) seviye taraması onu **çürüttü** ve `c = −½`'yi buldu (§3e);
  (iv) kalem karşılığı (Titchmarsh konvansiyonu) sonradan yazıldı.
  Bu sıra **saklanmıyor**: hipotez önce yanlıştı.
* **Gerçek gaz iki sadakatli gazın ARASINDA kalıyor** (`Hkeskin` fazı
  %5 altında ve δ'yı %8–15 aşıyor; `HA4` fazı %12 üstünde ve δ'yı
  %15–25 altında). **"Tam kapanış" iddiası YOKTUR**; ölçülen şey,
  gerçeğin sadakatli gazların kuşattığı bölgede olduğudur.
  σ_ds² her iki `c = −½` gazında da **%11–13 yüksektir** ve bu tek
  yönlü bir sapmadır.
* **`NKkeskin` kontrolünün sıralılığı BOZUK.** Sönümsüz Newton
  sıralamadan önce **872** sıra bozan çift bıraktı ve `maks|F| = 4.19`
  ile yakınsamadı (birkaç bin tekne). Bu kontrolün amacı "en yakın kök
  kuralında ne olur" sorusuydu ve o soruya cevap verdi (genlik 0.553,
  faz oranı 1.553); **sadakatli bir gaz olarak KULLANILAMAZ** ve
  hiçbir kapanış hükmüne girmedi.
* **A4/SA4'ün yüksek-τ sadakat oranları TANIMSIZ.** erfc penceresi
  τ > 0.65'te nominal genliği sıfıra sürdüğü için `|c|/b` bir
  gürültü/sıfır bölümüdür (2.27, 6.71, 20.3, 209 gibi). Tabloda ⚠ ile
  işaretli ve hiçbir hükme girmedi.
* **`τ = 0.815` bandı hiçbir hükme girmedi** — 152'de de sekiz koşunun
  altısında fazını sarıyordu; bu koşuda beş gazın dördünde sarıyor.
  Faz oranı ortalaması 152'nin konvansiyonuyla **yalnız b1–b4** üzerinden.
* **`SA1` (erfc 0.75/0.10) ve N5-tarzı yumuşatılmış gaz KOŞULMADI.**
  Süre bütçesi sadakatli keskin + A4 (iki konvansiyonda, dört gaz) +
  eski-inşa yeniden üretimi + kök-seçimi kontrolü + seviye taramasına
  gitti (her sadakatli inşa 9.5 dk, ızgara 10.4M nokta × 15450 çizgi).
  `164_insa.py SA1` hazırdır.
* **`HA4`'ün ilk sadakat koşusu penceresiz nominale karşı yapıldı**
  (`PENCERE` sözlüğünde `HA4` yoktu); hata fark edildi, sözlük
  düzeltildi ve koşu tekrarlandı. Rapordaki `HA4` sadakat sayıları
  **pencereli** (doğru) sürümdendir; tek etkilenen satır q = 1009
  (1.356 → 1.536).
* **`HA4`'ün τ = 0.74 bandında `|Γ| = 0.215` ve `±φ = 0.178`** — o
  hücrenin faz oranı (1.064) ±%10 belirsizdir ve ortalamaya bu
  belirsizlikle giriyor. Diğer üç bandı sağlıklıdır.
* **`orta`, `dusuk`, `J14`, `J26`, `N5`, `N5z`, `P1` YENİDEN İNŞA
  EDİLMEDİ.** Hepsi 152/153'ün kusurlu tabanından türüyor
  (`N5/N5z/J14/J26` doğrudan `z_taban.npy`'yi, `P1` kendi gap
  örnekleyicisini kullanıyor — sonuncusu bu kusurdan ETKİLENMEZ).
  161'in dokuz-gaz merdiveni bu yüzden **yeniden ölçülmelidir**;
  bu koşu onu yapmadı ve açık borç olarak yazıyor.
* **Figür üretilmedi;** bu koşuda görselleştirme istenmedi.
* **Süreler.** Sadakatli inşa **9.5 dk/gaz** (ızgara 483 s + Newton
  70–90 s, 7 süreç, 10.4M × 15450); 152-Newton yeniden üretimi
  9.2–9.7 dk/gaz; `NKkeskin` 13.4 dk (30 iterasyon, yakınsamadı);
  sadakat denetimi 39 s/gaz; artık tanısı 105 s/gaz; 152-konvansiyon
  faz koşusu 1.1 dk/gaz; 158-ızgarası δ koşusu ~3–4 dk/(gaz,taban);
  seviye taraması 5 dk (40k tekne × 7 `c`). Toplam ~2 saat.

---

## Sıradaki adım (bu ölçümün işaret ettiği)

1. **BÜTÜN SENTETİK GAZ KÜTÜPHANESİ `c = −½` + İLK-KÖK ile YENİDEN
   İNŞA EDİLMELİ.** Bu koşu iki basamağı (keskin, A4) onardı; kalan
   yedi (J14, J26, N5, N5z, P1, orta-pencere, dusuk-pencere) hâlâ
   152/153'ün kusurlu tabanından geliyor. 155/157/158/161'in bütün
   sıralama hükümleri bu yeniden inşadan SONRA okunmalıdır.
2. **163'ün üçlü-toplamı sadakatli gazda yeniden koşulmalı.** 163'ün
   `P0` sürümü keskin'de sistematik ×6 aşıyordu ve sebebi §1c'nin
   0.54'üydü (güç oranı 0.54² = 0.29 farkı tam açıklıyordu);
   sadakatli gazda `P0 ≡ P0m` olur ve ±%25 sınavı ilk kez **ADİL**
   koşulmuş olur. Ayrıca 163'ün "`arg x_q = πτ_q` belirlenimli merdiven"
   varsayımı `c = −½`'de doğrudan doğrulanmıştır (|arg| ≤ 8e−04).
3. **σ_ds²'nin kalan %13'ü.** `c = −½` gazlarının ikisi de gerçekten
   daha dağınık. Adaylar: (i) merdiven kesimi (τ ≤ 1.00 → daha dar
   bir pencere; `HA4`'ün erfc'si zaten bu yönde ama yetmiyor),
   (ii) merdivenin ikinci mertebe terimleri (`m ≥ 2` kule katkıları),
   (iii) `N̄`'nin daha yüksek RvM terimleri. Sadakatli erfc taraması
   (`SA1/SA2/SA3`'ün `c = −½` sürümleri) ucuzdur ve doğrudan sınar.
4. **`Hkeskin` ↔ `HA4` KUŞATMASI bir FİT KAPISIDIR.** Gerçek gaz
   ikisinin arasında; tek parametreli erfc ailesinde `τ_c`'yi
   δ tayfına oturtan değer aranmalıdır (152'nin A4'ü S3'e oturtmuştu;
   şimdi aynı şey **fazla** yapılabilir — ve S3 zaten oturuyor).
5. **Ψ ⋆ κ formülasyonu** (163 §Sıradaki adım 2) hâlâ açık; bu koşu
   ona dokunmadı.
6. **`P1` bu kusurdan ETKİLENMEZ** (kendi gap örnekleyicisi + tek
   geçiş merdiven boyası kullanıyor, Newton yok). 161'in
   "P1 gerçeğe en yakın sentetik" okuması bu yüzden ayakta olabilir
   ve yeniden inşa edilen merdivenle karşılaştırılmalıdır.

---

## Ek — dosyalar ve tekrar-üretim

| dosya | ne |
|---|---|
| `164_configs/164_insa.py` | teşhis (docstring'de tam), `merdiven`, `coz_sadakatli` (sıralı ilk-kök), `coz_enyakin` (kök-seçimi kontrolü), `coz_eski` (152'nin Newton'u, satır satır) |
| `164_configs/164_sadakat.py` | genlik sadakat denetimi (163'ün yöntemi + bant agregasyonu + gürültü çıkarımı) |
| `164_configs/164_tani_artik.py` | `ds = −ΔS + ΔF` ayrışımı; artık kanalının ölçümü |
| `164_configs/164_faz.py` | 152'nin ÖZGÜN faz-oranı ölçüsü (taban 0.52, cap 720) |
| `164_configs/164_kos.py` | 158 ızgarasında δ ölçümü (155_cekirdek.olc155, kopya değil) |
| `164_configs/164_analiz.py` | hüküm tabloları (161_cekirdek'in fit çekirdeğiyle) |
| `164_configs/164_tani_seviye.py` | seviye konvansiyonu `c` taraması — §3d'nin sınavı, §3e'nin bulgusu (ızgara `c`'den bağımsız, tek pasoyla yedi `c`) |

Ham çıktılar `scratchpad/164/`: `z_<gaz>.npy` (7), `insa_<gaz>.json` (7),
`sadakat_<gaz>.json` (9), `artik_<gaz>.json` (4), `faz152_<gaz>.json` (8),
`F_<gaz>_t<taban>.json` (5 gaz × 5 taban = 25), `seviye_taramasi.json`,
`analiz_cikti.txt`, `log_*.txt`.

**Gazların künyesi**

| ad | denklem | kök kuralı | pencere |
|---|---|---|---|
| `eski_keskin` / `eski_A4` | `N̄+S = n` | 152'nin sönümlü/kelepçeli Newton'u | — / erfc 0.68/0.125 |
| `NKkeskin` | `N̄+S = n` | **en yakın kök** (kontrol) | — |
| `Skeskin` / `SA4` | `N̄+S = n` | **sıralı ilk kök** | — / erfc 0.68/0.125 |
| **`Hkeskin` / `HA4`** | **`N̄+S = n − ½`** | **sıralı ilk kök** | — / erfc 0.68/0.125 |

```
for c in Skeskin SA4 Hkeskin HA4 eski_keskin eski_A4 NKkeskin; do
  .venv/bin/python 164_configs/164_insa.py $c
done
.venv/bin/python 164_configs/164_sadakat.py son keskin A4 Skeskin SA4 Hkeskin HA4 NKkeskin eski_keskin
.venv/bin/python 164_configs/164_tani_artik.py keskin Skeskin A4 SA4
for v in son keskin A4 Skeskin SA4 Hkeskin HA4 NKkeskin; do .venv/bin/python 164_configs/164_faz.py $v; done
for t in 0.40 0.28 0.34 0.46 0.52; do
  for v in Skeskin SA4 NKkeskin Hkeskin HA4; do .venv/bin/python 164_configs/164_kos.py $v $t; done
done
.venv/bin/python 164_configs/164_tani_seviye.py 40000
.venv/bin/python 164_configs/164_analiz.py
```
