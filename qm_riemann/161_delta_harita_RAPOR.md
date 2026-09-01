# 161 — korelasyon probunun YENİ-KİMLİKLİ haritası: δ(τ) tayfları ve parmak-izi üçlüsü

**Hedef.** 159, ölçülen fazın özdeş olarak bir kinematik omurga + mekanizma
fazı olduğunu gösterdi:

> **φ_Γ = (4π·τ_eff − 2π) + δ**,  a = 4π + dδ/dτ,  τ₀ = ½ − δ(½)/a

Fizik artık δ'da. Bu rapor korelasyon probunun **tamamını** yeni kimlikle
yeniden haritalandırıyor: dokuz gaz × beş taban × dört fit penceresi × iki
ağırlık konvansiyonunda δ(τ_eff) tayfları, **δ(½) · dδ/dτ|½ · b** parmak-izi
tablosu, ve probun asıl sorusu — *bu üçlü gerçeği sentetiklerden ayırıyor mu,
158'in a-uzayında ayıramadığı N5z dahil?*

**Figür:** `161_delta_harita.png` (6 panel).
P1 δ tayfları (9 gaz) · P2 taban konvansiyonunun δ'yı ötelemesi ·
**P3 ANA PANEL: (δ(½), dδ/dτ) düzlemi, konvansiyon bulutları + 2σ elipsleri** ·
P4 üçüncü eksen b · P5 tamamlayıcılık · P6 merdiven kesirleri.

**Veri kaynağı ve dürüstlük çerçevesi.** Bu raporda **hiçbir yeni ölçüm
yapılmadı.** Bütün sayılar 158'in `scratchpad/158/F_<veri>_t<taban>.json`
dosyalarındaki **45 kayıtlı koşudan** (bant düzeyinde `tau_eff`, `phi`,
`sPhi_jk` ve 8 jackknife replikası `phi_jk`/`teff_jk`) **yeniden
hesaplandı** — tablodan okunmadı. 155'in `tau0_*_genis.json` koşuları yalnız
**ızgara ekseni** için (T4e) kullanıldı. Tablodan-okunan tek şey 158'in W-A
(τ₀, a, b) satırları ve 159'un δ tablosudur; ikisi de yalnız **denetim
referansı** olarak duruyor ve her biri V1–V4'te bu koşunun kendi
hesabıyla karşılaştırılıyor.

---

> ## Kısa hüküm
>
> 1. **Rotor düzeltmesi bir yeniden-ölçüm değil, TAM bir koordinat
>    değişimidir — ve bu ölçüldü.** Aynı ağırlıkla yapılan fitler özdeş
>    olarak `c_φ = c_δ + (0, 4π, −2π)` ile bağlı: 178 fitte
>    maks |Δkatsayı| = **1.6e−10**. Sonuçları:
>    * **a_τ₀ = 4π + dδ/dτ|τ₀ ≡ a(158)**, maks fark **1.5e−12**;
>    * **τ₀* = ½ − δ(½)/a_½** 158'in τ₀'ını 40 satırda maks **0.00015** ile
>      veriyor — görevin 0.001 eşiğinin **0.15 katı**;
>    * **b(δ) ≡ b(φ) ÖZDEŞ.** b, rotor düzeltmesinin **dokunmadığı tek
>      sayıdır**: yeniden ölçülmedi, aynı sayıdır. (159'un gördüğü
>      |b(φ) − b(δ)| = 0.17…2.42 bir kimlik etkisi değil, **ağırlık
>      konvansiyonu** etkisidir — 159 δ'yı 1/σ_δ ile tartmıştı; §1c.)
> 2. **δ(½), τ₀'ın YENİDEN ETİKETLENMESİDİR — yeni eksen değil.**
>    Dokuz gazda korel(δ(½), τ₀) = **−0.997**, ve δ(½) ≈ a·(½ − τ₀)
>    özdeşliği ≤ 0.009 ile tutuyor. Aynı şekilde dδ/dτ, a'nın kendisidir
>    (korel +0.998; τ₀ çapasında ÖZDEŞ). **Yeni kimlik, yeni bilgi
>    üretmiyor; taşıdığı şey yorum ve hata bütçesidir.**
> 3. **Eksenler daha bağımsız DEĞİL.** korel(δ(½), dδ/dτ) = **+0.785**;
>    158'in (τ₀, a) çiftinde −0.762 idi. Yani "a ve τ₀ örtüşen ama özdeş
>    olmayan iki prob" hükmü δ-uzayında **aynen** duruyor, ne iyileşiyor
>    ne kötüleşiyor.
> 4. **YENİ VE ÖLÇÜLEN KAZANÇ — istatistik hata.** δ'nın kendi
>    jackknife'ı φ'ninkinden küçüktür: medyan σ_δ/σ_φ = **0.38 (gerçek-son)**,
>    **0.20 (orta)**, sentetiklerde 0.60 – 0.95. φ ile τ_eff aynı çizgilerden
>    geldiği için birlikte dalgalanıyorlar; omurga çıkınca **ortak gürültü
>    de çıkıyor**. Gerçek gazda **2.6 – 5 kat** daha hassas bir gözlenebilir.
> 5. **AMA δ(½) ÜÇLÜNÜN EN KIRILGAN AYAĞIDIR.** Gerçek gazda
>    σ_konv(δ(½)) = **0.057** ve bunun **0.0569'u tek başına taban
>    eksenidir** (bağıl %83). Karşılaştırma: dδ/dτ|½ **±0.105 (%5.8)**,
>    a_τ₀ **±0.098 (%0.9)**. Yani **½ çapası kavramsal olarak temiz ama
>    ampirik olarak kırılgan olanıdır**; τ₀ çapası hâlâ en dayanıklı
>    okumadır. 158 §6f'nin hükmü (a'nın farkı taban boyunca %6, τ₀'ınki
>    %83) δ-uzayında **aynen tekrarlanıyor**: Δδ(½)'nin menzil/ortalaması
>    **%74**, Δa_τ₀'ınki **%6.1**.
> 6. **N5z KARŞI-ÖRNEĞİ δ-UZAYINDA AYRIŞIYOR — ama zayıf.**
>    * dδ/dτ|½ ekseninde **0.9 σ_konv** (158'in a'daki 1.1σ'sının aynısı —
>      beklendiği gibi, aynı sayı);
>    * **δ(½) ekseninde 2.7 σ_konv**, ve **tek tek fitlerde ÖRTÜŞME YOK**:
>      N5z'nin 32 fiti [+0.0551, +0.0789], gerçeğin 64 fiti
>      [−0.1337, +0.0051]. (158, a'da örtüşme bulmuştu: N5z min 10.509 <
>      gerçek maks 10.824.)
>    * **2B (δ(½), dδ/dτ): 2.9 σ_konv**; **3B (+b): 2.9 σ_konv.**
>    * Bir-taban-dışarıda sınavında 2.5 – 5.5 arasında kalıyor.
>    **Hüküm: N5z artık ayırt ediliyor, ama 3σ'ın altında.**
> 7. **PARMAK İZİ ÜÇLÜSÜ, GERÇEĞİ YEDİ SENTETİĞİN HEPSİNDEN AYIRIYOR.**
>    3B Mahalanobis (konvansiyon bulutu birimi): **N5z 2.9 · P1 7.3 ·
>    J26 8.3 · J14 9.6 · keskin 9.9 · A4 10.8 · N5 11.1**. Hiçbir TEK
>    eksen bunu yapamıyor (δ(½) P1'de 1.4; dδ/dτ N5z'de 0.9; b N5z'de 0.6),
>    hiçbir İKİLİ de yapamıyor (en iyisi δ(½)+b, taban 2.8).
>    **Üçlü GEREKLİDİR**; yeterliliği **2.9σ ile sınırdadır** (3σ'ın
>    altında — mühür değil). Bu, 158'in "a gerçeği yedi sentetiğin
>    ALTISINDAN ayırır" hükmünü **kapatan ama mühürlemeyen** sayıdır.
> 8. **Üç eksen AYRI gazları yakalıyor (tamamlayıcılık ölçüldü).**
>    dδ/dτ beş gazı 5.2 – 10.3 ile yakalıyor ama N5z'yi (0.9) ve P1'i (2.2)
>    kaçırıyor; **δ(½) N5z'yi** (2.7) yakalıyor; **b P1'i** (3.0) yakalıyor
>    — P1'de 3B kazancı **×2.47**. Üçlünün gücü ortalamadan değil,
>    **ortogonal başarısızlıklardan** geliyor.
> 9. **b'nin işareti dördüncü bir sınıflandırma veriyor ve GERÇEK TEK
>    BAŞINA.** (sign δ(½), sign dδ/dτ, sign b) dörtlü sınıf üretiyor:
>    **(+,+,+)** A4/J14/J26 (erfc-yumuşatılmış aile) · **(+,−,+)** N5, P1 ·
>    **(+,−,−)** keskin, N5z · **(−,−,−)** son, orta. **δ(½) < 0 olan tek
>    gaz gerçektir** (τ₀ > ½ olan tek gaz olmasının aynısı).
> 10. **Taban konvansiyonu δ'yı gerçek gazda neredeyse SALT ÖTELİYOR**
>    (157'nin φ için bulduğunun δ'daki karşılığı): 0.28→0.40 adımında
>    ortalama Δδ = **−0.1021** (sd 0.0095), kalan τ-eğimi −0.132 ± 0.030 —
>    yani eğimin fit penceresi (0.18) boyunca katkısı ötelemenin **%23'ü**.
>    0.34→0.40'ta %15, 0.40→0.46'da %73. **Salt öteleme değil ama öteleme
>    baskın.** Sentetiklerde ötelemeler 10–30 kat küçük (0.003 – 0.06).
> 11. **L-DEĞİŞMEZLİK TAM.** son (L = 12.0296) ↔ orta (L = 11.4638):
>    |Δδ(½)| = 0.0013, |Δ(dδ/dτ)| = 0.0087, |Δb| = 0.077, |Δa_τ₀| = 0.0028,
>    |Δτ₀*| = 0.0001 — **hepsi ≤ 0.08 σ_konv**. Üçlünün üç ayağı da L'den
>    bağımsız.

---

## 1. Yöntem

### 1a. Girdi ve zincir

`161_configs/161_cekirdek.py` hiçbir ölçüm parçası yazmaz; 158'in kayıtlı
JSON'larını okur ve **158_analiz.py'nin fit çekirdeğini birebir yeniden
kurar** (aynı `np.polyfit(x, y, derece, w=1/σ)`, aynı kök seçimi, aynı
8-gruplu jackknife mimarisi, aynı bant pencereleri). Tek eklemesi
**okuma çapasıdır**.

| kaynak | ne | kullanım |
|---|---|---|
| `scratchpad/158/F_<veri>_t<taban>.json` | 45 koşu (9 gaz × 5 taban), bant düzeyi | **BİRİNCİL** |
| `scratchpad/155/tau0_<veri>_t0.4_genis.json` | 155'in ızgarası | yalnız T4e (ızgara ekseni) |
| 158 §4 W-A tablosu · 159 §2c/§2d | tablodan okunan referanslar | yalnız V1/V3/V4 denetimi |

`P0` (Poisson) 158'de "ÖLÇÜLEMEDİ" olduğu için buraya **hiç alınmadı**.

### 1b. Çapa: τ₀ → ½

158 üçlüyü **τ₀'da** okuyordu; τ₀ gaza bağlı bir noktadır, yani her gaz
kendi ekseninde farklı bir yerden okunuyordu. 161 çapayı **τ = ½**'ye alıyor
— omurganın kendi sıfırı, gazdan bağımsız, konvansiyonsuz:

```
δ(½)        = δ fitinin τ=½'deki değeri
dδ/dτ|½     = türevi
b           = ½·δ''(çapa)          (kuadratikte çapadan bağımsız: c₀)
```

Türetilmişler (158'e köprü): `a_½ = 4π + dδ/dτ|½`, `a_τ₀ = 4π + dδ/dτ|τ₀`,
`τ₀* = ½ − δ(½)/a_½`. Çapa farkı **ölçüldü ve özdeş**:
dδ/dτ|½ − dδ/dτ|τ₀ = 2b(½ − τ₀), dokuz gazda maks sapma **9.3e−16** (V6).

### 1c. Ağırlık: iki konvansiyon, ikisi de koşuldu

δ'nın kendi jackknife hatası σ_δ, φ'ninkinden farklıdır (§1d). Yani **iki
meşru ağırlık** var ve ikisi de bir konvansiyon eksenidir:

| | ağırlık | ne yapar |
|---|---|---|
| **w = 1/σ_φ** | 158/159'un ağırlığı | omurga model uzayında kalır ⇒ `c_φ = c_δ + (0,4π,−2π)` **ÖZDEŞ** ⇒ a_τ₀ ve b **tam** olarak 158'inkiler |
| **w = 1/σ_δ** | δ'yı ilkel gözlenebilir sayarsan doğru olan | fit farklılaşır; b'yi **−7.35'e kadar** kaydırıyor (N5z), δ(½)'yi ≤ 0.015, a_τ₀'ı ≤ 0.09 |

**159'un kullandığı ağırlık 1/σ_δ'dır** (bu koşuda tanımlandı: 1/σ_δ ile
δ(½) = −0.0958 / dδ/dτ = −1.998 çıkıyor, 159'un yazdığı −0.0955 / −1.990'a
karşı; 1/σ_φ ile −0.0936 / −1.951). **159'un raporladığı
|b(φ) − b(δ)| = 0.17 … 2.42 farkı bu ağırlık seçiminin sonucudur, kimliğin
değil** — aynı ağırlıkta fark özdeş sıfırdır (V2).

Her iki ağırlık da konvansiyon bütçesine dahil edildi (σ_ağırlık).

### 1d. Pencere ekseni ve ‡ kuralı

158'in dört penceresi aynen: **W-A** τ̄ ∈ [0.43, 0.61] kuadratik (BİRİNCİL) ·
**W-D** τ̄ ≤ 0.59 · **W-B** τ_eff − τ₀ ∈ [−0.075, +0.105] · **W-C** W-A üstünde
kübik.

**‡ kuralı (158'in kuralının δ hâli):** taban 0.52'de W-A'nın en alt bandı
τ̄ = 0.53 → **çapa τ = ½ pencerenin ALTINDA kalıyor**, δ(½) bir
ekstrapolasyondur. Bütün taban-0.52 satırları tablolarda duruyor ama
**hiçbir ortalamaya girmedi**. (Taban 0.52 + W-B beş koşuda bant yetersizliğinden
hiç fit vermiyor: A4/P1/keskin — 360 hedeften **355 fit** kuruldu, 288'inde
çapa pencere içinde.)

**Konvansiyon bütçesi:**
`σ_konv = σ_taban ⊕ σ_pencere ⊕ σ_ağırlık`, birincil değer 4 taban × 4
pencere × 2 ağırlık = **32 fitin ortalamasıdır** (gerçek için son+orta = 64).

---

## 2. T0 — δ(τ_eff) tayfları

### 2a. Dokuz gaz, taban 0.40, 158'in ızgarası

| gaz | 0.41 | 0.43 | 0.45 | 0.47 | 0.49 | 0.51 | 0.53 | 0.55 | 0.57 | 0.59 | 0.61 | 0.63 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| keskin | +0.192 | +0.181 | +0.168 | +0.164 | +0.155 | +0.113 | +0.102 | +0.071 | +0.053 | −0.001 | −0.019 | −0.050 |
| A4 | +0.098 | +0.110 | +0.119 | +0.135 | +0.132 | +0.145 | +0.159 | +0.174 | +0.196 | +0.237 | +0.276 | +0.453 |
| J14 | +0.100 | +0.115 | +0.132 | +0.131 | +0.141 | +0.151 | +0.143 | +0.170 | +0.206 | +0.248 | +0.267 | +0.463 |
| J26 | +0.134 | +0.151 | +0.165 | +0.157 | +0.180 | +0.185 | +0.141 | +0.196 | +0.208 | +0.288 | +0.310 | +0.585 |
| N5 | +0.063 | +0.071 | +0.073 | +0.081 | +0.070 | +0.075 | +0.072 | +0.075 | +0.084 | +0.101 | +0.142 | +0.268 |
| N5z | +0.130 | +0.139 | +0.123 | +0.119 | +0.079 | +0.070 | +0.025 | −0.002 | −0.026 | −0.007 | +0.130 | +0.723 |
| P1 | +0.109 | +0.071 | +0.051 | +0.041 | +0.017 | −0.006 | −0.048 | −0.006 | −0.017 | +0.086 | +0.100 | +0.119 |
| **son** | **+0.026** | **+0.005** | **−0.020** | **−0.046** | **−0.082** | **−0.115** | **−0.159** | **−0.209** | **−0.257** | **−0.312** | **−0.390** | **−0.454** |
| **orta** | +0.028 | +0.006 | −0.019 | −0.048 | −0.080 | −0.118 | −0.157 | −0.209 | −0.253 | −0.313 | −0.381 | −0.457 |

**Fark fitten önce, gözle görünür (figür P1).** Dört ayrı davranış:

* **gerçek (son/orta):** τ̄ = 0.43 ile 0.45 arasında sıfırı geçiyor ve
  **tekdüze** −0.454'e iniyor — δ'sı en derin gaz, farkla;
* **erfc ailesi (A4/J14/J26):** hiç sıfırı geçmiyor, tekdüze **yukarı**
  çıkıyor (+0.10 → +0.31);
* **keskin:** tekdüze iniyor ama sıfırı ancak τ̄ ≈ 0.58'de geçiyor ve
  yalnız −0.050'ye ulaşıyor — gerçeğin yolunun **%11'i**;
* **N5z ve P1:** salınıyor (N5z τ̄ = 0.54'te sıfırı geçip 0.59'da geri
  dönüyor; P1 0.52'de geçip 0.58'de dönüyor).

**Bu, 158'in "φ eğrileri sentetiklerde daha dik" ham gözleminin
kinematikten arındırılmış hâlidir** — ve arındırma, gözlemi tersine
çeviriyor: sentetikler daha dik değil, **gerçek daha derin**.

σ_δ (jackknife) aynı bantlarda: gerçek **0.0012 – 0.0059**, keskin
0.004 – 0.011, A4 0.002 – 0.048, P1 **0.013 – 0.141**. Sentetiklerin
τ̄ ≥ 0.61 bantları (158 §7'nin bant sağlığı kuralı) bozuk ve hiçbir fite
girmiyor (W-A üst sınırı 0.61).

### 2b. Aynı gaz, beş taban — konvansiyonun δ üzerindeki izi

Gerçek gaz (son), δ:

| taban | 0.47 | 0.49 | 0.51 | 0.53 | 0.55 | 0.57 | 0.59 | 0.61 |
|---|---|---|---|---|---|---|---|---|
| 0.28 | +0.049 | +0.022 | −0.016 | −0.059 | −0.105 | −0.151 | −0.208 | −0.264 |
| 0.34 | +0.002 | −0.032 | −0.066 | −0.109 | −0.158 | −0.208 | −0.261 | −0.331 |
| 0.40 | −0.046 | −0.082 | −0.115 | −0.159 | −0.209 | −0.257 | −0.312 | −0.390 |
| 0.46 | −0.081 | −0.115 | −0.149 | −0.190 | −0.239 | −0.283 | −0.334 | −0.408 |
| 0.52 ‡ | — | — | — | −0.189 | −0.230 | −0.273 | −0.322 | −0.379 |

Aynı tabanlarda A4: +0.133 / +0.128 / +0.132 / +0.132 / +0.135 (τ̄ = 0.49
satırı) — **taban ekseninde neredeyse hiç kımıldamıyor**. Nicel sınav §6a'da.

---

## 3. T1 — KESİN PARMAK-İZİ TABLOSU

### 3a. Birincil satır (W-A, taban 0.40, w = 1/σ_φ; jackknife hatalarıyla)

| gaz | **δ(½)** | ±jk | **dδ/dτ\|½** | ±jk | **b** | ±jk | a_½ | a_τ₀ | τ₀* |
|---|---|---|---|---|---|---|---|---|---|
| keskin | +0.1329 | 0.0048 | **−0.993** | 0.045 | **−3.811** | 0.853 | 11.573 | 11.660 | 0.4885 |
| A4 | +0.1347 | 0.0033 | **+0.602** | 0.039 | **+5.541** | 1.100 | 13.168 | 13.054 | 0.4898 |
| J14 | +0.1333 | 0.0023 | +0.548 | 0.081 | +6.212 | 1.508 | 13.114 | 12.987 | 0.4898 |
| J26 | +0.1601 | 0.0059 | +0.449 | 0.103 | +7.330 | 1.788 | 13.016 | 12.834 | 0.4877 |
| N5 | +0.0670 | 0.0027 | +0.076 | 0.041 | +4.092 | 1.080 | 12.642 | 12.599 | 0.4947 |
| N5z | +0.0580 | 0.0026 | −1.178 | 0.070 | +5.761 | 1.334 | 11.388 | 11.329 | 0.4949 |
| P1 | +0.0050 | 0.0127 | −0.573 | 0.167 | +6.835 | 4.948 | 11.993 | 11.988 | 0.4996 |
| **son** | **−0.0936** | 0.0012 | **−1.815** | 0.015 | **−7.759** | 0.372 | 10.751 | 10.615 | 0.5087 |
| **orta** | **−0.0936** | 0.0018 | **−1.774** | 0.043 | **−7.801** | 0.746 | 10.792 | 10.656 | 0.5087 |

### 3b. PARMAK İZİ — konvansiyon ortalaması (32 fit: 4 taban × 4 pencere × 2 ağırlık)

| gaz | **δ(½)** ± σ_konv | (σ_tab / σ_pen / σ_ağ) | **dδ/dτ\|½** ± σ_konv | (σ_tab / σ_pen / σ_ağ) | **b** ± σ_konv | işaret |
|---|---|---|---|---|---|---|
| keskin | **+0.1344** ± 0.0028 | (0.0026/0.0009/0.0001) | **−1.039** ± 0.091 | (0.089/0.015/0.013) | **−3.98** ± 0.48 | 31/32 **neg** |
| A4 | **+0.1374** ± 0.0093 | (0.0073/0.0029/0.0050) | **+0.473** ± 0.186 | (0.154/0.099/0.034) | **+4.05** ± 6.03 | 2/32 neg |
| J14 | **+0.1361** ± 0.0117 | (0.0098/0.0029/0.0057) | **+0.368** ± 0.200 | (0.179/0.089/0.012) | **+5.44** ± 6.27 | 0/32 neg |
| J26 | **+0.1663** ± 0.0071 | (0.0045/0.0037/0.0040) | **+0.164** ± 0.365 | (0.278/0.234/0.036) | **+4.66** ± 3.74 | 1/32 neg |
| N5 | **+0.0671** ± 0.0112 | (0.0091/0.0033/0.0057) | **−0.107** ± 0.175 | (0.121/0.123/0.023) | **+1.19** ± 3.58 | 9/32 neg |
| N5z | **+0.0691** ± 0.0141 | (0.0050/0.0078/0.0106) | **−1.571** ± 0.427 | (0.252/0.344/0.024) | **−2.99** ± 9.07 | 24/32 neg |
| P1 | **+0.0080** ± 0.0275 | (0.0274/0.0021/0.0011) | **−0.842** ± 0.588 | (0.513/0.287/0.007) | **+8.86** ± 5.31 | 0/32 neg |
| **son** | **−0.0685** ± 0.0570 | (**0.0569**/0.0027/0.0016) | **−1.801** ± 0.105 | (0.073/0.060/0.046) | **−5.98** ± 2.04 | **32/32 neg** |
| **orta** | **−0.0698** ± 0.0587 | (**0.0586**/0.0026/0.0031) | **−1.792** ± 0.121 | (0.085/0.040/0.075) | **−6.06** ± 2.17 | **32/32 neg** |

**Gerçek gaz havuzu (son + orta, 64 fit).** İki hata okuması ayrı ayrı
veriliyor: *havuz sd* = 64 fitin doğrudan yayılımı (Mahalanobis
hesaplarında kullanılan); *σ_konv* = §3b'nin kare toplamı (son satırı).

| | ort | havuz sd | bağıl | σ_konv | bağıl | menzil |
|---|---|---|---|---|---|---|
| **δ(½)** | **−0.0692** | 0.0508 | %73.5 | 0.0570 | %82 | −0.1337 … +0.0051 |
| **dδ/dτ\|½** | **−1.7962** | 0.0564 | %3.1 | 0.1049 | %5.8 | −1.880 … −1.612 |
| **b** | **−6.019** | 1.226 | %20.4 | 2.038 | %34 | −8.38 … −1.28 (64/64 neg) |
| a_τ₀ | +10.6925 | 0.0731 | %0.7 | 0.0978 | %0.9 | 10.568 … 10.824 |
| τ₀* | +0.5064 | 0.0047 | %0.9 | 0.0052 | %1.0 | 0.4995 … 0.5124 |

> türetilmiş anahtar sayı: **a_τ₀ − 4π = −1.8738 ± 0.0731**

(159'un aynı sayıyı 158'in 32 fitinden okuduğu hâli: a − 4π = −1.853 ± 0.059.
Fark, buraya eklenen **ağırlık ekseninden** gelir; sadece 158'in alt
kümesinde — 4 taban × W-A × 1/σ_φ, 8 fit — bu koşu **a_τ₀ = 10.7170 ± 0.0670**
ve **b = −7.2229 ± 0.7909** veriyor: 158'in yazdığı 10.717 ± 0.067 ve
−7.2229 ± 0.7909 ile **birebir**.)

### 3c. Türetilmişler ve 158'in ham değerleri

| gaz | a_½ = 4π+dδ/dτ\|½ | a_τ₀ = 4π+dδ/dτ\|τ₀ | τ₀* = ½ − δ(½)/a_½ | 158'in a'sı | 158'in τ₀'ı |
|---|---|---|---|---|---|
| keskin | 11.528 ± 0.091 | 11.621 ± 0.092 | 0.4883 ± 0.0002 | 11.635 | 0.4885 |
| A4 | 13.039 ± 0.186 | 12.961 ± 0.107 | 0.4895 ± 0.0008 | 13.073 | 0.4900 |
| J14 | 12.935 ± 0.200 | 12.832 ± 0.105 | 0.4895 ± 0.0010 | 12.984 | 0.4902 |
| J26 | 12.730 ± 0.365 | 12.632 ± 0.283 | 0.4869 ± 0.0008 | 12.747 | 0.4875 |
| N5 | 12.459 ± 0.175 | 12.448 ± 0.152 | 0.4946 ± 0.0009 | 12.528 | 0.4950 |
| N5z | 10.995 ± 0.427 | 11.042 ± 0.390 | 0.4937 ± 0.0014 | 11.046 | 0.4945 |
| P1 | 11.725 ± 0.588 | 11.731 ± 0.534 | 0.4994 ± 0.0024 | 11.890 | 0.4995 |
| **son** | 10.766 ± 0.105 | 10.691 ± 0.098 | 0.5064 ± 0.0052 | 10.704 | 0.5062 |
| **orta** | 10.775 ± 0.121 | 10.694 ± 0.097 | 0.5065 ± 0.0054 | 10.729 | 0.5063 |

**Kalan farklar tümüyle ağırlık ekseninden gelir** (T2 aynı satırları
1/σ_φ'de birebir veriyor).

---

## 4. T2 — TUTARLILIK SINAVI (eşik 0.001)

W-A, w = 1/σ_φ, 45 (gaz, taban) satırı. **Çapası pencere içinde olan
40 satırda:**

| sınav | sonuç | eşik |
|---|---|---|
| maks \|a_τ₀ − a(158)\| | **0.00049** | 158'in tablosu 3 haneye yuvarlı ⇒ tavan 0.0005 |
| maks \|τ₀* − τ₀(158)\| | **0.00015** | **0.001** ⇒ 0.15× eşik ✓ |

Satır satır en büyük sapmalar: J26 taban 0.34'te Δτ₀ = +0.00015, son taban
0.40'ta −0.00010. ‡ satırları (taban 0.52) dahil edilseydi en büyük sapma
P1'de +0.00068'e, N5z'de +0.00047'ye çıkardı — hâlâ eşiğin altında ama
sistematik olarak kötüleşiyor; **çapanın pencere dışında kalması ölçülebilir
bir bozulmadır.**

---

## 5. T3 — PROBE ANALİZİ

### 5a. Korelasyon merdiveni δ(½)'de

**δ(½) azalan sıra:**
`J26 (+0.166) > A4 (+0.137) > J14 (+0.136) > keskin (+0.134) > N5z (+0.069) > N5 (+0.067) > P1 (+0.008) > son (−0.069) ≈ orta (−0.070)`

**τ₀ artan sıra (155/158):**
`J26 (0.4875) < keskin (0.4885) < A4 (0.4900) < J14 (0.4902) < N5z (0.4945) < N5 (0.4950) < P1 (0.4995) < son (0.5062) ≈ orta (0.5063)`

**Merdiven ayakta ve şekli aynı.** İki sıralama arasındaki tek fark
keskin ↔ A4/J14 üçlüsünün iç sırasıdır (δ(½) merdiveninde keskin dördüncü,
τ₀'da ikinci); ayrım o üçlüde 0.003 mertebesindedir, σ_konv'un (0.003 – 0.012)
içinde. **Sorunun cevabı: EVET — merdiven δ(½)'de "katı/yapısız → korele
itme → GUE-taban → gerçek" sırasıyla görünüyor.**

**Merdiven kesirleri — A4 → gerçek yolunun kaçta kaçı:**

| nicelik | A4 | N5 | P1 | gerçek | **N5 %** | **P1 %** |
|---|---|---|---|---|---|---|
| **δ(½)** | +0.1374 | +0.0671 | +0.0080 | −0.0692 | **34.0** | **62.7** |
| **dδ/dτ\|½** | +0.473 | −0.107 | −0.842 | −1.796 | **25.6** | **57.9** |
| τ₀ (158) | 0.4900 | 0.4950 | 0.4995 | 0.5063 | 30.7 | 58.4 |
| a (158) | 13.073 | 12.528 | 11.890 | 10.717 | 23.1 | 50.2 |
| *155'in τ₀'ı* | | | | | *25* | *52* |

δ(½)'nin merdiveni τ₀'ınkinden **~3 puan önde**, dδ/dτ'ınki a'nınkinden
~3 puan önde (çapa farkı). **155'in "%25 / %52" şekli dört nicelikte de
korunuyor.**

### 5b. dδ/dτ ayrı bilgi taşıyor mu? — 158'in "bağımsız eksenler" bulgusunun yeni hâli

Dokuz gaz üzerinde:

| çift | korel |
|---|---|
| **δ(½) ↔ dδ/dτ\|½** | **+0.785** ← YENİ eksen çifti |
| **τ₀ ↔ a** | **−0.762** ← 158'in eksen çifti |
| δ(½) ↔ τ₀ | **−0.997** |
| dδ/dτ ↔ a | **+0.998** (τ₀ çapasında ÖZDEŞ) |
| δ(½) ↔ a | +0.796 |
| dδ/dτ ↔ τ₀ | −0.747 |
| b ↔ δ(½) | +0.501 |
| b ↔ dδ/dτ | +0.761 |

**Harita — tam-tersinir koordinat değişimi:**

```
dδ/dτ|τ₀ = a − 4π                 (ÖZDEŞ, 1.5e−12)
δ(½)     ≈ a·(½ − τ₀)             (doğrusallaştırılmış)
```

Kontrol (a·(½−τ₀) vs ölçülen δ(½), konvansiyon ortalamaları):

| gaz | δ(½) | a·(½−τ₀) | Δ |
|---|---|---|---|
| keskin | +0.1344 | +0.1341 | +0.0003 |
| A4 | +0.1374 | +0.1311 | +0.0064 |
| J14 | +0.1361 | +0.1272 | +0.0088 |
| N5 | +0.0671 | +0.0630 | +0.0042 |
| N5z | +0.0691 | +0.0605 | +0.0086 |
| P1 | +0.0080 | +0.0062 | +0.0017 |
| son | −0.0685 | −0.0669 | −0.0016 |

> **Hüküm.** dδ/dτ, a'nın **kendisidir** — yeni bilgi taşımıyor, taşıyamaz.
> δ(½), τ₀'ın a ile ölçeklenmiş hâlidir (korel −0.997) — o da yeni bilgi
> taşımıyor. **158'in "a ve τ₀ örtüşen ama özdeş olmayan iki prob" hükmü
> δ-uzayında aynen geçerli, ve eksenler daha bağımsız değil (|r| ≈ 0.78
> her iki koordinat sisteminde de).** Yeni kimliğin kazandırdığı şey bilgi
> değil, **yorum** (omurga ayrıldı) ve **hata bütçesi** (§5d).

### 5c. N5z KARŞI-ÖRNEĞİ — probun en önemli sorusu

158'in en ciddi niteliği: `a(N5z) = 11.046 ± 0.300` vs `a(gerçek) = 10.717 ±
0.075` → **1.1 σ_konv, ayırt edilemiyor**; tek tek fitlerde de örtüşme var.
δ-uzayında aynı sınav:

| eksen | Δ (N5z − gerçek) | σ_bileşik | **ayrım** |
|---|---|---|---|
| dδ/dτ\|½ | +0.225 | 0.259 | **0.9 σ** ← 158'in a'sının aynısı |
| **δ(½)** | **+0.1383** | 0.0512 | **2.7 σ** |
| b | — | — | 0.6 σ |
| **2B (δ(½), dδ/dτ)** | — | — | **2.9 σ** |
| **3B (+b)** | — | — | **2.9 σ** |

**Tek tek fitlerde örtüşme:**

| eksen | N5z menzili (32 fit) | gerçek menzili (64 fit) | örtüşüyor mu? |
|---|---|---|---|
| **δ(½)** | [+0.0551, +0.0789] | [−0.1337, +0.0051] | **HAYIR** |
| dδ/dτ\|½ | [−2.157, −1.178] | [−1.880, −1.612] | EVET |

**Bir-taban-dışarıda dayanıklılık (3B):**

| çıkarılan taban | — | 0.28 | 0.34 | 0.40 | 0.46 |
|---|---|---|---|---|---|
| N5z 3B ayrımı | 2.9 | 5.5 | 2.7 | 2.5 | 3.0 |

**Sürücü ölçüldü — δ(½), taban taban:**

| taban | gerçek-son | gerçek-orta | N5z | \|Δ\| |
|---|---|---|---|---|
| 0.28 | +0.0009 | +0.0051 | +0.0646 | 0.0637 |
| 0.34 | −0.0463 | −0.0507 | +0.0653 | 0.1116 |
| 0.40 | −0.0936 | −0.0936 | +0.0580 | 0.1516 |
| 0.46 | −0.1301 | −0.1312 | +0.0551 | 0.1852 |

> **Hüküm: N5z δ-uzayında AYRIŞIYOR — ama zayıf ve şartlı.**
> Ayrım **her tabanda pozitif** ve tek tek fitler δ(½)'de örtüşmüyor; bu,
> 158'in a-uzayındaki durumundan (örtüşme var) **kesin bir iyileşmedir**.
> Buna karşılık σ_konv birimindeki ayrım **2.5 – 5.5** arasında, yani
> **3σ eşiğinin sınırında**. Ve bu ayrımı taşıyan eksen (δ(½)) tam da
> gerçek gazda taban ekseninde 0.13 yürüyen eksendir (§5d) — yani
> **ayrımın büyüklüğü 155'in bilinen taban sorununa asılıdır**; 158'in
> a için kurabildiği "konvansiyondan bağımsız fark" cümlesi burada
> **kurulamaz**.

### 5d. Konvansiyon bütçesi — ½ çapasının bedeli

| gaz | δ(½): ort ± σ_konv (bağıl) | dδ/dτ\|½: ort ± σ_konv (bağıl) | a_τ₀: ort ± σ_konv (bağıl) |
|---|---|---|---|
| keskin | +0.1344 ± 0.0028 (2.1%) | −1.039 ± 0.091 (8.7%) | 11.621 ± 0.092 (0.8%) |
| A4 | +0.1374 ± 0.0093 (6.8%) | +0.473 ± 0.186 (39.3%) | 12.961 ± 0.107 (0.8%) |
| J14 | +0.1361 ± 0.0117 (8.6%) | +0.368 ± 0.200 (54.4%) | 12.832 ± 0.105 (0.8%) |
| J26 | +0.1663 ± 0.0071 (4.2%) | +0.164 ± 0.365 (222%) | 12.632 ± 0.283 (2.2%) |
| N5 | +0.0671 ± 0.0112 (16.7%) | −0.107 ± 0.175 (163%) | 12.448 ± 0.152 (1.2%) |
| N5z | +0.0691 ± 0.0141 (20.4%) | −1.571 ± 0.427 (27.2%) | 11.042 ± 0.390 (3.5%) |
| P1 | +0.0080 ± 0.0275 (345%) | −0.842 ± 0.588 (69.8%) | 11.731 ± 0.534 (4.6%) |
| **son** | **−0.0685 ± 0.0570 (83.2%)** | −1.801 ± 0.105 (5.8%) | 10.691 ± 0.098 (**0.9%**) |
| **orta** | **−0.0698 ± 0.0587 (84.1%)** | −1.792 ± 0.121 (6.7%) | 10.694 ± 0.097 (0.9%) |

**158 §6f'nin δ-hâli — "gerçek − A4" farkının taban ekseni boyunca kararlılığı:**

| taban | Δδ(½) | Δ(dδ/dτ\|½) | Δa (158) | Δτ₀ (158) |
|---|---|---|---|---|
| 0.28 | −0.1233 | −2.495 | −2.304 | +0.0093 |
| 0.34 | −0.1728 | −2.576 | −2.436 | +0.0139 |
| 0.40 | −0.2283 | −2.417 | −2.439 | +0.0191 |
| 0.46 | −0.2701 | −2.077 | −2.295 | +0.0228 |
| **ortalama** | **−0.1986** | **−2.392** | **−2.369** | **+0.0163** |
| **menzil / ortalama** | **%73.9** | **%20.9** | **%6.1** | **%82.9** |

> **δ(½), τ₀'ın kırılganlığını miras alıyor, a'nın dayanıklılığını değil.**
> ½ çapası kavramsal olarak temiz (gazdan bağımsız bir nokta) ama
> **ampirik olarak üçlünün en oynak ayağıdır**. Türetim hedefi olarak
> **hâlâ dδ/dτ|τ₀ = a − 4π** en temiz sayıdır (%6 taban kararlılığı).

### 5e. b'nin işareti — üçüncü eksen ne ayırıyor?

| gaz | b ortalama | neg/32 | \|ort\|/sd | sign δ(½) | sign dδ/dτ | sign b |
|---|---|---|---|---|---|---|
| keskin | −3.98 | 31 | 2.88 | + | − | **−** |
| A4 | +4.05 | 2 | 1.54 | + | + | + |
| J14 | +5.44 | 0 | 1.67 | + | + | + |
| J26 | +4.66 | 1 | 1.46 | + | + | + |
| N5 | +1.19 | 9 | 0.77 | + | − | + |
| N5z | −2.99 | 24 | 0.60 | + | − | − |
| P1 | +8.86 | 0 | 1.82 | + | − | + |
| **son** | **−5.98** | **32** | **4.62** | **−** | − | **−** |
| **orta** | **−6.06** | **32** | **5.16** | **−** | − | **−** |

**Üçlü işaret sınıfları:**

| sınıf (δ(½), dδ/dτ, b) | gazlar | ne |
|---|---|---|
| (+, +, +) | A4, J14, J26 | **erfc-yumuşatılmış aile** — tek tutarlı aile imzası |
| (+, −, +) | N5, P1 | korele itme (tek doz) + GUE aralık |
| (+, −, −) | keskin, N5z | sert kesim + korele itme ×8 |
| **(−, −, −)** | **son, orta** | **GERÇEK — δ(½) < 0 olan tek gaz** |

**b, 158'in bıraktığı yerde duruyor** (sentetiklerde |ort|/sd ≤ 1.9,
"ölçülemedi"); işareti gerçekte 64/64 negatif (|ort|/sd = 4.6 – 5.2),
P1'de 32/32 pozitif. **Ama ayırıcı olarak b, δ(½) ve dδ/dτ'ın ikisinin de
kaçırdığı P1'i yakalıyor (3.0 σ) — 3B'de P1'in ayrımını ×2.47 büyütüyor.**
Buna karşılık N5z'de (0.6) ve keskin'de (1.1) tek başına işe yaramıyor.

---

## 6. T4 — konvansiyon ekseni

### 6a. Taban öteleme sınavı — 157'nin "%7 salt öteleme"sinin δ hâli

Ortak bantlarda (W-A menzili, τ̄ ∈ [0.43, 0.61]) `δ_taban2 − δ_taban1`:

| gaz | çift | n | ort Δδ | sd | kalan τ-eğimi ± | (eğim×0.18)/öteleme |
|---|---|---|---|---|---|---|
| **son** | 0.28→0.40 | 10 | **−0.1021** | 0.0095 | −0.132 ± 0.030 | **0.23** |
| **son** | 0.34→0.40 | 10 | **−0.0503** | 0.0034 | −0.041 ± 0.013 | **0.15** |
| **son** | 0.40→0.46 | 8 | −0.0285 | 0.0059 | +0.115 ± 0.016 | 0.73 |
| **son** | 0.40→0.52 | 5 | −0.0131 | 0.0150 | +0.457 ± 0.080 | 6.30 |
| orta | 0.28→0.40 | 10 | −0.1051 | 0.0071 | −0.099 ± 0.022 | 0.17 |
| orta | 0.34→0.40 | 10 | −0.0449 | 0.0018 | −0.015 ± 0.009 | **0.06** |
| keskin | 0.28→0.40 | 10 | −0.0039 | 0.0035 | +0.016 ± 0.020 | 0.76 |
| A4 | 0.34→0.40 | 10 | −0.0207 | 0.0509 | −0.543 ± 0.230 | 4.71 |
| N5 | 0.34→0.40 | 10 | −0.0039 | 0.0245 | −0.220 ± 0.121 | 10.27 |
| P1 | 0.28→0.40 | 10 | −0.0450 | 0.0187 | −0.219 ± 0.078 | 0.88 |

> **Gerçek gazda taban konvansiyonu δ'yı BÜYÜK ve NEREDEYSE SALT bir
> ötelemeyle kaydırıyor** (0.28→0.46 toplam öteleme ≈ 0.13; eğim katkısı
> ötelemenin %15 – %73'ü, 0.34→0.40'ta yalnız %6 – %23). 157'nin φ için
> yazdığı "salt öteleme" hükmü δ'da **kalitatif olarak korunuyor ama
> ölçülebilir bir eğim kalıntısıyla** — yani tam değil.
>
> **Sentetiklerde ötelemeler 10 – 30 kat küçüktür** (0.003 – 0.06) ve
> oranlar (4 – 20) bu küçüklüğün yapaylığıdır, eğimin büyüklüğünün değil.
>
> **δ(½)'nin konvansiyon bütçesi budur:** gerçekte 0.13 (gerçek↔A4
> sinyalinin **%65'i**), sentetiklerde ≤ 0.02.

### 6a2. δ(½)'nin taban yürüyüşü, gaz gaz — dördüncü eksen adayı

W-A, w = 1/σ_φ:

| gaz | 0.28 | 0.34 | 0.40 | 0.46 | **menzil** | gerçeğin kaçta biri |
|---|---|---|---|---|---|---|
| **son** | +0.0009 | −0.0463 | −0.0936 | −0.1301 | **0.1310** | — |
| **orta** | +0.0051 | −0.0507 | −0.0936 | −0.1312 | **0.1363** | — |
| P1 | +0.0359 | +0.0171 | +0.0050 | −0.0293 | 0.0652 | 1/2.0 |
| J14 | +0.1156 | +0.1249 | +0.1333 | +0.1377 | 0.0221 | 1/6.0 |
| N5 | +0.0518 | +0.0602 | +0.0670 | +0.0730 | 0.0213 | 1/6.3 |
| A4 | +0.1243 | +0.1265 | +0.1347 | +0.1400 | 0.0157 | 1/8.5 |
| J26 | +0.1547 | +0.1610 | +0.1601 | +0.1656 | 0.0110 | 1/12.2 |
| N5z | +0.0646 | +0.0653 | +0.0580 | +0.0551 | 0.0102 | 1/13.1 |
| keskin | +0.1367 | +0.1347 | +0.1329 | +0.1306 | **0.0061** | 1/22.0 |

**Gerçek gaz, taban eksenine bütün sentetiklerden daha duyarlıdır** —
P1'in 2 katı, geri kalan altısının 6 – 22 katı. Ayrıca **işaret bile
gaza bağlı**: gerçek, keskin, J26 ve N5z tabanla aşağı yürüyor;
A4, J14, N5 yukarı. Bu, §5d'de "bütçe" diye sayılan şeyin kendi başına
bir yapı taşıdığının kaydıdır. **Bir eksen olarak sınanmadı.**

### 6b. L-değişmezlik (son L = 12.0296 ↔ orta L = 11.4638)

| nicelik | son | orta | \|Δ\| | \|Δ\|/σ_konv |
|---|---|---|---|---|
| δ(½) | −0.0685 | −0.0698 | 0.0013 | **0.02** |
| dδ/dτ\|½ | −1.8006 | −1.7918 | 0.0087 | **0.08** |
| b | −5.980 | −6.057 | 0.077 | **0.04** |
| a_τ₀ | 10.6912 | 10.6939 | 0.0028 | **0.03** |
| τ₀* | 0.50637 | 0.50647 | 0.0001 | **0.02** |

**Üçlünün üç ayağı da L'den bağımsız** — 157/158'in L-değişmezlik hükmü
δ-uzayında tam olarak korunuyor.

### 6c. Izgara ekseni (155'in `genis` ızgarası vs 158'inki, taban 0.40)

| gaz | δ(½): 158 → 155 | Δ | dδ/dτ: 158 → 155 | Δ |
|---|---|---|---|---|
| keskin | +0.1329 → +0.1349 | +0.0020 | −0.993 → −0.993 | −0.000 |
| A4 | +0.1347 → +0.1380 | +0.0033 | +0.602 → +0.591 | −0.011 |
| N5z | +0.0580 → +0.0700 | +0.0119 | −1.178 → −1.261 | −0.083 |
| P1 | +0.0050 → +0.0066 | +0.0016 | −0.573 → −0.688 | −0.115 |
| **son** | −0.0936 → −0.0976 | −0.0040 | −1.815 → −1.866 | −0.051 |
| **orta** | −0.0936 → −0.0967 | −0.0031 | −1.774 → −1.833 | −0.059 |

Izgara ekseni δ(½)'de ≤ 0.012, dδ/dτ'da ≤ 0.115 — **taban ekseninin
(0.13 / 0.25) onda biri kadar.** Bütçeye eklenmedi (155'in ızgarası bütün
gazlarda mevcut değil), ama kaydediliyor.

---

## 7. T5 — HÜKÜM ARİTMETİĞİ: hangi birleşim gerçeği tüm sentetiklerden ayırıyor?

### 7a. Bütün eksen kombinasyonları (Mahalanobis, σ_konv birimi)

| kombinasyon | keskin | A4 | J14 | J26 | N5 | N5z | P1 | **TABAN** |
|---|---|---|---|---|---|---|---|---|
| δ(½) | 4.0 | 4.0 | 4.0 | 4.6 | 2.7 | 2.7 | **1.4** | **1.4** (P1) |
| dδ/dτ\|½ | 6.9 | 9.5 | 7.8 | 5.2 | 10.3 | **0.9** | 2.2 | **0.9** (N5z) |
| b | 1.1 | 3.5 | 3.3 | 3.1 | 3.7 | **0.6** | 3.0 | **0.6** (N5z) |
| δ(½) + dδ/dτ | 7.9 | 10.4 | 8.8 | 6.9 | 10.6 | 2.9 | **2.2** | **2.2** (P1) |
| δ(½) + b | 4.3 | 5.7 | 5.6 | 5.9 | 4.8 | **2.8** | 4.0 | **2.8** (N5z) |
| dδ/dτ + b | 8.8 | 9.8 | 8.4 | 6.7 | 10.7 | **0.9** | 7.1 | **0.9** (N5z) |
| **ÜÇLÜ** | **9.9** | **10.8** | **9.6** | **8.3** | **11.1** | **2.9** | **7.3** | **2.9** (N5z) |

*(karşılaştırma: 158'in tek ekseni a — en zayıf ayrım **1.1 σ_konv**, N5z.)*

### 7b. Tamamlayıcılık — üçlünün gücü nereden geliyor

| gaz | δ(½) | dδ/dτ\|½ | b | en iyi tek | 3B | **kazanç** |
|---|---|---|---|---|---|---|
| keskin | 4.0 | **6.9** | 1.1 | dδ/dτ | 9.9 | ×1.43 |
| A4 | 4.0 | **9.5** | 3.5 | dδ/dτ | 10.8 | ×1.13 |
| J14 | 4.0 | **7.8** | 3.3 | dδ/dτ | 9.6 | ×1.23 |
| J26 | 4.6 | **5.2** | 3.1 | dδ/dτ | 8.3 | ×1.59 |
| N5 | 2.7 | **10.3** | 3.7 | dδ/dτ | 11.1 | ×1.08 |
| **N5z** | **2.7** | 0.9 | 0.6 | **δ(½)** | 2.9 | ×1.08 |
| **P1** | 1.4 | 2.2 | **3.0** | **b** | 7.3 | **×2.47** |

> **Üçlünün gücü ortalamadan değil, ORTOGONAL BAŞARISIZLIKLARDAN geliyor.**
> dδ/dτ (yani a) beş gazı güçlü ayırıyor ama N5z'de ve P1'de düşüyor;
> tam o iki gazı **δ(½)** ve **b** yakalıyor. Bu, 158'in "a ile τ₀
> birbirinin yerine geçmez, birlikte kullanılmalıdır" hükmünün **nicel
> kapanışıdır**.

---

# HÜKÜM

## (i) Rotor düzeltmesi yeni bir ölçüm değil, tam bir koordinat değişimidir

Aynı ağırlıkla `c_φ = c_δ + (0, 4π, −2π)` (1.6e−10). Bunun üç sonucu
ölçüldü ve hepsi tam:

* **a_τ₀ ≡ a(158)** (1.5e−12) — a "yeniden ölçülmedi", **aynı sayıdır**;
* **b(δ) ≡ b(φ)** — **b, rotor düzeltmesinin dokunmadığı tek sayıdır**;
  159'un gördüğü 0.17 – 2.42'lik fark **ağırlık konvansiyonundan** gelir;
* **τ₀* = ½ − δ(½)/a_½** 158'in τ₀'ını **0.00015** ile veriyor (eşik 0.001).

**Yani "fizik artık δ'da" cümlesi bir yorum kazancıdır, bir bilgi kazancı
değil.** δ(½) ile τ₀ arasında korel −0.997; dδ/dτ ile a arasında özdeşlik.
155/157/158'in bütün sayıları yeniden koşulmadan, tek bir cebirsel
dönüşümle δ-uzayına taşınabilir — ve bu raporda taşındı.

## (ii) Tek gerçek KAZANÇ: istatistik hata ve üçlünün tamamlayıcılığı

1. **σ_δ < σ_φ.** Medyan oran gerçek gazda **0.38 (son) / 0.20 (orta)**,
   sentetiklerde 0.60 – 0.95. φ ile τ_eff birlikte dalgalandığı için omurga
   çıkınca ortak gürültü de çıkıyor. **δ, φ'den 2.6 – 5 kat daha hassas
   ölçülen bir gözlenebilirdir** (gerçek gazda).
2. **Üçlü (δ(½), dδ/dτ|½, b) gerçeği YEDİ SENTETİĞİN HEPSİNDEN ayırıyor:**
   en zayıf ayrım **2.9 σ_konv** (N5z), en güçlüsü 11.1 (N5). 158'in tek
   ekseni a bunu yapamıyordu (N5z 1.1σ). **Hiçbir tek eksen ve hiçbir ikili
   de yapamıyor** (en iyi ikili δ(½)+b, taban 2.8). Üçlü **gereklidir**.

## (iii) En ciddi nitelendirme: δ(½) kırılgan, ve ayrımı ona asılı

δ(½)'nin gerçek gazdaki konvansiyon bütçesi **σ_konv = 0.057** ve bunun
**0.0569'u taban eksenidir** (bağıl %83). "Gerçek − A4" farkının taban
boyunca menzil/ortalaması **%74** — τ₀'ınkiyle (%83) aynı mertebe, a'nınkinin
(%6.1) **on iki katı**. Ve **N5z'yi yakalayan eksen tam olarak budur.**

> Kurulabilecek cümle: **"parmak-izi üçlüsü gerçeği yedi sentetiğin
> hepsinden ≥ 2.9 σ_konv ile ayırır."**
> Kurulamayacak cümle: *"bu ayrım konvansiyondan bağımsızdır."* Değildir —
> ayrımın büyüklüğü 155'in bilinen taban sorununu miras alıyor. 158'in a
> için kurabildiği "fark taban boyunca %6 içinde sabit" cümlesinin
> δ(½) karşılığı **%74**'tür.

**½ çapası kavramsal olarak temiz, ampirik olarak kırılgandır.** Türetim
hedefi olarak en temiz sayı hâlâ **a − 4π = dδ/dτ|τ₀ = −1.874 ± 0.073
(gerçek, 64 fit)**'tir — 159'un hükmü ayakta.

## (iv) Parmak izi üçlüsünün son hâli

> **GERÇEK GAZ** (son + orta, 4 taban × 4 pencere × 2 ağırlık = 64 fit;
> hata = havuz sd, parantezde σ_konv):
>
> | | değer | ± sd (σ_konv) | rol |
> |---|---|---|---|
> | **δ(½)** | **−0.0692** | 0.0508 (0.057) | **N5z'yi yakalayan eksen**; τ₀'ın a-ölçekli hâli; en kırılgan ayak |
> | **dδ/dτ\|½** | **−1.7962** | 0.0564 (0.105) | a − 4π'nin ½-çapalı hâli; **beş gazı** yakalıyor |
> | **b** | **−6.019** | 1.226 (2.038) | rotorun **dokunmadığı** sayı; **P1'i** yakalıyor |
>
> türetilmiş: **a_τ₀ − 4π = −1.8738 ± 0.0731** · **τ₀* = 0.5064 ± 0.0047**
> işaret imzası: **(−, −, −)** — dokuz gaz içinde **yalnız gerçek**.

Sentetiklerin tam tablosu §3b'de. Üç eksenin ayrı ayrı yakaladığı gazlar
§7b'de.

## Sıradaki adım (bu ölçümün işaret ettiği)

1. **δ(½)'nin TABAN DUYARLILIĞI dördüncü bir eksen adayıdır — ve ölçüldü
   (§6a2), ama bir eksen olarak SINANMADI.** δ(½)'nin taban 0.28→0.46
   menzili (W-A, w = 1/σ_φ):

   | gerçek | P1 | J14 | N5 | A4 | J26 | N5z | keskin |
   |---|---|---|---|---|---|---|---|
   | **0.134** | 0.065 | 0.022 | 0.021 | 0.016 | 0.011 | 0.010 | 0.006 |

   Gerçek gaz, P1'in **2 katı**, geri kalan altısının **6 – 22 katı**
   kadar taban-duyarlıdır. Yani §5d'de "kirlilik" diye kaydedilen şey
   kendi başına bir sinyal olabilir. Sınanması ucuz (bu koşunun
   çıktısında hazır), ama bir hata bütçesi kurmak yeni bir konvansiyon
   ekseni gerektirir (taban duyarlılığının duyarlılığı) — bu koşu onu
   kurmadı.
2. **N5z'yi 3σ'ın üstüne çıkarmanın yolu taban ekseninden geçiyor.**
   δ(½)'nin taban yürüyüşünü fiziksel olarak anlamak (η regresyonunun
   taban altındaki çizgileri çıkarmasının δ'ya etkisi) bütçeyi
   küçültebilir. Bu, 155'ten beri açık olan tek sistematiktir.
3. **b'nin sentetiklerdeki ölçülemezliği hâlâ kapalı** (158 §7). Ama
   161 yeni bir kayıt ekliyor: b, **ağırlık konvansiyonuna** çok duyarlı
   (N5z'de 1/σ_φ → 1/σ_δ geçişi b'yi +5.76'dan −1.59'a taşıyor). Yani
   b'nin "ölçülemezliği" yalnız taban ekseninden değil, ağırlık ekseninden
   de besleniyor — ve bu ikincisi 158'de ölçülmemişti.
4. **Bu raporun hiçbir sayısı yeni ölçüm gerektirmedi.** Aynı 45 koşu,
   `dusuk` ve `P0` dışındaki bütün gazlar için tam bir δ haritası veriyor.
   `dusuk` gazı 158'de koşulmadığı için burada da yok — koşulursa
   üçlü tablosuna doğrudan eklenebilir.

---

## Dürüstlük notları

* **TABLODAN-OKUNAN vs YENİDEN-HESAPLANAN, satır satır.**
  * *Yeniden hesaplanan (bu raporun bütün birincil sayıları):* δ tayfları,
    δ(½), dδ/dτ, b, a_½, a_τ₀, τ₀*, σ_δ, bütün konvansiyon bütçeleri,
    bütün Mahalanobis ayrımları. Kaynak: `scratchpad/158/F_*.json`'un
    bant düzeyi kayıtları (`tau_eff`, `phi`, `sPhi_jk`, `phi_jk`, `teff_jk`).
  * *Tablodan okunan (yalnız referans, hükme girmiyor):* 158 §4'ün W-A
    (τ₀, a, b) tablosu (45 satır, `161_cekirdek.REF158_WA`), 158 §6a'nın
    σ_konv sütunu, 159 §2c'nin δ tablosu ve §2d'nin beş gazlık özeti.
    Hepsi V1/V3/V4'te bu koşunun kendi hesabıyla karşılaştırıldı ve
    tutarlı çıktı.
  * **Yeniden ölçüm YOK.** Ne η zinciri, ne gaz üreticileri, ne çizgi
    döngüsü koşuldu. Bu, 158/159'dan farklı bir çalışma tarzıdır ve
    bilinçlidir: δ, φ'nin cebirsel bir dönüşümüdür; yeniden ölçmek aynı
    sayıyı verirdi (159 V1 bunu bit düzeyinde göstermişti).
* **V1 denetimi 45 satırın 45'inde geçiyor:** maks |Δτ₀| = 0.00005,
  |Δa| = 0.00049, |Δb| = 0.00048 — 158'in tablosu 3 haneye yuvarlı olduğu
  için bu **yuvarlama tavanıdır**, bir sapma değil.
* **159'un δ tablosuyla farklar açıklandı, gizlenmedi.** 159 δ'yı 1/σ_δ
  ile tartmıştı; bu koşu her iki ağırlığı da koşuyor. **1/σ_δ'de** bu
  raporun sayıları 159'unkilere δ(½)'de **≤ 0.006** (son 0.0003, orta
  0.0014, keskin 0.0024, A4 0.0053; **P1 istisna: 0.0121**), dδ/dτ'da
  **≤ 0.010** yaklaşıyor. **1/σ_φ'de** fark daha büyük (δ(½)'de 0.0019 –
  0.0137) ve **beklenen**: farklı ağırlık, farklı fit. Kalan mikro-farklar
  159'un tablosunun 3–4 haneye yuvarlı olmasındandır.
* **Taban 0.52 satırları bir ölçüm değildir (‡)** ve bütün ortalamalardan
  çıkarıldı. Gerekçesi mekanik: çapa τ = ½, W-A'nın en alt bandı τ̄ = 0.53
  → δ(½) bir ekstrapolasyondur. Satırlar T2'de duruyor; dahil edilselerdi
  τ₀* sapması 0.00015'ten 0.00068'e çıkardı.
* **355 / 360 fit kuruldu.** Kurulamayan beşi taban 0.52 + W-B
  (A4, P1, keskin — sıfır-merkezli pencere o tabanda yeterli bant
  bulamıyor). Hepsi zaten ‡'dir.
* **σ_konv'un üç bileşeni birbirinden bağımsız varsayılarak toplandı**
  (kare toplam). Bağımsız olmayabilirler; bu bir modeldir, ölçüm değil.
  Bileşenler tabloda ayrı ayrı duruyor (§3b), okuyucu kendi birleştirmesini
  yapabilir.
* **Mahalanobis ayrımları KONVANSİYON bulutunun kovaryansıyla hesaplandı,
  istatistik hatayla değil.** Bulut 32 (sentetik) / 64 (gerçek) noktadır ve
  noktalar bağımsız değildir (aynı ham koşuların farklı okumaları). Yani
  bu sayılar bir olasılık değil, **konvansiyon bütçesi biriminde bir
  mesafedir** — 158 §6a'nın "kaç σ_konv" ölçütüyle aynı cinsten ve aynı
  sınırlamalarla.
* **N5z'nin 2.9 σ'sı 3σ'ın altındadır ve öyle yazıldı.** Bir-taban-dışarıda
  sınavında 2.5 – 5.5; en küçüğü taban 0.40 çıkarılınca. **Mühür YOK.**
* **P1'in δ(½)'si sıfıra çok yakın** (+0.008), bu yüzden bağıl bütçesi
  %345 çıkıyor. Bu bir ölçüm bozukluğu değil, sıfıra bölmenin yapaylığıdır;
  mutlak bütçe 0.0275 ile sentetiklerin en büyüğü ama gerçeğinkinden
  (0.057) küçüktür.
* **b hakkında 158'in "ölçülemedi" hükmü DEĞİŞMEDİ.** 161 yalnız b'nin
  (a) rotor düzeltmesinden etkilenmediğini ve (b) ağırlık konvansiyonuna
  duyarlı olduğunu ekliyor. b'nin ayırıcı gücü (P1'de 3.0σ) bir
  **kombinasyon** ifadesidir, tek başına bir ölçüm iddiası değil.
* **`dusuk` gazı ve P0 bu haritada yok.** `dusuk` 158'de koşulmadı,
  P0 158'de "ÖLÇÜLEMEDİ" idi. İkisi de hiçbir hükme girmedi.
* **Süreler ölçüt değildi.** Analiz ~40 s, doğrulama ~35 s, figür ~15 s;
  hepsi tek makinede, yeniden ölçüm olmadığı için.

---

## Ek — dosyalar ve tekrar-üretim

| dosya | ne |
|---|---|
| `161_configs/161_cekirdek.py` | yükleyici + `delta_bant` (δ ve σ_δ) + `_fit`/`secim`/`olc` (158'in fit konvansiyonu + ½ çapası + iki ağırlık) + referans tabloları |
| `161_configs/161_dogrulama.py` | V1 (158'in W-A tablosu, 45 satır) · V2 (özdeşlik c_φ − c_δ) · V3 (159'un bant tablosu) · V4 (159 §2d) · V5 (σ_δ vs σ_φ) · V6 (çapa farkı) |
| `161_configs/161_analiz.py` | T0 (tayflar) · T1 (parmak izi) · T2 (tutarlılık) · T3 (probe) · T4 (konvansiyon) · T5 (hüküm aritmetiği) |
| `161_configs/161_figur.py` | 6 panelli figür |
| `161_delta_harita.png` | figür |

Ham çıktılar `scratchpad/161/`: `analiz_cikti.txt` (431 satır),
`dogrulama_cikti.txt` (121 satır), `ozet.json` (gaz başına parmak izi +
2B/3B Mahalanobis). Girdi: `scratchpad/158/F_*.json` (**değiştirilmedi**),
`scratchpad/155/tau0_*_genis.json` (yalnız okundu). Git'e dokunulmadı.

**Denetim özeti:**

| | sınav | sonuç |
|---|---|---|
| **V1** | 161'in fiti 158'in W-A tablosunu üretiyor mu? (45 satır) | maks \|Δτ₀\| = **0.00005**, \|Δa\| = **0.00049**, \|Δb\| = **0.00048** (158 3 haneye yuvarlı) |
| **V2** | `c_φ − c_δ = (0, 4π, −2π)` özdeşliği (178 fit) | maks \|Δkatsayı\| **1.6e−10** · \|b(δ)−b(φ)\| **3.1e−12** · \|a_τ₀−a(φ)\| **1.5e−12** |
| **V3** | 159 §2c bant tablosu (gerçek, taban 0.40) | maks \|Δδ\| = **0.00080** (159 4 haneye yuvarlı) |
| **V4** | 159 §2d beş gazlık özeti | §1c'de tablolandı; fark = ağırlık konvansiyonu |
| **V5** | σ_δ vs σ_φ | medyan oran **0.197 – 0.949**; gerçek 0.384/0.197 |
| **V6** | Çapa farkı = 2b(½−τ₀) mü? | 9 gazda maks sapma **9.3e−16** (ÖZDEŞ) |

Tekrar üretmek için (yeniden ölçüm gerekmez):

```
.venv/bin/python qm_riemann/161_configs/161_dogrulama.py
.venv/bin/python qm_riemann/161_configs/161_analiz.py
.venv/bin/python qm_riemann/161_configs/161_figur.py
```
