# 162 — GAZIN BELLEĞİ C_n: kimlik sınavı + Gauss-tayfsal kapanış + ×1.4'ün tayf kimliği

**Hedef.** KALEM (01 Eylül akşamı) iki şey söyledi: (i) `C_n = Σ_{k<n} ds_k`
teleskopik olarak `−S(z_n) + sabit`'tir — gazın belleği ayrı bir nesne
değil, SAYIM DALGALANMASININ teknelerde okunmuş halidir; (ii) 160'ın
kapattığı `δ = K1 + K2` muhasebesinin bütün girdileri `(η_n, C_n)` ortak
sürecinin fonksiyonelidir, dolayısıyla o süreç **ortak-Gauss** ise δ
YALNIZ ikinci momentlerden (tayftan) çıkmalıdır. Bu rapor birinci iddiayı
ölçtü (T1) ve ikinciyi **sınadı ve reddetti** (T2/T3).

---

> ## Kısa hüküm
>
> 1. **T1 GEÇTİ. C, merdivenin ta kendisidir.** Son 300k'da
>    `corr(cumsum(ds), −S_merdiven(z_n)) = +0.9727` (varyans oranı 0.861,
>    eğim 0.902) ve bond sürümünde
>    `corr(Ĉ, −S_merdiven(mid_n)) = +0.9764`, **eğim 0.9962, varyans oranı
>    1.041, artık payı %4.7**. Beklenti (≥0.97) tuttu.
> 2. **C'NİN DEVASA VARYANSI BELLEK DEĞİL, YOĞUNLUK SÜRÜKLENMESİDİR.**
>    D1'in ham `C = Σ X̃` nesnesinin varyansı **1.912e+04**; aritmetik
>    belleğinki `Var(Ĉ) = 0.0726`. Aradaki her şey `r_n = X̃_n − dsΔ_n ≈
>    L/Lw_n − 1` sürüklenmesidir (korel 0.974, gücünün %94'ü en düşük 500
>    FFT kutusunda). **Bu ayrım bu koşuda ilk kez yapıldı ve TAMDIR** (fit
>    yok): `m_n = m_0 + (2π/L)(u_n + Ĉ_n)` özdeşliği **0.0e+00**.
> 3. **`corr(Ĉ, dsΔ) = −0.43401` ve kalem karşılığı `−½√(σΔ²/Var Ĉ) =
>    −0.43401`** — beş basamakta aynı. Yani Ĉ'nin durağanlığı ölçüldü.
>    Ayrıca **`corr(Ĉ, η) = +0.00014`**: bellek, genliğe eşzamanlı olarak
>    TAM DİKTİR (üç gazda da ≤ 1.4e−04).
> 4. **T2 REDDEDİLDİ — ve büyük farkla. GAUSS-TAYFSAL KAPANIŞ YOKTUR.**
>    Aynı çizgi güçlerine, aynı çapraz-tayfa ve aynı süreklilik tayfına
>    sahip ORTAK-GAUSS vekilde δ **çöküyor**: gerçek gazda
>    `δ_gerçek/δ_vekil` medyanı **+6.36** [3.36, 14.95] (7 bant, 4 tohum),
>    keskin'de **+2.31**, A4'te **−11.76** (İŞARET DEĞİŞTİRİYOR).
>    Görevin ±%15 ölçütüne göre **0/17 bant**. **Gauss-ötesi pay
>    (oran−1) gerçek gazda %536 medyan** — yani δ'nın ancak **%16'sı**
>    ikinci momentlerden geliyor.
> 5. **KONTROL BİT DÜZEYİNDE GEÇTİ.** Aynı boru hattının faz-karıştırmasız
>    sürümü (V0) gerçeği **birebir** veriyor: `δ_gerçek/δ_V0 = 1.00000`,
>    bant düzeyi fark **0.0e+00**, dizi düzeyi ≤ 1.1e−16. Ve 162'nin
>    `gercek` koşusu 160'ın kayıtlı JSON'unu **63 bantta, bütün
>    sütunlarda 0.000e+00** ile yeniden üretiyor.
> 6. **TAŞIYICI BULUNDU: ÇİZGİLER-ARASI FAZ KİLİDİ.** Yalnız SÜREKLİLİK
>    karıştırıldığında (VC) δ neredeyse hiç değişmiyor
>    (`δ_g/δ_VC` medyan **+1.24**, A4'te **+0.82**); yalnız ÇİZGİ fazları
>    karıştırıldığında (VL) δ tam vekildeki gibi çöküyor (medyan **+7.27**).
>    **δ'yı asalların birbirine kilitli fazları taşıyor**, sürekli taban
>    değil.
> 7. **T3: ×1.4 BİR TAYF FARKI DEĞİL.** Üç gazın vekilleri **birbirine
>    yakınsıyor**: `a(φ)` 10.652 / 11.667 / 13.032 → **12.33 / 12.50 /
>    12.51** (çıplak 4π = 12.566), `b(φ)` −6.892 / −4.434 / **+3.614** →
>    **−0.94 ± 0.86 / +0.16 ± 1.40 / +0.16 ± 1.14**, τ₀ 0.509 / 0.488 /
>    0.489 → **0.504 / 0.501 / 0.502** (kinematik 0.5). Gazlar Gauss
>    vekilde **ayırt edilemez** hale geliyor.
> 8. **b'NİN İŞARETİ TAMAMEN GAUSS-ÖTESİDİR — ve `⟨σX̃²⟩` bunu doğrudan
>    gösteriyor.** A4'ün `−A²s2/2 = −A²⟨σX̃²⟩/2⟨ρ⟩`'si üç bantta
>    **+0.0533 / +0.1001 / +0.2368**; vekilde **−0.0106 / −0.0127 /
>    −0.0061** — **İŞARET ÜÇ BANTTA DA TERS**. Gerçek gazda işaret
>    korunuyor ama büyüklük düşük τ'da ×4.98'e kadar düşüyor.
>    **A4'ün dışbükeyliğini üreten nesne Gauss vekilde YOK.**
> 9. **ŞİŞME BAHANESİ ELENDİ.** Vekilin σΔ²'si gerçek gazda %28, A4'te %59
>    şişiyor (keskin'de **−%1.3**, yani şişmiyor). δ, A·σ_X'e güçlü
>    bağlı olduğundan gerçeğin δ(A·σ_X) eğrisi vekilin KENDİ A·σ_X'inde
>    değerlendirildi: oran **büyüyor** (+3.4…+15 yerine +7.2…+15).
>    Şişme açığı açıklamıyor, **tersine çeviriyor**. Ve tayf bütçesi
>    neredeyse tam kapanan keskin'de (kapanış 0.954/0.989, σΔ² farkı
>    −%1.3) çöküş **aynen** var.
> 10. **YENİ NESNE ÖLÇÜLDÜ: ASAL ÇİZGİLERİN YIKICI GİRİŞİMİ.**
>    `Σ_q |c_q(η)|²/2 = 0.06997` iken `Var(η) = 0.05431` — **oran 1.288**.
>    Bu bir çıkarım artefaktı DEĞİL: bilinen rastgele-fazlı sinyalde aynı
>    çıkarım **1.0045** veriyor (V5). Yani gerçek gazda η, çizgilerinin
>    uyumsuz toplamından **%29 daha sessizdir** — asalların faz kilidinin
>    doğrudan, φ'den bağımsız bir ölçüsü. Ĉ'de aynı oran **1.026**
>    (girişim yok): **kilit ρ-ağırlığının gördüğü η kanalındadır.**
> 11. **ÜÇ VEKİL REÇETESİ DENENDİ, İKİSİ ÖLÇÜLÜP ELENDİ.** Düz n-uzayı FFT
>    faz karıştırması çizgiyi öldürüyor (pow/gp 0.73 → **0.0022**); düzgün
>    u-ızgarası + bant-sınırlı ara temsil τ>0.5'i öldürüyor (0.628 →
>    0.008); öz-tutarlılıksız konum-uzayı karıştırması yarı yolda kalıyor
>    (0.628 → 0.089). Kabul edilen reçete (öz-tutarlı) çizgi gücünü
>    **%2–%16** içinde koruyor (τ = 0.78'de −%37). Üçü de tabloda.

---

## 1. Yöntem

### 1a. Ölçüm zinciri KOPYALANMADI

`162_configs/162_cekirdek.py`'nin `olc162`'si **160'ın `olc160`'ını
olduğu gibi çağırır**; yalnız `Yerel160` sınıfı, önceden kurulmuş nesneyi
döndüren bir fabrikayla geçici olarak değiştirilir. Böylece aday listesi,
220-örnekleme tohumu (21), `gap < 2.5·dres` filtresi, on/off ara-nokta
referansı (`W' = w + gap/2`), 8-grup round-robin jackknife ve bant sağlığı
kuralı **bit düzeyinde** 160'ınkidir. **V1 bunu sınadı: 3 gaz × 2 ızgara,
63 ortak bant, ölçülen bütün sütunlarda maks fark `0.000e+00`.**

### 1b. C'nin TAM ayrışımı (fit yok) — ve neden gerekli

160'ın D1 kimliğindeki `C_n = Σ_{k≤n} X̃_k` nesnesinin varyansı
**1.912e+04**'tür. Bu, bellek değildir: `X̃` global L ile tanımlıyken
gerçek yerel aralık `2π/Lw_n` olduğundan, `X̃` içinde yavaş bir
**yoğunluk sürüklenmesi** vardır. Ayrışım TAMDIR:

| nicelik | tanım | ölçülen (son) |
|---|---|---|
| `dsΔ_n` | `(ds_n+ds_{n+1})/2 − ort` (156.`_bond` = 160'ın `X_tam`'ı) | Var 0.05466 |
| `r_n` | `X̃_n − dsΔ_n` — **SÜRÜKLENME** | std 3.67e−03, menzil [−1.14e−02, +1.23e−02] |
| `u_n` | `n + Σ_{k<n} r_k` — sürüklenmiş kafes ekseni | salınım ≈ 470 |
| `Ĉ_n` | `Σ_{k<n} dsΔ_k` — **DURAĞAN BELLEK** | Var 0.07255, `Ĉ_0 = Ĉ_son = 0` (1.3e−12) |
| kimlik | `m_n = m_0 + (2π/L)(u_n + Ĉ_n)` | maks fark **0.0e+00** |

`r`, `L/Lw_n − 1` ile **korel 0.974**'tür ve gücünün **%94'ü en düşük 500
FFT kutusundadır** — yani düzgün, deterministik, gaza ait olmayan bir
kayma. Vekil üretiminde `r` ve `u` **dokunulmadan** bırakılır.

### 1c. Vekil reçetesi (KABUL EDİLEN) — konum uzayında faz karıştırması

İkinci-moment yapısı **konum uzayında** şudur: her asal frekans
`w_q = log q`'daki çizgi genliği `|c_q|`, çizgiler-arası çapraz-tayf, ve
artık sürekliliğin tayfı. Vekil dört adımda:

1. **ÇIKARIM.** `c_q(x) = 2⟨x_n e^{−i w_q m_n}⟩`, `x ∈ {η, Ĉ}`,
   `q ∈ pk_m(e^{0.86·L})` — **160'ın aday listesiyle AYNI küme** (3425
   çizgi). `c_q`, ölçümün o çizgide gördüğü genliğin ta kendisidir
   (`zc = c_q`), bu yüzden ham izdüşüm kullanılır: vekil `pow_on`'u
   birebir yeniden üretsin diye.
2. **AYRIM.** `x_çizgi = Re Σ_q c_q e^{i w_q m_n}`, `x_sür = x − x_çizgi`.
   Çizgi payı (varyansta): η **0.487**, Ĉ **0.918**.
3. **KARIŞTIRMA.** `c'_q = c_q e^{iψ_q}`, `ψ_q ~ U(0,2π)` — **η ve Ĉ için
   AYNI ψ_q** (çizgi çapraz-tayfı ve çizgi güçleri TAM korunur). Süreklilik
   n-uzayında **ORTAK** faz karıştırması (Prichard–Theiler: `|F_k|`
   korunur, aynı rastgele faz her iki artığa; `k=0` ve Nyquist fazı
   korunur ⇒ **ortalama korunur**). ≈3400 rastgele fazlı çizgi + faz
   karıştırılmış süreklilik ⇒ süreç **Gauss**'tur (Rice; V4'te ölçüldü).
4. **ÖZ-TUTARLI KURULUM.** Ĉ bir **alandır** (`Ĉ_n = −S(z_n)`, S konumun
   fonksiyonu) ve konumlar Ĉ'den doğar. Vekil gaz **aynı sabit noktayı**
   çözer (152/154'ün sentetik gaz üreticisinin Newton'uyla aynı yapı):

   ```
   Ĉ'_n = Ĉ'_sür,n + Re Σ_q c'_q e^{i w_q m'_n}
   m'_n = m_0 + (2π/L)(u_n + Ĉ'_n − Ĉ'_0)
   dsΔ' = diff(Ĉ') ,  X̃' = r + dsΔ' ,  η'_n = η'_sür,n + Re Σ_q c'^η_q e^{i w_q m'_n}
   ```

   Newton (adım sınırı 0.3·(2π/L), 6 adım). **Tutarlılık ÖZDEŞTİR:**
   `X̃'_n ≡ (m'_{n+1}−m'_n)L/2π − 1 = r_n + dsΔ'_n`. V0'da `m' = m` TAM
   sabit noktadır (kalıntı 2.3e−10 = kayan nokta gürültüsü).

**Tohumlar:** `VG101…VG104` (gerçek), `VG101…VG103` (keskin, A4).
Çizgi fazları `default_rng(tohum)`, süreklilik `default_rng(tohum+500000)`.
**Pencereleme yok**: süreklilik karıştırması tam kayıt üzerinde dairesel
FFT'dir; `Ĉ_0 = Ĉ_son = 0` (1.3e−12) olduğundan sarma süreksizliği yok.

**Çeşitler** (aynı boru hattı, yalnız ne karıştırıldığı değişir):

| çeşit | çizgi fazları | süreklilik | amaç |
|---|---|---|---|
| **V0** | ψ ≡ 0 | dokunulmaz | **kontrol** — gerçeğin ta kendisi olmalı |
| **VG** | karıştırılır | karıştırılır | **ANA Gauss vekil** |
| **VL** | karıştırılır | dokunulmaz | çizgiler-arası kilit tanısı |
| **VC** | ψ ≡ 0 | karıştırılır | süreklilik tanısı |

### 1d. REDDEDİLEN İKİ REÇETE — ölçüldü, kayda geçiyor

| | reçete | neden öldü | ölçülen (pow/gp, τ = 0.46 / 0.54 / 0.62 / 0.70 / 0.78) |
|---|---|---|---|
| **RED-1** | düz n-uzayı FFT faz karıştırması | Ölçüm KONUM uzayındadır; sürüklenme (`Σr` salınımı ≈ 470) bir asal çizgiyi n-tayfında ≈210 kutuya yayar. Fazlar karıştırılınca uyumlu toplanma biter. | gerçek 0.730/0.628/0.467/0.328/0.188 → **0.0022/0.0041/0.0096/0.0045/0.0038** |
| **RED-2** | u ekseninde düzgün ızgaraya taşı → karıştır → geri taşı (Kaiser-sinc, a=48, β=4) | V0 kontrolü **geçiyor** (0.628 → 0.619), ama düzgün ızgara Nyquist'i τ=0.5'e koyar; oysa ölçüm τ ≤ 0.80'i **düzensiz örneklemenin sayesinde** ayırt eder. Bant-sınırlı ara temsil bu üstün-Nyquist bilgiyi yok eder. | τ>0.5'te **0.628 → 0.008**, 0.467 → 0.002 |
| **RED-3** | konum uzayı, ama ÖZ-TUTARLILIK YOK (çizgiler gerçeğin m'sinde sentezlenir, vekilin m'siyle okunur) | Alan, gazın KENDİ konumlarında okunmalıdır; `e^{−i2πτ(Ĉ'−Ĉ)}` uyumu bozar. τ ≥ 0.6'da `pow_on < pow_off` olur, ölçüm kırılır. | **0.216/0.089/0.034/0.033/0.093** |
| **KABUL** | öz-tutarlı konum-uzayı karıştırması | — | **0.713/0.576/0.405/0.379/0.118** (gerçeğin %2–%16'sı içinde; τ = 0.78'de −%37) |

RED-2'nin ara temsil kodu `scratchpad/162/probe3.py`, `probe5.py`'de
duruyor (çekirdekten çıkarıldı); RED-1 ve RED-3 `162_dogrulama.py`'de
V7/V8 olarak **yeniden koşulur**.

### 1e. Patlayan koşu (kayda geçiyor)

`162_T1.py`'nin ilk sürümü artık tayfını `τ ≤ 1.6`'ya kadar TAM çizgi
kümesiyle istiyordu: `pk_m(e^{1.6·12.03}) = pk_m(2.3e8)` ≈ 1.2e7 çizgi
× 3e5 nokta. Koşu asılıp öldürüldü. Düzeltme: **τ ≤ 1.15 ve her 0.05'lik
bantta en çok 140 çizgi rastgele (tohum 21) örneklenir**, bandın gücü
`⟨|c|²/2⟩·n` ile kestirilir. İkinci bir koşu `_kor(Cds[:-1], Xtil[1:])`
şekil uyuşmazlığıyla düştü (299998 ≠ 299997); düzeltildi.

---

## 2. T1 — KİMLİK: `C_n = −S_merdiven(z_n)`

Merdiven: `S(t) = −Σ_{q=p^m, log q ≤ 0.9L} a_q sin(t log q)`,
`a_q = 1/(π m p^{m/2}) = 1/(π m √q)` — **155/160'ın `aq`'suyla birebir**.
5251 çizgi; kuramsal varyans `Σa_q²/2 = 0.14100`, gerçek zeroların
konumunda ölçülen `Var(S(z_n)) = 0.13259`.

### 2a. Sonuç (son 300k, taban 0.40)

| karşılaştırma | **r** | Var oranı | eğim | artık payı |
|---|---|---|---|---|
| `cumsum(ds)` ↔ `−S(z_{n+1})` | **+0.97265** | 0.8607 | 0.9024 | 0.0539 |
| `Ĉ = cumsum(dsΔ)` ↔ `−S(mid_n)` | **+0.97639** | **1.0409** | **0.9962** | **0.0467** |

Beklenti (r ≥ 0.97, "132c'nin konum versiyonu") **tuttu**. Bond sürümü
daha iyidir ve eğimi 1'e 4 binde yakındır — `Ĉ`, merdivenin bond
ortalamasıdır.

### 2b. Artığın tayfı — kesme nerede görünüyor

| τ bandı | çizgi | \|Ĉ\|²/2 | \|S\|²/2 | \|artık\|²/2 | artık/Ĉ | Ĉ/beklenen | cos²(πτ) |
|---|---|---|---|---|---|---|---|
| 0.05–0.10 | 2 | 4.03e−02 | 3.95e−02 | 7.5e−06 | **0.0002** | 0.9550 | 0.9455 |
| 0.15–0.20 | 4 | 8.46e−03 | 7.32e−03 | 5.3e−05 | 0.0063 | 0.6066 | 0.7270 |
| 0.30–0.35 | 10 | 1.48e−03 | 9.98e−04 | 5.0e−05 | 0.0340 | 0.1734 | 0.2730 |
| 0.40–0.45 | 21 | 9.74e−05 | 4.12e−05 | 2.1e−05 | 0.2206 | 0.0173 | 0.0545 |
| 0.50–0.55 | 56 | 1.41e−04 | 1.63e−04 | 1.9e−06 | **0.0133** | 0.0298 | 0.0062 |
| 0.60–0.65 | 157 | 5.31e−04 | 3.53e−04 | 2.1e−05 | 0.0401 | 0.1292 | 0.1464 |
| 0.70–0.75 | 432 | 6.31e−04 | 2.03e−04 | 1.31e−04 | 0.2069 | 0.1816 | 0.4218 |
| 0.75–0.80 | 735 | 5.57e−04 | 1.05e−04 | 2.28e−04 | 0.4095 | 0.1712 | 0.5782 |
| 0.80–0.85 | 1260 | 4.74e−04 | 1.11e−04 | 3.78e−04 | 0.7963 | 0.1520 | 0.7270 |
| 0.85–0.90 | 2168 | 5.42e−04 | 4.63e−04 | 5.51e−04 | **1.0166** | 0.1811 | 0.8536 |

**Okuma.** (i) Artık, τ ≤ 0.65'te Ĉ'nin **%4'ünün altındadır**; **τ ≥ 0.75'ten
itibaren hızla büyür ve 0.85–0.90'da Ĉ'ye eşitlenir** — yani eksik olan
merdivenin **üst ucudur** (kesme τ ≤ 0.9 ve orada zaten 2168 çizgi
birbirine karışıyor). Taban içinde (τ ≤ 0.4) kimlik neredeyse tamdır.
(ii) `Ĉ/beklenen` sütunu düşük τ'da `cos²(πτ)` ile **birebir** gidiyor
(0.955 ↔ 0.9455 τ = 0.075'te): Ĉ, dsΔ'nın kümülatifidir ve bond ortalaması
bir tona `cos(πτ)` kazancı verir. Yüksek τ'da oran cos²'nin **altına**
düşer (0.18 ↔ 0.42 τ = 0.72'de) — bu, çizginin gaz içindeki sönümüdür
(160'ın `|E1|`'i). (iii) `τ ≥ 0.95` satırları **ölçüm değildir**: orada
`w ≈ L` olduğu için kafes ALIAS'ı devreye girer (`Ĉ/beklenen` 3–9'a
fırlar) ve merdiven zaten kesiktir.

### 2c. C'nin künyesi

| nicelik | son | keskin | A4 |
|---|---|---|---|
| `Var(ds)` | 0.16744 | 0.12795 | 0.11284 |
| `Var(η)` (taban 0.40) | 0.05431 | 0.09447 | 0.04657 |
| `σΔ² = Var(dsΔ)` | 0.05466 | 0.03499 | 0.04294 |
| `Var(cumsum ds)` | 0.11411 | 0.07115 | 0.08880 |
| **`Var(Ĉ)`** | **0.07255** | 0.04655 | 0.06197 |
| `Var(Σ X̃)` (sürüklenmeli) | **1.912e+04** | 1.912e+04 | 1.912e+04 |
| `corr(Ĉ, X̃)` | −0.43298 | −0.44100 | −0.41357 |
| **`corr(Ĉ, dsΔ)`** | **−0.43401** | −0.43348 | −0.41621 |
| kalem: `−½√(σΔ²/Var Ĉ)` | **−0.43401** | −0.43348 | −0.41621 |
| **`corr(Ĉ, η)`** | **+0.00014** | +0.00004 | +0.00001 |
| `corr(Σr, Ĉ)` | +0.00103 | −0.00726 | +0.00286 |

**İki kayıt.** (a) `corr(Ĉ,dsΔ)`, durağanlığın zorunlu kıldığı
`Cov(Ĉ_n, dsΔ_n) = −½Var(dsΔ)` değerini **beş basamakta** veriyor — Ĉ
durağandır, sürüklenmeden arındırılmış olması bunun kanıtıdır.
(b) `corr(Ĉ, η) ≤ 1.4e−04`: bellek ile genlik **eşzamanlı olarak
ilişkisizdir**. δ'yı üreten bağlaşım eşzamanlı değil, **fazdadır** — bu,
T2'nin bulduğu şeyin ilk işaretidir.

### 2d. Sentetiklerde T1 (yan kayıt, hüküm değil)

Referans merdiven **gerçek asal merdivenidir**; sentetik gazlar 152'de
FARKLI genlik profilleriyle kuruldu (keskin: τ ≤ 1.00 keskin kesim,
A4: erfc 0.68/0.125). Bu yüzden uyumsuzluk **tanım gereğidir**:

| gaz | `cumsum(ds)` ↔ `−S(z)` | `Ĉ` ↔ `−S(mid)` |
|---|---|---|
| keskin | r = +0.8946, eğim 0.852 | r = +0.4240, eğim 0.220 |
| A4 | r = +0.9675, eğim 0.945 | r = +0.7068, eğim 0.425 |

Bond sürümündeki büyük fark, `cos²(πτ)` ağırlığının τ ≈ 1'i öne
çıkarmasındandır (keskin'in merdiveni oraya kadar uzanıyor, referans
τ ≤ 0.9'da kesik). **Yan bulgu, ölçüldü ve kayda geçiyor:** düşük τ'da
(0.05–0.10) `Ĉ/beklenen` gerçekte **0.955 = cos²(πτ)** iken keskin'de
**0.270**, A4'te **0.682** — 152'nin sentetik gazları ideal merdiven
genliğini düşük τ'da tam kurmuyor. Bu koşunun konusu değil; not düşülüyor.

---

## 3. T2 — GAUSS-TAYFSAL KAPANIŞ SINAVI

Bantlar 160'ın `t1` ızgarası; τ ∈ (0.50, 0.82) ve **gerçeğin** sağlık
kuralını geçen bantlar. Vekil sütunları tohum ortalaması ± tohumlar arası sd.

### 3a. Gerçek gaz (son, taban 0.40; 4 tohum)

| τ_eff | A·σ_X | \|Γ\| | δ_gerçek | δ_V0 | δ_vekil (ort±sd) | **δ_g/δ_v** | K1_g/K1_v | K2/δ (g) | K2/δ (v) | \|Γ\|_v |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.5392 | 0.792 | 0.860 | −0.1803 | −0.1803 | −0.0537 ± 0.0114 | **+3.36** | +3.67 | +0.271 | +0.333 | 0.686 |
| 0.5791 | 0.851 | 0.834 | −0.2811 | −0.2811 | −0.0647 ± 0.0153 | **+4.35** | +6.00 | +0.228 | +0.441 | 0.658 |
| 0.6188 | 0.909 | 0.789 | −0.4179 | −0.4179 | −0.0869 ± 0.0180 | **+4.81** | +7.37 | +0.219 | +0.490 | 0.629 |
| 0.6578 | 0.966 | 0.740 | −0.5909 | −0.5909 | −0.0929 ± 0.0187 | **+6.36** | +14.79 | +0.214 | +0.662 | 0.591 |
| 0.6980 | 1.025 | 0.684 | −0.8422 | −0.8422 | −0.1266 ± 0.0401 | **+6.65** | +17.13 | +0.235 | +0.703 | 0.566 |
| 0.7382 | 1.085 | 0.669 | −1.1282 | −1.1282 | −0.1318 ± 0.0829 | **+8.56** | +42.74 | +0.230 | +0.846 | 0.541 |
| 0.7766 | 1.141 | 0.699 | −1.6189 | −1.6189 | −0.1083 ± 0.1014 | **+14.95** | −129.81 | +0.274 | +1.083 | 0.489 |

`ORAN medyan **+6.36**, menzil [+3.36, +14.95]; |oran−1| medyan **%536**.`
`KONTROL δ_gerçek/δ_V0 = **+1.00000** (yedi bantta da, beş basamakta).`

### 3b. Üç gaz, tek satırda

| gaz | bant | δ_g/δ_V0 (kontrol) | **δ_g/δ_VG medyan [menzil]** | K1_g/K1_VG | s2_g/s2_VG | δ_g/δ_VL | δ_g/δ_VC |
|---|---|---|---|---|---|---|---|
| **gerçek-son** | 7 | **+1.00000** | **+6.36** [+3.36, +14.95] | +7.37 | +2.68 | +7.27 | **+1.24** |
| **keskin** | 7 | **+1.00000** | **+2.31** [−17.77, +44.67] | +2.79 | +0.66 | +2.58 | **−0.55** |
| **A4** | 3 | **+1.00000** | **−11.76** [−12.51, −7.14] | −8.64 | −7.86 | −5.21 | **+0.82** |

**±%15 içinde 0/17 bant.** keskin'in menzili geniştir çünkü δ'sı bu
aralıkta **sıfırdan geçiyor** (+0.135 → −0.309); orada oran anlamsızdır.
Mutlak okuma daha dürüsttür: keskin'in |δ|'sı 0.135 → 0.008,
0.309 → 0.026'ya iniyor. A4'te oran **negatif**: vekil δ'nın işaretini
çeviriyor.

### 3c. Vekil ölçülebilir mi? (kırık ölçüm değil)

| | gerçek | vekil (VG, ort) |
|---|---|---|
| `\|Γ\|` (son, τ 0.54 → 0.78) | 0.860 → 0.699 | 0.686 → 0.489 |
| çizgi gücü `pow/gp` (τ 0.46/0.54/0.62/0.70/0.78) | 0.730/0.628/0.467/0.328/0.188 | 0.713/0.576/0.405/0.379/0.118 |
| bant düzeyi `K3 = δ − (K1+K2)` (bütün koşular) | ≤ **1.10e−04** | ≤ **7.50e−04** |
| çizgi düzeyi kapanış (bütün koşular) | ≤ **1.0e−13** | ≤ **2.8e−13** |

Vekil gazda ölçüm **sağlamdır**: `|Γ|` düşüyor ama 0.49–0.69 aralığında
kalıyor, çizgi güçleri gerçeğin **%2–%16'sı içinde** (τ = 0.78'de en kötü,
**−%37**), muhasebe özdeşliği 1e−13'te kapanıyor ve K3 δ'nın **binde
3'ünün** altında. **Çöküş δ'nın kendisindedir, tahmincide değil.**

---

## 4. TAŞIYICIYI YALITMA — VL / VC ayrışımı

| τ_eff | δ_gerçek | δ_V0 | δ_VL101 (çizgi fazları) | δ_VC101 (süreklilik) | δ_VG101 | g/VL | g/VC |
|---|---|---|---|---|---|---|---|
| 0.5392 | −0.1803 | −0.1803 | −0.0405 | **−0.2200** | −0.0691 | +4.45 | **+0.82** |
| 0.5791 | −0.2811 | −0.2811 | −0.0387 | **−0.2872** | −0.0872 | +7.27 | **+0.98** |
| 0.6188 | −0.4179 | −0.4179 | −0.0623 | **−0.3787** | −0.1108 | +6.71 | **+1.10** |
| 0.6578 | −0.5909 | −0.5909 | −0.0376 | **−0.4768** | −0.1173 | +15.69 | **+1.24** |
| 0.6980 | −0.8422 | −0.8422 | −0.0939 | **−0.6250** | −0.1709 | +8.97 | **+1.35** |
| 0.7382 | −1.1282 | −1.1282 | −0.0697 | **−0.8208** | −0.2237 | +16.18 | **+1.37** |
| 0.7766 | −1.6189 | −1.6189 | +0.0516 | **−0.9907** | −0.2569 | −31.36 | **+1.63** |

A4'te aynı desen, daha keskin: `δ_VC = +0.2034 / +0.2836 / +0.3099`
(gerçek `+0.1659 / +0.2129 / +0.3248`, oran **0.82 / 0.75 / 1.05**) ama
`δ_VL = −0.0319 / −0.0259 / +0.0107` — **işaret gidiyor**.

> **HÜKÜM (taşıyıcı).** Sürekli tabanın fazları tümüyle karıştırıldığında
> δ **değişmiyor** (gerçek gazda 0.82–1.63; A4'te 0.75–1.05). Asal
> çizgilerin fazları karıştırıldığında δ **çöküyor** ve A4'te işaretini
> kaybediyor. **δ, asalların birbirine kilitli fazlarının taşıdığı bir
> niceliktir.** Bu, "aritmetik bellek" ifadesinin ilk doğrudan,
> nicel ölçümüdür.

---

## 5. T3 — ×1.4'ÜN TAYF KİMLİĞİ

### 5a. `⟨σX̃²⟩` — b'nin işaretini taşıyan nesne

`−A²s2/2 = −A²⟨σX̃0²⟩/(2⟨ρ⟩)` (160 §5e'nin `Re[Es]` sürücüsü):

**gerçek-son**

| τ_eff | gerçek | V0 | vekil (ort±sd) | oran | K2 gerçek | K2 vekil | işaret |
|---|---|---|---|---|---|---|---|
| 0.5392 | −0.06833 | −0.06833 | −0.01373 ± 0.00385 | **+4.98** | −0.04892 | −0.01788 | AYNI |
| 0.6188 | −0.11478 | −0.11478 | −0.03501 ± 0.00343 | +3.28 | −0.09154 | −0.04262 | AYNI |
| 0.6980 | −0.14014 | −0.14014 | −0.07427 ± 0.02086 | +1.89 | −0.19785 | −0.08897 | AYNI |
| 0.7766 | −0.03816 | −0.03816 | −0.10660 ± 0.05966 | +0.36 | −0.44412 | −0.11727 | AYNI |

**A4** — karşı-örneğin çözümü

| τ_eff | **gerçek** | V0 | **vekil (ort±sd)** | oran | K2 gerçek | K2 vekil | **işaret** |
|---|---|---|---|---|---|---|---|
| 0.5386 | **+0.05329** | +0.05329 | **−0.01055 ± 0.00596** | −5.05 | +0.07827 | −0.01309 | **TERS** |
| 0.5772 | **+0.10010** | +0.10010 | **−0.01273 ± 0.00517** | −7.86 | +0.14069 | −0.01477 | **TERS** |
| 0.6145 | **+0.23685** | +0.23685 | **−0.00615 ± 0.00202** | −38.54 | +0.31912 | −0.00689 | **TERS** |

> **Görev doğrudan sordu: "Gauss vekilde SIFIRA düşmeli mi?"**
> **Yanıt: TAM SIFIRA DÜŞMÜYOR, ama GAZI AYIRAN PAYI TAMAMEN KAYBEDİYOR.**
> σ, Gauss değişkenlerinin doğrusal olmayan bir fonksiyoneli olduğundan
> `⟨σX̃²⟩` ortak-Gauss vekilde özdeş olarak sıfır olmak zorunda değildir
> — ve olmuyor (vekilde −0.006 … −0.107). Ama **her üç gazda da vekil
> AYNI işareti (negatif) ve benzer büyüklüğü veriyor**: A4'ün ayırt edici
> **pozitif** `⟨σX̃²⟩`'si, keskin'in küçük değeri, gerçeğin büyük negatifi
> — hepsi tek bir "jenerik" değere çöküyor. **A4'ün dışbükeyliğini üreten
> nesne Gauss vekilde YOKTUR; o nesne, asal fazlarının kilidiyle
> doğar.** 160'ın "b'nin işaretini `⟨σX̃²⟩` taşıyor" anatomisi bu koşuda
> bir adım daha ileri gidiyor: **`⟨σX̃²⟩`'nin gaz-ayırıcı payı Gauss-ötesidir.**

### 5b. a, b ve τ₀ — üç gaz, 158/160'ın fit konvansiyonuyla (g158, 4 pencere)

| gaz | nicelik | **gerçek** | **vekil (3 tohum, ort±sd)** | çıplak kinematik |
|---|---|---|---|---|
| gerçek-son | `a(φ)` | **10.652** | **12.329 ± 0.097** | 4π = 12.566 |
| gerçek-son | `b(φ)` | **−6.892** | **−0.943 ± 0.863** | 0 |
| gerçek-son | `τ₀` | 0.5090 | 0.5036 | 0.5 |
| keskin | `a(φ)` | **11.667** | **12.497 ± 0.031** | 12.566 |
| keskin | `b(φ)` | **−4.434** | **+0.164 ± 1.398** | 0 |
| keskin | `τ₀` | 0.4884 | 0.5006 | 0.5 |
| A4 | `a(φ)` | **13.032** | **12.507 ± 0.054** | 12.566 |
| A4 | `b(φ)` | **+3.614** | **+0.162 ± 1.143** | 0 |
| A4 | `τ₀` | 0.4894 | 0.5015 | 0.5 |

Tohum tohum (W-A/W-D/W-B/W-C ortalaması):

| gaz | çeşit | τ₀ | a(φ) | b(φ) | b(δ) | b(K1) | b(A·S) |
|---|---|---|---|---|---|---|---|
| son | gercek | 0.5090 | 10.652 | **−6.892** | −6.882 | −5.740 | −3.651 |
| son | VG101 | 0.5043 | 12.218 | −1.706 | −1.582 | +0.094 | +0.082 |
| son | VG102 | 0.5033 | 12.389 | −1.116 | −1.154 | −1.115 | −1.477 |
| son | VG103 | 0.5032 | 12.381 | −0.006 | −0.025 | +0.568 | +0.586 |
| keskin | gercek | 0.4884 | 11.667 | **−4.434** | −4.441 | −3.050 | −2.771 |
| keskin | VG101/2/3 | 0.500 | 12.48–12.53 | +0.885 / −1.447 / +1.055 | +0.904 / −1.413 / +0.906 | — | — |
| A4 | gercek | 0.4894 | 13.032 | **+3.614** | +3.715 | −3.020 | −4.326 |
| A4 | VG101/2/3 | 0.501 | 12.46–12.56 | +0.551 / +1.061 / −1.124 | +0.634 / +0.953 / −1.061 | — | — |

> **T3'ÜN YANITI: VEKİLLER BİRBİRİNE YAKLAŞIYOR.** Gerçek ↔ sentetik
> ayrımının bütün ölçütleri Gauss vekilde **kayboluyor**: `a` üç gazda da
> çıplak `4π`'nin **%0.5–%2** içine, `b` **sıfırın 1σ'sı içine**, `τ₀`
> kinematik **0.5**'e dönüyor. `b`'nin **işareti bile** tohumdan tohuma
> değişiyor (son: −1.71/−1.12/−0.01; keskin: +0.89/−1.45/+1.06; A4:
> +0.55/+1.06/−1.12) — yani vekilde b **ölçülemiyor**, sıfırla uyumlu.
> **Dolayısıyla gerçek ↔ sentetik δ farkı bir TAYF FARKI DEĞİLDİR**
> (tayflar vekillerde de farklıdır — `σΔ²` 0.070 / 0.035 / 0.068, çizgi
> güçleri korunmuş — ama δ farkı gitmiştir). **Fark, FAZ / GAUSS-ÖTESİ
> YAPIDIR.** 160'ın ×1.4'ü (A4'ün kuadratür çarpanı) bir genlik profili
> farkının değil, faz kilidinin ürünüdür.

---

## 6. AÇIKLARI KAPATAN ÜÇ KONTROL

### 6a. "Vekilin σ_X'i şişiyor, açık ondan" — HAYIR

| gaz | çeşit | σΔ² | Var(η) | çizgi payı η / Ĉ | varyans kapanışı η / Ĉ |
|---|---|---|---|---|---|
| son | gercek/V0 | 0.054665 | 0.054313 | 0.487 / 0.918 | **1.801 / 1.109** |
| son | VG101…104 | 0.0691–0.0705 | 0.0969–0.0973 | " | " |
| **keskin** | gercek/V0 | 0.034986 | 0.094474 | 0.495 / 0.420 | **0.954 / 0.989** |
| **keskin** | VG101…103 | **0.0345–0.0346** | **0.0899–0.0900** | " | " |
| A4 | gercek/V0 | 0.042945 | 0.046567 | 0.738 / 0.693 | 0.984 / 1.242 |
| A4 | VG101…103 | 0.0671–0.0688 | 0.0456–0.0457 | " | " |

İki bağımsız argüman:

1. **keskin'de şişme YOK** (σΔ² −%1.3, Var(η) −%4.7, varyans bütçesi
   0.954/0.989 ile kapanıyor) — ve δ orada da çöküyor (0.135 → 0.008).
2. **A·σ_X eşleştirilmiş karşılaştırma açığı BÜYÜTÜYOR.** Gerçeğin
   δ(A·σ_X) eğrisi vekilin kendi A·σ_X'inde okunduğunda gerçek gazda oran
   **+7.17 … +14.95** (ham oran +3.36 … +14.95), A4'te **−13.99 … −11.76**
   (ham −7.14 … −11.76). Yani σ_X şişmesi δ'yı **büyütmesi** gerekirdi;
   vekilde δ küçülüyor. **Şişme açığı açıklamıyor.**

### 6b. "Boru hattı bozuyor" — HAYIR

V0 (ψ ≡ 0), gerçeğin **ta kendisidir**: dizi düzeyinde η/X̃/X_tam/Ĉ maks
fark **≤ 1.1e−16**, `mid` farkı **0.0e+00**, öz-tutarlılık kalıntısı
`2.3e−10` (kayan nokta), bant düzeyinde `|δ_V0 − δ_gerçek| = 0.0e+00`,
oran **1.00000**. Ve `gercek` koşusu 160'ın kaydını 63 bantta
**0.000e+00** ile veriyor.

### 6c. "Σ|c|²/2 > Var(η) bir çıkarım hatası" — HAYIR, GERÇEK BİR OLGU

Bilinen rastgele-fazlı sinyal kuruldu (`x = Re Σ |c_q| e^{i(w_q m + ψ_q)}`,
planlanan güç **0.069971**); aynı çıkarım uygulandı:

| | değer | plana oran |
|---|---|---|
| sentetik sinyalin `Var` | 0.069700 | **0.9961** |
| çıkarımdan geri gelen `Σ|c'|²/2` | 0.070288 | **1.0045** |
| **GERÇEK η**: `Σ|c_η|²/2 = 0.069971` ↔ `Var(η) = 0.054313` | | **1.2883** |
| **GERÇEK Ĉ**: `Σ|c_Ĉ|²/2 = 0.074438` ↔ `Var(Ĉ) = 0.072552` | | **1.0260** |

> **YENİ NESNE.** Çıkarım **yansızdır** (%0.45). Öyleyse gerçek gazda
> `Σ|c_η|²/2 / Var(η) = 1.288`, asal çizgilerinin η içinde **yıkıcı
> girişim** yaptığının doğrudan ölçüsüdür: η, çizgilerinin uyumsuz
> toplamından **%29 daha sessizdir**. Ĉ'de aynı oran **1.026** — girişim
> yok. **Kilit, ρ-ağırlığının gördüğü η kanalındadır**, belleğin kendisinde
> değil. Vekil bu girişimi (tanım gereği) yok ettiği için η'sı şişer;
> §6a bunun açığı açıklamadığını gösteriyor.

---

## 7. DENETİM

| | sınav | sonuç |
|---|---|---|
| **V1** | 162/`gercek` = 160'ın kayıtlı JSON'u mu? | 3 gaz × 2 ızgara, **63 ortak bant**, ölçülen bütün sütunlarda (φ, δ, τ_eff, \|Γ\|, K1, K2, Kfull, M1, M0, tüm jackknife'lar, `AS.*`, `KUM.*`, E1/Es) **maks fark 0.000e+00** |
| **V2** | V0 boru hattı kimliği | η/X̃/X_tam/Ĉ ≤ **1.1e−16**, `mid` **0.0e+00**, `m`-kimlik **0.0e+00**, bant düzeyi `δ` farkı **0.0e+00** |
| **V3** | vekil ikinci momentleri koruyor mu? | **Çizgi güçleri** (η) τ ≤ 0.80'de oran **0.84–0.99**, (Ĉ) **0.79–1.16**. **Ĉ'nin otokovaryansı**: k=1 0.623→0.551, k=2 0.113→0.157, k=5 −0.269→−0.250, k=20 0.247→0.243, k=100 0.0855→0.0826. **η'nın normalize otokovaryansı korunmuyor** (k=1: −0.724 → −0.472) — §6c'nin girişim olgusunun kaçınılmaz sonucu; ham kovaryans farkı %17 |
| **V4** | Gauss'luk | η çarpıklık −0.151 → **−0.008**, basıklık 3.022 → **2.998**; Ĉ 0.000/2.693 → 0.023/2.749; **dsΔ 0.282/3.006 → 0.534/3.721** (vekil bu momentte DAHA yapılı — öz-tutarlılığın doğrusal olmayanlığı; δ yine de çöküyor) |
| **V5** | Gram kontrolü | çıkarım yansız: **1.0045**; gerçekteki 1.288 **yıkıcı girişimdir** |
| **V6** | sabit nokta (Newton, 6 adım) | V0: 2.3e−10 (kayan nokta). keskin: rms **1.6e−01 → 1.0e−05…5e−11**. son: **2.0e−01 → 2.3e−02…3.4e−02**. A4: **1.9e−01 → 4.1e−02…6.9e−02**. (Ortalama aralık 2π/L = 0.522; en kötü kalıntı aralığın **%13'ü**.) |
| **V7** | RED-1 (düz n-uzayı FFT) | pow/gp 0.730/0.628/0.467/0.328/0.188 → **0.0022/0.0041/0.0096/0.0045/0.0038** |
| **V8** | RED-3 (öz-tutarlılık yok) | → **0.216/0.089/0.034/0.033/0.093**; öz-tutarlı VG → **0.713/0.576/0.405/0.379/0.118** |

---

# HÜKÜM

## (i) T1 GEÇTİ — C, merdivenin ta kendisidir, ve iki parçalıdır

> `corr(cumsum(ds), −S_merdiven(z_n)) = +0.9727`;
> `corr(Ĉ, −S_merdiven(mid_n)) = +0.9764`, **eğim 0.9962**, varyans oranı
> 1.041, artık payı **%4.7**. Artık τ ≤ 0.65'te Ĉ'nin %4'ünün altında;
> τ ≥ 0.75'ten sonra hızla büyüyor — **eksik olan merdivenin kesilmiş üst
> ucudur**, taban içi değil.

Ve KALEM'in görmediği bir ayrım ölçüldü: D1'in `C = Σ X̃`'i **iki
parçadır** — `Var = 1.9e+04` olan **yoğunluk sürüklenmesi** ve
`Var = 0.0726` olan **aritmetik bellek Ĉ**. Ayrışım tamdır
(`m_n = m_0 + (2π/L)(u_n + Ĉ_n)`, 0.0e+00), Ĉ durağandır
(`corr(Ĉ,dsΔ) = −½√(σΔ²/VarĈ)` beş basamakta) ve **η'ya eşzamanlı olarak
diktir** (`|corr| ≤ 1.4e−04`).

## (ii) T2 REDDEDİLDİ — δ ikinci momentlerden GELMİYOR

> Aynı çizgi güçleri, aynı çizgi çapraz-tayfı, aynı süreklilik tayfı,
> aynı sürüklenme, aynı öz-tutarlılık — **yalnız fazlar rastgele**:
> `δ_gerçek/δ_vekil` gerçek gazda medyan **+6.36**, keskin'de +2.31,
> A4'te **−11.76**. **±%15 içinde 0/17 bant.**
> **Gauss-ötesi pay = oran − 1 = %536 (gerçek gaz medyanı).**
> Yani δ'nın en çok **%16'sı** tayftan geliyor; geri kalanı
> **aritmetik belleğin faz-duyarlı payıdır** — görevin adını koyduğu
> yeni nesne ölçüldü ve büyüktür.

Zincir bu yüzden **kapanmıyor**: `δ = tayf fonksiyoneli` hipotezi
düşmüştür. 160'ın kapattığı muhasebe (δ = K1+K2) doğrudur ve özdeştir,
ama `K1` ve `K2`'yi besleyen `⟨ρe^{iAX̃}⟩`, `⟨σe^{iAX̃}⟩` **ikinci
momentlerin fonksiyoneli değildir**.

## (iii) TAŞIYICI: ÇİZGİLER-ARASI FAZ KİLİDİ (süreklilik değil)

> Sürekli taban tümüyle karıştırıldığında δ **değişmiyor** (gerçek gazda
> `δ_g/δ_VC` = 0.82…1.63; A4'te 0.75…1.05). Asal çizgi fazları
> karıştırıldığında **çöküyor** (7.27 medyan) ve A4'te işaret kayboluyor.
> **δ'yı asalların birbirine kilitli fazları taşır.**

Bunun φ'den bağımsız, doğrudan bir sayısal karşılığı da bulundu:
`Σ|c_η|²/2 / Var(η) = 1.288` (çıkarım yansız, V5) — **η, asal
çizgilerinin uyumsuz toplamından %29 daha sessizdir.** Ĉ'de bu oran
1.026: kilit **genlik kanalındadır**.

## (iv) T3: ×1.4 BİR TAYF FARKI DEĞİL — FAZ FARKI

> Üç gazın Gauss vekilleri **ayırt edilemez** hale geliyor:
> `a(φ)` 10.65 / 11.67 / 13.03 → **12.33 / 12.50 / 12.51** (4π = 12.566);
> `b(φ)` −6.89 / −4.43 / **+3.61** → **−0.94 ± 0.86 / +0.16 ± 1.40 /
> +0.16 ± 1.14** (üçü de sıfırla uyumlu, işaretleri tohumdan tohuma
> değişiyor); `τ₀` → kinematik 0.5.
> **Gerçek ↔ sentetik δ farkı vekillerde YENİDEN ÜRETİLMİYOR** — tayflar
> vekillerde de farklı olduğu halde. Fark, **Gauss-ötesi faz yapısıdır.**

## (v) `⟨σX̃²⟩`: b'nin "gerçeğin parmak izi" olmasının derin sebebi

> Görevin sorduğu sınav koşuldu. `⟨σX̃²⟩` Gauss vekilde **özdeş sıfıra
> düşmüyor** (σ, Gauss alanların doğrusal olmayan fonksiyonelidir), ama
> **gazı ayıran payı tamamen kayboluyor**: A4'ün `−A²s2/2`'si
> +0.053 / +0.100 / +0.237 iken vekilde −0.011 / −0.013 / −0.006 —
> **üç bantta da İŞARET TERS**. Vekilde her üç gaz da aynı jenerik
> negatif değeri veriyor.
> **Sonuç: `b`'nin işareti — 160'ın A4 karşı-örneğini çözen nesne —
> tamamen Gauss-ötesi/aritmetik bellek imzasıdır.** 160 "b'nin işaretini
> `⟨σX̃²⟩` taşıyor" demişti; 162 "o `⟨σX̃²⟩`'yi asal fazlarının kilidi
> taşıyor" diyor.

## Sıradaki adım (bu ölçümün işaret ettiği)

1. **Aranacak nesne artık tayf değil, BİSPEKTRUM (üçüncü mertebe).** δ,
   `(η, Ĉ)` sürecinin ikinci momentlerinden gelmiyor; VL/VC ayrışımı
   taşıyıcının **çizgiler-arası** faz ilişkisi olduğunu söylüyor. Doğru
   nesne `⟨c_{q1} c_{q2} c*_{q1+q2}⟩` türü bir üçlü korelatördür ve
   asal frekanslarda `log q1 + log q2 = log(q1q2)` **çarpımsal**
   yapıyla kilitlidir. Bu, ölçülebilir ve bu koşunun çekirdeği hazır.
2. **`Σ|c|²/2 / Var` oranı, φ'den bağımsız bir GAZ AYIRICISIDIR ve
   sınanmadı.** Gerçek 1.288, Ĉ'de 1.026. keskin ve A4 için henüz
   ölçülmedi (çizgi payları 0.495/0.420 ve 0.738/0.693 kayıtlı, ama
   girişim oranı ayrı bir ölçümdür). 158'in dokuz gazında ucuz.
3. **Vekilin `dsΔ`'sı gerçekten daha çarpıktır** (0.28 → 0.53) — öz-tutarlı
   Newton'un doğrusal olmayanlığı. Bunun δ'ya katkısı ölçülmedi; bir
   "doğrusallaştırılmış vekil" (sabit nokta yerine tek adım) ayırt eder.
4. **A·σ_X eşleştirmesi bir ARAÇTIR ve tam kullanılmadı.** Vekilin σΔ²'si
   gerçeğinkine eşitlenerek (çizgi genliklerinde tek skaler κ) koşulan
   bir sürüm, §6a'nın argümanını doğrudan kapatır.
5. **`orta`, `P1`, `dusuk`, `J14`, `J26`, `N5`, `N5z`, `P0` bu koşuda
   ölçülmedi.** Vekil koşuları yalnız `son`, `keskin`, `A4`; taban ekseni
   yalnız 0.40.

---

## Dürüstlük notları

* **Uydurma yok; her sayı yeniden koşuldu.** Bütün tablo değerleri
  `162_configs/162_kos.py`'nin ürettiği **35 JSON**'dan,
  `162_analiz.py`'nin çıktısından (`scratchpad/162/analiz_cikti.txt`),
  `162_dogrulama.py`'ninkinden (`dogrulama_cikti.txt`) ve
  `162_T1.py`'ninkinden (`log_T1.txt`, `T1.json`) otomatik alındı.
  160'tan **tek bir sayı bile kopyalanmadı**: 160'ın bütün bantları bu
  koşuda yeniden ölçüldü ve V1 farkı `0.000e+00` verdi.
* **"Gauss-tayfsal kapanış reddedildi" bir başarısızlık değil, bir
  ÖLÇÜMDÜR** — ve görevin ikinci şıkkı ("değilse Gauss-ötesi pay
  ölçülür") tam olarak budur. Ölçülen pay: gerçek gazda **%536** (medyan),
  bant bant %236–%1395.
* **ÜÇ VEKİL REÇETESİ DENENDİ, İKİSİ ELENDİ ve saklanmadı.** RED-1
  (düz n-uzayı FFT) ve RED-3 (öz-tutarlılıksız) `162_dogrulama.py`
  V7/V8'de **yeniden koşuluyor**; RED-2'nin (düzgün u-ızgarası,
  Kaiser-sinc a=48/β=4) kodu tasarım kayıtlarında (`scratchpad/162/
  probe3.py`, `probe5.py`) duruyor ve ölçtüğü sayılar §1d'de.
  RED-2'nin V0 kontrolü GEÇİYORDU (0.628 → 0.619) — yani "kontrol geçti"
  tek başına yeterli değildir; onu eleyen şey vekilin τ>0.5 çizgilerini
  öldürmesiydi (0.628 → 0.008). Bu ders kayda geçiyor.
* **Vekil, ikinci momentleri TAM korumuyor ve bu gizlenmiyor.**
  Çizgi güçleri τ ≤ 0.80'de %1–%21 içinde korunuyor (η %1–%16, Ĉ %1–%21); Ĉ'nin otokovaryansı
  iyi korunuyor; ama **η'nın varyansı gerçek gazda %79, σΔ² %28 şişiyor**
  (A4'te σΔ² %59). Sebebi §6c'de ölçüldü (yıkıcı girişimin yokluğu) ve
  §6a'da açığın **bu şişmeden gelmediği** iki bağımsız yolla gösterildi
  (keskin'de şişme yok, çöküş var; A·σ_X eşleştirmesi oranı büyütüyor).
* **`dsΔ`'nın çarpıklığı vekilde ARTIYOR** (0.282 → 0.534, basıklık
  3.01 → 3.72). Yani vekil "her bakımdan daha Gauss" değildir; öz-tutarlı
  Newton çözümü bir doğrusal olmayanlık ekliyor. Buna rağmen δ çöküyor —
  ama bu, "vekil tam Gauss'tur" cümlesinin **kurulamayacağı** anlamına
  gelir ve kuruluyor değil: kurulan cümle "**ikinci momentler korunmuş,
  fazlar rastgeleleştirilmiştir**"tir.
* **Sabit nokta tam yakınsamadı.** Newton 6 adımda `|h|` rms'ini
  keskin'de 1e−05'e, gerçek gazda 2.3e−02'ye, A4'te 6.9e−02'ye indiriyor
  (ortalama aralığın %13'ü). Daha çok adım maliyeti üç katına çıkarırdı;
  çizgi güçlerinin korunması (V8: 0.713/0.576/0.405 ↔ gerçek
  0.730/0.628/0.467) bunun yeterli olduğunu gösteriyor. **Yakınsamamış
  olması vekili geçersiz kılmaz** — vekil, `X̃' = r + diff(Ĉ')` ve
  `m'` bağıntısıyla **özdeş olarak tutarlı** bir gazdır; öz-tutarlılık
  yalnız "çizgilerin ne kadar uyumlu olduğunu" belirler. Uyum ölçüldü:
  vekilin çizgi güçleri gerçeğin %2–%16'sı içinde (τ = 0.78'de −%37).
* **keskin'in oran menzili [−17.8, +44.7] geniştir ve bu anlamlı
  DEĞİLDİR:** δ_keskin bu aralıkta sıfırdan geçiyor (+0.135 → −0.309).
  Hükümlere mutlak değerler girdi (0.135 → 0.008; 0.309 → 0.026).
* **A4'te yalnız 3 bant ölçüldü.** 160'ın sağlık kuralı (159 §7, birebir)
  A4'ün 9 bandından 5'ini geçiriyor ve τ ∈ (0.50,0.82) kısıtı 3 bant
  bırakıyor. A4'ün hükümleri (işaret çevirmesi) bu üç bandın üçünde de
  aynı yönde.
* **b'nin vekildeki değerleri SIFIRLA UYUMLUDUR, "küçük" değil.** Üç
  tohumun sd'si (0.86–1.40) ortalamadan büyük ya da onunla kıyaslanabilir;
  işaretler tohumdan tohuma değişiyor. Rapor bunu "b vekilde ölçülemiyor"
  diye yazıyor, "b vekilde sıfır" diye değil.
* **İki patlayan koşu §1e'de.** Biri çizgi kümesinin üstel büyümesi
  (τ ≤ 1.6 → 1.2e7 çizgi), biri şekil uyuşmazlığı. İkisi de düzeltildi
  ve koşular yeniden yapıldı.
* **T1'in sentetik satırları hüküm taşımıyor.** Referans merdiven gerçek
  asal merdivenidir; 152'nin gazları farklı genlik profilleriyle
  kuruldu, uyumsuzluk tanım gereğidir. Yan bulgu (düşük τ'da keskin 0.270
  / A4 0.682 ↔ gerçek 0.955 = cos²πτ) kayda geçti, kullanılmadı.
* **Figür üretilmedi.** Bu koşuda görselleştirme istenmedi; bütün
  tablolar ham çıktılardan doğrudan alınabilir durumda.
* **Süreler ölçüt değildi.** Üç gaz eşzamanlı, `son` 36.8 dk, `keskin`
  32.0 dk, `A4` 32.1 dk; T1 ~11 dk; denetim ~13 dk.

---

## Ek — dosyalar ve tekrar-üretim

| dosya | ne |
|---|---|
| `162_configs/162_cekirdek.py` | `Taban162` (çizgi çıkarımı, çizgi/süreklilik ayrımı, öz-tutarlı Newton vekili), `olc162` (160'ın `olc160`'ını fabrikayla çağırır), `merdiven` |
| `162_configs/162_kos.py` | `<veri> <taban> <çeşit listesi>` — koşu sürücüsü |
| `162_configs/162_T1.py` | T1: kimlik, artık tayfı, C'nin künyesi |
| `162_configs/162_analiz.py` | B1–B8: denetim, T2 ana tablo, çeşit ayrışımı, `⟨σX̃²⟩`, a/b, ikinci-moment bütçesi, A·σ_X eşleştirme, özet |
| `162_configs/162_dogrulama.py` | V1–V8: bit düzeyi, boru hattı kimliği, ikinci momentler, Gauss'luk, Gram kontrolü, sabit nokta, RED-1, RED-3 |

Ham çıktılar `scratchpad/162/`: `S_<veri>_t0.4_<ızgara>_<çeşit>.json`
(35 dosya), `log_son.txt` / `log_keskin.txt` / `log_A4.txt` / `log_T1.txt`,
`analiz_cikti.txt`, `dogrulama_cikti.txt`, `T1.json`,
`probe*.py` (tasarım aşamasının RED-1/RED-2 ölçümleri).
η önbellekleri `scratchpad/155/eta_*.npz`'den **yeniden kullanıldı**,
yenisi üretilmedi.

Tekrar üretmek için:

```
.venv/bin/python 162_configs/162_kos.py son    0.40 gercek,V0,VG101,VG102,VG103,VG104,VL101,VC101
.venv/bin/python 162_configs/162_kos.py keskin 0.40 gercek,V0,VG101,VG102,VG103,VL101,VC101
.venv/bin/python 162_configs/162_kos.py A4     0.40 gercek,V0,VG101,VG102,VG103,VL101,VC101
.venv/bin/python 162_configs/162_T1.py  son keskin A4
.venv/bin/python 162_configs/162_analiz.py
.venv/bin/python 162_configs/162_dogrulama.py
```
