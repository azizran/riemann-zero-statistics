# 159 — mekanizma denklemi φ_Γ(τ) = A·S(τ)'nin sınavı: T1–T3

**Hedef.** KALEM (01 Eylül, `KALEM_AB_TURETIM_01EYL2026.md`) a ve b'nin
türetim zincirini tek bir parametresiz denkleme bağlamıştı:

> Γ_rot = ⟨w·e^{−iA·dsΔ}⟩/⟨w⟩,  w = yerel bant gücü ⇒
> **φ_Γ(τ) = A·S(τ)**,  S = Cov(P_yerel, dsΔ)/⟨P_yerel⟩,  A = 2πτ

"İki tarafı bağımsız ölçülebilir bir denklem — T1 sınavı." Bu rapor o
sınavı koştu: sol taraf 158'in ölçüm zinciriyle (τ_eff apsisi dahil,
**bit düzeyinde** yeniden üretilerek), sağ taraf çizgi bazında yerel güç
ile, dokuz–on bantta, iki tabanda, beş gazda.

**Figür:** `159_phi_AS_sinavi.png` (6 panel).
P1 ham φ_Γ ve kinematik omurga · P2 omurga çıkınca kalan faz δ ·
P3 mekanizma merdiveni (gerçek gaz) · P4 T1 hükmü, bant bant oran ·
P5 T2 kanal ayrışımı · P6 T3 a öngörüsü.

---

> ## Kısa hüküm
>
> 1. **T1 (birebir okunuşuyla) DÜŞTÜ, ve nedeni ölçüldü: φ_Γ, A·S DEĞİL.**
>    Gerçek gazda dokuz bandın **0'ı** ±%25 içinde; oran −38.7 ile +4.3
>    arasında zıplıyor. Ama başarısızlık gürültü değil, **yapısal ve tam
>    olarak adreslenebilir**.
> 2. **BULGU (bu koşunun en önemli çıktısı): ölçülen φ_Γ, ÖZDEŞ olarak
>    bir KİNEMATİK OMURGA + mekanizma fazıdır.**
>    ```
>    φ_Γ = (4π·τ_eff − 2π) + δ ,   δ ≈ A·S
>    ```
>    Bu bir yorum değil, cebir: `zp·conj(zc)·e^{−iA}/4 = ⟨(ρ+iσ)e^{+iA·X̃}⟩`
>    özdeşliği bağıl **7.8e−09** ile doğrulandı (V2), sentetik kontrolde
>    (bilinen tek çizgi + bilinen adım dizisi) e^{−iA} fazı **+0.015**'e
>    indiriyor, e^{+iA} ise **+1.272 = 2A−2π**'ye taşıyor (V3).
>    148'den beri bütün zincir (148/149/150/151/152/153/154/155/156/157)
>    `cr = zp*conj(zc)*np.exp(+1j*TWO_PI*W/L)` yazıyor; 148'in kendi
>    docstring'i niyeti "**bilinen faz ilerlemesi 2πτ ile döndür**" diye
>    yazmış. Ham korelatör o ilerlemeyi **zaten taşıdığı** için çarpan
>    ilerlemeyi **çıkarmıyor, ikiye katlıyor**.
> 3. **Bunun a ve τ₀ için sonucu ölçüldü ve tam:**
>    * **a = 4π + dδ/dτ** — beş gazda hata ≤ **0.10** (korel = **+0.9999**).
>      Yani 158'in "en dayanıklı ölçülmüş sayısı" a = 10.713 ± 0.059'un
>      **12.566'sı (= 4π) kinematiktir**, gaza duyarlı kısım yalnız
>      **dδ/dτ = −1.99 (gerçek) … +0.58 (A4)**.
>    * **τ₀ = ½ − δ(½)/a** — beş gazda fark ≤ **0.0011**. 155/157/158'in
>      τ₀ merdiveni (0.4886 … 0.5088) **omurganın ½'si + δ(½)/a**'dır.
> 4. **Kinematikten arındırılınca mekanizma denklemi ANLAMLI ama MÜHÜRSÜZ.**
>    Gerçek gazda δ/(A·S) **τ ≥ 0.54'te 1.82 – 2.47**, taban 0.40 ve 0.52'de
>    aynı (medyan 2.00 vs 2.02) — yani **kararlı bir ×2 açık**. Açığın iki
>    parçası ayrı ayrı ölçüldü: **birinci-mertebe kesme ×1.53** ve
>    **kuadratür (σ) kanalı ×1.30**. Çarpımları 2.0'ı veriyor.
>    **keskin** gazda aynı oran, sıfır geçişinden uzak altı bantta
>    **1.01 – 1.83** (τ_eff = 0.46/0.50/0.54 → 1.30/1.21/1.01;
>    0.66/0.70/0.74 → 1.83/1.29/1.02); 8 bandın **3'ü** ±%25 içinde,
>    δ/M1'de **4'ü** — mekanizma denklemi orada neredeyse tutuyor.
> 5. **Görevin verdiği P_loc tanımı (kayan pencere |ĉ|², W=64) ÖLÇÜLEBİLİR
>    DEĞİL, ve nedeni sayısal: GÜRÜLTÜ PEDESTALI.** ⟨|ĉ|²⟩ = pow/4 + σ_η²/W;
>    ölçülen pedestal/koherent oranı gerçek gazda **38 – 6371**, A4'te
>    101 – 10708, keskin'de 1.4e6'ya kadar. A·S(|ĉ|²) bu yüzden δ'nın
>    yalnız **%0.02 – %1.0**'i kadar. Denklemi ÖZDEŞ yapan ağırlık **ρ_n = Re[c_{n+1}·conj⟨c⟩]**
>    (çizginin yerel KOHERENT gücü, ⟨ρ⟩ = pow/4 tam olarak); bütün
>    birincil sayılar onunla ölçüldü. Sözel okuma (⟨η²⟩_W) frekans-kör
>    olduğu için hiçbir bandı ayırmıyor (A·S ≈ 2e−4, τ'dan bağımsız).
> 6. **Kayan pencere sınavının kendisi bir sonuç verdi: bağlaşım BOND
>    ölçeğinde.** X'i pencere üzerinden ortalamak A·S'yi silip süpürüyor:
>    |A·S(X̄_W)/A·S(X_bond)| = **%0.01 – %0.74 (W=32)**, **%0.13 – %1.18
>    (W=64)**, **%0.26 – %2.72 (W=128)**. Yani yerel güç ile adım
>    arasındaki korelasyonun 32–128 bond ölçeğinde bileşeni yok;
>    KALEM'in önerdiği "~64 bond" penceresi 2–3 mertebe yanlış ölçek.
> 7. **T2 tuttu ve 156'yı doğruladı — ama sentetiklerde TERSİ.** Kovaryans
>    doğrusal olduğu için S_tam = S_lad + S_η + S_drift **ÖZDEŞ** (kapanış
>    ≤ 7e−16). Gerçek gazda merdiven payı τ ≥ 0.62'de **%79 → %106**'ya
>    çıkıyor, η payı %21 → %0'a iniyor; drift 1e−6 (ölçülemez). S'nin
>    τ ≈ 0.47'deki işaret değişimini **merdiven kanalı** yapıyor.
>    **A4 ve keskin'de tam tersi:** merdiven payı bütün menzilde POZİTİF
>    kalıyor (+0.11 → +0.21), işaret değişimini **η kanalı** yapıyor.
>    Bu, gerçek ↔ sentetik arasında ilk kez ölçülen bir **mekanizma
>    farkıdır** (moment farkı değil).
> 8. **T3 — KALEM'in birebir formülü ÖLÜ, düzeltilmişi ÇALIŞIYOR.**
>    a_pred = d(A·S)/dτ | (A·S)'nin sıfır geçişinde = **−0.92 … −1.44**,
>    ölçülen a'nın **%9'u** (tam başarısızlık; formül 4π omurgasını
>    saymıyor). Buna karşılık
>    * **a_pred = 4π + d(A·S)/dτ|τ₀ : hata +7.4 / +7.0 / +3.6 / −4.6 / +9.9 %**
>      (gerçek-son / orta / keskin / A4 / P1), korel = +0.746;
>    * **a_pred = 4π + dM1/dτ|τ₀ : hata +3.2 / +2.8 / +3.0 / −3.6 / +0.7 %**,
>      korel = **+0.975** ve **gaz SIRALAMASI ölçülenle BİREBİR AYNI**
>      (gerçek-son < orta < keskin < P1 < A4).
>    * çıplak 4π (mekanizmasız): hata +18.4 / +17.9 / +7.8 / −3.7 / +4.8 %.
> 9. **H-b SINAVI GEÇTİ (işaret düzeyinde).** Kinematik omurga τ'da tam
>    doğrusal olduğundan **b'ye katkısı SIFIRDIR — b saf mekanizmadır**
>    (ölçüldü: |b(φ) − b(δ)| = 0.17 … 2.42). A·S'nin eğriliği gerçek gazda
>    **−3.83**, GUE-boyalı P1'de **+1.42**; M1'inki −5.63 / +4.78; ölçülen
>    −7.76 / +6.83. **İşaret dönmesi üretiliyor.** Beş gazın dördünde
>    b(A·S) ve b(M1) doğru işareti veriyor; tek başarısızlık **A4**
>    (ölçülen +5.54, öngörü −6.09 / −3.87).
> 10. **Ölçüm zinciri bit düzeyinde denetlendi.** 159'un φ, τ_eff, σ_φ ve
>    |Γ|'sı **altı (gaz, taban) çiftinde, 66 bantta 158'inkiyle 0.000e+00**;
>    bu raporun fit kodu 158'in W-A tablosunu (τ₀, a, b) beş gazda
>    yazıldığı üç haneye kadar yeniden üretiyor.

---

## 1. Yöntem

`159_configs/159_cekirdek.py` hiçbir ölçüm parçasını kopyalamaz; **import
eder**:

| kaynak | import edilen | ne |
|---|---|---|
| `155_cekirdek` | `eta_onbellek`, `izgara` | η zinciri (154'ün `eta_zinciri`'i), önbellek |
| `155_kos` | `veri_yukle` | gaz/pencere yükleyici (152/154/155'in `z_*.npy`'leri) |
| `156_cekirdek` | `zincir3`, `_bond` | ds = drift + lad + η ayrışımı, bond adımı |
| `154_cekirdek` | `pk_m` | çizgi listesi |

Bant/çizgi döngüsü `olc155` ile **birebir aynı konvansiyondadır**: aynı
aday listesi, aynı 220-örnekleme tohumu (21), aynı `gap < 2.5·dres`
filtresi, aynı ara-nokta (off-line) referansı `W' = w + gap/2`, aynı
8-grup round-robin jackknife. Doğrulama V1 bu iddiayı bit düzeyinde
sınıyor.

### İki bant ızgarası

| ad | ızgara | niçin |
|---|---|---|
| **t1** | `izgara(0.44, 0.80, 0.04)` → 9 bant, merkez 0.46 … 0.78 | **T1'in birincil taraması**; görevin istediği τ∈(0.42,0.80) menzili ve 8–10 bant. Taban 0.52'de lo ≥ 0.52 kalan **7 bant**. |
| **g158** | `izgara(0.28, 0.64, 0.02)` → 158'in ızgarası | **denetim + a/b bağı**: 158'in φ tablosunu bit düzeyinde yeniden üretir ve a/b fitini 158'in W-A penceresinde (τ̄ ∈ [0.43, 0.61]) koşar. |

### Koşulan ölçümler

12 koşu (`S_<veri>_t<taban>_<ızgara>.json`), **115 ölçülen bant**:
`son` (0.40/0.52 × t1/g158), `orta` (0.40 g158), `keskin` (0.40 ×
t1/g158), `A4` (0.40 × t1/g158; 0.52 g158), `P1` (0.40 × t1/g158).
Tipik koşu 0.7–2.4 dk; η önbellekleri 155'in `scratchpad/155/eta_*.npz`
dosyalarından yeniden kullanıldı (yeni önbellek üretilmedi).

### Yerel gücün tanımı — ve neden görevin verdiği tanım kullanılamadı

Görev metni: *"u_n = η_n·cos(ω m_n), v_n = η_n·sin(ω m_n); P_loc,n =
kayan-pencere ortalaması (u²+v²)"*. İki okuma var ve **ikisi de ölçüldü**:

* **Sözel okuma.** u² + v² ≡ η_n² **özdeş olarak**; kayan ortalaması
  ⟨η²⟩_W frekanstan bağımsızdır, yani "bant çizgilerinin YEREL gücü"
  olamaz — bandın bütün çizgilerine aynı sayıyı verir. Ölçüldü:
  A·S(⟨η²⟩_W64) = **+1.5e−4 … +2.6e−4**, τ ile neredeyse hiç değişmiyor
  (δ = −0.03 … −1.62 iken).
* **Demodüle okuma.** ĉ = (η e^{−iωm})'nin W-bond kayan ortalaması,
  P_dem = |ĉ|². Bu frekans-seçicidir ama ⟨P_dem⟩ = pow/4 + σ_η²/W, yani
  W-bond penceresine sızan geniş-bant η gürültüsünün **pedestalını**
  taşır. Ölçülen pedestal/koherent oranı (W=64): gerçek gaz **37.9 …
  6371**, A4 **101 … 10708**, keskin **166 … 1.4e6**. Sonuç:
  A·S(|ĉ|²) = **−0.003 … +0.005**, δ'nın binde birkaçı.

**Kullanılan tanım (BİRİNCİL).** Tahmincinin cebrinden ÖZDEŞ olarak çıkan
ağırlık:

    c_n = η_n·e^{−iW·m_n},   ⟨c⟩ = zc/2
    ρ_n = Re[c_{n+1}·conj⟨c⟩]        ⟨ρ⟩ = |⟨c⟩|² = pow/4   (TAM)
    S   = Cov(ρ, dsΔ)/⟨ρ⟩

ρ, çizginin **yerel koherent gücüdür**: küresel faz yönüne izdüşmüş
bond-başına güç payı; toplamı tam olarak çizginin ölçülen gücüdür.
Pencere parametresi taşımaz — ve pencere duyarlılığı yine ölçülebilir,
çünkü simetrik çekirdek için **Cov(K*ρ, X) = Cov(ρ, K*X)**: "ağırlığı
pencerelemek" ile "X'i pencerelemek" aynı sınavdır (§4).

### X'in tanımı — bond başına, gerekçesiyle

Bond n'de fazın ilerlemesi **tam olarak** A·(1 + X̃_n)'dir,
X̃_n = (m_{n+1}−m_n)·L/2π − 1. Yani tahmincinin üstelinde duran nicelik
**bond başına** adımdır; pencere ortalaması değil. Birincil X = dsΔ
(150/154/155'in değişkeni); X̃ ile farkı ölçüldü ve **≤ %0.35**
(§4). Pencere-ortalamalı X ayrı bir satır olarak tablolandı.

---

## 2. T0 — ölçülen φ_Γ'nın cebri: kinematik omurga + mekanizma fazı

### 2a. Özdeşlik

Tanımlar (n = 0 … N−2, m = bond orta noktaları):

    zc = 2⟨η_n·e^{−iW m_n}⟩          (kodun pow_on = |zc|²)
    zp = 2⟨η_{n+1}·e^{−iW m_n}⟩
    cr = zp·conj(zc)·e^{+i·2πW/L}    ← 148'den beri bütün zincirdeki satır

m_{n+1} − m_n = (2π/L)(1 + X̃_n) olduğu için, **hiçbir yaklaşım
yapmadan**:

    zp·conj(zc)·e^{−iA}/4 = ⟨(ρ_n + i·σ_n)·e^{+iA·X̃_n}⟩
    ρ + iσ ≡ c_{n+1}·conj⟨c⟩ ,  A = 2πW/L

Yani **ham korelatör deterministik ilerlemeyi e^{+iA} olarak zaten
taşıyor**. Onu çıkarmak için e^{−iA} gerekir; kod e^{+iA} ile çarpıyor:

> **φ_Γ = arg Γ = (2A − 2π) + δ = (4π·τ_eff − 2π) + δ,
> δ = arg[zp·conj(zc)·e^{−iA}]**

### 2b. Üç bağımsız delil

| | sınav | sonuç |
|---|---|---|
| **V2** | Özdeşliğin sayısal doğrulaması, gerçek veride 12 çizgi (τ = 0.30 … 0.84) | maks **bağıl** fark **7.8e−09** |
| **V3** | **Sentetik kontrol**: bilinen tek çizgi + bilinen adım dizisi (L=12, σ_adım=0.20), τ=0.60 | arg(raw·e^{−iA}) = **+0.0152**; arg(raw·e^{+iA}) = **+1.2719**; 2A−2π = **+1.2566**. Gürültü eklenince aynı (+0.0150 / +1.2716) |
| **V5** | Bant düzeyinde φ − (4π·τ_eff − 2π) = δ mi? | maks kalıntı **8.3e−04** (gerçek, g158) … **1.5e−02** (P1, g158) — bant genişliği mertebesi |

### 2c. Ölçülen ayrışım (gerçek gaz, taban 0.40, 158'in ızgarası)

| τ̄ | τ_eff | φ_meas | 4πτ_eff−2π | **δ** | \|Γ\| |
|---|---|---|---|---|---|
| 0.41 | 0.4095 | −1.1115 | −1.1378 | **+0.0262** | 0.922 |
| 0.43 | 0.4299 | −0.8765 | −0.8813 | **+0.0049** | 0.916 |
| 0.45 | 0.4508 | −0.6379 | −0.6178 | **−0.0198** | 0.912 |
| 0.47 | 0.4698 | −0.4258 | −0.3793 | **−0.0463** | 0.903 |
| 0.49 | 0.4908 | −0.1972 | −0.1154 | **−0.0817** | 0.892 |
| 0.51 | 0.5095 | **+0.0040** | +0.1192 | **−0.1150** | 0.885 |
| 0.53 | 0.5306 | +0.2248 | +0.3841 | **−0.1592** | 0.870 |
| 0.55 | 0.5506 | +0.4267 | +0.6357 | **−0.2083** | 0.857 |
| 0.57 | 0.5703 | +0.6270 | +0.8837 | **−0.2561** | 0.844 |
| 0.59 | 0.5899 | +0.8169 | +1.1292 | **−0.3121** | 0.831 |
| 0.61 | 0.6098 | +0.9899 | +1.3794 | **−0.3887** | 0.797 |
| 0.63 | 0.6297 | +1.1760 | +1.6304 | **−0.4542** | 0.787 |

korel(φ, omurga) = **0.999246**; sd(φ) = 0.7188, sd(omurga) = 0.8674,
sd(δ) = 0.1515; medyan |δ|/|φ| = **0.389**.

### 2d. Sonuç: a ve τ₀'ın anatomisi

Omurga **τ'da tam doğrusaldır** (eğim 4π = 12.5664, eğrilik SIFIR).
Dolayısıyla

* **a = dφ/dτ|τ₀ = 4π + dδ/dτ|τ₀**, ve
* **b = eğrilik(φ) = eğrilik(δ)** — b **saf mekanizmadır**.

Ölçülen (g158, taban 0.40, 158'in W-A penceresi):

| gaz | a(φ) = 158'in sayısı | dδ/dτ | **4π + dδ/dτ** | fark |
|---|---|---|---|---|
| gerçek-son | 10.615 | −1.990 | **10.577** | 0.038 |
| gerçek-orta | 10.656 | −1.981 | **10.585** | 0.071 |
| keskin | 11.660 | −0.882 | **11.685** | 0.025 |
| A4 | 13.054 | +0.584 | **13.151** | 0.097 |
| P1 (GUE) | 11.988 | −0.559 | **12.007** | 0.019 |

korel(a_meas, 4π+dδ/dτ) = **+0.9999**, sıralama aynı.

**τ₀ AYRIŞIMI.** Omurga tek başına τ₀ = ½ verir; gaz-bağımlı kayma
−δ(½)/a olmalıdır:

| gaz | τ₀ (ölçülen) | δ(τ=½) | ½ − δ(½)/a | fark |
|---|---|---|---|---|
| gerçek-son | 0.5088 | −0.0955 | **0.5090** | −0.0002 |
| gerçek-orta | 0.5087 | −0.0966 | **0.5091** | −0.0003 |
| keskin | 0.4886 | +0.1352 | **0.4884** | +0.0002 |
| A4 | 0.4897 | +0.1364 | **0.4895** | +0.0002 |
| P1 (GUE) | 0.4996 | −0.0087 | **0.5007** | −0.0011 |

**155/157/158'in τ₀ merdiveni (0.4886 … 0.5088), ½ ile δ(½)/a'nın
toplamıdır.** Merdivenin toplam yüksekliği 0.0202; δ(½)/a'nın menzili
0.0207. Yani merdiven **tamamen** mekanizma fazının τ=½'deki değerinden
geliyor.

---

## 3. T1 — mekanizma denklemi, bant bant

### 3a. Gerçek gaz, taban 0.40, ızgara t1 (9 bant, hepsi sağlıklı)

`M1 = arg⟨ρ·e^{iA·dsΔ}⟩/⟨ρ⟩` — mekanizmanın **bütün mertebeli** hâli
(σ kanalı atılmış, birinci-mertebe kesme YOK). `ReΓ_arın` = kinematikten
arındırılmış korelatörün reel kısmı (sönüm).

| τ_eff | φ_meas | A·S | **φ/(A·S)** | δ_meas | **δ/(A·S)** | M1 | **δ/M1** | **M1/(A·S)** |
|---|---|---|---|---|---|---|---|---|
| 0.4607 | −0.5287 | +0.0137 | −38.71 | −0.0335 | −2.450 | −0.0120 | +2.795 | −0.88 |
| 0.4997 | −0.1019 | −0.0237 | +4.30 | −0.0975 | +4.118 | −0.0631 | +1.544 | +2.67 |
| 0.5392 | +0.3115 | −0.0731 | −4.26 | −0.1803 | **+2.465** | −0.1313 | **+1.373** | 1.80 |
| 0.5791 | +0.7116 | −0.1350 | −5.27 | −0.2811 | **+2.082** | −0.2170 | **+1.295** | 1.61 |
| 0.6188 | +1.0729 | −0.2064 | −5.20 | −0.4179 | **+2.025** | −0.3262 | **+1.281** | 1.58 |
| 0.6578 | +1.3875 | −0.3046 | −4.56 | −0.5909 | **+1.940** | −0.4640 | **+1.274** | 1.52 |
| 0.6980 | +1.6440 | −0.4212 | −3.90 | −0.8422 | **+1.999** | −0.6442 | **+1.307** | 1.53 |
| 0.7382 | +1.8626 | −0.6058 | −3.07 | −1.1282 | **+1.862** | −0.8683 | **+1.299** | 1.43 |
| 0.7766 | +1.8604 | −0.8874 | −2.10 | −1.6189 | **+1.824** | −1.1744 | **+1.378** | 1.32 |

**Oran özeti (ideal 1.00, MÜHÜR eşiği ±%25 ⇒ 0.75–1.25):**

| oran | ortalama | medyan | menzil | ±%25 içinde |
|---|---|---|---|---|
| **φ/(A·S)** — KALEM birebir | −6.973 | −4.259 | [−38.7, +4.30] | **0/9** |
| **δ/(A·S)** — kinematikten arınmış | +1.763 | +1.999 | [−2.45, +4.12] | 0/9 |
| **δ/M1** — ρ, bütün mertebeler | +1.505 | +1.307 | [+1.27, +2.80] | 0/9 |
| δ/(A·S/ReΓ) — sönüm düzeltmeli | +1.027 | +1.194 | [−2.23, +3.65] | 2/9 |

İlk iki bant (τ_eff = 0.46, 0.50) δ ve A·S'nin **ikisinin de sıfırdan
geçtiği** bölgededir; oran orada tanımsıza yakındır ve menzili tek başına
o iki bant açıyor. **τ_eff ≥ 0.54'ün yedi bandında δ/(A·S) = 1.82 – 2.47,
δ/M1 = 1.27 – 1.38** — yani açık gürültü değil, **kararlı bir çarpan**.

### 3b. Açığın anatomisi (ölçüldü, varsayılmadı)

δ = arg⟨(ρ+iσ)·e^{iAX}⟩ özdeşliğinden üç basamak:

| basamak | ne atılıyor | gerçek gazda etkisi (τ ≥ 0.54) |
|---|---|---|
| δ → M1 | **kuadratür (σ) kanalı** — c_{n+1}'in küresel faza dik bileşeni | **×1.27 – 1.38** |
| M1 → A·S | **birinci-mertebe kesme** (açılım parametresi A·σ_X = **0.68 … 1.14**, küçük DEĞİL) | **×1.32 – 1.80** |
| çarpım | | **×1.82 – 2.47** |

Yani mekanizma denklemi "biraz kaçırmıyor"; **iki tanımlı yerde**, her
biri ölçülebilir büyüklükte kaçırıyor. σ kanalı, KALEM'in reel-ağırlık
ansatzında **yeri olmayan** bir kanaldır: yerel genliğin küresel faza
göre **kayması** ile adımın korelasyonudur.

### 3c. Diğer gazlar (t1, taban 0.40; ✗ = sağlıksız bant, §7)

| gaz | sağlıklı bant | δ/(A·S) medyan (menzil) | ±%25 | δ/M1 medyan (menzil) | ±%25 |
|---|---|---|---|---|---|
| **gerçek-son** | 9/9 | +1.999 (−2.45 … +4.12) | 0/9 | +1.307 (+1.27 … +2.80) | 0/9 |
| **keskin** | 8/9 | **+1.116** (−7.83 … +1.83) | **3/8** | **+1.199** (+0.66 … +2.71) | **4/8** |
| **A4** | 5/9 | +1.816 (−5.25 … +4.92) | 0/5 | +1.893 (+1.45 … +57.3) | 0/5 |
| **P1 (GUE)** | 6/9 | +0.325 (−0.22 … +0.76) | 1/6 | +0.073 (−1.68 … +4.24) | 0/6 |
| **gerçek-son, taban 0.52** | 7/7 | +2.018 (+1.93 … +2.71) | 0/7 | **+1.267** (+1.23 … +1.38) | **3/7** |

**Keskin gazda mekanizma denklemi neredeyse tutuyor** (τ_eff = 0.4608 …
0.7393'te δ/(A·S) = 1.30/1.21/1.01/0.53/–/1.83/1.29/1.02). Keskin, |Γ|'sı
menzil boyunca **düz** kalan tek gazdır (0.78 – 0.96); yani yüksek
mertebeler orada küçüktür ve birinci-mertebe kesme daha az kaybettirir.
**P1 tersi uçta:** |Γ| = 0.40 – 0.52 (bütün gazların en sönümlüsü), faz
neredeyse hiç ilerlemiyor (δ = −0.03 … +0.17, altı sağlıklı bant boyunca)
ve A·S onu **1.3 – 10 kat aşıyor**. Yani P1'de mekanizma denklemi
fazı **fazla** öngörüyor; gerçek gazda ise **eksik**.

### 3d. Taban dayanıklılığı (0.40 ↔ 0.52, örtüşen bantlar)

Gerçek gaz, t1'in 7 örtüşen bandında δ/(A·S): taban 0.40'ta 2.465 /
2.082 / 2.025 / 1.940 / 1.999 / 1.862 / 1.824; taban 0.52'de 2.708 /
2.377 / 2.171 / 2.018 / 1.966 / 1.925 / 1.965. **Açığın büyüklüğü taban
konvansiyonuna asılı değil** (medyan 2.00 vs 2.02; en büyük fark %10).
g158 ızgarasında da aynı (1.96 – 2.67 vs 2.10 – 2.82).

---

## 4. T1b — tanım, pencere ve ara-nokta referansı sistematikleri

### 4a. P_loc tanımı (gerçek gaz, t1, taban 0.40)

| τ_eff | δ | **A·S(ρ)** | A·S(\|ĉ\|², W=32) | (W=64) | (W=128) | A·S(⟨η²⟩_64) | pedestal(W=64) |
|---|---|---|---|---|---|---|---|
| 0.4607 | −0.0335 | **+0.0137** | +0.00023 | −0.00021 | −0.00006 | +0.00015 | 37.9 |
| 0.5392 | −0.1803 | **−0.0731** | +0.00040 | +0.00101 | +0.00082 | +0.00018 | 52.0 |
| 0.6188 | −0.4179 | **−0.2064** | +0.00030 | +0.00328 | +0.00187 | +0.00021 | 6371.0 |
| 0.6980 | −0.8422 | **−0.4212** | +0.00165 | +0.00279 | +0.00347 | +0.00023 | 459.2 |
| 0.7766 | −1.6189 | **−0.8874** | −0.00319 | +0.00036 | +0.00482 | +0.00026 | 1193.0 |

Demodüle P_loc, δ'nın **%0.02 – %1.0**'i kadar; işareti bile kararlı
değil. ⟨η²⟩_W ise τ ile neredeyse hiç değişmiyor. **Pedestal** sütunu
nedenidir: ⟨|ĉ|²⟩/(pow/4) = 1 + (σ_η²/W)/(pow/4), ve tek bir çizginin
gücü W=64 penceresine sızan η gürültüsünün 1/38 … 1/6371'i kadardır.
Pedestalı yenmek için W ≳ 64·pedestal ≈ **2.4·10³ – 4.1·10⁵ bond**
gerekirdi; o ölçekte de X'in dalgalanması ortalamayla siliniyor (§4b) —
yani **pencere-tabanlı yerel güç ile bu bağlaşım ölçülemez**, ρ ile
ölçülür. (Demodüle A·S / δ oranı: W=32'de %0.02–1.00, W=64'te
%0.02–0.79, W=128'de %0.17–0.56.)

### 4b. Kayan pencere duyarlılığı (ağırlığı pencerelemek ≡ X'i pencerelemek)

| τ_eff | A·S(X bond) | A·S(X̄_32) | A·S(X̄_64) | A·S(X̄_128) | X̄_64 / bond |
|---|---|---|---|---|---|
| 0.4997 | −0.0237 | +0.00017 | +0.00006 | −0.00006 | −0.0025 |
| 0.5791 | −0.1350 | −0.00007 | +0.00053 | +0.00103 | −0.0039 |
| 0.6578 | −0.3046 | +0.00003 | +0.00107 | +0.00201 | −0.0035 |
| 0.7382 | −0.6058 | −0.00058 | +0.00203 | +0.00378 | −0.0033 |
| 0.7766 | −0.8874 | −0.00132 | +0.00340 | +0.00544 | −0.0038 |

|A·S(X̄_W)/A·S(X_bond)| menzili: **W=32: %0.01 – %0.74**, **W=64:
%0.13 – %1.18**, **W=128: %0.26 – %2.72** (τ_eff ≥ 0.54'te sırasıyla
%0.01–0.27, %0.13–0.39, %0.40–0.76). **Bağlaşım bond ölçeğindedir**;
32–128 bond bandında bileşeni yok. (Bu, KALEM'in "kayan pencere ~64
bond" önerisinin ölçülmüş reddidir — ve kendi başına bir sonuçtur.)

### 4c. Ara-nokta (off-line) referansının sistematiği

S için de S_on − S_off uygulandı (birincil). Referanssız (yalnız S_on)
sürümle farkı:

| τ_eff | A·S(on−off) | A·S(yalnız on) | fark | % |
|---|---|---|---|---|
| 0.4607 | +0.01366 | +0.01339 | −0.00027 | −2.0 |
| 0.4997 | −0.02368 | −0.02468 | −0.00100 | −4.2 |
| 0.5791 | −0.13502 | −0.14238 | −0.00736 | −5.5 |
| 0.6578 | −0.30460 | −0.32840 | −0.02380 | −7.8 |
| 0.6980 | −0.42122 | −0.45683 | −0.03561 | **−8.5** |
| 0.7766 | −0.88744 | −0.89592 | −0.00848 | −1.0 |

Gerçek gazda sistematik **%1 – %8.5**, hep aynı yönde (referanssız S
büyütüyor). δ/(A·S) ≈ 2 açığının **%8'inden azını** açıklıyor. A4'te
sistematik büyük (%0.5 – %27 sağlıklı bantlarda, bozuk bantlarda %40).

### 4d. X değişkeni: dsΔ mı, X̃ mi

| τ_eff | A·S(dsΔ) | A·S(X̃) | fark |
|---|---|---|---|
| 0.4607 | +0.01366 | +0.01361 | −0.35% |
| 0.5791 | −0.13502 | −0.13517 | −0.11% |
| 0.7766 | −0.88744 | −0.88786 | −0.05% |

İkisi **korel = 0.999877**, maks nokta farkı 1.2e−02. dsΔ (KALEM'in
değişkeni) birincil olarak kullanıldı; seçim hiçbir hükmü değiştirmiyor.

---

## 5. T2 — S'nin kanal ayrışımı

Kovaryans doğrusal, `_bond` doğrusal ve ds = drift + lad + η **özdeş**
olduğundan

    S_tam = S_lad + S_η + S_drift        (ÖZDEŞ, model seçimi DEĞİL)

Kapanış ölçüldü: maks |S_tam − Σ| = **6.7e−16** (115 bandın hepsinde
≤ 7e−16). Bu, 156'nın R-tabanlı ayrışımından **daha keskin bir araçtır**:
orada kanallar ayrı ayrı fit edilip aşırı-belirleme sınavına sokuluyordu;
burada ayrışım tamdır ve pay doğrudan okunur.

### 5a. Gerçek gaz (bond kovaryans payı: lad %92.79, η %7.21, drift %0.00)

| τ_eff | A·S_tam | A·S_lad | A·S_η | A·S_drift | lad % | η % |
|---|---|---|---|---|---|---|
| 0.4607 | +0.01366 | **+0.07240** | −0.05874 | +2.6e−07 | +530 | −430 |
| 0.4997 | −0.02368 | **+0.03521** | −0.05889 | +3.2e−07 | −149 | +249 |
| 0.5392 | −0.07313 | **−0.01782** | −0.05531 | +4.8e−07 | +24 | +76 |
| 0.5791 | −0.13502 | **−0.08243** | −0.05259 | +8.4e−07 | +61 | +39 |
| 0.6188 | −0.20643 | **−0.16376** | −0.04267 | +1.2e−06 | +79 | +21 |
| 0.6578 | −0.30460 | **−0.26967** | −0.03493 | +1.8e−06 | +89 | +12 |
| 0.6980 | −0.42122 | **−0.40922** | −0.01200 | +1.4e−06 | +97 | +3 |
| 0.7382 | −0.60579 | **−0.60390** | −0.00190 | +1.9e−06 | +100 | +0 |
| 0.7766 | −0.88744 | **−0.94242** | +0.05497 | +3.3e−06 | +106 | −6 |

**156 ile tam tutarlı:** yüksek τ'da bağlaşımı taşıyan tek kanal
**merdivendir** (τ ≥ 0.70'te η payı ≤ %3). Drift kanalı 1e−6 —
"drift-küçük" varsayımı burada da ölçüldü ve tuttu.

**Yeni ve 156'da olmayan kayıt:** S_tam'ın τ ≈ 0.47'deki **işaret
değişimini merdiven kanalı yapıyor** (S_lad: +0.0724 → +0.0352 →
−0.0178 → −0.0824), η kanalı ise o bölgede işaret değiştirmeden
−0.053 … −0.059 civarında duruyor. Yani **φ_Γ'nın sıfır geçişini —
dolayısıyla τ₀'ı — merdiven kanalı belirliyor.**

### 5b. Sentetikler — desen TERSİNE dönüyor

| gaz | bond payı lad / η | A·S_lad'ın davranışı | işaret değişimini yapan |
|---|---|---|---|
| **gerçek** | %92.8 / %7.2 | +0.072 → **−0.942** (işaret değiştiriyor) | **merdiven** |
| **A4** | %76.4 / %23.6 | +0.114 → +0.208; sağlıklı bantların **hepsinde pozitif** | **η** (−0.033 → −0.270) |
| **keskin** | %39.8 / %60.2 | +0.105 → +0.204 → +0.185; **hiç işaret değiştirmiyor** | **η** (+0.023 → −0.488) |

Bu, gerçek ile sentetik arasında **moment farkı değil, mekanizma farkı**
olan ilk ölçülmüş nesnedir: gerçek gazda merdivenin yerel gücü ile yerel
adım arasındaki korelasyon τ ile işaret değiştiriyor; sentetiklerde
değiştirmiyor ve sıfır geçişi η kanalından geliyor. (keskin'de bond
payının **çoğunluğu** η'dadır — taban 0.40'ta %60.2; 156 taban 0.52'de
%55.1 ölçmüştü, fark taban konvansiyonudur — ama 156'nın hükmü orada da
merdiveni gösteriyordu. Burada ölçülen ise keskin'in η kanalının S'de
baskın ve gerçeğinkinden **çok daha büyük** olduğudur.)

---

## 6. T3 — a_pred ve b'nin işareti

Fit konvansiyonu 158/157'nin W-A'sı: apsis τ_eff, ağırlık 1/σ_jk,
`polyfit(·,·,2)`, τ₀ = köklerden x ortalamasına en yakını, a = c₁+2c₀τ₀,
b = c₀; pencere τ̄ ∈ [0.43, 0.61] (g158, taban 0.40, 10 bant).

### 6a. Fit kodunun denetimi

| gaz | a(φ) bu koşu | a — 158 W-A | b(φ) bu koşu | b — 158 W-A |
|---|---|---|---|---|
| gerçek-son | **10.615** | 10.615 | **−7.76** | −7.759 |
| gerçek-orta | **10.656** | 10.656 | **−7.80** | −7.801 |
| keskin | **11.660** | 11.660 | **−3.81** | −3.811 |
| A4 | **13.054** | 13.054 | **+5.54** | +5.541 |
| P1 | **11.988** | 11.988 | **+6.83** | +6.835 |

### 6b. KALEM'in birebir T3'ü — ÖLÜ

a_pred = d(A·S)/dτ, **(A·S)'nin kendi sıfır geçişinde**:

| gaz | τ₀(A·S) | d(A·S)/dτ | a_meas | **a_pred/a_meas** |
|---|---|---|---|---|
| gerçek-son | 0.4767 | −0.916 | 10.615 | **−0.086** |
| gerçek-orta | 0.4751 | −0.922 | 10.656 | −0.087 |
| keskin | 0.6183 | −1.350 | 11.660 | −0.116 |
| A4 | 0.5982 | −1.435 | 13.054 | −0.110 |
| P1 | — (kök yok) | — | 11.988 | — |

Formül a'nın **%9'unu** veriyor ve işareti ters. Sebebi §2: formül 4π
omurgasını saymıyor. Ayrıca A·S'nin sıfır geçişi φ'nin τ₀'ından
sistematik olarak farklıdır (gerçekte 0.4767 vs 0.5088; keskin'de 0.6183
vs 0.4886) — yani "sıfır geçişinde türev" reçetesi gazdan gaza farklı
bir noktada okuyor.

### 6c. Düzeltilmiş T3 — a_pred = 4π + d(mekanizma)/dτ|τ₀(φ)

| gaz | **a_meas** | 4π+dδ/dτ | **4π+dM1/dτ** | hata % | **4π+d(A·S)/dτ** | hata % | çıplak 4π | hata % |
|---|---|---|---|---|---|---|---|---|
| gerçek-son | 10.615 | 10.577 | **10.956** | **+3.2** | 11.405 | +7.4 | 12.566 | +18.4 |
| gerçek-orta | 10.656 | 10.585 | **10.958** | **+2.8** | 11.406 | +7.0 | 12.566 | +17.9 |
| keskin | 11.660 | 11.685 | **12.008** | **+3.0** | 12.083 | +3.6 | 12.566 | +7.8 |
| A4 | 13.054 | 13.151 | **12.588** | **−3.6** | 12.453 | −4.6 | 12.566 | −3.7 |
| P1 (GUE) | 11.988 | 12.007 | **12.076** | **+0.7** | 13.170 | +9.9 | 12.566 | +4.8 |
| **korel(a_meas, ·)** | — | **+0.9999** | **+0.9748** | | +0.7463 | | — | |
| **gaz sıralaması** | — | **AYNI** | **AYNI** | | FARKLI (A4↔P1) | | — | |

> **Gerçek gazın referansı** (158, 32 fit): a = **10.713 ± 0.059**,
> σ_konv = 0.075. M1 öngörüsü (10.956) referanstan **+0.24 = 3.2 σ_konv**;
> A·S öngörüsü (11.405) **+0.69 = 9.2 σ_konv**; çıplak 4π **+1.85 =
> 24.7 σ_konv**. **Hiçbiri ölçüm hassasiyetine ulaşmıyor** — ama sıralama
> ve mertebe doğru.

**Sentetik izleme sınavı (görevin istediği):** üç sentetik gazda
(keskin, A4, P1) a_pred, a_meas'i izliyor mu? **M1 tabanlı öngörü beş
gazın sıralamasını BİREBİR üretiyor** (gerçek-son < orta < keskin < P1 <
A4) ve hepsinde hata ≤ %3.6. A·S tabanlı öngörü sıralamada A4 ile P1'i
takas ediyor.

### 6d. b'nin işareti — H-b sınavı

Kinematik omurga τ'da tam doğrusal olduğundan **b'ye hiç katkı vermez**;
ölçüldü: |b(φ) − b(δ)| = 0.17 (P1) … 2.42 (A4), tipik 0.35 – 1.8 (fark
tümüyle bant-düzeyi arındırmanın kalıntısıdır).

| gaz | b(φ) ölçülen | b(δ) | **b(M1)** | işaret | **b(A·S)** | işaret |
|---|---|---|---|---|---|---|
| **gerçek-son** | **−7.76** | −6.50 | **−5.63** | ✓ | **−3.83** | ✓ |
| gerçek-orta | −7.80 | −6.00 | −5.42 | ✓ | −3.54 | ✓ |
| keskin | −3.81 | −4.16 | −3.72 | ✓ | −3.34 | ✓ |
| A4 | +5.54 | +3.12 | −3.87 | ✗ | −6.09 | ✗ |
| **P1 (GUE)** | **+6.83** | +7.00 | **+4.78** | ✓ | **+1.42** | ✓ |

> **H-b: GEÇTİ (işaret düzeyinde).** KALEM'in hipotezi "gerçek gazın
> ters-nefesi ↔ GUE-boyalının adyabatik tepkisi ⇒ S'nin τ-eğriliği
> işaret değiştirir" idi. Ölçülen: **b(A·S) = −3.83 (gerçek) / +1.42
> (P1)**, **b(M1) = −5.63 / +4.78**, ölçülen b = −7.76 / +6.83.
> **S'nin eğriliği işareti gerçekten döndürüyor** ve dönüş yönü doğru.
> Büyüklük tutmuyor (%50 – %80 eksik) — beklenen, çünkü §3b'nin iki
> kesme faktörü b'ye de giriyor.
>
> **Tek karşı-örnek A4 ve saklanmıyor:** ölçülen b(A4) = +5.54 dışbükey,
> ama hem A·S hem M1 içbükey (−6.09 / −3.87) öngörüyor. 158 zaten
> sentetiklerde b'nin "ölçülemediğini" (taban yayılımı %84 – %474)
> yazmıştı; A4'ün b'si o listenin en oynak üyelerindendi (16 fitte
> |ort|/sd = 1.10). Yine de bu bir tutmama kaydıdır.

---

## 7. Denetim

| | sınav | sonuç |
|---|---|---|
| **V1** | 159'un φ, τ_eff, σ_φ, \|Γ\|'sı = 158'in kaydı mı? | **6 (gaz, taban) çifti, 66 ortak bant: maks \|Δφ\| = \|Δτ_eff\| = \|Δσ_φ\| = \|Δ\|Γ\|\| = 0.000e+00** (bit düzeyinde özdeş) |
| **V1b** | Bu raporun fit kodu 158'in W-A tablosunu üretiyor mu? | beş gazın beşinde de a ve b, 158'in yazdığı üç haneye kadar aynı (§6a) |
| **V2** | `zp·conj(zc)·e^{−iA}/4 = ⟨(ρ+iσ)e^{+iAX̃}⟩` özdeşliği | 12 çizgide maks **bağıl** fark **7.8e−09** |
| **V3** | Sentetik kontrol (bilinen çizgi, bilinen adım) | e^{−iA}: **+0.0152**; e^{+iA}: **+1.2719**; 2A−2π = +1.2566 |
| **V4** | 156'nın `zincir3` η'sı = 155'in önbelleklediği η mi? | maks fark **3.3e−16**; bond kapanışı (X_tam − ΣX_k) **3.3e−16** |
| **V5** | φ − (4πτ_eff − 2π) = δ mi (bant düzeyi)? | maks kalıntı 8.3e−04 (gerçek, g158) … 1.5e−02 (P1, g158) |
| **V6** | S_tam = S_lad + S_η + S_drift kapanışı | 115 bandın hepsinde ≤ **6.7e−16** |

### Bant sağlığı kuralı

Bir bant "sağlıklı" sayıldı: |Γ| ≤ 1.0 (Γ normalize korelasyondur),
|φ| ≤ 2.8 (π'ye sarma payı), |τ_eff − τ̄| ≤ yarım bant; ve **ilk
bozulmadan sonraki bantlar da düştü** (tahminci bir kez bozulunca üst
bantlar bozulmayı devralır). Sonuç: gerçek 9/9 ve 7/7, keskin 8/9,
P1 6/9, A4 5/9. **Sentetikler t1 ızgarasında τ ≈ 0.62–0.66'nın üstünde
bozuluyor** (A4: τ̄=0.66'da φ π'yi geçip sarıyor, τ̄=0.70'te |Γ| = 1.27;
P1: τ̄=0.70'te |Γ| = 2.71, τ̄=0.74'te 1.56). Bozuk bantlar tablolarda
duruyor, hiçbir özete girmedi.

---

# HÜKÜM

## (i) T1: mekanizma denklemi yazıldığı hâliyle YANLIŞ — ama eksiği ölçülebilir ve tam

`φ_Γ(τ) = A·S(τ)` gerçek gazda dokuz bandın hiçbirinde ±%25 içinde
değil (oran −38.7 … +4.3). Sebebi ölçüldü ve iki parçadır:

1. **Ölçülen φ_Γ, denklemin sol tarafı değildir.** φ_Γ = (4πτ_eff − 2π) + δ;
   denklem δ hakkındadır. Omurga, tahmincinin de-rotasyon çarpanının
   deterministik faz ilerlemesini çıkarmak yerine ikiye katlamasından
   gelir (V2/V3 ile kanıtlı, 148'den beri 11 dosyada aynı satır).
2. **δ ile A·S arasında bile kararlı bir ×2 açık var** ve ikiye ayrılıyor:
   kuadratür (σ) kanalı ×1.27–1.38, birinci-mertebe kesme ×1.32–1.80.
   İkisi de KALEM'in reel-ağırlık + küçük-faz ansatzının dışındadır.

**MÜHÜR YOK.** Ama denklem "kaba tahmin" seviyesinde değil: kinematikten
arındırıldıktan sonra A·S, δ'yı **her bantta doğru işaret ve 1.8–2.5 kat
içinde** veriyor; keskin gazda **1.0–1.3 kat** içinde.

## (ii) Bu koşunun asıl çıktısı: a ve τ₀'ın anatomisi

> **a = 4π + dδ/dτ|τ₀   (beş gazda hata ≤ 0.10, korel +0.9999)**
> **τ₀ = ½ − δ(½)/a     (beş gazda hata ≤ 0.0011)**

Bu, 155/157/158'in üç raporluk sayısal mirasını yeniden konumlandırıyor:

* **a = 10.713 ± 0.059'un 12.566'sı (4π) kinematiktir** ve bütün gazlarda
  aynıdır. 158'in "a makine-evrenseli DEĞİL" hükmü ayakta — ama
  evrensel-olmayan kısım a'nın kendisi değil, **a − 4π = −1.853 ± 0.059
  (gerçek)**'dir. Gazlar arası menzil (10.71 … 13.07) bu kalıntının
  menzilidir: −1.85 … +0.51. **Türetilecek sayı a değil, a − 4π'dir**;
  ve o sayının bağıl belirsizliği (%3.2) a'nınkinden (%0.55) **yaklaşık
  altı kat büyüktür** — hedef daha zor, ama tanımı çok daha temiz.
* **τ₀ merdiveni (0.4886 … 0.5088) ½ etrafında bir mekanizma kaymasıdır.**
  Merdivenin toplam yüksekliği 0.0202, −δ(½)/a'nın menzili 0.0207.
  155'in "τ₀ korelasyon merdivenini izliyor" bulgusu geçerliliğini
  koruyor; ama ölçtüğü şeyin adı artık belli: **δ(½)/a**.
* **b saf mekanizmadır** (omurga doğrusal, eğriliği sıfır). 158'in
  "b gerçekte 16/16 negatif, P1'de 16/16 pozitif" mührü, S'nin
  τ-eğriliğinin işaret değiştirmesinin doğrudan izidir (§6d).

## (iii) T2: ayrışım TAM ve 156'yı doğruluyor; sentetiklerde desen ters

S_tam = S_lad + S_η + S_drift **özdeş** (≤ 7e−16). Gerçek gazda
merdiven kanalı τ ≥ 0.70'te S'nin **%97–%106**'sını taşıyor, η kanalı
%0–3; drift ölçülemez (1e−6). **156'nın "bağlaşım merdivende yaşıyor"
hükmü, R-fitine hiç girmeden, doğrudan bir kovaryans ayrışımıyla
doğrulandı.**

Yeni: gerçek gazda S'nin işaret değişimini (dolayısıyla τ₀'ı) **merdiven
kanalı** yapıyor; A4 ve keskin'de merdiven kanalı hiç işaret değiştirmiyor
ve sıfır geçişi **η kanalından** geliyor. Bu, gerçek ↔ sentetik farkının
**mekanizma düzeyinde** ilk ölçümüdür.

## (iv) T3: KALEM'in reçetesi ölü, düzeltilmişi %3'te

`a_pred = d(A·S)/dτ|sıfır-geçiş` ölçülen a'nın %9'unu veriyor (işareti de
ters). Düzeltilmiş reçete `a_pred = 4π + d(mekanizma)/dτ|τ₀(φ)`:
**M1 ile beş gazda hata ≤ %3.6, korel +0.975, gaz sıralaması ölçülenle
birebir aynı**; A·S ile hata ≤ %9.9, korel +0.746, sıralamada A4↔P1
takası. Gerçek gazda M1 öngörüsü 10.956, referans 10.713 ± 0.059 —
**3.2 σ_konv**, yani **ölçüm hassasiyetine ulaşmıyor** ama mertebeyi
ve gaz sırasını yakalıyor.

**H-b geçti (işaret):** b(A·S) gerçekte −3.83, GUE-boyalı P1'de +1.42;
b(M1) −5.63 / +4.78. İşaret dönmesi üretiliyor. **A4 karşı-örnek**
(ölçülen +5.54, öngörü negatif) ve kayda geçiyor.

## Sıradaki adım (bu ölçümün işaret ettiği)

1. **De-rotasyon çarpanının işareti kalemle karara bağlanmalı.** Eğer
   niyet 148'in docstring'indeki gibi "bilinen faz ilerlemesini çıkarmak"
   ise, çarpan `np.exp(−1j*TWO_PI*W/L)` olmalı; o zaman ölçülen ilkel
   gözlenebilir **δ**'dır ve tayfı bu raporda ölçülüdür (a → dδ/dτ =
   −1.99 gerçek, −0.88 keskin, +0.58 A4, −0.56 P1; τ₀(δ) = 0.434 gerçek,
   0.597 keskin, kök yok A4, 0.507 P1). **Bu, 155/157/158'in bütün
   sayılarının yeniden koşulmasını gerektirir** ve bu raporda
   YAPILMAMIŞTIR — burada ölçülen yalnız ayrışımdır.
2. **Kuadratür (σ) kanalı yeni bir nesnedir ve adı yok.** δ'nın %23–30'unu
   taşıyor, gerçek gazda kararlı (δ/M1 = 1.27–1.38, iki tabanda da).
   Yerel genliğin küresel faza göre kayması ile adımın korelasyonudur;
   150'nin reel-ağırlık ailesinde karşılığı yok. Ölçülmesi ucuz (bu
   koşunun çıktılarında hazır), yorumu açık.
3. **Birinci-mertebe kesme kapatılabilir.** A·σ_X = 0.68 … 1.14; M1 (bütün
   mertebeler, ρ ağırlığı) zaten hesaplanıyor ve a'yı %3'te veriyor.
   Kapalı-form adayı artık `S` değil, **⟨ρ e^{iAX}⟩/⟨ρ⟩'nun kümülant
   açılımı**dır — ikinci ve üçüncü kümülantlar (⟨ρX²⟩, ⟨ρX³⟩) bu
   koşunun verisinden çıkarılabilir.
4. **A4'ün b'si ayrı bir kalemdir.** Beş gazın dördünde işaret tutuyor,
   A4'te tutmuyor — ve A4, 158'in b tablosunda da en oynak gazdı
   (|ort|/sd = 1.10). Sentetiklerde b'nin ölçülebilir olup olmadığı hâlâ
   açık; bu koşu onu çözmedi.

---

## Dürüstlük notları

* **Uydurma yok; her sayı yeniden koşuldu.** Bütün tablo değerleri
  `159_configs/159_kos.py`'nin ürettiği 12 JSON'dan ve
  `159_configs/159_analiz.py`'nin çıktısından
  (`scratchpad/159/analiz_cikti.txt`) otomatik alındı. **Tablodan okunan
  tek şey** 158'in referanslarıdır (a = 10.713 ± 0.059, b = −7.22 ± 0.79,
  σ_konv = 0.075) ve her biri kaynağıyla işaretli. 158'in W-A satırları
  ise okunmadı, **yeniden hesaplandı** ve aynı çıktı (§6a).
* **Görevin verdiği P_loc tanımı kullanılamadı ve bu bir seçim değil,
  bir ölçümdür.** Sözel okuma (⟨η²⟩_W) frekans-kördür — matematiksel
  olarak bandın çizgilerini ayıramaz. Demodüle okuma (|ĉ|²) frekans-
  seçicidir ama gürültü pedestalı koherent gücün 38–6371 katıdır
  (ölçüldü, §4a). Kullanılan ρ ağırlığı **keyfi bir alternatif değil**:
  mekanizma denklemini tahmincinin cebrinde ÖZDEŞ yapan tek ağırlıktır
  (V2, bağıl 7.8e−09) ve ortalaması tam olarak çizgi gücüdür.
* **"İşaret hatası" nitelemesi bir yorumdur; ÖLÇÜLEN şey ayrışımdır.**
  Ölçülen: (a) ham korelatör e^{+iA} taşır (V2, V3), (b) kod e^{+iA} ile
  çarpar, (c) sonuç φ = 2A−2π+δ (V5). 148'in docstring'i niyeti
  "bilinen faz ilerlemesi 2πτ ile döndür" diye yazar ve 150'nin KALEM'i
  Γ_rot = ⟨w e^{−iAdsΔ}⟩/⟨w⟩ (2A terimi YOK) der; uygulama bu tanımı
  gerçeklemiyor. **Ama "hangisinin doğru gözlenebilir olduğu" bu raporun
  kararı değildir** — her ikisinin tayfı da ölçülmüş olarak duruyor.
* **δ/(A·S) ≈ 2 açığının kaynakları ölçüldü, ama çarpanların
  ayrıştırılması bant başına yapıldı, global bir kanıt değil.** ×1.53
  (kesme) ve ×1.30 (σ) çarpımı 1.99; ölçülen oran 1.82–2.47. Kalan
  ±%15 bant-içi dağılımdır.
* **İlk iki t1 bandı (τ_eff ≈ 0.46, 0.50) oran özetlerini bozuyor** ve
  bu gizlenmedi: orada hem δ hem A·S sıfırdan geçiyor, oran tanımsıza
  yakın. Menzil sütunları o iki bandı içeriyor; medyan ve "τ ≥ 0.54"
  ifadeleri içermiyor. Her iki okuma da tabloda.
* **Sentetiklerin yüksek-τ bantları bozuk ve hükme girmedi.** A4 4 bant,
  P1 3 bant, keskin 1 bant elendi (§7 kuralı). Bozuk bantlar tablolarda
  ✗ ile duruyor.
* **σ_konv karşılaştırmaları yalnız gerçek gaz için anlamlıdır.**
  158, sentetiklerde b'nin ölçülemediğini yazmıştı; bu rapor onu
  değiştirmiyor — yalnız **işaretin** öngörülebilirliğini sınıyor.
* **Sentetik gaz üreticileri yeniden koşulmadı.** `z_keskin.npy` (152),
  `z_A4.npy` (154), `z_P1.npy` (155) önbellekten okundu; η zincirleri
  155'in önbelleğinden geldi. Momentlerin 158'in künyesiyle birebir
  aynı çıkması (σΔ² = 0.05466 gerçek, 0.04294 A4, 0.03499 keskin,
  0.14477 P1) bunu doğruluyor.
* **`dusuk`, `J14`, `J26`, `N5`, `N5z`, `P0` bu koşuda ölçülmedi.**
  Görev "en az bir sentetik gaz" istiyordu; üçü (keskin, A4, P1)
  koşuldu. Kalan dördünün S'si ölçülmedi ve hiçbir hükümde
  kullanılmadı.
* **Süreler ölçüt değildi.** 12 koşu iki akışta eşzamanlı; tipik koşu
  0.7–2.4 dk.

---

## Ek — dosyalar ve tekrar-üretim

| dosya | ne |
|---|---|
| `159_configs/159_cekirdek.py` | ölçüm çekirdeği: `Yerel` (kanal dizileri, pencere dilimleri, çizgi motoru), `olcS` (bant/çizgi döngüsü + mekanizma merdiveni) |
| `159_configs/159_kos.py` | `<veri> <taban> <t1\|g158>` — koşu sürücüsü |
| `159_configs/159_analiz.py` | T0–T4 bölümleri + fit konvansiyonu (`fit2`) + bant sağlığı (`saglik`) |
| `159_configs/159_dogrulama.py` | V1 (bit düzeyi), V2 (özdeşlik), V3 (sentetik kontrol), V5 (omurga) |
| `159_configs/159_figur.py` | 6 panelli figür |
| `159_phi_AS_sinavi.png` | figür |

Ham çıktılar `scratchpad/159/`: `S_<veri>_t<taban>_<ızgara>.json` (12),
`log_*.txt`, `analiz_cikti.txt`, `dogrulama_cikti.txt`, `kos_a.sh`,
`kos_b.sh`. η önbellekleri `scratchpad/155/eta_*.npz` (yeniden
kullanıldı, yenisi üretilmedi).

Tekrar üretmek için:

```
for v in son orta keskin A4 P1; do
  .venv/bin/python 159_configs/159_kos.py $v 0.40 g158
  .venv/bin/python 159_configs/159_kos.py $v 0.40 t1
done
.venv/bin/python 159_configs/159_kos.py son 0.52 g158
.venv/bin/python 159_configs/159_kos.py son 0.52 t1
.venv/bin/python 159_configs/159_kos.py A4  0.52 g158
.venv/bin/python 159_configs/159_dogrulama.py
.venv/bin/python 159_configs/159_analiz.py
.venv/bin/python 159_configs/159_figur.py
```
