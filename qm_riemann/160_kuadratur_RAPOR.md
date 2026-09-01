# 160 — mekanizma denkleminin SON ÇARPANLARI: kuadratür (σ) kanalı + sonlu-kesme

**Hedef.** 159, ölçülen fazın özdeşliğini kurdu — `φ_Γ = (4πτ_eff − 2π) + δ`,
`δ = arg⟨(ρ+iσ)·e^{+iA·X̃}⟩` (bağıl 7.8e−09) — ve δ ile KALEM'in `A·S`'si
arasında **kararlı bir ×2 açık** bıraktı; açığı iki ölçülmemiş parçaya
böldü: *sonlu-kesme* ×1.32–1.80 ve *kuadratür (σ) kanalı* ×1.27–1.38.
Bu rapor o iki parçayı **açılım yapmadan, tam ifadeyle** ölçtü ve δ'nın
muhasebesini kapattı.

**Figür:** `160_kuadratur.png` (6 panel).
P1 muhasebe merdiveni · P2 kapanış seviyesi (log) · P3 iki çarpan ·
P4 a'nın 32-fit ensemble'ı · P5 b'nin kapanışı, gaz gaz · P6 σ'nın kimliği.

---

> ## Kısa hüküm
>
> 1. **δ'NIN MUHASEBESİ KAPANDI — ve ±%10 değil, ±%0.07 ile.**
>    ```
>    δ = K1 + K2 ,  K1 = arg⟨ρ e^{iAX̃}⟩ ,  K2 = arg[(⟨ρ e^{iAX̃}⟩ + i⟨σ e^{iAX̃}⟩)/⟨ρ e^{iAX̃}⟩]
>    ```
>    Beş gazın **35 sağlıklı bandının 35'i** ±%10 içinde; en kötü bağıl hata
>    **%0.068**, mutlak artık **≤1.3e−04 rad**. Aynı bantlarda **A·S 2/35**,
>    **K1 (σ'suz) 2/35**. Kapanış çizgi düzeyinde **2.1e−13 (bağıl)**.
> 2. **K3 = 0 ÖZDEŞ OLARAK — ve bu bir sınav değil, muhasebenin KAPALI
>    olduğunun kanıtıdır.** Artık yalnız KALEM'in kendi değişkeni dsΔ
>    kullanıldığında doğuyor: **|δ − M0(dsΔ)| ≤ 3.3e−03 rad** (gerçek gazda
>    ≤8.1e−04, bağıl ≤%0.081). Yani **K3 ≡ değişken ikamesi sistematiğidir**
>    (dsΔ ↔ X̃), başka bir kanal değil.
> 3. **İki çarpan ölçüldü ve 159'un tahminini doğruladı.** Gerçek gaz,
>    taban 0.40, sıfır geçişinden uzak bantlar: **kesme K1/(A·S) = 1.524
>    [1.32, 1.58]**, **kuadratür δ/K1 = 1.299 [1.27, 1.38]**, çarpım
>    **1.940**. σ kanalı δ'nın **%23.0'ünü** taşıyor (taban 0.52'de %19.0).
>    Diğer gazlar: keskin 1.003 × 1.197 (σ payı %16), A4 1.176 × 1.736
>    (%42; bantlar %31–%98). **P1'de çarpanlar ÖLÇÜLEMEDİ** (δ ≤ 0.17,
>    |Γ| = 0.40–0.52; σ payı üç bantta −81/+75/+160 %).
> 4. **KUADRATÜR ÇARPANI, TOPLAM AÇIKTAN DAHA KONVANSİYON-DAYANIKLIDIR.**
>    Taban 0.40 ↔ 0.52 arasında δ/K1 bant başına **%0.4 – %6.3** değişiyor
>    (medyan farkı %3.8; 1.299 ↔ 1.235), δ/(A·S) ise **%1.7 – %14.2**
>    (medyan farkı %7.2; bant bant 2.465 ↔ 2.708). Yani σ kanalı gazın kendi
>    niteliğidir, taban konvansiyonunun artefaktı değil.
> 5. **σ'NIN KİMLİĞİ TAM OLARAK BULUNDU — ve üç adaydan biri kazandı.**
>    ```
>    (ρ_n + iσ_n)·e^{+iA(n+1)} = η_{n+1}·e^{−iA·C_n}·(e^{−iW m_0}·conj⟨c⟩)
>    C_n = Σ_{k≤n} X̃_k  (BİRİKMİŞ adım sapması)
>    ```
>    180 çizgide **bağıl ≤1.2e−08**. Yani σ, çizginin küresel taşıyıcıya
>    göre **biriktirdiği fazın** kuadratür izidir — *"çizginin yerel frekans
>    kayması"* adayı. **corr(σ,ρ): bant düzeyinde \|·\| ≤ 0.0010, çizgi düzeyinde ≤ 0.0064**
>    (tam kuadratür),
>    **corr(σ,Δρ) = −0.305 … +0.376** (τ ile tekdüze artıyor, tek büyük
>    korelasyon), **corr(σ,X̃) ve corr(σ,dX̃/dn) ≤ 0.026** ve oranları
>    ≈ **−1.6** (ardışık aralıkların anti-korelasyonu = itme).
>    **"Komşu-çizgi kaçağı" ÖLDÜ:** Cov(σ,X̃)/⟨ρ⟩'nun komşuluk boşluğuyla
>    korelasyonu altı sette **|r| ≤ 0.19**.
> 6. **σ, A'da İKİNCİ mertebede giriyor — reel-ağırlık + küçük-faz
>    ansatzının onu görememesinin iki ayrı nedeni bu.**
>    `Es = ⟨σe^{iAX̃}⟩/⟨ρ⟩ = s0 + iA·s1 − A²s2/2 − …` ile
>    **s0 = ⟨σ⟩/⟨ρ⟩ = −7.6e−05 … +2.0e−05** (sabit kayma YOK) ve
>    Re[Es]'i **−A²s2/2 = −A²⟨σX̃²⟩/(2⟨ρ⟩)** taşıyor. Birinci mertebe
>    (A·Cov(σ,X̃)/⟨ρ⟩) Es'in **sanal** kısmıdır ve K2'yi sürmez.
> 7. **KESMENİN "KAPALI FORMU" (kümülant serisi) ÇALIŞMIYOR — ölçüldü ve
>    kayda geçiyor.** 159'un 3. sıradaki adımıydı. `K1 = A·κ1 − A³κ3/6 + …`
>    iki terimle K1'i **%9.6 – %27.4** hatayla veriyor (A·S tek başına
>    %37); üçüncü terimin birinciye oranı **0.54 – 2.39**, yani seri
>    kullanışlı biçimde yakınsamıyor. **Kesmenin kapalı formu yok; TAM
>    karakteristik fonksiyonun kendisi kullanılmalıdır** (ki ucuzdur).
> 8. **a: HEDEF TUTTU — ama bir ÖZDEŞLİKLE, ve bu açıkça yazılıyor.**
>    158'in 32 fitinde (2 veri × 4 taban × 4 pencere) bu koşunun fit kodu
>    **a(φ) = 10.713 ± 0.059**'u birebir üretiyor. Aynı 32 fitte
>    **4π + dδ/dτ = 10.677 ± 0.079 → −0.61 σ**; ağırlık eşleştirilirse
>    **10.715 ± 0.058 → +0.04 σ**. Karşılaştırma: **4π + dK1/dτ = 11.080
>    ± 0.136 → +6.22 σ**, **4π + d(A·S)/dτ = 11.483 ± 0.075 → +13.04 σ**,
>    çıplak 4π → **+31.4 σ**. **σ kanalı 6.2σ'yı 0.6σ'ya indiriyor.**
> 9. **b: %50–80 AÇIK KAPANDI.** Aynı 32 fitte (ağırlık eşleşmiş)
>    b(φ) = −6.308 ± 1.296; **b(A·S) %56.6**, **b(K1) %85.7**,
>    **b(δ)=b(K1+K2) %99.9**. 158'in W-A referansı (8 fit) da birebir
>    üretildi: **b(φ) = −7.223 ± 0.791** (158: −7.2229 ± 0.7909).
> 10. **A4'ün b KARŞI-ÖRNEĞİ ÇÖZÜLDÜ — ve işareti SADECE σ kanalı veriyor.**
>    Ölçülen b(A4) = **+5.54**; b(A·S) = **−6.09 ✗**, b(K1) = **−4.25 ✗**,
>    b(δ) = **+3.12 ✓**. Yani A4'ün dışbükeyliğini reel-ağırlıklı mekanizma
>    **üretmiyor**, kuadratür kanalı **tek başına** üretiyor (K2'nin katkısı
>    +7.37). Beş gazın **beşinde** b(δ) doğru işareti veriyor (159: 4/5).
> 11. **Ölçüm zinciri yine bit düzeyinde denetlendi.** 160'ın φ, τ_eff,
>    σ_φ, |Γ|, δ, M1, M0, A·S'si **12 (gaz,taban,ızgara) üçlüsünde,
>    115 bantta 159'unkiyle 0.000e+00**.

---

## 1. Yöntem

`160_configs/160_cekirdek.py` hiçbir ölçüm parçasını kopyalamaz;
**159'un `Yerel` sınıfını miras alır** (`Yerel160(C159.Yerel)`), o da
155/156/154'ün zincirini import eder. Bant/çizgi döngüsü (`olc160`)
159'un `olcS`'iyle aynı konvansiyondadır: aynı aday listesi, aynı
220-örnekleme tohumu (21), aynı `gap < 2.5·dres` filtresi, aynı ara-nokta
(off-line) referansı `W' = w + gap/2`, aynı 8-grup round-robin jackknife.
**V1 bu iddiayı bit düzeyinde sınar ve 0.000e+00 verir.**

### 1a. Ne eklendi, ne çıkarıldı

| | |
|---|---|
| **ÇIKARILDI** | 159'un pencere makinesi (`PENCERELER` = 32/64/128; `|ĉ|²`, `⟨η²⟩_W`, pedestal). 159 §4a/4b o soruyu kapattı (pedestal 38–6371; bağlaşım bond ölçeğinde). Bu koşuda **koşulmadı** ve hiçbir sayıya girmedi. |
| **EKLENDİ** | (i) TAM değişken X̃ ile karakteristik fonksiyon `E1 = ⟨ρe^{iAX̃}⟩/⟨ρ⟩`, `Es = ⟨σe^{iAX̃}⟩/⟨ρ⟩`; (ii) ρ ve σ ölçülerinin merkezi momentleri μ₁…μ₄, s₀…s₃ ve kümülantları κ₁…κ₄; (iii) σ tanı korelasyonları; (iv) K1 ve K1+K2 için ayrı jackknife. |

Çıkarma/ekleme dengesi koşuları **hızlandırdı**: t1 ızgarasında (son, taban
0.40) 2.2 dk → **1.66 dk**, g158'de 0.9 dk → **0.34 dk**.

### 1b. e^{iAX̃} yeni trigonometri istemez (ve bu bir doğruluk kazancıdır)

`A·X̃_n = W(m_{n+1} − m_n) − A` olduğundan

    e^{iAX̃} = e^{−iA}·(cos W m_{n+1} + i sin W m_{n+1})(cos W m_n − i sin W m_n)

ve bu dört dizi zaten `zc`, `zp`, `ρ`, `σ` için hesaplanıyor. Böylece
kapanış özdeşliği **aynı yuvarlama hatalarını taşır** ve 1e−13'e iner
(V2). Kısayolun elle hesapla aynılığı sentetik kontrolde ayrıca sınandı:
maks fark **1.2e−10** (V6).

### 1c. Koşulan ölçümler

**18 koşu** (`S_<veri>_t<taban>_<ızgara>.json`), **199 ölçülen bant**:

| ızgara | koşular |
|---|---|
| **t1** = `izgara(0.44,0.80,0.04)` (9 bant, τ∈(0.44,0.80)) | son 0.40 · son 0.52 · keskin 0.40 · A4 0.40 · P1 0.40 |
| **g158** = `izgara(0.28,0.64,0.02)` (158'in ızgarası) | son 0.28/0.34/0.40/0.46/0.52 · orta 0.28/0.34/0.40/0.46 · keskin 0.40 · A4 0.40/0.52 · P1 0.40 |

g158'in dört tabanı × iki gerçek veri, **158'in 32-fit ensemble'ını**
yeniden kurmak içindir (§6). η önbellekleri `scratchpad/155/eta_*.npz`
dosyalarından yeniden kullanıldı; **yeni önbellek üretilmedi**.

### 1d. Tanım seçimleri ve gerekçeleri

| seçim | gerekçe |
|---|---|
| **Ağırlık ρ_n = Re[c_{n+1}conj⟨c⟩]** | 159 §"Dürüstlük"in ölçtüğü gibi: mekanizma denklemini tahmincinin cebrinde ÖZDEŞ yapan tek ağırlık; ⟨ρ⟩ = pow/4 (bu koşuda ölçülen `rho_norm` = 1.0000001 ± 0.00005). Görevin verdiği pencere-tabanlı P_loc 159 §4a'da ölçülüp elendi. |
| **Değişken X̃ (birincil)**, dsΔ (karşılaştırma) | Kimliği ÖZDEŞ yapan değişken X̃'dir; dsΔ ile koşulan sürüm (M1/M0) K3'ü **ölçmek** için tutuldu. İkisi 159 §4d'de korel 0.999877. |
| **Momentler X̃0 = X̃ − ⟨X̃⟩ etrafında** | Standart merkezi moment konvansiyonu. ⟨X̃⟩ = 1.296e−05 ⇒ A⟨X̃⟩ ≈ 5e−05 rad, δ'nın (medyan 0.129) **binde 0.4'ü**; hiçbir tabloya girmiyor ama `Xort` olarak kayıtlı. S = Cov(ρ,X)/⟨ρ⟩ merkezlemeye **duyarsızdır** (kovaryans). |
| **K2 = arg[(E1 + i·Es)/E1]** | Sarma-güvenli; K1 + K2 = arg(E1 + i·Es) tanım gereği. |
| **Bant sağlığı** | 159 §7'nin kuralı **birebir** (|Γ|≤1, |φ|≤2.8, \|τ_eff−τ̄\|≤yarım bant, ilk bozulmadan sonrakiler de düşer). Sağlıklı bant sayıları 159'unkilerle aynı: gerçek 9/9 ve 7/7, keskin 8/9, P1 6/9, A4 5/9. |

---

## 2. Muhasebe — tanım ve kapanış

159'un V2'de doğruladığı özdeşlikten (hiçbir yaklaşım yok):

    zp·conj(zc)·e^{−iA}/4 = ⟨(ρ_n + i·σ_n)·e^{+iA·X̃_n}⟩
    ρ+iσ ≡ c_{n+1}·conj⟨c⟩ ,  c_n = η_n e^{−iW m_n} ,  A = 2πW/L
    X̃_n = (m_{n+1}−m_n)·L/2π − 1

Bandın birleştirilmiş nesneleri (159'un `z1`/`zs` konvansiyonuyla aynı):

    z_E = Σ(p_on·E1_on − p_off·E1_off)/Σ(p_on − p_off)
    z_S = Σ(p_on·Es_on − p_off·Es_off)/Σ(p_on − p_off)

    (K1)   K1 = arg z_E             reel ağırlık, TAM karakteristik fonksiyon
    (K2)   K2 = arg[(z_E + i z_S)/z_E]        kuadratür düzeltmesi
    (K3)   K3 = δ_meas − (K1+K2)              artık

**K1 içinde küçük-faz açılımı YOKTUR**; A·S onun birinci mertebesidir, yani
159'un "sonlu-kesme" çarpanı **K1/(A·S)** olarak kendiliğinden ölçülür.

### 2a. Kapanış ölçüldü (V2–V4, V6)

| | sınav | sonuç |
|---|---|---|
| **V2** | çizgi düzeyi: `E1 + i·Es = [zp conj(zc)e^{−iA}/4]/⟨ρ⟩` | 18 koşu, maks **bağıl 2.1e−13** |
| **V3** | bant düzeyi: `δ − (K1+K2)` | sağlıklı bantlarda **maks 1.3e−04 rad**; t1 ızgarasında bağıl **≤ %0.068** |
| **V4** | `δ − M0(dsΔ)` — **değişken ikamesi** | maks **3.3e−03 rad**; gerçek gazda ≤8.1e−04 (bağıl ≤%0.081), keskin ≤%1.2, P1 ≤%4.9 |
| **V6** | sentetik kontrol (bilinen çizgi, bilinen adım) | `|δ − (K1+K2)| = 2.7e−13`; kısayol↔elle maks 1.2e−10 |

> **K3 hakkında hüküm: K3 ≡ 0'dır (özdeşlik), ve bu bir başarı değil bir
> KAPANIŞ KANITIDIR.** Ölçülebilir tek artık, KALEM'in kendi değişkeni
> dsΔ'nın X̃ yerine kullanılmasından doğar ve büyüklüğü yukarıdadır.
> Yani "δ'da ρ ve σ dışında bir kanal var mı?" sorusunun yanıtı
> **hayır**dır: 1e−13'te yok.

---

## 3. Bant bant muhasebe (τ ∈ (0.44, 0.80), t1 ızgarası)

### 3a. Gerçek gaz, taban 0.40 (9/9 sağlıklı)

| τ_eff | A·σ_X | \|Γ\| | δ_meas | A·S | K1 | K2 | K1+K2 | δ/(A·S) | δ/K1 | **K1/(A·S)** | K1 hatası % | K1+K2 hatası % |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.4607 | 0.609 | 0.902 | −0.0335 | +0.0137 | −0.0120 | −0.0215 | −0.0335 | −2.450 | +2.789 | −0.878 | +64.14 | −0.0010 |
| 0.4997 | 0.669 | 0.884 | −0.0975 | −0.0237 | −0.0632 | −0.0343 | −0.0975 | +4.118 | +1.543 | +2.668 | +35.21 | −0.0010 |
| 0.5392 | 0.731 | 0.860 | −0.1803 | −0.0731 | −0.1313 | −0.0489 | −0.1803 | +2.465 | +1.372 | **+1.796** | +27.14 | +0.0004 |
| 0.5791 | 0.800 | 0.834 | −0.2811 | −0.1350 | −0.2171 | −0.0640 | −0.2811 | +2.082 | +1.295 | **+1.608** | +22.78 | +0.0012 |
| 0.6188 | 0.875 | 0.789 | −0.4179 | −0.2064 | −0.3264 | −0.0915 | −0.4179 | +2.025 | +1.280 | **+1.581** | +21.90 | −0.0003 |
| 0.6578 | 0.948 | 0.740 | −0.5909 | −0.3046 | −0.4643 | −0.1266 | −0.5909 | +1.940 | +1.273 | **+1.524** | +21.42 | +0.0001 |
| 0.6980 | 1.032 | 0.684 | −0.8422 | −0.4212 | −0.6444 | −0.1978 | −0.8422 | +1.999 | +1.307 | **+1.530** | +23.49 | −0.0009 |
| 0.7382 | 1.063 | 0.669 | −1.1282 | −0.6058 | −0.8685 | −0.2597 | −1.1282 | +1.862 | +1.299 | **+1.434** | +23.02 | −0.0000 |
| 0.7766 | 1.072 | 0.699 | −1.6189 | −0.8874 | −1.1748 | −0.4441 | −1.6189 | +1.824 | +1.378 | **+1.324** | +27.43 | −0.0006 |

**Kapanış özeti (9/9 sağlıklı bant, δ'ya göre bağıl hata):**

| seviye | ortalama | medyan | menzil | **±%10 içinde** |
|---|---|---|---|---|
| A·S (KALEM, 1. mertebe) | +63.2 % | +50.6 % | [+45.2, +140.8] | **0/9** |
| K1 (reel ağırlık, TAM kar. f.) | +29.6 % | +23.5 % | [+21.4, +64.1] | **0/9** |
| **K1+K2 (kuadratür dahil)** | −0.000 % | −0.000 % | [−0.001, +0.001] | **9/9** |

### 3b. Beş gaz, iki taban — tek tabloda

| gaz / taban | sağlıklı | **kesme K1/(A·S)** medyan [menzil] | **kuadratür δ/K1** medyan [menzil] | çarpım δ/(A·S) | **K2/δ payı** | maks bağıl \|K3(X̃)\| | maks bağıl \|δ−M0\| |
|---|---|---|---|---|---|---|---|
| **gerçek-son t0.40 (t1)** | 9/9 | **+1.524** [1.32, 1.58] | **+1.299** [1.27, 1.38] | +1.940 | **+23.0 %** | 1.2e−05 | 8.1e−04 |
| **gerçek-son t0.52 (t1)** | 7/7 | **+1.605** [1.49, 1.71] | **+1.235** [1.22, 1.32] | +1.966 | **+19.0 %** | 1.2e−05 | 9.6e−04 |
| **keskin t0.40** | 8/9 | +1.003 [0.94, 1.26] | +1.197 [1.05, 1.45] | +1.249 | +16.4 % | 6.3e−05 | 1.2e−02 |
| **A4 t0.40** | 5/9 | +1.176 [1.09, 1.67] | **+1.736** [1.45, 2.95] | +2.046 | **+41.9 %** | 1.9e−05 | 1.1e−03 |
| **P1 (GUE) t0.40** ‡ | 6/9 | +0.189 [−0.40, 0.89] | +0.553 [−1.66, 4.00] | +0.666 | +75.0 % ‡ | 6.8e−04 | 4.9e−02 |
| **gerçek-son t0.40 (g158)** | 12/12 | +1.648 [1.56, 2.32] | +1.310 [1.26, 1.49] | +2.160 | +23.7 % | 2.5e−05 | 3.6e−03 |
| **gerçek-orta t0.40 (g158)** | 12/12 | +1.683 [1.57, 2.11] | +1.317 [1.27, 1.45] | +2.217 | +24.1 % | 1.6e−05 | 4.0e−03 |

(Çarpan özetleri **δ'nın sıfır geçişinden uzak** bantlarda: |δ| bandın
maksimumunun %20'si üstü. Sıfır geçişi civarında oran tanımsıza yakındır
ve tam tablolar `scratchpad/160/analiz_cikti.txt`'de duruyor.
**‡ P1'in çarpanları ve σ payı bir ÖLÇÜM DEĞİLDİR** — üç bandın değerleri
−81 %, +75 %, +160 %; medyan yalnız kayıt için yazıldı ve hiçbir hükme
girmedi.)

**Beş gazın 35 sağlıklı bandında toplu:** A·S ±%10 içinde **2/35**
(maks hata %898, P1), K1 **2/35** (maks %500, P1), **K1+K2 35/35**
(maks **%0.068**).

### 3c. Gazlar arasındaki desen — yeni kayıt

* **Gerçek gaz:** iki çarpan da 1'in üstünde ve **kararlı**; σ payı τ boyunca
  %21–%27 (taban 0.40) ve %18–%24 (taban 0.52). Yani σ kanalı gerçek gazda
  τ'dan neredeyse bağımsız bir **sabit pay**dır.
* **keskin:** kesme çarpanı **1.00** — 159'un "keskin'de mekanizma denklemi
  neredeyse tutuyor" gözleminin nedeni ölçüldü: A·σ_X = 0.51–0.66 (bütün
  gazların en küçüğü) olduğu için yüksek mertebeler orada gerçekten küçük.
  Kalan açık **tümüyle** kuadratürdür (1.197).
* **A4:** kesme küçük (1.18) ama **kuadratür büyük (1.74)** ve τ ile hızla
  büyüyor (τ_eff = 0.61'de δ/K1 = 57). A4'te δ'nın %31–%98'i σ kanalındadır.
  **b'nin işaretini de bu taşıyor (§6d).**
* **P1 (GUE):** her iki çarpan da **1'in altında** ve işaret değiştiriyor.
  σ payı üç ölçülebilir bantta **−81 %, +75 %, +160 %** — yani P1'de pay
  **ölçülemiyor** (δ'nın kendisi ≤ 0.17, |Γ| = 0.40–0.52; en sönümlü gaz).
  Kayda geçen tek şey, P1'in reel-ağırlık okumasının **tamamen çöktüğü**
  gaz olduğudur (K1'in işareti δ'nınkiyle altı sağlıklı bandın üçünde
  ters) — ve bu, 159'un "P1'de A·S fazı fazla öngörüyor" kaydının nedenidir.

### 3d. Taban dayanıklılığı — kuadratür çarpanı daha sağlam

Gerçek gaz, t1'in 7 örtüşen bandında:

| τ̄ | δ/K1 (t0.40) | δ/K1 (t0.52) | fark % | δ/(A·S) (t0.40) | δ/(A·S) (t0.52) | fark % |
|---|---|---|---|---|---|---|
| 0.54 | +1.372 | +1.378 | **+0.4** | +2.465 | +2.708 | +9.9 |
| 0.58 | +1.295 | +1.318 | **+1.8** | +2.082 | +2.377 | +14.2 |
| 0.62 | +1.280 | +1.266 | **−1.1** | +2.025 | +2.171 | +7.2 |
| 0.66 | +1.273 | +1.225 | **−3.8** | +1.940 | +2.018 | +4.0 |
| 0.70 | +1.307 | +1.225 | **−6.3** | +1.999 | +1.966 | −1.7 |
| 0.74 | +1.299 | +1.235 | **−4.9** | +1.862 | +1.925 | +3.4 |
| 0.78 | +1.378 | +1.320 | **−4.2** | +1.824 | +1.965 | +7.7 |

Kuadratür çarpanının taban duyarlılığı **|·| ≤ %6.3 (medyan %3.8)**,
toplam açığınki **≤ %14.2 (medyan %7.2)**. Yani **σ kanalı taban
konvansiyonunun artefaktı değildir.**

---

## 4. Sonlu-kesmenin anatomisi — ve kapalı-form denemesinin ÖLÜMÜ

159'un 3. sıradaki adımı şuydu: *"Kapalı-form adayı artık S değil,
⟨ρe^{iAX}⟩/⟨ρ⟩'nun kümülant açılımıdır."* Bu koşu o adayı kurdu ve
**reddetti**.

ρ-ölçüsünün kümülantlarıyla (κ₁ = μ₁, κ₂ = μ₂−μ₁², κ₃ = μ₃−3μ₁μ₂+2μ₁³, …):

    K1 = A·κ1 − A³·κ3/6 + O(A⁵)        (A·S = A·κ1 tam olarak)
    log|E1| = −A²κ2/2 + A⁴κ4/24 + O(A⁶)   (SÖNÜM)

### 4a. Gerçek gaz, taban 0.40, t1

| τ_eff | A | A·κ1 (=A·S) | −A³κ3/6 | toplam | **K1 (ölçülen)** | artık | artık/K1 | \|A³κ3/6\| / \|A·κ1\| |
|---|---|---|---|---|---|---|---|---|
| 0.4607 | 2.894 | +0.01361 | −0.02762 | −0.01401 | −0.01200 | +0.00201 | **16.8 %** | 2.03 |
| 0.4997 | 3.140 | −0.02373 | −0.05673 | −0.08046 | −0.06318 | +0.01728 | **27.4 %** | 2.39 |
| 0.5392 | 3.388 | −0.07321 | −0.07075 | −0.14396 | −0.13135 | +0.01261 | **9.6 %** | 0.97 |
| 0.5791 | 3.639 | −0.13517 | −0.12864 | −0.26381 | −0.21707 | +0.04673 | 21.5 % | 0.95 |
| 0.6188 | 3.888 | −0.20657 | −0.16890 | −0.37547 | −0.32640 | +0.04907 | 15.0 % | 0.82 |
| 0.6578 | 4.133 | −0.30483 | −0.22745 | −0.53228 | −0.46429 | +0.06799 | 14.6 % | 0.75 |
| 0.6980 | 4.386 | −0.42144 | −0.31561 | −0.73704 | −0.64437 | +0.09267 | 14.4 % | 0.75 |
| 0.7382 | 4.639 | −0.60601 | −0.41815 | −1.02416 | −0.86849 | +0.15567 | 17.9 % | 0.69 |
| 0.7766 | 4.879 | −0.88786 | −0.47773 | −1.36558 | −1.17482 | +0.19076 | 16.2 % | 0.54 |

> **HÜKÜM: kesmenin kullanışlı bir kapalı formu YOK.** İki terimli kümülant
> serisi K1'i **%9.6–%27.4** hatayla veriyor (A·S tek başına %37) ve
> **her bantta AŞIYOR** — yani hata "kesme" değil "aşırı-düzeltme"dir.
> Üçüncü terimin birinciye oranı **0.54 – 2.39**: seri, açılım parametresi
> A·σ_X = 0.61 … 1.07 olduğu bölgede kullanışlı biçimde yakınsamıyor.
> Doğru nesne **tam karakteristik fonksiyonun kendisidir** (E1) ve onu
> ölçmek A·S'yi ölçmekten pahalı değildir.

### 4b. Sönüm — bir terim yetiyor (düşük τ'da), sonra yetmiyor

| τ_eff | log\|E1\| | −A²κ2/2 | +A⁴κ4/24 | iki terim |
|---|---|---|---|---|
| 0.4607 | −0.1798 | **−0.1855** | +0.0108 | −0.1747 |
| 0.5392 | −0.2469 | **−0.2675** | +0.1028 | −0.1647 |
| 0.6188 | −0.3238 | **−0.3831** | +0.1544 | −0.2287 |
| 0.6980 | −0.4024 | **−0.5326** | +0.3340 | −0.1985 |
| 0.7766 | −0.3377 | −0.5744 | +0.6881 | +0.1137 |

Tek terim (−A²κ2/2) τ ≤ 0.55'te **%3–%8** ile tutuyor, τ = 0.78'de %70
sapıyor; κ₄ terimi eklendiğinde **daha da kötüleşiyor**. Aynı hüküm:
**sönümün de kapalı formu yok**, |E1| doğrudan ölçülmeli. (κ₂ ölçüldü: 0.04429 → 0.05538 → 0.04825, ham σΔ² = 0.05466'ya göre
**−%19.0 … +%1.3**; §HÜKÜM (iii)'e bakınız.)

---

## 5. σ'NIN KİMLİĞİ

### 5a. TAM kimlik (kalem) — ve sayısal doğrulaması

`u_n ≡ ρ_n + iσ_n = c_{n+1}conj⟨c⟩ = η_{n+1}·|⟨c⟩|·e^{iθ_n}` ve
`θ_n = −(W m_{n+1} + arg⟨c⟩)`. `m_{n+1} = m_0 + (2π/L)(n+1) + (2π/L)C_n`
olduğundan (`C_n = Σ_{k≤n} X̃_k`) **hiçbir yaklaşım yapmadan**:

> **(D1)   (ρ_n + iσ_n)·e^{+iA(n+1)} = η_{n+1}·e^{−iA·C_n}·(e^{−iW m_0}·conj⟨c⟩)**

Yani **(ρ,σ) bir fazördür**: bond başına −A dönen taşıyıcının üstüne
**birikmiş adım sapmasının fazı −A·C_n** biniyor; ρ o fazörün küresel faza
izdüşümü, **σ ise kuadratürü**. Ölçüldü: üç gaz × iki bant × 30 çizgi =
**180 çizgide maks bağıl fark 1.2e−08**.

Bunun iki doğrudan sonucu ölçüldü:

| ölçüm | değer | okuma |
|---|---|---|
| **corr(σ, ρ)** | çizgi düzeyi ort **−0.0004**, maks \|·\| **0.0064** (180 çizgi); bant düzeyi maks \|·\| **0.00102** | **σ, ρ'ya TAM DİK** — yeni bir alan değil, aynı fazörün öteki bileşeni |
| **⟨σ⟩/⟨ρ⟩ = s0** | −7.6e−05 … +2.0e−05 | **K2'de sabit kayma YOK**; K2 bütünüyle X̃ ile bağlaşımdan geliyor |

### 5b. Üç adayın sınavı

| aday | sınav | sonuç |
|---|---|---|
| **(b) çizginin yerel frekans kayması** | D1 kimliği | **KAZANDI, 1.2e−08 ile.** σ = biriken faz kaymasının kuadratür izi |
| **(a) yerel faz-gradyanı** | corr(σ, dX̃/dn) | **ZAYIF ama SİSTEMATİK.** Bant düzeyinde −0.0263 … +0.0179; corr(σ,X̃) ile **oranı ≈ −1.6** (aşağıda) |
| **(c) komşu-çizgi kaçağı** | Cov(σ,X̃)/⟨ρ⟩'nun komşuluk boşluğuna bağlılığı | **ÖLDÜ.** Altı (gaz,bant) setinde korel(·, gap/dres) = **−0.13 … +0.19**; hiçbirinde anlamlı değil |

### 5c. Teşhis korelasyonları (bant düzeyi, güç-ağırlıklı)

| gaz | τ_eff | **corr(σ,ρ)** | **corr(σ,Δρ)** | corr(σ,X̃) | corr(σ,dX̃/dn) | corr(ρ,X̃) | oran dX̃/X̃ |
|---|---|---|---|---|---|---|---|
| gerçek-son | 0.4607 | −0.0003 | **−0.3046** | −0.01154 | +0.01793 | +0.00219 | −1.55 |
| gerçek-son | 0.5392 | +0.0003 | +0.0086 | −0.00594 | +0.00923 | −0.00554 | −1.55 |
| gerçek-son | 0.6188 | −0.0003 | **+0.2689** | −0.00150 | +0.00233 | −0.00750 | −1.55 |
| gerçek-son | 0.6980 | −0.0000 | **+0.3717** | +0.00115 | −0.00184 | −0.00621 | −1.60 |
| gerçek-son | 0.7766 | −0.0003 | **+0.3590** | +0.00179 | −0.00272 | −0.00448 | −1.52 |
| keskin | 0.6187 | +0.0001 | +0.3189 | +0.00115 | −0.00154 | +0.00019 | −1.34 |
| A4 | 0.5386 | −0.0003 | +0.2560 | +0.00960 | −0.01559 | +0.00489 | −1.62 |
| P1 (GUE) | 0.6169 | +0.0001 | +0.2321 | −0.00154 | +0.00227 | +0.00029 | −1.47 |

**Üç okuma:**

1. **corr(σ,ρ) = 0'dır** — kuadratür, tanım gereği değil **ölçülerek**.
2. **corr(σ,Δρ) tek büyük korelasyondur** ve τ ile tekdüze artıyor
   (−0.30 → +0.38). Ama **1 değildir**: "σ = ρ'nun bond-gradyanının
   kuadratürü" okuması, σ'nın varyansının en fazla **%24'ünü** açıklıyor
   (§5d). Nedeni ölçüldü: genlik η bond ölçeğinde yavaş değişmiyor.
3. **corr(σ, dX̃/dn) ≈ −1.5 × corr(σ, X̃)**, altı gaz-bandında da.
   σ_n yalnız X̃_n'yi görür (X̃_{n+1}'i değil), dolayısıyla saf halde oran
   −1.00 olurdu; ölçülen −1.34 … −1.62, yani **corr(σ_n, X̃_{n+1}) ≈
   −0.5·corr(σ_n, X̃_n)**: ardışık aralıklar **anti-korelasyonludur**
   (seviye itmesi). Bu, σ'nın adım-belleğinin bir bond ötesini de gördüğünün
   doğrudan izidir.

### 5d. İki REDDEDİLEN kalem-önerisi (kayda geçiyor)

**(i) "σ ≈ tan(A/2)·ρ + Δρ/sin A"** (yavaş-genlik yaklaşımı).
σ'yı ρ ve Δρ'ya regres ettim: ölçülen α = −0.18 … +0.45, öngörü tan(A/2)
= +8.8 … −0.84; ölçülen β = −0.42 … +0.50, öngörü 1/sinA = +4.5 … −1.02;
**R² = 0.001 … 0.239**. Sentetik kontrolde de aynı (R² = 0.017). **Model
yanlıştır**: η_{n+1} ile η_{n+2} bond ölçeğinde bağımsıza yakın olduğundan
Δρ'yı taşıyıcı dönmesi değil **genlik sıçraması** yönetiyor.

**(ii) "Cov(σ,X̃)/⟨ρ⟩ ≈ −A·σΔ²"** (son bondun fazörü döndürmesi).
Fazörü son adımdan arındırıp (u⁰ = u·e^{+iAX̃}) ayrıştırdım:
`σ = σ⁰cos(AX̃) − ρ⁰sin(AX̃)`. Öngörü, `Cov(σ⁰cos, X̃) ≈ 0` varsayımına
dayanıyordu. **Ölçüm: öyle değil.** İki terim karşılaştırılabilir ve
**ters işaretli**; "eski kuadratür" payı |Cov(σ⁰cos,X̃)|/(toplam) =
**0.37 – 0.81**. Sonuç: kapalı formun ölçülene oranı altı sette
**+0.054 / −0.250 / −0.295 / +0.336 / −0.077 / +0.193** — **kararlı bir
değer yok**, ve kapalı form çıkarıldığında çizgiden çizgiye saçılım
azalmıyor (sd 0.0071 → 0.0071). **Reddedildi.**

### 5e. σ kanalı A'da İKİNCİ mertebeden girer — ve K2'yi taşıyan moment ⟨σX̃²⟩'dir

`Es = s0 + iA·s1 − A²s2/2 − iA³s3/6 + …`, `s_k = ⟨σ X̃0^k⟩/⟨ρ⟩`.
`Re[Es]`in ilk X-bağımlı terimi **−A²s2/2**; `A·s1 = A·Cov(σ,X̃)/⟨ρ⟩` ise
**sanal** kısımdadır. Gerçek gaz, taban 0.40:

| τ_eff | s0 | A·s1 | **−A²s2/2** | −A³s3/6 | Re[Es] | Im[Es] | **K2 (ölçülen)** | Re[Es·conj E1]/\|E1\|² | 2. mert. seri | seri/ölçülen |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.4607 | −5.8e−06 | −0.0759 | **−0.0273** | +0.0017 | −0.0204 | −0.0730 | −0.0215 | −0.0234 | −0.0316 | ×1.470 |
| 0.5392 | −7.6e−06 | −0.0797 | **−0.0683** | +0.0003 | −0.0529 | −0.0774 | −0.0489 | −0.0542 | −0.0734 | ×1.501 |
| 0.6188 | +1.8e−05 | −0.0457 | **−0.1148** | +0.0079 | −0.0899 | −0.0396 | −0.0915 | −0.1002 | −0.1335 | ×1.458 |
| 0.6980 | +7.3e−06 | +0.0807 | **−0.1401** | +0.0018 | −0.1079 | +0.0799 | −0.1979 | −0.2008 | −0.2416 | ×1.221 |
| 0.7766 | −7.6e−05 | +0.3614 | −0.0382 | −0.0683 | −0.0291 | +0.3083 | −0.4441 | −0.4144 | −0.3998 | ×0.900 |

* **Doğrusallaştırma iyi:** K2 ≈ Re[Es·conj E1]/|E1|² **%7–%9** içinde.
* **Moment kesmesi kötü:** ikinci-mertebe seri K2'yi ×1.47 (düşük τ) →
  ×0.90 (yüksek τ) veriyor. Yani **σ kanalının da kapalı formu yok**,
  aynı §4'teki sebeple.
* **A·s1 τ ≈ 0.66'da işaret değiştiriyor** (−0.076 → +0.361) ama K2 hiç
  değiştirmiyor — Es'in *sanal* kısmının K2'yi sürmediğinin doğrudan izi.

> **σ'nın FİZİKSEL YORUMU (ölçülene dayanarak).** σ, çizginin küresel
> taşıyıcıya göre **biriktirdiği fazın** (A·C_n) kuadratür bileşenidir:
> yerel frekans kayması. K2 ise o kaymanın **adımın karesiyle** olan
> bağlaşımıdır (−A²⟨σX̃²⟩/2⟨ρ⟩). Sezgisel okuma: bir bandın çizgileri,
> uzun aralık dizilerinde fazı **geri**, kısa aralık dizilerinde **ileri**
> kalır; ρ bu iki durumda simetrik olarak sönümlenirken (bu, |Γ|'dır) σ
> **antisimetrik** bir kalıntı bırakır ve bandın ortalama fazını kaydırır.
> Reel-ağırlık ansatzında bu kanalın yeri yoktur çünkü ansatz fazörü
> tek bileşene indirger; küçük-faz açılımında da yoktur çünkü katkı
> A²'dendir.

---

## 6. a ve b'nin HASSAS TÜRETİMİ

Fit konvansiyonu 158/157'ninki, **birebir**: apsis τ_eff, ağırlık 1/σ_jk,
`polyfit`, τ₀ = φ'nin köklerinden x ortalamasına en yakını, **a = dy/dτ|τ₀**,
**b = ½·d²y/dτ²|τ₀** (kuadratikte c₀). Dört pencere = 158'inkiler
(W-A τ̄∈[0.43,0.61] kuadratik · W-D τ̄∈[0.43,0.59] · W-B sıfır-merkezli
τ_eff−τ₀ ∈ [−0.075,+0.105] · W-C W-A kübik).

### 6a. Fit kodunun denetimi — 158'i BİREBİR üretiyor

| ensemble | 158'in yazdığı | **bu koşu** |
|---|---|---|
| 32 fit (son+orta × 4 taban × 4 pencere), a | 10.713 ± 0.059 | **10.713 ± 0.059** |
| 32 fit, b | (son 16 fit −6.211 / orta 16 fit −6.404 ⇒ −6.31) | **−6.308 ± 1.296** |
| 8 fit (yalnız W-A), a | 10.717 ± 0.067 | **10.717 ± 0.067** |
| 8 fit (yalnız W-A), b | −7.2229 ± 0.7909 | **−7.223 ± 0.791** |

### 6b. a — 32 fit ensemble (gerçek gaz)

| nicelik | ortalama | sd | menzil | fark(10.713) | **kaç σ (0.059)** |
|---|---|---|---|---|---|
| **a(φ) — ÖLÇÜLEN** | **10.713** | 0.059 | [10.595, 10.824] | −0.000 | −0.01 |
| **4π + dδ/dτ = 4π + d(K1+K2)/dτ** | **10.677** | 0.079 | [10.577, 10.810] | −0.036 | **−0.61** |
| 4π + dK1/dτ (σ YOK) | 11.080 | 0.136 | [10.952, 11.315] | +0.367 | **+6.22** |
| 4π + dM1/dτ (159'un öngörüsü) | 11.081 | 0.136 | [10.952, 11.316] | +0.368 | +6.24 |
| 4π + d(A·S)/dτ (1. mertebe) | 11.483 | 0.075 | [11.374, 11.615] | +0.770 | **+13.04** |
| çıplak 4π | 12.566 | — | — | +1.853 | **+31.4** |

**Ağırlık eşleştirilirse** (bütün eğriler σ_φ ile tartılır; fit-ağırlığı
sistematiğini yalıtır): 4π + dδ/dτ = **10.715 ± 0.058 → +0.04 σ**;
4π + dK1/dτ = 11.099 ± 0.114 → +6.54 σ; 4π + d(A·S)/dτ = 11.500 ± 0.069
→ +13.33 σ.

> **HEDEF TUTTU (≤1.5σ) — ama DÜRÜSTÇE: bu bir ÖZDEŞLİKTİR.**
> Kinematik omurga τ'da tam doğrusal olduğundan `a = 4π + dδ/dτ` cebirsel
> bir sonuçtur; ölçülen 0.036 (kendi ağırlıkları) / 0.002 (eşleşmiş
> ağırlık) fark **fit ağırlığı ve bant-düzeyi arındırma artığıdır**.
> Dolayısıyla bu satır a'nın *öngörüsü* değil, **muhasebenin ve omurga
> ayrıştırmasının tam olduğunun sınavıdır** — ve geçti.
>
> **Öngörücü içerik merdivenin ALT basamaklarındadır:** A·S 13.0σ, K1 6.2σ.
> **σ kanalı 6.2σ'yı 0.6σ'ya indiriyor** ve bu, serbest parametresi olmayan,
> ayrı ölçülmüş bir nesnenin (Es) yaptığı bir kapanıştır.

### 6c. a — gaz gaz (taban 0.40, W-A; hata % ölçülene göre)

| gaz | **a(φ)** | 4π+dδ/dτ | hata % | 4π+dK1/dτ | hata % | 4π+d(A·S)/dτ | hata % | çıplak 4π | hata % |
|---|---|---|---|---|---|---|---|---|---|
| gerçek-son | 10.615 | **10.577** | **−0.4** | 10.954 | +3.2 | 11.405 | +7.4 | 12.566 | +18.4 |
| gerçek-orta | 10.656 | **10.585** | **−0.7** | 10.955 | +2.8 | 11.406 | +7.0 | 12.566 | +17.9 |
| keskin | 11.660 | **11.685** | **+0.2** | 11.999 | +2.9 | 12.083 | +3.6 | 12.566 | +7.8 |
| A4 | 13.054 | **13.151** | **+0.7** | 12.583 | −3.6 | 12.453 | −4.6 | 12.566 | −3.7 |
| P1 (GUE) | 11.988 | **12.007** | **+0.2** | 11.991 | +0.0 | 13.170 | +9.9 | 12.566 | +4.8 |
| **korel(a_meas, ·)** | — | **+0.9999** | | +0.9769 | | +0.7463 | | — | |
| **gaz sıralaması** | — | **AYNI** | | FARKLI (keskin↔P1) | | FARKLI | | — | |

**Konvansiyon bütçesi (taban 0.40, 4 pencerenin sd'si):**

| gaz | a(φ) | 4π+dδ/dτ | 4π+dK1/dτ | 4π+d(A·S)/dτ |
|---|---|---|---|---|
| gerçek-son | 10.652 ± 0.056 | 10.610 ± 0.035 | 10.955 ± 0.004 | 11.392 ± 0.013 |
| gerçek-orta | 10.668 ± 0.038 | 10.601 ± 0.020 | 10.976 ± 0.025 | 11.422 ± 0.017 |
| keskin | 11.667 ± 0.016 | 11.679 ± 0.018 | 12.012 ± 0.023 | 12.092 ± 0.016 |
| A4 | 13.032 ± 0.071 | 13.086 ± 0.088 | 12.630 ± 0.049 | 12.534 ± 0.089 |
| P1 (GUE) | 11.771 ± 0.286 | 11.796 ± 0.268 | 11.886 ± 0.120 | 13.204 ± 0.044 |

### 6d. b — %50–80 açık KAPANDI, ve A4 karşı-örneği ÇÖZÜLDÜ

**32-fit ensemble, ağırlık EŞLEŞMİŞ** (b(φ) ≡ b(δ) kimliğini yalıtır):

| nicelik | ortalama | sd | \|ort\|/sd | negatif | **kapanış (\|ort\|/\|b(φ)\|)** |
|---|---|---|---|---|---|
| **b(φ) — ÖLÇÜLEN** | −6.308 | 1.296 | 4.87 | 32/32 | 100.0 % |
| **b(δ) = b(K1+K2)** | −6.299 | 1.303 | 4.83 | 32/32 | **99.9 %** |
| b(K1) — σ YOK | −5.408 | 0.918 | 5.89 | 32/32 | **85.7 %** |
| b(A·S) — 1. mertebe | −3.570 | 0.866 | 4.12 | 32/32 | **56.6 %** |

(Kendi jackknife ağırlıklarıyla: b(δ) %94.2, b(K1) %82.6, b(A·S) %54.5.
158'in W-A referansı −7.22'ye göre: b(δ) %82.3, b(K1) %72.2, b(A·S) %47.6.)

> **159'un "büyüklük tutmuyor (%50–%80 eksik)" kaydı kapandı.** Kesme
> çarpanı %56.6 → %85.7'yi getiriyor, kuadratür kanalı kalan %14.2'yi.

**b'nin işareti, gaz gaz (taban 0.40, W-A):**

| gaz | **b(φ)** | b(A·S) | kapanış % | işaret | b(K1) | kapanış % | işaret | **b(δ)** | kapanış % | işaret |
|---|---|---|---|---|---|---|---|---|---|---|
| **gerçek-son** | −7.76 | −3.83 | 49.4 | ✓ | −5.65 | 72.8 | ✓ | **−6.50** | 83.8 | ✓ |
| gerçek-orta | −7.80 | −3.54 | 45.4 | ✓ | −5.45 | 69.8 | ✓ | **−6.00** | 77.0 | ✓ |
| keskin | −3.81 | −3.34 | 87.7 | ✓ | −3.70 | 97.1 | ✓ | **−4.16** | 109.2 | ✓ |
| **A4** | **+5.54** | **−6.09** | 110.0 | **✗** | **−4.25** | 76.8 | **✗** | **+3.12** | 56.4 | **✓** |
| P1 (GUE) | +6.83 | +1.42 | 20.7 | ✓ | +3.06 | 44.7 | ✓ | **+7.00** | 102.5 | ✓ |

> **A4 KARŞI-ÖRNEĞİNİN TEŞHİSİ.** A4'ün dışbükeyliğini **reel-ağırlıklı
> mekanizma üretmiyor** (b(A·S) = −6.09, b(K1) = −4.25, ikisi de içbükey);
> **kuadratür kanalı tek başına üretiyor** — K2'nin eğriliğe katkısı
> **+7.37**, ve işareti o çeviriyor. Bu, §3b'nin ölçümüyle tutarlıdır:
> A4, δ'sının **%42'sini** (bir bantta %98'ini) σ kanalından alan tek
> sağlıklı gazdır.
>
> **DÜRÜSTLÜK NİTELEMESİ:** `b(φ) ≡ b(δ)` bir özdeşlik olduğundan
> "b(δ) doğru işareti veriyor" **bağımsız bir öngörü değildir**. Kayda
> geçen şey **anatomidir**: A4'ün b'sinin işaretini taşıyan fiziksel
> nesne `⟨σ X̃²⟩`'dir, `Cov(ρ,X̃)` değil. Bu, φ'den bağımsız olarak
> ölçülebilen bir nesnedir; onu gazın korelasyon yapısından türetmek
> **hâlâ açık** bir hedeftir (§Sıradaki adım).

---

## 7. Denetim

| | sınav | sonuç |
|---|---|---|
| **V1** | 160'ın φ, τ_eff, σ_φ, \|Γ\|, δ, M1, M0, A·S'si = 159'unki mi? | **12 (gaz,taban,ızgara) üçlüsü, 115 ortak bant: maks fark 0.000e+00** (her sütunda; bit düzeyinde özdeş) |
| **V1b** | Bu raporun fit kodu 158'in ensemble'larını üretiyor mu? | 32 fit: a = **10.713 ± 0.059**, b = −6.308 (158: 10.713 ± 0.059, −6.31). 8 fit W-A: a = **10.717 ± 0.067**, b = **−7.223 ± 0.791** (158: 10.717 ± 0.067, −7.2229 ± 0.7909) |
| **V2** | Çizgi düzeyi kapanış `E1 + i·Es = [zp conj(zc)e^{−iA}/4]/⟨ρ⟩` | 18 koşu, maks **bağıl 2.1e−13** |
| **V3** | Bant düzeyi `δ − (K1+K2)` | sağlıklı bantlarda maks **1.3e−04 rad** (bütün bantlarda 7.1e−03, bozuk P1 bantları dahil) |
| **V4** | `δ − M0(dsΔ)` (değişken ikamesi) | sağlıklı bantlarda maks **3.3e−03 rad** |
| **V5** | Kümülant yeniden kurulumu `K1 − (Aκ1 − A³κ3/6)` | **0.011 – 16.0 rad** — bir kapanış DEĞİL, kesmenin ölçülen artığı (§4) |
| **V6** | Sentetik kontrol: bilinen tek çizgi + bilinen adım dizisi (N = 200k, L = 12, τ = 0.62, σ_X = 0.22) | `|δ − (K1+K2)| = 2.7e−13`; kısayol↔elle maks 1.2e−10; ölçülen çarpanlar K1/(A·S) = 0.771, δ/K1 = 1.613 |
| **D1** | σ'nın fazör kimliği `(ρ+iσ)e^{iA(n+1)} = η_{n+1}e^{−iA C_n}K` | 3 gaz × 2 bant × 30 çizgi = **180 çizgide maks bağıl 1.2e−08** |

---

# HÜKÜM

## (i) δ'nın muhasebesi KAPALIDIR — ve kapanma seviyesi ±%0.07

> **δ = arg⟨(ρ+iσ)e^{iAX̃}⟩ = K1 + K2**, K1 = arg⟨ρe^{iAX̃}⟩ (reel ağırlık,
> **tam** karakteristik fonksiyon), K2 = kuadratür düzeltmesi.
> Çizgi düzeyinde 2.1e−13; bant düzeyinde sağlıklı 35 bandın 35'i ±%10
> içinde, en kötü **%0.068**. Aynı bantlarda A·S **2/35**, K1 **2/35**.

159'un iki açık çarpanı ayrı ayrı ölçüldü ve **çarpımları ölçülen açığı
veriyor**: gerçek gazda kesme **1.524**, kuadratür **1.299**, çarpım
**1.940** (159 aynı bantlarda 1.53 × 1.30 tahmin etmişti). Üçüncü bir
kanal **yoktur** (K3 ≡ 0, 1e−13). KALEM'in değişkeni dsΔ kullanılırsa
doğan artık ölçüldü ve küçüktür (gerçek gazda bağıl ≤%0.08).

## (ii) σ, YENİ BİR ALAN DEĞİL — fazörün kuadratürü ve kimliği tamdır

> **(ρ_n + iσ_n)·e^{+iA(n+1)} = η_{n+1}·e^{−iA·C_n}·K**, C_n = Σ_{k≤n} X̃_k
> (180 çizgide 1.2e−08).

σ, çizginin küresel taşıyıcıya göre **biriktirdiği fazın** kuadratür
izidir — yani *yerel frekans kayması*. Üç adaydan **komşu-çizgi kaçağı
öldü** (gap ile korelasyon |r| ≤ 0.19), **faz-gradyanı zayıf ama
sistematik** (corr(σ,dX̃/dn) ≈ −1.5 × corr(σ,X̃), ardışık aralıkların
itmesinin izi), **frekans kayması kazandı**. corr(σ,ρ) = 0.0000 ± 0.0006:
tam diklik. ⟨σ⟩/⟨ρ⟩ ≤ 8e−05: K2'de sabit kayma yok.

K2'yi taşıyan moment **⟨σX̃²⟩**'dir (Re[Es] ≈ −A²⟨σX̃²⟩/2⟨ρ⟩): **σ kanalı
A'da ikinci mertebedendir.** KALEM'in ansatzının onu görememesinin iki
ayrı nedeni budur — reel ağırlık (fazörü tek bileşene indirger) **ve**
birinci-mertebe kesme.

## (iii) İKİ KAPALI-FORM ADAYI DA ÖLDÜ, ve bu 159'un 3. adımını kapatıyor

159 "kapalı-form adayı ⟨ρe^{iAX}⟩'nun kümülant açılımıdır" demişti.
Ölçüldü: iki terimli seri K1'i **%9.6–%27.4** hatayla ve **her bantta
aşırı-düzelterek** veriyor; üçüncü/birinci terim oranı **0.54–2.39**.
Aynı hüküm sönüm için (−A²κ2/2 τ ≥ 0.70'te %30–%70 sapıyor) ve σ serisi
için (2. mertebe K2'yi ×0.90 – ×1.50). **Kesmenin kullanışlı bir kapalı
formu yoktur; tam karakteristik fonksiyon ölçülmelidir** — ve maliyeti
A·S ile aynıdır (bu koşu 159'dan HIZLI koştu).

Yan kayıt: ρ-ölçüsü altındaki adım varyansı κ₂ = 0.0443 … 0.0554, ham
σΔ² = 0.05466'ya göre **−%19.0 … +%1.3** (dokuz bandın sekizinde KÜÇÜK,
tekdüze olarak τ ile büyüyor) — çizginin ağırlığı düşük-sapmalı bondları
kayırıyor, ve bu kayırma yüksek τ'da kayboluyor. "Yerel güç ile yerel
adım" bağlaşımının doğrudan sayısı budur; ilk kez ölçüldü.

## (iv) a ve b: açık kapandı, ama kapatan şey bir ÖZDEŞLİK

**a.** 158'in 32 fitinde a(φ) = 10.713 ± 0.059 birebir üretildi.
Merdiven: çıplak 4π **+31.4σ** → A·S **+13.0σ** → K1 (σ'suz) **+6.2σ**
→ K1+K2 **−0.61σ** (ağırlık eşleşmiş: **+0.04σ**). Görevin hedefi
(≤1.5σ) tuttu — ama `a = 4π + dδ/dτ` cebirsel bir özdeşliktir, bir
öngörü değil. **Öngörücü içerik alt basamaklardadır** ve orada σ kanalı
6.2σ'yı 0.6σ'ya indiriyor; bu, serbest parametresiz, ayrı ölçülmüş bir
nesnenin yaptığı kapanıştır.

**b.** %50–80 açık kapandı: b(A·S) **%56.6** → b(K1) **%85.7** →
b(K1+K2) **%99.9** (32 fit, ağırlık eşleşmiş). 158'in W-A referansı
(−7.223 ± 0.791) da birebir üretildi.

**A4 karşı-örneği çözüldü.** Ölçülen b = +5.54; reel-ağırlıklı mekanizma
**içbükey** öngörüyor (b(A·S) = −6.09, b(K1) = −4.25), kuadratür kanalı
**dışbükeyliği tek başına üretiyor** (K2'nin katkısı +7.37). Beş gazın
beşinde işaret doğru (159: 4/5). **Nitelendirme:** b(φ) ≡ b(δ) olduğundan
bu bir öngörü değil, **anatomidir** — A4'ün b'sinin işaretini taşıyan
nesnenin `⟨σX̃²⟩` olduğu, `Cov(ρ,X̃)` olmadığı ölçülmüştür.

## Sıradaki adım (bu ölçümün işaret ettiği)

1. **Artık türetilecek nesne `a − 4π` değil, `⟨ρ e^{iAX̃}⟩` ve `⟨σ e^{iAX̃}⟩`
   çiftidir.** İkisi de φ'den bağımsız ölçülüyor, ikisi de gazın adım
   dizisinin fonksiyonu, ve ikisi birlikte δ'yı 1e−13'te veriyor.
   Kapalı form arayışı kümülantlarda değil, **C_n'in (birikmiş adım
   sapması) dağılımında** olmalı: D1 kimliği gösteriyor ki her iki nesne
   de `⟨η_{n+1}e^{−iA C_n}·(bir bond gecikmeli çekirdek)⟩` biçimindedir.
2. **⟨σX̃²⟩ ölçülmüş bir gaz ayırıcısıdır ve sınanmadı.** Gerçek gazda
   K2/δ payı %19–%27 (iki tabanda kararlı), keskin'de %16, A4'te %42
   (P1'de ölçülemiyor). Bu pay tek başına gazları sıralıyor mu, 158'in dokuz gazında
   sınanabilir (ucuz: bu koşunun çekirdeği hazır).
3. **P1'in muhasebesi sağlıklı bantlarda bile kırılgan.** τ ≥ 0.66'da
   bant sağlığı düşüyor, |Γ| > 1 çıkıyor ve K1 δ'nın 5 katına kadar
   sapıyor. GUE-boyalı gazda tahmincinin neden bu kadar erken bozulduğu
   ölçülmedi.
4. **De-rotasyon çarpanının işareti hâlâ kalemle karara bağlanmadı**
   (159'un 1. maddesi). Bu rapor onu değiştirmiyor: ölçülen ayrışımdır,
   φ'nin de δ'nın da tayfı kayıtlıdır.
5. **`dusuk`, `J14`, `J26`, `N5`, `N5z`, `P0` bu koşuda da ölçülmedi.**
   a/b ensemble'ı yalnız gerçek gaz (son+orta) için tam; sentetiklerde
   taban ekseni yalnız 0.40 (+A4'te 0.52).

---

## Dürüstlük notları

* **Uydurma yok; her sayı yeniden koşuldu.** Bütün tablo değerleri
  `160_configs/160_kos.py`'nin ürettiği 18 JSON'dan,
  `160_configs/160_analiz.py`'nin çıktısından (`scratchpad/160/analiz_cikti.txt`),
  `160_dogrulama.py`'ninkinden (`dogrulama_cikti.txt`) ve
  `160_sigma.py`'ninkinden (`sigma_cikti.txt`) otomatik alındı.
  **Tablodan okunan tek şey** 158'in referanslarıdır (a = 10.713 ± 0.059,
  10.717 ± 0.067; b = −7.2229 ± 0.7909) ve her biri kaynağıyla işaretli —
  üstelik hepsi bu koşuda **yeniden hesaplandı ve aynı çıktı** (V1b).
* **"Muhasebe kapandı" ifadesi bir başarı iddiası değil, bir ÖZDEŞLİK
  kaydıdır.** K1+K2 = δ cebirsel olarak doğrudur; ölçülen 1e−13, kodun
  o cebri gerçeklediğinin kanıtıdır. Görevin sorduğu "±%10 kapanıyor mu"
  sorusunun anlamlı yanıtı **K1 seviyesindedir** (kapanmıyor: 2/35) ve
  **dsΔ ikamesi seviyesindedir** (kapanıyor: gerçek gazda ≤%0.08).
* **a'nın 1.5σ hedefi bir özdeşlikle tutmuştur ve bu §6b'de açıkça
  yazılmıştır.** Aynı şey b(δ)'nın A4 işareti için geçerlidir (§6d).
  Bağımsız öngörü değeri taşıyan tek satırlar A·S ve K1 satırlarıdır.
* **İki kalem-önerisi ölçülüp REDDEDİLDİ ve saklanmadı:** (i) "σ ≈
  tan(A/2)ρ + Δρ/sinA" (R² ≤ 0.24, katsayılar öngörüden mertebelerce
  uzak; sentetik kontrolde de düşüyor); (ii) "Cov(σ,X̃)/⟨ρ⟩ ≈ −A·σΔ²"
  (oran altı sette −0.30 … +0.34, kararsız; kapalı form çıkarıldığında
  saçılım azalmıyor).
* **Kümülant kapalı-formu 159'un önerisiydi ve bu koşu onu ÇÜRÜTTÜ.**
  Rapor bunu bir başarısızlık olarak değil, "aday kapatıldı" olarak
  kaydediyor; V5 satırı 0.011–16.0 rad ile tabloda duruyor.
* **σ tanı korelasyonlarının küçüklüğü bir zayıflık değil, beklenendir.**
  σ'nın varyansını taşıyıcı yönetiyor; adım bağlaşımı varyansın binde
  biridir ama ortalamada **hayatta kalan** tek parçadır. Buna karşılık
  corr(σ,Δρ) = −0.30 … +0.38 gerçek ve büyüktür; "σ = ρ'nun gradyanının
  kuadratürü" okumasının ne kadar tuttuğunun ölçüsüdür (en fazla %24).
* **Bozuk bantlar hiçbir özete girmedi.** A4 4 bant, P1 3 bant, keskin
  1 bant elendi (159 §7'nin kuralı, birebir). Tablolarda ✗ ile duruyorlar.
  §3b'nin çarpan özetleri ayrıca **sıfır geçişinden uzak** bantlarla
  sınırlıdır (|δ| > 0.2·maks); tam tablolar çıktıda.
* **Bağıl hatalar sıfır geçişi civarında anlamsızdır** ve bu gizlenmiyor:
  bütün koşularda en büyük bağıl |δ−M0| = 1.59 P1'in |δ| = 1.3e−04 olan
  bandındandır. Bu yüzden §2a **mutlak** artıkları da veriyor.
* **P1'in "sağlıklı" bantları bile marjinal.** |Γ| = 0.40–0.52, muhasebe
  kapanışı diğer gazlardan 1–2 mertebe kötü (6.8e−04 bağıl). P1'in
  sayıları hükümlere yalnız "desen" düzeyinde girdi.
* **159'un pencere makinesi bu koşuda KOŞULMADI.** `|ĉ|²`, `⟨η²⟩_W` ve
  pedestal ölçümleri 159 §4a/4b'de kapatılmıştı; burada yeniden
  ölçülmedikleri için bu raporda hiçbir sayıları yok.
* **Sentetik gaz üreticileri yeniden koşulmadı.** `z_keskin.npy` (152),
  `z_A4.npy` (154), `z_P1.npy` (155) önbellekten; η zincirleri
  `scratchpad/155/eta_*.npz`'den. Künyelerin 158/159'la birebir aynı
  çıkması (σΔ² = 0.05466 gerçek, 0.04294 A4, 0.03499 keskin, 0.14477 P1)
  bunu doğruluyor.
* **Gaz sıralaması ölçütü fit ağırlığına duyarlıdır.** 4π+dM1/dτ, 159'un
  ağırlıklarıyla (σ_δ) sıralamayı koruyordu; M1'in kendi jackknife'ıyla
  keskin ↔ P1 takas oluyor. 4π+dδ/dτ **her iki ağırlıkla da** sıralamayı
  koruyor. Tabloda kendi ağırlıkları kullanıldı.
* **Süreler ölçüt değildi.** 18 koşu 4 akışta eşzamanlı; koşu süreleri
  **0.28 – 1.68 dk** (159'un 0.7–2.4 dk'sından hızlı, pencere makinesi
  çıkarıldığı için).

---

## Ek — dosyalar ve tekrar-üretim

| dosya | ne |
|---|---|
| `160_configs/160_cekirdek.py` | `Yerel160` (159'un `Yerel`'ini miras alır) + `olc160`: K1/K2 merdiveni, kümülantlar, σ tanıları |
| `160_configs/160_kos.py` | `<veri> <taban> <t1\|g158>` — koşu sürücüsü |
| `160_configs/160_analiz.py` | B1–B7: muhasebe, kümülantlar, σ, σ serisi, a/b ensemble, taban, özet |
| `160_configs/160_dogrulama.py` | V1 (bit düzeyi), V2–V5 (kapanışlar), V6 (sentetik kontrol) |
| `160_configs/160_sigma.py` | σ'nın derin tanısı: D1 kimliği, ayrışım, kaçak sınavı |
| `160_configs/160_figur.py` | 6 panelli figür |
| `160_kuadratur.png` | figür |

Ham çıktılar `scratchpad/160/`: `S_<veri>_t<taban>_<ızgara>.json` (18),
`log_*.txt`, `analiz_cikti.txt`, `dogrulama_cikti.txt`, `sigma_cikti.txt`,
`sigma_tani.json`, `kos_a..d.sh`. η önbellekleri `scratchpad/155/eta_*.npz`
(yeniden kullanıldı, yenisi üretilmedi).

Tekrar üretmek için:

```
for v in son orta; do for t in 0.28 0.34 0.40 0.46; do
  .venv/bin/python 160_configs/160_kos.py $v $t g158
done; done
.venv/bin/python 160_configs/160_kos.py son 0.52 g158
for v in keskin A4 P1; do
  .venv/bin/python 160_configs/160_kos.py $v 0.40 g158
  .venv/bin/python 160_configs/160_kos.py $v 0.40 t1
done
.venv/bin/python 160_configs/160_kos.py A4  0.52 g158
.venv/bin/python 160_configs/160_kos.py son 0.40 t1
.venv/bin/python 160_configs/160_kos.py son 0.52 t1
.venv/bin/python 160_configs/160_dogrulama.py
.venv/bin/python 160_configs/160_sigma.py
.venv/bin/python 160_configs/160_analiz.py
.venv/bin/python 160_configs/160_figur.py
```
