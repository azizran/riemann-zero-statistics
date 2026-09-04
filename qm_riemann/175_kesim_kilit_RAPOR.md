# 175 — FAZLANIN AYRIŞIMI: KESİM PAYI + KİLİT PAYI
### (4 Eylül 2026, Opus tayfası — KALEM_KESIM_VE_KILIT_04EYL2026.md)

**Soru.** 174, gerçek gazın `ΔM = +%4.97`'lik fazlasını λ ailesinin
**dışına** yerleştirdi (`R_η = 1.2883`, ailenin tavanı 1.2745) ve açığın
**%83'ünü τ > 0.70**'e adresledi — tam da 152'nin ölçtüğü erfc kuyruğunun
bölgesine. 175'in görevi bu fazlayı ikiye ayırmak: ne kadarı **KESİM
ŞEKLİ** (gerçek erfc-kuyruklu, ikiz keskin), ne kadarı **SAF FAZ-KİLİDİ**?

Betikler (`175_configs/`) — `175g` dışında hepsi ön-mühürlü (`175g`
ölçümden sonra yazıldı ve öyle etiketlidir); git'e dokunulmadı:

| betik | kapı | koşu | ham çıktı |
|---|---|---|---|
| `175a_onkayit.py` | **K2 ön-kaydı** + K1 inşa kapıları | 25 s | `175/ONKAYIT_K2.json` |
| `175b_K1_kesim.py` | **K1** — 162 makinesi, kesim ailesi (4 gaz) | 2.1-2.2 dk/gaz | `174/K1_<gaz>.json` |
| `175c_K3_kesim.py` | **K3** — DC kaçağı defteri, kesim ailesi (4 gaz) | 1.1-1.3 dk/gaz | `174/K3_<gaz>.json` |
| `175d_K2_yuzlesme.py` | **K2** — dondurulmuş kuralın uygulanışı + ayrışım | 2 s | `175/K2.json` |
| `175e_K3_carpimsal.py` | **K3** — çarpımsal κ (asal ↔ kule) + rakip sınav (**ön-mühürlü**) | 1.2 dk | `175/K3M.json` |
| `175f_K4_figur.py` | **K4** — θ defteri + figür | 25 s | `175/K4.json`, `175_kesim_kilit.png` |
| `175g_kappa_yasasi.py` | κ'nın düzeltilmiş kapalı yasası (**ön-mühürsüz ikinci tur**) | 1.5 dk | `175/K3G.json` |

Önbellek kökü `/private/tmp/claude-501/.../scratchpad/`; girdi
`155/eta_*.npz`, `164,155/z_HA4.npy`, `165/tayf_*.npz`, `167/C_*.json`,
`172/G1-G4.json`, `174/K1_*,K3_*.json`; çıktı `175/`.

---

> ## TEK CÜMLELİK HÜKÜM
>
> **174'ün "aranacak yeni eksen" borcu ödendi: EKSEN KESİMDİR** —
> λ ailesinin taşımadığı dört nicelik (`R_η` 1.2883, `m3` 0.2768,
> `skew x1` 0.2817, `skew ds` 0.4681) **kesim ailesinin içindedir**
> (`R_η` tavanı 1.2745 → **1.2963**), ve iki eksen **tamamlayıcıdır**:
> **X kanalı bir λ-gazı, η kanalı bir KESİM-gazıdır**.
> Ama KALEM'in literal mührü **ÖLDÜ** (ön-kayıt bunu önceden yazdı):
> τ>0.70 açığı erfc-ikizde küçülmedi, **%56 büyüdü ve işaret değiştirdi**
> — `τ_c = 0.68` fazla serttir; gerçek gaz iki ikizin **arasındadır**.
> O "ara" ölçüldü ve **açığın %83'ünü taşıyan bantlarda SABİTTİR**:
>
> > **KESİM PAYI = %39,  KİLİT PAYI = %61**  (τ > 0.70; `f_b` = 0.3975
> > ve 0.3871, %2.7 içinde aynı) — ve aynı %39 `g_E`'den (0.390) ve
> > `μ̂²_E`'den (0.397) bağımsız olarak çıkıyor.
>
> **θ ise kesim ekseninin dışındadır ve TERS yöndedir** (`f(θ) = −0.143`;
> gerçek 0.9219, kesim ailesinin tavanı 0.8933): ΔM'nin %51.5'lik θ payı
> kesim şekliyle kapanmıyor, **büyüyor** ⇒ **H-K4 ÖLDÜ**.
> Ve K3'te iki ölüm bir yasa doğurdu: aşımın çarpımsal olduğu hipotezi
> **öldü** (fazla asal çizgilerde de var, orada çarpımsal kanal özdeş
> sıfırdır), ön-kayıtlı **rakip kazandı** — 174e'nin `W_pos` çarpanı
> **hatalıydı**; doğru yasa `κ = −πA_Qτ_Q cos(πτ_Q)` sönümsüzdür ve
> 7000+ asal çizgide **±%10**'da tutar. Kulelerin kalan **×2.0**'ı ise
> gerçekten **çarpımsaldır**: `A^eff = A_Q + (πτ_Q/2)Σ_{q₁q₂=Q}A₁A₂`
> sıfır parametreyle oranı **1.07**'ye indiriyor.

---

## 0. ÖN-KAYIT (K2) — 4 Eylül 2026, 15:43:46 +03

`175a_onkayit.py` hiçbir 175 ölçümü koşmadan **önce** koştu; kendi
sha256'sını
(`13a5d1957f5983bf69add1d889241a9d34dbb1356eaa5178c97518018617cc8c`)
`175/ONKAYIT_K2.json`'a yazdı ve dosya bir daha yazılmaz (varlık
kontrolü).

### 0.1 DÜRÜSTLÜK BEYANI (ön-kayıtta aynen duruyor)

1. **Erfc-ikiz YENİDEN İNŞA EDİLMEDİ.** KALEM "`167_insa.main` ile erfc
   zarflı gaz kur" diyor; ama o gaz **zaten vardır**: `HA4`
   (`164_insa`, `tip = "sadakatli"`, `par = {tau_c: 0.68, delta: 0.125,
   c: −0.5}`). 167_insa aynı çözücüyü (`164_insa.coz_sadakatli`) çağırır
   ve **deterministiktir** — yeniden inşa bit-bit aynı diziyi verirdi.
   10 dakikalık artıksız bir tekrar yerine mevcut gazın kapıları
   **bağımsız olarak** `z_HA4.npy`'den yeniden denetlendi (§0.3).
2. **σ_ds ön-mühürü GEÇERSİZDİR ve öyle kaydedildi.** HA4'ün
   `σ_ds² = 0.18861`'i 164'te **yayımlanmıştır**. Ölçülmüş bir sayıya
   "geniş bant ön-mühür" atmak sahte olurdu. KALEM'in bu maddesi bu
   koşuda **uygulanamaz**; uzatmanın dört ıskası defterde kalır,
   beşincisi eklenmez.
3. **172'nin defteri OKUNMUŞTUR.** Ön-kayıt anında `172/G1.json` ve
   `G2.json`'un `HA4 / K090 / K070 / E060` satırları (`g_E`, `g_X`, `θ`,
   `KAL`, `μ̂²_E`, `π_E`, `ρ_E`, `σ`'lar) biliniyordu — **körlük iddiası
   yoktur**. Ön-kayıt edilen şey KURALDIR; öngörüler bu defterden
   TÜRETİLİR ve sınanan şey, defterin (172: site `s_n`, τ≤0.95, 8981
   çizgi) **başka iki makinede** — 162'nin girişim ölçümünde (site
   `m_n`, τ≤0.86, 3425 çizgi) ve 174d'nin DC-kaçağı defterinde — tutup
   tutmadığıdır.
4. **Ön-kayıt anında ÖLÇÜLMEMİŞ olanlar** (hepsi bu koşuda ilk kez
   ölçüldü): kesim ailesinin hiçbir üyesinde `R_η` / `R_Ĉ` / `R_X`,
   bant-bant defteri, üçüncü momentler; hiçbirinde `Σξ_q` DC-kaçağı
   defteri; **hiçbir gazda** κ'nın çizgi tipine (asal ↔ kule) göre
   ayrışması.

### 0.2 ERFC ZARFININ TAM FORMÜLÜ (KALEM'in istediği)

152'nin "S3'ü oturtan kesim" dediği (`152` raporu, A4 satırı) ve 164'ün
`HA4` gazını kuran zarf, **kaynak koddan birebir doğrulanarak**:

```
w_q   = ½ · erfc( (τ_q − τ_c) / Δ ) ,      τ_c = 0.68 ,  Δ = 0.125
τ_q   = ω_q / L_hedef ,  ω_q = log q ,  q = p^k  (asal kuvvet)
a_q   = 1 / (π · k · √q)
A_q   = λ · a_q · w_q ,                    λ = 1.00
S(t)  = − Σ_q A_q sin(ω_q t)
seviye:  N̄(z) + S(z) = n − ½               (c = −½ konvansiyonu, 164 §3d)
L_hedef = 12.030361949155486 ,  nline = 15450  (τ ≤ 1.00)
```

`175a` bu formülü `164_insa.merdiven`'in ürettiği ağırlık dizisiyle
karşılaştırdı: **maks |w_kod − w_formül| = 0.00e+00** (`175a`).
Zarfın değerleri: `w(0.50) = 0.9791`, `w(0.60) = 0.8173`,
`w(0.68) = 0.5000`, `w(0.75) = 0.2142`, `w(0.80) = 0.0873`,
`w(0.90) = 0.0064`.

### 0.3 K1 — ERFC-İKİZİN İNŞA KAPILARI (bağımsız denetim, `175a`)

| kapı | eşik | ölçülen | hüküm |
|---|---|---|---|
| ilk-kök hücre | 300000 / 300000 | **300000 / 300000** | ✓ |
| tekne sayısı | 300000 | **300000** | ✓ |
| maks \|F\| | ≤ 1e−8 | **1.863e−09** (aşan tekne 0) | ✓ |
| sıralılık | TAM | **TAM**, min Δz = 0.149853 | ✓ |
| σ_ds (z'den) | — | **0.43429** (σ_ds² = 0.18861; 164 yayımı 0.18861, fark 0.0e+00) | — |
| L (z'den) | — | **12.029593242** (164: 12.029593242) | — |
| `155/z_HA4` ≡ `164/z_HA4` | — | **bit-bit aynı** | ✓ |
| **SAĞLIK** `R_bant` (169 filtresi, lo ∈ [0.52,0.68]) | ≥ 0.98 | **1.4527 1.5522 1.7204 2.0149 2.5853** ⇒ min 1.4527 | ✓ |

> **K1 KAPILARI GEÇTİ.** Erfc-ikiz sadakatlidir ve ölçüme uygundur.

### 0.4 KESİM EKSENİ (inşa tarafı, ölçüm değil)

174'ün λ ailesinin yanına bu koşuda bir **ikinci eksen** kondu: aynı
sadakatli çözücüyle kurulmuş **beş kesim gazı** (hepsi λ = 1, c = −½).
τ̄_A := Σ A²τ / Σ A² (güç-ağırlıklı ortalama τ; yalnız merdivenden):

| gaz | kesim | τ̄_A | rms S | P(τ>0.5)/P |
|---|---|---|---|---|
| `Hkeskin` | keskin τ ≤ 1.00 | 0.31712 | 0.38253 | 0.23830 |
| `K090` | keskin τ ≤ 0.90 | 0.29323 | 0.37550 | 0.20950 |
| `K070` | keskin τ ≤ 0.70 | 0.24347 | 0.35819 | 0.13125 |
| **`HA4`** | **erfc 0.68/0.125** | **0.22605** | **0.35025** | **0.09259** |
| `E060` | erfc 0.60/0.125 | 0.20493 | 0.34025 | 0.04752 |

### 0.5 DONDURULMUŞ KURAL (ön-kayıtta aynen)

| # | kural | ölçüt / ölüm eşiği |
|---|---|---|
| **P1** | `R_η = π_E × r̄`, `r̄ = 0.98709 ± 0.00300` (174'ün 8 gazından) | dört kesim gazında `\|R_η/π_E − 1\| ≤ %2`; biri aşarsa köprü kesim ekseninde KIRIK |
| **P2** | gerçek gaz **kesim ailesinin içinde mi**? Öngörülen aralık **[1.2177, 1.2904]** ∋ 1.28829 | kapsarsa ⇒ *fazlanın ekseni KESİMDİR*; kapsamazsa H-K1'in η-kanalı ölür |
| **P3** | `\|ort E\| = √(μ̂²_E · K_E/π_E)` özdeşliği; öngörüler: `son` +0.046900, `Hkeskin` +0.062139, `K090` +0.064031, `K070` +0.036543, **`HA4` +0.031399**, `E060` +0.016825 | 174d'nin `Σ_q ξ_q`'su ‰5 içinde, işaret POZİTİF |
| **P4** | **MÜHÜR (KALEM, literal):** τ>0.70 açığı ≥ %50 küçülürse H-K1 yaşar. Keskin açık = **−1.270040e−02**, eşik **6.3502e−03** | ön-kayıtlı BEKLENTİ: **ölecek** — erfc-ikiz alttan ıskalayacak, işaret dönecek (sandviç) |
| **P5** | **kesim kesri** `f(y) := [y(son) − y(Hkeskin)] / [y(HA4) − y(Hkeskin)]` — koordinat yok, ara değer yok, parametre yok | bütün `f`'ler aynıysa açık SAF KESİMDİR (bantlarda ±0.15); yayılım = **KİLİT PAYI**; `f < 0` ⇒ kesim o nicelikte TERS yönde |
| **P6** | `M ≡ g_E·g_X²·θ` özdeşliği; çapa = keskin ikiz. `Δlog M = +0.05754 = +0.02455 + 0.00148 + 0.03152` (kalıntı 4.6e−16) | ayrışım: kesim = `f*·ΔlogM(erfc)`, kilit = E/X artığı, θ = θ artığı; toplam ÖZDEŞ |
| **P7** | çarpımsal κ: merdiven asal kuvvetler ⇒ `q₁q₂ = Q` yalnız **KULE**'de çözülür | **K3-a:** medyan(asal) ≥ 1.5 (τ>0.60) ⇒ H-K3 ÖLÜR. **K3-c (rakip):** τ>0.55'te medyan `\|κ\|/\|κ̂₀\|` ∈ [0.8,1.25] ⇒ rakip kazanır |
| **P8** | θ defteri | `f(θ) < 0` ⇒ H-K4 ÖLÜR |

**Ön-kayıtta yazılı 172-defteri türevi kesim kesirleri** (okunmuş
sayılardan; sınavı 162/174d makineleri yapacak):

```
f(π_E) = +0.8533   f(ort E) = +0.4957   f(g_E) = +0.3828   f(g_X) = +0.0414
f(θ)   = −0.1430   f(M)     = −0.5566   f(ρ_E) = +0.2051   f(σ_ds) = −7.5124
```

### 0.6 K3'ÜN TÜRETİMİ (ön-kayıtta, ölçümden önce)

`s_n = m̄_n + δ(m̄_n)`, `δ(m̄) = ḡ Σ_c A_c cos(πτ_c) sin(ω_c m̄)`
(143'ün G-yasası köşesi). `κ(ν) = ⟨e^{iνm̄} e^{iνδ}⟩` açılımında
`⟨e^{iνm̄} cos(Ωm̄)⟩ = ½·1{Ω = ±ν}` olduğundan

```
κ⁽¹⁾ = −π A_ν τ_ν cos(πτ_ν)                       (174e'nin yasası)
κ⁽²⁾ = (π²τ_ν²/2)·( S_ν − D_ν )
  S_ν = Σ_{q₁q₂=Q} A₁A₂ cos(πτ₁)cos(πτ₂)          (ÇARPIMSAL kanal)
  D_ν = 2 Σ_{q₁/q₂=Q} A₁A₂ cos(πτ₁)cos(πτ₂)       (fark kanalı)
```

**Merdiven asal kuvvetlerdir** (`q = p^k`) ⇒ `q₁q₂ = Q`'nun merdiven içi
çözümü **yalnız `Q = p^k`, `k ≥ 2`** için vardır: **kule çizgileri**.
`Q = p²` için kapalı, sıfır-parametreli oran:

```
|κ⁽²⁾/κ⁽¹⁾| = 2 τ_p cos²(πτ_p) / ( |cos 2πτ_p| · W_pos(2τ_p) )
```

| τ_p | 0.25 | 0.30 | 0.35 | 0.40 | 0.45 |
|---|---|---|---|---|---|
| τ_Q | 0.50 | 0.60 | 0.70 | 0.80 | 0.90 |
| oran | ∞ (cos 2πτ_p = 0) | 1.139 | 0.505 | 0.242 | **0.076** |

> Yani çarpımsal katkı **τ_Q ≈ 0.5 civarında baskın**, ama τ_Q → 0.9'da
> yalnız **%8**. 174e'nin ölçtüğü aşım ise τ ile **büyüyor** (×1.9 →
> ×4.2). Ön-kayıt bu gerilimi ölçümden önce yazdı ve bir **rakip**
> koydu (K3-c): fazla, `W_pos`'un aşırı sönümü olabilir.

---

## 1. K1 — GİRİŞİM ORANLARI, KESİM AİLESİNDE (`175b`, 162 makinesi)

174b'nin makinesi (= 162'nin `Taban162`'si: taban 0.40, cap 4000,
τ_çizgi 0.86, `NITER = 0`, çıkarım `c_q(x) = 2⟨x_n e^{−i w_q m_n}⟩`)
**hiç değiştirilmeden** import edildi. Dört kesim gazının hepsinde
`L = 12.029593242`, `nline = 3425`, `m`-kimliği `0.0e+00`
(**Ö-b1 ✓**), bant defteri bağıl kalıntısı ≤ **1.5e−13** (**Ö-b2 ✓**,
eşik 1e−12).

### 1a. KESİM AİLESİ — ilk kez ölçüldü

| gaz | kesim | τ̄_A | **`R_η`** | **`R_Ĉ`** | **`R_X`** | `P_η` | `Var(η)` | çizgi payı η |
|---|---|---|---|---|---|---|---|---|
| `E060` | erfc 0.60/0.125 | 0.20493 | **1.22124** | 1.01802 | 1.04751 | 0.06323 | 0.05178 | 0.5855 |
| **`HA4`** | **erfc 0.68/0.125** | 0.22605 | **1.29627** | 1.02270 | 1.06599 | 0.08163 | 0.06297 | 0.5319 |
| `K070` | keskin τ≤0.70 | 0.24347 | **1.29543** | 1.02460 | 1.05629 | 0.09300 | 0.07179 | 0.4992 |
| `K090` | keskin τ≤0.90 | 0.29323 | **1.28842** | 0.98261 | 1.08118 | 0.09189 | 0.07132 | 0.4391 |
| `Hkeskin` | keskin τ≤1.00 | 0.31712 | **1.26504** | 0.99086 | 1.07256 | 0.08593 | 0.06793 | 0.4445 |
| **`son` (gerçek)** | — | — | **1.28829** | 1.02599 | 1.09022 | 0.06997 | 0.05431 | 0.4869 |

### 1b. P1 — KÖPRÜ KESİM EKSENİNDE DE TUTUYOR (**✓ GEÇTİ**)

Ön-kayıt, 174'ün K2-F köprüsünden (`R_η = π_E × r̄`, `r̄ = 0.98709`,
8 λ-gazından) dört kesim gazının `R_η`'sını **ölçümden önce** yazdı:

| gaz | ön-kayıt `R_η^ön` | ölçülen `R_η` | sapma | `\|R_η/π_E − 1\|` |
|---|---|---|---|---|
| `HA4` | 1.29042 | **1.29627** | **+0.45%** | 0.84% |
| `K090` | 1.28915 | **1.28842** | **−0.06%** | 1.35% |
| `K070` | 1.28997 | **1.29543** | **+0.42%** | 0.87% |
| `E060` | 1.21766 | **1.22124** | **+0.29%** | 1.00% |

Ön-kayıt ölçütü %2; ölçülen maks **%1.35** ⇒ **✓ GEÇTİ**, ve dört
öngörünün hepsi **±%0.5 içinde** tuttu.

> **HÜKÜM (P1).** 174'ün "162'nin `R`'si = 172'nin `π`'sidir" köprüsü
> **kesim ekseninde de geçerlidir** — λ ailesinde ölçülen ‑%1.3'lük
> sistematik kayma (çizgi kümesi 3425↔8981, site `m_n`↔`s_n`)
> kesimden bağımsızdır. Köprü artık iki eksende de mühürlüdür.

### 1c. P2 — **AİLE-DIŞI FAZLANIN ADRESİ BULUNDU: KESİM EKSENİ**

174'ün en keskin negatif bulgusu şuydu: gerçek gazın beş niceliği
**λ ailesinin dışındadır** (`R_η`, `m3`, `skew(x1)`, `skew(ds)`,
`skew(e1)`) ve "yeni bir eksen aranmalı". O eksen ölçüldü:

| nicelik | gerçek | λ ailesi [min, max] | λ | KESİM ailesi [min, max] | kesim |
|---|---|---|---|---|---|
| **`R_η`** | 1.28829 | [1.1739, **1.2745**] | **DIŞINDA** | [1.2212, **1.2963**] | **İÇERİDE** |
| **`m3`** | +0.27675 | [0.1502, **0.2516**] | **DIŞINDA** | [0.2333, **0.2911**] | **İÇERİDE** |
| **`skew(x1)`** | +0.28174 | [0.1570, **0.2636**] | **DIŞINDA** | [0.2288, **0.3386**] | **İÇERİDE** |
| **`skew(ds)`** | +0.46806 | [0.2032, **0.4148**] | **DIŞINDA** | [0.3185, **0.7645**] | **İÇERİDE** |
| `skew(e1)` | −0.15129 | [−0.2048, −0.1322] | içeride¹ | [−0.2885, −0.1910] | dışında |
| `m3/m3_çizgi` | −4.20933 | [−4.3605, −3.0083] | içeride | [−38.007, −2.3329] | İÇERİDE |
| `R_Ĉ` | 1.02599 | [0.9408, 1.1073] | içeride | [0.9826, 1.0246] | dışında (+%0.14) |
| **`R_X`** | 1.09022 | [1.0124, 1.1633] | **içeride** | [1.0475, 1.0812] | **dışında** |
| `σ_ds` | 0.40919 | [0.2790, 0.4889] | içeride | [0.4217, 0.4443] | dışında |
| `P_η`, `Var(η)`, çizgi payı η/X̃, `Kov_arası`, `Var(η_çizgi)` | — | — | içeride | — | İÇERİDE |

¹ 174 `skew(e1)`'i "aralık dışı" saymıştı çünkü λ_eş ters çevirmesi
**tekdüze kolda** yapılıyor ve eğri tekdüze değil; [min,max] ölçütüyle
içeridedir. Fark burada kayda geçirilir.

> **BULGU (K1-P2) — 174'ÜN AÇIK BIRAKTIĞI EKSEN, KESİM EKSENİDİR.**
> λ ailesinin taşımadığı **dört** nicelik (`R_η`, `m3`, `skew x1`,
> `skew ds` — hepsi normalize oran ya da üçüncü moment) kesim ailesinin
> **içindedir**. `R_η(λ)` tümseğinin tavanı 1.2745 iken kesim ailesinin
> tavanı **1.2963**'tür ve gerçek gazın 1.28829'u tam ortasına düşer.
>
> **VE ASİMETRİ TERSİNE DÖNER.** 174'ün K2-G'si "X kanalı bir λ-gazıdır"
> demişti (`λ_eş(R_X) = 0.9026 ↔ defterin 0.9089'u`). Şimdi görülüyor ki
> **X kanalı bir λ-gazıdır ama bir kesim-gazı DEĞİLDİR** (`R_X(son) =
> 1.09022`, kesim ailesinin tavanı 1.08118'in üstünde); `σ_ds` de öyle.
> Buna karşılık **η kanalının oranı ve üçüncü momentleri bir kesim-gazıdır
> ama bir λ-gazı değildir.** İki eksen **tamamlayıcıdır** ve gerçek gaz
> her birinde ötekinin taşımadığını taşır.

### 1d. Bant defteri — kesim ailesinde de aynı mekanizma

`P_b/V_b` (bant içi) dört kesim gazında da **< 1** ve bantlar arası
kovaryans **artı** — 174'ün "yıkıcı girişim çizgiler arasında değil,
çizgi ↔ artık arasındadır" hükmü kesim ekseninde de geçerli:

| gaz | iç-bant Σ V_b | bantlar-arası Σ Kov | `Var(η_çizgi)` |
|---|---|---|---|
| `E060` | 0.06957 | +0.02658 | 0.09615 |
| `HA4` | 0.08971 | +0.04005 | 0.12976 |
| `K070` | 0.10236 | +0.04780 | 0.15016 |
| `K090` | 0.10009 | +0.05237 | 0.15246 |
| `Hkeskin` | 0.09350 | +0.04817 | 0.14167 |
| **`son`** | **0.07612** | **+0.03737** | **0.11350** |

`HA4`'ün en üst bantlarında (`τ > 0.75`) `P_b/V_b` **0.986 → 0.995**'e
çıkıyor: erfc zarfı orada çizgileri neredeyse tamamen sildiği için
kalan tek şey uyumlu bir kalıntıdır (`P_b` 2.7e−04 → 2.0e−04; keskin
ikizde 2.3e−03 → 1.5e−03, gerçekte 1.3e−03 → 8.1e−04). **Gerçek gaz
burada da tam ortada.**

---

## 2. K3 — DC KAÇAĞI, KESİM AİLESİNDE (`175c`)

174d'nin makinesi (`ort(E) = Σ_q Re[hp_q κ(ω_q)]`, `ξ_q := Re[hp_q κ_q]`)
değiştirilmeden koştu. **Ön-kayıtlı P3 öngörüsü altı gazda da tam
isabet**: `|ort E| = √(μ̂²_E · K_E/π_E)` özdeşliğinden yazılan değerler
ölçülenle **4.1e−13**'e kadar aynı (ön-kayıt eşiği ‰5):

| gaz | `ort(E)` ölçülen | ön-kayıt | sapma | `μ̂²(ξ)/μ̂²(G1) − 1` |
|---|---|---|---|---|
| `E060` | +0.016825 | +0.016825 | 4.1e−13 | 8.2e−13 |
| **`HA4`** | **+0.031399** | **+0.031399** | 3.1e−13 | 6.2e−13 |
| `K070` | +0.036543 | +0.036543 | 1.3e−13 | 2.6e−13 |
| `K090` | +0.064031 | +0.064031 | 9.3e−14 | 1.9e−13 |
| `Hkeskin` | +0.062139 | +0.062139 | 9.1e−14 | 1.8e−13 |
| **`son`** | **+0.046900** | +0.046900 | 6.4e−14 | 1.3e−13 |

> **HÜKÜM (P3 ✓).** Gerçek gazın DC kaçağı (0.0469) **kesim ailesinin
> içindedir** (0.0168 … 0.0640) ve tam `K070` (0.0365) ile `Hkeskin`
> (0.0621) arasına düşer. λ ailesinde `μ̂²_E` de içerideydi; ama
> **iki eksen gerçeği farklı yerlerden yakalıyor.**

### 2a. Bant-bant DC kaçağı — **SANDVİÇ ÖLÇÜLDÜ**

| τ bandı | n | `Σξ`(gerçek) | `Σξ`(keskin ikiz) | `Σξ`(**erfc ikiz**) | `Σξ`(K090) |
|---|---|---|---|---|---|
| 0.40–0.50 | 58 | −4.9350e−03 | −6.3898e−03 | −4.9456e−03 | −6.5814e−03 |
| 0.50–0.60 | 146 | +9.5090e−03 | +1.0125e−02 | +1.3421e−02 | +1.0963e−02 |
| 0.60–0.70 | 411 | +1.7202e−02 | +2.0580e−02 | +1.7556e−02 | +2.2396e−02 |
| **0.70–0.80** | 1167 | **+1.5384e−02** | **+2.0628e−02** | **+7.4349e−03** | +2.2829e−02 |
| **0.80–0.95** | 7158 | **+9.7396e−03** | **+1.7196e−02** | **−2.0669e−03** | +1.4424e−02 |
| **TOPLAM** | | **+0.046900** | **+0.062139** | **+0.031399** | +0.064031 |

`⟨cos Δφ⟩` (aynı bantlar): gerçek −0.849/+1.000/+0.999/+0.988/**+0.604**;
keskin −0.912/+1.000/+0.999/+0.992/**+0.772**; erfc −0.749/+1.000/
+0.999/+0.941/**−0.192**.

> **BULGU (K3-a) — GERÇEK GAZ İKİ İKİZİN TAM ARASINDADIR.**
> τ > 0.70'te: gerçek **+2.5124e−02**, keskin ikiz **+3.7824e−02**,
> erfc ikiz **+5.3680e−03**. Keskin ikiz gerçeği **aşıyor**, erfc ikiz
> gerçeğin **altında kalıyor** — 164'ün "gerçek `Hkeskin` ile `HA4`
> ARASINDA" sandviçi burada **nicelleşti**.
> Erfc zarfı en üst bantta (`w(0.90) = 0.0064`) çizgileri o kadar
> siliyor ki `Σξ` **işaret değiştiriyor** (−2.07e−03) ve `⟨cosΔφ⟩`
> negatife düşüyor (−0.19): orada ölçülen κ artık merdivenin tepkisi
> değil, bir **tabandır**. Bu, §4'te bağımsız olarak da görülecek.

### 2b. ÖN-MÜHÜR Ö-c3'ün KISMİ ISKASI (kayda geçiyor)

`τ = 0.5` işaret dönüşünün (`⟨cosΔφ⟩ < 0` altta, `> 0` üstte) kesim
ailesinin hepsinde görülmesi beklenmişti. `K090` ✓; ama `HA4`, `K070`
ve `E060`'ın **en üst bandında** `⟨cosΔφ⟩` negatife dönüyor
(−0.192 / −0.233 / −0.379) ⇒ **Ö-c3 kısmen ıskaladı**. Sebep
kurtarma değil ölçümdür: o bantlarda kesim çizgiyi sıfıra sürmüştür
(K070'te tam sıfır, HA4/E060'ta `w < 0.01`), dolayısıyla orada
tarak-çizgi bağıl fazı tanımlı değildir. Kurtarma yapılmadı.

---

## 3. K2 — MÜHÜR ve AYRIŞIM DEFTERİ (`175d`)

### 3a. P4 — MÜHÜR KURALI: **H-K1'in LİTERAL BİÇİMİ ÖLDÜ** (ön-kayıt bunu yazmıştı)

```
Σξ(τ>0.70):  gerçek +2.512385e−02 | keskin ikiz +3.782353e−02
                                   | erfc ikiz  +5.368003e−03
açık(gerçek − keskin) = −1.269968e−02      (174'ün sayısı)
açık(gerçek − erfc)   = +1.975585e−02
|açık_erfc| / |açık_keskin| = 1.556        (mühür eşiği ≤ 0.50)
```

> **HÜKÜM (P4).** Ön-kayıtlı mühür kuralı (**"τ>0.70 açığı ≥%50
> küçülürse H-K1 yaşar"**) **SAĞLANMADI**: açık küçülmedi, **%56
> BÜYÜDÜ ve İŞARET DEĞİŞTİRDİ**. Ön-kayıt bu ölümü ölçümden önce
> yazmıştı ("erfc-ikiz alttan ıskalayacak, işaret dönecek") ve ölçüm
> onu doğruladı — **kurtarma yapılmadı**.
> Ölümün bilgisi şudur: **açık tek bir kesim NOKTASI değildir**;
> `τ_c = 0.68` erfc kesimi **fazla serttir**. Gerçek gaz keskin ile
> erfc arasında bir yerdedir ve o "yer" bantlara göre değişir (§3b).

### 3b. P5 — KESİM KESRİ `f`: **KESİM PAYI ↔ KİLİT PAYI**

Parametresiz tanım (koordinat yok, ara değer yok):
`f(y) := [y(gerçek) − y(keskin ikiz)] / [y(erfc ikiz) − y(keskin ikiz)]`.

| τ bandı | n | `Σξ`(gerçek) | `Σξ`(keskin) | `Σξ`(erfc) | **`f_b`** |
|---|---|---|---|---|---|
| 0.40–0.50 | 58 | −4.9350e−03 | −6.3898e−03 | −4.9456e−03 | **+1.0073** |
| 0.50–0.60 | 146 | +9.5090e−03 | +1.0125e−02 | +1.3421e−02 | **−0.1871** |
| 0.60–0.70 | 411 | +1.7202e−02 | +2.0580e−02 | +1.7556e−02 | **+1.1169** |
| **0.70–0.80** | 1167 | +1.5384e−02 | +2.0628e−02 | +7.4349e−03 | **+0.3975** |
| **0.80–0.95** | 7158 | +9.7396e−03 | +1.7196e−02 | −2.0669e−03 | **+0.3871** |
| **TOPLAM ort(E)** | | +0.046900 | +0.062139 | +0.031399 | **+0.4957** |

Bantlar arası yayılım **1.3040** — ön-kayıtlı "saf kesim" eşiği ±0.15'in
çok üstünde ⇒ **KESİM TEK BAŞINA TAŞIMIYOR**.

> **BULGU (K2-P5) — AMA AÇIĞIN YAŞADIĞI YERDE `f` SABİTTİR.**
> Açığın **%83'ünü** taşıyan iki bantta (τ > 0.70) kesim kesri
> **0.3975 ve 0.3871** — birbirinden **%2.7** farkla aynı. Yani:
>
> > **τ > 0.70'te gerçek−keskin açığının payları:
> > KESİM PAYI = %39,  KİLİT PAYI = %61.**
>
> Alt bantlarda (τ < 0.70, açığın %17'si) `f` kararsızdır
> (+1.01 / −0.19 / +1.12): orada kesim zarfı zaten ≈ 1'dir
> (`w(0.60) = 0.817`), dolayısıyla payda küçüktür ve oran gürültülüdür.

Gözlenebilir bazında `f` (üç bağımsız makineden):

| nicelik | `f` | | nicelik | `f` |
|---|---|---|---|---|
| `R_η` (162 makinesi) | **+0.7445** | | `m3` | +0.8472 |
| `R_Ĉ` | +1.1033 | | `skew(x1)` | +0.7768 |
| `R_X` | **−2.6884** | | `skew(ds)` | +0.5550 |
| `ort(E)` (174d) | **+0.4957** | | `skew(e1)` | −0.6334 |
| `μ̂²_E` | +0.3968 | | `çizgi payı η` | +0.4845 |
| **`g_E`** | **+0.3828** | | `m3/m3_çizgi` | +0.2132 |
| **`g_X`** | **+0.0414** | | `P_η` | +3.7078 |
| **`θ`** | **−0.1430** | | `Var(η)` | +2.7466 |
| **`M = g_E g_X²θ`** | **−0.5566** | | `σ_ds` | −7.5124 |

> **BULGU (K2-P5b) — `g_E`'NİN KESİM PAYI ile DC-KAÇAĞININ KESİM PAYI
> AYNI SAYIDIR.** `f(g_E) = 0.383` (172 defteri, `s_n` sitesi, τ≤0.95)
> ↔ `f`(τ>0.70 DC kaçağı) = **0.387–0.398** (174d makinesi) ↔
> `f(μ̂²_E) = 0.397`. Üç ayrı yoldan **%39**. E kanalının kesim payı
> ölçülmüştür ve tektir.
> Buna karşılık `f(g_X) = 0.04` (X kanalı kesimden neredeyse hiç
> etkilenmiyor) ve **`f(θ) = −0.14`** (kesim θ'yı TERS yöne taşıyor).

### 3c. P6 — ΔM AYRIŞIM DEFTERİ (özdeş, parametresiz)

Özdeşlik `M ≡ g_E · g_X² · θ` (`175d`, kalıntı 3.1e−16); çapa **keskin ikiz**:

```
Δlog M(gerçek ← keskin) = +0.057542  (+5.92%)
   = Δlog g_E +0.024547 + 2Δlog g_X +0.001477 + Δlog θ +0.031518
Δlog M(erfc  ← keskin) = −0.112511
   =            +0.062890 +            +0.035358 +            −0.210759
```

Dondurulmuş formül (`f* := f(R_η) = +0.74447`, 162 makinesinin ölçtüğü):

| pay | Δlog | ΔM'nin yüzdesi |
|---|---|---|
| **KESİM payı** (`f*`·erfc yolu) | **−0.083761** | **−145.6%** |
| **KİLİT payı** (E, X artığı) | **−0.047119** | **−81.9%** |
| **θ payı** (θ artığı) | **+0.188422** | **+327.5%** |
| TOPLAM | +0.057542 | +100.0% (kalıntı 3.1e−16) |

Aynı defterin okunması kolay biçimi — her çarpanın **kendi** kesim kesri:

| çarpan | `Δlog`(gerçek←keskin) | `Δlog`(erfc←keskin) | kesim kesri | kesimin kapattığı |
|---|---|---|---|---|
| `g_E` | +0.024547 | +0.062890 | **0.390** | **%39.0** |
| `2·g_X` | +0.001477 | +0.035358 | **0.042** | **%4.2** |
| `θ` | +0.031518 | −0.210759 | **−0.150** | **−%15.0** |
| **`M`** | **+0.057542** | **−0.112511** | **−0.511** | **−%51.1** |

> **HÜKÜM (K2-P6) — AYRIŞIM CÜMLESİ.**
> **ΔM'nin E-kanalı payının %39'u KESİM ŞEKLİDİR** (üç bağımsız ölçüden
> aynı sayı), **%61'i KİLİTTİR**; X kanalının yalnız **%4**'ü kesimdir;
> **θ'nın kesim payı NEGATİFTİR (−%15)** — kesim ekseni θ'yı gerçeğin
> tersi yöne taşır. Net etki: **kesim ekseni tek başına M'yi −%10.6
> hareket ettirir, gerçek ise +%5.9'dur** ⇒ M kanalında kesim payı
> **−%51**. Yani "fazla kesim şeklindendir" cümlesi **E kanalında
> kısmen (%39) doğrudur ve M bütününde YANLIŞTIR.**

---

## 4. K3 — ÇARPIMSAL κ (`175e`, ön-mühürlü) ve DÜZELTİLMİŞ YASA (`175g`, ön-mühürsüz)

### 4a. Ön-mühürlü sınav: **H-K3 (tam-rezonans çarpımsal) ÖLDÜ**

Merdiven asal kuvvetlerdir ⇒ `q₁q₂ = Q`'nun çözümü yalnız **kule**
çizgilerinde vardır: 8981 çizginin **105'i** (%1.2). Ölçüm
(`175e`, gerçek gaz, `κ̂ = −πAτcos(πτ)·W_pos`, 174e'nin yasası):

| τ bandı | 0.60–0.65 | 0.65–0.70 | 0.70–0.75 | 0.75–0.80 | 0.80–0.85 | 0.85–0.90 | 0.90–0.95 |
|---|---|---|---|---|---|---|---|
| medyan **ASAL** | 1.903 | 2.010 | 2.166 | 2.380 | 2.692 | 3.046 | 3.781 |
| medyan **KULE** | 4.481 | 4.144 | 4.440 | 4.701 | 5.664 | 6.405 | 7.837 |
| `\|κ⁽²⁾/κ⁽¹⁾\|` (kule) | 0.474 | 0.266 | 0.193 | 0.124 | 0.084 | 0.043 | **0.019** |

> **HÜKÜM (K3-a, ön-mühürlü).** Ön-kayıt: *"medyan(ASAL) ≥ 1.5
> (τ>0.60) ⇒ H-K3 ÖLÜR"*. Ölçülen medyan(ASAL) = **2.380** (gerçek),
> **2.595** (keskin ikiz) ⇒ **H-K3 ÖLDÜ.** Fazla, çarpımsal kanalın
> **özdeş sıfır** olduğu asal çizgilerde de vardır; üstelik kulelerde
> hesaplanan çarpımsal terim τ = 0.9'da yalnız **%1.9**'dur.
> 174e'nin "fazlanın adayı çarpımsal çözümlerdir" cümlesi **yanlıştı**.

### 4b. Ön-kayıtlı RAKİP kazandı: fazla `W_pos`'un AŞIRI SÖNÜMÜDÜR

Ön-kayıt bir rakip koymuştu: `κ̂₀ := −π A_Q τ_Q cos(πτ_Q)` (`W_pos` YOK),
ölçüt *"τ>0.55'te medyan ∈ [0.8, 1.25] ⇒ rakip kazanır"*.

| gaz | medyan `\|κ\|/\|κ̂₀\|` (τ>0.55) | hüküm |
|---|---|---|
| **gerçek** | **0.982** | **RAKİP KAZANDI** |
| **keskin ikiz** | **1.019** | **RAKİP KAZANDI** |
| `HA4` (erfc) | 4.798 | rakip de ıska (§4d) |

Asal çizgilerde bant bant (gerçek gaz, `175g`): 1.044 (0.60–0.65),
1.000, 0.965, 0.945, 0.924, 0.899, **0.998** (0.90–0.95).

> **HÜKÜM (K3-c, ön-mühürlü).** 174e'nin ×1.9 → ×4.2 "aşımı"
> **tamamen `W_pos` çarpanının kendisidir**. Doğru yasa
> **`κ(ω_Q) = −π A_Q τ_Q cos(πτ_Q)`** — saf birinci-mertebe doğrusal
> tepki, **sönüm çarpanı YOK** — ve τ = 0.60 … 0.95 arasındaki
> **7000+ asal çizgide ±%10 içinde** tutuyor. 174e'nin `W_pos`
> düzeltmesi bir hataydı: `κ`'nın rezonans terimi tarağın faz
> seğirmesinden **bağımsız değildir** (o seğirmenin kendisi bu
> çizgidir), dolayısıyla `⟨e^{iωḡĈ}⟩` çarpanı oraya takılamaz.

### 4c. **ÖN-MÜHÜRSÜZ İKİNCİ TUR** (`175g`): kulelerin ×2'si ÇARPIMSAL kanaldır

*(Bu bölüm ölçümden SONRA yazıldı ve öyle etiketlenir — 174e'nin aynı
konumu. Yukarıdaki ön-mühürlü hükümler değişmez.)*

Kule çizgileri asallara göre tam **×2.03 – 2.33** fazla veriyor ve
kulelerin medyan çokluğu tam **2.0**. Kaynağı, seviye denkleminin
**ikinci mertebe tersidir**:

```
N̄(z_n) + S(z_n) = n − ½ ,  N̄(z̄_n) = n − ½ ,  d := z_n − z̄_n
N̄'d = −S − S'd + O(d³)   ⇒   d = −ḡS + ḡ² S S' + …
S S' = ½ Σ_{c₁c₂} A₁A₂ω₂ [ sin((ω₁+ω₂)t) + sin((ω₁−ω₂)t) ]
```

`ḡω_Q = 2πτ_Q` ile, yer değiştirmenin `ω_Q` bileşeninin **etkin genliği**

> ```
> A_Q^eff = A_Q + (πτ_Q/2)·Σ_{q₁q₂=Q} A₁A₂ − πτ_Q·Σ_{q₁/q₂=Q} A₁A₂
> κ(ω_Q) ≈ −π τ_Q cos(πτ_Q) · A_Q^eff                    (175g yasası)
> ```

**Sıfır parametre.** Kapalı kehanet: `A_p²/A_{p²} = 2/π` özdeşliğinden
`Q = p²` için ikinci terim birincinin tam **τ_Q** katıdır, `Q = p³` için
**1.5 τ_Q**. Sınav (gerçek gaz, kule çizgileri):

| τ bandı | 0.65–0.70 | 0.70–0.75 | 0.75–0.80 | 0.80–0.85 | 0.85–0.90 | 0.90–0.95 |
|---|---|---|---|---|---|---|
| `\|κ\|/\|κ̂₀\|` (kule) | 2.048 | 1.966 | 1.839 | 1.919 | 1.913 | 2.095 |
| **`\|κ\|/\|κ̂₂\|` (kule)** | **1.192** | **1.127** | **1.030** | **1.021** | **0.985** | **1.067** |
| `A^eff/A` (ölçülen) | 1.691 | 1.746 | 1.778 | 1.841 | 1.878 | 1.935 |
| öngörü `1 + c_k τ_Q` | 1.691 | 1.746 | 1.778 | 1.841 | 1.878 | 1.935 |
| `\|κ\|/\|κ̂₀\|` (**asal**) | 1.000 | 0.965 | 0.945 | 0.924 | 0.899 | 0.998 |

τ > 0.60 medyanları: **gerçek** kule 1.966 → **1.067**; **keskin ikiz**
2.040 → **1.026**; `K090` 2.076 → 1.133.

> **BULGU (K3-d, ön-mühürsüz).** **Kulelerin ×2 fazlası ÇARPIMSAL
> kanaldır ve sıfır parametreyle kapanır.** 143'ün G-yasası köşesindeki
> "toplam frekansı Σ ⇔ q₁q₂ asal-kuvvet — yalnız KULE çiftleri" kuralı,
> `κ`'da doğrudan ölçüldü: `A^eff/A = 1 + τ_Q` (k=2) / `1 + 1.5τ_Q`
> (k=3) kehaneti **altı bantta da birebir** tutuyor.
> Yani KALEM'in H-K3'ü **yanlış yerde** aranmıştı: çarpımsal kanal
> gerçektir ve ölçülür, ama τ>0.5'teki geniş aşımın kaynağı o değildir;
> o aşım `W_pos`'un yanlış takılmasıydı.

### 4d. Erfc gazında yasa neden çöküyor (ölçülmüş bir sınır)

`HA4`'te asal-çizgi oranı τ ile patlıyor: 1.64 (0.60–0.65) → 2.04 →
3.14 → 6.40 → **18.06** (0.80–0.85) → 35.35 (0.85–0.90). Sebep
mekaniktir: erfc zarfı `A_q`'yu sıfıra sürerken ölçülen `|κ|` sonlu bir
**tabana** oturuyor; oran o tabanın merdiven genliğine bölünmesidir.
Bu, §2a'daki iki gözlemin (Σξ'nin işaret değiştirmesi, `⟨cosΔφ⟩`'nin
negatife düşmesi) **bağımsız üçüncü kanıtıdır**: erfc-ikizin en üst
bandında ölçülen şey artık merdivenin tepkisi değildir.

---

## 5. K4 — θ DEFTERİ (`175d` K2-7, `175f`)

| gaz | λ / τ̄_A | `θ`(G1) | `θ`(G2) | `g_E` | `g_X` | `g_cal` | `M` |
|---|---|---|---|---|---|---|---|
| L050 | 0.50 | 0.97697 | 0.96362 | 0.64696 | 0.71134 | 0.32737 | 0.31983 |
| L070 | 0.70 | 0.92405 | 0.91249 | 0.60366 | 0.70331 | 0.29860 | 0.27592 |
| L085 | 0.85 | 0.90424 | 0.89589 | 0.58957 | 0.70279 | 0.29120 | 0.26331 |
| **Hkeskin** | 1.00 | **0.89333** | 0.88844 | 0.58370 | 0.70421 | 0.28947 | 0.25859 |
| L130 | 1.30 | 0.86369 | 0.86677 | 0.58766 | 0.70456 | 0.29172 | 0.25195 |
| — | | | | | | | |
| `E060` | 0.2049 | 0.74729 | 0.70589 | 0.65325 | 0.72494 | 0.34331 | 0.25655 |
| **`HA4`** | 0.2261 | **0.72357** | 0.70857 | 0.62159 | 0.71677 | 0.31935 | 0.23107 |
| `K070` | 0.2435 | 0.71970 | 0.71416 | 0.61144 | 0.71529 | 0.31283 | 0.22515 |
| `K090` | 0.2932 | 0.81167 | 0.80667 | 0.58120 | 0.70353 | 0.28767 | 0.23349 |
| **`son`** | — | **0.92194** | 0.91416 | 0.59821 | 0.70473 | 0.29710 | 0.27391 |

| nicelik | kesim ailesi [min,max] | gerçek | λ ailesi [min,max] |
|---|---|---|---|
| `θ`(G1) | [0.71970, 0.89333] | **DIŞINDA** (üstünde) | [0.86369, 0.97697] içeride |
| `θ`(G2) | [0.70589, 0.88844] | **DIŞINDA** (üstünde) | [0.86677, 0.96362] içeride |
| `g_cal` | [0.28767, 0.34331] | içeride | [0.28947, 0.32737] içeride |
| `M` | [0.22515, 0.25859] | **DIŞINDA** (üstünde) | [0.25195, 0.31983] içeride |

> **HÜKÜM (K4 / P8) — H-K4 ÖLDÜ, ve ölüm keskindir.**
> Ön-kayıt: *"`f(θ) < 0` ⇒ H-K4 ölür"*. Ölçülen **`f(θ) = −0.1430`**.
> Kesim ekseni θ'yı **monoton olmayan ama net biçimde AŞAĞI** taşır
> (0.8933 → 0.8117 → 0.7197/0.7236/0.7473), gerçek gaz ise **yukarıdadır**
> (0.9219). **Gerçek gaz θ'da kesim ailesinin TAMAMININ üstündedir.**
> ⇒ ΔM'nin %51.5'lik θ payı kesim şekliyle **kapanmıyor, BÜYÜYOR**:
> kesim ekseni boyunca gidilirse θ farkı +%3.2'den (gerçek ↔ keskin ikiz)
> +%27.4'e (gerçek ↔ erfc ikiz) çıkar.
> **θ'nın taşıyıcısı kesim şekli DEĞİLDİR.** (`g_cal = g_E g_X²` ise
> kesim ailesinin içindedir — yani ΔM'nin θ-DIŞI kısmı kesimle
> uyumludur, θ değil.)

---

## 6. FİGÜR

`175_kesim_kilit.png` (`175f_K4_figur.py`) — dört panel:
**(a)** `R_η`: λ ailesinin tümseği (tavan 1.2745, gerçek gaz üstünde) ve
KESİM ailesinin tümseği (tavan 1.2963, gerçek gaz **içinde**);
**(b)** DC kaçağının bant defteri — gerçek gazın keskin ve erfc ikizler
arasındaki sandviçi, erfc-ikizin en üst bantta işaret değiştirmesi;
**(c)** kesim kesri `f`'nin gözlenebilir bazında yayılımı (θ ve M'de
ters işaret; DC-kaçağı ağırlıklı `f` = 0.421 çizgisi);
**(d)** `κ` oranları: 174e'nin `W_pos`'lu yasası (τ ile büyüyen sahte
aşım), `W_pos`'suz yasa (asallarda ≈ 1), kulelerin ×2'si ve çarpımsal
düzeltmeyle kapanışı.

---

## 7. DENETİM

| | sınav | sonuç |
|---|---|---|
| **D1** | erfc zarfı 152/164 ile birebir mi? | `164_insa.merdiven` ile **maks fark 0.00e+00** (`175a`) |
| **D2** | erfc-ikizin inşa kapıları | ilk-kök 300000/300000, maks\|F\| **1.86e−09**, sıralılık TAM, `R_bant` min **1.4527** (`175a`) |
| **D3** | 162 makinesi bit düzeyinde mi? | dört kesim gazında `L = 12.029593242`, `nline = 3425`, `m`-kimliği **0.0e+00** (`175b`) |
| **D4** | bant defteri özdeş mi? | 4 gaz × 3 kanal, bağıl kalıntı **≤ 1.5e−13** (eşik 1e−12) |
| **D5** | DC kaçağı özdeşliği | `Σξ` ↔ `√(μ̂²K/π)`: **≤ 4.1e−13** altı gazda; `μ̂²(ξ)/μ̂²(G1) − 1 ≤ 8.2e−13` |
| **D6** | ΔM defteri özdeş mi? | `Δlog M = Δlog g_E + 2Δlog g_X + Δlog θ`, kalıntı **3.1e−16**; ayrışım kalıntısı **3.1e−16** |
| **D7** | ön-kayıt sonradan değişti mi? | `175/ONKAYIT_K2.json` bir kez yazıldı; betik varlık kontrolüyle üzerine yazmayı reddediyor; sha256 raporda |
| **D8** | git | **dokunulmadı** |

### 7a. Dürüstlük notları (kayda geçiyor)

* **Erfc-ikiz yeniden inşa EDİLMEDİ** (§0.1/1). KALEM'in K1 metni
  "167_insa.main ile kur" diyordu; gaz 164'te sadakatli çözücüyle zaten
  kurulmuştu ve çözücü deterministiktir. 10 dakikalık artıksız bir
  tekrar yerine kapılar `z_HA4.npy`'den bağımsız denetlendi. Bu bir
  kestirmedir ve burada açıkça yazılmıştır.
* **σ_ds ön-mühürü GEÇERSİZ ilan edildi** — HA4'ün σ_ds²'si 164'te
  yayımlıydı. "Geniş bant ön-mühür" atmak sahte olurdu. Uzatmanın dört
  ıskası defterde kalır; beşinci bir sahte isabet eklenmedi.
* **172 defteri ön-kayıt anında OKUNMUŞTU** (§0.1/3). Ön-kayıtlı
  öngörüler (P1, P3) o defterden TÜRETİLDİ; sınanan şey, defterin
  **başka iki makinede** tutup tutmadığıydı — ve tuttu (P1 ±%0.5,
  P3 1e−13).
* **`175e` ilk sürümünde `κ⁽²⁾`'ye `|W_pos|` çarpanı elle takılmıştı**;
  koşudan önce "her iki mertebeye AYNI sönüm" biçiminde düzeltildi.
  Ön-mühürlü hükümler (K3-a, K3-c) bu değişiklikten etkilenmez.
* **`175d` iki kez koştu.** İlk koşuda `ort(E)`'nin λ-ailesi kıyası
  eksik ölçümden dolayı yanlış "DIŞINDA" veriyordu (K3 makinesi λ
  merdiveninde yalnız `L085` ve `Hkeskin`'de koşmuştu). Düzeltme:
  `ort(E)` λ merdiveninde **özdeşlikten** (`√(μ̂²K/π)`, 175c'de altı
  gazda 4.1e−13'te doğrulanmış) hesaplandı. Bu bir kurtarma değil, bir
  **hata düzeltmesidir** ve `ort(E)`'yi "λ-dışı kurtarılanlar"
  listesinden **ÇIKARDI** (yani düzeltme hükmü zayıflattı, güçlendirmedi).
* **`175g` ÖLÇÜMDEN SONRA yazıldı** ve "ön-mühürsüz ikinci tur" olarak
  etiketlidir. 175a'nın ön-mühürlü K3 hükümleri (H-K3 öldü, rakip
  kazandı) değişmedi.
* **Ö-c3 kısmen ıskaladı** (§2b) — kurtarılmadı.
* **Bir koşu boşa gitti:** `175c`'yi kuyruğa alan sarmalayıcı
  `pgrep -f 175b_K1` ile bekliyordu ve kendi komut satırı o dizgeyi
  içerdiği için **kendini bekledi**; ~10 dakika kaybedildi, süreç
  öldürülüp `175c` doğrudan koşuldu. Bilimsel bir etkisi yoktur.

---

# HÜKÜM

## (i) AİLE-DIŞI FAZLANIN EKSENİ BULUNDU: **KESİM EKSENİ**

> 174 beş niceliğin "λ ailesinin dışında" olduğunu ölçmüş ve
> *"yeni bir eksen aranmalı"* demişti. 175 o ekseni ölçtü: aynı
> sadakatli çözücüyle kurulmuş **beş kesim gazı** (`E060`, `HA4`,
> `K070`, `K090`, `Hkeskin`). λ ailesinin taşımadığı **dördü** —
> `R_η` (1.2883, λ tavanı 1.2745), `m3` (0.2768, λ tavanı 0.2516),
> `skew(x1)` (0.2817 ↔ 0.2636), `skew(ds)` (0.4681 ↔ 0.4148) —
> **kesim ailesinin İÇİNDEDİR**. `R_η(kesim)` tavanı **1.2963**'tür.
> **İki eksen tamamlayıcıdır:** `R_X`, `R_Ĉ`, `σ_ds`, `θ`, `M`
> λ-ailesindedir ama kesim ailesinin **dışındadır**; η kanalının
> oranları ve üçüncü momentleri tam tersi. 174'ün K2-G'si
> ("X kanalı bir λ-gazıdır") ve K2-A'sı ("η kanalında λ_eş YOKTUR")
> şimdi tek bir cümledir: **X kanalı λ-gazı, η kanalı KESİM-gazıdır.**

## (ii) KÖPRÜ İKİNCİ EKSENDE DE MÜHÜRLENDİ

> `R_η = π_E × r̄` öngörüsü dört kesim gazında **±%0.5** içinde tuttu
> (1.29627/1.29042, 1.28842/1.28915, 1.29543/1.28997, 1.22124/1.21766)
> ve `|R_η/π_E − 1| ≤ %1.35`. 174'ün K2-F köprüsü (`K = ⟨E·e1⟩ ≡
> Σ|hp|²/2`, `g_E ≡ R_η/ρ_E`) kesimden bağımsızdır.

## (iii) MÜHÜR KURALI ÖLDÜ — ama SANDVİÇ nicelleşti

> Ön-kayıtlı kural (**τ>0.70 açığı ≥%50 küçülmeli**) **SAĞLANMADI**:
> açık %56 **büyüdü ve işaret değiştirdi** (−1.270e−02 → +1.976e−02).
> Ön-kayıt bu ölümü önceden yazmıştı. Bilgisi: `τ_c = 0.68` erfc kesimi
> **fazla serttir**; gerçek gaz keskin ile erfc ikizin **arasındadır**
> (164'ün sandviçi), ve o "ara" ölçülebilir bir sayıdır.

## (iv) AYRIŞIM: **KESİM %39, KİLİT %61** (E kanalında)

> Açığın **%83'ünü** taşıyan τ > 0.70 bantlarında kesim kesri
> **0.3975 / 0.3871** — birbirinin %2.7'si içinde SABİT. Aynı sayı
> `g_E`'den (0.390) ve `μ̂²_E`'den (0.397) bağımsız olarak çıkıyor.
> ⇒ **E kanalındaki fazlanın %39'u KESİM ŞEKLİ, %61'i SAF FAZ-KİLİDİ.**
> X kanalında kesim payı yalnız **%4**. Alt bantlarda (τ<0.70, açığın
> %17'si) tek bir kesim kesri yoktur (f: +1.01 / −0.19 / +1.12).

## (v) θ: KESİM DEĞİL — ve şimdi bu KESİN

> **`f(θ) = −0.143`**: kesim ekseni θ'yı gerçeğin **tersi** yöne taşır
> (0.8933 → 0.7236), gerçek gaz ise **kesim ailesinin tamamının
> üstündedir** (0.9219 > 0.8933). ΔM'nin %51.5'lik θ payı kesim şekliyle
> kapanmıyor, **büyüyor**. M bütününde kesim payı **−%51**.
> **H-K4 ÖLDÜ.** 174'ün "θ'nın taşıyıcısı köprüde yoktur" hükmü
> güçlendi: taşıyıcı **kesim ekseninde de yoktur**. (Buna karşılık
> `g_cal = g_E g_X²` kesim ailesinin İÇİNDEDİR — ΔM'nin θ-dışı bütün
> kısmı kesimle uyumludur.)

## (vi) κ'NIN KAPALI YASASI DÜZELTİLDİ — ve ÇARPIMSAL KANAL ÖLÇÜLDÜ

> **H-K3 (ön-mühürlü) ÖLDÜ:** τ>0.5 aşımı asal çizgilerde de var
> (medyan 2.38), oysa çarpımsal kanal orada **özdeş sıfırdır**.
> **Ön-kayıtlı rakip kazandı:** aşımın tamamı 174e'nin taktığı `W_pos`
> çarpanıdır. Doğru yasa
> > **`κ(ω_Q) = −π A_Q τ_Q cos(πτ_Q)`** — sönüm YOK —
>
> ve τ = 0.60…0.95'te **7000+ asal çizgide ±%10**'da tutuyor
> (medyan 0.965).
> **Ön-mühürsüz ikinci tur (`175g`):** kulelerin kalan **×2.0**'ı
> seviye denkleminin ikinci-mertebe tersinden (`d = −ḡS + ḡ²SS'`)
> gelen **ÇARPIMSAL** kanaldır; `A_Q^eff = A_Q + (πτ_Q/2)Σ_{q₁q₂=Q}A₁A₂
> − πτ_Q Σ_{q₁/q₂=Q}A₁A₂` sıfır parametreyle `A^eff/A = 1 + τ_Q`
> (k=2) kehanetini veriyor ve **altı bantta birebir** tutuyor
> (kule oranı 1.97 → **1.07**). 143'ün G-yasası köşesi `κ`'da
> doğrudan ölçüldü.

---

## Sıradaki adım (bu ölçümün işaret ettiği)

1. **θ artık iki eksende de dışarıdadır.** λ (174) ve kesim (175)
   ailelerinin ikisi de θ(gerçek) = 0.9219'u taşımıyor — biri altında
   kalıyor (kesim, tavan 0.8933), öteki içinde ama yanlış λ'da
   (λ_eş(θ) = 0.6916 ↔ λ_eş(σ_ds) = 0.9014). Üçüncü eksen adayı
   174'ün bıraktığı yerde: θ'nın kendi tanımındaki `W`-çözünürlüklü
   oran `⟨e1x1²e^{−iWs}⟩ / ⟨EX²e^{−iWs}⟩`, λ **ve** kesim merdivenlerinde.
2. **`κ`'nın düzeltilmiş yasası (175g) bütün defteri etkiler.**
   `W_pos` çarpanı 166/167'nin kalibrasyonunda da kullanılıyor;
   `κ`'da yanlış olduğu ölçüldüğüne göre, `KALIB`/`θ` zincirinde
   nerede doğru nerede yanlış olduğu yeniden sınanmalı — bu, θ'nın
   %51.5'lik payına doğrudan dokunabilir.
3. **Kesim ekseni yalnız 5 noktalı.** `R_η(kesim)` tümseğinin tepesi
   `HA4`/`K070` (τ̄_A ≈ 0.226–0.243) civarında; gerçek gaz 1.28829 ile
   `K090`–`Hkeskin` kolundadır. Bir erfc `τ_c = 0.80` gazı (tek inşa,
   ~10 dk) tepe ile gerçek arasını doldurur ve `τ_c^eş`'i tek anlamlı
   kılar.
4. **İki eksenli (λ, kesim) okuma yapılmadı.** Her nicelik için tek
   eksende "içeride/dışarıda" bakıldı; iki eksenli bir yüzeyde gerçek
   gazın konumu (ve artığı) parametreli ama bilgilendirici olurdu —
   özellikle `θ`'nın hangi yönde kaçtığını gösterirdi.
5. **Erfc gazının yüksek-τ TABANI ölçüldü ama açıklanmadı** (§4d):
   `A_q → 0` iken `|κ|` sonlu kalıyor. Bu taban, sonlu-N tarağın kendi
   gürültüsü müdür (∝ 1/√N) yoksa ikinci-mertebe kanalların toplamı mı?
   `175g`'nin `A^eff` formülü bunu doğrudan sınayabilir.
