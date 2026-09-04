# 173 — HAKEM GAZLARI: λ = 0.40 ve λ = 1.45
### (4 Eylül 2026, Opus tayfası — görev 173)

**Tek soru:** 172 §G4.4, iki gaz henüz İNŞA EDİLMEDEN, `M` ve `c` için
**iki bağımsız yol** mühürledi ve bunlar λ = 0.40'ta **%6.7** ayrışıyor.
Hangisi doğru? Tek ölçüm hakemdir. İkinci soru: genişletilmiş λ
merdiveninde 171/172'nin resimleri (c vadisi, α tümseği, Q_E tümseği,
ν(λ), M9 çıpası) uçlarda ayakta mı?

*(Merdiven dokuz nokta olarak planlandı; λ = 1.45'in ÖLÇÜM kapısı
düşünce ön-kayıtlı yedek λ = 1.40 da kuruldu ve merdiven **on** nokta
oldu — §H2, §H5.)*

Bu rapor tek başına okunur. Ham çıktılar `scratchpad/173/`, betikler
`173_configs/`.

| betik | kapı | koşu | ham çıktı |
|---|---|---|---|
| `173a_insa_uclar.py` | **H0/H5** — üç gazın inşası (ön-mühürlü) | 9.5 + 9.8 + 9.4 dk | `173/log_insa_L{040,145,140}.txt`, `167/z_L*.npy`, `167/insa_L*.json` |
| `173b_onkayit.py` | **H1** — ölçüm-öncesi ön kayıt + alan girdilerinin hakemi | 3.5-5.6 dk | `173/ONKAYIT_L{040,145,140}.json` |
| `173c_olcum.py` | **H2** — `167_olcum.kos` (ölçüm zinciri kopyalanmadı) | 3.7-3.8 dk | `167/C_L{040,145,140}.json` |
| `173d_hakem.py` | **H3** — 172/G4 mühürleriyle yüzleşme | 20 s | `173/HAKEM.json`, `173/log_173d.txt` |
| `173e_egri_uclari.py` | **H4** — on noktalı merdiven, vadi/tümsekler/ν | 20 s | `173/EGRI.json`, `173/log_173e.txt` |
| `173f_ek_tanilar.py` | **H6** — ön-mühürsüz koklamalar | 1 s | `173/EK.json`, `173/log_173f.txt` |
| `173g_figur.py` | figür | 3 s | `173_hakem.png` |

Toplam: **3 yeni gaz inşası (~29 dk), 3 ölçüm (~11 dk)**, gerisi
önbellek üstünde saniyeler. **Git'e dokunulmadı.**

---

## 0. TEK CÜMLELİK HÜKÜM

> **172'nin λ = 0.40 için mühürlediği iki yol arasındaki %6.7'lik
> ayrışmanın hakemi konuştu: DOĞRUDAN (empirik log-λ) yol +0.9σ ile
> ayakta, ÇARPAN yolu +4.9σ ile ölü — ve ıskanın %84'ü girdilerden
> değil, θ'nın ρ_X-taşıyıcısından geliyor** (G4'ün ekstrapole ettiği
> `g_E`, `g_X`, `ρ_X` ‰0.3-5 isabet etti). Aynı ölçümler `c(λ)`'nın
> vadisini **λ = 0.6486**'da (171'in bağımsız ölçtüğü 0.6487'den
> **0.0001** fark) yeniden bulurken üst uçta **iki yeni sınır** ortaya
> çıkardı: `c(λ)`'nın **λ ≈ 1.279'daki ikinci durgun noktası (bir tepe)**
> ve on gazın ölçülmüş sadakat oranından okunan **inşa sadakat sınırı
> λ_sad ≈ 1.38-1.42**, ki ötesinde (λ = 1.45) 169'un `R_bant ≥ 0.98`
> filtresi düşüyor ve hüküm verilemiyor. **M(λ), λ ∈ [0.40, 1.30]'da
> ölçüm gürültüsünün altında bir log-λ parabolüdür (rms %0.73 <
> %0.83) — ama hâlâ EMPİRİK bir eğridir, kimlik değil.**

---

## H0 — İNŞA: ÖN-MÜHÜR ve KAPILAR (`173a`)

`167_insa.main` AYNEN çağrıldı (170/171 ile birebir aynı yol); tek
değişiklik `KONFIG`'e `L040 = dict(lam=0.40)`, `L145 = dict(lam=1.45)`
satırlarının eklenmesidir. Merdiven, çözücü, seviye konvansiyonu
(c = −½), sıfır sayısı (300000), ızgara adımı (h = 0.015) ve çizgi
sayısı (15450) Hkeskin/L050…L130 ile AYNI.

### H0.1 ÖN-MÜHÜR (koşudan önce, 4 Eylül 2026 13:44 — `173a` docstring'i)

**(A) Merdivenden ÖZDEŞ gelen** (uyum yok, λ ile tam doğrusal):

| büyüklük | λ = 0.40 | λ = 1.45 |
|---|---|---|
| rms S′ = 1.8940·λ | **0.75760** | **2.74630** |
| Σa_qω_q = 259.9·λ | 103.96 | 376.86 |
| rms S = 0.3825·λ | 0.15300 | 0.55463 |
| N̄′ | 1.9147 | 1.9147 |

**(B) Uzatma — ve neden bant GENİŞ ve ÇARPIK.** 171 §T2c.0 "kuvvet-yasası
uzatması λ > 1'de üçüncü kez aşağıdan ıskaladı" demişti. 173a bunu
sayıyla ölçtü (**bir-adım-dışarı sınavı**, log-log kuadratik, yalnız
ölçülmüş yedi noktanın üstünde):

| sınav | σ_ds | σ_X̃ | σ_Ĉ | min Δz |
|---|---|---|---|---|
| 1.30'u öngör (0.50…1.15 ile) | **−1.87%** | **−2.48%** | **−2.62%** | −0.37% |
| 1.30'u öngör (0.50…1.00 ile) | −2.80% | −3.41% | −3.36% | +0.70% |
| 0.50'yi öngör (0.60…1.30 ile) | **+1.97%** | **+2.37%** | **+2.26%** | −0.57% |

Yani uzatma ÜST uçta aşağıdan, ALT uçta yukarıdan ıskalıyor —
log-log'daki eğrilik kuadratiğin taşıdığından fazla. Bantlar bu yüzden
hem geniş hem **çarpık** yazıldı (kurtarma değil: koşudan önce):

| büyüklük | ham kuadratik | çarpıklık düzeltmesi | **ÖN-KAYIT BANDI** |
|---|---|---|---|
| σ_ds(0.40) | 0.23378 | 0.2292 | **[0.224, 0.238]** |
| σ_ds(1.45) | 0.50583 | 0.5160 | **[0.503, 0.527]** |
| σ_X̃(0.40) | 0.13096 | 0.1273 | **[0.1240, 0.1330]** |
| σ_X̃(1.45) | 0.28731 | 0.2913 | **[0.2845, 0.2990]** |
| σ_Ĉ(0.40) | 0.15430 | 0.1478 | **[0.1440, 0.1545]** |
| σ_Ĉ(1.45) | 0.33184 | 0.3336 | **[0.3265, 0.3420]** |
| min Δz(0.40) | 0.22186 | — | **[0.214, 0.233]** |
| min Δz(1.45) | 0.07494 | — | **[0.069, 0.081]** |
| ΔG<0(0.40) | Gauss 0.00575 | — | **[0.0005, 0.005]** |
| ΔG<0(1.45) | Gauss 0.24284 | — | **[0.258, 0.288]** |

**(C) İnşa kapıları — hepsi geçilmeli, kurtarma yok:** (K1) ilk-kök
hücreleri benzersiz 300000/300000; (K2) maks|F| ≤ 1e−8; (K3) sıralılık
TAM. Biri düşerse gaz reddedilir, ölçüm yapılmaz, λ = 1.40'a çekilinir
(yedek ön-kayıt: rms S′ = 2.65160, σ_ds ≈ 0.5100 [0.497, 0.521],
σ_X̃ ≈ 0.2880 [0.2810, 0.2950], σ_Ĉ ≈ 0.3300, min Δz ≈ 0.0776, ΔG<0 ≈ 0.263).

**(D) Yazılı risk (λ = 1.45).** min Δz ≈ 0.075 ile en dar tekne yalnız
~5 ızgara hücresi; F′ zamanın ~%27'sinde negatif. Ölçek kestirimi
min Δz = h ancak λ ≈ 10'da olur ⇒ K1'in gerçek riski düşük; asıl bedel
Newton'un daha çok ikiye-bölme istemesi (süre, kapı değil).

### H0.2 SONUÇ — GERÇEK KOŞUDAN (`173/log_insa_L{040,145,140}.txt`)

**λ = 0.40 (L040, 9.5 dk):** üç kapı da geçti — ilk-kök hücreleri
benzersiz **300000/300000**, maks|F| = **1.863e−9**, sıralılık **TAM**,
ikiye-bölme adımı 0.

| büyüklük | ön-kayıt | ölçülen | hüküm |
|---|---|---|---|
| rms S′ | 0.75760 | **0.7576** | ✓ (özdeş) |
| Σa_qω_q | 103.96 | 103.9 | ✓ |
| rms S | 0.15300 | 0.1530 | ✓ |
| σ_ds | [0.224, 0.238] | **0.22728** | ✓ (alt kenarda) |
| σ_X̃ | [0.1240, 0.1330] | **0.12673** | ✓ |
| σ_Ĉ | [0.1440, 0.1545] | **0.15011** | ✓ |
| min Δz | [0.214, 0.233] | **0.223681** | ✓ (merkez 0.223!) |
| ΔG<0 kesri | [0.0005, 0.005] | **0.0002** | **✗ ISKA (altında)** |
| sıralılık | TAM | TAM | ✓ |

> **ÇARPIKLIK DÜZELTMESİ İŞE YARADI (ilk kez).** Ham log-log kuadratiği
> σ_ds = 0.23378, σ_X̃ = 0.13096, σ_Ĉ = 0.15430 diyordu; ölçüm sırasıyla
> **−%2.8, −%3.2, −%2.7** aşağıda. Ön-kaydın çarpıklık-düzeltmeli
> merkezi (0.2292 / 0.1273 / 0.1478) yalnız **+%0.8 / +%0.4 / −%1.5**
> ıskaladı. Yani "alt uçta uzatma yukarıdan ıskalar" kuralı **dördüncü
> kez** ve şimdi ön-kayıtlı olarak doğrulandı; bant hem geniş hem doğru
> yöne çarpık tutulduğu için üç marjinal de içeride kaldı.
> **Tek ıska ΔG<0 kesridir** (0.0002 vs [0.0005, 0.005]): Gauss
> kestiriminin gözlemle oranı λ = 0.50'de 0.49 idi, λ = 0.40'ta **0.035**
> çıktı — kuyruk kesri λ küçüldükçe Gauss'tan çok daha hızlı sönüyor.
> Kurtarma yok; bant düzeltilmedi.

**λ = 1.45 (L145, 9.8 dk):** üç kapı da geçti — ilk-kök hücreleri
benzersiz **300000/300000**, maks|F| = **1.863e−9**, sıralılık **TAM**;
**36 ikiye-bölme adımı** gerekti (L130'da 6, L050/L040'ta 0).

| büyüklük | ön-kayıt | ölçülen | hüküm |
|---|---|---|---|
| rms S′ | 2.74630 | **2.7463** | ✓ (özdeş) |
| Σa_qω_q | 376.86 | 376.8 | ✓ |
| rms S | 0.55463 | 0.5547 | ✓ |
| σ_ds | [0.503, 0.527] | **0.5444** | **✗ ISKA (+%3.3 ÜSTÜNDE)** |
| min Δz | [0.069, 0.081] | **0.075756** | ✓ (merkez 0.0749) |
| ΔG<0 kesri | [0.258, 0.288] | **0.2717** | ✓ |
| sıralılık | TAM (RİSKLİ) | **TAM** | ✓ (risk yalnız süre oldu) |

> **KURTARMASIZ NOT — "λ > 1'de uzatma aşağıdan ıskalar" DÖRDÜNCÜ KEZ,
> VE BU SEFER ÇOK DAHA SERT.** Ham kuadratik σ_ds(1.45) = 0.50583 dedi
> (**−%7.6**); çarpıklık-düzeltmeli ön-kayıt merkezi 0.5160 dedi
> (**−%5.2**); bandın üst ucu 0.527 bile **−%3.3** altında kaldı.
> **Bant yeterince geniş değildi ve bu bir ıskadır.** Sayının kendisi
> yeni bir olguya işaret ediyor: yerel kuvvet üsteli
> `p = dlogσ_ds/dlogλ` λ_orta = 1.2227'de **0.504** iken
> 1.30 → 1.45 adımında **0.985**'e fırlıyor — neredeyse iki katı ve
> neredeyse tam 1. Yani λ ≳ 1.3'te σ_ds merdiven genliğiyle **doğru
> orantılı** olmaya başlıyor (S'nin N̄'yi ezdiği rejim). 171'in
> "uzatma λ > 1'de daima aşağıdan ıskalar" hükmü artık bir yasa değil,
> bir **rejim değişimi** işareti gibi okunmalı. (min Δz ve ΔG<0 kesri
> aynı gazda tam isabet ettiği için bu bir inşa hatası değil.)

> **İKİ GAZ DA İNŞA KAPILARINDAN GEÇTİ ve kabul edildi.** λ = 1.40
> yedeği İNŞA kapısı düştüğü için değil, λ = 1.45'in **ÖLÇÜM** kapısı
> düştüğü için (§H2) ayrıca kuruldu — §H5.

**λ = 1.40 (L140, 9.4 dk; §H5'in gerekçesiyle sonradan kuruldu):**
üç kapı da geçti (300000/300000, maks|F| = 1.863e−9, TAM, 30
ikiye-bölme). rms S′ = **2.6516** ✓; σ_ds = **0.52231** (ön-kayıt
[0.497, 0.521] — **ıska, +%0.25**); min Δz = **0.078138** ✓
([0.072, 0.084]); ΔG<0 = **0.2654** ✓ ([0.250, 0.278]).

---

## H1 — ÖLÇÜM-ÖNCESİ ÖN KAYIT ve ALAN GİRDİLERİNİN HAKEMİ (`173b`)

`173b_onkayit.py` korelatöre BAKMAZ: yalnız marjinaller
(σ_ds, σ_X̃, σ_Ĉ, bant-başına W_X) ve `Model165.alanlar`'ın
korelatör-öncesi çıktıları (g_E, g_X, Q, ρ, μ̂²). Zaman damgaları
JSON'da. **172'nin mühürleri okundu, yeniden türetilmedi.**

### H1.1 λ = 0.40 — marjinaller ve alanlar (zaman damgası 14:00:34)

```
N = 299998   L = 12.02959   çizgi = 8981
σ_ds = 0.22728   σ_X̃ = 0.12673   σ_Ĉ = 0.15011
g_E = 0.68657   g_X = 0.72094   g_cal = 0.35685
ÖZDEŞLİK  |ΔK|/K = 1.6e−15 (E) / 3.3e−15 (X);  μ̂²_E = +0.01337,
          μ̂²_X = −3.5e−15   [172 §G1.1'in iki mührü de tuttu]
Q_E = 0.56708  ρ_E = 1.80928  rE = 0.31151
Q_X = 0.47940  ρ_X = 1.71790  rX = 0.24090
```

### H1.2 **P6 — 172/G4 ÇARPAN YOLUNUN GİRDİLERİ ÖLÇÜLDÜ** (λ = 0.40)

| girdi | G4 mührü (13:32, gaz yokken) | **ÖLÇÜLEN** | fark |
|---|---|---|---|
| g_E | 0.68415 | **0.68657** | **+0.35%** |
| g_X | 0.71754 | **0.72094** | **+0.47%** |
| ρ_X | 1.71842 | **1.71790** | **−0.03%** |
| Q_E | 0.54484 | **0.56708** | +4.08% |

> **HÜKÜM (H1a — ön-kayıt P6 TUTTU).** 172'nin çarpan yolunun
> **ekstrapole ettiği alan girdileri λ = 0.40'ta neredeyse tam
> isabettir**: g_E ‰3.5, g_X ‰4.7, ve `ρ_X` **‰0.3** — bir gaz
> kurulmadan, log-λ kuadratiğiyle. Ön-kayıt "girdi hataları ≤ %1.5"
> demişti ✓. **Bunun doğrudan sonucu:** λ = 0.40'taki %6.7'lik iki-yol
> ayrışması girdiden gelmiyor; kaynağı **model** olmak zorunda
> (θ = θ₀(ρ_X/ρ_X₀)^0.651 taşıyıcısı). Hata ayrıştırması §H3'te.
> `Q_E`'nin +%4.08'i de ön-kayıt bandı [0.54, 0.68] içindedir (P4 ✓).

*(σ_X̃/σ_Ĉ ön-kayıt bantları §H0.2'de; ikisi de içeride.)*

### H1.3 λ = 1.45 — marjinaller ve alanlar (zaman damgası 14:08:29)

```
σ_ds = 0.54445   σ_X̃ = 0.31782   σ_Ĉ = 0.37422
g_E = 0.60919   g_X = 0.69707   g_cal = 0.29601
ÖZDEŞLİK  |ΔK|/K = 9.9e−16 (E) / 0.0 (X);  μ̂²_E = +0.02067,
          μ̂²_X = −1.1e−14        [iki mühür yine tuttu]
Q_E = 0.55489  ρ_E = 1.41984  rE = 0.66928
Q_X = 0.42012  ρ_X = 1.38687  rX = 0.45338
```

| büyüklük | ön-kayıt | ölçülen | hüküm |
|---|---|---|---|
| σ_X̃ | [0.2845, 0.2990] | **0.31782** | **✗ ISKA (+%6.3)** |
| σ_Ĉ | [0.3265, 0.3420] | **0.37422** | **✗ ISKA (+%9.4)** |

> **REJİM DEĞİŞİMİ — ORANLAR DA KIRILDI.** λ ailesi boyunca
> `σ_X̃/σ_ds` 0.5582 → 0.5668 ve `σ_Ĉ/σ_ds` 0.6498 → 0.6530 arasında
> ‰-mertebesinde donuktu (171 §T2c bunu bir sabit gibi kullanmıştı).
> λ = 1.45'te **0.5837** ve **0.6873** çıkıyor — ikisi de yedi gazın
> tüm yayılımının DIŞINDA. σ'ların ıskası bu yüzden yalnız bir
> "uzatma hatası" değil: **λ ≈ 1.3–1.45 arasında marjinallerin
> iç oranları değişiyor.**

### H1.4 **P6 — G4 GİRDİLERİ, λ = 1.45**

| girdi | G4 mührü | **ÖLÇÜLEN** | fark | ön-kayıt ≤ %1.5 |
|---|---|---|---|---|
| g_E | 0.59391 | **0.60919** | **+2.57%** | **✗** |
| g_X | 0.70852 | **0.69707** | **−1.62%** | **✗** |
| ρ_X | 1.43463 | **1.38687** | −3.33% | — |
| **Q_E** | 0.80917 | **0.55489** | **−31.42%** | — |

> **HÜKÜM (H1b — ÜST UÇTA GİRDİLER DE ÖLDÜ; P6 λ=1.45'te ISKA).**
> λ = 0.40'ta ‰3-5 olan girdi hataları λ = 1.45'te **%1.6-2.6**'ya
> çıkıyor ve `Q_E` **%31 çöküyor**. `Q_E(1.45) = 0.5549` — yalnız
> L130'un 0.8390'ının değil, **λ = 0.40'ın 0.5671'inin bile ALTINDA**
> ve λ = 0.50'nin 0.6903'ünün çok altında. **P4'ün üst-uç bandı
> [0.70, 0.83] ÖLDÜ (kurtarma yok).** Aynı gazda `ρ_E` 2.0348 → 1.4198
> çöküyor ve `rE` 0.5854 → 0.6693'e fırlıyor; `g_E = 1 − Q_E/ρ_E`
> özdeşliği makine hassasiyetinde tutuyor (pay ve payda birlikte
> çöktüğü için g_E aslında YÜKSELİYOR: 0.5877 → 0.6092).
> Yani **Q_E tümseğinin sağ kolu bir yamaç değil, bir uçurumdur** ve
> 172'nin bütün uzatma araçları (5-orta parabol 0.779, 7-nokta 0.739,
> log-λ kuadratiği 0.809) bunu %33-45 ıskaladı.

---

## H2 — ÖLÇÜM ve **λ = 1.45'İN ÖLÇÜM-DÜZEYİ KAPI DÜŞÜŞÜ** (`173c`)

`167_olcum.kos(<gaz>, 0.40, 0.95, 0, 0)` AYNEN — 170/171 ile birebir
aynı çağrı. Her biri 3.7-3.8 dk.

**λ = 0.40 sağlıklıdır.** Hüküm penceresinin beş bandında sadakat
`R_bant` = 1.529 … 1.684, SNR = 245 … 29, hepsi 169'un değişmez
filtresini (R_bant ≥ 0.98, SNR ≥ 3, τ_eff < 0.85) geçti.

**λ = 1.45, 169'UN SAĞLIK FİLTRESİNİ BEŞ BANDIN HEPSİNDE DÜŞÜRDÜ:**

```
R_bant(L145) = 0.9209  0.9303  0.9373  0.9505  0.9598     (eşik 0.98)
```

> **HÜKÜM (H2 — KURTARMA YOK).** Filtre **gevşetilmedi**. λ = 1.45 için
> `c` ve `M` ailenin kendi konvansiyonuyla **hükme bağlanamaz**;
> raporun geri kalanında o satır **TEŞHİS** olarak, açık etiketle
> taşınır. Bu bir başarısızlık değil, bir **ölçümdür**: `R_bant`
> (gerçekleşen/nominal merdiven genliği) λ boyunca tekdüze düşüyor
> — L050 1.56-1.72, Hkeskin 1.29-1.39, L130 1.09-1.15, **L145
> 0.92-0.96** — ve **1'i λ ≈ 1.38 (alt bant) … 1.42 (üst bant)
> arasında kesiyor** (`173f` T3). Yani **164'ün SADAKATLİ çözücüsünün
> genlik sadakati λ ≈ 1.4 civarında biter**; inşa denklemi hâlâ tam
> çözülüyor (maks|F| = 1.9e−9, sıralılık TAM) ama ds'nin çizgi
> genlikleri nominalin altına düşüyor. Aynı sınırın ikinci bağımsız
> işareti §H0.2'deki σ_ds rejim değişimidir (p: 0.50 → 0.99).
> **Ön-kayıtlı yedek λ = 1.40 bu sınırı iki yandan sıkıştırmak için
> ayrıca kuruldu** (§H5).

---

## H3 — **HAKEM: 172/G4'ÜN İKİ YOLU** (`173d`, `173/HAKEM.json`)

### H3.1 λ = 0.40 — **ASIL HAKEM** (sağlıklı gaz)

```
c(L040) = 0.4098 ± 0.0010(jk) ± 0.0048(bant) = ±0.0049   [5 bant]
M(L040) = 1.4086 ± 0.0232      (bant bant 1.3428 1.3740 1.4203 1.4455 1.4639)
θ(orta bant) = 1.02923    g_cal = 0.35685    α_g = 0.6474
```

| yol (mühür 4 Eyl 13:32) | M(mühür) | M ölç−ön | σ | c(mühür) | c ölç−ön | σ | **HÜKÜM** |
|---|---|---|---|---|---|---|---|
| **ÇARPAN yolu** | 1.2947 | **+8.80%** | **+4.9σ** | 0.3790 | +8.12% | +6.3σ | **ÖLDÜ** |
| **DOĞRUDAN yol (M9′)** | 1.3881 | **+1.47%** | **+0.9σ** | 0.4052 | +1.15% | +0.9σ | **AYAKTA** |

> **HÜKÜM (H3a — HAKEM KONUŞTU: DOĞRUDAN YOL KAZANDI).**
> 172 §G4.4 iki yolu mühürlerken kendi tercihini de yazmıştı:
> *"çarpan yolu λ = 0.50'yi zaten −%4.4 ıskaladığı için **doğrudan yol
> daha güvenilirdir**."* **Bu yargı örneklem-dışı doğrulandı:**
> doğrudan yol **+0.9σ** ile ayakta, çarpan yolu **+4.9σ** ile öldü.
> %6.7'lik ayrışmanın hakemi tek bir gazda kapandı.

### H3.2 171'in bütün M ailesi, λ = 0.40'ta (parametreler donduruldu)

| aday | M(ön) | M ölç−ön | σ | hüküm |
|---|---|---|---|---|
| M0 düz | 1.0000 | +40.86% | +17.6σ | ölü |
| M2 kırık kuvvet (λ_c = 1.0109) | 1.1924 | +18.13% | +9.3σ | **ölü** |
| M3 varyans-açığı | 1.1674 | +20.66% | +10.4σ | **ölü** |
| M8 ölçülen çarpan | 1.2328 | +14.26% | +7.6σ | ölü |
| M1 Gram-doyum | 1.1476 | +22.74% | +11.2σ | ölü |
| M7 merdiven kesri | 1.2407 | +13.53% | +7.2σ | ölü |
| **M9 log-λ kuad (171, 5 gaz)** | **1.4112** | **−0.19%** | **−0.1σ** | **TAM İSABET** |
| M9′ log-λ kuad (172, 7 gaz) | 1.3881 | +1.47% | +0.9σ | ayakta |
| G4 çarpan yolu | 1.2947 | +8.80% | +4.9σ | ölü |
| M10 = ölçülen çarpan + θ-taşıyıcı | 1.3113 | +7.41% | +4.2σ | ölü |
| *(çapa) 4/π² — c dilinde* | *0.4053* | *+1.12%* | *+0.9σ* | — |

> **HÜKÜM (H3b — TÜRETİMLİ M AİLELERİNİN İKİNCİ TOPLU ÖLÜMÜ).**
> 171 §T2c-b düşük uçta (λ = 0.50) bütün türetimli adayları
> 9-14σ ile öldürmüştü. λ = 0.40 aynı adayları **7-11σ** ile bir kez
> daha, **daha uzak bir uçta** öldürüyor; eşik ailelerinin (M2/M3)
> yükselen kolu artık **+%18-21** geride. **Ayakta kalan hâlâ yalnız
> EMPİRİK log-λ kuadratiğidir** — ve bu sefer 171'in **beş** gazdan
> (0.60…1.15) kurulmuş M9'u, Δlogλ = 0.405'lik bir uzatmayla
> **−%0.19 (−0.1σ)** isabet ediyor. **Siren kuralı yürürlüktedir:
> bu bir KİMLİK DEĞİLDİR** — ama artık iki bağımsız örneklem-dışı
> noktada (λ = 0.50: +0.4σ; λ = 0.40: −0.1σ) yaşayan tek eğri odur.

### H3.3 **D3 — ÇARPAN YOLU NEREDE ÖLDÜ: GİRDİ Mİ, MODEL Mİ?**

| | MÜHÜR | ÖLÇÜLEN | fark |
|---|---|---|---|
| g_E (girdi) | 0.68415 | 0.68657 | **+0.35%** |
| g_X (girdi) | 0.71754 | 0.72094 | **+0.47%** |
| θ (model: `(ρ_X/ρ_X₀)^0.651`) | 0.94525 | **1.02923** | **+8.88%** |
| M: mühür → ölçülen girdilerle | 1.29469 | 1.31134 | +1.29% |
| M: ölçülen girdiler + model → gerçek | 1.31134 | **1.40857** | **+7.41%** |

```
toplam +8.80%  =  GİRDİ +1.29%  +  MODEL +7.41%
```

> **HÜKÜM (H3c — ÖLÜM ADRESİ: θ TAŞIYICISI).** Çarpan yolunun λ = 0.40
> ıskasının **%84'ü modelden**, yalnız %16'sı girdiden geliyor.
> 172 §G2.1'in "θ'nın en iyi tek-değişkenli taşıyıcısı ρ_X" bulgusu
> (üs 0.651) λ ailesinin İÇİNDE ±%1.4 tutuyordu; **λ = 0.40'ta −%8.18,
> λ = 1.45'te +%27.9 artık veriyor.** Ön-kayıt "|artık| ≤ %3 ⇒ ayakta"
> demişti: **P5 İKİ UÇTA DA ÖLDÜ (kurtarma yok).** 172'nin kendi
> dürüstlük notu ("bu bir KİMLİK DEĞİLDİR, en iyi tek-değişkenli
> taşıyıcıdır") tam olarak doğrulandı — taşıyıcı, uzatıldığında kırılır.

### H3.4 λ = 1.45 — TEŞHİS (hüküm değil; sağlık filtresi düştü)

```
[SAĞLIKSIZ]  c ≈ 0.4098,  M ≈ 0.7663,  θ ≈ 0.6426,  α_g ≈ 0.1866
```

Bu sayılarla **her aday** −16% … −26% (−5σ … −10σ) ıskalar; G4'ün iki
yolu da (çarpan −21.4%, doğrudan −23.8%) uzaktadır. **Ama bu bir hüküm
değildir:** R_bant ≈ 0.94, yani gazın ds çizgi genlikleri nominalin
%6 altında.

**Ancak §H5'in λ = 1.40 gazı bu okumayı düzeltiyor:** λ = 1.40
**SAĞLIKLIDIR** (R_bant = 0.983 … 1.037, beş bandın hepsi geçti) ve
orada bile M = 0.8634 ± 0.0310, yani λ ≤ 1.30 parabolünün **%14.1
altında (−4.6σ)**. **Yani M'nin çöküşü sadakat sınırından ÖNCE
başlıyor** — sadakat kaybı çöküşü büyütüyor ama sebebi değil.
(Nicel not, kurtarma değil: λ = 1.45'te genlik açığı −%6, M açığı −%24;
ilişkinin biçimi ölçülmedi, 174'e açık kapı.)

---


## H4 — EĞRİ UÇLARI: ON NOKTALI MERDİVEN (`173e`, `173/EGRI.json`)

*(Görev "dokuz noktalı merdiven" istemişti; §H5'in sadakat-sınırı gazı
λ = 1.40 ile onuncu nokta eklendi.)*

### H4.0 DEFTER

| gaz | λ | c ± σ | M | g_E | g_X | θ | Q_E | ρ_X | α_g | σ_X̃ | R_bant |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **L040** | **0.40** | **0.4098 ± 0.0049** | **1.4086** | 0.6866 | 0.7209 | **1.0292** | 0.5671 | 1.7179 | 0.6474 | 0.12673 | 1.53-1.68 |
| L050 | 0.50 | 0.3800 ± 0.0072 | 1.2266 | 0.6470 | 0.7113 | 0.9770 | 0.6903 | 1.7046 | 1.0242 | 0.15576 | 1.56-1.72 |
| L060 | 0.60 | 0.3689 ± 0.0072 | 1.1205 | 0.6205 | 0.7058 | 0.9452 | 0.7830 | 1.6834 | 1.1941 | 0.17951 | 1.54-1.70 |
| L070 | 0.70 | 0.3690 ± 0.0058 | 1.0595 | 0.6037 | 0.7033 | 0.9241 | 0.8463 | 1.6561 | 1.2264 | 0.19897 | 1.49-1.63 |
| L085 | 0.85 | 0.3817 ± 0.0028 | 1.0144 | 0.5896 | 0.7028 | 0.9042 | 0.8986 | 1.6094 | 1.1789 | 0.22272 | 1.39-1.51 |
| Hkeskin | 1.00 | 0.4035 ± 0.0018 | 1.0000 | 0.5837 | 0.7042 | 0.8933 | 0.9156 | 1.5624 | 1.0977 | 0.24204 | 1.29-1.39 |
| L115 | 1.15 | 0.4312 ± 0.0053 | 1.0022 | 0.5828 | 0.7066 | 0.8873 | 0.9098 | 1.5180 | 1.0158 | 0.25839 | 1.19-1.27 |
| L130 | 1.30 | **0.4540 ± 0.0136** | 0.9832 | 0.5877 | 0.7046 | 0.8637 | 0.8390 | 1.4767 | 0.6874 | 0.27713 | 1.09-1.15 |
| **L140** | **1.40** | **0.4338 ± 0.0240** | **0.8634** | 0.5998 | 0.6982 | 0.7395 | 0.6624 | 1.4263 | 0.1590 | 0.30150 | **0.98-1.04** |
| *L145* | *1.45* | *0.4098 ± 0.0246* | *0.7663* | *0.6092* | *0.6971* | *0.6426* | *0.5549* | *1.3869* | *0.1866* | *0.31782* | ***0.92-0.96 SAĞLIKSIZ*** |

`g_E = 1 − Q_E/ρ_E` özdeşliği (172 §G1) **on gazın hepsinde** |ΔK|/K ≤
5.3·10⁻¹⁵ ile kapanıyor; `μ̂²_X ≡ 0` mührü de (≤ 1.1·10⁻¹⁴).

### H4.1 **ν(λ) — VE İKİNCİ BİR DURGUN NOKTA** (P1)

| λ adımı | 0.50→0.40 | 0.60→0.50 | 0.70→0.60 | 0.85→0.70 | 1.00→0.85 | 1.15→1.00 | 1.30→1.15 | **1.40→1.30** | *1.45→1.40* |
|---|---|---|---|---|---|---|---|---|---|
| λ_orta | 0.450 | 0.550 | 0.650 | 0.775 | 0.925 | 1.075 | 1.225 | **1.350** | *1.425* |
| **ν (5 bant)** | **+2.203** | +1.489 | +0.993 | +0.563 | +0.205 | −0.034 | +0.271 | **+1.541** | *+1.909* |
| ± | 0.056 | 0.063 | 0.073 | 0.055 | 0.062 | 0.069 | 0.073 | **0.103** | *0.181* |
| ν−1 | **+21.7σ** | +7.7σ | −0.1σ | −8.0σ | −12.8σ | −15.0σ | −10.0σ | **+5.3σ** | *+5.0σ* |

| ön-kayıt (P1) | öngörü | ölçüm | hüküm |
|---|---|---|---|
| ν(0.50→0.40) ∈ [1.55, 2.45] | 1.99 | **+2.203 ± 0.056** | **✓ İÇERİDE** |
| c(0.40) > c(0.50) | — | 0.4098 > 0.3800 | **✓** |
| ν(1.30→1.45) ∈ [0.30, 0.90] | 0.58 | **+1.697 ± 0.063** | **✗ ISKA** |
| c(1.45) > c(1.30) | — | 0.4098 < 0.4540 | **✗** |

> **HÜKÜM (H4a — ALT UÇTA İSABET, ÜST UÇTA YENİ BİR OLGU).**
> **Alt uç:** ν(0.50→0.40) = 2.203, ön-kayıt bandının tam içinde.
> 171 §T2c-d'nin "c'nin tabanı yok, λ\* = 0.649'da bir MİNİMUMU var"
> hükmü **ikinci ve daha uzak bir örneklem-dışı noktayla** doğrulandı:
> c düşen λ ile **hızlanarak** artıyor.
> **Üst uç:** ön-kaydım ν'nün 1'in altında kalacağını söylüyordu;
> **ölçüm ν = 1.541 ± 0.103 verdi ve bu SAĞLIKLI bir adımdır**
> (L130 → L140, ikisi de 169'un filtresini geçiyor). Yani
> **c(λ)'nın λ ≈ 1.28'de bir MAKSİMUMU vardır** ve bu bir sadakat
> artefaktı değildir. Ön-kaydım kurtarılmadan öldü.

### H4.2 **VADİ / EĞİM-KESİŞMESİ — YER ‰0.15 İSABET, SAYI ISKA** (P2)

| λ_orta | dlogKALİB/dλ | dlogW_X/dλ | fark (= dlog c/dλ) |
|---|---|---|---|
| **0.450** | −1.3830 | −0.6277 | **−0.7552** |
| 0.550 | −0.9046 | −0.6077 | −0.2970 |
| **0.650** | −0.5603 | −0.5644 | **+0.0041** |
| 0.775 | −0.2899 | −0.5151 | +0.2253 |
| 0.925 | −0.0953 | −0.4662 | +0.3708 |
| 1.075 | +0.0144 | −0.4280 | +0.4425 |
| 1.225 | −0.1274 | −0.4704 | +0.3430 |
| **1.350** | −1.2989 | −0.8430 | **−0.4559** |
| *1.425* | *−2.3868* | *−1.2506* | *−1.1362* |

```
KESİŞİM(ler) = 0.6486   ve   1.2787       (ön-kayıt: SAYI 1, yer [0.63, 0.67])
```

| ön-kayıt (P2) | ölçüm | hüküm |
|---|---|---|
| kesişim yeri ∈ [0.63, 0.67] | **0.6486** | **✓✓ (171'in λ\* = 0.6487'sinden 0.0001 fark)** |
| kesişim SAYISI = 1 | **2** | **✗ ÖLDÜ** |
| fark(0.45) ∈ [−0.85, −0.40] | **−0.7552** | **✓ İÇERİDE** |
| fark(1.375) ∈ [+0.08, +0.42] | *(1.350'de −0.4559)* | **✗ ISKA** |

> **HÜKÜM (H4b — 172 §G4a DAHA DA KESKİNLEŞTİ; AMA RESİM İKİ VADİLİ
> DEĞİL, BİR VADİ + BİR TEPE).** Vadinin yeri, alt uca yeni bir nokta
> eklendiğinde **λ = 0.6486** çıkıyor; 171'in bağımsız ölçtüğü
> **λ\* = 0.6487** ile fark **0.0001**. 172'nin 0.6506'sı ‰3 idi, şimdi
> **‰0.15**. "Vadi = KALİB'in eğiminin W_X'in sabit eğimini geçtiği
> yer" hükmü artık dört haneli.
> **Ama "tek durgun nokta" ön-kaydım öldü:** λ ≈ **1.2787**'de ikinci
> bir kesişim var ve orada `c` **maksimum** yapıyor (c_max ≈ 0.4547,
> λ = 1.15/1.30/1.40 parabolünden tepe **λ = 1.2786**). Bu ikinci
> durgun nokta **sağlıklı gazlara dayanıyor** (L130 ve L140) — L145'e
> ihtiyaç yok. **c(λ) bir U değil, bir U-artı-tepedir.**
> Mekanizma tabloda açık: 172'nin "dlogW_X/dλ neredeyse sabit" hükmü
> λ ≤ 1.30'un yasasıdır; λ > 1.30'da **payda da çöküyor**
> (−0.47 → −0.84 → −1.25) ve pay ondan daha hızlı çöküyor.

### H4.3 **α(λ) ve A-ÇARPANLAŞMASININ PENCERESİ — 171 KAZANDI, ÖN-KAYDIM ÖLDÜ** (P3)

| gaz | λ | α_g | σ*/2 | S/A2Hk−1 (%, 5 bant) | rms% |
|---|---|---|---|---|---|
| **L040** | 0.40 | **0.6474** | 0.18110 | −5.92 −4.18 0 +1.72 +3.57 | **3.69** |
| L050 | 0.50 | 1.0242 | 0.22779 | −2.61 −2.10 0 −0.09 −0.29 | 1.50 |
| L060 | 0.60 | 1.1941 | 0.24596 | −0.99 −1.32 0 −0.77 −2.05 | 1.23 |
| L070 | 0.70 | **1.2264** | 0.24926 | −0.45 −1.06 0 −0.84 −2.17 | 1.16 |
| L085 | 0.85 | 1.1789 | 0.24438 | −0.33 −0.97 0 −0.53 −1.06 | 0.70 |
| Hkeskin | 1.00 | 1.0977 | 0.23582 | −0.49 −0.95 0 −0.05 +0.48 | **0.53** |
| L115 | 1.15 | 1.0158 | 0.22685 | −0.70 −0.94 0 +0.49 +1.99 | 1.06 |
| L130 | 1.30 | 0.6874 | 0.18661 | −2.01 −3.27 0 +2.45 +6.66 | 3.61 |
| **L140** | 1.40 | **0.1590** | 0.08974 | −3.76 −2.56 0 +5.15 +17.94 | **8.59** |
| *L145* | *1.45* | *0.1866* | *0.09723* | *−0.55 −3.83 0 +6.34 +19.40* | *9.29* |

| ön-kayıt (P3) | öngörü | ölçüm | hüküm |
|---|---|---|---|
| α(0.40) ∈ [0.85, 1.15] | 1.144 | **0.6474** | **✗ ISKA** |
| α(1.45) ∈ [0.35, 0.75] | 0.754 | *0.1866* | **✗ ISKA** |
| şekil rms(0.40) ≤ %2.5 | — | **3.69%** | **✗ ISKA** |
| şekil rms(1.45) > %4 | — | *9.29%* | **✓** |
| σ*/2(1.45) < 0.18661 | — | *0.09723* | **✓** |

> **HÜKÜM (H4c — 171'İN PENCERESİ İKİ YANDAN DA DOĞRULANDI).**
> 171 §T2c-c "**çarpanlaşma λ ∈ [0.50, 1.15] penceresinin yasasıdır**"
> demişti ve o zaman yalnız üst yandan sınanmıştı. **173 alt yanı da
> sınadı: λ = 0.40'ta şekil rms %3.69** — L130'un %3.61'iyle neredeyse
> aynı. Pencere ön-kayıtlı biçimde iki yandan da kapanıyor.
> **α(λ) bir tümsek değil, İKİ YANDAN UÇURUMLU BİR PLATODUR:**
> 0.647 ↗ 1.024 ↗ 1.194 ↗ **1.226** ↘ 1.179 ↘ 1.098 ↘ 1.016 ↘ 0.687
> ↘ **0.159**. Ön-kaydım α'nın uçlarda yumuşak ineceğini varsaydı ve
> **iki uçta da ıskaladı**; kurtarılmadı.
> **172 §G4.2'nin "α tepesi = c vadisi" tam isabeti dar pencereye
> BAĞIMLIDIR:** 5 orta nokta **0.6675** (ön-kayıt λ\* ± 0.03 ✓),
> 7 nokta 0.7877, 9 nokta 0.8058 — parabol plato+uçurum biçimine
> uydurulunca tepe sağa kayıyor. 172 "5 orta nokta" şartını yazmıştı;
> şimdi biliyoruz ki o şart **süs değil, gereklilik**tir.

### H4.4 **Q_E TÜMSEĞİ — SOL KOL DOĞRULANDI, SAĞ KOL UÇURUM** (P4)

| gaz | λ | Q_E | ρ_E | rE | μ̂²_E | g_E |
|---|---|---|---|---|---|---|
| **L040** | 0.40 | **0.56708** | 1.80928 | 0.31151 | +0.01337 | 0.68657 |
| L050 | 0.50 | 0.69027 | 1.95523 | 0.40515 | +0.02016 | 0.64696 |
| L060 | 0.60 | 0.78299 | 2.06334 | 0.47429 | +0.02834 | 0.62052 |
| L070 | 0.70 | 0.84627 | 2.13522 | 0.52052 | +0.03680 | 0.60366 |
| L085 | 0.85 | 0.89861 | 2.18942 | 0.55965 | +0.04815 | 0.58957 |
| Hkeskin | 1.00 | **0.91558** | 2.19932 | 0.57499 | +0.05684 | 0.58370 |
| L115 | 1.15 | 0.90981 | 2.18090 | 0.57603 | +0.06270 | 0.58283 |
| L130 | 1.30 | 0.83903 | 2.03480 | 0.58540 | +0.05786 | 0.58766 |
| **L140** | 1.40 | **0.66240** | 1.65501 | 0.63570 | +0.03410 | 0.59976 |
| *L145* | *1.45* | *0.55489* | *1.41984* | *0.66928* | *+0.02067* | *0.60919* |

| uyum penceresi | Q_E tepesi | λ_c = 1.0109 farkı |
|---|---|---|
| 5 orta (0.60-1.15) | **1.0191** | +0.0081 ✓ |
| 7 nokta (0.50-1.30) | **1.0009** | −0.0100 ✓ |
| 9 nokta (0.40-1.45) | 0.9407 | −0.0702 ✗ |

| ön-kayıt (P4) | ölçüm | hüküm |
|---|---|---|
| Q_E(0.40) ∈ [0.54, 0.68] | **0.56708** | **✓ İÇERİDE** |
| Q_E(1.45) ∈ [0.70, 0.83] | *0.55489* | **✗ ISKA (−%33)** |
| *(yedek)* Q_E(1.40) ∈ [0.71, 0.84] | **0.66240** | **✗ ISKA (−%21), SAĞLIKLI GAZDA** |

> **HÜKÜM (H4d).** 172 §G1.5'in tam isabeti — **Q_E tümseğinin tepesi =
> inşa doyum eşiği λ_c = 1.0109** — dar ve orta pencerede ayakta ve
> **yeni alt-uç noktası (λ = 0.40) tümseğin sol kolunu ön-kayıtlı
> bantta doğruladı**. Üst uçta ise Q_E çöküyor: **λ = 1.40 SAĞLIKLI
> bir gazdır ve orada Q_E zaten −%21 düşmüş.** Yani sağ koldaki çöküş
> bir ölçüm artefaktı değil, gerçek. Çöküş `Q_E` ile `ρ_E`'de birlikte
> oluyor (2.035 → 1.655 → 1.420) ve oranları `g_E`'yi tersine
> **yükseltiyor** (0.5877 → 0.5998 → 0.6092). **Yeni resim: Q_E de bir
> tümsek değil, sağ tarafı uçurumla biten bir tepedir** ve uçurum
> λ ≈ 1.3–1.4'te, tam sadakat sınırının yakınında başlıyor.

### H4.5 **M9 ÖRNEKLEM-DIŞI: M-KİMLİK AVI İÇİN NE SÖYLÜYOR?** (P6)

| uyum | λ = 0.40 öngörü | **ölçüm 1.4086** | λ = 1.40 öngörü | **ölçüm 0.8634** | λ = 1.45 öngörü | *ölçüm 0.7663* |
|---|---|---|---|---|---|---|
| **M9 (171, 5 gaz: 0.60-1.15)** | 1.4112 | **−0.19% (−0.1σ)** | 1.0324 | **−16.4% (−5.5σ)** | 1.0410 | *−26.4% (−9.6σ)* |
| **M9′ (172, 7 gaz: 0.50-1.30)** | 1.3881 | **+1.47% (+0.9σ)** | 1.0013 | **−13.8% (−4.4σ)** | 1.0060 | *−23.8% (−8.4σ)* |
| 8 sağlıklı nokta (0.40-1.30) | *(uyum içi)* | — | 1.0046 | **−14.1% (−4.6σ)** | — | — |

Sekiz **sağlıklı** noktayla (λ = 0.40 … 1.30) uyum (`173f` T1):

```
log M = +0.32810 (logλ)²  −0.07710 logλ  −0.00658
        durgun nokta λ = 1.1247,  M = 0.9890
artıklar (%): +0.30 −0.02 −0.46 −0.49 −0.03 +0.66 +1.32 −1.26
rms 0.73%      (ortalama ölçüm hatası 0.83%)
[karş.] saf kuvvet yasası M = 1.0104·λ^(−0.2839):  rms 4.35%
[karş.] dokuz sağlıklı nokta (1.40 dahil):         rms 3.61%  ⇒ 1.40 UYMUYOR
```

> **HÜKÜM (H4e — M-KİMLİK AVI İÇİN ASIL SONUÇ).**
> **(i) M(λ), λ ∈ [0.40, 1.30] aralığında (3.25 kat) iki parametreli bir
> log-λ PARABOLÜDÜR ve artığı ölçüm gürültüsünün ALTINDADIR
> (%0.73 < %0.83).** Saf kuvvet yasası aynı aralıkta %4.35 ile ölür;
> eğrilik gerçek ve zorunludur. Aranan kimliğin şartı artık çok dar:
> **düşük λ'da parabolik, λ ≈ 1.12'de yassı bir minimum (M ≈ 0.989),
> λ = 1.30'a kadar düz.**
> **(ii) Bu parabolün geçerlilik penceresi de KAPANDI: λ = 1.40, sağlık
> filtresini geçen bir gaz olmasına rağmen parabolün %14.1 altında
> (−4.6σ).** Yani M9 iki uçta değil, **yalnız alt uçta** örneklem-dışı
> yaşıyor. Doğru ifade: **M(λ)'nın log-λ parabolü λ ∈ [0.40, 1.30]'un
> yasasıdır — tıpkı A(τ) çarpanlaşmasının λ ∈ [0.50, 1.15]'in yasası
> olması gibi.** İki pencere de aynı yerlerde bitiyor ve üst uçta ikisi
> de sadakat sınırının (λ ≈ 1.38-1.42) hemen altında kırılıyor.
> **(iii) Siren kuralı yürürlükte: M9 hâlâ EMPİRİKTİR.** Ama üç
> örneklem-dışı sınavda (λ = 0.50: +0.4σ; λ = 0.40: −0.1σ; λ = 1.15/1.30:
> ±1σ) yaşayan tek eğri odur ve 174'ün kimlik avında hedef odur.
> **Türetimli aileler (M0/M1/M2/M3/M7/M8) düşük uçta ikinci kez
> 7-18σ ile öldü** ve hiçbiri kurtarılmadı.

---

## H5 — SADAKAT SINIRI: λ = 1.40 (ön-kayıtlı yedek, `ONMUHUR_L140_Rbant.txt`)

λ = 1.45 **inşa** kapılarını geçip **ölçüm** kapısını düşürünce
(§H2), 173a'nın yedeği devreye alındı ve **koşudan önce** ek bir
ön-mühür yazıldı (`scratchpad/173/ONMUHUR_L140_Rbant.txt`,
4 Eylül 14:14:49):

| ön-kayıt | öngörü | ölçüm | hüküm |
|---|---|---|---|
| **Y1** R_bant(L140): beş bandın **2-3'ü** geçer | 2-3 / 5 | **5 / 5** (0.9829 … 1.0370) | **✗ ISKA** |
| **Y2** 5/5 geçerse ⇒ sadakat sınırı 1.40 ile 1.45 ARASINDA | — | **öyle çıktı** | **✓ (Y2 kolu)** |
| **Y3** σ_ds ≈ 0.5100 [0.497, 0.521] *(ön-kayıt: "muhtemelen yine yukarıdan ıskalanacak")* | — | **0.52231** | **✗ ISKA (+%0.25)**, uyarı tuttu |
| Y3 σ_X̃ [0.2810, 0.2950] | — | **0.30150** | **✗ ISKA (+%2.2)** |
| Y3 σ_Ĉ [0.3230, 0.3380] | — | **0.35174** | **✗ ISKA (+%4.1)** |
| Y3 min Δz [0.072, 0.084] | 0.0776 | **0.078138** | **✓** |
| Y3 ΔG<0 [0.250, 0.278] | 0.263 | **0.2654** | **✓** |
| **Y4** λ = 1.40 için G4 mührü YOKTUR, uydurulmayacak | — | uydurulmadı | ✓ |

İnşa kapıları: hücre benzersizliği **300000/300000**, maks|F| =
**1.863e−9**, sıralılık **TAM**, 30 ikiye-bölme adımı (9.4 dk).

```
R_bant(λ):  0.40  1.53-1.68 | 0.50  1.56-1.72 | 0.60  1.54-1.70 | 0.70  1.49-1.63
            0.85  1.39-1.51 | 1.00  1.29-1.39 | 1.15  1.19-1.27 | 1.30  1.09-1.15
            1.40  0.98-1.04 | 1.45  0.92-0.96  ← 169 filtresi burada düşüyor
R_bant = 1 geçişi:  λ ≈ 1.3841 (alt bant)  …  1.4240 (üst bant)
169 filtresinin (0.98) geçişi: λ ≈ 1.4369
```

> **HÜKÜM (H5 — YENİ BİR SABİT: İNŞA SADAKAT SINIRI λ ≈ 1.38–1.42).**
> 164'ün SADAKATLİ çözücüsü inşa denklemini λ = 1.45'te de tam çözüyor
> (maks|F| = 1.9e−9, sıralılık TAM), **ama ürettiği gazın ds çizgi
> genlikleri λ ≈ 1.4'ün üstünde nominalin altına düşüyor.**
> `R_bant(λ)` on gaz boyunca tekdüze ve düzgün iniyor ve **1'i
> λ ≈ 1.38-1.42'de kesiyor**. Bu sayı hiçbir uyumdan gelmedi; on
> bağımsız gazın ölçülmüş sadakat oranından çıktı. **170'in inşa doyum
> eşiği λ_c = 1.0109'un (rms S′ = N̄′) yanına ikinci bir λ-ölçeği
> geliyor: λ_sad ≈ 1.40.** İkisinin oranı 1.385; ön-mühürsüz koklama:
> λ_sad/λ_c ≈ √2 = 1.414 (%2 fark) — **sınanmadı, iddia değil.**
> Pratik sonuç: **λ ailesinin ölçülebilir üst sınırı λ ≈ 1.40'tır**;
> bunun ötesindeki gazlar kurulabilir ama hüküm taşıyamaz.

---

## H6 — EK TANILAR (`173f`, **ÖN-MÜHÜRSÜZ KOKLAMALAR — iddia değil**)

172 §6.4'ün kuralı gereği açıkça etiketlenmiştir: aşağıdakiler
ön-kayıtlı DEĞİLDİR ve hüküm sayılmazlar; 174'e nottur.

**K1 — `θ(λ) = 1` GEÇİŞİ λ ≈ 0.456.** θ ≡ KALİB_orta/g_cal ve
172 §G2.0'a göre `θ = 1` ⟺ **beş artık-alan teriminin toplamı SIFIR**.

```
λ:  0.40    0.50    0.60    0.70    0.85    1.00    1.15    1.30    1.40   1.45
θ:  1.0292  0.9770  0.9452  0.9241  0.9042  0.8933  0.8873  0.8637  0.7395 0.6426
```

λ = 0.40, ailenin θ > 1 olan **ilk** gazıdır; doğrusal ara değer
**λ ≈ 0.4559**. Yani artık-alan eşleşmesi orada işaret değiştiriyor.
Bu, λ\* = 0.6487'den (c vadisi) ve λ_c = 1.0109'dan (inşa doyumu)
**ayrı, üçüncü bir λ-ölçeğidir**. Sınanmadı.

**K2 — `4/π²` İKİ KEZ KESİLİYOR.** c(0.40) = 0.4098 (+%1.12, +0.9σ) ve
c(1.00) = 0.4035 (−%0.43, −1.0σ); arada c minimumu −%9'a iniyor,
sonra λ = 1.30'da +%12'ye çıkıyor. **170 §K2′'nin "4/π² bir limit değil,
λ = 1.00'in tesadüfüdür" hükmü üçüncü kez ve şimdi çok daha güçlü
doğrulandı: `4/π²` bir limit değil, `c(λ)`'nın İKİ KEZ kestiği bir
seviye çizgisidir** (λ ≈ 0.43 ve λ ≈ 1.01; sağlıksız L145 üçüncü kez
oradan geçiyor).

**K3 — σ MARJİNALLERİNİN İÇ ORANLARI ve `p(λ)`.**

| λ | σ_ds | σ_X̃/σ_ds | σ_Ĉ/σ_ds | p = dlogσ_ds/dlogλ |
|---|---|---|---|---|
| 0.40 | 0.22728 | 0.55759 | **0.66049** | +0.9194 |
| 0.50 | 0.27903 | 0.55821 | 0.64980 | +0.7745 |
| 0.60 | 0.32135 | 0.55861 | 0.64509 | +0.6617 |
| 0.70 | 0.35586 | 0.55914 | 0.64328 | +0.5719 |
| 0.85 | 0.39765 | 0.56009 | **0.64290** | +0.5004 |
| 1.00 | 0.43134 | 0.56114 | 0.64376 | **+0.4538** |
| 1.15 | 0.45958 | 0.56223 | 0.64521 | +0.5044 |
| 1.30 | 0.48890 | 0.56684 | 0.65295 | +0.8921 |
| 1.40 | 0.52231 | 0.57723 | 0.67344 | +1.1828 |
| 1.45 | 0.54445 | **0.58374** | **0.68734** | — |

`p(λ)` bir U'dur ve **minimumu λ ≈ 1.07'dedir** — `ν(λ)`'nün minimumu
(λ_orta = 1.075, ν = −0.034) ve `dlogKALİB/dλ`'nın maksimumu
(λ_orta = 1.075, +0.0144) ile **aynı yerde**. Üç bağımsız büyüklüğün
durgun noktası çakışıyor; sınanmadı.

**K4 — M parabolünün katsayısı ve θ üssü.** 8 sağlıklı noktanın uyumu
`log M = 0.32810 (logλ)² − ...`; **2 × 0.32810 = 0.6562**, 172 §G2.1'in
θ taşıyıcı üssü **0.651** ile ‰8 içinde aynı sayı. Sınanmadı.

---

## 7. HÜKÜM TABLOSU — 173'ÜN BÜTÜN ÖN-KAYITLARI

| # | ön-kayıt | eşik/öngörü | ölçüm | hüküm |
|---|---|---|---|---|
| 173a-A | rms S′ = 1.8940·λ (üç gaz) | 0.7576 / 2.7463 / 2.6516 | aynen | ✓✓ |
| 173a-K1 | ilk-kök hücreleri benzersiz | 300000/300000 | üç gazda da | ✓ |
| 173a-K2 | maks\|F\| ≤ 1e−8 | — | 1.863e−9 (üçü de) | ✓ |
| 173a-K3 | sıralılık TAM (λ=1.45 RİSKLİ) | — | TAM (üçü de) | ✓ |
| 173a-B | σ_ds(0.40) ∈ [0.224, 0.238] | çarpık bant | 0.22728 | ✓ |
| 173a-B | σ_X̃(0.40) ∈ [0.1240, 0.1330] | — | 0.12673 | ✓ |
| 173a-B | σ_Ĉ(0.40) ∈ [0.1440, 0.1545] | — | 0.15011 | ✓ |
| 173a-B | min Δz(0.40) ∈ [0.214, 0.233] | 0.223 | 0.223681 | ✓✓ |
| 173a-B | ΔG<0(0.40) ∈ [0.0005, 0.005] | 0.002 | **0.0002** | ✗ |
| 173a-B | **σ_ds(1.45) ∈ [0.503, 0.527]** | 0.5160 | **0.5444** | **✗** |
| 173a-B | **σ_X̃(1.45) ∈ [0.2845, 0.2990]** | 0.2913 | **0.31782** | **✗** |
| 173a-B | **σ_Ĉ(1.45) ∈ [0.3265, 0.3420]** | 0.3336 | **0.37422** | **✗** |
| 173a-B | min Δz(1.45) ∈ [0.069, 0.081] | 0.0749 | 0.075756 | ✓ |
| 173a-B | ΔG<0(1.45) ∈ [0.258, 0.288] | 0.272 | 0.2717 | ✓✓ |
| 173a-D | λ=1.45 riski: yalnız süre mi? | — | 36 ikiye-bölme, kapı düşmedi | ✓ |
| 173b-P6 | **G4 girdi hataları ≤ %1.5 (λ=0.40)** | — | **g_E +0.35%, g_X +0.47%, ρ_X −0.03%** | **✓✓** |
| 173b-P6 | G4 girdi hataları ≤ %1.5 (λ=1.45) | — | g_E +2.57%, g_X −1.62% | ✗ |
| 173b-P4 | Q_E(0.40) ∈ [0.54, 0.68] | — | 0.56708 | ✓ |
| 173b-P4 | Q_E(1.45) ∈ [0.70, 0.83] | — | 0.55489 | ✗ |
| 173b-P4 | Q_E tepesi = λ_c ± 0.05 (9 nokta) | 1.0109 | 0.9407 *(sağlıksız kaldıraç)* | ✗ |
| 173b-P5 | **θ taşıyıcısı \|artık\| ≤ %3** | — | **−8.18% (0.40), +27.9% (1.45)** | **✗ İKİ UÇTA DA** |
| 173b-P1 | **ν(0.50→0.40) ∈ [1.55, 2.45]** | 1.99 | **+2.203** | **✓** |
| 173b-P1 | c(0.40) > c(0.50) | — | 0.4098 > 0.3800 | ✓ |
| 173b-P1 | ν(1.30→1.45) ∈ [0.30, 0.90] | 0.58 | **+1.697** | ✗ |
| 173b-P1 | c(1.45) > c(1.30) | — | 0.4098 < 0.4540 | ✗ |
| 173b-P2 | **kesişim yeri ∈ [0.63, 0.67]** | 0.6487 | **0.6486** | **✓✓✓** |
| 173b-P2 | kesişim SAYISI = 1 | 1 | **2** (0.6486, 1.2787) | ✗ |
| 173b-P2 | fark(0.45) ∈ [−0.85, −0.40] | −0.60 | −0.7552 | ✓ |
| 173b-P2 | fark(1.375) ∈ [+0.08, +0.42] | +0.24 | −0.4559 (1.350'de) | ✗ |
| 173b-P3 | α(0.40) ∈ [0.85, 1.15] | 1.144 | **0.6474** | ✗ |
| 173b-P3 | α(1.45) ∈ [0.35, 0.75] | 0.754 | 0.1866 | ✗ |
| 173b-P3 | şekil rms(0.40) ≤ %2.5 | — | 3.69% | ✗ |
| 173b-P3 | şekil rms(1.45) > %4 | — | 9.29% | ✓ |
| 173b-P3 | σ*/2(1.45) < 0.18661 | — | 0.09723 | ✓ |
| 173d-D1 | **172'nin tercihi ("doğrudan yol") doğru mu?** | — | **+0.9σ vs +4.9σ** | **✓✓✓** |
| 173d-D2 | hata kestirimi ±0.008 (M, λ=0.40) | — | ±0.0232 (3 kat büyük) | ✗ |
| 173d-D3 | **ıska modelden gelmeli, girdiden değil** | — | **%84 model, %16 girdi** | **✓✓** |
| 173-Y1 | R_bant(L140): 2-3 bant geçer | 2-3/5 | **5/5** | ✗ |
| 173-Y2 | 5/5 ⇒ sınır 1.40-1.45 arasında | — | öyle | ✓ |
| 173-Y3 | σ_ds(1.40) ∈ [0.497, 0.521] *(uyarılı)* | — | 0.52231 | ✗ |
| 173-Y3 | min Δz(1.40) ∈ [0.072, 0.084] | 0.0776 | 0.078138 | ✓ |
| 173-Y3 | ΔG<0(1.40) ∈ [0.250, 0.278] | 0.263 | 0.2654 | ✓ |
| 173-Y4 | λ=1.40 için G4 mührü uydurulmayacak | — | uydurulmadı | ✓ |

**Sayım (43 satır): 24 ✓, 19 ✗.** Hiçbir ölen kurtarılmadı, hiçbir eşik
koşudan sonra gevşetilmedi, 169'un sağlık filtresi (R_bant ≥ 0.98)
gevşetilmedi, 172'nin mühürleri değiştirilmeden yüzleştirildi.

---

## 8. DENETİM ve DÜRÜSTLÜK NOTLARI

1. **Hiçbir inşa/ölçüm parçası kopyalanmadı.** `173a` yalnız
   `167_insa.KONFIG`'e üç satır ekleyip `167_insa.main`'i çağırır;
   `173c` yalnız `167_ortak.KUNYE`'ye künye ekleyip
   `167_olcum.kos(<gaz>, 0.40, 0.95, 0, 0)` çağırır — 170/171 ile
   birebir aynı çağrı. `173d/e/f` `169_k1`'i aynen import eder.
   **Git'e dokunulmadı.**
2. **172'nin mühürleri okundu, yeniden türetilmedi.** `173b`
   `172/G4.json`'u dosyadan okur; bu raporda geçen 1.29469 / 1.38809 /
   0.97435 / 1.00599 sayıları 4 Eylül 13:32'de yazılmış olanlardır.
   (Denetim: `173b`'nin `ρ_X` ve `g_X` uzatmaları G4'ün sayılarını
   5 hanede yeniden üretti — iki kod yolu tutarlı.)
3. **Zaman damgaları:** ONKAYIT_L040 **14:00:34**, ONKAYIT_L145
   **14:08:29**, ONKAYIT_L140 **14:27:58**; ilgili ölçümler
   sırasıyla 14:05, 14:09 ve 14:28'de başladı. λ = 1.40'ın R_bant
   ön-mührü **14:14:49** (inşadan önce).
4. **169'un sağlık filtresi gevşetilmedi.** λ = 1.45'in sayıları
   raporun her yerinde *italik* ve "SAĞLIKSIZ" etiketiyle taşındı;
   hiçbir hüküm ona dayandırılmadı. §H4.1/H4.2'nin ikinci durgun
   noktası **sağlıklı L130-L140 çiftinden** okundu.
5. **Ön-kayıtlarımın kendi eşikleri iki kez fazla dar çıktı:**
   (a) σ(1.45) bandı (%3-9 ıska) — "çarpıklık düzeltmesi" alt uçta
   çalıştı, üst uçta yetmedi; (b) α ve şekil-rms bantları uçlarda
   yumuşak iniş varsaydı, gerçek uçurumdu. İkisi de yazıldı,
   düzeltilmedi.
6. **173d-D2'nin hata kestirimi 3 kat küçüktü** (±0.008 dedim,
   ±0.023 çıktı) çünkü bant-bant M yayılımı uçlarda büyüyor
   (λ=0.40'ta 1.343…1.464). Hüküm bundan etkilenmedi: çarpan yolu
   daha büyük hatayla bile +4.9σ ile ölüyor.
7. **Ön-mühürsüz koklamalar §H6'da toplandı ve etiketlendi**
   (θ = 1 geçişi, 4/π²'nin iki kesişimi, p(λ)/ν(λ) durgun noktalarının
   çakışması, 2a ≈ p_θ, λ_sad/λ_c ≈ √2). Hiçbiri iddia değildir.
8. **λ = 1.40 bir "kurtarma" değildir.** λ = 1.45'in ölçüm kapısı
   düştükten SONRA, 173a'nın kendi yedek ön-kaydı devreye alındı ve
   koşudan önce yeni bir ön-mühür (Y1-Y4) yazıldı. λ = 1.45'in ölümü
   geri alınmadı; λ = 1.40 onun yerine geçmiyor, sınırı sıkıştırıyor.

---

## 9. NE KAZANILDI, NE KALDI

**Kazanılanlar.**
1. **HAKEM KAPANDI.** 172'nin iki yolu λ = 0.40'ta **%6.7** ayrışıyordu;
   ölçüm **doğrudan yolu (+0.9σ)** seçti, **çarpan yolunu (+4.9σ)**
   öldürdü. 172'nin kendi yazılı tercihi doğrulandı.
2. **ÖLÜM ADRESİ ÖLÇÜLDÜ:** çarpan yolunun ıskasının **%84'ü model**
   (θ = θ₀(ρ_X/ρ_X₀)^0.651 taşıyıcısı), **%16'sı girdi**. G4'ün
   ekstrapole ettiği alan girdileri λ = 0.40'ta ‰0.3-5 isabetliydi.
3. **θ taşıyıcısı iki uçta da öldü** (−%8.2 ve +%27.9) — 172'nin
   "taşıyıcı, kimlik değil" dürüstlük notu doğrulandı.
4. **VADİ DÖRT HANEDE:** eğim-kesişmesi **λ = 0.6486**, 171'in bağımsız
   λ\* = **0.6487**'si ile fark **0.0001**.
5. **c(λ) bir U DEĞİL, U + TEPE:** λ ≈ 1.279'da bir maksimum var
   (sağlıklı gazlardan). ν(1.30→1.40) = 1.541 ± 0.103 (+5.3σ).
6. **171'in A-çarpanlaşma penceresi [0.50, 1.15] iki yandan da
   doğrulandı** (şekil rms: 0.40'ta %3.69, 1.30'da %3.61).
7. **Q_E tümseğinin tepesi = λ_c mührü sol kolda doğrulandı**
   (Q_E(0.40) ön-kayıtlı bantta); sağ kol bir uçurum çıktı.
8. **YENİ SABİT: İNŞA SADAKAT SINIRI λ_sad ≈ 1.38-1.42** — `R_bant(λ)`
   on gazda tekdüze inip 1'i orada kesiyor. λ ailesinin ölçülebilir
   üst sınırı budur.
9. **M(λ) = log-λ parabolü, λ ∈ [0.40, 1.30]'da gürültü altında
   (rms %0.73 < ölçüm %0.83)**; λ = 1.40'ta (sağlıklı gaz) −4.6σ ile
   kırılıyor. Kimlik avının hedefi ve penceresi netleşti.
10. Türetimli M aileleri düşük uçta **ikinci kez** toplu öldü (7-18σ).

**Kalanlar (174'ün kapıları).**
1. **M-KİMLİĞİ HÂLÂ YOK.** Aranan şey artık çok dar tanımlı:
   λ ∈ [0.40, 1.30]'da `log M ≈ 0.328 (logλ)² − 0.077 logλ`, λ ≈ 1.12'de
   yassı minimum (M ≈ 0.989), sonra λ_sad'a doğru çöküş. Türetimli
   hiçbir aile bunu vermiyor.
2. **λ_sad ≈ 1.40'ın ilk-ilke türetimi.** `R_bant = 1` koşulu neyin
   eşiği? λ_c = 1.0109 (rms S′ = N̄′) gibi kapalı bir koşul var mı?
   (Koklama: λ_sad/λ_c ≈ √2.)
3. **θ'nın gerçek kimliği.** ρ_X^0.651 taşıyıcısı uçlarda kırıldı;
   θ = 1 geçişi (λ ≈ 0.456) yeni ve sert bir çapa sunuyor: artık-alan
   eşleşmesinin işaret değiştirdiği λ.
4. **c'nin ikinci durgun noktası λ ≈ 1.279'un adresi.** Vadi (0.6486)
   "KALİB eğimi = W_X eğimi" ile açıklanmıştı; tepe için aynı cebir
   `dlogW_X/dλ`'nın SABİT OLMAKTAN ÇIKTIĞI yeri istiyor.
5. **172'nin asıl kapısı hâlâ açık:** `μ̂²_E` ve `Q_E`'yi κ cinsinden
   DOĞRUDAN ölçmek (asal-oran frekanslarında). 173 bu kapıya
   girmedi; onun yerine iki hakem gazı kurdu.

---

## 10. FİGÜR

`173_hakem.png` (`173g_figur.py`) — dört panel:
**(a)** `c(λ)`'nın U + tepe biçimi, vadi λ\* = 0.6486, 4/π² seviyesi ve
G4'ün iki `c` mührü; **(b)** hakem: `M(λ)`, sekiz sağlıklı noktanın
log-λ parabolü, G4'ün iki mührü (çarpan +4.9σ ÖLDÜ, doğrudan +0.9σ
AYAKTA), 171'in M ailesinin λ = 0.40 öngörüleri ve λ = 1.40'ın
parabolden −4.6σ kopuşu; **(c)** `α(λ)` platosu (171'in [0.50, 1.15]
penceresi taralı) ile `Q_E(λ)` tümseği ve λ_c; **(d)** `ν(λ)` merdiveni
(ν = 1 çizgisi) ve `R_bant(λ)` sadakat sınırı.

*(Tek cümlelik hüküm için §0'a bakınız.)*
