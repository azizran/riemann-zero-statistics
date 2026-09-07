# 181 — VS3/VS4 TOHUMLARI

### (7 Eylül 2026, Opus tayfası — 180'in nedensel defterini n = 4 tohuma büyütme)

> **Bu rapor tek başına okunur.** Her sayının yanında onu üreten betiğin
> adı vardır. Dar kapsamlı, tek amaçlı sefer: **yalnız VS ailesinin
> tohum sayısı** 2 → 4. Kırpma hakemi (K2 / VF ailesi / `169_k2b`)
> **bu görevin dışındadır** — dokunulmadı, tek sayısı yeniden
> hesaplanmadı, **DEĞİŞMEDİ**; 180'in H-180a **HÜKÜMSÜZ** hükmü aynen
> durur. Ölümler kurtarılmadan, tutmayan öngörüler süslenmeden yazılır.

**Betikler:** `181_configs/181a_onkayit_n4.py` (K0 ön-kayıt),
`180_configs/180b_vs_insa.py` + `180c_olcum.py` (**aynen çağrıldı,
kopyalanmadı**), `181_configs/181c_defter_n4.py` (n = 4 defteri; kestirici
ve nicelik listeleri `180f_defter.py`'den ithal, **üreme kapısıyla**).

---

## 0. K0 — n = 4 UZANTISININ DONMUŞ ÖN-KAYDI (`181a_onkayit_n4.py`)

**Zaman:** `Mon Sep  7 23:18:17 2026`
**Betik sha256:** `3b16cbdf31fc4e2030fb0e28e53768d02f189cf2e2f41a5bde8e7b6599515a4d`
**Dosya:** `181/ONKAYIT_181_K0.json`
**Devraldığı mühürler:** 180a ön-kaydı
`22af1b06c6d6cc26b5b19876646d99e3078d05af129e620d54363f06ee01126d`;
180'in `K3_DEFTER.json`'u (`Mon Sep 7 22:55:07 2026`).

Bu betik **VS3/VS4 doğmadan önce** koştu; hiçbir yeni yüzleşme niceliği
hesaplamadı. Üç iş yaptı: zarfın bit-bit aynılığını doğruladı,
okunurluk ölçütünü **sayısallaştırıp kalibre etti**, n = 4'ün hata
çubuklarını ve eşiklerini **önceden yazdı**.

### 0.1 Zarf — bit-bit aynı

`A_sonzarf.npy` sha256 =
`4ac2f82748a2a7ac8973cc38d514d4650bede00014a098c8963242455d0e52e1`
= 180a'nın ön-kayıtta dondurduğu dize. **VS3 ve VS4, VS1/VS2 ile
tamamen aynı genlik dizisini kullanır; tek fark φ_q'dur.**

### 0.2 Devralınan (GEVŞETİLMEYEN) eşikler — hepsi 180a'nın

* İki konvansiyon (doğrusal + log), `ZARF := ⟨N⟩_VS − ⟨N⟩_VF`,
  `KİLİT := ΔN − ZARF`;
* **hata çubuğu:** `se(⟨N⟩_A) = sd(N_A, ddof=1)/√n_A`,
  `se(ZARF) = √(se_VS² + se_VF²)`, `se(KİLİT) = se(ZARF)`,
  `se(p) = se(ZARF)/|ΔN|`;
* **işaret kuralı:** ZARF ile ΔN aynı işaretli değilse pay > 1 ya da
  < 0 çıkar; bu bir hata değil ölçümdür, aynen yazılır;
* **inşa kapıları V1–V6:** sha(A) = ön-kayıttaki, `maks|A − a_Hk·r(τ)| = 0`,
  φ≡0 sağlaması `0.0/0.0`, ilk-kök hücre **300000/300000**,
  `maks|F| ≤ 1e−8` (aşan 0), sıralılık **TAM**.

### 0.3 TEK YENİ ŞEY: okunurluk ölçütü, veriden önce ve kalibre edilmiş

180, "okunur ↔ HÜKÜMSÜZ" ayrımını raporunun §4.1'inde **sözel**
uyguladı. n = 4'e çıkarken bu ayrımın **sayısal** hâli gerekliydi.
Uydurulmadı: 180'in **kendi partisyonunu birebir üreten** ölçüt
donduruldu ve kalibrasyonu veri gelmeden gösterildi.

```
O1 (KARAR VEREN):  |ZARF_N| ≥ 2 · se(ZARF_N)
O2 (yalnız raporlanır, karar vermez):  se(p_zarf)
```

**Kalibrasyon (`181a`, VS3/VS4 yokken):** O1, 180 §4.1'in partisyonunu
— okunur = {M, Q_E, ρ_E, g_E}, HÜKÜMSÜZ = {Q_X, ρ_X, θ, g_X, g_cal} —
**18/18 satırda, iki konvansiyonda da birebir** üretti. En dar okunur
satır ρ_E `|ZARF|/se = 4.29`; en geniş HÜKÜMSÜZ satır g_cal `1.68`.
Ölçüt tutmasaydı betik ölecekti.

### 0.4 ÖN-KAYITLI ÖNGÖRÜLER (`181a`, sayılardan önce)

**H1** (koşul: tohum popülasyonunun sd'si n = 2 kestirimindeki değerinde
kalırsa): `se_VS(n=4) = sd_VS(n=2)/2 = se_VS(n=2)/√2`,
`se_ZARF(n=4) = √(se_VS(n=2)²/2 + se_VF²)` — **se_VF değişmez** (VF
ailesi zaten n = 4).

| satır | konv. | se_p (n=2) | **se_p (n=4) öngörü** | öngörülen hüküm |
|---|---|---|---|---|
| **M** | log | %11.3 | **%9.4** | OKUNUR |
| **M** | doğrusal | %18.7 | %15.3 | OKUNUR |
| **Q_E** | log | %7.0 | %5.7 | OKUNUR |
| **Q_E** | doğrusal | %2.3 | %1.9 | OKUNUR |
| **ρ_E** | log | %33.8 | %27.2 | OKUNUR |
| g_cal | doğrusal | %188 | %147 | **OKUNUR (yeni)** |
| g_cal | log | %108 | %85 | **OKUNUR (yeni)** |
| g_X | doğ./log | %2243 / %1910 | %1763 / %1503 | HÜKÜMSÜZ (sd 2.7× küçülmeli) |
| Q_X, ρ_X, θ | ikisi de | — | — | **ULAŞILAMAZ** |

**H2** (yapısal): bir satır ancak
`sd_VS(n=4) ≤ sd*_VS := 2√((|ZARF|/2)² − se_VF²)` olduğunda okunur
olabilir. Karekökün içi negatifse (**|ZARF| < 2·se_VF**) satır **VS
tohumu eklemekle ULAŞILAMAZ**dır — n ne olursa olsun, çünkü
`se_ZARF ≥ se_VF` her zaman:

> **Q_X, ρ_X ve θ satırları — iki konvansiyonda da — VS tohumu ekleyerek
> okunur hâle GETİRİLEMEZ.** Onları boğan şey VS'nin az tohumu değil,
> **VF ailesinin kendi saçılımıdır.**

---

## 1. K1 — İNŞA: VS3 ve VS4 (`180b_vs_insa.py`, aynen)

Makine 180'inkinin **aynısı**: alan değerlendiricisi
`176_vekil_cekirdek` (faz taşır), çözücü `164_insa.coz_sadakatli`
(h = 0.015, nz = 300000, c = −½). Yeni satır yazılmadı; betik olduğu
gibi `VS3 3` ve `VS4 4` ile çağrıldı.

### T1 — İnşa kapıları (eşikler 180a'dan, gevşetme yok)

| kapı | eşik | `VS3` (tohum 3) | `VS4` (tohum 4) |
|---|---|---|---|
| V1 `sha256(A_VS)` | ön-kayıttaki | `4ac2f827…52e1` ✓ | `4ac2f827…52e1` ✓ |
| V2 `maks\|A − a_Hk·r(τ)\|` | = 0.0 | **0.0e+00** ✓ | **0.0e+00** ✓ |
| V3 φ≡0 (ΔS / ΔS′) | = 0.0 | **0.0 / 0.0** ✓ | **0.0 / 0.0** ✓ |
| V4 ilk-kök hücre | 300000/300000 | **300000/300000** ✓ | **300000/300000** ✓ |
| V5 maks\|F\| | ≤ 1e−8 | **1.863e−09**, aşan 0 ✓ | **1.863e−09**, aşan 0 ✓ |
| V6 sıralılık | TAM | **TAM**, min Δz = 0.099232 ✓ | **TAM**, min Δz = 0.090580 ✓ |
| σ_ds | — | 0.41923 | 0.41808 |
| L | — | 12.029593229 | 12.029593229 |
| süre | — | 10.1 dk | 11.9 dk |

> **DÖRT TOHUMDA DA V1–V6 GEÇTİ.** Genlik dizisi dördünde de bit-bit
> aynı (tek sha), `L` dördünde de `12.029593229`; σ_ds dört tohumda
> 0.41808–0.41989 (VF: 0.4416/0.4427). Vekiller arasındaki **tek fark
> φ_q dizisidir**.

### Ölçüm (`180c_olcum.py`, aynen: `167_olcum.kos(ad, 0.40, 0.95, 0, 0)`)

Özdeşlik denetimi `|ΔK|/K`: VS3 E 3.97e−15 / X 4.30e−15; VS4 E 4.92e−15
/ X 1.42e−16 (172b'nin Ö-A özdeşliği kırılmadı). Süre: VS3 9.7 dk
(VS4 inşasıyla eşzamanlı koştu), VS4 5.1 dk.

---

## 2. K2 — DEFTER, n = 4 (`181c_defter_n4.py`)

### 2.0 ÜREME KAPISI (R-KAPISI) — makinenin 180f'inki olduğunun kanıtı

180f'in `main()`'i VS listesini gövdesinde sabitlediği için n = 4
döngüsü 181c'de kuruldu; kestirici (`istat`) ve nicelik listeleri
180f'den **ithal edildi**. Buna karşılık aynı döngü VS = (VS1, VS2) ile
koşuldu ve 180'in `K3_DEFTER.json`'uyla karşılaştırıldı:

> **270 alan (delta, ZARF, KİLİT, se, paylar, K4 dâhil) — 0 fark, TAM
> SIFIR.** Tutmasaydı betik ölecekti. n = 4 tablosu, 180'in makinesinin
> ürettiği tabloyla aynı koddan çıkmıştır.

### T2 — DEFTER (HAM esas | oran EK), dört VS tohumu

| gaz | M | Q_E | ρ_E | Q_X | ρ_X | \| g_E | g_X | θ | g_cal |
|---|---|---|---|---|---|---|---|---|---|
| `son` | 0.273907 | 0.875748 | 2.179597 | 0.467931 | 1.584774 | 0.598206 | 0.704733 | 0.921938 | 0.297099 |
| `Hkeskin` | 0.258590 | 0.915576 | 2.199323 | 0.462127 | 1.562364 | 0.583701 | 0.704213 | 0.893334 | 0.289467 |
| `VF1` | 0.422643 | 0.297426 | 1.130635 | 0.193999 | 1.040067 | 0.736939 | 0.813474 | 0.866671 | 0.487663 |
| `VF2` | 0.425711 | 0.297109 | 1.128368 | 0.164799 | 0.988779 | 0.736692 | 0.833331 | 0.832136 | 0.511588 |
| `VF3` | 0.430532 | 0.299440 | 1.133933 | 0.190991 | 1.035086 | 0.735928 | 0.815483 | 0.879713 | 0.489401 |
| `VF4` | 0.428004 | 0.298609 | 1.136974 | 0.204761 | 1.059888 | 0.737365 | 0.806809 | 0.891712 | 0.479981 |
| `VS1` | 0.444879 | 0.288527 | 1.150329 | 0.189736 | 1.051637 | 0.749179 | 0.819580 | 0.884044 | 0.503232 |
| `VS2` | 0.449528 | 0.287084 | 1.144479 | 0.160029 | 1.000126 | 0.749157 | 0.839991 | 0.850423 | 0.528594 |
| **`VS3`** | **0.450379** | **0.290417** | **1.152932** | **0.185700** | **1.045695** | 0.748105 | 0.822415 | 0.890089 | 0.505993 |
| **`VS4`** | **0.450341** | **0.289891** | **1.157151** | **0.199024** | **1.069210** | 0.749478 | 0.813859 | 0.907159 | 0.496429 |

### T3 — AYRIŞIM, **DOĞRUSAL** konvansiyon, n_VS = 4

| N | ΔN | ⟨N⟩_VS | ⟨N⟩_VF | **ZARF ± se** | **KİLİT ± se** | p_zarf ± se | p_kilit ± se | \|Z\|/se | hüküm |
|---|---|---|---|---|---|---|---|---|---|
| **M** | +0.015316 | 0.448782 | 0.426723 | **+0.022059 ± 0.002133** | **−0.006743 ± 0.002133** | **+1.4403 ± 0.1393** | **−0.4403 ± 0.1393** | 10.34 | OKUNUR |
| **Q_E** | −0.039828 | 0.288980 | 0.298146 | −0.009166 ± 0.000921 | −0.030661 ± 0.000921 | +0.2301 ± 0.0231 | +0.7699 ± 0.0231 | 9.95 | OKUNUR |
| **ρ_E** | −0.019726 | 1.151223 | 1.132478 | +0.018745 ± 0.003253 | −0.038471 ± 0.003253 | −0.9503 ± 0.1649 | +1.9503 ± 0.1649 | 5.76 | OKUNUR |
| **Q_X** | +0.005804 | 0.183622 | 0.188637 | −0.005015 ± 0.011896 | +0.010819 ± 0.011896 | −0.864 ± 2.050 | +1.864 ± 2.050 | 0.42 | **HÜKÜMSÜZ** |
| **ρ_X** | +0.022410 | 1.041667 | 1.030955 | +0.010712 ± 0.021047 | +0.011698 ± 0.021047 | +0.478 ± 0.939 | +0.522 ± 0.939 | 0.51 | **HÜKÜMSÜZ** |
| *(EK)* g_E | +0.014505 | 0.748980 | 0.736731 | +0.012249 ± 0.000426 | +0.002256 ± 0.000426 | +0.8445 ± 0.0294 | +0.1555 ± 0.0294 | 28.76 | OKUNUR |
| *(EK)* g_X | +0.000520 | 0.823961 | 0.817274 | +0.006687 ± 0.007987 | −0.006167 ± 0.007987 | +12.85 ± 15.35 | −11.85 ± 15.35 | 0.84 | **HÜKÜMSÜZ** |
| *(EK)* θ | +0.028604 | 0.882929 | 0.867558 | +0.015371 ± 0.017519 | +0.013233 ± 0.017519 | +0.537 ± 0.613 | +0.463 ± 0.613 | 0.88 | **HÜKÜMSÜZ** |
| *(EK)* g_cal | +0.007632 | 0.508562 | 0.492158 | +0.016404 ± 0.009734 | −0.008772 ± 0.009734 | +2.149 ± 1.276 | −1.149 ± 1.276 | 1.69 | **HÜKÜMSÜZ** |

### T4 — AYRIŞIM, **LOG** konvansiyon, n_VS = 4

| N | ΔlogN | **ZARFlog ± se** | **KİLİTlog ± se** | p_zarf ± se | p_kilit ± se | \|Z\|/se | hüküm |
|---|---|---|---|---|---|---|---|
| **M** | +0.057542 | **+0.050413 ± 0.004912** | **+0.007129 ± 0.004912** | **+0.8761 ± 0.0854** | **+0.1239 ± 0.0854** | 10.26 | OKUNUR |
| **Q_E** | −0.044475 | −0.031232 ± 0.003155 | −0.013243 ± 0.003155 | **+0.7022 ± 0.0709** | **+0.2978 ± 0.0709** | 9.90 | OKUNUR |
| **ρ_E** | −0.009010 | +0.016413 ± 0.002842 | −0.025423 ± 0.002842 | −1.8217 ± 0.3154 | +2.8217 ± 0.3154 | 5.78 | OKUNUR |
| **Q_X** | +0.012481 | −0.027025 ± 0.066103 | +0.039506 ± 0.066103 | −2.17 ± 5.30 | +3.17 ± 5.30 | 0.41 | **HÜKÜMSÜZ** |
| **ρ_X** | +0.014242 | +0.010357 ± 0.020507 | +0.003885 ± 0.020507 | +0.727 ± 1.440 | +0.273 ± 1.440 | 0.51 | **HÜKÜMSÜZ** |
| *(EK)* g_E | +0.024547 | +0.016490 ± 0.000574 | +0.008057 ± 0.000574 | +0.6718 ± 0.0234 | +0.3282 ± 0.0234 | 28.75 | OKUNUR |
| *(EK)* g_X | +0.000739 | +0.008151 ± 0.009689 | −0.007412 ± 0.009689 | +11.04 ± 13.12 | −10.04 ± 13.12 | 0.84 | **HÜKÜMSÜZ** |
| *(EK)* θ | +0.031518 | +0.017622 ± 0.020192 | +0.013896 ± 0.020192 | +0.559 ± 0.641 | +0.441 ± 0.641 | 0.87 | **HÜKÜMSÜZ** |
| *(EK)* g_cal | +0.026024 | +0.032791 ± 0.019266 | −0.006767 ± 0.019266 | +1.260 ± 0.740 | −0.260 ± 0.740 | 1.70 | **HÜKÜMSÜZ** |

### T8 — n = 2 → n = 4: **merkez mi kaydı, çubuk mu sıkıştı?**

| satır | konv. | p_zarf(2) | **p_zarf(4)** | Δmerkez | se_p(2) | se_p(4) öngörü | **se_p(4)** | çubuk× | hüküm(2) → (4) | öngörü |
|---|---|---|---|---|---|---|---|---|---|---|
| **M** | log | +0.8149 | **+0.8761** | **+0.0612** | 0.1133 | 0.0936 | **0.0854** | **0.75** | OKUNUR → OKUNUR | ✓ |
| **M** | doğrusal | +1.3372 | +1.4403 | +0.1030 | 0.1872 | 0.1534 | 0.1393 | 0.74 | OKUNUR → OKUNUR | ✓ |
| **Q_E** | log | +0.7936 | **+0.7022** | **−0.0914** | 0.0695 | 0.0569 | **0.0709** | **1.02** | OKUNUR → OKUNUR | ✓ |
| **Q_E** | doğrusal | +0.2596 | +0.2301 | −0.0295 | 0.0226 | 0.0186 | 0.0231 | 1.02 | OKUNUR → OKUNUR | ✓ |
| **ρ_E** | log | −1.4535 | −1.8217 | −0.3682 | 0.3379 | 0.2723 | 0.3154 | 0.93 | OKUNUR → OKUNUR | ✓ |
| **ρ_E** | doğrusal | −0.7567 | −0.9503 | −0.1936 | 0.1764 | 0.1418 | 0.1649 | 0.94 | OKUNUR → OKUNUR | ✓ |
| g_E | log | +0.6820 | +0.6718 | −0.0102 | 0.0167 | 0.0167 | 0.0234 | **1.40** | OKUNUR → OKUNUR | ✓ |
| g_E | doğrusal | +0.8574 | +0.8445 | −0.0130 | 0.0208 | 0.0208 | 0.0294 | **1.41** | OKUNUR → OKUNUR | ✓ |
| Q_X | log | −6.102 | −2.165 | +3.937 | 7.770 | 6.092 | 5.296 | 0.68 | HÜKÜMSÜZ → HÜKÜMSÜZ | ✓ |
| Q_X | doğrusal | −2.370 | −0.864 | +1.506 | 2.947 | 2.326 | 2.050 | 0.70 | HÜKÜMSÜZ → HÜKÜMSÜZ | ✓ |
| ρ_X | log | −0.346 | +0.727 | +1.073 | 2.044 | 1.620 | 1.440 | 0.70 | HÜKÜMSÜZ → HÜKÜMSÜZ | ✓ |
| ρ_X | doğrusal | −0.226 | +0.478 | +0.704 | 1.331 | 1.054 | 0.939 | 0.71 | HÜKÜMSÜZ → HÜKÜMSÜZ | ✓ |
| g_X | log | +20.57 | +11.04 | −9.53 | 19.10 | 15.03 | 13.12 | 0.69 | HÜKÜMSÜZ → HÜKÜMSÜZ | ✓ |
| g_X | doğrusal | +24.05 | +12.85 | −11.19 | 22.43 | 17.63 | 15.35 | 0.68 | HÜKÜMSÜZ → HÜKÜMSÜZ | ✓ |
| θ | log | −0.0072 | **+0.5591** | **+0.5663** | 0.7770 | 0.6438 | 0.6407 | 0.83 | HÜKÜMSÜZ → HÜKÜMSÜZ | ✓ |
| θ | doğrusal | −0.0113 | **+0.5374** | **+0.5487** | 0.7401 | 0.6124 | 0.6125 | 0.83 | HÜKÜMSÜZ → HÜKÜMSÜZ | ✓ |
| **g_cal** | log | +1.811 | +1.260 | −0.551 | 1.081 | 0.850 | 0.740 | 0.69 | HÜKÜMSÜZ → HÜKÜMSÜZ | **✗ ÖNGÖRÜ TUTMADI** |
| **g_cal** | doğrusal | +3.113 | +2.149 | −0.963 | 1.885 | 1.474 | 1.276 | 0.68 | HÜKÜMSÜZ → HÜKÜMSÜZ | **✗ ÖNGÖRÜ TUTMADI** |

### 2.1 Hangi boğulan satır okunur oldu? — **HİÇBİRİ. Tek tek:**

| satır | n=2 | **n=4** | neden |
|---|---|---|---|
| **Q_X** (doğ./log) | HÜKÜMSÜZ | **HÜKÜMSÜZ** | ön-kayıtta **ULAŞILAMAZ** yazılmıştı (\|ZARF\| < 2·se_VF); `\|Z\|/se` 0.80 → **0.42** ve 0.79 → **0.41** — *uzaklaştı* |
| **ρ_X** (doğ./log) | HÜKÜMSÜZ | **HÜKÜMSÜZ** | ön-kayıtta **ULAŞILAMAZ**; `\|Z\|/se` 0.17 → 0.51 |
| **θ** (doğ./log) | HÜKÜMSÜZ | **HÜKÜMSÜZ** | ön-kayıtta **ULAŞILAMAZ**; merkez çok kaydı (pay −0.01 → **+0.56**) ama `\|Z\|/se` yalnız 0.015 → **0.88** — eşiğin (2) hâlâ çok altında |
| **g_X** (doğ./log) | HÜKÜMSÜZ | **HÜKÜMSÜZ** | ön-kayıt "sd 2.7× küçülmeli" demişti; sd 0.78× küçüldü, yetmedi (`\|Z\|/se` 1.07 → 0.84) |
| **g_cal** (doğ./log) | HÜKÜMSÜZ | **HÜKÜMSÜZ** | **ön-kayıt OKUNUR öngörmüştü — TUTMADI.** Çubuk beklendiği gibi daraldı (×0.69), ama **ZARF merkezi de küçüldü** (log +0.0471 → +0.0328); `\|Z\|/se` 1.68 → **1.70**, yerinde saydı |

> **Ön-kayıtın iki öngörüsünden biri tuttu, biri tutmadı ve tutmadığı
> yazılıyor.** Tutan: Q_X/ρ_X/θ'nın **yapısal ulaşılamazlığı** (VF'nin
> kendi saçılımı boğuyor) — dördü de aynen HÜKÜMSÜZ kaldı. Tutmayan:
> `g_cal`'ın okunur olacağı — H1'in "ZARF merkezi sabit kalır"
> varsayımı yanlıştı; yeni tohumlar merkezi de aşağı çekti.

### T7 — K4: karıştırma çöküşleriyle çapraz tutarlılık (n = 4)

| N | K_N (Hk zarfı) | K′_N (son zarfı, n=4) | log K − log K′ | KİLİTlog | kalıntı | *(180 n=2: K′_N)* |
|---|---|---|---|---|---|---|
| **M** | +%65.01 | **+%63.84** | +0.007129 | +0.007129 | 0.0e+00 | +%63.27 |
| **Q_E** | −%67.44 | **−%67.00** | −0.013243 | −0.013243 | 4.2e−17 | −%67.14 |
| **ρ_E** | −%48.51 | −%47.18 | −0.025423 | −0.025423 | 2.8e−17 | −%47.36 |
| Q_X | −%59.31 | −%60.89 | +0.039506 | +0.039506 | 0.0e+00 | −%62.76 |
| ρ_X | −%34.03 | −%34.29 | +0.003885 | +0.003885 | −6.9e−18 | −%35.29 |
| *(EK)* g_E | +%26.22 | +%25.20 | +0.008057 | +0.008057 | 0.0e+00 | +%25.24 |
| *(EK)* θ | −%2.92 | −%4.26 | +0.013896 | +0.013896 | 0.0e+00 | −%5.95 |

Özdeşlik makine sıfırında duruyor. 176'nın manşet çöküşleri n = 4'te de
birebir: `Q_E` −%67, `ρ_E` −%48.5, `g_E` +%26, `M` +%65. **Kilit iki
zarfta da neredeyse aynı** (M: +%65.0 ↔ +%63.8; Q_E: −%67.4 ↔ −%67.0)
— 180'in (iii) hükmü ayakta.

### T9 — Zarf aktarımının sağlaması (tanı, kapı değil), n = 4

| τ_eff | ⟨R_bant⟩_VS | ⟨R_bant⟩_VF | ölçülen oran | tarif `r(τ)` | fark |
|---|---|---|---|---|---|
| 0.4607 | 0.9676 | 1.0098 | 0.9582 | 0.9460 | +1.3% |
| 0.5392 | 0.9446 | 0.9896 | 0.9545 | 0.9323 | +2.4% |
| 0.6185 | 0.9099 | 0.9605 | 0.9473 | 0.9178 | +3.2% |
| 0.6984 | 0.8611 | 0.9093 | 0.9470 | 0.9062 | +4.5% |
| 0.7776 | 0.7850 | 0.8345 | 0.9407 | 0.8917 | +5.5% |

180'in ölçülmüş sınırı **n = 4'te de aynen duruyor ve keskinleşiyor**:
ölçülen oran τ boyunca **düz** (0.941–0.958), tarif ise **eğimli**
(0.892–0.946); fark τ ile büyüyor (+%1.3 → +%5.5).

---

## 3. K3 — KISA DEĞERLENDİRME

### 3.1 Manşet cümlesi değişti mi? — **Merkez de kaydı, çubuk da sıkıştı**

| | 180 (n=2) | **181 (n=4)** | ne oldu |
|---|---|---|---|
| **M, log** | %81.5 ± %11.3 zarf | **%87.6 ± %8.5 zarf** | merkez **+6.1 puan** (eski se'nin 0.54'ü), çubuk **×0.75** |
| **M, log kilit** | %18.5 ± %11.3 | **%12.4 ± %8.5** | sıfırdan **1.45σ** (n=2'de 1.6σ) — **hâlâ sıfırdan ayırt edilemiyor** |
| **Q_E, log** | %79.4 ± %7.0 zarf | **%70.2 ± %7.1 zarf** | merkez **−9.1 puan** (eski se'nin 1.3'ü), çubuk **×1.02 — sıkışmadı** |
| **Q_E, log kilit** | %20.6 ± %7.0 | **%29.8 ± %7.1** | sıfırdan **4.20σ** |
| **M, doğrusal** | %133.7 ± %18.7 | **%144.0 ± %13.9** | çubuk ×0.74 |
| **Q_E, doğrusal** | %26.0 ± %2.3 | **%23.0 ± %2.3** | çubuk ×1.02 |

**Yalnız çubuk sıkışmadı; merkez de kaydı** — ve iki HAM nicelik
**ters yönlere** kaydı: `M` yukarı (zarf lehine), `Q_E` aşağı (kilit
lehine).

`Q_E`'nin çubuğunun sıkışmaması ölçülmüş bir olgudur: `sd_VS(Q_E, log)`
n = 2'de 0.003543 idi, n = 4'te **0.005174** (×1.46). Aynı şey `g_E`'de
daha çarpıcı: n = 2'nin iki tohumu neredeyse özdeş `g_E` vermişti
(sd = 0.000016), n = 4'te sd gerçekçi bir değere çıktı (0.000601,
**×38.8**) ve çubuk **×1.40 büyüdü**. **İki tohumluk se, bir yön garantisi değildi.**

### 3.2 "%80/%20" resmi ayakta mı? — **Ayakta ama gevşedi**

180'in gücü, "iki bağımsız HAM niceliğin aynı şeyi söylemesi"ydi
(%81.5 ↔ %79.4). n = 4'te ikisi **80'in iki yanına ayrıldı**:
`M` **%87.6 ± %8.5**, `Q_E` **%70.2 ± %7.1**.

* Aralarındaki kilit payı farkı: **0.174 ± 0.111 ⇒ 1.57σ** —
  istatistiksel olarak henüz çelişki değil, ama n = 2'deki
  yakınsamanın **tesadüf payı olduğunu** gösteriyor.
* Her ikisi de %80/%20'yi kendi çubuğu içinde barındırıyor
  (M: 80 merkezden 0.9σ; Q_E: 80 merkezden 1.4σ).
* **Ama kilidin sıfırdan farklı olup olmadığı para birimine bağlı:**
  `M`'de kilit payı %12.4 ± %8.5 (**1.45σ — sıfırla uyumlu**),
  `Q_E`'de %29.8 ± %7.1 (**4.20σ — sıfırdan uzak**). Aynı defterde
  iki HAM nicelik "kilit var mı yok mu"da ayrışıyor.

**Değişmeyen:** hiçbir okuma "%91.4 kilit"i desteklemiyor (en yüksek
log kilit payı %29.8 ± %7.1). 180 HÜKÜM(ii) **ayakta**.
**Değişmeyen:** konvansiyon uçurumu (M: log %87.6 ↔ doğrusal %144.0;
Q_E: log %70.2 ↔ doğrusal %23.0) — ikisi de basılıyor, aralarında
seçim yapılmıyor.

### 3.3 Okunur olan yeni satırlar ne söylüyor? — **Yeni okunur satır YOK**

Görevin sorduğu iki satırın açık yanıtı:

* **Δθ:** hâlâ **HÜKÜMSÜZ** (`|Z|/se` = 0.87 log, 0.88 doğrusal).
  Merkez ciddi kaydı — zarf payı −0.007 → **+0.559** — yani iki
  tohumdaki "θ tamamen kilittir" görüntüsü **tohum gürültüsüydü**;
  ama yeni merkez de okunur değil. n = 4 defteri θ hakkında **hiçbir
  şey söylemiyor** ve ön-kayıt bunun **VS tohumu eklenerek
  düzeltilemeyeceğini** (|ZARF| < 2·se_VF) sayılardan önce yazmıştı.
* **ΔQ_X:** hâlâ **HÜKÜMSÜZ** ve n = 2'dekinden **daha kötü**
  (`|Z|/se` 0.79 → 0.41 log). ZARF merkezi −0.076 → −0.027'ye küçüldü,
  yani "Q_X'te büyük bir zarf etkisi var" görüntüsü de **tohum
  gürültüsüydü**. Zarf/kilit bölünmesi **okunamaz**; süslenmiyor.

Okunur kalan satırların hepsi n = 2'dekilerle **aynı satırlardır**:
`M`, `Q_E`, `ρ_E`, `g_E`. `ρ_E`'nin işaret dönmesi n = 4'te
**güçlendi** (log pay −1.82 ± 0.32, kilit payı +2.82 ± 0.32, 8.95σ):
180a'nın işaret kuralı maddesi yine tam olarak bunu öngörmüştü;
kaydediliyor, kurtarılmıyor.

---

# HÜKÜM

## (i) İNŞA: dört tohumun dördü de temiz

V1–V6, VS3 ve VS4'te de geçti; zarf dört tohumda **bit-bit aynı**
(tek sha), `L` dördünde de aynı. n_VS = 4 defteri **kapı ihlali
olmadan** kuruldu.

## (ii) MANŞET: **%87.6 ± %8.5 zarf / %12.4 ± %8.5 kilit** (M, log)

180'in `%81.5 ± %11.3`'ü ile karşılaştırıldığında: **çubuk ×0.75
sıkıştı, merkez +6.1 puan kaydı** (eski çubuğun yarısı kadar).
`Q_E` (log) ise ters yöne kaydı ve çubuğu **sıkışmadı**:
**%70.2 ± %7.1 zarf / %29.8 ± %7.1 kilit**. Doğrusal konvansiyonda
`M` zarf payı %144.0 ± %13.9, `Q_E` kilit payı %77.0 ± %2.3.

## (iii) "%80/%20" AYAKTA, AMA n = 2'DEKİ MUTABAKAT TESADÜFTÜ

İki HAM niceliğin kilit payları 1.57σ ayrıldı (%12.4 ↔ %29.8) ve
"kilit sıfırdan farklı mı" sorusuna **farklı yanıt veriyorlar**
(1.45σ ↔ 4.20σ). Her ikisi de %80/%20'yi çubuğu içinde tutuyor;
hiçbiri "%91.4 kilit"e yaklaşmıyor.

## (iv) BOĞULAN SATIRLARIN HİÇBİRİ KURTULMADI

`Q_X`, `ρ_X`, `θ`, `g_X`, `g_cal` — **beşi de, iki konvansiyonda da
HÜKÜMSÜZ**. Bunlardan `Q_X`, `ρ_X`, `θ` için ön-kayıt bunu
**yapısal** gerekçeyle (|ZARF| < 2·se_VF) sayılardan önce yazmıştı:
onları boğan VS'nin tohum sayısı değil, **VF ailesinin saçılımı**.
`g_cal` için ön-kayıt "okunur olacak" demişti; **tutmadı** ve
tutmadığı burada yazılıdır.

## (v) İKİ TOHUMLUK HATA ÇUBUĞU GÜVENİLMEZDİ — İKİ YÖNDE DE

Çubuklar öngörüldüğü gibi ×0.71 civarı sıkışmadı: `M`'de ×0.75,
`Q_E`'de **×1.02**, `g_E`'de **×1.40**. Sebep ölçülü: `sd_VS(n=2)`
kimi nicelikte gerçek saçılımı **eksik** kestiriyordu (`sd_VS`: `Q_E` ×1.46,
`g_E` **×38.8**). Bu, 180'in sayılarına yönelik bir düzeltme
değil, **çözünürlüğün kendisinin ölçümüdür**.

## (vi) KAPSAM DIŞI KALAN, DEĞİŞMEDİ

Kırpma hakemi (K2 / `169_k2b` / ρ₃ / C / z_G / VF ailesi) **bu
görevde ele alınmadı, hiçbir sayısı yeniden hesaplanmadı**; 180'in
H-180a **HÜKÜMSÜZ** hükmü ve "DAL B öldü" alt-bulgusu aynen durur.
`W_pos` borcu (179 ŞART ③) 181'de de ele alınmadı.

## (vii) ÖLÇÜLMÜŞ SINIR — 180'inki, n = 4'te keskinleşerek duruyor

`r(τ)` kısmen kilit taşıyor: ölçülen ⟨R_bant⟩_VS/⟨R_bant⟩_VF oranı
τ boyunca **düz** (0.941–0.958), tarif **eğimli** (0.892–0.946), fark
τ ile **+%1.3'ten +%5.5'e** büyüyor. Zarf payının **üst sınır tadı**
n = 4'te de sürüyor; yönü belli, büyüklüğü hâlâ ölçülmedi.

---

# 🔒 MÜHÜR

> **Dört tohumlu VS ailesiyle, gerçeğin fazlasının çapasız zarf payı
> log konvansiyonunda `M`'de %87.6 ± %8.5 (kilit %12.4 ± %8.5, sıfırdan
> 1.45σ), `Q_E`'de %70.2 ± %7.1 (kilit %29.8 ± %7.1, sıfırdan 4.20σ);
> doğrusal konvansiyonda `M` zarf payı %144.0 ± %13.9'dur. n = 2'den
> n = 4'e geçişte yalnız çubuk sıkışmadı, merkez de kaydı ve iki HAM
> nicelik ters yönlere gitti; %80/%20 resmi her iki çubuğun içinde
> ayakta, ama n = 2'deki birebir mutabakat tesadüftü. Boğulan beş
> satırın hiçbiri okunur olmadı: `Q_X`, `ρ_X` ve `θ` için ön-kayıt
> bunu YAPISAL gerekçeyle önceden yazmıştı (onları VF'nin saçılımı
> boğuyor, VS'nin tohum sayısı değil); `g_cal` için ön-kayıt okunur
> olacağını yazmıştı ve TUTMADI.**

**Ne mühürlenmedi:** `Δθ` ve `ΔQ_X`'in zarf/kilit bölünmesi — n = 4'te
de okunamıyor ve VS tohumu eklemekle okunamayacağı ölçülü;
"%91.4 kilit" — bu defterde de hiçbir konvansiyonda yeniden
üretilemedi.

---

## Sıradaki adım (bu ölçümün işaret ettiği)

1. **VF ailesini büyütmek, VS'yi değil.** 181'in en sağlam bulgusu:
   `Q_X`, `ρ_X`, `θ` satırlarını boğan **VF saçılımıdır**
   (|ZARF| < 2·se_VF). VS5–VS8 bu satırları **ilkesel olarak**
   açamaz; VF5–VF8 açabilir. Aynı fiyata (~16 dk/tohum) doğru yere.
2. **`M` ↔ `Q_E` ayrışmasının kaynağı.** İki HAM nicelik kilit payında
   1.57σ ayrıldı ve "kilit sıfır mı" sorusuna farklı yanıt veriyor;
   bu, tek bir "kilit payı" sayısının var olup olmadığını sınayan
   doğrudan bir sorudur.
3. **E bacağının kırpma defteri** (180'in 2. maddesi) — hâlâ borç,
   yeni inşa gerektirmiyor.
4. **W_pos borcu** (179 ŞART ③) — hâlâ ele alınmadı.
