# 177 — ÜÇÜNCÜ TOHUM HAKEMLİĞİ
### (4 Eylül 2026, Opus tayfası — tek soruluk sefer)

**Soru (tek).** 176'nın `F3`'ü, θ'nın kilit-duyarlılığını **HÜKÜMSÜZ**
bıraktı: iki vekil tohumun θ'ları `0.8666712 / 0.8321358`, ortalama
`0.8494035`; dal (a) eşiğine uzaklık `0.0162134`, F7 saçılımı
`0.0345354` — gürültü marjı örtüyordu. 177 **yalnızca tohum sayısını
artırır** (2 → 3, gerekirse 4) ve aynı soruyu yeniden sorar:
**θ'nın kilit-duyarlılığı YAŞADI mı, ÖLDÜ mü?**

Kapsam dar: `R_η`, `g_E`, `M`, DC kaçağı, H-F2 — hepsi 176'da hükme
bağlandı; 177 onları **yeniden açmaz**. Git'e **dokunulmadı**.

| betik | kapı | koşu | ham çıktı |
|---|---|---|---|
| `177a_onkayit.py` | **K2 ön-kaydı** (kural donduruldu) | 0.1 s | `177/ONKAYIT_177.json` |
| `177b_guc_geometrisi.py` | ön-mühürlü: kuralın gücü (**ölçüm YOK**) | 3 s | `177/GUC_GEOMETRISI.json` |
| `176_configs/176b_vekil_insa.py` | **K1 (inşa)** — *değiştirilmeden çağrıldı* | ~10 dk/gaz | `176/z_VF*.npy`, `176/insa_VF*.json` |
| `176_configs/176c_olcum.py` | **K1 (ölçüm)** — *değiştirilmeden çağrıldı* | ~5–9 dk/gaz | `167/C_VF*.json`, `176/G_VF*.json` |
| `177c_hukum.py` | **K3 hüküm** — dondurulmuş kural + θ satırı | 1 s | `177/HUKUM_177_n3.json`, `…_n4.json` |
| `177d_pencere.py` | *ölçüm sonrası*: `n = 4` penceresi (θ₄ yokken) | 2 s | `177/PENCERE_n4.json` |
| `177e_tani.py` | *ölçüm sonrası*: gürültünün kaynağı + maliyet | 1 s | `177/TANI.json` |

`177a`, `177b`, `177c` **ön-mühürlüdür** (`177c`, VF3 kurulurken ve
hiçbir 177 θ'sı yokken yazıldı; sha `1e5cdd1a…`). `177d` ve `177e`
ölçümden sonra yazılmıştır ve **hiçbir eşiğe dokunmaz**.

> **Makine sadakati.** 177 inşa ya da ölçüm tarafında **tek satır yeni
> kod yazmadı**: 176'nın `176b_vekil_insa.py` ve `176c_olcum.py`
> dosyaları **aynı dosyalar olarak, alt-süreç halinde** çağrıldı. Bit-bit
> aynılık iddiası bu yüzden dosya kimliğine dayanır, kopyaya değil.

Önbellek kökü `/private/tmp/claude-501/.../scratchpad/`; vekil ham
çıktıları makinenin kendi yerine (`176/`, `167/`, `155/`) düşer, 177'nin
hüküm çıktıları `177/`'dedir.

---

> ## TEK CÜMLELİK HÜKÜM
>
> **H-F1b KALICI HÜKÜMSÜZ — ve 176'nın tek ayakta kalan faturası
> düştü.** Üçüncü ve dördüncü vekil tohumlar (zarf `Hkeskin`'inkiyle
> `sha256` özdeş, kapılar `G1–G5` tam) θ'yı **`0.8797` ve `0.8917`**
> verdi; dört-tohumlu ortalama `0.8494 → 0.8676` yükselerek dal (a)
> eşiğinin (`0.8656170`) **üstüne** çıktı ⇒ o dal artık **hiçbir tohum
> sayısıyla yaşayamaz**; dal (b) hâlâ dışlanıyor ama **kıl payıyla**
> (`0.0257759` vs `0.0257336`). Ödeme oranı
> **`ω_θ`: +0.716 (n=2) → +0.932 (n=3) → +1.2235 (n=4)** — ön-kayıtlı
> `(0,1]` bandının **dışında**, ve dışında olduğu da kesin değil; hatta
> `n = 3`'te KESİN olan **işaret** sınavı `n = 4`'te kıl payı kaçtı
> (`0.0296119 < 0.0299223`). **4/4 tohumda θ ikizin altındadır** ve
> kilidin θ'da taşıdığı toplam gerçeğin fazlasının **%94.0**'üdür — yani
> duyarlılığın **işareti** her tohumda aynı, **büyüklüğü** ölçülemez
> (tohum başına %5.8 – %225.2). Sebep ölçüldü ve cebridir:
> `θ ≡ M/(g_E g_X²)` ve tohum varyansının **%92.4'ü `g_X`'ten**,
> **%0.1'i `g_E`'den** gelir. ⇒ ΔM defterinin **%91.4 kilit** satırının
> sayıları değişmedi (vekilden bağımsızdırlar), ama onu ayakta tutan
> **tek nedensel dayanak kalmadı**. Eşik gevşetilmedi, beşinci tohum
> koşulmadı.

---

## 0. ÖN-KAYIT (K2) — 4 Eylül 2026, 19:42:56 +03

`177a_onkayit.py`, **`tohum = 3` ile hiçbir gaz kurulmadan önce** koştu;
kendi sha256'sını

```
03a9b9addedf738b4f4b03f0dedd0d54848fd11764cb4100e9490630cd12aebe
```

`177/ONKAYIT_177.json`'a yazdı; dosya varlık kontrolüyle bir daha
yazılmaz. 176'nın ön-kaydı (`18:53:39`, sha `6173943f…`) referans olarak
içine gömüldü.

### 0.1 DEĞİŞMEYEN — 176'dan bit-bit alındı

* `F3` eşikleri: **dal (a) `0.8656169530534596`** (H-F1b YAŞAR),
  **dal (b) `0.8933338132966645`** (H-F1b ÖLÜR).
* θ çapaları: `Hkeskin 0.8933338132966645`, `son 0.9219381611740145`,
  `HA4 0.7235729225877114`; `Δlog θ(son←Hk) = +0.031517828888066625`,
  `Δlog θ(erfc←Hk) = −0.2107589897030476`.
* `F9`'un ödeme-oranı bandı: **`ω_j ∈ (0,1]`**.
* İnşa kapıları `G1–G5` (176b'nin kendi kodunda; **dokunulmadı**).
* Vekil reçetesi: zarf birebir `Hkeskin`, `φ_q ~ U(0,2π)`,
  `164_insa.coz_sadakatli`, `h = 0.015`, `nz = 300000`, `c = −½`.
* Ölçüm: `167_olcum.kos(ad, 0.40, 0.95, kule=0, düz=0)`;
  **θ := KALİB_u2 / g_cal (lo = 0.60)** — başka kestirici yok.

### 0.2 DEĞİŞEN — **yalnız tohum sayısı**

`n = 2 → 3` (VF3, `tohum = 3`), gerekirse `→ 4` (VF4, `tohum = 4`).
Tohumlar 176'nın `1, 2` dizisinin doğal devamıdır ve ölçümden **önce**
donduruldu (tohum alışverişi yok). **176'nın hiçbir eşiği
değiştirilmedi.**

### 0.3 `F7`'NİN n-TOHUMLU BİÇİMİ — eşik gevşetme değil, **özdeş talep**

176'nın `F7`'si iki tohumda `|y₁−y₂| > |ȳ − eşik| ⇒ HÜKÜMSÜZ` der.
`n = 2` için, `s` örnek standart sapması (`ddof = 1`) olmak üzere

```
s = |y₁−y₂|/√2 ,   SEM = s/√n = |y₁−y₂|/2   ⇒   |y₁−y₂| ≡ 2·SEM
```

Yani `F7`'nin literal metni, **"ortalamanın eşiğe uzaklığı iki standart
hatadan büyük olmalı"** talebinin ta kendisidir. 177 bu talebi AYNEN
korur ve yalnız daha çok tohumla değerlendirir:

```
SAÇ_n := 2·std(y, ddof=1)/√n        HÜKÜMSÜZ ⟺ SAÇ_n > |ȳ − eşik|
```

**Ön-kayıtın kendi içindeki sayısal sağlaması** (176'nın iki tohumuyla):

| | değer |
|---|---|
| 176'nın literal metni `\|y₁−y₂\|` | `0.0345353851862731` |
| 177'nin biçimi `2·s/√n` | `0.0345353851862731` |
| katı okuma `max−min` | `0.0345353851862731` |
| **\|fark\|** | **`0.000e+00` — ÖZDEŞ** |

Betik bu özdeşlik tutmazsa ön-kaydı **yazmayı reddeder** (koşuldu,
tuttu). Ayrıca bağlayıcı olmayan **katı okuma** `SAÇ^kat_n := max−min`
her hükümde raporlanır; ikisi ayrışırsa rapor bunu açıkça yazar, hüküm
`SAÇ_n`'e göre verilir.

### 0.4 HÜKÜM MERDİVENİ (dondurulmuş)

| adım | hüküm |
|---|---|
| **H1 (F3 / H-F1b)** | `ȳ ≤ 0.8656169530534596` **ve** `\|ȳ−eşik\| ≥ SAÇ_n` ⇒ **YAŞADI**; `ȳ ≥ 0.8933338132966645` **ve** marj ≥ SAÇ_n ⇒ **ÖLDÜ**; aksi ⇒ `n ← n+1` |
| **H2 (F9 / ω_θ)** | `Δ_i := log θ(Hk) − log θ(VF_i)`, `Δ̄ = ort_i Δ_i`, `ω_θ = kilit_θ/Δ̄`; **`ω_θ ∈ (0,1] ⟺ Δ̄ ≥ kilit_θ`** (kilit_θ > 0 olduğu için özdeş); `\|Δ̄ − kilit_θ\| ≥ SAÇ_n(Δ)` ise KESİN |
| **üst sınır** | `n = 4`. Dörtte de kesinleşmezse **KALICI HÜKÜMSÜZ** yazılır. **Eşik gevşetme YOK.** |

`kilit_θ` elle girilmedi, ön-kayıtta çapalardan türetildi (176'nın F9
okumasının aynısı):

```
f_θ     = Δ(son←Hk)/Δ(erfc←Hk) = −0.1495444
kesim_θ = f_θ·Δ(son←Hk)        = −0.0047133
kilit_θ = (1−f_θ)·Δ(son←Hk)    = +0.0362311      ⇒  ω_θ ∈ (0,1] ⟺ Δ̄ ≥ 0.0362311
```

### 0.5 ÖLÜM MADDELERİ (ölçümden önce yazıldı)

* **İnşa kapısı ölümü.** Bir tohum `G1–G5`'ten birini tutturamazsa gaz
  **ölçülmez**, ortalamaya **girmez**, rapora yazılır; merdiven bir
  sonraki tohum indeksine geçer (en çok **2** kapı ölümü).
* **`F0`'ın `R_bant` alt-kapısı.** 176'da ıskalamıştı (`0.9125 / 0.8925`
  < `0.98`) ve kurtarılmadı. 177'de yeni tohum için de **raporlanır** ve
  **hiçbir hükme dayanak yapılmaz** — ne lehte ne aleyhte. Bu, ön-kayıt
  anında böyle yazılmıştır.

### 0.6 DÜRÜSTLÜK BEYANI (ön-kayıtta aynen duruyor)

1. **KÖRLÜK İDDİASI YOKTUR.** 176'nın bütün sayıları (θ: `0.8666712 /
   0.8321358`, `ω_θ = +0.716`, tohum başına `1.196 / 0.511`, saçılım
   `0.0345354`) bu tayfa tarafından okunmuştur. Ön-kayıt edilen şey
   **KURALDIR**.
2. **ÖN-KAYIT ANINDA ÖLÇÜLMEMİŞ OLAN:** `tohum = 3` (ve `4`) ile kurulan
   vekil gazın **hiçbir niceliği**. Böyle bir gaz kurulmamıştı.
3. **Ön-kayıt sonrası hiçbir eşik değiştirilmez.** Ölümler kurtarmasız.
4. Git'e dokunulmaz.

### 0.7 KURALIN GÜÇ GEOMETRİSİ — **ölçümden önce** (`177b`)

`177b_guc_geometrisi.py`, `tohum = 3` vekilinin hiçbir niceliği var
olmadan (VF3 inşası sürerken, `176c` hiç koşmamışken) koştu. **Yeni
ölçüm yoktur**; yalnız dondurulmuş kuralın cebri sorulmuştur: *176'nın
iki θ'sı sabitken üçüncü tohumun hangi değerleri hükmü KESİN kılar?*

| tarama | sonuç |
|---|---|
| `n = 3`, θ₃ ∈ [0.55, 0.98] — dal (a) KESİN kümesi | **BOŞ** |
| `n = 3` — dal (b) KESİN kümesi | **BOŞ** |
| `n = 3` — `H2` (ω_θ ≤ 1) KESİN kümesi | **BOŞ** |
| `n = 4`, θ₃ = θ₄ = c kesiti — dal (a) KESİN | `c ∈ [0.6510, 0.85325]` |
| `n = 4`, aynı kesit — `H2` KESİN | `c ∈ [0.71795, 0.84510]` |
| `n = 4`, (θ₃,θ₄) ∈ [0.70,0.95]² ızgarası | dal (a) **%26.0**, `H2` **%14.3** kesin |

> **ÖLÇÜMDEN ÖNCE KAYDA GEÇEN.** Dondurulmuş kural, iki mevcut tohumun
> saçılımı verildiğinde **hiçbir θ₃ değeri için `n = 3`'te
> kesinleşemez** — bu bir veri bulgusu değil, kuralın cebridir
> (`ȳ` eşiğe `Δ/3` hızıyla yaklaşırken `s` daha hızlı büyür).
> Dolayısıyla ön-kayıtlı merdiven **zorunlu olarak `n = 4`'e gider**;
> dördüncü tohumu kurma kararı **veriye değil kurala** bağlıdır ve bu
> satır θ₃ doğmadan yazılmıştır. `n = 4`'te kural gerçek bir güce
> sahiptir (kesin bölge boş değil), yani sefer "kurgulanmış hükümsüzlük"
> değildir.

---

## 1. K1 (İNŞA) — ÜÇÜNCÜ (ve DÖRDÜNCÜ) TOHUM

`176_configs/176b_vekil_insa.py` **değiştirilmeden**, alt-süreç olarak
`VF3 3` argümanlarıyla çağrıldı. Zarf yine `Hkeskin`'inkidir; değişen
tek şey `φ_q` dizisidir.

### 1a. İNŞA KAPILARI (`F0` / `G1–G5`) — **HEPSİ GEÇTİ**

| kapı | eşik | `VF1` (t=1) | `VF2` (t=2) | **`VF3` (t=3)** |
|---|---|---|---|---|
| G1 zarf farkı | `= 0.0` | 0.0e+00 ✓ | 0.0e+00 ✓ | **0.0e+00** ✓ |
| G1 `sha256(A)` | = Hkeskin | ✓ | ✓ | **`45af6fa5…9e43c267`** ✓ |
| G2 `φ≡0` sağlaması | `= 0.0` | 0.0 / 0.0 ✓ | 0.0 / 0.0 ✓ | **0.0 / 0.0** ✓ |
| G3 ilk-kök hücre | 300000/300000 | ✓ | ✓ | **300000/300000** ✓ |
| G4 `maks\|F\|` | ≤ 1e−8 | 1.863e−09 ✓ | 1.863e−09 ✓ | **1.863e−09** (aşan 0) ✓ |
| G5 sıralılık | TAM | ✓ | ✓ | **TAM**, min Δz = 0.091732 ✓ |
| ızgara nokta sayısı | — | 10 449 281 | 10 449 281 | **10 449 281** |
| `Σ\|A\|` / `ΣA²` | — | 27.815371002 / 0.292663602 | aynı | **aynı** |
| `rms S` / `rms S'` | — | 0.3825334 / 1.8939673 | aynı | **aynı** |
| `σ_ds` (z'den) | — | 0.44159 | 0.44273 | **0.44234** |
| `L` (z'den) | — | 12.029593224 | 12.029593224 | **12.029593224** |
| ikiye-bölme adımı | — | 26 | 23 | **21** |
| süre | — | 10.0 dk | 10.6 dk | **9.9 dk** |

> Üçüncü tohum, ilk ikisiyle **aynı zarfı** (bit-bit, aynı `sha256`),
> aynı ızgarayı ve aynı sadakat ölçütünü taşır. **Kapı ölümü yok.**

### 1b. Dördüncü tohum — kararın gerekçesi

`177b` (§0.7) **θ₃ doğmadan önce** kanıtladı ki dondurulmuş kural
`n = 3`'te hiçbir θ₃ için kesinleşemez. Ön-kayıtlı merdiven bu durumda
`n ← 4` der. Bu yüzden `VF4` (`tohum = 4`) **veriye bakılmadan**,
kuralın gereği olarak sıraya alındı: `176c VF3 → 176b VF4 → 176c VF4`.

### 1c. `VF4` (tohum 4) — inşa kapıları da **TAM**

| kapı | `VF4` |
|---|---|
| G1 zarf farkı / `sha256(A)` | **0.0e+00** / `45af6fa5…9e43c267` (= Hkeskin) ✓ |
| G2 `φ≡0` (S, S') | **0.0e+00 / 0.0e+00** ✓ |
| G3 ilk-kök hücre | **300000 / 300000** ✓ |
| G4 `maks\|F\|` | **1.863e−09** (aşan 0) ✓ |
| G5 sıralılık | **TAM**, min Δz = 0.083139 ✓ |
| ızgara / ikiye-bölme / süre | 10 449 281 / 25 / **9.8 dk** |
| `σ_ds` / `L` | 0.44108 / 12.029593224 |

> **Dört tohumun dördünde de `G1–G5` tam geçti; kapı ölümü yok.**
> Dört gazın genlik dizisi `sha256` düzeyinde **aynı**; aralarındaki tek
> fark `φ_q`'dur.

---

## 2. K1 (ÖLÇÜM) — ÜÇÜNCÜ TOHUMUN DEFTERİ

`176_configs/176c_olcum.py` **değiştirilmeden** çağrıldı
(`167_olcum.kos(VF3, 0.40, 0.95, kule=0, düz=0)`). `nline = 8981`,
özdeşlik denetimi `|ΔK|/K`: E `2.55e−15`, X `3.79e−15`.

| gaz | `θ` | `g_E` | `g_X` | `M` | `g_cal` | `Q_E` | `ρ_E` | `μ̂²_E` | `R_bant` min |
|---|---|---|---|---|---|---|---|---|---|
| `Hkeskin` (ikiz) | 0.893334 | 0.583701 | 0.704213 | 0.258590 | — | 0.91558 | 2.19932 | 0.056843 | 1.2857 |
| `son` (gerçek) | 0.921938 | 0.598206 | 0.704733 | 0.273907 | — | 0.87575 | 2.17960 | 0.040498 | 1.1986 |
| `VF1` (t=1) | 0.866671 | 0.736939 | 0.813474 | 0.422643 | 0.487663 | 0.29743 | 1.13064 | 0.007998 | 0.9125 |
| `VF2` (t=2) | 0.832136 | 0.736692 | 0.833331 | 0.425711 | 0.511588 | 0.29711 | 1.12837 | 0.007858 | 0.8925 |
| **`VF3` (t=3)** | **0.879713** | **0.735928** | **0.815483** | **0.430532** | 0.489401 | **0.29944** | **1.13393** | **0.007838** | **0.9182** |

> **Üçüncü tohum, 176'nın ÖLÜ maddelerini teyit etti** (yeniden
> açılmadı, yalnız kayda geçer): `g_E = 0.7359` — üç tohumda da
> `Hkeskin`'in `0.5837`'sinin **çok üstünde** (`F2` ölümü sağlam);
> `Q_E ≈ 0.2974` ve `μ̂²_E ≈ 0.0079` tohumdan tohuma **%1–3 içinde**
> aynı (kilit çöküşü tohuma duyarsız); `R_bant` min `0.9182` — yine
> `0.98`'in altında (ön-kayıt gereği **hiçbir hükme dayanak değil**).
> **Tohumdan tohuma ciddi biçimde oynayan tek nicelik θ'dır** (§5).

### 2a. `n = 3` HÜKMÜ (`177c`, dondurulmuş kural)

```
θ = 0.8666712 / 0.8321358 / 0.8797129     ȳ = 0.8595067
SAÇ_3 = 2s/√3 = 0.0283876                 [katı max−min = 0.0475771]

dal (a) ≤ 0.8656170 : taraf ✓  marj 0.0061103 < SAÇ_3       ⇒ HÜKÜMSÜZ
dal (b) ≥ 0.8933338 : taraf ✗  marj 0.0338272 > SAÇ_3 0.0283876
                       ⇒ (b) KESİN olarak DIŞLANDI (rahatlıkla)
⇒ H1 = H-F1b HÜKÜMSÜZ    (merdiven: n ← 4)

Δ_i = +0.0303006 / +0.0709646 / +0.0153647   Δ̄ = +0.0388767
SAÇ_3(Δ) = 0.0332265
ω_θ = 0.036231/Δ̄ = +0.9320   (tohum başına 1.1957 / 0.5106 / 2.3581)
  alt kenar (ω_θ > 0 ⟺ Δ̄ > 0):  marj 0.0388767 > 0.0332265  ⇒ **KESİN**
  üst kenar (ω_θ ≤ 1 ⟺ Δ̄ ≥ 0.0362311): marj 0.0026455 < SAÇ  ⇒ HÜKÜMSÜZ
⇒ H2 = BANDIN İÇİNDE ama kesin değil
```

Üçüncü tohum θ'yı **yukarı** çekti (0.8797 — ilk ikisinin de üstünde):
ortalama `0.8494 → 0.8595`, eşiğe marj `0.0162 → 0.0061`; `SAÇ` düştü
(`0.0345 → 0.0284`) ama marj daha hızlı eridi.

### 2b. `n = 4` PENCERESİ — **θ₄ doğmadan hesaplandı** (`177d`)

`177d_pencere.py`, `VF4` inşası sürerken (θ₄ yokken) koştu ve
dondurulmuş kuralın geometrisini ölçülmüş üç θ ile yeniden sordu:

| tarama | sonuç |
|---|---|
| `n = 4`, θ₄ ∈ [0.40, 0.99] — dal (a) KESİN kümesi | **BOŞ** |
| `n = 4` — dal (b) KESİN kümesi | **BOŞ** |
| `n = 4` — `H2` (ω_θ ≤ 1) KESİN kümesi | **BOŞ** |

> **KAYDA GEÇEN (θ₄ doğmadan).** Üç ölçülmüş θ'nın saçılımı verildiğinde
> **hiçbir θ₄ değeri** `n = 4`'te H1'i (ya da H2'nin üst kenarını)
> kesinleştiremez. Ön-kayıtlı üst sınır `n = 4` olduğuna ve **eşik
> gevşetme yasak** olduğuna göre hüküm, VF4 ne verirse versin,
> **KALICI HÜKÜMSÜZ** olacaktır. VF4 yine de kuruldu ve ölçüldü: ön-kayıt
> onu emrediyor ve dört-tohumlu ortalama ile `ω_θ` **teslim edilecek
> sayılardır**.

---

## 3. K3 — DÖRT TOHUMLU HÜKÜM (`177c`, dondurulmuş kural)

`VF4`'ün defteri: `θ = 0.891712`, `g_E = 0.737365`, `g_X = 0.806809`,
`M = 0.428004`, `Q_E = 0.298609`, `μ̂²_E = 0.008099`, `R_bant` min
`0.9141`, `nline = 8981`, `|ΔK|/K` E `3.44e−15`.

### 3a. H1 — `F3` / H-F1b

```
θ = 0.8666712 / 0.8321358 / 0.8797129 / 0.8917115     ȳ = 0.8675579
SAÇ_4 = 2s/√4 = 0.0257336            [katı max−min = 0.0595757]

dal (a) ≤ 0.8656170 : taraf ✗ (ȳ eşiğin ÜSTÜNDE)  marj 0.0019409  ⇒ HÜKÜMSÜZ
dal (b) ≥ 0.8933338 : taraf ✗ (ȳ eşiğin ALTINDA)  marj 0.0257759  ⇒ (b) DIŞLANDI
                       0.0257759 > SAÇ_4 0.0257336 — **kıl payı** (4.2e−05)
⇒ H1 = **H-F1b HÜKÜMSÜZ**   |   merdiven: n = 4 üst sınır ⇒ **KALICI HÜKÜMSÜZ**
```

> **DÖRDÜNCÜ TOHUM YÖNÜ ÇEVİRDİ.** `θ(VF4) = 0.8917115`, ikizin
> `0.8933338`'ine **%0.18 uzaklıkta**. Dört-tohumlu ortalama
> `0.8675579`, dal (a) eşiğinin **üstüne** çıktı (2 tohumda 0.8494,
> 3 tohumda 0.8595, 4 tohumda 0.8676 — her yeni tohum ortalamayı eşiğe
> doğru itti). Ortalama artık eşiğin **yanlış tarafında** olduğu için
> dal (a) **hiçbir tohum sayısıyla YAŞAyamaz** (§5).
> Dal (b) hâlâ dışlanıyor ama `n = 3`'teki rahatlığını
> (`0.0338 > 0.0284`) yitirdi; `n = 4`'te marj saçılımı **kıl payı**
> geçiyor. Bu satır ağırlık taşımaz.

### 3b. H2 — `F9` / `ω_θ` bandı

```
Δ_i = log θ(Hk) − log θ(VF_i)
    = +0.0303006 / +0.0709646 / +0.0153647 / +0.0018176      (4/4 pozitif)
Δ̄ = +0.0296119     SAÇ_4(Δ) = 0.0299223     [katı 0.0691470]
kilit_θ = +0.0362311

ω_θ = kilit_θ/Δ̄ = **+1.2235**            (tohum başına 1.1957 / 0.5106 /
                                            2.3581 / 19.9331)
üst kenar (ω_θ ≤ 1 ⟺ Δ̄ ≥ 0.0362311): taraf ✗  marj 0.0066192 < SAÇ ⇒ HÜKÜMSÜZ
alt kenar (ω_θ > 0 ⟺ Δ̄ > 0)        : taraf ✓  marj 0.0296119 < SAÇ 0.0299223
                                        ⇒ HÜKÜMSÜZ (kıl payı KAÇTI, 3.1e−04)
⇒ H2 = **BANDIN DIŞINDA — ama kesin değil**
```

> **176'NIN `ω_θ = +0.716`'SI İKİ TOHUMUN ŞANSIYDI.** Dört tohumla nokta
> kestirimi **`+1.2235`**'e çıktı, yani ön-kayıtlı `(0,1]` bandının
> **dışına**. Ne bandın içinde olduğu ne dışında olduğu kesinleşiyor.
> Üstelik `n = 3`'te KESİN olan **işaret sınavı** (`ω_θ > 0`) `n = 4`'te
> **kıl payı kaçtı** (`0.0296119 < 0.0299223`) — dördüncü tohumun
> `Δ₄ = +0.0018` gibi neredeyse sıfır bir değer vermesi yüzünden.
> **4/4 tohumda `Δ_i > 0`**'dır (yani her tohumda θ ikizin altına
> düşüyor); çöken şey işaretin kendisi değil, ortalamanın
> **2-standart-hata** sınavıdır.

---

## 4. NİHAİ AYRIŞIM TABLOSUNUN θ SATIRI — GÜNCEL (parametresiz)

`ΔlogM ≡ Δlog g_E + 2Δlog g_X + Δlog θ`. `KESİM_j` ve `KİLİT_j`
yalnız `son`/`Hkeskin`/`HA4` çapalarından çıkar; **vekilden bağımsızdır
ve 177'de değişmemiştir.** Değişen tek sütun `Δ(Hk←vekil)` ve ondan
türeyen `ω_j`'dir.

| çarpan | `Δ(son←Hk)` | `f_j` | **KESİM_j** | **KİLİT_j** | `Δ(Hk←vekil)` | **ω_j** | hüküm |
|---|---|---|---|---|---|---|---|
| `g_E` | +0.024547 | +0.3903 | +0.009581 | +0.014966 | −0.232948 *(176, n=2)* | **−0.064** | ✗ ÖLDÜ (176) |
| `g_X²` | +0.001477 | +0.0418 | +0.000062 | +0.001415 | −0.312583 *(176, n=2)* | **−0.005** | ✗ ÖLDÜ (176) |
| **`θ`** | **+0.031518** | **−0.1495** | **−0.004713** | **+0.036231** | **+0.029612** *(n = 4)* | **+1.2235** | **HÜKÜMSÜZ** |
| `M` | +0.057542 | −0.5114 | −0.029429 | +0.086970 | −0.494899 *(176, n=2)* | **−0.176** | ✗ ÖLDÜ (176) |
| `R_η` (doğrusal) | +0.023245 | +0.7445 | +0.017305 | +0.005940 | +0.456547 *(176, n=2)* | **+0.013** | ✓ (176) |

**θ satırının eski hâli (176, n = 2):** `Δ(Hk←vekil) = +0.050633`,
`ω_θ = +0.716` ✓ *(tek ödenen fatura)*.
**Yeni hâli (177, n = 4):** `Δ(Hk←vekil) = +0.029612`,
`ω_θ = +1.2235` **HÜKÜMSÜZ**.

Kilidin θ'da taşıdığı toplam / gerçeğin ikizden fazlası:
**%94.0** (tohum başına **%96.1 / %225.2 / %48.7 / %5.8**) — 176'daki
%160 yerine.

### 4a. ΔM'nin KESİM/KİLİT defteri — **sayıları değişmedi**

| | KESİM | KİLİT | toplam |
|---|---|---|---|
| E kanalı (`g_E`, `g_X²`) | +0.009643 | +0.016381 | +0.026024 |
| θ kanalı | −0.004713 | +0.036231 | +0.031518 |
| **TOPLAM (M)** | **+0.004929** | **+0.052612** | **+0.057542** |
| ΔM'nin yüzdesi | **+%8.6** | **+%91.4** | %100 |

Bu tablo `son`, `Hkeskin`, `HA4`'ün çapalarından çıkar; **hiçbir vekil
niceliği girmez**. 177 onu **rakamsal olarak değiştirmez** — ve
değiştiremez.

---

## 5. TANI — θ'nın tohum gürültüsü nereden geliyor? (`177e`, ölçüm sonrası)

Özdeşlik `θ ≡ M / (g_E · g_X²)` dört vekilde de **sayısal olarak tam**
(kalıntı `0.0e+00 / 0.0e+00 / 0.0e+00 / 1.1e−16`). Tohumlar arası
saçılım:

| nicelik | dört değer | CV | `Var(log θ)` payı |
|---|---|---|---|
| `M` | 0.422643 0.425711 0.430532 0.428004 | **%0.79** | **%7.5** |
| `g_E` | 0.736939 0.736692 0.735928 0.737365 | **%0.08** | **%0.1** |
| `g_X` | 0.813474 0.833331 0.815483 0.806809 | **%1.39** | **%92.4** (4·Var log g_X) |
| **`θ`** | 0.866671 0.832136 0.879713 0.891712 | **%2.97** | — |

> **BULGU (177-a).** **θ'nın tohum gürültüsünün %92'si `g_X`'ten
> gelir**; `g_E`'ninki **%0.1**'dir. Yani vekil gazlarda ölçülen "θ
> gürültüsü" θ'nın kendi üçüncü-moment içeriğinin gürültüsü değil,
> `θ = M/(g_E g_X²)` bölmesindeki **X-kanalı model kalitesinin**
> gürültüsüdür — ve `g_X²` olduğu için **dört katına** çıkarak girer.
> `g_E` ise tohumdan tohuma neredeyse sabittir (`0.7359–0.7374`), bu
> yüzden 176'nın `F2` ölümü dört tohumda da **sarsılmaz**.

### 5a. Kesinlik için gereken tohum sayısı (ölçülmüş `s` ile)

| sınav | ortalama | `s` | marj | gereken `n` |
|---|---|---|---|---|
| H1 dal (a): `θ̄ ≤ 0.8656170` | +0.8675579 | 0.0257336 | 0.0019409 | **∞** — ortalama eşiğin **yanlış tarafında** |
| H2 üst kenar: `Δ̄ ≥ 0.0362311` | +0.0296119 | 0.0299223 | 0.0066192 | **≳ 82** (≈ 20 saat inşa+ölçüm) |
| H2 alt kenar: `Δ̄ > 0` | +0.0296119 | 0.0299223 | 0.0296119 | **≳ 5** |

> Ön-kayıtlı tavan `n = 4`'tür ve **eşik gevşetme yasaktır**; beşinci
> tohum **koşulmadı**. Alt kenarın `n ≳ 5` ile kapanacak olması 178'in
> en ucuz adımıdır, 177'nin hükmü değildir.

---

## 6. HÜKÜM

```
F0 (G1–G5)   ✓✓✓✓  dört tohumda da TAM         (R_bant ✗, bağlayıcı değil)
H1 (F3)      dal (a) HÜKÜMSÜZ (ȳ eşiğin üstünde) ; dal (b) dışlandı (kıl payı)
H2 (ω_θ)     BANDIN DIŞINDA, kesin değil ; işaret sınavı da kıl payı kaçtı
merdiven     n = 4 üst sınır  ⇒  KALICI HÜKÜMSÜZ
```

> ## (i) H-F1b **KALICI HÜKÜMSÜZDÜR** — askıda değil, tavanda kapandı
>
> Ön-kayıtlı tavana (`n = 4`) kadar gidildi, eşiklere **dokunulmadı**.
> θ'nın kilit-duyarlılığı ne **YAŞADI** ne **ÖLDÜ**: dal (a)'nın literal
> eşiği (`θ̄ ≤ 0.8656170`) **kaçırıldı** — üstelik dört-tohumlu ortalama
> `0.8675579` eşiğin **üstünde** olduğu için bu dal artık hiçbir tohum
> sayısıyla yaşayamaz. Dal (b) (θ çökmez) dışlanmayı sürdürüyor ama
> `n = 4`'te **kıl payıyla** (`0.0257759` vs `0.0257336`).

> ## (ii) 176'NIN "TEK ÖDENEN FATURA: θ" OKUMASI **AYAKTA DEĞİL**
>
> `ω_θ`: **+0.716 (n = 2) → +0.932 (n = 3) → +1.2235 (n = 4)**.
> Nokta kestirimi ön-kayıtlı `(0,1]` bandının **dışına** çıktı; dışında
> olduğu da kesin değil. 176'nın "θ faturası ödenebilir" hükmü **iki
> tohumun şansıydı**: tohum başına `ω_θ` **0.51'den 19.93'e** kadar
> savruluyor.

> ## (iii) AMA KİLİT-DUYARLILIĞININ **İŞARETİ** DÖRT TOHUMDA DA AYNI
>
> `Δ_i = log θ(Hk) − log θ(VF_i) > 0` **4/4**. Her tohumda θ ikizin
> altına düşüyor (`0.8917`, `0.8797`, `0.8667`, `0.8321` — hepsi
> `0.8933`'ün altında). Kilidin θ'da taşıdığı toplam, gerçeğin ikizden
> fazlasının **%94.0**'ü (tohum başına %5.8 – %225.2). Çöken şey
> işaretin kendisi değil, **ortalamanın 2-standart-hata sınavıdır**
> (`n ≳ 5` yeterdi; tavan 4'tü).

> ## (iv) SEBEP ÖLÇÜLDÜ: GÜRÜLTÜ θ'NIN DEĞİL, `g_X`'İN
>
> `θ ≡ M/(g_E g_X²)` özdeşliği vekillerde tamdır ve tohum
> varyansının **%92.4'ü `g_X`'ten**, %7.5'i `M`'den, **%0.1'i
> `g_E`'den** gelir. Yani H-F1b'yi hükümsüz bırakan şey θ'nın fiziği
> değil, **böleninin ölçüm gürültüsüdür** — ve o bölen faz karıştırılmış
> gazda tek bantlık bir kestirimdir (`R_bant ≥ 0.98` filtresi vekilde
> yalnız en alt bandı geçiriyor; 176 §2c).

---

## 7. "%91.4 KİLİT" DEFTERİ DEĞİŞİYOR MU? — kısa değerlendirme

**Rakamlar: HAYIR. Gerekçe: EVET — ve kötüye.**

1. **Sayılar yerinde.** `KESİM = +0.004929 (%8.6)`,
   `KİLİT = +0.052612 (%91.4)` yalnız `son`, `Hkeskin`, `HA4`
   çapalarından çıkar; **hiçbir vekil niceliği girmez**. 177 bu
   tabloyu ne değiştirdi ne değiştirebilirdi.
2. **Ama tablo bir ATIFTIR, bir ölçüm değil.** 176'nın nedensel sınavı
   ona **tek bir dayanak** bırakmıştı: θ satırının `ω_θ = +0.716`
   ile ödenebilir olması. **177 o dayanağı kaldırdı** — `ω_θ` dört
   tohumda `+1.2235`, bandın dışında ve hükümsüz.
3. **Sonuç:** 177'den sonra ΔM'nin **hiçbir** çarpanında kilit faturası
   **kesin biçimde ödenebilir değildir**: `g_E` (−0.064), `g_X²`
   (−0.005), `M` (−0.176) işaret ters (176, kesin); `θ` hükümsüz.
   Kilit para biriminde ayakta kalan tek ödeme `R_η`'dır
   (`ω_R = +0.013`) — ve `R_η` ΔM defterinin çarpanı **değildir**.
   ⇒ **"Gerçeğin fazlasının %91.4'ü kilittir" cümlesi, 176'da olduğu
   gibi 177'de de bir defter satırıdır; nedensel sınavda karşılığı
   yoktur.** 175'in atfı çürütülmedi — **desteksiz** bırakıldı.
4. **Kaybedilmeyen.** Kilidin **varlığı** ve devasalığı (`R_η`
   1.265 → 0.808, `Q_E` −%68, `μ̂²_E` −%86, `m3_çizgi` −%77) 176'da
   ölçülmüştü ve 177 ona dokunmadı; üçüncü ve dördüncü tohum onu ayrıca
   teyit ediyor: `Q_E` dört tohumda **%0.8**, `μ̂²_E` **%3.3**, `g_E`
   **%0.2** içinde aynı. **Kilit ölçerleri tohuma duyarsızdır; duyarsız
   olmayan tek şey θ'nın böleni `g_X`'tir.**

---

## 8. DENETİM

| | sınav | sonuç |
|---|---|---|
| **D1** | makine aynı mı? | `176b`/`176c` **aynı dosyalar**, alt-süreç olarak çağrıldı; 177 inşa/ölçüm tarafında tek satır yazmadı |
| **D2** | zarf birebir mi? | `maks\|A−A_Hk\| = 0.0e+00`, `sha256` dört tohumda da `45af6fa5…9e43c267` (= `Hkeskin`) |
| **D3** | değerlendirici doğru mu? | `φ≡0`'da `maks\|ΔS\| = maks\|ΔS'\| = 0.0e+00` (VF3, VF4) |
| **D4** | inşa kapıları | ilk-kök **300000/300000**, `maks\|F\| = 1.863e−09` (aşan 0), sıralılık **TAM** — dört tohumda da |
| **D5** | ızgara aynı mı? | `ng = 10 449 281` dört tohumda ve `Hkeskin`'de aynı; `L = 12.0295932` |
| **D6** | ölçüm zinciri | `nline = 8981`, `\|ΔK\|/K` ≤ `5.2e−15` (VF3, VF4) |
| **D7** | çarpan özdeşliği | `θ − M/(g_E g_X²)` = `0.0 / 0.0 / 0.0 / 1.1e−16` |
| **D8** | ön-kayıt sonradan değişti mi? | `177/ONKAYIT_177.json` bir kez yazıldı (varlık kontrolü); sha `03a9b9ad…` §0'da |
| **D9** | `F7` genellemesi eşik gevşetiyor mu? | ön-kaydın kendi sağlaması: `n = 2`'de `2s/√n − \|y₁−y₂\| = 0.0e+00` |
| **D10** | eşik değişti mi? | **hayır** — 176'nın `F3(a)`, `F3(b)`, `ω ∈ (0,1]` eşikleri aynen kullanıldı |
| **D11** | tohum alışverişi | yok: `tohum = 3`, `4` ön-kayıtta donduruldu; beşinci tohum **koşulmadı** |
| **D12** | git | **dokunulmadı** |

### 8a. Dürüstlük notları

* **`177b` hiçbir 177 ölçümü yokken, `177d` ise `θ₄` doğmadan koştu**;
  ikisi de "bu merdiven kesinleşemez" dedi. Yani **KALICI HÜKÜMSÜZ**
  sonucu, θ₃/θ₄'ün ne çıktığına bakılarak seçilmiş bir dal değildir.
* **`177c` ön-mühürlüdür:** VF3 kurulurken, hiçbir 177 θ'sı yokken
  yazıldı (sha `1e5cdd1a…`) ve `n = 3` ile `n = 4`'te **değiştirilmeden**
  koşuldu.
* **`177e` ölçümden sonra yazıldı** (`176g` gibi) ve hiçbir hükmü
  değiştirmez; yalnız kayda geçmiş sayıların aritmetiğidir.
* **Dal (b)'nin dışlanması `n = 4`'te kıl payıdır** (`4.2e−05`).
  Rapor bunu ağırlık taşımayan bir satır olarak yazar.
* **`F0`'ın `R_bant` alt-kapısı VF3/VF4'te de ıskaladı** (`0.9182`,
  `0.9141`) ve kurtarılmadı; ön-kayıt gereği hiçbir hükme dayanak
  yapılmadı.
* **Boşa giden koşu yok**; her koşulan gaz rapora girdi.
* `176d` (girişim) ve `176e` (DC kaçağı) VF3/VF4 için **koşulmadı**:
  177'nin sorusu yalnız θ'dır ve `R_η`/DC hükümleri 176'da kapandı.
  Bu yüzden yukarıdaki tablonun `R_η` satırı `n = 2` etiketlidir.

---

## Sıradaki adım (bu ölçümün işaret ettiği)

1. **`g_X`'i düzelt, θ kendiliğinden kesinleşir.** θ gürültüsünün
   **%92'si `g_X`'ten** geliyor ve `g_X` vekilde tek bantlık bir
   kestirim (çünkü `R_bant ≥ 0.98` filtresi vekilde yalnız en alt bandı
   geçiriyor — 176 §2c). Filtre yerine **bant-ağırlıklı** bir `g_X`
   (ya da θ'yı `M`'den `g_X²`'yi bölmeden okuyan bir kanal) mevcut dört
   tohumla H-F1b'yi **yeni koşu olmadan** karara bağlayabilir. En ucuz
   adım.
2. **Beşinci tohum işaret sınavını kapatır** (`n ≳ 5`), ama `ω_θ ≤ 1`
   için `n ≳ 82` gerekir (~20 saat). İşaret sınavı ucuz, band sınavı
   pahalıdır — ve 177'nin tavanı 4'tü.
3. **`ω_θ`'nın tohum saçılımı bir OLGUDUR** (0.51 → 19.93). Bu, `Δ_i`
   dağılımının ağır kuyruklu olduğunu söylüyor; ortalama yerine
   **medyan/işaret** temelli bir ödeme oranı ön-kaydı 178'in konusu
   olabilir — ama 177'nin eşikleri **geriye dönük gevşetilemez**.
4. **ΔM'yi `g_E` yerine `(Q_E, ρ_E)` çiftinde ayrıştır** (176'nın 2.
   maddesi hâlâ ayakta ve şimdi daha da çekici): `Q_E` dört tohumda
   **%0.8**, `μ̂²_E` **%3.3** içinde aynı — kilidin tekdüze ve
   tohuma-duyarsız ölçerleri bunlar; θ'nınki değil.
