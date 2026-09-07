# 180 — NEDENSEL DEFTER
### (7 Eylül 2026, Opus tayfası — `KALEM_NEDENSEL_DEFTER_07EYL2026.md`)
### 179 teftişinin "Not 5 öncesi ŞART" işi: çapasız ayrışım + çifte-sayım hakemi

> **Bu rapor tek başına okunur.** Her sayının yanında onu üreten betiğin
> adı vardır. Ölümler kurtarılmadan, belirsizlikler süslenmeden yazılır;
> hüküm verilemeyen yere **HÜKÜMSÜZ** yazılır.

---

## 0. K0 — DONMUŞ ÖN-KAYIT (`180a_onkayit.py`)

**Zaman (`date`):** `Mon Sep  7 22:11:42 +03 2026`
**Betik sha256:** `22af1b06c6d6cc26b5b19876646d99e3078d05af129e620d54363f06ee01126d`
**KALEM sha256:** `af823bbcb66b72f87edb8ad0f403b04fbe285b48e8b4d41b9cb4fc65a4b5c03c`
**Dosya:** `180/ONKAYIT_K0.json` (bir daha değiştirilmedi)

Bu betik **hiçbir yüzleşme niceliği hesaplamadan** koştu: ZARF payı,
KİLİT payı, bunların oranları, VS ailesinin herhangi bir defter satırı ve
VF ailesinin kırpma aktarımı **burada yoktur** — hepsi ön-kayıt diske
yazıldıktan sonra ölçüldü.

### 0.1 Aileler

| aile | tanım |
|---|---|
| `son` | gerçek ζ gazı (son 300000 sıfır) |
| `Hk` | `Hkeskin` — keskin ikiz: nominal merdiven `a_q = 1/(πk√q)`, τ≤1.00, **φ ≡ 0** |
| `VF` | **Hk-zarflı karışık** vekiller `VF1…VF4` (176/177'de kurulu; φ_q ~ U(0,2π), tohum 1–4) |
| `VS` | **son-zarflı karışık** vekiller `VS1, VS2` (180'de kuruldu; φ_q ~ U(0,2π), tohum 1–2) |

### 0.2 ZARF tanımı (dondurulmuş, parametresiz)

```
r(τ) := R_bant(son, τ_eff) / R_bant(Hkeskin, τ_eff)         [167 bant defteri]
A_q(VS) := a_q · r(τ_q) ,   τ_q = ω_q / L_hedef
r(τ): τ_eff ızgarasında DOĞRUSAL interpolasyon, dışarıda SABİT TUTMA (clamp)
```

Ölçülen profil (`180a`):

| τ_eff | R_bant(son) | R_bant(Hkeskin) | r |
|---|---|---|---|
| 0.4607 | 1.15761 | 1.22372 | 0.94598 |
| 0.4997 | 1.17807 | 1.25583 | 0.93808 |
| 0.5393 | 1.19865 | 1.28572 | 0.93227 |
| 0.5792 | 1.21504 | 1.31574 | 0.92347 |
| 0.6188 | 1.23202 | 1.34244 | 0.91775 |
| 0.6580 | 1.24850 | 1.37193 | 0.91003 |
| 0.6981 | 1.25526 | 1.38511 | 0.90625 |
| 0.7384 | 1.25783 | 1.40381 | 0.89602 |
| 0.7770 | 1.25002 | 1.40189 | 0.89167 |

`Σ|A_VS| = 25.103417316` (Hk: 27.815371002), `ΣA²_VS = 0.256832569`
(Hk: 0.292663602; oran **0.877569**).
clamp'e giren çizgiler: 14015/15450 **ama** ΣA²'nin yalnız **%7.91**'i;
doğrusal-uzatma alternatifi ΣA²'yi **%0.31** değiştiriyor ⇒ clamp seçimi
zarf gücünün binde üçünü etkiliyor (tanı, `180a`).

> **BİLİNEN SINIR — ölçümden önce yazıldı, süslenmeyecek.** `R_bant`'ın
> kendisi 176-HÜKÜM(vi)'ye göre **kilit tarafından da beslenir** (karışık
> vekillerde 0.84–1.01'e düşüyor). Dolayısıyla `r(τ)` **saf bir zarf
> niceliği değildir**; "gerçeğin ölçülen çizgi-gücü profilinin nominal
> merdivene taşınması"dır. ZARF PAYI bu tanıma **GÖREdir**. Bu, çapa
> bağımsızlığını bozmaz (hiçbir üçüncü referans gaza başvurulmuyor), ama
> "saf zarf" iddiasını sınırlar.

### 0.3 Ayrışım formülleri — **İKİ KONVANSİYON, İKİSİ DE RAPORLANIR**

**Doğrusal:**
```
ΔN      := N(son) − N(Hk)
ZARF_N  := ⟨N⟩_VS − ⟨N⟩_VF          [ÇAPASIZ — iki karışık ailenin farkı]
KİLİT_N := ΔN − ZARF_N               [ARTIK: kilit + zarf×kilit etkileşimi]
p_zarf  := ZARF_N/ΔN ,  p_kilit := KİLİT_N/ΔN ,  p_zarf + p_kilit ≡ 1
```
**Log:**
```
ΔlogN      := log N(son) − log N(Hk)
ZARFlog_N  := ⟨log N⟩_VS − ⟨log N⟩_VF
KİLİTlog_N := ΔlogN − ZARFlog_N
p_zarf^log := ZARFlog/ΔlogN ,  p_kilit^log := 1 − p_zarf^log
```

**ESAS para birimi (HAM):** `M (=KALİB_u2, lo=0.60)`, `Q_E`, `ρ_E`,
`Q_X`, `ρ_X`.
**EK (yalnız ek, oran dili):** `g_E`, `g_X`, `θ`, `g_cal`.

**Hata çubuğu:** `se(⟨N⟩_A) = sd(N_A, ddof=1)/√n_A` (n_VS = 2 ⇒
`se = |N₁−N₂|/2`); `se(ZARF) = √(se_VS² + se_VF²)`;
`se(KİLİT) = se(ZARF)`; `se(p) = se(ZARF)/|ΔN|`.

**İşaret kuralı (dondurulmuş):** ZARF ile ΔN aynı işaretli değilse pay
`>1` veya `<0` çıkar; bu bir hata değil **ölçümdür** ve aynen yazılır
(176'nın `g_E` dersi).

### 0.4 H-180a — çifte-sayım hakeminin **kararı ölçümden önce**

Ölçülen: `169_k2b.main(g, 0.40, 0.95) → ortalama['111_hepsi'] =: ρ₃(g)`
(3-bacak kırpma aktarımı, lo ∈ [0.52, 0.68] bantlarının ortalaması).
Gazlar: `VF1…VF4`.

Çapalar: `ρ₃(Hkeskin) = 0.5284`, `ρ₃(son) = 0.5356`,
`GAUSS3 = (2/π)^{3/2} = 0.507949`.

```
C   := (ρ₃(Hk) − ⟨ρ₃⟩_VF) / (ρ₃(Hk) − GAUSS3)        ÇÖKÜŞ KESRİ
z_G := (⟨ρ₃⟩_VF − GAUSS3) / se_VF
```

| dal | koşul | hüküm |
|---|---|---|
| **A** | `C ≥ 0.70` **ve** `\|z_G\| ≤ 3` | kırpma aktarımı **kilidin bir GÖRÜNÜMÜ**dür; %93 kırpma-kapanışı ile %91.4 kilit defteri aynı mekanizmanın iki defteri; **çifte sayım YOKTUR** |
| **B** | `C ≤ 0.30` | **iki ayrı bileşen**; %184 gerçek çifte sayımdır; paylar yeniden bölünmelidir |
| **HÜKÜMSÜZ** | `0.30 < C < 0.70`, ya da `C ≥ 0.70` ama `\|z_G\| > 3`, ya da tohum saçılımı `\|ρ₃(Hk) − GAUSS3\|`'ün yarısını aşarsa | HÜKÜMSÜZ yazılır, kurtarma yok |

### 0.5 Ön-kayıt anında önbellekten okunan ÇAPALAR (`180a`)

| gaz | M | Q_E | ρ_E | Q_X | ρ_X | g_E | θ |
|---|---|---|---|---|---|---|---|
| `son` | 0.273907 | 0.875748 | 2.179597 | 0.467931 | 1.584774 | 0.598206 | 0.921938 |
| `Hkeskin` | 0.258590 | 0.915576 | 2.199323 | 0.462127 | 1.562364 | 0.583701 | 0.893334 |
| `HA4` | 0.231073 | 0.795861 | 2.103165 | 0.430387 | 1.519587 | 0.621589 | 0.723573 |
| `VF1` | 0.422643 | 0.297426 | 1.130635 | 0.193999 | 1.040067 | 0.736939 | 0.866671 |
| `VF2` | 0.425711 | 0.297109 | 1.128368 | 0.164799 | 0.988779 | 0.736692 | 0.832136 |
| `VF3` | 0.430532 | 0.299440 | 1.133933 | 0.190991 | 1.035086 | 0.735928 | 0.879713 |
| `VF4` | 0.428004 | 0.298609 | 1.136974 | 0.204761 | 1.059888 | 0.737365 | 0.891712 |

---

## 1. K1 — İNŞA: GERÇEK-ZARFLI KARIŞIK VEKİLLER (`180b_vs_insa.py`)

Makine **176'nınkidir, aynen**: alan değerlendiricisi
`176_configs/176_vekil_cekirdek.py` (faz taşır; φ≡0'da `164_insa.S_ve_Sp`
ile **bit-bit** aynı), çözücü `164_insa.coz_sadakatli` (ızgara braketi +
SIRALI İLK-KÖK + korumalı Newton, h = 0.015, nz = 300000, c = −½) —
**kopyalanan tek satır yok**, yalnız `S_par`/`SSp_par` isimleri koşudan
önce faz taşıyan sürümlerle değiştirildi. Tek fark genlik dizisidir.

### T1 — İnşa kapıları (hepsi 180a'da dondurulmuştu)

| kapı | eşik | `VS1` (tohum 1) | `VS2` (tohum 2) |
|---|---|---|---|
| V1 `sha256(A_VS)` = ön-kayıttaki | birebir | `4ac2f827…52e1` ✓ | `4ac2f827…52e1` ✓ |
| V2 `maks\|A_VS − a_Hk·r(τ)\|` | = 0.0 | **0.0e+00** ✓ | **0.0e+00** ✓ |
| V3 φ≡0 sağlaması (ΔS / ΔS′) | = 0.0 | **0.0e+00 / 0.0e+00** ✓ | **0.0e+00 / 0.0e+00** ✓ |
| V4 ilk-kök hücre (benzersiz) | 300000/300000 | **300000 / 300000** ✓ | *(aşağıda)* |
| V5 maks\|F\| | ≤ 1e−8 | **1.863e−09**, aşan 0 ✓ | *(aşağıda)* |
| V6 sıralılık | TAM | **TAM**, min Δz = 0.098913 ✓ | *(aşağıda)* |
| σ_ds | — | 0.41837 | *(aşağıda)* |
| L (z'den) | — | 12.029593229 | *(aşağıda)* |
| süre | — | 10.8 dk | *(aşağıda)* |

> **Genlik dizisi iki tohumda BİT-BİT aynıdır** — `sha256(A_VS1) =
> sha256(A_VS2) = 4ac2f82748a2a7ac8973cc38d514d4650bede00014a098c8963242455d0e52e1`,
> ve bu dize **ön-kayıtta**, herhangi bir inşa koşmadan önce yazılmıştı.
> Vekiller arasındaki **tek fark φ_q dizisidir**.

**VS2 (tohum 2):** V1 ✓ (aynı sha), V2 **0.0e+00** ✓, V3 **0.0/0.0** ✓,
V4 **300000/300000** ✓, V5 **1.863e−09**, aşan 0 ✓, V6 **TAM**,
min Δz = 0.100535 ✓; σ_ds = 0.41989; L = 12.029593229; süre 10.7 dk.

> **İNŞA KAPILARININ HEPSİ (V1–V6) İKİ TOHUMDA DA GEÇTİ.**
> Karşılaştırma için: `VF` ailesi σ_ds = 0.4416/0.4427, `Hkeskin`
> L = 12.029593242. Zarf zayıfladığı için (ΣA² oranı 0.8776) VS'nin
> σ_ds'i VF'ninkinin altında — beklenen yön.

---

## 2. K2 — H-180a: ÇİFTE-SAYIM HAKEMİ (`180d_kirpma_kos.py`, `180e_hakem.py`)

`169_k2b.main(g, 0.40, 0.95)` **aynen** koşuldu (174e'nin kullandığı
çağrının birebir aynısı); yeni makine yazılmadı. Hüküm, 180a'da
dondurulmuş kuralla `180e` tarafından verildi.

### T5 — 3-bacak kırpma aktarımı (lo ∈ [0.52, 0.68] bantlarının ortalaması)

| gaz | ρ(E) | ρ(Xa) | ρ(Xb) | ρ(XaXb) | ρ(EXa) | **ρ₃ (hepsi)** |
|---|---|---|---|---|---|---|
| `Hkeskin` | 0.7947 | 0.8463 | 0.8468 | 0.6220 | 0.6952 | **0.5284** |
| `son` | 0.8229 | 0.8366 | 0.8367 | 0.6023 | 0.7147 | **0.5356** |
| `VF1` | 0.6487 | 0.8290 | 0.8179 | 0.6693 | 0.5467 | **0.4474** |
| `VF2` | 0.6374 | 0.8283 | 0.8163 | 0.6692 | 0.5362 | **0.4477** |
| `VF3` | 0.6452 | 0.8244 | 0.8365 | 0.6778 | 0.5427 | **0.4635** |
| `VF4` | 0.6529 | 0.8206 | 0.8095 | 0.6575 | 0.5530 | **0.4553** |
| **GAUSS** | 0.7979 | 0.7979 | 0.7979 | 0.6366 | 0.6366 | **0.5079** |

`⟨ρ₃⟩_VF = 0.453484`, sd = 0.007615, se = 0.003807 (n = 4).

```
C   = (0.5284 − 0.453484)/(0.5284 − 0.507949) = +3.6632
z_G = (0.453484 − 0.507949)/0.003807         = −14.305
tohum saçılımı = 0.016091   (ön-kayıtlı HÜKÜMSÜZ eşiği 0.010225)
```

### ⚠ ÖN-KAYITLI HÜKÜM: **HÜKÜMSÜZ**

> Dondurulmuş kuralın ilk maddesi tetiklendi: **tohum saçılımı
> (0.016091), iki çapa arasındaki mesafenin yarısını (0.010225) aşıyor**
> ⇒ sınavın çözünürlüğü, "Gauss'a çöktü mü" sorusunu ön-kayıtta
> istenen kesinlikte yanıtlamaya yetmiyor. **Kurtarma yok.**

**Ne yine de kesin olarak ölçüldü** (kurtarma değil, aynı koşunun
kaçınılmaz okuması):

1. **DAL B (iki ayrı bileşen, "%184 gerçek çifte sayım") ÇÜRÜDÜ.**
   Koşulu `C ≤ 0.30` idi. Dört tohumun her biri **tek başına**:
   C(VF1) = 3.96, C(VF2) = 3.95, C(VF3) = 3.17, C(VF4) = 3.57.
   En yakın tohum bile eşiğin **10 katı** ötesinde. Faz kilidi
   söküldüğünde kırpma aktarımının Hkeskin'e göre fazlası **ayakta
   kalmıyor** — ortalamaya, konvansiyona, tohum seçimine bağlı değil.
2. **DAL A da mühürlenmedi.** Karışık vekiller Gauss'a çökmedi, Gauss'un
   **altına düştü**: `⟨ρ₃⟩_VF = 0.4535`, Gauss'un **%10.7 altında**,
   `z_G = −14.3`. "Karışık gaz = Gauss alanı" **yanlıştır**.
3. **Fazlanın ve açığın adresi tek bacak: E.** Xa ve Xb bacakları
   karıştırmaya **kör** (0.82–0.83; Gauss 0.7979, +3…+5% — kilitli
   gazlarda da aynı). Bütün hareket E bacağında: `son` 0.8229 →
   `Hkeskin` 0.7947 → karışık **0.6374–0.6529** (Gauss'un **%18.5
   altı**). Yani kırpma aktarımı anomalisi **η artık alanının kendi
   kırpma tepkisidir**, çizgi alanının değil.

> **Bu, ön-kayıtlı ikili dalın hiçbirinin öngörmediği ÜÇÜNCÜ bir
> olgudur** ve süslenmeden kaydediliyor: kilit sökülünce E bacağı
> Gauss'un üstünden değil **altına** geçiyor. Kilit, kırpma aktarımının
> Gauss-üstü fazlasının **gerekli koşuludur** (fazla kilitsiz hayatta
> kalmıyor) ama karışık gazın kırpma defteri Gauss'la da açıklanmıyor.

---

## 3. K1 (ÖLÇÜM) — VS'nin defteri (`180c_olcum.py`)

`167_olcum.kos(ad, 0.40, 0.95, kule=0, düz=0)` — VF ile **birebir aynı
çağrı**. Özdeşlik denetimi `|ΔK|/K`: E 4.2e−15 / 3.0e−15, X 8.6e−16 /
1.2e−15 (172b'nin Ö-A özdeşliği kırılmadı).

### T2 — DEFTER (HAM esas | oran EK)

| gaz | M | Q_E | ρ_E | Q_X | ρ_X | \| g_E | g_X | θ | g_cal |
|---|---|---|---|---|---|---|---|---|---|
| `son` | 0.273907 | 0.875748 | 2.179597 | 0.467931 | 1.584774 | 0.598206 | 0.704733 | 0.921938 | 0.297099 |
| `Hkeskin` | 0.258590 | 0.915576 | 2.199323 | 0.462127 | 1.562364 | 0.583701 | 0.704213 | 0.893334 | 0.289467 |
| `VF1` | 0.422643 | 0.297426 | 1.130635 | 0.193999 | 1.040067 | 0.736939 | 0.813474 | 0.866671 | 0.487663 |
| `VF2` | 0.425711 | 0.297109 | 1.128368 | 0.164799 | 0.988779 | 0.736692 | 0.833331 | 0.832136 | 0.511588 |
| `VF3` | 0.430532 | 0.299440 | 1.133933 | 0.190991 | 1.035086 | 0.735928 | 0.815483 | 0.879713 | 0.489401 |
| `VF4` | 0.428004 | 0.298609 | 1.136974 | 0.204761 | 1.059888 | 0.737365 | 0.806809 | 0.891712 | 0.479981 |
| **`VS1`** | **0.444879** | **0.288527** | **1.150329** | **0.189736** | **1.051637** | 0.749179 | 0.819580 | 0.884044 | 0.503232 |
| **`VS2`** | **0.449528** | **0.287084** | **1.144479** | **0.160029** | **1.000126** | 0.749157 | 0.839991 | 0.850423 | 0.528594 |

### Zarf aktarımının sağlaması (tanı, kapı değil — `180c` + `180a`)

| τ_eff | ⟨R_bant⟩_VS | ⟨R_bant⟩_VF | ölçülen oran | tarif `r(τ)` | fark |
|---|---|---|---|---|---|
| 0.4607 | 0.9651 | 1.0098 | 0.9557 | 0.9460 | +1.0% |
| 0.5394 | 0.9436 | 0.9896 | 0.9535 | 0.9323 | +2.3% |
| 0.6189 | 0.9125 | 0.9605 | 0.9500 | 0.9177 | +3.5% |
| 0.6986 | 0.8533 | 0.9093 | 0.9384 | 0.9063 | +3.6% |
| 0.7782 | 0.7886 | 0.8345 | 0.9449 | 0.8917 | +6.0% |

Aktarım çalıştı (oran her bantta 0.94–0.96), **ama ölçülen oran τ'da
tarifin eğimini taşımıyor**: gerçek gazın yüksek-τ çizgi gücü açığının
bir kısmı, kilitsiz bir gazda genlik ölçeklemesiyle **yeniden
üretilemiyor**. Bu, ön-kayıtta yazılmış sınırın ölçülmüş hâlidir:
`r(τ)`'nun τ-eğimi kısmen bir **kilit** etkisidir, saf zarf değil.

---

## 4. K3 — NİHAİ ÇAPASIZ DEFTER (`180f_defter.py`)

`ΔM(son↔Hk)`: doğrusal **+0.015316**, log **+0.057542** — 175/179'un
`Δlog M = +0.057542` sayısıyla **birebir**.

### T3 — AYRIŞIM, **DOĞRUSAL** konvansiyon

| N | ΔN | ⟨N⟩_VS | ⟨N⟩_VF | **ZARF ± se** | **KİLİT ± se** | p_zarf ± se | p_kilit ± se |
|---|---|---|---|---|---|---|---|
| **M** | +0.015316 | 0.447204 | 0.426723 | **+0.020481 ± 0.002867** | **−0.005165 ± 0.002867** | **+1.337 ± 0.187** | **−0.337 ± 0.187** |
| **Q_E** | −0.039828 | 0.287806 | 0.298146 | −0.010341 ± 0.000900 | −0.029487 ± 0.000900 | +0.260 ± 0.023 | +0.740 ± 0.023 |
| **ρ_E** | −0.019726 | 1.147404 | 1.132478 | +0.014927 ± 0.003479 | −0.034653 ± 0.003479 | −0.757 ± 0.176 | +1.757 ± 0.176 |
| **Q_X** | +0.005804 | 0.174883 | 0.188637 | −0.013755 ± 0.017103 | +0.019559 ± 0.017103 | −2.37 ± 2.95 | +3.37 ± 2.95 |
| **ρ_X** | +0.022410 | 1.025882 | 1.030955 | −0.005073 ± 0.029828 | +0.027483 ± 0.029828 | −0.226 ± 1.331 | +1.226 ± 1.331 |
| *(EK)* g_E | +0.014505 | 0.749168 | 0.736731 | +0.012437 ± 0.000302 | +0.002068 ± 0.000302 | +0.857 ± 0.021 | +0.143 ± 0.021 |
| *(EK)* θ | +0.028604 | 0.867234 | 0.867558 | −0.000324 ± 0.021170 | +0.028929 ± 0.021170 | −0.011 ± 0.740 | +1.011 ± 0.740 |
| *(EK)* g_X | +0.000520 | 0.829785 | 0.817274 | +0.012511 ± 0.011672 | −0.011991 ± 0.011672 | +24.0 ± 22.4 | −23.0 ± 22.4 |
| *(EK)* g_cal | +0.007632 | 0.515913 | 0.492158 | +0.023755 ± 0.014385 | −0.016123 ± 0.014385 | +3.11 ± 1.88 | −2.11 ± 1.88 |

### T4 — AYRIŞIM, **LOG** konvansiyon

| N | ΔlogN | **ZARFlog ± se** | **KİLİTlog ± se** | p_zarf ± se | p_kilit ± se |
|---|---|---|---|---|---|
| **M** | +0.057542 | **+0.046890 ± 0.006520** | **+0.010652 ± 0.006520** | **+0.815 ± 0.113** | **+0.185 ± 0.113** |
| **Q_E** | −0.044475 | −0.035297 ± 0.003089 | −0.009178 ± 0.003089 | **+0.794 ± 0.070** | **+0.206 ± 0.070** |
| **ρ_E** | −0.009010 | +0.013095 ± 0.003044 | −0.022105 ± 0.003044 | −1.454 ± 0.338 | +2.454 ± 0.338 |
| **Q_X** | +0.012481 | −0.076166 ± 0.096981 | +0.088647 ± 0.096981 | −6.10 ± 7.77 | +7.10 ± 7.77 |
| **ρ_X** | +0.014242 | −0.004925 ± 0.029113 | +0.019167 ± 0.029113 | −0.346 ± 2.044 | +1.346 ± 2.044 |
| *(EK)* g_E | +0.024547 | +0.016741 ± 0.000410 | +0.007806 ± 0.000410 | +0.682 ± 0.017 | +0.318 ± 0.017 |
| *(EK)* θ | +0.031518 | −0.000228 ± 0.024488 | +0.031746 ± 0.024488 | −0.007 ± 0.777 | +1.007 ± 0.777 |
| *(EK)* g_X | +0.000739 | +0.015188 ± 0.014102 | −0.014450 ± 0.014102 | +20.6 ± 19.1 | −19.6 ± 19.1 |
| *(EK)* g_cal | +0.026024 | +0.047118 ± 0.028129 | −0.021094 ± 0.028129 | +1.811 ± 1.081 | −0.811 ± 1.081 |

### 4.1 Hangi satır bir şey söylüyor, hangisi söylemiyor

**Hata çubuğu (tohum saçılımı) taşıyabilenler:**

| N | konvansiyon | ZARF payı | KİLİT payı | kaç σ |
|---|---|---|---|---|
| **M** | log | **%81.5 ± %11.3** | **%18.5 ± %11.3** | zarf 7.2σ, kilit 1.6σ |
| **Q_E** | log | **%79.4 ± %7.0** | **%20.6 ± %7.0** | zarf 11.4σ, kilit 3.0σ |
| **Q_E** | doğrusal | %26.0 ± %2.3 | %74.0 ± %2.3 | — |
| **M** | doğrusal | %133.7 ± %18.7 | −%33.7 ± %18.7 | — |
| *(EK)* g_E | log | %68.2 ± %1.7 | %31.8 ± %1.7 | — |

**HÜKÜMSÜZ satırlar (hata çubuğu payı yutuyor):** `Q_X`
(se ±%295 doğrusal, ±%777 log), `ρ_X` (±%133 / ±%204), `θ`
(±%74 / ±%78), `g_X` (±%2240 / ±%1910), `g_cal` (±%188 / ±%108).
Bu niceliklerde ΔN, iki tohumluk VS saçılımının altında kalıyor —
**hiçbir pay okunamaz**, süslenmez.

**`ρ_E` bir ölçümdür, bir hata değil:** Δρ_E = −0.0197 iken ZARF_ρ_E
**+0.0149** (ters işaret) ⇒ pay > 1. Ön-kayıtın "işaret kuralı"
maddesi tam olarak bunu öngörüp yazmıştı; kaydediliyor, kurtarılmıyor.

### 4.2 KONVANSİYON, ÇAPASIZ DÜNYADA DA KAÇMIYOR

> 179/S4'ün "konvansiyon karması" eleştirisi **çapa kaldırıldıktan sonra
> da ayakta**: `M`'de doğrusal zarf payı %133.7, log'da %81.5;
> `Q_E`'de doğrusal %26.0 **↔** log %79.4 — yani aynı ham nicelik iki
> konvansiyonda **birbirinin tersi** hikâye anlatıyor. Sebep ölçülüdür:
> karışık ailelerin mutlak seviyesi kilitli çiftinkinden çok uzak
> (`Q_E`: 0.876/0.916 ↔ 0.288/0.298; `M`: 0.274/0.259 ↔ 0.447/0.427),
> bu yüzden doğrusal FARK ile ORANSAL fark aynı şeyi ölçmüyor.
> **İki konvansiyon da ön-kayıtta dondurulmuştu ve ikisi de burada.**
> Aralarında seçim yapılmıyor.

---

## 5. K4 — 176'NIN KARIŞTIRMA ÇÖKÜŞLERİYLE ÇAPRAZ TUTARLILIK (`180f`)

Türetim (varsayımsız, ölçümden önce ön-kayıtta):

```
K_N  := exp⟨log N⟩_VF / N(Hk)      Hk zarfında karıştırma çöküşü  (=176'nınki)
K'_N := exp⟨log N⟩_VS / N(son)     son zarfında AYNI çöküş        (=180'in yenisi)
⇒   KİLİTlog_N  ≡  log K_N − log K'_N          (ÖZDEŞLİK)
```

| N | K_N (Hk zarfı) | K′_N (son zarfı) | log K − log K′ | KİLİTlog | kalıntı |
|---|---|---|---|---|---|
| **M** | **+%65.01** | **+%63.27** | +0.010652 | +0.010652 | 0.00e+00 |
| **Q_E** | **−%67.44** | **−%67.14** | −0.009178 | −0.009178 | 4.2e−17 |
| **ρ_E** | −%48.51 | −%47.36 | −0.022105 | −0.022105 | −2.8e−17 |
| Q_X | −%59.31 | −%62.76 | +0.088647 | +0.088647 | 0.00e+00 |
| ρ_X | −%34.03 | −%35.29 | +0.019167 | +0.019167 | 3.1e−17 |
| *(EK)* g_E | +%26.22 | +%25.24 | +0.007806 | +0.007806 | 0.00e+00 |
| *(EK)* θ | −%2.92 | −%5.95 | +0.031746 | +0.031746 | 0.00e+00 |

176'nın manşet çöküşleri **birebir yeniden üretildi** (`Q_E` −%67.4 ↔
176 HÜKÜM(ii)'nin "−%68"i; `ρ_E` −%48.5 ↔ "−%49"; `g_E` +%26.2 ↔
"+%26"; `M` +%65.0 ↔ "+%64").

> **K4'ÜN TEK CÜMLESİ:** Faz kilidinin çöküşü **iki zarfta da neredeyse
> aynıdır** (`M`: +%65.0 ↔ +%63.3; `Q_E`: −%67.4 ↔ −%67.1; `g_E`:
> +%26.2 ↔ +%25.2). Kilit **devasadır ama zarfa duyarsızdır**, bu
> yüzden `son`−`Hkeskin` farkından **büyük ölçüde sadeleşir**. Log
> dilindeki "kilit payı", ayrı bir kestirici değil, **bu iki çöküşün
> farkının ta kendisidir** — kalıntı makine sıfırı (≤ 4.2e−17).

---

## 6. T6 — ESKİ ÇAPA-BAĞIL SAYILAR ↔ 180'in ÇAPASIZ KARŞILIĞI

| kestirici | eski değer | konvansiyon / çapa | 180'in çapasız karşılığı | fark |
|---|---|---|---|---|
| ΔM ayrışımı | **%8.6 kesim / %91.4 kilit** | log, çapa = (Hkeskin, HA4) | **%81.5 zarf / %18.5 kilit** (log) | paylar **ters dönüyor** |
| ΔM ayrışımı | (yok) | doğrusal | **%133.7 zarf / −%33.7 kilit** | — |
| E-kanal ara-konumu f(g_E) | 0.3828 (doğrusal) / 0.3903 (log) | HA4 çapası | g_E log zarf payı **0.682 ± 0.017** | +%75 |
| f(μ̂²_E) | 0.3968 | HA4 çapası | *(karşılığı yok — μ̂²_E ayrışımı ölçülmedi)* | — |
| f_b (τ bantları) | 0.3975 / 0.3871 | HA4 çapası | *(karşılığı yok — bant kesri VS'de ölçülmedi)* | — |
| Q_E payı | (yok — 176'nın borcu) | — | **log %79.4 ± %7.0 zarf / %20.6 ± %7.0 kilit**; **doğrusal %26.0 / %74.0** | — |

> **Eski "%39"un nerede durduğu:** o sayı, `HA4` (erfc-kesim ikizi)
> çapasına göre `son`'un `Hkeskin`–`HA4` ekseninde nereye düştüğünü
> söylüyordu. 180'in çapasız sayısı **aynı soruyu sormuyor**: burada
> referans üçüncü bir gaz değil, **aynı zarfların kilitsiz sürümleri**.
> İki sayıyı doğrudan eşitlemek yanlış olur; tablo yalnızca **nerede
> durduklarını** yan yana koyar. Ortak nokta: her iki defterde de
> "kilit payı" **konvansiyona ve referansa aşırı duyarlıdır** —
> %91.4 ↔ %18.5 ↔ −%33.7 aralığı bunun ölçüsüdür.

---

## 7. FİGÜR

`180_nedensel_defter.png` (`180g_figur.py`) — dört panel:
**(a)** H-180a: ρ₃ çubukları (`son`, `Hkeskin`, VF1–VF4), Gauss çizgisi,
C ve z_G; **(b)** çapasız ayrışım, doğrusal konvansiyon (hata
çubuklarıyla; `Q_X` ölçeği taşıyor, değerler etikette);
**(c)** aynısı log konvansiyonunda + eski %91.4 çizgisi;
**(d)** zarf aktarımının sağlaması: tarif `r(τ)` ↔ ölçülen
⟨R_bant⟩_VS/⟨R_bant⟩_VF.

---

# HÜKÜM

## (i) K2 — ÇİFTE-SAYIM HAKEMİ: **HÜKÜMSÜZ** (ön-kayıtlı çözünürlük maddesi)

Tohum saçılımı 0.016091, iki çapa arasındaki mesafenin yarısını
(0.010225) aştı. **Kurtarma yok.** Ama iki alt-iddia kesin:
**DAL B öldü** (dört tohumun her biri tek başına C ≈ 3.2–4.0; eşik ≤ 0.30
— "iki ayrı bileşen / %184 gerçek çifte sayım" okuması **desteksiz**),
**DAL A da mühürlenmedi** (karışık gaz Gauss'a çökmedi, **%10.7 altına**
düştü, z_G = −14.3). Anomalinin adresi **tek bacaktır: E**.

## (ii) K3 — ÇAPASIZ DEFTER: KİLİT PAYI **KÜÇÜLDÜ**, ama konvansiyondan kurtulmadı

Log konvansiyonunda, iki bağımsız HAM nicelik aynı şeyi söylüyor:
`M` **%81.5 ± %11.3 zarf / %18.5 ± %11.3 kilit**, `Q_E` **%79.4 ± %7.0
zarf / %20.6 ± %7.0 kilit`. Doğrusal konvansiyonda `M` zarf payı %133.7,
`Q_E` kilit payı %74.0. **Ortak olan:** hiçbir okuma "%91.4 kilit"i
desteklemiyor — o sayı bu defterde **yeniden üretilemedi**.

## (iii) NEDEN: KİLİT DEVASADIR AMA ZARFA DUYARSIZDIR (K4)

`M` karıştırmayla iki zarfta da ~+%64 (65.01 ↔ 63.27), `Q_E` −%67
(67.44 ↔ 67.14). 176'nın "kilit nedensel ve devasa" bulgusu **ayakta**
— ama bu kilidin, `son` ile `Hkeskin`'i **birbirinden ayıran** şey
olduğu anlamına gelmiyor: aynı büyüklükte olduğu için farktan
sadeleşiyor. **Kilidin VARLIĞI ile kilidin FAZLAYA KATKISI ayrı
sorulardır**; 176 birincisini ölçtü, 180 ikincisinin küçük olduğunu
ölçüyor.

## (iv) ÖLÇÜLMÜŞ SINIR (kurtarma değil, sınırın kendisi)

Zarf tarifi `r(τ)` ölçülen `R_bant` oranından geldiği için, 176-HÜKÜM(vi)
uyarınca kısmen kilit taşır. Sağlama bunu **görünür kıldı**: ölçülen
⟨R_bant⟩_VS/⟨R_bant⟩_VF oranı 0.938–0.956 arasında **düz**, tarif ise
0.892–0.946 arasında **eğimli**. Yani gerçek gazın yüksek-τ açığının bir
kısmı kilitsiz bir gazda genlik ölçeklemesiyle yeniden üretilemiyor.
Bu, ZARF payını **kısmen şişirmiş** olabilir; yönü belli, büyüklüğü
ölçülmedi ⇒ %81.5'in bir **üst sınır tadı** vardır. Ölçülmeyen şey
ölçüldü sayılmaz.

## (v) NE İSTENDİ, NE VERİLDİ

179'un ŞART'ı üç ayaklıydı: ① (Q_E, ρ_E) para biriminde nedensel
ayrışım — **verildi** (Q_E'de hata çubuğuyla; ρ_E'de işaret dönüyor ve
bu yazıldı); ② kırpma ↔ kilit çifte-sayım yüzleşmesi — **HÜKÜMSÜZ**,
ama DAL B öldü; ③ W_pos'un KALİB/θ zincirindeki doğruluğu — **180'de
ele alınmadı**, borç duruyor.

---

# 🔒 MÜHÜR

> **Gerçeğin fazlasının çapasız zarf payı log konvansiyonunda
> %81.5 ± %11.3 (Q_E'de %79.4 ± %7.0), kilit payı %18.5 ± %11.3
> (Q_E'de %20.6 ± %7.0)'dir — doğrusal konvansiyonda zarf payı %133.7,
> kilit payı −%33.7'dir ve iki konvansiyon arasındaki bu uçurum
> ölçümün kendisidir; kırpma aktarımının kilidin görünümü mü ayrı
> bileşen mi olduğu ise HÜKÜMSÜZDÜR — ama "ayrı bileşen" dalı
> (çifte sayım) dört tohumun her birinde ayrı ayrı ÖLDÜ.**

**Ne mühürlenmedi:** "%91.4 kilit" — bu defterde hiçbir konvansiyonda
yeniden üretilemedi; "kırpma = kilidin görünümü" — çözünürlük yetmedi.

---

## Sıradaki adım (bu ölçümün işaret ettiği)

1. **VS'nin ÜÇÜNCÜ ve DÖRDÜNCÜ TOHUMU** (~2 × 16 dk). `M` ve `Q_E`
   dışındaki bütün satırlar iki tohumluk VS saçılımında boğuldu;
   n = 4'e çıkmak `θ`, `ρ_X`, `g_cal` satırlarını okunur yapabilir ve
   `M`'nin ±%11.3'ünü ~±%8'e indirir. **En ucuz, en yüksek getirili.**
2. **E BACAĞININ KIRPMA DEFTERİ.** K2'nin gerçek bulgusu şu: bütün
   hareket `ρ(E)`'de (`son` 0.823 → karışık 0.637; Xa/Xb kör). Bu, η
   artık alanının kırpma tepkisinin doğrudan ölçümünü ister — yeni
   inşa gerekmez, veri diskte (`169/K2b_VF*.json`, `167/C_VF*.json`).
3. **ARA KARIŞTIRMA (α taraması).** Hem 179/S4'ün "ara doz hiç
   örneklenmedi" eleştirisini kapatır hem K2'nin çözünürlük ölümünü
   giderir: `φ_q ~ U(0, 2πα)` ile `ρ₃(α)` eğrisi, gerçeğin α_eş'ini
   verir. Beş inşa (~50 dk).
4. **W_pos BORCU** (179 ŞART ③) — 180'de ele alınmadı.

