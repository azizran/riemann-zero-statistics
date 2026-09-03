# 170 — EJDERHA: SONLU-DOYUM DÜZELTME YASASI c(λ) (3 Eylül 2026, gece)

Kalem: `KALEM_SONLU_DOYUM_03EYL2026.md` (H-D1, H-D2; K0–K3 + bonus).

Betikler (`170_configs/`):

| betik | ne yapar | koşu | ham çıktı |
|---|---|---|---|
| `170a_muhasebe.py` | **K0** — muhasebe tablosu (yeni ölçüm yok) | 12 s | `170/K0.json`, `log_K0.txt` |
| `170l_uye_ussu.py` | **K0(d)** — üye üssü ν'nün iki yönde ayrıştırılması | 2 s | `170/UYE_USSU.json` |
| `170b_T_yasasi.py` | **K1** — T(σ) türetimi + MC doğrulaması + bacak-başına σ_b | 11 dk (4 gaz) | `170/K1_T.json`, `K1_T_kesim.json` |
| `170_insa60.py` | λ=0.60 gazının inşası (`167_insa.main` aynen) | 10.4 dk | `167/z_L060.npy`, `insa_L060.json` |
| `170c_onkayit60.py` | **K2 ÖN KAYIT** (korelatöre bakmadan) | 1.3 dk | `170/ONKAYIT_L060.json` |
| `170d_olcum60.py` | **K2 ölçüm** (`167_olcum.kos` aynen) | 3.7 dk | `167/C_L060.json` |
| `170g_yuzlesme.py` | **K2 yüzleşme** (ön kayıt ↔ ölçüm) | 20 s | `170/K2_yuzlesme.json` |
| `170e_bacak_aktarim.py` | **K2c/K3d** — bacak-bacak korelatör aktarımı (`169_k2b.main` aynen) | 4×5–7 dk | `169/K2b_*.json` |
| `170m_k3.py` | **K3** — gerçek gazın fazlası ve H-D2 | 2 s | `170/K3.json` |
| `170f_beta.py` | **BONUS** — (T5) açılımı, β, T/4 artığı | 3 s | `170/BONUS.json` |
| `170i_insa115.py` | EK KAPI: λ=1.15 gazının inşası | 14 dk | `167/z_L115.npy` |
| `170j_onkayit115.py` | EK KAPI ÖN KAYIT (λ=1.15) | 2 dk | `170/ONKAYIT_L115.json` |
| `170k_olcum115.py` | EK KAPI ölçüm (`167_olcum.kos` aynen) | 4.2 dk | `167/C_L115.json` |
| `170n_yuzlesme115.py` | EK KAPI yüzleşme | 25 s | `170/K2p_yuzlesme.json` |
| `170h_figur.py` | figür | 5 s | `170_sonlu_doyum.png` |

(`scratchpad` kökü: `/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/
71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad`)

Hiçbir ölçüm/çözüm parçası kopyalanmadı: `166_T1.bant_agg/_jk`,
`165_cekirdek.Model165/sentez/tayf_s/gaz`, `163_cekirdek.olc_cizgi/
bant_adaylari`, `166_bacak.karakteristik`, `167_insa.main`,
`167_olcum.kos`, `169_k1.*`, `169_k2.kuiper`, `169_k2b.main` **import
edildi**. Git'e dokunulmadı.

---

## 0. TEK CÜMLELİK HÜKÜM

> **H-D1 ÖLDÜ — hem büyüklükte hem İŞARETTE; ejderhanın adresi `W_X`
> PAYDASI çıktı.** Bacak-başına sarılmış-aktarım yasası
> `T(σ) = √(ρ/arcsin ρ)`, `ρ = 1−e^{−σ²}` türetildi (Monte-Carlo ile
> 5·10⁻⁵'te doğrulandı; iki-bacak biçimi `T₂(r) = (2/π)arcsin r/r`
> 169'un arcsine ara-değer noktalarını %0.8–3.5 içinde veriyor) ve
> bacak-başına `σ_b = 2πτ_bσ_Ĉ` ölçüldü. Ölçüm yasanın **temelini**
> yıktı: **tek bacaklar sarılmıyor** (`σ_b = 1.0–1.26 rad`, Kuiper
> `√N·V` = 154–240, düzgünlük eşiği 2.00); yalnız DÖRT-frekans toplamı
> sarılıyor (`σ_Σ = 3.6–5.1 rad`, `√N·V` = 1.2–1.8). Sonuç:
> `Π_b T(σ_b) = 0.771` — ölçülen `c = 0.4035`'in **1.90 katı** — ve
> yasa λ küçüldükçe `c`'nin **ARTMASINI** istiyor (λ = 1.15 → 0.60:
> 0.738 → 0.771 → 0.809 → 0.853 → 0.887), ölçülen ise düşüyor. **ÖN KAYITLI örneklem-dışı
> nokta λ = 0.60 (20:48:26'da yazıldı, 20:48:52'de ölçüldü):
> öngörü 0.4645, ölçüm 0.3689 ± 0.0072 ⇒ −13.2σ.** Korelatör
> düzeyinde bağımsız üçüncü bir ölüm: ölçülen üç-bacak kırpma aktarımı
> λ ile **ARTIYOR** (0.5284 → 0.5285 → 0.5542 → 0.5904), H-D1 ise
> azalmasını isterdi (→ 0.4940). Üstelik bu artışın işareti ve
> mertebesi tam olarak `T(σ)`'nın öngördüğüdür (dört bacağa uzatımda
> +%16.2 ölçülen ↔ +%15.0 yasa): **`T(σ)`
> kırpma aktarımını doğru tarif ediyor, `c` o aktarım DEĞİL.**
> **Buna karşılık K0 ejderhanın adresini buldu:** `c(λ)`'nın hareketinin
> TAMAMI `W_X` paydasındadır — ham `KALİB_u2` λ azaldıkça **BÜYÜYOR**
> (+%1.4/+%6.0/+%12.1) ve λ ≥ 1'de neredeyse SABİT (λ=1.15'te +%0.2);
> `W_X` ise sırasıyla +%7.1/+%15.6/+%22.3 ve −%6.1. Üye üssü
> `ν = d log KALİB/d log W_X` λ = 1.15 → 0.60 boyunca
> **−0.02 → 0.23 → 0.59 → 1.02** diye tırmanıyor ve **`c(λ)` eğrisi
> tam olarak ν(λ) eğrisidir**: ν = 1 olan yerde `c` düz
> (`c(0.60) = 0.3689` ≡ `c(0.70) = 0.3690`), ν ≈ 0 olan yerde `c`
> `1/W_X` hızıyla fırlıyor.
>
> **İKİNCİ ÖN-KAYITLI SINAV (λ = 1.15, kendi eklediğim kapı):**
> `c(1.15) = 0.4312 ± 0.0053`, yani `4/π² = 0.40528`'in **+%6.4
> ÜSTÜNDE (+4.9σ)**. `c` λ ile artmaya devam ediyor ve 4/π²'yi
> λ = 1.00'de **hiçbir özellik göstermeden kesip geçiyor**.
> **⇒ `c = 4/π²` bir doyum limiti DEĞİL, λ = 1.00'in tesadüfüdür.**
>
> **TAÇ (K3) ise POZİTİF:** gerçek ζ gazının sadakatli ikizine göre
> +3.00σ'lık fazlası, o gazın **ÖLÇÜLEN** bacak kırpma aktarımıyla
> **%86 kapanıyor** (+3.00σ → **+0.43σ**).

*(§4'teki tam hüküm tablosu, §K2/§K2′'deki ön-kayıt yüzleşmeleri ve
§K3 bağlayıcıdır. Figür: `170_sonlu_doyum.png`.)*

---

## K0 — MUHASEBE TABLOSU (`170a_muhasebe.py`, `log_K0.txt`)

Her `c` sayısının **koşul defteri**. Ortak bant penceresi her yerde
`lo ∈ [0.52, 0.68]`, `SNR ≥ 3`, `R_bant ≥ 0.98`, `τ_eff < 0.85` (beş bant);
`c_WX ≡ KALİB_u2/W_X` (168 §A2.8'in üyesi), `c_ampX ≡ KALİB_u2/(W_amp W_X)`
(JSON'daki `c_u2`, 167'nin üyesi), `KALİB_u2 = g_cal·θ`, `g_cal = g_E g_X²`.

| gaz | eksen | λ | kesim | φ | N | T/T₀ | σ_ds | σ_X̃ | σ_Ĉ | g_cal | θ | KALİB | W_amp | W_X | **c_WX ± σ_tot** | c_ampX |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **L115** ‡ | λ | 1.15 | τ≤1.00 | 0.0000 | 299998 | 1.00 | 0.4596 | 0.2584 | 0.2965 | 0.2910 | 0.8856 | 0.2577 | 0.6669 | 0.6003 | **0.4312 ± 0.0053** | 0.6485 |
| Hkeskin | TABAN | 1.00 | τ≤1.00 | 0.0000 | 299998 | 1.00 | 0.4313 | 0.2420 | 0.2777 | 0.2895 | 0.8884 | 0.2572 | 0.7004 | 0.6394 | **0.4035 ± 0.0018** | 0.5774 |
| **son** | TABAN (GERÇEK ζ) | — | gerçek | 0.0000 | 299998 | 1.00 | 0.4092 | 0.2338 | 0.2730 | 0.2971 | 0.9142 | 0.2716 | 0.7290 | 0.6607 | **0.4122 ± 0.0023** | 0.5664 |
| L085 | λ | 0.85 | τ≤1.00 | 0.0000 | 299998 | 1.00 | 0.3976 | 0.2227 | 0.2556 | 0.2912 | 0.8959 | 0.2609 | 0.7393 | 0.6851 | **0.3817 ± 0.0028** | 0.5171 |
| L070 | λ | 0.70 | τ≤1.00 | 0.0000 | 299998 | 1.00 | 0.3559 | 0.1990 | 0.2289 | 0.2986 | 0.9125 | 0.2725 | 0.7854 | 0.7395 | **0.3690 ± 0.0058** | 0.4703 |
| **L060** ‡ | λ | 0.60 | τ≤1.00 | 0.0000 | 299998 | 1.00 | 0.3213 | 0.1795 | 0.2073 | 0.3091 | 0.9323 | 0.2882 | 0.8212 | 0.7820 | **0.3689 ± 0.0072** | 0.4495 |
| K090 | kesim | 1.00 | τ≤0.90 | 0.4153 | 299998 | 1.00 | 0.4389 | 0.2446 | 0.2787 | 0.2877 | 0.8067 | 0.2321 | 0.6875 | 0.6297 | **0.3699 ± 0.0015** | 0.5393 |
| K070 | kesim | 1.00 | τ≤0.70 | 0.9270 | 299998 | 1.00 | 0.4443 | 0.2421 | 0.2766 | 0.3128 | 0.7142 | 0.2234 | 0.6769 | 0.6441 | **0.3480 ± 0.0031** | 0.5154 |
| HA4 | kesim | 1.00 | erfc 0.68/0.125 | 0.9098 | 299998 | 1.00 | 0.4343 | 0.2396 | 0.2755 | 0.3194 | 0.7086 | 0.2263 | 0.6916 | 0.6499 | **0.3492 ± 0.0073** | 0.5061 |
| E060 | kesim | 1.00 | erfc 0.60/0.125 | 0.9595 | 299998 | 1.00 | 0.4217 | 0.2369 | 0.2745 | 0.3433 | 0.7059 | 0.2423 | 0.7108 | 0.6567 | **0.3701 ± 0.0165** | 0.5216 |
| HkT2a | pencere | 1.00 | τ≤1.00 | 0.0000 | 149998 | 0.50 | 0.4299 | 0.2412 | 0.2722 | 0.2701 | 0.9166 | 0.2476 | 0.7022 | 0.6417 | **0.3871 ± 0.0030** | 0.5524 |
| HkT2b | pencere | 1.00 | τ≤1.00 | 0.0000 | 149998 | 0.50 | 0.4328 | 0.2428 | 0.2749 | 0.2652 | 0.8873 | 0.2353 | 0.6981 | 0.6366 | **0.3708 ± 0.0023** | 0.5324 |
| HkT4a | pencere | 1.00 | τ≤1.00 | 0.0000 | 74998 | 0.25 | 0.4297 | 0.2412 | 0.2720 | 0.2208 | 0.9595 | 0.2119 | 0.7028 | 0.6424 | **0.3309 ± 0.0046** | 0.4717 |
| HkT4b | pencere | 1.00 | τ≤1.00 | 0.0000 | 74998 | 0.25 | 0.4300 | 0.2412 | 0.2724 | 0.2157 | 0.9832 | 0.2121 | 0.7024 | 0.6420 | **0.3314 ± 0.0061** | 0.4728 |

‡ = 170'te İNŞA EDİLDİ ve ÖN KAYITLI olarak ölçüldü (§K2 / §K2′).

**Dokuz-bant (lo ≤ 0.80) karşılığı** — 168 §A2.8'in penceresi: Hkeskin
0.4051, son 0.4053, L085 0.3751, **L070 0.3642**, **L060 0.3623**,
**L115 0.4364** (son üçü 169'da yoktu; L070/L060/L115 altı sağlıklı bant
veriyor), K090 0.3709, K070 0.3480, HA4 0.3345, E060 0.3701. Hüküm hiçbir
yerde pencereye bağlı değildir: `c(1.15) > 4/π² > c(0.60)` her iki
pencerede de geçerlidir.

### K0.1 YENİDEN ÜRETİM SINAVI — TAM

Bağımsız kod yolundan (169_k1'in fonksiyonlarını import ederek, ama
kendi tablosuyla) 168 §A2.6/A2.7/A2.8 ve 169 §K1.1/K1.2'nin **her satırı
dört hanede** yeniden üretildi: `c_WX`(5 bant), `c_WX`(9 bant), `g_cal`,
`θ`, `c_ampX` — **on iki gazda tek bir sapma yok** (eşik %0.15).
`θ/θ₀` de birebir: son 1.0290, L085 1.0084, K090 0.9080, K070 0.8038,
HA4 0.7975, E060 0.7945 (168 §A2.7 ile aynı dört hane).
**Yeni:** `θ/θ₀`(L070) = **1.0271** (169 vermemişti).

> **HÜKÜM (K0a).** Defter temiz. 168 ve 169'un `c` sayıları arasında hiçbir
> tutarsızlık yoktur; iki farklı sayı (ör. Hkeskin 0.4035 ↔ 0.4051) yalnız
> **bant penceresidir** (5 bant lo≤0.68 ↔ 9 bant lo≤0.80), üye değil.
> `c_ampX` ile `c_WX` arasındaki fark ise **yalnız `W_amp` bölmesidir**
> (168 §A2.8: `W_amp` çift sayım).

### K0.2 ÇARPAN AYRIŞTIRMASI — **λ EKSENİNİN ADRESİ `W_X`'TİR**

`c_WX = (g_E · g_X² · θ)/W_X` özdeşliğinin çarpanları, Hkeskin'e oranla:

| gaz | λ | σ_ds/σ₀ | σ_X̃/σ₀ | σ_Ĉ/σ₀ | g_E/g₀ | g_X/g₀ | θ/θ₀ | **W_X/W₀** | c/c₀ (çarpım) | c/c₀ (ölçülen) |
|---|---|---|---|---|---|---|---|---|---|---|
| Hkeskin | 1.00 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | **1.0000** | 1.0000 | 1.0000 |
| L085 | 0.85 | 0.9219 | 0.9202 | 0.9207 | 1.0101 | 0.9980 | 1.0084 | **1.0714** | 0.9468 | 0.9459 |
| L070 | 0.70 | 0.8250 | 0.8221 | 0.8244 | 1.0342 | 0.9987 | 1.0271 | **1.1564** | 0.9162 | 0.9145 |
| **L060** | 0.60 | 0.7450 | 0.7416 | 0.7465 | 1.0631 | 1.0022 | 1.0493 | **1.2230** | 0.9162 | 0.9141 |
| **L115** | 1.15 | 1.0655 | 1.0676 | 1.0679 | 0.9985 | 1.0034 | 0.9968 | **0.9388** | 1.0675 | 1.0686 |

İki bulgu, ikisi de yeni:

1. **λ gazı λ ile ölçeklenmiyor.** `σ_ds`, `σ_X̃`, `σ_Ĉ` üçü de aynı
   çarpanla iniyor ama o çarpan λ değil: 0.9219 / 0.8250 / 0.7450
   (λ = 0.85 / 0.70 / 0.60). Kuvvet yasası `σ ∝ λ^p`: `p` = 0.500,
   0.539, **0.576** — yani merdiven genliğini λ ile çarpmak, gazın
   dalgalanmasını kabaca **√λ** kadar değiştiriyor ve üstel λ ile
   yavaşça büyüyor. İnşa denkleminin KENDİSİ doymuş: `rms S' = 1.894λ`
   ile `N̄' = 1.9147` aynı mertebede (λ=1'de tam sınırda; `F' = N̄'+S'`
   işaret değiştirme kesri λ=0.70'te %8.8, λ=0.60'ta %4.4).
   **"Sonlu doyum" en başta burada.**
2. **`c`'nin λ ile değişiminin tamamı `W_X` bölmesindendir.** Payın
   (`KALİB_u2 = g_cal·θ`) λ ile davranışı **TERS yönde**: +%1.4 / +%6.0 /
   **+%12.1** (λ = 0.85 / 0.70 / 0.60) ve λ = 1.15'te **+%0.2** (yani
   λ ≥ 1'de neredeyse SABİT). Payda `W_X` ise +%7.1 / +%15.6 / **+%22.3**
   ve λ = 1.15'te **−%6.1**. `g_X` λ-değişmez (‰3), `θ` −%0.3…+%4.9,
   `g_E` −%0.2…+%6.3 (Gram şişmesi λ ile zayıflıyor — 168 §A1.5'in doğru
   adresi). λ = 0.70 → 0.60'ta pay ve payda **eşit hızda** büyüyor
   (+%5.8 ↔ +%5.8) ⇒ `c` duruyor (§K2.3'ün tabanı); λ = 1.00 → 1.15'te
   pay durup payda düşüyor ⇒ `c` fırlıyor (§K2′).

> **HÜKÜM (K0b — 170'in ilk kazanımı).** `c(λ)` bir **doyum eğrisi değil**,
> bir **payda eğrisidir**. Ölçülen ham kalibrasyon λ düştükçe BÜYÜR;
> `c`'yi düşüren tek şey `W_X` ile bölmektir. Bu, 169 §K1.5'in üye
> sınavının (λ-değişmezliği `√W_X`'i seçiyordu) bağımsız ve **çarpan
> düzeyinde** doğrulanmasıdır. Herhangi bir "bacak-başına doyum" yasası
> bu çarpanın karşısına konmadan sınanamaz.

### K0.3 ÜYE ÜSSÜ ν: İKİ ÖLÇÜT NEDEN UZLAŞMIYOR (`170l_uye_ussu.py`)

169 §K1.5'in "hiçbir üye iki ölçütü birden sağlamıyor" hükmü tek bir
sayıya indirgenebilir. `c = KALİB_u2/W_X` üyesinin doğru olması demek,

```
KALİB_u2 ∝ W_X^ν   ile   ν = 1
```

demektir. ν'nün İKİ bağımsız yönü vardır ve 169 bunları karıştırıyordu:

**(a) BANT yönü** (aynı gaz, τ arttıkça W_X düşer):

| gaz | λ | **ν_bant** | artık rms |
|---|---|---|---|
| L115 | 1.15 | 0.7515 | 0.00371 |
| **Hkeskin** | **1.00** | **0.9302** | **0.00268** |
| L085 | 0.85 | 1.1849 | 0.00366 |
| L070 | 0.70 | 1.5483 | 0.00561 |
| L060 | 0.60 | **1.8515** | 0.00629 |

**(b) λ yönü** (aynı bant, gaz değişir) — komşu λ çiftlerinin yerel eğimi:

| λ çifti | lo=0.52 | 0.56 | 0.60 | 0.64 | 0.68 | **ort** |
|---|---|---|---|---|---|---|
| 1.15 → 1.00 | 0.0758 | 0.0249 | 0.0229 | −0.0535 | −0.1643 | **−0.0188** |
| 1.00 → 0.85 | 0.3760 | 0.2960 | 0.2612 | 0.1679 | 0.0294 | **0.2261** |
| 0.85 → 0.70 | 0.7851 | 0.6849 | 0.6113 | 0.5028 | 0.3626 | **0.5893** |
| **0.70 → 0.60** | 1.2217 | 1.1165 | 1.0237 | 0.9123 | 0.8189 | **1.0186** |

> **HÜKÜM (K0c — 169 §K1d'nin nicel biçimi).** `W_X` üyesi
> * **bant yönünde yalnız λ = 1.00'de doğrudur** (ν = 0.930, artık rms
>   ‰2.7 — güç yasası beş bantta ‰3 içinde tutuyor; 168 §A2.8'in
>   "τ-eğimi ≈ 0" ölçütünün ta kendisi); iki yana da hızla bozulur
>   (λ=1.15'te 0.75, λ=0.60'ta 1.85);
> * **λ yönünde yalnız λ ≈ 0.65 civarında doğrudur** (ν = 1.019),
>   λ = 1 yakınında feci biçimde yanlıştır (ν = 0.226), λ ≥ 1'de ise
>   sıfırdır (ν = −0.019 ⇒ `KALİB` λ-değişmez).
>
> **İki ölçüt hiçbir λ'da uzlaşmıyor** ve `c(λ)`'nın bütün şekli budur:
> `ν = 1` olan yerde `c` düz, `ν < 1` olan yerde `c` λ ile düşer.
> λ = 0.70 → 0.60'ta `ν = 1.019` ⇒ **`c` düz** (0.3690 → 0.3689);
> λ = 1.00 → 0.85'te `ν = 0.226` ⇒ `c` sert düşer (0.4035 → 0.3817);
> λ = 1.15 → 1.00'de `ν ≈ 0` ⇒ `c` `1/W_X` gibi tam hızla değişir
> (0.4312 → 0.4035). **`c(λ)` eğrisi tamamen ν(λ) eğrisidir ve
> `ν(λ)` λ = 1.15 → 0.60 boyunca −0.02'den 1.02'ye TEKDÜZE tırmanır.**

---

## K1 — T(σ) AKTARIM YASASI (`170b_T_yasasi.py`, `log_K1_T.txt`)

### K1.1 TÜRETİM (kalemle)

**Tanım.** Bir bacağın "aktarımı", alanı `F` yerine doyurulmuş karşılığı
`g(F)` konduğunda korelatörün kazandığı çarpandır. Korelatör `F`'ye
doğrusal ve öngörü varyans-eşleştirilmiş olduğu için aktarım, `g(F)` ile
`F`'nin **korelasyon katsayısıdır**:

```
T ≡ ⟨g(F)·F⟩ / (σ_{g(F)} σ_F)                                        (T1)
```

TAM doyumda `g = sgn`: `⟨sgn(F)F⟩ = ⟨|F|⟩ = σ√(2/π)`, `σ_sgn = 1` ⇒

```
T_∞ = √(2/π) = 0.797885                    ← 169 §K2.3'ün bacak-başı sayısı
```

**Kısmi doyum ailesi.** `sgn`'in Gauss yumuşatması `erf`'tir:
`g_α(F) = erf(F/(α σ_F √2))`, `α→0` ⇒ `sgn`. `F ~ N(0,σ²)` için Stein
özdeşliği ve Gauss integralleriyle

```
⟨F g_α⟩ = σ √(2/π) √ρ ,   Var(g_α) = (2/π) arcsin ρ ,   ρ ≡ 1/(1+α²)
⇒   T(ρ) = √( ρ / arcsin ρ )                                          (T2)
    T(ρ→0) = 1   (doyum yok)          T(ρ=1) = √(2/π)   (TAM doyum)
```

**İki bacak ⇒ arcsine ailesi.** İki bacak birden kırpıldığında Gauss
çiftinin işaret korelasyonu `(2/π)arcsin r`, doğrusal karşılığı `r`:

```
T₂(r) = (2/π) arcsin(r)/r ,   T₂(0) = 2/π = T_∞² ,  T₂(1) = 1          (T3)
```

(T2) ve (T3) **aynı arcsine ailesidir** (her ikisi de `ρ/arcsin ρ`
oranıdır); (T3) 169 §K2.3'ün ölçtüğü nesnedir.

**σ ↔ ρ köprüsü (sarılmış Gauss).** Bacağın taşıdığı birikmiş faz
`ϑ_b = 2πτ_b Ĉ ~ N(0, σ_b²)`. Sarılmış Gauss'un koherent (temel-harmonik)
genliği `|⟨e^{iϑ}⟩| = e^{−σ_b²/2}`, koherent GÜÇ payı `e^{−σ_b²}`; geri
kalan güç sarılmıştır ⇒ **doymuş pay**

```
ρ(σ_b) = 1 − e^{−σ_b²}      ⇒   T(σ_b) = √( ρ(σ_b)/arcsin ρ(σ_b) )     (T4)
σ_b → ∞ ⇒ T → √(2/π)  (KALEM'in şartı)      σ_b → 0 ⇒ T → 1
```

**Türetimin kendi öngörüsü (ölçümden ÖNCE, §A5 docstring'inde):** `T(σ)`
σ'da **monoton AZALAN** ve her yerde `≥ √(2/π)`'dir. λ genlikleri
küçültür ⇒ `σ_b` küçülür ⇒ `T` BÜYÜR ⇒ **H-D1, `c`'nin λ azaldıkça
ARTMASINI ve 4/π²'nin ÜSTÜNDE kalmasını ister.**

### K1.2 (T2)'NİN MONTE-CARLO DOĞRULAMASI — TAM (N = 4·10⁶)

| α | ρ = 1/(1+α²) | T (Monte-Carlo) | T(ρ) formül | fark |
|---|---|---|---|---|
| 0.00 (sgn) | 1.00000 | 0.797834 | **0.797885** | −5.1e−5 |
| 0.25 | 0.94118 | 0.876153 | 0.876140 | +1.3e−5 |
| 0.50 | 0.80000 | 0.928827 | 0.928829 | −2e−6 |
| 1.00 | 0.50000 | 0.977194 | 0.977205 | −1.1e−5 |
| 2.00 | 0.20000 | 0.996620 | 0.996622 | −3e−6 |
| 4.00 | 0.05882 | 0.999711 | 0.999711 | −4e−7 |

> **(T2) sayısal olarak MÜHÜRLÜ** (bütün ailede `< 6·10⁻⁵`), `ρ=1`
> ucunda `√(2/π)`'ye yakınsıyor.

### K1.3 T-EĞRİSİ 169'UN ARCSINE ARA-DEĞER NOKTALARINI VURUYOR

KALEM'in K1 şartı: "T-eğrisi 169'un arcsine ara-değer noktalarını
vurmalı". (T3) ile 169 §K2.3'ün beş ölçülen noktası
(`scratchpad/169/K2b_Hkeskin.json`):

| bant lo | n_çizgi | r | **⟨sgn·sgn⟩/r ÖLÇÜLEN** | **T₂(r) = (2/π)arcsin r/r** | fark |
|---|---|---|---|---|---|
| 0.52 | 53 | +0.39118 | 0.67691 | 0.65409 | +3.49% |
| 0.56 | 75 | +0.35468 | 0.66934 | 0.65078 | +2.85% |
| 0.60 | 116 | +0.32068 | 0.65625 | 0.64807 | +1.26% |
| 0.64 | 176 | +0.28116 | 0.66214 | 0.64532 | +2.61% |
| 0.68 | 264 | +0.23418 | 0.64767 | 0.64259 | +0.79% |

> **HÜKÜM (K1a — GEÇTİ).** Türetilen `T₂` eğrisi, 169'un beş ara-değer
> noktasını **%0.8–3.5** içinde veriyor (169'un kendi sapmalarıyla
> BİREBİR aynı sayılar: aynı +%3 mertebeli Gauss'tan-sapma kayması,
> alan-düzeyinde ölçülen `γ₂ < 0` ile tutarlı — §K1.5). **T yasasının
> ŞEKLİ doğrulandı.**

### K1.4 BACAK-BAŞINA σ_b ÖLÇÜMÜ — **TEK BACAKLAR SARILMIYOR**

`σ_b = 2πτ_b σ_Ĉ`; `τ_b` model çizgi genlikleriyle ağırlıklı ortalama τ
(169 §K2.1 ile aynı tanım); sarılma ölçütü `R(τ) = |⟨e^{2πiτĈ}⟩|` ve
Kuiper `√N·V` (düzgünde medyan 1.62, %1 red eşiği 2.00).

**Hkeskin (λ = 1.00, σ_Ĉ = 0.27768):**

| bacak | τ_b | **σ_b (rad)** | ρ = 1−e^{−σ²} | **T(σ_b)** | R(τ_b) | **√N·V** | mod 2π düzgün mü? |
|---|---|---|---|---|---|---|---|
| taşıyıcı h_Q | 0.6200 | 1.0817 | 0.6897 | 0.95196 | 0.54668 | 192.2 | **HAYIR** |
| E (h'_1) | 0.7054 | 1.2307 | 0.7801 | 0.93370 | 0.45401 | 158.7 | **HAYIR** |
| X (y_2) | 0.7166 | 1.2503 | 0.7906 | 0.93118 | 0.44209 | 154.4 | **HAYIR** |
| X (y_3) | 0.7166 | 1.2503 | 0.7906 | 0.93118 | 0.44209 | 154.4 | **HAYIR** |
| **Ç4 TOPLAM** | **2.7587** | **4.8131** | **1.00000** | **0.79789** | 0.00245 | **1.55** | **EVET** |

> **HÜKÜM (K1b — H-D1'in temeli ÇÖKTÜ).** H-D1 "her bacağın kendi
> birikmiş-faz dağılımı" der; ölçüm **tek bacakların sarılMADIĞINI**
> gösteriyor (`σ_b ≈ 1.08–1.25 rad`, `√N·V` = 154–240 ⇒ düzgünlük EZİCİ
> biçimde reddediliyor). Sarılan tek nesne DÖRT frekansın TOPLAMIDIR
> (`σ_Σ = 4.81 rad`, `√N·V = 1.55`) — 169 §K2.1'in bulgusu bacak
> düzeyinde tekrarlandı. `√(2/π)` **bacak başına değil, yalnız toplamda**
> kazanılıyor; bacaklara çarpanlara ayırmak fizikî olarak yanlıştır.

**λ ekseninde (aynı tablo, dört gaz):**

| gaz | λ | σ_Ĉ | σ_E | σ_X | σ_Σ (Ç4) | **Π_b T(σ_b)** | ΠT/(4/π²) | **ölçülen c** |
|---|---|---|---|---|---|---|---|---|
| L115 † | 1.15 | 0.29653 | 1.3189 | 1.3282 | 5.1306 | **0.73806** | 1.8211 | 0.4312 |
| Hkeskin | 1.00 | 0.27768 | 1.2307 | 1.2503 | 4.8131 | **0.77072** | 1.9017 | 0.4035 |
| L085 | 0.85 | 0.25565 | 1.1274 | 1.1594 | 4.4420 | **0.80875** | 1.9955 | 0.3817 |
| L070 | 0.70 | 0.22892 | 1.0027 | 1.0474 | 3.9893 | **0.85336** | 2.1056 | 0.3690 |
| L060 † | 0.60 | 0.20730 | 0.9029 | 0.9549 | 3.6202 | **0.88654** | 2.1875 | 0.3689 |
| son (GERÇEK ζ) | — | 0.27303 | 1.1932 | 1.2266 | 4.7100 | **0.78140** | 1.9280 | 0.4122 |

† = `Π T` ön-kayıt betiklerinden (`ONKAYIT_L060/L115.json`), yani
ölçümden ÖNCE hesaplandı; σ_E/σ_X/σ_Σ aynı marjinallerden.

> **HÜKÜM (K1c).** `Π_b T(σ_b)` ölçülen `c`'nin **1.71–2.40 KATI** ve
> λ ile **TERS yönde** gidiyor: yasa 0.738 → 0.771 → 0.809 → 0.853 →
> 0.887 (λ küçüldükçe artıyor), ölçülen 0.4312 → 0.4035 → 0.3817 →
> 0.3690 → 0.3689 (azalıyor). H-D1 hem büyüklükte (×1.9) hem işarette
> düşüyor. Resmî ön-kayıtlı sınav §K2 ve §K2′'de.

### K1.5 ALAN-DÜZEYİ AKTARIM ve BASIKLIK (T'nin şekil ekseni)

Aynı `T = ⟨|F|⟩/σ_F` alan düzeyinde doğrudan ölçüldü; Edgeworth açılımı
`T = √(2/π)(1 − γ₂/24)` (γ₂ = basıklık fazlası) ile karşılaştırıldı:

| gaz | λ | T_ölç(E) | γ₂(E) | Edgeworth(E) | T_ölç(X) | γ₂(X) | Edgeworth(X) |
|---|---|---|---|---|---|---|---|
| Hkeskin | 1.00 | 0.82961 | −0.5615 | 0.81655 (+1.60%) | 0.82144 | −0.6193 | 0.81847 (+0.36%) |
| L085 | 0.85 | 0.82944 | −0.5586 | 0.81645 (+1.59%) | 0.82287 | −0.6511 | 0.81953 (+0.41%) |
| L070 | 0.70 | 0.82850 | −0.5430 | 0.81594 (+1.54%) | 0.82532 | −0.6943 | 0.82097 (+0.53%) |
| son | — | 0.82680 | −0.5257 | 0.81536 (+1.40%) | 0.82254 | −0.6299 | 0.81883 (+0.45%) |

> İki gözlem:
> 1. **Alan-düzeyi aktarım λ-DEĞİŞMEZ** (E: 0.8296 → 0.8285, ‰1.3;
>    X: 0.8214 → 0.8253, +‰5). λ yalnız genlikleri ölçekler, dağılımın
>    ŞEKLİNİ değiştirmez — bu, `Π T`'nin λ-bağımlılığının **yalnız**
>    `σ_b`'den gelebileceğini, o kanalın da yanlış işarette olduğunu
>    (K1c) bağımsız olarak gösteriyor.
> 2. `γ₂ < 0` (basıklık EKSİĞİ) ⇒ `T > √(2/π)`; ölçülen +%3.8'lik kayma
>    (169 §K2.3) **tam olarak budur** ve Edgeworth ile %0.4–1.6 içinde
>    açıklanıyor. Ancak `γ₂ = −0.56` "etkin katılım sayısı"
>    `n_p = −1.5/γ₂ = 2.7` demektir, oysa genlik katılım oranı
>    `n_p(genlik) = 217`'dir: **çizgiler bağımsız fazlı değildir**
>    (çözülmemiş yüksek-τ çizgileri — 168 §A2.2'nin Gram sorunu).

---

## K2 — ÖN-KAYITLI ÖNGÖRÜ ve YÜZLEŞME

### K2.1 YENİ GAZ: λ = 0.60 (`170_insa60.py`)

164'ün sadakatli çözücüsüyle (`167_insa.main`, kod kopyalanmadı) kuruldu,
**10.4 dk**: ilk-kök hücreleri benzersiz **300000/300000**, `maks|F| =
1.86e−9` (eşik 1e−8), ikiye-bölme adımı 0, **sıralılık TAM**,
`min Δz = 0.1640`, `σ_ds = 0.32135`, `L = 12.02959`.
Merdiven: `rms S' = 1.1364`, `Σa_qω_q = 155.9` (ön-mühürdeki kalem
değerleriyle birebir).

**Ön-mührün ıskası (açık yazılıyor).** `170_insa60.py`'nin docstring'i
`σ_ds(0.60) ∈ [0.327, 0.340]` demişti; ölçülen **0.32135**, aralığın
%1.7 altında. Kuvvet yasası `σ_ds ∝ λ^p`'nin üsteli λ ile hafifçe
kayıyor: `p` = 0.500 (0.85'ten), 0.539 (0.70'ten), **0.576** (0.60'tan).
`σ_X̃(0.60)` için ön-kayıt betiğinin tahmini 0.1803 idi, ölçülen
**0.17951** (−%0.4); `W_X`(bant 0.52) tahmini 0.8285, ölçülen **0.8304**
(+%0.2) — marjinal öngörüler tuttu.

### K2.2 ÖN KAYIT (`170c_onkayit60.py`, saat **20:48:26**)

Korelatöre (J₂ / KALİB_u2) **hiç bakılmadan**, yalnız L060'ın
marjinallerinden (`σ_ds, σ_X̃, σ_Ĉ, ⟨τ⟩_E, ⟨τ⟩_X, W_X(τ)`) altı rakip
yasa sayısallaştırıldı; ölçüm (`170d_olcum60.py = 167_olcum.kos`)
**20:48:52**'de başlatıldı (sıra dosya damgalarıyla kanıtlı).
L060 marjinalleri: `σ_Ĉ = 0.20730`, `⟨τ⟩_E = 0.6932`, `⟨τ⟩_X = 0.7331`.

```
** ÖN KAYIT — λ = 0.60 (ÖLÇÜMDEN ÖNCE YAZILDI, 20:48:26) **
   P1a  H-D1 MUTLAK   Π_b T(σ_b)                [İLK-İLKE] : 0.8865
   P1b  H-D1 ORANLI   c(Hk)·ΠT(0.60)/ΠT(1.00)   [İLK-İLKE] : 0.4645
   P2   H-C1          c = 4/π² evrensel         [İLK-İLKE] : 0.4035
   P3   λ-kuadratik   3 ölçülen λ noktasından   [EMPİRİK]  : 0.3657
   P4   W_X^ν üyesi   ν = 0.823 (log-log EKK)   [EMPİRİK]  : 0.3757
   P5   KALİB λ-değişmez (168 §A3)              [ölmüştü]  : 0.3292
   [betiğin bastığı satır, aynen:]
   "H-D1 (P1b) TEK BAŞINA c(0.60) > c(0.70) diyor;
    P2/P3/P4/P5 hepsi < 0.3690 diyor."
```

*(Ön-kayıt metnindeki bu son cümlenin bir hatası vardır ve düzeltilerek
yazılıyor: **P2 = 0.4035 de 0.3690'ın ÜSTÜNDEDİR**; `c(0.70)`'in altını
söyleyen yasalar yalnız P3, P4, P5'tir. Ayrım yine de keskindir: P1b
`c(0.60) = 0.4645`, yani ölçülen dizinin her noktasından yukarıda.)*

### K2.3 ÖLÇÜM (`170d_olcum60.py`, 3.7 dk) ve YÜZLEŞME (`170g_yuzlesme.py`)

| bant lo | τ_eff | KALİB_u2 | W_X | **c_WX ± jk** |
|---|---|---|---|---|
| 0.52 | 0.5393 | 0.3218 | 0.8304 | 0.3875 ± 0.0016 |
| 0.56 | 0.5793 | 0.3046 | 0.8071 | 0.3774 ± 0.0022 |
| 0.60 | 0.6189 | 0.2922 | 0.7829 | 0.3732 ± 0.0019 |
| 0.64 | 0.6578 | 0.2736 | 0.7576 | 0.3611 ± 0.0026 |
| 0.68 | 0.6981 | 0.2536 | 0.7318 | 0.3465 ± 0.0029 |

```
c(L060) = 0.3689 ± 0.0010(jk) ± 0.0071(bant) = ±0.0072      [5 bant]
```

| ön-kayıtlı öngörü | değer | ölçüm−öngörü | **σ_tot** | hüküm |
|---|---|---|---|---|
| **P1a H-D1 MUTLAK** [ilk-ilke] | 0.8865 | −58.39% | **−71.9σ** | **ÖLDÜ** |
| **P1b H-D1 ORANLI** [ilk-ilke] | 0.4645 | −20.58% | **−13.3σ** | **ÖLDÜ** |
| P2 H-C1 (4/π²) [ilk-ilke] | 0.4035 | −8.59% | **−4.8σ** | ÖLDÜ |
| P3 λ-kuadratik [EMPİRİK] | 0.3657 | +0.87% | **+0.44σ** | **tuttu** |
| P4 W_X^ν üyesi [EMPİRİK] | 0.3757 | −1.81% | **−0.95σ** | **tuttu** |
| P5 KALİB λ-değişmez (168 §A3) | 0.3292 | +12.05% | **+5.5σ** | ÖLDÜ |
| *(referans)* 4/π² = 0.40528 | 0.4053 | −8.98% | **−5.1σ** | — |

**λ dizisinin tamamında H-D1 (oranlı biçim, en iyi şans):**

| gaz | λ | σ_Ĉ | Π_b T(σ_b) | **H-D1 öngörüsü** | **ÖLÇÜLEN c** | ±σ_tot | fark | σ |
|---|---|---|---|---|---|---|---|---|
| Hkeskin | 1.00 | 0.27768 | 0.77072 | 0.4035 (çapa) | 0.4035 | 0.0018 | — | — |
| L085 | 0.85 | 0.25565 | 0.80875 | **0.4235** | 0.3817 | 0.0028 | −9.86% | **−15.0σ** |
| L070 | 0.70 | 0.22892 | 0.85336 | **0.4468** | 0.3690 | 0.0058 | −17.41% | **−13.5σ** |
| **L060** ‡ | 0.60 | 0.20730 | 0.88654 | **0.4642** | **0.3689** | 0.0072 | −20.53% | **−13.2σ** |

‡ = ÖN KAYITLI, ÖRNEKLEM-DIŞI.

> **HÜKÜM (K2 — H-D1 ÖLDÜ, MÜHÜR KURALI UYGULANDI).** KALEM'in kuralı:
> "K2'de tek nokta bile >3σ ıskalarsa H-D1 ÖLÜR ve ölümü aynen
> raporlanır." **Üç noktanın ÜÇÜ de ıskalıyor** (−15.0σ, −13.5σ,
> −13.2σ); mutlak biçim −71.9σ. Kurtarma bahanesi yok:
> * **büyüklük**: `Π_b T(σ_b)` λ=1.00'de bile ölçülenin 1.90 katı;
> * **işaret**: yasa λ düştükçe `c`'nin artmasını ister, ölçülen düşüyor;
> * **temel**: tek bacaklar sarılmıyor (§K1.4), yani `√(2/π)`'nin
>   bacaklara çarpanlara ayrılması ölçülmüş bir olguya aykırı.
>
> **YENİ OLGU (ön-kayıtlı ölçümün asıl kazancı):** `c(0.60) = 0.3689 ±
> 0.0072` ile `c(0.70) = 0.3690 ± 0.0058` **AYNI** (fark +0.0%,
> +0.01σ). λ inişi λ ≈ 0.65 civarında **DURUYOR**. Üç noktadan kurulan
> kuadratik uyumun tepe noktası λ = 0.566 idi ve örneklem-dışı nokta
> onu **+0.44σ**'da vurdu. Yani `c(λ)` tekdüze inen bir eğri değil,
> **λ ≈ 0.6'da tabanı olan bir çukurdur** — 169'un "tekdüze iniyor"
> okuması λ=0.60 ile düzeltildi. §K0.3'ün diliyle: `ν(0.70→0.60) =
> 1.019` ⇒ `c` düz.

### K2.4 KORELATÖR DÜZEYİNDE DOĞRUDAN SINAV (`170e_bacak_aktarim.py`)

169 §K2.3'ün BÖLÜNMÜŞ-MERDİVEN kurulumu (`169_k2b.main`, kod
kopyalanmadı) λ gazlarında da koşturuldu: `G = E·X_a·X_b`, üç BAĞIMSIZ
alan bacağı, `clip(F) = σ_F sgn(F)`, `ρ = Pu2(kırpılmış)/Pu2(tam)`.
H-D1 tutuyorsa ölçülen ÜÇ-bacak aktarımı `c(λ)^{3/4}` gibi DÜŞMELİ.

| gaz | λ | 1 bacak | 2 bacak | **3 bacak (ölçülen)** | **H-D1'in istediği 3 bacak** | fark |
|---|---|---|---|---|---|---|
| Hkeskin | 1.00 | 0.8293 | 0.6586 | **0.5284** | 0.5284 (çapa) | — |
| L085 | 0.85 | 0.8306 | 0.6502 | **0.5285** | 0.5068 | **+4.3%** |
| L070 | 0.70 | 0.8448 | 0.6641 | **0.5542** | 0.4945 | **+12.1%** |
| **L060** | 0.60 | 0.8611 | 0.6894 | **0.5904** | 0.4940 | **+19.5%** |
| *son* (GERÇEK ζ) | — | 0.8321 | 0.6585 | *0.5356* | *0.5369* (c(son)'dan) | *−0.2%* |

> **HÜKÜM (K2c — H-D1'in üçüncü ve bağımsız ölümü).** Korelatörün
> ölçülen kırpma aktarımı λ ile **ARTIYOR** (0.5284 → 0.5285 → 0.5542 →
> 0.5904), H-D1 ise **AZALMASINI** istiyordu (→ 0.5068 → 0.4945 →
> 0.4942). İşaret ters, fark λ=0.60'ta +%19.5.
>
> **Aynı ölçüm T-yasasını KISMEN AKLIYOR (dürüstlük):** kırpma
> aktarımının λ ile ARTMASI tam olarak `T(σ)`'nın öngördüğü yöndür
> (σ_Ĉ düşer ⇒ T artar). Dört-bacak uzatmasında ölçülen artış
> 0.4730 → 0.4760 → 0.5092 → 0.5497 (**+%16.2**, λ = 1.00 → 0.60),
> T-yasasının öngördüğü `Π T` artışı 0.7707 → 0.8088 → 0.8534 →
> 0.8865 (**+%15.0**) — **aynı işaret, neredeyse aynı büyüklük.**
> Yani `T(σ)` KIRPMA AKTARIMINI doğru tarif ediyor; yanlış olan,
> **`c`'nin o aktarım olduğu** varsayımıdır. `c` λ ile DÜŞÜYOR, kırpma
> aktarımı ARTIYOR: **`c` bir bacak-aktarımı çarpımı DEĞİLDİR.**
>
> **AMA (K3'ün tohumu):** aynı λ'daki İKİ gaz karşılaştırıldığında
> (`son` ↔ `Hkeskin`) ölçülen kırpma aktarımı farkı `c` farkını
> **doğru veriyor** (satırın son sütunu: −%0.2). Yani kırpma aktarımı
> `c`'nin **λ-SABİT kesitinde** doğru bir tahmincidir, λ boyunca
> değildir — çünkü λ, kırpma fiziğinden bağımsız olarak `W_X`
> paydasını da değiştiriyor (§K0.2).

---

## K3 — TAÇ: GERÇEK GAZIN FAZLASI ve H-D2 (`170m_k3.py`)

**Ölçülen fazla.** `son` (gerçek ζ) ile `Hkeskin` (aynı ölçekte sadakatli
sentetik) arasındaki fark:

```
c(son) − c(Hkeskin) = 0.4122 − 0.4035 = +0.0087 ± 0.0029  =  +3.0σ
(169'un '+3.06σ'sı 4/π²'ye göreydi: c(son) − 4/π² = +0.0069 = +3.06σ)
```

### K3.1 162'nin %29'u σ_b'DE GÖRÜNMÜYOR — ve 162 bunu ZATEN SÖYLÜYOR

H-D2'nin zinciri: kilitli asal fazlar η'yı %29 sessizleştirir ⇒ birikmiş
faz σ_b'ler küçülür ⇒ T büyür ⇒ c yükselir. **İlk halka kopuk.**

| büyüklük | son/Hkeskin |
|---|---|
| σ_ds | 0.9488 (−%5.1) |
| σ_X̃ | 0.9662 (−%3.4) |
| **σ_Ĉ** (σ_b'yi belirleyen) | **0.9833 (−%1.7)** |

162 §6'nın kendi ölçümü bunu zaten söylüyordu: yıkıcı girişim oranı
`Σ|c_η|²/2 ÷ Var(η) = 1.288` (η kanalı, %29 sessizleşme) ama
**`Ĉ` kanalında aynı oran 1.026 — girişim YOK.** `σ_b = 2πτ_bσ_Ĉ`
tamamen `Ĉ` kanalından gelir; dolayısıyla **H-D2'nin mekanizması
σ_b üzerinden çalışamaz.** Ölçülen −%1.7, 162'nin Ĉ-kanalı ölçümüyle
(≈%0) tutarlı, %29 ile değil.

### K3.2 FAZLA NE KADAR KAPANIYOR — ÜÇ AYRI KURAL

| kural | öngörülen c(son) | ölçülen | **artık** |
|---|---|---|---|
| (yok) — sadakatli ikizle doğrudan fark | 0.4035 | 0.4122 | **+3.00σ** |
| **(b)** T-yasası: `c(Hk)·ΠT(son)/ΠT(Hk)` = 0.4035·1.0139 | **0.4091** | 0.4122 | **+1.07σ** |
| **(e)** **ÖLÇÜLEN** 3-bacak kırpma aktarımı, `(ρ₃/ρ₃⁰)^{4/3}` = 1.0182 | **0.4110** | 0.4122 | **+0.43σ** |
| (c) λ ekseninden ölçülen `dc/dσ_Ĉ` doğrusu (`c = 0.4911σ_Ĉ + 0.2618`) | 0.3958 | 0.4122 | **+7.22σ** |

`son`un bacak kırpma aktarımları (`169_k2b.main` ile, Hkeskin'le aynı
kurulum): 1 bacak **0.8321** (Hk: 0.8293), 2 bacak **0.6585** (0.6586),
3 bacak **0.5356** (0.5284).

> **HÜKÜM (K3 — H-D2'nin İŞARETİ ve BÜYÜKLÜĞÜ TUTTU, ama MEKANİZMASI
> DEĞİL.** Gerçek gazın sadakatli ikizine göre +3.0σ'lık fazlası,
> **ölçülen bacak kırpma aktarımıyla %86 kapanıyor** (+3.00σ → +0.43σ);
> T-yasasının `σ_Ĉ` üzerinden verdiği kapanış %64'tür (+1.07σ). İkisi de
> **doğru yönde**. Ama H-D2'nin zinciri (162'nin %29'u ⇒ küçük σ_b)
> **kopuktur**: σ_Ĉ yalnız −%1.7 küçülüyor ve 162 zaten "girişim η
> kanalında, Ĉ kanalında YOK (1.026)" demişti.
>
> **En sert kısıt burada:** aynı "ölçülen bacak aktarımı" kuralı λ
> ekseninde **çöküyor** (L085 −6.7σ, L070 −10.1σ, L060 −13.4σ). Yani
> kırpma aktarımı, **aynı λ'daki iki gaz arasındaki** farkı doğru
> veriyor, ama **λ'yı değiştirince** vermiyor — çünkü λ'yı değiştirmek
> `W_X` paydasını da değiştiriyor (§K0.2) ve o kanal kırpma fiziğinden
> bağımsızdır.

### K3.3 BÜTÜN λ DİZİSİNİN TOPLU DEFTERİ

| gaz | λ | σ_Ĉ | Π_b T(σ_b) | ΠT/ΠT₀ | KALİB | W_X | **c_WX** | ±σ_tot |
|---|---|---|---|---|---|---|---|---|
| L115 | 1.15 | 0.29653 | 0.73806 | 0.9576 | 0.2577 | 0.6003 | **0.4312** | 0.0053 |
| Hkeskin | 1.00 | 0.27768 | 0.77072 | 1.0000 | 0.2572 | 0.6394 | **0.4035** | 0.0018 |
| L085 | 0.85 | 0.25565 | 0.80875 | 1.0493 | 0.2609 | 0.6851 | **0.3817** | 0.0028 |
| L070 | 0.70 | 0.22892 | 0.85336 | 1.1072 | 0.2725 | 0.7395 | **0.3690** | 0.0058 |
| L060 | 0.60 | 0.20730 | 0.88654 | 1.1503 | 0.2882 | 0.7820 | **0.3689** | 0.0072 |
| **son** (GERÇEK ζ) | — | 0.27303 | 0.78140 | 1.0139 | 0.2716 | 0.6607 | **0.4122** | 0.0023 |

**Kırpma aktarımları (`170e`, tam beş gaz):**

| gaz | 1 bacak | 2 bacak | 3 bacak | (1 bacak)⁴ | ölçülen c |
|---|---|---|---|---|---|
| Hkeskin | 0.8293 | 0.6586 | 0.5284 | 0.4730 | 0.4035 |
| L085 | 0.8306 | 0.6502 | 0.5285 | 0.4760 | 0.3817 |
| L070 | 0.8448 | 0.6641 | 0.5542 | 0.5092 | 0.3690 |
| L060 | 0.8611 | 0.6894 | 0.5904 | 0.5497 | 0.3689 |
| **son** | 0.8321 | 0.6585 | 0.5356 | 0.4794 | 0.4122 |
| *Gauss* | *0.7979* | *0.6366* | *0.5079* | *0.4053* | *(4/π²)* |

---

## K2′ — EK KAPI (kendi kapım): λ = 1.15, DERİN UÇ

K2, `c(λ)`'nın λ ≤ 0.70'te bir **tabana** oturduğunu gösterdi. Geriye
şu ayrım kaldı: `c` λ ile ARTARKEN neye yaklaşıyor? İki okuma:

* **Y1** — derin (büyük λ) doyum limitine, yani `4/π²`'ye **ALTTAN**
  yaklaşıyor ⇒ `c(1.15) ∈ (0.4035, 0.4053]`;
* **Y2** — H-D1'in T-yasası: λ büyüyünce `σ_Ĉ` büyür ⇒ `Π T` KÜÇÜLÜR ⇒
  `c(1.15) < 0.4035`.

**Gaz (`170i_insa115.py`, 11.2 dk).** 167'de yarım kalan L115 inşası
tamamlandı: ilk-kök hücreleri benzersiz 300000/300000, `maks|F| =
1.86e−9`, ikiye-bölme adımı 0, **sıralılık TAM**, `min Δz = 0.0934`,
`σ_ds = 0.45958`. Ön-mühür `σ_ds(1.15) ≈ 0.4646` demişti — **%1.1 iska**
(doğru yönde, kuvvet üsteli beklenenden küçük). Marjinaller:
`σ_X̃ = 0.25839`, `σ_Ĉ = 0.29653`, `⟨τ⟩_E = 0.7079`, `⟨τ⟩_X = 0.7129`.

**ÖN KAYIT (`170j_onkayit115.py`, saat 21:16:11; ölçüm 21:16:21):**

```
   Q1  H-D1 (T-yasası, oranlı)         : c(L115) = 0.3865   ← TEK 'AŞAĞI'
   Q2  H-C1 (c = 4/π²)                 : c(L115) = 0.4053
   Q3  derin limite ALTTAN yaklaşma    : c(L115) ∈ (0.4035, 0.4053]
   Q4  λ-kuadratik (4 nokta, EMPİRİK)  : c(L115) = 0.4376
   Q5  W_X^ν üyesi (ν = 0.830, EMPİRİK): c(L115) = 0.3917
   Q6  KALİB λ-değişmez (168 §A3)      : c(L115) = 0.4303
   çapa: c(1.00) = 0.4035 ± 0.0018
```

**ÖLÇÜM (`170k_olcum115.py`, 4.2 dk) ve yüzleşme (`170n_yuzlesme115.py`):**

| bant lo | τ_eff | KALİB_u2 | W_X | **c_WX ± jk** |
|---|---|---|---|---|
| 0.52 | 0.5393 | 0.2852 | 0.6794 | 0.4198 ± 0.0033 |
| 0.56 | 0.5792 | 0.2702 | 0.6401 | 0.4222 ± 0.0033 |
| 0.60 | 0.6189 | 0.2582 | 0.6005 | 0.4300 ± 0.0016 |
| 0.64 | 0.6582 | 0.2449 | 0.5604 | 0.4369 ± 0.0030 |
| 0.68 | 0.6983 | 0.2334 | 0.5210 | 0.4479 ± 0.0037 |

```
c(L115) = 0.4312 ± 0.0014(jk) ± 0.0051(bant) = ±0.0053     [5 bant]
```

| ön-kayıtlı öngörü | değer | ölçüm−öngörü | **σ_tot** | hüküm |
|---|---|---|---|---|
| **Q1 H-D1 (T-yasası)** [ilk-ilke] | 0.3865 | +11.59% | **+8.5σ** | **ÖLDÜ** |
| **Q2 H-C1 (4/π²)** [ilk-ilke] | 0.4053 | +6.40% | **+4.9σ** | **ÖLDÜ** |
| **Q3 "derin limite ALTTAN"** [ilk-ilke] | (0.4035, 0.4053] | — | ölçüm **ARALIĞIN ÜSTÜNDE** | **ÖLDÜ** |
| Q4 λ-kuadratik [EMPİRİK] | 0.4376 | −1.45% | **−1.20σ** | tuttu |
| Q5 W_X^ν üyesi [EMPİRİK] | 0.3917 | +10.10% | **+7.5σ** | ÖLDÜ |
| **Q6 KALİB λ-değişmez** (168 §A3) | 0.4303 | +0.22% | **+0.18σ** | **TUTTU** |

**Q6'nın içeriği doğrudan:** `KALİB_u2` λ = 1.00 → 1.15'te
−0.36 / −0.14 / −0.15 / +0.39 / +1.35 % (bant bant), ortalama **+%0.22**;
`W_X` aynı aralıkta **−%6.21**. Yerel üye üssü
**ν(1.00→1.15) = −0.034 ≈ 0**.

> **HÜKÜM (K2′ — ÜÇ İLK-İLKE YASA DA ÖLDÜ; 4/π² BİR LİMİT DEĞİL).**
> `c(1.15) = 0.4312 ± 0.0053`, `4/π² = 0.40528`'in **+%6.4 ÜSTÜNDE
> (+4.9σ)**. Yani `c` λ ile artmaya devam ediyor ve `4/π²`'yi λ ≈ 1.00'de
> **hiçbir özellik göstermeden kesip geçiyor**. Bu, "λ = 1'de gördüğümüz
> 4/π² uyumu bir LİMİT'tir" okumasını (Y1/Q3) örneklem-dışı biçimde
> öldürür: **`c = 4/π²` bir yakınsama noktası değil, λ = 1.00'in
> tesadüfüdür.** 169'un `c₀ = 0.4057 ± 0.0016 = 4/π² + 0.23σ` kesişimi
> de bu ışıkta yeniden okunmalıdır (o kesişim λ = 1.00 gazlarından
> kuruluydu).
>
> **Tersine, 168 §A3'ün "KALİB_u2 λ-DEĞİŞMEZ" yasası λ ≥ 1'de DİRİ**
> (+%0.2, +0.18σ) — 169 onu λ = 0.70'te +3.6σ ile reddetmişti. Şimdi
> tablo tam: **`ν(λ)` λ = 1.15 → 0.60 boyunca −0.03 → 0.23 → 0.59 →
> 1.02 diye tırmanıyor**; yasa λ ≥ 1'de doğru, λ ≤ 0.7'de yanlış.

---

## BONUS (`170f_beta.py`, `log_BONUS.txt`)

### BONUS-i — SONLU-DOYUM DÜZELTMESİNİN KAPALI BİÇİMİ (T5)

`ρ = 1 − ε`, `ε = e^{−σ²}`; `arcsin(1−ε) = π/2 − √(2ε) + O(ε^{3/2})` ⇒

```
T² = (1−ε)/(π/2 − √(2ε)) = (2/π)[1 + (2√2/π)√ε − ε + O(ε^{3/2})]
⇒  T(σ) = √(2/π)·[ 1 + (√2/π) e^{−σ²/2} + O(e^{−σ²}) ],   √2/π = 0.45016   (T5)
```

Yani **sonlu-doyum düzeltmesi, sarılmış Gauss'un KOHERENT genliği
`e^{−σ²/2}` ile doğrusaldır** ve katsayısı `√2/π`'dir. Dört bacak için
`c ≈ (4/π²)[1 + (√2/π) Σ_b e^{−σ_b²/2}]` (T6). Açılımın kalitesi:

| σ | 0.80 | 1.00 | 1.25 | 1.50 | 2.00 | 3.00 | 4.80 |
|---|---|---|---|---|---|---|---|
| T(σ) tam | 0.9798 | 0.9611 | 0.9312 | 0.8987 | 0.8437 | 0.8019 | 0.79789 |
| (T5) | 1.0587 | 1.0157 | 0.9623 | 0.9145 | 0.8465 | 0.8019 | 0.79789 |
| fark | +8.0% | +5.7% | +3.3% | +1.8% | +0.33% | +0.002% | +0.000% |

> (T5) `σ ≥ 2` için ‰3, `σ ≥ 3` için 2·10⁻⁵ içinde. Ölçülen bacak
> σ'ları 1.0–1.3 olduğu için orada açılım %3–6 hatalıdır (tam form
> kullanıldı); Ç4 toplamında (σ_Σ ≈ 4–4.8) makinede ayırt edilemez.

### BONUS-ii — β = 0.2175 T'DEN TÜRÜYOR MU? **HAYIR**

| gaz | φ | σ_Ĉ | σ_Ĉ/σ₀ | Π_b T(σ_b) | **ΠT/ΠT₀** | **θ/θ₀ ÖLÇÜLEN** | 1 − 0.2175φ |
|---|---|---|---|---|---|---|---|
| Hkeskin | 0.0000 | 0.27768 | 1.0000 | 0.77072 | 1.0000 | 1.0000 | 1.0000 |
| K090 | 0.4153 | 0.27872 | 1.0037 | 0.76885 | **0.9976** | **0.9080** | 0.9097 |
| K070 | 0.9270 | 0.27663 | 0.9962 | 0.78496 | **1.0185** | **0.8038** | 0.7984 |
| HA4 | 0.9098 | 0.27547 | 0.9921 | 0.78575 | **1.0195** | **0.7975** | 0.8021 |
| E060 | 0.9595 | 0.27450 | 0.9886 | 0.78667 | **1.0207** | **0.7945** | 0.7913 |

```
(T6)'nın ürettiği eğim:  β_T = −0.0245     ölçülen β = +0.2175
                          oran = −0.113   (İŞARET TERS, mertebe 1/9)
```

> **HÜKÜM (BONUS-ii DÜŞTÜ).** Kesim ekseninde `σ_Ĉ` neredeyse sabittir
> (±%1.1), dolayısıyla `Π_b T(σ_b)` de sabittir (±%2) — üstelik φ
> arttıkça **artıyor**, `θ` ise **azalıyor**. **T çerçevesi β'yı
> üretemez**: ne büyüklüğü (1/9) ne de işareti. `β`'nın türetimi 168 →
> 169 → 170 boyunca açık kalmıştır ve T-yasası onun adresi DEĞİLDİR.

### BONUS-iii — T/4 ARTIĞI: (T6) ERİTMİYOR; ama artığın ADRESİ DARALDI

| gaz | N | σ_Ĉ/σ₀ | g_E/g₀ | g_X/g₀ | (g_X/g₀)² | g_cal/g₀ | (1+2n/N)₀/(1+2n/N) | **artık** | θ/θ₀ |
|---|---|---|---|---|---|---|---|---|---|
| Hkeskin | 299998 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| HkT2a | 149998 | 0.9803 | 0.9904 | 0.9706 | 0.9421 | 0.9330 | 0.9498 | 0.9824 | 1.0317 |
| HkT2b | 149998 | 0.9899 | 0.9841 | 0.9648 | 0.9308 | 0.9160 | 0.9432 | 0.9712 | 0.9987 |
| HkT4a | 74998 | 0.9795 | **0.9549** | **0.8938** | 0.7988 | 0.7628 | 0.8629 | **0.8840** | 1.0800 |
| HkT4b | 74998 | 0.9810 | **0.9514** | **0.8851** | 0.7833 | 0.7452 | 0.8575 | **0.8691** | 1.1066 |

> **HÜKÜM (BONUS-iii DÜŞTÜ, ama bir adres kazanıldı).** `σ_Ĉ` pencereyle
> ‰20 içinde sabittir ⇒ `Π_b T(σ_b)` sabittir ⇒ **(T6) T/4 artığına
> hiçbir katkı vermiyor.** YENİ BİLGİ: `(1+2n/N)^{−1}` sızıntı
> düzeltmesinden artan −%12/−%13'lük açık **neredeyse tamamen `g_X`
> bacağındadır**: T/4'te `g_E` yalnız −%4.5/−%4.9 düşerken `g_X`
> −%10.6/−%11.5 düşüyor (kareye girdiği için `g_cal`'daki payı iki
> katı). Yani T/4 borcunun adresi "genel sonlu-örneklem sızıntısı"
> değil, **X̃0 kanalının çizgi-uyum sızıntısıdır.**

---

## 4. HÜKÜM

| kapı | sonuç | dayanak |
|---|---|---|
| **K0 — MUHASEBE** | **GEÇTİ (yeni bulguyla)** | 12 gazın `c_WX`(5 ve 9 bant), `g_cal`, `θ`, `c_ampX` değerleri 168/169 ile **dört hanede** aynı, sapma yok. YENİ: λ ekseninin adresi `W_X` PAYDASI (§K0.2); üye üssü `ν` bant yönünde yalnız λ=1.00'de, λ yönünde yalnız λ≈0.65'te 1'dir (§K0.3). |
| **K1 — T(σ) TÜRETİMİ** | **GEÇTİ (yasa), ÇÖKTÜ (uygulama)** | `T(ρ)=√(ρ/arcsin ρ)`, `ρ=1−e^{−σ²}` türetildi; Monte-Carlo `<6·10⁻⁵`; iki-bacak biçimi 169'un beş arcsine noktasını %0.8–3.5'te veriyor. AMA `σ_b` ölçüldü ve **tek bacaklar sarılmıyor** (Kuiper 154–240 ↔ eşik 2.00), `Π_b T(σ_b)` ölçülen `c`'nin 1.90–2.11 katı. |
| **K2 — ÖN-KAYITLI ÖNGÖRÜ** | **H-D1 ÖLDÜ** | λ=0.85/0.70/0.60'ta **−15.0σ / −13.5σ / −13.2σ** (mutlak biçim −71.9σ); λ=0.60 ÖN KAYITLI ve ÖRNEKLEM-DIŞI. Korelatör düzeyinde bağımsız üçüncü ölüm (§K2.4). YENİ OLGU: `c(0.60) = c(0.70)` — λ inişi duruyor. |
| **K3 — TAÇ (H-D2)** | **İŞARET+BÜYÜKLÜK TUTTU, MEKANİZMA DÜŞTÜ** | Gerçek gazın ikizine göre +3.00σ fazlası, ÖLÇÜLEN bacak kırpma aktarımıyla **%86 kapanıyor** (+3.00σ → **+0.43σ**); T-yasasının σ_Ĉ üzerinden kapanışı %64 (+1.07σ). Ama H-D2'nin zinciri kopuk: `σ_Ĉ(son)/σ_Ĉ(Hk) = 0.983` (−%1.7), 162'nin %29'u değil — ve 162 zaten "girişim η'da, Ĉ'de YOK (1.026)" demişti. Aynı kural λ ekseninde çöküyor (−6.7σ…−13.4σ). |
| **K2′ — EK KAPI λ=1.15** | **ÜÇ İLK-İLKE YASA DA ÖLDÜ** | ÖN KAYITLI: `c(1.15) = 0.4312 ± 0.0053`. H-D1 (0.3865) **+8.5σ**; H-C1 4/π² **+4.9σ**; "derin limite alttan yaklaşma" — ölçüm aralığın ÜSTÜNDE. **`c` 4/π²'yi λ=1.00'de özelliksiz kesip geçiyor ⇒ 4/π² bir LİMİT DEĞİL.** 168 §A3'ün "KALİB λ-değişmez" yasası λ≥1'de DİRİ (+0.18σ). |
| **BONUS — β, T/4** | **İKİSİ DE DÜŞTÜ (adres daraldı)** | (T5) kapalı biçimi türetildi: `T = √(2/π)[1+(√2/π)e^{−σ²/2}]`, σ≥2'de ‰3. β: T'nin ürettiği eğim **−0.0245** ↔ ölçülen **+0.2175** (işaret ters, mertebe 1/9). T/4: σ_Ĉ pencereyle sabit ⇒ katkı yok; YENİ: artığın neredeyse tamamı **`g_X` bacağında** (T/4'te g_E −%4.7, g_X −%11). |

### MÜHÜR KURALININ UYGULANMASI (açık)

KALEM: *"K1 + K2 zorunlu; K3 taç. K2'de tek nokta bile >3σ ıskalarsa
H-D1 ÖLÜR ve ölümü aynen raporlanır (kurtarma bahanesi yok)."*

**K2'de ÜÇ noktanın ÜÇÜ de ıskaladı (−15.0σ / −13.5σ / −13.2σ), ek
kapıda dördüncü nokta da (+8.5σ). H-D1 ÖLDÜ ve ölümü aynen
raporlanmıştır.** Hiçbir kurtarma denenmemiş, hiçbir serbest parametre
eklenmemiştir. `T(σ)` yasasının KENDİSİ (K1'in türetimi + Monte-Carlo +
arcsine noktaları + korelatör kırpma aktarımının λ-eğilimi) doğrudur ve
ayrıca mühürlenmiştir; ölen şey, **`c`'nin o yasanın bacak çarpımı
olduğu** varsayımıdır.

### 0'ın cevabı — TEK CÜMLELİK HÜKÜM

> **`c` BİR BACAK-AKTARIMI ÇARPIMI DEĞİL, BİR BÖLME ARTIĞIDIR.**
> Bacak-başına sarılmış-aktarım yasası `T(σ) = √(ρ/arcsin ρ)`,
> `ρ = 1−e^{−σ²}` türetildi, sayısal olarak mühürlendi ve korelatörün
> ÖLÇÜLEN kırpma aktarımını λ boyunca **doğru işaret ve doğru
> mertebede** (+%16.2 ↔ +%15.0) tarif ettiği gösterildi — ama `c` o
> aktarımın tersine gidiyor. `c(λ)`'nın bütün şekli tek bir sayının,
> **üye üssü `ν(λ) = d log KALİB_u2 / d log W_X`'in** şeklidir
> (λ = 1.15 → 0.60 boyunca −0.02 → 0.23 → 0.59 → 1.02) ve iki ÖN
> KAYITLI örneklem-dışı nokta bunu mühürlüyor: λ = 0.60'ta `c` düz
> kalıyor (ν = 1.02), λ = 1.15'te `4/π²`'nin **+4.9σ üstüne** çıkıyor
> (ν ≈ 0). **Dolayısıyla `4/π²` bir doyum limiti değil, λ = 1.00'in
> tesadüfüdür; ejderha bir doyum eğrisi değil, bir NORMALİZASYON
> eğrisidir.** Tek olumlu taç: gerçek ζ gazının sadakatli ikizine göre
> +3.00σ'lık fazlası, o gazın ölçülen bacak kırpma aktarımıyla **%86
> kapanıyor** (artık +0.43σ) — yani kırpma resmi λ-SABİT kesitte
> DOĞRU, λ boyunca DEĞİL.

---

## 5. DENETİM ve DÜRÜSTLÜK NOTLARI

* **Yeniden üretim tam.** `170a_muhasebe.py` 168 §A2.6/A2.7/A2.8 ve
  169 §K1.1/K1.2'nin on iki gazdaki bütün sayılarını (`c_WX` 5 bant ve
  9 bant, `g_cal`, `θ`, `c_ampX`, `θ/θ₀`) **dört hanede** yeniden üretti;
  eşik %0.15, sapma sıfır (`log_K0.txt`).
* **Ölçüm/çözüm parçası kopyalanmadı.** İnşa `167_insa.main`
  (→ `164_insa.coz_sadakatli`), ölçüm `167_olcum.kos`
  (→ `166_T1.bant_agg/_jk`, `165_cekirdek.Model165/tayf_s/sentez/s_den_J`,
  `163_cekirdek.olc_cizgi/bant_adaylari`, `166_bacak.karakteristik`),
  bacak aktarımı `169_k2b.main`, bant seçimi/jackknife `169_k1.*`,
  Kuiper `169_k2.kuiper` — hepsi **import**.
* **ÖN KAYIT sırası dosya damgalarıyla kanıtlı.**
  `ONKAYIT_L060.json` **20:48:26**, `170d_olcum60.py` **20:48:52**'de
  başlatıldı (`log_C_L060.txt`); `ONKAYIT_L115.json` **21:16:11**,
  `170k_olcum115.py` **21:16:21**'de başlatıldı.
  Ön kayıtlar korelatöre (J₂ / KALİB_u2 / c) HİÇ bakmaz; yalnız gazın
  marjinallerini (`σ_ds, σ_X̃, σ_Ĉ, ⟨τ⟩_E, ⟨τ⟩_X, W_X(τ)`) ve ÖNCEDEN
  ölçülmüş gazların bant tablolarını kullanır (169 §K3.2'nin protokolü).
* **Ön-mühürlerin ıskaları gizlenmedi.** (i) `170_insa60.py`'nin
  `σ_ds(0.60) ∈ [0.327, 0.340]` tahmini ISKA (ölçülen 0.32135, −%1.7);
  (ii) `170i_insa115.py`'nin `σ_ds(1.15) ≈ 0.4646` tahmini ISKA
  (ölçülen 0.45958, −%1.1); (iii) `170e_bacak_aktarim.py`'nin
  "korelatör aktarımı λ ile ~sabit kalır (|Δ| < %2)" beklentisi
  λ = 0.85'te tuttu (‰0.2) ama λ = 0.70 ve 0.60'ta ISKA (+%4.9, +%11.7);
  (iv) aynı betiğin "`son`un fazlası bacak aktarımıyla KAPANMAMALI"
  beklentisi de ISKA — kapandı (§K3.2, +0.43σ). Dördü de ilgili
  bölümlerde yazılıdır. Buna karşılık `170c`'nin marjinal öngörüleri
  (σ_X̃(0.60) = 0.1803 ↔ 0.17951; W_X = 0.8285 ↔ 0.8304) tuttu.
* **Ön-kayıt metninde bir hata vardı ve düzeltilerek yazıldı:**
  `170c_onkayit60.py`'nin son satırı "P2/P3/P4/P5 hepsi < 0.3690" der;
  P2 = 0.4035 aslında ÜSTTEDİR (§K2.2'nin notu). Hata yalnız o cümlenin
  özetindedir; öngörü SAYILARI etkilenmez.
* **İki betik hatası yapıldı ve düzeltildi.** (i) `170e`'nin ilk
  sürümünde bir parantez eksikti (SyntaxError, hiçbir şey koşmadı; 3 dk
  kayıp). (ii) `170b`'nin ilk `_kesim` koşusu `L060` künyesi
  `167_ortak.KUNYE`'de olmadığı için `KeyError` ile bitti — DÖRT kesim
  gazının ekran çıktısı geçerlidir ama JSON yazılmamıştı; künye araması
  güvenli hâle getirilip dört kesim gazıyla yeniden koşuldu. Rapordaki
  kesim sayıları ikinci koşudandır.
* **`son` (gerçek ζ) gazının λ'sı yoktur** (merdiveni inşa edilmemiştir);
  K0 tablosunda `λ = —` yazılıdır ve λ-ekseni karşılaştırmalarına
  KATILMAZ; yalnız K3'te Hkeskin'e karşı konumlandırılır.
* **σ tanımı 169 ile aynı**: `σ_jk` bant-içi 8-grup jackknife'ların gaz
  düzeyine taşınması, `σ_bant` bantlar arası saçılımın standart hatası,
  `σ_tot = hipot(σ_jk, σ_bant)`; bütün hükümler `σ_tot` ile.
* **Bant penceresi** her yerde ortak: `lo ∈ [0.52, 0.68]`, `SNR ≥ 3`,
  `R_bant ≥ 0.98`, `τ_eff < 0.85` (beş bant). Dokuz-bant (lo ≤ 0.80)
  sayıları K0'da ayrıca verildi.
* **Uydurma yok.** Bu raporda hesaplanmamış tek sayı yoktur; her tablo
  `scratchpad/170`'teki bir log/JSON'a karşılık gelir (`K0.json`,
  `UYE_USSU.json`, `K1_T.json`, `K1_T_kesim.json`, `ONKAYIT_L060.json`,
  `K2_yuzlesme.json`, `K3.json`, `BONUS.json`, `ONKAYIT_L115.json`,
  `K2p_yuzlesme.json`, `log_*.txt`) ya da `scratchpad/167`
  (`C_L060.json`, `C_L115.json`, `insa_L060/L115.json`,
  `z_L060/L115.npy`) / `scratchpad/169` (`K2b_*.json`) dosyalarına.
* **169'un `K2b_Hkeskin.json`'ı yeniden koşulmadı**, aynen okundu;
  L085/L070/L060/son için aynı betik (`169_k2b.main`) 170'te koşuldu ve
  çıktıları 169'un dizinine yazıldı (aynı ad şeması).
* **Git'e dokunulmadı** (commit/add/push yok).

---

## 6. NE KAZANILDI, NE KALDI

**KAZANILDI (170'in net çıktısı)**

1. **`T(σ) = √(ρ/arcsin ρ)`, `ρ = 1−e^{−σ²}` türetildi ve mühürlendi**
   (Monte-Carlo `<6·10⁻⁵`; iki-bacak biçimi 169'un beş arcsine noktasını
   %0.8–3.5'te veriyor); derin açılımı **`T = √(2/π)[1+(√2/π)e^{−σ²/2}]`**
   (T5) — istenen "sonlu-doyum düzeltme yasası"nın kapalı biçimi.
2. **H-D1 ÖN KAYITLI olarak öldü** (−13…−15σ, üç λ noktasında; korelatör
   düzeyinde de bağımsız olarak), ve ölümünün NEDENİ ölçüldü: **tek
   bacaklar sarılmıyor**, yalnız dörtlü toplam sarılıyor.
3. **λ ekseninin ADRESİ bulundu: `W_X` paydası.** `c(λ)` = `ν(λ)` eğrisi
   (§K0.3). Ham `KALİB_u2` λ ≥ 1'de λ-DEĞİŞMEZ (168 §A3 orada diri),
   λ ≤ 0.7'de `W_X` ile birlikte gidiyor (ν → 1).
4. **`c = 4/π²` bir limit değil, λ = 1.00'in tesadüfü** — ÖN KAYITLI
   λ = 1.15 noktası 4/π²'nin +4.9σ üstünde (§K2′).
5. **`c(λ)`'nın λ ≈ 0.6'daki tabanı ölçüldü** (`c(0.60) ≡ c(0.70)`,
   +0.01σ) ve ön kayıtla öngörüldü (+0.44σ).
6. **K3 POZİTİF:** gerçek gazın +3.00σ fazlasının %86'sı, o gazın
   ölçülen bacak kırpma aktarımıyla kapanıyor (artık +0.43σ).
7. **λ gazının kendisi doymuş:** `σ_ds ∝ λ^p`, p = 0.50 → 0.58 (λ ile
   yavaşça artıyor); λ ile ölçeklemek gazı λ ile ölçeklemiyor.
8. **İki yeni gaz** (L060, L115) inşa edildi ve ölçüldü; L115 167'den
   beri yarım kalan borçtu.

**KALDI — AÇIK BORÇLAR (170'ten sonra)**

1. **`c` hangi üyenin sayısıdır?** 170'in en sert bulgusu: `W_X` üyesi
   BANT yönünde yalnız λ = 1.00'de (ν = 0.930), λ yönünde yalnız
   λ ≈ 0.65'te (ν = 1.019) doğrudur (§K0.3). İki ölçütü birden sağlayan
   bir üye ya YOKTUR ya da `KALİB_u2 = A·W_X^{ν(λ,τ)}` biçiminde
   ν'nün kendisi türetilmelidir. **170'in bıraktığı ASIL borç budur** ve
   169'un bıraktığı "λ ekseni" borcunun kesin biçimidir.
2. **`Π_b T(σ_b)` ile ölçülen kırpma aktarımı arasındaki ×1.6.**
   T-yasası bacak başına 0.93 verirken korelatörde ölçülen SERT kırpma
   0.84 veriyor; ikisi aynı yönde hareket ediyor (§K2.4) ama mutlak
   ölçekleri farklı. Kırpmanın "sertliği" ile `σ_b` arasındaki doğru
   köprü (T4) değildir; kapalı formu açık.
3. **β = 0.2175'in türetimi** — 168'den 169'a, 169'dan 170'e devroldu;
   T çerçevesi onu ÜRETEMİYOR (§BONUS-ii, işaret ters, mertebe 1/10).
4. **T/4 dilimlerindeki +%8–11 artık** — T çerçevesinde de erimiyor
   (§BONUS-iii); yeni bulgu: artığın büyük kısmı `g_X` bacağındadır.
5. **Dördüncü bacağın doğrudan kırpılması** (taşıyıcının kendi
   η-izdüşümü `h_Q`) — 169'dan devraldı, hâlâ bir adım uzatma.
6. **`c`'nin λ ≲ 0.6'daki tabanı** iki noktayla (0.70, 0.60) ölçüldü;
   λ = 0.50 ve 0.40 gazları tabanın gerçek olup olmadığını (ve varsa
   değerini) mühürler. Simetrik olarak λ = 1.30/1.50 `c`'nin üst uçta
   sınırsız mı arttığını sınar (169'un `c₀ = 0.4057` kesişimi de bu
   ışıkta yeniden okunmalıdır: o kesişim yalnız λ = 1.00 gazlarından
   kuruluydu).
7. **169'un `c₀ = 4/π² + 0.23σ` kesişimi ne anlama geliyor?** K2′
   4/π²'nin bir limit olmadığını gösterdi; kesim ekseninin φ→0
   kesişiminin λ = 1.00'e kilitli olduğu artık AÇIKTIR ve bu, 169'un
   §K1.3 hükmünün kapsamını daraltır.
