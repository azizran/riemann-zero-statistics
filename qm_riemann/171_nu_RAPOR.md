# 171 — ν ÜYESİ: A(τ) ve M(λ)'NIN KİMLİĞİ (3–4 Eylül 2026, gece-2)

Kalem: `KALEM_NU_UYE_03EYL2026.md` (T1 kaptan kaleminde bitti — `171a`;
bu rapor T2a/T2b/T2c/T3/T4).

Betikler (`171_configs/`):

| betik | ne yapar | koşu | ham çıktı |
|---|---|---|---|
| `171a_nu_carpanlasma.py` ‡ | **T1** (KAPTAN) — çarpanlaşma sınaması | 2 s | (docstring mührü) |
| `171b_insa_uclar.py` | **T2c** — λ=0.50 ve λ=1.30 gazlarının inşası | 9.3 + 9.4 dk | `167/z_L050.npy`, `z_L130.npy` |
| `171c_A_kimlik.py` | **T2a-1** — A(τ) hakemi + aday uyumu + 9-bant ÖN KAYDI | 25 s | `171/A_ONKAYIT.json` |
| `171d_A_yuzlesme.py` | **T2a-2** — 9-bant örneklem-dışı yüzleşme | 30 s | `171/A_YUZLESME.json` |
| `171e_M_kimlik.py` | **T2b** — M(λ) tek-değişkenlik + aday aileler + A-genişliği | 30 s | `171/M_ONKAYIT.json` |
| `171f_onkayit_uclar.py` | **T2c ÖN KAYIT** (korelatöre bakmadan) | 3.3 + 3.4 dk | `171/ONKAYIT_L050.json`, `..._L130.json` |
| `171g_olcum.py` | **T2c ölçüm** (`167_olcum.kos` aynen) | 2 × 3.8 dk | `167/C_L050.json`, `C_L130.json` |
| `171h_yuzlesme_uclar.py` | **T2c yüzleşme** | 30 s | `171/T2C_YUZLESME.json` |
| `171i_T3.py` | **T3** — gerçek gazın dili (ΔM'nin adresi) | 10 s | `171/T3.json` |
| `171j_T4.py` | **T4** — β köprüsü (kesim gazları aynı A üstünde mi) | 10 s | `171/T4.json` |
| `171k_figur.py` | figür + ν merdiveni + c'nin minimumu | 20 s | `171_nu.png`, `171/NU_MERDIVEN.json` |

‡ = depo kökünde (kaptan kalemi); geri kalanı `171_configs/` altında.
Her betiğin docstring'inde ÖN-MÜHÜR (koşudan önce) ve SONUÇ (yalnız gerçek
koşudan) blokları vardır.

(`scratchpad` kökü: `/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/
71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad`; bu görevin dosyaları
`scratchpad/171/` altında: `A_ONKAYIT.json`, `A_YUZLESME.json`,
`M_ONKAYIT.json`, `ONKAYIT_L050.json`, `ONKAYIT_L130.json`,
`T2C_YUZLESME.json`, `T3.json`, `T4.json`, `NU_MERDIVEN.json` + log'lar;
yeni gazlar `scratchpad/167/z_L050.npy`, `z_L130.npy`, `C_L050.json`,
`C_L130.json`.)

Hiçbir ölçüm/çözüm/inşa parçası kopyalanmadı: `167_insa.main`,
`167_olcum.kos`, `167_ortak`, `169_k1.yukle/saglikli/band_jk`,
`166_T1.bant_agg/_jk`, `165_cekirdek.gaz/Model165`,
`163_cekirdek.bant_adaylari`, `166_bacak.karakteristik`
**import edildi**. Git'e dokunulmadı.

---

## 0. TEK CÜMLELİK HÜKÜM

> **A(τ)'nun kimliği bulundu ve SIFIR serbest parametreyle mühürlendi:**
> `A(τ) = W_X(τ; σ_X̃ = σ_X̃(λ=1.00)) = exp(−2π²τ²·0.24204²)` — yani
> λ-değişmez bant şekli, **fiziksel (λ=1) merdivenin X̃ Debye–Waller
> çarpanıdır**, gazın kendi σ_X̃'i değil. Beş bantta artık rms **%0.62**,
> ÖRNEKLEM-DIŞI dokuz-bant uzantısında (τ = 0.7386 ve 0.7774, ön kayıt
> 23:44:17) **−%0.26 ve −%3.69**, dokuz-bant rms **%1.49** — bir serbest
> parametreli en iyi uyumdan (saf Gauss, α = 1.146, rms %1.58) DAHA İYİ.
> Rakiplerin hepsi ön-kayıtlı biçimde öldü: `W_amp·W_X` (DW-tam) +%16/+%18,
> 165'in yarı-analitik F ailesi +%12/+%24, doymuş tarak-sayımı +%9/+%17,
> çıplak tarak-sayımı −%77/−%86 (işaret), λ-değişmez genişlik taban −%11/−%18.
>
> **Buna karşılık M(λ)'nın kimliği BULUNAMADI ve `c`'nin tabanı YOK.**
> `M` hiçbir marjinalin tek-değişkenli fonksiyonu değildir (aynı σ_X̃'te
> dokuz gaz −%17.8 ile +%5.0 arasına yayılıyor, −17.7σ … +8.4σ); `M`
> λ EKSENİNİN ETİKETİDİR ve çarpan adresi `g_E·θ`'dır (yaklaşık eşit
> paylı; `g_X` λ-değişmez). ÖN-KAYITLI iki yeni gaz (λ = 0.50 ve 1.30,
> ön kayıtlar **00:00:05** ve **00:03:56**) şunu verdi:
> `c(0.50) = 0.3800 ± 0.0072`, `M(0.50) = 1.2266 ± 0.0069`;
> `c(1.30) = 0.4540 ± 0.0136`, `M(1.30) = 0.9832 ± 0.0157`.
> **λ = 1.30 tek-taraflılığı doğruladı** (M, 1'den −1.1σ; Gram-doyum
> adayı +2.8σ ile öldü), ama **λ = 0.50 bütün türetimli M ailelerini
> öldürdü** (+9.4σ … +15.3σ); ayakta kalan tek eğri EMPİRİK log-λ
> kuadratiğidir (+0.4σ, kimlik değil). Ve en sert sürpriz:
> **`c(λ)` bir tabana oturmuyor — λ* = 0.6487'de bir MİNİMUMU var**
> (ν(0.60→0.50) = 1.489 ± 0.063, ν = 1'i **+7.7σ** aşıyor), yani 170
> §K2'nin "λ ≤ 0.70'te taban" okuması örneklem-dışı ÖLDÜ.
> Gerçek ζ gazının fazlası ise **saf seviyedir**: şekli A'nın tam
> ortasında (+%1.19), ΔM = +%5.61'in **%93'ünü** ölçülen 3-bacak kırpma
> aktarımı kapatıyor. β = 0.2175 ise **seviyenin (θ ⊂ M) yasasıdır,
> şeklin (A) değil** — kesim gazlarından yalnız K090 aynı A üstündedir.

---

## T2a — A(τ)'NUN KİMLİĞİ

### T2a.0 Hakem: A(τ) gerçek bant KALİB'inden (`171c_A_kimlik.py`)

171a, A(τ)'yu `KALİB_b = c_bant_b · W_X_b(parametrik)` yeniden kurulumuyla
mühürlemişti. Burada A doğrudan **ölçülen bant KALİB_u2'sünden**
(`169_k1.band_jk`) kuruldu; fark ‰1–2.5:

| gaz | λ | S₁ | S₂ | S₃ | S₄ | S₅ |
|---|---|---|---|---|---|---|
| L115 | 1.15 | 1.1046 | 1.0465 | 1.0000 | 0.9483 | 0.9037 |
| Hkeskin | 1.00 | 1.1070 | 1.0464 | 1.0000 | 0.9433 | 0.8904 |
| L085 | 0.85 | 1.1088 | 1.0462 | 1.0000 | 0.9387 | 0.8767 |
| L070 | 0.70 | 1.1074 | 1.0453 | 1.0000 | 0.9358 | 0.8669 |
| L060 | 0.60 | 1.1014 | 1.0426 | 1.0000 | 0.9364 | 0.8679 |
| **A(τ)** (5-gaz geo. ort.) | — | **1.1058** | **1.0454** | **1.0000** | **0.9405** | **0.8810** |
| τ_ort | | 0.5393 | 0.5792 | 0.6189 | 0.6580 | 0.6982 |
| gazlar arası s.h. | | 0.12% | 0.07% | — | 0.25% | 0.80% |
| bant-içi jackknife | | 0.60% | 0.68% | 0.41% | 0.77% | 0.81% |

171a'nın A'sı (1.1048/1.0448/1.0000/0.9423/0.8832) ile fark
+0.09/+0.06/0.00/−0.19/−0.25 % — **çarpanlaşma yeniden-kurulumdan
bağımsızdır.**

### T2a.1 Adaylar ve ön-kayıtlı ölümler (5 bant, `171c`)

Bütün adaylar **YALNIZ beş bantta** (τ ∈ [0.539, 0.698]) uyduruldu;
dokuz-bant öngörüleri ölçümden önce `A_ONKAYIT.json`'a yazıldı
(zaman damgası **2026-09-03 23:44:17**).

| aday | türetim | serbest p | param | 5-bant rms | 5-bant en kötü |
|---|---|---|---|---|---|
| A0 düz | NUL çıpa | 0 | — | 7.87% | 11.90% |
| **A1 DW-tam (Hkeskin)** `W_amp·W_X` | c_ampX üyesi bant-düz olsaydı | 0 | — | **6.48%** | 9.43% |
| A1λ DW-tam (λ-ailesi geo.) | aynı, λ-ortalama σ | 0 | — | 3.83% | 5.38% |
| **A2Hk `W_X`-tek, σ_X̃(λ=1.00)** | ν_bant = 1, λ=1 gazında | **0** | — | **0.62%** | **1.05%** |
| A2sn `W_X`-tek, σ_X̃(son) | gerçek gazın genişliği | 0 | — | 0.77% | 1.38% |
| A2iv `W_X`-tek, λ-değişmez taban | σ_X̃² = A_Xλ²+B_X ⇒ √B_X = 0.1475 | 0 | — | 4.95% | 7.86% |
| A3f 165 §6b yarı-analitik F | 0.546sin²(2πτ)+0.666sin⁴(πτ) | 0 | — | 13.68% | 22.63% |
| A4 tarak-sayım (çıplak) | 168 §A1.3: `n_Δ·w_eff ∝ e^{τL}/(τL)` | 0 | — | **79.15%** | 150.93% |
| A2s saf Gauss `e^{−ατ²}` | **≡ DW-tam, σ* serbest** | 1 | α = 1.146 | 0.61% | 1.00% |
| A3s 165-F ailesi (c₂/c₁) | yarım-gap ailesinin iki üyesi | 1 | 2.457 | 2.52% | 4.49% |
| A4s doymuş tarak-sayım | `1/(1+κe^{τL})` | 1 | 6.45e−5 | 1.97% | 3.81% |
| A5 kuvvet `τ^{−p}` | **EMPİRİK ÇIPA** (türetimsiz) | 1 | 0.861 | 1.45% | 2.26% |

**Yapısal not (yeni).** Kalemin ilk iki aday ailesi — "`W_amp·W_X`
sabit-σ*" ile "saf Gauss `e^{−ατ²}`" — **AYNI TEK-PARAMETRELİ AİLEDİR**:

```
W_amp(τ)·W_X(τ) = exp[−½(πτ)²σ_ds²]·exp[−½(2πτ)²σ_X̃²]
                 = exp[−½π²τ²(σ_ds² + 4σ_X̃²)] ≡ e^{−ατ²},
                 α = ½π²σ*²,  σ*² = σ_ds² + 4σ_X̃².
```

Yarışan şey ailenin kendisi değil, **σ*'ın KİMLİĞİDİR**.

### T2a.2 σ*'ın kimliği — `A(τ)` λ=1 gazının `W_X`'idir

Serbest uyum α = **1.146** ⇒ σ*_eff = √(2α)/π = **0.48190** ⇒ `W_X`
dilinde **σ_X̃-eşdeğeri = 0.24095**. Bunun kimlik defteri:

| aday kimlik | değer | σ_X̃-eşdeğeri / bu |
|---|---|---|
| **σ_X̃(Hkeskin, λ = 1.00)** | **0.24204** | **−0.46%** |
| σ_X̃(son, GERÇEK ζ) | 0.23384 | +3.03% |
| σ_X̃(λ-ailesi geo. ort.) | 0.21845 | +10.29% |
| σ_X̃(λ-değişmez taban √B_X) | 0.14749 | +63.35% |
| σ_X̃(L115) | 0.25839 | −6.76% |
| σ_X̃(L060) | 0.17951 | +34.21% |
| σ_Ĉ(Hkeskin) | 0.27768 | −13.24% |
| σ_ds(Hkeskin)/2 | 0.21567 | +11.71% |

σ_X̃-eşdeğerini veren λ (log-log kuadratik ara değer): **λ_eff = 0.9886**.

> **HÜKÜM (T2a-a).** `A(τ)`'nun genişliği **λ = 1.00 gazının σ_X̃'idir**
> (%0.5 içinde; λ_eff = 0.989). Gazın KENDİ σ_X̃'i değildir (L115'te −%6.8,
> L060'ta +%34 uzak) — λ-değişmezliğin bütün içeriği budur. `W_amp`
> çarpanı **yoktur**: onu eklemek σ*'ı 0.482'den 0.648'e çıkarır ve
> beş bantta ±%9.4 artık bırakır. **168 §A2.8'in "`W_amp` çift sayımdır"
> hükmü burada bant şekli düzeyinde bağımsız olarak doğrulanmıştır.**

### T2a.3 ÖRNEKLEM-DIŞI SINAV: dokuz bant (`171d_A_yuzlesme.py`)

Dokuz-bant penceresi (lo ≤ 0.80) iki yeni bant ekler: τ ≈ **0.7386**
(beş λ-gazının hepsinde) ve τ ≈ **0.7774** (Hkeskin + L085). Uyum
penceresinin **dışındadır**; her adayın oradaki değeri 23:44:17'de
yazılmıştı.

**Ölçülen:** A(0.7386) = **0.8265 ± 0.0122**,
A(0.7774) = **0.7456 ± 0.0223** (gazlar arası s.h.).

| aday | ön-kayıt A(0.7386) | ön-kayıt A(0.7774) | sapma₆ | sapma₇ | **9-bant rms** |
|---|---|---|---|---|---|
| **A2Hk (0 param)** | **0.8287** | **0.7742** | **−0.26%** | **−3.69%** | **1.49%** |
| A2s (1 param, α=1.146) | 0.8301 | 0.7760 | −0.43% | −3.92% | 1.58% |
| A2sn (0 param) | 0.8391 | 0.7875 | −1.50% | −5.32% | 2.19% |
| A5 kuvvet (EMPİRİK) | 0.8588 | 0.8218 | −3.76% | −9.27% | 3.97% |
| A1λ DW-tam λ-geo | 0.7595 | 0.6876 | +8.82% | +8.44% | 5.64% |
| A4s doymuş tarak-sayım | 0.7575 | 0.6370 | +9.12% | +17.04% | 7.49% |
| A2iv λ-değişmez taban | 0.9326 | 0.9093 | −11.37% | −18.01% | 9.07% |
| A3s 165-F ailesi | 0.7367 | 0.6034 | +12.19% | +23.57% | 10.25% |
| A1 DW-tam Hkeskin | 0.7138 | 0.6318 | +15.79% | +18.01% | 10.58% |
| A0 düz | 1.0000 | 1.0000 | −17.35% | −25.44% | 13.40% |
| A3f 165-F sabit katsayı | 0.9759 | 0.8562 | −15.31% | −12.92% | 13.82% |
| A4 tarak-sayım çıplak | 3.5373 | 5.3597 | −76.63% | −86.09% | 79.83% |

> **HÜKÜM (T2a-b — MÜHÜR).** `A(τ) = exp(−2π²τ²σ_X̃(λ=1)²)` **sıfır
> serbest parametreyle** dokuz bantta rms %1.49 verir ve **bir serbest
> parametreli en iyi uyumu (A2s, %1.58) YENER** — çünkü serbest α, beş
> bantta ‰1 kazanıp uzantıda ‰2 kaybediyor. Bütün türetimli rakipler
> örneklem-dışı noktada **%9–86 arası** ıskalayarak ÖLDÜ.
>
> **Tarak-sayımının ölümü İŞARETTEDİR** (kalemde böyle ön-kayıtlıydı):
> merdiven yoğunluğu `dn/dω = e^ω/ω` τ = 0.539 → 0.698 boyunca **×6.0
> ARTIŞ** ister, ölçülen ×0.797 **DÜŞÜŞ**tür. Doyumlu biçimi (1 param)
> beş bantta %2'ye iner ama uzantıda +%17 ile ölür. **Rezonans doluluk
> oranı `n_Δ·w_eff` bant şeklinin adresi DEĞİLDİR.**
>
> **165'in yarı-analitik F'i de ölür.** Sabit katsayılı biçim (0.546/0.666)
> düşük τ'da +%22.6 ıskalar; serbest c₂/c₁ = 2.46 beş bantta %2.5'e iner
> ama uzantıda +%12/+%24 ile ölür. `F(τ_q)` çizgi-başına yoğunluktur,
> **bant şekli değildir** — 165 §6a'nın "Ç1'in F'i τ_band'dan bağımsızdır"
> hükmüyle tutarlı.

### T2a.4 DÜRÜSTLÜK: çarpanlaşma τ > 0.70'te ZAYIFLIYOR

6. bantta (τ = 0.7386) gazların S değerleri:

| gaz | L115 | Hkeskin | L085 | L070 | L060 |
|---|---|---|---|---|---|
| S₆ | 0.8634 | 0.8449 | 0.8235 | 0.8039 | 0.7988 |

**Yayılım %8.08** ve λ ile **tekdüze sıralı** — beş bantlık pencerede
kolaps rms %0.81 iken. Yani `KALİB = A(τ)·M(λ)` çarpanlaşması
τ ≲ 0.70'te ‰8, τ ≈ 0.74'te ancak %8'dir. **A(τ) tam λ-değişmez
değildir**; gaz-başına Gauss genişliği (§T2b.4) λ ile ±%5 kayar ve bu
kayma yüksek τ'da büyütülür. A(0.7386)'nın ±%1.5'lik s.h.'si bu
yayılımdan gelir; hüküm tablosundaki −%0.26 bu belirsizliğin içindedir.
**Bu, hükmü değil GEÇERLİLİK PENCERESİNİ sınırlar: A·M ayrışması
τ ≤ 0.70'in yasasıdır.**

---

## T2b — M(λ)'NIN KİMLİĞİ (`171e_M_kimlik.py`)

M burada A'yı **özdeş olarak eleyen** biçimde tanımlandı:
`M(g) = geo.ort_b [KALİB_b(g) / KALİB_b(Hkeskin)]` (beş bant). Bu, 171a'nın
orta-bant M'sinden ‰3–8 farklıdır (bant-ortalaması gürültüyü √5 azaltır).

| gaz | λ | bant bant oran | **M** | ± |
|---|---|---|---|---|
| L115 | 1.15 | 0.9964 0.9986 0.9985 1.0039 1.0135 | **1.0022** | 0.0054 |
| Hkeskin | 1.00 | 1 1 1 1 1 | 1.0000 | 0.0044 |
| L085 | 0.85 | 1.0199 1.0181 1.0182 1.0133 1.0026 | **1.0144** | 0.0054 |
| L070 | 0.70 | 1.0674 1.0658 1.0670 1.0586 1.0389 | **1.0595** | 0.0071 |
| L060 | 0.60 | 1.1242 1.1257 1.1299 1.1217 1.1014 | **1.1205** | 0.0069 |
| **son** | (1.00) | 1.0662 1.0549 1.0592 1.0564 1.0438 | **1.0561** | 0.0060 |

*(ÖN-MÜHÜR N1 orta-bant M'sini 0.999/1.000/1.018/1.067/1.129 demişti;
bant-ortalamalı M sistematik olarak %0.3–0.8 aşağıda çıktı — 1σ içinde,
fark tanımdandır. Kurtarma yok: N1 tam isabet değil.)*

### T2b.1 M'nin çarpan ayrışması — ÖZDEŞLİK (‰0'da kapanıyor)

`KALİB_u2 = g_E·g_X²·θ` olduğundan `M = (g_E/g₀)(g_X/g₀)²(θ/θ₀)`
**özdeştir** ve dört hanede kapanır (fark ≤ 0.0005%). İçerik, hangi
çarpanın taşıdığındadır:

| gaz | λ | g_E/g₀ | (g_X/g₀)² | θ/θ₀ | M |
|---|---|---|---|---|---|
| L115 | 1.15 | 0.9985 | 1.0068 | 0.9968 | 1.0022 |
| L085 | 0.85 | 1.0101 | 0.9960 | 1.0084 | 1.0144 |
| L070 | 0.70 | 1.0342 | 0.9974 | 1.0271 | 1.0595 |
| **L060** | 0.60 | **1.0631** | 1.0045 | **1.0493** | 1.1205 |
| **son** | — | **1.0249** | 1.0015 | **1.0290** | **1.0561** |
| K090 | (kesim) | 0.9957 | 0.9981 | **0.9080** | 0.9023 |
| K070 | (kesim) | 1.0475 | 1.0317 | **0.8038** | 0.8687 |
| HA4 | (kesim) | 1.0649 | 1.0360 | **0.7975** | 0.8799 |
| E060 | (kesim) | 1.1191 | 1.0597 | **0.7945** | 0.9423 |
| HkT2a | (T/2) | 0.9904 | 0.9421 | 1.0317 | 0.9626 |
| HkT4b | (T/4) | 0.9514 | **0.7833** | 1.1066 | 0.8247 |

> **HÜKÜM (T2b-a).** **λ ekseninde M'yi İKİ çarpan neredeyse EŞİT
> paylaşır:** Gram şişmesi `g_E` (λ=0.60'ta +%6.3) ve faz tutarlılığı `θ`
> (+%4.9); `g_X` λ-değişmezdir (≤‰5). Üç eksen ÜÇ AYRI çarpanı
> hareket ettirir: **λ → g_E+θ**, **kesim → yalnız θ** (g_E ters yönde
> ARTAR), **pencere → yalnız g_X²** (168 §BONUS-iii'ün adresi).
> M(λ) "bir doyum fonksiyonu" değil, `g_E·θ` çarpımının λ-profilidir.

### T2b.2 TEK-DEĞİŞKENLİK SINAVI — **M hiçbir marjinalin tek-değişkenli
fonksiyonu DEĞİL**

λ-ailesi içinde λ, σ_ds, σ_X̃, σ_Ĉ **birebir dejeneredir** (hepsi
birbirinin tekdüze fonksiyonu) — "hangi değişken" sorusu bu aileden
cevaplanamaz. Dejenerasyonu kıran gazlar λ = 1.00 kurulumundadır ama
σ'ları farklıdır. `M(σ_X̃)` λ-ailesinden log-log kuadratikle kurulup
onlara uygulanınca:

| gaz | σ_X̃ | λ_eş(σ_X̃) | M öngörü | M ölçülen | sapma |
|---|---|---|---|---|---|
| **son** | 0.2338 | 0.9335 | 1.0056 | 1.0561 | **+5.02% (+8.4σ)** |
| HkT2a | 0.2412 | 0.9934 | 1.0020 | 0.9626 | −3.93% (−6.0σ) |
| HkT2b | 0.2428 | 1.0073 | 1.0015 | 0.9148 | −8.65% (−16.8σ) |
| HkT4a | 0.2412 | 0.9931 | 1.0020 | 0.8238 | −17.78% (−16.3σ) |
| HkT4b | 0.2412 | 0.9937 | 1.0020 | 0.8247 | −17.69% (−12.9σ) |
| K090 | 0.2446 | 1.0228 | 1.0010 | 0.9023 | −9.86% (−17.7σ) |
| K070 | 0.2421 | 1.0014 | 1.0017 | 0.8687 | −13.27% (−14.8σ) |
| HA4 | 0.2396 | 0.9800 | 1.0026 | 0.8799 | −12.24% (−6.5σ) |
| E060 | 0.2369 | 0.9576 | 1.0039 | 0.9423 | −6.13% (−1.5σ) |

(σ_ds ve σ_Ĉ için λ_eş sütunları da hesaplandı; hüküm değişmiyor.)

> **HÜKÜM (T2b-b — kalemin ilk sorusunun cevabı).** **M yalnız σ_X̃
> üzerinden DEĞİLDİR; σ_ds ya da σ_Ĉ üzerinden de değildir.** Aynı
> σ_X̃'te dokuz gaz −%17.8 ile +%5.0 arasına yayılıyor (en küçük
> ayrışma −1.5σ, en büyüğü −17.7σ). `M`, bir MARJİNALİN fonksiyonu
> değil, **λ EKSENİNİN ETİKETİDİR**; başka eksenlere taşınamaz.
> Çarpanlaşma `KALİB = A(τ)·M` λ ekseninin yasasıdır — evrensel bir
> ayrışma değil. (Bunun T3 için doğrudan sonucu: gerçek gazın ΔM'si
> σ_X̃'ten öngörülemez; §T3.)

### T2b.3 Aday M aileleri — 1. sınav TEK-TARAFLILIK

Eşik ailelerinin eşiği **serbest bırakılmadı**: inşa denkleminin kendi
doyum noktası `rms S'(λ) = N̄'` ⇒ **λ_c = 1.9147/1.894 = 1.0109**
(170 §K0.2). Bu, sıfır ek parametredir.

| aday | türetim | p | param | rms | λ=1.15 öngörüsü (ölç 1.0022) | hüküm |
|---|---|---|---|---|---|---|
| M0 düz (168 §A3) | KALİB λ-değişmez | 0 | — | 6.05% | 1.0000 | λ<1'de ölü |
| M8 ölçülen `(g_E/g₀)(g_X/g₀)²` | θ λ-değişmez | 0 | — | 2.55% | 1.0053 | tek-taraflı ✓, %5 eksik |
| M1 Gram-doyum `1/g_E = 1+Cλ` | κ-dışı ∝ 𝒢 ∝ λ | 1 | C=0.273 | 2.09% | **0.9689 (−3.3%)** | **ÖLDÜ** |
| M1b Gram-doyum `1+Cλ²` | ikinci mertebe | 1 | C=0.140 | 2.55% | **0.9620 (−4.0%)** | **ÖLDÜ** |
| **M2 kırık kuvvet** `(λ_c/λ)^p`, λ_c sabit | inşa doyum eşiği | 1 | p=0.190 | **1.19%** | 1.0000 | **AYAKTA** |
| **M3 varyans-açığı** `1+κ(u_c−u)_+/u_c` | aynı eşik, u=σ_X̃² | 1 | κ=0.232 | **1.24%** | 1.0000 | **AYAKTA** |
| M7 merdiven-sürüşlü kesir `f_L^{−q}` | `f_L = A_Xλ²/(A_Xλ²+B_X)` | 1 | q=0.197 | 1.29% | **0.9811 (−2.1%)** | **ÖLDÜ** |
| M9 log-λ kuadratik | **EMPİRİK ÇIPA** | 2 | — | 0.04% | 1.0025 | (siren) |
| *[ref] ΠT/ΠT₀ (170 §K3.3)* | *170'te öldü, aday değil* | — | — | 3.41% | *0.9576 (−4.5%)* | *ölü* |

> **HÜKÜM (T2b-c).** **TEK-TARAFLILIK üç türetimli adayı öldürdü**
> (M1, M1b, M7 — hepsi λ ≥ 1'de düşme istiyor, ölçüm düz). Ayakta kalan
> iki aday da **AYNI eşiği** kullanıyor: inşa denkleminin doyum noktası
> `λ_c = 1.0109`. Bu eşik uydurulmadı, 170 §K0.2'den geldi ve
> M(1.15) = 1.0022 ± 0.0054 (1'den **+0.4σ**) onu doğruluyor.
> **λ = 1.30 gazı bu hükmün örneklem-dışı hakemidir** (§T2c): M2/M3
> tam 1.000 der, M1 0.940, M7 0.967 derdi.

### T2b.4 A(τ)'nun gaz-başına genişliği — λ-değişmezlik ne kadar tam?

Her gazın kendi bant KALİB'ine `e^{−ατ²}` uydurulup `σ*_g/2` okundu:

| gaz | λ | α_g | σ*_g/2 | /σ_X̃(Hk) | /σ_X̃(kendi) |
|---|---|---|---|---|---|
| L115 | 1.15 | 1.0158 | 0.22685 | −6.27% | −12.21% |
| Hkeskin | 1.00 | 1.0977 | 0.23582 | −2.57% | −2.57% |
| L085 | 0.85 | 1.1789 | 0.24438 | +0.97% | +9.73% |
| L070 | 0.70 | 1.2264 | 0.24926 | +2.98% | +25.27% |
| L060 | 0.60 | 1.1941 | 0.24596 | +1.62% | +37.02% |
| son | — | 1.1841 | 0.24493 | **+1.19%** | +4.74% |
| K090 | kesim | 1.2195 | 0.24855 | +2.69% | +1.60% |
| K070 | kesim | 1.4213 | 0.26834 | **+10.86%** | +10.81% |
| HA4 | kesim | 1.7038 | 0.29380 | **+21.38%** | +22.62% |
| E060 | kesim | 2.3001 | 0.34136 | **+41.03%** | +44.12% |
| HkT2a | T/2 | 0.9892 | 0.22386 | −7.51% | −7.19% |
| HkT4b | T/4 | 0.7087 | 0.18948 | **−21.71%** | −21.45% |

λ-ailesi ortalaması **0.24045** = σ_X̃(Hkeskin) − %0.66; yayılım **%9.9**.

> **HÜKÜM (T2b-d — T2a'nın geçerlilik sınırı).** A(τ)'nun genişliği
> λ ekseninde ±%5 içinde sabittir (ortalaması σ_X̃(λ=1)) ama **tam sabit
> değildir** ve kayma tekdüzedir (λ büyüdükçe daha SIĞ şekil). §T2a.4'ün
> 6. bant %8'lik yayılımının kaynağı budur. **Kesim ve pencere
> eksenlerinde A hiç sabit değildir** (E060 +%41, HkT4b −%22) — bu, T4'ün
> ilk cevabıdır: **kesim gazları AYNI A(τ) üstünde DEĞİLDİR.**

---

## T2c — YENİ GAZLAR: λ = 0.50 ve λ = 1.30 (ÖN-KAYITLI)

### T2c.0 İnşa (`171b_insa_uclar.py`) ve inşa ön-mührünün yüzleşmesi

Her iki gaz da `167_insa.main` ile, Hkeskin/L060/L115 ile **AYNI**
merdiven ve çözücüyle kuruldu (`S(t) = −Σ λ a_q w_q sin ω_q t`, keskin
τ ≤ 1.00, c = −½, 300000 sıfır, h = 0.015).

| büyüklük | ön-mühür (171b docstring) | L050 ölçülen | hüküm | ön-mühür | L130 ölçülen | hüküm |
|---|---|---|---|---|---|---|
| rms S′ | 0.947 / 2.462 | **0.9470** | ✓ | | **2.4622** | ✓ |
| σ_ds | [0.279, 0.285] / [0.477, 0.487] | **0.27903** | ✓ (alt sınırda) | | **0.48890** | **✗ (+%0.4 ÜSTÜNDE)** |
| σ_X̃ | [0.1555, 0.1590] / [0.2685, 0.2742] | **0.15576** | ✓ | | **0.27713** | **✗ (+%1.1)** |
| σ_Ĉ | [0.1795, 0.1835] / [0.3073, 0.3138] | **0.18131** | ✓ | | **0.31923** | **✗ (+%1.7)** |
| ΔG<0 kesri | 0.019 ± 0.005 / 0.26 ± 0.03 | **0.0106** | **✗ (düşük)** | | **0.2511** | ✓ |
| min Δz | ~0.21 / ~0.055 (0.03–0.08) | **0.18913** | ✗ (−%10) | | **0.08348** | ✗ (üstünde) |
| sıralılık | TAM (λ=1.30 RİSKLİ) | **TAM** | ✓ | | **TAM** | ✓ (risk gerçekleşmedi) |

maks|F| = 1.86·10⁻⁹ ikisinde de; ilk-kök hücreleri benzersiz 300000/300000
ikisinde de; L130'da 6 ikiye-bölme adımı gerekti (L050'de 0).

> **Kurtarmasız not.** λ = 1.30'un σ'ları ön-mührün ÜSTÜNDE çıktı (+%0.4
> … +%1.7). Kuvvet-yasası uzatması yüksek λ'da sistematik olarak AZ
> veriyor — 170 §K2′ aynı ıskayı L115'te (%1.1) yapmıştı. **Bu üçüncü
> tekrardır ve artık bir yasa gibi okunmalıdır: σ(λ) uzatması λ > 1'de
> daima aşağıdan ıskalar.**

### T2c.1 ÖN KAYIT (`171f_onkayit_uclar.py`) — korelatöre bakılmadan

Zaman damgaları: **L050 = 2026-09-04 00:00:05**, **L130 = 00:03:56**
(ölçümler 00:05 ve 00:08'de başladı). Öngörü zinciri:
`KALİB_b(yeni) = KALİB_b(Hkeskin)·M` (A özdeş elenir),
`c_b = KALİB_b/W_X_b` — `W_X_b` yeni gazın **marjinalinden** ölçüldü.
M adaylarının parametreleri 171e'de donduruldu.

| aday | M(0.50) | c(0.50) | M(1.30) | c(1.30) |
|---|---|---|---|---|
| M0 düz | 1.0000 | 0.3098 | 1.0000 | 0.4618 |
| **M2 kırık kuvvet** (λ_c = 1.0109) | 1.1430 | 0.3541 | **1.0000** | 0.4618 |
| **M3 varyans-açığı** | 1.1347 | 0.3515 | **1.0000** | 0.4618 |
| M8 ölçülen `(g_E/g₀)(g_X/g₀)²` | 1.1309 | 0.3504 | 1.0078 | 0.4654 |
| M1 Gram-doyum *(171e'de ölü)* | 1.1200 | 0.3470 | 0.9396 | 0.4339 |
| M7 merdiven kesri *(171e'de ölü)* | 1.1614 | 0.3598 | 0.9674 | 0.4467 |
| M9 log-λ kuadratik **[EMPİRİK]** | 1.2238 | 0.3791 | 1.0175 | 0.4699 |

### T2c.2 ÖLÇÜM (`171g_olcum.py`, 3.8 dk her biri) ve YÜZLEŞME (`171h`)

```
c(L050) = 0.3800 ± 0.0010(jk) ± 0.0072(bant) = ±0.0072   [5 bant]
M(L050) = 1.2266 ± 0.0069     (bant bant 1.2105 1.2225 1.2368 1.2363 1.2273)

c(L130) = 0.4540 ± 0.0019(jk) ± 0.0135(bant) = ±0.0136   [5 bant]
M(L130) = 0.9832 ± 0.0157     (bant bant 0.9594 0.9515 0.9743 0.9987 1.0343)
```

| ön-kayıtlı aday | **λ = 0.50** M sapması | **λ = 1.30** M sapması | hüküm |
|---|---|---|---|
| M0 düz | +22.66% (**+32.6σ**) | −1.68% (−1.1σ) | λ<1'de ölü (bilinen) |
| **M2 kırık kuvvet** | **+7.32% (+12.0σ)** | −1.68% (**−1.1σ**) | **DÜŞÜK UÇTA ÖLDÜ** |
| **M3 varyans-açığı** | **+8.11% (+13.2σ)** | −1.68% (**−1.1σ**) | **DÜŞÜK UÇTA ÖLDÜ** |
| M8 ölçülen çarpan | +8.46% (+13.8σ) | −2.44% (−1.6σ) | düşük uçta öldü |
| M1 Gram-doyum | +9.52% (+15.3σ) | +4.64% (+2.8σ) | iki uçta da ölü |
| M7 merdiven kesri | +5.62% (+9.4σ) | +1.64% (+1.0σ) | düşük uçta öldü |
| **M9 log-λ kuadratik [EMPİRİK]** | **+0.23% (+0.4σ)** | −3.37% (−2.2σ) | **tuttu (ama çıpa)** |
| (çapa) 4/π² — c dilinde | −6.24% (−3.5σ) | **+12.02% (+3.6σ)** | — |

> **HÜKÜM (T2c-a — TEK-TARAFLILIK AYAKTA).** `M(1.30) = 0.9832 ± 0.0157`
> 1'den **−1.1σ**: `M`, λ ≥ 1'de gerçekten düzdür ve bu artık İKİ
> örneklem-dışı noktayla (λ = 1.15 ve 1.30) desteklenmektedir.
> λ = 1.30 Gram-doyum adayını (**M1, öngörü 0.9396**) +2.8σ ile bir kez
> daha öldürdü.
>
> **HÜKÜM (T2c-b — TÜRETİMLİ M AİLELERİNİN TOPLU ÖLÜMÜ).** Düşük uçta
> **her türetimli aday öldü**: M2 +12.0σ, M3 +13.2σ, M8 +13.8σ, M7 +9.4σ.
> Ölçülen `M(0.50) = 1.2266`, inşa-doyum eşiğinin verdiği 1.1430'un
> **%7.3 üstünde**. Yani `M`'nin yükselen kolu, eşik ailelerinin
> öngördüğünden **daha hızlı** dikleşiyor. Ayakta kalan tek eğri,
> beş λ noktasından kurulmuş **EMPİRİK** log-λ kuadratiğidir
> (öngörü 1.2238, ölçüm 1.2266, **+0.4σ**) — bu bir KİMLİK DEĞİLDİR ve
> siren kuralı gereği öyle sayılmaz. **T2b'nin kimlik sorusu AÇIK
> KALMIŞTIR.**
>
> **HÜKÜM (T2c-c — A'NIN GEÇERLİLİK PENCERESİ λ ≤ 1.15).** λ = 0.50
> şekilde A(τ)'nun üyesidir (S/A2Hk rms **%1.50**, σ*/2 = **0.22779**,
> λ-ailesi bandı [0.22685, 0.24926] içinde). λ = 1.30 **DEĞİLDİR**:
> rms **%3.62**, σ*/2 = **0.18661** — bandın **%18 ALTINDA**, bant-bant
> oran 0.959'dan 1.034'e tekdüze tırmanıyor. **Çarpanlaşma λ ∈ [0.50, 1.15]
> penceresinin yasasıdır; λ = 1.30'da şekil kanadı kırılır.** (M(1.30)'un
> ±0.0157'lik hatasının tamamı bu şekil sürüklenmesidir.)

### T2c.3 EN BÜYÜK SÜRPRİZ: `c(λ)`'NIN TABANI YOK (`171k_figur.py`)

Yedi noktalı defter:

| gaz | λ | **c ± σ_tot** | M | σ_X̃ |
|---|---|---|---|---|
| L130 | 1.30 | **0.4540 ± 0.0136** | 0.9832 | 0.27713 |
| L115 | 1.15 | 0.4312 ± 0.0053 | 1.0022 | 0.25839 |
| Hkeskin | 1.00 | 0.4035 ± 0.0018 | 1.0000 | 0.24204 |
| L085 | 0.85 | 0.3817 ± 0.0028 | 1.0144 | 0.22272 |
| L070 | 0.70 | 0.3690 ± 0.0058 | 1.0595 | 0.19897 |
| L060 | 0.60 | 0.3689 ± 0.0072 | 1.1205 | 0.17951 |
| **L050** | **0.50** | **0.3800 ± 0.0072** | **1.2266** | 0.15576 |
| son | (1.00) | 0.4122 ± 0.0023 | 1.0561 | 0.23384 |

ν merdiveni (ν_λ = Δlog KALİB/Δlog W_X, 5-bant ortalaması, jackknife):

| λ adımı | 0.60→0.50 | 0.70→0.60 | 0.85→0.70 | 1.00→0.85 | 1.15→1.00 | 1.30→1.15 |
|---|---|---|---|---|---|---|
| **ν** | **+1.489 ± 0.063** | +0.993 ± 0.073 | +0.563 ± 0.055 | +0.205 ± 0.062 | −0.034 ± 0.069 | **+0.271 ± 0.073** |
| ν − 1 | **+7.7σ** | −0.1σ | −8.0σ | −12.8σ | −15.0σ | −10.0σ |

> **HÜKÜM (T2c-d — 170 §K2'NİN 'TABAN' OKUMASI ÖLDÜ).** `c(λ)` λ ≤ 0.70'te
> bir tabana oturmuyor; **λ* = 0.649'da bir MİNİMUMU var** (c_min = 0.3675)
> ve λ = 0.50'de yeniden **0.3800**'e çıkıyor. Doğrudan c farkı yalnız
> +1.1σ'dır, ama ayrımı yapan büyüklük ν'dür: **ν(0.60→0.50) = 1.489 ±
> 0.063, yani 1'i +7.7σ ile aşıyor** (ν = 1 ⟺ dc/dλ = 0, özdeşlik).
> Aynı özdeşlik minimumun yerini de veriyor: ν'nün 1'i kestiği λ ≈ 0.65,
> parabolün verdiği λ* = 0.6487 ile aynıdır.
>
> **Ve ν(λ) TEKDÜZE DEĞİL:** λ = 1.15 → 1.30'da −0.03'ten **+0.27**'ye
> geri tırmanıyor (+3.0σ). 170 §K0.3'ün "ν λ = 1.15 → 0.60 boyunca
> tekdüze tırmanır" resmi, **ν'nün λ ≈ 1.08 civarındaki MİNİMUMUNUN
> sağ tarafını görmemiş olmasıdır.** `c(λ)` eğrisi bir U'dur, ν(λ) ise
> ters bir U'nun aynası: `c` λ* = 0.649'da minimum, ν orada 1'den geçiyor.
>
> **4/π² üstüne:** `c(1.30) = 0.4540 ± 0.0136` = 4/π² + **%12.0 (+3.6σ)**.
> `c` λ ile artmaya devam ediyor; 170 §K2′'nin "4/π² bir limit değil,
> λ = 1.00'in tesadüfüdür" hükmü ikinci kez örneklem-dışı doğrulandı.

---

## T3 — GERÇEK GAZIN DİLİ: ΔM = +%5.6'NIN ADRESİ (`171i_T3.py`)

```
M(son) = 1.0561 ± 0.0060   ⇒   ΔM = +5.61%  =  +9.4σ
ayrışma:  g_E/g₀ = 1.0249 (+2.49%)   (g_X/g₀)² = 1.0015 (+0.15%)
          θ/θ₀   = 1.0290 (+2.90%)
log-paylar:  g_E %45   g_X² %3   θ %52
```

### T3.1 Gerçek gaz A(τ)'nun ÜYESİDİR (şekil), M'de DEĞİLDİR (seviye)

`S/A2Hk − 1` (%) — dokuz bant, ilk beşi uyum penceresi:

| gaz | b1 | b2 | b3 | b4 | b5 | **b6 (ÖD)** | **b7 (ÖD)** | 5-bant rms |
|---|---|---|---|---|---|---|---|---|
| L115 | −0.70 | −0.94 | 0 | +0.51 | +2.00 | +4.20 | — | 1.06% |
| Hkeskin | −0.49 | −0.95 | 0 | −0.05 | +0.48 | +1.96 | −0.78 | 0.53% |
| L085 | −0.33 | −0.97 | 0 | −0.55 | −1.07 | −0.63 | −6.53 | 0.71% |
| L070 | −0.45 | −1.06 | 0 | −0.87 | −2.19 | −3.01 | — | 1.17% |
| L060 | −0.99 | −1.31 | 0 | −0.81 | −2.08 | −3.62 | — | 1.24% |
| **son** | +0.18 | −1.35 | 0 | −0.34 | −1.01 | **−2.85** | **−7.12** | **0.77%** |

Gaz-başına genişlik: **σ*/2(son) = 0.24493**, σ_X̃(λ=1)'den **+%1.19** —
Hkeskin'in kendisinden (−%2.57) daha yakın, λ-ailesi bandının tam ortası.

*(Kurtarmasız: ön-mühür S2 örneklem-dışı sapmayı "≤ %4" demişti; 7. bantta
−%7.12 çıktı. Ama bu son'a özgü değil — L085 aynı bantta −%6.53. τ > 0.77
bütün gazlarda A'nın dışıdır; §T2a.4'ün penceresi burada da geçerli.)*

### T3.2 ΔM'yi hangi ölçülür büyüklük öngörür? — **HİÇBİR MARJİNAL**

λ-ailesinden kurulan `M(x)` eğrileri gerçek gaza uygulanınca:

| öngörücü x | x(son) | M öngörü | ölçülen − öngörü |
|---|---|---|---|
| σ_ds | 0.40919 | 1.0090 | **+4.67% (+7.9σ)** |
| σ_X̃ | 0.23384 | 1.0056 | **+5.02% (+8.4σ)** |
| σ_Ĉ | 0.27303 | 1.0031 | **+5.29% (+8.9σ)** |
| *g_E* | *0.59821* | *1.0408* | *+1.47% (+2.6σ)* |
| *g_X* | *0.70473* | *1.0626* | *−0.61% (−1.1σ)* |

> **Dürüstlük uyarısı:** `g_X`'in λ-ailesindeki tüm yayılımı ±%0.3'tür ve
> λ ile TEKDÜZE bile değildir; log g_X'te kuadratik uyum **kötü koşullu**
> bir ekstrapolasyondur. `g_X` satırının "başarısı" anlamsızdır ve
> kimlik sayılmaz. `g_E` satırı da (+2.6σ) tam kapatmıyor.

### T3.3 170 §K3'ün ölçülen kırpma kapanışı, A·M dilinde

170 §K3.2 kuralı (e): ölçülen 3-bacak kırpma aktarımı
`ρ₃(son)/ρ₃(Hk) = 0.5356/0.5284 = 1.0136`, `(…)^{4/3} = 1.0182`. Bu kural
`c` üzerinde yazılmıştı; **M diline** çevirisi payda çarpanını geri verir:

```
M_kırpma = (ρ₃/ρ₃⁰)^{4/3} · W_X(son)/W_X(Hk) = 1.0182 · 1.0333 = 1.0521
```

| | değer | ΔM'nin kapanan payı | kalan açık |
|---|---|---|---|
| ölçülen | M = 1.0561 ± 0.0060 | — | — |
| kırpma kuralı | M_kırpma = **1.0521** | **%93** | +0.37% = **+0.7σ** |

*(c dilinde aynı kural %86 kapatıyordu — 170 §K3. M dilinde daha iyi
kapanıyor çünkü `W_X` payı ayrıştırılmış oluyor.)*

**Ama aynı kural λ ekseninde M dilinde DAHA SERT çöküyor:**

| gaz | M_kırpma | M ölçülen | sapma |
|---|---|---|---|
| L085 | 1.0717 | 1.0144 | −5.34% (**−10.6σ**) |
| L070 | 1.2323 | 1.0595 | −14.03% (**−24.4σ**) |
| L060 | 1.4180 | 1.1205 | −20.97% (**−43.3σ**) |

*(170 §K3 aynı çöküşü c dilinde −6.7σ / −10.1σ / −13.4σ olarak görmüştü;
M dili `W_X`'in kısmi telafisini kaldırdığı için çöküş 3 kat büyüyor.)*

> **HÜKÜM (T3).** Gerçek ζ gazının fazlası **saf bir SEVİYE olayıdır**:
> şekli (A) λ-ailesinin tam ortasında (σ*/2 = 0.24493, +%1.19), seviyesi
> (M) ise **hiçbir marjinalden öngörülemez** (σ_ds/σ_X̃/σ_Ĉ ile +7.9…+8.9σ).
> ΔM'nin %93'ünü kapatan tek büyüklük, **korelatör düzeyinde ÖLÇÜLEN**
> 3-bacak kırpma aktarımıdır — ama o kural λ ekseninde M dilinde
> **−10.6σ'dan −43.3σ'ya** çöker. **Yani gerçek gazın fazlası, aynı
> λ'daki iki gaz arasındaki kırpma farkıyla açıklanıyor; λ'yı
> değiştirmek bambaşka bir kanaldır (M'nin kendisi) ve o kanal hâlâ
> kimliksizdir.** ΔM'nin çarpan adresi ise nettir: **%45 Gram şişmesi
> (g_E) + %52 faz tutarlılığı (θ)** — ikisi de yaklaşık eşit paylı.

---

## T4 (bonus) — β KÖPRÜSÜ (`171j_T4.py`)

### T4.1 Kesim gazları AYNI A(τ) üstünde **DEĞİL**

`S/A2Hk − 1` (%) beş bantta, ve gaz-başına şekil genişliği:

| gaz | φ | b1 | b2 | b3 | b4 | b5 | rms | τ-eğimi | σ*/2 |
|---|---|---|---|---|---|---|---|---|---|
| λ-ailesi | 0 | — | — | — | — | — | 0.53–1.24% | −8…+17 | 0.2269–0.2493 |
| **son** | 0 | +0.18 | −1.35 | 0 | −0.34 | −1.01 | **0.77%** | −3.4 | 0.24493 |
| **K090** | 0.4153 | +0.90 | −0.75 | 0 | −0.81 | −0.66 | **0.70%** | −8.0 | 0.24855 |
| K070 | 0.9270 | +2.19 | −0.49 | 0 | −1.78 | −3.33 | 1.97% | −32.4 | 0.26834 |
| HA4 | 0.9098 | +3.37 | −0.35 | 0 | −3.64 | −7.86 | 4.16% | −65.7 | 0.29380 |
| E060 | 0.9595 | +3.97 | −2.02 | 0 | −8.71 | −17.99 | **9.16%** | −129.6 | 0.34136 |

> **HÜKÜM (T4-a).** **K090 aynı A(τ) üstündedir** (rms %0.70 — λ-ailesinin
> içinde), K070 sınırdadır (%1.97), **HA4 ve E060 kesinlikle değildir**
> (%4.2 ve %9.2; σ*/2 +%21 ve +%41). Kesim, bant şeklini **dikleştirir**
> ve dikleşme kesim şiddetiyle büyür. `A(τ)`'nun λ-değişmezliği
> **kesim ekseninde geçerli değildir.**

### T4.2 β = 0.2175 hangi çarpanın işi — **SEVİYENİN, ŞEKLİN DEĞİL**

Sınav: `KALİB_b(kesim) = KALİB_b(Hk)·(g_cal/g₀)·(1 − βφ)` (yani
`θ/θ₀ = 1 − βφ`, 168 §A2.7). Artıklar (%):

| gaz | φ | 1−βφ | θ/θ₀ ölçülen | b1 | b2 | b3 | b4 | b5 | rms | τ-eğimi |
|---|---|---|---|---|---|---|---|---|---|---|
| K090 | 0.4153 | 0.9097 | **0.9080** | +1.26 | +0.09 | −0.12 | −0.89 | −1.26 | **0.89%** | −15.2 |
| K070 | 0.9270 | 0.7984 | **0.8038** | +3.62 | +1.39 | +0.91 | −0.86 | −1.56 | 1.95% | −33.3 |
| HA4 | 0.9098 | 0.8021 | **0.7975** | +4.80 | +1.54 | +0.98 | −2.53 | −7.20 | 4.11% | −71.6 |
| E060 | 0.9595 | 0.7913 | **0.7945** | +10.32 | +4.51 | +5.71 | −3.23 | −13.47 | 8.38% | −141.7 |
| *(kıyas)* L060 | 0 | 1.0000 | *1.0493* | +5.28 | +5.42 | +5.81 | +5.04 | +3.14 | 5.02% | *(DÜZ)* |
| *(kıyas)* son | 0 | 1.0000 | *1.0290* | +3.88 | +2.78 | +3.20 | +2.92 | +1.70 | 2.98% | *(DÜZ)* |

Bu tablo iki ekseni **temiz biçimde ayırıyor**:

* **λ ekseni**: `1−βφ = 1` (φ = 0); artık **DÜZ** (+5.3/+5.4/+5.8/+5.0/+3.1)
  — yani saf bir **SEVİYE** borcu, tam olarak `θ`'nın λ-hareketi (ve o da
  M'nin yarısı). Şekil borcu yok.
* **Kesim ekseni**: `1−βφ` seviyeyi ‰2–7 içinde **kapatıyor**
  (θ/θ₀ ölçülen ile 1−βφ farkı K090 −%0.19, K070 +%0.68, HA4 −%0.57,
  E060 +%0.40) ama geriye **tekdüze bir τ-EĞİMİ** kalıyor
  (−15 → −142 per birim τ). Bu, β yasasının **görmediği** şekil borcudur.

> **HÜKÜM (T4-b).** **β = 0.2175 SEVİYE çarpanının (θ, yani M'nin
> içindeki faz tutarlılığı) yasasıdır; ŞEKİL çarpanının (A) yasası
> değildir.** Kesim ekseninde `θ = 1 − βφ` seviyeyi ‰7 içinde verir ama
> bant şekli aynı anda dikleşir ve bu dikleşmeyi β hiç tarif etmez.
> Yani kesim ekseni **İKİ boyutludur**: bir seviye (β·φ) ve bir şekil
> (A'nın dikleşmesi). 168 §A2.7'nin kesim yasası bunun yalnız birinci
> boyutudur.

### T4.3 Kesim ekseninde ortak şekil var mı? — **φ tek değişken DEĞİL**

| gaz | φ | α_g | α − α_λ | σ*/2 | kesim biçimi |
|---|---|---|---|---|---|
| Hkeskin | 0.0000 | 1.0977 | −0.0449 | 0.23582 | keskin τ≤1.00 |
| K090 | 0.4153 | 1.2195 | +0.0769 | 0.24855 | keskin τ≤0.90 |
| K070 | 0.9270 | 1.4213 | +0.2787 | 0.26834 | keskin τ≤0.70 |
| HA4 | 0.9098 | 1.7038 | +0.5613 | 0.29380 | erfc 0.68/0.125 |
| E060 | 0.9595 | 2.3001 | +1.1575 | 0.34136 | erfc 0.60/0.125 |

Δα ile φ korelasyonu (4 kesim gazı) **r = +0.673**; ama sıralama bozuluyor:
φ(K070) = 0.9270 > φ(HA4) = 0.9098 iken α(K070) = 1.4213 **<** α(HA4) = 1.7038.

> **HÜKÜM (T4-c).** Kesim ekseninde de **tek-değişkenlik yok**: şeklin
> dikleşmesi φ'nin (boş çizgi kesri) tek-değişkenli fonksiyonu değildir;
> kesimin BİÇİMİ (keskin ↔ erfc yumuşaklığı) bağımsız olarak giriyor.
> Bu, §T2b.2'nin λ ekseninde bulduğu tek-değişkensizliğin kesim
> eksenindeki ikizidir: **hiçbir eksende `c` tek bir skaler etikete
> indirgenmiyor.**

---

## 4. HÜKÜM TABLOSU

| kapı | soru | **hüküm** |
|---|---|---|
| **T2a** | A(τ)'nun kimliği | **BULUNDU ve MÜHÜRLENDİ**: `A(τ) = exp(−2π²τ²σ_X̃(λ=1)²)`, **sıfır serbest parametre**, 9-bant rms %1.49; 1-parametreli en iyi uyumu yener. Bütün rakipler örneklem-dışı öldü. |
| T2a′ | `W_amp` var mı | **YOK.** DW-tam (W_amp·W_X) beş bantta ±%9.4, uzantıda +%16/+%18 ile öldü — 168 §A2.8'in "W_amp çift sayım" hükmü şekil düzeyinde doğrulandı. |
| T2a″ | 165'in F'i, tarak-sayımı | **İKİSİ DE ÖLDÜ.** Tarak-sayımı **işarette** (×6 artış ister, ×0.80 düşüş var); 165'in F ailesi uzantıda +%12/+%24. |
| **T2b** | M yalnız σ_X̃ (ya da σ_ds/σ_Ĉ) üzerinden mi | **HAYIR — hiçbiri.** Aynı σ_X̃'te dokuz gaz −%17.8 ile +%5.0 arasına yayılıyor (−17.7σ … +8.4σ). `M`, marjinalin değil **λ EKSENİNİN ETİKETİDİR.** |
| T2b′ | M'nin çarpan adresi | **g_E (Gram şişmesi) + θ (faz tutarlılığı), yaklaşık eşit paylı**; g_X λ-değişmez. Üç eksen üç ayrı çarpanı sürüyor (λ→g_E+θ, kesim→θ, pencere→g_X²). |
| T2b″ | tek-taraflılık | **AYAKTA, iki örneklem-dışı noktayla**: M(1.15) = 1.0022 ± 0.0054, **M(1.30) = 0.9832 ± 0.0157** (1'den −1.1σ). Eşik λ_c = 1.0109 uydurulmadı (inşa doyumu). |
| **T2b‴** | M(λ)'nın KİMLİĞİ | **BULUNAMADI.** Türetimli beş ailenin **hepsi** λ = 0.50'de öldü (+9.4σ … +15.3σ). Ayakta kalan tek eğri EMPİRİK log-λ kuadratiğidir (+0.4σ) — kimlik değildir. **AÇIK BORÇ.** |
| **T2c** | λ = 0.50 ve λ = 1.30 | ön-kayıtlı ölçüldü; c(0.50) = 0.3800 ± 0.0072, c(1.30) = 0.4540 ± 0.0136. |
| **T2c′** | `c`'nin tabanı var mı | **YOK.** `c(λ)`'nın λ* = **0.6487**'de bir **MİNİMUMU** var (ν(0.60→0.50) = 1.489 ± 0.063, 1'i **+7.7σ** aşıyor). 170 §K2'nin "taban" okuması öldü. |
| T2c″ | 4/π² | `c(1.30) = 4/π² + %12.0 (+3.6σ)`. Limit değil; λ=1.00'in tesadüfü — ikinci örneklem-dışı doğrulama. |
| T2c‴ | çarpanlaşmanın penceresi | **λ ∈ [0.50, 1.15] ve τ ≲ 0.70.** λ = 1.30'da şekil kanadı kırılıyor (σ*/2 = 0.1866, bandın %18 altında); τ ≈ 0.74'te gazlar arası yayılım %8. |
| **T3** | gerçek gazın ΔM'si | **Saf SEVİYE olayı.** Şekli A'nın tam ortasında (+%1.19), seviyesi hiçbir marjinalden öngörülemiyor (+7.9σ…+8.9σ). ΔM'nin **%93'ünü** ölçülen 3-bacak kırpma aktarımı kapatıyor (+0.7σ artık) — ama aynı kural λ ekseninde M dilinde −10.6σ…−43.3σ çöküyor. Çarpan payları: g_E %45, θ %52. |
| **T4** | kesim gazları aynı A'da mı | **K090 EVET** (rms %0.70), K070 sınırda, **HA4 ve E060 HAYIR** (%4.2, %9.2). |
| T4′ | β = 0.2175 kimin işi | **SEVİYENİN (θ ⊂ M), ŞEKLİN (A) DEĞİL.** `θ = 1−βφ` seviyeyi ‰7 içinde kapatıyor; geriye β'nın hiç görmediği tekdüze bir τ-eğimi kalıyor (−15 → −142). Kesim ekseni **iki boyutlu**. |
| T4″ | kesimde tek-değişkenlik | **YOK.** Δα ile φ korelasyonu r = +0.673 ve sıralama bozuluyor (K070 ↔ HA4): kesimin BİÇİMİ bağımsız değişken. |

### 4.1 Tek cümlelik hüküm

> **A(τ) çözüldü, M(λ) çözülmedi ve `c(λ)`'nın tabanı yok.** λ-değişmez
> bant şekli, **fiziksel (λ = 1) merdivenin X̃ Debye–Waller çarpanının ta
> kendisidir** — sıfır serbest parametreyle, örneklem-dışı iki bantta
> mühürlü. Seviye çarpanı `M(λ)` ise ne bir marjinalin fonksiyonudur
> (aynı σ_X̃'te dokuz gaz %23'lük bir yelpazeye yayılıyor) ne de beş
> türetimli aileden birine uyar (hepsi λ = 0.50'de 9–15σ ile öldü);
> yalnız `λ ≥ 1'de düz` özelliği iki örneklem-dışı noktayla ayakta.
> Ve λ = 0.50 gazı 170'in son resmini de bozdu: `c` bir tabana oturmuyor,
> **λ* = 0.649'da bir minimumdan geçip yeniden yükseliyor** (+7.7σ).

---

## 5. DENETİM ve DÜRÜSTLÜK NOTLARI

1. **Ön kayıtların zaman damgaları.** A(τ)'nun 9-bant öngörüleri
   `A_ONKAYIT.json`'a **23:44:17**'de yazıldı, hüküm 171d'de verildi.
   Yeni gazların c/M öngörüleri `ONKAYIT_L050.json` **00:00:05** ve
   `ONKAYIT_L130.json` **00:03:56**'da yazıldı; ölçümler 00:05 ve
   00:08'de başladı. Bu dosyalara sonradan dokunulmadı.
2. **`171h`in ön-mührü KÖR DEĞİLDİR.** `171g`in ekran çıktısında L050'nin
   gaz-düzeyi `c_ampX`'i (0.4369, 167'nin kendi 6-bant penceresi)
   görüldü ve V1 ön-mührü buna dayanarak yazıldı; betiğin docstring'inde
   böyle işaretlidir. Hükmü belirleyen ön kayıt (00:00:05) bundan önce
   mühürlenmiştir.
3. **Iskalayan ön-mühürler (kurtarmasız).**
   * `171b`: σ_ds(1.30) = 0.48890 aralığın (**0.477–0.487**) üstünde;
     σ_X̃(1.30) +%1.1, σ_Ĉ(1.30) +%1.7 üstünde. ΔG<0 kesri(0.50) = 0.0106,
     öngörü 0.019 ± 0.005 idi. min Δz iki gazda da ıskaladı.
   * `171c` P4/P6: A1'in artığı öngörülenden büyük (±%9.4 ↔ ±%4-5),
     A3s'inki de (%2.5 ↔ %0.5-1.5).
   * `171d` Q3: 6. bantta gaz-gaz yayılımı **%8.08**, öngörü ≤%3 idi.
     Q4: A2Hk'nın 9-bant rms'i %1.49, öngörü ≤%1.2 idi.
   * `171e` N1: bant-ortalamalı M değerleri orta-bant öngörüsünün %0.3–0.8
     altında (tanım farkı; 1σ içinde ama tam isabet değil).
   * `171f` R1: λ=1.30'un "5–10σ'lık ayrım" öngörüsü **ıskaladı** —
     gerçekleşen ayrım 1.1σ (M2/M3) ↔ 2.8σ (M1), çünkü σ_M = 0.0157
     (öngörü ~0.006) ve fazlası **şekil sürüklenmesinden** geliyor.
     R3: c(0.50) öngörü aralığı 0.360–0.372 idi, ölçüm **0.3800**.
     R4: σ_X̃(1.30) +%2.3 iskaladı.
   * `171i` S2: son'un 7. bant sapması −%7.12, öngörü ≤%4 idi (ama bu
     bant bütün gazlarda A'nın dışında).
   * `171h` V3: L130'un şekil rms'i %3.62 (öngörü ≤%2.5) ve σ*/2 = 0.18661
     (öngörü ≈0.220).
4. **ν'nün hata çubukları** yalnız bant-içi jackknife'tandır; bantlar arası
   saçılım iki gaz arasında büyük ölçüde ORTAKTIR (aynı bant tanımı, aynı
   çizgiler) ve oranda sadeleşir. Bu yüzden ν'nün hatası c farkınınkinden
   küçüktür ve "c(0.50) > c(0.60)" hükmü doğrudan c farkında +1.1σ,
   ν'de **+7.7σ**'dır. İki sayı çelişmez; ikincisi doğru istatistiktir.
5. **9-bant hakeminin zayıf noktası:** 7. bant (τ ≈ 0.777) yalnız İKİ
   λ-gazında (Hkeskin, L085) sağlıklıdır; oradaki A'nın hata çubuğu
   (±%3) iki noktalık bir saçılımdan gelir. 6. bant beş gazın hepsinde
   vardır ve asıl hakem odur.
6. **E060 gürültülüdür** (σ_tot = 0.0165, M hatası 0.0415). T4'ün
   E060 satırları yön gösterir, tek başına hüküm taşımaz. Kesim
   hükümleri K090/K070/HA4 üçlüsüne dayanmaktadır.
7. **R_bant** (sadakat ölçütü) bütün gazlarda 1'in üstündedir ve λ ile
   düşer: L130 1.05–1.16, L115 1.14–1.28, Hkeskin 1.22–1.40,
   L060 1.47–1.74, **L050 1.50–1.82**. L050'nin değerleri ailenin
   düzgün devamıdır, aykırı değildir; yine de not edilir.
8. **Hiçbir ölçüm/çözüm parçası kopyalanmadı** (yalnız import).
   **Git'e dokunulmadı.**

---

## 6. NE KAZANILDI, NE KALDI

**Kazanıldı**

* `A(τ)`'nun kimliği: `W_X(τ; σ_X̃(λ=1))`, sıfır parametre, örneklem-dışı
  mühürlü. `c`'nin τ-bağımlılığı artık bir yasadır, uyum değil.
* `M(λ)`'nın çarpan adresi: `g_E·θ` (yaklaşık eşit paylı); `g_X` λ-değişmez.
* Üç eksenin üç ayrı çarpanı sürdüğünün tablosu (λ→g_E+θ, kesim→θ,
  pencere→g_X²) — 168 §BONUS-iii'ün "T/4 borcu g_X'tedir" hükmü
  bu tabloda kendiliğinden yerine oturuyor.
* β = 0.2175'in yeri: **seviyenin (θ) yasası**, şeklin değil; kesim
  ekseninin ikinci (şekil) boyutu ölçüldü ve β'nın onu görmediği gösterildi.
* Gerçek gazın fazlasının **saf seviye** olduğu; şeklinin λ-ailesinin
  tam ortasında durduğu.
* `c(λ)`'nın **U şekli** ve minimumunun yeri (λ* = 0.6487) — ve ν(λ)'nın
  tekdüze olmadığı (λ ≈ 1.08'de minimum).

**Kaldı (açık borçlar)**

1. **`M(λ)`'nın kimliği.** Beş türetimli aile de öldü. Yeni fikir gerekiyor:
   `M`'nin yükselen kolu λ = 0.60 → 0.50'de eşik ailelerinin izin
   verdiğinden **%7 daha hızlı** dikleşiyor. Bir sonraki doğal hamle
   `g_E` ve `θ`'yı **ayrı ayrı** modellemektir (M'nin kendisini değil):
   `g_E` Gram matrisinin ölçülür bir fonksiyonelidir (168 §A1.5),
   `θ` ise faz tutarlılığıdır — ikisinin λ-yasaları farklı olabilir.
2. **λ = 0.40 (ya da 0.35) gazı**: `c`'nin U'sunun sol kolu ne kadar
   yükseliyor? `M`'nin kimliği için en bilgilendirici tek nokta budur
   (M(0.40) eşik ailelerinden ~%15 sapmalı).
3. **A'nın λ = 1.30'da kırılması**: şekil genişliği neden λ ≥ 1.2'de
   çöküyor (σ*/2: 0.2269 → 0.1866)? Bu, çarpanlaşmanın kendi sınırıdır
   ve mekanizması bilinmiyor.
4. **τ > 0.70 bandı**: A·M ayrışması orada zayıflıyor (%8 yayılım).
   9-bant penceresinin hükümleri bu yüzden 5-bant hükümlerinden zayıftır.
5. **β'nın türetimi** 168 → 171 boyunca hâlâ açıktır; 171 yalnız
   **nerede yaşadığını** (θ, seviye) ve **neyi görmediğini** (şekil)
   söyledi.

*(Figür: `171_nu.png` — (a) A(τ) yarışı ve örneklem-dışı noktalar,
(b) M(λ) ve ön-kayıtlı öngörüler, (c) c(λ)'nın U'su, (d) şekil
genişliğinin geçerlilik penceresi.)*
