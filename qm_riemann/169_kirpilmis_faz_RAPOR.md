# 169 — DOYMUŞ ÇARPANIN KAPALI FORMU: c = (2/π)² = 4/π² ÜÇ-KAPILI SINAV
(3 Eylül 2026)

Hipotez (KALEM_KIRPILMIS_FAZ_03EYL2026.md):

> **H-C1: c = (2/π)² = 4/π² = 0.4052847**, doymuş (tam sarılmış) fazlı
> korelatörün İŞARET-korelasyonuna inmesinden gelen evrensel geometrik sabit.

Scriptler: `169_configs/169_k1.py` (K1 — yalnız `C_*.json`, yeni ölçüm yok),
`169_configs/169_uye.py` (üye sınavı: iki ölçüt),
`169_configs/169_k2.py` (K2a sarılma + alan-düzeyi arcsine),
`169_configs/169_k2b.py` (K2b — BÖLÜNMÜŞ MERDİVEN, çok-bacaklı kırpma),
`169_configs/169_k3_onkayit.py` (K3b ÖN KAYIT — ölçümden ÖNCE),
`169_configs/169_k3_bos.py` (K3a — β'nın adresi: boş-çizgi deneyi).
Yeni ölçümler: `167_insa.py L070` (gaz inşası) + `167_olcum.py L070`.
Ham çıktılar `scratchpad/169`, girdiler `scratchpad/167`.

---

## 0. TEK CÜMLELİK HÜKÜM

> **MEKANİZMA DOĞRULANDI, EVRENSELLİK ÇÜRÜTÜLDÜ — H-C1 MÜHÜRSÜZ.**
> Ç4 rejiminde birikmiş faz `2πτĈ` mod 2π **ölçülerek** tam sarılıyor
> (Kuiper `√N·V` = 1.10–1.68, düzgünlük eşiği 2.00; aynı ölçüt tek
> merdiven çizgilerinde 78–500 veriyor — **beş mertebe kontrast**), ve
> korelatörün işaret-kırpma aktarımı Gauss ailesini `(2/π)^{n/2}` bacak
> bacak %3.8 içinde izliyor (1/2/3 bacak: 0.8293 / **0.6586** / 0.5284
> ↔ 0.7979 / **0.6366 = 2/π** / 0.5079). Kesim yasasının MUTLAK φ→0
> kesişimi `c₀ = 0.4057 ± 0.0016` = **4/π² + 0.23σ**. **Ama** `c` evrensel
> bir sayı değil: aynı merdivenin genlik ölçeği λ = 1.00 → 0.85 → **0.70**
> yapılınca `c` = 0.4035 → 0.3817 → **0.3690** ve **ÖN KAYITLI**
> örneklem-dışı öngörü **−5.9σ**'da düşüyor (rakip 168-λ-yasası da +3.6σ'da).
> Üstelik üye seçimi (`W_X`) — 4-hanelik uyumun tümüyle dayandığı seçim —
> λ ekseninde çürüyor (§K1.5). Sabitin **ailesi** kazanıldı, **kendisi**
> değil.

---

## K1 — HASSASİYET

### K1.0 Üye ve pencere (tanım)

168 §A2.8'in hükmü: `c` **`W_X` üyesinin** sayısıdır,

```
c_WX ≡ KALİB_u2 / W_X ,      W_X(τ) = ⟨e^{−2πiτ X̃0}⟩   (TAM karakteristik fn)
```

(`C_*.json`'daki `c_u2` alanı `KALİB_u2/(W_amp·W_X)`'tir; `W_amp` çift
sayımdır — 168 §A2.8.) Jackknife: bant içi 8 grup (160'ın round-robin'i),
silme-1; her silmede `W_amp`/`W_X` de `gp`-ağırlığıyla yeniden hesaplanır
(oran tahmincilerinde doğrusal hata yayılımı YOK). Bant birleştirmesi
`166_T1.bant_agg`, jackknife `166_T1._jk` — **kopyalanmadı, import edildi**.

Ortak pencere (168 §A2.6 ile aynı): `lo ∈ [0.52, 0.68]`, `SNR ≥ 3`,
`R_bant ≥ 0.98`, `τ_eff < 0.85` — beş bant.

### K1.1 168'İN TABLOSUNUN YENİDEN ÜRETİMİ (dürüstlük şartı)

Bağımsız kod yolundan (169_k1.py), 168 §A2.6/§A2.8'in satırları:

| gaz | c_ampX (bu) | 168 | g_cal (bu) | 168 | θ (bu) | 168 | c_WX 9-bant (bu) | 168 §A2.8 |
|---|---|---|---|---|---|---|---|---|
| Hkeskin | 0.5774 | 0.5774 | 0.2895 | 0.2895 | 0.8884 | 0.8884 | **0.4051** | 0.4051 |
| son | 0.5664 | 0.5664 | 0.2971 | 0.2971 | 0.9142 | 0.9142 | **0.4053** | 0.4053 |
| L085 | 0.5171 | 0.5171 | 0.2912 | 0.2912 | 0.8959 | 0.8959 | **0.3751** | 0.3751 |
| K090 | 0.5393 | 0.5393 | 0.2877 | 0.2877 | 0.8067 | 0.8067 | **0.3709** | 0.3709 |
| K070 | 0.5154 | 0.5154 | 0.3128 | 0.3128 | 0.7142 | 0.7142 | **0.3480** | 0.3480 |
| HA4 | 0.5061 | 0.5061 | 0.3194 | 0.3194 | 0.7086 | 0.7086 | **0.3345** | 0.3345 |
| E060 | 0.5216 | 0.5216 | 0.3433 | 0.3433 | 0.7059 | 0.7059 | **0.3701** | 0.3701 |

> **YENİDEN ÜRETİM TAM** (dört hane). 168 §A2.8'in `c_WX` sayıları DOKUZ
> (lo ≤ 0.80) banttan, §A2.6/§A2.7'ninkiler ORTAK PENCEREDEN (lo ≤ 0.68)
> geliyor; ikisi de yeniden üretildi. Aşağıda **ortak pencere** kullanılır
> (kesim/pencere eksenleri orada karşılaştırılabilir); fark Hkeskin'de
> 0.4035 ↔ 0.4051'dir (%0.4) ve raporda hep ikisi de yazılıdır.

### K1.2 c'NİN JACKKNIFE'LI DEĞERİ ve 4/π²'YE σ-UZAKLIĞI

Ortak pencerede, beş bant (bant-bant tablolar `scratchpad/169/K1.json`):

| gaz | eksen | **c_WX ± σ_jk ± σ_bant** | %fark (4/π²) | **σ_jk** | **σ_tot** |
|---|---|---|---|---|---|
| **Hkeskin** | taban | **0.4035 ± 0.0012 ± 0.0013** | **−0.43%** | −1.39 | **−0.97** |
| **son** (GERÇEK ζ) | taban | **0.4122 ± 0.0013 ± 0.0018** | **+1.71%** | +5.15 | **+3.06** |
| L085 (λ=0.85) | λ | 0.3817 ± 0.0012 ± 0.0025 | −5.82% | −20.5 | −8.46 |
| **L070 (λ=0.70)** ‡ | λ | 0.3690 ± 0.0011 ± 0.0057 | **−8.95%** | −33.5 | **−6.29** |
| K090 | kesim | 0.3699 ± 0.0011 ± 0.0009 | −8.74% | −31.1 | −23.9 |
| K070 | kesim | 0.3480 ± 0.0014 ± 0.0028 | −14.15% | −41.3 | −18.3 |
| HA4 | kesim | 0.3492 ± 0.0017 ± 0.0071 | −13.83% | −33.4 | −7.68 |
| E060 | kesim | 0.3701 ± 0.0037 ± 0.0160 | −8.69% | −9.54 | −2.14 |
| HkT2a/b | pencere | 0.3871 / 0.3708 | −4.5% / −8.5% | | −6.1 / −15.1 |
| HkT4a/b | pencere | 0.3309 / 0.3314 | −18.4% / −18.2% | | −16.3 / −12.1 |

‡ = 169'da İNŞA EDİLDİ ve ÖRNEKLEM-DIŞI ölçüldü (§K3.2).

Dokuz bantta (168 §A2.8'in penceresi): `c(Hkeskin) = 0.4051` (**−0.04%**),
`c(son) = 0.4053` (**+0.00%**).

> **HÜKÜM (K1a).** İki TABAN gazında (biri GERÇEK zeta sıfırları) ölçülen
> `c` 4/π²'ye **‰0.4 ile ‰0.0** uzaklıkta; ortak pencerede sapma
> `−0.97σ_tot` (Hkeskin) ve `+3.06σ_tot` (son). **Hassasiyet sınamaya
> yetiyor ve taban gazları sınavı geçiyor.** Ama `c` ölçülen hâliyle
> EVRENSEL DEĞİL: kesim/pencere/λ eksenlerinde %5–18 kayıyor. Sınav bu
> yüzden K1(b)'ye (kesim yasasının kesişimi) taşınır.

### K1.3 KESİM YASASININ φ→0 KESİŞİMİ (mutlak normalizasyon)

168 §A2.7'nin yasası `θ/θ₀ = 1 − β φ` (φ = BOŞ ÇİZGİ KESRİ). Aynı yasayı
**`c` üzerinde ve MUTLAK normalizasyonla** (referans oranı yok, kesişim
serbest parametre) yazıyoruz:

```
c(φ) = c₀ (1 − β φ)     — ağırlıklı en-küçük-kareler, ağırlık 1/σ_tot²
```

| gaz | φ | c_WX ± σ | c_WX·(g₀/g_cal) |
|---|---|---|---|
| Hkeskin | 0.0000 | 0.4035 ± 0.0018 | 0.4035 |
| K090 | 0.4153 | 0.3699 ± 0.0015 | 0.3722 |
| HA4 | 0.9098 | 0.3492 ± 0.0073 | 0.3165 |
| K070 | 0.9270 | 0.3480 ± 0.0031 | 0.3220 |
| E060 | 0.9595 | 0.3701 ± 0.0165 | 0.3120 |

| uyum | **c₀ (φ→0 kesişimi)** | β | 4/π²'ye uzaklık |
|---|---|---|---|
| ham `c_WX` | 0.4002 ± 0.0016 | 0.1573 | −1.25% = **−3.19σ** |
| **`g_cal`-düzeltmeli** (168 §A2.6'nın ayrıştırması) | **0.4057 ± 0.0016** | **0.2143** | **+0.09% = +0.23σ** |

> **HÜKÜM (K1b).** Kesim ekseninin dört gazından kurulan doğrunun φ→0
> kesişimi, `g_cal` sızıntısı 168 §A2.6'nın ölçülmüş ayrıştırmasıyla
> çıkarıldığında **c₀ = 0.4057 ± 0.0016 = 4/π² + 0.23σ** — yani sıfır
> serbest parametreli hedefi **‰1 içinde** vuruyor. Aynı uyumun eğimi
> `β = 0.2143`, 168'in bağımsız olarak ölçtüğü `β = 0.2175 ± 0.0054` ile
> **0.6σ** uyumlu. `g_cal` düzeltmesi YAPILMAZSA kesişim 0.4002 (−3.2σ):
> **düzeltmenin gerekli olduğu ve 168'in ayrıştırmasına bağımlı olduğumuz
> AÇIK bir borçtur.**

### K1.4 TAM DÜZELTMEDEN SONRA ÇÖKME (H-C1'in gerçek sınavı)

`c* = c_WX · (g₀/g_cal) / (1 − βφ)`, β = 0.2143. H-C1: **her gazda
c\* = 4/π².**

| gaz | eksen | φ | g₀/g_cal | c_WX | **c\*** | c\*/(4/π²)−1 |
|---|---|---|---|---|---|---|
| Hkeskin | taban | 0.000 | 1.0000 | 0.4035 | 0.4035 | **−0.43%** |
| son | taban | 0.000 | 0.9743 | 0.4122 | 0.4016 | **−0.90%** |
| K070 | kesim | 0.927 | 0.9253 | 0.3480 | 0.4018 | **−0.86%** |
| K090 | kesim | 0.415 | 1.0063 | 0.3699 | 0.4085 | **+0.80%** |
| HA4 | kesim | 0.910 | 0.9064 | 0.3492 | 0.3932 | **−2.98%** |
| E060 | kesim | 0.960 | 0.8432 | 0.3701 | 0.3928 | **−3.07%** |
| HkT2a | pencere | 0.000 | 1.0718 | 0.3871 | 0.4148 | +2.36% |
| HkT2b | pencere | 0.000 | 1.0917 | 0.3708 | 0.4048 | −0.11% |
| HkT4a | pencere | 0.000 | 1.3110 | 0.3309 | 0.4337 | **+7.02%** |
| HkT4b | pencere | 0.000 | 1.3418 | 0.3314 | 0.4447 | **+9.73%** |
| **L085** | **λ** | 0.000 | 0.9941 | 0.3817 | 0.3794 | **−6.37%** |
| **L070** | **λ** | 0.000 | 0.9694 | 0.3690 | 0.3577 | **−11.73%** |

> **HÜKÜM (K1c).** Oniki gazdan **yedisi ±3% içinde** 4/π²'ye çöküyor
> (taban 2, kesim 4, pencere 1). Kalan beşi iki bilinen borcu tekrar
> ediyor: **T/4 dilimleri** (+%7/+%10 — 168 §6.2'nin açık borcu; `(1+2n/N)`
> düzeltmesi tam değil) ve **λ ekseni** (−%6.4 / **−%11.7**). λ ekseni
> sadece bir sapma değil, H-C1'in içinde bir ÇELİŞKİDİR ve K3(b) onu
> örneklem-dışı sınava sokar (§K3.2).

### K1.5 ÜYE SINAVI — İKİ ÖLÇÜT ÇATIŞIYOR (169'un yeni bulgusu)

168 §A2.8 üyeyi TEK ölçütle seçti (bant τ_eff eğimi + `W_amp`'ın çift
sayım olduğu türetim argümanı). 169 ikinci bir ölçüt ekliyor: **λ
ekseninde değişmezlik** (üç gaz aynı merdiven, yalnız genlik ölçeği
farklı ⇒ aynı fizik). `169_configs/169_uye.py`:

| üye | τ-eğimi (7 bant) | artık rms | c(λ=1.00) | c(0.85) | c(0.70) | **λ-menzil** | 4/π²'ye |
|---|---|---|---|---|---|---|---|
| K0 (sabit) | −1.474 | 0.0121 | 0.2572 | 0.2609 | 0.2725 | 5.95% | −36.5% |
| W_amp | −0.228 | 0.0085 | 0.3680 | 0.3534 | 0.3473 | 5.96% | −9.2% |
| **√W_X** | −0.695 | 0.0097 | 0.3222 | 0.3156 | 0.3171 | **2.09%** | −20.5% |
| **W_X** (168'in üyesi) | **+0.085** | 0.0080 | 0.4035 | 0.3817 | 0.3690 | 9.35% | **−0.43%** |
| √(W_amp·W_X) | **−0.072** | 0.0082 | 0.3853 | 0.3673 | 0.3580 | 7.64% | −4.9% |
| W_pos | +0.698 | 0.0076 | 0.4724 | 0.4353 | 0.4097 | 15.3% | +16.6% |
| W_amp·W_X | +1.331 | 0.0081 | 0.5774 | 0.5171 | 0.4703 | 22.8% | +42.5% |
| W_X^{3/2} | +0.865 | 0.0076 | 0.5055 | 0.4617 | 0.4295 | 17.7% | +24.7% |

> **HÜKÜM (K1d — YENİ SERT KISIT).** Hiçbir üye iki ölçütü BİRDEN
> sağlamıyor. τ-eğimi ölçütü `W_X`/`√(W_amp W_X)` çiftini seçiyor
> (λ-menzili %7.6–9.4); λ-değişmezliği `√W_X` seçiyor (τ-eğimi −0.695,
> yani bant bağımlılığı 8 kat büyük) ve o üyenin mutlak değeri
> **0.3222 — 4/π²'nin %20 altında.** **`c = 4/π²`'nin 4-hanelik uyumu
> tamamen `W_X` üye seçimine bağlıdır ve bu seçim λ ekseninde çürüyor.**

---

## K2 — MEKANİZMA

### K2.1 SARILMA (tam-kırpılmış faz) — TANIM ve SINAV

165 §7 / 168 §A1.5'in geometrisi:

```
s_n = s̄ + ḡ(n−n̄) + ḡ·Ĉ_n ,   ḡ = 2π/L ,   Ĉ = kübik-trendsiz Σ_{k<n} X̃0_k
ν frekansının taşıdığı SİTE-SEĞİRME fazı:  ϑ_n(ν) = ν ḡ Ĉ_n = 2π τ_ν Ĉ_n
```

Ölçütler (N = 299 998 örnek): `R(τ) = |⟨e^{2πiτĈ}⟩|` (Rayleigh bileşke
uzunluğu; düzgünde ≈ √π/2/√N = 0.0016) ve **Kuiper `√N·V`** (düzgünde
medyan ≈ 1.62; %1 red eşiği 2.001).

`σ_Ĉ = 0.27768`. Model bacaklarının genlik-ağırlıklı τ'su: `⟨τ⟩_E = 0.7054`,
`⟨τ⟩_X = 0.7166` ⇒ Ç4'ün dört-frekans toplamı `τ_Q + ⟨τ⟩_E + 2⟨τ⟩_X`.

| etiket | τ | **R(τ)** | **√N·V (Kuiper)** | düzgün mü? |
|---|---|---|---|---|
| KONTROL q=2 | 0.0576 | 0.99496 | 499.5 | **HAYIR** |
| KONTROL q=7 | 0.1618 | 0.96088 | 427.5 | **HAYIR** |
| KONTROL q=1747 | 0.6206 | 0.54601 | 192.0 | **HAYIR** |
| merdiven ucu (tek çizgi) | 0.9500 | 0.22515 | 78.4 | **HAYIR** |
| bant τ_Q=0.54 (TEK taşıyıcı) | 0.5400 | 0.63480 | 226.4 | **HAYIR** |
| bant τ_Q=0.58 | 0.5800 | 0.59080 | 209.0 | **HAYIR** |
| bant τ_Q=0.62 | 0.6200 | 0.54668 | 192.2 | **HAYIR** |
| bant τ_Q=0.66 | 0.6600 | 0.50286 | 176.2 | **HAYIR** |
| bant τ_Q=0.70 | 0.7000 | 0.45974 | 160.7 | **HAYIR** |
| **bant τ_Q=0.54, Ç4 (4 frekans)** | **2.6787** | **0.00247** | **1.60** | **EVET** |
| **bant τ_Q=0.58, Ç4** | 2.7187 | 0.00246 | 1.58 | **EVET** |
| **bant τ_Q=0.62, Ç4** | 2.7587 | 0.00245 | 1.55 | **EVET** |
| **bant τ_Q=0.66, Ç4** | 2.7987 | 0.00244 | 1.65 | **EVET** |
| **bant τ_Q=0.70, Ç4** | 2.8387 | 0.00243 | 1.68 | **EVET** |
| Ç4 alt sınır (τ_Q+3·0.5) | 2.0400 | 0.00107 | 1.10 | **EVET** |

> **HÜKÜM (K2a GEÇTİ, KONTRASTLA).** Ç4'ün yaşadığı dört-frekans
> rejiminde (τ_eff = 2.0–2.8) birikmiş faz **mod 2π TAM DÜZGÜN**:
> Kuiper `√N·V = 1.10–1.68`, düzgünlüğün medyanının (1.62) etrafında,
> %1 eşiğinin (2.00) ALTINDA; bileşke uzunluk `R ≈ 0.0011–0.0025`, yani
> düzgün-null'un (0.0016) mertebesinde. **Aynı ölçüt aynı veride tek
> çizgilerde EZİCİ biçimde düzgünlüğü reddediyor** (`√N·V = 78–500`,
> R = 0.23–0.99). Kontrast beş mertebe: **"faz tam sarılıyor" bir varsayım
> değil, ölçülmüş bir olgudur** — ve yalnız Ç4 rejiminde.

### K2.2 İŞARET-AKTARIMI — TANIMLAR

`clip(F) ≡ σ_F · sgn(F)` (VARYANS-EŞLEŞMİŞ kırpma). Gauss alanı için
`⟨clip(F)·F⟩/⟨F²⟩ = √(2/π) = 0.79788`; saf sinüs için `√2·(2/π) = 0.90032`.
Ölçülen (Hkeskin, model alanları):

```
E bacağı: 0.82775      X bacağı: 0.82144
```

— alanlar Gauss'a (0.798) sinüsten (0.900) çok daha yakın, ama tam Gauss
değil (%3.7/%3.0 üstünde).

**Arcsine ailesi** (model alanı ↔ ÖLÇÜLEN alan çifti; Gauss çiftinde
`⟨sgn x sgn y⟩ = (2/π) arcsin r`):

| bacak | r | ⟨sgn·sgn⟩ | (2/π)arcsin r | fark | ⟨sgn sgn⟩/r |
|---|---|---|---|---|---|
| **E** (E_mod ↔ η) | 0.87704 | 0.68053 | 0.68097 | **−0.07%** | 0.77594 |
| **X** (X_mod ↔ X̃0) | 0.88023 | 0.74898 | 0.68522 | **+9.31%** | 0.85090 |

> E bacağı arcsine yasasını **‰0.7** içinde sağlıyor (çift Gauss); X bacağı
> %9.3 sapıyor (X̃0 Gauss değil). Yani "Gauss çiftin işaret-korelasyonu"
> resmi E'de TAM, X'te YAKLAŞIK.

### K2.3 İŞARET-AKTARIMI — KORELATÖR ÜZERİNDE (bant bant)

**Bir tuzak ve düzeltmesi (dürüstlük).** İlk deneme (`169_k2.py`) `G₂ = E·X·X`
içinde İKİ X'i birden kırptı; ama `clip(X)·clip(X) = σ_X²` SABİTTİR —
rezonans yok olur, o iki sütun (`ρ(X1X2)`, `ρ(E+X1X2)`) **kullanılamaz**.
`169_k2b.py` bacakları GERÇEKTEN ayırıyor: model çizgileri τ'ya göre
sıralanıp bir atlamalı iki AYRIK yarıya bölünüyor (|A| = 4491, |B| = 4490,
σ_Xa = 0.19897, σ_Xb = 0.19548, korr(Xa,Xb) = +0.177) ve sınav nesnesi

```
G = E · X_a · X_b        (Ç4 ile aynı yapı: üç BAĞIMSIZ alan bacağı + taşıyıcı)
```

Alan düzeyinde tek-bacak kırpma aktarımı `⟨clip(F)F⟩/⟨F²⟩`:
`E: 0.82775`, `X_a: 0.80854`, `X_b: 0.80915` (Gauss √(2/π) = 0.79788).

Korelatör düzeyinde `ρ = Pu2(varyant)/Pu2(tam)`, bant bant:

| bant τ_eff | ρ(E) | ρ(X_a) | ρ(X_b) | ρ(X_aX_b) | ρ(EX_a) | ρ(hepsi) |
|---|---|---|---|---|---|---|
| 0.5393 | 0.7142 | 0.8163 | 0.8348 | 0.5906 | 0.5873 | 0.4248 |
| 0.5792 | 0.7613 | 0.8180 | 0.8316 | 0.5866 | 0.6342 | 0.4561 |
| 0.6189 | 0.8042 | 0.8508 | 0.8368 | 0.6135 | 0.7077 | 0.5297 |
| 0.6581 | 0.8271 | 0.8652 | 0.8542 | 0.6471 | 0.7496 | 0.5912 |
| 0.6983 | 0.8669 | 0.8814 | 0.8768 | 0.6721 | 0.7973 | 0.6400 |
| **ORTALAMA** | 0.7947 | 0.8463 | 0.8468 | 0.6220 | 0.6952 | **0.5284** |

| kırpılan bacak sayısı | **ölçülen** | Gauss ailesi (2/π)^{n/2} | fark |
|---|---|---|---|
| 1 | **0.8293** | 0.79788 | **+3.94%** |
| 2 | **0.6586** | 0.63662 (= **2/π**) | **+3.45%** |
| 3 | **0.5284** | 0.50791 | **+4.03%** |
| *4 (uzatma, ort. kayma ×1.0381)* | *0.4207* | *0.40528 (= **4/π²**)* | *+3.81%* |

> **HÜKÜM (K2b — BÜYÜK ÖLÇÜDE GEÇTİ).** Korelatörün işaret-kırpma
> aktarımı **Gauss ailesini `(2/π)^{n/2}` bacak bacak, SABİT +%3.8'lik
> bir kaymayla izliyor** (kayma bacak sayısıyla neredeyse hiç değişmiyor:
> 1.0394 / 1.0345 / 1.0403) (alanlar tam Gauss değil, sinüse doğru %4 kayık —
> §K2.2'nin alan-düzeyi ölçümüyle birebir tutarlı). Üç bacaktan çıkan
> ETKİN bacak-başı aktarım `0.5284^{1/3} = 0.8085` = √(2/π) + %1.3.
> **İKİ bacakta ölçülen 0.6586, `2/π = 0.63662`'nin %3.5 üstünde.**
>
> **AMA:** KALEM'in "bacak başına 2/π" muhasebesi YANLIŞ. Ölçülen
> bacak-başı aktarım `√(2/π) = 0.798`'dir; `2/π` İKİ bacakta çıkar.
> Aynı sabit `4/π²` için **DÖRT** kırpılmış bacak gerekir. Bu, Ç4'ün
> yapısıyla (dört merdiven frekansı: `h_Q`, `h'_1`, `y_2`, `y_3` —
> 165 §2, "dördüncü kat") tam örtüşür; ama dördüncü bacak (taşıyıcının
> kendi η-izdüşümü `h_Q`) bu kurulumda kırpılamaz, **dolayısıyla
> 4 bacaklı kapanış ÜÇ ÖLÇÜLEN bacaktan bir adım UZATMADIR.**

**Arcsine ARA-DEĞER eğrisi (kısmi doyum).** Model alanı yalnız bandın
τ penceresindeki çizgilerden kurulunca `r` düşer (kısmi doyum); Gauss
çiftinde `⟨sgn x sgn y⟩ = (2/π) arcsin r` beklenir:

| bant lo | n_çizgi | r | ⟨sgn·sgn⟩ | (2/π)arcsin r | fark |
|---|---|---|---|---|---|
| 0.52 | 53 | +0.39118 | +0.26480 | +0.25587 | +3.49% |
| 0.56 | 75 | +0.35468 | +0.23740 | +0.23082 | +2.85% |
| 0.60 | 116 | +0.32068 | +0.21045 | +0.20782 | +1.26% |
| 0.64 | 176 | +0.28116 | +0.18617 | +0.18144 | +2.61% |
| 0.68 | 264 | +0.23418 | +0.15167 | +0.15048 | +0.79% |

> Arcsine ailesi **beş bantta da %0.8–3.5 içinde** tutuyor (aynı +%3
> mertebesindeki Gauss'tan sapma). Kısmi-doyum ara-değerleri isteniyordu:
> `r = 0.234 … 0.391` aralığında eğri ölçülmüş ve doğrulanmıştır.

---

## K3 — YENİ ÖNGÖRÜ

### K3.1 (a) β = 0.2175'in TÜRETME DENEMESİ — **DÜŞTÜ** (açık yazılıyor)

168'in açık borcu: `θ/θ₀ = 1 − βφ` yasasının katsayısı. Kırpılmış modelden
üç deneme yapıldı; **üçü de düştü.**

**Deneme 1 — KUVVET (katılımcı sayımı).** Ç4 dört merdiven frekansı
istiyor; her bacağın "gerçek" olma olasılığı `1−φ` ise `θ/θ₀ = (1−φ)^k`
olmalı, `k` = gereken gerçek bacak sayısı. Ölçülen `k`:

| gaz | φ | θ/θ₀ | **k (kuvvet modeli)** | β (doğrusal) |
|---|---|---|---|---|
| K090 | 0.4153 | 0.9080 | **0.1798** | 0.2215 |
| K070 | 0.9270 | 0.8038 | **0.0834** | 0.2117 |
| HA4 | 0.9098 | 0.7975 | **0.0941** | 0.2226 |
| E060 | 0.9595 | 0.7945 | **0.0717** | 0.2142 |

> `k` gazdan gaza 2.5 kat değişiyor (0.072–0.180) — **kuvvet/olasılık
> sayımı REDDEDİLDİ.** `β` ise %2.5 içinde sabit: yasa gerçekten
> DOĞRUSAL, üstel değil. (Bu, 168'in şeklini bağımsız olarak doğruluyor.)

**Deneme 2 — SEYRELME (varyans muhasebesi).** Kesim gazında model hâlâ
8981 çizgi uydurur; φ kadarı gazın merdiveninde YOKTUR ("boş çizgi") ve
gürültüye uyar. `E = E_g + E_b`, `X = X_g + X_b` (gerçek + boş) yazılıp
boş parçaların sinyalle ilişkisiz olduğu varsayılırsa

```
g_cal = g_E g_X² ∝ (1−ρ_E)(1−ρ_X)² ,   ρ = boş çizgilerin varyans payı
KALİB'in koherent parçası değişmez        ⇒  θ = KALİB/g_cal ∝ 1/[(1−ρ_E)(1−ρ_X)²]
```

yani θ φ ile **BÜYÜMELİ**. Ölçülen θ **KÜÇÜLÜYOR** ⇒ **saf seyrelme
İŞARETİYLE reddedildi.**

**Deneme 3 — HİPOTEZ B (boş çizgilerin aşırı-uyumu), DOĞRUDAN SINANDI.**
Eğer φ-bağımlılığı boş çizgilerden geliyorsa, **model alanından boş
çizgileri çıkarınca φ-bağımlılığı KAYBOLMALIDIR.** `169_configs/169_k3_bos.py`
bunu ölçüyor (`E_g`, `X_g` yalnız GERÇEK çizgilerden sentezleniyor;
regresyon ve öngörü baştan yeniden):

| gaz | φ | boş çizginin Var payı (E / X) | Π = Pu2^g/Pu2 | g_cal → g_cal^g | θ/θ₀ | **θ^g/θ₀** |
|---|---|---|---|---|---|---|
| K070 | 0.9270 | %2.9 / %16.3 | 0.885 | 0.3128 → 0.3862 | 0.8039 | **0.7322** |
| K090 | 0.4153 | %0.9 / %2.7 | 0.978 | 0.2877 → 0.2966 | 0.9080 | **0.8991** |

(Her iki gazda `θ/θ₀` 168 §A2.7'nin sayısını dört hanede yeniden üretti:
0.8039 ↔ 0.8038 ve 0.9080 ↔ 0.9080.)

> **HÜKÜM (K3a DÜŞTÜ).** Boş çizgileri modelden ATINCA sapma
> KAPANMIYOR, **BÜYÜYOR** (K070: −%19.6 → −%26.8; K090: −%9.2 → −%10.1).
> Yani `β` bir
> "boş-çizgi aşırı-uyumu" artefaktı DEĞİL; kesim yasası gazın GERÇEK
> merdiveninin bir özelliğidir (az sayıda gerçek çizgi ⇒ zayıf Ç4) —
> 168'in "katkı ŞİDDETTEN değil SAYIDAN geliyor" hükmüyle aynı yönde,
> ama **katsayı hâlâ TÜRETİLMEMİŞTİR.**
>
> **Sayısal tesadüf (KAZANILMAMIŞ, sonuç değildir).** Ölçülen
> `β = 0.2175 ± 0.0054` (168, dört gaz) ve `β = 0.2143` (169, K1.3'ün
> mutlak uyumu) değerleri `1 − π/4 = 0.21460`'a 0.5σ, `π²/45 = 0.21932`'ye
> 0.8σ uzaklıkta. **Hiçbirinin türetimi yoktur; buraya yalnız kayıt
> için yazılmıştır.**

### K3.2 (b) ÖRNEKLEM-DIŞI ÖNGÖRÜ — λ = 0.70 GAZI (ÖN KAYIT)

168 §A3, λ ekseninde tek nokta (L085) ile şunu ölçmüştü: **ham kalibrasyon
`KALİB_u2` λ-DEĞİŞMEZ** (L085/Hkeskin = 1.014), ama `W_X` λ ile büyüdüğü
için `c = KALİB/W_X` λ ile DÜŞÜYOR (0.946). H-C1 ise `c`'nin evrensel bir
geometrik sabit olduğunu söyler. **İkisi aynı anda doğru olamaz.** λ = 0.70
gazı bu ikisini ayırır.

`z_L070.npy` 169'da İNŞA EDİLDİ (`167_insa.py L070`, 9.6 dk, ilk-kök
hücreleri benzersiz 300000/300000, maks|F| < 2e−9, sıralılık TAM,
σ_ds = 0.35586). Gazın X̃ marjinalinden `W_X` — **korelatörden bağımsız,
ölçümden ÖNCE** (`169_k3_onkayit.py`, `scratchpad/169/ONKAYIT_L070.json`,
yazılma saati **14:14:14**; `167_olcum.py L070` **14:14:38**'de başlatıldı):

σ_X̃(L070) = 0.19897 (Hkeskin: 0.24204).

| bant lo | W_X(L070) | W_X(Hkeskin) | c(Hkeskin) | **P_A** | **P_B** |
|---|---|---|---|---|---|
| 0.52 | 0.7959 | 0.7128 | 0.4016 | 0.3597 | 0.4016 |
| 0.56 | 0.7685 | 0.6766 | 0.3999 | 0.3521 | 0.3999 |
| 0.60 | 0.7403 | 0.6399 | 0.4041 | 0.3493 | 0.4041 |
| 0.64 | 0.7110 | 0.6025 | 0.4049 | 0.3430 | 0.4049 |
| 0.68 | 0.6814 | 0.5653 | 0.4073 | 0.3379 | 0.4073 |

```
** ÖN KAYIT (ölçümden ÖNCE yazıldı) **
   P_A  (168 §A3: KALİB_u2 λ-değişmez)  :  c(L070) = 0.3483
   P_B  (H-C1: c evrensel = 4/π²)       :  c(L070) = 0.4035   [4/π² = 0.4053]
   ayrım +15.9%  ≈ 18σ
```

**ÖLÇÜM** (`167_olcum.py L070 0.40 0.95 0 0`, 5.0 dk; ölçüm zinciri
167'nin kendisi, kopyalanmış kod yok):

| τ_eff | KALİB_u2 ± jk | W_X | **c_WX ± jk** |
|---|---|---|---|
| 0.5393 | 0.3056 ± 0.0015 | 0.7959 | 0.3839 ± 0.0020 |
| 0.5792 | 0.2884 ± 0.0018 | 0.7685 | 0.3753 ± 0.0024 |
| 0.6189 | 0.2759 ± 0.0012 | 0.7403 | 0.3727 ± 0.0016 |
| 0.6579 | 0.2582 ± 0.0021 | 0.7110 | 0.3631 ± 0.0028 |
| 0.6981 | 0.2392 ± 0.0020 | 0.6814 | 0.3510 ± 0.0029 |

```
c(L070) = 0.3690 ± 0.0011(jk) ± 0.0057(bant) = ±0.0058     [5 bant]

   P_A öngörüsü 0.3483  →  ölçüm  +5.9%  ( +3.6σ_tot )     P_A DÜŞTÜ
   P_B öngörüsü 0.4035  →  ölçüm  −8.6%  ( −5.9σ_tot )     P_B DÜŞTÜ
   4/π² = 0.4053        →  ölçüm  −8.95% ( −6.29σ_tot )
```

> **HÜKÜM (K3b DÜŞTÜ — ama iki yönde birden).** Örneklem-dışı ölçüm İKİ
> rakip öngörünün de ARASINA düştü. (i) **H-C1'in evrensellik öngörüsü
> −5.9σ ile reddedildi**: `c`, λ = 1.00 → 0.85 → 0.70 boyunca
> **0.4035 → 0.3817 → 0.3690** diye tekdüze iniyor; bu bir gürültü değil,
> ölçülen bir eksen. (ii) **168 §A3'ün "KALİB_u2 λ-değişmez" yasası da
> düştü**: `KALİB_u2` (bant 1) λ ile 0.2863 → 0.2920 (+%2.0) → 0.3056
> (+%6.7) diye BÜYÜYOR; 168 bu yasayı yalnız λ=0.85'te (+%1.4) ölçmüştü
> ve orada gerçekten küçüktü — λ=0.70 onu kırıyor. **λ ekseni ne saf DW
> (W_X) ne de saf doyum resmiyle açıklanıyor.**
>
> Not (dürüstlük): P_B'nin öngörüsü `g_cal` düzeltmesi İÇERMİYORDU;
> düzeltme dahil edilirse (g₀/g_cal = 0.9694) beklenen `c*` = 0.3577,
> yani 4/π²'den −%11.7 — **düzeltme sapmayı BÜYÜTÜYOR**, küçültmüyor.

---

## 4. HÜKÜM

| kapı | sonuç | dayanak |
|---|---|---|
| **K1 — HASSASİYET** | **KOŞULLU GEÇTİ** | taban gazları 4/π²'ye ‰0.4 / ‰0.0 (dokuz bant: 0.4051 / 0.4053); ortak pencerede −0.97σ / +3.06σ; **kesim yasasının mutlak φ→0 kesişimi 0.4057 ± 0.0016 = 4/π² + 0.23σ**. **AMA** λ ekseni −6.3σ ve hiçbir üye iki seçim ölçütünü birden sağlamıyor (§K1.5). |
| **K2 — MEKANİZMA** | **GEÇTİ** | (a) Ç4 rejiminde faz mod 2π TAM DÜZGÜN (Kuiper √N·V = 1.10–1.68, eşik 2.00) — aynı ölçüt tek çizgilerde 78–500 veriyor: **beş mertebe kontrast**. (b) korelatörün kırpma aktarımı Gauss ailesini `(2/π)^{n/2}` sabit +%3.8 kaymayla izliyor; 2 bacak = 0.6586 ↔ 2/π; arcsine ara-değer eğrisi beş bantta %0.8–3.5. |
| **K3 — YENİ ÖNGÖRÜ** | **DÜŞTÜ** | (a) β'nın üç türetme denemesi de düştü; Hipotez B (boş-çizgi aşırı-uyumu) DOĞRUDAN ölçülüp ÇÜRÜTÜLDÜ. (b) ÖN KAYITLI örneklem-dışı öngörü λ=0.70'te **−5.9σ ile reddedildi** (rakip 168-yasası da +3.6σ ile reddedildi). |

### 0'ın cevabı — TEK CÜMLELİK HÜKÜM

> **KIRPILMIŞ-FAZ MEKANİZMASI DOĞRULANDI, SABİTİN EVRENSELLİĞİ ÇÜRÜTÜLDÜ:**
> Ç4 rejiminde birikmiş faz mod 2π **ölçülerek** tam sarılıyor (Kuiper
> √N·V = 1.1–1.7, tek çizgilerde 78–500) ve korelatörün işaret-kırpma
> aktarımı Gauss `(2/π)^{n/2}` ailesini bacak bacak %3.8 içinde izliyor —
> yani `c`'nin `4/π²` **ailesinden** olduğu artık ölçülmüş bir olgudur;
> ama `c`'nin KENDİSİ evrensel bir sayı DEĞİL: aynı merdivenin yalnız
> genlik ölçeği λ = 1.00 → 0.85 → 0.70 yapılınca `c` = 0.4035 → 0.3817 →
> **0.3690** diye tekdüze iniyor ve ÖN KAYITLI öngörü **−5.9σ**'da
> düşüyor. **H-C1 MÜHÜRSÜZ: "güçlü aday", ama λ ekseninde çürük.**

**MÜHÜR KURALININ UYGULANMASI (açık).** KALEM'in kuralı "K1 + (K2 veya
K3)" idi; harfiyen K1 (koşullu) + K2 geçiyor. **Mühür yine de VERİLMİYOR**,
çünkü K1'in geçtiği eksen (kesim) ile K3'ün düştüğü eksen (λ) AYNI
büyüklüğün iki ölçümüdür ve K3 ÖN KAYITLI, örneklem-dışı ve
6σ'lıktır. Bir sabitin "evrensel" ilan edilmesi, evrenselliğin
sınandığı ve düştüğü bir eksen varken dürüst değildir.

### Ne kazanıldı, ne kaldı

**KAZANILDI (169'un net çıktısı)**
1. **Sarılma bir varsayım değil, ölçüm** (§K2.1) — 168 §A1.5'in
   `2πτσ_Ĉ ≫ 1` uyarısı artık dağılım düzeyinde doğrulandı.
2. **Kırpma aktarımı Gauss ailesini izliyor** (§K2.3) ve KALEM'in
   muhasebesi düzeltildi: **bacak başına √(2/π), 2/π İKİ bacakta,
   4/π² DÖRT bacakta** — Ç4'ün dört merdiven frekansıyla birebir.
3. **Kesim yasasının MUTLAK kesişimi** `c₀ = 0.4057 ± 0.0016` (§K1.3):
   sıfır serbest parametreli 4/π² hedefini ‰1 içinde veriyor; aynı
   uyumun eğimi β = 0.2143, 168'in bağımsız β = 0.2175'iyle 0.6σ.
4. **λ ekseni ikinci noktasını kazandı** (168'in 4. açık borcu): L070
   inşa edildi ve ölçüldü; hem H-C1 hem 168 §A3'ün λ-yasası düştü.
5. **β'nın adresi daraldı**: boş-çizgi artefaktı DEĞİL (§K3.1).

**KALDI (169'dan sonraki borçlar)**
1. **λ ekseni.** `c(λ)` neden düşüyor? Ne `W_X` (DW) ne de doyum resmi
   veriyor. Üye sınavı (§K1.5) `√W_X`'i işaret ediyor ama o üye bant
   eğimini ve mutlak değeri bozuyor. **Bu, 169'un bıraktığı ASIL borç.**
2. **Dördüncü bacağın kırpılması** — taşıyıcının kendi η-izdüşümü `h_Q`;
   4/π² kapanışı hâlâ bir adım uzatma (§K2.3).
3. **β = 0.2175'in türetimi** (168'den devraldı, kapanmadı).
4. **T/4 dilimlerindeki +%7…%10** (168 §6.2, kapanmadı).
5. λ = 1.15 (L115) inşası yarım kaldı (167'den beri).

---

## 5. DENETİM ve DÜRÜSTLÜK NOTLARI

* **168'in tablosu yeniden üretildi** — bağımsız kod yolundan (`169_k1.py`),
  yedi gazda `c_ampX`, `g_cal`, `θ` ve (dokuz bantta) `c_WX` **dört
  hanede birebir** (§K1.1). Ayrıca `169_k3_bos.py` bağımsız olarak
  `θ/θ₀`(K070) = 0.8039 ↔ 168'in 0.8038 ve `θ/θ₀`(K090) = 0.9080 ↔ 0.9080
  verdi.
* **Ölçüm parçası kopyalanmadı**: `166_T1.bant_agg`/`_jk`,
  `165_cekirdek.Model165/tayf_s/sentez/s_den_J`, `163_cekirdek.olc_cizgi/
  bant_adaylari`, `166_bacak.karakteristik`, `167_insa/167_olcum` aynen
  import edildi.
* **Bir hata yapıldı ve düzeltildi:** `169_k2.py`'nin `ρ(X1X2)` ve
  `ρ(E+X1X2)` sütunları GEÇERSİZDİR — aynı alanı iki kez kırpmak
  `clip(X)² = σ_X²` sabitini verir ve rezonansı yok eder. Doğru ölçüm
  `169_k2b.py`'nin BÖLÜNMÜŞ MERDİVEN kurulumudur (§K2.3); rapordaki
  bütün çok-bacak sayıları oradan gelir. (`169_k2.py`'nin sarılma tablosu
  ve arcsine satırları geçerlidir.) `169_k2.py` ayrıca sonunda bir JSON
  serileştirme hatasıyla bitiyor (numpy bool); ekran çıktısı tamdır ve
  `scratchpad/169/log_K2_Hkeskin.txt`'de saklıdır.
* **ÖN KAYIT sırası kanıtlı:** `ONKAYIT_L070.json` **14:14:14**'te,
  `167_olcum.py L070` **14:14:38**'de başlatıldı; öngörüler yalnız gazın
  X̃ MARJİNALİNDEN (`W_X`) ve Hkeskin'in ölçülmüş `c`'sinden kuruldu,
  L070'in korelatörüne hiç bakılmadı.
* **L070 gazı 164'ün sadakatli çözücüsüyle kuruldu** (`167_insa.py L070`,
  9.6 dk): ilk-kök hücreleri benzersiz 300000/300000, maks|F| = 1.9e−9
  (eşik 1e−8), sıralılık TAM, min Δz = 0.13, σ_ds = 0.35586.
* **Kuiper eşikleri**: düzgün dağılımda `√N·V`'nin asimptotik medyanı
  ≈1.62, %5 = 1.747, %1 = 2.001. Ç4 satırları 1.10–1.68 ⇒ düzgünlük
  REDDEDİLEMİYOR; tek-çizgi satırları 78–500 ⇒ ezici RED. `R`'nin düzgün
  beklentisi `√π/2/√N = 0.00161`.
* **σ tanımı**: `σ_jk` bant-içi 8-grup jackknife'ların gaz düzeyine
  taşınması, `σ_bant` bantlar arası saçılımın standart hatası,
  `σ_tot = hipot(σ_jk, σ_bant)`. Raporda hüküm hep `σ_tot` ile verildi
  (daha muhafazakâr).
* **`g_cal` düzeltmesine bağımlılık açıktır**: K1.3'ün kesişimi `g_cal`
  düzeltmesi olmadan 0.4002 ± 0.0016 (−3.2σ), düzeltmeyle 0.4057 ± 0.0016
  (+0.23σ). Düzeltme 168 §A2.6'nın ÖLÇÜLEN ayrıştırmasıdır, 169'un
  uydurması değil — ama sonuç ona bağlıdır.
* **Bant penceresi**: bütün eksen karşılaştırmaları ortak pencerede
  (lo ∈ [0.52, 0.68], SNR ≥ 3, R_bant ≥ 0.98, τ_eff < 0.85 — beş bant);
  168 §A2.8'in dokuz-bant sayıları da ayrıca verildi (§K1.1/K1.5) ve
  ikisi arasındaki fark (Hkeskin: 0.4035 ↔ 0.4051) hiçbir yerde
  gizlenmedi.
* **Uydurma yok**: bu raporda hesaplanmamış tek sayı yoktur; her tablo
  `scratchpad/169`'daki bir log/JSON'a karşılık gelir
  (`K1.json`, `K2b_Hkeskin.json`, `UYE.json`, `ONKAYIT_L070.json`,
  `BOS_K070.json`, `BOS_K090.json`, `log_*.txt`).
