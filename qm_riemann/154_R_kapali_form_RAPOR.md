# 154 — R(τ) tayfının kapalı-form avı

**Hedef:** Γ_rot(τ) = ⟨w·e^{−iA·dsΔ}⟩/⟨w⟩, w = e^{−A·R(τ)·dsΔ} mekanizmasının kalan tek türetme dişlisi R(τ) için kapalı form. Görev dört adım istedi: (1) veriyi ikinci pencereyle güçlendir, (2) sentetik tanık, (3) az-parametreli yarış, (4) hüküm.

> ## Kısa hüküm
>
> 1. **R(τ) gerçek, L-değişmez bir fonksiyondur.** İki bağımsız pencere (ΔL = 0.566) dört sağlam bantta %0.5–2.0 içinde aynı R'yi veriyor. 150'nin tayfı birebir yeniden üretildi.
> 2. **Sıfır-parametreli aday (d) — R = −dlnρ_erfc/dA — ÖLDÜ.** τ=0.5375'te tesadüfen tutuyor (0.414 vs 0.391) ama τ=0.74'te 4.59 öngörüp 2.12 buluyor: %117 sapma.
> 3. **Ölçülmüş ρ ile kurulan sıfır-parametreli sürüm (f) de ÖLDÜ.** R/(−dlnρ/dA) oranı bantlar boyunca 0.57 → 2.52 → 1.23 geziniyor; sabit değil.
> 4. **Sentetik tanık bu aileyi mekanizma düzeyinde kesiyor:** keskin merdiven gazında ρ(τ) 0.53–0.74'te PRATİK OLARAK DÜZ (0.415→0.392) olduğu hâlde R yine 1.57 → 4.57 tırmanıyor. R'nin tırmanışı soğurmadan gelmiyor. **Şekil "makine"dendir**; ama sayılar (τ₀ ve ölçek) gerçeğe özgüdür.
> 5. **"Taban-kenarı" motifi (aday c'nin gerekçesi) ÖLDÜ:** taban 0.46/0.52/0.58 taramasında τ=0.60'ta R = 1.305/1.351/1.462 — motif doğru olsaydı taban 0.58 sütunu R≈0.4 vermeliydi. R, (τ−taban)'ın fonksiyonu değil. (Aday c ampirik ŞEKİL olarak hayatta; gerekçesi değil.)
> 6. **Hayatta kalanlar (ikisi de 2 parametre, ayırt edilemez):** (c) R = R∞(1−e^{−(τ−0.52)/w}) ve (g) R = c(τ−τ₀)/τ². Ölçülmüş sistematiklerle χ²/dof = 0.34 ve 0.33.
> 7. **En sağlam yeni sayı:** R(τ) sıfırını τ₀ = **0.5153 ± 0.0008(ist.) ± 0.002(pencere/konvansiyon)**'de kesiyor — DOĞRUDAN ölçüldü (taban 0.46 koşusunda τ=0.48/0.51 bantlarında R = −0.706 / −0.099), tabandan ve L'den bağımsız. **½ değildir** (muhafazakâr hatayla bile ~7σ). (g) adayı bu sayıyı yalnız τ≥0.5375 bantlarından fit ederek 0.513–0.519 buluyor.
> 8. **Kapalı form MÜHÜRLENMEDİ.** İki 2-parametreli form ayırt edilemiyor, ve sıfır-parametreli hiçbir türetim ayakta değil. AÇIK.

---

## 1. Yöntem — tek kod yolu, üç denetim

`154_configs/154_cekirdek.py`: 150'nin zinciri (ds → 0.52-taban regresyonu → η → bant Γ_rot → R faz-eşlemesi) tek fonksiyona (`olc`) kapatıldı. **Gerçek pencereler, taban taraması ve sentetik gazlar aynı fonksiyonu çağırır** — 152 raporunun uyardığı kopya-kayması riski sıfır.

150'ye göre kasıtlı eklemeler (ölçümü değiştirmez, denetim/hata üretir):

| | ekleme |
|---|---|
| E1 | R taraması hem 150'nin `linspace(0,3,121)` ızgarasıyla (sütun `R_150`, tekrar-üretim kontrolü) hem tavanı 8'e çekilmiş ızgara + bisection ile (`R_ince`) yapılır |
| E2 | faz artığı \|arg M(R) − arg Γ\| kaydedilir; >0.02 rad ise R eşleşmesi TUTMAMIŞTIR |
| E3 | jackknife: aday çizgiler 8 gruba bölünüp birer grup dışarıda bırakılarak R yeniden ölçülür → bant-içi istatistiksel hata |
| E4 | n_eff = (Σw)²/Σw² — ağırlık birkaç uç noktaya çökmüşse ortalama anlamsızdır |
| E5 | w = exp(x − max x) (oran alındığından sabit sadeleşir; taşma yok) |
| E6 | **aynı aday çizgilerden ρ(τ) da ölçülür** (144'ün tanımı birebir: ρ = (Σon−Σoff)/Σ(2a·sin(πω/L))², a = 1/(πm√q)). Böylece R ve ρ aynı bantlardan, aynı koşudan çıkar |

**Denetim 1 — tekrar-üretim.** Son-300k penceresinde `R_150` sütunu: **0.40 / 1.07 / 1.85 / 2.12 / 1.50** — 150'nin yayımlanmış tayfıyla birebir. σ_ds² = 0.1674, σ_η² = 0.0227, c₁ = −0.01158: birebir.

**Denetim 2 — sentetik kopya kayması.** 152_gaz.py'nin merdiven+Newton bloğu 154_sentetik.py'ye kopyalandı; koşu sonunda 152'nin S3 üçlüsüyle karşılaştırıldı:

| | σ_ds² | σ_η² | c₁ |
|---|---|---|---|
| keskin — bu koşu / 152 | 0.1280 / **0.1280** | 0.0767 / **0.0767** | −0.03535 / **−0.03535** |
| A4 — bu koşu / 152 | 0.1128 / **0.1128** | 0.0230 / **0.0230** | −0.00918 / **−0.00918** |

Dört anlamlı basamakta özdeş — kopya sapmamıştır.

**Denetim 3 — kök seçimi.** Negatif R dalı açıldığında (φ<0 bantları için gerekli) arg M(R)'nin büyük \|R\|'de monoton OLMADIĞI, sahte kök ürettiği görüldü (τ=0.69 bandı bir ara −2.87'ye yapıştı, artık 1.56). Düzeltildi: kök R=0'dan **dışa doğru** aranır, ilk (en küçük \|R\|) kök alınır. Rmin=0 ile davranış 150'ninkiyle özdeş; düzeltmeden sonra bütün pozitif-R bantları birebir aynı çıktı (kontrol koşusu: 1.351 / 2.135 → 1.351 / 2.135).

---

## 2. İki pencere — L-değişmezlik ve hata

| pencere | t aralığı | n | L | σΔ² | σ_ds² | σ_η² | c₁ |
|---|---|---|---|---|---|---|---|
| son | 975 796 – 1 132 491 | 300 000 | 12.0296 | 0.0547 | 0.1674 | 0.0227 | −0.01158 |
| orta | 517 561 – 681 994 | 300 000 | 11.4638 | 0.0543 | 0.1666 | 0.0225 | −0.01154 |

Standart bantlar (150 ile aynı):

| τ̄ | R_son | ±jk | R_orta | ±jk | \|fark\| | % | R_gauss son | R_gauss orta | ρ_son | ρ_orta | n_eff |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.5375 | 0.391 | 0.013 | 0.399 | 0.019 | 0.008 | 2.0 | 0.410 | 0.418 | 0.4217 | 0.4237 | 274 950 |
| 0.5850 | 1.086 | 0.015 | 1.081 | 0.017 | 0.005 | 0.5 | 0.982 | 0.976 | 0.3463 | 0.3440 | 153 942 |
| 0.6600 | 1.852 | 0.019 | 1.882 | 0.029 | 0.030 | 1.6 | 1.469 | 1.481 | 0.2464 | 0.2460 | 41 855 |
| 0.7400 | 2.120 | 0.028 | 2.144 | 0.019 | 0.024 | 1.1 | 1.574 | 1.579 | 0.1577 | 0.1564 | 17 361 |
| 0.8150 | 1.507 | 0.049 | 1.462 | 0.044 | 0.046 | 3.1 | 1.205 | 1.173 | 0.0787 | 0.0789 | 41 212 |

**Hüküm:** ΔL = 0.566 (çizgi kümeleri tümüyle farklı) olmasına rağmen dört sağlam bantta uyuşma %0.5–2.0. **R(τ) gerçek bir τ fonksiyonudur; L'ye asılı değildir** — 138'in L-değişmezliğiyle tutarlı. τ=0.815 bandı iki pencerede %3.1 ayrışıyor; görev notunun uyarısı doğrulanıyor, fite katılmadı.

**İki tahminci ayrışıyor.** `R_gauss = φ/(A²σΔ²)` (Gauss limiti) ile `R_emp` (150'nin ampirik-ağırlıklı tanımı) düşük τ'da uyuşuyor (0.410 vs 0.391) ama τ=0.74'te %35 ayrışıyor (1.574 vs 2.120). Ayrışma Gauss-dışı düzeltmenin büyüklüğüdür. **İlkel gözlenebilir ikisi de değil, φ_Γ'dır;** R hangi tanımla okunduğuna bağlıdır ve bu bir kapalı-form iddiasında belirtilmek zorundadır. Bu raporda birincil nesne, mekanizmanın kendi tanımı olan **R_emp**'tir.

---

## 3. İnce ızgara — R(τ) ve ρ(τ) aynı koşudan

0.03 genişlikli 9 bant, iki pencere:

| τ̄ | R_son | ±jk | R_orta | ±jk | ρ_son | ± | ρ_orta | ± | φ_son |
|---|---|---|---|---|---|---|---|---|---|
| 0.540 | 0.423 | 0.012 | 0.438 | 0.011 | 0.4170 | 0.0019 | 0.4155 | 0.0025 | +0.276 |
| 0.570 | 0.921 | 0.005 | 0.910 | 0.006 | 0.3622 | 0.0015 | 0.3618 | 0.0015 | +0.601 |
| 0.600 | 1.351 | 0.007 | 1.340 | 0.008 | 0.3254 | 0.0012 | 0.3247 | 0.0013 | +0.910 |
| 0.630 | 1.678 | 0.010 | 1.720 | 0.007 | 0.2838 | 0.0026 | 0.2855 | 0.0012 | +1.179 |
| 0.660 | 1.945 | 0.009 | 2.013 | 0.014 | 0.2449 | 0.0019 | 0.2374 | 0.0022 | +1.432 |
| 0.690 | 2.135 | 0.010 | 2.192 | 0.012 | 0.2068 | 0.0020 | 0.2061 | 0.0023 | +1.655 |
| 0.720 | 2.213 | 0.020 | 2.194 | 0.007 | 0.1720 | 0.0030 | 0.1727 | 0.0016 | +1.827 |
| 0.750 | 2.166 | 0.031 | 2.175 | 0.023 | 0.1488 | 0.0023 | 0.1444 | 0.0017 | +1.931 |
| 0.780 | 1.799 | 0.048 | 1.867 | 0.025 | 0.1079 | 0.0037 | 0.1095 | 0.0026 | +1.923 |

Şekil: R hızla tırmanıyor, τ≈0.72'de ~2.2'de düzleşiyor, 0.78'de düşüyor (çözünmeyen-çizgi bölgesinin başlangıcı — fitte kullanılmadı).

Not: geniş bantla ölçülen R, dar bantla ölçülenden farklıdır (τ=0.66'da geniş bant 1.852, dar bant 1.945) — R bandın içinde değiştiği için bant-ortalaması bant genişliğine bağlıdır. Aynı ızgara her karşılaştırmada kullanıldı.

### Aday (f): R = −(1/π)·dlnρ/dτ = 2·(−dlnρ/dA)

2 çarpanının gerekçesi ciddi: bond İKİ adımlıdır (dsΔ = (ds_n+ds_{n+1})/2), adım başına e^{−A·R₁·ds} iki kez çarpılırsa ölçülen R = 2R₁ olur. lnρ'ya kübik fitin analitik türeviyle:

| τ | R_ölç(son) | −dlnρ/dA | ×2 | **oran R/(−dlnρ/dA)** |
|---|---|---|---|---|
| 0.540 | 0.423 | 0.745 | 1.490 | **0.57** |
| 0.570 | 0.921 | 0.689 | 1.379 | **1.34** |
| 0.600 | 1.351 | 0.675 | 1.351 | **2.00** |
| 0.630 | 1.678 | 0.703 | 1.406 | **2.39** |
| 0.660 | 1.945 | 0.773 | 1.545 | **2.52** |
| 0.690 | 2.135 | 0.884 | 1.768 | **2.42** |
| 0.720 | 2.213 | 1.037 | 2.074 | **2.13** |
| 0.750 | 2.166 | 1.231 | 2.463 | **1.76** |
| 0.780 | 1.799 | 1.468 | 2.935 | **1.23** |

Oran 0.57'den 2.52'ye çıkıp 1.23'e iniyor. Orta bantlarda 2'ye yakın olması cezbedici ama **sabit değil** — iki-adım yorumu bu veriyle taşınmıyor. Sebebi de açık: ρ'nun log-eğimi 0.54–0.63'te neredeyse SABİT (−4.3), R ise aynı aralıkta 4 kat artıyor.

---

## 4. Sentetik tanık — şekil "makine"den mi?

152'nin iki konfigürasyonunda AYNI R-çıkarımı (`keskin`: 152'nin kayıtlı Newton çözümü; `A4`: erfc 0.68/0.125, Newton bu koşuda yeniden çözüldü, 21 dk).

| τ̄ | R_ger | R_kes | R_A4 | φ_ger | φ_kes | φ_A4 | φ/A² ger | φ/A² kes | φ/A² A4 | ρ_ger | ρ_kes | ρ_A4 | n_eff kes | n_eff A4 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.5375 | 0.391 | 1.574 | 1.215 | 0.256 | 0.539 | 0.595 | 0.0224 | 0.0473 | 0.0521 | 0.4217 | 0.4146 | 0.3969 | 137 461 | 157 053 |
| 0.5850 | 1.086 | 2.721 | 2.034 | 0.725 | 1.042 | 1.091 | 0.0537 | 0.0771 | 0.0808 | 0.3463 | 0.4248 | 0.2700 | 18 590 | 59 808 |
| 0.6600 | 1.852 | 3.626 | 3.867 | 1.381 | 1.822 | 1.821 | 0.0803 | 0.1059 | 0.1059 | 0.2464 | 0.4132 | 0.0534 | 1 824 | 9 029 |
| 0.7400 | 2.120 | 4.574 | ‡ | 1.860 | 2.636 | ‡ | 0.0860 | 0.1220 | ‡ | 0.1577 | 0.3924 | 0.0020 | 362 | 1 786 |
| 0.8150 | 1.507 | ‡ | ‡ | 1.727 | ‡ | ‡ | 0.0659 | ‡ | ‡ | 0.0787 | 0.3419 | 0.0080 | 183 | 1 974 |

‡ = faz sarmış (φ monoton artarken aniden −2.9/−2.7/−2.5'e düşmüş) ya da \|Γ\|>1.5 (ayrışım paydası patlamış) — ölçüm değildir. σΔ²: gerçek 0.0547, keskin 0.0350, A4 0.0429; bu yüzden ölçekten arınmış karşılaştırma φ/A² sütunlarındadır. **n_eff uyarısı:** keskin gazın τ=0.74 bandında ağırlığın etkin örneklemi 362/300 000'e düşüyor; o satır yön gösterir, ölçü değil.

**Bulgu 1 — R'nin tırmanışı soğurmadan gelmiyor (KESİN).** `keskin` gazında ρ(τ) 0.53–0.74 aralığında pratik olarak DÜZDÜR (0.415 → 0.392, %5); merdivenin tek kesimi τ=1.00'de ve keskindir, ölçülen bölgede soğurma gradyanı yoktur. Buna rağmen R 1.574 → 4.574 tırmanıyor (2.9 kat) — n_eff'i güvenilir olan τ=0.66'ya kadar bakılsa bile 1.574 → 3.626 (2.3 kat). Aynı τ aralığında A4'ün ρ'su 200 kat çöküyor (0.397 → 0.002) ve R yine tırmanıyor (1.215 → 3.867). Üç gazın ρ profilleri (düz / orta / çöken) taban tabana zıtken R hepsinde benzer oranda tırmanıyor. **−dlnρ/dA ailesi (d ve f) fit ile değil, mekanizma ile reddedilmiştir.**

**Bulgu 2 — şekil makine, sayılar değil.** "Tırmanan + doygunlaşan R, neredeyse doğrusal φ" örüntüsü soğurmasız saf asal-merdiven gazında da çıkıyor; yani örüntü Γ_rot tahmincisinin kendi geometrisinden doğuyor. Ama sayılar farklı: φ'yi doğrusal ekstrapole edince keskin gazın sıfırı τ₀ ≈ 0.487, gerçeğinki 0.516; φ/A² ölçeği keskinde 2.1 kat büyük. **Kapalı form makine analizinden aranmalı, soğurma teorisinden değil** — ama makinenin bu gaza (gerçek sıfırların ikili korelasyonuna, GUE itmesine) nasıl bağlandığı ayrı bir dişli.

**Bulgu 3 — A4'ün ρ'su gerçeğe uymuyor.** 152, A4'ün σ_η² ve c₁'i tam oturttuğunu bulmuştu; ama ρ(τ) bakımından A4 gerçekten çok uzak (τ=0.74'te 0.0020 vs 0.1577). Yani "S3'ü oturtan kesim", soğurma profilini oturtan kesim DEĞİLDİR. 152'nin A4 seçimi tek bir gözlenebilir çiftine dayanıyordu; ρ üçüncü bir kısıt olarak onu reddediyor. (Yeni kayıt.)

---

## 5. Taban taraması — R, (τ−taban)'ın fonksiyonu mu?

Aday (c)'nin gerekçesi ve ilk bandın düşüklüğüne dair "taban-kenarı seyrelmesi" şüphesi doğrudan sınandı: regresyon tabanı 0.46 / 0.52 / 0.58, üçünde de ortak 0.03'lük bant ızgarası (cap 720 → 4000; aksi hâlde taban>0.55'te sessizce devreye girip "taban değişti" iddiasını sahtelerdi). Taban 0.52 sütunu standart koşuyla birebir aynı çıktı — cap değişikliği zinciri bozmuyor.

| τ | R@0.46 | R@0.52 | R@0.58 | ρ@0.46 | ρ@0.52 | ρ@0.58 | R@0.58/R@0.46 |
|---|---|---|---|---|---|---|---|
| 0.60 | 1.305 | 1.351 | 1.462 | 0.3987 | 0.3254 | 0.2796 | 1.12 |
| 0.63 | 1.611 | 1.678 | 1.830 | 0.3425 | 0.2838 | 0.2479 | 1.14 |
| 0.66 | 1.863 | 1.945 | 2.122 | 0.3087 | 0.2449 | 0.2184 | 1.14 |
| 0.69 | 2.010 | 2.135 | 2.377 | 0.2621 | 0.2068 | 0.1873 | 1.18 |
| 0.72 | 2.049 | 2.213 | 2.543 | 0.2113 | 0.1720 | 0.1528 | 1.24 |
| 0.75 | 2.051 | 2.166 | 2.528 | 0.1783 | 0.1488 | 0.1259 | 1.23 |
| 0.78 | 1.729 | 1.799 | 2.202 | 0.1328 | 0.1079 | 0.0951 | 1.27 |

**Hüküm 5a — motif ÖLDÜ.** (τ−taban) motifi doğru olsaydı taban 0.58 sütunu τ=0.60'ta (τ−taban = 0.02) R ≈ 0.4 vermeliydi; 1.462 veriyor — 3.5 kat fazla. Üç eğri yatayda kaymıyor, yalnız çarpımsal olarak %12–27 yükseliyor. **R (τ−taban)'ın değil τ'nun fonksiyonudur; ilk bandın düşüklüğü taban-kenarı artefaktı değil, gerçek.** Aday (c) ampirik şekil olarak hayatta ama gerekçesi çürüdü.

**Hüküm 5b — ölçülmüş bir konvansiyon sistematiği.** R, tabanla çarpımsal olarak kayıyor (142'nin "taban-akışı"nın R'deki yüzü) ve kayma τ ile büyüyor: Δtaban = 0.12 için %12 (τ=0.60) → %27 (τ=0.78). Bu sistematik, jackknife hatasından (%0.5–3) bir mertebe büyüktür ve kapalı-form yarışında hata bütçesine katılmıştır.

---

## 6. İlkel gözlenebilir: φ_Γ(τ) ve sıfırı

R türetilmiş bir niceliktir; doğrudan ölçülen φ_Γ = arg Γ_rot'tur. Taban 0.46 koşusu, standart konvansiyonun ulaşamadığı τ<0.52 bantlarını açıyor ve orada **φ NEGATİF**:

| τ | 0.48 | 0.51 | 0.54 | 0.57 |
|---|---|---|---|---|
| φ_Γ (taban 0.46) | **−0.360** | **−0.036** | +0.275 | +0.592 |
| R (taban 0.46) | **−0.706** | **−0.099** | +0.421 | +0.904 |

Yerel parabol ile sıfırlar:

| nicelik / koşu | τ₀ |
|---|---|
| φ sıfırı — taban 0.46 (DOĞRUDAN, φ<0 bantlarıyla) | 0.5134 ± 0.0009 |
| φ sıfırı — taban 0.52, son pencere (ekstrapolasyon) | 0.5157 |
| φ sıfırı — taban 0.52, orta pencere (ekstrapolasyon) | 0.5126 |
| **R sıfırı — taban 0.46 (DOĞRUDAN)** | **0.5153 ± 0.0008** |

(İstatistiksel hatalar 20 000 örneklik Monte Carlo ile jackknife σ_R'lerden yayıldı.) R'nin sıfırı φ'nin sıfırından biraz yukarıdadır çünkü R=0'da M = M_emp'tir ve arg M_emp ≈ +0.03 (149: Im M_emp = 0.017–0.042).

**τ₀ = 0.5153 ± 0.0008(ist.)**; pencereler arası yayılım (0.5126 ↔ 0.5157) ~0.003 ekliyor, muhafazakâr toplam **±0.002**. Tabandan (0.46 ↔ 0.52) ve pencereden (L: 12.03 ↔ 11.46) bağımsız. **½ değildir**: 0.015 sapma, muhafazakâr hatayla ~7σ (yalnız istatistikle 18σ). 0.52 tabanı da değildir — taban 0.46 ile ölçüldü. Şimdilik açıklanmamış bir sabit.

---

## 7. Kapalı-form yarışı

Hata modeli: σ² = σ_jk² + σ_pencere² + σ_taban², σ_taban(τ) = \|R@0.58 − R@0.46\|/2 (§5'te ölçülmüş konvansiyon duyarlılığı). τ>0.76 bantları elendi (152'nin çözünmeyen-çizgi uyarısı). Skor: χ² = Σ((R−f)/σ)², AIC = χ² + 2k. n = 16 (8 τ × 2 pencere).

| aday | k | χ² | AIC | χ²/dof | parametreler |
|---|---|---|---|---|---|
| **g** R = c(τ−τ₀)/τ² | 2 | 4.66 | **8.66** | 0.33 | c=5.87, τ₀=0.5187 |
| **c** R = R∞(1−e^{−(τ−0.52)/w}) | 2 | 4.75 | **8.75** | 0.34 | R∞=2.74, w=0.1185 |
| c′ R∞(1−e^{−(τ−τ₀)/w}), τ₀=0.5139 ölçülen | 2 | 13.0 | 17.0 | 0.93 | R∞=3.65, w=0.201 |
| gτ R = c(τ−τ₀)/τ², τ₀=0.5139 ölçülen | **1** | 31.6 | 33.6 | 2.10 | c=5.23 |
| b R = c(τ−τ₀) doğrusal | 2 | 64.5 | 68.5 | 4.6 | c=12.3, τ₀=0.5036 |
| d″ erfc, τ_c ve Δ serbest (p=2) | 2 | 205 | 209 | 14.7 | τ_c=0.692, Δ=0.157 |
| g½ R = c(τ−½)/τ² (τ₀≡½ SABİT) | 1 | 237 | 239 | 15.8 | c=3.84 |
| d′ s·(−dlnρ_erfc/dA) | 1 | 409 | 411 | 27 | s=0.883 |
| **d −dlnρ_erfc/dA (0.68/0.125, p=2)** | **0** | **494** | **494** | 31 | — |
| e R = k·A²σΔ² | 1 | 899 | 901 | 60 | k=0.908 |
| a R = c (sabit) | 1 | 1319 | 1321 | 88 | c=0.567 |
| f1 −(1/2π)·dlnρ_ölçülen/dτ | 0 | 1552 | 1552 | 97 | — |
| f′ s·(−(1/2π)dlnρ_ölçülen/dτ) | 1 | 1272 | 1274 | 85 | s=0.790 |
| **f −(1/π)·dlnρ_ölçülen/dτ (iki-adım)** | **0** | **10541** | **10541** | 659 | — |
| d2 2·(−dlnρ_erfc/dA) | 0 | 8125 | 8125 | 508 | — |

Aynı sıralama, hata modeli değişse de bozulmuyor: yalnız istatistiksel hatalarla 8 nokta üzerinde (c) AIC 30.4 / (g) 46.6 / (d) 23 839 / (f) 6 417; 18 noktalık ince-ızgara çapraz sınavında (c) 394.5 / (g) 605 / (d) 61 554 / (f) 8 723. `R_gauss` tahmincisiyle koşulduğunda da aynı iki form önde ((c′) 19.7, (c) 24.6, (g) 82.7) ve aynı adaylar dipte.

### Artık yapıları (sistematikli koşu, son pencere)

| aday | 0.540 | 0.570 | 0.600 | 0.630 | 0.660 | 0.690 | 0.720 | 0.750 |
|---|---|---|---|---|---|---|---|---|
| c | −0.003 | −0.022 | +0.006 | +0.020 | +0.045 | +0.047 | −0.020 | −0.180 |
| g | −0.007 | −0.007 | +0.024 | +0.030 | +0.039 | +0.021 | −0.069 | −0.249 |
| d [0 par] | −0.012 | +0.180 | +0.184 | −0.037 | −0.431 | −1.003 | −1.772 | −2.736 |
| f [0 par] | −1.024 | −0.450 | −0.015 | +0.244 | +0.374 | +0.355 | +0.154 | −0.243 |

(c) ve (g)'nin artıkları yapısızdır (yalnız en üst bantta −0.18/−0.25'lik ortak kayma, çözünürlük bölgesinin başlangıcı). (d)'nin artığı τ ile tek yönlü ve patlayarak büyüyor — şekil hatası, ölçek hatası değil. (f)'nin artığı bir kemer çiziyor (−1.02 → +0.37 → −0.24).

### (g)'nin bağımsız öngörüsü

(g)'nin τ₀ parametresi yalnız τ ≥ 0.5375 bantlarından fit ediliyor, yani sıfır noktasının ÖLÇÜLDÜĞÜ bölgeyi (τ = 0.48–0.51, taban 0.46) hiç görmüyor. Fit sonuçları: **0.5155** (8 nokta, R_emp), **0.5187 ± 0.0011** (16 nokta, sistematikli), **0.5128** (18 nokta, ince). Doğrudan ölçüm: **0.5153 ± 0.0008**. Üç fit değeri ölçüleni kuşatıyor ve hepsi 0.003 içinde. Dürüst nitelendirme: **uyum niteliksel olarak çarpıcı, nicel olarak sınırda** — sistematikli fitin 0.5187'si, iki formal hata birleştirildiğinde ölçülen 0.5153'ten 2.4σ uzak. Yine de bu, herhangi bir adayın gösterebildiği tek bağımsız öngörüdür.

Buna karşılık aynı formun 1-parametreli τ₀ ≡ ½ sürümü (g½) çöküyor (χ²/dof = 15.8): **τ₀ ½ değil.**

---

## 8. HÜKÜM

**Hayatta olan:** iki 2-parametreli ampirik form, birbirinden ayırt edilemez (ΔAIC = 0.09).

- **(c)** R = R∞·(1 − e^{−(τ−0.52)/w}), R∞ = 2.74 ± 0.17, w = 0.1185 ± 0.0099
- **(g)** R = c·(τ−τ₀)/τ², c = 5.874 ± 0.147, τ₀ = 0.5187 ± 0.0011

(hatalar Δχ² = 1 profilinden, diğer parametre her adımda yeniden optimize edilerek)

(g)'nin lehindeki fazladan iki şey: (i) τ₀'ı, bağımsız bir bölgede doğrudan ölçülen 0.5153 ± 0.0008 ile 2.4σ içinde uyuşuyor; (ii) yorumu var — φ_Γ derinlikte (τ−τ₀) ile doğrusal, R'nin eğriliği ise tümüyle kinematik A⁻² = (2πτ)⁻² çarpanı. (c)'nin lehindeki tek şey: %5 daha küçük artık; motifi (taban-kenarı) ise §5'te çürüdü — yani (c) artık gerekçesiz bir eğri uydurmasıdır.

**Ölen:**
- **(d) sıfır-parametreli erfc türevi — ÖLDÜ.** τ=0.5375'teki isabeti (0.414 vs 0.391) tesadüftür; τ=0.74'te 4.59 öngörüyor, 2.12 ölçülüyor. τ_c ve Δ serbest bırakılsa bile (14.7 χ²/dof) şekil tutmuyor — aile reddedildi. *Not:* serbest fit τ_c = 0.69 (16 nokta) / 0.68 (8 nokta) buluyor, yani bağımsız olarak bilinen etkin kesimi geri veriyor; ama Δ ≈ 0.16–0.25 gerekiyor (merdiven kesiminin 1.3–2 katı) ve bununla bile uymuyor.
- **(f) ölçülmüş ρ'nun log-türevi (tek adım ve iki adım) — ÖLDÜ.** Oran sabit değil (0.57–2.52).
- **(a) sabit R = 6.20/2π — ÖLDÜ** (zaten "ölçek doğru şekil yanlış" biliniyordu; sayısal olarak χ²/dof = 88).
- **(e) R ∝ A²σΔ² — ÖLDÜ** (τ² çok yavaş büyüyor: τ 0.54→0.75'te A² 1.93 kat artarken R 5.1 kat artıyor).
- **(b) doğrusal — ÖLDÜ** (R içbükey; artık kemeri belirgin).
- **(τ−taban) motifi — ÖLDÜ** (taban taraması).

**Mühürlenmedi.** Görev "sıfır-parametreli aday tutarsa vurgula, hiçbiri tatmin etmiyorsa dürüstçe AÇIK yaz" diyordu: **sıfır-parametreli adayların hiçbiri tutmuyor; kapalı form AÇIKTIR.** Elde iki eşdeğer ampirik form var, birini seçecek delil yok.

**Avın yönü değişti.** Sentetik tanık (§4), R(τ)'nun tırmanan-doygunlaşan şeklinin soğurmasız bir asal-merdiven gazında da doğduğunu gösteriyor. Yani R'nin kapalı formu **soğurma profilinden türetilemez**; Γ_rot tahmincisinin kendi yapısından (η'nın bant-bant lag-1 çizgi korelatörü + ampirik ağırlıklı faz eşlemesi) türetilmelidir. Sıradaki analitik dişli budur; ve gerçek gazın sentetikten ayrıldığı iki sayı — τ₀ = 0.515 (sentetikte ≈0.486) ve φ/A² ölçeği (sentetikte 2.1 kat büyük) — o türetimin gerçek girdisi olmalıdır.

*Çapraz not (bu koşunun bulgusu değil):* 152'nin işaret ettiği "sıradaki aday: GUE seviye itmesi" bu kalemle eşzamanlı yürüyen **153** koşusunda sınanmış ve orada da kapanmamış görünüyor (`153_itmeli_gaz_RAPOR.md`, kısa hüküm: "hiçbiri 1.00'e inmiyor"). İki kalem aynı yöne bakıyor: sentetik gazın hangi istatistiğini düzeltirsen düzelt faz kapanmıyor, dolayısıyla eksik olan malzeme tahmincinin gördüğü *ikili yapının* kendisi. 153'ün sayıları bu raporda doğrulanmadı; yalnız yön uyumu kayda geçiriliyor.

---

## 9. Açıklar

1. **τ₀ = 0.516 ± 0.002 nedir?** ½ değil, taban değil, L'den bağımsız. En sağlam yeni sayı ve tümüyle açık. (Sentetik keskin gazda ≈0.487 — fark gerçek gazın ikili korelasyonundan geliyor olmalı.)
2. **Hangi R?** R_emp ve R_gauss τ=0.74'te %35 ayrışıyor. Kapalı form iddiası tahminciyi belirtmek zorunda. R_gauss ile en iyi form (c′) oluyor, R_emp ile (g)/(c) — sıralama tahminciye duyarlı.
3. **Konvansiyon sistematiği.** R, taban 0.46→0.58 ile %12–27 kayıyor. Hiçbir kapalı form bundan daha iyi konuşamaz; bir "yasa" iddiası için taban-bağımsız bir R tanımı gerekiyor.
4. **Yüksek τ güvenilmez.** τ=0.74'te ağırlığın etkin örneklemi n_eff = 17 361/300 000 (%5.8); τ=0.78'de R iki pencerede de düşüyor. R'nin doygunluğu kısmen tahmincinin derin-kuyruk ekstrapolasyonu olabilir — sentetik gazda aynı tahminci n_eff'i 362'ye düşürüp anlamsız değerler üretiyor. Doygunluk platosunun (R∞ ≈ 2.2–2.7) fiziksel mi tahminci kaynaklı mı olduğu AÇIK.
5. **A4'ün ρ'su gerçeğe uymuyor** (τ=0.74'te 0.0020 vs 0.1577) — 152'nin "S3'ü oturtan kesim" seçimi üçüncü bir gözlenebilir tarafından reddediliyor. Etkin merdiven kesiminin yeniden kalibrasyonu ayrı bir kalem.
6. **ρ'nun ince yapısı.** lnρ'nun 0.03'lük ızgaradaki yerel eğimi tekdüze değil (−4.70 / −3.57 / −4.57 …) ve bu düzensizlik İKİ PENCEREDE DE aynı — gürültü değil, yapı. Kaynağı açık.

---

## Ek — dosyalar ve tekrar-üretim

| dosya | ne |
|---|---|
| `154_configs/154_cekirdek.py` | tek ölçüm çekirdeği (150 + E1–E6); gerçek/taban/sentetik hepsi bunu çağırır |
| `154_configs/154_gercek.py` | `<son\|orta> <std\|ince>` — iki pencere |
| `154_configs/154_sentetik.py` | `<keskin\|A4> [std\|ince]` — sentetik tanık (152'nin merdiven+Newton'u birebir, S3 denetimli) |
| `154_configs/154_taban_taramasi.py` | `<0.46\|0.52\|0.58>` — taban taraması |
| `154_configs/154_yarisma.py` | kapalı-form yarışı + φ analizi + taban tablosu |
| `154_configs/154_figur.py` | 6 panelli figür |
| `154_R_kapali_form.png` | figür |

Ham çıktılar ve JSON'lar: scratchpad `154/` (`R_gercek_*.json`, `R_sentetik_*.json`, `R_taban_*.json`, `yarisma.json`, `yarisma_cikti.txt`, `log_*.txt`). Bu rapordaki her sayı o JSON/loglardan alınmıştır; elle girilen tek şey 150/152'nin referans değerleridir (denetim sütunlarında karşılaştırma için).

Koşu süreleri: gerçek pencereler 4.0 dk (std) / 6.6–6.8 dk (ince); taban taraması 3.7–6.3 dk; sentetik keskin 2.5 dk (152'nin z'si okundu), A4 21 dk (Newton yeniden çözüldü, `z_A4.npy` önbelleğe alındı — yeniden ölçüm 1.7 dk).
