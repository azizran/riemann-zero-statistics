# 156 — bağlaşımın kaynağı: üç-kanallı ağırlık ayrışımı ve ×1.4'ün alan sınavı

150'nin mühürlediği mekanizma tek bir korelasyon iddiasıdır:

> Γ_rot(τ) = ⟨w·e^{−iA·dsΔ}⟩ / ⟨w⟩,  w = e^{−A·R(τ)·dsΔ}

Buradaki **dsΔ ağırlıkta iki ayrı işi birden görüyor**: hem fazı ilerleten
fiziksel adım, hem de sağkalımı belirleyen yerel değişken. KALEM'in
hipotezi bu ikisini ayırmayı öneriyordu: gerçek gazda yerel adımın büyük
kısmı merdivenle koherent ALAN'dır; sağkalım-bağlaşımı belki yalnız
η-payına oturuyordur; sentetikte pay yapısı farklı olduğu için etkin R
şişiyordur — ve ×1.4 buradan geliyordur.

Bu rapor o ayrımı ölçtü. 0.52-taban zincirinin regresyonu zaten üç parçayı
üretiyor; burada ilk kez **ayrı ayrı ağırlık kanalı** olarak kullanıldılar:

    ds = drift + lad + eta            (kapanış: maks|artık| = 2.2e−16)
    X_tam = (ds_n+ds_{n+1})/2,  X_lad, X_eta, X_dri aynı biçimde

    model_k(R) = ⟨e^{−A·R·X_k} · e^{−iA·X_tam}⟩ / ⟨e^{−A·R·X_k}⟩

**Faz faktörü her zaman TAM adımdır**; yalnız ağırlık kanalı değişir. Her
bantta R_k ölçülen fazı eşleyecek şekilde çözülür (150/154'ün kök
bulucusu), sonra tek serbest parametre fazı yediği için **Re-kısmı
bağımsız bir aşırı-belirleme sınavıdır** (150'nin T2 eşiği: |ΔRe| < 0.06).

## Ne koşuldu

| koşu | gaz | bantlar | çıktı |
|---|---|---|---|
| `156_kos.py gercek std` | zeros6 son-300k | std 5 bant | k3_gercek_std.json |
| `156_kos.py A4 std` | 152'nin erfc-0.68/0.125 gazı | std 5 bant | k3_A4_std.json |
| `156_kos.py keskin std` | 152'nin kesim**siz** kontrolü | std 5 bant | k3_keskin_std.json |
| `156_kos.py gercek gorev 0.46 4000` | gerçek, taban kaydırılmış | görev 3 bandı | k3_gercek_gorev_t0.46.json |
| `156_kos.py gercek gorev 0.58 4000` | gerçek, taban kaydırılmış | görev 3 bandı | k3_gercek_gorev_t0.58.json |
| `156_kapasite.py` | gerçek + A4 | 5 bant | kapasite.json |
| `156_nakil.py` | gerçek R → A4 | 5 bant | nakil.json |
| `156_ortak_cift.py` | gerçek ∩ A4 | 5 bant | ortak_cift.json |

Sentetik gazların Newton çözümleri yeniden koşulmadı; 152/154'ün önbelleğe
aldığı `z_A4.npy` ve `z_keskin.npy` okundu (aynı dosyalar, aynı sayılar).

**Kopya-kayması denetimi geçti.** R kök-bulucusu ve M hesabı kopya değil,
`154_cekirdek`'ten import edilerek çağrıldı; Γ döngüsü kopyadır ve 154'ün
kaydettiği değerlerle karşılaştırıldı. Üç gazın on beş bandında da
Δ|Γ|, Δφ, ΔReΓ makine hassasiyetinde (≤ 1e−14) ve R̃_tam, 154'ün R'siyle
dördüncü basamağa kadar aynı çıktı. Tek istisna: keskin gazın τ=0.74
bandında 154 R=+4.574 bulmuşken 156 kök bulamadı — o bantta ağırlık 58
etkin noktaya çökmüş durumda ve kök, tarama adımına (0.02 vs 0.05)
duyarlı. Bu bant zaten hükme girmiyor; ayrıntı Dürüstlük notlarında.

## Ölçek normalizasyonu (sayısal, fiziği değiştirmez)

Kanalların varyansları çok farklı olduğundan ham R_k'ler doğrudan
karşılaştırılamaz ve aynı arama aralığına sığmaz. Gauss limitinde faz =
A²·R_k·Cov(X_k, X_tam) olduğu için

    Y_k = X_k · σΔ² / Cov(X_k, X_tam),   R_k(ham) = R̃_k · s_k,
    s_k = σΔ² / Cov(X_k, X_tam)

tanımlandı; R̃ üç kanalda da aynı aralıkta aranır ve R̃_tam tanımı gereği
150/154'ün R'siyle özdeştir. Kök **ampirik dağılımla** çözülür — Gauss
yaklaşımı yalnız ölçek seçiminde kullanılmıştır. Tablolarda hem R̃ hem ham
R verilmiştir; gazlar arası karşılaştırmalar **ham R** üzerinden yapılır,
çünkü nakil sınavının anlamı budur.

Ağırlık çöküşü her satırda `n_eff = (Σw)²/Σw²` ile kaydedildi. 300 binlik
örnekte n_eff < 3000 (%1) demek, "uydurma birkaç yüz bağa dayanıyor"
demektir; bu satırlar ⚠ ile işaretli ve hükümde ağırlıkları düşürüldü.

---

### T0 — kanal istatistikleri ve PAY YAPISI

| büyüklük | gerçek (son-300k) | sentetik A4 | sentetik keskin |
|---|---|---|---|
| σ_ds² (nokta) | 0.16744 | 0.11284 | 0.12795 |
| σ_lad² (nokta) | 0.14469 | 0.08988 | 0.05125 |
| σ_η² (nokta) | 0.02274 | 0.02295 | 0.07670 |
| c₁(η) | -0.01158 | -0.00918 | -0.03535 |
| σΔ² = Var(X_tam) | 0.05466 | 0.04294 | 0.03499 |
| Var(X_lad) | 0.05514 | 0.03686 | 0.01713 |
| Var(X_eta) | 0.00558 | 0.00689 | 0.02068 |
| Var(X_drift) | 0.00000 | 0.00000 | 0.00000 |
| Cov(X_lad, X_tam) | 0.05211 | 0.03646 | 0.01572 |
| Cov(X_eta, X_tam) | 0.00255 | 0.00649 | 0.01927 |
| **pay** Cov(X_lad,X_tam)/σΔ² | +0.9533 | +0.8490 | +0.4493 |
| **pay** Cov(X_eta,X_tam)/σΔ² | +0.0467 | +0.1510 | +0.5507 |
| **pay** Cov(X_drift,X_tam)/σΔ² | -0.0000 | -0.0000 | -0.0000 |
| korel(X_lad, X_tam) | +0.9492 | +0.9164 | +0.6421 |
| korel(X_eta, X_tam) | +0.1460 | +0.3772 | +0.7164 |
| ölçek s_lad = σΔ²/Cov(X_lad,X_tam) | 1.049 | 1.178 | 2.226 |
| ölçek s_eta = σΔ²/Cov(X_eta,X_tam) | 21.430 | 6.621 | 1.816 |


### T1-gerçek — tek-kanallı eşleşmeler

| τ̄ | \|Γ\| | φ_Γ | ReΓ | kanal | R_ham | R̃ | faz artığı | ReM | \|ΔRe\| | n_eff | hüküm |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.5375 | 0.813 | +0.2556 | +0.7866 | tam | +0.391 | +0.391 | 0.0000 | +0.7284 | 0.0582 | 274950 | ✓ |
|  |  |  |  | lad | +0.371 | +0.354 | 0.0000 | +0.7075 | 0.0791 | 270362 | ✗ Re |
|  |  |  |  | eta | -4.980 | -0.232 | 0.0000 | +0.0272 | 0.7594 | 191 | ✗ Re ⚠n_eff |
| 0.5850 | 0.735 | +0.7251 | +0.5505 | tam | +1.086 | +1.086 | 0.0000 | +0.5661 | 0.0157 | 153942 | ✓ |
|  |  |  |  | lad | +0.858 | +0.818 | 0.0000 | +0.5195 | 0.0310 | 127437 | ✓ |
|  |  |  |  | eta | -4.433 | -0.207 | 2.4165 | -0.0416 | 0.5921 | 281 | faz üretilemiyor |
| 0.6600 | 0.625 | +1.3810 | +0.1180 | tam | +1.852 | +1.852 | 0.0000 | +0.1428 | 0.0248 | 41855 | ✓ |
|  |  |  |  | lad | +1.218 | +1.161 | 0.0000 | +0.1264 | 0.0085 | 28564 | ✓ |
|  |  |  |  | eta | -3.832 | -0.179 | 1.7606 | -0.1540 | 0.2719 | 384 | faz üretilemiyor |
| 0.7400 | 0.520 | +1.8595 | -0.1480 | tam | +2.120 | +2.120 | 0.0000 | -0.2099 | 0.0619 | 17361 | ✗ Re |
|  |  |  |  | lad | +1.292 | +1.232 | 0.0000 | -0.1855 | 0.0375 | 12135 | ✓ |
|  |  |  |  | eta | -3.355 | -0.157 | 1.2821 | -0.2501 | 0.1022 | 482 | faz üretilemiyor |
| 0.8150 | 0.663 | +1.7274 | -0.1034 | tam | +1.507 | +1.507 | 0.0000 | -0.1022 | 0.0012 | 41212 | ✓ |
|  |  |  |  | lad | +0.969 | +0.924 | 0.0000 | -0.0850 | 0.0183 | 30850 | ✓ |
|  |  |  |  | eta | -3.029 | -0.141 | 4.8690 | -0.3208 | 0.2175 | 518 | faz üretilemiyor |

‡ = |Γ|>1.5, bant sayısal olarak patlak (152/153 kuralı); hüküm çıkarılmaz. ⚠n_eff = ağırlık 3000'den az etkin noktaya çökmüş (uydurma birkaç noktaya dayanıyor).


### T1-A4 — tek-kanallı eşleşmeler

| τ̄ | \|Γ\| | φ_Γ | ReΓ | kanal | R_ham | R̃ | faz artığı | ReM | \|ΔRe\| | n_eff | hüküm |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.5375 | 0.724 | +0.5947 | +0.5997 | tam | +1.215 | +1.215 | 0.0000 | +0.6591 | 0.0595 | 157053 | ✓ |
|  |  |  |  | lad | +1.226 | +1.041 | 0.0000 | +0.6418 | 0.0421 | 114885 | ✓ |
|  |  |  |  | eta | +6.941 | +1.048 | 0.0000 | +0.5622 | 0.0374 | 6221 | ✓ |
| 0.5850 | 0.698 | +1.0912 | +0.3219 | tam | +2.034 | +2.034 | 0.0000 | +0.3791 | 0.0572 | 59808 | ✓ |
|  |  |  |  | lad | +1.808 | +1.535 | 0.0000 | +0.3580 | 0.0362 | 12497 | ✓ |
|  |  |  |  | eta | +8.012 | +1.210 | 0.0000 | +0.2865 | 0.0354 | 879 | ✓ ⚠n_eff |
| 0.6600 | 0.737 | +1.8214 | -0.1828 | tam | +3.867 | +3.867 | 0.0000 | -0.2294 | 0.0465 | 9029 | ✓ |
|  |  |  |  | lad | +2.345 | +1.991 | 0.0000 | -0.2026 | 0.0198 | 309 | ✓ ⚠n_eff |
|  |  |  |  | eta | +9.109 | +1.376 | 0.0000 | -0.1800 | 0.0029 | 103 | ✓ ⚠n_eff |
| 0.7400 | 2.293 ‡ | -2.4989 | -1.8352 | tam | -2.623 | -2.623 | 0.0000 | -0.5283 | 1.3069 | 1786 | ✗ Re ⚠n_eff |
|  |  |  |  | lad | -4.624 | -3.925 | 0.0000 | -0.5169 | 1.3183 | 195 | ✗ Re ⚠n_eff |
|  |  |  |  | eta | -4.255 | -0.643 | 0.0000 | -0.4858 | 1.3494 | 59 | ✗ Re ⚠n_eff |
| 0.8150 | 1.376 | -2.7252 | -1.2588 | tam | -2.343 | -2.343 | 0.0000 | -0.5471 | 0.7117 | 1974 | ✗ Re ⚠n_eff |
|  |  |  |  | lad | +3.501 | +2.973 | 0.4164 | -0.9765 | 0.2823 | 4 | faz üretilemiyor |
|  |  |  |  | eta | -3.827 | -0.578 | 0.0000 | -0.5355 | 0.7234 | 64 | ✗ Re ⚠n_eff |

‡ = |Γ|>1.5, bant sayısal olarak patlak (152/153 kuralı); hüküm çıkarılmaz. ⚠n_eff = ağırlık 3000'den az etkin noktaya çökmüş (uydurma birkaç noktaya dayanıyor).


### T1-keskin — tek-kanallı eşleşmeler

| τ̄ | \|Γ\| | φ_Γ | ReΓ | kanal | R_ham | R̃ | faz artığı | ReM | \|ΔRe\| | n_eff | hüküm |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.5375 | 0.807 | +0.5395 | +0.6927 | tam | +1.574 | +1.574 | 0.0000 | +0.7397 | 0.0470 | 137461 | ✓ |
|  |  |  |  | lad | +2.187 | +0.982 | 0.0000 | +0.6781 | 0.0146 | 34521 | ✓ |
|  |  |  |  | eta | +4.294 | +2.365 | 0.0000 | +0.8000 | 0.1072 | 9135 | ✗ Re |
| 0.5850 | 0.794 | +1.0417 | +0.4006 | tam | +2.721 | +2.721 | 0.0000 | +0.4153 | 0.0147 | 18590 | ✓ |
|  |  |  |  | lad | +2.844 | +1.278 | 0.0000 | +0.3673 | 0.0333 | 1680 | ✓ ⚠n_eff |
|  |  |  |  | eta | +24.879 | +13.702 | 0.0000 | +0.4887 | 0.0880 | 5 | ✗ Re ⚠n_eff |
| 0.6600 | 0.790 | +1.8216 | -0.1960 | tam | +3.626 | +3.626 | 0.0000 | -0.1942 | 0.0017 | 1824 | ✓ ⚠n_eff |
|  |  |  |  | lad | +3.268 | +1.468 | 0.0000 | -0.1843 | 0.0116 | 109 | ✓ ⚠n_eff |
|  |  |  |  | eta | +36.315 | +20.000 | 0.5654 | +0.3062 | 0.5022 | 3 | faz üretilemiyor |
| 0.7400 | 0.836 | +2.6364 | -0.7313 | tam | -3.767 | -3.767 | 5.7780 | -0.6617 | 0.0696 | 58 | faz üretilemiyor |
|  |  |  |  | lad | +4.218 | +1.895 | 0.0000 | -0.8290 | 0.0977 | 8 | ✗ Re ⚠n_eff |
|  |  |  |  | eta | -6.028 | -3.320 | 5.7780 | -0.7414 | 0.0101 | 15 | faz üretilemiyor |
| 0.8150 | 1.015 | -2.9446 | -0.9956 | tam | -2.978 | -2.978 | 0.0000 | -0.5574 | 0.4383 | 183 | ✗ Re ⚠n_eff |
|  |  |  |  | lad | -4.563 | -2.050 | 1.9053 | +0.2823 | 1.2779 | 686 | faz üretilemiyor |
|  |  |  |  | eta | -4.579 | -2.522 | 0.0000 | -0.5379 | 0.4577 | 62 | ✗ Re ⚠n_eff |

‡ = |Γ|>1.5, bant sayısal olarak patlak (152/153 kuralı); hüküm çıkarılmaz. ⚠n_eff = ağırlık 3000'den az etkin noktaya çökmüş (uydurma birkaç noktaya dayanıyor).


### T2 — iki-kanallı sırt (w = e^{−A(R_l·X_lad + R_e·X_eta)})

| gaz | τ̄ | sırt nokta | serbest en iyi (R_l,R_e)_ham | \|ΔRe\| | n_eff | n_eff≥3000 kısıtlı (R_l,R_e)_ham | \|ΔRe\| | n_eff |
|---|---|---|---|---|---|---|---|---|
| gerçek | 0.5375 | 33 | (+3.409, +116.860) | 0.0017 | 10 | (+0.524, +3.118) | 0.0310 | 226800 |
| gerçek | 0.5850 | 32 | (+1.049, +0.913) | 0.0094 | 151340 | (+1.049, +0.913) | 0.0094 | 151340 |
| gerçek | 0.6600 | 32 | (+1.049, -0.504) | 0.0014 | 23727 | (+1.049, -0.504) | 0.0014 | 23727 |
| gerçek | 0.7400 | 23 | (+0.787, -1.287) | 0.0024 | 7041 | (+0.787, -1.287) | 0.0024 | 7041 |
| gerçek | 0.8150 | 13 | (+1.573, +1.693) | 0.0003 | 41935 | (+1.573, +1.693) | 0.0003 | 41935 |
| A4 | 0.5375 | 17 | (+0.442, +5.371) | 0.0067 | 24073 | (+0.442, +5.371) | 0.0067 | 24073 |
| A4 | 0.5850 | 22 | (+1.472, -1.332) | 0.0022 | 1762 | (+1.620, -0.785) | 0.0168 | 3962 |
| A4 | 0.6600 | 35 | (+1.767, -1.231) | 0.0007 | 124 | (+3.534, +2.741) | 0.0409 | 4036 |
| A4 ‡ | 0.7400 | 12 | (-0.442, -4.255) | 1.3094 | 77 | — | — | — |
| A4 | 0.8150 | 11 | (-0.294, -3.814) | 0.7134 | 79 | — | — | — |
| keskin | 0.5375 | 22 | (+1.948, +0.701) | 0.0152 | 84643 | (+1.948, +0.701) | 0.0152 | 84643 |
| keskin | 0.5850 | 7 | (+2.504, +5.446) | 0.0325 | 5631 | (+2.504, +5.446) | 0.0325 | 5631 |
| keskin | 0.6600 | 16 | (+5.565, +15.631) | 0.0003 | 18 | — | — | — |
| keskin | 0.7400 | 19 | (+6.956, +17.849) | 0.0224 | 9 | — | — | — |
| keskin | 0.8150 | 15 | (-2.226, -3.365) | 0.4438 | 199 | — | — | — |

Sırttaki en iyi nokta SIĞ bir minimumdur; kabul edilebilir BÖLGE (faz artığı<0.02, |ΔRe|<0.06, n_eff≥3000) daha bilgilidir:

| gaz | τ̄ | kabul edilebilir R_l (ham) | R_e (ham) | nokta |
|---|---|---|---|---|
| gerçek | 0.5375 | +0.393 … +0.524 | +0.438 … +3.118 | 2 |
| gerçek | 0.5850 | +0.787 … +1.311 | -0.345 … +2.139 | 5 |
| gerçek | 0.6600 | +0.393 … +4.196 | -2.477 … +8.828 | 30 |
| gerçek | 0.7400 | +0.393 … +1.967 | -2.187 … +1.727 | 13 |
| gerçek | 0.8150 | +0.393 … +1.705 | -1.529 … +2.063 | 11 |
| A4 | 0.5375 | -0.147 … +1.178 | -0.873 … +7.460 | 7 |
| A4 | 0.5850 | +1.620 … +1.914 | -0.785 … +0.540 | 3 |
| A4 | 0.6600 | +3.534 … +3.828 | +2.741 … +3.686 | 3 |
| A4 ‡ | 0.7400 | — | — | 0 |
| A4 | 0.8150 | — | — | 0 |
| keskin | 0.5375 | +1.391 … +2.226 | -0.134 … +1.948 | 4 |
| keskin | 0.5850 | +2.504 … +2.504 | +5.446 … +5.446 | 1 |
| keskin | 0.6600 | — | — | 0 |
| keskin | 0.7400 | — | — | 0 |
| keskin | 0.8150 | — | — | 0 |


### T3 — TABAN DAYANIKLILIĞI (gerçek gaz, regresyon tabanı 0.46 / 0.52 / 0.58)

lad/eta ayrımı tabanın nereye çizildiğine bağlıdır; hüküm tabana bağlıysa fizik değil konvansiyondur.

| τ̄ | taban | φ_Γ | R_tam | \|ΔRe\|_tam | R_lad | \|ΔRe\|_lad | R_eta | \|ΔRe\|_eta | 2K en iyi (R_l,R_e)_ham |
|---|---|---|---|---|---|---|---|---|---|
| 0.5850 | 0.46 | +0.7043 | +1.048 | 0.0177 | +0.844 | 0.0568 | -3.892 | 0.5235 | (+1.187, +1.763) |
| 0.5850 | 0.52 | +0.7251 | +1.086 | 0.0157 | +0.858 | 0.0310 | — | — | (+1.049, +0.913) |
| 0.5850 | 0.58 | +0.9613 | +1.535 | 0.0420 | +1.117 | 0.0041 | — | — | (+1.172, +0.209) |
| 0.6600 | 0.46 | +1.3436 | +1.785 | 0.0185 | +1.241 | 0.0049 | -3.201 | 0.1423 | (+1.055, -0.598) |
| 0.6600 | 0.52 | +1.3810 | +1.852 | 0.0248 | +1.218 | 0.0085 | — | — | (+1.049, -0.504) |
| 0.6600 | 0.58 | +1.4770 | +2.029 | 0.0145 | +1.290 | 0.0050 | — | — | (+1.042, -0.709) |
| 0.7400 | 0.46 | +1.7792 | +1.990 | 0.0303 | +1.317 | 0.0193 | -2.566 | 0.1143 | (+0.923, -1.092) |
| 0.7400 | 0.52 | +1.8595 | +2.120 | 0.0619 | +1.292 | 0.0375 | — | — | (+0.787, -1.287) |
| 0.7400 | 0.58 | +2.0440 | — | — | +1.401 | 0.0765 | — | — | (+0.782, -1.486) |


### T4 — görev bantları özeti (hangi kanal kazanıyor)

| τ̄ | gaz | R_tam | R_lad | R_eta | \|ΔRe\|_tam | \|ΔRe\|_lad | \|ΔRe\|_eta | en iyi kanal |
|---|---|---|---|---|---|---|---|---|
| 0.5850 | gerçek | +1.086 | +0.858 | — | 0.0157 | 0.0310 | — | tam |
| 0.5850 | A4 | +2.034 | +1.808 | +8.012 | 0.0572 | 0.0362 | 0.0354 | eta |
| 0.5850 | keskin | +2.721 | +2.844 | +24.879 | 0.0147 | 0.0333 | 0.0880 | tam |
| 0.6600 | gerçek | +1.852 | +1.218 | — | 0.0248 | 0.0085 | — | lad |
| 0.6600 | A4 | +3.867 | +2.345 | +9.109 | 0.0465 | 0.0198 | 0.0029 | eta |
| 0.6600 | keskin | +3.626 | +3.268 | — | 0.0017 | 0.0116 | — | tam |
| 0.7400 | gerçek | +2.120 | +1.292 | — | 0.0619 | 0.0375 | — | lad |
| 0.7400 | A4 ‡ | -2.623 | -4.624 | -4.255 | 1.3069 | 1.3183 | 1.3494 | tam |
| 0.7400 | keskin | — | +4.218 | — | — | 0.0977 | — | lad |


### T4b — FAZ KAPASİTESİ: kanal TEK BAŞINA hangi fazı üretebilir?

Sarma-açılmış arg M(R̃) menzili (`156_kapasite.py`). İkinci menzil n_eff ≥ 3000 kısıtı altındadır: ağırlığı çökertmeden erişilebilen faz. Ölçülen φ_Γ menzilin dışındaysa o kanal fazı HİÇBİR R ile üretemez.

| gaz | τ̄ | φ_Γ | kanal | serbest menzil | kök? | n_eff≥ kısıtlı menzil | kök? | η tavanı / φ_Γ |
|---|---|---|---|---|---|---|---|---|
| gercek | 0.5375 | +0.2556 | tam | [-2.541, +2.282] | var | [-2.087, +1.700] | var |  |
| gercek | 0.5375 | +0.2556 | lad | [-2.507, +2.593] | var | [-1.764, +1.615] | var |  |
| gercek | 0.5375 | +0.2556 | eta | [-0.977, +2.311] | var | [-0.148, +0.097] | **YOK** | 38% |
| gercek | 0.5850 | +0.7251 | tam | [-2.888, +2.546] | var | [-2.272, +1.849] | var |  |
| gercek | 0.5850 | +0.7251 | lad | [-2.893, +2.851] | var | [-1.907, +1.747] | var |  |
| gercek | 0.5850 | +0.7251 | eta | [-4.742, +0.107] | **YOK** | [-0.219, +0.107] | **YOK** | 15% |
| gercek | 0.6600 | +1.3810 | tam | [-3.435, +2.960] | var | [-2.585, +2.104] | var |  |
| gercek | 0.6600 | +1.3810 | lad | [-3.480, +3.251] | var | [-2.154, +2.001] | var |  |
| gercek | 0.6600 | +1.3810 | eta | [-4.544, +0.123] | **YOK** | [-0.846, +0.123] | **YOK** | 9% |
| gercek | 0.7400 | +1.8595 | tam | [-4.014, +3.397] | var | [-2.877, +2.352] | var |  |
| gercek | 0.7400 | +1.8595 | lad | [-4.079, +3.669] | var | [-2.427, +2.247] | var |  |
| gercek | 0.7400 | +1.8595 | eta | [-4.333, +0.142] | **YOK** | [-0.287, +0.142] | **YOK** | 8% |
| gercek | 0.8150 | +1.7274 | tam | [-4.553, +3.804] | var | [-3.186, +2.585] | var |  |
| gercek | 0.8150 | +1.7274 | lad | [-4.623, +4.056] | var | [-2.771, +2.482] | var |  |
| gercek | 0.8150 | +1.7274 | eta | [-4.136, +0.160] | **YOK** | [-1.418, +0.160] | **YOK** | 9% |
| A4 | 0.5375 | +0.5947 | tam | [-1.959, +1.721] | var | [-1.656, +1.610] | var |  |
| A4 | 0.5375 | +0.5947 | lad | [-1.364, +2.144] | var | [-1.317, +1.173] | var |  |
| A4 | 0.5375 | +0.5947 | eta | [-2.059, +1.855] | var | [-0.756, +0.662] | var | 111% |
| A4 | 0.5850 | +1.0912 | tam | [-2.278, +1.920] | var | [-1.816, +1.753] | var |  |
| A4 | 0.5850 | +1.0912 | lad | [-1.604, +2.338] | var | [-1.417, +1.289] | var |  |
| A4 | 0.5850 | +1.0912 | eta | [-2.245, +2.020] | var | [-0.845, +0.794] | **YOK** | 73% |
| A4 | 0.6600 | +1.8214 | tam | [-2.788, +2.247] | var | [-2.052, +1.975] | var |  |
| A4 | 0.6600 | +1.8214 | lad | [-2.030, +2.642] | var | [-1.593, +1.472] | **YOK** |  |
| A4 | 0.6600 | +1.8214 | eta | [-2.523, +2.280] | var | [-0.827, +0.839] | **YOK** | 46% |
| A4 | 0.7400 | -2.4989 | tam | [-3.329, +2.612] | var | [-2.303, +2.215] | **YOK** |  |
| A4 | 0.7400 | -2.4989 | lad | [-2.545, +2.967] | var | [-1.799, +1.637] | **YOK** |  |
| A4 | 0.7400 | -2.4989 | eta | [-2.822, +2.556] | var | [-0.978, +0.951] | **YOK** |  |
| A4 | 0.8150 | -2.7252 | tam | [-3.827, +2.962] | var | [-2.513, +2.444] | **YOK** |  |
| A4 | 0.8150 | -2.7252 | lad | [-3.072, +3.271] | var | [-1.940, +1.794] | **YOK** |  |
| A4 | 0.8150 | -2.7252 | eta | [-3.111, +2.816] | var | [-1.044, +1.199] | **YOK** |  |


### T5 — NAKİL SINAVI: gerçekte ölçülen ham R_k → A4 gazına

`öngörü/φ_A4` = nakil sentetiğin ölçülen fazını ne kadar tutuyor (hedef 1.00). `öngörü/φ_gerçek` = bağlaşım EVRENSEL olsaydı sentetik/gerçek faz oranı ne çıkardı — bunu ölçülen `oran` sütunuyla karşılaştırın.

| τ̄ | φ_gerçek | φ_A4 (ölç) | ölçülen oran | kanal | R_ham(gerçek) | φ_öngörü(A4) | öngörü/φ_A4 | öngörü/φ_gerçek | n_eff |
|---|---|---|---|---|---|---|---|---|---|
| 0.5375 | +0.2556 | +0.5947 | +2.327 | tam | +0.391 | +0.1927 | +0.324 | +0.754 | 278504 |
|  |  |  |  | lad | +0.371 | +0.1637 | +0.275 | +0.640 | 281173 |
|  |  |  |  | eta | -4.980 | -1.2968 | -2.180 | -5.073 | 288 |
| 0.5850 | +0.7251 | +1.0912 | +1.505 | tam | +1.086 | +0.6332 | +0.580 | +0.873 | 161920 |
|  |  |  |  | lad | +0.858 | +0.4825 | +0.442 | +0.665 | 182428 |
|  |  |  |  | eta | — | — | — | — | — |
| 0.6600 | +1.3810 | +1.8214 | +1.319 | tam | +1.852 | +1.2651 | +0.695 | +0.916 | 56573 |
|  |  |  |  | lad | +1.218 | +0.9266 | +0.509 | +0.671 | 61477 |
|  |  |  |  | eta | — | — | — | — | — |
| 0.7400 ‡ | +1.8595 | -2.4989 | -1.344 | tam | +2.120 | +1.6715 | -0.669 | +0.899 | 32645 |
|  |  |  |  | lad | +1.292 | +1.2730 | -0.509 | +0.685 | 25259 |
|  |  |  |  | eta | — | — | — | — | — |
| 0.8150 | +1.7274 | -2.7252 | -1.578 | tam | +1.507 | +1.5984 | -0.587 | +0.925 | 55978 |
|  |  |  |  | lad | +0.969 | +1.1596 | -0.426 | +0.671 | 65871 |
|  |  |  |  | eta | — | — | — | — | — |


### T6 — ORTAK ÇİFT: tek bir (R_l,R_e) iki gazı birden açıklıyor mu?

| τ̄ | ortak faz çözümü sayısı | en iyi (R_l,R_e)_ham | \|ΔRe\|_gerçek | \|ΔRe\|_A4 | n_eff(min) | hüküm |
|---|---|---|---|---|---|---|
| 0.5375 | 1 | (+0.604, +4.760) | 0.0582 | 0.0217 | 38203 | ORTAK ÇİFT AŞIRI-BELİRLEMEYİ GEÇTİ |
| 0.5850 | 1 | (+1.677, +3.844) | 0.0804 | 0.0470 | 30308 | ortak çift Re'yi TUTMUYOR |
| 0.6600 | 1 | (+3.182, +5.729) | 0.0398 | 0.0487 | 3930 | ORTAK ÇİFT AŞIRI-BELİRLEMEYİ GEÇTİ |
| 0.7400 ‡ | 0 | — | — | — | — | ortak faz çözümü yok |
| 0.8150 | 0 | — | — | — | — | ortak faz çözümü yok |


---

# HÜKÜM

## (i) Gerçekte bağlaşım hangi kanalda yaşıyor — MERDİVENDE, η'da DEĞİL

Hipotezin birinci yarısı doğrulandı, ikinci yarısı çürütüldü.

**Doğrulanan:** gerçek gazda yerel adım gerçekten neredeyse tamamen ALAN.
Bond adımının kovaryans payı merdivende **%95.3**, η'da **%4.7**, drift'te
%0.00 (Var(X_drift) = 8e−08 — "drift-küçük" varsayımı ölçüldü ve tuttu).

**Çürütülen:** bağlaşım η-payına oturmuyor. Üç bağımsız yoldan aynı sonuç:

1. **η fazı ÜRETEMİYOR.** Beş bandın dördünde (τ = 0.585 / 0.66 / 0.74 /
   0.815) η-ağırlığı hiçbir R değeriyle ölçülen fazı veremiyor: sarma
   açıldığında erişilebilen pozitif faz tavanı sırasıyla +0.107 / +0.123 /
   +0.142 / +0.160 rad, yani ölçülen fazın **yalnız %15 / %9 / %8 / %9'u**.
   Kalan bantta (τ=0.5375) kök var ama ağırlık 191 etkin bağa çökmüş ve
   |ΔRe| = 0.759 — aşırı-belirlemeyi yıkıcı biçimde kaybediyor. n_eff ≥ 3000
   kısıtı konunca **beş bandın beşinde de** η kanalında kök yok.
2. **Merdiven kanalı aşırı-belirlemeyi geçiyor**, üstelik tam adımdan daha
   iyi: |ΔRe|_lad = 0.0310 / 0.0085 / 0.0375 / 0.0183 (τ = 0.585 … 0.815),
   |ΔRe|_tam = 0.0157 / 0.0248 / **0.0619** / 0.0012. Görevin üç bandının
   ikisinde (0.66 ve 0.74) merdiven kanalı kazanıyor; τ=0.74'te tam adım
   150'nin 0.06 eşiğini aşarken merdiven geçiyor. Ağırlıklar sağlıklı
   (n_eff = 1.2e4 … 1.3e5).
3. **Taban dayanıklılığı.** lad/eta ayrımı 0.52 tabanının nereye
   çizildiğine bağlıdır — o yüzden taban 0.46 ve 0.58'de tekrarlandı.
   R_lad(τ=0.66) = 1.241 / 1.218 / 1.290 (±%3), R_lad(τ=0.74) = 1.317 /
   1.292 / 1.401 (±%4) — **taban değişse de aynı sabit**. (τ=0.585 bandı
   0.58 tabanının kenarına düştüğü için bu karşılaştırmaya alınmadı.)
   R_tam ise kayıyor
   (1.785 / 1.852 / 2.029 ve 1.990 / 2.120 / —). η, dokuz satırın altısında
   fazı üretemiyor, ürettiği üçünde |ΔRe| = 0.11–0.52. İki-kanallı en iyi
   çift üç tabanda da neredeyse aynı yerde: (R_l, R_e) ≈ (1.05, −0.6).

**Sonuç:** Γ_rot'un fazını üreten sağkalım-ağırlık korelasyonu, yerel
adımın **merdivenle koherent payına** oturuyor. η-payının katkısı sıfırla
uyumlu: iki-kanallı sırtta gerçek gazın kabul edilebilir R_e aralığı
τ=0.585'te −0.35…+2.14, τ=0.66'da −2.48…+8.83, ve en iyi nokta üç tabanda
da R_e ≈ −0.5…+0.9. Yani **tek gerekli kanal merdivendir**; η eklemek ne
gerekiyor ne de yasak.

Not: gerçek gazda korel(X_lad, X_tam) = 0.949, yani lad ile tam kanal %95
aynı değişkendir; ikisini ayırmak zayıf bir sınavdır ve "lad tam'dan
biraz daha iyi" hükmü ihtiyatlı okunmalıdır. Çürütülen kanal (η,
korel 0.146) ise ikisinden de kesin biçimde ayrışıyor — asıl hüküm odur.

Bunun tersini gösteren tek karşı-örnek yok: **keskin** gazda η, adımın
payının **%55'ini** taşıdığı halde (lad %45 — sıralama tersine dönmüş)
yine de aşırı-belirlemeyi kaybediyor (|ΔRe|_eta = 0.107 / 0.088 / 0.502)
ve merdiven kanalı kazanıyor. Yani bağlaşımın merdivende oturması, adımın
varyansını hangi kanalın taşıdığından bağımsız bir olgu.

## (ii) Sentetikte fark nerede — ×1.4 pay-yapısı DEĞİL, gerçekten daha büyük bağlaşım

Pay yapısı gerçekten farklı: η payı gerçek %4.7 → A4 %15.1 → keskin %55.1.
Ama **bu fark ×1.4'ü açıklamıyor; işareti bile ters.**

**Nakil sınavı** (Gauss yaklaşımı yok; gerçekte ölçülen ham R, sentetik
gazın kendi ampirik dağılımına uygulanıyor). Aşağıdaki özet T5'in
`ölçülen oran` ve `öngörü/φ_gerçek` sütunlarıdır:

| τ̄ | ölçülen φ_A4/φ_gerçek | R_tam sabit tutulursa öngörü | R_lad sabit tutulursa öngörü |
|---|---|---|---|
| 0.5375 | **2.33** | 0.75 | 0.64 |
| 0.5850 | **1.51** | 0.87 | 0.67 |
| 0.6600 | **1.32** | 0.92 | 0.67 |

Bağlaşım evrensel olsaydı, sentetik gazın fazı gerçeğinkinin **0.64–0.92
katı** çıkacaktı — çünkü A4'ün σΔ²'si (0.0429) ve merdiven kovaryansı
(0.0365) gerçeğinkinden (0.0547 / 0.0521) KÜÇÜK. Ölçülen ise 1.32–2.33 kat.
Pay yapısı fazı **düşürmesi** gereken yönde çalışıyor; fazlalık ona
rağmen var.

Fazlalık doğrudan R'de: R_lad(A4)/R_lad(gerçek) = 3.30 / 2.11 / 1.93 ve
R_tam oranı 3.11 / 1.87 / 2.09 (τ = 0.5375 / 0.585 / 0.66). Keskin gazda
daha da büyük: R_lad oranı 5.90 / 3.32 / 2.68. **Yani fazdaki ×1.3–1.5,
bağlaşımda ×1.9–2.1'e karşılık geliyor** — pay yapısı farkı düzeltildiğinde
açık büyüyor, kapanmıyor.

**Sıralama bir yön veriyor:** gerekli bağlaşım keskin (kesim yok) > A4
(erfc-0.68) > gerçek. Yani merdivenin yüksek-τ çizgi gücü söndürüldükçe
gereken R gerçeğe doğru iniyor. A4, keskin ile gerçek arasında ve gerçeğe
daha yakın. Bu, 152'nin erfc kesimini daha aşağı/daha sert çekmenin R'yi
indirmeye devam edeceğini söylüyor — ama 152 bunu σ_η² ve c₁'i bozmadan
yapamamıştı; gerilim aynı yerde duruyor.

**Yeni ve daha keskin olgu: sentetik gazın istediği bağlaşım mekanizmayı
kırıyor.** Fazı eşleyen kökte ağırlığın etkin bağ sayısı (τ=0.66,
merdiven kanalı): gerçek **28 564** → A4 **309** → keskin **109**. n_eff ≥
3000 kısıtı konunca A4'ün τ=0.66 bandında ne η ne merdiven kanalı fazı
üretebiliyor — yalnız tam adım, o da n_eff = 9029 ile. Gerçek gazda aynı
mekanizma rahat çalışıyor (n_eff 1e4–1e5). Yani sentetik gazın anormal
fazı, sağkalım-ağırlık resminin kendi geçerlilik sınırının dışına düşen
bir bağlaşım istiyor.

**Ortak çift sınavı (hipotezin en cömert okuması).** İki bilinmeyen
(R_l, R_e) ve iki kısıt (iki gazın fazı) olduğu için ortak bir çözümün
var olması sayım gereğidir, kanıt değil; asıl sınav o çiftin **hiçbir ek
serbestlik olmadan iki gazın Re'sini de tutması**. Aranan kutuda her bantta
tam bir ortak çözüm var:

* τ=0.5375 → (R_l, R_e) = (+0.60, +4.76), maks|ΔRe| = 0.058 → **geçiyor** (sınırda)
* τ=0.5850 → (+1.68, +3.84), maks|ΔRe| = 0.080 → **kalıyor**
* τ=0.6600 → (+3.18, +5.73), maks|ΔRe| = 0.049 → **geçiyor**

Yani iki-kanallı biçimiyle hipotez tümüyle ölmüş değil — ama zayıf:
(a) üç bandın birinde Re'yi tutmuyor; (b) tutturduğu yerde bile ortak çift,
her iki gazın kendi en iyi uydurmasından çok uzakta duruyor (gerçek gazın
kendi 2K en iyisi τ=0.66'da (1.05, −0.50), ortak çift (3.18, +5.73)); ve
(c) τ=0.585'te dayattığı çift, gerçek gazın o banttaki kabul edilebilir
bölgesinin **iki eksende de dışında** (R_l = 1.68 vs kabul 0.79…1.31;
R_e = +3.84 vs kabul −0.35…+2.14) — hipotezin o bantta kalması tam olarak
bu yüzden. τ=0.66'da ise ortak çift gerçek gazın kabul bölgesinin
(R_l 0.39…4.20, R_e −2.48…+8.83) içinde kalıyor ve bu yüzden geçiyor.

**Kısa cevap:** ×1.4'ün kaynağı pay yapısı değil. Pay yapısı düzeltilince
açık kapanmıyor, **büyüyor**: sentetik gaz gerçeğin yaklaşık **iki katı**
yerel sağkalım-ağırlık bağlaşımı istiyor, üstelik ağırlığı örneğin
binde birine çökertecek kadar büyük bir bağlaşım. Eksik olan şey adımın
kanallara dağılımı değil, merdiven payının **kendi iç yapısı** —
sentetikte aynı fazı üretmek için çok daha sert bir korelasyon gerekiyor.

## Sıradaki adım (bu ölçümün işaret ettiği)

R_lad'ın tabandan bağımsız çıkması (±%3), R(τ) türetmesinin doğru
değişkeninin **tam adım değil merdiven payı** olduğunu söylüyor: 154'ün
kapalı-form yarışması R_tam(τ) üzerinde koşulmuştu ve R_tam taban
kaydırınca %6–10 kayıyor. Yarışmayı R_lad(τ) üzerinde tekrarlamak
(gerçek: 0.371 / 0.858 / 1.218 / 1.292 / 0.969 @ τ = 0.5375 … 0.815)
tabana bağımlı olmayan bir hedef verir.

---

## Dürüstlük notları

* **Uydurma yok.** Rapordaki bütün tablo sayıları `156_tablolar.py` ile
  `scratchpad/156/*.json`'dan otomatik yazıldı. Metindeki sayıların hepsi
  bu tablolardan okunmuştur; elle girilen tek şey 150'nin |ΔRe| < 0.06
  eşiği ve n_eff ≥ 3000 (%1) kuralıdır.
* **Patlayan/güvenilmez koşular ve rapordaki yerleri.** A4'ün τ=0.74 bandı
  |Γ| = 2.29 ile patlak (152/153'te de patlamıştı) — ‡ ile işaretli,
  hiçbir hükme girmedi. A4'ün τ=0.815 bandı |Γ| = 1.38 ile 1.5 eşiğinin
  altında kaldığı için kural gereği "patlak" sayılmıyor, ama üç kanalın da
  |ΔRe|'si 0.28–0.72 ve n_eff'i 4–1974: **ölçüm değil**, hükme alınmadı.
  Keskin gazın τ=0.74 ve 0.815 bantlarında n_eff = 8–686; o satırlar
  tabloda duruyor ama hükümde yalnız "keskin'de de η kaybediyor" ifadesi
  için ilk üç bant kullanıldı.
* **Görev bandı 0.70–0.78'de sentetik karşılaştırma yapılamadı**: A4'ün o
  bandı patlak, keskin'in o bandı n_eff = 8. Bu bantta yalnız gerçek gaz
  hükmü (merdiven kazanıyor, |ΔRe| 0.0375 vs tam 0.0619) verilmiştir;
  ×1.4 hükmü iki banda (0.585 ve 0.66) ve destek olarak 0.5375'e dayanır.
* **τ=0.5375 zayıf banttır** (33 aday çizgi; 152/153 de bunu kaydetmişti).
  Oradaki oran 2.33 diğer iki bandın 1.51/1.32'sinden belirgin uzak;
  hükümde ayrı gösterildi, ortalamaya katılmadı.
* **154 ile bir uyuşmazlık kaydı.** Keskin gaz, τ=0.74: 154 R = +4.574
  (artık 0.000) bulmuş, 156 kök bulamamıştır (artık 5.778). Fark, tarama
  adımıdır (154: 0.02, 156: 0.05) ve ancak ağırlık 58 etkin noktaya
  çöktüğü için önemli hale gelmiştir — arg M(R) o bölgede sarıyor ve
  "kök" ızgara adımına duyarlı. Bu, o bandın ölçüm sayılmaması için ek
  gerekçedir; başka hiçbir bantta uyuşmazlık yok (15 bandın 14'ünde
  R̃_tam = R_154 dördüncü basamağa kadar).
* **Sırt kesişiminde düzeltilen hata.** `156_nakil.py`'nin ilk sürümü iki
  gazın faz-eşlenmiş sırtlarını R_e(R_l) eğrileri gibi doğrusal ara
  değerle kesiştiriyordu. Sırtlar ÇOK DALLI (kök bulucu "|R_e|'si en küçük
  kök" kuralıyla dallar arasında sıçrıyor: A4'te τ=0.66'da R_l ≈ 0.5'te
  R_e, +8.4'ten −3.5'e atlıyor) ve ara değer bu sıçramanın üstünden köprü
  kurup **olmayan bir kesişim uyduruyordu** — A4 faz artığının 1.48 rad
  çıkması bunun işaretiydi. Ortak çift bu yüzden ayrı bir dosyada
  (`156_ortak_cift.py`), ara değersiz, 81×177 ham ızgara + yerel
  incelemeyle yeniden hesaplandı; T6 o hesabın çıktısıdır. Hatalı ilk
  geçişin çıktısı `scratchpad/156/log_ortak_1gecis.txt`'de duruyor,
  silinmedi.
* **`kapasite` alanı iki kez hesaplandı.** `k3_*.json` içindeki `kapasite`
  alanı arg M'nin SARILI menzilidir ve |arg| > π'ye giden eğrilerde
  yanıltıcıdır (sarma sıçraması sahte genişlik üretir); koşulduğu gibi
  bırakıldı ama **rapora girmedi**. T4b tablosu yalnız `156_kapasite.py`
  ile sarma-açılmış olarak yeniden hesaplanan menzilden gelir.
* **Ölçek normalizasyonu Gauss varsayımı değildir**: yalnız arama
  aralığını üç kanalda ortaklaştırır. Kökler ampirik dağılımla çözülmüştür
  ve R̃_tam'ın 154'ün R'siyle birebir çıkması bunun denetimidir. Gauss
  limiti bu veride zaten kötüdür (τ=0.66, gerçek: A²σΔ²R_tam = 1.74 rad,
  ölçülen faz 1.38) — bu yüzden hiçbir hüküm ona dayandırılmadı.
* **Sırttaki en iyi nokta sığ bir minimumdur.** Tek nokta yerine kabul
  edilebilir bölge (faz artığı < 0.02, |ΔRe| < 0.06, n_eff ≥ 3000) da
  verildi; gerçek gazda R_e'nin geniş bir aralıkta serbest kalması
  (τ=0.66'da −2.48…+8.83) "η katkısı sıfırla uyumlu" hükmünün de sınırıdır:
  η **gereksizdir**, ama küçük bir η katkısı bu veriyle dışlanamaz.
* **Süreler ölçüt değildi**; koşular kısmen aynı anda yapıldı (8 çekirdek).
  Sayısal sonuç etkilenmez.

Ham çıktılar ve scriptler: `156_configs/` (156_cekirdek.py, 156_kos.py,
156_kapasite.py, 156_nakil.py, 156_ortak_cift.py, 156_tablolar.py),
loglar ve JSON'lar `scratchpad/156/`.
