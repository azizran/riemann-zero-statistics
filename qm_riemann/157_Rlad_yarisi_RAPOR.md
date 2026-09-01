# 157 — R_lad(τ): taban-bağımsız nesne sınavı ve kapalı-form yarışı

**Hedef.** 156, Γ_rot'un fazını üreten sağkalım-ağırlık bağlaşımının **merdiven
payında** yaşadığını buldu ve R_lad'ın regresyon tabanına (0.46 / 0.52 / 0.58)
±%3 içinde duyarsız olduğunu ölçtü; buradan "154'ün kapalı-form yarışını
taban-bağımlı R_tam yerine R_lad üstünde tekrarla" önerisi çıktı. Bu rapor o
görevi yerine getirdi — ve önerinin dayanağını çürüttü.

**Figür:** `157_Rlad_yarisi.png` (6 panel).

> ## Kısa hüküm
>
> 1. **R_lad(τ) tam tayfı ölçüldü:** iki gerçek pencere (son-300k L=12.030,
>    orta-300k L=11.464) × **altı taban** (0.28 / 0.34 / 0.40 / 0.46 / 0.52 /
>    0.58) × iki bant ızgarası (0.03 ve 0.02), τ ∈ (0.42, 0.80), apsis τ_eff,
>    negatif dal dahil. 22 koşu. İki pencere sağlam bantlarda %0–3 içinde
>    aynı R_lad'ı veriyor — **R_lad L-değişmez.**
> 2. **τ₀_lad, τ₀_tam'dan FARKLI BİR SAYI DEĞİL — KİMLİKTİR.** R_k = 0 ⇒
>    ağırlık ≡ 1 ⇒ M_k(0) = ⟨e^{−iA·X_tam}⟩ = M_emp, ve M_emp **kanaldan
>    bağımsızdır** (dört kanalda bit düzeyinde özdeş, fark 0.00e+00). Her
>    kanalın R'si tam olarak φ_Γ(τ) = arg M_emp(τ) noktasında sıfırlanır.
>    **Kanal değiştirerek τ₀'ın taban kayması kaldırılamaz** — ve
>    kaldırılmıyor: taban ekseni boyunca yayılım τ₀(φ) 0.0132, τ₀(R_tam)
>    0.0131, **τ₀(R_lad) 0.0140** (son). Eğimler: −0.174 / −0.176 / −0.163
>    (155'in −0.189'u ile uyumlu). Üstelik ikinci bir konvansiyon ekseninde
>    (fit derecesi/aralığı/ızgara, taban sabit) yayılım τ₀(φ) 0.0012,
>    τ₀(R_tam) 0.0038, **τ₀(R_lad) 0.0126** — R_lad'ın sıfırı iki eksende
>    birden en oynak olan. **Konvansiyonsuz τ₀ TANIMLANMADI.**
> 3. **156'nın "R_lad taban-bağımsız" bulgusu bir DURAĞAN NOKTA yanılsaması.**
>    R_lad(taban) eğrisinin minimumu tam 156'nın taradığı {0.46, 0.52, 0.58}
>    penceresinde. Eksen 0.28–0.58'e açıldığında (τ̄=0.625, n_eff ≥ 3000
>    süzgeçli): **R_lad yayılımı %98.7 (son) / %98.8 (orta)**, R_tam yayılımı
>    **%16.9 / %16.9**. 156'nın penceresinde ise sırasıyla %6.1 / %5.5 ve
>    %12.1 / %11.9 — yani dar pencerede sıralama TERS. Bu, 155'in 154'ün
>    "τ₀ tabandan bağımsız" hükmü için bulduğu hatanın birebir aynısıdır.
> 4. **Sebep analitik ve bit düzeyinde ölçüldü.** X_tam = bond(ds) tabandan
>    **bağımsızdır** (ds regresyondan önce tanımlı): arg M_tam(R) ızgarası ve
>    arg M_emp altı tabanda **maks fark 0.000e+00**. Yani R_tam, ilkel φ_Γ'nın
>    SABİT ve tersinir bir yeniden-parametrizasyonudur; bütün konvansiyon
>    bağımlılığı φ'nin kendisindedir. X_lad ise **tanımı gereği** tabana
>    bağlıdır (hangi çizgiler merdivene sayılıyor?): arg M_lad(R̃) ızgarası
>    tabanlar arası **5.3–6.3 rad** oynuyor. **R_lad, φ'nin konvansiyonuna
>    KENDİ konvansiyonunu ekler; azaltmaz.**
> 5. **Yarış koşuldu; kapalı form MÜHÜRLENMEDİ — ama 154'ün beraberliği
>    bozuldu.** Negatif dal (τ < τ₀) fite katıldığında (c) R∞(1−e^{−(τ−τ₀)/w})
>    ile (g) c(τ−τ₀)/τ² arasındaki 154 beraberliği R_lad'da **kesin biçimde
>    (c) lehine** çözülüyor (ΔAIC = 3436 istatistik bütçesinde, 9.8 tam
>    konvansiyon bütçesinde); (g)'nin artığı negatif dalda −0.61'e çıkıyor.
>    R_tam'da ise 154'ün beraberliği aynen sürüyor (ΔAIC = 1.0).
>    **Hangi formun kazandığı KANALA bağlı** — bu, "R'nin kapalı formu"
>    ifadesinin henüz konvansiyondan bağımsız bir anlamı olmadığının ölçüsü.
> 6. **φ-DOĞRUSALLIĞI ÖLDÜ; φ'nin eğriliği MÜHÜRLENDİ.** 155 b ≈ −6.5'i tek
>    pencerede marjinal bulmuştu; burada **10 bağımsız (pencere × taban)
>    koşusunun 9'unda** kuadratik doğrusalı yeniyor (ΔAIC = −7.1 … −134.6;
>    onuncusu, yalnız 5 bantlı orta/0.52, −1.4 ile ayrılamaz),
>    b = **−7.12 ± 0.78**. (h) ailesinin doğrusal üyeleri (h1, h1f, h1½)
>    yarışın dibinde (χ²/dof 227–370, istatistik bütçesi).
> 7. **YENİ ve KONVANSİYONSUZ bir sayı: φ'nin sıfırdaki eğimi.**
>    **a = dφ_Γ/dτ|_{τ₀} = 10.759 ± 0.114** — on koşuda (iki pencere × beş
>    taban) yayılım yalnız **%3.7**, oysa aynı koşularda τ₀ 0.0131 kayıyor.
>    Taban konvansiyonu ilkel φ_Γ'yı pratikte **SALT ÖTELİYOR**: a·Δτ₀ = 0.141
>    rad, sabit τ̄=0.53'te ölçülen taban aralığı 0.132 / 0.134 rad.
>    Konvansiyondan bağımsız içerik φ'nin **ŞEKLİ**, τ₀'ın **YERİ** değil.
> 8. **τ₀ ≡ ½ dört hata modelinin üçünde ölü.** (g½) ve (h1½) istatistik
>    bütçesinde χ²/dof = 202 / 226; yalnız tam konvansiyon bütçesinde
>    (⟨σ⟩ = 0.35, her şeyin uyduğu bütçe) hayatta. 154'ün "τ₀ ≠ ½" hükmü
>    R_lad'da da aynı statüde duruyor: **konvansiyon bütçesine bakılırsa
>    karar verilemez** (155'in R1/R2 açmazı aynen geçerli).

---

## 1. Yöntem — üç raporun zinciri, tek kod yolu

`157_configs/157_cekirdek.py` hiçbir ölçüm parçasını kopyalamaz; **import
eder**:

| kaynak | import edilen | ne |
|---|---|---|
| `154_cekirdek` | `pk`, `pk_m`, `_M`, `_R_150grid` | Γ/M çekirdeği ve R kök mantığı |
| `156_cekirdek` | `zincir3`, `_bond` | ds = drift + lad + eta ayrışımı, bond adımı |

Kendi yazdığı tek şey **bant döngüsü**dür; o da 155'in çizgi-bazında kayıt
(F1) + **τ_eff apsisi** (F2) eklentisiyle. 157'nin kendi eklentileri:

| | ekleme | niçin |
|---|---|---|
| G1 | apsis her yerde **τ_eff** = Σ τ_i(pow_on−pow_off)_i / Σ(pow_on−pow_off)_i | 155, bant-ortası apsisin τ₀'a SAHTE bir L-bağımlılığı imal ettiğini ölçmüştü |
| G2 | **arg M_k(R) ızgarası JSON'a yazılır** (R̃ ∈ [−8, 10], adım 0.05) | (h) "φ-ilkel" adayı — φ modellenip R'nin ondan TÜRETİLMESİ — ölçüm yeniden koşulmadan, aynı ampirik ağırlık dağılımıyla sınanabilsin |
| G3 | jackknife φ, τ_eff, ρ, R_tam, R_lad, R_eta için **eşzamanlı** (aynı 8 grup, bütün bantlarda aynı k indisi) | τ₀ fiti replika-replika tekrarlanabilsin (155'in hata mimarisi) |
| G4 | **sarma (phase-wrap) bayrağı**: R=0'dan köke giden ızgara parçasında ardışık \|Δ arg M\| > 1 rad ise kök sahtedir | 154'ün Denetim 3'ünün otomatikleştirilmesi |

Ölçek normalizasyonu 156'nınkinin aynısıdır: s_k = σΔ²/Cov(X_k, X_tam),
Y_k = X_k·s_k, R_k(ham) = R̃_k·s_k. **Ham R doğrudan X_k'nin ağırlıktaki
katsayısıdır** (w = e^{−A·R_ham·X_k}); tablolar ham R'dir.

### Kopya-kayması denetimi — dört düzeyde geçti

`157_configs/157_dogrulama.py` çıktısı (`scratchpad/157/dogrulama.txt`):

| | sınav | sonuç |
|---|---|---|
| **V1** | 156.zincir3'ün η'sı = 154.eta_zinciri'nin (155'in önbelleklediği) η'sı? | 8 (pencere,taban) çiftinde **maks\|Δη\| = 3.3e−16**, Δσ_η² = 0, Δc₁ ≤ 3e−18; ds = drift+lad+eta kapanışı ≤ 4.4e−16 |
| **V2/V3** | 157'nin Γ, φ_Γ, τ_eff, ρ ve R_tam'ı = 155'in kayıtlı JSON'ları? | beş koşuda 5–8 ortak bant: maks Δ\|Γ\| ≤ 3.3e−16, Δφ ≤ 1.1e−16, Δτ_eff ≤ 1.7e−16, Δρ ≤ 2.2e−16, **ΔR_tam ≤ 1.3e−15** |
| **V5** | M_k(R=0) kanaldan bağımsız mı? | τ = 0.47 / 0.51 / 0.55 / 0.66'da dört kanalda (tam, lad, eta, drift) **maks fark = 0.00e+00** |
| **V6** | arg M_tam(R) ve arg M_emp tabandan bağımsız mı? | 0.28↔0.52, …, 0.58↔0.52: **maks\|Δ arg M_tam(R)\| = 0.000e+00**, **maks\|Δ arg M_emp\| = 0.000e+00**; buna karşılık **maks\|Δ arg M_lad(R̃)\| = 5.30 … 6.28 rad** |
| **V7** | (h) adayı için JSON ızgarasının tersi = koşudaki bisection kökü? | maks\|Δ R̃\| = 8e−5 (tam) / 2.2e−4 (lad) — σ_jk ≈ 5e−3 yanında ihmal edilebilir |

**V5 ve V6 bu raporun analitik omurgasıdır** ve aşağıda §3–§4'te kullanılıyor.

### Koşulan ölçümler

| pencere | tabanlar | ızgaralar | bant |
|---|---|---|---|
| son (Z[−300000:], t ∈ [975 796, 1 132 491], L = 12.0296) | 0.28, 0.34, 0.40, 0.46, 0.52, 0.58 | `kaba` = 0.03, τ ∈ (0.43, 0.79]; `ince` = 0.02, τ ∈ (0.42, 0.80]; `tau0` = 0.02, τ ∈ (0.42, 0.62] | taban altı bantlar elenir |
| orta (Z[850000:1150000], L = 11.4638) | aynı | aynı | aynı |

`kaba` ızgarası 0.43'e demirlidir: kenarları **hem 0.46'ya hem 0.52'ye
düşer**, yani taban değişince ızgara kaymaz, yalnız alt bantlar elenir
(taban 0.28/0.34/0.40 → 12 bant, 0.46 → 11, 0.52 → 9, 0.58 → 7). Toplam 22
koşu; cap = 4000 (155 ile aynı). Bütün tablolar
`scratchpad/157/tablolar.md` ve `analiz_cikti.txt`'ten otomatik üretildi.

**Sağlamlık kuralları (156/150'den devralındı, her tabloda uygulanır):** faz
artığı < 0.02 rad · sarma yok · \|Γ\| ≤ 1.5 · **n_eff ≥ 3000** (ağırlığın
%1'i). Aşırı-belirleme eşiği \|ΔRe\| < 0.06 ayrıca raporlanır ama bant
elemesi için yalnız bir robustluk koşusunda kullanılır (§5).

---

## 2. T0 — kanal künyeleri ve tabanın ne yaptığı

| pencere | taban | regresör | σ_ds² | σ_lad² | σ_η² | c₁(η) | **σΔ² = Var(X_tam)** | pay Cov(X_lad,X_tam)/σΔ² | korel(X_lad,X_tam) | s_lad |
|---|---|---|---|---|---|---|---|---|---|---|
| son | 0.28 | 16 | 0.16744 | 0.06588 | 0.10156 | −0.06900 | **0.05466** | +0.7457 | +0.8394 | 1.3409 |
| son | 0.34 | 25 | 0.16744 | 0.08982 | 0.07762 | −0.05659 | **0.05466** | +0.8652 | +0.9007 | 1.1559 |
| son | 0.40 | 41 | 0.16744 | 0.11312 | 0.05431 | −0.03933 | **0.05466** | +0.9279 | +0.9312 | 1.0777 |
| son | 0.46 | 69 | 0.16744 | 0.13186 | 0.03558 | −0.02299 | **0.05466** | +0.9480 | +0.9427 | 1.0549 |
| son | 0.52 | 117 | 0.16744 | 0.14469 | 0.02274 | −0.01158 | **0.05466** | +0.9533 | +0.9492 | 1.0489 |
| son | 0.58 | 206 | 0.16744 | 0.15282 | 0.01461 | −0.00514 | **0.05466** | +0.9595 | +0.9568 | 1.0422 |
| orta | 0.28 | 13 | 0.16659 | 0.06444 | 0.10215 | −0.06931 | **0.05431** | +0.7419 | +0.8367 | 1.3478 |
| orta | 0.40 | 35 | 0.16659 | 0.11209 | 0.05451 | −0.03963 | **0.05431** | +0.9275 | +0.9312 | 1.0782 |
| orta | 0.46 | 58 | 0.16659 | 0.13148 | 0.03511 | −0.02276 | **0.05431** | +0.9480 | +0.9434 | 1.0548 |
| orta | 0.52 | 95 | 0.16659 | 0.14408 | 0.02251 | −0.01154 | **0.05431** | +0.9532 | +0.9497 | 1.0491 |
| orta | 0.58 | 159 | 0.16659 | 0.15209 | 0.01451 | −0.00515 | **0.05431** | +0.9593 | +0.9570 | 1.0424 |

(0.34 satırları tam tabloda; `tablolar.md`. Taban 0.52 sütunu 156'nın T0'ıyla
birebir: σ_ds² 0.16744, σ_lad² 0.14469, σ_η² 0.02274, c₁ −0.01158, pay
+0.9533, s_lad 1.049.)

**Kalın sütun raporun anahtarıdır.** σΔ² = Var(X_tam) altı tabanda **aynıdır**,
çünkü ds regresyondan önce tanımlıdır ve taban yalnız onun *ayrışımını*
değiştirir. σ_lad² ise 0.0659 → 0.1528 (2.3 kat) değişiyor: **merdiven kanalı
tabanın kendisidir**, bağımsız bir fiziksel nesne değil.

---

## 3. T1 — R tayfı (taban 0.46, 0.03'lük ızgara, apsis τ_eff)

| τ̄ | τ_eff | pencere | \|Γ\| | φ_Γ | arg M_emp | R_tam | ±jk | \|ΔRe\|_tam | **R_lad** | ±jk | \|ΔRe\|_lad | n_eff(lad) | ρ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.4750 | 0.4746 | son | 0.889 | −0.4095 | +0.0163 | −0.807 | 0.012 | 0.1236 | **−1.016** | 0.017 | 0.0960 | 207 535 | 0.6039 |
| 0.4750 | 0.4753 | orta | 0.892 | −0.4019 | +0.0164 | −0.797 | 0.045 | 0.1256 | **−0.999** | 0.065 | 0.0984 | 209 965 | 0.5963 |
| 0.5050 | 0.5042 | son | 0.873 | −0.0864 | +0.0197 | −0.189 | 0.014 | 0.1215 | **−0.207** | 0.016 | 0.1122 | 293 538 | 0.5558 |
| 0.5050 | 0.5045 | orta | 0.869 | −0.0857 | +0.0198 | −0.189 | 0.024 | 0.1165 | **−0.207** | 0.027 | 0.1073 | 293 580 | 0.5630 |
| 0.5350 | 0.5344 | son | 0.845 | +0.2315 | +0.0234 | +0.352 | 0.010 | 0.0897 | **+0.336** | 0.008 | 0.1080 | 276 176 | 0.5174 |
| 0.5350 | 0.5341 | orta | 0.847 | +0.2276 | +0.0236 | +0.348 | 0.015 | 0.0899 | **+0.332** | 0.013 | 0.1077 | 277 022 | 0.5147 |
| 0.5650 | 0.5651 | son | 0.811 | +0.5465 | +0.0277 | +0.840 | 0.005 | 0.0467 | **+0.710** | 0.004 | 0.0831 | 183 853 | 0.4674 |
| 0.5650 | 0.5647 | orta | 0.814 | +0.5406 | +0.0279 | +0.837 | 0.010 | 0.0473 | **+0.710** | 0.007 | 0.0832 | 185 594 | 0.4648 |
| 0.5950 | 0.5944 | son | 0.776 | +0.8325 | +0.0324 | +1.235 | 0.011 | 0.0122 | **+0.957** | 0.007 | 0.0499 | 103 818 | 0.4122 |
| 0.5950 | 0.5945 | orta | 0.775 | +0.8302 | +0.0327 | +1.244 | 0.011 | 0.0094 | **+0.966** | 0.007 | 0.0469 | 103 835 | 0.4147 |
| 0.6250 | 0.6238 | son | 0.725 | +1.0957 | +0.0377 | +1.555 | 0.012 | 0.0154 | **+1.131** | 0.007 | 0.0124 | 56 781 | 0.3494 |
| 0.6250 | 0.6243 | orta | 0.721 | +1.1032 | +0.0381 | +1.589 | 0.011 | 0.0189 | **+1.154** | 0.006 | 0.0082 | 54 916 | 0.3461 |
| 0.6550 | 0.6539 | son | 0.680 | +1.3465 | +0.0436 | +1.823 | 0.011 | 0.0172 | **+1.265** | 0.006 | 0.0041 | 31 109 | 0.3153 |
| 0.6550 | 0.6543 | orta | 0.678 | +1.3463 | +0.0440 | +1.848 | 0.007 | 0.0185 | **+1.283** | 0.004 | 0.0054 | 30 784 | 0.3057 |
| 0.6850 | 0.6836 | son | 0.644 | +1.5499 | +0.0500 | +1.982 | 0.012 | 0.0023 | **+1.337** | 0.007 | 0.0011 | 19 778 | 0.2720 |
| 0.6850 | 0.6843 | orta | 0.640 | +1.5710 | +0.0505 | +2.049 | 0.015 | 0.0000 | **+1.378** | 0.008 | 0.0000 | 18 234 | 0.2714 |
| 0.7150 | 0.7144 | son | 0.595 | +1.7157 | +0.0571 | +2.057 | 0.013 | 0.0215 | **+1.365** | 0.007 | 0.0141 | 14 415 | 0.2200 |
| 0.7150 | 0.7144 | orta | 0.585 | +1.7070 | +0.0577 | +2.074 | 0.017 | 0.0219 | **+1.377** | 0.009 | 0.0148 | 14 812 | 0.2143 |
| 0.7450 | 0.7445 | son | 0.580 | +1.8304 | +0.0648 | +2.038 | 0.029 | 0.0385 | **+1.342** | 0.016 | 0.0255 | 12 724 | 0.1819 |
| 0.7450 | 0.7442 | orta | 0.562 | +1.8084 | +0.0655 | +2.034 | 0.008 | 0.0399 | **+1.343** | 0.005 | 0.0276 | 13 697 | 0.1785 |
| 0.7750 ‡ | 0.7731 | son | 0.615 | +1.7933 | +0.0732 | +1.796 | 0.037 | 0.0184 | **+1.195** | 0.021 | 0.0049 | 18 730 | 0.1401 |
| 0.7750 ‡ | 0.7736 | orta | 0.603 | +1.7634 | +0.0741 | +1.781 | 0.016 | 0.0186 | **+1.190** | 0.009 | 0.0066 | 20 339 | 0.1377 |

‡ = τ > 0.76, 154'ün çözünmeyen-çizgi bölgesi; hiçbir fite girmedi.

**L-değişmezlik.** İki pencere sağlam bantlarda %0–3 içinde aynı: R_lad =
0.710/0.710 (τ̄=0.565, %0.0), 1.131/1.154 (0.625, %2.0), 1.265/1.283 (0.655,
%1.4), 1.337/1.378 (0.685, %3.0), 1.365/1.377 (0.715, %0.9), 1.342/1.343
(0.745, %0.1). ΔL = 0.566 ve çizgi kümeleri tümüyle farklı. **R_lad(τ) gerçek
bir τ fonksiyonudur** — konvansiyona asılı olması onu gürültü yapmıyor,
yalnız konvansiyonsuz yapmıyor.

**Negatif dal ölçüldü ve sağlam.** τ̄ = 0.475 ve 0.505'te R_lad = −1.016 ve
−0.207; faz artığı 0, ağırlık çökmemiş (n_eff = 2.1e5 / 2.9e5). Taban 0.40 ve
0.34/0.28 ile dal τ̄ = 0.445'e kadar açılıyor (R_lad = −2.207 / −2.092 /
−2.045). 155'in "derin negatif dalda R tanımsızlaşır" uyarısı bu aralıkta
devreye girmiyor.

**Dürüst kayıt — aşırı-belirleme τ < 0.58'de KAYBEDİLİYOR.** \|ΔRe\| eşiği
(0.06) τ̄ ≤ 0.565 bantlarında **her iki kanalda da** aşılıyor (lad: 0.083–0.112;
tam: 0.047–0.126). Yani sıfır bölgesinde ve negatif dalda tek-kanallı
sağkalım-ağırlık modeli bağımsız sınavı geçmiyor; oradaki R bir **faz-eşleme
sayısıdır**, doğrulanmış bir model değil. Bu, yarışın en ayırt edici
noktalarının (negatif dal) aynı zamanda modelin en zayıf olduğu yer olması
demektir ve §5'te bir robustluk koşusuyla ayrıca sınandı.

---

## 4. τ₀ — kimlik, taban yasası ve konvansiyon duyarlılığı

### 4a. KİMLİK: bütün kanalların R'si AYNI τ'da sıfırlanır

R_k = 0 ⇒ w_k ≡ 1 ⇒ M_k(0) = ⟨e^{−iA·X_tam}⟩ = M_emp. **M_emp kanaldan
bağımsızdır** (V5: dört kanalda bit düzeyinde özdeş). Dolayısıyla

> **τ₀ (her kanal) = φ_Γ(τ) = arg M_emp(τ) denkleminin kökü.**

Bu, görevin "en önemli soru"suna doğrudan yanıttır: **R_tam'daki taban kayması
R_lad'a geçerken kaybolamaz, çünkü sıfır-geçişi kanal değiştirmiyor.**
Aşağıda hem φ'nin sıfırı, hem bu **konvansiyonsuz kök τ₀\*** (= φ − arg M_emp
sıfırı), hem de iki kanalın kuadratik-fit kökleri veriliyor.

| pencere | taban | σ_η² | τ₀(φ_Γ) | **τ₀\*** = (φ−argM_emp) sıfırı | τ₀(R_tam) | τ₀(R_lad) | Δ(lad−tam) | 155'in τ₀(φ)'si |
|---|---|---|---|---|---|---|---|---|
| son | 0.28 | 0.1016 | 0.4998 ± 0.0001 | **0.5015 ± 0.0001** | 0.5013 ± 0.0001 | 0.5006 ± 0.0002 | −0.0007 | 0.4998 |
| son | 0.34 | 0.0776 | 0.5045 ± 0.0001 | **0.5063 ± 0.0001** | 0.5060 ± 0.0001 | 0.5044 ± 0.0002 | −0.0015 | 0.5044 |
| son | 0.40 | 0.0543 | 0.5092 ± 0.0001 | **0.5111 ± 0.0001** | 0.5105 ± 0.0001 | 0.5090 ± 0.0002 | −0.0015 | 0.5091 |
| son | 0.46 | 0.0356 | 0.5123 ± 0.0001 | **0.5143 ± 0.0001** | 0.5143 ± 0.0001 | 0.5146 ± 0.0001 | +0.0003 | 0.5122 |
| son | 0.52 | 0.0227 | 0.5129 ± 0.0005 | **0.5148 ± 0.0004** | 0.5144 ± 0.0004 | 0.5112 ± 0.0004 | −0.0032 | 0.5126 |
| orta | 0.28 | 0.1022 | 0.4997 ± 0.0002 | **0.5015 ± 0.0002** | 0.5012 ± 0.0002 | 0.5006 ± 0.0003 | −0.0007 | 0.4996 |
| orta | 0.34 | 0.0758 | 0.5050 ± 0.0001 | **0.5068 ± 0.0001** | 0.5064 ± 0.0001 | 0.5048 ± 0.0003 | −0.0017 | 0.5049 |
| orta | 0.40 | 0.0545 | 0.5093 ± 0.0001 | **0.5112 ± 0.0001** | 0.5106 ± 0.0001 | 0.5090 ± 0.0001 | −0.0016 | 0.5091 |
| orta | 0.46 | 0.0351 | 0.5124 ± 0.0001 | **0.5144 ± 0.0001** | 0.5144 ± 0.0001 | 0.5147 ± 0.0001 | +0.0003 | 0.5122 |
| orta | 0.52 | 0.0225 | 0.5123 ± 0.0004 | **0.5142 ± 0.0003** | 0.5138 ± 0.0003 | 0.5109 ± 0.0002 | −0.0029 | 0.5125 |

**155'in τ₀(φ) tablosu birebir yeniden üretildi** (10 satırın 10'unda fark
≤ 0.0003). Δ(lad−tam) sütununun sıfır olmaması kimliğin ihlali **değildir**:
kimlik *tam* sıfır-geçişi içindir, tablodaki kanal sütunları ise iki farklı
R(τ) eğrisine atılan **kuadratik fitin** kökleridir; iki eğrinin eğriliği
farklı olduğu için fit kökleri ~0.003 ayrışır. Kimliğin kendisi τ₀\* sütunudur.

### 4b. Taban yasası — kayma kaybolmuyor, hafifçe BÜYÜYOR

| pencere | nicelik | dτ₀/dσ_η² | kesişim (σ_η² → 0) | yayılım (taban 0.28↔0.52) |
|---|---|---|---|---|
| son | τ₀(φ) | −0.174 | 0.5179 | 0.0132 |
| son | **τ₀\*** | −0.175 | 0.5198 | 0.0132 |
| son | τ₀(R_tam) | −0.176 | 0.5196 | 0.0131 |
| son | **τ₀(R_lad)** | −0.163 | 0.5175 | **0.0140** |
| orta | τ₀(φ) | −0.167 | 0.5174 | 0.0127 |
| orta | **τ₀\*** | −0.169 | 0.5194 | 0.0129 |
| orta | τ₀(R_tam) | −0.169 | 0.5191 | 0.0132 |
| orta | **τ₀(R_lad)** | −0.159 | 0.5172 | **0.0141** |

155 (üç pencere, aralık (0.46,0.56]): eğim −0.189 / −0.189 / −0.195, kesişim
0.5191 / 0.5191 / 0.5196. Bu koşu (aralık τ̄ ≤ 0.57, yani 0.42'den itibaren)
−0.174 / −0.167 ve 0.5198 / 0.5194 veriyor — aynı yasa, aralık farkı kadar
kayık. **Hüküm: R_lad'a geçmek 0.19σ_η² kaymasını kaldırmıyor; R_lad'ın
yayılımı R_tam'ınkinden %7 DAHA BÜYÜK.**

### 4c. İkinci konvansiyon ekseni: fit derecesi / aralığı / ızgara

Taban **sabit 0.46**, yalnız fit konvansiyonu değişiyor (ızgara 0.02 ↔ 0.03,
aralık τ̄ ≤ 0.55 / 0.57 / 0.60 / 0.63, derece 1 ↔ 2):

n = 32 konvansiyon (2 pencere × 2 ızgara × 4 aralık × 2 derece, tanımsız
kombinasyonlar hariç):

| nicelik | ortalama ± sd | yayılım |
|---|---|---|
| τ₀(φ_Γ) | 0.5126 ± 0.0003 | **0.0012** |
| **τ₀\*** = (φ−argM_emp) sıfırı | 0.5146 ± 0.0004 | **0.0014** |
| τ₀(R_tam) | 0.5157 ± 0.0015 | 0.0038 |
| **τ₀(R_lad)** | 0.5191 ± 0.0045 | **0.0126** |

**R_lad'ın sıfırı bu eksende de en oynak olandır** — konvansiyon yayılımı
(0.0126) neredeyse taban yayılımı (0.0140) kadar büyük, yani R_lad'ın τ₀'ı
**iki ayrı eksende birden** belirsiz. En dayanıklısı τ₀\* (0.0014), sonra
φ'nin kendi sıfırı (0.0012). Ayrıntı (32 satır) `analiz_cikti.txt` B2
bölümünde.

Ayrıca: **154/155'in 0.5153'ü bu taramanın içinde yer alıyor.** τ₀(R_tam)
doğrusal fitle 0.5154–0.5180, kuadratik fitle 0.5142–0.5147; tarama ortalaması
**0.5157 ± 0.0015**. Yani 0.5153 ölçüm olarak doğrudur; "sabit" olarak
alıntılanması bu ±0.0038'lik konvansiyon yayılımını gizler.

---

## 5. Taban çöküşü — 156'nın ±%3'ü nereden geldi

### 5a. R(taban) sabit τ'da (pencere son, 0.03'lük ızgara)

**R_tam (ham)**

| taban | σ_η² | σ_lad² | τ̄=0.595 | 0.625 | 0.655 | 0.685 | 0.715 | 0.745 |
|---|---|---|---|---|---|---|---|---|
| 0.28 | 0.1016 | 0.0659 | +1.495 | +1.839 | +2.123 | +2.371 | +2.513 | — |
| 0.34 | 0.0776 | 0.0898 | +1.382 | +1.706 | +1.973 | +2.175 | +2.289 | +2.254 |
| 0.40 | 0.0543 | 0.1131 | +1.282 | +1.606 | +1.850 | +1.996 | +2.099 | +2.114 |
| 0.46 | 0.0356 | 0.1319 | +1.235 | +1.555 | +1.823 | +1.982 | +2.057 | +2.038 |
| 0.52 | 0.0227 | 0.1447 | +1.276 | +1.611 | +1.905 | +2.095 | +2.216 | +2.169 |
| 0.58 | 0.0146 | 0.1528 | +1.381 | +1.753 | +2.074 | +2.327 | +2.535 | — |
| **yayılım% — 156 penceresi {0.46,0.52,0.58}** | | | 11.2 | 12.1 | 13.0 | 16.2 | 21.1 | 6.2 |
| **yayılım% — tam eksen {0.28…0.58}** | | | **19.4** | **16.9** | **15.3** | **18.0** | **20.9** | 10.1 |

**R_lad (ham)**

| taban | σ_η² | σ_lad² | τ̄=0.595 | 0.625 | 0.655 | 0.685 | 0.715 | 0.745 |
|---|---|---|---|---|---|---|---|---|
| 0.28 | 0.1016 | 0.0659 | +1.994 | +2.602 | +3.260 | +4.071 ⚠ | +4.818 ⚠ | +5.960 ⚠ |
| 0.34 | 0.0776 | 0.0898 | +1.351 | +1.660 | +1.937 | +2.172 | +2.329 ⚠ | +2.310 ⚠ |
| 0.40 | 0.0543 | 0.1131 | +1.063 | +1.279 | +1.440 | +1.536 | +1.607 | +1.616 |
| 0.46 | 0.0356 | 0.1319 | +0.957 | +1.131 | +1.265 | +1.337 | +1.365 | +1.342 |
| 0.52 | 0.0227 | 0.1447 | +0.963 | +1.123 | +1.246 | +1.315 | +1.351 | +1.312 |
| 0.58 | 0.0146 | 0.1528 | +1.030 | +1.194 | +1.313 | +1.395 | +1.459 | +1.438 |
| **yayılım% — 156 penceresi {0.46,0.52,0.58}** | | | 7.4 | 6.1 | 5.2 | 6.0 | 7.7 | 9.2 |
| **yayılım% — tam eksen {0.28…0.58}** | | | **84.6** | **98.7** | **115.5** | 55.3 (5) | 17.7 (4) | 21.3 (4) |

⚠ = n_eff < 3000; yayılım satırlarında sayılmadı (parantez: kaç taban kaldı).
`orta` penceresi aynı örüntüyü veriyor: R_lad yayılımı 85.8 / 98.8 / 46.4 %
(tam eksen) vs 7.7 / 5.5 / 5.0 % (156 penceresi); R_tam 20.9 / 16.9 / 16.3 %
vs 11.8 / 11.9 / 14.3 %.

**Hüküm.** R_lad(taban) eğrisinin **minimumu tam 156'nın taradığı aralıkta**
(0.46–0.52). Orada eğri düz olduğu için ±%3–7 çıkıyor; eksen açıldığında
%85–116'ya fırlıyor. R_tam ise bütün eksende %15–21 ile **beş kat daha
kararlı**. **156'nın "R_lad taban-bağımsız" hükmü, dar bir taban penceresinin
yerel düzlüğünün yanlış okunmasıdır** — 155'in 154'ün "τ₀ tabandan bağımsız"
hükmünde bulduğu hatanın tam aynısı.

### 5b. Kaydırılmış apsis eğrileri çökertmiyor

Apsis (τ_eff − τ₀(taban)) yapılırsa tabanlar arası yayılım:

| pencere | kanal | ham apsis τ_eff | kaydırılmış τ_eff − τ₀ |
|---|---|---|---|
| son | tam | ⟨yayılım⟩ = 0.332 (17.6%) | 0.244 (12.3%) |
| son | **lad** | ⟨yayılım⟩ = 2.529 (**125.2%**) | 2.151 (**111.6%**) |
| orta | tam | 0.351 (18.2%) | 0.263 (13.0%) |
| orta | **lad** | 2.651 (127.8%) | 2.283 (114.6%) |

Yani R_lad'ın taban bağımlılığı **salt bir τ₀ ötelemesi değil**; ölçek/şekil
bağımlılığıdır ve kaydırmayla kurtarılmıyor.

### 5c. Neden — ve tersine, İLKEL nesnede işe yarayan çöküş

V6 (§1) bunun sebebini bit düzeyinde veriyor: X_tam tabandan bağımsız
olduğundan **arg M_tam(R) eğrisi de M_emp de altı tabanda özdeştir** (fark
0.000e+00). Dolayısıyla R_tam = argM_tam⁻¹(φ_Γ) — ilkel φ_Γ'nın sabit,
tersinir bir yeniden-parametrizasyonu; **bütün konvansiyon bağımlılığı φ'de.**
X_lad ise tabanın tanımladığı bir değişkendir; arg M_lad(R̃) tabanlar arası
5.3–6.3 rad oynuyor. **R_lad = argM_lad(taban)⁻¹(φ_Γ(taban))** — iki
konvansiyon kaynağı.

Buna karşılık **ilkel φ_Γ'da çöküş ÇALIŞIYOR.** φ_Γ(τ) = a(τ−τ₀) + b(τ−τ₀)²
fiti, on (pencere × taban) koşusunda:

| pencere | taban | τ₀(φ) | a = dφ/dτ\|_τ₀ | b |
|---|---|---|---|---|
| son | 0.28 | 0.4999 | 10.790 | −5.890 |
| son | 0.34 | 0.5043 | 10.689 | −6.877 |
| son | 0.40 | 0.5088 | 10.615 | −7.759 |
| son | 0.46 | 0.5120 | 10.724 | −7.502 |
| son | 0.52 | 0.5127 | 11.013 | −7.272 |
| orta | 0.28 | 0.4995 | 10.814 | −6.495 |
| orta | 0.34 | 0.5047 | 10.695 | −7.123 |
| orta | 0.40 | 0.5087 | 10.656 | −7.801 |
| orta | 0.46 | 0.5121 | 10.753 | −8.338 |
| orta | 0.52 | 0.5124 | 10.845 | −6.152 |

> **a = 10.759 ± 0.114** (yayılım %3.7) · **b = −7.121 ± 0.776** (yayılım
> %34.4) · **τ₀ = 0.5075 ± 0.0051** (yayılım 0.0131)

Denetim: a·Δτ₀ = **0.141 rad**; sabit τ̄ = 0.53'te φ_Γ'nın taban aralığı
**0.132 rad** (son) / **0.134 rad** (orta). **Taban konvansiyonu ilkel φ_Γ'yı
%7 içinde SALT ÖTELİYOR.** Konvansiyondan bağımsız içerik φ'nin *şeklidir*
(a, ve daha zayıf olarak b), τ₀'ın *yeri* değil.

---

## 6. İlkel gözlenebilir φ_Γ: doğrusal mı, eğri mi? (155'in b ≈ −6.5 sınavı)

τ̄ ≤ 0.61, apsis τ_eff, ağırlık = jackknife σ_φ:

| pencere | taban | bant | doğrusal τ₀ | eğim | χ²/dof | kuadratik τ₀ | a | b | χ²/dof | ΔAIC | hüküm |
|---|---|---|---|---|---|---|---|---|---|---|---|
| son | 0.28 | 10 | 0.5001 | 10.326 | 7.76 | 0.4999 | 10.790 | −5.890 | 0.11 | −59.4 | **EĞRİ** |
| son | 0.34 | 10 | 0.5044 | 10.109 | 14.93 | 0.5043 | 10.689 | −6.877 | 0.24 | −115.8 | **EĞRİ** |
| son | 0.40 | 10 | 0.5093 | 10.059 | 17.92 | 0.5088 | 10.615 | −7.759 | 0.96 | −134.6 | **EĞRİ** |
| son | 0.46 | 8 | 0.5118 | 10.168 | 8.03 | 0.5120 | 10.724 | −7.502 | 0.83 | −42.0 | **EĞRİ** |
| son | 0.52 | 5 | 0.5107 | 10.201 | 3.06 | 0.5127 | 11.013 | −7.272 | 0.03 | −7.1 | **EĞRİ** |
| orta | 0.28 | 10 | 0.4996 | 10.164 | 5.13 | 0.4995 | 10.814 | −6.495 | 0.21 | −37.6 | **EĞRİ** |
| orta | 0.34 | 10 | 0.5049 | 10.066 | 7.02 | 0.5047 | 10.695 | −7.123 | 0.08 | −53.6 | **EĞRİ** |
| orta | 0.40 | 10 | 0.5090 | 10.041 | 8.77 | 0.5087 | 10.656 | −7.801 | 0.29 | −66.2 | **EĞRİ** |
| orta | 0.46 | 8 | 0.5116 | 10.050 | 5.62 | 0.5121 | 10.753 | −8.338 | 0.24 | −30.5 | **EĞRİ** |
| orta | 0.52 | 5 | 0.5102 | 10.052 | 1.16 | 0.5124 | 10.845 | −6.152 | 0.05 | −1.4 | ayrılamaz |

**Not:** buradaki `a`, φ'nin τ₀'daki **türevidir**; 155'in kuadratik
tablosundaki "a" ham polinom katsayısı c₁'dir — iki sayı doğrudan
karşılaştırılmaz. 155'in DOĞRUSAL a'sı (10.686 / 10.687 / 10.665, taban 0.40)
ise bu tablonun doğrusal eğim sütunuyla (10.059 / 10.041) aynı mertebede;
kalan fark fit aralığındandır (155: τ ∈ (0.42,0.56], burada τ̄ ≤ 0.61).

**Hüküm:** 155'in "tek pencerede marjinal anlamlı" eğriliği **mühürlendi**.
On bağımsız (pencere × taban) koşusunun **dokuzunda** kuadratik doğrusalı
ΔAIC = −7.1 … −134.6 ile yeniyor; kalan biri (orta/0.52, yalnız 5 bant)
ayrılamaz diyor. b = −5.89 … −8.34, ortalama **−7.12 ± 0.78**.
Kaynağı hâlâ bilinmiyor (155'in Açık 5'i açık kalıyor).

---

## 7. KAPALI-FORM YARIŞI

**Kural (155/154 disiplini):** apsis τ_eff · iki pencere **ortak** fit ·
τ > 0.76 elenir · faz artığı ≥ 0.02, sarma, \|Γ\| > 1.5 ve n_eff < 3000
elenir · skor **AIC = χ² + 2k**.

**Adaylar** (≤ 2 parametreli olanlar kalın; 3 parametreliler bütçe dışıdır,
karşılaştırma için koşuldu):

| ad | tanım | k |
|---|---|---|
| **a** | R = c (sabit) | 1 |
| **b** | R = c(τ−τ₀) | 2 |
| **c** | R = R∞(1−e^{−(τ−τ₀)/w}), **τ₀ ölçülen değere sabit** | 2 |
| c3 | aynı, τ₀ serbest | 3 |
| **g** | R = c(τ−τ₀)/τ² | 2 |
| **g½** | R = c(τ−½)/τ² (τ₀ ≡ ½) | 1 |
| **e** | R = k·A²σΔ² | 1 |
| **h1** | **(h) φ-İLKEL:** φ = a(τ−τ₀), R = argM_lad⁻¹(φ) | 2 |
| **h1f** | aynı, τ₀ ölçülen değere sabit | 1 |
| **h1½** | aynı, τ₀ ≡ ½ | 1 |
| **h2b** | φ = a(τ−τ₀) − 6.5(τ−τ₀)² (b = 155'in değeri, SABİT) | 2 |
| h2 | φ = a(τ−τ₀) + b(τ−τ₀)² | 3 |

(h) ailesinde R **hesaplanmaz, TÜRETİLİR**: model φ'yi verir, R o bandın
**ölçülmüş ampirik** arg M_lad(R̃) ızgarası ters çevrilerek bulunur (V7:
ızgara tersi bisection kökünü 2e−4 içinde geri veriyor). Yani (h), "R'nin
bütün şekli tahmincinin kinematiğidir, fizik yalnız φ'dedir" hipotezidir.

**Hata bütçeleri (dört katman, 154'ün usulü).** Her katman bir öncekine bir
sistematik ekler; sıralamanın hangi katmanda çöktüğü hükmün kendisidir.

| model | bileşenler | ⟨σ⟩ (R_lad, taban 0.46) |
|---|---|---|
| **I** | σ_jk (yalnız istatistik) | 0.0123 |
| **II** | ⊕ σ_pencere = \|R_son − R_orta\|/2 | 0.0155 |
| **III** | ⊕ σ_taban, **156'nın penceresi** {0.46, 0.52, 0.58} yarı-menzili | 0.0347 |
| **IV** | ⊕ σ_taban, **tam eksen** {0.28…0.58} yarı-menzili | 0.3504 |

### 7a. BİRİNCİL: R_lad, taban 0.46, 0.03'lük ızgara, n = 20

AIC (parantezde χ²/dof):

| aday | k | **I** istatistik | **II** +pencere | **III** +156 tabanı | **IV** +tam taban |
|---|---|---|---|---|---|
| **c** | 2 | 380.0 (20.9) | 154.3 (8.35) | **11.7 (0.43)** | **4.27 (0.02)** |
| h2 | 3 | **276.4 (15.9)** | **130.0 (7.29)** | 19.0 (0.77) | 6.27 (0.02) |
| c3 | 3 | 381.4 (22.1) | 156.1 (8.83) | 13.6 (0.45) | 6.26 (0.02) |
| **h2b** | 2 | 1599.1 (88.6) | 918.2 (50.8) | 62.6 (3.25) | 5.25 (0.07) |
| **h1** | 2 | 4306.3 (239.0) | 2716.5 (150.7) | 302.4 (16.6) | 8.88 (0.27) |
| **h1f** | 1 | 7037.0 (370.3) | 4970.0 (261.5) | 432.3 (22.7) | 8.60 (0.35) |
| **g** | 2 | 3815.5 (211.8) | 2836.8 (157.4) | 523.2 (28.8) | 14.03 (0.56) |
| **h1½** | 1 | 4304.8 (226.5) | 2715.9 (142.8) | 1174.7 (61.7) | 42.80 (2.15) |
| **g½** | 1 | 3832.7 (201.6) | 2839.5 (149.3) | 1637.7 (86.1) | 62.39 (3.18) |
| **b** | 2 | 15159.6 (842.0) | 10444.7 (580.0) | 3551.4 (197.1) | 52.28 (2.68) |
| **e** | 1 | 22212.0 (1169) | 15718.0 (827) | 10919.4 (575) | 341.5 (17.9) |
| **a** | 1 | 60606.8 (3190) | 44503.8 (2342) | 15358.5 (808) | 494.2 (25.9) |

**En iyi parametreler (hata modeli IV, Δχ²=1 profili — diğer parametre her
adımda yeniden optimize edilerek):**

| aday | parametreler |
|---|---|
| **c** | R∞ = **1.4513 ± 0.0885**, w = **0.0751 ± 0.0069** (τ₀ = 0.5146 sabit) |
| h2b | a = 9.761 ± 0.322, τ₀ = 0.5143 ± 0.0031 |
| c3 | R∞ = 1.4464 ± 0.0984, τ₀ = 0.5141 ± 0.0046, w = 0.0740 ± 0.0117 |
| h2 | a = 10.788 ± 1.122, τ₀ = 0.5114 ± 0.0041, b = −12.13 ± 5.73 |
| h1f | a = 8.612 ± 0.310 |

**Artık yapıları (hata modeli I fitleri, pencere son):**

| aday | 0.475 | 0.504 | 0.534 | 0.565 | 0.594 | 0.624 | 0.654 | 0.684 | 0.714 | 0.745 |
|---|---|---|---|---|---|---|---|---|---|---|
| h2 | +0.125 | +0.019 | −0.020 | −0.020 | −0.008 | +0.005 | +0.025 | +0.027 | +0.014 | −0.016 |
| c | +0.043 | +0.016 | −0.009 | −0.017 | −0.013 | −0.002 | +0.021 | +0.019 | −0.003 | −0.059 |
| c3 | +0.053 | +0.021 | −0.007 | −0.016 | −0.012 | −0.002 | +0.020 | +0.019 | −0.003 | −0.059 |
| h2b | −0.231 | −0.121 | −0.048 | +0.003 | +0.036 | +0.050 | +0.059 | +0.035 | −0.018 | −0.102 |
| **g** | **−0.607** | **−0.244** | −0.059 | +0.028 | +0.056 | +0.053 | +0.043 | +0.001 | −0.065 | −0.161 |
| g½ | −0.633 | −0.263 | −0.072 | +0.019 | +0.051 | +0.052 | +0.044 | +0.004 | −0.060 | −0.154 |

### 7b. KONTROL: aynı yarış R_tam üstünde (n = 20)

| aday | k | I | II | III | **IV** |
|---|---|---|---|---|---|
| h2 | 3 | 275.6 (15.9) | 150.4 (8.49) | 15.6 (0.57) | **7.45 (0.09)** |
| **g** | 2 | 989.1 (54.7) | 518.7 (28.6) | 48.7 (2.48) | **9.80 (0.32)** |
| **c** | 2 | 991.4 (54.9) | 539.6 (29.8) | 50.4 (2.58) | **10.79 (0.38)** |
| c3 | 3 | 921.6 (53.9) | 470.5 (27.3) | 36.0 (1.77) | 11.24 (0.31) |
| **h2b** | 2 | 1595.6 (88.4) | 1025.9 (56.8) | 82.0 (4.33) | 14.85 (0.60) |
| **g½** | 1 | 3403.0 (179.0) | 2537.3 (133.4) | 1585.1 (83.3) | 32.11 (1.58) |

**Burası raporun en ayırt edici karşılaştırmasıdır.** R_tam'da (c) ile (g)
arasındaki fark ΔAIC = 1.0 — **154'ün "ayırt edilemez" beraberliği aynen
duruyor** (154: ΔAIC = 0.09). R_lad'da ise (c), (g)'yi ΔAIC = 9.8 (model IV)
ve 3436 (model I) ile açık ara geçiyor; (g)'nin artığı negatif dalda
−0.607'ye çıkıyor. Fark, **negatif dalın fite katılmasından** geliyor: 154'ün
yarışı yalnız τ ≥ 0.54'te koşulmuştu, (c) ile (g) orada gerçekten
ayrılmıyorlar.

### 7c. Robustluk koşuları

| koşu | n | model I kazananı (χ²/dof) | model IV kazananı (ΔAIC ikinciye) |
|---|---|---|---|
| birincil (lad, taban 0.46, 0.03) | 20 | h2 (15.9); ≤2 par: **c** (20.9) | **c** (0.98) |
| lad, taban 0.46, **0.02 ince ızgara** | 30 | h2 (18.8); ≤2 par: **c** (31.0) | **c** (1.95) |
| lad, **ana taban 0.52** | 16 | h2 (7.70); ≤2 par: **c** (22.4) | g½ (0.03) ⚠ |
| lad, **ana taban 0.40** | 22 | h2 (17.3); ≤2 par: **c** (36.5) | **c** (1.13) |
| lad, **yalnız \|ΔRe\| < 0.06 bantları** (τ ≥ 0.59) | 12 | h2 (2.09); ≤2 par: **c** (28.3) | g½ (0.29) ⚠ |

⚠ = model IV'te ⟨σ⟩ = 0.40–0.45'e çıkıyor ve **her aday uyuyor** — "R = c(τ−½)/τ²"
bile χ²/dof = 0.07–0.08 veriyor. Bu satırlar bir kazanan bildirmiyor, hata
bütçesinin ayırt etme gücünün bittiğini bildiriyor.

**Aşırı-belirlemeyi geçen bantlarla koşu (son satır) ayrıca önemlidir:** o
kısıt negatif dalı ve sıfır bölgesini tümüyle siliyor (12 nokta, τ_eff ∈
[0.594, 0.745]) ve tam konvansiyon bütçesinde **R = sabit** bile uyuyor
(χ²/dof = 0.12, AIC 3.34). Yani (c)'nin (g)'yi yenmesi tümüyle
**aşırı-belirlemeyi GEÇEMEYEN** bantlara dayanıyor. Bu, §7a hükmünün en
ciddi niteliğidir.

---

# HÜKÜM

## (i) R_lad konvansiyonsuz nesne statüsünü HAK ETMİYOR

Üç bağımsız gerekçe, ikisi analitik:

1. **Sıfır-geçişi kanal değiştirmiyor (kimlik).** M_k(0) = M_emp bütün
   kanallarda bit düzeyinde aynı; her R_k tam olarak φ_Γ = arg M_emp
   noktasında sıfırlanır. 155'in 0.19σ_η² taban kayması R_lad'da
   **kaybolmuyor** — ölçülen yayılım R_tam'da 0.0131, R_lad'da 0.0140.
2. **Ağırlık operatörü tabana bağlı.** arg M_tam(R) ve M_emp altı tabanda
   **özdeş** (0.000e+00); arg M_lad(R̃) tabanlar arası **5.3–6.3 rad**
   oynuyor. R_tam ilkel φ_Γ'nın sabit yeniden-parametrizasyonudur; R_lad ise
   φ'nin konvansiyonuna X_lad'ın konvansiyonunu **ekler**.
3. **Ölçüm bunu doğruluyor.** Taban ekseni 0.28–0.58'e açıldığında R_lad'ın
   bant yayılımı %85–116, R_tam'ınki %15–21. 156'nın gördüğü ±%3, R_lad(taban)
   eğrisinin minimumunun tam o pencerede olmasıdır. Apsis kaydırması da
   kurtarmıyor (125% → 112%).

**Yani 156'nın "sıradaki adım" önerisinin dayanağı çürüdü.** Türetmenin doğru
değişkeni merdiven payı DEĞİL. 156'nın *fiziksel* hükmü — bağlaşımın
merdivende yaşadığı, η'nın fazı üretemediği — bu koşuyla çelişmiyor ve
tartışılmıyor; çürüyen, R_lad'ın **konvansiyonsuz bir tahminci** olduğu
iddiasıdır.

## (ii) Kapalı form MÜHÜRLENEMEDİ — ama tablo değişti

* **İstatistik hassasiyetinde (model I) hiçbir ≤2 parametreli form ayakta
  değil:** en iyisi (c), χ²/dof = 20.9. Üç parametreli φ-ilkel (h2) bile 15.9.
  Yani R_lad(τ)'nun şekli, ölçüm hassasiyetinde hiçbir basit kapalı formla
  uyuşmuyor.
* **Dürüst konvansiyon bütçesinde (model IV) hepsi uyuyor:** (c) 4.27, h2b
  5.25, c3 6.26, h2 6.27, h1f 8.60, h1 8.88 — ΔAIC ≤ 4.6, yani seçilemez.
  154'ün "iki eşdeğer form var, seçecek delil yok" hükmü aynen geçerli, ama
  şimdi **beş** eşdeğer form var.
* **Yeni ve gerçek bir ayrım:** negatif dal fite katıldığında (c) ile (g)
  arasındaki 154 beraberliği **R_lad'da** çözülüyor ((c) lehine, ΔAIC 9.8/3436),
  **R_tam'da çözülmüyor** (ΔAIC 1.0). Hangi formun kazandığı kanala bağlı —
  bu da (i)'nin bir başka yüzü. *Nitelendirme:* ayrım tümüyle
  aşırı-belirlemeyi geçemeyen (\|ΔRe\| > 0.06) bantlardan geliyor; o bantlar
  atıldığında hiçbir ayrım kalmıyor.
* **φ-doğrusallığı ÖLDÜ.** (h1/h1f/h1½) yarışın dibinde (χ²/dof 227–370,
  model I) ve φ'nin kendi fitinde kuadratik doğrusalı 10 koşunun 9'unda
  yeniyor. **φ-ilkel yaklaşımın kuadratik sürümü (h2) ise model I ve II'de
  BİRİNCİ** — yani "R'nin şekli tahmincinin kinematiğidir, fizik φ'dedir"
  hipotezi, cebirsel R formlarının hepsinden iyi çalışıyor. Bu, 154'ün
  "kapalı form makine analizinden aranmalı" yönünü destekliyor.

## (iii) τ₀_lad ve ½ ile ilişkisi

**τ₀_lad ayrı bir sayı değildir.** Konvansiyonsuz tanım şudur:

> τ₀\* = φ_Γ(τ) − arg M_emp(τ) denkleminin kökü — **her kanal için aynı**.

Ölçülen değerler (τ₀\*, kuadratik fit, τ_eff apsisi):

| taban | 0.28 | 0.34 | 0.40 | 0.46 | 0.52 |
|---|---|---|---|---|---|
| son | 0.5015 | 0.5063 | 0.5111 | 0.5143 | 0.5148 |
| orta | 0.5015 | 0.5068 | 0.5112 | 0.5144 | 0.5142 |

İstatistik hata ±0.0001–0.0004; **taban sistematiği 0.0132**, yani 33–130 kat
büyük (fit-konvansiyonu sistematiği ayrıca 0.0014, §4c).
σ_η² → 0 ekstrapolasyonu (155'in R2 okuması) 0.5198 / 0.5194 veriyor.

**½ ile ilişki, 155'in R1/R2 açmazı içinde AYNEN duruyor:** taban 0.28'de
τ₀\* = 0.5015 (½'den 0.0015 uzakta: **15σ istatistik, ama yalnız 0.11
σ_taban**), taban 0.52'de 0.5148. Yarış tarafından: τ₀ ≡ ½ sürümleri (g½, h1½) model
I–III'te ölü (χ²/dof 62–226), model IV'te canlı. **Konvansiyon ekseni
kapatılmadan "τ₀ ≠ ½" cümlesi kurulamaz** — 157 bu ekseni kapatmadı, yalnız
kanal değiştirmenin onu kapatmadığını **kanıtladı**.

## (iv) Bu koşunun tek KAZANIMI: konvansiyonsuz bir sayı

Taban konvansiyonu ilkel φ_Γ'yı %7 içinde salt öteliyor. Bu yüzden
φ'nin sıfırdaki **eğimi** tabandan (ve L'den) bağımsızdır:

> **a = dφ_Γ/dτ \|_{τ₀} = 10.759 ± 0.114**
> (10 koşu = 2 pencere × 5 taban; yayılım %3.7; aynı koşularda τ₀ 0.0131
> kayıyor. İki pencere ortalaması 10.766 vs 10.753 — ΔL = 0.566 ile
> ayırt edilemez.)

Eğrilik de aynı ailede ama daha zayıf: **b = −7.121 ± 0.776** (yayılım %34).
τ₀'ın *yeri* konvansiyona asılı, φ'nin *şekli* değil. Bir sonraki kalemin
hedefi bu olmalı: **a ≈ 10.76'nın (ve b ≈ −7.1'in) türetilmesi** — çünkü
τ₀'ın kendisinden farklı olarak bunlar konvansiyondan bağımsız ölçülmüş
sayılardır. (Sentetik gazlarda a'nın ne olduğu bu koşuda ölçülmedi; 155'in
φ tablolarından okunabilir ve doğrudan bir sonraki sınavdır.)

## Sıradaki adım (bu ölçümün işaret ettiği)

154 "kapalı form makine analizinden aranmalı" demişti; 156 "R_lad'a geç"
demişti. 157 ikincisini kapatıyor ve birincisini keskinleştiriyor:

1. **Hedef nesne R değil, φ_Γ'nın SIFIRDAKİ EĞİMİ olmalı.** a = 10.759 ±
   0.114, altı tabanda ve iki pencerede %3.7 içinde sabit; τ₀ ise 0.013
   kayıyor. Türetilecek sayı a'dır. İlk sınav ucuz: 155'in kayıtlı φ
   tablolarından **sentetik gazların a'sı** okunabilir (keskin, A4, N5, J14,
   P1). Gerçek/sentetik a oranı, 152–156'nın "×1.4" ve "×2" açıklarıyla aynı
   dişli mi, farklı mı — bu tek tabloyla anlaşılır.
2. **τ₀'ın konvansiyon ekseni ancak analitik olarak kapanır.** 155'in Açık
   1'i aynen duruyor: Γ_rot'un çıkarılmamış asal çizgi içeriğine
   duyarlılığının hesaplanması. 157 bunun için tek ek bilgi veriyor —
   duyarlılık, φ'yi ötelemekten ibaret (eğimi bozmuyor), yani hesabın
   üreteceği şey bir **kayma** olmalı, bir yeniden ölçekleme değil.
3. **Aşırı-belirlemenin τ ≲ 0.58'de kaybolması ayrı bir kalemdir.** Bu, yarışı
   ayıran bölgenin aynı zamanda tek-kanallı sağkalım-ağırlık resminin
   çalışmadığı bölge olması demek. İki-kanallı (156'nın 2K sırtı) ya da
   Gauss-dışı düzeltmeli bir ağırlık orada Re'yi kurtarıyor mu — sınanmadı.

---

## Dürüstlük notları

* **Uydurma yok.** Bütün tablo sayıları `157_configs/157_analiz.py` ve
  `157_configs/157_tablolar.py` ile `scratchpad/157/R_*.json`'lardan otomatik
  üretildi (`analiz_cikti.txt`, `tablolar.md`). Elle girilen tek şey
  154/155/156'nın referans değerleridir; her biri kaynağıyla işaretli.
* **156 ile çelişki nerede ve nerede DEĞİL.** 156'nın merdiven-kanalı fiziği
  (η fazı üretemiyor, merdiven aşırı-belirlemeyi geçiyor) bu koşuda
  tartışılmadı ve çürütülmedi. Çürütülen tek şey **"R_lad tabandan bağımsız
  (±%3)"** cümlesidir; 156 o cümleyi {0.46, 0.52, 0.58} penceresinde kurmuştu
  ve o pencerede bu koşu da ±%5–9 buluyor (T2'nin "156 penceresi" satırı).
  Fark yalnızca eksenin genişletilmesinden geliyor. 156'nın kendi T3'ünde de
  τ=0.585 bandında taban 0.46 → 0.58 arasında +%32'lik bir sıçrama vardı; o
  satır "0.58 tabanının kenarına düştüğü için" karşılaştırmaya alınmamıştı.
* **Aşırı-belirleme τ ≲ 0.58'de kayboluyor** (\|ΔRe\| = 0.05–0.13, eşik 0.06)
  — hem tam hem merdiven kanalında. Bu, yarışın en ayırt edici bölgesidir ve
  aynı zamanda modelin en zayıf olduğu bölge. §7c'nin son satırı bu kısıtla
  koşuldu ve **hiçbir ayrım kalmadı**; §7a'nın (c) > (g) hükmü bu nitelendirme
  ile okunmalıdır.
* **n_eff çöküşü düşük tabanlarda ciddi.** Taban 0.28'de merdiven kanalı τ̄ ≥
  0.685'te n_eff = 770 / 260 / 77'ye iniyor; taban 0.34'te τ̄ ≥ 0.715'te
  2223 / 1784. Bu satırlar tablolarda ⚠ ile duruyor, yayılım hesaplarına ve
  hata bütçelerine **girmedi**. Yayılım hükümleri (τ̄ = 0.595 / 0.625 / 0.655)
  bütün tabanların sağlıklı olduğu bantlardan alındı.
* **Sarma (phase-wrap) bayrağı iki satır yakaladı** (son, taban 0.28 ve 0.58,
  τ̄ = 0.745, R_tam kökü −2.333); ikisinde de faz artığı zaten 5.3 rad'dı ve
  artık süzgeci onları elemiş olurdu. Bayrak fazladan bir güvenlik ağıdır,
  hiçbir hükmü değiştirmedi.
* **154'ün 0.5153'ü ile bu koşunun 0.5143'ü.** Fark fit konvansiyonundandır,
  ölçümden değil: aynı veriye doğrusal fit **0.5154–0.5180**, kuadratik fit
  **0.5142–0.5147** veriyor; 32 konvansiyonun ortalaması 0.5157 ± 0.0015 (§4c). 155 bu sayıyı 0.5153 ± 0.0004 olarak yeniden
  üretmişti; 157 farklı bir (kuadratik, τ̄ ≤ 0.57, 0.02 ızgara) konvansiyonla
  koşuyor ve o konvansiyonu her yerde aynı tutuyor. Konvansiyonlar arası
  yayılım §4c'de tablolu.
* **τ₀(R_lad) fit aralığına asılı.** τ̄ ≤ 0.55 → 0.5141, τ̄ ≤ 0.63 → 0.5169
  (taban 0.46, son). Aynı taramada τ₀(φ) 0.5123 ↔ 0.5120 ve τ₀\* 0.5142 ↔
  0.5139 içinde kalıyor. R_lad'ın eğrisi kuadratikten en çok sapan olduğu için
  kökü de en oynak olan. Bu, (i)'nin dördüncü (bağımsız) gerekçesidir.
* **(h) ailesinin sayısal temeli.** R, JSON'a yazılmış arg M ızgarasının
  doğrusal ara değerli tersinden geliyor; V7 bu tersin koşudaki bisection
  kökünü 2.2e−4 içinde geri verdiğini gösteriyor (σ_jk ≈ 5e−3). (h)
  modellerinin χ²'sinde bu kaynaktan gelen hata ihmal edilebilir.
* **Model IV'ün ⟨σ⟩'sı çok büyük (0.35), χ²/dof'ları 0.02–0.07.** Bu, hata
  bütçesinin abartılı olduğunun değil, **R_lad'ın taban sistematiğinin gerçekten
  o kadar büyük olduğunun** işaretidir (§5a). Model IV bir "fit kalitesi"
  ölçüsü değil, "bu konvansiyon belirsizliğiyle hiçbir şey ayırt edilemez"
  ifadesinin niceliğidir.
* **Sentetik gaz koşulmadı.** Görev gerçek gazın iki penceresini istedi; 152/153/154'ün
  sentetik gazları bu kalemde kullanılmadı. (iv)'ün işaret ettiği a-ölçümü
  sentetikte de yapılmalıdır — yapılmadı.
* **Süreler ölçüt değildi.** 22 koşu 8 çekirdekte kısmen eşzamanlı; tipik koşu
  2.7–7.0 dk (`kaba`), 8–14 dk (`ince`). Sayısal sonuç etkilenmez.

---

## Ek — dosyalar ve tekrar-üretim

| dosya | ne |
|---|---|
| `157_configs/157_cekirdek.py` | ölçüm çekirdeği: 154+156'yı import eder, çizgi-bazında kayıt + τ_eff + arg M ızgarası + sarma bayrağı |
| `157_configs/157_kos.py` | `<son\|orta> <taban> <kaba\|ince\|tau0>` — koşu sürücüsü |
| `157_configs/157_analiz.py` | A/A2/B/B2/C/D/E bölümleri: tayf, τ₀ kimliği, taban çöküşü, yarış, φ şekli |
| `157_configs/157_tablolar.py` | rapor tablolarını markdown olarak üretir |
| `157_configs/157_dogrulama.py` | V1–V7 kopya-kayması ve kimlik denetimleri |
| `157_configs/157_figur.py` | 6 panelli figür |
| `157_Rlad_yarisi.png` | figür |

Ham çıktılar `scratchpad/157/`: `R_<pencere>_t<taban>_<ızgara>.json` (22 adet,
çizgi-bazında kayıt + arg M ızgarasıyla), `kanal_*.npz` (drift/lad/eta
önbellekleri), `log_*.txt`, `analiz_cikti.txt`, `tablolar.md`,
`dogrulama.txt`, `kos_a.sh`, `kos_b.sh`.
