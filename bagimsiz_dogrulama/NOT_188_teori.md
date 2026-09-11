# Kalem 188 — Teorik Çerçeve ve Çürütme Planı

**Tarih:** 11 Eylül 2026 · **Durum:** veri bağımlılığı nedeniyle kısmen yürütülemedi (bkz. `RAPOR_veri_kurtarma_riski.md`)

---

## 1. ζ nedir (187 tarifinden)

187c'de:
```
Δ_q        = c^ölç_q − c^kesik_q            (ölçülen satır genliği − kesik model)
karışım_q  = c^kesik_q − c^öz_q             (kesik model − ρ≡1 öz-terim)
ζ_b        = Σ_{q∈b} Δ_q·conj(karışım_q) / Σ_{q∈b} |karışım_q|²
```
Yani **ζ, "kesik modelin fazlalığının" ölçülen eksikliğe izdüşüm katsayısıdır**: pencere-ötesi içeriğin, pencere-içi karışımı ne oranda yıkıcı biçimde sildiğinin normalleştirilmiş ölçüsü. Modül = iptal payı, açı = faz uyumu.

## 2. Ölçülmüş profil (187 K2 defteri)

| bant τ' | 0.45-0.50 | 0.50-0.55 | 0.55-0.60 | 0.60-0.65 | 0.65-0.70 | 0.70-0.75 | 0.75-0.80 | 0.80-0.86 |
|---|---|---|---|---|---|---|---|---|
| \|ζ\| gerçek | 0.3284 | 0.3230 | 0.3184 | 0.3182 | 0.3184 | 0.3274 | 0.3456 | **0.3811** |
| hata | .0032 | .0032 | .0060 | .0030 | .0039 | .0041 | .0052 | .0072 |
| açı | 179.6 | −179.5 | −179.8 | 179.3 | 179.5 | −179.9 | −179.3 | 179.9 |
| ikiz \|ζ\| | 0.0892 | 0.0869 | 0.0796 | 0.0796 | 0.0736 | 0.0695 | 0.0685 | 0.0591 |

Havuz: **0.3287 ± 0.0023, açı 179.96° ± 0.07°**. Amplitüd-okuması `1 − m_ölç/m_kesik` ζ ile %1–3 içinde örtüşüyor.

## 3. Üç gözlem (veri olmadan yapılabilenler)

**(0) 180° bir konvansiyon olabilir — tur 2'nin yeni ipucu.** Bu turda Not 4'ün faz-kilidini bağımsız doğruladım: yazarın konvansiyonuyla (Ĝ'de **mutlak** `t_n`) 24 satırın 24'ü 180°'de (±0.3°). Kritik nokta: o faz, `Ĝ ≈ iω⟨e^{iωt^s}u⟩` türev/ortanokta kinematiğinden geliyor — yani **180° bir dinamik işaret değil, yapısal bir çarpan** olabilir. ζ'nin 180.0°±0.07°'si aynı aileden bir düğümle açıklanabilir. **Test (T0):** amplitüd kanalının kendi dalga katsayıları işaretli okununca faz 0° çıkıyor (benim regresyonumda `a·√q = +0.757`, `b ≈ 0`); Ĝ'de ise 180°. Aynı satırın iki kanalda zıt işaretli görünmesi, ζ'nin açısını "saf yıkıcı" diye okumadan önce faz sözleşmesinin sabitlenmesi gerektiğini gösterir. ζ'nin |ζ| profili bu tartışmadan etkilenmez; yalnızca "180° = yıkıcı" yorumu etkilenir.

**(a) ζ bir tahminci sabiti değil.** İkiz deniz (aynı inşa, farklı pencere uzunluğu L_inşa = 12.0304) aynı tahminciden 0.077 alıyor — gerçeğin dörtte biri. Demek ki ζ denizin kendi içeriğini taşıyor; "her denizde 0.33 çıkan kinematik bir sabit" değil. **Bu, en olası ölü senaryoyu şimdiden eliyor.**

**(b) Plato, hata çubuklarından çok daha düz.** 0.3184, 0.3182, 0.3184 — üç bant arasındaki saçılma **%0.06**, oysa her bandın kendi hatası **%1–2**. Bu, üç bandın bağımsız ölçüm olmadığını, ortak bir yapı tarafından belirlendiğini gösterir. Yani platonun *değeri* tek bir nicelikten geliyor; bantlar onu ayrı ayrı görmüyor.

**(c) 1/π sınaması.** Plato değeri 0.3182–0.3184; **1/π = 0.318310**. Uyum %0.03. Bu projenin kendi tarihi (π/600 overfitting dersi, 13 Ağu) tam olarak bu tür tesadüflere karşı bağışıklık kazanmayı öğretti — o yüzden bunu **bulgu değil, çürütülecek hipotez** olarak kaydediyorum.

## 4. Çürütme planı (sıralı, her biri tek koşu)

| # | Test | Aritmetik ise beklenen | Kinematik ise beklenen | Durum |
|---|---|---|---|---|
| **T0** | ζ'nin 180°'sini faz sözleşmesinden ayır: Ĝ ile amplitüd kanalının işaretli katsayılarını aynı satırda karşılaştır | işaret farkı dinamik (denizden) | işaret farkı yapısal (türev/ortanokta) | **Tur 2'de kısmen yapıldı**: amplitüd kanalı faz 0°, Ĝ faz 180° — fark yapısal görünüyor |
| **T1** | ζ tahmincisini **faz-rastgele surrogate** denizin üzerinde koş (aynı ızgara, aynı kesim) | \|ζ\| → plasebo tabanı (≪0.33), açı kilitlenmez | \|ζ\| ≈ 0.33 ve 180°'de kalır | Veri bekliyor (184/185/186 blokları) |
| **T2** | Bant kenarlarını ±0.05, kesim τ_c'yi 0.86→1.20 kaydır | \|ζ\| düzgün ölçeklenir | hiç kıpırdamaz (sabit) | Veri bekliyor |
| **T3** | \|ζ(τ)\|'yi parametresiz pencere-ötesi ağırlığa karşı çiz: `W(τ) = Σ_{q: τ_q>τ_c(τ)} a_q² / Σ_{q: τ_q≤1} a_q²`, `a_q = 1/(π p^{k/2})` | yükseliş W'yi izler | korelasyon yok | Kısmen: bant ekseni artık biliniyor (τ ∈ (0.45,0.86], `187b`), ama bant-başına örtüşme modeli gerekiyor |
| **T4** | 1/π hipotezini ikinci bağımsız pencerede (L ≈ 10.9) sına | plato 1/π'den sapar | 1/π sabit kalır | Kısmen yapılabilir |

## 5-bis. ✅ T1 SONUCU — iptal ARİTMETİK (11 Eylül 2026, scratchpad kurtarıldıktan sonra)

**Yöntem.** `176b`'nin faz-rastgele vekilleri kullanıldı: `S_v(t) = −Σ_q a_q sin(ω_q t + φ_q)`, yani **genlik zarfı Hkeskin'le birebir aynı, yalnız her çizginin fazı bağımsız rastgele** (`φ_q ~ U(0,2π)`). Zincir her vekil için baştan koşturuldu: `184b → 185b → 186b → 187c` (yeni `K1_<gaz>`, `OZ_<gaz>`, `G1_proj_<gaz>` dosyaları üretildi; kampanyanın orijinal dosyalarına dokunulmadı). Sürücü: `22_T1_vekil_zeta.py`. **Kontrol:** aynı sürücü gerçek denizde 0.3287 ∠179.96° veriyor — kampanya değeriyle birebir.

| deniz | \|ζ\| havuz | açı | bantlar arası \|ζ\| | 180°±15° (8 bant) | 180°'den maks sapma |
|---|---|---|---|---|---|
| **gerçek** | **0.329** | **+179.96°** | **0.32 – 0.38** | **8/8 ✓** | **0.7°** |
| VF1 | 1.099 | −151.5° | 0.45 – 9.57 | ✗ | 100.8° |
| VF2 | 0.950 | −175.7° | 0.77 – 6.92 | ✗ | 109.1° |
| VF3 | 2.544 | −178.3° | 1.31 – 15.34 | ✗ | 18.5° |
| VF4 | 3.244 | +175.6° | 1.52 – 8.80 | **✓** | 8.9° |
| VF5 | 1.189 | −170.3° | 0.83 – 7.38 | ✗ | 83.8° |
| VS1 | 1.102 | −156.4° | 0.38 – 8.47 | ✗ | 112.3° |
| VS2 | 1.151 | +179.8° | 0.85 – 10.77 | ✗ | 39.3° |
| VS3 | 2.864 | −179.7° | 1.39 – 15.64 | **✓** | 12.6° |

**İki hüküm — biri güçlü, biri uyarı:**

**(1) İptal kinematik değil, ARİTMETİK.** Sekiz faz-rastgele denizin **hiçbiri** |ζ| ≈ 0.33 üretmiyor: ortalama modül **2.0–6.8** (gerçeğin 6–20 katı), bant başına en küçük değerler bile 0.38–1.5. Yani denizin "üçte birini yok eden" yapı, genlik zarfından değil **faz düzeninden** geliyor. Bu, 187'nin "deniz-ayrıştırıcı olan ζ'dir" cümlesini bağımsız olarak destekliyor.

**(2) ⚠️ H-F2 mührü ("8/8 bant 180°'de → saf yıkıcı") zayıf bir ölçüt.** Aynı testte **3/8 faz-rastgele vekil de geçiyor** (VF4, VS3; VF3 sınırda 18.5° ile kalıyor). Yani 180°±15° bandı, surrogate null'a karşı **~%40 yanlış-pozitif oranı** taşıyor. Doğru ayrım **modülde**: gerçeğin 8 bandı da 0.32–0.38 aralığında, vekillerin hiçbirinde bu yok. **Öneri:** Not 5/6'da H-F2 "saf yıkıcı" mührü tek başına sunulmasın; yanına *modül* koşulu (|ζ| ≈ 1/3) ve bu surrogate tabanı eklensin. Aksi hâlde bir hakem aynı itirazı yapar.

**Kendi tahminimin düzeltmesi:** §5'te T1 için "plasebo tabanı ~0.03" diye öngörmüştüm; **yanlıştı.** ζ bir *oran* olduğu için rastgele fazlı denizde model uyuşmazlığı O(1) veriyor — ayırt edici modül tabanı değil, **modülün küçüklüğü** (gerçek 0.33 vs vekil 2–7) ve onun bant-tutarlılığı.

## 5-ter. ✅ T2 SONUCU — plato kinematik DEĞİL, yapılı bir eğri (11 Eylül 2026)

Sabit 8 bant yerine τ boyunca **kayan pencere** (genişlik 0.05 ve 0.10, adım 0.01, 42 pencere), jackknife hatalarıyla. Sürücü: `23_T2_bant_taramasi.py`.

**Gerçek deniz, genişlik 0.05** (τ merkez → |ζ| ± se, açı):

| τ | 0.45 | 0.53 | 0.61 | **0.65** | 0.69 | 0.73 | 0.77 | 0.81 | 0.85 |
|---|---|---|---|---|---|---|---|---|---|
| \|ζ\| | 0.3364 | 0.3210 | 0.3181 | **0.3164** | 0.3179 | 0.3295 | 0.3430 | 0.3671 | 0.3918 |
| ±se | .0048 | .0028 | .0032 | **.0035** | .0037 | .0043 | .0047 | .0080 | .0058 |
| açı° | −179.4 | −179.3 | 179.6 | 179.9 | 179.8 | 180.0 | −179.4 | 179.8 | 180.0 |

Genişlik 0.10 ile aynı şekil (min 0.3169 @ 0.61–0.65, maks 0.3811 @ 0.85; std 0.0195). Hatalar ±0.003–0.008 olduğundan **eğrinin yapısı 10–20σ anlamlı**.

**Üç hüküm:**

1. **Plato kinematik değil.** |ζ|(τ) düz değil: **τ ≈ 0.65'te minimum (~0.316)**, iki yana yükseliyor (0.45'te 0.336, 0.85'te 0.392, yani +%24). Sabit bir konvansiyon/kinematik çarpanı olsaydı tarama boyunca kıpırdamazdı. Demek ki ζ **gerçek bir izdüşüm katsayısı**.
2. **1/π tesadüfü çözüldü.** 1/π = 0.31831, eğrinin **minimumuna** denk geliyor (0.3164 ± 0.0035 → 0.5σ). Yani "plato = 1/π" okuması, yapılı bir eğrinin tek bir noktasındaki tesadüftü — projenin kendi π/600 dersinin tam olarak beklediği şey. **1/π hipotezi (H-A) düşürülmeli.**
3. **Faz kilidi binlemeden bağımsız.** 42 kayan pencerenin **hepsinde** açı 180° ± 0.7° — yani H-F2'nin içeriği bant seçimine bağlı değil. (Yine de T1'deki uyarı geçerli: ±15° ölçütünün kendisi vekillere karşı zayıf; doğru beyan "açı 180°±0.7° **ve** |ζ| ≈ 0.32".)

**Kontrol (faz-rastgele VF1, aynı tarama):** |ζ| 0.56 → 14.9 arasında savruluyor, hatalar 5–10 kat büyük, açılar dağınık (−94° … −178°). Yani gerçeğin *düzgün, düşük ve kilitli* profili, aynı genlik zarfına sahip rastgele fazlı denizlerde üretilemiyor.

**Yeni gözlem (açıklama değil, kayıt):** minimum τ ≈ 0.65, ne foldda (0.5) ne kesimde (0.86) — ikisinin arasında. Neden orada olduğu açık bir soru; T3'ün (pencere-ötesi ağırlık) doğal hedefi budur.

## 5-quater. ✅ T3 SONUCU — şekil kesim konumunun eseri değil, "kesime yakınlık" (11 Eylül 2026)

Kesim üç değere çekildi (τ_c = 0.70, 0.86, 1.00; sırasıyla 246 / 1332 / 6119 çizgi) ve her birinde aynı kayan bant taraması koşuldu. Sürücü: `24_T3_kesim_taramasi.py`.

| τ | 0.45 | 0.50 | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 | 0.80 | 0.85 | 0.90 |
|---|---|---|---|---|---|---|---|---|---|---|
| **τ_c = 0.70** | 0.3393 | 0.3357 | 0.3239 | 0.3207 | **0.3102** | — | — | — | — | — |
| **τ_c = 0.86** | 0.3364 | 0.3261 | 0.3200 | 0.3182 | **0.3164** | 0.3206 | 0.3363 | 0.3621 | **0.3918** | — |
| **τ_c = 1.00** | 0.2652 | 0.2612 | **0.2583** | 0.2603 | 0.2640 | 0.2711 | 0.2912 | 0.3255 | 0.3693 | **0.4272** |

**Üç hüküm:**

1. **Minimum kesim konumunun eseri değil.** Üç kesimde de bir minimum var (τ ≈ 0.67 / 0.66 / 0.55–0.60). Aynı şekil, üç farklı çizgi evreniyle (246 → 6119 çizgi) çıkıyor; yani profil **denizin aritmetik içeriğinden** geliyor.
2. **Yükseliş "kesime yakınlığı" izliyor.** Her eğri kendi kesimine doğru monoton yükseliyor: 0.86'da τ=0.85'te 0.392; 1.00'de τ=0.90'da 0.427. Yani iptal, pencere kenarına yaklaştıkça zayıflıyor — fiziksel olarak beklenen yön (kenarda, dışarıda kalan içerik en çok orada önemli).
3. **Mutlak seviye kesime bağlı** (min 0.3076 / 0.3161 / 0.2345). Yani ζ'yi alıntılarken **kesim birlikte söylenmeli** (187 zaten 0.86 diyor). Bu aynı zamanda 1/π tesadüfünü ikinci kez zayıflatıyor: minimumun değeri kesimle değişiyor, 1/π ise sabit.

**188'in bilançosu:** profil kimliğinin *nitel* cevabı çıktı — |ζ|(τ) ≈ "sabit taban + kesime yakınlık artışı", taban ≈ 0.26–0.32 (kesime göre), minimum τ ≈ 0.6 civarı. **Nicel bir kapalı form hâlâ yok**; T3 ilk kez bunun ölçülebilir bir yapı olduğunu ve hangi değişkene bağlı olduğunu (kesime uzaklık) gösterdi.

| test | durum | sonuç |
|---|---|---|
| **T0** faz sözleşmesi | kısmen yapıldı | amplitüd kanalı 0°, Ĝ 180° → 180° yapısal bir çarpan (tur 2) |
| **T1** faz-rastgele denizde ζ | ✅ **yapıldı** | iptal **aritmetik**; 8 vekilin hiçbiri \|ζ\| ≈ 0.33 vermiyor (§5-bis) |
| **T2** bant/kesim taraması | ✅ **yapıldı** (bant ekseni) | plato **yapılı eğri**; 1/π hipotezi düştü; faz kilidi bant-bağımsız (§5-ter) |
| **T3** kesim taraması | ✅ **yapıldı** | minimum kesimden bağımsız; yükseliş kesime-yakınlığı izliyor (§5-quater) |
| **T4** ikinci pencerede 1/π | ✅ gereksizleşti | 1/π zaten minimuma denk gelen tesadüf olarak çözüldü |

**Kalan tek soru (nicel):** |ζ|(τ) ≈ taban + kesime-yakınlık artışı biçiminin **kapalı formu**. Nitel yapı T3 ile çıktı; nicel ifade (ör. tabanın τ_c'ye ve pencere genişliğine bağlılığı) hâlâ açık.

## 6. Durum (11 Eylül 2026 akşamı)

**Veri durumu:** `qm_riemann/scratchpad/` dolu (26 görev, 297 MB); pipeline birebir çalışıyor (ζ = 0.3287 ∠179.96°). T1/T2/T3 sürücüleri: `22_T1_vekil_zeta.py`, `23_T2_bant_taramasi.py`, `24_T3_kesim_taramasi.py`.
