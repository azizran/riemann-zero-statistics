# 168 — T3: c'NİN ANALİTİK TÜRETİMİ (Ç4 REZONANS-İNTEGRALİ) (3 Eylül 2026)

Durum: **A1 BİTTİ · A2 BİTTİ · A3 BİTTİ** (E060 çapraz sınavı ayrı bölümde).
167'nin bıraktığı iki sert kısıt (c(kesim) ve c(pencere)) türetimin karşısına
kondu ve **türetim bu kısıtları ÜRETEMEDİ**; onun yerine kısıtların hangi
çarpandan geldiği ölçüldü ve **kesim ekseni için sayısal olarak kapanan bir
yasa** bulundu. Aşağıda hangi adımın tuttuğu, hangisinin düştüğü açık yazılı.

Scriptler: `168_configs/168_olcek.py` (ucuz ayrıştırma, yalnız C_*.json),
`168_configs/168_rezonans.py` (★'nın sayısal değerlendirmesi + yerleşim
oranı), `168_configs/168_profil.py` (κ'nın tepe profili, merkez fazı
çıkarılmış). Ham çıktılar `scratchpad/168`, `scratchpad/167`.

---

## 0. TEK CÜMLELİK HÜKÜM

> **Ç4'ün ortak kazancı, birinci-mertebe rezonans integraliyle (★)
> AÇIKLANMIYOR: ★ ölçülen Ç4'ün ancak %5'ini veriyor, çünkü ★ `u`'da
> doğrusal ve geçerlilik penceresi `τ ≤ 1/(2πσ_Ĉ) = 0.573` iken Ç4 dörtlü
> kombinasyon frekanslarında (τ_etkin ≈ 2–3) yaşıyor. Aynı yerde ÜÇ bağımsız
> DOYUM imzası ölçüldü (★'nın ×18 açığı; kesim yasasının genlik değil ÇİZGİ
> SAYISI ile gitmesi; λ-ablasyonunun neredeyse boş çıkması). c bir rezonans
> integrali değil, DOYMUŞ bir Gram/Debye–Waller çarpanıdır.**

Bu düşüşün yanında ÜÇ şey kapandı:
1. **c(pencere) ölçek yasası rezonans fiziği DEĞİL** — tamamı çizgi-uyum
   sızıntısı `g_cal = g_E g_X²`'nin sonlu-örneklem yasasıdır (§A2.5).
2. **c(kesim) ölçek yasası SAYISAL OLARAK KAPANDI:**
   `θ/θ₀ = 1 − 0.2175·φ`, φ = BOŞ ÇİZGİ KESRİ; DÖRT kesim gazında ±%0.7,
   biri **ÖRNEKLEM-DIŞI** (E060 önce öngörüldü 0.7908, sonra ölçüldü
   0.7945; §A2.7/§5). Bu bir ÖLÇÜLEN yasa; katsayı türetilmedi (açık borç).
3. **c hangi üyeye ait: `W_X` üyesi** (c ≈ 0.4051/0.4053), türetimden
   gerekçelendirildi (§A2.8). `0.5647` aynı sayının `W_amp` ile bir kez
   fazla bölünmüş hâli.

---

## A1 — KALEM: Ç4 TERİMİ ve c'NİN KAPALI İFADESİ

### A1.1 Ç4'ün tek-site çekirdeğindeki yeri

165 §2a/2b'nin indirgemesinde bütün alanlar tek sitede (`s_n ≡ m_{n+1}`)
yazılır ve ölçüm nesnesi budamasız dörtlü toplamdır:

```
J_2(W) = ⟨E_mod(s) X_mod(s)² e^{−iW s}⟩_n
       = (1/8) Σ_{(1,ε₁)(2,ε₂)(3,ε₃)} h'^{(ε₁)}_1 y^{(ε₂)}_2 y^{(ε₃)}_3 · κ(Δ)
Δ = ε₁ω₁ + ε₂ω₂ + ε₃ω₃ − W ,    κ(ν) = (1/N) Σ_n e^{iν s_n}
```

`κ`'nın üç yapısı vardır (165 §4d): (i) `Δ ≈ 0`'da Dirichlet zarfı,
(ii) merkez fazı `arg κ ≈ ν s̄`, (iii) **`Δ ≈ ±ω_r`'de TARAK REZONANSI**.
Ç4, (iii)'ün taşıdığı paydır ve ölçülen toplamın %66–78'idir.

### A1.2 Rezonans integralinin kapalı formu (★)

Tarağın kendisi bir merdivendir: `N̄(t)+S(t) = n−½` ⇒

```
dN/dt = N̄'(t)[1 + u(t)] ,  u(t) = 2 Σ_r 𝒢_r cos(ω_r t) ,
𝒢_r = −π τ_r a_r cos(π τ_r)                       (143'ün yasası)
```

(Türetim: `S(t) = −Σ_r a_r sin(ω_r t)` ⇒ `S'/N̄' = −2π Σ_r τ_r a_r cos(ω_r t)`;
`cos(πτ_r)` çarpanı `s_n`'in yarım-gap ORTA NOKTASINDA örneklenmesinden.)

Her pencere ortalaması ikiye ayrılır ve `f = E X² e^{−iWt}` konunca

```
★   J_2(W) = Ĵ_D(W) + Σ_r 𝒢_r [ Ĵ_D(W−ω_r) + Ĵ_D(W+ω_r) ] + O(u²)
    Ç4(W)  = Σ_r 𝒢_r [ Ĵ_D(W−ω_r) + Ĵ_D(W+ω_r) ]
```

`Ĵ_D` = DÜZGÜN tarak tayfı (165 §4b'nin kontrolü). Yani **Ç4, düzgün-tarak
tayfının merdiven frekanslarında ÖRNEKLENMESİ ve 𝒢 ile tartılmasıdır.**
(★ 167_rezonans.py'de kurulmuştu; burada aynen kullanılıyor.)

### A1.3 c'nin kapalı ifadesi: ayrık örnekleme ↔ profil integrali

`Ĵ_D` kendisi bir tepe dizisidir:

```
Ĵ_D(ν) = Σ_j A_j κ_D(ν − ν_j) ,   ν_j = ε₁ω₁+ε₂ω₂+ε₃ω₃  (kombinasyon frekansları)
κ_D(ν) = (1/T)∫_pencere e^{iνt}dt = e^{iν t̄} · sinc(νT/2)
```

**Sinc ailesinin sabiti:** `∫κ_D dν = ∫sinc²... = 2π/T`, yani

```
w_eff ≡ ∫κ_tepe dν / κ_tepe(0) = 2π/T = dres            (★★)
```

Model, ★'yı **ÇIPLAK AYRIK** olarak, tam `ν = W ∓ ω_r` noktalarında
değerlendirir. Fiziksel karşılığı ise tepe profilinin, kombinasyon-Δ
yoğunluğu `n_Δ` ile ağırlıklı integralidir. Görevin istediği oran:

```
          n_Δ(ω_r) × ∫κ_tepe        n_Δ(ω_r) · 𝒢_r · w_eff
c_pred = ────────────────────── = ────────────────────────── = n_Δ(ω_r) · w_eff
          çıplak ayrık değer            𝒢_r · κ_D(0)
```

yani **c_pred = bir çözünürlük genişliğine düşen (genlik-ağırlıklı)
kombinasyon sayısı** — rezonansın "doluluk oranı". Boyutsuzdur
(`n_Δ` ~ 1/frekans, `w_eff` ~ frekans) ve doygun olmayan rejimde ≤ 1'dir.

### A1.4 SEMBOLİK ÖLÇEK YASALARI (A1'in üç öngörüsü)

**(i) Pencere ekseni.** `n_Δ` yalnız merdivene bağlıdır (T'den bağımsız);
`w_eff = 2π/T`. Dolayısıyla

```
c ∝ T^{−1}    ⇒   c(T/2)/c(T) = 2.00 ,   c(T/4)/c(T) = 4.00
```

**T-bağımlılığı BİRİNCİ kuvvettedir** ve işaret POZİTİF'tir (pencere
kısaldıkça yerleşim BÜYÜR).

**(ii) Kesim ekseni.** Kesim hem `n_Δ`'yı (kombinasyon üreten çizgiler) hem
ağırlıkları (`𝒢_r`) değiştirir. Merdiven yoğunluğu `dn/dω = e^ω/ω` ve
`𝒢 = −τ cos(πτ)/(m√q)` olduğundan `1/q = e^{−τL}` TAM sadeleşir ve

```
Σ_r 𝒢_r² ≈ ∫_0^{τ_c} τ cos²(πτ) dτ
         = τ_c²/4 + τ_c sin(2πτ_c)/(4π) + [cos(2πτ_c) − 1]/(8π²)
```

(sayısal merdivenle uyum %2.5–14; §4'ün notu). Kesim bağımlılığı bu
TOPLAMLARIN ORANINDADIR:

```
c(τ_c)/c(0.95) = Σ𝒢²(τ_c)/Σ𝒢²(0.95)   [genlik-kare ağırlıklı]  ya da
                 Σ|𝒢|(τ_c)/Σ|𝒢|(0.95)  [genlik ağırlıklı]
```

Sayıları: `τ_c = 0.70` → **0.245** (kare) / **0.105** (mutlak);
`τ_c = 0.90` → **0.778** / **0.669**; `erfc(0.68)` → **0.207** / **0.125**.

**(iii) λ ekseni.** `𝒢_r ∝ a_r ∝ λ` ve `Ĵ_D ∝ λ³` ⇒ `Ç4/Ĵ_D ∝ λ`. Ölçülen
`Ç4/tot ≈ 0.73` ⇒ `Ç4/Ĵ_D ≈ 2.7`; λ: 1.00 → 0.85 kazancı `3.7 → 3.30`'a
indirir, yani **KALİB_u2 +%12 artmalı**.

### A1.5 ★'NIN GEÇERLİLİK SINIRI (türetimin kendi iç uyarısı)

★ `u`'da BİRİNCİ mertebedir. Oysa sitelerin gerçek yeri
`s_n = s̄ + ḡ(n−n̄) + ḡ Ĉ_n` (σ_Ĉ = 0.278) ve `κ`'nın tam çarpanı
`e^{iν ḡ Ĉ}`'dir. Doğrusallaştırma ancak

```
ν ḡ σ_Ĉ = 2π τ σ_Ĉ ≲ 1   ⇒   τ ≤ 1/(2π σ_Ĉ) = 0.573
```

iken geçerlidir — bu **165 §7'nin FAZ SEĞİRMESİ ölçütünün ta kendisidir.**
Oysa Ç4, DÖRT merdiven frekansının kombinasyonlarında yaşar
(`τ_etkin ≈ 2–3`) ⇒ `2πτσ_Ĉ ≈ 3.5–5 rad`. **Yani ★ Ç4 için geçerlilik
penceresinin tamamen dışındadır ve ciddi biçimde AZ vermek zorundadır.**
Bu, A1'in kendi öngörüsüdür; A2 onu ölçer.

**Doğru (doyumlu) çerçeve:** üstel faz korunursa model alanının gerçek
tarak üzerindeki Gram matrisi

```
G_{qq'} = ⟨e^{i(ω_q−ω_{q'}) s}⟩ = κ(ω_q − ω_{q'})
```

olur: köşegen 1, köşegen-dışı ise **`log(q/q')` bir merdiven frekansına
denk düştüğünde tarak rezonansı `𝒢`, `|log(q/q')| ≲ dres` olduğunda ise
Dirichlet tepesi**. Model alanının şişmesi tam olarak
`h'†G h'/h'†h'`'dir ve ölçülen karşılığı `1/g_E`, `1/g_X`'tir. Rezonans
integralinin doğru adresi budur — A2 bu adresi doğruluyor.

---

## A2 — SAYISAL YÜZLEŞME (yeni gaz kurulmadı)

### A2.1 κ'nın TEPE PROFİLİ ölçüldü — (★★) DOĞRULANDI

`168_profil.py Hkeskin` (merkez fazı `e^{iνs̄}` çıkarılmış; ham `Re κ` ile
ölçmek yanlış — faz ızgara adımından hızlı döner, ilk deneme `w_eff/dres =
0.027` verdi, bu bir artefakttır).

**Δ ≈ 0 çapası (165 §4d'nin yeniden üretimi):**

| Δ/dres | 0.00 | 0.25 | 0.50 | 0.75 | 1.00 | 2.00 | 10.0 |
|---|---|---|---|---|---|---|---|
| \|κ̃\| ölçülen | 1.0000 | 0.9003 | 0.6366 | 0.3002 | 0.0020 | 0.0010 | 0.0002 |
| Dirichlet | 1.0000 | 0.9003 | 0.6366 | 0.3001 | 0 | 0 | 0 |

**`ω_r` çevresindeki tepe profili (adım dres/8, ±5·dres):**

| q | τ_r | 𝒢_r | \|κ̃(ω_r)\| | **∫κ̃dx/κ̃(0) [dres]** | ∫\|κ̃\|dx/\|κ̃\|(0) | \|κ̃\|(±dres)/\|κ̃\|(0) |
|---|---|---|---|---|---|---|
| 2 | 0.0576 | −0.040078 | 0.039963 | **1.0396** | 1.815 | 0.0031 |
| 3 | 0.0913 | −0.050572 | 0.050460 | **1.0392** | 1.815 | 0.0028 |
| 7 | 0.1618 | −0.053413 | 0.052992 | **1.0396** | 1.815 | 0.0030 |
| 43 | 0.3127 | −0.026470 | 0.025313 | **1.0405** | 1.823 | 0.0037 |
| 211 | 0.4449 | −0.005276 | 0.004207 | **1.0344** | 1.885 | 0.0309 |
| 1747 | 0.6206 | +0.005492 | 0.005724 | 1.1365 | 1.901 | 0.0375 |
| 12073 | 0.7813 | +0.005497 | 0.005364 | 1.2148 | 3.233 | 0.0845 |
| 26711 | 0.8473 | +0.004599 | 0.005082 | 0.2975 | 3.477 | 0.1382 |
| 58391 | 0.9123 | +0.003633 | 0.003776 | 0.1931 | 5.800 | 0.6385 |
| 80687 | 0.9392 | +0.003246 | 0.004159 | 0.6945 | 6.768 | 0.4416 |

> **HÜKÜM (★★ TUTTU).** `τ_r ≤ 0.45` çizgilerinde tepe integrali
> **`w_eff = (1.039 ± 0.002)·2π/T`** — A1'in sinc-ailesi sabiti **%4
> içinde doğrulandı**. `|κ̃(ω_r)| ↔ |𝒢_r|` uyumu: `τ_r ≤ 0.32`'de
> **%0.2–4.4**, `τ_r ≥ 0.44`'te %2–28 (143'ün yasası düşük τ'da bir kez
> daha; yüksek τ'da tepeler örtüştüğü için bozuluyor — A2.2).
>
> **Ama `τ_r ≳ 0.78` üstünde profil AYRIŞMIYOR:** `|κ̃|(±dres)/|κ̃|(0)`
> 0.003'ten 0.64'e çıkıyor, `∫κ̃` kararsızlaşıyor (1.21 → 0.19 → 0.69).
> **Tepe kavramı merdivenin üst ucunda ÖLÜYOR.**

### A2.2 MERDİVEN ÇÖZÜNÜRLÜĞÜ (A2.1'in nedeni)

Ardışık merdiven çizgilerinin aralığı `δω`, çözünürlük `dres = 2π/T`
biriminde (Hkeskin):

| τ kutusu | 0.3–0.4 | 0.4–0.5 | 0.5–0.6 | 0.6–0.7 | 0.7–0.8 | 0.8–0.9 | **0.9–1.0** |
|---|---|---|---|---|---|---|---|
| n_çizgi | 23 | 58 | 146 | 411 | 1167 | 3428 | **3730** |
| ⟨δω⟩/dres | 1442 | 524 | 205 | 72.9 | 25.8 | 8.76 | **4.02** |
| %(δω < dres) | 0 | 0 | 0 | 0 | 0 | 0.2% | **11.6%** |

Merdivenin üst %41'i pencerede ancak **4 çözünürlük genişliği** aralıkla
duruyor, %11.6'sı hiç ayrışmıyor. Gram köşegen-dışısının Dirichlet
kolunun (A1.5) yaşadığı yer burasıdır.

### A2.3 ★ SAYISAL OLARAK DEĞERLENDİRİLDİ — ve ×18 AZ VERİYOR

`168_rezonans.py Hkeskin 3`: `Ĵ_D` düzgün tarakta DTFT (OVS=8) ile TAM,
`S_exact = Σ_r 𝒢_r[Ĵ_D(W−ω_r)+Ĵ_D(W+ω_r)]` 8981 çizgi × 3 bant × 758
frekansta, 160'ın bant birleştirmesiyle `A²s2` birimine getirildi.
(GERÇEK-tarak `J_2` hiç hesaplanmadı — 167'de bekçiyi patlatan parça oydu;
karşılaştırma 165 §4b'nin kayıtlı Ç4'üne karşı yapıldı.)

| τ_eff | **S_exact = ★'nın Ç4'ü** | 165 §4b'nin ÖLÇÜLEN Ç4 (Hkeskin) | **★ / ölçülen** |
|---|---|---|---|
| 0.5393 | +0.015237 | +0.279046 | **0.055** |
| 0.5792 | +0.010729 | — | — |
| 0.6189 | +0.008913 | +0.609409 | **0.015** |

> **HÜKÜM (★ DÜŞTÜ).** Birinci-mertebe rezonans integrali ölçülen Ç4'ün
> **%1.5–5.5'ini** veriyor. A1.5'in kendi uyarısı (geçerlilik `τ ≤ 0.573`,
> Ç4 ise `τ_etkin ≈ 2–3`) **ölçüldü ve doğrulandı**: ★, Ç4 için
> geçerlilik penceresinin dışındadır.
>
> Not: 167'nin (R1) sınavında ★ `κ(ω_r)`'yi %0.2–6 içinde veriyordu —
> orada `ν = ω_r` TEK bir merdiven frekansıdır (`2πτσ_Ĉ ≤ 1.6`); burada
> `ν` DÖRT frekansın toplamıdır. **Aynı formül, farklı ν rejimi.**

### A2.4 YERLEŞİM ORANI ölçüldü — kararsız

A1.3'ün oranı doğrudan ölçüldü: `Ĵ_D`'nin `θ`-ızgarasında `±K·dres`
kutu ortalaması (≡ pencereyi `T/K`'ya kısaltmak) ve merdiven noktalarının
rastgele kaydırılması (`ω_r → ω_r + δ_r`, `δ_r ∈ ±3dres` / `±30dres`):

| τ_eff | exact | K=1 | K=2 | K=4 | K=8 | K=32 | shift±3 | shift±30 |
|---|---|---|---|---|---|---|---|---|
| 0.5393 | +0.015237 | **0.515** | 0.199 | 0.133 | 0.089 | −0.014 | 0.049 | −0.115 |
| 0.5792 | +0.010729 | **0.247** | 0.161 | 0.069 | −0.032 | −0.011 | 0.313 | 0.082 |
| 0.6189 | +0.008913 | **0.474** | 0.014 | −0.192 | −0.291 | −0.232 | 1.112 | −0.261 |

> Merdiven noktaları gerçekten ÖZELDİR (kaydırma toplamı 3–20 kat küçük ve
> işaret değiştiriyor) — yani A1.3'ün "ayrık örnekleme tepelerin üstüne
> düşüyor" resmi NİTEL olarak doğru. **Ama oranın kendisi banttan banda
> 0.25–0.52 arasında oynuyor ve K ≥ 2'de işaret değiştiriyor**; 167'nin
> açıklamamızı istediği %10–14'lük etkileri taşıyacak kararlılıkta DEĞİL.

### A2.5 İKİ SERT KISITLA YÜZLEŞME — A1'in İKİ ÖNGÖRÜSÜ DE DÜŞTÜ

`168_olcek.py` (ortak pencere `lo ∈ [0.52, 0.68]`, beş bant, log-ortalama):

| eksen | gaz | **c/c₀ ÖLÇÜLEN** | **A1 yerleşim öngörüsü** | fark |
|---|---|---|---|---|
| pencere | HkT2a (T/2) | 0.9566 | **2.00** | ×2.1, TERS YÖN |
| pencere | HkT2b (T/2) | 0.9220 | **2.00** | ×2.2, TERS YÖN |
| pencere | HkT4a (T/4) | 0.8170 | **4.00** | ×4.9, TERS YÖN |
| pencere | HkT4b (T/4) | 0.8188 | **4.00** | ×4.9, TERS YÖN |
| kesim | K090 (≤0.90) | 0.9340 | 0.778 (Σ𝒢²) / 0.669 (Σ\|𝒢\|) | ×1.20 / ×1.40 |
| kesim | K070 (≤0.70) | 0.8927 | **0.245** / **0.105** | ×3.6 / ×8.5 |
| kesim | HA4 (erfc 0.68) | 0.8766 | **0.207** / **0.125** | ×4.2 / ×7.0 |
| kesim | E060 (erfc 0.60) | 0.9034 | **0.174** / **0.063** | ×5.2 / ×14 |

> **HÜKÜM.** Yerleşim formu (`c = n_Δ·w_eff`) **iki eksende de düştü**:
> pencere ekseninde işaret bile ters (`c ∝ T^{−1}` öngörüyor, ölçülen
> `c ∝ T^{+0.06…+0.15}`), kesim ekseninde 3.6–14 kat fazla düşüş
> öngörüyor. **Saf geometrik sabit adaylarından sonra (167 H1), şimdi de
> BİRİNCİ-MERTEBE REZONANS İNTEGRALİ ölüyor.**

### A2.6 ÖLÇEK YASALARI HANGİ ÇARPANDAN GELİYOR — AYRIŞTIRMA

`c = KALİB_u2/(W_amp·W_X)` ve `KALİB_u2 = g_cal · θ`, `g_cal = g_E g_X²`
(165 §7'nin çizgi-uyum sızıntısı; `g_E = ⟨η E_mod⟩/⟨E_mod²⟩`).

| gaz | c | g_cal | ĉ = c/g_cal | **θ = KALİB/g_cal** | KALİB | W_amp·W_X | n/N |
|---|---|---|---|---|---|---|---|
| Hkeskin | 0.5774 | 0.2895 | 1.9947 | **0.8884** | 0.2572 | 0.4454 | 0.0299 |
| son | 0.5664 | 0.2971 | 1.9065 | **0.9142** | 0.2716 | 0.4795 | 0.0299 |
| L085 (λ=0.85) | 0.5171 | 0.2912 | 1.7758 | **0.8959** | 0.2609 | 0.5045 | 0.0299 |
| K090 | 0.5393 | 0.2877 | 1.8746 | **0.8067** | 0.2321 | 0.4303 | 0.0299 |
| K070 | 0.5154 | 0.3128 | 1.6476 | **0.7142** | 0.2234 | 0.4334 | 0.0299 |
| HA4 | 0.5061 | 0.3194 | 1.5849 | **0.7086** | 0.2263 | 0.4471 | 0.0299 |
| E060 | 0.5216 | 0.3433 | 1.5195 | **0.7059** | 0.2423 | 0.4646 | 0.0299 |
| HkT2a | 0.5524 | 0.2701 | 2.0452 | **0.9166** | 0.2476 | 0.4482 | 0.0580 |
| HkT2b | 0.5324 | 0.2652 | 2.0078 | **0.8873** | 0.2353 | 0.4419 | 0.0619 |
| HkT4a | 0.4717 | 0.2208 | 2.1365 | **0.9595** | 0.2119 | 0.4491 | 0.1141 |
| HkT4b | 0.4728 | 0.2157 | 2.1917 | **0.9832** | 0.2121 | 0.4486 | 0.1180 |

**Eksen eksen oranlar (referans Hkeskin):**

| eksen | gaz | c/c₀ | **g_cal/g₀** | **θ/θ₀** | (1+2n/N)₀/(1+2n/N) |
|---|---|---|---|---|---|
| pencere | HkT2a | 0.9566 | **0.9330** | 1.0317 | **0.9498** |
| pencere | HkT2b | 0.9220 | **0.9160** | 0.9987 | **0.9432** |
| pencere | HkT4a | 0.8170 | **0.7628** | 1.0800 | **0.8629** |
| pencere | HkT4b | 0.8188 | **0.7452** | 1.1066 | **0.8575** |
| kesim | K090 | 0.9340 | 0.9938 | **0.9080** | 1.0000 |
| kesim | K070 | 0.8927 | 1.0807 | **0.8038** | 1.0000 |
| kesim | HA4 | 0.8766 | 1.1032 | **0.7975** | 1.0000 |
| kesim | E060 | 0.9034 | 1.1860 | **0.7945** | 1.0000 |
| kesim | son | 0.9810 | 1.0264 | 1.0290 | 1.0000 |
| **λ** | **L085** | 0.8956 | **1.0060** | **1.0084** | 1.0000 |

> **HÜKÜM (c(T) REZONANS DEĞİL).** Pencere ekseninde etkinin TAMAMI
> `g_cal`'da: `g/g₀ = 0.933 / 0.916 / 0.763 / 0.745`, ve bunlar
> **örneklem-içi çizgi-uyum sızıntısının** `(1+2n_çizgi/N)^{−1}`
> yasasıyla (0.950 / 0.943 / 0.863 / 0.858) aynı işaret ve aynı
> mertebede. Geriye kalan `θ` pencereye neredeyse duyarsız:
> 1.032 / 0.999 / 1.080 / 1.107 — yani `c`'deki −%18'lik etki `θ`'da
> +%0…+%11'e düşüyor (ve İŞARETİ dönüyor).
> **167'nin "c(pencere): T/2 → −%10, T/4 → −%25" kısıtı bir rezonans-sayım
> etkisi DEĞİL, sonlu-örneklem uydurma sızıntısıdır.** (Aşırı düzeltme
> payı T/4'te +%8–11; açık borç.)
>
> **HÜKÜM (c(kesim) GERÇEK).** Kesim ekseninde `g_cal` ≈ sabit
> (0.99–1.10) ve etki bütünüyle `θ`'da kalıyor — hatta güçleniyor
> (0.934 → 0.908, 0.893 → 0.804). **Kesim ekseni gerçek merdiven fiziğidir.**

### A2.7 KESİM YASASI KAPANDI: θ = 1 − β·φ (β = 0.2175 ± 0.0054)

Kesim gazlarının merdiveni `τ_ust` (ya da erfc) ile kesilmiştir, ama MODEL
hâlâ `τ_c = 0.95`'e kadar 8981 çizgi uydurur. **BOŞ ÇİZGİ KESRİ**

```
φ = 1 − (1/N_r) Σ_q w_q ,   w_q = 1{τ_q ≤ τ_ust}  ya da  ½erfc((τ_q−τ_c)/Δ)
```

| gaz | kesim | φ | **θ/θ₀ ölçülen** | β = (1−θ/θ₀)/φ | 1 − 0.218φ | fark |
|---|---|---|---|---|---|---|
| Hkeskin | ≤1.00 | 0.0000 | 1.0000 | — | 1.0000 | 0.00% |
| son | gerçek | 0.0000 | 1.0290 | — | 1.0000 | −2.81% |
| L085 | ≤1.00 | 0.0000 | 1.0084 | — | 1.0000 | −0.83% |
| **K090** | ≤0.90 | 0.4153 | **0.9080** | **0.2216** | 0.9095 | **+0.16%** |
| **K070** | ≤0.70 | 0.9270 | **0.8038** | **0.2116** | 0.7979 | **−0.74%** |
| **HA4** | erfc 0.68 | 0.9098 | **0.7975** | **0.2225** | 0.8017 | **+0.51%** |
| **E060** ‡ | erfc 0.60 | 0.9595 | **0.7945** | **0.2141** | **0.7908** | **−0.47%** |
| HkT2a/b | ≤1.00 | 0.0000 | 1.032 / 0.999 | — | 1.0000 | −3.1 / +0.1% |
| HkT4a/b | ≤1.00 | 0.0000 | 1.080 / 1.107 | — | 1.0000 | −7.4 / −9.6% |

‡ = **ÖRNEKLEM-DIŞI**: β = 0.218 üç gazdan (K090/K070/HA4) kuruldu, E060
öngörüsü (0.7908) YAZILDI, ölçüm sonra alındı (§5).

> **β = 0.2175 ± 0.0054 (dört kesim gazı, %2.5 tutarlılık); yasa kesim
> gazlarını ±%0.75 içinde veriyor.**
>
> **YASANIN ŞEKLİ TÜRETİM İÇİN BELİRLEYİCİ: `φ` bir ÇİZGİ SAYISI kesridir,
> bir GENLİK kesri değil.** Aynı gazlar için genlik-ağırlıklı formlar
> (`Σ𝒢²`: 0.778 / 0.245 / 0.207 / 0.174; `Σ|𝒢|`: 0.669 / 0.105 / 0.125 /
> 0.063) ölçülenden (0.908 / 0.804 / 0.798 / 0.795) 3.6–14 kat uzaktır.
> **Rezonanslar ŞİDDETLERİYLE değil, SAYILARIYLA katkı veriyor — yani
> tepki DOYMUŞ.** E060 bunu keskinleştiriyor: `Σ|𝒢|` HA4'ün yarısı kadar
> (0.063 ↔ 0.125) ama `θ` neredeyse aynı (0.7945 ↔ 0.7975), çünkü çizgi
> sayıları yakın (φ = 0.960 ↔ 0.910).

### A2.8 MUTLAK DEĞER: c HANGİ ÜYENİN SAYISI?

`Hkeskin`, dokuz sağlıklı bantta, `log KALİB/üye` doğrusu:

| üye | **τ-eğimi** | artık rms | c (log-ort) |
|---|---|---|---|
| K0 (sabit) | −1.371 | 0.0148 | 0.2548 |
| W_amp | −0.205 | 0.0079 | 0.3673 |
| **W_X** | **+0.091** | **0.0071** | **0.4036** |
| √(W_amp W_X) | −0.057 | 0.0074 | 0.3850 |
| W_pos | +0.654 | 0.0080 | 0.4747 |
| W_amp·W_X | **+1.257** | 0.0102 | 0.5818 |

**Türetimden gerekçe.** 163, `J_k`'yı `m_n` sitesinde yazıp
`m_{n+1} = m_n + ḡ(1+X̃_n)` kaymasını `Γ` çekirdeğiyle taşımak zorunda
kalmıştı; 165 §2a bunu tek siteye geçerek kaldırdı. O kaymanın taşıyıcı
üstündeki TAM çarpanı

```
e^{−iω_Q (m_{n+1}−m_n)} = e^{−iω_Q ḡ(1+X̃_n)} = e^{−2πiτ_Q} · e^{−2πiτ_Q X̃_n}
⇒  ⟨e^{−2πiτ_Q X̃}⟩ = W_X(τ_Q)                        (TAM karakteristik fonksiyon)
```

yani **`W_X` tam olarak tek-gap site kaymasının Debye–Waller çarpanıdır ve
öngörüde BİR KEZ görünmelidir.** Buna karşılık `W_amp = ⟨cos(πτ·ds)⟩`
yarım-gap GENLİK modülasyonudur ve **zaten ÖLÇÜLEN çizgi genliklerinin
(`h'_q = 2⟨η e^{−iω_q s}⟩`) içindedir**; bir kez daha bölmek çift sayımdır.

Ölçüm bunu doğruluyor: `W_X` ile bölünce bant bağımlılığı **kayboluyor**
(eğim −1.371 → **+0.091**, artık rms bütün üyelerin en küçüğü), `W_amp` de
bölününce **ters yönde açılıyor** (+1.257).

> **HÜKÜM. c, `W_X` üyesinin sayısıdır:**
> ```
> c(Hkeskin) = 0.4051 ,  c(son) = 0.4053   (%0.05!)
> c(HA4) = 0.3345 , c(K070) = 0.3480 , c(K090) = 0.3709 , c(L085) = 0.3751 ,
> c(E060) = 0.3701
> ```
> `0.5647` (W_amp·W_X, 167 T4) aynı büyüklüğün `W_amp` ile bir kez fazla
> bölünmüş hâlidir; 167'nin 13/20'lik kapanışını bozmaz (tek küresel
> ölçek), ama TÜRETİMİN hedefi 0.405'tir.

### A2.9 ÜÇ BAĞIMSIZ DOYUM İMZASI (§0'ın hükmünün dayanağı)

1. **★ ×18–68 az veriyor** (A2.3) — doğrusal-`𝒢` tahmini çok küçük.
2. **Kesim yasası ÇİZGİ SAYISI ile gidiyor, `Σ𝒢²`/`Σ|𝒢|` ile değil**
   (A2.7) — katkı şiddetten bağımsız.
3. **λ ekseni neredeyse BOŞ** (A3): λ 1.00 → 0.85 (genlikte −%15)
   `KALİB_u2`'yi yalnız **+%1.4** değiştiriyor; A1.4(iii)'ün doğrusal
   öngörüsü **+%12** idi (×8.6 fazla).

Üçü de aynı şeyi söylüyor: `2πτσ_Ĉ ≫ 1` rejiminde tepki genlikte
DOYMUŞ. Rezonans integralinin doğru biçimi üstel (Gram/DW) formdur:
`G_{qq'} = κ(ω_q−ω_{q'})`, ve ölçülen karşılığı `g_E, g_X`'tir (A1.5).

---

## A3 — λ-ABLASYONU (167'nin açık borcu KAPANDI)

`z_L085.npy` hazırdı; ölçüm `167_configs/167_olcum.py L085 0.40 0.95 0 0`
ile tek koşuda alındı (5.1 dk, arka planda). `C_L085.json` yazıldı.

```
**c(L085) = 0.5305**  (W_amp·W_X üyesi, log-sd 0.0589, 7 bant, eğim +0.663)
**c(L085) = 0.3751**  (W_X üyesi — A2.8'in üyesi)
```

Ortak pencerede (lo ≤ 0.68), Hkeskin'e oranla:

| büyüklük | L085/Hkeskin | yorum |
|---|---|---|
| **KALİB_u2** | **1.014** | **λ-DEĞİŞMEZ (%1.4)** |
| g_cal | 1.006 | λ-değişmez |
| θ | 1.008 | λ-değişmez |
| W_amp·W_X | 1.133 | λ ile büyüyor (σ_ds, σ_X̃ küçüldü ⇒ DW zayıfladı) |
| c (W_amp·W_X üyesi) | 0.896 | yalnız W çarpanından |
| c (W_X üyesi) | 0.946 | yalnız W_X'ten |

> **HÜKÜM (λ ekseni).** 167'nin dördüncü ekseni **NULL çıktı**: ölçülen
> ham kalibrasyon λ'ya duyarsız (%1.4). `c`'nin λ ile düşmesi bir fizik
> değil, W çarpanlarının λ ile değişmesinin artefaktıdır — A2.8'in
> "W_amp çift sayım" hükmüyle aynı yönde ikinci bir kanıt.
>
> A1.4(iii)'ün doğrusal rezonans öngörüsü (+%12) **düştü** (ölçülen +%1.4)
> — A2.9'un üçüncü doyum imzası.

---

## 4. DENETİM ve DÜRÜSTLÜK NOTLARI

* **Yeni gaz kurulmadı.** Bütün sayılar `scratchpad/167`'nin hazır
  `C_*.json` / `z_*.npy` dosyalarından; tek yeni ÖLÇÜM L085 (ve E060,
  §5) — ikisi de 167'nin kendi zinciriyle (`167_olcum.py`), kopyalanmış
  kod yok.
* **İlk (P1) denemesi yanlıştı ve düzeltildi:** ham `Re κ` ile ölçülen
  `w_eff/dres = 0.027` bir MERKEZ FAZI artefaktıdır (`arg κ ≈ νs̄`,
  `s̄ = 1.05e6`, ızgara adımı `dres/6` ⇒ faz > 2π döner). Düzeltilmiş
  ölçüm (`168_profil.py`) `1.0396` verdi. Yanlış sayı `168_rezonans.py`'nin
  (P1) bloğunda hâlâ üretiliyor — **oradaki (P1) tablosu kullanılmamalıdır.**
* **★'nın 165'in Ç4'üyle karşılaştırması bant eşleşmelidir:** her ikisi de
  aynı `agg` (160'ın `A²s2` birleştirmesi) ile; `τ_eff` 0.5393 ve 0.6189
  bantları birebir örtüşüyor.
* **A2.7'nin yasası ÜÇ noktalı bir uyumla kuruldu** (K090/K070/HA4) ve
  katsayısı (0.218) TÜRETİLMEMİŞTİR; §5'teki E060 bunun ÖRNEKLEM-DIŞI
  sınavıdır (öngörü ölçümden önce yazıldı, %0.5 tuttu). Yasa dört noktada
  duruyor, ama hâlâ tek parametrelidir ve φ ∈ (0, 0.41) aralığında hiç
  veri yoktur.
* **E060'ın gaz-düzeyi `c`'si `nan`** (yüksek bantlarda `KALİB` işaret
  değiştiriyor, SNR < 3). Bütün 168 karşılaştırmaları ORTAK PENCEREDE
  (lo ≤ 0.68, SNR ≥ 3, R ≥ 0.98) yapıldığı için etkilenmiyor; E060 orada
  beş sağlıklı bant veriyor.
* `θ`'nın bant-içi τ-eğimi düzeltilmemiştir (`g_cal` bant-bağımsız bir
  gaz sayısı olduğu için `θ`, `KALİB`'in eğimini aynen taşır); bütün
  eksen karşılaştırmaları AYNI bant kümesinde (lo ∈ [0.52, 0.68]) yapıldı.
* `Σ𝒢² ≈ ∫τcos²(πτ)dτ` yaklaşımı sayısal merdivenle %2.5 (τ_c=0.95) —
  %14 (τ_c=0.60) arası uyuşuyor (fark asal-kuvvet `m ≥ 2` terimleri ve
  π(x) ≈ x/log x yaklaşımı); kesim ORANLARINDA her yerde SAYISAL merdiven
  kullanıldı, kapalı integral yalnız yasanın şeklini göstermek için.

---

## 5. ÖRNEKLEM-DIŞI SINAV — E060 (kesim yasasının tek gerçek testi)

A2.7'nin `β = 0.218` katsayısı ÜÇ gazdan (K090, K070, HA4) kuruldu.
Dördüncü kesim gazı `E060` (erfc 0.60/0.125) o sırada ölçülmemişti
(`z_E060.npy` hazırdı, `log_C_E060.txt` boştu — 167'nin kesilmiş
koşusundan). Yasanın öngörüsü ÖLÇÜMDEN ÖNCE yazıldı:

```
φ(E060) = 1 − Σw_q/N_r = 0.9595   ⇒   θ/θ₀ = 1 − 0.218·0.9595 = 0.7908
```

Ölçüm (`167_olcum.py E060 0.40 0.95 0 0`, 5.0 dk, arka planda):

```
θ/θ₀ (E060) = 0.7945                    ⇒  ÖNGÖRÜ − ÖLÇÜM = −0.47%
```

> **HÜKÜM. Kesim yasası örneklem-dışı GEÇTİ (%0.5).** Karşılaştırma için
> aynı noktada A1'in genlik-ağırlıklı öngörüleri **0.174** (`Σ𝒢²`) ve
> **0.063** (`Σ|𝒢|`) — ölçülenin 4.6 ve 12.6 katı uzağında.
>
> E060 ayrıca `g_cal` ekseninin ayrı yürüdüğünü gösteriyor: `g/g₀ = 1.186`
> (bütün gazların en büyüğü) iken `θ` yasanın üstünde duruyor. Yani
> ayrıştırma (`KALİB = g_cal · θ`) sadece tanım değil, **iki çarpanın
> bağımsız değiştiği ölçülmüş bir ayrıştırma.**

---

## 6. AÇIK BORÇLAR (168'den sonra)

1. **β = 0.2175'in türetimi.** Çizgi-sayısı-orantılı kesim yasasının
   katsayısı; doymuş Gram resminde `h'†Gh'` köşegen-dışı sayımından
   çıkması beklenir. (Yasanın KENDİSİ örneklem-dışı geçti — §5; eksik
   olan yalnız katsayının kalemle çıkarılması.)
2. **T/4 dilimlerinde θ'nın +%8–11'lik artığı** (`(1+2n/N)^{−1}`
   düzeltmesi tam değil).
3. **Doymuş rezonans integralinin kapalı formu**: `κ(ω_q−ω_{q'})` Gram
   matrisi üzerinden `1/g_E`, `1/g_X`'in analitik hesabı — 168'in
   bıraktığı asıl T3 borcu.
4. L070 / L115 gazları inşa edilmedi (λ ekseninde tek nokta var).
