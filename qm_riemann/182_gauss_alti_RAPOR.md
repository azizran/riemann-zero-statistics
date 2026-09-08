# 182 — GAUSS-ALTININ ANATOMİSİ + VF AİLESİ BÜYÜTMESİ (8 Eylül 2026)

*(ARTIMLI yazıldı, K3 ile tamamlandı. Her sayının yanında onu üreten
betiğin adı vardır. Sonuçlar YALNIZ gerçek koşudan; eşik gevşetme yok;
ölümler ve hükümsüzlükler kurtarmasız.)*

## Sorulan soru (180-K2'nin sahipsiz keşfi)

Fazları karıştırılmış (KİLİTSİZ) vekil gaz, 3-bacak kırpma aktarımında
Gauss referansına çökmüyor — **%10.7 ALTINA** iniyor. Hareketin tamamı
E-bacağında (`ρ(E)`: `son` 0.8229 → karışık 0.637); Xa/Xb karıştırmaya
kör. Aynı vekillerde `R_η = 0.809 < 1`: karışık gaz **anti-girişimli**.
*Kilitsiz bir gaz Gauss'tan neden ve nasıl sapar?*

---

## 0. K0 — DONMUŞ ÖN-KAYIT (`182a_onkayit.py`)

**Mühür:** `Tue Sep 8 16:57:40 +03 2026`,
betik sha256 `a221cc242f55daa2ed5bc6e535088884228a166a8d23af8d389952785a8f3bbe`,
KALEM sha256 `e2c2c91dd9e9e8c8ea5b2dbccab1b6ec70f219d5b4138a8a8bf6cab450c1a7d4`.
Dosya: `182/ONKAYIT_182.json`. Bu betik koştuğunda VF5..VF8 daha
**kurulmamıştı** (inşa 16:57'de, ön-kayıt diske yazıldıktan sonra
başladı — `182/insa_durum.txt` zaman damgaları).

### 0.1 Devralınan eşikler — GEVŞETME YOK

| eşik | değer | kaynak |
|---|---|---|
| `GAUSS1 = √(2/π)` | 0.7978845608028654 | 180a (bit-bit devralındı) |
| `GAUSS3 = (2/π)^{3/2}` | 0.5079490874739278 | 180a (bit-bit devralındı) |
| `ρ₃(Hkeskin)` / `ρ₃(son)` | 0.5284 / 0.5356 | 180a |
| 180e ÇÖZÜNÜRLÜK eşiği | **0.0102254562630361** | 180a, `0.5·|ρ₃(Hk)−GAUSS3|` |
| 181a okunurluk O1 | **`|ZARF| ≥ 2·se(ZARF)`** | 181a, `k_sigma = 2.0` |
| inşa kapıları G1..G5 | sha(A) bit-bit, φ≡0 → 0.0/0.0, ilk-kök 300000/300000, maks\|F\| ≤ 1e−8 (aşan 0), sıralılık TAM | 176a/180a |

> `182a` eşiği yeniden türetmiyor, **devralıyor**: `√(2/π)³` ile 180'in
> `GAUSS3`'ü son bitte 1 ULP farklı (0.5079490874739279 ↔ …78). 180'in
> sayısı esas alındı ve çözünürlük eşiği **0.0e+00 farkla** yeniden
> üretildi (`182a` çıktısı). Bu bir denetimdir, seçim değil.

### 0.2 Gauss referansının TÜRETİM VARSAYIMLARI (H-G1'in sınav nesnesi)

Gauss referansı bir ölçüm değil bir **türetimdir**: `V1 ∧ V2 ∧ V3 ∧ V4`.

* **V1 (marjinal):** her bacak `F` ayrı ayrı Gauss ⇒
  `κ_F := ⟨clip(F)·F⟩/⟨F²⟩ = E|F|/σ_F = √(2/π)`.
  *(169_k2b bu sayıyı zaten basıyor — sınav bedava.)*
* **V2 (ortak Gaussluk):** `E, X_a, X_b` **birlikte** Gauss ⇒ Price
  teoremi: `⟨sgn(F)·G⟩ = √(2/π)·⟨F·G⟩/σ_F` her `G` için; kırpma çarpımın
  her çapraz momentini aynı çarpanla ölçekler.
* **V3 (çarpımsallık):** bacaklar **bağımsız** ⇒ k bacakta `(√(2/π))^k`.
* **V4 (bant yansızlığı):** `P`'nin bant ağırlıklandırması
  (`pN·A_on² − pO·A_off²`) ve `s_den_J` zinciri, kırpmanın yol açtığı
  spektral yeniden dağılıma kör değildir; referans bu ağırlığın ρ'yu
  ölçeklemediğini varsayar.

### 0.3 Ölüm koşulları (sayısal, veriden önce)

| hipotez | ÖLÜM | YAŞAMA |
|---|---|---|
| **H-G1** (referansın türetimi) | 8/8 tohumda `\|κ_E/√(2/π)−1\| ≤ 0.05` iken açık `(√(2/π)−ρ(E))/√(2/π) ≥ 0.10` **VE** `\|⟨ρ°(E)⟩/√(2/π)−1\| ≤ 0.02` | `\|⟨ρ°(E)⟩/√(2/π)−1\| > 0.05`, ya da `κ_E`'nin sapması açığın en az yarısını tek başına açıklıyor |
| **H-G2** (inşa öz-tutarlılığı) | `⟨korr(X_a,X_b)⟩` sıfırla uyumlu (≤2se) **veya** işaret 8/8 tekdüze değil, **veya** `Rec ≤ 0.30` | `\|⟨r⟩\| ≥ 5·se` ∧ işaret 8/8 tekdüze ∧ `Rec ≥ 0.70` (≥7/8 tohumda) |
| **H-G3** (kestirimci artefaktı) | 8/8 tohumda `\|ρ_J(E)/ρ(E)−1\| ≤ 0.05` **ve** `ρ_J`'nin Gauss açığı ≥ 0.10 | `⟨\|ρ_J/ρ−1\|⟩ > 0.10` |

Aralarda kalan bölgeler için **HÜKÜMSÜZ** yazılır (kurtarma yok).

Tanımlar (182a'da donduruldu):

```
ρ°(E)  : ÜÇ bacak da faz-rastgeleleştirilmiş surrogate (güç tayfı bit-bit
         korunur, fazlar bağımsız U(0,2π)) — referansın KENDİ varsayımları
         altında sınavı.
ρ*(E)  : YALNIZ E surrogate (X_a, X_b aynen)  ⇒
         Rec := (ρ*(E) − ρ_J(E)) / (√(2/π) − ρ_J(E))   [geri kazanım]
ρ_J    : BAĞIMSIZ KESTİRİMCİ, yalnız model tarafı —
         ρ_J(var) = Σ_k Re[J_var(W_k)·conj(J_000(W_k))] / Σ_k |J_000(W_k)|²
         (ölçülen gazın hiçbir niceliği girmez; bant ağırlığı yok;
          s_den_J kullanılmaz; W_k kümesi 169_k2b'nin Wall'ı, aynen)
ρ(E) ≡ κ_E · β_E ;  PAY_ÇİZGİ := (log κ_E − log√(2/π))/(log ρ(E) − log√(2/π))
```

### 0.4 ÖN-KAYITLI ÖNGÖRÜ P1 — kırpma hakemi n=8'de KURTULAMAZ (yapısal)

180e'nin çözünürlük istatistiği **RANGE**'dir
(`saçılım := max ρ₃ − min ρ₃`) ve range **n'de monoton artandır**.
`n = 4`'te ölçülen 0.016091 zaten eşiğin (0.010225) 1.57 katıydı.
Dolayısıyla:

> **VF5..VF8 eklenince saçılım ≥ 0.016091 KESİNDİR ⇒ HÜKÜMSÜZ hükmü
> n = 8'de de değişemez.** Bu, ölçümden önce yazılmış yapısal bir
> öngörüdür ve şu anlama gelir: *180'in çözünürlük maddesi tohum
> EKLEYEREK aşılamaz.* Eşik gevşetilmiyor; hüküm aynen bırakılıyor.
> KALEM'in "(b) saçılım eşiğin altına indi mi" sorusunun yanıtı, veriden
> önce, **hayır — ve inemez**dir.

### 0.5 ÖN-KAYITLI ÖNGÖRÜ P2 — O1'in n_VF = 8 hâli (`182a`)

`se_VF(8) = sd_VF(4)/√8` varsayımıyla (merkezler sabit):

| satır | konv. | `|Z|/se` (n=4) | **öngörü (n=8)** | öngörülen hüküm |
|---|---|---|---|---|
| M | doğ./log | 10.34 / 10.26 | 12.45 / 12.45 | OKUNUR |
| Q_E | doğ./log | 9.95 / 9.90 | 10.93 / 10.83 | OKUNUR |
| ρ_E | doğ./log | 5.76 / 5.78 | 6.32 / 6.35 | OKUNUR |
| g_E | doğ./log | 28.76 / 28.75 | 33.23 / 33.31 | OKUNUR |
| **Q_X** | doğ./log | 0.42 / 0.41 | 0.49 / 0.47 | HÜKÜMSÜZ |
| **ρ_X** | doğ./log | 0.51 / 0.51 | 0.59 / 0.59 | HÜKÜMSÜZ |
| **θ** | doğ./log | 0.88 / 0.87 | 1.03 / 1.03 | HÜKÜMSÜZ |
| **g_X** | doğ./log | 0.84 / 0.84 | 0.97 / 0.97 | HÜKÜMSÜZ |
| **g_cal** | doğ./log | 1.69 / 1.70 | 1.94 / 1.97 | HÜKÜMSÜZ |

> **ÖNGÖRÜLEN YENİ OKUNUR SATIR: YOK.** Üstelik `Q_X`, `ρ_X`, `θ`, `g_X`
> için ön-kayıt **yapısal** bir gerekçe de yazdı: `se_ZARF ≥ se_VS` her
> zaman, ve bu dört satırda `|ZARF| < 2·se_VS`. Yani **VF tohumu
> ekleyerek de ULAŞILAMAZ**lar — 181 onları "VF'nin saçılımı boğuyor"
> diye HÜKÜMSÜZ bırakmıştı; 182'nin ön-kaydı bunun **yarısının doğru
> olmadığını** ölçümden önce söylüyor: VF tarafı sıfıra gitse bile
> `Q_X, ρ_X, θ, g_X` açılmaz, çünkü onları **VS tarafı da** boğuyor.
> `g_cal` tek "yakın" satır (1.94 öngörü, eşik 2.0) ama yine HÜKÜMSÜZ.

---

## 1. K1 — VF5..VF8 İNŞASI (`176b_vekil_insa.py`, AYNEN)

Makine 176'nınkinin **aynısı**: zarf `Hkeskin`'in genlik dizisi
(`a_q = 1/(πk√q)`, `τ_q ≤ 1.00`, 15450 çizgi), her çizgiye bağımsız
`φ_q ~ U(0,2π)`, sonra merdiven `164_insa.coz_sadakatli` ile
(ızgara braketi + SIRALI ilk-kök + korumalı Newton, `h = 0.015`,
`nz = 300000`, `c = −½`) sıfırdan çözülür. Yeni satır yazılmadı;
betik `VF5 5`, `VF6 6`, `VF7 7`, `VF8 8` ile çağrıldı.

### T1 — İnşa kapıları (eşikler 176a/180a'dan; gevşetme yok)

**ZARFIN BİT-BİT SHA KANITI** — dört tohumda da:

```
sha256(A)  vekil   = 45af6fa572bdd99c3a39c1ff75ccb8693038f733bf32d4ef7f921cae9e43c267
sha256(A) Hkeskin  = 45af6fa572bdd99c3a39c1ff75ccb8693038f733bf32d4ef7f921cae9e43c267
maks|A(vekil) − A(Hkeskin)| = 0.0e+00     maks|ω − ω| = 0.0e+00
```

| gaz | tohum | G1 sha(A) | G2 φ≡0 (ΔS/ΔS′) | G3 ilk-kök | G4 maks\|F\| (aşan) | G5 sıralılık (min Δz) | σ_ds | L | süre |
|---|---|---|---|---|---|---|---|---|---|
| `VF5` | 5 | ✓ eşit | 0.0e+00 / 0.0e+00 | **300000/300000** | 1.863e−09 (0) | **TAM** (0.093720) | 0.44362 | 12.029593224 | 9.8 dk |
| `VF6` | 6 | ✓ eşit | 0.0e+00 / 0.0e+00 | **300000/300000** | 1.863e−09 (0) | **TAM** (0.092361) | 0.44150 | 12.029593224 | 11.4 dk |
| `VF7` | 7 | ✓ eşit | 0.0e+00 / 0.0e+00 | **300000/300000** | 1.863e−09 (0) | **TAM** (0.094447) | 0.44259 | 12.029593224 | 11.0 dk |
| `VF8` | 8 | ✓ eşit | 0.0e+00 / 0.0e+00 | **300000/300000** | 1.863e−09 (0) | **TAM** (0.091496) | 0.44190 | 12.029593224 | 11.6 dk |

> **DÖRT TOHUMDA DA BEŞ KAPININ HEPSİ GEÇTİ** (`KAPILAR: G1..G5 = True`).
> `maks|F|` dördünde de birebir aynı (1.863e−09): çözücü tohumdan
> bağımsız olarak aynı yakınsama tabanına oturuyor.
>
> Karşılaştırma: `VF1..VF4`'ün σ_ds'i 0.4416/0.4427; `Hkeskin`
> L = 12.029593242. Yeni dört tohum aynı ailede.

---

### T1b — Ölçüm zinciri (`182c_olcum.py` → `176c_olcum.kos` → `167_olcum.kos(ad, 0.40, 0.95, 0, 0)` AYNEN)

Özdeşlik denetimi `|ΔK|/K` (172b'nin Ö-A özdeşliği): E 7.4e−16 … 5.8e−15,
X 7.8e−16 … 4.2e−15 — dört tohumda da makine sıfırında.

**176a'nın `F0` alt-kapısı (`R_bant ≥ 0.98`, lo ∈ [0.52,0.68]) DÖRT YENİ
TOHUMDA DA ISKALADI** ve kurtarılmadı:

| gaz | `R_bant` (lo = 0.52 … 0.68) | min |
|---|---|---|
| `VF5` | 0.9793 0.9826 0.9563 0.9181 0.8975 | **0.8975** ✗ |
| `VF6` | 0.9959 0.9827 0.9660 0.9404 0.9159 | **0.9159** ✗ |
| `VF7` | 0.9957 0.9833 0.9521 0.9355 0.8972 | **0.8972** ✗ |
| `VF8` | 0.9908 0.9763 0.9612 0.9331 0.9051 | **0.9051** ✗ |

> Bu bir **devralınan ıskadır**, yeni bir kırılma değil: `VF1..VF4`'ün
> min'leri de 0.9125 / 0.8925 / 0.9182 / 0.9141 idi ve 176 bunu §2c'de
> "ıskanın kendisi bir ölçümdür" diye kaydetmişti. Yeni dört tohum
> **aynı aralıkta**. Eşik gevşetilmedi; ıska aynen yazılıyor. (Bu kapı
> 176'nın F0 yüzleşmesini kapatır; defter satırını değil.)

---

## 2. K1(a) — KIRPMA HAKEMİ, n = 8 (`182e_hakem8.py`)

**ÜREME KAPISI:** aynı kod `VF1..VF4` ile koşulduğunda
`180/K2_HAKEM.json`'un 8 karar alanını (`dal, ⟨ρ₃⟩, sd, se, C, z_G,
saçılım, eşik`) **0 farkla** üretti ⇒ makine 180e'ninkinin bit-bit aynısı.

### T5 — 3-bacak kırpma aktarımı, n = 8 (`169_k2b.main(g, 0.40, 0.95)` aynen)

| gaz | ρ(E) | ρ(Xa) | ρ(Xb) | ρ(XaXb) | ρ(EXa) | **ρ₃** |
|---|---|---|---|---|---|---|
| `VF1` | 0.6487 | 0.8290 | 0.8179 | 0.6693 | 0.5467 | 0.4474 |
| `VF2` | 0.6374 | 0.8283 | 0.8163 | 0.6692 | 0.5362 | 0.4477 |
| `VF3` | 0.6452 | 0.8244 | 0.8365 | 0.6778 | 0.5427 | 0.4635 |
| `VF4` | 0.6529 | 0.8206 | 0.8095 | 0.6575 | 0.5530 | 0.4553 |
| **`VF5`** | 0.6360 | 0.8279 | 0.8283 | 0.6789 | 0.5403 | **0.4476** |
| **`VF6`** | 0.6443 | 0.8228 | 0.8213 | 0.6731 | 0.5480 | **0.4620** |
| **`VF7`** | 0.6424 | 0.8256 | 0.8199 | 0.6715 | 0.5453 | **0.4506** |
| **`VF8`** | 0.6460 | 0.8214 | 0.8174 | 0.6644 | 0.5372 | **0.4500** |
| **GAUSS** | 0.7979 | 0.7979 | 0.7979 | 0.6366 | 0.6366 | **0.5079** |
| *(180-T5)* `Hkeskin` | 0.7947 | 0.8463 | 0.8468 | 0.6220 | 0.6952 | *0.5284* |
| *(180-T5)* `son` | 0.8229 | 0.8366 | 0.8367 | 0.6023 | 0.7147 | *0.5356* |

VF ailesinin bacak bacak n=8 özeti (ort ± sd):

| bacak | ort ± sd | Gauss | fark |
|---|---|---|---|
| ρ(E) | **0.6441 ± 0.0056** | 0.7979 | **−19.27%** |
| ρ(Xa) | 0.8250 ± 0.0032 | 0.7979 | +3.40% |
| ρ(Xb) | 0.8209 ± 0.0082 | 0.7979 | +2.88% |
| ρ(XaXb) | 0.6702 ± 0.0070 | 0.6366 | +5.27% |
| ρ(EXa) | 0.5437 ± 0.0057 | 0.6366 | −14.60% |
| **ρ₃** | **0.4530 ± 0.0065** | 0.5079 | **−10.82%** |

```
⟨ρ₃⟩_VF = 0.453014   sd = 0.006544   se = 0.002314   (n = 8)
C   = +3.6862        z_G = −23.745
saçılım = 0.016091   (ön-kayıtlı HÜKÜMSÜZ eşiği 0.010225)
```

### (a) GAUSS-ALTI HÜKMÜ: **VERİLDİ** (`182a`/KARAR.A, n = 8)

> `z_G = −23.7` (eşik ≤ −5 ✓) ve **8/8 tohumun her biri tek başına**
> `ρ₃ < GAUSS3` ✓ ⇒ **KARIŞIK (KİLİTSİZ) VEKİL GAZ GAUSS'UN ALTINDADIR.**
> Açık `−%10.82`; `n = 4`'teki `−%10.7` ile birebir aynı yerde, ama
> hata çubuğu yarıya indi (se 0.003807 → 0.002314) ve z 14.3'ten
> **23.7**'ye çıktı. Bu bir eğilim değil, **ailenin özelliği**.

### (b) KIRPMA HAKEMİ: **HÜKÜMSÜZ — AYNEN** (P1 tuttu)

> saçılım(n=8) = **0.016091** > eşik 0.010225 ⇒ ön-kayıtlı dal
> **HÜKÜMSÜZ**. Ön-kayıtta yazılan P1 öngörüsü doğrulandı
> (`saçılım(8) ≥ saçılım(4)` ✓). **Üstelik tam eşitlikle:** dört yeni
> tohumun hiçbiri `[0.4474, 0.4635]` aralığının dışına çıkmadı, yani
> range **birebir aynı kaldı**. Saçılım azalmadı, azalamazdı ve
> azalmadı.
>
> **Bu bir başarısızlık değil, bir sınırın ölçümüdür:** 180'in
> çözünürlük maddesi bir RANGE üstüne kurulmuştu; range tohum sayısıyla
> monotondur. Dolayısıyla **H-180a hakemi, VF ailesini büyüterek
> KAZANILAMAZ** — hüküm 180'de nasıl bırakıldıysa öyle kalır.
> (Eşiği ölçek-değişmez bir istatistiğe — sd ya da se'ye — çevirmek
> hükmü değiştirirdi; bu bir **eşik gevşetmesi** olurdu ve
> YAPILMAMIŞTIR. `sd = 0.006544`, `se = 0.002314` yalnız kayıt için.)
>
> `DAL B`'nin 180'deki ölümü n = 8'de de sağlam: `C = +3.69`, koşulu
> `C ≤ 0.30` idi; sekiz tohumun her biri tek başına eşiğin 10 katı
> ötesinde.

---

## 3. K2 — E-BACAĞININ ANATOMİSİ (`182g_anatomi.py`, `182i_hipotez.py`)

`182g`, `169_k2b.main`'in ilk yarısını **birebir** yeniden kurar
(`K.gaz` → `Model165(...,'olculen',τ_c=0.95)` → `.alanlar()` →
τ-sıralı bir-atlamalı bölünme). Doğrulama: yeniden kurulan zincir
`VF1` için `korr(Xa,Xb) = +0.10440` ve `κ_E = 0.80990` verdi —
`180d`'nin logundaki `+0.1044` / `0.80990` ile **aynı**.

### 3.1 T-A1 — ÇİZGİ / BANT AYRIŞIMI: açığın tamamı BANT'ta

`ρ(E) ≡ κ_E · β_E`, `κ_E = ⟨clip(E)·E⟩/⟨E²⟩ = E|E|/σ_E`.

| gaz | κ_E | κ_Xa | κ_Xb | kurt(E) | ρ(E) | β_E | PAY_ÇİZGİ |
|---|---|---|---|---|---|---|---|
| VF1 | 0.80990 | 0.80401 | 0.80576 | −0.26207 | 0.64869 | 0.80095 | −0.0722 |
| VF2 | 0.81016 | 0.80471 | 0.80520 | −0.26153 | 0.63738 | 0.78673 | −0.0680 |
| VF3 | 0.81013 | 0.80435 | 0.80615 | −0.25935 | 0.64515 | 0.79636 | −0.0717 |
| VF4 | 0.80998 | 0.80480 | 0.80603 | −0.26632 | 0.65286 | 0.80603 | −0.0750 |
| VF5 | 0.81007 | 0.80413 | 0.80512 | −0.26301 | 0.63604 | 0.78516 | −0.0669 |
| VF6 | 0.81023 | 0.80440 | 0.80625 | −0.26250 | 0.64425 | 0.79515 | −0.0718 |
| VF7 | 0.81009 | 0.80424 | 0.80446 | −0.26991 | 0.64239 | 0.79298 | −0.0700 |
| VF8 | 0.81019 | 0.80423 | 0.80538 | −0.26265 | 0.64605 | 0.79740 | −0.0725 |
| **⟨⟩ ± sd** | **0.810093 ± 0.000110** | | | **−0.26342 ± 0.00325** | **0.644101 ± 0.005553** | **0.795096 ± 0.006906** | **−0.0710 ± 0.0026** |

> **`PAY_ÇİZGİ = −0.071`.** Yani ρ(E)'nin Gauss açığının **%107'si
> BANT'tan, %−7'si çizgiden** geliyor: alanın marjinal biçim çarpanı
> κ_E, Gauss'un **üstünde** (0.8101 vs 0.7979, +%1.5) — açığı
> kapatmaya çalışıyor, açmaya değil. E alanı **basıklığı eksik**
> (kurt = −0.263, yani platikurtik), ama bu Gauss'tan sapmanın
> **işareti yanlış ve büyüklüğü 6 kat küçük**.
>
> **ρ(E)'nin %19'luk açığı bir genlik-dağılımı olgusu DEĞİLDİR.**

### 3.2 T-A2 — τ eğimleri: E bacağı bant-bağımlı, X bacakları değil

```
⟨s_E⟩  = +0.38243 ± 0.09330   (se 0.03299)
⟨s_Xa⟩ = +0.06920 ± 0.12943
⟨s_Xb⟩ = +0.04434 ± 0.11629
ÖN-KAYIT: |⟨s_E⟩| ≥ ⟨|s_X|⟩ + 3·se(s_E) ?   0.38243 ≥ 0.21615   ⇒ ✓
```

> `ρ(E)` beş bantta 0.60 → 0.69'a **tırmanıyor**; `ρ(Xa)`, `ρ(Xb)` düz.
> Açık en dar bantta (τ_eff ≈ 0.54) en derin. **Kırpma açığı τ'ya
> bağlıdır ⇒ spektral (BANT) karakterlidir.**

### 3.3 T-A3 — Surrogate merdiveni ve bağımsız kestirimci

| gaz | ρ(E) (169) | ρ_J(E) (bağımsız) | \|ρ_J/ρ−1\| | ρ°(E) (üçü de surr.) | ρ*(E) (yalnız E surr.) | Rec |
|---|---|---|---|---|---|---|
| VF1 | 0.64869 | 0.67774 | 0.04479 | 0.78584 | 0.80149 | +1.0300 |
| VF2 | 0.63738 | 0.66922 | 0.04996 | 0.80108 | 0.79536 | +0.9804 |
| VF3 | 0.64515 | 0.67491 | 0.04613 | 0.79761 | 0.79397 | +0.9681 |
| VF4 | 0.65286 | 0.68124 | 0.04347 | 0.80074 | 0.79253 | +0.9541 |
| VF5 | 0.63604 | 0.66738 | 0.04928 | 0.80890 | 0.79964 | +1.0135 |
| VF6 | 0.64425 | 0.67506 | 0.04781 | 0.79086 | 0.79778 | +0.9991 |
| VF7 | 0.64239 | 0.67673 | **0.05346** | 0.80431 | 0.79690 | +0.9919 |
| VF8 | 0.64605 | 0.67652 | 0.04716 | 0.76855 | 0.80179 | +1.0322 |
| **⟨⟩** | 0.6441 | 0.6743 | 0.04776 | **0.79474** | **0.79743** | **+0.9962** |

### 3.4 T-A4 — İnşa öz-tutarlılığı: bağımsızlık ÜRETİLMİYOR

| gaz | korr(Xa,Xb) | korr(Xa°,Xb°) | R_η | Kov(η_artık,η)/Var(η) ≡ 1−R_η | Var(çizgi)/P |
|---|---|---|---|---|---|
| VF1 | +0.10440 | −0.00266 | 0.80989 | +0.19011 | 1.32170 |
| VF2 | +0.09000 | −0.00054 | 0.80710 | +0.19290 | 1.32221 |
| VF3 | +0.10391 | +0.00430 | 0.81118 | +0.18882 | 1.32388 |
| VF4 | +0.10458 | +0.00774 | 0.81539 | +0.18461 | 1.32122 |
| VF5 | +0.09820 | −0.00790 | 0.80372 | +0.19628 | 1.31983 |
| VF6 | +0.10713 | −0.00801 | 0.81244 | +0.18756 | 1.32129 |
| VF7 | +0.09607 | +0.00440 | 0.80962 | +0.19038 | 1.32168 |
| VF8 | +0.10016 | −0.00078 | 0.81291 | +0.18709 | 1.32297 |
| **⟨⟩ ± sd** | **+0.10056 ± 0.00564** | **−0.00043 ± 0.00574** | **0.81028 ± 0.00363** | **+0.18972** | **1.32185** |

*(R_η: `176d_K1.py` → `174b_K1_girisim.kos`, 162'nin makinesi aynen;
VF3..VF8 bu görevde ilk kez ölçüldü.)*

### 3.5 ÖN-KAYITLI HÜKÜMLER (`182i_hipotez.py`)

#### H-G1 (Gauss referansının türetimi): **ÖLDÜ** — kurtarma yok

```
(i)  maks|κ_E/√(2/π) − 1| = 0.01547 ≤ 0.05          ✓
     min açık (√(2/π) − ρ(E))/√(2/π) = 0.18176 ≥ 0.10 ✓
(ii) ⟨ρ°(E)⟩ = 0.79474 ,  |⟨ρ°(E)⟩/√(2/π) − 1| = 0.00394 ≤ 0.02 ✓
     κ_E'nin açığın yarısını taşıması: 0.01234 vs 0.07251 ⇒ TAŞIMIYOR
```

> Üç bacak da faz-rastgeleleştirilmiş surrogate'la (güç tayfı bit-bit
> korunarak) değiştirilince **kestirimci Gauss referansını %0.4 içinde
> geri veriyor** (0.79474 ↔ 0.79788). Yani **referansın türetimi
> doğrudur ve varsayımları altında bu makineden doğru çıkar.**
> Açığı referans taşımıyor. **V1 (marjinal Gaussluk) sağlanıyor,
> V4 (bant yansızlığı) sağlanıyor.** Geriye V2/V3 kalıyor.

#### H-G3 (kestirimci artefaktı): **HÜKÜMSÜZ** — gevşetme yok

```
maks|ρ_J(E)/ρ(E) − 1| = 0.05346  ≤ 0.05 ?  ✗  (ort 0.04776; 7/8 tohum geçti)
bağımsız kestirimcinin Gauss açığı: min 0.14619 ≥ 0.10  ✓
```

> Ölüm koşulu **8/8** istiyordu; `VF7` 0.05346 ile eşiği (0.05) 0.0035
> aştı. **Eşik gevşetilmedi ⇒ H-G3 ÖLMEDİ, HÜKÜMSÜZ yazıldı.**
> Kayda geçen: ölçülen sapma tek yönlü ve küçük (`ρ_J` her tohumda
> `ρ`'dan %4.3–5.3 **büyük**), ve **bağımsız kestirimci de Gauss'un
> %14.6–%16.4 altındadır** — yani açığın kendisi kestirimciye bağlı
> değildir; hüküm bir **eşik meselesidir**, olgu meselesi değil.
> *(Bunu ölüme çevirecek olan şey bir eşik değişikliğidir ve
> yapılmamıştır; ileride ölçülecekse ön-kayıt yeniden yazılmalıdır.)*

#### H-G2 (inşa öz-tutarlılığı): **YAŞIYOR**

```
⟨korr(Xa,Xb)⟩ = +0.10056 ± 0.00564 (sd), se = 0.00199  ⇒ |⟨r⟩|/se = 50.4  (≥5 ✓)
işaret 8/8 tekdüze ✓            sıfırla uyumlu ? HAYIR ✓
surrogate kontrolü ⟨korr(Xa°,Xb°)⟩ = −0.00043 ± 0.00574   (sıfır)
⟨Rec⟩ = +0.9962 ,  Rec ≥ 0.70 olan tohum: 8/8  ✓
```

> **İki bağımsız parmak izi, aynı yeri gösteriyor:**
>
> 1. **Merdivenin iki AYRIK yarısı bağımsız değil.** `X_a` ve `X_b`
>    çizgileri paylaşmaz (τ'ya göre sıralı, bir atlamalı bölünme), ama
>    ölçüm sitelerinde sentezlenmiş alanları `+0.1006 ± 0.0056`
>    ilişkilidir — sıfırdan **50 σ** uzakta. Aynı genlik tayfının
>    faz-rastgeleleştirilmiş surrogate'ında bu ilişki **−0.0004**'e
>    (yani sıfıra) düşüyor. Demek ki bağımlılığı üreten şey genlik
>    tayfı değil, **sadakatli zincirin ürettiği sitelerin kendisidir**
>    (ızgara braketi + SIRALI ilk-kök + n−½).
> 2. **E'nin fazı silinince açık kapanıyor.** `X_a`, `X_b` aynen
>    bırakılıp yalnız `E`, kendi güç tayfını bit-bit koruyan bir Gauss
>    surrogate'la değiştirildiğinde `ρ(E)` **0.6441 → 0.79743**'e
>    çıkıyor: `Rec = +0.9962`, yani **açığın %99.6'sı geri geliyor**,
>    sekiz tohumun sekizinde de.
>
> **Açığı taşıyan şey E alanının FAZ YAPISI'dır** — genliği değil,
> bant ağırlığı değil, kestirimci değil.

#### Anti-girişimin adresi (`R_η < 1`)

174-K1d'nin özdeşliği `P ≡ Kov(η_çizgi, η)` ⇒
`Var(η) = P + Kov(η_artık, η)` ⇒ **`1 − R_η ≡ Kov(η_artık,η)/Var(η)`**.
Tabloda bu özdeşlik sekiz tohumda da basamak basamak tutuyor.

> `⟨R_η⟩ = 0.81028 ± 0.00363` (8/8 tohumda < 1) ⇒
> **`Kov(η_artık, η)/Var(η) = +0.190`**: karışık gazda çizgi çıkarımı
> alanı **eksik açıklıyor**, artık η ile **pozitif** ilişkili kalıyor.
> Gerçek/sadakatli gazda işaret terstir (`R_η(Hkeskin) = 1.265` ⇒
> `Kov(artık,η) = −0.0557`).
>
> Ve aynı vekillerde **`Var(η_çizgi)/P = 1.3219`**: çizgi alanının
> varyansı, çizgi güçlerinin toplamını **%32 aşıyor**. Rastgele,
> birbirinden bağımsız fazlı çizgilerde bu oran 1 olurdu. **Çizgiler
> birbirinden bağımsız değil** — bantlar arası kovaryans/iç-bant oranı
> `+0.268` (VF1). Anti-girişim, çizgi alanının kendi içindeki bu
> yapıdan geliyor: aynı yapı hem `R_η`'yı 1'in altına, hem `ρ(E)`'yi
> Gauss'un altına indiriyor.

### 3.6 K2'nin HÜKMÜ, tek cümlede

> **Gauss-altı, referansın hatası değil (H-G1 ÖLDÜ), kestirimcinin
> artefaktı da görünmüyor (H-G3 HÜKÜMSÜZ; bağımsız kestirimci aynı
> %15'lik açığı görüyor); açığı taşıyan şey İNŞANIN KENDİSİDİR
> (H-G2 YAŞIYOR): sadakatli zincir, fazlar U(0,2π)'den bağımsız
> çekilse bile, ölçüm sitelerinde bağımsız OLMAYAN bir alan üretiyor —
> ve E bacağının bu artık faz yapısı silindiğinde açığın %99.6'sı
> kapanıyor.**

---

## 4. K1(c)/K3 — DEFTER, n_VF = 8 (`182h_defter8.py`; O1 AYNEN, k = 2.0)

**ÜREME KAPISI:** aynı kod `VF1..VF4` ile koşulduğunda
`181/K3_DEFTER_n4.json`'un **270 alanını 0 farkla** üretti.

### T3 — AYRIŞIM, DOĞRUSAL konvansiyon (n_VS = 4, n_VF = 8)

| N | Δ | ⟨N⟩_VS | ⟨N⟩_VF | ZARF ± se | KİLİT ± se | p_zarf ± se | \|Z\|/se | hüküm |
|---|---|---|---|---|---|---|---|---|
| **M** | +0.015316 | 0.448782 | 0.429529 | **+0.019253 ± 0.002310** | **−0.003937 ± 0.002310** | +1.2570 ± 0.1508 | 8.33 | **OKUNUR** |
| **Q_E** | −0.039828 | 0.288980 | 0.297532 | −0.008552 ± 0.000944 | −0.031276 ± 0.000944 | +0.2147 ± 0.0237 | 9.06 | **OKUNUR** |
| **ρ_E** | −0.019726 | 1.151223 | 1.131230 | +0.019993 ± 0.003113 | −0.039719 ± 0.003113 | −1.0135 ± 0.1578 | 6.42 | **OKUNUR** |
| **Q_X** | +0.005804 | 0.183622 | 0.186968 | −0.003346 ± 0.009888 | +0.009150 ± 0.009888 | −0.577 ± 1.704 | 0.34 | HÜKÜMSÜZ |
| **ρ_X** | +0.022410 | 1.041667 | 1.027760 | +0.013907 ± 0.017447 | +0.008503 ± 0.017447 | +0.621 ± 0.779 | 0.80 | HÜKÜMSÜZ |
| *(EK)* g_E | +0.014505 | 0.748980 | 0.736985 | +0.011995 ± 0.000366 | +0.002510 ± 0.000366 | +0.8270 ± 0.0252 | 32.75 | **OKUNUR** |
| *(EK)* g_X | +0.000520 | 0.823961 | 0.818306 | +0.005655 ± 0.006640 | −0.005135 ± 0.006640 | +10.87 ± 12.76 | 0.85 | HÜKÜMSÜZ |
| *(EK)* θ | +0.028604 | 0.882929 | 0.870626 | +0.012303 ± 0.013850 | +0.016301 ± 0.013850 | +0.430 ± 0.484 | 0.89 | HÜKÜMSÜZ |
| *(EK)* g_cal | +0.007632 | 0.508562 | 0.493567 | +0.014995 ± 0.008164 | −0.007364 ± 0.008164 | +1.965 ± 1.070 | 1.84 | HÜKÜMSÜZ |

### T4 — AYRIŞIM, LOG konvansiyon (n_VS = 4, n_VF = 8)

| N | ΔlogN | ZARFlog ± se | KİLİTlog ± se | p_zarf ± se | \|Z\|/se | hüküm |
|---|---|---|---|---|---|---|
| **M** | +0.057542 | **+0.043903 ± 0.005298** | **+0.013639 ± 0.005298** | +0.7630 ± 0.0921 | 8.29 | **OKUNUR** |
| **Q_E** | −0.044475 | −0.029161 ± 0.003238 | −0.015313 ± 0.003238 | +0.6557 ± 0.0728 | 9.01 | **OKUNUR** |
| **ρ_E** | −0.009010 | +0.017519 ± 0.002719 | −0.026528 ± 0.002719 | −1.9444 ± 0.3018 | 6.44 | **OKUNUR** |
| **Q_X** | +0.012481 | −0.018457 ± 0.055050 | +0.030938 ± 0.055050 | −1.479 ± 4.411 | 0.34 | HÜKÜMSÜZ |
| **ρ_X** | +0.014242 | +0.013429 ± 0.016938 | +0.000813 ± 0.016938 | +0.943 ± 1.189 | 0.79 | HÜKÜMSÜZ |
| *(EK)* g_E | +0.024547 | +0.016145 ± 0.000492 | +0.008401 ± 0.000492 | +0.6577 ± 0.0200 | 32.84 | **OKUNUR** |
| *(EK)* g_X | +0.000739 | +0.006882 ± 0.008045 | −0.006144 ± 0.008045 | +9.32 ± 10.89 | 0.86 | HÜKÜMSÜZ |
| *(EK)* θ | +0.031518 | +0.013993 ± 0.015862 | +0.017524 ± 0.015862 | +0.444 ± 0.503 | 0.88 | HÜKÜMSÜZ |
| *(EK)* g_cal | +0.026024 | +0.029910 ± 0.016064 | −0.003886 ± 0.016064 | +1.149 ± 0.617 | 1.86 | HÜKÜMSÜZ |

### T6 — Boğulan satırlar: **HİÇBİRİ AÇILMADI** (P2 öngörüsü 18/18 tuttu)

| satır | konv. | \|Z\|/se (4) | öngörü (8) | **ölçülen (8)** | sd_VF(4) → sd_VF(8) | hüküm(4) → hüküm(8) |
|---|---|---|---|---|---|---|
| M | doğ. | 10.34 | 12.45 | **8.33** | 0.00336 → **0.00537** | OKUNUR → OKUNUR |
| M | log | 10.26 | 12.45 | **8.29** | 0.00787 → **0.01247** | OKUNUR → OKUNUR |
| Q_E | doğ. | 9.95 | 10.93 | **9.06** | 0.00108 → **0.00164** | OKUNUR → OKUNUR |
| Q_E | log | 9.90 | 10.83 | **9.01** | 0.00361 → **0.00551** | OKUNUR → OKUNUR |
| ρ_E | doğ./log | 5.76 / 5.78 | 6.32 / 6.35 | 6.42 / 6.44 | 0.00377 → 0.00461 | OKUNUR → OKUNUR |
| g_E | doğ./log | 28.76 / 28.75 | 33.23 / 33.31 | 32.75 / 32.84 | 0.00060 → 0.00059 | OKUNUR → OKUNUR |
| **Q_X** | doğ./log | 0.42 / 0.41 | 0.49 / 0.47 | **0.34 / 0.34** | 0.01696 → 0.01500 | HÜKÜMSÜZ → HÜKÜMSÜZ |
| **ρ_X** | doğ./log | 0.51 / 0.51 | 0.59 / 0.59 | **0.80 / 0.79** | 0.03009 → 0.02650 | HÜKÜMSÜZ → HÜKÜMSÜZ |
| **θ** | doğ./log | 0.88 / 0.87 | 1.03 / 1.03 | **0.89 / 0.88** | 0.02573 → 0.02009 | HÜKÜMSÜZ → HÜKÜMSÜZ |
| **g_X** | doğ./log | 0.84 / 0.84 | 0.97 / 0.97 | **0.85 / 0.86** | 0.01133 → 0.00995 | HÜKÜMSÜZ → HÜKÜMSÜZ |
| **g_cal** | doğ./log | 1.69 / 1.70 | 1.94 / 1.97 | **1.84 / 1.86** | 0.01359 → 0.01201 | HÜKÜMSÜZ → HÜKÜMSÜZ |

```
n=4 → n=8 AÇILAN satır : HİÇBİRİ
n=4 → n=8 KAPANAN satır: HİÇBİRİ
P2 öngörüsünden sapma  : 0 / 18
```

> **Ön-kaydın HÜKÜM öngörüsü 18/18 tuttu; MEKANİZMA öngörüsü tutmadı**
> ve tutmadığı burada yazılıdır. P2, `sd_VF` sabit kalırsa `|Z|/se`nin
> M'de 10.34 → 12.45'e çıkacağını söylemişti; **çıkmadı, 8.33'e
> düştü**, çünkü dört yeni tohum `sd_VF`'i M'de ×1.60, Q_E'de ×1.52
> **büyüttü**. Yani **n = 4'ün VF hata çubukları — tıpkı 181'in VS
> tarafında bulduğu gibi — gerçek saçılımı EKSİK kestiriyordu.**
> Bu 180/181'in sayılarına düzeltme değil, **çözünürlüğün kendisinin
> ölçümüdür**: 180-T3'ün `M` ZARF payı `+0.020481 ± 0.002867` idi,
> n = 8'de `+0.019253 ± 0.002310` — merkez 0.43σ kaydı, çubuk
> ×0.81 daraldı; hüküm değişmedi.
>
> `ρ_X`'in `|Z|/se`'si 0.51 → 0.80'e çıktı (merkez kaydı), `Q_X`'inki
> 0.42 → 0.34'e düştü — ikisi de eşiğin (2.0) çok altında.
> **`Q_X`, `ρ_X`, `θ`, `g_X` için ön-kayıt bunu YAPISAL gerekçeyle
> söylemişti:** `se_ZARF ≥ se_VS` her zaman ve bu dört satırda
> `|ZARF| < 2·se_VS`; **VF tarafı sıfıra gitse bile açılmazlar.**
> 181, onları "VF ailesinin saçılımı boğuyor" diye bırakmıştı; 182
> bunu düzeltiyor: **onları VS tarafı da boğuyor, ve VF tohumu
> eklemek yetmez.**
>
> `g_cal`: 181'in ön-kaydı n=4'te "OKUNUR olacak" demiş, tutmamıştı;
> 182'nin ön-kaydı n=8'de "HÜKÜMSÜZ" dedi ve **tuttu** (1.84 / 1.86;
> eşik 2.0). İki görevde iki farklı yönde sınandı, ikisinde de
> HÜKÜMSÜZ.

### T7 — K4: karıştırma çöküşleriyle çapraz tutarlılık (n_VF = 8)

| N | K_N (Hk zarfı) n=8 | *(n=4)* | K′_N (son zarfı) | logK − logK′ | KİLİTlog | kalıntı |
|---|---|---|---|---|---|---|
| **M** | **+%66.09** | +%65.01 | +%63.84 | +0.013639 | +0.013639 | 0.0e+00 |
| **Q_E** | **−%67.50** | −%67.44 | −%67.00 | −0.015313 | −0.015313 | 4.2e−17 |
| **ρ_E** | −%48.56 | −%48.51 | −%47.18 | −0.026528 | −0.026528 | 2.8e−17 |
| Q_X | −%59.66 | −%59.31 | −%60.89 | +0.030938 | +0.030938 | 0.0e+00 |
| ρ_X | −%34.24 | −%34.03 | −%34.29 | +0.000813 | +0.000813 | −6.9e−18 |
| *(EK)* g_E | +%26.26 | +%26.22 | +%25.20 | +0.008401 | +0.008401 | 0.0e+00 |
| *(EK)* θ | −%2.56 | −%2.92 | −%4.26 | +0.017524 | +0.017524 | 0.0e+00 |

Özdeşlik makine sıfırında duruyor. 176'nın manşet çöküşleri n = 8'de de
birebir: `Q_E` −%67.5, `ρ_E` −%48.6, `g_E` +%26.3, `M` +%66.1.

---

## 5. HÜKÜMLER (özet)

### (i) GAUSS-ALTI: **HÜKÜM VERİLDİ, n = 8** (`182e`)

`⟨ρ₃⟩_VF = 0.453014 ± 0.002314` (se), Gauss'un **%10.82 altı**,
`z_G = −23.7`, **8/8 tohum tek tek altında**. 180'in n=4'te sahipsiz
bıraktığı olgu artık **adlı ve sağlam**: *kilitsiz gaz Gauss'a
çökmüyor, Gauss'un altına geçiyor.*

### (ii) KIRPMA HAKEMİ: **HÜKÜMSÜZ — aynen, ve YAPISAL olarak öyle** (`182e`)

Saçılım 0.016091 (eşik 0.010225); yeni dört tohum range'i **hiç
değiştirmedi**. Ön-kayıtlı P1 öngörüsü tuttu. **H-180a hakemi VF
ailesini büyüterek kazanılamaz** — range istatistiği monotondur.
Eşik gevşetilmedi. `DAL B` 180'deki gibi çürük (`C = +3.69`).

### (iii) E-BACAĞININ ANATOMİSİ: açık **BANT**'ta, adres **İNŞA** (`182g`/`182i`)

* `PAY_ÇİZGİ = −0.071` ⇒ açığın %107'si bant, %−7'si çizgi.
* `κ_E = 0.8101` Gauss'un **üstünde**; `kurt(E) = −0.263` (platikurtik).
* `s_E = +0.382` vs `s_Xa/s_Xb ≈ +0.07/+0.04` ⇒ τ-bağımlı, spektral.
* **H-G1 ÖLDÜ** (ρ°(E) = 0.79474, referans %0.4 içinde geri geliyor).
* **H-G3 HÜKÜMSÜZ** (7/8 tohum geçti, VF7 0.05346 > 0.05; gevşetme yok).
* **H-G2 YAŞIYOR**: `korr(Xa,Xb) = +0.1006 ± 0.0056` (50σ; surrogate'ta
  −0.0004), `Rec = +0.9962` (8/8).
* **Anti-girişimin adresi:** `1 − R_η ≡ Kov(η_artık,η)/Var(η) = +0.190`
  (8/8 tohumda pozitif; Hkeskin'de negatif), ve `Var(η_çizgi)/P = 1.322`
  — çizgiler birbirinden bağımsız değil.

### (iv) BOĞULAN SATIRLAR: **HİÇBİRİ AÇILMADI** (`182h`)

`Q_X`, `ρ_X`, `θ`, `g_X`, `g_cal` — beşi de, iki konvansiyonda da
HÜKÜMSÜZ. Ön-kayıt bunu 18/18 doğru öngördü ve dördü için (`Q_X`,
`ρ_X`, `θ`, `g_X`) **yapısal ulaşılamazlığı** VF tarafında da gösterdi.
**180 defterine işlenecek yeni satır yoktur; çapasız paylar
GÜNCELLENMEDİ** (yalnız `M`, `Q_E`, `ρ_E`, `g_E`'nin merkez/çubukları
n=8 ile yenilendi — hükümleri değişmedi).

### (v) n = 4'ÜN VF HATA ÇUBUKLARI GÜVENİLMEZDİ

`sd_VF` M'de ×1.60, Q_E'de ×1.52 **büyüdü**; `|Z|/se` M'de 10.34 →
8.33'e **düştü** (öngörü 12.45 idi). Bu, 181'in VS tarafında bulduğu
olgunun VF tarafındaki eşidir ve aynı sonucu verir: **iki-dört
tohumluk saçılım kestirimleri iyimserdir.**

### (vi) KAPSAM DIŞI, DEĞİŞMEDİ

`W_pos` borcu (179 ŞART ③) bu görevde de ele alınmadı. VS ailesi n = 4'te
bırakıldı. Gerçek gazın (`son`) hiçbir niceliği yeniden ölçülmedi.
176a'nın `F0`/`R_bant` ıskası devralındı ve kurtarılmadı.

---

## 6. DOSYALAR

| betik | ne yapar | çıktı |
|---|---|---|
| `182_configs/182a_onkayit.py` | K0 ön-kayıt (mühür, eşikler, ölüm koşulları, P1/P2) | `182/ONKAYIT_182.json` |
| `176_configs/176b_vekil_insa.py` *(aynen)* | VF5..VF8 inşası | `182/log_182b_VF*.txt`, `176/insa_VF*.json` |
| `182_configs/182c_olcum.py` | `176c_olcum.kos` + `180a.defter` | `182/G_VF5..8.json` |
| `180_configs/180d_kirpma_kos.py` *(aynen)* | `169_k2b.main(g,0.40,0.95)` | `169/K2b_VF5..8.json` |
| `176_configs/176d_K1.py` *(aynen)* | `174b` girişim makinesi (R_η) | `176/K1_VF3..8.json` |
| `182_configs/182e_hakem8.py` | hakem n=8 + Gauss-altı hükmü (R-KAPISI) | `182/K2_HAKEM_n8.json` |
| `182_configs/182g_anatomi.py` | κ/β, τ-eğimi, ρ_J, surrogate'lar | `182/ANATOMI_VF1..8.json` |
| `182_configs/182i_hipotez.py` | H-G1/G2/G3 yarışı (ön-kayıtlı ölümler) | `182/K2_HIPOTEZ.json` |
| `182_configs/182h_defter8.py` | defter n_VF=8, O1 aynen (R-KAPISI) | `182/K3_DEFTER_n8.json` |
| `182_configs/182j_figur.py` | 6 panel | `182_gauss_alti.png` |
