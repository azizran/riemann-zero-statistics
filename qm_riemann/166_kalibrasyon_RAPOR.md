# 166 — 165'İN KALİBRASYONUNUN TÜRETİMİ: ÜÇÜNCÜ-MOMENT ÖNGÖRÜSÜNÜ PARAMETRESİZ YAPMAK

**Hedef.** KALEM (2 Eylül akşam) 165'in `n_u2 = 0.15–0.31` kalibrasyonunun
kimliğini sordu. Üç aday: **H-K2** `kalib = ρ(τ_bant)` (tarak-sağkalımı),
**H-K1** bacak-sönümleri çarpımı, **H-K3** karışımlar. Katil test: **aynı
yasa üç gazın (`son`, `Hkeskin`, `HA4`) kalibrasyonunu, her gazın KENDİ
ölçülü girdileriyle, gaz-başına serbest ölçek OLMADAN vurmalı.**

> **Bu koşunun bir cümlelik özeti: H-K2 (ρ) DÜŞTÜ — bant şeklinde
> %44–73 artık bırakıyor ve `HA4` gazında ÖLDÜ (orada ρ, τ > 0.68'de
> sıfıra iner, kalibrasyon ise 0.15–0.19'da durur). H-K1 (BACAK
> SÖNÜMÜ) KAZANDI: kalibrasyonun BANT ŞEKLİ, X-bacağının ölçülü
> Debye–Waller çarpanları `W_amp(τ)·W_X(τ)` ile — sıfır uydurma
> parametreyle — %11 rms içinde çıkıyor, ve bu yasa 165'in
> u2-sabitlemesini tek bir EVRENSEL sayıyla (c = 0.5647) değiştirince
> ±%25 sayımı 12/20 → **13/20** (sağlıklı pencerede 13/15) oluyor.
> AMA c'nin kendisi TÜRETİLEMEDİ ve üç gazda −%12 … +%8 oynuyor: öngörü
> "parametresiz" DEĞİL, **"tek evrensel sayı + parametresiz bant
> şekli"** mertebesine geldi.**

---

## Kısa hüküm

1. **T1 ZİNCİRİ BİT DÜZEYİNDE YENİDEN ÜRETİLDİ.** `166_T1.py` 165'in
   ölçüm+öngörü zincirini (163'ün `bant_adaylari`/`olc_cizgi`'si,
   165'in `Model165`'i) import ederek koşar. Üç gazın **27 ölçülen
   bandında**, 13 büyüklükte (`A²s_k`, `A²u_k`, öngörüleri, `norm_u2`,
   `τ_eff`, `Ç1`, `Ç2`, `Ç3`) 165'in kayıtlı JSON'una karşı maks fark
   **0.0e+00** (istenen 1e−9'un çok ötesi). Yeni olan: **çizgi bazında
   kayıt** ⇒ `norm_u2` ve ölçüm/öngörü oranı için de bant-içi 8-grup
   jackknife.
2. **KALİB'İN İKİ TANIMI AYRIŞTIRILDI.**
   `KALİB_u2 = A²u2_ölç/A²u2_öng` (165'in KULLANDIĞI, yalnız ikinci
   momentler, jackknife %0.35–6.5, τ_eff ≤ 0.70'te ≤ %1.4) ve `KALİB_s2 = A²s2_ölç/A²s2_öng`
   (üçüncü momentin İSTEDİĞİ). İkisi τ_eff ≈ 0.62'de kesişiyor;
   `KALİB_s2` daha dik iniyor (τ ≥ 0.74'te 165'in bilinen çöküşü).
3. **ρ ÜÇ GAZDA ÖLÇÜLDÜ (144'ün yöntemiyle) ve 144 DOĞRULANDI.**
   `son`/0.52-taban/144-bantları: **0.3806 / 0.2664 / 0.1589 / 0.0796**
   ↔ 144'ün 0.381/0.266/0.159/0.080. Aynı estimatör iki tabanda ve iki
   bant kümesinde koşuldu (146'nın taban-bağımlılık dersi açıkta).
4. **H-K2 DÜŞTÜ — İKİ BAĞIMSIZ SEBEPLE.** (i) *Bant şekli*: 0.52-tabanda
   ρ, hüküm penceresinde 3.9 kat inerken kalibrasyon yalnız 1.55 kat
   iniyor; sıfır-ölçek rms **%56**, tek-ölçek **%56**, 0.40-tabanda
   **%73 / %44**. (ii) *ÇAPRAZ-GAZ*: `HA4`'ün erfc penceresi merdiveni
   τ > 0.68'de sıfırlıyor ⇒ ρ(0.72–0.76] = **+0.026**, ρ(0.76–0.80] =
   **−0.002**; oysa aynı bantlarda `KALİB_u2` = **0.189 / 0.146**.
   **Kalibrasyon, ρ'nun sıfırlandığı yerde sıfırlanmıyor.**
5. **H-K1'İN BACAK YAPISI TÜRETİLDİ ve ÖLÇÜLDÜ.** 165'in tek-site
   çekirdeğinde üç alan bacağı üç FARKLI faz referansı taşır:
   **E-bacağı** `A = πτ_q`, dalgalanma `ds` ⇒ `W_amp(τ)=⟨cos πτ·ds⟩`;
   **X-bacağı** ayrıca `A = 2πτ_q`, dalgalanma `X̃` ⇒
   `W_X(τ)=⟨e^{−2πiτX̃}⟩`; **taşıyıcı** ölçümde de öngörüde de AYNI
   sitede olduğu için doğrudan bir sönüm taşımaz. Kalem sonucu
   `y_q = ½ b_q W_amp(τ_q)[W_X(τ_q)e^{−2πiτ_q} + 1]` çıplak limitte
   165 §2a'yı (`B_q e^{−iπτ_q}`) verir ve ölçümle τ ≤ 0.30'da **%1–5**
   içinde uyuşur (§3b).
6. **YARIŞIN GALİBİ: `kalib(τ_bant) = c · W_amp(τ)·W_X(τ)`.**
   Hüküm bantlarında (20) log-artık **rms %11.4, maks %23.0** (bütün
   yasaların EN KÜÇÜK maksimumu), çapraz-gaz ölçek yayılımı **1.225**.
   Kardeş üyeler `c·W_amp` (%9.8), `c·W_X` (%10.6), `c·W_pos` (%10.8)
   sayısal olarak **ayırt edilemez**; hepsi aynı aileden. Karşılaştırma:
   bant-bağımsız sabit **%17.0**, ρ aileleri **%44–73**.
7. **T3 — KAPANIŞ: 13/20 (165: 12/20).** Ölçek `c` YALNIZ `KALİB_u2`ye
   (ikinci momentlere) oturtuldu, üçüncü momentten hiçbir girdi yok.
   `c = 0.5647` (üç gaz ortak). Sayım: `son` 5/7, `Hkeskin` 5/7,
   `HA4` 3/6 ⇒ **13/20**; τ_eff ∈ 0.54–0.70 penceresinde **13/15**
   (165: 12/15). Kaçan yedi bandın **beşi** τ_eff ≥ 0.73'tedir
   (165 §3c'nin bilinen çöküşü); kalan ikisi `HA4`'ün 0.657 (oran
   **1.252** — sınırı binde ikiyle kaçırıyor) ve 0.696 bantlarıdır.
8. **AMA SIFIR ÖLÇEK OLMADI.** Denenen bütün SIFIR-ölçek yasaları
   (ρ, `g·W_pos`, `1/λ_E`, `1/λ_X`, `1/(λ_Eλ_X²)`, `1/(λ_XΛ_EΛ_X)`)
   %27–79 artık bırakıyor. `c`'nin gazdan gaza kalan yayılımı
   **`son` 0.5822 / `Hkeskin` 0.6100 / `HA4` 0.4980** (maks/min = 1.225;
   ortak `c = 0.5647`'ye göre +%3 / +%8 / −%12).
   `W_X` üyesinde `son` ve `Hkeskin` **%0.06 içinde** aynı
   (0.4053 / 0.4051), `HA4` **%17.5 aşağıda**. Bu, açık borç olarak yazılıyor.
9. **T4 — ŞİŞME 1/ρ TİPİ DEĞİL, GRAM TİPİ; ve KAPALI FORMU BULUNDU.**
   Özdeşlik türetildi ve sayısal olarak doğrulandı:
   `Var(X_mod)/Var(X̃) = Λ_X · Γ_X` **üç gazda da ÜÇ HANEDE TAM**
   (1.585 / 1.562 / 1.520), `Var(E_mod)/Var(η) = Λ_E·Γ_E` %2 içinde
   (2.139 ↔ 2.180). Burada `Λ = Σ|h|²Reλ/Σ|h|²` merdivenin
   **Gram sızıntısı** (ve `Λ_E ≡ 1/g_E` özdeş çıkıyor), `Γ = ½Σ|h|²/Var`
   ise **Gram atfetmesi**. 165 §7'nin 2.1–3.7'si böylece iki ölçülen
   çarpana ayrıldı; **ρ ile ilgisi yok.**
10. **Ç4'ÜN KAPALI FORMU HÂLÂ YOK.** Bu koşu kalibrasyonun ŞEKLİNİ
    türetti, mekanizmasını (Ç4'ün ortak kazancı) türetmedi: naif
    bacak-DW muhasebesi Ç1(B) ailesinde bir **kazanç** (`e^{+A₁²σ²}`)
    öngörüyor, ölçülen ise bir **sönüm**. Bu çelişki §3e'de açık
    yazılmıştır ve 167'nin işidir.

---

## 1. T1 — KALİB'İN HAM TABLOSU (u2-sabitlemesinden ÖNCE)

### 1a. Ne ölçüldü

`166_configs/166_T1.py`, 165'in T1 sürücüsünün aynısını koşar
(`kaynak = olculen`, `τ_c = 0.95`, `taban = 0.40`, `IZGARA_T1`), ama
her bandın her çizgisini saklar. Bant birleştirmesi 160'ınkiyle
birebirdir:
`A²v = Σ(p_on A_on² v_on − p_off A_off² v_off)/Σ(p_on−p_off)`.
Hata: **bant içi 8-grup jackknife**, ORAN tahmincileri için silme-1
yeniden hesaplanarak (doğrusal hata yayılımı yapılmadı).

| | tanım | ne söyler |
|---|---|---|
| `KALİB_s2` | `A²s2_ölç / A²s2_öng` | üçüncü momentin İSTEDİĞİ kalibrasyon |
| `KALİB_u2` | `A²u2_ölç / A²u2_öng` (= 165'in `norm_u2`) | 165'in KULLANDIĞI (yalnız ikinci momentler) |

### 1b. Tablo (hüküm bantları: `lo ≥ 0.52`, 160'ın sağlık kuralı)

| gaz | τ_eff | **KALİB_s2 ± jk** | **KALİB_u2 ± jk** | ρ(0.52-tb) | ρ(0.40-tb) | W_amp | W_X | W_pos |
|---|---|---|---|---|---|---|---|---|
| **son** | 0.5392 | +0.3909 ± 0.0094 | 0.3052 ± 0.0018 | +0.4177 | +0.6347 | 0.7869 | 0.7303 | 0.6454 |
| son | 0.5791 | +0.3198 ± 0.0055 | 0.2855 ± 0.0018 | +0.3484 | +0.5757 | 0.7587 | 0.6961 | 0.6026 |
| son | 0.6188 | +0.2785 ± 0.0036 | 0.2739 ± 0.0011 | +0.3024 | +0.5020 | 0.7298 | 0.6613 | 0.5597 |
| son | 0.6578 | +0.2313 ± 0.0072 | 0.2577 ± 0.0022 | +0.2442 | +0.4243 | 0.6998 | 0.6257 | 0.5165 |
| son | 0.6980 | +0.1989 ± 0.0070 | 0.2403 ± 0.0025 | +0.1957 | +0.3361 | 0.6696 | 0.5903 | 0.4742 |
| son | 0.7382 | +0.1318 ± 0.0055 | 0.2206 ± 0.0064 | +0.1543 | +0.2773 | 0.6382 | 0.5542 | 0.4320 |
| son | 0.7766 | +0.0423 ± 0.0208 | 0.1972 ± 0.0063 | +0.1083 | +0.1837 | 0.6082 | 0.5200 | 0.3928 |
| **Hkeskin** | 0.5393 | +0.3998 ± 0.0116 | 0.2863 ± 0.0020 | +0.5091 | +0.7640 | 0.7641 | 0.7128 | 0.6352 |
| Hkeskin | 0.5792 | +0.3196 ± 0.0072 | 0.2706 ± 0.0020 | +0.4472 | +0.7201 | 0.7331 | 0.6766 | 0.5915 |
| Hkeskin | 0.6189 | +0.2646 ± 0.0032 | 0.2586 ± 0.0009 | +0.4056 | +0.6546 | 0.7013 | 0.6399 | 0.5478 |
| Hkeskin | 0.6581 | +0.2153 ± 0.0066 | 0.2439 ± 0.0019 | +0.3516 | +0.5865 | 0.6683 | 0.6025 | 0.5039 |
| Hkeskin | 0.6983 | +0.1776 ± 0.0071 | 0.2302 ± 0.0018 | +0.3018 | +0.4955 | 0.6352 | 0.5653 | 0.4611 |
| Hkeskin | 0.7386 | +0.1154 ± 0.0080 | 0.2185 ± 0.0048 | +0.2701 | +0.4496 | 0.6010 | 0.5275 | 0.4183 |
| Hkeskin | 0.7774 | +0.0073 ± 0.0188 | 0.1987 ± 0.0045 | +0.2205 | +0.3477 | 0.5681 | 0.4919 | 0.3788 |
| **HA4** | 0.5392 | +0.3312 ± 0.0101 | 0.2655 ± 0.0026 | +0.5364 | +0.8377 | 0.7582 | 0.7205 | 0.6403 |
| HA4 | 0.5788 | +0.2763 ± 0.0044 | 0.2431 ± 0.0022 | +0.4195 | +0.7188 | 0.7259 | 0.6856 | 0.5970 |
| HA4 | 0.6181 | +0.2336 ± 0.0047 | 0.2311 ± 0.0012 | +0.3116 | +0.5468 | 0.6927 | 0.6503 | 0.5538 |
| HA4 | 0.6566 | +0.1824 ± 0.0045 | 0.2104 ± 0.0020 | **+0.1790** | +0.3522 | 0.6581 | 0.6144 | 0.5103 |
| HA4 | 0.6962 | +0.1294 ± 0.0035 | 0.1891 ± 0.0032 | **+0.0857** | +0.1795 | 0.6231 | 0.5789 | 0.4679 |
| HA4 | 0.7349 | +0.0482 ± 0.0129 | 0.1463 ± 0.0095 | **+0.0256** | +0.0750 | 0.5868 | 0.5428 | 0.4254 |

(Hüküm dışı iki alt bant da ölçüldü ve şekil tanısında kullanıldı:
`son` τ_eff = 0.4607/0.4997'de `KALİB_u2` = 0.3290/0.3175. `HA4`'ün
τ_eff = 0.9217 bandı 165'te olduğu gibi ‡ ve hiçbir hükme girmedi.)

> **Okuma.** Hüküm penceresinde `KALİB_u2` yalnız **1.55 kat** iniyor
> (0.305 → 0.197). `ρ` aynı pencerede **3.86 kat** (0.52-taban) ya da
> **3.45 kat** (0.40-taban) iniyor. Bu tek başına H-K2'nin bant-şekli
> sınavını kaybettiğini söylüyor; §2'de sayıya dökülüyor.

---

## 2. ρ(τ_bant)'ın ÖLÇÜMÜ — 144'ün YÖNTEMİ, ÜÇ GAZDA

`166_configs/166_rho.py` 144'ün estimatörünü **birebir** koşar:

```
ρ(bant) = Σ_çizgi [ |ĥ(w)|² − |ĥ(w+gap/2)|² ] / Σ_çizgi b_q²
ĥ(W) = 2⟨η e^{−iW m_n}⟩ ,  b_q = 2a_q sin(πτ_q)
```

İki η konvansiyonu (146'nın dersi gizlenmiyor) ve iki bant kümesi:

| | taban / cap | ne |
|---|---|---|
| `t0.52_c720` | 144/145/146'nın standardı | ρ'nun tarihsel değerleri burada |
| `t0.40_c4000` | **165'in KENDİ konvansiyonu** | kalibrasyon burada doğdu |
| `b144` | 144'ün bantları, tohum 11, 260 örnek | **doğrulama** |
| `b165` | 163/165'in `bant_adaylari`'sı | kalib ile AYNI çizgiler — yarışın girdisi |

### 2a. V-R1: 144 doğrulandı

| bant | (0.525,0.60] | (0.60,0.70] | (0.70,0.78] | (0.78,0.84] |
|---|---|---|---|---|
| **166, `son`, 0.52-taban** | **+0.3806 ± 0.0015** | **+0.2664 ± 0.0031** | **+0.1589 ± 0.0030** | **+0.0796 ± 0.0038** |
| 144'ün kaydı | 0.381 | 0.266 | 0.159 | 0.080 |
| `Hkeskin`, 0.52-tb | +0.4753 | +0.3704 | +0.2712 | +0.1885 |
| `HA4`, 0.52-tb | +0.4743 | +0.2331 | **+0.0353** | **−0.0101** |
| `son`, 0.40-tb | +0.6037 | +0.4522 | +0.2799 | +0.1508 |
| `Hkeskin`, 0.40-tb | +0.7401 | +0.6077 | +0.4476 | +0.3107 |
| `HA4`, 0.40-tb | +0.7759 | +0.4295 | **+0.0888** | **−0.0167** |

146'nın E2'si (ρ taban-bağımlıdır, 0.40-tabanda ~1.5×) burada da aynen
görülüyor ve **her satırda hangi konvansiyon olduğu yazılıdır.**

### 2b. Çözünürlük / örtüşme kısıtları (açık)

* Aday çizgiler `gap ≥ 2.5·dres` filtresinden geçiyor (`dres = 4.01e−05`);
  ara-nokta referansı `W' = w + gap/2`. τ > 0.80'de merdiven çizgi
  aralığı `dres`'e yaklaşıyor, bu yüzden hüküm `lo ≤ 0.80` ile sınırlı.
* ρ **Gram-atfetme gölgesidir** (146'nın dersi): konvansiyonsuz bir
  "mutlak soğurma eğrisi" yoktur. Bu koşu ρ'yu iki konvansiyonda birden
  koşarak yarışı konvansiyon seçimine bağımlı olmaktan çıkardı —
  **H-K2 her iki konvansiyonda da kaybediyor.**
* `HA4`'te τ > 0.68'de `ρ ≈ 0` bir ölçüm gürültüsü değil, **inşa
  gerçeğidir**: gazın erfc penceresi (0.68/0.125) merdiveni orada
  söndürüyor, yani o bantlarda "bant çizgisi" fiziksel olarak yok.
  Kalibrasyonun orada 0.15–0.19'da durması bu yüzden H-K2 için
  **katil**dir.

---

## 3. BACAK SÖNÜMLERİNİN TÜRETİMİ (H-K1)

### 3a. Hangi bacak hangi `A` ile

Merdiven (163 §1, 165 §2a):
```
ds_n = Σ_q 2a_q sin(ω_q g_n/2)·cos(ω_q m_n) ,  g_n = ḡ(1+ds_n),
m_n = (z_n+z_{n+1})/2 ,  ω_q ḡ = 2πτ_q ,  165'in sitesi  s_n ≡ m_{n+1}
```

**(E-bacağı)** `η_{n+1}`'in q bileşeni SİTENİN KENDİSİNDE oturur
(`cos(ω_q m_{n+1}) = cos(ω_q s_n)`), ama genliği yerel boşlukla
modülelidir:
```
η_{n+1} ⊃ 2a_q sin(πτ_q(1+ds_{n+1}))·cos(ω_q s_n)
h'_q = 2⟨η_{n+1}e^{−iω_q s_n}⟩ = b_q·⟨cos(πτ_q ds)⟩ + 2a_q cos(πτ_q)⟨sin(πτ_q ds)⟩
     ≈ b_q · W_amp(τ_q) ,   **W_amp(τ) = ⟨cos(πτ·ds)⟩ ≈ e^{−½π²τ²σ_ds²}**
```
⇒ **E-bacağının `A`'sı `πτ_q`, dalgalanması `ds`.** (Bu, 165 §6b'nin
`sin⁴(πτ)` üyesinin kaynağıyla aynı Debye–Waller ailesidir.)

**(X-bacağı)** `X̃_n = ½(ds_n + ds_{n+1})`; ikinci yarı sitede,
BİRİNCİ yarı bir bond geride:
`ω_q(m_{n+1}−m_n) = ω_q ḡ(1+X̃_n) = 2πτ_q(1+X̃_n)` ⇒
```
y_q = 2⟨X̃0_n e^{−iω_q s_n}⟩ = ½ b_q W_amp(τ_q)·[ W_X(τ_q) e^{−2πiτ_q} + 1 ]
      **W_X(τ) = ⟨e^{−2πiτ·X̃}⟩ ≈ e^{−2π²τ²σ_X̃²}**
```
⇒ **X-bacağı İKİ çarpan taşır: `A = πτ_q` (ds) ve `A = 2πτ_q` (X̃).**
Çıplak limit `W_amp = W_X = 1`: `y_q = b_q cos(πτ_q)e^{−iπτ_q} =
B_q e^{−iπτ_q}` — **165 §2a AYNEN**. `W_X < 1` iken `arg y_q + πτ_q`
sıfır değildir ve τ = ½'te π/2'ye fırlar: **165 V3'ün "τ = 0.58'de
2.36 rad" kaydının kalem karşılığı budur.**

**(TAŞIYICI)** `e^{−iω_Q s_n}` hem ölçümde hem öngörüde AYNI sitede
değerlendirilir (`olc_cizgi`'de `mid[1:]`, `ongor`'da `Mo.s`), ve
normalizasyon `h_Q`/`⟨ρ⟩` ölçülendir. **Taşıyıcı doğrudan bir sönüm
çarpanı taşımaz;** τ_Q bağımlılığı ancak alanların ω_Q'daki
gösterim hatasından gelebilir.

**(DÖRDÜNCÜ, ÖLÇÜLEN ÇARPAN — GRAM SIZINTISI)** Merdiven frekansları
sonlu tarakta dik değildir:
```
λ_E(ω) ≡ 2⟨E_mod e^{−iωs}⟩ / h'_ω ,   λ_X(ω) ≡ 2⟨X_mod e^{−iωs}⟩ / y_ω  (≥ 1)
```
163'ün V4'ü tam olarak `u0_pred ≈ Re λ_E(ω_Q)`'yu 1'e çeker; 165'in
`u2`-normalizasyonu aynı ailenin ikinci-momentli üyesidir.

### 3b. Türetimin ÖLÇÜMLE sınanması (`166_bacak.py`)

Karakteristik fonksiyonlar TAM alındı (Gauss varsayımı yok);
`σ_ds / σ_X̃ / σ_Ĉ` = `son` 0.4092/0.2338/0.2730, `Hkeskin`
0.4313/0.2420/0.2777, `HA4` 0.4343/0.2396/0.2755.

**`son` (gerçek gaz) — E-bacağı sınavı `|h'_q|/b_q =? W_amp(τ)`**

| τ | 0.421 | 0.461 | 0.500 | 0.540 | 0.579 | 0.619 | 0.659 | 0.699 |
|---|---|---|---|---|---|---|---|---|
| `R_h = h'/b` | 0.8854 | 0.8574 | 0.8297 | 0.7971 | 0.7601 | 0.7122 | 0.6562 | 0.5901 |
| `W_amp` (öngörü) | 0.8642 | 0.8395 | 0.8140 | 0.7869 | 0.7587 | 0.7298 | 0.6998 | 0.6694 |
| fark | +2.5% | +2.1% | +1.9% | +1.3% | **+0.2%** | −2.4% | −6.2% | −11.8% |

> **E-bacağının Debye–Waller çarpanı gerçek gazda τ ∈ 0.42–0.60'ta
> ±%2.5 içinde ölçüldü.** τ ≥ 0.62'de ölçülen daha hızlı iniyor
> (çizgi gösteriminin bilinen sınırı, 165 §7).
> **`τ ≤ taban = 0.40` çizgilerinde `h' ≡ 0`** (η zinciri onları
> siliyor) — tabloda `R_h = 0.0000` olarak görünüyor ve 165'in
> `|hp|/b = 0.0000` kaydını doğruluyor.

**X-bacağı sınavı `|y_q|/B_q =? |½b W_amp(W_X e^{−2πiτ}+1)|/B`**
(düşük τ, sızıntının ihmal edilebildiği yer):

| τ | 0.180 | 0.224 | 0.261 | 0.301 | 0.340 | 0.380 |
|---|---|---|---|---|---|---|
| `son` ölçülen | 0.9498 | 0.9301 | 0.9005 | 0.8847 | 0.8521 | 0.8308 |
| `son` kalem | 0.9570 | 0.9350 | 0.9135 | 0.8894 | 0.8663 | 0.8487 |
| `Hkeskin` ölçülen | 0.9654 | 0.9531 | 0.9292 | 0.9216 | 0.9017 | 0.8911 |
| `Hkeskin` kalem | 0.9529 | 0.9288 | 0.9052 | 0.8789 | 0.8534 | 0.8339 |

Fark **%1–5**; ve **her iki taraf da τ = ½'te aynı kutbu** üretiyor
(`B_q = a sin 2πτ → 0`; ölçülen 28.3, kalem 15.9 — işaret ve yer
aynı, büyüklük kutup yakınında anlamsız). **Kalem formu doğrulandı.**

### 3c. Sentetik gazların E-bacağı DAHA AZ sönümlü (yeni sistematik)

τ ≈ 0.54'te `R_h`: `son` **0.797**, `Hkeskin` **0.875**, `HA4`
**0.972**; `W_amp` üçünde de 0.758–0.787. Yani 164'ün inşası çizgileri
yarım-gap modülasyon formülünün öngördüğünden daha **koherent**
gerçekleştiriyor. Bu, `HA4`'ün kalibrasyonunun neden sistematik olarak
küçük olduğunu **niteliksel olarak** açıklar (daha güçlü ölçülen
`h'` ⇒ daha büyük öngörü ⇒ daha küçük kalib) ama **niceliksel olarak
kapatmıyor** (§4c). `HA4`'te `R_h > 1` (τ ≥ 0.66) bir inşa
artefaktıdır: `b` erfc penceresini içerir ve orada sıfıra gider.

### 3d. Ölçülen Gram sızıntıları

| gaz | Λ_E (merdiven ort.) | Λ_X | λ_E(τ_Q) hüküm bantlarında | λ_X(τ_Q) |
|---|---|---|---|---|
| son | 1.6717 | 1.4190 | 1.641 → 2.406 | 2.436 → 1.924 |
| Hkeskin | 1.7132 | 1.4200 | 1.666 → 2.294 | 2.536 → 1.902 |
| HA4 | 1.6088 | 1.3951 | 1.613 → 2.654 | 2.219 → 1.867 |

(`τ ≤ taban`'da `h' = 0` olduğu için `λ_E` orada tanımsızdır;
merdiven ortalamaları `|h'|²`/`|y|²` ağırlığıyla ve yalnız
`h', y ≠ 0` çizgileri üzerinden alınmıştır.)

### 3e. TÜRETİMİN AÇIK BOŞLUĞU (gizlenmiyor)

Naif bacak-DW muhasebesi Ç1(B) ailesinde (`q₂ = Q, ε₂ = +`;
`q₃ = q₁, ε₃ = −ε₁`) şunu söyler: gerçek üçlü çarpımın
X̃-ortalaması `Φ₀(ΣA)`, modelinki ise `ΠΦ₀(A_i)`; Gauss limitinde
oran `exp(+A₁²σ²) > 1`, yani **bir KAZANÇ**. Ölçülen ise bir
**SÖNÜM** (`kalib ≈ 0.2–0.3`). Dolayısıyla:

> **Kalibrasyonun BÜYÜKLÜĞÜ bacak-DW muhasebesinden ÇIKMIYOR.**
> Çıkan şey BANT ŞEKLİDİR: taşıyıcı bir X-bacağında oturduğu için
> (165 §2c'nin (B)+(C) ailesi — ve `T_a ≡ 0` olduğundan Ç1'in
> TAMAMI odur) bandın τ_Q'suna bağlı çarpan X-bacağının kendi
> sönümüdür. Büyüklüğü taşıyan şey Ç4'ün ortak kazancıdır ve o
> kapalı formda hâlâ yok (165'in 1. sıradaki açık borcu).

---

## 4. T2 — HİPOTEZ YARIŞI

### 4a. Kurallar

* **Hedef:** `KALİB_u2` (birincil; bağıl jackknife %0.35–6.5) ve
  `KALİB_s2` (bağıl jackknife %1.2–257).
* **Ölçek sınıfı:** `S0` = sıfır ölçek (yasa doğrudan kalib);
  `S1` = ÜÇ GAZ İÇİN TEK global `c`. **Gaz-başına serbest ölçek YOK.**
* **Hakemler:** (i) log-artık rms/maks ve `χ²/dof` (jackknife);
  (ii) çapraz-gaz — her gazın tek başına isteyeceği `c_g`'nin
  maks/min yayılımı; (iii) artık yapısı — log-artığın τ eğimi.

### 4b. Sonuç (20 hüküm bandı, hedef `KALİB_u2`)

| yasa | sınıf | c | **rms %** | maks % | χ²/dof | **çapraz-gaz** | eğim |
|---|---|---|---|---|---|---|---|
| **K1d** `c·W_amp(τ)` | S1 | 0.3475 | **9.78** | 33.2 | 86 | **1.164** | −0.57 |
| **K1p** `c·√(W_X W_amp)` | S1 | 0.3646 | 10.10 | 34.1 | 89 | 1.187 | −0.41 |
| K1s `c·e^{−2π²σ_X̃²τ²}` (Gauss) | S1 | 0.3818 | 10.44 | 34.7 | 94 | 1.208 | −0.28 |
| K1t `c·e^{−½π²σ_ds²τ²}` (Gauss) | S1 | 0.3453 | 10.50 | 35.6 | 100 | 1.176 | −0.60 |
| **K1c** `c·W_X(τ)` | S1 | 0.3825 | 10.63 | 35.0 | 97 | 1.212 | **−0.26** |
| K1r `c·e^{−2π²σ_Ĉ²τ²}` (Gauss) | S1 | 0.4488 | 10.63 | 31.2 | 88 | 1.236 | +0.22 |
| **K1b** `c·W_pos(τ)` | S1 | 0.4606 | 10.77 | 29.2 | 89 | 1.236 | +0.38 |
| K1q `c·W_X^{3/2}` | S1 | 0.4877 | 11.14 | 28.7 | 93 | 1.242 | +0.49 |
| **K1o** `c·W_X·W_amp` **(X-bacağın tamamı)** | S1 | **0.5647** | 11.42 | **23.0** | 103 | 1.225 | +0.92 |
| K1n `c/(λ_Xλ_E)` | S1 | 0.9035 | 14.05 | 30.0 | 348 | 1.245 | −0.99 |
| K3f `c·W_pos/λ_X` | S1 | 0.9180 | 14.33 | 35.7 | 270 | 1.299 | −0.50 |
| K1e `c·W_pos·W_amp` | S1 | 0.6800 | 14.79 | 30.6 | 210 | 1.233 | +1.55 |
| K2f `c·√ρ` (0.40-tb) | S1 | 0.3618 | 16.11 | 39.0 | 304 | 1.244 | +1.09 |
| **K0 sabit (bant bağımsız)** | S1 | 0.2354 | **16.98** | 47.5 | 312 | 1.197 | −1.74 |
| K2e `c·√ρ` (0.52-tb) | S1 | 0.4764 | 21.53 | 65.2 | 422 | 1.271 | +1.48 |
| K1m `1/(λ_X Λ_E Λ_X)` | **S0** | 1 | 27.49 | 56.9 | 1085 | 1.330 | −2.57 |
| K0 `g = g_E g_X²` | **S0** | 1 | 31.00 | 78.0 | 882 | 1.287 | −1.68 |
| **K2d** `c·ρ(τ)` **(0.40-tb)** | S1 | 0.5561 | **43.90** | 125.5 | 2240 | **1.488** | +3.92 |
| K1g `g·W_pos` | **S0** | 1 | 44.94 | 59.4 | 2433 | 1.345 | +0.45 |
| **K2b** `c·ρ(τ)` **(0.52-tb)** | S1 | 0.9641 | **55.99** | 177.8 | 3471 | **1.726** | +4.69 |
| **K2a** `ρ(τ)` **(0.52-tb)** | **S0** | 1 | **56.11** | 174.2 | 4020 | 1.726 | +4.69 |
| K1k `1/(λ_Eλ_X²)` | **S0** | 1 | 62.99 | 112.0 | 4273 | 1.307 | −1.88 |
| **K2c** `ρ(τ)` **(0.40-tb)** | **S0** | 1 | **73.28** | 114.9 | 31670 | 1.488 | +3.92 |
| K1a `W_pos(τ)` | **S0** | 1 | 78.26 | 106.7 | 25353 | 1.236 | +0.38 |
| K1h `1/λ_X(ω_Q)` | **S0** | 1 | 79.47 | 129.8 | 18212 | 1.257 | −2.63 |

(Yasaların tamamı — 34 aday — `scratchpad/166/yaris_cikti.txt` içinde;
27 bantlık geniş şekil tanısı ve `KALİB_s2` hedefli iki koşu da orada.)

> **HÜKÜM (T2).**
> **(a) H-K2 DÜŞTÜ.** ρ'nun her dört sürümü de bant şeklinde
> **%44–73** artık bırakıyor — bant-bağımsız SABİTTEN (%17) bile
> KÖTÜ. Çapraz-gaz yayılımı 1.49–1.73 (K1 ailesi 1.16–1.25).
> Artık eğimi +3.9 … +4.7: ρ **çok dik** iniyor. KALEM'in kokusu
> (0.15–0.31 ↔ 0.38/0.27/0.16/0.08) bir **aralık çakışmasıydı**,
> yasa değil.
> **(b) H-K1 KAZANDI.** Bütün Debye–Waller üyeleri %9.8–11.4 ile
> önde; hepsi aynı aileden ve hüküm penceresinde **birbirinden
> ayırt edilemiyor** (τ ∈ 0.52–0.78'de üç çekirdek de yaklaşık aynı
> Gauss'a benziyor). Türetim (§3a) **X-bacağının tamamını**
> (`W_amp·W_X`) işaret ediyor: taşıyıcı Ç1'in (B)+(C) ailesinde bir
> X-bacağında oturuyor. Bu üye aynı zamanda **en küçük maksimum
> sapmaya** sahip (%23.0).
> **(c) H-K3 (çarpımlar) HİÇBİR YERDE saf K1'i geçmedi.** ρ'yu
> katmak daima artığı büyüttü (K3a/K3b: %56–67).
> **(d) SIFIR ÖLÇEK OLMADI.** En iyi S0 yasası %27.5 (K1m).

### 4c. Çapraz-gaz: kalan borç

| yasa | c(`son`) | c(`Hkeskin`) | c(`HA4`) | yayılım |
|---|---|---|---|---|
| `c·W_X` | **0.4053** | **0.4051** | 0.3345 | 1.212 |
| `c·W_amp` | 0.3618 | 0.3647 | 0.3133 | 1.164 |
| `c·W_X·W_amp` | 0.5822 | 0.6100 | 0.4980 | 1.225 |
| `c·ρ` (0.52-tb) | 1.090 | 0.702 | 1.211 | 1.726 |

> `W_X` üyesinde **gerçek gaz ile sadakatli-keskin gaz aynı ölçeği
> istiyor: 0.4053 ↔ 0.4051, yani %0.06 içinde.** Yani yasa gerçek ↔ sentetik
> arasında SERBEST PARAMETRESİZ taşınıyor. Kalan kırılma yalnız
> `HA4`'tedir (−%17.5) ve §3c'nin ölçtüğü inşa farkıyla aynı yöne
> bakıyor; ama denenen düzeltmeler (`R_h/W_amp` ve tersi) `Hkeskin`'i
> bozduğu için **kapanmadı**.

---

## 5. T3 — KAPANIŞ: 165'in T1 TABLOSU, u2-SABİTLEMESİZ

**Kural:** yasanın tek ölçeği `c` YALNIZ `KALİB_u2`ye (ikinci
momentlere) oturtuldu; üçüncü momentten hiçbir girdi yok. Sonra
`kalib_yasa · A²s2_öngörü` ile ölçülen karşılaştırıldı.

### 5a. Sayım

| yasa | sınıf | c | **±%25** | son | Hkeskin | HA4 |
|---|---|---|---|---|---|---|
| **K1o `c·W_X·W_amp`** | S1 | 0.5647 | **13/20** | 5/7 | 5/7 | 3/6 |
| K1e `c·W_pos·W_amp` | S1 | 0.6800 | **13/20** | 5/7 | 5/7 | 3/6 |
| K1f `c·W_pos²` | S1 | 0.9014 | **13/20** | 5/7 | 5/7 | 3/6 |
| K1u `c·W_pos·W_X` | S1 | 0.7486 | **13/20** | 5/7 | 5/7 | 3/6 |
| K2e `c·√ρ` (0.52-tb) | S1 | 0.4764 | 13/20 | 5/7 | 3/7 | 5/6 |
| **[165] u2-normalizasyonu** | (bant bant ölçülen) | — | **12/20** | 5/7 | 3/7 | 4/6 |
| K1b `c·W_pos` | S1 | 0.4606 | 12/20 | 5/7 | 4/7 | 3/6 |
| K1q `c·W_X^{3/2}` | S1 | 0.4877 | 12/20 | 5/7 | 4/7 | 3/6 |
| K2f `c·√ρ` (0.40-tb) | S1 | 0.3618 | 12/20 | 4/7 | 3/7 | 5/6 |
| K1c `c·W_X` / K1d `c·W_amp` | S1 | 0.3825 / 0.3475 | 11/20 | 4/7 | 4/7 | 3/6 |
| K2d `c·ρ` (0.40-tb) | S1 | 0.5561 | 10/20 | 6/7 | 1/7 | 3/6 |
| K2b `c·ρ` (0.52-tb) | S1 | 0.9641 | 8/20 | 6/7 | 1/7 | 1/6 |
| K0 sabit | S1 | 0.2354 | 7/20 | 3/7 | 2/7 | 2/6 |
| K2a `ρ` (0.52-tb) | **S0** | 1 | 7/20 | 6/7 | 0/7 | 1/6 |
| K2c `ρ` (0.40-tb) | **S0** | 1 | **0/20** | 0/7 | 0/7 | 0/6 |
| K1a `W_pos` | **S0** | 1 | **0/20** | 0/7 | 0/7 | 0/6 |

### 5b. Kazanan yasayla bant bant (`kalib = 0.5647·W_X·W_amp`)

| gaz | τ_eff | kalib_yasa | öngörü×kalib | ölçüm | **oran** | ±%25 |
|---|---|---|---|---|---|---|
| **son** | 0.5392 | 0.3245 | +0.113456 | +0.136666 | **0.830** | ✓ |
| son | 0.5791 | 0.2982 | +0.166368 | +0.178412 | **0.932** | ✓ |
| son | 0.6188 | 0.2725 | +0.224633 | +0.229560 | **0.979** | ✓ |
| son | 0.6578 | 0.2473 | +0.272760 | +0.255157 | **1.069** | ✓ |
| son | 0.6980 | 0.2232 | +0.314567 | +0.280278 | **1.122** | ✓ |
| son | 0.7382 | 0.1997 | +0.321100 | +0.211839 | 1.516 | |
| son | 0.7766 | 0.1786 | +0.322418 | +0.076314 | 4.225 | |
| **Hkeskin** | 0.5393 | 0.3076 | +0.110817 | +0.144047 | **0.769** | ✓ |
| Hkeskin | 0.5792 | 0.2801 | +0.161513 | +0.184278 | **0.876** | ✓ |
| Hkeskin | 0.6189 | 0.2534 | +0.212918 | +0.222319 | **0.958** | ✓ |
| Hkeskin | 0.6581 | 0.2274 | +0.251809 | +0.238466 | **1.056** | ✓ |
| Hkeskin | 0.6983 | 0.2028 | +0.276817 | +0.242369 | **1.142** | ✓ |
| Hkeskin | 0.7386 | 0.1790 | +0.263129 | +0.169534 | 1.552 | |
| Hkeskin | 0.7774 | 0.1578 | +0.238456 | +0.011062 | 21.56 | |
| **HA4** | 0.5392 | 0.3085 | +0.117731 | +0.126393 | **0.931** | ✓ |
| HA4 | 0.5788 | 0.2810 | +0.172677 | +0.169768 | **1.017** | ✓ |
| HA4 | 0.6181 | 0.2544 | +0.239087 | +0.219593 | **1.089** | ✓ |
| HA4 | 0.6566 | 0.2283 | +0.314952 | +0.251604 | 1.252 | *(binde 2 kaçtı)* |
| HA4 | 0.6962 | 0.2037 | +0.420744 | +0.267287 | 1.574 | |
| HA4 | 0.7349 | 0.1799 | +0.577710 | +0.154810 | 3.732 | |

**Sayım: 13/20; τ_eff ∈ 0.54–0.70 penceresinde 13/15** (165: 12/20 ve
12/15). Kaçan yedi bandın **beşi** τ_eff ≥ 0.73'tedir — 165 §3c'nin
zaten kayıtlı çöküşü, bu koşu onu düzeltmedi. Kalan iki kaçak
`HA4`'ün 0.657 (1.252, sınırı binde ikiyle) ve 0.696 (1.574)
bantlarıdır; ikisi de §4c'nin `HA4` ölçek kırılmasıyla aynı yöne
bakıyor.

### 5c. Ne ilan edilebilir, ne edilemez

> **EDİLEBİLİR:** Üçüncü-moment öngörüsünün **bant bant değişen
> kalibrasyonu ortadan kalktı.** 165 her bant için ayrı bir
> `n_u2(τ_bant)` ölçüyordu (20 sayı); 166 bunu **1 sayı + ölçülen bir
> şekil fonksiyonuyla** değiştiriyor. Şekil fonksiyonunun girdileri
> (`σ_ds`/`σ_X̃`, daha doğrusu `ds` ve `X̃`'nin TAM karakteristik
> fonksiyonları) gazın kendi ikinci-moment/marjinal ölçümleridir ve
> 142–150 zincirinde makineden gelir. **Üçüncü momentten hiçbir girdi
> yoktur** ve sonuç 165'ten İYİDİR (13/20 ↔ 12/20).
>
> **EDİLEMEZ:** "Tamamen parametresiz." Geriye **bir evrensel sayı**
> (`c = 0.5647`) kaldı, türetilmedi, ve üç gazda −%12 … +%8 oynuyor
> (`son`/`Hkeskin` %4.8 içinde, `HA4` `son`'dan −%14; `W_X` üyesinde
> `son`/`Hkeskin` farkı %0.06'ya düşüyor, `HA4` −%17.5 kalıyor). Dürüst çerçeve:
> **"tek evrensel ölçekli, bant-şekli parametresiz."**

---

## 6. T4 (bonus) — τ ≳ 0.57 GÖSTERİM ŞİŞMESİ

165 §7: `Var(E_mod)/Var(η) = 2.1–3.7`. Soru: kazanan yasadan çıkıyor mu,
`1/ρ` tipi mi?

**Türetilen özdeşlik.** `h'_q = 2⟨η e^{−iω_q s}⟩` olduğundan
`⟨E_mod·η⟩ = ½Σ_q|h'_q|²` ÖZDEŞ; ve `⟨E_mod²⟩ = ½Σ_q|h'_q|²Re λ_E(ω_q)`.
Dolayısıyla
```
Var(E_mod)/Var(η) = Λ_E · Γ_E ,
Λ_E = Σ|h'|²Reλ_E / Σ|h'|²   (GRAM SIZINTISI, sentez tarafı)
Γ_E = ½Σ|h'|² / Var(η)       (GRAM ATFETMESİ, ölçüm tarafı)
ve  Λ_E ≡ 1/g_E   (165'in alan-regresyonu kazancı)
```

| gaz | şişme (E) | `Λ_E·Γ_E` | `Λ_E` = `1/g_E` | `Γ_E` | şişme (X) | `Λ_X·Γ_X` | `1/W̄_pos²` |
|---|---|---|---|---|---|---|---|
| son | **2.139** | 2.180 | 1.672 = 1.672 | 1.304 | **1.585** | **1.585** | 1.788 |
| Hkeskin | **2.142** | 2.199 | 1.713 = 1.713 | 1.284 | **1.562** | **1.562** | 1.823 |
| HA4 | **2.088** | 2.103 | 1.609 = 1.609 | 1.307 | **1.520** | **1.520** | 1.570 |

> **HÜKÜM (T4).** `X` kanalında özdeşlik **üç hanede TAM** kapanıyor
> (`Λ_X·Γ_X` = şişme, üç gazda da). `E` kanalında %2 fark var ve
> sebebi biliniyor: `g_E` ortalaması çıkarılmamış `E` ile hesaplanıyor.
> **Şişme `1/ρ` tipi DEĞİLDİR:** aynı bantlarda `1/ρ` 1.86–39.0 arasında
> ve bantla hızla değişiyor, şişme ise merdiven-genelinde TEK bir sayı
> (2.09–2.14) ve üç gazda **%2.5 içinde aynı** — oysa aynı gazların
> ρ'ları birbirinden 2–20 kat farklı. Şişmenin kimliği
> **Gram (dik olmayan merdiven) sızıntısı + atfetmesidir**;
> 146'nın "ρ da r gibi Gram-atfetme gölgesidir" dersinin nicel karşılığı.
> `1/W̄_pos²` aynı mahallede (1.57–1.82 ↔ 1.52–1.61 X kanalı için)
> ama **mühürsüz** — kapalı form Λ ve Γ'nın kendisidir.

---

## 7. DENETİM

| | sınav | sonuç |
|---|---|---|
| **V1** | 166'nın T1 zinciri = 165'in kayıtlı JSON'u mu? | üç gaz, **27 ölçülen bant**, 13 büyüklük (`A²s0..s3`, `A²u2`, öngörüleri, `norm_u2`, `τ_eff`, `Ç1`, `Ç2`, `Ç3`): maks fark **0.0e+00** (istenen 1e−9'un ötesinde: **bit düzeyi**) |
| **V-R1** | 166'nın ρ estimatörü = 144'ün sayıları mı? | `son`, 0.52-taban, 144-bantları: **0.3806 / 0.2664 / 0.1589 / 0.0796** ↔ 144: 0.381/0.266/0.159/0.080 |
| **V-B1** | E-bacağı kalem formu (`h'/b = W_amp`) | `son`, τ ∈ 0.42–0.60: **±%2.5** (τ = 0.579'da %0.2) |
| **V-B2** | X-bacağı kalem formu (`y = ½bW_amp(W_X e^{−2πiτ}+1)`) | `son`/`Hkeskin`, τ ≤ 0.38: **%1–5**; çıplak limitte 165 §2a'yı özdeş veriyor; τ = ½ kutbunu her iki taraf da üretiyor |
| **V-B3** | `h' ≡ 0`, `τ ≤ taban` | üç gazda `R_h = 0.0000` (165'in `|hp|/b = 0.0000` kaydı) |
| **V-T4** | şişme özdeşliği `Var(X_mod)/Var(X̃) = Λ_X Γ_X` | 1.585/1.562/1.520 ↔ 1.585/1.562/1.520 — **üç hanede tam**; `Λ_E ≡ 1/g_E` de özdeş |
| **V-J** | jackknife tutarlılığı | `KALİB_u2` bağıl jackknife hatası **%0.35–6.5** (τ_eff ≤ 0.70'te ≤ %1.4), `KALİB_s2` **%1.2–257** (yüksek τ'da patlıyor — 165 §3c ile tutarlı) |

---

# HÜKÜM

## (i) H-K2 (ρ) DÜŞTÜ — ve nasıl düştüğü kayıtlı

> ρ üç gazda, iki η konvansiyonunda, iki bant kümesinde ölçüldü ve
> 144 birebir doğrulandı. Kalibrasyon ile bant şekli **uyuşmuyor**
> (%44–73 artık, sabit yasadan bile kötü), çapraz-gaz yayılımı 1.49–1.73,
> artık eğimi +3.9…+4.7. **Katil kanıt `HA4`'tür:** erfc penceresi
> merdiveni τ > 0.68'de sıfırlıyor, ρ oraya **0.026 → −0.002** iniyor,
> kalibrasyon ise **0.189 / 0.146**'da duruyor. KALEM'in kokusu bir
> aralık çakışmasıydı.

## (ii) H-K1 (BACAK SÖNÜMÜ) KAZANDI — bant şekli artık PARAMETRESİZ

> Türetim açık: **E-bacağı `A = πτ` ile `ds`'nin**, **X-bacağı ayrıca
> `A = 2πτ` ile `X̃`'nin** Debye–Waller çarpanını taşır; taşıyıcı
> ölçüm ve öngörüde aynı sitede olduğu için doğrudan bir çarpan
> taşımaz. Kalem formu (`y_q = ½b_q W_amp[W_X e^{−2πiτ}+1]`) çıplak
> limitte 165 §2a'yı verir ve düşük τ'da ölçümle %1–5 uyuşur.
> Ç1'in taşıyıcısı (B)+(C) ailesinde bir X-bacağında oturduğundan
> bandın çarpanı **X-bacağının kendi sönümü** `W_amp(τ)W_X(τ)`'dir:
> 20 bantta rms **%11.4**, maks **%23.0**, çapraz-gaz **1.225**.

## (iii) T3: 12/20 → **13/20** — VE BANT-BANT KALİBRASYON KALKTI

> `kalib(τ) = 0.5647 · W_amp(τ)·W_X(τ)`, ölçeği yalnız ikinci
> momentlere oturtulmuş, üçüncü momentte **13/20** (sağlıklı
> pencerede 13/15). 165'in yirmi ayrı `n_u2` ölçümü **bir sayıya**
> indi ve sonuç kötüleşmedi, iyileşti.

## (iv) T4: ŞİŞME 1/ρ DEĞİL, GRAM

> `Var(E_mod)/Var(η) = Λ_E·Γ_E` (sızıntı × atfetme); X kanalında
> özdeşlik üç hanede tam, `Λ_E ≡ 1/g_E`. 2.1'in kimliği budur.

## (v) AÇIK BORÇ — "PARAMETRESİZ" HENÜZ İLAN EDİLEMEZ

> Bir evrensel sayı (`c ≈ 0.56`) kaldı; türetilmedi ve `HA4`'te
> %14 kırılıyor. Naif bacak-DW muhasebesi Ç1(B)'de bir kazanç
> öngörüyor, ölçülen ise sönüm — yani **büyüklük Ç4'ün ortak
> kazancından gelmeli** ve o kapalı form hâlâ yok (165'in 1. sıradaki
> borcu). Bu koşu şekli kazandı, ölçeği kazanamadı.

---

## Sıradaki adım (bu ölçümün işaret ettiği)

1. **`c ≈ 0.56`'nın kaynağı = Ç4'ün ortak kazancı.** 165'in "Sıradaki
   adım 1"i (beşli rezonans `Δ = ±ω_r`'nin asal-kuvvet çözümlerini
   sayıp `Σ𝒢_r h'yy`'yi kapalı yazmak) artık YALNIZ BİR SAYI için
   gerekiyor — bant şekli çözüldü. Hedef daralınca iş kolaylaştı.
2. **`HA4`'ün −%15'i.** §3c ölçtü: 164'ün sadakatli inşası çizgileri
   yarım-gap formülünden daha koherent gerçekleştiriyor
   (`R_h`: son 0.80 / Hkeskin 0.87 / HA4 0.97 aynı τ'da). Bu farkın
   `c`'ye nasıl geçtiği, gaz inşasının kendi sınavıdır (164'e geri
   besleme).
3. **`W_amp` / `W_X` / `W_pos` ayrımı.** Üçü hüküm penceresinde
   ayırt edilemiyor. Ayırıcı, daha GENİŞ bir τ kolu: `taban = 0.28`
   ya da `0.34` η'sıyla (155'te önbellekleri var) 0.30–0.50 bantları
   açılırsa üç çekirdeğin eğrilikleri ayrışır.
4. **`KALİB_s2` ile `KALİB_u2`'nin farkı bir SİSTEMATİKTİR.**
   `KALİB_s2` daha dik iniyor (τ ≥ 0.74'te çöküyor). Sağlıklı
   pencerede bu farkı kapatan üye `c·W_pos²` (rms %12.4) — yani
   üçüncü moment ikinciden **bir çarpan fazla** sönüm istiyor.
   165'in "Sıradaki adım 3"üyle (τ ≥ 0.74'te ölçülen sönerken öngörü
   büyüyor) aynı olgunun iki yüzü olabilir.
5. **Koşulmayanlar:** `orta`, `dusuk`, `P1` ve 164'ün yedi açık gazı;
   `taban` ekseni yalnız 0.40 (ρ için ayrıca 0.52); `τ_c` yalnız 0.95;
   ızgara yalnız `t1`; `kaynak` yalnız `olculen`.

---

## Dürüstlük notları

* **Uydurma yok; her sayı bu koşuda üretildi.** Tablolar
  `scratchpad/166/`'daki `T1_<gaz>_olculen_tc0.95.json` (3),
  `RHO_<gaz>.json` (3), `BACAK_<gaz>.json` + `.npz` (3),
  `yaris_cikti.txt`, `artik_tablo.txt`, `girdiler.txt`, `log_*.txt`
  çıktılarından alındı. 165/144'ten alıntılanan tek şey referans
  satırlarıdır ve her biri V1 / V-R1'de bu koşunun kendi ölçümüyle
  karşılaştırılmıştır.
* **ÖLÇEK BİR SERBESTLİK DERECESİDİR ve saklanmıyor.** `c` ÜÇ GAZIN
  BÜTÜN bantlarına ortak tek sayıdır (gaz-başına ölçek YOK) ve
  YALNIZ `KALİB_u2`ye — yani ikinci momentlere — oturtulmuştur;
  üçüncü moment sınav tarafındadır. Yine de bir parametredir ve
  hükümde "parametresiz" değil **"tek evrensel ölçek"** diye anılır.
* **`χ²/dof` değerleri 86–31670 arasındadır ve hiçbiri "kabul"
  seviyesinde değildir.** Jackknife hataları %0.3–1 mertebesinde
  olduğundan sistematikler istatistik hatayı ezmektedir; bu yüzden
  hüküm log-artık rms / maks ve ÇAPRAZ-GAZ yayılımı üzerinden
  verilmiştir. `χ²` yalnız sıralama için kayda geçmiştir.
* **ρ TABAN-BAĞIMLIDIR (146 E2).** Bu yüzden iki tabanda birden
  koşuldu; **H-K2 her ikisinde de kaybediyor**, yani hüküm konvansiyon
  seçimine bağlı değildir. ρ'nun "mutlak" bir soğurma eğrisi olmadığı
  (Gram-atfetme gölgesi olduğu) kaydı ayakta.
* **`HA4`'te ρ ≈ 0 bir ölçüm kusuru değil inşa gerçeğidir**
  (erfc 0.68/0.125 merdiveni söndürüyor). Bu, H-K2 aleyhine bir kanıt
  olarak kullanıldı ve nedeni açıkça yazıldı.
* **Bacak türetimi TAMAMLANMADI (§3e).** Naif DW muhasebesi Ç1(B)'de
  kazanç veriyor, ölçülen sönüm. Yani **§3a bant ŞEKLİNİ açıklıyor,
  BÜYÜKLÜĞÜ açıklamıyor**; §3a'nın "taşıyıcı doğrudan sönüm taşımaz"
  hükmüyle "ama band çarpanı X-bacağının sönümüdür" ifadesi arasındaki
  köprü ÖLÇÜLMÜŞTÜR, TÜRETİLMEMİŞTİR.
* **`W_amp` / `W_X` / `W_pos` hüküm penceresinde ayırt edilemiyor**
  (rms %9.8 / %10.6 / %10.8). Galip üye türetim gerekçesiyle
  (`W_amp·W_X` = X-bacağının tamamı) seçildi, sayısal üstünlükle
  değil; `c·W_amp` tek başına rms'de biraz daha iyidir. Bu
  saklanmıyor.
* **`λ_E` `τ ≤ taban`'da tanımsızdır** (`h' = 0`); `BACAK_*.json`'daki
  `ozet` satırlarında orada 10³–10⁴ mertebesinde anlamsız sayılar
  görünür ve hiçbir hesaba girmemiştir. Merdiven ortalamaları
  `|h'|²`/`|y|²` ağırlığıyla ve yalnız sıfır-olmayan çizgiler
  üzerinden alınmıştır (`166_yaris.yukle` bunu yeniden hesaplar).
* **`K2e c·√ρ` de 13/20 veriyor** ama bant şekli rms'i %21.5 ve
  çapraz-gaz yayılımı 1.271'dir; ±%25 sayımı kaba bir ölçüttür ve
  tek başına hakem değildir. Hüküm üç hakemin birlikte okunmasıyla
  verilmiştir.
* **Figür üretilmedi;** bu koşuda görselleştirme istenmedi.
* **git'e dokunulmadı.**
* **Süreler.** ρ 1.0–1.1 dk/gaz; T1 (2276 frekans, ölçüm+öngörü+
  jackknife) 4.2–4.3 dk/gaz; bacak (karakteristik fonksiyonlar
  8981 çizgide + λ tayfı) 6.0–7.5 dk/gaz; yarış < 5 s. Üç gaz
  eşzamanlı; toplam ~25 dk duvar saati.

---

## Ek — dosyalar ve tekrar-üretim

| dosya | ne |
|---|---|
| `166_configs/166_T1.py` | 165'in T1 zinciri, ÇİZGİ BAZINDA kayıtla; `KALİB_s2`/`KALİB_u2` + bant-içi jackknife |
| `166_configs/166_rho.py` | 144'ün ρ estimatörü, iki taban × iki bant kümesi, üç gaz (V-R1 burada) |
| `166_configs/166_bacak.py` | bacak sönümlerinin TÜRETİMİ (docstring'de tam) + `W_amp`, `W_X`, `W_pos` tam karakteristik fonksiyonları + Gram sızıntıları `λ_E`, `λ_X` |
| `166_configs/166_yaris.py` | 34 aday yasa, üç hakem, T3 sayımı, T4 özdeşliği |

```
cd /Users/ugursezen/Desktop/arin/deney
for v in son Hkeskin HA4; do .venv/bin/python qm_riemann/166_configs/166_rho.py   $v; done
for v in son Hkeskin HA4; do .venv/bin/python qm_riemann/166_configs/166_T1.py    $v 0.40 olculen 0.95; done
for v in son Hkeskin HA4; do .venv/bin/python qm_riemann/166_configs/166_bacak.py $v; done
.venv/bin/python qm_riemann/166_configs/166_yaris.py
```

Ham çıktılar `scratchpad/166/`: `T1_<gaz>_olculen_tc0.95.json` (3),
`RHO_<gaz>.json` (3), `BACAK_<gaz>.json` + `BACAK_<gaz>.npz` (3),
`yaris_cikti.txt`, `artik_tablo.txt`, `girdiler.txt`,
`log_rho_*.txt`, `log_T1_*.txt`, `log_bacak_*.txt`.
