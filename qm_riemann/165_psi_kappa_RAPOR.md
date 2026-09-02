# 165 — Ψ⋆κ SINAVI: FAZ-KİLİDİ KATMANININ KANAL-ETİKETLİ, PENCERE-FAZLI, ADİL-ZEMİNLİ HESABI

**Hedef.** KALEM (2 Eylül) 163'ün patlamasını teşhis etmişti: kutu-pencere
yakın-rezonansları birim ağırlıkla topladı; doğru ağırlık pencere dönüşümü
`κ(ν)` × `cos(ν t̄ + faz)`'dır ve ölçülen moment `(Ψ⋆κ)(0)`'dır. Sınav
164'ün SADAKATLİ gazlarında (ilk-kök + `n−½`) koşuldu: `Hkeskin`, `HA4`
ve gerçek gaz `son`. Her `(q₁,q₂,q₃)` kombinasyonu **BUDAMASIZ** toplandı
ve **Ç1 / Ç2 / Ç3** etiketlerine ayrıldı — ve ölçüm görevin
hiyerarşisinde OLMAYAN dördüncü bir kanalı (Ç4, tarak-rezonansı) buldu.

> **Bu koşunun bir cümlelik özeti: SAYISAL HEDEF (±%25) BÜYÜK ÖLÇÜDE
> TUTTU (163: 4/27 → 165: 12/20 bant; τ_eff ∈ 0.54–0.70 penceresinde
> 12/15), ama KANAL HİYERARŞİSİ (H-Ç) DÜŞTÜ: Ç1 baskın DEĞİL (%21–41),
> Ç2 ölü (≤ 4×10⁻⁴), ve KALEM'in çare olarak gösterdiği Ç3 (pencere-fazlı
> yakın-rezonans) NEREDEYSE SIFIR (%0.2–1.4). Payın %66–78'ini
> hiyerarşide OLMAYAN dördüncü bir kanal taşıyor: TARAK-REZONANSI (Ç4),
> yani `Δ ≈ ±ω_r`'de `κ`'nın 143'ün 𝒢-tepeleri.**

---

## Kısa hüküm

1. **ÖLÇÜM ZİNCİRİ BİT DÜZEYİNDE (V1).** 165'in bant birleştirmesi 160'ın
   kayıtlı JSON'unu `son`'un dokuz bandında `A²s2`, `A·s1`, `A³s3`,
   `A·κ1`, `τ_eff`, çizgi sayısında **0.0e+00** farkla veriyor
   (`s0`'da 2.6e−16).
2. **HESAP TEK SİTEYE İNDİRGENDİ ve BUDAMASIZ KAPANDI.** Bütün alanlar
   `s_n ≡ m_{n+1}` sitesinde yazılınca 163'ün `Γ`/`Φ` kaydırma çekirdeği
   (birinci-mertebe kesmenin yeri) **tamamen kalkıyor**. Kalan özdeşlik
   `κ(Δ) = ⟨e^{iΔs_n}⟩` sayesinde toplam ve ortalama yer değiştiriyor:
   **(2·8981)³ ≈ 5.8×10¹² terim `O(N·n_çizgi)`'de TAM toplanıyor.**
   Doğrudan terim-terim sayımla özdeşlik **2.4e−09** bağıl farkla
   doğrulandı (V4; 19 044 çift × 138 işaretli çizgi, 3000 site).
3. **T1 SAYISAL: HEDEF TUTTU.** İkinci-moment normalizasyonuyla
   (`u2 = ⟨ρX̃²⟩/⟨ρ⟩` ölçüleni tutacak tek skaler; üçüncü momentten
   HİÇBİR girdi yok) `A²s2` öngörü/ölçüm oranı — **ölçülü çizgilerle**:
   `son` 0.78 / 0.89 / **0.98** / **1.11** / **1.21**;
   `Hkeskin` 0.72 / 0.85 / **0.98** / **1.13** / 1.30;
   `HA4` 0.80 / 0.88 / **0.99** / **1.15** / 1.46
   (τ_eff = 0.539 … 0.698). **±%25 içinde 12/20 bant** (163: 4/27),
   τ_eff ∈ 0.54–0.70'te **12/15**. Nominal (tamamen sıfır-parametre)
   sürüm: **11/20**, `HA4`'te 5/6.
4. **KANAL PAYLARI — H-Ç SINAVI DÜŞTÜ.** Üç gazda, yedi bantta:
   **Ç1 = %21–41**, **|Ç2|/tot ≤ 4×10⁻⁴** (yani ölü), **Ç3+Ç4 = %59–79**.
   Ç1'in payı bantla ARTIYOR (0.21 → 0.41).
   > "Ç1 çift-iptalli koherent kanal baskındır" hipotezi **düştü.**
5. **KALEM'in Ç3'Ü DE ÖLÜ ÇIKTI — ve yerine DÖRDÜNCÜ KANAL ÖLÇÜLDÜ.**
   Düzgün-tarak kontrolü (κ = Dirichlet ⇒ yalnız Ç1+Ç2+Ç3 yaşar) toplamın
   **%66–78'inin** kaybolduğunu gösteriyor: **Ç4 (tarak-rezonansı)**.
   Alt-merdiven terim sayımı bunu bağımsız doğruluyor: `|Δ| ≤ 200·dres`
   olan BÜTÜN terimlerin toplamı alt-merdiven toplamının yalnız
   **%3'ü**dür. **Ç3 ≈ %0.2–1.4.**
6. **κ ÖLÇÜLDÜ ve 143'ÜN 𝒢'SİNİ VERDİ.** `κ(ω_q)` ile
   `𝒢 = −πτ_q a_q cos(πτ_q)` sekiz düşük çizgide **%0.3–6** içinde
   uyuşuyor (ör. q=2: −0.03996 ↔ −0.04008). `κ̃`'nin merkez-fazsız
   zarfı ise Dirichlet: `Δ = k·dres`'te ÖZDEŞ 0, `k+½`'de 0.6366 = 2/π.
   **163'ün kutu penceresi tam olarak bu yapıyı 1 sayıyordu.**
7. **ÇİZGİ GÖSTERİMİNİN SINIRI ÖLÇÜLDÜ (yeni sistematik).** Sabit
   yarım-gap genlikli çizgi modeli `τ ≳ 0.57`'de bozuluyor: `Var(E_mod)`
   ölçülen `Var(η)`'yı τ≤0.95'te **2.1–3.7 kat** aşıyor. İki bağımsız
   ölçüt aynı yeri veriyor — faz seğirmesi `τ_c = 1/(2πσ_Ĉ) = 0.573 /
   0.578 / 0.583` (Hkeskin/HA4/son; 163 §6d'nin σ_Ĉ ≈ 0.27'si) ve
   yarım-gap modülasyonu `τ_c = 1/(πσ_ds) = 0.73–0.78`. **DÜZGÜN
   tarakta aynı model `u0 = 1.000–1.06` veriyor**, yani şişme bir
   sızıntı kusuru değil TARAK etkisidir.
8. **T3 — 163'ÜN A4 İŞARET ANOMALİSİ YOK OLDU.** Sadakatli gazlarda
   `Hkeskin`, `HA4` ve `son`'un τ-dilim matrisleri **hücre hücre aynı
   işarette ve ±0.01 içinde**. 163 §6b'nin "aynı hücrede son +0.0216,
   A4 −0.0554" karşıtlığı 152'nin inşa kusurundan geliyormuş (164'ün
   hükmüyle tutarlı). İşaret hâlâ **alçak-τ çift bloğunda** doğuyor:
   `(0–0.2)² + (0–0.2)×(0.2–0.4) + (0.2–0.4)²` toplamın **%40–51**'i.
9. **T2 — F(τ_q) ÇIKARILDI ve YARI-ANALİTİK FORMU İKİ KATSAYIYLA
   OTURDU.** Ç1'in çizgi-başına yoğunluğu KAPALI FORMDA
   `F(τ_q) = Re[h'_q ȳ_q]` ve **τ_band'dan BAĞIMSIZ**. Çıplak
   `F/a² = sin²(2πτ)` tek başına artığın %44'ünü bırakıyor;
   **`F/a² = 0.546 sin²(2πτ) + 0.666 sin⁴(πτ)` artığı %3.9'a indiriyor**
   (aynı yarım-gap ailesi). İkinci terim 163'ün `τ_Q = ½`'teki SAHTE
   sıfırını kaldırıyor — 163 §5(iii)'ün "kinematik kök yanlış yerde"
   kaydının çözümü budur.
10. **T4 — GAZ HL'Yİ DUYMUYOR (üst sınırla).** Gerçek gazın Ç3+Ç4 payı,
    saf-asal-merdiven gazlarınınkiyle **yedi bandın altısında
    KUŞATILMIŞ** (ör. τ_eff=0.619: HA4 0.7402 > son 0.7342 > Hkeskin
    0.7277). Fark τ_eff ≤ 0.70'te **≤ 0.014**, yani **⟨σX̃²⟩'nin
    %1.4'ünden az**. Bu, aritmetik-incelik (HL) katmanının bu ölçüdeki
    **üst sınırıdır**; ölçülen bir HL sinyali YOKTUR.

---

## 1. Yöntem

### 1a. Ölçüm zinciri kopyalanmadı

`165_configs/165_cekirdek.py`, `163_cekirdek`'i (o da 160/159/156/155/154'ü)
**import eder**. `gaz`, `bant_adaylari` (160'ın aday listesi, 220-örnekleme
tohumu 21, `gap < 2.5·dres` filtresi, ara-nokta `W' = w + gap/2`, 8-grup
round-robin) ve `olc_cizgi` **aynen** kullanılır. Değişen yalnız ÖNGÖRÜ
tarafıdır. **V1 bunu 160'ın kayıtlı JSON'una karşı `0.0e+00` ile sınar.**

### 1b. Koşulan gazlar, ızgara, taban

| gaz | ne | pencere |
|---|---|---|
| **`Hkeskin`** | 164'ün SADAKATLİ keskin gazı (`N̄+S = n−½`, sıralı ilk kök) | — |
| **`HA4`** | 164'ün SADAKATLİ erfc gazı (`n−½`, ilk kök) | erfc 0.68/0.125 |
| **`son`** | gerçek zeta sıfırlarının son 300k'sı | — |

`taban = 0.40`, `cap = 4000`, ızgara `IZGARA_T1 = izgara(0.44,0.80,0.04)`
(160'ın bantları); hüküm **`lo ≥ 0.52`** olan yedi bant üzerinden verilir
(görevin τ ∈ 0.52–0.80 kısıtı). Model merdiveni `τ ≤ τ_max = 0.95`
(8981 çizgi; 162 §2b'nin kafes-alias uyarısı korunuyor). Üç gazın
`N = 299 998`, `L = 12.0296`, `T = 156 693.7`, `dres = 4.010e−05`,
`s̄ = 1 054 305.2` değerleri özdeştir.

### 1c. Öngörünün iki sürümü (ikisi de üçüncü-momentten BAĞIMSIZ)

| sürüm | çizgi genlikleri | "sıfır parametre" mi? |
|---|---|---|
| **P0 (nominal)** | `h'_q = b_q·1{q>q_lim}`, `y_q = B_q e^{−iπτ_q}` (yalnız merdiven `a_q` + gazın penceresi) | **evet, tamamen** |
| **P0m (ölçülen)** | `h'_q = 2⟨η_{n+1}e^{−iω_q s}⟩`, `y_q = 2⟨X̃0 e^{−iω_q s}⟩` | evet (yalnız İKİNCİ momentler) |

`q_lim = min(e^{taban·L}, cap) = 122`: η zinciri (154.eta_zinciri)
`τ ≤ taban` çizgilerini regresyonla SİLER, bu yüzden nominal `h'` orada
sıfırdır (ölçülen `h'` de ölçülerek sıfır çıkıyor: `|hp|/b = 0.0000`).

---

## 2. TÜRETİM

### 2a. TEK SİTEYE İNDİRGEME (165'in kilit basitleştirmesi)

Ölçüm nesnesi 160 §2 / 163 §2b ile ÖZDEŞ:

```
c_{n+1} = η_{n+1} e^{−iω_Q m_{n+1}} ,  ⟨c⟩ = h_Q/2 ,  ⟨ρ⟩ = rm (ÖLÇÜLEN)
s_k = ⟨σ X̃0^k⟩/⟨ρ⟩ = Im[ h̄_Q J_k ] / (2 rm) ,   J_k = ⟨ c_{n+1} X̃0_n^k ⟩
```

163 `J_k`'yı `m_n` sitesinde yazdı ve `m_{n+1} = m_n + ḡ(1+X̃_n)` farkını
`Γ(θ,k)` (X̃'nin marjinali) ile taşımak zorunda kaldı — **birinci
mertebede kesilen yer orasıydı** (163 V4: `u0_pred` 1'den %1.6–10.5
sapıyordu). 165 bütün alanları **tek sitede, `s_n ≡ m_{n+1}`**'de yazar;
o zaman hiçbir kaydırma çekirdeği kalmaz:

```
η_{n+1} = Σ_q Re[h'_q e^{iω_q s_n}] + artık ,   h'_q = 2⟨η_{n+1}e^{−iω_q s_n}⟩
X̃0_n   = Σ_q Re[y_q  e^{iω_q s_n}] + artık ,   y_q  = 2⟨X̃0_n e^{−iω_q s_n}⟩
```

**ÇIPLAK DEĞERLER (163 §2a'dan tek satırda).** `ds_n = Σ_q b_q cos(ω_q m_n)`,
`b_q = 2a_q sin(πτ_q)` ⇒ `η_{n+1} ⊃ b_q cos(ω_q s_n)` ⇒ **`h'_q = b_q`
(REEL)**. `X̃0_n = Σ_q B_q cos(ω_q m_n + πτ_q)`, `B_q = a_q sin(2πτ_q)`;
`m_n = s_n − ḡ`, `ω_q ḡ = 2πτ_q` ⇒ **`y_q = B_q e^{−iπτ_q}`**.
Ölçüm bunu doğruluyor: `|y_q|/B_q = 0.995 / 0.993 / 0.982` ve
`arg y_q + πτ_q = −0.021 / −0.022 / −0.018` (q = 2,3,5; Hkeskin).

### 2b. BELİRLENİMLİ DÖRTLÜ TOPLAM — TAM, BUDAMASIZ

```
J_2 = (1/8) Σ_{(1,ε₁)} Σ_{(2,ε₂)} Σ_{(3,ε₃)}
          h'^{(ε₁)}_1 y^{(ε₂)}_2 y^{(ε₃)}_3 · κ(Δ)
Δ = ε₁ω₁ − ω_Q + ε₂ω₂ + ε₃ω₃      (x^{(+)}=x, x^{(−)}=x̄)
κ(ν) ≡ (1/N) Σ_n e^{iν s_n}       (PENCERE DÖNÜŞÜMÜ — TAM)
```

Dört merdiven frekansı (taşıyıcı `ω_Q` dahil) = "dördüncü kat". `κ` üç
şeyi birden taşır: (i) sonlu pencere zarfı, (ii) pencere MERKEZİ fazı
(`arg κ ≈ ν s̄`, `s̄ ≈ 1.05e6` ⇒ `ν ~ 1e−6` bile tam tur döner —
KALEM'in `cos(νt̄)`'si), (iii) tarağın kendi rezonans yapısı.

> **HESAP HİLESİ (ve 165'in neden budama yapmadığı).**
> `κ(Δ) = (1/N)Σ_n e^{iΔ s_n}` olduğundan toplam ile n-ortalaması yer
> değiştirir:
> ```
> J_2 = ⟨ [Σ_q Re(h'_q e^{iω_q s})]·e^{−iω_Q s}·[Σ_q Re(y_q e^{iω_q s})]² ⟩
> ```
> yani **bütün kombinasyonların κ-ağırlıklı toplamı = model alanlarının
> gerçek tarak üzerindeki ortalaması.** 5.8×10¹² terim `O(N·n_çizgi)`'de,
> **hiçbir kesme, hiçbir budama olmadan** toplanır. Sınavın "her
> kombinasyon" şartı böylece tam olarak yerine gelir.

### 2c. KANAL ETİKETLERİ

`q₁^{ε₁}q₂^{ε₂}q₃^{ε₃} = Q` denklemi asal kuvvetlerde (163 §2c'nin
teklik lemması) yalnız iki tür çözüm verir:

* **Ç1 — çift-iptal, ν ≡ 0 ARİTMETİKSİZ.** Dört frekans ikişerli birebir
  sadeleşir. Üç aile ve tam kapalı toplamları:
  ```
  (A) q₁=Q,ε₁=+ ; q₂=q₃=q, ε₂=−ε₃   → (h'_Q/4)·Σ_q |y_q|²         [= T_a]
  (B) q₂=Q,ε₂=+ ; q₃=q₁, ε₃=−ε₁     → (y_Q/4)·Σ_q Re[h'_q ȳ_q]
  (C) q₃=Q,ε₃=+ ; q₂=q₁, ε₂=−ε₁     → (B) ile aynı
  Ç1 = (A) + 2(B) − (üç çakışma terimi, mertebe a_Q³)
  ```
* **Ç2 — ν ≡ 0 ARİTMETİKLİ.** Yalnız aynı asalın kuleleri
  (`q_i = p^{m_i}`, `ε₁m₁+ε₂m₂+ε₃m₃ = m`), çift-iptal olanlar hariç.
* **Ç3 — ν ≠ 0, pencere-fazlı.** Geri kalan her şey;
  `Ç3 = TOPLAM − Ç1 − Ç2` (özdeş kapanış).

**OFF (ara-nokta) frekansında `W` merdiven çizgisi DEĞİLDİR ⇒ `Δ = 0`
çözümü yoktur ⇒ Ç1 = Ç2 = 0 ÖZDEŞ**; ara-nokta çıkarımının tamamı Ç3'e
gider. Bant birleştirmesi 160'ınkiyle birebirdir:
`A²s2 = Σ(p_on A_on² s2_on − p_off A_off² s2_off)/Σ(p_on−p_off)`.

### 2d. ÇIPLAK LİMİT = 163'ÜN KAPALI FORMU (V2)

`h' = b`, `y = B e^{−iπτ}`, `h_Q = b_Q`, `⟨ρ⟩ = b_Q²/4` koyunca:

```
Re[h'_q ȳ_q] = b_q B_q cos(πτ_q) = a_q² sin²(2πτ_q)
Im[h̄_Q y_Q]/(4⟨ρ⟩) = −½ sin(2πτ_Q) / Σ …
⇒  s2^{Ç1} = −½ sin(2πτ_Q) · Σ_q a_q² sin²(2πτ_q)      ← 163 §2e AYNEN
```

Sayısal doğrulama (V2, yedi τ_Q): kod ile kalem arasında maks fark
**4.8e−06** (bağıl ~3e−04; kalan çakışma terimlerinden).

---

## 3. T1 — ANA SINAV: (b)/(a) BANT BANT

Karşılaştırma nesnesi 160'ın `A²s2`'sidir. Öngörü hem on hem off
frekansta aynı formülle hesaplanır. **İKİNCİ-MOMENT NORMALİZASYONU
(NRM):** aynı dörtlü toplamın ρ-kanalı `u2 = ⟨ρX̃0²⟩/⟨ρ⟩` ölçüleni
tutacak biçimde tek skalerle ölçeklenir (`n_u2 = A²u2_ölç / A²u2_öng`).
Bu, 163'ün V4'ünün (`u0_pred ≡ 1`) aynı ailesinden bir kalibrasyondur ve
**üçüncü momentten hiçbir girdi taşımaz**; bütün kanallar aynı çarpanla
ölçeklendiği için **kanal payları normalizasyondan bağımsızdır.**

### 3a. ÖLÇÜLEN ÇİZGİLERLE (P0m), τ_c = 0.95

| gaz | τ_eff | **A²s2 ölç** | HAM öngörü | ham/ölç | `n_u2` | **NRM öngörü** | **nrm/ölç** | ±%25 |
|---|---|---|---|---|---|---|---|---|
| **son** | 0.5392 | +0.136666 | +0.349579 | 2.56 | 0.305 | +0.106696 | **0.781** | ✓ |
| son | 0.5791 | +0.178412 | +0.557829 | 3.13 | 0.285 | +0.159237 | **0.893** | ✓ |
| son | 0.6188 | +0.229560 | +0.824247 | 3.59 | 0.274 | +0.225767 | **0.983** | ✓ |
| son | 0.6578 | +0.255157 | +1.103086 | 4.32 | 0.258 | +0.284236 | **1.114** | ✓ |
| son | 0.6980 | +0.280278 | +1.409414 | 5.03 | 0.240 | +0.338720 | **1.209** | ✓ |
| son | 0.7382 | +0.211839 | +1.607653 | 7.59 | 0.221 | +0.354651 | 1.674 | |
| son | 0.7766 | +0.076314 | +1.805495 | 23.7 | 0.197 | +0.356046 | 4.666 | |
| **Hkeskin** | 0.5393 | +0.144047 | +0.360286 | 2.50 | 0.286 | +0.103136 | 0.716 | |
| Hkeskin | 0.5792 | +0.184278 | +0.576639 | 3.13 | 0.271 | +0.156036 | **0.847** | ✓ |
| Hkeskin | 0.6189 | +0.222319 | +0.840207 | 3.78 | 0.259 | +0.217269 | **0.977** | ✓ |
| Hkeskin | 0.6581 | +0.238466 | +1.107388 | 4.64 | 0.244 | +0.270115 | **1.133** | ✓ |
| Hkeskin | 0.6983 | +0.242369 | +1.365068 | 5.63 | 0.230 | +0.314299 | 1.297 | |
| Hkeskin | 0.7386 | +0.169534 | +1.469730 | 8.67 | 0.218 | +0.321109 | 1.894 | |
| Hkeskin | 0.7774 | +0.011062 | +1.510953 | 137 | 0.199 | +0.300156 | 27.1 | |
| **HA4** | 0.5392 | +0.126393 | +0.381631 | 3.02 | 0.265 | +0.101313 | **0.802** | ✓ |
| HA4 | 0.5788 | +0.169768 | +0.614420 | 3.62 | 0.243 | +0.149391 | **0.880** | ✓ |
| HA4 | 0.6181 | +0.219593 | +0.939921 | 4.28 | 0.231 | +0.217190 | **0.989** | ✓ |
| HA4 | 0.6566 | +0.251604 | +1.379333 | 5.48 | 0.210 | +0.290182 | **1.153** | ✓ |
| HA4 | 0.6962 | +0.267287 | +2.065514 | 7.73 | 0.189 | +0.390535 | 1.461 | |
| HA4 | 0.7349 | +0.154810 | +3.211518 | 20.7 | 0.146 | +0.469972 | 3.036 | |
| HA4 ‡ | 0.9217 | +12.1348 | −62.28 | −5.13 | 0.038 | −2.336 | −0.193 | ‡ |

‡ = 160'ın sağlık kuralını geçmeyen bant (HA4'ün erfc penceresi τ > 0.68'de
merdiveni sıfırlıyor; `τ_eff = 0.92`, `|Γ|` patlak). **Hiçbir hükme
girmedi.**

**Sayım:** `son` **5/7**, `Hkeskin` **3/7**, `HA4` **4/6** ⇒ **12/20**.
τ_eff ∈ 0.54–0.70 penceresinde (üç gazda beşer bant) **12/15**
(kaçanlar: Hkeskin 0.539 → 0.716 ve 0.698 → 1.297; HA4 0.696 → 1.461).

### 3b. NOMİNAL MERDİVENLE (P0 — tamamen sıfır parametre)

| gaz | 0.539 | 0.579 | 0.619 | 0.658 | 0.698 | 0.739 | 0.777 | ±%25 |
|---|---|---|---|---|---|---|---|---|
| **son** | 0.810 ✓ | 0.998 ✓ | 1.109 ✓ | 1.231 ✓ | 1.306 | 1.772 | 4.729 | **4/7** |
| **Hkeskin** | 0.739 | 0.958 ✓ | 1.117 ✓ | 1.276 | 1.436 | 2.058 | 28.3 | **2/7** |
| **HA4** | 0.783 ✓ | 0.882 ✓ | 0.969 ✓ | 1.087 ✓ | 1.190 ✓ | 2.011 | ‡ | **5/6** |

Toplam **11/20**; τ_eff ∈ 0.54–0.70 penceresinde **11/15**. Nominal ile ölçülü çizgiler arasındaki fark küçüktür —
164'ün sadakat onarımının doğrudan meyvesi (152'nin gazında P0/P0m farkı
×6 idi, 163 §3a).

### 3c. Nerede tutmuyor (ve neden)

* **τ_eff ≤ 0.50:** öngörü işaret bile veremiyor (`son`: −0.005).
  163 §5(iii) ile aynı: çıplak Ç1'in `sin(2πτ_Q)` kökü τ = ½'te;
  ölçülen orada sıfırdan geçmiyor. **T2 bu kökün neden yanlış yerde
  olduğunu çözüyor (§6).**
* **τ_eff ≥ 0.74:** ölçülen `A²s2` sönerken (0.212 → 0.076 → 0.011)
  öngörü büyümeye devam ediyor. Oran patlıyor; bu bantlar hükme
  sayıldı ama başarısızlık olarak yazılıyor.

---

## 4. KANAL PAYLARI — H-Ç SINAVI ve DÖRDÜNCÜ KANALIN KEŞFİ

### 4a. Ç1 / Ç2 / Ç3+Ç4 (τ_c = 0.95, ölçülü çizgiler)

| gaz | τ_eff | Ç1 | Ç1(T_a) | Ç1(B+C) | Ç2 | Ç3+Ç4 | **Ç1/tot** | **(Ç3+Ç4)/tot** |
|---|---|---|---|---|---|---|---|---|
| son | 0.5392 | +0.076834 | −5.0e−06 | +0.076919 | −6.5e−05 | +0.272810 | **0.220** | 0.780 |
| son | 0.5791 | +0.140978 | −1.1e−05 | +0.141095 | −2.0e−05 | +0.416871 | 0.253 | 0.747 |
| son | 0.6188 | +0.219125 | +1.6e−05 | +0.219216 | −4.5e−05 | +0.605168 | 0.266 | 0.734 |
| son | 0.6578 | +0.308387 | +1.1e−06 | +0.308480 | −6.5e−05 | +0.794764 | 0.280 | 0.720 |
| son | 0.6980 | +0.416387 | +1.4e−05 | +0.416445 | +1.6e−05 | +0.993011 | 0.295 | 0.705 |
| son | 0.7382 | +0.522678 | +2.2e−05 | +0.522703 | −2.3e−04 | +1.085205 | 0.325 | 0.675 |
| son | 0.7766 | +0.689266 | −1.5e−04 | +0.689444 | +1.7e−05 | +1.116212 | **0.382** | 0.618 |
| Hkeskin | 0.5393…0.7774 | | | | | | **0.211 → 0.411** | 0.789 → 0.589 |
| HA4 | 0.5392…0.7349 | | | | | | **0.227 → 0.309** | 0.773 → 0.691 |

> **HÜKÜM (H-Ç).** Ç1 baskın **DEĞİL**: payı **%21–41**. Ç2 **ölü**
> (|Ç2/tot| ≤ 4×10⁻⁴, sağlıklı bantların çoğunda ≤ 6×10⁻⁵; kule terimleri asal-kuvvet teklik lemmasının
> öngördüğü gibi ihmal edilebilir). `T_a` (çizginin kendi terimi) 163'ün
> ispatladığı gibi **özdeş olarak sıfır** — ölçülen karşılığı
> `|Ç1(T_a)| ≤ 1.5×10⁻⁴`, yani Ç1'in tamamı (B)+(C) ailesidir.

### 4b. DÜZGÜN TARAK KONTROLÜ — Ç3 ile Ç4'ü AYIRAN ÖLÇÜ

Aynı çizgi genlikleriyle ama siteler düzgün (`s_n → s_0 + ḡn`): orada
`κ` Dirichlet çekirdeğidir, `Δ = k·dres`'te ÖZDEŞ 0'dır ⇒ **düzgün-tarak
toplamı = Ç1 + Ç2 + Ç3 (yakın-rezonans)**; gerçek tarakla farkı
**Ç4 = tarak-rezonansı** kanalıdır.

| gaz | τ_eff | A²s2 ölç | **GERÇEK tarak** | **DÜZGÜN tarak** | **Ç4 = G−D** | **Ç4/G** | u0(G) | u0(D) |
|---|---|---|---|---|---|---|---|---|
| son | 0.5392 | +0.136666 | +0.349579 | +0.080074 | +0.269505 | **0.771** | 1.642 | **1.000** |
| son | 0.5791 | +0.178412 | +0.557829 | +0.156985 | +0.400844 | 0.719 | 1.713 | 1.004 |
| son | 0.6188 | +0.229560 | +0.824247 | +0.222455 | +0.601792 | 0.730 | 1.792 | 1.013 |
| son | 0.6578 | +0.255157 | +1.103086 | +0.328608 | +0.774478 | 0.702 | 1.889 | 1.023 |
| son | 0.6980 | +0.280278 | +1.409414 | +0.421013 | +0.988401 | 0.701 | 2.041 | 1.047 |
| Hkeskin | 0.5393 | +0.144047 | +0.360286 | +0.081241 | +0.279046 | **0.775** | 1.666 | **1.000** |
| Hkeskin | 0.6189 | +0.222319 | +0.840207 | +0.230797 | +0.609409 | 0.725 | 1.810 | 1.012 |
| Hkeskin | 0.6983 | +0.242369 | +1.365068 | +0.424731 | +0.940337 | 0.689 | 2.034 | 1.043 |
| HA4 | 0.5392 | +0.126393 | +0.381631 | +0.091138 | +0.290493 | 0.761 | 1.614 | 1.001 |
| HA4 | 0.6181 | +0.219593 | +0.939921 | +0.251008 | +0.688913 | 0.733 | 1.784 | 1.014 |
| HA4 | 0.6962 | +0.267287 | +2.065514 | +0.571772 | +1.493741 | 0.723 | 2.257 | 1.092 |

> **Ç4/G = 0.66–0.78, ÜÇ GAZIN YEDİ BANDINDA DA.** Ç1 (%21–41) ile
> birlikte alındığında geriye kalan **Ç3 payı yalnız %0.2–1.4**'tür
> (ör. Hkeskin τ_eff = 0.619: `(0.2308 − 0.2288)/0.8402 = 0.0024`).
>
> `u0` sütunu bağımsız bir teşhistir: DÜZGÜN tarakta model alanının
> kendi normalizasyonu **1.000–1.09** (çizgi gösterimi orada TAM),
> gerçek tarakta **1.61–2.26**. Yani §7'nin "alan şişmesi" bir sızıntı
> kusuru değil, **tarağın kendi yapı çarpanıdır** — 163 §6d'nin
> ×10–×40'lık çöküşünün 165'teki tam ölçüsü budur.

### 4c. ALT-MERDİVENDE TERİM SAYIMI — |Δ| dağılımı (bağımsız sınav)

`165_kanal.py` alt-merdivende (τ ≤ 0.60, 245 çizgi) terimleri **tek tek**
sayar; `κ(Δ) = e^{iΔs̄}·κ̃(Δ)` ayrıştırmasıyla merkez fazı TAM, `κ̃`
`dres/20` aralıklı tablodan alınır.

**Hkeskin, alt-merdiven toplamı (τ_s = 0.60), `s2` biriminde**

| Q (τ_Q) | alt-merd. TOPLAM | Δ ≡ 0 (n₀) | \|Δ\|≤1·dres | ≤5 | ≤20 | ≤100 | **≤200·dres** | **oran** |
|---|---|---|---|---|---|---|---|---|
| 521 (0.520) | +0.01397 | +0.00197 (1467) | +0.00194 | +0.00233 | +0.00221 | +0.00245 | +0.00248 | **0.18** |
| 853 (0.561) | +0.01933 | +0.00429 (1467) | +0.00426 | +0.00391 | +0.00400 | +0.00407 | +0.00402 | **0.21** |
| 1367 (0.600)† | +0.01483 | 0 (0) | +0.00047 | −0.00029 | −0.00039 | −0.00042 | −0.00040 | −0.03 |
| 2207 (0.640)† | +0.01701 | 0 (0) | +0.00099 | +0.00057 | +0.00116 | +0.00091 | +0.00087 | 0.05 |

† taşıyıcı çizgi alt-merdivenin DIŞINDA (τ_Q > τ_s) ⇒ Ç1 sıfırdır ve
tablo yalnız Ç3'ü gösterir.

> **HÜKÜM (Ç3'ün doğrudan sayımı).** `|Δ| ≤ 200·dres` olan BÜTÜN terimler
> (ki sayıları yüz binlerce) toplamda Δ ≡ 0'ın üstüne **hiçbir şey
> eklemiyor**: 0.00197 → 0.00248 (±%3 salınım, işareti bile değişiyor).
> Toplamın kalan **%80'i** `|Δ| > 200·dres` bölgesindedir.
> **KALEM'in "Ç3 — HL yakın-rezonansları, Ψ⋆κ ile tartılmış" katmanı
> ölçüldü ve NEREDEYSE SIFIR çıktı.**

### 4d. κ ÖLÇÜLDÜ: neden Ç4 var, Ç3 yok

| Δ/dres | \|κ̃(Δ)\| | Dirichlet (düzgün tarak) |
|---|---|---|
| 0.00 | 1.00000 | 1.00000 |
| 0.50 | 0.63664 | 0.63663 |
| 1.00 | 0.00197 | **0** (özdeş) |
| 2.00 | 0.00098 | 0 |
| 5.00 | 0.00039 | 0 |
| 10.0 | 0.00020 | 0 |

**Tarak rezonans tepeleri** `κ(ω_q)` ↔ 143'ün `𝒢 = −πτ_q a_q cos(πτ_q)`'si:

| q | 2 | 3 | 4 | 5 | 7 | 8 | 9 | 11 |
|---|---|---|---|---|---|---|---|---|
| `Re κ(ω_q)` | −0.03996 | −0.05046 | −0.02648 | −0.05416 | −0.05299 | −0.01639 | −0.02399 | −0.04809 |
| `𝒢` | −0.04008 | −0.05057 | −0.02694 | −0.05462 | −0.05341 | −0.01744 | −0.02557 | −0.04870 |
| fark | %0.3 | %0.2 | %1.7 | %0.8 | %0.8 | %6.0 | %6.2 | %1.3 |
| `Im κ` | ≤ 3e−05 | | | | | | | |

> Resim tamamlandı: `κ`'nın **`Δ ≈ 0`** civarındaki zarfı Dirichlet'tir ve
> `Δ ≥ dres`'te `2×10⁻³` mertebesine düşer (Ç3'ün küçüklüğü); buna
> karşılık `Δ ≈ ±ω_r`'de `κ` **`|𝒢| ≈ 0.02–0.055`**'e yükselir — yani
> yakın-rezonanslardan **yirmi–otuz kat** güçlüdür ve terim sayısı da
> kat kat fazladır. **Üçüncü momenti besleyen şey, dört merdiven
> frekansının ARTIĞINI beşinci bir merdiven çizgisiyle kapatan
> tarak-rezonansıdır.** 163 §6d'nin "taşıyıcı taraktır" hükmü, burada
> `κ`'nın kendi tayfıyla nicelendi.

---

## 5. T3 — b'NİN İŞARETİ: HANGİ KANAL, HANGİ τ_q

`165_T3.py` X̃'yi τ-dilimlerine ayırıp `A²⟨σX^{(j)}X^{(k)}⟩/⟨ρ⟩`
matrisini bant-toplu çıkarır — **iki kez**: (M) ölçülen η ile (163 §6b'nin
aynısı; dilimler + artık ⇒ kapanış `Σ_{jk} = A²s2` ÖZDEŞ, ölçülen fark
≤ 5e−06), (P) model alanı `E_mod` ile (yani dörtlü toplamın kendisi).

### 5a. 163'ün A4 ANOMALİSİ YOK OLDU

**Bant τ_eff ≈ 0.619, ölçülen (M) matrisi** (dilimler
0–0.2 / 0.2–0.4 / 0.4–0.55 / 0.55–0.7 / 0.7–0.85 / 0.85–0.95 / artık):

| hücre | **son** | **Hkeskin** | **HA4** |
|---|---|---|---|
| (0–0.2)² | +0.0247 | +0.0195 | +0.0256 |
| (0–0.2)×(0.2–0.4) | **+0.0612** | **+0.0557** | **+0.0567** |
| (0.2–0.4)² | +0.0281 | +0.0280 | +0.0200 |
| (0.2–0.4)×(0.4–0.55) | +0.0406 | +0.0417 | +0.0474 |
| (0.4–0.55)×(0.55–0.7) | +0.0678 | +0.0687 | +0.0920 |
| (0.55–0.7)² | +0.0641 | +0.0655 | +0.0771 |
| (0.55–0.7)×artık | −0.1527 | −0.1479 | −0.1689 |
| artık² | +0.0311 | +0.0253 | +0.0384 |
| **kapanış** | **+0.229559** | **+0.222319** | **+0.219592** |

> **HÜKÜM.** Üç gazın matrisleri **hücre hücre aynı işarette ve ±0.01
> içinde**. 163 §6b'nin "aynı hücrede son **+0.0216**, A4 **−0.0554**"
> karşıtlığı, **152'nin kusurlu inşasının artefaktıymış**; 164 kök
> seçimini ve `n−½` konvansiyonunu onarınca yok oldu. İşaret sayımı
> (|hücre| > 0.004, yedi bant): son **+120/−42**, Hkeskin **+123/−46**,
> HA4 **+118/−49** — üçü de aynı.
> **"b'nin işaret farkı" artık AÇIKLANACAK BİR OLGU DEĞİLDİR.**

### 5b. İŞARET HÂLÂ ALÇAK-τ ÇİFT BLOĞUNDA DOĞUYOR

`(0–0.2)² + (0–0.2)×(0.2–0.4) + (0.2–0.4)²` bloğunun toplam içindeki payı:

| τ_eff | 0.539 | 0.579 | 0.619 | 0.658 | 0.698 | 0.739 |
|---|---|---|---|---|---|---|
| **son** | 0.424 | 0.459 | 0.496 | 0.507 | 0.492 | 0.509 |
| **Hkeskin** | 0.402 | 0.427 | 0.464 | 0.488 | 0.474 | 0.451 |
| **HA4** | 0.395 | 0.418 | 0.466 | 0.524 | 0.509 | 0.626 |

**Toplamın %40–51'i, X̃'nin `q ≲ 123` (τ ≤ 0.4) çizgilerinin çift
bloğundan geliyor** — 163 §6b'nin lokalizasyonu ayakta; ama artık
GAZ-BAĞIMSIZ. Belirlenimli (P) matrisi de aynı hücreleri aynı işaretle
üretiyor (ör. τ_eff = 0.619: (0–0.2)×(0.2–0.4) = +0.0908 / +0.0917 /
+0.0915), yalnız toplu kazancı ~3.7 kat büyük (§4b).

**Kanal cevabı:** işaret Ç1'de de, Ç3+Ç4'te de aynıdır (her ikisi de
pozitif ve aynı bantta aynı yönde); yani `b`'nin işaretini **tek bir
kanal değil, alçak-τ çift bloğunun ortak kinematiği** taşıyor.

---

## 6. T2 — F(τ_q; τ_band): ÇİZGİ-BAŞINA KATKI YOĞUNLUĞU

### 6a. Ç1'in yoğunluğu KAPALI FORMDA ve BANT-BAĞIMSIZ

§2c'den `J_2^{Ç1(B+C)} = (y_Q/2)·Σ_q Re[h'_q ȳ_q]` ⇒

```
F(τ_q) = Re[ h'_q ȳ_q ]        (bant yalnız ortak taşıyıcı çarpanıyla girer)
```

**Yani Ç1'in F'i `τ_band`'dan BAĞIMSIZDIR** — sınavın sorduğu bant
bağımlılığı Ç1'de YOKTUR ve bu bir öngörü olarak kayda geçer.

### 6b. Empirik eğri ve yarı-analitik form (Hkeskin)

`Σ F_emp = 0.053564` ↔ `Σ F_çıplak = Σ a²sin²(2πτ) = 0.036784`
(**oran 1.456**). τ-kutularında:

| τ-bant | 0.40–45 | 0.45–50 | **0.50–55** | 0.55–60 | 0.60–65 | 0.65–70 | 0.70–75 | 0.75–80 | 0.85–90 |
|---|---|---|---|---|---|---|---|---|---|
| ΣF_emp | +0.00668 | +0.00647 | **+0.00654** | +0.00696 | +0.00729 | +0.00643 | +0.00541 | +0.00384 | +0.00121 |
| ΣF_çıplak | +0.00244 | +0.00039 | **+0.00030** | +0.00180 | +0.00407 | +0.00584 | +0.00671 | +0.00632 | +0.00291 |
| oran | 2.74 | **16.8** | **21.7** | 3.88 | 1.79 | 1.10 | 0.81 | 0.61 | 0.42 |

> **Kritik gözlem: `F_emp` τ = ½'te SIFIRLANMIYOR**, çıplak
> `a²sin²(2πτ)` ise sıfırlanıyor. 163'ün `τ_Q = ½`'teki sahte kökünün
> kaynağı budur.

**Yarı-analitik form (yarım-gap ailesi, a²-ağırlıklı en küçük kareler):**

| gaz | form | katsayı(lar) | bağıl artık |
|---|---|---|---|
| Hkeskin | `F/a² = c₁ sin²(2πτ)` (163'ün çıplak formu, c₁ ≡ 1) | c₁ = **0.910** | **0.443** |
| Hkeskin | `F/a² = c₁ sin²(2πτ) + c₂ sin⁴(πτ)` | c₁ = **0.546**, c₂ = **0.666** | **0.039** |
| **son** | `F/a² = c₁ sin²(2πτ)` | c₁ = **0.710** | **0.511** |
| **son** | `F/a² = c₁ sin²(2πτ) + c₂ sin⁴(πτ)` | c₁ = **0.384**, c₂ = **0.596** | **0.045** |

(gerçek gazda `Σ F_emp = 0.043905` ↔ `Σ F_çıplak = 0.036784`, oran **1.194**;
iki gazda da `c₂/c₁ ≈ 1.2–1.6` ve iki-terimli form artığı **%4–4.5**'e
indiriyor.)

Eşdeğer biçim: `F/a² = sin²(πτ)·[4c₁cos²(πτ) + c₂ sin²(πτ)]` — yani
**aynı yarım-gap ailesinin İKİ üyesi**, iki katsayı, artık %3.9.

**Türetimle destek.** Gerçek çizgi genliği sabit değil, yarım-gap
modülasyonludur: `2a_q sin(ω_q g_n/2) = 2a_q sin(πτ_q(1+ds_n))`.
Açılım
`sin(πτ(1+ds)) = sin(πτ)cos(πτ ds) + cos(πτ)sin(πτ ds)`
ve `X̃`'nin bond-ortalama çarpanı
`cos(πτ(1+X̃)) = cos(πτ)cos(πτX̃) − sin(πτ)sin(πτX̃)`
ikinci terimleri getirir; `Re[h'ȳ]`'de birinci terimlerin çarpımı
`4a²sin²πτ cos²πτ = a²sin²(2πτ)` (çıplak), İKİNCİ terimlerin çarpımı ise
`4a² sin⁴(πτ)·⟨…⟩`'dır. **`sin⁴(πτ)` üyesi bu yüzden ailenin doğal
ikinci elemanıdır ve katsayısı Debye–Waller tipi bir ortalamadır.**
`πτ σ_ds = π·0.5·0.431 = 0.68 rad` olduğundan bu terim τ ≈ ½'te
ihmal edilemez — ölçülen `c₂/c₁ = 1.22` ile tutarlı.

### 6c. TOPLAMIN yoğunluğu — bant bağımlılığı BURADA

`J_2 = Σ_q ⟨E X^{(q)} X e^{−iω_Q s}⟩` ile çizgi-başına yoğunluk
(Hkeskin, `s2` biriminde, τ_q kutuları 0.1 genişlikte):

| Q (τ_Q) | 0.0–0.1 | 0.1–0.2 | 0.2–0.3 | 0.3–0.4 | 0.4–0.5 | 0.5–0.6 | **0.6–0.7** | 0.7–0.8 | Σ |
|---|---|---|---|---|---|---|---|---|---|
| 521 (0.520) | +0.0014 | +0.0043 | +0.0054 | +0.0071 | +0.0028 | +0.0051 | −0.0004 | +0.0000 | +0.0262 |
| 1367 (0.600) | +0.0019 | +0.0069 | +0.0076 | +0.0074 | +0.0056 | +0.0036 | **+0.0217** | +0.0024 | +0.0582 |
| 3571 (0.680) | +0.0042 | +0.0063 | +0.0090 | +0.0109 | +0.0052 | +0.0033 | **+0.0246** | −0.0025 | +0.0624 |

> Toplamın yoğunluğu iki parçalıdır: (i) τ_q ≲ 0.5'te bant-bağımsız,
> Ç1'in F'ini izleyen geniş bir taban; (ii) **τ_q ≈ τ_Q'da keskin bir
> tepe** (Q=3571 için 0.6–0.7 kutusu tek başına toplamın %39'u).
> **Bant bağımlılığı Ç3+Ç4'ün özelliğidir; Ç1'in değil.** `son`'da
> aynı desen, aynı sayılarla (0.0259 / 0.0522 / 0.0700).

---

## 7. ÇİZGİ GÖSTERİMİNİN GEÇERLİLİK SINIRI (yeni sistematik)

Sabit yarım-gap genlikli çizgi modeli yüksek τ'da alanı **şişiriyor**:

| gaz | kaynak | τ_c | çizgi | `Var(E_mod)/Var(η)` | `Var(X_mod)/Var(X̃)` | `g_E` | `g_X` | ⟨nrm/ölç⟩ (0.54–0.70) |
|---|---|---|---|---|---|---|---|---|
| son | ölçülen | 0.60 | 245 | 1.377 | 0.974 | 0.721 | 0.921 | +1.482 |
| son | ölçülen | 0.70 | 656 | 1.838 | 1.198 | 0.651 | 0.819 | +1.050 |
| **son** | **ölçülen** | **0.95** | **8981** | **2.139** | **1.585** | 0.598 | 0.705 | **+0.996** |
| son | nominal | 0.95 | 8981 | 3.733 | 1.553 | 0.429 | 0.738 | +1.091 |
| Hkeskin | ölçülen | 0.95 | 8981 | 2.142 | 1.562 | 0.584 | 0.704 | +0.994 |
| HA4 | ölçülen | 0.95 | 8981 | 2.088 | 1.520 | 0.622 | 0.717 | +1.057 |
| HA4 | nominal | 0.95 | 8981 | 2.101 | 1.023 | 0.618 | 0.936 | +0.982 |

İki **bağımsız, sıfır-parametreli** ölçüt aynı sınırı veriyor:

```
faz seğirmesi:   τ_c = 1/(2π σ_Ĉ)  = 0.573 (Hkeskin) / 0.578 (HA4) / 0.583 (son)
yarım-gap mod.:  τ_c = 1/(π σ_ds)  = 0.738 / 0.733 / 0.778
```

(`σ_Ĉ` = birikmiş adım sapmasının TRENDSİZ standart sapması = 0.273–0.278;
163 §6d'nin `σ_Ĉ ≈ 0.27`'siyle aynı sayı. `Ĉ_n = (m_n−m_0)L/2π − n` sabit
`L` yüzünden pencere boyunca yavaş bir kayma taşır ve ondan arındırılmalıdır.)

Üçüncü, tamamen deneysel ölçüt aynı yeri gösteriyor: `Var(E_mod)`,
ölçülen `Var(η)`'yı **τ_c ≈ 0.57**'de tutuyor.

> **Ama şişme bir MODEL kusuru değil, TARAK etkisidir:** §4b'nin düzgün
> tarağında aynı model `u0 = 1.000–1.09` veriyor. Yani çizgi gösterimi
> düzgün tarakta TAM; gerçek tarakta ise çizgiler tarak-rezonansı
> üzerinden birbirine bağlanıyor. **Kalibrasyonun (n_u2 ≈ 0.15–0.31)
> fiziksel karşılığı tam olarak Ç4'ün ortak kazancıdır** — ve o kazanç
> ρ-kanalıyla σ-kanalında AYNI olduğu için ikinci momentten okunup
> üçüncüye taşınabiliyor (§3'ün başarısının sebebi budur).

---

## 8. T4 — GERÇEKTE: Ç3+Ç4 PAYI ve HL'nin ÜST SINIRI

Sadakatli gazlar **saf asal merdivendir**: içlerinde HL korelasyonu
YOKTUR (yalnız `N̄+S = n−½`). Gerçek gaz aynı hesapla koşulduğunda
kanal payı farkı, aritmetik-incelik katmanının payına **doğrudan üst
sınır** verir.

| τ_eff | **(Ç3+Ç4)/tot `son`** | `Hkeskin` | `HA4` | \|fark\| |
|---|---|---|---|---|
| 0.5392 | **+0.7804** | +0.7890 | +0.7732 | 0.0086 / 0.0072 |
| 0.5791 | **+0.7473** | +0.7460 | +0.7488 | 0.0013 / 0.0015 |
| 0.6188 | **+0.7342** | +0.7277 | +0.7402 | 0.0065 / 0.0060 |
| 0.6578 | **+0.7205** | +0.7118 | +0.7334 | 0.0087 / 0.0129 |
| 0.6980 | **+0.7046** | +0.6911 | +0.7189 | 0.0134 / 0.0143 |
| 0.7382 | **+0.6750** | +0.6546 | +0.6911 | 0.0204 / 0.0160 |
| 0.7766 | +0.6182 | +0.5886 | ‡ | 0.0296 / — |

> **HÜKÜM (T4).** Gerçek gazın kanal payı **altı sağlıklı bandın
> altısında da iki saf-merdiven gazının ARASINDA** (`HA4 > son >
> Hkeskin` ya da tersi). Fark τ_eff ≤ 0.70'te **≤ 0.014**, yani
> `⟨σX̃²⟩`'nin **%1.4'ünden azdır**; τ_eff ≤ 0.62'de **≤ %0.9**.
>
> **"Gaz Hardy–Littlewood'u duyuyor mu?" sorusuna ilk nicel cevap:
> BU ÖLÇÜDE DUYMUYOR. Ç3/HL payının ölçülen üst sınırı %1.4'tür**
> (ve Ç3'ün kendisi zaten %0.2–1.4 — §4c). Bu bir SINIRDIR, bir
> keşif değil; gerçek gaz ile saf asal merdiven gazları bu gözlenebilirde
> ayırt edilemiyor.
>
> Ayrıca (§3a) `son`'un ±%25 sayımı (5/7) `Hkeskin`'inkinden (3/7) daha
> iyidir — belirlenimli dörtlü toplam **gerçek gazda sentetiklerden daha
> iyi çalışıyor.**

---

## 9. DENETİM

| | sınav | sonuç |
|---|---|---|
| **V1** | 165'in bant birleştirmesi = 160'ın kayıtlı JSON'u mu? (`son`, t1 ızgarası, taban 0.40) | 9 bant; `A²s2`, `A·s1`, `A³s3`, `A·κ1`, `τ_eff`, `kul` **maks fark 0.0e+00**; `s0` **2.6e−16** |
| **V2** | Ç1'in ÇIPLAK limiti = 163 §2e'nin kapalı formu mu? | yedi τ_Q'da maks fark **4.8e−06** (bağıl 3e−04) |
| **V3** | s-sitesi tayfı 164 §3g'nin çizgi oranlarını veriyor mu? | `arg y_q + πτ_q` = −0.021 … +0.005 (τ ≤ 0.2; çıplak öngörü 0), `|y_q|/B_q` = 0.995…0.968 ↔ 164'ün `|c_q|/b_q` = 1.005…1.052; **fark %1–8** ve τ ile büyüyor (yarım-gap modülasyonunun Debye–Waller sönümü, §6b) |
| **V4** | SAYIM = SENTEZ özdeşliği (dörtlü toplamın kapanışı) | Hkeskin, Q=521, alt-merdiven τ ≤ 0.44 (138 işaretli çizgi, 19 044 çift), 3000 site: sentez `9.938680063721e−05 − 7.155902980906e−05j`, sayım `9.938680067489e−05 − 7.155903010417e−05j`, **bağıl fark 2.4e−09**; `son` için 1000 sitede **3.5e−10** |
| **V5** | τ-dilim matrisinin kapanışı (T3) | `Σ_{jk} M_{jk} − A²s2` yedi bantta **≤ 5.0e−06** (bağıl ≤ 5e−05) |
| **V6** | κ'nın Dirichlet limiti ve 𝒢 tepeleri | `κ̃(½·dres) = 0.63664` ↔ 2/π = 0.63662; `κ(ω_q)` ↔ `𝒢` sekiz çizgide **%0.3–6.2**, `Im κ ≤ 3e−05` |
| **V7** | düzgün tarakta normalizasyon | `u0(düzgün)` = **1.000 – 1.09** (τ_eff ≤ 0.70), gerçek tarakta 1.61–2.26 |

---

# HÜKÜM

## (i) HESAP KURULDU ve BUDAMASIZ KAPANDI

> Tek-site indirgemesi 163'ün birinci-mertebe kesmesini **tamamen
> ortadan kaldırdı**: `κ(Δ)`'nın tanımı sayesinde 5.8×10¹² kombinasyonun
> κ-ağırlıklı toplamı, model alanlarının gerçek tarak üzerindeki
> ortalamasına ÖZDEŞ. Hiçbir kesme, hiçbir budama yapılmadı; özdeşlik
> doğrudan terim sayımıyla **2.4e−09**'da doğrulandı (V4).

## (ii) T1 HEDEFİ (±%25) BÜYÜK ÖLÇÜDE TUTTU

> 163: **4/27** bant. 165: **12/20** (ölçülü çizgiler) ve **11/20**
> (nominal); τ_eff ∈ 0.54–0.70 penceresinde **12/15**. Tek kalibrasyon
> ikinci-moment normalizasyonudur (`u2`), üçüncü momentten girdi yoktur
> ve **kanal paylarını değiştirmez**. Düşük τ (≤0.50) ve yüksek τ
> (≥0.74) hâlâ tutmuyor ve bu böyle yazılıyor.

## (iii) H-Ç HİYERARŞİSİ DÜŞTÜ — ve yerine ÖLÇÜLMÜŞ bir hiyerarşi geldi

> | kanal | KALEM'in beklentisi | **165'in ölçtüğü** |
> |---|---|---|
> | **Ç1** çift-iptalli koherent | **baş şüpheli / baskın** | **%21–41** — baskın DEĞİL |
> | **Ç2** tam çarpımsal üçlüler | ikincil | **≤ 4×10⁻⁴** — ÖLÜ |
> | **Ç3** HL yakın-rezonans (Ψ⋆κ) | "asıl ev" | **%0.2–1.4** — NEREDEYSE SIFIR |
> | **Ç4** TARAK-REZONANSI (`Δ ≈ ±ω_r`) | *hiyerarşide yok* | **%66–78 — BASKIN** |
>
> Ç4'ün varlığı iki bağımsız ölçümle sabittir: düzgün-tarak kontrolünde
> toplamın %66–78'i kayboluyor (§4b), ve alt-merdiven terim sayımında
> `|Δ| ≤ 200·dres` bölgesi toplamın yalnız %18–21'ini (Ç1 dahil)
> topluyor (§4c). Mekanizması ölçüldü: `κ(±ω_r) = 𝒢 = −πτ a cos(πτ)`
> (%0.3–6), yani **dört merdiven frekansının artığını KAPATAN şey
> beşinci bir merdiven çizgisi, ama X̃'de değil TARAĞIN KENDİSİNDE.**

## (iv) 163'ÜN A4 İŞARET ANOMALİSİ BİR İNŞA ARTEFAKTIYMIŞ

> Sadakatli gazlarda `Hkeskin`, `HA4` ve `son`'un τ-dilim matrisleri
> hücre hücre aynı işarette ve ±0.01 içinde. 163 §6b'nin "gerçek(+) ↔
> A4(−)" karşıtlığı 152'nin kök-seçimi + yarım-seviye kusurundan
> geliyormuş. İşaret hâlâ **alçak-τ çift bloğunda** (toplamın %40–51'i)
> ama artık gaz-bağımsızdır.

## (v) T2: YARI-ANALİTİK FORM İKİ KATSAYIYLA OTURDU

> `F(τ_q) = Re[h'_q ȳ_q]`, Ç1'in bant-BAĞIMSIZ yoğunluğu.
> `F/a² = 0.546 sin²(2πτ) + 0.666 sin⁴(πτ)` (Hkeskin) artığı **%3.9**'a,
> `0.384 sin²(2πτ) + 0.596 sin⁴(πτ)` (gerçek gaz) **%4.5**'e indiriyor;
> çıplak tek-terimli form (163'ün `sin²(2πτ)`'si) %44–51 artık bırakıyor.
> İkinci üye yarım-gap açılımının ikinci teriminden geliyor ve
> **163'ün `τ_Q = ½`'teki sahte sıfırını kaldırıyor.**

## (vi) T4: HL PAYININ ÜST SINIRI %1.4

> Gerçek gazın kanal payı iki saf-asal-merdiven gazının arasında;
> fark τ_eff ≤ 0.70'te ≤ 0.014. **Gaz, bu gözlenebilirde HL'yi
> duymuyor.** Kayda geçen sonuç bir üst sınırdır.

---

## Sıradaki adım (bu ölçümün işaret ettiği)

1. **Ç4'ün KAPALI FORMU aranmalı.** `κ(±ω_r) = 𝒢_r` ölçüldü ve 143'ün
   formülünü verdi. Geriye kalan iş, `Δ = ε₁ω₁−ω_Q+ε₂ω₂+ε₃ω₃ = ±ω_r`
   BEŞLİ rezonans koşulunun asal-kuvvet çözümlerini saymak ve
   `Σ 𝒢_r · h'y y` toplamını kapalı yazmaktır. Asal-kuvvet teklik
   lemması burada da çalışmalı (ve Ç2'nin ölü olması bunun ön habercisidir);
   eğer çözerse **s2'nin %70'i kaleme geçer.**
2. **Şişme çarpanının (`1/n_u2 ≈ 3.4`) kapalı formu.** Düzgün ↔ gerçek
   tarak oranı üç gazda ve yedi bantta 0.66–0.78 aralığında ÇOK
   kararlı; muhtemelen `⟨e^{iωḡĈ}⟩` karakteristik fonksiyonunun bir
   fonksiyoneli. Ölçülmesi ucuz (bir O(N) geçiş).
3. **τ_eff ≥ 0.74'te ölçülen sönerken öngörü büyüyor.** Bu, dörtlü
   toplamın kapsamadığı bir sönüm mekanizmasıdır (muhtemelen `|Γ|`'nın
   kendi sönümü); 160'ın `E1`'iyle bağı kurulmalı.
4. **`orta`, `dusuk`, `P1` ve 164'ün açık borcu olan yedi gaz bu koşuda
   KOŞULMADI.** τ_c taraması yalnız 0.60/0.70/0.95; taban ekseni yalnız
   0.40; ızgara yalnız `t1`.
5. **Ç3'ün ölü çıkması, 163'ün T2 negatif sonucunu KAPATIYOR.** "HL
   katmanının büyüklüğü kutu pencereyle ölçülemez" (163 §4) doğruydu;
   doğru pencereyle ölçülünce **katman zaten yok**. Bu, aramanın
   yerini değiştirir: HL varsa `Ĉ`'nin/`κ`'nın kendi tayfında olmalıdır,
   `Ψ`'de değil.

---

## Dürüstlük notları

* **Uydurma yok; her sayı bu koşuda üretildi.** Tablolar
  `scratchpad/165/`'teki `K_<gaz>_t0.4.json` (3), `TAR_<gaz>_*.json` (3),
  `T3_<gaz>_*.json` (3), `T2_<gaz>_*.json` (2), `KAN_<gaz>_*.json` (2),
  `analiz_cikti.txt`, `dogrulama_cikti.txt`, `log_*.txt` çıktılarından
  **otomatik** alındı. 160/163/164'ten alıntılanan tek şey referans
  satırlarıdır ve her biri V1–V3'te bu koşunun kendi ölçümüyle
  karşılaştırılmıştır.
* **±%25 hedefi TEK BİR KALİBRASYONLA tutuyor ve bu saklanmıyor.**
  HAM (kalibrasyonsuz) oran 2.5–8.7'dir. Kalibrasyon ikinci-moment
  normalizasyonudur (`n_u2 = A²u2_ölç/A²u2_öng`, 0.15–0.31): üçüncü
  momentten hiçbir girdi taşımaz, ama **bir serbestlik derecesidir** ve
  hükümde "sıfır parametre" DEĞİL, "ikinci momentlerle sabitlenmiş"
  diye anılmalıdır. Alan-regresyonu kalibrasyonu (`g = g_E g_X²`)
  bağımsız olarak %1–20 içinde aynı sayıyı veriyor (ör. Hkeskin
  τ_eff=0.539: 0.2895 ↔ 0.286).
* **Çizgi gösteriminin şişmesi bir SİSTEMATİKTİR ve ölçülmüştür**
  (§7): `Var(E_mod)/Var(η)` = 2.1–3.7. Bu koşu onu düzeltmedi, ölçtü
  ve nedenini (tarak-rezonansı) yalıttı. τ_c taraması gösteriyor ki
  `τ_c = 0.60`'ta ham oran 1'e yaklaşıyor (0.77–1.36) ama o zaman
  τ_band > τ_c olan bantlarda taşıyıcı çizgi model kümesinin dışında
  kalıyor ve **Ç1 özdeş olarak sıfırlanıyor** — bu yüzden hüküm
  `τ_c = 0.95` (budamasız) sürümünden verildi.
* **Ç3'ün ölçümü İKİ farklı yöntemle yapıldı** (düzgün-tarak kontrolü ve
  alt-merdiven terim sayımı) ve ikisi de %0.2–3 veriyor. Alt-merdiven
  sayımı `τ_s = 0.60` (245 çizgi) ile sınırlıdır; tam merdivende terim
  sayımı yapılmadı, çünkü `|Δ| ≤ 200·dres` kesmesi orada 10⁹ terim
  üretiyor. **Budama gerekçesi:** `|κ̃(Δ)|` `Δ ≥ dres`'te 2e−03'e
  düşüyor ve `1/Δ` ile sönüyor (§4d); kesilen bölgenin katkısı
  düzgün-tarak kontrolüyle bağımsız olarak sınırlanmıştır (o kontrol
  BÜTÜN merdiveni kapsar).
* **HA4'ün `τ_eff = 0.92` bandı (nominal ızgarada (0.76,0.80))** 160'ın
  sağlık kuralını geçmiyor: erfc penceresi τ > 0.68'de merdiveni
  sıfırladığından bant güç farkı işaret değiştiriyor (`A²s2 = +12.13`,
  `u0p = −57`). `‡` ile işaretli ve **hiçbir hükme girmedi.**
* **`son`'un "nominal" sürümü penceresiz merdiveni kullanır** — gerçek
  gaz için doğru olan budur, ama gerçek gazın merdiveni τ ≤ 0.95'te
  kesilmiş sayılmıştır (kesim etkisi τ_c taramasında görünüyor).
* **V3'te 164'ün çizgi oranlarıyla %1–8 fark var ve gizlenmiyor.**
  164 `ds`'nin `m_n` sitesindeki çizgisini ölçtü, 165 `X̃0`'ın `s_n`
  sitesindeki çizgisini. Aradaki fark yarım-gap modülasyonunun
  Debye–Waller sönümüdür ve τ ile büyüyor — **§6b'nin `sin⁴(πτ)`
  teriminin bağımsız bir işareti.**
* **`arg y_q + πτ_q` düşük τ'da 0.02 rad, ama τ = 0.38'de 0.32 rad,
  τ = 0.58'de 2.36 rad.** Yani "belirlenimli merdiven, faz `πτ`"
  varsayımı `τ ≲ 0.3`'te doğrudur; ötesinde bozulur. 163 §1c'nin
  "bütün çizgilerin fazı ≈ 0" kaydı `ds`'nin düşük-τ çizgileri
  içindi ve orada hâlâ geçerlidir.
* **Ç1'in `T_a` (kendi-çizgi) ailesi ölçüldü ve özdeş sıfır çıktı**
  (`|A²s2_Ta| ≤ 1.5e−04`): 163'ün "köşegen tam ölür" hükmü 165'in
  bağımsız cebrinde de ayakta.
* **Figür üretilmedi;** bu koşuda görselleştirme istenmedi.
* **git'e dokunulmadı.**
* **Süreler.** Model tayfı (8981 çizgi × 3e5 site) 66 s/gaz (önbelleğe
  alınıyor); T1 sürücüsü (9 bant, 2276 frekans, 2 kaynak × 3 τ_c)
  **8.5–9.6 dk/gaz**; τ-dilim matrisi (49 alan) 6.0 dk/gaz;
  düzgün-tarak kontrolü 7.9 dk/gaz; T2 9.1 dk/gaz; kanal taraması
  5.7 dk/gaz; denetim 4 dk. Üç gaz eşzamanlı; toplam ~75 dk duvar saati.

---

## Ek — dosyalar ve tekrar-üretim

| dosya | ne |
|---|---|
| `165_configs/165_cekirdek.py` | türetim (docstring'de tam), `Model165` (s-sitesi tayfı, nominal/ölçülen genlikler, koherans ölçütleri, tek-skaler dekonvolüsyon), `kanal1` (Ç1 kapalı toplamı), `kule_cozumleri` (Ç2), `sentez`, `tayf_s`, `kappa`, `c1_capla` (163'ün çıplak formu) |
| `165_configs/165_kos.py` | T1/T3 sürücüsü: bant bant ölçüm + budamasız dörtlü toplam + kanal etiketleri + `u2` normalizasyonu + τ_c taraması |
| `165_configs/165_tarak.py` | DÜZGÜN tarak kontrolü — Ç4'ün payı |
| `165_configs/165_kanal.py` | alt-merdivende terim-terim sayım, `\|Δ\|` kümülatifi, `κ` tanısı (Dirichlet + 𝒢 tepeleri), V4 |
| `165_configs/165_T2.py` | F(τ_q) empirik eğrisi + yarı-analitik fit + toplamın çizgi-başına yoğunluğu |
| `165_configs/165_T3.py` | τ-dilim matrisi (ölçülen M + belirlenimli P) |
| `165_configs/165_dogrulama.py` | V1 (160 bit düzeyi), V2 (163 çıplak form), V3 (164 sadakat), V4 (sayım=sentez) |
| `165_configs/165_analiz.py` | hüküm tabloları |

Ham çıktılar `scratchpad/165/`: `K_<gaz>_t0.4.json` (3),
`TAR_<gaz>_olculen_tc0.95.json` (3), `T3_<gaz>_olculen_tc0.95.json` (3),
`T2_<gaz>_olculen_tc0.95.json` (2), `KAN_<gaz>_ts0.6_olculen.json` (2),
`tayf_<gaz>_t0.4_tm0.95.npz` (3), `analiz_cikti.txt`,
`dogrulama_cikti.txt`, `log_*.txt`.

```
cd /Users/ugursezen/Desktop/arin/deney
for v in Hkeskin HA4 son; do
  .venv/bin/python qm_riemann/165_configs/165_kos.py $v 0.40 olculen,nominal 0.95,0.70,0.60
done
for v in Hkeskin HA4 son; do .venv/bin/python qm_riemann/165_configs/165_T3.py   $v 0.40 olculen 0.95; done
for v in Hkeskin HA4 son; do .venv/bin/python qm_riemann/165_configs/165_tarak.py $v 0.40 olculen 0.95; done
for v in Hkeskin son;     do .venv/bin/python qm_riemann/165_configs/165_T2.py    $v 0.40 olculen 0.95; done
for v in Hkeskin son;     do .venv/bin/python qm_riemann/165_configs/165_kanal.py $v 0.60 0.40 olculen; done
.venv/bin/python qm_riemann/165_configs/165_dogrulama.py
.venv/bin/python qm_riemann/165_configs/165_analiz.py
```
