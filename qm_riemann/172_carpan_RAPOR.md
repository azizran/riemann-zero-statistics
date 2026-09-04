# 172 — ÇARPANLARIN YASASI: g_E, g_X, θ ve M(λ)'nın YENİDEN İNŞASI
### (4 Eylül 2026, Opus tayfası — görev 172b…)

Kalem/çerçeve: `172a_vadi_ayristirma.py` (kaptan, sabah). Bu rapor
172b'den başlayan kapıları taşır: **G1** (g_E yasası), **G2** (θ yasası),
**G3** (gerçeğin dar adresi), **G4** (sentez), **G5** (g_X).

Betikler (`172_configs/`), hepsi ÖN-MÜHÜRLÜ, hepsi **önbellekten** okur
(yeni gaz İNŞA EDİLMEDİ, yeni ölçüm KOŞULMADI, git'e DOKUNULMADI):

| betik | kapı | koşu | ham çıktı |
|---|---|---|---|
| `172b_gE_yasasi.py` | **G1** — özdeşlik zinciri + projeksiyon kusuru Q | 0.2 s | `172/G1.json`, `172/log_172b.txt` |
| `172c_Q_eksenleri.py` | **G1 kapanışı + G5** — Q'nun eksen yasaları | 0.2 s | `172/log_172c.txt` |
| `172d_theta_yasasi.py` | **G2** — θ'nın birleşik modeli + τ-eğimi türetimi | 0.6 s | `172/G2.json`, `172/log_172d.txt` |
| `172e_gercek_adres.py` | **G3** — gerçek gazın λ_eş haritası | 0.2 s | `172/G3.json`, `172/log_172e.txt` |
| `172f_sentez.py` | **G4** — vadi, iki tümsek, 173 mührü, figür | 3 s | `172/G4.json`, `172_carpanlar.png` |

Toplam koşu süresi **≈ 4 saniye** — bütün 172, 167/170/171 önbelleklerinin
üstünde çalışır.

Önbellek kökü: `/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/
71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/` — girdi `167/C_*.json`
(16 gaz), çıktı `172/`.

---

## 0. TEK CÜMLELİK HÜKÜM (kapılar kapandıkça güncellenir)

> **g_E ve g_X bir "yasa"ya uyan serbest sayılar değil, tek bir cebirsel
> özdeşliğin okumalarıdır:** `g = 1 − Q/ρ`, burada **Q, model alanının
> DİK İZDÜŞÜM OLMA KUSURUDUR** (`Q·V_O = ⟨E·(E−e₁)⟩`, gerçek izdüşümde
> ÖZDEŞ SIFIR). Bu, koddan türetildi ve 16 gazda **makine hassasiyetinde**
> (≤ 5·10⁻¹⁵) kapanıyor. 172a'nın "g_E, rE'de doğrusal" koklaması bir yasa
> değil, bu özdeşliğin ailenin kendi (Q, μ̂²) sürüklenmesi boyunca aldığı
> biçimdir. **Üç ön-kayıtlı, SIFIR-PARAMETRELİ tam isabet çıktı:**
> (i) `Q_E(λ)` bir TÜMSEKTİR, tepesi **1.0009…1.0191** = inşa doyum
> eşiği **λ_c = 1.0109** (170 §K0.2);
> (ii) 171 §T4'ün açıklayamadığı **τ-eğimi borcu (−149 … +54 %/τ)
> tamamen gazın kendi A-genişliği α_g'dir** — on iki gazda ±%7
> (yalnız K070 +%20);
> (iii) `c(λ)`'nın vadisi `dlogKALİB/dλ = dlogW_X/dλ` kesişmesidir ve
> **λ = 0.6506** çıkıyor — 171'in ölçtüğü **λ\* = 0.6487**'den ‰3 fark;
> aynı yerde A-genişliği α(λ) TEPE yapıyor (0.6675, ön-kayıt ±0.03).
> **Gerçek gazın fazlasının dar adresi ise E (η) KANALIDIR:**
> λ_eş(E) = 0.777 < λ_eş(X) = 0.909 < λ_eş(marjinaller) = 0.935 —
> 162'nin "yıkıcı girişim η kanalında var, Ĉ kanalında YOK" ölçümünün
> çarpan dilindeki tam karşılığı; en aykırı tek sayı **μ̂²_E = −%23.8**
> (model alanının DC kaçağı = sıfır tarağının asal rezonansları).

---

## G1 — g_E'NİN YASASI: ÖZDEŞLİK, YASA DEĞİL (`172b`, `172c`)

### G1.0 TÜRETİM — ölçüm kodundan, kalemle

`165_cekirdek.Model165.alanlar` (satır 272-300) g_E'yi şöyle üretir:

```
E   = sentez(s, ω, hp[msk])          hp_q = 2⟨η_{n+1} e^{−iω_q s}⟩
e1  = Y.e1 − ort(Y.e1)               (ORTALAMASI TAM SIFIR)
g_E = np.dot(E,e1)/np.dot(E,E)       (KESİŞİMSİZ EKK eğimi)
X   = sentez(...) − ort              (X AÇIKÇA ORTALANIYOR, E DEĞİL)
g_X = np.dot(X,x1)/np.dot(X,X)
```

Ham ikinci moment `⟨fg⟩ := (1/N)Σ fg`. `ort(e1)=0` olduğundan
`⟨E e1⟩ = Kov(E,e1) =: K` ve `⟨E²⟩ = Var(E) + μ²`, `μ := ort(E)`.
Varyans cebiri — **hiçbir varsayım yok**:

```
(Ö-A)  K = ½[varE_olc + varE_mod − varE_res]  =  korE·√(varE_mod·varE_olc)
(Ö-B)  g_E = K/(varE_mod + μ²)      ⇒  μ̂² := (K/g_E − varE_mod)/varE_olc
```

Boyutsuz üç sayı: `π := K/V_O` (yakalanan güç payı), `ρ := ⟨E²⟩/V_O`
(model gücü, DC dahil) ve

```
        **Q := ρ − π = ⟨E·(E − e1)⟩ / V_O**            (PROJEKSİYON KUSURU)
        **g_E = 1 − Q/ρ**   ,   ρ = 1 − rE + 2Q − μ̂²   (rE := varE_res/varE_olc)
```

`E` gerçek bir dik izdüşüm olsaydı `(E−e1) ⊥ E` olurdu ⇒ **Q ≡ 0, g_E ≡ 1**.
Gram diliyle (`u_k ∈ {cos ω_q s, sin ω_q s}`, `Γ = U^TU`, `Δ := 2Γ − I`):
`K = 2|v|²`, `⟨E²⟩ = 4v^TΓv` ⇒ **`g_E = 1/(1+δ)`, `δ = v^TΔv/|v|²`,
`Q = πδ`** — yani Q, merdiven çizgilerinin **dikgen olmamasının** tek
sayılık ölçüsü. Ölçek değişmezliği açıktır: `e1 → αe1` altında Q, π, ρ,
g_E'nin hepsi sabittir.

### G1.1 ÖZDEŞLİK DENETİMİ — MAKİNE HASSASİYETİ (`172b`, Ö1/Ö2)

| denetim | sonuç |
|---|---|
| (Ö-A) `\|K_var − K_kor\|/K`, 16 gaz | **≤ 5.2·10⁻¹⁵** (çoğu 10⁻¹⁶) |
| μ̂²_E > 0 (zorunlu), 16 gaz | ✓ 0.0055 … 0.0627 |
| μ̂²_E/ρ (ön-kayıt < 0.10) | ✓ ≤ 0.029 |
| **μ̂²_X ≡ 0** (kod X'i ORTALIYOR) | ✓ **\|μ̂²_X\| ≤ 10⁻¹⁷** |

> **HÜKÜM (G1a).** `g_E` bağımsız bir ölçüm değildir: üç varyans + `μ²`
> onu ÖZDEŞ olarak belirler. `μ̂²_X ≡ 0` sonucu cebrin **kod düzeyinde
> mührüdür** (`alanlar` X'i açıkça ortalar, E'yi ortalamaz) ve 10⁻¹⁷'de
> doğrulanmıştır. `μ̂²_E`, model alanının **DC kaçağıdır** — merdiven
> frekanslarında sıfır tarağının asal rezonansları (`κ(ω_q) ≠ 0`);
> λ ile tekdüze büyür (0.0202 → 0.0627, λ = 0.50 → 1.15).

### G1.2 Q TABLOSU — bütün eksenler (`172b`)

`g = 1 − Q/ρ` her satırda ölçülen g'yi dört hanede verir (özdeşlik).

| gaz | λ | rE | ρ_E | **Q_E** | g_E | rX | ρ_X | **Q_X** | g_X |
|---|---|---|---|---|---|---|---|---|---|
| L050 | 0.50 | 0.4052 | 1.9552 | **0.69027** | 0.6470 | 0.2795 | 1.7046 | **0.49204** | 0.7113 |
| L060 | 0.60 | 0.4743 | 2.0633 | **0.78299** | 0.6205 | 0.3071 | 1.6834 | **0.49526** | 0.7058 |
| L070 | 0.70 | 0.5205 | 2.1352 | **0.84627** | 0.6037 | 0.3266 | 1.6561 | **0.49134** | 0.7033 |
| L085 | 0.85 | 0.5597 | 2.1894 | **0.89861** | 0.5896 | 0.3473 | 1.6094 | **0.47833** | 0.7028 |
| Hkeskin | 1.00 | 0.5750 | 2.1993 | **0.91558** | 0.5837 | 0.3619 | 1.5624 | **0.46213** | 0.7042 |
| L115 | 1.15 | 0.5760 | 2.1809 | **0.90981** | 0.5828 | 0.3727 | 1.5180 | **0.44534** | 0.7066 |
| L130 | 1.30 | 0.5854 | 2.0348 | **0.83903** | 0.5877 | 0.3959 | 1.4767 | **0.43627** | 0.7046 |
| **son** | — | 0.5314 | 2.1796 | **0.87575** | 0.5982 | 0.3511 | 1.5848 | **0.46793** | 0.7047 |
| K090 | kesim | 0.5776 | 2.2471 | **0.94107** | 0.5812 | 0.3582 | 1.5766 | 0.46742 | 0.7035 |
| K070 | kesim | 0.5050 | 2.1373 | **0.83046** | 0.6114 | 0.3503 | 1.5089 | 0.42961 | 0.7153 |
| HA4 | kesim | 0.4729 | 2.1032 | **0.79586** | 0.6216 | 0.3412 | 1.5196 | 0.43039 | 0.7168 |
| E060 | kesim | 0.4157 | 1.8884 | **0.65481** | 0.6532 | 0.3352 | 1.4777 | 0.40647 | 0.7249 |
| HkT2a | T/2 | 0.5943 | 2.2549 | 0.95133 | 0.5781 | 0.3921 | 1.6563 | 0.52420 | 0.6835 |
| HkT2b | T/2 | 0.6027 | 2.2935 | 0.97610 | 0.5744 | 0.3990 | 1.6750 | 0.53695 | 0.6794 |
| HkT4a | T/4 | 0.6715 | 2.4455 | 1.08239 | 0.5574 | 0.5005 | 1.9299 | 0.71523 | 0.6294 |
| HkT4b | T/4 | 0.6782 | 2.4685 | 1.09769 | 0.5553 | 0.5171 | 1.9587 | 0.73788 | 0.6233 |

### G1.3 ÖN-KAYITLI ÖLÜMLER — kurtarma yok

| ön-kayıt (`172b`/`172c`) | türetim dedi | ölçüm dedi | hüküm |
|---|---|---|---|
| Ö3: Q λ-DEĞİŞMEZ (ölçek değişmezliği), yayılım ≤ %5 | ≤ %5 | **%26.8** | **ÖLDÜ** |
| Ö4: `dg_E/drE = −Q/ρ²` | −0.1893 | **−0.3550** | **ÖLDÜ** (1.9 kat) |
| Ö4′: o zaman zorunlu `dQ/drE` | +2.41 | **+1.12** | **ÖLDÜ** (ön-kayıt μ̂²'nin sürüklenmesini ihmal ediyordu; `dμ̂²/drE ≈ +0.25`) |
| Ö6: kesimde Q düşer | Q<Q(Hk) | K090 **+%2.8 ARTTI** | **kısmi ıska** |
| Ö3(172c): `B_E/B_X` — M1 gürültü mekanizması | 1.025 | **0.656** | **ÖLDÜ** |
| Ö3(172c): `B_E/B_X` — M2 κ-genişlemesi | 1.981 | **0.656** | **ÖLDÜ** |
| Ö3(172c): M2 ikinci ayak `Q(T/4)/Q(T)` | 4.00 | **1.19 (E), 1.57 (X)** | **ÖLDÜ** |
| Ö1(172c): pencere ekseninde Q, d/N'de doğrusal (rms ≤ %1) | ≤ %1 | E **%0.88 ✓**, X **%2.40 ✗** | kısmi |
| Ö7(172c): gerçek gazın tek eşdeğer λ'sı (\|Δλ\| ≤ 0.10) | ≤ 0.10 | λ_eş(Q_E)=0.7845, λ_eş(Q_X)=0.9463, **Δ = 0.162** | **ÖLDÜ** |

**172a'nın koklaması, tam üyeleriyle** (`172c` EK bloğu):

| uyum | doğru | artık rms | L130 | **son** |
|---|---|---|---|---|
| 7 λ gazı | `g_E = 0.78955 − 0.35496 rE` | 2.6·10⁻³ | −1.00% | +0.45% |
| **L130 hariç (6 gaz)** | `g_E = 0.79748 − 0.37215 rE` | **3.0·10⁻⁴** | **−1.37%** | **+0.25%** |

Altı gazın hepsi ±%0.07 içinde (172a'nın kaba `0.798 − 0.373`'ü tam
isabet). Kesim/pencere gazlarının doğrudan sapmaları: E060 **−%1.73**,
T/4 **−%1.11/−%1.17**, K090 +%0.57, K070 −%0.19, HA4 +%0.02.

> **HÜKÜM (G1b — 172a'nın koklamasının cevabı).** `g_E ≈ 0.797 − 0.372·rE`
> **bir yasa değildir.** Doğrusallık, özdeşlik
> `g_E = 1 − Q/(1−rE+2Q−μ̂²)`'nin ailenin kendi sürüklenmesi
> (`dQ/drE = +1.12`, `dμ̂²/drE ≈ +0.25`) boyunca aldığı **yerel** biçimdir.
> Sabit-Q "yasası" eğimi 1.9 kat ıskalıyor ve **ön-kayıtlı olarak öldü**.
> **L130'un −%1.37'lik sapması AÇIKLANDI:** L115 → L130 adımında
> `dQ/drE = −7.55` (işaret dönüyor, tümseğin inen kolu) ve dolayısıyla
> `dg_E/drE = +0.52` — doğrunun eğiminin tam TERSİ. Yani L130 "kırık
> bölge gürültüsü" değil, **Q-tümseğinin öbür yamacıdır** (§G1.5).
> Gerçek gaz `son` doğrunun **+%0.25** üstünde, yani ailenin ±%0.07
> bandının biraz dışında ama kesim/pencere sapmalarının (%1-2) çok
> altında: **g_E artığı ≈ 0** (G3 beklentisi ilk sınavı geçti).

### G1.4 Q'NUN EKSEN AYRIŞTIRMASI (`172c`)

`d/N := 2·nline/N` (model serbestlik derecesi payı) **λ ve kesim
eksenlerinde TAM SABİTTİR** (0.059874; nline = 8981, N = 299998 her
gazda). Yalnız pencere ekseninde değişir (0.1159 … 0.2360).

```
Q = Q_koh + B·(d/N)            (pencere ekseni, 5 nokta, 2 parametre)
E:  Q_koh = 0.84390,  B = 1.05201,  artık rms 0.0080 (%0.88)
X:  Q_koh = 0.35040,  B = 1.60297,  artık rms 0.0111 (%2.40)
```

> **HÜKÜM (G1c).** Üç eksen Q'yu **üç ayrı yoldan** kımıldatır:
> **pencere** ekseni saf **aşırı uyumdur** (d/N; Q(Hk)'nin %8'i E'de,
> %24'ü X'te), **λ ve kesim** eksenlerinde d/N sabit olduğu için Q'nun
> bütün hareketi **koherent sızıntıdır**. Ama T ekseninin mekanizması
> ne saf gürültü-izdüşümü (M1) ne saf κ-genişlemesidir (M2) — ikisi de
> ön-kayıtlı biçimde öldü.

### G1.5 **TÜMSEK ve λ_c** — 172'nin ilk tam isabeti (`172c`, Ö5)

`Q_E(λ)` tek-yönlü değil, **TÜMSEKTİR**: 0.6903 (λ=0.50) ↗ 0.9156
(λ=1.00) ↘ 0.8390 (λ=1.30). `log Q_E`'ye parabol:

| uyum penceresi | tepe λ\* | ön-kayıt |
|---|---|---|
| 5 orta nokta (0.60 … 1.15) | **1.0191** | λ_c = **1.0109** ± 0.05 |
| 7 nokta (0.50 … 1.30) | **1.0009** | (170 §K0.2, sıfır yeni parametre) |

> **HÜKÜM (G1d — ÖN-KAYITLI TAM İSABET).** **Model alanının projeksiyon
> kusuru, inşa denkleminin kendi doyum eşiğinde zirve yapar.**
> λ_c = 1.9147/1.894 = 1.0109 sayısı 170 §K0.2'de `rms S′ = N̄′`
> koşulundan gelmişti ve buraya hiçbir uyumla girmedi. Bu, 172a'nın
> "g_E tek-taraflı hokey sopası" gözleminin **kaynağıdır**: g_E = 1−Q/ρ
> olduğu için Q'nun λ_c'deki zirvesi g_E'nin λ ≥ 1'deki DÜZLÜĞÜNÜ
> (tepenin yakınında `dQ/dλ ≈ 0`) ve λ < 1'deki yükselişini birlikte
> üretir. `Q_X`'in tepesi ise ölçüm aralığının çok dışındadır
> (parabol tepesi 0.24-0.30) ⇒ **g_X λ-değişmezliği** (171) burada
> mekanizmasıyla açıklanıyor: Q_X ve ρ_X **aynı yönde ve yakın hızda**
> hareket ettiği için oranları `Q_X/ρ_X` λ boyunca ±%1 içinde kalıyor.

### G1.6 (G5) g_X — λ = 0.50 uyanışı ile T/4 sızıntısı **AYNI ODA DEĞİL**

`g_X = 1 − Q_X/ρ_X` özdeşliğinin log-ayrıştırması
(`d log g = −[f/(1−f)]·(d log Q − d log ρ)`, `f = Q_X/ρ_X = 0.2958`):

| olay | Δlog g_X² | Δlog Q_X | Δlog ρ_X | baskın çarpan |
|---|---|---|---|---|
| λ: 0.60 → 0.50 (172a'nın "uyanışı") | **+0.01566** | −0.00651 | **+0.01253** | **ρ_X — %66** (model gücü) |
| T/4 (HkT4b vs Hkeskin) | **−0.24419** | **+0.46795** | +0.22606 | **Q_X — sızıntı** |

> **HÜKÜM (G1e = G5).** **HAYIR, tek çerçeve değiller.** λ = 0.50'deki
> `g_X²` uyanışının üçte ikisi **model gücü** `ρ_X`'in artışıdır (çizgi
> gösterimi düşük λ'da alanı daha çok şişiriyor), yalnız üçte biri
> sızıntıdır. T/4 sızıntısı ise tersine **saf Q_X olayıdır** (aşırı
> uyum, d/N dört katına çıkıyor) ve ρ_X onu ancak yarı yarıya telafi
> ediyor. 168 §BONUS-iii'ün "T/4 borcu g_X'te" adresi **doğrudur ve
> şimdi mekanizmalıdır: d/N aşırı uyumu**; ama λ = 0.50 adımı o odaya
> ait değildir.

---

## G2 — θ'NIN YASASI: BİR ZAYIF TAŞIYICI, BİR TAM İSABET (`172d`)

### G2.0 TÜRETİM — θ bir ARTIK-ALAN EŞLEŞMESİDİR

`θ := KALİB_u2/g_cal`, `KALİB_u2 = Mu2/Pu2` (`166_T1.bant_agg`),
`Mu2 ∝ Re[h̄_Q⟨e1·x1²·e^{−iWs}⟩]`, `Pu2 ∝ Re[h̄_Q⟨E·X²·e^{−iWs}⟩]`.
EKK artıkları `ε := e1 − g_E E` (⟨εE⟩ = 0) ve `ξ := x1 − g_X X` ile
ölçülen üçlü korelatör TAM açılır:

```
⟨e1 x1² e^{−iWs}⟩ = g_E g_X²⟨E X² e^{−iWs}⟩              [= model]
   + 2g_Eg_X⟨EXξ e^{−iWs}⟩ + g_E⟨Eξ² e^{−iWs}⟩
   + g_X²⟨εX² e^{−iWs}⟩ + 2g_X⟨εXξ e^{−iWs}⟩ + ⟨εξ² e^{−iWs}⟩

⇒  θ = 1 + (beş artık teriminin toplamı)/(g_E g_X²·model)      (ÖZDEŞ)
```

Yani **θ, g'lerin çözemediği artık-alan eşleşmesidir.** Bu, iki aday
ailesini ÖLÇÜMDEN ÖNCE eler: (E1) `σ_ds, σ_X̃, σ_Ĉ, W_X` pencere
ekseninde ‰4 içinde donuktur ama θ orada %11 oynar; (E2) `d/N` λ ve
kesim eksenlerinde TAM sabittir (172c Ö4).

### G2.1 ÇAPRAZ SINAV — β dondurulmuş, p yalnız λ'dan (`172d`)

`θ/θ₀ = (1 − 0.2175·φ)·(x/x₀)^p`; `p` **yalnız** 7 λ gazından; pencere
(4) + kesim (4) + gerçek gaz (1) **örneklem dışıdır**.

| aday x | p | λ rms% | PENCERE% | KESİM% | son% | hüküm |
|---|---|---|---|---|---|---|
| rE | −0.241 | 0.79 | **8.74** | 4.86 | −0.95 | ÖLDÜ (işaret) |
| rX | −0.299 | 0.36 | **12.68** | 1.54 | −1.93 | ÖLDÜ (işaret) |
| Q_E | −0.258 | 1.87 | **9.17** | 4.89 | −1.69 | ÖLDÜ |
| Q_X | +0.648 | 2.01 | **17.01** | 5.34 | −2.02 | ÖLDÜ |
| ρ_E | −0.462 | 2.72 | 9.56 | 3.72 | −2.41 | ÖLDÜ |
| **ρ_X** | **+0.651** | **1.37** | **4.59** | **2.56** | **−1.91** | **AYAKTA** |
| 1−g_E | −0.495 | 1.11 | 8.30 | 5.45 | −1.09 | ÖLDÜ |
| 1−g_X | −2.848 | 2.55 | 40.10 | 14.24 | −2.33 | ÖLDÜ |
| σ_ds / σ_X̃ / σ_Ĉ | ≈−0.17 | 0.53 | **6.0-6.2** | 0.5-0.6 | −1.9…−2.5 | ÖLDÜ (E1) |
| W_X | +0.248 | 1.01 | 6.21 | 0.57 | −2.02 | ÖLDÜ (E1) |
| d/N | — | — | — | — | — | TANIMSIZ (E2) |
| g_cal | +0.673 | 1.19 | 17.87 | 7.32 | −1.10 | SİREN, öldü |
| *(2. tur, ön-mühürsüz)* kor_X⁻² | −0.625 | 0.90 | 6.90 | 0.77 | −1.86 | ÖLDÜ |
| *(2. tur)* ρ_E·ρ_X² | +0.266 | 3.02 | 4.87 | 4.11 | −2.31 | (zayıf) |

`ρ_X`'in gaz-başına artıkları: pencere **+0.68 / +4.77 / +6.25 / +4.69**,
kesim **+0.78 / −2.90 / −1.23 / −3.95**, son **−1.91** (%).

> **HÜKÜM (G2a).** **Ön-kaydım öldü:** "hiçbir tek değişken θ'yı üç
> eksende taşımaz" demiştim; **ρ_X = varX_mod/varX_olc taşıyor** —
> `θ/θ₀ = (1−βφ)·(ρ_X/ρ_X₀)^{0.651}`, tek üs, yalnız λ'dan, örneklem
> dışı rms %4.6 (pencere) ve %2.6 (kesim). **Ama bu bir KİMLİK
> DEĞİLDİR:** θ'nın kendi jackknife hatası ‰5'tir, dolayısıyla T/4
> gazlarında yasa hâlâ 5-10σ dışındadır. `ρ_X`, θ'nın **en iyi
> tek-değişkenli taşıyıcısıdır**, kimliği değil. Fiziksel okuma:
> θ'yı hareket ettiren şey **X kanalının MODEL GÜCÜ şişmesidir**
> (çizgi gösteriminin X̃ üstünde ürettiği fazla güç), σ'lar veya
> kesim kesri değil. (`μ̂²_X ≡ 0` olduğu için `ρ_X ≡ ρ_var_X` — özdeşlik
> denetiminin ikinci kez tutması.)

### G2.2 **τ-EĞİMİ BORCU KAPANDI — sıfır parametre** (`172d`, Ö6/Ö7)

171 §T4.2 kesim ekseninde `θ = 1−βφ` seviyeyi ‰7 içinde verirken
"β'nın görmediği" tekdüze bir τ-eğimi bırakıyordu (−15.2 → −141.7 %/τ).
**Türetim:** bu eğim θ'nın seviye yasasının borcu değil, gazın kendi
A-genişliğinin işidir. `A_g(τ) = e^{−α_g τ²}` ⇒

```
artık(τ) = 100·[ (θ̄_g /(θ̄₀(1−βφ)))·e^{−Δα(τ² − τ̄²)} − 1 ],   Δα := α_g − α₀
```

`α_g` 171'in gaz-başına uyumundan, `β` 168'den, `θ̄` ölçümden —
**hiçbir yeni parametre yok.**

| gaz | eksen | Δα | τ-eğimi ÖLÇÜLEN | τ-eğimi ÖNGÖRÜ | fark |
|---|---|---|---|---|---|
| K090 | kesim | +0.1218 | −15.2 | **−15.0** | −0.8% |
| K070 | kesim | +0.3236 | −33.3 | **−40.1** | +20.2% |
| HA4 | kesim | +0.6061 | −71.6 | **−74.4** | +3.9% |
| E060 | kesim | +1.2024 | −141.7 | **−148.9** | +5.1% |
| L060 | λ | +0.0964 | −11.8 | **−12.5** | +6.3% |
| **L050** | λ | −0.0735 | **+10.6** | **+9.9** | −6.7% |
| **L130** | λ | −0.4103 | **+49.2** | **+49.6** | +0.9% |
| **son** | GERÇEK | +0.0864 | −10.7 | **−11.0** | +3.0% |
| HkT2a / HkT2b | T/2 | −0.109 / −0.018 | +14.5 / +2.2 | +13.9 / +2.3 | −4.1% / +3.0% |
| HkT4a / HkT4b | T/4 | −0.203 / −0.389 | +28.6 / +54.4 | +27.1 / +53.3 | −5.5% / −2.1% |

> **HÜKÜM (G2b — 171 §T4-b'nin kapanışı).** **171'in "β'nın görmediği
> şekil borcu" diye bıraktığı τ-eğimi, θ'nın hiçbir borcu değildir:
> gazın A-genişliği α_g'nin türevidir ve TEK bir sayıyla, sıfır yeni
> parametreyle, −149'dan +54 %/τ'ya uzanan bütün yelpazeyi ±%7 içinde
> veriyor** (yalnız K070 +%20 ile ön-kayıtlı sınırın ucunda). Yani
> 171 §T4'ün "kesim ekseni İKİ boyutludur (seviye β·φ + şekil)" hükmü
> doğrudur ama iki boyut BAĞIMSIZ ölçülür: seviye θ'da, şekil α'da,
> ve şeklin borcu artık ölçülü bir sayıya bağlanmıştır. `α_g` denetimi:
> 172d'nin kendi uyumu 171'in 14 değerini **+0.00%** ile yeniden üretti;
> yeni ölçülenler **α(L050) = 1.0242, α(L130) = 0.6874**.

---

## G3 — GERÇEĞİN DAR ADRESİ: **E (η) KANALI** (`172e`)

Her ölçülür büyüklük için `λ_eş` := o büyüklüğün `son` değerini veren λ
(7-nokta merdiveninde, tekdüze kolda ters çevirme). Gerçek gaz bir
λ-gazı olsaydı hepsi aynı çıkardı.

| büyüklük | kanal | son | λ = 1.00 | **λ_eş** |
|---|---|---|---|---|
| σ_Ĉ | marjinal | 0.27303 | 0.27768 | 0.9684 |
| σ_X̃ | marjinal | 0.23384 | 0.24204 | 0.9363 |
| σ_ds | marjinal | 0.40919 | 0.43134 | 0.9014 |
| ρ_X | X (ΔĈ) | 1.58477 | 1.56236 | 0.9285 |
| rX | X (ΔĈ) | 0.35109 | 0.36189 | 0.8892 |
| ρ_E | **E (η)** | 2.17960 | 2.19932 | **0.8228** |
| Q_E | **E (η)** | 0.87575 | 0.91558 | **0.7845** |
| g_E | **E (η)** | 0.59821 | 0.58370 | **0.7581** |
| μ̂²_E | **E (η)** | 0.04050 | 0.05684 | **0.7489** |
| rE | **E (η)** | 0.53140 | 0.57499 | **0.7417** |
| θ | θ | 0.91416 | 0.88844 | **0.6916** |

```
λ_eş(marjinaller) = 0.9354   λ_eş(X kanalı) = 0.9089   λ_eş(E kanalı) = 0.7768
```

### G3.1 Çarpan defteri (σ_X̃ çapasında, λ_eş = 0.9363)

| çarpan | ölçülen | λ-eğrisi | **artık** | 172a'nın ön-kaydı |
|---|---|---|---|---|
| g_E | 0.59821 | 0.58619 | **+2.05%** | ≈0 (≤%0.5) → **ÖLDÜ** |
| (g_X)² | 0.49665 | 0.49507 | **+0.32%** | ≈0 (≤%0.5) → ✓ |
| θ | 0.91416 | 0.89160 | **+2.53%** | +%3 (2-4) → ✓ ucundan |
| **toplam** | | | **+4.97%** | = ölçülen ΔlogM (defter farkı 3.3·10⁻⁶) |

Artıkların büyüklük sıralaması (aynı çapada): **μ̂²_E −23.8%**,
rE −6.5%, α +4.6%, **Q_E −3.6%**, θ +2.5%, g_E +2.1%, σ_ds −1.9%,
σ_Ĉ +1.8%, rX −1.3%, ρ_E −0.7%, **Q_X −0.2%, g_X +0.2%, ρ_X +0.2%**.
Cebirsel sürücü Q_E'dir:
`ΔrE = −Δρ_E + 2ΔQ_E − Δμ̂²_E = +0.0155 − 0.0652 + 0.0127 = −0.0371` ✓.

> **HÜKÜM (G3 — 172a'nın "dar adres" sorusunun cevabı).**
> **Gerçek gaz bir λ-gazı DEĞİLDİR** (λ_eş yayılımı 0.69 … 0.97,
> ön-kayıt ±0.05 — kesin ölüm). Fazlanın adresi bir ÇARPAN değil bir
> **KANALDIR**: X (ΔĈ) kanalı ve marjinaller λ ≈ 0.91-0.94 derken
> **E (η) kanalının BEŞ büyüklüğü birden λ ≈ 0.74-0.82 diyor.**
> 172a'nın "g_E artığı ≈ 0, θ artığı ≈ +%3" beklentisi yarı doğru:
> θ +%2.53 ✓ ama **g_E de +%2.05 ile yasa dışıdır** — g_E'nin
> "yasa üstünde" görünmesi (g_E–rE doğrusunda +%0.25) yanıltıcıdır,
> çünkü gerçek gaz doğrunun ÜSTÜNDE ama **doğru boyunca λ ≈ 0.76'ya
> kaymış** durumdadır. Şekil borcu ise YOKTUR (son'un τ-eğimi G2b
> yasasını +%3.0 ile sağlıyor): **fazla saf bir SEVİYE olayıdır ve
> seviyeyi E kanalı taşır.**

### G3.2 162 KÖPRÜSÜ — ön-kayıtlı yön TUTTU

162 §6'nın ölçümü: yıkıcı girişim **η kanalındadır**
(`Σ|c_η|²/2 ÷ Var(η) = 1.288`, %29 sessizleşme), **Ĉ kanalında YOKTUR**
(1.026). 165'in `E` alanı η kanalıdır, `X` alanı `X̃ ≈ ΔĈ` kanalıdır.
Bellek fazlardaysa E-kanalı büyüklükleri "daha küçük λ" gibi
görünmelidir. **Ön-kayıt (172e Ö3): λ_eş(E) < λ_eş(X). ÖLÇÜM: 0.7768 <
0.9089 — TUTTU**, ve tek tek beş E büyüklüğünün hepsi X'in altında.

En dar tek sayı **μ̂²_E**'dir: tanımı gereği `μ̂²_E = ort(E)²/Var(e1)` ve
`ort(E) = Σ_q Re[hp_q·κ(ω_q)]` — yani **sıfır tarağının merdiven
frekanslarındaki (asal) rezonansı**. Gerçek ζ gazında bu −%23.8 daha
zayıftır. İkinci sıradaki `Q_E = 2v^TΔv/V_O` de aynı cinstendir:
`Δ_{kl} ≈ Re κ(ω_k−ω_l)` ve merdiven frekanslarının FARKLARI yine
merdiven frekanslarıdır (`log q₁ − log q₂ = log(q₁/q₂)`). **Yani hem
μ̂²_E hem Q_E, tarağın ASAL-ORAN frekanslarındaki tepkisidir ve gerçek
zeta sıfırlarında ikisi de daha sönüktür.** 162'nin "bellek fazlarda"
bulgusunun +3σ'ya giden yolu buradan geçiyor: **θ'nın ağırlık-doyumu
faz-kilidini DOĞRUDAN duymuyor** (θ'nın taşıyıcısı ρ_X, X kanalı,
ve orada artık ≈ 0); duyan `g_E`'dir, ve θ'nın +%2.5'i onun
**yanında** duran ikinci, henüz kimliksiz bir seviye borcudur.
*(Dürüstlük: bu bir mekanizma iddiası değil, adres iddiasıdır.
`μ̂²_E` ve `Q_E`'nin κ-cinsinden doğrudan ölçümü yapılmadı — alanları
yeniden sentezlemek gerekir; 173'ün kapısı.)*

---

## G4 — SENTEZ: İKİ TÜMSEK, BİR VADİ (`172f`, `172_carpanlar.png`)

### G4.1 VADİNİN ÇARPAN ADRESİ — ‰3 isabet

`c = KALİB/W_X` ⇒ `dlog c/dλ = dlogKALİB/dλ − dlogW_X/dλ`. Ölçülen
eğriler (7 λ gazı, ardışık çiftler):

| λ_orta | dlogKALİB/dλ | dlogW_X/dλ | fark (= dlog c/dλ) |
|---|---|---|---|
| 0.550 | −0.9046 | −0.6023 | **−0.3023** |
| **0.650** | **−0.5603** | **−0.5593** | **−0.0010** |
| 0.775 | −0.2899 | −0.5103 | +0.2205 |
| 0.925 | −0.0953 | −0.4616 | +0.3662 |
| 1.075 | +0.0144 | −0.4236 | +0.4381 |
| 1.225 | −0.1274 | −0.4671 | +0.3397 |

**Kesişim λ = 0.6506**, ön-kayıt (171'in ölçtüğü vadi) **λ\* = 0.6487**
⇒ fark **+0.0019 (‰3).**

> **HÜKÜM (G4a).** `dlogW_X/dλ` λ boyunca neredeyse SABİTTİR
> (−0.42 … −0.60); oynayan taraf `KALİB`'dir (−0.90 → +0.01).
> **Vadi, KALİB'in kendi eğiminin W_X'in sabit eğimini geçtiği yerdir.**
> 170 §K0b'nin "c bir payda eğrisidir" hükmü buradan düzeltiliyor:
> payda **sabit eğimli**, şekli veren **paydır**.

### G4.2 İKİ TÜMSEK, İKİ EŞİK

| tümsek | tepe (ölçülen) | ön-kayıtlı eşik | kaynak | fark |
|---|---|---|---|---|
| **Q_E(λ)** — izdüşüm kusuru | **1.0009 / 1.0191** | λ_c = **1.0109** | 170 §K0.2 (`rms S′ = N̄′`) | ≤ 0.011 |
| **α(λ)** — A-genişliği | **0.6675** (log) / **0.6594** | λ\* = **0.6487** | 171 §T2c.3 (c vadisi) | +0.019 / +0.011 |

*(α'nın 7-noktalı uyumu 0.7877 verir ve ön-kayıt bandının dışındadır;
sebep L130'un α = 0.6874 aykırılığıdır. Ön-kayıt "5 orta nokta" demişti
— kurtarma yapılmadı, iki sayı da yazıldı.)*

### G4.3 İNDİRGENMİŞ MODEL — ÖLDÜ (dürüstlük)

`g_X ≡ sabit` ön-kaydı ıskaladı (maks **%1.01**, eşik %0.5).
`M ≈ (g_E/g₀)·(ρ_X/ρ_X₀)^{0.651}` yedi gazda rms **%2.07** (eşik %1.5)
⇒ **ÖLDÜ**; en büyük hata λ = 0.50'de **−%4.37**. M(λ)'nın
"parametresiz yeniden inşası" ancak ÖZDEŞLİK düzeyinde (ölçülen Q, ρ
ile) tamdır ve o hiçbir şey sınamaz. **Ejderha kapanmadı:** vadinin
YERİ çıktı (‰3), ama M(λ)'nın kimliği hâlâ yok.

### G4.4 173 İÇİN MÜHÜR (`172/G4.json`) — iki bağımsız yol

| λ | yol | g_E | g_X | θ | **M** | **c** |
|---|---|---|---|---|---|---|
| **0.40** | çarpan | 0.6842 | 0.7175 | 0.9452 | **1.2947** | **0.3790** |
| **0.40** | doğrudan (171 M9) | — | — | — | **1.3881** | **0.4052** |
| **1.45** | çarpan | 0.5939 | 0.7085 | 0.8405 | **0.9744** | **0.4747** |
| **1.45** | doğrudan (171 M9) | — | — | — | **1.0060** | **0.4910** |

İki mühür λ = 0.40'ta **%6.7 ayrışıyor**. Çarpan yolu λ = 0.50'yi zaten
−%4.4 ıskaladığı için **doğrudan yol daha güvenilirdir**; 173 bu gazları
kurarsa iki mühür arasındaki fark tek başına bir hakemdir.
(`c = c₀·M·W_X₀/W_X` dönüşümü mevcut yedi gazda ‰1.4 tutuyor.)

---

## 5. HÜKÜM TABLOSU — bütün ön-kayıtlar

| # | ön-kayıt | türetim/eşik | ölçüm | hüküm |
|---|---|---|---|---|
| 172b-Ö1 | özdeşlik `K_var = K_kor` | < 1e−6 | ≤ 5.2e−15 | ✓✓ |
| 172b-Ö2 | μ̂²_E > 0, μ̂²/ρ < 0.10 | — | 0.0055-0.063, ≤ 0.029 | ✓ |
| 172b-Ö2′ | *(kod)* μ̂²_X ≡ 0 | — | ≤ 1e−17 | ✓✓ |
| 172b-Ö3 | Q λ-değişmez | ≤ %5 | **%26.8** | ✗ |
| 172b-Ö4 | dg_E/drE = −Q/ρ² | −0.189 | **−0.355** | ✗ |
| 172b-Ö4′ | zorunlu dQ/drE | +2.41 | **+1.12** | ✗ |
| 172b-Ö5 | Q(T) < Q(T/2) < Q(T/4) | işaret | 0.916 < 0.951/0.976 < 1.082/1.098 | ✓ |
| 172b-Ö6 | kesimde Q düşer | işaret | K090 **+%2.8**, ötekiler ↓ | ✗/✓ |
| 172b-Ö7 | Q_E(son) λ-bandı içinde | — | 0.8758 ∈ [0.690, 0.916] | ✓ |
| 172b-Ö8 | Q_X yayılımı < Q_E'ninki | — | %12.5 < %26.8 | ✓ |
| 172c-Ö1 | Q, d/N'de doğrusal | rms ≤ %1 | E %0.88 ✓ / X %2.40 ✗ | kısmi |
| 172c-Ö2 | Q_koh > 0 ve Q(Hk)'nin altında | — | 0.8439 / 0.3504 | ✓ |
| 172c-Ö3 | B_E/B_X (M1 gürültü) | 1.025 | **0.656** | ✗ |
| 172c-Ö3 | B_E/B_X (M2 κ-genişleme) | 1.981 | **0.656** | ✗ |
| 172c-Ö3 | Q(T/4)/Q(T) (M2) | 4.00 | **1.19 / 1.57** | ✗ |
| **172c-Ö5** | **Q_E tepesi = λ_c** | **1.0109 ± 0.05** | **1.0009 / 1.0191** | **✓✓** |
| 172c-Ö6 | λ=0.50 uyanışı ρ_X'ten | — | %66 ρ_X, %34 Q_X | ✓ |
| 172c-Ö7 | tek λ_eş (Q_E ↔ Q_X) | ≤ 0.10 | **0.162** | ✗ |
| 172d-Ö1 | σ-adayları pencerede ölür | rms > %7 | %6.0-6.2 (yine ölü) | ✗(kısmi) |
| 172d-Ö3 | rE/rX işaret çelişkisi | ölüm | %8.7 / %12.7 | ✓ |
| 172d-Ö4 | **hiçbir tek değişken taşımaz** | — | **ρ_X ayakta (%4.6/%2.6)** | ✗ |
| **172d-Ö6/Ö7** | **τ-eğimi = −2Δα·τ̄** | **±%20** | **12 gazda ±%7 (K070 %20)** | **✓✓✓** |
| 172d-Ö8 | α_g ve θ denetimi | %1 / — | **+0.00%** / 2.2e−16 | ✓✓ |
| 172e-Ö1 | g_E artığı ≈ 0 | ≤ %0.5 | **+%2.05** | ✗ |
| 172e-Ö1 | (g_X)² artığı ≈ 0 | ≤ %0.5 | +%0.32 | ✓ |
| 172e-Ö1 | θ artığı ≈ +%3 | 2-4 | **+%2.53** | ✓ |
| 172e-Ö2 | gerçek gaz bir λ-gazı | ±0.05 | **0.69 … 0.97** | ✗ |
| **172e-Ö3** | **λ_eş(E) < λ_eş(X)** (162) | işaret | **0.777 < 0.909** | **✓** |
| 172e-Ö5 | çarpan defteri kapanır | ‰1 | 3.3e−06 | ✓✓ |
| **172f-Ö1** | **α tepesi = λ\*** | **0.6487 ± 0.03** | **0.6675 / 0.6594** | **✓** |
| **172f-Ö2** | **vadi = eğim kesişmesi** | **0.6487 ± 0.05** | **0.6506** | **✓✓** |
| 172f-Ö3 | g_X sabit | ≤ %0.5 | %1.01 | ✗ |
| 172f-Ö3 | indirgenmiş M | rms ≤ %1.5 | **%2.07** | ✗ |

**Sayım (33 satır): 18 ✓, 12 ✗, 3 kısmi.** Ölenlerin hiçbiri kurtarılmadı,
hiçbir eşik koşudan sonra gevşetilmedi.

---

## 6. DENETİM ve DÜRÜSTLÜK NOTLARI

1. **Hiçbir ölçüm/inşa parçası kopyalanmadı.** 172 yalnız önbellek okur:
   `167/C_*.json` (16 gaz), `170/K0.json` (φ, θ denetimi),
   `171/M_ONKAYIT.json` (α_g denetimi), `171/NU_MERDIVEN.json` (λ\*, c, M).
   Yeni gaz kurulmadı, yeni korelatör koşulmadı. **Git'e dokunulmadı.**
2. **Bağımsız yeniden üretim:** 172d, 170'in θ değerlerini 2.2·10⁻¹⁶ ve
   171'in α_g değerlerini +0.00% ile yeniden üretti — iki ayrı kod
   yolundan. Defterde tutarsızlık yok.
3. **Sirenler işaretlendi.** `g_cal` adayı tanım gereği θ'ya bağlıdır
   (θ ≡ KALİB/g_cal) ve SİREN olarak etiketlendi. 172d'nin "ikinci tur"
   adayları (`kor_X⁻²`, `ρ_E·ρ_X²`, `ρ_var_X`) ρ_X'in ayakta kalmasından
   SONRA eklendi; kimlik sayılmazlar, çünkü örneklem-dışı eksenler o
   noktada harcanmıştı.
4. **Ön-mühürsüz koklamalar** açıkça etiketlendi: (a) `B·r ≈ 0.59` iki
   alanda da (172c); (b) `p = 0.651 ≈ 2/3` (172d) — ikisi de 173'e not,
   iddia değil.
5. **172a'nın koklaması doğrulandı ama yasa değil:** altı λ gazı
   `g_E = 0.79748 − 0.37215·rE` doğrusunda ±%0.07 (rms 3.0·10⁻⁴);
   L130 −%1.37 ile dışarıda ve bunun sebebi ölçüldü (Q tümseğinin inen
   kolu, `dQ/drE = −7.55`, işaret dönüyor).
6. **Ön-kayıtların ikisi kendi eşik sayısını ıskaladı** (172d-Ö1 "> %7"
   dedi, %6.2 çıktı; 172b-Ö4′ "+2.41" dedi, +1.12 çıktı çünkü μ̂²'nin
   sürüklenmesini ihmal etmişti). İkisi de yazıldı, düzeltilmedi.
7. **`son`un mutlak seviyesi hiçbir yerde uydurulmadı:** bütün gerçek-gaz
   sayıları λ ailesinin eğrilerine göre okundu ve çapa `σ_X̃` (171'in
   konvansiyonu) idi.

---

## 7. NE KAZANILDI, NE KALDI

**Kazanılanlar.**
1. `g = 1 − Q/ρ` — iki alan için de ÖZDEŞ, makine hassasiyetinde
   kapanan, **Q = izdüşüm kusuru** cinsinden tek satırlık cebir. g_E ve
   g_X artık serbest sayı değil, türetilmiş okumalar.
2. `μ̂²_X ≡ 0` (kod X'i ortalar, E'yi ortalamaz) — 10⁻¹⁷'de mühür.
3. Üç eksenin **üç ayrı Q mekanizması**: pencere = aşırı uyum (d/N),
   λ ve kesim = koherent sızıntı (d/N orada TAM sabit).
4. **Q_E'nin tümseği λ_c = 1.0109'da** — 170'in inşa doyum eşiği.
5. **171 §T4'ün τ-eğimi borcu kapandı**: α_g'nin türevi, sıfır parametre,
   −149 … +54 %/τ yelpazesinde ±%7.
6. **c vadisinin yeri türetildi**: `dlogKALİB/dλ = dlogW_X/dλ`,
   λ = 0.6506 ↔ ölçülen 0.6487 (‰3); ve aynı yerde α tepe yapıyor.
7. **Gerçeğin dar adresi: E (η) kanalı** — beş büyüklük birden
   λ ≈ 0.74-0.82, X kanalı ve marjinaller λ ≈ 0.91-0.94. 162'nin
   kanal ayrımının çarpan dilindeki tam karşılığı.
8. θ'nın en iyi tek-değişkenli taşıyıcısı: `ρ_X` (üs 0.651).

**Kalanlar (173'ün kapıları).**
1. **`μ̂²_E` ve `Q_E`'yi κ cinsinden DOĞRUDAN ölçmek.** İkisi de
   `κ(ω_q)` ve `κ(ω_k−ω_l)` toplamlarıdır; alanları yeniden sentezleyip
   (3-4 dk/gaz) gerçek gazın −%23.8'lik DC kaçağı eksikliğini
   asal-frekans bazında haritalamak. **+3σ'nın ilk-ilke kaynağına en
   dar kapı budur.**
2. `θ`'nın kimliği: `ρ_X^{0.651}` bir taşıyıcı, kimlik değil; üssün
   (2/3?) türetimi ve T/4'teki 5-10σ artığın adresi.
3. `Q`'nun T-ekseni mekanizması: M1 ve M2 öldü, `B·r ≈ 0.59` koklaması
   sınanmadı.
4. 173 mühürleri: λ = 0.40 ve λ = 1.45 gazları (her biri ≈ 9 dk inşa +
   4 dk ölçüm) iki tahmin yolunu ayırır.
