# 174 — GERÇEĞİN İMZASI: ΔM(son)'un ilk-ilke kaynağı
### (4 Eylül 2026, Opus tayfası — KALEM_GERCEK_IMZA_04EYL2026.md)

**Soru.** Gerçek ζ gazı, sadakatli ikizinden (`Hkeskin`) **ΔM = +%4.97**
fazla (defter 3.3e−6'da kapalı: `g_E` +%2.05, `g_X²` +%0.32, `θ` +%2.53).
172-G3 fazlanın adresini bir KANALA indirmişti (E/η kanalı, λ_eş = 0.777).
174'ün AV'ı: **"+3σ fazlası, asal fazlarının kilidinin üçüncü-moment
faturasıdır"** cümlesini nicel mühürlemek — 162 arkını (faz kilidi) ile
165-173 arkını (çarpan defteri) birleştirmek.

Betikler (`174_configs/`), hepsi ön-mühürlü, hiçbiri yeni gaz inşa etmedi,
git'e dokunulmadı:

| betik | kapı | koşu | ham çıktı |
|---|---|---|---|
| `174a_onkayit.py` | **K2 ön-kaydı** (kural donduruldu) | 0.1 s | `174/ONKAYIT_K2.json` |
| `174b_K1_girisim.py` | **K1** — girişim oranları + bant defteri + 3. momentler, 8 gaz | 2.2-2.4 dk/gaz | `174/K1_<gaz>.json`, `174/K1c_<gaz>.npz` |
| `174c_K2_yuzlesme.py` | **K2** — dondurulmuş kuralın uygulanışı ve yüzleşme | 0.3 s | `174/K2.json` |
| `174d_K3_dc.py` | **K3** — μ̂²_E'nin çizgi-çizgi anatomisi | 1.3 dk/gaz | `174/K3_<gaz>.json` |
| `174e_K3b_kappa.py` | **K3b** — κ'nın kapalı yasası + gerçek−ikiz ayrıştırması (**ön-mühürsüz**) | 1.7 dk | `174/K3b.json` |
| `174f_ek_figur.py` | **E1/E2** — λ_eş atlası, θ'nın oran vekili (**ön-mühürsüz**) + figür | 25 s | `174/ATLAS.json`, `174_gercek_imza.png` |

Önbellek kökü `/private/tmp/claude-501/.../scratchpad/`; girdi
`155/eta_*.npz`, `155,167/z_*.npy`, `165/tayf_*.npz`, `172/G1-G3.json`,
`167/C_*.json`; çıktı `174/`. Toplam koşu ≈ 24 dk (8 gaz × K1, 3 gaz × K3).

---

> ## TEK CÜMLELİK HÜKÜM
>
> **162'nin "girişim oranı" ile 172'nin `π = K/V_O`'su AYNI NESNEDİR**
> (8 gazda ≤ %1.9 fark, `174c`) ve `g_E ≡ R_η/ρ_E` cebri köprüyü kurar.
> Ama köprü, AV'ı **mühürlemiyor, ikiye bölüyor**:
> **X kanalında** girişim oranı λ ile tekdüze akar ve gerçek gaz tam
> yerine oturur — `λ_eş(R_X) = 0.9026` ↔ defterin `0.9089`'u
> (**ön-kayıtlı tam isabet**, `174c` K2-G).
> **η kanalında ise λ_eş YOKTUR:** `R_η(λ)` bir TÜMSEKTİR (tavan 1.2745,
> λ = 0.70) ve gerçek gazın `R_η = 1.2883`'ü **bütün ailenin üstündedir**
> — yani gerçeğin η-girişimi bir λ kayması DEĞİL, ailenin dışında bir
> FAZLADIR (`174b`/`174c` K2-A, ÖLDÜ ama bilgilendirici).
> Parametresiz köprü ΔM'nin yalnız **%31.4'ünü** kapatıyor (ön-kayıtlı
> eşik %70) ⇒ **H-G3 MÜHÜRLENMEDİ** — ön-kayıt bunu önceden yazmıştı.
> Ayakta kalan: **H-G1'in yönü** (`R_η(gerçek) > R_η(ikiz)`, 1.2883 >
> 1.2650) ve **H-G2'nin anatomisi**: DC kaçağı τ = 0.5'te işaret
> değiştiren bir ORTA-NOKTA tarak rezonansıdır ve gerçek−ikiz açığının
> **%83'ü τ > 0.70'te**, yarısı en üst bantta (τ > 0.80) oturur; orada
> gerçek gazın faz uyumu ⟨cos Δφ⟩ 0.60, ikizinki 0.77.

---

## 0. ÖN-KAYIT (K2) — 4 Eylül 2026, 14:59:18 +03

`174a_onkayit.py` K1 koşmadan **önce** koştu; kendi sha256'sını
(`b643e6c73a2513347503d317b5db284d5d47f25cafa844a5d0921170c4cb9192`)
`174/ONKAYIT_K2.json`'a yazdı ve dosya bir daha yazılmaz.

**Dürüstlük beyanı (ön-kayıtta aynen duruyor).** Bu tayfa 172'nin bütün
çarpan defterini ve hedef sayıları (0.7768; +%2.05/+%2.53) okumuş
durumdaydı — **körlük iddiası yoktur**. Ön-kayıt edilen şey KURALDIR:
hangi ölçümden hangi öngörünün, hangi formülle ve hangi ölüm eşiğiyle
çıkacağı. Ön-kayıt anında **ölçülmemiş** olanlar: `R_η`/`R_Ĉ`'nin
Hkeskin ve L085 (ve bütün λ merdiveni) değerleri — 162 bunları hiç
ölçmedi ve 162'nin "keskin"i 152'nin gazıdır, `Hkeskin` değildir;
`R_X` hiçbir gazda; bant-bant defteri hiçbir gazda; üçüncü momentler
hiçbir gazda.

### 0.1 KÖPRÜNÜN CEBRİ — ölçümden önce türetildi (`174a` başlığı)

`165_cekirdek.Model165.alanlar`'da `E = Σ_q Re[hp_q e^{iω_q s}]`,
`hp_q = 2⟨e1 e^{−iω_q s}⟩` olduğundan, **hiçbir varsayım olmadan**

```
K := ⟨E·e1⟩ = Σ_q Re[ conj(hp_q)·⟨e1 e^{iω_q s}⟩ ] = Σ_q |hp_q|²/2
```

yani 172b'nin `π := K/V_O`'su **162'nin girişim oranının ta kendisidir**.
172b'nin özdeşliği `g_E = 1 − Q/ρ = π/ρ` olduğuna göre

> ### `g_E = R_η / ρ_E` ,  `R_η := Σ_q|c_q(η)|²/2 ÷ Var(η)`

Bu bir yasa değil **cebirdir**. Ampirik içerik iki yerdedir: (i) 162'nin
`R`'si başka bir çizgi kümesinde (τ ≤ 0.86, 3425 çizgi) ve başka bir
sitede (`m_n`, `η_n`) ölçülür, 172'nin `π`'si τ ≤ 0.95'te 8981 çizgide
`s_n = m_{n+1}`'de — aynı sayıyı vermeleri zorunlu değildir (K2-F);
(ii) `R`'nin λ merdiveni boyunca nasıl aktığı hiç ölçülmemişti (K1).

---

## 1. K1 — GİRİŞİM ORANLARI (162 makinesi, aynı pencere/taban)

`174b_K1_girisim.py`, `162_configs/162_cekirdek.Taban162`'yi **aynen**
kurar: taban 0.40, cap 4000, çizgi evreni `pk_m(e^{0.86L})`, `NITER = 0`,
`BLOK = 2000`, çıkarım `c_q(x) = 2⟨x_n e^{−i w_q m_n}⟩`. Hiçbir ölçüm
parçası kopyalanmadı. Sekiz gazın hepsinde `L = 12.029593242`,
`nline = 3425`, `m`-kimliği `0.0e+00`.

### 1a. GERÇEK GAZ — 162'nin sayıları BİREBİR yeniden çıktı

| nicelik | 162 (yayımlanan) | **174b (yeniden ölçüm)** |
|---|---|---|
| `Σ\|c_η\|²/2` | 0.069971 | **0.069971** |
| `Var(η)` | 0.054313 | **0.054313** |
| **`R_η`** | **1.2883** | **1.2883** |
| `Σ\|c_Ĉ\|²/2` | 0.074438 | **0.074438** |
| `Var(Ĉ)` | 0.072552 | **0.072552** |
| **`R_Ĉ`** | **1.0260** | **1.0260** |
| çizgi payı η / Ĉ | 0.487 / 0.918 | **0.4869 / 0.9175** |
| varyans kapanışı η / Ĉ | 1.801 / 1.109 | **1.8014 / 1.1085** |

Ö1 (ön-mühür ±%1) **✓ sıfır farkla**; Ö5 (çizgi payları ±0.005) **✓**.
KALEM'in "alıntıyla yetinilmez" şartı karşılandı.

### 1b. ÜÇ GAZ + BÜTÜN λ MERDİVENİ (`174b`)

KALEM üç gaz istedi (gerçek, ikiz, L085); 172e'nin `λ_eş` tersi
7-noktalı merdiven üzerinde tanımlı olduğu için merdivenin tamamı
koşuldu (5 gaz daha, +12 dk).

| gaz | λ | **`R_η`** | **`R_Ĉ`** | **`R_X`** | `P_η` | `Var(η)` | çizgi payı η |
|---|---|---|---|---|---|---|---|
| L050 | 0.50 | 1.2541 | 1.1073 | 1.1633 | 0.03718 | 0.02965 | 0.6044 |
| L060 | 0.60 | 1.2677 | 1.0830 | 1.1459 | 0.04888 | 0.03856 | 0.5390 |
| L070 | 0.70 | **1.2745** | 1.0567 | 1.1277 | 0.05958 | 0.04675 | 0.4956 |
| L085 | 0.85 | 1.2740 | 1.0205 | 1.0998 | 0.07382 | 0.05795 | 0.4590 |
| **Hkeskin** | 1.00 | **1.2650** | 0.9909 | 1.0726 | 0.08593 | 0.06793 | 0.4445 |
| L115 | 1.15 | 1.2510 | 0.9676 | 1.0470 | 0.09628 | 0.07696 | 0.4433 |
| L130 | 1.30 | 1.1739 | 0.9408 | 1.0124 | 0.10268 | 0.08747 | 0.4315 |
| **son (gerçek)** | — | **1.2883** | 1.0260 | 1.0902 | 0.06997 | 0.05431 | 0.4869 |

> **BULGU (K1-a) — Ö2 GEÇTİ, H-G1'İN YÖNÜ AYAKTA.**
> `R_η(gerçek) = 1.2883 > R_η(ikiz) = 1.2650`. Gerçek gazın η'sı,
> çizgilerinin uyumsuz toplamından ikizinkinden **daha çok** sessizdir.
>
> **BULGU (K1-b) — ama gerçek gaz AİLENİN DIŞINDADIR.**
> `R_η(λ)` tekdüze değil, bir **TÜMSEKTİR** (tavan **1.2745**, λ = 0.70)
> ve `R_η(son) = 1.2883` bu tavanın **%1.1 üstündedir**. Yani η-kanalı
> girişim oranının bir `λ_eş`'i **YOKTUR**. Aynı şey defterin kendi
> `π_E`'sinde de doğrudur: merdiven 1.2650 / 1.2804 / 1.2890 / **1.2908**
> / 1.2837 / 1.2711 / 1.1958, `π_E(son) = 1.30385` — **aralık dışı**
> (`174c`, "defter denetimi" satırı).
>
> **BULGU (K1-c) — X ve Ĉ kanalları TEKDÜZE.** `R_X` 1.1633 → 1.0124 ve
> `R_Ĉ` 1.1073 → 0.9408, ikisi de λ ile tekdüze azalıyor ve gerçek gaz
> **içeride**. 162'nin "Ĉ kanalında girişim YOK (1.026 ≈ 1)" cümlesi bir
> yapı değil bir **rastlantıdır**: `R_Ĉ(λ)` ailesi tam da λ ≈ 0.83'te
> 1'i kesiyor ve gerçek gaz oraya düşüyor. Ailede yayılım **%17.7**.

### 1c. Bant defteri ÖZDEŞ kapanıyor (`174b`, Ö4)

`c_q` ham izdüşüm olduğu için `⟨x·Re(c_q e^{iw_q m})⟩ ≡ |c_q|²/2`;
dolayısıyla **fit içermeyen özdeş defter**

```
⟨x²⟩ = Σ_b P_b + ⟨x·x_artık⟩ ,   P_b := Σ_{q∈b}|c_q|²/2
```

Bağıl kalıntı sekiz gazda ve üç kanalda **≤ 4.1e−13** (gerçek gazda
η için **2.6e−16**). Ö4 eşiği 1e−12 ✓.

### 1d. BANT DEFTERİ — η kanalı (gerçek gaz / sadakatli ikiz)

| τ bandı | çizgi | `P_b`(son) | `V_b`(son) | `P/V` | `P_b`(ikiz) | `V_b`(ikiz) | `P/V` |
|---|---|---|---|---|---|---|---|
| 0.00–0.40 | 41 | 0.000000 | 0.000000 | 0.994 | 0.000000 | 0.000000 | 0.998 |
| 0.40–0.45 | 21 | 0.016518 | 0.017802 | 0.928 | 0.018476 | 0.019906 | 0.928 |
| 0.45–0.50 | 37 | 0.015596 | 0.017134 | 0.910 | 0.017939 | 0.019699 | 0.911 |
| 0.50–0.55 | 56 | 0.012320 | 0.013721 | 0.898 | 0.014677 | 0.016346 | 0.898 |
| 0.55–0.60 | 90 | 0.009540 | 0.010482 | 0.910 | 0.011865 | 0.013063 | 0.908 |
| 0.60–0.65 | 157 | 0.007017 | 0.007535 | 0.931 | 0.009208 | 0.009931 | 0.927 |
| 0.65–0.70 | 254 | 0.004255 | 0.004507 | 0.944 | 0.006008 | 0.006393 | 0.940 |
| 0.70–0.75 | 432 | 0.002595 | 0.002730 | 0.951 | 0.003977 | 0.004204 | 0.946 |
| 0.75–0.80 | 735 | 0.001322 | 0.001376 | 0.960 | 0.002266 | 0.002372 | 0.955 |
| 0.80–0.86 | 1602 | 0.000808 | 0.000837 | 0.966 | 0.001516 | 0.001586 | 0.956 |
| **iç-bant Σ V_b** | | **0.076124** | | | **0.093500** | | |
| **bantlar-arası Σ Kov** | | **+0.037374** | | | **+0.048168** | | |
| **Var(η_çizgi)** | | **0.113498** | | | **0.141668** | | |

> **BULGU (K1-d) — "yıkıcı girişim" ÇİZGİLER ARASINDA DEĞİL.**
> Her bantta `P_b/V_b < 1` (0.90–0.97): bant içindeki çizgiler **yapıcı**
> toplanıyor. Bantlar arası kovaryans da **artı** (+0.0374). Çizgi
> alanının kendisi uyumsuz toplamdan **%62 daha gürültülüdür**
> (`Var(η_çizgi)/P = 0.1135/0.0700 = 1.621`).
> `R_η > 1`'in gerçek anlamı şudur: `P ≡ Kov(η_çizgi, η)` özdeş olduğu
> için `R > 1` ⟺ **çizgi alanı η'yı AŞIYOR** ve artık onu geri kesiyor:
> `Kov(η_çizgi, η_artık) = P − Var(η_çizgi) = −0.04353` (ikizde
> −0.05574). Yani 162'nin "η, çizgilerinin uyumsuz toplamından %29 daha
> sessizdir" cümlesi doğrudur, ama mekanizması **çizgi ↔ artık yıkıcı
> girişimidir**, çizgi ↔ çizgi değil. *(162'nin VL/VC tanısıyla
> çelişmez: faz karıştırması ikisini birden yok ediyordu.)*

Aynı defter Ĉ ve X̃ kanallarında da koşuldu (`174/K1_*.json`). Öne çıkan:
Ĉ'nin gücünün **%91'i τ ≤ 0.40 bandındadır** (η regresyonunun η'dan
sildiği bölge Ĉ'de duruyor); X̃'de bant-içi oranlar 0.92–0.99, bantlar
arası kovaryans +0.0240 (gerçek) / +0.0252 (ikiz).

### 1e. ÜÇÜNCÜ MOMENTLER (162'de ölçülmemişti)

Rastgele fazlı çizgi alanının bütün üçüncü momentleri ÖZDEŞ SIFIRDIR
(Rice); dolayısıyla bunlar saf faz-kilidi ölçüleridir.

| gaz | `m3 = ⟨e1x1²⟩/σ_{e1}σ_{x1}²` | `m3_çizgi` (model) | `m3/m3_çizgi` | `skew(e1)` | `skew(η_çizgi)` | `skew(ds)` |
|---|---|---|---|---|---|---|
| L050 | +0.15025 | −0.03446 | −4.3605 | −0.1554 | −0.1975 | +0.2032 |
| L060 | +0.18492 | −0.04872 | −3.7954 | −0.1567 | −0.2175 | +0.2443 |
| L070 | +0.20665 | −0.06077 | −3.4007 | −0.1629 | −0.2300 | +0.2691 |
| L085 | +0.22751 | −0.07327 | −3.1050 | −0.1764 | −0.2420 | +0.2960 |
| Hkeskin | +0.24133 | −0.08022 | −3.0083 | −0.1910 | −0.2501 | +0.3185 |
| L115 | +0.25163 | −0.08316 | −3.0257 | −0.2048 | −0.2563 | +0.3396 |
| L130 | +0.25007 | −0.08083 | −3.0938 | −0.1322 | −0.2532 | +0.4148 |
| **son** | **+0.27675** | −0.06575 | **−4.2093** | −0.1513 | −0.2216 | **+0.4681** |

Gerçek gaz `m3`, `skew(x1)` ve `skew(ds)`'te **ailenin üstünde**, `skew(e1)`'de
ailenin (tekdüze kolunun) **dışında**. `skew(ds) = 0.4681` λ = 1.30'unkini
(0.4148) bile aşıyor.

---

## 2. K2 — DONDURULMUŞ KURALIN UYGULANMASI (`174c`)

**Denetim D1.** 172e'nin `lam_es` ters çevirmesi birebir kopyalandı;
172e'nin yayımlanmış λ_eş tablosunu `G1/G2.json`'dan yeniden üretiyor,
**maks fark 0.00e+00**.

### K2-F — KÖPRÜ KİMLİĞİ: 162'nin `R`'si = 172'nin `π`'si mi? **✓ GEÇTİ**

| gaz | `R_η` (162 makinesi, τ≤0.86, site `m_n`) | `π_E` (172 defteri, τ≤0.95, site `s_n`) | fark | `R_X` | `π_X` | fark |
|---|---|---|---|---|---|---|
| son | 1.28829 | 1.30385 | −1.19% | 1.09022 | 1.11684 | −2.38% |
| Hkeskin | 1.26504 | 1.28375 | −1.46% | 1.07256 | 1.10024 | −2.52% |
| L050 | 1.25410 | 1.26496 | −0.86% | 1.16334 | 1.21256 | −4.06% |
| L060 | 1.26768 | 1.28035 | −0.99% | 1.14587 | 1.18812 | −3.56% |
| L070 | 1.27449 | 1.28896 | −1.12% | 1.12775 | 1.16474 | −3.18% |
| L085 | 1.27397 | 1.29081 | −1.30% | 1.09975 | 1.13106 | −2.77% |
| L115 | 1.25103 | 1.27108 | −1.58% | 1.04705 | 1.07263 | −2.38% |
| L130 | 1.17388 | 1.19577 | −1.83% | 1.01243 | 1.04041 | −2.69% |

Ön-kayıt: ≤ %5. Ölçülen maks **%1.83** (E), **%4.06** (X) ⇒ **✓ GEÇTİ**.

> **HÜKÜM (K2-F).** 162 ile 172 **aynı sayıyı ölçüyor**. İki ark
> arasındaki köprü kurulmuştur ve farkı yalnız çizgi kümesi (3425 ↔
> 8981) ile site (m_n ↔ m_{n+1}) açıklar: sistematik olarak %1-2
> (E), %2-4 (X), işareti hep aynı.

### K2-A — λ_eş ÖNGÖRÜSÜ (E kanalı): **ÖLDÜ**

Birincil (7-nokta, 172e ile aynı ters çevirme): **λ_eş^ön(E) YOK** —
`R_η(son) = 1.2883` ailenin tavanının (1.2745) üstünde.
Yedek (2-nokta, L085↔Hkeskin log-doğrusal): **0.6114**.
Hedef 0.7768, ölçüt ≤ 0.05. **Her iki yoldan da ÖLDÜ** (kurtarma yok).

> **HÜKÜM (K2-A).** Ölüm bilgilendiricidir: 172e'nin `λ_eş(E) = 0.777`'si
> **tekdüze** E-büyüklüklerinden (rE, ρ_E, g_E, μ̂², Q_E) çıkar; onların
> altındaki **oran** `π = R_η` tekdüze değildir ve gerçek gaz onun
> aralığının dışındadır. Yani "gerçek gaz E kanalında λ ≈ 0.78 gibi
> davranıyor" ifadesi **güçlerin** dilinde doğrudur
> (`λ_eş(Var η) = 0.8012`, `λ_eş(P_η) = 0.8093`, `λ_eş(çizgi payı η)
> = 0.7359`; `174f` E1) ama **girişim oranının** dilinde YANLIŞTIR:
> orada gerçek gaz λ ailesinin DIŞINDADIR. Fazlanın adresi bir λ
> kayması değil, **ailenin taşımadığı bir fazladır.**

### K2-G (= K4) — X KANALI SAĞLAMASI: **TAM İSABET**

`λ_eş^ön(X) = 0.9026` (7-nokta) / 0.9022 (2-nokta) ↔ hedef **0.9089**;
fark **−0.0063**, ölçüt ≤ 0.05 ⇒ **TAM İSABET**.

> **HÜKÜM (K2-G).** KALEM'in bonus sorusu (**"aynı eşleme λ_eş(X) =
> 0.909'u da veriyor mu?"**) **EVET**. X (ΔĈ) kanalında girişim oranı
> λ ile tekdüze akar, gerçek gaz aralığın içindedir ve tam defterin
> söylediği yere düşer. **Asimetri budur ve şimdi mekanizmalıdır:**
> X kanalı bir λ-gazı gibi davranır, η kanalı davranmaz.

*(Yan sonuç, aynı kuralla: `λ_eş(R_Ĉ) = 0.8272` — 162'nin "Ĉ'de girişim
yok" okuması bir aile-rastlantısıdır, §1b.)*

### K2-B — g_E PAYI (parametresiz): **KISMİ**

```
Δlog g_E^ön = log R_η(son) − log R_η(λ_çapa=0.9363)
            = log(1.28829/1.26883) = +0.01522   ⇒  +1.53%
```
Hedef +2.05%; fark **−0.52 puan**; ölçüt TAM ≤ 0.5, KISMİ ≤ 1.0 ⇒
**KISMİ** (tam isabetin 0.02 puan dışında).

ÖZDEŞ ayrıştırma (172 defteriyle, `174c`):
```
Δlog g_E = Δlog π_E − Δlog ρ_E :  +0.02029 = +0.01320 − (−0.00710)
```
⇒ **g_E artığının %65'i GİRİŞİM (π), %35'i MODEL GÜCÜ (ρ)**.
(ρ_E düzeltmesi eklenirse öngörü +2.26% olur — hedefin 0.21 puan
üstünde — ama bu artık parametresiz değildir ve hükme girmez.)

### K2-C — θ PAYI: **ÖLDÜ, ve ön-kaydın PREMİSİ de yanlış çıktı**

Ön-kayıt "R_Ĉ üç gazda ±%2 içinde aynı" demişti; ölçülen yayılım
**%17.7** (1.1073 … 0.9408, tekdüze). **Premise yanlış.** Buna rağmen
sonuç değişmiyor: `Δθ^ön = 0.00%` ↔ hedef **+2.53%** ⇒ **ÖLDÜ**.
(`λ_eş(R_Ĉ) = 0.8272`, `λ_eş(θ) = 0.6916` — Ĉ-kanalı girişimi θ'nın
yerini vermiyor.)

### K2-E — ÜÇÜNCÜ MOMENT KÖPRÜSÜ: **ÖLDÜ, ve TERS YÖNDE**

Ön-kayıt tek yönlüydü: `λ_eş(m3) < λ_eş(R_η)`. Ölçüm:
`m3(son) = +0.27675` **ailenin üstünde** (tavan 0.25163, λ = 1.15) —
7-noktalı ters çevirmede λ_eş **YOK**; 2-noktalı yedek **1.3481**.
Yani gerçek gaz üçüncü momentte λ ≈ 1.35 gibi, θ'nın istediği 0.69'un
tam **tersi** yönde. **K2-E ÖLDÜ.**

### K2-D — ΔM KAPANIŞI ve H-G3: **MÜHÜRLENMEDİ**

```
Δlog M(ölçülen)         = +0.04850  (+4.97%)
köprünün kapattığı      = +0.01522  ⇒ %31.4      (eşik %70)
```
Ön-kayıtlı beklenti %41 idi; ölçülen **%31.4**. Defterin kendi π_E'siyle
%27.2, ρ düzeltmeli sürümle %46.1 — **hiçbiri %70'e yaklaşmıyor.**

> **HÜKÜM (K2-D / H-G3).** **H-G3 (taç) MÜHÜRLENMEDİ.** Ön-kayıt bunu
> önceden yazmıştı ve ölçüm onu doğruladı. `ΔM`'nin girişim eksikliğine
> parametresiz olarak bağlanabilen payı **%27–31**; geri kalanın tam
> defteri (log paylar, toplam %100): girişim `π_E` **%27.2**, model gücü
> `ρ_E` **%14.6**, `g_X²` **%6.6**, **`θ` %51.5** — ve θ'nın taşıyıcısı
> bu köprüde **yoktur**.

### K2 karnesi

| madde | ön-kayıt | ölçüm | hüküm |
|---|---|---|---|
| **K2-F** köprü kimliği | \|R/π−1\| ≤ %5 | %1.83 (E), %4.06 (X) | **✓ GEÇTİ** |
| **K2-G** λ_eş(X) | 0.9089 ± 0.05 | **0.9026** | **✓ TAM İSABET** |
| **K2-B** g_E payı | +2.05% ± 0.5 | **+1.53%** | **KISMİ** (0.52 pt) |
| **K2-A** λ_eş(E) | 0.7768 ± 0.05 | **YOK / 0.6114** | **✗ ÖLDÜ** |
| **K2-C** θ payı | 0 ↔ +2.53% | 0 | **✗ ÖLDÜ** (premis de yanlış) |
| **K2-E** λ_eş(m3) < λ_eş(R_η) | yön | **ters** (1.348 > 0.611) | **✗ ÖLDÜ** |
| **K2-D** ΔM kapanışı | ≥%70 mühür | **%31.4** | **✗ MÜHÜR YOK** |
| **Ö1** (174b) R_η(son) = 162 | ±%1 | 0.00% | ✓ |
| **Ö2** (174b) R_η(ikiz) < R_η(son) | yön | 1.2650 < 1.2883 | ✓ |
| **Ö3** (174b) R_Ĉ 1.00-1.05, ±%2 | — | 0.94-1.11, %17.7 | ✗ |
| **Ö4** (174b) defter özdeş | ≤1e−12 | ≤4.1e−13 | ✓ |
| **Ö5** (174b) çizgi payları | ±0.005 | 0.0001 | ✓ |

---

## 3. K3 — μ̂²_E'NİN ANATOMİSİ: DC KAÇAĞI (`174d`)

**Özdeşlik.** `ort(E) = Σ_q Re[hp_q κ(ω_q)]`, `κ(ν) = ⟨e^{iν s_n}⟩`,
`μ̂²_E = ort(E)²/Var(e1)`. Çizgi-çizgi ayrışım:
`ξ_q := Re[hp_q κ(ω_q)] = |hp_q||κ(ω_q)| cos Δφ_q`, üç çarpan: **çizgi
genliği**, **tarak rezonansı**, **bağıl faz**.

**Ö2 denetimi (bağımsız yol).** Σξ_q'dan hesaplanan μ̂², 172/G1.json'un
μ̂²_E'sini **1.3e−13** (son) ve **1.8e−13** (ikiz) bağıl farkla veriyor.
Bu, 165 → 172 → 174 zincirinin kapandığının mührüdür.

| gaz | `ort(E)` | `μ̂²_E` | `Var(e1)` |
|---|---|---|---|
| son | **+0.046900** | 0.040498 | 0.054313 |
| Hkeskin | **+0.062139** | 0.056843 | 0.067929 |
| L085 | **+0.052819** | 0.048145 | 0.057947 |

`ort(E)` oranı gerçek/ikiz = **0.7548** ⇒ μ̂² oranı **0.7125**
(172'nin σ-çapasındaki −%23.8'i, λ=1 ikizine karşı −%28.7).

### 3a. ÖN-MÜHÜRLÜ ÇIPLAK LİMİT — İKİ MADDESİ ÖLDÜ

Türetim (`174d` başlığı, koşudan önce): `dN = (N̄'+S')dt`,
`S = −Σ A_q sin(ω_q t)` ⇒ `κ(ω_Q) ≈ −A_Qω_Q/(2N̄') = −πA_Qτ_Q`
(reel, negatif) ⇒ `ξ_q^çıplak = −2πA_q²τ_q sin(πτ_q) < 0`.

| ön-mühür | beklenti | ölçüm | hüküm |
|---|---|---|---|
| **Ö1** `ort(E) < 0` | negatif | **+0.0469 / +0.0621 / +0.0528** | **✗ ÖLDÜ** |
| **Ö2** μ̂²(ξ) = μ̂²(G1) | ‰5 | **1.3e−13** | ✓ |
| **Ö3** en büyük 10 çizgi ≥ %50 | ≥%50 | **%5.2 / %4.8** | **✗ ÖLDÜ** |
| **Ö4** \|κ\|/\|κ^çıplak\| ∈ [0.3, 3] | — | 0.045 … 0.97 (11 bandın 7'si içeride) | **kısmi ıska** |
| **Ö5** \|ort E\|(gerçek) < (ikiz) | evet | 0.0469 < 0.0621 | ✓ |

### 3b. NEDEN — ve (−) işaretin gerçek kökeni

`174d`'nin bant defteri (gerçek gaz):

| τ bandı | çizgi | `Σξ_q` | `Σ\|hp\|\|κ\|` | `⟨cos Δφ⟩` |
|---|---|---|---|---|
| 0.40–0.50 | 58 | **−4.935e−03** | 5.815e−03 | **−0.8486** |
| 0.50–0.60 | 146 | +9.509e−03 | 9.511e−03 | **+0.9998** |
| 0.60–0.70 | 411 | +1.720e−02 | 1.722e−02 | +0.9987 |
| 0.70–0.80 | 1167 | +1.538e−02 | 1.557e−02 | +0.9880 |
| 0.80–0.95 | 7158 | +9.740e−03 | 1.612e−02 | **+0.6043** |

> **BULGU (K3) — DC KAÇAĞI τ = 0.5'TE İŞARET DEĞİŞTİRİR.**
> Çıplak türetim `κ`'yı **sıfır** tarağı için yazmıştı; oysa ölçüm sitesi
> `s_n = mid_{n+1}` — **ORTA NOKTA** tarağıdır. Yarım-gap kayması κ'ya
> `cos(πτ)` çarpanı takar (`B_q = b_q cos πτ_q` ile aynı çarpan) ve
> `cos(πτ)` τ = 0.5'te işaret değiştirir. Ölçüm bunu birebir gösteriyor:
> `⟨cos Δφ⟩` τ < 0.5'te **−0.85**, τ > 0.5'te **+1.00**; `|κ|`'nın
> medyanı (0.45, 0.50] bandında **1.17e−03**'e çöküyor (komşu bantlarda
> 5.1e−03 ve 3.1e−03) — **cos(πτ) çukuru ölçüldü**.
> Kaçağın toplamı, çizgi sayısı τ ile patladığı için τ > 0.5 tarafından
> yönetilir; bu yüzden `ort(E) > 0`. **Ö1 bu yüzden öldü ve ölümün
> sebebi tam olarak budur.** Ö3 de aynı sebeple öldü: kaçak birkaç büyük
> çizgide değil, **binlerce çözülmemiş yüksek-τ çizgisindedir.**

### 3c. GERÇEK − İKİZ FARKI NEREDE OTURUYOR (`174e`)

| τ bandı | `Σξ`(gerçek) | `Σξ`(ikiz) | Δ | Δ payı | `Σ\|hp\|` oranı | `Σ\|κ\|` oranı | `⟨cosΔφ⟩` oranı |
|---|---|---|---|---|---|---|---|
| 0.40–0.50 | −4.9350e−03 | −6.3898e−03 | +1.4548e−03 | **−9.5%** | 0.9386 | 0.8891 | 0.9302 |
| 0.50–0.60 | +9.5090e−03 | +1.0125e−02 | −6.1643e−04 | +4.0% | 0.9057 | 1.0312 | 1.0002 |
| 0.60–0.70 | +1.7202e−02 | +2.0580e−02 | −3.3780e−03 | +22.2% | 0.8555 | 0.9758 | 0.9999 |
| 0.70–0.80 | +1.5384e−02 | +2.0628e−02 | −5.2437e−03 | **+34.4%** | 0.7820 | 0.9547 | 0.9964 |
| 0.80–0.95 | +9.7396e−03 | +1.7196e−02 | −7.4560e−03 | **+48.9%** | 0.7745 | 0.9415 | **0.7825** |
| **TOPLAM** | **+0.046900** | **+0.062139** | **−0.015239** | 100% | | | |

> **BULGU (K3-b) — AÇIK, MERDİVENİN TEPESİNDEDİR.**
> Gerçek gazın DC kaçağı eksikliğinin **%83.3'ü τ > 0.70'te**, **%48.9'u
> tek başına τ > 0.80'de** oturuyor. En üst bantta iki mekanizma birden
> çalışıyor: çizgi genlikleri `Σ|hp|` ikizin **%77.5'i** ve faz uyumu
> `⟨cos Δφ⟩` ikizin **%78.3'ü** (0.6043 ↔ 0.7723). Alt bantta (τ < 0.5)
> fark **ters işaretlidir** (−%9.5): orada gerçek gaz ikizden daha
> "kaçaklı"dır.
> Yani 172'nin `μ̂²_E = −%23.8`'i, **çözülmemiş yüksek-τ çizgilerinin
> gerçek ζ tarağında hem daha zayıf hem daha dağınık fazlı olmasıdır** —
> tam da 165'in `τ_c ≈ 0.58` faz-seğirme sınırının ötesindeki bölge.

### 3d. `κ`'NIN KAPALI YASASI — **ön-mühürsüz ikinci tur** (`174e`)

Ölümden sonra sorulan soru: doğru çıplak limit ne? İki düzeltme,
**ikisi de yeni parametre içermiyor**:

```
κ(ω_q) ≈ −π A_q τ_q · cos(πτ_q) · W_pos(τ_q) ,
   A_q = λ a_q w_q  (gazın İNŞA genliği),
   W_pos(τ) = ⟨e^{2πiτ Ĉ_n}⟩  (166/167'nin ZATEN ölçtüğü nicelik)
```

| τ | 0.40–0.45 | 0.45–0.50 | 0.50–0.55 | 0.55–0.60 | 0.60–0.65 | 0.65–0.70 | 0.70–0.75 | 0.75–0.80 | 0.80–0.85 | 0.85–0.90 | 0.90–0.95 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `\|κ\|/\|κ̂\|` (son) | **0.966** | **0.807** | 2.345 | 1.902 | 1.906 | 2.011 | 2.167 | 2.383 | 2.694 | 3.050 | 3.784 |
| `\|κ\|/\|κ̂\|` (ikiz) | **1.109** | **0.969** | 2.192 | 1.922 | 1.981 | 2.130 | 2.330 | 2.599 | 2.979 | 3.414 | 4.241 |
| işaret uyumu (son) | 0.905 | 0.514 | 1.000 | 1.000 | 1.000 | 1.000 | 0.995 | 0.997 | 0.983 | 0.861 | 0.603 |

> **BULGU (K3-c, ön-mühürsüz).** Yasa **τ ≤ 0.50'de ±%20 içinde
> tutuyor** (0.81–1.11, üç gazda) — yani orada `κ` **saf birinci-mertebe
> doğrusal tepkidir** ve işaret yasası (`sign κ = −sign cos πτ`)
> çizgilerin %90'ında doğrudur. τ > 0.50'de ölçülen `|κ|` yasayı
> **×1.9'dan ×4.2'ye** kadar aşıyor ve fazla τ ile tekdüze büyüyor.
> Bu fazlanın adayı belli: `κ(ω_Q)`'ya yalnız `ω_Q` çizgisi değil,
> **çarpımsal çözümler** (`log q₁ + log q₂ = log(q₁q₂) = ω_Q`) de
> katkı verir ve q büyüdükçe çarpanlara ayırma sayısı patlar. Bu bir
> iddia değil, **ölçülmüş bir fazlanın adı konmuş adayıdır** — sınavı
> 175'e borç (alt-merdiven sayımıyla doğrudan sayılabilir).

---

## 4. K4 — X KANALININ SAĞLAMASI

KALEM'in bonus sorusu §K2-G'de yanıtlandı: **EVET**, aynı eşleme
`λ_eş(R_X) = 0.9026` veriyor (hedef 0.9089, fark 0.0063). Ek olarak:

* `R_X` merdiveni tekdüze (1.1633 → 1.0124), gerçek gaz **içeride**;
* `λ_eş(çizgi payı X̃) = 0.8757`, `λ_eş(σ_ds) = 0.9014`
  (172e'nin σ_ds'i için verdiği **0.9014** ile birebir),
  `λ_eş(σ_Ĉ) = 0.9654` (172e: 0.9684, tanım farkı: `174b` Ĉ'yi
  trendsizleştirmiyor);
* X kanalında bant-içi oranlar `P_b/V_b` 0.92–0.99 — η'daki 0.90–0.97'ye
  çok yakın; farkı yaratan bant sayısı değil, **artık payıdır**
  (çizgi payı X̃ 0.68 ↔ η 0.49).

> **HÜKÜM (K4).** X kanalı, girişim oranı dilinde de bir λ-gazıdır ve
> gerçek gaz oraya **marjinallerin biraz altına** (0.903 ↔ 0.935)
> oturur. 172-G3'ün "X kanalı marjinale yakın" hükmü **girişim
> ölçüsüyle bağımsız olarak doğrulandı.**

---

## 5. EK TANILAR — **ön-mühürsüz ikinci tur** (`174f`)

### 5a. λ_eş ATLASI (K1'in bütün nicelikleri, 7-nokta merdiveni)

| nicelik | son | λ_eş | tekdüze | not |
|---|---|---|---|---|
| `R_η` (girişim, η) | 1.28829 | **YOK** | ∩ | **aralık dışı (üstünde)** |
| `R_Ĉ` (girişim, Ĉ) | 1.02599 | 0.8272 | ↓ | |
| `R_X` (girişim, X̃) | 1.09022 | **0.9026** | ↓ | hedef 0.9089 |
| `P_η = Σ\|c\|²/2` | 0.06997 | 0.8093 | ↑ | |
| `Var(η)` | 0.05431 | 0.8012 | ↑ | |
| `Var(η_çizgi)` | 0.11350 | 0.8047 | ↑ | |
| çizgi payı η | 0.48688 | 0.7359 | ↓ | |
| çizgi payı X̃ | 0.68038 | 0.8757 | ↓ | |
| bantlar-arası Kov(η) | 0.03737 | 0.7938 | ↑ | |
| iç-bant Σ V_b(η) | 0.07612 | 0.8111 | ↑ | |
| `Kov(η_çiz, η_artık)` | −0.04353 | 0.7986 | ↓ | |
| `m3 = ⟨e1x1²⟩/σσ²` | 0.27675 | **YOK** | ↑ | **aralık dışı (üstünde)** |
| `m3_çizgi` (model) | −0.06575 | 0.7597 | ↓ | |
| **`m3/m3_çizgi`** | −4.20933 | **0.5267** | ∩ | θ'nın oran vekili |
| `skew(e1)` | −0.15129 | **YOK** | ↓ | aralık dışı |
| `skew(η_çizgi)` | −0.22159 | 0.6328 | ↓ | |
| `skew(x1)` | 0.28174 | **YOK** | ↑ | **aralık dışı (üstünde)** |
| `skew(ds)` | 0.46806 | **YOK** | ↑ | **aralık dışı (üstünde)** |
| `σ_ds` | 0.40919 | 0.9014 | ↑ | 172e ile birebir |
| `σ_Ĉ` | 0.26935 | 0.9654 | ↑ | 172e 0.9684 |

> **BULGU (E1).** Gerçek gazın λ ailesine göre konumu **niceliğe göre
> üç kümeye ayrılıyor**: (i) **marjinaller ve X kanalı** λ ≈ 0.88–0.97
> (aile içi, düzgün); (ii) **η kanalının güçleri** λ ≈ 0.74–0.81 (aile
> içi); (iii) **oranlar ve üçüncü momentler** (`R_η`, `m3`, `skew(x1)`,
> `skew(ds)`, `skew(e1)`) — **AİLE DIŞI**. Yani gerçek gazın imzası,
> hiçbir λ'nın veremeyeceği bir bölgede duruyor ve orası tam olarak
> **normalize edilmiş oranlar ve üçüncü momentlerdir** — 162'nin
> "aranacak nesne bispektrumdur" hükmünün 174'teki karşılığı.

### 5b. θ'nın ORAN vekili — K2-E'nin ölümünden sonra

θ tanımı gereği bir ORANDIR (`KALİB = ölçülen 3. moment / model 3.
moment`). K2-E'nin ön-kayıtlı vekili yalnız **ölçüleni** kullanıyordu ve
öldü. Oran vekili `m3/m3_çizgi`:

```
merdiven: −4.3605 −3.7954 −3.4007 −3.1050 −3.0083 −3.0257 −3.0938
son     : −4.2093        ⇒  λ_eş = 0.5267      (θ'nın hedefi 0.6916)
```

θ merdiveni: 0.9636 0.9323 0.9125 0.8959 0.8884 0.8856 0.8668;
θ(son) = 0.9142.

> **NOT (E2, ön-mühürsüz).** Oran vekili θ ile **aynı yönde** (λ çapanın
> çok altında) ama nicel ıska (0.53 ↔ 0.69). Bu hiçbir ölümü kurtarmaz;
> yalnız 175'e bir yön verir: **θ'yı arayan doğru nesne "ölçülen üçüncü
> moment" değil, "ölçülen/model üçüncü moment ORANI"dır.**

---

## 6. FİGÜR

`174_gercek_imza.png` (`174f_ek_figur.py`) — dört panel:
**(a)** `R_η(λ)` tümseği ve gerçek gazın tavanın üstündeki çizgisi
(K2-A'nın ölümü); **(b)** `R_X(λ)` ve `R_Ĉ(λ)`'nın tekdüze inişi,
gerçek gazın kesişimi ve `λ_eş = 0.9026 ↔ 0.9089` (K2-G'nin isabeti);
**(c)** DC kaçağının bant defteri (gerçek ↔ ikiz), τ = 0.5 işaret dönüşü
ve açığın τ > 0.70'te toplanması; **(d)** tarak rezonansının çıplak
doğrusal tepkiye oranı, τ < 0.5'te `cos πτ` çukuru, τ > 0.5'te
çarpımsal kanalın fazlası.

---

## 7. DENETİM

| | sınav | sonuç |
|---|---|---|
| **D1** | 172e'nin `lam_es`'i birebir mi? | 172e'nin λ_eş tablosunu G1/G2'den yeniden üretiyor, **maks fark 0.00e+00** (`174c`) |
| **D2** | 162 makinesi bit düzeyinde mi? | `L`, `nline = 3425`, `m`-kimliği `0.0e+00` sekiz gazda özdeş; gerçek gazda 162'nin altı sayısı **sıfır farkla** (`174b`) |
| **D3** | bant defteri özdeş mi? | 8 gaz × 3 kanal, bağıl kalıntı **≤ 4.1e−13** (gerçek gazda 2.6e−16) |
| **D4** | K3'ün μ̂²'si 172'nin defteriyle uyuşuyor mu? | `Σ_q ξ_q`'dan `μ̂² = ort(E)²/Var(e1)`, G1.json ile bağıl fark **1.3e−13** (son) / **1.8e−13** (Hkeskin) / **5.2e−14** (L085) |
| **D5** | ön-kayıt sonradan değişti mi? | `ONKAYIT_K2.json` yalnız bir kez yazıldı; betik varlık kontrolüyle üzerine yazmayı reddediyor; sha256 raporda |
| **D6** | git | **dokunulmadı** |

### 7a. Dürüstlük notları (kayda geçiyor)

* **`174e` iki kez koştu.** İlk koşuda κ-oranı bant ORTALAMASI olarak
  yazılmıştı; `cos(πτ)` çarpanı τ = 0.5'te sıfırlandığı için ortalama
  patlıyordu (bant değerleri 2.7 … 6.4). İkinci koşuda istatistik ince
  bantlarda MEDYANA çevrildi. Bu bir tanı düzeltmesidir, bir kurtarma
  değil: `174e` zaten **ön-mühürsüz ikinci tur** olarak etiketlidir ve
  `174d`'nin ölen ön-mühürleri (Ö1, Ö3) değişmedi.
* **172'nin defteri yeniden ölçülmedi**, önbellekten okundu
  (`172/G1.json`, `172/G2.json`, `167/C_*.json`). Yeniden ölçülen tek
  şey 162'nin girişim oranıdır (K1) ve o **sıfır farkla** çıktı.
* **`174b`'nin `σ_Ĉ`'si trendsizleştirilmemiştir** (0.26935), 172e/167
  ise kübik trendsiz `σ_Ĉ` kullanır (0.27303). λ_eş'lerdeki 0.003'lük
  fark (0.9654 ↔ 0.9684) bundandır.
* **Hedef sayılar ön-kayıt anında biliniyordu**; ön-kayıt körlük değil
  KURAL sabitlemesidir (§0). Ölçülmemiş olanlar §0'da tek tek sayılıdır
  ve hepsi bu koşuda ilk kez ölçüldü.

---

# HÜKÜM

## (i) KÖPRÜ KURULDU — 162'nin oranı 172'nin π'sidir

> `K = ⟨E·e1⟩ ≡ Σ_q|hp_q|²/2` özdeşliği (türetim `174a`'da, ölçümden
> önce) 162 ile 172'yi birleştiriyor ve ölçüm **8 gazda ≤ %1.9** (E),
> **≤ %4.1** (X) farkla doğruluyor. Bundan sonra 162'nin "yıkıcı
> girişim"i ile 172'nin "projeksiyon kusuru Q"su ayrı iki nesne değildir:
> **`g_E ≡ R_η/ρ_E`.**

## (ii) H-G1 YARIM AYAKTA — yönü doğru, ama "λ kayması" değil

> Yön tuttu: `R_η(gerçek) = 1.2883 > R_η(ikiz) = 1.2650` — sadakatli
> ikiz gerçekten de girişimin bir kısmını taşımıyor. **Ama** `R_η(λ)`
> tümsektir (tavan 1.2745) ve gerçek gaz **bütün ailenin üstündedir**:
> KALEM'in beklediği "fark λ_eş(E) = 0.777'yi öngörmeli" cümlesi
> **ÖLDÜ**. Doğrusu şudur: **η kanalındaki fazla, ailenin herhangi bir
> üyesinde bulunmayan bir fazladır** — bir parametre kayması değil,
> ailenin dışı.

## (iii) H-G3 MÜHÜRLENMEDİ — %31.4, eşik %70

> Parametresiz köprü ΔM'nin **%31.4'ünü** kapatıyor (defterin kendi
> π'siyle %27.2). g_E artığının %65'i girişim, %35'i model gücü; ΔM'nin
> **%52'si θ'dadır ve θ bu köprüde YOKTUR**. Ön-kayıt bu ölümü önceden
> yazmıştı (beklenti %41) ve ölçüm onu doğruladı — kurtarma yapılmadı.

## (iv) K4 TAM İSABET — asimetri gerçek ve ölçülür

> `λ_eş(R_X) = 0.9026` ↔ defterin `0.9089`'u (fark 0.0063). **X kanalı
> bir λ-gazıdır, η kanalı değildir.** 172-G3'ün "adres bir çarpan değil
> bir KANAL" hükmü, tamamen bağımsız bir ölçüyle (girişim oranı)
> doğrulandı — ve şimdi asimetrinin keskin biçimi var: X'te λ_eş VAR ve
> tutuyor, η'da λ_eş YOK.

## (v) K3: DC KAÇAĞI ORTA-NOKTA TARAĞININ REZONANSIDIR

> `ort(E) = Σ_q Re[hp_q κ(ω_q)]`, `μ̂²_E = ort(E)²/Var(e1)` özdeşliği
> 1.3e−13'te kapanıyor. İşaret ön-mühürlü türetimin dediği gibi negatif
> **değil pozitif** çıktı (Ö1 öldü) ve sebebi ölçüldü: site **orta
> noktadır**, bu da κ'ya `cos(πτ)` takar ve **τ = 0.5'te işaret döner**
> (`⟨cos Δφ⟩`: −0.85 → +1.00; `|κ|` medyanı (0.45,0.50] bandında
> 1.17e−03'e çöküyor). Kaçak birkaç büyük çizgide değil (Ö3 öldü),
> **binlerce çözülmemiş yüksek-τ çizgisindedir.**
> Gerçek−ikiz açığının **%83'ü τ > 0.70'te**, **%49'u τ > 0.80'de**;
> orada gerçek gazın çizgi genlikleri ikizin %77'si ve faz uyumu ikizin
> %78'i. **172'nin −%23.8'i budur: merdivenin çözülmemiş tepesinde
> hem daha zayıf hem daha dağınık fazlı bir tarak.**

## (vi) 162'NİN "YIKICI GİRİŞİM"İ ÇİZGİLER ARASINDA DEĞİLDİR

> Bant defteri (fit yok, kalıntı 2.6e−16): bant içi `P_b/V_b` 0.90–0.97
> (**yapıcı**), bantlar arası kovaryans **+0.0374** (**yapıcı**),
> `Var(η_çizgi)/P = 1.62`. `R_η > 1`'in mekanizması **çizgi ↔ artık**
> yıkıcı girişimidir: `Kov(η_çizgi, η_artık) = −0.0435`. 162'nin
> ölçtüğü olgu doğrudur; adı bu koşuda düzeltildi.

---

## Sıradaki adım (bu ölçümün işaret ettiği)

1. **θ'nın taşıyıcısı hâlâ bulunmadı ve ΔM'nin %52'si orada.** 174'ün
   ön-kayıtlı üçüncü-moment vekili (m3) ters yönde öldü; ön-mühürsüz
   ORAN vekili (`m3/m3_çizgi`, λ_eş = 0.527) doğru yönde ama ıska.
   Doğru nesne θ'nın kendi tanımındaki `W`-çözünürlüklü orandır:
   `⟨e1x1²e^{−iWs}⟩ / ⟨EX²e^{−iWs}⟩`. Bu, 165'in `ongor`'uyla ucuzdur
   (gaz başına ~1 dk) ve λ merdiveninde koşulmadı.
2. **`κ`'nın τ > 0.5 fazlası (×1.9 → ×4.2) ölçüldü, adı kondu,
   sayılmadı.** `κ(ω_Q)`'daki çarpımsal (`q₁q₂ = Q`) katkı 165'in
   `kule_cozumleri`/alt-merdiven sayımıyla doğrudan sayılabilir. Eğer
   fazlayı o veriyorsa, DC kaçağı **doğrudan aritmetik** bir nesnedir.
3. **Gerçek gazın "aile dışı" nicelikleri bir aile oluşturuyor:**
   `R_η`, `m3`, `skew(x1)`, `skew(ds)`, `skew(e1)`. Hepsi normalize
   oran/üçüncü moment. Yeni bir eksen (λ değil) aranmalı — 162'nin
   önerdiği bispektrum ekseni bu beş sayıyı birden içeriye alabilir mi?
4. **`R_η(λ)` tümseğinin tepesi λ ≈ 0.70-0.85'te.** 172-G4'ün iki
   makine-ölçeği (λ_c = 1.011, λ* = 0.6487) ile ilişkisi ölçülmedi;
   `Q_E(λ)` tümseğinin tepesi 1.0009-1.019 idi, `R_η = π_E` tümseğininki
   ise 0.70-0.85 — **ikisi AYNI eğrinin iki parçası** (`Q = ρ − π`), ve
   tepe yerlerinin farkı doğrudan `ρ_E(λ)`'nın eğimini veriyor. Ucuz.
5. **Ĉ kanalının λ_eş'i (0.8272) bu koşuda ilk kez ölçüldü** ve 162'nin
   "Ĉ'de girişim yok" okumasının bir aile-rastlantısı olduğunu
   gösterdi. 158'in dokuz gazında sınanmadı.
