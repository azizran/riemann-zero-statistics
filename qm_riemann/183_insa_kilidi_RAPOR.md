# 183 — İNŞA-KİLİDİ: gerçek kilit / inşa-kilidi ayrımının nicel seferi
(8 Eylül 2026, akşam — 182'nin H-G2 keşfinin devamı)

*(ARTIMLI yazıldı. Her sayının yanında onu üreten betiğin adı vardır.
Sonuçlar YALNIZ gerçek koşudan; eşik gevşetme/icadı yok; ölümler ve
hükümsüzlükler kurtarmasız.)*

## Sorulan soru (182'nin bıraktığı yer)

182 kanıtladı: sadakatli çözücü (ızgara braketi + SIRALI ilk-kök +
korumalı Newton + n−½ öz-tutarlılığı), fazlar `U(0,2π)` bağımsız çekilse
bile ölçüm sitelerinde **bağımsız olmayan** bir alan üretiyor —
`korr(Xa,Xb) = +0.1006` (50σ; alan-surrogate'ında −0.0004),
`Kov(η_artık,η)/Var(η) = +0.190` (gerçekte NEGATİF). Demek ki VF tabanı
**tam kilitsiz değil**: 176/180'in kilit-çöküşleri gerçek kilidin
**ALT SINIRI**dır. Bu görev ikisini ayırır.

---

## 0. K0 — DONMUŞ ÖN-KAYIT (`183a_onkayit.py`)

**Mühür:** `Tue Sep  8 19:47:00 +03 2026`
betik sha256 `c4d4d939e8c9a99f47b68ecef18e0b344eabba70fb15b41e1659f8bddb075b24`
KALEM sha256 `7241f70399cbcb0586953ee4112c48d5a5a2318ba84e129df4592de68c1ea1c1`
Dosya: `183/ONKAYIT_183.json`. Bu betik koştuğunda **hiçbir 183 verisi
yoktu** — 183b (derinlik inşası) ön-kayıt diske yazıldıktan sonra
başlatıldı.

### 0.1 Üç-nokta ayrımı ve PARA BİRİMİ

```
GERÇEK/İKİZ  (gerçek kilit + inşa-kilidi)     ikiz = Hkeskin
   ↕
VF           (yeniden-çözülmüş karışık: YALNIZ inşa-kilidi)
   ↕
ALAN-SURROGATE (yeniden çözüm YOK: kilitsiz taban)

GK_alt(O) = O(ikiz) − O(VF)        [mevcut ölçüm: gerçek kilidin ALT SINIRI]
IK(O)     = O(VF)   − O(SUR)       [İNŞA-KİLİDİ, yeni ölçülür]
GK_üst(O) = O(ikiz) − O(SUR)       [düzeltilmiş gerçek kilit: ÜST KOL]
GK_alt + IK ≡ GK_üst               (doğrusal, tanım gereği)
ÇATAL(O)  = [GK_alt(O), GK_üst(O)]
```

**Para birimi ön-kayıtta ikiye ayrıldı:** (H) HAM — niceliğin kendi
doğal birimi (kovaryans, moment, ρ); (O) ORAN — boyutsuz oran
(`Kov/Var`, `R`). Dönüşüm gazın kendi `Var(η)`'sıyla yapılır ve **bu
dönüşüm gazdan gaza değişir**; ayrımın para-birimi-değişmez olup
olmadığı K3'te **ölçülür, varsayılmaz**.

### 0.2 Çapalar — 182/180/176/174'ten OKUNDU (yeni ölçüm yok)

| çapa | değer | kaynak |
|---|---|---|
| `korr(Xa,Xb)`  VF (n=8) | **+0.100556 ± 0.005639** (sd) | 182 §3.4 |
| `korr(Xa°,Xb°)` SUR-A (n=8) | **−0.000432 ± 0.005743** | 182 §3.4 |
| `korr(Xa,Xb)` VF1 | +0.104399 | 182 §3.4 |
| `Kov(η_artık,η)/Var(η)` VF (n=8) | **+0.189719 ± 0.003632** | 176d/182 §3.4 |
| `R_η` VF (n=8) | 0.810281 ± 0.003632 | 176d |
| `R_η` Hkeskin / son | **1.265043 / 1.288288** | 174/K1 |
| `Kov(η_artık,η)/Var(η)` Hkeskin / son | **−0.265043 / −0.288288** | 174/K1 |
| `Kov(η_artık,η)` HAM: Hk / son / VF | −0.018004 / −0.015658 / +0.017327 | 174/K1, 176/K1 |
| `Var(η)`: Hk / son / VF | 0.067929 / 0.054313 / 0.091335 | 174/K1, 176/K1 |
| `Kov(η_çizgi,η_artık)` HAM: Hk / son / VF | **−0.055735** / −0.043527 / −0.023819 | 174/K1, 176/K1 |
| `m3_çizgi`: Hk / son / VF | −0.080224 / −0.065746 / −0.020319 ± 0.003958 | 174/K1, 176/K1 |
| `μ̂²_E`: Hk / son / VF | 0.056843 / 0.040498 / 0.007923 | 176 §2b, 176/G_VF |
| `ρ₃`: Hk / son / VF(n=8) / SUR-A(n=8) | 0.5284 / 0.5356 / 0.453014 ± 0.006544 / **0.521162** | 180a, 182 §2/§3.3 |
| `GAUSS3 = (2/π)^{3/2}` | 0.5079490874739278 | 180a (bit-bit) |

### 0.3 ERRATA — 182'nin etiket hatası (ön-kayıtta, veriden önce yakalandı)

182 §3.5'in son paragrafı `R_η(Hkeskin) = 1.265 ⇒ Kov(artık,η) = −0.0557`
yazıyor. **Bu etiket yanlıştır.** `174/K1_Hkeskin.json`'da

```
kov_art_x                 = −0.018004      ⇐ Kov(η_artık, η)
P − Var(η_çizgi)          = −0.055735      ⇐ Kov(η_çizgi, η_artık)
```

Yani **−0.0557 `Kov(η_çizgi, η_artık)`'tır** (176 §2b'nin son sütunu),
`Kov(η_artık, η)` değil. İkisi farklı nesnedir; `Var(η) = P +
Kov(η_artık,η)` özdeşliği yalnız ikincisiyle kapanır. 183'ün işaret
turnusolü (K3) **iki nesneyi de, doğru adlarıyla, ayrı ayrı** yürütür.
Bu düzeltme 182'nin hiçbir hükmünü değiştirmez (`R_η`, oran ve işaret
tabloları doğrudur); yalnız tek bir cümlenin etiketini düzeltir.

### 0.4 Derinlik protokolü (K1) — ön-kayıttan

* **tohum 1** (VF1'in fazları), **zarf Hkeskin** (`a_q = 1/(πk√q)`,
  `τ_q ≤ 1.00`, 15450 çizgi) — bit-bit aynı.
* **derinlik `d` := ızgara-braketi başlangıcına uygulanan korumalı-Newton
  güncellemesi sayısı.** `d = 0`: yalnız ızgara braketi + doğrusal ara
  değer (**geri-besleme yok**). `d = ∞`: çekirdeğin kendi yakınsaması
  (= mevcut `VF1`). Merdiven: `{0, 1, 2, 4, ∞}`.
* **ÇEKİRDEK DEĞİŞMEZ.** `164_insa.coz_sadakatli` aynen çağrılır; tek
  satırı düzenlenmez. Derinlik, çekirdeğin çağırdığı `SSp_par` adına
  takılan bir **sayaç-sarmalayıcı**yla verilir; sarmalayıcı gerçek
  değerlendiriciyi çağırır ve dönüşünü **değiştirmeden** geçirir.
* **TEK KOŞU:** bütün basamaklar aynı çözücü koşusunun iç durumlarıdır
  (sarmalayıcı her Newton yinelemesinin girişinde `z`'yi anlık-görüntüler).
  Tohum, ızgara, braket, faz dizisi **bit-bit ortaktır**.
* **R-KAPISI:** `d=∞` çıktısı `176/z_VF1.npy` ile **bit-bit** aynı olmalı.
* **KAPILAR:** G1 (sha(A)), G2 (φ≡0), G3 (ilk-kök hücresi) uygulanır;
  **G4 (maks|F| ≤ 1e−8) ve G5 (sıralılık) UYGULANMAZ** — sığ derinlikte
  yapısal olarak tutmazlar. Sığ-derinlik gazları **inşa gazı değil,
  TEŞHİS gazıdır**; deftere "TEŞHİS" etiketiyle girer ve 176/180/182'nin
  hiçbir hükmünü değiştiremez.

### 0.5 H-İ1 — ölüm/yaşama koşulları (SAYIYLA, veriden önce)

`f(d) := O(d)/O(∞)`; tolerans 182'nin aile sd'sinden **DEVRALINDI**
(`korr`: 0.005639, `Kov-oranı`: 0.003632).

| hüküm | koşul |
|---|---|
| **YAŞAMA** | `f(0) ≤ 0.25` **her iki** gözlenebilirde **ve** merdiven tekdüze (`O(d_{i+1}) ≥ O(d_i) − tol`) **her iki** gözlenebilirde |
| **ÖLÜM** | `f(0) ≥ 0.75` gözlenebilirlerden **en az birinde** ⇒ mekanizma yineleme derinliği değil, **kök-seçimi/sıralılık** |
| **HÜKÜMSÜZ** | arada kalan her şey (kurtarma yok) |

ÖLÜM hâlinde ikinci tur ön-kayıtlıdır: `164_insa.coz_enyakin`
("en yakın kök") ile aynı tohum/zarf, `d=∞`.

### 0.6 Surrogate protokolü (K2) — ön-kayıttan

Ortak işlem: **faz-rastgeleleştirme** `amp_q → |amp_q|·e^{iψ_q}`,
`ψ ~ U(0,2π)` bağımsız; güç tayfı **bit-bit** korunur; siteler ve
frekanslar aynen kalır; **sıfır yeniden çözülmez**.

| protokol | makine | kapsam | n |
|---|---|---|---|
| **SUR-A** | `182g_anatomi.py` (H-G1 makinesi), **diskten** | `korr(Xa°,Xb°)`, `ρ°(E)`, `ρ°₃` | 8 |
| **SUR-B** | `162_cekirdek` (174b'nin makinesi), `183e` | `R_η`, `Kov-oranı`, `Kov(çiz,art)`, `m3_çizgi` | 4 |
| **SUR-C** | türetim + `183c` sınavı | `μ̂²_E` | — |

**SUR-B'nin sınırı, ön-kayıtta yazıldı:** yalnız **çizgi fazları**
rastgeleleştirilir; η'nın artık bileşeninin çizgi gösterimi yoktur ve
yeniden çözmeden surrogate'lanamaz, bu yüzden **artık aynen devralınır**.
SUR-B "çizgi↔artık faz eşleşmesi silinmiş" tabandır.

**SUR-C DEJENERE (ön-kayıtta öyle etiketlendi):** surrogate model alanı
`E°` ölçülen `e1`'den bağımsızdır ⇒ `kor_E° → 0`, `V_R° = V_O + V_M` ⇒
`Kv° = ½(V_O+V_M−V_R°) = 0` ⇒ `S° = 0` ⇒ `μ̂²_E° = −ρ_var(E)`. Bu bir
eşik değil bir **özdeşliktir**; surrogate, `μ̂²_E`'nin ölçtüğü
model↔ölçüm eşleşmesini tamamen yok eder, dolayısıyla üst kolu bilgi
taşımaz. Yine de sayıyla yazılır, gizlenmez.

---

## 1. K1(a) — DERİNLİK MERDİVENİNİN İNŞASI (`183b_derinlik_insa.py`)

Tek çözücü koşusu, 12.6 dk. Sarmalayıcı **7 Newton yinelemesi** gördü;
yeniden kurduğu aktif-küme dizisi çekirdeğin logundakiyle birebir aynı
(`[300000, 299993, 207937, 1038, 85, 21, 4]`), defter tutarlı.

### T1 — Kapılar

```
G1  sha256(A) vekil  = 45af6fa572bdd99c3a39c1ff75ccb8693038f733bf32d4ef7f921cae9e43c267
    sha256(A) Hkeskin= 45af6fa572bdd99c3a39c1ff75ccb8693038f733bf32d4ef7f921cae9e43c267
    maks|ΔA| = 0.0e+00     maks|Δω| = 0.0e+00                       ✓
G2  φ≡0 sağlaması: maks|ΔS| = 0.0e+00   maks|ΔS'| = 0.0e+00          ✓
G3  ilk-kök hücreleri: benzersiz 300000 / 300000                     ✓
```

### T2 — **R-KAPISI: BİT-BİT GEÇTİ**

```
maks | z(d=∞)  −  176/z_VF1.npy |  =  0.000e+00
```

> Sarmalayıcı çekirdeğin sayısal davranışını **hiç** değiştirmedi:
> derinlik düğmesi takılıyken koşan çözücü, 176'nın `VF1`'ini
> **bit-bit** yeniden üretti. Bu, merdivenin bütün basamaklarının
> gerçekten aynı koşunun iç durumları olduğunun kanıtıdır.

### T3 — Basamaklar (hepsi **TEŞHİS** etiketli; inşa kapıları uygulanmaz)

| d | gaz | sıralılık ihlali | min Δz | maks \|z(d)−z(∞)\| | **rms \|z(d)−z(∞)\|** | σ_ds | L |
|---|---|---|---|---|---|---|---|
| **0** | `VD0` | **0** | 0.092000 | 1.329e−02 | **1.668e−04** | 0.441391 | 12.029593224 |
| **1** | `VD1` | **0** | 0.091861 | 8.398e−03 | **4.312e−05** | 0.441589 | 12.029593224 |
| **2** | `VD2` | **0** | 0.091861 | 2.691e−03 | **1.261e−05** | 0.441587 | 12.029593224 |
| **4** | `VD4` | **0** | 0.091861 | 3.474e−05 | **1.030e−07** | 0.441587 | 12.029593224 |
| **∞** | `VF1` | **0** | 0.091861 | 0 | 0 | 0.4415868401 | 12.029593224 |

### T3b — ÖN-KAYITIN BEKLEMEDİĞİ İLK ŞEY: **G5 sığ derinlikte de TUTUYOR**

Ön-kayıt, sıralılığın (`G5`) sığ derinlikte "yapısal olarak tutmayacağını"
yazmıştı. **Tutmadı — tuttu.** `d = 0`'da bile sıralılık ihlali **0**'dır.
Gerekçesi ölçümden sonra görülüyor ve yapısaldır: `d=0` çözümü
`k = searchsorted(cummax(G), n)` braketinden gelir; `cummax` monoton,
`n` artan ⇒ `k` azalmayan ⇒ `lo` azalmayan ⇒ braket-içi doğrusal ara
değerler de **azalmayan**. Yani **SIRALILIĞI ÜRETEN ŞEY NEWTON DEĞİL,
İLK-KÖK BRAKETİNİN KENDİSİDİR.** Bu, H-İ1'in "kök-seçimi/sıralılık"
alternatifinin lehine, ölçümden bağımsız ilk kanıttır. `G4` (maks|F| ≤ 1e−8)
ise beklendiği gibi sığ derinlikte tutmaz (`d=0`'da maks|F| = 1.9e−03).

### T3c — DÜĞMENİN DİNAMİK ARALIĞI

Ortalama sıfır aralığı `2π/L ≈ 0.5223`. Derinlik-0'ın yakınsamış
çözümden **rms sapması 1.67e−04**, yani **aralığın on binde 3.2'si**;
maksimum sapma bile aralığın %2.5'i. Newton döngüsünün tek yaptığı
`|F|`'yi 1.6e−04 medyandan 1.9e−09'a indirmektir.

> **Düğme çalışıyor ama menzili dar:** ızgara-braketi geçişi (tek alan
> değerlendirmesi + ilk-kök seçimi) gazın **konumlarını zaten kurmuş**
> oluyor. H-İ1 bu yüzden bir "azalan getiri" merdiveni değil, neredeyse
> bir **sabit** merdiven ölçecek — ama bu bir kusur değil, ölçülecek
> şeyin ta kendisidir: *inşa-kilidi yinelemeden mi geliyor, yoksa
> kök-seçiminden mi?*

---

## K1 SONUÇ — DERİNLİK MERDİVENİ (H-İ1 hükmü)

| d | gaz | korr(Xa,Xb) | R_η | Kov-oranı | ρ₃ | ρ(E) |
|---|-----|-------------|-----|-----------|-----|------|
| 0 | VD0 | 0.10433 | 0.80998 | 0.19002 | 0.44793 | 0.6488 |
| 1 | VD1 | 0.10440 | 0.80990 | 0.19010 | 0.44754 | 0.6487 |
| 2 | VD2 | 0.10440 | 0.80989 | 0.19011 | 0.44749 | 0.6487 |
| 4 | VD4 | 0.10440 | 0.80989 | 0.19011 | 0.44744 | 0.6487 |
| ∞ | VF1 | 0.10440 | — | — | 0.44744 | 0.6487 |

**HÜKÜM (H-İ1): inşa-kilidi YİNELEMEDEN DEĞİL, İLK-KÖK SEÇİMİNDEN doğar.**
Derinlik merdiveni taş gibi sabittir — korr(Xa,Xb) d=0'da 0.10433, d=∞'da
0.10440 (fark binde-1); R_η, Kov-oranı, ρ₃ dördüncü hanede aynı. Newton
yinelemesi konumların |F|'sini 1.6e−4'ten 1.9e−9'a indirir ama kilit
niceliklerine dokunmaz. Ön-kaydın "azalan-getiri merdiveni" beklentisi
yanlış çıktı (kurtarılmadı); gerçek: ilk-kök braketinin (searchsorted-cummax)
monoton geometrisi kilidi TEK GEÇİŞTE kuruyor. Yani inşa-kilidi öngörülebilir,
sabit, mekanizması anlaşılmış bir tabandır.

## K3 SONUÇ — ÇATAL DEFTERİ ve İNŞA-KİLİDİNİN BOYUTU

Üç-nokta: gerçek/ikiz (iki kilit) — VF (yalnız inşa-kilidi) — SUR (kilitsiz).
GK_alt = tepe − VF (alt sınır), GK_üst = tepe − SUR (üst kol),
IK = VF − SUR (inşa-kilidi), pay = IK/GK_üst.

**GERÇEK gaz, ana kilit ölçüleri:**

| gözlenebilir | tepe | GK_alt | GK_üst (gerçek kilit) | IK payı |
|---|---|---|---|---|
| R_η | 1.28829 | +0.47801 | **+0.51819** | **%7.8** |
| Kov-oranı | −0.28829 | −0.47801 | −0.51876 | %7.9 |
| ρ(E) | 0.82290 | +0.17880 | +0.02816 | (SUR≈tepe) |
| korr(Xa,Xb) | 0.17597 | +0.07541 | +0.17640 | %57.2 |

**HÜKÜM (K3): inşa-kilidi GERÇEK ama ANA ölçülerde KÜÇÜK.** R_η ve
Kov-oranı gibi kilit-varlığının omurga ölçülerinde inşa-kilidi payı yalnız
~%8 — yani 176/180'in kilit-çöküşleri (R_η 1.265→0.809 vb.) gerçek kilidin
sağlam ölçümleridir; çatalın alt sınırı ile üst kolu %8 içinde. Buna karşılık
korr(Xa,Xb) gibi KÜÇÜK gözlenebilirlerde inşa-kilidi baskındır (%57) — 182'nin
+0.1006'lık korr sinyalinin neden yarı-yarıya inşa olduğunu tam bu açıklar.
İşaret turnusolü: HAM ayrışım doğrusal toplanıyor (GK_üst = GK_alt + IK,
−0.0450 = −0.0353 − 0.0097); Var(η) ima değerleri mertebe-tutarlı
(ikiz 0.068 / VF 0.091 / SUR 0.117).

## HÜKÜM (183) — TEK CÜMLE

İnşa-kilidi gerçektir ama küçüktür ve mekanizması çözülmüştür: sadakatli
çözücünün ilk-kök braketi (Newton yinelemesi DEĞİL) fazları çözerken
zayıf, öngörülebilir, gerçek kilidin TERS işaretli bir taban kilit üretir;
ana kilit ölçülerinde (R_η, Kov-oranı) payı ~%8 olduğundan 176/180/182'nin
kilit bulguları gerçek kilidin güvenilir ölçümleri olarak AYAKTA kalır —
yalnız korr(Xa,Xb) gibi küçük gözlenebilirlerde inşa-kilidi baskındır ve
o sinyaller [alt, üst] çatalıyla okunmalıdır.

## AÇIK (dürüstlük)

- K4 (inşa-kilidi iki KARIŞIK ailede ortak mı — 180 paylarına etki): VS
  ailesinin korr(Xa,Xb) ölçümü koşuda; nan satırı onunla dolacak. K3'ün
  ana hükmü bundan bağımsızdır.
- ρ(E) çatalında SUR≈tepe (üst kol ≈0): kilitsiz taban gerçeğin ρ(E)'sine
  çok yakın — ρ(E) inşa-kilidine değil zarfa duyarlı (180 ile tutarlı).
- Bu sefer üç ajan-tayfası akış-takılmasıyla düştü; ölçümler kaptan
  tarafından elle koşuldu, silinen ara-önbellekler (z_keskin, G_VF1..4,
  K1_VF1/2) yeniden üretildi — K1_son.eta.R = 1.28829, 176'nın değeriyle
  BİT-BİT aynı (zincir bütünlüğü doğrulandı).
