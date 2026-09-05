# 178 — BANT-AĞIRLIKLI `g_X` KESTİRİCİSİ: HAKEM RAPORU

**Sefer:** 178 · **Hipotez:** `H-F1b^ba` (θ'nın kilit-duyarlılığı, bant-ağırlıklı
`g_X` para biriminde) · **Bakış sayacı:** k = 2 · **Tohumlar:** VF1..VF4 (n = 4,
tavan) · **Yeni gaz koşusu:** SIFIR · **Git:** dokunulmadı.

> # HÜKÜM: **HÜKÜMSÜZ — ve KALICI HÜKÜMSÜZ**
>
> Donmuş kuralın mekanik çıktısı `OLDU|ODENEMEZ`'dir ve iki bağımsız hesap bunda
> 3 × 10⁻¹² göreli fark içinde birleşir. **Ama hakem bu çıktıyı hükme
> çeviremez:** çürütücünün ölümcül kusurları ayaktadır ve hakemin kendi
> hesabı (Bölüm 6) bunu kesinleştirir — para birimi değişikliği, ölçtüğü her
> gazı gaza-özgü bir miktarda kaydırmıştır ve bu kaydırmanın tohumlar arası
> farkı, **sınanan etkinin 10.9–12.3 katıdır**. `ÖLDÜ` sonucu hipotezin değil,
> para biriminin ölçümüdür.
>
> 177'nin eski-para hükmü (**H-F1b KALICI HÜKÜMSÜZ**) ayakta kalır ve geri
> alınmaz. Ş2.6 + n = 4 tavanı gereği bu dört tohumda **ÜÇÜNCÜ kestirici
> DENENMEZ**; soru ancak yeni veriyle (yeni tohum / yeni gaz, ayrı sefer)
> açılır.

---

## 0. KÜNYE VE DENETİM ZİNCİRİ

Bütün sha256'lar hakem tarafından diskte yeniden hesaplandı
(`shasum -a 256`, 2026-09-05).

| dosya | sha256 | mtime (yerel) |
|---|---|---|
| `178_configs/178a_onkayit.py` | `f81dac6d8c3acebfc61d49d698ab2bed92ead407c479531489b40db7c6b7a1f6` | 17:55:38 |
| `scratchpad/178/ONKAYIT.json` | `0d084308a454428ce7478eba8090d7e5d15cdf691f06ccf85cb6c1ca86a7542f` | **18:03:00** |
| `178_configs/178b_hesap_A.py` | `63d6012e9ff3d359cd209626c1c374cdb078cffa1c5d4b6bffe7219decf47357` | 22:35:43 |
| `178_configs/178b_hesap_B.py` | `ee5684ac73aa6de9196608b149e9f5a8a395b23f196958c3e001065dc549fd42` | 22:41:01 |
| `scratchpad/178/HESAP_A.json` | `431eddf548724a61302fe979d4d2e72d1367f7daee913c8d0890a3ee6c924a28` | 22:44:55 |
| `scratchpad/178/HESAP_B.json` | `a4d0ffd124b35404edbdfce11c41a4abeca1caa1db773aa0cc2c7091efd110bd` | 22:46:26 |

**Ön-kayıt zaman damgası:** `2026-09-05 18:03:00 +0300`
(`ONKAYIT.json → "zaman"`, dosya mtime'ı ile birebir).
Ön-kayıt içindeki betik sha'sı (`f81dac6d…b7a1f6`) diskteki
`178a_onkayit.py` ile **eşleşir**. Ölçüm (`HESAP_A/B`) mühürden **4 saat 42
dakika sonradır** — sıra doğrudur.

Ön-kayıt kendinden öncekileri de gömüyor:
`onkayit_176 = 6173943f0f4ef0488dced150f9e5a767e13c6bb24af2837527d2bd3149690244`,
`onkayit_177 = 03a9b9addedf738b4f4b03f0dedd0d54848fd11764cb4100e9490630cd12aebe`
(`ONKAYIT.json → "onceki_onkayitlar"`).

**Hakemin kendi betikleri** (yalnız okur, hiçbir şeye yazmaz):
`scratchpad/178/hakem_capraz.py`, `scratchpad/178/hakem_mekanizma.py`.

---

## 1. KESTİRİCİ BEYANI — ne değişti, ne değişmedi

`ONKAYIT.json → "kestirici"`, `"beyan"`:

```
theta_bw := M / (g_E · g_X,bw²)
g_X,bw   := Σ_b W_b · g_X,b
g_X,b    := ⟨X_b·x1⟩ / ⟨X_b·X_b⟩          (165.alanlar'ın g_X biçimi, banda kısıtlı)
X_b      := sentez(s, w[lo<τ≤hi], y[lo<τ≤hi]) − ⟨·⟩ ,   x1 := Y.Xtil0
M (=KALİB_u2, lo=0.60) ve g_E  DEĞİŞMEDİ.
```

* **Bant kümesi (sabit):** `lo ∈ {0.52, 0.56, 0.60, 0.64, 0.68}`, genişlik 0.04
  ⇒ `τ ∈ (0.52, 0.72]`. Üyelik `lo < τ_q ≤ hi`.
* **Ağırlıklar:** `W_b ∝ 1/s_b²`, `Σ W_b = 1`, doğrusal uzayda, floor/kırpma yok.
  `W = [0.026030662647056793, 0.10131948125676836, 0.5462804304510797,
  0.2475890197938992, 0.07878040585119582]`
  (`ONKAYIT.json → kestirici.agirlik.W`; hakem `1/s_b²`'den yeniden türetti,
  **en kötü fark 0.000e+00**, `hakem_capraz.py` Bölüm 2).
* **Ağırlıkların kaynağı:** YALNIZ ikiz `Hkeskin`. Vekil dosyaları
  `178a_onkayit.py:256 _bekci()` körlük bekçisiyle kilitli — hakem kodu okudu,
  bekçi gerçek ve etkindir (vekil adı geçen her yol `SystemExit` atar).
  **Bu yönde sızıntı yoktur.**
* **Beyan:** *"BU BİR VARYANS-KÜÇÜLTME KESTİRİCİ DEĞİŞİKLİĞİDİR; EŞİK GEVŞETME
  DEĞİLDİR."* — Bölüm 5/K1'de bu beyanın **kendi sayılarıyla yanlışlandığı**
  gösterilecektir.

---

## 2. ÇAPALAR VE EŞİKLER (yeni para birimi)

`ONKAYIT.json → "capalar"`, `"esikler"`. Hakem her birini bağımsız yeniden
türetti (`hakem_capraz.py` Bölüm 2); farklar aşağıda.

| gaz | `g_X` global (donmuş) | `g_X,bw` | `θ_eski` | `θ_bw` |
|---|---|---|---|---|
| `Hkeskin` (ikiz) | 0.7042132063265797 | 0.9822421389493694 | 0.8933338132966645 | **0.4591821314567635** |
| `son` (gerçek) | 0.7047334771883751 | 0.9863866175580883 | 0.9219381611740145 | **0.4706057764995118** |
| `HA4` (erfc) | 0.7167735777866932 | 0.9942174293830159 | 0.7235729225877114 | **0.37608285664653684** |

| nicelik | ön-kayıt | hakemin yeniden türetmesi | fark |
|---|---|---|---|
| `D = Δlogθ_bw(son←Hk)` | +0.024573819057101542 | +0.024573819057101504 | 3.8e−17 |
| `E = Δlogθ_bw(erfc←Hk)` | −0.19963744938033756 | −0.19963744938033748 | 8.3e−17 |
| `f_θ = D/E` | −0.1230922311088284 | −0.12309223110882826 | — |
| `kesim_θ = f·D` | −0.0030248462146032747 | — | — |
| `kilit_θ = (1−f)·D` | +0.02759866527170482 | +0.02759866527170477 | 4.9e−17 |
| **`eşik_a` = θ_bw(Hk)·e^{−D}** | **0.4480357878679697** | 0.44803578786796977 | 5.6e−17 |
| **`eşik_b` = θ_bw(Hk)** | **0.4591821314567635** | özdeş (`==` doğru) | 0 |

`eşik_a`'nın ikinci yolu (`θ_bw(Hk)²/θ_bw(son)`) hakemde birinci yolla **birebir
aynı** çıktı (fark 0.000e+00); ön-kayıttaki 5.6e−17'lik ayrım yalnız numpy/math
toplama sırası farkıdır, maddi değildir.

**Formül sağlaması (eski para):** `0.8933338133²/0.9219381612 = 0.8656169530535`
= 176-F3'ün eşiği (`ONKAYIT.json → esikler.saglama_eski_para`, özdeş). Eşik
**formülleri** 176-F3 ile harf harf aynıdır — bu doğrudur ve teyit edilmiştir.

**Ön koşullar:** P1 (`D > 0`) ✓, P2 (`kilit_θ > 0`) ✓.

**Özdeşlik denetimi:** `|θ_bw − M/(g_E·g_X,bw²)|` üç çapada da hakemde
**0.000e+00** (`hakem_capraz.py` Bölüm 2). *(Çürütücü Kusur 13 haklıdır: bu
denetim totolojiktir — `θ_bw` zaten bu ifadeyle hesaplanıyor. Geçmesi hiçbir
hata sınıfını dışlamaz.)*

---

## 3. GÖREV (1) — A ↔ B ÇAPRAZ DOĞRULAMA: **GEÇTİ**

İki hesap tamamen ayrı makinelerdir: **A** numpy yolu
(`178b_hesap_A.py`), **B** saf Python (`178b_hesap_B.py`; betik başında
`if "numpy" in sys.modules: SystemExit` denetimi ile numpy hiç yüklenmez,
`.npy/.npz`'ler `zipfile+struct+array` ile elle çözülür, toplamlar
`math.fsum`).

`hakem_capraz.py` Bölüm 1 — **286 alan** karşılaştırıldı (dört vekilin
`θ_bw`, `g_X,bw`, `M`, `g_E`, `θ_eski`, `g_X` global; 20 bandın her birinde
`n`, `lo`, `g_X,b`, `pay`, `payda`, `s_b` ve **8 jackknife değeri**; `ȳ`, `s`,
`SAÇ_4`, `SAÇ^kat`, marjlar, `Δ_i`, `Δ̄`, `SAÇ_4(Δ)`, `ω_θ`, tohum-başına `ω`,
`R_bant_min`):

```
karşılaştırılan alan sayısı : 286
1e-9'u aşan alan sayısı     : 0
en kötü göreli fark         : 2.963e-12   (s_b[VF4, bant 3])
```

| nicelik | A (`HESAP_A.json`) | B (`HESAP_B.json`) | göreli fark |
|---|---|---|---|
| `θ_bw(VF1)` | 0.5991599552430427 | 0.5991599552430432 | 9.3e−16 |
| `θ_bw(VF2)` | 0.6019514482870469 | 0.6019514482870489 | 3.3e−15 |
| `θ_bw(VF3)` | 0.6090846668719566 | 0.6090846668719580 | 2.2e−15 |
| `θ_bw(VF4)` | 0.6018892959067726 | 0.6018892959067732 | 9.2e−16 |
| `ȳ` | 0.6030213415772048 | 0.6030213415772059 | 1.8e−15 |
| `SAÇ_4 = s(ddof=1)` | 0.004246583611324138 | 0.004246583611324387 | 5.9e−14 |
| `Δ̄` | −0.27248712049360957 | −0.27248712049361146 | 6.9e−15 |
| `ω_θ` | −0.10128429270972486 | −0.10128429270972417 | 6.9e−15 |
| H1 / H2 / bileşim | `OLDU` / `ODENEMEZ` / `OLDU\|ODENEMEZ` | aynı | — |

En büyük fark (`s_b`'lerde ~3e−12) jackknife farklarının kayan-nokta iptalinden
gelir ve hiçbir hükme girmez (`s_b` yalnız tanıdır; ağırlıklar ön-kayıttan
**sayı olarak** okunur).

⇒ **Kestirici muğlaklığı YOKTUR.** İcra düzeyinde ön-kayıt metni iki bağımsız
uygulayıcıyı aynı sayıya götürmüştür. *(Çürütücü Kusur 12'nin muğlaklık
listesi meşruiyet/yeniden-üretilebilirlik düzeyinde geçerlidir — ama bu iki
uygulayıcı arasında fiilen ayrışma üretmemiştir; ikisi de aynı ön-kayıt
metnini aynı biçimde okumuştur.)*

---

## 4. DONMUŞ KURALIN MEKANİK ÇIKTISI (hakem yeniden koştu)

`hakem_capraz.py` Bölüm 3, hem A'nın hem B'nin `θ_bw`'leriyle, ön-kayıt
kuralını sıfırdan uygulayarak:

```
ȳ        = 0.6030213415772048          (A) / …2059 (B)
s(ddof=1)= 0.004246583611324138   =  SAÇ_4   (n=4 ⇒ SAÇ_4 = s)
SAÇ^kat  = 0.009924711628913951        (max−min; BAĞLAYICI DEĞİL)
göreli saçılma s/ȳ = 0.704218 %
```

**H1 — dal (a):** `ȳ ≤ eşik_a` mı? `0.6030213 > 0.4480358` ⇒ `taraf_a = False`.
Marj `|ȳ−eşik_a| = 0.15498555370923506 ≥ SAÇ_4` ⇒ kesinlik var **ama yanlış
tarafta**. **H-F1b^ba YAŞAMADI.**

**H1 — dal (b):** `ȳ ≥ eşik_b` mı? `0.6030213 ≥ 0.4591821` ⇒ `taraf_b = True`.
Marj `|ȳ−eşik_b| = 0.14383921012044126 ≥ SAÇ_4` ⇒ **KESİN**. Marj SAÇ_4'ün
**33.87 katı**; bağlayıcı olmayan katı okumayla (`SAÇ^kat`) bile 14.5 katı.
⇒ **H1 = `OLDU`.**

**H2:** `Δ_i = logθ_bw(Hk) − logθ_bw(VF_i)` =
`−0.26608166768558705 / −0.27072985945499095 / −0.28251035219357560 /
−0.27062660264028470`; `Δ̄ = −0.27248712049360957`;
`SAÇ_4(Δ) = 0.007024824847814011`; `SAÇ^kat(Δ) = 0.016428684507988567`.
`ω_θ = kilit_θ/Δ̄ = −0.10128429270972486` — **(0,1] bandının dışında ve
NEGATİF**. Üst kenar (`Δ̄ ≥ kilit_θ`) sağlanmadı, marj 0.30008578576531440 ≥
SAÇ ⇒ **kesin ihlâl**; alt kenar (`Δ̄ > 0`) sağlanmadı, marj 0.27248712049360957
≥ SAÇ ⇒ **kesin ihlâl**. ⇒ **H2 = `ODENEMEZ`.**
Tohum başına `ω`: −0.10372254 / −0.10194172 / −0.09769081 / −0.10198061
(dördü de negatif, dağılım dar).

**Kapı ölümü yok:** 20 bandın hepsinde `payda > 0`, `θ_bw` sonlu; hiçbir vekil
düşürülmedi, dördü de ortalamaya girdi (`HESAP_A.json → kapilar.KAPI_OLUMU`,
`HESAP_B.json → kapi_olumu`).

**K-KİMLİK kapısı:** dört vekilde de `|Δg_E| = |Δg_X| = |Δθ_eski| = 0.000e+00`
(tol 1e−12) — ölçüm makinesi 176c/167 zinciriyle **bit düzeyinde aynı
makinedir**. B ayrıca üç çapayı saf-Python'la yeniden ölçtü: en kötü
`|Δg_X,b| = 6.33e−15`, bant çizgi sayıları (53/75/116/176/264) birebir. **Bu
kapı gerçektir ve gerçekten geçmiştir — raporun olumlu bulgusudur.**

**Ön-kayıtta ÖNCEDEN yazılı bileşim cümlesi** (`bilesim_tablosu`, anahtar
`OLDU|ODENEMEZ`):
> *"H-F1b^ba ÖLDÜ; F9-θ satırı ✗ fatura ödenemez — KESİN."* · mühür adayı:
> **HAYIR**.

**Bu, kuralın mekanik çıktısıdır. Hükmün kendisi değildir** — bkz. Bölüm 5–7.

---

## 5. GÖREV (2) — ÇÜRÜTÜCÜNÜN ÖLÜMCÜL KUSURLARI: HAKEM DENETİMİ

Çürütücü dört kusuru ölümcül ilan etti. Hakem **dördünü de ön-kaydın kendi
sayılarıyla yeniden hesapladı** (`hakem_capraz.py` Bölüm 4). **Dördü de
maddi olarak doğrudur.**

### K1 — "Çıta fiilen indi; beyanı kendi sayıları yanlışlıyor" ✅ DOĞRU

| | eski para (176/177) | yeni para (178) | değişim |
|---|---|---|---|
| `D` (dal-a için gereken düşüş) | +0.031517828888066625 | +0.024573819057101504 | **−%22.032** |
| `kilit_θ` (H2 çıtası) | +0.036231144031667745 | +0.027598665271704770 | **−%23.826** |

* İnen miktar `0.006944009830965121`; 177'nin kıl payı kaçırdığı açık
  (`D_eski − Δ̄_eski = 0.031517829 − 0.029611899 = 0.001905930`) ile oranı
  **3.643×**.
* `Δ̄_eski` hakem tarafından iki yoldan doğrulandı: `kilit_eski/ω_eski =
  0.029611898674445112` ve dört vekilin eski `θ`'sından doğrudan
  `0.029611898674445110` (fark 3.5e−18).
* **Karşı-olgu:** 177'nin `Δ̄`'sı değişmeden kalsa, yeni `kilit_θ` ile
  `ω = 0.027598665/0.029611899 = 0.9320126877079399 ≤ 1` — yani **177'de
  ihlâlde olan H2 üst kenarı, vekile hiç bakılmadan, yalnız kestirici
  seçilerek kendiliğinden geçerdi.**

⇒ Ön-kaydın *"EŞİK GEVŞETME DEĞİLDİR"* beyanı, ön-kaydın **kendi sayılarıyla
yanlışlanmıştır**. Formülün harf harf aynı olması çıtanın taşınmadığı anlamına
gelmez.

**Hakemin dürüstlük notu:** bu kusur, gözlenen dala göre *muhafazakârdır* —
indirilen çıta `YAŞADI` ve `ÖDENEBİLİR` yönünde çalışır, oysa sonuç `ÖLDÜ` +
`ÖDENEMEZ` çıkmıştır. K1 tek başına bir `ÖLDÜ`'yü çürütmez; **ön-kaydın
beyanını** çürütür. Ölümcüllüğü K2 ve K4 taşır.

### K2 — "İnişin kaynağı tam da dokunulan bileşen" ✅ DOĞRU (yön-bağımsız)

`2Δlog g_X`, iki çapa çiftinde iki farklı okumayla:

| çapa çifti | global okuma | bant-ağırlıklı okuma | oran |
|---|---|---|---|
| `son ← Hkeskin` | +0.0014770491787372366 | +0.0084210590097021360 | **5.701×** |
| `HA4 ← Hkeskin` | +0.0353576788613120800 | +0.0242361385386021100 | **0.686×** |

Özdeşlik tam tutuyor: `D_eski − (0.0084211 − 0.0014770) = 0.024573819057101723`
= `D_yeni` (fark 2.2e−16). Yani `D`'deki **%22'lik inişin tamamı `g_X`'in
okunuş biçiminden** gelir — dokunulan tek bileşenden.

Ve etki bir yeniden ölçekleme değildir: çapaların ölçek oranı
`c = g_X,bw/g_X,global` = 1.3948078935825203 / 1.3996590902611930 /
1.3870732127892250 — **yayılımı yalnızca 1.0091**, yani neredeyse sabit. Saf
yeniden ölçekleme olsaydı hüküm birebir 177'nin hükmü olurdu. Buna karşılık
`D` %22 kaydı: **hükmü belirleyen şey, `c`'nin o %0.9'luk çapa-içi
oynamasıdır.**

### K3 — "G-GÜÇ kapısı totolojik" ✅ DOĞRU

`178a_onkayit.py` (G-GÜÇ bölümü) tanıkları elde üretiyor:

```python
ta = [esik_a * (0.90 + 0.001 * k) for k in range(4)]   # dal (a) tanığı
tb = [esik_b * (1.10 + 0.001 * k) for k in range(4)]   # dal (b) tanığı
```

Hakem `ONKAYIT.json → kapilar.G_GUC.tanik` içindeki dörtlüleri bu formülle
yeniden üretti: **en kötü fark 0.000e+00** (her iki dalda). Herhangi bir eşik
için KESİN küme bu yolla "boş değil" gösterilebilir ⇒ kapının bloklama gücü
**sıfırdır**; Ş2.7'nin amacı (kurgulanmış hükümsüzlüğü yakalamak) hiç
karşılanmamıştır.

**Hakemin dürüstlük notu:** çürütücünün K3 içindeki *sayısal öngörüsü*
tutmadı — vekillerin `g_X,bw`'sinin 1.107–1.154 şeridine düşeceği öngörüsü
yanlış çıktı (ölçülen 0.9784–0.9820). Ama çürütücünün **adını koyduğu ikinci
senaryo** ("vekil bant-eğimleri çapalarınki gibi ≈0.98 kalırsa `ȳ ≈ 0.60` ile
kesin ÖLDÜ çıkar") **birebir gerçekleşti**. Kapının totolojik olduğu bulgusu
öngörüden bağımsızdır ve ayaktadır.

### K4 — "ÖLÇÜT-T silahsızlandırılmış; tek sayısal tanı kapının reddetmesi
gereken durumu gösteriyor" ✅ DOĞRU (yön-bağımsız)

`ONKAYIT.json → kapilar.OLCUT_T.tanim`:
> *"ikizin beş g_X,b'si sonlu + payda>0 + işaret = işaret(g_X global);
> **sayısal yayılım eşiği KASITLI OLARAK YOK** (Ş3.5)"*

Aynı kaydın raporladığı tanı — hakem yeniden hesapladı:
`χ² = 57.48384311692952` (ön-kayıt 57.48384311692951), `sd = 4`,
**`χ²/sd = 14.371`**, yayılım `max/min = 1.0342`.

Ters-varyans birleşimi, bantların **aynı** niceliği ölçtüğü varsayımına
dayanır. `χ² = 57.5` / sd = 4 bu varsayımı ezici biçimde reddeder. İkilem her
iki uçta da aleyhtedir: ya `s_b` geçerli bir hatadır (o hâlde bantlar ~14σ
düzeyinde tutarsızdır, ağırlıklama meşru değildir), ya değildir (o hâlde
`1/s_b²` "ters varyans" değildir, varyans-küçültme gerekçesi çöker).

**Bu kusur yön-bağımsızdır:** `g_X,bw`'yi tanımlayan ağırlıkları vurur,
dolayısıyla hem çapaları (⇒ `eşik_a`, `eşik_b`) hem vekilleri (⇒ `ȳ`) aynı
anda geçersizleştirir — hangi dalın ateşlendiğinden bağımsız olarak.

### Ölümcül olmayan ama teyit edilen kusurlar

* **Kusur 8 (jackknife bölüntüsü 178'in icadı, 163'e yanlış atfediliyor)** ✅
  DOĞRU. `163_cekirdek.bant_adaylari` `gap ≥ 2.5·dres` filtresinden **geçen
  aday** çizgilere, bant başına ≤220 (fazlası `rng(21)` ile rastgele altküme),
  `grup = kul % njack` atar ve çizgi kaydında dondurur. `178a_onkayit.py:303-307`
  ise `Mo.msk & (tau > lo) & (tau <= hi + 1e-12)` **merdiven** çizgilerini alır,
  `q`'ya göre sıralar ve kendi `np.arange(len(idx)) % 8`'ini uygular — farklı
  popülasyon, `grup` alanı yok. Bant sayıları da ele veriyor: `lo=0.68`'de
  **n = 264 > 220**, yani 163'ün aday listesi olamaz. Ş3.5 "gerekçesi yalnız
  178'in içinde doğan hiçbir sabit meşru değildir" der; `W` vektörünün tamamı
  böyle bir sabite dayanır.
* **Kusur 9 (`s_b` örnekleme hatası gibi davranmıyor)** ✅ DOĞRU.
  `n = [53, 75, 116, 176, 264]` ↔ `s_b = [5.046e−3, 2.558e−3, 1.101e−3,
  1.636e−3, 2.900e−3]`: en kalabalık bant ikinci en büyük hataya sahip,
  monotonluk yok.
* **Kusur 7 (destek daralması)** ✅ DOĞRU (hakem doğruladı): bantlardaki çizgi
  sayısı `Σ n_b = 684` / toplam merdiven `8981` = **%7.616**;
  `Σ⟨X_b,X_b⟩(Hkeskin) = 2433.197524449418`; ağırlığın **%87.265'i** `lo ≥ 0.60`
  bantlarındadır. *(Çürütücünün "alan gücünün %8.9'u" ve `τ_c_faz ≈ 0.577`
  sayıları hakem tarafından yeniden hesaplanmadı; çürütücünün ölçümü olarak
  kayda geçer.)*
* **Kusur 6 (indirgeme belgesi indirgemiyor)** ✅ DOĞRU: `Hkeskin`'de global
  `g_X = 0.7042132`, güç-ağırlıklı birleşim `0.9790545`, bant-ağırlıklı
  `0.9822421` — **`bw/global = 1.3948`**. Sayısal indirgeme yoktur; belge yalnız
  karşı-olgusal bir cümledir.
* **Kusur 5 (pay/payda aynı banda bağlı)**: ağırlığın **%54.628'i** `lo = 0.60`
  bandındadır ve `M = KALİB_u2` da o banttan okunur. Mekanizma olarak
  makuldür; hakem ayrıştırıcı bir ölçüm yapmadı.
* **Kusur 10/11 (mühürün kendine göndermesi; sabit-dizge denetim alanları)**
  ✅ DOĞRU: `178a_onkayit.py:627-629` `vekil_dosyasi_acilmadi=True`,
  `git="GİT'E DOKUNULMADI"`, `bosa_giden_kosu="YOK…"` hard-coded'dur; sayaç
  yoktur. JSON'un bağımsız özeti / dış zaman damgası yoktur; koruma yalnız
  `OUT.exists()`'tir (satır 336).
* **Kusur 13 (özdeşlik kapısı totolojik)** ✅ DOĞRU (Bölüm 2), ve çürütücünün
  K-KİMLİK'i gerçek sayması da doğrudur.

**Çürütücünün 2. ve 3. sorulara verdiği olumlu yanıtlar da teyit edilir:**
para-birimi tutarlılığı biçimsel olarak tamdır; ağırlıklar gerçekten vekile
bakılmadan sabitlenmiştir (`_bekci` gerçek ve etkin); eşik türetim formülleri
176-F3 ile birebirdir. **Sızıntı ağırlıklardan değil, PENCERE SEÇİMİNDEN
gelmektedir.**

---

## 6. HAKEMİN KESİCİ HESABI — para biriminin tohum-başına kaydırması

Çürütücünün K2/Kusur 6'sını hakem tam muhasebeye çevirdi
(`hakem_mekanizma.py`). Özdeşlik, her gaz için tam:

```
θ_bw / θ_eski  =  (g_X,global / g_X,bw)²          ⇒  kaydırma := 2·ln(g_X,global / g_X,bw)
```

| gaz | `g_X` global | `g_X,bw` | çarpan `k_g` | **kaydırma** = 2·ln(gg/gw) |
|---|---|---|---|---|
| `Hkeskin` | 0.7042132063 | 0.9822421389 | 0.5140095725 | **−0.665513390** |
| `son` | 0.7047334772 | 0.9863866176 | 0.5104526489 | −0.672457400 |
| `HA4` | 0.7167735778 | 0.9942174294 | 0.5197580574 | −0.654391850 |
| `VF1` | 0.8134744638 | 0.9783623314 | 0.6913347736 | **−0.369131096** |
| `VF2` | 0.8333309067 | 0.9797919595 | 0.7233812310 | **−0.323818905** |
| `VF3` | 0.8154829076 | 0.9800461192 | 0.6923675401 | **−0.367638337** |
| `VF4` | 0.8068089541 | 0.9820296264 | 0.6749820737 | **−0.393069146** |

(`k_g` ölçülen ile `(gg/gw)²` arasındaki fark yedi gazın hepsinde ≤ 2.2e−16.)

**İşaret dönüşünün tam muhasebesi** — `Δ_yeni(i) − Δ_eski(i)` tam olarak
`kaydırma(Hkeskin) − kaydırma(VF_i)`'dir, kalıntı ≤ 2.8e−16:

| tohum | `Δ_eski` | `Δ_yeni` | fark | öngörü (kaydırma farkı) | kalıntı |
|---|---|---|---|---|---|
| VF1 | **+0.030300627** | −0.266081668 | −0.296382295 | −0.296382295 | 2.8e−16 |
| VF2 | **+0.070964625** | −0.270729859 | −0.341694485 | −0.341694485 | 0.0 |
| VF3 | **+0.015364701** | −0.282510352 | −0.297875053 | −0.297875053 | 2.8e−16 |
| VF4 | **+0.001817642** | −0.270626603 | −0.272444244 | −0.272444244 | 2.8e−16 |

**Ortalama kaydırma farkı (ikiz − vekil) = −0.302099019.** Sınanan etkinin
ölçeğiyle karşılaştırın:

| | değer | oran |
|---|---|---|
| kaydırma farkı | 0.302099019 | — |
| `D` (dal-a'nın ölçtüğü etki) | 0.024573819 | **12.294 ×** |
| `kilit_θ` (H2'nin ölçtüğü etki) | 0.027598665 | **10.946 ×** |
| `SAÇ_4(Δ)` | 0.007024825 | **43.004 ×** |

Ve `D`'nin kendisi de aynı mekanizmanın çapa-içi artığıdır:
`kaydırma(son) − kaydırma(Hkeskin) = −0.006944010` ⇒
`D_yeni = D_eski + (−0.006944010)`, yani `D_eski`'nin **%22.032'si**.

### Bunun anlamı

Ölçüm makinesi bit düzeyinde aynıdır (K-KİMLİK: `|Δg_E| = |Δg_X| =
|Δθ_eski| = 0.000e+00`), gaz koşuları aynıdır, `M` ve `g_E` bit düzeyinde
değişmemiştir. Değişen tek şey `g_X`'in okunuşudur — ve o okunuş, **her gaza
kendi özel miktarını uygulamıştır**: ikize −0.6655, vekillere −0.32…−0.39.
Aradaki −0.3021, dört vekili ikizin **altından üstüne** taşımıştır.

Bir nicelik ki, çapasına göre sıralaması yalnızca kestirici değiştirilerek
tersine dönüyor ve bu dönüşün büyüklüğü sınanan etkinin **12 katı** —
o nicelik `θ`'nın yeniden ifade edilmiş hâli **değildir**; başka bir
büyüklüktür. `θ_bw`'de ölçülen `ÖLDÜ`, `θ`'nın kilit-duyarlılığı hakkında
değil, **bant penceresinin ikizle vekilleri farklı taşıması** hakkındadır.

Bu, çürütücünün Kusur 2 + Kusur 6'sının sayısal kesinleşmesidir ve
**yön-bağımsızdır**: hangi dal ateşlenirse ateşlensin, hüküm `θ` defterine
taşınamaz.

---

## 7. HÜKÜM

**Görev (1):** A ↔ B çapraz doğrulama **GEÇTİ** — en kötü göreli fark
`2.963e-12 < 1e-9`, 286 alanın hiçbiri toleransı aşmadı, iki hesap aynı H1,
aynı H2, aynı bileşim cümlesini veriyor. **"HÜKÜMSÜZ — kestirici muğlaklığı"
şıkkı kapalıdır.**

**Görev (2):** Çürütücünün dört ölümcül kusuru da **ayaktadır** ve hakem
tarafından ön-kaydın kendi sayılarıyla teyit edilmiştir (K1–K4). Bunlardan
**K2 ve K4 yön-bağımsızdır** — para biriminin kendisini geçersizleştirirler.
Ön-kayıt kuralı gereği **kurtarma yoktur**.

**Görev (3) — H-F1b HÜKMÜ:**

> ## **HÜKÜMSÜZ** — ve n = 4 tavanı + Ş2.6 gereği **KALICI HÜKÜMSÜZ**

**Gerekçe (üç ayak):**

1. **Ölçülen nicelik `θ` değildir.** (Bölüm 6, hakemin kendi hesabı.) Para
   birimi değişikliği her gaza gaza-özgü bir kaydırma uygulamıştır; ikiz ile
   vekiller arasındaki fark −0.3021, sınanan `D`'nin **12.29 katı**,
   `kilit_θ`'nın **10.95 katıdır**. Dört tohum eski parada ikizin **altında**
   (`Δ_i > 0`), yeni parada **üstündedir** (`Δ_i < 0`) — sıralama tersine
   dönmüştür. `θ_bw`, `θ`'nın daha sessiz kestirimi değil, **başka bir
   büyüklüktür** (`g_X,bw/g_X,global = 1.395`, %39 fark; sayısal indirgeme
   yoktur). Dolayısıyla `θ_bw`'deki `ÖLDÜ`, ΔM defterinin `θ` satırına
   taşınamaz.
2. **Ağırlıkların kendisi meşru değildir.** `χ²/sd = 14.371` (sd = 4),
   ters-varyans birleşiminin ön şartını reddeder; `s_b` bir örnekleme hatası
   gibi davranmaz (`n` ile monoton değil); jackknife bölüntüsü 178'in
   icadıdır ve 163'e yanlış atfedilmiştir (Ş3.5 ihlâli). Kapı bu durumu
   yakalamak için vardı ve **kasıtlı olarak silahsızlandırılmıştır**. Bu
   kusur `eşik_a`, `eşik_b` ve `ȳ`'yi aynı anda vurur.
3. **Ön-kaydın kendi beyanı yanlıştır ve güç kapısı boştur.** Çıta
   `D`'de −%22.03, `kilit_θ`'da −%23.83 inmiştir (inen miktar 177'nin açığının
   3.64 katı); G-GÜÇ tanıkları `eşik×(0.90+0.001k)` biçiminde elde
   uydurulmuştur (fark 0.000e+00). Bu ikisi tek başına bir `ÖLDÜ`'yü
   çürütmez — çünkü indirilen çıta `YAŞADI` yönünde çalışır — ama sefer
   kayıtlarının ön-kayıt disiplini iddiasını taşıyamaz hâle getirir.

**Ne değişmedi:** 177'nin eski-para hükmü — **H-F1b KALICI HÜKÜMSÜZ** —
ayaktadır ve 178 onu geri almaz. K-KİMLİK kapısının üç çapada ve dört vekilde
tam sıfır kalıntıyla geçmesi, ölçüm zincirinin bütünlüğünü kanıtlar ve
kayıtta olumlu olarak kalır.

**Ne kapandı:** Ş2.6 + n = 4 tavanı. Bakış sayacı `k = 2` tüketilmiştir. Bu
dört tohumda **ÜÇÜNCÜ kestirici DENENMEZ** — hükümsüzlüğe hangi yoldan
varıldığından bağımsız olarak. Soru ancak **yeni veriyle** (yeni tohum / yeni
gaz, ayrı sefer, ayrı ön-kayıt) açılır.

**Mühür adayı:** HAYIR. Ön-kaydın kendi tablosunda da mühür adayı yalnız
`YASADI|ODENEBILIR`'dır ve hiçbir okumada oraya varılmamıştır.

---

## 8. "%91.4 KİLİT" DEFTERİNE ETKİSİ

ΔM'nin KESİM/KİLİT defteri (`176_kilit_faturasi_RAPOR.md:495`,
`177_ucuncu_tohum_RAPOR.md:399`) — `KESİM = +0.004929 (%8.6)`,
`KİLİT = +0.052612 (%91.4)`, θ kanalının kilidi `+0.036231` — **sayısal olarak
zerre değişmez**: bu tablo yalnız `son`, `Hkeskin`, `HA4` çapalarından çıkar ve
hiçbir vekil niceliği içine girmez; 178 de tıpkı 177 gibi onu ne değiştirdi ne
değiştirebilirdi. Değişen, tablonun **dayanağıdır**. 176 defterin θ satırına
tek bir nedensel dayanak bırakmıştı (`ω_θ = +0.716`, iki tohum); 177 o dayanağı
kaldırdı (`ω_θ = +1.2235`, dört tohum, bandın dışında, hükümsüz); 178 onu
**geri getiremedi** — ve getiremediği gibi, defterin θ satırının kırılganlığını
görünür kıldı: aynı dört tohumda, aynı gaz koşularında, bit düzeyinde aynı
ölçüm makinesiyle, yalnız `g_X` başka okunduğunda `ω_θ` **−0.1013**'e düşüyor,
yani **işaret bile değiştiriyor** (+0.716 → +1.2235 → −0.1013). Bu bir hüküm
değil bir **tanıdır** — ama defterin ne kadar ince bir zeminde durduğunu
söyler: `θ` satırının ödenebilirliği, `g_X`'in nasıl okunduğuna karşı
sağlam değildir ve 178 bunu para birimi kaymasının kendi büyüklüğüyle
(sınanan etkinin 12 katı) belgelemiştir. Sonuç olarak: **"Gerçeğin fazlasının
%91.4'ü kilittir" cümlesi, 176'da ve 177'de olduğu gibi 178'den sonra da bir
defter satırı olarak durur; nedensel sınavda karşılığı yoktur — ve artık bu
dört tohumla aranamaz** (Ş2.6). Kilidin **varlığı** ve devasalığı (`R_η`
1.265 → 0.808, `Q_E` −%68, `μ̂²_E` −%86) 176'da ölçülmüştü; ona ne 177 ne 178
dokundu, dört tohum onu ayrıca teyit ediyor. Kaybedilen atıf, ölçüm değil.

---

## 9. DENETİM — hakemin yaptığı ve yapmadığı

| | |
|---|---|
| Yeni gaz koşusu | **YOK** — hakem yalnız JSON'ları ve kaynak betikleri okudu |
| Git | **DOKUNULMADI** (commit/add/status değişikliği yok) |
| Eşik gevşetme | **YOK** — hiçbir eşik, tolerans veya kural değiştirilmedi |
| Uydurma sayı | **YOK** — her sayı `ONKAYIT.json` / `HESAP_A.json` / `HESAP_B.json` / kaynak `.py`'den okundu ya da bunlardan hakemin betikleriyle yeniden türetildi |
| Değiştirilen sefer dosyası | **YOK** — `ONKAYIT.json`, `HESAP_*.json` ve `178_configs/*.py` okundu, yazılmadı |
| Hakemin yazdığı dosyalar | bu rapor + `scratchpad/178/hakem_capraz.py` + `scratchpad/178/hakem_mekanizma.py` |
| Hangisinin doğru olduğu seçimi (A mı B mi) | **YAPILMADI** — gerekmedi, ikisi 3e−12 içinde aynı |
| Kurtarma girişimi | **YAPILMADI** — ölümcül kusur karşısında ön-kayıt kuralı gereği yasaktır |

**Bilinen kirlilik (parent'ın bilmesi için):**
`scratchpad/178/log_178b.txt` iki koşunun (A ve B) aynı yola eşzamanlı
yönlendirilmiş çıktısının karışımıdır; B'nin K-KİMLİK^B bölümü üzerine
yazılmıştır. **Hesap dosyaları (`HESAP_A.json`, `HESAP_B.json`) bozulmamıştır**
ve kapı sonuçları dahil her şey onlardadır; hakem hükmünü yalnız bu iki
JSON'dan ve kaynak betiklerden kurmuştur. Log dosyası hükümde kullanılmadı ve
hakem tarafından silinmedi.

---

*Hakem raporu — 178, `/Users/ugursezen/Desktop/arin/deney/qm_riemann/178_gX_hakem_RAPOR.md`*
