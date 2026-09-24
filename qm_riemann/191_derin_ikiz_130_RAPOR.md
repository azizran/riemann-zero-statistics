# 191 — DERİN İKİZ τ≤1.30: ZARF TAM MI KAPANIR, YOKSA GERÇEK BİR KALINTI MI KALIR?

**Kalem (KALEM_DERIN_IKIZ_130_23EYL2026.md, commit `4efc958`, Wed Sep 23
18:31:28 +0300 2026):** 189, derin ikizi τ≤1.10 ve τ≤1.20'ye kadar
kurup zarfın (4. demir, w(τ)=1−0.149·τ^1.30) HAVUZ'da %32 → %51'ini
kapattığını ölçmüştü (H-189a MÜHÜR), ama D→∞ modelleri (A, B) — ikisi de
**dış-değerleme** — zarfın "tamamının" derinlik olduğunu iddia ediyordu
(P_der A %134, B %102). Kaptan, 189'un ÜÇ noktasından (D=1.00 orijin,
1.10, 1.20) veri-öncesi bir **harita-doğrusal yasa** fark etti:
f_D/Δz_D ≈ 3.38–3.39 sabit (Δz_D = 188'in GERÇEK-kinematikli harita
toplamı − z_1.00). Bu yasa **veri görülmeden** τ≤1.30'a **f_1.30 = 0.612
[0.58, 0.64]** ve **σ_ε(1.30) ∈ [0.4140, 0.4160]** öngörüyordu; yasanın
sonsuz-derinlik iması **f_∞ ≈ 0.86** idi (zarfın ~%14'ü derinlik-dışı
gerçek bir kalıntı) — 189'un modelleriyle (~%100) AÇIKÇA ÇELİŞEN bir
öngörü. İnşa (K1, 191b_insa.py) kalem commitinden SONRA, gece koşuldu
(425 863 çizgi) ve kapılardan geçti. Bu görev K2 (zincir) → K3 (hüküm +
kimlik kaydı) → K4 (kilit güç analizi) adımlarını yürütür.

Her sayının yanında onu üreten betik adı vardır. Tek dalga koşuldu.
Ölümler kurtarılmadı. Git'e dokunulmadı. Sonuç ORTAK TEFTİŞE sunulur,
commit sonra.

---

## K0 — ÖN-KAYIT [`191a_onkayit.py`]

`scratchpad/191/ONKAYIT_191.json` yazıldı, **ZİNCİRDEN (K2) ÖNCE**:

- **sha256 = `e8f07b9da9f8c15a…`**
- **damga** `Wed Sep 23 23:45:04 +0300 2026`
- KALEM commit (git log, salt-okunur doğrulama): **`4efc9586…`**,
  `Wed Sep 23 18:31:28 +0300 2026` — inşa (grid-sınavı logu 18:45, inşa
  bitişi 23:36) bu commit'ten SONRA, bu JSON'dan ÖNCE koştu; sıralama
  KALEM'in istediği gibi.
- **makine kimliği K0'da dondu:** sha256(`189c_zincir.py`) ve
  sha256(`189d_hukum.py`) kaydedildi — K2/K3 bu sha'ları yeniden
  hesaplayıp EŞİTLİĞİ doğrular (değişmemiş makine kanıtı, aşağıda).
- **Kaynaktan yeniden hesap (K0'ın kendi tutarlılık denetimi, YENİ bir
  öngörü değil):** 188b'nin `harita_K_gercek.npz` (HAVUZ) ve 187c'nin
  gerçek |ζ_g| değerinden z_D dörtlüsü ve harita-doğrusal yasa sabitleri
  yeniden türetildi:

  | D | z_D (kalem) | z_D (yeniden-hesap) | fark |
  |---|---|---|---|
  | 1.00 | 0.0758 | 0.0758 | 4.3e-05 |
  | 1.10 | 0.1717 | 0.1717 | 2.8e-05 |
  | 1.20 | 0.2261 | 0.2261 | 4.2e-05 |
  | 1.30 | 0.2566 | 0.2566 | 1.9e-05 |
  | tam (187c) | 0.3287 | 0.3287 | — |

  yasa eğimi: 1.10→3.3836 (kalem 3.38), 1.20→3.3881 (kalem 3.39),
  merkez 3.3858 (kalem 3.385) → **f_1.30 kaynaktan = 0.6119** (kalem
  0.612). **Dört sayı da EŞİT** (±5e-5, ±0.01). Ön-mühür tıpatıp
  doğrulandı.

---

## K1 — İNŞA (YAPILDI, tekrarlanmadı) [`191b_insa.py`, kaptan]

`scratchpad/189/insa_kapi_Hderin130.json`: **GEÇTİ** — maks|F| =
1.863e-09, sıralılık TAM (min Δz = 0.0905), L = 12.02959324, 425 863
çizgi, ızgara sınavı 0/5000 fark (h=0.015 seçildi). Bu görev bu adımı
tekrar ETMEDİ.

---

## K2 — ZİNCİR [`191c_zincir.py`]

189c_zincir.py (ve 189d_hukum.py) **DOSYA OLARAK DEĞİŞTİRİLMEDEN**
importlib ile yüklendi; tek yeni şey `b189c.zincir("Hderin130")`
çağrısı (184b2 → 185b → 186b → 187c). Mevcut Hkeskin/Hderin110/
Hderin120/gerçek dosyalarının HİÇBİRİ üzerine yazılmadı — `zincir()` ve
`gercek_zeta()` fonksiyonları Hkeskin/gerçek için TEKRAR ÇAĞRILMADI
(deterministik/değişmemiş kod olduğu için tekrarı bilgi katmazdı, yalnız
üzerine-yazma riski katardı); bunun yerine 189'un zaten geçmiş sonucu
referans alındı.

**MAKİNE MÜHRÜ:**

| kontrol | sonuç |
|---|---|
| sha256(189c_zincir.py) K0 ile bit-bit eşit | **EVET** |
| sha256(189d_hukum.py) K0 ile bit-bit eşit | **EVET** |
| 189/muhur_Hkeskin.json (184-187, 6/6, bit-bit; TEKRAR KOŞULMADI) | **6/6 TUTTU** (189'dan referans) |
| 189/zincir_gercek.npz mevcut (TEKRAR KOŞULMADI) | **EVET** |

Aynı, değişmemiş sarmalayıcı zinciri → Hderin130 için 184b2→187c koştu
(110 s). Çıktı: `scratchpad/189/{K1,OZ,G1_proj,zincir}_Hderin130.*`
(YENİ dosyalar), `scratchpad/191/ZINCIR_191_ozet.json`.

**HAVUZ (191c ekran defteri):** r_1.30 = 0.9665±0.0044, σ_ε(1.30) =
0.41564±0.00009, |ζ_H130| = 0.2334±0.0019.

---

## BANT TABLOSU — ölçülen (dört derinlik) [`191d_hukum.py`]

**r_D = w_g/w_HD (± σ_r, 184 konvansiyonu):**

| bant | r_1.00 | r_1.10 | r_1.20 | r_1.30 | f_1.30 |
|---|---|---|---|---|---|
| 0.45-0.50 | 0.9445 | 0.9678 | 0.9819 | 0.9892±.0015 | 0.806±.014 |
| 0.50-0.55 | 0.9351 | 0.9608 | 0.9756 | 0.9850±.0013 | 0.768±.010 |
| 0.55-0.60 | 0.9253 | 0.9522 | 0.9697 | 0.9803±.0033 | 0.737±.011 |
| 0.60-0.65 | 0.9169 | 0.9465 | 0.9646 | 0.9763±.0032 | 0.714±.011 |
| 0.65-0.70 | 0.9096 | 0.9388 | 0.9575 | 0.9692±.0039 | 0.659±.008 |
| 0.70-0.75 | 0.9021 | 0.9327 | 0.9513 | 0.9631±.0040 | 0.623±.012 |
| 0.75-0.80 | 0.8978 | 0.9293 | 0.9458 | 0.9580±.0078 | 0.589±.005 |
| 0.80-0.86 | 0.8931 | 0.9237 | 0.9382 | 0.9481±.0131 | 0.515±.012 |
| **HAVUZ** | **0.9096** | **0.9389** | **0.9556** | **0.9665±.0044** | **0.629±.008** |

f_1.10 = 0.3244±0.0105, f_1.20 = 0.5092±0.0092, f_1.30 = **0.6295±0.0080**
(HAVUZ; se ortak-blok loo jackknife). Kesime doğru kapanış payı düşük-τ
bantlarında en yüksek (%81 @ 0.45-0.50), kesimde en düşük (%52 @
0.80-0.86) — 189'un gözlemlediği eğilim sürüyor.

---

## HÜKÜM — H-191a/b/c [`191d_hukum.py`]

### H-191a — ZARF SÜRER (birincil, yön): **MÜHÜR**

- Δf(1.30−1.20) = **+0.1203 ± 0.0026** (ortak-blok loo) → HAVUZ yönü
  **>2se** (46× se) — açık ara ile geçti.
- 8/8 bantta r_1.20 ≤ r_1.30 + 2σ (184 konv., bağımsız yayılım).
- Ölüm koşulu (f_1.30 ≤ f_1.20+2se) **ateşlemedi**.

### H-191b — HARİTA-DOĞRUSAL YASA (nicel): **MÜHÜR**

- **f_1.30 = 0.6295 ± 0.0080**, bant **[0.58, 0.64]** — bandın **İÇİNDE**.
- Ön-mühür merkezinin (0.612) **+2.17 se üstünde** — bandı aşmadı ama
  merkeze göre yüksek uçta; kapanış yasanın öngördüğünden biraz daha
  hızlı (aşağıda K3'te bu, yeniden-fit edilen eğimin 3.385'ten 3.44'e
  çıkmasıyla aynı gözlem).
- Yasa üç derinlikte (1.10, 1.20, 1.30) tutarlı MÜHÜR aldığı için
  f_∞ ≈ 0.86 okuması (K3'te doğrudan ölçülür) **güç kazanıyor**.

### H-191c — σ_ε YAKINSAMASI: **MÜHÜR**

- σ_ε: Hk 0.43135 → 110 0.42248 → 120 0.41788 → **130 0.41564±0.00009**
  (gerçek 0.40921±0.00011).
- Δ(120−130) = **+0.002235±0.000076** → azalıyor, **29.5σ**.
- 130 − gerçek = **+0.006434±0.000099** → **aşma yok, 65σ** güvenle.
- Nicel bant **[0.4140, 0.4160]** içinde (0.41564, bandın üst-orta
  kesiminde).
- Gerçeğe doğru kapanan pay: **%70.9** (110: %40.1, 120: %60.8, 130:
  %70.9 — artış hızı yavaşlıyor, beklenen asimptotik yaklaşma).

**Üç hipotezin ÜÇÜ DE MÜHÜR.** Ön-mühür bandının hem yönü hem nicel
aralığı — veri görülmeden 12 saat önce donmuş — tutuldu.

---

## K3 — KİMLİK KAYDI (dört derinlik) [`191d_hukum.py`]

### r_D(τ) = 1 − c·τ^α (184c makinesi AYNEN)

| D | c | α | χ²/dof |
|---|---|---|---|
| 1.00 | 0.1492 ± 0.0080 / ± 0.0026 | 1.304 ± 0.092 / ± 0.044 | 0.343 |
| 1.10 | 0.1155 ± 0.0084 / ± 0.0040 | 1.681 ± 0.131 / ± 0.075 | 0.291 |
| 1.20 | 0.0997 ± 0.0105 / ± 0.0038 | 2.225 ± 0.193 / ± 0.083 | 0.140 |
| 1.30 | 0.0899 ± 0.0131 / ± 0.0036 | 2.802 ± 0.281 / ± 0.111 | **0.045** |

Trend sürüyor ve GÜÇLENİYOR: c 0.149→0.100→0.090 (küçülüyor), α
1.30→2.22→2.80 (büyüyor, hızlanarak), χ²/dof 0.34→0.14→**0.045** (fit
her derinlikte DAHA İYİ oturuyor — zarf kesime çekildikçe basit
güç-yasası formu daha az gerilimle uyuyor).

### f_D ve g_D'nin Δz_D'ye karşı doğrusallığı (dört nokta)

D=1.00 orijin (f≡g≡0 TANIM gereği, ölçüm değil) hariç, ÜÇ ölçülen
noktayla (1.10, 1.20, 1.30):

| nicelik | orijinden-geçen eğim ± se | serbest eğim ± se | serbest kesişim ± se |
|---|---|---|---|
| f_D | **3.443 ± 0.034** | 3.593 ± 0.156 | −0.023 ± 0.024 |
| g_D | **4.018 ± 0.090** | 3.674 ± 0.359 | +0.050 ± 0.051 |

Serbest fitin kesişimi ikisinde de **sıfırla tutarlı** (~1σ) — yasa
(orijinden geçme varsayımı) üç derinlikte de makul kalıyor; ama
orijinden-geçen eğim (3.443), kalemin veri-öncesi tahmininden (yalnız
1.10/1.20 çiftinden, 3.385) hafifçe **yukarı kaydı** — f_1.30'un bant
merkezinin üstünde çıkmasıyla aynı sinyal.

**f_∞ = eğim(orijin) × Δz_∞(0.2528) = 0.8706 ± 0.0086** ⇒ **zarfın
derinlik payı %87.1 ± %0.9, kalıntı payı %12.9 ± %0.9.**
**g_∞ = 1.016 ± 0.023** (σ_ε aralığı derinlikle gerçeğe pratik olarak
TAM yakınsıyor — ön-mühürün g_∞∈[1.02,1.06] imasıyla uyumlu, ufak
tefek altında).

### 189'un D→∞ modelleriyle yan yana

| okuma | P_der (HAVUZ) | doğa |
|---|---|---|
| 189 Model A (kinematik 1.20'de donuk + Γ_g) | **1.335 ± 0.008** (%134) | dış-değerleme, TEK sıçrama |
| 189 Model B (kinematik eğilimi Γ'ye doğrusal) | **1.015 ± 0.003** (%102) | dış-değerleme, TEK sıçrama |
| **191 harita-doğrusal yasa** (3 ölçülen nokta, 1.10/1.20/1.30) | **0.871 ± 0.009** (%87) | **iç-değerleme + kısa dış-değerleme**, üç bağımsız derinlikte sınandı |

**Bu üç okuma ANLAŞMIYOR** ve 191'in okuması ikisinin de ALTINDA: 189'un
A/B modelleri sıfır veya bir veri noktasından (D=1.20) kinematiği
dondurup tek sıçramayla D→∞'a taşıyordu; 191'in yasası ÜÇ bağımsız
derinlikte ÖLÇÜLEN f_D'nin kendisinden geliyor ve zarfın büyük
çoğunluğunun (~%87) derinlikle kapandığını ama gerçek, ölçülebilir bir
**~%13'lük kalıntının** (K0'ın veri-öncesi "f_∞≈0.86" imasına çok
yakın) ayakta kaldığını söylüyor. 189'un "zarfın tamamı derinlik"
okuması (A/B) bu üç-noktalı ampirik yasayla **DOĞRULANMADI**.

---

## K4 — KİLİT GÜÇ ANALİZİ (ölçümden ÖNCE) [`191e_guc.py`]

180'in HAM gerçek↔Hkeskin farkları (log konvansiyonu): Δlog M =
+0.057542, Δlog Q_E = −0.044475; kilit payı s_kilit: M %18.5±11.3,
Q_E %20.6±7.0. İki rakip okumanın (Δ(D)=Δ(1.00)(1−f_D) "kilit de
derinlik" vs Δ(D)=Δ(1.00)[(1−s_kilit)(1−f_D)+s_kilit] "kilit tabanı")
ayrımı, cebirsel özdeşlik **Ayrım(D) = Δ(1.00)·s_kilit·f_D** ile
D=1.30'da, ölçülen f_1.30=0.6295±0.0080 kullanılarak:

| para birimi | Ayrım | σ_birleşik | **Ayrım/σ_birleşik** | yeterli mi (≥3) |
|---|---|---|---|---|
| M (log) | +0.006701 | 0.004094 | **1.64** | HAYIR |
| Q_E (log) | −0.005767 | 0.001961 | **2.94** | HAYIR (eşiğe yakın) |

**Sınırlama (dürüstçe kaydedilir):** Δ(1.00)'ın kendi jackknife se'si
180 raporunda YOK (180a'nın formülü: ΔN "tohumsuz"), ve kaynak
scratchpad (166–183) bu oturumda temizlenmiş olduğu için yeniden
hesaplanamadı; σ_Δ(1.00)=0 varsayıldı. Bu varsayım σ_birleşik'i
KÜÇÜLTÜR (oranı büyütür) — yani gerçek oran muhtemelen 1.64 ve 2.94'ten
DAHA DA DÜŞÜK olurdu. Sonuç bu sınırlamadan zarar görmüyor, tersine
sağlamlaşıyor.

**KARAR: ÖLÇME YAPILMAZ — "güç yetersiz."** Ne M ne Q_E 3σ eşiğini
geçiyor (M 1.64σ, Q_E 2.94σ). İkincil kapı (180'in `167_olcum→172b`
hattının 60 dakikada yeni bir gaz için koşabilmesi) bilgi amaçlı
denetlendi: **scratchpad/166–183'ün TAMAMI bu oturumda temizlenmiş**
(167/C_*.json, 172b'nin band-defterleri, 176'nın VF vekilleri — hiçbiri
yok); bu hat şu an yeniden kurulmadan koşamaz, dolayısıyla ikinci kapı
da (varsayımsal olarak) kapanırdı. Karar tek başına GÜÇ testinden
geliyor; altyapı durumu yalnız destekleyici bilgi.

---

## TÜRETİLEN vs ÖLÇÜLEN

- **Türetilen (veri-öncesi, KALEM/K0):** harita-doğrusal yasa (yalnız
  189'un 1.10/1.20 çiftinden), f_1.30 = 0.612 [0.58,0.64], σ_ε(1.30) ∈
  [0.4140,0.4160], f_∞≈0.86 iması, K4'ün cebirsel özdeşliği ve karar
  kuralı.
- **Ölçülen (191c/191d):** r_D (4 derinlik, 8 bant + HAVUZ), f_1.30 =
  0.6295±0.0080, σ_ε(1.30) = 0.41564±0.00009, (c,α)(D) dört derinlikte,
  f_D/g_D'nin Δz_D'ye karşı doğrusal eğimi (3 nokta) ve bunun
  ekstrapolasyonu f_∞ = 0.871±0.009, K4'ün Ayrım/σ_birleşik oranları.
- **Doğrulanan:** yön (H-191a), nicel bant (H-191b, merkeze göre yüksek
  uçta ama içeride), σ_ε yakınsaması (H-191c) — **üçü de MÜHÜR**; yasa
  üç derinlikte tutarlı.
- **Kısmen sapan (kayıt, ölüm değil):** ölçülen f_1.30 ve buradan
  yeniden-fit edilen eğim (3.44), kalemin 2-noktalı tahmininin (3.385)
  hafif üstünde — kapanış, yasanın en saf halinden biraz daha hızlı.
  Bu f_∞'u da hafifçe yukarı taşıyor (0.871 vs kalemin informal "~0.86"
  iması) ama İKİSİ DE 189'un A/B modellerinin (1.335 / 1.015)
  ALTINDA kalıyor — 191'in temel bulgusu (gerçek bir derinlik-dışı
  kalıntı var) değişmiyor.
- **Açık kalanlar:** α'nın derinlikle BÜYÜME yasası (1.30→2.80,
  hızlanarak) hâlâ türetilmedi; K4'ün "güç yetersiz" sonucu bu görevi
  kapatıyor — kilit payının derin ikizlerde DOĞRUDAN ölçümü (180 hattının
  yeniden kurulmasını gerektirir) borç olarak kalıyor.

---

## MANŞET (aday cümle)

> **Zarf çoğunlukla derinliktir, ama gerçek bir kalıntı da vardır.**
> Kaptanın 189'un yalnız İKİ noktasından (τ≤1.10, τ≤1.20) veri-öncesi
> fark ettiği harita-doğrusal yasa (f_D/Δz_D ≈ 3.39 sabit), üçüncü,
> bağımsız bir derinlikte (τ≤1.30) SINANDI ve DOĞRULANDI: f_1.30 =
> **0.630 ± 0.008**, ön-mühürlü bandın [0.58,0.64] İÇİNDE (H-191b
> MÜHÜR); zarf sürmeye devam ediyor, 8/8 bantta (H-191a MÜHÜR); σ_ε
> öngörülen bandın [0.4140,0.4160] içinde gerçeğe yaklaşmaya devam etti,
> hâlâ aşmadan (H-191c MÜHÜR). Üç derinlikten doğrudan ölçülen eğim
> ekstrapole edildiğinde **f_∞ = 0.871 ± 0.009**: zarfın **%87'si**
> merdiven derinliğinin bir eseri, ama **%13'ü** derinlikle KAPANMAYAN,
> gerçek bir kinematik-fark kalıntısıdır. Bu, 189'un D→∞ modellerinin
> (A: %134, B: %102 — "zarfın TAMAMI derinlik") **İKİSİNİ DE
> desteklemiyor**: onlar tek bir dondurulmuş kinematik noktadan (D=1.20)
> sıçrayan dış-değerlemelerdi, 191'in yasası ise üç bağımsız, ölçülmüş
> derinlikten gelen bir iç-değerleme + kısa bir dış-değerleme. Kilit
> payının bu kalıntıya katkısını doğrudan ölçmek için gereken güç
> analizi (K4) YETERSİZ çıktı (M 1.64σ, Q_E 2.94σ, eşik 3σ) — ham para
> birimi ölçümü YAPILMADI, "güç yetersiz" diye kaydedildi.

---

Teslim:
- **Rapor:** bu dosya.
- **Betikler** (`191_configs/`):
  - `191a_onkayit.py`: K0
  - `191b_insa.py`: (K1, önceden koşuldu) ızgara sınavı + inşa sarmalayıcısı
  - `191c_zincir.py`: K2 — 189c ince sarmalayıcı (Hderin130)
  - `191d_hukum.py`: K3 — H-191a/b/c + kimlik kaydı
  - `191e_guc.py`: K4 — kilit güç analizi
  - `191f_figur.py`: figür
- **Figür:** `191_derin_ikiz_130.png`. Sol: r_D(τ) dört derinlik + 1−cτ^α
  fitleri. Orta: üst f_D, alt g_D — Δz_D'ye karşı dört nokta + yasa
  çizgisi + ∞ imaları. Sağ: üst σ_ε(D) + ön-mühür bandı, alt c(D)/α(D)
  (twin eksen).
- **`scratchpad/191/`:**
  - ONKAYIT_191.json, ZINCIR_191_ozet.json, HUKUM_191.json, GUC_191.json
  - log_izgara.txt, log_insa.txt (K1, önceden), log_zincir.txt,
    log_hukum.txt, log_guc.txt
- **`scratchpad/189/`:** `{K1,OZ,G1_proj,zincir}_Hderin130.{npz,json}`
  (YENİ dosyalar; mevcut Hkeskin/Hderin110/Hderin120/gerçek dosyalarının
  hiçbiri değiştirilmedi).
- **`scratchpad/164/`, `scratchpad/155/`:** `z_Hderin130.npy`,
  `insa_Hderin130.json` (K1, önceden koşuldu).
