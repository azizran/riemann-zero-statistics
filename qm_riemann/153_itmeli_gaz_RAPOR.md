# 153 — itmeli gaz: seviye itmesi anormal fazı gerçeğe indiriyor mu?

152'nin ölçüm zinciri **birebir** korunarak 15 koşu yapıldı: taban kontrolü (152'nin A4'ü), görevin istediği iki itme yolu (İ1, İ2) ve tarama sırasında zorunlu hale gelen iki ek seri (İ1b ve J — gerekçeleri aşağıda). Soru: saf-merdiven gazının gerçeğin ~1.4 katı olan anormal dispersiyon fazını GUE seviye itmesi kapatıyor mu?

> **Kısa hüküm:** hayır — ve sebebi görevin varsaydığının tersi. Sentetik gaz seviye itmesinden yoksun DEĞİL; kısa menzilde gerçeğin ~190 katı KATI (P(s<0.3): 0.00013 vs 0.02420). Görevde yazıldığı haliyle itme (ε>0) fazı 1.72'den 2.24'e **kötüleştiriyor**. Kısa menzili YUMUŞATMAK fazı kısmen indiriyor (1.72 → ~1.42) ama orada doyuyor ve S3'ü bozuyor. İ2/GUE tabanı biraz daha iyi (1.36) ama kendi kurgusundan gelen bir yapaylık taşıyor. **Hiçbiri 1.00'e inmiyor.**

**Koşulan kod:** `153_configs/153_gaz.py` (152 gibi tek kod yolu; zincir `g = np.diff(z)`den itibaren 152'nin karakteri karakterine aynısıdır). Kontrol bunu kanıtlıyor: `taban`, 152'nin A4'ünün her Newton iterasyonunda aynı maks|F|'yi ve beş bandın dördüncü basamağına kadar aynı Γ'sını üretti.


## 1. Kalibrasyon hedefi — ve taramayı baştan çeviren ölçüm

Görev, sentetik P(s)'in küçük-s kuyruğunu gerçeğinkine yaklaştırmayı kalibrasyon hedefi koydu. Hedef önce ÖLÇÜLDÜ (zeros6 son-300k, zincirin kendi açılması `s = g·log(mid/2π)/2π`), sonra taban gazında aynı büyüklüğe bakıldı:

| | P(s<0.1) | P(s<0.2) | P(s<0.3) | P(s<0.5) | σ_ds² |
|---|---|---|---|---|---|
| **gerçek (zeros6 son-300k)** | 0.00096 | 0.00733 | **0.02420** | 0.10385 | 0.1674 |
| taban gazı (erfc-0.68) | 0.00000 | 0.00000 | **0.00013** | 0.03550 | 0.1128 |
| GUE Wigner surmise (teorik) | 0.00107 | 0.00839 | 0.02725 | 0.11200 | 0.17810 |
| Poisson (teorik) | 0.09516 | 0.18127 | 0.25918 | 0.39347 | 1.0 |

Bu tablo taramanın yönünü değiştirdi. Sentetik gazda 300 bin aralığın **hiçbiri** 0.2ḡ'nin altında değil; gerçekte 2198 tanesi var. Yani sentetik gaz seviye itmesinden YOKSUN değil — kısa menzilde gerçeğin ~190 katı KATI. Birinci-mertebe merdiven gazı pürüzsüz ve deterministik bir haritadır; aralık dağılımı sınırlı bir fonksiyonun dağılımıdır, iki ucu da ince. Gerçek GUE'nin küçük-s kuyruğu bunun çok üstünde.

Sonuç: görevde yazıldığı haliyle itme (ε>0) kalibrasyon hedefinden UZAKLAŞTIRIR. Bu yüzden ε'nun **iki işareti de** koşuldu ve tabloda ayrı gruplar olarak duruyor — İ1a (ε>0, görevin harfi) ve İ1b (ε<0, görevin kalibrasyon hedefi). Bu tek sapmadır, sebebi de yukarıdaki ölçümdür.


## 2. Konfigürasyonlar

| konfig | grup | değişiklik | küçük-s: P(s<0.3) | σ_ds² | σ_η² | c₁ |
|---|---|---|---|---|---|---|
| **gerçek** | — | — | **0.02420** | **0.1674** | **0.0227** | **−0.01158** |
| taban erfc-0.68 (kontrol) | kontrol | değişiklik yok — 152'nin A4'ü, öz-tutarlı Newton | 0.00013 | 0.1128 | 0.0230 | -0.00918 |
| R1  İ1 ε=+0.003 ×5 | İ1a | itme (görevde yazıldığı gibi), 5 süpürme | 0.00000 | 0.0945 | 0.0200 | -0.00710 |
| R2  İ1 ε=+0.010 ×5 | İ1a | itme, 5 süpürme | 0.00000 | 0.0674 | 0.0159 | -0.00437 |
| R3  İ1 ε=+0.030 ×5 | İ1a | itme, 5 süpürme | 0.00000 | 0.0328 | 0.0091 | -0.00095 |
| N3  İ1 ε=−0.0070 ×3 | İ1b | yumuşatma, P(s<0.3)'e kalibre, 3 süpürme | 0.02539 | 0.1559 | 0.0337 | -0.01653 |
| N5  İ1 ε=−0.0037 ×5 | İ1b | yumuşatma, P(s<0.3)'e kalibre, 5 süpürme | 0.02324 | 0.1525 | 0.0330 | -0.01619 |
| N10 İ1 ε=−0.0018 ×10 | İ1b | yumuşatma, P(s<0.3)'e kalibre, 10 süpürme | 0.02426 | 0.1538 | 0.0338 | -0.01682 |
| N5x İ1 ε=−0.010 ×5 | İ1b | yumuşatma, kalibre dozun ~2.7 katı | 0.09516 | 0.2485 | 0.0675 | -0.03768 |
| N5y İ1 ε=−0.020 ×5 | İ1b | yumuşatma, ~5.4 kat | 0.16887 | 0.3672 | 0.1118 | -0.06020 |
| N5z İ1 ε=−0.030 ×5 | İ1b | yumuşatma, ~8.1 kat | 0.21955 | 0.4584 | 0.1464 | -0.07426 |
| J14 titreşim σ_j=0.141 | J | yapısız titreşim — σ_ds² N5 ile eşleşir | 0.02196 | 0.1527 | 0.0663 | -0.03009 |
| J26 titreşim σ_j=0.260 | J | yapısız titreşim — σ_ds² N5x ile eşleşir | 0.07170 | 0.2337 | 0.1570 | -0.07194 |
| Pd  İ2 taban=düzgün | İ2 | s≡1 taban + tek geçiş merdiven boyası | 0.05473 | 0.1627 | 0.0874 | -0.05636 |
| P1  İ2 taban=GUE | İ2 | GUE-surmise gap tabanı + tek geçiş boya | 0.10752 | 0.3367 | 0.2118 | -0.02410 |
| P0  İ2 taban=Poisson | İ2 | üstel (Poisson) gap tabanı + tek geçiş boya | 0.30972 | 1.1216 | 0.9854 | -0.00070 |

İ1'de yer değiştirme `Δz_n = ε·ḡ·[(ḡ/g_{n−1})² − (ḡ/g_n)²]`, ±0.4·min(g_{n−1},g_n) kelepçeli. Kelepçe sıralamanın bozulmamasını GARANTİ eder (komşular en çok 0.4g yaklaşabildiğinden yeni aralık ≥0.2g): dokuz İ1 koşusunun hiçbirinde sıra bozulmadı (sıra bozan çift = 0). Pencere ortalama-yoğunluğu (30 pencere) en kötü koşuda 6.4e−05·ḡ kaydı — z yeniden ölçeklenmedi/ötelenmedi, hareketin teleskopik olması yetti. J koşuları (yapısız titreşim) 152'nin H-B'siyle aynıdır ve sırayı BOZAR (J14: 245 çift %0.08, J26: 5130 çift %1.71); 152'de olduğu gibi sonrasında sıralanır.

**Kelepçenin bağlama oranı doz arttıkça patlıyor** ve bu, aşağıdaki doygunluğun büyük ihtimalle sebebidir — son süpürmede bağlanan hamlelerin oranı: R1 %0.00, R2 %0.01, R3 %1.05, N3 %1.59, N5 %1.66, N10 %2.09, N5x %13.5, **N5y %29.3, N5z %39.2**. Yani N5y ve N5z artık saf 1/s² gevşetmesi değil, hareketlerin üçte birinde kelepçenin belirlediği bir dinamiktir; bu ikisinden fizik hükmü çıkarılmamalıdır, yalnızca 'doz büyütmek daha fazla kazanç getirmiyor' gözlemi için duruyorlar.


## 3. Bantlar — Γ_rot(Re, Im), |Γ|, R

Hücre: `(Re,Im) |Γ| R`. Kurallar 152 ile AYNI: `‡` = |Γ|>1.5, bant sayısal olarak patlamıştır (ayrışım paydası sıfıra gitmiş), ölçüm değildir, hiçbir ortalamaya girmez. `*` = R uydurması fazı tutturamadı (faz-artığı > 0.05 rad).

| τ̄ | gerçek | taban | R1 | R2 | R3 | N3 | N5 | N10 | N5x | N5y | N5z | J14 | J26 | Pd | P1 | P0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.5375 | (+0.787,+0.206) 0.81 0.40 | (+0.600,+0.406) 0.72 1.23 | (+0.581,+0.429) 0.72 1.53 | (+0.551,+0.474) 0.73 2.60 | (+0.497,+0.591) 0.77 (3.00)* | (+0.640,+0.362) 0.74 0.83 | (+0.641,+0.363) 0.74 0.85 | (+0.645,+0.361) 0.74 0.83 | (+0.657,+0.319) 0.73 0.55 | (+0.627,+0.304) 0.70 0.43 | (+0.598,+0.304) 0.67 0.40 | (+0.592,+0.391) 0.71 0.98 | (+0.570,+0.381) 0.69 0.70 | (+1.076,+0.971) 1.45 1.25 | (+0.419,+0.180) 0.46 0.20 | (+0.217,-0.075) 0.23 (0.00)* |
| 0.585 | (+0.550,+0.488) 0.74 1.07 | (+0.322,+0.619) 0.70 2.02 | (+0.290,+0.637) 0.70 2.78 | (+0.234,+0.677) 0.72 (3.00)* | (+0.118,+0.780) 0.79 (3.00)* | (+0.382,+0.587) 0.70 1.33 | (+0.382,+0.589) 0.70 1.35 | (+0.386,+0.588) 0.70 1.33 | (+0.422,+0.540) 0.68 0.85 | (+0.404,+0.503) 0.64 0.62 | (+0.382,+0.481) 0.61 0.55 | (+0.320,+0.612) 0.69 1.55 | (+0.310,+0.603) 0.68 1.10 | (+0.67,-2.28) **‡** | (+0.244,+0.318) 0.40 0.45 | (+0.114,-0.085) 0.14 (3.00)* |
| 0.66 | (+0.118,+0.614) 0.63 1.85 | (-0.183,+0.714) 0.74 (3.00)* | (-0.221,+0.721) 0.75 (3.00)* | (-0.290,+0.762) 0.82 (3.00)* | (-0.346,+1.066) 1.12 (3.00)* | (-0.118,+0.716) 0.73 1.98 | (-0.118,+0.718) 0.73 2.02 | (-0.115,+0.720) 0.73 2.00 | (-0.096,+0.693) 0.70 1.20 | (-0.128,+0.632) 0.64 0.85 | (-0.127,+0.580) 0.59 0.70 | (-0.218,+0.744) 0.78 2.33 | (-0.344,+1.029) 1.08 1.58 | (-0.062,+0.670) 0.67 1.90 | (-0.063,+0.406) 0.41 0.75 | (+0.339,+0.149) 0.37 (0.00)* |
| 0.74 | (-0.148,+0.498) 0.52 2.12 | (-1.84,-1.37) **‡** | (+0.258,+1.298) 1.32 1.98 | (-0.448,+0.521) 0.69 (3.00)* | (-0.745,+0.257) 0.79 (3.00)* | (-0.566,-0.005) 0.57 (0.00)* | (-0.584,-0.034) 0.58 (0.00)* | (-0.558,-0.020) 0.56 (0.00)* | (-0.137,+0.357) 0.38 1.05 | (+0.025,+0.481) 0.48 0.60 | (+0.104,+0.507) 0.52 0.50 | (+0.033,-0.046) 0.06 (0.00)* | (+0.297,+0.205) 0.36 0.38 | (-3.97,-1.20) **‡** | (-0.050,+1.341) 1.34 0.50 | (+0.604,+0.616) 0.86 (0.00)* |
| 0.815 | (-0.103,+0.655) 0.66 1.50 | (-1.259,-0.557) 1.38 (0.00)* | (+5.19,+4.29) **‡** | (-0.013,+0.362) 0.36 (3.00)* | (-0.411,+0.181) 0.45 (3.00)* | (-0.523,+0.043) 0.52 (3.00)* | (-0.538,+0.023) 0.54 (3.00)* | (-0.518,+0.042) 0.52 (3.00)* | (-0.178,+0.439) 0.47 0.85 | (-0.033,+0.584) 0.59 (0.53)* | (+0.090,+0.678) 0.68 (0.48)* | (+0.471,+0.569) 0.74 0.62 | (+2.98,+1.17) **‡** | (+1.84,+4.49) **‡** | (+1.097,-0.792) 1.35 (0.00)* | (+0.49,-4.23) **‡** |

## 4. Anormal fazın dikliği — arg Γ_rot(sentetik) / arg Γ_rot(gerçek)

Hedef 1.000. Kapanış ölçütü budur.

`ort(b1–4)` 152'nin sütunuyla doğrudan karşılaştırılabilir olsun diye aynı kuralla hesaplandı. Yanına `ort(b1–3)` eklendi: 152, bant 4 ve 5'in her koşuda gürültülü olduğunu gösterdi, ve aşağıda görüleceği gibi bu taramadaki iyileşmenin büyük kısmı tam o iki bantta oturuyor — ayrı okunabilmeli. `b3` ise tek başına en sağlam banttır (220 aday, her koşuda R-artığı küçük).

| konfig | τ=0.5375 | τ=0.585 | τ=0.66 | τ=0.74 | τ=0.815 | **ort(b1–4)** | ort(b1–3) | b3 |
|---|---|---|---|---|---|---|---|---|
| **gerçek** | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | **1.00** | 1.00 | 1.00 |
| taban erfc-0.68 (kontrol) | 2.32 | 1.50 | 1.32 | ‡ | 2.06 | **1.72** (3b) | 1.72 (3b) | 1.32 |
| R1  İ1 ε=+0.003 ×5 | 2.48 | 1.58 | 1.35 | 0.74 | ‡ | **1.54** (4b) | 1.80 (3b) | 1.35 |
| R2  İ1 ε=+0.010 ×5 | 2.78 | 1.71 | 1.40 | 1.23 | 0.93 | **1.78** (4b) | 1.96 (3b) | 1.40 |
| R3  İ1 ε=+0.030 ×5 | 3.41 | 1.96 | 1.36 | 1.51 | 1.58 | **2.06** (4b) | 2.24 (3b) | 1.36 |
| N3  İ1 ε=−0.0070 ×3 | 2.01 | 1.37 | 1.26 | 1.69 | 1.77 | **1.58** (4b) | 1.55 (3b) | 1.26 |
| N5  İ1 ε=−0.0037 ×5 | 2.01 | 1.37 | 1.26 | 1.72 | 1.80 | **1.59** (4b) | 1.55 (3b) | 1.26 |
| N10 İ1 ε=−0.0018 ×10 | 1.99 | 1.36 | 1.25 | 1.71 | 1.77 | **1.58** (4b) | 1.54 (3b) | 1.25 |
| N5x İ1 ε=−0.010 ×5 | 1.76 | 1.25 | 1.24 | 1.04 | 1.13 | **1.32** (4b) | 1.42 (3b) | 1.24 |
| N5y İ1 ε=−0.020 ×5 | 1.76 | 1.23 | 1.28 | 0.82 | 0.94 | **1.27** (4b) | 1.43 (3b) | 1.28 |
| N5z İ1 ε=−0.030 ×5 | 1.84 | 1.24 | 1.29 | 0.74 | 0.83 | **1.28** (4b) | 1.46 (3b) | 1.29 |
| J14 titreşim σ_j=0.141 | 2.28 | 1.50 | 1.34 | -0.51 | 0.51 | **1.15** (4b) | 1.71 (3b) | 1.34 |
| J26 titreşim σ_j=0.260 | 2.30 | 1.51 | 1.37 | 0.32 | ‡ | **1.38** (4b) | 1.73 (3b) | 1.37 |
| Pd  İ2 taban=düzgün | 2.87 | ‡ | 1.20 | ‡ | ‡ | **2.04** (2b) | 2.04 (2b) | 1.20 |
| P1  İ2 taban=GUE | 1.58 | 1.26 | 1.25 | 0.86 | -0.36 | **1.24** (4b) | 1.36 (3b) | 1.25 |
| P0  İ2 taban=Poisson | -1.30 | -0.88 | 0.30 | 0.43 | ‡ | **-0.36** (4b) | -0.63 (3b) | 0.30 |

## 5. Skor

`s_Γ` = bant BAŞINA ortalama (|ΔRe| + |ΔIm|), bant 1–4 içinden patlak olmayanlar üzerinden — 152'deki düzeltmenin aynısı (toplam değil ortalama; yoksa en çok bandı patlayan konfigürasyon sahte biçimde 'en iyi' görünür). Kaç bant üzerinden alındığı parantezde. Küçük = gerçeğe yakın.

| konfig | s_Γ (bant başına) | faz oranı b1–4 | faz oranı b1–3 | s_R | patlak ‡ |
|---|---|---|---|---|---|
| taban erfc-0.68 (kontrol) | 0.382 (3b) | 1.72 | 1.72 | 0.89 (2b) | 0.74 |
| R1  İ1 ε=+0.003 ×5 | 0.623 (4b) | 1.54 | 1.80 | 0.99 (3b) | 0.815 |
| R2  İ1 ε=+0.010 ×5 | 0.472 (4b) | 1.78 | 1.96 | 2.20 (1b) | — |
| R3  İ1 ε=+0.030 ×5 | 0.788 (4b) | 2.06 | 2.24 | — | — |
| N3  İ1 ε=−0.0070 ×3 | 0.457 (4b) | 1.58 | 1.55 | 0.27 (3b) | — |
| N5  İ1 ε=−0.0037 ×5 | 0.470 (4b) | 1.59 | 1.55 | 0.30 (3b) | — |
| N10 İ1 ε=−0.0018 ×10 | 0.457 (4b) | 1.58 | 1.54 | 0.28 (3b) | — |
| N5x İ1 ε=−0.010 ×5 | 0.217 (4b) | 1.32 | 1.42 | 0.55 (5b) | — |
| N5y İ1 ε=−0.020 ×5 | 0.218 (4b) | 1.27 | 1.43 | 0.75 (4b) | — |
| N5z İ1 ε=−0.030 ×5 | 0.251 (4b) | 1.28 | 1.46 | 0.82 (4b) | — |
| J14 titreşim σ_j=0.141 | 0.481 (4b) | 1.15 | 1.71 | 0.60 (4b) | — |
| J26 titreşim σ_j=0.260 | 0.591 (4b) | 1.38 | 1.73 | 0.59 (4b) | 0.815 |
| Pd  İ2 taban=düzgün | 0.645 (2b) | 2.04 | 2.04 | 0.45 (2b) | 0.585, 0.74, 0.815 |
| P1  İ2 taban=GUE | 0.550 (4b) | 1.24 | 1.36 | 0.89 (4b) | — |
| P0  İ2 taban=Poisson | 0.854 (4b) | -0.36 | -0.63 | — | 0.815 |

---


## Hüküm


### 0. Kontrol tuttu

`taban`, 152'nin A4'ünü yeniden üretti: yirmi Newton iterasyonunun her birinde aynı maks|F| (1.199e+00 → 8.761e-02), aynı σ_ds²=0.1128 / σ_η²=0.0230 / c₁=−0.00918, ve beş bandın Γ'sı dördüncü basamağa kadar aynı. Aşağıdaki farklar konfigürasyondan geliyor.


### 1. İ1a (ε>0, görevde yazıldığı gibi itme) — fazı KÖTÜLEŞTİRİYOR

Doz arttıkça faz oranı monoton biçimde uzaklaşıyor: taban 1.72 → R1 1.80 → R2 1.96 → R3 2.24 (b1–3). Küçük-s göstergesi üç koşuda da tam **0.0000** — itme, zaten aşırı katı olan gazı büsbütün kristalleştiriyor; σ_ds² 0.1128'den 0.0328'e (gerçek 0.1674) çöküyor. Görevin öngördüğü malzeme bu değil: gazda eksik olan itme değil, itmenin fazlası var.


### 2. İ1b (ε<0, kalibrasyon hedefinin istediği yön) — kısmi, sonra doyuyor

Kalibre doz üç farklı süpürme sayısında aynı hedefi tutturuyor (P(s<0.3) = 0.02539 / 0.02324 / 0.02426, gerçek 0.02420) ve **aynı fazı** veriyor: b1–3 oranı 1.546 / 1.547 / 1.536. Süpürme sayısı sonucu değiştirmiyor — sonuç gerçekten dozun kendisine bağlı, sayısal ayrıntıya değil.

Kalibre dozda kazanç mütevazı: 1.72 → 1.55 (b1–3). Dozu 2.7 / 5.4 / 8.1 kat büyütünce oran 1.42 / 1.43 / 1.46'de **doyuyor** — daha fazla yumuşatma daha fazla kazanç getirmiyor, N5z'de geri bile dönüyor. 1.00'e gitmiyor. (Uyarı: N5y/N5z'de kelepçe hamlelerin %29–39'unda bağlıyor — bölüm 2 — yani doygunluğun bir kısmı fizik değil, kelepçe olabilir. Kelepçesiz daha büyük doz sıralamayı bozardı ve gaz olmaktan çıkardı; bu yüzden bu yönde daha ileri gidilmedi. Güvenilir en iyi İ1b koşusu N5x'tir: kelepçe %13.5.)

Bedeli ağır: aşırı dozda σ_η² 0.0230 → 0.1464 (gerçek 0.0227), c₁ -0.00918 → -0.07426 (gerçek −0.01158), P(s<0.3) 0.21955 (gerçeğin 9 katı). 152'nin bulduğu gerilim aynen sürüyor: fazı iyileştiren her hamle S3'ü bozuyor.


### 3. J kontrolü — kazanç σ_ds² artefaktı DEĞİL

İ1b fazı düşürürken σ_ds²'yi de büyütüyor. Faz hangisini görüyor: kısa-menzil YAPISINI mı, yalnız σ_ds²'yi mi? J koşuları aynı σ_ds²'yi yapısız (bağımsız, inkoherent) titreşimle üretir. Eşleşmeler:

| çift | σ_ds² | faz oranı b1–3 |
|---|---|---|
| N5  İ1 ε=−0.0037 ×5 | 0.1525 | **1.547** |
| J14 titreşim σ_j=0.141 | 0.1527 | 1.709 |
| N5x İ1 ε=−0.010 ×5 | 0.2485 | **1.417** |
| J26 titreşim σ_j=0.260 | 0.2337 | 1.728 |

Aynı σ_ds²'de yapılı yumuşatma yapısız titreşimden belirgin biçimde daha iyi (1.55 vs 1.71; 1.42 vs 1.73). Yani kazanç varyans artefaktı değil — kısa-menzil yapısı fazı gerçekten kıpırdatıyor. 152'nin 'titreşim fazı hiç kıpırdatmıyor' bulgusu da doğrulanıyor: J14/J26 tabana göre b1–3'te hiç iyileşme vermiyor.


### 4. İ2 (itmeli taban + boyalı merdiven) — en iyi faz, ama kurgudan gelen yapaylıkla

Üç koşu aynı boyayı paylaşır, yalnız taban istatistiği değişir:

| taban | var(s) girdide | faz oranı b1–3 | ölçülen σ_ds² | ölçülen P(s<0.3) |
|---|---|---|---|---|
| s≡1 taban | 0.000 | 2.036 | 0.1627 | 0.05473 |
| GUE-surmise gap tabanı | 0.178 | 1.364 | 0.3367 | 0.10752 |
| üstel (Poisson) gap tabanı | 1.004 | -0.628 | 1.1216 | 0.30972 |

Düzgün taban (Pd) tabandan **daha kötü** (2.04 vs 1.72) — yani tek-geçiş boyanın kendisi bir kazanç getirmiyor; öz-tutarlılığı atmak zarar veriyor. Poisson tabanı (P0) sinyali tümüyle yok ediyor (faz oranları negatif, R uydurmaları ölçülebilen dört bantta 1.76–2.33 rad artıkla düşüyor). Kazancı veren yalnız GUE tabanıdır: P1, tüm taramanın en düşük oranını veriyor (1.36 b1–3, 1.24 b1–4).

**Ama P1 temiz bir ölçüm değil.** Tek geçiş boya |Δz| en fazla 0.725 üretiyor (ḡ=0.522) — bir ortalama aralıktan büyük yer değiştirme. Sonuç: komşu çiftlerin %1.71'i sıralamayı bozuyor, zincirin gerektirdiği sıralama bunları sıfıra yakın aralıklara çeviriyor. Bunun kanıtı Pd'dir: tabanında s≡1 olduğu için P(s<0.3)=0.00000 girmesi gerekirken ölçümde 0.05473 çıkıyor — bu sayının **tamamı** çaprazlama yapaylığıdır. Aynı sebeple P1'in P(s<0.3)'ü tabanındaki 0.02709'dan 0.10752'e şişiyor (gerçeğin 4.4 katı) ve σ_ds² 0.3367'e (gerçeğin 2 katı) çıkıyor. P1'in düşük faz oranı bu yüzden GUE itmesine değil, kısmen o yapaylığa da yazılabilir; ayrıca |Γ| genlikleri gerçeğin çok altında (0.40–0.46 vs gerçek 0.63–0.81), s_Γ skoru 0.550 ile N5x'in 0.217'inden belirgin kötü.

Ek fizik notu: gerçek σ_ds²=0.1674 ile GUE surmise'ın 0.1781'i zaten aynı yerdedir. Yani gerçek sıfırların aralık varyansını GUE tek başına açıklıyor; üstüne BAĞIMSIZ bir merdiven yer değiştirmesi eklenecek yer yok — eklenince P1'de olduğu gibi 0.34'e fırlıyor. Gerçekte merdiven ile GUE dalgalanması bağımsız iki katman değil, aynı şeyin iki yazımıdır; İ2'nin kurgusu (bağımsız taban + üstüne boya) bunu ihlal ediyor.


### 5. Hangi yol daha iyi

Faz oranına tek başına bakılırsa İ2/GUE (P1, b1–3 1.36), İ1b'nin en iyisinden (N5x 1.42) bir tık önde. Ama İ1b her başka ölçütte kazanıyor: s_Γ 0.217 vs 0.550; |Γ| genlikleri gerçeğin aralığında; beş bandın BEŞİ DE ölçülebilir ve R uydurmalarının hepsi tutuyor (R-artık 0.002–0.015) — 152'nin sekiz ve 153'ün on beş koşusu içinde bunu başaran tek koşu N5x'tir; ve kurgusunda çaprazlama yapaylığı yok (sıra bozan çift 0, pencere yoğunluğu 4e−05·ḡ içinde korunmuş). **İ1b daha iyi yoldur**; İ2'nin sayısı daha parlak görünse de güvenilirliği düşüktür.


### 6. Faz ne kadar indi — dürüst muhasebe

Görevin sorusu: 1.4'ten 1.0'a çekiyor mu — kısmen mi, hiç mi?

| ölçüt | taban | en iyi İ1b (N5x) | en iyi İ2 (P1) | kapanan pay |
|---|---|---|---|---|
| ort b1–4 | 1.72 | 1.32 | 1.24 | %55 (N5x) |
| ort b1–3 | 1.72 | 1.42 | 1.36 | %42 (N5x) |
| b3 tek başına | 1.32 | 1.24 | 1.25 | %26 (N5x) |

**Kısmen — ve payın büyük kısmı zayıf bantlardan geliyor.** Bant bazında bakınca: τ=0.5375 2.32 → 1.76 (ama bu bantta yalnız 33 aday var), τ=0.585 1.50 → 1.25 (gerçek bir hareket, 144 aday), τ=0.66 1.32 → 1.24 — **en sağlam bantta %32'lik fazlalığın yalnız dörtte biri kapanıyor.** τ=0.74 ve 0.815 tabanda zaten patlamış/sarmıştı; oradaki 'iyileşme' esas olarak bantların ölçülebilir hale gelmesidir, ki bu da bir kazanç ama faz kapanışı değil.


### 7. Ne kaldı

- Eksik malzeme, GUE seviye itmesinin sentetik gaza EKLENMESİ değildir; ölçüm bunun tersini söylüyor (sentetik gaz zaten aşırı katı). Kısa-menzil yumuşatması doğru yön ama kapanış b3'te ~%25, b1–3'te ~%42 kalıp doyuyor.
- Yumuşatmanın yapısı yanlış: kalibre dozda P(s<0.3) gerçeği tutturuyor (0.02324 vs 0.02420) ama P(s<0.1) 0.00950 ile gerçeğin (0.00096) 10 katı — gevşetme GUE'nin s² yükselişini değil, bir kümelenme kuyruğu üretiyor. Doğru kısa-menzil ŞEKLİNİ kuran bir gevşetme (ör. Dyson log-gazı, 1/s kuvveti, ve ÖZ-TUTARLI merdivenle birlikte) sıradaki denemedir.
- İ2'yi kurtarmanın yolu, boyayı tek geçişte değil öz-tutarlı çözmektir; öyle olunca taban gap'leri ile merdiven artık bağımsız olmaz ve çaprazlama yapaylığı da kalkar. Bu, 153'ün kurgusunda kasten yoktu (görev tek geçiş istedi) ama sonuç onu işaret ediyor.
- 152'nin 'faz KATI' hükmü tam olarak yıkılmadı, yumuşatıldı: faz kıpırdıyor, ama en güvenilir bantta %6, ve 1.00'e giden yol görünmüyor.


---


## Dürüstlük notları

- **Sapma ve sebebi:** ε<0 (İ1b) ve J kontrolü görevde istenmedi. İ1b, görevin kendi kalibrasyon hedefinin (küçük-s kuyruğunu gerçeğe yaklaştırmak) ε>0 ile sağlanamaz olduğu ÖLÇÜLDÜĞÜ için eklendi (bölüm 1). J, İ1b'nin kazancının bir σ_ds² artefaktı olup olmadığını ayırmak için eklendi (bölüm 3). İkisi de tabloda ayrı grup.
- **Patlayan bantlar ve sebepleri.** Γ_rot'un paydası (on₀−off₀), η'nın asal çizgilerindeki gücün çizgi-dışı güçten farkıdır. Merdiven yapısı zayıfladığında ya da düzensizlik onu bastırdığında bu fark sıfıra gider ve Γ patlar. Bu taramada patlayanlar: `taban` τ=0.74 (|Γ|=2.29; 152'de de aynı bant A4'te patlamıştı), `R1` τ=0.815 (|Γ|=6.73 — en hafif itme dozunda bile; bu bandın niçin bu kadar kırılgan olduğu ayrıca sınanmadı), `J26` τ=0.815 (|Γ|=3.20), `Pd` τ=0.585/0.74/0.815, `P0` τ=0.815 (|Γ|=4.26). Hiçbiri ortalamalara katılmadı.
- **İ2'nin çaprazlama yapaylığı** (bölüm 4) raporun en zayıf halkasıdır ve P1'in 'en iyi faz' sonucunu doğrudan gölgeler; gizlenmedi, hükümde de ağırlığı düşürüldü.
- **|Γ| alt sınırı yok.** 152'nin kuralı yalnız üstten eler (|Γ|>1.5). J14'ün τ=0.74 bandı |Γ|=0.06 ile çökmüş olmasına rağmen kurala göre sayılıyor ve oranı −0.51 çıkıp J14'ün b1–4 ortalamasını yapay biçimde 1.15'e indiriyor. Bu yüzden hükümde b1–4 değil b1–3 kullanıldı; kural 152 ile karşılaştırılabilirlik için değiştirilmedi.
- **Süreler ölçüt değil.** Koşular sırasında aynı makinede başka bir iş (`154_sentetik.py`) 6 süreçle koştu; 8 çekirdeğe 12 süreç düştü. `taban` bu yüzden 20.4 dk sürdü (152'de aynı hesap 10.2 dk). Sayısal sonuç etkilenmez (aynı maks|F| dizisi), yalnız süre sütunu anlamsızdır — o yüzden tabloya konmadı.
- **Uydurma yok:** yukarıdaki tüm tablo sayıları `scratchpad/153/ozet_*.json`'dan otomatik yazıldı. Elle girilen tek sayılar gerçek referans bantları, gerçek S3 üçlüsü, ölçülmüş gerçek P(s) hedefleri ve GUE/Poisson teorik değerleridir.


Ham çıktılar: `153_configs/153_gaz.py <konfig>`; her koşunun tam log'u ve `ozet_*.json`'u üretilmiştir.

