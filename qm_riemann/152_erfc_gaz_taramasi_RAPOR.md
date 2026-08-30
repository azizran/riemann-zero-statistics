# 152 — erfc-kesimli sentetik gaz taraması

151b'nin (`151b_sentetik_gaz2.py`) ölçüm zinciri birebir korunarak **8 koşu** yapıldı: görevin istediği 5 konfigürasyon (A1–A3, B1–B2), bir kontrol (`keskin`) ve tarama sırasında ortaya çıkan iki ek koşu (A4, C1 — gerekçeleri aşağıda). Amaç: keskin-kesimli sentetik gazın gerçeğe göre ~1.45× dik anormal fazını hangi değişikliğin kapattığını bulmak.

> **Kısa hüküm:** hiçbiri kapatmıyor. H-A (erfc) *başka* bir açığı — σ_η² ve c₁'i — neredeyse tam kapatıyor; H-B (titreşim) yalnız σ_ds²'yi kıpırdatıyor; **anormal fazın dikliği sekiz koşunun hepsinde değişmeden duruyor.** Ayrıntı en altta.

**Koşulan kod:** `152_configs/152_gaz.py` (tek kod yolu, konfigürasyon `KONFIG` sözlüğünde tek satır — sekiz ayrı kopyada zincirin kazara ayrışma riskini almamak için).

**Farklar (151b'ye göre, kasıtlı ve tek):** Newton 20 iterasyon (151b: 40 — görev notu 20'yi yeterli sayıyor); S(z) değerlendirmesi 6 süreçli havuzla paralel (seri sonuca göre maks fark ~1e-15, yalnız toplama sırası). `keskin` kontrolü tam da bu iki farkın bantları kaydırmadığını sınamak için var.

**erfc:** stdlib `math.erfc` (tam; 1−erf'in iptalinden kaçınır). τ_q = ω_q/L_hedef, L_hedef = 12.030, merdiven 15450 çizgi (τ≤1.00).


## Bantlar — Γ_rot(Re, Im), |Γ|, R

Hücre: `(Re,Im) |Γ| R`. `*` = R uydurması fazı tutturamadı (faz-artığı > 0.05 rad) — o R anlamsız, skora katılmadı. `‡` = |Γ| > 1.5: ayrışım paydası sıfıra gitmiş, bant SAYISAL OLARAK PATLAMIŞTIR; ölçüm değildir, hiçbir ortalamaya girmez. `—` = bant ölçülemedi.


| τ̄ | gerçek | 151b(40it) | keskin (kontrol) | A3 erfc 0.90/0.05 | A1 erfc 0.75/0.10 | A4 erfc 0.68/0.125 | A2 erfc 0.60/0.15 | B1 titreşim 0.05 | B2 titreşim 0.10 | C1 = A4 + titreşim 0.10 |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.5375 | (+0.787,+0.206) 0.81 0.40 | (+0.691,+0.414) 0.81 1.58 | (+0.693,+0.415) 0.81 1.58 | (+0.681,+0.410) 0.79 1.48 | (+0.649,+0.408) 0.77 1.23 | (+0.600,+0.406) 0.72 1.23 | (+0.568,+0.313) 0.65 1.00 | (+0.694,+0.412) 0.81 1.48 | (+0.696,+0.407) 0.81 1.30 | (+0.596,+0.397) 0.72 1.07 |
| 0.5850 | (+0.550,+0.488) 0.74 1.07 | (+0.401,+0.685) 0.79 2.70 | (+0.401,+0.685) 0.79 2.73 | (+0.399,+0.679) 0.79 2.38 | (+0.343,+0.663) 0.75 2.08 | (+0.322,+0.619) 0.70 2.02 | (+0.433,+0.418) 0.60 1.35 | (+0.404,+0.686) 0.80 2.58 | (+0.409,+0.683) 0.80 2.27 | (+0.320,+0.616) 0.69 1.75 |
| 0.6600 | (+0.118,+0.614) 0.63 1.85 | (-0.195,+0.765) 0.79 3.00 | (-0.196,+0.765) 0.79 (3.00)* | (-0.215,+0.763) 0.79 (3.00)* | (-0.242,+0.721) 0.76 (3.00)* | (-0.183,+0.714) 0.74 (3.00)* | (+4.74,-0.94) **‡** | (-0.192,+0.766) 0.79 (3.00)* | (-0.180,+0.765) 0.79 3.00 | (-0.199,+0.724) 0.75 2.75 |
| 0.7400 | (-0.148,+0.498) 0.52 2.12 | (-0.732,+0.405) 0.84 3.00 | (-0.731,+0.404) 0.84 (3.00)* | (-0.744,+0.375) 0.83 (3.00)* | (-0.928,+0.555) 1.08 (3.00)* | (-1.84,-1.37) **‡** | (-1.60,-0.26) **‡** | (-0.731,+0.403) 0.83 (3.00)* | (-0.723,+0.403) 0.83 (3.00)* | (-0.219,-0.267) 0.35 (0.00)* |
| 0.8150 | (-0.103,+0.655) 0.66 1.50 | (-1.001,-0.195) 1.02 — | (-0.996,-0.199) 1.02 (0.00)* | (-0.962,-0.249) 0.99 (0.00)* | (-0.122,-0.042) 0.13 (0.00)* | (-1.259,-0.557) 1.38 (0.00)* | (+0.291,+0.136) 0.32 0.33 | (-0.968,-0.169) 0.98 (0.00)* | (-0.922,-0.133) 0.93 (0.00)* | (+0.019,+0.366) 0.37 1.20 |

## S3 — varyanslar


| konfig | değişiklik | σ_ds² | σ_η² | c₁ | süre |
|---|---|---|---|---|---|
| **gerçek (zeros6)** | — | **0.1674** | **0.0227** | **−0.01158** | — |
| 151b (40 it, rapor) | keskin | 0.1284 | 0.0769 | −0.03545 | — |
| keskin (kontrol) | değişiklik yok — 151b, 20 iterasyon | 0.1280 | 0.0767 | -0.03535 | 11.7 dk |
| A3 erfc 0.90/0.05 | a_q → a_q·½erfc((τ_q−0.90)/0.05) | 0.1208 | 0.0636 | -0.02904 | 10.5 dk |
| A1 erfc 0.75/0.10 | a_q → a_q·½erfc((τ_q−0.75)/0.10) | 0.1152 | 0.0380 | -0.01771 | 10.7 dk |
| A4 erfc 0.68/0.125 | a_q → a_q·½erfc((τ_q−0.68)/0.125) — S3'ü oturtan kesim | 0.1128 | 0.0230 | -0.00918 | 10.2 dk |
| A2 erfc 0.60/0.15 | a_q → a_q·½erfc((τ_q−0.60)/0.15) | 0.1085 | 0.0117 | -0.00232 | 10.6 dk |
| B1 titreşim 0.05 | Newton sonrası z += N(0, 0.05·ḡ_t), rng(7) | 0.1329 | 0.0820 | -0.03795 | 0.3 dk |
| B2 titreşim 0.10 | Newton sonrası z += N(0, 0.10·ḡ_t), rng(7) | 0.1478 | 0.0979 | -0.04580 | 0.3 dk |
| C1 = A4 + titreşim 0.10 | erfc(0.68/0.125) VE titreşim 0.10 — birlikte | 0.1331 | 0.0449 | -0.01979 | 10.3 dk |

## Anormal fazın dikliği — arg Γ_rot(sentetik) / arg Γ_rot(gerçek)

Hedef 1.000. Kapanış ölçütü budur: 151b'nin sorunu fazın ~1.45× dik olmasıydı.


| konfig | τ=0.5375 | τ=0.585 | τ=0.66 | τ=0.74 | τ=0.815 | **ort (b1–4)** |
|---|---|---|---|---|---|---|
| **gerçek** | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | **1.00** |
| 151b (40 it) | 2.11 | 1.43 | 1.32 | 1.42 | 1.93 | **1.57** |
| keskin (kontrol) | 2.11 | 1.44 | 1.32 | 1.42 | 1.93 | **1.57** |
| A3 erfc 0.90/0.05 | 2.12 | 1.43 | 1.34 | 1.44 | 1.97 | **1.58** |
| A1 erfc 0.75/0.10 | 2.19 | 1.51 | 1.37 | 1.40 | 2.01 | **1.62** |
| A4 erfc 0.68/0.125 | 2.32 | 1.50 | 1.32 | ‡ | 2.06 | **1.72** |
| A2 erfc 0.60/0.15 | 1.97 | 1.06 | ‡ | ‡ | 0.25 | **1.51** |
| B1 titreşim 0.05 | 2.09 | 1.43 | 1.32 | 1.42 | 1.92 | **1.56** |
| B2 titreşim 0.10 | 2.07 | 1.42 | 1.30 | 1.42 | 1.90 | **1.55** |
| C1 = A4 + titreşim 0.10 | 2.30 | 1.50 | 1.33 | 2.16 | 0.88 | **1.82** |

## Skor

`s_Γ` = bant BAŞINA ortalama (|ΔRe| + |ΔIm|), bant 1–4 içinden patlak olmayanlar üzerinden. **Toplam değil ortalama** — patlak bantları eleyince terim sayısı konfigürasyondan konfigürasyona değişiyor; toplam kullanılsaydı en çok bandı patlayan konfigürasyon (A2) sahte biçimde 'en iyi' görünürdü. Kaç bant üzerinden alındığı parantezde; 2 bantlık bir ortalama 4 bantlıkla aynı ağırlıkta okunmamalı. `s_R` aynı mantıkla R için. Küçük = gerçeğe yakın.


| konfig | s_Γ (bant başına) | faz oranı (b1–4) | s_R | patlak ‡ |
|---|---|---|---|---|
| keskin (kontrol) | 0.448 (4b) | 1.57 | 1.42 (2b) | — |
| A3 erfc 0.90/0.05 | 0.463 (4b) | 1.58 | 1.19 (2b) | — |
| A1 erfc 0.75/0.10 | 0.507 (4b) | 1.62 | 0.92 (2b) | — |
| A4 erfc 0.68/0.125 | 0.382 (3b) | 1.72 | 0.89 (2b) | 0.74 |
| A2 erfc 0.60/0.15 | 0.256 (2b) | 1.51 | 0.69 (3b) | 0.66, 0.74 |
| B1 titreşim 0.05 | 0.446 (4b) | 1.56 | 1.29 (2b) | — |
| B2 titreşim 0.10 | 0.437 (4b) | 1.55 | 1.08 (3b) | — |
| C1 = A4 + titreşim 0.10 | 0.501 (4b) | 1.82 | 0.64 (4b) | — |

---

## Hüküm

### 0. Kontrol tuttu — kıyas meşru

`keskin` (20 iterasyon, paralel çözücü) 151b'nin (40 iterasyon, seri) rapor edilen sayılarını 3–4 anlamlı basamak yeniden üretti: bantlar (+0.693,+0.415)/(+0.401,+0.685)/(−0.196,+0.765)/(−0.731,+0.404) vs 151b'nin (+0.691,+0.414)/(+0.401,+0.685)/(−0.195,+0.765)/(−0.732,+0.405); σ_ds² 0.1280 vs 0.1284, σ_η² 0.0767 vs 0.0769, c₁ −0.03535 vs −0.03545. Faz oranı ortalaması iki koşuda da 1.57. Yani aşağıdaki farklar konfigürasyondan geliyor, iterasyon sayısından ya da paralelleştirmeden değil.

### 1. H-B (inkoherent seyreltme) — fazı HİÇ kıpırdatmıyor

Titreşim bantları neredeyse aynen bırakıyor: σ_j=0.10'da bile τ=0.5375 bandı (+0.693,+0.415) → (+0.696,+0.407). Faz oranı 1.570 → 1.552 (%1). Buna karşılık σ_ds²'yi gerçeğe DOĞRU itiyor (0.1280 → 0.1478, gerçek 0.1674) ama σ_η²'yi (0.0767 → 0.0979) ve c₁'i (−0.0354 → −0.0458) gerçekten UZAKLAŞTIRIYOR. Konum seyreltmesi dispersiyonun eksik malzemesi değil.

### 2. H-A (erfc-kesim) — S3'ü kapatıyor, fazı kapatmıyor

σ_η² ve c₁ τ_c ile monoton: keskin 0.0767/−0.0354 → A3(0.90) 0.0636/−0.0290 → A1(0.75) 0.0380/−0.0177 → A2(0.60) 0.0117/−0.0023. Gerçek (0.0227/−0.01158) A1 ile A2 ARASINDA kalıyor; doğrusal ara değer her iki nicelikte de τ_c≈0.66–0.69 veriyor. **A4 (τ_c=0.68, Δ=0.125) tam o noktadır ve gerçeği vuruyor: σ_η²=0.0230 (gerçek 0.0227), c₁=−0.00918 (gerçek −0.01158).** Sentetik gazın şişkin η-varyansı ve fazla negatif c₁'i, merdivenin yüksek-τ keskin kesiminin artığıymış — erfc yumuşatması onu birebir kaldırıyor.

Ama faz kapanmıyor: A4'te oran 1.72 (keskin 1.57'den daha da dik), τ=0.66 bandı 1.32 → 1.32. Üstelik σ_ds² her erfc kesiminde gerçeğin TERSİNE gidiyor (0.1280 → 0.1128; gerçek 0.1674).

### 3. A2'nin 0.98'lik faz oranı sahte — okumayın

A2 ham ortalamada 0.98 (mükemmele yakın) gösteriyordu. Sebebi fizik değil: en ağır kesim ayrışım paydasını (on₀−off₀) sıfıra sürüyor, τ=0.66 bandı |Γ|=4.83'e patlıyor ve **negatif** faz oranı (−0.14) üretip ortalamayı yapay olarak 1'e çekiyor. |Γ|>1.5 elemesinden sonra A2'den geriye 4 bandın 2'si kalıyor (oranlar 1.97 ve 1.06) — kapanış değil, sinyal çöküşü. Aynı sebeple A2'nin toplam-skoru da aldatıcıydı; rapordaki skorlar bu yüzden bant BAŞINA ortalamadır.

### 4. İkisi birlikte de kapatmıyor (C1)

"Yoksa ikisi de mi kısmi?" sorusunu doğrudan sınamak için A4'ün kesimi + en güçlü titreşim birlikte koşuldu. Sinerji yok, çelişki var: titreşim σ_η²'yi A4'ün tam oturttuğu yerden geri kaldırıyor (0.0230 → 0.0449, gerçek 0.0227'nin üstüne), σ_ds² 0.1331'de yine eksik kalıyor, faz oranı 1.82 ile en kötülerden biri. İki knob σ_η² üzerinde birbirini iptal ediyor.

### 5. Asıl bulgu: faz KATI

Bant bazında faz oranı, iyi ölçülen orta bantlarda sekiz koşu boyunca pratik olarak sabit (τ=0.66'da 7 koşuda 1.30–1.37; τ=0.74'te 5 koşuda 1.40–1.44 — %5'ten dar). Kenar bantlar (0.5375 ve 0.815) her koşuda gürültülü; 0.815 zaten sekiz koşunun altısında fazını sarıyor, yani ondan hüküm çıkarılamaz.


| bant | aralık (sağlam bantlar) | kaç koşu | elenen |
|---|---|---|---|
| τ=0.5375 | **1.97 – 2.32** | 8 | — |
| τ=0.585 | **1.06 – 1.51** | 8 | — |
| τ=0.66 | **1.30 – 1.37** | 7 | A2(‡) |
| τ=0.74 | **1.40 – 1.44** | 5 | A4(‡), A2(‡), C1(faz sarmış) |
| τ=0.815 | **0.25 – 0.88** | 2 | keskin(faz sarmış), A3(faz sarmış), A1(faz sarmış), A4(faz sarmış), B1(faz sarmış), B2(faz sarmış) |

Merdiven genliklerine uygulanan hiçbir pencere ve konumlara eklenen hiçbir inkoherent titreşim bu sayıları kıpırdatmıyor — ve hiçbir bant, hiçbir konfigürasyonda 1.00'e yaklaşmıyor (A2'nin 1.06'sı hariç, o da bandın komşularını patlatan koşudan geliyor). Dispersiyon fazı, S3 istatistiklerinden BAĞIMSIZ bir kanaldan geliyor: A4 σ_η²'yi 3.3 kat düşürüp gerçeğe oturturken faz oranını 1.57'den 1.72'ye çıkardı — iki gözlenebilir ters yönde hareket etti, yani ortak bir sebeple açıklanamazlar.

### 6. Ne kaldı

- Eksik malzeme, merdiven genlik tayfında ya da konum koheransında DEĞİL. İkisi de tarandı, ikisi de fazı bırakmadı.
- Sentetik gaz aynı anda hem az dağınık (σ_ds² 0.11–0.15 vs 0.167) hem fazla dik (faz ×1.3–1.4). Tek bir 'yumuşatma' knobu bu ikisini birden düzeltemez; ters yönlere çekiyorlar.
- Sıradaki aday, saf asal-merdiven gazında BULUNMAYAN çift korelasyonu: gerçek sıfırların GUE seviye itmesi. Bu, konumları bağımsızca titretmekle (H-B) taklit edilemez — H-B'nin fazı hiç kıpırdatmaması da tam bunu söylüyor. Ölçüm zincirinin kendi geometrisi (0.52-taban regresyon, bant kenarları) ikinci aday.


---

Ham çıktılar: `152_configs/152_gaz.py <konfig>`; her koşunun tam log'u ve `ozet_*.json`'u üretilmiştir. Bu rapordaki tüm sayılar o JSON'lardan otomatik yazılmıştır (gerçek referans ve 151b(40it) referansı dışında elle sayı girilmemiştir).
