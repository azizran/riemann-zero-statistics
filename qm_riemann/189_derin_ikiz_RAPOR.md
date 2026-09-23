# 189 — DERİN İKİZ: ZARFIN KİMLİĞİ SINAVI

**Soru (KALEM_DERIN_IKIZ_23EYL2026):** 188, iptal çekirdeğinin kinematik-genel
olduğunu (corr 0.999) ve derinlik-eşli ζ'nin iki denizde aynı olduğunu
(τ'≤1.00: gerçek 0.0758, ikiz 0.0772) ölçtü. İkizle gerçeğin ζ farkı
(0.077'ye karşı 0.33), ikizin merdiveninin kurgu gereği τ=1'de bitmesinden
geliyor. 4. demirin manşeti "gerçeğin ikiz-fazlasının ~%80'i nazik bir zarf,
w(τ) = 1 − 0.149·τ^1.30" idi. **Bu zarf, ikizin merdiven kesiminin parmak izi
mi?** Sınav: ikizi daha derin merdivenle (τ≤1.10, τ≤1.20) AYNI sadakatli
çözücüyle yeniden kur, 184→187 zincirinden geçir, zarf ön-mühürlü hızla
kapanıyor mu, bak.

Her sayının yanında onu üreten betik adı vardır. Tek dalga koşuldu. Ölümler
kurtarılmadı. Git'e dokunulmadı. Sonuç ORTAK TEFTİŞE sunulur, commit sonra.

---

## Ölçüm tanımı (K0'da donmuş; 164 çözücüsü ve 184-187 makinesi AYNEN)

- **İkizler:** Hderin110 ve Hderin120, 164'ün sadakatli çözücüsüyle kuruldu
  (ızgara braketi + sıralı ilk-kök + korumalı Newton). Seviye c = −½,
  nz = 300 000, zeros6 son-300k penceresi. TEK fark merdiven kesimi:
  lim = exp(D·L_hedef), D = 1.10 ve 1.20 (Hkeskin'de D = 1.00).
  `164_insa.py` düzenlenmedi. `189b_insa.py` modülü importlib ile yükler,
  `merdiven` modül özniteliğini `partial(merdiven, tau_ust=D)` ile değiştirir,
  `KONFIG`'e iki ad ekler ve `main`'i AYNEN çağırır.
- **Zincir:** 184b2.kos → 185b.oz_hesapla → 186b.merdiven_izdusum (ρ≡1) →
  187c ζ defteri (`189c_zincir.py`). Modüller importlib ile yüklenir. Tek
  çalışma-anı yaması: `b185.kinematik` yeni adları 185b'nin Hkeskin koluyla
  aynı satırlarla okur. Çıktılar `scratchpad/189/`'a gider; 184-188
  önbellekleri yalnız okunur.
- **Nicelikler** (bant = τ∈[lo,hi), HAVUZ = [0.45, 0.86)):
  - w = Σ|c^ölç|/Σae
  - w^öz = Σ|c^öz|/Σae (ρ≡1)
  - m = w − w^öz
  - **r_D = w_g/w_HD**, σ_r = r·√((se_wg/wg)²+(se_wh/wh)²) (184 konvansiyonu)
  - σ_ε = √(Σ(g−ḡ)²/N)/ḡ (185b)
  - ζ = Σ Δ·conj(karışım)/Σ|karışım|² (187c; modül + açı). Genlik-okuma
    Γ = 1 − m_ölç/m^kesik.
  - Jackknife: 8-blok loo, se = √(7/8·Σ(θ_i−θ̄)²). Bloklar dört gazda aynı
    (N = 299 999).
- **Öngörü (TÜRETİLEN, kaptan, 188 haritası):** m^kesik_Hk = m_Hk/(1−z_1.00),
  m_HD = m^kesik_Hk·(1−z_D), r_D = w_g/(w^öz_Hk + m_HD), z_D = −Re Σ_{τ'≤D} K_Hk.
  Bu aritmetik kinematiği (w^öz) Hkeskin'de dondurur.

---

## K0 — ÖN-KAYIT  [189a_onkayit.py]

`189/ONKAYIT_189.json` yazıldı:
- **sha256 = 2ff314298ec49843…** (tam:
  2ff314298ec49843e21a894b2bc91453a66bcac591fec0500cc62e8538c4db7c)
- **damga Wed Sep 23 13:16:17 +0300 2026**
- kalem sha256 = 5a8d6b56…

Sıralama dosya doğum zamanlarından: ön-kayıt 13:16:17, ızgara sınavı logu
13:16:48, ilk inşa 13:23:44. Yani ön-kayıt her ölçümden ÖNCE yazıldı.

Donanlar:
- **Öngörü tablosu kalemden AYNEN.** Kaynağından yeniden hesaplandı
  (185/K1_faktorler + 188/harita_K_Hkeskin). 8 bantta 6 sayısal sütun 4
  haneye kadar EŞİT, kapanan-% sütunları tam sayıda EŞİT (maks fark 0.0).
- HAVUZ için aynı aritmetik (kalemde yok, bilgi amaçlı): r_1.00 = 0.9096,
  z_D = 0.0772 / 0.1794 / 0.2430, **r_1.10^öng = 0.9643, r_1.20^öng = 1.0017**
  (kapanan %61 / %102).
- H-189a..d eşikleri ve hüküm kuralları, ikiz adları, ızgara-sınavı kuralı,
  inşa kapıları, zincir ve makine mührü.
- K3 fit makinesi (184c AYNEN) ve D→∞ modelleri (A birincil, B ikincil).

Veriden önce JSON'a yazılan operasyonel netleştirmeler:
- H-189a'daki 2σ, σ_Δ = √(σ_r(1.10)²+σ_r(1.20)²) üzerinden okunur.
- "r HAVUZ'da azalır" ölümü: r_1.10 < r_1.00 ya da r_1.20 < r_1.10 − 2σ_Δ.
- Ne mühür ne ölüm çıkarsa hüküm KAYIT olur.
- H-189b derinlik başına 8/8 bant ister; iki derinlik de tutarsa MÜHÜR.
- H-189d'de ζ_HD = 187c modülü |ζ| (−Re ζ ikincil).

---

## K1 — IZGARA SINAVI + İNŞA  [189b_insa.py]

**Izgara sınavı.** İlk 5000 seviye (c = −½), aynı derin merdiven, h = 0.015
ve h = 0.0075 karşılaştırıldı. "İlk-kök farkı" |Δz| > 1e-6 olarak tanımlandı.

| ikiz | çizgi | ω_max → en kısa λ | fark / 5000 | maks \|Δz\| | karar |
|---|---|---|---|---|---|
| Hderin110 | 46 168 | 13.233 → 0.4748 | **0** | 9.3e-10 | **h = 0.015** |
| Hderin120 | 139 590 | 14.436 → 0.4352 | **0** | 9.3e-10 | **h = 0.015** |

(Hkeskin: 15 450 çizgi, en kısa λ = 0.522. h = 0.015, 1.20'de de en kısa dalga
boyunun 1/29'u. Farklar |F|≤1.9e-9 çözüm tabanında kalıyor.)

**İnşa kapıları** (nohup + ≤30 s yoklama; iki inşa ASLA aynı anda koşmadı):

| ikiz | ızgara noktası | S-ızgara | Newton it. / ikiye-bölme | **maks\|F\|** | **sıralılık** (min Δz) | **L** | σ_ds² | süre |
|---|---|---|---|---|---|---|---|---|
| Hkeskin (188-K0) | 10 449 281 | 544 s | 3 / 0 | 1.863e-9 | TAM (0.1067) | 12.02959324 | 0.1861 | 10.6 dk |
| **Hderin110** | 10 450 469 | 1582 s | 3 / 0 | **1.863e-9** | **TAM** (0.0994) | **12.02959324** (…241719) | 0.1785 | 30.8 dk |
| **Hderin120** | 10 452 451 | 4958 s | 3 / 0 | **1.863e-9** | **TAM** (0.0955) | **12.02959324** (…241723) | 0.1746 | 95.4 dk |

İki ikizin de üç kapısı GEÇTİ (`189/insa_kapi_*.json`): |F|>1e-8 olan tekne
yok, sıra bozan çift yok. Tanı notları: Σa 27.8 → 44.8 → 73.0; rms S′ 1.894 →
2.087 → 2.280; ΔG<0 kesri 0.192 → 0.225 → 0.253 (F daha sık monoton değil).
Sıralı ilk-kök kuralı tüm tekneleri benzersiz hücrede çözdü (300 000/300 000).
σ_ds² derinlikle gerçeğin 0.1674'üne doğru iniyor.

---

## K2 — ZİNCİR + MAKİNE MÜHRÜ  [189c_zincir.py]

**Makine mührü.** AYNI sarmalayıcı zinciri Hkeskin'e uygulandı (çıktılar
189/ altına yazıldı). Sonuçlar 184/185/186/187'nin kayıtlı Hkeskin
sayılarıyla **BİT-BİT** karşılaştırıldı (`189/muhur_Hkeskin.json`):

| kontrol | sonuç |
|---|---|
| 189/K1_Hkeskin = 184/K1_Hkeskin (tüm alanlar) | TUTTU |
| 189/OZ_Hkeskin = 185/OZ_Hkeskin | TUTTU |
| 189/G1_proj re_bir/im_bir/momentler = 186 | TUTTU |
| 185b defteri 9 satır (r, σ_r, F'ler, w'ler, se_tablo, τ̄) = 185/K1_faktorler | TUTTU |
| 187c ζ defteri Hkeskin (modül/açı/se/genlik/m) = 187/K2_zeta.json | TUTTU |
| 187c ζ defteri gerçek = 187/K2_zeta.json | TUTTU |

Ek mühürler:
- 184c fit makinesi (K3) Hkeskin'de **c = 0.1492, α = 1.304, χ²/dof = 0.343**
  verdi; 184'ün sayıları 0.149 / 1.30 / 0.34 [189d].
- σ_ε(Hk) = 0.43135, σ_ε(gerçek) = 0.40921; ikisi de 185 ile aynı [189c].
- Model A kimlik sınaması: Γ = Γ_HD konunca ölçülen r_HD'yi maks|Δ| = 2.2e-16
  ile geri veriyor [189d].

---

## BANT TABLOLARI — ölçülen vs ön-mühür  [189c → 189d]

**r_D = w_g/w_HD (± σ_r) ve ön-mühür öngörüsü:**

| bant | r_1.00 (Hk) | r_1.10 ölç | öng | r_1.20 ölç | öng | f_1.10 ölç (öng) | f_1.20 ölç (öng) |
|---|---|---|---|---|---|---|---|
| 0.45-0.50 | 0.9445±.0017 | 0.9678±.0017 | 0.9802 | 0.9819±.0015 | 1.0037 | 0.420±.021 (0.64) | 0.674±.016 (1.07) |
| 0.50-0.55 | 0.9351±.0019 | 0.9608±.0016 | 0.9747 | 0.9756±.0015 | 1.0028 | 0.396±.016 (0.61) | 0.625±.012 (1.04) |
| 0.55-0.60 | 0.9253±.0033 | 0.9522±.0031 | 0.9708 | 0.9697±.0032 | 1.0013 | 0.361±.012 (0.61) | 0.594±.012 (1.02) |
| 0.60-0.65 | 0.9169±.0040 | 0.9465±.0033 | 0.9666 | 0.9646±.0034 | 1.0003 | 0.356±.010 (0.60) | 0.574±.006 (1.00) |
| 0.65-0.70 | 0.9096±.0043 | 0.9388±.0036 | 0.9624 | 0.9575±.0041 | 0.9987 | 0.324±.013 (0.58) | 0.530±.012 (0.99) |
| 0.70-0.75 | 0.9021±.0044 | 0.9327±.0039 | 0.9628 | 0.9513±.0039 | 0.9985 | 0.313±.011 (0.62) | 0.502±.015 (0.99) |
| 0.75-0.80 | 0.8978±.0081 | 0.9293±.0074 | 0.9677 | 0.9458±.0077 | 1.0171 | 0.308±.012 (0.68) | 0.470±.008 (1.17) |
| 0.80-0.86 | 0.8931±.0135 | 0.9237±.0130 | 0.9638 | 0.9382±.0129 | 1.0203 | 0.286±.011 (0.66) | 0.422±.013 (1.19) |
| **HAVUZ** | **0.9096±.0051** | **0.9389±.0043** | 0.9643 | **0.9556±.0043** | 1.0017 | **0.324±.010** (0.61) | **0.509±.009** (1.02) |

f_D = (r_D − r_1.00)/(1 − r_1.00). f'nin se'si ortak-blok loo ile hesaplandı.

**Kinematik (w^öz) ve karışım (m) derinlikle** (HAVUZ ± jk; bant uçları):

| nicelik | D = 1.00 | 1.10 | 1.20 | gerçek |
|---|---|---|---|---|
| w_HD (HAVUZ) | 1.3675±.0062 | 1.3248±.0042 | 1.3017±.0041 | w_g = 1.2439 |
| w^öz_HD (HAVUZ) | 0.6672±.0013 | 0.6867±.0019 | 0.6976±.0019 | 0.7144 |
| w^öz_HD (0.45-0.50 … 0.80-0.86) | 0.8135 … 0.5820 | 0.8215 … 0.6130 | 0.8253 … 0.6322 | 0.8314 … 0.6619 |
| m_HD (HAVUZ) | 0.7003±.0053 | 0.6381±.0026 | 0.6040±.0024 | 0.5294 |
| m^kesik_HD (HAVUZ, 187c) | 0.7485 | 0.7598 | 0.7649 | 0.7944 |
| F_KİN = w^öz_g/w^öz_HD (HAVUZ) | 1.0708 | 1.0403 | 1.0241 | 1 |

İkizin öz-terimi w^öz derinlikle gerçeğinkine doğru büyüyor. Bu, aralık
saçılımının küçülmesinin (H-189c) doğrudan izi. Öngörü aritmetiği w^öz'yü
Hkeskin'de sabit tutuyordu.

---

## H-189a — ZARF = MERDİVEN DERİNLİĞİ (birincil)  [189d]

- **HAVUZ r:** 0.9096±0.0051 → 0.9389±0.0043 → 0.9556±0.0043 (monoton artıyor;
  c1 ve c2 ikisi de EVET).
- **Bantlar:** 8/8 bantta r_1.00 < r_1.10 ≤ r_1.20.
- **Kapanan pay:** f_1.10 = **0.324±0.010** (eşik ≥ 0.30), f_1.20 =
  **0.509±0.009**.
- **Hüküm: MÜHÜR.** Ölüm koşullarının hiçbiri ateşlemedi.
- **Dürüst not:** f_1.10, 0.30 eşiğinin yalnız 2.3 jk-se üstünde. Kapanış
  düşük-τ bantlarında daha hızlı (f_1.20 0.67), kesimde daha yavaş (0.42).

## H-189b — SABİT-KİNEMATİK ARİTMETİĞİ (ikincil, nicel)  [189d]

- **1.10'da 0/8 bant, 1.20'de 0/8 bant** max(2σ, 0.01) içinde.
- **Hüküm: ÖLDÜ** (kurtarmasız).
- Ölçülen kapanış öngörünün kabaca **yarı hızında**: HAVUZ f 0.32'ye karşı
  0.61, 0.51'e karşı 1.02. Sapma hep aynı yönde (ölçülen r < öngörü) ve kesime
  doğru büyüyor.

**Ayrışım tablosu (KAYIT):** Δw = w_HD^ölç − w_HD^öng = Δw^öz (kinematik) + Δm
(karışım).

| D | bant | sapma r | tol | Δw | Δw^öz (kin) | Δm (karış.) | kin payı | r^kin |
|---|---|---|---|---|---|---|---|---|
| 1.10 | 0.45-0.50 | −0.0124 | 0.0100 | +0.0153 | +0.0080 | +0.0074 | 0.52 | 0.9737 |
| 1.10 | 0.65-0.70 | −0.0236 | 0.0100 | +0.0327 | +0.0171 | +0.0156 | 0.52 | 0.9499 |
| 1.10 | 0.80-0.86 | −0.0401 | 0.0260 | +0.0566 | +0.0310 | +0.0256 | 0.55 | 0.9414 |
| 1.10 | **HAVUZ** | −0.0254 | — | +0.0349 | +0.0195 | +0.0154 | **0.56** | 0.9499 |
| 1.20 | 0.45-0.50 | −0.0218 | 0.0100 | +0.0259 | +0.0118 | +0.0141 | 0.46 | 0.9936 |
| 1.20 | 0.65-0.70 | −0.0412 | 0.0100 | +0.0542 | +0.0260 | +0.0281 | 0.48 | 0.9785 |
| 1.20 | 0.80-0.86 | −0.0821 | 0.0258 | +0.1077 | +0.0502 | +0.0575 | 0.47 | 0.9803 |
| 1.20 | **HAVUZ** | −0.0461 | — | +0.0599 | +0.0304 | +0.0295 | **0.51** | 0.9778 |

(Tüm 8 bant `189/HUKUM_189.json`'da; kinematik payı her iki derinlikte
0.41–0.57.)

**Okuma:** Sapmanın **~yarısı kinematik değişimden** geliyor: w^öz derinlikle
büyüyor, F_KİN telafisi 1.071 → 1.024'e eriyor. Kalem de bunu "derin ikizin
kinematiği de değişecek" diye uyarmıştı. Öbür yarısı karışımdan geliyor ve
bunun büyük kısmı derin ikizin kendi iptalinin haritadan zayıf olması (bkz.
H-189d ve EK (a)). Yalnız kinematiği değiştirilmiş öngörü (r^kin) bile
ölçülenin üstünde kalıyor.

## H-189c — GERİ-BESLEME KİNEMATİĞİ  [189c, 189d]

| | Hkeskin | Hderin110 | Hderin120 | gerçek |
|---|---|---|---|---|
| σ_ε | 0.43135±0.00057 | **0.42248±0.00014** | **0.41788±0.00011** | 0.40921±0.00011 |

- Δ(1.10−1.00) = −0.00887±0.00044; Δ(1.20−1.10) = −0.00460±0.00006. Farklar
  ortak-blok jackknife'la hesaplandı ve ikisi de ≫ 2 se.
- Gerçeğe doğru kapanan pay: %40.1 (1.10), %60.8 (1.20). Gerçek aşılmadı.
- **Hüküm: MÜHÜR.** σ_ε derinlikle monoton ve gerçeğe doğru azalıyor.
  Kalemin τ'≈1 negatif geri-besleme öngörüsü (sin(πτ'(1+ds)) ≈ −π·ds)
  yönüyle doğrulandı.

## H-189d — HARİTA GEÇERLİLİĞİ (kapı-niteliğinde)  [189c, 189d]

| D | HAVUZ \|ζ_HD\| | açı | z_D (ikiz haritası) | δ | eşik %10 |
|---|---|---|---|---|---|
| 1.00 (Hk) | 0.0768±0.0064 | +179.76°±0.34 | 0.0772±0.0065 | %0.5 | (kontrol) |
| 1.10 | **0.1656±0.0017** | −179.86°±0.22 | 0.1794±0.0030 | **%7.7** | geçti |
| 1.20 | **0.2121±0.0015** | +179.95°±0.19 | 0.2430±0.0028 | **%12.7** | **aştı** |

- **Hüküm: ÖLDÜ** (1.20'de %12.7 > %10).
- İkincil −Re ζ aynı oranları veriyor. Açı her bantta 180°±2.5° içinde
  (iptal saf yıkıcı kalıyor).
- **Bant bazında δ** (1.10 / 1.20):
  - 0.45-0.50: 9.1 / 14.4
  - 0.50-0.55: 6.0 / 14.0
  - 0.55-0.60: 10.3 / 15.1
  - 0.60-0.65: 7.9 / 13.2
  - 0.65-0.70: 6.6 / 11.8
  - 0.70-0.75: 7.2 / 9.1
  - 0.75-0.80: 9.7 / 15.2
  - 0.80-0.86: 2.8 / 8.0 (%)
- Ölçülen iptal haritanın her bantta ALTINDA kalıyor. Açık derinlikle
  büyüyor (HAVUZ oranı 0.995 → 0.923 → 0.873).

---

## K3 — ZARFIN KİMLİĞİ (KAYIT)  [189d]

**r_D(τ) = 1 − c·τ^α** (184c makinesi AYNEN; ± fit-perr / ± loo-jk yeniden fit):

| D | c | α | χ²/dof |
|---|---|---|---|
| 1.00 | 0.1492 ± 0.0080 / ± 0.0026 | 1.304 ± 0.092 / ± 0.044 | 0.343 |
| 1.10 | 0.1155 ± 0.0084 / ± 0.0040 | 1.681 ± 0.131 / ± 0.075 | 0.291 |
| 1.20 | 0.0997 ± 0.0105 / ± 0.0038 | 2.225 ± 0.193 / ± 0.083 | 0.140 |

- Form üç derinlikte de iyi uyuyor; sınıra dayanan fit yok. Ek işaret-serbest
  fit aynı sonucu veriyor.
- **Zarf derinlikle yalnız küçülmüyor, KESİME ÇEKİLİYOR:** c 0.149 → 0.100,
  α 1.30 → 2.22. Merdiven derinleşince düşük-τ zarfı hızla kapanıyor; kalan
  bastırma pencere kesimine yığılıyor.
- **4. demirin α = 1.30'u sabit bir üs değil.** τ≤1.00 merdiveninin değeri
  bu; α derinliğe bağlı.

**D → ∞ kestirimi (MODEL; K0'da donmuş tanımlar):**

- **Model A (birincil):** Kinematik 1.20'de dondurulur. Derinlik limiti
  gerçeğin tam-derinlik iptal payıdır: Γ_g = 0.3335, 187c genlik-okuması;
  188 derinlik-eşli ζ eşitliğine dayanır.
  r_∞ = w_g / (w^öz_H120 + m^kesik_H120·(1 − Γ_g)).
- **Model B (ikincil):** Kinematiğin derinlik eğilimi sürer. w^öz ve m^kesik,
  üç derinlikte Γ_HD'ye (0.064 / 0.160 / 0.210) karşı doğrusal fit edilir ve
  Γ_g'de değerlendirilir. Bu uzun bir dış-değerlemedir; jk se model
  belirsizliğini içermez.

| bant | r_1.00 | A: r_∞ | A: P_der | B: r_∞ | B: P_der |
|---|---|---|---|---|---|
| 0.45-0.50 | 0.9445 | 1.0276±.0008 | 1.50 | 1.0185±.0011 | 1.33 |
| 0.60-0.65 | 0.9169 | 1.0261±.0007 | 1.31 | 1.0090±.0010 | 1.11 |
| 0.70-0.75 | 0.9021 | 1.0263±.0008 | 1.27 | 0.9983±.0004 | 0.98 |
| 0.80-0.86 | 0.8931 | 1.0415±.0009 | 1.39 | 0.9896±.0008 | 0.90 |
| **HAVUZ** | 0.9096 | **1.0302±.0004** | **1.335±.008** | **1.0014±.0002** | **1.015±.003** |

- Duyarlılık: A'da Γ_g yerine |ζ_g| = 0.3287 konunca HAVUZ r_∞ = 1.0271,
  P_der = 1.30.
- **Payların tanımı:** DERİNLİK payı P_der = (r_∞ − r_1.00)/(1 − r_1.00);
  KİNEMATİK payı = 1 − P_der.
- **HAVUZ:**
  - A: derinlik %133.5, kinematik **−%33.5**. Kinematik fark ters işaretli:
    sonsuz-derin, 1.20-kinematikli ikiz gerçeği %3 AŞARDI.
  - B: derinlik %101.5, kinematik −%1.5 (~sıfır).
  - B'de kesim bandı (0.80-0.86) %90 derinlik, +%10 kinematik.
- **A'nın kinematik fark ayrışımı** (HAVUZ): w_H∞ − w_g = −0.0365. Bunun
  −0.0168'i öz-terim (w^öz_H120 0.6976 < gerçek 0.7144), −0.0197'si
  pencere-içi karışım (m^kesik_H120 0.7649 < gerçek 0.7944; ×(1−Γ_g)) payı.
  Model B'de bu iki fark, kinematiğin gerçeğe yürümesiyle kapanıyor.

**"4. demirin %80 zarfının ne kadarı derinlik":**
- **Ölçülen (modelsiz):** τ≤1.20'ye inen ikiz zarfın (HAVUZ 1 − r_1.00 =
  0.0904) **%50.9 ± 0.9'unu** kapattı (τ≤1.10'da %32.4 ± 1.0). 1.10→1.20
  artımı (+%18.5) hâlâ sürüyor.
- **Model:** D→∞'da zarfın **tamamı** merdiven derinliği (A %134, B %102).
  Gerçek-kinematik payı A'da ters işaretli (−%33), B'de ~sıfır (−%1.5).
- **Sonuç:** Hiçbir okumada zarfın gerçek bir kinematik gerçek-ikiz
  farkından geldiğine dair pay yok. 4. demirin "%80 zarf" bileşeni ikizin
  τ = 1 merdiven kesiminin parmak izidir. Ölçülen alt sınır zarfın yarısı,
  model kestirimi tamamı.

---

## ÖN-KAYITSIZ EK (veri görüldükten SONRA; hüküm DIŞI, KAYIT)  [189f_ek.py]

**(a) H-189b sapmasının karışım payı iki parçaya ayrıldı:**
Δm = (m^kesik_HD − m^kesik_öng)(1 − z_D) + m^kesik_HD·(z_D − Γ_HD). HAVUZ:

- **1.10:** Δm = +0.0154 = +0.0007 (pencere-içi karışım) + **+0.0146 (iptal
  açığı)**
- **1.20:** Δm = +0.0295 = +0.0045 + **+0.0250**
- **Uyarı:** D = 1.00'da da iki terim ±0.0096 (lehçe tabanı: Γ_Hk = 0.0644,
  z_1.00 = 0.0772); iki terim burada tam iptal eder. Bu taban derin
  satırlardaki iptal-açığı teriminin bir kısmını açıklar.
- **Kaba bütçe (HAVUZ, 1.20):** sapmanın ~%51'i öz-terim kinematiği, ~%7'si
  pencere-içi karışım, ~%42'si iptal açığı.

**(b) Derinlik-eşli kıyas.** Derin ikizin kendi ζ'sı, GERÇEK haritanın
derinlik-eşli toplamının da altında:

| D | \|ζ_HD\| | Σ K gerçek | oran |
|---|---|---|---|
| 1.10 | 0.1656 | 0.1717±0.0028 | 0.965 |
| 1.20 | 0.2121 | 0.2261±0.0026 | 0.938 |

D = 1.00'da oran 1.013 (188'in eşitliği). Okuma (sınanmadı): τ'>1 çekirdeği,
onu taşıyan kinematik daraldıkça (σ_ε ↓) zayıflıyor. 188'in "kinematik-genel"
çekirdeği ŞEKİLDE genel (corr 0.999), GENLİKTE kinematiğe duyarlı. Bu, H-189d
ölümünün ve H-189b karışım sapmasının ortak kaynağı olabilir; yeni ön-kayıt
ister.

---

## HÜKÜM (eşikler K0'da donmuş; kurtarma yok)

| hipotez/kapı | hüküm | dayanak |
|---|---|---|
| K0 ön-kayıt | **TUTTU** | sha 2ff31429…, 13:16:17; öngörü yeniden hesabı kalemle 4 hane eşit [189a] |
| Izgara sınavı | **GEÇTİ** | 0/5000 (1.10 ve 1.20), maks\|Δz\| 9.3e-10 → h = 0.015 [189b] |
| İnşa kapıları | **GEÇTİ** | maks\|F\| 1.863e-9 / sıralılık TAM / L = 12.02959324, iki ikizde [189b] |
| Makine mührü | **TUTTU** | 6/6 bit-bit; fit 0.1492/1.304/0.343 [189c, 189d] |
| **H-189a** zarf = derinlik | **MÜHÜR** | HAVUZ 0.9096 → 0.9389 → 0.9556; 8/8 bant; f_1.10 = 0.324±0.010 ≥ 0.30 [189d] |
| **H-189b** sabit-kinematik aritmetiği | **ÖLDÜ** | 0/8 (1.10), 0/8 (1.20); kapanış öngörünün ~yarı hızı; sapmanın ~yarısı kinematik (w^öz ↑) [189d] |
| **H-189c** geri-besleme kinematiği | **MÜHÜR** | σ_ε 0.4313 → 0.4225 → 0.4179 (gerçek 0.4092), farklar ≫ 2se [189c, 189d] |
| **H-189d** harita geçerliliği | **ÖLDÜ** | HAVUZ δ = %7.7 (1.10) ama %12.7 (1.20) > %10; açı 180°±0.2° [189d] |
| **K3** zarfın kimliği | **KAYIT** | c 0.149 → 0.116 → 0.100; α 1.30 → 1.68 → 2.22; D→∞ P_der A %134, B %102; ölçülen %51 [189d] |

**H-189d ölümünün dürüst okuması:** Ölen şey, 188 ikiz haritasının derin
ikize NİCEL olarak taşınabilirliğidir; yön ve açı tutuyor. Harita, kendi
kinematiğindeki (Hkeskin) iptali 1.00'da birebir veriyordu. Derin ikizin
kinematiği değişince aynı derinlikteki iptal 1.10'da %8, 1.20'de %13 daha
zayıf. Bu yüzden kalemin sabit-kinematik öngörü tablosu iki yerden birden
aşırı iyimserdi: w^öz ve iptal genliği.

---

## TÜRETİLEN vs ÖLÇÜLEN

- **Türetilen:**
  - Kalemin öngörü tablosu (sabit-kinematik aritmetiği + 188 haritası).
    YÖNÜ doğru çıktı (H-189a), HIZI yanlış (H-189b ÖLDÜ).
  - Mekanizma öngörüsü (τ'≈1 negatif geri-beslemesi ⇒ σ_ε ↓). YÖNÜ doğru
    (H-189c); büyüklüğü türetilmemişti.
- **Ölçülen:** r_D, w^öz_HD, m_HD, σ_ε, ζ_HD, (c, α).
- **Model:** D→∞ (A/B, K0'da donmuş). "Zarfın tamamı derinlik" cümlesi
  modele dayanır; "en az yarısı" cümlesi ölçümdür.
- **Açık kalanlar:**
  - α'nın derinlikle büyüme yasası (1.30 → 2.22) türetilmedi.
  - İptal genliğinin kinematiğe duyarlılığı (EK b) sınanmadı.

---

## MANŞET (aday cümle)

> **4. demirin zarfı derinliktir.** İkiz aynı sadakatli çözücüyle τ≤1.10 ve
> τ≤1.20 merdiveniyle yeniden kurulunca gerçek/ikiz genlik oranı HAVUZ'da
> 0.910 → 0.939 → 0.956'ya yükseldi. Zarfın %32'si, sonra %51'i kapandı,
> 8/8 bantta (H-189a MÜHÜR). İkizin aralık saçılımı öngörülen negatif
> geri-beslemeyle gerçeğe doğru indi: σ_ε 0.4313 → 0.4225 → 0.4179, gerçek
> 0.4092 (H-189c MÜHÜR). Kapanış, sabit-kinematik öngörünün ancak yarı
> hızında (H-189b ÖLDÜ). Bunun iki sebebi var: derinleşen ikizin öz-terimi
> gerçeğe yürüyor (sapmanın ~yarısı), ve derin ikizin kendi iptali 188
> haritasının %13 altında kalıyor (H-189d ÖLDÜ). Zarf küçülürken kesime
> çekiliyor (c 0.149 → 0.100, α 1.30 → 2.22): α = 1.30, τ≤1.00 merdiveninin
> parmak izi. D→∞ modelleri zarfın tamamını derinliğe veriyor (A %134, B
> %102); gerçek-kinematik fark zarfa katkı vermiyor.

---

Teslim:
- **Rapor:** bu dosya.
- **Betikler** (`189_configs/`):
  - `189a_onkayit.py`: K0
  - `189b_insa.py`: ızgara sınavı + inşa sarmalayıcısı
  - `189c_zincir.py`: 184b2→187c sürücüsü + makine mührü
  - `189d_hukum.py`: H-189a..d + K3
  - `189e_figur.py`
  - `189f_ek.py`: ön-kayıtsız ek
- **Figür:** `189_derin_ikiz.png`. Sol: r_D(τ) üç derinlik + ön-mühür
  öngörüleri kesikli + 1−cτ^α + r = 1. Orta: σ_ε(D) ve |ζ_HD|(D) + harita
  z_D + gerçek çizgileri. Sağ: c(D), α(D).
- **`scratchpad/189/`:**
  - ONKAYIT_189.json
  - izgara_Hderin{110,120}.json, insa_kapi_Hderin{110,120}.json
  - K1_/OZ_/G1_proj_{Hkeskin,Hderin110,Hderin120}.npz
  - zincir_{Hkeskin,Hderin110,Hderin120}.{json,npz}, zincir_gercek.npz
  - muhur_Hkeskin.json, HUKUM_189.json, EK_189f.json
  - tüm loglar
- **Yeni ikiz dosyaları** (164 main'in AYNEN yazdığı yerler):
  - `scratchpad/164/z_Hderin{110,120}.npy` ve `insa_Hderin{110,120}.json`
  - `scratchpad/155/z_Hderin{110,120}.npy`
- Mevcut Hkeskin/gerçek önbelleklerinin üzerine yazılmadı.
