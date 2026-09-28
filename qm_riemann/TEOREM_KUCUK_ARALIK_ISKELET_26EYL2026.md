# TEOREM İSKELETİ — Yakın iki sıfır arasındaki tepe (küçük-aralık seferi, KALEM 200 öncesi)

**Tarih:** 26 Eylül 2026
**Durum:** TASLAK v2 — TÜRETİM TEFTİŞİNDEN GEÇTİ (Sonnet, 26 Eyl): matematikte hata yok;
Teorem 1–2, Sonuç 3–4, merdiven tablosu ve Gauss rakibi bağımsız olarak yeniden hesaplandı.
İşlenen düzeltmeler: (1) aritmetik sabitler yakınsamış değerlere çekildi (−0.088124 /
+0.233653); (2) W1 ana risk olarak öne alındı; (3) Gauss rakibinin χ₃ türetimi eklendi;
(4) limit değişimi (k→0 türevi) varsayım olarak yazıldı (W4); (5) Teorem 1'in Palm-ölçüsü
olarak bilinme ihtimali §9'a eklendi; (6) Lemma B'de çok-noktalı durum açıkça yazıldı.
Zeta verisine
DOKUNULMADI (körlük: hiçbir küçük-aralık M_n/δ_n² değeri hesaplanmadı). Aşağıdaki
bütün sayılar ya kapalı formdur ya da CUE Monte Carlo'sudur.
**Kaynak soru:** Not 1'in açık problemi (1) + bir okurun sorusu: "asimptotik küçük
aralıklarda (δ_n, M_n²) ortak dağılımı bulunabilir mi?"
**Betikler:** `200_configs/200t_tahmin.py` (kapalı formlar), `200t_mc.py` (CUE MC),
`200t_carpim.py` (çarpım yasası örneklemesi), `200t_sabitler.py`, `200t_merdiven.py`,
`200t_aritmetik.py`; MC çıktıları `200_configs/200t_mc*_N*.json`.

---

## 0. Hikâye (önce basit)

İki sıfır birbirine çok yaklaşınca aralarındaki tepe iki şeyle belirlenir: aralığın
karesi (tepe, aralığın karesiyle ezilir) ve **geri kalan bütün sıfırların** o noktada
ürettiği yükseklik. Rastgele matriste bunun tam cevabı var: geri kalan N−2 özdeğer,
iki özdeğerin birleştiği noktadan **|Λ|⁴ ağırlığıyla uzaklaşmış** bir CUE gibi dağılır
(dairesel Jacobi, parametre 2). Bu ağırlıklı dünyada |Λ| bağımsız çarpanların
çarpımıdır, dolayısıyla bütün momentler ve bütün log-kümülantlar kapalı formdadır.

Zeta'ya geçerken asallar ekleniyor. Sürpriz şu: eğilmiş RMT log-tepeyi **sola**
çarpık bırakıyor, asalların (Keating–Snaith aritmetik çarpanı a_k) katkısı ise
**sağa** çarpıtıyor ve yakın çiftlerde bu katkı baskın hâle geliyor. Sanı, işaret
düzeyinde sınanabilir bir tahmin veriyor: **zeta'da yakın çift tepesinin log-dağılımı
sağa çarpıktır** (L≈10'da çarpıklık ≈ +1.7; salt RMT −0.42; Gauss −0.92).

---

## 1. Kurulum ve gösterim

- U ~ Haar, U(N); özaçılar θ_1,…,θ_N. |Λ_N(e^{iθ})| := Π_j |e^{iθ} − e^{iθ_j}|
  = Π_j |2 sin((θ−θ_j)/2)| (yalnız modül kullanılıyor).
- Döngüsel sıralama θ_(1)<…<θ_(N), θ_(N+1) := θ_(1)+2π; aralıklar
  s_n = θ_(n+1) − θ_(n); tepe M_n := max_{θ∈[θ_(n),θ_(n+1)]} |Λ_N(e^{iθ})|.
- Keating–Snaith: M_n(k) := E|Λ_n(1)|^{2k} = Π_{j=1}^n Γ(j)Γ(j+2k)/Γ(j+k)².
- **Eğilmiş yasa** (dairesel Jacobi, δ=2): P̃_n(dU) := |Λ_n(1)|⁴ P_n(dU)/M_n(2) on U(n);
  X_n := |Λ_n(1)| altında P̃_n. (n=0 için X_0 ≡ 1.)
- Açılmış (unfolded) birimler: s̃ = sN/2π; Y_N := (2π/N)²·X_{N−2}/4 (M/s̃²'nin
  küçük-aralık limiti).

## 2. Teorem 1 — sabit N'de küçük-aralık koşullu ortak yasa

**Teorem 1.** N ≥ 2 sabit. Her sınırlı sürekli F: [0,1]×(0,∞) → ℝ için

  lim_{ε→0} ε^{−3} · E[ Σ_{n=1}^{N} 1{s_n<ε} F(s_n/ε, M_n/s_n²) ]
   = (1/2π) ∫_0^1 u² · E_{N−2}[ |Λ_{N−2}(1)|⁴ F(u, |Λ_{N−2}(1)|/4) ] du.

Sonuçlar:
 (a) E#{n: s_n<ε} = N²(N²−1) ε³/(72π) · (1+o(1))  [F≡1; M_{N−2}(2) = N²(N²−1)/12].
 (b) Küçük aralıklar üzerinden (Palm tipi) koşullu yasa:
     (s_n/ε, M_n/s_n²) ⇒ (U^{1/3}, X_{N−2}/4), U ~ Düzgün(0,1), **U ve X bağımsız**.
     Yani küçük aralıkta **aralık ile normalize tepe asimptotik olarak bağımsızdır.**
 (c) Momentler: her k > −5/2 için
     E[(M_n/s_n²)^{2k} | s_n<ε] → 4^{−2k} M_{N−2}(k+2)/M_{N−2}(2);
     log-momentlerin hepsi de yakınsar.
 (d) Tepenin göreli konumu (θ_max−θ_(n))/s_n → 1/2.

**İspat iskeleti.**
 1. *Weyl + simetri:* E Σ_n G = N(N−1) ∫ p_N(θ) 1{θ_2−θ_1∈(0,ε), yay içinde başka θ_j
    yok} G dθ; p_N = |Δ_N|²/((2π)^N N!). Her ardışık aralık sıralı çift (sol, sağ)
    olarak tam bir kez sayılır.
 2. *Değişken değişimi:* θ_1=φ, θ_2=φ+εu (dθ_1dθ_2 = ε dφ du). Vandermonde ayrışır:
    |Δ_N|² = 4sin²(εu/2) · Π_{j≥3}|e^{iθ_j}−e^{iφ}|²|e^{iθ_j}−e^{i(φ+εu)}|² · |Δ_{N−2}(θ′)|².
    4sin²(εu/2) = ε²u²(1+O(ε²)) ⇒ ε³ çarpanı.
 3. *Noktasal limitler* (θ′ için h.h., yani hiçbir θ_j = φ değilken): ağırlık →
    |Λ_{N−2}(e^{iφ})|⁴; yay göstergesi → 1; M/s² → |Λ_{N−2}(e^{iφ})|/4 (Lemma A).
 4. *Baskın yakınsama:* ağırlık ≤ 4^{2(N−2)}, F sınırlı ⇒ DCT. Sabit:
    N(N−1)(2π)^{N−2}(N−2)!/((2π)^N N!) = 1/(4π²); ∫dφ = 2π ⇒ 1/(2π).
    Dönme değişmezliği e^{iφ} → 1.

**Lemma A (tepe limiti).** d := min_j |e^{iθ_j} − e^{iφ}| > 0 ise, yay I=[φ, φ+s] üzerinde
 |Λ_{N−2}(mid)|·(4sin²(s/4))/s² ≤ M/s² ≤ sup_I|Λ_{N−2}| · sup_I |2sin((θ−φ)/2)||2sin((φ+s−θ)/2)|/s²,
 iki taraf da |Λ_{N−2}(e^{iφ})|/4'e gider (süreklilik; |2sin(a/2)| ≤ |a|; x(1−x) ≤ 1/4).
 Ayrıca deterministik üst sınır: **M/s² ≤ 2^{N−2}/4** (her zaman). ⇒ k ≥ 0 momentleri
 doğrudan (sınırlı sürekli fonksiyon).

**Lemma B (negatif momentler ve log için düzgün integrallenebilirlik; teftiş edilmeli).**
 s ≤ 1 için M/s² ≥ |Λ_{N−2}(mid)|/5. x_j := θ_j − mid için yay dışı |x_j| ≥ s/2 ve
 ön-limit ağırlıkta |e^{iθ_j}−e^{iθ_1}||e^{iθ_j}−e^{iθ_2}| ≲ x_j² ⇒ tekil çarpan
 |x_j|^{2k}·|x_j|⁴ ⇒ 2k+4 > −1 ⇔ k > −5/2'de integrallenebilir; η>0 payıyla de la
 Vallée-Poussin ⇒ düzgün integrallenebilirlik. Kalan |Δ_{N−2}|² sınırlı; tekillikler
 ayrı değişkenlerde. (Limit yasasında E_tilt|Λ|^{2k} < ∞ ⇔ k > −5/2 ile tutarlı; eşik Γ(2k+5)'in kutbuyla da
 örtüşür.) Birden çok θ_j aynı anda mid'e yaklaşırsa aynı sınır indis indis uygulanır:
 tekil çarpanlar ayrı değişkenlerde olduğundan çarpımları da integrallenebilir (standart,
 yazımda ayrıca gösterilecek).

## 3. Teorem 2 — limit tepenin kesin çarpım yasası

**Teorem 2.** n ≥ 1 için X_n ≐ Π_{j=1}^n |1−ξ_j|, ξ_j bağımsız:
 - ξ_1 birim çemberde, yoğunluk |1−e^{iα}|⁴/6 (dα/2π'ye göre);
 - ξ_j (j≥2) birim diskte, yoğunluk c_j (1−|z|²)^{j−2} |1−z|⁴ d²z,
   c_j = ((j−1)/π)·j(j+1)/((j+2)(j+3)).

**İspat.** Bourgade–Hughes–Nikeghbali–Yor (Duke 2008, Prop. 2.2): CUE(n) altında
Λ_n(1) ≐ Π_{j=1}^n (1−γ_j), γ_j bağımsız, |γ_j|² ~ Beta(1, j−1) (γ_1 çemberde), dönmeyle
değişmez. |Λ|⁴ ile eğmek, bağımsız çarpanların çarpımına çarpım-biçimli ağırlık
vermektir ⇒ bağımsızlık korunur, her çarpan |1−γ_j|⁴ ile eğilir. (Eşdeğer olarak
Bourgade–Nikeghbali–Rouault, IMRN 2009, Prop. 4.3, δ=2.) ∎

**Sonuç (kapalı formlar).**
 - E X_n^{2k} = M_n(k+2)/M_n(2) = Π_{j=1}^n Γ(j+2)²Γ(j+4+2k)/(Γ(j+4)Γ(j+2+k)²), k > −5/2.
 - log X_n'in kümülantları: κ_r = Σ_{j=1}^n [ψ^{(r−1)}(j+4) − 2^{1−r} ψ^{(r−1)}(j+2)].
   κ_1 = H_{n+2}+H_{n+3} − 10/3;  κ_2 = Σ[ψ′(j+4) − ψ′(j+2)/2];  κ_3 = Σ[ψ″(j+4) − ψ″(j+2)/4].
 - Genel eğim b (|Λ|^{2b}, n = N−b kalan): κ_r = Σ_{j=1}^{n}[ψ^{(r−1)}(j+2b) − 2^{1−r}ψ^{(r−1)}(j+b)]
   (b=0: sıradan CUE; b=1: bir özdeğere koşullama = ζ′(ρ) modeli; b=2: bu teorem).
   Tam-sayı olmayan N için ψ-toplamlarının analitik devamı `200t_sabitler.py`'de.

## 4. Büyük N ve okurun gözlemi

**Sonuç 3 (N→∞, açılmış birimler; Y_N = (2π/N)² X_{N−2}/4).**
 (a) E Y_N^{2k} ~ π^{4k} · 12·G(k+3)²/G(2k+5) · N^{k²}   (k > −5/2).
     k=1: E Y_N² ~ π⁴N/720.
 (b) E log Y_N → **2 log π + 2γ − 10/3 = 0.110558** (sonlu N: H_N+H_{N+1}−10/3−log4+2log(2π/N)).
 (c) Var log Y_N = ½ log N + c₂ + o(1), **c₂ = 65/18 + γ/2 − π²/2 = −1.035083**.
 (d) κ_3(log Y_N) → −2[ζ(2)−H₄^{(2)} − 4(ζ(3)−H₄^{(3)})] + ½[ζ(2)−H₂^{(2)} − 2(ζ(3)−H₂^{(3)})]
     = **−0.127077** (her N'de negatif: eğilmiş RMT log-tepesi sola çarpık).
 (e) (log Y_N − E)/√(½ log N) ⇒ N(0,1) (bağımsız toplamlar, Lindeberg).
 Yorum: tipik yakın-çift tepesi/δ̃² O(1) kalır, ama log-normal yayılımı ½ log N ile
 büyür; ikinci momentin N ile doğrusal büyümesi bu kuyruktandır.

**Sonuç 4 (küçük aralıkların ikinci momente katkısı).**
 lim_{ε→0} ε^{−7} E[ Σ_{s_n<ε} M_n² ] = M_{N−2}(3)/(224π).
 (∫u⁶du = 1/7; E X²/16 = M(3)/(16 M(2)).) Toplam E Σ M_n² ~ N² olduğundan açılmış eşik
 ε̃ altındaki aralıkların payı O(ε̃⁷): küçük aralıklar ikinci momente pratikte hiç
 katkı vermez. Okurun gözleminin kesin biçimi budur.

## 5. Sayısal doğrulama (CUE, zeta YOK)

Haar QR ile CUE; her küçük aralıkta 41 noktalı ızgara + parabolik rafinasyon.
ε̃ = açılmış eşik. "tahmin" = Teorem 1–2 kapalı formu (θ birimleri).

| N | matris | tahmin E[M/s²] | ε̃<0.2 | ε̃<0.1 | tahmin E log | ε̃<0.1 | tahmin Var log | ε̃<0.1 | oran tahmin | ε̃<0.1 |
|---|---|---|---|---|---|---|---|---|---|---|
| 6 | 6·10⁵ | 1.4784 | 1.4828±.0030 | 1.4751±.0083 | 0.3232 | 0.3229±.0062 | 0.1486 | 0.1448 | 5.570 | 5.40 |
| 10 | 4·10⁵ | 3.9279 | 3.9354±.0111 | 3.9186±.0308 | 1.2292 | 1.2330±.0081 | 0.2987 | 0.2834 | 43.77 | 43.36 |
| 16 | 3·10⁵ | 10.194 | 10.178±.035 | 10.059±.099 | 2.1007 | 2.0822±.0095 | 0.4691 | 0.4781 | 288.6 | 289.1 |
| 24 | 2·10⁵ | 23.903 | 24.130±.102 | 24.103±.286 | 2.8723 | 2.8804±.0110 | 0.6340 | 0.6326 | 1464 | 1457 |

- Aralık–tepe korelasyonu ε̃ ↓ ile → 0 (N=10: 0.021, 0.013, 0.007, 0.000). ✔ (b)
- Sonlu-ε sapması **ε² ile** gidiyor (N=10, E[M/s²]: 0.4'te +0.076, 0.3'te +0.043; oran
  1.77 ≈ (4/3)² = 1.78). Simetri (s → −s) tek dereceli terimi öldürüyor. KALEM için:
  ε̃ ≤ 0.2'de log-ortalama sapması ≈ +0.003, varyans ≈ −0.002, κ₃ ≈ +0.005.
- κ₃ (N=10): tahmin −0.0688; MC ε̃<0.1: −0.064 (8526 olay). (N=16: −0.086 vs ≈ −0.085.)
- Teorem 2 örneklemesi (ξ_j reddetme örneklemesi, 4·10⁵): N=10 E log 1.2286 (1.2292),
  Var 0.2976 (0.2987), κ₃ −0.0675 (−0.0688), E X 3.924 (3.928); N=24 hepsi ≤%0.3 içinde.

## 6. Zeta sanısı

Zeta tarafı: Z(t) Hardy fonksiyonu; ardışık sıfırlar γ_n < γ_{n+1}, δ_n, M_n (Not 1);
L = log(t/2π); δ̃_n = δ_n L/2π; **Y_n := M_n/δ̃_n²**. Küçük δ'da M_n ≈ |Z″(m_n)| δ_n²/8.

**Sanı 5.** T → ∞ (sonra ε → 0), [T, 2T]'deki sıfırlar için, {δ̃_n < ε} üzerinde:
 (i) (δ̃_n/ε, Y_n) asimptotik bağımsız; δ̃_n/ε ≐ U^{1/3}.
 (ii) E[Y_n^{2k} | δ̃_n<ε] ~ a_k · E[Y_L^{2k}] ~ a_k · π^{4k}·12·G(k+3)²/G(2k+5) · L^{k²}
      (a_k: Keating–Snaith aritmetik çarpanı; k=1'de a_1 = 1 ⇒ E[Y²|küçük] ~ π⁴L/720).
 (iii) log-kümülantlar: κ_r(log Y_n | küçük) = κ_r^{eğik}(N=L) + κ_r^{arit} + o(1),
      κ_1^{arit} = 0 ((log a)′(0) = 0),
      κ_2^{arit} = ½ Σ_p [log(1−1/p) + Li₂(1/p)] = **−0.088124**,
      κ_3^{arit} = (3/2) Σ_p S_{1,2}(1/p) = **+0.233653**
      (S_{1,2}(x) = Σ_{m≥2} H_{m−1}x^m/m², Nielsen; asallar 3·10⁶'ya kadar, yakınsamış:
      `200t_aritmetik_kesin.py`; teftiş aynı değerleri bağımsız buldu. İlk tahmin
      `200t_aritmetik.py` h=0.02 sonlu farkla −0.08814 / +0.2335 vermişti.)

**Sezgisel türetim (teftiş edilecek).**
 1. Yerel: Z(t) ≈ (Z″/2)(t−γ_n)(t−γ_{n+1}) ⇒ M_n = |Z″(m_n)|δ_n²/8·(1+O((δL)²)).
 2. Hibrit Euler–Hadamard çarpımı (Gonek–Hughes–Keating 2007): ζ ≈ P_X·Z_X. Yakın çifte
    koşullama Z_X'e Teorem 1'deki gibi etki eder (|Λ|⁴ eğimi); P_X önde gelen mertebede
    eğilmez (Landau–Gonek kaymaları O(1/L); Bui–Gonek–Milinovich 2015'in ζ′(ρ) için
    gösterdiği mekanizmanın aynısı).
 3. Bölünme: momentler çarpanlara ayrılır ⇒ a_k × eğik-CUE. **Tutarlılık:** açılmış
    birimlerde E|P_X|^{2k} ~ a_k(e^γ log X)^{k²} ve eğik-CUE(N_X = L/(e^γ log X)) katkısı
    π^{4k}·12G(k+3)²/G(2k+5)·N_X^{k²}; çarpımda X sadeleşir.
 4. Analoji: Hughes–Keating–O'Connell (2000) ζ′(ρ) sanısı = a_k × (|Λ|²-eğik CUE), aynı a_k.

**Zayıf noktalar (teftişe açık):**
 W1. **ANA RİSK.** P_X'in yakın-çift koşullamasından bağımsızlığı. Programın kendi bulgusu
     (asal dalgaları aralıklara işliyor, Not 2–4/7) O(1/L) bağlaşım olduğunu söylüyor.
     Üstelik b=2'de aritmetik terim küçük bir düzeltme DEĞİL: L≈10'da κ₃'te RMT tabanı
     −0.069, aritmetik kayma +0.234. "Önde gelen mertebe + küçük pertürbasyon" mantığının en
     az güvenilir olduğu rejim tam burası. Bu yüzden işaret tahmini bir sınavdır, garanti
     değil; b=0 ve b=1 basamakları bu riski ölçmek için kalibrasyon olarak girer (§7).
 W2. Sonlu L: Not 1 N_eff = L + c buldu (c gözlenebilire bağlı). Burada c serbest mi,
     sabit mi? (κ_1 N'e zayıf bağlı: dκ_1/dN ≈ −0.013; κ_2: +0.03/birim N.)
 W3. Limit sırası (önce T, sonra ε); sonlu-ε düzeltmesi zeta'da CUE'dekiyle aynı mı?
 W4. **Varsayım:** (ii)'den (iii)'e geçerken L→∞ limiti ile k=0'da türev almanın yer
     değiştirdiği kabul ediliyor (momentlerin k'da düzgün asimptotiği). Bu, Keating–Snaith
     moment sanısı çerçevesinden miras; burada kanıtlanmıyor.

## 6b. Sonlu L'de hangi model? (26 Eyl, teftiş sonrası)

**Bulgu 1 — a_k biçimi L≈11'de olasılık yasası değil.** φ_A(ω) = φ_eğik(ω)·a_{iω/2} ters
Fourier ile yoğunluk verilmeye çalışıldığında patlıyor (`200t_tam_yasa.py`,
`scratchpad/k200/cf_profil.py`): log|a_{iω/2}| ω ile süper-karesel büyüyor (ω=10'da +30,
ω=20'de +163; her p ≲ ω² asalı (1−1/p)^{ω²/4} ile büyür, ₂F₁ çarpanı yalnız polinom söner),
eğik CUE'nun karakteristik fonksiyonu ise N≈11'de ancak ~e^{−ω²/4·log N} kadar söner. Yani
a_k × CUE(L) yalnız sabit k'da, L→∞ önde gelen mertebede anlamlı; sonlu L'de tam dağılım
tahmini olarak KULLANILAMAZ. Düşük kümülantlar (k=0'da türevler) yine tanımlı.

**Bulgu 2 — Sonlu-L model yelpazesi.** Gerçek bir olasılık modeli: hibrit(X) = bağımsız
rastgele Euler çarpımı (p ≤ X) ⊗ eğik CUE(N_X), N_X = L/(e^γ log X) (Gonek–Hughes–Keating
eşlemesi). `200t_hibrit.py`:

| b | L | salt CUE κ₂ / κ₃ | a_k κ₂ / κ₃ | hibrit X=2 | X=3 | X=5 | X=7 |
|---|---|---|---|---|---|---|---|
| 0 | 10 | 1.940 / −2.393 | 1.852 / −2.159 | 2.126 / −2.233 | 2.080 / −2.126 | 1.996 / −2.042 | 1.977 / −1.991 |
| 1 | 10 | 0.568 / −0.194 | 0.480 / +0.039 | 0.775 / −0.039 | 0.790 / +0.050 | 0.778 / +0.103 | 0.803 / +0.131 |
| 2 | 10 | 0.299 / −0.069 | 0.211 / **+0.165** | 0.523 / **+0.083** | 0.585 / **+0.160** | 0.622 / **+0.197** | 0.674 / **+0.213** |
| 2 | 12 | 0.362 / −0.077 | 0.274 / **+0.157** | 0.580 / **+0.075** | 0.628 / **+0.151** | 0.651 / **+0.188** | 0.695 / **+0.205** |

Sonuç:
 - **b=2'de κ₃'ün işareti modelden bağımsız:** bütün aritmetik modeller +0.075…+0.213;
   salt RMT −0.07, Gauss −0.10. Fark D₃ := κ₃ − κ₃^{CUE}(L) aritmetik ailede [+0.15, +0.28],
   salt RMT'de 0, Gauss'ta ≈ −0.03. **Birincil sınav adayı: D₃.**
 - **b=2'de κ₂ modele çok bağlı** (0.21…0.70): doğrulayıcı sınav olamaz; model ayırıcıdır.
 - b=0 basamağında modeller κ₂'de ~0.27 aralığa yayılıyor; 5·10⁴ rastgele t ile
   SE(κ₂) ≈ 0.013 ⇒ **alt basamaklar sonlu-L modelini seçebilir**, seçilen model b=2'ye
   kör tahmin olarak taşınır (merdiven sanısı).
 - Çarpıklık (κ₃/κ₂^{3/2}) modele bağlı κ₂ yüzünden birincil olmamalı; κ₃ (ya da D₃) olmalı.

## 7. Eğim merdiveni (kalibrasyon fikri)

Koşullanan sıfır sayısı b = 0, 1, 2 ⇒ geri kalan |Λ|^{2b} ile eğilir:
 - b=0: log|Z(t)|, rastgele t (Keating–Snaith 2000; bilinen fizik),
 - b=1: log|Z′(γ_n)| (HKO 2000),
 - b=2: log(|Z″(m_n)|/2) yakın çiftlerde (bu sefer).
**Merdiven sanısı:** κ_r(zeta, b) − κ_r(CUE_b, N=L) = κ_r^{arit}, **üç basamakta aynı.**
Alt iki basamak O(1/L) düzeltmelerini kalibre eder; yeni olan b=2.

Tahmin tablosu (`200t_merdiven.py`; κ₂, κ₃ ve çarpıklık κ₃/κ₂^{3/2}):

| b | L | CUE κ₂ | CUE κ₃ | sanı κ₂ | sanı κ₃ | çarpıklık CUE → sanı |
|---|---|---|---|---|---|---|
| 0 | 10 | 1.9403 | −2.3925 | 1.8522 | −2.1589 | −0.885 → −0.856 |
| 0 | 24.47 | 2.3874 | −2.4368 | 2.2993 | −2.2031 | −0.661 → −0.632 |
| 1 | 10 | 0.5681 | −0.1944 | 0.4800 | +0.0393 | −0.454 → +0.118 |
| 1 | 24.47 | 0.9600 | −0.2312 | 0.8719 | +0.0025 | −0.246 → +0.003 |
| 2 | 10 | 0.2987 | −0.0688 | 0.2106 | +0.1649 | −0.421 → **+1.706** |
| 2 | 12 | 0.3617 | −0.0765 | 0.2735 | +0.1571 | −0.352 → **+1.098** |
| 2 | 24.47 | 0.6423 | −0.0995 | 0.5541 | +0.1341 | −0.193 → **+0.325** |

Mesaj: aritmetik katkı b=0'da %10'luk bir düzeltme, b=1'de çarpıklığı sıfıra çeker,
**b=2'de işaretini çevirir.** Yakın çiftler asalları görünür kılıyor.

## 8. Rakip modeller

 - **R_RMT** (salt eğik CUE, κ^{arit}=0): b=2'de κ₃ < 0 (−0.07…−0.13).
 - **R_Gauss** (aynı güç spektrumlu durağan Gauss yüzeyi, Not 1'in Null 1'i).
   Türetim: iki noktalı Kac–Rice yoğunluğu E[|X′(t₁)||X′(t₂)| | X(t₁)=X(t₂)=0]·p(0,0).
   s→0'da iki sıfır koşulu X(mid)≈0, X′(mid)≈0'a dönüşür ve X′(t₁,₂) ≈ ∓X″·s/2, yani
   Jacobien |X′(t₁)X′(t₂)| ≈ X″²s²/4. Durağan süreçte X″ koşullu olarak Gauss:
   Var(X″ | X=0, X′=0) = σ² = λ₄ − λ₂²/λ₀ (X′ ⟂ X, X″). Merkezli Gauss'u x² ile eğmek
   χ₃ verir ⇒ |X″|/σ ~ χ₃; tepe M/s² = |X″|/8·(1+O(s²)) olduğundan log-tepenin şekli
   log χ₃'ünkü (σ yalnız κ₁'i kaydırır). χ_k için κ_r(log) = ψ^{(r−1)}(k/2)/2^r ⇒
   Var log = ψ′(3/2)/4 = 0.2337,
   κ₃ = ψ″(3/2)/8 = −0.1036, **L'den bağımsız**. (Dikkat: L≈10–12'de varyansı sanıyla
   neredeyse aynı: 0.234 vs 0.21–0.27. Ayırt edici olan κ₃'ün işareti ve varyansın
   ½ log L ile büyümesi.)
 - **R_1eğim** (yanlış koşullama, |Λ|²): κ_1 ve κ_2 farklı.

Güç kaba tahmini: n olayda SE(κ₃) ≈ √(6κ₂³/n); κ₂ = 0.25, n = 2000 ⇒ 0.007;
n = 14 000 (zeros6, δ̃<0.2) ⇒ 0.003. Sanı ile rakipler arasındaki κ₃ farkı ≈ 0.23.
**Düzeltme (200-A teftişi):** bu dar formül κ₄, κ₆'yı atıyor; bu çarpık ailede gerçek SE
~2.7× büyük (ürün yasası simülasyonu). Karar hesapları jackknife SE ile yapılır; ayrım yine
geniş (0.23 fark vs ~0.01–0.02 SE).

## 9. Literatür konumu (26 Eyl ajan taraması + 25 Eyl taraması)

 - Ben Arous–Bourgade (Ann. Probab. 41 (2013) 2648–2681): CUE en küçük aralıklar,
   N^{−4/3} ölçeği, Poisson, yoğunluk u²/(24π). Yöntem korelasyon fonksiyonları;
   **işaret (tepe/polinom değeri) yok, |Λ|⁴ koşullama yok.** Teorem 1(a) onların sabitiyle
   tutarlı (N⁴ε³/(72π) = ∫u²/(24π)).
 - Feng–Wei (Ann. Probab. 49 (2021) 997–1032), Feng–Tian–Wei (GAFA 29 (2019)),
   Bourgade (JEMS 24 (2022)): marjinal küçük-aralık yasaları; işaretli süreç yok.
 - Bourgade–Hughes–Nikeghbali–Yor (Duke 145 (2008) 45–69), Bourgade–Nikeghbali–Rouault
   (IMRN 2009(23) 4357–4394): çarpım yasası ve dairesel Jacobi — **yapı taşı**; birleşen
   özdeğer koşullamasıyla bağ kurulmamış.
 - Hughes–Keating–O'Connell (Proc. R. Soc. A 456 (2000) 2611–2627), Bui–Gonek–Milinovich
   (Forum Math. 27 (2015) 1799–1828, arXiv:1302.5032): b=1 basamağı, aynı a_k.
 - 25 Eyl taraması: Csordas–Smith–Varga, Stopple, Rodgers–Tao, Farmer–Gonek–Lee,
   Assiotis–Keating–Warren, Dehaye, Charlier–Claeys, Conrey–Ghosh, HLPC24, PC24 — hiçbiri
   küçük-aralık koşullu (δ, M) yasasını vermiyor.

**Dürüst yenilik iddiası:** Teorem 1 temel (Weyl + DCT); |Λ|⁴ Palm eğimi uzmanlarca
bilinen bir olgu olabilir (belirleyici süreçlerin Palm ölçüleri; Forrester, *Log-Gases and
Random Matrices*, 2010 çerçevesi — sayfa düzeyinde henüz doğrulanmadı, makaleden önce
bakılacak), ama (aralık, tepe) yasası olarak yazılmış bulunamadı —
"bildiğimiz kadarıyla" denmeli. Teorem 2, BHNY/BNR'nin doğrudan sonucu (öyle sunulmalı).
Yeni olan: Sanı 5, merdiven sanısı ve **b=2'de çarpıklığın işaret dönmesi** tahmini.

## 10. KALEM 200 için girdiler (tasarım, teftişten sonra)

 - Veri: `128_odl_zeros6_2e6_zeros.npz` (δ̃<0.1: 1824; <0.2: 14 279; <0.3: 47 458; L≤12.1),
   `55_win_1e+08…11` (L = 16.6–23.5, M zaten hesaplı, δ̃<0.2'de ~200–260/pencere),
   `53` (L=24.47, 83). Motor: yakın çiftte Z orta noktada + küçük ızgara (mpmath, yüksek
   hassasiyet; Motor B mutlak hatası 7.6e−6 δ̃<0.05'te sınırda).
 - Birincil gözlenebilirler: log Y'nin κ₁, κ₂, κ₃'ü (ε̃<0.2), L-pencerelerinde; bağımsızlık
   (aralık–tepe Spearman'ı → 0).
 - Kalibrasyon basamakları: b=0 (rastgele t'de log|Z|), b=1 (log|Z′(γ_n)|).
 - Sonlu-ε: aynı ε'de CUE MC ile model-veri karşılaştırması (ε² düzeltmesi küçük).
 - N_eff: κ₁'den (aritmetiksiz) N_eff, sonra κ₂ ve κ₃ tahmin (bir parametre, iki sınav).
 - Rakipler: R_RMT, R_Gauss (faz-rastgeleleştirilmiş yüzeyde aynı analiz), R_1eğim.
 - M6: sentetik güç (CUE N_eff + bağımsız Euler çarpanı örneklemesi ⇒ sanı altında sahte
   veri; salt eğik CUE ⇒ rakip altında sahte veri).
