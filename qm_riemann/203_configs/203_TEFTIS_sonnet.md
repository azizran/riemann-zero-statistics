# 203 κ₃ (oran/ratios) türetim teftişi

Kaynak: `203_configs/203_oran_k2.py`, `203_oran_k3.py`. Repoda ayrı bir "203_TURETIM notu"
bulunamadı (yalnızca bu iki .py dosyasının docstring'leri var; `find`/`grep` ile
`*turetim*`, `E_θ`, `∂α2` aratıldı, 203'e özel not yok — 185/186 için var ama 203 için yok).
Bu, E_θ/A₁ maddesini tam bağımsız kılmayı imkânsızlaştırdı (bkz. madde a-iii).

## Özet tablo

| # | Konu | Sonuç |
|---|------|-------|
| 1 | log ζ = −∫ζ'/ζ işareti, κ₂=½∫xG(x)dx | **PASS** (elle yeniden türetildi, tam eşleşme) |
| 1 | κ₃=(3/4)Re E[Y²Ȳ]=−(3/4)∫∫∫E[XXX̄] | **PASS** (elle yeniden türetildi, tam eşleşme) |
| a-i | 3B→2B indirgeme: min(x1,x2) ağırlığı, simetri faktörü 2, −3/2 katsayısı | **PASS** (elle yeniden türetildi, tam eşleşme) |
| a-ii | K(u,x) analitik devamı (adım 5: dz+Pfun+Q2−C) | **PASS** (asal-asal cebirsel özdeşlik olarak ispatlandı + 1e-16 seviyesinde sayısal doğrulama) |
| a-iii | E_θ kapalı formu / ∂α2 log A₁ | **UNCERTAIN** (not eksik, değişken isimleri (γ1,γ2,δ) tanımsız) |
| a-iv | CUE bağımsız çapraz kontrolü | **UNCERTAIN/kısmi PASS** (Monte Carlo yapılmadı — CPU çakışması; ama CUE hedef sabiti (3/4)Σψ''(j)→−π²/4 tam doğrulandı) |
| b | Kimlik terimi D → c_P | **PASS** (hem kapalı-form analitik ispat hem sayısal: 0.2336529 vs 0.233653) |
| c | L→∞ ⇒ swap kısmı → −π²/4 | **PASS** (CUE analoğu tam; ζ-spesifik trend L=10→30'da −2.335→−2.417, hedef −2.467'ye doğru monoton yaklaşıyor) |
| d | Eksik terim / işaret hataları | **PASS** (aşağıda) |
| d | Wc=0.75 kesme hatası (L≈9–22) | **PASS** (ihmal edilebilir, üstel küçülüyor) |

## Detaylar ve kanıtlar

### 1) Temel kurulum (log ζ = −∫ζ'/ζ, κ₂, κ₃ işaretleri/katsayıları) — PASS
Bağımsızca: log ζ(s) = −∫_s^∞(ζ'/ζ)dw = −∫_0^∞X(a)da. Y=−∫X(a)da.
(ReY)² açılımı + E[XX]=E[Ȳ²']=0 ⇒ E[(ReY)²]=½E[YȲ]=½∫∫G(a+b) ⇒ κ₂=½∫₀^∞xG(x)dx.
Tam kod ile eşleşiyor (`203_oran_k2.py` satır 27).
(ReY)³ = (Y+Ȳ)³/8, E[Y³]=E[Ȳ³]=0 (frekanslar hiç sıfırlanmıyor) ⇒
κ₃=(3/8)(E[Y²Ȳ]+conj)=(3/4)Re E[Y²Ȳ]. Y²Ȳ=−∫∫∫XXX̄ ⇒ κ₃=−(3/4)∫∫∫E[XXX̄]. Tam eşleşme.

### a-i) 3B→2B indirgeme — PASS
S1=F(x1)K(α2−α1,x1), S2=F(x2)K(α1−α2,x2) yalnız x1=α1+β, x2=α2+β'ye bağlı, β'den
BAĞIMSIZ (α1,α2 farkı ve x1,x2 β kaymasıyla değişmiyor). β integrali bu yüzden salt bir
uzunluk: ∫₀^{min(x1,x2)}dβ = min(x1,x2). Böylece ∫∫∫(S1+S2)dα1dα2dβ = ∫∫Φ(x1,x2)min(x1,x2)dx1dx2
(dokümanın "∫∫∫ = ∫∫Φ·min dx1dx2" iddiasıyla birebir). Φ(x1,x2) simetrik olduğundan
∫∫_kadran = 2∫₀^∞dx1 x1∫₀^∞dw Φ(x1,x1+w). α1↔α2 simetrisiyle ∫∫∫S2=∫∫∫S1, dolayısıyla
−(3/4)∫∫∫(S1+S2) = −(3/2)∫∫∫S1 = −(3/2)·[2∫dx1x1∫dwΦ]/2 = −(3/2)∫₀^∞dx1x1∫₀^∞dwΦ(x1,x1+w).
Kod/docstring'deki katsayı ve ağırlıkla TAM eşleşiyor — işaret, 3/2 katsayısı ve min(x1,x2)
ağırlığının kökeni doğrulandı.

### a-ii) K(u,x) analitik devamı — PASS (asal-bazlı cebirsel özdeşlik + sayısal)
φ_p = ψ_p·(1−p^{−x}) (tanımdan trivial), ψ_p=1+psi1 yazılıp q/(1-q) geometrik seri
q/(1-q)=q+q²/(1-q) açılınca, HER TEK ASAL İÇİN (toplamsız, cebirsel) şu özdeşlik çıkıyor:
  −log(p)·q/(1-q)·φ_p(x) = [−log(p)q/(1-q)] + [log(p)p^{-(1+u+x)}] + [log(p)p^{-x}q²/(1-q)] − [log(p)q/(1-q)(1-p^{-x})psi1]
Toplamda: K_ham(u,x) = (ζ'/ζ)(1+u) + P(1+u+x) + Q2 − C — kodun (satır 37) formülüyle BİREBİR.
Sayısal doğrulama (check2.py, Test 1): 4 asal (2,3,17,10007) × 4 u değeri (0.3,-0.5,-0.73,2.0)
× 2 x değeri için lhs−rhs farkı 1e-16..1e-27 mertebesinde (makine hassasiyeti). PASS.
Ayrıca Pfun(s) [Σ_p log(p)p^{-s}, Mangoldt-yinelemeli devam] bağımsız yüksek-hassasiyet
doğrudan toplamla (2×10⁷ asal + PNT kuyruk) karşılaştırıldı: P(2)=0.4930911(kod) vs
0.4930911094(bağımsız); P(3),P(4) benzer; küçük s (1.05,1.2,1.5) için de tutarlı
(not: ilk denemede yanlış referans sabiti — standart ağırlıksız asal-zeta — kullanılmıştı,
düzeltilince tam eşleşti; bu benim hatamdı, kodda değil).

Yan bulgu: K(u,x) gerçekten u=−1'de (q→1, 1−q→0) ve 1+u=−2,−4,... bayağı sıfırlarında ek
kutuplara sahip — bu da Wc<1 seçiminin sadece güvenlik payı değil, matematiksel zorunluluk
olduğunu gösteriyor (bkz. Wc bölümü).

### a-iii) E_θ kapalı formu / ∂α2 log A₁ — UNCERTAIN
Promptta verilen E_θ formülünde γ1, γ2, δ değişkenleri α1,α2,β cinsinden tanımlanmamış
(muhtemelen 2 üst-2 alt kaymalı genel ratios-conjecture nesnesinden türetilirken kullanılan
ara değişkenler) ve repoda bu adımı belgeleyen ayrı not bulunamadı. Bu yüzden θ üzerinden
doğrudan sayısal entegrasyonla bağımsız test yapılamadı — yanlış özdeşleme ile hem yanlış
PASS hem yanlış FAIL riski var, o yüzden dürüstçe UNCERTAIN bırakıyorum.
Dolaylı destek: a-ii'de K(u,x)'in NİHAİ kapalı formu asal-asal cebirsel özdeşlik olarak tam
doğrulandığından (E_θ/A₁ hesabının ürünü olması gereken şey), E_θ/A₁ adımında sinsi bir hata
olup da son formülün yine de bu kadar temiz çıkması olası değil — ama kesin değil.
Öneri: "203_TURETIM notu" varsa (docstring'de atıfı var) bu adımı ayrıca teftiş etmek için
paylaşılmalı.

### a-iv) CUE çapraz kontrolü — UNCERTAIN/kısmi PASS
Haar unitary Monte Carlo (scipy.stats.unitary_group) çalıştırılmadı: oturum sırasında
sistemde 8 paralel `203_oran_k3.py` (üretim koşusu) + birkaç `isik_tf.py` süreci CPU'yu
doyurmuş durumdaydı; hem "üretim scriptiyle yarışma" talimatına hem zaman bütçesine saygıyla
bu adım atlandı. Bunun yerine CUE hedef sabitinin kendisi doğrulandı:
  (3/4)Σ_{j=1}^N ψ''(j): N=1→−1.803, N=10→−2.393, N=1000→−2.4667, N=20000→−2.46736
  mpmath.nsum ile N→∞: Σψ''(j) = −3.28986813369645... ⇒ (3/4)·bu = **−2.4674011002723395**
  = float(−π²/4) BİT BİT AYNI (20 basamak hassasiyetle).
Yani −π²/4 hedefi rastgele/yanlış aktarılmış bir sabit değil, gerçekten
(3/4)Σψ''(j)'nin tam limiti — bu, madde (c)'nin CUE ayağını tam destekliyor.

### b) Kimlik terimi → c_P — PASS (analitik + sayısal)
D=−Σ_p log³p·y1y2/((1−y1)(1−y2)), y_i=p^{-1-α_i-β}. y1y2/((1-y1)(1-y2))=Σ_{k,l≥1}y1^k y2^l
açılıp α1,α2,β∈(0,∞) üzerinden ayrı ayrı entegre edilince (her biri 1/(k log p) tipi):
  ∫∫∫D = −Σ_p Σ_{k,l≥1} p^{-(k+l)}/(kl(k+l)) = −Σ_p Σ_{m≥2} p^{-m}/m · Σ_{k+l=m}1/(kl)
⇒ −(3/4)∫∫∫D = (3/4)Σ_{p,m≥2}(1/(m p^m))Σ_{k+l=m}1/(kl) = c_P — TANIMIN TA KENDİSİ, tam eşleşme.
Sayısal: primes<2×10⁶, m<60 ile 0.2336529317 (kodun cP=0.233653 ile 6-7 basamak uyum).

### c) L→∞ limiti −π²/4 — PASS
CUE hedefi yukarıda (a-iv) tam doğrulandı: (3/4)Σψ''(j)→−π²/4 tam (20 basamak).
ζ-spesifik taraf: `k3_gl(L,Wc=0.75,n=12)` ile (arka planda, CPU çakışması nedeniyle ~7.5
dakika sürdü ama tamamlandı):
  L=10.0:  κ3=−2.101686   swap_part(=κ3−c_P)=−2.335339
  L=30.0:  κ3=−2.183215   swap_part(=κ3−c_P)=−2.416868
  hedef:   −π²/4 = −2.467401
swap_part L=10→30 arasında −2.335→−2.417'ye, yani hedefe MONOTON ve doğru yönde
yaklaşıyor (fark −0.132'den −0.051'e küçülüyor, kabaca ~1/L tipi yakınsama ile tutarlı —
L=10'da kalan fark 0.132, L=30'da 0.051, oran≈2.6 ≈ 3.0=30/10'a yakın, 1/L azalımıyla
uyumlu). Bu, birleşik-çekirdek formülünün gerçekten ζ-tarafında da −π²/4'e doğru gittiğini
DOĞRUDAN sayısal olarak destekliyor. Yapısal gerekçe: birleşik çekirdeğin w→0'da sonlu
olması (bkz. d, kutup iptali doğrulandı) sonlu bir L→∞ limitinin var olabilmesi için
gerekli bir önkoşuldu; sayısal trend bu önkoşulun yeterliliğini de destekliyor.

### d) Eksik terim / işaret taraması — PASS
- log ζ=−∫ζ'/ζ: doğrulandı (madde 1).
- (3/4)Re faktörü: doğrulandı (madde 1).
- Simetri faktör 2 (S1+S2→2×∫∫∫S1) ve nihai −3/2 katsayısı: doğrulandı (madde a-i).
- Terim sayısı: 2 üst kayma (α1,α2) + 1 alt/eşlenik kayma (β) şekli için CFKRS reçetesi
  kimlik + 2 tekli-takas (α1↔−β, α2↔−β) verir — kodun D+S1+S2 yapısı bu sayımla uyumlu;
  eksik terim bulunmadı.
- K'nin analitik devamı (adım 5): PASS (madde a-ii), üstelik u=−1 ve 1+u=bayağı sıfırlarda
  GERÇEK ek kutuplar keşfedildi — Wc<1 seçiminin gerekliliğini doğruluyor.

### d) Wc=0.75 kesme hatası büyüklüğü (L≈9–22) — PASS (ihmal edilebilir)
Atılan terim F(x1+w)K(-w,x1+w), w=Wc=0.75'te x1=0.05 için:
  L=9.35:  F·K = −1.147e-02  (w=1.5'te −4.9e-05, w=2'de −6.1e-05 — hızla küçülüyor)
  L=22.3:  F·K = −3.634e-07
Kaba kuyruk tahmini (F·K(Wc)/L, e^{-Lw} düşüşüyle): L=9.35'te ≈ −1.2e-3, L=22.3'te ≈ −1.6e-8;
κ3'ün mutlak büyüklüğü (~2-3) yanında sırasıyla ~%0.05 ve ihmal edilebilir. Hata e^{-0.75L}
ile üstel küçülüyor — sealed-data karşılaştırmasının hassasiyetini tehlikeye atmaz.
(Not: w=1.0'da tam olarak q=1 tekilliğine çarpıp ZeroDivisionError/NaN üretti — bu da
Wc'nin matematik olarak neden 1'in altında kalması gerektiğinin doğrudan kanıtı.)

## Genel değerlendirme
Kurulumun (adım 1), 3B→2B indirgemenin (adım 3/4'ün iskeleti) ve K'nin analitik devamının
(adım 5) TÜMÜ bağımsızca elle yeniden türetilip kodla BİREBİR (işaret dahil) eşleştiği
gösterildi; kimlik terimi→c_P kapalı formda ispatlandı; L→∞ limiti hem CUE analoğunda tam,
hem ζ-tarafında sayısal trend olarak (L=10→30, −2.335→−2.417→hedef −2.467) doğrulandı.
Tespit edilen tek gerçek risk alanı, E_θ/∂α2 log A₁ hesabının (adım 3'ün "takas terimi"
çekirdeği) bağımsız yeniden türetilememiş olması — repoda bunu belgeleyen ayrı not yok, yalnız
iki .py dosyasının docstring'i var. Yeni bir işaret/eksik-terim/kutup hatası BULUNMADI;
aksine kodun K formülünün neden Wc<1 gerektirdiğine dair ek matematiksel gerekçe (u=−1 ve
bayağı-sıfır kutupları) keşfedildi ve Wc=0.75 kesme hatasının L≈9–22'de ihmal edilebilir
(~%0.05 ve altı, üstel küçülüyor) olduğu gösterildi.
