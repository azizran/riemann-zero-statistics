# Bogomolny–Keating aritmetik sanısı: ölçülen "uydu" yapısını öngörüyor mu?

Tarih: 24 Eylül 2026 — Hazırlayan: literatür teftişi (tam-metin doğrulamalı)

---

## 0. Erişim durumu (dürüstlük notu)

| # | Kaynak | Erişim |
|---|---|---|
| 1 | E. Bogomolny, "Riemann zeta function and quantum chaos", arXiv:0708.4223 | **TAM METİN OKUNDU** (PDF → pypdf, 18 sayfa) |
| 2 | Bogomolny & Keating, PRL 77, 1472 (1996) | **ERİŞİLEMEDİ** — APS paywall; sadece özet (abstract) ve künye doğrulandı |
| 3 | Bogomolny & Keating, Nonlinearity 8 (1995) 1115 / Nonlinearity 9 (1996) 911 | **ERİŞİLEMEDİ** — IOPscience bot-koruması (Radware CAPTCHA) tam metni engelledi; talimat gereği CAPTCHA aşılmadı. Sadece künye (CrossRef DOI kaydı + Bristol kurumsal sayfası) doğrulandı |
| 4 | Conrey & Snaith, "Applications of the L-functions ratios conjectures", arXiv:math/0509480 | **TAM METİN OKUNDU** (PDF → pypdf, 59 sayfa) |
| 5 | Hardy–Littlewood tekil serisinin Ramanujan açılımı | Ayrı bir kaynağa gerek kalmadı — **kaynak 1 içinde (eş. 6.8) birebir bu formda mevcut** ve Hardy–Littlewood 1923'e atıflı; orijinal Hardy–Littlewood (1923, *Acta Math.* 44, 1–70, "Partitio Numerorum III") makalesinin kendisine erişilemedi (çok eski, açık erişim değil) |

PRL 1996 ve Nonlinearity I/II'nin **tam metnine** erişilemediği için, o üç kaynağa atfedilen formüller bu raporda **doğrudan alıntılanmıyor**. Bunun yerine, kaynak 1 (Bogomolny 2007 derlemesi — PRL 1996'nın bizzat yazarlarından biri tarafından yazılmış, ilgili denklemleri açıkça "ref. 7) = PRL 1996"ya atfederek yeniden sunan bir derleme) kullanıldı; alttaki tüm formüller bu tam-metinden birebir aktarılmıştır, denklem numaralarıyla.

---

## 1. Doğrulanmış künyeler

1. **E. Bogomolny**, "*Riemann zeta function and quantum chaos*", arXiv:0708.4223 [nlin.CD] (30 Ağustos 2007); **Prog. Theor. Phys. Suppl. 166 (2007), 19–44**. Tek yazarlı. (Osaka'daki "Quantum Mechanics and Chaos" konferansı, Eylül 2006, davetli konuşma metni.) — arXiv özet sayfasından ve PDF'nin ilk sayfasından doğrulandı.
2. **E. B. Bogomolny, J. P. Keating**, "*Gutzwiller's Trace Formula and Spectral Statistics: Beyond the Diagonal Approximation*", **Phys. Rev. Lett. 77 (1996), 1472–1475**. DOI: 10.1103/PhysRevLett.77.1472. Yayın tarihi 19 Ağustos 1996. — APS özet sayfasından doğrulandı.
3. **E. B. Bogomolny, J. P. Keating**, "*Random matrix theory and the Riemann zeros I: three- and four-point correlations*", **Nonlinearity 8 (1995), 1115–1131**. DOI: 10.1088/0951-7715/8/6/013. — CrossRef kaydından ve Bristol Üniversitesi kurumsal araştırma sayfasından doğrulandı (yazarlar: E B Bogomolny, J P Keating).
   **Part II**: "*...II: n-point correlations*", **Nonlinearity 9 (1996), 911–935**. DOI: 10.1088/0951-7715/9/4/006 (published-print: 1996-07-01). — CrossRef kaydından doğrulandı. (Not: Bogomolny'nin 2007 derlemesindeki kendi referans listesi Part II'nin yılını yanlışlıkla "1995" olarak yazmış; CrossRef'e göre doğru yıl **1996**'dır.)
4. **J. B. Conrey, N. C. Snaith**, "*Applications of the L-functions ratios conjectures*", arXiv:math/0509480 [math.NT] (21 Eylül 2005); **Proc. London Math. Soc. 94:3 (2007), 594–646**. DOI: 10.1112/plms/pdl021. — arXiv ve CrossRef'ten doğrulandı.
5. **Hardy–Littlewood tekil serisi**: G. H. Hardy, J. E. Littlewood, "*Some problems of 'Partitio numerorum'; III: On the expression of a number as a sum of primes*", **Acta Mathematica 44 (1923), 1–70**. — Doğrudan erişilemedi; künye Bogomolny (2007)'nin kendi referans listesinden (ref. 6) alınıp aktarılmıştır. Formülün kendisi (eş. 6.8, aşağıda) Bogomolny (2007)'nin tam metninden doğrudan okundu.

---

## 2. (a) BK'nın R₂(ε) formülü — tam metinden, denklem numaralarıyla

Kaynak: Bogomolny (2007), §5–§7 (arXiv:0708.4223).

**Trace formülü / osilasyon terimi** (eş. 3.6):
$$d^{(osc)}(E) = -\frac{1}{\pi}\sum_{n=1}^\infty \frac{\Lambda(n)}{\sqrt n}\cos(E\ln n),\qquad \bar d(E)=\frac{1}{2\pi}\ln\frac{E}{2\pi}\ \ (3.5)$$

**Bağlantılı iki-nokta korelasyonu** (eş. 6.1) üç parçaya ayrılır (eş. 7.7):
$$R_2(\epsilon) = \bar d^2(E) + R_2^{(diag)}(\epsilon) + R_2^{(off)}(\epsilon)$$

**Köşegen (diagonal) terim** (eş. 6.2–6.4):
$$R_2^{(diag)}(\epsilon) = -\frac{1}{4\pi^2}\frac{\partial^2}{\partial\epsilon^2}\ln\!\Big(|\zeta(1+i\epsilon)|^2\,\Phi^{(diag)}(\epsilon)\Big),\qquad
\Phi^{(diag)}(\epsilon)=\exp\Big(\sum_p\sum_{m=1}^\infty \frac{1-m}{m^2p^m}e^{im\ln p\,\epsilon}+c.c.\Big)$$

**Köşegen-dışı (off-diagonal) terim — Hardy–Littlewood sanısı üzerinden** (eş. 6.5–7.6): $\Lambda(n_1)\Lambda(n_2)$ çarpımının ortalaması
$$\alpha(r)=\lim_{N\to\infty}\frac1N\sum_{n=1}^N \Lambda(n)\Lambda(n+r)$$
Hardy–Littlewood **tekil serisi** ile verilir (eş. 6.8):
$$\alpha(r)=\sum_{(p,q)=1}e^{2\pi i pr/q}\Big(\frac{\mu(q)}{\phi(q)}\Big)^2$$
(Bu toplam, $p<q$, $(p,q)=1$ üzerinden alınıyor — tanım gereği bu ifade tam olarak **Ramanujan toplamı** $c_q(r)=\sum_{(p,q)=1,0<p<q} e^{2\pi i pr/q}$'dir; metin bu terimi açıkça "Ramanujan toplamı" diye adlandırmıyor ama ifade birebir bu standart tanıma denk düşüyor.) Çift $r$ için kapalı form (eş. 6.11–6.12):
$$\alpha(r)=C_2\prod_{p\mid r,\,p>2}\frac{p-1}{p-2},\qquad C_2=2\prod_{p>2}\Big(1-\frac{1}{(p-1)^2}\Big)\approx 1.32032$$

Bunlar toplanınca (eş. 7.1–7.5) **kapalı form**:
$$R_2^{(off)}(\epsilon)=\frac{1}{4\pi^2}|\zeta(1+i\epsilon)|^2\,e^{2\pi i\bar d\epsilon}\,\Phi^{(off)}(\epsilon)+c.c.,\qquad
\Phi^{(off)}(\epsilon)=\prod_p\Big(1-\frac{(1-p^{i\epsilon})^2}{(p-1)^2}\Big)\ \ (7.5\text{–}7.6)$$

**Kırpılmış (unfolded) hâli** (eş. 7.8): $R_2(\varepsilon)=\bar d(E)^{-2}R_2(\varepsilon/\bar d(E))$, ve $E\to\infty$ limitinde köşegen-dışı terim GUE osilasyonuna indirgenir (eş. 7.9): $R_2^{(off)}(\varepsilon)\to (2\pi\varepsilon)^{-2}(e^{2\pi i\varepsilon}+e^{-2\pi i\varepsilon})$.

**İlk-mertebe düzeltmeler dahil açık asimptotik** (eş. 7.10–7.11):
$$R_2(\varepsilon)=1-\frac{\sin^2\pi\varepsilon}{\pi^2\varepsilon^2}-\frac{\beta}{\pi^2\bar d^2}\sin^2\pi\varepsilon-\frac{\delta}{2\pi^2\bar d^3\varepsilon}\sin 2\pi\varepsilon+O(\bar d^{-4})$$
$$\beta=\gamma_0^2+2\gamma_1+\sum_p\frac{\ln^2p}{(p-1)^2}\approx1.57314,\qquad \delta=\sum_p\frac{\ln^3p}{(p-1)^2}\approx2.3157$$

Bu formül Odlyzko'nun $2\cdot10^8$ sıfırlık verisiyle karşılaştırılmış ve fark saf istatistiksel (yapısız) bulunmuştur (Bogomolny 2007, Şekil 5).

---

## 3. (b) Fourier ikilisinde (form faktör / F(α)) ince yapı var mı?

**Doğrudan metinde AÇIKÇA YOK.** Kaynak 1 ve kaynak 4'te $R_2(\epsilon)$ (ya da onun ratios-conjecture eşdeğeri, aşağıda) **her zaman $\epsilon$ (enerji/seviye-farkı) uzayında, sürekli bir fonksiyon olarak** yazılıyor. Ne "form faktörü $K(\tau)$" ne de Montgomery'nin $F(\alpha)$ notasyonu, kaynak 1'in tam metninde (18 sayfa, arandı: "form factor", "$K(E_i,E_j)$" sadece GUE çekirdeği eş. 4.4 için kullanılmış, $F(\alpha)$ hiç geçmiyor) bir kez bile geçmiyor. Conrey–Snaith (kaynak 4) de aynı şekilde $\epsilon$-eşdeğeri değişkenle ($r$, sonra ölçeklenmiş $y=rL/2\pi$) çalışıyor (Teorem 4.1, eş. 4.25–4.30) ve sonucu doğrudan Montgomery'nin $1-(\sin\pi y/\pi y)^2$ limitine bağlıyor (satır 1929: *"The expression $1-(\sin 2\pi y)/(\pi y)^2$ is exactly the limiting two-point correlation function predicted by Montgomery"*) — yine $\alpha$-uzayında değil.

**Ancak formülden türetilebilir** — kısa taslak (bu benim türetimim, metinde yazılı değil):

$G(\omega)=\langle e^{i\omega t_n}\rangle$ tipi bir "yapı çarpanı" ile $|G(\omega)|^2$ arasındaki standart ilişki, tam olarak $\sum_{n,m} e^{i\omega(t_n-t_m)}$'dir — yani **çift toplamı**, ki bu Montgomery'nin $F(\alpha)=\sum_{\gamma,\gamma'}T^{i\alpha(\gamma-\gamma')}w(\gamma-\gamma')$ tanımıyla (ölçek/normalizasyon farkıyla) aynı nesnedir; $F(\alpha)$ zaten $R_2(\epsilon)$'un Fourier dönüşümüdür. Şu adımlar BK'nın (7.5)–(7.6) formülünün $\epsilon\to$ Fourier-ikili uzayına taşındığında neden ayrık uydu yapısı vereceğini gösteriyor:

1. $\zeta(1+i\epsilon)$'nın kendisi bir Dirichlet serisidir: $\zeta(1+i\epsilon)=\sum_{n\ge1} n^{-1}e^{-i\epsilon\ln n}$. Dolayısıyla
 $$|\zeta(1+i\epsilon)|^2=\zeta(1+i\epsilon)\zeta(1-i\epsilon)=\sum_{n,m\ge1}\frac{1}{nm}e^{i\epsilon\ln(m/n)}.$$
 $m/n$'yi en sade haliyle $a/b$ ($\gcd(a,b)=1$) yazıp $n=bk,\,m=ak$ ile toplarsak:
 $$|\zeta(1+i\epsilon)|^2=\sum_{(a,b)=1}\frac{1}{ab}\,e^{i\epsilon\ln(a/b)}\sum_{k\ge1}\frac1{k^2}=\frac{\pi^2}{6}\sum_{(a,b)=1}\frac{1}{ab}\,e^{i\epsilon\ln(a/b)}.$$
 Yani $|\zeta(1+i\epsilon)|^2$ zaten, $\epsilon$'un Fourier eşleniğinde, **tam olarak her asal-olmayan koprim çift $(a,b)$'de bir çizgiye** (frekans $\ln(a/b)$'de) karşılık gelen ayrık bir toplamdır — sadece $\epsilon$-uzayında yazıldığında sürekli-görünen kapalı formu ($|\zeta(1+i\epsilon)|^2$) bunu gizliyor.
2. (7.5)'teki taşıyıcı $e^{2\pi i\bar d\epsilon}$ ile çarpılınca, her çizgi $\bar d$ (yani $L=2\pi\bar d=\ln(E/2\pi)$ mertebesinden) merkezli, konumu $\ln(a/b)$ kadar kaydırılmış bir "uydu" hâline gelir — kullanıcının sorduğu $\omega\approx L+\log(a/b)$ konumuyla **birebir örtüşüyor**.
3. Genlik/koherens: $\Phi^{(off)}(\epsilon)=\prod_p(1-(1-p^{i\epsilon})^2/(p-1)^2)$ ifadesi (7.1)–(7.6) türetiminde tam olarak Hardy–Littlewood tekil serisinin (6.8) resummasyonundan geliyor — yani $\mu(q)/\phi(q)$ ağırlıklı Ramanujan toplamı $c_q(r)$'nin $r$ üzerinden toplanmasıyla elde edilmiş bir üretici fonksiyon. Dolayısıyla $\Phi^{(off)}$'un $|\zeta(1+i\epsilon)|^2$'nin verdiği $(a,b)$-çizgileriyle çarpımı, her uydunun genliğine **aynı $\mu(\cdot)/\phi(\cdot)$ tipi katsayıyı** taşıyacak şekilde iç içe geçmiş oluyor.

Bu üç adım, ölçülen $\mu(a)/\phi(a)$-koherensli, $\log(a/b)$-konumlu uydu resmini **niteliksel olarak** üretiyor, ama BK metninde bu adımlar açıkça yürütülmüş/yazılmış değil — bu bir **türetim**, alıntı değil.

---

## 4. (c) Ölçülen uydu yapısıyla örtüşme derecesi

**Kısmen — yapısal olarak evet, ama BK metninde açık formda değil.**

- **Konum $\log(a/b)$**: BK'nın $R_2^{(off)}$ formülünün taşıyıcısı $|\zeta(1+i\epsilon)|^2$'den geldiği için (§3'teki türetim), ilke olarak birebir örtüşüyor — ama bu BK tarafından yazılı olarak Fourier-dual uzayda ifade edilmemiş.
- **Koherens $\mu(a)/\phi(a)$**: BK'nın $\Phi^{(off)}$'u DOĞRUDAN Hardy–Littlewood tekil serisinin ($\mu(q)/\phi(q))^2$ ağırlıklı Ramanujan toplamından (eş. 6.8) inşa edilmiş — bu, kullanıcının tanımladığı asal-kuvvet başına $\mu(a)/\phi(a)$ ağırlığıyla **aynı aritmetik iskelet**. Ama BK'nın nihai kapalı formu (7.6) $r$ üzerinden ZATEN toplanmış (resummed) bir Euler çarpımıdır — yani $\mu(a)/\phi(a)$ orada görünmez hâlde, örtük olarak gömülü.
- **Sınıf ağırlıkları $\cos(2\pi r b/a)$**: Ramanujan toplamı $c_a(rb)=\sum_{(p,a)=1}e^{2\pi i p r b/a}$ gerçel bir sayı olduğundan, $\cos$ terimlerinin toplamı olarak yazılabilir — bu standart bir özdeşlik (Ramanujan toplamının reelliği), BK metninde bu şekilde açık yazılmıyor, ama eş. (6.8)'in cebirsel bir sonucu.
- BK'nın kendi vurgusu (7.10) $\epsilon$-uzayında **sürekli** bir $\bar d^{-2}, \bar d^{-3}$ açılımı — yani "belirli $n$'lere karşılık gelen ayrık uydular" değil, tüm asalların/asal kuvvetlerin toplamının ürettiği **düzgün (smooth)** bir zarf/düzeltme eğrisi olarak sunuluyor. Ayrık uydu tablosu (her $a/b$'ye bir çizgi) BK'nın kendi anlatımının odağı değil; bu, formülün $\epsilon$'dan $\omega/\alpha$'ya Fourier-dönüşümü alındığında ortaya çıkan (ama BK'nın yazmadığı) bir sonuç.

**Sonuç**: Konum yapısı (log(a/b)) ve koherens iskeleti (μ/φ, Ramanujan toplamı) BK'nın (7.5)–(6.8) formüllerinde **matematiksel olarak mevcut/örtük**, ama "uydu" dili ve $\cos(2\pi r b/a)$ sınıf-ağırlıklı açık gösterim BK metninde **yazılı değil** — bu, okunan tam metinlerde doğrulanamayan, sadece türetilebilir bir eşdeğerliktir.

---

## 5. (d) Dirichlet L-fonksiyonları için BK/ratios genellemesi var mı?

**Kısmen — ama BK'nın kendi nesnesinin (tek bir L-fonksiyonun kendi sıfırlarının çift korelasyonu) birebir analoğu değil.**

Conrey–Snaith (math/0509480) §2.2 ve §3'te, gerçel kuadratik karakterli Dirichlet L-fonksiyonları ailesi $L(s,\chi_d)$ üzerinde çalışıyorlar (eş. 2.14–2.19), fakat bu bir **aile ortalaması** (χ_d, farklı $d$'ler üzerinden, sempletik simetri sınıfı) çerçevesinde **tek-seviye yoğunluğu (one-level density)** ve momentler için — BK'nın yaptığı gibi **tek, sabit bir L-fonksiyonun kendi sıfırlarının iki-nokta korelasyonu** değil. Aritmetik ağırlık burada $a(n)=\prod_{p|n} p/(p+1)$ (eş. 2.19) — μ(q)/φ(q) tipi Ramanujan ağırlığından **farklı bir aritmetik obje**. Tam metinde ne "Gauss sum" ne "conductor" ne de "epsilon factor" ifadeleri geçiyor (arandı, sıfır eşleşme) — yani karakterin Gauss-toplamı/iletken yapısına dayanan açık bir BK-tipi arithmetik-faktör genellemesi bu makalede **yok**.

Metnin girişinde (§1, satır 61, 109) BK açıkça "ilk bulan" olarak anılıyor ve Conrey-Snaith'in ratios-sanısı yönteminin BK'nın Hardy–Littlewood temelli ayrıntılı analizine göre daha kolay bir yol sunduğu vurgulanıyor ("The strength of our method is that it allows us to avoid such detailed considerations") — ama bu, zeta'nın KENDİ pair-correlation'ı için, Dirichlet L'ler için değil.

**Sonuç**: Okunan 5 kaynak içinde, sabit bir Dirichlet $L(s,\chi)$'nin kendi sıfırlarının BK-tipi (Gauss toplamı/iletken ağırlıklı) çift korelasyonu için doğrudan bir formül **yok / erişilemedi**. Var olan genelleme (Conrey-Snaith) ailelerin ortalamasına dair, farklı bir aritmetik-faktör ailesiyle.

---

## 6. HÜKÜM

Bogomolny–Keating'in aritmetik sanısı, ölçülen "uydu" yapısını **açıkça (explicit olarak, o dilde yazılı şekilde) öngörmüyor, ama formülünde örtük olarak var**: okunan tam metin (Bogomolny 2007, arXiv:0708.4223, eş. 6.1–7.10), $R_2(\epsilon)$'u her zaman enerji-farkı ($\epsilon$) uzayında, $|\zeta(1+i\epsilon)|^2$ ile Hardy–Littlewood tekil serisinin (eş. 6.8, $\sum_q(\mu(q)/\phi(q))^2 c_q(r)$ biçiminde — bu Ramanujan toplamının standart tanımı) resummasyonundan üretilen bir Euler-çarpımı $\Phi^{(off)}(\epsilon)$'un (eş. 7.6) çarpımı olarak kapalı formda veriyor; form faktörü $K(\tau)$ ya da Montgomery'nin $F(\alpha)$'sı tam metinde bir kez bile geçmiyor, dolayısıyla "α ≈ 1+log(n)/L civarında ayrık uydular" ifadesi BK'nın kendi diliyle hiçbir yerde yazılı değil. Ne var ki, $|\zeta(1+i\epsilon)|^2=\sum_{n,m}(nm)^{-1}e^{i\epsilon\ln(m/n)}$ açılımı (temel bir Dirichlet-serisi özdeşliği) ve $\Phi^{(off)}$'un Ramanujan-toplamı kökeni bir araya getirildiğinde, $R_2(\epsilon)$'un $\epsilon\to$Fourier-ikili dönüşümünün konumu $\log(a/b)$'de, ağırlığı $\mu(\cdot)/\phi(\cdot)$-tipi katsayılarla çizgiler üretmesi **matematiksel olarak zorunlu** görünüyor — yani yapı BK'nın formülünün DNA'sında var, sadece o dilde (form faktör/α uzayında) ifade edilmiş değil. Conrey–Snaith (2007) BK'nın sonucunu ratios-sanısıyla yeniden üretiyor (Teorem 4.1) ama yine $\epsilon$-uzayında kalıyor ve Dirichlet L-fonksiyonları için sunduğu genelleme (§2.2/§3) BK'nın kendi nesnesinin değil, aile-ortalamalı tek-seviye yoğunluğunun bir genellemesi — dolayısıyla (d) sorusuna net bir "yok" ile yakın bir cevap veriyor. PRL 1996 ve Nonlinearity I/II'nin tam metnine erişilemediği için, BK'nın ORİJİNAL makalelerinde form-faktör dilinde bir tartışma olup olmadığı kesin olarak **belirlenemedi**; sadece 2007 derlemesindeki yeniden-sunumun içeriği doğrulanabildi.

---

## 7. Sayısal karşılaştırma için tam formül (programcı için)

Aşağıdaki blok, kaynak 1'in tam metninden (eş. 3.5, 6.3–6.4, 6.8, 6.11–6.12, 7.5–7.6, 7.7–7.11) birebir derlenmiştir; bir programcının $R_2(\epsilon)$'u (ya da onun Fourier dönüşümünü, F(α) karşılığı için — bkz. §3'teki türetim, DOĞRULANMAMIŞ ek adım) sayısal hesaplaması için yeterlidir.

```
# Girdi: E (yükseklik), epsilon (enerji farkı, kırpılmamış) ya da varepsilon (kırpılmış)
# Asal listesi p = 2,3,5,7,11,... (yakınsaklık için birkaç bin asal yeterli)

d_bar(E) = (1/(2*pi)) * ln(E / (2*pi))                      # (3.5)

# --- Diagonal terim ---
# Phi_diag(eps) = exp( sum_p sum_{m=1}^inf [(1-m)/(m^2 p^m)] * exp(i*m*ln(p)*eps) + c.c. )
# R2_diag(eps) = -(1/(4*pi^2)) * d^2/d(eps)^2 [ ln( |zeta(1+i*eps)|^2 * Phi_diag(eps) ) ]
#   (ikinci türevi nümerik/simgesel al ya da eps->0 limitinde -1/(2*pi^2*eps^2) kullan)

# --- Off-diagonal terim (kapali form) ---
Phi_off(eps) = product_p ( 1 - (1 - p**(i*eps))**2 / (p-1)**2 )        # (7.6)
R2_off(eps)  = (1/(4*pi**2)) * abs(zeta(1+i*eps))**2 * exp(2j*pi*d_bar*eps) * Phi_off(eps)
R2_off(eps) += conj(R2_off(eps))                                       # "+ c.c."

R2(eps) = d_bar(E)**2 + R2_diag(eps) + R2_off(eps)                     # (7.7)

# --- Kirpilmis (unfolded), sonlu epsilon icin ---
R2_unfolded(varepsilon) = R2(varepsilon * d_bar(E)) / d_bar(E)**2      # (7.8)

# --- Ilk-mertebe acik asimptotik (E->inf, dogrudan kontrol icin) ---
beta  = gamma0**2 + 2*gamma1 + sum_p (ln(p)**2 / (p-1)**2)   # ~= 1.57314   (7.11)
delta = sum_p (ln(p)**3 / (p-1)**2)                          # ~= 2.3157   (7.11)
R2(varepsilon) ~= 1 - sin(pi*varepsilon)**2/(pi*varepsilon)**2 \
                    - (beta/(pi**2*d_bar**2)) * sin(pi*varepsilon)**2 \
                    - (delta/(2*pi**2*d_bar**3*varepsilon)) * sin(2*pi*varepsilon)   # (7.10)

# --- Hardy-Littlewood tekil serisi / Ramanujan acilimi (eps. 6.8, 6.11-6.12) ---
def ramanujan_sum(q, r):
    return sum(exp(2j*pi*p*r/q) for p in range(1, q) if gcd(p, q) == 1)

def alpha_HL(r, Q_max):
    return sum( (mobius(q)/euler_phi(q))**2 * ramanujan_sum(q, r) for q in range(1, Q_max) )
# r cift ise kapali form: alpha(r) = C2 * product_{p|r, p odd} (p-1)/(p-2),  C2 ~= 1.32032

# --- (TÜRETİLMİŞ, BK METNİNDE YAZILI DEĞİL) satellite/form-factor eşdeğeri ---
# |zeta(1+i*eps)|^2 = (pi^2/6) * sum_{(a,b)=1, coprime} (1/(a*b)) * exp(i*eps*ln(a/b))
# Bu yüzden R2_off(eps)'in Fourier donusumu (omega ekseninde), tasiyici e^{2i*pi*d_bar*eps}
# ile birlikte, omega = 2*pi*d_bar + ln(a/b)  [L = 2*pi*d_bar = ln(E/(2*pi)) ile omega = L + ln(a/b)]
# konumlarinda cizgiler / "uydular" verir; her cizginin genligi Phi_off'un ayni (a,b) bileseni
# ile carpilir (Phi_off, alpha_HL'nin r uzerinden resummasyonundan turetilmisti, eş. 7.1-7.6).
```

**Not**: `zeta(1+i*eps)` standart bir kütüphane fonksiyonuyla (örn. mpmath `zeta`) doğrudan hesaplanabilir; `Phi_diag`/`Phi_off` ürünleri birkaç bin asal ile (p < 10^4–10^6) pratik yakınsaklığa ulaşır (yakınsaklık oranı $\sim p^{-2}$'dir çünkü paydalar $(p-1)^2$).
