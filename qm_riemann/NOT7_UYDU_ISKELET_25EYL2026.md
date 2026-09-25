# NOT 7 — TARAK UYDULARI: iskelet (25 Eylül 2026)

(Numara düzeltmesi: 8 Eylül planına göre (NOT5_ISKELET_08EYL2026.md) Not 5 = "Arithmetic
spectroscopy of L-function zero lattices" (96-118), Not 6 = "The locked phases of the Riemann
zero gas" (137-181). Bu not 188-198 arkını taşır ⇒ NOT 7. Not 5 ile temas: §5 (L-fonksiyonu
adaları, 195) Not 5'in ada verisini kullanır — çapraz atıf; Not 4 ile temas: G(ω) ve tarak.)

Durum: iskelet. Kaynak raporlar 188-198 (mühürlü). Her sayı yanında kaynağı var; LaTeX'e
dökülürken sayılar raporlardan yeniden okunacak. Atıflar [✓ doğrulandı] / [? doğrula].
Hedef uzunluk ~10-12 sayfa, 5 figür. Kategori math.NT (çapraz: math-ph).

## Başlık adayları

1. Arithmetic satellites of the Bragg comb of the Riemann zero lattice: residue classes,
   characters, and a blind test of the Bogomolny–Keating mirror law
2. Satellites of the Riemann zero comb and the Euler product of the Bogomolny–Keating pair
   correlation
3. Where the zero lattice hears the primes: comb satellites at log(a/b)

## Özet taslağı (EN)

Companion notes found that the per-gap response of Riemann zeros to the prime waves of the
explicit formula is cancelled beyond the solver window by a kernel that lives on the Bragg
comb of the zero-midpoint lattice. Here we study the satellites of that comb. (i) Their
positions form a height-universal grid Δω = ω − L_loc = log(a/b), L_loc = log(t/2π): a
pre-registered out-of-sample test at a new height places the peaks at 0.025, 0.725, 1.113,
1.800 against 0, log 2, log 3, log 6, and the profiles of two windows overlap with
correlation 0.98 on the Δω axis (0.002 on a τ-scaled axis). (ii) A stationary-phase
(Jacobi–Anger) analysis predicts that a satellite reads the residue class of the line
q' mod a; the naive rotation prediction dies, but a blind test confirms the refined law:
the class shares of +log 10 follow cos(2πr/10) within ±0.11. (iii) On Dirichlet
L-function islands the satellites whose prime divides the conductor vanish (21/21), the
cancellation keeps the zeta sign in both parities. (iv) The Euler product of the
Bogomolny–Keating pair correlation gives the satellite amplitudes in closed form,
c(r) = Π_p f_p(k_p) with f_p(k<0) = p^k, f_p(1) = p/(p−1)², f_p(k≥2) = 0 — a mirror law
(amplitude exactly 1/x at Δω = −log x) and a prime-power asymmetry. In a band never
computed before, a pre-registered blind profile test finds the mirror envelope at
1.023 ± 0.023 of the prediction (a line-density rival excluded at >12σ) and −log 8 on the
envelope. (v) Positive-exponent coefficients are transferred with a sub-unit factor κ_p
that relaxes toward 1 with height (1 − κ_p ∝ L^{−γ}, γ = 1.36 ± 0.21, three windows,
constant excluded at 6.4σ); since the ratios-conjecture pair correlation has no such
term, this is a property of the mixed zero–prime observable, whose derivation we leave
open. Nothing here bears on a proof of RH.

## Bölümler

### 1. Introduction
- Notlar 1-4 bağlamı (DOI 10.5281/zenodo.22942475): gap-genlik yasası, asal dalga
  anatomisi, tepki kuramı, sıcak kristal (orta-nokta örgüsünün yapı çarpanı G(ω):
  aritmetik noktalar log p^k, Bragg tarakları, karanlık alan).
- Yeni olan: tarak ω = L'nin UYDULARI ve aritmetik yapıları; ön-kayıt (KALEM) usulü, kör
  sınavlar ve ölümler dahil.
- Literatürdeki yeri: Landau–Gonek (tek-nokta, x^ρ toplamları) [? Landau 1911, Gonek 1993,
  DHPC 2601.18025 ✓]; asalların yapı çarpanı μ²/φ² (Zhang–Martelli–Torquato 1801.01541,
  Torquato–Zhang–de Courcy-Ireland 1802.10498) [✓ künye düzeltilmişti — yeniden doğrula];
  α>1 rejiminde BK konjektürel (Rodgers) [? doğrula]; Kanivets 2026 (Zenodo 20766728) —
  dipnot: "no arithmetic frequencies" iddiası büyük ihtimalle birim hatası (ln p yerine
  2π ln p / log(T/2π)).

### 2. The observable (tanım bölümü)
- Gap-başı özdeşlik ds_n = Σ_Q 2a_Q sin(ω_Q g_n/2) cos(ω_Q m_n), a_q = Λ(q)/(π√q log q).
- Pencere-ötesi iptal çekirdeği: C_{b,j}(q) = 2Σ_{n∈b} dds_{b,j}(n) e^{−iω_q m_n}/n_b,
  K_j = Σ_{q∈HAVUZ} C_j conj(mix_q)/Σ|mix_q|², κ = −Re K (188b/190b/197b).
- Blok-yerel eksen Δω = log q' − L_b; çizgi biçimi (blok-içi L dağılımı) ve neden 32 blok
  / eşit-ΔL (197/198). K_düz (dilimin kendi katkısı çıkarılmış).
- Kısa: iptalin büyüklüğü ζ_g = 0.333∠180° (düşük), 0.329∠180° (son) [190].

### 3. Positions: a height-universal grid (190; 188)
- 188: çekirdek bant-bağımsız τ'-tarağı, %77 rang-1; tepeler Bragg + uydular [188].
- 190 (kör, örneklem-dışı): tepeler 0.025, 0.725, 1.113, 1.800 (hedef 0, log2, log3, log6;
  ±0.05); τ'-sabit rakip (0.604/0.957/1.562) dışı; blok eğimleri 0; profiller Δω'da corr
  0.980, τ'-ölçekli eksende 0.002 [190].
- FİGÜR 1: iki pencerenin κ profili Δω ekseninde üst üste + uydu işaretleri.

### 4. Selection rules and residue classes (192-194)
- Türetim: seviye koşulu N(t_n) = n − ½, Jacobi–Anger e^{−2πiS} = Π_q Σ_k J_k(2πa_q)
  e^{ikω_q t}; durağan nokta L(t*) = ω' − log n; faz 2πq'b/a − 7π/4 − π/4 ≡ 0 (orta
  noktalar N(m_n) = n); koherens μ(a)/φ(a) (Ramanujan toplamı) [192, 193].
- Ölümler dürüstçe: dönme öngörüsü öldü (+log3 sınıfları −8.6° ± 3.9, 120° değil; H-192b)
  [192]; birleşme hipotezi (194) ve boş-pencere tabanı öldü [194].
- Kör geçiş: +log10'un dört sınıfı cos(2πr/10) işaretiyle, ±0.11 içinde, z = 11-52;
  +log7 yön doğru ama büyüklük 1.5-1.9× (KAYIT) [193].
- μ(a) = 0 konumları sönük ama küçük NEGATİF artık (L'den bağımsız, ~−0.02 s) [192, 197,
  198] — açık.
- FİGÜR 2: +log10 sınıf payları vs cos(2πr/10).

### 5. Dirichlet L-functions (195)
- L_χ = log(kt/2π), t* = 2πq'b/(ka), Gauss toplamı; k'yı bölen asallı uydular söner
  (yapı çarpanında 21/21; izinliler +9.8…+12.1σ); izinli uydular beş adada, iki paritede
  zeta ile aynı işaretli ve gerçel (∠180° ± 1.3°, z = −55…−82); χ₈ₑ/χ₈ₒ = 1.009 ± 0.024.
- Sapmalar: yasak artık izinlilerin %3-5'i, zıt işaretli (zeta μ=0 çukuruyla aynı düzey);
  √k/φ(k) ölçeği öldü (18-44σ) [195].
- FİGÜR 3: ada başına izinli/yasak uydular.
- NOT: 118b_chi3 / 118b_chi8o eksik sıfır çifti onarımı yayından önce (düşük öncelik ama
  bu bölüm için gerekli olabilir) [195].

### 6. The Bogomolny–Keating mirror law (196-197)
- BK/Conrey–Snaith çift korelasyonu (Bogomolny 2007 eş. 7.5-7.6 [✓ 0708.4223];
  Conrey–Snaith Thm 4.1 [✓ math/0509480]); Euler çarpımı ⇒ c(r) kapalı biçim (türetim:
  g_p'nin Fourier katsayıları; k≥2 iptali) [196, KALEM 197].
- Sonuçlar: ayna yasası (−log x zarfı 1/x; tam sayı ve buçuklu aileler), çeyrek aileler
  (0.25/x), asal-kuvveti asimetrisi.
- 196 veri-sonrası kıyas + şerh (8-blok çizgi biçimi karışık oranları bozuyordu).
- 197 kör sınav: bant [−2.12, −1.58], profil uyumu, iki-model öz-tutarlılık; R̄ = 1.023 ±
  0.023; rakip (çizgi yoğunluğu) >12σ dışı; ψ(−log8) = 0.993 ± 0.039; buçuklu aile 0.79 ±
  0.04 (bant içinde, ama 1'in altında — dürüstçe yazılır) [197].
- 32-blok ikincil okuma: +log(5/3) −0.44 → +0.42 vb. [197].
- FİGÜR 4: kör bant profili + R_x paneli (197 figürünün sadeleştirilmiş hâli).

### 7. Positive exponents and height (198)
- κ_p tablosu (üç pencere); γ̂ = 1.36 ± 0.21; sabit 6.4σ (çapasız 7.9σ) dışı; kör W_alt
  tek başına belirsiz (dürüstçe).
- Yorum (198 şerhi): CS/BK'da 1/L terimi yok ⇒ gözlenebilirin aktarım özelliği; ayna
  yanı birebir, pozitif üsler yükseklikle birebire. p^{−1/3} yükseklik tesadüfü.
- FİGÜR 5: κ_p vs L.

### 8. Open problems
1. Aktarım fonksiyonu T_p(L): çekirdek tanımından (gap-başı özdeşlik + stationary phase)
   türetim; negatif/pozitif üs asimetrisinin kaynağı.
2. μ=0 negatif artıkları (zeta ve L-adaları, ~%2-5): BK-ötesi terim mi?
3. Büyüklük bulmacaları: 193 +log7 fazlalığı (1.5-1.9×), 195 L-eşli ölçek ((k/φ(k))/ΠJ₀
   ipucu).
4. Buçuklu ailenin 0.8'i ile κ_2'nin ilişkisi (198'de ikisi de L ile 1'e gidiyor).

### Reproducibility / AI disclosure
- Yayın reposu + DOI; KALEM dosyaları ve sha256 damgaları; ölümler kayıtta.
- Claude (Anthropic) kapsamlı yardım; türetim teftişleri Sonnet, makine inşası Opus
  ajanları; yazar sorumluluğu.

## Kaynakça (doğrulama durumu)

- Bogomolny, E., Keating, J. P. (1996) — BK PRL / Nonlinearity 1995-96 [? tam künye; 196
  literatür raporunda künyeler doğrulandı — oradan al]
- Bogomolny, E. (2007) arXiv:0708.4223 [✓ eş. 7.5-7.6 okundu]
- Conrey, J. B., Snaith, N. C. (2007) Proc. LMS, math/0509480 [✓ Thm 4.1 okundu]
- Bogomolny, Bohigas, Leboeuf, Monastra (2006) math/0602270 [✓]
- Landau (1911); Gonek (1993) [? birincil künye]
- Durkan, Hughes, Pearce-Crump (2026) arXiv:2601.18025 [✓ tam metin okundu]
- Zhang, Martelli, Torquato (2018) arXiv:1801.01541 [✓ yazar hatası düzeltilmişti — yeniden
  bak]; Torquato, Zhang, de Courcy-Ireland (2018) arXiv:1802.10498 [✓]
- Rodgers (2012/2013) [? künye]
- Montgomery (1973); Odlyzko (1987 + tablolar); Keating–Snaith (2000) [✓ Not 1'de var]
- Ramanujan toplamı / Gauss toplamı standart kaynak (Apostol ya da Iwaniec–Kowalski) [?]
- Kanivets (2026) Zenodo 10.5281/zenodo.20766728 [✓ okundu] — dipnot
- Notlar 1-4 (Sezen 2026) — Zenodo DOI

## Açık kararlar (kullanıcıyla)

- Başlık seçimi.
- Ayrı not (Not 7) mi, yoksa Not 4'ün devamı mı (Not 4 "sıcak kristal" G(ω)'yı zaten kuruyor)? Not 5/6 sırası: yayın sırası Not 5 → 6 → 7 mi, yoksa Not 7 önce mi (en taze, kör sınavlı)?
- 118b onarımı notu bekletir mi?
- Figür dili: EN (notlar gibi).

## Sıradaki adımlar

1. Kaynakçayı birincil kaynaklardan doğrula (Sonnet, ucuz).
2. Figürler: raporların mevcut png'lerinden sadeleştirilmiş EN sürümleri (her biri bir
   betik, 1xx_fig_en.py deseninde).
3. LaTeX: bölüm bölüm; her sayı rapordan okunarak.
4. Türetim teftişi (Sonnet): §4 ve §6'daki cebir.
5. İç teftiş + ortak teftiş, sonra yayın reposu + Zenodo yeni sürüm + arXiv.
