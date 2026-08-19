# arXiv GÖNDERİM REHBERİ — dörtleme (21 Ağustos 2026)

Bu dosya: adım-adım plan, endorsement mektubu taslağı, ve dört notun
gönderim-formu metadata'sı (kısaltılmış özetler dahil — form alanı
~1.920 karakter sınırlı; makale-içi özetler olduğu gibi kalır).

---

## A. YOL HARİTASI (sıra önemli)

1. **Hesap** (sen): https://arxiv.org/user/register — gmail ile olur.
   ORCID bağlamak önerilir (ücretsiz, orcid.org).
2. **Endorsement kodu al** (sen): giriş yaptıktan sonra "Start new
   submission" → arşiv `math` / kategori `math.NT` seç. Sistem
   endorsement gerekiyorsa söyler ve sana özel bir **endorsement kodu**
   verir (endorsement durumu e-posta alanına ve geçmişe göre değişir;
   kurumsal e-postası olmayan bağımsız araştırmacıda genelde gerekir).
3. **Endorser bul ve yaz** (sen; taslak mektup aşağıda, bölüm B).
   Endorser kodu https://arxiv.org/auth/endorse adresine girer; onay
   e-postası sana düşer. Süre: günler mertebesi olabilir — o yüzden
   ERKEN başlatıyoruz.
4. **Gönderim sırası** (endorsement gelince):
   - Not 1'i gönder → gönderim anında **arXiv kimliği** (arXiv:YYMM.NNNNN)
     atanır (duyurudan önce görünür).
   - Ben Not 2-3-4'ün `companion` bib girdilerine bu kimliği işlerim,
     yeniden derlerim.
   - Aynı oturumda Not 2 → kimlik → Not 3-4'e işle → Not 3 → Not 4.
     (Hepsi aynı gün gönderilebilir; duyuru bir sonraki iş günü,
     yeni gönderici için moderasyon 1-2 gün uzatabilir.)
5. **Teknik**: TeX kaynağı yüklenir (PDF değil) — dosya + 3-4 PNG figür.
   amsart + amsmath/hyperref/graphicx standarttır, arXiv'in TeX Live'ı
   derler. Lisans: varsayılan "arXiv non-exclusive license" yeterli.
6. **v2 planı**: duyuru sonrası çapraz-ID'ler zaten içeride olur;
   düzeltme gerekirse v2 ucuzdur.

Kategori önerisi: birincil **math.NT**; çapraz-liste (opsiyonel)
**math-ph** (Not 3-4 için uygun). İkisini de formda seçebilirsin.

---

## B. ENDORSEMENT — NE, KİM, NASIL

**Ne:** arXiv, yeni göndericilerin ilgili kategoriye ilk gönderiminden
önce o alanda yayın yapmış birinin "bu kişi bilimsel içerik gönderiyor"
onayını ister. İçerik hakemliği DEĞİLDİR — yalnızca "spam/saçma değil"
onayıdır; endorser makaleyi savunmuş sayılmaz.

**Kim endorser olabilir:** son yıllarda `math.NT` içinde yeterli sayıda
arXiv makalesi olan yazarlar (arXiv sayfalarında "Which authors of this
paper are endorsers?" bağlantısından tek tek kontrol edilebilir).

**Kime yazmalı:** doğal adaylar, notların zaten atıf verdiği ve konusu
en yakın olan araştırmacılar — ör. between-zero maxima literatürünün
yazarları (Hughes/Pearce-Crump çizgisi), sıfır-dağılımı ölçümleri
(Odlyzko), γ log p dağılımları (Ford). Nazik, kısa, PDF ekli tek
e-posta; 3-4 kişiye ayrı ayrı yazmak normaldir (ilk cevap veren yeter).

**Mektup taslağı (İngilizce, kişiselleştir):**

> Subject: arXiv endorsement request (math.NT) — numerical study of
> gap–amplitude statistics of Riemann zeros
>
> Dear Professor X,
>
> I am an independent researcher writing to ask whether you would be
> willing to endorse my first arXiv submission in math.NT (endorsement
> code: XXXXXX, entered at https://arxiv.org/auth/endorse).
>
> The submission is a short numerical note measuring the joint law of
> consecutive-zero gaps and between-zero maxima of Hardy's Z-function,
> validated against the Conrey–Ghosh and Hughes–Lugmayer–Pearce-Crump
> mean-value results [and building on your work on ...]. Three companion
> notes dissect the deviation from CUE into explicit-formula channels.
> I attach the PDF; code and data are public at
> https://github.com/azizran/deney.
>
> Endorsement only confirms the submission is appropriate for arXiv;
> it does not imply refereeing. Thank you for considering it.
>
> Best regards, Uğur Sezen (Istanbul)

---

## C. GÖNDERİM FORMU METADATA'SI

Ortak `Comments` kalıbı:
`N pages, M figures. Companion to [arXiv IDs once assigned]. Code and
data: https://github.com/azizran/deney (qm_riemann/).`

### Not 1 — arxiv_gap_amplitude.tex
- Title: The joint law of zero gaps and local maxima of Hardy's
  Z-function: a numerical study and a CUE effective-dimension anomaly
- Abstract: makale-içi özet OLDUĞU GİBİ (1.667 karakter ✓ sınıra uyar).
- Categories: math.NT (primary)

### Not 2 — arxiv_prime_wave_anatomy.tex (form özeti, ~1.600 kr)
> In a companion note we measured the joint law of consecutive-zero
> gaps and between-zero maxima of Hardy's Z-function, finding that its
> correlation matches CUE characteristic polynomials only at a shifted,
> metric-dependent effective dimension. Here we dissect the deviation
> by regressing log-maxima and log-gaps on the explicit-formula waves
> cos(t log p^k), with placebo controls throughout. This revision also
> identifies a sampling-grid systematic — the gap midpoints are
> displaced by the primes, biasing restricted-basis regressions
> (ground-truth-validated) — and quotes all channel numbers in a
> corrected complete-basis convention. We find: (i) a sum rule — under
> the complete basis all fifteen primes carry unit content within one
> percent, with no frequency drift, and the prime powers fall a few
> percent below; (ii) the content splits into a direct amplitude
> channel w and a zero-displacement channel v, each collapsing onto
> one-variable laws in tau = log p / log(t/2pi); (iii) stripping the
> measured waves leaves a core pair correlation r* ~ 0.977, insensitive
> to the basis and tighter than any CUE sampled; (iv) an exhaustion
> analysis yields a role-crossover at tau* ~ 0.40 and an upper bound on
> the intrinsic variance of the maxima sequence. The laws, fitted at
> t <= 1.6e6, are re-tested out of sample at the 10^12-th zero and the
> gap channel at the 10^21-st and 10^22-nd zeros, extending the
> validated range to eighteen orders of magnitude. At single-gap scale,
> two thirds or more of the apparent randomness of the maxima sequence
> is deterministic prime signal.
- Categories: math.NT (primary)

### Not 3 — arxiv_response_theory.tex (form özeti, ~1.700 kr)
> Two companion notes measured how the local statistics of Riemann
> zeros respond to the prime waves of the explicit formula: an
> amplitude channel w(tau) and a zero-displacement channel v(tau),
> tau = log p / L. Here we organize the measurements into a response
> theory — after identifying and correcting a sampling-grid systematic
> (the gap midpoints are displaced by the primes; a ground-truth test
> shows restricted bases bias w downward, while complete bases over all
> prime powers recover the truth to ~1%). Under the corrected
> convention: (i) the two channels obey a tight linear constraint
> sqrt(w) = 1.017 - 0.884 v (RMS 0.0075), extending to prime powers,
> surviving an out-of-sample test on unseen primes; intensity-unitarity
> fails by a factor 6, while the sqrt(w)-linear and w-linear forms are
> currently indistinguishable. (ii) w changes sign at
> tau_0 = 0.447 +/- 0.005 with response phases that jump from 0 to pi
> without rotating: consistent with a lossless medium. (iii) A thermal
> Bragg model reproduces the crossing: the ideal mirror sits at the
> Riemann-Siegel fold tau = 1/2; dressed with the measured jitter and
> plateau it places the zero at 0.443, within 1 sigma of measurement.
> (iv) The structure factor (1 - 2 tau) is derived as the pair count of
> the folded main sum; pointwise |Z|^2 regressions follow the
> shifted-second-moment family through the fold, so the mirror channel
> is a tent peaked at the fold — the approximate functional equation as
> an interference channel. (v) The stripped core correlation 0.977 is
> basis-insensitive; remaining previous-convention results are flagged
> for re-measurement.
- Categories: math.NT (primary), cross-list math-ph

### Not 4 — arxiv_warm_crystal.tex (form özeti, ~1.750 kr)
> Companion notes identified a sampling-grid systematic in per-gap
> measurements of Riemann-zero statistics: the gap midpoints are
> displaced by the primes. Here the grid becomes the object. Its
> structure factor G(omega) — the diffraction pattern of the zero
> lattice — has three components: arithmetic spots at log p^k, Bragg
> combs damped by the lattice temperature, and a dark field a factor
> ~25 below shot noise (hyperuniformity seen directly). The spots obey
> an absolute, parameter-free law from the explicit formula,
> |G(log q)| = Lambda(q) (L sqrt(q))^{-1} cos(pi tau) DW, verified to
> 0.1-4% across fourteen lines with phases consistent with pure real;
> the comb tracks the measured window temperature; and binned over the
> spots the pattern rebuilds the Montgomery ramp in the first
> line-bearing bin — Berry's decomposition observed directly. All Gram
> couplings between wave columns are values of G, and the regression
> dressing is consistent with the crystal's thermal diffuse scattering.
> The same displacement field yields a derived gap-channel sum rule —
> rigid onset 2 tau and scale 2/pi — filtered through a measured
> screening function D(tau), with Lambda-weighting confirmed for prime
> powers (k v(p^k) on the prime curve to 1-5%). The conditional gap
> response splits into two branches: to its own off-line fluctuations
> the gas responds like CUE for tau > 0.3 (within 3%, 0.3% at Nyquist),
> while below the first prime line it is effectively frozen; the
> coherent prime branch starts near unity and is screened, the branches
> crossing at tau* ~ 0.14. Two mechanisms for the measured prime-power
> weight deficit are refuted and the surviving constraints stated as
> open problems.
- Categories: math.NT (primary), cross-list math-ph

---

## D. GÖNDERİM GÜNÜ KONTROL LİSTESİ

- [ ] Endorsement onayı geldi (math.NT)
- [ ] Not 1 kaynak paketi: .tex (figürsüz) — yükle, derlemeyi önizle
- [ ] Not 1 kimliğini al → 2/3/4 bib'lerine işle (Claude), yeniden derle
- [ ] Not 2 paketi: .tex + fig_tau_en.png + fig_v_en.png
- [ ] Not 3 paketi: .tex + fig_law_en.png + fig_crossing_en.png + fig_tent_en.png
- [ ] Not 4 paketi: .tex + fig_diffraction_en.png + fig_geometry_en.png
      + fig_bridge_en.png + fig_twobranch_en.png
- [ ] Her formda: kısaltılmış özet (yukarıdan), Comments kalıbı,
      kategori(ler), lisans varsayılan
- [ ] Duyuru sonrası: GitHub README'ye arXiv rozetleri (Claude)
