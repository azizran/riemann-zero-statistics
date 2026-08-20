# Not 4 Etki Listesi — 102 Taperli Yeniden-Denetim Seferi (20 Ağustos 2026)

Kaynak: 102a/b/c scriptleri (Opus tayfası, denetim Fable). Bu liste
`arxiv_warm_crystal.tex` revizyon oturumu içindir; TEX HENÜZ DEĞİŞMEDİ.

## MERKEZİ YENİ KİMLİK (102a-A5)

    Ĝ_orta(ω) = (1/n) Σ_n e^{iωz_n} · e^{iωg_n/2}

Orta noktalar sıfır örgüsünü *dalgalanan yarım-gap kaymasıyla* örnekler;
⟨e^{iωg/2}⟩ Debye-Waller çarpanı, dalgalanması difüz saçılmadır.
Aynı ω'da SIFIRLAR 8-10 kademe daha karanlık: sıfır tayfı çözülmüş
asal-kuvvet çizgi ormanı + GERÇEK boşluk; "karanlık alan" ve "TDS"
orta-nokta ızgarasının malıdır. Örnekleme-fazı eğrisi (yeni gözlemlenebilir):
t(c) = z_n + ḡ/2 + c(g_n−ḡ)/2 → I: c=0: 2.4e-14 | c=0.5: 2.3e-4 |
c=1 (orta): 4.2e-4 TEPE | c=2: 2.4e-14 (yine örgü).

## DEĞİŞMESİ GEREKENLER (satır numaraları mevcut tex'e göre)

1. **33-35 (özet):** "dark field ~25 below shot noise (hyperuniformity
   seen directly)" → taperli çizgi-dışı düzey banda göre 15-490× shot-
   noise altı; alan ORTA-NOKTA örneklemesinin; sıfırlarda veri-sınırına
   dek karanlık.
2. **46-48 (özet):** "dressing lives in the thermal diffuse scattering"
   → "...lives in the off-line continuum of the midpoint grid — the
   non-linear half-gap sampling sidebands of the arithmetic lines,
   whose level a noiseless explicit-formula lattice reproduces to 6%."
3. **42-45 (özet, Berry):** güçlendirilebilir: ilk çizgili binde rampa,
   mutlak benek yasasıyla %2 içinde parametresiz; çizgi-üstü pay %99.9.
4. **145-152 (§diff, Dark field parag.):** ekle: (a) ω<1 tapersiz taban
   kusursuz örgüden ayırt edilemez (1.20×); (b) Hann: 48×/490×/15×
   (hepsi/ω<1/ω>3); (c) doğru ifade = bağımsız-gap vekilinin 40-1150×
   altı; (d) sıfırlarda ≤7e-16 (veri-sınırlı). "Hyperuniformity ...
   seen directly as darkness" cümlesi orta-nokta ızgarası için yeniden.
5. **154-162:** fig_diffraction'a taperli + sıfır-ızgarası eğrileri
   (102a_karanlik_alan.png malzeme).
6. **238-247 (§ruler, TDS parag.):** korunacak: çizgi-dışılık, 4-18×,
   mertebe. Değişecek: mesafeler 0.03-0.44 → çözünürlüğün 7-199 katı
   (tam 7327'lik p^k listesi); 1B sıcak-kristal kestirimi kaldır
   (rastlantısal uyum); kimlik → örnekleme ızgarasının difüz saçılması;
   lab örgüsü medyanda %6. Tek cümle: "the bending of the regressions
   is arithmetic in origin and lives in the diffuse scattering of the
   sampling grid."
7. **249-260 (Montgomery/Berry):** "18 dark" → aritmetiksiz RvM kontrol
   aynı sayıyı verir (17.2×); yenisi ≥6.5e3× (orta) / ≥2.1e12× (sıfır),
   alet-sınırlı. "0.102 vs 0.089" AYAKTA, güçlendir (taper-değişmez
   0.1013, %99.9 çizgi-üstü, mutlak yasa 0.0997/0.1086 → 1.016/1.038).
   "10-20 times below" → 924× (orta, taperli) / 2.5e5× (sıfır); 88'in
   sayısı sinc kuyruklarıyla ~90× şişikti. "At larger α falls behind
   ramp" → sıfırlarda rampa ±%25 izlenir; açık orta-nokta cos(πτ)+DW
   sönümüdür.
8. **356-361 (§branches):** yalnız "factor 18" ibaresi → "at least
   three orders on the midpoint grid, at least twelve on the zeros";
   donma iddiası etkilenmez.
9. **420-426 (§Synthesis):** "thermal diffuse scattering" → "diffuse
   scattering of the sampling grid — arithmetic in origin"; dark field
   → "two to three orders below shot noise once window leakage removed".
10. **439-441 (Açık problem 5):** yeni sayılarla.
11. **446 (Reproducibility):** scripts listesine 101h-j, 102a-c ekle.

**Yeni eklenebilir kazançlar:** (i) örnekleme-fazı eğrisi; (ii) mutlak
benek yasasının bin düzeyinde rampayı %2 ile öngörmesi; (iii) "sıfır
tayfı atomik, orta-nokta tayfı atomik+sürekli" ayrımı (Not 5 köprüsü).

## AYAKTA KALANLAR (dokunma)

Gram=Ĝ kimliği (5e-8); mutlak benek yasası (çizgi-üstü; 102c bin
düzeyinde bağımsız doğruladı); faz kilidi 180°; πτcot(πτ) ve k₀; tarak/
Edgeworth/0.885; kolektif transfer; ilk çizgili binde rampa inşası
(GÜÇLENDİ); v=(2/π)sin(πτ), D(τ), Λ-ağırlığı, 1/k; iki-dal + donma
(101i tek-frekans kontrastıyla bağımsız: ζ 1.93× vs β 0.97×); kuvvet
açığı + iki ölü mekanizma; 87'nin "taşıyıcılar çizgi-dışı, 4-18×"
olgusu (yalnız kimlik ve mesafe sayıları değişir).

## DÜRÜST KAYITLAR (102 raporundan)

- Sıfır tabanı 7.3e-16 VERİ-SINIRLI (bisection kuantizasyonu σ=3.3e-8'e
  1.3-1.5×) → "≥10¹²×" sayıları üst sınır.
- Lab örgüsü uyumu DÜZEY uyumu (medyan 1.06; nokta-nokta 0.48-4.7×) —
  mekanizma+mertebe kanıtlı, frekans-frekans kimliği değil; 120k orman
  kontrastında 0.41× ters işaret (ince yapıyı üretemiyor).
- I(ω) üstel dağılımlı; medyan/ortalama ln2 → mutlak sayılar ~1.44×
  oynayabilir. DW konvansiyonu (σ_t sıfırlardan mı ortalardan mı)
  yüksek-α'yı ~1.3× oynatır. 102a çizgi maskesi 86'dan sıkı.
