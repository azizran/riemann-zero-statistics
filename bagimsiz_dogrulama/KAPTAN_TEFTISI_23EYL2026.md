# KAPTAN TEFTİŞİ — Bağımsız Doğrulama Paketi (23 Eylül 2026)

Bu klasördeki raporlar (11 Eylül, DeepSeek oturumu) birincil kayıttır ve
DEĞİŞTİRİLMEDİ. Aşağıdaki şerhler, raporları okuyacak herkes için teftiş
sonucudur. Ana defter kaydı: `qm_riemann/LITERATUR_TARAMASI_16AGU2026.md`
("ARA 11-23 Eyl" girdisi).

## Alınanlar (teyitli, programa katkı)

- **Bağımsız çoğaltma 21/21:** Not 1 sayıları (r tablosu 12/12, surrogate
  0.494, CG/HLPC), toplam kuralı (13 asal u≈1.00; plasebo 44-112σ), v(τ)
  kolapsı 10²²'ye dek (16.1 mertebe), benek fazları 24/24 @180°, tarak oranı
  0.8847±0.0050, benek yasası asallarda 0.987±0.054 (DW dahil). Yazarın
  betikleri kullanılmadan, SHA-teyitli ham Odlyzko verisinden.
- **BBLM c₀ erratumu doğru:** c₀ = Σ_p (log p)²/(p−1)² = 1.3855389
  (Σ_{r≥1}(r−1)p^{−r} = 1/(p−1)²); Λ = 1.573085 ✓. Değer doğruydu, yazım yanlıştı.
- **Metin önerileri (kabul):** "own-gap kontrolü"nün tanımı Not 2/3'e;
  N_eff'in BBLM'ninkiyle karışmaması için "korelasyon-kalibreli" etiketi;
  ζ alıntısında kesim (τ_c=0.86) yazılsın; iptal açısı daima modülle birlikte
  raporlansın.

## Düzeltilenler (teftişte yakalanan)

1. **"Çadır kanalı çoğalmadı (✗)" hükmü YANLIŞ — normalizasyon hatası.**
   §4-quater tablosunda ölçülen genlik 2|c|/DC q^{−1/2} çarpanını taşıyor,
   karşılaştırılan R_tam = 2(L−log q+2γ−1)/(L+2γ−1) taşımıyor. Raporun KENDİ
   sayıları √q ile düzeltilince: ölçülen/R_tam = 1.006, 1.011, 1.017, 1.027,
   1.062, 1.082, 1.058, 1.063 (q = 2…157) — 2(1−τ) ailesi %1-8 içinde TUTUYOR.
   Ayna kalıntısı (ölçülen−R_pred)·√q = 0.12, 0.20, 0.30, 0.43, 0.74 (fold
   altı, ≈2τ) → 1.07, 1.05 (fold üstü, ≈2(1−τ)): Not 3'ün "foldda tepe yapan
   çadır"ı raporun kendi verisiyle ÇOĞALIYOR (fold kenarı q=137'de yumuşama:
   0.74 vs 1.00 — kayıt).
2. **T1 (faz-rastgele vekiller) kurgusu geçersiz.** 185b/186b modeli
   (c^öz, c^kesik) nominal asal fazıyla cos(ω_q m) kurulur; φ_q-rastgele
   vekillere uygulanınca model denizin kendi çizgi fazlarıyla uyuşmaz, ζ
   O(1)-O(10) çıkar (gözlenen 2-7). Bu "iptal aritmetik" hükmünü sınamaz;
   o hükmün dayanağı 187'nin ikiz karşılaştırmasıdır (ζ 0.33 vs 0.077). Aynı
   sebeple "±15° ölçütünün ~%40 yanlış-pozitifi" sayısı anlamsızdır (ama
   "açıyı modülle raporla" önerisi yerinde).
3. **"188 kapalı form" türetim DEĞİL, ampirik uydurma.** 39 kayan pencere
   (genişlik 0.05, adım 0.01) ~%80 örtüşür → χ²/dof bağımsız-nokta
   varsayımıyla anlamsız; altı adaydan veri-sonrası seçim; τ_c=1.00'da
   tutmuyor. Değerli olan NİTEL bulgu: |ζ|(τ) ≈ taban + kesime-yakınlık
   artışı, minimum τ≈0.6-0.65 (T2/T3). 188'in asıl sorusu (kimlik + tarak
   köprüsü) açık kalıyor. "1/π" hipotezi paketin kendi hipoteziydi; düşmesi doğru.
4. **τ₀ = 0.397 vs 0.447 — açıklı.** 0.397 programın ESKİ dar-taban
   değeriyle (τ*≈0.40, 18 Ağu) örtüşür; tam-tabana geçişle 0.447'ye kaydı
   (görev 82-83). Aynı taban-duyarlılığı (Gap 1). Ayrıca dışbükey a(τ)'ya
   doğrusal fit sıfır geçişini sola çeker. Makalede taban-bağımlılığı açık yazılsın.

## Risk kaydı

152-187 ara verisi (~300 MB) yalnız Mac mini'de (`qm_riemann/scratchpad/`,
.gitignore'da); MacBook'taki geçici kopya sistemce silinmiş. 188, gerçek
sıfırlardan (repoda) + 164_insa.py'den yeniden kurulumla ve 187
sayılarının yeniden üretimi kapısıyla yürütülür.
