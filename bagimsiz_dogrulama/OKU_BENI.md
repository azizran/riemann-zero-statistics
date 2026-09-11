# Bağımsız Doğrulama Paketi — OKU BENİ

**Ne bu?** Uğur Sezen'in Riemann sıfır programındaki manşet iddiaların, **yazarın scriptleri kullanılmadan**, ham üçüncü-taraf veriden yeniden ölçüldüğü paket. 11 Eylül 2026, DSH ajanı tarafından, Uğur'un isteğiyle.

**Tek komutla hepsi:**
```bash
cd bagimsiz_dogrulama && python3 10_dogrula.py      # ~25 s → 21/21 PASS
```
```

---

## 1. Bağımsızlık beyanı

| | |
|---|---|
| **Kod** | Bu klasördeki tüm scriptler sıfırdan yazıldı. Yazarın `41_bigT_scan.py`, `74_kat_turetimi.py` vb. kodundan **satır alınmadı**; yalnızca fiziksel tanımlar (Riemann–Siegel formülü, `R_pred` özdeşliği) okundu. |
| **Motor doğrulaması** | Kendi `Z(t)` motoru mpmath'a karşı sınandı: `t ≥ 10⁵`'te maks hata ≤ 5×10⁻⁶. (İlk sürümümde θ(t)'de `−(t/2)ln π` terimi eksikti; hata kendi testimle yakalandı.) |
| **Veri** | Odlyzko'nun resmî tabloları. `zeros1` indirilip SHA-256'sı teyit edildi: `3436c916a7878261ac183fd7b9448c9a4736b8bbccf1356874a6ce1788541632` — repodaki `qm_riemann/veri_odlyzko/KAYNAK.txt` ile birebir. |
| **Kapsam** | 2.001.052 sıfır (t = 14…1,13×10⁶) + 10¹², 10²¹, 10²² tabloları; 200 surrogate; N = 5…22 CUE; 16,1 mertebe t aralığı. |

---

## 2. Sonuç özeti (2026-09-11 itibarıyla: sertifika **21/21 PASS**)

### ✅ Bağımsız olarak çoğalanlar

| İddia | Ölçülen |
|---|---|
| M² ortalaması Conrey–Ghosh eğimi | +%0.05 (makale −%0.05) |
| M² sabit terimi (HLPC) | −%0.14 (makale +%0.97) |
| r tablosu, 12 pencere | 12/12, sapma ≤0.004 |
| Gaussian surrogate null | 0.4938 ± 0.0059 (makale 0.494 ± 0.007) |
| Fazlalık | +51σ (makale +43σ — onların σ'sı 40 örneklemle ihtiyatlı) |
| N_eff anomalisi (Pearson kayması, Spearman sabitliği) | trend ve büyüklük tuttu |
| **Toplam kural** (on üç asalda u = 1) | 1.001–1.008 |
| — plasebo tabanı (aynı tayflı Gaussian) | 0.030 → **44–112σ ayrım** |
| Kanallar saf reel (kuadratür ~0) | faz 0.0–0.3°, kuadratür ~10⁻⁴ |
| `w` işaret değişimi + dönmeden 0→π atlama | var; kuadratür ≤%5 |
| **`v(τ)` kolapsı 10²¹ ve 10²²'de** | 10²¹/10¹² = 0.995 ± 0.147, 10²²/10¹² = 1.000 ± 0.142 → **16,1 mertebe** |
| Not 4 benek fazları 180° | 24/24 satır, sapma ≤0.3° |
| Not 3 (iv) doğrudan \|Z\|² kanalı `R_pred` | r = 0.9989 (21 asal, fold altı) |

### ⚠️ Konvansiyon teyidi bekleyenler

| Konu | Bende | Makalede |
|---|---|---|
| `√w–v` kısıtı eğimi | −1.85 (tüm asallar) / −1.26 (6 asal) | −0.884 |
| `w` tablosu seviyesi (p = 2…13) | 0.781…0.285 | 0.83…0.38 |
| `w` sıfır geçişi τ₀ (model-bağımsız) | **0.397** | 0.447 ± 0.005 |
| Tarak parlaklığı / DW oranı | 0.045–0.131 | 0.88–0.89 sabit |
| Benek genliği yasası (asal çizgiler) | %1–20 bant | %0.1–4 |

**Sorulacak dört şey:** (i) "own-gap kontrolü"nün tam tanımı, (ii) `w` tablosu τ-kolapsı mı pencere ortalaması mı, (iii) tarak parlaklığının normalizasyonu (`|Ĝ|` / `I(ω)` / ideal kafese göre), (iv) benek genliğinde DW/v normu.

### ✗ Şekil olarak çoğalmayan

- **"Çadır" (ayna) kanalı:** makale `R_tam − R_pred`'in foldda ~1.0'a çıkmasını bekliyor; benim ölçtüğüm kalıntı **~0.06–0.14** ve tepe yapmıyor. Doğrudan kanal kimliği (`R_pred`) ise r = 0.9989 ile tutuyor. → `RAPOR_acik_uyusmazliklar.md` §4-quater.

### ⛔ Veri bekleyen

- **Kalem 188** (ζ(τ) profil kimliği) ve **152–187 arku**: girdi blokları (`184/185/186` scratchpad) repoda ve bu diskte yok. → `RAPOR_veri_kurtarma_riski.md`.

---

## 3. Senin için eylem listesi

> **GÜNCELLEME 11 Eyl 2026:** Scratchpad **kurtarıldı** (277 MB zip → `qm_riemann/scratchpad/`, 26 görev klasörü).
> Pipeline bu makinede birebir çalışıyor (ζ = 0.3287 ∠179.96°). Kalem 188'in T1 testi **yapıldı** —
> sonuç: iptal aritmetik; ama H-F2 mührü surrogate null'a karşı zayıf (3/8 vekil de geçiyor).
> Ayrıntı: `NOT_188_teori.md` §5-bis ve `RAPOR_acik_uyusmazliklar.md` §4-quinquies.

1. ~~**Scratchpad'i kurtar (tek gerçek engel).**~~ ✅ yapıldı. Diğer makinede (içinde `/Users/ugursezen/Desktop/arin/deney` olan):
   `/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad` → zip'leyip getir, `qm_riemann/scratchpad/` altına aç. Köprü zaten kurulu (`09_yol_koprular.sh`), zincir çalışır hale gelir. 188 dahil her şey açılır.
2. **`QM` köprüsü için tek `sudo`** (yalnız gerektiğinde):
   ```bash
   sudo mkdir -p /Users/ugursezen/Desktop/arin
   sudo ln -sfn "/Users/ugur/Desktop/Deney" /Users/ugursezen/Desktop/arin/deney
   ```
   Bu, 189 script'in mutlak yollarını **hash'leri bozmadan** karşılar. (Toplu dosya düzenlemesi yapılmadı: 13 ön-kayıt scripti kendi sha256'sını donduruyor, düzenlemek mühürleri kırardı.)
3. **Dört konvansiyon sorusunu** (yukarıda) Claude'a sor — cevaplar gelince üç ⚠️ satırı ya ✅ ya da gerçek uyuşmazlık olur.

---

## 4. Teknik notlar (tekrar edenler için)

- **float64 tuzağı:** 10²¹ tablosunda `γ = taban + ofset`, `taban ≈ 1,44×10²⁰`. float64'te ULP ≈ 3×10⁴ olduğu için `taban + ofset` **ofseti yutar** (aralıklar 0 çıkar). Çözüm: faz sütunlarını **ofset** üzerinden kur — sabit faz kayması `cos/sin` çiftince yutulur, genlik değişmez; `L`'yi gerçek `t`'den hesapla. (`20_derin_zincir.py`)
- **`k` faktörü:** açık formülde asal-kuvvet satır genliği `1/(πk√q)`; norm `k√q` olmalı. `√q` ile asal-kuvvetler birimin %20–80 altında görünür (benim ilk hatam). (`15_k_faktoru_ve_tarak.py`)
- **Benek fazı:** `Ĝ(ω) = (1/n)Σ e^{iωt_n}`'de **mutlak** `t_n` kullanılmalı; ortalamayı çıkarmak `e^{−iωt̄}` bindirir ve 180° kilidini gizler (benim ikinci hatam). (`13_derin_pencere.py`)
- **Torak `Ψ(p)`:** `Ψ(p) = cos(2π(p²−p−1/16))/cos(2πp)` (p², p²/2 değil), `p = √(t/2π) − ⌊·⌋`. (`01_Z_motoru.py`)

---

## 5. Deponun kendi belgesinde bulunan erratum (yaması hazır, uygulanmadı)

`qm_riemann/KESIF_SEFERI_KUANTUM_KAOS_29AGU2026.md` satır 114:

```
- ... `c₀ = Σ_p (log p)⁴ Σ_{r≥1} (r−1)r²/p^r`, ayrıca ...
+ ... `c₀ = Σ_p (log p)²/(p−1)²`, ayrıca ...
```

Gerekçe: BBLM (arXiv math/0602270) tanımı `c_n = [(−1)ⁿ/(2n)!] Σ_p (log p)^{2(n+1)} Σ_r (r−1)r^{2n}/p^r`; `n=0`'da üs 2, `r^{2n}=1`. Doğru c₀ = **1.3855389** (repo'daki formül 33.81 verir). Λ = 1.573085 (makale 1.57314 ✓), C = Q/Λ = 1.47161 (makale 1.4720 ✓). Ayrıntı: `18_bblm_dogrulama.py`, `RAPOR_acik_uyusmazliklar.md` §4-ter.

---

## 6. Git'e gönderme (ÇALIŞTIRILMADI — onay bekliyor)

Depo `origin/main` ile senkron (`2c3ad16`). Bu paket **takipsiz**: `bagimsiz_dogrulama/`. Tek yerel değişiklik: `polygon_music/serve.js` (`/tunnel` rotası, korundu). Yedek: `_yedek_20260911/` (302 MB tar + çakışan dosyalar).

```bash
cd /Users/ugur/Desktop/Deney
git add bagimsiz_dogrulama qm_riemann/_yollar.py          # _yollar.py: yeni scriptler için taşınabilir QM/SCR
git status --short                                         # gözden geçir
git commit -m "Bağımsız doğrulama paketi: 21/21 PASS — Not 1/2/3/4 manşet sayıları, v(tau) zinciri 16.1 mertebe,
plasebo (toplam kural aritmetik), BBLM c0 erratumu, yol köprüleri; sertifika tek komut"
# git push origin main      # <-- SENİN ONAYINLA
```

---

## 7. Dosya envanteri

| Script | İş |
|---|---|
| `01_Z_motoru.py` | Z(t) motoru + mpmath doğrulaması |
| `02_olcum.py` | gap ve M_n, M²/asimptotik, r tablosu |
| `03_surrogate.py` | Theiler Gaussian surrogate null |
| `04_benek.py` | benek yasası ilk denemesi |
| `05_cue.py` | CUE eğrisi r(N) |
| `06_neff.py` | N_eff ters çevirme |
| `07_wv.py` | tam-taban regresyonu, w/v kanalları |
| `10_dogrula.py` | **sertifika — 21 kontrol, tek komut** |
| `11_placebo.py` | toplam kuralın plasebo testi |
| `12_w_konvansiyon.py` | `w` gap-koşullama varyantları |
| `13_derin_pencere.py` | 10¹²'de v kanalı |
| `14_kisit_fiti.py` | `√w–v` kısıtının yeniden fiti |
| `15_k_faktoru_ve_tarak.py` | `k` faktörü + Bragg tarağı |
| `16_isaret_degisimi.py` | `w` işaret geçişi (τ₀) |
| `18_bblm_dogrulama.py` | BBLM sabitleri + erratum |
| `19_ayna_cadir.py` | \|Z\|² modülasyonu ve çadır kanalı |
| `20_derin_zincir.py` | v(τ) zinciri 10⁵→10²² |
| `09_yol_koprular.sh` | mutlak yol köprüleri (hash'siz müdahale) |

**Raporlar:** `RAPOR_bagimsiz_cogaltma.md` (ana çoğaltma), `RAPOR_acik_uyusmazliklar.md` (hakem için açık noktalar), `RAPOR_veri_kurtarma_riski.md` (152–187 verisi), `NOT_188_teori.md` (ζ(τ) çerçevesi + T0–T4 testleri).
