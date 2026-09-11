# Reproducibility Riski — 152–187 Arku'nun Ara Verisi Kayıp

**Tarih:** 11 Eylül 2026 · **Durum:** doğrulandı, kanıtlı

---

## Bulgu

`qm_riemann/` içindeki 152–187 görev zinciri (Not 5/6 malzemesinin tamamı) **repodan yeniden çalıştırılamaz.** Zincirin ara çıktıları yalnızca **başka bir makinenin geçici scratchpad'inde** yaşıyor ve repoya hiç commit edilmemiş.

### Kanıtlar

| # | Bulgu | Ölçüm |
|---|---|---|
| 1 | Scriptler yabancı makine yoluna sabit kodlu | `qm_riemann/*.py` içinde **189/196** dosyada `/Users/ugursezen/...` veya `/private/tmp/claude-501/...` |
| 2 | 152–187 config'leri scratchpad'e bağımlı | `*_configs/*.py` içinde **191/210** |
| 3 | O kullanıcı bu makinede yok | `/Users/ugursezen` → **No such file or directory** |
| 4 | Temp'te o projenin girdisi yok | `/private/tmp/claude-501/` içinde `-Users-ugur`, `-Users-ugur-Desktop-3m`, `-Users-ugur-Desktop-belestepe` var; **`-Users-ugursezen-Desktop-arin` yok** |
| 5 | Hiç commit edilmemiş | `git rev-list --all --objects \| grep -i scratchpad\|K2_zeta\|ONKAYIT_18` → **boş** |
| 6 | 155–187 için tek bir `.npz` repoda yok | `git ls-tree -r origin/main \| grep -E '(15[5-9]\|1[6-8][0-9]).*\.npz'` → **0** |
| 7 | İkinci yerel klon da içermiyor | `~/Desktop/riemann` (576 MB, commit `72784cb`, 22 Ağu) → 152–187 config ve scratchpad **yok** |
| 8 | Tüm diskte kurtarılabilir kopya yok | `~/` altında `K2_zeta.json`, `ONKAYIT_18*.json`, `arin`, `K1_katman_defteri.json` araması → **sonuç yok** |

### Bağımlılık zinciri (config'lerin `SCR/'NNN'` referansları)

`155 → 156 → 157 → 158 → 161 → 164 → 165 → 166 → 167 → 169 → 172 → 173 → 174 → 175 → 176 → 177 → 178 → 180 → 181 → 182 → 183 → 184 → 185 → 186 → 187`

Her halka bir öncekinin `.npz`/`.json`'unu scratchpad'den okuyor (ör. 187, `SCR/155/z_Hkeskin.npy`, `SCR/184/K1_*.npz`, `SCR/185/OZ_*.npz`, `SCR/186/G1_proj_*.npz` dosyalarını istiyor). Zincirin kökü yoksa **hiçbir halka başlamaz** — ve dosyalar yerelde yok.

### Neden önemli

1. **Kalem 188 çoğaltılamaz.** ζ(τ) profili 187'nin çıktısına, o da 186→185→184→…→155 girdilerine dayanıyor. Bu veri olmadan `ζ = 0.329∠180°` bağımsız olarak yeniden ölçülemez — yalnız raporlardaki sayı okunabilir.
2. **Makalelerin Reproducibility bölümü fazla söz veriyor.** `arxiv_warm_crystal.tex` §Reproducibility: *"All measurements are scripts 79–95 … runnable from the stored datasets"* diyor ve GitHub'a işaret ediyor. 79–95 ve 101 serisi için veri kısmen var (`33/36/41/53/55/84/99/101c/101e/101f/102b/102c/117a/118b` npz'leri repoda); **152–187 için yok.**
3. **Hakem ilk denemede duvara çarpar.** Bir hakem `187_configs/187c_zeta_defteri.py`'yi açtığında `/Users/ugursezen/...` görüyor.

---

## Önerilen düzeltmeler (sırayla)

1. **Kurtarma (acil):** `/Users/ugursezen/Desktop/arin/deney` klonunun bulunduğu makine hâlâ duruyorsa, oradaki `scratchpad/15[2-9]`, `scratchpad/1[6-8][0-9]` klasörlerini (npz + json, tahminen onlarca MB) olduğu gibi arşivleyip repoya almak. Bu, 188 dahil her şeyi kurtarır.
2. **Taşınabilirlik:** 189 script'teki `QM`/`SCR` tanımlarını göreli hale getirmek:
   `QM = Path(__file__).resolve().parents[1]`, `SCR = Path(os.environ.get("QM_SCRATCH", QM/"scratchpad"))`.
   Mekanik bir refactor; istersen tek commit'te yapabilirim.
3. **Arşiv kuralı:** bundan sonra her görev çıktısını (`.npz`/`.json`) `qm_riemann/scratchpad/NNN/` altına yazıp commit etmek — `.gitignore`'a `scratchpad/*.tmp` gibi istisnalar koyup küçük artifact'ları takip etmek.
4. **Makale metni:** Reproducibility bölümlerini "şu dosyalar repoda, şunlar için şu girdi gerekir" diye kesinleştirmek; `Not 5/6` malzemesi için henüz yaşayan repo iddiası kurmamak.
5. **Sertifika:** her manşet sayıyı saklı veriden yeniden üretip karşılaştıran tek bir `dogrula.py` (ben yazabilirim) — böylece "runnable" iddiası test edilebilir olur.

---

## Bu turda yapılanlar (bağlam)

- Klon `origin/main`'e hizalandı (`2c3ad16`), takipsiz 1185 dosya yedeklendi (`_yedek_20260911/`, 302 MB).
- Bağımsız çoğaltma yapıldı ve raporlandı: `bagimsiz_dogrulama/RAPOR_bagimsiz_cogaltma.md`. Not 1 tam, Not 2 toplam kuralı tam, Not 3 saf-reel kanal tam, `√w–v` kısıtı kısmi, Not 4 kısmi.
- Kalem 188 için girdi verisi **yok** — bu rapor onun ön koşulunu belgeliyor.
