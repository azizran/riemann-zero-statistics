# İki Mac arası senkron (MacBook + Mac mini)

**İlke:** Kod, defter ve raporlar GitHub üzerinden senkronlanır. Büyük ara veri
(`qm_riemann/scratchpad/`, git dışı) TAŞINMAZ, gerektiğinde gerçek sıfırlardan
yeniden üretilir (~40 dk, 187'nin mühürlü sayılarıyla kapılı).

## Her oturum

| ne zaman | komut |
|---|---|
| oturum başı | `./araclar/senkron.sh` — GitHub'dan çeker, veri köprüsünü kurar, durumu yazar |
| oturum sonu | `./araclar/senkron.sh --gonder` — yerel commit'leri GitHub'a gönderir |

Komutlar repo kökünden (`deney/`) çalıştırılır. Claude ile çalışırken bunu Claude
da yapar (her mühürde push standart).

## Mac mini'de ilk kurulum (bir kez)

1. `cd /Users/ugur/Desktop/Deney && git pull` (araçlar gelir)
2. `./araclar/senkron.sh` — "eski yol yok" derse yazdığı iki `sudo` satırını bir kez
   çalıştır (eski betiklerin mutlak yolu `/Users/ugursezen/...` için köprü; mühürlü
   ön-kayıt betikleri düzenlenmediği için sha'lar korunur).
3. Veri eksikse: `./araclar/yeniden_kur.sh` (numpy/scipy/mpmath'li bir Python ister;
   repo içinde `.venv` yoksa `python3` kullanılır).

## Kurallar

- Aynı anda iki makinede aynı işe dokunma. Bir makineden çıkmadan `--gonder`.
- Betik "AYRIŞMA" derse otomatik birleştirme yapılmaz; Claude'a sor.
- Force-push asla.
