"""
200a_olcum.py — KALEM 200-A: ÖLÇÜM (b=0, b=1) → yalnız blok güç toplamları
=================================================================================

*** BU BETİK GERÇEK VERİ ÜZERİNDE ÇALIŞTIRILMADI *** (görev kuralı: makine ajanı
yalnız İNŞA EDER ve SENTETİK girdiyle uçtan uca test eder; gerçek ölçüm ortak teftiş
sonrası, eşikler donduktan ve ONKAYIT_200A.json push edildikten sonra başlar).

Tasarım (KALEM "Ölçüm tanımı" + Kapı M9):
  - rung b=0: pencere başına n0 nokta, t ~ Uniform[t_a, t_b] (varsayılan tohum 200),
    x = log|Z(t)|.
  - rung b=1: penceredeki BÜTÜN sıfırlar, x = log(|Z'(gamma_n)| * 2*pi / L(gamma_n)).
  - Ortak t-blokları: 128 eşit-t-uzunluklu blok / pencere (sınırlar t'de; b=0 örnekleri
    ve b=1 sıfırları AYNI blok sınırlarına göre kutulanır).
  - Diske YALNIZ blok başına güç toplamları yazılır: (n, Sum x, Sum x^2, Sum x^3, Sum x^4).
  - Hiçbir sıfır-başına / örnek-başına dizi diske YAZILMAZ; bellekteki büyük diziler her
    parça (chunk) sonunda atılır (M9). Betik sonunda np.savez'e verilen her dizinin
    boyutu (n_pencere * 2 * n_blok * 5) mertebesini AŞAMAYACAĞI assert edilir.

--fake bayrağı: motoru ve sıfır konumlarını TAMAMEN SENTETİK, gerçek veriyle hiçbir
ilişkisi olmayan bir kaynakla değiştirir (ayrı tohum, ayrı sahte t aralıkları) — şekil/IO
testi için. Gerçek Motor B / gerçek zeros6 dosyası --fake modunda HİÇ İÇE AKTARILMAZ.
"""
import argparse
import importlib.util
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
N_BLOCKS_DEFAULT = 128


def _load(modname):
    spec = importlib.util.spec_from_file_location(modname, HERE / f"{modname}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ---------------------------------------------------------------------------
# Blok biriktirici
# ---------------------------------------------------------------------------
def block_index(t, t_a, t_b, nblocks):
    idx = ((t - t_a) / (t_b - t_a) * nblocks).astype(np.int64)
    return np.clip(idx, 0, nblocks - 1)


def accumulate_into(acc, blk_idx, x, nblocks):
    """acc: (nblocks,5) float64, yerinde güncellenir. x, blk_idx bu çağrıdan sonra
    çağıran tarafından atılabilir (burada saklanmaz)."""
    counts = np.bincount(blk_idx, minlength=nblocks).astype(np.float64)
    s1 = np.bincount(blk_idx, weights=x, minlength=nblocks)
    s2 = np.bincount(blk_idx, weights=x ** 2, minlength=nblocks)
    s3 = np.bincount(blk_idx, weights=x ** 3, minlength=nblocks)
    s4 = np.bincount(blk_idx, weights=x ** 4, minlength=nblocks)
    acc[:, 0] += counts
    acc[:, 1] += s1
    acc[:, 2] += s2
    acc[:, 3] += s3
    acc[:, 4] += s4


# ---------------------------------------------------------------------------
# Gerçek motor / gerçek pencereler (yalnız --fake OLMADIĞINDA içe aktarılır)
# ---------------------------------------------------------------------------
def load_real_windows():
    ortak = _load("200a_ortak")
    tahmin = ortak.load_tahmin()
    zeros = np.load(HERE.parent / "128_odl_zeros6_2e6_zeros.npz")["zeros"]
    wins = []
    for w in ortak.WINDOWS:
        idx0, idx1 = tahmin[w]["idx"]
        t_a, t_b = tahmin[w]["t"]
        wins.append({"name": w, "t_a": float(t_a), "t_b": float(t_b),
                     "zero_t": zeros[idx0:idx1]})
    return wins


def real_engine_Z(t, chunk=None):
    motor = _load("200a_motor")
    return motor.Z_vec(t, chunk=chunk)


def real_engine_Zprime(t, chunk=None):
    motor = _load("200a_motor")
    return motor.Zprime_vec(t, chunk=chunk)


# ---------------------------------------------------------------------------
# SENTETİK sahte motor / sahte pencereler (yalnız --fake): gerçek veriyle SIFIR
# ilişki. Ayrı, açıkça "sahte" bir tohumla üretilir.
# ---------------------------------------------------------------------------
def load_fake_windows(n_fake_zeros):
    rng = np.random.default_rng(999999)  # SENTETİK test tohumu (gerçek 200/2001/2002 değil)
    wins = []
    bounds = [(1000.0, 2000.0), (2000.0, 3200.0), (3200.0, 5000.0)]
    for i, (t_a, t_b) in enumerate(bounds):
        zt = np.sort(rng.uniform(t_a, t_b, n_fake_zeros))
        wins.append({"name": f"FAKE{i+1}", "t_a": t_a, "t_b": t_b, "zero_t": zt})
    return wins


def fake_engine_Z(t, rng):
    # t'den bağımsız, saf sentetik gürültü — gerçek RS formülüyle hiçbir ilgisi yok.
    return rng.standard_normal(t.shape) * 1.3 + 0.05


def fake_engine_Zprime(t, rng):
    return rng.standard_normal(t.shape) * 4.0 + 6.0


# ---------------------------------------------------------------------------
# Ana ölçüm döngüsü
# ---------------------------------------------------------------------------
def measure(windows, n0, nblocks, seed, chunk, fake, out_path):
    n_win = len(windows)
    acc = np.zeros((n_win, 2, nblocks, 5), dtype=np.float64)
    edges = np.zeros((n_win, nblocks + 1), dtype=np.float64)

    rng_t = np.random.default_rng(seed)
    rng_fake = np.random.default_rng(seed + 555000) if fake else None

    for wi, win in enumerate(windows):
        t_a, t_b = win["t_a"], win["t_b"]
        edges[wi] = np.linspace(t_a, t_b, nblocks + 1)

        # ---- rung 0: n0 uniform t, parçalar hâlinde (M9: bellek sınırlı, disk YOK) ----
        remaining = n0
        while remaining > 0:
            take = min(chunk, remaining)
            t_chunk = rng_t.uniform(t_a, t_b, take)
            if fake:
                Z_chunk = fake_engine_Z(t_chunk, rng_fake)
            else:
                Z_chunk = real_engine_Z(t_chunk)
            x_chunk = np.log(np.abs(Z_chunk))
            blk = block_index(t_chunk, t_a, t_b, nblocks)
            accumulate_into(acc[wi, 0], blk, x_chunk, nblocks)
            del t_chunk, Z_chunk, x_chunk, blk
            remaining -= take

        # ---- rung 1: penceredeki TÜM sıfırlar, parçalar hâlinde ----
        zero_t = win["zero_t"]
        nz = len(zero_t)
        for s in range(0, nz, chunk):
            zt_chunk = zero_t[s:s + chunk]
            if fake:
                Zp_chunk = fake_engine_Zprime(zt_chunk, rng_fake)
            else:
                Zp_chunk = real_engine_Zprime(zt_chunk)
            L_chunk = np.log(zt_chunk / TWO_PI)
            x_chunk = np.log(np.abs(Zp_chunk) * TWO_PI / L_chunk)
            blk = block_index(zt_chunk, t_a, t_b, nblocks)
            accumulate_into(acc[wi, 1], blk, x_chunk, nblocks)
            del zt_chunk, Zp_chunk, L_chunk, x_chunk, blk

        print(f"  {win['name']}: rung0 n={n0}, rung1 n={nz} -> bloklandı (nblocks={nblocks})")

    # ---------------- M9: körlük bariyeri assert'i ----------------
    assert acc.shape == (n_win, 2, nblocks, 5), f"beklenmeyen şekil: {acc.shape}"
    to_save = {"blok_guc": acc, "t_kenar": edges,
               "pencere_adlari": np.array([w["name"] for w in windows])}
    for key, arr in to_save.items():
        if key == "pencere_adlari":
            continue
        assert arr.size <= n_win * 2 * (nblocks + 1) * 5, (
            f"M9 IHLALI: '{key}' dizisi beklenenden büyük (olası örnek-başına sızıntı): "
            f"{arr.size} elemanlı")
    print(f"M9 kontrolü GEÇTİ: blok_guc.shape={acc.shape}, örnek başına hiçbir dizi diske yazılmıyor.")

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(out_path, **to_save)
    print(f"Kaydedildi: {out_path}  (blok_guc {acc.nbytes/1024:.1f} KB)")
    return acc, edges


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fake", action="store_true",
                     help="SENTETİK stand-in motor + sentetik sahte pencereler (gerçek veri YOK)")
    ap.add_argument("--n0", type=int, default=1_000_000, help="rung-0 örnek sayısı / pencere")
    ap.add_argument("--fake-nzero", type=int, default=4000,
                     help="--fake modunda sentetik 'sıfır' sayısı / pencere")
    ap.add_argument("--nblocks", type=int, default=N_BLOCKS_DEFAULT)
    ap.add_argument("--seed", type=int, default=200)
    ap.add_argument("--chunk", type=int, default=200_000)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    if args.fake:
        windows = load_fake_windows(args.fake_nzero)
        out_path = args.out or (HERE.parent / "scratchpad" / "k200a_makine" / "test_bloklar_fake.npz")
        print("*** --fake modu: SENTETİK motor + SENTETİK pencereler, gerçek veri kullanılmıyor ***")
    else:
        windows = load_real_windows()
        out_path = args.out or (HERE / "200a_bloklar.npz")
        print("!!! GERÇEK MOD: bu, gerçek Z/Z' değerlerini hesaplar ve pencere başına blok güç "
              "toplamlarını yazar. Görev kuralı gereği bu ajan tarafından ÇALIŞTIRILMAMALIDIR "
              "(yalnız inşa edilip --fake ile test edilmiştir).")

    measure(windows, n0=args.n0, nblocks=args.nblocks, seed=args.seed,
            chunk=args.chunk, fake=args.fake, out_path=out_path)


if __name__ == "__main__":
    main()
