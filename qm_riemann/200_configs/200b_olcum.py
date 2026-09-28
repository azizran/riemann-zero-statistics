"""
200b_olcum.py — KALEM 200-B: ÖLÇÜM (b=2, yakın çift) → olay başına (pencere, n, m_n,
L_n, δ̃_n, M_n) → 200b_olaylar.npz
=================================================================================

*** BU BETİK GERÇEK VERİ ÜZERİNDE ÇALIŞTIRILMADI *** (görev kuralı: makine ajanı
yalnız İNŞA EDER ve --fake ile SENTETİK girdiyle uçtan uca test eder; gerçek ölçüm
ortak teftiş sonrası, eşikler donduktan ve ONKAYIT_200B push edildikten sonra başlar
— KALEM "Ölçüm tanımı" M9b: "ölçüm betiği yalnız ONKAYIT push'landıktan sonra çalışır").

Tasarım (KALEM "Gözlenebilir" + "Ölçüm tanımı" + Görev 1):
  - Olay: ardışık çift (γ_n, γ_{n+1}), HER İKİ indeks de pencere aralığında ve
    n+1 < pencere_bitiş (200b_ortak.events_from_zeros zaten bunu sağlıyor: zero_t
    tam pencere dilimidir, taşan çift üretilmez).
  - m_n=(γ_n+γ_{n+1})/2, L_n=log(m_n/2π), δ̃_n=(γ_{n+1}−γ_n)L_n/2π.
  - δ̃_n < 0.3 (ikincil ε=0.3'ü de kapsayacak şekilde EN GENİŞ eşik; 200b_analiz.py
    daha sonra ε∈{0.1,0.2,0.3} alt-kümelerini kesecek).
  - M_n = max_{[γ_n,γ_{n+1}]} |Z(t)| — 200a_motor Z, 33 iç nokta ızgara + 3-nokta
    parabolik rafinasyon (|Z| üzerinde; 200b_ortak.peak_grid_parabola).
  - Diske YALNIZ olay başına (pencere, n, m_n, L_n, δ̃_n, M_n) yazılır — hiçbir
    dağılım özeti (ortalama/varyans/k-istatistiği) burada HESAPLANMAZ (körlük
    200b_analiz.py'ye bırakılır; bu betik yalnız ÖLÇER ve KAYDEDER).

--fake bayrağı: motoru ve sıfır konumlarını TAMAMEN SENTETİK bir kaynakla değiştirir
(ayrı tohum, ayrı sahte t aralıkları) — şekil/IO testi için. Gerçek Motor B / gerçek
zeros6 dosyası --fake modunda HİÇ İÇE AKTARILMAZ.
"""
import argparse
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ortak_b = None  # gecikmeli import (aşağıda _load ile)


def _load(modname):
    import importlib.util
    spec = importlib.util.spec_from_file_location(modname, HERE / f"{modname}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


EPS_MAX_SAVE = 0.3  # KALEM: ikincil ε=0.3'ü de kapsayan en geniş eşik


# ---------------------------------------------------------------------------
# SENTETİK sahte pencereler / sahte motor (yalnız --fake): gerçek veriyle SIFIR ilişki.
# ---------------------------------------------------------------------------
def load_fake_windows(n_fake_zeros, seed=999111):
    rng = np.random.default_rng(seed)  # SENTETİK test tohumu (gerçek değil)
    bounds = [(1000.0, 3000.0), (3000.0, 6000.0), (6000.0, 10000.0)]
    wins = {}
    for i, (t_a, t_b) in enumerate(bounds):
        name = f"FAKE{i+1}"
        zt = np.sort(rng.uniform(t_a, t_b, n_fake_zeros))
        wins[name] = {"idx": (0, len(zt)), "t_a": t_a, "t_b": t_b, "zero_t": zt}
    return wins


def fake_Z_func(t, rng):
    """Sahte motor: gerçek Riemann-Siegel ile hiçbir ilgisi yok, ama sıfır-benzeri
    bir şekli var (deterministik salınım + küçük gürültü) — grid+parabol testinin
    en azından 'düzgün bir tepe bulma' davranışını sınaması için."""
    t = np.asarray(t, dtype=np.float64)
    return np.sin(0.37 * t) * (1.0 + 0.05 * rng.standard_normal(t.shape)) + 0.3 * np.cos(0.071 * t)


# ---------------------------------------------------------------------------
# Ana ölçüm
# ---------------------------------------------------------------------------
def measure(windows, Z_func, n_grid, chunk, out_path, save_max_events=None):
    """windows: dict(name -> {idx:(i0,i1), t_a, t_b, zero_t}). Z_func: t(1B)->Z(1B).
    save_max_events: --fake testinde/güvenlik için toplam olay sayısını sınırlamak
    isteyen çağıran (None ise sınırsız — gerçek ölçümde kullanılacak mod)."""
    ortak_b = _load("200b_ortak")

    all_win, all_n, all_m, all_L, all_dt, all_M = [], [], [], [], [], []

    for wname, win in windows.items():
        i0, _ = win["idx"]
        ev = ortak_b.events_from_zeros(win["zero_t"], i0, eps_max=EPS_MAX_SAVE)
        n_ev = len(ev["n"])
        if save_max_events is not None and n_ev > save_max_events:
            # yalnız --fake/test modunda: rastgele alt-örnek (gerçek ölçümde KULLANILMAZ)
            rng = np.random.default_rng(12321)
            keep = rng.choice(n_ev, size=save_max_events, replace=False)
            for k in ev:
                ev[k] = ev[k][keep]
            n_ev = save_max_events

        M = ortak_b.peak_grid_parabola(Z_func, ev["ga"], ev["gb"], n_grid=n_grid, chunk=chunk)

        all_win.append(np.full(n_ev, wname))
        all_n.append(ev["n"])
        all_m.append(ev["m"])
        all_L.append(ev["L"])
        all_dt.append(ev["delta_tilde"])
        all_M.append(M)
        print(f"  {wname}: δ̃<{EPS_MAX_SAVE} olay sayısı={n_ev} -> M hesaplandı (n_grid={n_grid})")

    out = {
        "window": np.concatenate(all_win) if all_win else np.array([], dtype="<U8"),
        "n": np.concatenate(all_n) if all_n else np.array([], dtype=np.int64),
        "m": np.concatenate(all_m) if all_m else np.array([]),
        "L": np.concatenate(all_L) if all_L else np.array([]),
        "delta_tilde": np.concatenate(all_dt) if all_dt else np.array([]),
        "M": np.concatenate(all_M) if all_M else np.array([]),
    }
    n_total = len(out["n"])
    # Diske yalnız olay-başına 6 alan yazılır (M9 ruhu: hiçbir dağılım özeti burada YOK)
    assert all(len(v) == n_total for v in out.values()), "alan uzunlukları eşleşmiyor"

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(out_path, **out)
    print(f"Kaydedildi: {out_path}  (toplam olay={n_total}, hiçbir dağılım özeti hesaplanmadı)")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fake", action="store_true",
                     help="SENTETİK stand-in motor + sentetik sahte pencereler (gerçek veri YOK)")
    ap.add_argument("--fake-nzero", type=int, default=6000,
                     help="--fake modunda sentetik 'sıfır' sayısı / pencere")
    ap.add_argument("--fake-max-events", type=int, default=None,
                     help="--fake modunda olay sayısını sınırlamak için (yalnız test)")
    ap.add_argument("--n-grid", type=int, default=129)  # M2b: 33-nokta 1e-6 kapısını geçemedi, 129 geçti (26 Eyl)
    ap.add_argument("--chunk", type=int, default=2_000_000)
    ap.add_argument("--seed", type=int, default=200200)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    t0 = time.time()
    if args.fake:
        ortak_b = _load("200b_ortak")
        windows = load_fake_windows(args.fake_nzero, seed=args.seed)
        rng_fake = np.random.default_rng(args.seed + 1)
        Z_func = lambda t: fake_Z_func(t, rng_fake)
        out_path = args.out or (HERE.parent / "scratchpad" / "k200b_makine" / "test_olaylar_fake.npz")
        print("*** --fake modu: SENTETİK motor + SENTETİK pencereler, gerçek veri kullanılmıyor ***")
        measure(windows, Z_func, args.n_grid, args.chunk, out_path,
                save_max_events=args.fake_max_events)
    else:
        ortak_b = _load("200b_ortak")
        motor = _load("200a_motor")
        real_windows = ortak_b.load_real_windows()
        zeros = np.load(HERE.parent / "128_odl_zeros6_2e6_zeros.npz")["zeros"]
        windows = {}
        for w, info in real_windows.items():
            i0, i1 = info["idx"]
            windows[w] = {"idx": (i0, i1), "t_a": info["t_a"], "t_b": info["t_b"],
                           "zero_t": zeros[i0:i1]}
        out_path = args.out or (HERE / "200b_olaylar.npz")
        print("!!! GERÇEK MOD: bu, gerçek Z değerlerini ve M_n'yi hesaplar. Görev kuralı "
              "gereği bu ajan tarafından ÇALIŞTIRILMAMALIDIR (yalnız inşa edilip --fake ile "
              "test edilmiştir; gerçek koşu ONKAYIT_200B push'landıktan SONRA yapılır).")
        measure(windows, motor.Z_vec, args.n_grid, args.chunk, out_path)

    print(f"Toplam süre: {time.time()-t0:.1f} s")


if __name__ == "__main__":
    main()
