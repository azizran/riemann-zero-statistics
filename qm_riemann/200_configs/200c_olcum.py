"""
200c_olcum.py — KALEM 200-C: OLCUM (b0,b1,b2) -> 200c_olcum.npz
=================================================================================

*** BU BETIK GERCEK VERI UZERINDE CALISTIRILMADI *** (gorev kurali: makine ajani
yalniz INSA EDER ve --fake ile SENTETIK girdiyle uctan uca test eder; gercek
olcum ortak teftis sonrasi, esikler donduktan ve ONKAYIT_200C push edildikten
sonra baslar).

Tasarim (KALEM Gorev 6):
  - b0: pencere basina 10^6 duzgun rastgele t (tohum 300), x=log|Z(t)|.
  - b1: pencere basina TUM 1.5e6 sifir EGER projekte edilen calisma suresi
    <=30 dk, aksi halde 10^6 rastgele sifir (tohum 301). x=log(|Z'(gamma_n)|*2pi/L_n).
    KARAR (gercek motor hiz olcumunden -- bkz asagida "B1_DECISION" ve
    scratchpad/k200c_makine/ hiz notlari, M3c'nin GERCEK 200-olay x (129+513)
    izgara kosusundan cikarilan per-nokta maliyet C1=0.0195ms C2=0.081ms
    C3=0.29ms C4=1.9ms; 1e6/1.5e6 nokta icin projeksiyon:
      C1: b0~19.5s b1~29s (30 dk'nin COK altinda) -> TUMU
      C2: b0~81s b1~122s -> TUMU
      C3: b0~290s(4.8dk) b1~435s(7.3dk) -> TUMU
      C4: b0~1900s(31.7dk) b1~2850s(47.5dk) -> b1 30 dk'yi ASIYOR -> ALT-ORNEK (1e6, tohum 301)
    Bu tahminler MOTOR HIZI uzerine (Z/Z' hesap maliyeti), hicbir Z/Z'
    DAGILIM OZETI icermez -- korluk ihlali degil.).
  - b2: pencere basina TUM olaylar (delta_tilde < 0.3): n (dosya-ici sifir
    indeksi), m_n (ofset), L_n, delta_tilde_n, M_n (max|Z| tepe, 129-nokta
    izgara+parabol -- M2b/M3c GEREGI 33 DEGIL 129 kullanilir, KALEM 200-B
    deneyiminin tasidigi ders).
  - Ortak 128 esit-t-uzunluklu blok/pencere (t sinirlarinda; b0 orneklerinin VE
    b1 sifirlarinin VE b2 olaylarinin (m_n ile) AYNI blok sinirlarina gore
    kutulanmasi -- 200c_analiz.py'nin ortak-blok capraz-basamak kovaryansi
    (H-200C-1) icin GEREKLI).
  - Diske b0/b1 icin YALNIZ blok basina guc toplamlari (n,Sx..Sx4) yazilir;
    b2 icin olay basina diziler (n,m,L,delta_tilde,M) yazilir -- hicbir ornek/
    olay uzerinden DAGILIM OZETI burada HESAPLANMAZ (korluk 200c_analiz.py'ye
    birakilir).

--fake bayragi: motoru VE sifir konumlarini TAMAMEN SENTETIK, gercek veriyle
hicbir iliskisi olmayan bir kaynakla degistirir (ayri tohum, ayri sahte t
araliklari) -- sekil/IO testi icin. Gercek Motor C / gercek zeros dosyalari
--fake modunda HIC ICE AKTARILMAZ.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import importlib.util


def _load(modname):
    if modname in sys.modules:
        return sys.modules[modname]
    spec = importlib.util.spec_from_file_location(modname, HERE / f"{modname}.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[modname] = m
    spec.loader.exec_module(m)
    return m


NBLOCKS = 128
SEED_B0 = 300
SEED_B1 = 301
EPS_B2_SAVE = 0.3
N_GRID_B2 = 129  # KALEM 200-B dersi: 33 nokta 1e-6 kapisini gecemedi, 129 gecti

# KARAR (bkz modul dokstring ust-kismindaki hiz projeksiyonu): b1 icin
# pencere basina "all" (1.5e6 sifir) ya da "subsample" (1e6, tohum 301).
B1_DECISION = {"C1": "all", "C2": "all", "C3": "all", "C4": "subsample"}


# ---------------------------------------------------------------------------
# SENTETIK sahte pencereler / sahte motor (yalniz --fake): gercek veriyle SIFIR
# iliski.
# ---------------------------------------------------------------------------
class FakeWindow:
    __slots__ = ("name", "T0", "off", "n_read", "t_a_off", "t_b_off")


def load_fake_windows(n_fake_zeros, seed=999222):
    rng = np.random.default_rng(seed)
    bounds = [(1000.0, 3000.0), (3000.0, 6000.0), (6000.0, 10000.0), (10000.0, 15000.0)]
    wins = []
    for i, (a, b) in enumerate(bounds):
        w = FakeWindow()
        w.name = f"FAKE{i+1}"
        w.T0 = 1_000_000  # sahte tamsayi taban
        w.off = np.sort(rng.uniform(a, b, n_fake_zeros))
        w.n_read = len(w.off)
        w.t_a_off, w.t_b_off = float(w.off[0]), float(w.off[-1])
        wins.append(w)
    return wins


def fake_Z_func(t_off, rng):
    t_off = np.asarray(t_off, dtype=np.float64)
    return np.sin(0.31 * t_off) * (1.0 + 0.05 * rng.standard_normal(t_off.shape)) + 0.2 * np.cos(0.053 * t_off)


def fake_Zprime_func(t_off, rng):
    t_off = np.asarray(t_off, dtype=np.float64)
    return 3.0 * np.cos(0.31 * t_off) * (1.0 + 0.05 * rng.standard_normal(t_off.shape)) + 1.5


# ---------------------------------------------------------------------------
# Blok yardimcilari (ortak, --fake ve gercek modda ayni)
# ---------------------------------------------------------------------------
def block_index(t_off, t_a, t_b, nblocks):
    idx = ((t_off - t_a) / (t_b - t_a) * nblocks).astype(np.int64)
    return np.clip(idx, 0, nblocks - 1)


def accumulate_into(acc, blk_idx, x, nblocks):
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


TWO_PI = 2 * np.pi


def L_of(T0, t_off):
    return np.log((T0 + np.asarray(t_off, dtype=np.float64)) / TWO_PI)


# ---------------------------------------------------------------------------
# Ana olcum
# ---------------------------------------------------------------------------
def measure(windows, get_Z_and_Zprime, get_events_peak, n0, nblocks, seed_b0, seed_b1,
            chunk, out_path, fake=False, b1_decision=None):
    """get_Z_and_Zprime(win, t_off) -> (Z, Zp). get_events_peak(win) -> dict with
    n,ga,gb,m,L,delta_tilde (KONUM, --fake'de de saglanir)."""
    n_win = len(windows)
    acc = np.zeros((n_win, 2, nblocks, 5), dtype=np.float64)  # b0,b1 guc toplamlari
    edges = np.zeros((n_win, nblocks + 1), dtype=np.float64)

    rng_b0 = np.random.default_rng(seed_b0)
    rng_b1 = np.random.default_rng(seed_b1)

    all_win2, all_n2, all_m2, all_L2, all_dt2, all_M2 = [], [], [], [], [], []

    for wi, win in enumerate(windows):
        t_a, t_b = win.t_a_off, win.t_b_off
        edges[wi] = np.linspace(t_a, t_b, nblocks + 1)

        # ---- b0: n0 uniform t ----
        remaining = n0
        while remaining > 0:
            take = min(chunk, remaining)
            t_chunk = rng_b0.uniform(t_a, t_b, take)
            Z_chunk, _ = get_Z_and_Zprime(win, t_chunk)
            x_chunk = np.log(np.abs(Z_chunk))
            blk = block_index(t_chunk, t_a, t_b, nblocks)
            accumulate_into(acc[wi, 0], blk, x_chunk, nblocks)
            del t_chunk, Z_chunk, x_chunk, blk
            remaining -= take

        # ---- b1: tum sifirlar ya da alt-ornek ----
        decision = (b1_decision or {}).get(win.name, "all")
        if decision == "all":
            zt_all = win.off
        else:
            n_sub = min(1_000_000, win.n_read)
            idx_sub = rng_b1.choice(win.n_read, size=n_sub, replace=False)
            zt_all = np.sort(win.off[idx_sub])

        nz = len(zt_all)
        for s in range(0, nz, chunk):
            zt_chunk = zt_all[s:s + chunk]
            _, Zp_chunk = get_Z_and_Zprime(win, zt_chunk)
            L_chunk = L_of(win.T0, zt_chunk)
            x_chunk = np.log(np.abs(Zp_chunk) * TWO_PI / L_chunk)
            blk = block_index(zt_chunk, t_a, t_b, nblocks)
            accumulate_into(acc[wi, 1], blk, x_chunk, nblocks)
            del zt_chunk, Zp_chunk, L_chunk, x_chunk, blk

        # ---- b2: TUM olaylar (delta_tilde<0.3), M_n = 129-izgara+parabol tepe ----
        ev = get_events_peak(win)
        n_ev = len(ev["n"])
        M = ortak_b_peak(win, ev, get_Z_and_Zprime) if not fake else ev["M_fake"]

        all_win2.append(np.full(n_ev, win.name))
        all_n2.append(ev["n"])
        all_m2.append(ev["m"])
        all_L2.append(ev["L"])
        all_dt2.append(ev["delta_tilde"])
        all_M2.append(M)

        print(f"  {win.name}: b0 n={n0}  b1 n={nz} ({decision})  b2 olay={n_ev} (eps<{EPS_B2_SAVE})")

    b2_out = {
        "window": np.concatenate(all_win2) if all_win2 else np.array([], dtype="<U8"),
        "n": np.concatenate(all_n2) if all_n2 else np.array([], dtype=np.int64),
        "m": np.concatenate(all_m2) if all_m2 else np.array([]),
        "L": np.concatenate(all_L2) if all_L2 else np.array([]),
        "delta_tilde": np.concatenate(all_dt2) if all_dt2 else np.array([]),
        "M": np.concatenate(all_M2) if all_M2 else np.array([]),
    }

    # ---------------- M9: korluk bariyeri assert'i ----------------
    assert acc.shape == (n_win, 2, nblocks, 5), f"beklenmeyen sekil: {acc.shape}"
    n_events_total = len(b2_out["n"])
    for k, v in b2_out.items():
        assert len(v) == n_events_total, f"b2 alan uzunluklari eslesmiyor: {k}"
    print(f"M9 kontrolu GECTI: blok_guc.shape={acc.shape}, b2 olay-basina 6 alan (toplam {n_events_total} olay) "
          f"-- hicbir ornek-basina dagilim ozeti diske yazilmiyor.")

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    to_save = {"blok_guc": acc, "t_kenar": edges,
               "pencere_adlari": np.array([w.name for w in windows]),
               "b1_decision": np.array([(b1_decision or {}).get(w.name, "all") for w in windows]),
               **b2_out}
    np.savez(out_path, **to_save)
    print(f"Kaydedildi: {out_path}  (blok_guc {acc.nbytes/1024:.1f} KB, b2 {n_events_total} olay)")
    return acc, edges, b2_out


def ortak_b_peak(win, ev, get_Z_and_Zprime):
    ortak_b = _load("200b_ortak")
    Z_func = lambda t_off: get_Z_and_Zprime(win, t_off)[0]
    return ortak_b.peak_grid_parabola(Z_func, ev["ga"], ev["gb"], n_grid=N_GRID_B2)


# ---------------------------------------------------------------------------
# GERCEK mod baglamalari
# ---------------------------------------------------------------------------
def real_get_Z_and_Zprime_factory():
    ortak = _load("200c_ortak")
    motor = _load("200c_motor")
    tables = {}

    def get(win, t_off):
        if win.name not in tables:
            tables[win.name] = motor.build_anchor_table(win.T0, win.t_a_off, win.t_b_off, n_workers=8)
        return motor.Z_and_Zprime(t_off, tables[win.name])
    return get


def real_get_events_peak(win):
    ortak = _load("200c_ortak")
    ev = ortak.events_from_window(win, eps_max=EPS_B2_SAVE)
    return ev


def fake_get_events_peak_factory(rng):
    def get(win):
        ga = win.off[:-1]
        gb = win.off[1:]
        n_idx = np.arange(len(win.off) - 1, dtype=np.int64)
        m = 0.5 * (ga + gb)
        L = L_of(win.T0, m)
        delta = gb - ga
        delta_tilde = delta * L / TWO_PI
        mask = delta_tilde < EPS_B2_SAVE
        M_fake = np.abs(fake_Z_func(m[mask], rng)) + 0.5
        return {"n": n_idx[mask], "ga": ga[mask], "gb": gb[mask], "m": m[mask], "L": L[mask],
                "delta_tilde": delta_tilde[mask], "M_fake": M_fake}
    return get


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fake", action="store_true")
    ap.add_argument("--fake-nzero", type=int, default=8000)
    ap.add_argument("--n0", type=int, default=1_000_000)
    ap.add_argument("--nblocks", type=int, default=NBLOCKS)
    ap.add_argument("--seed-b0", type=int, default=SEED_B0)
    ap.add_argument("--seed-b1", type=int, default=SEED_B1)
    ap.add_argument("--chunk", type=int, default=300_000)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    t0 = time.time()
    if args.fake:
        windows = load_fake_windows(args.fake_nzero)
        rng_fake = np.random.default_rng(555222)

        def fake_get_ZZp(win, t_off):
            return fake_Z_func(t_off, rng_fake), fake_Zprime_func(t_off, rng_fake)

        out_path = args.out or (HERE.parent / "scratchpad" / "k200c_makine" / "test_olcum_fake.npz")
        print("*** --fake modu: SENTETIK motor + SENTETIK pencereler, gercek veri kullanilmiyor ***")
        measure(windows, fake_get_ZZp, fake_get_events_peak_factory(rng_fake),
                n0=min(args.n0, 20000), nblocks=args.nblocks, seed_b0=args.seed_b0,
                seed_b1=args.seed_b1, chunk=args.chunk, out_path=out_path, fake=True,
                b1_decision={w.name: "all" if i < 2 else "subsample" for i, w in enumerate(windows)})
    else:
        ortak = _load("200c_ortak")
        windows = [ortak.load_window(w) for w in ortak.WINDOWS]
        out_path = args.out or (HERE / "200c_olcum.npz")
        print("!!! GERCEK MOD: bu, gercek Z/Z' degerlerini hesaplar ve pencere basina blok "
              "guc toplamlarini / b2 olay dizilerini yazar. Gorev kurali geregi bu ajan "
              "tarafindan CALISTIRILMAMALIDIR (yalniz insa edilip --fake ile test edilmistir).")
        measure(windows, real_get_Z_and_Zprime_factory(), real_get_events_peak,
                n0=args.n0, nblocks=args.nblocks, seed_b0=args.seed_b0, seed_b1=args.seed_b1,
                chunk=args.chunk, out_path=out_path, fake=False, b1_decision=B1_DECISION)

    print(f"Toplam sure: {time.time()-t0:.1f} s")


if __name__ == "__main__":
    main()
