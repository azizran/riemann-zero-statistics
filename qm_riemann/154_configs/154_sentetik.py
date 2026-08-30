"""
154 — SENTETİK GAZDA AYNI R(τ) ÇIKARIMI
=======================================
Kullanım: 154_sentetik.py <keskin|A4> [std|ince]
  keskin : 152'nin kontrolü — SAF asal merdiveni, τ≤1.00'de KESKİN kesim,
           yani 0.53–0.79 aralığında soğurma YOK. z, 152'nin koştuğu
           Newton çözümünden okunur (z_keskin.npy).
  A4     : erfc-kesim τ_c=0.68, Δ=0.125 — 152'de gerçeğin σ_η² ve c₁'ini
           oturtan konfigürasyon. Newton burada YENİDEN koşulur.

KOPYA-KAYMASI DENETİMİ: merdiven + Newton kodu 152_gaz.py'nin main()'inden
BİREBİR kopyadır. Kopyanın sapmadığı, koşu sonunda 152'nin rapor ettiği
S3 üçlüsüyle (A4: σ_ds²=0.1128 σ_η²=0.0230 c₁=−0.00918; keskin: 0.1280/
0.0767/−0.03535) karşılaştırılarak SINANIR. R çıkarımı ise kopya değil:
gerçek pencerelerle AYNI 154_cekirdek.olc çağrısıdır.
"""
import sys, math, json, time
import numpy as np
from pathlib import Path
import multiprocessing as mp
import importlib

sys.path.insert(0, str(Path(__file__).resolve().parent))
CEK = importlib.import_module("154_cekirdek")

HERE = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR152 = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
              "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/152")
OUT = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/154")
TWO_PI = 2 * np.pi
NEWTON_IT = 20
NWORK = 6
NPT = 25000
S3_BEKLENEN = {"keskin": (0.1280, 0.0767, -0.03535),
               "A4": (0.1128, 0.0230, -0.00918)}


# ---- aşağıdaki 4 blok 152_gaz.py'den BİREBİR ----
def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8


def rvm_d(t):
    return np.log(t / TWO_PI) / TWO_PI


def S_ve_Sp(z, om, a, blok=800, npt=40000):
    S = np.zeros_like(z); Sp = np.zeros_like(z)
    for b0 in range(0, len(om), blok):
        w = om[b0:b0 + blok]; aa = a[b0:b0 + blok]
        for s0 in range(0, len(z), npt):
            sl = slice(s0, min(s0 + npt, len(z)))
            arg = np.outer(z[sl], w)
            S[sl] += -np.sin(arg) @ aa
            Sp[sl] += -np.cos(arg) @ (aa * w)
            del arg
    return S, Sp


_G = {}


def _init(om, a):
    _G["om"] = om; _G["a"] = a


def _work(zc):
    return S_ve_Sp(zc, _G["om"], _G["a"], npt=NPT)
# ---- birebir kopya sonu ----


def gaz_kur(ad, t_bas):
    """152_gaz.main()'in merdiven + (erfc) + Newton bölümü, birebir."""
    d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
    Z = np.sort(np.asarray(d["zeros"], dtype=float))
    zr = Z[len(Z) - 300000:]
    t0, t1 = float(zr[0]), float(zr[-1])
    L_hedef = float(np.log(0.5 * (t0 + t1) / TWO_PI))

    from sympy import primerange
    lim = int(np.exp(1.00 * L_hedef))
    lad = []
    for p in primerange(2, lim + 1):
        q, m = p, 1
        while q <= lim:
            lad.append((np.log(q), 1.0 / (np.pi * m * np.sqrt(q))))
            q *= p; m += 1
    om_l = np.array([w for w, _ in lad])
    a_l = np.array([a for _, a in lad])
    print(f"merdiven: {len(lad)} çizgi (τ≤1.00, L={L_hedef:.3f})", flush=True)

    if ad == "A4":
        tc, dl = 0.68, 0.125
        tau = om_l / L_hedef
        wgt = np.array([0.5 * math.erfc((x - tc) / dl) for x in tau])
        a_ham = a_l.copy()
        a_l = a_l * wgt
        print(f"  erfc-kesim τ_c={tc} Δ={dl}: Σa·w/Σa="
              f"{a_l.sum()/a_ham.sum():.4f}  (w>0.01 çizgi: "
              f"{int((wgt>0.01).sum())})", flush=True)

    n0 = int(np.ceil(rvm_N(t0))); n1 = n0 + 300000
    ns = np.arange(n0, n1, dtype=float)
    z = np.full_like(ns, 0.5 * (t0 + t1))
    for _ in range(30):
        F = rvm_N(z) - ns
        z = np.clip(z - F / rvm_d(z), 100.0, None)
        if np.max(np.abs(F)) < 1e-9:
            break
    print(f"  pürüzsüz çözüm: maks|F|={np.max(np.abs(rvm_N(z)-ns)):.2e}",
          flush=True)
    gbar_t = 1.0 / rvm_d(z.mean())
    ctx = mp.get_context("spawn")
    pool = ctx.Pool(NWORK, initializer=_init, initargs=(om_l, a_l))
    try:
        for it in range(NEWTON_IT):
            parts = np.array_split(z, NWORK)
            res = pool.map(_work, parts)
            S = np.concatenate([r[0] for r in res])
            Sp = np.concatenate([r[1] for r in res])
            F = rvm_N(z) + S - ns
            payda = np.maximum(rvm_d(z) + Sp, 0.3 * rvm_d(z))
            adim = np.clip(0.8 * F / payda, -1.0 * gbar_t, 1.0 * gbar_t)
            z = np.clip(z - adim, 100.0, None)
            aF = np.abs(F); mf = float(aF.max())
            if it % 4 == 0 or it == NEWTON_IT - 1 or mf < 1e-3:
                print(f"  Newton {it:2d}: maks|F|={mf:.3e} "
                      f"medyan|F|={np.median(aF):.3e} "
                      f"({time.time()-t_bas:.0f}s)", flush=True)
            if mf < 1e-3:
                break
    finally:
        pool.close(); pool.join()
    return np.sort(z)


def main(ad, bset):
    t_bas = time.time()
    print(f"=== SENTETİK {ad} / bant {bset} ===", flush=True)
    onb = OUT / f"z_{ad}.npy"
    if ad == "keskin":
        z = np.sort(np.load(SCR152 / "z_keskin.npy"))
        print(f"  152'nin keskin Newton çözümü okundu ({len(z)} nokta)",
              flush=True)
    elif onb.exists():
        z = np.load(onb)
        print(f"  önbellekten okundu: {onb}", flush=True)
    else:
        z = gaz_kur(ad, t_bas)
        np.save(onb, z)
        print(f"  z kaydedildi -> {onb}", flush=True)

    bant = CEK.BANTLAR_STD if bset == "std" else CEK.BANTLAR_INCE
    r = CEK.olc(z, f"sentetik-{ad}-{bset}", bantlar=bant)
    e = S3_BEKLENEN[ad]
    print(f"\nKOPYA DENETİMİ [{ad}] σ_ds² {r['s_ds']:.4f} vs {e[0]:.4f} | "
          f"σ_η² {r['s_eta']:.4f} vs {e[1]:.4f} | c₁ {r['c1']:+.5f} vs "
          f"{e[2]:+.5f}", flush=True)
    r["ad"] = ad; r["bantset"] = bset; r["s3_beklenen"] = e
    r["sure_s"] = time.time() - t_bas
    (OUT / f"R_sentetik_{ad}_{bset}.json").write_text(json.dumps(r, indent=1))
    print(f"\nbitti — {(time.time()-t_bas)/60:.1f} dk", flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "std")
