"""
152 — ERFC-KESİMLİ SENTETİK GAZ TARAMASI  (151b'nin parametrik kopyası)
=======================================================================
151b_sentetik_gaz2.py'nin ölçüm zinciri BİREBİR korunmuştur. Tek fark:
merdiven genliklerine bir pencere (H-A) ya da Newton çözümünden sonra
konumlara bağımsız titreşim (H-B) uygulanabilmesi.

NEDEN TEK DOSYA, 6 KOPYA DEĞİL: altı ayrı kopyada zincirin bir satırının
kazara ayrışması (kopya kayması) bulguyu sessizce bozar. Burada zincir
TEK kod yolundadır; her konfigürasyon KONFIG sözlüğünde tek bir satırdır
— "tek değişiklik" böylece denetlenebilir kalır.

KONFİGÜRASYONLAR
  keskin   : değişiklik yok — 151b kontrolü (20 iterasyon).
             H-B'nin Newton tabanı da budur; z diske yazılır.
  A1/A2/A3 : H-A, erfc-kesim (Berry–Keating 5.24 ailesi)
             a_q → a_q · ½erfc((τ_q − τ_c)/Δ),  τ_q = ω_q / L_hedef
             A1=(0.75,0.10)  A2=(0.60,0.15)  A3=(0.90,0.05)
  B1/B2    : H-B, inkoherent seyreltme — keskin merdiven AYNEN, Newton'dan
             SONRA z += N(0, σ_j·ḡ_t),  ḡ_t = 2π/L,  rng(7)
             B1: σ_j=0.05   B2: σ_j=0.10
             (Newton'u tekrarlamaz; 'keskin' koşusunun z'sini okur —
              merdiven özdeş olduğu için hesap da özdeştir.)

erfc: stdlib math.erfc (tam; 1−erf'in büyük x'teki iptalinden kaçınır).
S(z) değerlendirmesi 6 süreçli havuzla paralel — seri sonuca göre maks
fark 1e-15 (yalnız toplama sırası), yani sayısal olarak özdeş.

  A4       : H-A, S3'u gercege oturtan kesim (0.68/0.125) —
             taramadan sonra, A1/A2 kusatmasindan turetildi.
  C1       : A4 + titresim 0.10 birlikte ('ikisi de kismi mi?').

Kullanım:  python 152_gaz.py <keskin|A1|A2|A3|A4|B1|B2|C1>
"""

import sys
import time
import math
import json
import numpy as np
from pathlib import Path
import multiprocessing as mp

HERE = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCRATCH = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
               "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/152")
TWO_PI = 2 * np.pi
NEWTON_IT = 20          # görev notu: 20 iterasyon yeter
NWORK = 6               # S(z) havuzu
NPT = 25000             # süreç başına nokta bloğu (bellek tavanı)

# konfig: (tip, parametreler)
KONFIG = {
    "keskin": ("kontrol", {}),
    "A1":     ("erfc", {"tau_c": 0.75, "delta": 0.10}),
    "A2":     ("erfc", {"tau_c": 0.60, "delta": 0.15}),
    "A3":     ("erfc", {"tau_c": 0.90, "delta": 0.05}),
    # A4: taramadan SONRA eklendi. A1 ve A2, gerçek σ_η² ve c₁'i iki yandan
    # kuşattı; ikisinde de doğrusal ara değer τ_c≈0.66–0.69 veriyor. A4 tam
    # o noktadır: S3 istatistiklerini GERÇEĞE OTURTAN erfc kesimi. Soru —
    # oraya oturunca faz da kapanıyor mu?
    "A4":     ("erfc", {"tau_c": 0.68, "delta": 0.125}),
    "B1":     ("jitter", {"sigma_j": 0.05}),
    "B2":     ("jitter", {"sigma_j": 0.10}),
    # C1: taramadan SONRA eklendi. "İkisi de kısmi mi?" sorusunu doğrudan
    # sınar: H-A'nın S3'ü oturtan kesimi (A4) + H-B'nin en güçlü titreşimi
    # birlikte. Newton'u erfc merdiveniyle kendi çözer, sonra titreşimi
    # ekler.
    "C1":     ("erfc+jitter", {"tau_c": 0.68, "delta": 0.125,
                               "sigma_j": 0.10}),
}

GERCEK = {0.5375: (0.787, 0.206, 0.40), 0.585: (0.550, 0.488, 1.07),
          0.660: (0.118, 0.614, 1.85), 0.740: (-0.148, 0.498, 2.12),
          0.815: (-0.103, 0.655, 1.50)}
BANTLAR = [(0.525, 0.55), (0.55, 0.62), (0.62, 0.70),
           (0.70, 0.78), (0.78, 0.85)]


def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8


def rvm_d(t):
    return np.log(t / TWO_PI) / TWO_PI


def pk(lim):
    from sympy import primerange
    out = []
    for p in primerange(2, lim + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))


# ---------- S(z) ve S'(z): 151b ile aynı formül, bloklu ----------
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


def main(ad):
    tip, par = KONFIG[ad]
    t_bas = time.time()
    print(f"=== KONFIG {ad}  ({tip} {par}) ===", flush=True)

    d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
    Z = np.sort(np.asarray(d["zeros"], dtype=float))
    zr = Z[len(Z) - 300000:]
    t0, t1 = float(zr[0]), float(zr[-1])
    L_hedef = float(np.log(0.5 * (t0 + t1) / TWO_PI))

    # ---- merdiven (151b: τ ≤ 1.00) ----
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

    # ---- H-A: erfc penceresi (TEK DEĞİŞİKLİK) ----
    if tip in ("erfc", "erfc+jitter"):
        tc, dl = par["tau_c"], par["delta"]
        tau = om_l / L_hedef
        wgt = np.array([0.5 * math.erfc((x - tc) / dl) for x in tau])
        a_ham = a_l.copy()
        a_l = a_l * wgt
        print(f"  erfc-kesim τ_c={tc} Δ={dl}: Σa·w/Σa="
              f"{a_l.sum()/a_ham.sum():.4f}  (w>0.01 olan çizgi: "
              f"{int((wgt>0.01).sum())})", flush=True)

    # ---- H-B: kayıtlı keskin çözümü oku, Newton'u atla ----
    z_onbellek = SCRATCH / "z_keskin.npy"
    if tip == "jitter":
        if not z_onbellek.exists():
            raise SystemExit("HATA: önce 'keskin' konfigürasyonu koşulmalı "
                             f"({z_onbellek} yok).")
        z = np.load(z_onbellek)
        print(f"  keskin Newton çözümü okundu ({len(z)} nokta) — merdiven "
              "özdeş olduğu için Newton tekrarlanmadı.", flush=True)
        ns = None
    else:
        n0 = int(np.ceil(rvm_N(t0))); n1 = n0 + 300000
        ns = np.arange(n0, n1, dtype=float)
        # 1) pürüzsüz sayım (S'siz Newton)
        z = np.full_like(ns, 0.5 * (t0 + t1))
        for _ in range(30):
            F = rvm_N(z) - ns
            z = np.clip(z - F / rvm_d(z), 100.0, None)
            if np.max(np.abs(F)) < 1e-9:
                break
        print(f"  pürüzsüz çözüm: maks|F|="
              f"{np.max(np.abs(rvm_N(z)-ns)):.2e}", flush=True)
        # 2) öz-tutarlı Newton (sönümlü, kelepçeli) — havuzla
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
                aF = np.abs(F)
                mf = float(aF.max())
                if it % 2 == 0 or it == NEWTON_IT - 1 or mf < 1e-3:
                    print(f"  Newton {it:2d}: maks|F|={mf:.3e}  "
                          f"medyan|F|={np.median(aF):.3e}  "
                          f"%99|F|={np.percentile(aF, 99):.3e}  "
                          f"({time.time()-t_bas:.0f}s)", flush=True)
                if mf < 1e-3:
                    break
        finally:
            pool.close(); pool.join()
    z = np.sort(z)
    if tip == "kontrol":
        np.save(z_onbellek, z)
        print(f"  z kaydedildi -> {z_onbellek}", flush=True)

    # ---- H-B: inkoherent seyreltme (TEK DEĞİŞİKLİK) ----
    if tip in ("jitter", "erfc+jitter"):
        sj = par["sigma_j"]
        gbar_t = 1.0 / rvm_d(z.mean())      # = 2π/L
        rj = np.random.default_rng(7)
        z = z + rj.normal(0.0, sj * gbar_t, size=len(z))
        capraz = int(np.sum(np.diff(z) < 0))
        z = np.sort(z)
        print(f"  titreşim σ_j={sj} (σ={sj*gbar_t:.4f}, ḡ_t={gbar_t:.4f}), "
              f"rng(7); sıra bozan çift: {capraz} "
              f"({100.0*capraz/len(z):.2f}%)", flush=True)

    # ================= GERÇEK ZİNCİRİN AYNISI (151b) =================
    g = np.diff(z)
    mid = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(mid / TWO_PI)
    L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    Nn = len(ds)
    tt = (mid - mid.mean()) / (mid[-1] - mid[0])
    qs = pk(min(int(np.exp(0.52 * L)), 720))
    fr = np.array([np.log(q) for q in qs])
    C = 3 + 2 * len(fr)
    XtX = np.zeros((C, C)); Xty = np.zeros(C)
    for s0 in range(0, Nn, 40000):
        sl = slice(s0, min(s0 + 40000, Nn))
        arg = np.outer(mid[sl], fr)
        Xc = np.empty((sl.stop - sl.start, C))
        Xc[:, 0] = 1; Xc[:, 1] = tt[sl]; Xc[:, 2] = tt[sl]**2
        Xc[:, 3::2] = np.cos(arg); Xc[:, 4::2] = np.sin(arg)
        XtX += Xc.T @ Xc; Xty += Xc.T @ ds[sl]
        del Xc, arg
    b = np.linalg.solve(XtX, Xty)
    eta = np.empty(Nn)
    for s0 in range(0, Nn, 40000):
        sl = slice(s0, min(s0 + 40000, Nn))
        arg = np.outer(mid[sl], fr)
        eta[sl] = ds[sl] - (b[0] + b[1]*tt[sl] + b[2]*tt[sl]**2 +
                            np.cos(arg) @ b[3::2] + np.sin(arg) @ b[4::2])
        del arg
    c1 = float(np.mean(eta[:-1] * eta[1:]))
    s_ds = float(np.var(ds)); s_eta = float(np.var(eta))
    print(f"\nS3 — {ad}: σ_ds²={s_ds:.4f}  σ_η²={s_eta:.4f}  c₁={c1:+.5f}"
          f"   [gerçek: 0.1674 / 0.0227 / −0.01158]", flush=True)

    dsA = 0.5 * (ds[:-1] + ds[1:]); dsA -= dsA.mean()
    allq = pk(int(np.exp(0.86 * L)))
    allw = np.array([np.log(q) for q in allq])
    T = mid[-1] - mid[0]; dres = TWO_PI / T
    e0, e1 = eta[:-1], eta[1:]
    m0 = mid[:-1]
    rng = np.random.default_rng(21)
    print("\nS1/S2 — bant   sentetik Γrot(Re,Im)  R_sent  | gerçek Γrot  R_ger",
          flush=True)
    satirlar = []
    for lo, hi in BANTLAR:
        tumu = [(q, w) for q, w in zip(allq, allw) if lo < w / L <= hi]
        cand = tumu
        if len(cand) > 220:
            idx = rng.choice(len(cand), 220, replace=False)
            cand = [cand[i] for i in idx]
        crN = 0j; crO = 0j; on0 = off0 = 0.0
        nkul = 0
        for q, w in cand:
            j = np.searchsorted(allw, w)
            koms = [allw[k] for k in (j - 1, j + 1)
                    if 0 <= k < len(allw) and abs(allw[k] - w) > 1e-12]
            gap = min(abs(w - k) for k in koms)
            if gap < 2.5 * dres:
                continue
            nkul += 1
            for W, hedef in ((w, True), (w + gap / 2, False)):
                cw, sw = np.cos(W * m0), np.sin(W * m0)
                zc = complex(2 * np.mean(e0 * cw), -2 * np.mean(e0 * sw))
                zp = complex(2 * np.mean(e1 * cw), -2 * np.mean(e1 * sw))
                cr = zp * np.conj(zc) * np.exp(1j * TWO_PI * W / L)
                if hedef:
                    crN += cr; on0 += abs(zc)**2
                else:
                    crO += cr; off0 += abs(zc)**2
        tb = 0.5 * (lo + hi)
        gr = GERCEK[round(tb, 4)]
        # Ayrışım paydası sıfırsa (çözünürlük eleği bandı boşaltmışsa) bant
        # ÖLÇÜLEMEZ — koşuyu düşürmeden '—' olarak kaydet.
        if nkul == 0 or abs(on0 - off0) < 1e-300:
            print(f"  {tb:.4f}  (—,—)        —      "
                  f"| ({gr[0]:+.3f},{gr[1]:+.3f})  {gr[2]:.2f}    "
                  f"[ölçülemedi: kullanılan aday {nkul}, payda "
                  f"{on0-off0:.2e}]", flush=True)
            satirlar.append({"tau": round(tb, 4), "re": None, "im": None,
                             "R": None, "R_artik": None, "n_aday": nkul})
            continue
        G = (crN - crO) / (on0 - off0)
        A = TWO_PI * tb
        phm = float(np.angle(G))
        best = None
        for R in np.linspace(0.0, 3.0, 121):
            wg = np.exp(-A * R * dsA)
            M = np.mean(wg * np.exp(-1j * A * dsA)) / np.mean(wg)
            if best is None or abs(np.angle(M) - phm) < best[0]:
                best = (abs(np.angle(M) - phm), R)
        # R eşleşmesinin gerçekten tutup tutmadığını da kaydet (faz sarması
        # durumunda 'en iyi' R yalnız tavana/tabana yapışır)
        artik = best[0]
        print(f"  {tb:.4f}  ({G.real:+.3f},{G.imag:+.3f})   {best[1]:4.2f}   "
              f"| ({gr[0]:+.3f},{gr[1]:+.3f})  {gr[2]:.2f}    "
              f"[R-artık {artik:.3f}, aday {nkul}]", flush=True)
        satirlar.append({"tau": round(tb, 4), "re": G.real, "im": G.imag,
                         "R": best[1], "R_artik": artik, "n_aday": nkul})

    ozet = {"ad": ad, "tip": tip, "par": par, "newton_it": NEWTON_IT,
            "sigma_ds2": s_ds, "sigma_eta2": s_eta, "c1": c1,
            "bantlar": satirlar, "sure_s": time.time() - t_bas}
    (SCRATCH / f"ozet_{ad}.json").write_text(json.dumps(ozet, indent=1))
    print(f"\n[{ad}] bitti — {(time.time()-t_bas)/60:.1f} dk", flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in KONFIG:
        raise SystemExit(f"kullanım: 152_gaz.py <{'|'.join(KONFIG)}>")
    main(sys.argv[1])
