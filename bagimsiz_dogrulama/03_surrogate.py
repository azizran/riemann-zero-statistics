"""
BAĞIMSIZ SURROGATE NULL (Theiler yontemi):
Ayni guc tayfli Gaussian surecler uret, ayni izgaradan sifir+maksimum cikar,
ayni korelasyonu olc. Zeta de ayni izgaradan okunur (esit muamele).
"""
import numpy as np, json, sys, time
from importlib import import_module
m2 = import_module("02_olcum") if False else None
import importlib.util
spec = importlib.util.spec_from_file_location("olcum", "02_olcum.py"); olcum = importlib.util.module_from_spec(spec); spec.loader.exec_module(olcum)
Z, TWO_PI, A_CG = olcum.Z, olcum.TWO_PI, olcum.A_CG

def sifir_ve_maks(ts, ys, esik=0.0):
    """Izgaradan sifir gecisleri (lineer) + aralik ici |y| maksimumu (vektorize)."""
    a = ys[:-1]; b = ys[1:]
    idx = np.where(((a <= esik) & (b > esik)) | ((a > esik) & (b <= esik)))[0]
    if idx.size < 3: return None, None, None
    t0, t1 = ts[idx], ts[idx+1]; y0, y1 = ys[idx], ys[idx+1]
    zpos = t0 - y0*(t1-t0)/(y1-y0)                      # lineer interpolasyon
    absy = np.abs(ys)
    seg = idx[:-1] + 1                                  # her araligin ilk izgara noktasi
    M = np.maximum.reduceat(absy, seg)[:seg.size]       # son parca kirpilir
    d = np.diff(zpos)
    tmid = 0.5*(zpos[:-1] + zpos[1:])
    return d[:M.size], M, tmid[:M.size]

def korelasyon(d, M, tmid):
    L = np.log(np.mean(tmid)/TWO_PI)
    dt = d*L/TWO_PI
    Mt = M/np.sqrt(A_CG*L + 2.758 - 0.054/L)
    return float(np.corrcoef(dt, Mt)[0, 1]), float(L)

if __name__ == "__main__":
    n_sur = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    t0, t1 = 107252.0, 132748.0
    L_ref = np.log(0.5*(t0+t1)/TWO_PI); ort_gap = TWO_PI/L_ref
    h = ort_gap/32.0
    n = int((t1-t0)/h)
    ts = t0 + h*np.arange(n)
    yaz = lambda *a: print(*a, flush=True)
    yaz(f"izgara: {n} nokta, adim {h:.5f} (ort. aralik/{ort_gap/h:.0f}), L~{L_ref:.3f}")
    t = time.time(); yz = Z(ts); yaz(f"Z(t) hesaplandi [{time.time()-t:.1f}s]")

    d, M, tmid = sifir_ve_maks(ts, yz)
    r_zeta, L = korelasyon(d, M, tmid)
    yaz(f"ZETA (ayni izgaradan): {len(d)} gap, L={L:.3f}, r={r_zeta:.4f}")

    rng = np.random.default_rng(20260911)
    F = np.fft.rfft(yz); amp = np.abs(F)
    rs = []
    t = time.time()
    for k in range(n_sur):
        faz = rng.uniform(0, 2*np.pi, amp.size)
        faz[0] = 0.0
        if n % 2 == 0: faz[-1] = 0.0
        ys = np.fft.irfft(amp*np.exp(1j*faz), n=n)
        dd, MM, tt = sifir_ve_maks(ts, ys)
        rk, _ = korelasyon(dd, MM, tt)
        rs.append(rk)
        if (k+1) % 10 == 0: yaz(f"  surrogate {k+1}/{n_sur}: r={rk:.4f}  [{time.time()-t:.1f}s]")
    rs = np.array(rs)
    z = (r_zeta - rs.mean())/rs.std(ddof=1)
    yaz(f"\nSONUC  r_zeta = {r_zeta:.4f}   r_null = {rs.mean():.4f} ± {rs.std(ddof=1):.4f}  (n={n_sur})")
    yaz(f"       fark = {r_zeta-rs.mean():+.4f}  =>  {z:+.1f} sigma")
    json.dump(dict(t0=t0, t1=t1, h=h, n_izgara=n, r_zeta=r_zeta, r_null_ort=float(rs.mean()),
                   r_null_std=float(rs.std(ddof=1)), n_sur=n_sur, z=float(z), rs=rs.tolist()),
              open("03_sonuc.json", "w"), indent=1)
