"""
BAĞIMSIZ ÖLÇÜM: gap (δ_n) ve aralık-içi maksimum (M_n) — kendi motorumla.
Veri: qm_riemann/128_odl_zeros6_2e6_zeros.npz (Odlyzko zeros6, ham γ).
"""
import numpy as np, json, sys, time

TWO_PI = 2*np.pi
A_CG = (np.e**2 - 5)/2.0                 # 1.1945280...
ALFA_M1 = 5 - np.e**2 - 10*0.5772156649015329 + 2*np.e**2*0.5772156649015329
ALFA_0  = -0.423285

def theta(t):
    t = np.asarray(t, float)
    return t/2*np.log(t/TWO_PI) - t/2 - np.pi/8 + 1/(48*t) + 7/(5760*t**3) + 31/(80640*t**5)

def Z(t, blok=40000):
    """Vektörize Riemann-Siegel Z(t), N-gruplu; kendi yazimim."""
    t = np.atleast_1d(np.asarray(t, float))
    a = np.sqrt(t/TWO_PI); N = a.astype(np.int64); p = a - N
    th = theta(t); out = np.empty_like(t)
    for i in range(0, t.size, blok):
        sl = slice(i, min(i+blok, t.size))
        tt, NN, tt_th = t[sl], N[sl], th[sl]
        zz = np.zeros_like(tt)
        for Nv in np.unique(NN):
            m = NN == Nv
            n = np.arange(1, Nv+1)
            zz[m] = 2.0*(np.cos(tt_th[m, None] - tt[m, None]*np.log(n)[None, :]) @ (n**-0.5))
        cp = np.cos(TWO_PI*p[sl])
        cp = np.where(np.abs(cp) < 1e-8, 1e-8, cp)
        zz += ((-1.0)**((NN-1) % 2))*(tt/TWO_PI)**-0.25*np.cos(TWO_PI*(p[sl]**2 - p[sl] - 1/16))/cp
        out[sl] = zz
    return out

def M_aralik(g0, g1, k=24):
    """[g0,g1] aralığında |Z| maksimumu: k iç nokta + parabolik düzeltme."""
    h = (g1-g0)/(k+1)
    ts = g0 + h*np.arange(1, k+1)
    zv = np.abs(Z(ts))
    j = int(np.argmax(zv)); M = float(zv[j])
    if 0 < j < k-1:                       # parabolik tepe düzeltmesi
        y0, y1, y2 = zv[j-1], zv[j], zv[j+1]
        den = (y0 - 2*y1 + y2)
        if den != 0:
            d = 0.5*(y0 - y2)/den
            if abs(d) <= 1:
                M = float(y1 - 0.25*(y0-y2)*d)
    return M

def pencere(zeros, merkez, n_gap):
    i0 = int(np.argmin(np.abs(zeros - merkez)))
    lo = max(0, i0 - n_gap//2); hi = min(len(zeros)-1, lo + n_gap + 1)
    g = zeros[lo:hi]
    return g

def olc(zeros, merkez, n_gap, k=24):
    g = pencere(zeros, merkez, n_gap)
    d = np.diff(g)                                  # δ_n
    tmid = 0.5*(g[:-1] + g[1:])
    M = np.array([M_aralik(g[i], g[i+1], k) for i in range(len(d))])
    L = np.log(np.mean(tmid)/TWO_PI)
    dtilde = d*L/TWO_PI
    Mtilde = M/np.sqrt(A_CG*L + 2.758 - 0.054/L)
    r = float(np.corrcoef(dtilde, Mtilde)[0, 1])
    tahmin = A_CG*L + (ALFA_M1 + 2*A_CG) + (ALFA_0 + ALFA_M1)/L
    return dict(L=float(L), n=int(len(d)), r=r, meanM2=float(np.mean(M**2)),
                meanM2_tahmin=float(tahmin), sapma_yuzde=float(100*(np.mean(M**2)/tahmin-1)),
                ort_gap=float(np.mean(d)), tmid0=float(tmid[0]), tmid1=float(tmid[-1]))

if __name__ == "__main__":
    n_gap = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
    dz = np.load("../qm_riemann/128_odl_zeros6_2e6_zeros.npz")
    zeros = np.sort(np.asarray(dz["zeros"], float))
    print(f"yuklenen sifir: {zeros.size}, t: {zeros[0]:.3f} .. {zeros[-1]:.1f}")
    merkezler = [1.2e5, 2.0e5, 3.5e5, 6.0e5, 1.0e6]
    sonuc = []
    for mc in merkezler:
        t0 = time.time(); s = olc(zeros, mc, n_gap); s["merkez"] = mc; s["sure_s"] = round(time.time()-t0, 1)
        sonuc.append(s)
        print(f"  t~{mc:.0e}  L={s['L']:.3f}  n={s['n']}  r={s['r']:.4f}  "
              f"meanM2={s['meanM2']:.4f} (tahmin {s['meanM2_tahmin']:.4f}, {s['sapma_yuzde']:+.2f}%)  [{s['sure_s']}s]")
    Ls = np.array([s["L"] for s in sonuc]); m2 = np.array([s["meanM2"] for s in sonuc])
    a_fit, b_fit = np.polyfit(Ls, m2, 1)
    print(f"\nfit: meanM2 = {a_fit:.4f}*L + {b_fit:.4f}   (Conrey-Ghosh A = {A_CG:.5f}, sapma {100*(a_fit/A_CG-1):+.2f}%)")
    print(f"     b tahmini (HLPC): {ALFA_M1+2*A_CG:.4f}  |  b fit: {b_fit:.4f}")
    json.dump(dict(n_gap=n_gap, sonuc=sonuc, a_fit=float(a_fit), b_fit=float(b_fit), A_CG=float(A_CG)),
              open(f"02_sonuc_{n_gap}.json", "w"), indent=1)
