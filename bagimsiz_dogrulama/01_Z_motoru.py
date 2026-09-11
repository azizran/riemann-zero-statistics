"""
BAĞIMSIZ Z(t) MOTORU — Uğur Sezen'in scriptlerinden bağımsız, sıfırdan yazıldı.
Riemann-Siegel: ana toplam + Phi_0 düzeltmesi.
Doğrulama: mpmath ile Z(t) = Re(e^{i theta} zeta(1/2+it)).
"""
import numpy as np

# --- Riemann-Siegel theta (asimptotik, t>=100 icin ~1e-13) ---
def theta(t):
    t = np.asarray(t, dtype=float)
    return (t/2.0)*np.log(t/(2.0*np.pi)) - t/2.0 - np.pi/8.0 \
           + 1.0/(48.0*t) + 7.0/(5760.0*t**3) + 31.0/(80640.0*t**5)

# --- Phi_0(z) = cos(pi(z^2/2+3/8))/cos(pi z), z=2p ; z=1/2 civarinda kararli ---
def phi0(z):
    z = np.asarray(z, dtype=float)
    A = np.pi*(z*z/2.0 + 3.0/8.0)
    B = np.pi*z
    cA, cB = np.cos(A), np.cos(B)
    out = np.empty_like(z)
    safe = np.abs(cB) > 1e-7
    out[safe] = cA[safe]/cB[safe]
    if (~safe).any():                      # z ~ 1/2 noktalari: Taylor
        e = z[~safe] - 0.5
        out[~safe] = 0.5 + e*(1.0 - np.pi*np.pi/4.0) + e*e*(np.pi*np.pi/2.0)
    return out

def Z_rs(t, blok=200000):
    """Vektorize Riemann-Siegel Z(t): ana toplam + (-1)^(N-1) (t/2pi)^(-1/4) Phi_0(2p)."""
    t = np.atleast_1d(np.asarray(t, dtype=float))
    a = np.sqrt(t/(2.0*np.pi))
    N = np.floor(a).astype(np.int64)
    p = a - N
    th = theta(t)
    out = np.empty_like(t)
    for i in range(0, t.size, blok):
        sl = slice(i, min(i+blok, t.size))
        tc, thc, Nc = t[sl], th[sl], N[sl]
        acc = np.zeros_like(tc)
        nmax = int(Nc.max())
        for n in range(1, nmax+1):
            m = n <= Nc
            if not m.any():
                break
            acc[m] += np.cos(thc[m] - tc[m]*np.log(n))/np.sqrt(n)
        r = 2.0*acc
        corr = ((-1.0)**((Nc-1) % 2)) * (tc/(2.0*np.pi))**(-0.25) * phi0(2.0*(a[sl]-Nc))
        out[sl] = r + corr
    return out

if __name__ == "__main__":
    import mpmath as mp
    mp.mp.dps = 30
    def Z_mm(t):
        return float(mp.re(mp.e**(1j*mp.im(mp.loggamma(0.25+0.5j*mp.mpf(t)))) * mp.zeta(mp.mpf('0.5')+1j*mp.mpf(t))))
    testler = [14.134725142, 100.0, 999.999, 12345.678, 100000.5, 500000.25, 1130000.75, 1000.0]
    print(f"{'t':>14} {'Z_rs':>16} {'Z_mpmath':>16} {'|fark|':>12}")
    maks = 0.0
    for t in testler:
        z1 = float(Z_rs(np.array([t]))[0]); z2 = Z_mm(t)
        maks = max(maks, abs(z1-z2))
        print(f"{t:14.3f} {z1:16.9f} {z2:16.9f} {abs(z1-z2):12.2e}")
    # rastgele 40 noktada tarama
    rng = np.random.default_rng(7)
    ts = rng.uniform(100, 1.13e6, 40)
    zr = Z_rs(ts); zm = np.array([Z_mm(float(x)) for x in ts])
    print(f"\n40 rastgele nokta (100..1.13e6): maks|fark| = {np.abs(zr-zm).max():.3e}, ort = {np.abs(zr-zm).mean():.3e}")
