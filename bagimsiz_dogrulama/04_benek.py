"""
BAGIMSIZ BENEK YASASI TESTI (Not 4):
G_hat(omega) = <exp(i omega t_n)>  (t_n = aralik orta noktalari)
Yasa: |G_hat(log q)| = Lambda(q)/(L sqrt(q)) * |cos(pi tau)| * DW,  tau=log q/L
"""
import numpy as np, json

dz = np.load("../qm_riemann/128_odl_zeros6_2e6_zeros.npz")
zeros = np.sort(np.asarray(dz["zeros"], float))
TWO_PI = 2*np.pi

def asal_kuvvetleri(limit):
    """(q, p, k, Lambda(q)=log p) listesi, q<=limit."""
    asallar = []
    for n in range(2, limit+1):
        if all(n % d for d in range(2, int(n**0.5)+1)): asallar.append(n)
    out = []
    for p in asallar:
        q = p; k = 1
        while q <= limit:
            out.append((q, p, k, np.log(p))); q *= p; k += 1
    return sorted(out)

merkez, n_gap = 1.2e5, 40000
i0 = int(np.argmin(np.abs(zeros-merkez)))
g = zeros[max(0, i0-n_gap//2): i0-n_gap//2+n_gap+1]
mid = 0.5*(g[:-1]+g[1:]); mid = mid - mid.mean()
L = np.log(np.mean(g)/TWO_PI)
print(f"pencere: {len(mid)} orta nokta, L={L:.4f}\n")
print(f"{'q':>6} {'tau':>6} {'|G|':>10} {'yasa(cos yok)':>14} {'oran':>8} {'faz(deg)':>9}")
sat = []
for (q,p,k,Lam) in asal_kuvvetleri(60):
    tau = np.log(q)/L
    G = np.mean(np.exp(1j*np.log(q)*mid))
    yasa = Lam/(L*np.sqrt(q))
    oran = np.abs(G)/yasa
    sat.append((q, tau, abs(G), yasa, oran, np.degrees(np.angle(G))))
    print(f"{q:6d} {tau:6.3f} {abs(G):10.5f} {yasa:14.5f} {oran:8.3f} {np.degrees(np.angle(G)):9.2f}")
r = np.array([s[4] for s in sat]); tau = np.array([s[1] for s in sat])
print(f"\ncos(pi tau) ile duzeltilmis oran (DW adayi):")
for (q,t_,a,y,o,f) in sat:
    c = np.cos(np.pi*t_)
    if abs(c) > 0.05: print(f"  q={q:3d} tau={t_:.3f}  oran/cos = {o/c:7.4f}")
json.dump(dict(L=float(L), n=len(mid), satir=[dict(q=int(a),tau=float(b),G=float(c),yasa=float(d),oran=float(e),faz=float(f)) for a,b,c,d,e,f in sat]), open("04_sonuc.json","w"), indent=1)
