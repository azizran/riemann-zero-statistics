"""
BAGIMSIZ CUE TESTI: Haar-rastgele U(N) karakteristik polinomu
|Lambda(theta)| = prod_j 2|sin((theta-phi_j)/2)| ; aralik ici maksimum
ile aralik arasindaki Pearson/Spearman korelasyonu.
"""
import numpy as np, json, sys, time

def cue_ornek(N, mat_sayisi, k=24, rng=None):
    rng = rng or np.random.default_rng(1)
    G = np.zeros(mat_sayisi*N); Mg = np.zeros(mat_sayisi*N); idx = 0
    for _ in range(mat_sayisi):
        A = (rng.normal(size=(N, N)) + 1j*rng.normal(size=(N, N)))/np.sqrt(2)
        Q, R = np.linalg.qr(A); d = np.diagonal(R); Q = Q*(d/np.abs(d))
        phi = np.sort(np.angle(np.linalg.eigvals(Q)))
        # aralik ici maksimum: her aralikta k nokta
        gaps = np.diff(np.concatenate([phi, [phi[0]+2*np.pi]]))
        # log|Lambda| izgara uzerinde: theta matrisi (N, k)
        for j in range(N):
            th = phi[j] + gaps[j]*(np.arange(1, k+1)/(k+1))
            D = 2*np.abs(np.sin((th[:, None] - phi[None, :])/2))
            with np.errstate(divide='ignore'):
                L = np.log(np.maximum(D, 1e-300)).sum(axis=1)
            Mg[idx+j] = np.exp(L.max())
        G[idx:idx+N] = gaps; idx += N
    return G, Mg

if __name__ == "__main__":
    mat = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    Ns = list(range(5, 23))
    sonuc = {}
    t0 = time.time()
    print(f"{'N':>3} {'gap':>9} {'r_Pearson':>10} {'r_Spearman':>11}")
    for N in Ns:
        rng = np.random.default_rng(1000+N)
        g, M = cue_ornek(N, mat, rng=rng)
        gt = g*N/(2*np.pi)                       # unfold: ortalama 1
        rp = float(np.corrcoef(gt, M)[0, 1])
        rg = np.argsort(np.argsort(gt)); rm = np.argsort(np.argsort(M))
        rs = float(np.corrcoef(rg, rm)[0, 1])
        sonuc[N] = dict(n=int(g.size), r_pearson=rp, r_spearman=rs)
        print(f"{N:3d} {g.size:9d} {rp:10.4f} {rs:11.4f}   [{time.time()-t0:.0f}s]")
    json.dump(sonuc, open(f"05_cue_sonuc_{mat}.json", "w"), indent=1)
