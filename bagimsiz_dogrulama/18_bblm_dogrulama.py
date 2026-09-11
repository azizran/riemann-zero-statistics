# -*- coding: utf-8 -*-
"""
BBLM SABITLERI — bagimsiz dogrulama + repo transkripsiyon hatasinin gosterimi.
BBLM (arXiv math/0602270):  c_n = [(-1)^n/(2n)!] sum_p (log p)^(2(n+1)) sum_r (r-1) r^(2n)/p^r
  => c_0 = sum_p (log p)^2 / (p-1)^2
  Lambda = g0^2 + 2 g1 + c0 = 1.57314... ;  Q = sum_p log^3 p/(p-1)^2 ;  C = Q/Lambda = 1.4720...
"""
from mpmath import mp, mpf, log
mp.dps = 25
from sympy import primerange

PRIMES = list(primerange(2, 200000))
def topla(f): return sum(f(p) for p in PRIMES)

g0 = mpf('0.577215664901532860606512090082')
g1 = mpf('-0.072815845483676724860586375874')

c0_dogru = topla(lambda p: log(p)**2/(p-1)**2)
c0_repo  = topla(lambda p: log(p)**4*sum(mpf((r-1)*r*r)/mpf(p)**r for r in range(1, 40)))
Q        = topla(lambda p: log(p)**3/(p-1)**2)
Lam      = g0**2 + 2*g1 + c0_dogru
C        = Q/Lam

if __name__ == "__main__":
    import json
    print(f"c0 (doğru)      = {mp.nstr(c0_dogru, 12)}")
    print(f"c0 (repo hatası)= {mp.nstr(c0_repo, 12)}")
    print(f"Λ  = {mp.nstr(Lam, 10)}   (makale 1.57314)")
    print(f"√(12Λ) = {mp.nstr(mp.sqrt(12*Lam), 10)}")
    print(f"Q  = {mp.nstr(Q, 8)}")
    print(f"C  = {mp.nstr(C, 8)}   (makale 1.4720)")
    olcum = {5.60:6.52, 6.99:7.97, 8.39:9.45, 9.86:10.91, 10.93:12.08, 11.98:13.40}
    print(f"\n{'L':>7} {'BBLM Neff':>11} {'bizim Neff':>12} {'oran':>7}")
    for L, v in olcum.items():
        b = float(mpf(L)/mp.sqrt(12*Lam)); print(f"{L:7.2f} {b:11.3f} {v:12.2f} {v/b:7.2f}")
    json.dump(dict(c0=float(c0_dogru), c0_repo_transkripsiyon=float(c0_repo), Lambda=float(Lam),
                   Q=float(Q), C=float(C), sqrt12L=float(mp.sqrt(12*Lam))), open("18_bblm.json", "w"), indent=1)
