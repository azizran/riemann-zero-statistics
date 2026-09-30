# 203-b1 kappa3 (log|zeta'(rho)| third cumulant) — numerics/derivation audit
(auditor: Sonnet 5.5; no repo files touched, no git; scratch code in scratchpad/a/jobs/*.py, logs *.out)

## Verdict table
| # | Item | Verdict | Key numbers |
|---|------|---------|-------------|
| 1 | 4-pt ratios recipe bookkeeping (Q, Q2, DS) | PASS | hand derivation + per-prime theta-quadrature agree to 1e-16 (K_j, f_jk, f_i'j', conj-side K, E_p closed forms, DS E/Y incl. (1-1/p)^4); A_DS -> 1 at zero shifts (1-A ~ scale^2: 2.4e-6 at 0.003x); CUE formulas vs Weyl N=1..4 (incl. shifts ~0.02): 8-9 digits; Zeta-class assembly fed CUE ingredients == cue.integrand to 6e-14 |
| 2 | b1 assembly, CUE Palm kappa3 | PASS | k3_b1_cue (amin_f=.001,n=12) vs exact sum: N=2 -3e-7, 3 -2e-7, 5 +9e-7, 10 +1e-7; N=9.35/11.66/22.3 (n=8) vs kap_vec: +7e-7 (strip); Palm formula independently re-derived from Keating-Snaith (d^3 log M_{N-1}(1+s/2)) |
| 3 | numerical scheme (kappa3^CUE + Int Delta + slab) | scheme PASS; error-budget claim "~1e-4" UNCERTAIN/optimistic | slab/edge/quadrature fine (<=5e-5); integrand implementation errors sum to -2.0e-4 (L=11.66) with gross terms up to 4.6e-4; see below |
| 4 | A_DS interpolation | limits PASS; interpolation-region size UNCERTAIN | a1->0: A_DS->A(x21) (0.920745 vs 0.920776 at a1=1e-4), a2->0: A(x11); interp is +1.1..1.4% above true product for x_m in .15-.44 -> bias +3e-4 (L=9.35), +7e-5 (11.66); DSVAR1-2 spread 6e-5 underestimates |
| 5 | asymptotic tails | f_jk tail: FAIL (bug, -3.1e-4); others PASS | see below |
| 6 | comparison to sealed data | PASS | window mixture <= 2e-4 (<0.1 sigma); systematic corrections do not change conclusion (chi2 11.9 -> 11.7) |

## 1. Recipe bookkeeping (PASS)
Derived from the CFZ ratios theorem (A'=S-complement + T^-, B' likewise, Z, A=Prod E_p/Y_p):
* 3+1 swap i: vanishing denominators 1/zeta(1-beta+delta), 1/zeta(1-alpha_i+gamma_i) each absorb exactly one of d/dbeta, d/dalpha_i (factor (-1) each, product +1); the remaining d/dalpha_j d/dalpha_k act on zeta(1+alpha_j-alpha_i)/zeta(1-alpha_i+gamma_j) * zeta(1+gamma_j+beta)/zeta(1+alpha_j+beta) * A. Product of first derivatives gives K K (K = (z'/z)(1+u) - (z'/z)(1+u+x) + d log A); cross-derivative of the zeta part is zero, only A contributes f_jk. F(x)=e^{-Lx}zeta(1+x)zeta(1-x)A(x): I verified E_p=(1-2/p+p^{-1-x})/(1-p^{-1+x}) by hand and Y_p=(1-1/p)^2 z(x)z(-x), hence A(x) exactly as coded.
* 2+2 single swap: the alpha_{i'} derivative hits zeta(1+alpha_{i'}-alpha_i) [num.] and zeta(1+alpha_{i'}+beta_j) [den.] (the zeta(1+alpha_{i'}+beta_{j'}) num./den. pair cancels), beta_{j'} likewise -> both K(u,x) with x=alpha_i+beta_j, u=alpha_{i'}-alpha_i / beta_{j'}-beta_j; E_p is invariant under (alpha<->beta, w<->wbar) so conj-side K is the same function (checked numerically: dbeta log E_p == Kcf(b1-b0,x) to 1e-12). Mixed derivative: only zeta(1+alpha_{i'}+beta_{j'}) -> (z'/z)'(1+x2) = H(x2)+Sum ey; f_{i'j'} = H(x2)+Sum_p (d_a d_b log E_p - ey) — consistent.
* Double swap: all 4 derivatives each kill one of the 4 vanishing denominators -> e^{-L Sigma} Prod z(1+-x)/[z(1+-Db)z(1+-Da)] * A_DS with Y_p containing (1-1/p)^4. Coded E_p assignment (Un,Ud,Cn,Cd) matches the recipe (numerators of zeta -> (1-w a)^{-1}).
* Numerics (own theta-quadrature, 4096 midpoints, independent of the authors' closed forms): p=2..997, random shifts: dj logE, d_j d_k logE vs closed forms 2e-16; single-swap (i=j=0): fab, d_alpha, d_beta at p=2,3,7,29,211 identical to 12 digits with the vectorised Zeta.fij (monkey-patched to one prime); DS: E_theta/Y_p vs code A_DS single prime 1e-12.
* Identity terms: 3+1: kappa4_p = l^4 Prod y_i/(1-y_i) exactly (all lower cumulants vanish) — no tail needed (>1e4: 4e-6/L); 2+2: H11H22+H12H21+Sum kappa4 with kappa4_p = l^4 E4 - e e - e e; E4 closed form verified by the authors, leading tail p^{-2-Sigma} verified below.
* CUE: reran 203_cue4.py (N=1,2,3: 8 digits; N=3 Q2 shows 1e-6 = grid aliasing of the 40^3 brute force; with M=80: 1e-9) and N=4 (M=56): Q 1e-9, Q2 3e-9; small shifts (0.02-0.06) at N=1,2 with M=2e4/1800: exact to 9 digits.
* Extra: the vectorised Zeta class (T,Q,Q2,DS with masks) fed with CUE functions (F,H,K=zp(u)-zp(u+x), DS with z-functions, all prime sums off) reproduces cue.integrand to 6e-14 at 8 random triples -> index mapping ok.
Caveat: the independent real-zero test validates 4-point pieces only at few-% level (they enter kappa_z at 1/L weight); exactness of Q/Q2 rests on the above derivation+per-prime checks.

## 2. b1 assembly (PASS)
E_z moments checked by hand: m(u)=G/L; U=T(u,v,0)/L; V=G(u+w)+(T(u,0,w)+T(w,0,u))/L; M3=T(a,b,c)+(Q(a,b,0;c)+Q2(a,b;c,0))/L; M3u=Q(a,b,c;0)/L; kappa(A,B,C)=E[ABC]-E[AB]E[C]-...+2E[A]E[B]E[C]; kappa3(Re Y)=1/4 Re[kappa(YYY)+3 kappa(YYYbar)]. Numbers above. Exact Palm: X=log|L'| under |Lambda_{N-1}(1)|^2 weight => cumulants of log M_{N-1}(1+s/2), kappa3 = Sum_{j=1}^{N-1}[psi''(j+2)-psi''(j+1)/4] (matches kap_vec incl. non-integer N).

## 3. Scheme, slab, edges, quadrature
* Quadrature: with the same integrand structure (prime cutoff 1e3, tails consistent) n=4,6,8 gives kappa3 = -0.047351, -0.047296, -0.047300 (L=11.66): 6->8 change 4e-6; (real code n=4->6: 5.5e-5). Quadrature error <= 5e-5.
* Slab/edges: with converged integrand (1e4 primes + continuous tails, exact K): the UNsymmetrised Delta does NOT vanish at a shift -> 0 in general (e.g. (b,c)=(0.3,0.5): Delta(0+)=-0.10; c->0 at (a,b)=(0.1,0.05): +1.7; edges Delta(s,s,c)=-3.6, Delta(s,c,s)=+2.6 at s=.0075,c=.05). The symmetrised Delta_sym vanishes ~linearly: Delta_sym(s,0.1,0.05)/s = 28, 40, 45 for s=.02,.0075,.004 (Delta ~ 70 s - 2175 s^2). Because the integration domain [s0,oo)^3 and the three slab grids are permutation-symmetric, the omitted s0/2*Sum_perm Delta(0,..) terms cancel exactly in the continuum; edge (two-small) terms are Delta_sym ~ 66 s -> total edge omission ~<3e-6; curvature error of trapezoid slab (c2 s0^3/6) ~ 8% of slab (4.5e-4) = 4e-5 worst. So slab treatment is JUSTIFIED (for the symmetrised integrand), but the statement "Delta -> 0 linearly as one shift -> 0" is only true for Delta_sym; unsymmetrised slabs are individually O(1e-3) (-3.8e-4,-3.8e-4,+1.2e-3 at L=11.66) and only their sum is meaningful.
* Near s0 the CODE's integrand is NOT accurate: at (0.0075,0.1,0.05): code Delta-integrand error ~ +0.21 (converged I = -2.274, code -2.486; Delta_exact=+0.40); at a=0.0015 the code integrand blows up (-11.8 vs converged -0.51). Sources (all separated numerically): (i) DS Euler product truncated at p<1000 (dominant, error ~ a^-3), (ii) f_jk tail bug, (iii) K spline tables. s0=0.0075 hides it because the integrated effect is small (below).
* Integrated systematic corrections to kappa3 (shift = corrected - code; n=6, s0=.0075; interior+slab):
  | L | f_jk tail bug | DS p<1e4+tail | f_ij rem. p>1e4 | K tables | Wc cut (est) | DS hi-region interp (est) | net |
  | 9.35 | -4.63e-4 | +1.84e-4 | +1.8e-5 | -7.9e-5 | -5e-5 | +3.1e-4 (+-1.5e-4) | -8e-5 |
  | 10.59| -3.70e-4 | +1.33e-4 | +1.6e-5 | -6.7e-5 | -2e-5 | +1.6e-4 | -1.5e-4 |
  | 11.66| -3.13e-4 | +1.02e-4 | +1.8e-5 | -5.9e-5 | -9e-6 | +6.6e-5 | -1.98e-4 |
  | 14.19| -2.25e-4 | +5.8e-5 | ~+1e-5 | -4.5e-5 | -1e-6 | ~+1e-5 | -1.9e-4 |
  | 16.58| ~-1.7e-4 | ~+3e-5 | | ~-3.5e-5 | 0 | 0 | -1.65e-4 |
  | 18.88| -1.41e-4 | ~+2e-5 | | ~-2.9e-5 | | | -1.4e-4 |
  | 22.31| -1.08e-4 | ~+1e-5 | | -2.3e-5 | | | -1.2e-4 |
  DIRECT CHECK (my one full 3D run, L=11.66, n=6, s0=0.0075, integrand with 1e4 primes + continuous prime tails for all sums, corrected f_jk tail, DS product to 1e4+tail, refined K tables): kappa3 = -0.049262 (interior +0.152998, slab +0.000450) vs authors' -0.049007 => shift -2.55e-4; additive estimate of the same four items = -2.52e-4 (additivity confirmed).
  => "systematic ~1e-4" is optimistic: gross errors are 3-5e-4 (they mostly cancel to -1e-4..-2e-4). All << SE (2.7e-3 ... 1.0e-3 at L=9.35..11.66).
* Wc=0.75 cut in zeta only (CUE uncut): CUE-only test of the cut: +6.8e-5 (9.35), +1.25e-5 (11.66), +1.5e-6 (14.19), 5e-9 (22.3); zeta effect ~0.7x.

## 4. A_DS interpolation
* Limits (b2=0, exact Euler product to 1e6): a1=1e-4: A_DS=0.920745 vs A(x21)=0.920776; a1=1e-5: 0.899227 vs 0.899231; mirror a2->0 -> A(x11) identically. PASS.
* Accuracy: Euler product truncated at p<1000 is biased HIGH: vs p<1e6 (rel.) 3e-5 (x=.18), 5e-4 (.35), 0.7-1.7% (x_pair-sum ~.7-.85, x_m .40-.42); not even converged at 1e6 for pair-sums >0.8. (Convergence needs max pair-sum 2a1+b1+b2 etc. <1, not x_lk<.45.)
* DS-only integrals (zeta DS contribution to kappa3, interior): L=9.35: lo(x_m<.3) -1.9033, mid(.3-.45) -0.1025, hi(>=.45) -2.35e-2; L=11.66: -1.129, -0.0396, -5.07e-3. DSVAR 1 vs 2 changes hi by 6.1e-5 (9.35), 1.3e-5 (11.66). But measured interp/true ratio in x_m<.44: 1.005 (x .1), 1.010 (.2), 1.014 (.3), 1.011 (.4) (max 2.5%) i.e. interpolation too high by ~1.3% => hi-region correction +3.1e-4 (9.35), +6.6e-5 (11.66), +-100%. UNCERTAIN (no independent A_DS beyond pair-sum 1 without analytic continuation).
* DS p<1e3 truncation in mid/lo region: +1.84e-4 (9.35) / +1.02e-4 (11.66) as in table.

## 5. Tails / p=2
* f_jk: `tail_log(3,sj+x)-tail_log(3,sj+2x)` in Zeta.fjk is Sum log^3 p p^-s but the summand is log^2 p*om*qj*qk: needs tail_log(2,.). Exact Sum_{1e4<p<1e6} vs code: 1.14e-5 vs 1.16e-4; 8.3e-5 vs 8.6e-4; 4.6e-5 vs 4.6e-4; (x=.1,u=-.3): .127 vs 1.79; k=2 gives 1.16e-5, 8.4e-5, 4.6e-5, .152. Effect on kappa3: -3.1e-4 (L=11.66), -4.6e-4 (9.35) (table). BUG.
* f_ij remainder (p>1e4, x_ij<0.5): R_ij=3e-5..2e-4 at moderate shifts (R/f ~1e-7), kappa3 effect +1.8e-5. For x_ij>0.5 the truncated prime sum is not the analytic continuation (drifts with cutoff) but F(x)~e^{-Lx} makes F*psum<=1e-4 there. PASS.
* c4 tail: exact tail(997..1e6) vs code -tail4: ratio 0.987 (1.3% error of a ~1e-3..4e-3 correction); effect ~1e-6. PASS.
* D tail (tail_log(3) correct there): ratio 0.974-0.988; effect 2e-7. Q identity beyond 1e4: 2e-10. PASS.
* K: p=2 exact (R2) fine; but table accuracy: absolute error up to 6e-6 (u~.01) at all x (R spline; refined tables: rms 8e-8), which F~1/x^2 amplifies to integrand errors 0.02-0.06 at a~s0; integrated -5.9e-5 (11.66). (Their 2.6e-6 "relative to max(1,|K|)" test hides this.)

## 6. Sealed-data comparison
* x=log(|Z'|*2pi/L_n): mixture kappa3 = E[k3(L)] + 3 Cov(mu(L),k2(L)) + kappa3(mu). Computed mu(L) (ratios, exact G): mu' = -0.0076/unit L (near 9-11), k2' ~ +0.038, k3'' ~ +2e-4: window corrections 3 mu' k2' Var(L) + 1/2 k3'' Var(L) = -2.0e-4 (W1, Var .26), -7e-5 (W2), -8e-5 (W3): < 0.1 sigma. Zero-weighted Lbar (9.3427 vs 9.35) irrelevant here.
* With the systematic corrections above (all 7 windows) pulls: 9.35 +0.84->+0.87, 10.59 -1.31->-1.24, 11.66 -1.89->-1.70, 14.19 -0.57->-0.42, 16.58 +2.24->+2.38, 18.88 +0.63->+0.72, 22.31 -0.45->-0.39; chi2 (rounded observed) 11.9 -> 11.7. Comparison is fair; conclusions unchanged.
* Note: tested at L in table only; W-window L spread assumed Var L as in earlier audit.

## Recommendations
1. Fix `tail_log(3,..)` -> `tail_log(2,..)` in Zeta.fjk (and mention the -2e-4..-5e-4 shift).
2. DS Euler product: use p<1e4-1e5 (or analytic tails) — p<1000 makes the integrand near s0 wrong by 10% of Delta.
3. Rebuild R tables denser (u-nodes 0.005, x log 45) — costs 4 min.
4. State that the slab argument holds for the symmetrised integrand; error budget: -2e-4 +- 1.5e-4 at L=11.66 (dominated by interpolation region at L<=10.6).
5. Optional: pre-registered value stands; corrected values: -0.04885, -0.04931, -0.04920, -0.04815, -0.04685, -0.04558, -0.04386 (L=9.35..22.31).
