# 203-b1k3: ratios-conjecture formula vs REAL Riemann zeros (independent numerical check)

## Verdict
- The 203_b1k3_zeta.py formula (Zeta(L), L = window mean 14.19425) reproduces the third-order joint zero-cumulants
  kappa_z(X(a),X(b),conj X(c)) and kappa_z(X(a),X(b),X(c)) of the LMFDB zeros within statistical errors at all 7 shift triples:
  all 14 pulls satisfy |pull| <= 1.1 (XXB: 0.0..0.6; XXX: 0.0..1.1). Relative data-formula difference 0.004-0.7 % (jackknife SE 0.2-1.8 %).
  Second-order warm-up (14 pairs): |pull| <= 0.7, relative agreement 1e-4..1e-3.
- CUE(N=L) is clearly worse: off by 9-26 % on the third cumulants, pulls 5 to 118 sigma (second order 27 to 394 sigma).
  => the formula (arithmetic-corrected) captures the real-zero value; CUE does not.
- Imaginary parts of all third cumulants are consistent with 0 (|Im/SE| <= 1.2).
- The mean E_z Re X(a) is reproduced by G(a)/L to ~1e-4 relative (checks the far-zero tail correction; table at end).

## Setup
- Data: qm_riemann/veri_lmfdb/zeros_8846000.dat (LMFDB file, not redistributed), first 1.5e6 zeros (window C1, read with 200c_veri.oku, base 8846000, offsets float64).
  Evaluated zeros: file index 800 .. 1499168 (n = 1,498,368, so every evaluated zero has a full +-800 zero neighbourhood). Window mean L = 14.194254 (matches HUKUM Lbar1).
- X_rho(a) = 1/a + sum_{0<|d|<=K} 1/(a + i(gamma - gamma')) - L_i/2 + (L_i/pi)*atan(a/W_i)   (last term = smooth-density far-zero real-part correction, W_i = half-width of the window),
  L_i = log(gamma_i/2pi) per zero, a = alpha/Lbar (alpha = shift in units 1/L). O(1/t) (t ~ 9e6) terms dropped. Window = +-K zeros (index-symmetric), K = 50,100,200,400,800 (W = K*0.4426 in gamma, i.e. K mean spacings; K=800 -> W = 354).
- Cumulants: exact joint-cumulant formula from the task, complex arithmetic (X-bar = complex conjugate), real part compared; delete-one-block jackknife, 32 blocks of consecutive zeros (also 64, 128 as SE check).
- Formula: Zeta(14.194254).S.G/T/Q/Q2 combined exactly as in Zeta.integrand (E_z moments with density weight); CUE: 203_b1k3_cue.funcs(N=14.194254) combined identically.
- Coincident nominal shifts: the formula (T, Q, Q2, DS) is singular at exact coincidence (a=b or a=b=c), so for (1,1,1),(3,3,2),(2,2,4) the two equal shifts are split by +-2 % (actual shifts printed in the table); data and formula are evaluated at exactly the same shifts.
  Formula sensitivity to this splitting is tiny and smooth (last section: delta 0.005..0.08 changes the formula by <= 0.002 rel. at delta=0.02).

## IMPORTANT: W (window) stability -- there IS a 1/K systematic in the data cumulants
The truncated-window data cumulants converge like c/K (NOT stable at a few hundred zeros): e.g. XXB(3,3) 3.1797, 3.2218, 3.2430, 3.2538, 3.2592 for K = 50,100,200,400,800
(increments halve each doubling). Origin: the imaginary (Hilbert-kernel, conditionally convergent) part -- the far-zero remainder is correlated with the near sum at O(1/K)
(cross term Cov(near, tail), not the tail variance, which is O(1/K^2)). Constant shifts (far-zero mean, real-part tail) do not matter (cumulants are shift invariant), so the
effect is genuinely the far-zero fluctuation. Size: K=800 is still low by 0.3-0.6 % (third order 0.04-0.13 abs; larger than the SE at several triples!). Therefore the data numbers
in the main table are Richardson-extrapolated to K -> infinity (quadratic in 1/K from K=200,400,800). Extrapolation is very stable: linear (400,800), quadratic (100,200,400),
cubic (all four) agree to <= 1e-4 abs. The extrapolation is validated independently by the second-order warm-up (formula known exactly there): extrapolated data agree with the formula
to 1e-4..1e-3 rel while raw K=800 values are off by 0.1-0.2 %. SE of the extrapolated value = SE at K=800 (the K-increments are almost perfectly correlated; SE(K) constant to 1e-3 rel across K).
Effect on pulls if NOT extrapolated: raw K=800 third-order pulls are -0.5..+1.6 (XXX all +1.0..+1.6), second-order up to -6.3 sigma ((3,3) XXB) and +4.3; raw K=200: third-order up to +4.6, second-order up to 24 sigma. So the raw finite-window numbers alone would look like a (spurious) few-sigma discrepancy at second order; K->infinity removes it.

## Tables (data = K->inf extrapolated; SE = 32-block jackknife at K=800)
Lbar = 14.194254; data: zeros_8846000.dat, first 1.5e6 zeros (window C1); evaluated zeros idx 800..1499168 (n=1498368 used); 32-block jackknife

## W-stability (window = ±K zeros, K=50..800; mean spacing 2π/L=0.4426; W=K·0.4426) – cumulant vs K
K-dependence ~ c/K (systematic). Richardson (200,400,800; quadratic in 1/K) used as K→∞ estimate.
| cumulant | K=50 | K=100 | K=200 | K=400 | K=800 | extrap(200,400,800) | extrap(400,800; lin) | extrap(100,200,400) | extrap(all4; cubic) |
|---|---|---|---|---|---|---|---|---|---|
| (1,1,1) XXB | 17.0528 | 17.1114 | 17.1409 | 17.1558 | 17.1632 | 17.1706 | 17.1706 | 17.1706 | 17.1706 |
| (1,1,1) XXX | -22.1661 | -22.3416 | -22.4303 | -22.4748 | -22.4972 | -22.5195 | -22.5195 | -22.5195 | -22.5195 |
| (2,1,0.5) XXB | 20.7515 | 20.8450 | 20.8922 | 20.9158 | 20.9277 | 20.9396 | 20.9396 | 20.9396 | 20.9396 |
| (2,1,0.5) XXX | -17.6417 | -17.7965 | -17.8746 | -17.9139 | -17.9335 | -17.9533 | -17.9532 | -17.9532 | -17.9533 |
| (0.5,2,1) XXB | 14.1555 | 14.1950 | 14.2149 | 14.2249 | 14.2299 | 14.2349 | 14.2349 | 14.2349 | 14.2349 |
| (0.5,2,1) XXX | -17.6417 | -17.7965 | -17.8746 | -17.9139 | -17.9335 | -17.9533 | -17.9532 | -17.9532 | -17.9533 |
| (3,3,2) XXB | 2.6767 | 2.6954 | 2.7049 | 2.7096 | 2.7120 | 2.7144 | 2.7144 | 2.7144 | 2.7144 |
| (3,3,2) XXX | -1.4482 | -1.4945 | -1.5179 | -1.5298 | -1.5357 | -1.5417 | -1.5416 | -1.5417 | -1.5417 |
| (1,4,0.5) XXB | 9.6110 | 9.6708 | 9.7010 | 9.7162 | 9.7238 | 9.7314 | 9.7314 | 9.7314 | 9.7314 |
| (1,4,0.5) XXX | -7.3298 | -7.4372 | -7.4914 | -7.5187 | -7.5323 | -7.5460 | -7.5460 | -7.5460 | -7.5460 |
| (0.3,0.6,0.4) XXB | 35.4090 | 35.4844 | 35.5224 | 35.5415 | 35.5511 | 35.5607 | 35.5607 | 35.5607 | 35.5607 |
| (0.3,0.6,0.4) XXX | -54.0609 | -54.2731 | -54.3802 | -54.4339 | -54.4608 | -54.4878 | -54.4878 | -54.4878 | -54.4878 |
| (2,2,4) XXB | 1.1182 | 1.1370 | 1.1465 | 1.1513 | 1.1537 | 1.1561 | 1.1561 | 1.1561 | 1.1561 |
| (2,2,4) XXX | -1.5630 | -1.6132 | -1.6386 | -1.6514 | -1.6579 | -1.6643 | -1.6643 | -1.6643 | -1.6643 |
| XXB2 (0.5,0.5) | 11.8200 | 11.8995 | 11.9396 | 11.9598 | 11.9700 | 11.9801 | 11.9801 | 11.9801 | 11.9802 |
| XXB2 (1.0,1.0) | 8.9990 | 9.0698 | 9.1055 | 9.1235 | 9.1326 | 9.1416 | 9.1416 | 9.1416 | 9.1416 |
| XXB2 (2.0,2.0) | 5.2424 | 5.2968 | 5.3243 | 5.3382 | 5.3452 | 5.3522 | 5.3521 | 5.3521 | 5.3522 |
| XXB2 (1.0,2.0) | 6.5637 | 6.6263 | 6.6579 | 6.6739 | 6.6819 | 6.6899 | 6.6899 | 6.6899 | 6.6899 |
| XXB2 (0.5,2.0) | 6.8648 | 6.9318 | 6.9656 | 6.9826 | 6.9912 | 6.9998 | 6.9997 | 6.9997 | 6.9998 |
| XXB2 (3.0,3.0) | 3.1797 | 3.2218 | 3.2430 | 3.2538 | 3.2592 | 3.2646 | 3.2646 | 3.2645 | 3.2646 |
| XXB2 (0.3,0.6) | 11.9035 | 11.9837 | 12.0242 | 12.0445 | 12.0548 | 12.0651 | 12.0650 | 12.0650 | 12.0651 |
| XXB2 (1.0,4.0) | 3.4235 | 3.4754 | 3.5016 | 3.5149 | 3.5215 | 3.5282 | 3.5282 | 3.5282 | 3.5282 |
| XX2 (1.0,2.0) | -3.7419 | -3.8045 | -3.8361 | -3.8520 | -3.8600 | -3.8680 | -3.8680 | -3.8680 | -3.8681 |
| XX2 (0.5,2.0) | -4.7842 | -4.8511 | -4.8849 | -4.9020 | -4.9105 | -4.9191 | -4.9191 | -4.9191 | -4.9191 |
| XX2 (1.0,4.0) | -1.7085 | -1.7603 | -1.7866 | -1.7998 | -1.8064 | -1.8131 | -1.8131 | -1.8131 | -1.8131 |
| XX2 (0.3,0.6) | -9.8859 | -9.9661 | -10.0065 | -10.0269 | -10.0372 | -10.0474 | -10.0474 | -10.0474 | -10.0474 |
| XX2 (0.5,1.0) | -7.4388 | -7.5139 | -7.5518 | -7.5709 | -7.5805 | -7.5901 | -7.5901 | -7.5901 | -7.5901 |
| XX2 (2.0,4.0) | -1.0477 | -1.0913 | -1.1134 | -1.1246 | -1.1302 | -1.1359 | -1.1358 | -1.1358 | -1.1359 |

## Main table: third-order joint zero-cumulants, real part (data = K→∞ extrapolated, SE = jackknife at K=800)
Actual shifts used (units 1/L) in brackets; coincident nominal shifts perturbed by ±2% (formula singular at exact coincidence), data and formula evaluated at the SAME shifts.
| triple, type | data ± SE | formula | pull_F | CUE(N=L) | pull_CUE |
|---|---|---|---|---|---|
| (1,1,1) [0.98,1.02,1] X X X̄ | +17.1706 ± 0.1049 | +17.1469 | +0.2 | +18.8611 | -16.1 |
| (1,1,1) [0.98,1.02,1] X X X | -22.5195 ± 0.0796 | -22.5988 | +1.0 | -24.3245 | +22.7 |
| (2,1,0.5) [2,1,0.5] X X X̄ | +20.9396 ± 0.0973 | +20.9161 | +0.2 | +23.0248 | -21.4 |
| (2,1,0.5) [2,1,0.5] X X X | -17.9533 ± 0.0574 | -18.0168 | +1.1 | -19.4974 | +26.9 |
| (0.5,2,1) [0.5,2,1] X X X̄ | +14.2349 ± 0.0790 | +14.2191 | +0.2 | +15.7223 | -18.8 |
| (0.5,2,1) [0.5,2,1] X X X | -17.9533 ± 0.0574 | -18.0168 | +1.1 | -19.4974 | +26.9 |
| (3,3,2) [2.94,3.06,2] X X X̄ | +2.7144 ± 0.0046 | +2.7145 | -0.0 | +3.2573 | -117.5 |
| (3,3,2) [2.94,3.06,2] X X X | -1.5417 ± 0.0057 | -1.5416 | -0.0 | -1.7077 | +29.1 |
| (1,4,0.5) [1,4,0.5] X X X̄ | +9.7314 ± 0.0358 | +9.7216 | +0.3 | +11.0170 | -35.9 |
| (1,4,0.5) [1,4,0.5] X X X | -7.5460 ± 0.0205 | -7.5654 | +0.9 | -8.4166 | +42.5 |
| (0.3,0.6,0.4) [0.3,0.6,0.4] X X X̄ | +35.5607 ± 0.6457 | +35.3115 | +0.4 | +38.9092 | -5.2 |
| (0.3,0.6,0.4) [0.3,0.6,0.4] X X X | -54.4878 ± 0.5370 | -54.7999 | +0.6 | -59.6118 | +9.5 |
| (2,2,4) [1.96,2.04,4] X X X̄ | +1.1561 ± 0.0055 | +1.1526 | +0.6 | +1.4535 | -54.2 |
| (2,2,4) [1.96,2.04,4] X X X | -1.6643 ± 0.0058 | -1.6652 | +0.2 | -1.8579 | +33.6 |

## Imaginary parts of the same third cumulants (K=800; expected ≈ 0 by symmetry)
| triple, type | Im part ± SE | Im/SE |
|---|---|---|
| (1,1,1) XXB | -0.0008 ± 0.0467 | -0.0 |
| (1,1,1) XXX | -0.0381 ± 0.0342 | -1.1 |
| (2,1,0.5) XXB | -0.0160 ± 0.0407 | -0.4 |
| (2,1,0.5) XXX | -0.0306 ± 0.0249 | -1.2 |
| (0.5,2,1) XXB | -0.0078 ± 0.0395 | -0.2 |
| (0.5,2,1) XXX | -0.0306 ± 0.0249 | -1.2 |
| (3,3,2) XXB | -0.0002 ± 0.0042 | -0.0 |
| (3,3,2) XXX | -0.0011 ± 0.0088 | -0.1 |
| (1,4,0.5) XXB | -0.0126 ± 0.0234 | -0.5 |
| (1,4,0.5) XXX | -0.0156 ± 0.0131 | -1.2 |
| (0.3,0.6,0.4) XXB | -0.0748 ± 0.1407 | -0.5 |
| (0.3,0.6,0.4) XXX | -0.0588 ± 0.0897 | -0.7 |
| (2,2,4) XXB | +0.0009 ± 0.0028 | +0.3 |
| (2,2,4) XXX | -0.0012 ± 0.0086 | -0.1 |

## Second-order warm-up (real part; data K→∞ extrapolated; SE at K=800)
| pair, type | data ± SE | formula | pull_F | CUE(N=L) | pull_CUE |
|---|---|---|---|---|---|
| X X̄ (0.5,0.5) | +11.9801 ± 0.0308 | +11.9795 | +0.0 | +12.8742 | -29.0 |
| X X̄ (1.0,1.0) | +9.1416 ± 0.0117 | +9.1441 | -0.2 | +9.8358 | -59.1 |
| X X̄ (2.0,2.0) | +5.3522 ± 0.0028 | +5.3532 | -0.4 | +5.8292 | -169.2 |
| X X̄ (1.0,2.0) | +6.6899 ± 0.0055 | +6.6914 | -0.3 | +7.2341 | -98.9 |
| X X̄ (0.5,2.0) | +6.9998 ± 0.0075 | +7.0008 | -0.1 | +7.5439 | -72.1 |
| X X̄ (3.0,3.0) | +3.2646 ± 0.0009 | +3.2648 | -0.3 | +3.6202 | -393.6 |
| X X̄ (0.3,0.6) | +12.0651 ± 0.0333 | +12.0633 | +0.1 | +12.9623 | -26.9 |
| X X̄ (1.0,4.0) | +3.5282 ± 0.0018 | +3.5283 | -0.0 | +3.8854 | -195.4 |
| X X (1.0,2.0) | -3.8680 ± 0.0047 | -3.8697 | +0.4 | -4.0765 | +44.6 |
| X X (0.5,2.0) | -4.9191 ± 0.0055 | -4.9223 | +0.6 | -5.2209 | +54.4 |
| X X (1.0,4.0) | -1.8131 ± 0.0020 | -1.8134 | +0.2 | -1.9217 | +54.6 |
| X X (0.3,0.6) | -10.0474 ± 0.0160 | -10.0573 | +0.6 | -10.7389 | +43.2 |
| X X (0.5,1.0) | -7.5901 ± 0.0101 | -7.5973 | +0.7 | -8.0706 | +47.7 |
| X X (2.0,4.0) | -1.1359 ± 0.0012 | -1.1355 | -0.3 | -1.1839 | +39.2 |

## Subset consistency (L-mixing / statistical): halves and middle 50% of the window, extrapolated data
| item | full | half1 | half2 | mid50 | formula |
|---|---|---|---|---|---|
| (1,1,1) XXB | 17.171±0.105 | 16.980±0.126 | 17.362±0.120 | 17.076±0.132 | 17.147 |
| (1,1,1) XXX | -22.520±0.080 | -22.553±0.103 | -22.486±0.118 | -22.461±0.101 | -22.599 |
| (2,1,0.5) XXB | 20.940±0.097 | 20.765±0.120 | 21.115±0.113 | 20.846±0.126 | 20.916 |
| (2,1,0.5) XXX | -17.953±0.057 | -17.971±0.075 | -17.936±0.087 | -17.900±0.071 | -18.017 |
| (0.5,2,1) XXB | 14.235±0.079 | 14.095±0.097 | 14.375±0.091 | 14.166±0.103 | 14.219 |
| (0.5,2,1) XXX | -17.953±0.057 | -17.971±0.075 | -17.936±0.087 | -17.900±0.071 | -18.017 |
| (3,3,2) XXB | 2.714±0.005 | 2.714±0.009 | 2.715±0.007 | 2.708±0.008 | 2.714 |
| (3,3,2) XXX | -1.542±0.006 | -1.550±0.011 | -1.533±0.006 | -1.543±0.011 | -1.542 |
| (1,4,0.5) XXB | 9.731±0.036 | 9.681±0.048 | 9.782±0.043 | 9.688±0.048 | 9.722 |
| (1,4,0.5) XXX | -7.546±0.021 | -7.558±0.028 | -7.534±0.032 | -7.529±0.028 | -7.565 |
| (0.3,0.6,0.4) XXB | 35.561±0.646 | 34.390±0.773 | 36.733±0.771 | 35.454±0.877 | 35.312 |
| (0.3,0.6,0.4) XXX | -54.488±0.537 | -53.648±0.627 | -55.329±0.719 | -53.918±0.597 | -54.800 |
| (2,2,4) XXB | 1.156±0.005 | 1.153±0.008 | 1.159±0.007 | 1.146±0.007 | 1.153 |
| (2,2,4) XXX | -1.664±0.006 | -1.672±0.011 | -1.656±0.006 | -1.666±0.011 | -1.665 |

## Jackknife block-count check of SE (K=800): NB=32 / 64 / 128
| item | SE32 | SE64 | SE128 |
|---|---|---|---|
| (1,1,1) XXB | 0.1049 | 0.0898 | 0.1005 |
| (1,1,1) XXX | 0.0796 | 0.0778 | 0.0844 |
| (2,1,0.5) XXB | 0.0973 | 0.0847 | 0.0951 |
| (2,1,0.5) XXX | 0.0574 | 0.0570 | 0.0616 |
| (0.5,2,1) XXB | 0.0790 | 0.0681 | 0.0754 |
| (0.5,2,1) XXX | 0.0574 | 0.0570 | 0.0616 |
| (3,3,2) XXB | 0.0046 | 0.0055 | 0.0066 |
| (3,3,2) XXX | 0.0057 | 0.0065 | 0.0074 |
| (1,4,0.5) XXB | 0.0358 | 0.0325 | 0.0362 |
| (1,4,0.5) XXX | 0.0205 | 0.0212 | 0.0236 |
| (0.3,0.6,0.4) XXB | 0.6457 | 0.5615 | 0.5603 |
| (0.3,0.6,0.4) XXX | 0.5370 | 0.4848 | 0.5027 |
| (2,2,4) XXB | 0.0055 | 0.0053 | 0.0059 |
| (2,2,4) XXX | 0.0058 | 0.0065 | 0.0075 |

## Formula near-coincidence stability (formula value vs perturbation δ of coincident shifts; XXB / XXX)
(1,1,1): δ=0.005: zeta 17.1492/-22.5999; δ=0.01: zeta 17.1487/-22.5984; δ=0.02: zeta 17.1469/-22.5988; δ=0.04: zeta 17.1394/-22.6014; δ=0.08: zeta 17.1089/-22.6095
(3,3,2): δ=0.005: zeta 2.7144/-1.5412; δ=0.01: zeta 2.7144/-1.5413; δ=0.02: zeta 2.7145/-1.5416; δ=0.04: zeta 2.7148/-1.5429; δ=0.08: zeta 2.7160/-1.5481
(2,2,4): δ=0.005: zeta 1.1529/-1.6650; δ=0.01: zeta 1.1528/-1.6651; δ=0.02: zeta 1.1526/-1.6652; δ=0.04: zeta 1.1517/-1.6659; δ=0.08: zeta 1.1481/-1.6686

## Mean check (far-zero tail correction), E_z Re X(alpha/L): formula G(a)/L vs data (K=800, with tail corr.)
alpha: formula, data -> 0.3: 40.853819, 40.853982 | 1: 8.923349, 8.923129 | 2: 3.012073, 3.011812 | 3: 1.446845, 1.446693 | 4: 0.825826, 0.825749  (all 14 shifts agree to <= 2.7e-4 abs, <= 1e-4 rel)

## L-variation across the window and subset check
Window L runs 14.158 -> 14.230 (+-0.036 about the mean). Formula slopes dkappa/dL (XXB): (1,1,1) 3.9, (2,1,0.5) 4.7, (0.5,2,1) 3.2, (3,3,2) 0.63, (1,4,0.5) 2.2, (0.3,0.6,0.4) 8.0, (2,2,4) 0.27 per unit L
(XXX: (1,1,1) -5.0, (2,1,0.5) -4.0, (3,3,2) -0.35). Using the window mean L is accurate to O(Var L) (negligible, <= 0.005). Half2-half1 differences in the subset table are consistent with this trend
(expected +0.14 for (1,1,1) XXB, observed +0.38 +- 0.17) plus noise. Middle-50 % subset (L range +-0.018) agrees with the full sample within 1 sigma.
SE robustness: jackknife SE with 32/64/128 blocks agree within about +-25 % (no systematic growth beyond that; (3,3,2)XXB 0.0046/0.0055/0.0066 is the largest spread -- pulls stay <= 0.3 sigma even with the largest).

## Caveats
- Pulls are highly correlated (all cumulants come from the same 1.5e6 zeros); the XXX pulls are all >= 0 (mean +0.6): data slightly less negative than formula by 0.06-0.08 (0.3-0.4 %) at the precise triples, still <= 1.1 sigma.
- 1/K extrapolation is an empirical (validated) correction; without it the K=800 numbers are 0.3-0.6 % low (third order: XXB up to 0.08, XXX up to 0.09 abs).
- Exact-coincidence shifts not tested directly (formula singular there); split by +-2 %.
- No repo files modified; scratch code/data in this directory: s203b1k3/{data_run.py, formula_run.py, extra_run.py, analyse.py, cumlib.py, shifts.py, data_out.json, formula_out.json, extra_out.json, report_tables.md, X400.npy}.
