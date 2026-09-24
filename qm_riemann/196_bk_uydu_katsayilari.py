# -*- coding: utf-8 -*-
"""
196 — Bogomolny–Keating çift-korelasyon sanısının TARAK UYDUSU katsayıları
(24 Eylül 2026; KEŞİF — veri-SONRASI, parametresiz karşılaştırma; ön-kayıtlı DEĞİL)

Kaynak formül (Bogomolny 2007, arXiv:0708.4223, eş. 7.5–7.6; BK 1996):
  R2_off(ε) = (1/4π²) |ζ(1+iε)|² e^{2πi d̄ ε} Φ_off(ε) + c.c.,
  Φ_off(ε)  = Π_p [1 − (1 − p^{iε})²/(p−1)²],  2π d̄ = L = log(E/2π).
Taşıyıcı e^{iLε}, yapı çarpanında tarağı (ω = L) verir; f(ε) = |ζ(1+iε)|² Φ_off(ε)
hemen-periyodiktir ve Euler çarpımıdır: f = Π_p g_p(θ_p), θ_p = ε log p,
  g_p(θ) = |1 − e^{−iθ}/p|^{−2} · [1 − (1 − e^{iθ})²/(p−1)²].
Bohr katsayısı r = Π p^{k_p} frekansında c(r) = Π_p ĝ_p(k_p) ⇒ uydu Δω = log r'nin
göreli genliği c(r)/c(1) = Π_{p | r} ĝ_p(k_p)/ĝ_p(0). Serbest parametre YOK.
"""
from fractions import Fraction

import numpy as np
from sympy import factorint

M = 1 << 14
TH = 2 * np.pi * np.arange(M) / M


def ghat(p):
    z = np.exp(1j * TH)
    g = np.abs(1 - np.conj(z) / p) ** -2 * (1 - (1 - z) ** 2 / (p - 1) ** 2)
    return np.fft.fft(g) / M          # [k] = e^{ikθ} katsayısı


G = {p: ghat(p) for p in (2, 3, 5, 7, 11, 13)}


def c(r):
    r = Fraction(r)
    k_p = {**factorint(r.numerator),
           **{q: -e for q, e in factorint(r.denominator).items()}}
    out = 1 + 0j
    for p, k in k_p.items():
        out *= G[p][k % M] / G[p][0]
    return out


if __name__ == "__main__":
    hedef = {"+log2": 2, "+log3": 3, "+log5": 5, "+log6": 6, "+log7": 7, "+log10": 10,
             "+log4": 4, "+log(3/2)": Fraction(3, 2), "-log2": Fraction(1, 2),
             "-log3": Fraction(1, 3)}
    cc = {k: c(v) for k, v in hedef.items()}
    for k, v in cc.items():
        print(f"{k:10s} c/c(1) = {v.real:+.4f} {v.imag:+.1e}i")
    # 192 (H-192c) düşük pencere çekirdek oranları — KARŞILAŞTIRMA (veri-sonrası)
    olc = {"+log3": 0.310, "+log5": 0.089, "+log6": 0.602, "+log7": 0.049, "+log10": 0.217}
    for k, m in olc.items():
        print(f"{k:8s} ölçülen/log2 {m:.3f}  BK/log2 {abs(cc[k]) / abs(cc['+log2']):.3f}")
    print(f"-log2/+log2 ölçülen 0.240  BK {abs(cc['-log2']) / abs(cc['+log2']):.3f}")
    print(f"-log3/+log3 ölçülen 0.466  BK {abs(cc['-log3']) / abs(cc['+log3']):.3f}")
