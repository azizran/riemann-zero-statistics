# -*- coding: utf-8 -*-
"""
176 — VEKİL ÇEKİRDEĞİ: FAZ TAŞIYAN ALAN DEĞERLENDİRİCİSİ
========================================================
Bu modülde ÇÖZÜCÜ YOKTUR. Yalnızca merdiven alanının

    S_v(t)  = − Σ_q A_q sin(ω_q t + φ_q)
    S_v'(t) = − Σ_q A_q ω_q cos(ω_q t + φ_q)

değerlendiricisi ve onun çok-süreçli sarmalayıcıları vardır.
`164_insa.S_ve_Sp`'nin BİREBİR aynı blok yapısı (blok=800 çizgi,
npt=10000 nokta), tek farkı `arg`'a φ_q'nun eklenmesidir. φ ≡ 0 iken
fonksiyon 164'ünkiyle bit-bit AYNI sonucu verir (176b bunu sınar).

Çözüm tarafında hiçbir şey kopyalanmaz: `164_insa.coz_sadakatli`
(ızgara braketi + SIRALI İLK-KÖK + korumalı Newton) AYNEN çağrılır;
yalnız onun kullandığı `S_par` / `SSp_par` isimleri buradaki faz taşıyan
sürümlerle değiştirilir (176b'de, koşudan önce, açıkça).
"""
import numpy as np

BLOK = 800          # 164_insa.BLOK ile aynı
NPT = 10000         # 164_insa.NPT ile aynı
NWORK = 7           # 164_insa.NWORK ile aynı

_GV = {}


def _init_v(om, A, phi):
    _GV["om"] = om
    _GV["A"] = A
    _GV["phi"] = phi


def S_ve_Sp_fazli(z, om, A, phi, deriv=True, blok=BLOK, npt=NPT):
    """S = −Σ A_q sin(ω_q z + φ_q),  S' = −Σ A_q ω_q cos(ω_q z + φ_q)."""
    S = np.zeros_like(z)
    Sp = np.zeros_like(z) if deriv else None
    for b0 in range(0, len(om), blok):
        w = om[b0:b0 + blok]
        aa = A[b0:b0 + blok]
        ph = phi[b0:b0 + blok]
        aw = aa * w if deriv else None
        for s0 in range(0, len(z), npt):
            sl = slice(s0, min(s0 + npt, len(z)))
            arg = np.outer(z[sl], w)
            arg += ph                       # ω_q z + φ_q  (satıra yayım)
            S[sl] += -np.sin(arg) @ aa
            if deriv:
                Sp[sl] += -np.cos(arg) @ aw
            del arg
    return (S, Sp) if deriv else (S, None)


def _work_Sv(zc):
    return S_ve_Sp_fazli(zc, _GV["om"], _GV["A"], _GV["phi"], deriv=False)[0]


def _work_SPv(zc):
    return S_ve_Sp_fazli(zc, _GV["om"], _GV["A"], _GV["phi"], deriv=True)


def _par_v(pool, fn, z, nw=NWORK):
    if len(z) == 0:
        return None
    parts = np.array_split(z, min(nw, max(1, len(z) // 200 + 1)))
    return pool.map(fn, parts)


def S_par_v(pool, z):
    return np.concatenate(_par_v(pool, _work_Sv, z))


def SSp_par_v(pool, z):
    r = _par_v(pool, _work_SPv, z)
    return (np.concatenate([x[0] for x in r]),
            np.concatenate([x[1] for x in r]))
