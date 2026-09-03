"""
167 — TEK GAZDA `c` ÖLÇÜMÜ (166'nın T1 zinciri + bacak sönümleri, yalın)
=======================================================================
166_T1'in bant birleştirmesi (`166_T1.bant_agg`) ve 165'in ölçüm/öngörü
zinciri (`165_cekirdek.Model165`, `163_cekirdek.olc_cizgi/bant_adaylari`)
AYNEN import edilir — hiçbir ölçüm parçası kopyalanmaz. Eklenen:

  * W_amp(τ)=⟨cos πτ·ds⟩, W_X(τ)=⟨e^{−2πiτX̃}⟩, W_pos(τ)=⟨e^{2πiτĈ}⟩
    aday çizgilerde TAM (166_bacak.karakteristik ile aynı tanım),
    bantta `gp`-ağırlıklı ortalama (166_bacak §5 ile birebir);
  * **c_bant = KALİB_u2 / (W_amp·W_X)** ve bant-içi 8-grup jackknife;
  * SADAKAT DENETİMİ (164 ölçütü) AYNI aday çizgilerde:
        R_bant = sqrt[(Σ|c_on|² − Σ|c_off|²)/Σ b_nom²],
        c_W = 2⟨(ds−⟨ds⟩)e^{−iWm}⟩ , b_nom = 2·λ·a_q·w_q·sin(πτ_q)
  * isteğe bağlı DÜZGÜN TARAK kontrolü (165_tarak ile aynı) → Ç4 payı.

Kullanım: 167_olcum.py <gaz> [taban] [tau_c] [kule=0|1] [duz=0|1]
Çıktı:    scratchpad/167/C_<gaz>.json
"""
import importlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("167_configs", "166_configs", "165_configs", "163_configs",
           "160_configs", "159_configs", "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
K = importlib.import_module("165_cekirdek")
C163 = importlib.import_module("163_cekirdek")
T166 = importlib.import_module("166_T1")
B166 = importlib.import_module("166_bacak")
ORT = importlib.import_module("167_ortak")

TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/167")
KMAX = 3
LO_MIN = 0.52          # hüküm penceresi (166 ile aynı)

# 165'in erfc pencereleri: 167'nin gazlarını ekle
K.PENCERE.update(ORT.pencere_dict())


def _proj(f, mid, Wler, blok=20000, fblok=256):
    """2⟨f e^{−iW m}⟩ — 163.tayf ile aynı cebir, tek alan."""
    N = len(mid)
    Wler = np.asarray(Wler, float)
    acc = np.zeros(len(Wler), dtype=complex)
    for f0 in range(0, len(Wler), fblok):
        fs = slice(f0, min(f0 + fblok, len(Wler)))
        Wc = Wler[fs]
        for s in range(0, N, blok):
            sl = slice(s, min(s + blok, N))
            arg = np.outer(mid[sl], Wc)
            acc[fs] += f[sl] @ np.cos(arg) - 1j * (f[sl] @ np.sin(arg))
            del arg
    return 2 * acc / N


def kos(veri, taban=0.40, tau_c=0.95, kule=1, duz=0, njack=8):
    t0 = time.time()
    kun = ORT.KUNYE.get(veri, dict(lam=1.0, tau_ust=1.00, pen=None,
                                   aile="?", gercek=False))
    lam = kun.get("lam") or 1.0
    tau_ust = kun.get("tau_ust") or 1.00
    pen = kun.get("pen")
    Y = K.gaz(veri, taban, 4000)
    Mo = K.Model165(veri, taban, 4000, 0.95, "olculen", Y=Y, tau_c=tau_c)
    Mo.sec("olculen", tau_c).alanlar(kmax=KMAX)
    A = Mo.artik
    print(f"=== 167-ÖLÇÜM {veri} (aile={kun['aile']} λ={lam} τ_ust={tau_ust} "
          f"pen={pen}) taban={taban} τ_c={tau_c} ===", flush=True)
    print(f"  N={len(Y.m0)} L={Y.L:.5f} T={Mo.T:.1f} dres={Mo.dres:.4e} "
          f"σ_Ĉ={Mo.sigC:.5f} σ_ds={Mo.sigds:.5f} σ_X̃={np.std(Y.Xtil0):.5f}",
          flush=True)
    print(f"  g_E={A['gE']:.4f} g_X={A['gX']:.4f} g={A['gcal']:.4f} "
          f"Var(E_mod)/Var(η)={A['varE_mod']/A['varE_olc']:.3f} "
          f"Var(X_mod)/Var(X̃)={A['varX_mod']/A['varX_olc']:.3f} "
          f"çizgi={A['nline']}", flush=True)

    ban = C163.bant_adaylari(Y, K.IZGARA_T1)
    Wall, bas = [], []
    for b in ban:
        bas.append(len(Wall))
        for r in b["cizgi"]:
            Wall.append(r["w"])
            Wall.append(r["w"] + r["gap"] / 2)
    Wall = np.array(Wall)
    print(f"  {len(ban)} bant, {len(Wall)//2} çizgi ({len(Wall)} frekans)",
          flush=True)

    tm = time.time()
    olc = [C163.olc_cizgi(Y, float(w), kmax=KMAX) for w in Wall]
    print(f"  ölçüm {time.time()-tm:.0f}s", flush=True)
    tm = time.time()
    J = Mo.ongor(Wall, kmax=KMAX)
    print(f"  öngörü {time.time()-tm:.0f}s", flush=True)

    # --- DÜZGÜN TARAK (Ç4 payı) ---------------------------------------
    Jd = None
    if duz:
        tm = time.time()
        gbar = TWO_PI / Y.L
        su = Mo.s[0] + gbar * np.arange(len(Mo.s))
        w_ = Mo.M["w"][Mo.msk]
        Ed = K.sentez(su, w_, Mo.hp[Mo.msk])
        Xd = K.sentez(su, w_, Mo.y[Mo.msk])
        Xd = Xd - Xd.mean()
        Gd = [Ed, Ed * Xd, Ed * Xd * Xd]
        Jd = [x / 2.0 for x in K.tayf_s(su, Gd, Wall)]
        print(f"  düzgün tarak {time.time()-tm:.0f}s "
              f"(Var Ed={np.var(Ed):.5f} Var Xd={np.var(Xd):.5f})", flush=True)

    # --- SADAKAT: ds projeksiyonları + nominal b ------------------------
    tm = time.time()
    dsm = Y.ds - Y.ds.mean()
    cds = _proj(dsm, Y.mid, Wall)
    print(f"  ds projeksiyonu {time.time()-tm:.0f}s", flush=True)

    def b_nom(q, mult, tau):
        if tau > tau_ust + 1e-12:
            return 0.0
        a = lam / (np.pi * mult * math.sqrt(q))
        if pen is not None:
            a *= 0.5 * math.erfc((tau - pen[0]) / pen[1])
        return 2 * a * math.sin(np.pi * tau)

    # --- W çarpanları aday çizgilerde (TAM karakteristik fonksiyon) ----
    tm = time.time()
    tauc = np.array([r["tau"] for b in ban for r in b["cizgi"]])
    W_amp_c = B166.karakteristik(Y.ds, np.pi * tauc).real
    W_X_c = B166.karakteristik(Y.Xtil0, -TWO_PI * tauc)
    C_ = np.cumsum(Y.Xtil0)
    nn = np.arange(len(C_), dtype=float)
    nn = (nn - nn.mean()) / (nn[-1] if len(nn) > 1 else 1.0)
    V = np.vstack([np.ones_like(nn), nn, nn * nn, nn ** 3]).T
    cc, *_ = np.linalg.lstsq(V, C_, rcond=None)
    Chat = C_ - V @ cc
    W_pos_c = B166.karakteristik(Chat, TWO_PI * tauc)
    print(f"  W çarpanları {time.time()-tm:.0f}s", flush=True)

    qidx = {int(q): i for i, q in enumerate(Mo.M["q"])}
    bantlar = []
    ci = 0
    for bi, b in enumerate(ban):
        L = []
        for li, r in enumerate(b["cizgi"]):
            d = dict(q=r["q"], w=r["w"], tau=r["tau"], grup=r["grup"],
                     gap=r["gap"], gp=r["gp"])
            k0 = bas[bi] + 2 * li
            for et, off in (("on", 0), ("off", 1)):
                o = olc[k0 + off]
                d[f"A_{et}"] = o["A"]
                d[f"pow_{et}"] = o["pow"]
                for kk in range(KMAX + 1):
                    d[f"s{kk}_{et}"] = o[f"s{kk}"]
                    d[f"u{kk}_{et}"] = o[f"u{kk}"]
                    sp, up = K.s_den_J(o["h"], J[kk][k0 + off], o["rho_ort"])
                    d[f"s{kk}p_{et}"] = sp
                    d[f"u{kk}p_{et}"] = up
                d[f"cds_{et}"] = float(abs(cds[k0 + off]) ** 2)
            if Jd is not None:
                for et, off in (("on", 0), ("off", 1)):
                    o = olc[k0 + off]
                    sp, up = K.s_den_J(o["h"], Jd[2][k0 + off], o["rho_ort"])
                    d[f"s2d_{et}"] = sp
                    d[f"u2d_{et}"] = up
            iQ = qidx[int(r["q"])]
            d["s2_C1"] = K.s_den_J(olc[k0]["h"], Mo.kanal1(iQ)[0],
                                   olc[k0]["rho_ort"])[0]
            d["u2_C1"] = K.s_den_J(olc[k0]["h"], Mo.kanal1(iQ)[0],
                                   olc[k0]["rho_ort"])[1]
            d["s2_C2"] = (K.s_den_J(olc[k0]["h"], Mo.kule_cozumleri(iQ)[0],
                                    olc[k0]["rho_ort"])[0] if kule else 0.0)
            d["s2_capla"] = K.c1_capla(r["tau"], Mo.M["tau"], Mo.M["a"],
                                       taban, tau_c)
            d["bnom"] = b_nom(r["q"], Mo.M["mult"][iQ], r["tau"])
            d["W_amp"] = float(W_amp_c[ci])
            d["W_X"] = float(W_X_c[ci].real)
            d["W_pos"] = float(W_pos_c[ci].real)
            ci += 1
            L.append(d)
        if not L:
            bantlar.append(dict(tau=b["tau"], lo=b["lo"], hi=b["hi"],
                                olculdu=False))
            continue
        rec = dict(tau=b["tau"], lo=b["lo"], hi=b["hi"], N=b["N"],
                   kul=b["kul"], olculdu=True)
        rec.update(_agg(L))
        grp = np.array([r["grup"] for r in L])
        jk = {k: [] for k in ("Ms2", "Ps2", "Mu2", "Pu2", "KALIB_s2",
                              "KALIB_u2", "c_u2", "c_s2", "C1", "C3")}
        for g in range(njack):
            m = grp != g
            if not m.any():
                continue
            try:
                a = _agg(L, m)
            except (ZeroDivisionError, FloatingPointError):
                continue
            for k in jk:
                if k in a:
                    jk[k].append(a[k])
        for k, v in jk.items():
            if v:
                rec[f"s{k}"] = T166._jk(v)
        rec["cizgi"] = L
        bantlar.append(rec)

    out = dict(veri=veri, kunye=kun, lam=lam, tau_ust=tau_ust, pen=pen,
               taban=taban, tau_c=tau_c, kule=int(kule), duz=int(duz),
               L=float(Y.L), N=int(len(Y.m0)), T=Mo.T, dres=Mo.dres,
               sbar=Mo.sbar, sigC=Mo.sigC, sigds=Mo.sigds,
               sigX=float(np.std(Y.Xtil0)), sigChat=float(np.std(Chat)),
               tau_c_faz=Mo.tau_c_faz, tau_c_gap=Mo.tau_c_gap,
               artik=dict(A), bant=bantlar, sure_s=time.time() - t0)
    # gaz-düzeyi c (hüküm bantlarında log-ortalama)
    H = [b for b in bantlar if b.get("olculdu") and b["lo"] >= LO_MIN - 1e-9
         and b["tau_eff"] < 0.85 and b["Ms2"] > 0 and b["R_bant"] >= 0.98]
    if H:
        lg = np.log([b["c_u2"] for b in H])
        out["c_gaz"] = float(np.exp(lg.mean()))
        out["c_gaz_sd"] = float(np.std(lg, ddof=1)) if len(H) > 1 else 0.0
        out["c_nbant"] = len(H)
        te = np.array([b["tau_eff"] for b in H])
        out["c_egim"] = float(np.polyfit(te, lg - lg.mean(), 1)[0]) \
            if len(H) > 2 else float("nan")
    SCR.mkdir(parents=True, exist_ok=True)
    p = SCR / f"C_{veri}.json"
    p.write_text(json.dumps(out, indent=1))

    print("\n  τ_eff  KALİB_u2±jk      W_amp   W_X    W_amp·W_X   "
          "**c_u2±jk**       R_bant  SNR    KALİB_s2")
    for b in bantlar:
        if not b.get("olculdu"):
            continue
        yz = "*" if (b["lo"] >= LO_MIN - 1e-9 and b["tau_eff"] < 0.85
                     and b["Ms2"] > 0 and b["R_bant"] >= 0.98) else " "
        print(f" {yz}{b['tau_eff']:.4f} {b['KALIB_u2']:.4f}±"
              f"{b.get('sKALIB_u2', float('nan')):.4f}  {b['W_amp']:.4f} "
              f"{b['W_X']:.4f} {b['W_amp']*b['W_X']:.4f}   "
              f"{b['c_u2']:.4f}±{b.get('sc_u2', float('nan')):.4f}   "
              f"{b['R_bant']:.4f} {b['SNR']:7.1f} {b['KALIB_s2']:+.4f}",
              flush=True)
    if H:
        print(f"\n  **c({veri}) = {out['c_gaz']:.4f}**  "
              f"(log-sd {out['c_gaz_sd']:.4f}, {out['c_nbant']} bant, "
              f"eğim {out['c_egim']:+.3f})", flush=True)
    print(f"-> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)
    return out


def _agg(L, msk=None):
    """166_T1.bant_agg + W çarpanları + sadakat + c."""
    out = T166.bant_agg(L, msk)
    if msk is None:
        msk = np.ones(len(L), bool)
    gp = np.array([r["gp"] for r in L])[msk]
    bn = np.array([r["bnom"] for r in L])[msk]
    on = np.array([r["cds_on"] for r in L])[msk]
    of = np.array([r["cds_off"] for r in L])[msk]
    pN = np.array([r["pow_on"] for r in L])[msk]
    pO = np.array([r["pow_off"] for r in L])[msk]
    gw = gp / gp.sum()
    for k in ("W_amp", "W_X", "W_pos"):
        out[k] = float((np.array([r[k] for r in L])[msk] * gw).sum())
    Bp = float(np.sum(bn ** 2))
    net = float(on.sum() - of.sum())
    out["R_bant"] = float(np.sqrt(max(net, 0.0) / Bp)) if Bp > 0 else 0.0
    out["SNR"] = float(pN.sum() / pO.sum()) if pO.sum() > 0 else float("inf")
    out["c_u2"] = out["KALIB_u2"] / (out["W_amp"] * out["W_X"])
    out["c_s2"] = out["KALIB_s2"] / (out["W_amp"] * out["W_X"])
    if "s2d_on" in L[0]:
        Ao = np.array([r["A_on"] for r in L])[msk]
        Af = np.array([r["A_off"] for r in L])[msk]
        d = float((pN - pO).sum())
        for nm, key in (("Ps2d", "s2d"), ("Pu2d", "u2d")):
            vN = np.array([r[f"{key}_on"] for r in L])[msk]
            vO = np.array([r[f"{key}_off"] for r in L])[msk]
            p = 2 if nm[1] == "s" or nm[1] == "u" else 0
            out[nm] = float((pN * Ao ** p * vN - pO * Af ** p * vO).sum() / d)
        out["C4"] = out["Ps2"] - out["Ps2d"]
    return out


if __name__ == "__main__":
    kos(sys.argv[1],
        float(sys.argv[2]) if len(sys.argv) > 2 else 0.40,
        float(sys.argv[3]) if len(sys.argv) > 3 else 0.95,
        int(sys.argv[4]) if len(sys.argv) > 4 else 1,
        int(sys.argv[5]) if len(sys.argv) > 5 else 0)
