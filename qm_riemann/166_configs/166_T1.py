"""
166 — T1: KALİB(τ_bant)'ın HAM TABLOSU (u2-sabitlemeden ÖNCE)
=============================================================
165'in T1 zincirini AYNEN koşar (165_cekirdek.Model165 + 163.olc_cizgi +
163.bant_adaylari; hiçbir ölçüm parçası kopyalanmaz) ama ÇİZGİ BAZINDA
kayıt tutar; böylece 165'in yalnız bant-toplamı saklanan büyüklükleri
(özellikle `norm_u2` ve ölçüm/çıplak-öngörü oranı) için de JACKKNIFE
hatası çıkarılabilir.

TANIMLAR (hepsi bant-birleştirmesinden SONRA, 160'ın formülüyle):
    Ms2  = A²s2 (ÖLÇÜLEN üçüncü moment)        Ps2  = A²s2 ÇIPLAK ÖNGÖRÜ
    Mu2  = A²u2 (ÖLÇÜLEN ikinci moment)        Pu2  = A²u2 ÇIPLAK ÖNGÖRÜ
    KALİB_s2 ≡ Ms2/Ps2   — üçüncü momentin İSTEDİĞİ kalibrasyon (hedef)
    KALİB_u2 ≡ Mu2/Pu2   — 165'in KULLANDIĞI kalibrasyon (`norm_u2`)
165'in hükmü `KALİB_u2`yi öngörüye çarpıp `KALİB_u2/KALİB_s2` oranına
bakar (raporun `nrm/ölç` sütunu).

Jackknife: bant içi 8-grup (`grup` = 160'ın round-robin'i); silme-1
tahmincisi her BİLEŞİK büyüklük için yeniden hesaplanır (oran
tahmincilerinde doğrusal hata yayılımı YAPILMAZ).

Kullanım:  166_T1.py <gaz> [taban] [kaynak] [tau_c]
Çıktı:     scratchpad/166/T1_<gaz>_<kaynak>_tc<tau_c>.json
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("166_configs", "165_configs", "163_configs", "160_configs",
           "159_configs", "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
K = importlib.import_module("165_cekirdek")
C163 = importlib.import_module("163_cekirdek")

TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/166")
KMAX = 3


def _jk(v):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if len(v) < 4:
        return float("nan")
    return float(np.sqrt((len(v) - 1) / len(v) * np.sum((v - v.mean()) ** 2)))


def bant_agg(L, msk=None):
    """160'ın bant birleştirmesi (birebir) — çizgi kayıtlarından."""
    if msk is None:
        msk = np.ones(len(L), bool)
    pN = np.array([r["pow_on"] for r in L])[msk]
    pO = np.array([r["pow_off"] for r in L])[msk]
    Ao = np.array([r["A_on"] for r in L])[msk]
    Af = np.array([r["A_off"] for r in L])[msk]
    tv = np.array([r["tau"] for r in L])[msk]
    d = float((pN - pO).sum())
    out = {"tau_eff": float((tv * (pN - pO)).sum() / d)}

    def AG(key, p):
        vN = np.array([r[f"{key}_on"] for r in L])[msk]
        vO = np.array([r[f"{key}_off"] for r in L])[msk]
        return float((pN * Ao ** p * vN - pO * Af ** p * vO).sum() / d)

    for k in range(KMAX + 1):
        out[f"Ms{k}"] = AG(f"s{k}", k)
        out[f"Ps{k}"] = AG(f"s{k}p", k)
        out[f"Mu{k}"] = AG(f"u{k}", k)
        out[f"Pu{k}"] = AG(f"u{k}p", k)
    for key, nm in (("s2_C1", "C1"), ("s2_C2", "C2"), ("s2_capla", "capla")):
        v = np.array([r[key] for r in L])[msk]
        out[nm] = float((pN * Ao ** 2 * v).sum() / d)
    out["C3"] = out["Ps2"] - out["C1"] - out["C2"]
    # ρ (144 estimatörü) AYNI çizgilerde, iki tabanda
    for et in ("r52", "r40"):
        if f"{et}_on" in L[0]:
            on = np.array([r[f"{et}_on"] for r in L])[msk]
            of = np.array([r[f"{et}_off"] for r in L])[msk]
            gp = np.array([r["gp"] for r in L])[msk]
            out[f"rho_{et}"] = float((on - of).sum() / gp.sum())
    out["KALIB_s2"] = out["Ms2"] / out["Ps2"] if out["Ps2"] else float("nan")
    out["KALIB_u2"] = out["Mu2"] / out["Pu2"] if out["Pu2"] else float("nan")
    return out


def kos(veri, taban=0.40, kaynak="olculen", tau_c=0.95, njack=8):
    t0 = time.time()
    Y = K.gaz(veri, taban, 4000)
    Mo = K.Model165(veri, taban, 4000, 0.95, kaynak, Y=Y, tau_c=tau_c)
    Mo.sec(kaynak, tau_c).alanlar(kmax=KMAX)
    A = Mo.artik
    print(f"=== 166-T1 {veri} taban={taban} kaynak={kaynak} τ_c={tau_c} ===")
    print(f"  N={len(Y.m0)} L={Y.L:.5f} σ_Ĉ={Mo.sigC:.5f} σ_ds={Mo.sigds:.5f}"
          f"  g_E={A['gE']:.4f} g_X={A['gX']:.4f} g={A['gcal']:.4f} "
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

    # ρ (144 estimatörü) AYNI frekanslarda, iki η konvansiyonunda
    rho_proj = {}
    try:
        R = importlib.import_module("166_rho")
        z = importlib.import_module("155_kos").veri_yukle(veri)
        mid = 0.5 * (z[:-1] + z[1:])
        K155 = importlib.import_module("155_cekirdek")
        for et, (tb, cp) in (("r52", (0.52, 720)), ("r40", (0.40, 4000))):
            Ç = K155.eta_onbellek(z, veri, tb, cp)
            rho_proj[et] = np.abs(R.proj(Ç["eta"], mid, Wall)) ** 2
        print(f"  ρ projeksiyonları hazır", flush=True)
    except Exception as e:                                # pragma: no cover
        print(f"  [ρ atlandı: {e}]", flush=True)

    qidx = {int(q): i for i, q in enumerate(Mo.M["q"])}
    bantlar = []
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
                for et2, P in rho_proj.items():
                    d[f"{et2}_{et}"] = float(P[k0 + off])
            iQ = qidx[int(r["q"])]
            Jc1 = Mo.kanal1(iQ)[0]
            Jc2 = Mo.kule_cozumleri(iQ)[0]
            hQ, rm = olc[k0]["h"], olc[k0]["rho_ort"]
            d["s2_C1"] = K.s_den_J(hQ, Jc1, rm)[0]
            d["s2_C2"] = K.s_den_J(hQ, Jc2, rm)[0]
            d["s2_capla"] = K.c1_capla(r["tau"], Mo.M["tau"], Mo.M["a"],
                                       taban, tau_c)
            L.append(d)
        if not L:
            bantlar.append(dict(tau=b["tau"], lo=b["lo"], hi=b["hi"],
                                olculdu=False))
            continue
        rec = dict(tau=b["tau"], lo=b["lo"], hi=b["hi"], N=b["N"],
                   kul=b["kul"], olculdu=True)
        rec.update(bant_agg(L))
        grp = np.array([r["grup"] for r in L])
        jk = {k: [] for k in ("Ms2", "Ps2", "Mu2", "Pu2", "KALIB_s2",
                              "KALIB_u2", "C1", "C3", "rho_r52", "rho_r40")}
        for g in range(njack):
            m = grp != g
            if not m.any():
                continue
            try:
                a = bant_agg(L, m)
            except (ZeroDivisionError, FloatingPointError):
                continue
            for k in jk:
                if k in a:
                    jk[k].append(a[k])
        for k, v in jk.items():
            if v:
                rec[f"s{k}"] = _jk(v)
        rec["cizgi"] = L
        bantlar.append(rec)

    out = dict(veri=veri, taban=taban, kaynak=kaynak, tau_c=tau_c,
               L=float(Y.L), N=int(len(Y.m0)), T=Mo.T, dres=Mo.dres,
               sbar=Mo.sbar, sigC=Mo.sigC, sigds=Mo.sigds,
               tau_c_faz=Mo.tau_c_faz, tau_c_gap=Mo.tau_c_gap,
               artik=dict(A), bant=bantlar, sure_s=time.time() - t0)
    SCR.mkdir(parents=True, exist_ok=True)
    p = SCR / f"T1_{veri}_{kaynak}_tc{tau_c}.json"
    p.write_text(json.dumps(out, indent=1))

    print("\n  τ_eff   Ms2        Ps2        KALİB_s2±jk     "
          "Mu2      Pu2       KALİB_u2±jk    ρ(0.52)  ρ(0.40)")
    for b in bantlar:
        if not b.get("olculdu"):
            continue
        print(f"  {b['tau_eff']:.4f} {b['Ms2']:+.6f} {b['Ps2']:+.6f} "
              f"{b['KALIB_s2']:+8.4f}±{b.get('sKALIB_s2', float('nan')):.4f} "
              f"{b['Mu2']:+.4f} {b['Pu2']:+8.4f} "
              f"{b['KALIB_u2']:.4f}±{b.get('sKALIB_u2', float('nan')):.4f} "
              f"{b.get('rho_r52', float('nan')):+.4f} "
              f"{b.get('rho_r40', float('nan')):+.4f}", flush=True)
    print(f"-> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)
    return out


if __name__ == "__main__":
    kos(sys.argv[1],
        float(sys.argv[2]) if len(sys.argv) > 2 else 0.40,
        sys.argv[3] if len(sys.argv) > 3 else "olculen",
        float(sys.argv[4]) if len(sys.argv) > 4 else 0.95)
