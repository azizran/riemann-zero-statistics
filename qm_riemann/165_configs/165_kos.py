"""
165 — T1/T3 SÜRÜCÜSÜ: bant bant ÖLÇÜM (referans) ↔ KANAL-ETİKETLİ
belirlenimli dörtlü toplam (pencere dönüşümü κ TAM, budama YOK).

Kullanım:  165_kos.py <veri> [taban] [kaynaklar] [tau_c listesi]
   ör:  165_kos.py Hkeskin 0.40 olculen,nominal 0.95,0.80,0.70,0.60,0.50

Çıktı: scratchpad/165/K_<veri>_t<taban>.json
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "165_configs"))
sys.path.insert(0, str(QM / "163_configs"))
K = importlib.import_module("165_cekirdek")
C163 = importlib.import_module("163_cekirdek")

TWO_PI = 2 * np.pi
SCR = K.SCR


def jkerr(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if len(v) < 4:
        return float("nan")
    return float(np.sqrt((len(v) - 1) / len(v) * np.sum((v - v.mean()) ** 2)))


def kos(veri, taban=0.40, kaynaklar=("olculen", "nominal"),
        tau_cs=(0.95,), bantlar=None, kmax=3, tau_max=0.95):
    t0 = time.time()
    bantlar = bantlar or K.IZGARA_T1
    Y = K.gaz(veri, taban, 4000)
    Mo = K.Model165(veri, taban, 4000, tau_max, "olculen", Y=Y)
    print(f"=== 165 {veri} taban={taban} tau_max={tau_max} ===")
    print(f"  N={len(Y.m0)} L={Y.L:.4f} T={Mo.T:.1f} dres={Mo.dres:.3e} "
          f"s̄={Mo.sbar:.1f} çizgi={len(Mo.M['q'])} qlim={Mo.qlim} "
          f"pencere={Mo.M['pen']}", flush=True)
    print(f"  KOHERANS: σ_Ĉ={Mo.sigC:.4f} → τ_c(faz)=1/(2πσ_Ĉ)="
          f"{Mo.tau_c_faz:.4f} ; σ_ds={Mo.sigds:.4f} → τ_c(gap)=1/(πσ_ds)="
          f"{Mo.tau_c_gap:.4f}", flush=True)

    ban = C163.bant_adaylari(Y, bantlar)
    qidx = {int(q): i for i, q in enumerate(Mo.M["q"])}
    Wall, bas = [], []
    for bi, b in enumerate(ban):
        bas.append(len(Wall))
        for r in b["cizgi"]:
            Wall.append(r["w"])
            Wall.append(r["w"] + r["gap"] / 2)
    Wall = np.array(Wall)
    print(f"  {len(ban)} bant, {len(Wall)//2} çizgi ({len(Wall)} frekans)",
          flush=True)

    tm = time.time()
    olc = [C163.olc_cizgi(Y, float(w), kmax=kmax) for w in Wall]
    print(f"  ölçüm {time.time()-tm:.0f}s", flush=True)

    sonuc = {}
    for kaynak in kaynaklar:
        for tc in tau_cs:
            etiket = f"{kaynak}_tc{tc}"
            tm = time.time()
            Mo.sec(kaynak, tc).alanlar(kmax=kmax)
            A = Mo.artik
            print(f"  [{etiket}] {A['nline']} çizgi, alanlar "
                  f"{time.time()-tm:.0f}s | Var(E) mod {A['varE_mod']:.5f} "
                  f"ölç {A['varE_olc']:.5f} kor {A['korE']:.3f} g_E "
                  f"{A['gE']:.3f} | Var(X) mod {A['varX_mod']:.5f} ölç "
                  f"{A['varX_olc']:.5f} kor {A['korX']:.3f} g_X "
                  f"{A['gX']:.3f} | g={A['gcal']:.4f}", flush=True)
            tm = time.time()
            J = Mo.ongor(Wall, kmax=kmax)
            print(f"  [{etiket}] öngörü {time.time()-tm:.0f}s", flush=True)

            sat = []
            for bi, b in enumerate(ban):
                recs = []
                for li, r in enumerate(b["cizgi"]):
                    d = dict(q=r["q"], w=r["w"], tau=r["tau"],
                             grup=r["grup"], gap=r["gap"], gp=r["gp"])
                    k0 = bas[bi] + 2 * li
                    for et, off in (("on", 0), ("off", 1)):
                        o = olc[k0 + off]
                        d[f"A_{et}"] = o["A"]
                        d[f"pow_{et}"] = o["pow"]
                        hQ, rm = o["h"], o["rho_ort"]
                        for kk in range(kmax + 1):
                            d[f"s{kk}_{et}"] = o[f"s{kk}"]
                            d[f"u{kk}_{et}"] = o[f"u{kk}"]
                            sp, up = K.s_den_J(hQ, J[kk][k0 + off], rm)
                            d[f"s{kk}p_{et}"] = sp
                            d[f"u{kk}p_{et}"] = up
                    iQ = qidx[int(r["q"])]
                    Jc1, JA, JBC, Jcak = Mo.kanal1(iQ)
                    Jc2, nc2 = Mo.kule_cozumleri(iQ)
                    hQ, rm = olc[k0]["h"], olc[k0]["rho_ort"]
                    d["s2_C1"] = K.s_den_J(hQ, Jc1, rm)[0]
                    d["s2_C1_Ta"] = K.s_den_J(hQ, JA, rm)[0]
                    d["s2_C1_BC"] = K.s_den_J(hQ, JBC, rm)[0]
                    d["s2_C2"] = K.s_den_J(hQ, Jc2, rm)[0]
                    d["n_C2"] = nc2
                    d["s2_capla"] = K.c1_capla(
                        r["tau"], Mo.M["tau"], Mo.M["a"], taban, tc)
                    recs.append(d)
                sat.append(dict(tau=b["tau"], lo=b["lo"], hi=b["hi"],
                                N=b["N"], kul=b["kul"], cizgi=recs))
            sonuc[etiket] = dict(bant=birlestir(sat, kmax, Mo.gcal),
                                 artik=dict(A), kaynak=kaynak, tau_c=tc)
            yaz_ozet(veri, etiket, sonuc[etiket]["bant"], Mo.gcal)

    out = dict(veri=veri, taban=taban, tau_max=tau_max, L=float(Y.L),
               N=int(len(Y.m0)), T=Mo.T, dres=Mo.dres, sbar=Mo.sbar,
               qlim=int(Mo.qlim), nline=int(len(Mo.M["q"])),
               sigC=Mo.sigC, sigds=Mo.sigds, tau_c_faz=Mo.tau_c_faz,
               tau_c_gap=Mo.tau_c_gap, kosum=sonuc, sure_s=time.time() - t0)
    SCR.mkdir(parents=True, exist_ok=True)
    p = SCR / f"K_{veri}_t{taban}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"-> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)
    return out


def birlestir(sat, kmax=3, gcal=1.0, njack=8):
    """160'ın bant birleştirmesi (birebir): Σ(pN A_on^p v_on
       − pO A_off^p v_off)/Σ(pN−pO)."""
    out = []
    for b in sat:
        L = b["cizgi"]
        if not L:
            out.append(dict(tau=b["tau"], olculdu=False))
            continue
        pN = np.array([r["pow_on"] for r in L])
        pO = np.array([r["pow_off"] for r in L])
        Ao = np.array([r["A_on"] for r in L])
        Af = np.array([r["A_off"] for r in L])
        tv = np.array([r["tau"] for r in L])
        grp = np.array([r["grup"] for r in L])
        u = pN - pO
        payda = float(u.sum())

        def AG(key, p=1, msk=None):
            m = slice(None) if msk is None else msk
            d = payda if msk is None else float(u[msk].sum())
            vN = np.array([r[f"{key}_on"] for r in L])
            vO = np.array([r[f"{key}_off"] for r in L])
            return float((pN[m] * Ao[m] ** p * vN[m]
                          - pO[m] * Af[m] ** p * vO[m]).sum() / d)

        def AC(key, p=2, msk=None):
            m = slice(None) if msk is None else msk
            d = payda if msk is None else float(u[msk].sum())
            v = np.array([r[key] for r in L])
            return float((pN[m] * Ao[m] ** p * v[m]).sum() / d)

        rec = dict(tau=b["tau"], lo=b["lo"], hi=b["hi"], N=b["N"],
                   kul=b["kul"], olculdu=True, gcal=gcal,
                   tau_eff=float((tv * u).sum() / payda))
        for k in range(kmax + 1):
            rec[f"A{k}s{k}"] = AG(f"s{k}", p=k)
            rec[f"A{k}s{k}p"] = AG(f"s{k}p", p=k)
            rec[f"A{k}u{k}"] = AG(f"u{k}", p=k)
            rec[f"A{k}u{k}p"] = AG(f"u{k}p", p=k)
        rec["A2s2_C1"] = AC("s2_C1")
        rec["A2s2_C1_Ta"] = AC("s2_C1_Ta")
        rec["A2s2_C1_BC"] = AC("s2_C1_BC")
        rec["A2s2_C2"] = AC("s2_C2")
        rec["A2s2_capla"] = AC("s2_capla")
        rec["A2s2_C3"] = rec["A2s2p"] - rec["A2s2_C1"] - rec["A2s2_C2"]
        rec["A2s2p_cal"] = rec["A2s2p"] * gcal
        rec["A2s2_C1_cal"] = rec["A2s2_C1"] * gcal
        rec["A2s2_C3_cal"] = rec["A2s2_C3"] * gcal
        # İKİNCİ-MOMENT NORMALİZASYONU: aynı dörtlü toplamın ρ-kanalı
        # (u2 = ⟨ρX̃²⟩/⟨ρ⟩) ölçüleni tutacak şekilde ölçeklenir.
        # (163'ün V4'ü "u0_pred ≡ 1" ile aynı ailedendir; üçüncü
        # momentten HİÇBİR girdi yoktur.)
        nrm = (rec["A2u2"] / rec["A2u2p"]) if rec["A2u2p"] else float("nan")
        rec["norm_u2"] = nrm
        rec["A2s2p_nrm"] = rec["A2s2p"] * nrm
        rec["A2s2_C1_nrm"] = rec["A2s2_C1"] * nrm
        rec["A2s2_C3_nrm"] = rec["A2s2_C3"] * nrm
        jm, jp, jc1, jc3 = [], [], [], []
        for g in range(njack):
            msk = grp != g
            if not msk.any() or abs(float(u[msk].sum())) < 1e-300:
                continue
            jm.append(AG("s2", p=2, msk=msk))
            jp.append(AG("s2p", p=2, msk=msk))
            jc1.append(AC("s2_C1", msk=msk))
            jc3.append(AG("s2p", p=2, msk=msk) - AC("s2_C1", msk=msk)
                       - AC("s2_C2", msk=msk))
        rec["sA2s2_jk"] = jkerr(jm)
        rec["sA2s2p_jk"] = jkerr(jp)
        rec["sA2s2_C1_jk"] = jkerr(jc1)
        rec["sA2s2_C3_jk"] = jkerr(jc3)
        out.append(rec)
    return out


def yaz_ozet(veri, etiket, bant, gcal):
    print(f"  --- {veri} / {etiket}  (g_alan={gcal:.4f}) ---")
    print("   τ_eff   A²s2 ölç   A²s2 HAM  ham/ölç  n_u2   A²s2 NRM  nrm/ölç"
          "   Ç1/tot Ç3/tot  Ç2/tot   u0p  A²u2ölç A²u2öng")
    for b in bant:
        if not b.get("olculdu"):
            continue
        m, p = b["A2s2"], b["A2s2p"]
        tot = p if abs(p) > 1e-300 else float("nan")
        pn = b["A2s2p_nrm"]
        print(f"   {b['tau_eff']:.4f} {m:+.6f} {p:+.6f} "
              f"{p/m if m else float('nan'):+7.3f} {b['norm_u2']:6.3f} "
              f"{pn:+.6f} {pn/m if m else float('nan'):+7.3f} "
              f"{b['A2s2_C1']/tot:+6.3f} {b['A2s2_C3']/tot:+6.3f} "
              f"{b['A2s2_C2']/tot:+8.1e} {b['A0u0p']:+.3f} "
              f"{b['A2u2']:+.4f} {b['A2u2p']:+.4f}", flush=True)


if __name__ == "__main__":
    veri = sys.argv[1]
    taban = float(sys.argv[2]) if len(sys.argv) > 2 else 0.40
    kk = (sys.argv[3].split(",") if len(sys.argv) > 3
          else ["olculen", "nominal"])
    tcs = ([float(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4
           else [0.95])
    kos(veri, taban, tuple(kk), tuple(tcs))
