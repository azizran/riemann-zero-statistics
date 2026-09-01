"""
162 — DENETİM
V1  BİT DÜZEYİ: `gercek` koşusu 160'ın kayıtlı JSON'unu birebir veriyor mu
V2  BORU HATTI KİMLİĞİ: V0 (faz karıştırmasız vekil) = gerçek mi (dizi düzeyi)
V3  İKİNCİ MOMENTLER: vekil, oto-/çapraz-kovaryansı ve çizgi güçlerini
    koruyor mu (otokovaryans profili + |c_q| bant bant)
V4  GAUSS'LUK: çarpıklık/basıklık gerçek ↔ vekil
V5  GRAM KONTROLÜ: Σ|c|²/2 > Var(η) bir çıkarım artefaktı mı, yoksa
    gerçeğin YIKICI GİRİŞİMİ mi? (bilinen rastgele-fazlı sinyalde sına)
V6  SABİT NOKTA: Newton izleri (koşu JSON'larından)
V7  RED-1: düz n-uzayı FFT faz karıştırması çizgiyi öldürüyor mu (ölç)
V8  RED-3: öz-tutarlılık olmadan çizgi gücü ne oluyor (ölç)
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
C = importlib.import_module("162_cekirdek")
sys.path.insert(0, str(C.QM / "160_configs"))
A160 = importlib.import_module("160_analiz")
SCR, SCR160 = C.SCR, A160.SCR
TWO_PI = 2 * np.pi
GAZ = ("son", "keskin", "A4")


def yukle(v, c, g="t1", t=0.4):
    p = SCR / f"S_{v}_t{t}_{g}_{c}.json"
    return json.load(open(p)) if p.exists() else None


def V1():
    print("\n" + "=" * 78)
    print("V1. BİT DÜZEYİ — 162/gercek  vs  160'ın kayıtlı JSON'u")
    print("=" * 78)
    kol = ["phi", "delta", "tau_eff", "absG", "absGc", "K1", "K2", "Kfull",
           "M1", "M0", "sPhi_jk", "sDelta_jk", "sK1_jk", "sKf_jk", "sAS_jk",
           "rho", "rho_norm", "absE1", "absEf", "E1_re", "E1_im", "Es_re",
           "Es_im", "kinematik", "delta_b", "sA2", "L"]
    top, tn = 0.0, 0
    for v in GAZ:
        for g in ("t1", "g158"):
            d = yukle(v, "gercek", g)
            p = SCR160 / f"S_{v}_t0.4_{g}.json"
            if d is None or not p.exists():
                continue
            r = json.load(open(p))
            ix = {b["tau"]: b for b in r["bantlar"] if b.get("olculdu")}
            mx, n, kt = 0.0, 0, None
            for b in A160.bantlar(d):
                o = ix.get(b["tau"])
                if o is None:
                    continue
                n += 1
                for k in kol + ["Xort", "sTeff_jk"]:
                    if k in b and k in o and np.isfinite(b[k]) and \
                            np.isfinite(o[k]) and abs(b[k] - o[k]) > mx:
                        mx, kt = abs(b[k] - o[k]), k
                for k in b["AS"]:
                    if abs(b["AS"][k] - o["AS"][k]) > mx:
                        mx, kt = abs(b["AS"][k] - o["AS"][k]), "AS." + k
                for k in b["KUM"]:
                    if abs(b["KUM"][k] - o["KUM"][k]) > mx:
                        mx, kt = abs(b["KUM"][k] - o["KUM"][k]), "KUM." + k
            top = max(top, mx); tn += n
            print(f"   {v:<8} {g:<5} {n:>3} ortak bant  maks fark {mx:.3e}"
                  f"   ({kt})")
    print(f"   → {tn} bant, BÜTÜN sütunlarda maks fark {top:.3e}")


def V2(veri="son"):
    print("\n" + "=" * 78)
    print("V2. BORU HATTI KİMLİĞİ — V0 (ψ≡0) gerçeğin ta kendisi mi?")
    print("=" * 78)
    T = C.Taban162(veri, 0.40, ayrinti=False)
    Y0 = T.vekil_yerel(None, "V0", ayrinti=False)
    for ad, a, b in (("η", T.Y.eta, Y0.eta), ("mid", T.Y.mid, Y0.mid),
                     ("X̃", T.Y.Xtil, Y0.Xtil),
                     ("X_tam", T.Y.X["tam"], Y0.X["tam"]),
                     ("Ĉ", T.Chat, Y0.Chat)):
        print(f"   {ad:<7} maks |fark| = {np.max(np.abs(a-b)):.3e}   "
              f"(rms sinyal {np.std(a):.3e})")
    print(f"   m-kimlik  m_n − [m_0+(2π/L)(u_n+Ĉ_n)] maks = {T.m_kimlik:.3e}")
    d0, dg = yukle(veri, "V0"), yukle(veri, "gercek")
    if d0 and dg:
        ix = {b["tau"]: b for b in A160.bantlar(d0)}
        mx = max(abs(b["delta"] - ix[b["tau"]]["delta"])
                 for b in A160.bantlar(dg) if b["tau"] in ix)
        print(f"   bant düzeyi maks |δ_V0 − δ_gerçek| = {mx:.3e}")
    return T


def V3(T, tohum=101):
    print("\n" + "=" * 78)
    print("V3. İKİNCİ MOMENTLER — vekil oto-/çapraz-kovaryansı koruyor mu?")
    print("=" * 78)
    Yv = T.vekil_yerel(tohum, "VG", ayrinti=False)
    Y = T.Y
    print("   otokovaryans C_xx(k)/C_xx(0) — η (indeks gecikmesi)")
    print("     k      gerçek      vekil")
    e0 = Y.eta - Y.eta.mean(); e1 = Yv.eta - Yv.eta.mean()
    for k in (1, 2, 3, 5, 10, 50):
        print(f"    {k:>3}   {np.mean(e0[:-k]*e0[k:])/e0.var():+9.5f}   "
              f"{np.mean(e1[:-k]*e1[k:])/e1.var():+9.5f}")
    print("   otokovaryans — Ĉ")
    c0 = T.Chat - T.Chat.mean(); c1 = Yv.Chat - Yv.Chat.mean()
    for k in (1, 2, 5, 20, 100):
        print(f"    {k:>3}   {np.mean(c0[:-k]*c0[k:])/c0.var():+9.5f}   "
              f"{np.mean(c1[:-k]*c1[k:])/c1.var():+9.5f}")
    print(f"   çapraz  ⟨ηĈ⟩/σ_ησ_C  gerçek {np.mean(e0*c0)/np.sqrt(e0.var()*c0.var()):+.5f}"
          f"   vekil {np.mean(e1*c1)/np.sqrt(e1.var()*c1.var()):+.5f}")
    ce, cc = C.cizgi_cikar(Yv.mid, [Yv.eta, Yv.Chat], T.w)
    tau = T.w / T.Y.L
    print("\n   ÇİZGİ GÜÇLERİ (|c_q|²/2, bant bant) gerçek → vekil")
    print("    τ bandı    n     η gerçek     η vekil    oran   "
          "Ĉ gerçek     Ĉ vekil    oran")
    for lo in np.arange(0.40, 0.86, 0.05):
        m = (tau > lo) & (tau <= lo + 0.05)
        if not m.any():
            continue
        pe = np.sum(np.abs(T.c_eta[m]) ** 2) / 2
        qe = np.sum(np.abs(ce[m]) ** 2) / 2
        pc = np.sum(np.abs(T.c_C[m]) ** 2) / 2
        qc = np.sum(np.abs(cc[m]) ** 2) / 2
        print(f"    {lo:.2f}-{lo+0.05:.2f} {int(m.sum()):5d}  {pe:.4e}  "
              f"{qe:.4e}  {qe/pe:5.2f}   {pc:.4e}  {qc:.4e}  {qc/pc:5.2f}")
    return Yv


def V4(T, Yv):
    print("\n" + "=" * 78)
    print("V4. GAUSS'LUK — çarpıklık / basıklık")
    print("=" * 78)
    print("   dizi      gerçek çarp/bas        vekil çarp/bas")
    for ad, a, b in (("η", T.Y.eta, Yv.eta), ("Ĉ", T.Chat, Yv.Chat),
                     ("dsΔ", T.Y.X["tam"], Yv.X["tam"])):
        def cb(x):
            x = x - x.mean()
            v = x.var()
            return np.mean(x ** 3) / v ** 1.5, np.mean(x ** 4) / v ** 2
        s0, k0 = cb(a); s1, k1 = cb(b)
        print(f"   {ad:<8}  {s0:+8.4f} / {k0:7.4f}      "
              f"{s1:+8.4f} / {k1:7.4f}")


def V5(T, tohum=901):
    print("\n" + "=" * 78)
    print("V5. GRAM KONTROLÜ — Σ|c|²/2 > Var(η) bir artefakt mı?")
    print("    Bilinen rastgele-fazlı sinyal kurulur (planlanan güç BİLİNİYOR),")
    print("    aynı çıkarım uygulanır. Geri gelen güç planlananla eşitse")
    print("    çıkarım YANSIZDIR ⇒ gerçekteki fazla, YIKICI GİRİŞİMDİR.")
    print("=" * 78)
    rng = np.random.default_rng(tohum)
    rot = np.exp(1j * rng.uniform(0, TWO_PI, len(T.w)))
    x = C.cizgi_kur(T.Y.mid, [T.c_eta * rot], T.w)[0]
    plan = float(np.sum(np.abs(T.c_eta) ** 2) / 2)
    c2 = C.cizgi_cikar(T.Y.mid, [x], T.w)[0]
    geri = float(np.sum(np.abs(c2) ** 2) / 2)
    print(f"   planlanan Σ|c|²/2 = {plan:.6f}   sinyalin Var = {np.var(x):.6f}"
          f"   (oran {np.var(x)/plan:.4f})")
    print(f"   çıkarımdan geri gelen Σ|c'|²/2 = {geri:.6f}  "
          f"(plana oran {geri/plan:.4f})")
    print(f"   GERÇEK: Σ|c_η|²/2 = {plan:.6f}  Var(η) = {np.var(T.Y.eta):.6f}"
          f"  (oran {plan/np.var(T.Y.eta):.4f})")
    print(f"   GERÇEK Ĉ: Σ|c_Ĉ|²/2 = {np.sum(np.abs(T.c_C)**2)/2:.6f}  "
          f"Var(Ĉ) = {np.var(T.Chat):.6f}  "
          f"(oran {float(np.sum(np.abs(T.c_C)**2)/2)/np.var(T.Chat):.4f})")


def V6():
    print("\n" + "=" * 78)
    print("V6. SABİT NOKTA — Newton izleri (|h| maks / rms, adım adım)")
    print("=" * 78)
    for v in GAZ:
        for c in ("V0", "VG101", "VG102", "VG103", "VG104", "VL101", "VC101"):
            d = yukle(v, c)
            if d is None or "v162" not in d:
                continue
            iz = d["v162"].get("fix_iz", [])
            if not iz:
                continue
            print(f"   {v:<8} {c:<7} {len(iz)} adım   "
                  f"maks {iz[0][0]:.2e}→{iz[-1][0]:.2e}   "
                  f"rms {iz[0][1]:.2e}→{iz[-1][1]:.2e}   "
                  f"(ortalama aralık 2π/L = {TWO_PI/d['L']:.3f})")


def V7(T, tohum=101):
    print("\n" + "=" * 78)
    print("V7. RED-1 — düz n-uzayı FFT faz karıştırması (ÇALIŞMAYAN reçete)")
    print("    η ve Ĉ, indeks uzayında ortak faz karıştırılır; konumlar")
    print("    yeniden kurulur. Sürüklenme çizgiyi n-tayfında yaydığı için")
    print("    uyumlu toplanma biter.")
    print("=" * 78)
    e, c = C.faz_karistir([T.Y.eta, T.Chat], tohum)
    Y = T._kur(e, c, tohum, "RED1")
    _cizgi_kiyas(T, {"gerçek": T.Y, "RED-1": Y})


def V8(T, tohum=101):
    print("\n" + "=" * 78)
    print("V8. RED-3 — öz-tutarlılık YOK (çizgiler gerçeğin m'sinde)")
    print("=" * 78)
    rng = np.random.default_rng(tohum)
    rot = np.exp(1j * rng.uniform(0, TWO_PI, len(T.w)))
    es, cs = C.faz_karistir([T.eta_sur, T.C_sur], tohum + 500000)
    le, lc = C.cizgi_kur(T.Y.mid, [T.c_eta * rot, T.c_C * rot], T.w)
    Y = T._kur(es + le, cs + lc, tohum, "RED3")
    YG = T.vekil_yerel(tohum, "VG", ayrinti=False)
    _cizgi_kiyas(T, {"gerçek": T.Y, "RED-3": Y, "VG (öz-tutarlı)": YG})


def _cizgi_kiyas(T, gazlar):
    L = T.Y.L
    qm = C.C154.pk_m(int(np.exp(0.86 * L)))
    aq_ = np.array(sorted(qm)); ww = np.log(aq_.astype(float))
    print("    τ       q      " + "".join(f"{k:>16}" for k in gazlar)
          + "     (pow/gp)")
    for tt in (0.46, 0.54, 0.62, 0.70, 0.78):
        j = int(np.argmin(np.abs(ww / L - tt)))
        w = ww[j]
        a = 1.0 / (np.pi * qm[int(aq_[j])] * np.sqrt(float(aq_[j])))
        gp = (2 * a * np.sin(np.pi * w / L)) ** 2
        vals = [Y.cizgi(w)["pow"] / gp for Y in gazlar.values()]
        print(f"    {w/L:.4f} {int(aq_[j]):7d}  " +
              "".join(f"{x:16.4f}" for x in vals))


if __name__ == "__main__":
    V1()
    T = V2("son")
    Yv = V3(T)
    V4(T, Yv)
    V5(T)
    V6()
    V7(T)
    V8(T)
