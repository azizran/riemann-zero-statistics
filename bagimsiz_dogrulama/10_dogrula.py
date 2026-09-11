# -*- coding: utf-8 -*-
"""
SERTİFİKA — bağımsız çoğaltılan tüm manşet iddiaları saklı veriden yeniden üretir.
Kullanım:
   python3 10_dogrula.py           # tam (~1-2 dk)
   python3 10_dogrula.py --hizli   # yalnız Not 1 çekirdeği (~30 s)
Çıkış: PASS/FAIL tablosu + 10_sertifika.json
"""
import json, sys, time, importlib.util
import numpy as np

def yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

o  = yukle("o",  "02_olcum.py")
s3 = yukle("s3", "03_surrogate.py")
s7 = yukle("s7", "07_wv.py")
hizli = "--hizli" in sys.argv
satirlar = []

def kontrol(grup, ad, olculen, beklenen, tol, birim=""):
    ok = abs(olculen - beklenen) <= tol
    satirlar.append(dict(grup=grup, kontrol=ad, olculen=float(olculen), beklenen=float(beklenen),
                         tolerans=float(tol), birim=birim, sonuc="PASS" if ok else "FAIL"))
    print(f"  [{'PASS' if ok else 'FAIL'}] {ad:<48} ölçülen={olculen:9.5f}  hedef={beklenen:9.5f}  ±{tol}")

t0 = time.time()
print("=" * 84)
print("SERTİFİKA — Uğur Sezen / Riemann sıfır programı, bağımsız doğrulama")
print("=" * 84)
zeros = np.sort(np.asarray(np.load("../qm_riemann/128_odl_zeros6_2e6_zeros.npz")["zeros"], float))
print(f"veri: {zeros.size} sıfır (Odlyzko zeros6), t = {zeros[0]:.2f} … {zeros[-1]:.1f}\n")

# ---------- A) Not 1 ----------
print("A) Not 1 — gap–max ortak yasası")
Ls, M2 = [], []
for mc in [1.2e5, 2.0e5, 3.5e5, 6.0e5, 1.0e6]:
    s = o.olc(zeros, mc, 40000); Ls.append(s["L"]); M2.append(s["meanM2"])
a_fit, b_fit = np.polyfit(Ls, M2, 1)
kontrol("A", "Conrey–Ghosh eğimi (A=1.19453)", a_fit, o.A_CG, 0.006)
kontrol("A", "HLPC sabit terimi (2.758)", b_fit, o.ALFA_M1 + 2*o.A_CG, 0.055)
kontrol("A", "r (L=9.86)", o.olc(zeros, 1.2e5, 40000)["r"], 0.7947, 0.005)
kontrol("A", "r (L=11.98)", o.olc(zeros, 1.0e6, 40000)["r"], 0.7589, 0.005)

print("\nB) Gaussian surrogate null")
import subprocess
subprocess.run([sys.executable, "03_surrogate.py", "40"], capture_output=True, text=True)
d = json.load(open("03_sonuc.json"))
kontrol("B", "zeta r (aynı ızgara)", d["r_zeta"], 0.7943, 0.005)
kontrol("B", "Gaussian null seviyesi", d["r_null_ort"], 0.494, 0.010)
kontrol("B", "fazlalık (sigma)", (d["r_zeta"]-d["r_null_ort"])/d["r_null_std"], 43.0, 15.0)

if not hizli:
    # ---------- C) Not 2/3 ----------
    print("\nC) Not 2/3 — toplam kural, kanal saflığı, işaret geçişi")
    d1, M1, tm1 = s7.pencere_verisi(zeros, 1.2e5, 40000)
    L1 = float(np.log(np.mean(tm1)/s7.TWO_PI)); Q1 = s7.taban_olustur(L1)
    _, w1, v1, cM1, cd1 = s7.regresyon(d1, M1, tm1, Q1, gap_kontrol=True)
    _, u1, _, _, _ = s7.regresyon(d1, M1, tm1, Q1, gap_kontrol=False)
    for p in (2, 3, 5, 7):
        kontrol("C", f"toplam kural u(p={p})", u1[p], 1.0, 0.02)
    idx = {q: i for i, (q, p) in enumerate(Q1)}
    k2 = abs(cM1[2+2*idx[2]])/max(abs(cM1[1+2*idx[2]]), 1e-12)
    kontrol("C", "kuadratür/kosinüs (p=2)", k2, 0.0, 0.05)
    kontrol("C", "v(p=2) L=9.86", v1[2], 0.11, 0.03)
    kontrol("C", "v(p=13) L=9.86", v1[13], 0.37, 0.07)
    # k-faktörü: asal kuvvetleri birime dönüyor mu
    k_oku = []
    for (q, p) in Q1:
        k = int(round(np.log(q)/np.log(p)))
        k_oku.append(u1[q]*k)
    kontrol("C", "u·k ortalaması (tüm satırlar)", float(np.mean(k_oku)), 1.0, 0.10)

    # ---------- D) 10^12 penceresi: v(τ) kolapsı ----------
    print("\nD) Örneklem-dışı: v(τ) kolapsı (t oranı ~2,2×10⁶)")
    sat = [l.strip() for l in open("../zeros3.txt") if l.strip()]
    vals = np.array([float(x) for x in sat if x.replace('.', '', 1).replace('-', '', 1).isdigit()])
    zdeep = 267653395647.0 + vals
    dd = np.diff(zdeep); tmd = 0.5*(zdeep[:-1]+zdeep[1:])
    Ld = float(np.log(np.mean(tmd)/s7.TWO_PI)); Qd = s7.taban_olustur(Ld)
    cc = [np.ones_like(dd)]
    for (q, p) in Qd: cc += [np.cos(tmd*np.log(q)), np.sin(tmd*np.log(q))]
    cd, *_ = np.linalg.lstsq(np.column_stack(cc), np.log(dd*Ld/s7.TWO_PI), rcond=None)
    vd = {q: float(np.hypot(cd[1+2*i], cd[2+2*i])*np.sqrt(q)) for i, (q, p) in enumerate(Qd)}
    def asal(n): return n > 1 and all(n % x for x in range(2, int(n**0.5)+1))
    oran = []
    for qs in sorted(v1):
        if not asal(qs): continue
        ts = np.log(qs)/L1
        aday = [q for q in vd if asal(q) and abs(np.log(q)/Ld - ts) < 0.0025]
        if aday: oran.append(vd[min(aday, key=lambda q: abs(np.log(q)/Ld-ts))]/v1[qs])
    kontrol("D", "v(10^12)/v(10^5) eşleşen τ'lerde", float(np.mean(oran)), 1.0, 0.15)

    # ---------- D2) 10^21 ve 10^22: zincirin tamamı ----------
    def asal(n): return n > 1 and all(n % x for x in range(2, int(n**0.5)+1))
    print("\nD2) Derin zincir: v(τ) kolapsı 10^12 → 10^21 → 10^22 (t oranı ~1,1×10^16)")
    def _ofset(dosya):
        sat = [l.strip() for l in open(dosya) if l.strip()]
        return np.array([float(x) for x in sat if x.replace('.','',1).replace('-','',1).isdigit()])
    def _v(ofset, taban):
        dd = np.diff(ofset); t_true = taban + 0.5*(ofset[:-1]+ofset[1:])
        LL = float(np.log(np.mean(t_true)/s7.TWO_PI)); dtl = dd*LL/s7.TWO_PI
        QQ = s7.taban_olustur(LL); zz = 0.5*(ofset[:-1]+ofset[1:]); zz = zz - zz.mean()
        cc = [np.ones_like(dtl)]
        for (q, p) in QQ: cc += [np.cos(zz*np.log(q)), np.sin(zz*np.log(q))]
        cf, *_ = np.linalg.lstsq(np.column_stack(cc), np.log(dtl), rcond=None)
        vv = {q: float(np.hypot(cf[1+2*i], cf[2+2*i])*np.sqrt(q)) for i, (q, p) in enumerate(QQ)}
        return LL, {q: x for q, x in vv.items() if asal(q)}   # yalnız asallar: karışık eğri interpole edilmesin
    B21, B22 = 144176897509546973000.0, 1370919909931995300000.0
    a21, a22 = _ofset("../zeros4.txt"), _ofset("../zeros5.txt")
    L21, V21 = _v(a21, B21); L22, V22 = _v(a22, B22)
    _, V12 = _v(_ofset("../zeros3.txt"), 267653395647.0)
    def _oran(La, Va, Lb, Vb):
        ta = np.array([np.log(q)/La for q in sorted(Va)]); xa = np.array([Va[q] for q in sorted(Va)])
        o_ = np.argsort(ta); ta, xa = ta[o_], xa[o_]
        r = [Vb[q]/np.interp(np.log(q)/Lb, ta, xa) for q in sorted(Vb)
             if ta.min() <= np.log(q)/Lb <= ta.max() and asal(q)]
        return float(np.mean(r)), len(r)
    o21, n21 = _oran(24.475, V12, L21, V21)
    o22, n22 = _oran(24.475, V12, L22, V22)
    kontrol("D2", f"v(10^21)/v(10^12), eşleşen τ (n={n21})", o21, 1.0, 0.20)
    kontrol("D2", f"v(10^22)/v(10^12), eşleşen τ (n={n22})", o22, 1.0, 0.20)

    # ---------- E) Not 4 faz kilidi ----------
    print("\nE) Not 4 — benek faz kilidi (mutlak t_n konvansiyonu)")
    g = zeros[np.argmin(np.abs(zeros-1.2e5))-20000: np.argmin(np.abs(zeros-1.2e5))+20001]
    mid = 0.5*(g[:-1]+g[1:])
    fazlar = []
    for q in [2, 3, 4, 5, 7, 8, 9, 11, 13, 17, 19, 23, 29, 37, 41, 43, 53, 59]:
        G = np.exp(1j*np.log(q)*mid).mean(); faz = abs(abs(np.degrees(np.angle(G)))-180.0)
        fazlar.append(faz)
    kontrol("E", "maks |faz − 180°| (18 satır)", float(np.max(fazlar)), 0.0, 1.0)

    # ---------- F) BBLM sabitleri ----------
    print("\nF) BBLM sabitleri (bağımsız asal toplamı)")
    from mpmath import mp, mpf, log
    mp.dps = 20
    from sympy import primerange
    g0 = mpf('0.577215664901532860606512090082'); g1 = mpf('-0.072815845483676724860586375874')
    P = list(primerange(2, 200000))
    c0 = sum(log(p)**2/(p-1)**2 for p in P); Q_ = sum(log(p)**3/(p-1)**2 for p in P)
    Lam = g0**2 + 2*g1 + c0
    kontrol("F", "Λ = γ₀²+2γ₁+c₀", float(Lam), 1.57314, 0.001)
    kontrol("F", "C = Q/Λ", float(Q_/Lam), 1.4720, 0.001)

n_pass = sum(1 for s in satirlar if s["sonuc"] == "PASS")
print("\n" + "=" * 84)
print(f"SONUÇ: {n_pass}/{len(satirlar)} PASS   ({time.time()-t0:.0f} s)")
if n_pass < len(satirlar):
    print("BAŞARISIZ: " + ", ".join(s["kontrol"] for s in satirlar if s["sonuc"] == "FAIL"))
json.dump(dict(tarih=time.strftime("%Y-%m-%d %H:%M"), satirlar=satirlar, pass_=n_pass, toplam=len(satirlar)),
          open("10_sertifika.json", "w"), indent=1, ensure_ascii=False)
print("yazıldı: 10_sertifika.json")
print("=" * 84)
