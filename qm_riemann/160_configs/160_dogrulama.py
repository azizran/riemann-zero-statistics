"""
160 — DENETİM
V1  160'ın φ, τ_eff, σ_φ, |Γ|, δ, M1, M0, A·S'si = 159'unki mi? (BİT DÜZEYİ)
V2  Çizgi düzeyinde kapanış: E1 + i·Es = [zp conj(zc) e^{−iA}/4]/⟨ρ⟩
V3  Bant düzeyinde muhasebe kapanışı: δ − (K1+K2)          [X̃ değişkeni]
V4  Değişken ikamesi artığı: δ − M0                         [dsΔ değişkeni]
V5  Kümülant yeniden kurulumu: K1 vs A·κ1 − A³κ3/6 (ve sönüm)
V6  Muhasebenin ANALİTİK kimliği — sentetik kontrol: bilinen tek çizgi,
    bilinen adım dizisi; K1, K2 ve δ elle hesaplananla karşılaştırılır.
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
C = importlib.import_module("160_cekirdek")
SCR160 = C.SCR
SCR159 = Path(str(C.SCR).replace("/160", "/159"))
TWO_PI = 2 * np.pi


def yuk(scr, veri, taban, g):
    p = scr / f"S_{veri}_t{taban}_{g}.json"
    return json.load(open(p)) if p.exists() else None


CIFT = [(v, t, g) for g in ("g158", "t1")
        for v, t in (("son", 0.4), ("son", 0.52), ("orta", 0.4),
                     ("keskin", 0.4), ("A4", 0.4), ("A4", 0.52), ("P1", 0.4))]


def V1():
    print("=" * 78)
    print("V1 — 160'ın ölçüm zinciri = 159'unki mi? (bit düzeyi)")
    print("=" * 78)
    alan = ["phi", "tau_eff", "sPhi_jk", "absG", "delta", "sDelta_jk",
            "M1", "M0", "absGc", "rho", "kinematik"]
    top = {a: 0.0 for a in alan}
    tS = 0.0
    nb = nc = 0
    for v, t, g in CIFT:
        a = yuk(SCR160, v, t, g)
        b = yuk(SCR159, v, t, g)
        if a is None or b is None:
            continue
        A = {x["tau"]: x for x in a["bantlar"] if x.get("olculdu")}
        B = {x["tau"]: x for x in b["bantlar"] if x.get("olculdu")}
        ort = sorted(set(A) & set(B))
        if not ort:
            continue
        nc += 1
        nb += len(ort)
        d = {al: max(abs(A[k][al] - B[k][al]) for k in ort) for al in alan}
        dS = max(max(abs(A[k]["AS"][s] - B[k]["AS"][s])
                     for s in ("S_tam", "S_lad", "S_eta", "S_dri", "S_xtil"))
                 for k in ort)
        for al in alan:
            top[al] = max(top[al], d[al])
        tS = max(tS, dS)
        print(f"  {v:7s} t{t:<5} {g:5s} {len(ort):2d} bant  "
              f"|Δφ|={d['phi']:.1e} |Δτe|={d['tau_eff']:.1e} "
              f"|Δδ|={d['delta']:.1e} |ΔM1|={d['M1']:.1e} "
              f"|ΔM0|={d['M0']:.1e} |ΔA·S|={dS:.1e}")
    print(f"\n  → {nc} (gaz,taban,ızgara) üçlüsü, {nb} bant. MAKS fark:")
    for al in alan:
        print(f"      {al:10s} {top[al]:.3e}")
    print(f"      {'AS[*]':10s} {tS:.3e}")


def V2_V5():
    print()
    print("=" * 78)
    print("V2–V5 — muhasebenin kapanışı ve kümülant yeniden kurulumu")
    print("=" * 78)
    print("  koşu                 bant  V2 çizgi-kapanış  V3 |δ−(K1+K2)|  "
          "V4 |δ−M0|   V5 |K1−(Aκ1−A³κ3/6)|")
    g2 = g3 = g4 = 0.0
    for p in sorted(SCR160.glob("S_*.json")):
        d = json.load(open(p))
        B = [b for b in d["bantlar"] if b.get("olculdu")]
        if not B:
            continue
        v2 = max(b["kapanis"] for b in B)
        v3 = max(abs(b["delta"] - b["Kfull"]) for b in B)
        v4 = max(abs(b["delta"] - b["M0"]) for b in B)
        v5 = max(abs(b["K1"] - (b["KUM"]["A1k1"] - b["KUM"]["A3k3"] / 6.0))
                 for b in B)
        g2 = max(g2, v2); g3 = max(g3, v3); g4 = max(g4, v4)
        print(f"  {p.stem:20s} {len(B):3d}   {v2:14.2e}  {v3:14.2e}  "
              f"{v4:10.2e}  {v5:16.4f}")
    print(f"\n  → V2 maks {g2:.2e} (bağıl) · V3 maks {g3:.2e} · "
          f"V4 maks {g4:.2e}")
    print("  (V5 bir KAPANIŞ değil, üçüncü-mertebe KESMENİN artığıdır — "
          "büyüklüğü §2'de yorumlanır.)")


def V6():
    """Sentetik kontrol: TEK bilinen çizgi + bilinen adım dizisi.
    ρ, σ, K1, K2 ve δ elle hesaplanır; ölçüm cebriyle karşılaştırılır."""
    print()
    print("=" * 78)
    print("V6 — SENTETİK KONTROL (bilinen çizgi, bilinen adım dizisi)")
    print("=" * 78)
    rng = np.random.default_rng(7)
    N = 200000
    L = 12.0
    tau = 0.62
    Wf = tau * L
    A = TWO_PI * Wf / L
    Xt = rng.normal(0.0, 0.22, N)          # bilinen adım sapması
    Xt -= Xt.mean()
    dm = (TWO_PI / L) * (1.0 + Xt)
    mid = np.concatenate(([1000.0], 1000.0 + np.cumsum(dm)))
    eta = np.cos(Wf * mid) + 0.35 * rng.normal(0, 1, N + 1)
    e0, e1, m0 = eta[:-1], eta[1:], mid[:-1]
    cw, sw = np.cos(Wf * m0), np.sin(Wf * m0)
    cw1, sw1 = np.cos(Wf * mid[1:]), np.sin(Wf * mid[1:])
    zc = complex(2 * np.mean(e0 * cw), -2 * np.mean(e0 * sw))
    zp = complex(2 * np.mean(e1 * cw), -2 * np.mean(e1 * sw))
    zb = zc / 2.0
    c1r, c1i = e1 * cw1, -e1 * sw1
    rho = c1r * zb.real + c1i * zb.imag
    sig = c1i * zb.real - c1r * zb.imag
    rm = rho.mean()
    # ELLE: e^{iA·X̃} doğrudan X̃'den (ölçümdeki kısayol DEĞİL)
    Xtil = (mid[1:] - mid[:-1]) * L / TWO_PI - 1.0
    ph = np.exp(1j * A * Xtil)
    E1 = complex(np.mean(rho * ph.real), np.mean(rho * ph.imag)) / rm
    Es = complex(np.mean(sig * ph.real), np.mean(sig * ph.imag)) / rm
    K1 = np.angle(E1)
    Kf = np.angle(E1 + 1j * Es)
    dlt = np.angle(zp * np.conj(zc) * np.exp(-1j * A))
    print(f"  τ={tau}  A={A:.4f}  N={N}  σ_X={Xt.std():.3f}  A·σ_X={A*Xt.std():.3f}")
    print(f"  ELLE   : K1 = {K1:+.8f}   K1+K2 = {Kf:+.8f}")
    print(f"  ÖLÇÜM  : δ  = {dlt:+.8f}   (arg[zp conj(zc) e^{{−iA}}])")
    print(f"  |δ − (K1+K2)| = {abs(dlt-Kf):.3e}    K2 = {Kf-K1:+.8f}")
    # kısayolun (yeni trigonometri yok) elle hesapla aynılığı
    dr = cw1 * cw + sw1 * sw
    di = sw1 * cw - cw1 * sw
    ca, sa = np.cos(A), np.sin(A)
    Cx, Sx = dr * ca + di * sa, di * ca - dr * sa
    print(f"  KISAYOL vs ELLE  maks|Δcos| = {np.max(np.abs(Cx-ph.real)):.2e}  "
          f"maks|Δsin| = {np.max(np.abs(Sx-ph.imag)):.2e}")
    # birinci mertebe ve kümülant
    S = float((np.dot(rho, Xt) / N - rm * Xt.mean()) / rm)
    u1 = float(np.dot(rho, Xt) / N / rm); u2 = float(np.dot(rho, Xt**2) / N / rm)
    u3 = float(np.dot(rho, Xt**3) / N / rm)
    k3 = u3 - 3 * u1 * u2 + 2 * u1**3
    print(f"  A·S = {A*S:+.6f}   A·κ1 − A³κ3/6 = {A*u1 - A**3*k3/6:+.6f}   "
          f"K1 = {K1:+.6f}   δ = {dlt:+.6f}")
    print(f"  → kesme çarpanı K1/(A·S) = {K1/(A*S):+.3f} ; "
          f"kuadratür çarpanı δ/K1 = {dlt/K1:+.3f}")
    # σ'nın kuadratür kimliği
    dro = np.concatenate((np.diff(rho), np.zeros(1)))
    a11 = np.dot(rho, rho); a12 = np.dot(rho, dro); a22 = np.dot(dro, dro)
    b1 = np.dot(rho, sig); b2 = np.dot(dro, sig)
    det = a11 * a22 - a12 * a12
    al = (a22 * b1 - a12 * b2) / det; be = (a11 * b2 - a12 * b1) / det
    R2 = (al * b1 + be * b2) / np.dot(sig, sig)
    print(f"  σ ≈ α·ρ + β·Δρ :  α = {al:+.4f} (öngörü tan(A/2) = "
          f"{np.tan(A/2):+.4f})  β = {be:+.4f} (öngörü 1/sinA = "
          f"{1/sa:+.4f})  R² = {R2:.4f}")


if __name__ == "__main__":
    V1()
    V2_V5()
    V6()
