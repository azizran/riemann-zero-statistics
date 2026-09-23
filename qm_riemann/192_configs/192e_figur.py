# -*- coding: utf-8 -*-
"""
192e — FİGÜR: 192_uydu_teorisi.png
 Sol: DÜŞÜK pencere Δω havuz profili κ(Δω) (ω-dilim 0.025; ± jk se bandı) —
      yanmalı konumlar yeşil (düz çizgi, ▲), sönmeli konumlar kırmızı (kesikli, ▼);
      işaret içi dolu = ön-kayıtlı hükümle uyumlu, boş = uyumsuz (belirsiz/söner/…).
      Alt şerit: SON penceresi profili (KAYIT; yoğunluk, normalize).
 Sağ: A1 kutupsal çizimler (SON; 5 uydu): her sınıf vektörü K_r / |K_top|, açı
      toplam yönüne göre (toplam = 0° yönünde, uzunluk 1); ± jk açı se yayı.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S192 = SCR / "192"
ONK = json.load(open(S192 / "ONKAYIT_192.json"))
B = json.load(open(S192 / "K1_katalog.json"))
A = json.load(open(S192 / "K2_sinif.json"))
H = json.load(open(S192 / "HUKUM_192.json"))
PF = np.load(S192 / "profiller_192.npz")
KAT = ONK["B"]["katalog"]

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#a3a29c"
YESIL, KIRMIZI = "#008300", "#e34948"
SINIF_RENK = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.titlesize": 10})

fig = plt.figure(figsize=(15.5, 7.6))
gs = fig.add_gridspec(2, 5, width_ratios=[1.35, 1.35, 1, 1, 1], height_ratios=[1, 1],
                      wspace=0.38, hspace=0.62, left=0.05, right=0.985, top=0.88, bottom=0.1)
axL = fig.add_subplot(gs[0, 0:2])
axS = fig.add_subplot(gs[1, 0:2], sharex=axL)

m = PF["merkez"]
ok = (m >= -1.35) & (m <= 2.35) & (PF["d_nkap"] >= 4)
k = PF["d_kap"]
se = np.sqrt(7 / 8 * np.sum((PF["d_kap_reps"] - PF["d_kap_reps"].mean(0)) ** 2, 0))
axL.fill_between(m[ok], (k - se)[ok], (k + se)[ok], color=MUTED, alpha=0.35, lw=0,
                 label="± jk se")
axL.plot(m[ok], k[ok], color=INK, lw=1.3, label="κ(Δω) düşük (L=10.48)")
axL.axhline(0, color=MUTED, lw=0.7)
ymax = np.nanmax(k[ok])
for ad, kk in KAT.items():
    h = kk["delta_omega"]
    r = B["dusuk"][ad]
    yan = kk["liste"] == "yanmali"
    renk = YESIL if yan else KIRMIZI
    beklenen = "yanar" if yan else "söner"
    uyum = r["hukum"] == beklenen
    axL.axvline(h, color=renk, lw=0.9, ls="-" if yan else "--", alpha=0.55, zorder=0)
    y = r["kappa_tepe"]
    axL.plot([r["tepe_merkez"]], [y], marker="^" if yan else "v", ms=7, color=renk,
             mfc=renk if uyum else "white", mew=1.4, zorder=5)
    if yan and ad in ONK["B"]["ana6"]:
        axL.annotate(ad.replace("log", "log "), (r["tepe_merkez"], y), xytext=(0, 8),
                     textcoords="offset points", ha="center", va="bottom",
                     fontsize=7.5, color=INK2)
    if r.get("cukur"):
        axL.annotate(f"{ad}\nçukur {r['kappa_bolu_se']:.1f}σ", (r["tepe_merkez"], y),
                     xytext=(0, -24), textcoords="offset points", ha="center",
                     fontsize=7, color=KIRMIZI)
axL.set_ylim(np.nanmin(k[ok]) * 3.6, ymax * 1.18)
axL.set_ylabel("κ = −Re K_HAVUZ (ω-dilim 0.025)")
axL.set_title("B — seçim kuralı kataloğu, DÜŞÜK pencere Δω profili (birincil)", loc="left")
from matplotlib.lines import Line2D
leg = [Line2D([], [], color=INK, lw=1.3, label="κ(Δω), düşük"),
       plt.Rectangle((0, 0), 1, 1, color=MUTED, alpha=0.35, label="± jk se"),
       Line2D([], [], color=YESIL, marker="^", ls="-", label="yanmalı (μ(a)≠0)"),
       Line2D([], [], color=KIRMIZI, marker="v", ls="--", label="sönmeli (μ(a)=0)"),
       Line2D([], [], color=INK2, marker="o", mfc="white", ls="", label="boş işaret: hüküm ≠ öngörü")]
axL.legend(handles=leg, loc="upper left", fontsize=7.5, frameon=True, framealpha=0.9,
           edgecolor="none", ncol=1)

s = PF["s_kap"]
oks = (m >= -1.35) & (m <= 2.35) & (PF["s_nkap"] >= 4)
ses = np.sqrt(7 / 8 * np.sum((PF["s_kap_reps"] - PF["s_kap_reps"].mean(0)) ** 2, 0))
nd = np.nanmax(k[ok] / 0.025)
ns = np.nanmax(s[oks])
axS.fill_between(m[oks], ((s - ses) / ns)[oks], ((s + ses) / ns)[oks], color=MUTED,
                 alpha=0.35, lw=0)
axS.plot(m[oks], s[oks] / ns, color=INK, lw=1.2, label="son (L=12.03), 188 τ'-dilim → Δω")
axS.plot(m[ok], (k / 0.025)[ok] / nd, color=SINIF_RENK[0], lw=0.9, alpha=0.8,
         label="düşük (üst panel), normalize")
axS.axhline(0, color=MUTED, lw=0.7)
for ad, kk in KAT.items():
    yan = kk["liste"] == "yanmali"
    r = B["son"][ad]
    renk = YESIL if yan else KIRMIZI
    axS.axvline(kk["delta_omega"], color=renk, lw=0.9, ls="-" if yan else "--",
                alpha=0.55, zorder=0)
    uyum = r["hukum"] == ("yanar" if yan else "söner")
    axS.plot([r["tepe_merkez"]], [r["kappa_tepe"] / ns], marker="^" if yan else "v", ms=6,
             color=renk, mfc=renk if uyum else "white", mew=1.3, zorder=5)
axS.set_xlabel("Δω = ω' − L_yerel")
axS.set_ylabel("normalize yoğunluk")
axS.set_title("SON penceresi (ikincil, KAYIT) — aynı kurallar", loc="left")
axS.legend(loc="upper left", fontsize=7.5, frameon=False)

# ---------------- sağ: A1 kutupsal ----------------
U = A["uydular"]
yer = {"+log3": gs[0, 2], "-log3": gs[0, 3], "+log2": gs[0, 4],
       "+log5": gs[1, 2], "+log6": gs[1, 3]}
ONG_VEK = {"+log3": [(-60, 1.0), (60, 1.0)], "-log3": [(0, 0.5), (0, 0.5)],
           "+log5": [(-108, 1.0), (-36, 1.0), (36, 1.0), (108, 1.0)],
           "+log6": [(60, 1.0), (-60, 1.0)], "+log2": [(0, 0.5), (0, 0.5)]}
ong = {"+log3": "öngörü: r1,r2 ±60° (120° ayrık), oran 1",
       "-log3": "öngörü: aynı yön, oran ≈ 0.5",
       "+log5": "öngörü: 72° adımlar, oran 1",
       "+log6": "öngörü: r1,r5 ∓60°, oran 1",
       "+log2": "kontrol: mod-4 sınıfları aynı yön"}
for ad, g in yer.items():
    ax = fig.add_subplot(g, projection="polar")
    u = U[ad]
    ax.set_theta_zero_location("E")
    ax.set_rlim(0, 1.6)
    ax.set_rticks([0.5, 1.0, 1.5])
    ax.tick_params(labelsize=7)
    ax.annotate("", xy=(0, 1.0), xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=2.2))
    for th_o, r_o in ONG_VEK[ad]:
        ax.plot([0, np.radians(th_o)], [0, r_o], color=MUTED, lw=1.2, ls=(0, (3, 2)))
    ci = 0
    satir = []
    for et, c in u["siniflar"].items():
        if et == "bolunen":
            continue
        th = np.radians(c["aci_goreli_top"])
        r_ = c["oran_top"]
        col = SINIF_RENK[ci % 4]
        ax.annotate("", xy=(th, r_), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.8))
        sd = np.radians(c["se_aci_goreli"])
        tt = np.linspace(th - sd, th + sd, 30)
        ax.plot(tt, np.full_like(tt, r_), color=col, lw=3, alpha=0.35)
        satir.append((col, f"{et.replace('mod4', '₄')}: {r_:.2f}∠{c['aci_goreli_top']:+.0f}°"))
        ci += 1
    for kk_, (col, tx) in enumerate(satir):
        x0 = 0.02 + (kk_ % 2) * 0.52
        y0 = -0.16 - 0.09 * (kk_ // 2)
        ax.text(x0, y0, "■", transform=ax.transAxes, ha="left", va="top", fontsize=8,
                color=col)
        ax.text(x0 + 0.08, y0, tx, transform=ax.transAxes, ha="left", va="top",
                fontsize=7.2, color=INK2)
    ff = u["farklar"]
    fs = "; ".join(f"Δ{a_.replace('mod4', '').replace('->', '→')}={v['fark']:+.0f}°±{v['se']:.0f}"
                   for a_, v in ff.items() if a_ != "1->4")
    ax.set_title(f"{ad}  |K_top|={u['mod_top']:.4f}\n{fs}\n{ong[ad]}", fontsize=7.5,
                 color=INK, pad=14)
axT = fig.add_subplot(gs[1, 4])
axT.axis("off")
hk = H["hukum"]
axT.text(0, 0.95, "HÜKÜM (ön-kayıtlı)", fontsize=9.5, fontweight="bold", color=INK,
         va="top")
import textwrap
axT.text(0, 0.84, textwrap.fill("H-192a: " + hk["H-192a"], 34), fontsize=7.5, color=INK,
         va="top")
axT.text(0, 0.50, "H-192b: " + hk["H-192b"] + " (+log3 sınıfları\naynı yönde; dönme yok)",
         fontsize=7.5, color=INK, va="top")
axT.text(0, 0.34, "H-192c: " + hk["H-192c"], fontsize=7.5, color=INK, va="top")
axT.text(0, 0.20, "Kutupsal: siyah ok = toplam K_top (0° yönü,\nuzunluk 1); renkli ok = sınıf "
         "K_r/|K_top|,\naçı toplama göre; yay = ± jk açı se;\ngri kesikli = KALEM öngörüsü.",
         fontsize=7, color=INK2, va="top")
fig.suptitle("192 — Uyduların teorisi: seçim kuralı μ(a)/φ(a) (B) ve kalıntı-sınıfı "
             "dönmesi (A1)", fontsize=12, color=INK, x=0.05, ha="left")
out = QM / "192_uydu_teorisi.png"
fig.savefig(out, dpi=150)
print(f"-> {out}")
