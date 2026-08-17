"""
50b — SON PARÇA: VARYANS EŞLEMESİYLE N_eff (17 Ağustos 2026)
==============================================================

50'nin resmi: ζ = sıkı çekirdek (r≈0.94) + tam-ağırlıklı asal dalgası.
Bu resim doğruysa N_eff = L + c sabitinin AÇIKLAMASI şu olmalı:
iletilen asal varyansı, CUE'nun N = L + c boyutundaki İÇ gürültüsünü
taklit ediyor. Bağımsız test: r'den değil VARYANSTAN N_eff hesapla.

  Var(log ã_ζ) = 1.3765 (L=12.45)  →  CUE'da hangi N bu varyansı verir?
  Var(log g̃_ζ) = 0.2176           →  aynı soru boşluk kanalı için.

r-tabanlı N_eff (42b): 13.52 (= L + 1.07, o pencerede).
Varyans-tabanlı N_eff de ≈ 13-14 çıkarsa → iki bağımsız gözlemlenebilir
aynı efektif boyutu gösteriyor → kompanzasyon hikayesi sayısal olarak
kapanıyor.
"""

import numpy as np
from pathlib import Path
from scipy.stats import pearsonr

rng = np.random.default_rng(50)
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

def haar_unitary_batch(M, N):
    G = (rng.standard_normal((M, N, N)) + 1j * rng.standard_normal((M, N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(G)
    diag = np.einsum("mii->mi", R)
    return Q * (diag / np.abs(diag))[:, None, :]

def cue_stats(N, n_gaps=200000, grid=24, chunk=1200):
    M = int(np.ceil(n_gaps / N))
    gs, ms = [], []
    for s in range(0, M, chunk):
        m = min(chunk, M - s)
        U = haar_unitary_batch(m, N)
        ph = np.sort(np.angle(np.linalg.eigvals(U)), axis=1)
        gp = np.diff(np.concatenate([ph, ph[:, :1] + TWO_PI], axis=1), axis=1)
        uu = np.arange(1, grid + 1) / (grid + 1)
        th = ph[:, :, None] + gp[:, :, None] * uu[None, None, :]
        df = th[:, :, :, None] - ph[:, None, None, :]
        amp = np.prod(2 * np.abs(np.sin(df / 2)), axis=-1)
        gs.append(gp.ravel()); ms.append(amp.max(axis=-1).ravel())
    g = np.concatenate(gs) * N / TWO_PI
    a = np.concatenate(ms); a /= np.sqrt((a**2).mean())
    return float(np.log(a).var()), float(np.log(g).var())

# ζ hedefleri (L=12.45, 50'den)
V_A_ZETA, V_G_ZETA, L = 1.3765, 0.2176, 12.45

print("CUE log-varyans eğrileri:")
Ns = np.arange(11, 18)
VA, VG = [], []
for N in Ns:
    va, vg = cue_stats(N)
    VA.append(va); VG.append(vg)
    print(f"  N={N:>2}: Var(log ã) = {va:.4f}   Var(log g̃) = {vg:.4f}")

VA, VG = np.array(VA), np.array(VG)
print(f"\nζ hedefleri (L={L}): Var(log ã)={V_A_ZETA}, Var(log g̃)={V_G_ZETA}")
print(f"\n=== N_eff KARŞILAŞTIRMASI (L={L} penceresi) ===")
print(f"  r-tabanlı (42b, Pearson): N_eff = 13.52 (L+1.07)")

# genlik: VA N ile ARTIYOR → küçük varyans küçük N ister
slope = np.polyfit(Ns, VA, 1)[0]
NA_ext = Ns[0] + (V_A_ZETA - VA[0]) / slope
print(f"  GENLİK varyansı: ζ değeri {V_A_ZETA} < tüm CUE aralığı "
      f"(min {VA.min():.3f} @N=11)")
print(f"    → lineer ekstrapolasyon N ≈ {NA_ext:.1f} (L{NA_ext-L:+.1f}) — TERS YÖN!")

# boşluk: eğri düz → eşleme anlamsız
print(f"  BOŞLUK varyansı: CUE eğrisi DÜZ ({VG.min():.4f}–{VG.max():.4f}, "
      f"N'e duyarsız); ζ ({V_G_ZETA}) tümünün altında → eşleme YOK")

print("""
SONUÇ: kompanzasyon hikayesi REDDEDİLDİ. Üç gözlemlenebilir üç farklı yer
gösteriyor (r: L+1; genlik varyansı: ~L−4 yönü; boşluk varyansı: hiçbiri).
ζ hiçbir CUE değil — N_eff = L + c yalnız r-gözlemlenebiliri için geçerli
bir EFEKTİF tanım. Gerçek yapı 50'nin bulduğu: sıkı çekirdek + tam-ağırlıklı
asal dalgası. CUE benzerliği yüzeysel bir denk gelme.""")
