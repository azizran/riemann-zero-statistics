"""
135 — ÇEKİRDEĞİN KAPALI FORMU: SON TERİM TÜRETİLDİ (27 Ağustos akşamı)
==========================================================================
KALEM TÜRETİMİ (bu akşam): kesin özdeşlik ds=ΔS'ten, kuyruk çizgisinin
p-dalgasına iki modülasyonu:
  (a) genlik: yarım-gap fazından sin(πτ(1+W)) → oran πτ·cot(πτ)
      (üç kez sahneye çıkan geometri fonksiyonunun KAYNAĞI);
  (b) lag-fazı: cos(2πτ(1+W)) → oran −2πτ·sin(2πτ).
Toplamda trigonometrik sadeleşme:
  πτ[cot(πτ)cos(2πτ) − sin(2πτ)] = πτ·cos(3πτ)/sin(πτ).
cos(πτ_p) kinematik çarpanlaşır (mm∓ḡ/2 ortalaması); emilim oranda
sadeleşir. KAPALI FORM:
  κ_çekirdek = Σ_Q A_Q²·πτ_Q·cos(3πτ_Q)/sin(πτ_Q)
             / Σ_Q A_Q²·cos(2πτ_Q)          [kuyruk üzerinden]
Ölçü mucizesi: sürekli limitte a²ρ·dω = dτ/τ (ölçek-değişmez).

SONUÇ:
  AYRIK (D3 kuyruğu, 2100 çizgi): κ_çekirdek = −3.182 — 132c'nin
  SIFIR-PARAMETRE ölçümü −3.22'nin %1.2 İÇİNDE. ✓✓✓
  YASANIN TAM MONTAJI: R_nn = cos(πτ_p)·[κ_ad + κ_çekirdek]
    = cos(πτ_p)·[+1.2 − 3.18] = −1.98·cos(πτ_p) ≈ −2cos(πτ_p) ✓
  DÜRÜST KAYIT: sürekli-limit (dτ/τ) değeri kesim-duyarlı (−1.6;
  tmax'la +0.1..−0.9 salınır) — bu derinlikte ayrık toplam ilk birkaç
  kuyruk çizgisince domine, sürekli yaklaşım henüz geçersiz; gerçek
  gazda evrensel −2'nin kesim-kararlılığı (DW-ağırlıklı kuyrukta
  κ'nın L-değişmezliği) SON açık soru — yarının testi: κ_çekirdek(L)
  ayrık toplamını ζ/ada parametreleriyle L-taraması. κ_ad (+1.2)'nin
  kendi türetimi de açık (adyabatik dal; πτcotπτ-ailesi aday).
"""

import numpy as np
from sympy import primerange, factorint
from scipy.integrate import quad

L0, Lb = 7.0, 8.111

def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

QS = [q for q in pk_list(int(np.exp(1.4 * L0))) if np.log(q) / Lb > 0.52]
num = den = 0.0
for q in QS:
    (pp, kk), = factorint(q).items()
    a = np.log(pp) / (np.pi * np.sqrt(q) * np.log(q))
    tau = np.log(q) / Lb
    A2 = (2 * a * np.sin(np.pi * tau))**2
    num += A2 * np.pi * tau * np.cos(3 * np.pi * tau) / np.sin(np.pi * tau)
    den += A2 * np.cos(2 * np.pi * tau)
print(f"AYRIK (D3, {len(QS)} çizgi): κ_çekirdek = {num/den:+.3f} "
      f"[132c ölçümü: −3.22]")
print(f"MONTAJ: κ_ad(+1.2) + κ_çekirdek({num/den:+.2f}) = "
      f"{1.2 + num/den:+.2f}  [yasa: −2]")
tc, tm = 0.52 * L0 / Lb, 1.4 * L0 / Lb
N, _ = quad(lambda t: np.pi * np.sin(np.pi * t) * np.cos(3 * np.pi * t), tc, tm)
D, _ = quad(lambda t: np.sin(np.pi * t)**2 * np.cos(2 * np.pi * t) / t, tc, tm)
print(f"sürekli limit [{tc:.3f},{tm:.3f}] (dτ/τ ölçüsü): κ = {N/D:+.3f} "
      f"(kesim-duyarlı — dürüst kayıt)")
