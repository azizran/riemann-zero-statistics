# -*- coding: utf-8 -*-
"""
ERRATUM YAMASI — repo'nun literatür notundaki BBLM c0 formülünü düzeltir.
Varsayılan KURU ÇALIŞTIRMA; uygulamak için:  python3 21_erratum_yamasi.py --uygula
(Yedek: aynı dizine <dosya>.bak yazılır.)
"""
import sys, shutil
from pathlib import Path

DOSYA = Path("../qm_riemann/KESIF_SEFERI_KUANTUM_KAOS_29AGU2026.md")
YANLIS = "`c₀ = Σ_p (log p)⁴ Σ_{r≥1} (r−1)r²/p^r`"
DOGRU  = "`c₀ = Σ_p (log p)²/(p−1)²`"
NOT    = ("`c₀ = Σ_p (log p)²/(p−1)²` (düzeltme 11 Eyl 2026: BBLM tanımı "
          "`c_n = [(−1)ⁿ/(2n)!] Σ_p (log p)^{2(n+1)} Σ_r (r−1)r^{2n}/p^r`; "
          "n=0'da `c₀ = Σ_p (log p)²/(p−1)² = 1.3855389…` — eski yazımdaki "
          "`(log p)⁴` ve `r²` transkripsiyon hatasıydı, 33.81 veriyordu; "
          "Λ = 1.573085 ve C = Q/Λ = 1.47161 doğrulandı)")

metin = DOSYA.read_text(encoding="utf-8")
if DOGRU in metin and YANLIS not in metin:
    print("Yama zaten uygulanmış — değişiklik yok."); sys.exit(0)
if YANLIS not in metin:
    print("HATA: beklenen yanlış ifade bulunamadı. Dosya değişmiş olabilir.")
    print(f"  aranan: {YANLIS}"); sys.exit(1)
yeni = metin.replace(YANLIS, NOT, 1)
print(f"dosya : {DOSYA}")
print(f"  −  {YANLIS}")
print(f"  +  {NOT[:120]}…")
if "--uygula" in sys.argv:
    shutil.copy2(DOSYA, str(DOSYA) + ".bak")
    DOSYA.write_text(yeni, encoding="utf-8")
    print(f"UYGULANDI. Yedek: {DOSYA}.bak")
else:
    print("\n(kuru çalıştırma — uygulamak için --uygula ekle)")
