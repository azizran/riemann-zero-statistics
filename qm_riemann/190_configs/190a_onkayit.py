# -*- coding: utf-8 -*-
"""
190a — K0b: KURAL-ÖNCE ÖN-KAYIT (TARAĞIN EVRENSELLİĞİ — düşük pencere)
=====================================================================
KALEM_TARAK_EVRENSELLIK_23EYL2026 AYNEN. K0 zinciri ('dusuk') ve makine mührü
('son' → 187 ζ_g hane hane) GEÇTİKTEN SONRA, HARİTA ÖLÇÜLMEDEN donan:
hedef/rakip Δω konumları, ω-dilim ızgarası ve Δω aralığı, 8 blok ve L_b
(kinematik — mid'den; sonuç değil), tepe konumu tanımı, H-190a/b/c eşikleri
ve hüküm kuralları, KAYIT tanımları, makine mühürleri, girdi dosyalarının sha'sı.

Çıktı: scratchpad/190/ONKAYIT_190.json (sha256 = bu betiğin sha'sı + damga).
"""
import hashlib
import json
import subprocess
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S190 = SCR / "190"
ZD = S190 / "zincir_dusuk"
TWO_PI = 2 * np.pi

# ---------------- kinematik (sonuç DEĞİL: mid'den L, L_b) ----------------
E = np.load(ZD / "eta_dusuk_t0.4_c4000.npz")
MID = np.asarray(E["mid"], float)
L = float(E["L"])
N = len(MID)
KJ = np.linspace(0, N, 9).astype(int)            # 184-188 jk blokları AYNEN
LW = np.log(MID / TWO_PI)
L_B = [float(LW[KJ[b]:KJ[b + 1]].mean()) for b in range(8)]     # 188f(c) AYNEN
L_B_YAYILIM = [float(LW[KJ[b + 1] - 1] - LW[KJ[b]]) for b in range(8)]
L_SON = 12.029593241726252                       # 188 ONKAYIT L

# ---------------- hedef ve rakip konumlar ----------------
LOG2, LOG3, LOG6 = float(np.log(2)), float(np.log(3)), float(np.log(6))
HEDEF_HUKUM = {"0": 0.0, "+log2": LOG2, "+log3": LOG3, "+log6": LOG6}
HEDEF_KAYIT = {"-log2": -LOG2, "-log3": -LOG3}
# τ'-sabit rakip: tepe son penceredeki τ' = 1 ± log n / L_son'da kalır →
# düşük pencerede Δω = ± log n · L_dusuk / L_son (KALEM: 0.604 / 0.957 / 1.562)
RAKIP = {"+log2": LOG2 * L / L_SON, "+log3": LOG3 * L / L_SON,
         "+log6": LOG6 * L / L_SON, "-log2": -LOG2 * L / L_SON,
         "-log3": -LOG3 * L / L_SON, "0": None}
RAKIP_KALEM = {"+log2": 0.604, "+log3": 0.957, "+log6": 1.562}

DW = 0.025          # ω-dilim genişliği
PENCERE = 0.15      # tepe arama yarı-genişliği
TOL = 0.05          # H-190a konum toleransı


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


ONKAYIT = {
    "gorev": "190 — TARAĞIN EVRENSELLİĞİ: örneklem-dışı sınav (düşük pencere)",
    "kalem": "KALEM_TARAK_EVRENSELLIK_23EYL2026.md AYNEN",
    "pencere": {"ad": "dusuk", "tanim": "zeros6 Z[200000:500000] (155_kos.PENCERE)",
                "L": L, "N": N, "L_min": float(LW.min()), "L_max": float(LW.max()),
                "t_lo": float(MID[0]), "t_hi": float(MID[-1])},
    "L_son": L_SON,
    "girdi_sha256": {
        "eta_dusuk_t0.4_c4000.npz": sha(ZD / "eta_dusuk_t0.4_c4000.npz"),
        "K1_gercek_dusuk.npz": sha(ZD / "K1_gercek_dusuk.npz"),
        "OZ_gercek_dusuk.npz": sha(ZD / "OZ_gercek_dusuk.npz"),
        "G1_proj_gercek_dusuk.npz": sha(ZD / "G1_proj_gercek_dusuk.npz"),
        "188b_harita.py": sha(QM / "188_configs" / "188b_harita.py"),
        "188d_analiz.py": sha(QM / "188_configs" / "188d_analiz.py"),
        "188f_kesif.py": sha(QM / "188_configs" / "188f_kesif.py")},
    # ---------------- bloklar ----------------
    "bloklar": {
        "tanim": ("8 eşit bitişik blok, kenar = linspace(0, N, 9).astype(int) "
                  "(184-188 jk blokları AYNEN); L_b = blok içi mean log(m_n/2π) "
                  "(188f(c) AYNEN)."),
        "kenar": KJ.tolist(), "L_b": L_B, "blok_ici_L_yayilimi": L_B_YAYILIM},
    # ---------------- τ'-dilim haritası (H-190b) ----------------
    "dilim_izgara": [round(0.86 + 0.005 * k, 3) for k in range(89)],
    "dilim_tanim": ("188 AYNEN: dilim s = τ' ∈ (e_s, e_{s+1}], 88 × 0.005, "
                    "(0.86, 1.30]; τ' = log q'/L_dusuk; 187b.asal_kuvvetler "
                    "(q ≤ e^{1.30·L_dusuk}); a = Λ/(π√q log q)."),
    "ince_bant_kenar": [round(0.45 + 0.01 * k, 2) for k in range(41)] + [0.86],
    "bant8_kenar": [0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.86],
    "havuz": "HAVUZ = pencere çizgileri τ ∈ [0.45, 0.86) (τ = w/L_dusuk; K1 evreni)",
    "makine": ("188b AYNEN (importlib ile 188b_harita modülü: seri_ve_G, izdusum, "
               "c_proj, K_matris, G_loo, bant_maskeleri, karisim, dilimleri_kur); "
               "karışım = c^kesik − c^öz (187c AYNEN, 'dusuk' zinciri, tam örneklem, "
               "loo); Ĝ_mid MUTLAK m_n; S = Σ a'|Ĝ|; W = Σ a'; "
               "S_Re = Σ a'·πτ'·cos(πτ')·Re Ĝ_mid(ω')."),
    "orneklem": ("TAM örneklem (ALT = 1; N = 299 999) — görev: alt-örneklem "
                 "gerekmez; 188'den tek sapma budur (188: her 3. nokta). jk: 8-blok loo."),
    "analiz": ("188d.analiz AYNEN: 41×88 karmaşık K, ağırlıksız SVD; v = conj(V1), "
               "faz Σv reel-pozitif; corr = Pearson (88 dilim); her loo replikası "
               "baştan (SVD dahil)."),
    # ---------------- ω-dilim (Δω) blok profili (H-190a) ----------------
    "omega_dilim": {
        "genislik": DW,
        "tanim": ("Her blok b için, pencere-ötesi çizgi q' (τ' ∈ (0.86, 1.30], "
                  "L_dusuk ile) Δω = log q' − L_b ekseninde j = floor(Δω/0.025 + 0.5) "
                  "dilimine atanır (dilim merkezi j·0.025; dilim [ (j−½)·0.025, "
                  "(j+½)·0.025 ) ). Dilim serisi dds_{b,j}(n) YALNIZ blok b "
                  "noktalarında 188b.seri_ve_G ile (AYNEN); izdüşüm blok b "
                  "noktalarında: C_{b,j}(q) = 2·Σ_{n∈b} dds e^{−iω_q m_n}/n_b "
                  "(188f(c) tek-blok deseni); K_HAVUZ(b,j) = Σ_{q∈HAVUZ} C conj(mix)"
                  "/Σ|mix|², mix TAM örneklem (188f(c) AYNEN); κ_b(j) = −Re K_HAVUZ(b,j)."),
        "delta_omega_araligi": {
            f"blok{b}": [0.86 * L - L_B[b], 1.30 * L - L_B[b]] for b in range(8)},
        "dilim_kapsami": "en az bir çizgi içeren dilimler (uçlarda kısmi dilimler olabilir)"},
    "havuz_profili": ("8 bloğun HAVUZ κ profili: κ_b(j) (blok başına, Δω dilimleri "
                      "bloklar arasında hizalı); şekil/figür için bloklar toplamı "
                      "κ_Σ(j) = Σ_b κ_b(j)·n_b/N."),
    # ---------------- tepe konumu ----------------
    "tepe_tanim": ("Blok başına: merkezi |j·0.025 − hedef| ≤ 0.15 (+1e-9) olan "
                   "dilimler arasında κ_b(j) maksimumu → tepe_b = j·0.025 (dilim "
                   "merkezi; alt-dilim ara değerleme YOK). Blok-medyan = "
                   "np.median(tepe_b, 8 blok). se: blok-loo jackknife (7 bloğun "
                   "medyanı; se = √(7/8·Σ(m_i − m̄)²)) — medyanın jk se'si "
                   "kaba bir ölçüdür, hükme girmez."),
    "hedef_hukum": HEDEF_HUKUM,
    "hedef_kayit": HEDEF_KAYIT,
    "rakip_tau_sabit": RAKIP,
    "rakip_kalem_yuvarlak": RAKIP_KALEM,
    "rakip_not": ("Bragg (Δω=0) için rakip sayısal konum yok (τ'-sabit okuma "
                  "blok-medyanda ≈ L − medyan L_b ≈ 0): yalnız ±0.05 koşulu. "
                  "+log6 için rakip 1.562 arama penceresinin [1.642, 1.942] "
                  "DIŞINDA: rakip-yakınlık koşulu tepe > 1.677 ile eşdeğer. "
                  "+log2 (ölüm kolu) ayırıcı sınavdır: arama penceresi [0.543, 0.843] "
                  "rakibi (0.604) içerir; orta nokta 0.6485."),
    # ---------------- hipotezler ----------------
    "H_190a": {
        "hipotez": ("Δω EVRENSELLİĞİ: düşük pencerede tepeler Δω = 0, +log2, +log3, "
                    "+log6'da (yükseklikten bağımsız)."),
        "kosul": ("her hedef için: |blok-medyan tepe − hedef| ≤ 0.05 VE (rakip "
                  "varsa) |medyan − hedef| < |medyan − rakip|."),
        "olum": ("+log2 tepesi: |medyan − 0.6931| > 0.05 YA DA |medyan − 0.6931| ≥ "
                 "|medyan − 0.6041|."),
        "hukum": ("dört hedefin hepsi koşulu sağlarsa MÜHÜR; ölüm koşulu sağlanırsa "
                  "ÖLDÜ; aksi halde (log2 geçer, diğerlerinden biri geçmez) KAYIT "
                  "(kısmi). H-190c tutarsa → ÖLDÜ (tarak resmi ölür).")},
    "H_190b": {
        "hipotez": "İMZALI ÖNGÖRÜCÜ: S_Re, Re v'yi izler (88 τ'-dilim, 188 ızgarası).",
        "kosul": "corr(Re v, S_Re) ≤ −0.80 → MÜHÜR",
        "olum": "|corr(Re v, S_Re)| < 0.60 → ÖLDÜ",
        "ara": ("−0.80 < corr ve |corr| ≥ 0.60 → KAYIT (pozitif işaretse 'işaret "
                "ters' notu). H-190c tutarsa → ÖLDÜ.")},
    "H_190c": {
        "hipotez": "YAPISIZ SIFIR: |corr(Re v, W)| ≥ |corr(Re v, S_Re)|.",
        "hukum": ("koşul sağlanırsa H-190c MÜHÜR ve tarak resmi ÖLÜR (H-190a ve "
                  "H-190b ÖLDÜ yazılır); sağlanmazsa H-190c ÖLDÜ.")},
    # ---------------- KAYIT (eşiksiz) ----------------
    "KAYIT": {
        "zeta_g": ("düşük pencere ζ_g (τ_c = 0.86; 187c AYNEN; 8 bant + HAVUZ; "
                   "modül + açı ± jk) — son (187 K2) ile yan yana."),
        "zeta_130": "HAVUZ Σ_{τ'≤1.20} K ve Σ_{τ'≤1.30} K (188b AYNEN) + ζ_g'ye oranı.",
        "rang1": "rang-1 payı, corr(|v|,S), corr(|v|,W) (188d AYNEN).",
        "negatif_hedefler": "−log2, −log3 blok-medyan tepeleri (rakip −0.604/−0.957).",
        "egim": ("blok tepelerinin L_b'ye eğimi (en küçük kareler, 8 blok): "
                 "Δω-evrensel → 0; τ'-sabit → −1 (eşiksiz, yorum)."),
        "son_profili": ("SON penceresi blok-blok Δω profili 188'in dosyalarından "
                        "(harita_proj_gercek.npz, 184/185/186 son zinciri; yeniden "
                        "koşu YOK): κ_b(s) = −Re K_HAVUZ(b,s), Δω_b(s) = τ'_s·L_son − "
                        "L_b,son (τ'-dilim çözünürlüğü 0.005·L_son = 0.060); aynı "
                        "tepe kuralı (merkez ±0.15)."),
        "dusuk_tau_dilim_profili": ("düşük pencerenin τ'-dilim (0.005·L = 0.052) "
                                    "blok profili — ω-dilim profilinin çapraz "
                                    "kontrolü (aynı tepe kuralı).")},
    # ---------------- makine mühürleri ----------------
    "makine_muhurleri": {
        "M1": "seri_ve_G ≡ 187b.katman_serisi (düşük, dilim 0, tam örneklem): maks|Δ| = 0",
        "M2": ("188'in dilim_0.npz'i (son, alt-örneklem) 190b kod yoluyla bit-bit "
               "yeniden üretilir (dds, Gre, Gim)"),
        "M3": ("doğrusallık: her blokta Σ_j K_HAVUZ(b,j) (ω-dilimleri) = Σ_s "
               "K_HAVUZ(b,s) (τ'-dilimleri, blok izdüşümü) — göreli fark ≤ 1e-10"),
        "K0": "190k0 'son' zinciri 155/184/185/186/187 dosyalarını bit-bit verdi"},
    "jackknife": "8-blok loo (184-188 AYNEN), se = √(7/8·Σ(θ_i−θ̄)²)",
    "kurallar": ("TEK DALGA; ölümler kurtarmasız; eşik gevşetme yok; ölçülemeyen "
                 "'erişilemedi'; git'e dokunulmaz; sonuç ORTAK TEFTİŞE."),
    "zaman": None,
    "sha256": None,
}

if __name__ == "__main__":
    try:
        ONKAYIT["zaman"] = subprocess.check_output(["date"], text=True).strip()
    except Exception:
        ONKAYIT["zaman"] = time.strftime("%a %b %d %H:%M:%S %Z %Y")
    s = hashlib.sha256(Path(__file__).resolve().read_bytes()).hexdigest()
    ONKAYIT["sha256"] = s
    yol = S190 / "ONKAYIT_190.json"
    if yol.exists():
        raise SystemExit(f"ONKAYIT_190 ZATEN VAR — üzerine yazılmaz: {yol}")
    json.dump(ONKAYIT, open(yol, "w"), indent=1, ensure_ascii=False)
    print("=" * 78)
    print("190a / K0b ÖN-KAYIT yazıldı")
    print("=" * 78)
    print(f"  sha256 = {s}")
    print(f"  damga  = {ONKAYIT['zaman']}")
    print(f"  L_dusuk = {L:.9f}  N = {N}  L∈[{LW.min():.4f},{LW.max():.4f}]")
    print(f"  L_b = {np.round(L_B, 4).tolist()}")
    print(f"  blok-içi L yayılımı = {np.round(L_B_YAYILIM, 4).tolist()}")
    print(f"  hedef = {HEDEF_HUKUM}")
    print(f"  rakip = { {k: (round(v, 4) if v is not None else None) for k, v in RAKIP.items()} }")
    print(f"  -> {yol}")
