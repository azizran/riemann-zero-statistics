"""
167 — ORTAK TANIMLAR: ablasyon ailesi tablosu ve gaz künyeleri
==============================================================
Her gazın merdiveni:   S(t) = −Σ_q (λ·a_q·w_q) sin(ω_q t),  a_q = 1/(πm√q)
    λ      = GENLİK ÖLÇEĞİ ekseni
    τ_ust  = KESKİN kesim (w_q = 1{τ_q ≤ τ_ust})
    (τ_c,Δ)= ERFC kesimi   (w_q = ½erfc((τ_q−τ_c)/Δ))
    pencere= kullanılan sıfır penceresi (T ekseni): (bas, bit) dilimi

`Hkeskin` (λ=1, keskin τ≤1.00) ve `HA4` (λ=1, erfc 0.68/0.125) 164'ün
gazlarıdır; `son` gerçek zeta sıfırlarının son 300k'sıdır (merdiven
NOMİNALDİR, inşa edilmemiştir — λ/kesim tanımsız).
"""

# ad -> künye
KUNYE = {
    # --- taban (164/165/166'nın gazları) ---
    "son":      dict(gercek=True, lam=None, tau_ust=None, pen=None,
                     aile="taban", T=1.0),
    "Hkeskin":  dict(gercek=False, lam=1.00, tau_ust=1.00, pen=None,
                     aile="taban", T=1.0),
    "HA4":      dict(gercek=False, lam=1.00, tau_ust=1.00, pen=(0.68, 0.125),
                     aile="taban", T=1.0),
    # --- (a) GENLİK ÖLÇEĞİ ekseni (keskin τ≤1.00) ---
    "L070":     dict(gercek=False, lam=0.70, tau_ust=1.00, pen=None,
                     aile="lam", T=1.0),
    "L085":     dict(gercek=False, lam=0.85, tau_ust=1.00, pen=None,
                     aile="lam", T=1.0),
    "L115":     dict(gercek=False, lam=1.15, tau_ust=1.00, pen=None,
                     aile="lam", T=1.0),
    # --- (c) KESİM ekseni (λ=1) ---
    "K090":     dict(gercek=False, lam=1.00, tau_ust=0.90, pen=None,
                     aile="kesim", T=1.0),
    "K070":     dict(gercek=False, lam=1.00, tau_ust=0.70, pen=None,
                     aile="kesim", T=1.0),
    "E060":     dict(gercek=False, lam=1.00, tau_ust=1.00, pen=(0.60, 0.125),
                     aile="kesim", T=1.0),
}

# --- (b) PENCERE ekseni: aynı gazın alt-pencereleri -------------------
#     ad -> (kaynak gaz, bas_kesir, bit_kesir)
DILIM = {
    "HkT2a": ("Hkeskin", 0.0, 0.5),   "HkT2b": ("Hkeskin", 0.5, 1.0),
    "HkT4a": ("Hkeskin", 0.0, 0.25),  "HkT4b": ("Hkeskin", 0.25, 0.5),
    "HkT4c": ("Hkeskin", 0.5, 0.75),  "HkT4d": ("Hkeskin", 0.75, 1.0),
    "snT2a": ("son", 0.0, 0.5),       "snT2b": ("son", 0.5, 1.0),
    "snT4a": ("son", 0.0, 0.25),      "snT4d": ("son", 0.75, 1.0),
}
for _ad, (_k, _a, _b) in DILIM.items():
    KUNYE[_ad] = dict(KUNYE[_k], aile="pencere", T=round(_b - _a, 4),
                      dilim=(_k, _a, _b))


def pencere_dict():
    """165_cekirdek.PENCERE'ye eklenecek erfc pencereleri."""
    return {ad: k["pen"] for ad, k in KUNYE.items() if k.get("pen")}


def lam_of(ad):
    return KUNYE[ad].get("lam") or 1.0


def tau_ust_of(ad):
    return KUNYE[ad].get("tau_ust") or 1.00
