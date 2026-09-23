#!/bin/bash
# yeniden_kur.sh — 184-187 ara önbelleklerini GERÇEK SIFIRLARDAN yeniden üretir
# (188-K0 reçetesi). Betikler AYNEN koşulur; sonunda 187'nin mühürlü sayılarıyla kapı.
# Önce: ./araclar/senkron.sh (köprüler kurulmuş olmalı). Süre ~40 dk.
set -e
REPO="$(cd "$(dirname "$0")/.." && pwd)"
QM="$REPO/qm_riemann"
S="/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad"
P="$REPO/.venv/bin/python"; [ -x "$P" ] || P=python3
LOG="$S/yeniden_kur"; mkdir -p "$S"/{155,164,184,185,186,187} "$LOG"
cd "$QM"
adim() { echo "=== $1 ($(date +%H:%M)) ==="; }
adim "1/7 eta_son (155 üreticisi)";   $P -u "$REPO/araclar/k0_eta_son.py"      > "$LOG/eta_son.txt" 2>&1; tail -2 "$LOG/eta_son.txt"
adim "2/7 184b gercek";               $P -u "$REPO/araclar/k0_184b_gercek.py"  > "$LOG/184b.txt" 2>&1;    tail -1 "$LOG/184b.txt"
adim "3/7 Hkeskin inşası (~11 dk)";   $P -u 164_configs/164_insa.py Hkeskin    > "$LOG/164.txt" 2>&1;     tail -1 "$LOG/164.txt"
adim "4/7 184b2 Hkeskin";             $P -u 184_configs/184b2_zarf_ikiz.py Hkeskin="$S/155/z_Hkeskin.npy" > "$LOG/184b2.txt" 2>&1; tail -1 "$LOG/184b2.txt"
adim "5/7 185b";                      $P -u 185_configs/185b_oz_muhasebe.py    > "$LOG/185b.txt" 2>&1;    tail -1 "$LOG/185b.txt"
adim "6/7 186b";                      $P -u 186_configs/186b_g1_oztutarlilik.py > "$LOG/186b.txt" 2>&1;   tail -1 "$LOG/186b.txt"
adim "7/7 187c + KAPI";               $P -u 187_configs/187c_zeta_defteri.py   > "$LOG/187c.txt" 2>&1
grep -i "HAVUZ" "$LOG/187c.txt" | head -3
echo "KAPI: HAVUZ gerçek ζ = 0.3287±0.0023 ∠179.96°, ikiz 0.0768±0.0064 olmalı (187 mührü)."
