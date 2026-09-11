#!/bin/bash
# -*- coding: utf-8 -*-
# YOL KÖPRÜLERİ — 152-187 arkunun mutlak yollarını yerel karşılıklarına bağlar.
# Hiçbir git-takipli dosyayı DEĞİŞTİRMEZ; bu yüzden sha256 ön-kayıt mühürleri korunur.
# Geri almak: script sonundaki "GERİ ALMA" bölümüne bak.
set -u

DENEY="${DENEY:-/Users/ugur/Desktop/Deney}"
QM="$DENEY/qm_riemann"
UUID="71b922d7-9b7c-485d-9927-3b0db425e2a2"
ESKI_TMP="/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/$UUID"
ESKI_HOME="/Users/ugursezen/Desktop/arin/deney"

echo "== 1) scratchpad köprüsü (bu makinede yazılabilir) =="
mkdir -p "$ESKI_TMP"
if [ -e "$ESKI_TMP/scratchpad" ] && [ ! -L "$ESKI_TMP/scratchpad" ]; then
  echo "  UYARI: $ESKI_TMP/scratchpad zaten var ve symlink değil — dokunulmadı."
else
  ln -sfn "$QM/scratchpad" "$ESKI_TMP/scratchpad"
  echo "  $ESKI_TMP/scratchpad -> $QM/scratchpad"
fi
mkdir -p "$QM/scratchpad"
echo "  yerel scratchpad: $QM/scratchpad (şu an $(ls -1 "$QM/scratchpad" 2>/dev/null | wc -l | tr -d ' ') girdi)"

echo
echo "== 2) QM köprüsü (yönetici izni gerekir) =="
if [ -d "$ESKI_HOME" ]; then
  echo "  $ESKI_HOME zaten var — dokunulmadı."
else
  echo "  Şu komutu bir kez çalıştırman gerekiyor:"
  echo "    sudo mkdir -p /Users/ugursezen/Desktop/arin"
  echo "    sudo ln -sfn \"$DENEY\" /Users/ugursezen/Desktop/arin/deney"
fi

echo
echo "== 3) sınama =="
python3 - <<'PY'
from pathlib import Path
p = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
         "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
print("  scratchpad köprüsü:", "TAMAM" if p.exists() else "YOK")
q = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
print("  QM köprüsü       :", "TAMAM" if q.exists() else "YOK (2. adım)")
PY

# GERİ ALMA:
#   rm -f "$ESKI_TMP/scratchpad"
#   sudo rm -f /Users/ugursezen/Desktop/arin/deney         (ve boşsa üst dizinleri)
