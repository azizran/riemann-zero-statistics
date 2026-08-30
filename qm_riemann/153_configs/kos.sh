#!/bin/zsh
S=/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/153
PY=/Users/ugursezen/Desktop/arin/deney/.venv/bin/python
G=/Users/ugursezen/Desktop/arin/deney/qm_riemann/153_configs/153_gaz.py
for c in R1 R2 R3 N3 N5 N10 N5x Pd P1 P0; do
  echo "@@@@ $c basladi $(date +%H:%M:%S)"
  $PY $G $c > $S/log_$c.txt 2>&1
  echo "@@@@ $c bitti rc=$? $(date +%H:%M:%S)"
done
echo "@@@@ TARAMA TAMAM $(date +%H:%M:%S)"
