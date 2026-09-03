#!/bin/zsh
# 167 — PENCERE EKSENİ ölçümleri (alt-pencereler; N küçük)
cd /Users/ugursezen/Desktop/arin/deney || exit 1
SCR=/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/167
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
for grup in "HkT2a HkT2b" "HkT4a HkT4b" "HkT4c HkT4d" "snT2a snT2b" "snT4a snT4d"; do
  for v in ${=grup}; do
    .venv/bin/python qm_riemann/167_configs/167_olcum.py $v 0.40 0.95 0 0 \
        > $SCR/log_C_$v.txt 2>&1 &
  done
  wait
done
echo "PENCERE DONE"
