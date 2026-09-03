#!/bin/zsh
# 167 — λ ve kesim gazları hazır olunca ölç (inşa bitişini bekler)
cd /Users/ugursezen/Desktop/arin/deney || exit 1
SCR=/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/167
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
for v in E060 L085 L115 L070; do
  while [ ! -f $SCR/z_$v.npy ]; do sleep 20; done
  sleep 5
  .venv/bin/python qm_riemann/167_configs/167_olcum.py $v 0.40 0.95 0 0 \
      > $SCR/log_C_$v.txt 2>&1
  echo "OLCULDU $v"
done
echo "LAM+KESIM DONE"
