#!/bin/bash

cd "$(dirname "$0")/.." || exit

OUT="scripts/direct_time.csv"
echo "method,N,time_ms" > "$OUT"

METHODS=(0 1 2)
N_LIST=(100 1000 10000 100000)

for m in "${METHODS[@]}"; do
  for n in "${N_LIST[@]}"; do
    echo "Running method $m with N=$n"
    result=$(./bin/tpPoisson1D_direct "$m" "$n")
    time_ms=$(echo "$result" | grep "Elapsed time" | awk '{print $(NF-1)}')
    echo "$m,$n,$time_ms" >> "$OUT"
  done
done

