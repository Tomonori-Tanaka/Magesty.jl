#!/usr/bin/env bash
# Run thread-scaling benchmark with 1, 2, and 4 threads.
# Usage: bash test/run_thread_bench.sh [julia-executable]
set -euo pipefail

JULIA="${1:-julia}"
SCRIPT="$(cd "$(dirname "$0")" && pwd)/benchmark_threads.jl"
LOGDIR="$(cd "$(dirname "$0")" && pwd)/thread_bench_results"
mkdir -p "$LOGDIR"

echo "Julia: $("$JULIA" --version)"
echo

for NT in 1 2 4; do
    LOG="$LOGDIR/threads_${NT}.txt"
    echo "========== -t $NT =========="
    "$JULIA" -t "$NT" "$SCRIPT" 2>&1 | tee "$LOG"
    echo
done

echo "========== SCALING SUMMARY =========="
echo "build_design_matrix_torque median times:"
for NT in 1 2 4; do
    LOG="$LOGDIR/threads_${NT}.txt"
    VAL=$(grep "build_design_matrix_torque" "$LOG" | awk '{print $NF, $(NF-1)}' | head -1)
    echo "  threads=$NT : $VAL"
done
