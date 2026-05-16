#!/usr/bin/env bash
# Run thread-scaling benchmark for build_design_matrix_energy and
# build_design_matrix_torque at multiple thread counts and print a
# scaling summary.
#
# Usage:
#   bash bench/run_threads_scaling.sh
#
# Environment overrides:
#   THREAD_COUNTS   space-separated thread counts (default: "1 2 4 8")
#   JULIA           julia executable (default: julia)
#   BENCH_INPUT     forwarded to the bench script
#   BENCH_SAMPLES   forwarded to the bench script
#   BENCH_SECONDS   forwarded to the bench script
#
# Per-thread-count logs land in bench/results/threads_N.txt.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPT="$REPO_ROOT/bench/benchmark_threads_scaling.jl"
LOGDIR="$REPO_ROOT/bench/results"
mkdir -p "$LOGDIR"

JULIA="${JULIA:-julia}"
THREAD_COUNTS="${THREAD_COUNTS:-1 2 4 8}"

echo "Julia: $("$JULIA" --version)"
echo

for NT in $THREAD_COUNTS; do
    LOG="$LOGDIR/threads_${NT}.txt"
    echo "========== threads=$NT =========="
    "$JULIA" --project="$REPO_ROOT/bench" -t "$NT" "$SCRIPT" 2>&1 | tee "$LOG"
    echo
done

echo "========== SCALING SUMMARY =========="
for KERNEL in build_design_matrix_energy build_design_matrix_torque; do
    echo "$KERNEL:"
    BASE_MS=""
    for NT in $THREAD_COUNTS; do
        LOG="$LOGDIR/threads_${NT}.txt"
        LINE="$(grep "^$KERNEL:" "$LOG" || true)"
        if [ -z "$LINE" ]; then
            printf "  threads=%-3s : (no data)\n" "$NT"
            continue
        fi
        MS="$(printf '%s' "$LINE" | awk -F'median=' '{print $2}' | awk '{print $1}')"
        if [ -z "$BASE_MS" ]; then
            BASE_MS="$MS"
            SPEEDUP="1.00"
        else
            SPEEDUP="$(awk -v b="$BASE_MS" -v m="$MS" 'BEGIN { if (m > 0) printf "%.2f", b/m; else print "n/a" }')"
        fi
        printf "  threads=%-3s : median=%s ms  speedup=%sx\n" "$NT" "$MS" "$SPEEDUP"
    done
done
