#!/usr/bin/env bash
# Build/run GJK-iterations benchmark and optional plots.
# Usage: ./run_gjk_iterations.sh [--trials N] [--seed S] [--no-build] [--no-plot]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

DATA_DIR="${PACIFIC_BUILDDIR_ABS}/Grains3DGPU/Validations/GJKIterations/data"
CSV_FILE="${DATA_DIR}/gjk_iterations.csv"
PYTHON="${PACIFIC_VENV_PY_ABS}"
EXECUTABLE="GJKIterationsTest"

TRIALS=100
SEED=42
DO_BUILD=1
DO_PLOT=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --trials)  TRIALS="$2";  shift 2 ;;
        --seed)    SEED="$2";    shift 2 ;;
        --no-build) DO_BUILD=0;  shift   ;;
        --no-plot)  DO_PLOT=0;   shift   ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: $0 [--trials N] [--seed S] [--no-build] [--no-plot]"
            exit 1
            ;;
    esac
done


if [[ $DO_BUILD -eq 1 ]]; then
    echo "Building GJKIterationsTest..."
    make all
fi

echo "Running GJK iterations benchmark (trials=$TRIALS seed=$SEED -> $CSV_FILE)"
mkdir -p ${DATA_DIR}
"$EXECUTABLE" "$TRIALS" "$SEED" "$CSV_FILE"

if [[ $DO_PLOT -eq 1 ]]; then
    echo "Generating plots..."
    ${PYTHON} plot.py --csv "$CSV_FILE" --plot-dir ${DATA_DIR}
    echo "Plots written to ${DATA_DIR}"
fi

echo "Done."
