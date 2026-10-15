#!/bin/bash
# Run ContactTablePerformance across capacities and load factors.

set -e

PYTHON=$PACIFIC_VENV_PY_ABS
EXECUTABLE="ContactForces"
DATADIR="${PACIFIC_BUILDDIR_ABS}/Grains3DGPU/Validations/ContactForces/data"
CSV_FILE="${DATADIR}/contact_table_bench.csv"

TRIALS=3
SEED=42

mkdir -p ${DATADIR}
rm -f "$CSV_FILE"

echo "================================================================================"
echo "Running Contact Table benchmarks"
echo "================================================================================"
echo "Output file: $CSV_FILE"
echo "Trials per configuration: $TRIALS"
echo "Random seed: $SEED"
echo

capacities=(2048 8192 32768 131072 524288 2097152)
load_factors=(0.25 0.5 0.75)

total=$(( ${#capacities[@]} * ${#load_factors[@]} ))
current=0

for cap in "${capacities[@]}"; do
  for lf in "${load_factors[@]}"; do
    current=$((current + 1))
    echo "--------------------------------------------------------------------------------"
    echo "Configuration $current/$total"
    echo "  capacity: $cap"
    echo "  loadFactor: $lf"
    echo "  trials: $TRIALS"
    echo "--------------------------------------------------------------------------------"

    # Executable args: <capacity> <loadFactor> <numLookups> <numTrials> <csv>
    numLookups=$((cap / 2))
    $EXECUTABLE "$cap" "$lf" "$numLookups" "$TRIALS" "$CSV_FILE"

    echo
  done
done

echo "================================================================================"
echo "All benchmarks completed!"
echo "Results saved to: $CSV_FILE"
echo "================================================================================"
