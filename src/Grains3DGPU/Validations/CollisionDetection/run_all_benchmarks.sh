#!/bin/bash
# Run all collision detection benchmarks (shape families below).

set -e

PYTHON=$PACIFIC_VENV_PY_ABS
EXECUTABLE="CDBenchmark"
DATADIR="${PACIFIC_BUILDDIR_ABS}/Grains3DGPU/Validations/CollectionDetection/data"
CSV_FILE="${DATADIR}/collision_benchmark_comprehensive.csv"
TRIALS=3
SEED=42

r=0.05

BOX_SIZE=$(echo "scale=10; 1 / sqrt(3) * $r" | bc -l)

BOX4_X=$(echo "scale=10; $r / 2" | bc -l)
BOX4_Y=$(echo "scale=10; $r / 2" | bc -l)
BOX4_Z=$(echo "scale=10; 4 * $r" | bc -l)

SQ4_X=$(echo "scale=10; $r / 2" | bc -l)
SQ4_Y=$(echo "scale=10; $r / 2" | bc -l)
SQ4_Z=$(echo "scale=10; 4 * $r" | bc -l)

mkdir -p ${DATADIR}
rm -f "$CSV_FILE"

echo "================================================================================"
echo "Running Comprehensive Collision Detection Benchmark"
echo "================================================================================"
echo "Output file: $CSV_FILE"
echo "Trials per configuration: $TRIALS"
echo "Random seed: $SEED"
echo ""

total_configs=20
current=0

run_benchmark() {
    local particles=$1
    local domain=$2
    local shape=$3
    local size_x=$4
    local size_y=$5
    local size_z=$6
    local aspect=$7
    local append_flag=$8
    
    current=$((current + 1))
    echo "================================================================================"
    echo "Configuration $current/$total_configs"
    echo "  Particles: $particles"
    echo "  Domain: [$domain, $domain, $domain]"
    echo "  Shape: $shape"
    echo "  Size: [$size_x, $size_y, $size_z]"
    echo "  Aspect: $aspect"
    echo "================================================================================"
    
    # compute-sanitizer --tool memcheck 
    # cuda-gdb --args \
    $EXECUTABLE \
        --particles "$particles" \
        --domain "$domain" \
        --shape "$shape" \
        --size "$size_x" "$size_y" "$size_z" \
        --aspect "$aspect" \
        --precision both \
        --platform both \
        --trials "$TRIALS" \
        --seed "$SEED" \
        --csv "$CSV_FILE" \
        $append_flag
    
    echo ""
}

# Spheres - S1 (aspect ratio 1.0)
echo "Running Sphere benchmarks (S1)..."

run_benchmark 512 $(echo "32 * $r" | bc -l) sphere $r $r $r 1.0 ""
run_benchmark 2048 $(echo "48 * $r" | bc -l) sphere $r $r $r 1.0 "--append"
run_benchmark 8192 $(echo "80 * $r" | bc -l) sphere $r $r $r 1.0 "--append"
run_benchmark 32768 $(echo "112 * $r" | bc -l) sphere $r $r $r 1.0 "--append"
run_benchmark 131072 $(echo "176 * $r" | bc -l) sphere $r $r $r 1.0 "--append"

# Boxes - B1 (aspect ratio 1.0)
echo "Running Box benchmarks (B1)..."

run_benchmark 512 $(echo "32 * $r" | bc -l) box $BOX_SIZE $BOX_SIZE $BOX_SIZE 1.0 "--append"
run_benchmark 2048 $(echo "48 * $r" | bc -l) box $BOX_SIZE $BOX_SIZE $BOX_SIZE 1.0 "--append"
run_benchmark 8192 $(echo "80 * $r" | bc -l) box $BOX_SIZE $BOX_SIZE $BOX_SIZE 1.0 "--append"
run_benchmark 32768 $(echo "112 * $r" | bc -l) box $BOX_SIZE $BOX_SIZE $BOX_SIZE 1.0 "--append"
run_benchmark 131072 $(echo "176 * $r" | bc -l) box $BOX_SIZE $BOX_SIZE $BOX_SIZE 1.0 "--append"

# Superquadrics - S4 (aspect ratio 4.0)
echo "Running Superquadric benchmarks (S4)..."

run_benchmark 512 $(echo "31 * $r" | bc -l) superquadric $SQ4_X $SQ4_Y $SQ4_Z 4.0 "--append"
run_benchmark 2048 $(echo "48 * $r" | bc -l) superquadric $SQ4_X $SQ4_Y $SQ4_Z 4.0 "--append"
run_benchmark 8192 $(echo "80 * $r" | bc -l) superquadric $SQ4_X $SQ4_Y $SQ4_Z 4.0 "--append"
run_benchmark 32768 $(echo "128 * $r" | bc -l) superquadric $SQ4_X $SQ4_Y $SQ4_Z 4.0 "--append"
run_benchmark 131072 $(echo "200 * $r" | bc -l) superquadric $SQ4_X $SQ4_Y $SQ4_Z 4.0 "--append"

# Boxes - B4 (aspect ratio 4.0)
echo "Running Box benchmarks (B4)..."

run_benchmark 512 $(echo "31 * $r" | bc -l) box $BOX4_X $BOX4_Y $BOX4_Z 4.0 "--append"
run_benchmark 2048 $(echo "48 * $r" | bc -l) box $BOX4_X $BOX4_Y $BOX4_Z 4.0 "--append"
run_benchmark 8192 $(echo "80 * $r" | bc -l) box $BOX4_X $BOX4_Y $BOX4_Z 4.0 "--append"
run_benchmark 32768 $(echo "128 * $r" | bc -l) box $BOX4_X $BOX4_Y $BOX4_Z 4.0 "--append"
run_benchmark 131072 $(echo "200 * $r" | bc -l) box $BOX4_X $BOX4_Y $BOX4_Z 4.0 "--append"

echo "================================================================================"
echo "All benchmarks completed!"
echo "Results saved to: $CSV_FILE"
echo "================================================================================"