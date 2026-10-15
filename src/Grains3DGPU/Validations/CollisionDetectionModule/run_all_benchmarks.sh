PYTHON=$PACIFIC_VENV_PY_ABS
DATADIR="${PACIFIC_BUILDDIR_ABS}/Grains3DGPU/Validations/CollectionDetectionModule/data"
EXECUTABLE="CDModuleBenchmark"
CSV="${DATADIR}/benchmark_comprehensive.csv"

# Simple Cartesian-product benchmark: vary N x shape x aspect ratio
set -e

SEED=42
NPLATP="--platform both --precision both"

R=0.05
SP_S="$R $R $R"
BOX1=$(python3 -c "import math; a=0.05/math.sqrt(3); print(f'{a:.6f} {a:.6f} {a:.6f}')")
BOX4_X=$(python3 -c "import math; a=0.05/math.sqrt(18); b=4*a; print(f'{a:.6f}')")
BOX4_Y=$BOX4_X
BOX4_Z=$(python3 -c "import math; a=0.05/math.sqrt(18); b=4*a; print(f'{b:.6f}')")
SQ4_X=$BOX4_X
SQ4_Y=$BOX4_Y
SQ4_Z=$BOX4_Z

mkdir -p ${DATADIR}
rm -f "$CSV"

N_LIST=(512 2048 8192 32768 131072)
SHAPES=(sphere box_ar1 box_ar4 superquadric_ar4)

# Domain multipliers matching the attached script
SPHERE_DOM=(32 48 80 112 176)
BOX1_DOM=(32 48 80 112 176)
BOX4_DOM=(31 48 80 128 200)
SQ4_DOM=(31 48 80 128 200)

total_configs=$(( ${#N_LIST[@]} * ${#SHAPES[@]} ))
current=0

run_benchmark() {
    local particles=$1
    local domain=$2
    local shape=$3
    local size_x=$4
    local size_y=$5
    local size_z=$6
    local aspect=$7

    current=$((current + 1))
    echo "================================================================================"
    echo "Configuration $current/$total_configs"
    echo "  Particles: $particles"
    echo "  Domain: [$domain, $domain, $domain]"
    echo "  Shape: $shape"
    echo "  Size: [$size_x, $size_y, $size_z]"
    echo "  Aspect: $aspect"
    echo "================================================================================"

    local append_flag=""
    [ $current -gt 1 ] && append_flag="--append"

    "$EXECUTABLE" \
        --particles "$particles" \
        --domain "$domain" \
        --shape "$shape" \
        --size "$size_x" "$size_y" "$size_z" \
        --aspect "$aspect" \
        $NPLATP \
        --warmup 2 --measure 5 \
        --csv "$CSV" --seed "$SEED" \
        $append_flag

    echo ""
}

for shape_key in "${SHAPES[@]}"; do
    for idx in "${!N_LIST[@]}"; do
        N=${N_LIST[$idx]}
        case "$shape_key" in
            sphere)
                dom_mul=${SPHERE_DOM[$idx]}
                domain=$(python3 -c "print(${dom_mul} * ${R})")
                run_benchmark $N "$domain" sphere $R $R $R 1.0
                ;;
            box_ar1)
                dom_mul=${BOX1_DOM[$idx]}
                domain=$(python3 -c "print(${dom_mul} * ${R})")
                run_benchmark $N "$domain" box $BOX1 $BOX1 $BOX1 1.0
                ;;
            box_ar4)
                dom_mul=${BOX4_DOM[$idx]}
                domain=$(python3 -c "print(${dom_mul} * ${R})")
                run_benchmark $N "$domain" box $BOX4_X $BOX4_Y $BOX4_Z 4.0
                ;;
            superquadric_ar4)
                dom_mul=${SQ4_DOM[$idx]}
                domain=$(python3 -c "print(${dom_mul} * ${R})")
                run_benchmark $N "$domain" superquadric $SQ4_X $SQ4_Y $SQ4_Z 4.0
                ;;
        esac
    done
done

echo "================================================================================"
echo "All benchmarks completed."
echo "Results: $CSV"
echo "================================================================================"
