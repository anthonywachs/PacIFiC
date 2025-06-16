#!/bin/bash
set -e  # exit on any error

# Activate the virtual environment
source /arc/project/st-wachs-1/lmj44/GNN/environment/env_py3.8.10/bin/activate

# Loop over shape, Re, phi, hidden_channels, and N_neigh
for shape in Icosa Cube Octa Dodeca Tetra
do
    for Re in 1 10 100
    do
        for phi in 005 01 02
        do
            for hidden_channels in 84
            do
                for N_neigh in 5 10 20 30 40 50
                do
                    # Create output directory name
                    directory_name="/scratch/st-wachs-1/lmj44/GNN/results/Phi${phi}/Re${Re}/${shape}/"

                    # Create SLURM job script filename
                    job_script="job_${shape}_Re${Re}_Phi${phi}_Neighbour${N_neigh}_Hidden${hidden_channels}.sh"

                    # Write the job script
                    cat << EOF > "$job_script"
#!/bin/bash
#SBATCH --job-name=GNN_Fy_${shape}_Re${Re}_Phi${phi}_Neighbour${N_neigh}_Hidden${hidden_channels}
#SBATCH --account=st-wachs-1-gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --gpus-per-node=1
#SBATCH --output=GNN_Fy_${shape}_Re${Re}_Phi${phi}_Neighbour${N_neigh}_Hidden${hidden_channels}.log.%j
#SBATCH --error=GNN_Fy_${shape}_Re${Re}_Phi${phi}_Neighbour${N_neigh}_Hidden${hidden_channels}.err.%j
#SBATCH --constraint=gpu_mem_32

# Activate virtual environment
source /arc/project/st-wachs-1/lmj44/GNN/environment/env_py3.8.10/bin/activate

# Set CUDA workspace config for determinism
export CUBLAS_WORKSPACE_CONFIG=:4096:8

# Run the training script
python Train_Fx.py --directory_name "$directory_name" --N_neigh $N_neigh --hidden_channels $hidden_channels
EOF

                    # Submit the job
                    sbatch "$job_script"

                    # Optionally remove the script after submission
                    # rm "$job_script"
                done
            done
        done
    done
done
