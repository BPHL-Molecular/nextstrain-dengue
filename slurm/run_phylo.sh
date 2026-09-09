#!/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=dengue-phylo
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128gb
#SBATCH --time=48:00:00
#SBATCH --output=nextstrain_%j.out
#SBATCH --error=nextstrain_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arnold.rodriguezhilario@flhealth.gov

module load conda
conda activate nextstrain

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

nextstrain build phylogenetic \
    --configfile build-configs/florida/config.yaml \
    --rerun-incomplete \
    --printshellcmds \
    --show-failed-logs
