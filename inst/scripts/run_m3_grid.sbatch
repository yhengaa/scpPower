#!/bin/bash
#SBATCH --job-name=scpPower-M3
#SBATCH --account=def-xcc_cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=reports/M3/cluster-logs/%A_%a.log
set -euo pipefail
module load r/4.6.1
export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
Rscript inst/scripts/run_m3_grid.R "$SLURM_ARRAY_TASK_ID"
