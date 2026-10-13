#!/bin/bash
#SBATCH --job-name=scpPower-M4b
#SBATCH --account=def-xcc_cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=reports/M4b/cluster-logs/%j.log
set -euo pipefail
module load r/4.6.1
export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
Rscript inst/scripts/run_m4b_cluster.R
