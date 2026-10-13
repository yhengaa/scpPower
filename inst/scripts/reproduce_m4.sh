#!/bin/bash
# Run from the repository root with the registered renv dependencies restored.
# Public-source downloads and hashes precede M4; released presets are not written.
set -euo pipefail
export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
Rscript -e 'dir.create(".local/R-library", recursive=TRUE, showWarnings=FALSE); if (!requireNamespace("renv",quietly=TRUE)) install.packages("renv",repos="https://cloud.r-project.org"); renv::restore(library=".local/R-library",prompt=FALSE)'
if ! test -f reports/M1-final/calibration-results.rds; then
  mkdir -p .local/datasets reports/M1-final
  python inst/scripts/fetch_bulk_metadata.py
  python inst/scripts/fetch_resumed_sources.py
  Rscript inst/scripts/inspect_resumed.R
  Rscript inst/scripts/calibrate_resumed.R
  Rscript inst/scripts/extract_tmt_psms.R
  Rscript inst/scripts/close_m1.R
fi
Rscript inst/scripts/prepare_m4.R
# Requires Slurm on the target cluster. Wait for this job to finish before report.
m4_batch_id=$(sbatch --parsable --wait inst/scripts/run_m4.sbatch)
sacct -j "${m4_batch_id%%;*}" --format=JobID,JobName,State,Elapsed,AllocCPUS,ReqMem,NodeList -P > reports/M4/cluster-accounting.txt
Rscript inst/scripts/report_m4.R
Rscript inst/scripts/read_m3_grid.R
python inst/scripts/render_reading_grid.py
python inst/scripts/consolidate_m4_layouts.py
python inst/scripts/finalize_m4_manifest.py
