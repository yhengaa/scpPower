#!/bin/bash
# Exact scientific order: public inputs -> J report -> I -> M2b -> M4b.
set -euo pipefail
export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
Rscript -e 'dir.create(".local/R-library",recursive=TRUE,showWarnings=FALSE);if(!requireNamespace("renv",quietly=TRUE))install.packages("renv",repos="https://cloud.r-project.org");renv::restore(library=".local/R-library",prompt=FALSE)'
if ! test -f reports/M1-final/calibration-results.rds; then
  mkdir -p .local/datasets reports/M1-final
  python inst/scripts/fetch_bulk_metadata.py
  python inst/scripts/fetch_resumed_sources.py
  Rscript inst/scripts/inspect_resumed.R
  Rscript inst/scripts/calibrate_resumed.R
  Rscript inst/scripts/extract_tmt_psms.R
  Rscript inst/scripts/close_m1.R
fi
Rscript inst/scripts/diagnostic_j.R
Rscript inst/scripts/audit_j_pipeline.R
Rscript inst/scripts/report_j.R
python inst/scripts/report_tmt_closure.py
Rscript inst/scripts/apply_i_presets.R
m2b_job=$(sbatch --parsable --wait --job-name=scpPower-M2b --cpus-per-task=2 --mem=8G --time=00:15:00 --output=reports/M2b/cluster-reproduction.log --wrap='Rscript inst/scripts/fidelity_m2b.R')
sacct -j "${m2b_job%%;*}" --format=JobID,JobName,State,Elapsed,AllocCPUS,ReqMem,NodeList -P > reports/M2b/cluster-accounting.txt
Rscript inst/scripts/report_m2b.R
Rscript -e 'd<-read.csv("reports/M2b/fidelity-summary.csv");if(!all(d$status=="PASS"))stop("M2b failed: report retained; M4b unauthorised")'
Rscript inst/scripts/prepare_m4b.R
m4b_job=$(sbatch --parsable --wait inst/scripts/run_m4b.sbatch)
sacct -j "${m4b_job%%;*}" --format=JobID,JobName,State,Elapsed,AllocCPUS,ReqMem,NodeList -P > reports/M4b/cluster-accounting.txt
Rscript inst/scripts/report_m4b.R
python inst/scripts/consolidate_m4b_layouts.py
python inst/scripts/finalize_m4b_manifest.py
