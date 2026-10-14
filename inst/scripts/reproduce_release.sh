#!/bin/bash
# Public input reconstruction followed by the unchanged scientific scripts.
# Use a fresh checkout: existing completed point files are restart receipts.
set -euo pipefail
export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
Rscript -e 'dir.create(".local/R-library",recursive=TRUE,showWarnings=FALSE);if(!requireNamespace("renv",quietly=TRUE))install.packages("renv",repos="https://cloud.r-project.org");renv::restore(library=".local/R-library",prompt=FALSE)'
mkdir -p .local/datasets reports/M1-final reports/M3b/points reports/M5
python inst/scripts/fetch_bulk_metadata.py
python inst/scripts/fetch_resumed_sources.py
Rscript inst/scripts/inspect_resumed.R
Rscript inst/scripts/calibrate_resumed.R
Rscript inst/scripts/extract_tmt_psms.R
Rscript inst/scripts/close_m1.R
# Regenerate G priors, then apply the accepted endpoint-mask correction.
Rscript inst/scripts/generate_presets_revised.R
python inst/scripts/extract_instrument_times.py
Rscript inst/scripts/add_instrument_times.R
Rscript inst/scripts/apply_i_presets.R
# Recorded M2/M4/M4b reports have their own registered reproduction scripts;
# this release regenerates M3b and examples. Unchanged non-TMT results are archived.
Rscript inst/scripts/register_m3b.R
m3b_job=$(sbatch --parsable --wait --account=def-xcc_cpu --cpus-per-task=32 --mem=64G --time=01:00:00 --output=reports/M3b/cluster-reproduction.log --wrap='Rscript inst/scripts/run_m3b_parallel.R')
sacct -j "${m3b_job%%;*}" --format=JobIDRaw,State,ElapsedRaw,TotalCPU,MaxRSS,Submit,Start,End -P > reports/M3b/cluster-accounting.txt
Rscript inst/scripts/report_m3b.R
Rscript inst/scripts/read_m3b_grid.R
python inst/scripts/render_reading_m3b.py
python inst/scripts/finalize_release_reports.py
example_job=$(sbatch --parsable --wait --account=def-xcc_cpu --cpus-per-task=4 --mem=8G --time=00:20:00 --output=reports/M5/examples-reproduction.log --wrap='Rscript inst/scripts/run_release_examples.R')
python inst/scripts/write_release_readme.py
R CMD build . --no-manual
R CMD check scpPower_1.0.0.tar.gz --no-manual --as-cran
