# Apply only approved E/F; accepted TMT and Derks variance estimates are retained.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
source('R/scp_calibrate.R');started<-proc.time()[['elapsed']];seed<-20261013L
out<-'reports/M1-final';dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(!file.exists('reports/M1-resumed/calibration-results.rds'))source('inst/scripts/calibrate_resumed.R')
old<-readRDS('reports/M1-resumed/calibration-results.rds');manifest<-jsonlite::read_json('reports/M1-resumed/run-manifest.json',simplifyVector=FALSE)
versions<-.scp_versions();vt<-paste(names(versions),versions,sep='=',collapse='; ')
write_csv<-function(d,file){d$software_versions<-rep(vt,nrow(d));write.csv(d,file.path(out,file),row.names=FALSE,na='NA')}
results<-list();summary<-read.csv('reports/M1-resumed/dataset-summary.csv');inputs<-jsonlite::read_json(file.path(out,'tmt-input-manifest.json'),simplifyVector=FALSE)
for(n in names(old)) {
 t0<-proc.time()[['elapsed']];r<-old[[n]];variance_input<-readRDS(paste0('.local/datasets/',n,'-aggregated.rds'))
 missing_input<-variance_input
 if(n%in%c('specht2019v3','leduc2022_pSCoPE')) {
  missing_input<-readRDS(paste0('.local/datasets/',n,'-psm-proteins.rds'))
  z<-missing_input$matrix;fit<-.scp_missingness(z,rep('all',ncol(z)))
  r$missingness<-fit$parameters;r$curve<-fit$curve
  r$manifest$missingness_assays<-paste('Raw reporter-ion PSM assays:',inputs[[n]]$psm_assays)
  r$manifest$raw_psm_input<-inputs[[n]]
  r$manifest$variance_scale<-'Accepted ratio-based peptide-to-protein scale retained; not refitted'
  r$manifest$missingness_scale<-'Mean observed log2 raw reporter intensity; no normalisation'
  if(!fit$parameters$converged)stop(paste('M1 dataset curve unavailable:',n))
 } else if(n=='brunner2022') {
  z<-variance_input$matrix;md<-variance_input$cell_metadata
  fresh<-scp_calibrate(z,md,sample_col='SampleType',condition_col='SampleType',batch_col='Run',biological_unit='culture',normalize_cells=TRUE,absorb_cell_batches=TRUE,seed=seed)
  r$proteins$sigma2_e<-fresh$proteins$sigma2_e;r$proteins$sigma2_b<-fresh$proteins$sigma2_b
  r$proteins$sigma2_e_reason<-fresh$proteins$sigma2_e_reason;r$proteins$sigma2_b_reason<-fresh$proteins$sigma2_b_reason;r$proteins$fit_status<-fresh$proteins$fit_status
  r$line_components<-fresh$line_components;r$manifest$absorb_cell_batches<-TRUE
  r$manifest$variance_scale<-'Total within-line log2 cell variance after per-cell median shift; run and cell noise are one unit'
  for(q in 1:5){j<-r$proteins$tier==q;r$tiers$sigma2_e[q]<-mean(r$proteins$sigma2_e[j],na.rm=TRUE);r$tiers$sigma2_b[q]<-0;r$tiers$sigma2_e_proteins[q]<-sum(is.finite(r$proteins$sigma2_e[j]));r$tiers$sigma2_b_proteins[q]<-sum(j)}
 }
 z<-missing_input$matrix
 r$missingness_proteins<-data.frame(protein=rownames(z),mean_log2=rowMeans(z,na.rm=TRUE),quantified_cells=rowSums(!is.na(z)),missing_fraction=rowMeans(is.na(z)))
 r$missingness_proteins$mean_log2[!is.finite(r$missingness_proteins$mean_log2)]<-NA_real_
 ix<-match(r$proteins$protein,r$missingness_proteins$protein)
 r$proteins$missingness_mean_log2<-r$missingness_proteins$mean_log2[ix]
 r$proteins$missingness_missing_fraction<-r$missingness_proteins$missing_fraction[ix]
 r$manifest$versions<-as.list(versions);r$manifest$seed<-seed;r$manifest$status<-'CALIBRATED';r$manifest$flags<-if(is.finite(r$missingness$slope_high)&&r$missingness$slope_high<0)'FLAG: fitted raw-intensity quantification slope is negative; retain the unconstrained registered fit under E and assess M2 fidelity without retuning' else character();r$manifest$mandate_addendum<-'1.1.2'
 r$manifest$missingness_input_hash<-digest::digest(z,algo='sha256');r$manifest$missingness_cells<-ncol(z)
 r$manifest$variance_input_hash<-digest::digest(variance_input$matrix,algo='sha256');r$manifest$wall_seconds<-proc.time()[['elapsed']]-t0
 r$manifest$retained_abundance_tiers<-'Accepted variance-calibration quintile membership retained; raw missingness means mapped into those tiers by exact protein ID'
 r$manifest$gate_observation_scale<-'Missingness and coverage: raw PSM medians for TMT, quantitative protein inputs for DIA. Cell and batch SD: accepted variance assay after median shifts. Same paired cells and fixed protein universe.'
 for(q in 1:5) {
  j<-r$proteins$tier==q;vals<-r$proteins$missingness_mean_log2[j];r$tiers$mean_log2[q]<-mean(vals,na.rm=TRUE)
  r$tiers$missing_fraction[q]<-mean(r$proteins$missingness_missing_fraction[j],na.rm=TRUE)
 }
 r$tiers$status<-'CALIBRATED';r$tiers$missingness_intercept<-r$missingness$intercept;r$tiers$missingness_slope<-r$missingness$slope
 r$tiers$intensity_center<-r$missingness$intensity_center;r$tiers$effective_slope<-r$missingness$effective_slope;r$tiers$missingness_raw_intercept<-r$missingness$uncentered_intercept;r$tiers$missingness_slope_low<-r$missingness$slope_low;r$tiers$missingness_slope_high<-r$missingness$slope_high
 r$tiers$batch_reason<-paste(unique(r$proteins$sigma2_b_reason),collapse='; ');r$tiers$cell_reason<-paste(unique(r$proteins$sigma2_e_reason),collapse='; ')
 r$tiers$missingness_assay<-if(n%in%c('specht2019v3','leduc2022_pSCoPE'))'Raw PSM reporter intensities; paired accepted cells' else r$manifest$assay
 k<-match(n,summary$dataset);summary$missing_fraction[k]<-mean(is.na(z));summary$mean_proteins_per_cell[k]<-mean(colSums(!is.na(z)));summary$median_proteins_per_cell[k]<-median(colSums(!is.na(z)));summary$minimum_proteins_per_cell[k]<-min(colSums(!is.na(z)))
 summary$quantification_slope[k]<-r$missingness$slope;summary$slope_low[k]<-r$missingness$slope_low;summary$slope_high[k]<-r$missingness$slope_high;summary$status[k]<-'CALIBRATED';summary$sigma2_e_proteins[k]<-sum(is.finite(r$proteins$sigma2_e));summary$sigma2_b_proteins[k]<-sum(is.finite(r$proteins$sigma2_b))
 write_csv(r$proteins,paste0(n,'-proteins.csv'));write_csv(r$missingness_proteins,paste0(n,'-missingness-proteins.csv'));write_csv(r$curve,paste0(n,'-missingness-curve.csv'))
 results[[n]]<-r;cat(n,'CALIBRATED slope',r$missingness$slope,'missing',summary$missing_fraction[k],'\n');flush.console()
}
write_csv(do.call(rbind,lapply(results,`[[`,'tiers')),'calibration-table.csv');write_csv(summary,'dataset-summary.csv')
saveRDS(results,file.path(out,'calibration-results.rds'))
manifest$datasets<-lapply(results,`[[`,'manifest');manifest$versions<-as.list(versions);manifest$mandate<-'1.1 with addenda 1.1.1 and 1.1.2';manifest$seed<-seed
manifest$decision<-list(M1_closed=TRUE,calibrated_datasets=names(results),amendments=c('E','F'),variance_not_retuned=TRUE,M4_criterion_unchanged=TRUE)
manifest$wall_seconds<-proc.time()[['elapsed']]-started;manifest$input_hashes<-list(accepted_results=digest::digest(file='reports/M1-resumed/calibration-results.rds',algo='sha256'),mandate=digest::digest(file='scpPower_mandate_v1.1.md',algo='sha256'))
jsonlite::write_json(manifest,file.path(out,'run-manifest.json'),pretty=TRUE,auto_unbox=TRUE,na='null',digits=NA)
