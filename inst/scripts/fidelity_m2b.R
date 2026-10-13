# Registered M2 gate. Definitions and tolerances are fixed before simulation.
# Twenty replicates at beta=0 in each dataset's exact sample-by-batch counts.
# Missingness/coverage use the paired raw mask; SD uses the accepted variance
# scale, with the same per-cell median shift applied to observed and simulated.
# Tier tests are individual, without averaging away a failing tier.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
source('R/scp_calibrate.R');source('R/scp_simulate.R')
seed<-20261013L;started<-proc.time()[['elapsed']];out<-'reports/M2b';dir.create(out,recursive=TRUE,showWarnings=FALSE)
rs<-readRDS('reports/M1-final/calibration-results.rds')[c('specht2019v3','leduc2022_pSCoPE')];vs<-.scp_versions();vt<-paste(names(vs),vs,sep='=',collapse='; ')
write_csv<-function(d,path){d$software_versions<-rep(vt,nrow(d));write.csv(d,file.path(out,path),row.names=FALSE,na='NA')}
normalize_cells<-function(z){med<-apply(z,2,median,na.rm=TRUE);sweep(z,2,med-median(med,na.rm=TRUE),'-')}
row_sd<-function(z)apply(z,1,function(v){v<-v[is.finite(v)];if(length(v)>1)sd(v) else NA_real_})
safe_median<-function(v)if(any(is.finite(v)))median(v[is.finite(v)]) else NA_real_
measure<-function(missing_matrix,variance_matrix,md,tier,lf){
 z<-normalize_cells(variance_matrix);lines<-unique(md$sample)
 cell_sds<-vapply(lines,function(line)row_sd(z[,md$sample==line,drop=FALSE]),numeric(nrow(z)))
 if(!is.matrix(cell_sds))cell_sds<-matrix(cell_sds,nrow=nrow(z))
 batch_sds<-rep(NA_real_,nrow(z))
 if(!lf) {
  line_mean<-vapply(lines,function(line)rowMeans(z[,md$sample==line,drop=FALSE],na.rm=TRUE),numeric(nrow(z)))
  if(!is.matrix(line_mean))line_mean<-matrix(line_mean,nrow=nrow(z))
  residual<-z-line_mean[,match(md$sample,lines),drop=FALSE]
  batches<-split(seq_len(nrow(md)),md$batch)
  mixed<-batches[vapply(batches,function(j)length(unique(md$sample[j]))>1L,logical(1))]
  if(length(mixed)>1L) {
   means<-vapply(mixed,function(j)rowMeans(residual[,j,drop=FALSE],na.rm=TRUE),numeric(nrow(z)))
   means[!is.finite(means)]<-NA_real_;batch_sds<-row_sd(means)
  }
 }
 do.call(rbind,lapply(sort(unique(tier[!is.na(tier)])),function(t){j<-tier==t&!is.na(tier);data.frame(tier=t,missing_fraction=mean(is.na(missing_matrix[j,,drop=FALSE])),mean_proteins_per_cell=mean(colSums(!is.na(missing_matrix))),half_cells_fraction=mean(rowMeans(!is.na(missing_matrix))>=.5),cell_sd=safe_median(cell_sds[j,,drop=FALSE]),batch_sd=if(lf)NA_real_ else safe_median(batch_sds[j]),cell_sd_pairs=sum(is.finite(cell_sds[j,,drop=FALSE])),batch_sd_proteins=sum(is.finite(batch_sds[j])))}))
}
intensity_bins<-function(z,breaks){mu<-rowMeans(z,na.rm=TRUE);miss<-rowMeans(is.na(z));bins<-cut(mu,breaks=breaks,include.lowest=TRUE,labels=FALSE);do.call(rbind,lapply(sort(unique(bins[!is.na(bins)])),function(b){j<-which(bins==b);data.frame(bin=b,mean_log2=mean(mu[j]),missing_fraction=mean(miss[j]),proteins=length(j))}))}
all_checks<-list();all_reps<-list();all_observed<-list();all_summary<-list();all_bins<-list();all_half<-list();manifests<-list()
for(di in seq_along(rs)) {
 dataset_started<-proc.time()[['elapsed']];n<-names(rs)[di];r<-rs[[di]];cat('M2 native design',n,'\n');flush.console()
 obs<-readRDS(paste0('.local/datasets/',n,'-aggregated.rds'));raw<-obs # Amendment I: endpoint quantitative mask, same accepted cells/proteins
 # Both observed quantities use identical cell/protein identifiers. No new
 # filtering is introduced and all fixed accepted tier labels are retained.
 ids<-r$proteins$protein;ov<-obs$matrix[ids,,drop=FALSE];om<-raw$matrix[ids,colnames(ov),drop=FALSE]
 md<-data.frame(sample=as.character(obs$cell_metadata[[r$manifest$sample_col]]),batch=as.character(obs$cell_metadata[[r$manifest$batch_col]]),row.names=colnames(ov))
 samples<-unique(md$sample);batches<-unique(md$batch)
 layout<-unclass(table(factor(md$sample,levels=samples),factor(md$batch,levels=batches)));storage.mode(layout)<-'double'
 ncounts<-c(ceiling(length(samples)/2),floor(length(samples)/2));cellcounts<-as.integer(rowSums(layout))
 lf<-n=='brunner2022';observed<-measure(om,ov,md,r$proteins$tier,lf)
 observed$dataset<-n;all_observed[[n]]<-observed
 means<-rowMeans(om,na.rm=TRUE);breaks<-unique(c(-Inf,as.numeric(quantile(means[is.finite(means)],seq(0,1,length.out=21),names=FALSE)),Inf))
 obin<-intensity_bins(om,breaks);obin$replicate<-0L;obin$dataset<-n;obin$source<-'observed';bins<-list(obin)
 replicates<-list();run_records<-list();layout_df<-data.frame(sample=rep(rownames(layout),ncol(layout)),batch=rep(colnames(layout),each=nrow(layout)),cells=as.vector(layout));write_csv(layout_df,paste0(n,'-native-layout.csv'))
 # Hash the preset before the first draw; the gate never writes calibration,
 # preset, or simulator code in response to outcomes.
 preset_path<-paste0('inst/presets/',n,'.yaml');preset_hash<-digest::digest(file=preset_path,algo='sha256')
 for(rep in seq_len(20L)) {
  simseed<-as.integer(seed+1000L*di+rep)
  sim<-scp_simulate(r,ncounts,cellcounts,beta=0,preset=preset_path,layout=layout,k=1,seed=simseed)
  sm<-sim$cell_metadata;metrics<-measure(sim$matrix,sim$matrix,sm,sim$protein_metadata$tier,lf)
  metrics$dataset<-n;metrics$replicate<-rep;metrics$seed<-simseed;replicates[[rep]]<-metrics;run_records[[rep]]<-sim$manifest;run_records[[rep]]$input_hashes<-as.list(sim$manifest$input_hashes)
  bin<-intensity_bins(sim$matrix,breaks);bin$replicate<-rep;bin$dataset<-n;bin$source<-'simulated';bins[[rep+1L]]<-bin
  cat(n,'replicate',rep,'missing',mean(is.na(sim$matrix)),'\n');flush.console()
 }
 if(digest::digest(file=preset_path,algo='sha256')!=preset_hash)stop('Preset changed during the fidelity gate.')
 d<-do.call(rbind,replicates);all_reps[[n]]<-d;all_bins[[n]]<-do.call(rbind,bins)
 average<-do.call(rbind,lapply(observed$tier,function(t){j<-d$tier==t;data.frame(tier=t,missing_fraction=mean(d$missing_fraction[j]),mean_proteins_per_cell=mean(d$mean_proteins_per_cell[j]),cell_sd=mean(d$cell_sd[j],na.rm=TRUE),batch_sd=if(lf)NA_real_ else mean(d$batch_sd[j],na.rm=TRUE))}))
 checks<-list();for(metric in c('missing_fraction','mean_proteins_per_cell','cell_sd','batch_sd')) {
  tt<-if(metric=='mean_proteins_per_cell')0L else observed$tier
  for(t in tt){idx<-if(t==0)1L else match(t,observed$tier);actual<-observed[[metric]][idx];predicted<-average[[metric]][idx]
   tolerance<-switch(metric,missing_fraction=.05,mean_proteins_per_cell=.1,cell_sd=.2,batch_sd=.25)
   applicable<-!(metric=='batch_sd'&&lf)
   error<-if(metric=='missing_fraction')abs(predicted-actual) else if(is.finite(actual)&&actual>0)abs(predicted/actual-1) else if(is.finite(actual)&&actual==0&&predicted==0)0 else NA_real_
   passed<-if(!applicable)NA else is.finite(error)&&error<=tolerance
   checks[[length(checks)+1L]]<-data.frame(dataset=n,tier=t,metric=metric,observed=actual,simulated_mean=predicted,error=error,tolerance=tolerance,error_units=if(metric=='missing_fraction')'absolute fraction' else 'relative fraction',applicable=applicable,pass=passed,reason=if(!applicable)'Not applicable: run variance absorbed into sigma2_e under F' else if(passed)'within registered tolerance' else 'FAIL: unchanged registered model/preset')
  }
 }
 checks<-do.call(rbind,checks);all_checks[[n]]<-checks
 all_half[[n]]<-data.frame(dataset=n,observed=observed$half_cells_fraction[1],simulated_mean=mean(d$half_cells_fraction[d$tier==observed$tier[1]]),simulated_sd=sd(d$half_cells_fraction[d$tier==observed$tier[1]]),gate=FALSE)
 passed<-all(checks$pass[checks$applicable]);summary<-data.frame(dataset=n,status=if(passed)'PASS' else 'FAIL',checks_passed=sum(checks$pass,na.rm=TRUE),checks_applicable=sum(checks$applicable),maximum_missingness_error=max(checks$error[checks$metric=='missing_fraction']),coverage_relative_error=checks$error[checks$metric=='mean_proteins_per_cell'],maximum_cell_sd_error=max(checks$error[checks$metric=='cell_sd']),maximum_batch_sd_error=if(lf)NA_real_ else max(checks$error[checks$metric=='batch_sd']),proteins=nrow(om),cells=ncol(om),samples=length(samples),batches=length(batches),replicates=20)
 all_summary[[n]]<-summary
 manifests[[n]]<-list(seed=seed,versions=as.list(vs),native_n_per_condition=ncounts,native_cells_per_sample=cellcounts,input_hashes=list(observed_missingness=digest::digest(om,algo='sha256'),observed_variance=digest::digest(ov,algo='sha256'),metadata=digest::digest(md,algo='sha256'),layout=digest::digest(layout,algo='sha256'),preset=preset_hash),replicates=run_records,status=summary$status,wall_seconds=proc.time()[['elapsed']]-dataset_started)
 cat(n,summary$status,'\n');flush.console()
}
write_csv(do.call(rbind,all_half),'half-cells-detection.csv');write_csv(do.call(rbind,all_checks),'fidelity-table.csv');write_csv(do.call(rbind,all_summary),'fidelity-summary.csv');write_csv(do.call(rbind,all_reps),'replicate-metrics.csv');write_csv(do.call(rbind,all_observed),'observed-metrics.csv');write_csv(do.call(rbind,all_bins),'missingness-intensity-bins.csv')
jsonlite::write_json(list(mandate_addendum='1.1.5',seed=seed,versions=as.list(vs),package_version='0.3.0.9000',registered_replicates=20,tolerances=list(missingness_absolute=.05,coverage_relative=.1,cell_sd_relative=.2,batch_sd_relative=.25),metric_definitions=list(cell_sd='Median of measured per-protein/per-line SDs after cell median shifts, within each retained tier',batch_sd='Median across proteins of SD of mixed-batch means after subtracting observed within-line protein means; same statistic applied to simulations',missingness='Endpoint aggregated quantitative mask for TMT under I, fixed accepted cell/protein population',aggregation='Compare means over twenty independent simulations; every applicable tier must meet its tolerance',native_design='Cell lines are samples. Two condition sample counts can be unequal or include a zero count at beta=0. Exact source sample-by-batch occupancy retained.'),datasets=manifests,input_hashes=list(final_calibration=digest::digest(file='reports/M1-final/calibration-results.rds',algo='sha256'),simulator=digest::digest(file='R/scp_simulate.R',algo='sha256'),gate_script=digest::digest(file='inst/scripts/fidelity_m2b.R',algo='sha256'),renv_lock=digest::digest(file='renv.lock',algo='sha256')),wall_seconds=proc.time()[['elapsed']]-started,no_retuning=TRUE,between_donor_transfer_not_validated=TRUE,M4_registration_unchanged=TRUE),file.path(out,'run-manifest.json'),pretty=TRUE,auto_unbox=TRUE,na='null',digits=NA)
cat('Registered fidelity gate completed for every dataset; failures preserved.\n')
