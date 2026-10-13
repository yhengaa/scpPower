# Reproduce one registered YAML preset for each final calibrated series.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
source('R/scp_calibrate.R');source('R/scp_simulate.R')
rs<-readRDS('reports/M1-final/calibration-results.rds')
b<-jsonlite::read_json('reports/M1-final/paper-benchmarks.json',simplifyVector=TRUE)$benchmarks
started<-proc.time()[['elapsed']];dir.create('inst/presets',recursive=TRUE,showWarnings=FALSE);records<-list()
for(n in names(rs)) {
 r<-rs[[n]];config<-.scp_simulation_config(r,NULL)
 md<-readRDS(paste0('.local/datasets/',n,'-aggregated.rds'))$cell_metadata
 batch<-md[[r$manifest$batch_col]]
 config$schema_version<-1L;config$dataset<-n
 config$plex_size<-if(config$label_free)1L else as.integer(max(table(batch)))
 config$carrier_ratio<-switch(n,specht2019v3=b$specht$carrier_ratio,leduc2022_pSCoPE=b$leduc$carrier_ratio,NULL)
 config$carrier_ratio_source<-if(n=='specht2019v3')b$specht$source else if(n=='leduc2022_pSCoPE')b$leduc$source else 'No isobaric carrier'
 config$k<-1;config$k_status<-'Registered user assumption, not estimated';config$k_sensitivity<-c(.75,.5)
 config$provenance<-list(mandate='1.1.2',calibration_version='0.1.0',calibration_source_sha256=digest::digest(file='reports/M1-final/calibration-results.rds',algo='sha256'),variance_scale=r$manifest$variance_scale,missingness_scale=r$manifest$missingness_scale,bulk_proxy='PXD058573 whole PBMC; donor-transfer assumption not validated',flags=r$manifest$flags,seed=20261013L,versions=as.list(.scp_versions()))
 columns<-c('tier','proteins','mean_log2','sigma2_s','sigma2_b','sigma2_e','missing_fraction')
 config$tiers<-lapply(seq_len(nrow(config$tiers)),function(i)as.list(config$tiers[i,columns,drop=FALSE]))
 path<-paste0('inst/presets/',n,'.yaml');yaml::write_yaml(config,path,precision=17)
 records[[n]]<-list(path=path,sha256=digest::digest(file=path,algo='sha256'),seed=20261013L,versions=as.list(.scp_versions()))
}
jsonlite::write_json(list(seed=20261013L,versions=as.list(.scp_versions()),input_hashes=list(final_calibration=digest::digest(file='reports/M1-final/calibration-results.rds',algo='sha256'),paper_benchmarks=digest::digest(file='reports/M1-final/paper-benchmarks.json',algo='sha256')),presets=records,wall_seconds=proc.time()[['elapsed']]-started),'reports/M2/preset-manifest.json',pretty=TRUE,auto_unbox=TRUE,digits=NA)
cat('Generated all five presets from fixed final calibration.\n')
