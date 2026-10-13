# Amendment G: regenerate complete empirical protein tables from frozen M1.
# No calibration, tier assignment or variance estimates are changed here.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
source('R/scp_calibrate.R');source('R/scp_simulate.R')
rs<-readRDS('reports/M1-final/calibration-results.rds')
b<-jsonlite::read_json('reports/M1-final/paper-benchmarks.json',simplifyVector=TRUE)$benchmarks
started<-proc.time()[['elapsed']];out<-'reports/M2-revised'
dir.create(out,recursive=TRUE,showWarnings=FALSE);records<-list()
for(n in names(rs)) {
 r<-rs[[n]];config<-.scp_simulation_config(r,NULL)
 obs<-readRDS(paste0('.local/datasets/',n,'-aggregated.rds'))
 md<-obs$cell_metadata;batch<-as.character(md[[r$manifest$batch_col]])
 tmt<-n%in%c('specht2019v3','leduc2022_pSCoPE')
 raw<-if(tmt)readRDS(paste0('.local/datasets/',n,'-psm-proteins.rds')) else obs
 z<-raw$matrix[config$protein_parameters$protein,colnames(obs$matrix),drop=FALSE]
 mask<-is.finite(z);mu<-rowMeans(z,na.rm=TRUE)
 pp<-config$protein_parameters
 pp$mean_fallback<-!is.finite(mu)
 pp$mean_log2<-mu
 pp$mean_log2[pp$mean_fallback]<-config$tiers$mean_log2[match(pp$tier[pp$mean_fallback],config$tiers$tier)]
 pp$r_g<-rowMeans(mask)
 if(tmt) {
  runs<-split(seq_len(ncol(z)),batch)
  identified<-vapply(runs,function(j)rowSums(mask[,j,drop=FALSE])>0,logical(nrow(z)))
  denominator<-as.vector(identified%*%as.integer(lengths(runs)))
  pp$p_g<-rowMeans(identified)
  pp$q_g<-ifelse(denominator>0,rowSums(mask)/denominator,0)
  pp$identified_plexes<-rowSums(identified);pp$observed_plexes<-length(runs)
  pp$channels_in_identified_plexes<-denominator
  config$detection_mode<-'tmt_two_level'
 } else config$detection_mode<-'protein_rate'
 config$protein_parameters<-lapply(seq_len(nrow(pp)),function(i)as.list(pp[i,,drop=FALSE]))
 config$schema_version<-2L;config$dataset<-n
 config$plex_size<-if(config$label_free)1L else as.integer(max(table(batch)))
 config$carrier_ratio<-switch(n,specht2019v3=b$specht$carrier_ratio,leduc2022_pSCoPE=b$leduc$carrier_ratio,NULL)
 config$carrier_ratio_source<-if(n=='specht2019v3')b$specht$source else if(n=='leduc2022_pSCoPE')b$leduc$source else 'No isobaric carrier'
 config$k<-1;config$k_status<-'Registered user assumption, not estimated';config$k_sensitivity<-c(.75,.5)
 config$provenance<-list(mandate='1.1.3',calibration_version='0.1.0',calibration_source_sha256=digest::digest(file='reports/M1-final/calibration-results.rds',algo='sha256'),variance_scale=r$manifest$variance_scale,missingness_scale=r$manifest$missingness_scale,bulk_proxy='PXD058573 whole PBMC; donor-transfer assumption not validated',flags=r$manifest$flags,seed=20261013L,versions=as.list(.scp_versions()),detection_population='Unchanged accepted M1 cells and proteins, grouped by native batch; finite raw quantities indicate detection',mean_fallback='Never-quantified proteins retain zero detection and use their tier mean only as an undefined-intensity fallback',G2_consequence='For carrier TMT, p_g and q_g are effect-independent because carrier pools both conditions; testability depends on plex count and barely on effect size')
 columns<-c('tier','proteins','mean_log2','sigma2_s','sigma2_b','sigma2_e','missing_fraction')
 config$tiers<-lapply(seq_len(nrow(config$tiers)),function(i)as.list(config$tiers[i,columns,drop=FALSE]))
 path<-paste0('inst/presets/',n,'.yaml');yaml::write_yaml(config,path,precision=17)
 csv<-pp;csv$software_versions<-paste(names(.scp_versions()),.scp_versions(),sep='=',collapse='; ')
 write.csv(csv,file.path(out,paste0(n,'-protein-parameters.csv')),row.names=FALSE,na='NA')
 records[[n]]<-list(path=path,sha256=digest::digest(file=path,algo='sha256'),mask_sha256=digest::digest(mask,algo='sha256'),metadata_sha256=digest::digest(md,algo='sha256'),seed=20261013L,versions=as.list(.scp_versions()),mean_fallback_count=sum(pp$mean_fallback),component_fallback_counts=as.list(colSums(pp[,grep('_tier_fallback$',names(pp)),drop=FALSE])))
}
jsonlite::write_json(list(seed=20261013L,versions=as.list(.scp_versions()),input_hashes=list(final_calibration=digest::digest(file='reports/M1-final/calibration-results.rds',algo='sha256'),paper_benchmarks=digest::digest(file='reports/M1-final/paper-benchmarks.json',algo='sha256'),generator=digest::digest(file='inst/scripts/generate_presets_revised.R',algo='sha256')),presets=records,wall_seconds=proc.time()[['elapsed']]-started),file.path(out,'preset-manifest.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
cat('Generated five per-protein presets under G; calibration unchanged.\n')
