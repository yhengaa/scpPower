# Evaluation helpers only; no simulator or calibration model changes.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
source('R/scp_calibrate.R');source('R/scp_simulate.R');source('R/scp_power.R')
m4_out <- 'reports/M4b'
m4_versions <- .scp_versions()
m4_vtext <- paste(names(m4_versions),m4_versions,sep='=',collapse='; ')
m4_hash <- function(p)digest::digest(file=p,algo='sha256')
m4_json <- function(x,p)jsonlite::write_json(x,p,pretty=TRUE,auto_unbox=TRUE,na='null',digits=NA)
m4_csv <- function(x,p){x$software_versions<-rep(m4_vtext,nrow(x));write.csv(x,p,row.names=FALSE,na='NA')}
m4_assignment <- function(cells,n,m,seed){
 RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(seed)
 chosen<-sample.int(nrow(cells),2L*n*m,replace=FALSE)
 data.frame(cell=cells$cell[chosen],sample=rep(paste0('S',seq_len(2L*n)),each=m),condition=rep(c('A','B'),each=n*m),batch=cells$batch[chosen],stringsAsFactors=FALSE)
}
m4_layout <- function(md,n){
 z<-unclass(table(factor(md$sample,levels=paste0('S',seq_len(2L*n))),factor(md$batch,levels=unique(md$batch))))
 storage.mode(z)<-'double';z
}
m4_map_bulk <- function(ids,bulk){
 ix<-match(ids,bulk$protein);v<-bulk$sigma2_s[ix]
 v[grepl(';',ids)]<-NA_real_;v
}
m4_fit_half <- function(dataset){
 start<-proc.time()[['elapsed']]
 rs<-readRDS('reports/M1-final/calibration-results.rds');old<-rs[[dataset]]
 obs<-readRDS(paste0('.local/datasets/',dataset,'-aggregated.rds'))
 tmt<-dataset%in%c('specht2019v3','leduc2022_pSCoPE')
 raw<-if(tmt)readRDS(paste0('.local/datasets/',dataset,'-psm-proteins.rds')) else obs
 split<-read.csv(file.path(m4_out,paste0(dataset,'-protein-split.csv')))
 ids<-split$protein[split$half=='calibration'];val<-split[split$half=='validation',]
 z<-obs$matrix[ids,,drop=FALSE];mask_z<-raw$matrix[ids,colnames(z),drop=FALSE]
 md<-obs$cell_metadata
 # Protein-held-out calibration uses only calibration rows, over the native
 # accepted calibration cell population. H1 restricts validation source cells.
 fresh<-scp_calibrate(z,md,sample_col=old$manifest$sample_col,
  condition_col=old$manifest$sample_col,batch_col=old$manifest$batch_col,
  biological_unit='culture',normalize_cells=TRUE,missingness_input=mask_z,
  absorb_cell_batches=dataset=='brunner2022',seed=split$split_seed[1],
  reference_versions=unlist(old$manifest$versions))
 # The bulk variance is external to both SCP protein halves. No SCP value
 # from validation proteins enters the calibration or any fallback estimate.
 bulk<-old$bulk_proteins;fresh$bulk_proteins<-bulk
 fresh$proteins$sigma2_s<-m4_map_bulk(ids,bulk)
 fresh$proteins$sigma2_s_reason<-ifelse(is.finite(fresh$proteins$sigma2_s),'External Walzik mapped variance','Unmapped: calibration-half tier fallback')
 raw_mu<-rowMeans(mask_z,na.rm=TRUE)
 for(t in fresh$tiers$tier){j<-which(fresh$proteins$tier==t)
  fresh$tiers$sigma2_s[fresh$tiers$tier==t]<-mean(fresh$proteins$sigma2_s[j],na.rm=TRUE)
  fresh$tiers$mean_log2[fresh$tiers$tier==t]<-mean(raw_mu[j],na.rm=TRUE)
  fresh$tiers$missing_fraction[fresh$tiers$tier==t]<-mean(is.na(z[j,,drop=FALSE]))
 }
 fresh$manifest$dataset<-dataset;fresh$manifest$mandate_addendum<-'1.1.5'
 config<-.scp_simulation_config(fresh,NULL);pp<-config$protein_parameters
 qz<-z[pp$protein,,drop=FALSE] # I: detection mask of the endpoint assay;mask<-is.finite(qz)
 pp$r_g<-rowMeans(mask)
 if(tmt){
  runs<-split(seq_len(ncol(qz)),as.character(md[[old$manifest$batch_col]]))
  identified<-vapply(runs,function(j)rowSums(mask[,j,drop=FALSE])>0,logical(nrow(qz)))
  den<-as.vector(identified%*%as.integer(lengths(runs)))
  pp$p_g<-rowMeans(identified);pp$q_g<-ifelse(den>0,rowSums(mask)/den,0)
  pp$identified_plexes<-rowSums(identified);pp$observed_plexes<-length(runs)
  config$detection_mode<-'tmt_two_level'
 }else config$detection_mode<-'protein_rate'
 release<-.scp_simulation_config(NULL,paste0('inst/presets/',dataset,'.yaml'))
 for(field in c('plex_size','cells_per_run','minutes_per_run','minutes_per_run_source','carrier_ratio','carrier_ratio_source','k_status','k_sensitivity'))config[[field]]<-release[[field]]
 config$schema_version<-2L;config$k<-1;config$dataset<-dataset
 config$versions<-as.list(m4_versions)
 # Preserve the accepted source cell transformation. Its shift is estimated
 # using calibration proteins only and applied to validation rows before any
 # injected signal; injected donor/condition effects are never normalised out.
 med<-apply(z,2,median,na.rm=TRUE);shifts<-med-median(med)
 stopifnot(all(is.finite(shifts)))
 shift_df<-data.frame(cell=colnames(z),log2_shift=shifts)
 m4_csv(shift_df,file.path(m4_out,paste0(dataset,'-calibration-cell-shifts.csv')))
 val$sigma2_s<-m4_map_bulk(val$protein,bulk)
 val$sigma2_s_tier_fallback<-!is.finite(val$sigma2_s)
 val$sigma2_s[val$sigma2_s_tier_fallback]<-config$tiers$sigma2_s[match(val$released_tier[val$sigma2_s_tier_fallback],config$tiers$tier)]
 stopifnot(all(is.finite(val$sigma2_s)),all(val$sigma2_s>=0),!any(pp$protein%in%val$protein))
 m4_csv(val,file.path(m4_out,paste0(dataset,'-validation-donor-variances.csv')))
 config$provenance<-list(mandate='1.1.5',seed=split$split_seed[1],versions=as.list(m4_versions),
  calibration_proteins=ids,calibration_cell_population='Native accepted M1 calibration cells, calibration protein rows only',
  detection_mask_source='Endpoint aggregated assay under I; raw PSM means and curve retained',
  validation_source='H1 largest line only; no validation protein intensities used to fit this preset',
  validation_tier_assignment='Frozen M1 tier labels used only to select calibration-half donor fallback for unmapped validation proteins',
  bulk_proxy='PXD058573 whole PBMC; injected variance cannot validate transfer',
  variance_normalization='Per-cell median shift estimated from calibration half only',
  input_hashes=list(protein_split=m4_hash(file.path(m4_out,paste0(dataset,'-protein-split.csv'))),
   quantitative_assay=m4_hash(paste0('.local/datasets/',dataset,'-aggregated.rds')),
   missingness_assay=m4_hash(paste0('.local/datasets/',dataset,if(tmt)'-psm-proteins.rds' else '-aggregated.rds')),
   external_bulk=digest::digest(bulk,algo='sha256')),
  wall_seconds=proc.time()[['elapsed']]-start)
 m4_csv(pp,file.path(m4_out,paste0(dataset,'-calibration-protein-parameters.csv')))
 # Serialise full precision; the released preset files are never written.
 config$protein_parameters<-lapply(seq_len(nrow(pp)),function(i)as.list(pp[i,,drop=FALSE]))
 config$tiers<-lapply(seq_len(nrow(config$tiers)),function(i)as.list(config$tiers[i,,drop=FALSE]))
 path<-paste0('inst/presets/m4b/',dataset,'.yaml');yaml::write_yaml(config,path,precision=17)
 m4_csv(fresh$proteins,file.path(m4_out,paste0(dataset,'-half-calibration.csv')))
 m4_csv(fresh$tiers,file.path(m4_out,paste0(dataset,'-half-tier-calibration.csv')))
 m4_json(list(seed=split$split_seed[1],versions=as.list(m4_versions),input_hashes=config$provenance$input_hashes,
  outputs=list(preset=m4_hash(path),donor_map=m4_hash(file.path(m4_out,paste0(dataset,'-validation-donor-variances.csv')))),
  calibration_proteins=length(ids),validation_proteins=nrow(val),generating_proteins=nrow(pp),
  excluded_undefined_tier_calibration_proteins=sum(is.na(fresh$proteins$tier)),
  calibration=fresh$manifest,wall_seconds=proc.time()[['elapsed']]-start),file.path(m4_out,paste0(dataset,'-half-preset-manifest.json')))
 cat('Calibration half fixed:',dataset,nrow(pp),'proteins\n');flush.console()
 dataset
}
