# Apply the conditional decision recorded in Diagnostic J, changing mask only.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
source('R/scp_calibrate.R');source('R/scp_simulate.R')
start<-proc.time()[['elapsed']];out<-'reports/M2b';dir.create(out,recursive=TRUE,showWarnings=FALSE)
j<-read.csv('reports/J/decision-table.csv');stopifnot(all(j$amendment_I_applies))
vs<-.scp_versions();vt<-paste(names(vs),vs,sep='=',collapse='; ');records<-list()
for(d in c('specht2019v3','leduc2022_pSCoPE')){
 t0<-proc.time()[['elapsed']];path<-paste0('inst/presets/',d,'.yaml')
 baseline<-file.path(out,paste0(d,'-preset-v0.3.0.yaml'))
 if(!file.exists(baseline))file.copy(path,baseline)
 old<-.scp_simulation_config(NULL,baseline)
 before<-digest::digest(file=baseline,algo='sha256');config<-old
 obs<-readRDS(paste0('.local/datasets/',d,'-aggregated.rds'));rs<-readRDS('reports/M1-final/calibration-results.rds')[[d]]
 pp<-config$protein_parameters;mask<-is.finite(obs$matrix[pp$protein,,drop=FALSE]);batch<-as.character(obs$cell_metadata[[rs$manifest$batch_col]])
 runs<-split(seq_len(ncol(mask)),batch);identified<-vapply(runs,function(k)rowSums(mask[,k,drop=FALSE])>0,logical(nrow(mask)))
 den<-as.vector(identified%*%as.integer(lengths(runs)))
 pp$p_g<-rowMeans(identified);pp$q_g<-ifelse(den>0,rowSums(mask)/den,0);pp$r_g<-rowMeans(mask)
 pp$identified_plexes<-rowSums(identified);pp$observed_plexes<-length(runs);pp$channels_in_identified_plexes<-den
 for(t in config$tiers$tier){ix<-pp$tier==t;config$tiers$missing_fraction[config$tiers$tier==t]<-mean(!mask[ix,,drop=FALSE])}
 # No generating intensity or component estimate is changed under I.
 fields<-c('protein','tier','mean_log2','sigma2_s','sigma2_b','sigma2_e')
 stopifnot(identical(pp[,fields],old$protein_parameters[,fields]),identical(config$curve,old$curve))
 config$protein_parameters<-lapply(seq_len(nrow(pp)),function(i)as.list(pp[i,,drop=FALSE]))
 config$tiers<-lapply(seq_len(nrow(config$tiers)),function(i)as.list(config$tiers[i,,drop=FALSE]))
 config$provenance$mandate<-'1.1.5';config$provenance$detection_mask_source<-'Finite ratio-based aggregated assay used by the primary endpoint; amendment I'
 config$provenance$diagnostic_J_sha256<-digest::digest(file='reports/J/diagnostic-report.md',algo='sha256')
 config$provenance$pre_I_preset_sha256<-before
 yaml::write_yaml(config,path,precision=17)
 pp$software_versions<-vt;write.csv(pp,file.path(out,paste0(d,'-protein-parameters.csv')),row.names=FALSE,na='NA')
 records[[d]]<-list(seed=20261013L,versions=as.list(vs),input_hashes=list(endpoint_file=digest::digest(file=paste0('.local/datasets/',d,'-aggregated.rds'),algo='sha256'),pre_I_preset=before,diagnostic_J=digest::digest(file='reports/J/diagnostic-report.md',algo='sha256')),output_preset_sha256=digest::digest(file=path,algo='sha256'),unchanged_fields=fields,unchanged_curve=TRUE,wall_seconds=proc.time()[['elapsed']]-t0)
}
jsonlite::write_json(list(seed=20261013L,versions=as.list(vs),amendment_I=TRUE,datasets=records,wall_seconds=proc.time()[['elapsed']]-start),file.path(out,'preset-manifest.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
cat('Amendment I applied to two TMT masks only, after Diagnostic J report.\n')
