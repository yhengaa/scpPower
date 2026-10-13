# Name only pipeline losses traceable from the deposited peptide assay.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
source('R/scp_calibrate.R')
start<-proc.time()[['elapsed']];vs<-.scp_versions();vt<-paste(names(vs),vs,sep='=',collapse='; ');out<-'reports/J';rows<-list();records<-list()
for(d in c('specht2019v3','leduc2022_pSCoPE')){
 t0<-proc.time()[['elapsed']];raw<-readRDS(paste0('.local/datasets/',d,'-raw.rds'))$assays$peptides
 agg<-readRDS(paste0('.local/datasets/',d,'-aggregated.rds'));psm<-readRDS(paste0('.local/datasets/',d,'-psm-proteins.rds'))
 split<-read.csv(paste0('reports/M4/',d,'-protein-split.csv'));ids<-split$protein[split$half=='validation'];cells<-read.csv(paste0('reports/M4/',d,'-source-cells.csv'))$cell
 z<-raw$matrix;rd<-raw$row_metadata
 if(d=='specht2019v3'){mapping<-as.character(rd$protein);feature_log<-z}else{
  bad<-grepl('REV__|CON__',rd$Leading.razor.protein);z<-z[!bad,,drop=FALSE];mapping<-as.character(rd$Leading.razor.protein.id[!bad]);feature_log<-.scp_matrix(z,FALSE,'deposited peptides')
 }
 reproduced<-.scp_aggregate(feature_log,mapping)[ids,cells,drop=FALSE]
 endpoint<-agg$matrix[ids,cells,drop=FALSE];stopifnot(identical(is.finite(reproduced),is.finite(endpoint)))
 upstream_stage<-if(d=='specht2019v3')'Deposited Peptides-raw.csv after producer reference normalisation, quality/cell and missing-data peptide filtering; exported from cr_norm_log(t3)' else 'Deposited processed peptide assay imported by scpdata; reference-ratio quantitative values and retained peptide/protein-group mapping'
 rows[[d]]<-data.frame(dataset=d,proteins=length(ids),cells=length(cells),psm_finite_endpoint_NA_pairs=sum(is.finite(psm$matrix[ids,cells,drop=FALSE]) & !is.finite(endpoint)),endpoint_mask_reproduced_from_peptides=TRUE,
  newly_lost_finite_pairs_at_local_protein_median=0L,paired_cells_removed_between_masks=0L,stage=upstream_stage,
  attribution_limit='For a PSM-positive/endpoint-NA pair every included deposited peptide quantity is unavailable. The local median does not discard any finite included peptide. Exact attribution to a particular upstream filter requires per-PSM retention lineage, not available in these two processed objects.',software_versions=vt)
 records[[d]]<-list(seed=20261013L,versions=as.list(vs),input_hashes=list(raw_assay=digest::digest(file=paste0('.local/datasets/',d,'-raw.rds'),algo='sha256'),quantitative_assay=digest::digest(file=paste0('.local/datasets/',d,'-aggregated.rds'),algo='sha256'),psm_assay=digest::digest(file=paste0('.local/datasets/',d,'-psm-proteins.rds'),algo='sha256')),wall_seconds=proc.time()[['elapsed']]-t0)
}
write.csv(do.call(rbind,rows),file.path(out,'pipeline-audit.csv'),row.names=FALSE)
jsonlite::write_json(list(seed=20261013L,versions=as.list(vs),datasets=records,input_hashes=list(script=digest::digest(file='inst/scripts/audit_j_pipeline.R',algo='sha256')),wall_seconds=proc.time()[['elapsed']]-start),file.path(out,'pipeline-audit-manifest.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
