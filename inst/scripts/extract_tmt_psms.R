# Download registered raw PSM assays and pair their cell/protein identifiers
# with the accepted variance population. No new QC threshold is applied.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
options(timeout=1800);source('R/scp_calibrate.R')
records<-list();dir.create('reports/M1-final',recursive=TRUE,showWarnings=FALSE)
for(dataset in c('specht2019v3','leduc2022_pSCoPE')) {
 t0<-proc.time()[['elapsed']];cat('Extracting raw PSMs',dataset,'\n');flush.console()
 qf<-getExportedValue('scpdata',dataset)();nr<-if(dataset=='specht2019v3')177L else 134L
 old<-readRDS(paste0('.local/datasets/',dataset,'-aggregated.rds'))
 global<-as.data.frame(SummarizedExperiment::colData(qf));types<-unique(as.character(old$cell_metadata$SampleType))
 all_sc<-rownames(global)[global$SampleType%in%types]
 z<-matrix(NA_real_,nrow(old$matrix),ncol(old$matrix),dimnames=dimnames(old$matrix))
 seen_cells<-character();seen_proteins<-character();ratios<-numeric();runs<-list();zero_count<-0L
 for(a in names(qf)[seq_len(nr)]) {
  se<-qf[[a]];mat<-as.matrix(SummarizedExperiment::assay(se));md<-global[colnames(mat),,drop=FALSE]
  if(anyNA(md$SampleType))stop('PSM channel identity cannot be established.')
  rd<-as.data.frame(SummarizedExperiment::rowData(se))
  original_ids<-if(dataset=='specht2019v3')as.character(rd$protein) else as.character(rd$Leading.razor.protein)
  ids<-gsub('^.*\\|([^|]+)\\|.*$','\\1',original_ids)
  good<-!is.na(ids)&nzchar(ids)&!grepl('REV__|CON__',original_ids)
  ids<-ids[good];mat<-mat[good,,drop=FALSE]
  raw_sc<-md$SampleType%in%types
  sc<-which(raw_sc)
  cells<-which(raw_sc & colnames(mat)%in%colnames(z))
  seen_proteins<-union(seen_proteins,ids);seen_cells<-union(seen_cells,colnames(mat)[sc])
  if(!length(cells))next
  mat<-mat[,cells,drop=FALSE];zero_count<-zero_count+sum(!is.finite(mat)|mat<=0)
  mat[!is.finite(mat)|mat<=0]<-NA_real_;mat<-log2(mat)
  aggregated<-.scp_aggregate(mat,ids)
  common<-intersect(rownames(z),rownames(aggregated));z[common,colnames(aggregated)]<-aggregated[common,,drop=FALSE]
  runs[[a]]<-list(quantitative_sha256=digest::digest(SummarizedExperiment::assay(se),algo='sha256'),proteins=nrow(aggregated),matched_proteins=length(common),matched_cells=length(cells))
 }
 if(!setequal(colnames(z),intersect(colnames(z),seen_cells)))stop('An accepted cell is absent from raw PSM channels.')
 path<-paste0('.local/datasets/',dataset,'-psm-proteins.rds');saveRDS(list(matrix=z,cell_metadata=old$cell_metadata),path)
 m<-list(dataset=dataset,seed=20261013L,versions=as.list(.scp_versions()),scpdata_version=as.character(packageVersion('scpdata')),psm_assays=nr,raw_single_cell_channels=length(all_sc),paired_accepted_cells=ncol(z),accepted_protein_groups=nrow(z),raw_target_protein_ids=length(seen_proteins),exact_protein_overlap=sum(rownames(z)%in%seen_proteins),never_quantified_proteins=sum(rowSums(!is.na(z))==0),processing=c('SampleType excludes carrier, reference, blank, negative and unused channels','Pair exact cell/protein IDs to the accepted variance population; no new QC threshold','Quantified if any PSM is finite and >0; log2 positive raw intensities; per-run protein median; no normalisation','Reverse/contaminant PSM annotations excluded as non-target identities; ambiguous group IDs never invented'),matched_population_hash=digest::digest(old$matrix,algo='sha256'),serialized_sha256=digest::digest(file=path,algo='sha256'),runs=runs,wall_seconds=proc.time()[['elapsed']]-t0)
 records[[dataset]]<-m;jsonlite::write_json(records,'reports/M1-final/tmt-input-manifest.json',pretty=TRUE,auto_unbox=TRUE,digits=NA,na='null')
 cat(dataset,'raw SC',length(all_sc),'paired',ncol(z),'matched proteins',m$exact_protein_overlap,'missing',mean(is.na(z)),'\n');flush.console()
 rm(qf,z);gc()
}
