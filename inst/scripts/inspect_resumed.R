# Download the registered quantitative assays through scpdata/ExperimentHub.
if (dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'), .libPaths()))
options(timeout=1800)
source('R/scp_calibrate.R')
dir.create('reports/M1-resumed', recursive=TRUE, showWarnings=FALSE)
records <- list()
for (name in c('specht2019v3','leduc2022_pSCoPE','derks2022','brunner2022')) {
 t0 <- proc.time()[['elapsed']]
 cat('Loading',name,'\n'); flush.console()
 x <- getExportedValue('scpdata',name)()
 assays <- switch(name, specht2019v3='peptides', leduc2022_pSCoPE='peptides', derks2022=c('tims_prec_extracted','qe_prec_extracted'), brunner2022='proteins')
 object <- list(global_metadata=as.data.frame(SummarizedExperiment::colData(x)), assays=setNames(lapply(assays,function(a) list(matrix=as.matrix(SummarizedExperiment::assay(x[[a]])), cell_metadata=as.data.frame(SummarizedExperiment::colData(x[[a]])), row_metadata=as.data.frame(SummarizedExperiment::rowData(x[[a]])))),assays))
 path <- paste0('.local/datasets/',name,'-raw.rds'); saveRDS(object,path)
 records[[name]] <- list(seed=20261013L,versions=as.list(.scp_versions()),scpdata_version=as.character(packageVersion('scpdata')),object_hash=digest::digest(x,algo='sha256'),extracted_file_hash=digest::digest(file=path,algo='sha256'),assays=assays,wall_seconds=proc.time()[['elapsed']]-t0)
 for (a in assays) {z <- object$assays[[a]]; cat(name,a,dim(z$matrix),'missing',mean(is.na(z$matrix)),'range',range(z$matrix,na.rm=TRUE),'\n');print(head(names(z$row_metadata),10));print(z$row_metadata[1:2,1:min(5,ncol(z$row_metadata)),drop=FALSE]);print(names(z$cell_metadata));print(head(z$cell_metadata,2)); print(table(object$global_metadata$SampleType));print(table(object$global_metadata$Celltype))}
 jsonlite::write_json(records,'reports/M1-resumed/input-inspection.json',pretty=TRUE,auto_unbox=TRUE,digits=NA)
 rm(x,object);gc()
}
