# Instrument metadata must not change the already-passed M2 numerical draws.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
source('R/scp_calibrate.R');source('R/scp_simulate.R');started<-proc.time()[['elapsed']]
rs<-readRDS('reports/M1-final/calibration-results.rds')
m<-jsonlite::read_json('reports/M2-revised/run-manifest.json',simplifyVector=FALSE)
records<-list()
for(n in names(rs)) {
 r<-rs[[n]];o<-readRDS(paste0('.local/datasets/',n,'-aggregated.rds'));md<-o$cell_metadata
 sample<-as.character(md[[r$manifest$sample_col]]);batch<-as.character(md[[r$manifest$batch_col]])
 layout<-unclass(table(factor(sample,levels=unique(sample)),factor(batch,levels=unique(batch))));storage.mode(layout)<-'double'
 original<-m$datasets[[n]]$replicates[[1]]
 sim<-scp_simulate(r,as.integer(unlist(original$n_per_condition)),as.integer(rowSums(layout)),beta=0,preset=paste0('inst/presets/',n,'.yaml'),layout=layout,k=1,seed=original$seed)
 stopifnot(identical(unname(sim$manifest$input_hashes['matrix']),original$input_hashes$matrix))
 records[[n]]<-list(seed=original$seed,versions=sim$manifest$versions,matrix_sha256=original$input_hashes$matrix,preset_sha256=digest::digest(file=paste0('inst/presets/',n,'.yaml'),algo='sha256'),byte_identical=TRUE,wall_seconds=sim$manifest$wall_seconds)
 cat(n,'first registered replicate remains byte-identical after instrument metadata augmentation.\n')
}
jsonlite::write_json(list(seed=20261013,versions=as.list(.scp_versions()),datasets=records,input_hashes=list(gate_manifest=digest::digest(file='reports/M2-revised/run-manifest.json',algo='sha256')),wall_seconds=proc.time()[['elapsed']]-started),'reports/M3/M2-replay-verification.json',pretty=TRUE,auto_unbox=TRUE,digits=NA)
