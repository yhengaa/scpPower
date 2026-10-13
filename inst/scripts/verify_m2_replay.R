.libPaths(c(normalizePath('.local/R-library'),.libPaths()));source('R/scp_calibrate.R');source('R/scp_simulate.R')
rs<-readRDS('reports/M1-final/calibration-results.rds');n<-names(rs)[1];r<-rs[[1]];o<-readRDS(paste0('.local/datasets/',n,'-aggregated.rds'));md<-o$cell_metadata
sample<-as.character(md[[r$manifest$sample_col]]);batch<-as.character(md[[r$manifest$batch_col]])
layout<-unclass(table(factor(sample,levels=unique(sample)),factor(batch,levels=unique(batch))));storage.mode(layout)<-'double'
m<-jsonlite::read_json('reports/M2/run-manifest.json',simplifyVector=FALSE);original<-m$datasets[[n]]$replicates[[1]]
sim<-scp_simulate(r,1,as.integer(rowSums(layout)),beta=0,preset=paste0('inst/presets/',n,'.yaml'),layout=layout,k=1,seed=original$seed)
stopifnot(identical(unname(sim$manifest$input_hashes['matrix']),original$input_hashes$matrix))
cat('The first registered replicate is byte-identical after namespace qualification and whitespace formatting.\n')
