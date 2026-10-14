# Worked examples, separate from the registered M3/M3b grid.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
source('R/scp_calibrate.R');source('R/scp_simulate.R');source('R/scp_power.R')
start<-proc.time()[['elapsed']]
names<-c('leduc2022_pSCoPE','derks2022_timsTOF','brunner2022')
results<-parallel::mclapply(names,function(name){
 p<-paste0('inst/presets/',name,'.yaml')
 x<-scp_power(p,n_per_condition=c(8,12),m_cells_per_sample=32,beta=log2(2),replicates=200,seed=20261013)
 x$results$dataset<-name
 write.csv(x$results,paste0('reports/M5/',name,'-worked-example.csv'),row.names=FALSE)
 jsonlite::write_json(x$manifest,paste0('reports/M5/',name,'-worked-example-manifest.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
 x$results
},mc.cores=3)
stopifnot(!any(vapply(results,inherits,logical(1),'try-error')))
write.csv(do.call(rbind,results),'reports/M5/worked-examples.csv',row.names=FALSE)
jsonlite::write_json(list(seed=20261013,versions=as.list(.scp_versions()),registered_non_TMT_grid_not_rerun=TRUE,demonstration_runs_only=TRUE,input_hashes=setNames(lapply(names,function(n)digest::digest(file=paste0('inst/presets/',n,'.yaml'),algo='sha256')),names),wall_seconds=proc.time()[['elapsed']]-start,cluster=list(job=Sys.getenv('SLURM_JOB_ID'),cores=Sys.getenv('SLURM_CPUS_PER_TASK'))),'reports/M5/worked-examples-manifest.json',pretty=TRUE,auto_unbox=TRUE,digits=NA)
