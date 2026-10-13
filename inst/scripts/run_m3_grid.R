# One fixed M3 grid point per cluster array task; no parameter fitting/tuning.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
source('R/scp_calibrate.R');source('R/scp_simulate.R');source('R/scp_power.R')
args<-commandArgs(trailingOnly=TRUE);task<-as.integer(Sys.getenv('SCPPOWER_GRID_TASK',if(length(args))args[1] else NA_character_));out<-'reports/M3'
grid<-read.csv(file.path(out,'registered-grid.csv'));g<-grid[task,,drop=FALSE]
if(nrow(g)!=1L||is.na(task))stop('Select one registered grid row.')
path<-paste0('inst/presets/',g$dataset,'.yaml')
expected<-jsonlite::read_json(file.path(out,'grid-registration.json'),simplifyVector=TRUE)$preset_hashes[[g$dataset]]
if(digest::digest(file=path,algo='sha256')!=expected)stop('Preset changed after grid registration.')
start<-proc.time()[['elapsed']]
result<-scp_power(path,g$n_per_condition,g$m_cells_per_sample,g$beta,
 layout='balanced',k=g$k,replicates=200,seed=g$seed,target_power=.8)
result$results$dataset<-g$dataset;result$results$task<-task
result$replicate_metrics$dataset<-g$dataset;result$replicate_metrics$task<-task
prefix<-file.path(out,'points',sprintf('%04d',task))
write.csv(result$results,paste0(prefix,'-results.csv'),row.names=FALSE,na='NA')
con<-gzfile(paste0(prefix,'-replicates.csv.gz'),'wt')
write.csv(result$replicate_metrics,con,row.names=FALSE,na='NA');close(con)
m<-result$manifest
m$grid_task<-task;m$input_hashes$preset_file<-expected
m$output_file_hashes<-list(results=digest::digest(file=paste0(prefix,'-results.csv'),algo='sha256'),replicates=digest::digest(file=paste0(prefix,'-replicates.csv.gz'),algo='sha256'))
m$cluster<-list(job=Sys.getenv('SLURM_ARRAY_JOB_ID',Sys.getenv('SLURM_JOB_ID')),array_task=Sys.getenv('SLURM_ARRAY_TASK_ID'),grid_task=task,host=Sys.info()[['nodename']])
m$worker_wall_seconds<-proc.time()[['elapsed']]-start
jsonlite::write_json(m,paste0(prefix,'-manifest.json'),pretty=TRUE,auto_unbox=TRUE,na='null',digits=NA)
cat('Completed registered task',task,g$dataset,'n',g$n_per_condition,'m',g$m_cells_per_sample,'beta',g$beta,'k',g$k,'\n')
