# Complete calibration, all predictions, then all observations in a live allocation.
source('inst/scripts/m4_helpers.R');source('inst/scripts/run_m4_point.R')
if(!nzchar(Sys.getenv('SLURM_JOB_ID')))stop('M4 requires a compute allocation.')
start<-proc.time()[['elapsed']];cores<-as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','1'))
grid<-read.csv(file.path(m4_out,'registered-grid.csv'));datasets<-unique(grid$dataset[grid$counted])
handle<-function(expr,dataset,task,stage){
 t0<-proc.time()[['elapsed']]
 tryCatch(expr,error=function(e){
  path<-if(is.na(task))file.path(m4_out,paste0(dataset,'-',stage,'-ERROR.json')) else file.path(m4_out,'points',paste0(sprintf('%04d',task),'-',stage,'-ERROR.json'))
  m4_json(list(dataset=dataset,task=task,stage=stage,seed=20261013L,versions=as.list(m4_versions),error=conditionMessage(e),
   input_hashes=list(registration=m4_hash(file.path(m4_out,'registration.json'))),wall_seconds=proc.time()[['elapsed']]-t0),path)
  cat('Reported error',dataset,stage,task,conditionMessage(e),'\n');flush.console();NA_integer_
 })
}
parallel::mclapply(datasets,function(d)handle(m4_fit_half(d),d,NA_integer_,'calibration'),mc.cores=min(length(datasets),cores),mc.preschedule=FALSE)
tasks<-grid$task[grid$counted & file.exists(paste0('inst/presets/m4/',grid$dataset,'.yaml'))]
parallel::mclapply(tasks,function(i)handle(m4_point(i,'prediction'),grid$dataset[grid$task==i],i,'prediction'),mc.cores=cores,mc.preschedule=FALSE)
# Freeze every successful prediction before executing any observed endpoint.
records<-list()
for(i in tasks){prefix<-file.path(m4_out,'points',sprintf('%04d',i));path<-paste0(prefix,'-prediction-manifest.json')
 if(file.exists(path)){
  m<-jsonlite::read_json(path,simplifyVector=FALSE)
  records[[as.character(i)]]<-list(results=m4_hash(paste0(prefix,'-predicted.csv')),manifest=m4_hash(path),preset=m$input_hashes$preset_file)
 }
}
m4_json(list(frozen_UTC=format(Sys.time(),tz='UTC',usetz=TRUE),seed=20261013L,versions=as.list(m4_versions),
 input_hashes=list(registration=m4_hash(file.path(m4_out,'registration.json'))),prediction_hashes=records,
 counted_points=sum(grid$counted),successful_predictions=length(records),wall_seconds=proc.time()[['elapsed']]-start),file.path(m4_out,'predictions-frozen.json'))
cat('All predictions finished and frozen before first observed endpoint:',length(records),'\n');flush.console()
parallel::mclapply(as.integer(names(records)),function(i)handle(m4_point(i,'observation'),grid$dataset[grid$task==i],i,'observation'),mc.cores=cores,mc.preschedule=FALSE)
m4_json(list(seed=20261013L,versions=as.list(m4_versions),job=Sys.getenv('SLURM_JOB_ID'),host=Sys.info()[['nodename']],cores=cores,
 input_hashes=list(registration=m4_hash(file.path(m4_out,'registration.json')),predictions_frozen=m4_hash(file.path(m4_out,'predictions-frozen.json'))),
 counted_points=sum(grid$counted),successful_predictions=length(records),
 successful_observations=sum(file.exists(file.path(m4_out,'points',paste0(sprintf('%04d',grid$task[grid$counted]),'-observation-manifest.json')))),
 wall_seconds=proc.time()[['elapsed']]-start),file.path(m4_out,'cluster-execution.json'))
cat('M4 cluster stages completed. Errors, if any, are retained without tuning.\n')
