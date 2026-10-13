# Execute the same registered points within an interactive compute allocation.
# Numerical inputs and seeds are identical to the batch-array worker.
if(!nzchar(Sys.getenv('SLURM_JOB_ID')))stop('This launcher requires a compute allocation.')
cores<-as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','1'))
paths<-file.path('reports/M3/points',sprintf('%04d-manifest.json',seq_len(945L)))
pending<-which(!file.exists(paths));cat('Registered points remaining:',length(pending),'cores:',cores,'\n');flush.console()
status<-parallel::mclapply(pending,function(i){
 Sys.setenv(SCPPOWER_GRID_TASK=as.character(i))
 source('inst/scripts/run_m3_grid.R',local=new.env(parent=globalenv()))
 i
},mc.cores=cores,mc.preschedule=FALSE)
if(any(vapply(status,inherits,logical(1),'try-error'))||any(!file.exists(paths)))stop('Some registered grid points did not complete; preserve and report the failed workers.')
cat('All registered points now complete.\n')
