# Register the full design before any M3 power result is observed.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
source('R/scp_calibrate.R');start<-proc.time()[['elapsed']];out<-'reports/M3'
s<-read.csv('reports/M2-revised/fidelity-summary.csv');stopifnot(nrow(s)==5L,all(s$status=='PASS'))
grid<-do.call(rbind,lapply(s$dataset,function(n){
 k<-if(n%in%c('specht2019v3','leduc2022_pSCoPE'))c(1,.75,.5) else 1
 d<-expand.grid(n_per_condition=c(3,4,6,8,12),m_cells_per_sample=c(1,2,4,8,16,32,64),beta=log2(c(1.25,1.5,2)),k=k,KEEP.OUT.ATTRS=FALSE)
 d$dataset<-n;d
}))
grid$task<-seq_len(nrow(grid));grid$seed<-as.integer(20261013+1000*grid$task)
grid$software_versions<-paste(names(.scp_versions()),.scp_versions(),sep='=',collapse='; ')
write.csv(grid,file.path(out,'registered-grid.csv'),row.names=FALSE)
dir.create(file.path(out,'points'),showWarnings=FALSE);dir.create(file.path(out,'cluster-logs'),showWarnings=FALSE)
hashes<-setNames(lapply(s$dataset,function(n)digest::digest(file=paste0('inst/presets/',n,'.yaml'),algo='sha256')),s$dataset)
jsonlite::write_json(list(mandate_addendum='1.1.3',seed=20261013,versions=as.list(.scp_versions()),replicates_per_point=200,grid_points=nrow(grid),grid_seed_formula='20261013 + 1000*one-based task; replicate seed = grid seed + replicate index',input_hashes=list(gate=digest::digest(file='reports/M2-revised/fidelity-summary.csv',algo='sha256'),grid=digest::digest(file=file.path(out,'registered-grid.csv'),algo='sha256'),simulator=digest::digest(file='R/scp_simulate.R',algo='sha256'),power=digest::digest(file='R/scp_power.R',algo='sha256'),renv_lock=digest::digest(file='renv.lock',algo='sha256')),preset_hashes=hashes,layout='balanced; actual occupied batches counted',target_power=.8,wall_seconds=proc.time()[['elapsed']]-start,M4_registration_unchanged=TRUE),file.path(out,'grid-registration.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
cat(nrow(grid),'registered M3 points; no outcomes have been read.\n')
