# Re-register exactly the original TMT tasks and seeds, with amendment-I inputs.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
source('R/scp_calibrate.R');start<-proc.time()[['elapsed']];out<-'reports/M3b'
grid<-read.csv('reports/M3/registered-grid.csv');grid<-grid[grid$dataset %in% c('specht2019v3','leduc2022_pSCoPE'),]
stopifnot(nrow(grid)==630L,identical(grid$task,seq_len(630L)))
write.csv(grid,file.path(out,'registered-grid.csv'),row.names=FALSE)
hashes<-setNames(lapply(unique(grid$dataset),function(n)digest::digest(file=paste0('inst/presets/',n,'.yaml'),algo='sha256')),unique(grid$dataset))
jsonlite::write_json(list(seed=20261013,versions=as.list(.scp_versions()),replicates_per_point=200,grid_points=630,grid_seed_formula='Original M3: 20261013 + 1000*one-based task; replicate seed = grid seed + replicate index',input_hashes=list(original_grid=digest::digest(file='reports/M3/registered-grid.csv',algo='sha256'),grid=digest::digest(file=file.path(out,'registered-grid.csv'),algo='sha256'),simulator=digest::digest(file='R/scp_simulate.R',algo='sha256'),power=digest::digest(file='R/scp_power.R',algo='sha256'),renv_lock=digest::digest(file='renv.lock',algo='sha256')),preset_hashes=hashes,layout='balanced; actual occupied batches counted',target_power=.8,wall_seconds=proc.time()[['elapsed']]-start,non_TMT_not_rerun=TRUE),file.path(out,'grid-registration.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
