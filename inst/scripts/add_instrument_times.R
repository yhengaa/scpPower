# Append cited acquisition metadata after M2 passes; numerical generators stay fixed.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
source('R/scp_calibrate.R');source('R/scp_simulate.R')
start<-proc.time()[['elapsed']];out<-'reports/M3'
ss<-read.csv('reports/M2-revised/fidelity-summary.csv')
if(nrow(ss)!=5L||any(ss$status!='PASS'))stop('All five M2 series must pass before M3.')
time<-jsonlite::read_json(file.path(out,'instrument-times.json'),simplifyVector=TRUE)$sources
records<-list()
for(n in ss$dataset){
 path<-paste0('inst/presets/',n,'.yaml');before<-digest::digest(file=path,algo='sha256')
 config<-yaml::read_yaml(path,eval.expr=FALSE)
 tmt<-n%in%c('specht2019v3','leduc2022_pSCoPE')
 if(tmt){
  raw<-readRDS(paste0('.local/datasets/',n,'-raw.rds'))$global_metadata
  md<-readRDS(paste0('.local/datasets/',n,'-aggregated.rds'))$cell_metadata
  types<-unique(md$SampleType)
  cap<-as.integer(max(table(raw$Set[raw$SampleType%in%types])))
  capacity_hash<-digest::digest(raw,algo='sha256')
 }else{cap<-if(n=='brunner2022')1L else 3L;capacity_hash<-NA_character_}
 key<-switch(n,specht2019v3='specht',leduc2022_pSCoPE='leduc',derks2022_timsTOF='derks_tims',derks2022_QExactive='derks_qe',brunner2022='brunner')
 config$cells_per_run<-cap;config$minutes_per_run<-time[[key]]$minutes_per_run
 config$minutes_per_run_source<-time[[key]]
 config$cells_per_run_source<-if(tmt)'Maximum native single-cell channel count in raw scpdata SampleType/Set metadata, excluding carrier/reference/empty/control' else if(cap==3L)'Derks 2022 plexDIA three-plex single-cell design' else 'One label-free cell per run under F'
 config$provenance$M2_gate_preset_sha256<-before
 yaml::write_yaml(config,path,precision=17)
 records[[n]]<-list(seed=20261013,versions=as.list(.scp_versions()),before_sha256=before,after_sha256=digest::digest(file=path,algo='sha256'),capacity_metadata_sha256=capacity_hash,cells_per_run=cap,minutes_per_run=config$minutes_per_run,source=time[[key]])
}
jsonlite::write_json(list(seed=20261013,versions=as.list(.scp_versions()),input_hashes=list(times=digest::digest(file=file.path(out,'instrument-times.json'),algo='sha256'),passed_gate=digest::digest(file='reports/M2-revised/fidelity-summary.csv',algo='sha256')),presets=records,generation_parameters_unchanged=TRUE,wall_seconds=proc.time()[['elapsed']]-start),file.path(out,'instrument-preset-manifest.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA,na='null')
