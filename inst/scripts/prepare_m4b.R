# Freeze eligibility, protein split and first source assignments before outcomes.
source('inst/scripts/m4b_helpers.R')
start<-proc.time()[['elapsed']];seed<-20261013L
for(p in c(m4_out,paste0(m4_out,'/points'),paste0(m4_out,'/assignments'),paste0(m4_out,'/layouts'),paste0(m4_out,'/cluster-logs'),'inst/presets/m4b'))dir.create(p,recursive=TRUE,showWarnings=FALSE)
gate<-read.csv('reports/M2b/fidelity-summary.csv');stopifnot(nrow(gate)==2L,all(gate$status=='PASS'));
rs<-readRDS('reports/M1-final/calibration-results.rds')[c('specht2019v3','leduc2022_pSCoPE')];grids<-list();sources<-list();hashes<-list();task<-0L
for(di in seq_along(rs)){
 dataset<-names(rs)[di];r<-rs[[di]];obs<-readRDS(paste0('.local/datasets/',dataset,'-aggregated.rds'))
 counts<-table(obs$cell_metadata[[r$manifest$sample_col]])
 line<-names(counts)[which.max(counts)];nc<-as.integer(max(counts))
 released<-.scp_simulation_config(NULL,paste0('inst/presets/',dataset,'.yaml'))
 ids<-released$protein_parameters$protein
 split_seed<-as.integer(seed+8000000L+100L*di);set.seed(split_seed)
 cal<-sample(ids,floor(length(ids)/2),replace=FALSE)
 split_df<-data.frame(protein=ids,half=ifelse(ids%in%cal,'calibration','validation'),
  released_tier=released$protein_parameters$tier,split_seed=split_seed)
 m4_csv(split_df,file.path(m4_out,paste0(dataset,'-protein-split.csv')))
 ix<-which(as.character(obs$cell_metadata[[r$manifest$sample_col]])==line)
 cells<-data.frame(cell=colnames(obs$matrix)[ix],batch=as.character(obs$cell_metadata[[r$manifest$batch_col]][ix]),source_line=line)
 m4_csv(cells,file.path(m4_out,paste0(dataset,'-source-cells.csv')))
 g<-expand.grid(n_per_condition=c(3L,4L,6L,8L,12L),m_cells_per_sample=c(1L,2L,4L,8L,16L,32L,64L),beta=log2(c(1.25,1.5,2)),KEEP.OUT.ATTRS=FALSE)
 g$dataset<-dataset;g$source_line<-line;g$available_cells<-nc;g$cells<-2*g$n_per_condition*g$m_cells_per_sample
 g$eligible<-g$cells<=nc;g$counted<-g$eligible & g$m_cells_per_sample>=3 & dataset!='derks2022_timsTOF'
 g$task<-seq_len(nrow(g))+task;task<-task+nrow(g);g$k<-1
 g$prediction_seed<-as.integer(seed+4000000+1000*g$task)
 g$assignment_seed<-as.integer(seed+5000000+1000*g$task)
 g$injection_seed<-as.integer(seed+6000000+1000*g$task)
 g$replicates<-200L
 g$status<-ifelse(!g$eligible,'INELIGIBLE: exceeds source-line cell count',ifelse(dataset=='derks2022_timsTOF','NOT EVALUABLE: no eligible m >= 3 point',ifelse(!g$counted,'EXCLUDED FROM GATE: m < 3, deterministic zero power','REGISTERED FOR GATE')))
 for(i in which(g$eligible)){
  a<-m4_assignment(cells,g$n_per_condition[i],g$m_cells_per_sample[i],g$assignment_seed[i])
  layout<-m4_layout(a,g$n_per_condition[i]);id<-sprintf('%04d',g$task[i])
  m4_csv(a,file.path(m4_out,'assignments',paste0(id,'-first-assignment.csv')))
  long<-data.frame(sample=rep(rownames(layout),ncol(layout)),batch=rep(colnames(layout),each=nrow(layout)),cells=as.vector(layout))
  m4_csv(long,file.path(m4_out,'layouts',paste0(id,'-layout.csv')))
 }
 grids[[dataset]]<-g
 sources[[dataset]]<-data.frame(dataset=dataset,source_line=line,cells=nc,calibration_proteins=length(cal),validation_proteins=length(ids)-length(cal),split_seed=split_seed,eligible_points=sum(g$eligible),counted_points=sum(g$counted),status=if(any(g$counted))'EVALUABLE' else 'NOT EVALUABLE')
 hashes[[dataset]]<-list(released_preset=m4_hash(paste0('inst/presets/',dataset,'.yaml')),quantitative_input=m4_hash(paste0('.local/datasets/',dataset,'-aggregated.rds')),split=m4_hash(file.path(m4_out,paste0(dataset,'-protein-split.csv'))),source_cells=m4_hash(file.path(m4_out,paste0(dataset,'-source-cells.csv'))))
}
grid<-do.call(rbind,grids);m4_csv(grid,file.path(m4_out,'registered-grid.csv'));m4_csv(do.call(rbind,sources),file.path(m4_out,'source-summary.csv'))
code<-c('R/scp_calibrate.R','R/scp_simulate.R','R/scp_power.R','inst/scripts/m4b_helpers.R','inst/scripts/prepare_m4b.R','inst/scripts/run_m4b_point.R','inst/scripts/run_m4b_cluster.R')
m4_json(list(mandate_addendum='1.1.5',registered_UTC=format(Sys.time(),tz='UTC',usetz=TRUE),seed=seed,versions=as.list(m4_versions),
 criterion=list(replicates=200L,absolute_power_gap=.10,dataset_pass_fraction=.80,overall_minimum_passing_datasets=2L,count_only_m_ge_3=TRUE),
 prediction_order='All counted-point prediction receipts frozen and hashed before the first observed endpoint is evaluated',
 calibration_population='Calibration protein half across native accepted SCP calibration cells; H1 source restriction applied to validation cells',
 truth='Validation protein half; native quantitative mask and batches retained; cell shifts from calibration proteins; log2 condition effect in 10% of validation proteins; external Walzik per-protein donor injection with half-tier fallback',
 k=1,amendment_I=TRUE,diagnostic_J_sha256=m4_hash('reports/J/diagnostic-report.md'),M2b_sha256=m4_hash('reports/M2b/fidelity-summary.csv'),non_TMT_results_not_rerun=TRUE,validation_tiers='Frozen M1 labels for selecting training-half donor fallback only',
 input_hashes=list(datasets=hashes,bulk_calibration=m4_hash('reports/M1-final/calibration-results.rds'),grid=m4_hash(file.path(m4_out,'registered-grid.csv')),mandate=m4_hash('scpPower_mandate_v1.1.md'),renv_lock=m4_hash('renv.lock'),code=as.list(setNames(vapply(code,m4_hash,character(1)),code))),
 wall_seconds=proc.time()[['elapsed']]-start),file.path(m4_out,'registration.json'))
cat('Registered',sum(grid$counted),'counted points and',sum(grid$eligible & !grid$counted),'constructional/non-evaluable points\n')
