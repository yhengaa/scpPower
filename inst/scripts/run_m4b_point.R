# One registered M4 point. All prediction receipts precede any observation.
source('inst/scripts/m4b_helpers.R')
m4_point <- function(task,stage){
 start<-proc.time()[['elapsed']];grid<-read.csv(file.path(m4_out,'registered-grid.csv'))
 g<-grid[grid$task==task,,drop=FALSE];stopifnot(nrow(g)==1L,g$counted)
 dataset<-g$dataset;prefix<-file.path(m4_out,'points',sprintf('%04d',task))
 reg<-jsonlite::read_json(file.path(m4_out,'registration.json'),simplifyVector=FALSE)
 for(path in names(reg$input_hashes$code))stopifnot(m4_hash(path)==reg$input_hashes$code[[path]])
 for(d in names(reg$input_hashes$datasets))stopifnot(m4_hash(paste0('inst/presets/',d,'.yaml'))==reg$input_hashes$datasets[[d]]$released_preset)
 assignment_path<-file.path(m4_out,'assignments',paste0(sprintf('%04d',task),'-first-assignment.csv'))
 first<-read.csv(assignment_path);layout<-m4_layout(first,g$n_per_condition)
 preset_path<-paste0('inst/presets/m4b/',dataset,'.yaml')
 if(stage=='prediction'){
  result<-scp_power(preset_path,g$n_per_condition,g$m_cells_per_sample,g$beta,
   layout=layout,k=1,replicates=200L,seed=g$prediction_seed)
  result$results$task<-task;result$results$dataset<-dataset
  m4_csv(result$results,paste0(prefix,'-predicted.csv'))
  r<-result$replicate_metrics;r$task<-task;r$dataset<-dataset
  con<-gzfile(paste0(prefix,'-predicted-replicates.csv.gz'),'wt');write.csv(r,con,row.names=FALSE,na='NA');close(con)
  receipt<-result$manifest
  receipt$stage<-'prediction';receipt$completed_UTC<-format(Sys.time(),tz='UTC',usetz=TRUE)
  receipt$input_hashes$preset_file<-m4_hash(preset_path);receipt$input_hashes$first_assignment<-m4_hash(assignment_path)
  receipt$output_file_hashes<-list(results=m4_hash(paste0(prefix,'-predicted.csv')),replicates=m4_hash(paste0(prefix,'-predicted-replicates.csv.gz')))
 }else{
  stopifnot(stage=='observation',file.exists(file.path(m4_out,'predictions-frozen.json')))
  frozen<-jsonlite::read_json(file.path(m4_out,'predictions-frozen.json'),simplifyVector=FALSE)
  expected<-frozen$prediction_hashes[[as.character(task)]]
  stopifnot(!is.null(expected),m4_hash(paste0(prefix,'-predicted.csv'))==expected$results,
   m4_hash(paste0(prefix,'-prediction-manifest.json'))==expected$manifest,
   m4_hash(preset_path)==expected$preset)
  obs<-readRDS(paste0('.local/datasets/',dataset,'-aggregated.rds'))
  val_path<-file.path(m4_out,paste0(dataset,'-validation-donor-variances.csv'));val<-read.csv(val_path)
  shifts<-read.csv(file.path(m4_out,paste0(dataset,'-calibration-cell-shifts.csv')))
  cells<-read.csv(file.path(m4_out,paste0(dataset,'-source-cells.csv')))
  base<-obs$matrix[val$protein,cells$cell,drop=FALSE]
  base<-sweep(base,2,shifts$log2_shift[match(colnames(base),shifts$cell)],'-')
  stopifnot(!anyNA(shifts$log2_shift[match(colnames(base),shifts$cell)]))
  records<-vector('list',200L);np<-nrow(base);ns<-2L*g$n_per_condition
  for(j in seq_len(200L)){
   t0<-proc.time()[['elapsed']]
   assignment_seed<-as.integer(g$assignment_seed+j-1L);injection_seed<-as.integer(g$injection_seed+j)
   md<-m4_assignment(cells,g$n_per_condition,g$m_cells_per_sample,assignment_seed)
   if(j==1L)stopifnot(identical(md[,c('cell','sample','condition','batch')],first[,c('cell','sample','condition','batch')]))
   set.seed(injection_seed);changed<-rep(FALSE,np);changed[sample.int(np,round(.1*np),replace=FALSE)]<-TRUE
   donor<-matrix(rnorm(np*ns),np)*sqrt(val$sigma2_s)
   y<-base[,md$cell,drop=FALSE]+donor[,match(md$sample,paste0('S',seq_len(ns))),drop=FALSE]+outer(changed*g$beta,as.integer(md$condition=='B'))
   # Native NAs remain NAs: semi-synthetic injection cannot create measurements.
   stopifnot(identical(is.na(y),is.na(base[,md$cell,drop=FALSE])))
   truth<-data.frame(changed=changed,beta_observed=changed*g$beta)
   e<-.scp_primary_endpoint(list(matrix=y,cell_metadata=md,protein_metadata=truth))
   records[[j]]<-data.frame(task=task,dataset=dataset,replicate=j,
    assignment_seed=assignment_seed,injection_seed=injection_seed,
    n_per_condition=g$n_per_condition,m_cells_per_sample=g$m_cells_per_sample,beta=g$beta,k=1,
    power=e$power,conditional_power=e$conditional_power,testable_changed=e$testable_changed,
    testability=e$testability,realised_fdr=e$realised_fdr,changed=e$changed,detected=e$detected,
    testable=e$testable,rejected=e$rejected,false_discoveries=e$false_discoveries,wrong_sign_changed=e$wrong_sign_changed,
    runs=length(unique(md$batch)),matrix_sha256=digest::digest(y,algo='sha256'),
    assignment_sha256=digest::digest(md,algo='sha256'),changed_proteins_sha256=digest::digest(val$protein[changed],algo='sha256'),
    donor_effects_sha256=digest::digest(donor,algo='sha256'),wall_seconds=proc.time()[['elapsed']]-t0,software_versions=m4_vtext)
  }
  r<-do.call(rbind,records);avg<-function(x)if(all(is.na(x)))NA_real_ else mean(x,na.rm=TRUE)
  se<-function(x)if(sum(is.finite(x))>1)sd(x,na.rm=TRUE)/sqrt(sum(is.finite(x))) else NA_real_
  summary<-data.frame(task=task,dataset=dataset,n_per_condition=g$n_per_condition,m_cells_per_sample=g$m_cells_per_sample,beta=g$beta,k=1,
   power=mean(r$power),conditional_power=avg(r$conditional_power),conditional_power_se=se(r$conditional_power),
   testability=mean(r$testability),realised_fdr=mean(r$realised_fdr),simulation_se=se(r$power),
   testability_se=se(r$testability),realised_fdr_se=se(r$realised_fdr),
   mean_rejections=mean(r$rejected),zero_rejection_fraction=mean(r$rejected==0),
   mean_runs=mean(r$runs),replicates=200L,wall_seconds=proc.time()[['elapsed']]-start)
  m4_csv(summary,paste0(prefix,'-observed.csv'))
  con<-gzfile(paste0(prefix,'-observed-replicates.csv.gz'),'wt');write.csv(r,con,row.names=FALSE,na='NA');close(con)
  receipt<-list(stage='observation',seed=g$injection_seed,assignment_seed=g$assignment_seed,versions=as.list(m4_versions),
   completed_UTC=format(Sys.time(),tz='UTC',usetz=TRUE),replicates=200L,
   input_hashes=list(predictions_frozen=m4_hash(file.path(m4_out,'predictions-frozen.json')),prediction=expected,
    original_quantitative_input=m4_hash(paste0('.local/datasets/',dataset,'-aggregated.rds')),
    donor_variances=m4_hash(val_path),source_cells=m4_hash(file.path(m4_out,paste0(dataset,'-source-cells.csv'))),
    shifts=m4_hash(file.path(m4_out,paste0(dataset,'-calibration-cell-shifts.csv'))),source_matrix=digest::digest(base,algo='sha256')),
   output_file_hashes=list(results=m4_hash(paste0(prefix,'-observed.csv')),replicates=m4_hash(paste0(prefix,'-observed-replicates.csv.gz'))),
   native_missingness_preserved=TRUE,between_donor_variance_injected_not_validated=TRUE,k_not_validated=TRUE)
 }
 receipt$task<-task;receipt$dataset<-dataset
 receipt$cluster<-list(job=Sys.getenv('SLURM_JOB_ID'),host=Sys.info()[['nodename']])
 receipt$wall_seconds<-proc.time()[['elapsed']]-start
 m4_json(receipt,paste0(prefix,'-',stage,'-manifest.json'))
 cat(stage,'completed',task,dataset,'n',g$n_per_condition,'m',g$m_cells_per_sample,'beta',g$beta,'\n');flush.console()
 task
}
