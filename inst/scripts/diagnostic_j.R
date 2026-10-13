# Diagnostic J precedes every amendment-I change. Uses v0.3.0 validation halves.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
source('R/scp_calibrate.R');source('R/scp_simulate.R')
start<-proc.time()[['elapsed']];seed<-20261013L;out<-'reports/J';dir.create(out,recursive=TRUE,showWarnings=FALSE)
vs<-.scp_versions();vt<-paste(names(vs),vs,sep='=',collapse='; ')
csv<-function(d,p){d$software_versions<-rep(vt,nrow(d));write.csv(d,file.path(out,p),row.names=FALSE,na='NA')}
pq<-function(mask,batch){runs<-split(seq_len(ncol(mask)),batch);identified<-vapply(runs,function(j)rowSums(mask[,j,drop=FALSE])>0,logical(nrow(mask)));den<-as.vector(identified%*%as.integer(lengths(runs)));data.frame(p_g=rowMeans(identified),q_g=ifelse(den>0,rowSums(mask)/den,0),r_g=rowMeans(mask))}
all_summary<-list();all_tiers<-list();all_depth<-list();all_reps<-list();records<-list()
for(di in seq_along(c('specht2019v3','leduc2022_pSCoPE'))){
 d<-c('specht2019v3','leduc2022_pSCoPE')[di];t0<-proc.time()[['elapsed']]
 psm<-readRDS(paste0('.local/datasets/',d,'-psm-proteins.rds'));agg<-readRDS(paste0('.local/datasets/',d,'-aggregated.rds'))
 split<-read.csv(paste0('reports/M4/',d,'-protein-split.csv'));split<-split[split$half=='validation',]
 source<-read.csv(paste0('reports/M4/',d,'-source-cells.csv'))
 ids<-split$protein;cells<-source$cell;batch<-source$batch
 a<-is.finite(psm$matrix[ids,cells,drop=FALSE]);b<-is.finite(agg$matrix[ids,cells,drop=FALSE]);pa<-pq(a,batch);pb<-pq(b,batch)
 differences<-data.frame(protein=ids,tier=split$released_tier,pa,aggregated_p_g=pb$p_g,aggregated_q_g=pb$q_g,aggregated_r_g=pb$r_g,p_difference=pb$p_g-pa$p_g,q_difference=pb$q_g-pa$q_g,psm_finite_aggregated_NA_fraction=rowMeans(a & !b),aggregated_finite_psm_NA_fraction=rowMeans(b & !a))
 csv(differences,paste0(d,'-protein-mask-differences.csv'))
 table<-do.call(rbind,lapply(c(0L,sort(unique(split$released_tier))),function(t){j<-if(t==0L)seq_along(ids) else which(split$released_tier==t);data.frame(dataset=d,tier=t,proteins=length(j),cells=length(cells),psm_finite_fraction=mean(a[j,,drop=FALSE]),aggregated_finite_fraction=mean(b[j,,drop=FALSE]),aggregated_sparser_gap=mean(a[j,,drop=FALSE])-mean(b[j,,drop=FALSE]),psm_finite_aggregated_NA_fraction=mean(a[j,,drop=FALSE]&!b[j,,drop=FALSE]),aggregated_finite_psm_NA_fraction=mean(b[j,,drop=FALSE]&!a[j,,drop=FALSE]),median_p_difference=median(pb$p_g[j]-pa$p_g[j]),median_q_difference=median(pb$q_g[j]-pa$q_g[j]))}))
 all_tiers[[d]]<-table;summary<-table[table$tier==0L,];summary$amendment_I_applies<-summary$aggregated_sparser_gap>=.05;all_summary[[d]]<-summary
 # Observe both masks on exactly the same cells and validation proteins.
 depth<-data.frame(dataset=d,cell=cells,batch=batch,psm_proteins=colSums(a),aggregated_proteins=colSums(b))
 csv(depth,paste0(d,'-cell-depth.csv'))
 # Detection-only G2 at the exact native source-line layout. Fix protein
 # identities and margins to distinguish independent p/q draws from native
 # cell/channel dependence; this diagnostic never writes a released preset.
 runs<-unique(batch);reps<-list()
 for(maskname in c('PSM','aggregated')){
  rates<-if(maskname=='PSM')pa else pb
  for(j in seq_len(20L)){
   simseed<-as.integer(seed+7000000L+10000L*di+1000L*(maskname=='aggregated')+j)
   RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(simseed)
   identified<-matrix(runif(length(ids)*length(runs)),length(ids))<rates$p_g
   quantified<-matrix(runif(length(ids)*length(cells)),length(ids))<rates$q_g
   sim<-identified[,match(batch,runs),drop=FALSE]&quantified
   reps[[length(reps)+1L]]<-data.frame(dataset=d,mask=maskname,replicate=j,seed=simseed,proteins=length(ids),cells=length(cells),depth_sd=sd(colSums(sim)),mean_depth=mean(colSums(sim)),matrix_sha256=digest::digest(sim,algo='sha256'))
  }
 }
 rep<-do.call(rbind,reps);all_reps[[d]]<-rep
 for(maskname in c('PSM','aggregated')){
  rr<-rep[rep$mask==maskname,];observed<-depth[[if(maskname=='PSM')'psm_proteins' else 'aggregated_proteins']]
  all_depth[[paste(d,maskname)]]<-data.frame(dataset=d,mask=maskname,observed_depth_sd=sd(observed),simulated_mean_depth_sd=mean(rr$depth_sd),simulated_replicate_sd=sd(rr$depth_sd),observed_mean_depth=mean(observed),simulated_mean_depth=mean(rr$mean_depth),replicates=20L,simulation='Native batch occupancy; fixed validation protein p/q; independent G2 Bernoulli draws, no protein resampling')
 }
 records[[d]]<-list(seed=seed,split_seed=unique(split$split_seed),versions=as.list(vs),source_line=unique(source$source_line),
  input_hashes=list(PSM_file=digest::digest(file=paste0('.local/datasets/',d,'-psm-proteins.rds'),algo='sha256'),aggregated_file=digest::digest(file=paste0('.local/datasets/',d,'-aggregated.rds'),algo='sha256'),split=digest::digest(file=paste0('reports/M4/',d,'-protein-split.csv'),algo='sha256'),source_cells=digest::digest(file=paste0('reports/M4/',d,'-source-cells.csv'),algo='sha256'),psm_mask=digest::digest(a,algo='sha256'),aggregated_mask=digest::digest(b,algo='sha256')),
  decision=if(summary$amendment_I_applies)'I APPLIES: overall aggregated sparsity gap >= 5 pp' else 'I DOES NOT APPLY: gap < 5 pp; document dependence under K',wall_seconds=proc.time()[['elapsed']]-t0)
}
summary<-do.call(rbind,all_summary);csv(summary,'decision-table.csv');csv(do.call(rbind,all_tiers),'mask-summary-by-tier.csv');csv(do.call(rbind,all_depth),'native-depth-summary.csv');csv(do.call(rbind,all_reps),'depth-replicate-metrics.csv')
jsonlite::write_json(list(seed=seed,versions=as.list(vs),decision_threshold=.05,decision_basis='Exact same v0.3.0 held-out proteins and H1 source cells; PSM finite fraction minus endpoint aggregated finite fraction',datasets=records,amendment_I_applies=all(summary$amendment_I_applies),preset_changes_before_diagnostic=FALSE,input_hashes=list(script=digest::digest(file='inst/scripts/diagnostic_j.R',algo='sha256'),simulator=digest::digest(file='R/scp_simulate.R',algo='sha256'),renv_lock=digest::digest(file='renv.lock',algo='sha256')),wall_seconds=proc.time()[['elapsed']]-start),file.path(out,'run-manifest.json'),pretty=TRUE,auto_unbox=TRUE,na='null',digits=NA)
print(summary[,c('dataset','psm_finite_fraction','aggregated_finite_fraction','aggregated_sparser_gap','amendment_I_applies')])
