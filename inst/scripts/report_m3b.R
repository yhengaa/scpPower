# Collect only complete, checksum-verified registered points; then render M3.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
source('R/scp_calibrate.R');started<-proc.time()[['elapsed']];out<-'reports/M3b'
grid<-read.csv('reports/M3/registered-grid.csv')
reg<-jsonlite::read_json(file.path(out,'grid-registration.json'),simplifyVector=FALSE)
old_results<-read.csv('reports/M3b/archive-v0.2.0/power-grid.csv')
old_replicates<-read.csv('reports/M3/replicate-metrics.csv.gz')
old_manifests<-jsonlite::read_json('reports/M3/grid-point-manifests.json.gz',simplifyVector=FALSE)
results<-list();manifests<-list();replicate_con<-gzfile(file.path(out,'replicate-metrics.csv.gz'),'wt')
for(i in seq_len(nrow(grid))) {
 if(i<=630L) {
 prefix<-file.path(out,'points',sprintf('%04d',i));mpath<-paste0(prefix,'-manifest.json')
 if(!file.exists(mpath))stop('Incomplete M3b grid: missing task ',i)
 m<-jsonlite::read_json(mpath,simplifyVector=FALSE)
 rp<-paste0(prefix,'-results.csv');qp<-paste0(prefix,'-replicates.csv.gz')
 stopifnot(digest::digest(file=rp,algo='sha256')==m$output_file_hashes$results,
  digest::digest(file=qp,algo='sha256')==m$output_file_hashes$replicates)
 r<-read.csv(rp);q<-read.csv(qp)
 } else {
 r<-old_results[old_results$task==i,,drop=FALSE]
 q<-old_replicates[old_replicates$task==i,,drop=FALSE]
 m<-old_manifests[[as.character(i)]]
 }
 stopifnot(nrow(r)==1L,nrow(q)==200L,identical(q$seed,as.integer(grid$seed[i]+seq_len(200L))),
  r$task==i,r$dataset==grid$dataset[i],r$n_per_condition==grid$n_per_condition[i],
  r$m_cells_per_sample==grid$m_cells_per_sample[i],r$beta==grid$beta[i],r$k==grid$k[i])
 results[[i]]<-r;manifests[[as.character(i)]]<-m
 for(field in setdiff(c('conditional_power','testable_changed'),names(q)))q[[field]]<-NA_real_
 if(i==1L)qfields<-names(q)
 for(field in setdiff(qfields,names(q)))q[[field]]<-NA_real_
 q<-q[,qfields,drop=FALSE]
 write.table(q,replicate_con,sep=',',row.names=FALSE,col.names=i==1L,quote=TRUE,na='NA')
}
close(replicate_con)
fields<-unique(unlist(lapply(results,names)));results<-lapply(results,function(x){for(f in setdiff(fields,names(x)))x[[f]]<-NA_real_;x[,fields,drop=FALSE]});d<-do.call(rbind,results);write.csv(d,file.path(out,'power-grid.csv'),row.names=FALSE,na='NA')
write.csv(d[,c('dataset','n_per_condition','m_cells_per_sample','beta','k','realised_fdr','realised_fdr_se','replicates','software_versions')],file.path(out,'realised-fdr.csv'),row.names=FALSE)
write.csv(d[,c('dataset','n_per_condition','m_cells_per_sample','beta','k','cells','runs','cells_per_run','minutes_per_run','instrument_hours','wall_seconds','software_versions')],file.path(out,'instrument-hours.csv'),row.names=FALSE)
cheapest<-list()
for(n in unique(d$dataset))for(b in unique(d$beta))for(kk in unique(d$k[d$dataset==n])) {
 x<-d[d$dataset==n&d$beta==b&d$k==kk&is.finite(d$power)&d$power>=.8,,drop=FALSE]
 if(nrow(x)) {
  x<-x[order(x$instrument_hours,x$cells,x$n_per_condition),,drop=FALSE][1,,drop=FALSE]
  row<-data.frame(dataset=n,beta=b,k=kk,target_power=.8,n_per_condition=x$n_per_condition,m_cells_per_sample=x$m_cells_per_sample,power=x$power,simulation_se=x$simulation_se,runs=x$runs,instrument_hours=x$instrument_hours,status='Target reached in registered grid',software_versions=x$software_versions)
 }else row<-data.frame(dataset=n,beta=b,k=kk,target_power=.8,n_per_condition=NA,m_cells_per_sample=NA,power=NA,simulation_se=NA,runs=NA,instrument_hours=NA,status='No registered design reaches target',software_versions=d$software_versions[1])
 cheapest[[length(cheapest)+1L]]<-row
}
cheapest<-do.call(rbind,cheapest);write.csv(cheapest,file.path(out,'cheapest-designs.csv'),row.names=FALSE,na='NA')
vt<-paste(names(.scp_versions()),.scp_versions(),sep='=',collapse='; ')
colors<-c('#2374ab','#188977','#a98600','#bf5a23','#8d4e9b');ns<-sort(unique(d$n_per_condition));betas<-sort(unique(d$beta))
render<-function(metric,path,sensitivity=FALSE) {
 groups<-if(sensitivity)expand.grid(dataset=c('specht2019v3','leduc2022_pSCoPE'),k=c(1,.75,.5),stringsAsFactors=FALSE) else data.frame(dataset=unique(d$dataset),k=1)
 png(file.path(out,path),width=3300,height=as.integer(850*nrow(groups)+180),res=180)
 par(mfrow=c(nrow(groups),length(betas)),mar=c(4.2,4.8,3.0,1.0),oma=c(4,0,1,0))
 for(i in seq_len(nrow(groups)))for(b in betas) {
  g<-groups[i,];x<-d[d$dataset==g$dataset&d$k==g$k&d$beta==b,]
  plot(NA,xlim=c(0,6),ylim=c(0,1),xlab='Cells per sample (m)',ylab=if(metric=='power')'Power' else 'Testability',xaxt='n',main=paste(g$dataset,'k=',g$k,'fold=',format(2^b,digits=3)))
  axis(1,at=0:6,labels=2^(0:6));if(metric=='power')abline(h=.8,lty=3,col='grey60')
  for(j in seq_along(ns)){z<-x[x$n_per_condition==ns[j],];z<-z[order(z$m_cells_per_sample),];lines(log2(z$m_cells_per_sample),z[[metric]],type='o',pch=16,col=colors[j],cex=.4,lwd=1.2)}
 }
 par(xpd=NA);legend(x=grconvertX(.5,from='ndc'),y=grconvertY(70/(850*nrow(groups)+180),from='ndc'),xjust=.5,yjust=.5,legend=paste('n=',ns,'per condition'),col=colors,lty=1,pch=16,horiz=TRUE,bty='n',cex=.75)
 mtext('Registered grid; balanced layouts; 200 replicates per point; k values are assumptions; accepted validation range remains limited.',side=1,outer=TRUE,line=2,cex=.7)
 mtext(vt,side=1,outer=TRUE,line=3.1,cex=.5);dev.off()
}
render('power','power-curves.png');render('testability','testability-curves.png')
render('power','power-tmt-sensitivity.png',TRUE);render('testability','testability-tmt-sensitivity.png',TRUE)
f<-function(x,digits=3)ifelse(is.finite(x),format(round(x,digits),nsmall=digits,trim=TRUE),'NA')
lines<-c('# M3b updated carrier-TMT simulation power report','13 October 2026. Registered M3 grid v1.1.3, with amendment-I masks under v1.1.5. All five M2 series passed before the original M3. Only the two TMT presets use updated endpoint masks; variances, parameter draws, layouts, k assumptions and registered seeds remain fixed. This is simulation-based planning; the accepted M4 and M4b verdicts are unchanged.','',sprintf('The combined table contains every registered grid point: %d points with %d independent replicates each. TMT uses k=1, 0.75 and 0.5; DIA and label-free use k=1.',nrow(d),reg$replicates_per_point),sprintf('Maximum power simulation standard error: %s percentage points. Grid-point wall times range from %s to %s seconds (sum %s worker wall hours; excludes scheduler waiting and source preparation).',f(100*max(d$simulation_se),2),f(min(d$wall_seconds),1),f(max(d$wall_seconds),1),f(sum(d$wall_seconds)/3600,2)),'','## Endpoint and interpretation','', 'Sample means require at least three quantified cells in every sample. Limma-trend fits design ~ condition; BH uses 0.05 among testable proteins. Changed proteins require sign agreement for detection. Power has all truly changed proteins in the denominator, including untestable proteins. Testability counts all proteins satisfying the cell minimum. FDR is mean replicate FDP among BH rejections, null rejections divided by all rejections, and zero with no rejections; wrong-sign changed rejections remain separately recorded. Standard errors are replicate SD divided by sqrt(replicates).','', 'Designs with fewer than three cells per sample cannot meet the endpoint and have zero testability and power. Every registered point is still simulated; no points are removed. Increasing cells can reduce cell noise and improve testability, but biological-sample variance remains. In TMT, identification/detection does not respond to effect size because the pooled carrier governs p_g and q_g.','', '## Realised FDR','', '| Preset | k | Maximum realised FDR over grid |','| --- | ---: | ---: |')
for(n in unique(d$dataset))for(kk in unique(d$k[d$dataset==n])) {
 x<-d[d$dataset==n&d$k==kk,];lines<-c(lines,sprintf('| %s | %s | %s |',n,kk,f(max(x$realised_fdr),4)))
}
if(any(d$realised_fdr>.05))lines<-c(lines,'','Some simulated designs exceed nominal BH 0.05 in realised FDR. The full table shows the extent by design; no variance, detection, test, or threshold was adjusted. Monte Carlo standard errors are reported alongside estimates; these are estimates under the registered model, not guarantees of experimental FDR.')
lines<-c(lines,'','## Cheapest tested design reaching target power','', '| Preset | Fold change | k | n per condition | m cells per sample | Power | Instrument hours | Status |','| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |')
for(i in seq_len(nrow(cheapest))) {x<-cheapest[i,];lines<-c(lines,sprintf('| %s | %s | %s | %s | %s | %s | %s | %s |',x$dataset,f(2^x$beta,2),x$k,x$n_per_condition,x$m_cells_per_sample,f(x$power),f(x$instrument_hours,2),x$status))}
lines<-c(lines,'','The target is mean simulated power 0.8 within this finite grid; it is not a lower confidence bound. An unavailable target is reported explicitly rather than extrapolated.','', '## Instrument time','', 'Runs count actual occupied batches in the balanced simulation layout. A three-channel plexDIA run contains two cells under the registered equal-A/B rule, leaving one slot empty; hours therefore use actual runs rather than dividing cells by nominal channel capacity. Preset durations are cited in YAML and parsed from the same registered source papers. Specht and Leduc use total published TMT method durations; Derks uses the full single-cell instrument-specific gradient; Brunner uses its published single-cell EvoSep gradient. Unreported handling and turnover overhead are excluded. These are nominal hours, not measured lab scheduling times.','', 'Whole-PBMC donor variance transfer remains an assumption, including donor composition differences. k values are registered compression assumptions. Passing M2 does not validate power or donor transfer, and the registered M4 criterion remains unchanged.','', 'Outputs: [full power grid](power-grid.csv), [FDR table](realised-fdr.csv), [instrument-hour table](instrument-hours.csv), [cheapest designs](cheapest-designs.csv), [compressed replicate metrics](replicate-metrics.csv.gz), [manifest](run-manifest.json).','',paste('Software:',vt))
writeLines(lines,file.path(out,'power-report.md'))
writeLines(c('# M3 evidence','', '[Power report](power-report.md), [power curves](power-curves.png), [testability curves](testability-curves.png), [TMT power sensitivity](power-tmt-sensitivity.png), [TMT testability sensitivity](testability-tmt-sensitivity.png), [FDR](realised-fdr.csv), [instrument hours](instrument-hours.csv), [cheapest designs](cheapest-designs.csv), [manifest](run-manifest.json).'),file.path(out,'README.md'))
# One compressed manifest retains every task's seeds, inputs, outputs and cluster host.
con<-gzfile(file.path(out,'grid-point-manifests.json.gz'),'wt')
writeLines(jsonlite::toJSON(manifests,pretty=FALSE,auto_unbox=TRUE,na='null',digits=NA),con);close(con)
m<-list(mandate_addendum='1.1.5',seed=20261013,versions=as.list(.scp_versions()),package_version='1.0.0',grid_points=nrow(d),replicates_per_point=200,total_replicates=nrow(d)*200,grid_registration=reg,input_hashes=list(original_non_TMT_grid=digest::digest(file='reports/M3b/archive-v0.2.0/power-grid.csv',algo='sha256'),grid=digest::digest(file=file.path(out,'registered-grid.csv'),algo='sha256'),instrument_metadata=digest::digest(file=file.path(out,'instrument-preset-manifest.json'),algo='sha256')),per_point_manifests='grid-point-manifests.json.gz',per_replicate_metrics='replicate-metrics.csv.gz',worker_wall_seconds=sum(d$wall_seconds[d$task<=630L]),reused_non_TMT_points=315L,rerun_TMT_points=630L,rerun_replicates=126000L,original_non_TMT_manifest=jsonlite::read_json('reports/M3b/archive-v0.2.0/run-manifest.json',simplifyVector=FALSE),report_wall_seconds=proc.time()[['elapsed']]-started,endpoint=manifests[[1]]$endpoint,instrument_time=manifests[[1]]$instrument_time_definition,decision=list(M2_all_five_passed=TRUE,M3_registered_grid_complete=TRUE,accepted_M4_unchanged=TRUE,M4_criterion_unchanged=TRUE,no_tuning=TRUE),wall_seconds=sum(d$wall_seconds[d$task<=630L])+proc.time()[['elapsed']]-started)
jsonlite::write_json(m,file.path(out,'run-manifest.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA,na='null')
cat('All registered M3 points verified and reported.\n')
