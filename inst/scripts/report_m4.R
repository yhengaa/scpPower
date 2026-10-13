# Report all registered results without changing simulation, calibration or gate.
source('inst/scripts/m4_helpers.R')
start<-proc.time()[['elapsed']];grid<-read.csv(file.path(m4_out,'registered-grid.csv'))
sources<-read.csv(file.path(m4_out,'source-summary.csv'));records<-list();verdicts<-list();manifests<-list();replicates<-list()
for(dataset in sources$dataset){
 d<-grid[grid$dataset==dataset & grid$eligible,,drop=FALSE];rows<-vector('list',nrow(d))
 for(i in seq_len(nrow(d))){g<-d[i,];prefix<-file.path(m4_out,'points',sprintf('%04d',g$task))
  row<-data.frame(g,predicted_power=NA_real_,observed_power=NA_real_,absolute_gap=NA_real_,within_10_pp=FALSE,
   predicted_testability=NA_real_,observed_testability=NA_real_,observed_realised_fdr=NA_real_,predicted_realised_fdr=NA_real_,
   predicted_simulation_se=NA_real_,observed_simulation_se=NA_real_,gap_simulation_se=NA_real_,
   predicted_conditional_power=NA_real_,observed_conditional_power=NA_real_,
   predicted_mean_rejections=NA_real_,observed_mean_rejections=NA_real_,predicted_zero_rejection_fraction=NA_real_,observed_zero_rejection_fraction=NA_real_,
   prediction_wall_seconds=0,observation_wall_seconds=0,actual_replicates=0L)
  if(!g$counted){
   row$predicted_power<-row$observed_power<-row$absolute_gap<-0
   row$predicted_testability<-row$observed_testability<-row$observed_realised_fdr<-row$predicted_realised_fdr<-0
   row$predicted_simulation_se<-row$observed_simulation_se<-row$gap_simulation_se<-0
   row$predicted_mean_rejections<-row$observed_mean_rejections<-0
   row$predicted_zero_rejection_fraction<-row$observed_zero_rejection_fraction<-1
   row$within_10_pp<-TRUE
  }else{
   paths<-c(predicted=paste0(prefix,'-predicted.csv'),observed=paste0(prefix,'-observed.csv'))
   if(all(file.exists(paths))){
    p<-read.csv(paths['predicted']);o<-read.csv(paths['observed'])
    row$predicted_power<-p$power;row$observed_power<-o$power;row$absolute_gap<-abs(p$power-o$power)
    row$within_10_pp<-is.finite(row$absolute_gap) & row$absolute_gap<=.10
    row$predicted_testability<-p$testability;row$observed_testability<-o$testability
    row$predicted_realised_fdr<-p$realised_fdr;row$observed_realised_fdr<-o$realised_fdr
    row$predicted_simulation_se<-p$simulation_se;row$observed_simulation_se<-o$simulation_se
    row$gap_simulation_se<-sqrt(p$simulation_se^2+o$simulation_se^2)
    row$predicted_conditional_power<-p$conditional_power;row$observed_conditional_power<-o$conditional_power
    row$predicted_mean_rejections<-p$mean_rejections;row$observed_mean_rejections<-o$mean_rejections
    row$predicted_zero_rejection_fraction<-p$zero_rejection_fraction;row$observed_zero_rejection_fraction<-o$zero_rejection_fraction
    row$prediction_wall_seconds<-p$wall_seconds;row$observation_wall_seconds<-o$wall_seconds
    row$actual_replicates<-200L
    row$status<-if(row$within_10_pp)'WITHIN 10 PP' else 'OUTSIDE 10 PP'
    for(stage in c('prediction','observation')){
     mp<-paste0(prefix,'-',stage,'-manifest.json');m<-jsonlite::read_json(mp,simplifyVector=FALSE)
     m$receipt_file_sha256<-m4_hash(mp);manifests[[paste(g$task,stage,sep='-')]]<-m
     rp<-paste0(prefix,if(stage=='prediction')'-predicted-replicates.csv.gz' else '-observed-replicates.csv.gz')
     con<-gzfile(rp,'rt');rr<-read.csv(con);close(con);rr$stage<-stage
     replicates[[paste(g$task,stage,sep='-')]]<-rr
    }
   }else row$status<-'ERROR: missing stage receipt, retained without tuning'
  }
  rows[[i]]<-row
 }
 table<-do.call(rbind,rows);records[[dataset]]<-table
 m4_csv(table,file.path(m4_out,paste0(dataset,'-validation-table.csv')))
 counted<-table[table$counted,,drop=FALSE];total<-nrow(counted);passed<-sum(counted$within_10_pp)
 complete<-sum(counted$actual_replicates==200L)
 fraction<-if(total)passed/total else NA_real_
 status<-if(!total)'NOT EVALUABLE' else if(complete<total)'FAIL: incomplete evaluation' else if(fraction>=.8)'PASS' else 'FAIL'
 verdicts[[dataset]]<-data.frame(dataset=dataset,source_line=unique(table$source_line),source_cells=unique(table$available_cells),
  eligible_points=nrow(table),excluded_m_lt_3=sum(!table$counted),counted_points=total,completed_points=complete,
  within_10_pp_points=passed,pass_fraction=fraction,status=status,
  maximum_absolute_gap=if(total && any(is.finite(counted$absolute_gap)))max(counted$absolute_gap,na.rm=TRUE) else NA_real_,
  maximum_observed_se=if(total && any(is.finite(counted$observed_simulation_se)))max(counted$observed_simulation_se,na.rm=TRUE) else NA_real_)
 png(file.path(m4_out,paste0(dataset,'-predicted-observed.png')),width=1500,height=1300,res=180)
 par(mar=c(5,5,4,1));plot(0,0,type='n',xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i',
  xlab='Predicted unconditional power',ylab='Observed unconditional power',main=paste(dataset,'— registered M4'))
 polygon(c(0,1,1,0),c(.1,1.1,.9,-.1),col='#E6EDF2',border=NA);abline(a=0,b=1,col='#788995');abline(a=c(-.1,.1),b=1,lty=2,col='#788995');grid()
 colors<-c('#147D92','#D47B2C','#7B55AA');effects<-sort(unique(grid$beta))
 if(total){
  points(counted$predicted_power,counted$observed_power,pch=19,col=adjustcolor(colors[match(counted$beta,effects)],.75),cex=1.1)
  legend('topleft',legend=paste0(c('1.25','1.5','2'),' fold'),col=colors,pch=19,bty='n',cex=.9)
  mtext(sprintf('%s: %d/%d counted points within 10 pp (%.1f%%)',status,passed,total,100*fraction),side=3,line=.3,cex=.8)
 }else text(.5,.5,'NOT EVALUABLE\n10 cells per line; no eligible m >= 3 design',cex=.9)
 mtext(m4_vtext,side=1,line=3.8,cex=.49);dev.off()
}
all<-do.call(rbind,records);summary<-do.call(rbind,verdicts)
m4_csv(all,file.path(m4_out,'validation-table.csv'));m4_csv(summary,file.path(m4_out,'dataset-verdicts.csv'))
# Prediction and observation metrics differ in their recorded seed names.
cols<-Reduce(union,lapply(replicates,names))
replicates<-lapply(replicates,function(d){for(k in setdiff(cols,names(d)))d[[k]]<-NA;d[,cols]})
con<-gzfile(file.path(m4_out,'replicate-metrics.csv.gz'),'wt');write.csv(do.call(rbind,replicates),con,row.names=FALSE,na='NA');close(con)
con<-gzfile(file.path(m4_out,'point-manifests.json.gz'),'wt');writeLines(jsonlite::toJSON(manifests,pretty=TRUE,auto_unbox=TRUE,na='null',digits=NA),con);close(con)
np<-sum(summary$status=='PASS');overall<-if(np>=2L)'PASS' else 'FAIL'
lines<-c('# M4 semi-synthetic validation',
 '\nM4 validates cell-level noise, missingness and the primary test under the registered semi-synthetic construction. It cannot validate the transfer of between-donor variance from Walzik whole PBMC to these cell lines, because that variance is injected, or carrier compression k, which remains a user assumption. All predictions were completed and hashed before any observed endpoint was evaluated.',
 paste0('\n13 October 2026. Frozen addendum v1.1.4. Overall verdict: **',overall,'** (',np,' datasets pass; at least two required).'),
 '\n## Registered gate',
 'A dataset passes when the absolute predicted–observed unconditional power gap is at most 10 percentage points at at least 80% of eligible points with m >= 3. Every counted point uses 200 prediction and 200 observation replicates. Points with m=1 or 2 are listed as deterministic zero-power exclusions; they are not simulated and do not enter the pass fraction. H1 makes timsTOF not evaluable, not a pass. No model, calibration, tolerance or criterion was changed in response to the results.',
 '\n| Dataset | Source line | Cells | Counted points | Within 10 pp | Pass fraction | Verdict |',
 '| --- | --- | ---: | ---: | ---: | ---: | --- |')
for(i in seq_len(nrow(summary))){r<-summary[i,];lines<-c(lines,sprintf('| %s | %s | %d | %d | %d | %s | %s |',r$dataset,r$source_line,r$source_cells,r$counted_points,r$within_10_pp_points,if(is.finite(r$pass_fraction))sprintf('%.1f%%',100*r$pass_fraction) else 'NA',r$status))}
lines<-c(lines,'\n## Construction and separation',
 'Protein identities from the accepted G presets were split once per dataset, with the split seed and both halves recorded. Half-presets were refitted from the calibration proteins across the native accepted SCP calibration cells. H1 restricts validation source cells to the largest single line; it does not remove other native calibration cells from the training protein half. The released presets remain byte-identical. All tier summaries, missingness curves, p_g/q_g or r_g, variance estimates and cell median shifts in each new preset use calibration proteins only.',
 'The validation half retains native quantitative missingness and plex/run labels. Calibration-half cell median shifts are applied before injection; neither injected condition effects nor donor effects are normalised out. No missing value is filled and no cell-level noise or batch effect is newly drawn in the observed data: those components are already present in real cells. For each replicate, cells are randomly sampled without replacement and assigned to 2n pseudo-samples of m cells; pseudo-samples are assigned to A/B. Ten percent of validation proteins are randomly changed and receive beta in B. Independent normal donor effects use each protein’s external mapped Walzik variance, with the half-preset tier value for unmapped IDs. Frozen M1 tier labels select that fallback for validation proteins; no validation intensities are used to fit the half-preset.',
 'The first assignment at every eligible point is recorded as cell IDs and a pseudo-sample-by-native-batch count matrix. The prediction uses that exact custom matrix. Later observed replicates use independent source assignments; their seeds and assignment hashes are recorded. Brunner has one cell per batch. k=1 is used throughout M4. Predictions jointly resample training-half protein rows under G1; observations evaluate held-out real protein rows, rather than copying or fitting their cell-noise estimates.',
 '\n## Endpoint and uncertainty',
 'The primary endpoint remains per-sample means of log2 values, at least three quantified cells in every sample, limma-trend with design ~ condition, BH 0.05 and injected-sign agreement. Unconditional power divides detected changed proteins by all changed proteins. Conditional power divides by testable changed proteins, is undefined when none are testable, and is averaged over defined replicates. FDR remains mean null-protein FDP among BH rejections, zero when nothing is rejected. Mean rejections and zero-rejection fractions are reported beside it.',
 'Tables provide prediction and observation power SE separately, plus their quadrature sum for the gap. The gate compares mean powers with the fixed 10-pp tolerance, not confidence intervals. The plot band is ±10 pp. Finite Monte Carlo error and repeated reuse of the finite source cell pool limit interpretation; this is not an independent-donor experiment.',
 '\n## Disagreements and suspected causes')
for(dataset in summary$dataset){r<-summary[summary$dataset==dataset,];d<-all[all$dataset==dataset & all$counted,,drop=FALSE]
 lines<-c(lines,paste0('\n### ',dataset))
 if(!nrow(d)){lines<-c(lines,'Not evaluable: ten cells per source line cannot supply the smallest counted design (2 × 3 × 4 = 24 cells). No degenerate timsTOF simulations were run.');next}
 lines<-c(lines,sprintf('%s: %d/%d counted points within 10 pp. Maximum gap %.1f pp. Maximum observed power SE %.2f pp.',r$status,r$within_10_pp_points,r$counted_points,100*r$maximum_absolute_gap,100*r$maximum_observed_se))
 lines<-c(lines,sprintf('Observed mean rejections range from %.3f to %.3f per replicate; the no-rejection fraction ranges from %.1f%% to %.1f%%. Realised mean FDP ranges from %.4f to %.4f; %d counted points exceed nominal 0.05, without adjustment.',min(d$observed_mean_rejections),max(d$observed_mean_rejections),100*min(d$observed_zero_rejection_fraction),100*max(d$observed_zero_rejection_fraction),min(d$observed_realised_fdr),max(d$observed_realised_fdr),sum(d$observed_realised_fdr>.05)))
 if(dataset=='derks2022_QExactive')lines<-c(lines,'Q Exactive passes at very low power: every observed mean is below 2%, with at most 1.115 mean rejections per replicate. The absolute 10-pp criterion provides limited evidence in this region and does not establish agreement at useful high power. This limitation does not change its registered verdict.')
 bad<-d[!d$within_10_pp,,drop=FALSE]
 if(nrow(bad)){
  finite<-bad[is.finite(bad$absolute_gap),]
  if(nrow(finite)){
   under<-sum(finite$observed_power>finite$predicted_power)
   lines<-c(lines,sprintf('%d discrepant points have observed power above prediction; %d have prediction above observation. Mean testability in the discrepant points is %.3f predicted versus %.3f observed.',under,nrow(finite)-under,mean(finite$predicted_testability),mean(finite$observed_testability)))
  }
  if(dataset%in%c('specht2019v3','leduc2022_pSCoPE'))lines<-c(lines,'Suspected causes: real source-line plex and channel masks can have dependencies beyond independent p_g/q_g draws; training across the native lines can differ from the selected validation line. TMT raw-PSM detection and the endpoint’s ratio-based quantitative assay can also have different finite-value masks. These are interpretations of held-out discrepancies, not fitted corrections.')
  else lines<-c(lines,'Suspected causes: native real masks remain fixed when beta or donor effects are injected, whereas G3 shifts predicted detection with expected intensity. Protein-rate dependence, finite source-line selection, and departures of native cell/run errors from the independent normal model can contribute. No mask or variance was altered to match observations.')
 }else lines<-c(lines,'All counted points meet the registered 10-pp tolerance. This does not validate the donor-transfer or compression assumptions, or every feature of cell-level dependence.')
 lines<-c(lines,paste0('[Eligible-point table](',dataset,'-validation-table.csv) · [Comparison plot](',dataset,'-predicted-observed.png)'))
}
lines<-c(lines,'\n## Rejection counts and stop decision',
 'Read the FDR column with mean rejections and the zero-rejection fraction. Constructional zero-power points have zero rejections and zero FDP by definition; their zero-rejection fraction is one. Low FDR alone does not demonstrate useful detection.',
 if(overall=='PASS')'The overall registered criterion passes. The release is authorised after engineering checks; dataset failures and all discrepant points remain in the evidence.' else 'The overall registered criterion fails. Work stops after this report; no v0.3.0 tag or release is authorised. Suspected causes above are not used to retune the simulator, calibration or gate.',
 '\nNo M5 work was performed. [Full eligible-point table](validation-table.csv), [dataset verdicts](dataset-verdicts.csv), [replicate metrics](replicate-metrics.csv.gz), [manifest](run-manifest.json).',paste0('\nSoftware: ',m4_vtext))
writeLines(lines,file.path(m4_out,'validation-report.md'))
if(overall=='PASS')writeLines(c('# Registered M4 semi-synthetic validation',
 'M4 passes overall under the unchanged registered criterion, with H2 excluding constructional zero-power designs. Four datasets pass: Specht 67/72, Leduc 66/69, Derks Q Exactive 6/6 and Brunner 57/57 counted points within 10 percentage points. TimsTOF is not evaluable. Every counted point has 200 predictions and 200 held-out observations; all predictions were frozen before observed endpoints.',
 'The evidence retains all eight discrepant TMT points, FDR exceedances, seeds, hashes, versions, timings and cluster accounting. Q Exactive passes at observed power below 2%, so its agreement provides limited evidence about high-power designs. M4 tests cell-level noise, missingness and the primary test; it cannot validate whole-PBMC donor variance transfer or compression k. Released presets and M3 numerical results remain unchanged.',
 'Adds documented, tested conditional power, mean rejection counts and zero-rejection fractions beside the existing measures. Includes the one-page Reading the grid section. No M5 work.',
 paste0('Software: ',m4_vtext)),file.path(m4_out,'release-notes.md'))
m4_json(list(mandate_addendum='1.1.4',seed=20261013L,versions=as.list(m4_versions),
 registered_replicates=200L,counted_points=sum(grid$counted),prediction_replicates=sum(all$actual_replicates),observation_replicates=sum(all$actual_replicates),
 input_hashes=list(registration=m4_hash(file.path(m4_out,'registration.json')),predictions_frozen=m4_hash(file.path(m4_out,'predictions-frozen.json'))),
 source_summary=sources,criterion=list(absolute_gap=.10,dataset_fraction=.80,minimum_passing_datasets=2L,count_only_m_ge_3=TRUE),
 datasets=summary,overall_verdict=overall,passing_datasets=np,
 between_donor_transfer_not_validated=TRUE,k_not_validated=TRUE,no_retuning=TRUE,M3_numeric_results_unchanged=TRUE,M5_not_executed=TRUE,
 point_manifest_file='point-manifests.json.gz',point_manifest_sha256=m4_hash(file.path(m4_out,'point-manifests.json.gz')),
 replicate_metric_file='replicate-metrics.csv.gz',replicate_metric_sha256=m4_hash(file.path(m4_out,'replicate-metrics.csv.gz')),
 cluster_execution=jsonlite::read_json(file.path(m4_out,'cluster-execution.json'),simplifyVector=FALSE),
 report_wall_seconds=proc.time()[['elapsed']]-start),file.path(m4_out,'run-manifest.json'))
cat('M4 overall',overall,'passing datasets',np,'\n')
