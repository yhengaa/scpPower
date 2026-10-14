# Copy registered evidence and package it; never call calibration/simulation.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
start<-proc.time()[['elapsed']];out<-'paper/SI';dir.create(out,recursive=TRUE,showWarnings=FALSE)
entries<-list();inputs<-list()
sha<-function(p)digest::digest(file=p,algo='sha256')
copy_evidence<-function(source,subdir,name=basename(source)){
 dest<-file.path(out,subdir,name);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
 stopifnot(file.copy(source,dest,overwrite=TRUE),identical(sha(source),sha(dest)))
 entries[[length(entries)+1L]]<<-list(file=sub(paste0('^',out,'/'),'',dest),sha256=sha(dest),source_file=source,
  source_sha256=sha(source),copy_identical=TRUE)
 inputs[[source]]<<-sha(source);dest
}
read_report<-function(p){inputs[[p]]<<-sha(p);read.csv(p,check.names=FALSE)}
for(stage in c('M2','M2-revised','M2b'))for(name in c('fidelity-table.csv','fidelity-summary.csv'))
 copy_evidence(paste0('reports/',stage,'/',name),'S1-fidelity',paste0(stage,'-',name))
for(stage in c('M4','M4b')){
 for(name in c('validation-table.csv','dataset-verdicts.csv','run-manifest.json','predictions-frozen.json'))
  copy_evidence(paste0('reports/',stage,'/',name),'S2-validation',paste0(stage,'-',name))
 for(p in list.files(paste0('reports/',stage),pattern='-predicted-observed\\.png$',full.names=TRUE))
  copy_evidence(p,'S2-validation',paste0(stage,'-',basename(p)))
}
for(name in c('decision-table.csv','mask-summary-by-tier.csv','native-depth-summary.csv','pipeline-audit.csv','diagnostic-report.md'))
 copy_evidence(paste0('reports/J/',name),'S3-diagnostic-J')
for(name in c('power-grid.csv','realised-fdr.csv','instrument-hours.csv','cheapest-designs.csv','replicate-metrics.csv.gz','run-manifest.json'))
 copy_evidence(paste0('reports/M3b/',name),'S4-power-grid')
for(n in c('leduc2022_pSCoPE','derks2022_timsTOF','brunner2022')){
 copy_evidence(paste0('reports/M5/',n,'-worked-example-manifest.json'),'S5-worked-examples')
 copy_evidence(paste0('reports/M5/',n,'-worked-example.csv'),'S5-worked-examples')
}
copy_evidence('inst/doc/pilot-protocol.md','S6-pilot-protocol')
copy_evidence('scpPower_mandate_v1.1.md','S7-registered-record')
# Complete the requested rejection-reporting columns from already recorded
# replicate counts. Original grid bytes are retained beside this derivative.
grid_path<-'reports/M3b/power-grid.csv';rep_path<-'reports/M3b/replicate-metrics.csv.gz'
g<-read_report(grid_path);r<-read_report(rep_path)
rejections<-aggregate(rejected~task,r,mean);zero<-aggregate(I(rejected==0)~task,r,mean)
mean_rejections<-rejections$rejected[match(g$task,rejections$task)]
zero_rejections<-zero[[2]][match(g$task,zero$task)]
known<-is.finite(g$mean_rejections)
stopifnot(all(abs(g$mean_rejections[known]-mean_rejections[known])<1e-12))
known_zero<-is.finite(g$zero_rejection_fraction)
stopifnot(all(abs(g$zero_rejection_fraction[known_zero]-zero_rejections[known_zero])<1e-12))
g$mean_rejections[!known]<-mean_rejections[!known]
g$zero_rejection_fraction[!known_zero]<-zero_rejections[!known_zero]
g$rejection_summary_source<-ifelse(known,'Recorded power-grid columns','Derived from recorded per-replicate rejection counts')
dest<-file.path(out,'S4-power-grid/power-grid-complete.csv');write.csv(g,dest,row.names=FALSE,na='NA')
entries[[length(entries)+1L]]<-list(file='S4-power-grid/power-grid-complete.csv',sha256=sha(dest),copy_identical=FALSE,
 sources=list(list(source_file=grid_path,source_sha256=sha(grid_path)),list(source_file=rep_path,source_sha256=sha(rep_path))),
 expression='Fill unavailable mean_rejections with mean(rejected) by task; fill zero_rejection_fraction with mean(rejected == 0) by task. All existing values, FDR and hours unchanged.')
# Supporting captions are copied from the supplied manuscript, not rewritten.
manuscript<-readLines('paper/manuscript.md',warn=FALSE,encoding='UTF-8')
a<-match('# Supporting Information',manuscript);b<-which(seq_along(manuscript)>a&grepl('^# ',manuscript))[1]
si<-paste(manuscript[(a+1):(b-1)],collapse=' ')
si<-trimws(si);starts<-gregexpr('S[1-7],',si,perl=TRUE)[[1]];ends<-c(starts[-1]-1,nchar(si));parts<-vapply(seq_along(starts),function(i)trimws(substr(si,starts[i],ends[i])),character(1))
stopifnot(length(parts)==7)
writeLines(parts,'paper/SI-captions.txt')
mdtable<-function(x){z<-x;z[is.na(z)]<-'—';c(paste0('| ',paste(names(z),collapse=' | '),' |'),
 paste0('| ',paste(rep('---',ncol(z)),collapse=' | '),' |'),apply(z,1,function(v)paste0('| ',paste(v,collapse=' | '),' |')))}
label<-c(specht2019v3='Specht TMT',leduc2022_pSCoPE='Leduc TMT',derks2022_timsTOF='Derks timsTOF',
 derks2022_QExactive='Derks Q Exactive',brunner2022='Brunner label-free')
# One overview table per section; full tables, plots and manifests remain in SI/.
fg<-read_report('reports/M2-revised/fidelity-summary.csv');fi<-read_report('reports/M2b/fidelity-summary.csv')
old<-read_report('reports/M2/fidelity-summary.csv');final<-rbind(fi,fg[!fg$dataset %in% fi$dataset,])
final<-final[match(names(label),final$dataset),];old<-old[match(names(label),old$dataset),]
t1<-data.frame(Series=unname(label[final$dataset]),Final=final$status,
 Checks=paste0(final$checks_passed,'/',final$checks_applicable),Earlier=old$status,
 Earlier_checks=paste0(old$checks_passed,'/',old$checks_applicable))
v0<-read_report('reports/M4/dataset-verdicts.csv');vb<-read_report('reports/M4b/dataset-verdicts.csv')
t2<-do.call(rbind,lapply(names(label),function(n){x<-v0[v0$dataset==n,];y<-vb[vb$dataset==n,]
 data.frame(Series=label[[n]],First=if(x$status=='NOT EVALUABLE')'Not evaluable' else paste0(x$within_10_pp_points,'/',x$counted_points),
  First_max_gap_pp=if(is.finite(x$maximum_absolute_gap))sprintf('%.2f',100*x$maximum_absolute_gap) else '—',
  Rerun=if(nrow(y))paste0(y$within_10_pp_points,'/',y$counted_points) else 'Not rerun',
  Rerun_max_gap_pp=if(nrow(y))sprintf('%.2f',100*y$maximum_absolute_gap) else '—')}))
j<-read_report('reports/J/decision-table.csv');t3<-data.frame(Series=unname(label[j$dataset]),
 PSM_finite=sprintf('%.4f',j$psm_finite_fraction),Endpoint_finite=sprintf('%.4f',j$aggregated_finite_fraction),
 Gap_pp=sprintf('%.2f',100*j$aggregated_sparser_gap),Amendment_I=j$amendment_I_applies)
t4<-do.call(rbind,lapply(names(label),function(n)do.call(rbind,lapply(sort(unique(g$beta)),function(beta)
 do.call(rbind,lapply(if(n %in% names(label)[1:2])c(1,.75,.5) else 1,function(k){
  x<-g[g$dataset==n&abs(g$beta-beta)<1e-12&g$k==k,];x<-x[order(-x$power,x$instrument_hours),][1,]
  data.frame(Series=label[[n]],Fold=2^beta,k=k,Max_power=sprintf('%.4f',x$power),n=x$n_per_condition,m=x$m_cells_per_sample,
   Hours=sprintf('%.1f',x$instrument_hours))}))))))
t5<-do.call(rbind,lapply(c('leduc2022_pSCoPE','derks2022_timsTOF','brunner2022'),function(n){
 p<-paste0('reports/M5/',n,'-worked-example-manifest.json');inputs[[p]]<<-sha(p)
 x<-jsonlite::read_json(p,simplifyVector=FALSE)
 data.frame(Series=label[[n]],Seed=x$seed,Manifest_file=basename(p),SHA256_prefix=substr(sha(p),1,16))}))
pilot<-readLines('inst/doc/pilot-protocol.md',warn=FALSE)
t6<-data.frame(Item=c('Default design','Donor variance','Carrier compression','Background'),
 Record=c('4 donors × 2 conditions × 20 cells (default, not a derived minimum)',
 'Own cell type: log2 donor means in the reference condition',
 'Known-ratio two-proteome standard; k supplied to calibration','Blank wells in every acquisition batch'))
mandate<-readLines('scpPower_mandate_v1.1.md',warn=FALSE)
headings<-mandate[grepl('^## Frozen addendum',mandate)]
t7<-data.frame(Record=c('Frozen mandate v1.1',sub('^## ','',headings)),
 Deposited_file='scpPower_mandate_v1.1.md',SHA256_prefix=substr(sha('scpPower_mandate_v1.1.md'),1,16))
tables<-list(t1,t2,t3,t4,t5,t6,t7)
head<-c('---','title: "scpPower Supporting Information"','author: "Youshan Heng"','geometry: margin=0.65in','fontsize: 9pt','---','',
 '# Contents','',paste0('- S',1:7,'. ',c('Fidelity','Validation','Diagnostic J','Power grid','Worked examples','Pilot protocol','Registered record')),
 '', 'The full CSV tables, original plots, manifests and registered record are deposited in the accompanying SI directory. This PDF contains one overview table per section; SHA-256 checksums and source paths are in SI/manifest.json.', '', '\\newpage')
body<-head
for(i in seq_len(7))body<-c(body,paste0('# S',i), '',parts[i],'',mdtable(tables[[i]]),'',
 if(i==4)'Rejection summaries absent from older non-TMT grid columns are computed from their existing replicate counts. The original grid is retained unchanged; no runs were repeated.' else '',
 if(i<7)c('','\\newpage','') else '')
mdpath<-tempfile(fileext='.md');writeLines(body,mdpath)
status<-system2('pandoc',c(shQuote(mdpath),'-o',shQuote('paper/SI.pdf'),'--pdf-engine=pdflatex','-V','colorlinks=true'),stdout='paper/SI-render.log',stderr='paper/SI-render.log')
if(status!=0)stop('Supporting PDF rendering failed; see paper/SI-render.log')
unlink(mdpath)
reference_manifest<-'reports/M3b/run-manifest.json';inputs[[reference_manifest]]<-sha(reference_manifest)
reference_versions<-jsonlite::read_json(reference_manifest,simplifyVector=FALSE)$versions
versions<-list(R=as.character(getRversion()),scp=as.character(utils::packageVersion('scp')),limma=as.character(utils::packageVersion('limma')),edgeR=as.character(utils::packageVersion('edgeR')),scplainer=paste0('workflow in scp ',utils::packageVersion('scp'),'; no separate package'))
for(n in names(versions))if(!identical(versions[[n]],reference_versions[[n]]))warning('Version mismatch for ',n,': run=',versions[[n]],'; registered=',reference_versions[[n]],'; continuing.')
manifest<-list(files=entries,input_hashes=inputs,manuscript_sha256=sha('paper/manuscript.md'),seed=NULL,
 versions=versions,registered_versions=reference_versions,wall_seconds=proc.time()[['elapsed']]-start,new_scientific_runs=FALSE,
 checksum_scope='Every SI payload file; manifest.json excludes its own recursive checksum.',
 exceptions='S6 source is inst/doc/pilot-protocol.md; S7 source is the deposited frozen mandate. Both are byte-identical copies. The explicitly labelled complete grid is a derivative of registered replicate counts.',
 combined_pdf=list(file='paper/SI.pdf',sha256=sha('paper/SI.pdf'),supporting_captions=parts))
jsonlite::write_json(manifest,file.path(out,'manifest.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
cat('Supporting Information copied, hashed and rendered without scientific reruns.\n')
