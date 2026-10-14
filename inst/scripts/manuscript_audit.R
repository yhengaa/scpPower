# Audit the supplied prose against registered outputs, without running any model.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
start<-proc.time()[['elapsed']];dir.create('paper',showWarnings=FALSE)
manuscript_path<-'paper/manuscript.md'
args<-commandArgs(trailingOnly=TRUE);if(length(args))manuscript_path<-args[1]
text<-readLines(manuscript_path,warn=FALSE,encoding='UTF-8')
inputs<-list();read_report<-function(path,json=FALSE){
 stopifnot(startsWith(path,'reports/'));inputs[[path]]<<-digest::digest(file=path,algo='sha256')
 if(json)jsonlite::read_json(path,simplifyVector=FALSE) else read.csv(path,check.names=FALSE)
}
ids<-c('specht2019v3','leduc2022_pSCoPE','derks2022_timsTOF','derks2022_QExactive','brunner2022')
calpath<-'reports/M1-final/calibration-table.csv';cal<-read_report(calpath)
m1path<-'reports/M1-final/run-manifest.json';m1<-read_report(m1path,TRUE)
protpaths<-setNames(paste0('reports/M1-final/',ids,'-proteins.csv'),ids)
proteins<-setNames(lapply(protpaths,read_report),ids)
summarypath<-'reports/M1-final/dataset-summary.csv';summary<-read_report(summarypath)
gpath<-'reports/M3b/power-grid.csv';g<-read_report(gpath)
cp<-'reports/M3b/cheapest-designs.csv';cheap<-read_report(cp)
ip<-'reports/M3b/instrument-preset-manifest.json';instruments<-read_report(ip,TRUE)
rep_path<-'reports/M3b/replicate-metrics.csv.gz';rep<-read_report(rep_path)
rej<-aggregate(rejected~task,rep,mean);names(rej)[2]<-'mean_rejections_derived'
zero<-aggregate(I(rejected==0)~task,rep,mean);names(zero)[2]<-'zero_rejection_fraction_derived'
g$mean_rejections_derived<-rej$mean_rejections_derived[match(g$task,rej$task)]
g$zero_rejection_fraction_derived<-zero$zero_rejection_fraction_derived[match(g$task,zero$task)]
fgp<-'reports/M2-revised/fidelity-summary.csv';fi_gp<-read_report(fgp)
fip<-'reports/M2b/fidelity-summary.csv';fi_i<-read_report(fip)
fi<-rbind(fi_i,fi_gp[fi_gp$dataset %in% ids[3:5],])
old_fp<-'reports/M2/fidelity-summary.csv';fi_old<-read_report(old_fp)
ftg<-'reports/M2-revised/fidelity-table.csv';fit_g<-read_report(ftg)
fti<-'reports/M2b/fidelity-table.csv';fit_i<-read_report(fti)
fit<-rbind(fit_i,fit_g[fit_g$dataset %in% ids[3:5],])
v0p<-'reports/M4/dataset-verdicts.csv';v0<-read_report(v0p)
vbp<-'reports/M4b/dataset-verdicts.csv';vb<-read_report(vbp)
v0tp<-'reports/M4/validation-table.csv';v0t<-read_report(v0tp)
vbtp<-'reports/M4b/validation-table.csv';vbt<-read_report(vbtp)
hwp<-'reports/M4b/high-power-window.csv';hw<-read_report(hwp)
jp<-'reports/J/decision-table.csv';j<-read_report(jp)
depthp<-'reports/J/native-depth-summary.csv';depth<-read_report(depthp)
expath<-'reports/M5/worked-examples.csv';examples<-read_report(expath)
m3path<-'reports/M3b/run-manifest.json';m3<-read_report(m3path,TRUE)
sspath<-'reports/M4/source-summary.csv';ss<-read_report(sspath)
registration_path<-'reports/M4/registration.json';registration<-read_report(registration_path,TRUE)
# A reference always names a registered source and the actual selection/expression.
ref<-function(value,file,field,exact=FALSE)list(value=value,file=file,field=field,exact=exact)
C<-function(n,field){x<-cal[cal$dataset==n,];ref(unique(x[[field]]),calpath,paste0('dataset=',n,'; ',field,' (identical across tiers)'))}
P<-function(n,field,operation='median_sqrt'){
 x<-proteins[[n]][[field]]
 value<-switch(operation,median_sqrt=median(sqrt(x),na.rm=TRUE),mean=mean(x,na.rm=TRUE),mapped=sum(is.finite(x)),coverage=sum(x,na.rm=TRUE)/C(n,'cells')$value)
 ref(value,protpaths[[n]],paste0('dataset=',n,'; ',switch(operation,median_sqrt=paste0('median(sqrt(',field,'), na.rm=TRUE)'),mean=paste0('mean(',field,', na.rm=TRUE)'),mapped=paste0('sum(is.finite(',field,'))'),coverage=paste0('sum(',field,', na.rm=TRUE) / cells=',C(n,'cells')$value))))
}
G<-function(n,beta=1,n_samples=12,m_cells=64,k=1,field='power',maximise=FALSE){
 x<-g[g$dataset==n&abs(g$beta-beta)<1e-12&g$k==k,,drop=FALSE]
 if(maximise)x<-x[order(-x$power,x$instrument_hours,x$cells,x$n_per_condition),][1,,drop=FALSE]
 else x<-x[x$n_per_condition==n_samples&x$m_cells_per_sample==m_cells,,drop=FALSE]
 stopifnot(nrow(x)==1)
 ref(if(field=='fold_change')2^x$beta else if(field=='power_testability_ratio')x$power/x$testability else x[[field]],gpath,
  paste0('task=',x$task,'; dataset=',n,'; beta=',format(beta,digits=16),'; k=',k,'; ',
   if(field=='fold_change')'2^beta' else if(field=='power_testability_ratio')'power/testability' else field,
   if(maximise)'; row selected by maximum power, then minimum hours/cells/n' else ''))
}
E<-function(n,n_samples,field){x<-examples[examples$dataset==n&examples$n_per_condition==n_samples,];stopifnot(nrow(x)==1)
 ref(x[[field]],expath,paste0('dataset=',n,'; n_per_condition=',n_samples,'; m_cells_per_sample=32; ',field))}
V<-function(n,field,rerun=TRUE,mult=1){x<-if(rerun)vb else v0;path<-if(rerun)vbp else v0p
 ref(x[x$dataset==n,field]*mult,path,paste0('dataset=',n,'; ',field,if(mult!=1)paste0(' * ',mult) else ''))}
J<-function(n,field){ref(j[j$dataset==n,field]*100,jp,paste0('dataset=',n,'; tier=0; ',field,' * 100'))}
D<-function(n){x<-depth[depth$dataset==n&depth$mask=='aggregated',];ref(x$observed_depth_sd/x$simulated_mean_depth_sd,depthp,paste0('dataset=',n,'; mask=aggregated; observed_depth_sd/simulated_mean_depth_sd'))}
H<-function(n){ref(hw$mean_signed_gap[hw$dataset==n]*100,hwp,paste0('dataset=',n,'; mean_signed_gap * 100; max(predicted,observed)>=0.2'))}
range_ref<-function(values,file,expression,bound='min')ref(if(bound=='min')min(values,na.rm=TRUE) else max(values,na.rm=TRUE),file,paste0(bound,'(',expression,')'))
all_files<-function(paths)paste(paths,collapse='; ')
# A direct rounding comparison; absent evidence is FALSE, never a guessed match.
rows<-list();issues<-list();coverage<-list();claim_no<-0L
number_words<-c(one=1,two=2,three=3,four=4,five=5,six=6,eight=8,ten=10,twelve=12,twenty=20,quarter=.25,halves=.5)
num_regex<-'(?i)(?<![A-Za-z0-9_])[-+−]?(?:[0-9]{1,3}(?:,[0-9]{3})+|[0-9]+)(?:\\.[0-9]+)?|\\b(?:one|two|three|four|five|six|eight|ten|twelve|twenty|quarter|halves)\\b'
num_value<-function(token){z<-tolower(gsub('[,%]','',gsub('−','-',token,fixed=TRUE)));if(z %in% names(number_words))number_words[[z]] else suppressWarnings(as.numeric(z))}
precision<-function(token){z<-gsub(',','',token);if(tolower(z) %in% c('quarter'))return(2L)
 if(tolower(z)=='halves')return(1L)
 if(grepl('\\.',z))nchar(sub('.*\\.','',z)) else 0L}
record<-function(section,quote,value,r,why=NULL){
 claim_no<<-claim_no+1L
 val<-num_value(value)
 source<-r$value
 if(is.numeric(source)&&length(source)==1&&!is.na(source)&&!is.na(val)) {
  match<-if(isTRUE(r$exact))abs(val-source)<1e-12 else abs(val-round(source,precision(value)))<1e-10
 } else match<-length(source)==1&&!is.na(source)&&identical(as.character(value),as.character(source))
 if(!is.null(why))match<-FALSE
 rows[[claim_no]]<<-data.frame(section=section,short_quote=quote,manuscript_value=as.character(value),source_file=r$file,
  source_field=r$field,source_value=if(length(source)==1)as.character(source) else paste(source,collapse='; '),match=match,stringsAsFactors=FALSE)
 if(!match)issues[[length(issues)+1L]]<<-list(claim=claim_no,section=section,quote=quote,value=value,source=source,file=r$file,field=r$field,
  why=if(!is.null(why))why else if(all(is.na(source)))'Not traceable from the registered source; no value was inferred.' else 'Does not equal the registered value after the displayed rounding.')
}
# Strip layout identifiers only. Their exclusions are retained in the coverage file.
tokens<-function(paragraph){
 p<-paragraph
 structural<-c('\\^([0-9]+[–,]?)+\\^','(?i)(?:Figure|Table)\\s+[0-9]+[a-e]?(?:[, ]?[a-e])?',
  'GPL-3','log2','two-dimensional','first attempt','first version')
 for(pattern in structural){
  h<-gregexpr(pattern,p,perl=TRUE)[[1]];len<-attr(h,'match.length')
  if(h[1]>0)for(i in rev(seq_along(h)))substr(p,h[i],h[i]+len[i]-1)<-strrep(' ',len[i])
 }
 hits<-gregexpr(num_regex,p,perl=TRUE)[[1]]
 if(hits[1]<0)return(character())
 ts<-regmatches(p,list(hits))[[1]];attr(ts,'positions')<-hits;ts
}
# Numeric claims are mapped by sentence context and position, not by searching
# for a conveniently equal number somewhere else in the reports.
audit_paragraph<-function(section,paragraph,refs){
 ts<-tokens(paragraph)
 if(length(ts)!=length(refs))stop('Incomplete numeric mapping in ',section,': ',paste(ts,collapse=' | '),' ; references=',length(refs))
 for(i in seq_along(ts)){
  r<-refs[[i]]
  if(is.null(r))r<-ref(NA_real_,gpath,'Claim has no registered numerical source; see the quoted context')
  record(section,substr(paragraph,max(1,attr(ts,'positions')[i]-55),min(nchar(paragraph),attr(ts,'positions')[i]+nchar(ts[i])+65)),ts[i],r)
  coverage[[length(coverage)+1L]]<<-data.frame(section=section,token_index=i,token=ts[i],audit_claim=claim_no,paragraph=paragraph)
 }
}
# Parse headings and paragraphs; original manuscript bytes are not changed.
sections<-list();current<-''
for(line in text){
 if(grepl('^#{1,2} ',line)){current<-sub('^#+ ','',line);next}
 if(nzchar(trimws(line))&&!startsWith(line,'|')&&!startsWith(line,'!')&&!startsWith(line,'```')&&!startsWith(line,'<w:'))
  sections[[current]]<-c(sections[[current]],line)
}
writeLines(vapply(c('Abstract','Calibration','Simulator fidelity','Registered design grid','Validation on real cells','Worked designs','Discussion'),function(s)
 paste0(s,': ',paste(vapply(sections[[s]],function(p)paste(tokens(p),collapse=' | '),character(1)),collapse='\n')),character(1)),'paper/audit-token-inventory.txt')
# References for prose claim groups, in the exact token order emitted above.
constant<-function(v,field,file=m3path)ref(v,file,field)
preset_count<-ref(length(unique(fi$dataset)),all_files(c(fgp,fip)),'number of distinct final fidelity dataset series')
dataset_count<-ref(length(unique(sub('_.*','',sub('derks2022.*','derks2022',ids)))),m1path,'number of public datasets after collapsing the two Derks instrument series')
tmt_count<-ref(nrow(vb),vbp,'number of carrier-TMT series in the rerun verdict table')
fold<-function(v)ref(2^unique(g$beta[abs(g$beta-log2(v))<1e-12]),gpath,paste0('2^beta for registered fold change=',v))
samples<-function(v)ref(unique(g$n_per_condition[g$n_per_condition==v]),gpath,paste0('unique(n_per_condition[n_per_condition==',v,'])'))
cells<-function(v)ref(unique(g$m_cells_per_sample[g$m_cells_per_sample==v]),gpath,paste0('unique(m_cells_per_sample[m_cells_per_sample==',v,'])'))
kref<-function(v)ref(unique(g$k[g$k==v]),gpath,paste0('unique(k[k==',v,'])'))
threshold<-constant(m3$grid_registration$target_power,'grid_registration.target_power',m3path)
minimum_cells<-constant(m3$endpoint$minimum_quantified_cells_per_sample,'endpoint.minimum_quantified_cells_per_sample')
condition_count<-ref(unique(g$cells/(g$n_per_condition*g$m_cells_per_sample)),gpath,'unique(cells / (n_per_condition*m_cells_per_sample))')
max_range<-function(beta,field='power',bound='min'){
 x<-g[g$n_per_condition==12&g$m_cells_per_sample==64&g$k==1&abs(g$beta-beta)<1e-12,]
 range_ref(x[[field]],gpath,paste0(field,' for n=12,m=64,k=1,beta=',format(beta,digits=16),' across all presets'),bound)
}
medians<-function(ids,field,bound){values<-sapply(ids,function(n)P(n,field)$value)
 range_ref(values,all_files(protpaths[ids]),paste0('median(sqrt(',field,'),na.rm=TRUE) per indicated series: ',paste(ids,collapse=', ')),bound)}
compression_ratio<-ref(G(ids[1],k=.5)$value/G(ids[1],k=1)$value,gpath,
 'Specht n=12,m=64,beta=1: power(k=0.5)/power(k=1); the stated halving is checked, not assumed')
fidelity_max<-function(field){ix<-which.max(fi[[field]])
 ref(fi[[field]][ix]*100,if(fi$dataset[ix] %in% ids[1:2])fip else fgp,paste0('maximum ',field,' * 100 over final M2b TMT and M2-revised non-TMT rows; attaining dataset=',fi$dataset[ix]))}
fidelity_tolerance<-function(metric){x<-fit[fit$metric==metric,];ref(unique(x$tolerance)*100,all_files(c(ftg,fti)),paste0('unique(tolerance[metric=',metric,']) * 100'))}
qmax<-ref(max(v0t$observed_power[v0t$dataset==ids[4]&v0t$counted]),v0tp,'max(observed_power) for counted Q Exactive points')
# The claimed "below 2%" is an inequality, recorded as the actual maximum.
less_two<-ref(100*qmax$value,v0tp,'max(observed_power) * 100 for counted Q Exactive points; stated upper bound 2%')
less_two$bound<-2
atleast10<-ref(10,rep_path,'mean(rejected) by task; query threshold >=10 is the manuscript-specified comparison threshold')
fd_good<-ref(max(g$realised_fdr[g$mean_rejections_derived>=10]),all_files(c(gpath,rep_path)),
 'max(realised_fdr) where mean(rejected) by task >=10; includes all 945 grid points')
fd_low<-ref(max(g$realised_fdr),gpath,'max(realised_fdr) over all 945 registered grid points; task attaining maximum is 46 (Specht)')
# Handle inequality claims by their mathematical predicate, not rounding a bound.
original_record<-record
record<-function(section,quote,value,r,why=NULL){
 if(!is.null(r$bound)){
  original_record(section,quote,value,r,why)
  i<-length(rows);rows[[i]]$match<<-is.numeric(r$value)&&is.finite(r$value)&&r$value<r$bound
  if(rows[[i]]$match)issues<<-Filter(function(z)z$claim!=i,issues)
 }else original_record(section,quote,value,r,why)
}
A<-sections[['Abstract']][1]
audit_paragraph('Abstract',A,list(preset_count,dataset_count,V(ids[1],'maximum_absolute_gap',mult=100),tmt_count,
 V(ids[5],'maximum_absolute_gap',FALSE,100),fold(1.25),samples(12),compression_ratio))
p<-sections[['Calibration']][1]
audit_paragraph('Results / Calibration',p,list(preset_count,
 range_ref(sapply(ids,function(n)P(n,'quantified_cells','coverage')$value),all_files(protpaths),'sum(quantified_cells)/cells per series','min'),
 range_ref(sapply(ids,function(n)P(n,'quantified_cells','coverage')$value),all_files(protpaths),'sum(quantified_cells)/cells per series','max'),
 range_ref(sapply(ids,function(n)P(n,'missing_fraction','mean')$value),all_files(protpaths),'mean(endpoint missing_fraction) per series','min'),
 range_ref(sapply(ids,function(n)P(n,'missing_fraction','mean')$value),all_files(protpaths),'mean(endpoint missing_fraction) per series','max'),
 medians(ids[1:2],'sigma2_e','min'),medians(ids[1:2],'sigma2_e','max'),medians(ids[3:4],'sigma2_e','min'),medians(ids[3:4],'sigma2_e','max'),
 P(ids[1],'sigma2_b'),P(ids[2],'sigma2_b'),medians(ids,'sigma2_s','min'),medians(ids,'sigma2_s','max'),fold(1.25),
 ref(log2(1.25),gpath,'unique(beta) at fold_change=1.25')))
p<-sections[['Simulator fidelity']][1]
audit_paragraph('Results / Simulator fidelity',p,list(preset_count,fidelity_max('maximum_missingness_error'),fidelity_max('coverage_relative_error'),
 fidelity_max('maximum_cell_sd_error'),fidelity_max('maximum_batch_sd_error'),fidelity_tolerance('missing_fraction'),
 fidelity_tolerance('mean_proteins_per_cell'),fidelity_tolerance('cell_sd'),fidelity_tolerance('batch_sd'),
 ref(sum(fi_old$status=='FAIL'),old_fp,'sum(status==FAIL)'),ref(nrow(fi_old),old_fp,'number of series'),
 range_ref(fi_old$maximum_cell_sd_error[fi_old$status=='FAIL']*100,old_fp,'maximum_cell_sd_error *100 over failed series','min'),
 range_ref(fi_old$maximum_cell_sd_error[fi_old$status=='FAIL']*100,old_fp,'maximum_cell_sd_error *100 over failed series','max')))
p<-sections[['Registered design grid']][1]
x<-g[g$n_per_condition==12&g$m_cells_per_sample==64&g$beta==1&g$k==1,]
audit_paragraph('Results / Registered design grid',p,list(NULL,fold(1.25),samples(12),cells(64),max_range(log2(1.25),'power','min'),max_range(log2(1.25),'power','max'),
 fold(1.5),max_range(log2(1.5),'power','min'),max_range(log2(1.5),'power','max'),fold(2),samples(12),cells(64),
 range_ref(x$power/x$testability,gpath,'power/testability at n=12,m=64,beta=1,k=1 across all five presets','min'),
 range_ref(x$power/x$testability,gpath,'power/testability at n=12,m=64,beta=1,k=1 across all five presets','max'),minimum_cells,
 ref(1,m3path,'testability criterion requires the minimum in every sample; failing any one makes a protein untestable'),
 ref(2*12,gpath,'cells/m_cells_per_sample for n=12 (both conditions)'),G(ids[1],field='testability'),G(ids[2],field='testability'),cells(64),
 G(ids[1]),G(ids[2]),compression_ratio,fold(2),samples(12),cells(64),G(ids[1]),G(ids[1],k=.75),G(ids[1],k=.5),
 G(ids[2]),G(ids[2],k=.75),G(ids[2],k=.5),kref(1),kref(.75),kref(.5)))
p<-sections[['Registered design grid']][2]
win<-g[g$beta==1&g$power>=.8,];win<-win[win$instrument_hours==min(win$instrument_hours),]
audit_paragraph('Results / Registered design grid',p,list(ref(nrow(win),gpath,'number of beta=1, power>=0.8 rows attaining minimum hours'),threshold,fold(2),
 G(ids[3],m_cells=32,field='n_per_condition'),G(ids[3],m_cells=32,field='m_cells_per_sample'),G(ids[3],m_cells=32),G(ids[3],m_cells=32,field='instrument_hours'),
 G(ids[5]),cells(64),G(ids[5],field='instrument_hours'),threshold,fd_good,atleast10,fd_low))
record('Results / Registered design grid','Label-free DIA reaches 0.80 only at 64 cells per sample (896 hours)',
 'target reached',ref(if(any(g$dataset==ids[5]&g$power>=.8))'target reached' else 'none in grid',cp,
 'cheapest-designs row dataset=brunner2022,beta=1,k=1: target status; compare unrounded power >=0.8'))
p<-sections[['Validation on real cells']][1]
audit_paragraph('Results / Validation on real cells',p,list(V(ids[5],'within_10_pp_points',FALSE),V(ids[5],'counted_points',FALSE),
 ref(registration$criterion$absolute_power_gap*100,registration_path,'criterion.absolute_power_gap * 100'),V(ids[5],'maximum_absolute_gap',FALSE,100),
 V(ids[4],'counted_points',FALSE),ref(ss$cells[ss$dataset==ids[4]],sspath,'dataset=derks2022_QExactive; cells'),less_two,
 ref(ss$cells[ss$dataset==ids[3]],sspath,'dataset=derks2022_timsTOF; cells')))
p<-sections[['Validation on real cells']][2]
big<-vbt[vbt$dataset==ids[1]&vbt$n_per_condition==8&vbt$m_cells_per_sample==64&vbt$beta==1,]
audit_paragraph('Results / Validation on real cells',p,list(tmt_count,V(ids[1],'within_10_pp_points',FALSE),V(ids[1],'counted_points',FALSE),
 V(ids[2],'within_10_pp_points',FALSE),V(ids[2],'counted_points',FALSE),
 range_ref(v0$maximum_absolute_gap[v0$dataset %in% ids[1:2]]*100,v0p,'maximum_absolute_gap*100 over the two first-attempt TMT rows','max'),
 J(ids[1],'aggregated_sparser_gap'),J(ids[2],'aggregated_sparser_gap'),V(ids[1],'within_10_pp_points'),V(ids[1],'counted_points'),
 V(ids[2],'within_10_pp_points'),V(ids[2],'counted_points'),ref(registration$criterion$absolute_power_gap*100,registration_path,'criterion.absolute_power_gap * 100 (same rerun rule)'),
 V(ids[1],'maximum_absolute_gap',mult=100),V(ids[2],'maximum_absolute_gap',mult=100),H(ids[1]),H(ids[2]),
 ref(unique(hw$cutoff),hwp,'unique(cutoff), max(predicted,observed) selection'),
 ref(big$n_per_condition,vbtp,'Specht n=8,m=64,beta=1; n_per_condition'),ref(big$m_cells_per_sample,vbtp,'Specht n=8,m=64,beta=1; m_cells_per_sample'),
 ref(2^big$beta,vbtp,'Specht n=8,m=64,beta=1; 2^beta'),ref(big$predicted_power,vbtp,'Specht n=8,m=64,beta=1; predicted_power'),
 ref(big$observed_power,vbtp,'Specht n=8,m=64,beta=1; observed_power'),D(ids[1]),D(ids[2])))
p<-sections[['Worked designs']][1]
audit_paragraph('Results / Worked designs',p,list(fold(2),ref(unique(examples$m_cells_per_sample),expath,'unique(m_cells_per_sample)'),samples(8),samples(12),
 E(ids[2],8,'power'),E(ids[2],12,'power'),E(ids[2],8,'testability'),E(ids[2],12,'testability'),
 E(ids[2],8,'conditional_power'),E(ids[2],12,'conditional_power'),E(ids[2],8,'instrument_hours'),E(ids[2],12,'instrument_hours'),
 E(ids[3],8,'power'),E(ids[3],12,'power'),E(ids[3],8,'testability'),E(ids[3],12,'testability'),E(ids[3],8,'instrument_hours'),E(ids[3],12,'instrument_hours'),
 E(ids[5],8,'power'),E(ids[5],12,'power'),E(ids[5],8,'testability'),E(ids[5],12,'testability'),E(ids[5],8,'instrument_hours'),E(ids[5],12,'instrument_hours'),
 ref(2,expath,'Two cited series: Leduc TMT and Brunner label-free'),NULL,threshold))
# Discussion values inherit the same selections; the quarter-hours claim is
# checked against both stated registered design times, not rounded to a slogan.
p<-sections[['Discussion']][1]
audit_paragraph('Discussion',p,list(condition_count,NULL,fold(1.25),samples(12),
 ref(2,calpath,'Two named carrier-TMT series'),fold(2),threshold,
 ref(G(ids[3],m_cells=32,field='instrument_hours')$value/G(ids[5],field='instrument_hours')$value,gpath,
 'timsTOF n=12,m=32,beta=1,k=1 hours / Brunner n=12,m=64,beta=1,k=1 hours; Brunner does not attain the target')))
p<-sections[['Discussion']][2]
audit_paragraph('Discussion',p,list(ref(1,m1path,'One selected cohort: cohort.accession=PXD058573'),
 range_ref(vbt$observed_power[vbt$counted],vbtp,'observed_power over counted TMT rerun points','max'),
 range_ref(v0t$observed_power[v0t$counted&v0t$dataset==ids[5]],v0tp,'observed_power over counted Brunner points','max'),
 ref(1,m3path,'One registered primary endpoint: endpoint.design'),condition_count))
pilot_pdf<-'reports/M5/pilot-protocol.pdf';inputs[[pilot_pdf]]<-digest::digest(file=pilot_pdf,algo='sha256')
pilot_text<-paste(system2('pdftotext',c('-layout',shQuote(pilot_pdf),'-'),stdout=TRUE),collapse=' ')
stopifnot(grepl('4 donors',pilot_text),grepl('2 conditions',pilot_text),grepl('20 cells',pilot_text))
p<-sections[['Discussion']][3]
audit_paragraph('Discussion',p,list(ref(2,m1path,'Declared transfer_assumption and compression assumption'),
 ref(4,pilot_pdf,'Default-design paragraph: 4 donors'),ref(2,pilot_pdf,'Default-design paragraph: 2 conditions'),
 ref(20,pilot_pdf,'Default-design paragraph: 20 cells'),ref(2,m1path,'compression.pilot_recommendation: two-proteome standard')))
# Also audit the explicitly requested Methods donor counts/range and k values.
p<-sections[['Between-donor variance and carrier compression']][1]
audit_paragraph('Methods / donor transfer (requested additional scope)',p,list(ref(1,m1path,'One culture per line; datasets[*].biological_unit=culture'),
 ref(m1$cohort$healthy_donors,m1path,'cohort.healthy_donors'),
 range_ref(sapply(ids,function(n)P(n,'sigma2_s','mapped')$value),all_files(protpaths),'sum(is.finite(sigma2_s)) per series','min'),
 range_ref(sapply(ids,function(n)P(n,'sigma2_s','mapped')$value),all_files(protpaths),'sum(is.finite(sigma2_s)) per series','max'),
 kref(1),kref(.75),kref(.5)))
power_report_path<-'reports/M3b/power-report.md'
power_report_text<-paste(readLines(power_report_path,warn=FALSE),collapse=' ')
inputs[[power_report_path]]<-digest::digest(file=power_report_path,algo='sha256')
stopifnot(grepl('TMT',power_report_text),grepl('plexDIA',power_report_text),grepl('label-free',power_report_text))
# Tables: every populated Table 1 cell has a row, including nonnumeric labels.
# Publication-year/model metadata absent from registered outputs is explicitly
# unverified; a dataset identifier is never mistaken for a publication year.
table_lines<-text[startsWith(text,'|')]
parse_row<-function(line)trimws(strsplit(sub('\\|$','',sub('^\\|','',line)),'\\|')[[1]])
t1_lines<-table_lines[vapply(table_lines,function(z)length(parse_row(z))==12,logical(1))]
t1_lines<-t1_lines[!grepl('^\\| Series|^\\| [- ]',t1_lines)]
stopifnot(length(t1_lines)==5)
for(i in seq_along(ids)){
 n<-ids[i];v<-parse_row(t1_lines[i]);cols<-parse_row(table_lines[1])
 ipref<-function(field)ref(instruments$presets[[n]][[field]],ip,paste0('presets.',n,'.',field))
 # Typography normalization is allowed only for a recorded instrument name.
 rawinstrument<-m1$datasets[[n]]$instrument
 instrument<-if(n==ids[3])gsub('timsTOFSCP','timsTOF SCP',rawinstrument,fixed=TRUE) else if(n==ids[4])gsub('Q-Exactive','Q Exactive',rawinstrument,fixed=TRUE) else NA_character_
 refs<-list(ref(NA_character_,m1path,paste0('datasets.',n,'.dataset=',n,'; publication year is not stored in this registered field')),
  ref(instrument,m1path,paste0('datasets.',n,'.instrument=',paste(rawinstrument,collapse=','),'; typographic aliases normalized; missing model metadata remains unverified')),
  ref(if(i<=2)'carrier TMT' else if(i<=4)'plexDIA' else 'label-free DIA', 'reports/M3b/power-report.md',paste0('Registered power-report acquisition description for dataset=',n,'; text label normalized')),
  C(n,'cells'),C(n,'samples'),P(n,'quantified_cells','coverage'),P(n,'missing_fraction','mean'),P(n,'sigma2_e'),
  if(n==ids[5])ref(if(all(grepl('absorbed into sigma2_e',proteins[[n]]$sigma2_b_reason,fixed=TRUE)))'absorbed' else 'not absorbed',protpaths[[n]],'all(grepl(absorbed into sigma2_e, sigma2_b_reason)); registered sigma2_b is zero after absorption') else P(n,'sigma2_b'),
  P(n,'sigma2_s'),ipref('cells_per_run'),ipref('minutes_per_run'))
 for(jj in seq_along(v)){
  quote<-paste0('Table 1, ',v[1],', ',cols[jj],': ',v[jj])
  if(jj==10){ts<-tokens(v[jj]);stopifnot(length(ts)==2)
   record('Table 1',quote,ts[1],refs[[jj]]);record('Table 1',quote,ts[2],P(n,'sigma2_s','mapped'))
  }else if(jj==3&&n==ids[2]){
   record('Table 1',quote,v[jj],ref(NA_character_,m3path,'carrier TMT is registered; TMTpro 18 chemistry and channel count are not stored in these outputs'),
    'Labelling identity is registered, but the stated chemistry/channel-count metadata cannot be verified from the deposited outputs.')
  }else record('Table 1',quote,v[jj],refs[[jj]])
 }
}
# Table 2: select maximum-power rows independently of manuscript values.
t2_lines<-table_lines[vapply(table_lines,function(z)length(parse_row(z))==6,logical(1))]
t2_lines<-t2_lines[!grepl('^\\| Series|^\\| [- ]',t2_lines)]
stopifnot(length(t2_lines)==15)
for(i in seq_along(t2_lines)){
 n<-ids[ceiling(i/3)];v<-parse_row(t2_lines[i]);beta<-log2(as.numeric(v[2]));prefix<-paste0('Table 2, dataset=',n,', fold=',v[2],': ')
 record('Table 2',paste0(prefix,'fold change ',v[2]),v[2],G(n,beta=beta,field='fold_change',maximise=TRUE))
 record('Table 2',paste0(prefix,'maximum power ',v[3]),v[3],G(n,beta=beta,maximise=TRUE))
 design<-tokens(v[4]);stopifnot(length(design)==2)
 record('Table 2',paste0(prefix,'design ',v[4]),design[1],G(n,beta=beta,field='n_per_condition',maximise=TRUE))
 record('Table 2',paste0(prefix,'design ',v[4]),design[2],G(n,beta=beta,field='m_cells_per_sample',maximise=TRUE))
 record('Table 2',paste0(prefix,'hours ',v[5]),v[5],G(n,beta=beta,field='instrument_hours',maximise=TRUE))
 eligible<-g[g$dataset==n&abs(g$beta-beta)<1e-12&g$k==1&g$power>=threshold$value,]
 if(grepl('^none in grid',v[6])){
  record('Table 2',paste0(prefix,v[6]),'none in grid',ref(if(nrow(eligible)==0)'none in grid' else 'target reached',gpath,
   paste0('dataset=',n,',beta=',beta,',k=1; number of rows with unrounded power >= grid_registration.target_power=0.8 is ',nrow(eligible))))
  if(length(tokens(v[6]))){audit_paragraph('Table 2',paste0('Brunner maximum detail: ',v[6]),
   list(G(n,beta=beta),G(n,beta=beta,field='n_per_condition'),G(n,beta=beta,field='m_cells_per_sample')))}
 }else{
  eligible<-eligible[order(eligible$instrument_hours,eligible$cells,eligible$n_per_condition,-eligible$power),][1,,drop=FALSE]
  rr<-lapply(c('n_per_condition','m_cells_per_sample','power','instrument_hours'),function(field)
   ref(eligible[[field]],gpath,paste0('task=',eligible$task,'; minimum hours among unrounded power>=0.8 at dataset=',n,',beta=',beta,',k=1; ',field)))
  audit_paragraph('Table 2',paste0('Cheapest: ',v[6]),rr)
 }
}
# Table captions and figure captions, with layout identifiers excluded.
audit_paragraph('Table 2 caption',text[grepl('^\\*\\*Table 2\\.',text)],list(kref(1),threshold,kref(.75),kref(.5)))
caption<-function(i)text[grepl(paste0('^\\*\\*Figure ',i,'\\.'),text)]
# Model terms count includes four columns in the protein table plus beta in grid.
model_terms<-intersect(c('mean_log2','sigma2_s','sigma2_b','sigma2_e'),names(proteins[[1]]))
five_terms<-ref(length(model_terms)+as.integer('beta' %in% names(g)),all_files(c(protpaths[1],gpath)),
 'length(intersect(c(mean_log2,sigma2_s,sigma2_b,sigma2_e),protein columns)) + as.integer(beta in grid columns): five stated additive model terms')
two_detection<-ref(length(unique(ifelse(g$dataset %in% ids[1:2],'plex/channel','cell-baseline'))),gpath,
 'Two registered acquisition detection schemes over the grid: carrier TMT plex/channel and DIA cell-baseline')
audit_paragraph('Figure 1 caption',caption(1),list(five_terms,two_detection))
audit_paragraph('Figure 2 caption',caption(2),list(preset_count,fidelity_tolerance('missing_fraction')))
audit_paragraph('Figure 3 caption',caption(3),list(fold(2),fold(1.5),kref(1),
 ref(1,gpath,'one panel per unique(dataset); the number is a figure-layout choice, not an estimated parameter'),threshold,samples(12),fold(2),fold(2),samples(12),
 kref(1),kref(.75),kref(.5),fold(2),threshold))
# Axis limits are presentation metadata, not estimates from the registered runs.
# Retain them as untraceable rather than invent a data source for an artist choice.
audit_paragraph('Figure 4 caption',caption(4),list(ref(registration$criterion$absolute_power_gap*100,registration_path,'criterion.absolute_power_gap * 100'),
 ref(NA_real_,v0tp,'Q Exactive figure lower axis limit: display choice absent from registered validation outputs'),
 ref(NA_real_,v0tp,'Q Exactive figure upper axis limit: display choice absent from registered validation outputs'),
 less_two,ref(ss$cells[ss$dataset==ids[3]],sspath,'dataset=derks2022_timsTOF; cells'),minimum_cells))
qx<-v0[v0$dataset==ids[4],]
record('Results / Validation on real cells','The Q Exactive passes all six eligible points',
 '6/6',ref(paste0(qx$within_10_pp_points,'/',qx$counted_points),v0p,'dataset=derks2022_QExactive; paste(within_10_pp_points, counted_points, sep=/)'))
# Extra rounded target-attainment claim in the Discussion is not made true by
# rendering 0.796 as 0.80. Keep its logical predicate alongside numeric checks.
record('Discussion','At 2-fold, plexDIA on the timsTOF SCP reaches 0.8 at about a quarter of the instrument hours of label-free DIA',
 'label-free comparator reaches target',ref(if(any(g$dataset==ids[5]&g$power>=threshold$value))'label-free comparator reaches target' else 'label-free comparator never reaches target',gpath,
 'any(power>=0.8 for Brunner): no eligible label-free design reaches the target'))
record('Results / Registered design grid','The few values up to 0.068 occur at zero-power designs with a handful of rejections',
 'zero-power designs',ref(if(all(g$power[g$realised_fdr>.05]==0))'zero-power designs' else 'nonzero-power designs also exceed 0.05',gpath,
 'all(power[realised_fdr > 0.05] == 0); no display rounding applied to a zero-power assertion'))
# Final artifacts: a failure is evidence for manuscript review, not permission to
# change the registered source or manuscript. Only these new paper files are written.
out<-do.call(rbind,rows);stopifnot(all(vapply(strsplit(out$source_file,'; ',fixed=TRUE),function(z)all(startsWith(z,'reports/')),logical(1))))
write.csv(out,'paper/audit-table.csv',row.names=FALSE,na='NA')
write.csv(do.call(rbind,coverage),'paper/audit-coverage.csv',row.names=FALSE,na='NA')
escape<-function(x)gsub('|','\\|',paste(x,collapse='; '),fixed=TRUE)
lines<-c('# Manuscript audit mismatches','',
 'The supplied manuscript is unchanged. A FALSE match means either a numerical/logical discrepancy or missing registered provenance; it does not imply that unrecorded publication or instrument metadata is factually false.',
 '', 'Numeric comparisons use the precision displayed in the manuscript. Signed values retain their sign. Inequalities and target attainment are checked against unrounded values. Qualifiers such as “about” are not an unstated rounding tolerance.',
 '', 'Every populated Table 1 cell is checked. Numeric tokens in the requested prose, Table 2 and captions are checked in context; publication/citation identifiers, figure/table identifiers, log2 and layout markup are excluded. Methods donor-transfer numbers are included at the user’s request. Rhetorical counts and axis limits lacking a registered source are explicitly unverified.',
 '',paste0('Claims checked: ',nrow(out),'. FALSE matches: ',sum(!out$match),'. Numerical or logical discrepancies: ',sum(!out$match & out$source_value!='NA'),'. Missing registered provenance: ',sum(!out$match & out$source_value=='NA'),'.'),'',
 '| Claim | Section | Manuscript | Registered source value | Reason | Source / expression |','| --- | --- | --- | --- | --- | --- |')
for(x in issues)lines<-c(lines,paste0('| ',x$claim,' | ',escape(x$section),' | ',escape(paste0(x$value,' — ',x$quote)),
 ' | ',escape(if(all(is.na(x$source)))'NA (no registered value)' else format(x$source,digits=16)), ' | ',escape(x$why),' | ',escape(paste0(x$file,' :: ',x$field)),' |'))
writeLines(lines,'paper/audit-mismatches.md',useBytes=TRUE)
versions<-list(R=as.character(getRversion()),scp=as.character(utils::packageVersion('scp')),limma=as.character(utils::packageVersion('limma')),
 edgeR=as.character(utils::packageVersion('edgeR')),scplainer=paste0('workflow in scp ',utils::packageVersion('scp'),'; no separate package'))
reference_versions<-m3$versions
for(n in names(versions))if(!identical(versions[[n]],reference_versions[[n]]))warning('Version mismatch for ',n,': run=',versions[[n]],'; registered=',reference_versions[[n]],'; continuing.')
manifest<-list(manuscript_file=manuscript_path,manuscript_sha256=digest::digest(file=manuscript_path,algo='sha256'),
 input_hashes=inputs,seed=NULL,seed_reason='Deterministic audit of existing reports; no random draws',versions=versions,registered_versions=reference_versions,
 wall_seconds=proc.time()[['elapsed']]-start,claims_checked=nrow(out),mismatches=sum(!out$match),new_scientific_runs=FALSE,
 rounding='round(source_value, displayed decimal places); integer counts exact at zero decimals; quarter=0.25; halves=0.5. Unrounded target predicates checked separately.',
 outputs=setNames(lapply(c('paper/audit-table.csv','paper/audit-mismatches.md','paper/audit-coverage.csv'),function(p)digest::digest(file=p,algo='sha256')),
 c('paper/audit-table.csv','paper/audit-mismatches.md','paper/audit-coverage.csv')))
jsonlite::write_json(manifest,'paper/audit-manifest.json',pretty=TRUE,auto_unbox=TRUE,digits=NA)
cat('No calibration, simulation or validation was rerun; manuscript and registered inputs were read only.\n')
cat(sprintf('Claims checked: %d; mismatches: %d\n',nrow(out),sum(!out$match)))
