# Append interpretation only. Frozen M3 numeric tables and plots stay unchanged.
if(dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
source('R/scp_calibrate.R')
start<-proc.time()[['elapsed']];v<-.scp_versions();vt<-paste(names(v),v,sep='=',collapse='; ')
path<-'reports/M3b/power-report.md';original<-'reports/M3b/power-report-before-reading.md'
if(!file.exists(original))file.copy(path,original)
d<-read.csv('reports/M3b/power-grid.csv');d$fold_change<-2^d$beta
large<-d[d$n_per_condition==12 & d$m_cells_per_sample==64 & d$beta==1,]
large$power_over_testability<-large$power/large$testability
write.csv(large,'reports/M3b/reading-grid-diagnostics.csv',row.names=FALSE,na='NA')
one<-d[abs(d$beta-log2(1.25))<1e-12,]
maxima<-aggregate(power~dataset+k,one,max)
write.csv(maxima,'reports/M3b/reading-grid-small-effects.csv',row.names=FALSE,na='NA')
lines<-c('\n## Reading the grid',
 'The grid measures unconditional power: a changed protein counts in the denominator even when missingness prevents it from being tested. A large effect can make the test reliable for the proteins it can test while leaving total power below the target.',
 '\n### Testability caps twofold power',
 'At n=12 samples per condition and m=64 cells per sample, k=1, the following estimates come directly from the combined M3b table:',
 '\n| Preset | Unconditional power | Testability | Power / testability |','| --- | ---: | ---: | ---: |')
for(i in which(large$k==1)){r<-large[i,];lines<-c(lines,sprintf('| %s | %.3f | %.3f | %.2f |',r$dataset,r$power,r$testability,r$power_over_testability))}
lines<-c(lines,
 '\nFor carrier TMT, pooled-carrier p_g and q_g do not respond to the injected effect; proteins failing the three-cell minimum in any sample remain untestable. At twofold and k=1, missingness therefore sets the main limit. Power / overall testability is a diagnostic ratio, not the new conditional-power column: power uses changed proteins and testability uses all proteins. Random changed-protein selection explains a ratio slightly above one. The new conditional power divides detected changed proteins by testable changed proteins within each replicate and cannot exceed one.',
 '\n### Donor variance removes 1.25-fold changes from this grid',
 'No registered n, m or k combination reaches 80% mean power for a 1.25-fold change under the fixed Walzik donor variances. Increasing cells reduces cell-level noise and can improve testability; it does not remove the independent donor effect. More biological samples or smaller donor variance would define a different design or assumption and were not evaluated here. The grid cannot establish that whole-PBMC donor variance transfers to these cell lines; M4 injects the same external variance and cannot validate its transfer.',
 '\n### Read the k sensitivity alongside testability',
 'At the same n=12, m=64 and twofold effect:',
 '\n| Carrier-TMT preset | k | Power | Testability |','| --- | ---: | ---: | ---: |')
for(i in which(large$dataset%in%c('specht2019v3','leduc2022_pSCoPE'))){r<-large[i,];lines<-c(lines,sprintf('| %s | %.2f | %.3f | %.3f |',r$dataset,r$k,r$power,r$testability))}
lines<-c(lines,
 '\nCompression leaves TMT detection assumptions unchanged but reduces the tested condition effect. Testability stays nearly constant while power falls as k moves from 1 to 0.75 or 0.5. This is a registered sensitivity analysis; k was not estimated and M4 cannot validate it. Read low realised FDR alongside mean rejection counts and the no-rejection fraction: a near-zero FDP can describe very few discoveries.',
 paste0('\nSoftware: ',vt),
 '\nTMT numeric rows use amendment I; non-TMT rows are reused from v0.2.0. The archived original grid and report remain available in archive-v0.2.0/. reading-grid-manifest.json records these regenerated tables.')
writeLines(c(readLines(original),lines),path)
jsonlite::write_json(list(seed=20261013L,versions=as.list(v),input_hashes=list(power_grid=digest::digest(file='reports/M3b/power-grid.csv',algo='sha256'),original_report=digest::digest(file=original,algo='sha256'),script=digest::digest(file='inst/scripts/read_m3_grid.R',algo='sha256')),output_hashes=list(report=digest::digest(file=path,algo='sha256'),diagnostics=digest::digest(file='reports/M3b/reading-grid-diagnostics.csv',algo='sha256'),small_effects=digest::digest(file='reports/M3b/reading-grid-small-effects.csv',algo='sha256')),non_TMT_results_unchanged=TRUE,TMT_rows_updated_under_I=TRUE,wall_seconds=proc.time()[['elapsed']]-start),'reports/M3b/reading-grid-manifest.json',pretty=TRUE,auto_unbox=TRUE,digits=NA)
