# Regenerate all resumed plots and the one-page note without changing fits.
if (dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'),.libPaths()))
source('R/scp_calibrate.R')
started<-proc.time()[['elapsed']];out<-'reports/M1-resumed'
results<-readRDS(file.path(out,'calibration-results.rds'))
table<-read.csv(file.path(out,'calibration-table.csv'))
s<-read.csv(file.path(out,'dataset-summary.csv'))
manifest<-jsonlite::read_json(file.path(out,'run-manifest.json'),simplifyVector=FALSE)
bench<-jsonlite::read_json(file.path(out,'paper-benchmarks.json'),simplifyVector=TRUE)$benchmarks
versions<-.scp_versions();version_text<-paste(names(versions),versions,sep='=',collapse='; ')
manifest$versions<-as.list(versions);manifest$cohort$versions<-as.list(versions)
manifest$cohort$sorted_subset_proteomics<-FALSE
manifest$cohort$sorted_subset_evidence<-'Walzik methods: unsorted PBMC pellet used for proteomics; lymphocyte subsets measured by spectral flow cytometry, stored separately on Figshare; no sorted subset protein quantities in PXD058573'
manifest$cohort$composition_effect<-'Whole-PBMC donor variance includes differences in immune-cell composition; not a lineage-specific donor estimate.'
manifest$cohort$lineage_transfer<-'U-937 is monocyte-like. Melanoma, PDAC and HeLa are not lymphoid cells; no matching lineage proteomics exists here. No sorted or lymphoid estimate is fabricated.'
for(n in names(manifest$datasets)) manifest$datasets[[n]]$versions<-as.list(versions)
for(f in c('resumed-source-downloads','resumed-bulk-downloads','resumed-assay-doc-downloads','specht-analysis-source','bulk-metadata-download-manifest'))manifest[[f]]<-jsonlite::read_json(paste0('.local/sources/',f,'.json'))
manifest$engineering_checks<-list(local_tests='77 passed; no failures, warnings, or skips; all existing tests retained',test_seeds=c(304,508),local_check='Status: OK (CRAN incoming remote lookup disabled after a public CRAN archive timeout; all package checks and tests executed)')
manifest$missingness_equation<-'q(mu)=plogis(alpha + beta*(mu-center)); missingness=1-q; raw alpha=alpha-beta*center. Effective beta is zero when its Wald interval includes zero, preserving centered alpha.'
manifest$variance_equation<-'Log2 y = cell-line fixed mean + batch random mean + residual; REML on structurally mixed batches, without donor random effects. Aliased batches excluded from batch component; no-shared-batch cell variance needs replicated cells within line/batch.'
manifest$report_generation<-list(seed=20261013,versions=as.list(versions), input_hashes=setNames(lapply(c('calibration-results.rds','calibration-table.csv','dataset-summary.csv','paper-benchmarks.json'),function(f)digest::digest(file=file.path(out,f),algo='sha256')),c('results','table','summary','papers')))
labels<-names(results)
# Report the alias rule for each line/batch in the actual acquisition design.
for(n in labels) {
  input <- readRDS(paste0('.local/datasets/',n,'-aggregated.rds'))
  md <- input$cell_metadata
  lc <- manifest$datasets[[n]]$sample_col; bc <- manifest$datasets[[n]]$batch_col
  design <- do.call(rbind,lapply(split(seq_len(nrow(md)),md[[bc]]),function(j) {
    lines <- unique(as.character(md[[lc]][j]))
    data.frame(batch=as.character(md[[bc]][j[1]]),cell_line=lines,
      cells=vapply(lines,function(v)sum(md[[lc]][j]==v),integer(1)),
      sigma2_b_design=if(length(lines)>1)'Mixed-line batch; shared random effect can be identified across batches' else 'NA for this batch/line: only one line; aliased under amendment B')
  }))
  design$software_versions<-version_text
  write.csv(design,file.path(out,paste0(n,'-batch-design.csv')),row.names=FALSE)
}
plot_footer<-function() {mtext('k=1 is an assumption. Walzik donor variance is a whole-PBMC proxy. Flagged fits are diagnostic.',side=1,outer=TRUE,line=0.8,cex=.65);mtext(version_text,side=1,outer=TRUE,line=2.0,cex=.48)}
png(file.path(out,'variance-components.png'),width=2100,height=1250,res=180)
par(mfrow=c(2,3),mar=c(4,4,4,1),oma=c(4,0,2,0))
for(n in labels){d<-table[table$dataset==n,];v<-as.matrix(d[,c('sigma2_s','sigma2_b','sigma2_e')]);vmax<-max(v,na.rm=TRUE);if(!is.finite(vmax))vmax<-1
 matplot(d$tier,v,type='b',pch=c(16,17,15),lty=1,col=c('#d17c00','#7954a1','#1680a1'),xlab='Abundance quintile',ylab='Variance (log2 squared)',main=gsub('derks2022_','Derks ',n),ylim=c(0,vmax*1.08))
 if(grepl('STOPPED',d$status[1]))mtext('FLAG: diagnostic input/fit',side=3,line=.2,col='#b22222',cex=.7)
 if(n=='brunner2022')mtext('Cell and batch NA: one cell/run',side=3,line=.2,cex=.7)
}
plot.new();legend('center',legend=c('Between-donor: whole PBMC proxy','Batch (where identified)','Cell residual (where identified)'),col=c('#d17c00','#7954a1','#1680a1'),pch=c(16,17,15),lty=1,bty='n',cex=.9);plot_footer();dev.off()
png(file.path(out,'missingness-fit.png'),width=2100,height=1250,res=180)
par(mfrow=c(2,3),mar=c(4,4,4,1),oma=c(4,0,2,0))
for(n in labels){r<-results[[n]];p<-r$proteins;plot(p$mean_log2,p$missing_fraction,pch=16,cex=.25,col='#88888855',xlab='Mean observed log2 intensity',ylab='Missing fraction',main=gsub('derks2022_','Derks ',n),ylim=c(0,1));lines(r$curve$expected_log2,r$curve$missing_probability,lwd=2,col=if(length(r$manifest$flags))'#b22222' else '#1680a1');mtext(sprintf('Quantification slope %.3f [%.3f, %.3f]',r$missingness$slope,r$missingness$slope_low,r$missingness$slope_high),side=3,line=.2,cex=.65)}
plot.new();text(.5,.65,'Points: observed protein missingness\nLine: registered fitted curve\nRed: flagged diagnostic dataset',cex=.95);text(.5,.35,'TMT input ratios are not absolute abundance.\nNo imputation or model adjustment.',cex=.85);plot_footer();dev.off()
fmt<-function(x,digits=3)format(round(x,digits),nsmall=digits,trim=TRUE)
pct<-function(x)paste0(fmt(100*x,1),'%')
note<-c('# scpPower M1 plausibility note — frozen addendum v1.1.1',
 '13 October 2026. Every specified assay was evaluated; both Derks instruments are separate. Values below are unchanged estimates, not targets for tuning. Cultures are cell lines, not donors. Missingness uses observed pre-shift protein means; only variances receive per-cell median shifts.',
 '',
 sprintf('**Specht — FLAG / dataset stopped.** Missing %s, mean %s proteins/cell and median peptide-agreement cell CV %s. The paper implies about %s missing and %s proteins/cell from %s cells/protein; its %s peptide missingness is conditional on within-set identification, not global completeness. Coverage is consistent. Its reconstructed CV cannot reproduce the pre-normalisation QC metric exactly (FLAG: processed peptide input). However, the specified `peptides` is already log2 and row/column normalised (source export `cr_norm_log(t3)`); raw abundance cannot be recovered. Quantification slope %s [%s, %s] is implausibly negative. These tier/curve/variance estimates remain diagnostic. [Specht](%s).',pct(s$missing_fraction[1]),fmt(s$mean_proteins_per_cell[1],1),fmt(s$median_cell_peptide_cv[1]),pct(bench$specht$global_missing_fraction),fmt(bench$specht$implied_mean_proteins_per_cell,1),fmt(bench$specht$mean_cells_per_protein,0),pct(bench$specht$within_identified_peptide_missing/100),fmt(s$quantification_slope[1]),fmt(s$slope_low[1]),fmt(s$slope_high[1]),bench$specht$source),
 '',
 sprintf('**Leduc — %s.** Missing %s; mean %s proteins/cell versus the paper\'s %s (implied missing %s); median peptide-agreement cell CV %s versus the source QC cutoff %s (median below cutoff; source does not state an overall central cell CV). FLAG: coverage differs after explicit reverse/contaminant exclusion and reaggregation, with no filtering added to force agreement. `peptides`, not `peptides_log`, was used: the latter changes features, cells and normalisation. Input is reference-channel ratios. Quantification slope %s [%s, %s]; %s. [Leduc](%s).',if(grepl('STOPPED',s$status[2]))'FLAG / dataset stopped' else 'calibrated',pct(s$missing_fraction[2]),fmt(s$mean_proteins_per_cell[2],1),fmt(bench$leduc$mean_proteins_per_cell,0),pct(bench$leduc$implied_missing_fraction),fmt(s$median_cell_peptide_cv[2]),fmt(bench$leduc$cv_threshold,2),fmt(s$quantification_slope[2]),fmt(s$slope_low[2]),fmt(s$slope_high[2]),if(s$slope_high[2]<0)'FLAG: significantly negative; no corrective refit' else 'positive abundance dependence',bench$leduc$source),
 '',
 sprintf('**Derks — timsTOF / Q Exactive.** Global missing %s / %s; mean proteins/cell %s / %s versus approximately %s. Within-plex completeness %s / %s: timsTOF agrees with the published >%s%% claim; global and within-plex denominators are distinct. FLAG: Q Exactive coverage is below the headline count; instrument/gradient and early-assay scope differ. Median observed protein CV across cells %s / %s; the paper\'s repeat-injection CV is not cell-to-cell CV, so exact cell-CV agreement is untestable. File/label matching corrected permuted metadata, with verified instrument and single-cell identities. [Derks](%s).',pct(s$missing_fraction[3]),pct(s$missing_fraction[4]),fmt(s$mean_proteins_per_cell[3],1),fmt(s$mean_proteins_per_cell[4],1),fmt(bench$derks$proteins_per_cell,0),pct(s$mean_within_plex_completeness[3]),pct(s$mean_within_plex_completeness[4]),fmt(bench$derks$within_plex_completeness,0),fmt(s$median_observed_protein_cv[3]),fmt(s$median_observed_protein_cv[4]),bench$derks$source),
 '',
 sprintf('**Brunner — calibrated missingness, unavailable variance.** Raw protein missing %s; mean %s proteins/cell, minimum %s versus the paper\'s %s proteins/cell cutoff (FLAG: below it). Median observed linear protein CV %s versus reported %s: FLAG, different upstream processing/QC. No missingness percentage matching this raw assay is stated in the paper; completeness is abundance-dependent. Each run contains one HeLa cell, so both residual-within-run and separate run variance are unidentifiable; explained NA, not zero. [Brunner](%s).',pct(s$missing_fraction[5]),fmt(s$mean_proteins_per_cell[5],1),fmt(s$minimum_proteins_per_cell[5],0),fmt(bench$brunner$min_proteins_per_cell,0),fmt(s$median_observed_protein_cv[5]),fmt(bench$brunner$median_cv,1),bench$brunner$source),
 '',
 sprintf('**Walzik / compression.** %s healthy donors, %s baseline observations, exact identifier overlap (Specht, Leduc, timsTOF, Q Exactive, Brunner): %s. Proteomics is whole PBMC; subsets are flow cytometry, not sorted-subset proteomics. Lineage-specific donor variance is unavailable. Whole-PBMC variance includes donor composition differences and remains an unvalidated proxy for these cell lines. [Walzik](https://doi.org/10.1038/s41467-025-68101-9). k=1 is a registered assumption, not estimated; sensitivity assumptions are 0.75 and 0.5. The lab pilot needs a two-proteome ratio standard. The compression plot is dropped. [Cheung, Nature Methods](https://doi.org/10.1038/s41592-020-01002-5).',manifest$cohort$healthy_donors,manifest$cohort$baseline_observations,paste(s$bulk_overlap,collapse=', ')),
 '',paste('Reproduction: inst/scripts/calibrate_resumed.R and report_resumed.R; seed 20261013. Software:',version_text))
writeLines(note,file.path(out,'plausibility-note.md'))
plain<-gsub('\\[([^]]+)\\]\\([^)]*\\)','\\1',note);plain<-gsub('[`*#]','',plain);plain<-gsub('—','-',plain,fixed=TRUE)
wrapped<-unlist(lapply(plain,function(x)if(nzchar(x))strwrap(x,width=113) else ''))
pdf(file.path(out,'plausibility-note.pdf'),width=8.27,height=11.69,title='scpPower M1 resumed plausibility note')
par(mar=c(0,0,0,0));plot.new();plot.window(xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i')
line_height<-min(.0148,.9/length(wrapped));for(i in seq_along(wrapped))text(.055,.95-(i-1)*line_height,wrapped[i],adj=c(0,1),cex=if(i==1)1 else .68)
dev.off()
manifest$report_generation$wall_seconds<-proc.time()[['elapsed']]-started
manifest$decision<-list(registered_datasets_evaluated=TRUE,compression_plot_dropped=TRUE,stopped_datasets=labels[vapply(results,function(r)length(r$manifest$flags)>0,logical(1))],accepted_datasets=labels[vapply(results,function(r)length(r$manifest$flags)==0,logical(1))],all_datasets_accepted=all(vapply(results,function(r)length(r$manifest$flags)==0,logical(1))),M1_deliverables_generated=TRUE,no_model_adjustments_after_flags=TRUE)
jsonlite::write_json(manifest,file.path(out,'run-manifest.json'),pretty=TRUE,auto_unbox=TRUE,na='null',digits=NA)
writeLines(c('# Resumed M1 evidence','', 'Frozen addendum v1.1.1 governs this run. Earlier stopped evidence is retained in ../M1.', '', '[Calibration table](calibration-table.csv), [dataset summary](dataset-summary.csv), [variance plot](variance-components.png), [missingness plot](missingness-fit.png), [one-page note](plausibility-note.pdf), [manifest](run-manifest.json).', '', 'STOPPED dataset rows are diagnostic and not accepted calibration. Brunner run/cell components are unidentifiable because one cell was measured per run. Whole-PBMC donor variance is an unvalidated proxy; no sorted-subset proteomics is available. k is assumed, not estimated. All reported numbers are regenerated by the scripts listed in the root README.'),file.path(out,'README.md'))
cat('Plots, one-page note, and manifest generated.\n')
