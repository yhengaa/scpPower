# Manuscript artwork from immutable registered outputs; no fitting or simulation.
# Run from the repository root with the recorded renv environment restored.
if (dir.exists('.local/R-library')) .libPaths(c(normalizePath('.local/R-library'), .libPaths()))
suppressPackageStartupMessages({library(ggplot2);library(grid)})
started <- proc.time()[['elapsed']]
out <- 'paper/jpr/figures';dir.create(out,recursive=TRUE,showWarnings=FALSE)
inputs <- list()
allowed <- c('reports/M1-final/','reports/M2-revised/','reports/M2b/',
             'reports/M3b/','reports/M4/','reports/M4b/','reports/J/')
read_input <- function(path,json=FALSE) {
    stopifnot(any(vapply(allowed,function(prefix)startsWith(path,prefix),logical(1))))
    inputs[[path]] <<- digest::digest(file=path,algo='sha256')
    if(json) jsonlite::read_json(path,simplifyVector=FALSE) else read.csv(path,check.names=FALSE)
}
write_data <- function(x,name) {
    x$render_software_versions<-paste(names(versions),versions,sep='=',collapse='; ')
    write.csv(x,file.path(out,name),row.names=FALSE,na='NA')
}
ids <- c('specht2019v3','leduc2022_pSCoPE','derks2022_timsTOF','derks2022_QExactive','brunner2022')
labels <- setNames(c('Specht 2019 TMT','Leduc 2022 TMT','Derks 2022 timsTOF plexDIA',
                     'Derks 2022 Q Exactive plexDIA','Brunner 2022 label-free'),ids)
# Okabe-Ito: five fixed preset colours; its remaining colours encode fold change.
palette <- setNames(c('#E69F00','#0072B2','#009E73','#CC79A7','#D55E00'),ids)
fold_palette <- setNames(c('#000000','#56B4E9','#F0E442'),c('1.25×','1.5×','2×'))
wrap_labels <- function(x) vapply(labels[as.character(x)],function(z)paste(strwrap(z,width=18),collapse='\n'),character(1))
tag_series <- function(x) {x$dataset<-factor(x$dataset,levels=ids);x}
style <- theme_bw(base_size=8,base_family='sans') + theme(
    text=element_text(size=8),axis.text=element_text(size=7),
    axis.title=element_text(size=8),strip.text=element_text(size=7),
    strip.background=element_blank(),panel.grid.minor=element_blank(),
    panel.grid.major=element_line(linewidth=.15,colour='grey92'),
    legend.text=element_text(size=7),legend.title=element_text(size=7),
    legend.position='bottom',legend.box='vertical',legend.key.width=unit(14,'pt'),
    plot.margin=margin(14,7,7,7,unit='pt'))
series_scale <- function(legend=TRUE) scale_colour_manual(values=palette,
    breaks=ids,labels=function(x)vapply(labels[x],function(z)paste(strwrap(z,27),collapse='\n'),character(1)),
    name=NULL,guide=if(legend)guide_legend(nrow=2,byrow=TRUE) else 'none',drop=FALSE)
facets <- function(ncol=5,free='fixed') facet_wrap(~dataset,ncol=ncol,scales=free,
    labeller=labeller(dataset=wrap_labels))
# All text, including custom diagrams, is >=7 pt at the exported physical size.
placements <- function(plot,row,col=1,colspan=1,rowspan=1,letter=NULL)
    list(plot=plot,row=row,col=col,colspan=colspan,rowspan=rowspan,letter=letter)
render <- function(name,width,height,draw) {
    grDevices::cairo_pdf(file.path(out,paste0(name,'.pdf')),width=width,height=height,
                        family='sans',pointsize=8,onefile=TRUE)
    draw();dev.off()
    grDevices::png(file.path(out,paste0(name,'.png')),width=width,height=height,
                  units='in',res=600,type='cairo',pointsize=8)
    draw();dev.off()
}
compose <- function(name,width,height,items,heights,widths=c(1,1)) {
    draw <- function() {
        grid.newpage();pushViewport(viewport(layout=grid.layout(length(heights),length(widths),
            heights=unit(heights,'null'),widths=unit(widths,'null'))))
        for(item in items) {
            pushViewport(viewport(layout.pos.row=item$row+seq_len(item$rowspan)-1,
                                  layout.pos.col=item$col+seq_len(item$colspan)-1))
            if(is.function(item$plot))item$plot() else grid.draw(ggplotGrob(item$plot))
            if(!is.null(item$letter))grid.text(item$letter,x=unit(3,'pt'),y=unit(1,'npc')-unit(2,'pt'),
                just=c('left','top'),gp=gpar(fontsize=10,fontface='bold'))
            popViewport()
        }
        popViewport()
    }
    render(name,width,height,draw)
}
cal <- read_input('reports/M1-final/calibration-table.csv')
m1 <- read_input('reports/M1-final/run-manifest.json',TRUE)
versions <- c(R=as.character(getRversion()),scp=as.character(utils::packageVersion('scp')),
    limma=as.character(utils::packageVersion('limma')),edgeR=as.character(utils::packageVersion('edgeR')),
    scplainer=paste0('workflow in scp ',utils::packageVersion('scp'),'; no separate package'))
ref <- unlist(m1$versions);common<-intersect(names(ref),names(versions))
if(any(ref[common]!=versions[common]))warning('Software version mismatch; continuing figure rendering.')
proteins <- setNames(lapply(ids,function(n)read_input(paste0('reports/M1-final/',n,'-proteins.csv'))),ids)
pp_list <- lapply(ids,function(n){
    path<-paste0(if(n %in% ids[1:2])'reports/M2b/' else 'reports/M2-revised/',n,'-protein-parameters.csv')
    x<-read_input(path);x$dataset<-n;x
})
pp_fields<-unique(unlist(lapply(pp_list,names)))
pp<-do.call(rbind,lapply(pp_list,function(x){for(f in setdiff(pp_fields,names(x)))x[[f]]<-NA;x[,pp_fields,drop=FALSE]}))
# Figure 1: symbolic model and API workflow; grids are explicitly schematic.
model_data <- data.frame(term=c('mu_g','k * beta_g * C_s','u_gs','b_gb','epsilon_gsbc'),
    meaning=c('protein mean','condition effect with assumed compression','biological sample effect',
              'acquisition batch effect','cell residual'),
    distribution=c('fixed','fixed','Normal(0, sigma2_s,g)','Normal(0, sigma2_b,g)','Normal(0, sigma2_e,g)'))
mechanisms <- data.frame(mechanism=c('carrier TMT','DIA / label-free'),
    identification=c('I_gb ~ Bernoulli(p_g)','per-cell quantification'),
    quantification=c('Q_gbc | I_gb=1 ~ Bernoulli(q_g)',
        'logit P(M=1) = logit(r_g) + slope * (expected_log2 - mean_log2_g)'),
    observation=c('I_gb * Q_gbc','M_gc'))
write_data(model_data,'Fig1a-data.csv');write_data(mechanisms,'Fig1a-missingness-data.csv')
schematic<-expand.grid(mechanism=c('carrier TMT','DIA / label-free'),plex=1:2,channel=1:4)
schematic$quantified<-with(schematic,ifelse(mechanism=='carrier TMT',plex==1 & channel!=3,(plex+channel)%%3!=0))
schematic$illustrative_only<-TRUE;write_data(schematic,'Fig1a-grid-data.csv')
flow<-data.frame(order=seq_len(7),node=c('Processed dataset','scp_calibrate()','Preset',
    'scp_simulate()','Simulated experiment','scp_power()','Power table + hours'))
write_data(flow,'Fig1b-data.csv')
textg <- function(text,x,y,size=8,...)grid.text(text,x=x,y=y,gp=gpar(fontsize=size,...))
lineg <- function(x,y,arrow=FALSE,dashed=FALSE)grid.lines(x=x,y=y,
    arrow=if(arrow)grid::arrow(length=unit(4,'pt'),type='closed') else NULL,
    gp=gpar(col='grey30',lwd=.7,lty=if(dashed)2 else 1))
model_draw <- function() {
    textg(expression(Y[g*","*s*","*b*","*c]==mu[g]+k*beta[g]*C[s]+u[g*","*s]+b[g*","*b]+epsilon[g*","*s*","*b*","*c]),.5,.89,12)
    textg('Protein mean     Condition effect     Sample effect     Batch effect     Cell noise',.5,.76,8)
    textg(expression(u[g*","*s]%~%N(0,sigma[s*","*g]^2)~';'~b[g*","*b]%~%N(0,sigma[b*","*g]^2)~';'~epsilon[g*","*s*","*b*","*c]%~%N(0,sigma[e*","*g]^2)),.5,.65,10)
    textg('Carrier TMT: identify plex, then quantify channels',.25,.52,8)
    textg('DIA / label-free: quantify each cell',.75,.52,8)
    for(side in c(.25,.75))for(row in 1:2)for(channel in 1:4) {
        quantified<-if(side==.25)row==1 && channel!=3 else (row+channel)%%3!=0
        grid.rect(x=side+(channel-2.5)*.037,y=.38-(row-1)*.10,
            width=.028,height=.068,gp=gpar(fill=if(quantified)'grey35' else 'white',col='grey35',lwd=.6))
    }
    textg(expression(P(I[g*","*b]==1)==p[g]),.25,.18,10)
    textg(expression(P(Q[g*","*b*","*c]==1~'|'~I[g*","*b]==1)==q[g]),.25,.07,10)
    textg(expression(logit(P(M[g*","*c]==1))==logit(r[g])+lambda*(x-mu[g])),.75,.12,10)
}
workflow_draw <- function() {
    xx<-seq(.065,.935,length.out=7)
    nodes<-c('Dataset','scp_calibrate()','Preset','scp_simulate()','Simulated\nexperiment','scp_power()','Power table\n+ hours')
    for(i in seq_along(xx)) {
        grid.roundrect(x=xx[i],y=.57,width=.127,height=.28,r=unit(3,'pt'),
            gp=gpar(fill=if(i %in% c(2,4,6))'grey94' else 'white',col='grey45',lwd=.6))
        textg(nodes[i],xx[i],.57,7.5)
        if(i<7)lineg(c(xx[i]+.065,xx[i+1]-.065),c(.57,.57),TRUE,dashed=i==5)
    }
    lineg(c(xx[3],xx[3],xx[6],xx[6]),c(.40,.22,.22,.40),TRUE)
    textg('Preset + design feed scp_power(); it simulates and tests internally',.57,.10,7.5)
    textg('SCP measurements + bulk donor variance',.16,.88,8)
}
compose('Fig1',7,4.1,list(placements(model_draw,1,colspan=2,letter='a'),
    placements(workflow_draw,2,colspan=2,letter='b')),heights=c(2.7,1.4))
# Figure 2a: square roots of the registered tier estimates; no variance re-fit.
sds<-do.call(rbind,lapply(c('sigma2_e','sigma2_b','sigma2_s'),function(component){
    x<-cal[,c('dataset','tier',component)];names(x)[3]<-'variance';x$component<-component
    x$sd_log2<-sqrt(x$variance);x
}))
sds<-tag_series(sds);write_data(sds,'Fig2a-data.csv')
component_labels<-c(sigma2_e='Cell',sigma2_b='Batch',sigma2_s='Mapped donor')
f2a<-ggplot(sds,aes(tier,sd_log2,colour=dataset,linetype=component,shape=component))+
    geom_line(linewidth=.4)+geom_point(size=1.4)+facets()+series_scale(FALSE)+
    scale_linetype_manual(values=c(sigma2_e=1,sigma2_b=2,sigma2_s=3),labels=component_labels,name=NULL)+
    scale_shape_manual(values=c(sigma2_e=16,sigma2_b=17,sigma2_s=15),labels=component_labels,name=NULL)+
    scale_x_continuous(breaks=1:5)+labs(x='Abundance tier',y='SD (log2 units)')+style
# Figure 2b: amendment-I fidelity for TMT, amendment-G fidelity for non-TMT.
fid_g<-read_input('reports/M2-revised/fidelity-table.csv')
fid_i<-read_input('reports/M2b/fidelity-table.csv')
fid<-rbind(fid_i[fid_i$metric=='missing_fraction',],
           fid_g[fid_g$metric=='missing_fraction' & fid_g$dataset %in% ids[3:5],])
fid$source_report<-ifelse(fid$dataset %in% ids[1:2],'M2b','M2-revised')
fid$band_low<-fid$observed-fid$tolerance;fid$band_high<-fid$observed+fid$tolerance
fid<-tag_series(fid);write_data(fid,'Fig2b-data.csv')
f2b<-ggplot(fid,aes(tier,colour=dataset))+
    geom_ribbon(aes(ymin=band_low,ymax=band_high),fill='grey90',colour=NA)+
    geom_line(aes(y=observed,linetype='Observed'),linewidth=.4)+
    geom_line(aes(y=simulated_mean,linetype='Simulated mean'),linewidth=.4)+
    geom_point(aes(y=observed),shape=1,size=1.3)+geom_point(aes(y=simulated_mean),shape=16,size=1)+
    facets()+series_scale(FALSE)+scale_linetype_manual(values=c(Observed=1,'Simulated mean'=2),name=NULL)+
    scale_x_continuous(breaks=1:5)+coord_cartesian(ylim=c(0,1))+
    labs(x='Abundance tier',y='Missing fraction')+style
# Figure 2c: deterministic 2-D bin counts, plus every source protein.
tmt<-pp[pp$dataset %in% ids[1:2],c('dataset','protein','tier','p_g','q_g')]
write_data(tmt,'Fig2c-protein-data.csv')
stopifnot(all(is.finite(tmt$p_g)),all(is.finite(tmt$q_g)))
tmt$xbin<-pmin(39L,floor(tmt$p_g*40));tmt$ybin<-pmin(39L,floor(tmt$q_g*40))
bins<-aggregate(list(proteins=tmt$protein),tmt[,c('dataset','xbin','ybin')],length)
bins$p_g_min<-bins$xbin/40;bins$p_g_max<-(bins$xbin+1)/40
bins$q_g_min<-bins$ybin/40;bins$q_g_max<-(bins$ybin+1)/40
bins$p_g_centre<-(bins$xbin+.5)/40;bins$q_g_centre<-(bins$ybin+.5)/40
bins<-tag_series(bins);write_data(bins,'Fig2c-data.csv')
f2c<-ggplot(bins,aes(p_g_centre,q_g_centre,fill=dataset,alpha=proteins))+
    geom_tile(width=1/40,height=1/40)+facets(2)+
    scale_fill_manual(values=palette,guide='none')+
    scale_alpha_continuous(transform='log10',range=c(.2,1),name='Proteins / bin')+
    coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+
    scale_x_continuous(breaks=c(0,.5,1),labels=c('0','0.5','1'))+
    labs(x='Plex identification p_g',y='Channel quantification q_g')+style+
    theme(panel.spacing.x=unit(10,'pt'))
# Figure 2d: r_g and mean intensity from registered G priors, with the M1 curve.
dia<-pp[pp$dataset %in% ids[3:5],c('dataset','protein','tier','mean_log2','r_g')]
dia<-tag_series(dia);write_data(dia,'Fig2d-data.csv')
curves<-do.call(rbind,lapply(ids[3:5],function(n){
    x<-read_input(paste0('reports/M1-final/',n,'-missingness-curve.csv'))
    x$dataset<-n;x$quantification_probability<-1-x$missing_probability;x
}))
curves<-tag_series(curves);write_data(curves,'Fig2d-curve-data.csv')
centres<-cal[cal$dataset %in% ids[3:5]&cal$tier==1,
    c('dataset','missingness_intercept','effective_slope','intensity_center')]
centres<-tag_series(centres);write_data(centres,'Fig2d-parameters.csv')
f2d<-ggplot(dia,aes(mean_log2,r_g,colour=dataset))+
    geom_point(alpha=.22,size=.4)+
    geom_line(data=curves,aes(expected_log2,quantification_probability),linewidth=.65)+
    geom_vline(data=centres,aes(xintercept=intensity_center),colour='grey55',linetype=3,linewidth=.3)+
    facets(3,'free_x')+series_scale(FALSE)+coord_cartesian(ylim=c(0,1))+
    scale_x_continuous(breaks=function(limits)pretty(limits,n=2))+
    labs(x='Mean log2 intensity',y='Quantification r_g')+style
compose('Fig2',7,8.2,list(placements(f2a,1,colspan=2,letter='a'),
    placements(f2b,2,colspan=2,letter='b'),placements(f2c,3,letter='c'),
    placements(f2d,3,col=2,letter='d')),heights=c(2.65,2.65,2.9))
# Figure 3 uses power-grid.csv only: all filters and time values come from it.
grid_power<-tag_series(read_input('reports/M3b/power-grid.csv'))
levels_n<-sort(unique(grid_power$n_per_condition))
cell_breaks<-sort(unique(grid_power$m_cells_per_sample))
curve_plot <- function(x) {
    x$n_label<-factor(x$n_per_condition,levels=levels_n)
    ggplot(x,aes(m_cells_per_sample,power,group=n_label,linetype=n_label,colour=dataset))+
        geom_hline(yintercept=.8,linetype=3,colour='grey45',linewidth=.3)+
        geom_line(linewidth=.5)+geom_point(size=.9)+facets()+series_scale(FALSE)+
        scale_linetype_manual(values=c('solid','dashed','dotdash','longdash','twodash'),name='Samples / condition')+
        scale_x_continuous(transform='log2',breaks=cell_breaks)+
        coord_cartesian(ylim=c(0,1))+labs(x='Cells / sample',y='Unconditional power')+style
}
f3a_data<-grid_power[grid_power$beta==log2(2)&grid_power$k==1,]
f3b_data<-grid_power[abs(grid_power$beta-log2(1.5))<1e-12&grid_power$k==1,]
write_data(f3a_data,'Fig3a-data.csv');write_data(f3b_data,'Fig3b-data.csv')
f3a<-curve_plot(f3a_data);f3b<-curve_plot(f3b_data)
# Match the 2x, k=1 rows; detection can respond to effect size in DIA/label-free.
f3c_data<-f3a_data[f3a_data$n_per_condition==12,]
write_data(f3c_data,'Fig3c-data.csv')
f3c<-ggplot(f3c_data,aes(m_cells_per_sample,testability,colour=dataset))+
    geom_line(linewidth=.55)+geom_point(size=1.1)+series_scale(TRUE)+
    scale_x_continuous(transform='log2',breaks=cell_breaks)+coord_cartesian(ylim=c(0,1))+
    labs(x='Cells / sample',y='Testability')+style+
    guides(colour=guide_legend(ncol=1,byrow=TRUE))+
    theme(legend.key.height=unit(7,'pt'),legend.margin=margin(0,0,0,0))
f3d_data<-grid_power[grid_power$dataset %in% ids[1:2]&grid_power$beta==1&grid_power$n_per_condition==12,]
f3d_data$k_label<-factor(f3d_data$k,levels=c(1,.75,.5))
write_data(f3d_data,'Fig3d-data.csv')
f3d<-ggplot(f3d_data,aes(m_cells_per_sample,power,colour=dataset,linetype=k_label,group=k_label))+
    geom_line(linewidth=.6)+geom_point(size=1.1)+facets(2)+series_scale(FALSE)+
    scale_linetype_manual(values=c('solid','dashed','dotdash'),name='Assumed k')+
    scale_x_continuous(transform='log2',breaks=cell_breaks)+coord_cartesian(ylim=c(0,1))+
    labs(x='Cells / sample',y='Unconditional power')+style
f3e_data<-grid_power[grid_power$beta==1,]
passing<-f3e_data[f3e_data$power>=.8,]
stopifnot(nrow(passing)>0)
winning<-passing[passing$instrument_hours==min(passing$instrument_hours),]
stopifnot(nrow(winning)==1)
f3e_data$cheapest_reaching_target<-f3e_data$task==winning$task
write_data(f3e_data,'Fig3e-data.csv');write_data(winning,'Fig3e-cheapest-point.csv')
time_plot<-function(legend=TRUE,compact=FALSE) {
    x<-ggplot(f3e_data,aes(instrument_hours,power,colour=dataset))+
        geom_hline(yintercept=.8,colour='grey45',linetype=3,linewidth=.35)+
        geom_point(alpha=.55,size=if(compact).95 else 1.2)+
        geom_point(data=winning,shape=21,fill=NA,colour='black',size=2.7,stroke=.55)+
        series_scale(legend)+scale_x_log10()+coord_cartesian(ylim=c(0,1))+
        labs(x='Instrument hours',y='Unconditional power')+style
    if(!compact)x<-x+annotate('text',x=winning$instrument_hours,y=.94,
        label=sprintf('n = %d; m = %d; %.0f h',winning$n_per_condition,winning$m_cells_per_sample,winning$instrument_hours),
        size=7.5/ggplot2::.pt,colour='black')
    x
}
f3e<-time_plot(TRUE)
compose('Fig3',7,9,list(placements(f3a,1,colspan=2,letter='a'),
    placements(f3b,2,colspan=2,letter='b'),placements(f3c,3,letter='c'),
    placements(f3d,3,col=2,letter='d'),placements(f3e,4,colspan=2,letter='e')),
    heights=c(2.1,2.1,2.6,2.2))
# Figure 4: counted real-cell validation points; exclude construction-zero points.
m4<-read_input('reports/M4/validation-table.csv');m4b<-read_input('reports/M4b/validation-table.csv')
source_sets<-list(a=m4b[m4b$dataset==ids[1]&m4b$counted,],
    b=m4b[m4b$dataset==ids[2]&m4b$counted,],
    c=m4[m4$dataset==ids[5]&m4$counted,],
    d=m4[m4$dataset==ids[4]&m4$counted,],
    e=m4[m4$dataset==ids[1]&m4$counted,])
extent<-ceiling(max(unlist(lapply(source_sets[c('a','b','c','e')],function(x)
    c(x$predicted_power,x$observed_power))))*10)/10
fold_label<-function(beta)factor(c('1.25×','1.5×','2×')[
    match(round(beta,12),round(log2(c(1.25,1.5,2)),12))],levels=c('1.25×','1.5×','2×'))
f4plots<-list()
for(letter in names(source_sets)) {
    x<-source_sets[[letter]];stopifnot(nrow(x)>0,all(x$counted),all(x$m_cells_per_sample>=3))
    x$source_report<-if(letter %in% c('a','b'))'M4b (v0.3.1)' else 'M4 (v0.3.0)'
    x$fold_change<-fold_label(x$beta);x$identity_low<-x$observed_power-.1;x$identity_high<-x$observed_power+.1
    write_data(x,paste0('Fig4',letter,'-data.csv'))
    limit<-if(letter=='d').05 else extent
    if(letter=='d')stopifnot(all(x$observed_power<.02),all(x$predicted_power<.02))
    band<-data.frame(observed_power=c(0,limit),low=c(-.1,limit-.1),high=c(.1,limit+.1))
    write_data(band,paste0('Fig4',letter,'-band-data.csv'))
    p<-ggplot(x,aes(observed_power,predicted_power))+
        geom_ribbon(data=band,aes(observed_power,ymin=low,ymax=high),inherit.aes=FALSE,
            fill='grey92',colour=NA)+geom_abline(slope=1,intercept=0,colour='grey45',linewidth=.4)+
        geom_point(aes(fill=fold_change),shape=21,colour='grey25',stroke=.2,size=1.5)+
        scale_fill_manual(values=fold_palette,name='Fold change',drop=FALSE)+
        coord_fixed(xlim=c(0,limit),ylim=c(0,limit),expand=FALSE)+
        labs(x='Observed power',y='Predicted power')+style+
        theme(legend.position='none',plot.margin=margin(24,8,7,8,unit='pt'))
    panel_name<-paste0(labels[unique(x$dataset)],if(letter=='e')'\nbefore amendment I' else '')
    f4plots[[letter]]<-local({plot<-p;name<-panel_name;dataset<-unique(x$dataset)
        function(){grid.draw(ggplotGrob(plot));grid.text(name,x=.53,y=unit(1,'npc')-unit(4,'pt'),
            just=c('centre','top'),gp=gpar(fontsize=8,col=palette[[dataset]]))}})
}
f4legend<-function(){
    textg('Fold change',.5,.68,8)
    for(i in seq_along(fold_palette)) {
        xx<-c(.19,.5,.81)[i]
        grid.circle(x=xx,y=.54,r=unit(2.2,'pt'),gp=gpar(fill=fold_palette[i],col='grey25',lwd=.5))
        textg(names(fold_palette)[i],xx,.42,8)
    }
    textg('Shading: ±10 percentage points\nIdentity: predicted = observed',.5,.26,8)
}
compose('Fig4',7,7.7,list(placements(f4plots$a,1,letter='a'),placements(f4plots$b,1,2,letter='b'),
    placements(f4plots$c,2,letter='c'),placements(f4plots$d,2,2,letter='d'),
    placements(f4plots$e,3,letter='e'),placements(f4legend,3,2)),heights=c(2.6,2.6,2.5))
# TOC: exactly the Figure 3e values, reduced to the requested physical dimensions.
write_data(f3e_data,'TOC-data.csv')
toc<-time_plot(FALSE,TRUE)+labs(y='Power')+
    theme(plot.margin=margin(2,5,2,2,unit='pt'),axis.title=element_text(size=7),axis.text=element_text(size=7))
render('TOC',3.25,1.75,function(){
    grid.newpage();pushViewport(viewport(layout=grid.layout(2,1,heights=unit(c(.34,1.41),'in'))))
    pushViewport(viewport(layout.pos.row=1));textg('How many samples and cells does a\nsingle-cell proteomics study need?',.5,.56,8);popViewport()
    pushViewport(viewport(layout.pos.row=2));grid.draw(ggplotGrob(toc));popViewport();popViewport()
})
# Tables: summaries of existing per-protein estimates and the registered grid.
info<-read_input('reports/M3b/instrument-preset-manifest.json',TRUE)
hours<-read_input('reports/M3b/instrument-hours.csv')
cheapest<-read_input('reports/M3b/cheapest-designs.csv')
table1<-do.call(rbind,lapply(ids,function(n){
    p<-proteins[[n]];d<-cal[cal$dataset==n,];meta<-m1$datasets[[n]]
    instrument<-unlist(meta$instrument)
    if(!length(instrument))instrument<-'not recorded in registered outputs'
    absorbed<-n==ids[5]
    data.frame(dataset=labels[[n]],preset=n,instrument=paste(instrument,collapse='; '),
        labelling=if(n %in% ids[1:2])'carrier TMT' else if(n %in% ids[3:4])'plexDIA' else 'label-free',
        source_assay=unique(d$assay),cells=unique(d$cells),cell_lines=unique(d$samples),
        proteins_per_cell=sum(p$quantified_cells,na.rm=TRUE)/unique(d$cells),
        missing_fraction=mean(p$missing_fraction,na.rm=TRUE),
        median_cell_sd=median(sqrt(p$sigma2_e),na.rm=TRUE),
        median_batch_sd=if(absorbed)'absorbed' else as.character(median(sqrt(p$sigma2_b),na.rm=TRUE)),
        mapped_donor_sd=median(sqrt(p$sigma2_s),na.rm=TRUE),
        mapped_donor_proteins=sum(is.finite(p$sigma2_s)),
        cells_per_run=info$presets[[n]]$cells_per_run,minutes_per_run=info$presets[[n]]$minutes_per_run,
        citation=info$presets[[n]]$source$source,
        scale='Endpoint protein matrix; mapped donor SD excludes tier fallbacks')
}))
# Require metadata to agree with the unchanged instrument-hour table.
for(n in ids){h<-hours[hours$dataset==n,];stopifnot(all(h$cells_per_run==info$presets[[n]]$cells_per_run),
    all(h$minutes_per_run==info$presets[[n]]$minutes_per_run))}
write_data(table1,'Table1.csv')
table2<-do.call(rbind,lapply(ids,function(n)do.call(rbind,lapply(sort(unique(grid_power$beta)),function(beta)
    do.call(rbind,lapply(if(n %in% ids[1:2])c(1,.75,.5) else 1,function(k){
        x<-grid_power[grid_power$dataset==n & abs(grid_power$beta-beta)<1e-12 & grid_power$k==k,]
        best<-x[order(-x$power,x$instrument_hours,x$cells,x$n_per_condition),][1,]
        c<-cheapest[cheapest$dataset==n & abs(cheapest$beta-beta)<1e-12 & cheapest$k==k,]
        stopifnot(nrow(c)==1)
        data.frame(dataset=labels[[n]],preset=n,fold_change=as.character(fold_label(beta)),k=k,
            maximum_power=best$power,n_at_maximum=best$n_per_condition,m_at_maximum=best$m_cells_per_sample,
            hours_at_maximum=best$instrument_hours,maximum_power_se=best$simulation_se,
            cheapest_status=if(is.na(c$n_per_condition))'none in grid' else 'reaches 0.8',
            cheapest_n=c$n_per_condition,cheapest_m=c$m_cells_per_sample,
            cheapest_power=c$power,cheapest_hours=c$instrument_hours)
    }))))))
write_data(table2,'Table2.csv')
markdown <- function(x,path) {
    z<-x
    for(n in names(z))if(is.numeric(z[[n]]))z[[n]]<-ifelse(is.na(z[[n]]),'—',format(round(z[[n]],3),trim=TRUE,nsmall=0))
    for(n in names(z))if(is.character(x[[n]])) {
        numeric_values<-suppressWarnings(as.numeric(x[[n]]));available<-is.finite(numeric_values)
        z[[n]][available]<-format(round(numeric_values[available],3),trim=TRUE,nsmall=3)
    }
    z[is.na(z)]<-'—'
    writeLines(c(paste0('| ',paste(names(z),collapse=' | '),' |'),
        paste0('| ',paste(rep('---',ncol(z)),collapse=' | '),' |'),
        apply(z,1,function(r)paste0('| ',paste(r,collapse=' | '),' |'))),file.path(out,path))
}
markdown(table1,'Table1.md');markdown(table2,'Table2.md')
# Provenance and captions: no new interpretation of the registered estimates.
j<-read_input('reports/J/decision-table.csv')
stopifnot(all(j$amendment_I_applies))
source_summary<-read_input('reports/M4/source-summary.csv')
tims_cells<-source_summary$cells[source_summary$dataset==ids[3]]
write_data(data.frame(dataset=j$dataset,aggregated_sparser_gap=j$aggregated_sparser_gap,
    amendment_I_applies=j$amendment_I_applies),'TMT-mask-provenance.csv')
captions<-c('# Caption drafts','',
 'Figure 1. Registered model and workflow. (a) The five terms are the protein log2 mean, the condition effect scaled by the assumed carrier compression k, the biological-sample effect, the acquisition-batch effect and the cell residual. The sample variance is transferred from whole-PBMC bulk data; k is supplied, not estimated. The schematic plex grids illustrate carrier-TMT identification per plex followed by quantification per channel, and DIA or label-free quantification per cell with the registered intensity slope. Filled and empty squares are illustrative, not measured counts. (b) Processed SCP and bulk donor data feed scp_calibrate(); the resulting preset feeds scp_simulate() and scp_power(). The latter accepts a preset and design, simulates internally, and applies the registered sample-mean limma-trend endpoint to produce a power table with nominal instrument hours. Model components and endpoint definitions are recorded in reports/M1-final and reports/M3b; the dotted experiment-to-power arrow denotes internal testing, not a public matrix argument.','',
 'Figure 2. Calibration and fidelity across the five series. (a) Cell, batch and mapped donor SDs are square roots of the registered per-tier variance estimates in reports/M1-final/calibration-table.csv. Brunner batch variance is absorbed into its cell/run component; its zero batch SD is not zero acquisition noise. (b) Observed missing fraction and simulated mean are shown by tier with the registered ±5-percentage-point tolerance band: the two TMT series use reports/M2b/fidelity-table.csv and the other series use reports/M2-revised/fidelity-table.csv. (c) Two-dimensional protein counts show plex identification p_g against within-identified-plex channel quantification q_g from the amendment-I endpoint masks in reports/M2b; bins have a common size and logarithmic count opacity. (d) Per-protein r_g against mean log2 intensity uses reports/M2-revised, with the unchanged registered logistic quantification curves and centre from reports/M1-final. The smooth lines are the recorded population fits, not new fits to r_g; the simulator anchors each protein at its own r_g. Colours identify presets throughout these panels.','',
 'Figure 3. Registered design grid from reports/M3b/power-grid.csv. (a, b) Unconditional power against cells per sample is shown at 2× and 1.5×, respectively, with k = 1, one facet per preset and line types for biological samples per condition. The cells axis is logarithmic in base two and the horizontal target is 0.8. (c) Testability is shown for 12 samples per condition at 2× and k = 1. (d) Carrier-TMT power is shown at 2× and 12 samples per condition for the assumed k values 1, 0.75 and 0.5. (e) Power is plotted vertically against nominal instrument hours on a logarithmic horizontal axis for every 2× grid point, including all registered TMT k values. The circled point is the sole minimum-hour design reaching 0.8; its n, m and hours are labelled from the grid. TMT rows are the M3b rerun and non-TMT rows are the unchanged registered M3 results. Each point summarizes the registered replicate run; no simulations were added for these figures.','',
 sprintf('Figure 4. Validation on real cells at counted grid points only. Predicted unconditional power is vertical and observed power horizontal; the line is identity and shading spans ±10 percentage points. Point fill denotes 1.25×, 1.5× or 2×, with preset identity shown by the panel label. (a, b) Specht and Leduc use reports/M4b/validation-table.csv after amendment I. (c, d) Brunner and Q Exactive use reports/M4/validation-table.csv. Q Exactive axes are limited to 0–0.05 and all its plotted predicted and observed powers are below 2%%. (e) Specht uses the accepted v0.3.0 M4 results before amendment I, on the same axes as (a). Points with fewer than three cells per sample are excluded under H2. timsTOF is not evaluable because it has %d source cells per line. These comparisons validate cell-level noise, missingness and the test within the public-cell power range; they cannot validate donor-variance transfer or k.',tims_cells),'',
 'Table 1. Calibrated series and acquisition metadata. Cells, cell-line counts and source assay come from reports/M1-final/calibration-table.csv. Mean quantified proteins per cell and missing fraction describe the endpoint protein matrix, calculated from the per-protein quantified-cell counts and missing fractions in reports/M1-final/*-proteins.csv; the more permissive raw-PSM TMT mask is not substituted. SDs are medians of square roots of measured per-protein variances, retaining boundary estimates and excluding unavailable estimates. Mapped donor SD includes only uniquely mapped measured donor variances, not tier fallbacks. Brunner batch noise is recorded as absorbed into the cell component. Cells per run, method minutes and citation reproduce the preset metadata recorded in reports/M3b/instrument-preset-manifest.json and instrument-hours.csv. Instrument metadata identifiers are retained; absent instrument model names are marked unrecorded rather than inferred. Markdown display rounds measured SDs and coverage to three decimals; CSV retains full precision.','',
 'Table 2. Registered design summaries from reports/M3b/power-grid.csv and cheapest-designs.csv. Each preset and fold change has a k = 1 row; carrier-TMT presets also have k = 0.75 and 0.5 rows. The maximum grid power is reported with its biological samples per condition, cells per sample, nominal instrument hours and simulation SE. Ties in maximum power select fewer hours, then fewer cells, then fewer samples. The cheapest registered design reaching mean power 0.8 is reported separately, or as none in grid. The table does not extrapolate beyond the registered grid or replace unconditional power with conditional power. Markdown values are displayed to three decimals; CSV retains full precision.','',
 'TOC graphic. Reduced Figure 3e using exactly the same registered 2× grid points, preset colours, logarithmic instrument-hour axis and circled minimum-hour point reaching 0.8. The question describes the package purpose; no additional design calculation is shown.')
writeLines(captions,file.path(out,'captions.md'))
omissions<-paste('No requested figure panel was omitted for lack of registered numerical outputs.',
    'The timsTOF real-cell comparison remains absent because the registered M4 report classifies it as not evaluable; the Figure 4 caption states the source-cell limitation.',
    'Table 1 marks the two TMT instrument model names as unrecorded because the allowed reports and copied preset metadata contain no model names; Brunner retains the recorded TIMS04 identifier without inferring a model.',
    'All plotted estimates come from existing reports, with only the requested SD conversion, tabular summaries and deterministic display binning. No calibration, simulation, model or criterion was changed.')
writeLines(omissions,file.path(out,'omissions.md'))
writeLines(c('# Manuscript figure bundle','',
    'Reproduce from the repository root: `Rscript inst/scripts/figures_manuscript.R`.',
    'The script reads only the allowed registered report directories. It writes no source data.',
    'PDFs are vector artwork. PNGs are 600 dpi. Figures are 7 inches wide; the TOC is 3.25 × 1.75 inches.',
    'Every plotted data panel has an adjacent CSV; extra protein, curve and band CSVs preserve its source coordinates.',
    '[Caption drafts](captions.md), [availability note](omissions.md), [Table 1](Table1.md), [Table 2](Table2.md), [run manifest](run-manifest.json).'),file.path(out,'README.md'))
formats<-data.frame(figure=c('Fig1','Fig2','Fig3','Fig4','TOC'),
    width_inches=c(7,7,7,7,3.25),height_inches=c(4.1,8.2,9,7.7,1.75),dpi=600)
write_data(formats,'figure-dimensions.csv')
# Check the actual PDF font sizes when Poppler is available; the draw sizes
# above also enforce the minimum independently of this inspection utility.
font_inspection<-list()
if(nzchar(Sys.which('pdftohtml')))for(name in formats$figure) {
    target<-tempfile('manuscript-fonts-')
    status<-system2('pdftohtml',c('-xml','-hidden','-zoom','1',
        shQuote(file.path(out,paste0(name,'.pdf'))),shQuote(target)),stdout=FALSE,stderr=FALSE)
    if(status!=0L)stop('PDF inspection failed for ',name)
    xml<-paste(readLines(paste0(target,'.xml'),warn=FALSE),collapse='\n')
    tags<-regmatches(xml,gregexpr('<fontspec[^>]*size="[0-9.]+"',xml,useBytes=TRUE))[[1]]
    sizes<-as.numeric(sub('.*size="([0-9.]+)".*','\\1',tags,useBytes=TRUE))
    stopifnot(length(sizes)>0,min(sizes)>=7)
    font_inspection[[name]]<-list(minimum_points=min(sizes),font_sizes=sort(unique(sizes)))
    unlink(paste0(target,'.xml'))
}
outputs<-list.files(out,pattern='\\.(pdf|png|csv|md)$',full.names=TRUE)
manifest<-list(input_hashes=inputs,seed=NULL,randomness_used=FALSE,
    versions=as.list(versions),rendering_versions=list(ggplot2=as.character(utils::packageVersion('ggplot2')),
        cairo=grDevices::grSoftVersion()[['cairo']]),
    wall_seconds=proc.time()[['elapsed']]-started,formats=formats,
    minimum_text_points=7,pdf_font_inspection=font_inspection,preset_palette=as.list(palette),fold_change_palette=as.list(fold_palette),
    input_scope=allowed,scientific_source_files_written=FALSE,new_calibration=FALSE,new_simulations=FALSE,
    numerical_transformations=c('square roots of recorded variances','requested medians and coverage summaries',
        'deterministic 2-D count bins','selection of registered grid maxima and cheapest designs',
        'identity tolerance bands and recorded logistic probability complement'),
    output_sha256=setNames(lapply(outputs,function(p)digest::digest(file=p,algo='sha256')),basename(outputs)),
    availability_note=omissions)
jsonlite::write_json(manifest,file.path(out,'run-manifest.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA,na='null')
cat('Manuscript figures, panel CSVs, tables, captions and manifest complete.\n')
