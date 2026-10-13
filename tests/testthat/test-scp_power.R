power_fixture <- function() {
 p<-list(tiers=data.frame(tier=1,proteins=100,mean_log2=10,
     sigma2_s=.05,sigma2_b=.03,sigma2_e=.1),
     curve=list(intercept=10,effective_slope=0,intensity_center=10),
     label_free=FALSE,plex_size=8,cells_per_run=8,minutes_per_run=60)
 p$protein_parameters<-data.frame(protein=paste0('real',1:100),tier=1,
     mean_log2=seq(8,12,length.out=100),sigma2_s=.05,sigma2_b=.03,sigma2_e=.1,r_g=1)
 p$detection_mode<-'protein_rate';p
}

test_that('M3 endpoint requires three quantified cells in every sample',{
 p<-power_fixture()
 r<-scp_power(p,3,c(1,2,4),beta=3,replicates=4,seed=119)
 expect_equal(r$results$testability,c(0,0,1))
 expect_equal(r$results$power[1:2],c(0,0))
 expect_gt(r$results$power[3],.95)
 expect_equal(r$results$realised_fdr[1:2],c(0,0))
 expect_true(all(r$replicate_metrics$testable[r$replicate_metrics$m_cells_per_sample<3]==0))
 p$protein_parameters$r_g[1:50]<-0
 r<-scp_power(p,3,4,3,replicates=2,seed=121)
 expect_equal(r$results$testability,.5,tolerance=.15)
})

test_that('M3 records reproducibility, Monte Carlo SE, FDR and actual run hours',{
 p<-power_fixture()
 a<-scp_power(p,c(3,4),c(4,8),beta=1,k=c(1,.5),replicates=5,seed=123)
 b<-scp_power(p,c(3,4),c(4,8),beta=1,k=c(1,.5),replicates=5,seed=123)
 expect_identical(a$replicate_metrics$matrix_sha256,b$replicate_metrics$matrix_sha256)
 expect_equal(a$results$power,b$results$power)
 expect_equal(nrow(a$results),8)
 r<-a$replicate_metrics[a$replicate_metrics$grid_index==1,]
 expect_equal(a$results$simulation_se[1],sd(r$power)/sqrt(5))
 expect_equal(a$results$realised_fdr[1],mean(r$realised_fdr))
 expect_equal(a$results$runs[1],3)
 expect_equal(a$results$instrument_hours,a$results$runs)
 expect_true(all(c('R','scp','limma','edgeR','scplainer')%in%names(a$manifest$versions)))
 expect_equal(a$replicate_metrics$seed[1:5],123L+1:5)
 expect_equal(a$replicate_metrics$seed[6:10],100123L+1:5)
 expect_true(all(nchar(a$replicate_metrics$matrix_sha256)==64))
 expect_true(all(a$replicate_metrics$wall_seconds>=0))
 expect_true(all(is.finite(a$results$simulation_se)))
 expect_true(all(a$results$realised_fdr>=0&a$results$realised_fdr<=1))
 expect_match(a$manifest$endpoint$inference,'limma-trend')
})

test_that('limma BH FDR counts null rejections and sign agreement affects power',{
 set.seed(127)
 y<-matrix(rnorm(12*24,sd=.02),12,24)
 md<-data.frame(sample=rep(paste0('S',1:6),each=4),
    condition=rep(c('A','B'),each=12),batch=rep(paste0('B',1:3),each=8))
 y[9:12,md$condition=='B']<-y[9:12,md$condition=='A']
 effect<-c(3,3,-3,-3,3,3,3,3,rep(0,4));y[,md$condition=='B']<-y[,md$condition=='B']+effect
 truth<-data.frame(changed=c(rep(TRUE,4),rep(FALSE,8)),beta_observed=c(rep(3,4),rep(0,8)))
 endpoint<-scpPower:::.scp_primary_endpoint(list(matrix=y,cell_metadata=md,protein_metadata=truth))
 expect_equal(endpoint$power,.5)
 expect_equal(endpoint$detected,2)
 expect_equal(endpoint$wrong_sign_changed,2)
 expect_equal(endpoint$false_discoveries,4)
 expect_equal(endpoint$rejected,8)
 expect_equal(endpoint$realised_fdr,.5)
})

test_that('instrument accounting preserves odd balanced plexes and label-free runs',{
 p<-power_fixture();p$plex_size<-3;p$cells_per_run<-3;p$minutes_per_run<-30
 r<-scp_power(p,3,4,1,replicates=1,seed=129)
 expect_equal(r$results$runs,12)
 expect_equal(r$results$instrument_hours,6)
 expect_true(is.na(r$results$simulation_se))
 p$label_free<-TRUE;p$plex_size<-1;p$cells_per_run<-1;p$tiers$sigma2_b<-0;p$protein_parameters$sigma2_b<-0
 r<-scp_power(p,3,4,1,replicates=1,seed=129)
 expect_equal(r$results$runs,24)
 expect_equal(r$results$instrument_hours,12)
})

test_that('target design selection reports unavailable targets and software mismatches',{
 p<-power_fixture()
 r<-scp_power(p,3,c(1,4,8),3,replicates=2,seed=131)
 expect_equal(r$cheapest$m_cells_per_sample,4)
 expect_equal(r$cheapest$target_power,.8)
 r<-scp_power(p,3,1,1,replicates=2,seed=131)
 expect_match(r$cheapest$status,'No tested design')
 expect_true(is.na(r$cheapest$instrument_hours))
 p$versions<-list(R='0.0')
 expect_warning(r<-scp_power(p,3,4,1,replicates=2,seed=133),'mismatch.*continuing')
 expect_equal(nrow(r$results),1)
 p$minutes_per_run<-NULL
 expect_error(scp_power(p,3,4,1),'minutes_per_run')
 expect_error(scp_power(power_fixture(),1,4,1),'valid integers')
})
