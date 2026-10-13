simulation_fixture <- function() list(
    tiers = data.frame(tier = c(1,2), proteins = c(20,20), mean_log2 = c(8,12),
        sigma2_s = c(.2,.3), sigma2_b = c(.1,.2), sigma2_e = c(.3,.4)),
    curve = list(intercept=1,effective_slope=.6,intensity_center=10),
    plex_size=8L,label_free=FALSE,k=1)

 test_that("seeded simulation is reproducible and preserves caller RNG state", {
    p <- simulation_fixture();set.seed(721);before <- .Random.seed;kind <- RNGkind()
    a <- scp_simulate(NULL,2,7,beta=1,preset=p,seed=91)
    b <- scp_simulate(NULL,2,7,beta=1,preset=p,seed=91)
    expect_identical(a$matrix,b$matrix)
    expect_identical(a$cell_metadata,b$cell_metadata)
    expect_identical(a$protein_metadata,b$protein_metadata)
    expect_identical(a$truth,b$truth)
    expect_identical(.Random.seed,before);expect_identical(RNGkind(),kind)
    expect_equal(dim(a$matrix),c(40,28))
    expect_equal(as.integer(table(a$cell_metadata$sample)),rep(7L,4))
    expect_equal(sum(a$protein_metadata$changed),4)
    expect_equal(a$manifest$seed,91L)
    expect_true(all(c('R','scp','limma','edgeR','scplainer')%in%names(a$manifest$versions)))
    expect_true(all(nchar(a$manifest$input_hashes)==64L))
 })

 test_that("missingness follows expected intensity, independently of cell residual", {
    p <- simulation_fixture();p$tiers$proteins<-c(100,100);p$tiers$mean_log2<-10
    p$tiers$sigma2_s<-p$tiers$sigma2_b<-0;p$tiers$sigma2_e<-25
    p$curve<-list(intercept=0,effective_slope=1,intensity_center=10)
    a<-scp_simulate(NULL,1,100,beta=0,preset=p,seed=93)
    expect_equal(range(a$expected_log2),c(10,10))
    expect_false(any(a$protein_metadata$changed))
    expect_equal(mean(is.na(a$matrix)),.5,tolerance=.02)
    q<-a$matrix;q[is.na(q)]<-a$expected_log2[is.na(q)]
    # Direct expected-intensity semantics are also checked through low/high
    # expected intensities while the realized residual SD remains large.
    p$tiers$mean_log2<-c(7,13)
    b<-scp_simulate(NULL,1,100,beta=0,preset=p,seed=93)
    expect_gt(mean(is.na(b$matrix[1:100,])),.9)
    expect_lt(mean(is.na(b$matrix[101:200,])),.1)
    p$curve$intercept<-4
    c<-scp_simulate(NULL,1,100,beta=0,preset=p,seed=93)
    expect_lt(mean(is.na(c$matrix)),mean(is.na(b$matrix)))
 })

 test_that("balanced and confounded layouts change batch-induced contrasts", {
    p<-simulation_fixture();p$plex_size<-20;p$tiers$proteins<-c(200,200)
    p$tiers$sigma2_s<-p$tiers$sigma2_e<-0;p$tiers$sigma2_b<-9
    p$curve<-list(intercept=40,effective_slope=0,intensity_center=10)
    balanced<-scp_simulate(NULL,3,20,beta=0,preset=p,layout='balanced',seed=94)
    confounded<-scp_simulate(NULL,3,20,beta=0,preset=p,layout='confounded',seed=94)
    contrast<-function(r)rowMeans(r$matrix[,r$cell_metadata$condition=='B'])-rowMeans(r$matrix[,r$cell_metadata$condition=='A'])
    expect_lt(max(abs(contrast(balanced))),1e-12)
    expect_gt(sd(contrast(confounded)),1)
    expect_true(all(vapply(split(balanced$cell_metadata$condition,balanced$cell_metadata$batch),function(v)sum(v=='A')==sum(v=='B'),logical(1))))
    expect_true(all(vapply(split(confounded$cell_metadata$condition,confounded$cell_metadata$batch),function(v)length(unique(v))==1L,logical(1))))
 })

 test_that("sample effects are shared across cells and k scales only the condition effect", {
    p<-simulation_fixture();p$tiers$sigma2_b<-p$tiers$sigma2_e<-0
    p$curve<-list(intercept=40,effective_slope=0,intensity_center=10)
    a<-scp_simulate(NULL,2,10,beta=2,preset=p,k=1,seed=97,changed_fraction=1)
    b<-scp_simulate(NULL,2,10,beta=2,preset=p,k=.5,seed=97,changed_fraction=1)
    expect_equal(dim(a$truth$sample_effects),c(40,4))
    expect_true(all(vapply(split(seq_len(nrow(a$cell_metadata)),a$cell_metadata$sample),function(j)all(apply(a$matrix[,j],1,sd)==0),logical(1))))
    conditionB<-a$cell_metadata$condition=='B'
    expect_equal(unname(a$matrix[,conditionB]-b$matrix[,conditionB]),matrix(1,40,20))
    expect_identical(a$matrix[,!conditionB],b$matrix[,!conditionB])
 })

 test_that("custom layouts retain native unequal designs and label-free absorption", {
    p<-simulation_fixture()
    layout<-matrix(c(2,0,1, 0,3,1),3,2,dimnames=list(c('L1','L2','L3'),c('run1','run2')))
    r<-scp_simulate(NULL,c(2,1),c(2,3,2),beta=0,preset=p,layout=layout,seed=99)
    expect_equal(dim(r$matrix),c(40,7))
    expect_equal(unclass(table(factor(r$cell_metadata$sample,levels=rownames(layout)),factor(r$cell_metadata$batch,levels=colnames(layout)))),unname(layout),ignore_attr=TRUE)
    p$label_free<-TRUE;p$tiers$sigma2_b<-0
    native<-matrix(1,1,8,dimnames=list('HeLa',paste0('run',1:8)))
    r<-scp_simulate(NULL,c(1,0),8,beta=0,preset=p,layout=native,seed=99)
    expect_equal(length(unique(r$cell_metadata$batch)),8)
    expect_true(all(r$truth$batch_effects==0))
    expect_error(scp_simulate(NULL,c(1,0),8,beta=1,preset=p,layout=native),'beta=0')
    expect_error(scp_simulate(NULL,c(1,0),8,beta=0,preset=p,layout=matrix(8,1,1)),'one cell per batch')
    expect_error(scp_simulate(NULL,1,c(2,3),0,preset=simulation_fixture()),'equal total')
 })

 test_that("YAML user presets and calibration objects warn on version mismatch", {
    p<-simulation_fixture();path<-tempfile(fileext='.yaml');yaml::write_yaml(p,path)
    a<-scp_simulate(NULL,1,3,0,preset=p,seed=101)
    b<-scp_simulate(NULL,1,3,0,preset=path,seed=101)
    expect_identical(a$matrix,b$matrix)
    p$versions<-list(R='0.0')
    expect_warning(r<-scp_simulate(NULL,1,3,0,preset=p,seed=101),'mismatch.*continuing')
    expect_equal(dim(r$matrix),c(40,6))
    set.seed(882)
    z<-matrix(rnorm(12*30),12,30)+seq_len(12)
    z[matrix(runif(length(z)),12,30)<.15]<-NA_real_
    dimnames(z)<-list(paste0('P',1:12),paste0('C',1:30))
    md<-data.frame(sample=rep(paste0('S',1:3),each=10),batch=rep(rep(paste0('B',1:5),each=2),3),row.names=colnames(z))
    cal<-scp_calibrate(z,md,sample_col='sample',batch_col='batch')
    cal$tiers$sigma2_s<-0
    r<-scp_simulate(cal,1,5,0,seed=101)
    expect_equal(dim(r$matrix),c(12,10))
    multi<-scp_simulate(list(custom=cal),1,5,0,preset='custom',seed=101)
    expect_identical(r$matrix,multi$matrix)
    p$versions<-NULL
    p$tiers$sigma2_s[1]<-NA_real_
    expect_error(scp_simulate(NULL,1,3,0,preset=p),'Unidentified')
 })
