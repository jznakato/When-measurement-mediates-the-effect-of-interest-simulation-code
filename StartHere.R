
## Full script of opal manuscript revised and resubmitted version to statmed. 

##code adapted from TwoStageTMLE/StartHere.R 
#https://github.com/LauraBalzer/TwoStageTMLE/blob/main/StartHere.R
rm(list=ls())


#install.packages('earth')
library(dplyr)
library(tidyverse)
library(usedist)
library('nbpMatching')
library('MASS')
library('SuperLearner')
library('ltmle')
library('kableExtra')
library('earth')


## load functions to generate the data
source('generate_data.R')
## load functions for stage 1 estimation
source('do.stage.1.R')
## load functions for doing stage 2 estimation
source('stage2_functions.R')
source('tmle_functions.R')
source('adapt_functions.R')
## load functions for getting performance of the estimator
source('estimator.performance.adapt.meta.R')

# set seed here.
set.seed(423)



##do stage 1
R <-1000
J <-50 # 100, 150
# LB june 6: set indv sample size here
N.mean <-  200
N.sd <- 10
effect <- T
pooled <- F

#  specify the DGP
dgp <- 'complex'  # "simple"
# difference or ratio
goal <- 'RD' # 'aRR'

# Specify target for estimation i.e cluster-level  or individual level effect
target="clust"
verbose=F

# Specifiy superlearner library
SL.library<- c("SL.glm", "SL.earth", "SL.mean","SL.gam")


#  calculate truth here
Jtruth <- 5000

this_date <- '21st_Nov_25_adj_all'

filename_truth <- paste("R_",R,'truth_dgp_',dgp, '_J', Jtruth, '_N', N.mean, '_effect', effect,  this_date, '.Rdata', sep='' )

if(TRUE){
  truth <- get.truth( dgp=dgp, J=Jtruth, N.mean=N.mean, N.sd=N.sd, effect=effect, verbose=F)
# also returns the cluster-specific truth
# could be useful later if we switch to sample (vs. population effects)
  truth <- truth$truth
  save(truth, file=filename_truth)
} 
load(filename_truth)

Psi.P0.RD<-truth$RD.clust

Psi.P0.RR<-truth$RR.clust
# Risk prevalence (1-true HIV prevalence)
Psi_den <-truth$Psi_den

## create variables for result estimates
raw_all_unadj <- raw_eligible_unadj <- unadj_unadj <-
  glm_unadj <- glm_glm <-tmle_unadj_rd<-tmle_tmle_rd<-tmle_unadj_rr<-tmle_tmle_rr<-aug_rr<-standard<-gee_std_rd<-gee_std_rr<-glmm_std_rr<- NULL
Q_bar_adj<-rep(NA,R)
pscore_adj<-rep(NA,R)
Q_form<-rep(NA,R)
G_form<-rep(NA,R)

## Iterate over R simulations
for(r in 1:R){
  print(r)
  OUT <- looper(dpg=dgp, J=J, N.mean=N.mean, N.sd=N.sd, effect=effect, goal=goal, target=target, 
                Psi.P0.RD=Psi.P0.RD,Psi.P0.RR=Psi.P0.RR, Psi_den=Psi_den, SL.library=SL.library, verbose=verbose)
  standard<-rbind(standard,OUT$standard)
  gee_std_rd<-rbind(gee_std_rd,OUT$GEE_std_RD)
  raw_all_unadj <- rbind(raw_all_unadj, OUT$raw_all_unadj)
  raw_eligible_unadj <- rbind(raw_eligible_unadj, OUT$raw_eligible_unadj)
  unadj_unadj <- rbind(unadj_unadj, OUT$unadj_unadj)
  glm_unadj <- rbind(glm_unadj, OUT$glm_unadj)
  glm_glm <- rbind(glm_glm, OUT$glm_glm)
  tmle_unadj_rd <- rbind(tmle_unadj_rd, OUT$tmle_unadj_RD)
  tmle_tmle_rd <- rbind(tmle_tmle_rd, OUT$tmle_tmle_RD)
  gee_std_rr<-rbind(gee_std_rr,OUT$GEE_std_RR)
  glmm_std_rr<-rbind(glmm_std_rr,OUT$GLMM_std_RR)
  aug_rr<-rbind(aug_rr,OUT$aug_RR)
  tmle_unadj_rr <- rbind(tmle_unadj_rr, OUT$tmle_unadj_RR)
  tmle_tmle_rr <- rbind(tmle_tmle_rr, OUT$tmle_tmle_RR)


# Stage 2 adjustment variables for Q and g
  Q_bar_adj[r]<-tmle_tmle_rd$QAdj
pscore_adj[r]<-tmle_tmle_rd$gAdj
# Stage 2  functional forms for modelling  Q and g
 Q_form[r]<-tmle_tmle_rd$Qform  
 G_form[r]<-tmle_tmle_rd$gform  
}


get_performance <- function(X){
  data.frame(cbind(
    est= mean(X$est),
    ci_lo = mean(X$CI.lo),
    ci_hi = mean(X$CI.hi),
    bias=mean(X$bias),
    ave_se = mean(X$se),
    sd_pt = sd(X$est),
    cover = mean(X$cover),
    power = mean(X$reject)
    
  ))
  
  
}

results <- data.frame( rbind(
  standard=get_performance(standard),
  gee_std_rd=get_performance(gee_std_rd),
  raw_all_unadj=get_performance(raw_all_unadj),
  raw_eligible_unadj=get_performance(raw_eligible_unadj),
  unadj_unadj=get_performance(unadj_unadj),
  glm_unadj=get_performance(glm_unadj),
  glm_glm=get_performance(glm_glm),
  tmle_unadj_rd=get_performance(tmle_unadj_rd),
  tmle_tmle_rd=get_performance(tmle_tmle_rd),
  gee_std_rr=get_performance(gee_std_rr),
  glmm_std_rr=get_performance(glmm_std_rr),
  aug_rr=get_performance(aug_rr),
  tmle_unadj_rr = get_performance(tmle_unadj_rr),
  tmle_tmle_rr = get_performance(tmle_tmle_rr)
  
  
))


file.name1 <- paste('dgp_',dgp, '_goal_',goal,'_R',R, '_J', J, '_N', N.mean, '_effect', effect,
                    '_pooled',pooled,'_Y1_type',Y1_type,'_date_', format(Sys.Date(), "%Y%m%d"), '.Rdata', sep='' )
save(results, raw_all_unadj, raw_eligible_unadj,standard,gee_std_rd,gee_std_rr,glmm_std_rr,
     unadj_unadj, 
     glm_unadj, glm_glm,
     tmle_unadj_rd, tmle_tmle_rd,  tmle_unadj_rr, tmle_tmle_rr,aug_rr,
     truth, SL.library, Q_bar_adj,pscore_adj, Q_form, G_form,file=file.name1)

