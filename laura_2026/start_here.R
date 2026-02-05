# Code by Joy Nakato and Laura Balzer

rm(list=ls())

library(dplyr)
library(tidyverse)
# library(usedist)
library('MASS')
library('SuperLearner')
library('ltmle')
# library('kableExtra')
library('earth')

library(MRStdCRT) # for GEE

## load functions to generate the data
source('generate_data_2feb.R')
## load functions for stage 1 estimation
source('stage1.R')
## load functions for doing stage 2 estimation
source('stage2/stage2.R')
source('stage2/tmle.R')
source('stage2/aps.R')
## load functions for single stage estimators
source('single_stage.R')
## load functions for getting performance of the estimator
source('run_estimators.R')

# set seed here.
set.seed(423)

N_mean <- 200
N_sd <- 10

#  specify the DGP
dgp <- "null"

calculate_causal_estimand <- F

J_truth <- 5000
this_date <- '_2feb2026'
file_truth <- paste0('truth_dgp_', dgp, '_J', J_truth, this_date, '.Rdata')
# 
if(calculate_causal_estimand){
  #  calculate truth here
  truth <- get_truth(dgp=dgp, J=J_truth, N_mean=N_mean, N_sd=N_sd, verbose=F)
  truth
  save(truth,  generate_cluster, get_delta, get_Y2, file=file_truth)
} else{
  load(file_truth)
}

psi <- as.numeric(truth[1] - truth[2])

verbose=F

# Specify Super Learner library
SL.library <- NULL
# SL.library<- c("SL.glm", "SL.earth", "SL.mean")

J <- 70

R <- 1000

file_out <- paste0(dgp, '_J', J, '_rep', R,  this_date, '.Rdata')

screened <- eligible <- unadj <- tmle <- tmle_tmle <-  
  gee <- stmle_eligible <-  stmle_ratio <-  NULL

for(r in 1:R){
  out <- looper(dpg=dgp, J=J, N_mean=N_mean, N_sd=N_sd, 
                psi=psi,
                SL.library=SL.library, verbose=F)
  screened <-  rbind(screened, out$screened)
  eligible <-  rbind(eligible, out$eligible)
  unadj <-  rbind(unadj, out$unadj)
  tmle <- rbind(tmle, out$tmle)
  tmle_tmle <- rbind(tmle_tmle, out$tmle_tmle)
  gee <- rbind(gee, out$gee)
  stmle_eligible <- rbind(stmle_eligible, out$stmle_eligible)
  stmle_ratio <- rbind(stmle_ratio, out$stmle_ratio)
  print(r)
}

# reformat TMLE-TMLE with CV variance
# these are the only columns that we really need
tmle_tmle_cv <- tmle_tmle[,c("CV.est","CV.CI.lo","CV.CI.hi","CV.se", "CV.pval")]
colnames(tmle_tmle_cv) <- c('est', 'CI.lo', 'CI.hi', 'se', 'pval')
save(screened,  eligible, unadj, tmle, tmle_tmle,  tmle_tmle_cv,
     gee, stmle_eligible, stmle_ratio,
     generate_cluster, get_delta, get_Y2,
     file=file_out)


these_cols <- c('bias','cover', 'reject')
print( colMeans(screened[,these_cols]))
print( colMeans(eligible[,these_cols]))
print( colMeans(unadj[,these_cols]))
print( colMeans(tmle[,these_cols]) )
print( colMeans(tmle_tmle[,these_cols]) )
print( colMeans(tmle_tmle[,c('CV.bias','CV.cover','CV.reject')]) )

print( colMeans(gee) )
print( colMeans(stmle_eligible))
print( colMeans(stmle_ratio) )



# these_cols <- c('bias','cover', 'reject')
# print( colMeans(screened[,these_cols]))
# print( colMeans(eligible[,these_cols]))
# print( colMeans(unadj[,these_cols]))
# print( colMeans(tmle[,these_cols]) )
# # print( colMeans(tmle_tmle[,these_cols]) )
# #print( colMeans(gee[,these_cols]) )
# print( colMeans(stmle_eligible[,these_cols])) 
# print( colMeans(stmle_ratio[,these_cols]) )


# get_performance <- function(X){
#   data.frame(cbind(
#     est= mean(X$est),
#     ci_lo = mean(X$CI.lo),
#     ci_hi = mean(X$CI.hi),
#     bias=mean(X$bias),
#     ave_se = mean(X$se),
#     sd_pt = sd(X$est),
#     cover = mean(X$cover),
#     power = mean(X$reject)
#     
#   ))
#   
#   
# }
# 
# results <- data.frame( rbind(
#   standard=get_performance(standard),
#   gee_std_rd=get_performance(gee_std_rd),
#   raw_all_unadj=get_performance(raw_all_unadj),
#   raw_eligible_unadj=get_performance(raw_eligible_unadj),
#   unadj_unadj=get_performance(unadj_unadj),
#   glm_unadj=get_performance(glm_unadj),
#   glm_glm=get_performance(glm_glm),
#   tmle_unadj_rd=get_performance(tmle_unadj_rd),
#   tmle_tmle_rd=get_performance(tmle_tmle_rd)
#   
#   
# ))
# 
# 
# file.name1 <- paste('dgp_',dgp, '_goal_',goal,'_R',R, '_J', J, '_N', N_mean, '_effect', effect,
#                     '_date_', format(Sys.Date(), "%Y%m%d"), '.Rdata', sep='' )
# save(results, raw_all_unadj, raw_eligible_unadj,standard,gee_std_rd,
#      unadj_unadj, 
#      glm_unadj, glm_glm,
#      tmle_unadj_rd, tmle_tmle_rd, 
#      truth, SL.library, Q_bar_adj,pscore_adj, Q_form, G_form,file=file.name1)
# 
