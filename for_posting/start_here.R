# Code by Joy Nakato & Laura B. Balzer
# jznakato@berkeley.edu & laura.balzer@berkeley.edu

rm(list=ls())

library(dplyr)
library(tidyverse)
library('MASS')
library('ltmle')

library(MRStdCRT) # for GEE

## load functions to generate the data
source('generate_data.R')
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
# dgp <- "main" # with an effect
dgp <- "null" # under the null

calculate_causal_estimand <- F

J_truth <- 5000
file_truth <- paste0('truth_dgp_', dgp, '_J', J_truth, '.Rdata')

if(calculate_causal_estimand){
  #  calculate truth here
  truth <- get_truth(dgp=dgp, J=J_truth, N_mean=N_mean, N_sd=N_sd, verbose=F)
  truth
  save(truth,  generate_cluster, get_delta, get_Y2, file=file_truth)
} else{
  load(file_truth)
}

# focusing on the absolute effect
# truth contains the expected cf outcomes: E[Yc(1)] & E[Yc(0)]
psi <- as.numeric(truth[1] - truth[2])

verbose=F

# vary the number of clusters
J <- 20

# number of iterations
R <- 1000

file_out <- paste0(dgp, '_J', J, '_rep', R, '.Rdata')

screened <- eligible <- unadj <- tmle <- tmle_tmle <-  
  gee <- stmle_eligible <- NULL

for(r in 1:R){
  out <- looper(dpg=dgp, J=J, N_mean=N_mean, N_sd=N_sd, 
                psi=psi, verbose=F)
  screened <-  rbind(screened, out$screened)
  eligible <-  rbind(eligible, out$eligible)
  unadj <-  rbind(unadj, out$unadj)
  tmle <- rbind(tmle, out$tmle)
  tmle_tmle <- rbind(tmle_tmle, out$tmle_tmle)
  gee <- rbind(gee, out$gee)
  stmle_eligible <- rbind(stmle_eligible, out$stmle_eligible)
  print(r)
}

# for TMLE-TMLE, grab the CV variance columns
tmle_tmle <- tmle_tmle[,c("CV.est","CV.CI.lo","CV.CI.hi","CV.se", "CV.pval")]
colnames(tmle_tmle) <- c('est', 'CI.lo', 'CI.hi', 'se', 'pval')

# save everything
save(screened, eligible, unadj, tmle, tmle_tmle, 
     gee, stmle_eligible,
     generate_cluster, get_delta, get_Y2,
     file=file_out)
