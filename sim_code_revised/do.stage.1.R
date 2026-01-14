
##=================================================================================
# getYc: Function to get cluster level endpoint in Stage 1 of TMLE
## Inputs:
# O- observed data
# ind.cov: Individual level covariates
# clust.var: Cluster level variables 
## output: Returns a cluster-level data set after iterating over all clusters

##=================================================================================
getYc<- function(O, ind.cov, clust.var){
  
  clusters <- unique(O$id)
  nClust <- length(clusters)
  
  data.clust <- NULL
  
  for(j in 1:nClust){ 
# for each cluster separately
    
# subset the data on the cluster of interest 
    these.units <- clusters[j]==O$id
    Oc <- O[these.units, ]
    
# doing the numerator first P(Y2=1)
    suppressMessages(
      tmle_num <-  ltmle(data=Oc[,c('U','Y2')], 
                         Anodes='U', Ynodes='Y2',abar = 1,
                         estimate.time=F)
    )
    Psi_num<- as.numeric(tmle_num$estimates['tmle'])

    
# get the influence function for the numerator P(Y2=1)
    IC_num<-tmle_num$IC$tmle 
    
# unadj = P(Y2=1)/P(Y1=0| Delta=1)
unadj <- get_denominator_ratio(Oc, estimator='unadj', Psi_num=Psi_num, IC_num=IC_num)
    
# glm and TMLE = P(Y2=1)/ E[ P(Y1=0|Delta=1,W)]~BW
glm <- get_denominator_ratio(Oc, estimator='glm', ind.cov=ind.cov, 
                                 Psi_num=Psi_num, IC_num=IC_num)
    
tmle <- get_denominator_ratio(Oc, estimator='tmle', ind.cov=ind.cov, SL.library=SL.library,
                                  Psi_num=Psi_num, IC_num=IC_num)
    
# "raw_all" ("Screened"): P(Y2=1 | Delta=1) - this does not have a numerator/denominator
# this estimates a different parameter - does not account for eligibility 
raw_all <- among_linkers(Oc, among_eligible=F)
    
# "raw_eligible"("eligible"): P(Y2=1 | Delta=1, Y1=0
raw_eligible <- among_linkers(Oc, among_eligible=T)
  
    temp <- cbind(Oc[1, clust.var], W1=mean(Oc$W1),  W2=mean(Oc$W2),  W3=mean(Oc$W3),
                  unadj, glm, tmle, raw_all,raw_eligible)
    data.clust <-  rbind(data.clust, temp)
    
    
  }  
  data.clust
  
}

##=================================================================================
# get_denominator_ratio: This function estimates the demonimator portion of the parameter 
# using an unadjusted estimator, a simple glm adjusting for indv covariates or a doubly
# robust estimator (TMLE)
## Inputs
# -estimator; unadajusted, glm or tmle
# ind.cov: individual level covariates to adjust for
# SL.library: SuperLearner library
# Psi_num: estimate of the numerator
# IC_num: Influence function of the numerator
##=================================================================================

get_denominator_ratio <- function(Oc, estimator, ind.cov, SL.library,
                                  Psi_num, IC_num){
  
  if(estimator=='unadj'){
    # unadjusted estimator (1-HIV prevalence): mean among those measured 
    # ## adjust for nothing
    
    # ## several different implementations; all give same point estimate, but want
    # ## 2nd implementation for parallel structure and preserve length of ICs
    
    # option1 - subset on measured
    # Oc.sub <- Oc[Oc$Delta==1,c('U','Y1')]
    # tmle_den<-  suppressMessages( ltmle(data=Oc.sub, Anodes="U", Ynodes='Y1',abar = 1,
    #                   estimate.time=F))
    # option 2 - set adjustment covariates = U and intervene on Delta
    Oc.sub <- Oc[,c('U','Delta','Y1')]
    tmle_den<-  suppressMessages(  ltmle(data=Oc.sub, Anodes='Delta', 
                                         Ynodes='Y1', abar=1,
                                         estimate.time=F, stratif=T) )
    
  } else if (estimator=='glm') {
    # in two-stage we dont adjust for cluster-level covariates; they are constant 
    # this will run a glm
    Oc.sub <- Oc[,c(ind.cov,'Delta','Y1')]
    tmle_den<-  suppressMessages(  ltmle(data=Oc.sub, Anodes='Delta', 
                                         Ynodes='Y1', abar=1,
                                         estimate.time=F, stratif=T) )
    
  } else if (estimator=='tmle') {
    # run TMLE with SL
    # individual-level TMLE to adjust for differential measurement (using superlearner)
    Oc.sub <- Oc[,c(ind.cov,'Delta','Y1')]
    
    tmle_den<- suppressMessages(ltmle(data=Oc.sub, Anodes='Delta', Ynodes='Y1', abar=1 ,
                                      SL.library=SL.library,
                                      estimate.time=F, stratif=T))
  } 
  
  # get inference for HIV risk prevalence
  den <- 1-tmle_den$estimates['tmle'] 
  den_lo <-  1-summary(tmle_den)$treatment$CI[2]
  den_hi <-  1-summary(tmle_den)$treatment$CI[1]
  
  IC_den<- tmle_den$IC$tmle
  
  # NOW DO THE RATIO FOR POINT ESTIMATE AND INFERENCE
  # get IC-based 95%CI for TMLEs
  Yc <- get.var.delta(psi_num=Psi_num, psi_den=den, 
                      IC1=IC_num, IC0=IC_den, alpha=0.05)$est
  yay <- data.frame(den, den_lo, den_hi, Yc)
  colnames(yay) <- paste0(estimator, '_',colnames(yay))
  yay 
  
}


##=================================================================================
# get.var.delta : Applies the delta method to get inference for the ratio
## Inputs:
# psi_num: estimate of the numerator
# psi_den: estimate of the denominator
# IC1: Influence Curve of the numerator
# IC0: Influence cure of the denominator
# alpha: level of significance
## output: A data frame with point estimates and confidence intervals
##=================================================================================

get.var.delta <- function(psi_num, psi_den, IC1, IC0, alpha=0.05){
  
  # get inference via the delta method on log scale
  psi <- log( psi_num/psi_den)
  IC <- (1/psi_num)*(IC1) - (1/psi_den)*(IC0)
  
  # variance of asymptotically  linear estimator  is var(IC)/n
  var<- var(IC)/length(IC)
  # testing and CI	
  cutoff <- qnorm(alpha/2, lower.tail=F)
  se <- sqrt(var)
  CI.lo <- psi - cutoff*se
  CI.hi <- psi + cutoff*se
  
  est <- data.frame(pt=exp(psi), lo=exp(CI.lo), hi=exp(CI.hi) ) 
  
  list(est=est, IC=IC)
}


# Jul8 - new function to get raw uptake among linkers
# overall and among those who are eligible

among_linkers <- function(Oc, among_eligible){
  
  if(among_eligible){
    Oc.sub <- Oc[Oc$Delta==1 & Oc$Y1==0,c('U','Y2')]
    col_start <- 'raw_eligible_'
  } else{
    Oc.sub <- Oc[Oc$Delta==1 ,c('U','Y2')]
    col_start <- 'raw_all_'
  }
  
  tmle_Yc<- suppressMessages( ltmle(data=Oc.sub, Anodes="U", Ynodes='Y2',abar = 1,
                                    estimate.time=F))
  summ<-summary( tmle_Yc)
  
  
  Yc<- data.frame(cbind(
    pt = as.numeric(tmle_Yc$estimates['tmle']),
    lo = summ$treatment$CI[1,1],
    hi = summ$treatment$CI[1,2]
  )) 
  colnames(Yc) <- paste0(col_start, colnames(Yc))
  Yc
}





