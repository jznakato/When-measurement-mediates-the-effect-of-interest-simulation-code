
##=======================================================================================================

# looper: function to implement all estimators
## Inputs:
# dgp: Data generating process
# J: Number of clusters
# N_mean: cluster size mean
# N_sd: cluster size sd
# SL.library: SuperLearner library
## Output: A list of estimates from all estimators with inference

##=======================================================================================================

looper<- function(dpg, J, N_mean, N_sd, psi=NULL, 
                  SL.library, verbose=F){
  
  
  # full data including counterfactuals
  O <- get_full_data(dgp=dgp, J=J, N_mean=N_mean, N_sd=N_sd, verbose=F)
  
  indv_cov<-c("W1","W2","W3")  
  clust_cov<-c("E1","E2")
  
  screened <- eligible <- unadj <- tmle <- tmle_tmle <-  
    gee <- stmle_eligible <-  stmle_ratio <-  NULL
  
  
  # TWO-STAGE ESTIMATORS
  
  # do stage1
  Yc <- getYc(O=O, indv_cov=indv_cov, clust_cov=clust_cov, SL.library=SL.library)
  
  # do stage2
  goal <- 'RD' # absolute scale
  target <- 'clust' # weighting clusters equally
  psi <- mean(Yc$Yc1) - mean(Yc$Yc0) #sample effect
  
  # screened-unadjusted
  data.input <- cbind(Yc[,c('id','alpha','U','A')], Y=Yc$screened)
  screened <- Stage2(goal=goal, target=target, data.input=data.input,
                     do.unadjusted=T,  psi=psi)
  
  # eligible-unadjusted
  data.input <- cbind(Yc[,c('id','alpha','U','A')], Y=Yc$eligible)
  eligible <- Stage2(goal=goal, target=target, data.input=data.input,
                     do.unadjusted=T,  psi=psi)
  
  # unadjusted-unadjusted
  data.input <- cbind(Yc[,c('id','alpha','U','A')], Y=Yc$unadj)
  unadj <- Stage2(goal=goal, target=target, data.input=data.input,
                  do.unadjusted=T,  psi=psi)
  # tmle-unadjusted
  data.input <- cbind(Yc[,c('id','alpha','U','A')], Y=Yc$tmle)
  tmle <- Stage2(goal=goal, target=target, data.input=data.input,
                 do.unadjusted=T,  psi=psi)
  
  # tmle-tmle
  cand_cov <- c('U','E1','E2')
  data.input <- cbind(Yc[,c('id','alpha','A', cand_cov)], Y=Yc$tmle)
  tmle_tmle <- Stage2(goal=goal, target=target, data.input=data.input,
                      do.data.adapt =T, do.cv.variance=T,
                      cand.QAdj=cand_cov, cand.Qform='glm',
                      cand.gAdj=cand_cov, cand.gform='glm',
                      verbose=F, psi=psi)
  
  # SINGLE-STAGE APPPROACHES
  gee <- run_gee(O, psi)
  
  stmle_eligible <- single_tmle_eligible(O=O,
                                         SL.library=SL.library,
                                         psi=psi)
  
  stmle_ratio <- single_tmle_ratio(O=O,
                                   SL.library=SL.library,
                                   psi=psi)
  
  return(list(screened=screened, eligible=eligible,
              unadj=unadj, tmle=tmle, tmle_tmle=tmle_tmle,
              gee=gee, 
              stmle_eligible=stmle_eligible,
              stmle_ratio=stmle_ratio
  ))
}



