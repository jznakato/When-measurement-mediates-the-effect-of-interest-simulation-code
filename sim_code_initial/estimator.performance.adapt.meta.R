
##=======================================================================================================

# looper: Function that iterates over a full data set, estimating the cluster-level end points in stage 1 
# and treatment effect in stage 2
## Inputs:
# dgp: Data egenrating mechanism for Delta
# J: Number of clusters
# N.mean: cluster size mean
# N.sd: cluster size sd
# goal: whether effect is a risk difference or ratio
# target: whether effect is at cluster-level or individual level
# Psi.P0.RD: True risk difference
# Psi.P0.RR: True risk ratio
# Psi_den" True estimate of the denominator
# SL.library: SuperLearner library
## Output: A list of estimates from all estimators with inference

##=======================================================================================================

looper<- function(dpg, J, N.mean, N.sd, effect,
                  goal=goal, target='clust',
                  Psi.P0.RD=NA,Psi.P0.RR=NA, Psi_den=NA,
                  SL.library, verbose=F){
  
# full data including counterfactuals
  full.data<-get.full.data(dgp=dgp, J=J, N.mean=N.mean, N.sd=N.sd, effect=effect, verbose=F)
  
  # observed data (drop counterfactuals)
  O <- subset(full.data, select=-c(Y2.1,Y2.0,Delta.0,Delta.1,Y1_star))
  
  ind.cov<-c("W1","W2","W3")  
  clust.cov<-c("E1","E2")
  # cluster-level variables
  clust.var<-c("id","n","alpha","U","A","E1","E2")

  # do stage1
  Yc <- getYc(O=O, ind.cov=ind.cov, clust.var=clust.var)
  
  if(F){print( colMeans(Yc[,c('unadj_pt', 'glm_pt', 'tmle_pt',
                              'raw_all_pt', 'raw_eligible_pt')]) )}
  
  # do stage2 
  data.input <- cbind(Yc[,c('id','alpha','U','A')], Y=Yc$raw_all_pt)
  raw_all_unadj<- Stage2(goal=goal, target=target, data.input=data.input,
                         do.unadjusted=T,  psi=Psi.P0.RD)
  
  data.input <- cbind(Yc[,c('id','alpha','U','A')], Y=Yc$raw_eligible_pt)
  raw_eligible_unadj <- Stage2(goal=goal, target=target, data.input=data.input,
                               do.unadjusted=T,  psi=Psi.P0.RD)
  ## complete-case estimation of HIV prevalence among those measured P(Y1=1|delta=1) and unadjusted estimation in stage2 
  data.input <- cbind(Yc[,c('id','alpha','U','A')], Y=Yc$unadj_pt)
  unadj_unadj<- Stage2(goal=goal, target=target, data.input=data.input,
                       do.unadjusted=T,  psi=Psi.P0.RD)
  
# # glm_unadjusted (TMLE with glm in stage 1 and unadjusted in stage 2)
#   data.input <- cbind(Yc[,c('id','alpha','U','A')], Y=Yc$glm_pt)
#   glm_unadj<- Stage2(goal=goal, target=target, data.input=data.input,
#                      do.unadjusted=T,  psi=Psi.P0.RD)
#   
# # glm_glm - forced adjustment for W3 in glm in stage 2
#   data.input <- cbind(Yc[,c('id','alpha','U','W3','A')], Y=Yc$glm_pt)
#   glm_glm<- Stage2(goal=goal, target=target, data.input=data.input,
#                    do.data.adapt =F, 
#                    QAdj='W3',  Qform='glm', gAdj='U', gform='glm',
#                    verbose=F, psi=Psi.P0.RD)
#   
# tmle-unadjusted (tmle  with SL in stage 1 and unadjusted estimation in stage 2)
  
  
  data.input <- cbind(Yc[,c('id','alpha','U','A')], Y=Yc$tmle_pt)
  
  tmle_unadj_RD<- Stage2(goal="RD", target=target, data.input=data.input,
                      do.unadjusted=T, psi=Psi.P0.RD)
 
 # tmle-tmle (tmle with SL in stage 1 and tmle with APS in stage 2)
  data.input <- cbind(Yc[,c(clust.var, 'W1','W2','W3')], Y=Yc$tmle_pt)
  # Need to give it the unadjusted as an option
  W <- c('U','E1','E2')#,'W1','W2','W3') ##add in extra adjustment variables. 
  tmle_tmle_RD<- Stage2(goal="RD", target=target, data.input=data.input,
                      do.data.adapt =T, remove.pscore=T, 
                      cand.QAdj=W, cand.Qform='glm', cand.gAdj=W, cand.gform='glm',
                      verbose=F, psi=Psi.P0.RD)





  return( list(raw_all_unadj=raw_all_unadj, 
               raw_eligible_unadj=raw_eligible_unadj,
               unadj_unadj=unadj_unadj, 
               tmle_unadj_RD=tmle_unadj_RD, tmle_tmle_RD=tmle_tmle_RD 
              ))
}

