
##=======================================================================================================
# get.inference.1: function to calculate two-sided confidence intervals
#     & test the null hypothesis with a one-sided test
#	input: 
#		goal (aRR= arithmetic risk ratio; otherwise RD)
#   psi (true value)
#   psi.hat (estimate)
#   se (standard error)
#		df (degrees of freedom if using a Student's t-dist ) 
#		sig.level (significance level)
#   one.sided (if one-sided test)
# output: 
#		variance, test statistic, confidence intervals, pval, indicator reject null
# 		note: if goal=aRR, variance & test stat are on log-scale
##=======================================================================================================
# function for getting inference for Risk ratio estimators
get.inference.1<- function(goal='RD', psi=NA, psi.hat, se, df=99, sig.level=0.05, 
                          one.sided=F, alt.smaller=NULL){
  
  # if doing a one-sided test, need to specify the alternative
  # alt.smaller=T if intervention reduces mean outcome
  # alt.smaller=F if intervention increases mean outcome
  if(one.sided & is.null(alt.smaller)){
    print('*****ERROR: For one-sided test, need to specify the direction of the hypo')
  }
  
  # test statistic (on the log-transformed scale if goal= aRR or OR )
  tstat <- psi.hat/se
  
  if(df>40){
    # assume normal distribution
    cutoff <- qnorm(sig.level/2, lower.tail=F)
    # one.sided hypothesis test 
    if(one.sided){
      pval<- pnorm(tstat, lower.tail=alt.smaller) 
    } else{
      pval<- 2*pnorm(abs(tstat), lower.tail=F) 
    }
  }else{
    # use Student's t-distribution
    # print('Using t-distribution')
    cutoff <- qt(sig.level/2, df=df, lower.tail=F)
    # one.sided hypothesis test 
    if(one.sided){
      pval <- pt(tstat, df=df, lower.tail= alt.smaller ) 
    } else{
      pval <- 2*pt(abs(tstat), df=df, lower.tail=F)
    }
  }
  
  
  
  # 95% confidence interval 
  CI.lo <- (psi.hat - cutoff*se)
  CI.hi <- (psi.hat + cutoff*se)
  
  # transform back 
  if(goal=='aRR'){
    psi.hat<- exp(psi.hat)
    CI.lo <- exp(CI.lo)
    CI.hi <- exp(CI.hi)
    bias<-(psi.hat/psi)
  }else{ 
  psi.hat<-psi.hat
  # bias
  bias <- (psi.hat - psi)
  
  # 95% confidence interval 
  CI.lo <- (psi.hat - cutoff*se)
  CI.hi <- (psi.hat + cutoff*se)
  }
  # confidence interval coverage
  cover<- ( CI.lo <= psi & psi <= CI.hi )
  # reject the null
  reject <- as.numeric( pval < sig.level  )
  
  data.frame(est=psi.hat,  CI.lo, CI.hi, se=se,  pval, bias, cover, reject)
  
}

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
  
# glm_unadjusted (TMLE with glm in stage 1 and unadjusted in stage 2)
  data.input <- cbind(Yc[,c('id','alpha','U','A')], Y=Yc$glm_pt)
  glm_unadj<- Stage2(goal=goal, target=target, data.input=data.input,
                     do.unadjusted=T,  psi=Psi.P0.RD)
  
# glm_glm - forced adjustment for W3 in glm in stage 2
  data.input <- cbind(Yc[,c('id','alpha','U','W3','A')], Y=Yc$glm_pt)
  glm_glm<- Stage2(goal=goal, target=target, data.input=data.input,
                   do.data.adapt =F, 
                   QAdj='W3',  Qform='glm', gAdj='U', gform='glm',
                   verbose=F, psi=Psi.P0.RD)
  
# tmle-unadjusted (tmle  with SL in stage 1 and unadjusted estimation in stage 2)
  
  
  data.input <- cbind(Yc[,c('id','alpha','U','A')], Y=Yc$tmle_pt)
  
  tmle_unadj_RD<- Stage2(goal="RD", target=target, data.input=data.input,
                      do.unadjusted=T, psi=Psi.P0.RD)
  
  tmle_unadj_RR<- Stage2(goal="aRR", target=target, data.input=data.input,
                      do.unadjusted=T, psi=Psi.P0.RR)
 
 # tmle-tmle (tmle with SL in stage 1 and tmle with APS in stage 2)
  data.input <- cbind(Yc[,c(clust.var, 'W1','W2','W3')], Y=Yc$tmle_pt)
  # Need to give it the unadjusted as an option
  W <- c('U','E1','E2')#,'W1','W2','W3') ##add in extra adjustment variables. 
  tmle_tmle_RD<- Stage2(goal="RD", target=target, data.input=data.input,
                      do.data.adapt =T, remove.pscore=T, 
                      cand.QAdj=W, cand.Qform='glm', cand.gAdj=W, cand.gform='glm',
                      verbose=F, psi=Psi.P0.RD)
  
  
  tmle_tmle_RR <- Stage2(goal="aRR", target=target, data.input=data.input,
                      do.data.adapt =T, remove.pscore=T, 
                      cand.QAdj=W, cand.Qform='glm', cand.gAdj=W, cand.gform='glm',
                      verbose=F, psi=Psi.P0.RR)
##do the vanilla approach
##Fit LTMLE among those who are measured and at risk. 
  data_sub<-O[O$Delta==1 & O$Y1==0,] #when you subset your obs this doesnot subset the id vector
  id_sub<-data_sub$id
  data_sub<-data_sub[,c("U","A","Y2")]
  standard_ltmle<-ltmle(data_sub, Anodes="A",  Ynodes="Y2",abar=list(1,0),SL.library="glm",id=id_sub,estimate.time = F )
  summ<-summary(standard_ltmle)
  pval<-summ$effect.measures$ATE$pvalue
psi.hat<-summ$effect.measures$ATE$estimate
  se<-summ$effect.measures$ATE$std.dev
  CI.lo<-summ$effect.measures$ATE$CI[1]
  CI.hi<-summ$effect.measures$ATE$CI[2]
  cover<- (CI.lo<Psi.P0.RD ) & (Psi.P0.RD <CI.hi)
 reject <- as.numeric( pval < 0.05  )
 bias <- (psi.hat - Psi.P0.RD)
 standard<-data.frame(est=psi.hat,  CI.lo, CI.hi, se=se,  pval, bias, cover, reject)

  
#########  ADDING  STANDARDIZED GEE ESTIMATORS (AMONG MEAUSRED/AMONG MEAUSRED ELIGIBLE)
 
 ##: We used the MRStdCRT package to implement a standardized GEE estimator for the RD and RR that 
 ## returns the CLUSTER-level effect
 #devtools::install_github("deckardt98/MRStdCRT")
# library(MRStdCRT)

 ##calculate treatment assignment probabilities for each cluster store in vec "prob"
 
 data_std<-O[O$Delta==1,]  ##among MEAUSURED
 data_std_prob<- data_std %>% group_by(id) %>% mutate(first_trt=first(A)) %>%
   ungroup() %>% mutate(prob_A_1=mean(first_trt==1,na.rm=TRUE),prob_A_0=mean(first_trt==0,na.rm=TRUE)) %>%
   mutate(assigned_value=ifelse(A==1,prob_A_1,prob_A_0))
 
 prob<- data_std_prob$assigned_value
 
data_std$cluster<-data_std$id
 # RISK DIFFERENCE 
GEE_mod <-MRStdCRT::MRStdCRT_fit(
   formula = Y2~ A+ W1 + W2+ W3 + E1 + E2 +
    + n,
   data =  data_std,
   cluster = "cluster",
   trt = "A",
   trtprob = prob,
   family=binomial(link = "logit"),
   method = "GEE",
   corstr = "independence",
   scale = "RD"
 )
psi.hat<-GEE_mod$estimate$Estimate[1]
pval<-GEE_mod$estimate$`p-value`[1]
se<-GEE_mod$estimate$`Std. Error`[1]
CI.lo<-GEE_mod$estimate$`CI lower`[1]

CI.hi<-GEE_mod$estimate$`CI upper`[1]
cover<- (CI.lo<Psi.P0.RD ) & (Psi.P0.RD <CI.hi)
reject <- as.numeric( pval < 0.05  )
bias <- (psi.hat - Psi.P0.RD)
GEE_std_RD<-data.frame(est=psi.hat,  CI.lo, CI.hi, se=se,  pval, bias, cover, reject)


## RISK RATIO
GEE_mod <-MRStdCRT::MRStdCRT_fit(
  formula = Y2~ A+ W1 + W2+ W3 + E1 + E2 +
    + n,
  data =  data_std,
  cluster = "cluster",
  trt = "A",
  trtprob = prob,
  family=binomial(link = "logit"),
  method = "GEE",
  corstr = "independence",
  scale = "RR"
)
psi.hat<-GEE_mod$estimate$Estimate[1]
pval<-GEE_mod$estimate$`p-value`[1]
se<-GEE_mod$estimate$`Std. Error`[1]
CI.lo<-GEE_mod$estimate$`CI lower`[1]
CI.hi<-GEE_mod$estimate$`CI upper`[1]
cover<- (CI.lo<Psi.P0.RR ) & (Psi.P0.RR <CI.hi)
reject <- as.numeric( pval < 0.05  )
bias <- (psi.hat / Psi.P0.RR)
GEE_std_RR<-data.frame(est=psi.hat,  CI.lo, CI.hi, se=se,  pval, bias, cover, reject)


# ADDING RISK RATIO ESTIMATORS

# GEE: Generalized estimating equations 
library(geepack)
do.gee <- function(train, ind.cov, psi, paired=F, link, dropM=F){
  
  N <- sum(!duplicated(train$id))
  id <- train$id
  pair <- train$pair
  if(paired){
    train$pair <- as.factor(train$pair)
    train <- train[, c('pair', ind.cov, 'A', 'Y')]
  }else{
    train <- train[, c(ind.cov, 'A', 'Y')]
  }
  m <- geeglm(Y~ ., data=train, family=link, id=id)
  psi.hat<- coef(m)['A']
  se <- summary(m)$coefficients['A', 'Std.err']
  get.inference.1(goal='aRR', psi=psi, psi.hat=psi.hat, se=se, df=NA)
}



# MIXED MODELS (a.k.a., random effect models)
library('lme4')
do.mixed.models <- function(train, paired=F, psi, link, do.complex, dropM=F){
  
  N <- sum(!duplicated(train$id))

   
      m <-  glmer(Y2~ A+ W1 +W2 + W3 + E1 + E2 + (1 | id),
                  data=train, family=link ) 
  
  psi.hat <- summary(m)$coef['A','Estimate']
  se <- summary(m)$coef['A','Std. Error']
  get.inference.1(goal='aRR', psi=psi, psi.hat=psi.hat, se=se, df=99)
}

GLMM_std_RR <- do.mixed.models(train=O, psi=Psi.P0.RR, link="poisson")


# DOUBLE-ROBUST GEE (DR-GEE)
library('CRTgeeDR')
do.gee.aug <- function(train, psi, link){
  N <- sum(!duplicated(train$id))
  train[train$Delta==0,'Y1']<- NA

    out <- geeDREstimation(formula=Y2~A, nameY='Y2', id='id', nameTRT='A', nameMISS='Delta',
                           data=train, family=link, 
                           model.augmentation.ctrl = Y2~ W1 + W2 + W3 + E1 + E2,
                           model.augmentation.trt = Y2 ~ W1 + W2 + W3 + E1 + E2,
                           model.weights = Delta ~ W1 + W2+ W3 ) ##removed A and Es
  
  psi.hat <- summary(out)$beta[2]
  se <- summary(out)$se.robust[2]    
  get.inference.1(goal='aRR', psi=psi, psi.hat=psi.hat, se=se, df=N ) 

}

aug_RR <- do.gee.aug(train=O, psi=Psi.P0.RR, link="poisson")



  return( list(standard=standard,GEE_std_RD=GEE_std_RD,
               GEE_std_RR=GEE_std_RR,
               GLMM_std_RR=GLMM_std_RR,aug_RR=aug_RR,raw_all_unadj=raw_all_unadj, 
               raw_eligible_unadj=raw_eligible_unadj,
               unadj_unadj=unadj_unadj, 
               glm_unadj=glm_unadj,  glm_glm=glm_glm,
               tmle_unadj_RD=tmle_unadj_RD, tmle_tmle_RD=tmle_tmle_RD, 
               tmle_unadj_RR=tmle_unadj_RR, tmle_tmle_RR=tmle_tmle_RR))
}

