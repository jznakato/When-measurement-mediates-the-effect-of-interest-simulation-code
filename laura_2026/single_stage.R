#########  ADDING  STANDARDIZED GEE ESTIMATORS (AMONG MEAUSRED)



# SINGLE STAGE ESTIMATORS

#-----
# run_gee: implement standardized GEE for the cluster effect among screened
# via the MRStdCRT package

run_gee <- function(O, psi){
  
  # subset on screened
  data_std<-O[O$Delta==1,] 
  
  # calculate treatment assignment probabilities for each cluster 
  data_std_prob<- data_std %>% group_by(id) %>% mutate(first_trt=first(A)) %>%
    ungroup() %>% mutate(prob_A_1=mean(first_trt==1,na.rm=TRUE),
                         prob_A_0=mean(first_trt==0,na.rm=TRUE)) %>%
    mutate(assigned_value=ifelse(A==1,prob_A_1,prob_A_0))
  
  # store in vec "prob"
  prob<- data_std_prob$assigned_value
  
  data_std$cluster <- data_std$id
  # RISK DIFFERENCE
  gee <-MRStdCRT::MRStdCRT_fit(
    formula = Y2~ A+ W1 + W2+ W3 + E1 + E2 +
      + n,
    data =  data_std,
    cluster = "cluster",
    trt = "A",
    trtprob = prob,
    family=binomial(link = "logit"),
    method = "GEE",
    corstr = "independence",
    scale = "RD",
    jack=1
  )
  psi.hat<-gee$estimate$Estimate[1]
  CI.lo<-gee$estimate$`CI lower`[1]
  CI.hi<-gee$estimate$`CI upper`[1]
  cover<- ( CI.lo <= psi & psi <= CI.hi )
  pval<- gee$estimate$`p-value`[1]
  reject <- as.numeric( pval < 0.05  )

  gee_out<-data.frame(est=psi.hat,  
                         CI.lo, 
                         CI.hi, 
                         se=gee$estimate$`Std. Error`[1],  
                         pval, 
                         bias=(psi.hat - psi), 
                         cover, 
                         reject)
  
  gee_out
  
}

# single-stage TMLE among eligible(Y1=1)
single_tmle_eligible <- function(O, SL.library=NULL, psi){
  
  data_sub <- O[O$Y1==1,]
  id_sub <- data_sub$id
  data_sub <- data_sub[,c("E1","E2","W1","W2","W3","A","Y2")]
  single_tmle <-suppressMessages(suppressWarnings(
          ltmle(data=data_sub, Anodes="A",  
          Ynodes="Y2", abar=list(1,0),
          SL.library=SL.library, id=id_sub,
          estimate.time = F) ))
  e <- summary(single_tmle)$effect.measures$ATE
  
  psi.hat <- e$estimate
  CI.lo <- e$CI[1]
  CI.hi <- e$CI[2]
  cover <- ( CI.lo <= psi & psi <= CI.hi )
  pval<- e$pvalue
  reject <- as.numeric( pval < 0.05  )
  
  e_out <-data.frame(est=psi.hat,  
                      CI.lo, 
                      CI.hi, 
                      se=e$std.dev,  
                      pval, 
                      bias=(psi.hat - psi), 
                      cover, 
                      reject)
  
  e_out
}



# single-stage TMLE for the ratio
single_tmle_ratio <- function(O, SL.library=NULL, psi){
  
  
  # numerator: Prob( Y2 =1 | A=a)
  # using ltmle to get back IC estimates
  data_num <- O[,c('A','Y2')]
  num_1 <-  suppressMessages(suppressWarnings(
                ltmle(data=data_num, Anodes="A",  
                Ynodes="Y2", abar=1,
                SL.library=NULL, id=O$id,
                estimate.time = F) ))
  num_0 <-  suppressMessages(suppressWarnings(
                ltmle(data=data_num, Anodes="A",  
                Ynodes="Y2", abar=0,
                SL.library=NULL, id=O$id,
                estimate.time = F) ))
  
  # denominator: Exp[ Prob(Y1=1 | Delta=1, A=a, E, W)]
  data_den <- O[,c("E1","E2","W1","W2","W3","A","Delta","Y1")]
  data_den$Delta <- BinaryToCensoring(is.uncensored = data_den$Delta )
  den_1 <-suppressMessages(suppressWarnings(
                ltmle(data=data_den, Anodes="A", Cnodes='Delta',
                Ynodes="Y1", abar=1,
                SL.library=SL.library, id=O$id,
                estimate.time = F) ))
  
  den_0 <-suppressMessages(suppressWarnings(
    ltmle(data=data_den, Anodes="A", Cnodes='Delta',
          Ynodes="Y1", abar=0,
          SL.library=SL.library, id=O$id,
          estimate.time = F) ))

  psi.hat <- as.numeric(num_1$estimates['tmle']/den_1$estimates['tmle'] -
                        num_0$estimates['tmle']/den_0$estimates['tmle']
                        )
 
  ic_1 <- get_ic_delta(num=num_1, den_1)
  ic_0 <- get_ic_delta(num=num_0, den_0)
  
  ic <- ic_1 - ic_0
  # print(paste0('***solve EIC: ', mean(ic)))
  
  J <- length(ic)
  se <- sqrt( var(ic)/J)
  # from Stage2 functions
  e_out <- get.inference(goal='RD', psi=psi, 
                             psi.hat=psi.hat, se=se, df=(J-2))
  
  
  e_out
}

get_ic_delta <- function(num, den){
  
  # doing directly via Moore 2007 https://biostats.bepress.com/ucbbiostat/paper215/
  # instead of usual log version
  # txt-specific mean: num/den
  
  num_pt <- as.numeric(num$estimates['tmle'])
  den_pt <- as.numeric(den$estimates['tmle'])
  ic <- 1/den_pt*num$IC$tmle - num_pt/(den_pt^2)*den$IC$tmle
  
  ic
}


