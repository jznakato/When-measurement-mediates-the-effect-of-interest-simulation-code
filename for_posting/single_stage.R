# Code by Joy Nakato & Laura B. Balzer
# jznakato@berkeley.edu & laura.balzer@berkeley.edu

# SINGLE STAGE ESTIMATORS

#-----
# run_gee: implement standardized GEE for the cluster effect among screened
# via the MRStdCRT package

run_gee <- function(O, psi){
  
  # subset on screened
  data_std <- O[O$Delta==1,] 
  
  # calculate treatment assignment probabilities for each cluster 
  data_std_prob<- data_std %>% group_by(id) %>% mutate(first_trt=first(A)) %>%
    ungroup() %>% mutate(prob_A_1=mean(first_trt==1,na.rm=TRUE),
                         prob_A_0=mean(first_trt==0,na.rm=TRUE)) %>%
    mutate(assigned_value=ifelse(A==1,prob_A_1,prob_A_0))
  
  # store in vec "prob"
  prob <- data_std_prob$assigned_value
  
  data_std$cluster <- data_std$id
  
  # RISK DIFFERENCE
  gee <-MRStdCRT::MRStdCRT_fit(
    formula = Y2~ A+ W1 + W2+ W3 + E1 + E2 + n,
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
  
  psi.hat <- gee$estimate$Estimate[1]
  CI.lo <- gee$estimate$`CI lower`[1]
  CI.hi <- gee$estimate$`CI upper`[1]
  cover <- ( CI.lo <= psi & psi <= CI.hi )
  pval <- gee$estimate$`p-value`[1]
  reject <- as.numeric(pval < 0.05)

  gee_out<-data.frame(est=psi.hat,  
                      CI.lo, 
                      CI.hi, 
                      se=gee$estimate$`Std. Error`[1],  
                      pval, 
                      bias=(psi.hat-psi), 
                      cover, 
                      reject)
  
  gee_out
  
}

# single-stage TMLE among eligible(Y1=1)
single_tmle_eligible <- function(O, SL.library=NULL, psi){
  
  # subset on eligible
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
  cover <- (CI.lo <= psi & psi <= CI.hi)
  pval<- e$pvalue
  reject <- as.numeric(pval < 0.05)
  
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
