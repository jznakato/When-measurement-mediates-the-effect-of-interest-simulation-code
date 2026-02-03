
##=================================================================================
# getYc: Function to get cluster level endpoint in Stage 1 
## Inputs:
# O: observed data
# ind_cov: Individual level covariates
# clust_cov: Cluster level covariates 
## Output: 
# Returns a cluster-level data set after iterating over all clusters
##=================================================================================
getYc<- function(O, indv_cov, clust_cov, SL.library=NULL){
  
  clusters <- unique(O$id)
  J <- length(clusters)
  
  stage1 <- data.frame(matrix(NA, nrow=J, ncol= 16))
  
  for(j in clusters){ 

    # subset the data on the cluster of interest 
    Oc <- O[O$id==j, ]
    
    # Prob(Y2=1 | Delta=1)
    screened <- mean(Oc[Oc$Delta==1, 'Y2'])
    
    # Prob(Y2=1 | Y1=0) 
    eligible <- mean(Oc[Oc$Y1==1, 'Y2'])
    
    # ratio based Prob(Y2=1) / E[Prob(Y2=1|Delta, W)]
    
    # doing the numerator first P(Y2=1)
    numerator <- mean(Oc$Y2)

    # unadjusted denominator  
    den_unadj <- mean(Oc[Oc$Delta==1, 'Y1'])
    
    # TMLE denominator
    den_tmle <- suppressMessages(suppressWarnings(
                                      ltmle(data=Oc[,c(indv_cov,'Delta','Y1')], 
                                      Anodes='Delta', Ynodes='Y1', abar=1,
                                      SL.library=SL.library,
                                      estimate.time=F, stratify=T)))
    den_tmle <- unname(den_tmle$estimates["tmle"])
  
    # for spot checking bias
    cf <- Oc[Oc$Y1_star==1,]
    Yc1 <- mean(cf$Y2_1)
    Yc0 <- mean(cf$Y2_0)
    A <- Oc[1,'A'] 
    Yc_true <- ifelse(A==1, Yc1, Yc0)
    temp <- cbind(id=j, U=1, alpha=1, Oc[1, c(clust_cov)], 
                  W1=mean(Oc$W1),  W2=mean(Oc$W2),  W3=mean(Oc$W3),
                  A,
                  screened, eligible, 
                  unadj=numerator/den_unadj,
                  tmle=numerator/den_tmle,
                  Yc1,Yc0,
                  Yc_true)
    stage1[j,] <- temp
    
  }  
  colnames(stage1) <- colnames(temp)
  stage1
  
}


